/* TissueFlu - the biological model
 *
 * One time step of the simulation is one minute and consists of three kernels,
 * launched from main.cu in this order:
 *
 *   1. tissue_update    - what happens inside and to each cell
 *   2. tissue_diffusion - how virions and IFN spread between cells
 *   3. update_Tsys      - the systemic (whole-animal) CTL response
 *
 * Every kernel maps one GPU thread to one cell.
 */
#include <curand_kernel.h>
#include "headers.h"

/*==========================================*/
/* Random numbers                           */
/*==========================================*/

/* Three independent uniform deviates in (0,1] for one cell at one time step.
 *
 * Philox is a counter-based generator: the deviates are computed directly from
 * (seed, cell, step) instead of being read from a stored stream.  Two useful
 * consequences: nothing has to be stored or transferred, and a given cell at a
 * given step always draws the same numbers, whatever the internal ordering of
 * the cells or the number of GPU threads.  A run is therefore reproducible
 * from its seed alone.
 */
__device__ inline float4 cell_random(unsigned long long seed, int cellId, int step)
{
	curandStatePhilox4_32_10_t rng;
	curand_init(seed, cellId, 4ull*step, &rng);
	return curand_uniform4(&rng);
}

/*==========================================*/
/* Virion release by an infected cell       */
/*==========================================*/

/* Cumulative release follows a Hill curve in time, V(t) = A t^n / (K^n + t^n);
   a cell releases dV/dt virions during the current minute. */
__device__ inline float virion_production(int minutesInfected)
{
	const float A  = 17703.6603f;   // total virions released by one cell
	const float n  = 2.3904f;       // steepness
	const float Kn = 13710739.0f;   // K^n, with K = 16.1275 h = 967.65 min
	                                // the time of half-maximal release

	float t = (float)minutesInfected;
	if (t <= 0.0f) return 0.0f;

	float tn  = powf(t, n);
	float den = Kn + tn;
	return A * n * Kn * (tn/t) / (den*den);   // tn/t is t^(n-1)
}

/*==========================================*/
/* 1. Cell state machine                    */
/*==========================================*/

__global__ void tissue_update(Tissue t, Params pars, const float *T_sys,
                              unsigned long long seed, int step)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	if (i >= t.numCells) return;

	CellState state = t.state[i];
	if (state == NONPERMISSIVE || state == REFRACTORY || state == DEAD) return;

	float4 ran = cell_random(seed, t.cellId[i], step);
	float2 field = t.field[i];          // .x = virions, .y = IFN around the cell

	switch (state)
	{
		case SUSCEPTIBLE:
		{
			// Refractory mechanism: IFN puts the cell into an antiviral state
			if (pars.flagRefrac && hill(pars.refrac, log10f(field.y + 1.0f)) > ran.x)
			{
				t.state[i] = REFRACTORY;
				break;
			}

			// Infection: driven by the local virion load, damped by local IFN
			float infecProb = hill(pars.infect, log10f(field.x + 1.0f));
			if (pars.flagSupp)
				infecProb *= 1.0f - hill(pars.suppress, log10f(field.y + 1.0f));

			if (infecProb > ran.y)
				t.state[i] = (ran.z < pars.IFNcellProb) ? INFECTED_PLUS   // makes IFN
				                                        : INFECTED_MINUS; // makes none
			break;
		}

		case INFECTED_PLUS:
		case INFECTED_MINUS:
		{
			// CTL killing
			if (pars.flagCTL && hill(pars.ctlKill, *T_sys) > ran.x)
			{
				t.state[i] = DEAD;
				break;
			}

			// Virion production, optionally blocked by IFN (BP mechanism)
			float produced = virion_production(t.internalTime[i]++);
			if (pars.flagBP) produced *= 1.0f - hill(pars.blockProd, field.y);
			field.x += produced;

			// Only IFN-competent cells run the dsRNA -> IFN circuit
			if (state == INFECTED_PLUS)
			{
				float dsRNA = t.dsRNA[i];
				dsRNA += pars.k_syn*produced - pars.k_deg*dsRNA;
				t.dsRNA[i] = dsRNA;

				float ifn = pars.pFmax*dsRNA;
				// Positive-feedback mechanism: IFN amplifies its own production
				if (pars.flagPF) ifn *= 1.0f + pars.alpha_pf*hill(pars.posFeed, field.y);
				field.y += ifn;
			}

			t.field[i] = field;

			// The cell dies once its infectious lifetime runs out
			if (--t.infectingTime[i] <= 0) t.state[i] = DEAD;
			return;
		}

		default:
			break;
	}
}

/*==========================================*/
/* 2. Extracellular transport               */
/*==========================================*/

/* Virions and IFN hop between cells that are within neighRadius of each other,
 * at a rate weighted by 1/d^2, and are cleared at a constant rate.  Reading
 * "field" and writing "fieldNext" keeps every cell seeing the same state of the
 * tissue, independently of the order in which the GPU happens to run threads.
 *
 * The kernel also accumulates the tissue-wide IFN total needed by the CTL
 * compartment: each block sums its own cells and leaves one partial total for
 * update_Tsys to add up.  Summing in a fixed order (rather than with atomics,
 * which arrive in whatever order the blocks happen to finish) keeps a run
 * reproducible from its seed.
 */
__global__ void tissue_diffusion(Tissue t, Params pars, float *IFNblock)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;

	float2 self = make_float2(0.0f, 0.0f);
	if (i < t.numCells)
	{
		self = t.field[i];

		// Porous-medium transport moves f^2 instead of f
		bool porous = (pars.flagPorousDiff != 0);
		float2 from = porous ? make_float2(self.x*self.x, self.y*self.y) : self;

		float2 flux = make_float2(0.0f, 0.0f);
		int begin = t.linkStart[i], end = t.linkStart[i+1];
		for (int k = begin; k < end; k++)
		{
			Neighbor nb = t.link[k];
			float2 to = t.field[nb.index];
			if (porous) { to.x *= to.x; to.y *= to.y; }

			flux.x += nb.weight*(to.x - from.x);
			flux.y += nb.weight*(to.y - from.y);
		}

		float aux = porous ? 0.5f : 1.0f;   // from d(f^2)/dx = 2 f df/dx
		float2 next;
		next.x = (1.0f - pars.virionClearance)*(self.x + aux*pars.virionDiffusion*flux.x);
		next.y = (1.0f - pars.IFNclearance)   *(self.y + aux*pars.IFNdiffusion   *flux.y);
		t.fieldNext[i] = next;
	}

	if (!pars.flagCTL) return;

	// Block-wide sum of the IFN present before this step's transport
	__shared__ float partial[THREADS_PER_BLOCK];
	partial[threadIdx.x] = self.y;
	__syncthreads();

	for (int stride = blockDim.x/2; stride > 0; stride >>= 1)
	{
		if (threadIdx.x < stride) partial[threadIdx.x] += partial[threadIdx.x + stride];
		__syncthreads();
	}
	if (threadIdx.x == 0) IFNblock[blockIdx.x] = partial[0];
}

/*==========================================*/
/* 3. Systemic CTL compartment              */
/*==========================================*/

/* A single well-mixed pool of cytotoxic T cells shared by the whole tissue:
   it expands when the mean IFN level is high and decays otherwise. */
__global__ void update_Tsys(Params pars, float *T_sys, const float *IFNblock,
                            int numBlocks, int numCells)
{
	// Add up the per-block IFN totals left by tissue_diffusion
	__shared__ double partial[THREADS_PER_BLOCK];
	double sum = 0.0;
	for (int k = threadIdx.x; k < numBlocks; k += THREADS_PER_BLOCK) sum += IFNblock[k];
	partial[threadIdx.x] = sum;
	__syncthreads();

	for (int stride = THREADS_PER_BLOCK/2; stride > 0; stride >>= 1)
	{
		if (threadIdx.x < stride) partial[threadIdx.x] += partial[threadIdx.x + stride];
		__syncthreads();
	}
	if (threadIdx.x != 0) return;

	float T = *T_sys;
	float IFNmean = (float)(partial[0]/numCells);

	float expansion = pars.rho_T * hill(pars.ctlGrowth, IFNmean) * T
	                * (1.0f - T/pars.T_max);        // logistic cap at T_max
	T += expansion - pars.delta_T*T;

	*T_sys = fmaxf(T, 0.0f);
}

/*==========================================*/
/* Tissue-wide measurement                  */
/*==========================================*/

/* Per-block totals of virions and IFN, and the number of cells in each state.
   Run only on measurement steps.  The host adds the block totals up in index
   order, so the reported numbers do not depend on GPU scheduling. */
__global__ void tissue_reduce(Tissue t, double *blockVirions, double *blockIFN,
                              int *counts)
{
	__shared__ double sumVirions[THREADS_PER_BLOCK];
	__shared__ double sumIFN[THREADS_PER_BLOCK];
	__shared__ int    tally[NUM_STATES];

	int i = threadIdx.x + blockIdx.x*blockDim.x;

	if (threadIdx.x < NUM_STATES) tally[threadIdx.x] = 0;
	sumVirions[threadIdx.x] = 0.0;
	sumIFN[threadIdx.x]     = 0.0;
	__syncthreads();

	if (i < t.numCells)
	{
		float2 field = t.field[i];
		sumVirions[threadIdx.x] = field.x;
		sumIFN[threadIdx.x]     = field.y;
		atomicAdd(&tally[t.state[i]], 1);
	}
	__syncthreads();

	for (int stride = blockDim.x/2; stride > 0; stride >>= 1)
	{
		if (threadIdx.x < stride)
		{
			sumVirions[threadIdx.x] += sumVirions[threadIdx.x + stride];
			sumIFN[threadIdx.x]     += sumIFN[threadIdx.x + stride];
		}
		__syncthreads();
	}

	if (threadIdx.x == 0)
	{
		blockVirions[blockIdx.x] = sumVirions[0];
		blockIFN[blockIdx.x]     = sumIFN[0];
	}
	if (threadIdx.x < NUM_STATES && tally[threadIdx.x])
		atomicAdd(&counts[threadIdx.x], tally[threadIdx.x]);
}
