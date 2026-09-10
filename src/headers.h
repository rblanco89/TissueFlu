/* TissueFlu - data types and function declarations
 *
 * Layout notes (why the code looks the way it does):
 *
 *  - The tissue is stored as a *structure of arrays*: one flat array per
 *    biological variable, all resident in GPU memory.  Every kernel then reads
 *    one variable for many neighbouring cells at once, which is what the GPU
 *    is fast at.
 *
 *  - Neighbour lists use the *compressed-sparse-row* (CSR) layout: a single
 *    flat array of links plus one offset per cell.  A cell may have any number
 *    of neighbours (no fixed maximum) and its links sit contiguously in memory.
 *
 *  - Cells are re-ordered by spatial proximity at start-up (see tissue.cu), so
 *    that neighbours in the tissue are also neighbours in memory.  Each cell
 *    keeps the line number it came from in the input file (cellId), and its
 *    random numbers are derived from that number, so the results of a run do
 *    not depend on the internal ordering.
 */
#ifndef TISSUEFLU_HEADERS_H
#define TISSUEFLU_HEADERS_H

#include <stdio.h>
#include <math.h>
#include <cuda_runtime.h>

#define THREADS_PER_BLOCK 256

/*==========================================*/
/* Saturating (Hill) response               */
/*==========================================*/

/* h(x) = x^n / (K^n + x^n):  0 at x = 0, 1/2 at x = K, -> 1 for large x.
   K^n is cached because K and n never change during a run. */
typedef struct
{
	float K;   // half-maximum: input level giving a half-maximal response
	float n;   // steepness (Hill coefficient)
	float Kn;  // cached K^n, filled in by parse_parameters()
}
Hill;

__host__ __device__ inline float hill(const Hill h, float x)
{
	if (x <= 0.0f) return 0.0f;
	// n == 2 is the usual choice and avoids a costly general power
	float xn = (h.n == 2.0f) ? x*x : powf(x, h.n);
	return xn / (h.Kn + xn);
}

__host__ Hill make_hill(float K, float n);

/*==========================================*/
/* Parameters                               */
/*==========================================*/

typedef struct
{
	/* Simulation control */
	int timeSteps;        // number of 1-minute steps to simulate
	int numReplicates;
	int measureInterval;  // steps between tissue-wide measurements
	int snapInterval;     // steps between snapshots
	int printSnap;        // 0 = none, 1 = whole tissue, 2 = infected cells only
	int ranSeed;          // seeds everything that differs between replicates
	int tissueSeed;       // if >= 0: fixed layout of non-permissive cells

	/* Tissue and initial condition */
	int numInfections;    // cells seeded with virions at t = 0
	int infectingPeriod;  // mean lifetime of an infected cell (min)
	float neighRadius;    // cells closer than this exchange virions and IFN
	float nonPermProb;    // fraction of cells that can never be infected
	float initialVirions; // virions deposited on each seeded cell
	float IFNcellProb;    // fraction of infections that produce IFN

	/* Mechanism switches (each turns one biological feedback on or off) */
	int flagRefrac;       // IFN turns susceptible cells refractory
	int flagSupp;         // IFN suppresses new infections
	int flagBP;           // IFN blocks virion production ("BP")
	int flagPF;           // IFN amplifies its own production ("PF")
	int flagPorousDiff;   // porous-medium (non-linear) transport
	int flagCTL;          // cytotoxic T lymphocytes kill infected cells

	/* Extracellular transport, per minute */
	float virionDiffusion;
	float virionClearance;
	float IFNdiffusion;
	float IFNclearance;

	/* Intracellular IFN circuit */
	float pFmax;          // IFN production rate     (IFN / min / dsRNA)
	float k_syn;          // dsRNA synthesis rate    (dsRNA / min / virion)
	float k_deg;          // dsRNA degradation rate  (1 / min)
	float alpha_pf;       // strength of the positive-feedback amplification

	/* Systemic CTL compartment */
	float rho_T;          // CTL expansion rate  (1 / min)
	float delta_T;        // CTL decay rate      (1 / min)
	float T0;             // CTL level before the infection
	float T_max;          // CTL carrying capacity

	/* Dose-response curves; K_* and nHill come from the configuration file */
	float nHill;
	Hill refrac;          // K_r   : log10(IFN+1) -> susceptible turns refractory
	Hill suppress;        // K_s   : log10(IFN+1) -> new infection is blocked
	Hill infect;          // K_v   : log10(virions+1) -> cell becomes infected
	Hill blockProd;       // K_bp  : IFN -> virion production is blocked
	Hill posFeed;         // K_pf  : IFN -> IFN production is amplified
	Hill ctlKill;         // K_T   : CTL level -> infected cell is killed
	Hill ctlGrowth;       // K_ifn : mean IFN -> CTL expansion
}
Params;

extern Params params;

/*==========================================*/
/* Tissue                                   */
/*==========================================*/

typedef enum
{
	NONPERMISSIVE,
	SUSCEPTIBLE,
	INFECTED_MINUS,
	INFECTED_PLUS,
	REFRACTORY,
	DEAD,
	NUM_STATES
}
CellState;

/* One entry of a neighbour list: who, and how strongly it is coupled. */
typedef struct __align__(8)
{
	int index;    // neighbouring cell
	float weight; // 1/d^2, d = distance between the two cells
}
Neighbor;

/* All simulation state, in GPU memory. */
typedef struct
{
	int numCells;
	int numLinks;           // total number of neighbour entries

	/* Extracellular fields.  field.x = virions, field.y = IFN.  The two
	   components travel together, so they are stored together and fetched
	   with a single memory transaction.  Diffusion reads "field" and writes
	   "fieldNext"; the two are swapped at the end of every step. */
	float2 *field;
	float2 *fieldNext;

	/* Per-cell state */
	CellState *state;
	int   *infectingTime;   // minutes left before the cell dies of infection
	int   *internalTime;    // minutes elapsed since the cell was infected
	float *dsRNA;           // viral double-stranded RNA, the IFN trigger
	int   *cellId;          // line in the input file; seeds this cell's RNG
	float3 *position;

	/* Neighbour lists, CSR layout: cell i owns link[linkStart[i] .. linkStart[i+1]) */
	int *linkStart;         // numCells + 1 entries
	Neighbor *link;         // numLinks entries
}
Tissue;

/* Host-side mirror, used only for input and output. */
typedef struct
{
	float3 *position;       // in simulation order
	int *cellId;            // simulation slot -> line in the input file
	int *slotOf;            // line in the input file -> simulation slot
	float2 *field;          // staging buffer for snapshots
	CellState *state;       // staging buffer for snapshots
}
TissueHost;

/* Tissue-wide observables at one point in time. */
typedef struct
{
	double virions;         // total extracellular virions
	double IFN;             // total interferon
	int counts[NUM_STATES]; // number of cells in each state
	float T_sys;            // systemic CTL level
}
TissueState;

/*==========================================*/
/* Functions                                */
/*==========================================*/

/* params.cu */
void parse_parameters(const char *filename);

/* tissue.cu */
int  tissue_create(const char *structureFile, const Params *pars,
                   Tissue *t, TissueHost *h);
void tissue_destroy(Tissue *t, TissueHost *h);

/* model.cu - the biological model */
__global__ void tissue_update(Tissue t, Params pars, const float *T_sys,
                              unsigned long long seed, int step);
__global__ void tissue_diffusion(Tissue t, Params pars, float *IFNblock);
__global__ void update_Tsys(Params pars, float *T_sys, const float *IFNblock,
                            int numBlocks, int numCells);
__global__ void tissue_reduce(Tissue t, double *blockVirions, double *blockIFN,
                              int *counts);

/* output.cu */
void measure_open(int numCells);
void measure_close(void);
void measure_tissue(const Tissue *t, const float *d_T_sys, TissueState *out);
void record_measurement(const TissueState *s, int measureIdx, int measureInterval,
                        double *aucVirus, double *aucIFN,
                        double *prevVirus, double *prevIFN, FILE *fRep);
void print_snapshot(const Tissue *t, TissueHost *h, FILE *fSnap, int onlyInfected);

#endif
