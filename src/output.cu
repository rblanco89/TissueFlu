/* TissueFlu - measurements and output files */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "headers.h"

/* Accumulators the reduction kernel writes into, plus a pinned host buffer to
   copy them back through.  Allocated once, reused at every measurement. */
static double *d_blockSum = NULL;   // 2*numBlocks: virions, then IFN
static int    *d_counts   = NULL;
static double *h_blockSum = NULL;
static int    *h_counts   = NULL;
static float  *h_T_sys    = NULL;
static int     numBlocks  = 0;

void measure_open(int numCells)
{
	numBlocks = (numCells + THREADS_PER_BLOCK - 1)/THREADS_PER_BLOCK;
	cudaMalloc(&d_blockSum, 2*numBlocks*sizeof(double));
	cudaMalloc(&d_counts, NUM_STATES*sizeof(int));
	cudaMallocHost(&h_blockSum, 2*numBlocks*sizeof(double));
	cudaMallocHost(&h_counts, NUM_STATES*sizeof(int));
	cudaMallocHost(&h_T_sys, sizeof(float));
}

void measure_close(void)
{
	cudaFree(d_blockSum);
	cudaFree(d_counts);
	cudaFreeHost(h_blockSum);
	cudaFreeHost(h_counts);
	cudaFreeHost(h_T_sys);
}

/* Totals over the whole tissue.  The sums are computed on the GPU: bringing
   the cell arrays back to the host just to add them up would cost more than a
   whole time step. */
void measure_tissue(const Tissue *t, const float *d_T_sys, TissueState *out)
{
	cudaMemset(d_counts, 0, NUM_STATES*sizeof(int));

	tissue_reduce<<<numBlocks, THREADS_PER_BLOCK>>>(*t, d_blockSum,
	                                                d_blockSum + numBlocks, d_counts);

	cudaMemcpy(h_blockSum, d_blockSum, 2*numBlocks*sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_counts, d_counts, NUM_STATES*sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_T_sys, d_T_sys, sizeof(float), cudaMemcpyDeviceToHost);

	out->virions = 0.0;
	out->IFN     = 0.0;
	for (int b = 0; b < numBlocks; b++)
	{
		out->virions += h_blockSum[b];
		out->IFN     += h_blockSum[numBlocks + b];
	}
	memcpy(out->counts, h_counts, NUM_STATES*sizeof(int));
	out->T_sys = *h_T_sys;
}

/* Store one measurement: one row of this replicate's time course, and the
   running areas under the viral-load and IFN curves. */
void record_measurement(const TissueState *s, int measureIdx, int measureInterval,
                        double *aucVirus, double *aucIFN,
                        double *prevVirus, double *prevIFN, FILE *fRep)
{
	double val[8];
	val[0] = s->virions;
	val[1] = s->IFN;
	val[2] = s->counts[SUSCEPTIBLE];
	val[3] = s->counts[REFRACTORY];
	val[4] = s->counts[INFECTED_PLUS] + s->counts[INFECTED_MINUS];
	val[5] = s->counts[DEAD];
	val[6] = s->counts[NONPERMISSIVE];
	val[7] = s->T_sys;

	if (measureIdx > 0)   // trapezoidal rule between consecutive measurements
	{
		*aucVirus += 0.5*(s->virions + *prevVirus)*measureInterval;
		*aucIFN   += 0.5*(s->IFN     + *prevIFN)  *measureInterval;
	}
	*prevVirus = s->virions;
	*prevIFN   = s->IFN;

	if (fRep)
	{
		fprintf(fRep, "%d", measureIdx*measureInterval);
		for (int m = 0; m < 8; m++) fprintf(fRep, ",%e", val[m]);
		fprintf(fRep, "\n");
	}
}

/* Extended-XYZ snapshot for viewers such as Ovito.  Cells are written in the
   order of the structure file, not in the internal (spatially sorted) order. */
void print_snapshot(const Tissue *t, TissueHost *h, FILE *fSnap, int onlyInfected)
{
	int numCells = t->numCells;
	cudaMemcpy(h->field, t->field, numCells*sizeof(float2), cudaMemcpyDeviceToHost);
	cudaMemcpy(h->state, t->state, numCells*sizeof(CellState), cudaMemcpyDeviceToHost);

	int written = numCells;
	if (onlyInfected)
	{
		written = 0;
		for (int c = 0; c < numCells; c++)
			if (h->state[c] == INFECTED_PLUS || h->state[c] == INFECTED_MINUS) written++;
	}

	fprintf(fSnap, "%d\n", written);
	fprintf(fSnap, "Properties=species:I:1:pos:R:3:virions:R:1:IFN:R:1\n");

	for (int o = 0; o < numCells; o++)
	{
		int s = h->slotOf[o];
		CellState state = h->state[s];
		if (onlyInfected && state != INFECTED_PLUS && state != INFECTED_MINUS) continue;

		float3 r = h->position[s];
		fprintf(fSnap, "%d %f %f %f %f %f\n",
		        (int)state, r.x, r.y, r.z, h->field[s].x, h->field[s].y);
	}
}
