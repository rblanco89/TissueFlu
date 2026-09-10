/* TissueFlu: simulation of influenza spread in the lungs of mice
   Author: Rodolfo Blanco
   Date: April 2026

   One time step is one minute.  The tissue lives in GPU memory for the whole
   run; the host only sets up the initial condition and collects measurements.
   The biological model itself is in model.cu. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <errno.h>

#include <cuda_runtime.h>

#include "headers.h"
#include "ranNumbers.h"

int main(int argc, char *argv[])
{
	const char *config_file = NULL;
	const char *structure_file = NULL;
	const char *output_dir = ".";

	/*==========================================*/
	// Parse command-line arguments
	/*==========================================*/

	int i = 1;
	while (i < argc)
	{
		if (strncmp(argv[i], "--config", 8) == 0)
		{
			config_file = argv[++i];
			i++;
		}
		else if (strncmp(argv[i], "--structure", 11) == 0)
		{
			structure_file = argv[++i];
			i++;
		}
		else if (strncmp(argv[i], "--output", 8) == 0)
		{
			output_dir = argv[++i];
			i++;
		}
		else
		{
			fprintf(stderr, "Unknown argument: %s\n", argv[i]);
			return 1;
		}
	}

	if (config_file == NULL || structure_file == NULL)
	{
		fprintf(stderr,
				"Usage: %s --config FILE --structure FILE [--output DIR]\n"
				"  --config     : path to configuration file\n"
				"  --structure  : path to cell positions file (e.g., CSV)\n"
				"  --output     : directory where results will be written (default is in-place)\n"
				"  Neighbor list is loaded automatically from <structure>_neighbors_r<R>.csv\n"
				"  next to the structure file; if absent it is computed and saved there.\n",
				argv[0]);
		return 1;
	}

	/*==========================================*/
	// Build the tissue
	/*==========================================*/

	parse_parameters(config_file);

	Tissue tissue;
	TissueHost host;
	if (!tissue_create(structure_file, &params, &tissue, &host)) return 1;

	int numCells = tissue.numCells;
	int blocks = (numCells + THREADS_PER_BLOCK - 1)/THREADS_PER_BLOCK;

	/* Systemic CTL level, and the per-block IFN totals that drive it.  Both
	   live on the GPU: the host never needs them between measurements. */
	float *d_T_sys, *d_IFNblock;
	cudaMalloc(&d_T_sys, sizeof(float));
	cudaMalloc(&d_IFNblock, blocks*sizeof(float));

	measure_open(numCells);

	/* Staging buffers for the initial condition, filled on the host once per
	   replicate and uploaded in one go (pinned memory makes the copy faster). */
	float2 *h_field;  CellState *h_state;
	int *h_infectingTime, *h_internalTime;
	float *h_dsRNA;
	cudaMallocHost(&h_field, numCells*sizeof(float2));
	cudaMallocHost(&h_state, numCells*sizeof(CellState));
	cudaMallocHost(&h_infectingTime, numCells*sizeof(int));
	cudaMallocHost(&h_internalTime, numCells*sizeof(int));
	cudaMallocHost(&h_dsRNA, numCells*sizeof(float));

	/*==========================================*/
	// Output files
	/*==========================================*/

	if (mkdir(output_dir, 0755) != 0 && errno != EEXIST)
	{
		fprintf(stderr, "Error: could not create output directory %s\n", output_dir);
		return 1;
	}

	char path_snapshots[512], path_cellState[512], path_auc[512];
	snprintf(path_snapshots, sizeof(path_snapshots), "%s/snapshots.xyz",   output_dir);
	snprintf(path_cellState, sizeof(path_cellState), "%s/cellState.csv",   output_dir);
	snprintf(path_auc,       sizeof(path_auc),       "%s/auc_results.csv", output_dir);

	const char *header = "Time,Virus,IFN,Susceptible,Refractory,Infected,Dead,nonPermissive,Tcells\n";

	if (params.tissueSeed >= 0)
		printf("Non-permissive cells: fixed layout from tissueSeed = %d\n", params.tissueSeed);
	else
		printf("Non-permissive cells: new layout for every replicate\n");

	FILE *fAUC = fopen(path_auc, "w");
	if (!fAUC) fprintf(stderr, "Warning: could not open %s for writing\n", path_auc);
	else fprintf(fAUC, "Replicate,AUC_Virus,AUC_IFN\n");

	FILE *fSnap = NULL;
	if (params.numReplicates == 1 && params.printSnap)
	{
		fSnap = fopen(path_snapshots, "w");
		if (!fSnap) fprintf(stderr, "Warning: could not open %s for writing\n", path_snapshots);
	}

	FILE *fCell = NULL;
	if (params.numReplicates == 1)
	{
		fCell = fopen(path_cellState, "w");
		if (!fCell) fprintf(stderr, "Warning: could not open %s for writing\n", path_cellState);
		else fprintf(fCell, "Time,ViralLoad,IFN\n");
	}

	/*==========================================*/
	// Replicates Loop
	/*==========================================*/

	for (int rep = 0; rep < params.numReplicates; rep++)
	{
		printf("\nRunning replicate %d/%d\n", rep+1, params.numReplicates);

		// Time course of this replicate (averaging is left to the analysis scripts)
		char path_rep[512];
		snprintf(path_rep, sizeof(path_rep), "%s/tissueState_%d.csv", output_dir, rep+1);
		FILE *fRep = fopen(path_rep, "w");
		if (!fRep) fprintf(stderr, "Warning: could not open %s for writing\n", path_rep);
		else fprintf(fRep, "%s", header);

		/*------------------------------------------*/
		// Initial condition
		/*------------------------------------------*/

		unsigned long long seed = (unsigned long long)params.ranSeed + rep;
		Ran ranUni(seed);
		Poissondev ranInfecting(params.infectingPeriod/60, seed^0x9E3779B9); // hours

		// Which cells are non-permissive.  With tissueSeed set, a generator
		// restarted from that seed gives the same layout in every replicate and
		// every run; otherwise the layout is part of the replicate's randomness.
		Ran ranTissue(params.tissueSeed >= 0 ? params.tissueSeed : 0);
		Ran &ranLayout = (params.tissueSeed >= 0) ? ranTissue : ranUni;

		// Drawn in the order of the structure file, so the initial condition
		// does not depend on the internal ordering of the cells.
		for (int o = 0; o < numCells; o++)
		{
			int c = host.slotOf[o];
			h_state[c] = (ranLayout.doub() < params.nonPermProb) ? NONPERMISSIVE : SUSCEPTIBLE;
			h_field[c] = make_float2(0.0f, 0.0f);
			h_dsRNA[c] = 0.0f;
			h_infectingTime[c] = 60*ranInfecting.dev();   // hours -> minutes
			h_internalTime[c] = 0;
		}

		// Deposit the inoculum on randomly chosen cells (differs per replicate)
		int tracked = 0;
		for (int n = 0; n < params.numInfections; n++)
		{
			int o;
			do o = numCells*ranUni.doub();
			while (h_field[host.slotOf[o]].x > 0.0f);
			h_field[host.slotOf[o]].x = params.initialVirions;
			tracked = host.slotOf[o];   // cellState.csv follows the last one
		}

		cudaMemcpy(tissue.field, h_field, numCells*sizeof(float2), cudaMemcpyHostToDevice);
		cudaMemcpy(tissue.state, h_state, numCells*sizeof(CellState), cudaMemcpyHostToDevice);
		cudaMemcpy(tissue.dsRNA, h_dsRNA, numCells*sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(tissue.infectingTime, h_infectingTime, numCells*sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(tissue.internalTime, h_internalTime, numCells*sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(d_T_sys, &params.T0, sizeof(float), cudaMemcpyHostToDevice);

		/*==========================================*/
		// Simulation Loop
		/*==========================================*/

		int measureIdx = 0;
		double aucVirus = 0.0, aucIFN = 0.0;
		double prevVirus = 0.0, prevIFN = 0.0;
		int progressInterval = params.timeSteps/10;
		if (progressInterval == 0) progressInterval = 1;

		for (int step = 0; step <= params.timeSteps; step++)
		{
			if (step % progressInterval == 0) { printf("."); fflush(stdout); }

			if (step % params.measureInterval == 0)
			{
				TissueState now;
				measure_tissue(&tissue, d_T_sys, &now);
				record_measurement(&now, measureIdx++, params.measureInterval,
				                   &aucVirus, &aucIFN, &prevVirus, &prevIFN, fRep);

				if (fCell)
				{
					float2 one;
					cudaMemcpy(&one, tissue.field + tracked, sizeof(float2),
					           cudaMemcpyDeviceToHost);
					fprintf(fCell, "%d,%f,%f\n", step, one.x, one.y);
				}
			}

			if (fSnap && step % params.snapInterval == 0)
				print_snapshot(&tissue, &host, fSnap, params.printSnap == 2);

			// 1. what happens inside each cell
			tissue_update<<<blocks, THREADS_PER_BLOCK>>>(tissue, params, d_T_sys, seed, step);

			// 2. transport of virions and IFN between cells
			tissue_diffusion<<<blocks, THREADS_PER_BLOCK>>>(tissue, params, d_IFNblock);

			// 3. the systemic CTL response
			if (params.flagCTL)
				update_Tsys<<<1, THREADS_PER_BLOCK>>>(params, d_T_sys, d_IFNblock,
				                                      blocks, numCells);

			// diffusion wrote the new fields into fieldNext; make them current
			float2 *swap = tissue.field; tissue.field = tissue.fieldNext; tissue.fieldNext = swap;
		}

		if (fRep) fclose(fRep);
		if (fAUC) fprintf(fAUC, "%d,%e,%e\n", rep+1, aucVirus, aucIFN);
	}

	if (fCell) fclose(fCell);
	if (fSnap) fclose(fSnap);
	if (fAUC)  fclose(fAUC);

	printf("\nCompleted\n");

	/*==========================================*/
	// Clean up
	/*==========================================*/

	cudaFreeHost(h_field); cudaFreeHost(h_state); cudaFreeHost(h_dsRNA);
	cudaFreeHost(h_infectingTime); cudaFreeHost(h_internalTime);
	cudaFree(d_T_sys); cudaFree(d_IFNblock);
	measure_close();
	tissue_destroy(&tissue, &host);

	return 0;
}
