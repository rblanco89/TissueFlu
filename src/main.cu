/* AeroFlue: Influenza simulation
   Author: Rodolfo Blanco
   Date: April 2026 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <errno.h>

#include <cuda_runtime.h>
#include <curand.h>

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

	// Check that both required arguments are set
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
	// Fetch structure and parameters
	/*==========================================*/

	parse_parameters(config_file);

	FILE *fp = fopen(structure_file, "r");
	if (!fp)
	{
		fprintf(stderr, "Error: could not open structure file %s\n", structure_file);
		return 1;
	}

	// Count number of lines (cells)
	int numCells = 0;
	char line[256];
	while (fgets(line, sizeof(line), fp)) numCells++;
	rewind(fp); // reset file pointer
	numCells--; // Skip the header line

	if (fgets(line, sizeof(line), fp) == NULL || numCells <= 0)
	{
		fprintf(stderr, "Error: structure file is empty: %s\n", structure_file);
		fclose(fp);
		return 1;
	}

	// Allocate memory for cells
	Cell *cells;
	cudaError_t err = cudaMallocManaged((void**)&cells, numCells*sizeof(Cell), cudaMemAttachGlobal);
	if (err != cudaSuccess || cells == NULL)
	{
		fprintf(stderr, "cudaMallocManaged failed for cells: %s\n", cudaGetErrorString(err));
		fclose(fp);
		return 1;
	}

	// Read positions and initialize fields
	for (int c=0; c<numCells; c++)
	{
		float3 r;
		if (fgets(line, sizeof(line), fp) == NULL)
		{
			fprintf(stderr, "Unexpected end of structure file at cell %d\n", c);
			fclose(fp);
			cudaFree(cells);
			return 1;
		}
		
		if (sscanf(line, "%f,%f,%f", &r.x, &r.y, &r.z) != 3)
		{
			fprintf(stderr, "Invalid format in structure file at line %d\n", c + 2);
			fclose(fp);
			cudaFree(cells);
			return 1;
		}

		cells[c].position = r;
		cells[c].numNeighbors = 0;
	}
	fclose(fp);

	printf("Loaded %d cells from %s\n", numCells, structure_file);

	/*==========================================*/
	// Pre-simulation setup
	/*==========================================*/

	// Set GPU random generator
	float *d_ranUni;
	curandGenerator_t gen;
	cudaMalloc(&d_ranUni, 3*numCells*sizeof(float)); // 3 non-overlapping slots per cell
	curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_MTGP32);

	// Estimate the number of threads and blocks for the GPU
	int ths = (numCells < THS_MAX) ? nextPow2(numCells) : THS_MAX;
	int blks = 1 + (numCells - 1)/ths;

	// Derive auto-save path for neighbor list: <structure_base>_neighbors<ext>
	char neigh_path_buf[512];
	{
		const char *dot = strrchr(structure_file, '.');
		if (dot && dot != structure_file)
		{
			int base_len = (int)(dot - structure_file);
			snprintf(neigh_path_buf, sizeof(neigh_path_buf),
					"%.*s_neighbors_r%d%s", base_len, structure_file,
					int(params.neighRadius), dot);
		}
		else
		{
			snprintf(neigh_path_buf, sizeof(neigh_path_buf), "%s_neighbors_r%d",
					 structure_file, int(params.neighRadius));
		}
	}

	FILE *fNeigh = fopen(neigh_path_buf, "r");
	if (fNeigh)
	{
		char nline[MAX_NEIGHBORS * 8 + 32];
		for (int c = 0; c < numCells; c++)
		{
			if (fgets(nline, sizeof(nline), fNeigh) == NULL)
			{
				fprintf(stderr, "Error: neighbor file too short at cell %d\n", c);
				fclose(fNeigh);
				cudaFree(cells);
				cudaFree(d_ranUni);
				curandDestroyGenerator(gen);
				return 1;
			}
			char *tok = strtok(nline, ",\n");
			cells[c].numNeighbors = atoi(tok);
			for (int n = 0; n < cells[c].numNeighbors; n++)
			{
				tok = strtok(NULL, ",\n");
				if (tok == NULL)
				{
					fprintf(stderr, "Error: malformed neighbor entry for cell %d\n", c);
					fclose(fNeigh);
					cudaFree(cells);
					cudaFree(d_ranUni);
					curandDestroyGenerator(gen);
					return 1;
				}
				cells[c].neighbors[n] = atoi(tok);
			}
		}
		fclose(fNeigh);
		printf("Loaded neighbor list from %s\n", neigh_path_buf);

		// Compute 1/d^2 weights from the loaded neighbor indices
		compute_weights<<<blks, ths>>>(cells, numCells);
		cudaDeviceSynchronize();
	}
	else
	{
		printf("Creating list of neighbors...\n");

		// Find neighbors, store indices and 1/d² weights in cell structure
		build_neighbors<<<blks, ths>>>(cells, numCells, params.neighRadius);
		cudaDeviceSynchronize();

		// Save neighbor list next to the structure file
		printf("Saving neighbor list to %s...\n", neigh_path_buf);
		FILE *fSave = fopen(neigh_path_buf, "w");
		if (!fSave)
		{
			fprintf(stderr, "Warning: could not save neighbor list to %s\n", neigh_path_buf);
		}
		else
		{
			for (int c = 0; c < numCells; c++)
			{
				fprintf(fSave, "%d", cells[c].numNeighbors);
				for (int n = 0; n < cells[c].numNeighbors; n++)
					fprintf(fSave, ",%d", cells[c].neighbors[n]);
				fprintf(fSave, "\n");
			}
			fclose(fSave);
		}
	}

	float maxDiffusion = stabilityCondition(cells, numCells);
	if (maxDiffusion < params.virionDiffusion || maxDiffusion < params.IFNdiffusion)
	{
		printf("Diffusion parameters must be less than %f\nStopping...\n", maxDiffusion);
		cudaFree(cells);
		cudaFree(d_ranUni);
		curandDestroyGenerator(gen);
		return 1;
	}

	// Create output directory if it does not exist
	if (mkdir(output_dir, 0755) != 0 && errno != EEXIST)
	{
		fprintf(stderr, "Error: could not create output directory %s\n", output_dir);
		cudaFree(cells);
		cudaFree(d_ranUni);
		curandDestroyGenerator(gen);
		return 1;
	}

	// Build output file paths
	char path_snapshots[512], path_cellState[512];
	char path_tissueState[512], path_tissueAvg[512], path_tissueStd[512], path_auc[512];
	snprintf(path_snapshots,  sizeof(path_snapshots),  "%s/snapshots.xyz",       output_dir);
	snprintf(path_cellState,  sizeof(path_cellState),  "%s/cellState.csv",       output_dir);
	snprintf(path_tissueState,sizeof(path_tissueState),"%s/tissueState.csv",     output_dir);
	snprintf(path_tissueAvg,  sizeof(path_tissueAvg),  "%s/tissueState_avg.csv", output_dir);
	snprintf(path_tissueStd,  sizeof(path_tissueStd),  "%s/tissueState_std.csv", output_dir);
	snprintf(path_auc,        sizeof(path_auc),        "%s/auc_results.csv",     output_dir);

	// Prepare for result accumulation
	int numMeasures = params.timeSteps / params.measureInterval + 1;

	// Arrays for averaging tissue state (7 variables: V, IFN, S, R, I, D, NP)
	double *tissueSum = (double*)calloc(numMeasures * 7, sizeof(double));
	double *tissueSqSum = (double*)calloc(numMeasures * 7, sizeof(double));

	// Cell state counters (indexed by CellState enum): updated atomically by GPU
	int *cellCounts;
	cudaMallocManaged(&cellCounts, 6 * sizeof(int));

	// Shadow arrays for race-free diffusion (read-only snapshots of virions and IFN)
	float *virions_old, *IFN_old;
	cudaMalloc(&virions_old, numCells*sizeof(float));
	cudaMalloc(&IFN_old,     numCells*sizeof(float));

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
		fprintf(fCell, "Time,ViralLoad,IFN\n");
	}

	/*==========================================*/
	// Replicates Loop
	/*==========================================*/

	for (int rep=0; rep<params.numReplicates; rep++)
	{
		printf("\nRunning replicate %d/%d\n", rep+1, params.numReplicates);

		// Initialize/Reset random numbers
		ulong seed = params.ranSeed + rep;
		Ran ranUni(seed);
		Poissondev ranInfecting(params.infectingPeriod/60, seed^0x9E3779B9);
		curandSetPseudoRandomGeneratorSeed(gen, seed);

		// Reset cells state
		memset(cellCounts, 0, 6 * sizeof(int));
		for (int i=0; i<numCells; i++)
		{
			if (ranUni.doub() < params.nonPermProb) cells[i].state = NONPERMISSIVE;
			else cells[i].state = SUSCEPTIBLE;
		
			cellCounts[cells[i].state]++;

			cells[i].virions = 0.0f;
			cells[i].dsRNA = 0.0f;
			cells[i].IFN = 0.0f;
			cells[i].infectingTime = 60*ranInfecting.dev(); // Convert to minutes
			cells[i].internalTime = 0;
		}

		int ind = 0;
		for (int i=0; i<params.numInfections; i++)
		{
			do ind = numCells*ranUni.doub();
			while (cells[ind].virions > 0.0f);
			cells[ind].virions = params.initialVirions;
		}

		// Infecting a central cell of a rectangle tissue
		// ind = numCells/2 + 149; // 300 x 300 square tissue
		// ind = 27332; // lung windows selection
		// cells[ind].state = INFECTED_PLUS;
		// cells[ind].virions = params.initialVirions;

		/*==========================================*/
		// Simulation Loop
		/*==========================================*/
		int measureIdx = 0;
		double aucVirus = 0.0, aucIFN = 0.0;
		double prevVirus = 0.0, prevIFN = 0.0;
		int progressInterval = params.timeSteps / 10;
		if (progressInterval == 0) progressInterval = 1;
		for (int step=0; step<=params.timeSteps; step++)
		{
			if (step % progressInterval == 0) printf("."); fflush(stdout);

			if (step % params.measureInterval == 0)
			{
				tissue_metrics(cells, numCells, cellCounts, tissueSum, tissueSqSum, measureIdx++,
							   &aucVirus, &aucIFN, &prevVirus, &prevIFN, params.measureInterval);

				if (fCell) fprintf(fCell, "%d,%f,%f\n", step, cells[ind].virions, cells[ind].IFN);
			}

			// Print snapshots only if 1 replicate and printSnap is ON
			if (fSnap && step % params.snapInterval == 0)
			{
				if (params.printSnap == 2) print_infectedSnapshots(cells, numCells, fSnap,
												cellCounts[INFECTED_PLUS] + cellCounts[INFECTED_MINUS]);
				else print_tissueSnapshots(cells, numCells, fSnap);
			}

			// Generate GPU random numbers (3 non-overlapping draws per cell)
			curandGenerateUniform(gen, d_ranUni, 3*numCells);

			tissue_update<<<blks, ths>>>(cells, numCells, cellCounts, &params, d_ranUni);
				
			copy_fields<<<blks, ths>>>(cells, virions_old, IFN_old, numCells);
			tissue_diffusion<<<blks, ths>>>(cells, virions_old, IFN_old, numCells, &params);

			cudaDeviceSynchronize();
		}

		if (fCell) fclose(fCell);
		if (fSnap) fclose(fSnap);

		if (fAUC) fprintf(fAUC, "%d,%e,%e\n", rep+1, aucVirus, aucIFN);
	}

	if (fAUC) fclose(fAUC);

	/*==========================================*/
	// Finalize results: Average and Std Dev
	/*==========================================*/

	if (params.numReplicates == 1)
	{
		FILE *fTissue = fopen(path_tissueState, "w");
		if (fTissue)
		{
			fprintf(fTissue, "Time,Virus,IFN,Susceptible,Refractory,Infected,Dead,nonPermissive\n");
			for (int s=0; s<numMeasures; s++)
			{
				fprintf(fTissue, "%d", s*params.measureInterval);
				int tBase = s*7;
				for (int m=0; m<7; m++)
					fprintf(fTissue, ",%e", tissueSum[tBase+m]);
				fprintf(fTissue, "\n");
			}
			fclose(fTissue);
		}
	}
	else
	{
		FILE *fTavg = fopen(path_tissueAvg, "w");
		FILE *fTstd = fopen(path_tissueStd, "w");
		if (fTavg && fTstd)
		{
			fprintf(fTavg, "Time,Virus,IFN,Susceptible,Refractory,Infected,Dead,nonPermissive\n");
			fprintf(fTstd, "Time,Virus,IFN,Susceptible,Refractory,Infected,Dead,nonPermissive\n");

			for (int s=0; s<numMeasures; s++)
			{
				fprintf(fTavg, "%d", s*params.measureInterval);
				fprintf(fTstd, "%d", s*params.measureInterval);
				int tBase = s*7;
				for (int m=0; m<7; m++)
				{
					double avg = tissueSum[tBase+m] / params.numReplicates;
					double var = (tissueSqSum[tBase+m] / params.numReplicates) - (avg*avg);
					double sdev = sqrt(fmax(0.0, var));
					fprintf(fTavg, ",%e", avg);
					fprintf(fTstd, ",%e", sdev);
				}
				fprintf(fTavg, "\n");
				fprintf(fTstd, "\n");
			}
			fclose(fTavg);
			fclose(fTstd);
		}
	}

	printf("\nCompleted\n");

	// Clean up
	cudaFree(cells);
	cudaFree(virions_old);
	cudaFree(IFN_old);
	cudaFree(d_ranUni);
	cudaFree(cellCounts);
	curandDestroyGenerator(gen);
	free(tissueSum);
	free(tissueSqSum);

	return 0;
}
