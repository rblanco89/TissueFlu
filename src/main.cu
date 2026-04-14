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
	const char *output_dir = "results";

	/*==========================================*/
	// Parse command-line arguments
	/*==========================================*/

	// for (int i = 1; i < argc; i++)
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
				"  --output     : directory where results will be written (default: results)\n",
				argv[0]);
		return 1;
	}

	/*==========================================*/
	// Fetch structure and parameters
	/*==========================================*/

	parse_options(config_file);

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
	cudaMalloc(&d_ranUni, numCells*sizeof(float)); // Array only for GPU
	curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_MTGP32);

	// Estimate the number of threads and blocks for the GPU
	int ths = (numCells < THS_MAX) ? nextPow2(numCells) : THS_MAX;
	int blks = 1 + (numCells - 1)/ths;

	printf("Creating list of neighbors...\n");

	// Find neighbors and store in cell structure
	build_neighbors<<<blks, ths>>>(cells, numCells, options.neighRadius);
	cudaDeviceSynchronize();

	float maxDiffusion = stabilityCondition(cells, numCells);
	if (maxDiffusion < options.virionDiffusion || maxDiffusion < options.IFNdiffusion)
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
	char path_tissueAvg[512], path_tissueStd[512];
	snprintf(path_snapshots, sizeof(path_snapshots), "%s/snapshots.xyz",       output_dir);
	snprintf(path_cellState, sizeof(path_cellState), "%s/cellState.csv",       output_dir);
	snprintf(path_tissueAvg, sizeof(path_tissueAvg), "%s/tissueState_avg.csv", output_dir);
	snprintf(path_tissueStd, sizeof(path_tissueStd), "%s/tissueState_std.csv", output_dir);

	// Prepare for result accumulation
	int numMeasures = options.timeSteps / options.measureInterval + 1;

	// Arrays for averaging tissue state (7 variables: V, IFN, S, R, I, D, NP)
	double *tissueSum = (double*)calloc(numMeasures * 7, sizeof(double));
	double *tissueSqSum = (double*)calloc(numMeasures * 7, sizeof(double));

	FILE *fSnap = NULL;
	if (options.numReplicates == 1 && options.printSnap)
	{
		fSnap = fopen(path_snapshots, "w");
		if (!fSnap) fprintf(stderr, "Warning: could not open %s for writing\n", path_snapshots);
	}

	FILE *fCell = NULL;
	if (options.numReplicates == 1)
	{
		fCell = fopen(path_cellState, "w");
		if (!fCell) fprintf(stderr, "Warning: could not open %s for writing\n", path_cellState);
		fprintf(fCell, "Time,ViralLoad,IFN\n");
	}

	/*==========================================*/
	// Replicates Loop
	/*==========================================*/

	for (int rep=0; rep<options.numReplicates; rep++)
	{
		printf("\nRunning replicate %d/%d\n", rep+1, options.numReplicates);

		// Initialize/Reset random numbers
		ulong seed = options.ranSeed + rep;
		Ran ranUni(seed);
		Poissondev ranInfecting(options.infectingPeriod/60, seed);
		curandSetPseudoRandomGeneratorSeed(gen, seed);

		// Reset cells state
		for (int i=0; i<numCells; i++)
		{
			if (ranUni.doub() < options.nonPermProb) cells[i].state = NONPERMISSIVE;
			else cells[i].state = SUSCEPTIBLE;

			cells[i].virions = 0.0f;
			cells[i].dsRNA = 0.0f;
			cells[i].IFN = 0.0f;
			cells[i].infectingTime = 60*ranInfecting.dev(); // Convert to minutes
			cells[i].internalTime = 0;
		}
	
		int ind = 0;
		for (int i=0; i<options.numInfections; i++)
		{
			do ind = numCells*ranUni.doub();
			while (cells[ind].virions > 0.0f);
			cells[ind].virions = options.initialVirions;
		}

		// Infecting a central cell of a rectangle tissue
		// ind = numCells/2 + 149; // 300 x 300 square tissue
		// ind = 27332; // lung windows selection
		// cells[ind].state = INFECTED_PLUS;
		// cells[ind].virions = options.initialVirions;

		/*==========================================*/
		// Simulation Loop
		/*==========================================*/

		int measureIdx = 0;
		double metrics[7];
		int progressInterval = options.timeSteps / 10;
		if (progressInterval == 0) progressInterval = 1;
		for (int step=0; step<=options.timeSteps; step++)
		{
			if (step % progressInterval == 0) printf("."); fflush(stdout);

			if (step % options.measureInterval == 0)
			{
				memset(metrics, 0, sizeof(metrics));
				tissue_metrics(cells, numCells, metrics);
				for (int m=0; m<7; m++)
				{
					tissueSum[measureIdx*7 + m] += metrics[m];
					tissueSqSum[measureIdx*7 + m] += metrics[m]*metrics[m];
				}
				measureIdx++;
				
				if (fCell) fprintf(fCell, "%d,%f,%f\n", step, cells[ind].virions, cells[ind].IFN);
			}

			// Print snapshots only if 1 replicate and printSnap is ON
			if (fSnap && step % options.snapInterval == 0)
			{
				if (options.printSnap == 2)
				{
					if (options.snapInterval != options.measureInterval)
					{
						memset(metrics, 0, sizeof(metrics));
						tissue_metrics(cells, numCells, metrics);
					}
					print_infectedSnapshots(cells, numCells, fSnap, int(metrics[4]));
				}
				else print_tissueSnapshots(cells, numCells, fSnap);
			}

			// Generate GPU random numbers
			curandGenerateUniform(gen, d_ranUni, numCells);

			tissue_update<<<blks, ths>>>(cells, numCells, options.IFNcellProb, d_ranUni);
			tissue_diffusion<<<blks, ths>>>(cells, numCells, options.virionDiffusion,
				options.virionClearance, options.IFNdiffusion, options.IFNclearance);

			cudaDeviceSynchronize();
		}

		if (fCell) fclose(fCell);
		if (fSnap) fclose(fSnap);
	}

	/*==========================================*/
	// Finalize results: Average and Std Dev
	/*==========================================*/

	FILE *fTavg = fopen(path_tissueAvg, "w");
	FILE *fTstd = fopen(path_tissueStd, "w");
	if (fTavg && fTstd)
	{
		fprintf(fTavg, "Time,Virus,IFN,Susceptible,Refractory,Infected,Dead,nonPermissive\n");
		fprintf(fTstd, "Time,Virus,IFN,Susceptible,Refractory,Infected,Dead,nonPermissive\n");

		for (int s=0; s<numMeasures; s++)
		{
			fprintf(fTavg, "%d", s*options.measureInterval);
			fprintf(fTstd, "%d", s*options.measureInterval);
			int tBase = s*7;
			for (int m=0; m<7; m++)
			{
				double avg = tissueSum[tBase+m] / options.numReplicates;
				double var = (tissueSqSum[tBase+m] / options.numReplicates) - (avg*avg);
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

	printf("\nCompleted\n");

	// Clean up
	cudaFree(cells);
	cudaFree(d_ranUni);
	curandDestroyGenerator(gen);
	free(tissueSum);
	free(tissueSqSum);

	return 0;
}
