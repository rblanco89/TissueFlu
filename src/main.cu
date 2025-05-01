/* AeroFlue: Influenza simulation
   Author: Rodolfo Blanco
   Date: April 2025 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cuda_runtime.h>
#include <curand.h>

#include "headers.h"
#include "ranNumbers.h"

int main(int argc, char *argv[])
{
	const char *config_file = NULL;
	const char *structure_file = NULL;

	/*==========================================*/
	// Parse command-line arguments
	/*==========================================*/

	for (int i = 1; i < argc; i++)
	{
		if (strncmp(argv[i], "--config=", 9) == 0)
		{
			config_file = argv[i] + 9;
		}
		else if (strncmp(argv[i], "--structure=", 12) == 0)
		{
			structure_file = argv[i] + 12;
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
				"Usage: %s --config=FILE --structure=FILE\n"
				"  --config     : path to configuration file\n"
				"  --structure  : path to cell positions file (e.g., CSV)\n",
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
	for (int i = 0; i < numCells; i++)
	{
		float3 r;
		if (fgets(line, sizeof(line), fp) == NULL)
		{
			fprintf(stderr, "Unexpected end of structure file at cell %d\n", i);
			break;
		}

		if (sscanf(line, "%f,%f,%f", &r.x, &r.y, &r.z) != 3)
		{
			fprintf(stderr, "Invalid format in structure file at line %d\n", i + 1);
			break;
		}

		cells[i].position = r;
		cells[i].state = HEALTHY;
		cells[i].virions = 0.0f;
		cells[i].incubationTime = options.incubationPeriod;
		cells[i].expressingTime = options.expressingPeriod;
	}

	fclose(fp);

	printf("Loaded %d cells from %s\n", numCells, structure_file);

	// Initialize random numbers
	ulong seed = options.ranSeed;
	Ran ranUni(seed);
	//Poissondev poisson(options.incubationPeriod, seed); // a = (aveTime/stdTime)^2; b = aveTime/stdTime^2

	//printf("%d\n", poisson.dev());
	
	// Initialize random numbers with a standard normal distribution
	float *d_ranUni;
	curandGenerator_t gen;
	cudaMalloc(&d_ranUni, numCells*sizeof(float)); // Array only for GPU
	curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_MTGP32);
	curandSetPseudoRandomGeneratorSeed(gen, seed);

	/*==========================================*/
	// Find neighbors and store in cell structure
	/*==========================================*/

	build_neighbors(cells, numCells, options.neighRadius);

	/*==========================================*/
	// Initial infection
	/*==========================================*/

	// For now I'm infecting a central cell of a rectangle tissue
	cells[315].state = INCUBATING;
	cells[315].virions = options.initialVirions;

	/*==========================================*/
	// Initialize files for results
	/*==========================================*/

	char filename[128];

	sprintf(filename, "results/snapshots.xyz");
	FILE *fSnap = fopen(filename, "w");
	if (!fSnap)
	{
		fprintf(stderr, "Error: could not open %s for writing\n", filename);
		return 1;
	}

	sprintf(filename, "results/viralLoad.csv");
	FILE *fVirions = fopen(filename, "w");
	if (!fVirions)
	{
		fprintf(stderr, "Error: could not open %s for writing\n", filename);
		return 1;
	}
	fprintf(fVirions, "Days,Viral Load\n");

	/*==========================================*/
	// Main loop: Simulation
	/*==========================================*/

	int ths, blks;

	// Estimate the number of threads and blocks for the GPU
	ths = (numCells < THS_MAX) ? nextPow2(numCells) : THS_MAX;
	blks = 1 + (numCells - 1)/ths;

	float viralLoad;

	short noPrintFlag;
	int printStep = 0.01*options.timeSteps;

	printf("Starting simulation...\n");

	for (int step=0; step<options.timeSteps; step++)
	{
		noPrintFlag = step%printStep;
		if (!noPrintFlag) printf("Step %d/%d\n", step, options.timeSteps);

		viralLoad = 0.0;
		for (int i=0; i<numCells; i++) viralLoad += cells[i].virions;
		fprintf(fVirions, "%d,%f\n", step, viralLoad);

		// GPU functions (kernels)

		// Generate random numbers and then update positions
		curandGenerateUniform(gen, d_ranUni, numCells);

		tissue_update<<<blks, ths>>>(cells, numCells, options.virionProduction, d_ranUni);
		tissue_infection<<<blks, ths>>>(cells, numCells, options.virionDiffusion,
								  options.virionClearance);
		cudaDeviceSynchronize();

		// Host functions
		//tissue_update(cells, numCells);
		//tissue_infection(cells, numCells);

		if (!noPrintFlag) tissue_snapshots(cells, fSnap, numCells);
	}

	fclose(fVirions);
	fclose(fSnap);

	printf("Simulation completed\n");


	// Clean up
	cudaFree(cells);

	return 0;
}
