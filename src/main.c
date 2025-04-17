/* AeroFlue: Influenza simulation
   Author: Rodolfo Blanco
   Date: April 2025 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cuda_runtime.h>
#include <curand.h>

#include "options.h"
#include "headers.h"

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
		cells[i].state = CELL_SUSCEPTIBLE;
		cells[i].virions = 0.0f;
		cells[i].incubationTime = options.incubationPeriod;
		cells[i].expressingTime = options.expressingPeriod;
	}

	fclose(fp);

	printf("Loaded %d cells from %s\n", numCells, structure_file);

	/*==========================================*/
	// Find neighbors and store in cell structure
	/*==========================================*/

	build_neighbors(cells, numCells, 1.5);

	/*==========================================*/
	// Initial infection
	/*==========================================*/

	cells[315].state = CELL_INCUBATING;

	/*==========================================*/
	// Initialize snapshot file
	/*==========================================*/

	char filename[128];
	sprintf(filename, "results/snapshots.xyz");
	FILE *fSnap = fopen(filename, "w");
	if (!fSnap)
	{
		fprintf(stderr, "Error: could not open %s for writing\n", filename);
		return 1;
	}

	/*==========================================*/
	// Main loop: Simulation
	/*==========================================*/

	printf("Starting simulation...\n");

	for (int step=0; step<options.timeSteps; step++)
	{
		printf("Step %d/%d\n", step, options.timeSteps);

		tissue_advance(cells, numCells, options.incubationPeriod, options.expressingPeriod);

		tissue_snapshots(cells, fSnap, numCells);

		cudaDeviceSynchronize();
	}
	fclose(fSnap);

	printf("Simulation completed\n");


	// Clean up
	cudaFree(cells);

	return 0;
}
