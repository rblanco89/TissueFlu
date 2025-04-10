/* AeroFlue: Influenza simulation
   Author: Rodolfo Blanco
   Date: April 2025 */

#include <stdio.h>
#include <stdlib.h>
#include <curand.h>

#include "options.h"
#include "tissue.h"

int main(int argc, char *argv[])
{
	const char *config_file = (argc > 1) ? argv[1] : "config.conf";
	parse_options(config_file);

	int gridWidth = options.gridWidth;
	int gridHeight = options.gridHeight;
	printf("Grid: %dx%d\n", gridWidth, gridHeight);
	
	// Load tissue
	Tissue *tissue = tissue_create(gridWidth, gridHeight);
	if (tissue == NULL)
	{
		fprintf(stderr, "Error: copuld not create tissue/n");
		return 1;
	}
	
	// Main simulation loop
	printf("Starting simulation...\n");

	//int step;
	//int max_steps = 100;

	//for (step=0; step<max_steps; step++)
	//{
		// simulate virus and cell dynamic
		// gather statistics
		// write to file

	//}

	printf("Simulation completed\n");

	return 0;
}
