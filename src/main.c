// AeroFlue: Influenza simulation
// Author: Rodolfo Blanco
// Date: April 2025

#include <stdio.h>
#include <stdlib.h>
#include <curand.h>

int main()
{
	// Load tissue data
	//Tissue tissue;
	
	// Main simulation loop
	printf("Starting simulation...\n");

	int step;
	int max_steps = 100;

	for (step=0; step<max_steps; step++)
	{
		// simulate virus and cell dynamic
		// gather statistics
		// write to file

	}

	printf("Simulation completed\n");

	exit (0);
}
