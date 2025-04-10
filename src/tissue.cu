/* Functions for create, initialize, and clean the tissue */
#include <stdio.h>
#include <stdlib.h>
#include "tissue.h"

// ========================================================================

Tissue* tissue_create(int width, int height)
{
	Tissue *tissue;
	cudaMallocManaged(&tissue, sizeof(Tissue));
	if (tissue == NULL)
	{
		fprintf(stderr, "cudaMalloc fail for Tissue struct\n");
		return NULL;
	}

	tissue->width = width;
	tissue->height = height;

	// Allocate unified memory for cells
	size_t totalCells = width*height;
	cudaMallocManaged(&tissue->cells, totalCells*sizeof(Cell));
	if (tissue->cells == NULL)
	{
		fprintf(stderr, "cudaMalloc fail for Cell struct/n");
		cudaFree(tissue);
		return NULL;
	}

	// Initialize cells on host
	for (int i=0; i<totalCells; i++)
	{
		tissue->cells[i].state = CELL_SUSCEPTIBLE;
		tissue->cells[i].virionCount = 0.0f;
	}

	return tissue;
}

// ========================================================================

void tissue_initialize(Tissue *tissue)
{
	int cx = tissue->width/2;
	int cy = tissue->height/2;
	int idx = cy*tissue->width + cx;

	tissue->cells[idx].state = CELL_INCUBATING;
	tissue->cells[idx].virionCount = 10.0f;
}

// ========================================================================

void tissue_free(Tissue *tissue)
{
	cudaFree(tissue->cells);
	cudaFree(tissue);
}
//void tissue_advance(Tissue *tissue, int step);
