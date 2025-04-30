/* Functions for create, initialize, and clean the tissue */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "headers.h"
#include "options.h"

// ==================================================================

__host__ long nextPow2(long x)
{
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return ++x;
}

// ==================================================================

__host__ void tissue_snapshots(Cell *cells, FILE *fSnap, int numCells)
{
	// Define the grid and declare properties
	fprintf(fSnap, "%d\n", numCells);
	//fprintf(fSnap, "Lattice=\"%.2f 0.00 0.00 ", width);
	//fprintf(fSnap, "0.00 %.2f 0.00 ", height);
	//fprintf(fSnap, "0.00 0.00 %.2f\" ", 2.0);
	fprintf(fSnap, "Properties=species:I:1:Radius:R:1:pos:R:3\n");

    for (int i=0; i<numCells; i++)
	{
		Cell cell = cells[i];

        // State as numeric for coloring in Ovito
        int state = (int)cell.state;
		float3 r = cell.position;

        // Format: state radius x y
        fprintf(fSnap, "%d %f %f %f %f\n", state, 0.5, r.x, r.y, r.z);
	}
}

// ==================================================================

__host__ void build_neighbors(Cell *cells, int numCells, float cutoff)
{
	for (int i = 0; i < numCells; i++)
		cells[i].numNeighbors = 0;

	for (int i=0; i<numCells-1; i++)
		for (int j=i+1; j<numCells; j++)
		{
			//if (i == j) continue;
			float3 ri = cells[i].position;
			float3 rj = cells[j].position;

			float dx = ri.x - rj.x;
			float dy = ri.y - rj.y;
			float dz = ri.z - rj.z;

			float dist2 = dx*dx + dy*dy + dz*dz;

			if (dist2 > cutoff*cutoff) continue;

			// Add j to i's neighbor list
			if (cells[i].numNeighbors < MAX_NEIGHBORS)
				cells[i].neighbors[cells[i].numNeighbors++] = j;
			else
			{
				fprintf(stderr, "Warning: cell %d neighbor list full\n", i);
				break;
			}

			// Add i to j's neighbor list
			if (cells[j].numNeighbors < MAX_NEIGHBORS)
				cells[j].neighbors[cells[j].numNeighbors++] = i;
			else
			{
				fprintf(stderr, "Warning: cell %d neighbor list full\n", j);
				break;
			}
		}
}

// ==================================================================

//__host__ void tissue_update(Cell *cells, int numCells)
//{
//	for (int i=0; i<numCells; i++)
//	{
//		Cell *cell = &cells[i];
//		switch (cell->state)
//		{
//			case HEALTHY:
//				if (cell->virions > 50)
//					cell->state = INCUBATING;
//				break;
//			case INCUBATING:
//				cell->incubationTime--;
//				if (cell->incubationTime <= 0)
//					cell->state = EXPRESSING;
//				break;
//			case EXPRESSING:
//				cell->expressingTime--;
//				cell->virions += options.virionProduction;
//				if (cell->expressingTime <= 0)
//					cell->state = DEAD;
//				break;
//			default:
//				break;
//		}
//	}
//}

__global__ void tissue_update(Cell *cells, int numCells, float virionProduction)
{
	int ind = threadIdx.x + blockIdx.x*blockDim.x;
	if (ind >= numCells) return;

	Cell *cell = &cells[ind];
	switch (cell->state)
	{
		case HEALTHY:
			if (cell->virions > 50)
				cell->state = INCUBATING;
			break;
		case INCUBATING:
			cell->incubationTime--;
			if (cell->incubationTime <= 0)
				cell->state = EXPRESSING;
			break;
		case EXPRESSING:
			cell->expressingTime--;
			cell->virions += virionProduction;
			if (cell->expressingTime <= 0)
				cell->state = DEAD;
			break;
		default:
			break;
	}
}

// ==================================================================

__host__ void tissue_infection(Cell *cells, int numCells)
{
	// Infect susceptible neighbors
	for (int i=0; i<numCells; i++)
	{
		float meanVirions = cells[i].virions;
		for (int j=0; j<cells[i].numNeighbors; j++)
			meanVirions += cells[cells[i].neighbors[j]].virions;
		meanVirions /= cells[i].numNeighbors + 1;

		// Update count of virions for each cell
		float diffusedVirions = options.virionDiffusion*(meanVirions - cells[i].virions);
		float updatedVirions = cells[i].virions + diffusedVirions;
		cells[i].virions = (1 - options.virionClearance)*updatedVirions;
	}
}
