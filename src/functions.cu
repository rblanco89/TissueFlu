/* Functions for create, initialize, and clean the tissue */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "headers.h"

// ==================================================================
// HOST FUNCTIONS
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
// DEVICE FUNCTIONS
// ==================================================================

__device__ float sigmoidFun(float x, float A, float K)
{
	return 1/(1 + exp(-A*(x-K)));
}

// ==================================================================

__global__ void tissue_update(Cell *cells, int numCells,
							  float virionProduction, float IFNproduction,
							  float *ranUni)
{
	int ind = threadIdx.x + blockIdx.x*blockDim.x;
	if (ind >= numCells) return;

	short infecFlag = 0, refracFlag = 0;
	float infecProb, refracProb;
	Cell *cell = &cells[ind];
	switch (cell->state)
	{
		case HEALTHY:
			//infecProb = 0.001*cell->virions;
			//refracProb = 0.001*cell->IFN;

			infecProb = sigmoidFun(log10(cell->virions), 2, 3);
			refracProb = sigmoidFun(log10(cell->IFN), 2, 4);

			if (infecProb > ranUni[ind]) infecFlag = 1;
			if (refracProb > ranUni[ind]) refracFlag = 1;

			if (infecFlag && refracFlag)
			{
				if (0.5 > ranUni[(ind+1)%ind]) infecFlag = 0;
				else refracFlag = 0;
			}

			if (infecFlag)
			{
				cell->state = INCUBATING;
			
				// Randomly protect a percentage of cell (fake refractory state)
				//if (0.2 > ranUni[(ind+1)%ind]) cell->state = REFRACTORY;
				//else cell->state = INCUBATING;

				break;
			}

			if (refracFlag) cell->state = REFRACTORY;
			break;

		case REFRACTORY:
			break;

		case INCUBATING:
			cell->incubationTime--;
			if (cell->incubationTime <= 0)
				cell->state = EXPRESSING;
			break;

		case EXPRESSING:
			cell->expressingTime--;
			cell->virions += virionProduction;
			cell->IFN += IFNproduction;
			if (cell->expressingTime <= 0)
				cell->state = DEAD;
			break;

		default:
			break;
	}
}

// ==================================================================

__global__ void tissue_infection(Cell *cells, int numCells,
								 float virionDiffusion, float virionClearance,
								 float IFNdiffusion, float IFNclearance)
{
	int ind = threadIdx.x + blockIdx.x*blockDim.x;
	if (ind >= numCells) return;

	Cell cell = cells[ind];

	float meanVirions = cell.virions;
	float meanIFN = cell.IFN;

	for (int j=0; j<cell.numNeighbors; j++)
	{
		meanVirions += cells[cell.neighbors[j]].virions;
		meanIFN += cells[cell.neighbors[j]].IFN;
	}

	meanVirions /= cell.numNeighbors + 1;
	meanIFN /= cell.numNeighbors + 1;

	// Update count of virions for each cell
	float diffusedVirions = virionDiffusion*(meanVirions - cell.virions);
	cells[ind].virions = (1.0 - virionClearance)*(cell.virions + diffusedVirions);


	// Update IFN for each cell
	float diffusedIFN= IFNdiffusion*(meanIFN - cell.IFN);
	cells[ind].IFN = (1.0 - IFNclearance)*(cell.IFN + diffusedIFN);
}
