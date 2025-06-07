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
__host__ void print_tissueStatus(Cell *cells, int numCells, int step, FILE *fStat)
{
	float viralLoad = 0.0, IFNlevel = 0.0;
	int healthyCells = 0, refractoryCells = 0,
		infectedCells = 0, deadCells = 0;

	for (int i=0; i<numCells; i++)
	{
		viralLoad += cells[i].virions;
		IFNlevel += cells[i].IFN;
		switch (cells[i].state)
		{
			case HEALTHY:
				healthyCells++;
				break;

			case REFRACTORY:
				refractoryCells++;
				break;

			case INCUBATING:
				infectedCells++;
				break;

			case EXPRESSING:
				infectedCells++;
				break;

			case DEAD:
				deadCells++;
				break;

			default:
				break;
		}
	}

	fprintf(fStat, "%d,%e,%e,%d,%d,%d,%d\n", step, viralLoad, IFNlevel,
		healthyCells, refractoryCells, infectedCells, deadCells);
}
// ==================================================================

__host__ void print_tissueSnapshots(Cell *cells, int numCells, FILE *fSnap)
{
	// Define the grid and declare properties
	fprintf(fSnap, "%d\n", numCells);
	//fprintf(fSnap, "Lattice=\"%.2f 0.00 0.00 ", width);
	//fprintf(fSnap, "0.00 %.2f 0.00 ", height);
	//fprintf(fSnap, "0.00 0.00 %.2f\" ", 2.0);
	fprintf(fSnap, "Properties=species:I:1:pos:R:3\n");

    for (int i=0; i<numCells; i++)
	{
		Cell cell = cells[i];

        // State as numeric for coloring in Ovito
        int state = (int)cell.state;
		float3 r = cell.position;

        // Format: state radius x y
        fprintf(fSnap, "%d %f %f %f\n", state, r.x, r.y, r.z);
	}
}

// ==================================================================

//__host__ void build_neighbors(Cell *cells, int numCells, float cutoff)
//{
//	for (int i = 0; i < numCells; i++)
//		cells[i].numNeighbors = 0;
//
//	for (int i=0; i<numCells-1; i++)
//		for (int j=i+1; j<numCells; j++)
//		{
//			float3 ri = cells[i].position;
//			float3 rj = cells[j].position;
//
//			float dx = ri.x - rj.x;
//			float dy = ri.y - rj.y;
//			float dz = ri.z - rj.z;
//
//			float dist2 = dx*dx + dy*dy + dz*dz;
//
//			if (dist2 > cutoff*cutoff) continue;
//
//			// Add j to i's neighbor list
//			if (cells[i].numNeighbors < MAX_NEIGHBORS)
//				cells[i].neighbors[cells[i].numNeighbors++] = j;
//			else
//			{
//				fprintf(stderr, "Warning: cell %d neighbor list full\n", i);
//				break;
//			}
//
//			// Add i to j's neighbor list
//			if (cells[j].numNeighbors < MAX_NEIGHBORS)
//				cells[j].neighbors[cells[j].numNeighbors++] = i;
//			else
//			{
//				fprintf(stderr, "Warning: cell %d neighbor list full\n", j);
//				break;
//			}
//		}
//}

// ==================================================================
// DEVICE FUNCTIONS
// ==================================================================

// Function to calculate the virion release rate at time t (in min)
__device__ float virionRelease(int time)
{
	// Parameters for the viral production rate (logistic)
	float K = 15000.0;
	float r = 0.004;
	float t0 = 315.0; //(13.25 hrs)
	float t = time; 
	
	// Logistic curve derivative
	float exponent = exp(-r*(t - t0));

	return (K*r*exponent) / pow(1 + exponent, 2);
}

// ==================================================================

__device__ float sigmoidFun(float x, float A, float K)
{
	return 1/(1 + exp(-A*(x-K)));
}

// ==================================================================

__global__ void build_neighbors(Cell *cells, int numCells, float cutoff)
{
	int ind = threadIdx.x + blockIdx.x*blockDim.x;
	if (ind >= numCells) return;

	float3 ri, rj, dr;
	float dist2;
	int numNeighbors = 0;

	for (int ind_j=0; ind_j<numCells; ind_j++)
	{
		if (ind == ind_j) continue;

		ri = cells[ind].position;
		rj = cells[ind_j].position;

		dr.x = ri.x - rj.x;
		dr.y = ri.y - rj.y;
		dr.z = ri.z - rj.z;

		dist2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;

		if (dist2 > cutoff*cutoff) continue;
		if (dist2 == 0.0f)
		{
			printf("Warning: cell %d and cell %d are at the same position\n", ind, ind_j);
			break;
		}

		// Add j to i's neighbor list
		if (numNeighbors < MAX_NEIGHBORS)
		{
			cells[ind].neighbors[numNeighbors] = ind_j;
			cells[ind].neighDist2[numNeighbors] = dist2;
			numNeighbors++;
		}
		else
		{
			printf("Warning: cell %d neighbor list full\n", ind);
			break;
		}
	}

	cells[ind].numNeighbors = numNeighbors;
}

// ==================================================================

__global__ void tissue_update(Cell *cells, int numCells,
							  float virionProduction, float IFNproduction,
							  float *ranUni)
{
	int ind = threadIdx.x + blockIdx.x*blockDim.x;
	if (ind >= numCells) return;

	short infecFlag = 0, refracFlag = 0;
	float virions, ifn;
	float infecProb, refracProb;
	Cell *cell = &cells[ind];
	switch (cell->state)
	{
		case HEALTHY:
			//infecProb = 0.001*cell->virions;
			//refracProb = 0.001*cell->IFN;

			virions = cell->virions <= 0.0f ? -38 : log10(cell->virions);
			ifn = cell->IFN <= 0.0f ? -38 : log10(cell->IFN);
			infecProb = sigmoidFun(virions, 2, 3);
			refracProb = sigmoidFun(ifn, 2, 3);

			//refracProb = 1 - sigmoidFun(ifn, 2, 3);
			if (infecProb*(1.0f-refracProb) > ranUni[ind]) infecFlag = 1;
			//if (infecProb > ranUni[ind]) infecFlag = 1;

			refracProb = sigmoidFun(ifn, 2, 4);
			if (refracProb > ranUni[ind]) refracFlag = 1;

			if (infecFlag && refracFlag)
			{
				if (infecProb < refracProb) infecFlag = 0;
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
			//cell->IFN += IFNproduction;
			if (cell->incubationTime <= 0)
				cell->state = EXPRESSING;
			break;

		case EXPRESSING:
			cell->expressingTime--;
			cell->virions += virionRelease(cell->internalTime++);
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

	// MODEL 1
	float diffVirions = 0.0;
	float diffIFN = 0.0;

	for (int j=0; j<cell.numNeighbors; j++)
	{
		int ind_j = cell.neighbors[j];
		float dist2 = cell.neighDist2[j];

		diffVirions += (cells[ind_j].virions - cell.virions) / dist2;
		diffIFN += (cells[ind_j].IFN - cell.IFN) / dist2;
	}

	// Update count of virions for each cell
	float diffusedVirions = virionDiffusion*diffVirions;
	cells[ind].virions = (1.0 - virionClearance)*(cell.virions + diffusedVirions);

	// Update IFN for each cell
	float diffusedIFN = IFNdiffusion*diffIFN;
	cells[ind].IFN = (1.0 - IFNclearance)*(cell.IFN + diffusedIFN);

	// MODEL 2 (New Mexico Approach)
	//float meanVirions = cell.virions;
	//float meanIFN = cell.IFN;

	//for (int j=0; j<cell.numNeighbors; j++)
	//{
	//	meanVirions += cells[cell.neighbors[j]].virions;
	//	meanIFN += cells[cell.neighbors[j]].IFN;
	//}

	//meanVirions /= cell.numNeighbors + 1;
	//meanIFN /= cell.numNeighbors + 1;

	//// Update count of virions for each cell
	//float diffusedVirions = virionDiffusion*(meanVirions - cell.virions);
	//cells[ind].virions = (1.0 - virionClearance)*(cell.virions + diffusedVirions);

	//// Update IFN for each cell
	//float diffusedIFN = IFNdiffusion*(meanIFN - cell.IFN);
	//cells[ind].IFN = (1.0 - IFNclearance)*(cell.IFN + diffusedIFN);
}
