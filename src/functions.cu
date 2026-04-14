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

__host__  float stabilityCondition(Cell *cells, int numCells)
{
	float maxWeightSum = 0.0f;
	for (int i = 0; i < numCells; i++)
	{
		float ws = 0.0f;
		for (int j = 0; j < cells[i].numNeighbors; j++) ws += cells[i].weights[j];
		maxWeightSum = fmaxf(maxWeightSum, ws);
	}

	return 1.0f/maxWeightSum;
}

// ==================================================================

__host__ void tissue_metrics(Cell *cells, int numCells, int *cellCounts,
                             double *tissueSum, double *tissueSqSum, int measureIdx)
{
	double val[7];

	double viralLoad = 0.0;
	double IFNlevel = 0.0;
	for (int i=0; i<numCells; i++)
	{
		viralLoad += cells[i].virions;
		IFNlevel += cells[i].IFN;
	}

	val[0] = viralLoad;
	val[1] = IFNlevel;
	val[2] = (double)cellCounts[SUSCEPTIBLE];
	val[3] = (double)cellCounts[REFRACTORY];
	val[4] = (double)(cellCounts[INFECTED_PLUS] + cellCounts[INFECTED_MINUS]);
	val[5] = (double)cellCounts[DEAD];
	val[6] = (double)cellCounts[NONPERMISSIVE];

	int base = measureIdx * 7;
	for (int m = 0; m < 7; m++)
	{
		tissueSum[base + m]   += val[m];
		tissueSqSum[base + m] += val[m] * val[m];
	}
}

// ==================================================================

__host__ void print_tissueSnapshots(Cell *cells, int numCells, FILE *fSnap)
{
	// Define the grid and declare properties
	fprintf(fSnap, "%d\n", numCells);
	fprintf(fSnap, "Properties=species:I:1:pos:R:3:virions:R:1:IFN:R:1\n");

    for (int i=0; i<numCells; i++)
	{
		Cell cell = cells[i];

        // State as numeric for coloring in Ovito
        int state = cell.state;
		float3 r = cell.position;
		float virions = cell.virions;
		float IFN = cell.IFN;

        // Format: state radius x y
        fprintf(fSnap, "%d %f %f %f %f %f\n", state, r.x, r.y, r.z, virions, IFN);
	}
}

// ==================================================================

__host__ void print_infectedSnapshots(Cell *cells, int numCells, FILE *fSnap, int numInf)
{
	// Define the grid and declare properties
	fprintf(fSnap, "%d\n", numInf);
	fprintf(fSnap, "Properties=species:I:1:pos:R:3:virions:R:1:IFN:R:1\n");

	if (numInf == 0) return;

    for (int i=0; i<numCells; i++)
	{
		Cell cell = cells[i];

        // State as numeric for coloring in Ovito
        int state = cell.state;
		if (state != INFECTED_PLUS && state != INFECTED_MINUS) continue;
		float3 r = cell.position;
		float virions = cell.virions;
		float IFN = cell.IFN;

        // Format: state radius x y
        fprintf(fSnap, "%d %f %f %f %f %f\n", state, r.x, r.y, r.z, virions, IFN);
	}
}
// ==================================================================
// DEVICE FUNCTIONS
// ==================================================================

__global__ void build_neighbors(Cell *cells, int numCells, float cutoff)
{
	int ind = threadIdx.x + blockIdx.x*blockDim.x;
	if (ind >= numCells) return;

	float3 ri, rj, dr;
	float dist2, w;
	// float weightSum = 0.0f;
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
			w = 1.0f / dist2;
			// weightSum += w;
			cells[ind].neighbors[numNeighbors] = ind_j;
			cells[ind].weights[numNeighbors] = w;
			numNeighbors++;
		}
		else
		{
			printf("Warning: cell %d neighbor list full\n", ind);
			break;
		}
	}

	cells[ind].numNeighbors = numNeighbors;
	// for (int j=0; j<numNeighbors; j++)
		// cells[ind].weights[j] /= weightSum;
}

// ==================================================================

__global__ void compute_weights(Cell *cells, int numCells)
{
	int ind = threadIdx.x + blockIdx.x*blockDim.x;
	if (ind >= numCells) return;

	float3 ri = cells[ind].position;
	int numNeighbors = cells[ind].numNeighbors;

	for (int j = 0; j < numNeighbors; j++)
	{
		int ind_j = cells[ind].neighbors[j];
		float3 rj = cells[ind_j].position;

		float dx = ri.x - rj.x;
		float dy = ri.y - rj.y;
		float dz = ri.z - rj.z;
		float dist2 = dx*dx + dy*dy + dz*dz;

		if (dist2 == 0.0f)
		{
			printf("Warning: cell %d and neighbor %d are at the same position\n", ind, ind_j);
			cells[ind].weights[j] = 0.0f;
		}
		else
		{
			cells[ind].weights[j] = 1.0f / dist2;
		}
	}
}

// ==================================================================

__device__ float hillFun(float x, float K, float n)
{
	if (x <= 0.0f) return 0.0f;
    float xn = powf(x, n);
    float Kn = powf(K, n);
    return xn / (Kn + xn);
}

// ==================================================================

__device__ float virionProduction(int time)
{
	const float A  = 17703.6603f;    // max cumulative virions
    const float K  = 16.1275f * 60; // half-max time (min)
    const float n  = 2.3904;        // Hill coefficient

    float t  = (float)time;
    float Kn = powf(K, n);
    float tn = powf(t, n);

    // Derivative of Hill function: dH/dt
    return A * n * Kn * powf(t, n - 1.0f) / powf(Kn + tn, 2.0f);
}

// ==================================================================

__global__ void tissue_update(Cell *cells, int numCells, int *cellCounts,
							  float IFNcellProb, float *ranUni)
{
	int ind = threadIdx.x + blockIdx.x*blockDim.x;
	if (ind >= numCells) return;

	float pFmax = 0.00025f; // max IFN production rate (IFN min^-1 dsRNA^-1)
	float k_syn = 1.0f; // dsRNA synthesis rate (dsRNA min^-1 virions^-1)
	float k_deg = 0.15f / 60.0f; // dsRNA degradation rate (min^-1)

	float logVirions, virions, refracProb, infecProb, suppProb, effInfProb;
	Cell *cell = &cells[ind];
	switch (cell->state)
	{
		// case NONPERMISSIVE:
		// 	refracProb = hillFun(cell->IFN, 10.0f, 3.0f);
		// 	if (refracProb > ranUni[ind]) cell->state = REFRACTORY;
		// 	break;

		case SUSCEPTIBLE:
			refracProb = hillFun(cell->IFN, 10.0f, 3.0f); // refractory mechanism
			// refracProb = 0.0f;
			if (refracProb > ranUni[ind])
			{
				cell->state = REFRACTORY;
				atomicAdd(&cellCounts[SUSCEPTIBLE], -1);
				atomicAdd(&cellCounts[REFRACTORY],   1);
				break;
			}

			logVirions = log10(cell->virions + 1.0f);
			infecProb = hillFun(logVirions, 3.0f, 3.0f); // infection mechanism
			suppProb = 1.0f - hillFun(cell->IFN, 5.0f, 3.0f); // suppression mechanism
			effInfProb = infecProb * suppProb;
			// effInfProb = infecProb;
			// effInfProb = 0.0f;
			if (effInfProb > ranUni[(ind+1)%numCells])
			{
				if (ranUni[(ind+2)%numCells] < IFNcellProb)
				{
					cell->state = INFECTED_PLUS;
					atomicAdd(&cellCounts[SUSCEPTIBLE],    -1);
					atomicAdd(&cellCounts[INFECTED_PLUS],   1);
				}
				else
				{
					cell->state = INFECTED_MINUS;
					atomicAdd(&cellCounts[SUSCEPTIBLE],    -1);
					atomicAdd(&cellCounts[INFECTED_MINUS],  1);
				}
			}
			break;

		case INFECTED_PLUS:
			cell->infectingTime--;
			virions = virionProduction(cell->internalTime++);
			cell->virions += virions; // virus field (virions * dt)

			cell->dsRNA += k_syn * virions - k_deg * cell->dsRNA;
			cell->IFN += pFmax * cell->dsRNA;
			if (cell->infectingTime <= 0)
			{
				cell->state = DEAD;
				atomicAdd(&cellCounts[INFECTED_PLUS], -1);
				atomicAdd(&cellCounts[DEAD],           1);
			}
			break;

		case INFECTED_MINUS:
			cell->infectingTime--;
			cell->virions += virionProduction(cell->internalTime++);
			// No dsRNA, no IFN production for minus-strand infected cells
			if (cell->infectingTime <= 0)
			{
				cell->state = DEAD;
				atomicAdd(&cellCounts[INFECTED_MINUS], -1);
				atomicAdd(&cellCounts[DEAD],            1);
			}
			break;

		default:
			break;
	}
}

// ==================================================================

__global__ void tissue_diffusion(Cell *cells, int numCells,
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
		float w = cell.weights[j];

		diffVirions += w*(cells[ind_j].virions - cell.virions);
		diffIFN += w*(cells[ind_j].IFN - cell.IFN);
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
