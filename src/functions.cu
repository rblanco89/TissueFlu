/* Functions for create, initialize, and clean the tissue */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "headers.h"

// ==================================================================

void tissue_snapshots(Cell *cells, FILE *fSnap, int numCells)
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
        fprintf(fSnap, "%d %f %f %f %f\n",
						state, 0.5, r.x, r.y, r.z);
	}
}

// ==================================================================

void build_neighbors(Cell *cells, int numCells, float cutoff)
{
	for (int i = 0; i < numCells; i++)
	{
        cells[i].numNeighbors = 0;
    }

	for (int i=0; i<numCells; i++)
	{
        cells[i].numNeighbors = 0;

		for (int j=i+1; j<numCells; j++)
		{
			float3 ri = cells[i].position;
			float3 rj = cells[j].position;

            float dx = ri.x - rj.x;
            float dy = ri.y - rj.y;
            float dz = ri.z - rj.z;

            float dist2 = dx*dx + dy*dy + dz*dz;

            if (dist2 <= cutoff*cutoff)
			{
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
    }
}

// ==================================================================

//__host__ void tissue_advance(Tissue *tissue)
//{
//	int width = tissue->width;
//	int height = tissue->height;
//    //int totalCells = width*height;
//
//    for (int y=0; y<height; y++)
//	{
//        for (int x=0; x<width; x++)

//		{
//            int idx = y*width + x;
//            Cell *cell = &tissue->cells[idx];
//
//            if (cell->state == CELL_INCUBATING) 
//			{
//				cell->state = CELL_EXPRESSING;
//				continue;
//			}
//
//            if (cell->state == CELL_EXPRESSING)
//			{
//                // Produce virions, infect neighbors
//                cell->virionCount += 3.0f;
//
//                int dx[4] = { -1, 1, 0, 0 };
//                int dy[4] = { 0, 0, -1, 1 };
//
//                for (int d=0; d<4; d++)
//				{
//                    int nx = x + dx[d];
//                    int ny = y + dy[d];
//                    if (nx >= 0 && nx < width && ny >= 0 && ny < height)
//					{
//                        int nidx = ny*width + nx;
//                        Cell *neighbor = &tissue->cells[nidx];
//                        if (neighbor->state == CELL_SUSCEPTIBLE && cell->virionCount > 3.0f)
//                            neighbor->state = CELL_INCUBATING;
//                    }
//                }
//            }
//        }
//    }
//}
