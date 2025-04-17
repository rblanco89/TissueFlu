#define MAX_NEIGHBORS 32

typedef enum
{
	CELL_SUSCEPTIBLE,
	CELL_INCUBATING,
	CELL_EXPRESSING,
	CELL_DEAD
}
CellState;

typedef struct
{
	float3 position;
	CellState state;

	float incubationTime;
	float expressingTime;
	float virions;

	int numNeighbors;
	int neighbors[MAX_NEIGHBORS];
}
Cell;

#ifdef __cplusplus
extern "C" {
#endif

void tissue_snapshots(Cell *cells, FILE *fSnap, int numCells);
void build_neighbors(Cell *cells, int numCells, float cutoff);
//Tissue* tissue_create(int width, int height);
//void tissue_initialize(Tissue *tissue);
//void tissue_free(Tissue *tissue);
//void tissue_advance(Tissue *tissue);
//void tissue_snapshots(Tissue *tissue, FILE *fSnap);

#ifdef __cplusplus
}
#endif
