#define MAX_NEIGHBORS 32

typedef enum
{
	SUSCEPTIBLE,
	INCUBATING,
	EXPRESSING,
	DEAD
}
CellState;

typedef struct
{
	float3 position;
	CellState state;

	int incubationTime;
	int expressingTime;
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
void tissue_update(Cell *cells, int numCells);
void tissue_infection(Cell *cells, int numCells);
//Tissue* tissue_create(int width, int height);
//void tissue_initialize(Tissue *tissue);
//void tissue_free(Tissue *tissue);
//void tissue_advance(Tissue *tissue);
//void tissue_snapshots(Tissue *tissue, FILE *fSnap);

#ifdef __cplusplus
}
#endif
