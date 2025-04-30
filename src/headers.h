#define THS_MAX 256
#define MAX_NEIGHBORS 32

typedef enum
{
	HEALTHY,
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

//#ifdef __cplusplus
//extern "C" {
//#endif

void parse_options(const char *filename);
__host__ long nextPow2(long x);
__host__ void tissue_snapshots(Cell *cells, FILE *fSnap, int numCells);
__host__ void build_neighbors(Cell *cells, int numCells, float cutoff);
__global__ void tissue_update(Cell *cells, int numCells, float virionProduction);
__host__ void tissue_infection(Cell *cells, int numCells);

//#ifdef __cplusplus
//}
//#endif
