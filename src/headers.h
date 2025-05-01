#define THS_MAX 256
#define MAX_NEIGHBORS 32

typedef struct
{
    int timeSteps;
    int numInfections;
    int incubationPeriod;
    int expressingPeriod;
	int ranSeed;

	float neighRadius;
	float initialVirions;
	float virionProduction;
	float virionDiffusion;
	float virionClearance;
}
Options;

extern Options options;

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

/*==========================================*/
// Functions
/*==========================================*/

void parse_options(const char *filename);
__host__ long nextPow2(long x);
__host__ void tissue_snapshots(Cell *cells, FILE *fSnap, int numCells);
__host__ void build_neighbors(Cell *cells, int numCells, float cutoff);

//__host__ void tissue_update(Cell *cells, int numCells);
//__host__ void tissue_infection(Cell *cells, int numCells);

__global__ void tissue_update(Cell *cells, int numCells, float virionProduction, float *d_ranUni);
__global__ void tissue_infection(Cell *cells, int numCells, float virionDiffusion,
		float virionClearance);
