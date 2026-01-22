#define THS_MAX 256
#define MAX_NEIGHBORS 64

typedef struct
{
    int timeSteps;
    int numInfections;
    int incubationPeriod;
    int expressingPeriod;
	int ranSeed;

	float neighRadius;
	float intrinRefracProb;
	float initialVirions;
	float virionDiffusion;
	float virionClearance;
	float IFNproduction;
	float IFNdiffusion;
	float IFNclearance;
}
Options;

extern Options options;

typedef enum
{
	HEALTHY,
	REFRACTORY,
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
	int internalTime;
	float virions;
	float IFN;

	int numNeighbors;
	int neighbors[MAX_NEIGHBORS];
	float neighDist2[MAX_NEIGHBORS];
}
Cell;

/*==========================================*/
// Functions
/*==========================================*/

void parse_options(const char *filename);
__host__ long nextPow2(long x);
__host__ void print_tissueSnapshots(Cell *cells, int numCells, FILE *fSnap);
__host__ void print_tissueStatus(Cell *cells, int numCells, int step, FILE *fStat);

//__host__ void build_neighbors(Cell *cells, int numCells, float cutoff);
//__host__ void tissue_update(Cell *cells, int numCells);
//__host__ void tissue_infection(Cell *cells, int numCells);

__global__ void build_neighbors(Cell *cells, int numCells, float cutoff);
__global__ void tissue_update(Cell *cells, int numCells, 
		float IFNproduction, float *d_ranUni);
__global__ void tissue_infection(Cell *cells, int numCells, float virionDiffusion,
		float virionClearance, float IFNdiffusion, float IFNclearance);
