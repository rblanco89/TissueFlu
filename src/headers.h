#define THS_MAX 256
#define MAX_NEIGHBORS 64

typedef struct
{
    int timeSteps;
    int numInfections;
    int infectingPeriod;
	int ranSeed;
	int printSnap;
	int snapInterval;
	int measureInterval;

	float neighRadius;
	float nonPermProb;
	float initialVirions;
	float virionDiffusion;
	float virionClearance;
	float IFNcellProb;
	float IFNdiffusion;
	float IFNclearance;
}
Options;

extern Options options;

typedef enum
{
	NONPERMISSIVE,
	SUSCEPTIBLE,
	INFECTED_MINUS,
	INFECTED_PLUS,
	REFRACTORY,
	DEAD
}
CellState;

typedef struct
{
	float3 position;
	CellState state;

	int infectingTime;
	int internalTime;
	float virions;
	float dsRNA;
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
__host__ void sustainability_check(Cell *cells, int numCells); 

__global__ void build_neighbors(Cell *cells, int numCells, float cutoff);
__global__ void tissue_update(Cell *cells, int numCells, float IFNcellProb, float *d_ranUni);
__global__ void tissue_diffusion(Cell *cells, int numCells, float virionDiffusion,
		float virionClearance, float IFNdiffusion, float IFNclearance);
