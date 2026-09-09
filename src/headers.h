#define THS_MAX 256
#define MAX_NEIGHBORS 512

typedef struct
{
    int timeSteps;
    int numInfections;
	int numReplicates;
    int infectingPeriod;
	int ranSeed;
	int printSnap;
	int snapInterval;
	int measureInterval;
	int printReplicates;

	int flagRefrac;
	int flagSupp;
	int flagBP;
	int flagPF;
	int flagPorousDiff;
	int flagCTL;

	float neighRadius;
	float nonPermProb;
	float initialVirions;
	float virionDiffusion;
	float virionClearance;
	float IFNcellProb;
	float IFNdiffusion;
	float IFNclearance;

	float pFmax;
	float k_syn;
	float k_deg;
	float K_r;
	float K_s;
	float K_v;
	float K_bp;
	float K_pf;
	float alpha_pf;
	float nHill;

	float rho_T;
	float delta_T;
	float K_ifn;
	float T0;
	float K_T;
	float T_max;

	float T_sys;
}
Params;

extern Params params;

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
	float weights[MAX_NEIGHBORS];
}
Cell;

/*==========================================*/
// Functions
/*==========================================*/

void parse_parameters(const char *filename);
__host__ long nextPow2(long x);
__host__ float stabilityCondition(Cell *cells, int numCells);
__host__ __device__ float hillFun(float x, float K, float n);
__host__ void tissue_metrics(Cell *cells, int numCells, int *cellCounts,
                             double *tissueSum, double *tissueSqSum, int measureIdx,
                             double *aucVirus, double *aucIFN,
                             double *prevVirus, double *prevIFN, int measureInterval,
                             FILE *fRep, float T_sys);
__host__ void print_tissueSnapshots(Cell *cells, int numCells, FILE *fSnap);
__host__ void print_infectedSnapshots(Cell *cells, int numCells, FILE *fSnap, int numInf);

__global__ void build_neighbors(Cell *cells, int numCells, float cutoff, int *overflowFlag);
__global__ void compute_weights(Cell *cells, int numCells);
__global__ void copy_fields(Cell *cells, float *virions_old, float *IFN_old, int numCells,
							Params *pars, float *d_IFNsum);
__global__ void tissue_update(Cell *cells, int numCells, int *cellCounts, Params *pars, float *d_ranUni);
__global__ void tissue_diffusion(Cell *cells, const float *virions_old, const float *IFN_old,
								int numCells, Params *pars);
__global__ void update_Tsys(Params *pars, const float *d_IFNsum, int numCells);
