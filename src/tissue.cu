/* TissueFlu - building the tissue
 *
 * Reads the cell positions, obtains the neighbour lists (from a cached file or
 * by searching), re-orders the cells for memory locality and uploads
 * everything to the GPU.  All of this happens once, before the simulation.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <vector>
#include "headers.h"

/*==========================================*/
/* Small helpers                            */
/*==========================================*/

#define CUDA_CHECK(call) do {                                                  \
	cudaError_t err_ = (call);                                                 \
	if (err_ != cudaSuccess) {                                                 \
		fprintf(stderr, "CUDA error at %s:%d: %s\n",                           \
		        __FILE__, __LINE__, cudaGetErrorString(err_));                 \
		return 0;                                                              \
	}                                                                          \
} while (0)

/* Read a whole file into memory.  The neighbour file can be tens of megabytes
   and is parsed on every run, so it is read in one call and scanned in place
   rather than line by line. */
static char *slurp(const char *path, size_t *size)
{
	FILE *f = fopen(path, "rb");
	if (!f) return NULL;

	fseek(f, 0, SEEK_END);
	long n = ftell(f);
	rewind(f);

	char *buf = (char*)malloc(n + 1);
	if (!buf) { fclose(f); return NULL; }

	size_t got = fread(buf, 1, n, f);
	fclose(f);
	buf[got] = '\0';
	*size = got;
	return buf;
}

static inline int parse_int(const char **s)
{
	const char *p = *s;
	while (*p && (*p < '0' || *p > '9') && *p != '-') p++;
	int sign = 1;
	if (*p == '-') { sign = -1; p++; }
	int v = 0;
	while (*p >= '0' && *p <= '9') v = v*10 + (*p++ - '0');
	*s = p;
	return sign*v;
}

/*==========================================*/
/* Cell positions                           */
/*==========================================*/

/* Structure file: one header line, then "x,y,z" per cell. */
static int read_positions(const char *path, std::vector<float3> &pos)
{
	size_t size;
	char *buf = slurp(path, &size);
	if (!buf)
	{
		fprintf(stderr, "Error: could not open structure file %s\n", path);
		return 0;
	}

	const char *p = buf;
	while (*p && *p != '\n') p++;      // skip the header line
	if (*p) p++;

	while (*p)
	{
		char *end;
		float3 r;
		r.x = strtof(p, &end); if (end == p) break; p = end; while (*p == ',' || *p == ' ') p++;
		r.y = strtof(p, &end); if (end == p) break; p = end; while (*p == ',' || *p == ' ') p++;
		r.z = strtof(p, &end); if (end == p) break; p = end;
		while (*p && *p != '\n') p++;
		if (*p) p++;
		pos.push_back(r);
	}

	free(buf);

	if (pos.empty())
	{
		fprintf(stderr, "Error: no cells found in %s\n", path);
		return 0;
	}
	return 1;
}

/*==========================================*/
/* Spatial ordering                         */
/*==========================================*/

/* Interleave the bits of x, y and z (a Morton, or Z-order, code).  Cells with
   nearby codes are nearby in space, so sorting by the code puts each cell's
   neighbours close to it in memory - which is what makes the diffusion kernel
   fast, since it reads every neighbour of every cell at every time step. */
static inline unsigned long long spread_bits(unsigned int v)
{
	unsigned long long x = v & 0x1fffff;              // keep 21 bits
	x = (x | x << 32) & 0x1f00000000ffffULL;
	x = (x | x << 16) & 0x1f0000ff0000ffULL;
	x = (x | x <<  8) & 0x100f00f00f00f00fULL;
	x = (x | x <<  4) & 0x10c30c30c30c30c3ULL;
	x = (x | x <<  2) & 0x1249249249249249ULL;
	return x;
}

/* slot[o] = position in simulation order of the cell on line o of the file. */
static void spatial_order(const std::vector<float3> &pos, float binSize,
                          std::vector<int> &cellId, std::vector<int> &slotOf)
{
	int numCells = (int)pos.size();

	float3 lo = pos[0];
	for (int i = 1; i < numCells; i++)
	{
		lo.x = fminf(lo.x, pos[i].x);
		lo.y = fminf(lo.y, pos[i].y);
		lo.z = fminf(lo.z, pos[i].z);
	}
	if (!(binSize > 0.0f)) binSize = 1.0f;

	std::vector<std::pair<unsigned long long,int> > key(numCells);
	for (int i = 0; i < numCells; i++)
	{
		unsigned int bx = (unsigned int)((pos[i].x - lo.x)/binSize);
		unsigned int by = (unsigned int)((pos[i].y - lo.y)/binSize);
		unsigned int bz = (unsigned int)((pos[i].z - lo.z)/binSize);
		key[i] = std::make_pair((spread_bits(bx) << 2) | (spread_bits(by) << 1)
		                                                | spread_bits(bz), i);
	}
	std::stable_sort(key.begin(), key.end());

	cellId.resize(numCells);
	slotOf.resize(numCells);
	for (int s = 0; s < numCells; s++)
	{
		cellId[s] = key[s].second;
		slotOf[key[s].second] = s;
	}
}

/*==========================================*/
/* Neighbour search                         */
/*==========================================*/

/* Every cell tests every other cell.  This runs once per structure (the result
   is cached on disk) and the positions stay in the GPU caches, so the simple
   all-pairs search is fast enough not to need a spatial index. */
__global__ void count_neighbors(const float3 *pos, int numCells, float cutoff2,
                                int *count)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	if (i >= numCells) return;

	float3 ri = pos[i];
	int n = 0;
	for (int j = 0; j < numCells; j++)
	{
		float dx = ri.x - pos[j].x, dy = ri.y - pos[j].y, dz = ri.z - pos[j].z;
		float d2 = dx*dx + dy*dy + dz*dz;
		if (d2 > 0.0f && d2 <= cutoff2) n++;
	}
	count[i] = n;
}

__global__ void fill_neighbors(const float3 *pos, int numCells, float cutoff2,
                               const int *start, int *neighbor)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	if (i >= numCells) return;

	float3 ri = pos[i];
	int at = start[i];
	for (int j = 0; j < numCells; j++)
	{
		float dx = ri.x - pos[j].x, dy = ri.y - pos[j].y, dz = ri.z - pos[j].z;
		float d2 = dx*dx + dy*dy + dz*dz;
		if (d2 > 0.0f && d2 <= cutoff2) neighbor[at++] = j;
	}
}

static int build_neighbors(const std::vector<float3> &pos, float radius,
                           std::vector<int> &start, std::vector<int> &neighbor)
{
	int numCells = (int)pos.size();
	int blocks = (numCells + THREADS_PER_BLOCK - 1)/THREADS_PER_BLOCK;

	float3 *d_pos; int *d_count;
	CUDA_CHECK(cudaMalloc(&d_pos, numCells*sizeof(float3)));
	CUDA_CHECK(cudaMalloc(&d_count, numCells*sizeof(int)));
	CUDA_CHECK(cudaMemcpy(d_pos, &pos[0], numCells*sizeof(float3), cudaMemcpyHostToDevice));

	count_neighbors<<<blocks, THREADS_PER_BLOCK>>>(d_pos, numCells, radius*radius, d_count);
	CUDA_CHECK(cudaDeviceSynchronize());

	start.resize(numCells + 1);
	std::vector<int> count(numCells);
	CUDA_CHECK(cudaMemcpy(&count[0], d_count, numCells*sizeof(int), cudaMemcpyDeviceToHost));

	start[0] = 0;
	for (int i = 0; i < numCells; i++) start[i+1] = start[i] + count[i];
	neighbor.resize(start[numCells]);

	int *d_start, *d_neighbor;
	CUDA_CHECK(cudaMalloc(&d_start, (numCells + 1)*sizeof(int)));
	CUDA_CHECK(cudaMalloc(&d_neighbor, (start[numCells] ? start[numCells] : 1)*sizeof(int)));
	CUDA_CHECK(cudaMemcpy(d_start, &start[0], (numCells + 1)*sizeof(int), cudaMemcpyHostToDevice));

	fill_neighbors<<<blocks, THREADS_PER_BLOCK>>>(d_pos, numCells, radius*radius,
	                                              d_start, d_neighbor);
	CUDA_CHECK(cudaDeviceSynchronize());
	if (start[numCells])
		CUDA_CHECK(cudaMemcpy(&neighbor[0], d_neighbor, start[numCells]*sizeof(int),
		                      cudaMemcpyDeviceToHost));

	cudaFree(d_pos); cudaFree(d_count); cudaFree(d_start); cudaFree(d_neighbor);
	return 1;
}

/*==========================================*/
/* Neighbour list cache file                */
/*==========================================*/

/* Format, one line per cell in the order of the structure file:
       <number of neighbours>,<index>,<index>,...                              */
static int read_neighbor_file(const char *path, int numCells,
                              std::vector<int> &start, std::vector<int> &neighbor)
{
	size_t size;
	char *buf = slurp(path, &size);
	if (!buf) return 0;

	start.assign(numCells + 1, 0);
	neighbor.clear();
	neighbor.reserve(size/6);

	const char *p = buf;
	for (int c = 0; c < numCells; c++)
	{
		while (*p == '\n' || *p == '\r') p++;
		if (!*p)
		{
			fprintf(stderr, "Error: neighbour file %s ends at cell %d\n", path, c);
			free(buf);
			return 0;
		}
		int n = parse_int(&p);
		for (int k = 0; k < n; k++) neighbor.push_back(parse_int(&p));
		start[c+1] = (int)neighbor.size();
		while (*p && *p != '\n') p++;
	}

	free(buf);
	return 1;
}

static void write_neighbor_file(const char *path, int numCells,
                                const std::vector<int> &start,
                                const std::vector<int> &neighbor,
                                const std::vector<int> &cellId,
                                const std::vector<int> &slotOf)
{
	FILE *f = fopen(path, "w");
	if (!f)
	{
		fprintf(stderr, "Warning: could not save neighbour list to %s\n", path);
		return;
	}

	// Written in the order of the structure file, so the cache does not depend
	// on the internal ordering used by the simulation.
	for (int o = 0; o < numCells; o++)
	{
		int s = slotOf[o];
		fprintf(f, "%d", start[s+1] - start[s]);
		for (int k = start[s]; k < start[s+1]; k++)
			fprintf(f, ",%d", cellId[neighbor[k]]);
		fputc('\n', f);
	}
	fclose(f);
}

/*==========================================*/
/* Assembling the tissue                    */
/*==========================================*/

/* Path of the cached neighbour list: <structure>_neighbors_r<R>.<ext> */
static void neighbor_path(const char *structureFile, int radius,
                          char *out, size_t outSize)
{
	const char *dot = strrchr(structureFile, '.');
	if (dot && dot != structureFile)
		snprintf(out, outSize, "%.*s_neighbors_r%d%s",
		         (int)(dot - structureFile), structureFile, radius, dot);
	else
		snprintf(out, outSize, "%s_neighbors_r%d", structureFile, radius);
}

static void report_tissue(const std::vector<float3> &pos,
                          const std::vector<int> &start,
                          const std::vector<Neighbor> &link, const Params *pars)
{
	int numCells = (int)pos.size();
	int minN = start[1] - start[0], maxN = minN;
	float maxWeightSum = 0.0f;

	for (int i = 0; i < numCells; i++)
	{
		int n = start[i+1] - start[i];
		minN = std::min(minN, n);
		maxN = std::max(maxN, n);

		float ws = 0.0f;
		for (int k = start[i]; k < start[i+1]; k++) ws += link[k].weight;
		maxWeightSum = fmaxf(maxWeightSum, ws);
	}

	printf("Neighbours per cell: min = %d, max = %d, avg = %.1f\n",
	       minN, maxN, (float)link.size()/numCells);

	// Explicit diffusion is stable while D * sum(weights) <= 1 for every cell
	if (maxWeightSum > 0.0f)
	{
		float limit = 1.0f/maxWeightSum;
		float used = fmaxf(pars->virionDiffusion, pars->IFNdiffusion);
		printf("Diffusion stability limit: D <= %.4g (largest D in use: %.4g)%s\n",
		       limit, used, (used > limit) ? "  <-- UNSTABLE" : "");
	}
}

int tissue_create(const char *structureFile, const Params *pars,
                  Tissue *t, TissueHost *h)
{
	/* ---- positions ---- */
	std::vector<float3> filePos;
	if (!read_positions(structureFile, filePos)) return 0;

	int numCells = (int)filePos.size();
	printf("Loaded %d cells from %s\n", numCells, structureFile);

	/* ---- re-order the cells for memory locality ---- */
	std::vector<int> cellId, slotOf;
	spatial_order(filePos, pars->neighRadius, cellId, slotOf);

	std::vector<float3> pos(numCells);
	for (int s = 0; s < numCells; s++) pos[s] = filePos[cellId[s]];

	/* ---- neighbour lists, in simulation order ---- */
	char cachePath[512];
	neighbor_path(structureFile, (int)pars->neighRadius, cachePath, sizeof(cachePath));

	std::vector<int> start, neighbor;
	std::vector<int> fileStart, fileNeighbor;

	if (read_neighbor_file(cachePath, numCells, fileStart, fileNeighbor))
	{
		printf("Loaded neighbour list from %s\n", cachePath);

		start.resize(numCells + 1);
		neighbor.resize(fileNeighbor.size());
		start[0] = 0;
		for (int s = 0; s < numCells; s++)
		{
			int o = cellId[s];
			int at = start[s];
			for (int k = fileStart[o]; k < fileStart[o+1]; k++)
			{
				int j = fileNeighbor[k];
				if (j < 0 || j >= numCells)
				{
					fprintf(stderr, "Error: %s refers to cell %d, outside the "
					        "structure file. Delete it and it will be rebuilt.\n",
					        cachePath, j);
					return 0;
				}
				neighbor[at++] = slotOf[j];
			}
			start[s+1] = at;
		}
	}
	else
	{
		printf("Creating list of neighbours (radius %.3g)...\n", pars->neighRadius);
		if (!build_neighbors(pos, pars->neighRadius, start, neighbor)) return 0;

		printf("Saving neighbour list to %s\n", cachePath);
		write_neighbor_file(cachePath, numCells, start, neighbor, cellId, slotOf);
	}

	/* Sorting each cell's list makes the run independent of where the list came
	   from, and makes the neighbour reads more cache-friendly. */
	for (int s = 0; s < numCells; s++)
		std::sort(neighbor.begin() + start[s], neighbor.begin() + start[s+1]);

	/* ---- coupling weights: 1/d^2 ---- */
	std::vector<Neighbor> link(neighbor.size());
	long coincident = 0;
	for (int s = 0; s < numCells; s++)
	{
		float3 ri = pos[s];
		for (int k = start[s]; k < start[s+1]; k++)
		{
			int j = neighbor[k];
			float dx = ri.x - pos[j].x, dy = ri.y - pos[j].y, dz = ri.z - pos[j].z;
			float d2 = dx*dx + dy*dy + dz*dz;

			link[k].index  = j;
			link[k].weight = (d2 > 0.0f) ? 1.0f/d2 : 0.0f;
			if (!(d2 > 0.0f)) coincident++;
		}
	}
	if (coincident)
		printf("Warning: %ld neighbour pairs sit at the same position "
		       "and were left uncoupled\n", coincident);

	report_tissue(pos, start, link, pars);

	/* ---- upload ---- */
	int numLinks = (int)link.size();
	t->numCells = numCells;
	t->numLinks = numLinks;

	CUDA_CHECK(cudaMalloc(&t->field,         numCells*sizeof(float2)));
	CUDA_CHECK(cudaMalloc(&t->fieldNext,     numCells*sizeof(float2)));
	CUDA_CHECK(cudaMalloc(&t->state,         numCells*sizeof(CellState)));
	CUDA_CHECK(cudaMalloc(&t->infectingTime, numCells*sizeof(int)));
	CUDA_CHECK(cudaMalloc(&t->internalTime,  numCells*sizeof(int)));
	CUDA_CHECK(cudaMalloc(&t->dsRNA,         numCells*sizeof(float)));
	CUDA_CHECK(cudaMalloc(&t->cellId,        numCells*sizeof(int)));
	CUDA_CHECK(cudaMalloc(&t->position,      numCells*sizeof(float3)));
	CUDA_CHECK(cudaMalloc(&t->linkStart,     (numCells + 1)*sizeof(int)));
	CUDA_CHECK(cudaMalloc(&t->link,          (numLinks ? numLinks : 1)*sizeof(Neighbor)));

	CUDA_CHECK(cudaMemcpy(t->cellId, &cellId[0], numCells*sizeof(int), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpy(t->position, &pos[0], numCells*sizeof(float3), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpy(t->linkStart, &start[0], (numCells + 1)*sizeof(int), cudaMemcpyHostToDevice));
	if (numLinks)
		CUDA_CHECK(cudaMemcpy(t->link, &link[0], numLinks*sizeof(Neighbor), cudaMemcpyHostToDevice));

	/* ---- host mirror, used for the initial condition and for output ---- */
	h->position = (float3*)malloc(numCells*sizeof(float3));
	h->cellId   = (int*)malloc(numCells*sizeof(int));
	h->slotOf   = (int*)malloc(numCells*sizeof(int));
	memcpy(h->position, &pos[0], numCells*sizeof(float3));
	memcpy(h->cellId, &cellId[0], numCells*sizeof(int));
	memcpy(h->slotOf, &slotOf[0], numCells*sizeof(int));
	CUDA_CHECK(cudaMallocHost(&h->field, numCells*sizeof(float2)));
	CUDA_CHECK(cudaMallocHost(&h->state, numCells*sizeof(CellState)));

	return 1;
}

void tissue_destroy(Tissue *t, TissueHost *h)
{
	cudaFree(t->field);         cudaFree(t->fieldNext);
	cudaFree(t->state);         cudaFree(t->infectingTime);
	cudaFree(t->internalTime);  cudaFree(t->dsRNA);
	cudaFree(t->cellId);        cudaFree(t->position);
	cudaFree(t->linkStart);     cudaFree(t->link);

	free(h->position); free(h->cellId); free(h->slotOf);
	cudaFreeHost(h->field); cudaFreeHost(h->state);
}
