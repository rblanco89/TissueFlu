#ifndef TISSUE_H
#define TISSUE_H

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
	CellState state;
	float virionCount;
}
Cell;

typedef struct
{
	int width, height;
	Cell *cells; // 1D array: cells[y*witdth + x]
}
Tissue;

// Tissue API
//Tissue* tissue_create(int width, int height);
//void tissue_initialize(Tissue *tissue);
//void tissue_advance(Tissue *tissue, int step);
//void tissue_free(Tissue *tissue);

#ifdef __cplusplus
extern "C" {
#endif

// Tissue API:
Tissue* tissue_create(int width, int height);
void tissue_initialize(Tissue *tissue);
void tissue_free(Tissue *tissue);

#ifdef __cplusplus
}
#endif

#endif
