#ifndef OPTIONS_H
#define OPTIONS_H

typedef struct
{
    int gridWidth;
    int gridHeight;
    int timeSteps;
}
SimOptions;

extern SimOptions options;

void parse_options(const char *filename);

#endif

