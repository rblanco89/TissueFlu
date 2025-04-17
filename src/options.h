#ifndef OPTIONS_H
#define OPTIONS_H

typedef struct
{
    int timeSteps;
    int incubationPeriod;
    int expressingPeriod;
}
SimOptions;

extern SimOptions options;

void parse_options(const char *filename);

#endif

