#ifndef OPTIONS_H
#define OPTIONS_H

typedef struct
{
    int timeSteps;
    int numInfections;
    int incubationPeriod;
    int expressingPeriod;

	float neighRadius;
	float initialVirions;
	float virionProduction;
	float virionDiffusion;
	float virionClearance;
}
Options;

extern Options options;

#endif
