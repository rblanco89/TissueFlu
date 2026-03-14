#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "headers.h"

Options options;  // global struct

void parse_options(const char *filename)
{
	FILE *file = fopen(filename, "r");
	if (!file)
	{
        	perror("Could not open config file");
        	exit(1);
	}

	char line[256];
	while (fgets(line, sizeof(line), file))
	{
		char key[64], value[128];

		if (line[0] == '#' || line[0] == '\n')
			continue; // skip comments and blank lines

		if (sscanf(line, "%[^=]=%s", key, value) == 2)
		{
			if (strcmp(key, "timeSteps ") == 0)
				options.timeSteps = atoi(value);
			else if (strcmp(key, "numInfections ") == 0)
				options.numInfections = atoi(value);
			else if (strcmp(key, "infectingPeriod ") == 0)
				options.infectingPeriod = atoi(value);
			else if (strcmp(key, "randSeed ") == 0)
				options.ranSeed = atoi(value);

			else if (strcmp(key, "neighRadius ") == 0)
				options.neighRadius = atof(value);
			else if (strcmp(key, "nonPermissiveProbability ") == 0)
				options.nonPermProb = atof(value);
			else if (strcmp(key, "initialVirions ") == 0)
				options.initialVirions = atof(value);
			else if (strcmp(key, "virionDiffusion ") == 0)
				options.virionDiffusion = atof(value);
			else if (strcmp(key, "virionClearance ") == 0)
				options.virionClearance = atof(value);
			else if (strcmp(key, "IFNcellProbability ") == 0)
				options.IFNcellProb = atof(value);
			else if (strcmp(key, "IFNdiffusion ") == 0)
				options.IFNdiffusion = atof(value);
			else if (strcmp(key, "IFNclearance ") == 0)
				options.IFNclearance = atof(value);
		}
	}

	fclose(file);
}

