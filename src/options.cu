#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "options.h"

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
			else if (strcmp(key, "incubationPeriod ") == 0)
				options.incubationPeriod = atoi(value);
			else if (strcmp(key, "expressingPeriod ") == 0)
				options.expressingPeriod = atoi(value);

			else if (strcmp(key, "neighRadius ") == 0)
				options.neighRadius = atof(value);
			else if (strcmp(key, "initialVirions ") == 0)
				options.initialVirions = atof(value);
			else if (strcmp(key, "virionProduction ") == 0)
				options.virionProduction = atof(value);
			else if (strcmp(key, "virionDiffusion ") == 0)
				options.virionDiffusion = atof(value);
			else if (strcmp(key, "virionClearance ") == 0)
				options.virionClearance = atof(value);
		}
	}

	fclose(file);
}

