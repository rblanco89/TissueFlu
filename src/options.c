#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "options.h"

SimOptions options;  // global struct

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
			if (strcmp(key, "gridWidth") == 0) options.gridWidth = atoi(value);
			else if (strcmp(key, "gridHeight") == 0) options.gridHeight = atoi(value);
			else if (strcmp(key, "timeSteps") == 0) options.timeSteps = atoi(value);
		}
	}

	fclose(file);
}

