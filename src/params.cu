#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "headers.h"

Params params;  // global struct

void parse_parameters(const char *filename)
{
	// Set default values
	params.timeSteps = 1000;
	params.numInfections = 1;
	params.initialVirions = 1.0f;
	params.numReplicates = 1;
	params.neighRadius = 2.0f;

	params.infectingPeriod = 1800; // 30 hours
	params.nonPermProb = 0.0f;
	params.IFNcellProb = 0.0f;

	params.pFmax = 1.0f; // max IFN production rate (IFN min^-1 dsRNA^-1)
	params.k_syn = 1.0f; // dsRNA synthesis rate (dsRNA min^-1 virions^-1)
	params.k_deg = 0.15f / 60.0f; // dsRNA degradation rate (min^-1)

	params.flagRefrac = 0; // default: refractory mechanism OFF
	params.flagSupp = 0; // default: infection suppression OFF
	params.flagBP = 0; // default: BP mechanism OFF
	params.flagPF = 0; // default: PF mechanism OFF
	params.flagPorousDiff = 0; // default: porous diffusion OFF
	params.flagCTL = 0; // default: CTL killing mechanism OFF

	params.K_r = 10.0f; // IFN half-max for refractory mechanism (IFN)
	params.K_s = 5.0f; // IFN half-max for suppression mechanism (IFN)
	params.K_v = 3.0f; // virion half-max for infection mechanism (log10(virions))
	params.K_bp = 5.0f; // IFN half-max for BP mechanism (IFN)
	params.K_pf = 5.0f; // IFN half-max for PF mechanism (IFN)
	params.alpha_pf = 1.0f; // PF mechanism enhancement factor (unitless)
	params.nHill = 2.0f; // Hill coefficient (unitless)

	params.rho_T = 0.01f; // CTL autocatalytic growth rate (min^-1)
	params.delta_T = 0.002f; // CTL decay rate (min^-1)
	params.K_ifn = 5.0f; // IFN half-max for CTL growth gating (IFN)
	params.T0 = 1e-4f; // initial systemic CTL pool
	params.K_T = 5.0f; // CTL half-max for killing probability
	params.T_max = 1000.0f; // CTL carrying capacity (logistic cap on T_sys growth)

	params.virionDiffusion = 0.01f; // virion diffusion coefficient (cell diam^2 min^-1)
	params.virionClearance = 0.0115f; // virion clearance rate (min^-1)
	params.IFNdiffusion = 0.04f; // IFN diffusion coefficient (cell diam^2 min^-1)
	params.IFNclearance = 0.005760f; // IFN clearance rate (min^-1)

	params.printSnap = 0;
	params.snapInterval = 60;
	params.measureInterval = 60;
	params.printReplicates = 0;

	params.ranSeed = 42; // default random seed
	
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
				params.timeSteps = atoi(value);
			else if (strcmp(key, "numInfections ") == 0)
				params.numInfections = atoi(value);
			else if (strcmp(key, "infectingPeriod ") == 0)
				params.infectingPeriod = atoi(value);
			else if (strcmp(key, "ranSeed ") == 0)
				params.ranSeed = atoi(value);
			else if (strcmp(key, "printSnapshots ") == 0)
				params.printSnap = atoi(value);
			else if (strcmp(key, "snapshotInterval ") == 0)
				params.snapInterval = atoi(value);
			else if (strcmp(key, "measureInterval ") == 0)
				params.measureInterval = atoi(value);
			else if (strcmp(key, "numReplicates ") == 0)
				params.numReplicates = atoi(value);
			else if (strcmp(key, "flagRefrac ") == 0)
				params.flagRefrac = atoi(value);
			else if (strcmp(key, "flagSupp ") == 0)
				params.flagSupp = atoi(value);
			else if (strcmp(key, "flagBP ") == 0)
				params.flagBP = atoi(value);
			else if (strcmp(key, "flagPF ") == 0)
				params.flagPF = atoi(value);
			else if (strcmp(key, "flagPorousDiff ") == 0)
				params.flagPorousDiff = atoi(value);
			else if (strcmp(key, "flagCTL ") == 0)
				params.flagCTL = atoi(value);
			else if (strcmp(key, "printReplicates ") == 0)
				params.printReplicates = atoi(value);

			else if (strcmp(key, "neighRadius ") == 0)
				params.neighRadius = atof(value);
			else if (strcmp(key, "nonPermissiveProbability ") == 0)
				params.nonPermProb = atof(value);
			else if (strcmp(key, "initialVirions ") == 0)
				params.initialVirions = atof(value);
			else if (strcmp(key, "virionDiffusion ") == 0)
				params.virionDiffusion = atof(value);
			else if (strcmp(key, "virionClearance ") == 0)
				params.virionClearance = atof(value);
			else if (strcmp(key, "IFNcellProbability ") == 0)
				params.IFNcellProb = atof(value);
			else if (strcmp(key, "IFNdiffusion ") == 0)
				params.IFNdiffusion = atof(value);
			else if (strcmp(key, "IFNclearance ") == 0)
				params.IFNclearance = atof(value);
			else if (strcmp(key, "K_r ") == 0)
				params.K_r = atof(value);
			else if (strcmp(key, "K_s ") == 0)
				params.K_s = atof(value);
			else if (strcmp(key, "K_v ") == 0)
				params.K_v = atof(value);
			else if (strcmp(key, "K_bp ") == 0)
				params.K_bp = atof(value);
			else if (strcmp(key, "K_pf ") == 0)
				params.K_pf = atof(value);
			else if (strcmp(key, "alpha_pf ") == 0)
				params.alpha_pf = atof(value);
			else if (strcmp(key, "nHill ") == 0)
				params.nHill = atof(value);
			else if (strcmp(key, "k_syn ") == 0)
				params.k_syn = atof(value);
			else if (strcmp(key, "k_deg ") == 0)
				params.k_deg = atof(value);
			else if (strcmp(key, "pFmax ") == 0)
				params.pFmax = atof(value);
			else if (strcmp(key, "rho_T ") == 0)
				params.rho_T = atof(value);
			else if (strcmp(key, "delta_T ") == 0)
				params.delta_T = atof(value);
			else if (strcmp(key, "K_ifn ") == 0)
				params.K_ifn = atof(value);
			else if (strcmp(key, "T0 ") == 0)
				params.T0 = atof(value);
			else if (strcmp(key, "K_T ") == 0)
				params.K_T = atof(value);
			else if (strcmp(key, "T_max ") == 0)
				params.T_max = atof(value);
			else
				fprintf(stderr, "Warning: Unknown parameter '%s' in config file.\n", key);
		}
	}

	fclose(file);
}

