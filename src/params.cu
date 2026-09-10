/* TissueFlu - reading the configuration file
 *
 * The file holds one "name = value" pair per line; '#' starts a comment,
 * either on its own line or after a value.
 * To add a new parameter: declare it in headers.h, give it a default below,
 * and add one line to the corresponding table.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "headers.h"

Params params;  // global struct

__host__ Hill make_hill(float K, float n)
{
	Hill h;
	h.K = K;
	h.n = n;
	h.Kn = powf(K, n);
	return h;
}

static void set_defaults(void)
{
	/* Simulation control */
	params.timeSteps = 1000;
	params.numReplicates = 1;
	params.measureInterval = 60;
	params.snapInterval = 60;
	params.printSnap = 0;
	params.ranSeed = 42;
	params.tissueSeed = -1;   // not set: non-permissive cells redrawn per replicate

	/* Tissue and initial condition */
	params.numInfections = 1;
	params.infectingPeriod = 1800;   // 30 hours
	params.neighRadius = 2.0f;
	params.nonPermProb = 0.0f;
	params.initialVirions = 1.0f;
	params.IFNcellProb = 0.0f;

	/* Mechanism switches: all feedbacks off by default */
	params.flagRefrac = 0;
	params.flagSupp = 0;
	params.flagBP = 0;
	params.flagPF = 0;
	params.flagPorousDiff = 0;
	params.flagCTL = 0;

	/* Extracellular transport */
	params.virionDiffusion = 0.01f;
	params.virionClearance = 0.0115f;
	params.IFNdiffusion = 0.04f;
	params.IFNclearance = 0.005760f;

	/* Intracellular IFN circuit */
	params.pFmax = 1.0f;
	params.k_syn = 1.0f;
	params.k_deg = 0.15f/60.0f;
	params.alpha_pf = 1.0f;

	/* Systemic CTL compartment */
	params.rho_T = 0.01f;
	params.delta_T = 0.002f;
	params.T0 = 1e-4f;
	params.T_max = 1000.0f;

	/* Dose-response curves (half-maxima; the exponent is shared) */
	params.nHill     = 2.0f;
	params.refrac    = make_hill(10.0f, 2.0f);
	params.suppress  = make_hill( 5.0f, 2.0f);
	params.infect    = make_hill( 3.0f, 2.0f);
	params.blockProd = make_hill( 5.0f, 2.0f);
	params.posFeed   = make_hill( 5.0f, 2.0f);
	params.ctlKill   = make_hill( 5.0f, 2.0f);
	params.ctlGrowth = make_hill( 5.0f, 2.0f);
}

/* Strip leading and trailing blanks in place. */
static char *trim(char *s)
{
	while (*s == ' ' || *s == '\t') s++;
	char *end = s + strlen(s);
	while (end > s && (end[-1] == ' ' || end[-1] == '\t' ||
	                   end[-1] == '\n' || end[-1] == '\r')) end--;
	*end = '\0';
	return s;
}

void parse_parameters(const char *filename)
{
	set_defaults();

	struct { const char *name; int *dst; } intParams[] =
	{
		{"timeSteps",        &params.timeSteps},
		{"numReplicates",    &params.numReplicates},
		{"measureInterval",  &params.measureInterval},
		{"snapshotInterval", &params.snapInterval},
		{"printSnapshots",   &params.printSnap},
		{"ranSeed",          &params.ranSeed},
		{"tissueSeed",       &params.tissueSeed},
		{"numInfections",    &params.numInfections},
		{"infectingPeriod",  &params.infectingPeriod},
		{"flagRefrac",       &params.flagRefrac},
		{"flagSupp",         &params.flagSupp},
		{"flagBP",           &params.flagBP},
		{"flagPF",           &params.flagPF},
		{"flagPorousDiff",   &params.flagPorousDiff},
		{"flagCTL",          &params.flagCTL},
	};

	struct { const char *name; float *dst; } floatParams[] =
	{
		{"neighRadius",              &params.neighRadius},
		{"nonPermissiveProbability", &params.nonPermProb},
		{"initialVirions",           &params.initialVirions},
		{"IFNcellProbability",       &params.IFNcellProb},
		{"virionDiffusion",          &params.virionDiffusion},
		{"virionClearance",          &params.virionClearance},
		{"IFNdiffusion",             &params.IFNdiffusion},
		{"IFNclearance",             &params.IFNclearance},
		{"pFmax",                    &params.pFmax},
		{"k_syn",                    &params.k_syn},
		{"k_deg",                    &params.k_deg},
		{"alpha_pf",                 &params.alpha_pf},
		{"rho_T",                    &params.rho_T},
		{"delta_T",                  &params.delta_T},
		{"T0",                       &params.T0},
		{"T_max",                    &params.T_max},
		{"nHill",                    &params.nHill},
		/* half-maxima of the dose-response curves */
		{"K_r",   &params.refrac.K},
		{"K_s",   &params.suppress.K},
		{"K_v",   &params.infect.K},
		{"K_bp",  &params.blockProd.K},
		{"K_pf",  &params.posFeed.K},
		{"K_T",   &params.ctlKill.K},
		{"K_ifn", &params.ctlGrowth.K},
	};

	const int numInts = sizeof(intParams)/sizeof(intParams[0]);
	const int numFloats = sizeof(floatParams)/sizeof(floatParams[0]);

	FILE *file = fopen(filename, "r");
	if (!file)
	{
		perror("Could not open config file");
		exit(1);
	}

	char line[256];
	while (fgets(line, sizeof(line), file))
	{
		char *comment = strchr(line, '#');
		if (comment) *comment = '\0';

		char *eq = strchr(line, '=');
		if (!eq) continue;
		*eq = '\0';

		char *key = trim(line);
		char *value = trim(eq + 1);
		if (!*key) continue;

		int found = 0;
		for (int i = 0; i < numInts && !found; i++)
			if (strcmp(key, intParams[i].name) == 0)
			{ *intParams[i].dst = atoi(value); found = 1; }
		for (int i = 0; i < numFloats && !found; i++)
			if (strcmp(key, floatParams[i].name) == 0)
			{ *floatParams[i].dst = (float)atof(value); found = 1; }

		if (!found)
			fprintf(stderr, "Warning: Unknown parameter '%s' in config file.\n", key);
	}
	fclose(file);

	/* All dose-response curves share the exponent nHill; cache K^n now that
	   both K and n are known. */
	params.refrac    = make_hill(params.refrac.K,    params.nHill);
	params.suppress  = make_hill(params.suppress.K,  params.nHill);
	params.infect    = make_hill(params.infect.K,    params.nHill);
	params.blockProd = make_hill(params.blockProd.K, params.nHill);
	params.posFeed   = make_hill(params.posFeed.K,   params.nHill);
	params.ctlKill   = make_hill(params.ctlKill.K,   params.nHill);
	params.ctlGrowth = make_hill(params.ctlGrowth.K, params.nHill);
}
