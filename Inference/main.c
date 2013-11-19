#include <unistd.h>
#include<math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_matrix_complex_float.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_vector_complex.h>

#include "gslworkspace.h"
#include "tree.h"
#include "matrix.h"
#include "mcmc.h"
#include "modellike.h"
#include "beliefprop.h"
#include "sampman.h"

void handle_treeset(FILE *logfp, mcmc_t *mcmc, node_t **trees, gslws_t *ws, sampman_t *sm, int burnin, int samples, int lag, int multitree) {
	int i, j;
	/* Burn in */
	for(i=0; i<burnin; i++) mcmc_iteration(logfp, mcmc, trees, ws, multitree);
	/* Take samples */
	for(i=0; i<samples; i++) {
		if(gsl_rng_uniform_int(mcmc->r, 10000) >= 9999) {
			/* Random restart! */
			random_restart(mcmc);
			compute_probabilities(mcmc, trees, ws, multitree);
			for(i=0; i<burnin; i++) mcmc_iteration(logfp, mcmc, trees, ws, multitree);
		}

		for(j=0; j<lag; j++) {
			mcmc_iteration(logfp, mcmc, trees, ws, multitree);
		}

		build_q(mcmc);
		upwards_belprop(logfp, trees, mcmc->Q, ws, multitree);

		/* Record sample */
		process_sample(sm, mcmc, trees);
	}
}

int main(int argc, char **argv) {

	// Variable setup, allocation, etc.
	FILE *logfp;
	int i, c;
	int logging = 0;
	int multitree, treeindex, treeclass;
	char *treefile, *leaffile, *outdir;
	char families[][16] = {"indo", "austro", "niger", "afro", "nilo", "sino"};
	char types[][16] = {"geographic", "genetic", "feature", "combination" };
	int burnin, lag, samples;
	node_t **trees = calloc(sizeof(node_t*), 6);
	mcmc_t mcmc;
	gslws_t ws;
	sampman_t sm;


	// Option parsing
	// defaults
	multitree = 0;
	burnin = 5000;
	lag = 100;
	samples = 1000;
	treefile = NULL;
	leaffile = NULL;
	outdir = ".";
	while((c = getopt(argc, argv, "b:c:i:l:ms:t:o:L")) != -1) {
		switch(c) {
			case 'b':
				burnin = atoi(optarg);
				break;
			case 'c':
				treeclass = atoi(optarg);
				break;
			case 'i':
				treeindex = atoi(optarg);
				break;
			case 'l':
				leaffile = optarg;
				break;
			case 'm':
				multitree = 1;
				break;
			case 's':
				samples = atoi(optarg);
				break;
			case 't':
				treefile = optarg;
				break;
			case 'o':
				outdir = optarg;
				break;
			case 'L':
				logging = 1;
				break;
		}
	}
	if(treefile == NULL) treefile = calloc(sizeof(char), 128);
	if(leaffile == NULL) leaffile = calloc(sizeof(char), 128);

	// TODO check for conflicting options

	// Open log file
	if(logging) {
		logfp = fopen("logfile", "w");
	} else {
		logfp = fopen("/dev/null", "w");
	}
	// Open other files

	initialise_sampman(&sm, outdir);
	initialise_mcmc(&mcmc);
	alloc_gslws(&ws);

	// Loop over trees...
	for(treeindex=0; treeindex<5; treeindex++) {
		/* Build tree(s) */
		if(multitree) {
			for(i=0; i<6; i++) {
				sprintf(treefile, "../TreeBuilder/generated_trees/%s/%s/tree_%d.simple", types[treeclass], families[i], treeindex+1);
				sprintf(leaffile, "../TreeBuilder/generated_trees/%s.leafdata", families[i]);
				trees[i] = build_tree(treefile, leaffile);
				reset_tree(trees[i]);
			}	
		} else {
			trees[0] = build_tree(treefile, leaffile);
			reset_tree(trees[0]);
		}

		/* Draw samples for this tree (set) */
		compute_probabilities(&mcmc, trees, &ws, multitree);
		handle_treeset(logfp, &mcmc, trees, &ws, &sm, burnin, samples, lag, multitree);
		/* Free up tree memory */
		free(trees[0]);
		if(multitree) for(i=1; i<6; i++) free(trees[i]);
	}

	// Finish up
	shout_out(&sm);
	fclose(logfp);
	return 0;
}
