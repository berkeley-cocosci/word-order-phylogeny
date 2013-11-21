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
		process_sample(sm, mcmc, trees[0]);
	}
}

void whole_shared_q(int method, int shuffle, int burnin, int samples, int lag, char *outdir, int logging) {
	FILE *logfp;
	int treeindex, i;
	int multitree = 1;  // Family variable in load_trees is not used if multitree=1
	int family = 0;  // Family variable in load_trees is not used if multitree=1
	node_t **trees = calloc(sizeof(node_t*), 6);
	mcmc_t mcmc;
	gslws_t ws;
	sampman_t sm;

	// Open log file
	if(logging) {
		logfp = fopen("logfile", "w");
	} else {
		logfp = fopen("/dev/null", "w");
	}

	initialise_sampman(&sm, outdir);
	initialise_mcmc(&mcmc);
	alloc_gslws(&ws);

	// Loop over trees...
	for(treeindex=0; treeindex<5; treeindex++) {
		/* Build tree(s) */
		load_trees(trees, "../TreeBuilder/generated_trees/whole/", method, family, treeindex, shuffle, multitree);
		/* Draw samples for this tree (set) */
		compute_probabilities(&mcmc, trees, &ws, multitree);
		handle_treeset(logfp, &mcmc, trees, &ws, &sm, burnin, samples, lag, multitree);
		/* Free up tree memory */
		free(trees[0]);
		if(multitree) for(i=1; i<6; i++) free(trees[i]);
	}

	// Finish up
	compute_means(&sm);
	save_common_q("results", &sm);
	fclose(logfp);
}

void whole_indiv_q(int method, int shuffle, int burnin, int samples, int lag, char *outdir, int logging) {
	FILE *logfp;
	int treeindex, family;
	node_t **trees = calloc(sizeof(node_t*), 6);
	mcmc_t mcmc;
	gslws_t ws;
	sampman_t sm;

	// Open log file
	if(logging) {
		logfp = fopen("logfile", "w");
	} else {
		logfp = fopen("/dev/null", "w");
	}

	initialise_sampman(&sm, outdir);
	initialise_mcmc(&mcmc);
	alloc_gslws(&ws);

	// Loop over families...
	for(family=0; family<6; family++) {
		// Loop over trees...
		for(treeindex=0; treeindex<5; treeindex++) {
			/* Build tree(s) */
			load_trees(trees, "../TreeBuilder/generated_trees/whole/", method, family, treeindex, shuffle, 0);
			/* Draw samples for this tree (set) */
			compute_probabilities(&mcmc, trees, &ws, 0);
			handle_treeset(logfp, &mcmc, trees, &ws, &sm, burnin, samples, lag, 0);
			/* Free up tree memory */
			free(trees[0]);
		}
		// Finish up
		finish(&sm);
	}

	fclose(logfp);
}

void split_shared_q(int method, int shuffle, int burnin, int samples, int lag, char *outdir, int logging) {
	FILE *logfp;
	int treeindex, i;
	int multitree = 1;  // Family variable in load_trees is not used if multitree=1
	int family = 0;  // Family variable in load_trees is not used if multitree=1
	node_t **trees1 = calloc(sizeof(node_t*), 6);
	node_t **trees2 = calloc(sizeof(node_t*), 6);
	mcmc_t mcmc1, mcmc2;
	gslws_t ws1, ws2;
	sampman_t sm1, sm2;

	// Open log file
	if(logging) {
		logfp = fopen("logfile", "w");
	} else {
		logfp = fopen("/dev/null", "w");
	}

	initialise_sampman(&sm1, outdir);
	initialise_sampman(&sm2, outdir);
	initialise_mcmc(&mcmc1);
	initialise_mcmc(&mcmc2);
	alloc_gslws(&ws1);
	alloc_gslws(&ws2);

	// Loop over trees...
	for(treeindex=0; treeindex<5; treeindex++) {
		/* Build tree(s) */
		load_trees(trees1, "../TreeBuilder/generated_trees/split/1/", method, family, treeindex, shuffle, multitree);
		load_trees(trees2, "../TreeBuilder/generated_trees/split/2/", method, family, treeindex, shuffle, multitree);
		/* Draw samples for this tree (set) */
		compute_probabilities(&mcmc1, trees1, &ws1, multitree);
		compute_probabilities(&mcmc2, trees2, &ws2, multitree);
		handle_treeset(logfp, &mcmc1, trees1, &ws1, &sm1, burnin, samples, lag, multitree);
		handle_treeset(logfp, &mcmc2, trees2, &ws2, &sm2, burnin, samples, lag, multitree);
		/* Free up tree memory */
		free(trees1[0]);
		free(trees2[0]);
		if(multitree) for(i=1; i<6; i++) free(trees1[i]);
		if(multitree) for(i=1; i<6; i++) free(trees2[i]);
	}

	// Finish up
	finish(&sm1);
	finish(&sm2);
	fclose(logfp);
}

int main(int argc, char **argv) {

	// Variable setup, allocation, etc.
	int c;
	int logging = 0;
	int multitree, treeclass, shuffle, split;
	char *outdir;
	int burnin, lag, samples;

	// Option parsing
	// defaults
	shuffle = 0;
	treeclass = 0;
	multitree = 0;
	split = 0;
	burnin = 5000;
	lag = 100;
	samples = 1000;
	outdir = ".";
	while((c = getopt(argc, argv, "b:c:i:l:ms:St:o:Lx")) != -1) {
		switch(c) {
			case 'b':
				burnin = atoi(optarg);
				break;
			case 'c':
				treeclass = atoi(optarg);
				break;
/*			case 'i':
				treeindex = atoi(optarg);
				break;
				*/
			case 'l':
				lag = atoi(optarg);
				break;
			case 'm':
				multitree = 1;
				break;
			case 's':
				samples = atoi(optarg);
				break;
			case 'S':
				shuffle = 1;
				break;
			case 'x':
				split = 1;
				break;
				/*
			case 't':
				treefile = optarg;
				break;
				*/
			case 'o':
				outdir = optarg;
				break;
			case 'L':
				logging = 1;
				break;
		}
	}

	if(multitree && !split) {
		whole_shared_q(treeclass, shuffle, burnin, samples, lag, outdir, logging);
	} else if(multitree && split) {
		split_shared_q(treeclass, shuffle, burnin, samples, lag, outdir, logging);
	} else if(!multitree && !split) {
		whole_indiv_q(treeclass, shuffle, burnin, samples, lag, outdir, logging);
	} else if(!multitree && split) {
		// Not implemented yet
	}

	return 0;
}
