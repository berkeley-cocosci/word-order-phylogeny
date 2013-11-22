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

void do_single_tree_inference(FILE *logfp, mcmc_t *mcmc, node_t *tree, gslws_t *ws, sampman_t *sm, int burnin, int samples, int lag) {
	int i, j;
	/* Burn in */
	for(i=0; i<burnin; i++) single_tree_mcmc_iteration(logfp, mcmc, tree, ws);
	/* Take samples */
	for(i=0; i<samples; i++) {
		if(gsl_rng_uniform_int(mcmc->r, 10000) >= 9999) {
			/* Random restart! */
			random_restart(mcmc);
			compute_single_tree_probabilities(mcmc, tree, ws);
			for(i=0; i<burnin; i++) single_tree_mcmc_iteration(logfp, mcmc, tree, ws);
		}

		for(j=0; j<lag; j++) {
			single_tree_mcmc_iteration(logfp, mcmc, tree, ws);
		}

		build_q(mcmc);
		upwards_belprop(logfp, tree, mcmc->Q, ws);

		/* Record sample */
		process_sample(sm, mcmc, tree);
	}
}

void do_multi_tree_inference(FILE *logfp, mcmc_t *mcmc, node_t **trees, gslws_t *wses, sampman_t *sms, int burnin, int samples, int lag) {
	int i, j;
	/* Burn in */
	for(i=0; i<burnin; i++) multi_tree_mcmc_iteration(logfp, mcmc, trees, wses);
	/* Take samples */
	for(i=0; i<samples; i++) {
		if(gsl_rng_uniform_int(mcmc->r, 10000) >= 9999) {
			/* Random restart! */
			random_restart(mcmc);
			compute_multi_tree_probabilities(mcmc, trees, wses);
			for(i=0; i<burnin; i++) multi_tree_mcmc_iteration(logfp, mcmc, trees, wses);
		}

		for(j=0; j<lag; j++) {
			multi_tree_mcmc_iteration(logfp, mcmc, trees, wses);
		}

		build_q(mcmc);
		for(j=0; j<6; j++) {
			upwards_belprop(logfp, trees[j], mcmc->Q, &wses[j]);
			process_sample(&sms[j], mcmc, trees[j]);
		}
	}
}

void whole_shared_q(int method, int shuffle, int burnin, int samples, int lag, char *outdir, int logging) {
	FILE *logfp;
	int treeindex, i;
	node_t **trees = calloc(6, sizeof(node_t*));
	mcmc_t mcmc;
	gslws_t *wses = calloc(6, sizeof(gslws_t));
	sampman_t *sms = calloc(6, sizeof(sampman_t));

	// Open log file
	if(logging) {
		logfp = fopen("logfile", "w");
	} else {
		logfp = fopen("/dev/null", "w");
	}

	initialise_mcmc(&mcmc);
	for(i=0; i<6; i++) {
		alloc_gslws(&wses[i]);
		initialise_sampman(&sms[i], outdir);
	}

	// Loop over trees...
	for(treeindex=0; treeindex<5; treeindex++) {
		/* Build tree(s) */
		for(i=0; i<6; i++) load_tree(&trees[i], "../TreeBuilder/generated_trees/whole/", method, i, treeindex, shuffle);
		/* Draw samples for this tree (set) */
		compute_multi_tree_probabilities(&mcmc, trees, wses);
		do_multi_tree_inference(logfp, &mcmc, trees, wses, sms, burnin, samples, lag);
		/* Free up tree memory */
		for(i=0; i<6; i++) free(trees[i]);
	}

	// Finish up
	for(i=0; i<6; i++) compute_means(&sms[i]);
	save_common_q("results", sms);
	fclose(logfp);
}

void whole_indiv_q(int method, int shuffle, int burnin, int samples, int lag, char *outdir, int logging) {
	FILE *logfp;
	uint8_t family;
	int treeindex;
	node_t *tree;
	mcmc_t mcmc;
	gslws_t ws;
	sampman_t sm;

	// Open log file
	if(logging) {
		logfp = fopen("logfile", "w");
	} else {
		logfp = fopen("/dev/null", "w");
	}

	initialise_mcmc(&mcmc);
	alloc_gslws(&ws);

	// Loop over families
	for(family=0; family<6; family++) {
		initialise_sampman(&sm, outdir);
		// Loop over trees
		for(treeindex=0; treeindex<5; treeindex++) {
			/* Build tree */
			load_tree(&tree, "../TreeBuilder/generated_trees/whole/", method, family, treeindex, shuffle);
			/* Draw samples for this tree (set) */
			compute_single_tree_probabilities(&mcmc, tree, &ws);
			do_single_tree_inference(logfp, &mcmc, tree, &ws, &sm, burnin, samples, lag);
			/* Free up tree memory */
			free(tree);

			// Finish up
			compute_means(&sm);
			save_indiv_q("results", &sm);
		}
	}

	fclose(logfp);
}

void split_shared_q(int method, int shuffle, int burnin, int samples, int lag, char *outdir, int logging) {
	FILE *logfp;
	int treeindex, i;
	node_t **trees1 = calloc(6, sizeof(node_t*));
	node_t **trees2 = calloc(6, sizeof(node_t*));
	mcmc_t mcmc1, mcmc2;
	gslws_t *wses1 = calloc(6, sizeof(gslws_t));
	gslws_t *wses2 = calloc(6, sizeof(gslws_t));
	sampman_t *sms1 = calloc(6, sizeof(sampman_t));
	sampman_t *sms2 = calloc(6, sizeof(sampman_t));

	// Open log file
	if(logging) {
		logfp = fopen("logfile", "w");
	} else {
		logfp = fopen("/dev/null", "w");
	}

	initialise_mcmc(&mcmc1);
	initialise_mcmc(&mcmc2);
	for(i=0; i<6; i++) {
		alloc_gslws(&wses1[i]);
		alloc_gslws(&wses2[i]);
		initialise_sampman(&sms1[i], outdir);
		initialise_sampman(&sms2[i], outdir);
	}

	// Loop over trees...
	for(treeindex=0; treeindex<5; treeindex++) {
		/* Build tree(s) */
		for(i=0; i<6; i++) {
			load_tree(&trees1[i], "../TreeBuilder/generated_trees/split/1/", method, i, treeindex, shuffle);
			load_tree(&trees2[i], "../TreeBuilder/generated_trees/split/2/", method, i, treeindex, shuffle);
		}
		/* Draw samples for this tree (set) */
		compute_multi_tree_probabilities(&mcmc1, trees1, wses1);
		compute_multi_tree_probabilities(&mcmc2, trees2, wses2);
		do_multi_tree_inference(logfp, &mcmc1, trees1, wses1, sms1, burnin, samples, lag);
		do_multi_tree_inference(logfp, &mcmc2, trees2, wses2, sms2, burnin, samples, lag);
		/* Free up tree memory */
		for(i=0; i<6; i++) {
			free(trees1[i]);
			free(trees2[i]);
		}
	}

	// Finish up
	for(i=0; i<6; i++) {
		compute_means(&sms1[i]);
		compute_means(&sms2[i]);
	}
	save_common_q("results", sms1);
	save_common_q("results", sms2);
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
