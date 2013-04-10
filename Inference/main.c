#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define BURNIN 5000
#define SAMPLES 100
#define LAG 10

#include "tree.h"
#include "matrix.h"
#include "beliefprop.h"
#include "modellike.h"
#include "mcmc.h"
#include "saveresults.h"

void count_transitions(FILE *fp, node_t *node, int transcount[][6]) {
	if(node->left_child != NULL) {
		transcount[node->ml_order][node->left_child->ml_order]++;
		count_transitions(fp, node->left_child, transcount);
	}
	if(node->right_child != NULL) {
		transcount[node->ml_order][node->right_child->ml_order]++;
		count_transitions(fp, node->right_child, transcount);
	}
}


int main(int argc, char **argv) {

	// Variable setup, allocation, etc.
	FILE *fp;
	//int fp;
	int i, j, c;
	int multitree, treeindex, treeclass;
	char *treefile, *leaffile, *outfile;
	char families[][16] = {"indo", "austro", "niger", "afro", "nilo", "sino"};
	char types[][16] = {"distance", "family", "feature" };
	int burnin, samples;
	unsigned long int seed;
	node_t **trees = calloc(sizeof(node_t*), 6);
	double likelihood, prior, posterior, max_posterior = -1000000;
	gsl_vector *ancestral_sum[6];
	gsl_vector *ancestral_max[6];

	gsl_vector *stabs = gsl_vector_alloc(6);
	gsl_vector *stabs_dash = gsl_vector_alloc(6);
	gsl_vector *stabs_sum = gsl_vector_alloc(6);
	gsl_vector *stabs_max = gsl_vector_alloc(6);
	gsl_matrix *trans = gsl_matrix_alloc(6, 6);
	gsl_matrix *trans_dash = gsl_matrix_alloc(6, 6);
	gsl_matrix *trans_sum = gsl_matrix_alloc(6, 6);
	gsl_matrix *trans_max = gsl_matrix_alloc(6, 6);
	gsl_matrix *Q = gsl_matrix_alloc(6, 6);
	gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);


	/* Seed random generator */
	fp = fopen("/dev/urandom", "r");
	fread(&seed, sizeof(seed), 1, fp);
	fclose(fp);
	gsl_rng_set(r, seed);

	/* Initialise some things */
	initialise_stabs(stabs, 1.0);
	initialise_trans(r, trans);
	build_q(Q, stabs, trans);

	// Option parsing
	// defaults
	multitree = 0;
	burnin = 5000;
	samples = 1000;
	treefile = NULL;
	leaffile = NULL;
	outfile = "output";
	while((c = getopt(argc, argv, "b:c:i:l:ms:t:o:")) != -1) {
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
				outfile = optarg;
				break;
		}
	}
	if(treefile == NULL) treefile = calloc(sizeof(char), 64);
	if(leaffile == NULL) leaffile = calloc(sizeof(char), 64);
	for(i=0; i<6; i++) {
		if(i>0 && multitree == 0) break;
		ancestral_sum[i] = gsl_vector_alloc(6);
		ancestral_max[i] = gsl_vector_alloc(6);
	}

	// TODO check for conflicting options

	// Open chatter file
	fp = fopen("chatter", "w");

	fprintf(fp, "Default Q:\n");
	fprint_matrix(fp, Q);

	/* Build trees */
	if(multitree) {
		for(i=0; i<6; i++) {
			sprintf(treefile, "../TreeBuilder/trees/%s%s%d.simple", families[i], types[treeclass], treeindex);
			sprintf(leaffile, "../TreeBuilder/trees/%s.leafdata", families[i]);
			trees[i] = build_tree(treefile, leaffile);
			reset_tree(trees[i]);
		}	
	} else {
		trees[0] = build_tree(treefile, leaffile);
		reset_tree(trees[0]);
	}

	/* Compute initial posteriors */
	prior = log(get_prior(stabs, trans));
	likelihood = get_model_likelihood(fp, trees, Q, multitree);
	posterior = prior + likelihood;
	printf("Initial log posterior: %f\n", posterior);

	/* Burn in */
	for(i=0; i<burnin; i++) posterior = mcmc_iteration(fp, r, trees, stabs, stabs_dash, trans, trans_dash, posterior, multitree);
	for(i=0; i<samples; i++) {
		if(gsl_rng_uniform_int(r, 10000) >= 9999) {
			/* Random restart! */
			fprintf(fp, "Random restart!\n");
			initialise_stabs(stabs, 1.0);
			initialise_trans(r, trans);

			// Pick a new random starting point by unconditionally accepting a few proposals
			for(j=0; j<25; j++) draw_proposal(fp, r, stabs, stabs, trans, trans);
			build_q(Q, stabs, trans);
			prior = log(get_prior(stabs, trans));
			likelihood = get_model_likelihood(fp, trees, Q, multitree);
			posterior = prior + likelihood;
			printf("New initial log posterior: %f\n", posterior);
			for(j=0; j<burnin; j++) posterior = mcmc_iteration(fp, r, trees, stabs, stabs_dash, trans, trans_dash, posterior, multitree);
		}

		for(j=0; j<LAG; j++) {
			posterior = mcmc_iteration(fp, r, trees, stabs, stabs_dash, trans, trans_dash, posterior, multitree);
			if(posterior > max_posterior) {
				max_posterior = posterior;
				printf("Max posterior: %e\n", posterior);
				vector_copy(stabs, stabs_max);
				matrix_copy(trans, trans_max);
			}
		}
		fprintf(fp, "POSTERTRACK: Sample number %d has log posterior %e.\n", i+1, posterior);
		fprintf(fp, "I've got %d of %d samples.\n", i, SAMPLES);
	        upwards_belprop(fp, trees, Q, multitree);
		record_sample(trees, ancestral_sum, stabs, trans, stabs_sum, trans_sum, multitree);

		fprintf(fp, "Here's the current sample:\n");
		build_q(Q, stabs, trans);
		fprint_matrix(fp, Q);
		fprintf(fp, "-----\n");
		fprintf(fp, "Here's the samples so far:\n");
		fprint_vector(fp, stabs_sum);
		fprint_matrix(fp, trans_sum);
	}

	normalise_samples(ancestral_sum, stabs_sum, trans_sum, samples, multitree);
	build_q(Q, stabs_sum, trans_sum);
	upwards_belprop(fp, trees, Q, multitree);

	fclose(fp);


	// Save results
	save_results(outfile, Q, trees, ancestral_sum, ancestral_max, stabs_sum, stabs_max, trans_sum, trans_max, max_posterior, multitree);


}

