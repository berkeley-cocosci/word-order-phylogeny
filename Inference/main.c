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
	FILE *logfp, *samplesfp, *ancestralsfp;
	int i, j, k, c;
	int logging = 0;
	int multitree, treeindex, treeclass;
	char *treefile, *leaffile, *outdir;
	char filename[1024];
	char mkcmd[1024];
	char families[][16] = {"indo", "austro", "niger", "afro", "nilo", "sino"};
	char types[][16] = {"geographic", "genetic", "feature", "combination" };
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
	gsl_vector_complex *evals = gsl_vector_complex_alloc(6);
        gsl_matrix_complex *evecs = gsl_matrix_complex_alloc(6,6);
        gsl_matrix_complex *evecs_inv = gsl_matrix_complex_alloc(6,6);
	gsl_matrix *Q = gsl_matrix_alloc(6, 6);
	gsl_matrix *P = gsl_matrix_alloc(6, 6);
	gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);


	/* Seed random generator */
	/* (Borrow the logfile handler to save variables) */
	logfp = fopen("/dev/urandom", "r");
	fread(&seed, sizeof(seed), 1, logfp);
	fclose(logfp);
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
	for(i=0; i<6; i++) {
		if(i>0 && multitree == 0) break;
		ancestral_sum[i] = gsl_vector_alloc(6);
		ancestral_max[i] = gsl_vector_alloc(6);
	}

	// TODO check for conflicting options

	// Open log file
	if(logging) {
		logfp = fopen("logfile", "w");
	} else {
		logfp = fopen("/dev/null", "w");
	}
	// Open other files
	sprintf(mkcmd, "mkdir -p %s", outdir);
	system(mkcmd);
	strcpy(filename, outdir);
	strcat(filename, "/samples");
	samplesfp = fopen(filename, "w");
	strcpy(filename, outdir);
	strcat(filename, "/ancestrals");
	ancestralsfp = fopen(filename, "w");
	
	fprintf(logfp, "Default Q:\n");
	fprint_matrix(logfp, Q);

	/* Build trees */
	if(multitree) {
		for(i=0; i<6; i++) {
			sprintf(treefile, "../TreeBuilder/generated_trees/%s/%s/tree_%d.simple", types[treeclass], families[i], treeindex);
			sprintf(leaffile, "../TreeBuilder/generated_trees/%s.leafdata", families[i]);
			trees[i] = build_tree(treefile, leaffile);
			reset_tree(trees[i]);
		}	
	} else {
		trees[0] = build_tree(treefile, leaffile);
		reset_tree(trees[0]);
	}

	/* Compute initial posteriors */
	prior = log(get_prior(stabs, trans));
	likelihood = get_model_likelihood(logfp, trees, Q, multitree);
	posterior = prior + likelihood;
	printf("Initial log posterior: %f\n", posterior);

	/* Burn in */
	for(i=0; i<burnin; i++) posterior = mcmc_iteration(logfp, r, trees, stabs, stabs_dash, trans, trans_dash, posterior, multitree);
	for(i=0; i<samples; i++) {
		if(gsl_rng_uniform_int(r, 10000) >= 9999) {
			/* Random restart! */
			fprintf(logfp, "Random restart!\n");
			initialise_stabs(stabs, 1.0);
			initialise_trans(r, trans);

			// Pick a new random starting point by unconditionally accepting a few proposals
			for(j=0; j<25; j++) draw_proposal(logfp, r, stabs, stabs, trans, trans);
			build_q(Q, stabs, trans);
			prior = log(get_prior(stabs, trans));
			likelihood = get_model_likelihood(logfp, trees, Q, multitree);
			posterior = prior + likelihood;
			printf("New initial log posterior: %f\n", posterior);
			for(j=0; j<burnin; j++) posterior = mcmc_iteration(logfp, r, trees, stabs, stabs_dash, trans, trans_dash, posterior, multitree);
		}

		for(j=0; j<LAG; j++) {
			posterior = mcmc_iteration(logfp, r, trees, stabs, stabs_dash, trans, trans_dash, posterior, multitree);
			if(posterior > max_posterior) {
				max_posterior = posterior;
				printf("Max posterior: %e\n", posterior);
				vector_copy(stabs, stabs_max);
				matrix_copy(trans, trans_max);
			}
		}
		fprintf(logfp, "POSTERTRACK: Sample number %d has log posterior %e.\n", i+1, posterior);
		fprintf(logfp, "I've got %d of %d samples.\n", i, SAMPLES);
		build_q(Q, stabs, trans);
	        upwards_belprop(logfp, trees, Q, multitree);
		decompose_q(Q, evals, evecs, evecs_inv);
		compute_p(evals, evecs, evecs_inv, 0.1, P);
		record_sample(trees, ancestral_sum, stabs, trans, stabs_sum, trans_sum, multitree);

		/* Record sample details */
		fprintf(samplesfp, "Log posterior: %f\n", posterior);
		fprintf(samplesfp, "Sample: %d\n", i+1);
		fprint_vector(samplesfp, stabs);
		fprint_matrix(samplesfp, trans);
		fprint_matrix(samplesfp, P);
		fprintf(samplesfp, "----------\n");

		/* Record ancestral distribution */
		if(multitree) {
			for(j=0; j<6; j++) {
				fprintf(ancestralsfp, families[j]);
				fprintf(ancestralsfp, ": ");
				for(k=0; k<5; k++) {
					fprintf(ancestralsfp, "%f ", trees[j]->dist[k]);
				}
				fprintf(ancestralsfp, "%f\n", trees[j]->dist[5]);
			}
			fprintf(ancestralsfp, "----------\n");
		} else {
			for(k=0; k<5; k++) {
				fprintf(ancestralsfp, "%f ", trees[0]->dist[k]);
			}
			fprintf(ancestralsfp, "%f\n", trees[0]->dist[5]);
		}
	}

	normalise_samples(ancestral_sum, stabs_sum, trans_sum, samples, multitree);
	build_q(Q, stabs_sum, trans_sum);
	//upwards_belprop(logfp, trees, Q, multitree);

	fclose(logfp);
	fclose(samplesfp);
	fclose(ancestralsfp);

	// Save results
	strcpy(filename, outdir);
	strcat(filename, "/summary");
	save_results(filename, Q, ancestral_sum, ancestral_max, stabs_sum, stabs_max, trans_sum, trans_max, max_posterior, multitree);


}

