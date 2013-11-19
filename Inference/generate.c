#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_matrix_complex_float.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_vector_complex.h>

#include "beliefprop.h"
#include "gslworkspace.h"
#include "matrix.h"
#include "mcmc.h"
#include "modellike.h"
#include "tree.h"
#include "wordorders.h"

void generate_random_mutation_model(mcmc_t *mcmc) {
	int i;

        initialise_stabs(mcmc, 1.0);
        initialise_trans(mcmc);
	for(i=0; i<25; i++) {
		draw_proposal(mcmc);
		gsl_vector_memcpy(mcmc->stabs, mcmc->stabs_dash);
		gsl_matrix_memcpy(mcmc->trans, mcmc->trans_dash);
	}
}

void read_mutation_model(char *filename, gsl_vector *stabs, gsl_matrix *trans) {
	FILE *fp;
	int i, j;
	double x1, x2, x3, x4, x5, x6;
	fp = fopen(filename, "r");
	printf("File pointer is: %p\n", fp);
	i = fscanf(fp, "%lf	%lf	%lf	%lf	%lf	%lf\n", &x1, &x2, &x3, &x4, &x5, &x6);
	printf("fscanf returned: %d\n", i);
	printf("Read in %f, %f, %f, %f, %f, %f\n", x1, x2, x3, x4, x5, x6);
	gsl_vector_set(stabs, 0, x1);
	gsl_vector_set(stabs, 1, x2);
	gsl_vector_set(stabs, 2, x3);
	gsl_vector_set(stabs, 3, x4);
	gsl_vector_set(stabs, 4, x5);
	gsl_vector_set(stabs, 5, x6);
	for(i=0; i<6; i++) {
		fscanf(fp, "%lf	%lf	%lf	%lf	%lf	%lf\n", &x1, &x2, &x3, &x4, &x5, &x6);
		printf("Read in %f, %f, %f, %f, %f, %f\n", x1, x2, x3, x4, x5, x6);
		gsl_matrix_set(trans, i, 0, x1);
		gsl_matrix_set(trans, i, 1, x2);
		gsl_matrix_set(trans, i, 2, x3);
		gsl_matrix_set(trans, i, 3, x4);
		gsl_matrix_set(trans, i, 4, x5);
		gsl_matrix_set(trans, i, 5, x6);
	}	
	fclose(fp);

	for(i=0; i<6; i++) {
		gsl_matrix_set(trans, i, i, 0.0);
	}
	for(i=0; i<6; i++) {
		x1 = 0;
		for(j=0; j<6; j++) {
			printf("Adding %f\n", gsl_matrix_get(trans, i, j));
			x1 += gsl_matrix_get(trans, i, j);
		}
		if(abs(x1 - 1.0) > 1e-10) {
			fprintf(stderr, "Invalid mutation model specified by file %s!\n", filename);
			fprintf(stderr, "Row %d of transition matrix sums to %f, not 1.0!\n", i, x1);
			fprintf(stderr, "Dying now.\n");
			exit(1);
		}
	}
}

void prepare_leaves_for_likelihood(node_t *node) {
	int i;
	if(node->left_child == NULL && node->right_child == NULL) {
		//printf("Did have dist: [%f, %f, %f, %f, %f, %f]\n", node->dist[0], node->dist[1], node->dist[2], node->dist[3], node->dist[4], node->dist[5]);
		for(i=0; i<6; i++) {
			node->dist[i] = node->l_message[i];
		}
		//printf("Now have dist: [%f, %f, %f, %f, %f, %f]\n", node->dist[0], node->dist[1], node->dist[2], node->dist[3], node->dist[4], node->dist[5]);
	} else {
		if(node->left_child != NULL) prepare_leaves_for_likelihood(node->left_child);
		if(node->right_child != NULL) prepare_leaves_for_likelihood(node->right_child);
	}
}

void save_leaf_data(FILE *fp, node_t *node) {
	int wordorder, i;
	if(node->left_child == NULL && node->right_child == NULL) {
		// I'm a leaf!
		//printf("Here's my l_message:\n");
		//printf("[%f, %f, %f, %f, %f, %f]\n", node->l_message[0], node->l_message[1], node->l_message[2], node->l_message[3], node->l_message[4], node->l_message[5]);
		for(i=0; i<6; i++) {
			if(node->l_message[i] != 0.0) wordorder = i;
		}
		//printf("I'm returning: %d\n", wordorder);
		fprintf(fp, "%s\t%d\n", node->nodename, wordorder+1);
	} else {
		if(node->left_child != NULL) save_leaf_data(fp, node->left_child);
		if(node->right_child != NULL) save_leaf_data(fp, node->right_child);
	}
}
	
int sample_6_dist(gsl_rng *r, double *dist) {
	double cumul = 0;
	double sample;
	int i;
	sample = gsl_rng_uniform(r);
	//printf("Sample is: %f\n", sample);
	for(i=0; i<6; i++) {
		cumul += dist[i];
		//printf("Cumul is: %f\n", cumul);
		if(cumul >= sample) return i;
	}	
	return 66;
}

void evolve_down(node_t *node, mcmc_t *mcmc, gslws_t *ws) {
	int i;
	int my_wordorder = 0;
	int child_wordorder = 0;
	while(!node->l_message[my_wordorder]) {
		my_wordorder++;
	}
	if(node->left_child != NULL) {
		compute_p(ws, node->left_branch);
		//fprint_matrix(stdout, P);
		for(i=0; i<6; i++) {
			node->left_child->dist[i] = gsl_matrix_get(ws->P, my_wordorder, i);
			node->left_child->l_message[i] = 0.0;
		}
		child_wordorder = sample_6_dist(mcmc->r, node->left_child->dist);
		//printf("Sampled word order %d\n", child_wordorder);
		node->left_child->l_message[child_wordorder] = 1.0;
		evolve_down(node->left_child, mcmc, ws);
	}
	if(node->right_child != NULL) {
		compute_p(ws, node->right_branch);
		for(i=0; i<6; i++) {
			node->right_child->dist[i] = gsl_matrix_get(ws->P, my_wordorder, i);
			node->right_child->l_message[i] = 0.0;
		}
		child_wordorder = sample_6_dist(mcmc->r, node->right_child->dist);
		//printf("Sampled word order %d\n", child_wordorder);
		node->right_child->l_message[child_wordorder] = 1.0;
		evolve_down(node->right_child, mcmc, ws);
	}
}

void evolve_on_tree(node_t *root, int rootorder, mcmc_t *mcmc, gslws_t *ws)  {
	int i;
        decompose_q(mcmc->Q, ws);
	for(i=0; i<6; i++) {
		root->l_message[i] = 0.0;
	}
	root->l_message[rootorder] = 1.0;
	evolve_down(root, mcmc, ws);
}

int main(int argc, char **argv) {

	// Variable setup, allocation, etc.
	FILE *fp;
	//int fp;
	int i, j, c;
	int rootorder, random_mutation;
	char *treefile, *leaffile, *mutfile, *outfile, *filename;
	char rootorderstr[64];
	unsigned long int seed;
	node_t *tree;
	double likelihood, prior, posterior, max_likelihood = 0;
	float cutoff, cumul;
	mcmc_t mcmc;
	gslws_t ws;
	filename = calloc(64, sizeof(char));
	gsl_vector *stabs = gsl_vector_alloc(6);
	gsl_matrix *trans = gsl_matrix_alloc(6, 6);
	gsl_matrix *Q = gsl_matrix_alloc(6, 6);
	gsl_matrix *P = gsl_matrix_alloc(6, 6);
	gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);
        gsl_vector_complex *evals = gsl_vector_complex_alloc(6);
        gsl_matrix_complex *evecs = gsl_matrix_complex_alloc(6,6);
        gsl_matrix_complex *evecs_inv = gsl_matrix_complex_alloc(6,6);

	fp = fopen("/dev/urandom", "r");
	fread(&seed, sizeof(seed), 1, fp);
	fclose(fp);
	gsl_rng_set(r, seed);


	// Option parsing
	// defaults
	treefile = "tree";
	leaffile = "leafdata";
	mutfile = "mutationmodel";
	outfile = "generated_";
	rootorder = SOV;
	random_mutation = 0;
	while((c = getopt(argc, argv, "a:r:t:l:m:o:x")) != -1) {
		switch(c) {
			case 'r':
				for(i=0; optarg[i]; i++) {
					optarg[i] = tolower(optarg[i]);
				}
				//strncpy(rootorderstr, optarg, 64), 
				//rootorderstr = tolower(rootorderstr);
				if(strcmp(optarg, "sov") == 0) {
					rootorder = SOV;
				} else if(strcmp(optarg, "svo") == 0) {
					rootorder = SVO;
				} else if(strcmp(optarg, "vso") == 0) {
					rootorder = VSO;
				} else if(strcmp(optarg, "vos") == 0) {
					rootorder = VOS;
				} else if(strcmp(optarg, "ovs") == 0) {
					rootorder = OVS;
				} else if(strcmp(optarg, "osv") == 0) {
					rootorder = OSV;
				} else {
					printf("Unrecognised root word order, using SOV\n");
					rootorder = SOV;
				}
				break;
			case 't':
				treefile = optarg;
				break;
			case 'l':
				leaffile = optarg;
				break;
			case 'm':
				mutfile = optarg;
				break;
			case 'o':
				outfile = optarg;
				break;
			case 'x':
				random_mutation = 1;
				break;
		}
	}

	initialise_mcmc(&mcmc);
	alloc_gslws(&ws);

	if(random_mutation) {
		printf("Generating random mutation model...\n");
		generate_random_mutation_model(&mcmc);
	} else {
		printf("Reading mutation model...\n");
		read_mutation_model(mutfile, mcmc.stabs, mcmc.trans);
	}
	printf("Building Q...\n");
	build_q(&mcmc);
	strcpy(filename, outfile);
	strcat(filename, "_matrix");

	// Find stationary and sample root
        decompose_q(mcmc.Q, &ws);
	compute_p(&ws, 1000000);
	cutoff = gsl_rng_uniform(mcmc.r);
	rootorder = 0;
	cumul = gsl_matrix_get(ws.P, 0, 0);
	while(cumul < cutoff) {
		printf("Cutoff: %f\n", cutoff);
		printf("Cumul: %f\n", cumul);
		printf("Root: %d\n", rootorder);
		rootorder++;
		cumul += gsl_matrix_get(ws.P,0, rootorder);
	}

	printf("Saving Q matrix...\n");
	fp = fopen(filename, "w");
	fprintf(fp, "%d\n", rootorder);
	fprintf(fp, "----------\n");
	fprint_matrix(fp, mcmc.Q);
	fclose(fp);

	printf("Building tree...\n");
	tree = build_tree(treefile, leaffile);

	printf("Evolving on tree...\n");
	evolve_on_tree(tree, rootorder, &mcmc, &ws); 

	strcpy(filename, outfile);
	strcat(filename, "_leaves");
	printf("Saving leaf data...\n");
	fp = fopen(filename, "w");
	save_leaf_data(fp, tree);
	fclose(fp);

	prepare_leaves_for_likelihood(tree);	
	fp = fopen("genchatter", "w");
//	likelihood = get_model_likelihood(fp, tree, Q);
	fclose(fp);
	printf("Generated data with:\n");
	prior = get_log_prior(mcmc.stabs, mcmc.trans);
	printf("Prior = %e\n", prior);
	printf("Likelihood = %e\n", likelihood);
	printf("Posterior = %e\n", prior*likelihood);

}

