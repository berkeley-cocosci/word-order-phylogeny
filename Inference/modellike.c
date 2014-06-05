#include <errno.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_matrix_complex_float.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_vector_complex.h>

#include "beliefprop.h"
#include "gslworkspace.h"
#include "matrix.h"
#include "modellike.h"
#include "tree.h"

void compute_likelihood(node_t *node, double *likelihood, gslws_t *ws);

double get_log_prior(gsl_vector *stabs, gsl_matrix *trans) {
	double log_prior = 0;
	int i;
	for(i=0; i<6; i++) {
		log_prior += log(gsl_ran_exponential_pdf(gsl_vector_get(stabs, i), 3.0));
	}
	return log_prior;
}

double get_model_loglh(node_t *tree, gsl_matrix *Q, gslws_t *ws) {
	uint8_t i;
	double likelihood;

	// Compute eigen-things of Q
	decompose_q(Q, ws);
	// Clean up from last time
	reset_tree(tree);
	// Set root prior
	likelihood = 0;
	for(i=0; i<6; i++) likelihood += get_tree_likelihood(NULL, tree, i, ws);
	likelihood *= (1.0/6.0);
//	likelihood = 1.0/6.0*get_tree_likelihood(NULL, tree, 0, ws);
//	for(i=0; i<6; i++) tree->dist[i] = 1.0/6.0;
//	compute_likelihood(tree, &likelihood, ws);
	return log(likelihood);
}

double get_balanced_model_loglh(node_t *tree, mcmc_t *mcmc, gslws_t *ws) {
	uint8_t i, j;
	double likelihood;

	// Compute eigen-things of Q
	decompose_q(mcmc->Q, ws);
	// Clean up from last time
	reset_tree(tree);
	// Set root priors
	compute_p(ws, mcmc->dummy_tree_age - tree->age);
	for(i=0; i<6; i++) {
		for(j=0; j<6; j++) {
			tree->dist[i] += gsl_matrix_get(ws->P, j, i);
		}
	}
	likelihood = 0;
	for(i=0; i<6; i++) likelihood += tree->dist[i];
	for(i=0; i<6; i++) tree->dist[i] /= likelihood;
	printf("Got balanced tree prior: ");
	for(i=0; i<6; i++) printf("%f, ", tree->dist[i]);
	printf("\n");
	likelihood = 0;
	for(i=0; i<6; i++) likelihood += get_tree_likelihood(NULL, tree, i, ws);
	likelihood *= (1.0/6.0);
//	likelihood = 1.0/6.0*get_tree_likelihood(NULL, tree, 0, ws);
//	for(i=0; i<6; i++) tree->dist[i] = 1.0/6.0;
//	compute_likelihood(tree, &likelihood, ws);
	return log(likelihood);
}
void compute_likelihood(node_t *node, double *likelihood, gslws_t *ws) {
	int i, j;
	double nodelikelihood = 0;
	if(node->left_child == NULL && node->right_child == NULL) {
		for(i=0; i<6; i++) {
			nodelikelihood += node->dist[i]*node->l_message[i];
		}
		/*
		if(nodelikelihood < 0 && nodelikelihood < 1e-15) {
			// If nodelikelihood is a very negative number, it's most likely supposed to be zero
			// That's bad because we can't return log(0)
			// In this situation, the current parameters are *terrible* and we don't want to accept them
			// So fudge it and set nodelikelihood to something insanely low
			fprintf(stderr, "Mmmm, fudge!\n");
			nodelikelihood = 1e-100;
		}
		if(nodelikelihood < 0) {
			fprintf(stderr, "compute_likelihood function has computed a negative likelihood!\n");
			fprintf(stderr, "Dying now.\n");
			exit(4);
		} else if(nodelikelihood > 1) {
			fprintf(stderr, "compute_likelihood function has computed a likelihood over 1!\n");
			fprintf(stderr, "Dying now.\n");
			exit(42);
		} else if(isnan(nodelikelihood)) {
			fprintf(stderr, "compute_likelihood function has computed a NAN likelihood!\n");
			fprintf(stderr, "Dying now.\n");
			exit(5);
		} else if(isinf(nodelikelihood)) {
			fprintf(stderr, "compute_likelihood function has computed an infinite likelihood!\n");
			fprintf(stderr, "Dying now.\n");
			exit(6);
		}
		*/
		
		*likelihood += log(nodelikelihood);
	} else {
		/* We're not a leaf!
		   So, update the node->dist[] parts of our children and then recurse
		*/

		// LEFT SIDE
		compute_p(ws, node->left_branch);
		memset(node->left_child->dist, 0, 6*sizeof(double));
		for(i=0; i<6; i++) {
			// Assume I am i
			for(j=0; j<6; j++) {
				// Adjust my child's probability of being j
				node->left_child->dist[j] += node->dist[i] * gsl_matrix_get(ws->P, i, j);
			}
		}
		compute_likelihood(node->left_child, likelihood, ws);
		// RIGHT SIDE
		compute_p(ws, node->right_branch);
		memset(node->right_child->dist, 0, 6*sizeof(double));
		for(i=0; i<6; i++) {
			// Assume I am i
			for(j=0; j<6; j++) {
				// Adjust my child's probability of being j
				node->right_child->dist[j] += node->dist[i] * gsl_matrix_get(ws->P, i, j);
			}
		}
		compute_likelihood(node->right_child, likelihood, ws);
	}
}

double dummy_tree_log_lh(mcmc_t *mcmc, gslws_t *ws) {
	int dummy_lang_counts[6] = {377, 110, 44, 15, 9, 2};
	int dummy_eurasia_counts[6] = {71, 34, 0, 1, 0, 0};
	int dummy_aust_counts[6] = {34, 13, 0, 2, 2, 1};
	int dummy_amer_counts[6] = {61+52, 14+18, 35+6, 8+1, 0+7, 0};
	int dummy_papu_counts[6] = {131, 19, 1, 0, 0, 0};
	int dummy_other_counts[6] = {28, 12, 3, 3, 0, 1};
	double eruasia_age = 7;
	double aust_age = 5;
	double amer_age = 1.55;
	double papu_age = 6;
	double other_age = 7;
	int i, j;
	double log_lh = 0;
	double dist1[6] = {0,0,0,0,0,0};
	double dist2[6] = {0,0,0,0,0,0};
	double norm;

	decompose_q(mcmc->Q_dash, ws);
	/* Eurasia and other */
	compute_p(ws, 6);
	for(i=0; i<6; i++) {
		for(j=0; j<6; j++) {
			dist1[i] += gsl_matrix_get(ws->P, j, i);
		}
	}
	norm = 0;
	for(i=0; i<6; i++) norm += dist1[i];
	for(i=0; i<6; i++) dist1[i] /= norm;
	return 200*log(dist1[0]);
	for(i=0; i<6; i++) log_lh += dummy_eurasia_counts[i]*log(dist1[i]);
	for(i=0; i<6; i++) log_lh += dummy_other_counts[i]*log(dist1[i]);
	/* Papuan */
	compute_p(ws, 6);
	for(i=0; i<6; i++) {
		for(j=0; j<6; j++) {
			dist1[i] += gsl_matrix_get(ws->P, j, i);
		}
	}
	norm = 0;
	for(i=0; i<6; i++) norm += dist1[i];
	for(i=0; i<6; i++) dist1[i] /= norm;
	for(i=0; i<6; i++) log_lh += dummy_papu_counts[i]*log(dist1[i]);
	/* Australian */
	compute_p(ws, 5);
	for(i=0; i<6; i++) {
		for(j=0; j<6; j++) {
			dist1[i] += gsl_matrix_get(ws->P, j, i);
		}
	}
	norm = 0;
	for(i=0; i<6; i++) norm += dist1[i];
	for(i=0; i<6; i++) dist1[i] /= norm;
	for(i=0; i<6; i++) log_lh += dummy_aust_counts[i]*log(dist1[i]);
	/* American */
	compute_p(ws, 1.5);
	for(i=0; i<6; i++) {
		for(j=0; j<6; j++) {
			dist1[i] += gsl_matrix_get(ws->P, j, i);
		}
	}
	norm = 0;
	for(i=0; i<6; i++) norm += dist1[i];
	for(i=0; i<6; i++) dist1[i] /= norm;
	for(i=0; i<6; i++) log_lh += dummy_amer_counts[i]*log(dist1[i]);

	return log_lh;
}
