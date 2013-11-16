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
#include "matrix.h"
#include "tree.h"

void compute_likelihood(node_t *node, double *likelihood, gsl_vector_complex *evals, gsl_matrix_complex *evecs, gsl_matrix_complex *evecs_inv, gsl_matrix *P);
double get_model_loglh(node_t **trees, gsl_matrix *Q, int multitree);

double get_log_prior(gsl_vector *stabs, gsl_matrix *trans) {
	double log_prior = 0;
	int i;
	for(i=0; i<6; i++) {
		log_prior += log(gsl_ran_exponential_pdf(gsl_vector_get(stabs, i), 3.0));
	}
	return log_prior;
}

double get_model_loglh(node_t **trees, gsl_matrix *Q, int multitree) {


	int i, j;
	double likelihood;
	node_t *root;
	gsl_vector_complex *evals = gsl_vector_complex_alloc(6);
        gsl_matrix_complex *evecs = gsl_matrix_complex_alloc(6,6);
        gsl_matrix_complex *evecs_inv = gsl_matrix_complex_alloc(6,6);
        gsl_matrix *P = gsl_matrix_alloc(6,6);

	// Compute eigen-things of Q
	decompose_q(Q, evals, evecs, evecs_inv);

	for(i=0; i<6; i++) {
		if(i>0 && multitree==0) break;
		root = trees[i];
		// Clean up from last time
		reset_tree(root);
		// Set root prior
		for(j=0; j<6; j++) {
			root->dist[j] = 1.0/6.0;
		}

		compute_likelihood(root, &likelihood, evals, evecs, evecs_inv, P);
	}
	gsl_vector_complex_free(evals);
	gsl_matrix_complex_free(evecs);
	gsl_matrix_complex_free(evecs_inv);
	return likelihood;
}

void compute_likelihood(node_t *node, double *likelihood, gsl_vector_complex *evals, gsl_matrix_complex *evecs, gsl_matrix_complex *evecs_inv, gsl_matrix *P) {
	int i, j;
	double nodelikelihood = 0;
	if(node->left_child == NULL && node->right_child == NULL) {
		for(i=0; i<6; i++) {
			nodelikelihood += node->dist[i]*node->l_message[i];
		}
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
		*likelihood += log(nodelikelihood);
	} else {
		/* We're not a leaf!
		   So, update the node->dist[] parts of our children and then recurse
		*/

		// LEFT SIDE
		compute_p(evals, evecs, evecs_inv, node->left_branch, P);
		for(i=0; i<6; i++) {
			node->left_child->dist[i] = 0;
		}
		for(i=0; i<6; i++) {
			// Assume I am i
			for(j=0; j<6; j++) {
				// Adjust my child's probability of being j
				node->left_child->dist[j] += node->dist[i] * gsl_matrix_get(P, i, j);
			}
		}
		compute_likelihood(node->left_child, likelihood, evals, evecs, evecs_inv, P);
		// RIGHT SIDE
		compute_p(evals, evecs, evecs_inv, node->right_branch, P);
		for(i=0; i<6; i++) {
			node->right_child->dist[i] = 0;
		}
		for(i=0; i<6; i++) {
			// Assume I am i
			for(j=0; j<6; j++) {
				// Adjust my child's probability of being j
				node->right_child->dist[j] += node->dist[i] * gsl_matrix_get(P, i, j);
			}
		}
		compute_likelihood(node->right_child, likelihood, evals, evecs, evecs_inv, P);
	}
}

