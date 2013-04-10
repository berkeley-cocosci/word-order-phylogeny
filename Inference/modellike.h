#include <errno.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


double get_model_likelihood(FILE *fp, node_t **trees, gsl_matrix *Q, int multitree) {


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

		compute_likelihood(fp, root, &likelihood, evals, evecs, evecs_inv, P);
	}

	gsl_vector_complex_free(evals);
	gsl_matrix_complex_free(evecs);
	gsl_matrix_complex_free(evecs_inv);
	return likelihood;
}

compute_likelihood(FILE *fp, node_t *node, double *likelihood, gsl_vector_complex *evals, gsl_matrix_complex *evecs, gsl_matrix_complex *evecs_inv, gsl_matrix *P) {
	int i, j;
	double nodelikelihood = 0;
	if(node->left_child == NULL && node->right_child == NULL) {
		/* We're a leaf! */
//		fprintf(fp, "DIST: [%f, %f, %f, %f, %f, %f]\n", node->dist[0], node->dist[1], node->dist[2], node->dist[3], node->dist[4], node->dist[5]);
//		fprintf(fp, "DATA: [%f, %f, %f, %f, %f, %f]\n", node->l_message[0], node->l_message[1], node->l_message[2], node->l_message[3], node->l_message[4], node->l_message[5]);
		for(i=0; i<6; i++) {
			nodelikelihood += node->dist[i]*node->l_message[i];
		}
		if(nodelikelihood < 0 && nodelikelihood < 1e-15) {
			// If nodelikelihood is negative a shit-tonne, it's most likely supposed to be zero
			// That's bad because we can't return log(0)
			// In this situation, the current parameters are *terrible* and we don't want to accept them
			// So fudge it and set nodelikelihood to something insanely low
			fprintf(stderr, "Mmmm, fudge!\n");
			fprintf(fp, "Mmmm, fudge!\n");
			//dist_printer(fp, node);
			nodelikelihood = 1e-100;
		}
		//fprintf(fp, "BELPROP: Node likelihood: %e\n", nodelikelihood);
		if(nodelikelihood < 0) {
			fprintf(stderr, "compute_likelihood function has computed a negative likelihood!\n");
			fprintf(stderr, "Dying now.\n");
			exit(4);
		}
		if(isnan(nodelikelihood)) {
			fprintf(stderr, "compute_likelihood function has computed a NAN likelihood!\n");
			fprintf(stderr, "Dying now.\n");
			exit(5);
		}
		if(isinf(nodelikelihood)) {
			fprintf(stderr, "compute_likelihood function has computed an infinite likelihood!\n");
			fprintf(stderr, "Dying now.\n");
			exit(6);
		}
//			fprintf(fp, "Here's what I've got:\n");
//			fprintf(fp, "Dist: [%f, %f, %f, %f, %f, %f]\n", node->dist[0], node->dist[1], node->dist[2], node->dist[3], node->dist[4], node->dist[5]);
//			fprintf(fp, "Lmsg: [%f, %f, %f, %f, %f, %f]\n", node->l_message[0], node->l_message[1], node->l_message[2], node->l_message[3], node->l_message[4], node->l_message[5]);
//			fprintf(fp, "Parent node branch legnths: %f (l) and %f (r)\n", node->parent->left_branch, node->parent->right_branch);
			
		*likelihood += log(nodelikelihood);
	} else {
		/* We're not a leaf!
		   So, update the node->dist[] parts of our children and then recurse
		*/

		// LEFT SIDE
		compute_p(evals, evecs, evecs_inv, node->left_branch, P);
//		fprint_matrix(fp, P);
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
		compute_likelihood(fp, node->left_child, likelihood, evals, evecs, evecs_inv, P);
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
		compute_likelihood(fp, node->right_child, likelihood, evals, evecs, evecs_inv, P);
	}
}

