#include <errno.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include<gsl/gsl_matrix.h>
#include<gsl/gsl_matrix_complex_float.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_vector_complex.h>

#include "tree.h"
#include "matrix.h"

void unpassed_detector(node_t *node) {
	if(!node->has_passed) {
		printf("I ain't passed!\n");
	}
	if(node->left_child != NULL) unpassed_detector(node->left_child);
	if(node->right_child != NULL) unpassed_detector(node->right_child);
}

void reset_tree(node_t *node) {
	node->ready_to_feed = 0;
	node->has_passed = 0;
	node->has_cached_matrices = 0;
	node->got_left = 0;
	node->got_right = 0;
	node->has_cached_matrices = 0;
	if(node->left_child == NULL && node->right_child == NULL) {
		/* Leaves are ready to feed, and need their dist reset */
		node->ready_to_feed = 1;
		memcpy(node->dist, node->l_message, 6*sizeof(double));
	} else {
		memset(node->dist, 0, 6*sizeof(double));
		memset(node->cache, 0, 6*sizeof(double));
		memset(node->has_cached, 0, 6*sizeof(int));
		memset(node->l_message, 0, 6*sizeof(double));
		memset(node->r_message, 0, 6*sizeof(double));
		memset(node->p_message, 0, 6*sizeof(double));
	}
	if(node->left_child != NULL) reset_tree(node->left_child);
	if(node->right_child != NULL) reset_tree(node->right_child);
}

void mark_ready_nodes(FILE *fp, node_t *node) {
	int i;
	float norm;
	if(!node->has_passed && node->got_left && node->got_right) {
		node->ready_to_feed = 1;
	}
	if(node->left_child != NULL) mark_ready_nodes(fp, node->left_child);
	if(node->right_child != NULL) mark_ready_nodes(fp, node->right_child);
	return;
		//printf("Marking as ready to feed a node with these messages:\n");
		//printf("LEFT		RIGHT\n");
		//for(i=0; i<6; i++) {
	//		printf("%f	%f\n", node->l_message[i], node->r_message[i]);
	//	}	

		if(node->left_child != NULL && node->right_child != NULL) {
			/* We don't want to do this for leaves, for whom dist is already set */
			norm = 0;
			for(i=0; i<6; i++) {
				node->dist[i] = node->l_message[i]*node->r_message[i];
				norm += node->dist[i];
			}	
			for(i=0; i<6; i++) {
				node->dist[i] /= norm;
			}

/*			for(i=0; i<6; i++) {
				fprintf(fp, "Left %d: %f\n", i, node->l_message[i]);
				fprintf(fp, "Right %d: %f\n", i, node->r_message[i]);
				fprintf(fp, "Dist %d: %f\n", i, node->dist[i]);
			}
*/

		}
//	}
}

void get_ready_nodes(node_t *root, node_t ***firstready, int *nodecount) {
	//printf("Pointer: %p\n", firstready);
	//printf("Nodecount: %d\n", *nodecount);
	if(root->ready_to_feed) {
		(*nodecount)++;
		*firstready = realloc(*firstready, (*nodecount) * sizeof(node_t*));
		(*firstready)[*nodecount-1] = root;
	}
	if(root->left_child != NULL) get_ready_nodes(root->left_child, firstready, nodecount);
	if(root->right_child != NULL) get_ready_nodes(root->right_child, firstready, nodecount);
}

void pass_up(FILE *fp, node_t *node, gslws_t *ws) {
	int i, j;
	double norm;
	if(node->parent == NULL) {
		node->ready_to_feed = 0;
		node->has_passed = 1;
		norm = 0;
		for(i=0; i<6; i++) {
			node->dist[i] = node->l_message[i]*node->r_message[i];
			norm += node->dist[i];
		}
		for(i=0; i<6; i++) {
			node->dist[i] /= norm;
		}
		return;
	}

	if(!(node->left_child == NULL && node->right_child == NULL)) {
		norm = 0;
		for(i=0; i<6; i++) {
			node->dist[i] = node->l_message[i]*node->r_message[i];
			norm += node->dist[i];
		}
		for(i=0; i<6; i++) {
			node->dist[i] /= norm;
		}
	}



	fprintf(fp, "Passing up from a node with this dist:\n");
	for(i=0; i<6; i++) {
		fprintf(fp, "Dist %d: %f\n", i, node->dist[i]);
	}

	if(node->parent->left_child == node) {
		compute_p(ws, node->parent->left_branch);
		for(i=0; i<6; i++) {
			node->parent->l_message[i] = 0;
			for(j=0; j<6; j++) {
				node->parent->l_message[i] += node->dist[j]*gsl_matrix_get(ws->P, i, j);
			}
		}
		norm = 0;
		for(i=0; i<6; i++) {
			norm += node->parent->l_message[i];
		}
		for(i=0; i<6; i++) {
			node->parent->l_message[i] /= norm;
		}	

		fprintf(fp, "GANGNAM: My dist is this: [");
		for(i = 0; i<6; i++) {
			fprintf(fp, "%f  ", node->dist[i]);
		}
		fprintf(fp, "\n");
		fprintf(fp, "GANGNAM: I just set an l message to this: [");
		for(i = 0; i<6; i++) {
			fprintf(fp, "%f  ", node->parent->l_message[i]);
		}
		fprintf(fp, "\n");

		node->parent->got_left = 1;
	} else {
		compute_p(ws, node->parent->right_branch);
		for(i=0; i<6; i++) {
			node->parent->r_message[i] = 0;
			for(j=0; j<6; j++) {
				node->parent->r_message[i] += node->dist[j]*gsl_matrix_get(ws->P, i, j);
			}
		}
		norm = 0;
		for(i=0; i<6; i++) {
			norm += node->parent->r_message[i];
		}
		for(i=0; i<6; i++) {
			node->parent->r_message[i] /= norm;
		}	

		node->parent->got_right = 1;
	}
	node->ready_to_feed = 0;
	node->has_passed = 1;
}

void pass_down(FILE *fp, node_t *node, gslws_t *ws) {
	int i, j;
	double norm;
	/* Update my dist by incorporating my parent's message */
	//fprintf(fp, "BELPROP: Here's my dist before incorporating parental message: [%f, %f, %f, %f, %f, %f]\n", node->dist[0], node->dist[1], node->dist[2], node->dist[3], node->dist[4], node->dist[5]);
	//fprintf(fp, "BELPROP: Here's my parental message: [%f, %f, %f, %f, %f, %f]\n", node->p_message[0], node->p_message[1], node->p_message[2], node->p_message[3], node->p_message[4], node->p_message[5]);
	norm = 0;
	if(node->left_child == NULL && node->right_child == NULL) {
		/* I'm a leaf!  My parent message is it. */
		for(i=0; i<6; i++) {
			node->dist[i] = node->p_message[i];
		}
		fprintf(fp, "BELPROP: Here's a leaf distribution: [%f, %f, %f, %f, %f, %f]\n", node->p_message[0], node->p_message[1], node->p_message[2], node->p_message[3], node->p_message[4], node->p_message[5]);
		return;
	}
	for(i=0; i<6; i++) {
		node->dist[i] *= node->p_message[i];
		norm += node->dist[i];
	}
	for(i=0; i<6; i++) {
		node->dist[i] /= norm;
	}
	//fprintf(fp, "BELPROP: Here's my final dist: [%f, %f, %f, %f, %f, %f]\n", node->dist[0], node->dist[1], node->dist[2], node->dist[3], node->dist[4], node->dist[5]);
	if(node->left_child != NULL) {
		/* Compute the transition matrix to my left child */
		compute_p(ws, node->left_branch);
		/* Compute my normalised pass down message, reusing my r_message space */
		norm = 0;
		for(i=0; i<6; i++) {
			node->r_message[i] = 0;
			for(j=0; j<6; j++) {
				node->r_message[i] += node->dist[j]*gsl_matrix_get(ws->P, j, i);
		}
			norm += node->r_message[i];
		}
		for(i=0; i<6; i++) {
			node->r_message[i] /= norm;
			node->left_child->p_message[i] = node->r_message[i];
		}
		//fprintf(fp, "BELPROP: Here's what I passed down: [%f, %f, %f, %f, %f, %f]\n", node->r_message[0], node->r_message[1], node->r_message[2], node->r_message[3], node->r_message[4], node->r_message[5]);
		pass_down(fp, node->left_child, ws);
	}
	if(node->right_child != NULL) {
		compute_p(ws, node->right_branch);
		/* Compute my normalised pass down message, reusing my r_message space */
		norm = 0;
		for(i=0; i<6; i++) {
			node->r_message[i] = 0;
			for(j=0; j<6; j++) {
				node->r_message[i] += node->dist[j]*gsl_matrix_get(ws->P, j, i);
			}
			norm += node->r_message[i];
		}
		for(i=0; i<6; i++) {
			node->r_message[i] /= norm;
			node->right_child->p_message[i] = node->r_message[i];
		}
		//fprintf(fp, "BELPROP: Here's what I passed down: [%f, %f, %f, %f, %f, %f]\n", node->r_message[0], node->r_message[1], node->r_message[2], node->r_message[3], node->r_message[4], node->r_message[5]);
		pass_down(fp, node->right_child, ws);
	}
}

void broken_compute_likelihood(FILE *fp, node_t *node, double *likelihood) {
	int i;
	double nodelikelihood = 0;
	if(node->left_child == NULL && node->right_child == NULL) {
		/* We're a leaf! */
		//fprintf(fp, "BELPROP: Here's my dist: [%e, %e, %e, %e, %e, %e]\n", node->dist[0], node->dist[1], node->dist[2], node->dist[3], node->dist[4], node->dist[5]);
		//fprintf(fp, "BELPROP: Here's my lmsg: [%e, %e, %e, %e, %e, %e]\n", node->l_message[0], node->l_message[1], node->l_message[2], node->l_message[3], node->l_message[4], node->l_message[5]);
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
		//fprintf(fp, "Running likelihood: %e\n", *likelihood);
	} else {
		if(node->left_child) broken_compute_likelihood(fp, node->left_child, likelihood);
		if(node->right_child) broken_compute_likelihood(fp, node->right_child, likelihood);
	}
}

void dist_printer(FILE *fp, node_t *node) {
	fprintf(fp, "Here's a dist trace from leaf to root:\n");
	fprintf(fp, "Left branch: %f\n", node->left_branch);
	fprintf(fp, "Left: [%f, %f, %f, %f, %f, %f]\n", node->l_message[0], node->l_message[1], node->l_message[2], node->l_message[3], node->l_message[4], node->l_message[5]);
	fprintf(fp, "Right branch: %f\n", node->right_branch);
	fprintf(fp, "Right: [%f, %f, %f, %f, %f, %f]\n", node->r_message[0], node->r_message[1], node->r_message[2], node->r_message[3], node->r_message[4], node->r_message[5]);
	fprintf(fp, "Dist: [%f, %f, %f, %f, %f, %f]\n", node->dist[0], node->dist[1], node->dist[2], node->dist[3], node->dist[4], node->dist[5]);
	fprintf(fp, "-----\n");
	while(node->parent != NULL) {
		node = node->parent;
		fprintf(fp, "Left branch: %f\n", node->left_branch);
		fprintf(fp, "Left: [%f, %f, %f, %f, %f, %f]\n", node->l_message[0], node->l_message[1], node->l_message[2], node->l_message[3], node->l_message[4], node->l_message[5]);
		fprintf(fp, "Right branch: %f\n", node->right_branch);
		fprintf(fp, "Right: [%f, %f, %f, %f, %f, %f]\n", node->r_message[0], node->r_message[1], node->r_message[2], node->r_message[3], node->r_message[4], node->r_message[5]);
		fprintf(fp, "Dist: [%f, %f, %f, %f, %f, %f]\n", node->dist[0], node->dist[1], node->dist[2], node->dist[3], node->dist[4], node->dist[5]);
	fprintf(fp, "-----\n");
	}
}

long double get_tree_likelihood(FILE *fp, node_t *node, int value, gslws_t *ws) {
	long double likelihood, left, right, a, b;
	int i;

	if(node->has_cached[value]) {
		return node->cache[value];
	}

	if(node->left_child == NULL && node->right_child == NULL) {
		likelihood = node->l_message[value] == 1 ? 1.0 : 0.0;
	} else {
		if(! node->has_cached_matrices) {
			compute_p(ws, node->left_branch);
			gsl_matrix_memcpy(node->left_P, ws->P);
			compute_p(ws, node->right_branch);
			gsl_matrix_memcpy(node->right_P, ws->P);
			node->has_cached_matrices = 1;
		}
		for(i=0; i<6; i++) node->left_child->dist[i] = gsl_matrix_get(node->left_P, value, i);
		for(i=0; i<6; i++) node->right_child->dist[i] = gsl_matrix_get(node->right_P, value, i);
		left = 0;
		right = 0;
		for(i=0; i<6; i++) {
			a = node->left_child->dist[i];
			b = get_tree_likelihood(fp, node->left_child, i, ws);
			left += a*b;
			a = node->right_child->dist[i];
			b = get_tree_likelihood(fp, node->right_child, i, ws);
			right += a*b;
		}
		likelihood = left*right;
	}

	node->has_cached[value] = 1;
	node->cache[value] = likelihood;
	return likelihood;
}

void cache_trace(node_t *node) {
	printf("Cache: [%e, %e, %e, %e, %e, %e]\n", node->cache[0], node->cache[1], node->cache[2], node->cache[3], node->cache[4], node->cache[5]);
	if(node->left_child != NULL) cache_trace(node->left_child);
}

void find_likely_interiors(FILE *fp, node_t *node) {
	int i;
	double maxprob = 0;
	for(i=0; i<6; i++) {
		if(node->dist[i] > maxprob) {
			maxprob = node->dist[i];
			node->ml_order = i;
		}	
	}

	if(node->left_child == NULL && node->right_child == NULL) {
		for(i=0; i<6; i++) {
			if(node->dist[i] == 1.0) {
				node->ml_order = i;
			}
		}
	}
	if(node->left_child != NULL) find_likely_interiors(fp, node->left_child);
	if(node->right_child != NULL) find_likely_interiors(fp, node->right_child);
}

void upwards_belprop(FILE *fp, node_t *tree, gsl_matrix *Q, gslws_t *ws) {
	node_t **listhead;
	int listsize = 0;
	int i;
	listhead = NULL;
	double norm;

	reset_tree(tree);

	decompose_q(Q, ws);

	get_ready_nodes(tree, &listhead, &listsize);
	while(listsize != 0) {
		fprintf(fp, "Found %d ready nodes!\n", listsize);
		for(i=0; i<listsize; i++) {
			pass_up(fp, listhead[i], ws);
		}
		mark_ready_nodes(fp, tree);
		listsize = 0;
		get_ready_nodes(tree, &listhead, &listsize);
	}

	if(!tree->has_passed) {
		fprintf(stderr, "Upwards belief propagation terminated before reaching root!\n");
		exit(22);
	}
	
	/* Set root poster */
	norm = 0;
	for(i=0; i<6; i++) {
		tree->dist[i] = tree->l_message[i] * tree->r_message[i];
		norm += tree->dist[i];
	}
	for(i=0; i<6; i++) {
		tree->dist[i] /= norm;
	}
}
