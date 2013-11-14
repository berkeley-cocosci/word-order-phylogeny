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

#define MIN_PROTOAGE 35000
#define MAX_PROTOAGE 100000
#define STEPS 100.0

#include "tree.h"

void load_mean_ancestrals_and_q(char *filename, float (*ancestrals)[6], gsl_matrix *Q) {
	FILE *fp;
	char junkbuffer[1024];
	float row[6];
	int i, j, k;
	fp = fopen(filename, "r");
	if(fp == NULL) {
		printf("Failed to open %s!\n", filename);
		return;
	}
	for(i=0; i<12; i++) fgets(junkbuffer, 1024, fp);
	for(i=0; i<6; i++) {
		fscanf(fp, "%f %f %f %f %f %f\n", &row[0], &row[1], &row[2], &row[3], &row[4], &row[5]);
		for(j=0; j<6; j++) gsl_matrix_set(Q, i, j, row[j]);
	}
	for(i=0; i<13; i++) fgets(junkbuffer, 1024, fp);
	for(i=0; i<6; i++) {
		fscanf(fp, "%f %f %f %f %f %f\n", &ancestrals[i][0], &ancestrals[i][1], &ancestrals[i][2], &ancestrals[i][3], &ancestrals[i][4], &ancestrals[i][5]);
		//printf("Loaded ancestrals: %f %f %f %f %f %f\n", ancestrals[i][0], ancestrals[i][1], ancestrals[i][2], ancestrals[i][3], ancestrals[i][4], ancestrals[i][5]);
	}
	fclose(fp);
}

float load_tree_age(char *filename) {
	FILE *fp;
	float age;

	fp = fopen(filename, "r");
	fscanf(fp, "%f\n", &age);
	fclose(fp);
	return age;
}

void compute_proto_dist(float *proto_dist, float (*sub_dists)[6], float *branch_lengths, gsl_matrix *Q, gsl_matrix *P, gsl_vector_complex *evals, gsl_matrix_complex *evecs, gsl_matrix_complex *evecs_inv) {
	int tree, proto_order, sub_order; 
	float message, norm;
	for(proto_order=0; proto_order<6; proto_order++) {
		proto_dist[proto_order] = 1.0;
	}
	decompose_q(Q, evals, evecs, evecs_inv);
	for(tree=0; tree<6; tree++) {
		compute_p(evals, evecs, evecs_inv, branch_lengths[tree]/10000.0, P);
		//printf("For branch length %f I computed this P:\n", branch_lengths[tree]);
		//printf("My dist: %f %f %f %f %f %f\n", sub_dists[tree][0], sub_dists[tree][1], sub_dists[tree][2], sub_dists[tree][3], sub_dists[tree][4], sub_dists[tree][5]);
		//fprint_matrix(stdout, P);
		
		for(proto_order=0; proto_order<6; proto_order++) {
			message = 0;
			for(sub_order=0; sub_order<6; sub_order++) {
				//printf("My prob %d is: %f\n", sub_order, sub_dists[tree][sub_order]);
				//printf("P from %d: %f\n", proto_order, sub_order, gsl_matrix_get(P, proto_order, sub_order));
				message += sub_dists[tree][sub_order] * gsl_matrix_get(P, proto_order, sub_order);
			}
			//printf("Computed message: %f\n", message);
			proto_dist[proto_order] *= message;
		}
		//printf("Here's the protodist: %f %f %f %f %f %f\n", proto_dist[0], proto_dist[1], proto_dist[2], proto_dist[3], proto_dist[4], proto_dist[5]);
	}
	//printf("Raw proto_dist: %f %f %f %f %f %f\n", proto_dist[0], proto_dist[1], proto_dist[2], proto_dist[3], proto_dist[4], proto_dist[5]);
	norm = 0;
	for(proto_order=0; proto_order<6; proto_order++) {
		norm += proto_dist[proto_order];
	}
	for(proto_order=0; proto_order<6; proto_order++) {
		proto_dist[proto_order] /= norm;
	}
}

int main(int argc, char **argv) {

	// Variable setup, allocation, etc.
	FILE *fp;
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

	float ancestrals[6][6];
	float protodist[6];
	float protodist_sum[6];
	float tree_ages[6];
	float branch_lengths[6];
	float protoage;

	fp = fopen("results/common-q/combination/common_ancestor", "w");

	/* For each age of language */
	for(j=0; j<STEPS; j++) {
		protoage = MIN_PROTOAGE + j*(MAX_PROTOAGE-MIN_PROTOAGE)/STEPS;
		for(i=0; i<6; i++) protodist_sum[i] = 0.0;
		/* For each tree */
		for(i=1; i<101; i++) {
			/* Load ages */
			for(k=0; k<6; k++) {
				sprintf(filename, "../TreeBuilder/generated_trees/combination/%s/tree_%d.age", families[k], i);
				tree_ages[k] = load_tree_age(filename);
				branch_lengths[k] = protoage - tree_ages[k];
			}

			sprintf(filename, "./results/common-q/combination/trees_%d/summary", i);
			load_mean_ancestrals_and_q(filename, ancestrals, Q);
			compute_proto_dist(protodist, ancestrals, branch_lengths, Q, P, evals, evecs, evecs_inv);
			for(k=0; k<6; k++) protodist_sum[k] += protodist[k];
		}
		for(i=0; i<6; i++) protodist_sum[i] /= 100.0;
		fprintf(fp, "%f, %f, %f, %f, %f, %f, %f\n", protoage, protodist_sum[0], protodist_sum[1], protodist_sum[2], protodist_sum[3], protodist_sum[4], protodist_sum[5]);
	}

	fclose(fp);
}
