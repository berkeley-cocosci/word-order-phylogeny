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

#define MIN_PROTOAGE 70000
#define MAX_PROTOAGE 150000
#define STEPS 100.0

#include "gslworkspace.h"
#include "matrix.h"
#include "tree.h"
#include "ubertree.h"

void initialise_ubertree(ubertree_t *ut) {
	ut->sample_count = 0;
	memset(ut->protodist_sum, 0, 100*6*sizeof(double));
}

void update_ubertree(ubertree_t *ut, node_t **trees, gsl_matrix *Q, gslws_t *ws) {
	uint8_t i, j;
	double protoage;
	for(i=0; i<STEPS; i++) {
		protoage = MIN_PROTOAGE + i*(MAX_PROTOAGE-MIN_PROTOAGE)/STEPS;
		compute_protodist(protoage, ut, trees, Q, ws);
		for(j=0;j<6;j++) {
			ut->protodist_sum[i][j] += ut->protodist[j];
		}
	}
	ut->sample_count++;
}
	
float load_tree_age(char *filename) {
	FILE *fp;
	float age;

	fp = fopen(filename, "r");
	fscanf(fp, "%f\n", &age);
	fclose(fp);
	return age;
}

void compute_protodist(double protoage, ubertree_t *ut, node_t **trees, gsl_matrix *Q, gslws_t *ws) {
	double branch_lengths[6];
	int i, tree, proto_order, sub_order; 
	double message, norm;
	for(i=0; i<6; i++) {
		branch_lengths[i] = protoage - trees[i]->age;
		ut->protodist[i] = 1.0;
	}
	for(tree=0; tree<6; tree++) {
		compute_p(ws, branch_lengths[tree]/10000.0);
		for(proto_order=0; proto_order<6; proto_order++) {
			message = 0;
			for(sub_order=0; sub_order<6; sub_order++) {
				message += trees[tree]->dist[sub_order] * gsl_matrix_get(ws->P, proto_order, sub_order);
			}
			ut->protodist[proto_order] *= message;
		}
	}
	norm = 0;
	for(proto_order=0; proto_order<6; proto_order++) {
		norm += ut->protodist[proto_order];
	}
	for(proto_order=0; proto_order<6; proto_order++) {
		ut->protodist[proto_order] /= norm;
	}
}

void save_ubertree(ubertree_t *ut, char *directory) {
	FILE *fp;
	char filename[1024];
	double protoage;
	uint8_t i, j;
	strcpy(filename, directory);
	strcat(filename, "/common_ancestor");
	fp = fopen(filename, "w");
	// For each age of language 
	for(i=0; i<STEPS; i++) {
		protoage = MIN_PROTOAGE + i*(MAX_PROTOAGE-MIN_PROTOAGE)/STEPS;
		for(j=0; j<6; j++) ut->protodist_sum[i][j] /= ut->sample_count;
		fprintf(fp, "%f, %f, %f, %f, %f, %f, %f\n", protoage, ut->protodist_sum[i][0], ut->protodist_sum[i][1], ut->protodist_sum[i][2], ut->protodist_sum[i][3], ut->protodist_sum[i][4], ut->protodist_sum[i][5]);
	}

	fclose(fp);
}
