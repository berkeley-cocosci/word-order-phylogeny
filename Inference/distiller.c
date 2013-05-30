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

int main(int argc, char **argv) {

	// Variable setup, allocation, etc.
	FILE *fp;
	int i, j, k, c;
	int multitree, treeindex, treeclass;
	char *directory;
	char filename[1024];
	char junkbuffer[1024];
	unsigned long int seed;
	node_t **trees = calloc(sizeof(node_t*), 6);
	double likelihood, prior, posterior, max_posterior = -1000000;
	int sample, maxsample;
	float logpost;
	float row[6];
	float maxlogpost = -1000000.0;

	gsl_vector *stabs = gsl_vector_alloc(6);
	gsl_vector *stabs_sum = gsl_vector_alloc(6);
	gsl_vector *stabs_map = gsl_vector_alloc(6);
	gsl_matrix *trans = gsl_matrix_alloc(6, 6);
	gsl_matrix *trans_sum = gsl_matrix_alloc(6, 6);
	gsl_matrix *trans_map = gsl_matrix_alloc(6, 6);
	gsl_vector_complex *evals = gsl_vector_complex_alloc(6);
        gsl_matrix_complex *evecs = gsl_matrix_complex_alloc(6,6);
        gsl_matrix_complex *evecs_inv = gsl_matrix_complex_alloc(6,6);
	gsl_matrix *Q = gsl_matrix_alloc(6, 6);
	gsl_matrix *Q_sum = gsl_matrix_alloc(6, 6);
	gsl_matrix *Q_map = gsl_matrix_alloc(6, 6);
	gsl_matrix *P = gsl_matrix_alloc(6, 6);
	gsl_matrix *P_sum = gsl_matrix_alloc(6, 6);
	gsl_matrix *P_map = gsl_matrix_alloc(6, 6);
	gsl_vector *stationary_sum = gsl_vector_alloc(6);
	gsl_vector *stationary_map = gsl_vector_alloc(6);
	gsl_matrix *ancestral = gsl_matrix_alloc(6, 6);
	gsl_matrix *ancestral_sum = gsl_matrix_alloc(6, 6);
	gsl_matrix *ancestral_map = gsl_matrix_alloc(6, 6);

	/* Handle options */
	directory = ".";
	while((c = getopt(argc, argv, "d:m")) != -1) {
		switch(c) {
			case 'd':
				directory = optarg;
				break;
			case 'm':
				multitree = 1;
				break;
		}
	}

	/* Read sample file */
	strcpy(filename, directory);
	strcat(filename, "/samples");
	fp = fopen(filename, "r");
	if(fp == NULL) {
		printf("Couldn't open %s!\n", filename);
		return;
	}
	while(1) {
		// Read sample number
		i = fscanf(fp, "Sample: %d\n", &sample);	
		if(i == EOF) {
			break;
		}
		// Read posterior prob
		i = fscanf(fp, "Log posterior: %f\n", &logpost);
		if(i == 0) {
			printf("Didn't match!\n");
			return;
		}
		if(logpost > maxlogpost) maxlogpost = logpost;
		// Read stabilities vector
		i = fscanf(fp, "%f %f %f %f %f %f\n", row, &row[1], &row[2], &row[3], &row[4], &row[5]);
		for(i=0; i<6; i++) {
			gsl_vector_set(stabs, i, row[i]);
		}
		gsl_vector_add(stabs_sum, stabs);
		if(logpost == maxlogpost) gsl_vector_memcpy(stabs_map, stabs);
		// Read six rows of trans
		for(i=0; i<6; i++) {
			fscanf(fp, "%f %f %f %f %f %f\n", &row[0], &row[1], &row[2], &row[3], &row[4], &row[5]);
			for(j=0; j<6; j++) {
				gsl_matrix_set(trans, i, j, row[j]);
			}
		}
		gsl_matrix_add(trans_sum, trans);
		if(logpost == maxlogpost) gsl_matrix_memcpy(trans_map, trans);
		// Eat divider line
		fgets(junkbuffer, 1024, fp);
		// Build Q
		build_q(Q, stabs, trans);
		gsl_matrix_add(Q_sum, Q);
		if(logpost == maxlogpost) gsl_matrix_memcpy(Q_map, Q);
		// Build example P
		decompose_q(Q, evals, evecs, evecs_inv);
		compute_p(evals, evecs, evecs_inv, 0.1, P);
		gsl_matrix_add(P_sum, P);
		if(logpost == maxlogpost) gsl_matrix_memcpy(P_map, P);
		// Find stationary
		compute_p(evals, evecs, evecs_inv, 100000000.0, P);
		for(i=0; i<6; i++) {
			gsl_vector_set(stationary_sum, i, gsl_vector_get(stationary_sum, i) + gsl_matrix_get(P, 0, i));
			if(logpost == maxlogpost) gsl_vector_set(stationary_map, i, gsl_matrix_get(P, 0, i));
		}
	}
	fclose(fp);

	/* Read ancestrals file */
	strcpy(filename, directory);
	strcat(filename, "/ancestrals");
	fp = fopen(filename, "r");
	maxsample = sample;

	for(sample=0; sample<maxsample; sample++) {
		if(multitree) {
			for(i=0; i<6; i++) {
				j = fscanf(fp, "%s%f %f %f %f %f %f\n", junkbuffer, &row[0], &row[1], &row[2], &row[3], &row[4], &row[5]);
				for(j=0; j<6; j++) gsl_matrix_set(ancestral, i, j, row[j]);
			}
			fgets(junkbuffer, 1024, fp);
		} else {
			fscanf(fp, "%f %f %f %f %f %f\n", &row[0], &row[1], &row[2], &row[3], &row[4], &row[5]);
			for(j=0; j<6; j++) gsl_matrix_set(ancestral, 0, j, row[j]);
		}
		gsl_matrix_add(ancestral_sum, ancestral);
	}

	/* Normalise sample means */
	for(i=0; i<6; i++) {
		gsl_vector_set(stabs_sum, i, gsl_vector_get(stabs_sum, i) / maxsample);
		gsl_vector_set(stationary_sum, i, gsl_vector_get(stationary_sum, i) / maxsample);
		for(j=0; j<6; j++) {
			gsl_matrix_set(trans_sum, i, j, gsl_matrix_get(trans_sum, i, j) / maxsample);
			gsl_matrix_set(Q_sum, i, j, gsl_matrix_get(Q_sum, i, j) / maxsample);
			gsl_matrix_set(P_sum, i, j, gsl_matrix_get(P_sum, i, j) / maxsample);
			gsl_matrix_set(ancestral_sum, i, j, gsl_matrix_get(ancestral_sum, i, j) / maxsample);
		}
	}

	/* Write summary file */
	strcpy(filename, directory);
	strcat(filename, "/summary");
	fp = fopen(filename, "w");

	fprintf(fp, "Posterior mean stabilities:\n");
	fprint_vector(fp, stabs_sum);
	fprintf(fp, "----------\n");
	fprintf(fp, "Posterior mean transitions:\n");
	fprint_matrix(fp, trans_sum);
	fprintf(fp, "----------\n");
	fprintf(fp, "Posterior mean Q:\n");
	fprint_matrix(fp, Q_sum);
	fprintf(fp, "----------\n");
	fprintf(fp, "P matrix over short branch:\n");
	fprint_matrix(fp, P_sum);
	fprintf(fp, "----------\n");
	fprintf(fp, "Posterior mean stationary:\n");
	fprint_vector(fp, stationary_sum);
	fprintf(fp, "----------\n");
	if(multitree) {
		fprintf(fp, "Posterior mean ancestrals:\n");
		fprint_matrix(fp, ancestral_sum);
		fprintf(fp, "----------\n");
	} else {
		fprintf(fp, "Posterior mean ancestral:\n");
		fprintf(fp, "%f", gsl_matrix_get(ancestral_sum, 0, 0));
		for(i=1;i<6;i++) fprintf(fp, " %f", gsl_matrix_get(ancestral_sum, 0, i));
		fprintf(fp, "\n");
		fprintf(fp, "----------\n");
	}
	fprintf(fp, "Maximum posteror: %f\n", maxlogpost);
	fprintf(fp, "----------\n");
	/*
	if(multitree) {
		fprintf(fp, "MAP ancestral word order distributions:\n");
		for(i=0; i<6; i++) fprint_vector(fp, ancestral_max[i]);
		fprintf(fp, "----------\n");
	} else {
		fprintf(fp, "MAP ancestral word order distribution:\n");
		fprint_vector(fp, ancestral_max[0]);
		fprintf(fp, "----------\n");
	}
	*/
	fprintf(fp, "MAP stabilities:\n");
	fprint_vector(fp, stabs_map);
	fprintf(fp, "MAP transitions:\n");
	fprint_matrix(fp, trans_map);
	fclose(fp);

}

