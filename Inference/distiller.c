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

#define SOV_TO_SVO 0
#define SVO_TO_SOV 1
#define VSO_TO_SOV 2
#define SVO_MOST_STAB 3
#define SOV_MOST_LIKE 4

#include "tree.h"
#include "matrix.h"
#include "beliefprop.h"
#include "modellike.h"
#include "mcmc.h"
#include "saveresults.h"

void handle_directory(char *directory, gsl_vector *stabs_sum, gsl_vector *stabs_map, gsl_matrix *trans_sum, gsl_matrix *trans_map, gsl_matrix *Q_sum, gsl_matrix *Q_map, gsl_matrix *P_sum, gsl_matrix *P_map, gsl_matrix *ancestral_sum, gsl_matrix *ancestral_map, gsl_vector *stationary_sum, gsl_vector *stationary_map, float *statistics, int multitree, int *samplecount) {
	FILE *fp;
	int i, j, k, c;
	int treeindex, treeclass;
	char filename[1024];
	char junkbuffer[1024];
	unsigned long int seed;
	int sample, maxsample;
	float logpost;
	float row[6];
	float maxlogpost = -1000000.0;
	gsl_vector_complex *evals = gsl_vector_complex_alloc(6);
        gsl_matrix_complex *evecs = gsl_matrix_complex_alloc(6,6);
        gsl_matrix_complex *evecs_inv = gsl_matrix_complex_alloc(6,6);
	gsl_vector *stabs = gsl_vector_alloc(6);
	gsl_matrix *trans = gsl_matrix_alloc(6, 6);
	gsl_matrix *Q = gsl_matrix_alloc(6, 6);
	gsl_matrix *P = gsl_matrix_alloc(6, 6);
	gsl_matrix *ancestral = gsl_matrix_alloc(6, 6);

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
		(*samplecount)++;
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
		//fprint_matrix(stdout, Q);
		//fprint_matrix(stdout, Q_sum);
		gsl_matrix_add(Q_sum, Q);
		//fprint_matrix(stdout, Q_sum);
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
		// Collect probabilities
		if(gsl_matrix_get(trans, 0, 1) > gsl_matrix_get(trans, 0, 2))  statistics[SOV_TO_SVO]++;
		if(gsl_matrix_get(trans, 1, 0) > gsl_matrix_get(trans, 1, 2))  statistics[SVO_TO_SOV]++;
		if(gsl_matrix_get(trans, 2, 0) > gsl_matrix_get(trans, 2, 1))  statistics[VSO_TO_SOV]++;
		if(gsl_vector_min_index(stabs) == 1)  statistics[SVO_MOST_STAB]++;
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
	fclose(fp);

	gsl_vector_complex_free(evals);
        gsl_matrix_complex_free(evecs);
        gsl_matrix_complex_free(evecs_inv);
	gsl_vector_free(stabs);
	gsl_matrix_free(trans);
	gsl_matrix_free(Q);
	gsl_matrix_free(P);

}



void write_summary(char *directory, gsl_vector *stabs_sum, gsl_vector *stabs_map, gsl_matrix *trans_sum, gsl_matrix *trans_map, gsl_matrix *Q_sum, gsl_matrix *Q_map, gsl_matrix *P_sum, gsl_matrix *P_map, gsl_matrix *ancestral_sum, gsl_matrix *ancestral_map, gsl_vector *stationary_sum, gsl_vector *stationary_map, float *statistics, int multitree) {
	char filename[1024];
	FILE *fp;
	int i;
	float maxlogpost;
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
	fprintf(fp, "Hypothesis probabilities:\n", maxlogpost);
	fprintf(fp, "SOV to SVO (over VSO): %f\n", statistics[SOV_TO_SVO]);
	fprintf(fp, "SVO to SOV (over VSO): %f\n", statistics[SVO_TO_SOV]);
	fprintf(fp, "VSO to SOV (over SVO): %f\n", statistics[VSO_TO_SOV]);
	fprintf(fp, "SVO most stable: %f\n", statistics[SVO_MOST_STAB]);
//	fprintf(fp, "SOV most likely ancestor: %f\n", statistics[SOV_MOST_LIKE]);
	fprintf(fp, "----------\n");
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

int main(int argc, char **argv) {

	// Variable setup, allocation, etc.
	FILE *fp;
	int i, j, k, c;
	int method, tree;
	int multitree, treeindex, treeclass;
	char *directory;
	char filename[1024];
	char junkbuffer[1024];
	unsigned long int seed;
	node_t **trees = calloc(sizeof(node_t*), 6);
	double likelihood, prior, posterior, max_posterior = -1000000;
	int sample, maxsample, local_samplecount, global_samplecount;
	float logpost;
	float row[6];
	float maxlogpost = -1000000.0;
	char families[][16] = {"indo", "austro", "niger", "afro", "nilo", "sino"};
	char types[][16] = {"geographic", "genetic", "feature", "combination" };

	gsl_vector *local_stabs_sum = gsl_vector_alloc(6);
	gsl_vector *local_stabs_map = gsl_vector_alloc(6);
	gsl_matrix *local_trans = gsl_matrix_alloc(6, 6);
	gsl_matrix *local_trans_sum = gsl_matrix_alloc(6, 6);
	gsl_matrix *local_trans_map = gsl_matrix_alloc(6, 6);
	gsl_matrix *local_Q_sum = gsl_matrix_alloc(6, 6);
	gsl_matrix *local_Q_map = gsl_matrix_alloc(6, 6);
	gsl_matrix *local_P_sum = gsl_matrix_alloc(6, 6);
	gsl_matrix *local_P_map = gsl_matrix_alloc(6, 6);
	gsl_vector *local_stationary_sum = gsl_vector_alloc(6);
	gsl_vector *local_stationary_map = gsl_vector_alloc(6);
	gsl_matrix *local_ancestral_sum = gsl_matrix_alloc(6, 6);
	gsl_matrix *local_ancestral_map = gsl_matrix_alloc(6, 6);
	float local_statistics[20];

	gsl_vector *global_stabs_sum = gsl_vector_alloc(6);
	gsl_vector *global_stabs_map = gsl_vector_alloc(6);
	gsl_matrix *global_trans = gsl_matrix_alloc(6, 6);
	gsl_matrix *global_trans_sum = gsl_matrix_alloc(6, 6);
	gsl_matrix *global_trans_map = gsl_matrix_alloc(6, 6);
	gsl_matrix *global_Q_sum = gsl_matrix_alloc(6, 6);
	gsl_matrix *global_Q_map = gsl_matrix_alloc(6, 6);
	gsl_matrix *global_P_sum = gsl_matrix_alloc(6, 6);
	gsl_matrix *global_P_map = gsl_matrix_alloc(6, 6);
	gsl_vector *global_stationary_sum = gsl_vector_alloc(6);
	gsl_vector *global_stationary_map = gsl_vector_alloc(6);
	gsl_matrix *global_ancestral_sum = gsl_matrix_alloc(6, 6);
	gsl_matrix *global_ancestral_map = gsl_matrix_alloc(6, 6);
	float global_statistics[10];

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

	/* Read all samples */
	for(method=0; method<4; method++) {

		/* New method, so reset all global variables */
		gsl_vector_set_zero(global_stabs_sum);
		gsl_vector_set_zero(global_stationary_sum);
		gsl_matrix_set_zero(global_trans_sum);
		gsl_matrix_set_zero(global_Q_sum);
		gsl_matrix_set_zero(global_P_sum);
		gsl_matrix_set_zero(global_ancestral_sum);
		for(i=0;i<10;i++) global_statistics[i] = 0.0;
		global_samplecount = 0;
		
		for(tree=1; tree<101; tree++) {

			/* New tree, so reset all local vars */

			gsl_vector_set_zero(local_stabs_sum);
			gsl_vector_set_zero(local_stationary_sum);
			gsl_matrix_set_zero(local_trans_sum);
			gsl_matrix_set_zero(local_Q_sum);
			gsl_matrix_set_zero(local_P_sum);
			gsl_matrix_set_zero(local_ancestral_sum);
			for(i=0;i<10;i++) local_statistics[i] = 0.0;
			local_samplecount = 0;
		
			/* Do it! */
			sprintf(filename, "results/individual-q/%s/tree_%d/", types[method], tree);
			//printf("Now handling %s...\n", filename);
			//fprint_matrix(stdout, Q_sum);
			handle_directory(filename, local_stabs_sum, local_stabs_map, local_trans_sum, local_trans_map, local_Q_sum, local_Q_map, local_P_sum, local_P_map, local_ancestral_sum, local_ancestral_map, local_stationary_sum, local_stationary_map, local_statistics, multitree, &local_samplecount);

			/* Normalise local variables */
			gsl_vector_scale(local_stabs_sum, 1.0 / local_samplecount);
			gsl_vector_scale(local_stationary_sum, 1.0 / local_samplecount);
			gsl_matrix_scale(local_trans_sum, 1.0 / local_samplecount);
			gsl_matrix_scale(local_Q_sum, 1.0 / local_samplecount);
			gsl_matrix_scale(local_P_sum, 1.0 / local_samplecount);
			gsl_matrix_scale(local_ancestral_sum, 1.0 / local_samplecount);
			for(i=0;i<10;i++) local_statistics[i] /= local_samplecount;

			/* Summarise local variables */
			write_summary(filename, local_stabs_sum, local_stabs_map, local_trans_sum, local_trans_map, local_Q_sum, local_Q_map, local_P_sum, local_P_map, local_ancestral_sum, local_ancestral_map, local_stationary_sum, local_stationary_map, local_statistics, multitree);

			/* Add normalised local vars to global vars */
			gsl_vector_add(global_stabs_sum, local_stabs_sum);
			gsl_vector_add(global_stationary_sum, local_stationary_sum);
			gsl_matrix_add(global_trans_sum, local_trans_sum);
			gsl_matrix_add(global_Q_sum, local_Q_sum);
			gsl_matrix_add(global_P_sum, local_P_sum);
			gsl_matrix_add(global_ancestral_sum, local_ancestral_sum);
			for(i=0;i<10;i++) global_statistics[i] += local_statistics[i];
			global_samplecount++;
		}

		/* Normalise global variables */
		gsl_vector_scale(global_stabs_sum, 1.0 / global_samplecount);
		gsl_vector_scale(global_stationary_sum, 1.0 / global_samplecount);
		gsl_matrix_scale(global_trans_sum, 1.0 / global_samplecount);
		gsl_matrix_scale(global_Q_sum, 1.0 / global_samplecount);
		gsl_matrix_scale(global_P_sum, 1.0 / global_samplecount);
		gsl_matrix_scale(global_ancestral_sum, 1.0 / global_samplecount);
		for(i=0;i<10;i++) global_statistics[i] /= global_samplecount;

		/* Summarise global variables */
		sprintf(filename, "results/individual-q/%s/", types[method]);
		write_summary(filename, global_stabs_sum, global_stabs_map, global_trans_sum, global_trans_map, global_Q_sum, global_Q_map, global_P_sum, global_P_map, global_ancestral_sum, global_ancestral_map, global_stationary_sum, global_stationary_map, global_statistics, multitree);

	}
}

