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

struct summary_context {
	gsl_vector *stabs_sum;
	gsl_vector *stabs_map;
	gsl_matrix *trans;
	gsl_matrix *trans_sum;
	gsl_matrix *trans_map;
	gsl_matrix *Q_sum;
	gsl_matrix *Q_map;
	gsl_matrix *P_sum;
	gsl_matrix *P_map;
	gsl_vector *stationary_sum;
	gsl_vector *stationary_map;
	gsl_matrix *ancestral_sum;
	gsl_matrix *ancestral_map;
	gsl_matrix *fuzz_prior_ancestral_sum;
	gsl_matrix *fuzz_prior_ancestral_map;
	gsl_matrix *stationary_prior_ancestral_sum;
	gsl_matrix *stationary_prior_ancestral_map;
	gsl_matrix **sliding_fuzz_ancestral_sum;
	float statistics[20];
	float sample_count;
	float max_log_post;
};
typedef struct summary_context summary_t;

void alloc_summary(summary_t *s) {
	int i;
	s->stabs_sum = gsl_vector_alloc(6);
	s->stabs_map = gsl_vector_alloc(6);
	s->trans = gsl_matrix_alloc(6, 6);
	s->trans_sum = gsl_matrix_alloc(6, 6);
	s->trans_map = gsl_matrix_alloc(6, 6);
	s->Q_sum = gsl_matrix_alloc(6, 6);
	s->Q_map = gsl_matrix_alloc(6, 6);
	s->P_sum = gsl_matrix_alloc(6, 6);
	s->P_map = gsl_matrix_alloc(6, 6);
	s->stationary_sum = gsl_vector_alloc(6);
	s->stationary_map = gsl_vector_alloc(6);
	s->ancestral_sum = gsl_matrix_alloc(6, 6);
	s->ancestral_map = gsl_matrix_alloc(6, 6);
	s->fuzz_prior_ancestral_sum = gsl_matrix_alloc(6, 6);
	s->fuzz_prior_ancestral_map = gsl_matrix_alloc(6, 6);
	s->stationary_prior_ancestral_sum = gsl_matrix_alloc(6, 6);
	s->stationary_prior_ancestral_map = gsl_matrix_alloc(6, 6);
	s->sliding_fuzz_ancestral_sum = calloc(100, sizeof(gsl_matrix*));
	for(i=0;i<100;i++) s->sliding_fuzz_ancestral_sum[i] = gsl_matrix_alloc(6, 6);
}

void reset_summary(summary_t *s) {
	int i;
	gsl_vector_set_zero(s->stabs_sum);
	gsl_vector_set_zero(s->stabs_map);
	gsl_vector_set_zero(s->stationary_sum);
	gsl_vector_set_zero(s->stationary_map);
	gsl_matrix_set_zero(s->trans_sum);
	gsl_matrix_set_zero(s->trans_map);
	gsl_matrix_set_zero(s->Q_sum);
	gsl_matrix_set_zero(s->Q_map);
	gsl_matrix_set_zero(s->P_sum);
	gsl_matrix_set_zero(s->P_map);
	gsl_matrix_set_zero(s->ancestral_sum);
	gsl_matrix_set_zero(s->ancestral_map);
	gsl_matrix_set_zero(s->fuzz_prior_ancestral_sum);
	gsl_matrix_set_zero(s->fuzz_prior_ancestral_map);
	gsl_matrix_set_zero(s->stationary_prior_ancestral_sum);
	gsl_matrix_set_zero(s->stationary_prior_ancestral_map);
	for(i=0;i<100;i++) gsl_matrix_set_zero(s->sliding_fuzz_ancestral_sum[i]);
	for(i=0;i<10;i++) s->statistics[i] = 0.0;
	s->sample_count = 0;
	s->max_log_post = -1000000000;
}

void normalise_summary(summary_t *s) {
	int i;
	gsl_vector_scale(s->stabs_sum, 1.0 / s->sample_count);
	gsl_vector_scale(s->stationary_sum, 1.0 / s->sample_count);
	gsl_matrix_scale(s->trans_sum, 1.0 / s->sample_count);
	gsl_matrix_scale(s->Q_sum, 1.0 / s->sample_count);
	gsl_matrix_scale(s->P_sum, 1.0 / s->sample_count);
	gsl_matrix_scale(s->ancestral_sum, 1.0 / s->sample_count);
	gsl_matrix_scale(s->fuzz_prior_ancestral_sum, 1.0 / s->sample_count);
	gsl_matrix_scale(s->stationary_prior_ancestral_sum, 1.0 / s->sample_count);
	for(i=0;i<100;i++) gsl_matrix_scale(s->sliding_fuzz_ancestral_sum[i], 1.0 / s->sample_count);
	for(i=0;i<10;i++) s->statistics[i] /= s->sample_count;
}

void add_summaries(summary_t *target, summary_t *source) {
	int i;
	gsl_vector_add(target->stabs_sum, source->stabs_sum);
	gsl_vector_add(target->stationary_sum, source->stationary_sum);
	gsl_matrix_add(target->trans_sum, source->trans_sum);
	gsl_matrix_add(target->Q_sum, source->Q_sum);
	gsl_matrix_add(target->P_sum, source->P_sum);
	gsl_matrix_add(target->ancestral_sum, source->ancestral_sum);
	gsl_matrix_add(target->fuzz_prior_ancestral_sum, source->fuzz_prior_ancestral_sum);
	gsl_matrix_add(target->stationary_prior_ancestral_sum, source->stationary_prior_ancestral_sum);
	for(i=0;i<100;i++) gsl_matrix_add(target->sliding_fuzz_ancestral_sum[i], source->sliding_fuzz_ancestral_sum[i]);
	for(i=0;i<10;i++) target->statistics[i] += source->statistics[i];
	target->sample_count++;
	if(source->max_log_post > target->max_log_post) {
		target->max_log_post = source->max_log_post;
		gsl_vector_memcpy(target->stabs_map, source->stabs_map);
		gsl_vector_memcpy(target->stationary_map, source->stationary_map);
		gsl_matrix_memcpy(target->trans_map, source->trans_map);
		gsl_matrix_memcpy(target->Q_map, source->Q_map);
		gsl_matrix_memcpy(target->P_map, source->P_map);
		gsl_matrix_memcpy(target->ancestral_map, source->ancestral_map);
		gsl_matrix_memcpy(target->ancestral_map, source->ancestral_map);
		gsl_matrix_memcpy(target->fuzz_prior_ancestral_map, source->fuzz_prior_ancestral_map);
		gsl_matrix_memcpy(target->stationary_prior_ancestral_map, source->stationary_prior_ancestral_map);
	}
}

void build_fuzz_prior(gsl_vector *fuzz_prior, gsl_vector_complex *evals, gsl_matrix_complex *evecs, gsl_matrix_complex *evecs_inv) {
	int i, j, k;
	double delta, norm;
	gsl_matrix *P = gsl_matrix_alloc(6, 6);
	gsl_vector *dist = gsl_vector_alloc(6);
	gsl_vector_set_all(fuzz_prior, 1.0/6.0);
	for(i=10; i<=250; i+= 10) {
		gsl_vector_set_zero(dist);
		delta = i / 10000.0;
		compute_p(evals, evecs, evecs_inv, delta, P);
		for(j=0; j<6; j++) {
			for(k=0; k<6; k++) {
				gsl_vector_set(dist, j, gsl_vector_get(dist, j) + gsl_matrix_get(P, k, j));
				if(gsl_vector_get(dist, j) < 0) {
					printf("Yikes, problem!\n");
					printf("Dist is: ");
					fprint_vector(stdout, dist);
					printf("P is:\n");
					fprint_matrix(stdout, P);
					exit(42);
				}
			}
		}
		norm = 0;
		for(j=0; j<6; j++) norm += gsl_vector_get(dist, j);
		gsl_vector_scale(dist, 1.0/norm);
		gsl_vector_add(fuzz_prior, dist);
	}
	norm = 0;
	for(j=0; j<6; j++) norm += gsl_vector_get(fuzz_prior, j);
	gsl_vector_scale(fuzz_prior, 1.0/norm);
	gsl_matrix_free(P);
	gsl_vector_free(dist);
	//fprint_vector(stdout, fuzz_prior);
}

void handle_directory(char *directory, summary_t *s, int multitree) {
	FILE *sample_fp, *ancestral_fp;
	int i, j, k, c;
	int treeindex, treeclass;
	char filename[1024];
	char junkbuffer[1024];
	unsigned long int seed;
	int sample, maxsample;
	float logpost, norm, t, years;
	float row[6];
	gsl_vector_complex *evals = gsl_vector_complex_alloc(6);
        gsl_matrix_complex *evecs = gsl_matrix_complex_alloc(6,6);
        gsl_matrix_complex *evecs_inv = gsl_matrix_complex_alloc(6,6);
	gsl_vector *stabs = gsl_vector_alloc(6);
	gsl_matrix *trans = gsl_matrix_alloc(6, 6);
	gsl_matrix *Q = gsl_matrix_alloc(6, 6);
	gsl_matrix *P = gsl_matrix_alloc(6, 6);
	gsl_matrix *ancestral = gsl_matrix_alloc(6, 6);
	gsl_matrix *prior_ancestral = gsl_matrix_alloc(6, 6);
	gsl_matrix *log_prior_ancestral = gsl_matrix_alloc(6, 6);
	gsl_vector *fuzz_prior = gsl_vector_alloc(6);
	gsl_vector *dist = gsl_vector_alloc(6);

	strcpy(filename, directory);
	strcat(filename, "/samples");
	sample_fp = fopen(filename, "r");
	if(sample_fp == NULL) {
		printf("Couldn't open %s!\n", filename);
		return;
	}

	strcpy(filename, directory);
	strcat(filename, "/ancestrals");
	ancestral_fp = fopen(filename, "r");
	if(ancestral_fp == NULL) {
		printf("Couldn't open %s!\n", filename);
		return;
	}

	while(1) {
		// Read sample number
		i = fscanf(sample_fp, "Sample: %d\n", &sample);	
		if(i == EOF) {
			break;
		}
		s->sample_count++;
		// Read posterior prob
		i = fscanf(sample_fp, "Log posterior: %f\n", &logpost);
		if(i == 0) {
			printf("Didn't match!\n");
			return;
		}
		// Read stabilities vector
		i = fscanf(sample_fp, "%f %f %f %f %f %f\n", row, &row[1], &row[2], &row[3], &row[4], &row[5]);
		for(i=0; i<6; i++) {
			gsl_vector_set(stabs, i, row[i]);
		}
		gsl_vector_add(s->stabs_sum, stabs);
		// Read six rows of trans
		for(i=0; i<6; i++) {
			fscanf(sample_fp, "%f %f %f %f %f %f\n", &row[0], &row[1], &row[2], &row[3], &row[4], &row[5]);
			for(j=0; j<6; j++) {
				gsl_matrix_set(trans, i, j, row[j]);
			}
		}
		gsl_matrix_add(s->trans_sum, trans);
		// Eat divider line
		fgets(junkbuffer, 1024, sample_fp);
		// Build Q
		build_q(Q, stabs, trans);
		gsl_matrix_add(s->Q_sum, Q);
		// Build example P
		decompose_q(Q, evals, evecs, evecs_inv);
		compute_p(evals, evecs, evecs_inv, 0.1, P);
		gsl_matrix_add(s->P_sum, P);
		// Find stationary
		compute_p(evals, evecs, evecs_inv, 100000000.0, P);
		for(i=0; i<6; i++) {
			gsl_vector_set(s->stationary_sum, i, gsl_vector_get(s->stationary_sum, i) + gsl_matrix_get(P, 0, i));
		}
		// Collect probabilities
		if(gsl_matrix_get(trans, 0, 1) > gsl_matrix_get(trans, 0, 2)) s->statistics[SOV_TO_SVO]++;
		if(gsl_matrix_get(trans, 1, 0) > gsl_matrix_get(trans, 1, 2)) s->statistics[SVO_TO_SOV]++;
		if(gsl_matrix_get(trans, 2, 0) > gsl_matrix_get(trans, 2, 1)) s->statistics[VSO_TO_SOV]++;
		if(gsl_vector_min_index(stabs) == 1)  s->statistics[SVO_MOST_STAB]++;
		// Uniform prior ancestrals
		if(multitree) {
			for(i=0; i<6; i++) {
				j = fscanf(ancestral_fp, "%s%f %f %f %f %f %f\n", junkbuffer, &row[0], &row[1], &row[2], &row[3], &row[4], &row[5]);
				for(j=0; j<6; j++) gsl_matrix_set(ancestral, i, j, row[j]);
				for(j=0; j<6; j++) gsl_matrix_set(prior_ancestral, i, j, row[j]);
			}
			fgets(junkbuffer, 1024, ancestral_fp);
		} else {
			fscanf(ancestral_fp, "%f %f %f %f %f %f\n", &row[0], &row[1], &row[2], &row[3], &row[4], &row[5]);
			for(j=0; j<6; j++) gsl_matrix_set(ancestral, 0, j, row[j]);
			for(j=0; j<6; j++) gsl_matrix_set(prior_ancestral, 0, j, row[j]);
		}
		gsl_matrix_add(s->ancestral_sum, ancestral);

		// Fuzzy prior ancestrals
		build_fuzz_prior(fuzz_prior, evals, evecs, evecs_inv);
		if(multitree) {
			for(i=0; i<6; i++) {
				for(j=0; j<6; j++) gsl_matrix_set(prior_ancestral, i, j, gsl_matrix_get(ancestral, i, j)*gsl_vector_get(fuzz_prior, j));
				norm = 0;
				for(j=0; j<6; j++) norm += gsl_matrix_get(prior_ancestral, i, j);
				for(j=0; j<6; j++) gsl_matrix_set(prior_ancestral, i, j, gsl_matrix_get(prior_ancestral, i, j) / norm);

			}
		} else {
			for(j=0; j<6; j++) gsl_matrix_set(prior_ancestral, 0, j, gsl_matrix_get(ancestral, 0, j)*gsl_vector_get(fuzz_prior, j));
			norm = 0;
			for(j=0; j<6; j++) norm += gsl_matrix_get(prior_ancestral, 0, j);
			for(j=0; j<6; j++) gsl_matrix_set(prior_ancestral, 0, j, gsl_matrix_get(prior_ancestral, 0, j) / norm);
		}
		gsl_matrix_add(s->fuzz_prior_ancestral_sum, prior_ancestral);

		// Stationary prior ancestrals
		if(multitree) {
			for(i=0; i<6; i++) {
				for(j=0; j<6; j++) gsl_matrix_set(prior_ancestral, i, j, gsl_matrix_get(ancestral, i, j)*gsl_matrix_get(P, 0, j));
				norm = 0;
				for(j=0; j<6; j++) norm += gsl_matrix_get(prior_ancestral, i, j);
				for(j=0; j<6; j++) gsl_matrix_set(prior_ancestral, i, j, gsl_matrix_get(prior_ancestral, i, j) / norm);

			}
		} else {
			for(j=0; j<6; j++) gsl_matrix_set(prior_ancestral, 0, j, gsl_matrix_get(ancestral, 0, j)*gsl_matrix_get(P, 0, j));
			norm = 0;
			for(j=0; j<6; j++) norm += gsl_matrix_get(prior_ancestral, 0, j);
			for(j=0; j<6; j++) gsl_matrix_set(prior_ancestral, 0, j, gsl_matrix_get(prior_ancestral, 0, j) / norm);
		}
		gsl_matrix_add(s->stationary_prior_ancestral_sum, prior_ancestral);

		// Sliding prior stuff
		//
		for(i=0;i<100;i++) {
			// Compute prior
			years = i*40000/99.0;
			t = years / 10000.0;
			compute_p(evals, evecs, evecs_inv, t, P);
			gsl_vector_set_zero(dist);
			for(j=0; j<6; j++) {
				for(k=0; k<6; k++) {
					gsl_vector_set(dist, j, gsl_vector_get(dist, j) + gsl_matrix_get(P, k, j));
				}
			}
			norm = 0.0;
			for(j=0;j<6;j++) norm += gsl_vector_get(dist, j);
			gsl_vector_scale(dist, 1.0/norm);

			// Multiply likelihood by prior
			if(multitree) {
				for(j=0; j<6; j++) {
					for(k=0; k<6; k++) gsl_matrix_set(prior_ancestral, j, k, gsl_matrix_get(ancestral, j, k)*gsl_vector_get(dist, k));
					norm = 0;
					for(k=0; k<6; k++) norm += gsl_matrix_get(prior_ancestral, j, k);
					for(k=0; k<6; k++) gsl_matrix_set(prior_ancestral, j, k, gsl_matrix_get(prior_ancestral, j, k) / norm);

				}
			} else {
				for(j=0; j<6; j++) gsl_matrix_set(prior_ancestral, 0, j, gsl_matrix_get(ancestral, 0, j)*gsl_vector_get(dist, j));
				norm = 0;
				for(j=0; j<6; j++) norm += gsl_matrix_get(prior_ancestral, 0, j);
				for(j=0; j<6; j++) gsl_matrix_set(prior_ancestral, 0, j, gsl_matrix_get(prior_ancestral, 0, j) / norm);
			}
			gsl_matrix_add(s->sliding_fuzz_ancestral_sum[i], prior_ancestral);
		}

		// MAP stuff
		if(logpost > s->max_log_post) {
			s->max_log_post = logpost;
			gsl_vector_memcpy(s->stabs_map, stabs);
			for(i=0; i<6; i++) gsl_vector_set(s->stationary_map, i, gsl_matrix_get(P, 0, i));
			gsl_matrix_memcpy(s->trans_map, trans);
			gsl_matrix_memcpy(s->Q_map, Q);
			gsl_matrix_memcpy(s->P_map, P);
		}
	}

	fclose(sample_fp);
	fclose(ancestral_fp);

	gsl_vector_complex_free(evals);
	gsl_matrix_complex_free(evecs);
	gsl_matrix_complex_free(evecs_inv);
	gsl_vector_free(stabs);
	gsl_vector_free(fuzz_prior);
	gsl_matrix_free(trans);
	gsl_matrix_free(Q);
	gsl_matrix_free(P);
	gsl_matrix_free(ancestral);
	gsl_matrix_free(prior_ancestral);

}

void write_summary(char *directory, summary_t *s, int multitree) {
	char filename[1024];
	FILE *fp;
	int i;
	float maxlogpost;
	/* Write summary file */
	strcpy(filename, directory);
	strcat(filename, "/summary");
	fp = fopen(filename, "w");

	fprintf(fp, "Posterior mean stabilities:\n");
	fprint_vector(fp, s->stabs_sum);
	fprintf(fp, "----------\n");
	fprintf(fp, "Posterior mean transitions:\n");
	fprint_matrix(fp, s->trans_sum);
	fprintf(fp, "----------\n");
	fprintf(fp, "Posterior mean Q:\n");
	fprint_matrix(fp, s->Q_sum);
	fprintf(fp, "----------\n");
	fprintf(fp, "P matrix over short branch:\n");
	fprint_matrix(fp, s->P_sum);
	fprintf(fp, "----------\n");
	fprintf(fp, "Posterior mean stationary:\n");
	fprint_vector(fp, s->stationary_sum);
	fprintf(fp, "----------\n");
	if(multitree) {
		fprintf(fp, "Posterior mean uniform ancestrals:\n");
		fprint_matrix(fp, s->ancestral_sum);
		fprintf(fp, "----------\n");
		fprintf(fp, "Posterior mean fuzzy ancestrals:\n");
		fprint_matrix(fp, s->fuzz_prior_ancestral_sum);
		fprintf(fp, "----------\n");
		fprintf(fp, "Posterior mean stationary ancestrals:\n");
		fprint_matrix(fp, s->stationary_prior_ancestral_sum);
		fprintf(fp, "----------\n");
	} else {
		fprintf(fp, "Posterior mean ancestral:\n");
		fprintf(fp, "%f", gsl_matrix_get(s->ancestral_sum, 0, 0));
		for(i=1;i<6;i++) fprintf(fp, " %f", gsl_matrix_get(s->ancestral_sum, 0, i));
		fprintf(fp, "\n");
		fprintf(fp, "----------\n");
		fprintf(fp, "Posterior mean fuzzy ancestrals:\n");
		fprintf(fp, "%f", gsl_matrix_get(s->fuzz_prior_ancestral_sum, 0, 0));
		for(i=1;i<6;i++) fprintf(fp, " %f", gsl_matrix_get(s->fuzz_prior_ancestral_sum, 0, i));
		fprintf(fp, "\n");
		fprintf(fp, "----------\n");
		fprintf(fp, "Posterior mean stationary ancestrals:\n");
		fprintf(fp, "%f", gsl_matrix_get(s->stationary_prior_ancestral_sum, 0, 0));
		for(i=1;i<6;i++) fprintf(fp, " %f", gsl_matrix_get(s->stationary_prior_ancestral_sum, 0, i));
		fprintf(fp, "\n");
		fprintf(fp, "----------\n");
	}
	fprintf(fp, "Hypothesis probabilities:\n", s->max_log_post);
	fprintf(fp, "SOV to SVO (over VSO): %f\n", s->statistics[SOV_TO_SVO]);
	fprintf(fp, "SVO to SOV (over VSO): %f\n", s->statistics[SVO_TO_SOV]);
	fprintf(fp, "VSO to SOV (over SVO): %f\n", s->statistics[VSO_TO_SOV]);
	fprintf(fp, "SVO most stable: %f\n", s->statistics[SVO_MOST_STAB]);
//	fprintf(fp, "SOV most likely ancestor: %f\n", statistics[SOV_MOST_LIKE]);
	fprintf(fp, "----------\n");
	fprintf(fp, "Maximum posteror: %f\n", s->max_log_post);
	fprintf(fp, "----------\n");
	fprintf(fp, "MAP stabilities:\n");
	fprint_vector(fp, s->stabs_map);
	fprintf(fp, "MAP transitions:\n");
	fprint_matrix(fp, s->trans_map);
	fprintf(fp, "MAP stationary:\n");
	fprint_vector(fp, s->stationary_map);
	if(multitree) {
		fprintf(fp, "MAP ancestral word order distributions:\n");
		fprint_matrix(fp, s->ancestral_sum);
		fprintf(fp, "----------\n");
	} else {
		fprintf(fp, "MAP ancestral word order distribution:\n");
		fprintf(fp, "%f", gsl_matrix_get(s->ancestral_sum, 0, 0));
		for(i=1;i<6;i++) fprintf(fp, " %f", gsl_matrix_get(s->ancestral_sum, 0, i));
		fprintf(fp, "\n");
		fprintf(fp, "----------\n");
	}

	fclose(fp);
}

void dump_long_prior(char *directory, summary_t *s, int multitree) {
	char filename[1024];
	FILE *fp;
	int i, j;

	/* Write summary file */
	strcpy(filename, directory);
	strcat(filename, "/sliding_prior");
	fp = fopen(filename, "w");

	if(multitree) {
		for(i=0;i<100;i++) {
			fprint_matrix(fp, s->sliding_fuzz_ancestral_sum[i]);
			fprintf(fp,"----------\n");
		}
	} else {
		for(i=0;i<100;i++) {
			fprintf(fp, "%f", gsl_matrix_get(s->sliding_fuzz_ancestral_sum[i], 0, 0));
			for(j=1;j<6;j++) fprintf(fp, " %f", gsl_matrix_get(s->sliding_fuzz_ancestral_sum[i], 0, j));
			fprintf(fp, "\n");
		}
	}

	fclose(fp);
}

int main(int argc, char **argv) {

	// Variable setup, allocation, etc.
	int family, method, tree;
	char filename[1024];
	char families[][16] = {"afro", "austro", "niger", "indo", "nilo", "sino"};
	char types[][16] = {"geographic", "genetic", "feature", "combination" };
	summary_t global_s;
	summary_t local_s;

	alloc_summary(&global_s);
	alloc_summary(&local_s);
	
	/* Common Q data */
	for(method=0; method<4; method++) {
		reset_summary(&global_s);
		for(tree=1; tree<101; tree++) {
			reset_summary(&local_s);
			sprintf(filename, "results/common-q/%s/trees_%d/", types[method], tree);
			printf("Handling %s\n", filename);
			handle_directory(filename, &local_s, 1);
			normalise_summary(&local_s);
			write_summary(filename, &local_s, 1);
			add_summaries(&global_s, &local_s);
		}
	normalise_summary(&global_s);
	sprintf(filename, "results/common-q/%s/", types[method]);
	write_summary(filename, &global_s, 1);
	dump_long_prior(filename, &global_s, 1);
	}

	/* Individual Q data */
	for(family=0; family<6; family++) {
		for(method=0; method<4; method++) {
			reset_summary(&global_s);
			for(tree=1; tree<101; tree++) {
				reset_summary(&local_s);
				sprintf(filename, "results/individual-q/%s/%s/tree_%d/", types[method], families[family], tree);
				printf("Handling %s\n", filename);
				handle_directory(filename, &local_s, 0);
				normalise_summary(&local_s);
				write_summary(filename, &local_s, 0);
				add_summaries(&global_s, &local_s);
			}
		normalise_summary(&global_s);
		sprintf(filename, "results/individual-q/%s/%s/", types[method], families[family]);
		write_summary(filename, &global_s, 0);
		dump_long_prior(filename, &global_s, 0);
		}
	}

}
