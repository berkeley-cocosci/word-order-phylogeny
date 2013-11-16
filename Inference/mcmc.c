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
#include<gsl/gsl_eigen.h>
#include<gsl/gsl_permutation.h>
#include<math.h>

#include "tree.h"
#include "matrix.h"
#include "mcmc.h"
#include "modellike.h"


void initialise_mcmc(mcmc_t *mcmc) {
	unsigned long int seed;
	FILE *fp;
	/* Allocate memory */
	mcmc->r = gsl_rng_alloc(gsl_rng_taus);
	mcmc->stabs = gsl_vector_alloc(6);
	mcmc->stabs_dash = gsl_vector_alloc(6);
	mcmc->stabs_max = gsl_vector_alloc(6);
	mcmc->trans = gsl_matrix_alloc(6, 6);
	mcmc->trans_dash = gsl_matrix_alloc(6, 6);
	mcmc->trans_max = gsl_matrix_alloc(6, 6);
	mcmc->Q = gsl_matrix_alloc(6, 6);
	mcmc->Q_dash = gsl_matrix_alloc(6, 6);

	/* Seed PRNG */
	fp = fopen("/dev/urandom", "r");
	fread(&seed, sizeof(seed), 1, fp);
	fclose(fp);
	gsl_rng_set(mcmc->r, seed);

	/* Initialise some variables */
	initialise_stabs(mcmc, 1.0);
	initialise_trans(mcmc);
	build_q(mcmc);
}

void compute_probabilities(mcmc_t *mcmc, node_t **trees, int multitree) {
	mcmc->log_prior = get_log_prior(mcmc->stabs, mcmc->trans);
	mcmc->log_lh = get_model_loglh(trees, mcmc->Q, multitree);
	mcmc->log_poster = mcmc->log_prior + mcmc->log_lh;
}

void random_restart(mcmc_t *mcmc) {
	uint8_t i;
	initialise_stabs(mcmc, 1.0);
	initialise_trans(mcmc);
	for(i=0; i<25; i++) draw_proposal(mcmc);
}

void initialise_stabs(mcmc_t *mcmc, double x) {
	uint8_t i;
	for(i=0; i<6; i++) {
		gsl_vector_set(mcmc->stabs, i, x);
	}
}

void initialise_trans(mcmc_t *mcmc) {
	int i, j;
	double norm;
	gsl_matrix_set_zero(mcmc->trans);
	for(i=0; i<6; i++) {
		for(j=0; j<6; j++) {
			if(i != j) gsl_matrix_set(mcmc->trans, i, j, 1.0/5);
		}
		norm = 0;
		for(j=0; j<6; j++) {
			norm+= gsl_matrix_get(mcmc->trans, i, j);
		}
		for(j=0; j<6; j++) {
			gsl_matrix_set(mcmc->trans, i, j, gsl_matrix_get(mcmc->trans, i, j)/norm);
		}
	}
}

void build_q(mcmc_t *mcmc) {
	int i, j;
	double diag, prod;
	for(i=0; i<6; i++) {
		diag = 0.0;
		for(j=0; j<6; j++) {
			if(i != j) {
				prod = gsl_vector_get(mcmc->stabs, i)*gsl_matrix_get(mcmc->trans, i, j);
				gsl_matrix_set(mcmc->Q, i, j, prod);
				diag -= prod;
			}
		}
		gsl_matrix_set(mcmc->Q, i, i, diag);
	}
}

void build_q_dash(mcmc_t *mcmc) {
	int i, j;
	double diag, prod;
	for(i=0; i<6; i++) {
		diag = 0.0;
		for(j=0; j<6; j++) {
			if(i != j) {
				prod = gsl_vector_get(mcmc->stabs_dash, i)*gsl_matrix_get(mcmc->trans_dash, i, j);
				gsl_matrix_set(mcmc->Q_dash, i, j, prod);
				diag -= prod;
			}
		}
		gsl_matrix_set(mcmc->Q_dash, i, i, diag);
	}
}

void swap_column_step(mcmc_t *mcmc) {
	double var;
	int i, j, k;
	i = gsl_rng_uniform_int(mcmc->r, 6);
	j = i;
	while(j == i) j = gsl_rng_uniform_int(mcmc->r, 6);
	k = i;
	while(k == i || k == j) k = gsl_rng_uniform_int(mcmc->r, 6);
	var = gsl_matrix_get(mcmc->Q_dash, j, i);
	gsl_matrix_set(mcmc->Q_dash, j, i, gsl_matrix_get(mcmc->Q_dash, k, i));
	gsl_matrix_set(mcmc->Q_dash, k, i, var);
}

void swap_row_step(mcmc_t *mcmc) {
	double var;
	int i, j, k;
	i = gsl_rng_uniform_int(mcmc->r, 6);
	j = i;
	while(j == i) j = gsl_rng_uniform_int(mcmc->r, 6);
	k = i;
	while(k == i || k == j) k = gsl_rng_uniform_int(mcmc->r, 6);
	var = gsl_matrix_get(mcmc->Q, i, j);
	gsl_matrix_set(mcmc->Q_dash, i, j, gsl_matrix_get(mcmc->Q_dash, i, k));
	gsl_matrix_set(mcmc->Q_dash, i, k, var);
}

void add_step(mcmc_t *mcmc, double variance) {
	double var;
	int i, j;
	i = gsl_rng_uniform_int(mcmc->r, 6);
	j = i;
	while(j == i) j = gsl_rng_uniform_int(mcmc->r, 6);
	var = gsl_ran_gaussian(mcmc->r, variance);
	gsl_matrix_set(mcmc->Q_dash, i, j, fabs(gsl_matrix_get(mcmc->Q_dash, i, j) + var));
}

void draw_proposal(mcmc_t *mcmc) {
	double var;
	int i, j;

	/* Build Q out of current stabs and trans */
	/* Copy Q into Q dash */
	build_q(mcmc);
	gsl_matrix_memcpy(mcmc->Q_dash, mcmc->Q);

	/* Change Q dash */
	i = gsl_rng_uniform_int(mcmc->r, 8);
	switch(i) {
		case 0:
			add_step(mcmc, 0.001);
			break;
		case 1:
			add_step(mcmc, 0.01);
			break;
		case 2:
			add_step(mcmc, 0.10);
			break;
		case 3:
			add_step(mcmc, 1.00);
			break;
		case 4:
		case 5:
			swap_row_step(mcmc);
			break;
		case 6:
		case 7:
			swap_column_step(mcmc);
			break;
	}

	/* Break Q_dash up into stabs_dash and trans_dash */
	for(i=0; i<6; i++) {
		var = 0;
		for(j=0; j<6; j++) {
			if(i != j) var += gsl_matrix_get(mcmc->Q_dash, i, j);
		}
		gsl_vector_set(mcmc->stabs_dash, i, var);
		for(j=0; j<6; j++) {
			if(i != j) gsl_matrix_set(mcmc->trans_dash, i, j, gsl_matrix_get(mcmc->Q_dash, i, j)/var);
		}
	}
}

void handle_acceptance_probability(double *a) {
	if(*a<0) {
		fprintf(stderr, "Acceptance probability is negative!\n");
		fprintf(stderr, "Dying now.\n");
		exit(7);
	}
	if(isnan(*a)) {
		fprintf(stderr, "Acceptance probability is NAN!\n");
		fprintf(stderr, "Dying now.\n");
		exit(8);
	}
	if(isinf(*a)) {
		fprintf(stderr, "Acceptance probability is infinite!\n");
		fprintf(stderr, "Dodgily accepting!\n");
		*a = 1.1;
	}
}

void mcmc_iteration(FILE *fp, mcmc_t *mcmc, node_t **trees, int multitree) {
	double new_log_prior, new_log_lh, new_log_poster, a, sample;

	build_q(mcmc);
	fprintf(fp, "Time to take a closer look...\n");
	fprintf(fp, "Current Q:\n");
	fprint_matrix(fp, mcmc->Q);
	fprintf(fp, "Current posterior: %e\n", mcmc->log_poster);

	draw_proposal(mcmc);

	fprintf(fp, "Proposed stabs:\n");
	fprint_vector(fp, mcmc->stabs_dash);
	fprintf(fp, "Proposed trans:\n");
	fprint_matrix(fp, mcmc->trans_dash);
	build_q_dash(mcmc);
	fprintf(fp, "Proposed Q:\n");
	fprint_matrix(fp, mcmc->Q_dash);

	new_log_prior = get_log_prior(mcmc->stabs_dash, mcmc->trans_dash);
	new_log_lh = get_model_loglh(trees, mcmc->Q_dash, multitree);
	new_log_poster = new_log_prior + new_log_lh;

	// Acceptance probability
	a = exp(new_log_poster - mcmc->log_poster);
	handle_acceptance_probability(&a);
	fprintf(fp, "Acceptance probability is: %f\n", a);
	if(a >= 1) {
		// Accept
		vector_copy(mcmc->stabs_dash, mcmc->stabs);
		matrix_copy(mcmc->trans_dash, mcmc->trans);
		mcmc->log_prior = new_log_prior;
		mcmc->log_lh = new_log_lh;
		mcmc->log_poster = new_log_poster;
	} else {
		sample = gsl_rng_uniform(mcmc->r);
		if(sample <= a) {
			vector_copy(mcmc->stabs_dash, mcmc->stabs);
			matrix_copy(mcmc->trans_dash, mcmc->trans);
		}
	}

}
