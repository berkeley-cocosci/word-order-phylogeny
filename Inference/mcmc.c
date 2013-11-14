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
#include "modellike.h"

void initialise_stabs(gsl_vector *stabs, double x) {
	int i;
	for(i=0; i<6; i++) {
		gsl_vector_set(stabs, i, x);
	}
}

void initialise_trans(gsl_rng *r, gsl_matrix *trans) {
	int i, j;
	double norm;
	gsl_matrix_set_zero(trans);
	for(i=0; i<6; i++) {
		for(j=0; j<6; j++) {
			if(i != j) gsl_matrix_set(trans, i, j, 1.0/5);
		}
		norm = 0;
		for(j=0; j<6; j++) {
			norm+= gsl_matrix_get(trans, i, j);
		}
		for(j=0; j<6; j++) {
			gsl_matrix_set(trans, i, j, gsl_matrix_get(trans, i, j)/norm);
		}
	}
}

void build_q(gsl_matrix *Q, gsl_vector *stabs, gsl_matrix *trans) {
	int i, j;
	double diag, prod;
	for(i=0; i<6; i++) {
		diag = 0.0;
		for(j=0; j<6; j++) {
			if(i != j) {
				prod = gsl_vector_get(stabs, i)*gsl_matrix_get(trans, i, j);
				gsl_matrix_set(Q, i, j, prod);
				diag -= prod;
			}
		}
		gsl_matrix_set(Q, i, i, diag);
	}
}
double get_prior(gsl_vector *stabs, gsl_matrix *trans) {
	double prior = 1;
	int i;
	for(i=0; i<6; i++) {
		prior *= gsl_ran_exponential_pdf(gsl_vector_get(stabs, i), 3.0);
	}
	return prior;
}

void swap_column_step(FILE *fp, gsl_rng *r, gsl_matrix *Q) {
	double var;
	int i, j, k;
	i = gsl_rng_uniform_int(r, 6);
	j = i;
	while(j == i) j = gsl_rng_uniform_int(r, 6);
	k = i;
	while(k == i || k == j) k = gsl_rng_uniform_int(r, 6);
	var = gsl_matrix_get(Q, j, i);
	gsl_matrix_set(Q, j, i, gsl_matrix_get(Q, k, i));
	gsl_matrix_set(Q, k, i, var);
}

void swap_row_step(FILE *fp, gsl_rng *r, gsl_matrix *Q) {
	double var;
	int i, j, k;
	i = gsl_rng_uniform_int(r, 6);
	j = i;
	while(j == i) j = gsl_rng_uniform_int(r, 6);
	k = i;
	while(k == i || k == j) k = gsl_rng_uniform_int(r, 6);
	var = gsl_matrix_get(Q, i, j);
	gsl_matrix_set(Q, i, j, gsl_matrix_get(Q, i, k));
	gsl_matrix_set(Q, i, k, var);
}

void add_step(FILE *fp, gsl_rng *r, gsl_matrix *Q, double variance) {
	double var;
	int i, j;
	i = gsl_rng_uniform_int(r, 6);
	j = i;
	while(j == i) j = gsl_rng_uniform_int(r, 6);
	var = gsl_ran_gaussian(r, variance);
	gsl_matrix_set(Q, i, j, fabs(gsl_matrix_get(Q, i, j) + var));
}

void draw_proposal(FILE *fp, gsl_rng *r, gsl_vector *stabs, gsl_vector *stabs_dash, gsl_matrix *trans, gsl_matrix *trans_dash) {
	gsl_matrix *Q = gsl_matrix_alloc(6, 6);
	double var;
	int i, j;

	/* Turn stabs and trans into Q */	
	build_q(Q, stabs, trans);
	vector_copy(stabs, stabs_dash);
	matrix_copy(trans, trans_dash);

	/* Change Q */
	i = gsl_rng_uniform_int(r, 8);
	switch(i) {
		case 0:
			add_step(fp, r, Q, 0.001);
			break;
		case 1:
			add_step(fp, r, Q, 0.01);
			break;
		case 2:
			add_step(fp, r, Q, 0.10);
			break;
		case 3:
			add_step(fp, r, Q, 1.00);
			break;
		case 4:
		case 5:
			swap_row_step(fp, r, Q);
			break;
		case 6:
		case 7:
			swap_column_step(fp, r, Q);
			break;
	}

	/* Break Q back into stabs and trans */
	for(i=0; i<6; i++) {
		var = 0;
		for(j=0; j<6; j++) {
			if(i != j) var += gsl_matrix_get(Q, i, j);
		}
		gsl_vector_set(stabs_dash, i, var);
		for(j=0; j<6; j++) {
			if(i != j) gsl_matrix_set(trans_dash, i, j, gsl_matrix_get(Q, i, j)/var);
		}
	}
	gsl_matrix_free(Q);
}

void broken_draw_proposal(FILE *fp, gsl_rng *r, gsl_vector *stabs, gsl_vector *stabs_dash, gsl_matrix *trans, gsl_matrix *trans_dash) {
	int i, j;
	double new, delta;
	vector_copy(stabs, stabs_dash);
	matrix_copy(trans, trans_dash);
	if(gsl_rng_uniform_int(r, 11) >= 8) {
		/* Change a stability */
		i = gsl_rng_uniform_int(r, 6);
		delta = gsl_ran_gaussian(r, 0.1);
		gsl_vector_set(stabs_dash, i, fabs(gsl_vector_get(stabs_dash, i) + delta));
	} else {
		/* Change a transition probability */
		fprintf(fp, "Changing transition probability!");
		i = gsl_rng_uniform_int(r, 6);
		j = i;
		while(j == i) {
			j = gsl_rng_uniform_int(r, 6);
		}
		delta = gsl_ran_gaussian(r, 0.3);
		new = fabs(gsl_matrix_get(trans_dash, i, j) + delta);
		gsl_matrix_set(trans_dash, i, j, new);
		/* Renormalise */
		delta = 0;
		for(j=0; j<6; j++) {
			delta += gsl_matrix_get(trans_dash, i, j);
		}
		for(j=0; j<6; j++) {
			gsl_matrix_set(trans_dash, i, j, gsl_matrix_get(trans_dash, i, j)/delta);
		}
	}
}


double mcmc_iteration(FILE *fp, gsl_rng* r, node_t **trees, gsl_vector *stabs, gsl_vector *stabs_dash, gsl_matrix *trans, gsl_matrix *trans_dash, double old_posterior, int multitree) {
	double new_likelihood, new_prior, new_posterior, a, sample;
	gsl_matrix *Q = gsl_matrix_alloc(6, 6);

	build_q(Q, stabs, trans);
	fprintf(fp, "Time to take a closer look...\n");
	fprintf(fp, "Current Q:\n");
	fprint_matrix(fp, Q);
	fprintf(fp, "Current posterior: %e\n", old_posterior);

	draw_proposal(fp, r, stabs, stabs_dash, trans, trans_dash);

	fprintf(fp, "Proposed stabs:\n");
	fprint_vector(fp, stabs_dash);
	fprintf(fp, "Proposed trans:\n");
	fprint_matrix(fp, trans_dash);
	build_q(Q, stabs_dash, trans_dash);
	fprintf(fp, "Proposed Q:\n");
	fprint_matrix(fp, Q);

	new_prior = log(get_prior(stabs_dash, trans_dash));
	fprintf(fp, "Proposal prior is: %e\n", new_prior);
	new_likelihood = get_model_loglh(fp, trees, Q, multitree);
	fprintf(fp, "ACHTUNG I received: %e\n", new_likelihood);
	fprintf(fp, "Proposal loglh is: %e\n", new_likelihood);
	new_posterior = new_prior + new_likelihood;
	fprintf(fp, "Proposal posterior: %e\n", new_posterior);

	if(new_likelihood == 0.0) {
		fprintf(fp, "Exactly zero likelihood!\n");
		fprintf(fp, "This could be a problem...\n");
	}
	a = exp(new_posterior - old_posterior);
	if(a<0) {
		fprintf(stderr, "Acceptance probability is negative!\n");
		fprintf(stderr, "Dying now.\n");
		exit(7);
	}
	if(isnan(a)) {
		fprintf(stderr, "Acceptance probability is NAN!\n");
		fprintf(stderr, "Dying now.\n");
		exit(8);
	}
	if(isinf(a)) {
		fprintf(stderr, "Acceptance probability is infinite!\n");
		fprintf(stderr, "Dodgily accepting!\n");
		a = 1.1;
	}

	fprintf(fp, "Acceptance probability is: %f\n", a);
	if(a >= 1) {
		fprintf(fp, "MCMC decision: accept\n");
		fprintf(fp, "Posterior probability of last accepted transition: %e\n", new_posterior);
		vector_copy(stabs_dash, stabs);
		matrix_copy(trans_dash, trans);
		return new_posterior;
	}
	sample = gsl_rng_uniform(r);
	if(sample <= a) {
		fprintf(fp, "MCMC decision: accept\n");
		vector_copy(stabs_dash, stabs);
		matrix_copy(trans_dash, trans);
		return new_posterior;
	} else {
		fprintf(fp, "MCMC decision: reject\n");
		return old_posterior;
	}

}

void record_sample(node_t **root, gsl_vector **ancestral_sum, gsl_vector *stabs, gsl_matrix *trans, gsl_vector *stabs_sum, gsl_matrix *trans_sum, int multitree) {
	int i, j;
	for(i=0; i<6; i++) {
		gsl_vector_set(ancestral_sum[0], i, gsl_vector_get(ancestral_sum[0], i) + root[0]->dist[i]);
		if(multitree) {
			for(j=1; j<6; j++) {
				gsl_vector_set(ancestral_sum[j], i, gsl_vector_get(ancestral_sum[j], i) + root[j]->dist[i]);
			}
		}
		gsl_vector_set(stabs_sum, i, gsl_vector_get(stabs_sum, i) + gsl_vector_get(stabs, i));
		for(j=0; j<6; j++) {
			gsl_matrix_set(trans_sum, i, j, gsl_matrix_get(trans_sum, i, j) + gsl_matrix_get(trans, i, j));
		}
	}
}

void normalise_samples(gsl_vector **ancestral_sum, gsl_vector *stabs_sum, gsl_matrix *trans_sum, int samples, int multitree) {
	int i, j;
	for(i=0; i<6; i++) {
		gsl_vector_set(ancestral_sum[0], i, gsl_vector_get(ancestral_sum[0], i) / samples);
		if(multitree) for(j=1; j<6; j++) gsl_vector_set(ancestral_sum[j], i, gsl_vector_get(ancestral_sum[j], i) / samples);
		gsl_vector_set(stabs_sum, i, gsl_vector_get(stabs_sum, i) / samples);
		for(j=0; j<6; j++) gsl_matrix_set(trans_sum, i, j, gsl_matrix_get(trans_sum, i, j) / samples);
	}
}
