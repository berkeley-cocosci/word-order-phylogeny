#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>

#include "matrix.h"
#include "mcmc.h"
#include "sampman.h"
#include "tree.h"

void initialise_sampman(sampman_t *sm, char *outdir) {
	uint8_t i;
	char mkcmd[1024];
	sprintf(mkcmd, "mkdir -p %s", outdir);
	system(mkcmd);

	sm->stabs_sum = gsl_vector_alloc(6);
	sm->stabs_map = gsl_vector_alloc(6);
	sm->trans = gsl_matrix_alloc(6, 6);
	sm->trans_sum = gsl_matrix_alloc(6, 6);
	sm->trans_map = gsl_matrix_alloc(6, 6);
	sm->Q_sum = gsl_matrix_alloc(6, 6);
	sm->Q_map = gsl_matrix_alloc(6, 6);
	sm->P_sum = gsl_matrix_alloc(6, 6);
	sm->P_map = gsl_matrix_alloc(6, 6);
	sm->stationary_sum = gsl_vector_alloc(6);
	sm->stationary_map = gsl_vector_alloc(6);
	sm->ancestral_sum = gsl_matrix_alloc(6, 6);
	sm->ancestral_map = gsl_matrix_alloc(6, 6);
	sm->fuzz_prior_ancestral_sum = gsl_matrix_alloc(6, 6);
	sm->fuzz_prior_ancestral_map = gsl_matrix_alloc(6, 6);
	sm->stationary_prior_ancestral_sum = gsl_matrix_alloc(6, 6);
	sm->stationary_prior_ancestral_map = gsl_matrix_alloc(6, 6);
	sm->sliding_fuzz_ancestral_sum = calloc(100, sizeof(gsl_matrix*));
	for(i=0;i<100;i++) sm->sliding_fuzz_ancestral_sum[i] = gsl_matrix_alloc(6, 6);
}

void compute_means(sampman_t *sm) {
	uint8_t i;
	gsl_vector_scale(sm->stabs_sum, 1.0 / sm->sample_count);
	gsl_vector_scale(sm->stationary_sum, 1.0 / sm->sample_count);
	gsl_matrix_scale(sm->trans_sum, 1.0 / sm->sample_count);
	gsl_matrix_scale(sm->Q_sum, 1.0 / sm->sample_count);
	gsl_matrix_scale(sm->P_sum, 1.0 / sm->sample_count);
	gsl_matrix_scale(sm->ancestral_sum, 1.0 / sm->sample_count);
	gsl_matrix_scale(sm->fuzz_prior_ancestral_sum, 1.0 / sm->sample_count);
	gsl_matrix_scale(sm->stationary_prior_ancestral_sum, 1.0 / sm->sample_count);
	for(i=0;i<100;i++) gsl_matrix_scale(sm->sliding_fuzz_ancestral_sum[i], 1.0 / sm->sample_count);
	for(i=0;i<10;i++) sm->statistics[i] /= sm->sample_count;
}

void reset_summary(sampman_t *sm) {
	uint8_t i;
	gsl_vector_set_zero(sm->stabs_sum);
	gsl_vector_set_zero(sm->stabs_map);
	gsl_vector_set_zero(sm->stationary_sum);
	gsl_vector_set_zero(sm->stationary_map);
	gsl_matrix_set_zero(sm->trans_sum);
	gsl_matrix_set_zero(sm->trans_map);
	gsl_matrix_set_zero(sm->Q_sum);
	gsl_matrix_set_zero(sm->Q_map);
	gsl_matrix_set_zero(sm->P_sum);
	gsl_matrix_set_zero(sm->P_map);
	gsl_matrix_set_zero(sm->ancestral_sum);
	gsl_matrix_set_zero(sm->ancestral_map);
	gsl_matrix_set_zero(sm->fuzz_prior_ancestral_sum);
	gsl_matrix_set_zero(sm->fuzz_prior_ancestral_map);
	gsl_matrix_set_zero(sm->stationary_prior_ancestral_sum);
	gsl_matrix_set_zero(sm->stationary_prior_ancestral_map);
	for(i=0;i<100;i++) gsl_matrix_set_zero(sm->sliding_fuzz_ancestral_sum[i]);
	for(i=0;i<10;i++) sm->statistics[i] = 0.0;
	sm->sample_count = 0;
	sm->max_log_poster = -1000000000;
}
void process_sample(sampman_t *sm, mcmc_t *mcmc, node_t **trees) {
	gsl_matrix_add(sm->Q_sum, mcmc->Q);
}

void shout_out(sampman_t *sm) {
	printf("Mean Q is:\n");
	fprint_matrix(stdout, sm->Q_sum);
}
