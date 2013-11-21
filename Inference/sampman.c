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
	sm->ancestral_sum = gsl_vector_alloc(6);
	sm->ancestral_map = gsl_vector_alloc(6);
	sm->fuzz_prior_ancestral_sum = gsl_matrix_alloc(6, 6);
	sm->fuzz_prior_ancestral_map = gsl_matrix_alloc(6, 6);
	sm->stationary_prior_ancestral_sum = gsl_matrix_alloc(6, 6);
	sm->stationary_prior_ancestral_map = gsl_matrix_alloc(6, 6);
	sm->sliding_fuzz_ancestral_sum = calloc(100, sizeof(gsl_matrix*));
	for(i=0;i<100;i++) sm->sliding_fuzz_ancestral_sum[i] = gsl_matrix_alloc(6, 6);
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
	gsl_vector_set_zero(sm->ancestral_sum);
	gsl_vector_set_zero(sm->ancestral_map);
	gsl_matrix_set_zero(sm->fuzz_prior_ancestral_sum);
	gsl_matrix_set_zero(sm->fuzz_prior_ancestral_map);
	gsl_matrix_set_zero(sm->stationary_prior_ancestral_sum);
	gsl_matrix_set_zero(sm->stationary_prior_ancestral_map);
	for(i=0;i<100;i++) gsl_matrix_set_zero(sm->sliding_fuzz_ancestral_sum[i]);
	for(i=0;i<10;i++) sm->statistics[i] = 0.0;
	sm->sample_count = 0;
	sm->max_log_poster = -1000000000;
}

void compute_means(sampman_t *sm) {
	uint8_t i;
	gsl_vector_scale(sm->stabs_sum, 1.0 / sm->sample_count);
	gsl_vector_scale(sm->stationary_sum, 1.0 / sm->sample_count);
	gsl_matrix_scale(sm->trans_sum, 1.0 / sm->sample_count);
	gsl_matrix_scale(sm->Q_sum, 1.0 / sm->sample_count);
	gsl_matrix_scale(sm->P_sum, 1.0 / sm->sample_count);
	gsl_vector_scale(sm->ancestral_sum, 1.0 / sm->sample_count);
	gsl_matrix_scale(sm->fuzz_prior_ancestral_sum, 1.0 / sm->sample_count);
	gsl_matrix_scale(sm->stationary_prior_ancestral_sum, 1.0 / sm->sample_count);
	for(i=0;i<100;i++) gsl_matrix_scale(sm->sliding_fuzz_ancestral_sum[i], 1.0 / sm->sample_count);
	for(i=0;i<10;i++) sm->statistics[i] /= sm->sample_count;
}

void process_sample(sampman_t *sm, mcmc_t *mcmc, node_t *tree) {
	uint8_t i;
	sm->sample_count++;
	gsl_vector_add(sm->stabs_sum, mcmc->stabs);
	gsl_matrix_add(sm->trans_sum, mcmc->trans);
	gsl_matrix_add(sm->Q_sum, mcmc->Q);
	for(i=0;i<6;i++) gsl_vector_set(sm->ancestral_sum, i, gsl_vector_get(sm->ancestral_sum, i) + tree->dist[i]);
	if(mcmc->log_poster > sm->max_log_poster) {
		sm->max_log_poster = mcmc->log_poster;
		gsl_vector_memcpy(sm->stabs_map, mcmc->stabs);
		gsl_matrix_memcpy(sm->trans_map, mcmc->trans);
		gsl_matrix_memcpy(sm->Q_map, mcmc->Q);
		gsl_matrix_add(sm->Q_sum, mcmc->Q);
		for(i=0;i<6;i++) gsl_vector_set(sm->ancestral_map, i, tree->dist[i]);
	}
}

void finish(sampman_t *sm) {
	printf("Mean Q is:\n");
	fprint_matrix(stdout, sm->Q_sum);
}

void save_common_q(char *directory, sampman_t *sm) {
	char filename[1024];
	FILE *fp;
	/* Write summary file */
	strcpy(filename, directory);
	strcat(filename, "/summary");
	fp = fopen(filename, "w");

	fprintf(fp, "Posterior mean stabilities:\n");
	fprint_vector(fp, sm->stabs_sum);
	fprintf(fp, "----------\n");
	fprintf(fp, "Posterior mean transitions:\n");
	fprint_matrix(fp, sm->trans_sum);
	fprintf(fp, "----------\n");
	fprintf(fp, "Posterior mean Q:\n");
	fprint_matrix(fp, sm->Q_sum);
	fprintf(fp, "----------\n");
	fprintf(fp, "P matrix over short branch:\n");
	fprint_matrix(fp, sm->P_sum);
	fprintf(fp, "----------\n");
	fprintf(fp, "Posterior mean stationary:\n");
	fprint_vector(fp, sm->stationary_sum);
	fprintf(fp, "----------\n");
	fprintf(fp, "Posterior mean uniform ancestrals:\n");
	fprint_vector(fp, sm->ancestral_sum);
	fprintf(fp, "----------\n");
	fprintf(fp, "Posterior mean fuzzy ancestrals:\n");
	fprint_matrix(fp, sm->fuzz_prior_ancestral_sum);
	fprintf(fp, "----------\n");
	fprintf(fp, "Posterior mean stationary ancestrals:\n");
	fprint_matrix(fp, sm->stationary_prior_ancestral_sum);
	fprintf(fp, "----------\n");
	fprintf(fp, "Hypothesis probabilities:\n");
	fprintf(fp, "SOV to SVO (over VSO): %f\n", sm->statistics[SOV_TO_SVO]);
	fprintf(fp, "SVO to SOV (over VSO): %f\n", sm->statistics[SVO_TO_SOV]);
	fprintf(fp, "VSO to SOV (over SVO): %f\n", sm->statistics[VSO_TO_SOV]);
	fprintf(fp, "SVO most stable: %f\n", sm->statistics[SVO_MOST_STAB]);
//	fprintf(fp, "SOV most likely ancestor: %f\n", statistics[SOV_MOST_LIKE]);
	fprintf(fp, "----------\n");
	fprintf(fp, "Maximum posterior: %f\n", sm->max_log_poster);
	fprintf(fp, "----------\n");
	fprintf(fp, "MAP stabilities:\n");
	fprint_vector(fp, sm->stabs_map);
	fprintf(fp, "MAP transitions:\n");
	fprint_matrix(fp, sm->trans_map);
	fprintf(fp, "MAP stationary:\n");
	fprint_vector(fp, sm->stationary_map);
	fprintf(fp, "MAP ancestral word order distributions:\n");
	fprint_vector(fp, sm->ancestral_sum);
	fprintf(fp, "----------\n");
	fclose(fp);
}

void save_indiv_q(char *directory, sampman_t **sm) {
	char filename[1024];
	FILE *fp;
	int i;
	/* Write summary file */
	strcpy(filename, directory);
	strcat(filename, "/summary");
	fp = fopen(filename, "w");

	/* This is all basically nonsense which makes no sense in an
	 * individual Q situation.  However, for backward compatibility
	 * with the code that parses up this output and turns it into
	 * MATLAB stuff, we need to generate it anyway.
	 */

	fprintf(fp, "Posterior mean stabilities:\n");
	fprint_vector(fp, sm[0]->stabs_sum);
	fprintf(fp, "----------\n");
	fprintf(fp, "Posterior mean transitions:\n");
	fprint_matrix(fp, sm[0]->trans_sum);
	fprintf(fp, "----------\n");
	fprintf(fp, "Posterior mean Q:\n");
	fprint_matrix(fp, sm[0]->Q_sum);
	fprintf(fp, "----------\n");
	fprintf(fp, "P matrix over short branch:\n");
	fprint_matrix(fp, sm[0]->P_sum);
	fprintf(fp, "----------\n");
	fprintf(fp, "Posterior mean stationary:\n");
	fprint_vector(fp, sm[0]->stationary_sum);
	fprintf(fp, "----------\n");
	fprintf(fp, "Posterior mean ancestral:\n");
	for(i=0; i<6; i++) fprint_vector(fp, sm[i]->ancestral_sum);
	fprintf(fp, "\n");
	fprintf(fp, "----------\n");
	fprintf(fp, "Posterior mean fuzzy ancestrals:\n");
	for(i=0; i<6; i++) fprint_matrix(fp, sm[i]->fuzz_prior_ancestral_sum);
	fprintf(fp, "\n");
	fprintf(fp, "----------\n");
	fprintf(fp, "Posterior mean stationary ancestrals:\n");
	for(i=0; i<6; i++) fprint_matrix(fp, sm[i]->stationary_prior_ancestral_sum);
	fprintf(fp, "\n");
	fprintf(fp, "----------\n");
	fprintf(fp, "Hypothesis probabilities:\n");
	fprintf(fp, "SOV to SVO (over VSO): %f\n", sm[0]->statistics[SOV_TO_SVO]);
	fprintf(fp, "SVO to SOV (over VSO): %f\n", sm[0]->statistics[SVO_TO_SOV]);
	fprintf(fp, "VSO to SOV (over SVO): %f\n", sm[0]->statistics[VSO_TO_SOV]);
	fprintf(fp, "SVO most stable: %f\n", sm[0]->statistics[SVO_MOST_STAB]);
//	fprintf(fp, "SOV most likely ancestor: %f\n", statistics[SOV_MOST_LIKE]);
	fprintf(fp, "----------\n");
	fprintf(fp, "Maximum posterior: %f\n", sm[0]->max_log_poster);
	fprintf(fp, "----------\n");
	fprintf(fp, "MAP stabilities:\n");
	fprint_vector(fp, sm[0]->stabs_map);
	fprintf(fp, "MAP transitions:\n");
	fprint_matrix(fp, sm[0]->trans_map);
	fprintf(fp, "MAP stationary:\n");
	fprint_vector(fp, sm[0]->stationary_map);
	fprintf(fp, "MAP ancestral word order distribution:\n");
	for(i=1;i<6;i++) fprint_vector(fp, sm[0]->ancestral_sum);
	fprintf(fp, "\n");
	fprintf(fp, "----------\n");

	fclose(fp);
}
