#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "params.h"
#include "matrix.h"
#include "mcmc.h"
#include "sampman.h"
#include "tree.h"

void initialise_sampman(sampman_t *sm, char *outdir) {
	int i;
	sm->sample_count = 0;
	sm->stabs_log_pointer = 0;
	sm->log_pointer = 0;
	sm->max_log_lh = -1000000000;
	sm->max_log_poster = -1000000000;
	sm->stabs_sum = gsl_vector_alloc(6);
	sm->stabs_map = gsl_vector_alloc(6);
	for(i=0; i<500; i++) sm->stabs_log[i] = gsl_vector_alloc(6);
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
	sm->stationary_prior_ancestral_sum = gsl_vector_alloc(6);
	sm->stationary_prior_ancestral_map = gsl_vector_alloc(6);
	for(i=0; i<100; i++) sm->sliding_prior_ancestral_sum[i] = gsl_vector_alloc(6);
//	sm->sliding_prior_ancestral_map = calloc(100, sizeof(gsl_matrix*));
	sm->evidence_sum = gsl_vector_alloc(6);
}

void reset_sampman(sampman_t *sm) {
	uint8_t i;
	sm->stabs_log_pointer = 0;
	sm->log_pointer = 0;
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
	gsl_vector_set_zero(sm->stationary_prior_ancestral_sum);
	gsl_vector_set_zero(sm->stationary_prior_ancestral_map);
	for(i=0;i<10;i++) sm->statistics[i] = 0.0;
	sm->sample_count = 0;
	sm->max_log_lh = -1000000000;
	sm->max_log_poster = -1000000000;
	for(i=0; i<100; i++) gsl_vector_set_zero(sm->sliding_prior_ancestral_sum[i]);
	//for(i=0; i<100; i++) gsl_vector_set_zero(sm->sliding_prior_ancestral_map[i]);
	gsl_vector_set_zero(sm->evidence_sum);
}

void compute_means(sampman_t *sm) {
	uint8_t i;
	gsl_vector_scale(sm->stabs_sum, 1.0 / sm->sample_count);
	gsl_vector_scale(sm->stationary_sum, 1.0 / sm->sample_count);
	gsl_matrix_scale(sm->trans_sum, 1.0 / sm->sample_count);
	gsl_matrix_scale(sm->Q_sum, 1.0 / sm->sample_count);
	gsl_matrix_scale(sm->P_sum, 1.0 / sm->sample_count);
	gsl_vector_scale(sm->ancestral_sum, 1.0 / sm->sample_count);
	for(i=0; i<100; i++) gsl_vector_scale(sm->sliding_prior_ancestral_sum[i], 1.0 / sm->sample_count);
	fprint_vector(stdout, sm->stationary_prior_ancestral_sum);
	printf("%d\n", sm->sample_count);
	gsl_vector_scale(sm->stationary_prior_ancestral_sum, 1.0 / sm->sample_count);
	fprint_vector(stdout, sm->stationary_prior_ancestral_sum);
	gsl_vector_scale(sm->evidence_sum, 1.0 / sm->sample_count);
	for(i=0;i<10;i++) sm->statistics[i] /= sm->sample_count;

}

void process_sample(sampman_t *sm, mcmc_t *mcmc, gslws_t *ws, node_t *tree) {
	uint8_t i, j, k;
	double years, t, norm;
	// Do simple stuff
	sm->sample_count++;
	gsl_vector_add(sm->stabs_sum, mcmc->stabs);
	gsl_matrix_add(sm->trans_sum, mcmc->trans);
	gsl_matrix_add(sm->Q_sum, mcmc->Q);
	for(i=0;i<6;i++) gsl_vector_set(sm->ancestral_sum, i, gsl_vector_get(sm->ancestral_sum, i) + tree->dist[i]);
	// Short P
	decompose_q(mcmc->Q, ws);
	compute_p(ws, 0.1);
	gsl_matrix_add(sm->P_sum, ws->P);
	// Stationary
	compute_p(ws, 100000000.0);
	for(i=0;i<6;i++) gsl_vector_set(sm->stationary_sum, i, gsl_vector_get(sm->stationary_sum, i) + gsl_matrix_get(ws->P, 0, i));
	// Update maxima
	if(mcmc->log_poster > sm->max_log_poster) {
		sm->max_log_poster = mcmc->log_poster;
		gsl_vector_memcpy(sm->stabs_map, mcmc->stabs);
		gsl_matrix_memcpy(sm->trans_map, mcmc->trans);
		gsl_matrix_memcpy(sm->Q_map, mcmc->Q);
		for(i=0;i<6;i++) gsl_vector_set(sm->ancestral_map, i, tree->dist[i]);
		for(i=0;i<6;i++) gsl_vector_set(sm->stationary_map, i, gsl_matrix_get(ws->P, 0, i));
		compute_p(ws, 0.1);
		gsl_matrix_memcpy(sm->P_map, ws->P);
	}
	if(mcmc->log_lh > sm->max_log_lh) {
		sm->max_log_lh = mcmc->log_lh;
	}
	// Take stability samples at appropriate rate
	// (Handle the case of stability_sampling_rate = 0 because sampling rate is
	// computed from tree and sample counts, and is sometimes zero for very small tasks)
	if(sm->stability_sampling_rate > 0) {
		if(sm->sample_count % sm->stability_sampling_rate == 0) {
			gsl_vector_memcpy(sm->stabs_log[sm->stabs_log_pointer], mcmc->stabs);
			sm->stabs_log_pointer++;
		}
	}
	if(sm->log_pointer < 100000) {
		sm->prior_log[sm->log_pointer] = mcmc->log_prior;
		sm->likelihood_log[sm->log_pointer] = mcmc->log_lh;
		sm->posterior_log[sm->log_pointer] = mcmc->log_poster;
		sm->log_pointer++;
	}

	// Stationary prior
	gsl_vector_set_zero(ws->tempvec1);
	for(i=0;i<6;i++) gsl_vector_set(ws->tempvec1, i, gsl_matrix_get(ws->P, 0, i) * tree->dist[i]);
	norm = 0.0;
	for(j=0;j<6;j++) norm += gsl_vector_get(ws->tempvec1, j);
	gsl_vector_scale(ws->tempvec1, 1.0/norm);
	gsl_vector_add(sm->stationary_prior_ancestral_sum, ws->tempvec1);

	// Sliding prior
	for(i=0;i<100;i++) {
		// Compute prior
		years = i*40000/99.0;
		t = years / 10000.0;
		compute_p(ws, t);
		gsl_vector_set_zero(ws->tempvec1);
		for(j=0; j<6; j++) {
			for(k=0; k<6; k++) {
				gsl_vector_set(ws->tempvec1, k, gsl_vector_get(ws->tempvec1, k) + gsl_matrix_get(ws->P, j, k));
			}
		}
		norm = 0.0;
		for(j=0;j<6;j++) norm += gsl_vector_get(ws->tempvec1, j);
		gsl_vector_scale(ws->tempvec1, 1.0/norm);

		// Multiply likelihood by prior
		gsl_vector_set_zero(ws->tempvec2);
		for(j=0; j<6; j++) gsl_vector_set(ws->tempvec2, j, gsl_vector_get(ws->tempvec1, j) * tree->dist[j]);
		norm = 0;
		for(j=0; j<6; j++) norm += gsl_vector_get(ws->tempvec2, j);
		gsl_vector_scale(ws->tempvec2, 1.0 / norm);
		gsl_vector_add(sm->sliding_prior_ancestral_sum[i], ws->tempvec2);
	}
	// Statistics
	if(gsl_matrix_get(mcmc->trans, 0, 1) > gsl_matrix_get(mcmc->trans, 0, 2)) sm->statistics[SOV_TO_SVO]++;
	if(gsl_matrix_get(mcmc->trans, 1, 0) > gsl_matrix_get(mcmc->trans, 1, 2)) sm->statistics[SVO_TO_SOV]++;
	if(gsl_matrix_get(mcmc->trans, 2, 0) > gsl_matrix_get(mcmc->trans, 2, 1)) sm->statistics[VSO_TO_SOV]++;
	if(gsl_vector_min_index(mcmc->stabs) == 0)  sm->statistics[SOV_MOST_STAB]++;
	if(gsl_vector_min_index(mcmc->stabs) == 1)  sm->statistics[SVO_MOST_STAB]++;
	if(gsl_vector_min_index(mcmc->stabs) == 2)  sm->statistics[VSO_MOST_STAB]++;
	// Evidence
	gsl_vector_set_zero(ws->tempvec1);
	for(i=0; i<6; i++) {
		norm = 0;
		for(j=0; j<6; j++) {
			if(j != i) norm += tree->dist[j]*(1.0/5.0);
		}
		gsl_vector_set(ws->tempvec1, i, log(tree->dist[i] / norm));
		gsl_vector_add(sm->evidence_sum, ws->tempvec1);
	}

}

void process_prior_sample(sampman_t *sm, mcmc_t *mcmc, gslws_t *ws) {
	if(sm->log_pointer < 200000) {
		sm->prior_log[sm->log_pointer] = mcmc->log_prior;
		sm->likelihood_log[sm->log_pointer] = mcmc->log_lh;
		sm->posterior_log[sm->log_pointer] = mcmc->log_poster;
		sm->log_pointer++;
	}
}

void finish(sampman_t *sm) {
	printf("Mean Q is:\n");
	fprint_matrix(stdout, sm->Q_sum);
}

void save_indiv_q(char *directory, sampman_t *sm) {
	char filename[1024];
	FILE *fp;
	char mkcmd[1024];
	int i;
	/* Make directory */
	sprintf(mkcmd, "mkdir -p %s", directory);
	system(mkcmd);
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
	fprint_vector(fp, sm->sliding_prior_ancestral_sum[0]);
	fprintf(fp, "----------\n");
	fprintf(fp, "Posterior mean stationary ancestrals:\n");
	fprint_vector(fp, sm->stationary_prior_ancestral_sum);
	fprintf(fp, "----------\n");
	fprintf(fp, "Posterior mean ancestry evidence:\n");
	fprint_vector(fp, sm->evidence_sum);
	fprintf(fp, "----------\n");
	fprintf(fp, "Hypothesis probabilities:\n");
	fprintf(fp, "SOV to SVO (over VSO): %f\n", sm->statistics[SOV_TO_SVO]);
	fprintf(fp, "SVO to SOV (over VSO): %f\n", sm->statistics[SVO_TO_SOV]);
	fprintf(fp, "VSO to SOV (over SVO): %f\n", sm->statistics[VSO_TO_SOV]);
	fprintf(fp, "SOV most stable: %f\n", sm->statistics[SOV_MOST_STAB]);
	fprintf(fp, "SVO most stable: %f\n", sm->statistics[SVO_MOST_STAB]);
	fprintf(fp, "VSO most stable: %f\n", sm->statistics[VSO_MOST_STAB]);
//	fprintf(fp, "SOV most likely ancestor: %f\n", statistics[SOV_MOST_LIKE]);
	fprintf(fp, "----------\n");
	fprintf(fp, "Maximum log likelihood: %f\n", sm->max_log_lh);
	fprintf(fp, "Maximum log posterior: %f\n", sm->max_log_poster);
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
	
	/* Save logged stability vectors */
	strcpy(filename, directory);
	strcat(filename, "/stabilities");
	fp = fopen(filename, "w");
	for(i=0; i<500; i++) fprint_vector(fp, sm->stabs_log[i]);
	fclose(fp);

	/* Save logged likelihoods */
	strcpy(filename, directory);
	strcat(filename, "/bayes_factor_posterior_samples");
	fp = fopen(filename, "w");
	for(i=0; i<sm->log_pointer; i++) fprintf(fp, "%f, %f, %f\n", sm->prior_log[i], sm->likelihood_log[i], sm->posterior_log[i]);
	fclose(fp);

	/* Save sliding prior ancestral distributions */
	strcpy(filename, directory);
	strcat(filename, "/sliding_prior");
	fp = fopen(filename, "w");
	for(i=0; i<100; i++) fprint_vector(fp, sm->sliding_prior_ancestral_sum[i]);
	fclose(fp);
}

void save_common_q(char *directory, sampman_t *sms) {
	char filename[1024];
	FILE *fp;
	int i, j;
	char mkcmd[1024];
	/* Make directory */
	sprintf(mkcmd, "mkdir -p %s", directory);
	system(mkcmd);
	/* Write summary file */
	strcpy(filename, directory);
	strcat(filename, "/summary");
	fp = fopen(filename, "w");

	/* Nonsense */
	fprintf(fp, "Posterior mean stabilities:\n");
	fprint_vector(fp, sms[0].stabs_sum);
	fprintf(fp, "----------\n");
	fprintf(fp, "Posterior mean transitions:\n");
	fprint_matrix(fp, sms[0].trans_sum);
	fprintf(fp, "----------\n");
	fprintf(fp, "Posterior mean Q:\n");
	fprint_matrix(fp, sms[0].Q_sum);
	fprintf(fp, "----------\n");
	fprintf(fp, "P matrix over short branch:\n");
	fprint_matrix(fp, sms[0].P_sum);
	fprintf(fp, "----------\n");
	fprintf(fp, "Posterior mean stationary:\n");
	fprint_vector(fp, sms[0].stationary_sum);
	fprintf(fp, "----------\n");

	fprintf(fp, "Posterior mean uniform ancestrals:\n");
	for(i=0; i<NUM_TREES; i++) fprint_vector(fp, sms[i].ancestral_sum);
	fprintf(fp, "----------\n");
	fprintf(fp, "Posterior mean fuzzy ancestrals:\n");
	for(i=0; i<NUM_TREES; i++) fprint_vector(fp, sms[i].sliding_prior_ancestral_sum[0]);
	fprintf(fp, "----------\n");
	fprintf(fp, "Posterior mean stationary ancestrals:\n");
	for(i=0; i<NUM_TREES; i++) fprint_vector(fp, sms[i].stationary_prior_ancestral_sum);
	fprintf(fp, "----------\n");
	fprintf(fp, "Posterior mean ancestry evidence:\n");
	for(i=0; i<NUM_TREES; i++) fprint_vector(fp, sms[i].evidence_sum);
	fprintf(fp, "----------\n");

	fprintf(fp, "Hypothesis probabilities:\n");
	fprintf(fp, "SOV to SVO (over VSO): %f\n", sms[0].statistics[SOV_TO_SVO]);
	fprintf(fp, "SVO to SOV (over VSO): %f\n", sms[0].statistics[SVO_TO_SOV]);
	fprintf(fp, "VSO to SOV (over SVO): %f\n", sms[0].statistics[VSO_TO_SOV]);
	fprintf(fp, "SOV most stable: %f\n", sms[0].statistics[SOV_MOST_STAB]);
	fprintf(fp, "SVO most stable: %f\n", sms[0].statistics[SVO_MOST_STAB]);
	fprintf(fp, "VSO most stable: %f\n", sms[0].statistics[VSO_MOST_STAB]);
//	fprintf(fp, "SOV most likely ancestor: %f\n", statistics[VSO_MOST_LIKE]);
	fprintf(fp, "----------\n");
	fprintf(fp, "Maximum log likelihood: %f\n", sms[0].max_log_lh);
	fprintf(fp, "Maximum log posteror: %f\n", sms[0].max_log_poster);
	fprintf(fp, "----------\n");
	fprintf(fp, "MAP stabilities:\n");
	fprint_vector(fp, sms[0].stabs_map);
	fprintf(fp, "MAP transitions:\n");
	fprint_matrix(fp, sms[0].trans_map);
	fprintf(fp, "MAP stationary:\n");
	fprint_vector(fp, sms[0].stationary_map);

	fprintf(fp, "MAP ancestral word order distributions:\n");
	fprint_vector(fp, sms[0].ancestral_map);
	fprintf(fp, "----------\n");
	fclose(fp);

	/* Save logged stability vectors */
	strcpy(filename, directory);
	strcat(filename, "/stabilities");
	fp = fopen(filename, "w");
	for(i=0; i<500; i++) fprint_vector(fp, sms[0].stabs_log[i]);
	fclose(fp);

	/* Save logged likelihoods */
	strcpy(filename, directory);
	strcat(filename, "/bayes_factor_posterior_samples");
	fp = fopen(filename, "w");
	for(i=0; i<sms[0].log_pointer; i++) {
		fprintf(fp, "%f, %f, %f\n", sms[0].prior_log[i], sms[0].likelihood_log[i], sms[0].posterior_log[i]);
	}
	fclose(fp);

	/* Save sliding prior ancestral distributions */
	strcpy(filename, directory);
	strcat(filename, "/sliding_prior");
	fp = fopen(filename, "w");
	for(i=0; i<100; i++) {
		for(j=0; j<NUM_TREES; j++) {
			fprint_vector(fp, sms[j].sliding_prior_ancestral_sum[i]);
		}
		fprintf(fp, "----------\n");
	}
	fclose(fp);

}
