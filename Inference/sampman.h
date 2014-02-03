#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include "gslworkspace.h"
#include "tree.h"

#define SOV_TO_SVO 0
#define SVO_TO_SOV 1
#define VSO_TO_SOV 2
#define SOV_MOST_STAB 3
#define SVO_MOST_STAB 4
#define VSO_MOST_STAB 5

struct sampman{
	int sample_count;
	int stabs_log_pointer;
	int log_pointer;
	int stability_sampling_rate;
	int multitree;
	FILE *samplesfp, *ancestralsfp;
	double max_log_lh;
	double max_log_poster;
	double prior_log[100000];
	double likelihood_log[100000];
	double posterior_log[100000];
	gsl_vector *stabs_sum;
	gsl_vector *stabs_map;
	gsl_vector *stabs_log[500];
	gsl_matrix *trans;
	gsl_matrix *trans_sum;
	gsl_matrix *trans_map;
	gsl_matrix *Q_sum;
	gsl_matrix *Q_map;
	gsl_matrix *P_sum;
	gsl_matrix *P_map;
	gsl_vector *stationary_sum;
	gsl_vector *stationary_map;
	gsl_vector *ancestral_sum;
	gsl_vector *ancestral_map;
	gsl_vector *sliding_prior_ancestral_sum[100];
	gsl_vector *sliding_prior_ancestral_map[100];
	gsl_vector *stationary_prior_ancestral_sum;
	gsl_vector *stationary_prior_ancestral_map;
	gsl_vector *evidence_sum;
	float statistics[20];
	char families[][16];
};
typedef struct sampman sampman_t;

void initialise_sampman(sampman_t *sm, char *outdir);
void compute_means(sampman_t *sm);
void reset_sampman(sampman_t *sm);
void process_sample(sampman_t *sm, mcmc_t *mcmc, gslws_t *ws, node_t *tree);
void finish(sampman_t *sm);
void save_common_q(char *directory, sampman_t *sms);
void save_indiv_q(char *directory, sampman_t *sm);
