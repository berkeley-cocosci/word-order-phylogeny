#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include "tree.h"

#define SOV_TO_SVO 0
#define SVO_TO_SOV 1
#define VSO_TO_SOV 2
#define SVO_MOST_STAB 3
#define SOV_MOST_LIKE 4

struct sampman{
	int sample_count;
	int multitree;
	FILE *samplesfp, *ancestralsfp;
	double max_log_poster;
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
	gsl_vector *ancestral_sum;
	gsl_vector *ancestral_map;
	gsl_matrix *fuzz_prior_ancestral_sum;
	gsl_matrix *fuzz_prior_ancestral_map;
	gsl_matrix *stationary_prior_ancestral_sum;
	gsl_matrix *stationary_prior_ancestral_map;
	gsl_matrix **sliding_fuzz_ancestral_sum;
	float statistics[20];
	char families[][16];
};
typedef struct sampman sampman_t;

void initialise_sampman(sampman_t *sm, char *outdir);
void compute_means(sampman_t *sm);
void reset_summary(sampman_t *sm);
void process_sample(sampman_t *sm, mcmc_t *mcmc, node_t *tree);
void finish(sampman_t *sm);
void save_common_q(char *directory, sampman_t *sm);
