#ifndef MCMC_H
#define MCMC_H
 
#include <gsl/gsl_rng.h>
#include "gslworkspace.h"
#include "tree.h"

struct mcmc {
	double log_prior;
	double log_lh;
	double log_poster;
	gsl_vector *stabs;
	gsl_vector *stabs_dash;
	gsl_vector *stabs_max;
	gsl_matrix *trans;
	gsl_matrix *trans_dash;
	gsl_matrix *trans_max;
	gsl_matrix *Q;
	gsl_matrix *Q_dash;
	gsl_rng *r;
};
typedef struct mcmc mcmc_t;

void initialise_mcmc(mcmc_t *mcmc);
void compute_single_tree_probabilities(mcmc_t *mcmc, node_t *tree, gslws_t *ws);
void compute_multi_tree_probabilities(mcmc_t *mcmc, node_t **trees, gslws_t *wes);
void random_restart(mcmc_t *mcmc);
void initialise_stabs(mcmc_t *mcmc, double x);
void initialise_trans(mcmc_t *mcmc);
void build_q(mcmc_t *mcmc);
void build_q_dash(mcmc_t *mcmc);
void swap_column_step(mcmc_t *mcmc);
void swap_row_step(mcmc_t *mcmc);
void add_step(mcmc_t *mcmc, double variance);
void draw_proposal(mcmc_t *mcmc);
void handle_acceptance_probability(double *a);
void single_tree_mcmc_iteration(FILE *fp, mcmc_t *mcmc, node_t *tree, gslws_t *ws);
void multi_tree_mcmc_iteration(FILE *fp, mcmc_t *mcmc, node_t **trees, gslws_t *wses);

#endif /* MCMC_H */
