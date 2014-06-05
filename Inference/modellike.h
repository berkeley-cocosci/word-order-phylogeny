#ifndef MODELLIKE_H
#define MODELLIKE_H

#include "gslworkspace.h"
#include "mcmc.h"

double get_log_prior(gsl_vector *stabs, gsl_matrix *trans);
double get_model_loglh(node_t *tree, gsl_matrix *Q, gslws_t *ws);
double get_balanced_model_loglh(node_t *tree, mcmc_t *mcmc, gslws_t *ws);
double dummy_tree_log_lh(mcmc_t *mcmc, gslws_t *ws);

#endif /* MODELLIKE_H */
