#ifndef MODELLIKE_H
#define MODELLIKE_H

#include "gslworkspace.h"

double get_log_prior(gsl_vector *stabs, gsl_matrix *trans);
double get_model_loglh(node_t *tree, gsl_matrix *Q, gslws_t *ws);

#endif /* MODELLIKE_H */
