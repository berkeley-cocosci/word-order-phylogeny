#ifndef MCMC_H
#define MCMC_H
 
void initialise_stabs(gsl_vector *stabs, double x);
void initialise_trans(gsl_rng *r, gsl_matrix *trans);
void build_q(gsl_matrix *Q, gsl_vector *stabs, gsl_matrix *trans);
double get_prior(gsl_vector *stabs, gsl_matrix *trans);
void draw_proposal(FILE *fp, gsl_rng *r, gsl_vector *stabs, gsl_vector *stabs_dash, gsl_matrix *trans, gsl_matrix *trans_dash);
double mcmc_iteration(FILE *fp, gsl_rng* r, node_t **trees, gsl_vector *stabs, gsl_vector *stabs_dash, gsl_matrix *trans, gsl_matrix *trans_dash, double old_posterior, int multitree);

#endif /* MCMC_H */
