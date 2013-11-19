#ifndef MATRIX_H
#define MATRIX_H

#include "gslworkspace.h"

void fprint_vector(FILE *fp, gsl_vector *v);
void fprint_matrix(FILE *fp, gsl_matrix *m);
void vector_copy(gsl_vector *in, gsl_vector *out);
void matrix_copy(gsl_matrix *in, gsl_matrix *out);
void compute_p(gslws_t *ws, double t);
void decompose_q(gsl_matrix *Q, gslws_t *ws);

#endif /* MATRIX_H */
