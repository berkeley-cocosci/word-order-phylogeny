#ifndef MATRIX_H
#define MATRIX_H
 
void fprint_vector(FILE *fp, gsl_vector *v);
void fprint_matrix(FILE *fp, gsl_matrix *m);
void vector_copy(gsl_vector *in, gsl_vector *out);
void matrix_copy(gsl_matrix *in, gsl_matrix *out);
void compute_p(gsl_vector_complex *evals, gsl_matrix_complex *evecs, gsl_matrix_complex *evecs_inv, double t, gsl_matrix *out);
void decompose_q(gsl_matrix *Q, gsl_vector_complex *evals, gsl_matrix_complex *evecs, gsl_matrix_complex *evecs_inv);

#endif /* MATRIX_H */
