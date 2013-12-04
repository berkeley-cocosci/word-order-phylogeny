#ifndef GSLWS_H
#define GSLWS_H
 
#include<gsl/gsl_eigen.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_matrix_complex_float.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_vector_complex.h>
#include<gsl/gsl_permutation.h>

struct gslws {
	gsl_vector_complex *evals;
        gsl_matrix_complex *evecs;
        gsl_matrix_complex *evecs_inv;
        gsl_matrix *P;
	gsl_matrix_complex *c_P;
	gsl_matrix_complex *temp1;
	gsl_matrix_complex *temp2;
	gsl_matrix *Q_copy;
	gsl_matrix_complex *evecs_copy;
	gsl_permutation *perm;
	gsl_eigen_nonsymmv_workspace *eigenws;
	gsl_vector *tempvec1, *tempvec2;
};
typedef struct gslws gslws_t;

void alloc_gslws(gslws_t *ws);
void free_gslws(gslws_t *ws);

#endif /* GSLWS_H */
