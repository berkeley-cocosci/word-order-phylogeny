#include<gsl/gsl_matrix.h>
#include<gsl/gsl_matrix_complex_float.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_vector_complex.h>

#include "gslworkspace.h"

void alloc_gslws(gslws_t *ws) {
	ws->evals = gsl_vector_complex_alloc(6);
        ws->evecs = gsl_matrix_complex_alloc(6,6);
        ws->evecs_inv = gsl_matrix_complex_alloc(6,6);
        ws->P = gsl_matrix_alloc(6,6);
	ws->c_P = gsl_matrix_complex_alloc(6,6);
	ws->temp1 = gsl_matrix_complex_alloc(6,6);
	ws->temp2 = gsl_matrix_complex_alloc(6,6);
	ws->Q_copy = gsl_matrix_alloc(6,6);
	ws->evecs_copy = gsl_matrix_complex_alloc(6,6);
	ws->perm = gsl_permutation_alloc(6);
	ws->eigenws = gsl_eigen_nonsymmv_alloc(6);
	ws->tempvec1 = gsl_vector_alloc(6);
	ws->tempvec2 = gsl_vector_alloc(6);
}

void free_gslws(gslws_t *ws) {
	gsl_vector_complex_free(ws->evals);
	gsl_matrix_complex_free(ws->evecs);
	gsl_matrix_complex_free(ws->evecs_inv);
	gsl_matrix_free(ws->P);
	gsl_matrix_complex_free(ws->c_P);
	gsl_matrix_complex_free(ws->temp1);
	gsl_matrix_complex_free(ws->temp2);
	gsl_matrix_free(ws->Q_copy);
	gsl_matrix_complex_free(ws->evecs_copy);
	gsl_permutation_free(ws->perm);
	gsl_eigen_nonsymmv_free(ws->eigenws);
	gsl_vector_free(ws->tempvec1);
	gsl_vector_free(ws->tempvec2);
}
