#include<stdio.h>
#include<math.h>

#include<gsl/gsl_blas.h>
#include<gsl/gsl_blas_types.h>
#include<gsl/gsl_complex.h>
#include<gsl/gsl_complex_math.h>
#include<gsl/gsl_eigen.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_matrix_complex_float.h>
#include<gsl/gsl_permutation.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_vector_complex.h>

#include "gslworkspace.h"

void fprint_vector(FILE *fp, gsl_vector *v) {
        int i;
        for(i=0; i<5; i++) {
                fprintf(fp, "%f ", gsl_vector_get(v, i));
        }
        fprintf(fp, "%f\n", gsl_vector_get(v, 5));
}

void fprint_matrix(FILE *fp, gsl_matrix *m) {
        int i, j;
        for(i=0; i<6; i++) {
                for(j=0; j<5; j++) {
                        fprintf(fp, "%f ", gsl_matrix_get(m, i, j));
                }
                fprintf(fp, "%f\n", gsl_matrix_get(m, i, 5));
        }
}

void initialise_q(gsl_matrix *Q, double q) {
	int i;
	gsl_matrix_set_all(Q, q);
	for(i=0; i<6; i++) {
		gsl_matrix_set(Q, i, i, -5.0*q);
	}
}

void realify_vector(gsl_vector_complex *complex, gsl_vector *real) {
	int i;
	for(i=0; i<6; i++) {
		gsl_vector_set(real, i, GSL_REAL(gsl_vector_complex_get(complex, i)));
	}
}

void realify_matrix(gsl_matrix_complex *complex, gsl_matrix *real) {
	int i;
	gsl_vector_complex *complex_row = gsl_vector_complex_alloc(6);
	gsl_vector real_row;
	for(i=0; i<6; i++) {
		gsl_matrix_complex_get_row(complex_row, complex, i);
		real_row = gsl_vector_complex_real(complex_row).vector;
		gsl_matrix_set_row(real, i, &real_row);
	}
	gsl_vector_complex_free(complex_row);
}

void invert(gsl_matrix_complex *in, gsl_matrix_complex *out, gsl_permutation *p, int *signum) {
	/* Why in the fuck is this not just part of gsl? */
	gsl_linalg_complex_LU_decomp(in, p, signum);
//	print_matrix(in);
	gsl_linalg_complex_LU_invert(in, p, out);

}

void eigensolve(gsl_matrix *Q, gslws_t *ws) {
	gsl_eigen_nonsymmv(Q, ws->evals, ws->evecs, ws->eigenws);
}

void compute_p(gslws_t *ws, double t) {
	int i;
	gsl_matrix_complex_set_zero(ws->temp1);
	gsl_matrix_complex_set_zero(ws->temp2);
	gsl_complex one, zero, z;
	GSL_SET_COMPLEX(&one, 1.0, 0.0);
	GSL_SET_COMPLEX(&zero, 0.0, 0.0);
	/* Fill temp1 with e^(t\lbamda_i) */
	for(i=0; i<6; i++) {
		z = gsl_complex_mul_real(gsl_vector_complex_get(ws->evals, i), t);
		gsl_matrix_complex_set(ws->temp1, i, i, gsl_complex_exp(z));
	}
	gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, one, ws->evecs, ws->temp1, zero, ws->temp2);
	/* Store temp2*evecs_inv in out */
	gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, one, ws->temp2, ws->evecs_inv, zero, ws->c_P);
	realify_matrix(ws->c_P, ws->P);
	//print_matrix(out);
//	printf("---\n");
}

void decompose_q(gsl_matrix *Q, gslws_t *ws) {
	int signum;
	int i,j;

	/* Make a copy of Q to give to the (destructive) eigensolver */
	gsl_matrix_memcpy(ws->Q_copy, Q);

	eigensolve(ws->Q_copy, ws);

	/* Make a copy of the eigenvectors to give to the (destructive) inverter */
	for(i=0; i<6; i++) {
		for(j=0; j<6; j++) {
			gsl_matrix_complex_set(ws->evecs_copy, i, j, gsl_matrix_complex_get(ws->evecs, i, j));
		}
	}
	
	/* Compute inverse eigenvector matrix */
	invert(ws->evecs_copy, ws->evecs_inv, ws->perm, &signum);

/*
	fp = fopen("eigendecomp.txt", "w");
	fprintf(fp, "Here's Q:\n");
	gsl_matrix_fprintf(fp, Q, "%f");
	fprintf(fp, "Here's the eigenvalues:\n");
	gsl_vector_fprintf(fp, evals, "%f");
	fprintf(fp, "Here're the eigenvectors:\n");
	gsl_matrix_fprintf(fp, evecs, "%f");
	fprintf(fp, "Here's the inverse of the eigenvector matrix\n");
	gsl_matrix_fprintf(fp, evecs_inv, "%f");
	fclose(fp);
*/
}
