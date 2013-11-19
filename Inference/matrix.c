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

void stationary_dist(gsl_matrix *Q, gsl_vector *stationary, gslws_t *ws) {
	gsl_matrix *D = gsl_matrix_alloc(6,6);
	gsl_matrix_complex *D_bullshit = gsl_matrix_complex_alloc(6,6);
	gsl_matrix_complex *D_inv = gsl_matrix_complex_alloc(6,6);
	gsl_matrix *S = gsl_matrix_alloc(6,6);
	gsl_matrix *realevecs = gsl_matrix_alloc(6,6);
	gsl_matrix_complex *evecs = gsl_matrix_complex_alloc(6,6);
	gsl_matrix_complex *evecs_inv = gsl_matrix_complex_alloc(6,6);
	gsl_vector_complex *evals = gsl_vector_complex_alloc(6);
	gsl_vector *phi = gsl_vector_alloc(6);
	gsl_permutation *p = gsl_permutation_alloc(6);
	int signum;
	double norm;
	int i, j;
	gsl_complex bullshit;
	gsl_complex zero;

	GSL_SET_COMPLEX(&zero, 0.0, 0.0);

	gsl_matrix_complex_set_all(D_bullshit, zero);
	for(i=0; i<6; i++) {
		GSL_SET_COMPLEX(&bullshit, gsl_matrix_get(Q, i, i), 0.0);
		gsl_matrix_complex_set(D_bullshit, i, i, bullshit);
	}
	invert(D_bullshit, D_inv, p, &signum);
	for(i=0; i<6; i++) {
		for(j=0; j<6; j++) {
			gsl_matrix_set(D, i, j, GSL_REAL(gsl_matrix_complex_get(D_inv, i, j)));
		}
	}

	// S = I - DinvQ
	gsl_matrix_set_identity(S);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.0, D, Q, 1.0, S);
	printf("S:\n");
	fprint_matrix(stdout, S);
	printf("----------\n");
	gsl_matrix_transpose(S);
	decompose_q(S, ws);
	for(i=0; i<6; i++) {
		for(j=0; j<6; j++) {
			gsl_matrix_complex_set(evecs, i, j, gsl_complex_conjugate(gsl_matrix_complex_get(evecs, i, j)));
		}
	}
	for(i=0; i<6; i++) {
		printf("%f + %fi\n", GSL_REAL(gsl_vector_complex_get(evals, i)), GSL_IMAG(gsl_vector_complex_get(evals, i)));
	}
	realify_matrix(evecs, realevecs);
	fprint_matrix(stdout, realevecs);
	for(i=0; i<6; i++) {
		if(abs(GSL_IMAG(gsl_vector_complex_get(evals, i))) < 0.00001 && abs(GSL_REAL(gsl_vector_complex_get(evals, i)) - 1) < 0.0001) {
			printf("MatcheD!\n");
			for(j=0; j<6; j++) {
				gsl_vector_set(phi, j, GSL_REAL(gsl_matrix_complex_get(evecs, j, i)));
			}
			break;
		}
	}
	norm = 0;
	for(i=0; i<6; i++) {
		norm += gsl_vector_get(phi, i);
	}
	printf("Norm: %f\n", norm);
	gsl_blas_dscal(1.0/norm, phi);

	gsl_matrix_transpose(D);
	gsl_blas_dgemv(CblasNoTrans, 1.0, D, phi, 1.0, stationary);

	norm = 0;
	for(i=0; i<6; i++) {
		norm += gsl_vector_get(stationary, i);
	}
	printf("Norm: %f\n", norm);
	gsl_blas_dscal(1.0/norm, stationary);
	fprint_vector(stdout, stationary);

	gsl_matrix_free(D);
	gsl_matrix_complex_free(D_bullshit);
	gsl_matrix_complex_free(D_inv);
	gsl_matrix_free(S);
	gsl_matrix_free(realevecs);
	gsl_matrix_complex_free(evecs);
	gsl_matrix_complex_free(evecs_inv);
	gsl_vector_complex_free(evals);
	gsl_vector_free(phi);
}
