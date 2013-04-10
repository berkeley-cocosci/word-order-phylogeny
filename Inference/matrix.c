#include<stdio.h>
#include<math.h>
#include<gsl/gsl_complex.h>
#include<gsl/gsl_complex_math.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_blas_types.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_matrix_complex_float.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_vector_complex.h>
#include<gsl/gsl_eigen.h>
#include<gsl/gsl_permutation.h>

//	gsl_vector_complex *evals = gsl_vector_complex_alloc(6);
//	gsl_matrix_complex *evecs = gsl_matrix_complex_alloc(6,6);

void save_matrix(gsl_matrix *M, char *filename);

void realify_vector(gsl_vector_complex *complex, gsl_vector *real) {
	int i;
	for(i=0; i<6; i++) {
		gsl_vector_set(real, i, GSL_REAL(gsl_vector_complex_get(complex, i)));
	}
}

void realify_matrix(gsl_matrix_complex *complex, gsl_matrix *real) {
	int i, j;
	for(i=0; i<6; i++) {
		for(j=0; j<6; j++) {
			gsl_matrix_set(real, i, j, GSL_REAL(gsl_matrix_complex_get(complex, i, j)));
		}
	}
}

void invert(gsl_matrix *in, gsl_matrix *out, gsl_permutation *p, int *signum) {
	/* Why in the fuck is this not just part of gsl? */
	gsl_linalg_LU_decomp(in, p, signum);
	gsl_linalg_LU_invert(in, p, out);

}

void eigensolve(gsl_matrix *Q, gsl_vector_complex *evals, gsl_matrix_complex *evecs) {
	gsl_eigen_nonsymmv_workspace *ws = gsl_eigen_nonsymmv_alloc(6);
	gsl_eigen_nonsymmv(Q, evals, evecs, ws);
}

void exponential(gsl_vector *evals, gsl_matrix *evecs, gsl_matrix *evecs_inv, double t, gsl_matrix *out) {
	int i;
	gsl_matrix *temp1 = gsl_matrix_alloc(6,6);
	gsl_matrix *temp2 = gsl_matrix_alloc(6,6);
	gsl_matrix_set_zero(temp1);
	gsl_matrix_set_zero(temp2);
	/* Fill temp1 with e^(t\lbamda_i) */
	for(i=0; i<6; i++) {
		gsl_matrix_set(temp1, i, i, exp(t * gsl_vector_get(evals, i)));
	}
	save_matrix(temp1, "diageigenv.mat");
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, evecs, temp1, 0.0, temp2);
	save_matrix(temp2, "firstproduct.mat");
	/* Store temp2*evecs_inv in out */
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, temp2, evecs_inv, 0.0, out);
}

void save_matrix(gsl_matrix *M, char *filename) {
	FILE *fp;
	int i, j;
	fp = fopen(filename, "w");
	for(i=0; i<6; i++) {
		for(j=0; j<5; j++) {
			fprintf(fp, "%f	", gsl_matrix_get(M, i,j));
		}
		fprintf(fp, "%f\n", gsl_matrix_get(M, i,5));
	}
	fclose(fp);
}

void save_complex_matrix(gsl_matrix_complex *M, char *filename) {
	int i, j;
	FILE *fp;
	fp = fopen(filename, "w");
	for(i=0; i<6; i++) {
		for(j=0; j<5; j++) {
			fprintf(fp, "%f + %fi	", GSL_REAL(gsl_matrix_complex_get(M, i,j)), GSL_IMAG(gsl_matrix_complex_get(M, i, j)));
		}
		fprintf(fp, "%f + %fi\n", GSL_REAL(gsl_matrix_complex_get(M, i,5)), GSL_IMAG(gsl_matrix_complex_get(M, i, 5)));
	}
	fclose(fp);
}

void gsl_save_complex(gsl_matrix_complex *M, char *filename) {
	FILE *fp;
	fp = fopen(filename, "w");
	gsl_matrix_complex_fprintf(fp, M, "%f");
	fclose(fp);
}

void save_vector(gsl_vector *V, char *filename) {
	FILE *fp;
	fp = fopen(filename, "w");
	gsl_vector_fprintf(fp, V, "%f");
	fclose(fp);
	
}

void initialise_q_matrix(gsl_matrix *Q) {
	int i;
	double q = 0.1;
	gsl_matrix_set_all(Q, q);
	for(i=0;i<6;i++) {
		gsl_matrix_set(Q,i,i,-5.0*q);
	}
}

void main() {
	gsl_matrix *Q = gsl_matrix_alloc(6,6);
	gsl_matrix *P = gsl_matrix_alloc(6,6);
	gsl_vector_complex *c_evals = gsl_vector_complex_alloc(6);
	gsl_matrix_complex *c_evecs = gsl_matrix_complex_alloc(6,6);
	gsl_vector *evals = gsl_vector_alloc(6);
	gsl_matrix *evecs = gsl_matrix_alloc(6,6);
	gsl_matrix *evecs_copy = gsl_matrix_alloc(6,6);
	gsl_matrix *evecs_inv = gsl_matrix_alloc(6,6);
	gsl_permutation *p = gsl_permutation_alloc(6);
	int signum;
	int i,j;

	initialise_q_matrix(Q);
	save_matrix(Q, "initialq.mat");

	/* Compute eigen-values and vectors and discard imaginary zeros */
	eigensolve(Q, c_evals, c_evecs);
	realify_vector(c_evals, evals);
	realify_matrix(c_evecs, evecs);

	/* Make a copy of the eigenvectors to give to the (destructive) inverter */
	for(i=0; i<6; i++) {
		for(j=0; j<6; j++) {
			gsl_matrix_set(evecs_copy, i, j, gsl_matrix_get(evecs, i, j));
		}
	}
	save_vector(evals, "eigenvals.vec");
	save_matrix(evecs, "eigenvecs_before.mat");
	
	/* Compute inverse eigenvector matrix */
	invert(evecs_copy, evecs_inv, p, &signum);
	save_matrix(evecs, "eigenvecs_after.mat");
	save_matrix(evecs_inv, "eigenvecs_inv.mat");

	/* Perform matrix exponentiation */
	exponential(evals, evecs, evecs_inv, 4.327, P);
	save_matrix(P, "p.mat");
}

