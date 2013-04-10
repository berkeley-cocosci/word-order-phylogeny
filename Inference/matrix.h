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

void vector_copy(gsl_vector *in, gsl_vector *out) {
        int i;
        for(i=0; i<6; i++) {
                gsl_vector_set(out, i, gsl_vector_get(in, i));
        }
}

int matrix_compare(gsl_matrix *a, gsl_matrix *b) {
        int i, j;
        for(i=0; i<6; i++) {
                for(j=0; j<6; j++) {
			if(gsl_matrix_get(a, i, j) != gsl_matrix_get(b, i, j)) return 0;
		}
	}
	return 1;
}

void matrix_copy(gsl_matrix *in, gsl_matrix *out) {
        int i, j;
        for(i=0; i<6; i++) {
                for(j=0; j<6; j++) {
                        gsl_matrix_set(out, i, j, gsl_matrix_get(in, i, j));
                }
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
	int i, j;
	for(i=0; i<6; i++) {
		for(j=0; j<6; j++) {
			if(GSL_IMAG(gsl_matrix_complex_get(complex, i, j)) != 0) {
//				printf("Shit!  Discarding imaginary stuff!");
			}
			gsl_matrix_set(real, i, j, GSL_REAL(gsl_matrix_complex_get(complex, i, j)));
		}
	}
}

void invert(gsl_matrix_complex *in, gsl_matrix_complex *out, gsl_permutation *p, int *signum) {
	/* Why in the fuck is this not just part of gsl? */
	gsl_linalg_complex_LU_decomp(in, p, signum);
//	print_matrix(in);
	gsl_linalg_complex_LU_invert(in, p, out);

}

void eigensolve(gsl_matrix *Q, gsl_vector_complex *evals, gsl_matrix_complex *evecs) {
	gsl_eigen_nonsymmv_workspace *ws = gsl_eigen_nonsymmv_alloc(6);
	gsl_eigen_nonsymmv(Q, evals, evecs, ws);
	gsl_eigen_nonsymmv_free(ws);
}

void compute_p(gsl_vector_complex *evals, gsl_matrix_complex *evecs, gsl_matrix_complex *evecs_inv, double t, gsl_matrix *out) {
	int i;
	gsl_matrix_complex *c_P = gsl_matrix_complex_alloc(6,6);
	gsl_matrix_complex *temp1 = gsl_matrix_complex_alloc(6,6);
	gsl_matrix_complex *temp2 = gsl_matrix_complex_alloc(6,6);
	gsl_matrix_complex_set_zero(temp1);
	gsl_matrix_complex_set_zero(temp2);
	gsl_complex one, zero, z;
	GSL_SET_COMPLEX(&one, 1.0, 0.0);
	GSL_SET_COMPLEX(&zero, 0.0, 0.0);
	/* Fill temp1 with e^(t\lbamda_i) */
	for(i=0; i<6; i++) {
		z = gsl_complex_mul_real(gsl_vector_complex_get(evals, i), t);
		gsl_matrix_complex_set(temp1, i, i, gsl_complex_exp(z));
	}
	gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, one, evecs, temp1, zero, temp2);
	/* Store temp2*evecs_inv in out */
	gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, one, temp2, evecs_inv, zero, c_P);
	realify_matrix(c_P, out);
	//print_matrix(out);
//	printf("---\n");

	gsl_matrix_complex_free(c_P);
	gsl_matrix_complex_free(temp1);
	gsl_matrix_complex_free(temp2);
}

void decompose_q(gsl_matrix *Q, gsl_vector_complex *evals, gsl_matrix_complex *evecs, gsl_matrix_complex *evecs_inv) {
	FILE *fp;
	gsl_matrix *Q_copy = gsl_matrix_alloc(6,6);
	gsl_matrix_complex *evecs_copy = gsl_matrix_complex_alloc(6,6);
	gsl_permutation *p = gsl_permutation_alloc(6);
	int signum;
	int i,j;

	/* Make a copy of Q to give to the (destructive) eigensolver */
	matrix_copy(Q, Q_copy);

	eigensolve(Q_copy, evals, evecs);
	//realify_vector(c_evals, evals);
	//realify_matrix(c_evecs, evecs);

	/* Make a copy of the eigenvectors to give to the (destructive) inverter */
	for(i=0; i<6; i++) {
		for(j=0; j<6; j++) {
			gsl_matrix_complex_set(evecs_copy, i, j, gsl_matrix_complex_get(evecs, i, j));
		}
	}
	
	/* Compute inverse eigenvector matrix */
	invert(evecs_copy, evecs_inv, p, &signum);

	gsl_matrix_free(Q_copy);
	gsl_matrix_complex_free(evecs_copy);
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

void stationary_dist(gsl_matrix *Q, gsl_vector *stationary) {
	gsl_matrix *D = gsl_matrix_alloc(6,6);
	gsl_matrix_complex *D_bullshit = gsl_matrix_complex_alloc(6,6);
	gsl_matrix_complex *D_inv = gsl_matrix_complex_alloc(6,6);
	gsl_matrix *S = gsl_matrix_alloc(6,6);
	gsl_matrix *realevecs = gsl_matrix_alloc(6,6);
	gsl_matrix_complex *evecs = gsl_matrix_complex_alloc(6,6);
	gsl_matrix_complex *evecs_inv = gsl_matrix_complex_alloc(6,6);
	gsl_vector_complex *evals = gsl_vector_complex_alloc(6);
	gsl_vector *phi = gsl_vector_alloc(6);
	gsl_vector *zerovec = gsl_vector_alloc(6);
	gsl_permutation *p = gsl_permutation_alloc(6);
	gsl_permutation *p2 = gsl_permutation_alloc(6);
	int signum;
	int signum2;
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
	decompose_q(S, evals, evecs, evecs_inv);
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
}
