#include<stdio.h>
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

void invert(gsl_matrix_complex in, gsl_matrix_complex *out, gsl_permutation *p, int *signum) {
	/* Why in the fuck is this not just part of gsl? */
	gsl_linalg_complex_LU_decomp(&in, p, signum);
	gsl_linalg_complex_LU_invert(&in, p, out);

}

void eigensolve(gsl_matrix *Q, gsl_vector_complex *evals, gsl_matrix_complex *evecs) {
	gsl_eigen_nonsymmv_workspace *ws = gsl_eigen_nonsymmv_alloc(6);
	gsl_eigen_nonsymmv(Q, evals, evecs, ws);
}

void exponential(gsl_matrix *in, double t, gsl_matrix_complex *out, gsl_vector_complex *evals, gsl_matrix_complex *evecs, gsl_matrix_complex *evecs_inv) {
	int i;
	gsl_complex one;
	gsl_complex zero;
	gsl_matrix_complex *temp1 = gsl_matrix_complex_alloc(6,6);
	gsl_matrix_complex *temp2 = gsl_matrix_complex_alloc(6,6);
	gsl_matrix_complex_set_zero(temp1);
	gsl_matrix_complex_set_zero(temp2);
	/* Fill temp1 with e^(t\lbamda_i) */
	for(i=0; i<6; i++) {
		one = gsl_complex_mul_real(gsl_vector_complex_get(evals, i), t);
		one = gsl_complex_exp(one);
		gsl_matrix_complex_set(temp1, i, i, one);
	}
	save_complex_matrix(temp1, "diageigenv.mat", 6, 6);
	/* Store evecs*temp1 in temp2 */
	GSL_SET_COMPLEX(&one, 1.0, 0.0);
	GSL_SET_COMPLEX(&zero, 0.0, 0.0);
	gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, one, evecs, temp1, zero, temp2);
	save_complex_matrix(temp2, "firstproduct.mat", 6, 6);
	/* Store temp2*evecs_inv in out */
	gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, one, temp2, evecs_inv, zero, out);
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

void save_complex_matrix(gsl_matrix_complex *M, char *filename, int x, int y) {
	int i, j;
	FILE *fp;
	fp = fopen(filename, "w");
	for(i=0; i<x; i++) {
		for(j=0; j<y-1; j++) {
			fprintf(fp, "%f + %fi	", GSL_REAL(gsl_matrix_complex_get(M, i,j)), GSL_IMAG(gsl_matrix_complex_get(M, i, j)));
		}
		fprintf(fp, "%f + %fi\n", GSL_REAL(gsl_matrix_complex_get(M, i,y-1)), GSL_IMAG(gsl_matrix_complex_get(M, i, y-1)));
	}
	fclose(fp);
}

void gsl_save_complex(gsl_matrix_complex *M, char *filename) {
	FILE *fp;
	fp = fopen(filename, "w");
	gsl_matrix_complex_fprintf(fp, M, "%f");
	fclose(fp);
}

void save_complex_vector(gsl_vector_complex *V, char *filename) {
	FILE *fp;
	fp = fopen(filename, "w");
	gsl_vector_complex_fprintf(fp, V, "%f");
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
	
	FILE *fp;
	int i, j, k;
	double a[] = {1,5,3,0,7,2,5,3,1};
	double b[] = {4,2,8,3,7,5,0,2,4};
	double c[] = {0,0,0,0,0,0,0,0,0};

	gsl_matrix_view A = gsl_matrix_view_array(a, 3, 3);
	gsl_matrix_view B = gsl_matrix_view_array(b, 3, 3);
	gsl_matrix_view C = gsl_matrix_view_array(c, 3, 3);

	gsl_matrix_complex *Ac = gsl_matrix_complex_alloc(3,3);
	gsl_matrix_complex *Bc = gsl_matrix_complex_alloc(3,3);
	gsl_matrix_complex *Cc = gsl_matrix_complex_alloc(3,3);

	gsl_complex z, z0, z1;
	gsl_matrix_complex_set_zero(Ac);
	gsl_matrix_complex_set_zero(Bc);
	gsl_matrix_complex_set_zero(Cc);

	GSL_SET_IMAG(&z, 0.0);
	GSL_SET_COMPLEX(&z0, 0.0, 0.0);
	GSL_SET_COMPLEX(&z1, 1.0, 0.0);
	for(i=0; i<3; i++) {
		for(j=0; j<3; j++) {
			GSL_SET_REAL(&z, gsl_matrix_get(&A.matrix, i, j));
			gsl_matrix_complex_set(Ac, i, j, z);
			GSL_SET_REAL(&z, gsl_matrix_get(&B.matrix, i, j));
			gsl_matrix_complex_set(Bc, i, j, z);
			GSL_SET_REAL(&z, gsl_matrix_get(&C.matrix, i, j));
			gsl_matrix_complex_set(Cc, i, j, z);
		}
	}
	
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &A.matrix, &B.matrix, 0.0, &C.matrix);
	
	gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, z1, Ac, Bc, z0, Cc);

	fp = fopen("testproduct", "w");
	gsl_matrix_fprintf(fp, &C.matrix, "%f");
	fclose(fp);	

	save_complex_matrix(Cc, "testproduccomp.mat", 3, 3);
}
