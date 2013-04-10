void save_results(char *outfile, gsl_matrix *Q, node_t **trees, gsl_vector **ancestral_sum, gsl_vector **ancestral_max, gsl_vector *stabs_sum, gsl_vector *stabs_max, gsl_matrix *trans_sum, gsl_matrix *trans_max, double max_posterior, int multitree) {

	int i;
	FILE *fp;
	gsl_vector_complex *evals = gsl_vector_complex_alloc(6);
	gsl_matrix_complex *evecs = gsl_matrix_complex_alloc(6, 6);
	gsl_matrix_complex *evecs_inv = gsl_matrix_complex_alloc(6, 6);
	gsl_matrix *P = gsl_matrix_alloc(6, 6);

	decompose_q(Q, evals, evecs, evecs_inv);
	compute_p(evals, evecs, evecs_inv, 1000000.0, P);

	fp = fopen(outfile, "w");
	if(multitree) {
		fprintf(fp, "Posterior ancestral word order distributions:\n");
		for(i=0; i<6; i++) fprint_vector(fp, ancestral_sum[i]);
		fprintf(fp, "----------\n");
	} else {
		fprintf(fp, "Posterior ancestral word order distribution:\n");
		fprint_vector(fp, ancestral_sum[0]);
		fprintf(fp, "----------\n");
	}
	fprintf(fp, "Posterior mean Q:\n");
	fprint_matrix(fp, Q);
	fprintf(fp, "----------\n");
	fprintf(fp, "Posterior mean stabilities:\n");
	fprint_vector(fp, stabs_sum);
	fprintf(fp, "----------\n");
	fprintf(fp, "Posterior mean transitions:\n");
	fprint_matrix(fp, trans_sum);
	fprintf(fp, "----------\n");
	compute_p(evals, evecs, evecs_inv, 0.1, P);
	fprintf(fp, "P matrix over short branch:\n");
	fprint_matrix(fp, P);
	fprintf(fp, "----------\n");
	compute_p(evals, evecs, evecs_inv, 1000000.0, P);
	fprintf(fp, "Posterior mean stationary:\n");
	fprint_matrix(fp, P);
	fprintf(fp, "----------\n");
	fprintf(fp, "Maximum posteror: %f\n", max_posterior);
	fprintf(fp, "----------\n");
	if(multitree) {
		fprintf(fp, "MAP ancestral word order distributions:\n");
		for(i=0; i<6; i++) fprint_vector(fp, ancestral_max[i]);
		fprintf(fp, "----------\n");
	} else {
		fprintf(fp, "MAP ancestral word order distribution:\n");
		fprint_vector(fp, ancestral_max[0]);
		fprintf(fp, "----------\n");
	}
	fprintf(fp, "MAP stabilities:\n");
	fprint_vector(fp, stabs_max);
	fprintf(fp, "MAP transitions:\n");
	fprint_matrix(fp, trans_max);
	fclose(fp);
}
