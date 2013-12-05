#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_matrix_complex_float.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_vector_complex.h>

#define MIN_PROTOAGE 35000
#define MAX_PROTOAGE 100000
#define STEPS 100.0

#include "tree.h"

struct ubertree {
	FILE *fp;
	int i, j, k, c;
	
	double protodist[6];
	double protodist_sum[100][6];
	float branch_lengths[6];
	unsigned long int sample_count;
};
typedef struct ubertree ubertree_t;

void initialise_ubertree(ubertree_t *ut);
void update_ubertree(ubertree_t *ut, node_t **trees, gsl_matrix *Q, gslws_t *ws);
float load_tree_age(char *filename);
void compute_protodist(double protoage, ubertree_t *ut, node_t **trees, gsl_matrix *Q, gslws_t *ws);
void save_ubertree(ubertree_t *ut, char *directory);
