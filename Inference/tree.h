#ifndef TREE_H
#define TREE_H
 
#include <errno.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_rng.h>
<<<<<<< HEAD
#include <gsl/gsl_matrix.h>
=======
#include<gsl/gsl_matrix.h>
>>>>>>> dfb46fa8349a0836e1f1f74e28984358b153f49a

struct leafdata{
	char *langname;
	uint8_t wordorder;
	struct leafdata* next;
};
typedef struct leafdata leafdata_t;

struct node{
	char nodename[16];
	float left_branch;
	float right_branch;
	gsl_matrix *left_P;
	gsl_matrix *right_P;
	float distance_from_root;
	double age;
	struct node *parent;
	struct node *left_child;
	struct node *right_child;
	double p_message[6];
	double l_message[6];
	double r_message[6];
	double dist[6];
	int ml_order;
	int got_left;
	int got_right;
	int ready_to_feed;
	int has_passed;
	int has_cached[6];
	int has_cached_matrices;
	double cache[6];
};
typedef struct node node_t;

leafdata_t* read_leaf_data(char* filename);
uint8_t get_leaf_wordorder(leafdata_t* node, char* langname);
int get_node_count(char *filename);
void populate_nodes(node_t *nodes, leafdata_t* leafdata, char *filename);
void load_tree(node_t **tree, char *dir, int method, int family, int treeindex, int shuffle);
void fprint_tree(FILE *fp, node_t *nodes, int length);
void die_due_to_neg();
void fix_negative_branches(node_t *node);
void compute_root_length(node_t *node);
void get_leaves(node_t *node, node_t ***listhead, int *nodecount);
void print_out_tree(node_t *node);
void double_branch_fixer(node_t *node);
void verify_tree_goodness(node_t *node);
void unknown_data_leafectomy(node_t *root);
node_t *build_tree(char *treefile, char *leaffile, int shuffle_leaves, gsl_rng *r);

#endif /* TREE_H */
