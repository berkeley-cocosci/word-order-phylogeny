#ifndef BELIEFPROP_H
#define BELIEFPROP_H
 
#include "gslworkspace.h"
#include "tree.h"

void reset_tree(node_t *node);
void upwards_belprop(FILE *fp, node_t *tree, gsl_matrix *Q, gslws_t *ws);
long double get_tree_likelihood(FILE *fp, node_t *node, int value, gslws_t *ws);

#endif /* BELIEFPROP_H */
