#ifndef BELIEFPROP_H
#define BELIEFPROP_H
 
#include "tree.h"

void reset_tree(node_t *node);
void upwards_belprop(FILE *fp, node_t **trees, gsl_matrix *Q, int multitree);

#endif /* BELIEFPROP_H */
