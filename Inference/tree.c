#include <errno.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "tree.h"
#include "beliefprop.h"

leafdata_t* read_leaf_data(char* filename) {
	FILE *fp;
	int wordorder;
	int ishead = 1;
	leafdata_t *head, *node;
	char langname[8];
	//printf("Initial langname pointer is %p\n", langname);
	//printf("Initial wordorder pointer is %p\n", &wordorder);
	fp = fopen(filename, "r");
	//printf("File pointer is: %p\n", fp);
	int count = 42;
	while(count != EOF) {
		count = fscanf(fp, "%s\t%d\n", langname, &wordorder);
		//printf("Fscanf returned: %d\n", count);
		if(ishead) {
			node = malloc(sizeof(leafdata_t));
			head = node;
			ishead = 0;
		} else {
			node->next = malloc(sizeof(leafdata_t));
			node = node->next;
		}
		node->langname = calloc(8, sizeof(char));
		//printf("Node langname pointer is %p\n", node->langname);
		strncpy(node->langname, langname, 8);
		//printf("File says langname is: %s\n", langname);
		//printf("I've copied in: %s\n", node->langname);
		node->wordorder = wordorder - 1;
//		printf("Assigning word order of: %d\n", node->wordorder);
		node->next = NULL;
		//printf("Now langname pointer is %p\n", langname);
		//printf("Now wordorder pointer is %p\n", &wordorder);
		//printf("File pointer is: %p\n", fp);
		
	}
	fclose(fp);
	return head;
}

void shuffle_leafdata(leafdata_t *head, gsl_rng *r) {
	int i;
	leafdata_t *node;
	int leafcount;
	int *wordorders;
	/* Figure out length of list */
	node = head;
	leafcount = 0;
	while(1) {
		if(node->next == NULL) {
			break;
		} else {
			leafcount++;
			node = node->next;
		}
	}
	/* Allocate array for word orders */
	wordorders = calloc(leafcount, sizeof(int));
	/* Populate array */
	node = head;
	for(i=0; i<leafcount; i++) { 
		wordorders[i] = node->wordorder;
		node = node-> next;
	}
	/* Shuffle array */
	gsl_ran_shuffle(r, wordorders, leafcount, sizeof(int));
	/* Change leaf word orders */
	node = head;
	for(i=0; i<leafcount; i++) { 
		node->wordorder = wordorders[i];
		node = node-> next;
	}
}

uint8_t get_leaf_wordorder(leafdata_t* node, char* langname) {
	//printf("Target langname: %s\n", langname);
	while(1) {
		//printf("Current node langname: %s\n", node->langname);
		if(strcmp(node->langname, langname) == 0) {
			return node->wordorder;
		}
		if(node->next == NULL) {
			return 66;
		} else {
			node = node->next;
		}
	}	
}

int get_node_count(char *filename) {
	FILE *fp;
	int to, from, max, dir;
	float length;
	char string[64];
	max = 0;
	fp = fopen(filename, "r");
	if(fp == NULL) {
		printf("Couldn't open %s!\n", filename);
	}
	while(fscanf(fp, "%d, %d, %d, %f, %s\n", &from, &dir, &to, &length, string) != EOF) {
		//printf("%d, %d, %d, %f, %s\n", from, dir, to, length, string);
		if(to > max) max = to;
		//printf("Max: %d\n");
	}
	fclose(fp);
	return max;
}

void populate_nodes(node_t *nodes, leafdata_t* leafdata, char *filename) {
	FILE *fp;
	int i;
	int to, from, dir, count;
	float length;
	char string[64];
	fp = fopen(filename, "r");
	if(fp == NULL) {
		printf("Couldn't open %s!\n", filename);
	}
	//printf("fp is initially: %d\n", fp);
	//printf("EOF is: %d\n", EOF);
	//printf("EBADF is: %d\n", EBADF);
	while(1) {
		count = fscanf(fp, "%d, %d, %d, %f, %s\n", &from, &dir, &to, &length, string);
/*		if(length <= 0.0) {
			fprintf(stderr, "Tree file is invalid:\n");
			fprintf(stderr, "Tree file specifies a branch of length 0!\n");
			fprintf(stderr, "Dying now.\n");
			exit(3);
		}
*/
		if(length == 0.0) length = 0.000001;
		//printf("%d\n", count);
		//printf("%d\n", ferror(fp));
		if(count == EOF) {
			//printf("fp is now: %d\n", fp);
			if(ferror(fp)) {
				perror("fscanf");
			}
			break;
		}
		else if(count == 0) {
			fprintf(stderr, "Error: fscanf early matching failure\n");
		}
		else if(count != 5) {
			fprintf(stderr, "Error: fscanf matched less items than expected\n");
		}

		if(ferror(fp)) {
			perror("fscanf");
		}


		//printf("%d, %d, %d, %f, %s\n", from, dir, to, length, string);
		nodes[to].parent = &nodes[from];
		nodes[to].ready_to_feed = 0;
		nodes[from].ready_to_feed = 0;
		if(dir) {
			/* dir = 1 means right */
			nodes[from].right_branch = length;
			nodes[from].right_child = &nodes[to];
		} else {
			/* dir = 0 means left */
			nodes[from].left_branch = length;
			nodes[from].left_child = &nodes[to];
		}
		//printf("String is: %s\n", string);
		//printf("String cmp NON-LEAF is: %d\n", strcmp(string, "NON-LEAF"));
		if(!strcmp(string, "NON-LEAF") == 0) {
			strcpy(nodes[to].nodename, string);
			nodes[to].left_branch = 0.0;
			nodes[to].right_branch = 0.0;
			nodes[to].left_child = NULL;
			nodes[to].right_child = NULL;
			for(i=0; i<6; i++) {
				nodes[to].dist[i] = 0;
				nodes[to].l_message[i] = 0;
			}
			i = get_leaf_wordorder(leafdata, string);
			if(i==66) {
				fprintf(stderr, "Leaf word order lookup failure!");
				fprintf(stderr, "Couldn't find leaf named: %s\n", string);
				fprintf(stderr, "Dying now.\n");
				exit(2);
			}
			if(i != 42) {
				nodes[to].dist[i] = 1.0;
				nodes[to].l_message[i] = 1.0;
				nodes[to].ready_to_feed = 1;
			}
		}
	}
	fclose(fp);
}

void fprint_tree(FILE *fp, node_t *nodes, int length) {
	int i;
	for(i=0; i<length; i++) {
		fprintf(fp, "Left child: %p\n", nodes[i].left_child);
		fprintf(fp, "Right child: %p\n", nodes[i].right_child);
		fprintf(fp, "Ready to feed: %d\n", nodes[i].ready_to_feed);
		fprintf(fp, "Has passed: %d\n", nodes[i].has_passed);
		fprintf(fp, "----------\n");
	
	}
}

void die_due_to_neg() {
	fprintf(stderr, "Couldn't remove all negative branches!\n");
	fprintf(stderr, "Dying now.\n");
	exit(42);
}

void fix_negative_branches(node_t *node) {
	if(node->left_branch < 0) {
		fprintf(stderr, "Okay, I have a left branch of: %f\n", node->left_branch);
		if(node->left_child != NULL) {
			// Try to push the negative distance down onto the child
			fprintf(stderr, "I'm going to push down\n");
			fprintf(stderr, "My left child has a left branch of: %f\n", node->left_child->left_branch);
			fprintf(stderr, "My left child has a right branch of: %f\n", node->left_child->left_branch);
			fprintf(stderr, "Doing some tinkering...\n");
			node->left_child->left_branch += node->left_branch;
			node->left_child->right_branch += node->left_branch;
			fprintf(stderr, "Now my left child has a left branch of: %f\n", node->left_child->left_branch);
			fprintf(stderr, "Now my left child has a right branch of: %f\n", node->left_child->left_branch);
			node->left_branch = 0.0;
		} else {
			// Try to push the negative distance up onto the parent.
			fprintf(stderr, "I'm going to push up\n");
			if(node->parent->left_child == node) {
				fprintf(stderr, "I'm a left child\n");
				fprintf(stderr, "The branch to me from my parent has length: %f\n", node->parent->left_branch);
				if(node->parent->left_branch > -1*node->left_branch) {
					node->parent->left_branch += node->left_branch;
					node->left_branch = 0.0;
				} else {
					// Shit.  Set to zero out of desparation
					node->left_branch = 0;
				}
			} else {
				fprintf(stderr, "I'm a right child\n");
				fprintf(stderr, "The branch to me from my parent has length: %f\n", node->parent->right_branch);
				if(node->parent->right_branch > -1*node->left_branch) {
					node->parent->right_branch += node->left_branch;
					node->left_branch = 0.0;
				} else {
					die_due_to_neg();
				}
			
			}
		}
	}
	if(node->right_branch < 0) {
		if(node->right_child != NULL) {
			node->right_child->left_branch += node->right_branch;
			node->right_child->right_branch += node->right_branch;
			node->right_branch = 0.0;
		} else {
			// Try to push the negative distance up onto the parent.
			if(node->parent->left_child == node) {
				if(node->parent->left_branch > -1*node->right_branch) {
					node->parent->left_branch += node->right_branch;
					node->right_branch = 0.0;
				} else {
					die_due_to_neg();
				}
			} else {
				if(node->parent->right_branch > -1*node->right_branch) {
					node->parent->right_branch += node->right_branch;
					node->right_branch = 0.0;
				} else {
					die_due_to_neg();
				}
			
			}
		}
	}
	if(node->left_child != NULL) fix_negative_branches(node->left_child);
	if(node->right_child != NULL) fix_negative_branches(node->right_child);
}

void compute_root_length(node_t *node) {
	if(node->left_child != NULL) {
		node->left_child->distance_from_root = node->distance_from_root + node->left_branch;
		compute_root_length(node->left_child);
	}
	if(node->right_child != NULL) {
		node->right_child->distance_from_root = node->distance_from_root + node->right_branch;
		compute_root_length(node->right_child);
	}
}

void get_leaves(node_t *node, node_t ***listhead, int *nodecount) {
        if(node->left_child == NULL && node->right_child == NULL) {
                (*nodecount)++;
                *listhead = realloc(*listhead, (*nodecount) * sizeof(node_t*));
                (*listhead)[*nodecount-1] = node;
        }
        if(node->left_child != NULL) get_leaves(node->left_child, listhead, nodecount);
        if(node->right_child != NULL) get_leaves(node->right_child, listhead, nodecount);
}

void print_out_tree(node_t *node) {
	printf("I'm the node at %p.\n", node);
	printf("My left child is %p, at a distance of %f\n", node->left_child, node->left_branch);
	printf("My right child is %p, at a distance of %f\n", node->right_child, node->right_branch);
	if(node->left_child != NULL) print_out_tree(node->left_child);
	if(node->right_child != NULL) print_out_tree(node->right_child);
}

void double_branch_fixer(node_t *node) {
	int childcount = 0;
	if(node->left_child != NULL) childcount++;
	if(node->right_child != NULL) childcount++;
	if(childcount == 1) {
		if(node->parent->left_child == node) {
			if(node->left_child != NULL) {
				node->parent->left_child = node->left_child;
				node->parent->left_branch += node->left_branch;
				node->left_child->parent = node->parent;
			} else {
				node->parent->left_child = node->right_child;
				node->parent->left_branch += node->right_branch;
				node->right_child->parent = node->parent;
			}
		}
		if(node->parent->right_child == node) {
			if(node->left_child != NULL) {
				node->parent->right_child = node->left_child;
				node->parent->right_branch += node->left_branch;
				node->left_child->parent = node->parent;
			} else {
				node->parent->right_child = node->right_child;
				node->parent->right_branch += node->right_branch;
				node->right_child->parent = node->parent;
			}
		}
	}
	if(node->left_child != NULL) double_branch_fixer(node->left_child);
	if(node->right_child != NULL) double_branch_fixer(node->right_child);
}

void verify_tree_goodness(node_t *node) {
	int childcount = 0;
	if(node->left_child != NULL) childcount++;
	if(node->right_child != NULL) childcount++;
	if(!(childcount == 0 || childcount == 2)) {
		printf("Tree is not binary!  Counted %d children.  Aborting.\n", childcount);
		exit(69);
	}
        if(node->left_child != NULL) verify_tree_goodness(node->left_child);
        if(node->right_child != NULL) verify_tree_goodness(node->right_child);
}

void unknown_data_leafectomy(node_t *root) {
	node_t **listhead=NULL;
	int i, j;
	int knowndata;
	int leafcount = 0;
	int perfection = 0;

	/* First, delete links to nodes with unknown data */
	while(1) {
		perfection = 1;
		leafcount = 0;
		get_leaves(root, &listhead, &leafcount);
		for(i=0; i<leafcount; i++) {
			knowndata = 0;
			for(j=0; j<6; j++) {
				if(listhead[i]->l_message[j]) knowndata = 1;
			}
			if(!knowndata) {
				perfection = 0;
				if(listhead[i]->parent->left_child == listhead[i])  listhead[i]->parent->left_child = NULL;
				if(listhead[i]->parent->right_child == listhead[i]) listhead[i]->parent->right_child = NULL;
				listhead[i]->parent = NULL;
			}
		}
		if(perfection) break;
	}

	/* Now, replace "double branches" with single branches */
	double_branch_fixer(root);
}

node_t *build_tree(char *treefile, char *leaffile, int shuffle_leaves, gsl_rng *r) {
	int i, j, k;
	leafdata_t *leafdata;
	node_t *nodes;
	float maxlength;
	leafdata = read_leaf_data(leaffile);
	if(shuffle_leaves) shuffle_leafdata(leafdata, r);
	i = get_node_count(treefile) + 1;
	nodes = calloc(i, sizeof(node_t));
	if(nodes == NULL) {
		fprintf(stderr, "Couldn't allocate memory for tree!\n");
		fprintf(stderr, "Dying now.");
		exit(1);
	}
	for(j=0; j<i; j++) {
		nodes[j].ready_to_feed = 0;
		nodes[j].has_passed = 0;
		nodes[j].got_left = 0;
		nodes[j].got_right = 0;
		for(k=0; k<6; k++) {
			nodes[j].has_cached[k] = 0;
			nodes[j].cache[k] = 0.0;
		}
	}
	populate_nodes(nodes, leafdata, treefile);
	unknown_data_leafectomy(&nodes[0]);
	verify_tree_goodness(&nodes[0]);

	compute_root_length(&nodes[0]);
	maxlength = 0;
	for(j=0; j<i; j++) {
		if(nodes[j].parent != NULL)  {
			if(nodes[j].left_child == NULL && nodes[j].right_child == NULL) {
				if(nodes[j].distance_from_root > maxlength) maxlength = nodes[j].distance_from_root;
			}
		}
	}
	printf("Maximum distance from root is %f\n", maxlength);
	printf("About to return this: %p\n", &nodes[0]);
	return &nodes[0];
}

void load_tree(node_t **tree, char *dir, int method, int family, int treeindex, int shuffle) {
	unsigned long int seed;
	FILE *fp;
	char treefile[1024];
	char leaffile[1024];
	char families[][16] = {"indo", "austro", "niger", "afro", "nilo", "sino"};
	char types[][16] = {"geographic", "genetic", "feature", "combination" };
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus);
	/* Seed PRNG */
	fp = fopen("/dev/urandom", "r");
	fread(&seed, sizeof(seed), 1, fp);
	fclose(fp);
	gsl_rng_set(rng, seed);

	sprintf(treefile, "%s/%s/%s/tree_%d.simple", dir, types[method], families[family], treeindex+1);
	sprintf(leaffile, "../TreeBuilder/generated_trees/%s.leafdata", families[family]);
	(*tree) = build_tree(treefile, leaffile, shuffle, rng);
	reset_tree(*tree);
}
