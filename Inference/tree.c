#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

struct node{
	float left_branch;
	float right_branch;
	struct node *parent;
	struct node *left_child;
	struct node *right_child;
};
typedef struct node node_t;

int get_node_count(char *filename) {
	FILE *fp;
	int to, from, max, dir;
	float length;
	char string[64];
	max = 0;
	fp = fopen(filename, "r");
	if(fp == NULL) {
		printf("Couldn't open that!\n");
	}
	while(fscanf(fp, "%d, %d, %d, %f, %s\n", &from, &dir, &to, &length, &string) != EOF) {
		printf("%d, %d, %d, %f, %s\n", from, dir, to, length, string);
		if(to > max) max = to;
	}
	fclose(fp);
	return max;
}

void populate_nodes(node_t *nodes, char *filename) {
	FILE *fp;
	int to, from, max, dir, count;
	float length;
	char string[64];
	fp = fopen(filename, "r");
	if(fp == NULL) {
		printf("Couldn't open that!\n");
	}
	printf("fp is initially: %d\n", fp);
	printf("EOF is: %d\n", EOF);
	printf("EBADF is: %d\n", EBADF);
	while(1) {
		count = fscanf(fp, "%d, %d, %d, %f, %s\n", &from, &dir, &to, &length, &string);
		printf("%d\n", count);
		printf("%d\n", ferror(fp));
		if(count == EOF) {
			printf("fp is now: %d\n", fp);
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
		printf("%d, %d, %d, %f, %s\n", from, dir, to, length, string);
		nodes[to].parent = &nodes[from];
		if(dir) {
			/* dir = 1 means right */
			nodes[from].right_branch = length;
			nodes[from].right_child = &nodes[to];
		} else {
			/* dir = 0 means left */
			nodes[from].left_branch = length;
			nodes[from].left_child = &nodes[to];
		}
		if(!strcmp(string, "NON-LEAF")) {
			nodes[to].left_branch = 0.0;
			nodes[to].right_branch = 0.0;
			nodes[to].left_child = NULL;
			nodes[to].right_child = NULL;
		}
	}
	fclose(fp);
}

node_t *build_tree(char *filename) {
	int i, j, k;
	node_t *nodes;
	printf("Countin'!\n");
	i = get_node_count(filename);
	printf("Found a total of %d nodes!\n", i);
	nodes = calloc(i, sizeof(node_t));
	printf("Callocin'!\n");
	if(nodes == NULL) {
		printf("Couldn't allocate memory for tree!\n");
	}
	printf("Populatin'!\n");
	populate_nodes(nodes, filename);
}

int main() {
	build_tree("easytree");
	fflush(stdout);
}

