#include <unistd.h>
#include<math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_matrix_complex_float.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_vector_complex.h>
#include<gsl/gsl_statistics_double.h>

#include "gslworkspace.h"
#include "tree.h"
#include "matrix.h"
#include "mcmc.h"
#include "modellike.h"
#include "beliefprop.h"
#include "sampman.h"
#include "ubertree.h"


void do_single_tree_inference(FILE *logfp, mcmc_t *mcmc, node_t *tree, gslws_t *ws, sampman_t *sm, int burnin, int samples, int lag) {
	int i, j;
	/* Burn in */
	for(i=0; i<burnin; i++) single_tree_mcmc_iteration(logfp, mcmc, tree, ws);
	/* Take samples */
	for(i=0; i<samples; i++) {
		if(gsl_rng_uniform_int(mcmc->r, 10000) >= 9999) {
			/* Random restart! */
			random_restart(mcmc);
			compute_single_tree_probabilities(mcmc, tree, ws);
			for(j=0; j<burnin; j++) single_tree_mcmc_iteration(logfp, mcmc, tree, ws);
		}

		for(j=0; j<lag; j++) {
			single_tree_mcmc_iteration(logfp, mcmc, tree, ws);
		}

		build_q(mcmc);
		upwards_belprop(logfp, tree, mcmc->Q, ws);

		/* Record sample */
		process_sample(sm, mcmc, ws, tree);
	}
}

void do_multi_tree_inference(FILE *logfp, mcmc_t *mcmc, node_t **subtrees, node_t **wholetrees, gslws_t *wses, sampman_t *sms, ubertree_t *ut, int burnin, int samples, int lag) {
	int i, j;
	/* Burn in */
	for(i=0; i<burnin; i++) multi_tree_mcmc_iteration(logfp, mcmc, subtrees, wses);
	/* Take samples */
	for(i=0; i<samples; i++) {
		if(gsl_rng_uniform_int(mcmc->r, 10000) >= 9999) {
			/* Random restart! */
			random_restart(mcmc);
			compute_multi_tree_probabilities(mcmc, subtrees, wses);
			for(j=0; j<burnin; j++) multi_tree_mcmc_iteration(logfp, mcmc, subtrees, wses);
		}

		for(j=0; j<lag; j++) {
			multi_tree_mcmc_iteration(logfp, mcmc, subtrees, wses);
		}

		build_q(mcmc);
		for(j=0; j<6; j++) {
			upwards_belprop(logfp, wholetrees[j], mcmc->Q, &wses[j]);
			process_sample(&sms[j], mcmc, &wses[j], wholetrees[j]);
		}
		if(ut != NULL) update_ubertree(ut, wholetrees, mcmc->Q, &wses[0]);
	}
}

void do_balanced_multi_tree_inference(FILE *logfp, mcmc_t *mcmc, node_t **subtrees, node_t **wholetrees, gslws_t *wses, sampman_t *sms, ubertree_t *ut, int burnin, int samples, int lag) {
	int i, j;
	/* Burn in */
	for(i=0; i<burnin; i++) balanced_multi_tree_mcmc_iteration(logfp, mcmc, subtrees, wses);
	/* Take samples */
	for(i=0; i<samples; i++) {
		if(gsl_rng_uniform_int(mcmc->r, 10000) >= 9999) {
			/* Random restart! */
			random_restart(mcmc);
			compute_balanced_multi_tree_probabilities(mcmc, subtrees, wses);
			for(j=0; j<burnin; j++) balanced_multi_tree_mcmc_iteration(logfp, mcmc, subtrees, wses);
		}

		for(j=0; j<lag; j++) {
			balanced_multi_tree_mcmc_iteration(logfp, mcmc, subtrees, wses);
		}

		build_q(mcmc);
		for(j=0; j<6; j++) {
			upwards_belprop(logfp, wholetrees[j], mcmc->Q, &wses[j]);
			process_sample(&sms[j], mcmc, &wses[j], wholetrees[j]);
		}
		if(ut != NULL) update_ubertree(ut, wholetrees, mcmc->Q, &wses[0]);
	}
}

void bespoke(char *treefile, char *leaffile, int burnin, int samples, int lag, char *outdir, int logging) {
	FILE *logfp;
	node_t *tree;
	mcmc_t mcmc;
	gslws_t ws;
	sampman_t sm;

	// Open log file
	if(logging) {
		logfp = fopen("logfile", "w");
	} else {
		logfp = fopen("/dev/null", "w");
	}

	initialise_mcmc(&mcmc);
	initialise_sampman(&sm, outdir);
	alloc_gslws(&ws);
	// Compute stability sampling rate
	sm.stability_sampling_rate = samples / 500.0;

	tree = build_tree(treefile, leaffile, 0, NULL);
	compute_single_tree_probabilities(&mcmc, tree, &ws);
	do_single_tree_inference(logfp, &mcmc, tree, &ws, &sm, burnin, samples, lag);
	free(tree);
	compute_means(&sm);
	save_indiv_q(outdir, &sm);
	fclose(logfp);
}

void whole_shared_q(int method, int shuffle, int burnin, int samples, int lag, int treecount, char *outdir, int logging) {
	FILE *logfp;
	int treeindex, i;
	char methods[][16] = {"geographic", "genetic", "feature", "combination" };
	char filename[1024];
	node_t **trees = calloc(6, sizeof(node_t*));
	node_t **subtrees = calloc(6, sizeof(node_t*));
	mcmc_t mcmc;
	gslws_t *wses = calloc(6, sizeof(gslws_t));
	sampman_t *sms = calloc(6, sizeof(sampman_t));
	ubertree_t ut;

	// Open log file
	if(logging) {
		logfp = fopen("logfile", "w");
	} else {
		logfp = fopen("/dev/null", "w");
	}

	initialise_mcmc(&mcmc);
	initialise_ubertree(&ut);

	for(i=0; i<6; i++) {
		alloc_gslws(&wses[i]);
		initialise_sampman(&sms[i], outdir);
		sms[i].stability_sampling_rate = (treecount*samples) / 500;
	}

	// Loop over trees...
	for(treeindex=0; treeindex<treecount; treeindex++) {
		/* Build tree(s) */
		for(i=0; i<6; i++) load_tree(&trees[i], "../TreeBuilder/generated_trees/whole/", method, i, treeindex, shuffle);
		for(i=0; i<6; i++) load_tree(&subtrees[i], "../TreeBuilder/generated_trees/whole/", method, i, treeindex, shuffle);
		/* Subsample languages to match global statistics */
		//subsample(&subtrees, mcmc.r);
		/* Draw samples for this tree (set) */
		compute_multi_tree_probabilities(&mcmc, subtrees, wses);
		do_multi_tree_inference(logfp, &mcmc, subtrees, trees, wses, sms, &ut, burnin, samples, lag);
		/* Free up tree memory */
		for(i=0; i<6; i++) free(trees[i]);
	}

	// Finish up
	for(i=0; i<6; i++) compute_means(&sms[i]);
	sprintf(filename, "%s/%s/", outdir, methods[method]);
	save_common_q(filename, sms);
	save_ubertree(&ut, filename);
	fclose(logfp);
}

void balanced_whole_shared_q(int method, int shuffle, int burnin, int samples, int lag, int treecount, char *outdir, int logging) {
	FILE *logfp;
	int treeindex, i;
	char methods[][16] = {"geographic", "genetic", "feature", "combination" };
	char filename[1024];
	node_t **trees = calloc(6, sizeof(node_t*));
	node_t **subtrees = calloc(6, sizeof(node_t*));
	mcmc_t mcmc;
	gslws_t *wses = calloc(6, sizeof(gslws_t));
	sampman_t *sms = calloc(6, sizeof(sampman_t));
	ubertree_t ut;

	// Open log file
	if(logging) {
		logfp = fopen("logfile", "w");
	} else {
		logfp = fopen("/dev/null", "w");
	}

	initialise_mcmc(&mcmc);
	initialise_ubertree(&ut);

	for(i=0; i<6; i++) {
		alloc_gslws(&wses[i]);
		initialise_sampman(&sms[i], outdir);
		sms[i].stability_sampling_rate = (treecount*samples) / 500;
	}

	// Loop over trees...
	for(treeindex=0; treeindex<treecount; treeindex++) {
		mcmc.dummy_tree_age = (55000+gsl_ran_gaussian(mcmc.r, 2500)) / 10000.0;
		/* Build tree(s) */
		for(i=0; i<6; i++) load_tree(&trees[i], "../TreeBuilder/generated_trees/whole/", method, i, treeindex, shuffle);
		for(i=0; i<6; i++) load_tree(&subtrees[i], "../TreeBuilder/generated_trees/whole/", method, i, treeindex, shuffle);
		/* Subsample languages to match global statistics */
		//subsample(&subtrees, mcmc.r);
		/* Draw samples for this tree (set) */
		compute_balanced_multi_tree_probabilities(&mcmc, subtrees, wses);
		do_balanced_multi_tree_inference(logfp, &mcmc, subtrees, trees, wses, sms, &ut, burnin, samples, lag);
		/* Free up tree memory */
		for(i=0; i<6; i++) free(trees[i]);
	}

	// Finish up
	for(i=0; i<6; i++) compute_means(&sms[i]);
	sprintf(filename, "%s/%s/", outdir, methods[method]);
	save_common_q(filename, sms);
	save_ubertree(&ut, filename);
	fclose(logfp);
}

void whole_indiv_q(int method, int shuffle, int burnin, int samples, int lag, int treecount, char *outdir, int logging) {
	FILE *logfp;
	uint8_t family;
	int treeindex;
	char filename[1024];
	char families[][16] = {"afro", "austro", "indo", "niger", "nilo", "sino"};
	char types[][16] = {"geographic", "genetic", "feature", "combination" };
	node_t *tree;
	mcmc_t mcmc;
	gslws_t ws;
	sampman_t sm;

	// Open log file
	if(logging) {
		logfp = fopen("logfile", "w");
	} else {
		logfp = fopen("/dev/null", "w");
	}

	initialise_mcmc(&mcmc);
	initialise_sampman(&sm, outdir);
	alloc_gslws(&ws);
	// Compute stability sampling rate
	sm.stability_sampling_rate = (treecount*samples) / 500;

	// Loop over families
	for(family=0; family<6; family++) {
		reset_sampman(&sm);
		// Loop over trees
		for(treeindex=0; treeindex<treecount; treeindex++) {
			/* Build tree */
			load_tree(&tree, "../TreeBuilder/generated_trees/whole/", method, family, treeindex, shuffle);
			/* Draw samples for this tree (set) */
			compute_single_tree_probabilities(&mcmc, tree, &ws);
			do_single_tree_inference(logfp, &mcmc, tree, &ws, &sm, burnin, samples, lag);
			/* Free up tree memory */
			free(tree);
		}

		// Finish up
		compute_means(&sm);
		sprintf(filename, "%s/%s/%s/", outdir, types[method], families[family]);
		save_indiv_q(filename, &sm);
	}

	fclose(logfp);
}

void split_shared_q(int method, int shuffle, int burnin, int samples, int lag, int treecount, char *outdir, int logging) {
	FILE *logfp;
	FILE *correlfp;
	int treeindex, i, j, k;
	char filename[1024];
	char methods[][16] = {"geographic", "genetic", "feature", "combination" };
	node_t **trees1 = calloc(6, sizeof(node_t*));
	node_t **subtrees1 = calloc(6, sizeof(node_t*));
	node_t **trees2 = calloc(6, sizeof(node_t*));
	node_t **subtrees2 = calloc(6, sizeof(node_t*));
	mcmc_t mcmc1, mcmc2;
	gslws_t *wses1 = calloc(6, sizeof(gslws_t));
	gslws_t *wses2 = calloc(6, sizeof(gslws_t));
	sampman_t *sms1 = calloc(6, sizeof(sampman_t));
	sampman_t *sms2 = calloc(6, sizeof(sampman_t));
	gsl_matrix *q1 = gsl_matrix_alloc(6, 6);
	gsl_matrix *q2 = gsl_matrix_alloc(6, 6);
	double qcorrel1[36], qcorrel2[36];
	double anccorrel1[36], anccorrel2[36];
	double qcorrelation = 0;
	double anccorrelation = 0;

	// Open log file
	if(logging) {
		logfp = fopen("logfile", "w");
	} else {
		logfp = fopen("/dev/null", "w");
	}

	initialise_mcmc(&mcmc1);
	initialise_mcmc(&mcmc2);
	for(i=0; i<6; i++) {
		alloc_gslws(&wses1[i]);
		alloc_gslws(&wses2[i]);
		initialise_sampman(&sms1[i], outdir);
		initialise_sampman(&sms2[i], outdir);
		// Compute stability sampling rate
		sms1[i].stability_sampling_rate = (treecount*samples) / 500.0;
		sms2[i].stability_sampling_rate = (treecount*samples) / 500.0;
	}

	// Loop over trees...
	for(treeindex=0; treeindex<treecount; treeindex++) {

		/* Build tree(s) */
		for(i=0; i<6; i++) {
			load_tree(&trees1[i], "../TreeBuilder/generated_trees/split/1/", method, i, treeindex, shuffle);
			load_tree(&subtrees1[i], "../TreeBuilder/generated_trees/split/1/", method, i, treeindex, shuffle);
			load_tree(&trees2[i], "../TreeBuilder/generated_trees/split/2/", method, i, treeindex, shuffle);
			load_tree(&subtrees2[i], "../TreeBuilder/generated_trees/split/2/", method, i, treeindex, shuffle);
		}
		/* Subsample languages to match global statistics */
		subsample(&subtrees1, mcmc1.r);
		subsample(&subtrees2, mcmc2.r);
		/* Draw samples for this tree (set) */
		compute_multi_tree_probabilities(&mcmc1, subtrees1, wses1);
		compute_multi_tree_probabilities(&mcmc2, subtrees2, wses2);

		/* Burn in */
		for(i=0; i<burnin; i++) {
			multi_tree_mcmc_iteration(logfp, &mcmc1, subtrees1, wses1);
			multi_tree_mcmc_iteration(logfp, &mcmc2, subtrees2, wses2);
		}

		gsl_matrix_set_zero(q1);
		gsl_matrix_set_zero(q2);
		memset(anccorrel1, 0, 6*sizeof(double));
		memset(anccorrel2, 0, 6*sizeof(double));
		/* Take samples */
		for(i=0; i<samples; i++) {
			if(gsl_rng_uniform_int(mcmc1.r, 10000) >= 9999) {
				/* Random restart! */
				random_restart(&mcmc1);
				random_restart(&mcmc2);
				compute_multi_tree_probabilities(&mcmc1, subtrees1, wses1);
				compute_multi_tree_probabilities(&mcmc2, subtrees2, wses2);
				for(j=0; j<burnin; j++) {
					multi_tree_mcmc_iteration(logfp, &mcmc1, subtrees1, wses1);
					multi_tree_mcmc_iteration(logfp, &mcmc2, subtrees2, wses2);
				}
			}

			for(j=0; j<lag; j++) {
				multi_tree_mcmc_iteration(logfp, &mcmc1, subtrees1, wses1);
				multi_tree_mcmc_iteration(logfp, &mcmc2, subtrees2, wses2);
			}

			build_q(&mcmc1);
			build_q(&mcmc2);

			for(j=0; j<6; j++) {
				upwards_belprop(logfp, trees1[j], mcmc1.Q, &wses1[j]);
				upwards_belprop(logfp, trees2[j], mcmc2.Q, &wses2[j]);
				process_sample(&sms1[j], &mcmc1, &wses1[j], trees1[j]);
				process_sample(&sms2[j], &mcmc2, &wses2[j], trees2[j]);

			}

			/* Combine Q samples for this tree pair */
			gsl_matrix_add(q1, mcmc1.Q);
			gsl_matrix_add(q2, mcmc2.Q);
			for(j=0; j<6; j++) {
				for(k=0; k<6; k++) {
					anccorrel1[j*6+k] += trees1[j]->dist[k];
					anccorrel2[j*6+k] += trees2[j]->dist[k];
				}
			}


		}
		
		gsl_matrix_scale(q1, 1.0 / samples);
		gsl_matrix_scale(q2, 1.0 / samples);
		for(i=0; i<6; i++) {
			for(j=0; j<6; j++) {
				qcorrel1[i*6+j] = gsl_matrix_get(q1, i, j);
				qcorrel2[i*6+j] = gsl_matrix_get(q2, i, j);
				anccorrel1[i*6+j] /= (1.0*samples);
				anccorrel2[i*6+j] /= (1.0*samples);
			}
		}
		qcorrelation += gsl_stats_correlation(qcorrel1, 1, qcorrel2, 1, 36);
		anccorrelation += gsl_stats_correlation(anccorrel1, 1, anccorrel2, 1, 36);

		/***************************************************/
		/* Free up tree memory */
		for(i=0; i<6; i++) {
			free(trees1[i]);
			free(subtrees1[i]);
			free(trees2[i]);
			free(subtrees2[i]);
		}
	}
	qcorrelation /= treecount;
	anccorrelation /= treecount;

	// Finish up
	for(i=0; i<6; i++) {
		compute_means(&sms1[i]);
		compute_means(&sms2[i]);
	}
	sprintf(filename, "%s/%s/", outdir, methods[method]);
	strcat(filename, "/left-half/");
	save_common_q(filename, sms1);
	sprintf(filename, "%s/%s/", outdir, methods[method]);
	strcat(filename, "/right-half/");
	save_common_q(filename, sms2);

	sprintf(filename, "%s/%s/correlations", outdir, methods[method]);
	correlfp = fopen(filename, "w");
	fprintf(correlfp, "Q: %f\n", qcorrelation);
	fprintf(correlfp, "Anc: %f\n", anccorrelation);
	fclose(correlfp);

	fclose(logfp);
}

void usage() {
	printf("Options:\n\n");
	printf("	-b BURNIN		Discard BURNIN samples at start of MCMC.\n");
	printf("	-c METHOD		Tree method: 0 is geo, 1 is gen, 2 is feat, 3 is combo.\n");
	printf("	-l LAG			Perform LAG MCMC iterations between samples.\n");
	printf("	-m 			Do common Q analysis (default is individual Q.\n");
	printf("	-s SAMPLES		Take SAMPLES MCMC samples per tree.\n");
	printf("	-S			Shuffle leaf word orders.\n");
	printf("	-x			Split trees.\n");
	printf("	-L			Be ridiculously verbose.\n");
}

int main(int argc, char **argv) {

	// Variable setup, allocation, etc.
	uint8_t ambiguous = 1;
	int c;
	int logging = 0;
	int balanced, multitree, treeclass, shuffle, split;
	char *outbase, *treefile, *leaffile;
	char outdir[1024];

	int burnin, lag, samples;
	int treecount;

	// Option parsing
	// defaults
	balanced = 0;
	shuffle = 0;
	treeclass = 0;
	multitree = 0;
	split = 0;
	burnin = 5000;
	lag = 100;
	samples = 1000;
	treecount = 100;
	outbase = "./results";
	treefile = NULL;
	leaffile = NULL;
	while((c = getopt(argc, argv, "b:Bc:i:l:ms:St:T:o:L:xv")) != -1) {
		switch(c) {
			case 'b':
				burnin = atoi(optarg);
				break;
			case 'B':
				balanced = 1;
				break;
			case 'c':
				treeclass = atoi(optarg);
				ambiguous = 0;
				break;
/*			case 'i':
				treeindex = atoi(optarg);
				break;
				*/
			case 'l':
				lag = atoi(optarg);
				break;
			case 'L':
				leaffile = optarg;
				ambiguous = 0;
				break;
			case 'm':
				multitree = 1;
				break;
			case 's':
				samples = atoi(optarg);
				break;
			case 'S':
				shuffle = 1;
				break;
			case 'x':
				split = 1;
				break;
			case 't':
				treecount = atoi(optarg);
				break;
			case 'T':
				treefile = optarg;
				ambiguous = 0;
				break;
			case 'o':
				outbase = optarg;
				break;
			case 'v':
				logging = 1;
				break;
		}
	}

	if(ambiguous) {
		usage();
		exit(0);
	}

	if(treefile && leaffile) {

		printf("Performing bespoke inference!");
		bespoke(treefile, leaffile, burnin, samples, lag, outbase, logging);

	} else if(balanced) {

		printf("Performing inference with one Q shared among all families, with counter-balancing languages.\n");
		sprintf(outdir, outbase);
		strcat(outdir, "/balanced-common-q/");
		balanced_whole_shared_q(treeclass, shuffle, burnin, samples, lag, treecount, outdir, logging);

	} else if(multitree && !split) {

		printf("Performing inference with one Q shared among all families.\n");
		sprintf(outdir, outbase);
		strcat(outdir, "/common-q/unsplit/");
		if(shuffle) {
			strcat(outdir, "/shuffled/");
		} else {
			strcat(outdir, "/unshuffled/");
		}
		whole_shared_q(treeclass, shuffle, burnin, samples, lag, treecount, outdir, logging);

	} else if(multitree && split) {

		printf("Performing inference with one Q shared among all families, with trees split in half.\n");
		sprintf(outdir, outbase);
		strcat(outdir, "/common-q/split/");
		if(shuffle) {
			strcat(outdir, "/shuffled/");
		} else {
			strcat(outdir, "/unshuffled/");
		}
		split_shared_q(treeclass, shuffle, burnin, samples, lag, treecount, outdir, logging);

	} else if(!multitree && !split) {
		printf("Performing inference with individual Qs per family.\n");
		sprintf(outdir, outbase);
		strcat(outdir, "/individual-q/unsplit/");
		if(shuffle) {
			strcat(outdir, "/shuffled/");
		} else {
			strcat(outdir, "/unshuffled/");
		}
		whole_indiv_q(treeclass, shuffle, burnin, samples, lag, treecount, outdir, logging);
	} else if(!multitree && split) {
		// Not implemented yet
		printf("Split individual Q not implemented yet.\n");
	} else {
		usage();
	}

	return 0;
}
