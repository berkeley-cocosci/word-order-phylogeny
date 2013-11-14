TreeBuilder
====================

This is a messy pile of scripts geared toward solving the problem of randomly
generating a large number of phylogenetic trees using data from WALS and
Ethnologue data and the Neighbour Joining algorithm.  It relies on a number
of third party libraries, some of which are also included here.

EXECUTABLE SCRIPTS:

	perform-calibration:
		Spends a lot of time trying to find optimal values of
		parameters for the various distance methods.  Should be run
		before running generate_trees.

	binarizer:
		Reads in a Newick formatted non-binary tree and produces as
		output a Newick formatted tree which is "equivalently binary",
		by adding in branches of length zero.  I don't believe this
		is currently being used by anything.  It was originally
		written to make the authoritative trees (which are non-binary)
		compatible with the inference code (which assumes a binary
		tree).

	simplify:
		Reads in a Newick formatted binary tree and produces as output
		a file describing that tree in a format which is much easier
		to parse in C, which is used by the inference code.

	generate_trees:
		Does about what you'd expect: generates a very large number
		of trees, for all combinations of language family and distance
		method.

	make_demo_trees:
		Copies the first 10 of the randomly generated IE trees into
		a separate directory and changes them slightly so they include
		human readable names instead of integers as branch labels.
		These files can then be used to create pretty tree images for
		use in papers etc.

LIBRARY MODULES:

	fileio:
		Boring file I/O stuff.

	distance:
		Implements the various distance method functions.

INPUT DATA:

	authoritative_trees:
		Authoritative trees for the Austronesian and Indo-European
		language families.  Don't touch anything here or else the
		calibration script will be totally broken.

OUTPUT DATA:

	generated_trees:
		Randomly generated trees.

	demo_trees:
		Pretty IE tree files.
