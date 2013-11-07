word-order-phylogeny
====================

Code for inferring ancestral word orders and word order change dynamics from
present day language data.

The code is relatively modular, with the following major parts:

EthnoScrape:
	A Python script to scrape the Ethnologue website and record genetic
	classifications.

	Requires Beautiful Soup, which it's up to you to install.

WALS2SQL:
	Python code that takes the WALS .csv data files and produces an SQLite
	database that has the information needed for the TreeBuilder code to
	build various classes of estimated phylogenetic trees.

TreeBuilder:
	Python code that uses the results of EthnoScrape and WALS2SQL to build
	pairwise distance matrices between languages.  Neighbour Joining is
	then used to turn these matrices into trees, which are later used as
	samples in an MCMC analysis.

	Requires:
		NINJA, a Java implementation of Neighbour Joining by Travis
		Wheeler (see http://nimbletwist.com/software/ninja/).
		
		Dendropy, a Python library for doing phylogenetic stuff (see
		http://pythonhosted.org/DendroPy/).

		Newick, a Python library for parsing Newick format tree files
		by Thomas Mailund (see
		http://www.daimi.au.dk/~mailund/newick.html).  This could
		probably be replaced with Dendropy to minimise dependencies.

	All the above three projects are included in this repository for
	convenience of distribution.  I believe this is done in compliance
	with the licenses for NINJA and Dendropy.  Strictly speaking, I
	shouldn't be doing this for Newick.  I'll either replace it with
	Dendropy or get the author's blessing in the near future...

Inference:
	Performs MCMC inference over ancestral word orders and word order
	change dynamics, using the trees produced by TreeBuilder.  All of the
	computational heavy lifting is done in C, with convenient Python
	scripts to handle basic logic.

MatlabAnalysis:
	MATLAB code to do some basic post-processing of the results produced
	by Inference and produce pretty pictures.  MATLAB is awful and I would
	love to replace this with Python stuff: Pandas, Matplotlib, etc.
