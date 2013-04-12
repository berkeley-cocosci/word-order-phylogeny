all: ethno wals trees inference

ethno:
	cd EthnoScrape; ./ethnoscrape.py

wals:
	cd WALS2SQL; ./wals2sql.py

trees:
	cd TreeBuilder; ./generate_trees.py

inference: 
	cd Inference; gcc -lgsl -lgslcblas -lm -g main.c -o inference.exe

