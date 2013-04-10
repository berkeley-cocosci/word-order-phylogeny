all: wals trees inference

wals:
	cd WALS2SQL; ./wals2sql.py

trees:
	cd TreeBuilder; ./generate_trees.py

inference: 
	cd Inference; gcc -lgsl -lgslcblas -lm -g main.c -o inference.exe

