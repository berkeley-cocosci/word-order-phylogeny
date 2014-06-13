all: inference generate 

inference:
	gcc -g -O3  main.c beliefprop.c matrix.c mcmc.c modellike.c tree.c gslworkspace.c sampman.c ubertree.c -o inference.exe -lgsl -lgslcblas -lm

generate:
	gcc -g  generate.c beliefprop.c matrix.c mcmc.c modellike.c tree.c gslworkspace.c -o generate.exe -lgsl -lgslcblas -lm

distill:
	gcc -g -lgsl -lgslcblas -lm distiller.c beliefprop.c matrix.c mcmc.c modellike.c -o distill.exe
