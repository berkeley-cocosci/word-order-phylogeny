all: inference generate distill ubertree

inference:
	gcc -Wall -g -lgsl -lgslcblas -lm main.c beliefprop.c matrix.c mcmc.c modellike.c saveresults.c tree.c sampman.c -o inference.exe

generate:
	gcc -g -lgsl -lgslcblas -lm generate.c beliefprop.c matrix.c mcmc.c modellike.c tree.c -o generate.exe

distill:
	gcc -g -lgsl -lgslcblas -lm distiller.c beliefprop.c matrix.c mcmc.c modellike.c -o distill.exe

ubertree:
	gcc -g -lgsl -lgslcblas -lm ubertree.c beliefprop.c matrix.c mcmc.c modellike.c -o ubertree.exe
