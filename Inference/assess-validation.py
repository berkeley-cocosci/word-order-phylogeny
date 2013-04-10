import random
import os
from math import sqrt

from meta import get_ancestral_dist

def load_true_root(filenameroot):
	filename = "validation/generated_%s_matrix" % filenameroot
	fp = open(filename, "r")
	root = int(fp.readline().strip())
	fp.close()
	return root

def load_recovered_root(filenameroot):	
	filename = "validation/recovered_%s" % filenameroot
	ancdist = get_ancestral_dist(filename)
	ancdist = zip(ancdist, range(0,6))
	ancdist.sort()
	ancdist.reverse()
	return ancdist[0][1]

def load_true_q(filenameroot):
	filename = "validation/generated_%s_matrix" % filenameroot
        fp = open(filename, "r")
        matrix = []
        for i in range(0,2):
                fp.readline()
        for i in range(0,6):
                matrix.append(map(float, fp.readline().strip().split()))
        fp.close()
        return matrix

def load_recovered_q(filenameroot):
	filename = "validation/recovered_%s" % filenameroot
        fp = open(filename, "r")
        matrix = []
        for i in range(0,6):
                fp.readline()
        for i in range(0,6):
                matrix.append(map(float, fp.readline().strip().split()))
        fp.close()
        return matrix

def correl(x, y):
	x = map(float, x)
	y = map(float, y)
	xbar = sum(x)/len(x)
	ybar = sum(y)/len(y)
	sumnumer = sum([(xi - xbar)*(yi -ybar) for (xi,yi) in zip(x,y)])
	sumxsqu = sum([(xi - xbar)**2 for xi in x])
	sumysqu = sum([(yi - ybar)**2 for yi in y])
	return sumnumer / (sqrt(sumxsqu)*sqrt(sumysqu))

avgr = 0
hits = 0
tries = 0
for treeindex in range(1,101):
	for family in "indo austro afro niger nilo sino".split():
		for type in "distance family feature".split():
			try:
				true = load_true_root("%s%s%d" % (family, type, treeindex))
				guess = load_recovered_root("%s%s%d" % (family, type, treeindex))
				if true == guess:
					hits += 1
				tries += 1
				print "Truth was %d, I guessed %d" % (true, guess)
				true = load_true_q("%s%s%d" % (family, type, treeindex))
				guess = load_recovered_q("%s%s%d" % (family, type, treeindex))
				truevec = []
				for row in true:
					truevec.extend(row)
				guessvec = []
				for row in guess:
					guessvec.extend(row)
				avgr += correl(truevec, guessvec)
			except:
				pass
print "Correct roots: ", 1.0*hits/tries
print "Mean Q correlation: ", avgr/tries
