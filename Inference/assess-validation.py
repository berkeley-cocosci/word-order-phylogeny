import random
import os
from math import sqrt


def get_ancestral_dist(filename):
    ancestral = [0, 0, 0, 0, 0, 0]
    n = 0
    fp = open(filename, "r")
    for line in fp:
        n += 1
        anc = map(float,line.strip().split())
        for i in range(0,6):
            ancestral[i] += anc[i]
    fp.close()
    ancestral = [p/n for p in ancestral]
    return ancestral

def load_true_root(filenameroot):
    filename = "validation/generated_%s_matrix" % filenameroot
    fp = open(filename, "r")
    root = int(fp.readline().strip())
    fp.close()
    return root

def load_recovered_root(filenameroot):    
    filename = "validation/recovered_%s/ancestrals" % filenameroot
    ancdist = get_ancestral_dist(filename)
    ancdist = zip(ancdist, range(0,6))
    ancdist.sort()
    ancdist.reverse()
    return [x[1] for x in ancdist[0:2]]

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
    filename = "validation/recovered_%s/summary" % filenameroot
    fp = open(filename, "r")
    matrix = []
    for i in range(0,12):
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
near = 0
tries = 0
for treeindex in range(1,5):
    for family in "indo austro afro niger nilo sino".split():
        for type in ("combination",):
                true = load_true_root("%s%s%d" % (family, type, treeindex))
                first, second = load_recovered_root("%s%s%d" % (family, type, treeindex))
                if true == first:
                    hits += 1
                elif true == second:
                    near += 1
                tries += 1
                print "Truth was %d, I guessed %d or %d" % (true, first, second)
                true = load_true_q("%s%s%d" % (family, type, treeindex))
                guess = load_recovered_q("%s%s%d" % (family, type, treeindex))
                truevec = []
                for row in true:
                    truevec.extend(row)
                guessvec = []
                for row in guess:
                    guessvec.extend(row)
                avgr += correl(truevec, guessvec)
print "Correct roots: ", 1.0*hits/tries
print "Nearly correct roots: ", 1.0*(hits+near)/tries
print "Mean Q correlation: ", avgr/tries
