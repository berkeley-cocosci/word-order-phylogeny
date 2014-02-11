#!/usr/bin/env python
from math import floor
import multiprocessing
import time
import pdb
import random
import os

families = "indo austro afro niger nilo sino".split()

rootstabs = [0.10,0.20,0.30,0.40,0.50,0.75,1.00,1.25,1.50,1.75,
             2.00,2.50,3.00,3.50,4.00,4.50,5.00,6.00,7.00,8.0]
log = []
stabs = []
trees = 20

def gen_and_inf(treeindex, family, type):
    # Generate leaf data
    rootorder = random.sample("SOV SVO VSO VOS OVS OSV".split(), 1)[0]
    rootstab = rootstabs[int(floor((treeindex-1)/(trees/20)))]
    os.system("./generate.exe -t ../TreeBuilder/generated_trees/whole/%s/%s/tree_%d.simple -l ../TreeBuilder/generated_trees/%s.leafdata -r %s -s %f -x -o validation/generated_%s%s%d" % (type, family, treeindex, family, rootorder, rootstab, family, type, treeindex))
    # Infer ancestors
    os.system("./inference.exe -b 50 -s 10 -T ../TreeBuilder/generated_trees/whole/%s/%s/tree_%d.simple -L validation/generated_%s%s%d_leaves -o validation/recovered_%s%s%d" % (type, family, treeindex, family, type, treeindex, family, type, treeindex))

# Build inference argument list
allargs = []
for treeindex in range(1,trees):
    random.shuffle(families)
    for family in families:
        for type in ("combination",):
            allargs.append((treeindex, family, type))

# Set 8 processes running
procs = []
#for i in range(0,8):
#    proc = multiprocessing.Process(target=gen_and_inf, args=allargs.pop())
#    procs.append(proc)
#    proc.start()

# Monitor those 8 and start a new proc each time one finishes
while procs:
    time.sleep(1)
    proc = procs.pop()
    if proc.exitcode != None:
        if allargs:
            proc = multiprocessing.Process(target=gen_and_inf, args=allargs.pop())
            proc.start()
            procs.insert(0, proc)
    else:
        procs.insert(0, proc)

# See how we did
for treeindex in range(1,trees):
    random.shuffle(families)
    for family in families:
        for type in ("combination",):

            fp = open("validation/generated_%s%s%d_matrix" % (family, type, treeindex))
            lines = fp.readlines()
            fp.close()
            rootorder = int(lines[0].strip())
            lines = lines[2:]
            line = lines[rootorder]
            line = map(float, line.split())
            stab = -1*line[rootorder]
            stabs.append(stab)
            # Compare
            fp = open("validation/recovered_%s%s%d/summary" % (family, type, treeindex))
            lines = fp.readlines()
            fp.close()
            ancestors = map(float, lines[31].split())
            #pdb.set_trace()

            #hokey_prior = [map(float, row.split())[i] for i, row in enumerate(lines[20:26])]
            #hokey_prior = [x**3 for x in hokey_prior]
            #hokey_prior = [x/sum(hokey_prior) for x in hokey_prior]
            #ancestors = [x*y for x,y in zip(hokey_prior, ancestors)]
            #ancestors = [x/sum(ancestors) for x in ancestors]
            print "Ancestors: ", ancestors
            print "True root: ", rootorder
            if ancestors.index(max(ancestors)) == rootorder:
                print "HIT"
                log.append((stab, 1.0))
            else:
                log.append((stab, 0.0))

log.sort()
print log
lower = 0.0
span = 0.1
while True:
    hits = 0
    tries = 0
    for stab, hit in log:
#        print stab, hit
        print "Is %f <= %f < %f?" % (lower, stab, lower+span)
        if lower <= stab < (lower+span):
            print "Yes!"
            hits += hit
            tries += 1.0
        else:
            print "No"
    if lower >= 7.0:
        break
    if tries:
        print lower+span, hits/tries, tries
    lower += span

#print min(stabs), max(stabs), sum(stabs)/len(stabs)
#print hits0 / tries0, hits1/tries1, hits2/tries2, hits3/tries3
