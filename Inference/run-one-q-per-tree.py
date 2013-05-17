#!/usr/bin/env python
import random
import os

for treeindex in range(1,101):
	for family in "austro afro indo niger nilo sino".split():
		for type in "geographic genetic feature combination".split():
			os.system("./inference.exe -b 500 -s 5000 -t ../TreeBuilder/generated_trees/%s/%s/tree_%d.simple -l ../TreeBuilder/generated_trees/%s.leafdata -o results/individual-q/%s/%s/tree_%d" % (type, family, treeindex, family, type, family, treeindex))
