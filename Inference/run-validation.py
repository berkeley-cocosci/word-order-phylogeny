#!/usr/bin/env python
import random
import os

for treeindex in range(1,21):
	for family in "indo austro afro niger nilo sino".split():
		for type in "distance family feature".split():
			rootorder = random.sample("SOV SVO VSO VOS OVS OSV".split(), 1)[0]
			os.system("./generate -t ../TreeBuilder/trees/%s%s%d.simple -l ../TreeBuilder/trees/%ssleafdata -r %s -x -o validation/generated_%s%s%d" % (family, type, treeindex, family, rootorder, family, type, treeindex))
			os.system("./go -b 500 -s 5000 -t ../TreeBuilder/trees/%s%s%d.simple -l validation/generated_%s%s%d_leaves -o validation/recovered_%s%s%d" % (family, type, treeindex, family, type, treeindex, family, type, treeindex))
