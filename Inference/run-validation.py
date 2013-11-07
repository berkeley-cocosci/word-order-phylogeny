#!/usr/bin/env python
import random
import os

for treeindex in range(1,5):
    for family in "indo austro afro niger nilo sino".split():
        for type in ("combination",):
            rootorder = random.sample("SOV SVO VSO VOS OVS OSV".split(), 1)[0]
            os.system("./generate.exe -t ../TreeBuilder/generated_trees/%s/%s/tree_%d.simple -l ../TreeBuilder/generated_trees/%s.leafdata -r %s -x -o validation/generated_%s%s%d" % (type, family, treeindex, family, rootorder, family, type, treeindex))
            os.system("./inference.exe -b 500 -s 5000 -t ../TreeBuilder/generated_trees/%s/%s/tree_%d.simple -l validation/generated_%s%s%d_leaves -o validation/recovered_%s%s%d" % (type, family, treeindex, family, type, treeindex, family, type, treeindex))
            os.system("./distill.exe -d validation/recovered_%s%s%d/" % (family, type, treeindex))
