import random
import os

for treeindex in range(1,101):
	for family in "indo austro afro niger nilo sino".split():
		for type in "distance family feature".split():
			os.system("./go -b 500 -s 5000 -t ../TreeBuilder/trees/%s%s%d.simple -l ../TreeBuilder/trees/%s.leafdata -o results/results_%s_%s_%d" % (family, type, treeindex, family, family, type, treeindex))
