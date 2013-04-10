#!/usr/bin/env python
import random
import os

for treeindex in range(65,101):
	for treeclass, type in enumerate("distance family feature".split()):
			os.system("./go -b 500 -s 5000 -m -i %d -c %d -o results/results_common_%s_%d" % (treeindex, treeclass, type, treeindex))
