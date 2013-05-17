#!/usr/bin/env python
import random
import os

for treeindex in range(65,101):
    for treeclass, type in enumerate("geographic genetic feature combination".split()):
			os.system("./inference.exe -b 500 -s 5000 -m -i %d -c %d -o results/common-q/%s/trees_%d" % (treeindex, treeclass, type, treeindex))
