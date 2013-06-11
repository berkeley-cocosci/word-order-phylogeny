#!/usr/bin/env python -u

import sys

values = []
for line in sys.stdin:
    try:
        values.append(float(line.strip()))
    except:
        continue

print "Min: ", min(values)
print "Mean: ", sum(values) / len(values)
print "Max: ", max(values)
