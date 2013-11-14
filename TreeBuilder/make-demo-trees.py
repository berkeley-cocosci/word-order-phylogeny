#!/usr/bin/env python

# Super gross script to make some demonstration IE trees with
# human friendly names

import glob
from fileio import load_translate_file

revtrans = load_translate_file("generated_trees/indo.translate")
names = revtrans.keys()
trans = {}
for name in names:
    trans[revtrans[name]] = name

for i in range(1,11):
    fp = open("generated_trees/combination/indo/tree_%d.tree" % i, "r")
    string = fp.read()
    fp.close()
    for number in trans:
        print "Replacing ", number, " with ", trans[number]
        string = string.replace("(%d:" % number, "(%s:" % trans[number])
        string = string.replace(",%d:" % number, ",%s:" % trans[number])
    fp = open("demo-trees/demo_indo_%d.tree" % i, "w")
    fp.write(string)
    fp.close()

backtrans = load_translate_file("generated_trees/indo.translate")
trans = {}
for key in backtrans:
    trans[backtrans[key]] = key

for filename in glob.glob("demo-trees/*.tree"):
    fp = open(filename, "r")
    text = fp.read()
    fp.close()

    for number in trans:
        text = text.replace("(%d:" % number, "(%s:" % trans[number])
        text = text.replace(",%d:" % number, ",%s:" % trans[number])

    fp = open(filename.replace(".tree","_named.tree"), "w")
    fp.write(text)
    fp.close()

