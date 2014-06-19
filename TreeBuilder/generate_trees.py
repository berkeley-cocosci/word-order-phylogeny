#!/usr/bin/env python

from codecs import open
import os
from sys import exit
import sys
import sqlite3
from random import shuffle, gauss, sample, normalvariate, lognormvariate
from copy import deepcopy

import fileio
from distance import *
import dendropy

sys.path.append("../WALS2SQL/")
import wals2sql

NOISES = {}
NOISES["geographic"] = 0.05 # was 001
NOISES["genetic"] = 0.05
NOISES["feature"] = 0.05

bwo = "Order of Subject, Object and Verb"

class Language:

    def __init__(self):
	self.code = "FOO"
	self.name = "BAR"
	self.data = {}
		
def has_negative_branches(tree):
    for edge in tree.get_edge_set():
        if edge.length and edge.length < 0:
            return True
    return False

def fix_negative_branches(tree):
    while(has_negative_branches(tree)):
        for edge in tree.get_edge_set():
            if edge.length and edge.length < 0:
                parent = edge.head_node
                negative_child = edge.tail_node
                for child in parent.child_nodes():
                    if child != negative_child:
                        other_child = child
                if other_child.edge_length > abs(edge.length):
                    # We can fix this!
                    other_child.edge.length -= edge.length
                    edge.length = 0.0
                    print "********** FIXED A NEGATIVE!!! **********"
                else:
                    # Can't fix this
                    raise Exception("Couldn't fix negative!")
        return False

def make_trees(languages, age_params, build_method, family_name, tree_count):

    # Build base matrix according to specified method
    base_matrix = build_matrix_by_method_name(languages, build_method)
    assert base_matrix.min() >=0

    # Build random variations of base matrix
    for index in range(0, tree_count):
        failures = 0
        filename = os.path.join("generated_trees", "whole", build_method, family_name, "tree_%d" % (index+1))
        while failures < 10:
            try:
                make_tree(base_matrix, age_params, filename)
                break
            except:
                failures += 1
                continue
        if failures == 10:
            print "Aborting after sustained failure to generated %d-th whole tree with the following parameters:" % index
            print build_method, family_name, tree_count
            exit(42)

def make_split_trees(languages, age_params, build_method, family_name, tree_count):
    for index in range(0, tree_count):
        failures = 0
        while failures < 10:
            try:
                # Split the languages up
                shuffle(languages)
                pivot = len(languages)/2
                langs1 = languages[0:pivot]
                langs2 = languages[pivot:]
                for half, langs in enumerate((langs1, langs2)):
                    filename = os.path.join("generated_trees", "split", str(half+1), build_method, family_name, "tree_%d" % (index+1,))
                    base = build_matrix_by_method_name(langs, build_method)
                    make_tree(base, age_params, filename, preserve_base=False)
                break
            except:
                failures += 1
                continue
        if failures == 10:
            print "Aborting after sustained failure to generated %d-th split tree with the following parameters:" % index
            print build_method, family_name, tree_count
            exit(42)

def make_tree(base_matrix, age_params, filename, preserve_base=True):

    if preserve_base:
        matrix = deepcopy(base_matrix)
    else:
        matrix = base_matrix

    # Add random noise
    # The biggest distance is now 1.0
    # Let's set our Gaussian parameters so that we very rarely add/take more than 0.2
    # to/from a distance.  3sigma = 0.2 => sigma = 0.0666
    for j in range(0, len(matrix)):
        for k in range(j+1, len(matrix)):
            matrix[j][k] = max(0.0, matrix[j][k] + gauss(0, 0.0666))
            matrix[k][j] = matrix[j][k]

    # Renormalise matrix
    norm = max([max(row) for row in matrix])
    for j in range(0, len(matrix)):
        for k in range(j+1, len(matrix)):
            matrix[j][k] /= norm
            matrix[k][j] = matrix[j][k]

    # Save the matrix and generate a tree from it
    directory = os.path.dirname(filename)
    if not os.path.exists(directory):
        os.makedirs(directory)
    fileio.save_matrix(matrix, filename+".distance")

    # Use NINJA to do Neighbour Joining
    os.system("java -jar ./ninja/Ninja.jar --in_type d ./%s.distance > %s.tree" % (filename, filename))

    # Read output of NINJA into Dendropy Tree
    fp = open("%s.tree" % filename, "r")
    tree_string = fp.read()
    fp.close()
    tree = dendropy.Tree.get_from_string(tree_string, "newick")

    # Fix any negative branches
    if has_negative_branches(tree):
        fix_negative_branches(tree)

    # Root at midpoint
    tree.reroot_at_midpoint()

    # Sample change range parameters, and adjust branch lengths
    # so they reflect *only* elapsed time, in preparation for age scaling
    change_rates = {}
    for edge in tree.get_edge_set():
        if edge.tail_node:
            rate = lognormvariate(0,0.25)
            change_rates[edge.oid] = rate
            edge.length /= rate

    # Sample tree age
    # age_params is a tuple of mean ages
    # Choose one mean at random (equally weighted mixture model)
    # Standard deviation is always 20% of the mean
    shuffle(age_params)
    age = gauss(age_params[0], age_params[0]*0.25/2)

    # Scale to fit age
    maxlength = max([leaf.distance_from_root() for leaf in tree.leaf_nodes()])
    desiredmax = age/10000.0
    scalefactor = desiredmax / maxlength 
    tree.scale_edges(scalefactor)

    # Reapply change rates
    for edge in tree.get_edge_set():
        if edge.oid in change_rates:
            edge.length *= change_rates[edge.oid]

    # Write newly scaled and rooted tree out to Newick file
    fp = open("%s.tree" % filename, "w")
    fp.write(tree.as_newick_string())
    fp.close()

    # Use simplification script to translate Newick file to Simple file
    os.system("python simplify.py -i %s.tree -o %s.simple" % (filename, filename))

    # Save tree age
    fp = open("%s.age" % filename, "w")
    fp.write("%f\n" % age)
    fp.close()
    
    # Clean up after ourselves...
    os.unlink("%s.distance" % filename)
    os.unlink("%s.tree" % filename)

def report_on_dense_langs(languages):
    family_counts = {}
    word_order_counts = {}
    for lang in languages:
        family_counts[lang.data["family"]] = family_counts.get(lang.data["family"], 0) + 1 
        word_order_counts[lang.data[bwo]] = word_order_counts.get(lang.data[bwo], 0) + 1 
    sortedfamily = [(family_counts[family], family) for family in family_counts]
    sortedfamily.sort()
    sortedfamily.reverse()
    fp = open("dense_lang_report", "w", "UTF-8")
    fp.write("Total languages: %d\n" % len(languages))
    fp.write("Family breakdown:\n")
    for count, family in sortedfamily:
        fp.write("%s: %d\n" % (family, count))
    fp.write("----------\n")
    fp.write("Word order breakdown:\n")
    for order in word_order_counts:
        fp.write("%d: %d\n" % (order, word_order_counts[order]))
    fp.write("----------\n")
    fp.close()

def load_ethnologue_classifications():
    ethnoclass = {}
    fp = open("../EthnoScrape/classifications.txt", "UTF-8")
    for line in fp:
        if not line or line.startswith("#"):
            continue
        code, classification = [x.strip() for x in line.split(":")]
        codes = code.split()
        for code in codes:
            ethnoclass[code] = classification
    fp.close()
    return ethnoclass

def apply_ethnoclass(lang, ethnoclasses):

    if lang.code == "nku":
        lang.data["iso_codes"] = "xnz"
    elif lang.code == "nbd":
        lang.data["iso_codes"] = "dgl"
    for isocode in lang.data["iso_codes"].split():
        # Patch WALS's outdated ISO codes...
        if isocode == "gmo":
            isocode = "gmv"
        elif isocode == "lav":
            isocode = "lvs"
        elif isocode == "nep":
            isocode = "npi"
        elif isocode == "nob":
            isocode = "nor"
        elif isocode == "ori":
            isocode = "ory"
        elif isocode == "daf":
            isocode = "dnj"
# I am pretty sure now (07/05/13) that changing dan to dnj is
# wrong (sends Denmark to Africa!), but I'm commenting rather
# than deleting because it's possible I meant to change something
# else to dnj and having this around as a clue may pay off...
#                elif isocode == "dan":
#                    isocode = "dnj"
        elif isocode == "izi":
            isocode = "izz"
        elif isocode == "jar":
            isocode = "anq"
        elif isocode == "baz":
            isocode = "tvu"
        elif isocode == "kln":
            isocode = "niq"
        # Australian / TNG languages
        elif isocode == "gbc":
            isocode = "wrk"
        elif isocode == "mwd":
            isocode = "dmw"
        elif isocode == "nbx":
            isocode = "xwk"
        elif isocode == "wiw":
            isocode = "wgu"
        elif isocode == "unp":
            isocode = "wro"
        # American languages
        elif isocode == "wit":
            if lang.code == "ptw":
                isocode = "pwi"
            elif lang.code == "win":
                isocode = "wnw"
        elif isocode == "tzc":
            isocode = "tzo"
        elif isocode == "ckf":
            isocode = "cak"
        elif isocode == "hsf":
            isocode = "hus"
        elif isocode == "ixi":
            isocode = "ixl"
        elif isocode == "mvc":
            isocode = "mam"
        elif isocode == "chb":
            isocode = "EXTINCT"

        if isocode in ethnoclasses:
            lang.data["ethnoclass"] = ethnoclasses[isocode].decode("UTF-8")
            if lang.data["family"] in ("Oto-Manguean", "Uto-Aztecan", "Mayan", "Algic", "Arawakan", "Salishan", "Penutian", "Na-Dene", "Tupian", "Hokan", "Cariban", "Tucanoan", "Eskimo-Aleut", "Macro-Ge", "Chibchan"):
                if lang.data["ethnoclass"].startswith(lang.data["family"]):
                    classification = "American, " + lang.data["ethnoclass"]
                else:
                    classification = ("American, %s, " % lang.data["family"]) + lang.data["ethnoclass"]
                lang.data["ethnoclass"] = classification

            return lang

    lang.data["ethnoclass"] = "Unclassified"
#    print lang.data.get("Order of Subject, Object and Verb",None)
    if lang.data.get("Order of Subject, Object and Verb",None) not in (None,7,'7'):
        fp = open("/home/luke/missing_codes.txt", "a")
        fp.write(lang.data["iso_codes"]+"\n")
        fp.close()
    lang.data["ethnoclass"] = lang.data["ethnoclass"].decode("UTF-8")
    return lang

def main():

    conn = sqlite3.connect("../WALS2SQL/wals.db")
    cursor = conn.cursor()
    cursor.execute('''PRAGMA cache_size = -25000''')
    cursor.execute('''CREATE TEMPORARY TABLE speedyfeatures AS SELECT name, value_id, wals_code FROM data_points INNER JOIN dense_features on data_points.feature_id = dense_features.id''')
    cursor.execute('''CREATE INDEX wals_code_index ON speedyfeatures(wals_code)''')


    afrolangs = wals2sql.get_dense_languages_by_family(conn, cursor, "Afro-Asiatic")
    austlangs = wals2sql.get_dense_languages_by_family(conn, cursor, "Australian")
    austrolangs = wals2sql.get_dense_languages_by_family(conn, cursor, "Austronesian")
    indolangs = wals2sql.get_dense_languages_by_family(conn, cursor, "Indo-European")
    nigerlangs = wals2sql.get_dense_languages_by_family(conn, cursor, "Niger-Congo")
    nilolangs = wals2sql.get_dense_languages_by_family(conn, cursor, "Nilo-Saharan")
    sinolangs = wals2sql.get_dense_languages_by_family(conn, cursor, "Sino-Tibetan")
    tnglangs = wals2sql.get_dense_languages_by_family(conn, cursor, "Trans-New Guinea")
    americanlangs = []
    for amfam in ("Oto-Manguean", "Uto-Aztecan", "Mayan", "Algic", "Arawakan", "Salishan", "Penutian", "Na-Dene", "Tupian", "Hokan", "Cariban", "Tucanoan", "Eskimo-Aleut", "Macro-Ge", "Chibchan"):
        americanlangs.extend(wals2sql.get_dense_languages_by_family(conn, cursor, amfam))
    cursor.close()
    conn.close()

#    report_on_dense_langs(languages)
    languages = (afrolangs, austrolangs, indolangs, nigerlangs, nilolangs, sinolangs, austlangs, tnglangs, americanlangs)
    ages = ([25000,], [7000,], [6000, 8750], [17500,], [17500,], [7500,], [50000,], [4250,], [16500,])
    names = ("afro", "austro", "indo", "niger", "nilo", "sino", "aust", "tng", "amer")

    if not os.path.exists("generated_trees"):
        os.mkdir("generated_trees")

    for langs, name in zip(languages, names):
        fileio.save_translate_file(langs, "generated_trees/" + name + ".translate")
        fileio.save_multistate_file(langs, "generated_trees/" + name + ".leafdata")

    for method in ("geographic", "genetic", "feature", "combination"):
        for langs, age, name in zip(languages, ages, names):
            make_trees(langs, age, method, name, 100)
            make_split_trees(langs, age, method, name, 100)

if __name__ == "__main__":
    main()
