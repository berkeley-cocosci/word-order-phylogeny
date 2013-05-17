#!/usr/bin/env python

from codecs import open
import os
from sys import exit
import sqlite3
from random import shuffle, gauss, sample, normalvariate, lognormvariate
from copy import deepcopy

import fileio
from distance import *
import dendropy

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
		
def get_languages_by_family(conn, cursor, family):
    ethnoclasses = load_ethnologue_classifications()
    cursor.execute("""SELECT wals_code FROM languages WHERE family=? AND iso_codes != ''""", (family,))
    codes = [code[0] for code in cursor.fetchall()]
    languages = map(lambda(x): language_from_wals_code(conn, cursor, x), codes)
    languages = map(lambda(x): apply_ethnoclass(x, ethnoclasses), languages)
    languages = filter(lambda(x): x.data["ethnoclass"].split(",")[0].strip() == x.data["family"], languages)
    languages = filter(lambda(x): x.data.get(bwo, None) not in (7,'7',None), languages)
    hierlengths = [len(x.data["ethnoclass"].split(",")[1:]) for x in languages ]
    return languages

def language_from_wals_code(conn, cursor, code):
    lang = Language()
    cursor.execute('''SELECT * FROM languages WHERE wals_code=?''',(code,))
    results = cursor.fetchone()
    lang.code = results[0]
    lang.name = results[1].replace(" ","_").replace("(","").replace(")","")
    lang.data = {}
    lang.data["location"] = (float(results[2]),float(results[3]))
    lang.data["genus"] = results[4]
    lang.data["family"] = results[5]
    lang.data["subfamily"] = results[6]
    lang.data["iso_codes"] = results[7]
    cursor.execute('''SELECT name, value_id FROM speedyfeatures WHERE wals_code=?''',(code,))
    for x in cursor.fetchall():
        name, value = x
        lang.data[name] = value
    return lang

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

    base_matrix = build_matrix_by_method_name(languages, build_method)
    for index in range(0, tree_count):
        make_tree(base_matrix, age_params, build_method, family_name, index)

def make_tree(base_matrix, age_params, build_method, family_name, index):

    # Age params is a tuple of mean ages
    # Choose one mean at random (equally weighted mixture model)
    # Standard deviation is always 20% of the mean
    shuffle(age_params)
    age = gauss(age_params[0], age_params[0]*0.25/2)

    matrix = deepcopy(base_matrix)
    failures = 0
    while failures < 5:
        # Add random noise
        for j in range(0, len(matrix)):
            for k in range(j+1, len(matrix)):
                matrix[j][k] = max(0.0, matrix[j][k] + gauss(0, matrix[j][k]*0.1))
                matrix[k][j] = matrix[j][k]

        # Randomly generate a ratio of highest pairwise distanc to lowest,
        # based on the apparent distribution in authoritative trees
        r = gauss(12.4, 1.04)
        print "MY RATIO IS ", r
        minn = min([min(row) for row in matrix])
        maxx = max([max(row) for row in matrix])
        offset = (maxx - r*minn) / (r-1)

        # Normalise the matrix
        norm = offset + max([max(row) for row in matrix])
        for j in range(0, len(matrix)):
            for k in range(j+1, len(matrix)):
                matrix[j][k] += offset
                matrix[j][k] /= norm
                matrix[k][j] = matrix[j][k]

        # Save the matrix and generate a tree from it
        filename = os.path.join("generated_trees", build_method, family_name, "tree_%d" % (index+1))
        directory = os.path.dirname(filename)
        if not os.path.exists(directory):
            os.makedirs(directory)
        fileio.save_matrix(matrix, filename+".distance")

        # Use NINJA to do Neighbour Joining
        os.system("java -jar ./ninja/Ninja.jar --in_type d ./%s.distance > %s.phylip" % (filename, filename))

        # Read output of NINJA into Dendropy Tree
        fp = open("%s.phylip" % filename, "r")
        tree_string = fp.read()
        fp.close()
        tree = dendropy.Tree.get_from_string(tree_string, "newick")

        # Die on negative branch length
        if has_negative_branches(tree):
            try:
                fix_negative_branches(tree)
            except:
                failures += 1
                continue

        break

    if failures == 5:
        # We tried 5 times to build a tree without negatives and failed
        print "Dying due to persistent negative branch problem!"
        exit(42)

    # Root at midpoint
    tree.reroot_at_midpoint()

    # Scale to fit age.
    maxlength = max([leaf.distance_from_root() for leaf in tree.leaf_nodes()])
    desiredmax = age/10000.0
    scalefactor = desiredmax / maxlength 
    tree.scale_edges(scalefactor)

    # Write newly scaled and rooted tree out to Newick file
    fp = open("%s.phylip" % filename, "w")
    fp.write(tree.as_newick_string())
    fp.close()

    # Use simplification script to translate Newick file to Simple file
    os.system("python simplify.py -i %s.phylip -o %s.simple" % (filename, filename))

    # Save tree age
    fp = open("%s.age" % filename, "w")
    fp.write("%f" % age)
    fp.close()
    
    # Clean up after ourselves...
    #os.unlink("%s.distance" % filename)
    #os.unlink("%s.phylip" % filename)

def report_on_dense_langs(languages):
    family_counts = {}
    word_order_counts = {}
    for lang in languages:
        family_counts[lang.data["family"]] = family_counts.get(lang.data["family"], 0) + 1 
        word_order_counts[lang.data[bwo]] = word_order_counts.get(lang.data[bwo], 0) + 1 
    sortedfamily = [(family_counts[family], family) for family in family_counts]
    sortedfamily.sort()
    sortedfamily.reverse()
    fp = codecs.open("dense_lang_report", "w", "UTF-8")
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
    fp = open("../EthnoScrape/classifications.txt")
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
        if isocode in ethnoclasses:
            lang.data["ethnoclass"] = ethnoclasses[isocode]
            return lang

    lang.data["ethnoclass"] = "Unclassified"
    return lang

def main():

    conn = sqlite3.connect("../WALS2SQL/wals.db")
    cursor = conn.cursor()
    cursor.execute('''PRAGMA cache_size = -25000''')
    cursor.execute('''CREATE TEMPORARY TABLE speedyfeatures AS SELECT name, value_id, wals_code FROM data_points INNER JOIN dense_features on data_points.feature_id = dense_features.id''')
    cursor.execute('''CREATE INDEX wals_code_index ON speedyfeatures(wals_code)''')


    afrolangs = get_languages_by_family(conn, cursor, "Afro-Asiatic")
    austrolangs = get_languages_by_family(conn, cursor, "Austronesian")
    indolangs = get_languages_by_family(conn, cursor, "Indo-European")
    nigerlangs = get_languages_by_family(conn, cursor, "Niger-Congo")
    nilolangs = get_languages_by_family(conn, cursor, "Nilo-Saharan")
    sinolangs = get_languages_by_family(conn, cursor, "Sino-Tibetan")
    cursor.close()
    conn.close()

#    report_on_dense_langs(languages)
    languages = (afrolangs, austrolangs, indolangs, nigerlangs, nilolangs, sinolangs)
    ages = ([25000,], [7000,], [6000, 8750], [17500,], [17500,], [7500,])
    names = ("afro", "austro", "indo", "niger", "nilo", "sino")


    if not os.path.exists("generated_trees"):
        os.mkdir("generated_trees")

    for langs, name in zip(languages, names):
        fileio.save_translate_file(langs, "generated_trees/" + name + ".translate")
        fileio.save_multistate_file(langs, "generated_trees/" + name + ".leafdata")

    for method in ("genetic", "geographic", "feature", "combination"):
        for langs, age, name in zip(languages, ages, names):
            make_trees(langs, age, method, name, 100)

if __name__ == "__main__":
    main()
