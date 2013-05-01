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
    cursor.execute('''SELECT wals_code FROM languages WHERE family=?''', (family,))
    codes = [code[0] for code in cursor.fetchall()]
    languages = map(lambda(x): language_from_wals_code(conn, cursor, x), codes)
    languages = filter(lambda(x): x.data.get(bwo, None) not in (7,'7',None), languages)
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
    lang.data["isocode"] = results[7]
    cursor.execute('''SELECT name, value_id FROM speedyfeatures WHERE wals_code=?''',(code,))
    for x in cursor.fetchall():
        name, value = x
        lang.data[name] = value
    return lang

def make_trees(languages, age_params, build_method, family_name, treecount):

    base_matrix = make_matrix_go_now(languages, build_method)

    for index in range(0, treecount):
        print index
        if type(age_params) is tuple:
            age = gauss(age_params[0], age_params[1])
        else:
            age = gauss(age_params, age_params*0.15/2)

        # Add random noise
        matrix = deepcopy(base_matrix)
        noisevar = NOISES[build_method]
        for j in range(0, len(languages)):
            for k in range(j+1, len(languages)):
                matrix[j][k] += gauss(0, noisevar)

        # Normalise the matrix
	norm = max([max(row) for row in matrix])
        for j in range(0, len(languages)):
            for k in range(j+1, len(languages)):
                matrix[j][k] /= norm
                matrix[k][j] = matrix[j][k]

        # Save the matrix and generate a tree from it
        filename = os.path.join("generated_trees", build_method, family_name, "tree_%d" % (index+1))
        directory = os.path.dirname(filename)
        if not os.path.exists(directory):
            os.makedirs(directory)
        fileio.save_matrix(matrix, languages, filename+".distance")

        # Use NINJA to do Neighbour Joining
        os.system("java -server -Xmx2G -jar ./ninja/Ninja.jar --in_type d ./%s.distance > %s.phylip" % (filename, filename))

        # Read output of NINJA into Dendropy Tree
        fp = open("%s.phylip" % filename, "r")
        tree_string = fp.read()
        fp.close()
        tree = dendropy.Tree.get_from_string(tree_string, "newick")

	# Die on negative branch length
	for edge in tree.get_edge_set():
		if edge.length and edge.length < 0:
			print "Negative branch length!  Dying!"
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

def main():

    conn = sqlite3.connect("../WALS2SQL/wals.db")
    cursor = conn.cursor()
    cursor.execute('''PRAGMA cache_size = -25000''')
    cursor.execute('''CREATE TEMPORARY TABLE speedyfeatures AS SELECT name, value_id, wals_code FROM data_points INNER JOIN dense_features on data_points.feature_id = dense_features.id''')
    cursor.execute('''CREATE INDEX wals_code_index ON speedyfeatures(wals_code)''')

    austrolangs = get_languages_by_family(conn, cursor, "Austronesian")
    afrolangs = get_languages_by_family(conn, cursor, "Afro-Asiatic")
    indolangs = get_languages_by_family(conn, cursor, "Indo-European")
    nigerlangs = get_languages_by_family(conn, cursor, "Niger-Congo")
    nilolangs = get_languages_by_family(conn, cursor, "Nilo-Saharan")
    sinolangs = get_languages_by_family(conn, cursor, "Sino-Tibetan")
    cursor.close()
    conn.close()

#    report_on_dense_langs(languages)
    languages = (afrolangs, austrolangs, indolangs, nigerlangs, nilolangs, sinolangs)
    ages = (155000, 6000, (7000, 0.2), 10000, 175000, 6000)
    names = ("afro", "austro", "indo", "niger", "nilo", "sino")

    if not os.path.exists("generated_trees"):
        os.mkdir("generated_trees")

    for langs, name in zip(languages, names):
        fileio.save_translate_file(langs, "generated_trees/" + name + ".translate")
        fileio.save_multistate_file(langs, "generated_trees/" + name + ".leafdata")

    for method in "geographic", "genetic", "feature":
        for langs, age, name in zip(languages, ages, names):
            make_trees(langs, age, method, name, 100)

if __name__ == "__main__":
    main()
