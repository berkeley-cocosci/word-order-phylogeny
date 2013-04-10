#!/usr/bin/env python

from codecs import open
from os import system, unlink
from sys import exit
import sqlite3
from random import shuffle, gauss, sample, normalvariate, lognormvariate
from copy import deepcopy
from distance import *
import dendropy

FUZZES = {}
FUZZES["family"] = 0.05
FUZZES["distance"] = 0.05 # was 001
FUZZES["feature"] = 0.05

bwo = "Order of Subject, Object and Verb"

class Language:

    def __init__(self):
	self.code = "FOO"
	self.name = "BAR"
	self.data = {}
		
def language_from_wals_code_factory(code, conn, cursor):
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
    cursor.execute('''SELECT name, value_id FROM speedyfeatures WHERE wals_code=?''',(code,))
    for x in cursor.fetchall():
        name, value = x
	#print name
        lang.data[name] = value
    return lang

def get_all_languages(conn, cursor):
    languages = []
    cursor.execute('''SELECT wals_code FROM languages''')
    for code in cursor.fetchall():
         lang = language_from_wals_code_factory(code[0], conn, cursor)
         if (bwo not in lang.data) or (lang.data[bwo] in (None, 7, '7')):
             continue
         else:
             languages.append(lang)
    return languages

def get_most_languages(conn, cursor):
    languages = []
    cursor.execute('''SELECT wals_code FROM dense_languages WHERE family IN ("Afro-Asiatic", "Austronesian", "Indo-European", "Niger-Congo", "Nilo-Saharan", "Sino-Tibetan")''')
    codes = cursor.fetchall()
    for index, code in enumerate(codes):
         print "Constructing lang %d of %d" %(index, len(codes))
         lang = language_from_wals_code_factory(code[0], conn, cursor)
         if (bwo not in lang.data) or (lang.data[bwo] in (None, 7, '7')):
             continue
         else:
             languages.append(lang)
    return languages

def instantiate_dense_language_objects(conn, cursor):
    languages = []
    cursor.execute('''CREATE TABLE IF NOT EXISTS dense_data AS
        SELECT wals_code, feature_id, value_id FROM data_points
        WHERE wals_code IN (SELECT wals_code FROM dense_languages) AND feature_id IN (SELECT
id FROM dense_features)''')

    cursor.execute('''SELECT wals_code FROM dense_languages''')
    dense_codes = [x[0] for x in cursor.fetchall()]
    for i, code in enumerate(dense_codes):
        language = language_from_wals_code_factory(code, conn, cursor)
        if "Order of Subject, Object and Verb" not in language.data:
            print "Dropping %s because it has no word order data!" % language.name.encode("UTF-8")
        elif language.data["Order of Subject, Object and Verb"] == 7:
            print "Dropping %s because it has no basic word order!" % language.name.encode("UTF-8")
        else:
            languages.append(language)
    return languages

def save_translate_file(languages, filename):
    fp = open(filename, "w", "utf8")
    fp.write("translate\n");
    fp.write("%d %s" % (0, languages[0].name)) 
    for i in range(1,len(languages)):
        fp.write(",\n%d %s" % (i, languages[i].name)) 
    fp.write(";\n")
    fp.close()

def save_multistate_file(languages, filename):
    fp = open(filename, "w", "utf8")
    for i, lang in enumerate(languages):
	fp.write("%d	%d\n" % (i, lang.data[u'Order of Subject, Object and Verb']))
    fp.close()

def subsample_langs(all_languages):
    sovlangs = []
    svolangs = []
    vsolangs = []
    voslangs = []
    ovslangs = []
    osvlangs = []
    for lang in all_languages:
	if lang.data[u'Order of Subject, Object and Verb'] == 1:
            sovlangs.append(lang)
	elif lang.data[u'Order of Subject, Object and Verb'] == 2:
            svolangs.append(lang)
	elif lang.data[u'Order of Subject, Object and Verb'] == 3:
            vsolangs.append(lang)
	elif lang.data[u'Order of Subject, Object and Verb'] == 4:
            voslangs.append(lang)
	elif lang.data[u'Order of Subject, Object and Verb'] == 5:
            ovslangs.append(lang)
	elif lang.data[u'Order of Subject, Object and Verb'] == 6:
            osvlangs.append(lang)
	else:
            print "SHIT!"
    print [len(x) for x in (sovlangs, svolangs, vsolangs, voslangs, ovslangs, osvlangs)]
    subsample = []
    subsample.extend(sample(sovlangs, 103))
    subsample.extend(sample(svolangs, 89))
    subsample.extend(sample(vsolangs, 17))
    subsample.extend(sample(voslangs, 5))
    subsample.extend(sample(ovslangs, 2))
    subsample.extend(sample(osvlangs, 0))
    return subsample

def make_trees(all_languages, ages, build_method, filename_basis, subsample=False):

    if subsample:
        languages = subsample_langs(all_languages)
    else:
        languages= all_languages
    matrix = make_matrix_go_now(languages, build_method)

    for index, age in enumerate(ages):
        print index
        # Fuzz the matrix up
        fuzzmatrix = deepcopy(matrix)
        minn = 999999
        maxx = 0
        fuzz_size = FUZZES[build_method]
        for j in range(0, len(languages)):
            for k in range(j, len(languages)):
                if j == k:
                    continue
                fuzzmatrix[j][k] += gauss(0, fuzz_size)
                fuzzmatrix[k][j] = fuzzmatrix[j][k]
                #print "Before fuzzing I had: ", matrix[j][k]
                #print "Now I have: ", fuzzmatrix[j][k]
                if fuzzmatrix[j][k] < minn:
                    minn = fuzzmatrix[j][k]
                if fuzzmatrix[j][k] > maxx:
                    maxx = fuzzmatrix[j][k]

        # Normalise the matrix
        # print "Got min, max: ", minn, maxx
        for j in range(0, len(languages)):
            for k in range(j, len(languages)):
                if j == k:
                    continue
                #print "Pre norm: ", fuzzmatrix[j][k]
                fuzzmatrix[j][k] /= (1.0*maxx)
                #print "Post norm: ", fuzzmatrix[j][k]
                fuzzmatrix[k][j] = fuzzmatrix[j][k]

        # Save the matrix and generate a tree from it
        fname = filename_basis + str(index+1)
        save_matrix(fuzzmatrix, languages, fname+".distance")
        run_diags(fuzzmatrix, languages, fname+".diag")

        # Use NINJA to do Neighbour Joining
        system("java -server -Xmx2G -jar ./ninja/Ninja.jar --in_type d ./%s.distance > %s.phylip" % (fname, fname))

        # Read output of NINJA into Dendropy Tree
        fp = open("%s.phylip" % fname, "r")
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
        fp = open("%s.phylip" % fname, "w")
        fp.write(tree.as_newick_string())
        fp.close()

        # Use simplification script to translate Newick file to Simple file
        system("python simplify.py -i %s.phylip -o %s.simple" % (fname, fname))

        # Clean up after ourselves...
        unlink("%s.distance" % fname)
        unlink("%s.phylip" % fname)

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

def gen_afro_age():
    return normalvariate(15500,1225)

def gen_austro_age():
    return normalvariate(6000,500)

def gen_indo_age():
    return normalvariate(7000,500)

def gen_niger_age():
    return normalvariate(10000,1000)

def gen_nilo_age():
    return normalvariate(17500,1225)

def gen_sino_age():
    return normalvariate(6000, 500)

def main():

    conn = sqlite3.connect("../WALS2SQL/wals.db")
    cursor = conn.cursor()
    cursor.execute('''PRAGMA cache_size = -25000''')
    cursor.execute('''CREATE TEMPORARY TABLE speedyfeatures AS SELECT name, value_id, wals_code FROM data_points INNER JOIN dense_features on data_points.feature_id = dense_features.id''')
    cursor.execute('''CREATE INDEX wals_code_index ON speedyfeatures(wals_code)''')
    #languages = get_all_languages(conn, cursor)
    languages = get_most_languages(conn, cursor)
    #return
    datasizes = []
    for lang in languages:
        datasizes.append(len(lang.data.keys()))
    datasizes = set(datasizes)
    cursor.close()
    conn.close()
    report_on_dense_langs(languages)

    nigerlangs = []
    sinolangs = []
    indolangs = [] 
    austrolangs = [] 
    afrolangs = [] 
    nilolangs = [] 

    for lang in languages:
        if lang.data["family"] == "Niger-Congo":
            nigerlangs.append(lang)
        elif lang.data["family"] == "Sino-Tibetan":
            sinolangs.append(lang)
        elif lang.data["family"] == "Indo-European":
            indolangs.append(lang)
        elif lang.data["family"] == "Austronesian":
            austrolangs.append(lang)
        elif lang.data["family"] == "Afro-Asiatic":
            afrolangs.append(lang)
        elif lang.data["family"] == "Nilo-Saharan":
            nilolangs.append(lang)

    for langs, name in zip((afrolangs, austrolangs, indolangs, nigerlangs, nilolangs, sinolangs),"afro austro indo niger nilo sino".split()):
        save_translate_file(langs, "trees/" + name + ".translate")
        save_multistate_file(langs, "trees/" + name + ".leafdata")

    treecount = 100
    afroages = [gen_afro_age() for i in range(0, treecount)]
    austroages = [gen_austro_age() for i in range(0, treecount)]
    indoages = [gen_indo_age() for i in range(0, treecount)]
    nigerages = [gen_niger_age() for i in range(0, treecount) ]
    niloages = [gen_nilo_age() for i in range(0, treecount)]
    sinoages = [gen_sino_age() for i in range(0, treecount)]

    make_trees(afrolangs, afroages, "feature", "trees/afrofeature")
    make_trees(afrolangs, afroages, "family", "trees/afrofamily")
    make_trees(afrolangs, afroages, "distance", "trees/afrodistance")
    make_trees(austrolangs, austroages, "feature", "trees/austrofeature")
    make_trees(austrolangs, austroages, "family", "trees/austrofamily")
    make_trees(austrolangs, austroages, "distance", "trees/austrodistance")
    make_trees(indolangs, indoages, "family", "trees/indofamily")
    make_trees(indolangs, indoages, "distance", "trees/indodistance")
    make_trees(indolangs, indoages, "feature", "trees/indofeature")
    make_trees(nigerlangs, nigerages, "family", "trees/nigerfamily")
    make_trees(nigerlangs, nigerages, "distance", "trees/nigerdistance")
    make_trees(nigerlangs, nigerages, "feature", "trees/nigerfeature")
    make_trees(nilolangs, niloages, "feature", "trees/nilofeature")
    make_trees(nilolangs, niloages, "family", "trees/nilofamily")
    make_trees(nilolangs, niloages, "distance", "trees/nilodistance")
    make_trees(sinolangs, sinoages, "family", "trees/sinofamily")
    make_trees(sinolangs, sinoages, "distance", "trees/sinodistance")
    make_trees(sinolangs, sinoages, "feature", "trees/sinofeature")
#        make_trees(languages, "family", "trees/family%d" % i, subsample=True)
#        make_trees(languages, "distance", "trees/distance%d" % i, subsample=True)
#        make_trees(languages, "feature", "trees/feature%d" % i, subsample=True)

if __name__ == "__main__":
    main()
