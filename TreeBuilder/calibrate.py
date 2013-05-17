#!/usr/bin/env python
import pdb
import getopt
import math
import sys
import sqlite3
import random

import dendropy

sys.path.append("../WALS2SQL/")
from wals2sql import create_dense_subset
from fileio import load_translate_file
from generate_trees import get_languages_by_family
import distance

def tree_from_line(line):
    line = line[line.index("("):]
    return dendropy.Tree.get_from_string(line, "newick")

def matrix_from_tree(tree):
    return dendropy.treecalc.PatristicDistanceMatrix(tree)

def load_indo_trees(filename, count):
    trees = []
    fp = open(filename)
    for i in range(0, count):
        line = fp.readline()
        if not line:
            break
        tree = tree_from_line(line)
        trees.append(tree)
    fp.close()
    return trees

def load_matrix(filename):
    fp = open(filename, "r")
    string = fp.read()
    fp.close()
    tree = dendropy.Tree.get_from_string(string, "newick")
    return dendropy.treecalc.PatristicDistanceMatrix(tree)

def average(x):
    assert len(x) > 0
    return float(sum(x)) / len(x)

def correlation(x, y):
    assert len(x) == len(y)
    n = len(x)
    assert n > 0
    avg_x = average(x)
    avg_y = average(y)
    diffprod = 0
    xdiff2 = 0
    ydiff2 = 0
    for idx in range(n):
        xdiff = x[idx] - avg_x
        ydiff = y[idx] - avg_y
        diffprod += xdiff * ydiff
        xdiff2 += xdiff * xdiff
        ydiff2 += ydiff * ydiff

    return diffprod / math.sqrt(xdiff2 * ydiff2)

def build_friendly_matrix_from_unfriendly_matrix(unfriendly, trans):
    friendly = {}
    name_pairs = []
    for lang1 in trans:
        friendly[lang1] = {}
        for lang2 in trans:
            if lang1 != lang2 and (lang2, lang1) not in name_pairs:
                name_pairs.append((lang1, lang2))
    for l1, l2 in name_pairs:
       friendly[l1][l2] = unfriendly[trans[l1]][trans[l2]]
       friendly[l2][l1] = friendly[l1][l2]
    return friendly


def build_friendly_matrix_from_tree(tree, trans):
    """
    Turn a tree and a translate dictionary into a "dictionary based matrix"
    showing distances between languages by name
    """
    matrix = dendropy.treecalc.PatristicDistanceMatrix(tree)
    friendly = {}
    name_pairs = []
    for lang1 in trans:
        friendly[lang1] = {}
        for lang2 in trans:
            if lang1 != lang2 and (lang2, lang1) not in name_pairs:
                name_pairs.append((lang1, lang2))

    for l1, l2 in name_pairs:
       t1 = tree.find_node_with_taxon_label(str(trans[l1])).taxon
       t2 = tree.find_node_with_taxon_label(str(trans[l2])).taxon
       friendly[l1][l2] = matrix(t1, t2)
       friendly[l2][l1] = friendly[l1][l2]
    return friendly

def get_correlation(matrix1, trans1, matrix2, trans2):
    common_langs = filter(lambda(x): x in trans1, trans2)
    pairs = []
    for lang1 in common_langs:
        for lang2 in common_langs:
            if lang1 != lang2 and (lang2, lang1) not in pairs:
                pairs.append((lang1, lang2))
    seq1 = [matrix1[t1][t2] for (t1, t2) in pairs]
    seq2 = [matrix2[t1][t2] for (t1, t2) in pairs]
    return correlation(seq1, seq2)

def compute_correlations(method_austro_matrix, wals_austro_trans, auth_austro_matrix, auth_austro_trans, method_indo_matrix, wals_indo_trans, auth_indo_matrices, auth_indo_trans):
    austro_correl = get_correlation(method_austro_matrix, wals_austro_trans, auth_austro_matrix, auth_austro_trans)
    indo_correls = [get_correlation(method_indo_matrix, wals_indo_trans, auth_indo_matrix, auth_indo_trans) for auth_indo_matrix in auth_indo_matrices]
    return austro_correl, indo_correls

def evaluate_method(long_name, matrix_builder, austrolangs, wals_austro_trans, auth_austro_matrix, auth_austro_trans, indolangs, wals_indo_trans, auth_indo_matrices, auth_indo_trans, vocal=True):
    method_austro_matrix = build_friendly_matrix_from_unfriendly_matrix(matrix_builder(austrolangs), wals_austro_trans)
    method_indo_matrix = build_friendly_matrix_from_unfriendly_matrix(matrix_builder(indolangs), wals_indo_trans)
    austro_correl, indo_correls = compute_correlations(method_austro_matrix, wals_austro_trans, auth_austro_matrix, auth_austro_trans, method_indo_matrix, wals_indo_trans, auth_indo_matrices, auth_indo_trans)
    if vocal:
        print long_name
        print austro_correl
        print sum(indo_correls)/len(indo_correls)
    return austro_correl, sum(indo_correls)/len(indo_correls)

def summarise(austrolangs, wals_austro_trans, auth_austro_matrix, auth_austro_trans, indolangs, wals_indo_trans, auth_indo_matrices, auth_indo_trans):

    evaluate_method("Linear geography", distance.linear_geography_matrix, austrolangs, wals_austro_trans, auth_austro_matrix, auth_austro_trans, indolangs, wals_indo_trans, auth_indo_matrices, auth_indo_trans)
    evaluate_method("Genetic", distance.genetic_matrix, austrolangs, wals_austro_trans, auth_austro_matrix, auth_austro_trans, indolangs, wals_indo_trans, auth_indo_matrices, auth_indo_trans)
    evaluate_method("Feature", distance.feature_matrix, austrolangs, wals_austro_trans, auth_austro_matrix, auth_austro_trans, indolangs, wals_indo_trans, auth_indo_matrices, auth_indo_trans)
    evaluate_method("Optimal combination", distance.optimal_combination_matrix, austrolangs, wals_austro_trans, auth_austro_matrix, auth_austro_trans, indolangs, wals_indo_trans, auth_indo_matrices, auth_indo_trans)

def optimise_n(wals_austro_trans, auth_austro_matrix, auth_austro_trans, wals_indo_trans, auth_indo_matrices, auth_indo_trans):

    conn = sqlite3.connect("../WALS2SQL/wals.db")
    cursor = conn.cursor()
    cursor.execute('''PRAGMA cache_size = -25000''')
    create_dense_subset(conn, cursor)
    cursor.execute('''CREATE TEMPORARY TABLE speedyfeatures AS SELECT name, value_id, wals_code FROM data_points INNER JOIN dense_features on data_points.feature_id = dense_features.id''')
    cursor.execute('''CREATE INDEX wals_code_index ON speedyfeatures(wals_code)''')

    bestn = 0
    bestc = 0
    for n in range(1,30):
        create_dense_subset(conn, cursor, n)
        cursor.execute('''DROP TABLE speedyfeatures''')
        cursor.execute('''CREATE TEMPORARY TABLE speedyfeatures AS SELECT name, value_id, wals_code FROM data_points INNER JOIN dense_features on data_points.feature_id = dense_features.id''')
        cursor.execute('''CREATE INDEX wals_code_index ON speedyfeatures(wals_code)''')
        austrolangs = get_languages_by_family(conn, cursor, "Austronesian")
        indolangs = get_languages_by_family(conn, cursor, "Indo-European")
        c1, c2 = evaluate_method("Feature", distance.feature_matrix, austrolangs, wals_austro_trans, auth_austro_matrix, auth_austro_trans, indolangs, wals_indo_trans, auth_indo_matrices, auth_indo_trans, vocal=False)
        meanc = 0.5*(c1+c2)
        if meanc > bestc:
            bestc = meanc
            bestn = n
            print "Best n so far is: ", bestn
            print "Best correlations: ", c1, c2

    cursor.close()
    conn.close()

def optimise_weights(austrolangs, wals_austro_trans, auth_austro_matrix, auth_austro_trans, indolangs, wals_indo_trans, auth_indo_matrices, auth_indo_trans):

    random.seed()
    bestweights = None
    bestc = 0
    weights = [random.random() for i in range(0,3)]
    func = distance.weighted_triple_factory(weights)
    c1, c2 = evaluate_method("Random weighted", func, austrolangs, wals_austro_trans, auth_austro_matrix, auth_austro_trans, indolangs, wals_indo_trans, auth_indo_matrices, auth_indo_trans, vocal=False)
    var = 0.1
    regress = 0.2
    for i in range(0,250):
        # Make a change
        if random.random() < 0.95:
            # Most of the time, just make a small change to one component
            newweights = weights[:]
            change_index = random.randint(0,2)
            newweights[change_index] += random.gauss(0, var)
            if newweights[change_index] < 0:
                newweights[change_index] = 0
        else:
            # Every now and then do a random restart
            newweights = [random.random() for i in range(0,3)]

        func = distance.weighted_triple_factory(newweights)
        newc1, newc2 = evaluate_method("Random weighted", func, austrolangs, wals_austro_trans, auth_austro_matrix, auth_austro_trans, indolangs, wals_indo_trans, auth_indo_matrices, auth_indo_trans)
        meanc = 0.5*(newc1+newc2)
        if meanc > 0.5*(c1+c2) or random.random() < regress:
            weights = newweights
            c1, c2 = newc1, newc2
        if meanc > bestc:
            bestc = meanc
            bestweights = weights
            print "New best: ", bestweights, " -=> ", c1, c2
        var *= 0.99
        regress *= 0.99
    print "Best weights: ", bestweights

def main():

    # Load WALS language data
    conn = sqlite3.connect("../WALS2SQL/wals.db")
    cursor = conn.cursor()
    cursor.execute('''PRAGMA cache_size = -25000''')
    create_dense_subset(conn, cursor)
    cursor.execute('''CREATE TEMPORARY TABLE speedyfeatures AS SELECT name, value_id, wals_code FROM data_points INNER JOIN dense_features on data_points.feature_id = dense_features.id''')
    cursor.execute('''CREATE INDEX wals_code_index ON speedyfeatures(wals_code)''')

    wals_austro_trans = load_translate_file("generated_trees/austro.translate")
    wals_indo_trans = load_translate_file("generated_trees/indo.translate")
    austrolangs = get_languages_by_family(conn, cursor, "Austronesian")
    indolangs = get_languages_by_family(conn, cursor, "Indo-European")

    cursor.close()
    conn.close()

    # Load authoritative trees
    auth_austro_tree = dendropy.Tree.get_from_path("authoritative_trees/austronesian/no_labels", "newick")
    auth_austro_trans = load_translate_file("authoritative_trees/austronesian/translate")
    common_austro_langs = filter(lambda(x): x in auth_austro_trans, wals_austro_trans)
    for lang in auth_austro_trans.keys():
        if lang not in common_austro_langs:
            auth_austro_trans.pop(lang)
    auth_austro_matrix = build_friendly_matrix_from_tree(auth_austro_tree, auth_austro_trans)

    auth_indo_trees = load_indo_trees("authoritative_trees/indo-european/no_labels", 10)
    auth_indo_trans = load_translate_file("authoritative_trees/indo-european/translate")
    common_indo_langs = filter(lambda(x): x in auth_indo_trans, wals_indo_trans)
    for lang in auth_indo_trans.keys():
        if lang not in common_indo_langs:
            auth_indo_trans.pop(lang)
    auth_indo_matrices = [build_friendly_matrix_from_tree(auth_indo_tree, auth_indo_trans) for auth_indo_tree in auth_indo_trees] 
 
    optimise_n(wals_austro_trans, auth_austro_matrix, auth_austro_trans, wals_indo_trans, auth_indo_matrices, auth_indo_trans)
    optimise_weights(austrolangs, wals_austro_trans, auth_austro_matrix, auth_austro_trans, indolangs, wals_indo_trans, auth_indo_matrices, auth_indo_trans)
    summarise(austrolangs, wals_austro_trans, auth_austro_matrix, auth_austro_trans, indolangs, wals_indo_trans, auth_indo_matrices, auth_indo_trans)
   

if __name__ == "__main__":
    main()
