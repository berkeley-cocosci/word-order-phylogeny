#!/usr/bin/env python
import pdb
import getopt
import itertools
import math
import numpy as np
import os
import pdb
import scipy.stats
import pandas as pd
import statsmodels.formula.api as smf
import sys
import sqlite3
import random

import dendropy

sys.path.append("../WALS2SQL/")
bwo = "Order of Subject, Object and Verb"
import wals2sql
from fileio import load_translate_file
import distance

class Calibrator:

    """
    The Calibrator's job is to find optimal values of various parameters
    relating to the different methods of computing between-language
    distances.  This is done by calibrating against some "authoritative"
    trees for the Austronesian and Indo-European families."
    """

    def __init__(self):
        self.load_wals_data()
        self.load_auth_trees()
        self.find_common_langs()
        self.compute_auth_vectors()

    def load_wals_data(self):

        conn = sqlite3.connect("../WALS2SQL/wals.db")
        cursor = conn.cursor()
        cursor.execute('''PRAGMA cache_size = -25000''')
        cursor.execute('''CREATE TEMPORARY TABLE speedyfeatures AS SELECT name, value_id, wals_code FROM data_points INNER JOIN dense_features on data_points.feature_id = dense_features.id''')
        cursor.execute('''CREATE INDEX wals_code_index ON speedyfeatures(wals_code)''')

        self.wals_austro_trans = load_translate_file("generated_trees/austro.translate")
        self.wals_indo_trans = load_translate_file("generated_trees/indo.translate")
        self.austrolangs = wals2sql.get_dense_languages_by_family(conn, cursor, "Austronesian")
        self.indolangs = wals2sql.get_dense_languages_by_family(conn, cursor, "Indo-European")

        cursor.close()
        conn.close()

    def load_indo_trees(self, filename, count):
        trees = []
        fp = open(filename)
        for i in range(0, count):
            line = fp.readline()
            if not line:
                break
            line = line[line.index("("):]
            tree = dendropy.Tree.get_from_string(line, "newick")
            trees.append(tree)
        fp.close()
        return trees

    def load_auth_trees(self):

        self.auth_austro_tree = dendropy.Tree.get_from_path("authoritative_trees/austronesian/no_labels", "newick")
        self.auth_austro_trans = load_translate_file("authoritative_trees/austronesian/wals_compatible_translate")
        self.common_austro_langs = filter(lambda(x): x in self.auth_austro_trans, self.wals_austro_trans)

        self.auth_indo_trees = self.load_indo_trees("authoritative_trees/indo-european/no_labels", 10)
        self.auth_indo_trans = load_translate_file("authoritative_trees/indo-european/wals_compatible_translate")
        self.common_indo_langs = filter(lambda(x): x in self.auth_indo_trans, self.wals_indo_trans)

    def find_common_langs(self):
        self.common_austro_langs = filter(lambda(x): x in self.auth_austro_trans, self.wals_austro_trans)
        self.common_indo_langs = filter(lambda(x): x in self.auth_indo_trans, self.wals_indo_trans)

    def compute_auth_vector(self, tree, common_langs, trans):
        matrix = dendropy.treecalc.PatristicDistanceMatrix(tree)
        auth_vector = []
        for l1, l2 in itertools.combinations(common_langs, 2):
            t1 = tree.find_node_with_taxon_label(str(trans[l1])).taxon
            t2 = tree.find_node_with_taxon_label(str(trans[l2])).taxon
            auth_vector.append(matrix(t1, t2))
        return np.array(auth_vector)

    def compute_auth_vectors(self):
        self.auth_austro_vector = self.compute_auth_vector(self.auth_austro_tree, self.common_austro_langs, self.auth_austro_trans)
        self.auth_indo_vectors = [self.compute_auth_vector(auth_indo_tree, self.common_indo_langs, self.auth_indo_trans) for auth_indo_tree in self.auth_indo_trees]

    def compute_method_vector(self, matrix, common_langs, trans):
        """
        Transform a matrix of pairwise distances into a vector of
        distances between those languages which are in our WALS data
        AND on the authoritative trees.
        """
        method_vector = [matrix[trans[l1]][trans[l2]] for l1, l2 in itertools.combinations(common_langs, 2)]
        return np.array(method_vector)

    def evaluate_method(self, long_name, matrix_builder):
        """
        Build matrices using matrix_builder function and compare the
        pairwise distances between vectors to those taken from the
        authoritative trees, computing the correlation between the two.
        """
        method_austro_vector = self.compute_method_vector(matrix_builder(self.austrolangs), self.common_austro_langs, self.wals_austro_trans)
        method_indo_vector = self.compute_method_vector(matrix_builder(self.indolangs), self.common_indo_langs, self.wals_indo_trans)
        austro_correl = scipy.stats.pearsonr(self.auth_austro_vector, method_austro_vector)[0]
        indo_correls = [scipy.stats.pearsonr(auth_indo_vector, method_indo_vector)[0] for auth_indo_vector in self.auth_indo_vectors]
        return austro_correl, sum(indo_correls)/len(indo_correls)

    def optimise_feature_count(self):
        """
        Find the best value of n, such that data about the n densest
        features yields the best correlation with the authoritative
        trees.
        """
        conn = sqlite3.connect("../WALS2SQL/wals.db")
        cursor = conn.cursor()
        cursor.execute('''PRAGMA cache_size = -25000''')

        bestn = 0
        bestc = 0
        for n in range(1,30):
            # Recreate all the WALS stuff for this value of N.
            # This potential changes the set of languages that are common between
            # WALS and our authoritative trees, so we need to recompute those
            wals2sql.create_dense_subset(conn, cursor, n)
            self.load_wals_data()
            self.find_common_langs()
            self.compute_auth_vectors()

            c1, c2 = self.evaluate_method("Feature", distance.feature_matrix)
            meanc = 0.5*(c1+c2)
            if meanc > bestc:
                bestc = meanc
                bestn = n
        cursor.close()
        conn.close()

        self.optimal_feature_count = bestn

    def optimise_feature_weights(self):
        """
        Use random search to find the optimal weights for the various
        language features to define a weighted Hamming distance.
        """
        random.seed()
        bestweights = None
        bestc = 0
        weights = {}

        for feature in self.austrolangs[0].data:
            if feature not in ["genus", "subfamily", "family", "location", "Order of Subject, Object and Verb", "Order of Subject and Verb", "Order of Object and Verb", "iso_codes", "ethnoclass"] and not feature.startswith("Relationship between the Order of Object and Verb"):
                weights[feature] = random.random()

        # Do a rough simulated annealing kind of thing
        func = distance.weighted_feature_factory(weights)
        c1, c2 = self.evaluate_method("Random weighted", func)
        var = 0.01
        regress = 0.05
        for i in range(0,25):
            newweights = weights.copy()
            changekey = random.sample(newweights.keys(),1)[0]
            if random.random() < 0.9:
                newweights[changekey] = max(0, newweights[changekey] + random.gauss(0, var))
            else:
                if newweights[changekey]:
                    newweights[changekey] = 0.0
                else:
                    newweights[changekey] = random.random()

            func = distance.weighted_feature_factory(weights)
            newc1, newc2 = self.evaluate_method("Random weighted", func)
            meanc = 0.5*(newc1+newc2)
            if meanc > 0.5*(c1+c2) or random.random() < regress:
                weights = newweights
                c1, c2 = newc1, newc2
            if meanc > bestc:
                bestc = meanc
                bestweights = weights
            var *= 0.99
            regress *= 0.99

        fp = open("calibration/optimal_feature_weights", "w")
        for feature in bestweights:
            fp.write("%f\t%s\n" % (bestweights[feature], feature))
        fp.close()

        self.optimal_feature_weights = bestweights
        self.maximum_feature_correlation = bestc

    def optimise_combination(self):
        """
        Use multiple linear regression to determine the optimal weighted
        combination of the GEOGRAPHIC, GENETIC and FEATUE methods.
        """
        # Austro
        geographic_austro_vector = self.compute_method_vector(distance.linear_geography_matrix(self.austrolangs), self.common_austro_langs, self.wals_austro_trans)
        genetic_austro_vector = self.compute_method_vector(distance.genetic_matrix(self.austrolangs), self.common_austro_langs, self.wals_austro_trans)
        feature_austro_vector = self.compute_method_vector(distance.optimal_feature_matrix(self.austrolangs), self.common_austro_langs, self.wals_austro_trans)

        df = {}
        df["auth"] = self.auth_austro_vector
        df["geo"] = geographic_austro_vector
        df["gen"] = genetic_austro_vector
        df["feat"] = feature_austro_vector
        df = pd.DataFrame(df)
        model = smf.ols('auth ~ geo + gen + feat', data=df).fit()
        weights1 = [model.params[x] for x in ("geo", "gen", "feat")]
        weights1 = [x/sum(weights1) for x in weights1]
        r21 = model.rsquared
        # Indo
        geographic_indo_vector = self.compute_method_vector(distance.linear_geography_matrix(self.indolangs), self.common_indo_langs, self.wals_indo_trans)
        genetic_indo_vector = self.compute_method_vector(distance.genetic_matrix(self.indolangs), self.common_indo_langs, self.wals_indo_trans)
        feature_indo_vector = self.compute_method_vector(distance.optimal_feature_matrix(self.indolangs), self.common_indo_langs, self.wals_indo_trans)
        df = {}
        df["auth"] = self.auth_indo_vectors[0]
        df["geo"] = geographic_indo_vector
        df["gen"] = genetic_indo_vector
        df["feat"] = feature_indo_vector
        df = pd.DataFrame(df)
        model = smf.ols('auth ~ geo + gen + feat', data=df).fit()
        weights2 = [model.params[x] for x in ("geo", "gen", "feat")]
        weights2 = [x/sum(weights2) for x in weights2]
        r22 = model.rsquared

        # TODO Come up with a more principled way of combining weights
        weights = [(x+y)/2.0 for x, y in zip(weights1, weights2)]

        self.optimal_combination_weights = weights
        self.maximum_combination_correlation = 0.5*(r21 + r22)

    def optimise_genetic(self):
        """
        Find the optimal rate at which to discount the importance of increasingly
        fine-grained genetic category matches.
        """
        N = 50
        bestparam = 0
        bestc = 0
        for i in range(0, N+1):
            param = 0.70 + i*(0.05/N)
            func = distance.param_genetic_matrix_factory(param)
            c1, c2 = self.evaluate_method("Random weighted", func)
            mean = 0.5*(c1+c2)
            if mean > bestc:
                bestc = mean
                bestparam = param

        self.optimal_genetic_param = bestparam
        self.maximum_genetic_correlation = bestc

    def summarise(self):

        print "SUMMARY:"
        print "-"*79
        print "Best feature count: %d" % self.optimal_feature_count
        print "Best GENETIC correlation: %f" % self.maximum_genetic_correlation
        print "Best FEATURE correlation: %f" % self.maximum_feature_correlation
        print "Best COMBINATION correlation: %f" % self.maximum_combination_correlation

if __name__ == "__main__":
    calibrator = Calibrator()
    calibrator.optimise_feature_count()
    calibrator.optimise_feature_weights()
    calibrator.optimise_genetic()
    calibrator.optimise_combination()
    calibrator.summarise()
