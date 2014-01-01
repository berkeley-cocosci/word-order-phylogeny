#!/usr/bin/env python
import gc
import pdb
import getopt
import itertools
import math
import numpy as np
import os
import pdb
import pandas as pd
import statsmodels.formula.api as smf
import sys
import sqlite3
import random

import dendropy
import prettytable

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
        self.austrolangs = wals2sql.get_languages_by_family(conn, cursor, "Austronesian")
        self.indolangs = wals2sql.get_languages_by_family(conn, cursor, "Indo-European")

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

        self.auth_indo_trees = self.load_indo_trees("authoritative_trees/indo-european/no_labels", 500)
        self.auth_indo_trans = load_translate_file("authoritative_trees/indo-european/wals_compatible_translate")

    def find_common_langs(self):
        self.common_austro_langs = filter(lambda(x): x in self.auth_austro_trans, self.wals_austro_trans)
        self.common_indo_langs = filter(lambda(x): x in self.auth_indo_trans, self.wals_indo_trans)
        if len(self.common_austro_langs) > len(self.common_indo_langs):
            self.all_common_austro_langs = self.common_austro_langs[:]
            self.all_common_indo_langs = None
            random.shuffle(self.common_austro_langs)
            self.common_austro_langs = self.common_austro_langs[0:len(self.common_indo_langs)]
        else:
            self.all_common_indo_langs = self.common_indo_langs[:]
            self.all_common_austro_langs = None
            random.shuffle(self.common_indo_langs)
            self.common_indo_langs = self.common_indo_langs[0:len(self.common_austro_langs)]
        assert len(self.common_austro_langs) == len(self.common_indo_langs)

    def compute_auth_vector(self, tree, common_langs, trans):
        matrix = dendropy.treecalc.PatristicDistanceMatrix(tree)
        auth_vector = []
        for l1, l2 in itertools.combinations(common_langs, 2):
            t1 = tree.find_node_with_taxon_label(str(trans[l1])).taxon
            t2 = tree.find_node_with_taxon_label(str(trans[l2])).taxon
            auth_vector.append(matrix(t1, t2))
        return np.array(auth_vector)

    def compute_auth_vectors(self):
        # Austronesian (simple)
        self.auth_austro_vector = self.compute_auth_vector(self.auth_austro_tree, self.common_austro_langs, self.auth_austro_trans)
        if self.all_common_austro_langs:
            self.auth_all_austro_vector = self.compute_auth_vector(self.auth_austro_tree, self.all_common_austro_langs, self.auth_austro_trans)
        # Indo-European (store mean of many vectors)
        auth_indo_vectors = np.array([self.compute_auth_vector(auth_indo_tree, self.common_indo_langs, self.auth_indo_trans) for auth_indo_tree in self.auth_indo_trees])
        self.auth_indo_vector = np.mean(auth_indo_vectors, axis=0)
        if self.all_common_indo_langs:
            auth_all_indo_vectors = np.array([self.compute_auth_vector(auth_indo_tree, self.all_common_indo_langs, self.auth_indo_trans) for auth_indo_tree in self.auth_indo_trees])
            self.auth_all_indo_vector = np.mean(auth_all_indo_vectors, axis=0)
        # Normalise
        self.auth_indo_vector = self.auth_indo_vector / self.auth_indo_vector.max()
        self.auth_austro_vector = self.auth_austro_vector / self.auth_austro_vector.max()
        # Combination
        self.auth_combo_vector = np.concatenate([self.auth_austro_vector, self.auth_indo_vector])

    def compute_method_vector(self, matrix, common_langs, trans):
        """
        Transform a matrix of pairwise distances into a vector of
        distances between those languages which are in our WALS data
        AND on the authoritative trees.
        """
        method_vector = [matrix[trans[l1]][trans[l2]] for l1, l2 in itertools.combinations(common_langs, 2)]
        return np.array(method_vector)

    def evaluate_method(self, matrix_builder, label=None):
        """
        Build matrices using matrix_builder function and compare the
        pairwise distances between vectors to those taken from the
        authoritative trees, computing the correlation between the two.
        """
        method_austro_vector = self.compute_method_vector(matrix_builder(self.austrolangs), self.common_austro_langs, self.wals_austro_trans)
        method_indo_vector = self.compute_method_vector(matrix_builder(self.indolangs), self.common_indo_langs, self.wals_indo_trans)
        return self.fit_models(method_austro_vector, method_indo_vector, label)

    def fit_models(self, method_austro_vector, method_indo_vector, label=None):
        """
        Fit a model to the authoritative data for both reference language
        families and return the model.
        """
        # Fit the model to the combined authoritative data
        df = {}
        df["auth"] = self.auth_combo_vector
        df["method"] = np.concatenate([method_austro_vector, method_indo_vector])
        df = pd.DataFrame(df)
        model = smf.ols('auth ~ method', data=df).fit()
        if label:
            df.to_csv("calibration_results/%s_data.csv" % label)

        return model

    def optimise_geographic(self):

        func = distance.geographic_matrix_factory()
        model = self.evaluate_method(func, "geo")
        with open("calibration_results/optimal_geographic_parameters", "w") as fp:
            fp.write("%f\n" % model.params["Intercept"])
            fp.write("%f\n" % model.params["method"])

    def optimise_genetic(self):
        """
        Find the optimal rate at which to discount the importance of increasingly
        fine-grained genetic category matches.
        """
        N = 100
        best_rsquared = 0
        for i in range(0, N):
            param = (i+1)*(1.0/N)
            func = distance.genetic_matrix_factory(param)
            model = self.evaluate_method(func, "gen")
            if model.rsquared > best_rsquared:
                best_model = model
                best_param = param

        with open("calibration_results/optimal_genetic_parameters", "w") as fp:
            fp.write("%f\n" % best_model.params["Intercept"])
            fp.write("%f\n" % best_model.params["method"])
            fp.write("%f\n" % best_param)

    def optimise_feature(self):

        conn = sqlite3.connect("../WALS2SQL/wals.db")
        cursor = conn.cursor()
        cursor.execute('''PRAGMA cache_size = -25000''')

        wals2sql.compute_dense_features(conn, cursor, 25)
        dense_features = wals2sql.get_dense_features(conn, cursor)
        cursor.close()
        conn.close()
        comparators = distance.build_comparators()

        # Ugly hack
        langs_by_name = {}
        for lang in self.austrolangs:
            langs_by_name[lang.name] = lang
        for lang in self.indolangs:
            langs_by_name[lang.name] = lang

        df = {}
        df["auth"] = self.auth_combo_vector
        good_features = []
        for index, feature in enumerate(dense_features):
            if feature == bwo:
                continue
            df["feat%d" % index] = []
            for l1, l2 in itertools.chain(itertools.combinations(self.common_austro_langs, 2), itertools.combinations(self.common_indo_langs, 2)):
                l1 = langs_by_name[l1]
                l2 = langs_by_name[l2]
                useful_points = 0
                if feature in l1.data and feature in l2.data:
                    df["feat%d" % index].append(comparators[feature](l1.data[feature], l2.data[feature]))
                    useful_points += 1
                else:
                    df["feat%d" % index].append(0.5)
            if useful_points > 0:
                good_features.append("feat%d" % index)
            else:
                # There's literally no pairwise data for this feature!
                df.pop("feat%d" % index)
        df["sum"] = []
        for i in range(0, len(df["auth"])):
            df["sum"].append(math.log(sum([df[x][i] for x in good_features])+0.0001))
        df = pd.DataFrame(df)
        df.to_csv("calibration_results/feature_data.csv")
        model_spec = "auth ~ " + " + ".join(good_features)
        model = smf.ols(model_spec, data=df).fit()

        #print "Linear: ", model.rsquared
        #df["fitted"] = [math.log(x) for x in list(model.fittedvalues)]
        #model2 = smf.ols("auth ~ fitted", data=df).fit()
        #print "Log: ", model2.rsquared
        #df["fitted"] = [math.exp(x) for x in list(model.fittedvalues)]
        #model3 = smf.ols("auth ~ fitted", data=df).fit()
        #print "Exp: ", model3.rsquared

        fp = open("calibration_results/optimal_feature_weights", "w")
        fp.write("%f\tintercept\n" % model.params["Intercept"])
        for index, feature in enumerate(dense_features):
            if "feat%d" % index in good_features:
                fp.write("%f\t%s\n" % (model.params["feat%d" % index], feature))
        fp.close()

    def optimise_combination(self):
        """
        Use multiple linear regression to determine the optimal weighted
        combination of the GEOGRAPHIC, GENETIC and FEATUE methods.
        """

        geo_austro = self.compute_method_vector(distance.build_optimal_geographic_matrix(self.austrolangs), self.common_austro_langs, self.wals_austro_trans)
        geo_indo =  self.compute_method_vector(distance.build_optimal_geographic_matrix(self.indolangs), self.common_indo_langs, self.wals_indo_trans)
        gen_austro = self.compute_method_vector(distance.build_optimal_genetic_matrix(self.austrolangs), self.common_austro_langs, self.wals_austro_trans)
        gen_indo =  self.compute_method_vector(distance.build_optimal_genetic_matrix(self.indolangs), self.common_indo_langs, self.wals_indo_trans)
        feat_austro = self.compute_method_vector(distance.build_optimal_feature_matrix(self.austrolangs), self.common_austro_langs, self.wals_austro_trans)
        feat_indo =  self.compute_method_vector(distance.build_optimal_feature_matrix(self.indolangs), self.common_indo_langs, self.wals_indo_trans)

        df = {}
        df["auth"] = self.auth_combo_vector
        df["geo"] = np.concatenate([geo_austro, geo_indo])
        df["gen"] = np.concatenate([gen_austro, gen_indo])
        df["feat"] = np.concatenate([feat_austro, feat_indo])
        df = pd.DataFrame(df)
        df.to_csv("calibration_results/combination_data.csv")
        model = smf.ols('auth ~ geo + gen + feat', data=df).fit()
        weights = [model.params[x] for x in ("geo", "gen", "feat")]

        fp = open("calibration_results/optimal_combination_weights", "w")
        fp.write("intercept\t%f\n" % model.params["Intercept"])
        fp.write("geo\t%f\n" % weights[0])
        fp.write("gen\t%f\n" % weights[1])
        fp.write("feat\t%f\n" % weights[2])
        fp.close()

        return (model.params["Intercept"], weights[0], weights[1], weights[2])

    def summarise(self):

        """
        Compute the correlation between each of the four methods and the
        two authoritative trees, displaying the results in a PrettyTable.
        """
        pt = prettytable.PrettyTable(["Method", "Austronesian correl", "Indo-European correl"])

        names = ("GEOGRAPHIC", "GENETIC", "FEATURE", "COMBINATION")
        funcs = (distance.build_optimal_geographic_matrix,
                    distance.build_optimal_genetic_matrix,
                    distance.build_optimal_feature_matrix,
                    distance.build_optimal_combination_matrix)
        for name, func in zip(names, funcs):
            if self.all_common_austro_langs:
                austro_data = self.compute_method_vector(func(self.austrolangs), self.all_common_austro_langs, self.wals_austro_trans)
                austro_correl = np.corrcoef(austro_data, self.auth_all_austro_vector)[0,1]
            else:
                austro_data = self.compute_method_vector(func(self.austrolangs), self.common_austro_langs, self.wals_austro_trans)
                austro_correl = np.corrcoef(austro_data, self.auth_austro_vector)[0,1]
            if self.all_common_indo_langs:
                indo_data =  self.compute_method_vector(func(self.indolangs), self.all_common_indo_langs, self.wals_indo_trans)
                indo_correl = np.corrcoef(indo_data, self.auth_all_indo_vector)[0,1]
            else:
                indo_data =  self.compute_method_vector(func(self.indolangs), self.common_indo_langs, self.wals_indo_trans)
                indo_correl = np.corrcoef(indo_data, self.auth_indo_vector)[0,1]
            pt.add_row([name, austro_correl, indo_correl])

        print pt

def main(samples=10):

    comb_params = []
    for i in range(0,samples):
        gc.collect()
        calibrator = Calibrator()
        calibrator.optimise_geographic()
        calibrator.optimise_genetic()
        calibrator.optimise_feature()
        comb_params.append(calibrator.optimise_combination())
        if i != (samples-1):
            del(calibrator)
    mean_params = [sum([params[index] for params in comb_params])/len(comb_params) for index in range(0,4)]
    fp = open("calibration_results/optimal_combination_weights", "w")
    fp.write("intercept\t%f\n" % mean_params[0])
    fp.write("geo\t%f\n" % mean_params[1])
    fp.write("gen\t%f\n" % mean_params[2])
    fp.write("feat\t%f\n" % mean_params[3])
    fp.close()

    calibrator.summarise()

if __name__ == "__main__":
    main()
