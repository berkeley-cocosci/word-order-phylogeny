#!/usr/bin/env python
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
            random.shuffle(self.common_austro_langs)
            self.common_austro_langs = self.common_austro_langs[0:len(self.common_indo_langs)]
        else:
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
        # Indo-European (store mean of many vectors)
        auth_indo_vectors = np.array([self.compute_auth_vector(auth_indo_tree, self.common_indo_langs, self.auth_indo_trans) for auth_indo_tree in self.auth_indo_trees])
        self.auth_indo_vector = np.mean(auth_indo_vectors, axis=0)
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
        families and return a 3-tuple of the model, the Austronesian correlation
        and the Indo-European correlation
        """
        # Fit the model to the combined authoritative data
        df = {}
        df["auth"] = self.auth_combo_vector
        df["method"] = np.concatenate([method_austro_vector, method_indo_vector])
        DF = pd.DataFrame(df)
        combined_model = smf.ols('auth ~ method', data=DF).fit()
        if label:
            df.to_csv("calibration_results/%s_data.csv" % label)

        # Get predictions for Austronesian, and compare to authoritative
        df.pop("auth")
        df["method"] = method_austro_vector
        DF = pd.DataFrame(df)
        combo_austro_predictions = combined_model.predict(DF)
        austro_correl = np.corrcoef(combo_austro_predictions, self.auth_austro_vector)[0,1]

        # Get predictions for Indonesian, and compare to authoritative
        df["method"] = method_indo_vector
        DF = pd.DataFrame(df)
        combo_indo_predictions = combined_model.predict(DF)
        indo_correl = np.corrcoef(combo_indo_predictions, self.auth_indo_vector)[0,1]

        return combined_model, austro_correl, indo_correl

    def optimise_geographic(self):

        func = distance.geographic_matrix_factory()
        combined_model, austro_correl, indo_correl = self.evaluate_method(func, "geo")
        self.max_geo_combined = math.sqrt(combined_model.rsquared)
        self.max_geo_austro = austro_correl
        self.max_geo_indo = indo_correl
        with open("calibration_results/optimal_geographic_parameters", "w") as fp:
            fp.write("%f\n" % combined_model.params["Intercept"])
            fp.write("%f\n" % combined_model.params["method"])

    def optimise_genetic(self):
        """
        Find the optimal rate at which to discount the importance of increasingly
        fine-grained genetic category matches.
        """
        N = 100
        bestparam = 0
        bestc = 0
        for i in range(0, N):
            param = (i+1)*(1.0/N)
            func = distance.genetic_matrix_factory(param)
            combined_model, austro_correl, indo_correl = self.evaluate_method(func, "gen")
            if combined_model.rsquared > bestc:
                best_model = combined_model
                best_combined = math.sqrt(combined_model.rsquared)
                best_austro = austro_correl
                best_indo = indo_correl
                best_param = param
        self.max_gen_combined = best_combined
        self.max_gen_austro = best_austro
        self.max_gen_indo = best_indo

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
        df = pd.DataFrame(df)
        df.to_csv("calibration_results/feature_data.csv")
        model_spec = "auth ~ " + " + ".join(good_features)
        model = smf.ols(model_spec, data=df).fit()

        fp = open("calibration_results/optimal_feature_weights", "w")
        for index, feature in enumerate(dense_features):
            fp.write("%f\tintercept\n" % model.params["Intercept"])
            if "feat%d" % index in good_features:
                fp.write("%f\t%s\n" % (model.params["feat%d" % index], feature))
        fp.close()

        self.max_feat_combined = model.rsquared

        olddf = df
        df = {}
        for key in olddf:
            df[key] = olddf[key][0:len(self.auth_austro_vector)]
        df = pd.DataFrame(df)
        model_spec = "auth ~ " + " + ".join(good_features)
        model = smf.ols(model_spec, data=df).fit()
        self.max_feat_austro = model.rsquared

        df = {}
        for key in olddf:
            df[key] = olddf[key][len(self.auth_austro_vector):]
            assert len(df[key]) == len(self.auth_indo_vector)
        df = pd.DataFrame(df)
        model_spec = "auth ~ " + " + ".join(good_features)
        model = smf.ols(model_spec, data=df).fit()
        self.max_feat_indo = model.rsquared

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

        self.max_combo_combined = model.rsquared

        df = {}
        df["auth"] = self.auth_austro_vector
        df["geo"] = geo_austro
        df["gen"] = gen_austro
        df["feat"] = feat_austro
        df = pd.DataFrame(df)
        model = smf.ols('auth ~ geo + gen + feat', data=df).fit()
        self.max_combo_austro = model.rsquared

        df = {}
        df["auth"] = self.auth_indo_vector
        df["geo"] = geo_indo
        df["gen"] = gen_indo
        df["feat"] = feat_indo
        df = pd.DataFrame(df)
        model = smf.ols('auth ~ geo + gen + feat', data=df).fit()
        self.max_combo_indo = model.rsquared

    def summarise(self):

        pt = prettytable.PrettyTable(["Method", "Combined correl", "Austronesian correl", "Indo-European correl"])
        pt.add_row(["GEOGRAPHIC", self.max_geo_combined, self.max_geo_austro, self.max_geo_indo])
        pt.add_row(["GENETIC", self.max_gen_combined, self.max_gen_austro, self.max_gen_indo])
        pt.add_row(["FEATURE", self.max_feat_combined, self.max_feat_austro, self.max_feat_indo])
        pt.add_row(["COMBINATION", self.max_combo_combined, self.max_combo_austro, self.max_combo_indo])
        print pt

if __name__ == "__main__":
    calibrator = Calibrator()
    calibrator.optimise_geographic()
    calibrator.optimise_genetic()
    calibrator.optimise_feature()
    calibrator.optimise_combination()
    calibrator.summarise()
