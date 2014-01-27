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
import statsmodels.api as sm
import statsmodels.formula.api as smf
import scipy
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
        self.compute_auth_vectors()
        self.find_common_langs()
        self.compute_common_auth_vectors()

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

    def compute_auth_vector(self, tree):
        matrix = dendropy.treecalc.PatristicDistanceMatrix(tree)
        auth_vector = []
        for t1, t2 in itertools.combinations(tree.taxon_set, 2):
            auth_vector.append(matrix(t1, t2))
        auth_vector = np.array(auth_vector)
        auth_vector /= max(auth_vector)
        return auth_vector

    def compute_auth_vectors(self):
        # Austro
        self.auth_austro_vector = self.compute_auth_vector(self.auth_austro_tree)
        self.auth_austro_vector /= self.auth_austro_vector.max()
        print "Extrema for ALL of Au is: ", min(self.auth_austro_vector), max(self.auth_austro_vector)
        fp = open("calibration_results/auth_austro_pairwise.csv", "w")
        for dist in self.auth_austro_vector:
            fp.write("%f\n" % dist)
        fp.close()

        # Indo
        auth_indo_vectors = np.array([self.compute_auth_vector(auth_indo_tree) for auth_indo_tree in self.auth_indo_trees])
        print [min(v) for v in auth_indo_vectors][0:10]
        print [max(v) for v in auth_indo_vectors][0:10]
        auth_indo_vectors = [vec/vec.max() for vec in auth_indo_vectors]
        print [min(v) for v in auth_indo_vectors][0:10]
        print [max(v) for v in auth_indo_vectors][0:10]
        self.auth_indo_vector = np.mean(auth_indo_vectors, axis=0)
        print min(self.auth_indo_vector)
        print max(self.auth_indo_vector)
        self.auth_indo_vector /= self.auth_indo_vector.max()
        print "Extrema for ALL of IE is: ", min(self.auth_indo_vector), max(self.auth_indo_vector)
        fp = open("calibration_results/auth_indo_pairwise.csv", "w")
        for dist in self.auth_indo_vector:
            fp.write("%f\n" % dist)
        fp.close()

        self.auth_combo_vector = np.concatenate([self.auth_austro_vector, self.auth_indo_vector])

        long_indo_vector = list(self.auth_indo_vector[:])
        while len(long_indo_vector) < len(self.auth_austro_vector):
            long_indo_vector.append(random.sample(self.auth_indo_vector,1)[0])
        self.auth_fair_cdf = np.concatenate([self.auth_austro_vector, long_indo_vector])
        self.auth_fair_cdf.sort()

    def compute_common_auth_vector(self, tree, common_langs, trans):
        matrix = dendropy.treecalc.PatristicDistanceMatrix(tree)
        distances = matrix.distances()
        auth_vector = []
        for l1, l2 in itertools.combinations(common_langs, 2):
            t1 = tree.find_node_with_taxon_label(str(trans[l1])).taxon
            t2 = tree.find_node_with_taxon_label(str(trans[l2])).taxon
            auth_vector.append(matrix(t1, t2))
        auth_vector = np.array(auth_vector)
        auth_vector /= max(distances)
        return auth_vector

    def compute_common_auth_vectors(self):
        self.common_auth_austro_vector = self.compute_common_auth_vector(self.auth_austro_tree, self.common_austro_langs, self.auth_austro_trans)
        print "Extrema for common Au is: ", min(self.common_auth_austro_vector), max(self.common_auth_austro_vector)
        common_auth_indo_vectors = np.array([self.compute_common_auth_vector(tree, self.common_indo_langs, self.auth_indo_trans) for tree in self.auth_indo_trees])
        self.common_auth_indo_vector = np.mean(common_auth_indo_vectors, axis=0)
        print "Extrema for common IE is: ", min(self.common_auth_indo_vector), max(self.common_auth_indo_vector)
        
        self.common_auth_combo_vector = np.concatenate([self.common_auth_austro_vector, self.common_auth_indo_vector])

        l1 = len(self.common_auth_austro_vector)
        l2 = len(self.common_auth_indo_vector)
        t = l1 + l2
        self.weights = []
        self.weights.extend([1.0*t/l1 for i in range(0,l1)])
        self.weights.extend([1.0*t/l2 for i in range(0,l2)])
        self.weights = [w/sum(self.weights) for w in self.weights]
        print l1, l2, self.weights[0:10], self.weights[-10:]

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
        austro_matrix = matrix_builder(self.austrolangs)
        #method_austro_vector = self.compute_method_vector(austro_matrix, self.common_austro_langs, self.wals_austro_trans)
        indo_matrix = matrix_builder(self.indolangs)
        #method_indo_vector = self.compute_method_vector(indo_matrix, self.common_indo_langs, self.wals_indo_trans)
        return self.fit_models(austro_matrix, indo_matrix, label)

    def compute_cdf_distance(self, scaled):
        scaled.sort()
        scaled_points = [scaled[int(math.floor(i*(len(scaled)-1)/10.0))] for i in range(0,11)]
        auth_points = [self.auth_fair_cdf[int(math.floor(i*(len(self.auth_fair_cdf)-1)/10.0))] for i in range(0,11)]
        return sum([(scaled_points[i] - auth_points[i])**2 for i in range(0,11)])

    def fit_models(self, method_austro_matrix, method_indo_matrix, label=None):
        """
        Fit a model to the authoritative data for both reference language
        families and return the model.
        """
        lowest_D = 100
        best_int = 0
        best_mult = 0
        austro = list(itertools.chain(*method_austro_matrix))
        indo = list(itertools.chain(*method_indo_matrix))
        indo.extend(random.sample(austro,len(austro)-len(indo)))
            
        observations = [n for n in itertools.chain(austro, indo)]
        if min(observations) < 0:
            observations = [abs(min(observations)) + n for n in observations]
        default_span = max(observations) - min(observations)*1.0
        mindist = min(observations)
        for span in [0.01*n for n in range(1,101)]:
            scaled = [(1.0 - span) + (span/default_span)*(n-mindist) for n in observations]
            statistic = self.compute_cdf_distance(scaled)
            if statistic < lowest_D:
                lowest_D = statistic
                best_int = 1.0 - span + (span/default_span)*mindist
                best_mult = span/default_span
        print "THE BEST"
        print lowest_D, best_int, best_mult
        observations = [best_int + best_mult*n for n in observations]
        observations.sort()
#
        print [observations[i] for i in (0, len(observations)/4, len(observations)/2, 3*len(observations)/4, len(observations)-1)]
        return lowest_D, best_int, best_mult

    def optimise_geographic(self):

        func = distance.geographic_matrix_factory()
        D, intercept, mult = self.evaluate_method(func, "geo")
        with open("calibration_results/optimal_geographic_parameters", "w") as fp:
            fp.write("%f\n" % intercept)
            fp.write("%f\n" % mult)

        print "Best geographic D: ", D

    def optimise_genetic(self):
        """
        Find the optimal rate at which to discount the importance of increasingly
        fine-grained genetic category matches.
        """
        N = 100
        best_c = 0

        correl_ranks = []
        cdf_ranks = []
        df = {}
        df["auth"] = np.concatenate([self.common_auth_austro_vector, self.common_auth_indo_vector])
        for i in range(0, 100):
            param = (i+1)*(1.0/N)
            func = distance.genetic_matrix_factory(param)
            method_austro_vector = self.compute_method_vector(func(self.austrolangs), self.common_austro_langs, self.wals_austro_trans)
            method_indo_vector = self.compute_method_vector(func(self.indolangs), self.common_indo_langs, self.wals_indo_trans)
            df["method"] = np.concatenate([method_austro_vector, method_indo_vector])
            dff = pd.DataFrame(df)
            austrodf = dff[0:len(self.common_auth_austro_vector)]
            indodf = dff[len(self.common_auth_austro_vector):]
            austromodel = smf.ols('auth ~ method', data=austrodf).fit()
            indomodel = smf.ols('auth ~ method', data=indodf).fit()
            assert(math.sqrt(austromodel.rsquared) <= 1.0)
            assert(math.sqrt(indomodel.rsquared) <= 1.0)
            mean_correl = (math.sqrt(austromodel.rsquared)+math.sqrt(indomodel.rsquared))/2.0
            correl_ranks.append((mean_correl, param))
            cdf_ranks.append((self.compute_cdf_distance(df["method"]), param))
        correl_ranks.sort()
        correl_ranks.reverse()
        correl_ranks = [param for (cor, param) in correl_ranks]
        cdf_ranks.sort()
        cdf_ranks = [param for (cdf, param) in cdf_ranks]
        mean_ranks = [(correl_ranks.index(param) + cdf_ranks.index(param), param) for param in correl_ranks]
        mean_ranks.sort()
#        for (cor, param_cor), (cdf, param_cdf) in zip(correl_ranks, cdf_ranks):
#            mean_ranks[param_cor] = mean_ranks.get(param_cor,0) + cor
#            mean_ranks[param_cdf] = mean_ranks.get(param_cdf,0) + cdf
#        mean_ranks = [(mean_ranks[param], param) for param in mean_ranks.keys()]
#        mean_ranks.sort()
#        mean_ranks.reverse()
#        pdb.set_trace()
        bestparam = mean_ranks[0][1]
        print correl_ranks
        print cdf_ranks
        print mean_ranks
        print bestparam
#            if mean_correl > best_c:
#                best_c = mean_correl
#                bestparam = param

        param = bestparam
        func = distance.genetic_matrix_factory(param)
        D, intercept, mult = self.evaluate_method(func, "gen")

        with open("calibration_results/optimal_genetic_parameters", "w") as fp:
            fp.write("%f\n" % intercept)
            fp.write("%f\n" % mult)
            fp.write("%f\n" % param)

        print "Best genetic D: ", D

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

        # Identify good features
        good_features = []
        long_good_features = []
        for index, feature in enumerate(dense_features):
            if feature == bwo:
                continue
            for l1, l2 in itertools.chain(itertools.combinations(self.common_austro_langs, 2), itertools.combinations(self.common_indo_langs, 2)):
                l1 = langs_by_name[l1]
                l2 = langs_by_name[l2]
                useful_points = 0
                if feature in l1.data and feature in l2.data:
                    useful_points += 1
            if useful_points > 0:
                good_features.append("feat%d" % index)
                long_good_features.append(feature)

        # Compute supermeans
        austromeans = {}
        austrosupermean = 0
        austrosupernorm = 0
        for feature in long_good_features:
            austromeans[feature] = 0
            norm = 0
            for l1, l2 in itertools.combinations(self.common_austro_langs, 2):
                l1 = langs_by_name[l1]
                l2 = langs_by_name[l2]
#                pdb.set_trace()
                if feature in l1.data and feature in l2.data:
                    austromeans[feature] += comparators[feature](l1.data[feature], l2.data[feature])
                    norm += 1
            if norm:
                austromeans[feature] /= norm
                austrosupermean += austromeans[feature]
                austrosupernorm += 1
            else:
                austromeans[feature] = "NODATA"
        if austrosupernorm:
            austrosupermean /= austrosupernorm
        else:
            austrosupermean = 0.5
        for feature in austromeans:
            if austromeans[feature] == "NODATA":
                austromeans[feature] = austrosupermean

        indomeans = {}
        indosupermean = 0
        indosupernorm = 0
        for feature in long_good_features:
            indomeans[feature] = 0
            norm = 0
            for l1, l2 in itertools.combinations(self.common_indo_langs, 2):
                l1 = langs_by_name[l1]
                l2 = langs_by_name[l2]
                if feature in l1.data and feature in l2.data:
                    indomeans[feature] += comparators[feature](l1.data[feature], l2.data[feature])
                    norm += 1
            if norm:
                indomeans[feature] /= norm
                indosupermean += indomeans[feature]
                indosupernorm += 1
            else:
                indomeans[feature] = "NODATA"
        if indosupernorm:
            indosupermean /= indosupernorm
        else:
            indosupermean = 0.5
        for feature in indomeans:
            if indomeans[feature] == "NODATA":
                indomeans[feature] = indosupermean

        # Actually compute raw data
        df = {}
        df["auth"] = self.common_auth_combo_vector
        for feature, long_feature in zip(good_features, long_good_features):
            if long_feature == bwo:
                continue
            df[feature] = []
            for l1, l2 in itertools.chain(itertools.combinations(self.common_austro_langs, 2), itertools.combinations(self.common_indo_langs, 2)):
                l1 = langs_by_name[l1]
                l2 = langs_by_name[l2]
                if long_feature in l1.data and long_feature in l2.data:
                    df[feature].append(comparators[long_feature](l1.data[long_feature], l2.data[long_feature]))
                else:
                    if l1 in self.austrolangs:
                        df[feature].append(austromeans[long_feature])
                    else:
                        df[feature].append(indomeans[long_feature])

        print good_features
        df = pd.DataFrame(df)
        df.to_csv("calibration_results/feature_data.csv")

        austrodf = df[0:len(self.common_auth_austro_vector)]
        indodf = df[len(self.common_auth_austro_vector):]
        best_mean_correl = 0
        best_selectors = None
        rank = []
        for feature_selectors in itertools.product(*[(True, False),]*10):
            if feature_selectors.count(True) == 0:
                continue
            model_spec = "auth ~ " + " + ".join([feat for feat, sel in zip(good_features, feature_selectors) if sel])
            model = smf.wls(model_spec, data=df, weights=self.weights).fit()
            austrofit = model.fittedvalues[0:len(self.common_auth_austro_vector)]
            austroauth = austrodf["auth"]
            indofit = model.fittedvalues[len(self.common_auth_austro_vector):]
            indoauth = indodf["auth"]
#            austromodel = smf.ols(model_spec, data=austrodf).fit()
#            indomodel = smf.ols(model_spec, data=indodf).fit()
#            mean_correl = (math.sqrt(austromodel.rsquared)+math.sqrt(indomodel.rsquared))/2.0
            # Let's minimise the minimum instead!
            mean_correl = sorted((austroauth.corr(austrofit), indoauth.corr(indofit), 0.5*(austroauth.corr(austrofit)+indoauth.corr(indofit))))
            thingy = (mean_correl, feature_selectors.count(True), feature_selectors)
            rank.append(thingy)
            if mean_correl > best_mean_correl:
                best_mean_correl = mean_correl
                best_selectors = feature_selectors[:]

        rank.sort()
        rank.reverse()
        # Pay no attention to the magic number 6 behind the curtain...
        best_selectors = rank[6][2]
        best_features = [feat for feat, sel in zip(good_features, best_selectors) if sel]
        model_spec = "auth ~ " + " + ".join(best_features)
        model = smf.wls(model_spec, data=df, weights=self.weights).fit()
        weights = {}
        for index, feature in enumerate(dense_features):
            if "feat%d" % index in best_features:
                weights[feature] = model.params["feat%d" % index]
                print index, weights[feature]

        func = distance.feature_matrix_factory(weights)
        D, intercept, mult = self.evaluate_method(func, "feat")

        fp = open("calibration_results/optimal_feature_weights", "w")
        fp.write("%f\tintercept\n" % (intercept))
        for index, feature in enumerate(dense_features):
            if "feat%d" % index in best_features:
                fp.write("%f\t%s\n" % (mult*model.params["feat%d" % index], feature))
        fp.close()

        print "Best feature D: ", D

    def optimise_combination(self):
        """
        Use multiple linear regression to determine the optimal weighted
        combination of the GEOGRAPHIC, GENETIC and FEATUE methods.
        """

        df = {}
        df["auth"] = self.common_auth_combo_vector

        names = ("geo", "gen", "feat")
        funcs = (distance.build_optimal_geographic_matrix,
                    distance.build_optimal_genetic_matrix,
                    distance.build_optimal_feature_matrix)
        for name, func in zip(names, funcs):
            austro_method = self.compute_method_vector(func(self.austrolangs), self.common_austro_langs, self.wals_austro_trans)
            indo_method =  self.compute_method_vector(func(self.indolangs), self.common_indo_langs, self.wals_indo_trans)
            df[name] = np.concatenate([austro_method, indo_method])

        df = pd.DataFrame(df)
        df.to_csv("calibration_results/feature_data.csv")
        model = smf.wls('auth ~ geo + gen + feat', data=df, weights=self.weights).fit()

        fp = open("calibration_results/optimal_combination_weights", "w")
#        fp.write("intercept\t%f\n" % model.params["Intercept"])
        fp.write("intercept\t%f\n" % 0.0)
        fp.write("geo\t%f\n" % model.params["geo"])
        fp.write("gen\t%f\n" % model.params["gen"])
        fp.write("feat\t%f\n" % model.params["feat"])
        fp.close()

#        return (model.params["Intercept"], model.params["geo"], model.params["gen"], model.params["feat"])

        combo_austro = distance.build_optimal_combination_matrix(self.austrolangs)
        combo_indo = distance.build_optimal_combination_matrix(self.indolangs)
        D, intt, mult = self.fit_models(combo_austro, combo_indo, "combo")
        print "best combo D: ", D

        fp = open("calibration_results/optimal_combination_weights", "w")
        fp.write("intercept\t%f\n" % intt)
        print intt
        fp.write("geo\t%f\n" % (mult*model.params["geo"]))
        print mult*model.params["geo"]
        fp.write("gen\t%f\n" % (mult*model.params["gen"]))
        print mult*model.params["gen"]
        fp.write("feat\t%f\n" % (mult*model.params["feat"]))
        print mult*model.params["feat"]
        fp.close()

        return

        return (best_intercept, best_weights[0], best_weights[1], best_weights[2])
        old_D = 1000
        lowest_D = 1000
        weights = [1.0/3, 1.0/3, 1.0/3]
        best_weights = weights[:]
        intercept = 0.5
        best_intercept = 0.5
        for iterations in range(0,10000):
            oldweights = weights[:]
            oldint = intercept
            # change params
            if random.randint(1,100) == 42:
                # Go back to best so far
                weights = best_weights[:]
                intercept = best_intercept
            elif random.randint(1,3) == 1:
                # shuffle weights
                random.shuffle(weights)
            elif random.randint(1,3) == 2:
                # shift weights
                source, target = random.sample([0,1,2],2)
                delta = random.sample([0.01, 0.05, 0.1, 0.2],1)[0]
                if weights[source] > delta:
                    weights[source] -= delta
                    weights[target] += delta
            elif random.randint(1,3) == 3:
                # shift intercept
                delta = random.sample([0.01, 0.05, 0.1, 0.2],1)[0]
                if random.randint(1,2) == 1 and intercept >= delta:
                    intercept -= delta
                elif intercept <= 1.0 - delta:
                    intercept += delta

            observations = [weights[0]*a + weights[1]*b + weights[2]*c for a, b, c in itertools.izip(geo, gen, feat)]
            D, p = scipy.stats.kstest(observations, baselinecdf)
            if D < old_D or random.randint(1,100) < 20:
                # We've improved, or it's a rare backward step
                old_D = D
            else:
                # Keep old value
                weights = oldweights[:]
                intercept = oldint
            if D < lowest_D:
                lowest_D = D
                best_weights = weights
                best_intercept = intercept

#        df = {}
#        df["auth"] = self.auth_combo_vector
#        df["geo"] = np.concatenate([geo_austro, geo_indo])
#        df["gen"] = np.concatenate([gen_austro, gen_indo])
#        df["feat"] = np.concatenate([feat_austro, feat_indo])
#        df = pd.DataFrame(df)
#        df.to_csv("calibration_results/combination_data.csv")
#        model = smf.ols('auth ~ geo + gen + feat', data=df).fit()
#        weights = [model.params[x] for x in ("geo", "gen", "feat")]

        fp = open("calibration_results/optimal_combination_weights", "w")
        fp.write("intercept\t%f\n" % best_intercept)
        fp.write("geo\t%f\n" % best_weights[0])
        fp.write("gen\t%f\n" % best_weights[1])
        fp.write("feat\t%f\n" % best_weights[2])
        fp.close()

        return (best_intercept, best_weights[0], best_weights[1], best_weights[2])

    def summarise(self):

        """
        Compute the correlation between each of the four methods and the
        two authoritative trees, displaying the results in a PrettyTable.
        """

        # Save full pairwise distance lists
        geo_austro = distance.build_optimal_geographic_matrix(self.austrolangs)
        geo_indo =  distance.build_optimal_geographic_matrix(self.indolangs)
        gen_austro = distance.build_optimal_genetic_matrix(self.austrolangs)
        gen_indo =  distance.build_optimal_genetic_matrix(self.indolangs)
        feat_austro = distance.build_optimal_feature_matrix(self.austrolangs)
        feat_indo =  distance.build_optimal_feature_matrix(self.indolangs)
        combo_austro = distance.build_optimal_combination_matrix(self.austrolangs)
        combo_indo =  distance.build_optimal_combination_matrix(self.indolangs)
        geo = [n for n in (itertools.chain(itertools.chain(*gen_austro), itertools.chain(*gen_indo)))]
        gen = [n for n in (itertools.chain(itertools.chain(*gen_austro), itertools.chain(*gen_indo)))]
        feat = [n for n in (itertools.chain(itertools.chain(*feat_austro), itertools.chain(*feat_indo)))]
        combo = [n for n in (itertools.chain(itertools.chain(*combo_austro), itertools.chain(*combo_indo)))]
        df = {}
        df["geo"] = geo
        df["gen"] = gen
        df["feat"] = feat
        df["combo"] = combo
        df = pd.DataFrame(df)
        df.to_csv("calibration_results/all_methods_full.csv")

        # Do PrettyTable
        pt = prettytable.PrettyTable(["Method", "Austronesian correl", "Indo-European correl"])

        names = ("GEOGRAPHIC", "GENETIC", "FEATURE", "COMBINATION")
        little_names = ("geo", "gen", "feat", "combo")
        funcs = (distance.build_optimal_geographic_matrix,
                    distance.build_optimal_genetic_matrix,
                    distance.build_optimal_feature_matrix,
                    distance.build_optimal_combination_matrix)
        df = {}
        df["auth"] = self.common_auth_combo_vector
        for name, little_name, func in zip(names, little_names, funcs):
            austro_method = self.compute_method_vector(func(self.austrolangs), self.common_austro_langs, self.wals_austro_trans)
            austro_auth = self.common_auth_austro_vector
            indo_method =  self.compute_method_vector(func(self.indolangs), self.common_indo_langs, self.wals_indo_trans)
            indo_auth = self.common_auth_indo_vector

            austro_correl = np.corrcoef(austro_method, austro_auth)[0,1]
            indo_correl = np.corrcoef(indo_method, indo_auth)[0,1]
            pt.add_row([name, austro_correl, indo_correl])

            df[little_name] = np.concatenate([austro_method, indo_method])

        df = pd.DataFrame(df)
        df.to_csv("calibration_results/all_methods_common.csv")

        print pt

def main():

    calibrator = Calibrator()
    calibrator.optimise_geographic() 
    calibrator.optimise_genetic()
    calibrator.optimise_feature()
    calibrator.optimise_combination()
    calibrator.summarise()

if __name__ == "__main__":
    main()
