import itertools
from math import cos, sin, tan, atan, sqrt, pi, radians
import math
import numpy as np

import codecs
from random import shuffle, random
from csv import DictReader
import sys
import pdb

import numpy as np

def build_matrix(languages, distance_function):
    n = len(languages)
    matrix = np.zeros(shape=(n, n))
    for i, j in itertools.combinations(range(0,n), 2):
        d = distance_function(languages[i], languages[j])
        matrix[i][j] = d
        matrix[j][i] = d
    maxd = matrix.max()
    mind = matrix.min()
    if mind < 0:
        matrix /= abs(mind)
    else:
        matrix /= maxd
    return matrix

def build_matrix_by_method_name(languages, method):
    if method == "geographic":
        return build_optimal_geographic_matrix(languages)
    elif method == "genetic":
        return build_optimal_genetic_matrix(languages)
    elif method == "feature":
        return build_optimal_feature_matrix(languages)
    elif method == "combination":
        return build_optimal_combination_matrix(languages)

##########
# GEOGRAPHIC STUFF
##########

def haversine_distance(val1, val2):
    R = 6371.0
    lat1, lon1 = val1
    lat2, lon2 = val2
    dlat = radians(lat2-lat1)
    dlon = radians(lon2-lon1)
    lat1 = radians(lat1)
    lat2 = radians(lat2)
    a = sin(dlat/2)**2 + cos(lat1)*cos(lat2)*sin(dlon/2)**2
    c = 2*atan(sqrt(a)/sqrt(1-a))
    d = R*c
    d = 2*d/(2*pi*R)
    return d

def geographic_function_factory(languages, intercept=0.0, slope=1.0):
    done = []
    maxdist = 0
    mindist = 1000
    for i in range(0,len(languages)):
        for j in range(i,len(languages)):
            if i != j and (j,i) not in done:
                done.append((i,j))
                dist = haversine_distance(languages[i].data["location"], languages[j].data["location"])
                maxdist = max(maxdist, dist)
                mindist = min(mindist, dist)
    def distance(lang1, lang2):
        linear = haversine_distance(lang1.data["location"], lang2.data["location"]) / maxdist
        log = np.log(linear)
        minlog = abs(np.log(mindist))
        log /= minlog
        return intercept + slope*log
    return distance

def geographic_matrix_factory():
    def matrix_builder(languages):
        distfunc = geographic_function_factory(languages)
        return build_matrix(languages, distfunc)
    return matrix_builder

def build_optimal_geographic_matrix(languages):
    with open("calibration_results/optimal_geographic_parameters", "r") as fp:
        intercept = float(fp.readline().strip())
        slope = float(fp.readline().strip())
    geographic_distance = geographic_function_factory(languages, intercept, slope)
    return build_matrix(languages, geographic_distance)

##########
# GENETIC STUFF
##########

def genetic_function_factory(param, intercept=0.0, slope=1.0):
    """
    Return a function which measures pairwise distances calculated
    using the GENETIC method with parameter param.
    """
    def genetic_distance(lang1, lang2):
        lang1hier = lang1.data["ethnoclass"].split(",")
        lang2hier = lang2.data["ethnoclass"].split(",")
        n = min(len(lang1hier), len(lang2hier))
        similarity = 0
        for i in range(0,n):
            if lang1hier[i] == lang2hier[i]:
                similarity += param**(i+1)
        maxsim = sum([param**(i+1) for i in range(0,n)])
        distance = (maxsim - similarity) / maxsim
        return intercept + slope*distance
    return genetic_distance

def genetic_matrix_factory(param):
    """
    Return a function which takes as input a list of languages and
    produces as output a matrix of pairwise language distances
    calcualted using the GENETIC method with parameter param.
    """
    def matrix_builder(languages):
        distfunc = genetic_function_factory(param)
        return build_matrix(languages, distfunc)
    return matrix_builder

def build_optimal_genetic_matrix(languages):
    with open("calibration_results/optimal_genetic_parameters", "r") as fp:
        intercept = float(fp.readline().strip())
        slope = float(fp.readline().strip())
        param = float(fp.readline().strip())
    genetic_distance = genetic_function_factory(param, intercept, slope)
    return build_matrix(languages, genetic_distance)

##########
# FEATURE STUFF
##########

def euclidean_factory(max):
	def func(x,y):
		if x == None or y == None:
			return 0.5
		else:
			x = int(x)
			y = int(y)
			return abs(x-y)/(max*1.0)
	return func

def hamming_distance(val1, val2):
	if val1 == None or val2 == None:
		return 0.5
	if val1 == val2:
		return 0.0
	return 1.0

def build_comparators():
	fp = open("comparators.csv","r")
	reader = DictReader(fp)
	comparators = {}
	for row in reader:
		if row["comparator"] == "euclidean":
			func = euclidean_factory(int(row["comparg"]))
		else:
			func = hamming_distance
		comparators[row["name"]] = func
	fp.close()
	for x in "genus family subfamily".split():
		comparators[x] = hamming_distance
	comparators["location"] = haversine_distance
	return comparators

def feature_function_factory(weights, intercept=0, slope=1.0):
    """
    Return a function which computes pairwise distances according
    to the FEATURE method, with feature weights as given by weights.
    """
    comparators = build_comparators()
    def feature_distance(lang1, lang2):
        dist = 0
        norm = 0
        for feature in weights:
            if feature not in ["genus", "subfamily", "family", "location", "Order of Subject, Object and Verb", "Order of Subject and Verb", "Order of Object and Verb", "iso_codes", "ethnoclass"] and not feature.startswith("Relationship between the Order of Object and Verb"):
                if feature in lang1.data and feature in lang2.data:
                    dist += weights[feature]*comparators[feature](lang1.data[feature], lang2.data[feature])
                else:
                    dist += weights[feature]*0.5
                norm += weights[feature]
                raw = dist/norm
        return intercept + slope*raw
    return feature_distance

def feature_matrix_factory(weights, intercept=0, slope=1.0):
    def matrix_builder(languages):
        distfunc = feature_function_factory(weights, intercept, slope)
        return build_matrix(languages, distfunc)
    return matrix_builder

def build_optimal_feature_matrix(languages):
    weights = {}
    with open("calibration_results/optimal_feature_weights", "r") as fp:
        for line in fp:
            weight, feature = line.split("\t")
            weight = float(weight.strip())
            feature = feature.strip()
            weights[feature] = weight
    feature_distance = feature_function_factory(weights, intercept=weights["intercept"])
    return build_matrix(languages, feature_distance)

##########
# COMBINATION STUFF
##########

def build_optimal_combination_matrix(languages):

    with open("calibration_results/optimal_geographic_parameters", "r") as fp:
        intercept = float(fp.readline().strip())
        slope = float(fp.readline().strip())
    geo = geographic_function_factory(languages, intercept, slope)

    with open("calibration_results/optimal_genetic_parameters", "r") as fp:
        intercept = float(fp.readline().strip())
        slope = float(fp.readline().strip())
        param = float(fp.readline().strip())
    gen = genetic_function_factory(param, intercept, slope)

    weights = {}
    with open("calibration_results/optimal_feature_weights", "r") as fp:
        for line in fp:
            weight, feature = line.split("\t")
            weight = float(weight.strip())
            feature = feature.strip()
            weights[feature] = weight
    feat = feature_function_factory(weights, intercept=weights["intercept"])

    fp = open("calibration_results/optimal_combination_weights", "r")
    intercept, geo_w, gen_w, feat_w = [float(line.split()[1]) for line in fp.readlines()]
    fp.close()
    def distfunc(l1, l2):
        return intercept + geo_w*geo(l1, l2) + gen_w*gen(l1, l2) + feat_w*feat(l1, l2)
    return build_matrix(languages, distfunc)

