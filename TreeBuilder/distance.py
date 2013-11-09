from math import cos, sin, tan, atan, sqrt, pi, radians
import math

import codecs
from random import shuffle, random
from csv import DictReader
import sys

MINDIST = {}
MINDIST["geographic"] = 0.15 # was 001
MINDIST["genetic"] = 0.30
MINDIST["feature"] = 0.45

def build_matrix(languages, distance_function):
    matrix = []
    for x in range(0,len(languages)):
        matrix.append([0]*len(languages))
    for i in range(0,len(languages)):
        for j in range(i,len(languages)):
            if i == j:
                continue
            matrix[i][j] = distance_function(languages[i],languages[j])
            matrix[j][i] = matrix[i][j]
    return matrix

def build_matrix_by_method_name(languages, method):
    if method == "geographic":
        return built_optimal_geographic_matrix(languages)
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
	if a == 1.0:
		print "Ack, cheating on: ", lat1, lon1, lat2, lon2
		return 0.5
        c = 2*atan(sqrt(a)/sqrt(1-a))
        d = R*c
        d = 2*d/(2*pi*R)
        return d

def geographic_function_factory(languages, param):
    done = []
    maxdist = 0
    for i in range(0,len(languages)):
        for j in range(i,len(languages)):
            if i != j and (j,i) not in done:
                done.append((i,j))
                dist = haversine_distance(languages[i].data["location"], languages[j].data["location"])
                if dist > maxdist:
                    maxdist = dist
    def distance(lang1, lang2):
        return haversine_distance(lang1.data["location"], lang2.data["location"]) / maxdist
    return distance

def geographic_matrix_factory(param):
    def matrix_builder(languages):
        distfunc = geographic_function_factory(languages, param)
        return build_matrix(languages, distfunc)
    return matrix_builder

def build_optimal_geographic_matrix(languages):
    distfunc = geographic_function_factory(languages, None)
    return build_matrix(languages, distfunc)

##########
# GENETIC STUFF
##########

def genetic_function_factory(param):
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
        return distance
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
    with open("calibration_results/optimal_genetic_parameter", "r") as fp:
        optimal_param = float(fp.readline().strip())
    genetic_distance = genetic_function_factory(optimal_param)
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

def feature_function_factory(weights=None):
    """
    Return a function which computes pairwise distances according
    to the FEATURE method, with feature weights as given by weights.
    """
    comparators = build_comparators()
    if not weights:
        weights = {}
    def feature_distance(lang1, lang2):
        dist = 0
        norm = 0
        for key in lang1.data:
            if key not in ["genus", "subfamily", "family", "location", "Order of Subject, Object and Verb", "Order of Subject and Verb", "Order of Object and Verb", "iso_codes", "ethnoclass"] and not key.startswith("Relationship between the Order of Object and Verb"):
               dist += weights.get(key, 1.0)*comparators[key](lang1.data[key], lang2.data[key])
               norm += weights.get(key, 1.0)
        return dist/norm
    return feature_distance

def feature_matrix_factory(weights=None):
    def matrix_builder(languages):
        distfunc = feature_function_factory(weights)
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
    feature_distance = feature_function_factory(weights)
    return build_matrix(languages, feature_distance)

##########
# COMBINATION STUFF
##########

def build_optimal_combination_matrix(languages):
    geo = geographic_function_factory()

    with open("calibration_results/optimal_genetic_parameter", "r") as fp:
        optimal_param = float(fp.readline().strip())
    gen = genetic_function_factory(optimal_param)

    weights = {}
    with open("calibration_results/optimal_feature_weights", "r") as fp:
        for line in fp:
            weight, feature = line.split("\t")
            weight = float(weight.strip())
            feature = feature.strip()
            weights[feature] = weight
    feat = feature_function_factory(weights)

    fp = open("calibration_results/optimal_combination_weights", "r")
    geo_w, gen_w, feat_w = [float(line.split()[1]) for line in fp.readlines()]
    fp.close()
    def distfunc(l1, l2):
        return geo_w*geo(l1, l2) + gen_w*gen(l1, l2) + feat_w*feat(l1, l2)
    return build_matrix(languages, distfunc)

