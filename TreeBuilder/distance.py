from math import cos, sin, tan, atan, sqrt, pi, radians
import math

import codecs
from random import shuffle, random
from csv import DictReader
import sys

from diagnostics import run_diags

MINDIST = {}
MINDIST["geographic"] = 0.15 # was 001
MINDIST["genetic"] = 0.30
MINDIST["feature"] = 0.45

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

def param_genetic_factory(param):
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

genetic_distance = param_genetic_factory(0.65)

#def geographic_distance(lang1, lang2):
#	return param_genetic_factory((0.1, 0.5,))(lang1, lang2) + haversine_distance(lang1.data["location"], lang2.data["location"])

#def weighted_gen_geo_factory(params):
#    assert len(params) == 2
#    gen_weight, gen_param = params
#    gen = param_genetic_factory((gen_param,))
#    geo = puregeo_lin
#    return lambda x,y: gen_weight*gen(x,y) + (1-gen_weight)*geo(x,y)

def feature_distance_factory():
    comparators = build_comparators()
    def feature_distance(lang1, lang2):
        feature_distance = 0
        norm = 0
        for key in lang1.data:
            if key not in ["genus", "subfamily", "family", "location", "Order of Subject, Object and Verb", "Order of Subject and Verb", "Order of Object and Verb", "iso_codes", "ethnoclass"] and not key.startswith("Relationship between the Order of Object and Verb"):
               feature_distance += comparators[key](lang1.data[key], lang2.data[key])
               norm += 1.0
        return feature_distance/norm
    return feature_distance

def build_matrix_by_method_name(languages, method):
    if method == "genetic":
        return genetic_matrix(languages)
    elif method == "geographic":
        return linear_geography_matrix(languages)
    elif method == "feature":
        return feature_matrix(languages)
    elif method == "combination":
        return optimal_combination_matrix(languages)

def linear_geography_factory(languages):
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

def linear_geography_matrix(languages):
    distfunc = linear_geography_factory(languages)
    return build_matrix(languages, distfunc)

def genetic_matrix(languages):
    return build_matrix(languages, genetic_distance)

#def genetic_matrix_factory(params):
#    distfunc = param_genetic_factory(params)
#    def matrix_builder(languages):
#        return build_matrix(languages, distfunc)
#    return matrix_builder

def feature_matrix(languages):
    distfunc = feature_distance_factory()
    return build_matrix(languages, distfunc)

def optimal_combination_matrix(languages):
    geo = linear_geography_factory(languages)
    gen = genetic_distance
    feat = feature_distance_factory()
    def distfunc(l1, l2):
#        return 0.1115*geo(l1, l2) + 0.7657*gen(l1, l2) + 0.1228*feat(l1, l2)
        return 0.1348*geo(l1, l2) + 0.6951*gen(l1, l2) + 0.1700*feat(l1, l2)
    return build_matrix(languages, distfunc)

def weighted_triple_factory(weights):
    assert len(weights) == 3
    def matrix_builder(languages):
        distfunc1 = linear_geography_factory(languages)
        distfunc2 = genetic_distance
        distfunc3 = feature_distance_factory()
        def distfunc(lang1, lang2):
            return weights[0]*distfunc1(lang1, lang2) + weights[1]*distfunc2(lang1, lang2) + weights[2]*distfunc3(lang1, lang2)
        return build_matrix(languages, distfunc)
    return matrix_builder

#def build_lang_data():
#	fp = open("languages.csv","r")
#	reader = DictReader(fp)
#	data = {}
#	for row in reader:
#		code = row["wals code"]
#		data[code] = {}
#		for key in row.keys():
#			data[code][key] = row[key] 
#	fp.close()
#	return data

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

