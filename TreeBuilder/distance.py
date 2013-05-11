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

def puregeo_lin(lang1, lang2):
	return haversine_distance(lang1.data["location"], lang2.data["location"])

def puregeo_exp(lang1, lang2):
	return exp(haversine_distance(lang1.data["location"], lang2.data["location"]))

def param_genetic_factory(params):
    assert len(params) == 1
    base = params[0]
    def genetic_distance(lang1, lang2):
        lang1hier = lang1.data["ethnoclass"].split(",")
        lang2hier = lang2.data["ethnoclass"].split(",")
        n = min(len(lang1hier), len(lang2hier))
        similarity = 0
        for i in range(0,n):
            if lang1hier[i] == lang2hier[i]:
                similarity += base**(i+1)
        maxsim = sum([base**(i+1) for i in range(0,n)])
        distance = (maxsim - similarity) / maxsim
        return 0.30 + distance
    return genetic_distance

def geographic_distance(lang1, lang2):
	return param_genetic_factory((0.5,))(lang1, lang2) + haversine_distance(lang1.data["location"], lang2.data["location"])

def feature_distance_factory(comparators):
    comparators = build_comparators()
    def feature_distance(lang1, lang2):
        feature_distance = 0
        norm = 0
        for key in lang1.data:
            if key not in ["genus", "subfamily", "family", "location", "Order of Subject, Object and Verb", "iso_codes", "ethnoclass"]:
                feature_distance += comparators[key](lang1.data[key], lang2.data[key])
                norm += 1.0
        return geographic_distance(lang1, lang2) + feature_distance/norm
    return feature_distance

def build_matrix_by_method_name(languages, method):
    if method == "genetic":
        distance = param_genetic_factory((0.5,))
    elif method == "puregeo_lin":
        distance = puregeo_lin
    elif method == "puregeo_exp":
        distance = puregeo_exp
    elif method == "geographic":
        distance = geographic_distance
    if method == "feature":
        distance = feature_distance_factory()
    return build_matrix(languages, distance)

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

def build_lang_data():
	fp = open("languages.csv","r")
	reader = DictReader(fp)
	data = {}
	for row in reader:
		code = row["wals code"]
		data[code] = {}
		for key in row.keys():
			data[code][key] = row[key] 
	fp.close()
	return data

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

