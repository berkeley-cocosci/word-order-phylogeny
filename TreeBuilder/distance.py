from math import cos, sin, tan, atan, sqrt, pi, radians
import codecs
from random import shuffle, random
from csv import DictReader

from diagnostics import run_diags

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

def genetic_distance(lang1, lang2, comparators):
	if(lang1.data["genus"] and lang2.data["genus"] and lang1.data["genus"] == lang2.data["genus"]):
		# Known matching genera - minimum distance
		return 1.0
	if(lang1.data["subfamily"] and lang2.data["subfamily"] and lang1.data["subfamily"] == lang2.data["subfamily"]):
		# Known matching subfamilies
		# What about genera?
		if(lang1.data["genus"] and lang2.data["genus"] and lang1.data["genus"] != lang2.data["genus"]):
			# Definitely not matching genera
			return 2.0
		else:
			# Uncertain about genera
			return 1.5
	if(lang1.data["family"] == lang2.data["family"]):
		# Known matching families
		# What about subfamilies?
		if(lang1.data["subfamily"] and lang2.data["subfamily"] and lang1.data["subfamily"] != lang2.data["subfamily"]):
			# Definitely not matching subfamilies
			return 3.0
		else:
			# Uncertain about subfamily
			return 2.5
	return 4.0

def geographic_distance(lang1, lang2, comparators):
	return genetic_distance(lang1, lang2, comparators) + comparators["location"](lang1.data["location"], lang2.data["location"])

def feature_distance(lang1, lang2, comparators):
	feature_distance = 0
	norm = 0
	for key in lang1.data:
		if key not in ["genus", "subfamily", "family", "location", "Order of Subject, Object and Verb", "isocode"]:
			feature_distance += comparators[key](lang1.data[key], lang2.data[key])
			norm += 1.0
	return geographic_distance(lang1, lang2, comparators) + feature_distance/norm

def build_matrix(languages, method):
	comparators = build_comparators()
	feat_data = build_feat_data()
	if method == "genetic":
		distance = genetic_distance
	elif method == "geographic":
		distance = geographic_distance
	if method == "feature":
		distance = feature_distance
	matrix = []
	for x in range(0,len(languages)):
		matrix.append([0]*len(languages))
	for i in range(0,len(languages)):
		for j in range(i,len(languages)):
			if i == j:
				continue
			matrix[i][j] = distance(languages[i],languages[j],comparators)
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


def build_feat_data():
	fp = open("../WALS2SQL/wals_data/features.csv","r")
	reader = DictReader(fp)
	data = {}
	for row in reader:
		id = row["id"]
		name = row["name"]
		data[id] = name
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

def make_matrix_go_now(languages, method):
	newlangs = []
	print "I had: ", len(languages)
	for lang in languages:
		if lang.name not in [u'Fiote', u'Jiarong']:
			newlangs.append(lang)
	print len(languages), len(newlangs)
	languages = newlangs
	matrix = build_matrix(languages, method)
	return matrix
