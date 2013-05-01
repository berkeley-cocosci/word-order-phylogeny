#!/usr/bin/env python
import os

import fileio

def run_diags(matrix, languages, filename):

    fp = open(filename, "w")
    hits = 0
    for i in range(0,len(languages)):
        if languages[i].data["family"] != "Indo-European":
            continue
        min_dist = 9999.0
        for j in range(0, len(matrix[i])):
            if j==i:
                continue
            #print matrix[i][j]
            if matrix[i][j] < min_dist:
               # print "In!"
                min_dist = matrix[i][j]
                min_index = j
                min_name = languages[j].name.encode("UTF-8")
        fp.write("Closest to %s is %s\n" % (languages[i].name.encode("UTF-8"), min_name))


#        if languages[i].data["family"] == languages[min_index].data["family"]:
#            hits += 1

    fp.close()
    #print "Got %d correct out of %d\n" % (hits, len(languages))

if __name__ == "__main__":

    translate = fileio.load_translate_file("generated_trees/indo.translate")
    print translate
    langs = ("Dutch", "Danish", "Swedish", "French", "Russian", "Ukranian", "Hindi")
    distances = {}
    for method in ("geographic", "genetic", "feature"):
        for lang in langs:
            distances[lang] = {}
        for index in range(0,100):
            filename = os.path.join("generated_trees", method, "indo", "tree_%d" % (index+1))
            matrix = fileio.load_matrix(filename)
            for lang in langs:
                distances[lang].append(matrix(translate["English"],translate[lang]))
        for lang in langs:
            fp = open(method+"_"+lang.lower()+"_diag","w")
            for dist in distances[lang]:
                fp.write("%f\n" % dist)
            fp.close()

    main()
