#!/usr/bin/env python
import random
import itertools

def load_stabilities(filename):
    # If we load ALL the stability samples, we get an absurd
    # amount.  So let's randomly sample from 10 trees and then
    # take a random subsample of 500 points
    fp = open(filename, "r")
    stabs = [map(float, line) for line in [x.strip().split() for x in fp.readlines()]]
    fp.close()
    return stabs
    
def parse_ubertree(filename):
    ages = []
    dists = []
    fp = open(filename, "r")
    for line in fp:
        age, dist = line.strip().split(",",1)
        ages.append(float(age))
        dists.append(map(float, dist.split(",")))
    fp.close()
    return ages, dists
   
def parse_sliding_prior_file(method="combination", multitree=False):
    posteriors = {}
    if multitree:
        for family in "afro austro indo niger nilo sino".split():
            posteriors[family] = []
        fp = open("../Inference/results/common-q/unsplit/unshuffled/%s/sliding_prior" % (method), "r")

        for family, line in zip(itertools.cycle("afro austro indo niger nilo sino junk".split()), fp.readlines()):
            if family != "junk":
                posteriors[family].append(map(float, line.strip().split()))
        fp.close()
    else:
        for family in "afro austro indo niger nilo sino".split():
            posteriors[family] = []
            fp = open("../Inference/results/individual-q/unsplit/unshuffled/%s/%s/sliding_prior" % (method, family), "r")
            for line in fp.readlines():
                posteriors[family].append(map(float, line.strip().split()))
            fp.close()

    return posteriors

def parse_summary_file(filename, multitree=False):
    summary = {}
    fp = open(filename, "r")
    lines = [line.strip().split() for line in fp.readlines()]
    fp.close()
    summary["mean_stabs"] = map(float, lines[1])
    summary["mean_trans"] = [map(float, line) for line in lines[4:10]]
    summary["mean_Q"] = [map(float, line) for line in lines[12:18]]
    summary["mean_P"] = [map(float, line) for line in lines[20:26]]
    summary["mean_stationary"] = map(float, lines[28])
    if multitree:
        summary["mean_uniform_ancestral"] = [map(float, line) for line in lines[31:37]]
        summary["mean_fuzzy_ancestral"] = [map(float, line) for line in lines[39:45]]
        summary["mean_stationary_ancestral"] = [map(float, line) for line in lines[47:53]]
    else:
        summary["mean_uniform_ancestral"] = map(float, lines[31])
        summary["mean_fuzzy_ancestral"] = map(float, lines[34])
        summary["mean_stationary_ancestral"] = map(float, lines[37])

    fp.close()
    return summary

def format_vector(vector):
    string = "["
    for v in vector:
        string += ("%f, " % v)
    string += "]"
    return string

def format_matrix(matrix):
    string = "["
    for row in matrix:
        for r in row[:-1]:
            string += ("%f, " % r)
        string += ("%f;\n" % row[-1])
    string += "]"
    return string

def format_summary(summary, method, family=None):
    # Note that 2D asociative arrays are too much for MATLAB's pretty little head,
    # so we simulate them using underscore separation of key strings
    lines = []
    if family:
        key = "%s_%s" % (method, family)
    else:
        key = method
    lines.append("stabs('%s') = %s\n" % (key, format_vector(summary["mean_stabs"])))
    lines.append("trans('%s') = %s\n" % (key, format_matrix(summary["mean_trans"])))
    if family:
        lines.append("indiv_uniform_ancestrals('%s') =  %s;\n" % (key, format_vector(summary["mean_uniform_ancestral"])))
        lines.append("indiv_fuzzy_ancestrals('%s') =  %s;\n" % (key, format_vector(summary["mean_fuzzy_ancestral"])))
        lines.append("indiv_stationary_ancestrals('%s') =  %s;\n" % (key, format_vector(summary["mean_stationary_ancestral"])))
    else:
        lines.append("Q('%s') =  %s\n" % (key, format_matrix(summary["mean_Q"])))
        for i, family in enumerate("afro austro indo niger nilo sino".split()):
            key = "%s_%s" % (method, family)
            lines.append("multi_uniform_ancestrals('%s') =  %s;\n" % (key, format_vector(summary["mean_uniform_ancestral"][i])))
            lines.append("multi_fuzzy_ancestrals('%s') =  %s;\n" % (key, format_vector(summary["mean_fuzzy_ancestral"][i])))
            lines.append("multi_stationary_ancestrals('%s') =  %s;\n" % (key, format_vector(summary["mean_stationary_ancestral"][i])))
    return "\n".join(lines)

def main():
    fp = open("results.m", "w")
    for var in "stabs all_stabs trans indiv_uniform_ancestrals indiv_fuzzy_ancestrals indiv_stationary_ancestrals multi_uniform_ancestrals multi_fuzzy_ancestrals multi_stationary_ancestrals indiv_sliding_priors multi_sliding_priors Q".split():
        fp.write("%s = containers.Map();\n" % (var,))

    for method in "geographic genetic feature combination".split():

        # SLIDING PRIOR STUFF

        sliding_priors = parse_sliding_prior_file(method, multitree=False)
        for family in "afro austro indo niger nilo sino".split():
            fp.write("indiv_sliding_priors('%s_%s') = %s;\n" % (method, family, format_matrix(sliding_priors[family])))

        sliding_priors = parse_sliding_prior_file(method, multitree=True)
        for family in "afro austro indo niger nilo sino".split():
            fp.write("multi_sliding_priors('%s_%s') = %s;\n" % (method, family, format_matrix(sliding_priors[family])))

        # INDIVIDUAL Q
        for family in "afro austro indo niger nilo sino".split():
            summary = parse_summary_file("../Inference/results/individual-q/unsplit/unshuffled/%s/%s/summary" % (method, family))
            fp.write(format_summary(summary, method, family))
            stabs = load_stabilities("../Inference/results/individual-q/unsplit/unshuffled/%s/%s/stabilities" % (method, family))
            fp.write("all_stabs('%s_%s') = %s;\n" % (method, family, format_matrix(stabs)))

        # COMMON Q

        summary = parse_summary_file("../Inference/results/common-q/unsplit/unshuffled/%s/summary" % (method,), multitree=True)
        fp.write(format_summary(summary, method))
        stabs = load_stabilities("../Inference/results/common-q/unsplit/unshuffled/%s/stabilities" % (method,))
        fp.write("all_stabs('%s') = %s;\n" % (method, format_matrix(stabs)))
        

    # UBER TREE
    ages, dists = parse_ubertree("../Inference/results/common-q/unsplit/unshuffled/combination/common_ancestor")
    fp.write("common_ancestor_ages = %s;\n" % format_vector(ages))
    fp.write("common_ancestor_dists = %s;\n" % format_matrix(dists))

    fp.close()

if __name__ == "__main__":
    main()

