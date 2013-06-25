import random
import itertools

def load_stabilities(basename, multitree=False):
    # If we load ALL the stability samples, we get an absurd
    # amount.  So let's randomly sample from 10 trees and then
    # take a random subsample of 500 points
    stabs = []
    trees = random.sample(range(1,101), 10)
    for i in trees:
        if multitree:
            filename = basename + "/trees_%d/samples" % i
        else:
            filename = basename + "/tree_%d/samples" % i
        fp = open(filename, "r")
        lines = fp.readlines()
        for line in [lines[x] for x in range(2,len(lines),10)]:
            stabs.append(map(float, line.strip().split()))
        fp.close()
    stabs = random.sample(stabs, 500)
    return stabs
    
def parse_ubertree():
    ages = []
    dists = []
    fp = open("../Inference/ubertree_dists", "r")
    for line in fp:
        age, dist = line.strip().split(",",1)
        ages.append(float(age))
        dists.append(map(float, dist.split(",")))
    fp.close()
    return ages, dists
   
def parse_sliding_prior_file():
    posteriors = {}
    for family in "afro austro indo niger nilo sino".split():
        posteriors[family] = []
    fp = open("../Inference/results/common-q/combination/sliding_prior", "r")

    for family, line in zip(itertools.cycle("afro austro indo niger nilo sino junk".split()), fp.readlines()):
        if family != "junk":
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
    lines.append("stabs(%s) = %s\n" % (key, format_vector(summary["mean_stabs"])))
    lines.append("trans(%s) = %s\n" % (key, format_matrix(summary["mean_trans"])))
    if family:
        lines.append("indiv_uniform_ancestrals('%s') =  %s\n" % (key, format_vector(summary["mean_uniform_ancestral"])))
        lines.append("indiv_fuzzy_ancestrals('%s') =  %s\n" % (key, format_vector(summary["mean_fuzzy_ancestral"])))
        lines.append("indiv_stationary_ancestrals('%s') =  %s\n" % (key, format_vector(summary["mean_stationary_ancestral"])))
    else:
        for i, family in enumerate("afro austro indo niger nilo sino".split()):
            key = "%s_%s" % (method, family)
            lines.append("multi_uniform_ancestrals('%s') =  %s\n" % (key, format_vector(summary["mean_uniform_ancestral"][i])))
            lines.append("multi_fuzzy_ancestrals('%s') =  %s\n" % (key, format_vector(summary["mean_fuzzy_ancestral"][i])))
            lines.append("multi_stationary_ancestrals('%s') =  %s\n" % (key, format_vector(summary["mean_stationary_ancestral"][i])))
    return "\n".join(lines)

def main():
    fp = open("results.m", "w")
    for var in "stabs all_stabs trans indiv_uniform_ancestrals indiv_fuzzy_ancestrals indiv_stationary_ancestrals multi_uniform_ancestrals multi_fuzzy_ancestrals mutli_stationary_ancestrals sliding_priors".split():
        fp.write("%s = containers.Map()\n" % (var,))

    for method in "geographic genetic feature combination".split():
        for family in "afro austro indo niger nilo sino".split():
            summary = parse_summary_file("../Inference/results/individual-q/%s/%s/summary" % (method, family))
            fp.write(format_summary(summary, method, family))
            stabs = load_stabilities("../Inference/results/individual-q/%s/%s/" % (method, family))
            fp.write("all_stabs('%s_%s') = %s\n" % (method, family, format_matrix(stabs)))

        summary = parse_summary_file("../Inference/results/common-q/%s/summary" % (method,), multitree=True)
        fp.write(format_summary(summary, method))
        stabs = load_stabilities("../Inference/results/common-q/%s/" % (method,), multitree=True)
        fp.write("all_stabs('%s') = %s\n" % (method, format_matrix(stabs)))

    ages, dists = parse_ubertree()
    fp.write("common_ancestor_ages = %s\n" % format_vector(ages))
    fp.write("common_ancestor_dists = %s\n" % format_matrix(dists))

    sliding_posteriors = parse_sliding_prior_file()
    for family in "afro austro indo niger nilo sino".split():
        fp.write("sliding_priors('%s') = %s\n" % (family, format_matrix(sliding_posteriors[family])))
    fp.close()

if __name__ == "__main__":
    main()

