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
    return "\n".join(lines)

def main():
    fp = open("results.m", "w")
    for var in "stabs trans uniform_ancestrals".split():
        fp.write("%s = containers.Map()\n" % (var,))

    for method in "geographic genetic feature combination".split():
        for family in "afro austro indo niger nilo sino".split():
            summary = parse_summary_file("results/individual-q/%s/%s/summary" % (method, family))
            fp.write(format_summary(summary, method, family))
        summary = parse_summary_file("results/common-q/%s/summary" % (method,), multitree=True)
    fp.close()

if __name__ == "__main__":
    main()

