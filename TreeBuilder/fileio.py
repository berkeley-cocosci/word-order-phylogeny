import codecs

def save_translate_file(languages, filename):
    fp = codecs.open(filename, "w", "utf8")
    fp.write("translate\n");
    fp.write("%d %s" % (0, languages[0].name)) 
    for i in range(1,len(languages)):
        fp.write(",\n%d %s" % (i, languages[i].name)) 
    fp.write(";\n")
    fp.close()

def load_translate_file(filename):
    langs = {}
    fp = codecs.open(filename, "r", "utf8")
    line = fp.readline().strip()
    for line in fp:
        index, lang = [x.strip() for x in line.split()]
        langs[lang[:-1]] = int(index)
    return langs

def save_multistate_file(languages, filename):
    fp = codecs.open(filename, "w", "utf8")
    for i, lang in enumerate(languages):
        fp.write("%d\t%d\n" % (i, lang.data[u'Order of Subject, Object and Verb']))
    fp.close()

def save_matrix(matrix, filename):
    fp = codecs.open(filename,encoding="utf-8",mode="w")
    fp.write("%d\n" % len(matrix))
    for i in range(0,len(matrix)):
        fp.write("%d " % i)
        for j in range(0, len(matrix)):
            if i == j:
                fp.write("0.0 ")
            else:
                fp.write(" %f " % (matrix[i][j]))
        fp.write("\n")
    fp.close()

def load_matrix(filename):
    fp = codecs.open(filename,encoding="utf-8",mode="r")
    n = int(fp.readline().strip())
    matrix = []
    for line in fp:
        matrix.append(map(float, line.strip().split())[1:])
    fp.close()
    return matrix
