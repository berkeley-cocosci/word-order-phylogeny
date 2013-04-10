import getopt, sys
import newick.tree

n = 0

# add_edge(subtree, bootstrap, length)
# I need to give the binarise constructor a list of tuples containing leaf ndoes AND distances!!!
def binarise(nodes):
    good_tree = newick.tree.Tree()
    n = len(nodes)
    if n == 1:
       if isinstance(nodes[0], newick.tree.Leaf):
           return nodes[0]
       else:
           return binarise(nodes[0].get_edges())
    elif n == 2:
       good_tree.add_edge((binarise([nodes[0][0],]), None, nodes[0][2]))
       good_tree.add_edge((binarise([nodes[1][0],]), None, nodes[1][2]))
    elif n == 3:
       # 3 branches - one on left, two on right
       good_tree.add_edge((binarise([nodes[0][0],]), None, nodes[0][2]))
       good_tree.add_edge((binarise(nodes[1:]), None, 0.0))
    elif n == 4:
       # 4 branches - two on left, two on right
       good_tree.add_edge((binarise(nodes[0:2]), None, 0.0))
       good_tree.add_edge((binarise(nodes[2:]), None, 0.0))
    elif n == 5:
       # 5 branches - two on left, three on right
       good_tree.add_edge((binarise(nodes[0:2]), None, 0.0))
       good_tree.add_edge((binarise(nodes[2:]), None, 0.0))
    elif n == 6:
       # 6 branches - three on left, three on right
       good_tree.add_edge((binarise(nodes[0:3]), None, 0.0))
       good_tree.add_edge((binarise(nodes[3:]), None, 0.0))
    elif n == 7:
       # 7 branches - three on left, four on right
       good_tree.add_edge((binarise(nodes[0:3]), None, 0.0))
       good_tree.add_edge((binarise(nodes[3:]), None, 0.0))
    elif n == 8:
       # 8 branches - four on left, four on right
       good_tree.add_edge((binarise(nodes[0:4]), None, 0.0))
       good_tree.add_edge((binarise(nodes[4:]), None, 0.0))
    else:
        print "WTF?!"
    return good_tree

def forkcounter(t, fp):
    edges = t.get_edges()
    #print len(edges)
    for edge in edges:
        if isinstance(edge[0], newick.tree.Tree):
            forkcounter(edge[0], fp)

def callback(t, fp):
    global n
    m = n
    #print len(t.get_edges())
    left, right = t.get_edges()
    n += 1
    if isinstance(left[0], newick.tree.Leaf):
        fp.write("%d, 0, %d, %f, %s\n" % (m, n, left[2], left[0].identifier))
    else:
        fp.write("%d, 0, %d, %f, %s\n" % (m, n, left[2], "NON-LEAF"))
        callback(left[0], fp)
    n += 1
    if isinstance(right[0], newick.tree.Leaf):
        fp.write("%d, 1, %d, %f, %s\n" % (m, n, right[2], right[0].identifier))
    else:
        fp.write("%d, 1, %d, %f, %s\n" % (m, n, right[2], "NON-LEAF"))
        callback(right[0], fp)
    
def format_tree(t, output):
    fp = open(output, "w")
    callback(t, fp)
    fp.close()

class BranchLengthSum(newick.AbstractHandler):
    def __init__(self):
        self.sum = 0.0

    def new_edge(self,b,l):
        if l:
            self.sum += l

    def get_result(self):
        return self.sum

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:o:")
    except getoptGetoptError, err:
        print str(err)
	sys.exit(2)
    input = None
    output = None
    print opts
    for o, a in opts:
	print o, a
        if o == "-i":
            input = a
        elif o == "-o":
            output = a
    
    print input, output   
    if not (input and output):
        print "I need to know what to do!"
        sys.exit(2)

    fp = open(input, "r")
    tree_string = fp.read()
    fp.close()
    raw_tree = newick.tree.parse_tree(tree_string)
    #for leaf in raw_tree.get_leaves():
    #   print leaf.identifier
    cooked_tree = binarise([raw_tree,])
    format_tree(cooked_tree, output)

if __name__ == "__main__":
    main()
