def get_ancestral_dist(filename):
	fp = open(filename, "r")
	for i in range(0,3):
		fp.readline()
	dist = map(float, fp.readline().strip().split())
	fp.close()
	return dist

def get_stationary_dist(filename):
	fp = open(filename, "r")
	for i in range(0,33):
		fp.readline()
	dist = map(float, fp.readline().strip().split())
	fp.close()
	return dist

def get_transition_matrix(filename):
	matrix = []
	fp = open(filename, "r")
	for i in range(0,25):
		fp.readline()
	for i in range(0,6):
		matrix.append(map(float, fp.readline().strip().split()))
	fp.close()
	return matrix

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

def main():
	SAMPLES = 100
	for family in "afro austro indo niger nilo sino".split():
		for method in "distance family feature".split():
			ancestral = [0, 0, 0, 0, 0, 0]
			stationary = [0, 0, 0, 0, 0, 0]
			transition = []
			for i in range(0,6):
				transition.append([0, 0, 0, 0, 0, 0])
			for n in range(1,SAMPLES + 1):
				sample = get_ancestral_dist("results/results_%s_%s_%d" % (family, method, n))
				for i in range(0, 6):
					ancestral[i] += sample[i] / (1.0*SAMPLES)
				sample = get_stationary_dist("results/results_%s_%s_%d" % (family, method, n))
				for i in range(0, 6):
					stationary[i] += sample[i] / (1.0*SAMPLES)
				sample = get_transition_matrix("results/results_%s_%s_%d" % (family, method, n))
				for i in range(0, 6):
					for j in range(0, 6):
						transition[i][j] += sample[i][j] / (1.0*SAMPLES)
			print "%s_%s_ancestral = %s\n\n" % (family, method, format_vector(ancestral))
			print "%s_%s_stationary = %s\n\n" % (family, method, format_vector(stationary))
			print "%s_%s_transition = %s\n\n" % (family, method, format_matrix(transition))


if __name__ == "__main__":
	main()

