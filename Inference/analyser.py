ORDERS = "SOV SVO VSO VOS OVS OSV".split()

def get_results(filename):
	fp = open(filename, "r")
	fp.readline()
	ancestral = map(float, fp.readline().strip().split())
	fp.readline()
	stabilities = map(float, fp.readline().strip().split())
	fp.readline()
	Q = []
	for i in range(0,6):
		Q.append(map(float, fp.readline().strip().split()))
	fp.close()
	return ancestral, stabilities, Q

def handle_analysis(filename):
	print "Analysing %s", filename

	ancestral, stabilities, Q = get_results(filename)
	ancestral = zip(ancestral, ORDERS)
	ancestral.sort()
	ancestral.reverse()
	print "Most likely ancestral word order: ", ancestral[0][1]

	stabs2 = stabilities[:]
	stabilities = zip(stabilities, ORDERS)
	stabilities.sort()
	print "Stability ranking: ", ", ".join([x[1] for x in stabilities])

	max = 0
	avg = 0
	for i in range(0,6):
		avg += Q[i][4] + Q[i][5]
		if Q[i][4] > max:
			max = Q[i][4]
		if Q[i][5] > max:
			max = Q[i][5]
	print "Highest probability of transitioning to OXX: ", max
	print "Mean probability of transitioning to OXX: ", avg/12.0

	SOV = Q[0][:]
	SOV = zip(SOV, ORDERS)
	SOV.sort()
	SOV.reverse()
	print "Most probable target of change from SOV: ", SOV[0][1]

	if Q[1][2] > Q[2][1]:
		print "Raw SVO -> VSO > VSO -> SVO"
	else:
		print "Raw VSO -> SVO > SVO -> VSO"

	if stabs2[1]*Q[1][2] > stabs2[2]*Q[2][1]:
		print "Raw SVO -> VSO > VSO -> SVO"
	else:
		print "Raw VSO -> SVO > SVO -> VSO"
	print "----------"

def main():
	for clazz in "family distance feature".split():
		for i in range(0, 5):
			handle_analysis(clazz + str(i) + ".out")

if __name__ == "__main__":
	main()
