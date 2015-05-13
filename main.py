import numpy as np

atoms = []
bindingTreshold = 1.1



def main():
	with open("data.txt", "r") as f:
		s = f.read()
		lines = parseRawInput(s)
		print("")
		global atoms
		atoms = makeAtoms(lines)
		crunchAtoms()
		result = textAtoms(atoms)
		print("\nResult:")
		print(result)

def crunchAtoms():
	pairs, threes = atomGroups()
	print("\nA total of ", len(pairs), "pairs; ", len(threes), " bound triplets")

	#WIP
	# for pair in pairs:
	# 	lj = lj_deriv(pair)
	# 	a, b = pair
	# 	aDir = b - a
	# 	bDir = -aDir




def atomGroups():
	def pairs(l):
		if(len(l) > 1):
			head, *tail = l
			firstPairs = [(head, other) for other in tail]
			return firstPairs + pairs(tail)
		else:
			return []

	def strPair(p):
		fst, snd = p
		return str(fst.name) + "=="  + str(snd.name) + " Dist: " + str(dist(fst, snd))

	def strThree(th):
		a, c, b = th
		return a.name + "--" + c.name + "--" + b.name

	def printPairs(ps):
		stringed = [strPair(pair) for pair in ps]
		print(*stringed, sep="\n")

	def printThrees(ths):
		stringed = [strThree(three) for three in ths]
		print(*stringed, sep="\n")

	def validTriple(t):
		(a, b), snd = t
		return a in snd or b in snd

	def three(t):
		(a, b), (c, d) = t
		if a == c:
			return (b, a, d)
		if a == d:
			return (b, a, c)
		if b == c:
			return (a, b, d)
		if b == d:
			return (a, b, c)

	def threes(ts):
		return [three(t) for t in ts]

	allPairs = pairs(atoms)

	print("\nPairs:")
	printPairs(allPairs)

	boundPairs = list(filter(isBound_p, allPairs))
	possibleBoundTriples = pairs(boundPairs)
	boundTriples = list(filter(validTriple, possibleBoundTriples))
	boundThrees = threes(boundTriples)

	print("\nBound Threes: ")
	printThrees(boundThrees)

	return allPairs, boundThrees

def isBound_p(p):
	return isBound(*p)

def isBound(a, b):
	return dist(a, b) < bindingTreshold

def dist(a, b):
	return np.linalg.norm(b.pos - a.pos)

def parseRawInput(s):
	lines = s.split("\n")
	lines = [line.split() for line in lines]
	lines = [words for words in lines if words[0][0] != '#']
	print(str(len(lines)) + " lines of input")
	return lines

def makeAtoms(lines):
	result = [Atom(words[0], vector3(words[1], words[2], words[3])) for words in lines]
	print(str(len(result)) + " Atoms registered")
	return result

def textAtoms(atoms):
	result = [str(atom) for atom in atoms]
	return "\n".join(result)

def lennard_jones(a, b):
	r = dist(a, b)
	rOpt, E = pairProps(a, b)

def vector3(x, y ,z):
	return np.array([float(x), float(y), float(z)])

class Atom:
	""" An atom, with its name and position """
	atoms = "CHO"
	names = { 
		"C": "Coal", 
		"H": "Hydrogen", 
		"O": "Oxygen", 
		"X": "Invalid" }
	name = "X"
	pos = np.array([0, 0, 0])
	forces = []
	def clearForces(self):
		self.forces = []

	def norm(self):
		return self.pos/np.linalg.norm(self.pos)

	def __init__(self, _name, _pos):
		_name = _name.upper()
		assert(_name in Atom.atoms)
		assert(len(_pos) == 3)
		
		self.name = _name.upper()
		self.pos = _pos

		print("Created Atom: " + Atom.names[self.name] + " at pos: " + str(self.pos))

	def __str__(self):
		return self.name + " " + str(self.pos);

	# Energi, i hartree
well_depth = {
	"CH": 0.1332,
	"CO": 0.0921590301,
	"CC": 0.1433381031,
	"HO": 0.1107685941
}

# Avstånd, i ångström
best_dist = {
	"HH": 0.731597,
	"CH": 1.080454,
	"CC": 1.513910667,
	"CO": 1.428491,
	"HO": 0.946065
}

def LJ_potential(a, r):
	rm = best_dist[a]
	E = well_depth[a]
	return E*((rm/r)**12 - 2*(rm/r)**6)

def LJ_deriv(a, r):
	rm = best_dist[a]
	E = well_depth[a]
	return (12 * E * rm**6 * (r**6 - rm**6)) / r**13


main()