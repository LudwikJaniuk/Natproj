from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import math

iterations = 1000
atoms = []
bindingTreshold = 1.1
step = 0.01
pairs = []
threes = []


def main():
	s = ""
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

	print("\nPairs:")
	printPairs(pairs)

	plotAtoms()

def plotAtoms():
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')

	ps = [a.pos for a in atoms]
	xs, ys, zs = zip(*ps)

	ax.scatter(xs, ys, zs, c='r', marker='o')

	# ax.set_xlabel('X Label')
	# ax.set_ylabel('Y Label')
	# ax.set_zlabel('Z Label')

	plt.show()

def crunchAtoms():
	atomGroups()
	print("\nA total of ", len(pairs), "pairs; ", len(threes), " bound triplets")

	for _ in range(iterations):
		simulate()

def simulate():
	clearForces()
	applyLJForces()
	makeThrees()
	applyAngularForces()
	sumForces()
	moveAtoms()
	# printPairs(pairs)
	# print("")
	# printInfo()

def logify(n):
	# Why does adding e not work better? SOLVED
	# Now I know. log(e) = 1, but log(1) = 0.
	# So in adding e, we'd be setting out lowest value at 1
	# But adding 1, we set it at 0.
	# Which is much better.
	return math.log(abs(n)+1)*(n/abs(n))

def printInfo():
	for atom in atoms:
		print("\n", atom)
		for force in atom.forces:
			print(force, "Mag: ", np.linalg.norm(force))
		print("Total: ", atom.force, np.linalg.norm(atom.force))

def clearForces():
	for atom in atoms:
		atom.forces = []

def sumForces():
	for atom in atoms:
		atom.force = sum(atom.forces)

def applyLJForces():
	global pairs
	for pair in pairs:
		lj = LJ_deriv(pair)
		a, b = pair
		aDir = b.pos - a.pos
		# A direct approach can cause crazy jumps
		# aDir *= lj 
		aDir *= logify(lj)
		bDir = -aDir

		# print(a, b)
		# print(len(a.forces), len(b.forces))
		a.forces.append(aDir)
		# print(len(a.forces), len(b.forces))
		b.forces.append(bDir)
		# print(len(a.forces), len(b.forces), "\n")

def applyAngularForces():
	for three in threes:
		a, _, b = three
		d = angle_deriv(three)
		d = logify(d)

		aDir = b.pos - a.pos
		aDir *= logify(d)
		bDir = -aDir

		a.forces.append(aDir)
		b.forces.append(bDir)



def moveAtoms():
	for atom in atoms:
		print(atom.pos, atom.force)
		atom.pos += atom.force * step

def pair(l):
	if(len(l) > 1):
		head, *tail = l
		firstPairs = [(head, other) for other in tail]
		return firstPairs + pair(tail)
	else:
		return []

def strPair(p):
	fst, snd = p
	return str(fst.name) + "=="  + str(snd.name) + " Dist: " + str(dist(fst, snd))

def strThree(th):
	a, c, b = th
	return a.name + "--" + c.name + "--" + b.name

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

def toThrees(ts):
	return [three(t) for t in ts]

def printPairs(ps):
	stringed = [strPair(pair) for pair in ps]
	print(*stringed, sep="\n")

def atomGroups():
	global pairs

	pairs = pair(atoms)
	makeThrees()

def makeThrees():
	global threes

	boundPairs = list(filter(isBound_p, pairs))
	possibleBoundTriples = pair(boundPairs)
	boundTriples = list(filter(validTriple, possibleBoundTriples))
	boundThrees = toThrees(boundTriples)

	print("\nBound Threes: ")
	printThrees(boundThrees)

	threes = boundThrees

def isBound_p(p):
	a, b = p
	return dist(a, b) < bind_treshold[pairName(p)]

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
	result = [str(atom) + " " + str(np.linalg.norm(atom.force)) for atom in atoms]
	return "\n".join(result)

# def lennard_jones(a, b):
# 	r = dist(a, b)
# 	rOpt, E = pairProps(a, b)


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
	force = None
	def clearForces(self):
		self.forces = []

	# def norm(self):
	# 	return self.pos/np.linalg.norm(self.pos)

	def __init__(self, _name, _pos):
		_name = _name.upper()
		assert(_name in Atom.atoms)
		assert(len(_pos) == 3)

		self.name = _name.upper()
		self.pos = _pos
		self.forces = []

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

# Avstånd, i ångström
bind_treshold = {
	"HH": -1.0973955,
	"CH": -1.620681,
	"CC": 2.270866,
	"CO": -2.1427365,
	"HO": -1.4190975
}

three_der_a = {
	"HCH": 0.000070322,
	"CCC": 0.000111,
	"HOH": 0.000073,
	"OCO": 0.00005400,
	"COH": 0.000096572,
	"HCO": 0.000096572
}

three_der_b = {
	"HCH": -0.0077106,
	"CCC": -0.012611,
	"HOH": -0.008059,
	"OCO": -0.009721,
	"COH": -0.011591,	
	"HCO": -0.011591
}

# UTIL FUNCTIONS

def vector3(x, y ,z):
	return np.array([float(x), float(y), float(z)])

def pairName(p):
	a, b = p
	nm = a.name + b.name
	result = "".join(sorted(nm))
	return result

def threeName(t):
	a, c, b = t
	sides = pairName((a, b))
	return sides[0] + c.name + sides[1]

def LJ_potential(p):
	a = pairName(p)

	if(a not in best_dist or a not in well_depth):
		return 0;

	rm = best_dist[a]
	E = well_depth[a]
	r = dist(*p)

	return E*((rm/r)**12 - 2*(rm/r)**6)

def LJ_deriv(p):
	a = pairName(p)

	if(a not in best_dist or a not in well_depth):
		return 0;

	rm = best_dist[a]
	E = well_depth[a]
	r = dist(*p)
	
	return (12 * E * rm**6 * (r**6 - rm**6)) / r**13

def three_angle(t):
	a, c, b = t
	a = a.pos
	b = b.pos
	c = c.pos
	v1 = a - c
	v2 = b - c
	return math.degrees(angle_between(v1, v2))

def angle_deriv(t):
	nm = threeName(t)
	a = three_der_a[nm]
	b = three_der_b[nm]
	x = three_angle(t)
	return a*x + b


# SRC: http://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python
def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

# SRC: http://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python
def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    angle = np.arccos(np.dot(v1_u, v2_u))
    if np.isnan(angle):
        if (v1_u == v2_u).all():
            return 0.0
        else:
            return np.pi
    return angle

main()