from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import argparse
import math
import sys

atoms = []
bindingTreshold = 1.1
maxForce = 0.4
step = 0.01
pairs = []
threes = []
args = {}

def main():
	setupArgs()

	s = ""
	with open("data.txt", "r") as f:
		s = f.read()
	makeAtoms(parseRawInput(s))
	
	result = crunchAtoms()

	print("\nResult:")
	print(result)

	plotAtoms()

def setupArgs():
	global args

	parser = argparse.ArgumentParser(description="Simulate atoms with regard to Lennard-Jones potentials and angular spring potentials")
	group = parser.add_mutually_exclusive_group()
	group.add_argument("-i", "--iterations", help="Manually set the amount of iterations.", type=int, default=1000)
	group.add_argument("-f", "--noiter", help="Don't iterate, only compute forces in current positions.", action="store_true", default=False)
	parser.add_argument("-a", "--noang", help="Don't compute angular forces.", action="store_true", default=False)
	parser.add_argument("-l", "--nolj", help="Don't compute lennard-jones forces.", action="store_true", default=False)
	parser.add_argument("-H", "--nohydrogenforces", help="Ignore forces between hydrogen atoms", action="store_true", default=False)
	args = parser.parse_args()

	if args.nohydrogenforces:
		del well_depth["HH"]
		del best_dist["HH"]

def plotAtoms():
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	# ax.set_aspect('equal', 'datalim')

	ps = [a.pos for a in atoms]
	xs, ys, zs = zip(*ps)
	X = list(xs)
	Y = list(ys)
	Z = list(zs)

	# Create cubic bounding box to simulate equal aspect ratio
	max_range = np.array([max(X)-min(X), max(Y)-min(Y), max(Z)-min(Z)]).max()
	Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(max(X)+min(X))
	Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(max(Y)+min(Y))
	Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(max(Z)+min(Z))
	# Comment or uncomment following both lines to test the fake bounding box:
	for xb, yb, zb in zip(Xb, Yb, Zb):
	   ax.plot([xb], [yb], [zb], 'w')

	ax.scatter(xs, ys, zs, c='r', marker='o')

	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.set_zlabel('Z')

	plt.show()

def crunchAtoms():
	makePairs()
	if args.noiter:
		applyAllForces()
		sumForces()
		limitForces()
		return textAtoms(atoms)
	else:
		print("\nA total of ", len(pairs), "pairs; ", len(threes), " bound triplets")

		for _ in range(args.iterations):
			simulate()

		return textAtoms(atoms)

def applyAllForces():
	if not args.nolj:
		applyLJForces()
	if not args.noang:
		applyAngularForces()


def applyLJForces():
	global pairs
	for pair in pairs:
		lj = LJ_deriv(pair)
		if(lj == 0): continue
		a, b = pair
		aDir = unit_vector(b.pos - a.pos)
		aDir *= lj 
		# aDir *= logify(lj)
		bDir = -aDir

		# print(a, b)
		# print(len(a.forces), len(b.forces))
		a.forces.append(aDir)
		# print(len(a.forces), len(b.forces))
		b.forces.append(bDir)
		# print(len(a.forces), len(b.forces), "\n")

def applyAngularForces():
	print("ng")
	makeThrees()
	for three in threes:
		a, c, b = three
		d = -angle_deriv(three)
		if(d == 0): continue
		#d = logify(d)

		aRPos = unit_vector(a.pos - c.pos)
		bRPos = unit_vector(b.pos - c.pos)

		# interm = unit_vector(aRPos + bRPos)

		# aInterm = interm * vlen(aRPos)
		# bInterm = interm * vlen(bRPos)

		# aDir = unit_vector(aRPos - aInterm)
		# bDir = unit_vector(bRPos - bInterm)

		o = np.cross(bRPos, aRPos)
		aDir = unit_vector(np.cross(o, aRPos))
		bDir = unit_vector(np.cross(bRPos, o))

		aDir *= d
		bDir *= d

		a.forces.append(aDir)
		b.forces.append(bDir)
		c.forces += [-aDir, -bDir]

def simulate():
	applyAllForces()
	sumForces()
	limitForces()
	moveAtoms()

def limitForces():
	global atoms
	largestForce = max([vlen(atom.force) for atom in atoms])
	if largestForce > maxForce:
		limitFactor = largestForce / maxForce
		for atom in atoms:
			atom.force /= limitFactor


def logify(n):
	# Why does adding e not work better? SOLVED
	# Now I know. log(e) = 1, but log(1) = 0.
	# So in adding e, we'd be setting out lowest value at 1
	# But adding 1, we set it at 0.
	# Which is much better.
	return math.log(abs(n)+1)*(n/abs(n)) if n != 0 else 0

def printInfo():
	for atom in atoms:
		print("\n", atom)
		for force in atom.forces:
			print(force, "Mag: ", vlen(force))
		print("Total: ", atom.force, vlen(atom.force))

def sumForces():
	for atom in atoms:
		atom.force = sum(atom.forces)
		atom.forces = []


def moveAtoms():
	for atom in atoms:
		# print(atom.pos, atom.force)
		atom.pos += atom.force# * step

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

def textThrees(ths):
	stringed = [strThree(three) for three in ths]
	#print(*stringed, sep="\n")
	return "\n".join(stringed)

def validTriple(t):
	(a, b), snd = t
	return a in snd or b in snd

def validThree(t):
	_, c, _ = t
	return c.name != "H"

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

def textPairs(ps):
	stringed = [strPair(pair) for pair in ps]
	#print(*stringed, sep="\n")
	return "\n".join(stringed)

def makePairs():
	global pairs
	pairs = pair(atoms)
	pairs = list(filter(validPair, pairs))

def validPair(p):
	if args.nohydrogenforces:
		return pairName(p) != "HH"
	else:
		return True

def makeThrees():
	global threes

	boundPairs = list(filter(isBound_p, pairs))
	possibleBoundTriples = pair(boundPairs)
	boundTriples = list(filter(validTriple, possibleBoundTriples))
	boundThrees = toThrees(boundTriples)
	threes = list(filter(validThree, boundThrees))

def isBound_p(p):
	a, b = p
	return dist(a, b) < bind_treshold[pairName(p)]

def dist(a, b):
	return vlen(b.pos - a.pos)

def parseRawInput(s):
	lines = s.split("\n")
	lines = [line.split() for line in lines]
	lines = [words for words in lines if words and words[0] and words[0][0] != '#']
	print(str(len(lines)) + " lines of input")
	return lines

def makeAtoms(lines):
	global atoms
	atoms = [Atom(words[0], vector3(words[1], words[2], words[3])) for words in lines]
	print(str(len(atoms)) + " Atoms registered")

def textAtoms(atoms):
	result = ["Atoms:"]
	for atom in atoms:
		thisAtom = []
		thisAtom.append(str(atom))
		if atom.forces:
			thisAtom.append("Forces:")
			thisAtom += [" " + str(force) + "\n\t" + str(vlen(force)) for force in atom.forces]
		elif atom.force != None:
			thisAtom.append("Force:\n  " + str(atom.force) + "\n  " + str(vlen(atom.force)))

		result += thisAtom
	result.append("\nPairs:")
	result.append(textPairs(pairs))

	result.append("\nThrees:")
	result.append(textThrees(threes))	
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
	# 	return self.pos/vlen(self.pos)

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
	"HH": 0.3501999999,
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
	"HH": 1.0973955,
	"CH": 1.620681,
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

def printList(os):
	head, *tail = os
	print("[ ", head, end="")
	for o in os:
		print("\n, ", o, end="")
	print(" ]")

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

	if(nm not in three_der_a or nm not in three_der_b):
		return 0;

	a = three_der_a[nm]
	b = three_der_b[nm]
	x = three_angle(t)
	return a*x + b

def vlen(v):
	return np.linalg.norm(v)

# SRC: http://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python
def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / vlen(vector)

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