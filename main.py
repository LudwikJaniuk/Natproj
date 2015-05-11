import numpy as np

def main():
	with open("data.txt", "r") as f:
		s = f.read()
		lines = parseRawInput(s)
		atoms = makeAtoms(lines)
		result = textAtoms(atoms)
		print(result)

def parseRawInput(s):
	lines = s.split("\n")
	lines = [line.split() for line in lines]
	lines = [words for words in lines if words[0][0] != '#']
	return lines

def makeAtoms(lines):
	return [Atom(words[0], vector3(words[1], words[2], words[3])) for words in lines]

def textAtoms(atoms):
	result = [str(atom) for atom in atoms]
	return "\n".join(result)

def vector3(x, y ,z):
	return np.array([float(x), float(y), float(z)])

class Atom:
	""" An atom, with its name and position """
	atoms = "CHO"
	name = "X"
	pos = np.array([0, 0, 0])

	def __init__(self, _name, _pos):
		_name = _name.upper()
		assert(_name in Atom.atoms)
		assert(len(_pos) == 3)
		
		self.name = _name.upper()
		self.pos = _pos

	def __str__(self):
		return self.name + " " + str(self.pos);







main()