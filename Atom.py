import numpy as np

class Atom:
    def __init__(self, dipoles, eiendiffs):
        self.dipoles = dipoles # Transition dipole moment matrix of the atom
        self.eiendiffs = eiendiffs # The matrix of the atom's valence electron's energy levels' differences

