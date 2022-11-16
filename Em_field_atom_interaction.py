import cmath
import math
import math_utilities
import math_utilities
import cmath
import numpy as np
from scipy.integrate import quad
import math_utilities as math_utils


class Em_field_atom_interaction():
    def __init__(self, atom, em_field):
        self.atom = atom
        self.em_field = em_field

    #The hamiltonian operator which is a time dependent matrix

    def hamiltonian(self, time, i, j):
        omega_ji = float(self.atom.eiendiffs.iloc[j, i])
        mu_ij = complex(self.atom.dipoles.iloc[i,j])
        return complex(0,1)*cmath.exp(complex(0,-time * omega_ji))*mu_ij*self.em_field.efield(time)


    #The peak value of the dinamic Stark shift between the 3s and the 4s states
    def DSS_4s_3s_calc(self):
        self.DSS4s_3s = (self.starkShift(1)-self.starkShift(0)) *self.em_field.E0**2

    def DSS_4s_3s(self, time):
        return self.DSS4s_3s*math_utilities.gaussian_function(time,self.em_field.tau)
    
    # A state's energy level's Stark shift
    # g is the ground state's label
    # e is the excited state's label
    # m is every other states' label

    def Stark_shift(self, k):
        res = 0.0
        state_labels = self.atom.dipoles.index.to_list()
        for m in state_labels:
            mu_km = complex(self.atom.dipoles.loc[k][m])
            omega_mk = complex(self.atom.eiendiffs.loc[m][k])
            omega_0 = self.em_field.omega

            res += math_utilities.cabssq(mu_km) * omega_mk / (omega_mk**2 - omega_0**2)
        return res*0.5

    # The Rabi osciallation between two states
    # g is the ground state's label
    # e is the excited state's label
    # m is every other states' label
    def Rabi_frequency(self, g, e):
        res = 0.0
        states = self.atom.dipoles.index.to_list()
        for m in states:

            mu_em = complex(self.atom.dipoles.loc[e][ m])
            mu_mg = complex(self.atom.dipoles.loc[ m] [g])
            omega_mg = float(self.atom.eiendiffs.loc[m][g])
            omega_0 = self.em_field.omega
            res +=-0.25*mu_em * mu_mg /( omega_mg - omega_0 )
        return res


    # The specific field strenght that with the rigth phase function (C. Trallero-Herrero 2005)
    def pi_pulse_field_strength(self, int_from, int_to):
    
        int_of_envelope_fun_sq, dev = quad(lambda t: math_utils.cabssq( self.em_field.envelope_fun(t)),int_from, int_to)
        rabi_freq = self.Rabi_frequency( 0,1 )
        
        E0 = math.sqrt(math.pi/(2*(abs(rabi_freq)*int_of_envelope_fun_sq)))
        return E0 