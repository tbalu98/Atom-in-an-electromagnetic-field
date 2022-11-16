import math
import cmath
import numpy as np
from scipy.integrate import quad

au_to_ev = 27.211324570273

#Different types of envelope functions to simulate a laser impulse
class Sin_sq_envelope:
    def __init__(self, DT, T):
        self.DT = DT
        self.T = T


    def evelope(self, t):
        if t<0 or t>self.T:
            return 0.0
        if 0<=t<=self.DT:
            return math.sin(math.pi * t /(2*self.DT))**2
        if self.DT <= t <= self.T-self.DT:
            return 1.0
        if self.T-self.DT<= t <= self.T:
            return math.sin(math.pi * t /(2*self.DT))**2
    
class Gaussian_envelope:
    def __init__(self, FWHM_fs):
        self.set_tau_from_FWHM_fs(FWHM_fs)

    def set_tau_from_FWHM_fs(self,FWHM_fs):
        FWHM = FWHM_fs*41.34137333656
        self.tau = FWHM/(2*math.sqrt(math.log(2)))

    def envelope(self, t):
        return gaussian_function(t, self.tau)

    def envelope_sq(self, time):
        return self.envelope(time)**2

# A function that cannot be expressed analitically
class numeric_fuction:
    def __init__(self, xs,ys):
        self.xs = xs
        self.dx = xs[1]-xs[0]
        self.lower = xs[0]
        self.higher = xs[-1]
        self.ys = ys

    def eval(self, x):
        if self.lower<=x<=self.higher:
            index = int((x-self.lower) / self.dx)
            #print(index)
            return self.ys[index]
        else:
            return None

#A few function for miscellous purposes

def wavelength_to_energy_au(wavelength):
    
    return 1239.8/(wavelength * au_to_ev) 


def cabssq(z):
    return math.pow(abs(z),2)

def gaussian_function(t, tau):
    return math.exp(-(pow(t,2)/(2*pow(tau,2))))


def cnst_envelope(t):
    s = t
    return 1

def abssq(x):
    return abs(x)**2

def simple_integrator(dx, ys):
    res = 0
    for i in range(1,len(ys)):
        trapez = (ys[i] + ys[i-1])/2
        res += trapez

    return res*dx

def definite_integral(dx, ys):
    res = 0
    vals = [0]

    for i in range(1,len(ys)):
        trapez = (ys[i] + ys[i-1])/2
        res += trapez*dx
        vals.append(res)
    return vals


def pi_pulse_field_strength( int_from, int_to, rabi_freq, envelope_fun):
    
    int_of_envelope_fun_sq, dev = quad(lambda t: cabssq( envelope_fun(t)),int_from, int_to)
    E0 = math.sqrt(math.pi/(2*(abs(rabi_freq)*int_of_envelope_fun_sq)))
    return E0        
