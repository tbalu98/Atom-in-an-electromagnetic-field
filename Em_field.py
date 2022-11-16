from math import cos

# A class for creating custom electromagnetic field. The user can specify the envelope, frequency and frequecy modulation. 
class Em_field:

    def __init__(self,  E0, omega, envelope_fun, phase_fun):
        
        self.phase_fun = phase_fun
        self.E0 = E0
        self.omega = omega
        self.envelope_fun = envelope_fun
        self.phase_fun = phase_fun
    

    def efield(self,time):
        return self.E0*cos( self.omega * time -self.phase_fun(time)*0.5 )*self.envelope_fun(time)





