import numpy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import RK4_algorithm
import Em_field
from Atom import Atom
from Em_field_atom_interaction import Em_field_atom_interaction
import RK4_algorithm
from Sparse_matrix import Trip_in_sparse_matrix
import winsound
import math_utilities as math_utils
res_folder = 'Na_results/'




#E0 = 0.005

E0 = 0.002025638142579677

#properties of Na atom:
data_folder = 'Na_data/'

dipolematrix = pd.read_csv(data_folder+'dipolemat67_df.csv', index_col = 'row_index')
eiendiffmat = pd.read_csv(data_folder+'eiendiffmat67_df.csv', index_col = 'row_index')
states_labels = ['3sb','4sb','5sb','6sb','7sb','8sb','9sb','10sb','11sb','12sb','13sb','14sb','15sb','16sb',
                '3pb','4pb','5pb','6pb','7pb','8pb','9pb','10pb','11pb','12pb','13pb','14pb','15pb','16pb',
                '3db','4db','5db','6db','7db','8db','9db','10db','11db','12db','13db','14db','15db','16db',
                '4fb','5fb','6fb','7fb','8fb','9fb','10fb','11fb','12fb','13fb','14fb','15fb','16fb',
                '5gb','6gb','7gb','8gb','9gb','10gb','11gb','12gb','13gb','14gb','15gb','16gb',                
                ]

trip = Trip_in_sparse_matrix(dipolematrix.to_numpy(dtype = np.complex64))


#Laser parameters

half_res_omega = 0.1167599709407699943*0.5

frequency = half_res_omega*1.0

print(frequency)



# The electromagnetic field which is a laser impulse
E0 = 0.002


FWHM_fs = 12.5

# The envelope function
gaussian_envelope = math_utils.Gaussian_envelope(FWHM_fs)

#The phase fun
temp_phase_fun = lambda x: 0 # This electromagnetic field od not frequency modulated
em_field = Em_field.Em_field(E0, frequency, gaussian_envelope.envelope , temp_phase_fun )

atom = Atom(dipolematrix, eiendiffmat)

interaction = Em_field_atom_interaction(atom, em_field)

print(interaction.Stark_shift('4sb'))
print(interaction.Rabi_frequency( '3sb','4sb'))

#The simulation:
iterations = 5000

start = -FWHM_fs*100
stop = FWHM_fs*100

# The initial state of the atom
eq_num = 67
c = numpy.zeros((67,1), dtype="complex")
c.itemset(0,complex(1,0))
comm_const = complex(1,0)


# The execution of the Runge-Kutta algorithm
rk4o = RK4_algorithm.RK4_algorithm(start,stop,iterations,comm_const, interaction.hamiltonian, c, eq_num, trip,'times', states_labels)

c_out =  rk4o.calculate()
c_out = pd.DataFrame(c_out)
c_out = c_out.set_index('times')
c_out.to_csv('Na_result_'+str(E0)+'_'+str(round(frequency,4))+'.csv')


plt.plot(c_out.index, list(map(math_utils.cabssq,c_out['3s'] )) )
plt.plot(c_out.index, list(map(math_utils.cabssq,c_out['4s'] )) )
plt.plot(c_out.index, list(map(math_utils.cabssq,c_out['4d'] )) )

plt.show()
winsound.Beep(440,10000)