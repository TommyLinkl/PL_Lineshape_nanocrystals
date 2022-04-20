import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.fft import fft, ifft
from scipy.integrate import simpson

from allFunc import *

T = 1

dt = 50   #1e-15                                         # in a.u.
t_range = np.arange(0.0, 4.134e5+dt, dt)             # in a.u.
# erf param
'''
cutoff = 5e-12                                     # in s
slope = 0.01
'''

domega = 0.1e-3                                    # in eV
w_range = np.arange(1., 3., domega)                # in eV

exc_E = 0.080311203446                             # in Hartree
# exc_E *= AUTOJ                                     # in J

phonon_filename = "w.dat"
coupling_filename = "Vklq-diabatic.dat"

BETA = 1./ (KB * T)

# Reading frequency and coupling matrix element data
nExc = 20
phonons = read_w(phonon_filename)
coupling = read_Vklq(coupling_filename, len(phonons), nExc)
phonons = phonons[6:]
V = coupling[6:, 0]
delta_range = (V) / phonons**2
print(delta_range)
print("DONE reading")

# Calculating the dephasing function
(F_t_Re, F_t_Im) = calc_dephasing_F_MC_classical(t_range, delta_range, phonons, T)
print("DONE calc_dephasing_F_MC_classical")
print(t_range)
print(F_t_Re)
plt.clf()
plt.plot(t_range*AUTOS, F_t_Re, label="Classical, Re")

f = open("AU_F_t.dat", "w")
for i in range(len(t_range)):
    f.write("%.20f      %.20f      %.20f\n" % (t_range[i]*AUTOS, F_t_Re[i], F_t_Im[i]))
f.close()

'''
(F_t_Re, F_t_Im) = calc_dephasing_F_MC(t_range, delta_range, phonons, T)
print("DONE calc_dephasing_F_MC_QM")
plt.plot(t_range, F_t_Re, label="QM, Re")
plt.plot(t_range, F_t_Im, label="QM, Im")
'''

plt.grid()
plt.xlabel("Time (s)")
plt.ylabel("Classical dephasing function (unitless)")
plt.legend()
plt.savefig("test_au.png")
