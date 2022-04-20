import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.fft import fft, ifft
from scipy.integrate import simpson

from allFunc import *

T = 1

dt = 1e-20                                         # in s
t_range = np.arange(0.0, 1e-16+dt, dt)             # in s
# erf param
cutoff = 5e-12                                     # in s
slope = 0.01

domega = 0.1e-3                                    # in eV
w_range = np.arange(1., 3., domega)                # in eV

exc_E = 0.080311203446                             # in Hartree
exc_E *= AUTOJ                                     # in J

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

reorg_E = reorg_E(delta_range, phonons)

# Calculating the dephasing function
(F_t_Re, F_t_Im) = calc_dephasing_F_MC_classical(t_range, delta_range, phonons, T)
print("DONE calc_dephasing_F_MC_classical")
print(t_range)
print(F_t_Re)
# plt.clf()
# plt.plot(t_range, F_t_Re, label="Classical, Re")

f = open("test_Fat1K.dat", "w")
for i in range(len(t_range)):
    f.write("%.20f      %.20f \n" % (t_range[i], F_t_Re[i]))
f.close()

'''
(F_t_Re, F_t_Im) = calc_dephasing_F_MC(t_range, delta_range, phonons, T)
print("DONE calc_dephasing_F_MC_QM")
plt.plot(t_range, F_t_Re, label="QM, MC, Re")
plt.plot(t_range, F_t_Im, ":", label="QM, MC, Im")

(F_t_Re, F_t_Im) = calc_dephasing_F_E17(t_range, delta_range, phonons, T)
print("DONE calc_dephasing_F_E17_QM")
plt.plot(t_range, F_t_Re, label="QM, E17, Re")
plt.plot(t_range, F_t_Im, ":", label="QM, E17, Im")

plt.grid()
plt.title("Dephasing function at 1K")
plt.xlabel("Time [s]")
plt.ylabel("Dephasing function [unitless]")
plt.legend()
plt.savefig("test.png")
'''
