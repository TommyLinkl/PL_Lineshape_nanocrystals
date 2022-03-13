import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.fft import fft, ifft
from scipy.integrate import simpson

from allFunc import *

global T
T = 1                            # unit: K 
global BETA
BETA = 1./ (KB * T)   

dt = 1e-15                                         # in s
t_range = np.arange(0.0, 8e-12+dt, dt)             # in s
# erf param
cutoff = 5e-12                                     # in s
slope = 0.01

domega = 0.1e-3                                    # in eV
w_range = np.arange(1., 3., domega)                # in eV

exc_E = 0.080311203446                             # in Hartree
exc_E *= AUTOJ                                     # in J

phonon_filename = "w.dat"
coupling_filename = "Vklq-diabatic.dat"

# Reading frequency and coupling matrix element data
nExc = 20
w = read_w(phonon_filename)
coupling = read_Vklq(coupling_filename, len(w), nExc)
w = w[6:]
V = coupling[6:, 0]
delta_range = (V) / w**2
print("DONE reading")

# Calculating the correlation function, and writing to file
(F_t_Re_E17, F_t_Im_E17) = calc_dephasing_F_E17(t_range, delta_range, w)
print("DONE calc_dephasing_F_E17")

(F_t_Re_MC, F_t_Im_MC) = calc_dephasing_F_MC(t_range, delta_range, w)
print("DONE calc_dephasing_F_MC")

(DC_t_Re_E17, DC_t_Im_E17) = calc_delta_corr_E17(t_range, delta_range, w)
DC_t_Re_E17 *= JTOEV**2
DC_t_Im_E17 *= JTOEV**2
print("DONE calc_delta_corr_E17")

(DC_t_Re_MC, DC_t_Im_MC) = calc_delta_corr_MC(t_range, delta_range, w)
DC_t_Re_MC *= JTOEV**2
DC_t_Im_MC *= JTOEV**2
print("DONE calc_delta_corr_MC")
write_dephasing_F("zeroCoupling_F.dat", t_range, F_t_Re_MC, F_t_Im_MC, F_t_Re_E17, F_t_Im_E17)
write_delta_corr("zeroCoupling_Delta_Corr.dat", t_range, DC_t_Re_MC, DC_t_Im_MC, DC_t_Re_E17, DC_t_Im_E17)

# erf correction to the dephasing functions
F_t_Re_E17_erf = multiply_erf(cutoff, slope, F_t_Re_E17, dt)
F_t_Im_E17_erf = multiply_erf(cutoff, slope, F_t_Im_E17, dt)
F_t_Re_MC_erf = multiply_erf(cutoff, slope, F_t_Re_MC, dt)
F_t_Im_MC_erf = multiply_erf(cutoff, slope, F_t_Im_MC, dt)
print("DONE erf_multiplication")

# FT the dephasing functions to obtain the spectrum from the correlation function
delta_E_E17 = excitonic_e + reorg_E(delta_range, w)    # in J
delta_E_MC = excitonic_e                               # in J
(spectrum_w_E17, spectrum_I_E17) = spectrum(delta_E_E17, t_range, F_t_Re_E17_erf, F_t_Im_E17_erf, w_range)
print("DONE spectrum_E17")
(spectrum_w_MC, spectrum_I_MC) = spectrum(delta_E_MC, t_range, F_t_Re_MC_erf, F_t_Im_MC_erf, w_range)
print("DONE spectrum_MC")
write_spectrum("zeroCoupling_spectrum_E17.dat", spectrum_w_E17, spectrum_I_E17)
write_spectrum("zeroCoupling_spectrum_MC.dat", spectrum_w_MC, spectrum_I_MC)

# plotting
fig, axs = plt.subplots(2,2, figsize=(8,8))
fig.suptitle(r"The dephasing function and spectra, MC and E17 models, 3nmCdSe, T=1K, zero sys-bath coupling")

for i in range(2): 
    for j in range(2): 
        axs[i,j].grid()
axs[0,0].set(xlabel="Time (s)", ylabel=r"$\langle F(t) \rangle$ (unitless)", title="MC dephasing function")
axs[0,0].plot(t_range, F_t_Re_MC, label=r"$\langle F(t) \rangle _{Re} $")
axs[0,0].plot(t_range, F_t_Im_MC, ":", label=r"$\langle F(t) \rangle _{Im}$")
axs[0,0].plot(t_range, F_t_Re_MC_erf, label=r"$\langle F(t) \rangle _{Re} \cdot erf$")
axs[0,0].plot(t_range, F_t_Im_MC_erf, ":", label=r"$\langle F(t) \rangle _{Im} \cdot erf$")

axs[1,0].set(xlabel="Time (s)", ylabel=r"$\langle F(t) \rangle$ (unitless)", title="E17 dephasing function")
axs[1,0].plot(t_range, F_t_Re_E17, label=r"$\langle F(t) \rangle _{Re} $")
axs[1,0].plot(t_range, F_t_Im_E17, ":", label=r"$\langle F(t) \rangle _{Im}$")
axs[1,0].plot(t_range, F_t_Re_E17_erf, label=r"$\langle F(t) \rangle _{Re} \cdot erf$")
axs[1,0].plot(t_range, F_t_Im_E17_erf, ":", label=r"$\langle F(t) \rangle _{Im} \cdot erf$")

axs[0,1].set(xlabel="Energey [meV]", title="Spectra", xlim=(2175, 2225))
axs[0,1].plot(spectrum_w_MC*1000, spectrum_I_MC, label="MC spectrum (w/ erf)")
axs[0,1].plot(spectrum_w_E17*1000, spectrum_I_E17, label="E17 spectrum (w/ erf)")

axs[1,1].set(xlabel="Energey [meV]", title="Spectra Diff", xlim=(2175, 2225))
axs[1,1].plot(spectrum_w_E17*1000, spectrum_I_MC-spectrum_I_E17, label="MC spectrum - E17 spectrum")


for i in range(2): 
    for j in range(2): 
        axs[i,j].legend()
fig.tight_layout()
fig.savefig("nanocrystal_zeroCoupling.png")
