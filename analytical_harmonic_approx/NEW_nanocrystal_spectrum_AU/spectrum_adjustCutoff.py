import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.fft import fft, ifft
from scipy.integrate import simpson

from allFunc import *

T = 1                            # unit: K 
BETA = 1./ (KB * T)   

dt = 1e-15
t_range = np.arange(0.0, 30e-12+dt, dt)
domega = 0.1e-3                  # in eV
w_range = np.arange(1., 3., domega)

# erf param
cutoff_1 = 5e-12                 # in s
cutoff_2 = 15e-12                # in s
cutoff_3 = 25e-12                # in s
slope = 0.01

excitonic_e = 0.080311203446 * AUTOJ

# Reading frequency and coupling matrix element data
nExc = 20
w = read_w("w.dat")
coupling = read_Vklq("Vklq-diabatic.dat", len(w), nExc)
w = w[6:]
V = coupling[6:, 0]
delta_range = V / w**2
print("DONE reading")

# Calculating the correlation function
(F_t_Re, F_t_Im) = calc_dephasing_F_E17(t_range, delta_range, w, T)
print("DONE calc_dephasing_F_E17")

'''
(DC_t_Re, DC_t_Im) = calc_delta_corr_E17(t_range, delta_range, w, T)
DC_t_Re *= JTOEV**2
DC_t_Im *= JTOEV**2
print("DONE calc_delta_corr_E17")
'''

# erf correction to the dephasing functions
F_t_Re_erf_1 = multiply_erf(cutoff_1, slope, F_t_Re-np.mean(F_t_Re[-1000:-1]), dt, t_range)
F_t_Im_erf_1 = multiply_erf(cutoff_1, slope, F_t_Im-np.mean(F_t_Im[-1000:-1]), dt, t_range)
F_t_Re_erf_2 = multiply_erf(cutoff_2, slope, F_t_Re-np.mean(F_t_Re[-1000:-1]), dt, t_range)
F_t_Im_erf_2 = multiply_erf(cutoff_2, slope, F_t_Im-np.mean(F_t_Im[-1000:-1]), dt, t_range)
F_t_Re_erf_3 = multiply_erf(cutoff_3, slope, F_t_Re-np.mean(F_t_Re[-1000:-1]), dt, t_range)
F_t_Im_erf_3 = multiply_erf(cutoff_3, slope, F_t_Im-np.mean(F_t_Im[-1000:-1]), dt, t_range)
print("DONE erf_multiplication")

# FT the dephasing functions to obtain the spectrum from the correlation function
delta_E_E17 = excitonic_e + reorg_E(delta_range, w)    # in J
(spectrum_w_cutoff1, spectrum_I_cutoff1) = spectrum(delta_E_E17, t_range, F_t_Re_erf_1, F_t_Im_erf_1, w_range)
(spectrum_w_cutoff2, spectrum_I_cutoff2) = spectrum(delta_E_E17, t_range, F_t_Re_erf_2, F_t_Im_erf_2, w_range)
(spectrum_w_cutoff3, spectrum_I_cutoff3) = spectrum(delta_E_E17, t_range, F_t_Re_erf_3, F_t_Im_erf_3, w_range)
print("DONE spectrum_E17")
write_spectrum("adjustCutoff_spectrum_E17_cutoff1_T=1.dat", spectrum_w_cutoff1, spectrum_I_cutoff1)
write_spectrum("adjustCutoff_spectrum_E17_cutoff2_T=1.dat", spectrum_w_cutoff2, spectrum_I_cutoff2)
write_spectrum("adjustCutoff_spectrum_E17_cutoff3_T=1.dat", spectrum_w_cutoff3, spectrum_I_cutoff3)

#FWHM
k = open("adjustCutoff_FWHMs.dat", "w")
k.write("# Cutoff Time (ps)       FWHM (eV)")
(linewidth1, l_index_1, r_index_1) = calc_FWHM(spectrum_w_cutoff1, spectrum_I_cutoff1)
k.write("5ps      %f"%linewidth1)
(linewidth2, l_index_2, r_index_2) = calc_FWHM(spectrum_w_cutoff2, spectrum_I_cutoff2)
k.write("15ps      %f"%linewidth2)
(linewidth3, l_index_3, r_index_3) = calc_FWHM(spectrum_w_cutoff3, spectrum_I_cutoff3)
k.write("25ps      %f"%linewidth3)

# plotting
fig, axs = plt.subplots(3,2, figsize=(10,10))
fig.suptitle(r"The dephasing function and spectra at different cutoff points, E17 model, 3nmCdSe at T=1")

for i in range(3): 
    for j in range(2): 
        axs[i,j].grid()

axs[0,0].set(xlabel="Time (s)", ylabel=r"$\langle F(t) \rangle$ (unitless)", title="E17 dephasing function cutoff at 5 ps", yscale="log")
axs[0,0].plot(t_range, F_t_Re, label=r"$\langle F(t) \rangle _{Re} $")
axs[0,0].plot(t_range, F_t_Im, ":", label=r"$\langle F(t) \rangle _{Im}$")
axs[0,0].plot(t_range, F_t_Re_erf_1, label=r"$\langle F(t) \rangle _{Re, shifted} \cdot erf$")
axs[0,0].plot(t_range, F_t_Im_erf_1, ":", label=r"$\langle F(t) \rangle _{Im, shifted} \cdot erf$")

axs[1,0].set(xlabel="Time (s)", ylabel=r"$\langle F(t) \rangle$ (unitless)", title="E17 dephasing function cutoff at 15 ps", yscale="log")
axs[1,0].plot(t_range, F_t_Re, label=r"$\langle F(t) \rangle _{Re} $")
axs[1,0].plot(t_range, F_t_Im, ":", label=r"$\langle F(t) \rangle _{Im}$")
axs[1,0].plot(t_range, F_t_Re_erf_2, label=r"$\langle F(t) \rangle _{Re, shifted} \cdot erf$")
axs[1,0].plot(t_range, F_t_Im_erf_2, ":", label=r"$\langle F(t) \rangle _{Im, shifted} \cdot erf$")

axs[2,0].set(xlabel="Time (s)", ylabel=r"$\langle F(t) \rangle$ (unitless)", title="E17 dephasing function cutoff at 25 ps", yscale="log")
axs[2,0].plot(t_range, F_t_Re, label=r"$\langle F(t) \rangle _{Re} $")
axs[2,0].plot(t_range, F_t_Im, ":", label=r"$\langle F(t) \rangle _{Im}$")
axs[2,0].plot(t_range, F_t_Re_erf_3, label=r"$\langle F(t) \rangle _{Re, shifted} \cdot erf$")
axs[2,0].plot(t_range, F_t_Im_erf_3, ":", label=r"$\langle F(t) \rangle _{Im, shifted} \cdot erf$")

# spectra
axs[0,1].set(xlabel="Energey [meV]", title="Spectra", xlim=(2100, 2250))
axs[0,1].plot(spectrum_w_cutoff1*1000, spectrum_I_cutoff1, linewidth=3, label="erf cutoff at 5ps")
axs[0,1].plot(spectrum_w_cutoff2*1000, spectrum_I_cutoff2, linewidth=2, label="erf cutoff at 15ps")
axs[0,1].plot(spectrum_w_cutoff3*1000, spectrum_I_cutoff3, linewidth=1, label="erf cutoff at 25ps")

axs[1,1].set(xlabel="Energey [meV]", title="Spectra Diff", xlim=(2100, 2250))
axs[1,1].plot(spectrum_w_cutoff1*1000, spectrum_I_cutoff2-spectrum_I_cutoff1, linewidth=3, label="15ps spectrum - 5ps spectrum")
axs[1,1].plot(spectrum_w_cutoff1*1000, spectrum_I_cutoff3-spectrum_I_cutoff1, linewidth=1, label="25ps spectrum - 5ps spectrum")

for i in range(3): 
    for j in range(2): 
        axs[i,j].legend()
fig.tight_layout()
fig.savefig("adjustCutoff.png")
