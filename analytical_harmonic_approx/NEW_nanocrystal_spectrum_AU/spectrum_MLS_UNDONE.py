import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.fft import fft, ifft
from scipy.integrate import simpson

from allFunc import *

T = 1                                      # unit: K 
BETA = 1./ (KB * T)   

dt = 1e-15                                 # in s
t_range = np.arange(0.0, 8e-12+dt, dt)     # in s
domega = 0.1e-3                            # in eV
w_range = np.arang(1., 3., domega)         # in eV

# erf param
cutoff = 5e-12
slope = 0.01

excitonic_e = np.array([0.080794824219, 0.081730217568, 0.081821304713, 0.081854660953, 0.082673096323, 0.082678266958, 0.083266678531, 0.083287746427, 0.084415605775, 0.084449844719, 0.084489040104, 0.085492049118, 0.085559936920, 0.085679382996, 0.085697756841, 0.086162572228, 0.086229061493, 0.086240691720, 0.086749006766, 0.087217358603])       # in Hartree
excitonic_e *= AUTOJ                     # in J


# Reading frequency and coupling matrix element data
nExc = 20
w = read_w("w.dat")
coupling = read_Vklq("Vklq-diabatic.dat", len(w), nExc)
w = w[6:]

fig, axs = plt.subplots(1,2, figsize=(20,10))
fig.suptitle(r"Absorption spectra for the MLS, E17 model, 3nmCdSe at T=1K")

axs[0].grid()
axs[1].grid()

total_spectrum = np.zeros(int((4000e12 - 2000e12)/domega))

for exc in range(nExc):
    V = coupling[6:, exc]
    delta_range = (V) / w**2
    print("DONE reading")

    # Calculating the correlation function, and writing to file
    (F_t_Re_E17, F_t_Im_E17) = calc_dephasing_F_E17(t_range, delta_range, w)
    print("DONE calc_dephasing_F_E17")
    
    '''
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
    write_dephasing_F("3nmCdSe_F.dat", t_range, F_t_Re_MC, F_t_Im_MC, F_t_Re_E17, F_t_Im_E17)
    write_delta_corr("3nmCdSe_Delta_Corr.dat", t_range, DC_t_Re_MC, DC_t_Im_MC, DC_t_Re_E17, DC_t_Im_E17)
    '''

    # erf correction to the dephasing functions
    F_t_Re_E17_erf = multiply_erf(cutoff, slope, F_t_Re_E17-np.mean(F_t_Re_E17[-1000:-1]), dt)
    F_t_Im_E17_erf = multiply_erf(cutoff, slope, F_t_Im_E17-np.mean(F_t_Im_E17[-1000:-1]), dt)
    # write_corr("3nmCdSe_T=16_shorttime_corr_erf.dat", t_range, corr_erf)
    print("DONE erf_multiplication")

    # FT the dephasing functions to obtain the spectrum from the correlation function
    delta_E_E17 = excitonic_e[exc] + reorg_E(delta_range, w)    # in J
    # delta_E_MC = excitonic_e                               # in J
    (this_exc_spectrum_w_E17, this_exc_spectrum_I_E17) = spectrum(delta_E_E17, t_range, F_t_Re_E17_erf, F_t_Im_E17_erf, domega)
    print("DONE spectrum_E17")
    # (spectrum_w_MC, spectrum_I_MC) = spectrum(delta_E_MC, t_range, F_t_Re_MC_erf, F_t_Im_MC_erf, domega)
    # print("DONE spectrum_MC")
    write_spectrum("3nmCdSe_spectrum_E17_Exc_%d.dat" % exc, this_exc_spectrum_w_E17, this_exc_spectrum_I_E17)

    total_spectrum += this_exc_spectrum_I_E17

    # plotting
    axs[0].set(xlabel="Energey [meV]", title="Spectra", xlim=(2100, 2500))
    axs[0].plot(this_exc_spectrum_w_E17*1000, this_exc_spectrum_I_E17, label="exciton #%d" % exc)
    axs[0].legend()

axs[1].set(xlabel="Energey [meV]", title="Total Spectra", xlim=(2100, 2500))
axs[1].plot(this_exc_spectrum_w_E17*1000, total_spectrum, label="Total spectra")
axs[1].legend()

fig.tight_layout()
fig.savefig("nanocrystal_MLS.png")
