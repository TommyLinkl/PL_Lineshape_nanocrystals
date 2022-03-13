import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.fft import fft, ifft
from scipy.integrate import simpson

from allFunc import *

# main
Temp_range = np.array([0, 1, 2])    # , 5, 10, 20, 30, 40])
# Temp_range = np.append(Temp_range, np.arange(50, 450, 50))

linewidth_range = np.zeros(0)
k = open("T_width_FWHM_vs_T.dat", "w")
k.write("# Temperature (K)         FWHM (eV)")

fig, axs = plt.subplots(3, figsize=(5,15))
fig.suptitle("Spectral linewidth and temperature")      
axs[0].grid()
axs[1].grid()
axs[2].grid()

for T in Temp_range: 
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
    I_file = "Spectrum_E17_5ps_T=%d.dat" % (T)
    FWHM_file = "FWHM_T=%d.dat" % (T)
    
    (spectrum_w, spectrum_I, linewidth, F_t_Re, F_t_Im, F_t_Re_erf, F_t_Im_erf, DC_t_Re, DC_t_Im) = calc_spectrum_E17(T, t_range, w_range, cutoff, slope, exc_E, phonon_filename, coupling_filename, I_file, FWHM_file)
    linewidth_range = np.append(linewidth_range, linewidth)
    k.write("%f         %f"%(T, linewidth))
    axs[0].plot(spectrum_w*1000, spectrum_I, label="T=%d K" % T)   


axs[0].set(xlabel="Energey [meV]", title="Normalized Spectra", xlim=(1800, 2600))
axs[1].set(xlabel="Temperature [K]", ylabel="Linewidth [meV]", title="FWHM vs. T")
axs[1].plot(Temp_range, linewidth_range*1000, "*-")
axs[2].set(xlabel="Inverse Temperature [1/K]", ylabel="Linewidth [meV]", title=r"FWHM vs. $\beta$")
axs[2].plot(1/Temp_range, linewidth_range*1000, "*-")
axs[0].legend()
axs[1].legend()
axs[2].legend()
fig.tight_layout()
fig.savefig("T_width.png")

