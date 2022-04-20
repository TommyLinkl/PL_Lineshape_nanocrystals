import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.fft import fft, ifft
from scipy.integrate import simpson

from allFunc import *

# main
Temp_range = np.array([0, 1, 2, 5, 10, 20, 30, 40])
Temp_range = np.append(Temp_range, np.arange(50, 450, 50))

linewidthQM_range = np.zeros(0)
linewidthCL_range = np.zeros(0)
freqShiftQM_range = np.zeros(0)
freqShiftCL_range = np.zeros(0)
k = open("T_width_FWHMFreqShift_vs_T.dat", "w")
k.write("# Temperature (K)         Quantum FWHM (eV)        Classical FWHM (eV)       Quantum FreqShift (eV)      Classicla FreqShift (eV)   \n")

fig, axs = plt.subplots(3,2, figsize=(12,18))
fig.suptitle("Spectral linewidth, frequency shifts, and temperature")      
for a in range(3): 
    for b in range(2): 
        axs[a,b].grid()

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
    I_file_QM = "T_width_QMSpectrum_E17_5ps_T=%d.dat" % (T)
    I_file_CL = "T_width_CLSpectrum_E17_5ps_T=%d.dat" % (T)
    FWHM_file = "T_width_FWHM_T=%d.dat" % (T)
    
    (spectrum_w, spectrum_I, linewidth_QM, freq_shift_QM, F_t_Re, F_t_Im, F_t_Re_erf, F_t_Im_erf, DC_t_Re, DC_t_Im) = calc_spectrum_E17(T, t_range, w_range, cutoff, slope, exc_E, phonon_filename, coupling_filename, I_file_QM, FWHM_file)
    linewidthQM_range = np.append(linewidthQM_range, linewidth_QM)
    freqShiftQM_range = np.append(freqShiftQM_range, freq_shift_QM)
    axs[0,0].plot(spectrum_w*1000, spectrum_I, label="T=%d K" % T)   

    (spectrum_w, spectrum_I, linewidth_CL, freq_shift_CL, F_t_Re, F_t_Im, F_t_Re_erf, F_t_Im_erf, DC_t_Re, DC_t_Im) = calc_spectrum_classical_MC(T, t_range, w_range, cutoff, slope, exc_E, phonon_filename, coupling_filename, I_file_CL, FWHM_file)
    linewidthCL_range = np.append(linewidthCL_range, linewidth_CL)
    freqShiftCL_range = np.append(freqShiftCL_range, freq_shift_CL)
    axs[0,1].plot(spectrum_w*1000, spectrum_I, label="T=%d K" % T)   

    k.write("%f     %f     %f     %f     %f  \n"%(T, linewidth_QM, linewidth_CL, freq_shift_QM, freq_shift_CL))


axs[0,0].set(xlabel="Energey [meV]", title="Normalized Spectra, QM", xlim=(1800, 2600))
axs[0,1].set(xlabel="Energey [meV]", title="Normalized Spectra, Classical", xlim=(1800, 2600))

axs[1,0].set(xlabel="Temperature [K]", ylabel="Linewidth [meV]", title="FWHM vs. T")
axs[1,0].plot(Temp_range, linewidthQM_range*1000, "*-", label="QM")
axs[1,0].plot(Temp_range, linewidthCL_range*1000, "*-", label="Classical")

axs[1,1].set(xlabel="Inverse Temperature [1/K]", ylabel="Linewidth [meV]", title=r"FWHM vs. $\beta$")
axs[1,1].plot(1/Temp_range, linewidthQM_range*1000, "*-", label="QM")
axs[1,1].plot(1/Temp_range, linewidthCL_range*1000, "*-", label="Classical")

axs[2,0].set(xlabel="Temperature [K]", ylabel="Frequency shift [meV]", title=r"$\Delta \omega$ vs. T")
axs[2,0].plot(Temp_range, freqShiftQM_range*1000, "*-", label="QM")
axs[2,0].plot(Temp_range, freqShiftCL_range*1000, "*-", label="Classical")

axs[2,1].set(xlabel="Inverse Temperature [1/K]", ylabel="Frequency shift [meV]", title=r"$\Delta \omega$ vs. $\beta$")
axs[2,1].plot(1/Temp_range, freqShiftQM_range*1000, "*-", label="QM")
axs[2,1].plot(1/Temp_range, freqShiftCL_range*1000, "*-", label="Classical")

for a in range(3): 
    for b in range(2): 
        axs[a,b].legend()

fig.tight_layout()
fig.savefig("T_width_freqShift.png")

