import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.fft import fft, ifft
from scipy.integrate import simpson

from allFunc import *

# main
# Temp_range = np.array([1, 2, 5, 10, 20, 50, 100, 200, 300, 400])
Temp_range = np.array([1, 2, 5, 10])

linewidthCL_range = np.zeros(0)
freqShiftCL_range = np.zeros(0)


k = open("NEW_WidthFreqShift_vs_T.dat", "w")
k.write("# Temperature (K)   WidthCL(eV)   FreqShiftCL(eV)\n")

'''
fig, axs = plt.subplots(3,2, figsize=(12,18))
fig.suptitle("Spectral linewidth, frequency shifts, and temperature")      
for a in range(3): 
    for b in range(2): 
        axs[a,b].grid()

'''

for T in Temp_range: 
    dt = 1e-15                                         # in s
    t_range = np.arange(0.0, 1e-11+dt, dt)             # in s
    # erf param
    cutoff = 8e-12                                     # in s
    slope = 0.01
    
    domega = 0.1e-3                                    # in eV
    w_range = np.arange(1.5, 2.5, domega)                # in eV
    
    exc_E = 0.080311203446                             # in Hartree
    exc_E *= AUTOJ                                     # in J
    
    phonon_filename = "w.dat"
    coupling_filename = "Vklq-diabatic.dat"
    # I_file_QM_E17 = "Spectrum_QM_E17_T=%d.dat" % (T)
    # I_file_QM_MC = "Spectrum_QM_MC_T=%d.dat" % (T)
    I_file_CL = "NEW_Spectrum_CL_T=%d.dat" % (T)
    F_file_CL = "NEW_F_CL_T=%d.dat" % (T)
    FWHM_file = "FWHM_T=%d.dat" % (T)
    
    '''
    (spectrum_w, spectrum_I, linewidth_QM_MC, freq_shift_QM_MC, F_t_Re_QM_MC, F_t_Im_QM_MC, F_t_Re_erf, F_t_Im_erf, DC_t_Re, DC_t_Im) = calc_spectrum_QM_MC(T, t_range, w_range, cutoff, slope, exc_E, phonon_filename, coupling_filename, I_file_QM_MC, FWHM_file)
    linewidthQM_MC_range = np.append(linewidthQM_MC_range, linewidth_QM_MC)
    freqShiftQM_MC_range = np.append(freqShiftQM_MC_range, freq_shift_QM_MC)
    # axs[0,0].plot(spectrum_w*1000, spectrum_I, label="T=%d K" % T)   
    
    (spectrum_w, spectrum_I, linewidth_QM_E17, freq_shift_QM_E17, F_t_Re_QM_E17, F_t_Im_QM_E17, F_t_Re_erf, F_t_Im_erf, DC_t_Re, DC_t_Im) = calc_spectrum_QM_E17(T, t_range, w_range, cutoff, slope, exc_E, phonon_filename, coupling_filename, I_file_QM_E17, FWHM_file)
    linewidthQM_E17_range = np.append(linewidthQM_E17_range, linewidth_QM_E17)
    freqShiftQM_E17_range = np.append(freqShiftQM_E17_range, freq_shift_QM_E17)
    axs[0,0].plot(spectrum_w*1000, spectrum_I, label="T=%d K" % T)   
    '''

    (spectrum_w, spectrum_I, linewidth_CL, freq_shift_CL, F_t_Re_CL, F_t_Im_CL, F_t_Re_erf, F_t_Im_erf, DC_t_Re, DC_t_Im) = calc_spectrum_classical_MC(T, t_range, w_range, cutoff, slope, exc_E, phonon_filename, coupling_filename, I_file_CL, FWHM_file)
    linewidthCL_range = np.append(linewidthCL_range, linewidth_CL)
    freqShiftCL_range = np.append(freqShiftCL_range, freq_shift_CL)
    # axs[0,1].plot(spectrum_w*1000, spectrum_I, label="T=%d K" % T)   

    g = open(F_file_CL, "w")
    g.write("# Time    F_t_Re (Classical)      F_t_Im (Classical)       F_t_Re_erf       F_t_Im_erf  \n")
    for i in range(len(t_range)):
        g.write("%.20f        %f        %f         %f         %f     \n" % (t_range[i], F_t_Re_CL[i], F_t_Im_CL[i], F_t_Re_erf[i], F_t_Im_erf[i]))
    g.close()


    k.write("%f     %f     %f\n"%(T, linewidth_CL, freq_shift_CL))

