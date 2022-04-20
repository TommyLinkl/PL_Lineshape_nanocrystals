import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.fft import fft, ifft
from scipy.integrate import simpson

from allFunc import *

# main
# Temp_range = np.array([1, 2, 5, 10, 20, 50, 100, 200, 300, 400])
# Temp_range = np.array([1, 2, 200, 300])
Temp_range = np.array([300])

linewidthQM_range = np.zeros(0)
freqShiftQM_range = np.zeros(0)

k = open("fitted_SuperOhmic_WidthFreqShift_vs_T.dat", "w")
k.write("# Temperature (K)   WidthQM(eV)   FreqShiftQM(eV)\n")

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
    I_file_QM = "fitted_Spectrum_T=%d.dat" % (T)
    F_file_QM = "fitted_F_T=%d.dat" % (T)
    FWHM_file = "FWHM_T=%d.dat" % (T)
    
    '''
    (spectrum_w, spectrum_I, linewidth_QM_MC, freq_shift_QM_MC, F_t_Re_QM_MC, F_t_Im_QM_MC, F_t_Re_erf, F_t_Im_erf, DC_t_Re, DC_t_Im) = calc_spectrum_QM_MC(T, t_range, w_range, cutoff, slope, exc_E, phonon_filename, coupling_filename, I_file_QM_MC, FWHM_file)
    linewidthQM_MC_range = np.append(linewidthQM_MC_range, linewidth_QM_MC)
    freqShiftQM_MC_range = np.append(freqShiftQM_MC_range, freq_shift_QM_MC)
    # axs[0,0].plot(spectrum_w*1000, spectrum_I, label="T=%d K" % T)   
    '''

    (spectrum_w, spectrum_I, linewidth_QM_E17, freq_shift_QM_E17, F_t_Re_QM_E17, F_t_Im_QM_E17, F_t_Re_erf, F_t_Im_erf, DC_t_Re, DC_t_Im) = calc_spectrum_QM_E17_fitted(T, t_range, w_range, cutoff, slope, exc_E, phonon_filename, coupling_filename, I_file_QM, FWHM_file)
    linewidthQM_range = np.append(linewidthQM_range, linewidth_QM_E17)
    freqShiftQM_range = np.append(freqShiftQM_range, freq_shift_QM_E17)
    # axs[0,0].plot(spectrum_w*1000, spectrum_I, label="T=%d K" % T)   

    '''
    (spectrum_w, spectrum_I, linewidth_CL, freq_shift_CL, F_t_Re_CL, F_t_Im_CL, F_t_Re_erf, F_t_Im_erf, DC_t_Re, DC_t_Im) = calc_spectrum_classical_MC(T, t_range, w_range, cutoff, slope, exc_E, phonon_filename, coupling_filename, I_file_CL, FWHM_file)
    linewidthCL_range = np.append(linewidthCL_range, linewidth_CL)
    freqShiftCL_range = np.append(freqShiftCL_range, freq_shift_CL)
    # axs[0,1].plot(spectrum_w*1000, spectrum_I, label="T=%d K" % T)   
    '''

    g = open(F_file_QM, "w")
    g.write("# Time    F_t_Re (QM,E17)      F_t_Im (QM,E17)       F_t_Re_erf       F_t_Im_erf  \n")
    for i in range(len(t_range)):
        g.write("%.20f        %f        %f         %f         %f     \n" % (t_range[i], F_t_Re_QM_E17[i], F_t_Im_QM_E17[i], F_t_Re_erf[i], F_t_Im_erf[i]))
    g.close()


    k.write("%f     %f     %f\n"%(T, linewidth_QM_E17, freq_shift_QM_E17))

