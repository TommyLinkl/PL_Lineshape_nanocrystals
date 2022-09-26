import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.fft import fft, ifft
from scipy.integrate import simpson

from allFunc import *

# main
Temp_range = np.array([0, 1, 2, 5, 10, 20, 50, 100, 150, 200, 250, 300, 350, 400])
# Temp_range = np.array([0])
# Shrink_range = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
Shrink_range = np.array([1.0])

linewidthQM_E17_range = np.zeros(0)
freqShiftQM_E17_range = np.zeros(0)

k = open("results/NEW_WidthFreqShift_vs_T.dat", "w")
k.write("# Temperature (K)   WidthQM_E17(eV)   FreqShiftQM_E17(eV)\n")

for T in Temp_range: 
    shrink = 1.0
    dt = 1e-15                                         # in s
    t_range = np.arange(0.0, 5e-11+dt, dt)               # in s
    
    # erf param
    cutoff = 8e-6                                     # in s
    slope = 0.0001
    
    domega = 0.1e-3                                    # in eV
    w_range = np.arange(1.5, 2.5, domega)                # in eV
    
    exc_E = 0.075068450260   # 0.080311203446                             # in Hartree
    exc_E *= AUTOJ                                     # in J
    
    phonon_filename = "w.dat"
    coupling_filename = "NEW_Vklq-diabatic.dat"
    I_file_QM_E17 = "results/NEW_Spectrum_QM_E17_T=%d.dat" % (T)
    F_file_QM_E17 = "results/NEW_F_QM_E17_T=%d.dat" % (T)
    FWHM_file = "FWHM_T=%d.dat" % (T)

    (spectrum_w, spectrum_I, linewidth_QM_E17, freq_shift_QM_E17, F_t_Re_QM_E17, F_t_Im_QM_E17, F_t_Re_erf, F_t_Im_erf, DC_t_Re, DC_t_Im) = calc_spectrum_QM_E17(T, t_range, w_range, cutoff, slope, exc_E, phonon_filename, coupling_filename, I_file_QM_E17, FWHM_file, shrink)
    linewidthQM_E17_range = np.append(linewidthQM_E17_range, linewidth_QM_E17)
    freqShiftQM_E17_range = np.append(freqShiftQM_E17_range, freq_shift_QM_E17)

    g = open(F_file_QM_E17, "w")
    g.write("# Time    F_t_Re (QM_E17)      F_t_Im (QM_E17)       F_t_Re_erf       F_t_Im_erf  \n")
    for i in range(len(t_range)):
        g.write("%.20f        %.20f        %.20f         %.20f         %.20f     \n" % (t_range[i], F_t_Re_QM_E17[i], F_t_Im_QM_E17[i], F_t_Re_erf[i], F_t_Im_erf[i]))
    g.close()

    k.write("%f     %f     %f\n"%(T, linewidth_QM_E17, freq_shift_QM_E17))

