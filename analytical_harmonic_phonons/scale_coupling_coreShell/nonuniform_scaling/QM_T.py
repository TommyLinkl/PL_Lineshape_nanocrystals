import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.fft import fft, ifft
from scipy.integrate import simpson

from allFunc import *

# main
# Temp_range = np.array([0, 1, 2, 5, 10, 20, 50, 100, 150, 200, 250, 300, 350, 400])
Temp_range = np.array([0])
ac_low = 0
ac_upper = 2.5
alpha_range = np.arange(0, 1.1, 0.1)
op_low = 6
op_upper = 12
beta_range = np.arange(1, 11, 1)
print("Acoustic range: %f %f" % (ac_low, ac_upper))
print(alpha_range)
print("Optical range: %f %f" % (op_low, op_upper))
print(beta_range)

dt = 1e-15                                         # in s
t_range = np.arange(0.0, 5e-11+dt, dt)               # in s
# erf param
cutoff = 8e-6                                     # in s
slope = 0.0001
domega = 0.1e-3                                    # in eV
w_range = np.arange(1.5, 2.5, domega)                # in eV
exc_E = 0.0367493    # 0.080311203446             # in Hartree
exc_E *= AUTOJ                                     # in J

k = open("WidthFreqShift_vs_T.dat", "w")
k.write("# Temperature (K)   alpha   beta     WidthQM_E17(eV)   FreqShiftQM_E17(eV)\n")

for T in Temp_range: 
    for alpha_index in range(len(alpha_range)): 
        for beta_index in range(len(beta_range)):
            print("alpha = %f, beta = %f" % (alpha_range[alpha_index], beta_range[beta_index]))
            phonon_filename = "w.dat"
            coupling_filename = "Vklq_trim.dat"

            nExc = 1
            phonons = read_w(phonon_filename)
            coupling = read_Vklq(coupling_filename, len(phonons), nExc, 1.0)
            phonons = phonons[6:]
            V = coupling[6:, 0]
            delta_range = (V) / phonons**2

            rescale_Vklq(phonons, V, ac_low, ac_upper, alpha_range[alpha_index], \
                op_low, op_upper, beta_range[beta_index], \
                "Vklq_aIndex%d_bIndex%d.dat"%(alpha_index, beta_index))
            print("DONE rescaling V")

            I_file_QM_E17 = "Spectrum_QM_E17_T=%d_aIndex%d_bIndex%d.dat" % (T, alpha_index, beta_index)
            F_file_QM_E17 = "F_QM_E17_T=%d_aIndex%d_bIndex%d.dat" % (T, alpha_index, beta_index)
            FWHM_file = "FWHM_T=%d.dat" % (T)

            (spectrum_w, spectrum_I, linewidth_QM_E17, freq_shift_QM_E17, F_t_Re_QM_E17, F_t_Im_QM_E17, \
                F_t_Re_erf, F_t_Im_erf, DC_t_Re, DC_t_Im) = calc_spectrum_QM_E17(T, t_range, \
                w_range, cutoff, slope, exc_E, phonon_filename, \
                "Vklq_aIndex%d_bIndex%d.dat"%(alpha_index, beta_index), I_file_QM_E17, FWHM_file, 1.0)

            g = open(F_file_QM_E17, "w")
            g.write("# Time    F_t_Re (QM_E17)      F_t_Im (QM_E17)       F_t_Re_erf       F_t_Im_erf  \n")
            for i in range(len(t_range)):
                g.write("%.20f        %.20f        %.20f         %.20f         %.20f     \n" % (t_range[i], F_t_Re_QM_E17[i], F_t_Im_QM_E17[i], F_t_Re_erf[i], F_t_Im_erf[i]))
            g.close()

            k.write("%f     %d      %d      %f     %f\n"%(T, alpha_range[alpha_index], \
                beta_range[beta_index], linewidth_QM_E17, freq_shift_QM_E17))

