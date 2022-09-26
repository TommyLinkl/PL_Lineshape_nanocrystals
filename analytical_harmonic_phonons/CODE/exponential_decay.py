import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.fft import fft, ifft
from scipy.integrate import simpson

from allFunc import *

# main
w_range = np.zeros(0)
exc_E = 0.073095124824                             # in Hartree
exc_E *= AUTOJ                                     # in J\
Temp_range = np.array([0, 1, 2, 4, 5, 10, 20, 30, 50, 60, 100, 150, 200, 220, 250, 290, 300, 350, 400])
# oneOverGamma_range = np.array([1, 10, 20, 30, 40, 50, 100, 500, 1000]) # in fs
oneOverGamma_range = np.array([200]) # in fs

for oneOverGamma in oneOverGamma_range: 
    print("1/gamma = %d fs" % oneOverGamma)
    gamma = 1/(oneOverGamma * 1e-15)  # in 1/s

    # g = open("results/NEW_WidthFreqShift_vs_T_exponential%dfs.dat" % (oneOverGamma), "w")
    g = open("results/NEW_WidthFreqShift_vs_T_exponential%dfs_gaussian500fs.dat" % (oneOverGamma), "w")
    g.write("# Temperature (K)  WidthQM_E17(eV)\n")
    
    for T in Temp_range:
        # gaussian_filename = "results/NEW_F_QM_E17_T=%d.dat" % (T)
        # exponential_F_filename = "results/NEW_exponential%dfs_F_QM_E17_T=%d.dat" % (oneOverGamma, T)
        # exponential_I_filename = "results/NEW_exponential%dfs_Spectrum_QM_E17_T=%d.dat" % (oneOverGamma, T)
        gaussian_filename = "results/NEW_gaussian500fs_F_QM_E17_T=%d.dat" % (T)
        exponential_F_filename = "results/NEW_exponential%dfs_gaussian500fs_F_QM_E17_T=%d.dat" % (oneOverGamma, T)
        exponential_I_filename = "results/NEW_exponential%dfs_gaussian500fs_Spectrum_QM_E17_T=%d.dat" % (oneOverGamma, T)
        
        time = np.array([float(line.strip().split()[0]) for line in open(gaussian_filename, 'r') if line[0]!="#"])
        F_Re = np.array([float(line.strip().split()[1]) for line in open(gaussian_filename, 'r') if line[0]!="#"])
        F_Im = np.array([float(line.strip().split()[2]) for line in open(gaussian_filename, 'r') if line[0]!="#"])
        
        exponential_factor = np.exp(-gamma * time)
        exponential_F_Re = F_Re * exponential_factor
        exponential_F_Im = F_Im * exponential_factor
        
        f = open(exponential_F_filename, "w")
        # f.write("# Time    F_t_Re*exponential (QM_E17)      F_t_Im*exponential (QM_E17)\n")
        f.write("# Time    F_t_Re*exponential*gaussian (QM_E17)      F_t_Im*exponential*gaussian (QM_E17)\n")
        for i in range(len(time)):
            f.write("%.20f    %.20f    %.20f\n" % (time[i], exponential_F_Re[i], exponential_F_Im[i]))
        f.close()
    
        # FT the dephasing functions to obtain the spectrum from the dipole-dipole correlation function
        delta_E_E17 = exc_E   #  - reorg_E(delta_range, phonons)    # in J
        (spectrum_w, spectrum_I) = spectrum(delta_E_E17, time, exponential_F_Re, exponential_F_Im, w_range)
        print("DONE spectrum")
        write_spectrum(exponential_I_filename, spectrum_w, spectrum_I)
        
        g.write("%d   %f\n"%(T, calc_FWHM(spectrum_w, spectrum_I)[0]))
    g.close()
