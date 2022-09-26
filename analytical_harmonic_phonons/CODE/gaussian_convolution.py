import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.fft import fft, ifft
from scipy.integrate import simpson

from allFunc import *

# main
w_range = np.zeros(0)
exc_E = 0.075068450260                             # in Hartree
exc_E *= AUTOJ                                     # in J\
Temp_range = np.array([0, 1, 2, 5, 10, 20, 50, 100, 150, 200, 250, 300, 350, 400])
gaussian_range = np.array([1, 2, 5, 10, 20, 30, 40, 50, 100, 200, 300, 400, 500, 600, 700, 800, 1000, 2000, 5000, 10000, 20000, 50000]) # in fs

for gaussian in gaussian_range: 
    print("%d fs" % gaussian)
    sigma_t = gaussian * 1e-15      # seconds
    g = open("results/NEW_WidthFreqShift_vs_T_gaussian_%dfs.dat" % (gaussian), "w")
    g.write("# Temperature (K)  WidthQM_E17(eV)\n")
    
    for T in Temp_range:
        highRes_filename = "results/NEW_F_QM_E17_T=%d.dat" % (T)
        gaussian_F_filename = "results/NEW_gaussian%dfs_F_QM_E17_T=%d.dat" % (gaussian, T)
        gaussian_I_filename = "results/NEW_gaussian%dfs_Spectrum_QM_E17_T=%d.dat" % (gaussian, T)
        
        time = np.array([float(line.strip().split()[0]) for line in open(highRes_filename, 'r') if line[0]!="#"])
        F_Re = np.array([float(line.strip().split()[1]) for line in open(highRes_filename, 'r') if line[0]!="#"])
        F_Im = np.array([float(line.strip().split()[2]) for line in open(highRes_filename, 'r') if line[0]!="#"])
        
        gaussian_factor = np.exp(-time**2/ (2*sigma_t**2))
        Gaussian_F_Re = F_Re * gaussian_factor
        Gaussian_F_Im = F_Im * gaussian_factor
        
        f = open(gaussian_F_filename, "w")
        f.write("# Time    F_t_Re*gaussian (QM_E17)      F_t_Im*gaussian (QM_E17)\n")
        for i in range(len(time)):
            f.write("%.20f    %.20f    %.20f\n" % (time[i], Gaussian_F_Re[i], Gaussian_F_Im[i]))
        f.close()
    
        # FT the dephasing functions to obtain the spectrum from the dipole-dipole correlation function
        delta_E_E17 = exc_E   #  - reorg_E(delta_range, phonons)    # in J
        (spectrum_w, spectrum_I) = spectrum(delta_E_E17, time, Gaussian_F_Re, Gaussian_F_Im, w_range)
        print("DONE spectrum")
        write_spectrum(gaussian_I_filename, spectrum_w, spectrum_I)
        
        g.write("%d   %f\n"%(T, calc_FWHM(spectrum_w, spectrum_I)[0]))
    g.close()
