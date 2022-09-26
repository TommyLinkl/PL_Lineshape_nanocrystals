import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.fft import fft, ifft
from scipy.integrate import simpson
from scipy.interpolate import CubicSpline

from allFunc import *

# main
# Temp_range = np.array([0, 1, 2, 4, 5, 10, 20, 30, 50, 60, 100, 150, 200, 220, 250, 290, 300, 350, 400])
Temp_range = np.array([300])
# gaussian_range = np.array([1, 2, 5, 10, 20, 30, 40, 50, 100, 200, 300, 400, 500, 600, 700, 800, 1000, 2000, 5000, 10000, 20000, 50000]) # in fs
gaussian_range = np.array([500]) # in fs

exc_energy_range = np.array([float(line.strip().split()[1]) for line in open("exciton.dat", 'r') if line[0]!="#"])[:100]
mu_x_range = np.array([float(line.strip().split()[4]) for line in open("OS.dat", 'r') if line[0]!="#"])[:100]
mu_y_range = np.array([float(line.strip().split()[5]) for line in open("OS.dat", 'r') if line[0]!="#"])[:100]
mu_z_range = np.array([float(line.strip().split()[6]) for line in open("OS.dat", 'r') if line[0]!="#"])[:100]
mu_squared_range = mu_x_range**2 + mu_y_range**2 + mu_z_range**2

k = open("results/WidthFreqShift_vs_T.dat", "w")
k.write("# Temperature (K)   ExcNumber    Gaussian_std(fs)    WidthQM_E17(eV)     scalingRatio\n")

for gaussian in gaussian_range: 
    print("\n%d fs" % gaussian)
    sigma_t = gaussian * 1e-15      # seconds

    for T in Temp_range: 
        shrink = 1.0
        dt = 1e-15                                         # in s
        t_range = np.arange(0.0, 5e-11+dt, dt)               # in s
        domega = 0.1e-3                                    # in eV
        w_range = np.arange(1.5, 2.5, domega)                # in eV
        
        # erf param
        cutoff = 8e-6                                     # in s
        slope = 0.0001
        
        manyExc_w_range = np.arange(-0.050, 4.300, 0.00001)
        manyExc_I_range = np.zeros(len(manyExc_w_range))

        for exc_index in range(20): 
            print("\nExciton #%d" % exc_index)
            exc_E = exc_energy_range[exc_index]                # in Hartree
            exc_E *= AUTOJ                                     # in J
            
            phonon_filename = "w.dat"
            coupling_filename = "Vklq_diabatic.dat"
            I_file_QM_E17 = "results/Exc%d_Spectrum_QM_E17_T=%d.dat" % (exc_index, T)
            F_file_QM_E17 = "results/Exc%d_F_QM_E17_T=%d.dat" % (exc_index, T)
            FWHM_file = "FWHM_T=%d.dat" % (T)

            (spectrum_w, spectrum_I, linewidth_QM_E17, freq_shift_QM_E17, F_t_Re_QM_E17, F_t_Im_QM_E17, F_t_Re_erf, F_t_Im_erf, DC_t_Re, DC_t_Im) = calc_spectrum_QM_E17(T, 20, exc_index, t_range, w_range, cutoff, slope, exc_E, phonon_filename, coupling_filename, I_file_QM_E17, FWHM_file, shrink)
            k.write("%d     %d      inf       %f     0\n"%(T, exc_index, linewidth_QM_E17))

            g = open(F_file_QM_E17, "w")
            g.write("# Time    F_t_Re (QM_E17)      F_t_Im (QM_E17)       F_t_Re_erf       F_t_Im_erf  \n")
            for i in range(len(t_range)):
               g.write("%.20f        %.20f        %.20f         %.20f         %.20f     \n" % (t_range[i], F_t_Re_QM_E17[i], F_t_Im_QM_E17[i], F_t_Re_erf[i], F_t_Im_erf[i]))
            g.close()

            # Get Gaussian-convolved spectrum
            gaussian_F_filename = "results/gaussian%dfs_Exc%d_F_QM_E17_T=%d.dat" % (gaussian, exc_index, T)
            gaussian_I_filename = "results/gaussian%dfs_Exc%d_Spectrum_QM_E17_T=%d.dat" % (gaussian, exc_index, T)
            time = t_range
            F_Re = F_t_Re_QM_E17
            F_Im = F_t_Im_QM_E17
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
            (Gaussian_spectrum_w, Gaussian_spectrum_I) = spectrum(delta_E_E17, time, Gaussian_F_Re, Gaussian_F_Im, w_range)
            print("DONE Gaussian spectrum")
            write_spectrum(gaussian_I_filename, Gaussian_spectrum_w, Gaussian_spectrum_I)

            # Interpolate the Gaussian-convolved spectrum for summation
            interpolate_cs = CubicSpline(Gaussian_spectrum_w, Gaussian_spectrum_I)
            interpolated_gaussian_I_filename = "results/interpolate_gaussian%dfs_Exc%d_Spectrum_QM_E17_T=%d.dat" % (gaussian, exc_index, T)
            write_spectrum(interpolated_gaussian_I_filename, manyExc_w_range, interpolate_cs(manyExc_w_range))
            
            # rescale manyExc spectra
            manyExc_I_range += interpolate_cs(manyExc_w_range) * np.exp(-(exc_E-exc_energy_range[0]*AUTOJ)/(KB*T)) * (mu_squared_range[exc_index]/mu_squared_range[0])**2
            print(np.exp(-(exc_E-exc_energy_range[0]*AUTOJ)/(KB*T)))
            print((mu_squared_range[exc_index]/mu_squared_range[0])**2)
            print(manyExc_I_range)

            k.write("%d     %d      %d       %f     %.20f\n"%(T, exc_index, gaussian, calc_FWHM(Gaussian_spectrum_w, Gaussian_spectrum_I)[0], np.exp(-(exc_E-exc_energy_range[0]*AUTOJ)/(KB*T)) * (mu_squared_range[exc_index]/mu_squared_range[0])**2))

        # write manyExc spectrum for this gaussian broadening
        manyExc_I_filename = "results/manyExc_Spectrum_QM_E17_T=%d.dat" % (T)
        write_spectrum(manyExc_I_filename, manyExc_w_range, manyExc_I_range/max(manyExc_I_range))

        # measure the linewidth of this manyExc spectrum
        k.write("%d     all    %d       %f     0\n"%(T, gaussian, calc_FWHM(manyExc_w_range, manyExc_I_range)[0]))

k.close()