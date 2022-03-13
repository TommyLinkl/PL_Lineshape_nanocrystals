import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.fft import fft, ifft
from scipy.integrate import simpson

def write_dephasing_F(filename, t_range, F_Re_MC, F_Im_MC, F_Re_E17, F_Im_E17):
    f = open(filename, "w")
    f.write("# time        F_Re_MC      F_Im_MC      F_Re_E17      F_Im_E17 \n")
    for i in range(len(t_range)):
        f.write("%.20f       %.10f       %.10f       %.10f       %.10f\n" % (t_range[i], F_Re_MC[i], F_Im_MC[i], F_Re_E17[i], F_Im_E17[i]))
    f.close()
    return

def write_delta_corr(filename, t_range, DC_Re_MC, DC_Im_MC, DC_Re_E17, DC_Im_E17):
    f = open(filename, "w")
    f.write("# time        DeltaCorr_Re_MC      DeltaCorr_Im_MC      DeltaCorr_Re_E17      DeltaCorr_Im_E17 \n")
    for i in range(len(t_range)):
        f.write("%.10f       %.10f       %.10f       %.10f       %.10f\n" % (t_range[i], DC_Re_MC[i], DC_Im_MC[i], DC_Re_E17[i], DC_Im_E17[i]))
    f.close()
    return

def write_spectrum(filename, w_range, spectrum):
    f = open(filename, "w")
    f.write("# Frequency w         I(w)\n")
    for i in range(len(w_range)):
        f.write("%f       %.20f\n" % (w_range[i], spectrum[i]))
    f.close()
    return

def write_FWHM(filename, w, I_w, l_index, r_index):
    f = open(filename, "w")
    f.write("# FWHM width       points_w      points_I\n")
    f.write("%f        %f       %.20f\n" % (w[r_index]-w[l_index], w[l_index], I_w[l_index]))
    f.write("%f        %f       %.20f\n" % (w[r_index]-w[l_index], w[r_index], I_w[r_index]))
    f.close()
    return
