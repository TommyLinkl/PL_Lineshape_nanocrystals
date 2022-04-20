import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.fft import fft, ifft
from scipy.integrate import simpson

# Settings
T = 8.
BETA = 1/T
spectral_density = "Ohmic"

# CONST
PI = 3.14159265359
HBAR = 1.     # 1.054571817e-34 Js

domega = 0.04
dt = 0.05
t_range = np.arange(0.0, 10.0+dt, dt)
omega_10 = 0.

corr_file = "pure_linear_shift_corr_Ohmic_T=8_dw0.04.dat"
spectrum_file = "pure_linear_shift_I_Ohmic_T=8_dw0.04.dat"

def coth(x):
    return (np.exp(x) + np.exp(-x))/(np.exp(x) - np.exp(-x))

def J(omega):
    if (spectral_density=="Ohmic"):
        return 0.25 * omega * np.exp(-omega)
    elif (spectral_density=="Super-Ohmic"):
        return 0.05 * omega**3 * np.exp(-omega)

def DCL_Re(omega, time):
    return J(omega) * (np.cos(omega*time)-1) / (HBAR**2 * omega**4 * BETA)

def Int_DCL_Re(time):
    w_range = np.arange(0.0+domega, 20.0+domega, domega)
    y_range = np.zeros(0)
    for w in w_range:
        y_range = np.append(y_range, DCL_Re(w, time))
    return simpson(y_range, w_range)

def DCL_Im(omega, time):
    return J(omega) * time / (HBAR* omega**2)

def Int_DCL_Im(time):
    w_range = np.arange(0.0+domega, 20.0+domega, domega)
    y_range = np.zeros(0)
    for w in w_range:
        y_range = np.append(y_range, DCL_Im(w, time))
    return simpson(y_range, w_range)

def Int_ACL_Re(t):
    return Int_DCL_Re(t)

def ACL_Im(omega, time):
    return J(omega) * (np.sin(omega*time)/(2*HBAR*omega**3) + time/(2*HBAR*omega**2))

def Int_ACL_Im(time):
    w_range = np.arange(0.0+domega, 20.0+domega, domega)
    y_range = np.zeros(0)
    for w in w_range:
        y_range = np.append(y_range, ACL_Im(w, time))
    return simpson(y_range, w_range)

def corr_FQM_Re(t):
    return 0.

def corr_FQM_Im(t):
    return 0.

def corr_DCL_Re(t):
    return np.exp(Int_DCL_Re(t)) * np.cos(Int_DCL_Im(t))

def corr_DCL_Im(t):
    return np.exp(Int_DCL_Re(t)) * np.sin(Int_DCL_Im(t))

def corr_ACL_Re(t):
    return np.exp(Int_ACL_Re(t)) * np.cos(Int_ACL_Im(t))

def corr_ACL_Im(t):
    return np.exp(Int_ACL_Re(t)) * np.sin(Int_ACL_Im(t))

def write_corr(filename, t_range, c_FQM_Re, c_FQM_Im, c_DCL_Re, c_DCL_Im, c_ACL_Re, c_ACL_Im): 
    f = open(filename, "w")
    f.write('# Analytical correlation function\n')
    f.write("# time   C_FQM_Re    C_FQM_Im    C_DCL_Re    C_DCL_Im    C_ACL_Re   C_ACL_Im\n")
    for i in range(len(t_range)): 
        f.write("%f      %f      %f      %f      %f      %f      %f\n" % (t_range[i], c_FQM_Re[i], c_FQM_Im[i], c_DCL_Re[i], c_DCL_Im[i], c_ACL_Re[i], c_ACL_Im[i]))
    f.close()
    return

def multiply_erf(cutoff, slope, corr_t, dt):
    corr_erf = np.zeros(0)
    for i in range(len(t_range)):
        corr_erf = np.append(corr_erf, corr_t[i]* (-0.5*math.erf(slope*(float(i)-cutoff/dt))+0.5))
    return corr_erf

def spectrum(omega_10, t_range, corr_Re, corr_Im, domega):
    I = np.zeros(0)
    w_range = np.arange(-20.0, 20.0, domega)

    for w in w_range:
        integrand = np.zeros(0)
        for i in range(len(t_range)):
            integrand = np.append(integrand, 1/PI * np.cos((omega_10 - w) * t_range[i]) * corr_Re[i] + 1/PI * np.sin((omega_10 - w) * t_range[i]) * corr_Im[i])
        I_w = simpson(integrand, x=t_range)
        I = np.append(I, I_w)
    return (w_range, I)

def write_spectrum(filename, FQM_w, FQM_I, DCL_I, ACL_I):
    f = open(filename, "w")
    f.write("# FQM_w,     FQM_I,     DCL_I,     ACL_I\n")
    for i in range(len(FQM_w)):
        f.write("%f       %f       %f       %f\n" % (FQM_w[i], FQM_I[i], DCL_I[i], ACL_I[i]))
    f.close()
    return

# main

# Calc FQM, DCL, ACL analytical correlation functions, and write to file
c_FQM_Re = np.zeros(0)
c_FQM_Im = np.zeros(0)
c_DCL_Re = np.zeros(0)
c_DCL_Im = np.zeros(0)
c_ACL_Re = np.zeros(0)
c_ACL_Im = np.zeros(0)

for t in t_range: 
    c_FQM_Re = np.append(c_FQM_Re, corr_FQM_Re(t))
    c_FQM_Im = np.append(c_FQM_Im, corr_FQM_Im(t))
    c_DCL_Re = np.append(c_DCL_Re, corr_DCL_Re(t))
    c_DCL_Im = np.append(c_DCL_Im, corr_DCL_Im(t))
    c_ACL_Re = np.append(c_ACL_Re, corr_ACL_Re(t))
    c_ACL_Im = np.append(c_ACL_Im, corr_ACL_Im(t))

write_corr(corr_file, t_range, c_FQM_Re, c_FQM_Im, c_DCL_Re, c_DCL_Im, c_ACL_Re, c_ACL_Im)

# FT -> Spectrum, and write to file

# (FQM_w, FQM_I) = spectrum(omega_10, t_range, c_FQM_Re-corr_FQM_Re(1e5), c_FQM_Im-corr_FQM_Im(1e5), domega)
# (DCL_w, DCL_I) = spectrum(omega_10, t_range, c_DCL_Re, c_DCL_Im, domega)
# (ACL_w, ACL_I) = spectrum(omega_10, t_range, c_ACL_Re-corr_ACL_Re(1e5), c_ACL_Im-corr_ACL_Im(1e5), domega)
# write_spectrum(spectrum_file, FQM_w, FQM_I, DCL_I, ACL_I)
