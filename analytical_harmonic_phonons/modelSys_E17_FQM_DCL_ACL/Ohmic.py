import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.fft import fft, ifft
from scipy.integrate import simpson

T = 1
BETA = 1./T
PI = 3.14159265359

domega = 0.01
dt = 0.1
t_range = np.arange(0.0, 30.0+dt, dt)
omega_10 = 0.

corr_file = "corr_Ohmic_T=1.dat"
spectrum_file = "I_Ohmic_T=1.dat"

def coth(x):
    return (np.exp(x) + np.exp(-x))/(np.exp(x) - np.exp(-x))

def J(omega):
    # Ohmic spectral density
    return 0.25 * omega * np.exp(-omega)

def FQM_Re(omega, time):
    return 0.5 * J(omega) * omega * coth(BETA * omega / 2) * (np.cos(omega*time) - 1)

def Int_FQM_Re(time): 
    w_range = np.arange(0.0+domega, 20.0+domega, domega)
    y_range = np.zeros(0)
    for w in w_range:
        y_range = np.append(y_range, FQM_Re(w, time))
    return simpson(y_range, w_range)

def Int_FQM_Im(t):
    return 0.25 * t * (t**2 -3) / (t**2 + 1)**3

def Int_DCL_Re(t):
    return -0.25 / BETA * t**2 * (t**2 + 3)/(t**2 + 1)**2

def Int_DCL_Im(t):
    return -0.75 * t

def Int_ACL_Re(t):
    return Int_DCL_Re(t)

def Int_ACL_Im(t):
    return Int_FQM_Im(t)

def corr_FQM_Re(t): 
    return np.exp(Int_FQM_Re(t)) * np.cos(Int_FQM_Im(t))

def corr_FQM_Im(t):
    return np.exp(Int_FQM_Re(t)) * np.sin(Int_FQM_Im(t))

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
(FQM_w, FQM_I) = spectrum(omega_10, t_range, c_FQM_Re-corr_FQM_Re(1e5), c_FQM_Im-corr_FQM_Im(1e5), domega)
(DCL_w, DCL_I) = spectrum(omega_10, t_range, c_DCL_Re, c_DCL_Im, domega)
(ACL_w, ACL_I) = spectrum(omega_10, t_range, c_ACL_Re-corr_ACL_Re(1e5), c_ACL_Im-corr_ACL_Im(1e5), domega)
write_spectrum(spectrum_file, FQM_w, FQM_I, DCL_I, ACL_I)
