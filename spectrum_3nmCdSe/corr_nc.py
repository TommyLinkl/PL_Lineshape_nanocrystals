import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.fft import fft, ifft
from scipy.integrate import simpson

T = 16
KB = 1.380649e-23   # J/K, SI units
BETA = 1./ (KB * T)
PI = 3.14159265359
HBAR = 1.054571817e-34  # Js, SI units
JTOEV = 6.24150907446076e18 

domega = 0.01e12 # in Hz
dt = 1e-15
t_range = np.arange(0.0, 5e-11+dt, dt)
# t_range = np.arange(0.0, 2e-12+dt, dt)

# erf param
cutoff = 4e-11
slope = 0.005

omega_10 = 6.154592518543e12  # in Hz

def read_w(filename):
    w = np.zeros(0)
    f = open(filename, "r")
    while True:
        line = f.readline()
        if not line:
            break
        if line[0]!= "#":
            w = np.append(w, float(line.split()[1]))
    f.close()
    w *= 1e12 * (2*PI) 
    return w

def read_Vklq(filename):
    V00 = np.zeros(0)
    V11 = np.zeros(0)

    f = open(filename, "r")
    while True:
        line = f.readline()
        if not line:
            break
        if line[0]!= "#":
            nums = line.split()
            if (int(nums[1])==0) and (int(nums[2])==0):
                V00 = np.append(V00, float(nums[3]))
            if (int(nums[1])==1) and (int(nums[2])==1):
                V11 = np.append(V11, float(nums[3]))
    f.close()
    return (V00, V11)

def calc_corr_t(t_range):
    corr_t = np.zeros(0)
    for t in t_range: 
        # print(t)
        func = 0
        for i in range(len(w)):
            if (i>5):
                func += (V11[i]-V00[i])**2 / (HBAR**2 * BETA * w[i]**4) * (np.cos(w[i] * t) - 1)
        corr_t = np.append(corr_t, np.exp(func))
    return corr_t

def write_corr(filename, t_range, corr_t):
    f = open(filename, "w")
    f.write("# time   C_Re\n")
    for i in range(len(t_range)):
        f.write("%.15f       %.15f\n" % (t_range[i], corr_t[i]))
    f.close()
    return

def multiply_erf(cutoff, slope, corr_t, dt):
    corr_erf = np.zeros(0)
    for i in range(len(t_range)):
        corr_erf = np.append(corr_erf, corr_t[i]* (-0.5*math.erf(slope*(float(i)-cutoff/dt))+0.5))
    return corr_erf

def spectrum(omega_10, t_range, corr_Re, corr_Im, domega):
    I = np.zeros(0)
    w_range = np.arange(-30e12, 30e12, domega)  # in Hz

    for w in w_range:
        integrand = np.zeros(0)
        for i in range(len(t_range)):
            integrand = np.append(integrand, 1/PI * np.cos((omega_10 - w) * t_range[i]) * corr_Re[i] + 1/PI * np.sin((omega_10 - w) * t_range[i]) * corr_Im[i])
        I_w = simpson(integrand, x=t_range) / (HBAR * JTOEV)
        I = np.append(I, I_w)
    return (w_range, I)

def write_spectrum(filename, w_range, spectrum):
    f = open(filename, "w")
    f.write("# Frequency w         I(w)\n")
    for i in range(len(w_range)):
        f.write("%f       %.15f\n" % (w_range[i], spectrum[i]))
    f.close()
    return

# main

# Reading frequency and coupling matrix element data
w = read_w("w.dat")
(V00, V11) = read_Vklq("Vklq-diabatic.dat")

# Calculating the correlation function, and writing to file
corr_t = calc_corr_t(t_range)
print(corr_t)
write_corr("3nmCdSe_T=16_shorttime_corr.dat", t_range, corr_t)

# erf correction to the correlation function, and writing to file
corr_erf = multiply_erf(cutoff, slope, corr_t-np.mean(corr_t[-200:-1]), dt)   
write_corr("3nmCdSe_T=16_shorttime_corr_erf.dat", t_range, corr_erf)

# FT to obtain the spectrum from the correlation function
corr_Im = np.zeros(len(corr_erf))
(spectrum_w, spectrum_I) = spectrum(omega_10, t_range, corr_erf, corr_Im, domega)
write_spectrum("3nmCdSe_T=16_shorttime_spectrum.dat", spectrum_w, spectrum_I)
