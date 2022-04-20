import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.fft import fft, ifft
from scipy.integrate import simpson

# T = 1
# BETA = 1./T
PI = 3.14159265359

omega_10 = 0.
domega = 0.01

corr_file = "corr.dat"
rate_file = "k.dat"

def read_corr(filename): 
    f = open(filename, "r")
    t_range = np.zeros(0)
    c_re = np.zeros(0)
    c_im = np.zeros(0)
    while True: 
        line = f.readline()
        if not line:
            break
        if (line[0]!="#") and (len(line)!=0):
            t_range = np.append(t_range, float(line.split()[0]))
            c_re = np.append(c_re, float(line.split()[1]))
            c_im = np.append(c_im, float(line.split()[2]))
    f.close()
    return (t_range, c_re, c_im)

def multiply_erf(cutoff, slope, t_range, corr_Re, corr_Im):
    corr_Re_erf = np.zeros(0)
    corr_Im_erf = np.zeros(0)
    dt = t_range[1] - t_range[0]
    for i in range(len(t_range)):
        corr_Re_erf = np.append(corr_Re_erf, corr_Re[i]* (-0.5*math.erf(slope*(float(i)-cutoff/dt))+0.5))
        corr_Im_erf = np.append(corr_Im_erf, corr_Im[i]* (-0.5*math.erf(slope*(float(i)-cutoff/dt))+0.5))
    return (corr_Re_erf, corr_Im_erf)

def rate(omega_10, t_range, corr_Re, corr_Im, domega):
    k = np.zeros(0)
    w_range = np.arange(0.0, 20.0, domega)

    for w in w_range:
        integrand = np.zeros(0)
        for i in range(len(t_range)):
            integrand = np.append(integrand, 0.02 * np.cos(w*t_range[i]) * corr_Re[i] - 0.02 * np.sin(w*t_range[i]) * corr_Im[i])
        k_w = simpson(integrand, x=t_range)
        k = np.append(k, k_w)
    return (w_range, k)

def write_rates(filename, w_range, k):
    f = open(filename, "w")
    f.write("# w,    rate k\n")
    for i in range(len(w_range)):
        f.write("%f       %f\n" % (w_range[i], k[i]))
    f.close()
    return

# main

(t_range, corr_Re, corr_Im) = read_corr(corr_file)
# multiply_erf
(w_range, k) = rate(omega_10, t_range, corr_Re, corr_Im, domega)
write_rates(rate_file, w_range, k)
