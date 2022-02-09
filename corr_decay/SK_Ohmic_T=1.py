import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.fft import fft, ifft
from scipy.integrate import simpson

T = 1
BETA = 1./T
PI = 3.14159265359

dt = 0.1
t_range = np.arange(0.0, 40.0+dt, dt)

# corr_file = "corr_Ohmic_T=1.dat"

def coth(x):
    return (np.exp(x) + np.exp(-x))/(np.exp(x) - np.exp(-x))

def J(omega):
    return omega**(-2) * np.exp(-omega)

def FQM_Re(omega, time):
    return 0.5 * J(omega) * omega * coth(BETA * omega / 2) * (np.cos(omega*time) - 1)

def Int_FQM_Re(time, omega_range): 
    y_range = np.zeros(0)
    for omega in omega_range:
        y_range = np.append(y_range, FQM_Re(omega, time))
    return simpson(y_range, omega_range)

def FQM_Im(omega, time):
    return 0.5 * J(omega) * omega * (-np.sin(omega * time) + omega*time)

def Int_FQM_Im(time, omega_range): 
    y_range = np.zeros(0)
    for omega in omega_range:
        y_range = np.append(y_range, FQM_Im(omega, time))
    return simpson(y_range, omega_range)

def corr_FQM_Re(t, omega_range): 
    return np.exp(Int_FQM_Re(t, omega_range)) * np.cos(Int_FQM_Im(t, omega_range))

def corr_FQM_Im(t, omega_range):
    return np.exp(Int_FQM_Re(t, omega_range)) * np.sin(Int_FQM_Im(t, omega_range))




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

fig, axs = plt.subplots(2)
fig.suptitle("Convergence for MC's equivalent of Skinner's Ohmic, T=1")

axs[0].grid()
axs[1].grid()
axs[0].set(xlabel="Time", ylabel="C(t)")
axs[1].set(xlabel="Time", ylabel="C(t)", yscale="log")
# Calc FQM, DCL, ACL analytical correlation functions, and write to file
for N in np.array([2, 10, 100, 1000, 10000]): 
    domega = 20. / N
    omega_range = np.arange(0.0+domega, 20.+domega, domega)

    c_FQM_Re = np.zeros(0)
    c_FQM_Im = np.zeros(0)
    for t in t_range: 
        c_FQM_Re = np.append(c_FQM_Re, corr_FQM_Re(t, omega_range))
        c_FQM_Im = np.append(c_FQM_Im, corr_FQM_Im(t, omega_range))

    axs[0].plot(t_range, c_FQM_Re, label="NB=%.1f" % N)
    axs[1].plot(t_range, c_FQM_Re, label="NB=%.1f" % N)


axs[0].legend()
axs[1].legend()
fig.tight_layout()
fig.savefig("Skinner_Ohmic_T=1.png")
