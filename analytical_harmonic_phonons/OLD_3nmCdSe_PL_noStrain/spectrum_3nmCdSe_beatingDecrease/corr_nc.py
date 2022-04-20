import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.fft import fft, ifft
from scipy.integrate import simpson

T = 1                            # unit: K 
KB = 1.380649e-23                # unit: J/K, SI units
BETA = 1./ (KB * T)   
PI = 3.14159265359
HBAR = 1.054571817e-34           # Unit: Js, SI units
JTOEV = 6.24150907446076e18 
AUTOEV = 27.211396641308 
AUTOJ = 4.359744e-18

dt1 = 1e-14
dt2 = 1e-12
dt3 = 1e-9
t_range = np.arange(0.0, 1e-12, dt1)
t_range = np.append(t_range, np.arange(1e-12, 1e-9, dt2))
t_range = np.append(t_range, np.arange(1e-9, 1e-7+dt3, dt3))
domega = 0.01e12                 # in Hz

# erf param
cutoff = 0.8e-11
slope = 0.01

excitonic_e = 0.080794824219     # in Hartree

def coth(x):
    return (np.exp(x) + np.exp(-x))/(np.exp(x) - np.exp(-x))

def func_n(x): 
    return 1 / (np.exp(BETA * HBAR * x) - 1)

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
    # print(w)
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
    # print(V00)
    # print(V11)
    return (V00, V11)

def reorg_E(delta_range, omega_range):
    E = 0
    for i in range(len(omega_range)):
        E += 0.5 * omega_range[i]**2 * delta_range[i]**2
    print("reorg_E = %f J, %f eV \n" % (E, E*JTOEV))
    return E   # In units of J

def calc_dephasing_F_E17(t_range, delta_range, omega_range): 
    F_t_Re = np.zeros(0)
    F_t_Im = np.zeros(0)
    for t in t_range: 
        # print(t)
        Re = 0
        Im = 0
        for i in range(len(omega_range)):
            Re += 1 / (2*HBAR) * delta_range[i]**2 * omega_range[i] * coth(BETA * HBAR * omega_range[i] / 2) * (np.cos(omega_range[i] * t) - 1)
            Im += 1 / (2*HBAR) * delta_range[i]**2 * omega_range[i] * (-np.sin(omega_range[i] * t))
        F_t_Re = np.append(F_t_Re, np.exp(Re) * np.cos(Im))
        F_t_Im = np.append(F_t_Im, np.exp(Re) * np.sin(Im))
    return (F_t_Re, F_t_Im)

def calc_dephasing_F_MC(t_range, delta_range, omega_range): 
    F_t_Re = np.zeros(0)
    F_t_Im = np.zeros(0)
    for t in t_range: 
        # print(t)
        Re = 0
        Im = 0
        for i in range(len(omega_range)):
            Re += 1 / (2*HBAR) * delta_range[i]**2 * omega_range[i] * coth(BETA * HBAR * omega_range[i] / 2) * (np.cos(omega_range[i] * t) - 1)
            Im += 1 / (2*HBAR) * delta_range[i]**2 * omega_range[i] * (-np.sin(omega_range[i] * t) + omega_range[i] * t)
        F_t_Re = np.append(F_t_Re, np.exp(Re) * np.cos(Im))
        F_t_Im = np.append(F_t_Im, np.exp(Re) * np.sin(Im))
    return (F_t_Re, F_t_Im)

def calc_delta_corr_E17(t_range, delta_range, omega_range): 
    (C_t_Re, C_t_Im) = calc_delta_corr_MC(t_range, delta_range, omega_range)
    C_t_Re += reorg_E(delta_range, omega_range)**2
    return (C_t_Re, C_t_Im)

def calc_delta_corr_MC(t_range, delta_range, omega_range): 
    C_t_Re = np.zeros(0)
    C_t_Im = np.zeros(0)
    for t in t_range: 
        Re = 0.
        Im = 0.
        for i in range(len(omega_range)): 
            w = omega_range[i]
            Re += HBAR / 2 * delta_range[i]**2 * w**3 * ((func_n(w) + 1)*np.cos(w*t) + func_n(w)*np.cos(w*t))
            Im += HBAR / 2 * delta_range[i]**2 * w**3 * ((func_n(w) + 1)*np.sin(-w*t) + func_n(w)*np.sin(w*t))
        C_t_Re = np.append(C_t_Re, Re)
        C_t_Im = np.append(C_t_Im, Im)
    return (C_t_Re, C_t_Im)

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

def multiply_erf(cutoff, slope, corr_t, dt):
    corr_erf = np.zeros(0)
    for i in range(len(t_range)):
        corr_erf = np.append(corr_erf, corr_t[i]* (-0.5*math.erf(slope*(float(i)-cutoff/dt))+0.5))
    return corr_erf

def spectrum(w_shift, t_range, corr_Re, corr_Im, domega):
    I = np.zeros(0)
    w_range = np.arange(-30e12, 30e12, domega)  # in Hz
    # reads in w_shift in Hz
    
    for w in w_range:
        integrand = np.zeros(0)
        for i in range(len(t_range)):
            integrand = np.append(integrand, 1/PI * np.cos((w_shift - w) * t_range[i]) * corr_Re[i] + 1/PI * np.sin((w_shift - w) * t_range[i]) * corr_Im[i])
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
w = w[6:]
(V00, V11) = read_Vklq("Vklq-diabatic.dat")
V00 = V00[6:]
V11 = V11[6:]
delta_range = (V00) / w**2
print("DONE reading")

# Calculating the correlation function, and writing to file
(DC_t_Re_MC, DC_t_Im_MC) = calc_delta_corr_MC(t_range, delta_range, w)
DC_t_Re_MC *= JTOEV**2
DC_t_Im_MC *= JTOEV**2
print("DONE calc_delta_corr_MC")


w = w[:500]
V00 = V00[:500]
delta_range = (V00) / w**2
print("DONE reading")

(DC_t_Re_MC_less, DC_t_Im_MC_less) = calc_delta_corr_MC(t_range, delta_range, w)
DC_t_Re_MC *= JTOEV**2
DC_t_Im_MC *= JTOEV**2
print("DONE calc_delta_corr_MC")

# plotting
plt.plot(t_range, DC_t_Re_MC, label=r"$C_{\Delta \Delta}(t)_{Re}$, 1332 modes")
plt.plot(t_range, DC_t_Re_MC_less, label=r"$C_{\Delta \Delta}(t)_{Re}$, 500 modes")
plt.xlabel("Time (s)")
plt.ylabel(r"$\langle \Delta (t) \Delta (0) \rangle$, ($eV^2$)")
plt.title(r"MC $\Delta - \Delta$ correlation function - semilog")
plt.yscale("log")
plt.grid()
plt.legend()

plt.tight_layout()
plt.savefig("nanocrystal_beating.png")
