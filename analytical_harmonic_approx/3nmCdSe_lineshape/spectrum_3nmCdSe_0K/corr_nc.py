import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.fft import fft, ifft
from scipy.integrate import simpson

T = 0                            # unit: K 
KB = 1.380649e-23                # unit: J/K, SI units
# BETA = 1./ (KB * T)   
PI = 3.14159265359
HBAR = 1.054571817e-34           # Unit: Js, SI units
JTOEV = 6.24150907446076e18 
AUTOEV = 27.211396641308 
AUTOJ = 4.359744e-18

dt = 1e-15
t_range = np.arange(0.0, 30e-12+dt, dt)
domega = 0.5e12    # 1e12                 # in Hz

# erf param
cutoff_1 = 5e-12
cutoff_2 = 15e-12
cutoff_3 = 25e-12
slope = 0.01

excitonic_e = 0.080794824219 * AUTOJ     # in J
# omega_10 = 6.154592518543e12     # in Hz

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
            Re += 1 / (2*HBAR) * delta_range[i]**2 * omega_range[i] * (np.cos(omega_range[i] * t) - 1)
            Im += 1 / (2*HBAR) * delta_range[i]**2 * omega_range[i] * (-np.sin(omega_range[i] * t))
        F_t_Re = np.append(F_t_Re, np.exp(Re) * np.cos(Im))
        F_t_Im = np.append(F_t_Im, np.exp(Re) * np.sin(Im))
    return (F_t_Re, F_t_Im)

'''
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
'''

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

def spectrum(E_shift, t_range, F_Re, F_Im, domega):
    I = np.zeros(0)
    w_range = np.arange(2000e12, 4000e12, domega)  # in Hz
    print("E_shift: %f J, %f eV" % (E_shift, E_shift*JTOEV))
    print("w_shift in Hz: %f" % (E_shift / HBAR))

    for w in w_range:
        integrand = np.zeros(0)
        for i in range(len(t_range)):
            integrand = np.append(integrand, 1/PI * np.cos((w-E_shift/HBAR) * t_range[i]) * F_Re[i] + 1/PI * np.sin((w-E_shift/HBAR) * t_range[i]) * F_Im[i])
        I_w = simpson(integrand, x=t_range)
        I = np.append(I, I_w)
    return (w_range* HBAR * JTOEV, I)    # return w_range in eV

def write_spectrum(filename, w_range, spectrum):
    f = open(filename, "w")
    f.write("# Frequency w         I(w)\n")
    for i in range(len(w_range)):
        f.write("%f       %.20f\n" % (w_range[i], spectrum[i]))
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
(F_t_Re_E17, F_t_Im_E17) = calc_dephasing_F_E17(t_range, delta_range, w)
print("DONE calc_dephasing_F_E17")

'''
(F_t_Re_MC, F_t_Im_MC) = calc_dephasing_F_MC(t_range, delta_range, w)
print("DONE calc_dephasing_F_MC")

(DC_t_Re_E17, DC_t_Im_E17) = calc_delta_corr_E17(t_range, delta_range, w)
DC_t_Re_E17 *= JTOEV**2
DC_t_Im_E17 *= JTOEV**2
print("DONE calc_delta_corr_E17")

(DC_t_Re_MC, DC_t_Im_MC) = calc_delta_corr_MC(t_range, delta_range, w)
DC_t_Re_MC *= JTOEV**2
DC_t_Im_MC *= JTOEV**2
print("DONE calc_delta_corr_MC")
write_dephasing_F("3nmCdSe_F.dat", t_range, F_t_Re_MC, F_t_Im_MC, F_t_Re_E17, F_t_Im_E17)
write_delta_corr("3nmCdSe_Delta_Corr.dat", t_range, DC_t_Re_MC, DC_t_Im_MC, DC_t_Re_E17, DC_t_Im_E17)
'''

# erf correction to the dephasing functions
F_t_Re_E17_erf_1 = multiply_erf(cutoff_1, slope, F_t_Re_E17-np.mean(F_t_Re_E17[-1000:-1]), dt)
F_t_Im_E17_erf_1 = multiply_erf(cutoff_1, slope, F_t_Im_E17-np.mean(F_t_Im_E17[-1000:-1]), dt)
F_t_Re_E17_erf_2 = multiply_erf(cutoff_2, slope, F_t_Re_E17-np.mean(F_t_Re_E17[-1000:-1]), dt)
F_t_Im_E17_erf_2 = multiply_erf(cutoff_2, slope, F_t_Im_E17-np.mean(F_t_Im_E17[-1000:-1]), dt)
F_t_Re_E17_erf_3 = multiply_erf(cutoff_3, slope, F_t_Re_E17-np.mean(F_t_Re_E17[-1000:-1]), dt)
F_t_Im_E17_erf_3 = multiply_erf(cutoff_3, slope, F_t_Im_E17-np.mean(F_t_Im_E17[-1000:-1]), dt)
# F_t_Re_MC_erf = multiply_erf(cutoff, slope, F_t_Re_MC-np.mean(F_t_Re_MC[-1000:-1]), dt)
# F_t_Im_MC_erf = multiply_erf(cutoff, slope, F_t_Im_MC-np.mean(F_t_Im_MC[-1000:-1]), dt)
# write_corr("3nmCdSe_T=16_shorttime_corr_erf.dat", t_range, corr_erf)
print("DONE erf_multiplication")

# FT the dephasing functions to obtain the spectrum from the correlation function
delta_E_E17 = excitonic_e + reorg_E(delta_range, w)    # in J
delta_E_MC = excitonic_e                               # in J
(spectrum_w_E17_cutoff1, spectrum_I_E17_cutoff1) = spectrum(delta_E_E17, t_range, F_t_Re_E17_erf_1, F_t_Im_E17_erf_1, domega)
(spectrum_w_E17_cutoff2, spectrum_I_E17_cutoff2) = spectrum(delta_E_E17, t_range, F_t_Re_E17_erf_2, F_t_Im_E17_erf_2, domega)
(spectrum_w_E17_cutoff3, spectrum_I_E17_cutoff3) = spectrum(delta_E_E17, t_range, F_t_Re_E17_erf_3, F_t_Im_E17_erf_3, domega)
print("DONE spectrum_E17")
# (spectrum_w_MC, spectrum_I_MC) = spectrum(delta_E_MC, t_range, F_t_Re_MC_erf, F_t_Im_MC_erf, domega)
# print("DONE spectrum_MC")
write_spectrum("3nmCdSe_spectrum_E17_cutoff1.dat", spectrum_w_E17_cutoff1, spectrum_I_E17_cutoff1)
write_spectrum("3nmCdSe_spectrum_E17_cutoff2.dat", spectrum_w_E17_cutoff2, spectrum_I_E17_cutoff2)
write_spectrum("3nmCdSe_spectrum_E17_cutoff3.dat", spectrum_w_E17_cutoff3, spectrum_I_E17_cutoff3)
# write_spectrum("3nmCdSe_spectrum_MC.dat", spectrum_w_MC, spectrum_I_MC)

# plotting
fig, axs = plt.subplots(3,3, figsize=(11,11))
fig.suptitle(r"The dephasing function and spectra at different cutoff points, E17 model, 3nmCdSe at T=300K")

for i in range(3): 
    for j in range(3): 
        axs[i,j].grid()

axs[0,0].set(xlabel="Time (s)", ylabel=r"$\langle F(t) \rangle$ (unitless)", title="E17 dephasing function cutoff at 5 ps", yscale="log")
axs[0,0].plot(t_range, F_t_Re_E17, label=r"$\langle F(t) \rangle _{Re} $")
axs[0,0].plot(t_range, F_t_Im_E17, ":", label=r"$\langle F(t) \rangle _{Im}$")
axs[0,0].plot(t_range, F_t_Re_E17_erf_1, label=r"$\langle F(t) \rangle _{Re, shifted} \cdot erf$")
axs[0,0].plot(t_range, F_t_Im_E17_erf_1, ":", label=r"$\langle F(t) \rangle _{Im, shifted} \cdot erf$")

axs[1,0].set(xlabel="Time (s)", ylabel=r"$\langle F(t) \rangle$ (unitless)", title="E17 dephasing function cutoff at 15 ps", yscale="log")
axs[1,0].plot(t_range, F_t_Re_E17, label=r"$\langle F(t) \rangle _{Re} $")
axs[1,0].plot(t_range, F_t_Im_E17, ":", label=r"$\langle F(t) \rangle _{Im}$")
axs[1,0].plot(t_range, F_t_Re_E17_erf_2, label=r"$\langle F(t) \rangle _{Re, shifted} \cdot erf$")
axs[1,0].plot(t_range, F_t_Im_E17_erf_2, ":", label=r"$\langle F(t) \rangle _{Im, shifted} \cdot erf$")

axs[2,0].set(xlabel="Time (s)", ylabel=r"$\langle F(t) \rangle$ (unitless)", title="E17 dephasing function cutoff at 25 ps", yscale="log")
axs[2,0].plot(t_range, F_t_Re_E17, label=r"$\langle F(t) \rangle _{Re} $")
axs[2,0].plot(t_range, F_t_Im_E17, ":", label=r"$\langle F(t) \rangle _{Im}$")
axs[2,0].plot(t_range, F_t_Re_E17_erf_3, label=r"$\langle F(t) \rangle _{Re, shifted} \cdot erf$")
axs[2,0].plot(t_range, F_t_Im_E17_erf_3, ":", label=r"$\langle F(t) \rangle _{Im, shifted} \cdot erf$")

# spectra
axs[0,1].set(xlabel="Energey [meV]", title="Spectra", xlim=(2000, 2400))
axs[0,1].plot(spectrum_w_E17_cutoff1*1000, spectrum_I_E17_cutoff1, linewidth=3, label="erf cutoff at 5ps")
axs[0,1].plot(spectrum_w_E17_cutoff2*1000, spectrum_I_E17_cutoff2, linewidth=2, label="erf cutoff at 15ps")
axs[0,1].plot(spectrum_w_E17_cutoff3*1000, spectrum_I_E17_cutoff3, linewidth=1, label="erf cutoff at 25ps")

axs[1,1].set(xlabel="Energey [meV]", title="Spectra Diff", xlim=(2000, 2400))
axs[1,1].plot(spectrum_w_E17_cutoff1*1000, spectrum_I_E17_cutoff2-spectrum_I_E17_cutoff1, linewidth=3, label="15ps spectrum - 5ps spectrum")
axs[1,1].plot(spectrum_w_E17_cutoff1*1000, spectrum_I_E17_cutoff3-spectrum_I_E17_cutoff1, linewidth=1, label="25ps spectrum - 5ps spectrum")

'''
# Delta-Delta
axs[0,2].set(xlabel="Time (s)", ylabel=r"$\langle \Delta (t) \Delta (0) \rangle$, ($eV^2$)", title=r"E17 $\Delta - \Delta$ correlation function")
axs[0,2].plot(t_range, DC_t_Re_E17, label=r"$C_{\Delta \Delta}(t)_{Re} $")
axs[0,2].plot(t_range, DC_t_Im_E17, ":", label=r"$C_{\Delta \Delta}(t)_{Im} $")

axs[1,2].set(xlabel="Time (s)", ylabel=r"$\langle \Delta (t) \Delta (0) \rangle$, ($eV^2$)", title=r"E17 $\Delta - \Delta$ correlation function", yscale="log")
axs[1,2].plot(t_range, DC_t_Re_E17, label=r"$C_{\Delta \Delta}(t)_{Re} $")
axs[1,2].plot(t_range, DC_t_Im_E17, ":", label=r"$C_{\Delta \Delta}(t)_{Im} $")
'''

for i in range(3): 
    for j in range(3): 
        axs[i,j].legend()
fig.tight_layout()
fig.savefig("nanocrystal_0K.png")
