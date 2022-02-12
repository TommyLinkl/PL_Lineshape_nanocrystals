import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.fft import fft, ifft
from scipy.integrate import simpson

T = 1
BETA = 1./T
PI = 3.14159265359

dt = 0.1
t_range = np.arange(0.0, 50.0+dt, dt)

omega_10 = 0.
cutoff = 40.
slope = 0.1

corr_file_MC = "MC_oscillating_T=1_corr.dat"
I_file_MC = "MC_oscillating_T=1_I.dat"
corr_file_E17 = "E17_nonzero_T=1_corr.dat"
I_file_E17 = "E17_nonzero_T=1_I.dat"

def coth(x):
    return (np.exp(x) + np.exp(-x))/(np.exp(x) - np.exp(-x))

def J(omega):
    return omega * np.exp(-omega)

def FQM_Re(omega, time):
    return 0.5 * J(omega) * omega * coth(BETA * omega / 2) * (np.cos(omega*time) - 1)

def Int_FQM_Re(time, omega_range): 
    y_range = np.zeros(0)
    for omega in omega_range:
        y_range = np.append(y_range, FQM_Re(omega, time))
    return simpson(y_range, omega_range)

# MC
def FQM_Im_MC(omega, time):
    return 0.5 * J(omega) * omega * (-np.sin(omega * time) + omega*time)

def Int_FQM_Im_MC(time, omega_range): 
    y_range = np.zeros(0)
    for omega in omega_range:
        y_range = np.append(y_range, FQM_Im_MC(omega, time))
    return simpson(y_range, omega_range)

def corr_FQM_Re_MC(t, omega_range): 
    return np.exp(Int_FQM_Re(t, omega_range)) * np.cos(Int_FQM_Im_MC(t, omega_range))

def corr_FQM_Im_MC(t, omega_range):
    return np.exp(Int_FQM_Re(t, omega_range)) * np.sin(Int_FQM_Im_MC(t, omega_range))

# E17
def FQM_Im_E17(omega, time):
    return 0.5 * J(omega) * omega * (-np.sin(omega * time))

def Int_FQM_Im_E17(time, omega_range): 
    y_range = np.zeros(0)
    for omega in omega_range:
        y_range = np.append(y_range, FQM_Im_E17(omega, time))
    return simpson(y_range, omega_range)

def corr_FQM_Re_E17(t, omega_range): 
    return np.exp(Int_FQM_Re(t, omega_range)) * np.cos(Int_FQM_Im_E17(t, omega_range))

def corr_FQM_Im_E17(t, omega_range):
    return np.exp(Int_FQM_Re(t, omega_range)) * np.sin(Int_FQM_Im_E17(t, omega_range))

def write_corr(filename, t_range, c_FQM_Re, c_FQM_Im, c_FQM_Re_erf, c_FQM_Im_erf): 
    f = open(filename, "w")
    f.write("# time   C_FQM_Re    C_FQM_Im    C_FQM_Re_erf    C_FQM_Im_erf \n")
    for i in range(len(t_range)): 
        f.write("%f       %f       %f       %f      %f \n" % (t_range[i], c_FQM_Re[i], c_FQM_Im[i], c_FQM_Re_erf[i], c_FQM_Im_erf[i]))
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

def write_spectrum(filename, FQM_w, FQM_I, FQM_I_erf):
    f = open(filename, "w")
    f.write("# FQM_w,     FQM_I,    FQM_I_of_C_after_erf  \n")
    for i in range(len(FQM_w)):
        f.write("%f       %f       %f \n" % (FQM_w[i], FQM_I[i], FQM_I_erf[i]))
    f.close()
    return


# main
N = 1000
domega = 20. / N
omega_range = np.arange(0.0+domega, 20.+domega, domega)

# MC
c_FQM_Re_MC = np.zeros(0)
c_FQM_Im_MC = np.zeros(0)
for t in t_range:
    c_FQM_Re_MC = np.append(c_FQM_Re_MC, corr_FQM_Re_MC(t, omega_range))
    c_FQM_Im_MC = np.append(c_FQM_Im_MC, corr_FQM_Im_MC(t, omega_range))
c_FQM_Re_MC_erf = multiply_erf(cutoff, slope, c_FQM_Re_MC, dt)
c_FQM_Im_MC_erf = multiply_erf(cutoff, slope, c_FQM_Im_MC, dt)

write_corr(corr_file_MC, t_range, c_FQM_Re_MC, c_FQM_Im_MC, c_FQM_Re_MC_erf, c_FQM_Im_MC_erf)
(w_range_MC, I_MC) = spectrum(omega_10, t_range, c_FQM_Re_MC, c_FQM_Im_MC, domega)
(w_range_MC_erf, I_MC_erf) = spectrum(omega_10, t_range, c_FQM_Re_MC_erf, c_FQM_Im_MC_erf, domega)

write_spectrum(I_file_MC, w_range_MC, I_MC, I_MC_erf)

# E17
c_FQM_Re_E17 = np.zeros(0)
c_FQM_Im_E17 = np.zeros(0)
for t in t_range:
    c_FQM_Re_E17 = np.append(c_FQM_Re_E17, corr_FQM_Re_E17(t, omega_range))
    c_FQM_Im_E17 = np.append(c_FQM_Im_E17, corr_FQM_Im_E17(t, omega_range))
c_FQM_Re_E17_erf = c_FQM_Re_E17 - np.mean(c_FQM_Re_E17[-100:])
c_FQM_Im_E17_erf = c_FQM_Im_E17 - np.mean(c_FQM_Im_E17[-100:])

write_corr(corr_file_E17, t_range, c_FQM_Re_E17, c_FQM_Im_E17, c_FQM_Re_E17_erf, c_FQM_Im_E17_erf)
(w_range_E17, I_E17) = spectrum(omega_10, t_range, c_FQM_Re_E17, c_FQM_Im_E17, domega)
(w_range_E17_erf, I_E17_erf) = spectrum(omega_10, t_range, c_FQM_Re_E17_erf, c_FQM_Im_E17_erf, domega)

write_spectrum(I_file_E17, w_range_E17, I_E17, I_E17_erf)


# plotting
fig, axs = plt.subplots(2,3, figsize=(12,8))
fig.suptitle(r"The MC and E17 cases for T=1, $\langle F(t) \rangle$ & its FT")

axs[0,0].grid()
axs[0,1].grid()
axs[0,2].grid()
axs[1,0].grid()
axs[1,1].grid()
axs[1,2].grid()
axs[0,0].set(xlabel="Time", ylabel=r"$\langle F(t) \rangle$", title="E17 dephasing function", ylim=(-0.5, 1.))
axs[1,0].set(xlabel="Time", ylabel=r"$\langle F(t) \rangle$", yscale="log")
axs[0,1].set(xlabel="Time", ylabel=r"$\langle F(t) \rangle$", title="MC dephasing function", ylim=(-0.5, 1.))
axs[1,1].set(xlabel="Time", ylabel=r"$\langle F(t) \rangle$", yscale="log")
axs[0,2].set(xlabel="Frequency", title=r"FT of $\langle F(t) \rangle$")
axs[1,2].set(xlabel="Frequency", xlim=(-8,8), ylim=(-0.1, 0.3))
 
axs[0,0].plot(t_range, c_FQM_Re_E17, label=r"$\langle F(t) \rangle _{Re} $")
axs[0,0].plot(t_range, c_FQM_Im_E17, ":", label=r"$\langle F(t) \rangle _{Im}$")
axs[0,0].plot(t_range, c_FQM_Re_E17_erf, label=r"$\langle F(t) \rangle _{Re} \cdot erf$")
axs[0,0].plot(t_range, c_FQM_Im_E17_erf, ":", label=r"$\langle F(t) \rangle _{Im} \cdot erf$")

axs[1,0].plot(t_range, c_FQM_Re_E17, color="C0", label=r"$\langle F(t) \rangle _{Re} $")
axs[1,0].plot(t_range, c_FQM_Re_E17_erf, color="C2", label=r"$\langle F(t) \rangle _{Re} \cdot erf$")

axs[0,1].plot(t_range, c_FQM_Re_MC, label=r"$\langle F(t) \rangle _{Re} $")
axs[0,1].plot(t_range, c_FQM_Im_MC, ":", label=r"$\langle F(t) \rangle _{Im}$")
axs[0,1].plot(t_range, c_FQM_Re_MC_erf, label=r"$\langle F(t) \rangle _{Re} \cdot erf$")
axs[0,1].plot(t_range, c_FQM_Im_MC_erf, ":", label=r"$\langle F(t) \rangle _{Im} \cdot erf$")

axs[1,1].plot(t_range, c_FQM_Re_MC, color="C0", label=r"$\langle F(t) \rangle _{Re} $")
axs[1,1].plot(t_range, c_FQM_Re_MC_erf, color="C2", label=r"$\langle F(t) \rangle _{Re} \cdot erf$")

axs[0,2].plot(w_range_MC_erf, I_MC_erf, label=r"FT of $\langle F(t) \rangle _{MC} \cdot erf$")
axs[0,2].plot(w_range_E17_erf, I_E17_erf, label=r"FT of $\langle F(t) \rangle _{E17} \cdot erf$")

axs[1,2].plot(w_range_MC_erf+3., I_MC_erf, label=r"FT of $\{ \langle F(t) \rangle _{MC} \cdot erf + \lambda \}$, $\lambda=3$")
axs[1,2].plot(w_range_E17_erf, I_E17_erf, label=r"FT of $\{ \langle F(t) \rangle _{E17} \cdot erf \}$")

axs[0,0].legend()
axs[0,1].legend()
axs[1,0].legend()
axs[1,1].legend()
axs[0,2].legend()
axs[1,2].legend()
fig.tight_layout()
fig.savefig("MC_E17_FT.png")

