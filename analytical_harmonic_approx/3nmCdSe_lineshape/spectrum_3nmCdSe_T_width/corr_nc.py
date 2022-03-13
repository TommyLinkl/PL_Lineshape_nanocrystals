import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.fft import fft, ifft
from scipy.integrate import simpson

KB = 1.380649e-23                # unit: J/K, SI units
PI = 3.14159265359
HBAR = 1.054571817e-34           # Unit: Js, SI units
JTOEV = 6.24150907446076e18 
AUTOEV = 27.211396641308 
AUTOJ = 4.359744e-18

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

def read_Vklq(filename, nPhonon, nExc):
    coupling = np.zeros((nPhonon, nExc))
    f = open(filename, "r")
    while True:
        line = f.readline()
        if not line:
            break
        if (line[0]!= "#") and (int(line.split()[1])==int(line.split()[2])):
            coupling[int(line.split()[0]), int(line.split()[1])] = float(line.split()[3])
    f.close()
    return (coupling)

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

def spectrum(E_shift, t_range, F_Re, F_Im, w_range):
    # E_shift  in J, w_range in eV
    I = np.zeros(0)
    print("E_shift: %f J, %f eV" % (E_shift, E_shift*JTOEV))
    print("w_shift in Hz: %f" % (E_shift / HBAR))

    for w in w_range:
        w *= 1/(JTOEV * HBAR)
        integrand = np.zeros(0)
        for i in range(len(t_range)):
            integrand = np.append(integrand, 1/PI * np.cos((w-E_shift/HBAR) * t_range[i]) * F_Re[i] + 1/PI * np.sin((w-E_shift/HBAR) * t_range[i]) * F_Im[i])
        I_w = simpson(integrand, x=t_range)
        I = np.append(I, I_w)
    return (w_range, I)    # return w_range in eV

def write_spectrum(filename, w_range, spectrum):
    f = open(filename, "w")
    f.write("# Frequency w         I(w)\n")
    for i in range(len(w_range)):
        f.write("%f       %.20f\n" % (w_range[i], spectrum[i]))
    f.close()
    return

def calc_FWHM(w, I_w): 
    half_max = I_w.max()/2
    max_index = np.argmax(I_w)
    for r_index in range(max_index, len(w), 1): 
        if I_w[r_index]<=half_max:
            break
    for l_index in range(max_index, -1, -1): 
        if I_w[l_index]<=half_max: 
            break
    return (w[r_index]-w[l_index], l_index, r_index)

def write_FWHM(filename, w, I_w, l_index, r_index):
    f = open(filename, "w")
    f.write("# FWHM width       points_w      points_I\n")
    f.write("%f        %f       %.20f\n" % (w[r_index]-w[l_index], w[l_index], I_w[l_index]))
    f.write("%f        %f       %.20f\n" % (w[r_index]-w[l_index], w[r_index], I_w[r_index]))
    f.close()
    return

def calc_spectrum_E17(T, t_range, w_range, cutoff, slope, exc_E, phonon_filename, coupling_filename, I_file, FWHM_file):
    # give T in K, t_range in s, w_range in eV, cutoff in s, exc_E in J

    global BETA
    BETA = 1./ (KB * T)
    dt = t_range[1] - t_range[0]      # assuming equal spacing

    # Reading frequency and coupling matrix element data
    nExc = 20
    phonons = read_w(phonon_filename)
    coupling = read_Vklq(coupling_filename, len(phonons), nExc)
    phonons = phonons[6:]
    V = coupling[6:, 0]
    delta_range = (V) / phonons**2
    print("DONE reading")

    # Calculating the dephasing function
    (F_t_Re, F_t_Im) = calc_dephasing_F_E17(t_range, delta_range, phonons)
    print("DONE calc_dephasing_F_E17")
    
    '''
    # Calculating the delta-delta correlation function
    (DC_t_Re, DC_t_Im) = calc_delta_corr_E17(t_range, delta_range, phonons)
    DC_t_Re *= JTOEV**2
    DC_t_Im *= JTOEV**2
    print("DONE calc_delta_corr_E17")
    '''
    DC_t_Re = np.zeros(0)
    DC_t_Im = np.zeros(0)

    # erf correction to the dephasing functions
    F_t_Re_erf = multiply_erf(cutoff, slope, F_t_Re-np.mean(F_t_Re[-1000:-1]), dt)
    F_t_Im_erf = multiply_erf(cutoff, slope, F_t_Im-np.mean(F_t_Im[-1000:-1]), dt)
    print("DONE erf_multiplication")

    # FT the dephasing functions to obtain the spectrum from the dipole-dipole correlation function
    delta_E_E17 = exc_E + reorg_E(delta_range, phonons)    # in J
    (spectrum_w, spectrum_I) = spectrum(delta_E_E17, t_range, F_t_Re_erf, F_t_Im_erf, w_range)
    print("DONE spectrum")
    write_spectrum(I_file, spectrum_w, spectrum_I)

    (FWHM, l_index, r_index) = calc_FWHM(spectrum_w, spectrum_I)
    print("DONE calculating FWHM")
    write_FWHM(FWHM_file, spectrum_w, spectrum_I, l_index, r_index)

    return (spectrum_w, spectrum_I, FWHM, F_t_Re, F_t_Im, F_t_Re_erf, F_t_Im_erf, DC_t_Re, DC_t_Im)



# main
Temp_range = np.array([1, 2, 5, 10, 20, 30, 40])
Temp_range = np.append(Temp_range, np.arange(50, 450, 50))

linewidth_range = np.zeros(0)

fig, axs = plt.subplots(3, figsize=(6,18))
fig.suptitle("Spectral linewidth and temperature")      
axs[0].grid()
axs[1].grid()
axs[2].grid()

for T in Temp_range: 
    dt = 1e-15                                         # in s
    t_range = np.arange(0.0, 8e-12+dt, dt)            # in s
    # erf param
    cutoff = 5e-12                                     # in s
    slope = 0.01
    
    domega = 0.1e-3                                    # in eV
    w_range = np.arange(1., 3., domega)                # in eV
    
    exc_E = 0.080794824219 * AUTOJ               # in J
    
    phonon_filename = "w.dat"
    coupling_filename = "Vklq-diabatic.dat"
    I_file = "Spectrum_E17_5ps_T=%d.dat" % (T)
    FWHM_file = "FWHM_T=%d.dat" % (T)
    
    (spectrum_w, spectrum_I, linewidth, F_t_Re, F_t_Im, F_t_Re_erf, F_t_Im_erf, DC_t_Re, DC_t_Im) = calc_spectrum_E17(T, t_range, w_range, cutoff, slope, exc_E, phonon_filename, coupling_filename, I_file, FWHM_file)
    linewidth_range = np.append(linewidth_range, linewidth)
    axs[0].plot(spectrum_w*1000, spectrum_I, label="T=%d K" % T)   


axs[0].set(xlabel="Energey [meV]", title="Spectra", xlim=(1800, 2600))
axs[1].set(xlabel="Temperature [K]", ylabel="Linewidth [meV]", title="FWHM vs. T")
axs[1].plot(Temp_range, linewidth_range*1000, "*-")
axs[2].set(xlabel="Inverse Temperature [1/K]", ylabel="Linewidth [meV]", title=r"FWHM vs. $\beta$")
axs[2].plot(1/Temp_range, linewidth_range*1000, "*-")
axs[0].legend()
axs[1].legend()
axs[2].legend()
fig.tight_layout()
fig.savefig("nanocrystal_T_width.png")

'''
# plotting
fig, axs = plt.subplots(2,2, figsize=(8,8))
fig.suptitle("The dephasing function and spectra, E17 model, 3nmCdSe at T= %d K" % (T))

for i in range(2): 
    for j in range(2): 
        axs[i,j].grid()

axs[0,0].set(xlabel="Time (s)", ylabel=r"$\langle F(t) \rangle$ (unitless)", title="E17 dephasing function, cutoff at %.2f ps" % (cutoff*1e12), yscale="log")
axs[0,0].plot(t_range, F_t_Re, label=r"$\langle F(t) \rangle _{Re} $")
axs[0,0].plot(t_range, F_t_Im, ":", label=r"$\langle F(t) \rangle _{Im}$")
axs[0,0].plot(t_range, F_t_Re_erf, label=r"$\langle F(t) \rangle _{Re, shifted} \cdot erf$")
axs[0,0].plot(t_range, F_t_Im_erf, ":", label=r"$\langle F(t) \rangle _{Im, shifted} \cdot erf$")

# spectra
axs[1,0].set(xlabel="Energey [meV]", title="Spectra", xlim=(2000, 2400))
axs[1,0].plot(spectrum_w*1000, spectrum_I) #, label="erf cutoff at 5ps")

# Delta-Delta
axs[0,1].set(xlabel="Time (s)", ylabel=r"$\langle \Delta (t) \Delta (0) \rangle$, ($eV^2$)", title=r"E17 $\Delta - \Delta$ correlation function", yscale="log")
axs[0,1].plot(t_range, DC_t_Re, label=r"$C_{\Delta \Delta}(t)_{Re} $")
axs[0,1].plot(t_range, DC_t_Im, ":", label=r"$C_{\Delta \Delta}(t)_{Im} $")

for i in range(2): 
    for j in range(2): 
        axs[i,j].legend()
fig.tight_layout()
fig.savefig("nanocrystal_spectrum.png")
'''
