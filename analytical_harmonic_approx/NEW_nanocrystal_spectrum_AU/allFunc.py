import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.fft import fft, ifft
from scipy.integrate import simpson

# Now all calculated in a.u.

PI = 3.14159265359
HBAR = 1.                                # 1 au = 1.054571817e-34 Js
JTOEV = 6.24150907446076e18 
AUTOEV = 27.211396641308 
AUTOJ = 4.359744e-18
AUTOS = 2.4188843265857e-17
KB = 1.380649e-23 / AUTOJ                # unit: hartree/K, a.u.

def coth(x):
    return (np.exp(x) + np.exp(-x))/(np.exp(x) - np.exp(-x))

def func_n(x, T): 
    BETA = 1 / (KB * T)
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
    w *= 1e12 * (2*PI)        # Now w in units of 1/s
    w *= AUTOS                # Now w in units of a.u.
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
    coupling *= AUTOS / np.sqrt(AUTOJ)        # Now coupling in units of a.u.
    return (coupling)

def reorg_E(delta_range, omega_range):
    E = 0
    for i in range(len(omega_range)):
        E += 0.5 * omega_range[i]**2 * delta_range[i]**2
    print("reorg_E = %f a.u., %f J, %f eV \n" % (E, E*AUTTOJ, E*AUTOEV))
    return E   # In units of a.u.

def calc_dephasing_F_E17(t_range, delta_range, omega_range, T): 
    if T==0:
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
    else:
        BETA = 1 / (KB * T)
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

def calc_dephasing_F_MC(t_range, delta_range, omega_range, T): 
    BETA = 1 / (KB * T)
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

def calc_dephasing_F_MC_classical(t_range, delta_range, omega_range, T): 
    BETA = 1 / (KB * T)
    F_t_Re = np.zeros(0)
    F_t_Im = np.zeros(0)
    for t in t_range: 
        Re = 0
        for i in range(len(omega_range)):
            Re += 1 / (HBAR**2 * BETA) * (delta_range[i])**2 * (np.cos(omega_range[i] * t) - 1)
        # print(Re)
        F_t_Re = np.append(F_t_Re, np.exp(Re))
        F_t_Im = np.append(F_t_Im, 0)
    return (F_t_Re, F_t_Im)

def calc_delta_corr_E17(t_range, delta_range, omega_range, T): 
    (C_t_Re, C_t_Im) = calc_delta_corr_MC(t_range, delta_range, omega_range, T)
    C_t_Re += reorg_E(delta_range, omega_range)**2
    return (C_t_Re, C_t_Im)       # in Hartree**2

def calc_delta_corr_MC(t_range, delta_range, omega_range, T): 
    BETA = 1 / (KB * T)
    C_t_Re = np.zeros(0)
    C_t_Im = np.zeros(0)
    for t in t_range: 
        Re = 0.
        Im = 0.
        for i in range(len(omega_range)): 
            w = omega_range[i]
            Re += HBAR / 2 * delta_range[i]**2 * w**3 * ((func_n(w, T) + 1)*np.cos(w*t) + func_n(w, T)*np.cos(w*t))
            Im += HBAR / 2 * delta_range[i]**2 * w**3 * ((func_n(w, T) + 1)*np.sin(-w*t) + func_n(w, T)*np.sin(w*t))
        C_t_Re = np.append(C_t_Re, Re)
        C_t_Im = np.append(C_t_Im, Im)
    return (C_t_Re, C_t_Im)       # in Hartree**2

'''
def write_dephasing_F(filename, t_range, F_Re_MC, F_Im_MC, F_Re_E17, F_Im_E17):
    f = open(filename, "w")
    f.write("# time(s)        F_Re_MC      F_Im_MC      F_Re_E17      F_Im_E17 \n")
    for i in range(len(t_range)):
        f.write("%.20f       %.10f       %.10f       %.10f       %.10f\n" % (t_range[i], F_Re_MC[i], F_Im_MC[i], F_Re_E17[i], F_Im_E17[i]))
    f.close()
    return

def write_delta_corr(filename, t_range, DC_Re_MC, DC_Im_MC, DC_Re_E17, DC_Im_E17):
    f = open(filename, "w")
    f.write("# time(s)        DeltaCorr_Re_MC(eV^2)      DeltaCorr_Im_MC(eV^2)      DeltaCorr_Re_E17(eV^2)      DeltaCorr_Im_E17(eV^2) \n")
    for i in range(len(t_range)):
        f.write("%.10f       %.10f       %.10f       %.10f       %.10f\n" % (t_range[i], DC_Re_MC[i], DC_Im_MC[i], DC_Re_E17[i], DC_Im_E17[i]))
    f.close()
    return

'''

def multiply_erf(cutoff, slope, corr_t, dt, t_range):
    corr_erf = np.zeros(0)
    for i in range(len(t_range)):
        corr_erf = np.append(corr_erf, corr_t[i]* (-0.5*math.erf(slope*(float(i)-cutoff/dt))+0.5))
    return corr_erf

'''
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
    
    I_norm = I / I.max()
    # return (w_range, I)    # return w_range in eV
    return (w_range, I_norm)    # return w_range in eV
    
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

def calc_freqShift(spec_w, spec_I, exciton_E):
    max_index = np.argmax(spec_I)
    max_freq = spec_w[max_index]
    return (max_freq - exciton_E)

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
    (F_t_Re, F_t_Im) = calc_dephasing_F_E17(t_range, delta_range, phonons, T)
    print("DONE calc_dephasing_F_E17")
    
    # Calculating the delta-delta correlation function
    (DC_t_Re, DC_t_Im) = calc_delta_corr_E17(t_range, delta_range, phonons)
    DC_t_Re *= JTOEV**2
    DC_t_Im *= JTOEV**2
    print("DONE calc_delta_corr_E17")
    DC_t_Re = np.zeros(0)
    DC_t_Im = np.zeros(0)

    # erf correction to the dephasing functions
    F_t_Re_erf = multiply_erf(cutoff, slope, F_t_Re-np.mean(F_t_Re[-1000:-1]), dt, t_range)
    F_t_Im_erf = multiply_erf(cutoff, slope, F_t_Im-np.mean(F_t_Im[-1000:-1]), dt, t_range)
    print("DONE erf_multiplication")

    # FT the dephasing functions to obtain the spectrum from the dipole-dipole correlation function
    delta_E_E17 = exc_E + reorg_E(delta_range, phonons)    # in J
    (spectrum_w, spectrum_I) = spectrum(delta_E_E17, t_range, F_t_Re_erf, F_t_Im_erf, w_range)
    print("DONE spectrum")
    write_spectrum(I_file, spectrum_w, spectrum_I)

    (FWHM, l_index, r_index) = calc_FWHM(spectrum_w, spectrum_I)
    print("DONE calculating FWHM")
    # write_FWHM(FWHM_file, spectrum_w, spectrum_I, l_index, r_index)

    freq_shift = calc_freqShift(spectrum_w, spectrum_I, exc_E)

    return (spectrum_w, spectrum_I, FWHM, freq_shift, F_t_Re, F_t_Im, F_t_Re_erf, F_t_Im_erf, DC_t_Re, DC_t_Im)


def calc_spectrum_classical_MC(T, t_range, w_range, cutoff, slope, exc_E, phonon_filename, coupling_filename, I_file, FWHM_file):
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
    (F_t_Re, F_t_Im) = calc_dephasing_F_MC_classical(t_range, delta_range, phonons, T)
    print("DONE calc_dephasing_F_MC_classical")
    
    DC_t_Re = np.zeros(0)
    DC_t_Im = np.zeros(0)

    # erf correction to the dephasing functions
    F_t_Re_erf = multiply_erf(cutoff, slope, F_t_Re-np.mean(F_t_Re[-1000:-1]), dt, t_range)
    F_t_Im_erf = multiply_erf(cutoff, slope, F_t_Im-np.mean(F_t_Im[-1000:-1]), dt, t_range)
    print("DONE erf_multiplication")

    # FT the dephasing functions to obtain the spectrum from the dipole-dipole correlation function
    (spectrum_w, spectrum_I) = spectrum(exc_E, t_range, F_t_Re_erf, F_t_Im_erf, w_range)
    print("DONE spectrum")
    write_spectrum(I_file, spectrum_w, spectrum_I)

    (FWHM, l_index, r_index) = calc_FWHM(spectrum_w, spectrum_I)
    # write_FWHM(FWHM_file, spectrum_w, spectrum_I, l_index, r_index)

    freq_shift = calc_freqShift(spectrum_w, spectrum_I, exc_E)

    return (spectrum_w, spectrum_I, FWHM, freq_shift, F_t_Re, F_t_Im, F_t_Re_erf, F_t_Im_erf, DC_t_Re, DC_t_Im)
'''
