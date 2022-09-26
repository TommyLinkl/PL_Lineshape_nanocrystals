import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.fft import fft, ifft
from scipy.integrate import simpson
from scipy.special import erf

KB = 1.380649e-23                # unit: J/K, SI units
PI = 3.14159265359
HBAR = 1.054571817e-34           # Unit: Js, SI units
JTOEV = 6.24150907446076e18 
AUTOEV = 27.211396641308 
AUTOJ = 4.359744e-18

def coth(x):
    return 1/np.tanh(x)

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
    w *= 1e12 * (2*PI)     # in 1/s
    # print(w)
    return w

def read_Vklq(filename, nPhonon, nExc, shrink):
    coupling = np.zeros((nPhonon, nExc))
    f = open(filename, "r")
    while True:
        line = f.readline()
        if not line:
            break
        if (line[0]!= "#") and (int(line.split()[1])==int(line.split()[2])):
            coupling[int(line.split()[0]), int(line.split()[1])] = float(line.split()[3])
    f.close()
    # print(coupling*shrink)
    return (coupling*shrink)      # in sqrt(J) / s

def rescale_Vklq(w_range, coupling, ac_low, ac_upper, alpha, op_low, op_upper, beta, \
    int_low, int_upper, gamma, new_V_filename):
#################################################
## BUG!!!
#################################################
    # w_range doesn't include the 6 zero-frequency phonons
    # coupling is a 1D array of length nPhonon. Only V_00 is kept in this input. Units: sqrt(J) / s
    rescaled_coupling = np.array(coupling)
    f = open(new_V_filename, "w")
    f.write("0 0 0 0.0\n1 0 0 0.0\n2 0 0 0.0\n3 0 0 0.0\n4 0 0 0.0\n5 0 0 0.0\n")
    for i in range(len(w_range)): 
        if (w_range[i]>=ac_low*1e12*2*PI) and (w_range[i]<=ac_upper*1e12*2*PI): 
            rescaled_coupling[i] *= alpha
        if (w_range[i]>=int_low*1e12*2*PI) and (w_range[i]<=int_upper*1e12*2*PI): 
            rescaled_coupling[i] *= gamma
        if (w_range[i]>=op_low*1e12*2*PI) and (w_range[i]<=op_upper*1e12*2*PI): 
            rescaled_coupling[i] *= beta
        f.write("%d 0 0 %.16f\n" % (i+6, rescaled_coupling[i]))
    f.close()
    return rescaled_coupling      # in sqrt(J) / s

def NEW_rescale_Vklq(w_range, coupling, a1_lower, a1_upper, a1, a2_lower, a2_upper, a2, a3_lower, a3_upper, a3, a4_lower, a4_upper, a4, a5_lower, a5_upper, a5, b_lower, b_upper, b, new_V_filename):
    # w_range doesn't include the 6 zero-frequency phonons
    # coupling is a 1D array of length nPhonon. Only V_00 is kept in this input. Units: sqrt(J) / s
    rescaled_coupling = np.reshape(coupling, -1)
    f = open(new_V_filename, "w")
    f.write("0 0 0 0.0\n1 0 0 0.0\n2 0 0 0.0\n3 0 0 0.0\n4 0 0 0.0\n5 0 0 0.0\n")
    for i in range(6, len(w_range)): 
        if (w_range[i]>=a1_lower*1e12*2*PI) and (w_range[i]<a1_upper*1e12*2*PI): 
            rescaled_coupling[i] *= a1
        if (w_range[i]>=a2_lower*1e12*2*PI) and (w_range[i]<a2_upper*1e12*2*PI): 
            rescaled_coupling[i] *= a2
        if (w_range[i]>=a3_lower*1e12*2*PI) and (w_range[i]<a3_upper*1e12*2*PI): 
            rescaled_coupling[i] *= a3
        if (w_range[i]>=a4_lower*1e12*2*PI) and (w_range[i]<a4_upper*1e12*2*PI): 
            rescaled_coupling[i] *= a4
        if (w_range[i]>=a5_lower*1e12*2*PI) and (w_range[i]<a5_upper*1e12*2*PI): 
            rescaled_coupling[i] *= a5
        if (w_range[i]>=b_lower*1e12*2*PI) and (w_range[i]<b_upper*1e12*2*PI): 
            rescaled_coupling[i] *= b
        f.write("%d 0 0 %.16f\n" % (i, rescaled_coupling[i]))
    f.close()
    return rescaled_coupling      # in sqrt(J) / s

def delete_modes(w_filename, V_filename, phonon_delete_range, new_w_filename, new_V_filename): 
    all_w = read_w(w_filename)
    all_coupling = read_Vklq(V_filename, len(all_w), 1, 1.0)
    
    new_w = np.delete(all_w, phonon_delete_range)
    new_coupling = np.delete(all_coupling, phonon_delete_range)
    print(np.shape(new_w))
    print(np.shape(new_coupling))

    f = open(new_w_filename, "w")
    for i in range(len(new_w)):
        f.write("%d %.16f\n" % (i,new_w[i]/(2*PI*1e12)))
    f.close()

    f = open(new_V_filename, "w")
    for i in range(len(new_coupling)): 
        f.write("%d 0 0 %.16f\n" % (i, new_coupling[i]))
    f.close()
    return




def reorg_E(delta_range, omega_range):
    E = np.sum(0.5 * omega_range**2 * delta_range**2)
    print("reorg_E = %f J, %f eV, %6f meV" % (E, E*JTOEV, E*JTOEV*1000))
    return E   # In units of J

def calc_dephasing_F_E17(t_range, delta_range, omega_range, T): 
    if T==0:
        Re = 1 / (2*HBAR) * np.sum((delta_range**2 * omega_range).reshape(-1, 1) * (np.cos(np.outer(omega_range, t_range))-1), axis=0)
        Im = -1 / (2*HBAR) * np.sum((delta_range**2 * omega_range).reshape(-1, 1) * (-np.sin(np.outer(omega_range, t_range))), axis=0)
        F_t_Re = np.exp(Re) * np.cos(Im)
        F_t_Im = np.exp(Re) * np.sin(Im)
    else:
        BETA = 1 / (KB * T)
        Re = 1 / (2*HBAR) * np.sum((delta_range**2 * omega_range * coth(BETA*HBAR*omega_range/2)).reshape(-1, 1) * (np.cos(np.outer(omega_range, t_range))-1), axis=0)
        Im = -1 / (2*HBAR) * np.sum((delta_range**2 * omega_range).reshape(-1, 1) * (-np.sin(np.outer(omega_range, t_range))), axis=0)
        F_t_Re = np.exp(Re) * np.cos(Im)
        F_t_Im = np.exp(Re) * np.sin(Im)
    return (F_t_Re, F_t_Im)

def calc_dephasing_F_MC(t_range, delta_range, omega_range, T): 
    if T==0:
        Re = 1 / (2*HBAR) * np.sum((delta_range**2 * omega_range).reshape(-1, 1) * (np.cos(np.outer(omega_range, t_range))-1), axis=0)
        Im = -1 / (2*HBAR) * np.sum((delta_range**2 * omega_range).reshape(-1, 1) * (-np.sin(np.outer(omega_range, t_range)) + np.outer(omega_range, t_range)), axis=0)
        F_t_Re = np.exp(Re) * np.cos(Im)
        F_t_Im = np.exp(Re) * np.sin(Im)
    else:
        BETA = 1 / (KB * T)
        Re = 1 / (2*HBAR) * np.sum((delta_range**2 * omega_range * coth(BETA*HBAR*omega_range/2)).reshape(-1, 1) * (np.cos(np.outer(omega_range, t_range))-1), axis=0)
        Im = -1 / (2*HBAR) * np.sum((delta_range**2 * omega_range).reshape(-1, 1) * (-np.sin(np.outer(omega_range, t_range)) + np.outer(omega_range, t_range)), axis=0)
        F_t_Re = np.exp(Re) * np.cos(Im)
        F_t_Im = np.exp(Re) * np.sin(Im)
    return (F_t_Re, F_t_Im)

def calc_dephasing_F_E17_classical(t_range, delta_range, omega_range, T): 
    Re = 1 / (HBAR**2) * KB * T * np.sum((delta_range**2).reshape(-1, 1) * (np.cos(np.outer(omega_range, t_range))-1), axis=0)
    Im = -1 / (HBAR) * reorg_E(delta_range, omega_range) * t_range
    F_t_Re = np.exp(Re) * np.cos(Im)
    F_t_Im = np.exp(Re) * np.sin(Im)
    return (F_t_Re, F_t_Im)

def calc_dephasing_F_MC_classical(t_range, delta_range, omega_range, T): 
    Re = 1 / (HBAR**2) * KB * T * np.sum((delta_range**2).reshape(-1, 1) * (np.cos(np.outer(omega_range, t_range))-1), axis=0)
    F_t_Re = np.exp(Re)
    F_t_Im = np.zeros(len(t_range))
    return (F_t_Re, F_t_Im)

def calc_delta_corr_E17(t_range, delta_range, omega_range, T): 
    (C_t_Re, C_t_Im) = calc_delta_corr_MC(t_range, delta_range, omega_range, T)
    C_t_Re += reorg_E(delta_range, omega_range)**2
    return (C_t_Re, C_t_Im)

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
    return (C_t_Re, C_t_Im)

def write_dephasing_F(filename, t_range, F_Re_MC, F_Im_MC, F_Re_E17, F_Im_E17, F_CL):
    f = open(filename, "w")
    f.write("# time        F_Re_MC      F_Im_MC      F_Re_E17      F_Im_E17     F_CL\n")
    for i in range(len(t_range)):
        f.write("%.20f       %.20f       %.20f       %.20f       %.20f     %.20f\n" % (t_range[i], F_Re_MC[i], F_Im_MC[i], F_Re_E17[i], F_Im_E17[i], F_CL[i]))
    f.close()
    return

def write_delta_corr(filename, t_range, DC_Re_MC, DC_Im_MC, DC_Re_E17, DC_Im_E17):
    f = open(filename, "w")
    f.write("# time        DeltaCorr_Re_MC      DeltaCorr_Im_MC      DeltaCorr_Re_E17      DeltaCorr_Im_E17 \n")
    for i in range(len(t_range)):
        f.write("%.10f       %.10f       %.10f       %.10f       %.10f\n" % (t_range[i], DC_Re_MC[i], DC_Im_MC[i], DC_Re_E17[i], DC_Im_E17[i]))
    f.close()
    return

def multiply_erf(cutoff, slope, corr_t, dt, t_range):
    corr_erf = corr_t * (-0.5*erf(slope*(t_range/dt-cutoff/dt))+0.5)
    return corr_erf

def spectrum(E_shift, t_range, F_Re, F_Im, E_range):
    # E_shift  in J, E_range in eV
    I = np.zeros(0)
    print("E_shift: %f J, %f eV, %f Hartree. This is equal to omega= %f Hz " % (E_shift, E_shift*JTOEV, E_shift/AUTOJ, E_shift / HBAR))

    F_t = F_Re + 1j * F_Im
    I = np.fft.hfft(F_t)
    I_norm =  I / I.max()

    N = 2 * (len(t_range) - 1) 
    dt = t_range[1] - t_range[0]
    E_range = np.fft.fftfreq(N, d=dt)
    E_range *= HBAR * JTOEV * (2 * PI)
    E_range += E_shift * JTOEV
    
    E_range = np.append(E_range[int(len(E_range)/2):], E_range[0:int(len(E_range)/2)])
    I_norm = np.append(I_norm[int(len(I_norm)/2):], I_norm[0:int(len(I_norm)/2)])

    return (E_range, I_norm)    # return E_range in eV
    
def write_spectrum(filename, w_range, spectrum):
    f = open(filename, "w")
    f.write("# Frequency w         I(w)\n")
    for i in range(len(w_range)):
        f.write("%.20f       %.20f\n" % (w_range[i], spectrum[i]))
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

def calc_spectrum_QM_E17(T, t_range, w_range, cutoff, slope, exc_E, phonon_filename, \
    coupling_filename, I_file, FWHM_file, shrink):
    # give T in K, t_range in s, w_range in eV, cutoff in s, exc_E in J

    dt = t_range[1] - t_range[0]      # assuming equal spacing

    # Reading frequency and coupling matrix element data
    nExc = 1
    phonons = read_w(phonon_filename)
    coupling = read_Vklq(coupling_filename, len(phonons), nExc, shrink)
    phonons = phonons[6:]
    V = coupling[6:, 0]
    delta_range = (V) / phonons**2
    print("DONE reading")

    # Calculating the dephasing function
    (F_t_Re, F_t_Im) = calc_dephasing_F_E17(t_range, delta_range, phonons, T)
    ########################################################################################
    # F_t_Re = np.array([1.0]*len(t_range))
    # F_t_Im = np.zeros(len(t_range))
    ########################################################################################
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
    # F_t_Re_erf = multiply_erf(cutoff, slope, F_t_Re-np.mean(F_t_Re[-1000:-1]), dt, t_range)
    # F_t_Im_erf = multiply_erf(cutoff, slope, F_t_Im-np.mean(F_t_Im[-1000:-1]), dt, t_range)
    F_t_Re_erf = multiply_erf(cutoff, slope, F_t_Re, dt, t_range)
    F_t_Im_erf = multiply_erf(cutoff, slope, F_t_Im, dt, t_range)
    print("DONE erf_multiplication")

    # FT the dephasing functions to obtain the spectrum from the dipole-dipole correlation function
    delta_E_E17 = exc_E   #  - reorg_E(delta_range, phonons)    # in J
    reorg_E(delta_range, phonons) 
    (spectrum_w, spectrum_I) = spectrum(delta_E_E17, t_range, F_t_Re_erf, F_t_Im_erf, w_range)
    print("DONE spectrum")
    write_spectrum(I_file, spectrum_w, spectrum_I)

    (FWHM, l_index, r_index) = calc_FWHM(spectrum_w, spectrum_I)
    print("DONE calculating FWHM\n")
    # write_FWHM(FWHM_file, spectrum_w, spectrum_I, l_index, r_index)

    freq_shift = calc_freqShift(spectrum_w, spectrum_I, exc_E)

    return (spectrum_w, spectrum_I, FWHM, freq_shift, F_t_Re, F_t_Im, F_t_Re_erf, F_t_Im_erf, DC_t_Re, DC_t_Im)






def calc_spectrum_QM_MC(T, t_range, w_range, cutoff, slope, exc_E, phonon_filename, coupling_filename, I_file, FWHM_file, shrink):
    # give T in K, t_range in s, w_range in eV, cutoff in s, exc_E in J

    dt = t_range[1] - t_range[0]      # assuming equal spacing

    # Reading frequency and coupling matrix element data
    nExc = 1
    phonons = read_w(phonon_filename)
    coupling = read_Vklq(coupling_filename, len(phonons), nExc, shrink)
    phonons = phonons[6:]
    V = coupling[6:, 0]
    delta_range = (V) / phonons**2
    print("DONE reading")

    # Calculating the dephasing function
    (F_t_Re, F_t_Im) = calc_dephasing_F_MC(t_range, delta_range, phonons, T)
    print("DONE calc_dephasing_F_MC")
    
    '''
    # Calculating the delta-delta correlation function
    (DC_t_Re, DC_t_Im) = calc_delta_corr_MC(t_range, delta_range, phonons)
    DC_t_Re *= JTOEV**2
    DC_t_Im *= JTOEV**2
    print("DONE calc_delta_corr_MC")
    '''
    DC_t_Re = np.zeros(0)
    DC_t_Im = np.zeros(0)

    # erf correction to the dephasing functions
    # F_t_Re_erf = multiply_erf(cutoff, slope, F_t_Re-np.mean(F_t_Re[-1000:-1]), dt, t_range)
    # F_t_Im_erf = multiply_erf(cutoff, slope, F_t_Im-np.mean(F_t_Im[-1000:-1]), dt, t_range)
    F_t_Re_erf = multiply_erf(cutoff, slope, F_t_Re, dt, t_range)
    F_t_Im_erf = multiply_erf(cutoff, slope, F_t_Im, dt, t_range)
    print("DONE erf_multiplication")

    # FT the dephasing functions to obtain the spectrum from the dipole-dipole correlation function
    delta_E_MC = exc_E    # in J
    (spectrum_w, spectrum_I) = spectrum(delta_E_MC, t_range, F_t_Re_erf, F_t_Im_erf, w_range)
    print("DONE spectrum")
    write_spectrum(I_file, spectrum_w, spectrum_I)

    (FWHM, l_index, r_index) = calc_FWHM(spectrum_w, spectrum_I)
    print("DONE calculating FWHM")
    # write_FWHM(FWHM_file, spectrum_w, spectrum_I, l_index, r_index)

    freq_shift = calc_freqShift(spectrum_w, spectrum_I, exc_E)

    return (spectrum_w, spectrum_I, FWHM, freq_shift, F_t_Re, F_t_Im, F_t_Re_erf, F_t_Im_erf, DC_t_Re, DC_t_Im)
 

def calc_spectrum_classical_MC(T, t_range, w_range, cutoff, slope, exc_E, phonon_filename, coupling_filename, I_file, FWHM_file, shrink):
    # give T in K, t_range in s, w_range in eV, cutoff in s, exc_E in J

    dt = t_range[1] - t_range[0]      # assuming equal spacing

    # Reading frequency and coupling matrix element data
    nExc = 1
    phonons = read_w(phonon_filename)
    coupling = read_Vklq(coupling_filename, len(phonons), nExc, shrink)
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
    F_t_Re_erf = multiply_erf(cutoff, slope, F_t_Re-np.mean(F_t_Re[-2000:-1]), dt, t_range)
    F_t_Im_erf = multiply_erf(cutoff, slope, F_t_Im-np.mean(F_t_Im[-2000:-1]), dt, t_range)
    print("DONE erf_multiplication")

    # FT the dephasing functions to obtain the spectrum from the dipole-dipole correlation function
    (spectrum_w, spectrum_I) = spectrum(exc_E, t_range, F_t_Re_erf, F_t_Im_erf, w_range)
    print("DONE spectrum")
    write_spectrum(I_file, spectrum_w, spectrum_I)

    (FWHM, l_index, r_index) = calc_FWHM(spectrum_w, spectrum_I)
    # write_FWHM(FWHM_file, spectrum_w, spectrum_I, l_index, r_index)

    freq_shift = calc_freqShift(spectrum_w, spectrum_I, exc_E*JTOEV)

    return (spectrum_w, spectrum_I, FWHM, freq_shift, F_t_Re, F_t_Im, F_t_Re_erf, F_t_Im_erf, DC_t_Re, DC_t_Im)


def calc_spectrum_QM_E17_noCoupling(T, t_range, w_range, cutoff, slope, exc_E, phonon_filename, coupling_filename, I_file, FWHM_file, shrink):
    # give T in K, t_range in s, w_range in eV, cutoff in s, exc_E in J

    dt = t_range[1] - t_range[0]      # assuming equal spacing

    # Reading frequency and coupling matrix element data
    nExc = 1
    phonons = read_w(phonon_filename)
    coupling = read_Vklq(coupling_filename, len(phonons), nExc, shrink)
    phonons = phonons[6:]
    V = coupling[6:, 0]
    delta_range = (V) / phonons**2
    delta_range = np.zeros(len(delta_range))
    print("DONE reading")

    # Calculating the dephasing function
    (F_t_Re, F_t_Im) = calc_dephasing_F_E17(t_range, delta_range, phonons, T)
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


def calc_spectrum_QM_E17_fitted(T, t_range, w_range, cutoff, slope, exc_E, phonon_filename, coupling_filename, I_file, FWHM_file, shrink):
    # give T in K, t_range in s, w_range in eV, cutoff in s, exc_E in J

    global BETA
    BETA = 1./ (KB * T)
    dt = t_range[1] - t_range[0]      # assuming equal spacing

    # Reading frequency and coupling matrix element data
    nExc = 1
    phonons = read_w(phonon_filename)
    coupling = read_Vklq(coupling_filename, len(phonons), nExc, shrink)
    phonons = phonons[6:]
    V = coupling[6:, 0]
    delta_range = (V) / phonons**2
    print(phonons)
    print(delta_range)
    plt.plot(phonons, delta_range, label="original")

    delta_range = np.zeros(0)
    A = 4427036.251
    B = 3.728
    C = -0.401
    for i in range(len(phonons)):
        delta_range = np.append(delta_range, np.sqrt(A * (phonons[i]/2/PI*1e-12)**(B) * np.exp((phonons[i]/2/PI*1e-12)/C) * 2 / PI / (phonons[i])**3))
    print(delta_range)
    print("DONE reading")
    # plt.plot(phonons, delta_range, label="fitted")
    plt.legend()
    plt.grid()
    plt.savefig("testing.png")


    # Calculating the dephasing function
    (F_t_Re, F_t_Im) = calc_dephasing_F_E17(t_range, delta_range, phonons, T)
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

    
