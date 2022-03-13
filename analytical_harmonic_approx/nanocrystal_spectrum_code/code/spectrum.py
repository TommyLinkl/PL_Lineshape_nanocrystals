import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.fft import fft, ifft
from scipy.integrate import simpson

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


