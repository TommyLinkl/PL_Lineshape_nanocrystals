import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.fft import fft, ifft
from scipy.integrate import simpson

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
