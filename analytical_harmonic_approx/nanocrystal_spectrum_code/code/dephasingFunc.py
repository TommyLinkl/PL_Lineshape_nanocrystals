import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.fft import fft, ifft
from scipy.integrate import simpson

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

