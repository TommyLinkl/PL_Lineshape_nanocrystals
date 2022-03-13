import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.fft import fft, ifft
from scipy.integrate import simpson

def reorg_E(delta_range, omega_range):
    E = 0
    for i in range(len(omega_range)):
        E += 0.5 * omega_range[i]**2 * delta_range[i]**2
    print("reorg_E = %f J, %f eV \n" % (E, E*JTOEV))
    return E   # In units of J

