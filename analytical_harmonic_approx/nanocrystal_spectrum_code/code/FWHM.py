import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.fft import fft, ifft
from scipy.integrate import simpson

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
