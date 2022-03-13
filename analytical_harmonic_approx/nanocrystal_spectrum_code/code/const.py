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
