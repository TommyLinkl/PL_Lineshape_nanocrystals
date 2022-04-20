import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.integrate import simpson

T = 1                            # unit: K 
KB = 1.380649e-23                # unit: J/K, SI units
BETA = 1./ (KB * T)   
PI = 3.14159265359
HBAR = 1.054571817e-34           # Unit: Js, SI units
JTOEV = 6.24150907446076e18 
AUTOEV = 27.211396641308 
AUTOJ = 4.359744e-18

dt = 1e-15
t_range = np.arange(0.0, 1e-11+dt, dt)
domega = 1e12                 # in Hz

def read_spectrum(filename):
    w = np.zeros(0)
    I_w = np.zeros(0)
    f = open(filename, "r")
    while True:
        line = f.readline()
        if not line:
            break
        if line[0]!= "#":
            w = np.append(w, float(line.split()[0]))
            I_w = np.append(I_w, float(line.split()[1]))
    f.close()
    return (w, I_w)

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

# main

# Reading frequency and coupling matrix element data
(w, I_w) = read_spectrum("3nmCdSe_spectrum_E17.dat")
print("DONE reading")

(width, l_index, r_index) = calc_FWHM(w, I_w)
print("DONE calculating FWHM")

write_FWHM("FWHM_E17.dat", w, I_w, l_index, r_index)
