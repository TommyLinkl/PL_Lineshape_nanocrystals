import numpy as np
import scipy as sp
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import argparse

PI = 3.14159265358979323846
HBAR = 1.054571817e-34         # Js
JTOEV = 6.242e18

nExc = 20

parser = argparse.ArgumentParser()
parser.add_argument('w_file', type=str)
parser.add_argument('Vklq_file', type=str)
parser.add_argument('output_file', type=str)
args = parser.parse_args()

def J(x, A, B, C):
    return A * x**B * np.exp(-x/C)

def delta_to_Gaussian(center, gaussian_std, x):
    y = 1/(gaussian_std * np.sqrt(2*PI)) * np.exp(-(x - center)**2 / (2*gaussian_std**2))
    return y

def read_ph(filename):
    w = np.zeros(0)
    f = open(filename, "r")
    while True:
        line = f.readline()
        if not line:
            break
        if line[0]!= "#":
            w = np.append(w, float(line.split()[1]))
    f.close()
    w *= 1e12 * 2 * PI
    return w      # in Hz

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
    return (coupling)      # in sqrt(J) / s

def PDOS(w_range, gaussian_std):
    # w_range given in Hz
    x_range = np.arange(0, max(w_range)+1e12, 1e8)
    y_range = np.zeros(len(x_range))
    for i in range(len(w_range)):
        w = w_range[i]
        y_range += delta_to_Gaussian(w, gaussian_std, x_range)
    return (x_range, y_range)

def SK_spectral_density(w_range, delta_range, gaussian_std):
    # w_range given in Hz, 
    x_range = np.arange(0, max(w_range)+1e12, 1e8)
    y_range = np.zeros(len(x_range))
    for i in range(len(w_range)):
        w = w_range[i]
        delta = delta_range[i]
        y_range += delta_to_Gaussian(w, gaussian_std, x_range) * PI/2*w**3*delta**2
    return (x_range, y_range)

def Dipti_spectral_density(w_range, reorg_e_range, gaussian_std):
    # w_range given in Hz, 
    x_range = np.arange(0, max(w_range)+1e12, 1e8)
    y_range = np.zeros(len(x_range))
    for i in range(len(w_range)):
        w = w_range[i]
        reorg_e = reorg_e_range[i]
        y_range += delta_to_Gaussian(w, gaussian_std, x_range) * reorg_e
    return (x_range, y_range)

# main
a = 2e11            # in units of Hz
phonons = read_ph(args.w_file)
V = read_Vklq(args.Vklq_file, len(phonons), nExc)
print("DONE reading. ")

fig, axs = plt.subplots(6, 1 , figsize=(6,18))
fig.suptitle(args.output_file)
(x, y) = PDOS(phonons, a)
axs[0].plot(x / (2*PI*1e12), y)
axs[0].grid()
axs[0].set(title="PDOS", xlabel="Phonon frequency [THz]", ylabel="PDOS")
axs[0].axes.yaxis.set_ticks([])

axs[1].plot(phonons/(2*PI*1e12), V[:,0], "o", markersize=1.2)
axs[1].grid()
axs[1].set(title="Vklq naCoupling", xlabel="Phonon frequency [THz]", ylabel=r"$V_{0,0} ^ {\alpha}$ $[\frac {\sqrt{J}} {s}]$")

phonons = phonons[6:]
delta_range = V[6:, 0] / phonons**2

(x,y) = SK_spectral_density(phonons, delta_range, a)
axs[2].plot(x/(2*PI*1e12), y)
axs[2].grid()
axs[2].set(title="Skinner's Spectral density", xlabel="Phonon frequency [THz]", ylabel="[J]")

reorg_e_range = 0.5 * (delta_range)**2 * phonons**2 * JTOEV * 1000
print("Reorganization energy is %f meV\n" % np.sum(reorg_e_range))
axs[3].plot(phonons/(2*PI*1e12), reorg_e_range, "o", markersize=1.2)
axs[3].grid()
axs[3].set(title="Reorganization energy", xlabel="Phonon frequency [THz]", ylabel=r"$\lambda$ [meV]") # , yscale="log")

huang_range = reorg_e_range / HBAR / phonons 
axs[4].plot(phonons/(2*PI*1e12), huang_range, "o", markersize=1.2)
axs[4].grid()
axs[4].set(title="Huang-Rhys parameter", xlabel="Phonon frequency [THz]", ylabel=r"$S_{\alpha}$ [dimensionless]") #, yscale="log")

(x,y) = Dipti_spectral_density(phonons, reorg_e_range, a)
axs[5].plot(x/(2*PI*1e12), y, "o-", markersize=1.2)
axs[5].grid()
axs[5].set(title="Dipti's Spectral density", xlabel="Phonon frequency [THz]", ylabel="[Js]")

fig.tight_layout()
fig.savefig(args.output_file)
