import numpy as np
import scipy as sp
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

PI = 3.14159265358979323846
a = 0.05                       # in units of THz

def J(x, A, B, C):
    return A * x**B * np.exp(-x/C)

def read_ph(filename):
    f = open(filename, "r")
    w_range = np.zeros(0)
    while True:
        line = f.readline()
        if not line:
            break
        if line[0]!= "#":
            w_range = np.append(w_range, float(line.split()[1]))
    f.close()
    # w_range *= 2*PI
    return w_range        # return \nu in units of THz

def read_nacoupling(filename):
    coupling = np.zeros((1332, 20))
    f = open(filename, "r")
    while True:
        line = f.readline()
        if not line:
            break
        if (line[0]!= "#") and (int(line.split()[1])==int(line.split()[2])):
            coupling[int(line.split()[0]), int(line.split()[1])] = float(line.split()[3])
    f.close()
    return (coupling)

def calc_y(w_range, V, gaussian_std, x_range):
    y_range = np.zeros(0)
    for x in x_range: 
        y = 0
        for i in range(len(w_range)):
            this_phonon = w_range[i]
            this_V = V[i] 
            y += np.sqrt(PI)/(2*gaussian_std) * this_V**2 / this_phonon * np.exp(-(x - this_phonon)**2/(a**2)) 
        y_range = np.append(y_range, y)
    return y_range

# main
phonon_w = read_ph("w.dat")[6:]
coupling = read_nacoupling("Vklq-diabatic.dat")
g = open("J_fitting.dat", "w")
h = open("3nmCdSe_spectralDensity.dat", "w")

n_count = 3

fig, axs = plt.subplots(n_count, figsize=(6,3*n_count))
fig.suptitle("3nm CdSe spectral density fitting")

# fit for Vnn
for n in range(20): 
    V = coupling[6:, n]
    x = np.arange(0, 10, 0.01)    # all in units of THz
    y = calc_y(phonon_w, V, a, x)

    if n==0:
        for l in range(len(x)):
            h.write("%f    %f \n" % (x[l], y[l]))

    popt, pcov = curve_fit(J, x, y) 
    print(popt)
    print(pcov)

    print("The form of the spectral density: %.3f * x**%.3f * np.exp(-x/%.3f)" % (popt[0], popt[1], popt[2]))
    g.write("# For the %d th exciton: J_{SK}(w) = %.3f * x^{%.3f} * exp(-x / %.3f)" % (n, popt[0], popt[1], popt[2]))
    g.write("\n")
    g.write("%d   %f   %f   %f \n" % (n, popt[0], popt[1], popt[2]))

    y_fitted = np.zeros(0)
    for this_x in x:
        y_fitted = np.append(y_fitted, J(this_x, popt[0], popt[1], popt[2]))

    if n<n_count:
        axs[n].grid()
        axs[n].set(xlabel=r"$\nu = \frac {\omega} {(2\pi)}$ , [THz]", ylabel=r"Spectral density, $J_{SK}(\omega)$, [arb units]")
        axs[n].set(title=r"For the %d th exciton: $J_{SK} (\omega) = %.3f \cdot x^{%.3f} \cdot exp(\frac {-x} {%.3f})$ " % (n, popt[0], popt[1], popt[2]))

        axs[n].plot(x, y, "+:", label="Spectral density")
        axs[n].plot(x, y_fitted, label="Fitted form")
        axs[n].legend()

g.close()
h.close()
fig.tight_layout()
fig.savefig("Spectral_density_fitting_3nmCdSe.png")