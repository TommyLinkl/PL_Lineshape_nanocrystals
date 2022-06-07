import numpy as np
import scipy as sp
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

PI = 3.14159265358979323846

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
    w_range *= 2*PI*1e12
    return w_range 

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

def calc_y(phonon_w_range, delta_range, gaussian_std, x_range):
    y_range = np.zeros(len(x_range))
    for i in range(len(phonon_w_range)):
        w = phonon_w_range[i]
        delta = delta_range[i]
        y_range += 1 / (gaussian_std*np.sqrt(2*PI)) * np.exp(-(x_range - w)**2/(2*gaussian_std**2)) * PI/2*w**3*delta**2 

    return y_range

# main
phonon_w = read_ph("w.dat")[6:]
V = read_nacoupling("Vklq-diabatic.dat")
coupling = V[6:, 0] / phonon_w**2
# g = open("J_fitting.dat", "w")

a = 1e11

fig, axs = plt.subplots(1, figsize=(6,3))
fig.suptitle("3nm CdSe spectral density fitting")

# fit for Vnn
x = np.arange(0, 50e12, 1e8)    # all in units of Hz for omega
y = calc_y(phonon_w, coupling, a, x)

popt, pcov = curve_fit(J, x, y) 
print(popt)
print(pcov)

print("The form of the spectral density: %.3f * x**%.3f * np.exp(-x/%.3f)" % (popt[0], popt[1], popt[2]))
# g.write("# For the %d th exciton: J_{SK}(w) = %.3f * x^{%.3f} * exp(-x / %.3f)" % (n, popt[0], popt[1], popt[2]))
# g.write("\n")
# g.write("%d   %f   %f   %f \n" % (n, popt[0], popt[1], popt[2]))

y_fitted = J(x, popt[0], popt[1], popt[2])

axs.grid()
axs.set(xlabel=r"$\omega$, [Hz]", ylabel=r"Spectral density, $J_{SK}(\omega)$, [arb units]")
axs.set(title=r"0th exciton: $J_{SK} (\omega) = %.3f \cdot x^{%.3f} \cdot exp(\frac {-x} {%.3f})$ " % (popt[0], popt[1], popt[2]))

axs.plot(x, y, ":", label="Spectral density")
axs.plot(x, y_fitted, label="Fitted form")
axs.legend()

# g.close()
fig.tight_layout()
fig.savefig("Spectral_density_fitting_3nmCdSe.png")
