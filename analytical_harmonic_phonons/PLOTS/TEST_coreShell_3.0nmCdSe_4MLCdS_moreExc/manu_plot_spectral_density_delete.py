import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib import rc
from allFunc import * 
matplotlib.rcParams['font.family'] = "serif"
matplotlib.rcParams['font.sans-serif'] = "Times New Roman"
rc('axes', linewidth=0.4)
matplotlib.rcParams['xtick.major.width'] = 0.4
matplotlib.rcParams['ytick.major.width'] = 0.4

gaussian_std = 0.03  # THz

def delta_to_Gaussian(center, gaussian_std, x):
    y = 1/(gaussian_std * np.sqrt(2*np.pi)) * np.exp(-(x - center)**2 / (2*gaussian_std**2))
    return y

def PDOS(w_range, gaussian_std):
    x_range = np.arange(0-gaussian_std, max(w_range)+3*gaussian_std, gaussian_std/100)
    y_range = np.zeros(len(x_range))
    for i in range(len(w_range)):
        w = w_range[i]
        y_range += delta_to_Gaussian(w, gaussian_std, x_range)
    return (x_range, y_range)

def SK_spectral_density(w_range, V_range, gaussian_std):
    x_range = np.arange(0-gaussian_std, max(w_range)+3*gaussian_std, gaussian_std/100)
    y_range = np.zeros(len(x_range))
    for i in range(len(w_range)):
        w = w_range[i]
        V = V_range[i]
        y_range += delta_to_Gaussian(w, gaussian_std, x_range) * np.pi/2*V**2 / w
    return (x_range, y_range)

def Dipti_spectral_density(w_range, V_range, gaussian_std):
    x_range = np.arange(0-3*gaussian_std, max(w_range)+3*gaussian_std, gaussian_std/100)
    y_range = np.zeros(len(x_range))
    for i in range(len(w_range)):
        w = w_range[i]
        V = V_range[i]
        y_range += delta_to_Gaussian(w, gaussian_std, x_range) * 1/2 * V**2 / w**2
    return (x_range, y_range)


# Data
phonons = np.array([float(line.strip().split()[1]) for line in open("../../CALCS/coreShell_3.0nmCdSe_4MLCdS/w.dat", 'r') if line[0]!="#"])
OLD_V = np.array([float(line.strip().split()[3]) for line in open("../../CALCS/coreShell_3.0nmCdSe_4MLCdS/OLD_Vklq-diabatic.dat", 'r') if line[0]!="#"])
NEW_V = np.array([float(line.strip().split()[3]) for line in open("../../CALCS/coreShell_3.0nmCdSe_4MLCdS/NEW_Vklq-diabatic.dat", 'r') if line[0]!="#"])
NEW_delete_V = np.array([float(line.strip().split()[3]) for line in open("../../CALCS/coreShell_3.0nmCdSe_4MLCdS/NEW_Vklq-diabatic_delete.dat", 'r') if line[0]!="#"])
phonons = phonons[6:]
OLD_V = OLD_V[6:]
NEW_V = NEW_V[6:]

# (x_OLD, y_OLD) = Dipti_spectral_density(phonons, OLD_V, gaussian_std)
(x_NEW, y_NEW) = Dipti_spectral_density(phonons, NEW_V, gaussian_std)
(x_NEW_delete, y_NEW_delete) = Dipti_spectral_density(phonons, NEW_delete_V, gaussian_std)
(x_inset, y_inset) = PDOS(phonons, gaussian_std)


fig, axs = plt.subplots(1,1, figsize=(6, 4))
fig.suptitle("3nm CdSe/4ML CdS core-shell quantum dot")

# axs.plot(x_OLD, y_OLD, "-r", linewidth=1, label="Old code")
# axs.plot(x_NEW, y_NEW, "-k", linewidth=1, label="New code")
axs.plot(x_NEW_delete, y_NEW_delete, "-k", linewidth=1, label="New code, delete 2 modes")
axs.get_yaxis().set_ticks([0])
axs.set(xlim=(0,11.4), xlabel="Frequency [THz]", ylabel="Spectral density", title=r"$J_{Dipti}\left( \omega \right) = \sum _\alpha \left( \frac {V^\alpha}{\omega_\alpha} \right)^2 \delta \left( \omega - \omega_\alpha \right)$")
axs.set_ylim(bottom=0, top=None)
axs.legend(loc="upper left")

axins = inset_axes(axs, width="100%", height="100%", bbox_to_anchor=(.2, .5, .4, .3), bbox_transform=axs.transAxes)
axins.plot(x_inset, y_inset, "-k", linewidth=1)
axins.get_yaxis().set_ticks([])
axins.set(xlim=(0,11.4), ylabel="PDOS")
axins.set_ylim(bottom=0, top=None)

'''
(x_OLD_SK, y_OLD_SK) = SK_spectral_density(phonons, OLD_V, gaussian_std)
(x_NEW_SK, y_NEW_SK) = SK_spectral_density(phonons, NEW_V, gaussian_std)

axs[1].plot(x_OLD_SK, y_OLD_SK, "-r", linewidth=1, label="Old code")
axs[1].plot(x_NEW_SK, y_NEW_SK, "-k", linewidth=1, label="New code")
axs[1].get_yaxis().set_ticks([0])
axs[1].set(xlim=(0,11.4), xlabel="Frequency [THz]", ylabel="Spectral density", title=r"$J_{Skinner}\left( \omega \right) = \sum _\alpha \frac {\left( V^\alpha \right)^2}{\omega_\alpha}  \delta \left( \omega - \omega_\alpha \right)$")
axs[1].set_ylim(bottom=0, top=None)
axs[1].legend(loc="upper left")
'''

fig.tight_layout()
fig.savefig("spectral_density_comparison_delete.pdf")
