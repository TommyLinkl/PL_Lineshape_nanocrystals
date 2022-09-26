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
# coupling = np.zeros((nPhonon, nExc))
all_V = read_Vklq("../../../../NC_exc_elph/StrainCube_filter_bse_naCoupling/coreShell_3.0nmCdSe_4MLCdS/Vklq_NEW_moreExc/Vklq-diabatic.dat", 10728, 10, 1.0)
V_0 = all_V[:,0]
V_1 = all_V[:,1]
V_2 = all_V[:,2]
V_3 = all_V[:,3]
V_4 = all_V[:,4]
print(V_0)

phonons = phonons[6:]
V_0 = V_0[6:]
V_1 = V_1[6:]
V_2 = V_2[6:]
V_3 = V_3[6:]
V_4 = V_4[6:]

(x_NEW_0, y_NEW_0) = Dipti_spectral_density(phonons, V_0, gaussian_std)
(x_NEW_1, y_NEW_1) = Dipti_spectral_density(phonons, V_1, gaussian_std)
(x_NEW_2, y_NEW_2) = Dipti_spectral_density(phonons, V_2, gaussian_std)
(x_NEW_3, y_NEW_3) = Dipti_spectral_density(phonons, V_3, gaussian_std)
(x_NEW_4, y_NEW_4) = Dipti_spectral_density(phonons, V_4, gaussian_std)
(x_inset, y_inset) = PDOS(phonons, gaussian_std)


fig, axs = plt.subplots(1,1, figsize=(6, 4))
fig.suptitle("3nm CdSe/4ML CdS core-shell quantum dot")

axs.plot(x_NEW_0, y_NEW_0, "-", linewidth=1, label="New code, Exc#0")
axs.plot(x_NEW_1, y_NEW_1, "-", linewidth=1, label="New code, Exc#1")
axs.plot(x_NEW_2, y_NEW_2, "-", linewidth=1, label="New code, Exc#2")
axs.plot(x_NEW_3, y_NEW_3, "-", linewidth=1, label="New code, Exc#3")
axs.plot(x_NEW_4, y_NEW_4, "-", linewidth=1, label="New code, Exc#4")
axs.get_yaxis().set_ticks([0])
axs.set(xlim=(0,11.4), xlabel="Frequency [THz]", ylabel="Spectral density", title=r"$J_{Dipti}\left( \omega \right) = \sum _\alpha \left( \frac {V^\alpha}{\omega_\alpha} \right)^2 \delta \left( \omega - \omega_\alpha \right)$")
axs.set_ylim(bottom=0, top=None)
axs.legend(loc="upper left")

'''
axins = inset_axes(axs, width="100%", height="100%", bbox_to_anchor=(.2, .5, .4, .3), bbox_transform=axs.transAxes)
axins.plot(x_inset, y_inset, "-k", linewidth=1)
axins.get_yaxis().set_ticks([])
axins.set(xlim=(0,11.4), ylabel="PDOS")
axins.set_ylim(bottom=0, top=None)
'''

fig.tight_layout()
fig.savefig("spectral_density_comparison_moreExc.pdf")
