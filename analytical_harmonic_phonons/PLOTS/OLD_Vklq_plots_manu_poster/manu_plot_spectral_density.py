import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib import rc
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
phonons = np.array([float(line.strip().split()[1]) for line in open("../CALCS/3nmCdSe-3MLCdS_a5/w.dat", 'r') if line[0]!="#"])
short_V = np.array([float(line.strip().split()[3]) for line in open("../CALCS/3nmCdSe-3MLCdS_a5/Vklq_trim.dat", 'r') if line[0]!="#"])
frohlich_V = np.array([float(line.strip().split()[3]) for line in \
 open("../scale_coupling_coreShell/nonuniform_scaling_trial7/results/Vklq_case77.dat", 'r') if line[0]!="#"])
phonons = phonons[6:]
short_V = short_V[6:]
frohlich_V = frohlich_V[6:]

(x_short, y_short) = Dipti_spectral_density(phonons, short_V, gaussian_std)
(x_frohlich, y_frohlich) = Dipti_spectral_density(phonons, frohlich_V, gaussian_std)
(x_inset, y_inset) = PDOS(phonons, gaussian_std)


fig, axs = plt.subplots(1,1, figsize=(6, 4))


axs.plot(x_short, y_short, "-.r", linewidth=1, label="Short-range model")
axs.plot(x_frohlich, y_frohlich, "-k", linewidth=1, label="Frohlich model")
axs.get_yaxis().set_ticks([0])
axs.set(xlim=(0,11.4), xlabel="Frequency [THz]", ylabel=r"Spectral density $J\left( \omega \right) = \sum _\alpha \left( \frac {V^\alpha}{\omega_\alpha} \right)^2 \delta \left( \omega - \omega_\alpha \right)$", title="3nm CdSe/3ML CdS core-shell quantum dot")
axs.set_ylim(bottom=0, top=None)
axs.legend(loc="upper left")

axins = inset_axes(axs, width="100%", height="100%", bbox_to_anchor=(.35, .5, .4, .3), bbox_transform=axs.transAxes)
axins.plot(x_inset, y_inset, "-k", linewidth=1)
axins.get_yaxis().set_ticks([])
axins.set(xlim=(0,11.4), ylabel="PDOS")
axins.set_ylim(bottom=0, top=None)

# fig.tight_layout()
fig.savefig("manu_plot_spectral_density.pdf")
