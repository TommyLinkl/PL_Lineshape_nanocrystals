import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib import rc
from allFunc import * 
matplotlib.rcParams['font.family'] = "serif"
matplotlib.rcParams['font.sans-serif'] = "Times New Roman"
matplotlib.rcParams['font.size'] = 12
matplotlib.rcParams['mathtext.fontset'] = "stix"
rc('axes', linewidth=0.4)
matplotlib.rcParams['xtick.major.width'] = 0.4
matplotlib.rcParams['ytick.major.width'] = 0.4

fig, axs = plt.subplots(1,2, figsize=(9, 4.5))

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_2MLCdS/results/NEW_gaussian40fs_Spectrum_QM_E17_T=4.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_2MLCdS/results/NEW_gaussian40fs_Spectrum_QM_E17_T=4.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("../experimental_digitized/ensemble/Meijerink_coreShell_4K_shape.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("../experimental_digitized/ensemble/Meijerink_coreShell_4K_shape.dat", 'r') if line[0]!="#"])
axs[0].plot(x1*1000, y1, "-k", linewidth=1, label="Model")
axs[0].plot(4149-x3, y3, ":b", linewidth=1, label="van der Bok et al.")
axs[0].set(xlim=(1900,2250), ylim=(0, 1), xlabel="Energy [meV]", ylabel="Abs")
axs[0].legend(loc="upper right", prop={"size": 10})
axs[0].text(1905, 0.95, "Ensemble, 4K")
axs[0].get_yaxis().set_ticks([0.0, 0.5, 1.0])
axs[0].get_xaxis().set_ticks([1900, 2000, 2100, 2200])

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_2MLCdS/results/NEW_exponential30fs_gaussian40fs_Spectrum_QM_E17_T=300.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_2MLCdS/results/NEW_exponential30fs_gaussian40fs_Spectrum_QM_E17_T=300.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("../experimental_digitized/ensemble/Meijerink_coreShell_295K_shape.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("../experimental_digitized/ensemble/Meijerink_coreShell_295K_shape.dat", 'r') if line[0]!="#"])
axs[1].plot(x1*1000, y1, "-k", linewidth=1, label="Model")
axs[1].plot(4097-x3, y3, ":b", linewidth=1, label="van der Bok et al.")
axs[1].set(xlim=(1900,2250), ylim=(0, 1), xlabel="Energy [meV]")
# axs[1].legend(loc="upper right", prop={"size": 10})
axs[1].text(1905, 0.95, "Ensemble, 300K")
axs[1].get_yaxis().set_ticks([0.0, 0.5, 1.0])
axs[1].get_xaxis().set_ticks([1900, 2000, 2100, 2200])


fig.tight_layout()
fig.savefig("manu_plot_ensemble.pdf")
