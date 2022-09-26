import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
matplotlib.rcParams['font.family'] = "serif"
matplotlib.rcParams['font.sans-serif'] = "Times New Roman"
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
axs[0].plot(x1*1000, y1, "-k", linewidth=1, label="Gauss40fs, T2=Inf")
axs[0].plot(4149-x3, y3, ":b", linewidth=1, label="Meijerink JPCC")
axs[0].set(xlim=(1900,2250), ylim=(0, 1), xlabel="Energy [meV]", ylabel="Abs", title="3nmCdSe/2MLCdS QD ensemble line shape at low T (4K)")
axs[0].legend(loc="upper right", prop={"size": 10})

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_2MLCdS/results/NEW_exponential30fs_gaussian40fs_Spectrum_QM_E17_T=300.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_2MLCdS/results/NEW_exponential30fs_gaussian40fs_Spectrum_QM_E17_T=300.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("../experimental_digitized/ensemble/Meijerink_coreShell_295K_shape.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("../experimental_digitized/ensemble/Meijerink_coreShell_295K_shape.dat", 'r') if line[0]!="#"])
axs[1].plot(x1*1000, y1, "-k", linewidth=1, label="Gauss40fs, T2=30fs")
axs[1].plot(4097-x3, y3, ":b", linewidth=1, label="Meijerink JPCC")
axs[1].set(xlim=(1900,2250), ylim=(0, 1), xlabel="Energy [meV]", ylabel="Abs", title="At high T (300K)")
axs[1].legend(loc="upper right", prop={"size": 10})


fig.tight_layout()
fig.savefig("TEST_manu_plot_ensemble.pdf")
