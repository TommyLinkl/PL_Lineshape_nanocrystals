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
    open("../../CALCS/coreShell_3.0nmCdSe_2MLCdS/results/NEW_gaussian30fs_Spectrum_QM_E17_T=300.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("../../CALCS/coreShell_3.0nmCdSe_2MLCdS/results/NEW_gaussian30fs_Spectrum_QM_E17_T=300.dat", 'r') if line[0]!="#"])
x2 = np.array([float(line.strip().split()[0]) for line in \
    open("../../CALCS/coreShell_3.0nmCdSe_2MLCdS/results/OLD_gaussian30fs_Spectrum_QM_E17_T=300.dat", 'r') if line[0]!="#"])
y2 = np.array([float(line.strip().split()[1]) for line in \
    open("../../CALCS/coreShell_3.0nmCdSe_2MLCdS/results/OLD_gaussian30fs_Spectrum_QM_E17_T=300.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("../../experimental_digitized/ensemble/Meijerink_coreShell_295K_shape.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("../../experimental_digitized/ensemble/Meijerink_coreShell_295K_shape.dat", 'r') if line[0]!="#"])
axs[0].plot(x1*1000, y1, "-k", linewidth=1, label="NEW code, 30fs")
axs[0].plot(x2*1000, y2, "-.r", linewidth=1, label="OLD code, 30fs")
axs[0].plot(4095-x3, y3, ":b", linewidth=1, label="Meijerink JPCC") # 295K
axs[0].set(xlim=(1900,2200), ylim=(0, 1), xlabel="Energy [meV]", ylabel="Abs", title="QD ensemble line shape at high T (300K)")
axs[0].legend(loc="upper right", prop={"size": 10})


x1 = np.array([float(line.strip().split()[0]) for line in \
    open("../../CALCS/coreShell_3.0nmCdSe_2MLCdS/results/NEW_gaussian30fs_Spectrum_QM_E17_T=5.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("../../CALCS/coreShell_3.0nmCdSe_2MLCdS/results/NEW_gaussian30fs_Spectrum_QM_E17_T=5.dat", 'r') if line[0]!="#"])
x2 = np.array([float(line.strip().split()[0]) for line in \
    open("../../CALCS/coreShell_3.0nmCdSe_2MLCdS/results/OLD_gaussian30fs_Spectrum_QM_E17_T=5.dat", 'r') if line[0]!="#"])
y2 = np.array([float(line.strip().split()[1]) for line in \
    open("../../CALCS/coreShell_3.0nmCdSe_2MLCdS/results/OLD_gaussian30fs_Spectrum_QM_E17_T=5.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("../../experimental_digitized/ensemble/Meijerink_coreShell_4K_shape.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("../../experimental_digitized/ensemble/Meijerink_coreShell_4K_shape.dat", 'r') if line[0]!="#"])
axs[1].plot(x1*1000, y1, "-k", linewidth=1, label="NEW code, 30fs")
axs[1].plot(x2*1000, y2, "-.r", linewidth=1, label="OLD code, 30fs")
axs[1].plot(4150-x3, y3, ":b", linewidth=1, label="Meijerink JPCC") # 4K
axs[1].set(xlim=(1900,2200), ylim=(0, 1), xlabel="Energy [meV]", ylabel="Abs", title="QD ensemble line shape at low T (5K)")
axs[1].legend(loc="upper right", prop={"size": 10})


fig.tight_layout()
fig.savefig("manu_plot_ensemble.pdf")
