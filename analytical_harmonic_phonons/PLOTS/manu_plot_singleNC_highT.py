import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import FormatStrFormatter
from matplotlib import rc
matplotlib.rcParams['font.family'] = "serif"
matplotlib.rcParams['font.sans-serif'] = "Times New Roman"
matplotlib.rcParams['mathtext.fontset'] = "stix"
rc('axes', linewidth=0.4)
matplotlib.rcParams['xtick.major.width'] = 0.4
matplotlib.rcParams['ytick.major.width'] = 0.4
fig, axs = plt.subplots(1,2, figsize=(8, 4))

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_exponential200fs_gaussian500fs_Spectrum_QM_E17_T=220.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_exponential200fs_gaussian500fs_Spectrum_QM_E17_T=220.dat", 'r') if line[0]!="#"])
x2 = np.array([float(line.strip().split()[0]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_exponential200fs_Spectrum_QM_E17_T=220.dat", 'r') if line[0]!="#"])
y2 = np.array([float(line.strip().split()[1]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_exponential200fs_Spectrum_QM_E17_T=220.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("../experimental_digitized/single_NC/Hendrik_220K_shape.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("../experimental_digitized/single_NC/Hendrik_220K_shape.dat", 'r') if line[0]!="#"])
axs[0].plot(x1*1000, y1, "-k", linewidth=1, label="Gaussian500fs, exp200fs")
axs[0].plot(x2*1000, y2, "-r", linewidth=0.5, label="Exp200fs")
axs[0].plot(4030-x3, y3, "g", linewidth=1, label="Utzat")
axs[0].set(xlim=(1900,2150), ylim=(0, 1), xlabel="Energy [meV]", ylabel="Abs", title="3nm/3ML single-NC line shape at 220K")
axs[0].legend(loc="upper right", prop={"size": 10})


x1 = np.array([float(line.strip().split()[0]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_exponential200fs_gaussian500fs_Spectrum_QM_E17_T=290.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_exponential200fs_gaussian500fs_Spectrum_QM_E17_T=290.dat", 'r') if line[0]!="#"])
x2 = np.array([float(line.strip().split()[0]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_exponential200fs_Spectrum_QM_E17_T=290.dat", 'r') if line[0]!="#"])
y2 = np.array([float(line.strip().split()[1]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_exponential200fs_Spectrum_QM_E17_T=290.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("../experimental_digitized/single_NC/Hendrik_290K_shape.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("../experimental_digitized/single_NC/Hendrik_290K_shape.dat", 'r') if line[0]!="#"])
axs[1].plot(x1*1000, y1, "-k", linewidth=1, label="Gaussian500fs, exp200fs")
axs[1].plot(x2*1000, y2, "-r", linewidth=0.5, label="Exp200fs")
axs[1].plot(4002-x3, y3, "g", linewidth=1, label="Utzat, CdSe/CdS/ZnS")
axs[1].set(xlim=(1900,2150), ylim=(0, 1), xlabel="Energy [meV]", title="290K")
axs[1].legend(loc="upper right", prop={"size": 10})

fig.tight_layout()
fig.savefig("manu_plot_singleNC_highT.pdf")
