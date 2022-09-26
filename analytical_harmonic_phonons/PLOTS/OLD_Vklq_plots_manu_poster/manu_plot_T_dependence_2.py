import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
matplotlib.rcParams['font.family'] = "serif"
matplotlib.rcParams['font.sans-serif'] = "Times New Roman"
rc('axes', linewidth=0.4)
matplotlib.rcParams['xtick.major.width'] = 0.4
matplotlib.rcParams['ytick.major.width'] = 0.4

fig, axs = plt.subplots(1,2, figsize=(8, 6))

x1 = np.array([float(line.strip().split()[0]) for line in open("../CALCS/coreShell_3.0nmCdSe_2MLCdS/results/NEW_WidthFreqShift_vs_T_gaussian_30fs.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in open("../CALCS/coreShell_3.0nmCdSe_2MLCdS/results/NEW_WidthFreqShift_vs_T_gaussian_30fs.dat", 'r') if line[0]!="#"])
x2 = np.array([float(line.strip().split()[0]) for line in open("../CALCS/coreShell_3.0nmCdSe_2MLCdS/results/NEW_WidthFreqShift_vs_T_gaussian_20fs.dat", 'r') if line[0]!="#"])
y2 = np.array([float(line.strip().split()[1]) for line in open("../CALCS/coreShell_3.0nmCdSe_2MLCdS/results/NEW_WidthFreqShift_vs_T_gaussian_20fs.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in open("../experimental_digitized/ensemble/Meijerink_coreShell_inhomo_width.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in open("../experimental_digitized/ensemble/Meijerink_coreShell_inhomo_width.dat", 'r') if line[0]!="#"])
axs[0].plot(x1, y1*1000, "-+k", linewidth=1, label="Gaussian 30fs")
# axs[0].plot(x2, y2*1000, "-.+r", label="Gaussian 20fs")
axs[0].plot(x3, y3, ":sb", linewidth=1, label="Meijerink JPCC")
axs[0].set(ylim=(0,110), xlabel="Temperature [K]", ylabel="FWHM [meV]", title="Inhomogeneous PL Linewidth")
axs[0].legend(loc="lower right", prop={"size": 10})

x1 = np.array([float(line.strip().split()[0]) for line in open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_WidthFreqShift_vs_T_gaussian_500fs.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_WidthFreqShift_vs_T_gaussian_500fs.dat", 'r') if line[0]!="#"])
x2 = np.array([float(line.strip().split()[0]) for line in open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_WidthFreqShift_vs_T_gaussian_30fs.dat", 'r') if line[0]!="#"])
y2 = np.array([float(line.strip().split()[1]) for line in open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_WidthFreqShift_vs_T_gaussian_30fs.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in open("../experimental_digitized/single_NC/Hendrik_Dot1_width.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in open("../experimental_digitized/single_NC/Hendrik_Dot1_width.dat", 'r') if line[0]!="#"])
axs[1].plot(x1, y1*1000, "-+k", linewidth=1, label="Gaussian 500fs")
axs[1].plot(300, 41.729, "+r", linewidth=1, label="Finite ph lifetime 10 ps")
axs[1].plot(300, 42.018, "+m", linewidth=1, label="Finite ph lifetime 10 ps + gaussian")
# axs[1].plot(x2, y2*1000, "-.+r", label="Gaussian 30fs")
axs[1].plot(x3, y3, ":og", linewidth=1, label="Utzat, CdSe/CdS/ZnS")
axs[1].set(ylim=(0,110), xlabel="Temperature [K]", ylabel="FWHM [meV]", title="Single-NC PL Linewidth vs. Temp")
axs[1].legend(loc="upper right", prop={"size": 10})

fig.tight_layout()
fig.savefig("manu_plot_T_dependence.pdf")
