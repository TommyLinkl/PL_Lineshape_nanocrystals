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

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("../scale_coupling_coreShell/nonuniform_scaling_trial7/case77_T/results/WidthFreqShift_vs_T_gaussian_50fs.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[8]) for line in 
    open("../scale_coupling_coreShell/nonuniform_scaling_trial7/case77_T/results/WidthFreqShift_vs_T_gaussian_50fs.dat", 'r') if line[0]!="#"])
x2 = np.array([float(line.strip().split()[0]) for line in \
    open("../scale_coupling_coreShell/no_scaling_gaussian/results/WidthFreqShift_vs_T_gaussianBroadening_50fs.dat", 'r') if line[0]!="#"])
y2 = np.array([float(line.strip().split()[4]) for line in \
    open("../scale_coupling_coreShell/no_scaling_gaussian/results/WidthFreqShift_vs_T_gaussianBroadening_50fs.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("../experimental_digitized/ensemble/Meijerink_coreShell_inhomo_width.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("../experimental_digitized/ensemble/Meijerink_coreShell_inhomo_width.dat", 'r') if line[0]!="#"])
axs[0].plot(x1, y1*1000, "-+k", linewidth=1, label="Frohlich model")
axs[0].plot(x2, y2*1000, "-.+r", label="Short-range model")
axs[0].plot(x3, y3, ":sb", linewidth=1, label="Meijerink JPCC")
axs[0].set(ylim=(0,110), xlabel="Temperature [K]", ylabel="FWHM [meV]", title="Inhomogeneous PL Linewidth")
axs[0].legend(loc="lower right", prop={"size": 10})


x1 = np.array([float(line.strip().split()[0]) for line in \
    open("../scale_coupling_coreShell/nonuniform_scaling_trial7/case77_T/results/WidthFreqShift_vs_T_gaussian_800fs.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[8]) for line in \
    open("../scale_coupling_coreShell/nonuniform_scaling_trial7/case77_T/results/WidthFreqShift_vs_T_gaussian_800fs.dat", 'r') if line[0]!="#"])
x2 = np.array([float(line.strip().split()[0]) for line in \
    open("../scale_coupling_coreShell/no_scaling_gaussian/results/WidthFreqShift_vs_T_gaussianBroadening_800fs.dat", 'r') if line[0]!="#"])
y2 = np.array([float(line.strip().split()[4]) for line in \
    open("../scale_coupling_coreShell/no_scaling_gaussian/results/WidthFreqShift_vs_T_gaussianBroadening_800fs.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("../experimental_digitized/single_NC/Hendrik_Dot1_width.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("../experimental_digitized/single_NC/Hendrik_Dot1_width.dat", 'r') if line[0]!="#"])
axs[1].plot(x1, y1*1000, "-+k", linewidth=1, label="Frohlich model")
axs[1].plot(x2, y2*1000, "-.+r", label="Short-range model")
axs[1].plot(x3, y3, ":og", linewidth=1, label="Utzat, CdSe/CdS/ZnS")
axs[1].set(ylim=(0,110), xlabel="Temperature [K]", ylabel="FWHM [meV]", title="Single-NC PL Linewidth vs. Temp")
axs[1].legend(loc="lower right", prop={"size": 10})

fig.tight_layout()
fig.savefig("manu_plot_T_dependence.pdf")
