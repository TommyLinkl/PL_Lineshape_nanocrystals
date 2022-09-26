import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.fft import fft, ifft
from scipy.integrate import simpson

fig, axs = plt.subplots(1,2, figsize=(9, 3.7))

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian800fs_Spectrum_QM_E17_T=300_aIndex0_bIndex0.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian800fs_Spectrum_QM_E17_T=300_aIndex0_bIndex0.dat", 'r') if line[0]!="#"])
x2 = np.array([float(line.strip().split()[0]) for line in \
    open("scale_coupling_coreShell/no_scaling_gaussian/results/gaussian800fs_Spectrum_QM_E17_T=300.dat", 'r') if line[0]!="#"])
y2 = np.array([float(line.strip().split()[1]) for line in \
    open("scale_coupling_coreShell/no_scaling_gaussian/results/gaussian800fs_Spectrum_QM_E17_T=300.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("experimental_digitized/single_NC/Hendrik_290K_shape.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("experimental_digitized/single_NC/Hendrik_290K_shape.dat", 'r') if line[0]!="#"])
axs[0].plot(x1*1000, y1, "-r", label="Lin adjusted, CdSe/CdS, 300K")
axs[0].plot(x2*1000, y2, "-k", label="Lin original, CdSe/CdS, 300K")
axs[0].plot(4050-x3, y3, "-", label="Hendrik, CdSe/CdS/ZnS, 290K")
axs[0].set(xlim=(1900,2200), ylim=(0, 1), xlabel="Energy [meV]", ylabel="Abs", title="PL line shape at high T")
axs[0].grid()
# axs[0].legend(bbox_to_anchor=(0, 1), loc="lower left", prop={"size": 6})
axs[0].legend(prop={"size": 6})

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/WidthFreqShift_vs_T_gaussian800fs.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[3]) for line in \
    open("scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/WidthFreqShift_vs_T_gaussian800fs.dat", 'r') if line[0]!="#"])
x2 = np.array([float(line.strip().split()[0]) for line in \
    open("scale_coupling_coreShell/no_scaling_gaussian/results/WidthFreqShift_vs_T_gaussianBroadening_800fs.dat", 'r') if line[0]!="#"])
y2 = np.array([float(line.strip().split()[4]) for line in \
    open("scale_coupling_coreShell/no_scaling_gaussian/results/WidthFreqShift_vs_T_gaussianBroadening_800fs.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("experimental_digitized/single_NC/Hendrik_Dot1_width.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("experimental_digitized/single_NC/Hendrik_Dot1_width.dat", 'r') if line[0]!="#"])
axs[1].plot(x1, y1*1000, "-+r", label="Lin adjusted, CdSe/CdS")
axs[1].plot(x2, y2*1000, "-ok", label="Lin original, CdSe/CdS")
axs[1].plot(x3, y3, "-s", label="Hendrik, CdSe/CdS/ZnS")
axs[1].set(xlabel="Temperature [K]", ylabel="Linewidth FWHM [meV]", title="PL Linewidth vs. Temp")
axs[1].grid()
# axs[1].legend(bbox_to_anchor=(0, 1), loc="lower left", prop={"size": 6})
axs[1].legend(prop={"size": 6})

fig.tight_layout()
fig.savefig("poster_plot_single_highT.pdf")
