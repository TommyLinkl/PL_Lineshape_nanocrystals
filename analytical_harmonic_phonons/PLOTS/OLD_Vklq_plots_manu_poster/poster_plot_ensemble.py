import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.fft import fft, ifft
from scipy.integrate import simpson

fig, axs = plt.subplots(1,3, figsize=(12, 2.5))

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian50fs_Spectrum_QM_E17_T=5_aIndex0_bIndex0.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian50fs_Spectrum_QM_E17_T=5_aIndex0_bIndex0.dat", 'r') if line[0]!="#"])
x2 = np.array([float(line.strip().split()[0]) for line in \
    open("scale_coupling_coreShell/no_scaling_gaussian/results/gaussian50fs_Spectrum_QM_E17_T=5.dat", 'r') if line[0]!="#"])
y2 = np.array([float(line.strip().split()[1]) for line in \
    open("scale_coupling_coreShell/no_scaling_gaussian/results/gaussian50fs_Spectrum_QM_E17_T=5.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("experimental_digitized/ensemble/Meijerink_coreShell_4K_shape.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("experimental_digitized/ensemble/Meijerink_coreShell_4K_shape.dat", 'r') if line[0]!="#"])
axs[0].plot(x1*1000, y1, "-r", label="Lin adjusted, 5K")
axs[0].plot(x2*1000, y2, "-k", label="Lin original, 5K")
axs[0].plot(4150-x3, y3, "-", label="Meijerink JPCC, 4K")
axs[0].set(xlim=(1900,2200), ylim=(0, 1), xlabel="Energy [meV]", ylabel="Abs", title="Low T line shape, ensemble")
axs[0].grid()
# axs[0].legend(bbox_to_anchor=(0, 1), loc="lower left", prop={"size": 6})
axs[0].legend(loc="upper right", prop={"size": 6})

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian50fs_Spectrum_QM_E17_T=300_aIndex0_bIndex0.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian50fs_Spectrum_QM_E17_T=300_aIndex0_bIndex0.dat", 'r') if line[0]!="#"])
x2 = np.array([float(line.strip().split()[0]) for line in \
    open("scale_coupling_coreShell/no_scaling_gaussian/results/gaussian50fs_Spectrum_QM_E17_T=300.dat", 'r') if line[0]!="#"])
y2 = np.array([float(line.strip().split()[1]) for line in \
    open("scale_coupling_coreShell/no_scaling_gaussian/results/gaussian50fs_Spectrum_QM_E17_T=300.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("experimental_digitized/ensemble/Meijerink_coreShell_295K_shape.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("experimental_digitized/ensemble/Meijerink_coreShell_295K_shape.dat", 'r') if line[0]!="#"])
axs[1].plot(x1*1000, y1, "-r", label="Lin adjusted, 300K")
axs[1].plot(x2*1000, y2, "-k", label="Lin original, 300K")
axs[1].plot(4100-x3, y3, "-", label="Meijerink JPCC, 295K")
axs[1].set(xlim=(1900,2200), ylim=(0, 1), xlabel="Energy [meV]", ylabel="Abs", title="High T line shape, ensemble")
axs[1].grid()
# axs[1].legend(bbox_to_anchor=(0, 1), loc="lower left", prop={"size": 6})
axs[1].legend(loc="upper right", prop={"size": 6})

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/WidthFreqShift_vs_T_gaussian50fs.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[3]) for line in \
    open("scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/WidthFreqShift_vs_T_gaussian50fs.dat", 'r') if line[0]!="#"])
x2 = np.array([float(line.strip().split()[0]) for line in \
    open("scale_coupling_coreShell/no_scaling_gaussian/results/WidthFreqShift_vs_T_gaussianBroadening_50fs.dat", 'r') if line[0]!="#"])
y2 = np.array([float(line.strip().split()[4]) for line in \
    open("scale_coupling_coreShell/no_scaling_gaussian/results/WidthFreqShift_vs_T_gaussianBroadening_50fs.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("experimental_digitized/ensemble/Meijerink_coreShell_inhomo_width.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("experimental_digitized/ensemble/Meijerink_coreShell_inhomo_width.dat", 'r') if line[0]!="#"])
axs[2].plot(x1, y1*1000, "-+r", label="Lin adjusted")
axs[2].plot(x2, y2*1000, "-ok", label="Lin original")
axs[2].plot(x3, y3, "-s", label="Meijerink JPCC, CdSe/CdS ensemble")
axs[2].set(xlabel="Temperature [K]", ylabel="Linewidth FWHM [meV]", title="Inhomogeneous PL Linewidth")
axs[2].grid()
# axs[2].legend(bbox_to_anchor=(0, 1), loc="lower left", prop={"size": 6})
axs[2].legend(loc="lower right", prop={"size": 6})

fig.tight_layout()
fig.savefig("poster_plot_ensemble.pdf")
