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
    open("../scale_coupling_coreShell/nonuniform_scaling_trial7/case77_T/results/gaussian50fs_Spectrum_QM_E17_T=300_case77.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("../scale_coupling_coreShell/nonuniform_scaling_trial7/case77_T/results/gaussian50fs_Spectrum_QM_E17_T=300_case77.dat", 'r') if line[0]!="#"])
x2 = np.array([float(line.strip().split()[0]) for line in \
    open("../scale_coupling_coreShell/no_scaling_gaussian/results/gaussian50fs_Spectrum_QM_E17_T=300.dat", 'r') if line[0]!="#"])
y2 = np.array([float(line.strip().split()[1]) for line in \
    open("../scale_coupling_coreShell/no_scaling_gaussian/results/gaussian50fs_Spectrum_QM_E17_T=300.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("../experimental_digitized/ensemble/Meijerink_coreShell_295K_shape.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("../experimental_digitized/ensemble/Meijerink_coreShell_295K_shape.dat", 'r') if line[0]!="#"])
axs[0].plot(x2*1000, y2, "-.r", linewidth=1, label="Short-range model")
axs[0].plot(x1*1000+10, y1, "-k", linewidth=1, label="Frohlich model")
axs[0].plot(4100-x3, y3, ":b", linewidth=1, label="Meijerink JPCC") # 295K
axs[0].set(xlim=(1900,2300), ylim=(0, 1), xlabel="Energy [meV]", ylabel="Abs", title="QD ensemble line shape at high T (300K)")
# axs[0].legend(bbox_to_anchor=(0, 1), loc="lower left", prop={"size": 6})
axs[0].legend(loc="upper right", prop={"size": 10})

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("../scale_coupling_coreShell/nonuniform_scaling_trial7/case77_T/results/gaussian50fs_Spectrum_QM_E17_T=5_case77.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("../scale_coupling_coreShell/nonuniform_scaling_trial7/case77_T/results/gaussian50fs_Spectrum_QM_E17_T=5_case77.dat", 'r') if line[0]!="#"])
x2 = np.array([float(line.strip().split()[0]) for line in \
    open("../scale_coupling_coreShell/no_scaling_gaussian/results/gaussian50fs_Spectrum_QM_E17_T=5.dat", 'r') if line[0]!="#"])
y2 = np.array([float(line.strip().split()[1]) for line in \
    open("../scale_coupling_coreShell/no_scaling_gaussian/results/gaussian50fs_Spectrum_QM_E17_T=5.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("../experimental_digitized/ensemble/Meijerink_coreShell_4K_shape.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("../experimental_digitized/ensemble/Meijerink_coreShell_4K_shape.dat", 'r') if line[0]!="#"])
axs[1].plot(x2*1000, y2, "-.r", linewidth=1, label="Short-range model")
axs[1].plot(x1*1000+10, y1, "-k", linewidth=1, label="Frohlich model")
axs[1].plot(4150-x3, y3, ":b", linewidth=1, label="Meijerink JPCC") # 4K
axs[1].set(xlim=(1900,2300), ylim=(0, 1), xlabel="Energy [meV]", ylabel="Abs", title="QD ensemble line shape at low T (5K)")
# axs[1].legend(bbox_to_anchor=(0, 1), loc="lower left", prop={"size": 6})
axs[1].legend(loc="upper right", prop={"size": 10})


'''
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
# axs[2].plot(x1, y1*1000, "-+r", label="Lin adjusted")
axs[2].plot(x2, y2*1000, "-ok", label="Lin original")
axs[2].plot(x3, y3, "-s", label="Meijerink JPCC, CdSe/CdS ensemble")
axs[2].set(ylim=(0,110), xlabel="Temperature [K]", ylabel="Linewidth FWHM [meV]", title="Inhomogeneous PL Linewidth")
axs[2].grid()
# axs[2].legend(bbox_to_anchor=(0, 1), loc="lower left", prop={"size": 6})
axs[2].legend(loc="lower right", prop={"size": 6})
'''

fig.tight_layout()
fig.savefig("manu_plot_ensemble.pdf")
