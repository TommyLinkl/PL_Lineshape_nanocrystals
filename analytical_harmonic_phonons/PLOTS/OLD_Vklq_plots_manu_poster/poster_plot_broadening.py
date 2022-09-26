import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.fft import fft, ifft
from scipy.integrate import simpson

fig, axs = plt.subplots(4, figsize=(5,12))

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("scale_coupling_coreShell/no_scaling_gaussian/results/gaussian200fs_Spectrum_QM_E17_T=0.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("scale_coupling_coreShell/no_scaling_gaussian/results/gaussian200fs_Spectrum_QM_E17_T=0.dat", 'r') if line[0]!="#"])
x2 = np.array([float(line.strip().split()[0]) for line in \
    open("scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian200fs_Spectrum_QM_E17_T=0_aIndex0_bIndex0.dat", 'r') if line[0]!="#"])
y2 = np.array([float(line.strip().split()[1]) for line in \
    open("scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian200fs_Spectrum_QM_E17_T=0_aIndex0_bIndex0.dat", 'r') if line[0]!="#"])
axs[0].plot(x1, y1, "-k", label="Original couplings")
axs[0].plot(x2, y2, "-r", label="Adjusted couplings")
axs[0].set(xlim=(2,2.15), ylim=(-0.1, 1.0), ylabel="Abs", title="Guassian 200fs, 0K, CdSe/CdS")
axs[0].legend()
axs[0].grid()

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("scale_coupling_coreShell/no_scaling_gaussian/results/gaussian800fs_Spectrum_QM_E17_T=0.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("scale_coupling_coreShell/no_scaling_gaussian/results/gaussian800fs_Spectrum_QM_E17_T=0.dat", 'r') if line[0]!="#"])
x2 = np.array([float(line.strip().split()[0]) for line in \
    open("scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian800fs_Spectrum_QM_E17_T=0_aIndex0_bIndex0.dat", 'r') if line[0]!="#"])
y2 = np.array([float(line.strip().split()[1]) for line in \
    open("scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian800fs_Spectrum_QM_E17_T=0_aIndex0_bIndex0.dat", 'r') if line[0]!="#"])
axs[1].plot(x1, y1, "-k", label="Original couplings")
axs[1].plot(x2, y2, "-r", label="Adjusted couplings")
axs[1].set(xlim=(2,2.15), ylim=(-0.1, 1.0), ylabel="Abs", title="Guassian 800fs, 0K")
axs[1].grid()

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("scale_coupling_coreShell/no_scaling_gaussian/results/gaussian4ps_Spectrum_QM_E17_T=0.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("scale_coupling_coreShell/no_scaling_gaussian/results/gaussian4ps_Spectrum_QM_E17_T=0.dat", 'r') if line[0]!="#"])
x2 = np.array([float(line.strip().split()[0]) for line in \
    open("scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian4ps_Spectrum_QM_E17_T=0_aIndex0_bIndex0.dat", 'r') if line[0]!="#"])
y2 = np.array([float(line.strip().split()[1]) for line in \
    open("scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian4ps_Spectrum_QM_E17_T=0_aIndex0_bIndex0.dat", 'r') if line[0]!="#"])
axs[2].plot(x1, y1, "-k", label="Original couplings")
axs[2].plot(x2, y2, "-r", label="Adjusted couplings")
axs[2].set(xlim=(2,2.15), ylim=(-0.1, 1.0), ylabel="Abs", title="Guassian 4ps, 0K")
axs[2].grid()

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("scale_coupling_coreShell/no_scaling_gaussian/results/Spectrum_QM_E17_T=0_shrink#1_highRes.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("scale_coupling_coreShell/no_scaling_gaussian/results/Spectrum_QM_E17_T=0_shrink#1_highRes.dat", 'r') if line[0]!="#"])
x2 = np.array([float(line.strip().split()[0]) for line in \
    open("scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/highRes_Spectrum_QM_E17_T=0_aIndex0_bIndex0.dat", 'r') if line[0]!="#"])
y2 = np.array([float(line.strip().split()[1]) for line in \
    open("scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/highRes_Spectrum_QM_E17_T=0_aIndex0_bIndex0.dat", 'r') if line[0]!="#"])
axs[3].plot(x1, y1, "-k", label="Original couplings")
axs[3].plot(x2, y2, "-r", label="Adjusted couplings")
axs[3].set(xlim=(2,2.15), ylim=(-0.1, 1.0), xlabel="Energy [eV]", ylabel="Abs", title="Guassian Inf, i.e. highRes, 0K")
axs[3].grid()
axs[3].legend()

fig.tight_layout()
fig.savefig("poster_plot_broadening.pdf")
