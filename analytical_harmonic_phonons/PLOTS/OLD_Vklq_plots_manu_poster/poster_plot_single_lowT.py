import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.fft import fft, ifft
from scipy.integrate import simpson

fig, axs = plt.subplots(1,4, figsize=(12, 3.85))

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian800fs_Spectrum_QM_E17_T=5_aIndex0_bIndex0.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian800fs_Spectrum_QM_E17_T=5_aIndex0_bIndex0.dat", 'r') if line[0]!="#"])
x2 = np.array([float(line.strip().split()[0]) for line in \
    open("experimental_digitized/single_NC/Hendrik_4K_reversedAndShifted_shape.dat", 'r') if line[0]!="#"])
y2 = np.array([float(line.strip().split()[1]) for line in \
    open("experimental_digitized/single_NC/Hendrik_4K_reversedAndShifted_shape.dat", 'r') if line[0]!="#"])
axs[0].plot(x1*1000, y1, "-r", label="Lin adjusted, CdSe/CdS, 5K")
axs[0].plot(x2*1000-2, y2, "-", label="Utzat, CdSe/CdS/ZnS, 4K")
axs[0].set(xlim=(2000,2120), ylim=(0, 1), xlabel="Energy [meV]", ylabel="Abs")
axs[0].grid()
axs[0].legend(bbox_to_anchor=(0, 1), loc="lower left", prop={"size": 6})

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian150fs_Spectrum_QM_E17_T=5_aIndex0_bIndex0.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian150fs_Spectrum_QM_E17_T=5_aIndex0_bIndex0.dat", 'r') if line[0]!="#"])
x2 = np.array([float(line.strip().split()[0]) for line in \
    open("experimental_digitized/single_NC/Bawendi_PRB76_2007_singleNC_shape.dat", 'r') if line[0]!="#"])
y2 = np.array([float(line.strip().split()[1]) for line in \
    open("experimental_digitized/single_NC/Bawendi_PRB76_2007_singleNC_shape.dat", 'r') if line[0]!="#"])
axs[1].plot(x1*1000, y1, "-r", label="Lin adjusted, CdSe/CdS, 5K")
axs[1].plot(4417-x2, y2, "-", label="Bawendi PRB76, CdSe/ZnS, 5K")
axs[1].set(xlim=(2000,2120), ylim=(0, 1), xlabel="Energy [meV]", ylabel="Abs")
axs[1].grid()
axs[1].legend(bbox_to_anchor=(0, 1), loc="lower left", prop={"size": 6})

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian1ps_Spectrum_QM_E17_T=10_aIndex0_bIndex0.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian1ps_Spectrum_QM_E17_T=10_aIndex0_bIndex0.dat", 'r') if line[0]!="#"])
x2 = np.array([float(line.strip().split()[0]) for line in \
    open("experimental_digitized/single_NC/Klimov_PRL93_2004_Dot1_shape.dat", 'r') if line[0]!="#"])
y2 = np.array([float(line.strip().split()[1]) for line in \
    open("experimental_digitized/single_NC/Klimov_PRL93_2004_Dot1_shape.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("experimental_digitized/single_NC/Klimov_PRL93_2004_Dot2_shape.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("experimental_digitized/single_NC/Klimov_PRL93_2004_Dot2_shape.dat", 'r') if line[0]!="#"])
axs[2].plot(x1*1000, y1, "-r", label="Lin adjusted, CdSe/CdS, 10K")
axs[2].plot(2015.0212+x2, y2-0.08, "-", label="Klimov PRL93, CdSe/ZnS Dot1, 10K")
axs[2].plot(2015.0212+x3, y3-0.04+0.6, "-", label="Klimov PRL93, CdSe/ZnS Dot2, 10K")
axs[2].set(xlim=(2000,2120), ylim=(0,1.5), xlabel="Energy [meV]", ylabel="Abs")
axs[2].grid()
axs[2].legend(bbox_to_anchor=(0, 1), loc="lower left", prop={"size": 6})

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian150fs_Spectrum_QM_E17_T=10_aIndex0_bIndex0.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian150fs_Spectrum_QM_E17_T=10_aIndex0_bIndex0.dat", 'r') if line[0]!="#"])
x2 = np.array([float(line.strip().split()[0]) for line in \
    open("experimental_digitized/single_NC/PRL85_Dot1_shape.dat", 'r') if line[0]!="#"])
y2 = np.array([float(line.strip().split()[1]) for line in \
    open("experimental_digitized/single_NC/PRL85_Dot1_shape.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("experimental_digitized/single_NC/PRL85_Dot2_shape.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("experimental_digitized/single_NC/PRL85_Dot2_shape.dat", 'r') if line[0]!="#"])
x4 = np.array([float(line.strip().split()[0]) for line in \
    open("experimental_digitized/single_NC/PRL85_Dot3_shape.dat", 'r') if line[0]!="#"])
y4 = np.array([float(line.strip().split()[1]) for line in \
    open("experimental_digitized/single_NC/PRL85_Dot3_shape.dat", 'r') if line[0]!="#"])
axs[3].plot(x1*1000, y1, "-r", label="Lin adjusted, CdSe/CdS, 10K")
axs[3].plot(4106-x2, y2/0.853094, "-", label="Bawendi PRL85, CdSe/ZnS Dot1, 10K")
axs[3].plot(4115-x3, y3/0.853094+0.5, "-", label="Bawendi PRL85, CdSe/ZnS Dot2, 10K")
axs[3].plot(4128-x4, y4/0.787948+1, "-", label="Bawendi PRL85, CdSe/ZnS Dot3, 10K")
axs[3].set(xlim=(2000,2120), ylim=(0,2), xlabel="Energy [meV]", ylabel="Abs")
axs[3].grid()
axs[3].legend(bbox_to_anchor=(0, 1), loc="lower left", prop={"size": 6})

fig.tight_layout()
fig.savefig("poster_plot_single_lowT.pdf")
