import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import FormatStrFormatter
from matplotlib import rc
matplotlib.rcParams['font.family'] = "serif"
matplotlib.rcParams['font.sans-serif'] = "Times New Roman"
rc('axes', linewidth=0.4)
matplotlib.rcParams['xtick.major.width'] = 0.4
matplotlib.rcParams['ytick.major.width'] = 0.4
fig, axs = plt.subplots(1,2, figsize=(12, 6))

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("../scale_coupling_coreShell/nonuniform_scaling_trial7/case77_T/results/gaussian800fs_Spectrum_QM_E17_T=300_case77.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("../scale_coupling_coreShell/nonuniform_scaling_trial7/case77_T/results/gaussian800fs_Spectrum_QM_E17_T=300_case77.dat", 'r') if line[0]!="#"])
x2 = np.array([float(line.strip().split()[0]) for line in \
    open("../scale_coupling_coreShell/no_scaling_gaussian/results/gaussian800fs_Spectrum_QM_E17_T=300.dat", 'r') if line[0]!="#"])
y2 = np.array([float(line.strip().split()[1]) for line in \
    open("../scale_coupling_coreShell/no_scaling_gaussian/results/gaussian800fs_Spectrum_QM_E17_T=300.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("../experimental_digitized/single_NC/Hendrik_290K_shape.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("../experimental_digitized/single_NC/Hendrik_290K_shape.dat", 'r') if line[0]!="#"])
axs[0].plot(x1*1000, y1, "-k", linewidth=1, label="Frohlich model")
axs[0].plot(x2*1000, y2, "-.r", linewidth=1, label="Short-range model")
axs[0].plot(4050-x3, y3, "g", linewidth=1, label="Utzat, CdSe/CdS/ZnS") # 290K
axs[0].set(xlim=(1900,2200), ylim=(0, 1), xlabel="Energy [meV]", ylabel="Abs", title="Single-NC line shape at high T (300K)")
# axs[0].legend(bbox_to_anchor=(0, 1), loc="lower left", prop={"size": 10})
axs[0].legend(loc="upper right", prop={"size": 10})


x1 = np.array([float(line.strip().split()[0]) for line in \
    open("../scale_coupling_coreShell/nonuniform_scaling_trial7/case77_T/results/gaussian800fs_Spectrum_QM_E17_T=5_case77.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("../scale_coupling_coreShell/nonuniform_scaling_trial7/case77_T/results/gaussian800fs_Spectrum_QM_E17_T=5_case77.dat", 'r') if line[0]!="#"])
x2 = np.array([float(line.strip().split()[0]) for line in \
    open("../scale_coupling_coreShell/no_scaling_gaussian/results/gaussian800fs_Spectrum_QM_E17_T=5.dat", 'r') if line[0]!="#"])
y2 = np.array([float(line.strip().split()[1]) for line in \
    open("../scale_coupling_coreShell/no_scaling_gaussian/results/gaussian800fs_Spectrum_QM_E17_T=5.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("../experimental_digitized/single_NC/Hendrik_4K_reversedAndShifted_shape.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("../experimental_digitized/single_NC/Hendrik_4K_reversedAndShifted_shape.dat", 'r') if line[0]!="#"])
x4 = np.array([float(line.strip().split()[0]) for line in \
    open("../scale_coupling_coreShell/nonuniform_scaling_trial7/case77_T/results/highRes_Spectrum_QM_E17_T=0_case77.dat", 'r') if line[0]!="#"])
y4 = np.array([float(line.strip().split()[1]) for line in \
    open("../scale_coupling_coreShell/nonuniform_scaling_trial7/case77_T/results/highRes_Spectrum_QM_E17_T=0_case77.dat", 'r') if line[0]!="#"])
axs[1].plot(x1*1000, y1, "-k", linewidth=1, label="Frohlich model")
axs[1].plot(x2*1000, y2, "-.r", linewidth=1, label="Short-range model")
axs[1].plot(x3*1000, y3, "g", linewidth=1, label="Utzat, CdSe/CdS/ZnS") # 4K
axs[1].set(xlim=(1975,2250), ylim=(0, 1), xlabel="Energy [meV]", ylabel="Abs", title="Single-NC line shape at low T (5K)")
# axs[1].legend(bbox_to_anchor=(0, 1), loc="lower left", prop={"size": 10})
axs[1].legend(loc="upper right", prop={"size": 10})
axins = inset_axes(axs[1], width="100%", height="100%", bbox_to_anchor=(.5, .3, .48, .45), bbox_transform=axs[1].transAxes)
axins.plot(x1*1000, y1, "-k", linewidth=1)
axins.plot(x4*1000, y4, ":", color="dimgray", linewidth=1, label="No broadening")
axins.tick_params(labelleft=False)
axins.get_xaxis().set_ticks([2031.021218130123, 2032.65683813])
axins.get_yaxis().set_ticks([0, 0.5, 1.0])
axins.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axins.set(xlim=(2030,2034.5))
axins.legend(loc="center right", prop={"size": 10})

'''
x1 = np.array([float(line.strip().split()[0]) for line in \
    open("../scale_coupling_coreShell/nonuniform_scaling_trial6/a0.6_b2.5_g0.479/results/gaussian150fs_Spectrum_QM_E17_T=5_aIndex6_bIndex3.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("../scale_coupling_coreShell/nonuniform_scaling_trial6/a0.6_b2.5_g0.479/results/gaussian150fs_Spectrum_QM_E17_T=5_aIndex6_bIndex3.dat", 'r') if line[0]!="#"])
x2 = np.array([float(line.strip().split()[0]) for line in \
    open("../experimental_digitized/single_NC/Bawendi_PRB76_2007_singleNC_shape.dat", 'r') if line[0]!="#"])
y2 = np.array([float(line.strip().split()[1]) for line in \
    open("../experimental_digitized/single_NC/Bawendi_PRB76_2007_singleNC_shape.dat", 'r') if line[0]!="#"])
axs[2].plot(x1*1000, y1, "-k", linewidth=1, label="Frohlich model, CdSe/CdS, 5K")
axs[2].plot(4417-x2, y2, "-.m", linewidth=1, label="Bawendi PRB76, CdSe/ZnS, 5K")
axs[2].set(xlim=(1900,2200), ylim=(0, 1), xlabel="Energy [meV]", ylabel="Abs")
axs[2].legend(loc="upper right", prop={"size": 10})

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
axs[2].legend(bbox_to_anchor=(0, 1), loc="lower left", prop={"size": 7})

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
axs[3].legend(bbox_to_anchor=(0, 1), loc="lower left", prop={"size": 7})
'''

fig.tight_layout()
fig.savefig("manu_plot_singleNC.pdf")
