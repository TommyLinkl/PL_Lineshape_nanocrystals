import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import FormatStrFormatter
from matplotlib import rc
matplotlib.rcParams['font.family'] = "serif"
matplotlib.rcParams['font.sans-serif'] = "Times New Roman"
matplotlib.rcParams['font.size'] = 12
matplotlib.rcParams['mathtext.fontset'] = "stix"
rc('axes', linewidth=0.4)
matplotlib.rcParams['xtick.major.width'] = 0.4
matplotlib.rcParams['ytick.major.width'] = 0.4

fig = plt.figure(figsize=(10, 13))
ax1 = plt.subplot2grid((4, 8), (0, 0), colspan=5)
ax2 = plt.subplot2grid((4, 8), (1, 0), colspan=4)
ax3 = plt.subplot2grid((4, 8), (1, 4), colspan=4)
ax4 = plt.subplot2grid((4, 8), (2, 0), colspan=4)
ax5 = plt.subplot2grid((4, 8), (2, 4), colspan=4)
ax6 = plt.subplot2grid((4, 8), (3, 0), colspan=4)
ax7 = plt.subplot2grid((4, 8), (3, 4), colspan=4)
ax8 = plt.subplot2grid((4, 8), (0, 5), colspan=3)

# fig.suptitle("Single-NC line shape")

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_gaussian500fs_Spectrum_QM_E17_T=4.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_gaussian500fs_Spectrum_QM_E17_T=4.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("../experimental_digitized/single_NC/Hendrik_4K_shape.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("../experimental_digitized/single_NC/Hendrik_4K_shape.dat", 'r') if line[0]!="#"])
x4 = np.array([float(line.strip().split()[0]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_Spectrum_QM_E17_T=4.dat", 'r') if line[0]!="#"])
y4 = np.array([float(line.strip().split()[1]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_Spectrum_QM_E17_T=4.dat", 'r') if line[0]!="#"])
ax1.plot(x1*1000, y1, "-k", linewidth=1, label="Model")   # , Gauss500fs, $T_2=\infty$
ax1.plot(x4*1000, y4*np.sign(y4), ":", color="dimgray", linewidth=1, label="Model w/o broadening")
ax1.plot(4108-x3, y3, "g", linewidth=1, label="Single-QD PL Expt") 
ax1.set(xlim=(1950,2150), ylim=(0, 1), xlabel="Energy [meV]", ylabel="Abs")
ax1.legend(loc="upper right", prop={"size": 10})
ax1.text(1955, 0.9, "4K")
ax1.get_yaxis().set_ticks([0.0, 0.5, 1.0])
axins = inset_axes(ax1, width="100%", height="100%", bbox_to_anchor=(.5, .24, .43, .40), bbox_transform=ax1.transAxes)
axins.plot(x1*1000, y1, "-k", linewidth=1)
axins.plot(x4*1000, y4*np.sign(y4), ":", color="dimgray", linewidth=1, label="No broadening")
axins.tick_params(labelleft=False)
axins.get_xaxis().set_ticks([1989.0195080144384, 1995])
axins.get_yaxis().set_ticks([0, 0.5, 1.0])
axins.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axins.set(xlim=(1988,1998))

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_gaussian500fs_Spectrum_QM_E17_T=30.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_gaussian500fs_Spectrum_QM_E17_T=30.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("../experimental_digitized/single_NC/Hendrik_30K_shape.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("../experimental_digitized/single_NC/Hendrik_30K_shape.dat", 'r') if line[0]!="#"])
ax2.plot(x1*1000, y1, "-k", linewidth=1, label=r"Model, $T_2=\infty$")
ax2.plot(4092-x3, y3, "g", linewidth=1, label="Expt")
ax2.set(xlim=(1950,2150), ylim=(0, 1))
ax2.text(2125, 0.9, "30K")
ax2.get_yaxis().set_ticks([0.0, 0.5, 1.0])
# ax2.legend(loc="upper right", prop={"size": 10})

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_gaussian500fs_Spectrum_QM_E17_T=60.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_gaussian500fs_Spectrum_QM_E17_T=60.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("../experimental_digitized/single_NC/Hendrik_60K_shape.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("../experimental_digitized/single_NC/Hendrik_60K_shape.dat", 'r') if line[0]!="#"])
ax3.plot(x1*1000, y1, "-k", linewidth=1, label=r"Model, $T_2=\infty$")
ax3.plot(4112-x3, y3, "g", linewidth=1, label="Expt")
ax3.set(xlim=(1950,2150), ylim=(0, 1))
ax3.text(2125, 0.9, "60K")
ax3.get_yaxis().set_ticks([0.0, 0.5, 1.0])
# ax3.legend(loc="upper right", prop={"size": 10})

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_gaussian500fs_Spectrum_QM_E17_T=100.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_gaussian500fs_Spectrum_QM_E17_T=100.dat", 'r') if line[0]!="#"])
x2 = np.array([float(line.strip().split()[0]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_exponential1000fs_gaussian500fs_Spectrum_QM_E17_T=100.dat", 'r') if line[0]!="#"])
y2 = np.array([float(line.strip().split()[1]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_exponential1000fs_gaussian500fs_Spectrum_QM_E17_T=100.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("../experimental_digitized/single_NC/Hendrik_100K_shape.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("../experimental_digitized/single_NC/Hendrik_100K_shape.dat", 'r') if line[0]!="#"])
# ax4.plot(x1*1000, y1, "-k", linewidth=1, label="Same Gauss, T2=Inf")
ax4.plot(x2*1000, y2, "-k", linewidth=1, label=r"Model, $T_2=1000fs$")
ax4.plot(4064-x3, y3, "g", linewidth=1, label="Expt")
ax4.set(xlim=(1950,2150), ylim=(0, 1))
ax4.text(2120, 0.9, "100K")
ax4.get_yaxis().set_ticks([0.0, 0.5, 1.0])
# ax4.legend(loc="upper right", prop={"size": 10})

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_gaussian500fs_Spectrum_QM_E17_T=150.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_gaussian500fs_Spectrum_QM_E17_T=150.dat", 'r') if line[0]!="#"])
x2 = np.array([float(line.strip().split()[0]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_exponential500fs_gaussian500fs_Spectrum_QM_E17_T=150.dat", 'r') if line[0]!="#"])
y2 = np.array([float(line.strip().split()[1]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_exponential500fs_gaussian500fs_Spectrum_QM_E17_T=150.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("../experimental_digitized/single_NC/Hendrik_150K_shape.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("../experimental_digitized/single_NC/Hendrik_150K_shape.dat", 'r') if line[0]!="#"])
# ax5.plot(x1*1000, y1, "-k", linewidth=1, label="Same Gauss, T2=Inf")
ax5.plot(x2*1000, y2, "-k", linewidth=1, label=r"Model, $T_2=500fs$")
ax5.plot(4076-x3, y3, "g", linewidth=1, label="Expt")
ax5.set(xlim=(1950,2150), ylim=(0, 1))
ax5.text(2120, 0.9, "150K")
ax5.get_yaxis().set_ticks([0.0, 0.5, 1.0])
# ax5.legend(loc="upper right", prop={"size": 10})

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_gaussian500fs_Spectrum_QM_E17_T=220.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_gaussian500fs_Spectrum_QM_E17_T=220.dat", 'r') if line[0]!="#"])
x2 = np.array([float(line.strip().split()[0]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_exponential50fs_gaussian500fs_Spectrum_QM_E17_T=220.dat", 'r') if line[0]!="#"])
y2 = np.array([float(line.strip().split()[1]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_exponential50fs_gaussian500fs_Spectrum_QM_E17_T=220.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("../experimental_digitized/single_NC/Hendrik_220K_shape.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("../experimental_digitized/single_NC/Hendrik_220K_shape.dat", 'r') if line[0]!="#"])
# ax6.plot(x1*1000, y1, "-k", linewidth=1, label="Same Gauss, T2=Inf")
ax6.plot(x2*1000, y2, "-k", linewidth=1, label=r"Model, $T_2=50fs$")
ax6.plot(4030-x3, y3, "g", linewidth=1, label="Expt")
ax6.set(xlim=(1900,2150), ylim=(0, 1))
ax6.text(2115, 0.9, "220K")
ax6.get_yaxis().set_ticks([0.0, 0.5, 1.0])
# ax6.legend(loc="upper right", prop={"size": 10})

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_gaussian500fs_Spectrum_QM_E17_T=290.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_gaussian500fs_Spectrum_QM_E17_T=290.dat", 'r') if line[0]!="#"])
x2 = np.array([float(line.strip().split()[0]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_exponential30fs_gaussian500fs_Spectrum_QM_E17_T=290.dat", 'r') if line[0]!="#"])
y2 = np.array([float(line.strip().split()[1]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_exponential30fs_gaussian500fs_Spectrum_QM_E17_T=290.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("../experimental_digitized/single_NC/Hendrik_290K_shape.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("../experimental_digitized/single_NC/Hendrik_290K_shape.dat", 'r') if line[0]!="#"])
# ax7.plot(x1*1000, y1, "-k", linewidth=1, label="Same Gauss, T2=Inf")
ax7.plot(x2*1000, y2, "-k", linewidth=1, label=r"Model, $T_2=30fs$")
ax7.plot(4000-x3, y3, "g", linewidth=1, label="Expt")
ax7.set(xlim=(1900,2150), ylim=(0, 1))
ax7.text(2115, 0.9, "290K")
ax7.get_yaxis().set_ticks([0.0, 0.5, 1.0])
# ax7.legend(loc="upper right", prop={"size": 10})

Temp_range = np.array([4, 30, 60, 100, 150, 220, 290])
oneOverT2_range = np.array([0, 0, 0, 1/1000, 1/500, 1/50, 1/30])
ax8.plot(Temp_range, oneOverT2_range, "-ok", linewidth=1, markersize=5)
ax8.set(xlabel="Temperature [K]", ylabel=r"Dephasing rate $1/T_2$ [1/fs]")  #, title="Best fit dephasing rate"
ax8.set_ylim(bottom=0, top=None)
# ax8.legend(loc="upper right", prop={"size": 10})
ax8.get_yaxis().set_ticks([0.0, 0.01, 0.02, 0.03])

fig.tight_layout()
fig.savefig("manu_plot_singleNC_allT.pdf")
