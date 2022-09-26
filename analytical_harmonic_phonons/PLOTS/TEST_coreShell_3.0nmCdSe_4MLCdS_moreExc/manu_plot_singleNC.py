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
    open("../../CALCS/coreShell_3.0nmCdSe_4MLCdS/results/NEW_gaussian500fs_Spectrum_QM_E17_T=5.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("../../CALCS/coreShell_3.0nmCdSe_4MLCdS/results/NEW_gaussian500fs_Spectrum_QM_E17_T=5.dat", 'r') if line[0]!="#"])
x2 = np.array([float(line.strip().split()[0]) for line in \
    open("../../CALCS/coreShell_3.0nmCdSe_4MLCdS/results/OLD_gaussian500fs_Spectrum_QM_E17_T=5.dat", 'r') if line[0]!="#"])
y2 = np.array([float(line.strip().split()[1]) for line in \
    open("../../CALCS/coreShell_3.0nmCdSe_4MLCdS/results/OLD_gaussian500fs_Spectrum_QM_E17_T=5.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("../../experimental_digitized/single_NC/Hendrik_4K_shape.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("../../experimental_digitized/single_NC/Hendrik_4K_shape.dat", 'r') if line[0]!="#"])
x4 = np.array([float(line.strip().split()[0]) for line in \
    open("../../CALCS/coreShell_3.0nmCdSe_4MLCdS/results/NEW_Spectrum_QM_E17_T=5.dat", 'r') if line[0]!="#"])
y4 = np.array([float(line.strip().split()[1]) for line in \
    open("../../CALCS/coreShell_3.0nmCdSe_4MLCdS/results/NEW_Spectrum_QM_E17_T=5.dat", 'r') if line[0]!="#"])
axs[0].plot(x1*1000, y1, "-k", linewidth=1, label="NEW code, 500fs")
axs[0].plot(x2*1000, y2, "-.r", linewidth=1, label="OLD code, 500fs")
axs[0].plot(4060-x3, y3, "g", linewidth=1, label="Utzat, CdSe/CdS/ZnS") # 4K
axs[0].set(xlim=(1900,2100), ylim=(0, 1), xlabel="Energy [meV]", ylabel="Abs", title="Single-NC line shape at low T (5K)")
axs[0].legend(loc="upper right", prop={"size": 10})
axins = inset_axes(axs[0], width="100%", height="100%", bbox_to_anchor=(.5, .3, .48, .45), bbox_transform=axs[0].transAxes)
axins.plot(x1*1000, y1, "-k", linewidth=1)
axins.plot(x4*1000, y4, ":", color="dimgray", linewidth=1, label="No broadening")
axins.tick_params(labelleft=False)
axins.get_xaxis().set_ticks([1940.6130346163348, 1942.07960462])
axins.get_yaxis().set_ticks([0, 0.5, 1.0])
axins.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axins.set(xlim=(1939.5,1944))
axins.legend(loc="center right", prop={"size": 10})


x1 = np.array([float(line.strip().split()[0]) for line in \
    open("../../CALCS/coreShell_3.0nmCdSe_4MLCdS/results/NEW_gaussian30fs_Spectrum_QM_E17_T=300.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("../../CALCS/coreShell_3.0nmCdSe_4MLCdS/results/NEW_gaussian30fs_Spectrum_QM_E17_T=300.dat", 'r') if line[0]!="#"])
x2 = np.array([float(line.strip().split()[0]) for line in \
    open("../../CALCS/coreShell_3.0nmCdSe_4MLCdS/results/OLD_gaussian30fs_Spectrum_QM_E17_T=300.dat", 'r') if line[0]!="#"])
y2 = np.array([float(line.strip().split()[1]) for line in \
    open("../../CALCS/coreShell_3.0nmCdSe_4MLCdS/results/OLD_gaussian30fs_Spectrum_QM_E17_T=300.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("../../experimental_digitized/single_NC/Hendrik_290K_shape.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("../../experimental_digitized/single_NC/Hendrik_290K_shape.dat", 'r') if line[0]!="#"])
axs[1].plot(x1*1000, y1, "-k", linewidth=1, label="NEW code, 30fs")
axs[1].plot(x2*1000, y2, "-.r", linewidth=1, label="OLD code, 30fs")
axs[1].plot(3950-x3, y3, "g", linewidth=1, label="Utzat, CdSe/CdS/ZnS") # 290K
axs[1].set(xlim=(1800,2100), ylim=(0, 1), xlabel="Energy [meV]", ylabel="Abs", title="Single-NC line shape at high T (300K)")
axs[1].legend(loc="upper right", prop={"size": 10})



fig.tight_layout()
fig.savefig("manu_plot_singleNC.pdf")
