import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib import rc
from allFunc import * 
matplotlib.rcParams['font.family'] = "serif"
matplotlib.rcParams['font.sans-serif'] = "Times New Roman"
matplotlib.rcParams['mathtext.fontset'] = "stix"
rc('axes', linewidth=0.4)
matplotlib.rcParams['xtick.major.width'] = 0.4
matplotlib.rcParams['ytick.major.width'] = 0.4
fig, axs = plt.subplots(1,3, figsize=(14, 4))

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
axs[0].plot(x1*1000, y1, "-k", linewidth=1, label="This work")
axs[0].plot(x4*1000, y4, ":", color="dimgray", linewidth=1, label="This work, no broadening")
axs[0].plot(4108-x3, y3, "g", linewidth=1, label="Utzat, CdSe/CdS/ZnS, shifted") # 4K
axs[0].set(xlim=(1925,2150), ylim=(0, 1), xlabel="Energy [meV]", ylabel="Abs", title="Single-NC line shape at 4K")
axs[0].legend(loc="upper right", prop={"size": 10})
axins = inset_axes(axs[0], width="100%", height="100%", bbox_to_anchor=(.5, .3, .45, .40), bbox_transform=axs[0].transAxes)
axins.plot(x1*1000, y1, "-k", linewidth=1)
axins.plot(x4*1000, y4, ":", color="dimgray", linewidth=1, label="No broadening")
axins.tick_params(labelleft=False)
axins.get_xaxis().set_ticks([1989.0195080144384, 1995])
axins.get_yaxis().set_ticks([0, 0.5, 1.0])
axins.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axins.set(xlim=(1988,1998))
# axins.legend(loc="center right", prop={"size": 10})


x1 = np.array([float(line.strip().split()[0]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_gaussian500fs_Spectrum_QM_E17_T=60.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_gaussian500fs_Spectrum_QM_E17_T=60.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("../experimental_digitized/single_NC/Hendrik_60K_shape.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("../experimental_digitized/single_NC/Hendrik_60K_shape.dat", 'r') if line[0]!="#"])
axs[1].plot(x1*1000, y1, "-k", linewidth=1, label="This work")
axs[1].plot(4112-x3, y3, "g", linewidth=1, label="Utzat, CdSe/CdS/ZnS")
axs[1].set(xlim=(1925,2150), ylim=(0, 1), xlabel="Energy [meV]", title="60K")
# axs[1].legend(loc="upper right", prop={"size": 10})

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_gaussian500fs_Spectrum_QM_E17_T=150.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_gaussian500fs_Spectrum_QM_E17_T=150.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("../experimental_digitized/single_NC/Hendrik_150K_shape.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("../experimental_digitized/single_NC/Hendrik_150K_shape.dat", 'r') if line[0]!="#"])
axs[2].plot(x1*1000, y1, "-k", linewidth=1, label="This work")
axs[2].plot(4077-x3, y3, "g", linewidth=1, label="Utzat, CdSe/CdS/ZnS")
axs[2].set(xlim=(1925,2150), ylim=(0, 1), xlabel="Energy [meV]", title="150K")
# axs[2].legend(loc="upper right", prop={"size": 10})

fig.tight_layout()
fig.savefig("manu_plot_singleNC_lowT.pdf")
