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
fig, axs = plt.subplots(4,2, figsize=(8, 12))
fig.suptitle("Single-NC line shape, gaussian 50fs")

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("../../CALCS/coreShell_3.0nmCdSe_2MLCdS/results/NEW_gaussian50fs_Spectrum_QM_E17_T=4.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("../../CALCS/coreShell_3.0nmCdSe_2MLCdS/results/NEW_gaussian50fs_Spectrum_QM_E17_T=4.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("../../experimental_digitized/single_NC/Hendrik_4K_shape.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("../../experimental_digitized/single_NC/Hendrik_4K_shape.dat", 'r') if line[0]!="#"])
axs[0,0].plot(x1*1000, y1, "-k", linewidth=1, label="Model")
axs[0,0].plot(4162-x3, y3, "g", linewidth=1, label="Utzat") 
axs[0,0].set(xlim=(2000,2200), ylim=(0, 1), xlabel="Energy [meV]", ylabel="Abs", title="4K")
axs[0,0].legend(loc="upper right", prop={"size": 10})

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("../../CALCS/coreShell_3.0nmCdSe_2MLCdS/results/NEW_gaussian50fs_Spectrum_QM_E17_T=30.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("../../CALCS/coreShell_3.0nmCdSe_2MLCdS/results/NEW_gaussian50fs_Spectrum_QM_E17_T=30.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("../../experimental_digitized/single_NC/Hendrik_30K_shape.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("../../experimental_digitized/single_NC/Hendrik_30K_shape.dat", 'r') if line[0]!="#"])
axs[0,1].plot(x1*1000, y1, "-k", linewidth=1, label="Model")
axs[0,1].plot(4146-x3, y3, "g", linewidth=1, label="Utzat")
axs[0,1].set(xlim=(2000,2200), ylim=(0, 1), xlabel="Energy [meV]", ylabel="Abs", title="30K")

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("../../CALCS/coreShell_3.0nmCdSe_2MLCdS/results/NEW_gaussian50fs_Spectrum_QM_E17_T=60.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("../../CALCS/coreShell_3.0nmCdSe_2MLCdS/results/NEW_gaussian50fs_Spectrum_QM_E17_T=60.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("../../experimental_digitized/single_NC/Hendrik_60K_shape.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("../../experimental_digitized/single_NC/Hendrik_60K_shape.dat", 'r') if line[0]!="#"])
axs[1,0].plot(x1*1000, y1, "-k", linewidth=1, label="Model")
axs[1,0].plot(4166-x3, y3, "g", linewidth=1, label="Utzat")
axs[1,0].set(xlim=(2000,2200), ylim=(0, 1), xlabel="Energy [meV]", ylabel="Abs", title="60K")

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("../../CALCS/coreShell_3.0nmCdSe_2MLCdS/results/NEW_gaussian50fs_Spectrum_QM_E17_T=100.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("../../CALCS/coreShell_3.0nmCdSe_2MLCdS/results/NEW_gaussian50fs_Spectrum_QM_E17_T=100.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("../../experimental_digitized/single_NC/Hendrik_100K_shape.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("../../experimental_digitized/single_NC/Hendrik_100K_shape.dat", 'r') if line[0]!="#"])
axs[1,1].plot(x1*1000, y1, "-k", linewidth=1, label="Model")
axs[1,1].plot(4120-x3, y3, "g", linewidth=1, label="Utzat")
axs[1,1].set(xlim=(2000,2200), ylim=(0, 1), xlabel="Energy [meV]", ylabel="Abs", title="100K")

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("../../CALCS/coreShell_3.0nmCdSe_2MLCdS/results/NEW_gaussian50fs_Spectrum_QM_E17_T=150.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("../../CALCS/coreShell_3.0nmCdSe_2MLCdS/results/NEW_gaussian50fs_Spectrum_QM_E17_T=150.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("../../experimental_digitized/single_NC/Hendrik_150K_shape.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("../../experimental_digitized/single_NC/Hendrik_150K_shape.dat", 'r') if line[0]!="#"])
axs[2,0].plot(x1*1000, y1, "-k", linewidth=1, label="Model")
axs[2,0].plot(4132-x3, y3, "g", linewidth=1, label="Utzat")
axs[2,0].set(xlim=(2000,2200), ylim=(0, 1), xlabel="Energy [meV]", ylabel="Abs", title="150K")

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("../../CALCS/coreShell_3.0nmCdSe_2MLCdS/results/NEW_gaussian50fs_Spectrum_QM_E17_T=220.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("../../CALCS/coreShell_3.0nmCdSe_2MLCdS/results/NEW_gaussian50fs_Spectrum_QM_E17_T=220.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("../../experimental_digitized/single_NC/Hendrik_220K_shape.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("../../experimental_digitized/single_NC/Hendrik_220K_shape.dat", 'r') if line[0]!="#"])
axs[2,1].plot(x1*1000, y1, "-k", linewidth=1, label="Model")
axs[2,1].plot(4087-x3, y3, "g", linewidth=1, label="Utzat")
axs[2,1].set(xlim=(2000,2200), ylim=(0, 1), xlabel="Energy [meV]", ylabel="Abs", title="220K")

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("../../CALCS/coreShell_3.0nmCdSe_2MLCdS/results/NEW_gaussian50fs_Spectrum_QM_E17_T=290.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("../../CALCS/coreShell_3.0nmCdSe_2MLCdS/results/NEW_gaussian50fs_Spectrum_QM_E17_T=290.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("../../experimental_digitized/single_NC/Hendrik_290K_shape.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("../../experimental_digitized/single_NC/Hendrik_290K_shape.dat", 'r') if line[0]!="#"])
axs[3,0].plot(x1*1000, y1, "-k", linewidth=1, label="Model")
axs[3,0].plot(4060-x3, y3, "g", linewidth=1, label="Utzat")
axs[3,0].set(xlim=(2000,2200), ylim=(0, 1), xlabel="Energy [meV]", ylabel="Abs", title="290K")


fig.tight_layout()
fig.savefig("manu_plot_singleNC_allT_50fs.pdf")
