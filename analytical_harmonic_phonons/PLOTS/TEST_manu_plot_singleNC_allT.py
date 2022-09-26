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
fig.suptitle("Single-NC line shape, gaussian 500fs")

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_gaussian500fs_Spectrum_QM_E17_T=4.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_gaussian500fs_Spectrum_QM_E17_T=4.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("../experimental_digitized/single_NC/Hendrik_4K_shape.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("../experimental_digitized/single_NC/Hendrik_4K_shape.dat", 'r') if line[0]!="#"])
axs[0,0].plot(x1*1000, y1, "-k", linewidth=1, label="Gauss500fs, T2=Inf")
axs[0,0].plot(4108-x3, y3, "g", linewidth=1, label="Utzat") 
axs[0,0].set(xlim=(1950,2150), ylim=(0, 1), xlabel="Energy [meV]", ylabel="Abs", title="4K")
axs[0,0].legend(loc="upper right", prop={"size": 10})

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_gaussian500fs_Spectrum_QM_E17_T=30.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_gaussian500fs_Spectrum_QM_E17_T=30.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("../experimental_digitized/single_NC/Hendrik_30K_shape.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("../experimental_digitized/single_NC/Hendrik_30K_shape.dat", 'r') if line[0]!="#"])
axs[0,1].plot(x1*1000, y1, "-k", linewidth=1, label="Same Gauss, T2=Inf")
axs[0,1].plot(4092-x3, y3, "g", linewidth=1, label="Utzat")
axs[0,1].set(xlim=(1950,2150), ylim=(0, 1), xlabel="Energy [meV]", ylabel="Abs", title="30K")
axs[0,1].legend(loc="upper right", prop={"size": 10})

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_gaussian500fs_Spectrum_QM_E17_T=60.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_gaussian500fs_Spectrum_QM_E17_T=60.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("../experimental_digitized/single_NC/Hendrik_60K_shape.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("../experimental_digitized/single_NC/Hendrik_60K_shape.dat", 'r') if line[0]!="#"])
axs[1,0].plot(x1*1000, y1, "-k", linewidth=1, label="Same Gauss, T2=Inf")
axs[1,0].plot(4112-x3, y3, "g", linewidth=1, label="Utzat")
axs[1,0].set(xlim=(1950,2150), ylim=(0, 1), xlabel="Energy [meV]", ylabel="Abs", title="60K")
axs[1,0].legend(loc="upper right", prop={"size": 10})

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
axs[1,1].plot(x1*1000, y1, "-k", linewidth=1, label="Same Gauss, T2=Inf")
axs[1,1].plot(x2*1000, y2, "-b", linewidth=1, label="Same Gauss, T2=1000fs")
axs[1,1].plot(4064-x3, y3, "g", linewidth=1, label="Utzat")
axs[1,1].set(xlim=(1950,2150), ylim=(0, 1), xlabel="Energy [meV]", ylabel="Abs", title="100K")
axs[1,1].legend(loc="upper right", prop={"size": 10})

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
axs[2,0].plot(x1*1000, y1, "-k", linewidth=1, label="Same Gauss, T2=Inf")
axs[2,0].plot(x2*1000, y2, "-b", linewidth=1, label="Same Gauss, T2=500fs")
axs[2,0].plot(4076-x3, y3, "g", linewidth=1, label="Utzat")
axs[2,0].set(xlim=(1950,2150), ylim=(0, 1), xlabel="Energy [meV]", ylabel="Abs", title="150K")
axs[2,0].legend(loc="upper right", prop={"size": 10})

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
axs[2,1].plot(x1*1000, y1, "-k", linewidth=1, label="Same Gauss, T2=Inf")
axs[2,1].plot(x2*1000, y2, "-b", linewidth=1, label="Same Gauss, T2=50fs")
axs[2,1].plot(4030-x3, y3, "g", linewidth=1, label="Utzat")
axs[2,1].set(xlim=(1900,2150), ylim=(0, 1), xlabel="Energy [meV]", ylabel="Abs", title="220K")
axs[2,1].legend(loc="upper right", prop={"size": 10})

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
axs[3,0].plot(x1*1000, y1, "-k", linewidth=1, label="Same Gauss, T2=Inf")
axs[3,0].plot(x2*1000, y2, "-b", linewidth=1, label="Same Gauss, T2=30fs")
axs[3,0].plot(4000-x3, y3, "g", linewidth=1, label="Utzat")
axs[3,0].set(xlim=(1900,2150), ylim=(0, 1), xlabel="Energy [meV]", ylabel="Abs", title="290K")
axs[3,0].legend(loc="upper right", prop={"size": 10})


fig.tight_layout()
fig.savefig("TEST_manu_plot_singleNC_allT.pdf")
