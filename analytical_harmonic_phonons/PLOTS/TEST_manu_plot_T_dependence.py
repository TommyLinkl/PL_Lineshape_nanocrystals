import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
matplotlib.rcParams['font.family'] = "serif"
matplotlib.rcParams['font.sans-serif'] = "Times New Roman"
rc('axes', linewidth=0.4)
matplotlib.rcParams['xtick.major.width'] = 0.4
matplotlib.rcParams['ytick.major.width'] = 0.4

KB = 1.380649e-23                # unit: J/K, SI units
PI = 3.14159265359
HBAR = 1.054571817e-34           # Unit: Js, SI units
JTOEV = 6.24150907446076e18 
AUTOEV = 27.211396641308 
AUTOJ = 4.359744e-18
AUTOS = 2.4188843265857e-17      # a.u. of time to seconds

fig, axs = plt.subplots(1,3, figsize=(12, 6))

x1 = np.array([float(line.strip().split()[0]) for line in open("../CALCS/coreShell_3.0nmCdSe_2MLCdS/results/NEW_WidthFreqShift_vs_T_gaussian_40fs.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in open("../CALCS/coreShell_3.0nmCdSe_2MLCdS/results/NEW_WidthFreqShift_vs_T_gaussian_40fs.dat", 'r') if line[0]!="#"])
x2 = np.array([float(line.strip().split()[0]) for line in open("../CALCS/coreShell_3.0nmCdSe_2MLCdS/results/NEW_WidthFreqShift_vs_T_exponential50fs_gaussian30fs.dat", 'r') if line[0]!="#"])
y2 = np.array([float(line.strip().split()[1]) for line in open("../CALCS/coreShell_3.0nmCdSe_2MLCdS/results/NEW_WidthFreqShift_vs_T_exponential50fs_gaussian30fs.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in open("../experimental_digitized/ensemble/Meijerink_coreShell_inhomo_width.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in open("../experimental_digitized/ensemble/Meijerink_coreShell_inhomo_width.dat", 'r') if line[0]!="#"])
axs[0].plot(x1, y1*1000, "-+k", linewidth=1, label="Gauss40fs (T2=Inf)")
axs[0].plot(100, 0.051406*1000, "-<r", label="Same Gauss, T2=1000fs")
axs[0].plot(150, 0.055418*1000, "-sr", label="Same Gauss, T2=500fs")
axs[0].plot(220, 0.076055*1000, "-or", label="Same Gauss, T2=50fs")
axs[0].plot(290, 0.085402*1000, "-^r", label="Same Gauss, T2=30fs")
axs[0].plot(x3, y3, ":sb", linewidth=1, label="Meijerink JPCC")
axs[0].set(ylim=(0,110), xlabel="Temperature [K]", ylabel="FWHM [meV]", title="3nm/2ML Inhomogeneous PL Linewidth")
axs[0].legend(loc="lower right", prop={"size": 10})

x1 = np.array([float(line.strip().split()[0]) for line in open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_WidthFreqShift_vs_T_gaussian_500fs.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_WidthFreqShift_vs_T_gaussian_500fs.dat", 'r') if line[0]!="#"])
x2 = np.array([float(line.strip().split()[0]) for line in open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_WidthFreqShift_vs_T_exponential50fs_gaussian500fs.dat", 'r') if line[0]!="#"])
y2 = np.array([float(line.strip().split()[1]) for line in open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_WidthFreqShift_vs_T_exponential50fs_gaussian500fs.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in open("../experimental_digitized/single_NC/Hendrik_Dot1_width.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in open("../experimental_digitized/single_NC/Hendrik_Dot1_width.dat", 'r') if line[0]!="#"])
axs[1].plot(x1, y1*1000, "-+k", linewidth=1, label="Gauss500fs (T2=Inf)")
axs[1].plot(100, 0.015343*1000, "-<r", label="Same Gauss, T2=1000fs")
axs[1].plot(150, 0.019934*1000, "-sr", label="Same Gauss, T2=500fs")
axs[1].plot(220, 0.049711*1000, "-or", label="Same Gauss, T2=50fs")
axs[1].plot(290, 0.071754*1000, "-^r", label="Same Gauss, T2=30fs")
axs[1].plot(x3, y3, ":og", linewidth=1, label="Utzat, CdSe/CdS/ZnS")
axs[1].set(ylim=(0,110), xlabel="Temperature [K]", ylabel="FWHM [meV]", title="3nm/3ML Single-NC PL Linewidth")
axs[1].legend(loc="upper right", prop={"size": 10})

Temp_range = np.array([100, 150, 220, 290])
T2_range = np.array([1000, 500, 50, 30])
axs[2].plot(Temp_range, 1/T2_range, "-ok", label="Best fit, this work")
axs[2].set(ylim=(0, 0.04), xlabel="Temperature [K]", ylabel=r"$1/T_2$ [1/fs]", title="Best fit population relaxation rate")
'''
axs[2].plot(1/(Temp_range*KB/AUTOJ), 1/(T2_range*1e-15/AUTOS), "-ok", label="Best fit, this work")
axs[2].set(xlabel=r"$\beta [au]", ylabel=r"$1/T_2$ [au]", title="Best fit population relaxation rate")
'''
axs[2].legend(loc="upper right", prop={"size": 10})

fig.tight_layout()
fig.savefig("TEST_manu_plot_T_dependence.pdf")
