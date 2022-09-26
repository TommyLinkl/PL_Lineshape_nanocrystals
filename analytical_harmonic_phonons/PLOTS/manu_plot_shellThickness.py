import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import FormatStrFormatter
from matplotlib import rc
matplotlib.rcParams['font.family'] = "serif"
matplotlib.rcParams['font.sans-serif'] = "Times New Roman"
matplotlib.rcParams['mathtext.fontset'] = "stix"
rc('axes', linewidth=0.4)
matplotlib.rcParams['xtick.major.width'] = 0.4
matplotlib.rcParams['ytick.major.width'] = 0.4
fig = plt.figure(figsize=(10, 6))
ax1 = plt.subplot2grid((2, 3), (0, 0), colspan=2)
ax2 = plt.subplot2grid((2, 3), (0, 2), colspan=1)
ax3 = plt.subplot2grid((2, 3), (1, 0), colspan=1)
ax4 = plt.subplot2grid((2, 3), (1, 1), colspan=1)
ax5 = plt.subplot2grid((2, 3), (1, 2), colspan=1)

x2 = np.array([float(line.strip().split()[0]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_2MLCdS/results/NEW_gaussian500fs_Spectrum_QM_E17_T=5.dat", 'r') if line[0]!="#"])
y2 = np.array([float(line.strip().split()[1]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_2MLCdS/results/NEW_gaussian500fs_Spectrum_QM_E17_T=5.dat", 'r') if line[0]!="#"])
x3 = np.array([float(line.strip().split()[0]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_gaussian500fs_Spectrum_QM_E17_T=5.dat", 'r') if line[0]!="#"])
y3 = np.array([float(line.strip().split()[1]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_3MLCdS/results/NEW_gaussian500fs_Spectrum_QM_E17_T=5.dat", 'r') if line[0]!="#"])
x4 = np.array([float(line.strip().split()[0]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_4MLCdS/results/NEW_gaussian500fs_Spectrum_QM_E17_T=5.dat", 'r') if line[0]!="#"])
y4 = np.array([float(line.strip().split()[1]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_4MLCdS/results/NEW_gaussian500fs_Spectrum_QM_E17_T=5.dat", 'r') if line[0]!="#"])
x5 = np.array([float(line.strip().split()[0]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_5MLCdS/results/NEW_gaussian500fs_Spectrum_QM_E17_T=5.dat", 'r') if line[0]!="#"])
y5 = np.array([float(line.strip().split()[1]) for line in \
    open("../CALCS/coreShell_3.0nmCdSe_5MLCdS/results/NEW_gaussian500fs_Spectrum_QM_E17_T=5.dat", 'r') if line[0]!="#"])
ax1.plot(x2*1000, y2, "-", linewidth=1, label="3nm core, 2ML shell")
ax1.plot(x3*1000, y3, "-", linewidth=1, label="3nm core, 3ML shell")
ax1.plot(x4*1000, y4, "-", linewidth=1, label="3nm core, 4ML shell")
ax1.plot(x5*1000, y5, "-", linewidth=1, label="3nm core, 5ML shell")
ax1.set(xlim=(1900,2150), ylim=(0, 1), xlabel="Energy [meV]", ylabel="Abs", title="Photoluminescence spectra at 5K")
ax1.legend(loc="upper right", prop={"size": 10})

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("3.0nm_shellComparison_FWHM.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[2]) for line in \
    open("3.0nm_shellComparison_FWHM.dat", 'r') if line[0]!="#"])
x2 = np.array([float(line.strip().split()[0]) for line in \
    open("3.9nm_shellComparison_FWHM.dat", 'r') if line[0]!="#"])
y2 = np.array([float(line.strip().split()[2]) for line in \
    open("3.9nm_shellComparison_FWHM.dat", 'r') if line[0]!="#"])
ax2.plot(x1, y1, "oc", markersize=5, label="3.0nm CdSe core")
ax2.plot(x2, y2, "*m", markersize=5, label="3.9nm CdSe core")
ax2.set(xlabel="Shell Layers [Monolayers]", ylabel="Reorganization Energy [meV]")
ax2.legend(loc="upper right", prop={"size": 10})

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("3.0nm_shellComparison_FWHM.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[1]) for line in \
    open("3.0nm_shellComparison_FWHM.dat", 'r') if line[0]!="#"])
x2 = np.array([float(line.strip().split()[0]) for line in \
    open("3.9nm_shellComparison_FWHM.dat", 'r') if line[0]!="#"])
y2 = np.array([float(line.strip().split()[1]) for line in \
    open("3.9nm_shellComparison_FWHM.dat", 'r') if line[0]!="#"])
ax3.plot(x1, y1*1000, "oc", markersize=5, label="3.0nm CdSe core")
ax3.plot(x2, y2*1000, "*m", markersize=5, label="3.9nm CdSe core")
ax3.set(ylim=(3,7), title="Width of the ZPL", xlabel="Shell Layers [Monolayers]", ylabel="FWHM [meV]")

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("3.0nm_shellComparison_FWHM.dat", 'r') if line[0]!="#"])
y1 = np.array([float(line.strip().split()[5]) for line in \
    open("3.0nm_shellComparison_FWHM.dat", 'r') if line[0]!="#"])
x2 = np.array([float(line.strip().split()[0]) for line in \
    open("3.9nm_shellComparison_FWHM.dat", 'r') if line[0]!="#"])
y2 = np.array([float(line.strip().split()[5]) for line in \
    open("3.9nm_shellComparison_FWHM.dat", 'r') if line[0]!="#"])
ax4.plot(x1, y1*1000, "oc", markersize=5, label="3.0nm CdSe core")
ax4.plot(x2, y2*1000, "*m", markersize=5, label="3.9nm CdSe core")
ax4.set(ylim=(3,7), title="Width of the phonon sideband", xlabel="Shell Layers [Monolayers]", ylabel="FWHM [meV]")

x1 = np.array([float(line.strip().split()[0]) for line in \
    open("3.0nm_shellComparison_FWHM.dat", 'r') if line[0]!="#"])
y1_main = np.array([float(line.strip().split()[3]) for line in \
    open("3.0nm_shellComparison_FWHM.dat", 'r') if line[0]!="#"])
y1_LO = np.array([float(line.strip().split()[4]) for line in \
    open("3.0nm_shellComparison_FWHM.dat", 'r') if line[0]!="#"])
x2 = np.array([float(line.strip().split()[0]) for line in \
    open("3.9nm_shellComparison_FWHM.dat", 'r') if line[0]!="#"])
y2_main = np.array([float(line.strip().split()[3]) for line in \
    open("3.9nm_shellComparison_FWHM.dat", 'r') if line[0]!="#"])
y2_LO = np.array([float(line.strip().split()[4]) for line in \
    open("3.9nm_shellComparison_FWHM.dat", 'r') if line[0]!="#"])
ax5.plot(x1, (y1_LO-y1_main)*1000, "oc", markersize=5, label="3.0nm CdSe core")
ax5.plot(x2, (y2_LO-y2_main)*1000, "*m", markersize=5, label="3.9nm CdSe core")
ax5.set(title="Sideband shift from main peak", xlabel="Shell Layers [Monolayers]", ylabel="Shift [meV]")

fig.tight_layout()
fig.savefig("manu_plot_shellThickness.pdf")
