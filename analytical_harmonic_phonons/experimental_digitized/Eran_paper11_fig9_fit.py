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

beta_au = np.array([float(line.strip().split()[0]) for line in open("Eran_paper11_fig9.dat", 'r') if line[0]!="#"])
oneOverT2_au = np.array([float(line.strip().split()[1]) for line in open("Eran_paper11_fig9.dat", 'r') if line[0]!="#"])

T_SI = 1/(KB*beta_au/AUTOJ)
oneOverT2_SI = oneOverT2_au / AUTOS

fig, axs = plt.subplots(1,1, figsize=(5, 5))
axs.plot(T_SI, oneOverT2_SI*1e15)
axs.set(xlabel="Temperature [K]", ylabel="1/T2 [1/fs]")

fig.tight_layout()
fig.savefig("Eran_paper11_fig9_plot.pdf")
