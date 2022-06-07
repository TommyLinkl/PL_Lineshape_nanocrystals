import numpy as np
import matplotlib.pyplot as plt

Hendrik_T = np.array([float(line.strip().split()[0]) for line in open("HendrikDot1_T_Width.dat", 'r')])
Hendrik_Width = np.array([float(line.strip().split()[1]) for line in open("HendrikDot1_T_Width.dat", 'r')])

Scholes_T = np.array([float(line.strip().split()[0]) for line in open("Scholes_T_Width.dat", 'r')])
Scholes_Width = np.array([float(line.strip().split()[1]) for line in open("Scholes_T_Width.dat", 'r')])

Meijerink_510_T = np.array([float(line.strip().split()[0]) for line in open("Meijerink_510QD.dat", 'r')])
Meijerink_510_Width = np.array([float(line.strip().split()[1]) for line in open("Meijerink_510QD.dat", 'r')])

Meijerink_580_T = np.array([float(line.strip().split()[0]) for line in open("Meijerink_580QD.dat", 'r')])
Meijerink_580_Width = np.array([float(line.strip().split()[1]) for line in open("Meijerink_580QD.dat", 'r')])

Meijerink_610_T = np.array([float(line.strip().split()[0]) for line in open("Meijerink_610QD.dat", 'r')])
Meijerink_610_Width = np.array([float(line.strip().split()[1]) for line in open("Meijerink_610QD.dat", 'r')])

Meijerink_coreShell_T = np.array([float(line.strip().split()[0]) for line in open("Meijerink_core-shellQD.dat", 'r')])
Meijerink_coreShell_Width = np.array([float(line.strip().split()[1]) for line in open("Meijerink_core-shellQD.dat", 'r')])

Bawendi_45A_T = np.array([float(line.strip().split()[0]) for line in open("Bawendi_45A.dat", 'r')])
Bawendi_45A_Width = np.array([float(line.strip().split()[1]) for line in open("Bawendi_45A.dat", 'r')])

Bawendi_coreShell_T = np.array([float(line.strip().split()[0]) for line in open("Bawendi_43ACoreShell.dat", 'r')])
Bawendi_coreShell_Width = np.array([float(line.strip().split()[1]) for line in open("Bawendi_43ACoreShell.dat", 'r')])

Lin_3nm_T = np.array([float(line.strip().split()[0]) for line in open("Lin_3nm.dat", 'r')])
Lin_3nm_Width = np.array([float(line.strip().split()[1]) for line in open("Lin_3nm.dat", 'r')])*1000

Lin_coreShell_T = np.array([float(line.strip().split()[0]) for line in open("Lin_3nm_3ML.dat", 'r')])
Lin_coreShell_Width = np.array([float(line.strip().split()[1]) for line in open("Lin_3nm_3ML.dat", 'r')])*1000

fig, axs = plt.subplots(1,2, figsize=(16,8))
fig.suptitle("Experimental and atomic pseudopotential-based PL results")
axs[0].grid()
axs[1].grid()

axs[0].set(xlabel="Temperature [K]", ylabel="Homogeneous Linewidth [meV]", title="Bare core QDs", xlim=(-10,450), ylim=(-15,160))
axs[0].plot(Lin_3nm_T, Lin_3nm_Width, "-", label="Lin, CdSe3nm, NEW a5 pp")
# axs[0].plot(Lin_3nm_T, Lin_3nm_Width-Lin_3nm_Width[0], "-", label="Lin, shifted, CdSe3nm")
axs[0].plot(Meijerink_510_T, Meijerink_510_Width, "s:", label="Meijerink, ensemble meas. deconvoluted, CdSe emitting at 510nm (?2.5nm)")
axs[0].plot(Meijerink_580_T, Meijerink_580_Width, "s:", label="Meijerink, ensemble meas. deconvoluted, CdSe3.3-3.9nm")
# axs[0].plot(Meijerink_610_T, Meijerink_610_Width, "s:", label="Meijerink, ensemble meas. deconvoluted, CdSe4.1-4.7nm")
axs[0].plot(Scholes_T, Scholes_Width, "<:", label="Scholes, 3PEPS, ?Maybe simulated? Not sure what QD")
axs[0].plot(Bawendi_45A_T, Bawendi_45A_Width, "d:", label="Bawendi, SDL, CdSe4.5nm")

axs[1].set(xlabel="Temperature [K]", ylabel="Homogeneous Linewidth [meV]", title="Core-shell QDs", xlim=(-10,450), ylim=(-15,160))
axs[1].plot(Lin_coreShell_T, Lin_coreShell_Width, "-", label="Lin, CdSe3nm/3MLCdS, NEW a5 pp")
# axs[1].plot(Lin_coreShell_T, Lin_coreShell_Width-Lin_coreShell_Width[0], "-", label="Lin, shifted, CdSe3nm/3?MLCdS")
axs[1].plot(Hendrik_T, Hendrik_Width, "o:", label="Hendrik manu., \'Dot1\', CdSe3nm/CdS3ML/ZnS1-2ML")
axs[1].plot(Meijerink_coreShell_T, Meijerink_coreShell_Width, "s:", label="Meijerink, ensemble meas. deconvoluted, CdSe3.6nm/CdS1-2layers, 4.6-5.4nm total")
axs[1].plot(Bawendi_coreShell_T, Bawendi_coreShell_Width, "d:", label="Bawendi, SDL, CdSe/ZnS, 4.3nm total")

axs[0].legend()
axs[1].legend()

fig.tight_layout()
fig.savefig("experimental_summary.pdf")