import numpy as np
import matplotlib.pyplot as plt

HBAR = 6.582119569e-16     # eV s
KB = 8.617333262145e-5     # eV K^-1

reorg_E = 50.192e-3        # eV
E_shift = 2.18538          # eV

T_range = np.array([1, 2, 10, 200, 300])

fig, axs = plt.subplots(len(T_range), 1, figsize=(12, 4*len(T_range)))

for i in range(len(T_range)): 
    T = T_range[i]
    F_file = "F_CL_T=%d.dat" % T
    I_file = "Spectrum_CL_T=%d.dat" % T

    t_range = np.zeros(0)
    F_t_Re_range = np.zeros(0)
    '''
    f = open(F_file, "r")
    while True:
        line = f.readline()
        if not line:
            break
        if line[0]!= "#":
            t_range = np.append(t_range, float(line.split()[0]))
            F_t_Re_range = np.append(F_t_Re_range, float(line.split()[1]))
    f.close()
    '''

    E_range = np.zeros(0)
    I_E_range = np.zeros(0)
    f = open(I_file, "r")
    while True:
        line = f.readline()
        if not line:
            break
        if line[0]!= "#":
            E_range = np.append(E_range, float(line.split()[0]))
            I_E_range = np.append(I_E_range, float(line.split()[1]))
    f.close()

    sigma = np.sqrt(HBAR**2 / (2 * KB * T * reorg_E))
    # print(1/(sigma**2))
    F_Gaussian = np.exp(-t_range**2 / (2 * sigma**2))
    I_Gaussian = np.exp(-((E_range-E_shift) / HBAR)**2 * sigma**2 / 2)
    
    '''
    axs[i,0].plot(t_range, F_t_Re_range, label="Classical")
    # axs[i,0].plot(t_range, F_Gaussian, label="Static limit")
    axs[i,0].grid()
    axs[i,0].legend()
    axs[i,0].set(title="Dephasing function at T=%d"%T, xlabel="time [s]", ylabel=r"$\langle F(t) \rangle$ [unitless]")

    axs[i,1].plot(t_range, F_t_Re_range, "*", label="Classical")
    axs[i,1].plot(t_range, F_Gaussian, label="Static limit")
    axs[i,1].grid()
    axs[i,1].legend()
    axs[i,1].set(title="Zooming in to short time", xlabel="time [s]", ylabel=r"$\langle F(t) \rangle$ [unitless]", xlim=[0,5e-14])
    '''

    axs[i].plot(E_range, I_E_range, label="Classical")
    # axs[i,1].plot(E_range, I_Gaussian, label="Static limit")
    axs[i].grid()
    axs[i].legend()
    axs[i].set(title="T=%d"%T, xlabel=r"Energy [eV]", ylabel=r"$I(\omega)$", xlim=[2.15,2.23])

'''
axs[3].set(xlim=[2.17, 2.20])
axs[4].set(xlim=[2.17, 2.20])
'''

fig.tight_layout()
fig.savefig("onePhononTest_Classical.png")
