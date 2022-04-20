import numpy as np
import matplotlib.pyplot as plt

HBAR = 6.582119569e-16     # eV s
KB = 8.617333262145e-5     # eV K^-1

reorg_E = 50.192e-3        # eV
E_shift = 2.18538          # eV

T_range = np.array([1, 2, 10, 200, 300])

fig, axs = plt.subplots(len(T_range), 3, figsize=(18, 4*len(T_range)))

for i in range(len(T_range)): 
    T = T_range[i]
    F_file = "F_QM_E17_T=%d.dat" % T
    I_file = "Spectrum_QM_E17_T=%d.dat" % T

    t_range = np.zeros(0)
    F_t_Re_range = np.zeros(0)
    f = open(F_file, "r")
    while True:
        line = f.readline()
        if not line:
            break
        if line[0]!= "#":
            t_range = np.append(t_range, float(line.split()[0]))
            F_t_Re_range = np.append(F_t_Re_range, float(line.split()[1]))
    f.close()

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
    
    # 2nd order derivative
    deriv_x = np.zeros(0)
    deriv_y = np.zeros(0)
    for j in range(len(E_range)):
        if (j!=0) and (j!=len(E_range)-1):
            deriv_x = np.append(deriv_x, E_range[j])
            deriv_y = np.append(deriv_y, (I_E_range[j+1]+I_E_range[j-1]-2*I_E_range[j])/((E_range[j+1]-E_range[j])**2))


    axs[i,0].plot(t_range, F_t_Re_range, label="QM_E17")
    # axs[i,0].plot(t_range, F_Gaussian, label="Static limit")
    axs[i,0].grid()
    axs[i,0].legend()
    axs[i,0].set(title="Dephasing function at T=%d K"%T, xlabel="time [s]", ylabel=r"$\langle F(t) \rangle$ [unitless]", xlim=[0, 1e-12])

    axs[i,1].plot(E_range, I_E_range, label="QM_E17")
    # axs[i,1].plot(E_range, I_Gaussian, label="Static limit")
    # axs[i,1].plot(deriv_x, deriv_y / np.max(deriv_y), ":", label="second derivative, normalized")
    axs[i,1].grid()
    axs[i,1].legend()
    axs[i,1].set(xlabel=r"Energy [eV]", ylabel=r"$I(\omega)$", xlim=[2.0,2.4])

    axs[i,2].plot(deriv_x, deriv_y)
    axs[i,2].grid()
    axs[i,2].set(title="Second derivative of spectrum", xlabel=r"Energy [eV]", ylabel=r"$ \frac {d^2 I}{dE^2}$", xlim=[2.23,2.30], ylim=[-10,10])


axs[3,2].set(xlim=[2.37, 2.50])
axs[4,2].set(xlim=[2.37, 2.50])


fig.tight_layout()
fig.savefig("allPhonon_QM_E17.pdf")
