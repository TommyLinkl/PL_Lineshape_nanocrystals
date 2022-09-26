import numpy as np
import math
from scipy.fft import fft, ifft
from scipy.integrate import simpson
from scipy.integrate import trapezoid
from scipy.special import erf
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib import rc
matplotlib.rcParams['font.family'] = "serif"
matplotlib.rcParams['font.sans-serif'] = "Times New Roman"
matplotlib.rcParams['font.size'] = 12
matplotlib.rcParams['mathtext.fontset'] = "stix"
rc('axes', linewidth=0.4)
matplotlib.rcParams['xtick.major.width'] = 0.4
matplotlib.rcParams['ytick.major.width'] = 0.4


KB = 1.380649e-23                # unit: J/K, SI units
PI = 3.14159265359
HBAR = 1.054571817e-34           # Unit: Js, SI units
JTOEV = 6.24150907446076e18 
AUTOEV = 27.211396641308 
AUTOJ = 4.359744e-18

def coth(x):
    return 1/np.tanh(x)

def func_n(x, T): 
    BETA = 1 / (KB * T)
    # print(x)
    return 1 / (np.exp(BETA * HBAR * x) - 1)

def read_w(filename):
    w = np.zeros(0)
    f = open(filename, "r")
    while True:
        line = f.readline()
        if not line:
            break
        if line[0]!= "#":
            w = np.append(w, float(line.split()[1]))
    f.close()
    w *= 1e12 * (2*PI)     # in 1/s
    # print(w)
    return w

def read_Vklq(filename, nPhonon, nExc, shrink):
    coupling = np.zeros((nPhonon, nExc))
    f = open(filename, "r")
    while True:
        line = f.readline()
        if not line:
            break
        if (line[0]!= "#") and (int(line.split()[1])==int(line.split()[2])):
            coupling[int(line.split()[0]), int(line.split()[1])] = float(line.split()[3])
    f.close()
    # print(coupling*shrink)
    return (coupling*shrink)      # in sqrt(J) / s

def reorg_E(delta_range, omega_range):
    E = np.sum(0.5 * omega_range**2 * delta_range**2)
    print("reorg_E = %f J, %f eV, %6f meV" % (E, E*JTOEV, E*JTOEV*1000))
    return E   # In units of J
    
def write_spectrum(filename, w_range, spectrum):
    f = open(filename, "w")
    f.write("# Frequency w         I(w)\n")
    for i in range(len(w_range)):
        f.write("%.20f       %.20f\n" % (w_range[i], spectrum[i]))
    f.close()
    return

def calc_spectrum_QM_E17(T, t_range, w_range, cutoff, slope, exc_E, phonon_filename, \
    coupling_filename, I_file, FWHM_file, shrink):
    # give T in K, t_range in s, w_range in eV, cutoff in s, exc_E in J

    dt = t_range[1] - t_range[0]      # assuming equal spacing

    # Reading frequency and coupling matrix element data
    nExc = 1
    phonons = read_w(phonon_filename)
    coupling = read_Vklq(coupling_filename, len(phonons), nExc, shrink)
    phonons = phonons[6:]
    V = coupling[6:, 0]
    delta_range = (V) / phonons**2
    print("DONE reading")

    # Calculating the dephasing function
    (F_t_Re, F_t_Im) = calc_dephasing_F_E17(t_range, delta_range, phonons, T)
    ########################################################################################
    # F_t_Re = np.array([1.0]*len(t_range))
    # F_t_Im = np.zeros(len(t_range))
    ########################################################################################
    print("DONE calc_dephasing_F_E17")
    
    '''
    # Calculating the delta-delta correlation function
    (DC_t_Re, DC_t_Im) = calc_delta_corr_E17(t_range, delta_range, phonons)
    DC_t_Re *= JTOEV**2
    DC_t_Im *= JTOEV**2
    print("DONE calc_delta_corr_E17")
    '''
    DC_t_Re = np.zeros(0)
    DC_t_Im = np.zeros(0)

    # erf correction to the dephasing functions
    # F_t_Re_erf = multiply_erf(cutoff, slope, F_t_Re-np.mean(F_t_Re[-1000:-1]), dt, t_range)
    # F_t_Im_erf = multiply_erf(cutoff, slope, F_t_Im-np.mean(F_t_Im[-1000:-1]), dt, t_range)
    F_t_Re_erf = multiply_erf(cutoff, slope, F_t_Re, dt, t_range)
    F_t_Im_erf = multiply_erf(cutoff, slope, F_t_Im, dt, t_range)
    print("DONE erf_multiplication")

    # FT the dephasing functions to obtain the spectrum from the dipole-dipole correlation function
    delta_E_E17 = exc_E   #  - reorg_E(delta_range, phonons)    # in J
    reorg_E(delta_range, phonons) 
    (spectrum_w, spectrum_I) = spectrum(delta_E_E17, t_range, F_t_Re_erf, F_t_Im_erf, w_range)
    print("DONE spectrum")
    write_spectrum(I_file, spectrum_w, spectrum_I)

    (FWHM, l_index, r_index) = calc_FWHM(spectrum_w, spectrum_I)
    print("DONE calculating FWHM\n")
    # write_FWHM(FWHM_file, spectrum_w, spectrum_I, l_index, r_index)

    freq_shift = calc_freqShift(spectrum_w, spectrum_I, exc_E)

    return (spectrum_w, spectrum_I, FWHM, freq_shift, F_t_Re, F_t_Im, F_t_Re_erf, F_t_Im_erf, DC_t_Re, DC_t_Im)



###########################################################################

def delta_to_Gaussian(center, gaussian_std, x):
    y = 1/(gaussian_std * np.sqrt(2*np.pi)) * np.exp(-(x - center)**2 / (2*gaussian_std**2))
    return y

def SK_spectral_density(w_range, h_range, gaussian_std):
    x_range = np.arange(0+gaussian_std/10, max(w_range)+3*gaussian_std, gaussian_std/10)
    y_range = np.zeros(len(x_range))
    for i in range(len(w_range)):
        w = w_range[i]
        h = h_range[i]
        y_range += delta_to_Gaussian(w, gaussian_std, x_range) * PI/HBAR * h**2
    return (x_range, y_range)

############################################################################################

def Skinner_dephasing_rate(spectral_density_w, spectral_density_J_w, Temp_range): 
    oneOverT2_range = np.zeros(0)
    for T in Temp_range: 
        oneOverT2 = 1/PI * trapezoid( func_n(spectral_density_w,T)* (func_n(spectral_density_w,T)+1) * spectral_density_J_w**2, spectral_density_J_w)
        oneOverT2_range = np.append(oneOverT2_range, oneOverT2)
    return (Temp_range, oneOverT2_range)


# main
shrink = 1.0
Temp_range = np.array([10, 50, 100, 150, 200, 220, 250, 290, 300, 350, 400])

# Reading frequency and coupling matrix element data
nExc = 1
phonons = read_w("w.dat")
coupling = read_Vklq("NEW_Vklq-diabatic.dat", len(phonons), nExc, shrink)
phonons = phonons[6:]
V = coupling[6:, 0]
delta_range = (V) / phonons**2
h_range = (V) * np.sqrt(HBAR) / np.sqrt(2*phonons) 
print("DONE reading")
# print(h_range)

fig, axs = plt.subplots(3,1, figsize=(6, 12))

gaussian_std_range = np.array([0.001, 0.005, 0.01, 0.03, 0.05, 0.1, 0.5]) * 2*PI*1e12  # Hz

for gaussian_std in gaussian_std_range: 
    print(gaussian_std/(2*PI*1e12))
    (spectral_density_w, spectral_density_J_w) = SK_spectral_density(phonons, h_range, gaussian_std)
    if gaussian_std==0.03* 2*PI*1e12: 
        axs[0].plot(spectral_density_w/(2*PI*1e12),spectral_density_J_w*JTOEV, "-k", linewidth=1, label="std=%f"%(gaussian_std/(2*PI*1e12)))
        axs[0].set(xlabel="Frequency [THz]", ylabel=r"Spectral density $J_{SK}\left( \omega \right)$ [eV]")
        axs[0].legend()

    (Temp_range, oneOverT2_range) = Skinner_dephasing_rate(spectral_density_w, spectral_density_J_w, Temp_range)
    axs[1].plot(Temp_range, oneOverT2_range, "-o", label="std=%.3f"%(gaussian_std/(2*PI*1e12)))
    axs[1].set(xlabel="Temperature [K]", ylabel=r"$1/T_2$ [arb.u.]")
    axs[1].legend()
    axs[2].plot(Temp_range, oneOverT2_range / oneOverT2_range[7] * 1/30, "-o", label="std=%.3f, rescaled"%(gaussian_std/(2*PI*1e12)))

Temp_range = np.array([4, 30, 60, 100, 150, 220, 290])
oneOverT2_range = np.array([0, 0, 0, 1/1000, 1/500, 1/50, 1/30])
axs[2].plot(Temp_range, oneOverT2_range, "-+k", label="Best fit dephasing rate")
axs[2].set(xlabel="Temperature [K]", ylabel=r"$1/T_2$ [1/fs]", )
axs[2].legend()

fig.tight_layout()
fig.savefig("SK_T2.pdf")
