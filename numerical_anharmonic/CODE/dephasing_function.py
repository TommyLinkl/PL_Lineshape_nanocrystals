import numpy as np
import matplotlib.pyplot as plt
import argparse

HBAR = 1.054571817e-34     # Js

parser = argparse.ArgumentParser()
parser.add_argument('eig_file', type=str)
parser.add_argument('w_file', type=str)
parser.add_argument('Vklq_file', type=str)
parser.add_argument('equil_file', type=str)
parser.add_argument('traj_file', type=str)
parser.add_argument('dephasing_file', type=str)
args = parser.parse_args()

######## parameters #############################################################
# number of atoms
natoms = 444
# number of dimensions
ndim = natoms*3
dt = 1e-15     # s

# # # # frequencies and modes # # # # # # # # # # # # # # # # # # # # # # # # # # #
freq = np.array([line.strip().split()[1] for line in open(args.w_file, 'r')]).astype(float)
freq *= 2 * np.pi * 1e12
eig = np.array([line.strip().split() for line in open(args.eig_file, 'r')]).astype(float).transpose()
print('Finished reading eigenvalues...\n')

####### coupling ##################################################################
V = np.zeros(len(freq))
f = open(args.Vklq_file, "r")
while True:
    line = f.readline()
    if not line:
        break
    if (line[0]!= "#") and (int(line.split()[1])==0) and (int(line.split()[2])==0):
        V[int(line.split()[0])] = float(line.split()[3])
f.close()
print('Finished reading couplings...\n')

####### positions ################################################################
####### equilibrium position ####################################################
x_0 = np.zeros(natoms)
y_0 = np.zeros(natoms)
z_0 = np.zeros(natoms)
count = 0
f = open(args.equil_file, "r")
line = f.readline()
line = f.readline()
while True:
    line = f.readline()
    if not line:
        break
    if (line[0]!= "#"):
        x_0[count] = float(line.split()[1])
        y_0[count] = float(line.split()[2])
        z_0[count] = float(line.split()[3])
        count += 1
f.close()
r_0 = np.dstack((x_0, y_0, z_0))
r_0 = r_0.reshape((r_0.shape[0], ndim))*1e10
print('Finished reading equilibrium positions...\n')

################################################################################
# read trajecory positions from file
lines = np.array([line.strip().split() for line in open(args.traj_file, 'r')])
x = []
y = []
z = []
count = 0
while count < len(lines):
    if len(lines[count]) > 1:
        if lines[count][1] == "ATOMS":
            count += 1
            x.append(np.array([l[2] for l in lines[count:count+natoms]]))
            y.append(np.array([l[3] for l in lines[count:count+natoms]]))
            z.append(np.array([l[4] for l in lines[count:count+natoms]]))
            count += natoms
        else:
            count += 1
    else:
        count += 1
print('Finished reading positions...\n')

# mass (kg)
at_idx = np.argsort(np.array([l[0] for l in lines[9:9+natoms]]).astype(int))
mass = ((np.array([l[1] for l in lines[9:9+natoms]])[at_idx]).astype(float)*1e-3/6.022140857e23)
mass3 = np.repeat(mass, 3)

# positions #(ang / ps)
x = np.array(x).astype(float)
y = np.array(y).astype(float)
z = np.array(z).astype(float)
r = np.dstack((x, y, z))
# positions in cartesian coordinates #(m)
r = r.reshape((r.shape[0], ndim))*1e10
# print(r.ndim)
# print(r.shape)
r = r - np.tile(r_0, (r.shape[0], 1))

# positions in normal mode coordinates (sqrt(kg)*m / s)
r_q = np.stack([np.matmul(eig, np.sqrt(mass3)*r[i,:]) for i in range(r.shape[0])])
print(r_q)
print('Finished computing normal mode positions...\n')

##########################################################################################
# # # # dephasing function # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
####### delete the first 6 phonon modes ##################################################
freq = freq[6:]
r_q = r_q[:,6:]
V = V[6:]
coupling = V / (freq**2)

delta = np.matmul(r_q, np.reshape(freq**2*coupling, (-1,1))).flatten()
print(delta.shape)
reorg_E = np.sum(0.5 * coupling**2 * freq**2)
print(reorg_E)

integral = np.zeros(len(delta))
F_t_Re = np.zeros(len(delta))
F_t_Im = np.zeros(len(delta))
for t in range(len(integral)):
    integral[t] = np.sum(delta[0:t]) * dt
    F_t_Re[t] = np.cos(-(integral[t]+ reorg_E*t*dt)/HBAR)
    F_t_Im[t] = np.sin(-(integral[t]+ reorg_E*t*dt)/HBAR)
print(integral)
print(F_t_Re)
print(F_t_Im)
print('Finished computing correlation functions...\n')

# # # # write to file
with open(args.dephasing_file, 'w') as f:
    f.write("# time (s)          F_t_Re             F_t_Im \n")
    for i in range(len(F_t_Re)):
        f.write("%.20f         %.20f         %.20f\n" % (i*dt, F_t_Re[i], F_t_Im[i]))

print('Finished writing correlation functions...\n')

