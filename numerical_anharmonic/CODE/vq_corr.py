import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('eig_file', type=str)
parser.add_argument('w_file', type=str)
parser.add_argument('vel_file', type=str)
parser.add_argument('vqc_file', type=str)
args = parser.parse_args()

# # # # parameters # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# number of atoms
natoms = 444
# number of dimensions
ndim = natoms*3

# # # # frequencies and modes # # # # # # # # # # # # # # # # # # # # # # # # # # #
freq = np.array([line.strip().split()[1] for line in open(args.w_file, 'r')]).astype(np.float)
eig = np.array([line.strip().split() for line in open(args.eig_file, 'r')]).astype(np.float).transpose()

print('Finished reading eigenvalues...\n')

# # # # velocities # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# read velocities from file
lines = np.array([line.strip().split() for line in open(args.vel_file, 'r')])
vx = []
vy = []
vz = []
count = 0
while count < len(lines):
    if len(lines[count]) > 1:
        if lines[count][1] == "ATOMS":
            count += 1
            vx.append(np.array([l[2] for l in lines[count:count+natoms]]))
            vy.append(np.array([l[3] for l in lines[count:count+natoms]]))
            vz.append(np.array([l[4] for l in lines[count:count+natoms]]))
            count += natoms
        else:
            count += 1
    else:
        count += 1

print('Finished reading velocities...\n')

# mass (kg)
at_idx = np.argsort(np.array([l[0] for l in lines[9:9+natoms]]).astype(np.int))
mass = ((np.array([l[1] for l in lines[9:9+natoms]])[at_idx]).astype(np.float)*1e-3/6.022140857e23)
mass3 = np.repeat(mass, 3)

# velocities (ang / ps)
vx = np.array(vx).astype(np.float)
vy = np.array(vy).astype(np.float)
vz = np.array(vz).astype(np.float)
v = np.dstack((vx, vy, vz))
# velocities in cartesian coordinates (m / s)
v = v.reshape((v.shape[0], ndim))*100

# velocities in normal mode coordinates (sqrt(kg)*m / s)
v_q = np.stack([np.matmul(eig, np.sqrt(mass3)*v[i,:]) for i in range(v.shape[0])])

print('Finished computing velocities...\n')

# # # # correlation function # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
corr_len = 20000 # 200ps
avg_len = int(1.5*corr_len)
corr = np.empty((3*natoms, corr_len))
for i in range(corr_len):
    corr[:,i] = np.mean(v_q[:avg_len,:]*v_q[i:i+avg_len,:], axis=0)

print('Finished computing correlation functions...\n')

# # # # write to file
with open(args.vqc_file, 'w') as f:
    for i in range(corr.shape[1]):
        for j in range(corr.shape[0]):
            f.write(str(corr[j,i]) + ' ')
        f.write('\n')

print('Finished writing correlation functions...\n')
