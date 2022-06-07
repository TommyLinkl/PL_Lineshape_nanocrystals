# # # # writes xyz to lammps conf file
# assumes orthorhombic cell and charges

import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename', type=str)
args = parser.parse_args()

triclinic = False

lines = np.array([line.strip() for line in open(args.filename, 'r')])
i = 0
while (i < len(lines)):
    if len(lines[i].split()) > 1:
        if lines[i].split()[1] == 'NUMBER':
            natoms = int(lines[i+1])
            i += 1
        elif lines[i].split()[1] == 'BOX':
            if lines[i].split()[3] == 'ff':
                xlo = float(lines[i+1].split()[0])
                xhi = float(lines[i+1].split()[1])
                ylo = float(lines[i+2].split()[0])
                yhi = float(lines[i+2].split()[1])
                zlo = float(lines[i+3].split()[0])
                zhi = float(lines[i+3].split()[1])
                i += 3
            elif lines[i].split()[3] == 'pp':
                xlo = float(lines[i+1].split()[0])
                xhi = float(lines[i+1].split()[1])
                ylo = float(lines[i+2].split()[0])
                yhi = float(lines[i+2].split()[1])
                zlo = float(lines[i+3].split()[0])
                zhi = float(lines[i+3].split()[1])
                i += 3
            elif lines[i].split()[3] == 'xy':
                xlo = float(lines[i+1].split()[0]) 
                xhi = float(lines[i+1].split()[1])-float(lines[i+1].split()[2])
                ylo = float(lines[i+2].split()[0])
                yhi = float(lines[i+2].split()[1])
                zlo = float(lines[i+3].split()[0])
                zhi = float(lines[i+3].split()[1])
                tilt = float(lines[i+1].split()[2])
                i += 3
                triclinic=True
        elif lines[i].split()[1] == 'ATOMS':
            atidx = np.argsort(np.array([l.split()[0] for l in lines[i+1:i+1+natoms]]).astype(int))
            mass = np.array([l.split()[1] for l in lines[i+1:i+1+natoms]]).astype(float)[atidx]
            x = np.array([l.split()[2] for l in lines[i+1:i+1+natoms]]).astype(float)[atidx]
            y = np.array([l.split()[3] for l in lines[i+1:i+1+natoms]]).astype(float)[atidx]
            z = np.array([l.split()[4] for l in lines[i+1:i+1+natoms]]).astype(float)[atidx]
            i += natoms+1
        else:
            i += 1
    else:
        i += 1

ntypes = len(np.unique(mass))
atypes = np.empty_like(mass)
chg = np.empty_like(mass)

# Cd
atypes[np.where(mass == 112.411)[0]] = 1
chg[np.where(mass == 112.411)[0]] = 1.18
# Se
atypes[np.where(mass == 78.960)[0]] = 2
chg[np.where(mass == 78.960)[0]] = -1.18
# S
if (np.where(mass == 78.960)[0]).shape[0] == 0:
    atypes[np.where(mass == 32.065)[0]] = 2
else:
    atypes[np.where(mass == 32.065)[0]] = 3
chg[np.where(mass == 32.065)[0]] = -1.18
# LJ fluid particle
atypes[np.where(mass == 40.)[0]] = 1
chg[np.where(mass == 40.)[0]] = -1000.


# # # write file
with open('lammpsconf_test.par', 'w') as f:
    f.write('LAMMPS configuration data file\n\n')
    f.write('  ' + str(natoms) + ' atoms\n\n')
    f.write('  ' + str(ntypes) + ' atom types\n\n')
    f.write('  ' + str(np.round(xlo, 6)) + '  ' + str(np.round(xhi, 6)) + '  xlo xhi\n')
    f.write('  ' + str(np.round(ylo, 6)) + '  ' + str(np.round(yhi, 6)) + '  ylo yhi\n')
    f.write('  ' + str(np.round(zlo, 6)) + '  ' + str(np.round(zhi, 6)) + '  zlo zhi\n')
    if triclinic:
        f.write('  ' + str(np.round(tilt, 6)) + ' 0.000 0.000 xy xz yz\n\n')
    else:
        f.write('\n')
    f.write(' Masses\n\n')
    for i in range(ntypes):
        f.write('  ' + str(i+1) + '   ' + str(mass[np.where(atypes == (i+1))[0]][0]) + '\n')
    f.write('\n')
    f.write(' Atoms\n\n')
    for i in range(natoms):
        #f.write('  ' + str(i+1) + ' ' + str(int(atypes[i])) + ' ' + str(chg[i]) + ' ' + str(x[i]) + ' ' + str(y[i]) + ' ' + str(z[i]) + '\n')
        f.write('  ' + str(i+1) + ' ' + str(int(atypes[i])) + ' ' + str(x[i]) + ' ' + str(y[i]) + ' ' + str(z[i]) + '\n')
