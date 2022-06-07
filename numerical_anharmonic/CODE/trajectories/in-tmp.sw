# Cd/Zn/Se/S/Te/Hg Nanostructure with Stillinger-Weber Potentials

# The file ~/Pseudopotentials/md/CdTeZnSeHgS0.sw has the parameters for the the interaction

# variables and units
units			metal # http://lammps.sandia.gov/doc/units.html A, ps, eV, K, 0.001ps tstep, 2.0 A skin
variable                atom_list string "Cd Se"
variable		time_step index 0.001           # 0.001 ps = 1.0 fs timestep - default = 0.001 ps for metal units
variable		seed index snum                 # random seed used in initializing the velocities
variable                t equal 300                     # temperature is 300.0 K
dimension		3  # If 2D then z must be periodic
boundary		f f f  # periodic (p), fixed (f) 

# set up atom positions and atom parameters (masses, charges, etc.)
atom_style		atomic	# using default values (coarse-grain liquids, solids, metals), charge specified in read_data file
newton			on # Calculated forces individually (less communication but more computing with it off) - on for Tersoff
read_data		lammpsconf.par # file with the atom charges, positions, etc.

# set up the interactions between the atoms
pair_style              sw                      # type of pair potential 
pair_coeff              * * ../CdTeZnSeHgS0.sw ${atom_list} # atom type i, atom type j, filename, atom symbols in order of their atom type

# set up neighbor list information
neighbor		5.0  bin 	# skin(A) type - make neighbor list containing atoms within cutoff+skin in A
neigh_modify	        every 1 delay 0 check yes 	# use this for minimization

# set up the md run 
neigh_modify	        delay 5
timestep		${time_step}  # 0.001 ps = 1.0 fs timestep - the default is 0.001 ps for metal units
run_style		verlet
velocity		all create ${t} ${seed} mom yes rot yes dist gaussian  # initial gaussian velocity distr centered @ start_temp

thermo 			10 	# outputs the thermo info every 3 steps
thermo_style	        custom step temp press pe ke etotal 	# prints thermo props at every thermo step

reset_timestep 0

# equilibrate  
fix                     1 all langevin ${t} ${t} 0.5 ${seed} zero yes
fix                     2 all nve
run                     5000
unfix                   1
unfix                   2

dump                    equil_dump all custom 2 equil-num.xyz id mass x y z
dump_modify             equil_dump sort id
run                     1
undump                  equil_dump

reset_timestep 0

# positions / velocities
fix                     1 all nve
velocity                all zero linear
velocity                all zero angular

# dump                  vel all custom 10 vel-nve-num.xyz id mass vx vy vz
# dump_modify           vel sort id
dump                    traj all custom 10 traj-nve-num.xyz id mass x y z
dump_modify             traj sort id
run                     10000
# undump                vel
undump                  traj
unfix                   1

# save final structure
dump                    final_dump all custom 2 final-num.xyz id mass x y z
dump_modify             final_dump sort id
run                     1
undump                  final_dump

