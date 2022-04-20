# Dipti Jasrasaria
# October 2020
#
# LAMMPS routine to minimize structure and compute dynamical matrix

# Cd/Zn/Se/S/Te/Hg nanostructure with Stillinger-Weber interaction potentials
# CdTeZnSeHgS0.sw file has the parameters

# variables and units
atom_modify             map yes
units			metal # http://lammps.sandia.gov/doc/units.html A, ps, eV, K, 0.001ps tstep, 2.0 A skin
variable		atom_list string "Cd Se"
variable		start_temp index 0.001    	# starting temperature of simulation
variable		time_step index 0.001       # 0.001 ps = 1.0 fs timestep - default = 0.001 ps for metal units
variable		seed index 359592          # random seed used in initializing the velocities
dimension		3  # If 2D then z must be periodic
boundary		f f f  # periodic (p), fixed (f) 

# set up atom positions and atom parameters (masses, charges, etc.)
atom_style		atomic	# using default values (coarse-grain liquids, solids, metals), charge specified in read_data file
newton			on # Calculated forces individually (less communication but more computing with it off) - on for Tersoff
read_data		lammpsconf.par # file with the atom charges, positions, etc.

# set up the interactions between the atoms
pair_style		sw   # type of pair potential 
pair_coeff		* * CdTeZnSeHgS0.sw ${atom_list} # atom type i, atom type j, filename, atom symbols in order of their atom type

# set up neighbor list information
neighbor		5.0  bin 	# skin(A) type - make neighbor list containing atoms within cutoff+skin in A
neigh_modify	        every 1 delay 0 check yes 	# use this for minimization

# set up the thermodynamic printing in the minimization
thermo 			3 	# outputs the thermo info every 3 steps
thermo_style	        custom step temp press pe ke etotal 	# prints thermo props at every thermo step

# set up and run the minimization
min_style		cg 	# conjugate-gradient minimization
dump			init_min all xyz 2 init_Mintraj.xyz # dumps the atom positions every 2 timesteps
minimize		0.0 1.0e-8 1000 100000
undump 			init_min 	# undumps the dumping of initial minimization

# set up the dummy md run 
neigh_modify	        delay 5
timestep		${time_step}  # 0.001 ps = 1.0 fs timestep - the default is 0.001 ps for metal units
run_style		verlet
velocity		all create ${start_temp} ${seed} dist gaussian  # initial gaussian velocity distr centered @ start_temp

# dump the initial minimized structure
dump			init_min_only all xyz 2 init-min.xyz 
run			1
undump			init_min_only
