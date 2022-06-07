for i in {0..99}; do
    mkdir $i
    echo $i
    #((j=${i}+10))
    #echo $j
    #mkdir $j

    # get snapshot from MD trajectory   
    ((a=10*453*$i+1))        # ((a=10*879*$i+1))
    ((b=$a+452))             # ((b=$a+878))
    ((c=$b+1))
    #((k=$i+1))
    #((a=5*195*$k+1))
    #((b=$a+194))
    #((c=$b+1))

    sed -n "${a},${b}p;${c}q" ../langevinMD/traj-langevin.xyz > $i/snapshot$i.xyz
    #sed -n "${a},${b}p;${c}q" traj-langevin.xyz > $j/snapshot$j.xyz

    # set up LAMMPS input file
    sed "s/snum/${RANDOM}/g" in-tmp.sw > tmp
    sed "s/num/${i}/g" tmp > $i/in-$i.sw
    #sed "s/num/${j}/g" tmp > $j/in-$j.sw

    cd $i
    #cd $j

    # make LAMMPS conf file
    python3 ~/tommylin/PL_Lineshape_nanocrystals/numerical_anharmonic/3.0nmCdSe/trajectories/xyzToLammpsConf.py snapshot$i.xyz
    #/home/djasras/anaconda/bin/python /home/djasras/scripts/xyzToLammpsConf.py snapshot$j.xyz
    mv lammpsconf_test.par lammpsconf.par

    # run LAMMPS
    lmp -in in-$i.sw > run_lammps-$i.dat
    #mpirun -np 16 lmp_mpi -in in-$j.sw > run_lammps-$j.dat
    cd ../
done

