for i in {0..9}; do
    echo $i
    
    python3 dephasing_function.py eig.dat w.dat Vklq-diabatic.dat unpassivated_conf.xyz trajectories/$i/traj-nve-$i.xyz dephasingFunc/$i.dat

done
