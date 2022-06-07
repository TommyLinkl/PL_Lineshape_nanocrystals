for i in {5..9}; do
    echo $i
    
    /home/tommylin/PL_Lineshape_nanocrystals/numerical_anharmonic/3.0nmCdSe vq_corr.py eig.dat w.dat trajectories/$i/vel-nve-$i.xyz vq_corr/vel-corr-$i.dat

done
