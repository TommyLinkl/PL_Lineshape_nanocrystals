set multiplot
set origin 0,0
set size 1,1 
unset grid
set title "CoreShell QD, time: 100ps"
plot [70:150] "Spectrum_QM_E17_T=0_shrink#1.dat" u ($1*1000):2 w l t "original", "Spectrum_QM_E17_T=0_shrink#2.dat" u ($1*1000):2 w l t "shrink by half"

set origin 0.4, 0.33  
set size 0.55, 0.5 
unset title
set grid
plot [74.5:78] "Spectrum_QM_E17_T=0_shrink#1.dat" u ($1*1000):2 w l t "original", "Spectrum_QM_E17_T=0_shrink#2.dat" u ($1*1000-18.9135):2 w l t "shrink by half, shifted"

unset multiplot
