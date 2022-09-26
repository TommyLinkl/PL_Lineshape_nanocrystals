set term qt size 600,500
set grid

set multiplot
set origin 0,0
set size 1,1
set xlabel "Energy [meV]"
set ylabel "Abs [Arb. unit]"
set title "PL spectrum for CdSe/CdS QD at low T. Adjusted couplings"
plot [2010:2130] "experimental_digitized/single_NC/Hendrik_4K_shape.dat" u (4150-$1):2 w l lt 7 t "Hendrik, CdSe/CdS/ZnS, 4K", \
"scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian800fs_Spectrum_QM_E17_T=5_aIndex0_bIndex0.dat" \
u ($1*1000):2 w l lc "black" t "Lin, Trial3BestFit, CdSe/CdS, 5K, w/ broadening of 0.8ps", \
"scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/highRes_Spectrum_QM_E17_T=5_aIndex0_bIndex0.dat" \
u ($1*1000):2 w l lt 6 t "Lin, Trial3BestFit, CdSe/CdS, 5K, high res"


set origin 0.4,0.4
set size 0.55,0.4
unset xlabel
unset ylabel
unset title
plot [2029:2037] [-0.2:1.2] "scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian800fs_Spectrum_QM_E17_T=5_aIndex0_bIndex0.dat" \
u ($1*1000):2 w l lc "black" t "", \
"scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/highRes_Spectrum_QM_E17_T=5_aIndex0_bIndex0.dat" \
u ($1*1000):($2*7) w l lt 6 t "x7"

unset multiplot