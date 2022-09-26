set term qt size 1200, 800
set multiplot
set grid

set origin 0,0.5
set size 0.5,0.5
set xlabel "Temperature [K]"
set ylabel "FWHM [meV]"
set title "Original couplings, core-shell"
plot "single_NC/Bawendi_PRB76_2007_singleNC_width.dat" t "Bawendi PRB76, CdSe/ZnS, 5K", \
"../scale_coupling_coreShell/no_scaling_gaussian/WidthFreqShift_vs_T_gaussianBroadening_5ps.dat" \
u 1:($5*1000) w lp t "Lin, original, CdSe/CdS, w/ broadening of 5ps"
set origin 0.25, 0.6
set size 0.2, 0.2
unset xlabel
unset ylabel
unset title
plot [0:20] [10:30] "single_NC/Bawendi_PRB76_2007_singleNC_width.dat" t "", \
"../scale_coupling_coreShell/no_scaling_gaussian/WidthFreqShift_vs_T_gaussianBroadening_5ps.dat" \
u 1:($5*1000) w lp t ""

set origin 0,0
set size 0.5,0.5
set xlabel "Energy [meV]"
set ylabel "Abs"
set title "Low T lineshape"
plot [2000:2120] "single_NC/Bawendi_PRB76_2007_singleNC_shape.dat" u (4417-$1):2 w l t "Bawendi PRB76, CdSe/ZnS, 5K", \
"../scale_coupling_coreShell/no_scaling_gaussian/gaussian5ps_Spectrum_QM_E17_T=5.dat"\
u ($1*1000):2 w l t "Lin, original, CdSe/CdS, 5K, w/ broadening of 5ps"

#############################################################################################

set origin 0.5,0.5
set size 0.5,0.5
set xlabel "Temperature [K]"
set ylabel "FWHM [meV]"
set title "Scaled couplings, Trial 3 best fit"
plot "single_NC/Bawendi_PRB76_2007_singleNC_width.dat" t "Bawendi PRB76, CdSe/ZnS, 5K", \
"../scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/WidthFreqShift_vs_T_gaussian150fs.dat" \
u 1:($4*1000) w lp t "Lin, Trial3BestFit, CdSe/CdS, w/ broadening of 150fs"
set origin 0.75, 0.6
set size 0.2, 0.2
unset xlabel
unset ylabel
unset title
plot [0:20] [10:30] "single_NC/Bawendi_PRB76_2007_singleNC_width.dat" t "", \
"../scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/WidthFreqShift_vs_T_gaussian150fs.dat" \
u 1:($4*1000) w lp t ""

set origin 0.5,0
set size 0.5,0.5
set xlabel "Energy [meV]"
set ylabel "Abs"
set title "Low T lineshape"
plot [2000:2120] [0:1] "single_NC/Bawendi_PRB76_2007_singleNC_shape.dat" u (4417-$1):2 w l t "Bawendi PRB76, CdSe/ZnS, 5K", \
"../scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian150fs_Spectrum_QM_E17_T=5_aIndex0_bIndex0.dat" \
u ($1*1000):2 w l t "Lin, Trial3BestFit, CdSe/CdS, 5K, w/ broadening of 150fs"


unset multiplot