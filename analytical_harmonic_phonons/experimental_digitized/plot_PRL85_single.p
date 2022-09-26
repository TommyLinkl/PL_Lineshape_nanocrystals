set term qt size 1200, 800
set multiplot
set grid

set origin 0,0.5
set size 0.5,0.5
set xlabel "Temperature [K]"
set ylabel "FWHM [meV]"
set title "Original couplings, core-shell"
plot "single_NC/PRL85_width.dat" t "PRL85, CdSe/ZnS, 10K", \
"../scale_coupling_coreShell/no_scaling_gaussian/WidthFreqShift_vs_T_gaussianBroadening_5ps.dat" \
u 1:($5*1000) w lp t "Lin, original, CdSe/CdS, w/ broadening of 5ps"
set origin 0.25, 0.6
set size 0.2, 0.2
unset xlabel
unset ylabel
unset title
plot [0:20] [5:30] "single_NC/PRL85_width.dat" t "", \
"../scale_coupling_coreShell/no_scaling_gaussian/WidthFreqShift_vs_T_gaussianBroadening_5ps.dat" \
u 1:($5*1000) w lp t ""

set origin 0,0
set size 0.5,0.5
set xlabel "Energy [meV]"
set ylabel "Abs"
set title "Low T lineshape"
plot [2000:2120] [0:1] "single_NC/PRL85_Dot1_shape.dat" u (4106-$1):($2/0.853094) w l t "PRL85, CdSe/ZnS Dot1, 10K", \
"single_NC/PRL85_Dot2_shape.dat" u (4115-$1):($2/0.950814) w l t "PRL85, CdSe/ZnS Dot2, 10K", \
"single_NC/PRL85_Dot3_shape.dat" u (4128-$1):($2/0.787948) w l t "PRL85, CdSe/ZnS Dot3, 10K", \
"../scale_coupling_coreShell/no_scaling_gaussian/gaussian5ps_Spectrum_QM_E17_T=10.dat"\
u ($1*1000):2 w l t "Lin, original, CdSe/CdS, 10K, w/ broadening of 5ps"

#############################################################################################

set origin 0.5,0.5
set size 0.5,0.5
set xlabel "Temperature [K]"
set ylabel "FWHM [meV]"
set title "Scaled couplings, Trial 3 best fit"
plot "single_NC/PRL85_width.dat" t "PRL85, CdSe/ZnS, 10K", \
"../scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/WidthFreqShift_vs_T_gaussian150fs.dat" \
u 1:($4*1000) w lp t "Lin, Trial3BestFit, CdSe/CdS, w/ broadening of 150fs"
set origin 0.75, 0.6
set size 0.2, 0.2
unset xlabel
unset ylabel
unset title
plot [0:20] [5:30] "single_NC/PRL85_width.dat" t "", \
"../scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/WidthFreqShift_vs_T_gaussian150fs.dat" \
u 1:($4*1000) w lp t ""

set origin 0.5,0
set size 0.5,0.5
set xlabel "Energy [meV]"
set ylabel "Abs"
set title "Low T lineshape"
plot [2000:2120] [0:1] "single_NC/PRL85_Dot1_shape.dat" u (4106-$1):($2/0.853094) w l t "PRL85, CdSe/ZnS Dot1, 10K", \
"single_NC/PRL85_Dot2_shape.dat" u (4115-$1):($2/0.950814) w l t "PRL85, CdSe/ZnS Dot2, 10K", \
"single_NC/PRL85_Dot3_shape.dat" u (4128-$1):($2/0.787948) w l t "PRL85, CdSe/ZnS Dot3, 10K", \
"../scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian150fs_Spectrum_QM_E17_T=10_aIndex0_bIndex0.dat" \
u ($1*1000):2 w lp t "Lin, Trial3BestFit, CdSe/CdS, 10K, w/ broadening of 150fs"


unset multiplot