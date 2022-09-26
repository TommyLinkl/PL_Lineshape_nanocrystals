set term qt size 1200, 800
set multiplot
set grid

set origin 0,0.5
set size 0.5,0.5
set xlabel "Temperature [K]"
set ylabel "FWHM [meV]"
set title "Original couplings, core-shell"
plot "<echo '10 3.19'" t "Klimov PRL93, CdSe/ZnS Dot1, 10K", \
"<echo '10 2.63'" t "Klimov PRL93, CdSe/ZnS Dot2, 10K", \
"../scale_coupling_coreShell/no_scaling_gaussian/WidthFreqShift_vs_T_gaussianBroadening_20ps.dat" \
u 1:($5*1000) w lp t "Lin, original, CdSe/CdS, w/ broadening of 20ps"
set origin 0.25, 0.6
set size 0.2, 0.2
unset xlabel
unset ylabel
unset title
plot [0:20] [-2:6] "<echo '10 3.19'" t "", \
"<echo '10 2.63'" t "", \
"../scale_coupling_coreShell/no_scaling_gaussian/WidthFreqShift_vs_T_gaussianBroadening_20ps.dat" \
u 1:($5*1000) w lp t ""

set origin 0,0
set size 0.5,0.5
set xlabel "Energy [meV]"
set ylabel "Abs"
set title "Low T lineshape"
plot [2000:2120] "single_NC/Klimov_PRL93_2004_Dot1_shape.dat" u (2015.0212+$1):($2-0.08) w l t "Klimov PRL93, CdSe/ZnS Dot1, 10K", \
"single_NC/Klimov_PRL93_2004_Dot2_shape.dat" u (2015.0212+$1):($2-0.04) w l t "Klimov PRL93, CdSe/ZnS Dot2, 10K", \
"../scale_coupling_coreShell/no_scaling_gaussian/gaussian20ps_Spectrum_QM_E17_T=10.dat"\
u ($1*1000):2 w l t "Lin, original, CdSe/CdS, 10K, w/ broadening of 20ps"

#####################################################################################

set origin 0.5,0.5
set size 0.5,0.5
set xlabel "Temperature [K]"
set ylabel "FWHM [meV]"
set title "Scaled couplings, Trial 3 best fit"
plot "<echo '10 3.19'" t "Klimov PRL93, CdSe/ZnS Dot1, 10K", \
"<echo '10 2.63'" t "Klimov PRL93, CdSe/ZnS Dot2, 10K", \
"../scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/WidthFreqShift_vs_T_gaussian1ps.dat" \
u 1:($4*1000) w lp t "Lin, Trial3BestFit, CdSe/CdS, w/ broadening of 1ps"
set origin 0.75, 0.6
set size 0.2, 0.2
unset xlabel
unset ylabel
unset title
plot [0:20] [-2:6] "<echo '10 3.19'" t "", \
"<echo '10 2.63'" t "", \
"../scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/WidthFreqShift_vs_T_gaussian1ps.dat" \
u 1:($4*1000) w lp t ""

set origin 0.5,0
set size 0.5,0.5
set xlabel "Energy [meV]"
set ylabel "Abs"
set title "Low T lineshape"
plot [2000:2120] [0:1] "single_NC/Klimov_PRL93_2004_Dot1_shape.dat" u (2015.0212+$1):($2-0.08) w l t "Klimov PRL93, CdSe/ZnS Dot1, 10K", \
"single_NC/Klimov_PRL93_2004_Dot2_shape.dat" u (2015.0212+$1):($2-0.04) w l t "Klimov PRL93, CdSe/ZnS Dot2, 10K", \
"../scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian1ps_Spectrum_QM_E17_T=10_aIndex0_bIndex0.dat" \
u ($1*1000):2 w l t "Lin, Trial3BestFit, CdSe/CdS, 10K, w/ broadening of 1ps"


unset multiplot