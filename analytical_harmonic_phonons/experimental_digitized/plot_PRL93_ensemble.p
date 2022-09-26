set term qt size 1200, 800
set multiplot
set grid

set origin 0,0.5
set size 0.5,0.5
set xlabel "Temperature [K]"
set ylabel "FWHM [meV]"
set title "Original couplings, core-shell"
plot "<echo '10 27.47'" w lp t "Klimov PRL93, CdSe/ZnS, ensemble", \
"../scale_coupling_coreShell/no_scaling_gaussian/WidthFreqShift_vs_T_gaussianBroadening_10ps.dat" \
u 1:($5*1000) w lp t "Lin, original, CdSe/CdS, w/ broadening of 10ps"

set origin 0,0
set size 0.5,0.5
set xlabel "Energy [meV]"
set ylabel "Abs"
set title "Low T lineshape"
plot [1900:2200] [0:1] "ensemble/Klimov_PRL93_2004_ensemble_shape.dat" u ($1+2004):2 w l t "Klimov PRL93, CdSe/ZnS, 10K, ensemble", \
"../scale_coupling_coreShell/no_scaling_gaussian/gaussian10ps_Spectrum_QM_E17_T=10.dat"\
u ($1*1000):2 w l t "Lin, original, CdSe/CdS, 10K, w/ broadening of 10ps"

##########################################################################################

set origin 0.5,0.5
set size 0.5,0.5
set xlabel "Temperature [K]"
set ylabel "FWHM [meV]"
set title "Scaled couplings, Trial 3 best fit"
plot "<echo '10 27.47'" w lp t "Klimov PRL93, CdSe/ZnS, ensemble", \
"../scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/WidthFreqShift_vs_T_gaussian200fs.dat" \
u 1:($4*1000) w lp t "Lin, Trial3BestFit, CdSe/CdS, w/ broadening of 200fs"

set origin 0.5,0
set size 0.5,0.5
set xlabel "Energy [meV]"
set ylabel "Abs"
set title "Low T lineshape"
plot [1900:2200] [0:1] "ensemble/Klimov_PRL93_2004_ensemble_shape.dat" u ($1+2005):2 w l t "Klimov PRL93, CdSe/ZnS, 10K, ensemble", \
"../scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian200fs_Spectrum_QM_E17_T=10_aIndex0_bIndex0.dat" \
u ($1*1000):2 w l t "Lin, Trial3BestFit, CdSe/CdS, 10K, w/ broadening of 200fs"



unset multiplot