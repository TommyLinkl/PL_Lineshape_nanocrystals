set term qt size 1200, 800
set multiplot
set grid

set origin 0,0.66
set size 0.5,0.33
set xlabel "Temperature [K]"
set ylabel "FWHM [meV]"
set title "Original couplings, core-shell"
plot "single_NC/Hendrik_Dot1_width.dat" w lp t "Hendrik, CdSe/CdS/ZnS", \
"../scale_coupling_coreShell/no_scaling_gaussian/WidthFreqShift_vs_T_gaussianBroadening_20ps.dat" \
u 1:($5*1000) w lp t "Lin, original, CdSe/CdS, w/ broadening of 20ps", \
"../scale_coupling_coreShell/no_scaling_gaussian/WidthFreqShift_vs_T_gaussianBroadening_800fs.dat" \
u 1:($5*1000) w lp t "Lin, original, CdSe/CdS, w/ broadening of 0.8ps"

set origin 0,0.33
set size 0.5,0.33
set xlabel "Energy [meV]"
set ylabel "Abs"
set title "Low T lineshape"
plot [2000:2120] "single_NC/Hendrik_4K_shape.dat" u (4150-$1):2 w l t "Hendrik, CdSe/CdS/ZnS, 4K", \
"../scale_coupling_coreShell/no_scaling_gaussian/gaussian20ps_Spectrum_QM_E17_T=5.dat"\
u ($1*1000):2 w l t "Lin, original, CdSe/CdS, 5K, w/ broadening of 20ps", \
"../scale_coupling_coreShell/no_scaling_gaussian/gaussian800fs_Spectrum_QM_E17_T=5.dat"\
u ($1*1000):2 w l t "Lin, original, CdSe/CdS, 5K, w/ broadening of 0.8ps"

set origin 0,0
set size 0.5,0.33
set xlabel "Energy [meV]"
set ylabel "Abs"
set title "High T lineshape"
plot [1900:2200] "single_NC/Hendrik_290K_shape.dat" u (4050-$1):2 w l t "Hendrik, CdSe/CdS/ZnS, 290K", \
"../scale_coupling_coreShell/no_scaling_gaussian/gaussian20ps_Spectrum_QM_E17_T=300.dat"\
u ($1*1000):2 w l t "Lin, original, CdSe/CdS, 300K, w/ broadening of 20ps", \
"../scale_coupling_coreShell/no_scaling_gaussian/gaussian800fs_Spectrum_QM_E17_T=300.dat"\
u ($1*1000):2 w l t "Lin, original, CdSe/CdS, 300K, w/ broadening of 0.8ps"



set origin 0.5,0.66
set size 0.5,0.33
set xlabel "Temperature [K]"
set ylabel "FWHM [meV]"
set title "Scaled couplings, Trial 3 best fit"
plot "single_NC/Hendrik_Dot1_width.dat" w lp t "Hendrik, CdSe/CdS/ZnS", \
"../scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/WidthFreqShift_vs_T_gaussian800fs.dat" \
u 1:($4*1000) w lp t "Lin, Trial3BestFit, CdSe/CdS, w/ broadening of 0.8ps"

set origin 0.5,0.33
set size 0.5,0.33
set xlabel "Energy [meV]"
set ylabel "Abs"
set title "Low T lineshape"
plot [2000:2120] "single_NC/Hendrik_4K_shape.dat" u (4150-$1):2 w l t "Hendrik, CdSe/CdS/ZnS, 4K", \
"../scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian800fs_Spectrum_QM_E17_T=5_aIndex0_bIndex0.dat" \
u ($1*1000):2 w l t "Lin, Trial3BestFit, CdSe/CdS, 5K, w/ broadening of 0.8ps"

set origin 0.5,0
set size 0.5,0.33
set xlabel "Energy [meV]"
set ylabel "Abs"
set title "High T lineshape"
plot [1900:2200] "single_NC/Hendrik_290K_shape.dat" u (4050-$1):2 w l t "Hendrik, CdSe/CdS/ZnS, 290K", \
"../scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian800fs_Spectrum_QM_E17_T=300_aIndex0_bIndex0.dat" \
u ($1*1000):2 w l t "Lin, Trial3BestFit, CdSe/CdS, 300K, w/ broadening of 0.8ps"



unset multiplot