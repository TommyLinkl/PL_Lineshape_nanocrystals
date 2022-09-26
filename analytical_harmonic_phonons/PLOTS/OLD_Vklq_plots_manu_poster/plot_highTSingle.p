set term qt size 800,400
set grid
set key outside right tmargin

set multiplot
set origin 0,0
set size 0.5,1
set xlabel "Energy [meV]"
set ylabel "Abs [Arb. unit]"
set title "PL line shape for CdSe/CdS QD at high T"
plot [1900:2200] "scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian800fs_Spectrum_QM_E17_T=300_aIndex0_bIndex0.dat" \
u ($1*1000):2 w l lt -1 t "Lin, adjusted coupling, CdSe/CdS, 300K", \
"scale_coupling_coreShell/no_scaling_gaussian/results/gaussian800fs_Spectrum_QM_E17_T=300.dat"\
u ($1*1000):2 w l lt 6 t "Lin, original, CdSe/CdS, 300K", \
"experimental_digitized/single_NC/Hendrik_290K_shape.dat" u (4050-$1):2 w l lt 7 t "Hendrik, CdSe/CdS/ZnS, 290K"


set origin 0.5,0
set size 0.5,1
set xlabel "Temperature [K]"
set ylabel "FWHM [meV]"
set title "PL linewidth"
plot "scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/WidthFreqShift_vs_T_gaussian800fs.dat" \
u 1:($4*1000) w lp lc "black" t "Lin, adjusted coupling, CdSe/CdS, w/ broadening of 0.8ps", \
"scale_coupling_coreShell/no_scaling_gaussian/results/WidthFreqShift_vs_T_gaussianBroadening_800fs.dat" \
u 1:($5*1000) w lp lt 6 t "Lin, original, CdSe/CdS, w/ broadening of 0.8ps", \
"experimental_digitized/single_NC/Hendrik_Dot1_width.dat" w lp lt 7 t "Hendrik, CdSe/CdS/ZnS", \



unset multiplot