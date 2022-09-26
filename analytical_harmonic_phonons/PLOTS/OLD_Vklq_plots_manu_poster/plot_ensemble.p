set term qt size 1000, 400
set multiplot
set grid
set key outside right tmargin

set origin 0,0
set size 0.33,1
set xlabel "Energy [meV]"
set ylabel "Abs [Arb. unit]"
set title "Low T lineshape, CdSe/CdS"
plot [1900:2200] [0:1] "experimental_digitized/ensemble/Meijerink_coreShell_4K_shape.dat" u (4150-$1):2 w l lt 7 t "Meijerink JPCC, 4K, ensemble", \
"scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian50fs_Spectrum_QM_E17_T=5_aIndex0_bIndex0.dat" \
u ($1*1000):2 w l lc "black" t "Lin, adjusted, 5K", \
"scale_coupling_coreShell/no_scaling_gaussian/results/gaussian50fs_Spectrum_QM_E17_T=5.dat"\
u ($1*1000):2 w l lt 6 t "Lin, original, 5K"

set origin 0.33,0
set size 0.33,1
set xlabel "Energy [meV]"
set ylabel "Abs [Arb. unit]"
set title "High T lineshape, CdSe/CdS"
plot [1900:2200] [0:1] "experimental_digitized/ensemble/Meijerink_coreShell_295K_shape.dat" u (4100-$1):2 w l lt 7 t "Meijerink JPCC, 295K, ensemble", \
"scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian50fs_Spectrum_QM_E17_T=300_aIndex0_bIndex0.dat" \
u ($1*1000):2 w l lc "black" t "Lin, adjusted, 300K", \
"scale_coupling_coreShell/no_scaling_gaussian/results/gaussian50fs_Spectrum_QM_E17_T=300.dat"\
u ($1*1000):2 w l lt 6 t "Lin, original, 300K"

set origin 0.66,0
set size 0.33,1
set xlabel "Temperature [K]"
set ylabel "FWHM [meV]"
set title "PL Linewidth"
plot "experimental_digitized/ensemble/Meijerink_coreShell_inhomo_width.dat" w lp lt 7 t "Meijerink JPCC, inhomo", \
"scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/WidthFreqShift_vs_T_gaussian50fs.dat" \
u 1:($4*1000) w lp lc "black" t "Lin, adjusted", \
"scale_coupling_coreShell/no_scaling_gaussian/results/WidthFreqShift_vs_T_gaussianBroadening_50fs.dat" \
u 1:($5*1000) w lp lt 6 t "Lin, original"



unset multiplot