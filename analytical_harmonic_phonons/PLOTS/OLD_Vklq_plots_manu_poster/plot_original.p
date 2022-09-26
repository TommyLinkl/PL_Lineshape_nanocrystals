set term qt size 600,500
set grid

set multiplot

set origin 0,0
set size 1,1
set title "PL spectrum for CdSe/CdS QD at 0K, original coupling"
set xlabel "Energy [meV]"
set ylabel "Abs [Arb. unit]"
plot [2020:2150] "scale_coupling_coreShell/no_scaling_gaussian/results/Spectrum_QM_E17_T=0_shrink#1_highRes.dat" u ($1*1000):2 w l t "High Res", \
"scale_coupling_coreShell/no_scaling_gaussian/results/gaussian5ps_Spectrum_QM_E17_T=0.dat" u ($1*1000):2 w l t "Broadened"

set origin 0.4,0.45
set size 0.55,0.4
unset xlabel
unset ylabel
unset title
plot [2030.9:2033.2] "scale_coupling_coreShell/no_scaling_gaussian/results/Spectrum_QM_E17_T=0_shrink#1_highRes.dat" u ($1*1000):2 w l t "", \
"scale_coupling_coreShell/no_scaling_gaussian/results/gaussian5ps_Spectrum_QM_E17_T=0.dat" u ($1*1000):2 w l t ""


unset multiplot