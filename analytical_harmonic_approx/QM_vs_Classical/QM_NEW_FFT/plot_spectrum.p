set multiplot

set size 1, 1
set origin 0, 0
set grid
set xlabel "Frequency [eV]"
set ylabel "Abs [Arb. unit]"
plot [2:2.5] "Spectrum_QM_E17_T=1.dat" w l, "Spectrum_QM_E17_T=10.dat" w l, "Spectrum_QM_E17_T=100.dat" w l, "Spectrum_QM_E17_T=200.dat" w l, "Spectrum_QM_E17_T=300.dat" w l

set size 0.48, 0.5
set origin 0.45, 0.25
set xtics (0, 100, 200, 300, 400)
set xlabel "Temp [K]"
set y2tics nomirror
set ytics nomirror
set ylabel "Width [meV]"
plot "WidthFreqShift_vs_T.dat" u 1:($2*1000) w lp t "Width" axis x1y1,"WidthFreqShift_vs_T.dat" u 1:($3*1000) w lp t "Freq shift" axis x1y2


unset multiplot