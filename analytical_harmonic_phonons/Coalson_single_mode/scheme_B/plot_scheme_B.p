set term qt size 1400, 700

set multiplot

set size 0.5, 1
set origin 0, 0
set grid
set xlabel "Time [Reduced units]"
set ylabel "Dephasing function [unitless]"
set title "Single mode correlation function results - Scheme B"
plot [0:1] "cutoff/F_QM_E17_T=0_extraPhase.dat" u 1:4 w l t "Lin - F_{Re}(t)*phase, damped", \
"cutoff/F_QM_E17_T=0_extraPhase.dat" u 1:5 w l t "Lin - F_{Im}(t)*phase, damped", \
"digitizedData_schemeB_g=0/corr_Re_fig4.dat" t "Coalson digitized fig4 - C_{Re}(t)", \
"digitizedData_schemeB_g=0/corr_Im_fig4.dat" t "Coalson digitized fig4 - C_{Im}(t)", \
"digitizedData_schemeB_g=0/corr_Re_fig5.dat" t "Coalson digitized fig5 - C_{Re}(t)", \
"digitizedData_schemeB_g=0/corr_Im_fig5.dat" t "Coalson digitized fig5 - C_{Im}(t)"
# exp(4.5*(cos(x)-1))*cos(-4.5*sin(x)+4*x) t "Coalson equation - C_{Re}(t)", \
# exp(4.5*(cos(x)-1))*sin(-4.5*sin(x)+4*x) t "Coalson equation - C_{Im}(t)"

set size 0.25, 0.48
set origin 0.24, 0.3
set xlabel "Time [Reduced units]"
set ylabel "Dephasing function [unitless]"
unset title
plot "F_QM_E17_T=0_extraPhase.dat" w l t "Lin - F_{Re}(t)*phase, no damping", \
"F_QM_E17_T=0_extraPhase.dat" u 1:3 w l t "Lin - F_{Im}(t)*phase, no damping"


set size 0.5, 1
set origin 0.5, 0
set xlabel "Energy [Reduced units]"
set ylabel "Abs"
set title "Single mode spectrum results - Scheme B"
plot [-30:5] "cutoff/Spectrum_QM_E17_T=0_extraPhase.dat" u (-$1/0.004136):($2*0.1136) w l t "Lin - Damped spectrum",\
"Spectrum_QM_E17_T=0_extraPhase.dat" u (-$1/0.004136):($2*0.1136) w l t "Lin - undamped spectrum",\
"digitizedData_schemeB_g=0/spectrum_fig6.dat" u ($1-12.5):2 t "Coalson digitized fig3 - Spectrum, manually shifted by 12.5[reduced Energy]"

unset multiplot
