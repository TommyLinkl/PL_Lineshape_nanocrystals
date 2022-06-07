set term x11 size 1400, 700

set multiplot

set size 0.5, 1
set origin 0, 0
set grid
set xlabel "Time [Reduced units]"
set ylabel "Dephasing function [unitless]"
set title "Single mode correlation function results - Scheme A"
plot [0:2] "cutoff/F_QM_E17_T=0_extraPhase.dat" u 1:4 w l t "Lin - F_{Re}(t)*phase, damped", \
"cutoff/F_QM_E17_T=0_extraPhase.dat" u 1:5 w l t "Lin - F_{Im}(t)*phase, damped", \
"digitizedData_schemeA_g=0/corr_Re_fig1.dat" t "Coalson digitized fig1 - C_{Re}(t)", \
"digitizedData_schemeA_g=0/corr_Im_fig1.dat" t "Coalson digitized fig1 - C_{Im}(t)", \
"digitizedData_schemeA_g=0/corr_Re_fig2.dat" t "Coalson digitized fig2 - C_{Re}(t)", \
"digitizedData_schemeA_g=0/corr_Im_fig2.dat" t "Coalson digitized fig2 - C_{Im}(t)"
# exp(4.5*(cos(x)-1))*cos(-4.5*sin(x)+4*x) t "Coalson equation - C_{Re}(t)", \
# exp(4.5*(cos(x)-1))*sin(-4.5*sin(x)+4*x) t "Coalson equation - C_{Im}(t)"

set size 0.27, 0.5
set origin 0.22, 0.26
set xlabel "Time [Reduced units]"
set ylabel "Dephasing function [unitless]"
unset title
plot [] [-0.2:1.2] "F_QM_E17_T=0_extraPhase.dat" w l t "Lin - F_{Re}(t)*phase, no damping", \
"F_QM_E17_T=0_extraPhase.dat" u 1:3 w l t "Lin - F_{Im}(t)*phase, no damping"


set size 0.5, 1
set origin 0.5, 0
set xlabel "Energy [Reduced units]"
set ylabel "Abs"
set title "Single mode spectrum results - Scheme A"
plot [-12:8] "cutoff/Spectrum_QM_E17_T=0_extraPhase.dat" u (-$1/0.004136):($2*0.1885) w l t "Lin - Damped spectrum",\
"Spectrum_QM_E17_T=0_extraPhase.dat" u (-$1/0.004136):($2*0.1885) w l t "Lin - undamped spectrum",\
"digitizedData_schemeA_g=0/spectrum_fig3.dat" u ($1-4.5):2 t "Coalson digitized fig3 - Spectrum, manually shifted by 4.5[reduced Energy]"

unset multiplot
