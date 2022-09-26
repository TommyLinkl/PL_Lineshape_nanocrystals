set term qt size 1000, 900
set grid

set multiplot
set origin 0,0.66
set size 0.33,0.33
unset xlabel 
unset ylabel
set title "Guassian 100fs, 0K, CdSe/CdS"
plot [2:2.15] [-0.1:1.0] "scale_coupling_coreShell/no_scaling_gaussian/results/gaussian100fs_Spectrum_QM_E17_T=0.dat"\
 w l lc "black" t "Original couplings", \
"scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian100fs_Spectrum_QM_E17_T=0_aIndex0_bIndex0.dat"\
 w l lc "red" t "Adjusted couplings"

set origin 0.33,0.66
set size 0.33,0.33
unset xlabel 
unset ylabel
set title "Guassian 200fs, 0K"
plot [2:2.15] [-0.1:1.0] "scale_coupling_coreShell/no_scaling_gaussian/results/gaussian200fs_Spectrum_QM_E17_T=0.dat" w l lc "black" t "", \
"scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian200fs_Spectrum_QM_E17_T=0_aIndex0_bIndex0.dat"\
 w l lc "red" t ""

set origin 0.66,0.66
set size 0.33,0.33
unset xlabel 
unset ylabel
set title "Guassian 800fs, 0K"
plot [2:2.15] [-0.1:1.0] "scale_coupling_coreShell/no_scaling_gaussian/results/gaussian800fs_Spectrum_QM_E17_T=0.dat" w l lc "black" t "", \
"scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian800fs_Spectrum_QM_E17_T=0_aIndex0_bIndex0.dat"\
 w l lc "red" t ""

set origin 0,0.33
set size 0.33,0.33
unset xlabel 
unset ylabel
set title "Guassian 2ps, 0K"
plot [2:2.15] [-0.1:1.0] "scale_coupling_coreShell/no_scaling_gaussian/results/gaussian2ps_Spectrum_QM_E17_T=0.dat" w l lc "black" t "", \
"scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian2ps_Spectrum_QM_E17_T=0_aIndex0_bIndex0.dat"\
 w l lc "red" t ""

set origin 0.33,0.33
set size 0.33,0.33
unset xlabel 
unset ylabel
set title "Guassian 4ps, 0K"
plot [2:2.15] [-0.1:1.0] "scale_coupling_coreShell/no_scaling_gaussian/results/gaussian4ps_Spectrum_QM_E17_T=0.dat" w l lc "black" t "", \
"scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian4ps_Spectrum_QM_E17_T=0_aIndex0_bIndex0.dat"\
 w l lc "red" t ""

set origin 0.66,0.33
set size 0.33,0.33
unset xlabel 
unset ylabel
set title "Guassian 10ps, 0K"
plot [2:2.15] [-0.1:1.0] "scale_coupling_coreShell/no_scaling_gaussian/results/gaussian10ps_Spectrum_QM_E17_T=0.dat" w l lc "black" t "", \
"scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian10ps_Spectrum_QM_E17_T=0_aIndex0_bIndex0.dat"\
 w l lc "red" t ""

set origin 0,0
set size 0.33,0.33
unset xlabel 
unset ylabel
set title "Guassian 20ps, 0K"
plot [2:2.15] [-0.1:1.0] "scale_coupling_coreShell/no_scaling_gaussian/results/gaussian20ps_Spectrum_QM_E17_T=0.dat" w l lc "black" t "", \
"scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian20ps_Spectrum_QM_E17_T=0_aIndex0_bIndex0.dat"\
 w l lc "red" t ""

set origin 0.33,0
set size 0.33,0.33
unset xlabel 
unset ylabel
set title "Guassian Inf, i.e. highRes, 0K"
plot [2:2.15] [-0.1:1.0] "scale_coupling_coreShell/no_scaling_gaussian/results/Spectrum_QM_E17_T=0_shrink#1_highRes.dat"\
 w l lc "black" t "Original couplings", \
"scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/highRes_Spectrum_QM_E17_T=0_aIndex0_bIndex0.dat"\
 w l lc "red" t "Adjusted couplings"

set origin 0.66,0
set size 0.33,0.33
unset xlabel
unset ylabel
set title "Comparison to Hendrik's experimental results"
plot [2:2.15] [-0.1:1.0] "scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian800fs_Spectrum_QM_E17_T=5_aIndex0_bIndex0.dat"\
 w l lc "red" t "Adjusted, 5K, CdSe/CdS", \
"experimental_digitized/single_NC/Hendrik_4K_reversedAndShifted_shape.dat" w l lc "blue" t "Hendrik, 4K, CdSe/CdS/ZnS"

unset multiplot