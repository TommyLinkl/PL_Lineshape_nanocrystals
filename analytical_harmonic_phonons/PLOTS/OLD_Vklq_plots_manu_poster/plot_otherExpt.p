set term qt size 1200,500
set grid
set key outside right tmargin

set multiplot
set origin 0,0
set size 0.33,1
set xlabel "Energy [meV]"
set ylabel "Abs [Arb. unit]"
unset title
plot [2000:2120] [0:1] "scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian150fs_Spectrum_QM_E17_T=5_aIndex0_bIndex0.dat" \
u ($1*1000):2 w l lc "black" t "Lin, CdSe/CdS, 5K, w/ broadening of 150fs", \
"experimental_digitized/single_NC/Bawendi_PRB76_2007_singleNC_shape.dat" u (4417-$1):2 w l lt 1 t "Bawendi PRB76, CdSe/ZnS, 5K"


set origin 0.33,0
set size 0.33,1
unset ylabel
plot [2000:2120] "scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian1ps_Spectrum_QM_E17_T=10_aIndex0_bIndex0.dat" \
u ($1*1000):2 w l lc "black" t "Lin, CdSe/CdS, 10K, w/ broadening of 1ps", \
"experimental_digitized/single_NC/Klimov_PRL93_2004_Dot1_shape.dat" u (2015.0212+$1):($2-0.08) w l lt 1 t "Klimov PRL93, CdSe/ZnS Dot1, 10K", \
"experimental_digitized/single_NC/Klimov_PRL93_2004_Dot2_shape.dat" u (2015.0212+$1):($2-0.04+0.6) w l lt 6 t "Klimov PRL93, CdSe/ZnS Dot2, 10K"

set origin 0.66,0
set size 0.33,1
unset ylabel
plot [2000:2120] "scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian150fs_Spectrum_QM_E17_T=10_aIndex0_bIndex0.dat" \
u ($1*1000):2 w l lc "black" t "Lin, CdSe/CdS, 10K, w/ broadening of 150fs", \
"experimental_digitized/single_NC/PRL85_Dot1_shape.dat" u (4106-$1):($2/0.853094) w l lt 1 t "PRL85, CdSe/ZnS Dot1, 10K", \
"experimental_digitized/single_NC/PRL85_Dot2_shape.dat" u (4115-$1):($2/0.950814+0.5) w l lt 6 t "PRL85, CdSe/ZnS Dot2, 10K", \
"experimental_digitized/single_NC/PRL85_Dot3_shape.dat" u (4128-$1):($2/0.787948+1) w l lt 7 t "PRL85, CdSe/ZnS Dot3, 10K"


unset multiplot