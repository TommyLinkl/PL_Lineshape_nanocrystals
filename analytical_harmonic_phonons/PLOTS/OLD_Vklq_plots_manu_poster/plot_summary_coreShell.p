set term qt size 1500, 900
set multiplot
set grid

set origin 0,0
set size 0.5, 0.5
set xlabel "Temperatrue [K]"
set ylabel "FWHM [meV]"
set title ""
set key bottom right
plot [0:400] [0:100]"scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/WidthFreqShift_vs_T_gaussian800fs.dat" u 1:($4*1000) w lp t "Lin Trial3BestFit, CdSe/CdS, w/ broadening of 0.8ps", \
"scale_coupling_coreShell/nonuniform_scaling_trial6/aIndex8_bIndex3_T/results/WidthFreqShift_vs_T_gaussian_800fs.dat" u 1:($5*1000) w lp t "Lin Trial6BestFit, CdSe/CdS, w/ broadening of 0.8ps", \
"experimental_digitized/HendrikDot1_T_Width.dat" w lp t "Hendrik manu, CdSe/CdS/ZnS, expt", \
"experimental_digitized/Bawendi_PRB76_2007_T_Width.dat" w lp t "Bawendi PRB76-2007, CdSe/ZnS, expt"


set origin 0,0.5
set size 0.5, 0.5
set xlabel "Energy [meV]"
set ylabel "Abs [arb. units]"
set title "Single NC results"
set key top right
plot [2000:2120] [0:1] "scale_coupling_coreShell/nonuniform_scaling_trial3/results/gaussian_Spectrum_QM_E17_T=0_aIndex0_bIndex0.dat" u ($1*1000):2 w lp t "Lin Trial3BestFit, CdSe/CdS, w/ broadening of 0.8ps", \
"scale_coupling_coreShell/nonuniform_scaling_trial6/results/gaussian800fs_Spectrum_QM_E17_T=0_aIndex8_bIndex3.dat" u ($1*1000):2 w lp t "Lin Trial6BestFit, CdSe/CdS, w/ broadening of 0.8ps", \
"experimental_digitized/Hendrik_4K.dat" u (4150-$1):2 w l t "Hendrik manu, CdSe/CdS/ZnS, 4K. Reversed and shifted", \
"experimental_digitized/Bawendi_PRB76_2007.dat" u (4416-$1):2 w l t "Bawendi PRB76-2007, CdSe/ZnS, 5K", \
"experimental_digitized/Klimov_PRL93_2004_Dot1.dat" u (2010+$1):2 w l t "Klimov PRL93-2004, CdSe/ZnS Dot#1, 10K", \
"experimental_digitized/Klimov_PRL93_2004_Dot2.dat" u (2010+$1):2 w l t "Klimov PRL93-2004, CdSe/ZnS Dot#2, 10K"


set origin 0.5,0
set size 0.5, 0.5
set xlabel "Temperatrue [K]"
set ylabel "FWHM [meV]"
set title ""
set key bottom right
plot [0:400] [40:100]"scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/WidthFreqShift_vs_T_gaussian40fs.dat" u 1:($4*1000) w lp t "Lin Trial3BestFit, CdSe/CdS, w/ broadening of 40fs", \
"scale_coupling_coreShell/nonuniform_scaling_trial6/aIndex8_bIndex3_T/results/WidthFreqShift_vs_T_gaussian_40fs.dat" u 1:($5*1000) w lp t "Lin Trial6BestFit, CdSe/CdS, w/ broadening of 40fs", \
"experimental_digitized/Meijerink_core-shellQD_inhomo.dat" w lp t "Meijerink JPCC-2020, CdSe/CdS, expt", \
"experimental_digitized/Bawendi_PRB76_2007_ensemble_24A.dat" w lp t "Bawendi PRB76-2007, CdSe/ZnS, R=2.4nm", \
"experimental_digitized/Bawendi_PRB76_2007_ensemble_41A.dat" w lp t "Bawendi PRB76-2007, CdSe/ZnS, R=4.1nm"


set origin 0.5,0.5
set size 0.5, 0.5
set xlabel "Energy [meV]"
set ylabel "Abs [arb. units]"
set title "NC ensemble results"
set key top left
plot [1900:2200] [0:1] "scale_coupling_coreShell/nonuniform_scaling_trial3/aIndex0_bIndex0_T/results/gaussian40fs_Spectrum_QM_E17_T=10_aIndex0_bIndex0.dat" u ($1*1000):2 w lp t "Lin Trial3BestFit, CdSe/CdS, w/ broadening of 40fs, 10K", \
"scale_coupling_coreShell/nonuniform_scaling_trial6/aIndex8_bIndex3_T/results/gaussian40fs_Spectrum_QM_E17_T=10_aIndex8_bIndex3.dat" u ($1*1000):2 w lp t "Lin Trial6BestFit, CdSe/CdS, w/ broadening of 40fs, 10K", \
"experimental_digitized/Klimov_PRL93_2004_ensemble.dat" u (2010+$1):2 w l t "Klimov PRL93-2004, CdSe/ZnS ensemble, 10K"


unset multiplot