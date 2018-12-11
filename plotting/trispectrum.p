reset

if(!exists('print')) {print=0}
if(print==0) {set term aqua}

set size square

set xrange [-3:2]
set xlabel 'log_{10}(k_1 / h Mpc^{-1})'

set yrange [-3:2]
set ylabel 'log_{10}(k_2 / h Mpc^{-1})'

set cblabel 'T_{1H}(k_1,k_2)'

set view map
#set pm3d

splot 'data/trispectrum.dat' u (log10($1)):(log10($2)):3 w image noti
