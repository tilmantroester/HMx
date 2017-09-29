reset

if(print==0) set term aqua dashed
if(print==1) set term post enh col; set output 'nz.eps'

nz='lensing/nz.dat'

set xlabel 'z'
set xrange [0:5]

set ylabel 'n(z)'
set yrange [0:1.5]

plot nz u 1:2 w l lw 3 noti
