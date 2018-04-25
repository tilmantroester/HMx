reset

file='diagnostics/properies_z0.0.dat'

if(print==0){sun='sun'}

set log x
set xlabel 'M / h^{-1} M_{'.sun.'}'
set format x '10^{%T}'

set lmargin 10
set rmargin 2

set multiplot layout 2,2

set log y
set ylabel 'R / h^{-1} Mpc'

set key top left

plot file u 1:2 w l lw 2 ti 'Lagrangian radius',\
     file u 1:3 w l lw 2 ti 'Virial radius'

set log y
set ylabel '{/Symbol n}'

plot file u 1:4 w l lw 2 lc -1 noti

unset log y
set ylabel 'c'

plot file u 1:5 w l lw 2 lc -1 noti

set log y
set ylabel '{/Symbol s}'

plot file u 1:6 w l lw 2 lc -1 noti

unset multiplot

