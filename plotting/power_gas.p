reset

if(!exists('print')) {print = 0}
if(print == 0) {set term qt}
if(print == 1) {set term post enh col; set output 'plots/power_gas.eps'}

file = 'data/power_gas.dat'

klab = 'k / h Mpc^{-1}'
plab = '{/Symbol D}^2(k)'

set log x
set xlabel klab

set log y
set ylabel plab
set format y '10^{%T}'

set key top left

plot file u 1:2 w l lw 3 ti 'gas - gas',\
     file u 1:3 w l lw 3 ti 'bound gas - bound gas',\
     file u 1:4 w l lw 3 ti 'ejected gas - ejected gas',\
     file u 1:5 w l lw 3 ti 'bound gas - ejected gas'

show output