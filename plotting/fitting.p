unset multiplot
reset

data='data/fitting.dat'

set log x

set lmargin 10
set rmargin 2

set multiplot layout 2,1

set key top left

set xlabel ''
set format x ''

set log y
set ylabel '{/Symbol D}^2(k)'
set format y '10^{%T}'

plot data u 1:3 w p lc -1 ti 'BAHAMAS',\
     data u 1:2 w l lc 1 lw 2 ti 'HMx'

set xlabel 'k / h^{-1} Mpc'
set format x

unset log y
set ylabel 'P(k) / P_{data}(k)'
set format y

plot 1 w l lt -1 noti,\
     data u 1:($2/$3) w l lc 1 lw 2 noti

unset multiplot
