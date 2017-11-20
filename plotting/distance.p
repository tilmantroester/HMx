reset

file='projection/distance.dat'

zmin=0.
zmax=4.
set xlabel 'z'
set xrange [zmin:zmax]

rmin=0.
rmax=6000.
set ylabel 'r / (h^{-1} Mpc)'
set yrange [rmin:*]

set key top left

plot file u 1:2 w p pt 7 ps 1 noti,\
     file u 1:2 w l lw 3 ti 'Comoving distance',\
     file u 1:3 w l lw 3 ti 'Comoving angular distance'

