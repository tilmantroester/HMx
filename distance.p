reset

file='projection/distance.dat'

set xlabel 'z'

set ylabel 'r / (h^{-1} Mpc)'

set key top left

plot file u 1:2 w l lw 3 ti 'Comoving distance',\
     file u 1:3 w l lw 3 ti 'Comoving angular distance'

