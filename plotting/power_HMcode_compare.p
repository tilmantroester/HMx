reset

file_HMcode = 'data/power_halomodel.dat'
file_halomodel = 'data/power_HMcode.dat'

kmin = 1e-2
kmax = 1e1
klab = 'k / h Mpc^{-1}'

pmin = 1e-3
pmax = 1e3
plab = '{/Symbol D}^2(k)'
set format y '10^{%T}'

set log x
set xlabel klab
set xrange [kmin:kmax]

set log y
set ylabel plab
set yrange [pmin:pmax]

set key top left

plot NaN w l lw 3 lc -1 dt 2 ti 'Two-halo term',\
   NaN w l lw 3 lc -1 dt 3 ti 'One-halo term',\
   file_halomodel u 1:3 w l lw 3 lc 1 dt 2 noti,\
   file_halomodel u 1:4 w l lw 3 lc 1 dt 3 noti,\
   file_halomodel u 1:5 w l lw 3 lc 1 dt 1 ti 'Halo model',\
   file_HMcode u 1:3 w l lw 3 lc 2 dt 2 noti,\
   file_HMcode u 1:4 w l lw 3 lc 2 dt 3 noti,\
   file_HMcode u 1:5 w l lw 3 lc 2 dt 1 ti 'HMcode'
