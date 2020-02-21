unset multiplot
reset

mode = 'data/power.dat'
base = 'data/power_baseline.dat'

c = 5

pmin = 1e-8
pmax = 1e4

rmin = 0.8
rmax = 1.2

set lmargin 10
set rmargin 2

set multiplot layout 2, 1

set key top left

set log x
set xlabel ''
set format x ''

set log y
set ylabel '{/Symbol D}^2(k)'
set format y '10^{%T}'
set yrange [pmin:pmax]

plot NaN w l lc -1 lw 3 ti 'Baseline',\
   NaN w l lc 1 lw 3 ti 'Comparison',\
   base u 1:2 w l lc -1 lw 3 dt 2 ti 'Linear',\
   base u 1:3 w l lc -1 lw 3 dt 3 ti '2-halo',\
   base u 1:4 w l lc -1 lw 3 dt 4 ti '1-halo',\
   base u 1:5 w l lc -1 lw 3 dt 1 ti 'Full',\
   mode u 1:2 w l lc 1 lw 2 dt 2 noti, \
   mode u 1:3 w l lc 1 lw 2 dt 3 noti, \
   mode u 1:4 w l lc 1 lw 2 dt 4 noti, \
   mode u 1:5 w l lc 1 lw 2 dt 1 noti

unset log y
set ylabel 'P(k) / P_{fid}(k)'
set format y
#set yrange [*:*]
set yrange [rmin:rmax]

set xlabel 'k / h Mpc^{-1}'
set format x

plot 1 w l lt -1 noti,\
   '<paste '.mode.' '.base.'' u 1:(column(2)/column(2+c)) w l lc 1 lw 2 dt 2 noti,\
   '<paste '.mode.' '.base.'' u 1:(column(3)/column(3+c)) w l lc 1 lw 2 dt 3 noti,\
   '<paste '.mode.' '.base.'' u 1:(column(4)/column(4+c)) w l lc 1 lw 2 dt 4 noti,\
   '<paste '.mode.' '.base.'' u 1:(column(5)/column(5+c)) w l lc 1 lw 2 dt 1 noti

unset multiplot