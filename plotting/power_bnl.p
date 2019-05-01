reset

set log x

hmo='data/power.dat'
bnl='data/power_bnl.dat'
emu='data/power_emu.dat'

set lmargin 10
set rmargin 2

set multiplot layout 2,1

set key top left

set xlabel ''
set format x ''

set log y
set ylabel '{/Symbol D}^2(k)'
set format y '10^{%T}'
pmin=1e-6
pmax=1e3
set yrange [pmin:pmax]

plot emu u 1:2 w p lc -1 pt 7 lw 3 ti 'Emulator',\
     hmo u 1:5 w l lc 1  dt 1 lw 3 ti 'Standard',\
     hmo u 1:3 w l lc 1  dt 2 lw 3 noti,\
     hmo u 1:4 w l lc 1  dt 3 lw 3 noti,\
     bnl u 1:5 w l lc 2  dt 1 lw 3 ti 'Non-linear bias',\
     bnl u 1:3 w l lc 2  dt 2 lw 3 noti,\
     bnl u 1:4 w l lc 2  dt 3 lw 3 noti

set xlabel 'k / h Mpc^{-1}'
set format x

unset log y
set format y
set yrange [*:*]

plot 1 w l lt -1 noti,\
     '<paste '.hmo.' '.emu.'' u 1:(column(5)/column(5+2)) w l lc 1 lw 3 noti,\
     '<paste '.bnl.' '.emu.'' u 1:(column(5)/column(5+2)) w l lc 2 lw 3 noti

unset multiplot
