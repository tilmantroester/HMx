unset multiplot
reset

if(print==0) set term aqua

if(thing eq 'radius') {set title 'Radius'; file='definitions/radius.dat'; vmin=5e-3; vmax=1e2; xlab='R_v / (h^{-1} Mpc)'; ylab='R_i / (h^{-1} Mpc)'; ylab2='R_i / R_v'; rmin=0.5; rmax=1.3; outfile='radius.eps'; set format x '10^{%T}'; set format y '10^{%T}'}
if(thing eq 'mass') {set title 'Mass'; file='definitions/mass.dat'; vmin=1e8; vmax=1e17; xlab='M_v / (h^{-1} M_{sun})'; ylab='M_i / (h^{-1} M_{sun})'; ylab2='M_i / M_v'; rmin=0.5; rmax=1.2; outfile='mass.eps'; set format x '10^{%T}'; set format y '10^{%T}'}
if(thing eq 'concentration') {set title 'Concentration'; file='definitions/concentration.dat'; vmin=2; vmax=50; xlab='c_v'; ylab='c_i'; ylab2='c_i / c_v'; rmin=0.5; rmax=1.3; outfile='concentration.eps'; set format x; set format y}

if(print==1) set term post enh col sol; set output outfile

set lmargin 10
set rmargin 2

set log x
set xrange [vmin:vmax]
set mxtics 10

set key top left

set multiplot layout 2,1

set xlabel ''
set format x ''

set log y
set yrange [vmin:vmax]
set ylabel ylab
set mytics 10

plot file u 1:1 w l lw 3 lc 1 ti '{/Symbol D}_v',\
     file u 1:2 w l lw 3 lc 2 ti '{/Symbol D}_m = 200',\
     file u 1:3 w l lw 3 lc 3 ti '{/Symbol D}_m = 500',\
     file u 1:4 w l lw 3 lc 4 ti '{/Symbol D}_c = 200',\
     file u 1:5 w l lw 3 lc 5 ti '{/Symbol D}_c = 500'

unset title

set xlabel xlab
if(thing eq 'radius') {set format x '10^{%T}'}
if(thing eq 'mass') {set format x '10^{%T}'}
if(thing eq 'concentration') {set format x}

unset log y
set yrange [rmin:rmax]
set ylabel ylab2
set format y

plot 1 w l lt -1 noti,\
     file u 1:($2/$1) w l lw 3 lc 2 noti,\
     file u 1:($3/$1) w l lw 3 lc 3 noti,\
     file u 1:($4/$1) w l lw 3 lc 4 noti,\
     file u 1:($5/$1) w l lw 3 lc 5 noti

unset multiplot

