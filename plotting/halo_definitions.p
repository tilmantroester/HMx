unset multiplot
reset

#Initial white space
print ''

print 'Set variable *thing* to change what is plotted (radius, mass, concentration)'
if(!exists('thing')){thing='radius'}
print 'thing: ', thing

#Set the redshift
z=0.

#Radius stuff
if(thing eq 'radius') {set title 'Radius'; file='diagnostics/radius_z0.0.dat'}
if(thing eq 'radius') {vmin=5e-3; vmax=1e1; xlab='R_v / (h^{-1} Mpc)'; ylab='R_i / (h^{-1} Mpc)'; ylab2='R_i / R_v'; rmin=0.5; rmax=1.3}
if(thing eq 'radius') {outfile='halo_radius.eps'; set format x '10^{%T}'; set format y '10^{%T}'}

#Mass stuff
if(thing eq 'mass') {set title 'Mass'; file='diagnostics/mass_z0.0.dat'}
if(thing eq 'mass') {vmin=1e8; vmax=1e17; xlab='M_v / (h^{-1} M_{sun})'; ylab='M_i / (h^{-1} M_{sun})'; ylab2='M_i / M_v'; rmin=0.5; rmax=1.2}
if(thing eq 'mass') {outfile='halo_mass.eps'; set format x '10^{%T}'; set format y '10^{%T}'}

#Concentration stuff
if(thing eq 'concentration') {set title 'Concentration'; file='diagnostics/concentration_z0.0.dat'}
if(thing eq 'concentration') {vmin=2; vmax=50; xlab='c_v'; ylab='c_i'; ylab2='c_i / c_v'; rmin=0.5; rmax=1.3}
if(thing eq 'concentration') {outfile='halo_concentration.eps'; set format x; set format y}

if(!exists('print')) {print=0}
if(print==0) {set term aqua}
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

#Final white space
print ''

unset multiplot

