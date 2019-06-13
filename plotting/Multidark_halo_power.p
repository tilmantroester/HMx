unset multiplot
reset

cmsy='/Users/Mead/Fonts/cmsy10.pfb'

# Initial white space
print ''

if(!exists('print')) {print=0}
print 'print = 0: aqua'
print 'print = 1: eps'
print 'print = 2: pdf'
if(print==0) {set term aqua dashed font ',10' size 1000,1200; sun='sun'}
if(print==1) {set term post enh col fontfile cmsy font ',8' size 10,12; sun='{/cmsy10 \014}'}
if(print==2) {set term pdfcairo font ',14' ps .25 size 10.24,10.24; set encoding utf8; sun='☉'}
print 'print = ', print
print ''

# File paths
power(sim,m,a,bin1,bin2)=sprintf('/Users/Mead/Physics/data/%s/power/M%d/halo_%1.5f_bin%d_bin%d_power.dat',sim,m,a,bin1,bin2)
hmod(a,bin1,bin2,opt)=sprintf('data/power_a%1.5f_hh_%d%d_%s.dat',a,bin1,bin2,opt)

if(print==1) {set output 'Multidark_halo_power.eps'}
if(print==2) {set output 'Multidark_halo_power.pdf'}

if(!exists('mesh')) {mesh=512}
print 'Simulation mesh size: mesh: ', mesh
print ''

if(!exists('sim')) {sim='Multidark'}
print 'Simulation: sim: ', sim
print ''

# k axis
kmin=7e-3
kmax=2e0
klab='k / h Mpc^{-1}'
set log x
set xrange [kmin:kmax]
set xlabel klab

# Displacement factor for residuals plots
disp=1.01
print 'Residual points being displaced in k by a small factor: ', disp
print ''

if(!exists('ptype')){ptype=1}
print 'ptype = 1: Delta^2(k)'
print 'ptype = 2: k^1.5 P(k)'
print 'ptype = 3: P(k)'
print 'ptype = ', ptype
print ''

# Axis stuff for power
if(ptype==1){
# Delta^2(K)
plab='{/Symbol D}^2_{uv}(k)'
pmin=1e-3
pmax=1e3
pow=0
con=1.
}
if(ptype==2){
# k^{1.5} P(k)
plab='{/Symbol D}^2_{uv}(k) / (k / h Mpc^{-1})^{3/2}'
pow=1.5
con=((2*pi)**3)/(4*pi)
pmin=1e0
pmax=1e2
}
if(ptype==3){
# P(k)
plab='P(k) / (h Mpc^{-1})^3'
pow=3
con=((2*pi)**3)/(4*pi)
pmin=1e1
pmax=3e3
}
pmin=pmin*1.1
pmax=pmax/1.1

# Colours
c1=7
c2=6

# Axis stuff for residual
rmin=0.5
rmax=1.5
rlab='P_{uv,mod}(k) / P_{uv,sim}(k)'

# Scale factor for plot
# a = 1.00109
# a = 0.68215
# a = 0.59103
# a = 0.49990
# a = 0.25690
if(!exists('a')) {a=1.00109}
print 'Scale factor: a: ', a
print ''

set key bottom right

print 'Example simulation file: ', power(sim,mesh,a,1,1)
print 'Example halo-model filed: ', hmod(a,1,1,'standard')
print ''

if(!exists('nf')) {nf=4}
print 'Number of fields: nf: ', nf
print ''

# Bins to use
array b[nf]
if(nf==3) {b[1]=0; b[2]=4; b[3]=5}
if(nf==4) {b[1]=0; b[2]=4; b[3]=5; b[4]=6}

# Types corresponding to bins
array type[nf]
if(nf==3){
type[1]='matter'
type[2]='haloes: 10^{12.5} -> 10^{13.0} M_{'.sun.'}/h'
type[3]='haloes: 10^{13.0} -> 10^{13.5} M_{'.sun.'}/h'
}
if(nf==4){
type[1]='matter'
type[2]='haloes: 10^{12.5} -> 10^{13.0} M_{'.sun.'}/h'
type[3]='haloes: 10^{13.0} -> 10^{13.5} M_{'.sun.'}/h'
type[4]='haloes: 10^{13.5} -> 10^{14.0} M_{'.sun.'}/h'
}

set multiplot layout nf+2,nf margins 0.06,0.98,0.05,0.98 spacing 0.005,0.005

do for [i=1:nf+2]{

do for [j=1:nf]{

okay=0

if(j>=i) {

b1=b[i]; type1=type[i]
b2=b[j]; type2=type[j]

print 'i = ', i, '; j = ', j, '; bin1 = ', b1, '; type 1 = ', type1
print 'i = ', i, '; j = ', j, '; bin2 = ', b2, '; type 2 = ', type2
print ''

set label 'u = '.type1.'' at graph 0.06,0.9
set label 'v = '.type2.'' at graph 0.06,0.8

# Set variable to say 'making a plot'
okay=1

# Power plots
set log y
ylab=plab
set yrange[pmin:pmax]

if(i==1 && j==1) {set key} else {unset key}

if(i==j) {
set xlabel klab; set format x
set ylabel ylab
set format y '10^{%T}'
} else {
set xlabel ''; set format x ''
set ylabel ''; set format y ''
}

# Do the plotting
plot power(sim,mesh,a,b1,b2) u 1:($2/$1**pow):($5/$1**pow) w e pt 7 ps .5 lw 2 lc -1 ti 'Multidark',\
     hmod(a,b1,b2,'standard')  u 1:($2/$1**pow) w l lw 3 lc -1 dt 1 ti 'Linear theory',\
     hmod(a,b1,b2,'standard')  u 1:($3/$1**pow) w l lw 3 lc c1 dt 2 noti 'Two-halo term',\
     hmod(a,b1,b2,'standard')  u 1:($4/$1**pow) w l lw 3 lc c1 dt 3 noti 'One-halo term',\
     hmod(a,b1,b2,'standard')  u 1:($5/$1**pow) w l lw 3 lc c1 dt 1 ti 'Standard halo model',\
     hmod(a,b1,b2,'bnl')       u 1:($3/$1**pow) w l lw 3 lc c2 dt 2 noti,\
     hmod(a,b1,b2,'bnl')       u 1:($5/$1**pow) w l lw 3 lc c2 dt 1 ti 'Non-linear bias halo model'

}

if(i>=j+2){

if(i==1+2 && j==1) {set key} else {unset key}

# Set variable to say 'making a plot'
okay=1

if(i==nf+2){
set format x
set xlabel klab
} else {
set xlabel ''
set format x ''
}

# Residual plots
unset log y
ylab=rlab
set yrange[rmin:rmax]
if(j==1){
set format y
set ylabel rlab
} else {
set format y ''
set ylabel ''
}

#b1=b[i-2]; type1=type[i-2]
#b2=b[j]; type2=type[j]

b1=b[j]; type1=type[j]
b2=b[i-2]; type2=type[i-2]

print 'i = ', i, '; j = ', j, '; bin1 = ', b1, '; type 1 = ', type1
print 'i = ', i, '; j = ', j, '; bin2 = ', b2, '; type 2 = ', type2
print ''

set label 'u = '.type1.'' at graph 0.94,0.9 right
set label 'v = '.type2.'' at graph 0.94,0.8 right

# Do the plotting
plot 1 w l lt -1 noti,\
     '<paste '.hmod(a,b1,b2,'standard').' '.power(sim,mesh,a,b1,b2).'' u ($1/disp):(column(5)/column(5+2)):(column(5+5)/column(5+2)) w e lw 2 pt 7 ps .5 lc c1 ti 'Standard halo model',\
     '<paste '.hmod(a,b1,b2,'bnl').' '.power(sim,mesh,a,b1,b2).''      u ($1*disp):(column(5)/column(5+2)):(column(5+5)/column(5+2)) w e lw 2 pt 7 ps .5 lc c2 ti 'Non-linear bias model'

}

if(okay==0){
set multiplot next
}

unset label

}
}

unset multiplot

print 'Scale factor: a: ', a

show output
unset output
