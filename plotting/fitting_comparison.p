unset multiplot
reset

if(!exists('print')) {print=0}
if(print==0) {set term aqua dashed dl 1}
if(print==1) {set term post enh col font ',9'; set output 'fitting_comparison.eps'}

# File
power(mod,z,n,m,f1,f2)=sprintf('fitting/%s_%s_n%d_m%d_best_cos1_%d%d_z1.dat',mod,z,n,m,f1,f2)

# Sets of names
mods="'AGN_7p6_nu0' 'AGN_TUNED_nu0' 'AGN_TUNED_nu0_v2' 'AGN_TUNED_nu0_v3' 'AGN_TUNED_nu0_v4' 'AGN_8p0_nu0'"
names="'AGN 7p6' 'AGN tuned v1' 'AGN tuned v2' 'AGN tuned v3' 'AGN tuned v4' 'AGN 8p0'"
cols="'orange' 'blue' 'light-blue' 'light-blue' 'light-blue' 'red'"
fields="'matter' 'CDM' 'gas' 'stars' 'electon pressure'"

print ''

# Number of points in chain
n=50000

# Fitting mode
m=16

# Fields
f1=4
f2=4

print 'field 1: ', f1
print 'field 2: ', f2
print ''

# Location of zlabel on plot
labx=0.8
laby=0.1

set log x
set xlabel 'k / h Mpc^{-1}'

#set log y
set ylabel '{/Symbol D}^2(k) / (k / h Mpc^{-1})^{1.5}'
#set format y '10^{%T}'

set key top left

set multiplot layout 2,2 title ''.word(fields,f1).'-'.word(fields,f2).''

set lmargin 10
set rmargin 4
set tmargin 2
set bmargin 3

do for [i=1:4] {

if(i==1) {z='z0.0'; zlab='z = 0.0'}
if(i==2) {z='z0.5'; zlab='z = 0.5'}
if(i==3) {z='z1.0'; zlab='z = 1.0'}
if(i==4) {z='z2.0'; zlab='z = 2.0'}

set label zlab at graph labx,laby

unset key
if(i==1) {set key top left}

plot for [i=1:words(mods)] NaN w l lw 2 lc rgb word(cols,i) ti word(names,i),\
     for [i=1:words(mods)] power(word(mods,i),z,n,m,f1,f2) u 1:($2/$1**1.5) w p pt 7 ps .5 lc rgb word(cols,i) noti,\
     for [i=1:words(mods)] power(word(mods,i),z,n,m,f1,f2) u 1:($3/$1**1.5) w l dt 1 lw 2  lc rgb word(cols,i) noti

unset label

}

unset multiplot
