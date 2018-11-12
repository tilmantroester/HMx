unset multiplot
reset

if(!exists('print')) {print=0}
if(print==0) {set term aqua dashed dl 1}
if(print==1) {set term post enh col font ',9'; set output 'fitting_comparison.eps'}

# File
power_zsep(dir,mod,z,n,m,f1,f2)=sprintf('%s/%s_%s_n%d_m%d_best_cos1_%d%d_z1.dat',dir,mod,z,n,m,f1,f2)
power_zsim(dir,mod,n,m,f1,f2,iz)=sprintf('%s/%s_n%d_m%d_best_cos1_%d%d_z%d.dat',dir,mod,n,m,f1,f2,iz)

# Sets of names
mods="'AGN_7p6_nu0' 'AGN_TUNED_nu0' 'AGN_TUNED_nu0_v2' 'AGN_TUNED_nu0_v3' 'AGN_TUNED_nu0_v4' 'AGN_8p0_nu0'"
names="'AGN 7p6' 'AGN tuned v1' 'AGN tuned v2' 'AGN tuned v3' 'AGN tuned v4' 'AGN 8p0'"
cols="'orange' 'blue' 'light-blue' 'light-blue' 'light-blue' 'red'"
fields="'matter' 'CDM' 'gas' 'stars' 'electon pressure'"

print ''

if(!exists('iplot')){iplot=1}
print 'iplot = 1: Fitting to z separately'
print 'iplot = 2: Fitting to z simultaneously'
print 'iplot = ', iplot
print ''

if(!exists('dir')) {dir='fitting'}
print 'Directory: ', dir

# Number of points in chain
if(!exists('n')) {n=50000}
print 'Number of points in chain: n: ', n


# Fitting mode
if(!exists('m')) {m=16}
print 'Fitting mode: m: ', m

# Fields
if(!exists('f1')) {f1=1}
if(!exists('f2')) {f2=1}
print 'field 1: f1:', f1
print 'field 2: f2:', f2

if(iplot==1) {print 'Example file: ', power_zsep(dir,word(mods,2),'z0.0',n,m,f1,f2)}
if(iplot==2) {print 'Example file: ', power_zsim(dir,word(mods,2),n,m,f1,f2,1)}
print ''

# Location of zlabel on plot
labx=0.88
laby=0.1

kmin=2e-2
kmax=1e1
set log x
set xlabel 'k / h Mpc^{-1}'
set xrange [kmin:kmax]

rmin=1
rmax=25
set ylabel '{/Symbol D}^2(k) / (k / h Mpc^{-1})^{1.5}'
#set yrange [rmin:rmax]

set key top left opaque

set multiplot layout 2,2 title ''.word(fields,f1).'-'.word(fields,f2).''

set lmargin 10
set rmargin 4
set tmargin 2
set bmargin 3

# Rectangle
set style rect fc lt -1 fs solid 0.15 noborder

do for [i=1:4] {

if(i==1) {z='z0.0'; iz=1; zlab='z = 0.0'; kmax_fit=10.}
if(i==2) {z='z0.5'; iz=2; zlab='z = 0.5'; kmax_fit=4.}
if(i==3) {z='z1.0'; iz=3; zlab='z = 1.0'; kmax_fit=2.}
if(i==4) {z='z2.0'; iz=4; zlab='z = 2.0'; kmax_fit=1.}

set label zlab at graph labx,laby

unset key
if(i==1) {set key top left}

set obj rect from kmax_fit,-1 to kmax,30

if(iplot==1){
plot for [i=1:words(mods)] NaN w l lw 2 lc rgb word(cols,i) ti word(names,i),\
     for [i=1:words(mods)] power_zsep(dir,word(mods,i),z,n,m,f1,f2) u 1:($2/$1**1.5) w p pt 7 ps .5 lc rgb word(cols,i) noti,\
     for [i=1:words(mods)] power_zsep(dir,word(mods,i),z,n,m,f1,f2) u 1:($3/$1**1.5) w l dt 1 lw 2  lc rgb word(cols,i) noti
}

if(iplot==2){
plot for [i=1:words(mods)] NaN w l lw 2 lc rgb word(cols,i) ti word(names,i),\
     for [i=1:words(mods)] power_zsim(dir,word(mods,i),n,m,f1,f2,iz) u 1:($2/$1**1.5) w p pt 7 ps .5 lc rgb word(cols,i) noti,\
     for [i=1:words(mods)] power_zsim(dir,word(mods,i),n,m,f1,f2,iz) u 1:($3/$1**1.5) w l dt 1 lw 2  lc rgb word(cols,i) noti
}

unset obj

unset label

}

unset multiplot
