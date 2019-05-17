unset multiplot
reset

if(!exists('print')){print=0}
if(print==0) {set term aqua dashed}
if(print==1) {set term post enh col font ',10'}

# Initial white space
print ''

# Choose plotting mode
if(!exists('imode')) {imode=1}
print 'imode = 1: Many cosmologies; no error blob'
print 'imode = 2: Many cosmologies; error blob'
print 'imode = 3: Single cosmology'
print 'imode = ', imode
print ''

# File name
file(n,z)=sprintf('data/cosmo%d_z%d.dat',n,z)

# Columns
c_li=2 # Linear
c_em=3 # Emulator
c_ql=4 # Quasi-linear halofit
c_oh=5 # Non-linear halofit
c_hf=6 # Halofit
c_2h=7 # Two-halo term
c_1h=8 # One-halo term
c_hm=9 # Halo model

# Label for k axis
klab='k / h Mpc^{-1}'

if(imode==1 || imode==2){

if(print==1) {set output 'power_emulator_many.eps'}

# k range
kmin=2e-3
kmax=9e0
set log x
set xrange [kmin:kmax]

# Number of redshifts
if(!exists('nz')) {nz=4}
print 'Number of redshifts: nz: ', nz
print ''

# Number of cosmologies
if(!exists('ncos')){ncos=10}
print 'Plotting this number of cosmolgies: ncos: ', ncos
print ''

# Residual range
dy=0.34
ddy=0.05
set yrange [1.-dy:1+dy]

# Label positions
xlab=0.05
ylab=0.90

# Final white space
print ''

# Error blob function
error_positive(k,A,b,knl,n)=1.+A*(1./(1.+(k/(b*knl))**(-2*n)))
error_negative(k,A,b,knl,n)=2.-error_positive(k,A,b,knl,n)

# Array of non-linear wave numbers at z = 0, 0.5, 1, 2
# These were calculated from the 'boring' LCDM cosmology
array knl[4]
knl[1]=0.400 # z=0.0
knl[2]=0.690 # z=0.5
knl[3]=1.206 # z=1.0
knl[4]=3.446 # z=2.0

set xlabel ''
set format x ''

set multiplot layout nz,2 margins 0.10,0.98,0.11,0.94 spacing 0.0,0.00

do for [j=1:nz]{

do for [i=1:2]{

if(i==1){c=c_hm; set format y;    set ylabel 'P(k) / P_{emu}(k)'}
if(i==2){c=c_hf; set format y ''; set ylabel ''}

unset title
if(i==1 && j==1){set title 'HMcode' offset 0,-0.8}
if(i==2 && j==1){set title 'HALOFIT' offset 0,-0.8}

if(j==1){set label 'z = 0.0' at graph xlab,ylab}
if(j==2){set label 'z = 0.5' at graph xlab,ylab}
if(j==3){set label 'z = 1.0' at graph xlab,ylab}
if(j==4){set label 'z = 2.0' at graph xlab,ylab}

set xlabel ''; set format x ''
if(j==nz){set xlabel klab; set format x}

if(imode==1){
plot 1 w l lt -1 noti,\
     1.+ddy w l lc -1 dt 2 noti,\
     1.-ddy w l lc -1 dt 2 noti,\
     for [icos=0:ncos] file(icos,j) u 1:(column(c)/column(c_em)) w l noti
}

if(imode==2) {
if(i==1) {A=0.04; b=0.1; n=1}
if(i==2) {A=0.08; b=0.1; n=1}
plot 1 w l lt -1 noti,\
     1.+ddy w l lc -1 dt 2 noti,\
     1.-ddy w l lc -1 dt 2 noti,\
     for [icos=0:ncos] file(icos,j) u 1:(column(c)/column(c_em)) w l noti,\
     error_positive(x,A,b,knl[j],n) w l lw 3 lc -1 noti,\
     error_negative(x,A,b,knl[j],n) w l lw 3 lc -1 noti
}

unset label

}

}

unset multiplot

}

if(imode==3){

if(print==1) {set output 'power_emulator.eps'}

# Set the cosmology
icos=0

# Set the redshift
# 1 - z = 0.0
# 2 - z = 0.5
# 3 - z = 1.0
# 4 - z = 2.0
iz=2

col_li=6
col_hm=7

# k range
kmin=0.01
kmax=10.
set log x
set xrange [kmin:kmax]
set format x ''

set multiplot layout 2,1 margins 0.08,0.97,0.13,0.97 spacing 0.01,0.03

# Power range
pmin=1e-3
pmax=2e2
set log y
set ylabel '{/Symbol D}^2(k)'
set yrange [pmin:pmax]
set format y '10^{%T}'

if(iz==1) {zlab='z = 0.0'}
if(iz==2) {zlab='z = 0.5'}
if(iz==3) {zlab='z = 1.0'}
if(iz==4) {zlab='z = 2.0'}

set label zlab at graph 0.9,0.1

set key top left

plot file(icos,iz) u 1:(column(c_em)) w p lw 3 lc -1     pt 1 ti 'Emulator',\
     file(icos,iz) u 1:(column(c_li)) w l lw 3 lc col_li dt 1 ti 'Linear',\
     file(icos,iz) u 1:(column(c_2h)) w l lw 3 lc col_hm dt 2 ti 'Two-halo term',\
     file(icos,iz) u 1:(column(c_1h)) w l lw 3 lc col_hm dt 3 ti 'One-halo term',\
     file(icos,iz) u 1:(column(c_hm)) w l lw 3 lc col_hm dt 1 ti 'Halo model'

unset label

set xlabel klab
set format x

rmin=0.5
rmax=1.2
unset log y
set format y
set yrange [rmin:rmax]
set ylabel 'P_{mod}(k) / P_{emu}(k)'
dy=0.1

plot 1 w l lt -1 noti,\
     1.-dy w l lc -1 dt 2 noti,\
     1.+dy w l lc -1 dt 2 noti,\
     '<paste '.file(icos,iz).' '.file(icos,iz).'' u 1:(column(c_li)/column(c_em)) w l lw 3 lc col_li dt 1 noti,\
     '<paste '.file(icos,iz).' '.file(icos,iz).'' u 1:(column(c_2h)/column(c_em)) w l lw 3 lc col_hm dt 2 noti,\
     '<paste '.file(icos,iz).' '.file(icos,iz).'' u 1:(column(c_1h)/column(c_em)) w l lw 3 lc col_hm dt 3 noti,\
     '<paste '.file(icos,iz).' '.file(icos,iz).'' u 1:(column(c_hm)/column(c_em)) w l lw 3 lc col_hm dt 1 noti

unset multiplot

}

show output
unset output
