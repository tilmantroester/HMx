unset multiplot
reset

if(!exists('print')){print=0}
if(print==0) {set term aqua dashed dl 1}
if(print==1) {set term post enh col}
if(print==2) {set term qt}

# Initial white space
print ''

# Data files
data(base,type,icosmo,if1,if2,iz)=sprintf('%s_%s_cos%d_%d%d_z%d.dat',base,type,icosmo,if1,if2,iz)

# Field types
field_names="'matter' 'CDM' 'gas' 'stars' 'electron pressure'"

array cols[5]
cols[1]=1
cols[2]=2
cols[3]=3
cols[4]=4
cols[5]=6

klab='k / h^{-1} Mpc'
plab='{/Symbol D}^2(k)'
rlab='R_{HM}(k) / R_{sim}(k)'

# Set number of cosmoloies
if(!exists('ncos')){ncos=1}
print 'Number of cosmologies: ncos: ', ncos

if(!exists('nf')){nf=1}
print 'Number of fields: nf: ', nf

# Set number of redshifts
if(!exists('nz')){nz=1}
print 'Number of redshifts: nz: ', nz
print ''

# Set what to plot
if(!exists('iplot')){iplot=1}
print 'iplot = 1 - Power and residual'
print 'iplot = 2 - Residual only'
print 'iplot = 3 - Many residuals M15 50,0000'
print 'iplot = 4 - Many residuals M11 50,0000'
print 'iplot = ', iplot
print ''

if(!exists('base')){base='fitting/test'}
print 'Data file base: base: ', base
print 'Example data file: ', data(base,'best',1,1,1,1)
print ''

off=1.5
kmin=0.1/off
kmax=10.*off
set xrange [kmin:kmax]
set log x

top=0.98
bot=0.10
lef=0.10
rig=0.98

if(iplot==1){

set lmargin 10
set rmargin 2

set multiplot layout 2,1

set key top left

set xlabel ''
set format x ''

set log y
set ylabel plab
set format y '10^{%T}'

col(i,ni,j,nj,k,nk,l,nl)=nj*nk*nl*(i-1)+nk*nl*(j-1)+nl*(k-1)+l

plot for [i=1:ncos] for [j=1:nz] for [j1=1:nf] for [j2=j1:nf] data(base,'best',i,j1,j2,j) u 1:3 w p lc col(i,ncos,j,nz,j1,nf,j2,nf) dt 1 lw 2 noti 'BAHAMAS',\
     for [i=1:ncos] for [j=1:nz] for [j1=1:nf] for [j2=j1:nf] data(base,'orig',i,j1,j2,j) u 1:2 w l lc 0                            dt j lw 2 noti 'HMx',\
     for [i=1:ncos] for [j=1:nz] for [j1=1:nf] for [j2=j1:nf] data(base,'best',i,j1,j2,j) u 1:2 w l lc col(i,ncos,j,nz,j1,nf,j2,nf) dt j lw 2 noti 'HMx'

set xlabel klab
set format x

rmin=0.8
rmax=1.2
unset log y
set yrange [rmin:rmax]
set ylabel rlab
set format y

plot 1 w l lt -1 noti,\
     0.95 w l lc -1 dt 2 noti,\
     1.05 w l lc -1 dt 2 noti,\
     0.99 w l lc -1 dt 2 noti,\
     1.01 w l lc -1 dt 2 noti,\
     for [i=1:ncos] for [j=1:nz] for [j1=1:nf] for [j2=j1:nf] data(base,'orig',i,j1,j2,j) u 1:($2/$3) w l lc 0                            dt j lw 2 noti,\
     for [i=1:ncos] for [j=1:nz] for [j1=1:nf] for [j2=j1:nf] data(base,'best',i,j1,j2,j) u 1:($2/$3) w l lc col(i,ncos,j,nz,j1,nf,j2,nf) dt j lw 2 noti

unset multiplot

}

if(iplot==2){

rmin=0.8
rmax=1.2
unset log y
set yrange [rmin:rmax]
set ylabel 'P_{HM}(k) / P_{sim}(k)'
set format y

set key top left

set title base noenhanced

plot 1 w l lt -1 noti,\
     0.95 w l lc -1 dt 2 noti,\
     1.05 w l lc -1 dt 2 noti,\
     0.99 w l lc -1 dt 2 noti,\
     1.01 w l lc -1 dt 2 noti,\
     for [j=1:nf] NaN w l lw 2 lc cols[j] ti word(field_names,j),\
     for [i=1:ncos] for [j=1:nz] for [j1=1:nf] for [j2=j1:nf] data(base,'orig',i,j1,j2,j) u 1:($2/$3) w l lc 0        dt j lw 2 noti,\
     for [i=1:ncos] for [j=1:nz] for [j1=1:nf] for [j2=j1:nf] data(base,'best',i,j1,j2,j) u 1:($2/$3) w l lc cols[j1] dt 1 lw 2 noti,\
     for [i=1:ncos] for [j=1:nz] for [j1=1:nf] for [j2=j1:nf] data(base,'best',i,j1,j2,j) u 1:($2/$3) w l lc cols[j2] dt 2 lw 2 noti

}

if(iplot==3 || iplot==4){

# iplot = 3: M15 matter fitting
# iplot = 4: M11 everything fitting
if(iplot==3 && print==1) {set output 'fitting_m15.eps'}
if(iplot==4 && print==1) {set output 'fitting_m11.eps'}

# models
mods="'AGN_7p6' 'AGN_TUNED' 'AGN_8p0'"
mod_names="'AGN 7p6' 'AGN tuned' 'AGN 8p0'"

# redshifts
zs="'z0.0' 'z0.5' 'z1.0' 'z2.0'"
z_names="'z = 0.0' 'z = 0.5' 'z = 1.0' 'z = 2.0'"

if(iplot==3) {nf=4; type='m15'}
if(iplot==4) {nf=5; type='m11'}

# y axis
dr=0.28
rmin=1.-dr
rmax=1.+dr
unset log y
set yrange [rmin:rmax]
set ylabel rlab
set format y

set key top left

set multiplot# layout 4,3

nm=words(mods)
nz=words(zs)

dy=(top-bot)/real(nz)
dx=(rig-lef)/real(nm)

do for [iz=1:nz]{

set tmargin at screen top-(iz-1)*dy
set bmargin at screen top-(iz-0)*dy

do for [im=1:nm]{

set lmargin at screen lef+(im-1)*dx
set rmargin at screen lef+(im-0)*dx

if(iz==1) {set label word(mod_names,im) at graph 0.65,0.9}
if(iz==1 || iz==2 || iz==3){set format x ''; set xlabel ''}
if(iz==4) {set format x; set xlabel klab}

if(im==1){set format y; set ylabel rlab; set label word(z_names,iz) at graph 0.05,0.1}
if(im==2 || im==3) {set format y ''; set ylabel ''}

unset key
if(iz==1 && im==1) {set key top left}

data(mod,z,type,i1,i2)=sprintf('fitting/%s_nu0_%s_n50000_%s_best_cos1_%d%d_z1.dat',mod,z,type,i1,i2)

plot 1 w l lt -1 noti,\
     0.95 w l lc -1 dt 2 noti,\
     1.05 w l lc -1 dt 2 noti,\
     for [j=1:nf] NaN w l lw 2 lc cols[j] ti word(field_names,j),\
     for [j1=1:nf] for [j2=j1:nf] data(word(mods,im),word(zs,iz),type,j1,j2) u 1:($2/$3) w l lc 0        dt j lw 1.5 noti,\
     for [j1=1:nf] for [j2=j1:nf] data(word(mods,im),word(zs,iz),type,j1,j2) u 1:($2/$3) w l lc cols[j1] dt 1 lw 1.5 noti,\
     for [j1=1:nf] for [j2=j1:nf] data(word(mods,im),word(zs,iz),type,j1,j2) u 1:($2/$3) w l lc cols[j2] dt 2 lw 1.5 noti

unset label

}

}

unset multiplot

}
