unset multiplot
reset

if(!exists('print')){print=0}
if(print==0) {set term aqua dashed font ',8' dl 1}
if(print==1) {set term post enh col font ',12'; set output 'fitting_power.eps'}
if(print==2) {set term qt}

# Initial white space
print ''

# Data files
data(base,type,icosmo,if1,if2,iz)=sprintf('%s_%s_cos%d_%d%d_z%d.dat',base,type,icosmo,if1,if2,iz)
#data(base,nchain,model,chain,type,icosmo,if1,if2,iz)=sprintf('%s_n%d_m%d_c%d_%s_cos%d_%d%d_z%d.dat',base,nchain,model,chain,type,icosmo,if1,if2,iz)

# Field types
field_names="'matter' 'CDM' 'gas' 'stars' 'pressure'"

array cols[5]
cols[1]=1
cols[2]=2
cols[3]=3
cols[4]=4
cols[5]=6

if(!exist('c')){c=1}
print 'Chain: ', c
print ''

klab='k / h^{-1} Mpc'

plab='{/Symbol D}^2(k)'

rlab='R_{HM}(k) / R_{sim}(k)'
rmin=0.8
rmax=1.2

labx=0.8
laby=0.9

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
print 'iplot = 1: Power and residual'
print 'iplot = 2: Residual only'
#print 'iplot = 3: Many residuals M15 50,000'
#print 'iplot = 4: Many residuals M11 50,000'
#print 'iplot = 5: Many residuals M16 50,000'
#print 'iplot = 6: Many residuals M17 50,000'
#print 'iplot = 7: Many residuals M19 50,000'
print 'iplot = 3: Separate z panels'
print 'iplot = ', iplot
print ''

if(!exists('base')){base='fitting/test'}
if(iplot==1 || iplot==2){
print 'Data file base: base: ', base
print 'Example data file: ', data(base,'best',1,1,1,1)
print ''
}

#off=1.5
#kmin=0.01/off
#kmax=10.*off
#set xrange [kmin:kmax]
set xrange [*:*]
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

plot for [i=1:ncos] for [j=1:nz] for [j1=1:nf] for [j2=j1:nf] data(base,'best',i,j1,j2,j) u 1:3 w p lc j1 dt 1 lw 2 noti 'BAHAMAS',\
     for [i=1:ncos] for [j=1:nz] for [j1=1:nf] for [j2=j1:nf] data(base,'orig',i,j1,j2,j) u 1:2 w l lc 0  dt 1 lw 2 noti 'HMx',\
     for [i=1:ncos] for [j=1:nz] for [j1=1:nf] for [j2=j1:nf] data(base,'best',i,j1,j2,j) u 1:2 w l lc j1 dt 1 lw 2 noti 'HMx',\
     for [i=1:ncos] for [j=1:nz] for [j1=1:nf] for [j2=j1:nf] data(base,'best',i,j1,j2,j) u 1:2 w l lc j2 dt 2 lw 2 noti

set xlabel klab
set format x

unset log y
set yrange [rmin:rmax]
set ylabel rlab
set format y

plot 1 w l lt -1 noti,\
     0.95 w l lc -1 dt 2 noti,\
     1.05 w l lc -1 dt 2 noti,\
     0.99 w l lc -1 dt 2 noti,\
     1.01 w l lc -1 dt 2 noti,\
     for [i=1:ncos] for [j=1:nz] for [j1=1:nf] for [j2=j1:nf] data(base,'orig',i,j1,j2,j) u 1:($2/$3) w l lc 0  dt j lw 2 noti,\
     for [i=1:ncos] for [j=1:nz] for [j1=1:nf] for [j2=j1:nf] data(base,'best',i,j1,j2,j) u 1:($2/$3) w l lc j1 dt 1 lw 2 noti,\
     for [i=1:ncos] for [j=1:nz] for [j1=1:nf] for [j2=j1:nf] data(base,'best',i,j1,j2,j) u 1:($2/$3) w l lc j2 dt 2 lw 2 noti

unset multiplot

}

if(iplot==2){

unset log y
set yrange [rmin:rmax]
if(nz==1 || nz==4) {set ylabel rlab}
if(nz==11) {set ylabel ''}
set format y

set key top left

print 'base: ', base
print ''

if(nz != 1 && nz != 4 && nz !=11) {print 'nz must equal either 1, 4 or 11'; print ''; exit}

set style rect fc lt -1 fs solid 0.25 noborder

if(nz==4)  {set multiplot layout 2,2}
if(nz==11) {set multiplot layout 4,3}

do for [j=1:nz]{

unset key
if(j==1){set key top left}

if(nz==1 || nz==4) {
if(j==1) {kmin=10.; zlab='z = 0.0'}
if(j==2) {kmin=4. ; zlab='z = 0.5'}
if(j==3) {kmin=2. ; zlab='z = 1.0'}
if(j==4) {kmin=1. ; zlab='z = 2.0'}
}

if(nz==11)  {
kmin=100.
if(j==1)  {zlab='z = 0.000'}
if(j==2)  {zlab='z = 0.125'}
if(j==3)  {zlab='z = 0.250'}
if(j==3)  {zlab='z = 0.375'}
if(j==5)  {zlab='z = 0.500'}
if(j==6)  {zlab='z = 0.750'}
if(j==7)  {zlab='z = 1.000'}
if(j==8)  {zlab='z = 1.250'}
if(j==9)  {zlab='z = 1.500'}
if(j==10) {zlab='z = 1.750'}
if(j==11) {zlab='z = 2.000'}
}

#set style fill transparent solid 0.5 noborder
set obj rect from kmin,rmin to 20.,rmax back
set label zlab at graph labx,laby front


plot 1 w l lt -1 noti,\
     0.95 w l lc -1 dt 2 noti,\
     1.05 w l lc -1 dt 2 noti,\
     0.99 w l lc -1 dt 2 noti,\
     1.01 w l lc -1 dt 2 noti,\
     for [i=1:ncos] for [j1=1:nf] for [j2=j1:nf] data(base,'orig',i,j1,j2,j) u 1:($2/$3) w l lc 0        dt 1 lw 2 noti,\
     for [i=1:ncos] for [j1=1:nf] for [j2=j1:nf] data(base,'best',i,j1,j2,j) u 1:($2/$3) w l lc cols[j1] dt 1 lw 2 noti,\
     for [i=1:ncos] for [j1=1:nf] for [j2=j1:nf] data(base,'best',i,j1,j2,j) u 1:($2/$3) w l lc cols[j2] dt 2 lw 2 noti

#for [j1=1:nf] NaN w l lw 2 lc cols[j1] ti word(field_names,j1),\

unset obj
unset label

}

if(nz==4 || nz==11) {unset multiplot}

}

if(iplot==3) {

unset log y
set yrange [rmin:rmax]
set ylabel rlab
set format y

set key top left

print 'base: ', base
print ''

if(nz != 11) {print 'nz must equal 11'; print ''; exit}

set style rect fc lt -1 fs solid 0.25 noborder

if(nz==4) {set multiplot layout 2,2}

do for [j=1:nz]{

unset key
if(j==1){set key top left}

if(j==1){kmin=10.; zlab='z = 0.0'}
if(j==2){kmin=4. ; zlab='z = 0.5'}
if(j==3){kmin=2. ; zlab='z = 1.0'}
if(j==4){kmin=1. ; zlab='z = 2.0'}

#set style fill transparent solid 0.5 noborder
set obj rect from kmin,rmin to 20.,rmax back
set label zlab at graph labx,laby front

plot 1 w l lt -1 noti,\
     0.95 w l lc -1 dt 2 noti,\
     1.05 w l lc -1 dt 2 noti,\
     0.99 w l lc -1 dt 2 noti,\
     1.01 w l lc -1 dt 2 noti,\
     for [j1=1:nf] NaN w l lw 2 lc cols[j1] ti word(field_names,j1),\
     for [i=1:ncos] for [j1=1:nf] for [j2=j1:nf] data(base,'orig',i,j1,j2,j) u 1:($2/$3) w l lc 0        dt 1 lw 2 noti,\
     for [i=1:ncos] for [j1=1:nf] for [j2=j1:nf] data(base,'best',i,j1,j2,j) u 1:($2/$3) w l lc cols[j1] dt 1 lw 2 noti,\
     for [i=1:ncos] for [j1=1:nf] for [j2=j1:nf] data(base,'best',i,j1,j2,j) u 1:($2/$3) w l lc cols[j2] dt 2 lw 2 noti

unset obj
unset label

}

}

#if(iplot==3 || iplot==4 || iplot==5 || iplot==6 || iplot==7){

# iplot = 3: M15 matter fitting
# iplot = 4: M11 everything fitting
# iplot = 4: M16 everything fitting (minus epressure-epressure)
#if(iplot==3 && print==1) {set output 'fitting_m15.eps'}
#if(iplot==4 && print==1) {set output 'fitting_m11.eps'}
#if(iplot==5 && print==1) {set output 'fitting_m16.eps'}
#if(iplot==6 && print==1) {set output 'fitting_m17.eps'}
#if(iplot==7 && print==1) {set output 'fitting_m19.eps'}

# models
#mods="'AGN_7p6_nu0' 'AGN_TUNED_nu0' 'AGN_8p0_nu0'"
#mod_names="'AGN 7p6' 'AGN tuned' 'AGN 8p0'"

# redshifts
#zs="'z0.0' 'z0.5' 'z1.0' 'z2.0'"
#z_names="'z = 0.0' 'z = 0.5' 'z = 1.0' 'z = 2.0'"

#if(iplot==3) {nf=4; type='m15'; n=50000}
#if(iplot==4) {nf=5; type='m11'; n=50000}
#if(iplot==5) {nf=5; type='m16'; n=50000}
#if(iplot==6) {nf=5; type='m17'; n=50000}
#if(iplot==7) {nf=1; type='m19'; n=50000}

# y axis
#dr=0.28
#rmin=1.-dr
#rmax=1.+dr
#unset log y
#set yrange [rmin:rmax]
#set ylabel rlab
#set format y

#set key top left

#set multiplot

#nm=words(mods)
#nz=words(zs)

#dy=(top-bot)/real(nz)
#dx=(rig-lef)/real(nm)

#do for [iz=1:nz]{

#set tmargin at screen top-(iz-1)*dy
#set bmargin at screen top-(iz-0)*dy

#do for [im=1:nm]{

#set lmargin at screen lef+(im-1)*dx
#set rmargin at screen lef+(im-0)*dx

#if(iz==1) {set label word(mod_names,im) at graph 0.65,0.9}
#if(iz==1 || iz==2 || iz==3){set format x ''; set xlabel ''}
#if(iz==4) {set format x; set xlabel klab}

#if(im==1){set format y; set ylabel rlab; set label word(z_names,iz) at graph 0.05,0.1}
#if(im==2 || im==3) {set format y ''; set ylabel ''}

#unset key
#if(iz==1 && im==1) {set key top left}

#if(iplot==3 || iplot==4 || iplot==5) {data(mod,z,n,type,c,best,i1,i2)=sprintf('fitting/%s_%s_n%d_%s_c%d_%s_cos1_%d%d_z1.dat',mod,z,n,type,c,best,i1,i2)}
#if(iplot==6 || iplot==7) {data(mod,n,type,c,best,i1,i2,iz)=sprintf('fitting/%s_n%d_%s_c%d_%s_cos1_%d%d_z%d.dat',mod,n,type,c,best,i1,i2,iz)}

#if(iplot==3 || iplot==4 || iplot==5){
#plot 1 w l lt -1 noti,\
#     0.95 w l lc -1 dt 2 noti,\
#     1.05 w l lc -1 dt 2 noti,\
#     for [j=1:nf] NaN w l lw 2 lc cols[j] ti word(field_names,j),\
#     for [j1=1:nf] for [j2=j1:nf] data(word(mods,im),word(zs,iz),n,type,c,'orig',j1,j2) u 1:($2/$3) w l lc rgb 'light-grey' dt 1 lw 1.5 noti,\
#     for [j1=1:nf] for [j2=j1:nf] data(word(mods,im),word(zs,iz),n,type,c,'best',j1,j2) u 1:($2/$3) w l lc cols[j1] dt 1 lw 1.5 noti,\
#     for [j1=1:nf] for [j2=j1:nf] data(word(mods,im),word(zs,iz),n,type,c,'best',j1,j2) u 1:($2/$3) w l lc cols[j2] dt 2 lw 1.5 noti
#}

#if(iplot==6 || iplot==7){
#plot 1 w l lt -1 noti,\
#     0.95 w l lc -1 dt 2 noti,\
#     1.05 w l lc -1 dt 2 noti,\
#     for [j=1:nf] NaN w l lw 2 lc cols[j] ti word(field_names,j),\
#     for [j1=1:nf] for [j2=j1:nf] data(word(mods,im),n,type,c,'orig',j1,j2,iz) u 1:($2/$3) w l lc rgb 'light-grey' dt 1 lw 1.5 noti,\
#     for [j1=1:nf] for [j2=j1:nf] data(word(mods,im),n,type,c,'best',j1,j2,iz) u 1:($2/$3) w l lc cols[j1] dt 1 lw 1.5 noti,\
#     for [j1=1:nf] for [j2=j1:nf] data(word(mods,im),n,type,c,'best',j1,j2,iz) u 1:($2/$3) w l lc cols[j2] dt 2 lw 1.5 noti
#}

#unset label

#}

#}

#unset multiplot

#}
