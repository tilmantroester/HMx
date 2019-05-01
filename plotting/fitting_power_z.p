reset
unset multiplot

if(!exists('print')) {print=0}
if(print==0) {set term aqua dashed}
if(print==1) {set term post enh col sol font ',10'}

# Initial white space
print ''

fname(base,z,n,m,best,f1,f2)=sprintf('%s_z%4.3f_n%i_m%i_c1_%s_cos1_%i%i_z1.dat',base,z,n,m,best,f1,f2)
sim(f1,f2)=sprintf('/Users/Mead/Physics/BAHAMAS/power/M1024/%s_L400N1024_WMAP9_snap32_%s_%s_power.dat',mod,f1,f2)
outfile(m,imode)=sprintf('fitting_m%d_type%d.eps',m,imode)
zlab(z)=sprintf('z = %4.3f',z)

# Set the mode
if(!exists('imode')) {imode=1}
print 'imode = 1: P_mod / P_sim'
print 'imode = 2: k^1.5*P(k)'
print 'imode = 3: P / P_mm'
print 'imode = ', imode
print ''

# Set the file base
if(!exists('base')) {base='fitting/AGN_TUNED_nu0'}
print 'File base: base: ', base

# Set the number of points in chain
if(!exists('n')) {n=3000}
print 'Number in chain: n: ', n

# Set the mode
if(!exists('m')) {m=34}
print 'Fitting mode: m: ', m
print ''

if(!exists('nf')) {nf=1}
print 'Number of fields: nf: ', nf
print ''

set log x

if(imode==1){
unset log y
set format y
dr=0.2
ddr=0.05
set yrange [1.-dr:1.+dr]
}
if(imode==2){
set log y
set format y '10^{%T}'
}
if(imode==3){
unset log y
}

if(print==1) {set output outfile(m,imode)}

set multiplot layout 3,4

do for [i=1:11]{

if(i==1)  {z=0.000}
if(i==2)  {z=0.125}
if(i==3)  {z=0.250}
if(i==4)  {z=0.375}
if(i==5)  {z=0.500}
if(i==6)  {z=0.750}
if(i==7)  {z=1.000}
if(i==8)  {z=1.250}
if(i==9)  {z=1.500}
if(i==10) {z=1.750}
if(i==11) {z=2.000}

set label zlab(z) at graph 0.1,0.9

if(imode==1){

plot 1 w l lt -1 noti,\
     1.-ddr w l lc -1 dt 2 noti,\
     1.+ddr w l lc -1 dt 2 noti,\
     for [f1=1:nf] for [f2=f1:nf] fname(base,z,n,m,'orig',f1,f2) u 1:($2/$3) w l lw 2 lc 0  dt 2 noti,\
     for [f1=1:nf] for [f2=f1:nf] fname(base,z,n,m,'best',f1,f2) u 1:($2/$3) w l lw 2 lc f1 dt 1 noti,\
     for [f1=1:nf] for [f2=f1:nf] fname(base,z,n,m,'best',f1,f2) u 1:($2/$3) w l lw 2 lc f2 dt 2 noti

}

if(imode==2){

plot for [f1=1:nf] for [f2=f1:nf] fname(base,z,n,m,'orig',f1,f2) u 1:($2/$1**1.5) w l lw 2 lc 0  dt 2 noti,\
     for [f1=1:nf] for [f2=f1:nf] fname(base,z,n,m,'best',f1,f2) u 1:($3/$1**1.5) w l lw 2 lc -1 dt 1 noti,\
     for [f1=1:nf] for [f2=f1:nf] fname(base,z,n,m,'best',f1,f2) u 1:($2/$1**1.5) w l lw 2 lc f1 dt 1 noti,\
     for [f1=1:nf] for [f2=f1:nf] fname(base,z,n,m,'best',f1,f2) u 1:($2/$1**1.5) w l lw 2 lc f2 dt 2 noti

}

if(imode==3){

plot for [f1=1:nf] for [f2=f1:nf] fname(base,z,n,m,'orig',f1,f2) u 1:($2/$4) w l lw 2 lc 0  dt 2 noti,\
     for [f1=1:nf] for [f2=f1:nf] fname(base,z,n,m,'best',f1,f2) u 1:($3/$4) w l lw 2 lc -1 dt 1 noti,\
     for [f1=1:nf] for [f2=f1:nf] fname(base,z,n,m,'best',f1,f2) u 1:($2/$4) w l lw 2 lc f1 dt 1 noti,\
     for [f1=1:nf] for [f2=f1:nf] fname(base,z,n,m,'best',f1,f2) u 1:($2/$4) w l lw 2 lc f2 dt 2 noti

}

unset label

}

unset multiplot

show output
