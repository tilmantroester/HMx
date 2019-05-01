unset multiplot
reset

if(!exists('print')){print=0}
if(print==0) {set term aqua dashed}
if(print==1) {set term post enh col font ',12'; set output 'emulator.eps'}

# Initial white space
print ''

if(!exists('imode')) {imode=1}
print 'imode = 1: No error blob'
print 'imode = 2: Error blob'
print 'imode = ', imode
print ''

file(n,z)=sprintf('data/cosmo%d_z%d.dat',n,z)

# Number of redshifts
nz=4

if(!exists('ncos')){ncos=10}
print 'Plotting this number of cosmolgies: ncos:', ncos

kmin=2e-3
kmax=9e0
set log x
set xlabel ''
set format x ''
set xrange [kmin:kmax]

dy=0.12
ddy=0.05
set yrange [1.-dy:1+dy]

# Graph positions - x
x1=0.1
x2=0.98
dx=(x2-x1)/2.

# Graph positions - y
y2=0.96
y1=0.1
dy=(y2-y1)/nz

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

set multiplot layout 4,2

do for [i=1:2]{

if(i==1){c=3; set lmargin at screen x1+0.*dx; set rmargin at screen x1+1.*dx; set format y;    set ylabel 'P(k) / P_{emu}(k)'}
if(i==2){c=2; set lmargin at screen x1+1.*dx; set rmargin at screen x1+2.*dx; set format y ''; set ylabel ''}

do for [j=1:nz]{

unset title
if(i==1 && j==1){set title 'HMcode' offset 0,-0.8}
if(i==2 && j==1){set title 'HALOFIT' offset 0,-0.8}

if(j==1){set tmargin at screen y2-0.*dy; set bmargin at screen y2-1.*dy; set label 'z = 0.0' at graph xlab,ylab}
if(j==2){set tmargin at screen y2-1.*dy; set bmargin at screen y2-2.*dy; set label 'z = 0.5' at graph xlab,ylab}
if(j==3){set tmargin at screen y2-2.*dy; set bmargin at screen y2-3.*dy; set label 'z = 1.0' at graph xlab,ylab}
if(j==4){set tmargin at screen y2-3.*dy; set bmargin at screen y2-4.*dy; set label 'z = 2.0' at graph xlab,ylab}

set xlabel ''; set format x ''
if(j==4){set xlabel 'k / h Mpc^{-1}'; set format x}

if(imode==1){
plot 1 w l lt -1 noti,\
     1.+ddy w l lc -1 dt 2 noti,\
     1.-ddy w l lc -1 dt 2 noti,\
     for [icos=0:ncos] file(icos,j) u 1:(column(c)/$4) w l noti
}

if(imode==2) {
if(i==1) {A=0.04; b=0.1; n=1}
if(i==2) {A=0.08; b=0.1; n=1}
plot 1 w l lt -1 noti,\
     1.+ddy w l lc -1 dt 2 noti,\
     1.-ddy w l lc -1 dt 2 noti,\
     for [icos=0:ncos] file(icos,j) u 1:(column(c)/$4) w l noti,\
     error_positive(x,A,b,knl[j],n) w l lw 3 lc -1 noti,\
     error_negative(x,A,b,knl[j],n) w l lw 3 lc -1 noti
}

unset label

}

}

unset multiplot

show output
