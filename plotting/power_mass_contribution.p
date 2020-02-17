unset multiplot
reset

cmsy='/Users/Mead/Fonts/cmsy10.pfb'

if(!exists('print')){print=0}
if(print==0){set term qt dashed size 1000,500; sun='sun'}
if(print==1){set term post enh col fontfile cmsy ',10' size 10,4; sun='{/cmsy10 \014}'}

# Initial white space
print ''

if(!exists('iplot')) {iplot=1}
print 'iplot = 1: Cumulative mass contraibution'
print 'iplot = 2: PAPER: Binned mass contribution'
print 'iplot = ', iplot
print ''

if(print==1 && iplot==1) {set output 'power_mass_contribution_cumulative.eps'}
if(print==1 && iplot==2) {set output 'paper/power_mass_contribution.eps'}

# data file
if(iplot==1) {power(f1,f2,m)=sprintf('data/power_%d%d_m%d.dat',f1,f2,m)}
if(iplot==2) {power(f1,f2,m1,m2)=sprintf('data/power_%d%d_m%d_m%d.dat',f1,f2,m1,m2); base(f1,f2)=sprintf('data/power_%d%d.dat',f1,f2)}

# Field integers
imatter=2
ipressure=8

# Mass range
if(iplot==1) {m1=10; m2=16}
if(iplot==2) {m1=10; m2=15}

# wavenumber axis
kmin=1e-2
kmax=1e2
set log x
set xrange [kmin:kmax]
set xlabel 'k / h Mpc^{-1}'

set multiplot layout 1,2

# power axis
dmin=1e-3
dmax=1e3
set log y
set yrange [dmin:dmax]
set ylabel '{/Symbol D}^2_{mm}(k)'
set format y '10^{%T}'

set palette defined ( 0 "light-blue", 1 "blue", 2 "black" )
set cblabel 'log_{10} (M / h^{-1} M_{'.sun.'})'
if(iplot==2) {set cbrange [m1:m2+1]}

set label 'matter-matter' at graph 0.1,0.93

if(iplot==1){
plot for [i=m1:m2] power(imatter,imatter,i) u 1:5:(i) w l lw 3 dt 1 lc palette noti#,\
     for [i=m1:m2] power(imatter,imatter,i) u 1:3:(i) w l lw 3 dt 2 lc palette noti,\
     for [i=m1:m2] power(imatter,imatter,i) u 1:4:(i) w l lw 3 dt 3 lc palette noti
}
if(iplot==2){
plot base(imatter,imatter) u 1:5 w l lw 5 dt 1 lc -1 noti,\
     for [i=m1:m2] power(imatter,imatter,i,i+1) u 1:5:(i) w l lw 3 dt 1 lc palette noti
}

unset label

# power axis
dmin=1e-7
dmax=1e-1
set log y
set yrange [dmin:dmax]
set ylabel '{/Symbol D}^2_{mp}(k) / eV cm^{-3}'
set format y '10^{%T}'

set palette defined ( 0 "pink", 1 "red", 2 "black" )
set cblabel 'log_{10} (M / h^{-1} M_{'.sun.'})'

set label 'matter-electron pressure' at graph 0.1,0.93

if(iplot==1){
plot for [i=m1:m2] power(imatter,ipressure,i) u 1:5:(i) w l lw 3 dt 1 lc palette noti#,\
     for [i=m1:m2] power(imatter,ipressure,i) u 1:3:(i) w l lw 3 dt 2 lc palette noti,\
     for [i=m1:m2] power(imatter,ipressure,i) u 1:4:(i) w l lw 3 dt 3 lc palette noti
}
if(iplot==2){
plot base(imatter,ipressure) u 1:5 w l lw 5 dt 1 lc -1 noti,\
     for [i=m1:m2] power(imatter,ipressure,i,i+1) u 1:5:(i) w l lw 3 dt 1 lc palette noti
}

unset label

unset multiplot

show output
