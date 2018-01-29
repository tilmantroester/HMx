reset
unset multiplot

cmsy='/Users/Mead/Fonts/cmsy10.pfb'

if(!exists("print")){print=0}
if(print==0){set term aqua dashed; sun='sun'}
if(print==1){set term post enh col fontfile cmsy; set output 'halo_pressure.eps'; sun='{/cmsy10 \014}'}

#Use a log r axis or not
ilog=0

#Constants
#JpeV=1.602e-19 # 1 = 1.602e-19 J/eV
#mpcm=100. #1 = 100 cm/m
#JpeV=1.
#mpcm=1.

if(ilog==0) {rmin=0.; rmax=3.; set xrange [rmin:rmax]}
if(ilog==1) {rmin=1e-3; rmax=1e1; set log x; set xrange [rmin:rmax]}
#if(ih==0)   {set xlabel 'r / Mpc'}
set xlabel 'r / h^{-1} Mpc'

if(ilog==0) {set yrange [*:*]}#{set yrange [0:3.5]}
if(ilog==1) {set log y; set yrange [1e-4:1e1]; set format y '10^{%T}'; set mytics 10}
set ylabel '4{/Symbol p}r^2 P_e(r,M) / (eV cm^{-3} h^{-2} Mpc^2)'

tits(m)=sprintf('M = 10^{%i} h^{-1} M_{'.sun.'}',m)
file(m)=sprintf('diagnostics/halo_profile_m%i.dat',m)
UPP(m)=sprintf('diagnostics/UPP/halo_profile_m%i.dat',m)

set multiplot layout 1,3

m1=13
m2=m1+2
do for [m=m1:m2] {

set title tits(m)

plot file(m) u 1:(4.*pi*$1*$1*$7) w l lw 3 dt 1 lc rgb 'red' ti 'Halo model',\
     UPP(m)  u 1:(4.*pi*$1*$1*$7) w l lw 3 dt 2 lc rgb 'black' ti 'UPP'

#if(ih==0 && i==1){tits=}
#if(ih==0 && i==2){tits='M = 10^{15} M_'.msun.''}

#if(ih==1 && i==1){tits='M = 10^{14}  M_'.msun.'; z = 0'}
#if(ih==1 && i==2){tits='M = 10^{15} h^{-1} M_'.msun.'; z = 0'}

#if(i==1 && icomp==0){file=file2}
#if(i==2 && icomp==0){file=file3}

#if(i==1 && icomp==1){filea=file2a; fileb=file2b}
#if(i==2 && icomp==1){filea=file3a; fileb=file3b}

#if(i==1 && icomp==2){filea=file2b; fileb=file2}
#if(i==2 && icomp==2){filea=file3b; fileb=file3}



#if(icomp==0){
#plot file u ($1):(($1)*($1)*$7) w l lw 3 dt 1 lc rgb 'black' noti
#}

#if(icomp==1){
#plot filea u ($1):(($1)*($1)*$7) w l lw 3 dt 1 lc rgb 'black' ti 'gas',\
#     fileb u ($1):(($1)*($1)*$7) w l lw 3 dt 1 lc rgb 'red'   ti 'UPP'
#}

#if(icomp==2){
#plot filea u ($1):(($1)*($1)*$7) w l lw 3 dt 1 lc rgb 'black' ti 'UPP',\
#     fileb u ($1):(($1)*($1)*$7) w l lw 2 dt 1 lc rgb 'red'   ti 'Halo model'
#}

}

unset multiplot
