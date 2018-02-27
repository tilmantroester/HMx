reset
unset multiplot

cmsy='/Users/Mead/Fonts/cmsy10.pfb'

if(print==0){set term aqua dashed; msun='sun'}
if(print==1){set term post enh col font ',10' fontfile cmsy; set output 'halo_profiles.eps'; msun='{/cmsy10 \014}'}

file(m)=sprintf('diagnostics/halo_profile_m%i.dat',m)

#Power for y-axis 4*pi * r**pow * rho(r)
#pow==2 is good for linear x axis
#pow==3 is good for log x axis
pow=2

#This is the range in the Fedeli paper
#set xrange [1e-2:1e1]

#set xrange [0:4]
rmin=1e-2
rmax=5
set log x
set xrange [rmin:rmax]
set xlabel 'r / h^{-1} Mpc'

#set yrange [0:3]
rhomin=1e-3
if(pow==2){rhomax=10.}
if(pow==3){rhomax=1.}
set log y
set yrange [rhomin:rhomax]
set ylabel '4{/Symbol p} r^'.pow.' {/Symbol r}(r) / M    [h Mpc^{-1}]' offset 2
set mytics 10

tits(m)=sprintf('M = 10^{%i} h^{-1} M_{'.msun.'}',m)

set multiplot layout 1,3

m1=13
m2=m1+2
if(pow==2){mtit=m2}
if(pow==3){mtit=m1}
do for [m=m1:m2] {

set title tits(m)

ti_CDM=''
ti_gas=''
ti_stars=''
ti_bound=''
ti_free=''
if(m==mtit) {ti_CDM='CDM'}
if(m==mtit) {ti_gas='Gas'}
if(m==mtit) {ti_stars='Stars'}
if(m==mtit) {ti_bound='Bound gas'}
if(m==mtit) {ti_free='Free gas'}

plot file(m) u 1:(4.*pi*($1**pow)*$2) w l lw 3 dt 1 lc rgb 'black' ti ti_CDM,\
     file(m) u 1:(4.*pi*($1**pow)*$3) w l lw 3 dt 1 lc rgb 'red'   ti ti_gas,\
     file(m) u 1:(4.*pi*($1**pow)*$4) w l lw 3 dt 1 lc rgb 'blue'  ti ti_stars,\
     file(m) u 1:(4.*pi*($1**pow)*$5) w l lw 3 dt 2 lc rgb 'red'   ti ti_bound,\
     file(m) u 1:(4.*pi*($1**pow)*$6) w l lw 3 dt 3 lc rgb 'red'   ti ti_free

}

unset multiplot
