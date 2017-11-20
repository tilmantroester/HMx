reset
unset multiplot

cmsy='/Users/Mead/Fonts/cmsy10.pfb'

if(print==0){set term aqua dashed; msun='sun'}
if(print==1){set term post enh col font ',10' fontfile cmsy; set output 'halo_profiles.eps'; msun='{/cmsy10 \014}'}

file(m)=sprintf('diagnostics/halo_profile_m%i.dat',m)

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
rhomax=3e1
set log y
set yrange [rhomin:rhomax]
set ylabel '4{/Symbol p} r^2 {/Symbol r}(r) / M' offset 2
set mytics 10

tits(m)=sprintf('M = 10^{%i} h^{-1} M_{'.msun.'}',m)

set multiplot layout 1,3

m1=13
m2=m1+2
do for [m=m1:m2] {

set title tits(m)

plot file(m) u 1:(4.*pi*$1*$1*$2) w l lw 3 dt 1 lc rgb 'black' ti 'CDM',\
     file(m) u 1:(4.*pi*$1*$1*$3) w l lw 3 dt 1 lc rgb 'red' ti 'Gas',\
     file(m) u 1:(4.*pi*$1*$1*$4) w l lw 3 dt 1 lc rgb 'blue' ti 'Stars',\
     file(m) u 1:(4.*pi*$1*$1*$5) w l lw 3 dt 2 lc rgb 'red' ti 'Bound gas',\
     file(m) u 1:(4.*pi*$1*$1*$6) w l lw 3 dt 3 lc rgb 'red' ti 'Free gas'

}

unset multiplot
