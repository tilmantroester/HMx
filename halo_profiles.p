reset
unset multiplot

cmsy='/Users/Mead/Fonts/cmsy10.pfb'

if(print==0) set term aqua dashed
if(print==1) set term post enh col fontfile cmsy; set output 'halo_profiles.eps'

file1='diagnostics/halo_profile_m13.dat'
file2='diagnostics/halo_profile_m14.dat'

msun='{/cmsy10 \014}'

#This is the range in the Fedeli paper
#set xrange [1e-2:1e1]

#set xrange [0:4]
set log x
set xrange [1e-3:1e1]
set xlabel 'r / (h^{-1} Mpc)'

#This is the range in the Fedeli paper
#set yrange [3.16e-7:3.16e0]

#set yrange [1e-2:1e6]
#set ylabel '{/Symbol r}(r,M)'
#set log y
#set format y '10^{%T}'

#set yrange [0:3]
set log y
set yrange [1e-3:3e1]
set ylabel '4{/Symbol p} (r / h^{-1} Mpc)^2 {/Symbol r}(r,M) / M'
set mytics 10

set multiplot layout 1,2

do for [i=1:2] {

if(i==1){file=file1; tits='M = 10^{13} h^{-1} M_'.msun.'; z = 0'}
if(i==2){file=file2; tits='M = 10^{14} h^{-1} M_'.msun.'; z = 0'}

set title tits

#plot file u 1:2 w l lw 3 dt 1 lc rgb 'black' ti 'CDM',\
     file u 1:3 w l lw 3 dt 1 lc rgb 'red' ti 'Gas',\
     file u 1:4 w l lw 3 dt 1 lc rgb 'blue' ti 'Stars',\
     file u 1:5 w l lw 3 dt 2 lc rgb 'red' ti 'Bound gas',\
     file u 1:6 w l lw 3 dt 3 lc rgb 'red' ti 'Free gas'

plot file u 1:(4.*pi*$1*$1*$2) w l lw 3 dt 1 lc rgb 'black' ti 'CDM',\
     file u 1:(4.*pi*$1*$1*$3) w l lw 3 dt 1 lc rgb 'red' ti 'Gas',\
     file u 1:(4.*pi*$1*$1*$4) w l lw 3 dt 1 lc rgb 'blue' ti 'Stars',\
     file u 1:(4.*pi*$1*$1*$5) w l lw 3 dt 2 lc rgb 'red' ti 'Bound gas',\
     file u 1:(4.*pi*$1*$1*$6) w l lw 3 dt 3 lc rgb 'red' ti 'Free gas'

}

unset multiplot
