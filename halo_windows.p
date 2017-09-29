reset
unset multiplot

cmsy='/Users/Mead/Fonts/cmsy10.pfb'

if(print==0) set term aqua dashed
if(print==1) set term post enh col fontfile cmsy; set output 'halo_windows.eps'

file1='diagnostics/halo_window_m13.dat'
file2='diagnostics/halo_window_m14.dat'

msun='{/cmsy10 \014}'

set xrange [1e-1:1e2]
set xlabel 'kr_v'
set log x

set yrange [3.16e-3:3.16e0]
set ylabel 'W(k,M)'
set log y

set multiplot layout 1,2

do for [i=1:2] {

if(i==1){file=file1; tits='M = 10^{13} h^{-1} M_'.msun.'; z = 0'}
if(i==2){file=file2; tits='M = 10^{14} h^{-1} M_'.msun.'; z = 0'}

set title tits

plot file u 1:2 w l dt 1 lw 3 lc rgb 'black' ti 'CDM',\
     file u 1:3 w l dt 1 lw 3 lc rgb 'red' ti 'Gas',\
     file u 1:4 w l dt 1 lw 3 lc rgb 'blue' ti 'Stars',\
     file u 1:5 w l dt 2 lw 3 lc rgb 'red' ti 'Bound gas',\
     file u 1:6 w l dt 3 lw 3 lc rgb 'red' ti 'Free gas'

}

unset multiplot
