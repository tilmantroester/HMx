reset
unset multiplot

cmsy='/Users/Mead/Fonts/cmsy10.pfb'

if(print==0){set term aqua dashed; sun='sun'}
if(print==1){set term post enh col fontfile cmsy; set output 'halo_windows.eps', sun='{/cmsy10 \014}'}

file(m)=sprintf('diagnostics/halo_window_m%i.dat',m)

kmin=1e-1
kmax=1e2
set xrange [kmin:kmax]
set xlabel 'kr_v'
set log x

wmin=3.16e-3
wmax=3.16e0
set yrange [wmin:wmax]
set ylabel 'W(k,M)'
set log y

set multiplot layout 1,3

tits(m)=sprintf('M = 10^{%i} h^{-1} M_{'.sun.'}; z = 0',m)

m1=13
m2=m1+2
do for [m=m1:m2] {

set title tits(m)

plot file(m) u 1:2 w l dt 1 lw 3 lc rgb 'black' ti 'CDM',\
     file(m) u 1:3 w l dt 1 lw 3 lc rgb 'red' ti 'Gas',\
     file(m) u 1:4 w l dt 1 lw 3 lc rgb 'blue' ti 'Stars',\
     file(m) u 1:5 w l dt 2 lw 3 lc rgb 'red' ti 'Bound gas',\
     file(m) u 1:6 w l dt 3 lw 3 lc rgb 'red' ti 'Free gas'

}

unset multiplot
