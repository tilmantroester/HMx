reset

cmsy='/Users/Mead/Fonts/cmsy10.pfb'

if(!exists('print')){print=0}
if(print==0){set term qt dashed; sun='sun'; width=3}
if(print==1){set term post enh col font ',18' fontfile cmsy; set output 'paper/halo_mass_fractions.eps'; sun='{/cmsy10 \014}'; width=5}

file='data/mass_fractions.dat'

#set size square

set log x
set xlabel 'M / h^{-1} M_{'.sun.'}'
set xrange [1e10:1e16]
set format x '10^{%T}'
set mxtics 10

set log y
set ylabel 'Halo mass fraction'
#set yrange [3.162e-3:3.162e0]
set yrange [3.162e-3:1e0]
set format y

#set key outside left box
set key top left box opaque

om_b=0.05
om_m=0.3

plot om_b/om_m  w l lw width dt 2 lc -1 ti 'Universal baryon',\
   file u 1:2 w l lw width dt 1 lc 2  ti 'CDM',\
   file u 1:3 w l lw width dt 1 lc 3  ti 'Gas',\
   file u 1:4 w l lw width dt 2 lc 3  ti 'Bound gas',\
   file u 1:5 w l lw width dt 3 lc 3  ti 'Unbound gas',\
   file u 1:6 w l lw width dt 1 lc 4  ti 'Stars',\
   file u 1:7 w l lw width dt 2 lc 4  ti 'Central stars',\
   file u 1:8 w l lw width dt 3 lc 4  ti 'Satellite stars',\
   file u 1:9 w l lw width dt 3 lc 5  ti 'Neutrinos'

show output

