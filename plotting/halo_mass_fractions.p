reset

# Font file
cmsy='/Users/Mead/Fonts/cmsy10.pfb'

# Output options
if(!exists('print')){print=0}
if(print==0){set term qt dashed font ',16'; sun='sun'; width=3}
if(print==1){set term post enh col font ',20' fontfile cmsy; set output 'paper/halo_mass_fractions.eps'; sun='{/cmsy10 \014}'; width=6}

# File to plot
file='data/mass_fractions.dat'

# Range for mass axis
mmin=1e10
mmax=3e15

# Range for mass-fraction axis
rmin=3e-3
rmax=2.

# X axis properties
set log x
set xlabel 'M / h^{-1} M_{'.sun.'}'
set xrange [mmin:mmax]
set format x '10^{%T}'
set mxtics 10

# Y axis properties
set log y
set ylabel 'Halo mass fraction'
set yrange [rmin:rmax]
set format y

#set key outside left box
set key top left box opaque font ',15'

# Cosmological paarameters taken from WMAP9
om_b=0.0463
om_m=0.2793

# Ploty plot plot
plot 1 w l lc -1 lw width noti 'Total',\
   file u 1:(1.-column(5)) w l lw width dt 1 lc 1 ti 'Total',\
   om_b/om_m  w l lw width dt 2 lc -1 ti 'Universal baryon',\
   file u 1:2 w l lw width dt 1 lc 2  ti 'CDM',\
   file u 1:3 w l lw width dt 2 lc 3  ti 'Total gas',\
   file u 1:4 w l lw width dt 1 lc 3  ti 'Bound gas',\
   file u 1:5 w l lw width dt 3 lc 3  ti 'Unbound gas',\
   file u 1:6 w l lw width dt 1 lc 4  ti 'Total stars',\
   file u 1:7 w l lw width dt 2 lc 4  ti 'Central stars',\
   file u 1:8 w l lw width dt 3 lc 4  ti 'Satellite stars'#,\
   file u 1:9 w l lw width dt 3 lc 5  ti 'Neutrinos'

show output

