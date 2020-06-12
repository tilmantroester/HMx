unset multiplot
reset

file = 'data/HMcode_baryon_model_power.dat'

mmin = 1e12
mmax = 1e16
nm = 9

klab = 'k / h Mpc^{-1}'

plab = '{/Symbol D}^2(k)'

set log x
set xlab klab

set log y
set ylab plab
set format y '10^{%T}'

unset key
set log cb
set cbrange [mmin:mmax]
set colorbox
set format cb '10^{%T}'

mbar(i, nm) = 10**(log10(mmin)+log10(mmax/mmin)*(i-1)/(nm-1))

set palette defined (1 'royalblue', 2 'grey', 3 'light-red')

set multiplot layout 2, 1 margins 0.06, 0.90, 0.06, 0.98

plot file u 1:2 w l lc -1 lw 2,\
   for [i=1:nm] file u 1:(column(i+2)):(mbar(i,nm)) w l lw 2 lc palette

unset log y
set format y

plot for [i=1:nm] file u 1:(column(i+2)/column(2)):(mbar(i, nm)) w l lw 2 lc palette

unset multiplot