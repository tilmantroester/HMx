reset

if(!exists("print")){print=0}
if(print==0){set term qt dashed font ',14'}
if(print==1){set term post enh col sol; set output 'power_halovoid.eps'}
if(print==2){set term pdfcairo; set output 'power.pdf'}

# Initial white space
print ''

if(!exists('delta')) {delta=1}
if(delta==0) {fac=4.*pi/(2.*pi)**3; p=3; pmin=1e-3; pmax=1e5; plab='P(k) / (h^{-1} Mpc)^3'}
if(delta==1) {fac=1.; p=0; pmin=1e-8; pmax=1e4; plab='{/Symbol D}^2(k)'}

kmin=1e-3
kmax=1e2

set log x
set xrange [kmin:kmax]
set xlabel 'k / h Mpc^{-1}'
set mxtics 10

set log y
set yrange [pmin:pmax]
set ylabel plab
set format y '10^{%T}'

col='red'

unset colorbox

power='data/power_halovoid.dat'

#Key stuff
if(delta==0) set key top right
if(delta==1) set key top left

# write
print 'delta = 0: Plot P(k)'
print 'delta = 1: Plot Delta^2(k)'
print 'delta = ', delta
print ''

# Now do the actual plotting
plot power u 1:($2/(fac*$1**p)) w l lc -1       dt 1 lw 3 ti 'Linear',\
     power u 1:($3/(fac*$1**p)) w l lc rgb col  dt 2 lw 3 ti '2-halo',\
     power u 1:($4/(fac*$1**p)) w l lc rgb col  dt 3 lw 3 ti '1-halo',\
     power u 1:($6/(fac*$1**p)) w l lc rgb col  dt 4 lw 3 ti '1-void',\
     power u 1:($5/(fac*$1**p)) w l lc rgb col  dt 1 lw 3 ti 'Total'





