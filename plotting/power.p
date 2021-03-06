reset

if(!exists("print")){print=0}
if(print==0){set term qt dashed font ',14'}
if(print==1){set term post enh col sol; set output 'power.eps'}
if(print==2){set term pdfcairo; set output 'power.pdf'}

# Initial white space
print ''

# Choose what to plot
if(!exists("imode")){imode=1}
print 'imode = 1: Plot P(k) breakdown'
print 'imode = 2: Plot Delta^2(k) breakdown'
print 'imode = 3: Plot P(k)'
print 'imode = 4: Plot Delta^2(k)'
print 'imode = 5: Ratio compared to linear'
print 'imode = 6: Ratio compared to non-linear'
print 'imode = ', imode
print ''

if(imode==1 || imode==3) {pmin=1e-2; pmax=1e5; plab='P(k) / (h^{-1} Mpc)^3'; fac=4.*pi/(2.*pi)**3; p=3}
if(imode==2 || imode==4) {pmin=1e-7; pmax=1e4; plab='{/Symbol D}^2(k)'; fac=1.; p=0}
if(imode==5) {pmin=0.9; pmax=1.1; plab='P(k) / P_{linear}(k)'}
if(imode==6) {pmin=0.9; pmax=1.1; plab='P(k) / P_{non-linear}(k)'}

# k axis
kmin=1e-3
kmax=1e2
set log x
set xrange [kmin:kmax]
set xlabel 'k / h Mpc^{-1}'
set mxtics 10

# power axis
set yrange [pmin:pmax]
set ylabel plab
if(imode==1 || imode==2 || imode==3 || imode==4) {set log y; set format y '10^{%T}'}
if(imode==5) {unset log y}

# File to plot
power='data/power.dat'

#Key stuff
if(imode==1 || imode==3 || imode==5 || imode==6) set key top right
if(imode==2 || imode==4) set key top left

#Now do the actual plotting
if(imode==1 || imode==2){
plot power u 1:($2/(fac*$1**p)) w l lc -1 dt 2 lw 5 ti 'Linear',\
     power u 1:($3/(fac*$1**p)) w l lc  7 dt 1 lw 4 ti '2-halo',\
     power u 1:($4/(fac*$1**p)) w l lc  6 dt 1 lw 4 ti '1-halo',\
     power u 1:($5/(fac*$1**p)) w l lc -1 dt 1 lw 4 ti 'Total'
}

#Now do the actual plotting
if(imode==3 || imode==4){
plot power u 1:($2/(fac*$1**p)) w l lc -1 dt 2 lw 4 ti 'Linear',\
     power u 1:($5/(fac*$1**p)) w l lc  7 dt 1 lw 4 ti 'Non-linear'
}

if(imode==5){
plot 1 w l lt -1 noti,\
   power u 1:($3/$2) w l lc 6 dt 2 lw 4 ti 'Two-halo term',\
   power u 1:($5/$2) w l lc 6 dt 1 lw 4 ti 'Halo model'
}

if(imode==6){
plot 1 w l lt -1 noti,\
   power u 1:($2/$5) w l lc 7 dt 1 lw 4 ti 'Linear',\
   power u 1:($3/$5) w l lc 7 dt 2 lw 4 ti 'Two-halo term',\
   power u 1:($4/$5) w l lc 7 dt 3 lw 4 ti 'One-halo term'
}

show output
unset output
