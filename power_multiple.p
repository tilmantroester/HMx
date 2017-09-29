reset

if(print==0) set term aqua dashed

#This gnuplot script plots the output file 'power.dat' that is spat out of HMcode.
#Simply load up gnuplot (type gnuplot in the terminal) and then type "gnuplot>load 'plot.p'"
#The plot should then be the non-linear spectrum at 16 redshifts

kmin=0.001
kmax=100.
set log x
set xrange [kmin:kmax]
set xlabel 'k / (h Mpc^{-1})'
set mxtics 10

set log y
#set yrange [1e-8:1e4]
set ylabel '{/Symbol D}^2(k)'
set format y '10^{%T}'

unset colorbox

file_linear='data/power_linear.dat'
file_full='data/power_full.dat'
file_1halo='data/power_1halo.dat'
file_2halo='data/power_2halo.dat'

#A=1e-11 #Appropriate for dd 1-halo term
A=1e-28 #Appropriate for dp 1-halo term
k0=kmin
ns=0.
kmax=2e-2
f(k) = k<kmax ? A*(k/k0)**(3+ns) : 1/0

#Key stuff
unset key

#Number of redshifts
n=16

# start value for H
#h1 = 117/360.0
# end value for H
#h2 = 227/360.0
# creating the palette by specifying H,S,V
#set palette model HSV functions (1-gray)*(h2-h1)+h1,1,0.68

set palette defined (1 'dark-red', 2 'gold')

#col(i)=sprintf("#%1x%1x0000",i-1,i-1)
#plot for[i=1:16] file u 1:(column(i+1)) w l lw 2 lc rgb col(i) noti

#Now do the actual plotting
plot file_2halo  u 1:(column(n+1)):(real(n-1)/real(n)) w l lw 2 dt 2 lc palette,\
     file_1halo  u 1:(column(n+1)):(real(n-1)/real(n)) w l lw 2 dt 2 lc palette,\
     for[i=1:n] file_full u 1:(column(i+1)):(real(i-1)/real(n)) w l lw 2 dt 1 lc palette noti,\
     f(x) w l lw 3 dt 2 lc -1 noti




