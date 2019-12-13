reset

if(!exists('print')) {print=0}
if(print==0) {set term qt dashed}

#Type of spectrum to plot
print('')
if(!exists('itype')) {itype=4}
print('Choose what power to print')
print('itype = 1: Linear')
print('itype = 2: Two-halo')
print('itype = 3: One-halo')
print('itype = 4: Full')
print('itype = '.itype.'')
print('')

# Files
file_li='data/power_linear.dat'
file_hm='data/power_hm.dat'
file_1h='data/power_1h.dat'
file_2h='data/power_2h.dat'

# Baseline files
file_base_li='data/power_baseline_linear.dat'
file_base_hm='data/power_baseline_hm.dat'
file_base_1h='data/power_baseline_1h.dat'
file_base_2h='data/power_baseline_2h.dat'

#Set k-range
kmin=1e-3
kmax=1e2


#Set power axis
pmin = 1e-8
pmax = 1e4


unset colorbox

#For low-k power-law function
#A=1e-11 #Appropriate for dd 1-halo term
#A=1e-28 #Appropriate for dp 1-halo term
#k0=kmin
#ns=0.
#kmax=2e-2
#f(k) = k<kmax ? A*(k/k0)**(3+ns) : 1/0

#Key stuff
unset key

#Number of redshifts
n=16

#Colour stuff
set palette defined (1 'dark-red', 2 'gold')

#Set the file type to plot
if(itype==1){file=file_li; file_base=file_base_li; tits='Linear power'}
if(itype==2){file=file_2h; file_base=file_base_2h; tits='Two-halo power'}
if(itype==3){file=file_1h; file_base=file_base_1h; tits='One-halo power'}
if(itype==4){file=file_hm; file_base=file_base_hm; tits='Halo-model power'}
print('File: '.file.'')
print('')
print('Title: '.tits.'')
print('')

set lmargin 10
set rmargin 2
set multiplot layout 2,1

set log y
#set yrange [pmin:pmax]
set ylabel '{/Symbol D}_{i,j}^2(k)'
set format y '10^{%T}'

set log x
set xrange [kmin:kmax]
set xlabel ''
set format x ''

#Actual plot
set title tits
plot for[i=1:n] file u 1:(column(i+1)):(real(i-1)/real(n)) w l lw 2 dt 1 lc palette noti,\
     for[i=1:n] file_base u 1:(column(i+1)):(real(i-1)/real(n)) w l lw 2 dt 1 lc palette noti

set xlabel 'k / (h Mpc^{-1})'
set mxtics 10
set format x

unset log y
set ylabel 'P(k) / P_{base}(k)'
set format y

plot 1 w l lt -1 noti,\
   for[i=1:n] '<paste '.file.' '.file_base.'' u 1:(column(i+1)/column(i+2+n)):(real(i-1)/real(n)) w l lw 2 dt 1 lc palette noti

unset multiplot

show output

