unset multiplot
reset

# Set terminal options
if(!exists('print')) {print=0}
if(print==0) {set term qt dashed}
if(print==1) {set term post enh col; set output 'power_multiple_comparison.eps'}

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
file_2h='data/power_2h.dat'
file_1h='data/power_1h.dat'
file_hm='data/power_hm.dat'

# Baseline files
file_base_li='data/power_baseline_linear.dat'
file_base_2h='data/power_baseline_2h.dat'
file_base_1h='data/power_baseline_1h.dat'
file_base_hm='data/power_baseline_hm.dat'

# Set k-range
kmin=1e-3
kmax=1e2

# Set power axis
pmin = 1e-8
pmax = 1e4

# Set ratio axis
rmin = 0.65
rmax = 1.05

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

# Number of redshifts
if(!exist('nz')){nz = 16}
print('Number of z being plotted: nz: ', nz)

# Colour stuff
set palette defined (1 'dark-red', 2 'gold')

# Set the file type to plot
if(itype==1){file=file_li; file_base=file_base_li; tits='Linear power'}
if(itype==2){file=file_2h; file_base=file_base_2h; tits='Two-halo power'}
if(itype==3){file=file_1h; file_base=file_base_1h; tits='One-halo power'}
if(itype==4){file=file_hm; file_base=file_base_hm; tits='Halo-model power'}
print('File: '.file.'')
print('')
print('Title: '.tits.'')
print('')

# k axis
set log x
set xrange [kmin:kmax]
set xlabel ''
set format x ''

# Margin options and multiplot
set lmargin 10
set rmargin 2
set multiplot layout 2,1

# Power axis
set log y
#set yrange [pmin:pmax]
set ylabel '{/Symbol D}@^2_{i,j}(k)'
set format y '10^{%T}'

# Plot power
plot for[i=1:nz] file u 1:(column(i+1)):(real(i-1)/real(nz)) w l lw 2 dt 1 lc palette not,\
     for[i=1:nz] file_base u 1:(column(i+1)):(real(i-1)/real(nz)) w l lw 1 dt 3 lc palette noti

# x axis
set xlabel 'k / (h Mpc^{-1})'
set mxtics 10
set format x

# ratio axis
unset log y
#set yrange [rmin:rmax]
set ylabel 'P(k) / P_{base}(k)'
set format y

# Plot response
plot 1 w l lt -1 noti,\
   for[i=nz:nz] '<paste '.file_li.' '.file_base_li.'' u 1:(column(i+1)/column(i+2+nz)):(real(i-1)/real(nz)) w l lw 2 dt 2 lc -1 noti,\
   for[i=1:nz] '<paste '.file.' '.file_base.'' u 1:(column(i+1)/column(i+2+nz)):(real(i-1)/real(nz)) w l lw 2 dt 1 lc palette noti
   

unset multiplot

show output
