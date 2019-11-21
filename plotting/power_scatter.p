reset

if(!exists('print')){print=0}
if(print==0) {set term qt dashed}
if(print==1) {set term post enh col; set output 'power_scatter.eps'}

pow_norm='data/power_hm.dat'
pow_scat='data/power_scatter_hm.dat'

set multiplot

set log x
set xlabel 'k / h^{-1} Mpc'

set ylabel 'P_{scatter}(k) / P_{normal}(k)'

na=16
amin=0.1
amax=1.
set cbrange [amin:amax]
set palette cubehelix negative
set cblabel 'a'

plot 1 w l lt -1 noti,\
     for [i=1:na] '<paste '.pow_scat.' '.pow_norm.'' u 1:(column(1+i)/column(2+na+i)):(amin+(amax-amin)*(real(i-1)/real(na-1))) w l lc palette lw 3 noti

top=0.9
bot=0.45
lef=0.2
rig=0.6

set tmargin at screen top
set bmargin at screen bot
set lmargin at screen lef
set rmargin at screen rig

unset log x
set xlabel 'c'

set ylabel 'p(c)'
unset ytics

plot 'data/p_conc.dat' u 1:2 w l lw 2 lc -1 noti

unset multiplot
