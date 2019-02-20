reset

file='data/Tinker_bias.dat'

numin=-0.525
numax=0.625

bmin=0
bmax=9

set xrange [numin:numax]
set xlabel 'log_{10}({/Symbol n})'

set yrange [bmin:bmax]
set ylabel 'b'

set title 'To be compared with Fig. 1 of Tinker et al. (2010)'

plot file u (log10($1)):2 w l lw 3 noti
