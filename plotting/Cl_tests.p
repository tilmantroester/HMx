unset multiplot
reset

benc(f1,f2)=sprintf('benchmarks/cl_%s_%s.txt',f1,f2)
data(f1,f2)=sprintf('data/cl_%s_%s.dat',f1,f2)

f1s="'RCSLenS' 'CMB' 'KiDS'"
f2s="'RCSLenS' 'y' 'y'"

set log x

set lmargin 10
set rmargin 3

set multiplot layout 2,1

set format x ''

set log y
set ylabel 'C(l)'
set format y '10^{%T}'

plot NaN w l lw 2 dt 1 lc -1 ti 'Benchmark',\
     NaN w l lw 2 dt 2 lc -1 ti 'HMx',\
     for [i=1:words(f1s)] benc(word(f1s,i),word(f2s,i)) u 1:2 w l lc i lw 2 dt 1 ti ''.word(f1s,i).'-'.word(f2s,i).'',\
     for [i=1:words(f1s)] data(word(f1s,i),word(f2s,i)) u 1:2 w l lc i lw 2 dt 2 noti

set xlabel 'l'
set format x

unset log y
set ylabel 'C(l) / C(l)'
set format y
dy=0.15
set yrange [1.-dy:1.+dy]

ddy=0.01
dddy=0.001

plot 1 w l lt -1 noti,\
     1.-ddy  w l lc -1 dt 3 noti,\
     1.+ddy  w l lc -1 dt 3 noti,\
     1.-dddy w l lc -1 dt 3 noti,\
     1.+dddy w l lc -1 dt 3 noti,\
     for [i=1:words(f1s)] '<paste '.data(word(f1s,i),word(f2s,i)).' '.benc(word(f1s,i),word(f2s,i)).'' u 1:($2/$5) w l lc i lw 2 dt 1 noti

unset multiplot
