unset multiplot
reset

set term qt dashed dl 1

benc(f1,f2)=sprintf('benchmarks/cl_%s_%s.txt',f1,f2)
data(f1,f2)=sprintf('data/cl_%s_%s.dat',f1,f2)

f1s="'RCSLenS' 'CMB'"
f2s="'RCSLenS' 'y'"

tracers="'RCSLenS' 'y' 'CMB'"

set log x

set lmargin 10
set rmargin 3

set multiplot layout 2,1

set format x ''

set log y
set ylabel 'C(l)'
set format y '10^{%T}'

plot for [i=1:words(tracers)] NaN w l lw 2 lc i ti word(tracers,i),\
     NaN w l lw 2 dt 3 lc -1 ti 'Benchmark',\
     NaN w l lw 2 dt 2 lc -1 ti 'HMx',\
     for [i=1:words(tracers)] for [j=1:words(tracers)] benc(word(tracers,i),word(tracers,j)) u 1:2 w l lc 0 lw 2 dt 1 noti,\
     for [i=1:words(tracers)] for [j=1:words(tracers)] data(word(tracers,i),word(tracers,j)) u 1:2 w l lc i lw 2 dt 1 noti,\
     for [i=1:words(tracers)] for [j=1:words(tracers)] data(word(tracers,i),word(tracers,j)) u 1:2 w l lc j lw 2 dt 2 noti

#for [i=1:words(f1s)] benc(word(f1s,i),word(f2s,i)) u 1:2 w l lc i lw 2 dt 3 noti,\
#for [i=1:words(f1s)] data(word(f1s,i),word(f2s,i)) u 1:2 w l lc i lw 2 dt 2 noti

set xlabel 'l'
set format x

unset log y
set ylabel 'C(l) / C(l)'
set format y
dy=0.05
set yrange [1.-dy:1.+dy]

ddy=0.01
dddy=0.001

plot 1 w l lt -1 noti,\
     1.-ddy  w l lc -1 dt 3 noti,\
     1.+ddy  w l lc -1 dt 3 noti,\
     1.-dddy w l lc -1 dt 3 noti,\
     1.+dddy w l lc -1 dt 3 noti,\
     for [i=1:words(tracers)] for [j=1:words(tracers)] '<paste '.data(word(tracers,i),word(tracers,j)).' '.benc(word(tracers,i),word(tracers,j)).'' u 1:($2/$5) w l lc i lw 2 dt 1 noti,\
     for [i=1:words(tracers)] for [j=1:words(tracers)] '<paste '.data(word(tracers,i),word(tracers,j)).' '.benc(word(tracers,i),word(tracers,j)).'' u 1:($2/$5) w l lc j lw 2 dt 2 noti

     #for [i=1:words(f1s)] '<paste '.data(word(f1s,i),word(f2s,i)).' '.benc(word(f1s,i),word(f2s,i)).'' u 1:($2/$5) w l lc i lw 2 dt 1 noti

unset multiplot
