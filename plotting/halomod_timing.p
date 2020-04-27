reset

file = 'data/halomod_timing.dat'

tmin = 0.05
tmax = 2.

set log x
set xlabel 'Number of redshifts'

set log y
set ylabel 'Runtime [s]'

set yrange [tmin:tmax]

# Rough model
t_cos = 0.05
t_hmx = 0.008
model(n) = t_cos+t_hmx*n

set key top left

plot model(x) w l lw 2 lc -1 ti 'Simple model',\
   file u 1:2 w p pt 7 ti 'Data'
   