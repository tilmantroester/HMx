reset

file='fitting/debug.dat'

set xlabel 'parameter value'

set ylabel 'figure of merit'

plot file u 1:2 w l lw 3 noti
