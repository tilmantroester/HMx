reset

CAMB='/Users/Mead/Physics/CAMB/boring/boring_scalCls.dat'
HMx='/Users/Mead/Physics/HMx/data/cl_linear.dat'

ell='l'

T0=2.725

fac=(1e6*T0)**2

set log x
set xlabel ''.ell.''
set format x '10^{%T}'
set xrange [2:2000]

set log y
set ylabel ''.ell.'^4 C_{{/Symbol F}{/Symbol F}}('.ell.') [{/Symbol m}K]^2'
set format y '10^{%T}'

plot CAMB u 1:5 w l lw 3 ti 'CAMB',\
     HMx u 1:(fac*$2) w l lw 3 ti 'HMx'
