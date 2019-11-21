unset multiplot
reset

# Font file
cmmi='/Users/Mead/Fonts/cmmi10.pfb'

# Print options
if(!exists('print')){print=0}
if(print==0){set term qt; ell='l'}
if(print==1){set term post enh col fontfile cmmi; set output 'CMB_lensing.eps'; ell='{/cmmi10 \140}'}

# File paths
#CAMB='/Users/Mead/Physics/CAMB_files/boring/boring_scalCls.dat'
CAMB='/Users/Mead/Physics/HMx/benchmarks/CAMB_boring_scalCls.txt'
HMx='/Users/Mead/Physics/HMx/data/CMBlensing_cl_linear.dat'

# CMB temperature in K
T0=2.725

# Factor to convert C(l) from dimensionless to uK^2
# Note T0 converts to K
# 1e6 converts to uK
# 2 converts between kappa and Phi (Del^2 Phi = 2 kappa)
fac=(1e6*T0)**2

# ell axis
ellmin=2
ellmax=10000
set log x
set xrange [ellmin:ellmax]

set lmargin 10
set rmargin 2

set multiplot layout 2,1

set xlabel ''
set format x ''

# C(ell) axis
set log y
set ylabel ''.ell.'^4 C_{'.ell.'}^{{/Symbol F}{/Symbol F}} [{/Symbol m}K]^2'
set format y '10^{%T}'

# In Fourier Space [l(l+1)]^2 Phi = 2 kappa
# CAMB files make l^4 C_l^PhiPhi
# So need to multiply my result by [2*l/(l+1)]^2
plot CAMB u 1:5 w l lw 2 lc -1 ti 'CAMB',\
     HMx u 1:(((2*$1/($1+1))**2)*(fac*$2)) w l lw 2 lc 1 ti 'HMx'

set xlabel ''.ell.''
set format x '10^{%T}'

rmin = 0.9
rmax = 1.1
unset log y
set format y
set yrange[rmin:rmax]
set ylabel 'C_{HMx}('.ell.')/C_{CAMB}('.ell.')'

plot 1 w l lt -1 noti,\
     0.99 w l dt 2 lc -1 noti,\
     1.01 w l dt 2 lc -1 noti,\
     0.95 w l dt 3 lc -1 noti,\
     1.05 w l dt 3 lc -1 noti,\
     '<paste '.HMx.' '.CAMB.'' u 1:((((2*$1/($1+1))**2)*(fac*$2))/$8) w l lc 1 lw 2 noti

unset multiplot


