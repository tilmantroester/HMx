reset

print ''

print 'Choose terminal'
print 'print = 0: aqua'
print 'print = 1: eps'
if(!exists('print')) {print=0}
if(print == 0) {set term aqua dashed}
if(print == 1) {set term post enh col; set output 'HMF_amplitude.eps'}
print 'print = ', print
print ''

power(i) = sprintf('data/power_HMFamp_%d.dat',i)
base = 'data/power_HMFamp_fid.dat'

n = 16

set log x
set xlabel 'k / h Mpc^{-1}'

set log y
set format y '10^{%T}'
set ylabel '{/Symbol D}^2(k)'

amp_min=0.1
amp_max=2.0
amp(i,n)=amp_min+(amp_max-amp_min)*(i-1.)/(n-1.)
set cbrange [amp_min:amp_max]
set cblabel '{/Symbol a}'
set palette cubehelix start 1.5 cycles -0.5 saturation 2.0

print 'Minimum amplitude: ', amp_min
print 'Maximum amplitude: ', amp_max
print 'Check that these are consistent with the output from code'
print ''

set key top left

plot for [i=1:n] power(n-i+1) u 1:5:(amp(n-i+1,n)) w l lw 3 lc palette noti#,\
   base u 1:5 w l lw 3 dt 2 lc -1 ti 'Fiducial'

show output
