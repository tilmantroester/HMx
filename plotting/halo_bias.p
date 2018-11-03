unset multiplot
reset

set term aqua dashed

power(f1,f2)=sprintf('data/power_%s%s_hm.dat',f1,f2)
simul(snap,f1,f2)=sprintf('/Users/Mead/Physics/BAHAMAS/power/M1024/DMONLY_nu0_L400N1024_WMAP9_%s_%s_%s_power.dat',snap,f1,f2)

print ''

if(!exists('iplot')){iplot=1}
print 'iplot = 1: Plot both raw power and ratio'
print 'iplot = 2: Plot only ratio '
print 'iplot = ', iplot

print ''

kmin=1e-2
kmax=1e1

pmin=1e-4
pmax=1e3

rmin=0.
rmax=10.

# 1 - z = 0.0
# 2 - z = 0.5
# 3 - z = 1.0
# 4 - z = 2.0
snaps="'snap32' 'snap28' 'snap26' 'snap22'"
names="'z = 0.0' 'z = 0.5' 'z = 1.0' 'z = 2.0'"

# File lengths
L=5

# Number of redshifts to plot (from zero)
nz=4

set log x
set xrange [kmin:kmax]

set lmargin 10
set rmargin 2

if(iplot==1){

set multiplot layout 2,1

set xlabel ''
set format x ''

set yrange [pmin:pmax]
set log y
set ylabel '{/Symbol D}^2_{uv}(k)'
set format y '10^{%T}'

set key top left

plot for [i=1:nz] power('m','m') u 1:(column(i+1)) w l lw 3 lc -1      noti 'matter-matter',\
     for [i=1:nz] power('m','f') u 1:(column(i+1)) w l lw 3 lc i  dt 2 noti 'matter-halo',\
     for [i=1:nz] power('f','f') u 1:(column(i+1)) w l lw 3 lc i  dt 1 noti 'halo-halo',\
     for [i=1:nz] simul(word(snaps,i),'all','all') u 1:2:5 w e pt 7 lc -1 noti,\
     for [i=1:nz] simul(word(snaps,i),'all','FOF') u 1:2:5 w e pt 6 lc i  noti,\
     for [i=1:nz] simul(word(snaps,i),'FOF','FOF') u 1:2:5 w e pt 7 lc i  noti

}

set key top left

set xlabel 'k / h Mpc^{-1}'
set format x

set ylabel 'b_h(k)'
unset log y
set yrange [rmin:rmax]

plot 1 w l lt -1 noti,\
     NaN w l lc -1 dt 2 lw 3 ti 'Cross-correlation bias',\
     NaN w l lc -1 dt 1 lw 3 ti 'Auto-correlation bias',\
     for [i=1:nz] '<paste '.power('m','f').' '.power('m','m').'' u 1:(column(i+1)/column(1+i+L))       w l lw 3 lc i dt 2 noti,\
     for [i=1:nz] '<paste '.power('f','f').' '.power('m','m').'' u 1:(sqrt(column(i+1)/column(1+i+L))) w l lw 3 lc i dt 1 ti word(names,i),\
     for [i=1:nz] '<paste '.simul(word(snaps,i),'all','FOF').' '.simul(word(snaps,i),'all','all').'' u 1:(column(2)/column(2+L))       w p pt 6 lc i dt 2 noti,\
     for [i=1:nz] '<paste '.simul(word(snaps,i),'FOF','FOF').' '.simul(word(snaps,i),'all','all').'' u 1:(sqrt(column(2)/column(2+L))) w p pt 7 lc i dt 2 noti

if(iplot==1) {unset multiplot}
