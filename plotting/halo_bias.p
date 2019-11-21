unset multiplot
reset

if(!exists('print')) {print=0}
if(print==0) {set term qt dashed}
if(print==1) {set term post enh col; set output 'halo_bias.eps'}

power(f1,f2)=sprintf('data/power_%s%s_hm.dat',f1,f2)
simul(snap,f1,f2)=sprintf('/Users/Mead/Physics/BAHAMAS/power/M1024/DMONLY_nu0_L400N1024_WMAP9_%s_%s_%s_power.dat',snap,f1,f2)

print ''

if(!exists('iplot')){iplot=1}
print 'iplot = 1: Plot both power and ratio'
print 'iplot = 2: Plot only ratio '
print 'iplot = 3: Plot only power'
print 'iplot = ', iplot

print ''

kmin=1e-2
kmax=1e1
klab='k / h Mpc^{-1}'

pmin=1e-6
pmax=1e3
plab='{/Symbol D}^2_{uv}(k)'

rmin=0.
rmax=10.
rlab='b_h(k)'

# Point size
size=0.5

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

if(iplot==1 || iplot==3){

if(iplot==1) {set multiplot layout 2,1}

if(iplot==1){set xlabel ''; set format x ''}
set xlabel klab

set yrange [pmin:pmax]
set log y
set ylabel plab
set format y '10^{%T}'

set key top left
if(iplot==1) {unset key}

plot for [i=1:nz] NaN lc i lw 3 ti word(names,i),\
     NaN lc -1 lw 3 dt 1 ti 'matter-matter',\
     NaN lc -1 lw 3 dt 2 ti 'matter-halo',\
     NaN lc -1 lw 3 dt 3 ti 'halo-halo',\
     for [i=1:nz] power('m','m') u 1:(column(i+1)/10**(i-1)) w l lw 3 lc i dt 1 noti,\
     for [i=1:nz] power('m','f') u 1:(column(i+1)/10**(i-1)) w l lw 3 lc i dt 2 noti,\
     for [i=1:nz] power('f','f') u 1:(column(i+1)/10**(i-1)) w l lw 3 lc i dt 3 noti,\
     for [i=1:nz] simul(word(snaps,i),'all','all') u 1:($2/10**(i-1)):($5/10**(i-1)) w e pt 7 ps size lc i noti,\
     for [i=1:nz] simul(word(snaps,i),'all','FOF') u 1:($2/10**(i-1)):($5/10**(i-1)) w e pt 6 ps size lc i noti,\
     for [i=1:nz] simul(word(snaps,i),'FOF','FOF') u 1:($2/10**(i-1)):($5/10**(i-1)) w e pt 5 ps size lc i noti

}

if(iplot==1 || iplot==2) {

set key top left

set xlabel klab
set format x

set ylabel rlab
unset log y
set yrange [rmin:rmax]
set format y

plot 1 w l lt -1 noti,\
     NaN w l lc -1 dt 2 lw 3 ti 'Cross-correlation bias',\
     NaN w l lc -1 dt 1 lw 3 ti 'Auto-correlation bias',\
     for [i=1:nz] '<paste '.power('m','f').' '.power('m','m').'' u 1:(column(i+1)/column(1+i+L))       w l lw 3 lc i dt 2 noti,\
     for [i=1:nz] '<paste '.power('f','f').' '.power('m','m').'' u 1:(sqrt(column(i+1)/column(1+i+L))) w l lw 3 lc i dt 1 ti word(names,i),\
     for [i=1:nz] '<paste '.simul(word(snaps,i),'all','FOF').' '.simul(word(snaps,i),'all','all').'' u 1:(column(2)/column(2+L))       w p pt 6 ps size lc i dt 2 noti,\
     for [i=1:nz] '<paste '.simul(word(snaps,i),'FOF','FOF').' '.simul(word(snaps,i),'all','all').'' u 1:(sqrt(column(2)/column(2+L))) w p pt 7 ps size lc i dt 2 noti

if(iplot==1) {unset multiplot}

}
