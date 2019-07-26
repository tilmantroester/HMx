reset

unset multiplot

load '/Users/Mead/Physics/library/gnuplot/mead.p'

if(!exists('print')) {print=0}
if(print==0) {set term aqua dashed}
if(print==1) {set term post enh col; set output 'HMcode_CAMB_comparison.eps'}

HMx='data/power_HMcode_HMx.dat'
CAMB='data/power_HMcode_CAMB.dat'

kmin=1e-4
kmax=1e2
#kmin=1e0
#kmax=1e1

pmin=1e-9
pmax=1e3
#pmin=1.
#pmax=2.

rmin=0.98
rmax=1.02
#rmin=0.9
#rmax=1.1
dr=0.01
ddr=0.001

amin=0.333333
amax=1.0
#amin=0.0
#amax=1.0
na=16
ia1=1
ia2=na

top = 0.98
bot = 0.1
lef = 0.1
rig = 0.88
spc = 0.02
cbs = 0.03

set log x
set xrange [kmin:kmax]

set multiplot layout 2,1 margins lef,rig,top,bot spacing spc

set log y
set format y '10^{%T}'
set ylabel '{/Symbol D}^2(k)'
set yrange [pmin:pmax]

set format x ''
set xlabel ''

set cbrange [amin:amax]
set colorbox vertical user origin rig+0.02, bot size cbs, top-bot
set cblabel 'a'
#set palette rgb color_scheme(rainbow[1], rainbow[2], rainbow[3])
#set palette rgb color_scheme(AFM_hot[1], AFM_hot[2], AFM_hot[3])
#set palette rgb 34,35,36
#set palette rgb 33,13,10

set key top left

# MEAD: Revert -1 -> Palette
plot\
   NaN w l lw 3 lc -1 dt 1 ti 'HMx',\
   NaN w l lw 3 lc -1 dt 2 ti 'CAMB',\
   for [i=ia1:ia2] HMx  u 1:(column(i+1)):(progression(amin,amax,i,na)) w l lw 3 lc palette dt 1 noti 'HMx',\
   for [i=ia1:ia2] CAMB u 1:(column(i+1)):(progression(amin,amax,i,na)) w l lw 3 lc palette dt 2 noti 'CAMB'


unset colorbox

set format x '10^{%T}'
set xlabel 'k / h Mpc^{-1}'

unset log y
set yrange[rmin:rmax]
set ylabel 'P_{HMx}(k) / P_{CAMB}(k)'
set format y

plot 1 w l lt -1 noti,\
   1.+dr w l dt 2 lc -1 noti,\
   1.-dr w l dt 2 lc -1 noti,\
   1.+ddr w l dt 3 lc -1 noti,\
   1.-ddr w l dt 3 lc -1 noti,\
   for [i=ia1:ia2] '<paste '.HMx.' '.CAMB.'' u 1:(column(i+1)/column(na+2+i)):(progression(amin,amax,i,na)) w l lw 3 lc palette noti

unset multiplot

show output