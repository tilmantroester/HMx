unset multiplot
reset

if(print==0) {set term aqua size 1000,400}
if(print==1) {set term post enh col font ',10' size 10,4; set output 'paper/power_cosmology.eps'}

power(i,f1,f2)=sprintf('data/cosmology_%d_%s-%s_power.dat',i,f1,f2)

set log x
set xlabel 'k / h Mpc^{-1}'

rmin=0.4
rmax=2.2
set yrange [rmin:rmax]
set ylabel 'P(k) / P_{mid}(k)'

cbmin=0.7
cbmax=0.9
set palette defined (1 'light-blue', 2 'black', 3 'light-red')
set cblabel '{/Symbol s}_8'
set cbrange [cbmin:cbmax]

labx=0.1
laby=0.95

ncos=15
f(i)=cbmin+(cbmax-cbmin)*real(i-1)/real(ncos-1)

set multiplot layout 1,2

do for [j=1:2]{

if(j==1) {f1='matter'; f2='matter'; f1_name='matter'; f2_name='matter'}
if(j==2) {f1='matter'; f2='epressure'; f1_name='matter'; f2_name='electron pressure'}

set label ''.f1_name.' - '.f2_name.'' at graph labx,laby
plot for [i=1:ncos] '<paste '.power(i,f1,f2).' '.power((ncos+1)/2,f1,f2).'' u 1:($5/$10):(f(i)) w l lw 2 lc palette noti

unset label

}

unset multiplot

show output
