reset
unset multiplot

cmmi='/Users/Mead/Fonts/cmmi10.pfb'

if(print==0){set term aqua dashed; ell='l'}
if(print==1){set term post enh col fontfile cmmi; set output 'Cl.eps'; ell='{/cmmi10 \140}'}

#Files to plot
twohalo='data/cl_2halo.dat'
onehalo='data/cl_1halo.dat'
full='data/cl_full.dat'

#Power-law plotting
ns=0.
#A=1e-10 #Appropriate for kk one-halo term
#A=1e-13 #Appropriate for ky one-halo term
A=1e3
lmax=20
f(l) = l<lmax ? A*(l+1)*l**(1+ns) : 1/0
tits(n)=sprintf('C('.ell.') \~ l^{%1.1f}',n)

icl=2
if(icl==1) c=2
if(icl==2) c=3

set log x
set xrange [1e0:1e5]
#set xrange [1e-3:1e5]
set mxtics 10

set lmargin 10
set rmargin 2

set multiplot layout 2,1

set format x ''
set xlabel ''

set log y
set mytics 10
#if(print==0){set ylabel 'l(l+1)C_{i,j}(l) / 2{/Symbol p}'}
#if(print==1){
if(icl==1){set ylabel 'C_{i,j}('.ell.')'}
if(icl==2){set ylabel ''.ell.'('.ell.'+1)C_{i,j}('.ell.') / 2{/Symbol p}'}
set format y '10^{%T}'
#set yrange [1e-8:1e-4]
#set yrange [1e-12:1e-8]

set key top left

plot twohalo  u 1:c w l lw 3 lc 1  dt 2 ti '2-halo term',\
     onehalo  u 1:c w l lw 3 lc 1  dt 3 ti '1-halo term',\
     full     u 1:c w l lw 3 lc 1  dt 1 ti 'Full'#,\
     f(x)           w l lw 2 lc -1 dt 2 ti tits(ns)

set format x '10^{%T}'
if(print==0) set xlabel 'l'
if(print==1) set xlabel ''.ell.''

unset log y
set yrange [0:1]
set format y
set mytics
if(print==0) set ylabel 'C_{n,halo}(l) / C_{full}(l)'
if(print==1) set ylabel 'C_{n,halo}('.ell.') / C_{full}('.ell.')'

plot '<paste '.twohalo.' '.full.'' u 1:($3/$6) w l lw 3 lc 1  dt 2 noti,\
     '<paste '.onehalo.' '.full.'' u 1:($3/$6) w l lw 3 lc 1  dt 3 noti

unset multiplot
