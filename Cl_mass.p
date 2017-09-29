reset

cmmi='/Users/Mead/Fonts/cmmi10.pfb'
cmsy='/Users/Mead/Fonts/cmsy10.pfb'

if(print==0){set term aqua}
if(print==1){set term post enh col sol fontfile cmmi fontfile cmsy; set output 'Cl_mass.eps'}

ell='{/cmmi10 \140}'
msun='{/cmsy10 \014}'

cl(m1,m2)=sprintf('data/mass_%i_%i_cl_full.dat',m1,m2)
#tits(m1,m2)=sprintf('10^{%2i} - 10^{%2i}',m1,m2)
#tits2(m)=sprintf('Minimum -> 10^{%2i}',m)

icumulative=1

set log x
if(print==0) set xlabel 'l'
if(print==1) set xlabel ''.ell.''
set mxtics 10
set xrange [1:1e5]
set format x '10^{%T}'

set log y
if(print==0) set ylabel 'l(l+1)C_{i,j}}(l) / 2{/Symbol p}'
if(print==1) set ylabel ''.ell.'('.ell.'+1)C_{i,j}}('.ell.') / 2{/Symbol p}'
set format y '10^{%T}'
set mytics 10
#set yrange [1e-13:1e-8]

mmin=1e11
mmax=1e16
nm=6
set palette defined (1 'pink', 2 'dark-red')
set log cb
set cbrange [mmin:mmax]
set format cb '10^{%T}'
if(print==0) set cblabel 'Integration maximum mass [M_'.msun.' / h]'
if(print==1) set cblabel 'Integration maximum mass [M_'.msun.' / h]'

set key top left

set title 'Cross correlation build-up as a function of halo mass'

if(icumulative==0){
plot for [i=10:15] cl(i,i+1) u 1:3 w l lw 3 lc i ti tits(i,i+1),\
     'data/cl_full.dat' u 1:3 w l lw 3 lc -1 ti 'Total'
}

if(icumulative==1){
plot for [i=1:nm] cl(7,i+10) u 1:3:(exp(log(mmin)+(log(mmax)-log(mmin))*real(i-1)/real(nm-1))) w l lw 3 lc palette noti,\
     'data/cl_full.dat' u 1:3 w l lw 3 lc -1 ti 'Total'
}

