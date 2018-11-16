reset

if(!exists('print')) {print=0}
if(print==0) {set term aqua dashed dl 1 font ',16'}
if(print==1) {set term post enh col sol font ',16'}

power(z,f1,f2)=sprintf('data/power_z%1.1f_%d%d.dat',z,f1,f2)
sim(mod,snap,f1,f2)=sprintf('/Users/Mead/Physics/BAHAMAS/power/M1536/%s_L400N1024_WMAP9_%s_%s_%s_power.dat',mod,snap,f1,f2)

fields="'all' 'dm' 'gas' 'stars' 'epressure'"
field_names="'all matter' 'dark matter' 'gas' 'stars' 'electron pressure'"
array ifield[5]
ifield[1]=0
ifield[2]=1
ifield[3]=2
ifield[4]=3
ifield[5]=6

snaps="'snap32' 'snap28' 'snap26' 'snap22'"
array zs[4]
zs[1]=0.0
zs[2]=0.5
zs[3]=1.0
zs[4]=2.0
z_names="'z = 0.0' 'z = 0.5' 'z = 1.0' 'z = 2.0'"
iz=1

print ''

if(!exists('iplot')){iplot=1}
print 'iplot = 1: Halo model compared to simulation'
print 'iplot = 2: Halo model compared to simulation with individual terms'
print 'iplot = 3: Halo model compared to simulations: CDM'
print 'iplot = 4: Halo model compared to simulations: gas'
print 'iplot = 5: Halo model compared to simulations: stars'
print 'iplot = ', iplot
print ''

kmin=3e-3
kmax=1e2
set log x
set xlabel 'k / h Mpc^{-1}'
set xrange [kmin:kmax]

dmin=1e-6
dmax=1e3
set log y
set ylabel '{/Symbol D}_{uv}^2(k)'
set format y '10^{%T}'
set yrange [dmin:dmax]

mod='AGN_TUNED_nu0'

set key top left

set label word(z_names,iz) at screen 0.9,0.2

if(iplot==1){

if(print==1){set output 'plots/power_talk_1.eps'}

plot power(zs[iz],ifield[1],ifield[1]) u 1:5 w l lc 1 dt 1 lw 3 ti 'Halo model',\
     sim(mod,word(snaps,iz),word(fields,1),word(fields,1)) u 1:($2-$3):5 w e lc -1 pt 7 ps .5 ti 'Simulation'

}

if(iplot==2){

if(print==1){set output 'plots/power_talk_2.eps'}

plot power(zs[iz],ifield[1],ifield[1]) u 1:5 w l lc 1 dt 1 lw 3 ti 'Halo model',\
     power(zs[iz],ifield[1],ifield[1]) u 1:3 w l lc 1 dt 2 lw 3 ti 'Two-halo term',\
     power(zs[iz],ifield[1],ifield[1]) u 1:4 w l lc 1 dt 3 lw 3 ti 'One-halo term',\
     sim(mod,word(snaps,iz),word(fields,1),word(fields,1)) u 1:($2-$3):5 w e lc -1 pt 7 ps .5 ti 'Simulation'

}

if(iplot==3 || iplot==4 || iplot==5){

if(print==1){
if(iplot==3){set output 'plots/power_talk_3.eps'}
if(iplot==4){set output 'plots/power_talk_4.eps'}
if(iplot==5){set output 'plots/power_talk_5.eps'}
}

if(iplot==3){c=2}
if(iplot==4){c=3}
if(iplot==5){c=4}

plot power(zs[iz],ifield[1],ifield[1]) u 1:5 w l lc 1 dt 1 lw 3 ti ''.word(field_names,1).'-'.word(field_names,1).'',\
     power(zs[iz],ifield[c],ifield[c]) u 1:5 w l lc c dt 1 lw 3 ti ''.word(field_names,c).'-'.word(field_names,c).'',\
     sim(mod,word(snaps,iz),word(fields,1),word(fields,1)) u 1:($2-$3):5 w e lc -1 pt 7 ps .5 ti 'Simulation',\
     sim(mod,word(snaps,iz),word(fields,c),word(fields,c)) u 1:($2-$3):5 w e lc -1 pt 7 ps .5 noti

}
