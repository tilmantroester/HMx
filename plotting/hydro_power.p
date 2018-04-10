unset multiplot
reset

if(!exists('print')){print=0}
if(print==0) {set term aqua dashed}
if(print==1) {set term post enh col dashed dl .5 font ',10'}

print ''

#Simulations to compare against
#1 - cosmo-OWLS
#2 - BAHAMAS
if(!exists('icomp')){icomp=2}
print 'icomp = 1: Compare to cosmo-OWLS'
print 'icomp = 2: Compare to BAHAMAS'
print 'icomp = 3: Generic hydro, no comparison'
print 'icomp = '.icomp.''
#print 'Variable *icomp* (1-2): '.icomp.''
if(icomp==1){sims='cosmo-OWLS'; Om_m=0.272; Om_b=0.0455}
if(icomp==2){sims='BAHAMAS'; Om_m=0.2793; Om_b=0.0463}
if(icomp==3){sims=''; Om_m=0.3; Om_b=0.05}
Om_c=Om_m-Om_b
print ''

#Redshift
if(!exists('iz')){iz=1}
print 'iz = 1: z = 0.0'
print 'iz = 2: z = 0.5'
print 'iz = 3: z = 1.0'
print 'iz = 4: z = 2.0'
print 'iz = '.iz.''
if(iz==1){z=0.0; snap='snap32'}
if(iz==2){z=0.5; snap='snap28'}
if(iz==3){z=1.0; snap='snap26'}
if(iz==4){z=2.0; snap='snap22'}
print 'z = ', z
print ''

#Plot to make
if(!exists('iplot')){iplot=1}
print 'iplot = 1: Power spectrum plot'
print 'iplot = 2: Power spectrum ratio plot'
print 'iplot = 3: Power spectrum suppression plot'
print 'iplot = 4: Power spectrum residual plot'
print 'iplot = '.iplot.''
#if(iplot==1){print 'Making power spectrum plot'}
#if(iplot==2){print 'Making power spectrum ratio plot'}
#if(iplot==3){print 'Making power spectrum suppression plot'}
#if(iplot==4){print 'Making power spectrum residual plot'}
print ''

#File names - cosmo-OWLS
if(icomp==1){
data(sim,type1,type2)=sprintf('/Users/Mead/Physics/cosmo-OWLS/power/N800/%s_%s_%s_power.dat',sim,type1,type2)
hmpk(sim,i,j)=sprintf('cosmo-OWLS/power_%s_%i%i.dat',sim,i,j)
}

#File names - BAHAMAS
if(icomp==2){
data(sim,snap,type1,type2)=sprintf('/Users/Mead/Physics/BAHAMAS/power/M1024/%s_nu0_L400N1024_WMAP9_%s_%s_%s_power.dat',sim,snap,type1,type2)
hmpk(sim,z,i,j)=sprintf('BAHAMAS/power_%s_z%1.1f_%i%i.dat',sim,z,i,j)
name(sim,z)=sprintf('BAHAMAS comarison of %s at z = %1.1f', sim, z)
}

if(icomp==3){
hmpk(sim,z,i,j)=sprintf('hydro/power_z%1.1f_%i%i.dat',z,i,j)
hmdm(z)=sprintf('hydro/power_z%1.1f.dat',z)
}

#Columns for simulation power
c=2
L=4

#Columns for halo-model power
d=5
M=5

#cosmo-OWLS simulation names
if(icomp==1){
hm_names="'DMONLY' 'REF' 'NOCOOL' 'AGN' 'AGN8p5' 'AGN8p7'"
owls_names="'DMONLY' 'REF' 'NOCOOL_UVB' 'AGN' 'AGN_Theat_8p5' 'AGN_Theat_8p7'"
}

#BAHAMAS simulation names
if(icomp==2){
hm_names="'DMONLY' 'AGN-lo' 'AGN-hi' 'AGN'"
owls_names="'DMONLY_2fluid' 'AGN_7p6' 'AGN_8p0' 'AGN_TUNED'"
}

if(icomp==3){
hm_names="''"
owls_names="''"
}

#Set the comparison model
if(!exists('nsim')){nsim=4}
hm_name=word(hm_names,nsim)
print 'Variable *nsim* '.nsim.''
print 'Simuation name: '.hm_name.''
if(icomp==1 || icomp==2){
owls_name=word(owls_names,nsim)
print 'Simuation file: '.owls_name.''
}
print ''

#All different fields
thing0='all'
thing1='dm'
thing2='gas'
thing3='stars'
thing6='pressure'

if(icomp==1 || icomp==2){print 'Example simulation file: ', data(owls_name,snap,thing0,thing0)}
print 'Example halo-model file: ', hmpk(hm_name,z,0,0)
print ''

#Set the files to compare against (DMONLY)
if(icomp==1){owl0='DMONLY'; hm0=hmpk('DMONLY',z,0,0)}
if(icomp==2){owl0='DMONLY_2fluid'; hm0=hmpk('DMONLY',z,0,0)}
if(icomp==3){owl0=''; hm0=hmdm(z)}

#Set colours
col0=0
col1=1
col2=2
col3=3
col4=4
col5=5
col6=6

#k range
if(icomp==1 || icomp==2){kmin=1e-2; kmax=1e1}
if(icomp==3){kmin=1e-2; kmax=1e2}
set xlabel 'k / h Mpc^{-1}'
set format x
set log x
set xrange [kmin:kmax]

#Delta^2(k) range
pmin=1e-5
pmax=1e3
set log y
set yrange [pmin:pmax]
set format y '10^{%T}'
set ylabel '{/Symbol D}^2(k)'
set mytics 10

set title name(owls_name,z)

if(iplot==1){

#set title 'Comparison of '.sims.' '.hm_name.' simulation to halo-model predictions'

set key bottom right

plot NaN w l lw 3 dt 1 lc -1 ti 'Autospectra',\
     NaN w l lw 3 dt 2 lc -1 ti 'Cross with matter',\
     data(owls_name,snap,thing0,thing0) u 1:(column(c)) w p pt 7 lc col1 noti,\
     data(owls_name,snap,thing1,thing1) u 1:(column(c)) w p pt 7 lc col2 noti,\
     data(owls_name,snap,thing2,thing2) u 1:(column(c)) w p pt 7 lc col3 noti,\
     data(owls_name,snap,thing3,thing3) u 1:(column(c)) w p pt 7 lc col4 noti,\
     data(owls_name,snap,thing0,thing1) u 1:(column(c)) w p pt 6 lc col2 noti,\
     data(owls_name,snap,thing0,thing2) u 1:(column(c)) w p pt 6 lc col3 noti,\
     data(owls_name,snap,thing0,thing3) u 1:(column(c)) w p pt 6 lc col4 noti,\
     hmpk(hm_name,z,0,0) u 1:(column(d)) w l lw 3 dt 1 lc col1 ti 'All matter',\
     hmpk(hm_name,z,1,1) u 1:(column(d)) w l lw 3 dt 1 lc col2 ti 'CDM',\
     hmpk(hm_name,z,2,2) u 1:(column(d)) w l lw 3 dt 1 lc col3 ti 'Gas',\
     hmpk(hm_name,z,3,3) u 1:(column(d)) w l lw 3 dt 1 lc col4 ti 'Stars',\
     hmpk(hm_name,z,0,1) u 1:(column(d)) w l lw 3 dt 2 lc col2 noti,\
     hmpk(hm_name,z,0,2) u 1:(column(d)) w l lw 3 dt 2 lc col3 noti,\
     hmpk(hm_name,z,0,3) u 1:(column(d)) w l lw 3 dt 2 lc col4 noti

}

if(iplot==2){

set multiplot layout 1,2

#Set the ratio axis
rmin=2e-3
rmax=1.5
set log y
set yrange [rmin:rmax]
set ylabel 'P_{OWL}/P_{DMONLY}'
set mytics 10

plot 1 w l lt -1 noti,\
     NaN w l lw 3 dt 1 lc -1 ti 'Autospectra',\
     NaN w l lw 3 dt 2 lc -1 ti 'Cross with matter',\
     Om_b/Om_m      w l lc -1 dt 2 noti,\
     Om_c/Om_m      w l lc -1 dt 2 noti,\
     (Om_b/Om_m)**2 w l lc -1 dt 2 noti,\
     (Om_c/Om_m)**2 w l lc -1 dt 2 noti,\
     '<paste '.data(owls_name,snap,thing0,thing0).' '.data(owl0,snap,thing0,thing0).'' u 1:(column(c)/column(c+L)) w p pt 7 lc col1 noti,\
     '<paste '.data(owls_name,snap,thing1,thing1).' '.data(owl0,snap,thing0,thing0).'' u 1:(column(c)/column(c+L)) w p pt 7 lc col2 noti,\
     '<paste '.data(owls_name,snap,thing2,thing2).' '.data(owl0,snap,thing0,thing0).'' u 1:(column(c)/column(c+L)) w p pt 7 lc col3 noti,\
     '<paste '.data(owls_name,snap,thing3,thing3).' '.data(owl0,snap,thing0,thing0).'' u 1:(column(c)/column(c+L)) w p pt 7 lc col4 noti,\
     '<paste '.data(owls_name,snap,thing0,thing0).' '.data(owl0,snap,thing0,thing0).'' u 1:(column(c)/column(c+L)) w p pt 6 lc col1 noti,\
     '<paste '.data(owls_name,snap,thing0,thing1).' '.data(owl0,snap,thing0,thing0).'' u 1:(column(c)/column(c+L)) w p pt 6 lc col2 noti,\
     '<paste '.data(owls_name,snap,thing0,thing2).' '.data(owl0,snap,thing0,thing0).'' u 1:(column(c)/column(c+L)) w p pt 6 lc col3 noti,\
     '<paste '.data(owls_name,snap,thing0,thing3).' '.data(owl0,snap,thing0,thing0).'' u 1:(column(c)/column(c+L)) w p pt 6 lc col4 noti,\
     '<paste '.hmpk(hm_name,z,0,0).' '.hm0.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 1 lc col1 ti 'All matter',\
     '<paste '.hmpk(hm_name,z,1,1).' '.hm0.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 1 lc col2 ti 'CDM',\
     '<paste '.hmpk(hm_name,z,2,2).' '.hm0.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 1 lc col3 ti 'Gas',\
     '<paste '.hmpk(hm_name,z,3,3).' '.hm0.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 1 lc col4 ti 'Stars',\
     '<paste '.hmpk(hm_name,z,0,0).' '.hm0.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 2 lc col1 noti,\
     '<paste '.hmpk(hm_name,z,0,1).' '.hm0.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 2 lc col2 noti,\
     '<paste '.hmpk(hm_name,z,0,2).' '.hm0.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 2 lc col3 noti,\
     '<paste '.hmpk(hm_name,z,0,3).' '.hm0.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 2 lc col4 noti

rmin=1e-9
rmax=2.
set yrange [rmin:rmax]
set format y '10^{%T}'

plot 1 w l lt -1 noti,\
     NaN w l lw 3 dt 1 lc -1 ti 'Autospectra',\
     NaN w l lw 3 dt 2 lc -1 ti 'Cross with matter',\
     '<paste '.data(owls_name,snap,thing0,thing0).' '.data(owl0,snap,thing0,thing0).'' u 1:(column(c)/column(c+L)) w p pt 7 lc col1 noti,\
     '<paste '.data(owls_name,snap,thing6,thing6).' '.data(owl0,snap,thing0,thing0).'' u 1:(column(c)/column(c+L)) w p pt 7 lc col6 noti,\
     '<paste '.data(owls_name,snap,thing0,thing6).' '.data(owl0,snap,thing0,thing0).'' u 1:(column(c)/column(c+L)) w p pt 7 lc col6 noti,\
     '<paste '.hmpk(hm_name,z,0,0).' '.hm0.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 1 lc col1 ti 'All matter',\
     '<paste '.hmpk(hm_name,z,6,6).' '.hm0.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 1 lc col6 ti 'Pressure',\
     '<paste '.hmpk(hm_name,z,0,6).' '.hm0.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 2 lc col6 noti

unset multiplot

}

if(iplot==3){

rmin=0.6
rmax=1.4
unset log y
set format y
set yrange [rmin:rmax]
set ylabel 'P(k) / P_{DMONLY}(k)'

if(icomp==1 || icomp==2){
plot 1 w l lt -1 noti,\
     for [i=1:words(owls_names)] '<paste '.data(word(owls_names,i),snap,thing0,thing0).' '.data(owl0,snap,thing0,thing0).'' u 1:(column(c)/column(c+L)) w p pt 7 dt 1 lc i noti,\
     for [i=1:words(hm_names)]   '<paste '.hmpk(word(hm_names,i),z,0,0).' '.hm0.'' u 1:(column(d)/column(d+M)) w l lw 3 lc i ti word(hm_names,i)
}

if(icomp==3){
plot 1 w l lt -1 noti,\
     for [i=1:words(hm_names)] '<paste '.hmpk(word(hm_names,i),z,0,0).' '.hm0.'' u 1:(column(d)/column(d+M)) w l lw 3 lc i ti word(hm_names,i)
}

}

if(iplot==4){

#Delta^2(k) range
rmin=0.5
rmax=1.5
unset log y
set format y
set yrange [rmin:rmax]
set ylabel 'P_{HM}(k) / P_{OWLS}(k)'

plot NaN w l lw 2 dt 1 lc -1 ti 'Autospectra',\
     NaN w l lw 2 dt 2 lc -1 ti 'Cross with matter',\
     1 w l lt -1 noti,\
     '<paste '.data(owls_name,snap,thing0,thing0).' '.hmpk(hm_name,z,0,0).'' u 1:(column(L+d)/column(c)) w l lw 2 dt 1 lc col1 ti 'All matter',\
     '<paste '.data(owls_name,snap,thing1,thing1).' '.hmpk(hm_name,z,1,1).'' u 1:(column(L+d)/column(c)) w l lw 2 dt 1 lc col2 ti 'CDM',\
     '<paste '.data(owls_name,snap,thing2,thing2).' '.hmpk(hm_name,z,2,2).'' u 1:(column(L+d)/column(c)) w l lw 2 dt 1 lc col3 ti 'Gas',\
     '<paste '.data(owls_name,snap,thing3,thing3).' '.hmpk(hm_name,z,3,3).'' u 1:(column(L+d)/column(c)) w l lw 2 dt 1 lc col4 ti 'Stars',\
     '<paste '.data(owls_name,snap,thing6,thing6).' '.hmpk(hm_name,z,6,6).'' u 1:(column(L+d)/column(c)) w l lw 2 dt 1 lc col6 ti 'Pressure',\
     '<paste '.data(owls_name,snap,thing0,thing1).' '.hmpk(hm_name,z,0,1).'' u 1:(column(L+d)/column(c)) w l lw 2 dt 2 lc col2 noti,\
     '<paste '.data(owls_name,snap,thing0,thing2).' '.hmpk(hm_name,z,0,2).'' u 1:(column(L+d)/column(c)) w l lw 2 dt 2 lc col3 noti,\
     '<paste '.data(owls_name,snap,thing0,thing3).' '.hmpk(hm_name,z,0,3).'' u 1:(column(L+d)/column(c)) w l lw 2 dt 2 lc col4 noti,\
     '<paste '.data(owls_name,snap,thing0,thing6).' '.hmpk(hm_name,z,0,6).'' u 1:(column(L+d)/column(c)) w l lw 2 dt 2 lc col6 noti

}
