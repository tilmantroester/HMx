unset multiplot
reset

if(!exists('print')){print=0}
if(print==0) set term aqua dashed
if(print==1) set term post enh col dashed dl .5 font ',10'#; set output 'allpower.eps'

print ''

#Simulations to compare against
#1 - cosmo-OWLS
#2 - BAHAMAS
if(!exists('icomp')){icomp=1}
print 'Variable *icomp* (1-2): '.icomp.''
if(icomp==1){print 'Comparing to cosmo-OWLS simulations'; Om_m=0.272; Om_b=0.0455}
if(icomp==2){print 'Comparing to BAHAMAS simulations'; Om_m=0.2793; Om_b=0.0463}
Om_c=Om_m-Om_b
print ''

#Plot to make
#1 - Power spectrum
#2 - Power spectrum ratio
#3 - Power suppression
#4 - Power residual
if(!exists('iplot')){iplot=1}
print 'Variable *iplot* (1-4): '.iplot.''
if(iplot==1){print 'Making power spectrum plot'}
if(iplot==2){print 'Making power spectrum ratio plot'}
if(iplot==3){print 'Making power spectrum suppression plot'}
if(iplot==4){print 'Making power spectrum residual plot'}
print ''

#File names - cosmo-OWLS
if(icomp==1){
data(sim,type1,type2)=sprintf('/Users/Mead/Physics/cosmo-OWLS/power/N800/%s_%s_%s_power.dat',sim,type1,type2)
hmpk(sim,i,j)=sprintf('cosmo-OWLS/data/power_%s_%i%i.dat',sim,i,j)
}

#File names - BAHAMAS
if(icomp==2){
data(sim,type1,type2)=sprintf('/Users/Mead/Physics/BAHAMAS/power/M512/%s_nu0_L400N1024_WMAP9_snap32_%s_%s_power.dat',sim,type1,type2)
hmpk(sim,i,j)=sprintf('BAHAMAS/data/power_%s_%i%i.dat',sim,i,j)
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
hm_names="'DMONLY' 'AGN' 'AGN-lo' 'AGN-hi'"
owls_names="'DMONLY_2fluid' 'AGN_TUNED' 'AGN_7p6' 'AGN_8p0'"
}

#Set the comparison model
if(!exists('nsim')){nsim=4}
hm_name=word(hm_names,nsim)
owls_name=word(owls_names,nsim)
print 'Variable *nsim* '.nsim.''
print 'Simuation name: '.hm_name.''
print 'Simuation file: '.owls_name.''
print ''

#Set the files to compare against (DMONLY)
if(icomp==1){owl0='DMONLY'}
if(icomp==2){owl0='DMONLY_2fluid'}
hm0=hmpk('DMONLY',0,0)

#Snapshot
#22: z=2.0
#26: z=1.0
#38: z=0.5
#32: z=0.0
#zs=32

#All different fields
thing0='all'
thing1='dm'
thing2='gas'
thing3='stars'
thing6='pressure'

#Set colours
col0=0
col1=1
col2=2
col3=3
col4=4
col5=5
col6=6

#k range
kmin=1e-2
kmax=1e1
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

if(iplot==1){

set title 'Comparison of '.hm_name.' simulation to halo-model predictions'

set key bottom right

plot NaN w l lw 3 dt 1 lc -1 ti 'Autospectra',\
     NaN w l lw 3 dt 2 lc -1 ti 'Cross with matter',\
     data(owls_name,thing0,thing0) u 1:(column(c)) w p pt 7 lc col1 noti,\
     data(owls_name,thing1,thing1) u 1:(column(c)) w p pt 7 lc col2 noti,\
     data(owls_name,thing2,thing2) u 1:(column(c)) w p pt 7 lc col3 noti,\
     data(owls_name,thing3,thing3) u 1:(column(c)) w p pt 7 lc col4 noti,\
     data(owls_name,thing0,thing1) u 1:(column(c)) w p pt 6 lc col2 noti,\
     data(owls_name,thing0,thing2) u 1:(column(c)) w p pt 6 lc col3 noti,\
     data(owls_name,thing0,thing3) u 1:(column(c)) w p pt 6 lc col4 noti,\
     hmpk(hm_name,0,0) u 1:(column(d)) w l lw 3 dt 1 lc col1 ti 'All matter',\
     hmpk(hm_name,1,1) u 1:(column(d)) w l lw 3 dt 1 lc col2 ti 'CDM',\
     hmpk(hm_name,2,2) u 1:(column(d)) w l lw 3 dt 1 lc col3 ti 'Gas',\
     hmpk(hm_name,3,3) u 1:(column(d)) w l lw 3 dt 1 lc col4 ti 'Stars',\
     hmpk(hm_name,0,1) u 1:(column(d)) w l lw 3 dt 2 lc col2 noti,\
     hmpk(hm_name,0,2) u 1:(column(d)) w l lw 3 dt 2 lc col3 noti,\
     hmpk(hm_name,0,3) u 1:(column(d)) w l lw 3 dt 2 lc col4 noti

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
     '<paste '.data(owls_name,thing0,thing0).' '.data(owl0,thing0,thing0).'' u 1:(column(c)/column(c+L)) w p pt 7 lc col1 noti,\
     '<paste '.data(owls_name,thing1,thing1).' '.data(owl0,thing0,thing0).'' u 1:(column(c)/column(c+L)) w p pt 7 lc col2 noti,\
     '<paste '.data(owls_name,thing2,thing2).' '.data(owl0,thing0,thing0).'' u 1:(column(c)/column(c+L)) w p pt 7 lc col3 noti,\
     '<paste '.data(owls_name,thing3,thing3).' '.data(owl0,thing0,thing0).'' u 1:(column(c)/column(c+L)) w p pt 7 lc col4 noti,\
     '<paste '.data(owls_name,thing0,thing0).' '.data(owl0,thing0,thing0).'' u 1:(column(c)/column(c+L)) w p pt 6 lc col1 noti,\
     '<paste '.data(owls_name,thing0,thing1).' '.data(owl0,thing0,thing0).'' u 1:(column(c)/column(c+L)) w p pt 6 lc col2 noti,\
     '<paste '.data(owls_name,thing0,thing2).' '.data(owl0,thing0,thing0).'' u 1:(column(c)/column(c+L)) w p pt 6 lc col3 noti,\
     '<paste '.data(owls_name,thing0,thing3).' '.data(owl0,thing0,thing0).'' u 1:(column(c)/column(c+L)) w p pt 6 lc col4 noti,\
     '<paste '.hmpk(hm_name,0,0).' '.hm0.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 1 lc col1 ti 'All matter',\
     '<paste '.hmpk(hm_name,1,1).' '.hm0.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 1 lc col2 ti 'CDM',\
     '<paste '.hmpk(hm_name,2,2).' '.hm0.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 1 lc col3 ti 'Gas',\
     '<paste '.hmpk(hm_name,3,3).' '.hm0.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 1 lc col4 ti 'Stars',\
     '<paste '.hmpk(hm_name,0,0).' '.hm0.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 2 lc col1 noti,\
     '<paste '.hmpk(hm_name,0,1).' '.hm0.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 2 lc col2 noti,\
     '<paste '.hmpk(hm_name,0,2).' '.hm0.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 2 lc col3 noti,\
     '<paste '.hmpk(hm_name,0,3).' '.hm0.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 2 lc col4 noti

rmin=1e-9
rmax=2.
set yrange [rmin:rmax]
set format y '10^{%T}'

plot 1 w l lt -1 noti,\
     NaN w l lw 3 dt 1 lc -1 ti 'Autospectra',\
     NaN w l lw 3 dt 2 lc -1 ti 'Cross with matter',\
     '<paste '.data(owls_name,thing0,thing0).' '.data(owl0,thing0,thing0).'' u 1:(column(c)/column(c+L)) w p pt 7 lc col1 noti,\
     '<paste '.data(owls_name,thing6,thing6).' '.data(owl0,thing0,thing0).'' u 1:(column(c)/column(c+L)) w p pt 7 lc col6 noti,\
     '<paste '.data(owls_name,thing0,thing6).' '.data(owl0,thing0,thing0).'' u 1:(column(c)/column(c+L)) w p pt 7 lc col6 noti,\
     '<paste '.hmpk(hm_name,0,0).' '.hm0.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 1 lc col1 ti 'All matter',\
     '<paste '.hmpk(hm_name,6,6).' '.hm0.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 1 lc col6 ti 'Pressure',\
     '<paste '.hmpk(hm_name,0,6).' '.hm0.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 2 lc col6 noti

unset multiplot

}

if(iplot==3){

rmin=0.6
rmax=1.4
unset log y
set format y
set yrange [rmin:rmax]
set ylabel 'P(k) / P_{DMONLY}(k)'

plot 1 w l lt -1 noti,\
     for [i=1:words(owls_names)] '<paste '.data(word(owls_names,i),thing0,thing0).' '.data(owl0,thing0,thing0).'' u 1:(column(c)/column(c+L)) w p pt 7 dt 1 lc i noti,\
     for [i=1:words(hm_names)] '<paste '.hmpk(word(hm_names,i),0,0).' '.hm0.'' u 1:(column(d)/column(d+M)) w l lw 3 lc i ti word(hm_names,i)

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
     '<paste '.data(owls_name,thing0,thing0).' '.hmpk(hm_name,0,0).'' u 1:(column(L+d)/column(c)) w l lw 2 dt 1 lc col1 ti 'All matter',\
     '<paste '.data(owls_name,thing1,thing1).' '.hmpk(hm_name,1,1).'' u 1:(column(L+d)/column(c)) w l lw 2 dt 1 lc col2 ti 'CDM',\
     '<paste '.data(owls_name,thing2,thing2).' '.hmpk(hm_name,2,2).'' u 1:(column(L+d)/column(c)) w l lw 2 dt 1 lc col3 ti 'Gas',\
     '<paste '.data(owls_name,thing3,thing3).' '.hmpk(hm_name,3,3).'' u 1:(column(L+d)/column(c)) w l lw 2 dt 1 lc col4 ti 'Stars',\
     '<paste '.data(owls_name,thing6,thing6).' '.hmpk(hm_name,6,6).'' u 1:(column(L+d)/column(c)) w l lw 2 dt 1 lc col6 ti 'Pressure',\
     '<paste '.data(owls_name,thing0,thing1).' '.hmpk(hm_name,0,1).'' u 1:(column(L+d)/column(c)) w l lw 2 dt 2 lc col2 noti,\
     '<paste '.data(owls_name,thing0,thing2).' '.hmpk(hm_name,0,2).'' u 1:(column(L+d)/column(c)) w l lw 2 dt 2 lc col3 noti,\
     '<paste '.data(owls_name,thing0,thing3).' '.hmpk(hm_name,0,3).'' u 1:(column(L+d)/column(c)) w l lw 2 dt 2 lc col4 noti,\
     '<paste '.data(owls_name,thing0,thing6).' '.hmpk(hm_name,0,6).'' u 1:(column(L+d)/column(c)) w l lw 2 dt 2 lc col6 noti

}
