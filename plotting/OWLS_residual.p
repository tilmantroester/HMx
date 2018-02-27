unset multiplot
reset

if(!exists('print')){print=0}
if(print==0) set term aqua dashed
if(print==1) set term post enh col dashed dl .5 font ',10'; set output 'allpower.eps'

#File names
data(sim,type1,type2)=sprintf('/Users/Mead/Physics/cosmo-OWLS/power/N800/%s_%s_%s_power.dat',sim,type1,type2)
hmpk(sim,i,j)=sprintf('cosmo-OWLS/data/power_%s_%i%i.dat',sim,i,j)

#All different simualtions (names is my names, while OWLS are the actual simulation names)
names="'DMONLY' 'REF' 'NOCOOL' 'AGN' 'AGN8p5' 'AGN8p7'"
owls="'DMONLY' 'REF' 'NOCOOL_UVB' 'AGN' 'AGN_Theat_8p5' 'AGN_Theat_8p7'"

#Set the comparison model
n=4
owl=word(owls,n)
mod=word(names,n)

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
rmin=0.5
rmax=1.5
set yrange [rmin:rmax]
#set format y '10^{%T}'
set ylabel 'P_{HM}(k) / P_{OWLS}(k)'
#set mytics 10

set title 'Residual of halo-model to '.owl.' simulation'

set key bottom right

#Columns for simulation power
c=2
L=4

#Columns for halo-model power
d=5
M=5

plot NaN w l lw 2 dt 1 lc -1 ti 'Autospectra',\
     NaN w l lw 2 dt 2 lc -1 ti 'Cross with matter',\
     1 w l lt -1 noti,\
     '<paste '.data(owl,thing0,thing0).' '.hmpk(mod,0,0).'' u 1:(column(L+d)/column(c)) w l lw 2 dt 1 lc col1 ti 'All matter',\
     '<paste '.data(owl,thing1,thing1).' '.hmpk(mod,1,1).'' u 1:(column(L+d)/column(c)) w l lw 2 dt 1 lc col2 ti 'CDM',\
     '<paste '.data(owl,thing2,thing2).' '.hmpk(mod,2,2).'' u 1:(column(L+d)/column(c)) w l lw 2 dt 1 lc col3 ti 'Gas',\
     '<paste '.data(owl,thing3,thing3).' '.hmpk(mod,3,3).'' u 1:(column(L+d)/column(c)) w l lw 2 dt 1 lc col4 ti 'Stars',\
     '<paste '.data(owl,thing6,thing6).' '.hmpk(mod,6,6).'' u 1:(column(L+d)/column(c)) w l lw 2 dt 1 lc col6 ti 'Pressure',\
     '<paste '.data(owl,thing0,thing1).' '.hmpk(mod,0,1).'' u 1:(column(L+d)/column(c)) w l lw 2 dt 2 lc col2 noti,\
     '<paste '.data(owl,thing0,thing2).' '.hmpk(mod,0,2).'' u 1:(column(L+d)/column(c)) w l lw 2 dt 2 lc col3 noti,\
     '<paste '.data(owl,thing0,thing3).' '.hmpk(mod,0,3).'' u 1:(column(L+d)/column(c)) w l lw 2 dt 2 lc col4 noti,\
     '<paste '.data(owl,thing0,thing6).' '.hmpk(mod,0,6).'' u 1:(column(L+d)/column(c)) w l lw 2 dt 2 lc col6 noti

unset multiplot
