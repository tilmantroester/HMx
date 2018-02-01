unset multiplot
reset

if(print==0) set term aqua dashed
if(print==1) set term post enh col dashed dl .5 font ',10'; set output 'ratio_pressure.eps'

set size square

N=800

data(N,sim,type1,type2)=sprintf('/Users/Mead/Physics/cosmo-OWLS/power/N%i/%s_%s_%s_power.dat',N,sim,type1,type2)
hmpk(i,j)=sprintf('data/power_%s%s.dat',i,j)
dmonly='data/power.dat'

#All different simualtions
#mod0='DMONLY'
#mod1='REF'
#mod2='NOCOOL_UVB'
#mod3='AGN'
#mod4='AGN_Theat_8p5'
#mod5='AGN_Theat_8p7'
models='DMONLY REF NOCOOL_UVB AGN AGN_Theat_8p5 AGN_Theat_8p7'
titles="'DMONLY' 'REF' 'NO COOL' 'AGN' 'AGN 8.5' 'AGN 8.7'"

#Set the comparison model
mod=word(models,3)

#All different fields
#thing0='overdensity_grid'
#thing1='overdensity_grid_DM'
#thing2='overdensity_grid_gas'
#thing3='overdensity_grid_stars'
#thing4='pressure_grid'
fields="'all' 'dm' 'gas' 'stars' 'pressure'"

#Set colours
col0=-1
col1=1
col2=2
col3=3
col4=4
col5=5

#Cosmological parameters
om_m=0.272
om_b=0.0455

#Conversion of CGS to SI pressure units (Pascals per barye)
#cgs=0.1

#Multiplication factor for pressure
fac=1.e-4

#set lmargin 10
#set rmargin 2

kmin=1e-2
kmax=1e1
set xlabel 'k / (h Mpc^{-1})'
set format x
set log x
set xrange [kmin:kmax]

rmin=1e-5
rmax=1e1
set log y
#set yrange [1e-3:1e1]
set yrange [*:*]
set ylabel 'P_{ij,OWL}/P_{DMONLY}'
set format y '10^{%T}'
set mytics 10

set title 'Comparison of '.mod.' simulation to DMONLY at z = 0 with pressure'

set key outside left box

#Columns for simulation power
c=2
L=4

#Columns for halo-model power
d=5
M=5

plot 1 w l lt -1 noti,\
     '<paste '.data(N,mod,'all','all').' '.data(N,'DMONLY','all','all').'' u 1:(column(c)/column(c+L)) w p pt 7 lc col1 noti,\
     '<paste '.data(N,mod,'all','pressure').' '.data(N,'DMONLY','all','all').'' u 1:(fac*column(c)/column(c+L)) w p pt 7 lc col2 noti,\
     '<paste '.data(N,mod,'pressure','pressure').' '.data(N,'DMONLY','all','all').'' u 1:(fac*fac*column(c)/column(c+L)) w p pt 7 lc col3 noti,\
     '<paste '.hmpk('d','d').' '.dmonly.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 1 lc col1 ti '{/Symbol d}{/Symbol d}',\
     '<paste '.hmpk('d','p').' '.dmonly.'' u 1:(fac*column(d)/column(d+M)) w l lw 3 dt 1 lc col2 ti '{/Symbol d}p',\
     '<paste '.hmpk('p','p').' '.dmonly.'' u 1:(fac*fac*column(d)/column(d+M)) w l lw 3 dt 1 lc col3 ti 'pp',\
     '<paste '.hmpk('d','d').' '.dmonly.'' u 1:(column(d-2)/column(d-2+M)) w l lw 3 dt 2 lc col1 noti,\
     '<paste '.hmpk('d','p').' '.dmonly.'' u 1:(fac*column(d-2)/column(d-2+M)) w l lw 3 dt 2 lc col2 noti,\
     '<paste '.hmpk('p','p').' '.dmonly.'' u 1:(fac*fac*column(d-2)/column(d-2+M)) w l lw 3 dt 2 lc col3 noti,\
     '<paste '.hmpk('d','d').' '.dmonly.'' u 1:(column(d-1)/column(d-1+M)) w l lw 3 dt 3 lc col1 noti,\
     '<paste '.hmpk('d','p').' '.dmonly.'' u 1:(fac*column(d-1)/column(d-1+M)) w l lw 3 dt 3 lc col2 noti,\
     '<paste '.hmpk('p','p').' '.dmonly.'' u 1:(fac*fac*column(d-1)/column(d-1+M)) w l lw 3 dt 3 lc col3 noti

unset multiplot
