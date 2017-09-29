unset multiplot
reset

if(print==0) set term aqua dashed
if(print==1) set term post enh col dashed dl .5 font ',10'; set output 'allpower.eps'

#set size square

data(sim,type1,type2)=sprintf('/Users/Mead/Physics/cosmo-OWLS/power/N400/%s_%s_%s_power.dat',sim,type1,type2)
hmpk(i,j)=sprintf('cosmo-OWLS/data/power_%i%i.dat',i,j)
dmonly='cosmo-OWLS/data/DMONLY.dat'

#All different simualtions
mod0='DMONLY'
mod1='REF'
mod2='NOCOOL_UVB'
mod3='AGN'
mod4='AGN_Theat_8p5'
mod5='AGN_Theat_8p7'

#Set the comparison model
mod='AGN'

#All different fields
thing0='overdensity_grid'
thing1='overdensity_grid_DM'
thing2='overdensity_grid_gas'
thing3='overdensity_grid_stars'

#Set colours
col0=0
col1=1
col2=2
col3=3
col4=4
col5=5

#Cosmological parameters
om_m=0.272
om_b=0.0455

bf=om_m/om_b
cf=om_m/(om_m-om_b)

bf2=bf**2.
cf2=cf**2.

set xlabel 'k / (h Mpc^{-1})'
set format x
set log x
set xrange [1e-2:1e1]

set log y
set yrange [1e-3:1.e3]
set format y '10^{%T}'
set ylabel '{/Symbol D}^2(k)'
set mytics 10

set title 'Comparison of '.mod.' simulation to DMONLY at z = 0'

set key bottom right

#Columns for simulation power
c=2
L=3

#Columns for halo-model power
d=5
M=5

#A small number
small=1e-22

plot small w l lw 3 dt 1 lc -1 ti 'Autospectra',\
     small w l lw 3 dt 2 lc -1 ti 'Cross with matter',\
     data(mod,thing0,thing0) u 1:(column(c)) w p pt 7 lc col1 noti,\
     data(mod,thing1,thing1) u 1:(cf2*column(c)) w p pt 7 lc col2 noti,\
     data(mod,thing2,thing2) u 1:(bf2*column(c)) w p pt 7 lc col3 noti,\
     data(mod,thing3,thing3) u 1:(bf2*column(c)) w p pt 7 lc col4 noti,\
     data(mod,thing0,thing1) u 1:(cf*column(c)) w p pt 6 lc col2 noti,\
     data(mod,thing0,thing2) u 1:(bf*column(c)) w p pt 6 lc col3 noti,\
     data(mod,thing0,thing3) u 1:(bf*column(c)) w p pt 6 lc col4 noti,\
     hmpk(0,0) u 1:(column(d)) w l lw 3 dt 1 lc col1 ti 'All matter',\
     hmpk(1,1) u 1:(cf2*column(d)) w l lw 3 dt 1 lc col2 ti 'CDM',\
     hmpk(2,2) u 1:(bf2*column(d)) w l lw 3 dt 1 lc col3 ti 'Gas',\
     hmpk(3,3) u 1:(bf2*column(d)) w l lw 3 dt 1 lc col4 ti 'Stars',\
     hmpk(0,1) u 1:(cf*column(d)) w l lw 3 dt 2 lc col2 noti,\
     hmpk(0,2) u 1:(bf*column(d)) w l lw 3 dt 2 lc col3 noti,\
     hmpk(0,3) u 1:(bf*column(d)) w l lw 3 dt 2 lc col4 noti

unset multiplot
