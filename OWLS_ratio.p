unset multiplot
reset

if(print==0) set term aqua dashed
if(print==1) set term post enh col dashed dl .5 font ',10'; set output 'allpower.eps'

set size square

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

#set lmargin 10
#set rmargin 2

set xlabel 'k / (h Mpc^{-1})'
set format x
set log x
set xrange [1e-2:1e1]

set log y
set yrange [2e-3:1.5e0]
#set yrange [0.5:1.1]
set ylabel 'P_{OWL}/P_{DMONLY}'
#set format y '10^{%T}'
set mytics 10

set title 'Comparison of '.mod.' simulation to DMONLY at z = 0'

set key outside left box

#Columns for simulation power
c=2
L=3

#Columns for halo-model power
d=5
M=5

#A small number
small=1e-22

plot 1 w l lt -1 noti,\
     small w l lw 3 dt 1 lc -1 ti 'Autospectra',\
     small w l lw 3 dt 2 lc -1 ti 'Cross with matter',\
     om_b/om_m w l lc -1 dt 2 noti,\
     (om_m-om_b)/om_m w l lc -1 dt 2 noti,\
     (om_b/om_m)**2. w l lc -1 dt 2 noti,\
     ((om_m-om_b)/om_m)**2. w l lc -1 dt 2 noti,\
     '<paste '.data(mod,thing0,thing0).' '.data(mod0,thing0,thing0).'' u 1:(column(c)/column(c+L)) w p pt 7 lc col1 noti,\
     '<paste '.data(mod,thing1,thing1).' '.data(mod0,thing0,thing0).'' u 1:(column(c)/column(c+L)) w p pt 7 lc col2 noti,\
     '<paste '.data(mod,thing2,thing2).' '.data(mod0,thing0,thing0).'' u 1:(column(c)/column(c+L)) w p pt 7 lc col3 noti,\
     '<paste '.data(mod,thing3,thing3).' '.data(mod0,thing0,thing0).'' u 1:(column(c)/column(c+L)) w p pt 7 lc col4 noti,\
     '<paste '.data(mod,thing0,thing0).' '.data(mod0,thing0,thing0).'' u 1:(column(c)/column(c+L)) w p pt 6 lc col1 noti,\
     '<paste '.data(mod,thing0,thing1).' '.data(mod0,thing0,thing0).'' u 1:(column(c)/column(c+L)) w p pt 6 lc col2 noti,\
     '<paste '.data(mod,thing0,thing2).' '.data(mod0,thing0,thing0).'' u 1:(column(c)/column(c+L)) w p pt 6 lc col3 noti,\
     '<paste '.data(mod,thing0,thing3).' '.data(mod0,thing0,thing0).'' u 1:(column(c)/column(c+L)) w p pt 6 lc col4 noti,\
     '<paste '.hmpk(0,0).' '.dmonly.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 1 lc col1 ti 'All matter',\
     '<paste '.hmpk(1,1).' '.dmonly.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 1 lc col2 ti 'CDM',\
     '<paste '.hmpk(2,2).' '.dmonly.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 1 lc col3 ti 'Gas',\
     '<paste '.hmpk(3,3).' '.dmonly.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 1 lc col4 ti 'Stars',\
     '<paste '.hmpk(0,0).' '.dmonly.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 2 lc col1 noti,\
     '<paste '.hmpk(0,1).' '.dmonly.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 2 lc col2 noti,\
     '<paste '.hmpk(0,2).' '.dmonly.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 2 lc col3 noti,\
     '<paste '.hmpk(0,3).' '.dmonly.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 2 lc col4 noti

#Checked that the auto-spectra residual was not simply the square of the cross-spectra residual
#'<paste '.data(mod,thing0,thing2).' '.data(mod0,thing0,thing0).'' u 1:((column(c)/column(c+L))**2.) w p pt 6 lc -1 noti,\
#'<paste '.hmpk(0,2).' '.dmonly.'' u 1:(column(d)/column(d+M))**2. w l lw 3 dt 2 lc -1 noti

unset multiplot
