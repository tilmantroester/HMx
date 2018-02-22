unset multiplot
reset

if(!exists("print")){print=0}
if(print==0) set term aqua dashed
if(print==1) set term post enh col dashed dl .5 font ',10'; set output 'allpower.eps'

set size square

data(sim,type1,type2)=sprintf('/Users/Mead/Physics/cosmo-OWLS/power/N400/%s_%s_%s_power.dat',sim,type1,type2)
hmpk(i,j)=sprintf('data/power_%i%i.dat',i,j)
dmonly='data/power.dat'

#Set the comparison model
if(!exists("sim")){sim=4; print 'Setting sim: ', sim}
simulation_names='DMONLY REF NOCOOL_UVB AGN AGN_Theat_8p5 AGN_Theat_8p7'
simulation_titles="'DMONLY' 'REF' 'NO COOL' 'AGN' 'AGN 8.5' 'AGN 8.7'"
mod0=word(simulation_names,1)
mod=word(simulation_names,sim)
tit=word(simulation_names,sim)

print ''
print 'Comparing to simulation: ', mod
print ''

#All different fields
thing0='all'
thing1='dm'
thing2='gas'
thing3='stars'
thing4='pressure'

#Set colours
col0=0
col1=1
col2=2
col3=3
col4=4
col5=6

#Cosmological parameters
Om_m=0.272
Om_b=0.0455

#set lmargin 10
#set rmargin 2

kmin=1e-2
kmax=1e1
set xlabel 'k / h Mpc^{-1}'
set format x
set log x
set xrange [kmin:kmax]

rmin=2e-3
rmax=1.5
set log y
set yrange [rmin:rmax]
set ylabel 'P_{OWL}/P_{DMONLY}'
set mytics 10

set title 'Comparison of '.word(simulation_titles,sim).' simulation to DMONLY at z = 0'

#set key outside left box
set key bottom right

#Columns for simulation power
c=2
L=4

#Columns for halo-model power
d=5
M=5

set multiplot layout 1,2

plot 1 w l lt -1 noti,\
     NaN w l lw 3 dt 1 lc -1 ti 'Autospectra',\
     NaN w l lw 3 dt 2 lc -1 ti 'Cross with matter',\
     Om_b/Om_m             w l lc -1 dt 2 noti,\
     (Om_m-Om_b)/Om_m      w l lc -1 dt 2 noti,\
     (Om_b/Om_m)**2        w l lc -1 dt 2 noti,\
     ((Om_m-Om_b)/Om_m)**2 w l lc -1 dt 2 noti,\
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

p1=1.
p2=p1**2

rmin=1e-9
rmax=2.
set yrange [rmin:rmax]
set format y '10^{%T}'

plot 1 w l lt -1 noti,\
     NaN w l lw 3 dt 1 lc -1 ti 'Autospectra',\
     NaN w l lw 3 dt 2 lc -1 ti 'Cross with matter',\
     '<paste '.data(mod,thing0,thing0).' '.data(mod0,thing0,thing0).'' u 1:(column(c)/column(c+L)) w p pt 7 lc col1 noti,\
     '<paste '.data(mod,thing4,thing4).' '.data(mod0,thing0,thing0).'' u 1:(p1*column(c)/column(c+L)) w p pt 7 lc col5 noti,\
     '<paste '.data(mod,thing0,thing4).' '.data(mod0,thing0,thing0).'' u 1:(p2*column(c)/column(c+L)) w p pt 7 lc col5 noti,\
     '<paste '.hmpk(0,0).' '.dmonly.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 1 lc col1 ti 'All matter',\
     '<paste '.hmpk(6,6).' '.dmonly.'' u 1:(p1*column(d)/column(d+M)) w l lw 3 dt 1 lc col5 ti 'Pressure',\
     '<paste '.hmpk(0,6).' '.dmonly.'' u 1:(p2*column(d)/column(d+M)) w l lw 3 dt 2 lc col5 noti

unset multiplot

#Checked that the auto-spectra residual was not simply the square of the cross-spectra residual
#'<paste '.data(mod,thing0,thing2).' '.data(mod0,thing0,thing0).'' u 1:((column(c)/column(c+L))**2.) w p pt 6 lc -1 noti,\
#'<paste '.hmpk(0,2).' '.dmonly.'' u 1:(column(d)/column(d+M))**2. w l lw 3 dt 2 lc -1 noti

unset multiplot
