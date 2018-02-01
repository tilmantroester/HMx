unset multiplot
reset

if(!exists('print')) {print=0}
if(print==0) {set term aqua dashed}
if(print==1) {set term post enh col dashed dl .5 font ',10'; set output 'pressure_power.eps'}

#Mesh size
N=800

#field(N,sim,type1,type2)=sprintf('/Users/Mead/Physics/cosmo-OWLS/power_fields/N%i/%s_%s_%s_power.dat',N,sim,type1,type2)
field(N,sim,type1,type2)=sprintf('/Users/Mead/Physics/cosmo-OWLS/power/N%i/%s_%s_%s_power.dat',N,sim,type1,type2)
hmpk(i,j)=sprintf('data/power_%s%s.dat',i,j)
dmonly='data/power.dat'

#Different simualtions
models='DMONLY REF NOCOOL_UVB AGN AGN_Theat_8p5 AGN_Theat_8p7'
titles="'DMONLY' 'REF' 'NO COOL' 'AGN' 'AGN 8.5' 'AGN 8.7'"

#Set the simulation here
m=6
mod=word(models,m)
print 'Comparing to simulation: '.mod.''

#All different fields
#fields="'overdensity_grid' 'overdensity_grid_DM' 'overdensity_grid_gas' overdensity_grid_stars' 'pressure_grid'"
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

#Factor to multiply by for pressure to get nice comparison
fm=1.
fh=1.
print 'Simulation factor: ', fm
print 'Model factor: ', fh
print ''

kmin=1e-2
kmax=1e1
set log x
set xrange [kmin:kmax]

set title 'Halo model power spectra for matter and pressure fields compared to simulation: '.word(titles,m).''

#Columns for simulation power
c=2
#L=3
L=4

#Columns for halo-model power
d=5
M=5

set lmargin 10
set rmargin 2

set multiplot layout 2,1

set key top left

set xlabel ''
set format x ''

pmin=1e-11
pmax=1e3
set log y
set yrange [pmin:pmax]
set ylabel '{/Symbol D}^2(k)'
set format y '10^{%T}'
set mytics 10

plot field(N,mod,word(fields,1),word(fields,1)) u 1:(column(c)) w p pt 7 lc col1 noti,\
     field(N,mod,word(fields,1),word(fields,5)) u 1:(fm*column(c)) w p pt 7 lc col2 noti,\
     field(N,mod,word(fields,5),word(fields,5)) u 1:(fm*fm*column(c)) w p pt 7 lc col3 noti,\
     dmonly u 1:2           w l lw 3 dt 2 lc col0 ti 'Linear',\
     dmonly u 1:(column(d)) w l lw 3 lc col0 ti 'DMONLY',\
     hmpk('d','d') u 1:(column(d))         w l dt 1 lw 3 lc col1 ti '{/Symbol d}{/Symbol d}',\
     hmpk('d','d') u 1:(column(d-1))       w l dt 2 lw 3 lc col1 noti,\
     hmpk('d','d') u 1:(column(d-2))       w l dt 3 lw 3 lc col1 noti,\
     hmpk('d','p') u 1:(fh*column(d))      w l dt 1 lw 3 lc col2 ti '{/Symbol d}p',\
     hmpk('d','p') u 1:(fh*column(d-1))    w l dt 2 lw 3 lc col2 noti,\
     hmpk('d','p') u 1:(fh*column(d-2))    w l dt 3 lw 3 lc col2 noti,\
     hmpk('p','p') u 1:(fh*fh*column(d))   w l dt 1 lw 3 lc col3 ti 'pp',\
     hmpk('p','p') u 1:(fh*fh*column(d-1)) w l dt 2 lw 3 lc col3 noti,\
     hmpk('p','p') u 1:(fh*fh*column(d-2)) w l dt 3 lw 3 lc col3 noti

unset title

set format x
set xlabel 'k / (h Mpc^{-1})'

pmin=1e-8
pmax=1e1
set yrange [pmin:pmax]
set ylabel 'P(k) / P_{DMONLY}(k)'

plot 1 w l lt -1 noti,\
     '<paste '.field(N,mod,word(fields,1),word(fields,1)).' '.field(N,'DMONLY',word(fields,1),word(fields,1)).'' u 1:(column(c)/column(c+L)) w p pt 7 lc col1 noti,\
     '<paste '.field(N,mod,word(fields,1),word(fields,5)).' '.field(N,'DMONLY',word(fields,1),word(fields,1)).'' u 1:(fm*column(c)/column(c+L)) w p pt 7 lc col2 noti,\
     '<paste '.field(N,mod,word(fields,5),word(fields,5)).' '.field(N,'DMONLY',word(fields,1),word(fields,1)).'' u 1:(fm*fm*column(c)/column(c+L)) w p pt 7 lc col3 noti,\
     '<paste '.hmpk('d','d').' '.dmonly.'' u 1:(column(d)/column(d+M))           w l lw 3 dt 1 lc col1 noti '{/Symbol d}{/Symbol d}',\
     '<paste '.hmpk('d','p').' '.dmonly.'' u 1:(fh*column(d)/column(d+M))        w l lw 3 dt 1 lc col2 noti '{/Symbol d}p',\
     '<paste '.hmpk('p','p').' '.dmonly.'' u 1:(fh*fh*column(d)/column(d+M))     w l lw 3 dt 1 lc col3 noti 'pp'

#One-halo and two-halo stuff (unnecessary?)
#'<paste '.hmpk('d','d').' '.dmonly.'' u 1:(column(d-2)/column(d-2+M))       w l lw 3 dt 2 lc col1 noti,\
#'<paste '.hmpk('d','p').' '.dmonly.'' u 1:(fh*column(d-2)/column(d-2+M))    w l lw 3 dt 2 lc col2 noti,\
#'<paste '.hmpk('p','p').' '.dmonly.'' u 1:(fh*fh*column(d-2)/column(d-2+M)) w l lw 3 dt 2 lc col3 noti,\
#'<paste '.hmpk('d','d').' '.dmonly.'' u 1:(column(d-1)/column(d-1+M))       w l lw 3 dt 3 lc col1 noti,\
#'<paste '.hmpk('d','p').' '.dmonly.'' u 1:(fh*column(d-1)/column(d-1+M))    w l lw 3 dt 3 lc col2 noti,\
#'<paste '.hmpk('p','p').' '.dmonly.'' u 1:(fh*fh*column(d-1)/column(d-1+M)) w l lw 3 dt 3 lc col3 noti

unset multiplot
