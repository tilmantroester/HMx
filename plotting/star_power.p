reset

owl(mod,type1,type2)=sprintf('/Users/Mead/Physics/cosmo-OWLS/power/N800/%s_%s_%s_power.dat',mod,type1,type2)
hmpk(i,j)=sprintf('data/power_%d%d.dat',i,j)

kmin=1e-2
kmax=1e1
set xrange [kmin:kmax]
set log x

set yrange [*:*]
set log y
set format y '10^{%T}'

mod='AGN'
type1='stars'
type2='stars'

plot owl(mod,type1,type2) u 1:2 w p,\
     hmpk(3,3) u 1:5 w l lw 3
