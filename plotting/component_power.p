reset

set term aqua dashed

owl(mod,type1,type2)=sprintf('/Users/Mead/Physics/cosmo-OWLS/power/N800/%s_%s_%s_power.dat',mod,type1,type2)
hmpk(i,j)=sprintf('data/power_%d%d.dat',i,j)

kmin=1e-2
kmax=1e1
set xrange [kmin:kmax]
#set xrange [*:*]
set log x
set xlabel 'k / h Mpc^{-1}'

pmin=1e-5
pmax=1e1
#set yrange [pmin:pmax]
set yrange [*:*]
set log y
set format y '10^{%T}'
set ylabel '{/Symbol D}_{i,j}^2(k)'

#mod='AGN'
#type1='stars'
#type2='stars'

if(!exists('mod')){mod='AGN'}
if(!exists('type1')){type1='all'}
if(!exists('type2')){type2='all'}

if(type1 eq 'all')     {i1=0}
if(type1 eq 'dm')      {i1=1}
if(type1 eq 'gas')     {i1=2}
if(type1 eq 'stars')   {i1=3}
if(type1 eq 'pressure'){i1=6}

if(type2 eq 'all')     {i2=0}
if(type2 eq 'dm')      {i2=1}
if(type2 eq 'gas')     {i2=2}
if(type2 eq 'stars')   {i2=3}
if(type2 eq 'pressure'){i2=6}

print ''
print 'Model *mod*: ', mod
print 'field 1 *type1*: ', type1
print 'field 2 *type2*: ', type2
print 'cosmo-OWLS file: ', owl(mod,type1,type2)
print 'halo-model file: ', hmpk(i1,i2)
print ''

set title 'Power: '.type1.' x '.type2.''

unset key

plot owl(mod,type1,type2) u 1:2 w p lc -1,\
     hmpk(i1,i2) u 1:3 w l lc -1 dt 2 lw 3,\
     hmpk(i1,i2) u 1:4 w l lc -1 dt 3 lw 3,\
     hmpk(i1,i2) u 1:5 w l lc -1 dt 1 lw 3
