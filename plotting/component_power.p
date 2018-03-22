unset multiplot
reset

if(!exists('print')){print=0}
if(print==0){set term aqua dashed}
if(print==1){set term post enh col font ',10'}

#Initial white space
print ''

#Compare to
print 'icomp = 1: Compare to cosmo-OWLS'
print 'icomp = 2: Compare to BAHAMAS'
if(!exists('icomp')){icomp=2}
print 'icomp = '.icomp.''
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

#Files for comparison
if(icomp==1){data(mod,field1,field2)=sprintf('/Users/Mead/Physics/cosmo-OWLS/power/N800/%s_%s_%s_power.dat',mod,field1,field2)}
if(icomp==2){data(sim,snap,field1,field2)=sprintf('/Users/Mead/Physics/BAHAMAS/power/M1024/%s_nu0_L400N1024_WMAP9_%s_%s_%s_power.dat',sim,snap,field1,field2)}
hmpk(z,i,j)=sprintf('hydro/power_z%1.1f_%i%i.dat',z,i,j)

#x axis
klab='k / h Mpc^{-1}'
kmin=1e-2
kmax=1e1
set xrange [kmin*1.1:kmax/1.1]
#set xrange [*:*]
set log x
set xlabel klab

#y axis
plab='{/Symbol D}_{i,j}^2(k)'
pmin=1e-7; pmax=1e3 #Suitable for matter spectra
if(field1 eq 'pressure'){pmin=pmin/1e3; pmax=pmax/1e3}
if(field2 eq 'pressure'){pmin=pmin/1e3; pmax=pmax/1e3}
set yrange [pmin*1.1:pmax/1.1]
#set yrange [*:*]
set log y
set format y '10^{%T}'
set ylabel plab

if(!exists('mod'))  {mod='AGN_TUNED'}
if(!exists('field1')){field1='all'}
if(!exists('field2')){field2='all'}

if(field1 eq 'all')     {i1=0}
if(field1 eq 'dm')      {i1=1}
if(field1 eq 'gas')     {i1=2}
if(field1 eq 'stars')   {i1=3}
if(field1 eq 'pressure'){i1=6}

if(field2 eq 'all')     {i2=0}
if(field2 eq 'dm')      {i2=1}
if(field2 eq 'gas')     {i2=2}
if(field2 eq 'stars')   {i2=3}
if(field2 eq 'pressure'){i2=6}

print 'Model *mod*: ', mod
print ''
print 'Filed types:'
print 'all - Matter'
print 'dm - CDM'
print 'gas - Gas'
print 'stars - Stars'
print 'pressure - Pressure'
print 'field 1 *field1*: ', i1, ' ', field1
print 'field 2 *field2*: ', i2, ' ', field2
print ''

print 'Simulation file: ', data(mod,snap,field1,field2)
print 'Halo-model file: ', hmpk(z,i1,i2)
print ''

#set title 'Simulation: '.mod.' || Power: '.field1.' x '.field2.' || z = '.z.''
#title_function(mod,field1,field2,z)=sprintf('Simulation: %s || Power: %s x %s || z = %1.1f',mod,field1,field2,z)
#set title title_function(mod,field1,field2,z)

unset key

#plot data(mod,snap,field1,field2) u 1:2 w p lc 1,\
#     hmpk(z,i1,i2) u 1:3 w l lc -1 dt 2 lw 3,\
#     hmpk(z,i1,i2) u 1:4 w l lc -1 dt 3 lw 3,\
#     hmpk(z,i1,i2) u 1:5 w l lc -1 dt 1 lw 3

title_function(mod,field1,field2)=sprintf('Simulation: %s || Power: %s x %s',mod,field1,field2)
#set title title_function(mod,field1,field2)

outfile(field1,field2)=sprintf('component_power_%s_%s.eps',field1,field2)
if(print==1){set output outfile(field1,field2); print 'output: ', outfile(field1,field2); print ''}

y2=0.9
y1=0.10
my=(y2+y1)/2

x1=0.10
x2=0.98
mx=(x2+x1)/2

set multiplot layout 2,2

do for [i=1:4]{

if(i==1){z=0.0; snap='snap32'; set xlabel ''; set format x ''; set ylabel plab; set format y '10^{%T}';
set tmargin at screen y2; set lmargin at screen x1; set rmargin at screen mx; set bmargin at screen my;
set label title_function(mod,field1,field2) at screen 0.37,0.95}

if(i==2){z=0.5; snap='snap28'; set xlabel ''; set format x ''; set ylabel ''; set format y '';
set tmargin at screen y2; set lmargin at screen mx; set rmargin at screen x2; set bmargin at screen my}

if(i==3){z=1.0; snap='snap26'; set xlabel klab; set format x; set ylabel plab; set format y '10^{%T}';
set tmargin at screen my; set lmargin at screen x1; set rmargin at screen mx; set bmargin at screen y1}

if(i==4){z=2.0; snap='snap22'; set xlabel klab; set format x; set ylabel ''; set format y '';
set tmargin at screen my; set lmargin at screen mx; set rmargin at screen x2; set bmargin at screen y1}

zlab(z)=sprintf('z = %1.1f', z)
set label zlab(z) at graph 0.05,0.9

plot data(mod,snap,field1,field2) u 1:2 w p lc 1,\
     hmpk(z,i1,i2) u 1:3 w l lc -1 dt 2 lw 3,\
     hmpk(z,i1,i2) u 1:4 w l lc -1 dt 3 lw 3,\
     hmpk(z,i1,i2) u 1:5 w l lc -1 dt 1 lw 3

unset label

}

unset multiplot
