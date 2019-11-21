unset multiplot
reset

# Terminal options
if(!exists('print')){print=0}
if(print==0) {set term qt dashed dl 1 size 1000,800}
if(print==1) {set term post enh col dashed dl .5 font ',10'}

# Initial white space
print ''

# Plot to make
if(!exists('iplot')){iplot=10}
print 'iplot = 1:  Power spectrum plot'
print 'iplot = 2:  Power spectrum ratio plot'
print 'iplot = 3:  Power spectrum suppression plot'
print 'iplot = 4:  Power spectrum residual plot'
print 'iplot = 5:  Power spectrum components'
print 'iplot = 6:  Power spectrum of electron pressure'
print 'iplot = 7:  Power spectrum with k^1.5 units'
print 'iplot = 8:  PAPER: Combination of iplot=1 and 2'
print 'iplot = 9:  Response residual'
print 'iplot = 10: Combination of iplot=1 and 2'
print 'iplot = 11: Matter-electron pressure spectrum variance demonstration'
print 'iplot = 12: electron pressure-electron pressure spectrum variance demonstration'
print 'iplot = 13: Halo model compared to simulation'
print 'iplot = 14: Halo model comparison: matter-matter'
print 'iplot = 15: Halo model comparison: CDM-CDM'
print 'iplot = 16: Halo model comparison: gas-gas'
print 'iplot = 17: Halo model comparison: stars-stars'
print 'iplot = 18: Halo model comparison: electron pressure-electron pressure'
print 'iplot = 19: Halo model comparison: matter-electron pressure'
print 'iplot = 20: Power as a functon of AGN strength: matter-matter'
print 'iplot = 21: Power as a functon of AGN strength: CDM-CDM'
print 'iplot = 22: Power as a functon of AGN strength: gas-gas'
print 'iplot = 23: Power as a functon of AGN strength: stars-stars'
print 'iplot = 24: Power as a functon of AGN strength: matter-electron pressure'
print 'iplot = 25: Response as a functon of AGN strength: matter-matter'
print 'iplot = 26: Response as a functon of AGN strength: CDM-CDM'
print 'iplot = 27: Response as a functon of AGN strength: gas-gas'
print 'iplot = 28: Response as a functon of AGN strength: stars-stars'
print 'iplot = 29: Response as a functon of AGN strength: matter-electron pressure'
print 'iplot = 30: All power as a function of AGN strength'
print 'iplot = 31: All responses as a function of AGN strength'
print 'iplot = 32: PAPER: Combination of iplot=31 and 32'
print 'iplot = '.iplot.''
print ''

# Simulations to compare against
# 1 - cosmo-OWLS
# 2 - BAHAMAS
# 3 - Generic hydro comparison
if(!exists('icomp')){icomp=3}
if(iplot==8){icomp=3}
print 'icomp = 1: Compare to cosmo-OWLS (NO LONGER SUPPORTED)'
print 'icomp = 2: Compare to BAHAMAS'
print 'icomp = 3: Generic hydro, no comparison simulation'
print 'icomp = '.icomp.''
print ''

if(icomp==1){print 'Twat; icomp=1 does not work any more'; exit}
#if(icomp==1){sims='cosmo-OWLS'; Om_m=0.272; Om_b=0.0455}
#if(icomp==2){sims='BAHAMAS'}
#if(icomp==3){sims=''}

# Cosmological parameters (only used for plotting Om_b/Om_m lines)
Om_m=0.2793
Om_b=0.0463
Om_c=Om_m-Om_b
print 'Omega_m: ', Om_m
print 'Omega_b: ', Om_b
print 'Omega_c: ', Om_c
print ''

# Redshift
if(!exists('z')){z=0.0}
if(iplot==8){z=0.0}
if(z==0.0){snap='snap32'}
if(z==0.5){snap='snap28'}
if(z==1.0){snap='snap26'}
if(z==2.0){snap='snap22'}
snaps="'snap32' 'snap28' 'snap26' 'snap22'"
z_names="'z = 0.0' 'z = 0.5' 'z = 1.0' 'z = 2.0'"
array zs[4]
zs[1]=0.0
zs[2]=0.5
zs[3]=1.0
zs[4]=2.0
print 'Redshift: z = ', z
print ''

array gs[4]
gs[1]=1.000
gs[2]=0.779
gs[3]=0.619
gs[4]=0.428
do for [iz=1:4] {print 'Redshift: z, g(z): ', zs[iz], gs[iz]}
print ''

# File names - cosmo-OWLS
#if(icomp==1){
#data(sim,type1,type2)=sprintf('/Users/Mead/Physics/cosmo-OWLS/power/N800/%s_%s_%s_power.dat',sim,type1,type2)
#hmpk(sim,i,j)=sprintf('cosmo-OWLS/power_%s_%i%i.dat',sim,i,j)
#print 'Error: cosmo-OWLS comparison needs to be updated'
#exit
#}

# Simulation P(k) mesh-size options
if(!exists('mesh')){mesh=1024}
if(iplot==8){mesh=1024}
print 'Mesh size: mesh: ', mesh
print ''

# Simulation data files
plot_title_name_z(sim,z)=sprintf('BAHAMAS comarison of %s at z = %1.1f', sim, z)
plot_title_z(z)=sprintf('BAHAMAS comarison at z = %1.1f', z)
simpk(sim,mesh,snap,type1,type2)=sprintf('/Users/Mead/Physics/BAHAMAS/power/M%d/%s_L400N1024_WMAP9_%s_%s_%s_power.dat',mesh,sim,snap,type1,type2)
sim_dmonly='DMONLY_2fluid_nu0'

# File names - BAHAMAS
if(icomp==2){
hmpk(sim,z,i,j)=sprintf('data/power_%s_z%1.1f_%i%i.dat',sim,z,i,j)
hmpk_dmonly=hmpk('DMONLY',z,0,0)
}

# File names - generic hydro
if(icomp==3){
hmpk(sim,z,i,j)=sprintf('data/power_z%1.1f_%i%i.dat',z,i,j)
hmdm(z)=sprintf('data/power_z%1.1f.dat',z)
hmpk_dmonly=hmdm(z)
}

# Number of columns and pertinent column for simulation power
c=2
s=3
L=5

# Number of columns and pertinent column for halo-model power
d=5
M=5

#cosmo-OWLS simulation names
#if(icomp==1){
#hm_names="'DMONLY' 'REF' 'NOCOOL' 'AGN' 'AGN8p5' 'AGN8p7'"
#owls_names="'DMONLY' 'REF' 'NOCOOL_UVB' 'AGN' 'AGN_Theat_8p5' 'AGN_Theat_8p7'"
#}

# BAHAMAS simulation names
if(icomp==2){hmpk_names="'DMONLY' 'AGN-lo' 'AGN' 'AGN-hi'"; hmpk_cols="'black' 'dark-yellow' 'blue' 'dark-plum'"}
if(icomp==3){hmpk_names="''"; hmpk_cols="'black'"}
sims="'DMONLY_2fluid_nu0' 'AGN_7p6_nu0' 'AGN_TUNED_nu0' 'AGN_8p0_nu0' 'AGN_TUNED_nu0_v2' 'AGN_TUNED_nu0_v3' 'AGN_TUNED_nu0_v4'"
sim_names="'DMonly' 'AGN-lo' 'AGN' 'AGN-hi' 'AGN v2' 'AGN v3' 'AGN v4'"
sim_cols="'black' 'dark-yellow' 'blue' 'dark-plum' 'light-blue' 'light-blue' 'light-blue'"

sims_few="'DMONLY_2fluid_nu0' 'AGN_7p6_nu0' 'AGN_TUNED_nu0' 'AGN_8p0_nu0'"
sim_few_names="'DMonly' 'AGN-lo' 'AGN' 'AGN-hi'"
sim_few_cols="'black' 'dark-yellow' 'blue' 'dark-plum'"

# Set the comparison model
if(!exists('nsim')){nsim=3}       # Default to AGN
if(iplot==8 || iplot==11){nsim=3} # Default to AGN
hmpk_name=word(hmpk_names,nsim)
print 'Simulation model number: nsim: '.nsim.''
if(icomp==2) print 'Halo model power file name: '.hmpk_name.''
sim=word(sims,nsim)
print 'Simulation power file: '.sim.''
print ''

# All different fields for power spectra
fields="'all' 'dm' 'gas' 'stars' 'epressure'"
field_names="'matter' 'CDM' 'gas' 'stars' 'electron pressure'"

# Field integers
array ifield[5]
ifield[1]=0
ifield[2]=1
ifield[3]=2
ifield[4]=3
ifield[5]=6

# Field colours
array icol[5]
icol[1]=1
icol[2]=2
icol[3]=3
icol[4]=4
icol[5]=6

# Few fields
fields_few="'all' 'epressure'"
field_few_names="'matter' 'electron pressure [keV cm^{-3}]'"

# Few field integers
array ifield_few[2]
ifield_few[1]=0
ifield_few[2]=6

# Few field colours
array icol_few[2]
icol_few[1]=1
icol_few[2]=6

# Fractions to multiply pressure spectra by
array fac[2]
fac[1]=1.
fac[2]=1e3 # Units of pressure become [keV/cm^3]

# Write to screen
print 'Pressure field multiplied by: ', fac[2]
print ''

# Write useful things to screen
print 'Example simulation file: ', simpk(sim,mesh,snap,word(fields,1),word(fields,1))
print 'Example halo-model file: ', hmpk(hmpk_name,z,ifield[1],ifield[1])
print ''

# k range
if(icomp==1 || icomp==2){kmin=1e-2; kmax=2e1}
if(icomp==3){kmin=1e-2; kmax=2e1}
klab='k / h Mpc^{-1}'
set xlabel klab
set format x
set log x
set xrange [kmin:kmax]

# Delta^2(k) range
pmin=1e-7
pmax=1e3
set log y
set yrange [pmin:pmax]
set format y '10^{%T}'
set ylabel '{/Symbol D}_{uv}^2(k)'
set mytics 10

# Set the overall plot titles
set title plot_title_name_z(sim,z) noenh

## ##

if(iplot==1){

if(print==1){
outfile(name,z)=sprintf('%s_z%1.1f_power.eps',name,z)
set output outfile(sim,z)
print 'Outfile: ', outfile(sim,z)
print ''
}

set multiplot layout 1,2

unset title

set key bottom right

plot NaN w l lw 3 dt 1 lc -1 ti 'Autospectra',\
     NaN w l lw 3 dt 2 lc -1 ti 'Cross with matter',\
     for [i=1:4] simpk(sim,mesh,snap,word(fields,i),word(fields,i)) u 1:(column(c)-column(s)):5 w e pt 7 ps .5 lc icol[i] noti,\
     for [i=1:4] simpk(sim,mesh,snap,word(fields,1),word(fields,i)) u 1:(column(c)-column(s)):5 w e pt 6 ps .5 lc icol[i] noti,\
     for [i=1:4] hmpk(hmpk_name,z,ifield[i],ifield[i]) u 1:(column(d)) w l lw 3 dt 1 lc icol[i] ti word(field_names,i),\
     for [i=1:4] hmpk(hmpk_name,z,ifield[1],ifield[i]) u 1:(column(d)) w l lw 3 dt 2 lc icol[i] noti

plot for [i=1:2] NaN w l lw 3 dt 1 lc icol_few[i] ti word(field_few_names,i),\
     for [i=1:2] for [j=i:2]   simpk(sim,mesh,snap,word(fields_few,i),word(fields_few,j)) u 1:(fac[i]*fac[j]*(column(c)-column(s))):(fac[i]*fac[j]*column(5)) w e pt 7 ps .5 lc icol_few[i] noti,\
     for [i=1:2] for [j=i:2]   hmpk(hmpk_name,z,ifield_few[i],ifield_few[j]) u 1:(fac[i]*fac[j]*column(d)) w l lw 3 dt 1 lc icol_few[i] noti,\
     for [i=1:2] for [j=i+1:2] hmpk(hmpk_name,z,ifield_few[i],ifield_few[j]) u 1:(fac[i]*fac[j]*column(d)) w l lw 3 dt 2 lc icol_few[j] noti

unset multiplot

}

## ##

## ##

if(iplot==2){

if(print==1){
outfile(name,z)=sprintf('%s_z%1.1f_ratio.eps',name,z)
set output outfile(sim,z)
}

set multiplot layout 1,2

# Set the ratio axis
rmin=2e-3
rmax=1.5
set log y
set yrange [rmin:rmax]
set ylabel 'P_{OWL}/P_{DMONLY}'
set mytics 10

# left - matter repsonse
plot 1 w l lt -1 noti,\
     NaN w l lw 3 dt 1 lc -1 ti 'Autospectra',\
     NaN w l lw 3 dt 2 lc -1 ti 'Cross with matter',\
     Om_b/Om_m      w l lc -1 dt 2 noti,\
     Om_c/Om_m      w l lc -1 dt 2 noti,\
     (Om_b/Om_m)**2 w l lc -1 dt 2 noti,\
     (Om_c/Om_m)**2 w l lc -1 dt 2 noti,\
     for [i=1:4] '<paste '.simpk(sim,mesh,snap,word(fields,i),word(fields,i)).' '.simpk(sim_dmonly,mesh,snap,'all','all').'' u 1:((column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 7 lc icol[i] noti,\
     for [i=1:4] '<paste '.simpk(sim,mesh,snap,'all',word(fields,i)).' '.simpk(sim_dmonly,mesh,snap,'all','all').'' u 1:((column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 6 lc icol[i] noti,\
     for [i=1:4] '<paste '.hmpk(hmpk_name,z,ifield[i],ifield[i]).' '.hmpk_dmonly.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 1 lc icol[i] ti word(field_names,i),\
     for [i=1:4] '<paste '.hmpk(hmpk_name,z,ifield[1],ifield[i]).' '.hmpk_dmonly.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 2 lc icol[i] noti

rmin=1e-6
rmax=2.
set yrange [rmin:rmax]
set format y '10^{%T}'

# right - pressure response
plot 1 w l lt -1 noti,\
     for [i=1:2] for [j=i:2] '<paste '.simpk(sim,mesh,snap,word(fields_few,i),word(fields_few,j)).' '.simpk(sim_dmonly,mesh,snap,'all','all').'' u 1:(fac[i]*fac[j]*(column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 7 ps .5 lc icol_few[i] noti,\
     for [i=1:2] for [j=i:2]   '<paste '.hmpk(hmpk_name,z,ifield_few[i],ifield_few[j]).' '.hmpk_dmonly.'' u 1:(fac[i]*fac[j]*column(d)/column(d+M)) w l lw 3 dt 1 lc icol_few[i] noti,\
     for [i=1:2] for [j=i+1:2] '<paste '.hmpk(hmpk_name,z,ifield_few[i],ifield_few[j]).' '.hmpk_dmonly.'' u 1:(fac[i]*fac[j]*column(d)/column(d+M)) w l lw 3 dt 2 lc icol_few[j] noti

unset multiplot

}

## ##

## ##

if(iplot==3){

set title plot_title_z(z) noenh

if(print==1){
outfile='paper/suppression.eps'
set output outfile
}

kmin=0.011
kmax=9.9
set xrange [kmin:kmax]

rmin=0.76
rmax=1.06
unset log y
set format y
set yrange [rmin:rmax]

top=0.98
bot=0.10
lef=0.10
rig=0.98
miy=(top+bot)/2.
mix=(lef+rig)/2.

unset title

set multiplot layout 2,2

do for [i=1:4] {

set key bottom left

if(i==1){snap='snap32'; zlab='z = 0.0'; set format x ''; set xlabel ''; set format y; set ylabel 'P(k) / P_{DMONLY}(k)'}
if(i==2){snap='snap28'; zlab='z = 0.5'; set format x ''; set xlabel ''; set format y ''; set ylabel ''}
if(i==3){snap='snap26'; zlab='z = 1.0'; set format x; set xlabel 'k / h^{-1} Mpc'; set format y; set ylabel 'P(k) / P_{DMONLY}(k)'}
if(i==4){snap='snap22'; zlab='z = 2.0'; set format x; set xlabel 'k / h^{-1} Mpc'; set format y ''; set ylabel ''}

if(i==1){set tmargin at screen top; set bmargin at screen miy; set lmargin at screen lef; set rmargin at screen mix}
if(i==2){set tmargin at screen top; set bmargin at screen miy; set lmargin at screen mix; set rmargin at screen rig}
if(i==3){set tmargin at screen miy; set bmargin at screen bot; set lmargin at screen lef; set rmargin at screen mix}
if(i==4){set tmargin at screen miy; set bmargin at screen bot; set lmargin at screen mix; set rmargin at screen rig}

names=hmpk_names
if(i==2 || i==3 || i==4){names="''''''''"}

set label zlab at graph 0.1,0.9

if(icomp==1 || icomp==2){
plot 1 w l lt -1 noti,\
     for [i=1:words(sims_few)]   '<paste '.simpk(word(sims_few,i),mesh,snap,'all','all').' '.simpk(sim_dmonly,mesh,snap,'all','all').'' u 1:((column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 7 dt 1 lc rgb word(sim_few_cols,i) noti,\
     for [i=1:words(hmpk_names)] '<paste '.hmpk(word(hmpk_names,i),z,0,0).' '.hmpk_dmonly.'' u 1:(column(d)/column(d+M)) w l lw 3 lc rgb word(sim_few_cols,i) ti word(names,i)
}

if(icomp==3){
plot 1 w l lt -1 noti,\
     for [i=1:words(sims_few)]    '<paste '.simpk(word(sims,i),mesh,snap,'all','all').' '.simpk(sim_dmonly,mesh,snap,'all','all').'' u 1:((column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 7 dt 1 lc rgb word(sim_few_cols,i) noti,\
     for [i=1:words(hmpk_names)] '<paste '.hmpk(word(hmpk_names,i),z,0,0).' '.hmpk_dmonly.'' u 1:(column(d)/column(d+M)) w l lw 3 lc -1 ti word(names,i)
}

unset label

}

unset multiplot

}

## ##

## ##

if(iplot==4){

if(print==1){
outfile(name,z)=sprintf('%s_z%1.1f_residual.eps',name,z)
set output outfile(sim,z)
}

print 'Note that k ranges must align here!'
print ''
#if(icomp==3){print 'iplot=4 does not work with icomp=3 because the k axis do not align'; print ''; exit}

# Delta^2(k) range
rmin=0.5
rmax=1.5
unset log y
set format y
set yrange [rmin:rmax]
set ylabel 'P_{HM}(k) / P_{OWLS}(k)'

plot NaN w l lw 2 dt 1 lc -1 ti 'Autospectra',\
     NaN w l lw 2 dt 2 lc -1 ti 'Cross with matter',\
     1 w l lt -1 noti,\
     for [i=1:1] '<paste '.simpk(sim,mesh,snap,word(fields,i),word(fields,i)).' '.hmpk(hmpk_name,z,ifield[i],ifield[i]).'' u 1:(column(L+d)/(column(c)-column(s))) w l lw 2 dt 1 lc icol[i] ti word(field_names,i),\
     for [i=5:5] '<paste '.simpk(sim,mesh,snap,word(fields,1),word(fields,i)).' '.hmpk(hmpk_name,z,ifield[1],ifield[i]).'' u 1:(column(L+d)/(column(c)-column(s))) w l lw 2 dt 2 lc icol[i] noti

#'<paste '.simpk(sim,mesh,snap,word(fields,1),word(fields,1)).' '.hmpk(hmpk_name,z,ifield[1],ifield[1]).'' u 1:(column(L+d)/(column(c)-column(s))) w l lw 2 dt 1 lc icol[1] ti word(field_names,i),\
#'<paste '.simpk(sim,mesh,snap,word(fields,1),word(fields,5)).' '.hmpk(hmpk_name,z,ifield[1],ifield[5]).'' u 1:(column(L+d)/(column(c)-column(s))) w l lw 2 dt 2 lc icol[5] noti

}

## ##

## ##

if(iplot==5){

# Set fields if none are specified
if(!exists('field1')){field1='all'}
if(!exists('field2')){field2='all'}

# File name for output
if(print==1){
outfile(name,field1,field2)=sprintf('%s_%s_%s_components.eps',name,field1,field2)
set output outfile(sim,field1,field2)
}

# x axis
set xrange [kmin*1.1:kmax/1.1]
set log x
set xlabel klab

# y axis
plab='{/Symbol D}_{uv}^2(k)'
pmin=1e-7; pmax=1e3 # Suitable for matter spectra
if(field1 eq 'pressure'){pmin=pmin/1e3; pmax=pmax/1e3}
if(field2 eq 'pressure'){pmin=pmin/1e3; pmax=pmax/1e3}
set yrange [pmin*1.1:pmax/1.1]
set log y
set format y '10^{%T}'
set ylabel plab

# Set field 1
if(field1 eq 'all')     {i1=0}
if(field1 eq 'dm')      {i1=1}
if(field1 eq 'gas')     {i1=2}
if(field1 eq 'stars')   {i1=3}
if(field1 eq 'pressure'){i1=6}

# Set field 2
if(field2 eq 'all')     {i2=0}
if(field2 eq 'dm')      {i2=1}
if(field2 eq 'gas')     {i2=2}
if(field2 eq 'stars')   {i2=3}
if(field2 eq 'pressure'){i2=6}

# Print field information to the screen
print 'Field types:'
print 'all - Matter'
print 'dm - CDM'
print 'gas - Gas'
print 'stars - Stars'
print 'pressure - Pressure'
print 'field 1 *field1*: ', i1, ' ', field1
print 'field 2 *field2*: ', i2, ' ', field2
print ''

unset key

# Title function for plots
title_function(name,field1,field2)=sprintf('Simulation: %s || Power: %s x %s',name,field1,field2)

# y scale info for multiplot
y2=0.9
y1=0.10
my=(y2+y1)/2

# x scale info for multiplot
x1=0.10
x2=0.98
mx=(x2+x1)/2

# Get rid of the generic title
unset title

# Do the actual plotting
set multiplot layout 2,2

do for [i=1:4]{

if(i==1){zz=0.0; ss='snap32'; set xlabel ''; set format x ''; set ylabel plab; set format y '10^{%T}';
set tmargin at screen y2; set lmargin at screen x1; set rmargin at screen mx; set bmargin at screen my;
set label title_function(sim,field1,field2) at screen 0.37,0.95}

if(i==2){zz=0.5; ss='snap28'; set xlabel ''; set format x ''; set ylabel ''; set format y '';
set tmargin at screen y2; set lmargin at screen mx; set rmargin at screen x2; set bmargin at screen my}

if(i==3){zz=1.0; ss='snap26'; set xlabel klab; set format x; set ylabel plab; set format y '10^{%T}';
set tmargin at screen my; set lmargin at screen x1; set rmargin at screen mx; set bmargin at screen y1}

if(i==4){zz=2.0; ss='snap22'; set xlabel klab; set format x; set ylabel ''; set format y '';
set tmargin at screen my; set lmargin at screen mx; set rmargin at screen x2; set bmargin at screen y1}

zlab(z)=sprintf('z = %1.1f', z)
set label zlab(zz) at graph 0.05,0.9

plot simpk(sim,mesh,ss,field1,field2) u 1:(column(c)-column(s)):5 w e lc 1,\
     hmpk(hmpk_name,zz,i1,i2) u 1:3 w l lc -1 dt 2 lw 3,\
     hmpk(hmpk_name,zz,i1,i2) u 1:4 w l lc -1 dt 3 lw 3,\
     hmpk(hmpk_name,zz,i1,i2) u 1:5 w l lc -1 dt 1 lw 3

unset label

}

unset multiplot

}

## ##

## ##

if(iplot==6){

# Pressure spectrum

if(print==1){
outfile(name,z)=sprintf('%s_z%1.1f_pressure.eps',name,z)
set output outfile(sim,z)
}

set multiplot layout 2,1

unset xlabel
set format x ''

set key top left

plot for [i=1:2] NaN w l lw 3 lc icol_few[i] ti word(field_few_names,i),\
     for [i=1:2] for [j=i:2] simpk(sim,mesh,snap,word(fields_few,i),word(fields_few,j)) u 1:(fac[i]*fac[j]*(column(c)-column(s))):(fac[i]*fac[j]*column(5)) w e pt 7 lc icol_few[i] noti,\
     for [k=1:3] for [i=1:2] for [j=i:2] hmpk(hmpk_name,z,ifield_few[i],ifield_few[j]) u 1:(fac[i]*fac[j]*column(d-k+1)) w l lw 3 dt k lc icol_few[i] noti,\
     for [i=1:2] for [j=i+1:2] hmpk(hmpk_name,z,ifield_few[i],ifield_few[j]) u 1:(fac[i]*fac[j]*column(d)) w l lw 3 dt 2 lc icol_few[j] noti

unset title

set xlabel klab
set format x

# Set the y axis for P(k) / P_{DMONLY}(k)
rmin=1e-5
rmax=2
set yrange [rmin:rmax]
set ylabel 'P(k) / P_{DMONLY}(k)'

plot 1 w l lt -1 noti,\
     for [i=1:2] for [j=i:2] '<paste '.simpk(sim,mesh,snap,word(fields_few,i),word(fields_few,j)).' '.simpk(sim_dmonly,mesh,snap,'all','all').'' u 1:(fac[i]*fac[j]*(column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 7 lc icol_few[i] noti,\
     for [i=1:2] for [j=i:2]   '<paste '.hmpk(hmpk_name,z,ifield_few[i],ifield_few[j]).' '.hmpk_dmonly.'' u 1:(fac[i]*fac[j]*column(d)/column(d+M)) w l lw 3 dt 1 lc icol_few[i] noti,\
     for [i=1:2] for [j=i+1:2] '<paste '.hmpk(hmpk_name,z,ifield_few[i],ifield_few[j]).' '.hmpk_dmonly.'' u 1:(fac[i]*fac[j]*column(d)/column(d+M)) w l lw 3 dt 2 lc icol_few[j] noti

unset multiplot

}

## ##

## ##

if(iplot==7){

pow=1.5

# Delta^2(k)/k^1.5 range
pmin=1e-4
pmax=1e2
set log y
set yrange [pmin:pmax]
set format y '10^{%T}'
set ylabel '{/Symbol D}_{uv}^2(k) / [k / h^{-1} Mpc]^{1.5}'
set mytics 10

if(print==1){
outfile(name,z)=sprintf('%s_z%1.1f_power.eps',name,z)
set output outfile(sim,z)
}

set key bottom right

plot NaN w l lw 3 dt 1 lc -1 ti 'Autospectra',\
     NaN w l lw 3 dt 2 lc -1 ti 'Cross with matter',\
     for [i=1:4] simpk(sim,mesh,snap,word(fields,i),word(fields,i)) u 1:((column(c)-column(s))/(column(1)**pow)):($5/(column(1)**pow)) w e pt 7 lc icol[i] noti,\
     for [i=1:4] simpk(sim,mesh,snap,word(fields,1),word(fields,i)) u 1:((column(c)-column(s))/(column(1)**pow)):($5/(column(1)**pow)) w e pt 6 lc icol[i] noti,\
     for [i=1:4] hmpk(hmpk_name,z,ifield[i],ifield[i]) u 1:((column(d)/column(1)**pow)) w l lw 3 dt 1 lc icol[i] ti word(fields,i),\
     for [i=1:4] hmpk(hmpk_name,z,ifield[1],ifield[i]) u 1:((column(d)/column(1)**pow)) w l lw 3 dt 2 lc icol[i] ti word(fields,i)

}

## ##

## PAPER FIGURE: big ##

if(iplot==8 || iplot==10){

if(print==0){set term qt dashed size 1200,800}

if(print==1){
#outfile(name,z)=sprintf('%s_z%1.1f_power.eps',name,z)
outfile='paper/hydro.eps'
set output outfile
print 'Outfile: ', outfile
print ''
}

set multiplot layout 2,2

set label ''.word(sim_names,nsim).'; z = '.sprintf('%1.1f', z).'' at graph 0.05,0.9

set xlabel ''
set format x ''

pmin=1e-7
pmax=1e3
set yrange [pmin:pmax]

set key bottom right

unset title

plot NaN w l lw 3 dt 1 lc -1 ti 'Autospectra',\
     NaN w l lw 3 dt 2 lc -1 ti 'Cross with matter',\
     for [i=1:4] simpk(sim,mesh,snap,word(fields,i),word(fields,i)) u 1:(column(c)-column(s)):5 w e pt 7 ps .5 lc icol[i] noti,\
     for [i=1:4] simpk(sim,mesh,snap,word(fields,1),word(fields,i)) u 1:(column(c)-column(s)):5 w e pt 6 ps .5 lc icol[i] noti,\
     for [i=1:4] hmpk(hmpk_name,z,ifield[i],ifield[i]) u 1:(column(d)) w l lw 3 dt 1 lc icol[i] ti word(field_names,i),\
     for [i=1:4] hmpk(hmpk_name,z,ifield[1],ifield[i]) u 1:(column(d)) w l lw 3 dt 2 lc icol[i] noti

unset label

pmin=1e-5
pmax=1e3
set yrange [pmin:pmax]

plot for [i=1:2] NaN w l lw 3 dt 1 lc icol_few[i] ti word(field_few_names,i),\
     for [i=1:2] for [j=i:i] simpk(sim,mesh,snap,word(fields_few,i),word(fields_few,j)) u 1:(fac[i]*fac[j]*(column(c)-column(s))):(fac[i]*fac[j]*column(5)) w e pt 7 ps .5 lc icol_few[i] noti,\
     for [i=1:1] for [j=2:2] simpk(sim,mesh,snap,word(fields_few,i),word(fields_few,j)) u 1:(fac[i]*fac[j]*(column(c)-column(s))):(fac[i]*fac[j]*column(5)) w e pt 6 ps .5 lc icol_few[i] noti,\
     for [i=1:2] for [j=i:2] hmpk(hmpk_name,z,ifield_few[i],ifield_few[j]) u 1:(fac[i]*fac[j]*column(d)) w l lw 3 dt 1 lc icol_few[i] noti,\
     for [i=1:2] for [j=i:2] hmpk(hmpk_name,z,ifield_few[i],ifield_few[j]) u 1:(fac[i]*fac[j]*column(d)) w l lw 3 dt 2 lc icol_few[j] noti

set xlabel 'k / h Mpc^{-1}'
set format x

rmin=1e-4
rmax=2
set ylabel 'P_{uv}(k) / P_{no-hydro}(k)'
set yrange [rmin:rmax]

# Bottom left - matter response
plot 1 w l lt -1 noti,\
     Om_b/Om_m      w l lc -1 dt 2 noti,\
     Om_c/Om_m      w l lc -1 dt 2 noti,\
     (Om_b/Om_m)**2 w l lc -1 dt 2 noti,\
     (Om_c/Om_m)**2 w l lc -1 dt 2 noti,\
     for [i=1:4] '<paste '.simpk(sim,mesh,snap,word(fields,i),word(fields,i)).' '.simpk(sim_dmonly,mesh,snap,word(fields,1),word(fields,1)).'' u 1:((column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 7 ps .5 lc icol[i] noti,\
     for [i=1:4] '<paste '.simpk(sim,mesh,snap,word(fields,1),word(fields,i)).' '.simpk(sim_dmonly,mesh,snap,word(fields,1),word(fields,1)).'' u 1:((column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 6 ps .5 lc icol[i] noti,\
     for [i=1:4] '<paste '.hmpk(hmpk_name,z,ifield[i],ifield[i]).' '.hmpk_dmonly.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 1 lc icol[i] noti,\
     for [i=1:4] '<paste '.hmpk(hmpk_name,z,ifield[1],ifield[i]).' '.hmpk_dmonly.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 2 lc icol[i] noti

rmin=1e-2
rmax=2
set ylabel 'P_{uv}(k) / P_{no-hydro}(k)'
set yrange [rmin:rmax]

# Bottom right - pressure response
plot 1 w l lt -1 noti,\
     for [i=1:2] for [j=i:2] '<paste '.simpk(sim,mesh,snap,word(fields_few,i),word(fields_few,j)).' '.simpk(sim_dmonly,mesh,snap,'all','all').'' u 1:(fac[i]*fac[j]*(column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 7 ps .5 lc icol_few[i] noti,\
     for [i=1:2] for [j=i:2] '<paste '.hmpk(hmpk_name,z,ifield_few[i],ifield_few[j]).' '.hmpk_dmonly.'' u 1:(fac[i]*fac[j]*column(d)/column(d+M)) w l lw 3 dt 1 lc icol_few[i] noti,\
     for [i=1:2] for [j=i:2] '<paste '.hmpk(hmpk_name,z,ifield_few[i],ifield_few[j]).' '.hmpk_dmonly.'' u 1:(fac[i]*fac[j]*column(d)/column(d+M)) w l lw 3 dt 2 lc icol_few[j] noti

unset multiplot

}

## ##

## Response ratio ##

if(iplot==9){

# Set to do BAHAMAS comparison only
#icomp=2
print 'Note that k ranges mus align here!'
print ''

# Terminal commands
if(print==0){set term qt dashed size 1200,800}
if(print==1){set term post enh col dashed; set output 'response_ratio.eps'}

# Label position
labx=0.25
laby=0.9

# Plot limits in x
lef=0.08
rig=0.98
nx=3.
dx=(rig-lef)/nx

# k axis
kmin=0.013
kmax=9
set xrange [kmin:kmax]
set log x

# Plot limits in y
top=0.99
bot=0.08
ny=4.
dy=(top-bot)/ny

# ratio axis
rmin=0.5
rmax=1.5
unset log y
set yrange [rmin:rmax]

# No title for this plot
unset title

set key top left

# What the fuck is this?
isim=3
iz=2

set multiplot layout 4,3

do for [isim=2:4]{

if(isim==2){set lmargin at screen lef+0*dx; set rmargin at screen lef+1*dx}
if(isim==3){set lmargin at screen lef+1*dx; set rmargin at screen lef+2*dx}
if(isim==4){set lmargin at screen lef+2*dx; set rmargin at screen lef+3*dx}

if(isim==2){set format y; set ylabel 'R_{HM} / R_{sim}'}
if(isim==3 || isim==4){set format y ''; set ylabel ''}

do for [iz=1:4]{

if(icomp==2){
combi(isim,iz,i,j)=sprintf('<paste '.hmpk(word(hmpk_names,isim),zs[iz],ifield[i],ifield[j]).' '.hmpk(word(hmpk_names,1),zs[iz],ifield[1],ifield[1]).' '.simpk(word(sims,isim),mesh,word(snaps,iz),word(fields,i),word(fields,j)).' '.simpk(sim_dmonly,mesh,word(snaps,iz),'all','all').'',isim,iz,i,j)
}

if(icomp==3){
combi(isim,iz,i,j)=sprintf('<paste '.hmpk(word(hmpk_names,isim),zs[iz],ifield[i],ifield[j]).' '.hmdm(zs[iz]).' '.simpk(word(sims,isim),mesh,word(snaps,iz),word(fields,i),word(fields,j)).' '.simpk(sim_dmonly,mesh,word(snaps,iz),'all','all').'',isim,iz,i,j)
}

if(iz==1){set tmargin at screen top-0*dy; set bmargin at screen top-1*dy}
if(iz==2){set tmargin at screen top-1*dy; set bmargin at screen top-2*dy}
if(iz==3){set tmargin at screen top-2*dy; set bmargin at screen top-3*dy}
if(iz==4){set tmargin at screen top-3*dy; set bmargin at screen top-4*dy}

if(iz==1 || iz==2 || iz==3){set xlabel ''; set format x ''; unset key}
if(iz==4){set xlabel 'k / h^{-1} Mpc'; set format x}

set label ''.word(hmpk_names,isim).'; '.word(z_names,iz).'' at graph labx,laby

if(iz==1 && isim==2){set key; unset label}

plot 1 w l lt -1 noti,\
     NaN w l lw 3 dt 1 lc -1 ti 'Autospectra',\
     NaN w l lw 3 dt 2 lc -1 ti 'Cross with matter',\
     for [i=1:5] combi(isim,iz,i,i) u 1:((column(d)/column(d+M))/((column(2*M+c)-column(2*M+s))/(column(2*M+c+L)-column(2*M+s+L)))) w l lw 2 dt 1 lc icol[i] ti word(field_names,i),\
     for [i=1:5] combi(isim,iz,1,i) u 1:((column(d)/column(d+M))/((column(2*M+c)-column(2*M+s))/(column(2*M+c+L)-column(2*M+s+L)))) w l lw 2 dt 2 lc icol[i] noti

     #combi(isim,iz,0,0) u 1:((column(d)/column(d+M))/((column(2*M+c)-column(2*M+s))/(column(2*M+c+L)-column(2*M+s+L)))) w l lw 3 dt 1 lc cols[1] ti word(field_names,1+0),\
     combi(isim,iz,1,1) u 1:((column(d)/column(d+M))/((column(2*M+c)-column(2*M+s))/(column(2*M+c+L)-column(2*M+s+L)))) w l lw 3 dt 1 lc cols[2] ti word(field_names,1+1),\
     combi(isim,iz,2,2) u 1:((column(d)/column(d+M))/((column(2*M+c)-column(2*M+s))/(column(2*M+c+L)-column(2*M+s+L)))) w l lw 3 dt 1 lc cols[3] ti word(field_names,1+2),\
     combi(isim,iz,3,3) u 1:((column(d)/column(d+M))/((column(2*M+c)-column(2*M+s))/(column(2*M+c+L)-column(2*M+s+L)))) w l lw 3 dt 1 lc cols[4] ti word(field_names,1+3),\
     combi(isim,iz,6,6) u 1:((column(d)/column(d+M))/((column(2*M+c)-column(2*M+s))/(column(2*M+c+L)-column(2*M+s+L)))) w l lw 3 dt 1 lc cols[6] ti word(field_names,1+6),\
     combi(isim,iz,0,1) u 1:((column(d)/column(d+M))/((column(2*M+c)-column(2*M+s))/(column(2*M+c+L)-column(2*M+s+L)))) w l lw 3 dt 2 lc cols[2] noti 'CDM',\
     combi(isim,iz,0,2) u 1:((column(d)/column(d+M))/((column(2*M+c)-column(2*M+s))/(column(2*M+c+L)-column(2*M+s+L)))) w l lw 3 dt 2 lc cols[3] noti 'gas',\
     combi(isim,iz,0,3) u 1:((column(d)/column(d+M))/((column(2*M+c)-column(2*M+s))/(column(2*M+c+L)-column(2*M+s+L)))) w l lw 3 dt 2 lc cols[4] noti 'stars',\
     combi(isim,iz,0,6) u 1:((column(d)/column(d+M))/((column(2*M+c)-column(2*M+s))/(column(2*M+c+L)-column(2*M+s+L)))) w l lw 3 dt 2 lc cols[6] noti 'electron_pressure'

unset label

}

}

unset multiplot

}

if(iplot==11 || iplot==12){

# Pressure spectrum variance demonstration

if(print==1){
if(iplot==11) {outfile(name,z)=sprintf('%s_z%1.1f_matter_epressure_variance.eps',name,z)}
if(iplot==12) {outfile(name,z)=sprintf('%s_z%1.1f_epressure_epressure_variance.eps',name,z)}
set output outfile(sim,z)
}

if(iplot==11) {pmin=1e-4; pmax=1e-2}
if(iplot==12) {pmin=1e-8; pmax=1e-5}
set yrange[pmin:pmax]
set ylabel '{/Symbol D}^2(k) / (k/h^{-1} Mpc)^{1.5}'

if(iplot==11) {f1='all';       f2='epressure'; i1=0; i2=6; set title 'matter-electron pressure BAHAMAS spectra'}
if(iplot==12) {f1='epressure'; f2='epressure'; i1=6; i2=6; set title 'electron pressure-electron pressure BAHAMAS spectra'}

set key top left

plot simpk('AGN_TUNED_nu0_v2',mesh,snap,f1,f2) u 1:(($2-$3)/$1**1.5) w p pt 6 lc rgb 'light-blue' noti,\
     simpk('AGN_TUNED_nu0_v3',mesh,snap,f1,f2) u 1:(($2-$3)/$1**1.5) w p pt 6 lc rgb 'light-blue' noti,\
     simpk('AGN_TUNED_nu0_v4',mesh,snap,f1,f2) u 1:(($2-$3)/$1**1.5) w p pt 6 lc rgb 'light-blue' noti,\
     simpk('AGN_TUNED_nu0',mesh,snap,f1,f2)    u 1:(($2-$3)/$1**1.5):($5/$1**1.5) w e pt 6 lc icol[5] ti 'AGN 7p8',\
     simpk('AGN_7p6_nu0',mesh,snap,f1,f2)      u 1:(($2-$3)/$1**1.5) w p pt 6 lc 5 ti 'AGN 7p6',\
     simpk('AGN_8p0_nu0',mesh,snap,f1,f2)      u 1:(($2-$3)/$1**1.5) w p pt 6 lc 7 ti 'AGN 8p0',\
     hmpk(hmpk_name,z,i1,i2)                   u 1:($5/$1**1.5)      w l lw 3 dt 1 lc icol[5] ti 'HMx'

}

if(iplot==13 || iplot==14 || iplot==15 || iplot==16 || iplot==17 || iplot==18 || iplot==19){

unset title

kmin=3e-3
kmax=1e2
set log x
set xlabel 'k / h Mpc^{-1}'
set xrange [kmin:kmax]

dmin=1e-6
dmax=1e3
set log y
set ylabel '{/Symbol D}_{uv}^2(k)'
set format y '10^{%T}'
set yrange [dmin:dmax]

mod='AGN_TUNED_nu0'

labx=0.85
laby=0.2

set key top left

#set label word(z_names,iz) at screen labx,laby

if(iplot==13){

if(print==1){set output 'plots/power_example.eps'}

plot hmpk(mod,z,0,0) u 1:5 w l lc 1 dt 1 lw 3 ti 'Halo model',\
     simpk(mod,mesh,snap,'all','all') u 1:($2-$3):5 w e lc -1 pt 7 ps .5 ti 'Simulation'

}

if(iplot==14 || iplot==15 || iplot==16 || iplot==17 || iplot==18){

if(print==1){
if(iplot==15){set output 'plots/power_example_matter-matter.eps'}
if(iplot==15){set output 'plots/power_example_cdm-cdm.eps'}
if(iplot==16){set output 'plots/power_example_gas-gas.eps'}
if(iplot==17){set output 'plots/power_example_stars-stars.eps'}
if(iplot==18){set output 'plots/power_example_epressure-epressure.eps'}
}

if(iplot==14){c1=1; c2=1; col=1; f=1.} # Matter - Matter
if(iplot==15){c1=2; c2=2; col=2; f=1.} # CDM - CDM
if(iplot==16){c1=3; c2=3; col=3; f=1.} # Gas - Gas
if(iplot==17){c1=4; c2=4; col=4; f=1.} # Stars - Stars
if(iplot==18){c1=5; c2=5; col=6; f=f2} # Electron pressure - Electron pressure


plot hmpk(mod,z,0,0) u 1:5 w l lc 1 dt 1 lw 3 ti ''.word(field_names,1).'',\
     hmpk(mod,z,ifield[c1],ifield[c2]) u 1:(f*$5) w l lc col dt 1 lw 3 ti ''.word(field_names,c1).'',\
     hmpk(mod,z,ifield[c1],ifield[c2]) u 1:(f*$3) w l lc col dt 2 lw 3 noti,\
     hmpk(mod,z,ifield[c1],ifield[c2]) u 1:(f*$4) w l lc col dt 3 lw 3 noti,\
     simpk(mod,mesh,snap,word(fields,1),word(fields,1))   u 1:($2-$3):5 w e lc -1 pt 7 ps .5 ti 'Simulation',\
     simpk(mod,mesh,snap,word(fields,c1),word(fields,c2)) u 1:(f*($2-$3)):(f*$5) w e lc -1 pt 7 ps .5 noti

}

if(iplot==19){

if(print==1){
set output 'plots/power_example_matter-epressure.eps'
show output
}

c1=1; c2=5; col=6; f=fac[1] # Matter - Electron pressure

plot hmpk(mod,z,ifield[c1],ifield[c2]) u 1:(f*$3) w l lc col dt 2 lw 3 noti,\
     hmpk(mod,z,ifield[c1],ifield[c2]) u 1:(f*$4) w l lc col dt 3 lw 3 noti,\
     hmpk(mod,z,0,0) u 1:5 w l lc 1 dt 1 lw 3 ti ''.word(field_names,1).'',\
     hmpk(mod,z,ifield[c1],ifield[c2]) u 1:(f*$5) w l lc col dt 1 lw 3 ti ''.word(field_names,c2).'',\
     hmpk(mod,z,ifield[1], ifield[c2]) u 1:(f*$5) w l lc 1   dt 2 lw 3 noti,\
     simpk(mod,mesh,snap,word(fields,1),word(fields,1))   u 1:($2-$3):5 w e lc -1 pt 7 ps .5 ti 'Simulation',\
     simpk(mod,mesh,snap,word(fields,c1),word(fields,c2)) u 1:(f*($2-$3)):(f*$5) w e lc -1 pt 7 ps .5 noti

}

}

if(iplot==20 || iplot==21 || iplot==22 || iplot==23 || iplot==24 || iplot==25 || iplot==26 || iplot==27 || iplot==28 || iplot==29){

unset title

pow=1.5

icomp=2
print 'icomp automatically set to 2 here'

if(iplot==20 || iplot==25) {f1=1; f2=1} # matter-matter
if(iplot==21 || iplot==26) {f1=2; f2=2} # CDM-CDM
if(iplot==22 || iplot==27) {f1=3; f2=3} # gas-gas
if(iplot==23 || iplot==28) {f1=4; f2=4} # stars-stars
if(iplot==24 || iplot==29) {f1=1; f2=5} # matter-electron pressure

if(print==1){
if(iplot==20 || iplot==25) {set output 'hydro_matter-matter.eps'}
if(iplot==21 || iplot==26) {set output 'hydro_cdm-cdm.eps'}
if(iplot==22 || iplot==27) {set output 'hydro_gas-gas.eps'}
if(iplot==23 || iplot==28) {set output 'hydro_stars-stars.eps'}
if(iplot==24 || iplot==29) {set output 'hydro_matter-epressure.eps'}
}

kmin=0.02
kmax=10.
set xrange [kmin:kmax]

# Delta^2(k)/k^1.5 range
unset log y
set format y
if(iplot==20 || iplot==21 || iplot==22 || iplot==23 || iplot==25) {set ylabel '{/Symbol D}_{uv}^2(k) / [k / h^{-1} Mpc]^{1.5}'}
if(iplot==25 || iplot==26 || iplot==27 || iplot==28 || iplot==29) {set ylabel 'P_{uv}(k) / P_{mm-dmony}(k)'}
set mytics 10

set multiplot layout 2,2

do for [iz=1:4]{

# matter-matter
if(iplot==20 && iz==1) {pmax=25.}
if(iplot==20 && iz==2) {pmax=12.}
if(iplot==20 && iz==3) {pmax=6.}
if(iplot==20 && iz==4) {pmax=2.}

# CDM-CDM
if(iplot==21 && iz==1) {pmax=20.}
if(iplot==21 && iz==2) {pmax=10.}
if(iplot==21 && iz==3) {pmax=5.}
if(iplot==21 && iz==4) {pmax=1.5}

# gas-gas
if(iplot==22 && iz==1) {pmax=0.4}
if(iplot==22 && iz==2) {pmax=0.2}
if(iplot==22 && iz==3) {pmax=0.14}
if(iplot==22 && iz==4) {pmax=0.06}

# stars-stars
if(iplot==23 && iz==1) {pmax=0.03}
if(iplot==23 && iz==2) {pmax=0.014}
if(iplot==23 && iz==3) {pmax=0.007}
if(iplot==23 && iz==4) {pmax=0.0016}

# matter-epressure
if(iplot==24 && iz==1) {pmax=0.006}
if(iplot==24 && iz==2) {pmax=0.0018}
if(iplot==24 && iz==3) {pmax=0.0006}
if(iplot==24 && iz==4) {pmax=8e-5}

if(iplot==25){pmin=0.75; pmax=1.02}
if(iplot==26){pmin=0.65; pmax=0.75}
if(iplot==27){pmin=0.00; pmax=0.030}
if(iplot==28){pmin=0.00; pmax=0.002}
if(iplot==29){pmin=0.00; pmax=3e-4}

if(iz==1) {zlab='z = 0.0'}
if(iz==2) {zlab='z = 0.5'}
if(iz==3) {zlab='z = 1.0'}
if(iz==4) {zlab='z = 2.0'}

if(iplot==20 || iplot==21 || iplot==22 || iplot==23 || iplot==24){
unset key
if(iz==1) {set key top left}
set yrange [0.:pmax]
plot NaN w p pt 7 ps .5 lc -1 ti 'Simulations',\
     NaN w l dt 1 lw  2 lc -1 ti 'Halo model',\
     for [i=2:words(sims)] simpk(word(sims,i),mesh,word(snaps,iz),word(fields,f1),word(fields,f2)) u 1:((column(c)-column(s))/(column(1)**pow)) w p pt 7 ps .5 lc rgb word(sim_cols,i) noti,\
     for [i=2:words(hmpk_names)] hmpk(word(hmpk_names,i),zs[iz],ifield[f1],ifield[f2]) u 1:((column(d)/column(1)**pow)) w l lw 2 dt 1 lc rgb word(hmpk_cols,i) ti word(hmpk_names,i)
}


if(iplot==25 || iplot==26 || iplot==27 || iplot==28 || iplot==29){
unset key
if(iz==1) {set key bottom left}
set yrange [pmin:pmax]
set label zlab at graph 0.1,0.9
plot NaN w p pt 7 ps .5 lc -1 ti 'Simulations',\
     NaN w l dt 1 lw  2 lc -1 ti 'Halo model',\
     for [i=2:words(sims)] '<paste '.simpk(word(sims,i),mesh,word(snaps,iz),word(fields,f1),word(fields,f2)).' '.simpk(word(sims,1),mesh,word(snaps,iz),word(fields,1),word(fields,1)).'' u 1:(column(c)-column(s))/(column(c+5)-column(s+5)) w p pt 7 ps .5 lc rgb word(sim_cols,i) noti,\
     for [i=2:words(hmpk_names)] '<paste '.hmpk(word(hmpk_names,i),zs[iz],ifield[f1],ifield[f2]).' '.hmpk(word(hmpk_names,1),zs[iz],ifield[1],ifield[1]).'' u 1:(column(d)/column(d+5)) w l lw 2 dt 1 lc rgb word(hmpk_cols,i) ti word(hmpk_names,i)
unset label
}

}

unset multiplot

}

if(iplot==30 || iplot==31){

if(print==1 && iplot==30) {set output 'power_AGN_model.eps'}
if(print==1 && iplot==31) {set output 'residuals_AGN_model.eps'}

unset title

set multiplot layout 5,4 margins 0.1,0.98,0.1,0.98 spacing 0.01,0.01

do for [ifi=1:5]{

if(ifi==1){f1=1; f2=1; rmin=0.75; rmax=1.05;   pmin=5e0;  pmax=3e1;    plab='k^{3/2}P_{mm}^2(k)'; rlab='P_{mm}(k) / P_{mm-dmony}(k)'}
if(ifi==2){f1=2; f2=2; rmin=0.65; rmax=0.80;   pmin=4e0;  pmax=2e1;    plab='k^{3/2}P_{cc}^2(k)'; rlab='P_{cc}(k) / P_{mm-dmony}(k)'}
if(ifi==3){f1=3; f2=3; rmin=0.00; rmax=0.03;   pmin=1.1e-2; pmax=1e0;  plab='k^{3/2}P_{gg}^2(k)'; rlab='P_{gg}(k) / P_{mm-dmony}(k)'}
if(ifi==4){f1=4; f2=4; rmin=0.00; rmax=0.0025; pmin=1.1e-4; pmax=1e-1; plab='k^{3/2}P_{**}^2(k)'; rlab='P_{**}(k) / P_{mm-dmony}(k)'}
if(ifi==5){f1=1; f2=5; rmin=0.00; rmax=2.5e-4; pmin=1.1e-4; pmax=1e-2; plab='k^{3/2}P_{mP}^2(k)'; rlab='P_{mP}(k) / P_{mm-dmony}(k)'}

set xlabel ''
set format x ''
if(ifi==5){set xlabel klab; set format x}

do for [iz=1:4]{

# Growth factors
g2=gs[iz]**2

if(iplot==30) {set yrange [pmin:pmax]; set log y}
if(iplot==31) {set yrange [rmin:rmax]; unset log y}

set format y ''
set ylabel ''
if(iplot==30) {ylab=plab}
if(iplot==31) {ylab=rlab}
if(iplot==30 && iz==1){set format y '10^{%T}'; set ylabel ylab}
if(iplot==31 && iz==1){set format y; set ylabel ylab}

zlab=''
if(ifi==1 && iz==1){zlab='z = 0.0'}
if(ifi==1 && iz==2){zlab='z = 0.5'}
if(ifi==1 && iz==3){zlab='z = 1.0'}
if(ifi==1 && iz==4){zlab='z = 2.0'}

simti=''
haloti=''
if(iz==2 && ifi==1) {simti='Simulations'; haloti='Halomodel'}

unset key
if(iplot==30 && iz==1 && ifi==1) {set key bottom right}
if(iplot==31 && iz==1 && ifi==1) {set key bottom left}

set label zlab at graph 0.1,0.9

if(iplot==30){
unset key
plot NaN w p pt 7 ps .5 lc -1 ti simti,\
     NaN w l dt 1 lw  2 lc -1 ti haloti,\
     for [i=2:words(sims_few)] simpk(word(sims,i),mesh,word(snaps,iz),word(fields,f1),word(fields,f2)) u 1:((column(c)-column(s))/(g2*column(1)**pow)) w p pt 7 ps .25 lc rgb word(sim_cols,i) noti,\
     for [i=2:words(hmpk_names)] hmpk(word(hmpk_names,i),zs[iz],ifield[f1],ifield[f2]) u 1:((column(d)/(g2*column(1)**pow))) w l lw 2 dt 1 lc rgb word(hmpk_cols,i) ti word(hmpk_names,i)
}

if(iplot==31){
plot NaN w p pt 7 ps .5 lc -1 ti simti,\
     NaN w l dt 1 lw  2 lc -1 ti haloti,\
     for [i=2:words(sims_few)] '<paste '.simpk(word(sims,i),mesh,word(snaps,iz),word(fields,f1),word(fields,f2)).' '.simpk(word(sims,1),mesh,word(snaps,iz),word(fields,1),word(fields,1)).'' u 1:(column(c)-column(s))/(column(c+5)-column(s+5)) w p pt 7 ps .5 lc rgb word(sim_cols,i) noti,\
     for [i=2:words(hmpk_names)] '<paste '.hmpk(word(hmpk_names,i),zs[iz],ifield[f1],ifield[f2]).' '.hmpk(word(hmpk_names,1),zs[iz],ifield[1],ifield[1]).'' u 1:(column(d)/column(d+5)) w l lw 2 dt 1 lc rgb word(hmpk_cols,i) ti word(hmpk_names,i)
}

unset label

}}

unset multiplot

}

if(iplot==32){

if(print==1) {set term post enh col size 10,12 font ',8'; set output 'paper/big_AGN_model.eps'}

unset title
unset mytics

set multiplot layout 10,4 margins 0.06,0.98,0.04,0.98 spacing 0.01,0.01

do for [ifi=1:5]{

if(ifi==1){f1=1; f2=1; rmin=0.75; rmax=1.05;   pmin=0; pmax=30;    plab='k^{3/2} P_{mm}^2(k)'; rlab='P_{mm}(k) / P_{mm-dmony}(k)'}
if(ifi==2){f1=2; f2=2; rmin=0.65; rmax=0.80;   pmin=0; pmax=20;    plab='k^{3/2} P_{cc}^2(k)'; rlab='P_{cc}(k) / P_{mm-dmony}(k)'}
if(ifi==3){f1=3; f2=3; rmin=0.00; rmax=0.03;   pmin=0; pmax=0.4;   plab='k^{3/2} P_{gg}^2(k)'; rlab='P_{gg}(k) / P_{mm-dmony}(k)'}
if(ifi==4){f1=4; f2=4; rmin=0.00; rmax=0.0025; pmin=0; pmax=0.04;  plab='k^{3/2} P_{**}^2(k)'; rlab='P_{**}(k) / P_{mm-dmony}(k)'}
if(ifi==5){f1=1; f2=5; rmin=0.00; rmax=2.5e-4; pmin=0; pmax=0.006; plab='k^{3/2} P_{mP}^2(k)'; rlab='P_{mP}(k) / P_{mm-dmony}(k)'}

do for [ir=1:2]{

set xlabel ''
set format x ''
if(ir==2 && ifi==5){set xlabel klab; set format x}

do for [iz=1:4]{

# Growth factors
g2=gs[iz]**2

if(ir==1) {set yrange [pmin:pmax]; unset log y}
if(ir==2) {set yrange [rmin:rmax]; unset log y}

set format y ''
set ylabel ''
if(ir==1) {ylab=plab}
if(ir==2) {ylab=rlab}
if(ir==1 && iz==1){set format y; set ylabel ylab}
if(ir==2 && iz==1){set format y; set ylabel ylab}

zlab=''
if(ir==1 && ifi==1 && iz==1){zlab='z = 0.0'}
if(ir==1 && ifi==1 && iz==2){zlab='z = 0.5'}
if(ir==1 && ifi==1 && iz==3){zlab='z = 1.0'}
if(ir==1 && ifi==1 && iz==4){zlab='z = 2.0'}

simti=''
haloti=''
if(iz==2 && ifi==1) {simti='Simulations'; haloti='Halomodel'}

unset key
if(ir==1 && iz==1 && ifi==1) {set key bottom right}
if(ir==2 && iz==1 && ifi==1) {set key bottom left}

set label zlab at graph 0.05,0.86

if(ir==1){
unset key
plot NaN w p pt 7 ps .5 lc -1 ti simti,\
     NaN w l dt 1 lw  2 lc -1 ti haloti,\
     for [i=2:words(sims_few)] simpk(word(sims,i),mesh,word(snaps,iz),word(fields,f1),word(fields,f2)) u 1:((column(c)-column(s))/(g2*column(1)**pow)):((column(5))/(g2*column(1)**pow)) w e pt 7 ps .25 lc rgb word(sim_cols,i) noti,\
     for [i=2:words(hmpk_names)] hmpk(word(hmpk_names,i),zs[iz],ifield[f1],ifield[f2]) u 1:((column(d)/(g2*column(1)**pow))) w l lw 2 dt 1 lc rgb word(hmpk_cols,i) ti word(hmpk_names,i)
}

if(ir==2){
plot NaN w p pt 7 ps .5 lc -1 ti simti,\
     NaN w l dt 1 lw  2 lc -1 ti haloti,\
     for [i=2:words(sims_few)] '<paste '.simpk(word(sims,i),mesh,word(snaps,iz),word(fields,f1),word(fields,f2)).' '.simpk(word(sims,1),mesh,word(snaps,iz),word(fields,1),word(fields,1)).'' u 1:(column(c)-column(s))/(column(c+5)-column(s+5)) w p pt 7 ps .25 lc rgb word(sim_cols,i) noti,\
     for [i=2:words(hmpk_names)] '<paste '.hmpk(word(hmpk_names,i),zs[iz],ifield[f1],ifield[f2]).' '.hmpk(word(hmpk_names,1),zs[iz],ifield[1],ifield[1]).'' u 1:(column(d)/column(d+5)) w l lw 2 dt 1 lc rgb word(hmpk_cols,i) ti word(hmpk_names,i)
}

unset label

}}}

unset multiplot

}

show output
