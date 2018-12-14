unset multiplot
reset

cmmi='/Users/Mead/Fonts/cmmi10.pfb'

if(!exists('print')){print=0}
if(print==0){set term aqua dashed font ',10'; ell='l'}
if(print==1){set term post enh col dashed fontfile cmmi font ',10'; ell='{/cmmi10 \140}'}

# File name and location functions
hm_file(mod,x,y)=sprintf('data/%s_%s-%s.dat',mod,x,y)
mod='triad_Cl_AGN'

print ''
if(!exists('iplot')){iplot=1}
print 'iplot = 1: Triad 3: Cross spectra; separate plots'
print 'iplot = 2: Triad 2: Cross spectra; same plot'
print 'iplot = 3: Triad 3: Cross spectra and residuals'
print 'iplot = 4: Triad 2: Cross spectra; separate plots'
print 'iplot = 5: Triad 3: Residuals'
print 'iplot = 6: Triad 3: Directly KiDS(0.1->0.9)-y'
if(!exists('iplot')){iplot=2}
print 'iplot = ', iplot
print ''

if(iplot==1){set output 'triad.eps'}
if(iplot==2){set output 'triad.eps'}
if(iplot==3){set output 'paper/triad.eps'}
if(iplot==4){set output 'triad.eps'}
if(iplot==5){set output 'triad_residuals.eps'}
if(iplot==6){set output 'triad_direct.eps'}

# Set units for y axis
print 'icl = 1: C(l)'
print 'icl = 2: l(l+1)C(l)/2pi'
print 'icl = 3: (l+1)C(l)'
if(!exists('icl')){icl=2}
print 'icl = ', icl
print ''

# x axis
set log x
set xlabel ''.ell.''

# y axis
if(icl==1) {DC(l,Cl)=Cl; ylab='C_{uv}('.ell.')'; set key top right}
if(icl==2) {DC(l,Cl)=l*(l+1)*Cl/(2.*pi); ylab=''.ell.'('.ell.'+1)C_{uv}('.ell.') / 2{/Symbol p}'; set key top left}
if(icl==3) {DC(l,Cl)=(l+1)*Cl; ylab='('.ell.'+1)C_{uv}('.ell.')'; set key top left}
set log y
set format y '10^{%T}'
set ylabel ylab
#set yrange [clmin:clmax]

# Halo model Cl file bases
mods="'triad_Cl_AGN-lo' 'triad_Cl_AGN' 'triad_Cl_AGN-hi'"

# BAHAMAS simulation names
sims="'LOW' 'TUNED' 'HIGH'"
sim_names="'AGN-lo' 'AGN' 'AGN-hi'"

# Colour scheme for feedback models
cols="'dark-yellow' 'blue' 'dark-plum'"

# Different tracer types and names
hms="'CMB' 'y' 'gal_z0.1-0.9' 'gal_z0.1-0.5' 'gal_z0.5-0.9'"
fields="'CMBkappa' 'tSZ' 'shear_z0.1-0.9' 'shear_z0.1-0.5' 'shear_z0.5-0.9'"
field_names="'{/Symbol f}' 'y' '{/Symbol g} (z = 0.1 -> 0.9)' '{/Symbol g} (z = 0.1 -> 0.5)' '{/Symbol g} (z = 0.5 -> 0.9)'"

# Planck beam
fwhm = 10. # FWHM [arcmin]
fwhm = fwhm/60. # FWHM [deg]
fwhm = fwhm*pi/180. # FWHM [rad]
sigma = fwhm/(2.*sqrt(2.*log(2.))) # sigma [rad]
beam(l)=exp(-0.5*l*(l+1.)*sigma**2)

# Distance and function to shift simulation points for clarity
# Keep the central simulation points unchanged (i=2)
disp(i,n)=1.+0.03*real(i-2)

# l axis
lmin=90.
lmax=3000.
set log x
set xlabel ''.ell.''
set xrange [lmin:lmax]

# y axis
set log y
set format y '10^{%T}'
set ylabel ylab
set yrange [*:*]

nsim=words(sims)
print 'Number of simulations: ', nsim
print ''

## TRIAD 3 ##

if(iplot==1 || iplot==3 || iplot==5 || iplot==6){

# Location of C(ell) measured from the simulations
sim_file(x,y,mod)=sprintf('/Users/Mead/Physics/people/Tilman/BAHAMAS_triad_3/mean_Cl_%s-%s_%s.txt', x, y, mod)

if(iplot==1 || iplot==5){

set multiplot layout 3,3

do for [i=1:9]{

if(i==1) {f1=3; f2=3; clmin=2e-6;  clmax=1e-4; sim_names="'AGN-lo' 'AGN' 'AGN-hi'"; t1=''; t2=''}
if(i==2) {f1=4; f2=4; clmin=2e-6;  clmax=1e-4; sim_names=""; t1='Halo model'; t2='BAHAMAS data'}
if(i==3) {f1=5; f2=5; clmin=2e-6;  clmax=1e-4; t1=''; t2=''}
if(i==4) {f1=4; f2=3; clmin=2e-6;  clmax=1e-4; t1=''; t2=''}
if(i==5) {f1=5; f2=3; clmin=2e-6;  clmax=1e-4; t1=''; t2=''}
if(i==6) {f1=5; f2=4; clmin=2e-6;  clmax=1e-4; t1=''; t2=''}
if(i==7) {f1=2; f2=3; clmin=3e-10; clmax=8e-9; t1=''; t2=''}
if(i==8) {f1=2; f2=4; clmin=3e-10; clmax=8e-9; t1=''; t2=''}
if(i==9) {f1=2; f2=5; clmin=3e-10; clmax=8e-9; t1=''; t2=''}

if(iplot==1){

set yrange [*:*]
if(icl==2) {set yrange [clmin:clmax]}

set title ''.word(field_names,f1).'-'.word(field_names,f2).'' enh
plot NaN w l lw 2 dt 1  lc -1 ti t1,\
     NaN w p pt 7 ps .5 lc -1 ti t2,\
     for [j=1:nsim] hm_file(word(mods,j),word(hms,f1),word(hms,f2)) u 1:(DC($1,$2)) w l lw 2 lc rgb word(cols,j) dt 1 noti,\
     for [j=1:nsim] sim_file(word(fields,f1),word(fields,f2),word(sims,j)) u ($1*disp(j,nsim)):(DC($1,$2)):(DC($1,$3)) w errorbars lc rgb word(cols,j) pt 7 ps .5 ti word(sim_names,j)

}

if(iplot==5){

set yrange [0.5:1.5]
set format y
unset log y
set ylabel 'C_{mod}('.ell.') / C_{sim}('.ell.')'

set title ''.word(field_names,f1).'-'.word(field_names,f2).'' enh
plot 1 w l lt -1 noti,\
     for [j=1:nsim] '<paste '.hm_file(word(mods,j),word(hms,f1),word(hms,f2)).' '.sim_file(word(fields,f1),word(fields,f2),word(sims,j)).'' u ($1*disp(j,nsim)):(column(2)/column(2+3)):(column(6)/column(5)) w e ps .5 pt 7 lc rgb word(cols,j) dt 1 noti

}

}

unset multiplot

}

if(iplot==3){

DC(l,Cl)=l*(l+1)*Cl

if(print==1) {set term post enh col font ',8'}

set multiplot layout 6,3 margins 0.06,0.98,0.06,0.98 spacing 0.01,0.01

set key bottom right

# Loop over rows of plots (gamma-gamma auto; gamma-gamma cross; gamma-y))
do for [i=1:3]{

# Loop over actual C(l) and residuals
do for [j=1:2]{

# Loop over columns
do for [k=1:3]{

if(i==1 && k==1) {f1=3; f2=3; clmin=0; clmax=40; sim_names="'AGN-lo' 'AGN' 'AGN-hi'"; t1=''; t2=''}
if(i==1 && k==2) {f1=4; f2=4; clmin=0; clmax=40; sim_names=""; t1='Halo model'; t2='BAHAMAS data'}
if(i==1 && k==3) {f1=5; f2=5; clmin=0; clmax=40; t1=''; t2=''}
if(i==2 && k==1) {f1=3; f2=4; clmin=0; clmax=40; t1=''; t2=''}
if(i==2 && k==2) {f1=4; f2=5; clmin=0; clmax=40; t1=''; t2=''}
if(i==2 && k==3) {f1=5; f2=3; clmin=0; clmax=40; t1=''; t2=''}
if(i==3 && k==1) {f1=2; f2=3; clmin=0; clmax=40; t1=''; t2=''}
if(i==3 && k==2) {f1=2; f2=4; clmin=0; clmax=40; t1=''; t2=''}
if(i==3 && k==3) {f1=2; f2=5; clmin=0; clmax=40; t1=''; t2=''}

set xlabel ''; set format x ''
if(i==3 && j==2) {set xlabel ''.ell.''; set format x}

set ylabel ''
set format y ''

# C(ell)
if(j==1){
if(i==1) {ylab='10^{5} '.ell.'('.ell.'+1)C_{{/Symbol g}{/Symbol g}}('.ell.')'; b=1e5}
if(i==2) {ylab='10^{5} '.ell.'('.ell.'+1)C_{{/Symbol g}{/Symbol g}}('.ell.')'; b=1e5}
if(i==3) {ylab='10^{9} '.ell.'('.ell.'+1)C_{y{/Symbol g}}('.ell.')'; b=1e9}
if(k==1) {set ylabel ylab; set format y}
set yrange [clmin:clmax]; unset log y
}

# Residuals
if(j==2){
if(i==1) {ylab='C_{{/Symbol g}{/Symbol g}}('.ell.') / C_{sim}('.ell.')'}
if(i==2) {ylab='C_{{/Symbol g}{/Symbol g}}('.ell.') / C_{sim}('.ell.')'}
if(i==3) {ylab='C_{y{/Symbol g}}('.ell.') / C_{sim}('.ell.')'}
if(k==1) {set ylabel ylab; set format y}
set yrange [0.5:1.5]; unset log y
}

set key top right

if(j==1){

set label ''.word(field_names,f1).'-'.word(field_names,f2).'' enh at graph 0.05,0.9
plot NaN w l lw 2 dt 1  lc -1 ti t1,\
     NaN w p pt 7 ps .5 lc -1 ti t2,\
     for [j=1:nsim] hm_file(word(mods,j),word(hms,f1),word(hms,f2)) u 1:(b*DC($1,$2)) w l lw 2 lc rgb word(cols,j) dt 1 noti,\
     for [j=1:nsim] sim_file(word(fields,f1),word(fields,f2),word(sims,j)) u ($1*disp(j,nsim)):(b*DC($1,$2)):(b*DC($1,$3)) w errorbars lc rgb word(cols,j) pt 7 ps .5 ti word(sim_names,j)
unset label
}

if(j==2){

plot 1 w l lt -1 noti,\
     for [j=1:nsim] '<paste '.hm_file(word(mods,j),word(hms,f1),word(hms,f2)).' '.sim_file(word(fields,f1),word(fields,f2),word(sims,j)).'' u ($1*disp(j,nsim)):(column(2)/column(2+3)):(column(6)/column(5)) w e ps .5 pt 7 lc rgb word(cols,j) dt 1 noti

}

}}}

unset multiplot

}

if(iplot==6){

f1=2
f2=3
print 'Field 1: ', word(field_names,f1)
print 'Field 2: ', word(field_names,f2)
print ''

directs="'data/triad_Cl_direct_AGN-lo_y-gal_z0.1-0.9.dat' 'data/triad_Cl_direct_AGN_y-gal_z0.1-0.9.dat' 'data/triad_Cl_direct_AGN-hi_y-gal_z0.1-0.9.dat'"

set log x

#clmin=2e-10
#clmax=5e-9
clmin=0.
clmax=4.5e-9
unset log y
set format y

set yrange [*:*]
if(icl==2) {set yrange [clmin:clmax]}

plot NaN w l dt 1 lw  2 lc -1 ti 'Halo model',\
     NaN w l dt 2 lw  2 lc -1 ti 'Direct integration of 3D spectra',\
     NaN w p pt 7 ps .5 lc -1 ti 'Direction simulation measurements',\
     for [i=1:nsim] hm_file(word(mods,i),word(hms,f1),word(hms,f2)) u 1:(DC($1,$2)) w l lw 2 lc rgb word(cols,i) dt 1 noti,\
     for [i=1:nsim] sim_file(word(fields,f1),word(fields,f2),word(sims,i)) u ($1*disp(i,nsim)):(DC($1,$2)):(DC($1,$3)) w errorbars lc rgb word(cols,i) pt 7 ps .5 ti word(sim_names,i),\
     for [i=1:nsim] word(directs,i) u 1:(DC($1,$2)) w l dt 2 lc rgb word(cols,i) lw 2 noti

}

}

## ##

## TRIAD 1 and 2 ##

if(iplot==2 || iplot==4 || iplot==6){

# Location of C(ell) measured from the simulations
sim_file(x,y,mod)=sprintf('/Users/Mead/Physics/people/Tilman/BAHAMAS_triad_1_and_2/mean_Cl_%s-%s_%s.txt', x, y, mod)

if(iplot==2){

set multiplot layout 1,2

set yrange [*:*]
if(icl==2) {clmin=1e-5; clmax=2e-4; set yrange [1e-5:2e-4]}

array ifield1[3]
array ifield2[3]

ifield1[1]=3 # gal_z0.1-0.9
ifield2[1]=1 # CMB

ifield1[2]=4 # gal_z0.1-0.5
ifield2[2]=1 # CMB

ifield1[3]=5 # gal_z0.5-0.5
ifield2[3]=1 # CMB

# Make the plot
plot for [i=1:3] hm_file(mod,word(hms,ifield1[i]),word(hms,ifield2[i])) u 1:(DC($1,$2)) w l lw 3 lc -1 dt i ti ''.word(field_names,ifield1[i]).'-'.word(field_names,ifield2[i]).'',\
     for [i=1:3] for [j=1:nsim] sim_file(word(fields,ifield1[i]),word(fields,ifield2[i]),word(sims,j)) u ($1*disp(j,nsim)):(DC($1,$2)):(DC($1,$3)) w errorbars lc rgb word(cols,j) pt 7 noti word(sim_names,j),\
     for [j=1:3] NaN w p lc rgb word(cols,j) pt 7 ti word(sim_names,j)

#for [i=1:3] hm_file(mod,word(hms,ifield1[i]),word(hms,ifield2[i])) u 1:(column(c)) w l lw 3 lc -1 dt i ti ''.word(field_names,ifield1[i]).'-'.word(field_names,ifield2[i]).'',\
#for [i=1:3] for [j=1:nsim] sim_file(word(fields,ifield1[i]),word(fields,ifield2[i]),word(sims,j)) u ($1*disp(j,nsim)):(($1**p1)*(($1+1)**p2)*$2/(2.*pi)**p3):(($1**p1)*(($1+1)**p2)*$3/(2.*pi)**p3) w errorbars lc rgb word(cols,j) pt 7 noti word(sim_names,j),\
#for [j=1:3] NaN w p lc rgb word(cols,j) pt 7 ti word(sim_names,j)

set yrange [*:*]
if(icl==2) {clmin=2e-10; clmax=2e-8; set yrange [clmin:clmax]}

array ifield1[4]
array ifield2[4]

ifield1[1]=1 # CMB
ifield2[1]=2 # y

ifield1[2]=2 # y
ifield2[2]=3 # gal_z0.1-0.9

ifield1[3]=2 # y
ifield2[3]=4 # gal_z0.1-0.5

ifield1[4]=2 # y
ifield2[4]=5 # gal_z0.5-0.5

# Make the plot
plot for [i=1:4] hm_file(mod,word(hms,ifield1[i]),word(hms,ifield2[i])) u 1:(beam($1)*DC($1,$2)) w l lw 3 lc -1 dt i ti ''.word(field_names,ifield1[i]).'-'.word(field_names,ifield2[i]).'',\
     for [i=1:4] for [j=1:nsim] sim_file(word(fields,ifield1[i]),word(fields,ifield2[i]),word(sims,j)) u ($1*disp(j,nsim)):(DC($1,$2)):(DC($1,$3)) w errorbars lc rgb word(cols,j) pt 7 noti word(sim_names,j)

#plot for [i=1:4] hm_file(mod,word(hms,ifield1[i]),word(hms,ifield2[i])) u 1:(beam($1)*column(c)) w l lw 3 lc -1 dt i ti ''.word(field_names,ifield1[i]).'-'.word(field_names,ifield2[i]).'',\
#for [i=1:4] for [j=1:nsim] sim_file(word(fields,ifield1[i]),word(fields,ifield2[i]),word(sims,j)) u ($1*disp(j,nsim)):(($1**p1)*(($1+1)**p2)*$2/(2.*pi)**p3):(($1**p1)*(($1+1)**p2)*$3/(2.*pi)**p3) w errorbars lc rgb word(cols,j) pt 7 noti word(sim_names,j)

unset multiplot

}

if(iplot==4){

set multiplot layout 3,3

do for [i=1:7]{

if(i==1) {f1=3; f2=1; clmin=1e-5;  clmax=2e-4; b(x)=1.; sim_names="'AGN-lo' 'AGN' 'AGN-hi'"; t1=''; t2=''}
if(i==2) {f1=4; f2=1; clmin=1e-5;  clmax=2e-4; b(x)=1.; sim_names=""; t1='Halo model'; t2='BAHAMAS data'}
if(i==3) {f1=5; f2=1; clmin=1e-5;  clmax=2e-4; b(x)=1.; t1=''; t2=''}
if(i==4) {f1=2; f2=3; clmin=3e-10; clmax=3e-9; b(x)=beam(x)}
if(i==5) {f1=2; f2=4; clmin=3e-10; clmax=3e-9; b(x)=beam(x)}
if(i==6) {f1=2; f2=5; clmin=3e-10; clmax=3e-9; b(x)=beam(x)}
if(i==7) {f1=1; f2=2; clmin=1e-9;  clmax=1e-8; b(x)=beam(x)}

set yrange [*:*]
if(icl==2) {set yrange [clmin:clmax]}

set title ''.word(field_names,f1).'-'.word(field_names,f2).'' enh
plot NaN w l lw 2 dt 1  lc -1 ti t1,\
     NaN w p pt 7 ps .5 lc -1 ti t2,\
     for [i=1:nsim] hm_file(word(mods,i),word(hms,f1),word(hms,f2)) u 1:(b($1)*DC($1,$2)) w l lw 2 lc rgb word(cols,i) dt 1 noti,\
     for [i=1:nsim] sim_file(word(fields,f1),word(fields,f2),word(sims,i)) u ($1*disp(i,nsim)):(DC($1,$2)):(DC($1,$3)) w errorbars lc rgb word(cols,i) pt 7 ps .5 ti word(sim_names,i)

#for [i=1:nsim] hm_file(word(mods,i),word(hms,f1),word(hms,f2)) u 1:(column(c)*b($1)) w l lw 2 lc rgb word(cols,i) dt 1 noti,\
#for [i=1:nsim] sim_file(word(fields,f1),word(fields,f2),word(sims,i)) u ($1*disp(i,nsim)):(($1**p1)*(($1+1)**p2)*$2/(2.*pi)**p3):(($1**p1)*(($1+1)**p2)*$3/(2.*pi)**p3) w errorbars lc rgb word(cols,i) pt 7 ps .5 ti word(sim_names,i)

}

unset multiplot

}

}

## ##

show output
