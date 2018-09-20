reset

cmmi='/Users/Mead/Fonts/cmmi10.pfb'

if(!exists('print')){print=0}
if(print==0){set term aqua dashed; ell='l'}
if(print==1){set term post enh col dashed fontfile cmmi font ',16'; set output 'paper/triad.eps'; ell='{/cmmi10 \140}'}

# File name and location functions
hm_file(x,y)=sprintf('data/triad_Cl_%s-%s.dat',x,y)

print ''
print 'icomp = 1: cosmo-OWLS'
print 'icomp = 2: BAHAMAS'
if(!exists('icomp')){icomp=2}
print 'icomp = ', icomp
print ''

# 1 - C(l)
# 2 - l(l+1)C(l)/2pi
print 'icl = 1: C(l)'
print 'icl = 2: l(l+1)C(l)/2pi'
if(!exists('icl')){icl=2}
print 'icl = ', icl
print ''

if(icomp==1){

# cosmo-OWLS feedback simulations
#sims="'REF' 'NOCOOL' 'AGN_8.0' 'AGN_8.5' 'AGN_8.7'"
#sim_names="'REF' 'NOCOOL' 'AGN 8.0' 'AGN 8.5' 'AGN 8.7'"
sims="'AGN_8.0' 'AGN_8.5' 'AGN_8.7'"
sim_names="'AGN 8.0' 'AGN 8.5' 'AGN 8.7'"

cols="'gold' 'orange' 'red'"

#Types of fields and names
fields="'CMB' 'y' 'gal'"
field_names="'{/Symbol f}' 'y' '{/Symbol g}'"

sim_file(x,y,mod)=sprintf('/Users/Mead/Physics/people/Tilman/cosmo-OWLS/Cl/Cl_%s-%s_%s.txt', x, y, mod)

}

if(icomp==2){

# Neutrino simulations (BAHAMAS)
#sims="'0.00_eV' '0.06_eV' '0.12_eV' '0.24_eV' '0.48_eV'"
#names="'0.00 eV' '0.06 eV' '0.12 eV' '0.24 eV' '0.48 eV'"

sims="'LOW' 'TUNED' 'HIGH'"
sim_names="'AGN-lo' 'AGN' 'AGN-hi'"

cols="'gold' 'orange' 'red'"

hms="'CMB' 'y' 'gal_z0.1-0.9' 'gal_z0.1-0.5' 'gal_z0.5-0.9'"
fields="'CMBkappa' 'tSZ' 'shear_z0.1-0.9' 'shear_z0.1-0.5' 'shear_z0.5-0.9'"
field_names="'{/Symbol f}' 'y' '{/Symbol g} (z = 0.1 -> 0.9)' '{/Symbol g} (z = 0.1 -> 0.5)' '{/Symbol g} (z = 0.5 -> 0.9)'"

# Location of C(ell) measured from the simulations
sim_file(x,y,mod)=sprintf('/Users/Mead/Physics/people/Tilman/BAHAMAS/mean_Cl_%s-%s_%s.txt', x, y, mod)

# Beam
fwhm = 10. # FWHM [arcmin]
fwhm = fwhm/60. # FWHM [deg]
fwhm = fwhm*pi/180. # FWHM [rad]
sigma = fwhm/(2.*sqrt(2.*log(2.))) # sigma [rad]
beam(l)=exp(-0.5*l*(l+1.)*sigma**2.)

}

# x axis
ellmin=90.
ellmax=3500.
set log x
set xlabel ''.ell.''
set xrange [ellmin:ellmax]

# y axis
set log y
set format y '10^{%T}'
set yrange [*:*]

# Distance and function to shift simulation points for clarity
disp(i,n)=1.+0.03*real(i-1)/real(n-1)

# Factor to shift down CMB-gal curves
f=5e-4
print 'Note that CMB-galaxy power has been multiplied by a factor to bring it down:'
print 'Factor: ', f
print ''

if(icl==1){p1=0; p2=0; p3=0; ylab='C_{ij}('.ell.')'                                 ; c=2; clmin=1e-14; clmax=2e-12; set key top right}
if(icl==2){p1=1; p2=1; p3=1; ylab=''.ell.'('.ell.'+1)C_{ij}('.ell.') / 2{/Symbol p}'; c=3; clmin=1e-10; clmax=3e-7;  set key top left}

set ylabel ylab
set yrange [clmin:clmax]

nsim=words(sims)

# Make the plot
plot hm_file(word(hms,3),word(hms,1)) u 1:(f*column(c)) w l lw 3 lc -1 dt 1 ti ''.word(field_names,3).'-'.word(field_names,1).'',\
     hm_file(word(hms,4),word(hms,1)) u 1:(f*column(c)) w l lw 3 lc -1 dt 4 ti ''.word(field_names,4).'-'.word(field_names,1).'',\
     hm_file(word(hms,5),word(hms,1)) u 1:(f*column(c)) w l lw 3 lc -1 dt 5 ti ''.word(field_names,5).'-'.word(field_names,1).'',\
     hm_file(word(hms,1),word(hms,2)) u 1:(column(c)*beam($1)) w l lw 3 lc -1 dt 2 ti ''.word(field_names,1).'-'.word(field_names,2).'',\
     hm_file(word(hms,2),word(hms,3)) u 1:(column(c)*beam($1)) w l lw 3 lc -1 dt 3 ti ''.word(field_names,2).'-'.word(field_names,3).'',\
     hm_file(word(hms,2),word(hms,4)) u 1:(column(c)*beam($1)) w l lw 3 lc -1 dt 6 ti ''.word(field_names,2).'-'.word(field_names,4).'',\
     hm_file(word(hms,2),word(hms,5)) u 1:(column(c)*beam($1)) w l lw 3 lc -1 dt 7 ti ''.word(field_names,2).'-'.word(field_names,5).'',\
     for [i=1:nsim] sim_file(word(fields,3),word(fields,1),word(sims,i)) u ($1*disp(i,nsim)):(f*($1**p1)*(($1+1)**p2)*$2/(2.*pi)**p3):(f*($1**p1)*(($1+1)**p2)*$3/(2.*pi)**p3) w errorbars lc rgb word(cols,i) pt 7 ti word(sim_names,i),\
     for [i=1:nsim] sim_file(word(fields,4),word(fields,1),word(sims,i)) u ($1*disp(i,nsim)):(f*($1**p1)*(($1+1)**p2)*$2/(2.*pi)**p3):(f*($1**p1)*(($1+1)**p2)*$3/(2.*pi)**p3) w errorbars lc rgb word(cols,i) pt 7 noti,\
     for [i=1:nsim] sim_file(word(fields,5),word(fields,1),word(sims,i)) u ($1*disp(i,nsim)):(f*($1**p1)*(($1+1)**p2)*$2/(2.*pi)**p3):(f*($1**p1)*(($1+1)**p2)*$3/(2.*pi)**p3) w errorbars lc rgb word(cols,i) pt 7 noti,\
     for [i=1:nsim] sim_file(word(fields,1),word(fields,2),word(sims,i)) u ($1*disp(i,nsim)):(($1**p1)*(($1+1)**p2)*$2/(2.*pi)**p3):(($1**p1)*(($1+1)**p2)*$3/(2.*pi)**p3)     w errorbars lc rgb word(cols,i) pt 7 noti,\
     for [i=1:nsim] sim_file(word(fields,2),word(fields,3),word(sims,i)) u ($1*disp(i,nsim)):(($1**p1)*(($1+1)**p2)*$2/(2.*pi)**p3):(($1**p1)*(($1+1)**p2)*$3/(2.*pi)**p3)     w errorbars lc rgb word(cols,i) pt 7 noti,\
     for [i=1:nsim] sim_file(word(fields,2),word(fields,4),word(sims,i)) u ($1*disp(i,nsim)):(($1**p1)*(($1+1)**p2)*$2/(2.*pi)**p3):(($1**p1)*(($1+1)**p2)*$3/(2.*pi)**p3)     w errorbars lc rgb word(cols,i) pt 7 noti,\
     for [i=1:nsim] sim_file(word(fields,2),word(fields,5),word(sims,i)) u ($1*disp(i,nsim)):(($1**p1)*(($1+1)**p2)*$2/(2.*pi)**p3):(($1**p1)*(($1+1)**p2)*$3/(2.*pi)**p3)     w errorbars lc rgb word(cols,i) pt 7 noti
