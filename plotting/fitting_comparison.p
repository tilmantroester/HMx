unset multiplot
reset

# Set terminal
if(!exists('print')) {print=0}
if(print==0) {set term aqua dashed dl 1 font ',10'}
if(print==1) {set term post enh col font ',9'}

# File
power(dir,mod,z,n,m,c,f1,f2,iz)=sprintf('%s/%s%s_n%d_m%d_c%d_best_cos1_%d%d_z%d.dat',dir,mod,z,n,m,c,f1,f2,iz)
dmonly(z)=sprintf('data/power%s_HMcode_BAHAMAS_M1536_k.dat',z) # This should be an HMcode prediction

# Sets of names of all models
mods="'AGN_7p6_nu0' 'AGN_TUNED_nu0_v2' 'AGN_TUNED_nu0_v3' 'AGN_TUNED_nu0_v4' 'AGN_TUNED_nu0' 'AGN_8p0_nu0'"
names="'AGN 7p6' 'AGN tuned v2' 'AGN tuned v3' 'AGN tuned v4' 'AGN tuned' 'AGN 8p0'"
cols="'dark-yellow' 'light-blue' 'light-blue' 'light-blue' 'blue' 'dark-plum'"

# Names of the subset of few models
mods_few="'AGN_7p6_nu0' 'AGN_TUNED_nu0' 'AGN_8p0_nu0'"
names_few="'AGN 7p6' 'AGN tuned' 'AGN 8p0'"
cols_few="'dark-yellow' 'blue' 'dark-plum'"

# Different fields
fields="'matter' 'cdm' 'gas' 'stars' 'epressure'"
field_names="'matter' 'CDM' 'gas' 'stars' 'electon pressure'"

# Initial white space
print ''

# Choose plot to make
if(!exists('iplot')){iplot=1}
print 'iplot = 1: Power comparison'
print 'iplot = 2: Response comparison'
print 'iplot = ', iplot
print ''

# Set the directory
if(!exists('dir')) {dir='fitting'}
print 'Directory: ', dir

# Number of points in chain
if(!exists('n')) {n=100000}
print 'Number of points in chain: n: ', n

# Fitting mode
if(!exists('m')) {m=17}
print 'Fitting mode: m: ', m

if(m==16) {zsep=1}
if(m==17 || m==18 || m==19) {zsep=0}

# Chain
if(!exists('c')) {c=1}
print 'Chain: c: ', c

# Fields
if(!exists('f1')) {f1=1}
if(!exists('f2')) {f2=1}
if(iplot==1){
print 'field 1: f1: ', f1, ' ', word(fields,f1), ' ', word(field_names,f1)
print 'field 2: f2: ', f2, ' ', word(fields,f2), ' ', word(field_names,f2)
}

# Set ouput file names
if(print==1 && iplot==1) {outfile(m,f1,f2)=sprintf('fitting_power_m%d_%s-%s.eps',m,f1,f2)}
if(print==1 && iplot==2) {outfile(m,f1,f2)=sprintf('fitting_response_m%d_%s-%s.eps',m,f1,f2)}

# Print an example file to screen to check
print 'Example file: ', power(dir,word(mods,2),'_z0.0',n,m,c,f1,f2,1)
print ''

# Location of zlabel on plot
labx=0.88
laby=0.9

# Line and point style and size
point_size=0.4
line_width=2.0
point_type=7
line_type=1

# Set x range
kmin=2e-2
kmax=1e1
set log x
set xlabel 'k / h Mpc^{-1}'
set xrange [kmin:kmax]

# Power plots (Delta^2(k)/k^1.5)
if(iplot==1){
rmin=1
rmax=25
set ylabel '{/Symbol D}^2(k) / (k / h Mpc^{-1})^{1.5}'
#set yrange [rmin:rmax]
tit=''.word(field_names,f1).'-'.word(field_names,f2).''
}

# Power suppression plots
if(iplot==2){
#rmin=0.75
#rmax=1.02
set ylabel 'P(k) / P_{dmonly}(k)'
#set yrange [rmin:rmax]
set yrange [*:*]
}

if(print==1 && iplot==1) {set output outfile(m,word(fields,f1),word(fields,f2))}
if(print==1 && iplot==2) {set output outfile(m,word(fields,f1),word(fields,f2))}
show output

# Key position
set key top left opaque

# Multiplot over z
set multiplot layout 2,2 title tit

# Margins
set lmargin 10
set rmargin 4
set tmargin 2
set bmargin 3

# Rectangle
set style rect fc lt -1 fs solid 0.15 noborder

# Loop over redshifts
do for [i=1:4] {

kmin_fit=0.15
if(i==1) {z='_z0.0'; iz=1; zlab='z = 0.0'; kmax_fit=10.}
if(i==2) {z='_z0.5'; iz=2; zlab='z = 0.5'; kmax_fit=4.}
if(i==3) {z='_z1.0'; iz=3; zlab='z = 1.0'; kmax_fit=2.}
if(i==4) {z='_z2.0'; iz=4; zlab='z = 2.0'; kmax_fit=1.}
za=z
if(zsep==1){izz=1}
if(zsep==0){za=''; izz=iz}

# If I need to do cross spectra I could probably do pmax=sqrt(pmax1*pmax2) or similar

# Matter
if(f1==1 && f2==1 && iz==1) {pmax=25}
if(f1==1 && f2==1 && iz==2) {pmax=12}
if(f1==1 && f2==1 && iz==3) {pmax=6}
if(f1==1 && f2==1 && iz==4) {pmax=2}

# CDM
if(f1==2 && f2==2 && iz==1) {pmax=20}
if(f1==2 && f2==2 && iz==2) {pmax=10}
if(f1==2 && f2==2 && iz==3) {pmax=5}
if(f1==2 && f2==2 && iz==4) {pmax=1.5}

# Gas
if(f1==3 && f2==3 && iz==1) {pmax=0.4}
if(f1==3 && f2==3 && iz==2) {pmax=0.2}
if(f1==3 && f2==3 && iz==3) {pmax=0.14}
if(f1==3 && f2==3 && iz==4) {pmax=0.06}

# Stars
if(f1==4 && f2==4 && iz==1) {pmax=0.03}
if(f1==4 && f2==4 && iz==2) {pmax=0.014}
if(f1==4 && f2==4 && iz==3) {pmax=0.007}
if(f1==4 && f2==4 && iz==4) {pmax=0.0016}

# Electron pressure
if(f1==1 && f2==5 && iz==1) {pmax=0.006}
if(f1==1 && f2==5 && iz==2) {pmax=0.0018}
if(f1==1 && f2==5 && iz==3) {pmax=0.0006}
if(f1==1 && f2==5 && iz==4) {pmax=8e-5}

set label zlab at graph labx,laby

unset key
if(i==1) {set key top left}

set obj rect from kmin,-1     to kmin_fit,pmax
set obj rect from kmax_fit,-1 to kmax,pmax

if(iplot==1){

# Set yrange (annoying large scatter at high k makes range crap)
set yrange [*:*]
if(m==16 || m==17 || m==18) {set yrange [0:pmax]}

plot for [i=1:words(mods)] NaN w l lw 2 lc rgb word(cols,i) ti word(names,i),\
     for [i=1:words(mods)] power(dir,word(mods,i),za,n,m,c,f1,f2,izz) u 1:($3/$1**1.5) w p pt point_type ps point_size lc rgb word(cols,i) noti,\
     for [i=1:words(mods)] power(dir,word(mods,i),za,n,m,c,f1,f2,izz) u 1:($2/$1**1.5) w l dt line_type  lw line_width lc rgb word(cols,i) noti

}

if(iplot==2){

print 'Note that you must have a DMONLY file for the correct k spacing'
print 'DMONLY file: ', dmonly(z)
print ''

plot for [i=1:words(mods_few)] '<paste '.power(dir,word(mods_few,i),za,n,m,c,f1,f2,izz).' '.dmonly(z).'' u 1:($3/$8) w p ps point_size pt point_type lc rgb word(cols_few,i) noti,\
     for [i=1:words(mods)]     '<paste '.power(dir,word(mods,i),za,n,m,c,f1,f2,izz).'     '.dmonly(z).'' u 1:($2/$8) w l lw line_width dt line_type  lc rgb word(cols,i) noti

}

#unset obj
unset obj

unset label

}

unset multiplot
