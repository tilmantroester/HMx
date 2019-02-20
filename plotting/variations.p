unset multiplot
reset

# Font file
cmsy='/Users/Mead/Fonts/cmsy10.pfb'

# Terminal information
if(!exists("print")){print=0}
if(print==0){set term aqua dashed font ',10' size 800,900; sun='sun'}
if(print==1){set term post enh col dashed dl 1 fontfile cmsy font ',7' size 8,9; sun='{/cmsy10 \014}'}

# Initial white space
print ''

# File paths
dmonly='data/DMONLY.dat'
power(param,value,type)=sprintf('data/power_param_%i_value_%i_%s.dat',param,value,type)
profile(mass,param,value)=sprintf('data/profile_mass_%i_param_%i_value_%i.dat',mass,param,value)
mass_fraction(param,value)=sprintf('data/mass_fractions_param_%i_value_%i.dat',param,value)
UPP(mass)=sprintf('/Users/Mead/Physics/HMx/diagnostics/UPP/halo_profile_m%i.dat',mass)

# Fix the parameter to plot
if(!exists("param")){param=1}
if(param==1)  {pname='{/Symbol a}';       min=0.0333; max=0.6333; ilog=0; ptype='mid'}
if(param==2)  {pname='{/Symbol e}';       min=0.5;    max=1.5;    ilog=0; ptype='mid'}
if(param==3)  {pname='{/Symbol G}';       min=1.12;   max=1.22;   ilog=0; ptype='mid'}
if(param==4)  {pname='M_g / M_{'.sun.'}'; min=1e13;   max=1e15;   ilog=1; ptype='mid'}
if(param==5)  {pname='A_*';               min=0.02;   max=0.04;   ilog=0; ptype='mid'}
if(param==6)  {pname='T_{WHIM} / K';      min=1e5;    max=1e7;    ilog=1; ptype='mid'}
if(param==7)  {pname='c_*';               min=5.;     max=20.;    ilog=1; ptype='mid'}
if(param==8)  {pname='f_{c}';             min=0.;     max=0.25;   ilog=0; ptype='low'}
if(param==9)  {pname='M_*';               min=5e11;   max=5e13;   ilog=1; ptype='mid'}
if(param==10) {pname='{/Symbol s}_*';     min=0.7;    max=1.7;    ilog=0; ptype='mid'}
if(param==11) {pname='{/Symbol a}^{P}';   min=-0.1;   max=0.1;    ilog=0; ptype='mid'}
if(param==12) {pname='{/Symbol G}^{P}';   min=-0.01;  max=0.01;   ilog=0; ptype='mid'}
if(param==13) {pname='c_*^{P}';           min=-0.1;   max=0.1;    ilog=0; ptype='mid'}
if(param==14) {pname='f_{h}';             min=0.;     max=0.25;   ilog=0; ptype='low'}
if(param==15) {pname='{/Symbol a}^{Z}';   min=-0.1;   max=0.1;    ilog=0; ptype='mid'}
if(param==16) {pname='{/Symbol G}^{Z}';   min=-0.1;   max=0.1;    ilog=0; ptype='mid'}
if(param==17) {pname='M_g^{Z}';           min=-0.1;   max=0.1;    ilog=0; ptype='mid'}
if(param==18) {pname='A_*^{Z}';           min=-0.1;   max=0.1;    ilog=0; ptype='mid'}
if(param==19) {pname='T_{WHIM}^{Z}';      min=-0.1;   max=0.1;    ilog=0; ptype='mid'}
if(param==20) {pname='{/Symbol h}';       min=-0.5;   max=0.0;    ilog=0; ptype='hig'}
if(param==33) {pname='{/Symbol e}^{Z}';   min=-0.1;   max=0.1;    ilog=0; ptype='mid'}
if(param==34) {pname='{/Symbol b}';       min=0.0333; max=0.6333; ilog=0; ptype='mid'}
if(param==35) {pname='{/Symbol b}^{P}';   min=-0.1;   max=0.1;    ilog=0; ptype='mid'}
if(param==36) {pname='{/Symbol b}^{Z}';   min=-0.1;   max=0.1;    ilog=0; ptype='mid'}
print 'Comparison for parameter (set with *param*): '.param.''

# Output figure
outfile(i)=sprintf('paper/variations_param_%i.eps',i)
if(print==1){set output outfile(param); print 'Output: ', outfile(param)}

# Number of different values for the parameters
n=9
if(ptype eq 'mid') {m=(n+1)/2}
if(ptype eq 'low') {m=1}
if(ptype eq 'hig') {m=n}

if(!exists('z')){z=0.0}
if(z==0.0) snap='snap32'
if(z==0.5) snap='snap28'
if(z==1.0) snap='snap26'
if(z==2.0) snap='snap22'
print 'Redshift: ', z
print 'Snapshot: ', snap

if(!exists('sim')){sim=4}
simulation_names="'DMONLY_2fluid' 'AGN_7p6' 'AGN_8p0' 'AGN_TUNED'"
simulation_titles="'DMONLY' 'AGN-lo' 'AGN-hi' 'AGN'"
simulation(sim,snap,type1,type2)=sprintf('/Users/Mead/Physics/BAHAMAS/power/M1024/%s_nu0_L400N1024_WMAP9_%s_%s_%s_power.dat',sim,snap,type1,type2)
sim_name=word(simulation_names,sim)
#sim_title(name,z)=''.word(simulation_titles,sim).'; z = '.z.''
name=word(simulation_titles,sim)
sim_title(name,z)=sprintf('%s; z = %1.1f',name,z)
print 'Comparing to simulation (set with *sim*): '.sim_title(name,z).''
dmsim=word(simulation_names,1)

# Names of simulation power spectrum types
dsim='all'
csim='dm'
gsim='gas'
ssim='stars'
psim='epressure'

# Name of halo-model power spectrum types
#types_density="'00' '01' '11' '02' '22' '03' '33''"
types_density="'00' '01' '02' '03'"
types_pressure="'00' '06' '66'"

# Label information
#ddlab = 'mm'
#dglab = 'mg'
#dslab = 'ms'
#dplab = 'mp'
#pplab = 'pp'
ddlab = 'i,j = m,m'
dglab = 'i,j = m,g'
dslab = 'i,j = m,s'
dplab = 'i,j = m,p'
pplab = 'i,j = p,p'

# power column for halo model power spectra files
cp=5
Lp=5

# power columns for simulation power spectra files
cs=2
ss=3
Ls=5

# Masses to plot
m1=13
m2=14
m3=15
mlab(m)=sprintf('10^{%i} M_{'.sun.'} / h',m)
mlabx=0.3
mlaby=0.94

# Offsets for axes from whole numbers
dx=1.05
dy=1.05

# Function for progressions in either linear or log space
if(ilog==0){prog(a,b,i,n)=a+(b-a)*real(i-1)/real(n-1)}
if(ilog==1){prog(a,b,i,n)=exp(log(a)+log(b/a)*real(i-1.)/real(n-1.))}

#Final white space
print ''

### Power spectrum plots ###

# k axis range
kmin=0.01
kmax=10.

# density power axis range
d_pmin=1e-5
d_pmax=1e3

# pressure power axis range
p_pmin=1e-10
p_pmax=1e3

# Density ratio axis
d_ratio_min=5e-3
d_ratio_max=2.

# Pressure ratio axis
p_ratio_min=1e-9
p_ratio_max=2.

# Axes labels
klab='k / h Mpc^{-1}'
plab='{/Symbol D}_{i,j}^2(k)'
slab='P_{i,j}(k) / P_{DMONLY}(k)'

### ###

### Density profile plots ###

# Density profile axis
rhomin=1e-3
rhomax=1e1
rholab='4{/Symbol p} r^2 {/Symbol r}(r) / M [Mpc/h]^{-1}'

# Pressure profile axis
premin=1e-4
premax=1e2
prelab='4{/Symbol p} r^2 P(r) [eV (Mpc/h)^{-1}]'

# Axis range for halo profiles
rmin=0.
rmax=1.1
rlab='r / r_v'

### ###

### Mass fraction plots ###

# Axis range for halo mass fractions
massmin=1e10
massmax=1e16
masslab='M / h^{-1} M_{'.sun.'}'

# Mass fraction range
fmin=1e-3
fmax=2
flab='Mass fraction'

### ###

### Plot positions ###

# Global boundaries of plot
all_top = 0.99
all_bottom = 0.06
all_left = 0.06
all_right = 0.90

# Local boundaries of plots
first_column_right = 0.47
second_column_left = 0.54
second_column_break_bottom = 0.41
second_column_break_top = 0.35
color_bar_size = 0.02

# Density power spectrum
power_density_top = all_top-0.*(all_top-all_bottom)/4.
power_density_bottom = all_top-1.*(all_top-all_bottom)/4.
power_density_left = all_left
power_density_right = first_column_right

# Density power spectrum ratio
ratio_density_top = power_density_bottom
ratio_density_bottom = all_top-2.*(all_top-all_bottom)/4.
ratio_density_left = power_density_left
ratio_density_right = power_density_right

# Pressure power spectrum
power_pressure_top = ratio_density_bottom
power_pressure_bottom = all_top-3.*(all_top-all_bottom)/4.
power_pressure_left = power_density_left
power_pressure_right = power_density_right

# Pressure power spectrum ratio
ratio_pressure_top = power_pressure_bottom
ratio_pressure_bottom = all_top-4.*(all_top-all_bottom)/4.
ratio_pressure_left = power_density_left
ratio_pressure_right = power_density_right

# Number of halo-mass profiles
nrho=3

# Left-hand density profile
rho1_top = all_top
rho1_bottom = all_top-(all_top-second_column_break_bottom)/2.
rho1_left = second_column_left
rho1_right = rho1_left+(all_right-rho1_left)/real(nrho)

# Middle density profile
rho2_top = rho1_top
rho2_bottom = rho1_bottom
rho2_left = rho1_right
rho2_right = rho1_left+2.*(all_right-rho1_left)/real(nrho)

# Right-hand density profile
rho3_top = rho1_top
rho3_bottom = rho1_bottom
rho3_left = rho2_right
rho3_right = rho1_left+3.*(all_right-rho1_left)/real(nrho)

# Left-hand pressure profile
pre1_top = rho1_bottom
pre1_bottom = second_column_break_bottom
pre1_left = rho1_left
pre1_right = rho1_right

# Middle pressure profile
pre2_top = rho2_bottom
pre2_bottom = pre1_bottom
pre2_left = rho2_left
pre2_right = rho2_right

# Right-hand pressure profile
pre3_top = rho3_bottom
pre3_bottom = pre1_bottom
pre3_left = rho3_left
pre3_right = rho3_right

# Mass fractions
mfrac_top = second_column_break_top
mfrac_bottom = all_bottom
mfrac_left = rho1_left
mfrac_right = rho3_right

### ###

# Colour palette - needs to be set after plot position stuff
if(ptype eq 'mid') {set palette defined (1 'light-blue', 2 'black', 3 'light-red')}
if(ptype eq 'low') {set palette defined (1 'black', 2 'light-red')}
if(ptype eq 'hig') {set palette defined (1 'light-blue', 2 'black')}
set cbrange [min:max]
if(ilog==0){unset log cb; set format cb}
if(ilog==1){set log cb; set format cb '10^{%T}'}
set colorbox vertical user origin all_right+0.02, all_bottom size color_bar_size, all_top-all_bottom
set cblabel pname

set multiplot

### Density: Power spectrum plots - Top left ###

set tmargin at screen power_density_top
set bmargin at screen power_density_bottom
set lmargin at screen power_density_left
set rmargin at screen power_density_right

set xlabel ''
set format x ''
set log x
set xrange [kmin:kmax]

set log y
set ylabel plab
set format y '10^{%T}'
set yrange[d_pmin*dy:d_pmax/dy]

set key top left

set label ddlab at graph 0.02,0.35
set label dglab at graph 0.02,0.25
set label dslab at graph 0.02,0.15

plot for [j=1:words(types_density)] for [i=1:n] power(param,i,word(types_density,j)) u 1:(column(cp)):(prog(min,max,i,n)) w l lw 2 lc palette noti,\
     for [j=1:words(types_density)] power(param,m,word(types_density,j)) u 1:(column(cp)):(prog(min,max,m,n)) w l lw 2 lc palette noti,\
     simulation(sim_name,snap,dsim,dsim) u 1:($2-$3):5 w e pt 7 ps .3 lc 'black' ti sim_title(name,z),\
     simulation(sim_name,snap,dsim,gsim) u 1:($2-$3):5 w e pt 7 ps .3 lc 'black' noti,\
     simulation(sim_name,snap,dsim,ssim) u 1:($2-$3):5 w e pt 7 ps .3 lc 'black' noti

unset colorbox
unset label

### ###

### Density: Power spectrum residual ###

set tmargin at screen ratio_density_top
set bmargin at screen ratio_density_bottom
set lmargin at screen ratio_density_left
set rmargin at screen ratio_density_right

#set xlabel klab
#set format x
set xlabel ''
set format x ''

#set ylabel 'P_{i,j}(k) / P_{DMONLY}(k)'
set ylabel slab
set yrange [d_ratio_min*dy:d_ratio_max/dy]

set key top left

set label ddlab at graph 0.03,0.93
set label dglab at graph 0.03,0.63
set label dslab at graph 0.03,0.27

plot for [j=1:words(types_density)] for [i=1:n] '<paste '.power(param,i,word(types_density,j)).' '.dmonly.'' u 1:(column(cp)/column(cp+Lp)):(prog(min,max,i,n)) w l lw 2 lc palette noti,\
     for [j=1:words(types_density)] '<paste '.power(param,m,word(types_density,j)).' '.dmonly.'' u 1:(column(cp)/column(cp+Lp)):(prog(min,max,m,n)) w l lw 2 lc palette noti,\
     '<paste '.simulation(sim_name,snap,dsim,dsim).' '.simulation(dmsim,snap,dsim,dsim).'' u 1:((column(cs)-column(ss))/(column(cs+Ls)-column(ss+Ls))) w p pt 7 ps .3 lc 'black' noti,\
     '<paste '.simulation(sim_name,snap,dsim,gsim).' '.simulation(dmsim,snap,dsim,dsim).'' u 1:((column(cs)-column(ss))/(column(cs+Ls)-column(ss+Ls))) w p pt 7 ps .3 lc 'black' noti,\
     '<paste '.simulation(sim_name,snap,dsim,ssim).' '.simulation(dmsim,snap,dsim,dsim).'' u 1:((column(cs)-column(ss))/(column(cs+Ls)-column(ss+Ls))) w p pt 7 ps .3 lc 'black' noti

unset label

### ###

### Pressure: power spectrum - Bottom right ###

set tmargin at screen power_pressure_top
set bmargin at screen power_pressure_bottom
set lmargin at screen power_pressure_left
set rmargin at screen power_pressure_right

set xlabel ''
set format x ''
set log x
set xrange [kmin:kmax]

set log y
set ylabel plab
set format y '10^{%T}'
set yrange[p_pmin*dy:p_pmax/dy]

set label ddlab at graph 0.03,0.45
set label dplab at graph 0.03,0.25
set label pplab at graph 0.03,0.10

set key top left

plot for [j=1:words(types_pressure)] for [i=1:n] power(param,i,word(types_pressure,j)) u 1:(column(cp)):(prog(min,max,i,n)) w l lw 2 lc palette noti,\
     for [j=1:words(types_pressure)] power(param,m,word(types_pressure,j)) u 1:(column(cp)):(prog(min,max,m,n)) w l lw 2 lc palette noti,\
     simulation(sim_name,snap,dsim,dsim) u 1:($2-$3):5 w e pt 7 ps .3 lc 'black' ti sim_title(name,z),\
     simulation(sim_name,snap,dsim,psim) u 1:($2-$3):5 w e pt 7 ps .3 lc 'black' noti,\
     simulation(sim_name,snap,psim,psim) u 1:($2-$3):5 w e pt 7 ps .3 lc 'black' noti

unset colorbox
unset label

### ###

### Pressure: power spectrum residual ###

set tmargin at screen ratio_pressure_top
set bmargin at screen ratio_pressure_bottom
set lmargin at screen ratio_pressure_left
set rmargin at screen ratio_pressure_right

set xlabel klab
set format x

#set ylabel 'P_{i,j}(k) / P_{DMONLY}(k)'
set ylabel slab
set yrange [p_ratio_min*dy:p_ratio_max/dy]

set label ddlab at graph 0.03,0.9
set label dplab at graph 0.03,0.6
set label pplab at graph 0.03,0.3

plot for [j=1:words(types_pressure)] for [i=1:n] '<paste '.power(param,i,word(types_pressure,j)).' '.dmonly.'' u 1:(column(cp)/column(cp+Lp)):(prog(min,max,i,n)) w l lw 2 lc palette noti,\
     for [j=1:words(types_pressure)] '<paste '.power(param,m,word(types_pressure,j)).' '.dmonly.'' u 1:(column(cp)/column(cp+Lp)):(prog(min,max,m,n)) w l lw 2 lc palette noti,\
     '<paste '.simulation(sim_name,snap,dsim,dsim).' '.simulation(dmsim,snap,dsim,dsim).'' u 1:((column(cs)-column(ss))/(column(cs+Ls)-column(ss+Ls))) w p pt 7 ps .3 lc 'black' noti,\
     '<paste '.simulation(sim_name,snap,dsim,psim).' '.simulation(dmsim,snap,dsim,dsim).'' u 1:((column(cs)-column(ss))/(column(cs+Ls)-column(ss+Ls))) w p pt 7 ps .3 lc 'black' noti,\
     '<paste '.simulation(sim_name,snap,psim,psim).' '.simulation(dmsim,snap,dsim,dsim).'' u 1:((column(cs)-column(ss))/(column(cs+Ls)-column(ss+Ls))) w p pt 7 ps .3 lc 'black' noti

unset label

### ###

### Density: Density profiles  - Top right ###

unset log x
set xrange[rmin:rmax/dx]
set xlabel ''
set format x ''

set yrange [rhomin*dy:rhomax/dy]

do for [j=1:nrho]{

if(j==1){
set tmargin at screen rho1_top
set bmargin at screen rho1_bottom
set lmargin at screen rho1_left
set rmargin at screen rho1_right
mass=m1
set ylabel rholab #offset 2
set label 'CDM' at graph 0.07,0.9
set label 'gas' at graph 0.07,0.5
}

if(j==2){
unset label
set tmargin at screen rho2_top
set bmargin at screen rho2_bottom
set lmargin at screen rho2_left
set rmargin at screen rho2_right
mass=m2
set ylabel ''
set format y ''
}

if(j==3){
set tmargin at screen rho3_top
set bmargin at screen rho3_bottom
set lmargin at screen rho3_left
set rmargin at screen rho3_right
mass=m3
set ylabel ''
set format y ''
}

set label mlab(mass) at graph mlabx,mlaby

plot for [c=2:5] for [i=1:n] profile(mass,param,i) u 1:(4.*pi*$1*$1*column(c)):(prog(min,max,i,n)) w l lw 2 dt 1 lc palette noti,\
     for [c=2:5] profile(mass,param,m) u 1:(4.*pi*$1*$1*column(c)):(prog(min,max,m,n)) w l lw 2 dt 1 lc palette noti

unset label

}

### Pressure profiles - middle right ###

### ###

set xrange[rmin:rmax/dx]
set xlabel rlab
#set format x '10^{%T}'
set format x

set yrange [premin*dy:premax/dy]

do for [j=1:nrho]{

if(j==1){
set tmargin at screen pre1_top
set bmargin at screen pre1_bottom
set lmargin at screen pre1_left
set rmargin at screen pre1_right
mass=m1
set format y '10^{%T}'
set ylabel prelab# offset 2
#set label 'CDM' at graph 0.03,0.9
}

if(j==2){
unset label
set tmargin at screen pre2_top
set bmargin at screen pre2_bottom
set lmargin at screen pre2_left
set rmargin at screen pre2_right
mass=m2
set ylabel ''
set format y ''
}

if(j==3){
set tmargin at screen pre3_top
set bmargin at screen pre3_bottom
set lmargin at screen pre3_left
set rmargin at screen pre3_right
mass=m3
set ylabel ''
set format y ''
}

#set label mlab(mass) at graph mlabx,mlaby

plot for [i=1:n] profile(mass,param,i) u 1:(4.*pi*$1*$1*column(6)):(prog(min,max,i,n)) w l lw 2 dt 1 lc palette noti,\
     profile(mass,param,m) u 1:(4.*pi*$1*$1*column(6)):(prog(min,max,m,n)) w l lw 2 dt 1 lc palette noti

unset label

}

###

### Halo mass fractions - bottom right ###

set tmargin at screen mfrac_top
set bmargin at screen mfrac_bottom
set lmargin at screen mfrac_left
set rmargin at screen mfrac_right

set xrange[massmin:massmax]
set xlabel masslab
set log x
set format x '10^{%T}'

set log y
set yrange [fmin:fmax]
set ylabel flab
set format y '10^{%T}'

set label 'CDM'   at graph 0.03,0.93
set label 'gas'   at graph 0.03,0.73
set label 'stars' at graph 0.03,0.35

plot for [c=2:6] for [i=1:n] mass_fraction(param,i) u 1:(column(c)):(prog(min,max,i,n)) w l lw 2 dt 1 lc palette noti,\
     for [c=2:6] mass_fraction(param,m) u 1:(column(c)):(prog(min,max,m,n)) w l lw 2 dt 1 lc palette noti

unset label

### ###

unset multiplot
