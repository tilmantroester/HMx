unset multiplot
reset

cmsy='/Users/Mead/Fonts/cmsy10.pfb'

if(!exists("print")){print=0}
if(print==0){set term aqua dashed font ',10'; sun='sun'}
if(print==1){set term post enh col dashed dl 1 fontfile cmsy font ',10'; sun='{/cmsy10 \014}'}

dmonly='variations/DMONLY.dat'
power(param,value,type)=sprintf('variations/power_param_%i_value_%i_%s.dat',param,value,type)
profile(mass,param,value)=sprintf('variations/profile_mass_%i_param_%i_value_%i.dat',mass,param,value)
mass_fraction(param,value)=sprintf('variations/mass_fractions_param_%i_value_%i.dat',param,value)
UPP(mass)=sprintf('/Users/Mead/Physics/HMx/diagnostics/UPP/halo_profile_m%i.dat',mass)

#Simulation pressure conversion factor
fm=6.242e11

#Fix the parameter to plot
if(!exists("param")){param=1}
if(param==1){pname='{/Symbol a}';       min=0.4;   max=2.0;   ilog=0; coll='light-blue'}
if(param==2){pname='{/Symbol b}';       min=0.0;   max=2.0;   ilog=0; coll='pink'}
if(param==3){pname='{/Symbol G}';       min=1.10;  max=1.26;  ilog=0; coll='orange'}
if(param==4){pname='M_B / M_{'.sun.'}'; min=1e13;  max=1e15;  ilog=1; coll='light-green'}
if(param==5){pname='A_*';               min=0.015; max=0.055; ilog=0; coll='gold'}

#Output figure
if(!exists("type")){type='matter'}
if(type eq 'matter'){outfile(i)=sprintf('variations_matter_param_%i.eps',i)}
if(type eq 'pressure'){outfile(i)=sprintf('variations_pressure_param_%i.eps',i)}
if(print==1){set output outfile(param); print 'Output: ', outfile(param)}

#Number of different values of the parameter
n=9

#Simulation conversion factor
#fm=6.242e11

#Cosmological parameters
#Om_m=0.272
#Om_b=0.0455

#bf=Om_m/Om_b
#cf=Om_m/(Om_m-Om_b)

#bf2=bf**2
#cf2=cf**2

if(!exists("isim")){sim=4}
simulation_names='DMONLY REF NOCOOL_UVB AGN AGN_Theat_8p5 AGN_Theat_8p7'
simulation_titles="'DMONLY' 'REF' 'NO COOL' 'AGN' 'AGN 8.5' 'AGN 8.7'"
simulation(n,sim,type1,type2)=sprintf('/Users/Mead/Physics/cosmo-OWLS/power/N%i/%s_%s_%s_power.dat',n,sim,type1,type2)
sim_name=word(simulation_names,sim)
sim_title=word(simulation_titles,sim)
dmsim='DMONLY'
mesh=800
dsim='overdensity_grid'
csim='overdensity_grid_DM'
gsim='overdensity_grid_gas'
ssim='overdensity_grid_stars'
psim='pressure_grid'

types_density='dd dc cc dg gg ds ss'
types_pressure='dd dp pp'

ddlab='{/Symbol d}{/Symbol d}'
dplab='{/Symbol d}P'
pplab='PP'

#power column for halo model power spectra files
cp=5
Lp=5

#power columns for simulation power spectra files
cs=2
Ls=3

#Masses to plot
m1=13
m2=14
m3=15
mlab(m)=sprintf('10^{%i} h^{-1} M_{'.sun.'}',m)
mlabx=0.3
mlaby=0.94

if(ilog==0){prog(a,b,i,n)=a+(b-a)*real(i-1)/real(n-1)}
if(ilog==1){prog(a,b,i,n)=exp(log(a)+log(b/a)*real(i-1.)/real(n-1.))}

#pfac=1e2

kmin=0.01
kmax=10.
klab='k / h Mpc^{-1}'

pmin=1e-5
pmax=1e3
plab='{/Symbol D}_{i,j}^2(k)'

d_ratio_min=5e-3
d_ratio_max=2.

p_ratio_min=1e-4
p_ratio_max=2.

rmin=0.01
rmax=1e1
rlab='r / h^{-1} Mpc'

rhomin=1e-3
rhomax=1e1
rholab='4{/Symbol p} r^2 {/Symbol r}(r) / M'

massmin=1e10
massmax=1e16
masslab='M / h^{-1} M_{'.sun.'}'

fmin=1e-3
fmax=2
flab='Mass fraction'

dx=1.05
dy=1.05

all_top=0.98
all_bottom=0.1
all_left=0.05
all_right=0.91

power_top=all_top
power_bottom=0.5
power_left=all_left
power_right=0.43

ratio_top=power_bottom
ratio_bottom=all_bottom
ratio_left=power_left
ratio_right=power_right

gap=0.07

nrho=3

rho1_top=all_top
rho1_bottom=0.5
rho1_left=power_right+gap
rho1_right=rho1_left+(all_right-rho1_left)/real(nrho)

rho2_top=rho1_top
rho2_bottom=rho1_bottom
rho2_left=rho1_right
rho2_right=rho1_left+2.*(all_right-rho1_left)/real(nrho)

rho3_top=rho1_top
rho3_bottom=rho1_bottom
rho3_left=rho2_right
rho3_right=rho1_left+3.*(all_right-rho1_left)/real(nrho)

mfrac_top=rho1_bottom-0.1
mfrac_bottom=all_bottom
mfrac_left=rho1_left
mfrac_right=rho3_right

set palette defined (1 coll, 2 'black')
set cbrange [min:max]
if(ilog==0){unset log cb; set format cb}
if(ilog==1){set log cb; set format cb '10^{%T}'}
set colorbox vertical user origin all_right+0.02, .1 size .02,all_top-all_bottom
set cblabel pname

if(type eq 'matter'){

set multiplot

###Power spectrum plots###

set tmargin at screen power_top
set bmargin at screen power_bottom
set lmargin at screen power_left
set rmargin at screen power_right

set xlabel ''
set format x ''
set log x
set xrange [kmin:kmax]

set log y
set ylabel plab
set format y '10^{%T}'
set yrange[pmin*dy:pmax/dy]

set key top left

set label '{/Symbol d}{/Symbol d}' at graph 0.02,0.35
set label '{/Symbol d}g'           at graph 0.02,0.25
set label '{/Symbol d}s'           at graph 0.02,0.15

plot for [i=1:n] power(param,i,'dd') u 1:(column(cp)):(prog(min,max,i,n)) w l lw 2 lc palette noti,\
     for [i=1:n] power(param,i,'dg') u 1:(column(cp)):(prog(min,max,i,n)) w l lw 2 lc palette noti,\
     for [i=1:n] power(param,i,'ds') u 1:(column(cp)):(prog(min,max,i,n)) w l lw 2 lc palette noti,\
     simulation(mesh,sim_name,dsim,dsim) u 1:2 w p pt 2 lc 'black' ti sim_title,\
     simulation(mesh,sim_name,dsim,gsim) u 1:2 w p pt 2 lc 'black' noti,\
     simulation(mesh,sim_name,dsim,ssim) u 1:2 w p pt 2 lc 'black' noti,\

unset label

set tmargin at screen ratio_top
set bmargin at screen ratio_bottom
set lmargin at screen ratio_left
set rmargin at screen ratio_right

set xlabel klab
set format x

set ylabel 'P_{i,j}(k) / P_{DMONLY}(k)'
set yrange [d_ratio_min*dy:d_ratio_max/dy]

set key top left

set label '{/Symbol d}{/Symbol d}' at graph 0.03,0.93
set label '{/Symbol d}g'           at graph 0.03,0.63
set label '{/Symbol d}s'           at graph 0.03,0.27

plot for [i=1:n] '<paste '.power(param,i,'dd').' '.dmonly.'' u 1:(column(cp)/column(cp+Lp)):(prog(min,max,i,n)) w l lw 2 lc palette noti,\
     for [i=1:n] '<paste '.power(param,i,'dg').' '.dmonly.'' u 1:(column(cp)/column(cp+Lp)):(prog(min,max,i,n)) w l lw 2 lc palette noti,\
     for [i=1:n] '<paste '.power(param,i,'ds').' '.dmonly.'' u 1:(column(cp)/column(cp+Lp)):(prog(min,max,i,n)) w l lw 2 lc palette noti,\
'<paste '.simulation(mesh,sim_name,dsim,dsim).' '.simulation(mesh,dmsim,dsim,dsim).'' u 1:(column(cs)/column(cs+Ls))              w p pt 2 lc 'black' noti,\
'<paste '.simulation(mesh,sim_name,dsim,gsim).' '.simulation(mesh,dmsim,dsim,dsim).'' u 1:(column(cs)/column(cs+Ls))      w p pt 2 lc 'black' noti,\
'<paste '.simulation(mesh,sim_name,dsim,ssim).' '.simulation(mesh,dmsim,dsim,dsim).'' u 1:(column(cs)/column(cs+Ls)) w p pt 2 lc 'black' noti,\

unset label

#NaN w l lt -1 ti '{/Symbol d}{/Symbol d}',\
#NaN w l lt -1 ti '{/Symbol d}g',\
#NaN w l lt -1 ti '{/Symbol d}s'

### ###

### Density profiles ###

set xrange[rmin:rmax/dx]
set xlabel rlab
set format x

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

plot for [i=1:n] profile(mass,param,i) u 1:(4.*pi*$1*$1*column(2)):(prog(min,max,i,n)) w l lw 2 dt 1 lc palette noti,\
     for [i=1:n] profile(mass,param,i) u 1:(4.*pi*$1*$1*column(3)):(prog(min,max,i,n)) w l lw 2 dt 1 lc palette noti

unset label

}

### ###

### Halo mass fractions ###

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

plot for [i=1:n] mass_fraction(param,i) u 1:2:(prog(min,max,i,n)) w l lw 2 dt 1 lc palette noti,\
     for [i=1:n] mass_fraction(param,i) u 1:3:(prog(min,max,i,n)) w l lw 2 dt 1 lc palette noti,\
     for [i=1:n] mass_fraction(param,i) u 1:4:(prog(min,max,i,n)) w l lw 2 dt 1 lc palette noti,\
     for [i=1:n] mass_fraction(param,i) u 1:5:(prog(min,max,i,n)) w l lw 2 dt 2 lc palette noti,\
     for [i=1:n] mass_fraction(param,i) u 1:6:(prog(min,max,i,n)) w l lw 2 dt 2 lc palette noti

unset label

### ###

}

if(type eq 'pressure'){

set multiplot

###Power spectrum plots###

set tmargin at screen power_top
set bmargin at screen power_bottom
set lmargin at screen power_left
set rmargin at screen power_right

set xlabel ''
set format x ''
set log x
set xrange [kmin:kmax]

set log y
set ylabel plab
set format y '10^{%T}'
set yrange[pmin*dy:pmax/dy]

set label ddlab at graph 0.03,0.45
set label dplab at graph 0.03,0.25
set label pplab at graph 0.03,0.10

set key top left

plot for [j=1:words(types_pressure)] for [i=1:n] power(param,i,word(types_pressure,j)) u 1:((pfac**(j-1))*column(cp)):(prog(min,max,i,n)) w l lw 2 lc palette noti,\
     simulation(mesh,sim_name,dsim,dsim) u 1:2                 w p pt 2 lc 'black' ti sim_title,\
     simulation(mesh,sim_name,dsim,psim) u 1:($2*fm*pfac)      w p pt 2 lc 'black' noti,\
     simulation(mesh,sim_name,psim,psim) u 1:($2*(fm*pfac)**2) w p pt 2 lc 'black' noti

#NaN ti '{/Symbol d}{/Symbol d}' lt -1,\
#NaN ti '{/Symbol d}P'           lt -1,\
#NaN ti 'PP'                     lt -1

unset label

set tmargin at screen ratio_top
set bmargin at screen ratio_bottom
set lmargin at screen ratio_left
set rmargin at screen ratio_right

set xlabel klab
set format x

set ylabel 'P_{i,j}(k) / P_{DMONLY}(k)'
set yrange [p_ratio_min*dy:p_ratio_max/dy]

set label ddlab at graph 0.03,0.9
set label dplab at graph 0.03,0.6
set label pplab at graph 0.03,0.3

plot for [j=1:words(types_pressure)] for [i=1:n] '<paste '.power(param,i,word(types_pressure,j)).' '.dmonly.'' u 1:((pfac**(j-1))*column(cp)/column(cp+Lp)):(prog(min,max,i,n)) w l lw 2 lc palette noti,\
     '<paste '.simulation(mesh,sim_name,dsim,dsim).' '.simulation(mesh,dmsim,dsim,dsim).'' u 1:(column(cs)/column(cs+Ls))              w p pt 2 lc 'black' noti,\
     '<paste '.simulation(mesh,sim_name,dsim,psim).' '.simulation(mesh,dmsim,dsim,dsim).'' u 1:(column(cs)*fm*pfac/column(cs+Ls))      w p pt 2 lc 'black' noti,\
     '<paste '.simulation(mesh,sim_name,psim,psim).' '.simulation(mesh,dmsim,dsim,dsim).'' u 1:(column(cs)*(fm*pfac)**2/column(cs+Ls)) w p pt 2 lc 'black' noti

unset label

### ###

### Density profiles ###

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
set ylabel rholab# offset 2
set label 'CDM' at graph 0.03,0.9
set label 'gas' at graph 0.03,0.5
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

plot for [i=1:n] profile(mass,param,i) u 1:(4.*pi*$1*$1*column(2)):(prog(min,max,i,n)) w l lw 2 dt 1 lc palette noti,\
     for [i=1:n] profile(mass,param,i) u 1:(4.*pi*$1*$1*column(3)):(prog(min,max,i,n)) w l lw 2 dt 1 lc palette noti

unset label

}

#set ylab ''
#set format y ''

#set tmargin at screen rho2_top
#set bmargin at screen rho2_bottom
#set lmargin at screen rho2_left
#set rmargin at screen rho2_right

#set label mlab(m2) at graph mlabx,mlaby

#plot for [i=1:n] profile(m2,param,i) u 1:(4.*pi*$1*$1*column(2)):(prog(min,max,i,n)) w l lw 2 dt 1 lc palette noti,\
#     for [i=1:n] profile(m2,param,i) u 1:(4.*pi*$1*$1*column(3)):(prog(min,max,i,n)) w l lw 2 dt 1 lc palette noti

#unset label

### ###

### Pressure profiles ###

set xrange[rmin:rmax/dx]
set xlabel rlab
set format x

set yrange [premin*dy:premax/dy]

do for [j=1:nrho]{

if(j==1){
set tmargin at screen pre1_top
set bmargin at screen pre1_bottom
set lmargin at screen pre1_left
set rmargin at screen pre1_right
set ylabel prelab# offset 2
set format y '10^{%T}'
mass=m1
}

if(j==2){
set tmargin at screen pre2_top
set bmargin at screen pre2_bottom
set lmargin at screen pre2_left
set rmargin at screen pre2_right
set ylab ''
set format y ''
mass=m2
}

if(j==3){
set tmargin at screen pre3_top
set bmargin at screen pre3_bottom
set lmargin at screen pre3_left
set rmargin at screen pre3_right
set ylab ''
set format y ''
mass=m3
}


plot for [i=1:n] profile(mass,param,i) u 1:(4.*pi*$1*$1*column(7)):(prog(min,max,i,n)) w l lw 2 dt 1 lc palette noti,\
     UPP(mass) u 1:(4.*pi*$1*$1*column(7)) w l lw 2 lc -1 dt 2 noti

}

#set ylab ''
#set format y ''

#set tmargin at screen pre2_top
#set bmargin at screen pre2_bottom
#set lmargin at screen pre2_left
#set rmargin at screen pre2_right

#plot for [i=1:n] profile(m2,param,i) u 1:(4.*pi*$1*$1*column(7)):(prog(min,max,i,n)) w l lw 2 dt 1 lc palette noti,\
#     UPP(m2) u 1:(4.*pi*$1*$1*column(7)) w l lw 2 lc -1 dt 2 noti

### ###

unset multiplot

}

unset multiplot
