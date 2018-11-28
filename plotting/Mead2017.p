reset

if(!exists('print')){print=0}
if(print==0) {set term aqua dashed}
if(print==1) {set term post enh col; set output 'Mead2017.eps'}

power(name,z)=sprintf('data/power_%s_z%d.dat',name,z)
sim_power(z,r,name)=sprintf('/Users/Mead/Physics/fixed_sigma/power/L200/z%d/R%d/%s_power.dat',z,r,name)
sim2_power(z,r,name1,name2)=sprintf('/Users/Mead/Physics/fixed_sigma/power/L200/z%d/R%d/%s_%s_power.dat',z,r,name1,name2)
sim_ratio(z,name)=sprintf('/Users/Mead/Physics/fixed_sigma/power/L200/z%d/%s_ratio.dat',z,name)


# Initial white space
print ''

# Plot choice
print 'iplot = 1: z = 0'
print 'iplot = 2: z = 1'
if(!exists('iplot')){iplot=1}
print 'iplot = ', iplot
print ''

# x range
kmin=0.02
kmax=10
set log x
set xlabel 'k / h Mpc^{-1}'
set xrange [kmin:kmax]

# y range
dr=0.2
rmin=0.9
rmax=1.16
set yrange[rmin:rmax]
set ylabel 'P(k) / P_{{/Symbol L}CDM}(k)'

base='LCDM'
models="'OCDM' 'w-0.7' 'w-1.3' 'wa0.5' 'wa-0.5' 'w0-0.7wa-1.5' 'w0-1.3wa0.5' 'SCDM'"
bases="'LCDM_OCDM' 'LCDM_w-0.7' 'LCDM_w-1.3' 'LCDM_wa0.5' 'LCDM_wa-0.5' 'LCDM_w0-0.7wa-1.5' 'LCDM_w0-1.3wa0.5' 'LCDM_SCDM'"
model_names="'Open' 'w = -0.7' 'w = -1.3' 'w_a = 0.5' 'w_a = -0.5' 'w = -0.7; w_a = -1.5' 'w = -1.3; w_a = 0.5' 'EdS'"

set key top left

if(iplot==1){

z=0

plot 1 w l lt -1 noti,\
     for [i=1:words(models)] '<paste '.power(word(models,i),z).' '.power(base,z).'' u 1:($5/$10) w l lc i lw 2 ti word(model_names,i),\
     for [i=1:words(models)] sim_ratio(z,word(models,i)) u 1:2:3 w errorbars pt 7 ps .5 dt 1 lc i noti

}

if(iplot==2){

r=1
z=1

plot 1 w l lt -1 noti,\
     for [i=1:words(models)] '<paste '.power(word(models,i),z).' '.power(word(bases,i),z).''                    u 1:($5/$10) w l lc i lw 2 ti word(model_names,i),\
     for [i=1:words(models)] '<paste '.sim_power(z,r,word(models,i)).' '.sim2_power(z,r,base,word(models,i)).'' u 1:($2/$6) pt 7 ps .5 dt 1 lc i noti

}
