reset
unset multiplot

set term aqua size 1000,1000

# Initial white space
print ''

# Data file
if(!exists('chain')){chain='data/mcmc_good.dat'}
print 'Plotting data file: chain: ', chain
print ''

# Chain pruning
if(!exists('m')){m=0}
print 'Skipping entries: m: ', m
print ''

# Number of parameters
if(!exists('np')){np=8}
print 'Number of parameters: np: ', np
print ''

# Binning function
min=0.
n=1000
width = 1./real(n) # binwidth; evaluates to 1.0
bin(x) = min+width*(floor((x-min)/width)+0.5)

# Everything
if(im==11) {
np=8
labels="'log_{10}{/Symbol a}' 'log_{10}{/Symbol e}' '{/Symbol G}-1' 'log_{10} M_0' 'log_{10} A_*' 'log_{10} T_{w}' 'log_{10} c_*' 'log_{10} f_c'"
}

# Gas
if(im==12) {
np=6
labels="'log_{10}{/Symbol e}' '{/Symbol G}-1' 'log_{10} M_0' 'log_{10} A_*' 'log_{10} f_c' '{/Symbol G}^p'"
}

# Stars
if(im==13) {
np=5
labels="'log_{10}{A_*}' 'log_{10}{c_*}' 'c_*^p' 'log_{10}M_*' '{/Symbol s}_*'"
}

# Gas and stars
if(im==14) {

}

# Matter
if(im==15) {
np=3
labels="'log_{10}{/Symbol e}' '{/Symbol G}-1' 'log_{10} M_0' 'log_{10} A_*' 'log_{10} c_*' 'log_{10} f_c' '{/Symbol G}^p' 'c_*^p' 'log_{10}M_*' '{/Symbol s}_*'"
}

# Write parameters to screen
print 'Parameters to plot:'
do for [i=1:np] {print 'Parameter: ', i, ' ', word(labels,i)}
print ''

# Plot boundary stuff
top=0.98
bot=0.10
lef=0.10
rig=0.98
#nx=np-1
#ny=np-1
nx=np
ny=np
dx=(rig-lef)/ny
dy=(top-bot)/nx

# Now actually do the plotting
#set multiplot layout np-1,np-1
set multiplot# layout np,np

#do for [i=1:np-1]{
do for [i=1:np]{

do for [j=i:np]{

set lmargin at screen lef+(i-1)*dx
set rmargin at screen lef+(i-0)*dx

set tmargin at screen top-(j-1)*dy
set bmargin at screen top-(j-0)*dy

if(i==1)  {set format y; set ylabel word(labels,j); set ytics} else {set format y ''; set ylabel ''}
if(j==np) {set format x; set xlabel word(labels,i); set ytics} else {set format x ''; set xlabel ''}
if(i==j)  {set format y ''; set ylabel ''; unset ytics}

if(i==j){
plot chain every ::m u (bin(column(i+1))):(1.0) lc -1 smooth freq with boxes noti
}

if(i!=j){
plot chain every ::m u (column(i+1)):(column(j+1)) w p lc -1 noti
}

}

}

unset multiplot
