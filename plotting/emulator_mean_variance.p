reset
unset multiplot

# Set terminal
if(!exists('print')) {print=0}
if(print == 0) {set term qt dashed}
if(print == 1) {set term post enh col}

# Initial white space
print('')

# Set mode
if(!exists('imode')) {imode = 1}
print('imode = 1: Compare HMcode versions')
print('imode = 2: Compare two-halo terms')
print('imode = 3: Compare mass functions')
print('imode = 4: Compare concentration relations')
print('imode = 5: Compare halo overdensity definitions')
print('imode = 6: Compare delta_c and Delta_v')
print('imode = 7: Neutrino choices')
print('imode = ', imode)
print('')

# Function for file path
file(emu, ihm, zlab) = sprintf('data/%s_mean_variance_hm%d_z%s.dat', emu, ihm, zlab)

# Redshift label
ztitle(z) = sprintf('z = %1.1f', z)
labx = 0.05
laby = 0.90

# Set redshift stuff
nz = 4
array z[nz]
z[1] = 0.0
z[2] = 0.5
z[3] = 1.0
z[4] = 2.0
array zlab[nz]
zlab[1] = '0p0'
zlab[2] = '0p5'
zlab[3] = '1p0'
zlab[4] = '2p0'

# Set halomodels to compare
if(imode == 1) {nhm = 3; dr = 0.10; ds = 0.05}
if(imode == 2) {nhm = 3; dr = 0.30; ds = 0.20}
if(imode == 3) {nhm = 4; dr = 0.50; ds = 0.20}
if(imode == 4) {nhm = 6; dr = 0.30; ds = 0.20}
if(imode == 5) {nhm = 4; dr = 0.50; ds = 0.20}
if(imode == 6) {nhm = 4; dr = 0.50; ds = 0.20}
if(imode == 7) {nhm = 6; dr = 0.50; ds = 0.20}
array ihms[nhm]
array hmlabel[nhm]
if (imode == 1){
   emu = 'FrankenEmu'
   ihms[1] = 7
   ihms[2] = 1
   ihms[3] = 15
   hmlabel[1] = 'HMcode (2015)'
   hmlabel[2] = 'HMcode (2016)'
   hmlabel[3] = 'HMcode (2020)'
   outfile = 'plots/emulator_mean_variance_HMcode.eps'
}
if (imode == 2){
   emu = 'FrankenEmu'
   ihms[1] = 3
   ihms[2] = 72
   ihms[3] = 73
   hmlabel[1] = 'Standard two-halo term'
   hmlabel[2] = 'Linear two-halo term'
   hmlabel[3] = 'Linear two-halo term with damped wiggles'
   outfile = 'plots/emulator_mean_variance_twohalo.eps'
}
if (imode == 3){
   emu = 'FrankenEmu'
   ihms[1] = 3
   ihms[2] = 27
   ihms[3] = 23
   #ihms[4] = 80
   ihms[4] = 87
   hmlabel[1] = 'Sheth & Tormen: virial'
   hmlabel[2] = 'Press & Schecter: virial'
   hmlabel[3] = 'Tinker mass function: virial'
   #hmlabel[4] = 'Jenkins: m178'
   hmlabel[4] = 'Despali: virial'
   outfile = 'plots/emulator_mean_variance_massfunctions.eps'
}
if (imode == 4){
   emu = 'FrankenEmu'
   ihms[1] = 3
   ihms[2] = 68
   ihms[3] = 69
   ihms[4] = 70
   ihms[5] = 71
   ihms[6] = 90
   hmlabel[1] = 'Duffy et al. (2008)'
   hmlabel[2] = 'Bullock et al. (2001)'
   hmlabel[3] = 'Simple Bullock et al. (2001)'
   hmlabel[4] = 'Duffy et al. with no Dolag correction'
   hmlabel[5] = 'Duffy et al. with Dolag with 1.5 exponent'
   hmlabel[6] = 'Duffy et al. with z-dependent Dolag correction'
   outfile = 'plots/emulator_mean_variance_concentration.eps'
}
if (imode == 5){
   emu = 'FrankenEmu'
   hmlabel[1] = 'Sheth & Tormen: virial'
   hmlabel[2] = 'Tinker: virial'
   hmlabel[3] = 'Tinker: m200'
   hmlabel[4] = 'Tinker: m200c'
   ihms[1] = 3
   ihms[2] = 23
   ihms[3] = 44
   ihms[4] = 42
   outfile = 'plots/emulator_mean_variance_halodefinition.eps'
}
if (imode == 6){
   emu = 'FrankenEmu'
   hmlabel[1] = 'Sheth & Tormen: {/Symbol d}_c: Nakamura & Suto; {/Symbol D}_v: Bryan & Norman'
   hmlabel[2] = 'Sheth & Tormen: {/Symbol d}_c = 1.686; {/Symbol D}_v: Bryan & Norman'
   hmlabel[3] = 'Sheth & Tormen: {/Symbol d}_c: Mead; {/Symbol D}_v: Bryan & Norman'
   hmlabel[4] = 'Sheth & Tormen: {/Symbol d}_c: Nakamura & Suto; {/Symbol D}_v: Mead'
   ihms[1] = 3
   ihms[2] = 74
   ihms[3] = 75
   ihms[4] = 76
   outfile = 'plots/emulator_mean_variance_dcDv.eps'
}
if (imode == 7){
   emu = 'MiraTitan_neutrinos'
   hmlabel[1] = 'Sigma from un-normalised cold'
   hmlabel[2] = 'Sigma from total matter'
   hmlabel[3] = 'No correction for halo mass for neutrinos'
   hmlabel[4] = 'Sigma from normalised cold'
   hmlabel[5] = 'Neutrinos affect halo virial radius'
   hmlabel[6] = 'Mead (2017) spherical collapse calculation'
   ihms[1] = 3
   ihms[2] = 97
   ihms[3] = 98
   ihms[4] = 99
   ihms[5] = 100
   ihms[6] = 52
   outfile = 'plots/emulator_mean_variance_neutrinos.eps'
}

if (print == 1) {set output outfile}

# x axis stuff
rmin = 1.-dr
rmax = 1.+dr
rlab = 'Mean deviation'
klab = 'k / h Mpc^{-1}'
set log x
set xlabel klab

# y axis stuff
smin = 0.0
smax = ds
slab = 'Standard deviation'

set multiplot layout nz, 2 margins 0.09, 0.99, 0.07, 0.99 spacing 0.05, 0.02

   do for [iz = 1:nz] {

      if(iz == nz) {set xlabel klab; set format x} else {unset xlabel; set format x ''}

      set yrange [rmin:rmax]
      set ylabel rlab

      if(iz == 1) {set key top left} else {unset key}
      set label ztitle(z[iz]) at graph labx, laby

      plot 1 w l lt -1 noti,\
         for [ihm = 1:nhm] file(emu, ihms[ihm], zlab[iz]) u 1:2 w l lw 3 lc ihm dt 1 noti  

      unset label

      set yrange [smin:smax]
      set ylabel slab

      plot for [ihm = 1:nhm] NaN lc ihm lw 3 ti ''.hmlabel[ihm].'',\
         for [ihm = 1:nhm] file(emu, ihms[ihm], zlab[iz]) u 1:3 w l lw 3 lc ihm dt 1 noti

   }

unset multiplot

show output