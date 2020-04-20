reset
unset multiplot

# Set terminal
set term qt dashed

# Initial white space
print('')

# Set mode
if(!exists('imode')) {imode = 1}
print('imode = 1: Compare HMcode versions')
print('imode = 2: Compare two-halo terms')
print('imode = 3: Compare mass functions')
print('imode = 4: Compare concentration relations')
print('imode = ', imode)
prin('')

# Function for file path
file(ihm, zlab) = sprintf('data/emulator_mean_variance_hm%d_z%s.dat', ihm, zlab)

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
if(imode == 3) {nhm = 9; dr = 0.50; ds = 0.20}
if(imode == 4) {nhm = 5; dr = 0.30; ds = 0.20}
array ihms[nhm]
array hmlabel[nhm]
if (imode == 1){
   ihms[1] = 7
   ihms[2] = 1
   ihms[3] = 15
   hmlabel[1] = 'HMcode (2015)'
   hmlabel[2] = 'HMcode (2016)'
   hmlabel[3] = 'HMcode (2020)'
}
if (imode == 2){
   ihms[1] = 3
   ihms[2] = 72
   ihms[3] = 73
   hmlabel[1] = 'Standard two-halo term'
   hmlabel[2] = 'Linear two-halo term'
   hmlabel[3] = 'Linear two-halo term with damped wiggles'
}
if (imode == 3){
   ihms[1] = 3
   ihms[2] = 27
   ihms[3] = 23
   ihms[4] = 44
   ihms[5] = 42
   ihms[6] = 74
   ihms[7] = 75
   ihms[8] = 76
   ihms[9] = 80
   hmlabel[1] = 'Sheth & Tormen: virial'
   hmlabel[2] = 'Press & Schecter: virial'
   hmlabel[3] = 'Tinker mass function: virial'
   hmlabel[4] = 'Tinker mass function: m200'
   hmlabel[5] = 'Tinker mass function: m200c'
   hmlabel[6] = 'Sheth & Tormen: {/Symbol d}_c = 1.686 fixed'
   hmlabel[7] = 'Sheth & Tormen: {/Symbol d}_c from Mead fit'
   hmlabel[8] = 'Sheth & Tormen: {/Symbol D}_v from Mead fit'
   hmlabel[9] = 'Jenkins: m178'
   #hmlabel[9] = 'Sheth & Tormen: {/Symbol d}_c and {/Symbol D}_v from Mead fit'
}
if (imode == 4){
   ihms[1] = 3
   ihms[2] = 68
   ihms[3] = 69
   ihms[4] = 70
   ihms[5] = 71
   hmlabel[1] = 'Duffy'
   hmlabel[2] = 'Bullock'
   hmlabel[3] = 'Simple Bullock'
   hmlabel[4] = 'No Dolag correction'
   hmlabel[5] = 'Dolag with 1.5 exponent'
}

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

      if(iz==1) {set key top left} else {unset key}
      set label ztitle(z[iz]) at graph labx, laby

      plot 1 w l lt -1 noti,\
         for [ihm = 1:nhm] file(ihms[ihm], zlab[iz]) u 1:2 w l lw 3 lc ihm dt 1 noti  

      unset label

      set yrange [smin:smax]
      set ylabel slab

      plot for [ihm = 1:nhm] NaN lc ihm lw 3 ti 'Halo model: '.hmlabel[ihm].'',\
         for [ihm = 1:nhm] file(ihms[ihm], zlab[iz]) u 1:3 w l lw 3 lc ihm dt 1 noti

   }

unset multiplot