reset
unset multiplot

# Set terminal
if(!exists('print')) {print=0}
if(print == 0) {set term qt dashed}
if(print == 1) {set term post enh col font ',8'; set output 'paper/emulator_mean_variance_all.eps'}

# Initial white space
print('')

# Function for file path
file(ihm, zlab) = sprintf('data/emulator_mean_variance_hm%d_z%s.dat', ihm, zlab)

# Redshift label
ztitle(z) = sprintf('z = %1.1f', z)
labx = 0.05
laby = 0.90

# Set redshift
z = 0.5
if(z == 0.0) {zlab = '0p0'}
if(z == 0.5) {zlab = '0p5'}
if(z == 1.0) {zlab = '1p0'}
if(z == 2.0) {zlab = '2p0'}

# Axis ranges
rmin = 0.75
rmax = 1.15
ds = 0.12

# Total number of different comparison plot pairs
nmode = 5

set multiplot layout nmode, 2 margins 0.09, 0.99, 0.07, 0.99 spacing 0.10, 0.02

do for [imode = 1:nmode] {

   # Set halomodels to compare
   #if(imode == 1) {nhm = 3; dr = 0.10; ds = 0.05}
   if(imode == 1) {nhm = 3}#; dr = 0.30; ds = 0.20}
   if(imode == 2) {nhm = 5}#; dr = 0.50; ds = 0.20} 
   if(imode == 3) {nhm = 4}#; dr = 0.50; ds = 0.20}
   if(imode == 4) {nhm = 4}#; dr = 0.50; ds = 0.20}
   if(imode == 5) {nhm = 5}#; dr = 0.30; ds = 0.20}
   array ihms[nhm]
   array hmlabel[nhm]
   #if (imode == 1){
   #   ihms[1] = 7
   #   ihms[2] = 1
   #   ihms[3] = 15
   #   hmlabel[1] = 'HMcode (2015)'
   #   hmlabel[2] = 'HMcode (2016)'
   #   hmlabel[3] = 'HMcode (2020)'
   #   outfile = 'plots/emulator_mean_variance_HMcode.eps'
   #}
   if (imode == 1){
      # Two-halo terms
      ihms[1] = 3
      ihms[2] = 72
      ihms[3] = 73
      hmlabel[1] = 'Standard two-halo term'
      hmlabel[2] = 'Linear two-halo term'
      hmlabel[3] = 'Linear two-halo term with damped wiggles'
   }
   if (imode == 2){
      # Mass functions with virial definition
      ihms[1] = 3
      ihms[2] = 27
      ihms[3] = 91
      ihms[4] = 23
      ihms[5] = 87
      hmlabel[1] = 'Sheth (1999): virial'
      hmlabel[2] = 'Press (1974): virial'
      hmlabel[3] = 'Tinker (2008): virial'
      hmlabel[4] = 'Tinker (2010): virial'
      hmlabel[5] = 'Despali (2016): virial'
   }
   if (imode == 3){
      # Mass functions with non-virial definitions
      hmlabel[1] = 'Sheth (1999): virial'
      hmlabel[2] = 'Tinker (2010): virial'
      hmlabel[3] = 'Tinker (2010): M_{200}'
      hmlabel[4] = 'Tinker (2010): M_{200,c}'
      ihms[1] = 3
      ihms[2] = 23
      ihms[3] = 44
      ihms[4] = 42
   }
   if (imode == 4){
      # Halo definitions
      hmlabel[1] = '{/Symbol d}_c: Nakamura (1997); {/Symbol D}_v: Bryan (1998)'
      hmlabel[2] = '{/Symbol d}_c = 1.686; {/Symbol D}_v: Bryan (1998)'
      hmlabel[3] = '{/Symbol d}_c: Mead (2017); {/Symbol D}_v: Bryan (1998)'
      hmlabel[4] = '{/Symbol d}_c: Nakamura (1997); {/Symbol D}_v: Mead (2017)'
      ihms[1] = 3
      ihms[2] = 74
      ihms[3] = 75
      ihms[4] = 76
   }
   if (imode == 5){
      # Concentration-mass relations
      ihms[1] = 3
      ihms[2] = 68
      ihms[3] = 69
      ihms[4] = 70
      #ihms[5] = 71
      #ihms[5] = 88
      ihms[5] = 90
      hmlabel[1] = 'Duffy (2008)'
      hmlabel[2] = 'Bullock (2001)'
      hmlabel[3] = 'Simple Bullock (2001)'
      hmlabel[4] = 'Duffy (2008) with no Dolag (2004) correction'
      #hmlabel[5] = 'Duffy et al. with Dolag with 1.5 exponent'
      #hmlabel[5] = 'Child (2018)'
      hmlabel[5] = 'Duffy (2008) with z-dependent Dolag (2004) correction'
   }

   # x axis stuff
   rlab = 'Mean deviation'
   klab = 'k / h Mpc^{-1}'
   set log x

   # y axis stuff
   smin = 0.0
   smax = ds
   slab = 'Standard deviation'

   if(imode == nmode) {set xlabel klab; set format x} else {unset xlabel; set format x ''}

   set yrange [rmin:rmax]
   set ylabel rlab

   set key top left
   if(imode == 1) {set label ztitle(z) at graph labx, laby}

   plot 1 w l lt -1 noti,\
      for [ihm = 1:nhm] file(ihms[nhm-ihm+1], zlab) u 1:2 w l lw 3 lc nhm-ihm+1 dt 1 noti  

   unset label

   set yrange [smin:smax]
   set ylabel slab

   plot for [ihm = 1:nhm] NaN lc nhm-ihm+1 lw 3 ti ''.hmlabel[nhm-ihm+1].'',\
      for [ihm = 1:nhm] file(ihms[nhm-ihm+1], zlab) u 1:3 w l lw 3 lc nhm-ihm+1 dt 1 noti

}

unset multiplot

show output