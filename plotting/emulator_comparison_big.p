unset multiplot
reset

# Terminal options
if(!exists('print')){print=0}
if(print==0) {set term qt dashed}
if(print==1) {set term post enh col font ',8' size 8,8}

# Initial white space
print ''

# Choose the mode
if(!exists('imode')) {imode = 2}
print('imode = 1: Cosmic Emu (all)')
print('imode = 2: Franken Emu (all)')
print('imode = 3: Mira Titan (all)')
print('imode = 4: Mira Titan (all; colour by neutrino mass)')
print('imode = 5: Cosmic Emu')
print('imode = 6: PAPER: Franken Emu')
print('imode = 7: Mira Titan')
print('imode = 8: Mira Titan (neutrino nodes; colour by neutrino mass)')
print('imode = 9: PAPER: Mira Titan (all nodes; colour by neutrino mass)')
print('imode = ', imode)
print('')

# Set the emulator
if(imode==1 || imode==5) {emulator = 'CosmicEmu';  icos1=1;  icos2=37}
if(imode==2 || imode==6) {emulator = 'FrankenEmu'; icos1=1;  icos2=37}
if(imode==3 || imode==7) {emulator = 'MiraTitan';  icos1=1;  icos2=36}
if(imode==4 || imode==8) {emulator = 'MiraTitan';  icos1=11; icos2=36}
if(imode==9) {emulator = 'MiraTitan';  icos1=1; icos2=36}

# Set file if writing to disk
if(print == 1){
   if(imode == 1) {set output 'plots/CosmicEmu_comparison.eps'}
   if(imode == 2) {set output 'plots/FrankenEmu_comparison.eps'}
   if(imode == 3) {set output 'plots/MiraTitan_comparison.eps'}
   if(imode == 4) {set output 'plots/MiraTitan_comparison_neutrinos.eps'}
   if(imode == 5) {set output 'plots/CosmicEmu_comparison.eps'}
   if(imode == 6) {set output 'paper/FrankenEmu_comparison.eps'}
   if(imode == 7) {set output 'plots/MiraTitan_comparison.eps'}
   if(imode == 8) {set output 'plots/MiraTitan_comparison_neutrinos.eps'}
   if(imode == 9) {set output 'paper/MiraTitan_comparison_neutrinos.eps'}
}

# Input file name
file(name, emulator, icos, zlab) = sprintf('data/power_%s_%s_cos%d_z%s.dat', name, emulator, icos, zlab)

# Function for redshift label
zlabel(z) = sprintf('z = %1.1f', z)
zlab_x = 0.75
zlab_y = 0.88

# Comparison label
clab_x = 0.05
clab_y = zlab_y

# k range
kmin = 2e-3
kmax = 12e0
klab = 'k / h Mpc^{-1}'
set log x
set xrange [kmin:kmax]

# Number of redshifts
if(!exists('nz')) {nz = 4}
print 'Number of redshifts: nz: ', nz
print ''

# Number of cosmologies
print 'First cosmology to plot: icos1: ', icos1
print 'Last cosmology to plot: icos2: ', icos2
print 'Total number of cosmologies to plot: ', icos2-icos1+1
print ''

# Residual range
dy = 0.17
ddy = 0.05
set yrange [1.-dy:1+dy]

#numax = 
set cblabel 'm_{/Symbol n} / eV'
#set cbrange [0.:numax]

# Margins
x1 = 0.09
if (imode == 1 || imode == 2 || imode == 3 || imode == 5 || imode == 6 || imode == 7) {x2 = 0.99}
if (imode == 4 || imode == 8 || imode == 9) {x2 = 0.91}
y1 = 0.07
y2 = 0.99
cbsep = 0.01
cbwid = 0.015

# Number of comparisons
if(imode == 1 || imode == 2 || imode == 3 || imode == 4) {ncom = 9}
if(imode == 5 || imode == 6 || imode == 7 || imode == 8 || imode == 9) {ncom = 6}

# Setup multiplot
set multiplot layout ncom, nz margins x1, x2, y1, y2 spacing 0.0, 0.0

# Colourbox
#set palette defined ( 0 'royalblue', 1 'orchid' )
set palette defined (1 'royalblue', 2 'grey', 3 'light-red')
set colorbox vertical user origin x2+cbsep, y1 size cbwid, y2-y1

# Loop over the comparisons (rows)
do for [icom = 1:ncom]{

   # Columns for comparison 
   if(imode == 1 || imode == 2 || imode == 3 || imode == 4) {
      if(icom == 1) {comp = 'HMcode_2020';             comp_title = 'HMcode (2020)'}
      if(icom == 2) {comp = 'HMcode_2019';             comp_title = 'HMcode (2019)'}
      if(icom == 3) {comp = 'HMcode_2018';             comp_title = 'HMcode (2018)'}
      if(icom == 4) {comp = 'HMcode_2016_neutrinofix'; comp_title = 'HMcode (2016) with neutrino fix'}
      if(icom == 5) {comp = 'HMcode_2016';             comp_title = 'Mead et al. (2016)'}
      if(icom == 6) {comp = 'HMcode_2015';             comp_title = 'Mead et al. (2015)'}
      if(icom == 7) {comp = 'HALOFIT_Takahashi';       comp_title = 'Takahashi et al. (2012)'}
      if(icom == 8) {comp = 'HALOFIT_Smith';           comp_title = 'Smith et al. (2003)'}
      if(icom == 9) {comp = 'linear';                  comp_title = 'Linear theory'}
   }
   if(imode == 5 || imode == 6 || imode == 7 || imode == 8 || imode == 9) {
      if(icom == 1) {comp = 'HMcode_2020';       comp_title = 'This paper'}
      if(icom == 2) {comp = 'HMcode_2016';       comp_title = 'Mead et al. (2016)'}
      if(icom == 3) {comp = 'HMcode_2015';       comp_title = 'Mead et al. (2015)'}
      if(icom == 4) {comp = 'HALOFIT_Takahashi'; comp_title = 'Takahashi et al. (2012)'}
      if(icom == 5) {comp = 'HALOFIT_Smith';     comp_title = 'Smith et al. (2003)'}
      if(icom == 6) {comp = 'linear';            comp_title = 'Linear theory'}
   }

   # Loop over redshifts (columns)
   do for [iz = 1:nz]{

      if(iz == 1) {z = 0.0; zfile = '0p0'}
      if(iz == 2) {z = 0.5; zfile = '0p5'}
      if(iz == 3) {z = 1.0; zfile = '1p0'}
      if(iz == 4) {z = 2.0; zfile = '2p0'}
      if(iz == 5) {z = 3.0; zfile = '3p0'}
      if(iz == 6) {z = 4.0; zfile = '4p0'}

      # x axis options
      set xlabel ''; set format x ''
      if(icom == ncom) {set xlabel klab; set format x}

      # y axis options
      set format y ''
      set ylabel ''
      if(iz == 1) {set format y; set ylabel 'P(k) / P_{emu}(k)'}

      # Plot titles
      unset title
      if(iz == 1) {set label comp_title at graph clab_x, clab_y}
   
      # Redshift label
      if(icom == 1) {set label zlabel(z) at graph zlab_x, zlab_y}

      # Plot
      if (imode == 1 || imode == 2 || imode == 3 || imode == 5 || imode == 6 || imode == 7){
         plot 1 w l lt -1 noti,\
            1.+ddy w l lc -1 dt 2 noti,\
            1.-ddy w l lc -1 dt 2 noti,\
            for [icos = icos1:icos2] file(comp, emulator, icos, zfile) u 1:(column(2)/column(3)) w l lc icos noti
      }
      if(imode == 4 || imode == 8 || imode == 9){
         plot 1 w l lt -1 noti,\
            1.+ddy w l lc -1 dt 2 noti,\
            1.-ddy w l lc -1 dt 2 noti,\
            for [icos = icos1:icos2] file(comp, emulator, icos, zfile) u 1:(column(2)/column(3)):4 w l lc palette noti
         unset colorbox
      }

      # Remove label
      unset label

   }

}

unset multiplot

show output