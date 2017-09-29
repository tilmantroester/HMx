reset
unset multiplot

cmsy='/Users/Mead/Fonts/cmsy10.pfb'

if(print==0) {set term aqua dashed}
if(print==1) {set term post enh col fontfile cmsy; set output 'halo_pressure.eps'}

msun='{/cmsy10 \014}'

#Use a log r axis or not
ilog=0

#So that units include the Hubble parameter factors
ih=1

#Do comparison or not
icomp=1

if(ih==0){
file1='diagnostics/halo_profile_m14noh.dat'
file2='diagnostics/halo_profile_m15noh.dat'
}
if(ih==1){
file1='diagnostics/halo_profile_m13.dat'
file2='diagnostics/halo_profile_m14.dat'
file3='diagnostics/halo_profile_m15.dat'
file1a='diagnostics/pressure_comparison/halo_profile_m13_gas.dat'
file1b='diagnostics/pressure_comparison/halo_profile_m13_UPP.dat'
file2a='diagnostics/pressure_comparison/halo_profile_m14_gas.dat'
file2b='diagnostics/pressure_comparison/halo_profile_m14_UPP.dat'
file3a='diagnostics/pressure_comparison/halo_profile_m15_gas.dat'
file3b='diagnostics/pressure_comparison/halo_profile_m15_UPP.dat'
}

#Constants
JpeV=1.602e-19 # 1 = 1.602e-19 J/eV
mpcm=100. #1 = 100 cm/m

if(ilog==0) {set xrange [0:8]}
if(ilog==1) {set log x; set xrange [1e-3:1e1]}
if(ih==0)   {set xlabel 'r / Mpc'}
if(ih==1)   {set xlabel 'r / (h^{-1} Mpc)'}

if(ilog==0) {set yrange [*:*]}#{set yrange [0:3.5]}
if(ilog==1) {set log y; set yrange [1e-4:1e1]; set format y '10^{%T}'; set mytics 10}
set ylabel '(r / Mpc)^2 P_e(r,M) / (eV cm^{-3})'

set multiplot layout 1,2

#print 'Wanker'

do for [i=1:2] {

#print 'Extreme wanker'

if(ih==0 && i==1){tits='M = 10^{14} M_'.msun.''}
if(ih==0 && i==2){tits='M = 10^{15} M_'.msun.''}
if(ih==1 && i==1){tits='M = 10^{14} h^{-1} M_'.msun.'; z = 0'}
if(ih==1 && i==2){tits='M = 10^{15} h^{-1} M_'.msun.'; z = 0'}
if(i==1 && icomp==0){file=file2}
if(i==2 && icomp==0){file=file3}
if(i==1 && icomp==1){filea=file2a; fileb=file2b}
if(i==2 && icomp==1){filea=file3a; fileb=file3b}

set title tits

if(icomp==0){
plot file u ($1):(($1)*($1)*$7/(JpeV*mpcm**3)) w l lw 3 dt 1 lc rgb 'black' noti
}

if(icomp==1){
plot filea u ($1):(($1)*($1)*$7/(JpeV*mpcm**3)) w l lw 3 dt 1 lc rgb 'black' ti 'gas',\
     fileb u ($1):(($1)*($1)*$7/(JpeV*mpcm**3)) w l lw 3 dt 1 lc rgb 'red'   ti 'UPP'
}

}

unset multiplot
