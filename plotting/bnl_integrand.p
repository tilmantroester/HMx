reset

print ''

if(!exists('print')) {print=0}
print 'print = 0: aqua'
print 'print = 1: eps'
if(print==0) {set term aqua font ',8' size 1000,600}
if(print==1) {set term post enh col font ',8' size 10,6; set output 'bnl.eps'}
print 'print = ', print
print ''

set size square
set view map

infile(k)=sprintf('data/bnl_k%1.1f.dat',k)

if(!exists('imode')) {imode=1}
print 'imode = 1 - Plot B_NL'
print 'imode = 2 - Plot integrand'
print 'imode = 3 - Plot both'
print 'imode = ', imode
print ''

if(imode==1 || imode==2){
if(!exists('iplot')) {iplot=1}
print 'iplot = 1 - Plot a single k value'
print 'iplot = 2 - Plot 9 k values from 0.1 to 0.9'
print 'iplot = ', iplot
print ''
}

lab_bnl='B^{NL}({/Symbol n}_1,{/Symbol n}_2,k,z)'
lab_Inl='I^{NL}_{mm}(k,z) integrand'
#lab_Inl='g({/Symbol n}_1)g({/Symbol n}_2)b({/Symbol n}_1)b({/Symbol n}_2)B_{NL}({/Symbol n}_1,{/Symbol n}_2,k)'

if(imode==1) {c=3; lab=lab_bnl}
if(imode==2) {c=4; lab=lab_Inl}

# nu range on x and y
numin=1e-1
numax=3.
set xrange[numin:numax]
set yrange[numin:numax]
set xlabel '{/Symbol n}_1'
set ylabel '{/Symbol n}_2'

set cblabel lab
#set cbrange [cbmin:cbmax]
ds=0.2

# Title for plots
tits(k)=sprintf('k = %1.1f h/Mpc', k)

if(imode==1 || imode==2){

if(iplot==1) {nk=1}
if(iplot==2) {set multiplot layout 3,3 margins 0.02,0.98,0.1,0.98 spacing 0.1,0.1; nk=9}

do for [i=1:nk] {

if(iplot==1){
if(!exists('k')) {k=0.3}
print 'k [h/Mpc]: ', k
print ''
}

if(iplot==2){
if(i==1) {k=0.1}
if(i==2) {k=0.2}
if(i==3) {k=0.3}
if(i==4) {k=0.4}
if(i==5) {k=0.5}
if(i==6) {k=0.6}
if(i==7) {k=0.7}
if(i==8) {k=0.8}
if(i==9) {k=0.9}
#set xlabel ''
#set ylabel ''
#set cblabel ''
}

set title tits(k)
splot infile(k) u 1:2:(column(c)) with image noti

}

if(iplot==2) {unset multiplot}

}

if(imode==3) {

nk=3
set multiplot layout 2,nk margins 0.05,0.93,0.06,0.98 spacing 0.08,0.01

do for [j=1:2] {

do for [i=1:nk] {

if(i==1) {k=0.1}
if(i==2) {k=0.3}
if(i==3) {k=1.0}

if(j==1) {c=3; set cblab lab_bnl; set xlabel ''; set format x ''}
if(j==2) {c=4; set cblab lab_Inl; set xlabel '{/Symbol n}_1'; set format x}

if(i==1 && j==1) {dcb=0.2}
if(i==2 && j==1) {dcb=0.6}
if(i==3 && j==1) {dcb=4.0}
if(i==1 && j==2) {dcb=0.015}
if(i==2 && j==2) {dcb=0.25}
if(i==3 && j==2) {dcb=0.9}
cbmin=-dcb
cbmax=dcb
#set cbrange [cbmin:cbmax]

if(j==1) {set title tits(k)} else {set title ''}
splot infile(k) u 1:2:(column(c)) with image noti

}
}

unset multiplot

}

show output


