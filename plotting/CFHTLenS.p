reset

# Initial white space
print ''

# Terminal
if(!exists('print')){print=0}
if(print==0) set term qt dashed
if(print==1) set term post enh col; set output 'CFHTLenS.eps'

ell='l'

# File names
types="'linear' '2h' '1h' 'hm'"
xi(type)=sprintf('data/CFHTLenS_xi_%s.dat',type)
cl(type)=sprintf('data/CFHTLenS_cl_%s.dat',type)

# CFHTLenS xi data
data='/Users/Mead/Physics/data/CFHTLenS/corr.txt'

# Colours
xip_col='red'
xim_col='orange'
cl_col='blue'

# Labels
xip_lab='{/Symbol x}_+'
xim_lab='{/Symbol x}_-'

#NaN=1e-18
# Arcminute to degree conversion
arcmin=60.

# Offset for xi plus and minus
offset=0.975

# Theta range
th_min=0.01
th_max=10.

# Xi range
xi_min=3e-8
xi_max=3e-4

print 'Select plot type'
if(!exists('iplot')){iplot=1}
print 'iplot = 1: Correlation function'
print 'iplot = 2: Power spectrum'
print 'iplot = ', iplot
print ''

# xi
if(iplot==1){

set log x
set xlabel '{/Symbol q} [degrees]'
set xrange [th_min:th_max]
set xtics nomirror

#Set x2 axis for arcminutes
set log x2
set format x2
set x2range [th_min*arcmin:th_max*arcmin]
set x2label '{/Symbol q} [arcminutes]'
set x2tics

set log y
set ylabel '{/Symbol x}({/Symbol q})'
set format y '10^{%T}'
set yrange [xi_min:xi_max]
set mytics 10

plot NaN w l lw 3 lc rgb xip_col dt 1 ti xip_lab,\
     NaN w l lw 3 lc rgb xim_col dt 1 ti xim_lab,\
     NaN w l lw 3 lc -1 dt 1 ti 'Full',\
     NaN w l lw 3 lc -1 dt 2 ti 'One halo',\
     NaN w l lw 3 lc -1 dt 3 ti 'Two halo',\
     NaN w l lw 3 lc -1 dt 4 ti 'Linear',\
     for [i=1:words(types)] xi(word(types,i)) u 1:2 w l lw 3 lc rgb xip_col dt 5-i noti,\
     for [i=1:words(types)] xi(word(types,i)) u 1:4 w l lw 3 lc rgb xim_col dt 5-i noti,\
     data u ($1/arcmin/offset):2 w p ps 1. pt 7 lc rgb 'black' ti 'CFHTLenS {/Symbol x}_+',\
     data u ($1/arcmin/offset):2:(sqrt($3)) w errorbars ps 0 lt 1 lc rgb 'black' noti,\
     data u (offset*$1/arcmin):4 w p ps 1. pt 7 lc rgb 'dark-grey' ti 'CFHTLenS {/Symbol x}_-',\
     data u (offset*$1/arcmin):4:(sqrt($5)) w errorbars ps 0 lt 1 lc rgb 'dark-grey' noti

}

# Cl
if(iplot==2){

   cl_min=1e-14
   cl_max=1e-7

   set log x
   set xlabel ''.ell.''
   set format x '10^{%T}'
   
   set log y
   set ylabel 'C('.ell.')'
   set yrange [cl_min:cl_max]

   plot \
      NaN w l lw 3 lc -1 dt 1 ti 'Full',\
      NaN w l lw 3 lc -1 dt 2 ti 'One halo',\
      NaN w l lw 3 lc -1 dt 3 ti 'Two halo',\
      NaN w l lw 3 lc -1 dt 4 ti 'Linear',\
      for [i=1:words(types)] cl(word(types,i)) u 1:2 w l lw 3 lc rgb cl_col dt 5-i noti

}

show output
