unset multiplot
reset session

bias='data/bias_functions.dat'
mass='data/mass_functions.dat'

set xlabel '{/Symbol n}'

set log y

set multiplot layout 2,1

do for [i=1:2]{

if(i==1){file=mass; set ylabel 'g({/Symbol n})'}
if(i==2){file=bias; set key top left; set ylabel 'b({/Symbol n})'}

plot file u 1:2 w l lw 3 ti 'Press & Schecter (1974)',\
     file u 1:3 w l lw 3 ti 'Sheth & Tormen (1999)',\
     file u 1:4 w l lw 3 ti 'Tinker et al. (2010)'

}

unset multiplot


