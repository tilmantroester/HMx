reset

fiducial='cosmo-OWLS/data/DMONLY.dat'
power(i,j,mod)=sprintf('cosmo-OWLS/data/power_%s_%d%d.dat',mod,i,j)
sim(m,mod)=sprintf('/Users/Mead/Physics/cosmo-OWLS/power/N%d/%s_all_all_power.dat',m,mod)

mods="'REF' 'NOCOOL' 'AGN' 'AGN8p5' 'AGN8p7'"
owls="'REF' 'NOCOOL_UVB' 'AGN' 'AGN_Theat_8p5' 'AGN_Theat_8p7'"

set log x
set xlabel 'k / h Mpc^{-1}'

set ylabel 'P(k) / P_{DMONLY}(k)'
set yrange [0.6:1.4]

set key top left

c1=2
L1=4
c2=5
L2=5
m=800

plot 1 w l lt -1 noti,\
     for [i=1:words(owls)] '<paste '.sim(m,word(owls,i)).' '.sim(m,'DMONLY').'' u 1:(column(c1)/column(c1+L1)) w p pt 7 dt 1 lc i noti,\
     for [i=1:words(mods)] '<paste '.power(0,0,word(mods,i)).' '.fiducial.''    u 1:(column(c2)/column(c2+L2)) w l lw 3 lc i ti word(mods,i)
