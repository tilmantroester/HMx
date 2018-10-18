reset

chain='fitting/mcmc.dat'

np=120

plot for [i=1:np] chain u (column(0)):(column(i+1)) w p noti
