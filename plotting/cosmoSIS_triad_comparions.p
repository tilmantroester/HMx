reset

hmx(f1,f2)=sprintf('data/triad_Cl_%s-%s.dat',f1,f2)

cs(f1,f2,thing)=sprintf('/Users/Mead/Physics/people/Tilman/cosmoSIS/%s_%s_cl/%s.txt',f1,f2,thing)

#fields="'gal_z0.1-0.9' 'gal_z0.1-0.5' 'gal_z0.5-0.9' 'CMB' 'y'"
#fields="''"

set log x

set log y

f=10**5

plot hmx('gal_z0.1-0.9','CMB') u 1:3 w l lw 3 noti 'HMx: {/Symbol k}-{}',\
     '<paste '.cs('shear','cmbkappa','ell').' '.cs('shear','cmbkappa','bin_1_1').'' u 1:2 w l lw 3 noti 'cosmoSIS: {/Symbol k}-y',\
     '<paste '.cs('shear','cmbkappa','ell').' '.cs('shear','cmbkappa','bin_2_1').'' u 1:2 w l lw 3 noti 'cosmoSIS: {/Symbol k}-y',\
     '<paste '.cs('shear','cmbkappa','ell').' '.cs('shear','cmbkappa','bin_3_1').'' u 1:2 w l lw 3 noti 'cosmoSIS: {/Symbol k}-y',\
     '<paste '.cs('shear','cmbkappa','ell').' '.cs('shear','cmbkappa','bin_4_1').'' u 1:2 w l lw 3 noti 'cosmoSIS: {/Symbol k}-y'
