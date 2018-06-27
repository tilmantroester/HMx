PROGRAM HMx_driver

  USE HMx
  USE Limber
  USE file_info
  USE random_numbers

  !Parameter definitions
  IMPLICIT NONE
  REAL, ALLOCATABLE :: k(:), a(:)
  REAL, ALLOCATABLE :: pow_lin(:), pow_2h(:), pow_1h(:), pow_full(:)
  REAL, ALLOCATABLE :: powa(:,:), powa_lin(:,:), powa_2h(:,:), powa_1h(:,:), powa_full(:,:)
  REAL, ALLOCATABLE :: ell(:), Cell(:), theta(:), xi(:,:)
  REAL, ALLOCATABLE :: z_tab(:)
  INTEGER :: i, j, nk, na, j1, j2, n, nl, nz, nth, nnz, m, ipa, npa, ii
  INTEGER :: ip(2), ix(2), ixx(2), ihalo
  REAL :: kmin, kmax, amin, amax, lmin, lmax, thmin, thmax, zmin, zmax
  REAL :: z, z1, z2, r1, r2, a1, a2
  TYPE(cosmology) :: cosm
  TYPE(halomod) :: hmod
  TYPE(projection) :: proj(2)
  TYPE(lensing) :: lens
  CHARACTER(len=256) :: infile, outfile, base, mid, ext, dir, name, fname
  CHARACTER(len=256) :: mode, halomodel, red, cosmo
  INTEGER :: imode, icosmo, iowl, ihm, irho, imeth
  REAL :: sig8min, sig8max
  INTEGER :: ncos
  REAL :: m1, m2, mass
  REAL :: c, rmin, rmax, rv, rs, p1, p2
  REAL :: t1, t2, w, w0
  LOGICAL :: timing
  
  !Baryon stuff
  REAL :: param_min, param_max, param
  LOGICAL :: ilog

  !Halo-model Parameters
  LOGICAL, PARAMETER :: verbose=.TRUE. !Verbosity
  REAL, PARAMETER :: mmin=1e7 !Minimum halo mass for the calculation
  REAL, PARAMETER :: mmax=1e17 !Maximum halo mass for the calculation

  !Output choices
  LOGICAL, PARAMETER :: icumulative=.TRUE. !Do cumlative distributions for breakdown
  LOGICAL, PARAMETER :: ixi=.TRUE. !Do correlation functions from C(l)
  LOGICAL, PARAMETER :: ifull=.FALSE. !Do only full halo model C(l), xi(theta) calculations (quicker, no breakdown ...)

  CALL get_command_argument(1,mode)
  IF(mode=='') THEN
     imode=-1
  ELSE
     READ(mode,*) imode
  END IF

  CALL get_command_argument(2,cosmo)
  IF(cosmo=='') THEN
     icosmo=-1
  ELSE
     READ(cosmo,*) icosmo
  END IF

  CALL get_command_argument(3,halomodel)
  IF(halomodel=='') THEN
     ihm=-1
  ELSE
     READ(halomodel,*) ihm
  END IF

  !Initial space
  WRITE(*,*)

  !Choose mode
  IF(imode==-1) THEN
     WRITE(*,*) 'HMx_DRIVER: Choose what to do'
     WRITE(*,*) '============================='
     WRITE(*,*) ' 0 - Matter power spectrum at z = 0'
     WRITE(*,*) ' 1 - Matter power spectrum over multiple z'
     WRITE(*,*) ' 2 - Produce all halo components cross and auto spectra'
     WRITE(*,*) ' 3 - Run diagnostics'
     WRITE(*,*) ' 4 - Do random cosmologies for bug testing'
     WRITE(*,*) ' 5 - NOT SUPPORTED: Pressure field comparison'
     WRITE(*,*) ' 6 - n(z) check'
     WRITE(*,*) ' 7 - Do cross correlation'
     WRITE(*,*) ' 8 - Cross correlation as a function of cosmology'
     WRITE(*,*) ' 9 - Breakdown correlations in halo mass'
     WRITE(*,*) '10 - Breakdown correlations in redshift'
     WRITE(*,*) '11 - Breakdown correlations in halo radius'
     WRITE(*,*) '12 - Project triad'
     WRITE(*,*) '13 - Cross-correlation coefficient'
     WRITE(*,*) '14 - 3D spectra as baryon parameters vary'
     WRITE(*,*) '15 - 3D spectra for cosmo-OWLS models'
     WRITE(*,*) '16 - 3D spectra for BAHAMAS models'
     WRITE(*,*) '17 - 3D spectra for user choice of fields'
     WRITE(*,*) '18 - 3D bias'
     WRITE(*,*) '19 - CCL comparison'
     WRITE(*,*) '20 - Make Ma et al. (2015) Fig. 1'
     WRITE(*,*) '21 - W(k) integrand diagnostics'
     WRITE(*,*) '22 - Time W(k) integration methods'
     WRITE(*,*) '23 - Produce results from Mead (2017)'
     READ(*,*) imode
     WRITE(*,*) '============================'
     WRITE(*,*)
  END IF

  IF(imode==0) THEN

     !Set number of k points and k range (log spaced)
     nk=128
     kmin=1e-3
     kmax=1e2
     CALL fill_array(log(kmin),log(kmax),k,nk)
     k=exp(k)
     ALLOCATE(pow_lin(nk),pow_2h(nk),pow_1h(nk),pow_full(nk))

     !Assigns the cosmological model
     CALL assign_cosmology(icosmo,cosm,verbose)

     !Normalises power spectrum (via sigma_8) and fills sigma(R) look-up tables
     CALL print_cosmology(cosm)

     !Sets the redshift
     z=0.

     !Initiliasation for the halomodel calcualtion
     CALL init_halomod(ihm,mmin,mmax,z,hmod,cosm,verbose)

     !Do the halo-model calculation
     CALL calculate_halomod(-1,-1,k,nk,z,pow_lin,pow_2h,pow_1h,pow_full,hmod,cosm,verbose)

     !Write out the results
     outfile='data/power.dat'
     CALL write_power(k,pow_lin,pow_2h,pow_1h,pow_full,nk,outfile,verbose)

     !Write the one-void term if necessary
     IF(hmod%voids) THEN
        OPEN(8,file='data/power_1void.dat')
        DO i=1,nk     
           WRITE(8,*) k(i), p_1v(k(i),hmod)
        END DO
        CLOSE(8)
     END IF

  ELSE IF(imode==19) THEN

     !Set number of k points and k range (log spaced)
     nk=128
     kmin=1e-3
     kmax=1e1
     CALL fill_array(log(kmin),log(kmax),k,nk)
     k=exp(k)
     ALLOCATE(pow_lin(nk),pow_2h(nk),pow_1h(nk),pow_full(nk))

     !Set the halo model to that for CCL tests
     ihm=9

     DO j=1,3

        !Assigns the cosmological model
        IF(j==1) icosmo=1
        IF(j==2) icosmo=2
        IF(j==3) icosmo=3
        CALL assign_cosmology(icosmo,cosm,verbose)

        !Normalises power spectrum (via sigma_8) and fills sigma(R) look-up tables
        CALL print_cosmology(cosm)

        DO i=1,2

           !Sets the redshift
           IF(i==1) THEN
              z=0.
              IF(j==1) outfile='/Users/Mead/Physics/HMx/benchmarks/HMx_power_model1_z0.txt'
              IF(j==2) outfile='/Users/Mead/Physics/HMx/benchmarks/HMx_power_model2_z0.txt'
              IF(j==3) outfile='/Users/Mead/Physics/HMx/benchmarks/HMx_power_model3_z0.txt'
           ELSE IF(i==2) THEN
              z=1.
              IF(j==1) outfile='/Users/Mead/Physics/HMx/benchmarks/HMx_power_model1_z1.txt'
              IF(j==2) outfile='/Users/Mead/Physics/HMx/benchmarks/HMx_power_model2_z1.txt'
              IF(j==3) outfile='/Users/Mead/Physics/HMx/benchmarks/HMx_power_model3_z1.txt'
           END IF

           !Initialise the halo-model calculation
           CALL init_halomod(ihm,mmin,mmax,z,hmod,cosm,verbose)

           !Do the halo-model calculation
           CALL calculate_halomod(-1,-1,k,nk,z,pow_lin,pow_2h,pow_1h,pow_full,hmod,cosm,verbose)

           !Write out the results
           CALL write_power(k,pow_lin,pow_2h,pow_1h,pow_full,nk,outfile,verbose)

        END DO

     END DO

  ELSE IF(imode==1) THEN

     !Assigns the cosmological model
     !icosmo=-1
     CALL assign_cosmology(icosmo,cosm,verbose)
     CALL print_cosmology(cosm)

     !Set number of k points and k range (log spaced)
     !The range kmin=1e-3 to kmax=1e4 is necessary to compare to HMcode
     nk=128
     kmin=1e-3
     kmax=1e2
     CALL fill_array(log(kmin),log(kmax),k,nk)
     k=exp(k)

     !Set the number of redshifts and range (linearly spaced) and convert z -> a
     nz=16
     zmin=0.
     zmax=4.
     CALL fill_array(zmin,zmax,a,nz)
     a=1./(1.+a)
     na=nz

     ip=-1 !Set DMONLY profiles 
     CALL calculate_HMx(ihm,ip,mmin,mmax,k,nk,a,na,powa_lin,powa_2h,powa_1h,powa_full,cosm,verbose)

     base='data/power'
     CALL write_power_a_multiple(k,a,powa_lin,powa_2h,powa_1h,powa_full,nk,na,base,verbose)

  ELSE IF(imode==17) THEN

     !Assigns the cosmological model
     icosmo=1
     CALL assign_cosmology(icosmo,cosm,verbose)

     !Normalises power spectrum (via sigma_8) and fills sigma(R) look-up tables
     CALL print_cosmology(cosm)

     !Set number of k points and k range (log spaced)
     !The range kmin=1e-3 to kmax=1e4 is necessary to compare to HMcode
     nk=128
     kmin=1e-3
     kmax=1e2
     CALL fill_array(log(kmin),log(kmax),k,nk)
     k=exp(k)

     !Set the number of redshifts and range (linearly spaced) and convert z -> a
     nz=16
     zmin=0.
     zmax=4.
     CALL fill_array(zmin,zmax,a,nz)
     a=1./(1.+a)
     na=nz

     !Choose the field types
     DO i=1,2
        WRITE(*,*) 'HMx_driver: Choose halo', i
        CALL set_halo_type(ip(i))
     END DO
     !WRITE(*,*) 'IP:', ip
     !STOP

     !User chooses halo model   
     CALL calculate_HMx(ihm,ip,mmin,mmax,k,nk,a,na,powa_lin,powa_2h,powa_1h,powa_full,cosm,verbose)

     base='data/power'
     CALL write_power_a_multiple(k,a,powa_lin,powa_2h,powa_1h,powa_full,nk,na,base,verbose)

  ELSE IF(imode==18) THEN

     !Assigns the cosmological model
     icosmo=1
     CALL assign_cosmology(icosmo,cosm,verbose)

     !Normalises power spectrum (via sigma_8) and fills sigma(R) look-up tables
     !CALL initialise_cosmology(verbose,cosm)
     CALL print_cosmology(cosm)

     !Set number of k points and k range (log spaced)
     !The range kmin=1e-3 to kmax=1e4 is necessary to compare to HMcode
     nk=128
     kmin=1e-3
     kmax=1e2
     CALL fill_array(log(kmin),log(kmax),k,nk)
     k=exp(k)

     !Set the number of redshifts and range (linearly spaced) and convert z -> a
     nz=16
     zmin=0.
     zmax=4.
     CALL fill_array(zmin,zmax,a,nz)
     a=1./(1.+a)
     na=nz

     !Select field type for bias study
     CALL set_halo_type(ihalo)

     DO j=1,3

        IF(j==1) THEN
           !DMONLY-DMONLY
           ip(1)=-1
           ip(2)=-1
           base='bias/power_mm'
        ELSE IF(j==2) THEN
           !DMONLY-field
           ip(1)=ihalo
           ip(2)=-1
           base='bias/power_mf'
        ELSE IF(j==3) THEN
           !field-field
           ip(1)=ihalo
           ip(2)=ihalo
           base='bias/power_ff'
        END IF
        
        CALL calculate_HMx(ihm,ip,mmin,mmax,k,nk,a,na,powa_lin,powa_2h,powa_1h,powa_full,cosm,verbose)
        CALL write_power_a_multiple(k,a,powa_lin,powa_2h,powa_1h,powa_full,nk,na,base,verbose)

     END DO


  ELSE IF(imode==2 .OR. imode==15 .OR. imode==16) THEN

     !Make cross/auto power spectra of all different components of haloes as well as pressure

     !Generic hydro
     IF(imode==2) THEN     
        
        !Only do one 'model' here
        n=1

        !Set the redshift
        nz=4
        ALLOCATE(z_tab(nz))
        z_tab(1)=0.0
        z_tab(2)=0.5
        z_tab(3)=1.0
        z_tab(4)=2.0
     
        !Set number of k points and k range (log spaced)
        nk=128
        kmin=1e-3
        kmax=1e2
        CALL fill_array(log(kmin),log(kmax),k,nk)
        k=exp(k)

     !cosmo-OWLS
     ELSE IF(imode==15) THEN

        !Do from REF, NOCOOL, AGN, AGN 8.5, AGN 8.7
        n=5

        !Set the redshift
        nz=1
        ALLOCATE(z_tab(nz))
        z_tab(1)=0.

        !Get the k values from the simulation measured P(k)
        infile='/Users/Mead/Physics/cosmo-OWLS/power/N800/DMONLY_all_all_power.dat'
        CALL get_k_values(infile,k,nk)

     !BAHAMAS
     ELSE IF(imode==16) THEN
       
        !Do AGN, AGN-lo and AGN-hi
        n=3

        !Set the redshift
        nz=4
        ALLOCATE(z_tab(nz))
        z_tab(1)=0.0
        z_tab(2)=0.5
        z_tab(3)=1.0
        z_tab(4)=2.0

        !Get the k values from the simulation measured P(k)
        infile='/Users/Mead/Physics/BAHAMAS/power/M1024/DMONLY_nu0_L400N1024_WMAP9_snap32_all_all_power.dat'        
        CALL get_k_values(infile,k,nk)
        
     END IF

     !Allocate the arrays for P(k)
     ALLOCATE(pow_lin(nk),pow_2h(nk),pow_1h(nk),pow_full(nk))

     DO iowl=1,n

        DO j=1,nz

           z=z_tab(j)

           !Assigns the cosmological model
           IF(imode==2)  icosmo=4
           IF(imode==15) icosmo=2
           IF(imode==16) icosmo=4
           CALL assign_cosmology(icosmo,cosm,verbose)

           IF(imode==15) THEN
           
              !cosmo-OWLS
              IF(imode==15 .AND. iowl==1) THEN
                 name='REF'
                 fname=name
                 !From my fitting by eye
                 cosm%alpha=2.
                 cosm%eps=1.
                 cosm%Gamma=1.24
                 cosm%M0=1e13
                 cosm%Astar=0.055
              ELSE IF(imode==15 .AND. iowl==2) THEN
                 name='NOCOOL'
                 fname=name
                 !From my fitting by eye
                 cosm%alpha=2.
                 cosm%eps=1.
                 cosm%Gamma=1.1
                 cosm%M0=0.
                 cosm%Astar=0.
              ELSE IF(imode==15 .AND. iowl==3) THEN
                 name='AGN'
                 fname=name
                 !From Tilman's preliminary results
                 cosm%alpha=0.52
                 cosm%eps=1.
                 cosm%Gamma=1.17
                 cosm%M0=1.047e14
                 cosm%Astar=0.02
              ELSE IF(imode==15 .AND. iowl==4) THEN
                 name='AGN 8.5'
                 fname='AGN8p5'
                 !From Tilman's preliminary results
                 cosm%alpha=0.56
                 cosm%eps=1.
                 cosm%Gamma=1.19
                 cosm%M0=3.548e14
                 cosm%Astar=0.01
              ELSE IF(imode==15 .AND. iowl==5) THEN
                 name='AGN 8.7'
                 fname='AGN8p7'
                 !From Tilman's preliminary results
                 cosm%alpha=0.53
                 cosm%eps=1.
                 cosm%Gamma=1.21
                 cosm%M0=7.586e14
                 cosm%Astar=0.01
              END IF

              !Need to do this after cosmological parameters are changed
              CALL init_cosmology(cosm)
              CALL print_cosmology(cosm)

           END IF

           !BAHAMAS
           IF(imode==16) THEN

              ! Set the halo model
              ihm=14

              IF(iowl==1) THEN

                 ! Simulation name and file name
                 name='AGN'
                 fname='AGN'
  
                 ! Best z=0 fit on 21/06/2018
                 IF(ihm==4) THEN
                    
                    IF(z==0.) THEN
                       cosm%alpha=0.379
                       cosm%eps=10**(-0.064)
                       cosm%Gamma=1.206
                       cosm%M0=10**(13.818)
                       cosm%Astar=0.029
                       cosm%whim=10**(5.752)
                    ELSE IF(z==0.5) THEN
                       cosm%alpha=0.537
                       cosm%eps=10**(-0.213)
                       cosm%Gamma=1.163
                       cosm%M0=10**(13.963)
                       cosm%Astar=0.024
                       cosm%whim=10**(5.749)
                    ELSE IF(z==1. .OR. z==2.) THEN
                       cosm%alpha=0.664
                       cosm%eps=10**(-0.394)
                       cosm%Gamma=1.152
                       cosm%M0=10**(13.973)
                       cosm%Astar=0.007
                       cosm%whim=10**(5.648)
                    END IF
                    
                 ELSE IF(ihm==6) THEN
                    
                    IF(z==0.) THEN
                       cosm%alpha=0.195
                       cosm%eps=10**(-0.484)
                       cosm%Gamma=1.399
                       cosm%M0=10**(13.807)
                       cosm%Astar=0.020
                       cosm%whim=10**(6.049)
                    ELSE IF(z==0.5) THEN
                       cosm%alpha=0.619
                       cosm%eps=10**(-0.309)
                       cosm%Gamma=1.507
                       cosm%M0=10**(14.937)
                       cosm%Astar=0.021
                       cosm%whim=10**(5.987)
                    ELSE IF(z==1. .OR. z==2.) THEN
                       cosm%alpha=0.384
                       cosm%eps=10**(-0.475)
                       cosm%Gamma=1.183
                       cosm%M0=10**(14.562)
                       cosm%Astar=0.017
                       cosm%whim=10**(5.796)
                    END IF
                    
                 ELSE IF(ihm==3 .OR. ihm==14) THEN
                    
                     IF(z==0.) THEN
                       cosm%alpha=0.428
                       cosm%eps=10**(0.015)
                       cosm%Gamma=1.287
                       cosm%M0=10**(13.233)
                       cosm%Astar=0.030
                       cosm%whim=10**(5.404)
                    ELSE IF(z==0.5) THEN
                       cosm%alpha=0.742
                       cosm%eps=10**(0.148)
                       cosm%Gamma=1.516
                       cosm%M0=10**(12.688)
                       cosm%Astar=0.026
                       cosm%whim=10**(5.531)
                    ELSE IF(z==1. .OR. z==2.) THEN
                       cosm%alpha=0.830
                       cosm%eps=10**(0.169)
                       cosm%Gamma=1.487
                       cosm%M0=10**(12.004)
                       cosm%Astar=0.024
                       cosm%whim=10**(5.643)
                    END IF
                    
                 END IF
                 
              ELSE IF(iowl==3) THEN

                 ! Simulation name and file name
                 name='AGN high'
                 fname='AGN-hi'

                 ! Best z=0 fit on 21/06/2018
                 IF(ihm==4) THEN

                    IF(z==0.) THEN
                       cosm%alpha=0.421
                       cosm%eps=10**(-0.160)
                       cosm%Gamma=1.211
                       cosm%M0=10**(14.313)
                       cosm%Astar=0.026
                       cosm%whim=10**(5.800)
                    ELSE IF(z==0.5) THEN
                       cosm%alpha=0.617
                       cosm%eps=10**(-0.491)
                       cosm%Gamma=1.170
                       cosm%M0=10**(14.451)
                       cosm%Astar=0.022
                       cosm%whim=10**(5.848)
                    ELSE IF(z==1. .OR. z==2.) THEN
                       cosm%alpha=1.038
                       cosm%eps=10**(-0.455)
                       cosm%Gamma=1.397
                       cosm%M0=10**(14.560)
                       cosm%Astar=0.007
                       cosm%whim=10**(5.780)
                    END IF

                 ELSE IF(ihm==6) THEN

                    IF(z==0.) THEN
                       cosm%alpha=0.302
                       cosm%eps=10**(-0.429)
                       cosm%Gamma=1.339
                       cosm%M0=10**(14.662)
                       cosm%Astar=0.026
                       cosm%whim=10**(6.080)
                    ELSE IF(z==0.5) THEN
                       cosm%alpha=0.751
                       cosm%eps=10**(-0.013)
                       cosm%Gamma=1.502
                       cosm%M0=10**(14.958)
                       cosm%Astar=0.014
                       cosm%whim=10**(5.959)
                    ELSE IF(z==1. .OR. z==2.) THEN
                       cosm%alpha=0.371
                       cosm%eps=10**(-0.127)
                       cosm%Gamma=1.101
                       cosm%M0=10**(14.966)
                       cosm%Astar=0.018
                       cosm%whim=10**(6.040)
                    END IF

                 ELSE IF(ihm==3 .OR. ihm==14) THEN

                    IF(z==0.) THEN
                       cosm%alpha=0.528
                       cosm%eps=10**(0.038)
                       cosm%Gamma=1.505
                       cosm%M0=10**(13.638)
                       cosm%Astar=0.027
                       cosm%whim=10**(5.078)
                    ELSE IF(z==0.5) THEN
                       cosm%alpha=0.742
                       cosm%eps=10**(0.125)
                       cosm%Gamma=1.547
                       cosm%M0=10**(13.481)
                       cosm%Astar=0.024
                       cosm%whim=10**(5.786)
                    ELSE IF(z==1. .OR. z==2.) THEN
                       cosm%alpha=0.918
                       cosm%eps=10**(0.289)
                       cosm%Gamma=1.996
                       cosm%M0=10**(13.022)
                       cosm%Astar=0.022
                       cosm%whim=10**(5.849)
                    END IF

                 END IF
                 
              ELSE IF(iowl==2) THEN

                 ! Simulation name and file name
                 name='AGN low'
                 fname='AGN-lo'
                 
                 ! Best z=0 21/06/2018
                 IF(ihm==4) THEN

                    IF(z==0.) THEN
                       cosm%alpha=0.353
                       cosm%eps=10**(-0.020)
                       cosm%Gamma=1.199
                       cosm%M0=10**(13.511)
                       cosm%Astar=0.031
                       cosm%whim=10**(5.764)
                    ELSE IF(z==0.5) THEN
                       cosm%alpha=0.492
                       cosm%eps=10**(-0.107)
                       cosm%Gamma=1.159
                       cosm%M0=10**(13.662)
                       cosm%Astar=0.025
                       cosm%whim=10**(5.720)
                    ELSE IF(z==1. .OR. z==2.) THEN
                       cosm%alpha=0.432
                       cosm%eps=10**(-0.424)
                       cosm%Gamma=1.089
                       cosm%M0=10**(13.687)
                       cosm%Astar=0.010
                       cosm%whim=10**(5.786)
                    END IF

                 ELSE IF(ihm==6) THEN

                    IF(z==0.) THEN
                       cosm%alpha=0.181
                       cosm%eps=10**(-0.489)
                       cosm%Gamma=1.432
                       cosm%M0=10**(13.632)
                       cosm%Astar=0.023
                       cosm%whim=10**(6.099)
                    ELSE IF(z==0.5) THEN
                       cosm%alpha=0.420
                       cosm%eps=10**(-0.413)
                       cosm%Gamma=1.467
                       cosm%M0=10**(14.606)
                       cosm%Astar=0.021
                       cosm%whim=10**(5.874)
                    ELSE IF(z==1. .OR. z==2.) THEN
                       cosm%alpha=0.715
                       cosm%eps=10**(-0.405)
                       cosm%Gamma=1.154
                       cosm%M0=10**(14.786)
                       cosm%Astar=0.019
                       cosm%whim=10**(5.434)
                    END IF

                 ELSE IF(ihm==3 .OR. ihm==14) THEN

                    IF(z==0.) THEN
                       cosm%alpha=0.409
                       cosm%eps=10**(0.045)
                       cosm%Gamma=1.275
                       cosm%M0=10**(12.737)
                       cosm%Astar=0.032
                       cosm%whim=10**(5.357)
                    ELSE IF(z==0.5) THEN
                       cosm%alpha=0.698
                       cosm%eps=10**(0.159)
                       cosm%Gamma=1.393
                       cosm%M0=10**(12.012)
                       cosm%Astar=0.028
                       cosm%whim=10**(5.491)
                    ELSE IF(z==1. .OR. z==2.) THEN
                       cosm%alpha=0.832
                       cosm%eps=10**(0.312)
                       cosm%Gamma=1.344
                       cosm%M0=10**(12.020)
                       cosm%Astar=0.025
                       cosm%whim=10**(5.723)
                    END IF

                 END IF

              END IF

              !Need to do this after cosmological parameters are changed  
              CALL init_cosmology(cosm)
              CALL print_cosmology(cosm)

           END IF

           IF(imode==15) WRITE(*,*) 'Comparing to OWLS model: ', TRIM(name)
           IF(imode==16) WRITE(*,*) 'Comparing to BAHAMAS model: ', TRIM(name)
            
           !Initiliasation for the halomodel calcualtion
           CALL init_halomod(ihm,mmin,mmax,z,hmod,cosm,verbose)

           !Runs the diagnostics
           IF(imode==2) THEN
              dir='diagnostics'
              CALL halo_diagnostics(z,hmod,cosm,dir)
              CALL halo_definitions(z,hmod,cosm,dir)
              CALL halo_properties(z,hmod,dir)
           END IF

           IF(imode==2) THEN
              !File base and extension
              IF(j==1) base='hydro/power_z0.0_'
              IF(j==2) base='hydro/power_z0.5_'
              IF(j==3) base='hydro/power_z1.0_'
              IF(j==4) base='hydro/power_z2.0_'
              mid=''
              ext='.dat'
           ELSE IF(imode==15) THEN
              base='cosmo-OWLS/power_'//TRIM(fname)//'_'
              mid=''
              ext='.dat'
           ELSE IF(imode==16) THEN
              IF(j==1) base='BAHAMAS/power_'//TRIM(fname)//'_z0.0_'
              IF(j==2) base='BAHAMAS/power_'//TRIM(fname)//'_z0.5_'
              IF(j==3) base='BAHAMAS/power_'//TRIM(fname)//'_z1.0_'
              IF(j==4) base='BAHAMAS/power_'//TRIM(fname)//'_z2.0_'
              mid=''
              ext='.dat'
           END IF

           !Dark-matter only
           IF(imode==2) THEN
              IF(j==1) outfile='hydro/power_z0.0.dat'
              IF(j==2) outfile='hydro/power_z0.5.dat'
              IF(j==3) outfile='hydro/power_z1.0.dat'
              IF(j==4) outfile='hydro/power_z2.0.dat'
           ELSE IF(imode==15) THEN
              outfile='cosmo-OWLS/power_DMONLY_00.dat'
           ELSE IF(imode==16) THEN
              IF(j==1) outfile='BAHAMAS/power_DMONLY_z0.0_00.dat'
              IF(j==2) outfile='BAHAMAS/power_DMONLY_z0.5_00.dat'
              IF(j==3) outfile='BAHAMAS/power_DMONLY_z1.0_00.dat'
              IF(j==4) outfile='BAHAMAS/power_DMONLY_z2.0_00.dat'
           END IF
           WRITE(*,*) -1, -1, TRIM(outfile)
           CALL calculate_halomod(-1,-1,k,nk,z,pow_lin,pow_2h,pow_1h,pow_full,hmod,cosm,verbose)
           CALL write_power(k,pow_lin,pow_2h,pow_1h,pow_full,nk,outfile,verbose)

           !Loop over matter types and do auto- and cross-spectra
           DO j1=0,6
              DO j2=j1,6

                 !Skip for the bound- and free-gas spectra, fuck 'em
                 IF(j1==4 .OR. j1==5) CYCLE
                 IF(j2==4 .OR. j2==5) CYCLE

                 !Fix output file and write to screen
                 outfile=number_file2(base,j1,mid,j2,ext)
                 WRITE(*,*) j1, j2, TRIM(outfile)

                 !Do the calculation and write P(k) to disk
                 CALL calculate_halomod(j1,j2,k,nk,z,pow_lin,pow_2h,pow_1h,pow_full,hmod,cosm,.FALSE.)
                 CALL write_power(k,pow_lin,pow_2h,pow_1h,pow_full,nk,outfile,.FALSE.)

              END DO
           END DO

        END DO

     END DO

  ELSE IF(imode==3) THEN

     !Assigns the cosmological model
     CALL assign_cosmology(icosmo,cosm,verbose)

     !TODO: Remove this
     !cosm%Gamma=1.2

     !Normalises power spectrum (via sigma_8) and fills sigma(R) look-up tables
     CALL print_cosmology(cosm)

     !Loop over redshifts
     DO j=1,4

        !Set the redshift
        IF(j==1) z=0.0
        IF(j==2) z=0.5
        IF(j==3) z=1.0
        IF(j==4) z=2.0

        !Initiliasation for the halomodel calcualtion
        CALL init_halomod(ihm,mmin,mmax,z,hmod,cosm,verbose)

        !Runs the diagnostics
        dir='diagnostics'
        CALL halo_diagnostics(z,hmod,cosm,dir)
        CALL halo_definitions(z,hmod,cosm,dir)
        CALL halo_properties(z,hmod,dir)

     END DO

  ELSE IF(imode==4) THEN

     !Ignore this, only useful for bug tests
     CALL RNG_set(0)

     !Set number of k points and k range (log spaced)
     nk=128
     kmin=1e-3
     kmax=1e2
     CALL fill_array(log(kmin),log(kmax),k,nk)
     k=exp(k)
     ALLOCATE(pow_lin(nk),pow_2h(nk),pow_1h(nk),pow_full(nk))

     !Assigns the cosmological model
     CALL assign_cosmology(icosmo,cosm,verbose)

     !Sets the redshift
     z=0.

     DO

        CALL random_baryon_parameters(cosm)

        CALL init_cosmology(cosm)

        !Normalises power spectrum (via sigma_8) and fills sigma(R) look-up tables
        CALL print_cosmology(cosm)

        !Initiliasation for the halomodel calcualtion
        CALL init_halomod(ihm,mmin,mmax,z,hmod,cosm,verbose)

        WRITE(*,fmt='(A12,6F12.6)') 'Parameters:', cosm%alpha, cosm%eps, cosm%Gamma, log10(cosm%M0), cosm%Astar, log10(cosm%whim)

        !Do the halo-model calculation
        DO j1=0,6
           DO j2=j1,6
              IF(j1==1 .OR. j1==4 .OR. j1==5) CYCLE
              IF(j2==1 .OR. j2==4 .OR. j2==5) CYCLE
              WRITE(*,*) 'Halo types:', j1, j2
              CALL calculate_halomod(j1,j2,k,nk,z,pow_lin,pow_2h,pow_1h,pow_full,hmod,cosm,.FALSE.)
              DO i=1,nk
                 IF(ISNAN(pow_full(i))) STOP 'HMX_DRIVER: Error, NaN found in pow_full array'
              END DO
           END DO
        END DO

        WRITE(*,*)

     END DO

  ELSE IF(imode==5) THEN

     STOP 'HMX_DRIVER: Error, imode=5 not supported any more'

     !Create spectra for 'all matter' and 'electron pressure' and their cross

     !Set number of k points and k range (log spaced)
     nk=128
     kmin=1e-3
     kmax=1e2
     CALL fill_array(log(kmin),log(kmax),k,nk)
     k=exp(k)
     ALLOCATE(pow_lin(nk),pow_2h(nk),pow_1h(nk),pow_full(nk))

     !Assigns the cosmological model
     icosmo=2
     CALL assign_cosmology(icosmo,cosm,verbose)

     !Normalises power spectrum (via sigma_8) and fills sigma(R) look-up tables
     !CALL initialise_cosmology(verbose,cosm)
     CALL print_cosmology(cosm)

     !File base and extension
     base='data/power_'
     ext='.dat'

     DO i=1,4

        !Set the redshift
        IF(i==1) THEN
           z=0.0
           red='z0.0'
        ELSE IF(i==2) THEN
           z=0.5
           red='z0.5'
        ELSE IF(i==3) THEN
           z=1.0
           red='z1.0'
        ELSE IF(i==4) THEN
           z=2.0
           red='z2.0'
        END IF

        !Initiliasation for the halomodel calcualtion
        CALL init_halomod(ihm,mmin,mmax,z,hmod,cosm,verbose)

        !Runs the diagnostics
        dir='diagnostics'
        CALL halo_diagnostics(z,hmod,cosm,dir)
        CALL halo_definitions(z,hmod,cosm,dir)
        CALL halo_properties(z,hmod,dir)

        !Do the calculation
        DO j=0,3

           IF(j==0) THEN
              !DMONLY
              j1=-1
              j2=-1
              outfile='data/power_'//TRIM(red)//TRIM(ext)
           ELSE IF(j==1) THEN
              !matter - matter
              j1=0
              j2=0
              outfile='dd'
           ELSE IF(j==2) THEN
              !matter - electron pressure
              j1=0
              j2=6
              outfile='dp'
           ELSE IF(j==3) THEN
              !electron pressure - electron pressure
              j1=6
              j2=6
              outfile='pp'
           END IF

           IF(j .NE. 0) outfile=TRIM(base)//TRIM(outfile)//'_'//TRIM(red)//TRIM(ext)

           WRITE(*,fmt='(3I5,A30)') j, j1, j2, TRIM(outfile)

           CALL calculate_halomod(j1,j2,k,nk,z,pow_lin,pow_2h,pow_1h,pow_full,hmod,cosm,.FALSE.)
           CALL write_power(k,pow_lin,pow_2h,pow_1h,pow_full,nk,outfile,.FALSE.)

        END DO

        WRITE(*,*)

     END DO

  ELSE IF(imode==6) THEN

     !n(z) normalisation check

     WRITE(*,*) 'HMx_DRIVER: Checking n(z) functions'
     WRITE(*,*)

     nnz=7
     DO i=1,nnz
        IF(i==1) nz=1
        IF(i==2) nz=4
        IF(i==3) nz=5
        IF(i==4) nz=6
        IF(i==5) nz=7
        IF(i==6) nz=8
        IF(i==7) nz=9
        WRITE(*,*) 'HMx_DRIVER: n(z) number:', nz
        WRITE(*,*)
        CALL get_nz(nz,lens)
        WRITE(*,*) 'HMx_DRIVER: n(z) integral (linear):', integrate_table(lens%z_nz,lens%nz,lens%nnz,1,lens%nnz,1)
        WRITE(*,*) 'HMx_DRIVER: n(z) integral (quadratic):', integrate_table(lens%z_nz,lens%nz,lens%nnz,1,lens%nnz,2)
        WRITE(*,*) 'HMx_DRIVER: n(z) integral (cubic):', integrate_table(lens%z_nz,lens%nz,lens%nnz,2,lens%nnz,3)
        WRITE(*,*)
     END DO

  ELSE IF(imode==7 .OR. imode==8 .OR. imode==9 .OR. imode==10 .OR. imode==11) THEN

     !General stuff for all cross correlations

     !Set the fields
     ix=-1
     CALL set_xcorr_type(ix,ip)

     !Assign the cosmological model
     CALL assign_cosmology(icosmo,cosm,verbose)

     !Set the k range
     kmin=1e-3
     kmax=1e2
     nk=128

     !Set the z range
     !amin=scale_factor_z(cosm%z_cmb) !Problems with one-halo term if amin is less than 0.1
     amin=0.1
     amax=1.
     na=16

     !Set number of k points and k range (log spaced)
     !Also z points and z range (linear)
     !Also P(k,z)
     CALL fill_array(log(kmin),log(kmax),k,nk)
     k=exp(k)
     CALL fill_array(amin,amax,a,na)
     ALLOCATE(powa_lin(nk,na),powa_1h(nk,na),powa_2h(nk,na),powa_full(nk,na),powa(nk,na))

     !Set the ell range
     lmin=1
     lmax=1e4 !Problems if this is pushed up to 10^5
     nl=128

     !Allocate arrays for l and C(l)
     CALL fill_array(log(lmin),log(lmax),ell,nl)
     ell=exp(ell)
     ALLOCATE(Cell(nl))

     !Set the angular arrays in degrees
     !Allocate arrays for theta and xi(theta)
     IF(ixi) THEN
        thmin=0.01
        thmax=10.
        nth=128
        CALL fill_array(log(thmin),log(thmax),theta,nth)
        theta=exp(theta)
        ALLOCATE(xi(3,nth))
     END IF

     WRITE(*,*) 'HMx_DRIVER: Cross-correlation information'
     WRITE(*,*) 'HMx_DRIVER: output directiory: ', TRIM(dir)
     WRITE(*,*) 'HMx_DRIVER: Profile type 1: ', TRIM(halo_type(ip(1)))
     WRITE(*,*) 'HMx_DRIVER: Profile type 2: ', TRIM(halo_type(ip(2)))
     WRITE(*,*) 'HMx_DRIVER: cross-correlation type 1: ', TRIM(xcorr_type(ix(1)))
     WRITE(*,*) 'HMx_DRIVER: cross-correlation type 2: ', TRIM(xcorr_type(ix(2)))
     WRITE(*,*) 'HMx_DRIVER: P(k) minimum k [h/Mpc]:', REAL(kmin)
     WRITE(*,*) 'HMx_DRIVER: P(k) maximum k [h/Mpc]:', REAL(kmax)
     WRITE(*,*) 'HMx_DRIVER: minimum a:', REAL(amin)
     WRITE(*,*) 'HMx_DRIVER: maximum a:', REAL(amax)
     WRITE(*,*) 'HMx_DRIVER: minimum ell:', REAL(lmin)
     WRITE(*,*) 'HMx_DRIVER: maximum ell:', REAL(lmax)
     WRITE(*,*)     

     IF(imode==7) THEN

        !Normalises power spectrum (via sigma_8) and fills sigma(R) look-up tables
        !CALL initialise_cosmology(verbose,cosm)
        CALL print_cosmology(cosm)

        !Initialise the lensing part of the calculation
        !CALL initialise_distances(verbose,cosm)
        CALL write_distances(cosm)

        !Write out diagnostics
        CALL init_halomod(ihm,mmin,mmax,z,hmod,cosm,verbose)

        !dir='diagnostics'
        !CALL halo_diagnostics(z,hmod,cosm,dir)
        !CALL halo_definitions(z,hmod,dir)
        !CALL halo_properties(z,hmod,dir)

        CALL calculate_HMx(ihm,ip,mmin,mmax,k,nk,a,na,powa_lin,powa_2h,powa_1h,powa_full,cosm,verbose)

!!$        !Fix the one-halo term P(k) to be a constant
!!$        DO j=1,na
!!$           DO i=1,nk
!!$              !IF(j==1) WRITE(*,*) i, k(i), powz(3,i,j)
!!$              powa(3,i,j)=powa(3,1,j)*(k(i)/k(1))**3
!!$              !powz(3,i,j)=(k(i)/k(1))**3
!!$              !IF(j==1) WRITE(*,*) i, k(i), powz(3,i,j)
!!$           END DO
!!$        END DO

        !Output directory
        dir='data/'
        base=TRIM(dir)//'power'
        CALL write_power_a_multiple(k,a,powa_lin,powa_2h,powa_1h,powa_full,nk,na,base,verbose)

        !Fill out the projection kernels
        CALL fill_projection_kernels(ix,proj,cosm)
        CALL write_projection_kernels(proj,cosm)

        !Set the distance range for the Limber integral
        r1=0. !100.
        r2=maxdist(proj)!proj%rs

        !Write to screen
        WRITE(*,*) 'HMx_DRIVER: Computing C(l)'
        WRITE(*,*) 'HMx_DRIVER: ell min:', REAL(ell(1))
        WRITE(*,*) 'HMx_DRIVER: ell max:', REAL(ell(nl))
        WRITE(*,*) 'HMx_DRIVER: number of ell:', nl
        WRITE(*,*) 'HMx_DRIVER: lower limit of Limber integral [Mpc/h]:', REAL(r1)
        WRITE(*,*) 'HMx_DRIVER: upper limit of Limber integral [Mpc/h]:', REAL(r2)
        WRITE(*,*)

        !Loop over all types of C(l) to create
        DO j=1,4

           IF(ifull .AND. (j .NE. 4)) CYCLE
           !IF(j==3) CYCLE !Skip the fucking one-halo term

           !Write information to screen
           IF(j==1) THEN
              WRITE(*,*) 'HMx_DRIVER: Doing linear'
              outfile=TRIM(dir)//'cl_linear.dat'
              powa=powa_lin
           ELSE IF(j==2) THEN
              WRITE(*,*) 'HMx_DRIVER: Doing 2-halo'
              outfile=TRIM(dir)//'cl_2halo.dat'
              powa=powa_2h
           ELSE IF(j==3) THEN
              WRITE(*,*) 'HMx_DRIVER: Doing 1-halo'
              outfile=TRIM(dir)//'cl_1halo.dat'
              powa=powa_1h
           ELSE IF(j==4) THEN
              WRITE(*,*) 'HMx_DRIVER: Doing full'
              outfile=TRIM(dir)//'cl_full.dat'
              powa=powa_full
           END IF

           WRITE(*,*) 'HMx_DRIVER: Output: ', TRIM(outfile)

           !Actually calculate the C(ell)
           CALL calculate_Cell(r1,r2,ell,Cell,nl,k,a,powa,nk,na,proj,cosm)
           CALL write_Cell(ell,Cell,nl,outfile)

           IF(j==4) CALL Cell_contribution(r1,r2,k,a,powa,nk,na,proj,cosm)

           IF(ixi) THEN

              !Set xi output files
              IF(j==1) outfile=TRIM(dir)//'xi_linear.dat'
              IF(j==2) outfile=TRIM(dir)//'xi_2halo.dat'
              IF(j==3) outfile=TRIM(dir)//'xi_1halo.dat'
              IF(j==4) outfile=TRIM(dir)//'xi_full.dat'
              WRITE(*,*) 'HMx_DRIVER: Output: ', TRIM(outfile)

              !Actually calculate the xi(theta)
              CALL calculate_xi(theta,xi,nth,ell,Cell,nl,NINT(lmax))
              CALL write_xi(theta,xi,nth,outfile)

           END IF

        END DO
        WRITE(*,*) 'HMx_DRIVER: Done'
        WRITE(*,*)

     ELSE IF(imode==8) THEN

        !Assess cross-correlation as a function of cosmology

        !Loop over cosmology
        sig8min=0.7
        sig8max=0.9
        ncos=5
        DO i=1,ncos

           !cosm%sig8=sig8min+(sig8max-sig8min)*float(i-1)/float(ncos-1)
           cosm%sig8=progression(sig8min,sig8max,i,ncos)
           CALL init_cosmology(cosm)

           CALL print_cosmology(cosm)
           !CALL initialise_distances(verbose,cosm)

           CALL calculate_HMx(ihm,ip,mmin,mmax,k,nk,a,na,powa_lin,powa_2h,powa_1h,powa_full,cosm,verbose)

           !Fill out the projection kernels
           CALL fill_projection_kernels(ix,proj,cosm)
           CALL write_projection_kernels(proj,cosm)

           !Now do the C(l) calculations
           !Set l range, note that using Limber and flat-sky for sensible results lmin to ~10
           CALL fill_array(log(lmin),log(lmax),ell,nl)
           ell=exp(ell)
           IF(ALLOCATED(Cell)) DEALLOCATE(Cell)
           ALLOCATE(Cell(nl))

           !Write to screen
           WRITE(*,*) 'HMx_DRIVER: Computing C(l)'
           WRITE(*,*) 'HMx_DRIVER: ell min:', REAL(ell(1))
           WRITE(*,*) 'HMx_DRIVER: ell max:', REAL(ell(nl))
           WRITE(*,*) 'HMx_DRIVER: number of ell:', nl
           WRITE(*,*)

           !Loop over all types of C(l) to create
           dir='data/'
           base=TRIM(dir)//'cosmology_'
           DO j=1,4
              IF(j==1) THEN
                 WRITE(*,*) 'HMx_DRIVER: Doing C(l) linear'
                 ext='_cl_linear.dat'
                 powa=powa_lin
              ELSE IF(j==2) THEN
                 WRITE(*,*) 'HMx_DRIVER: Doing C(l) 2-halo'
                 ext='_cl_2halo.dat'
                 powa=powa_2h
              ELSE IF(j==3) THEN
                 WRITE(*,*) 'HMx_DRIVER: Doing C(l) 1-halo'
                 ext='_cl_1halo.dat'
                 powa=powa_1h
              ELSE IF(j==4) THEN
                 WRITE(*,*) 'HMx_DRIVER: Doing C(l) full'
                 ext='_cl_full.dat'
                 powa=powa_full
              END IF           
              outfile=number_file(base,i,ext)
              WRITE(*,*) 'HMx_DRIVER: Output: ', TRIM(outfile)
              !Actually calculate the C(l)
              CALL calculate_Cell(0.,maxdist(proj),ell,Cell,nl,k,a,powa,nk,na,proj,cosm)
              CALL write_Cell(ell,Cell,nl,outfile)
           END DO
           WRITE(*,*) 'HMx_DRIVER: Done'
           WRITE(*,*)

        END DO

     ELSE IF(imode==9) THEN

        !Breakdown cross-correlation in terms of mass

        !Normalises power spectrum (via sigma_8) and fills sigma(R) look-up tables
        !CALL initialise_cosmology(verbose,cosm)
        CALL print_cosmology(cosm)

        !Initialise the lensing part of the calculation
        !CALL initialise_distances(verbose,cosm)
        CALL write_distances(cosm)

        !Fill out the projection kernels
        CALL fill_projection_kernels(ix,proj,cosm)
        CALL write_projection_kernels(proj,cosm)

        DO i=0,6
           IF(icumulative .EQV. .FALSE.) THEN
              !Set the mass intervals
              IF(i==0) THEN
                 m1=mmin
                 m2=mmax
              ELSE IF(i==1) THEN
                 m1=mmin
                 m2=1e11
              ELSE IF(i==2) THEN
                 m1=1e11
                 m2=1e12
              ELSE IF(i==3) THEN
                 m1=1e12
                 m2=1e13
              ELSE IF(i==4) THEN
                 m1=1e13
                 m2=1e14
              ELSE IF(i==5) THEN
                 m1=1e14
                 m2=1e15
              ELSE IF(i==6) THEN
                 m1=1e15
                 m2=1e16
              END IF
           ELSE
              !Set the mass intervals
              IF(i==0) THEN
                 m1=mmin
                 m2=mmax
              ELSE IF(i==1) THEN
                 m1=mmin
                 m2=1e11
              ELSE IF(i==2) THEN
                 m1=mmin
                 m2=1e12
              ELSE IF(i==3) THEN
                 m1=mmin
                 m2=1e13
              ELSE IF(i==4) THEN
                 m1=mmin
                 m2=1e14
              ELSE IF(i==5) THEN
                 m1=mmin
                 m2=1e15
              ELSE IF(i==6) THEN
                 m1=mmin
                 m2=1e16
              END IF
           END IF

           !Set the code to not 'correct' the two-halo power for missing
           !mass when doing the calcultion binned in halo mass
           !STOP 'HMx_DRIVER: Extreme caution here, need to set ip2h=0, but it is defined as parameter in HMx.f90'
           !IF((icumulative .EQV. .FALSE.) .AND. i>1) hmod%ip2h=0
           IF((icumulative .EQV. .TRUE.) .AND. i>0) hmod%ip2h=0
           
           WRITE(*,fmt='(A16)') 'HMx_DRIVER: Mass range'
           WRITE(*,fmt='(A16,I5)') 'HMx_DRIVER: Iteration:', i
           WRITE(*,fmt='(A21,2ES15.7)') 'HMx_DRIVER: M_min [Msun/h]:', m1
           WRITE(*,fmt='(A21,2ES15.7)') 'HMx_DRIVER: M_max [Msun/h]:', m2
           WRITE(*,*)

           !Loop over redshifts
           DO j=1,na

              z=redshift_a(a(j))

              !Initiliasation for the halomodel calcualtion
              CALL init_halomod(ihm,m1,m2,z,hmod,cosm,verbose)
              CALL calculate_halomod(ip(1),ip(2),k,nk,z,powa_lin(:,j),powa_2h(:,j),powa_1h(:,j),powa_full(:,j),hmod,cosm,verbose)

              !Write progress to screen
              IF(j==1) THEN
                 WRITE(*,fmt='(A5,A7)') 'i', 'a'
                 WRITE(*,fmt='(A13)') '   ============'
              END IF
              WRITE(*,fmt='(I5,F8.3)') j, a(j)

           END DO
           WRITE(*,fmt='(A13)') '   ============'
           WRITE(*,*)

           dir='data/'
           IF(i==0) THEN
              outfile=TRIM(dir)//'power'
           ELSE
              base=TRIM(dir)//'mass_'
              mid='_'
              ext='_power'
              outfile=number_file2(base,NINT(log10(m1)),mid,NINT(log10(m2)),ext)
           END IF
           WRITE(*,*) 'HMx_DRIVER: File: ', TRIM(outfile)
           !CALL write_power_a(k,a,powa,nk,na,output)

           !Loop over all types of C(l) to create
           base=TRIM(dir)//'mass_'
           mid='_' 
           DO j=1,4

              !Skip the 1-halo C(l) because it takes ages (2017/02/06)
              IF(j==3) CYCLE

              !Set output files
              IF(j==1) THEN
                 powa=powa_lin
                 outfile=TRIM(dir)//'cl_linear.dat'
                 ext='_cl_linear.dat'
              ELSE IF(j==2) THEN
                 powa=powa_2h
                 outfile=TRIM(dir)//'cl_2halo.dat'
                 ext='_cl_2halo.dat'
              ELSE IF(j==3) THEN
                 powa=powa_1h
                 outfile=TRIM(dir)//'cl_1halo.dat'
                 ext='_cl_1halo.dat'
              ELSE IF(j==4) THEN
                 powa=powa_full
                 outfile=TRIM(dir)//'cl_full.dat'
                 ext='_cl_full.dat'
              END IF

              IF(i>0) outfile=number_file2(base,NINT(log10(m1)),mid,NINT(log10(m2)),ext)

              WRITE(*,*) 'HMx_DRIVER: File: ', TRIM(outfile)

              CALL calculate_Cell(0.,maxdist(proj),ell,Cell,nl,k,a,powa,nk,na,proj,cosm)
              CALL write_Cell(ell,Cell,nl,outfile)

           END DO
           WRITE(*,*) 'HMx_DRIVER: Done'
           WRITE(*,*)

        END DO

     ELSE IF(imode==10) THEN

        !Break down cross-correlation in terms of redshift

        !Normalises power spectrum (via sigma_8) and fills sigma(R) look-up tables
        !CALL initialise_cosmology(verbose,cosm)
        CALL print_cosmology(cosm)

        CALL calculate_HMx(ihm,ip,mmin,mmax,k,nk,a,na,powa_lin,powa_2h,powa_1h,powa_full,cosm,verbose)

        !Initialise the lensing part of the calculation
        !CALL initialise_distances(verbose,cosm)
        CALL write_distances(cosm)

        !Fill out the projection kernels
        CALL fill_projection_kernels(ix,proj,cosm)
        CALL write_projection_kernels(proj,cosm)

        !Write to screen
        WRITE(*,*) 'HMx_DRIVER: Computing C(l)'
        WRITE(*,*) 'HMx_DRIVER: ell min:', REAL(ell(1))
        WRITE(*,*) 'HMx_DRIVER: ell max:', REAL(ell(nl))
        WRITE(*,*) 'HMx_DRIVER: number of ell:', nl
        WRITE(*,*)

        zmin=0.
        zmax=1.
        nz=8

        DO i=0,nz

           IF(i==0) THEN
              r1=0.
              r2=maxdist(proj)
           ELSE
              IF(icumulative .EQV. .FALSE.) THEN
                 z1=progression(zmin,zmax,i,nz)
              ELSE
                 z1=zmin
              END IF
              z2=zmin+(zmax-zmin)*float(i)/float(nz)
              a1=scale_factor_z(z1)
              a2=scale_factor_z(z2)
              r1=comoving_distance(a1,cosm)
              r2=comoving_distance(a2,cosm)
           END IF

           WRITE(*,*) 'HMx_DRIVER:', i
           IF(i>0) THEN
              WRITE(*,*) 'HMx_DRIVER: z1:', REAL(z1)
              WRITE(*,*) 'HMx_DRIVER: z2:', REAL(z2)
           END IF
           WRITE(*,*) 'HMx_DRIVER: r1 [Mpc/h]:', REAL(r1)
           WRITE(*,*) 'HMx_DRIVER: r2 [Mpc/h]:', REAL(r2)

           !Loop over all types of C(l) to create
           dir='data/'
           base=TRIM(dir)//'redshift_'
           mid='_'
           DO j=1,4

              !Set output files
              IF(j==1) THEN
                 ext='_cl_linear.dat'
                 outfile=TRIM(dir)//'cl_linear.dat'
                 powa=powa_lin
              ELSE IF(j==2) THEN
                 ext='_cl_2halo.dat'
                 outfile=TRIM(dir)//'cl_2halo.dat'
                 powa=powa_2h
              ELSE IF(j==3) THEN
                 ext='_cl_1halo.dat'
                 outfile=TRIM(dir)//'cl_1halo.dat'
                 powa=powa_1h
              ELSE IF(j==4) THEN
                 ext='_cl_full.dat'
                 outfile=TRIM(dir)//'cl_full.dat'
                 powa=powa_full
              END IF

              IF(i>0 .AND. (icumulative .EQV. .FALSE.)) THEN
                 outfile=number_file2(base,i-1,mid,i,ext)
              ELSE IF(i>0 .AND. icumulative) THEN
                 outfile=number_file2(base,0,mid,i,ext)
              END IF
              WRITE(*,*) 'HMx_DRIVER: Output: ', TRIM(outfile)

              !This crashes for the low r2 values for some reason
              !Only a problem if lmax ~ 10^5
              !STOP 'This crashes for the low r2 values for high ell for some reason - should debug'
              CALL calculate_Cell(r1,r2,ell,Cell,nl,k,a,powa,nk,na,proj,cosm)
              CALL write_Cell(ell,Cell,nl,outfile)

           END DO
           WRITE(*,*)

        END DO

        WRITE(*,*) 'HMx_DRIVER: Done'
        WRITE(*,*)

     ELSE IF(imode==11) THEN

        STOP 'HMx_DRIVER: Error, breakdown in radius is not supported yet'

     ELSE

        STOP 'HMx_DRIVER: Error, you have specified the mode incorrectly'

     END IF

  ELSE IF(imode==12) THEN

     !Project triad

     dir='data'

     !Assigns the cosmological model
     icosmo=3
     CALL assign_cosmology(icosmo,cosm,verbose)
     !cosm%alpha=0.485
     cosm%alpha=2.
     cosm%eps=10**0.103
     cosm%Gamma=1.212
     cosm%M0=10**13.836
     cosm%Astar=0.029
     cosm%whim=10**6.346

     !Normalises power spectrum (via sigma_8) and fills sigma(R) look-up tables
     CALL init_cosmology(cosm)
     CALL print_cosmology(cosm)

     !Initialise the lensing part of the calculation
     !CALL initialise_distances(verbose,cosm)
     CALL write_distances(cosm)

     !Set the ell range
     lmin=100.
     lmax=3000.
     nl=64

     !Allocate arrays for l and C(l)
     CALL fill_array(log(lmin),log(lmax),ell,nl)
     ell=exp(ell)
     ALLOCATE(Cell(nl))

     WRITE(*,*) 'HMx_DRIVER: Cross-correlation information'
     WRITE(*,*) 'HMx_DRIVER: output directiory: ', TRIM(dir)
     WRITE(*,*) 'HMx_DRIVER: minimum ell:', REAL(lmin)
     WRITE(*,*) 'HMx_DRIVER: maximum ell:', REAL(lmax)
     WRITE(*,*) 'HMx_DRIVER: number of ell:', nl
     WRITE(*,*)

     !Loop over the triad
     DO i=1,3

        IF(i==1) THEN
           !ix(1)=4 !CFHTLenS
           ix(1)=5 !KiDS
           ix(2)=3 !CMB
           outfile=TRIM(dir)//'/triad_Cl_gal-CMB.dat'
        ELSE IF(i==2) THEN
           ix(1)=3 !CMB
           ix(2)=2 !y
           outfile=TRIM(dir)//'/triad_Cl_CMB-y.dat'
        ELSE IF(i==3) THEN
           ix(1)=2 !y
           !ix(2)=4 !CFHTLenS
           ix(2)=5 !KiDS
           outfile=TRIM(dir)//'/triad_Cl_y-gal.dat'
        END IF

        CALL xcorr(ihm,ix,mmin,mmax,ell,Cell,nl,cosm,verbose)
        CALL write_Cell(ell,Cell,nl,outfile)

        WRITE(*,*) 'HMx_DRIVER: Done'
        WRITE(*,*)

     END DO

  ELSE IF(imode==13) THEN

     !Calculate the cross-correlation coefficient

     !Assign the cosmology
     CALL assign_cosmology(icosmo,cosm,verbose)

     !Normalises power spectrum (via sigma_8) and fills sigma(R) look-up tables
     !CALL initialise_cosmology(verbose,cosm)
     CALL print_cosmology(cosm)

     !Initialise the lensing part of the calculation
     !CALL initialise_distances(verbose,cosm)
     CALL write_distances(cosm)

     !Set the ell range and allocate arrays for l and C(l)
     lmin=1e0
     lmax=1e4 !Errors if this is increased to 10^5
     nl=64 
     CALL fill_array(log(lmin),log(lmax),ell,nl)
     ell=exp(ell)
     ALLOCATE(Cell(nl))

     dir='data'

     ixx=-1
     CALL set_xcorr_type(ixx,ip)

     DO i=1,3
        IF(i==1) THEN
           ix(1)=ixx(1)
           ix(2)=ixx(1)
           outfile=TRIM(dir)//'/cl_first.dat'
        ELSE IF(i==2) THEN
           ix(1)=ixx(2)
           ix(2)=ixx(2)
           outfile=TRIM(dir)//'/cl_second.dat'
        ELSE IF(i==3) THEN
           ix(1)=ixx(1)
           ix(2)=ixx(2)
           outfile=TRIM(dir)//'/cl_full.dat'
        END IF
        CALL xcorr(ihm,ix,mmin,mmax,ell,Cell,nl,cosm,verbose)
        CALL write_Cell(ell,Cell,nl,outfile)
     END DO

  ELSE IF(imode==14) THEN

     !Make power spectra as a function of baryon parameter variations

     !Number of values to try for each parameter
     m=9

     !Set number of k points and k range (log spaced)
     kmin=1e-3
     kmax=1e1
     nk=128
     CALL fill_array(log(kmin),log(kmax),k,nk)
     k=exp(k)
     ALLOCATE(pow_lin(nk),pow_2h(nk),pow_1h(nk),pow_full(nk))

     !Set the redshift
     z=0.

     !Assigns the cosmological model
     icosmo=4
     CALL assign_cosmology(icosmo,cosm,verbose)

     !Normalises power spectrum (via sigma_8) and fills sigma(R) look-up tables
     !CALL initialise_cosmology(verbose,cosm)
     CALL print_cosmology(cosm)

     !Initiliasation for the halo-model calcualtion
     CALL init_halomod(ihm,mmin,mmax,z,hmod,cosm,verbose)  

     !DMONLY
     j1=-1
     j2=-1
     outfile='variations/DMONLY.dat'
     CALL calculate_halomod(j1,j2,k,nk,z,pow_lin,pow_2h,pow_1h,pow_full,hmod,cosm,verbose)
     CALL write_power(k,pow_lin,pow_2h,pow_1h,pow_full,nk,outfile,verbose)

     !Prevents warning
     ilog=.FALSE.

     !Number of parameters
     npa=6 

     !Loop over parameters     
     DO ipa=1,npa

        !DO NOT DELETE - needs to be here to restore default cosmology on each loop
        !Reassigns the cosmological model
        CALL assign_cosmology(icosmo,cosm,verbose)

        !DO NOT DELETE - needs to be here to restore default cosmology on each loop
        !Normalises power spectrum (via sigma_8) and fills sigma(R) look-up tables
        !CALL initialise_cosmology(verbose,cosm)
        CALL print_cosmology(cosm) 
        
        !Set maximum and minimum parameter values and linear or log range
        IF(ipa==1) THEN
           !alpha - virial temperature pre factor
           param_min=0.1
           param_max=1.1
           ilog=.FALSE.
        ELSE IF(ipa==2) THEN
           !epsilon - concentration change due to gas
           param_min=0.5
           param_max=2.
           ilog=.TRUE.
        ELSE IF(ipa==3) THEN
           !Gamma - KS polytropic index
           param_min=1.15
           param_max=1.25
           ilog=.FALSE.
        ELSE IF(ipa==4) THEN
           !M0 - bound gas transition in Msun/h
           param_min=1e13
           param_max=1e15
           ilog=.TRUE.
        ELSE IF(ipa==5) THEN
           !A* - Stellar mass fraction
           param_min=0.01
           param_max=0.03
           ilog=.FALSE.
        ELSE IF(ipa==6) THEN
           !WHIM temperature in K
           param_min=1e5
           param_max=1e7
           ilog=.TRUE.
        END IF

        !Loop over parameter values
        DO i=1,m

           !Set the parameter value that is being varied
           IF(ilog) THEN
              param=progression_log(param_min,param_max,i,m)
           ELSE
              param=progression(param_min,param_max,i,m)
           END IF

           IF(ipa==1) cosm%alpha=param
           IF(ipa==2) cosm%eps=param
           IF(ipa==3) cosm%Gamma=param
           IF(ipa==4) cosm%M0=param
           IF(ipa==5) cosm%Astar=param
           IF(ipa==6) cosm%whim=param

           !DO NOT DELETE - needs to be here to restore default cosmology on each loop
           !Initiliasation for the halo-model calcualtion
           CALL init_cosmology(cosm)
           CALL init_halomod(ihm,mmin,mmax,z,hmod,cosm,verbose) 

           !DO NOT DELETE THIS
           !It is only used to print values to the screen later
           !For example, mass, which is inconvenient if written out in full
           IF(ilog) param=log10(param)

           !Write out halo matter and electron-pressure profile information
           !All the string crap is in the loop for a reason
           DO j=10,16
              base='variations/profile_mass_'
              ext='_param_'
              base=number_file(base,j,ext)
              mid='_value_'
              ext='.dat'
              outfile=number_file2(base,ipa,mid,i,ext)
              mass=10.**j
              CALL write_halo_profiles(mass,z,hmod,cosm,outfile)
           END DO

           !Write out halo mass fraction information
           base='variations/mass_fractions_param_'
           outfile=number_file2(base,ipa,mid,i,ext)
           CALL write_mass_fractions(cosm,outfile)

           !File base and extension
           base='variations/power_param_'
           mid='_value_'

           !Do the calculation
           DO j=1,9

              IF(j==1) THEN
                 !matter - matter
                 j1=0
                 j2=0
                 ext='_dd.dat'
              ELSE IF(j==2) THEN
                 !matter - electron pressure
                 j1=0
                 j2=6
                 ext='_dp.dat'
              ELSE IF(j==3) THEN
                 !electron pressure - electron pressure
                 j1=6
                 j2=6
                 ext='_pp.dat'
              ELSE IF(j==4) THEN
                 !matter-CDM
                 j1=0
                 j2=1
                 ext='_dc.dat'
              ELSE IF(j==5) THEN
                 !CDM-CDM
                 j1=1
                 j2=1
                 ext='_cc.dat'
              ELSE IF(j==6) THEN
                 !matter-gas
                 j1=0
                 j2=2
                 ext='_dg.dat'
              ELSE IF(j==7) THEN
                 !gas-gas
                 j1=2
                 j2=2
                 ext='_gg.dat'
              ELSE IF(j==8) THEN
                 !Matter-star
                 j1=0
                 j2=3
                 ext='_ds.dat'
              ELSE IF(j==9) THEN
                 !Star-star
                 j1=3
                 j2=3
                 ext='_ss.dat'
              END IF

              !Set output file
              outfile=number_file2(base,ipa,mid,i,ext)

              !Write progress to screen
              WRITE(*,fmt='(4I5,F14.7,A50)') ipa, i, j1, j2, param, TRIM(outfile)

              !Do the halo-model calculation and write to file
              CALL calculate_halomod(j1,j2,k,nk,z,pow_lin,pow_2h,pow_1h,pow_full,hmod,cosm,.FALSE.)
              CALL write_power(k,pow_lin,pow_2h,pow_1h,pow_full,nk,outfile,.FALSE.)

           END DO

        END DO

     END DO

  ELSE IF(imode==20) THEN

     !Set the cosmology
     icosmo=3
     CALL assign_cosmology(icosmo,cosm,verbose)
     CALL print_cosmology(cosm)

     !Set the halo model
     z=0.
     ihm=4
     CALL init_halomod(ihm,mmin,mmax,z,hmod,cosm,verbose)

     !Make the Figure
     CALL YinZhe_Fig1(z,hmod,cosm)

  ELSE IF(imode==21) THEN

     !Stuff for diagnosing problems with the window function integrand
     outfile='winint/integrand.dat'
     irho=11
     rv=1.
     c=4.
     rs=rv/c
     p1=1.18
     p2=0.
     rmin=0.
     rmax=rv
     CALL winint_diagnostics(rmin,rmax,rv,rs,p1,p2,irho,outfile)

  ELSE IF(imode==22) THEN

     !k range
     kmin=1e-2
     kmax=1e3
     nk=512
     CALL fill_array(log(kmin),log(kmax),k,nk)
     k=exp(k)

     !Halo parameters
     rv=1.
     c=4.
     rs=rv/c
     p1=1.5
     p2=0.
     irho=11
     rmin=0.
     rmax=rv

     CALL winint_speed_tests(k,nk,rmin,rmax,rv,rs,p1,p2,irho)

  ELSE IF(imode==23) THEN

     ! Set number of k points and k range (log spaced)
     nk=128
     kmin=1e-3
     kmax=1e1
     CALL fill_array(log(kmin),log(kmax),k,nk)
     k=exp(k)
     ALLOCATE(pow_lin(nk),pow_2h(nk),pow_1h(nk),pow_full(nk))

     ! Set the redshift
     z=0.

     ! Directory for output
     dir='Mead2017'

     ! Set the halo model
     ihm=12

     ! Loop over cosmologies
     DO i=1,9

        IF(i==1) THEN
           ! LCDM
           icosmo=1
           outfile='LCDM'
        ELSE IF(i==2) THEN
           ! OCDM
           icosmo=5
           outfile='OCDM'
        ELSE IF(i==3) THEN
           ! EdS
           icosmo=15
           outfile='SCDM'
        ELSE IF(i==4) THEN
           ! w = -0.7
           icosmo=16
           outfile='w-0.7'
        ELSE IF(i==5) THEN
           ! w = -1.3
           icosmo=17
           outfile='w-1.3'
        ELSE IF(i==6) THEN
           ! wa = 0.5
           icosmo=18
           outfile='wa0.5'
        ELSE IF(i==7) THEN
           ! wa = -0.5
           icosmo=19
           outfile='wa-0.5'
        ELSE IF(i==8) THEN
           ! w = -0.7; wa = -1.5
           icosmo=20
           outfile='w0-0.7wa-1.5'
        ELSE IF(i==9) THEN
           ! w = -1.3; wa = 0.5
           icosmo=21
           outfile='w0-1.3wa0.5'
        END IF

        ! Assigns the cosmological model
        CALL assign_cosmology(icosmo,cosm,verbose)

        ! Normalises power spectrum (via sigma_8) and fills sigma(R) look-up tables
        CALL print_cosmology(cosm)

        ! Initiliasation for the halomodel calcualtion
        CALL init_halomod(ihm,mmin,mmax,z,hmod,cosm,verbose)

        ! Do the halo-model calculation
        CALL calculate_halomod(-1,-1,k,nk,z,pow_lin,pow_2h,pow_1h,pow_full,hmod,cosm,verbose)

        ! Write out the results
        outfile=TRIM(dir)//'/power_'//TRIM(outfile)//'.dat'
        CALL write_power(k,pow_lin,pow_2h,pow_1h,pow_full,nk,outfile,verbose)

     END DO
     
  ELSE

     STOP 'HMx_DRIVER: Error, you have specified the mode incorrectly'

  END IF

CONTAINS

  SUBROUTINE get_k_values(infile,k,nk)

    IMPLICIT NONE
    CHARACTER(len=256), INTENT(IN) :: infile
    REAL, ALLOCATABLE, INTENT(OUT) :: k(:)
    INTEGER, INTENT(OUT) :: nk

    !Get the number of k values
    nk=file_length(infile)

    !Allocate the array in k
    ALLOCATE(k(nk))

    !Read in the k values
    OPEN(7,file=infile)
    DO i=1,nk
       READ(7,*) k(i)
    END DO
    CLOSE(7)
    
  END SUBROUTINE get_k_values

  SUBROUTINE YinZhe_Fig1(z,hmod,cosm)

    IMPLICIT NONE
    REAL, INTENT(IN) :: z
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: i
    REAL :: r, rs, rv, c, Mh, rh, r500c, m500c

    REAL, PARAMETER :: M=1e15 !Halo virial? mass [Msun]
    REAL, PARAMETER :: rmin=1e-3!Minimum radius [Mpc]
    REAL, PARAMETER :: rmax=8 !Maximum radius [Mpc] 
    INTEGER, PARAMETER :: nr=512 !Number of points in radius

    IF(hmod%has_mass_conversions .EQV. .FALSE.) CALL convert_mass_definitions(z,hmod,cosm)

    Mh=M*cosm%h !This virial mass is now [Msun/h]

    rv=exp(find(log(Mh),hmod%log_m,log(hmod%rv),hmod%n,3,3,2)) ![Mpc/h]
    c=find(log(Mh),hmod%log_m,hmod%c,hmod%n,3,3,2)
    rs=rv/c ![Mpc/h]

    m500c=exp(find(log(Mh),hmod%log_m,log(hmod%m500c),hmod%n,3,3,2)) ![Mpc/h]
    r500c=exp(find(log(Mh),hmod%log_m,log(hmod%r500c),hmod%n,3,3,2)) ![Mpc/h]

    WRITE(*,*) 'YINZHE_FIG1: Making data for this figure'
    WRITE(*,*) 'YINZHE_FIG1: Redshift:', z
    WRITE(*,*) 'YINZHE_FIG1: Virial radius [Mpc]:', rv/cosm%h
    WRITE(*,*) 'YINZHE_FIG1: Virial radius [Mpc/h]:', rv
    WRITE(*,*) 'YINZHE_FIG1: r_500,c [Mpc]:', r500c/cosm%h
    WRITE(*,*) 'YINZHE_FIG1: r_500,c [Mpc/h]:', r500c
    WRITE(*,*) 'YINZHE_FIG1: r_500,c / r_v:', r500c/rv
    WRITE(*,*) 'YINZHE_FIG1: Virial halo mass [log10 Msun]:', log10(M)
    WRITE(*,*) 'YINZHE_FIG1: Virial halo mass [log10 Msun/h]:', log10(Mh)    
    WRITE(*,*) 'YINZHE_FIG1: M_500,c [log10 Msun]:', log10(M500c/cosm%h)
    WRITE(*,*) 'YINZHE_FIG1: M_500,c [log10 Msun/h]:', log10(M500c)
    WRITE(*,*) 'YINZHE_FIG1: M_500,c / M_v:', M500c/Mh
    WRITE(*,*) 'YINZHE_FIG1: Halo concentraiton:', c

    OPEN(7,file='diagnostics/YinZhe_Fig1.dat')
    DO i=1,nr
       r=progression(rmin,rmax,i,nr) !Radius [Mpc]
       rh=r*cosm%h !Convert [Mpc/h]
       WRITE(7,*) r, UPP(.TRUE.,rh,z,Mh,rv,rs,hmod,cosm)*r**2, win_electron_pressure(.TRUE.,1,rh,z,Mh,rv,rs,hmod,cosm)*r**2
    END DO
    CLOSE(7)

    WRITE(*,*) 'YINZHE_FIG1: Done'
    WRITE(*,*)
    
  END SUBROUTINE YinZhe_Fig1

   SUBROUTINE write_power(k,pow_lin,pow_2h,pow_1h,pow,nk,output,verbose)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: output
    INTEGER, INTENT(IN) :: nk
    REAL, INTENT(IN) :: k(nk), pow_lin(nk), pow_2h(nk), pow_1h(nk), pow(nk)
    LOGICAL, INTENT(IN) :: verbose
    INTEGER :: i

    IF(verbose) WRITE(*,*) 'WRITE_POWER: Writing power to ', TRIM(output)

    !Loop over k values
    !Fill the tables with one- and two-halo terms as well as total
    OPEN(7,file=output)
    DO i=1,nk       
       WRITE(7,fmt='(5ES20.10)') k(i), pow_lin(i), pow_2h(i), pow_1h(i), pow(i)
    END DO
    CLOSE(7)

    IF(verbose) THEN
       WRITE(*,*) 'WRITE_POWER: Done'
       WRITE(*,*)
    END IF

  END SUBROUTINE write_power

  SUBROUTINE write_power_a_multiple(k,a,pow_lin,pow_2h,pow_1h,pow_full,nk,na,base,verbose)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: base
    INTEGER, INTENT(IN) :: nk, na
    REAL, INTENT(IN) :: k(nk), a(na), pow_lin(nk,na), pow_2h(nk,na), pow_1h(nk,na), pow_full(nk,na)
    LOGICAL, INTENT(IN) :: verbose
    REAL :: pow(nk,na)
    INTEGER :: i
    CHARACTER(len=512) :: output
    LOGICAL :: verbose2

    DO i=1,4
       IF(i==1) THEN
          output=TRIM(base)//'_linear.dat'
          pow=pow_lin
       ELSE IF(i==2) THEN
          output=TRIM(base)//'_2halo.dat'
          pow=pow_2h
       ELSE IF(i==3) THEN
          output=TRIM(base)//'_1halo.dat'
          pow=pow_1h
       ELSE IF(i==4) THEN
          output=TRIM(base)//'_full.dat'
          pow=pow_full
       ELSE
          STOP 'WRITE_POWER_A_MULTIPLE: Error, something went FUBAR'
       END IF
       IF(i==1) THEN
          verbose2=verbose
       ELSE
          verbose2=.FALSE.
       END IF
       CALL write_power_a(k,a,pow,nk,na,output,verbose2)
    END DO

  END SUBROUTINE write_power_a_multiple

  SUBROUTINE write_power_a(k,a,pow,nk,na,output,verbose)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: output
    INTEGER, INTENT(IN) :: nk, na
    REAL, INTENT(IN) :: k(nk), a(na), pow(nk,na)
    LOGICAL, INTENT(IN) :: verbose
    INTEGER :: i, j

    !Print to screen
    IF(verbose) THEN
       WRITE(*,*) 'WRITE_POWER_A: The first entry of the file is hashes - #####'
       WRITE(*,*) 'WRITE_POWER_A: The remainder of the first row are the scale factors - a'
       WRITE(*,*) 'WRITE_POWER_A: The remainder of the first column are the wave numbers - k'
       WRITE(*,*) 'WRITE_POWER_A: Each row then gives the power at that k and a'
       WRITE(*,*) 'WRITE_POWER_A: Output:', TRIM(output)
    END IF

    !Write out data to files
    OPEN(7,file=output)
    DO i=0,nk
       IF(i==0) THEN
          WRITE(7,fmt='(A20,40F20.10)') '#####', (a(j), j=1,na)
       ELSE
          WRITE(7,fmt='(F20.10,40E20.10)') k(i), (pow(i,j), j=1,na)
       END IF
    END DO
    CLOSE(7)

    !Print to screen
    IF(verbose) THEN
       WRITE(*,*) 'WRITE_POWER_A: Done'
       WRITE(*,*)
    END IF

  END SUBROUTINE write_power_a

  SUBROUTINE write_distances(cosm)

    ! Write file of z vs. r(z)
    IMPLICIT NONE
    TYPE(cosmology), INTENT(INOUT) :: cosm
    CHARACTER(len=256) :: output
    INTEGER :: i
    REAL :: z

    ! Now write the results of r(z) calculation
    output='projection/distance.dat'
    WRITE(*,*) 'WRITE_DISTANCE: Writing r(a): ', TRIM(output)
    OPEN(7,file=output)
    DO i=1,cosm%n_r
       z=redshift_a(cosm%a_r(i))
       WRITE(7,*) z, cosm%r(i), f_k(cosm%r(i),cosm)
    END DO
    CLOSE(7)
    WRITE(*,*) 'WRITE_DISTANCE: Done'
    WRITE(*,*)

  END SUBROUTINE write_distances

  SUBROUTINE random_cosmology(cosm)

    !Generate some random cosmological parameter
    USE random_numbers
    IMPLICIT NONE
    TYPE(cosmology), INTENT(INOUT) :: cosm

    REAL, PARAMETER :: Om_m_min=0.1
    REAL, PARAMETER :: Om_m_max=1.

    REAL, PARAMETER :: Om_b_on_Om_m_min=0.05
    REAL, PARAMETER :: Om_b_on_Om_m_max=0.5

    REAL, PARAMETER :: n_min=0.5
    REAL, PARAMETER :: n_max=1.5

    REAL, PARAMETER :: h_min=0.4
    REAL, PARAMETER :: h_max=1.2

    REAL, PARAMETER :: w_min=-1.5
    REAL, PARAMETER :: w_max=-0.5

    REAL, PARAMETER :: wa_min=0.
    REAL, PARAMETER :: wa_max=0.

    REAL, PARAMETER :: sig8_min=0.2
    REAL, PARAMETER :: sig8_max=1.5

    cosm%Om_m=random_uniform(Om_m_min,Om_m_max)

    !Enforce flatness
    cosm%Om_v=1.-cosm%Om_m

    cosm%Om_b=cosm%Om_m*random_uniform(Om_b_on_Om_m_min,Om_b_on_Om_m_max)

    cosm%n=random_uniform(n_min,n_max)

    cosm%h=random_uniform(h_min,h_max)

    cosm%w=random_uniform(w_min,w_max)

    cosm%wa=random_uniform(wa_min,wa_max)

    cosm%sig8=random_uniform(sig8_min,sig8_max)

  END SUBROUTINE random_cosmology

  SUBROUTINE random_baryon_parameters(cosm)

    IMPLICIT NONE
    TYPE(cosmology), INTENT(INOUT) :: cosm

    REAL, PARAMETER :: alpha_min=0.05
    REAL, PARAMETER :: alpha_max=2.5

    !REAL, PARAMETER :: eps_min=10**(-1.5)
    !REAL, PARAMETER :: eps_max=10**(1.5)
    REAL, PARAMETER :: eps_min=0.5
    REAL, PARAMETER :: eps_max=2.0
    !REAL, PARAMETER :: eps_min=1.
    !REAL, PARAMETER :: eps_max=1.

    REAL, PARAMETER :: Gamma_min=1.05
    REAL, PARAMETER :: Gamma_max=2.00
    !REAL, PARAMETER :: Gamma_min=1.15
    !REAL, PARAMETER :: Gamma_max=1.20

    REAL, PARAMETER :: M0_min=10**(12.)
    REAL, PARAMETER :: M0_max=10**(15.)

    REAL, PARAMETER :: Astar_min=0.002
    REAL, PARAMETER :: Astar_max=0.2

    REAL, PARAMETER :: whim_min=10**(5.)
    REAL, PARAMETER :: whim_max=10**(7.)
    
    cosm%alpha=random_uniform(alpha_min,alpha_max)

    cosm%eps=random_uniform(log(eps_min),log(eps_max))
    cosm%eps=exp(cosm%eps)

    cosm%Gamma=random_uniform(Gamma_min,Gamma_max)

    cosm%M0=random_uniform(log(M0_min),log(M0_max))
    cosm%M0=exp(cosm%M0)
    
    cosm%Astar=random_uniform(Astar_min,Astar_max)

    cosm%whim=random_uniform(log(whim_min),log(whim_max))
    cosm%whim=exp(cosm%whim)
    
  END SUBROUTINE random_baryon_parameters

END PROGRAM HMx_driver
