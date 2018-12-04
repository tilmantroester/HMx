PROGRAM HMx_fitting
  
  USE HMx
  USE cosmology_functions
  USE cosmic_emu_stuff
  USE string_operations
  USE random_numbers
  USE special_functions

  IMPLICIT NONE
  INTEGER :: im
  INTEGER :: ncos, nf, nz, nk
  REAL, ALLOCATABLE :: k(:), z(:), pow_bm(:,:,:,:,:), pow_hm(:,:,:,:,:), weight(:,:,:,:,:)
  INTEGER, ALLOCATABLE :: fields(:)
  TYPE(halomod), ALLOCATABLE :: hmod(:)
  TYPE(cosmology), ALLOCATABLE :: cosm(:)
  REAL :: kmin, kmax
  CHARACTER(len=256) :: name, base, mode, zin, outbase, outfile, nchain, maxtime, accuracy
  REAL, ALLOCATABLE :: pow_sim(:), k_sim(:)
  REAL, ALLOCATABLE :: p_min(:), p_max(:), p_bst(:), p_rge(:), p_new(:), p_old(:), p_ori(:)
  REAL, ALLOCATABLE :: p_lim(:), p_lam(:)
  CHARACTER(len=256), ALLOCATABLE :: p_nme(:)
  LOGICAL, ALLOCATABLE :: p_log(:), p_set(:), p_cov(:)
  REAL :: delta, fom_bst, fom_new, fom_old, fom_ori
  REAL :: t1, t2, tmax
  LOGICAL :: accept
  INTEGER :: icosmo, ihm, i_bst, np, ip(2)
  INTEGER :: i, j, l, j1, j2, n
  INTEGER :: i_bet, i_wor, i_acc, i_fai, i_tot
  LOGICAL :: verbose2
  INTEGER :: out

  ! Fitting parameters
  INTEGER, PARAMETER :: param_alpha=1
  INTEGER, PARAMETER :: param_eps=2
  INTEGER, PARAMETER :: param_gamma=3
  INTEGER, PARAMETER :: param_M0=4
  INTEGER, PARAMETER :: param_Astar=5
  INTEGER, PARAMETER :: param_Twhim=6
  INTEGER, PARAMETER :: param_cstar=7
  INTEGER, PARAMETER :: param_fcold=8
  INTEGER, PARAMETER :: param_mstar=9
  INTEGER, PARAMETER :: param_sstar=10
  INTEGER, PARAMETER :: param_alphap=11
  INTEGER, PARAMETER :: param_Gammap=12
  INTEGER, PARAMETER :: param_cstarp=13
  INTEGER, PARAMETER :: param_fhot=14
  INTEGER, PARAMETER :: param_alphaz=15
  INTEGER, PARAMETER :: param_Gammaz=16
  INTEGER, PARAMETER :: param_M0z=17
  INTEGER, PARAMETER :: param_Astarz=18
  INTEGER, PARAMETER :: param_Twhimz=19
  INTEGER, PARAMETER :: param_eta=20
  INTEGER, PARAMETER :: param_HMcode_Dv0=21
  INTEGER, PARAMETER :: param_HMcode_Dvp=22
  INTEGER, PARAMETER :: param_HMcode_dc0=23
  INTEGER, PARAMETER :: param_HMcode_dcp=24
  INTEGER, PARAMETER :: param_HMcode_eta0=25
  INTEGER, PARAMETER :: param_HMcode_eta1=26
  INTEGER, PARAMETER :: param_HMcode_f0=27
  INTEGER, PARAMETER :: param_HMcode_fp=28
  INTEGER, PARAMETER :: param_HMcode_kstar=29
  INTEGER, PARAMETER :: param_HMcode_As=30
  INTEGER, PARAMETER :: param_HMcode_alpha0=31
  INTEGER, PARAMETER :: param_HMcode_alpha1=32
  INTEGER, PARAMETER :: param_n=32

  ! Halo model calculation parameters
  REAL, PARAMETER :: mmin=mmin_HMx ! Minimum halo mass for the calculation
  REAL, PARAMETER :: mmax=mmax_HMx ! Maximum halo mass for the calculation

  INTEGER, PARAMETER :: m=huge(m)            ! Re-evaluate range every 'm' points
  INTEGER, PARAMETER :: seed=0               ! Random-number seed
  LOGICAL, PARAMETER :: random_start=.FALSE. ! Start from a random point within the prior range
  LOGICAL, PARAMETER :: mcmc=.TRUE.          ! Accept worse figure of merit with some probability
  INTEGER, PARAMETER :: computer=1           ! Which computer are you on?
  REAL, PARAMETER :: tmax_default=10000.*60  ! Default maximum time for run, should not be huge
  REAL, PARAMETER :: delta_default=1e-3      ! Default accuracy

  ! k cuts for the BAHAMAS power spectra
  REAL, PARAMETER :: kmin_BAHAMAS=0.15
  REAL, PARAMETER :: kmax_BAHAMAS_z0p0=10.
  REAL, PARAMETER :: kmax_BAHAMAS_z0p5=4.
  REAL, PARAMETER :: kmax_BAHAMAS_z1p0=2.
  REAL, PARAMETER :: kmax_BAHAMAS_z2p0=1.

  ! Read in starting option
  CALL get_command_argument(1,mode)
  IF(mode=='') THEN
     !STOP 'HMx_FITTING: Error, please specify mode'
     im=-1
  ELSE
     READ(mode,*) im
  END IF

  ! Decide what to do
  IF(im==-1) THEN
     WRITE(*,*)
     WRITE(*,*) 'HMx_FITTING: Choose what to do'
     WRITE(*,*) '=============================='
     WRITE(*,*) ' 1 - Mira Titan nodes'
     WRITE(*,*) ' 2 - FrankenEmu nodes'
     WRITE(*,*) ' 3 - Random Mira Titan'
     WRITE(*,*) ' 4 - Random FrankenEmu'
     WRITE(*,*) '11 - Hydro: fixed z; everything'
     WRITE(*,*) '12 - Hydro: fixed z; gas'
     WRITE(*,*) '13 - Hydro: fixed z; stars'
     WRITE(*,*) '14 - Hydro: fixed z; gas and stars'
     WRITE(*,*) '15 - Hydro: fixed z; matter'
     WRITE(*,*) '16 - Hydro: fixed z; everything but pressure-pressure'
     WRITE(*,*) '17 - Hydro: all z; everything but pressure-pressure'
     WRITE(*,*) '18 - Hydro: all z; only matter-matter and matter-pressure'
     WRITE(*,*) '19 - Hydro: all z; stars'
     WRITE(*,*) '20 - Hydro: all z; finessed version of 17: DOES NOT WORK'
     WRITE(*,*) '21 - Hydro: all z; no stars or pressure-pressure'
     READ(*,*) im
     WRITE(*,*)
  END IF

  ! Read in chain length
  CALL get_command_argument(2,nchain)
  IF(nchain=='') THEN
     WRITE(*,*) 'HMx_FITTING: Specify points in fitting chain'
     READ(*,*) n
     WRITE(*,*)
  ELSE
     READ(nchain,*) n
  END IF

  ! Read in maximum time
  CALL get_command_argument(3,maxtime)
  IF(maxtime=='') THEN
     tmax=tmax_default
  ELSE
     READ(maxtime,*) tmax
  END IF

  ! Accuracy
  CALL get_command_argument(4,accuracy)
  IF(accuracy=='') THEN
     delta=delta_default
  ELSE
     READ(accuracy,*) delta
  END IF

  ! Read in outfile
  CALL get_command_argument(5,outbase)
  IF(outbase=='') outbase='fitting/test'

  ! Read in BAHAMAS simulation name
  CALL get_command_argument(6,name)

  ! Read in BAHAMAS simulation redshift if doing fixed z
  CALL get_command_argument(7,zin)
  
  ! Set the random-number generator
  CALL RNG_set(seed)

  WRITE(*,*) 'HMx_FITTING: Fitting routine'
  WRITE(*,*) 'HMx_FITTING: Mode:', im
  WRITE(*,*) 'HMx_FITTING: Number of points in chain:', n
  WRITE(*,*) 'HMx_FITTING: Maximum run time [mins]:', tmax
  WRITE(*,*) 'HMx_FITTING: Maximum run time [hours]:', tmax/60.
  WRITE(*,*) 'HMx_FITTING: Accuracy:', delta
  WRITE(*,*) 'HMx_FITTING: Output file base: ', TRIM(outbase)
  WRITE(*,*) 'HMx_FITTING: Random number seed:', seed
  WRITE(*,*) 'HMx_FITTING: Random start:', random_start
  WRITE(*,*) 'HMx_FITTING: MCMC mode:', mcmc
  WRITE(*,*) 'HMx_FITTING: Re-evaluate covariance:', m
  WRITE(*,*)
  
  ! SET: number of cosmologies, fields and delta
  IF(im==1 .OR. im==2 .OR. im==3 .OR. im==4) THEN
     
     IF(im==1 .OR. im==3) THEN
        ncos=9 ! Number of Mita Titan nodes - only 10 have Omega_nu = 0. (ignore 1 because it is weird)
        nf=1
     ELSE IF(im==2 .OR. im==4) THEN
        ncos=37 ! Number of FrankenEmu nodes
        nf=1
     ELSE
        STOP 'HMx_FITTING: Error, something went wrong with setting fields'
     END IF
     
  ELSE IF(im>=11) THEN

     ! Set to the number of different cosmologies
     ncos=1 

     ! Set the number of different fields
     IF(im==11 .OR. im==16 .OR. im==17 .OR. im==18 .OR. im==20 .OR. im==21) THEN
        nf=5   
     ELSE IF(im==12 .OR. im==13 .OR. im==19) THEN
        nf=1
     ELSE IF(im==14) THEN
        nf=2
     ELSE IF(im==15) THEN
        nf=4
     ELSE
        STOP 'HMx_FITTING: Error, something went wrong with setting fields'
     END IF
     
  ELSE
     STOP 'HMx_FITTING: Error, something went wrong with setting fields'
  END IF

  ! Allocate arrays for cosmology and fields
  ALLOCATE(cosm(ncos),fields(nf))

  ! SET: redshifts and halo-profiles here     
  IF(im==1 .OR. im==2 .OR. im==3 .OR. im==4) THEN
     
     ! Mira Titan or FrankenEmu
     nz=4
     ALLOCATE(z(nz))
     z(1)=0.0
     z(2)=0.5
     z(3)=1.0
     z(4)=2.0
     fields(1)=field_dmonly
     
  ELSE IF(im>=11) THEN
     
     ! Hydro simulations
     IF(name=='') name='AGN_TUNED_nu0'
     
     IF(im==11 .OR. im==12 .OR. im==13 .OR. im==14 .OR. im==15 .OR. im==16) THEN
        nz=1
     ELSE IF(im==17 .OR. im==18 .OR. im==19 .OR. im==20 .OR. im==21) THEN
        nz=4
     ELSE
        STOP 'HMx_FITTING: Error, im not specified correctly'
     END IF
     
     ALLOCATE(z(nz))

     ! Set the redshifts
     IF(im==11 .OR. im==12 .OR. im==13 .OR. im==14 .OR. im==15 .OR. im==16) THEN        
        IF((zin)=='') THEN
           z(1)=0.0 
        ELSE
           READ(zin,*) z(1)
        END IF        
     ELSE IF(im==17 .OR. im==18 .OR. im==19 .OR. im==20 .OR. im==21) THEN
        z(1)=0.0
        z(2)=0.5
        z(3)=1.0
        z(4)=2.0
     ELSE
        STOP 'HMx_FITTING: Error, im not specified correctly'
     END IF

     ! Set the fields
     IF(im==11 .OR. im==16 .OR. im==17 .OR. im==18 .OR. im==20 .OR. im==21) THEN
        fields(1)=field_matter
        fields(2)=field_cdm
        fields(3)=field_gas
        fields(4)=field_star
        fields(5)=field_electron_pressure
     ELSE IF(im==12) THEN
        fields(1)=field_gas
     ELSE IF(im==13 .OR. im==19) THEN
        fields(1)=field_star
     ELSE IF(im==14) THEN
        fields(1)=field_gas
        fields(2)=field_star
     ELSE IF(im==15) THEN
        fields(1)=field_matter
        fields(2)=field_cdm
        fields(3)=field_gas
        fields(4)=field_star
     END IF
     
  ELSE
     STOP 'HMx_FITTING: Error, something went wrong'
  END IF

  ! Assign the cosmological models
  DO i=1,ncos
     
     IF(im==1) THEN
        icosmo=101+i ! Set Mira Titan node (note that we are skipping node 1)
     ELSE IF(im==2) THEN
        icosmo=200+i ! Set set FrankenEmu node
     ELSE IF(im==3) THEN
        icosmo=24    ! Random Mira Titan cosmology
     ELSE IF(im==4) THEN
        icosmo=25    ! Random FrankenEmu cosmology
     ELSE IF(im>=11) THEN
        icosmo=4     ! WMAP9
     ELSE
        STOP 'HMx_FITTING: Error, im not specified correctly'
     END IF
     
     CALL assign_cosmology(icosmo,cosm(i),verbose=.TRUE.)
     CALL init_cosmology(cosm(i))
     CALL print_cosmology(cosm(i))
     
  END DO

  ! SET: halo model here
  ALLOCATE(hmod(ncos))
  IF(im==1 .OR. im==2 .OR. im==3 .OR. im==4) THEN
     ihm=15 ! HMcode 2018
  ELSE IF(im>=11) THEN
     ihm=20 ! Standard halo-model response
  ELSE
     STOP 'HMx_FITTING: Error, im not specified correctly'
  END IF
     
  DO i=1,ncos
     CALL assign_halomod(ihm,hmod(i),verbose=.FALSE.)
  END DO

  ! Just print one halo model to screen to check
  CALL init_halomod(mmin,mmax,scale_factor_z(z(1)),hmod(1),cosm(1),verbose=.TRUE.)
  CALL print_halomod(hmod(1),cosm(1),verbose=.TRUE.)

  !! Read in the simulation power spectra for the models

  ! Loop over cosmologies and redshifts
  DO i=1,ncos
     DO j=1,nz

        ! Loop over fields
        DO j1=1,nf
           DO j2=j1,nf

              ! Read in power spectra
              IF(im==1 .OR. im==3) THEN
                 CALL read_Mira_Titan_power(k_sim,pow_sim,nk,z(j),cosm(i),rebin=.TRUE.)
              ELSE IF(im==2 .OR. im==4) THEN
                 CALL read_FrankenEmu_power(k_sim,pow_sim,nk,z(j),cosm(i),rebin=.TRUE.)
              ELSE IF(im>=11) THEN
                 ip(1)=fields(j1)
                 ip(2)=fields(j2)
                 CALL read_BAHAMAS_power(k_sim,pow_sim,nk,z(j),name,ip,cosm(i))!,kmin,kmax)
              ELSE
                 STOP 'HMx_FITTING: Error, something went wrong reading data'
              END IF

              ! Allocate big arrays for P(k,z,cosm)
              IF(.NOT. ALLOCATED(k))      ALLOCATE(k(nk))
              IF(.NOT. ALLOCATED(pow_bm)) ALLOCATE(pow_bm(ncos,nf,nf,nk,nz))
              k=k_sim
              pow_bm(i,j1,j2,:,j)=pow_sim
              
           END DO
        END DO

     END DO
  END DO

  ! Standard to give equal weight to everything
  ALLOCATE(weight(ncos,nf,nf,nk,nz))
  weight=1.

  ! Only weight to matter-pressure and matter-matter
  IF(im==18) THEN
     weight=0.
     weight(:,1,1,:,:)=1.
     weight(:,1,5,:,:)=1.
     weight(:,5,1,:,:)=1.
  END IF

  ! No weight to pressure-pressure
  IF(im==16 .OR. im==17 .OR. im==20) THEN
     weight(:,5,5,:,:)=0. 
  END IF

  ! No weight to pressure-pressure or to stars
  IF(im==21) THEN
     weight(:,4,:,:,:)=0.
     weight(:,:,4,:,:)=0.
     weight(:,5,5,:,:)=0.
  END IF 

  ! k range for multi-z
  IF(im>=11) THEN

     DO j=1,nz

        kmin=kmin_BAHAMAS
        IF(z(j)==0.0) THEN
           kmax=kmax_BAHAMAS_z0p0
        ELSE IF(z(j)==0.5) THEN
           kmax=kmax_BAHAMAS_z0p5
        ELSE IF(z(j)==1.0) THEN
           kmax=kmax_BAHAMAS_z1p0
        ELSE IF(z(j)==2.0) THEN
           kmax=kmax_BAHAMAS_z2p0
        ELSE
           STOP 'HMx_FITTING: Error, something went wrong setting z'
        END IF

        DO i=1,nk
           IF(k(i)<kmin .OR. k(i)>kmax) weight(:,:,:,i,j)=0.
        END DO

     END DO

  END IF

  !!

  ! Allocate arrays for halo-model power
  ALLOCATE(pow_hm(ncos,nf,nf,nk,nz))
  
  ! Allocate arrays for the parameters
  np=param_n
  ALLOCATE(p_nme(np),p_ori(np),p_min(np),p_max(np),p_log(np))
  ALLOCATE(p_lim(np),p_lam(np))
  ALLOCATE(p_bst(np),p_new(np),p_old(np))
  ALLOCATE(p_rge(np),p_set(np),p_cov(np))

  ! Set the initial parameter values, ranges, names, etc.
  CALL set_parameters(p_nme,p_ori,p_lim,p_lam,p_min,p_max,p_log,np)

  ! Initially all parameters default to not being varied
  p_set=.FALSE.

  ! Initially set the range for all parameter to be zero
  p_rge=0.

  ! Initially assume we calculate the step size for each parameter
  p_cov=.TRUE.

  IF(im==1 .OR. im==2 .OR. im==3 .OR. im==4) THEN

     ! HMcode
     p_set(param_HMcode_Dv0)=.TRUE.
     p_set(param_HMcode_Dvp)=.TRUE.
     p_set(param_HMcode_dc0)=.TRUE.
     p_set(param_HMcode_dcp)=.TRUE.
     p_set(param_HMcode_eta0)=.TRUE.
     p_set(param_HMcode_eta1)=.TRUE.
     p_set(param_HMcode_f0)=.TRUE.
     p_set(param_HMcode_fp)=.TRUE.
     p_set(param_HMcode_kstar)=.TRUE.
     p_set(param_HMcode_As)=.TRUE.
     p_set(param_HMcode_alpha0)=.TRUE.
     p_set(param_HMcode_alpha1)=.TRUE.

  ELSE IF(im>=11) THEN
     
     IF(im==11 .OR. im==16) THEN

        ! 11 - everything
        ! 16 - everything minus pressure-pressure
        p_set(param_alpha)=.TRUE.
        p_set(param_eps)=.TRUE.
        p_set(param_Gamma)=.TRUE.
        p_set(param_M0)=.TRUE.
        p_set(param_Astar)=.TRUE.
        p_set(param_Twhim)=.TRUE.
        p_set(param_cstar)=.TRUE.
        p_set(param_fcold)=.TRUE.
        p_set(param_Mstar)=.TRUE.
        p_set(param_sstar)=.TRUE.
        p_set(param_alphap)=.TRUE.
        p_set(param_Gammap)=.TRUE.
        p_set(param_cstarp)=.TRUE.
        p_set(param_fhot)=.TRUE.
        p_set(param_eta)=.TRUE.

     ELSE IF(im==17 .OR. im==18 .OR. im==20) THEN

        ! redshift dependent everything minus pressure-pressure
        p_set(param_alpha)=.TRUE.
        p_set(param_eps)=.TRUE.
        p_set(param_Gamma)=.TRUE.
        p_set(param_M0)=.TRUE.
        p_set(param_Astar)=.TRUE.
        p_set(param_Twhim)=.TRUE.
        p_set(param_cstar)=.TRUE.
        p_set(param_fcold)=.FALSE.
        p_set(param_Mstar)=.TRUE.
        p_set(param_sstar)=.TRUE.
        p_set(param_fhot)=.TRUE.
        IF(im .NE. 18) p_set(param_eta)=.TRUE.
        
        p_set(param_alphap)=.TRUE.
        p_set(param_Gammap)=.TRUE.
        p_set(param_cstarp)=.TRUE.        
        
        p_set(param_alphaz)=.TRUE.
        p_set(param_Gammaz)=.TRUE.
        p_set(param_M0z)=.TRUE.
        p_set(param_Astarz)=.TRUE.
        p_set(param_Twhimz)=.TRUE.

     ELSE IF(im==12) THEN

        ! gas
        p_set(param_eps)=.TRUE.
        p_set(param_Gamma)=.TRUE.
        p_set(param_M0)=.TRUE.
        p_set(param_Astar)=.TRUE.
        p_set(param_fcold)=.TRUE.
        p_set(param_Gammap)=.TRUE.

     ELSE IF(im==13) THEN

        ! stars
        p_set(param_Astar)=.TRUE.
        p_set(param_cstar)=.TRUE.
        p_set(param_cstarp)=.TRUE.
        p_set(param_Mstar)=.TRUE.
        p_set(param_sstar)=.TRUE.
        p_set(param_eta)=.TRUE.

     ELSE IF(im==19) THEN

        ! stars - z
        p_set(param_Astar)=.TRUE.
        p_set(param_cstar)=.TRUE.
        p_set(param_cstarp)=.TRUE.
        p_set(param_Mstar)=.TRUE.
        p_set(param_sstar)=.TRUE.
        p_set(param_eta)=.TRUE.
        p_set(param_Astarz)=.TRUE.

     ELSE IF(im==14 .OR. im==15) THEN

        ! 14 - gas and stars
        ! 15 - matter
        p_set(param_eps)=.TRUE.
        p_set(param_Gamma)=.TRUE.
        p_set(param_M0)=.TRUE.
        p_set(param_Astar)=.TRUE.
        p_set(param_cstar)=.TRUE.
        p_set(param_fcold)=.TRUE.
        p_set(param_Mstar)=.TRUE.
        p_set(param_sstar)=.TRUE.
        p_set(param_Gammap)=.TRUE.
        p_set(param_cstarp)=.TRUE.

     ELSE IF(im==21) THEN

        ! redshift dependent everything minus pressure-pressure
        p_set(param_alpha)=.TRUE.
        p_set(param_eps)=.TRUE.
        p_set(param_Gamma)=.TRUE.
        p_set(param_M0)=.TRUE.
        p_set(param_Twhim)=.TRUE.
        p_set(param_fhot)=.TRUE.
        
        p_set(param_alphap)=.TRUE.
        p_set(param_Gammap)=.TRUE.       
        
        p_set(param_alphaz)=.TRUE.
        p_set(param_Gammaz)=.TRUE.
        p_set(param_M0z)=.TRUE.
        p_set(param_Twhimz)=.TRUE.
        
     ELSE

        STOP 'HMx_FITTING: Something went wrong with setting parameters'
        
     END IF

     IF(im==20) CALL finesse_parameters(p_ori,n,name)

     IF(im==21) CALL fix_star_parameters(p_ori,n,name)
     
  ELSE

     STOP 'HMx_FITTING: Something went wrong with setting parameters'

  END IF

  ! Start from a random place within the prior
  IF(random_start) THEN
     DO i=1,np
        IF(p_set(i)) p_ori(i)=random_uniform(p_lim(i),p_lam(i))
     END DO
  END IF

  ! Take the log of those parameters that are explored in log
  DO i=1,np
     IF(p_log(i)) THEN
        p_ori(i)=log10(p_ori(i))
        p_min(i)=log10(p_min(i))
        p_max(i)=log10(p_max(i))
     END IF
  END DO

  ! Set the new parameters
  p_old=p_ori
  p_new=p_ori
  p_bst=p_ori

  ! Set the best figures-of-merit to some huge value
  fom_old=HUGE(fom_old)
  fom_new=HUGE(fom_new)
  fom_bst=HUGE(fom_bst)

  ! Loop over number of runs
  WRITE(*,*) 'HMx_FITTING: Starting fitting'
  WRITE(*,*) 'HMx_FITTING: Number of points in chain:', n
  WRITE(*,*)

  ! Set counting variables to zero
  i_bet=0
  i_wor=0
  i_acc=0
  i_fai=0
  i_tot=0

  ! Get the starting time
  CALL cpu_time(t1)
  CALL cpu_time(t2)

  ! Do the chain
  DO l=1,n+1

     IF(l==1 .OR. mod(l,m)==0) THEN
        IF(l==1) THEN
           verbose2=.TRUE.
        ELSE
           verbose2=.TRUE.
        END IF
        CALL set_parameter_ranges(delta,fields,nf,p_rge,p_set,p_cov,p_old,p_log,p_nme,np,k,nk,z,nz,pow_bm,weight,hmod,cosm,ncos,verbose2)
        IF(l==1) THEN
           outfile=trim(outbase)//'_chain.dat'
           OPEN(10,file=outfile)
        END IF

     END IF

     IF(l==1) THEN
        ! Do nothing
     ELSE IF(l==n+1 .OR. (t2-t1)>60.*tmax) THEN
        ! Set to best-fitting parameters on last go
        p_new=p_bst           
     ELSE
        ! Randomly jump parameters
        DO i=1,np
           IF(p_set(i)) THEN
              p_new(i)=random_Gaussian(p_old(i),p_rge(i))
           ELSE
              p_new(i)=p_old(i)
           END IF
           IF(p_new(i)<p_min(i)) p_new(i)=p_min(i)
           IF(p_new(i)>p_max(i)) p_new(i)=p_max(i)
        END DO
     END IF

     ! Calculate the figure-of-merit
     CALL fom_multiple(fields,nf,fom_new,p_new,p_log,np,k,nk,z,nz,pow_hm,pow_bm,weight,hmod,cosm,ncos)

!!$     WRITE(*,*) 'Go:', l
!!$     WRITE(*,*) 'FOM:', fom_new
!!$     DO i=1,np
!!$        IF(p_set(i)) WRITE(*,*) i, TRIM(p_nme(i)), p_new(i)
!!$     END DO
!!$     WRITE(*,*)
!!$     IF(l==2) STOP    

     ! Write original power spectra to disk
     IF(l==1) THEN

        ! Set the original figure-of-merit to the new figure-of-merit
        fom_old=fom_new
        fom_ori=fom_new
        fom_bst=fom_new

        ! Write out the original data
        base=trim(outbase)//'_orig_cos'
        CALL write_fitting_power(base,k,pow_hm,pow_bm,ncos,nf,nk,nz)

        accept=.TRUE.

        i_bst=1
        i_bet=i_bet+1
        i_tot=i_tot+1

     ELSE IF(l==n+1 .OR. (t2-t1)>60.*tmax) THEN

        WRITE(*,*)
        
        ! Output the best-fitting model
        base=trim(outbase)//'_best_cos'
        CALL write_fitting_power(base,k,pow_hm,pow_bm,ncos,nf,nk,nz)

        accept=.TRUE.
        EXIT

     ELSE

        i_tot=i_tot+1
        IF(fom_new < fom_bst) THEN
           ! If fit is the best then always accept...
           p_bst=p_new
           i_bst=l
           fom_bst=fom_new
           accept=.TRUE.
           i_bet=i_bet+1
        ELSE IF(fom_new <= fom_old) THEN
           ! ...also accept if fom is better than previous...
           accept=.TRUE.
           i_bet=i_bet+1
        ELSE IF(mcmc .AND. (fom_old/fom_new)**(1./delta) > random_uniform(0.,1.)) THEN
           ! ...otherwise accept poorer fit with some probability...
           accept=.TRUE.
           i_wor=i_wor+1
        ELSE
           ! ...otherwise, do not accept.
           accept=.FALSE.
           i_fai=i_fai+1
        END IF

     END IF

     IF(l .NE. n+1) WRITE(*,fmt='(I10,3F14.7,L3)') l, fom_bst, fom_old, fom_new, accept
     !IF(l .NE. n+1) WRITE(*,fmt='(I10,3F14.7,L3,F10.5)') l, fom_bst, fom_old, fom_new, accept, p_new(param_eta)

     IF(accept) THEN
        i_acc=i_acc+1
        p_old=p_new
        fom_old=fom_new
        WRITE(10,*) l, fom_old, (p_old(j), j=1,np)
     END IF

     CALL cpu_time(t2)

  END DO
  CLOSE(10)
  WRITE(*,*) 'HMx_FITTING: Done'
  WRITE(*,*)

  ! Write useful information to screen and file
  DO j=1,2

     IF(j==1) THEN
        out=6
     ELSE IF(j==2) THEN
        out=77
        outfile=trim(outbase)//'_params.dat'
        OPEN(out,file=outfile)
     ELSE
        STOP 'HMx_FITTING: Error, output fucked up badly'
     END IF

     WRITE(out,*) 'HMx_FITTING: Best location:', i_bst
     WRITE(out,*) 'HMx_FITTING: Total attempts:', i_tot
     WRITE(out,*) 'HMx_FITTING: Accepted steps:', i_acc
     WRITE(out,*) 'HMx_FITTING: Fraction accepted steps:', REAL(i_acc)/REAL(i_tot)
     WRITE(out,*) 'HMx_FITTING: Better steps:', i_bet
     WRITE(out,*) 'HMx_FITTING: Fraction better steps:', REAL(i_bet)/REAL(i_tot)
     WRITE(out,*) 'HMx_FITTING: Accepted worse steps:', i_wor
     WRITE(out,*) 'HMx_FITTING: Fraction accepted worse steps:', REAL(i_wor)/REAL(i_tot)
     WRITE(out,*) 'HMx_FITTING: Failed steps:', i_fai
     WRITE(out,*) 'HMx_FITTING: Fraction failed steps:', REAL(i_fai)/REAL(i_tot)
     WRITE(out,*) 'HMx_FITTING: Original figure-of-merit:', fom_ori
     WRITE(out,*) 'HMx_FITTING: Final figure-of-merit:', fom_new
     WRITE(out,*)

     WRITE(out,*) 'HMx_FITTING: Best-fitting parameters'
     WRITE(out,*) '================================================================================'
     WRITE(out,*) 'Parameter           Name      Original          Best       Minimum       Maximum'
     WRITE(out,*) '================================================================================'
     DO i=1,np
        IF(j==2 .OR. (j==1 .AND. p_set(i))) THEN
           IF(p_log(i)) THEN
              WRITE(out,fmt='(I10,A15,4F14.7)') i, trim(p_nme(i)), 10**p_ori(i), 10**p_bst(i), 10**p_min(i), 10**p_max(i)
              !WRITE(out,fmt='(I10,A15,A5,4F14.7)') i, trim(p_nme(i)), 'log', p_ori(i), p_bst(i), p_min(i), p_max(i)
           ELSE
              WRITE(out,fmt='(I10,A15,4F14.7)') i, trim(p_nme(i)), p_ori(i), p_bst(i), p_min(i), p_max(i)
           END IF
        END IF
     END DO
     WRITE(out,*) '====================================================================================='
     WRITE(out,*)

     IF(j==2) THEN
        CLOSE(out)
     END IF

  END DO

CONTAINS

  SUBROUTINE fom_multiple(fields,nf,fom,p,p_log,np,k,nk,z,nz,pow_hm,pow_sim,weight,hmod,cosm,n)

    USE special_functions

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: fields(nf) ! Field types
    INTEGER, INTENT(IN) :: nf ! Number of fields
    REAL, INTENT(OUT) :: fom ! Output figure of merit
    REAL, INTENT(IN) :: p(np) ! Array of varying parameters
    LOGICAL, INTENT(IN) :: p_log(np) ! Array of which parameters are explored in log
    INTEGER, INTENT(IN) :: np ! Number of varying parameters
    REAL, INTENT(IN) :: k(nk) ! Array of k values for comparison data
    INTEGER, INTENT(IN) :: nk ! Number of k values for comparison data
    REAL, INTENT(IN) :: z(nz) ! Array of z values for comparison data
    INTEGER, INTENT(IN) :: nz ! Number of z values
    REAL, INTENT(OUT) :: pow_hm(n,nf,nf,nk,nz) ! Output array for power as a function of k, z, cosmology
    REAL, INTENT(IN) :: pow_sim(n,nf,nf,nk,nz) ! Comparison power as a function of k, z, cosmology
    REAL, INTENT(IN) :: weight(n,nf,nf,nk,nz) ! Weight array
    TYPE(halomod), INTENT(INOUT) :: hmod(n) ! Array of halo models for each comparison
    TYPE(cosmology), INTENT(INOUT) :: cosm(n) ! Array of cosmological models for each comparison
    INTEGER, INTENT(IN) :: n ! Number of cosmological models being compared
    INTEGER :: i, j, ik
    REAL :: pow_li(n,nk,nz), pow_2h(n,nf,nf,nk,nz), pow_1h(n,nf,nf,nk,nz)
    REAL :: pexp(np), neff

    ! Set this counting output variable to zero
    fom=0.
    neff=0.

    ! Set this to zero too, for the banter
    pow_hm=0.

    ! Loop over cosmologies
    DO i=1,n

       ! SET: add parameters here

       ! Exponentiate those parameters that need it
       DO j=1,np
          IF(p_log(j)) THEN
             pexp(j)=10**p(j)
          ELSE
             pexp(j)=p(j)
          END IF
       END DO

       ! HMcode
       hmod(i)%Dv0=pexp(param_HMcode_Dv0)
       hmod(i)%Dv1=pexp(param_HMcode_Dvp)
       hmod(i)%dc0=pexp(param_HMcode_dc0)
       hmod(i)%dc1=pexp(param_HMcode_dcp)
       hmod(i)%eta0=pexp(param_HMcode_eta0)
       hmod(i)%eta1=pexp(param_HMcode_eta1)
       hmod(i)%f0=pexp(param_HMcode_f0)
       hmod(i)%f1=pexp(param_HMcode_fp)
       hmod(i)%ks=pexp(param_HMcode_kstar)
       hmod(i)%As=pexp(param_HMcode_As)
       hmod(i)%alp0=pexp(param_HMcode_alpha0)
       hmod(i)%alp1=pexp(param_HMcode_alpha1)

       ! Hydro
       hmod(i)%alpha=pexp(param_alpha)
       hmod(i)%eps=pexp(param_eps)
       hmod(i)%Gamma=1.+pexp(param_Gamma)
       hmod(i)%M0=10.**pexp(param_M0)
       hmod(i)%Astar=pexp(param_Astar)
       hmod(i)%Twhim=10.**pexp(param_Twhim)
       hmod(i)%cstar=pexp(param_cstar)
       hmod(i)%fcold=pexp(param_fcold)
       hmod(i)%Mstar=10.**pexp(param_Mstar)
       hmod(i)%sstar=pexp(param_sstar)
       hmod(i)%fhot=pexp(param_fhot)
       hmod(i)%eta=pexp(param_eta)
       hmod(i)%alphap=pexp(param_alphap)
       hmod(i)%Gammap=pexp(param_Gammap)
       hmod(i)%cstarp=pexp(param_cstarp)     
       hmod(i)%alphaz=pexp(param_alphaz)
       hmod(i)%Gammaz=pexp(param_Gammaz)
       hmod(i)%M0z=pexp(param_M0z)
       hmod(i)%Astarz=pexp(param_Astarz)
       hmod(i)%Twhimz=pexp(param_Twhimz)

       ! Loop over redshifts
       DO j=1,nz

          ! Initialise the halo-model calculation
          CALL init_halomod(mmin,mmax,scale_factor_z(z(j)),hmod(i),cosm(i),verbose=.FALSE.)
          CALL print_halomod(hmod(i),cosm(i),verbose=.FALSE.)

          ! Calculate the halo-model power spectrum
          CALL calculate_HMx_a(fields,nf,k,nk,pow_li(i,:,j),pow_2h(i,:,:,:,j),pow_1h(i,:,:,:,j),pow_hm(i,:,:,:,j),hmod(i),cosm(i),verbose=.FALSE.,response=.FALSE.)

          ! Calculate figure of merit and add to total
          DO j1=1,nf
             DO j2=j1,nf
                DO ik=1,nk
                   fom=fom+weight(i,j1,j2,ik,j)*(pow_hm(i,j1,j2,ik,j)/pow_sim(i,j1,j2,ik,j)-1.)**2
                   neff=neff+weight(i,j1,j2,ik,j)
                END DO
             END DO
          END DO

       END DO

    END DO

    ! Calculate the final figure-of-merit by dividing by the effective number of data points and sqrt
    fom=sqrt(fom/neff)

  END SUBROUTINE fom_multiple

  SUBROUTINE set_parameter_ranges(delta,fields,nf,p_rge,p_set,p_cov,p,p_log,p_nme,np,k,nk,z,nz,pow_sim,weight,hmod,cosm,n,verbose)

    IMPLICIT NONE
    REAL, INTENT(IN) :: delta
    INTEGER, INTENT(IN) :: fields(nf)
    INTEGER, INTENT(IN) :: nf
    REAL, INTENT(INOUT) :: p_rge(np)
    LOGICAL, INTENT(IN) :: p_set(np)
    LOGICAL, INTENT(IN) :: p_cov(np)
    REAL, INTENT(IN) :: p(np)
    LOGICAL, INTENT(IN) :: p_log(np)
    CHARACTER(len=*), INTENT(IN) :: p_nme(np)
    INTEGER, INTENT(IN) :: np
    REAL, INTENT(IN) :: k(nk)
    INTEGER, INTENT(IN) :: nk
    REAL, INTENT(IN) :: z(nz)
    INTEGER, INTENT(IN) :: nz
    REAL, INTENT(IN) :: pow_sim(n,nf,nf,nk,nz)
    REAL, INTENT(IN) :: weight(n,nf,nf,nz)
    TYPE(halomod), INTENT(INOUT) :: hmod(n)
    TYPE(cosmology), INTENT(INOUT) :: cosm(n)
    INTEGER, INTENT(IN) :: n
    LOGICAL, INTENT(IN) :: verbose
    INTEGER :: i, nv
    REAL :: fom_base, fom_diff, fom, df, p2(np), pow(n,nf,nf,nk,nz), dp
    
    REAL, PARAMETER :: eps=2.0    ! Tolerated error in fom difference when setting range
    REAL, PARAMETER :: deriv=1e-4 ! How much smaller is the derivative than delta

    ! Get the figure of merit for the base set of parameters
    CALL fom_multiple(fields,nf,fom_base,p,p_log,np,k,nk,z,nz,pow,pow_sim,weight,hmod,cosm,n)

    ! Calculate an initial value for the paramter change
    dp=deriv*delta

    ! Count the number of parameters being varied
    nv=0
    DO i=1,np
       IF(p_set(i)) nv=nv+1
    END DO

    ! Initial parameter perturbation
    DO i=1,np
       IF(p_set(i)) THEN
          IF(p_cov(i)) THEN
             IF(p(i) .NE. 0.) THEN
                p_rge(i)=p(i)*dp
             ELSE
                p_rge(i)=dp
             END IF
          END IF
       ELSE
          p_rge(i)=0.
       END IF
    END DO

    ! Write to screen
    IF(verbose) THEN
       WRITE(*,*) 'SET_PARAMETER_RANGES: Setting parameter step sizes'
       WRITE(*,*) 'SET_PARAMETER_RANGES: Number of varied parameters:', nv
       WRITE(*,*) 'SET_PARAMETER_RANGES: Total number of parameters:', np
       WRITE(*,*) 'SET_PARAMETER_RANGES: Number of cosmologies:', n
       WRITE(*,*) 'SET_PARAMETER_RANGES: Number of fields:', nf
       WRITE(*,*) 'SET_PARAMETER_RANGES: Number of wavenumbers:', nk
       WRITE(*,*) 'SET_PARAMETER_RANGES: Number of redshifts:', nz
       WRITE(*,*) 'SET_PARAMETER_RANGES: Derivatives being calculated with:', dp
       WRITE(*,*) 'SET_PARAMETER_RANGES: Fixing sigma to give change in fom:', delta
       WRITE(*,*) 'SET_PARAMETER_RANGES: Baseline fom:', fom_base
       WRITE(*,*) '================================================================================'
       WRITE(*,*) 'Parameter           Name    Base value     New Value         Sigma        Change'
       WRITE(*,*) '================================================================================'
    END IF

    ! Loop over parameters
    DO i=1,np

       IF(p_set(i)) THEN

          ! Set the range of p to take the derivative over
          p2=p ! Reset all
          p2(i)=p(i)+p_rge(i) ! Perturb parameter i

          ! Get the figure of merit for the updated parameter
          CALL fom_multiple(fields,nf,fom,p2,p_log,np,k,nk,z,nz,pow,pow_sim,weight,hmod,cosm,n)

          ! Calculate the change in the figure of merit for this parameter
          df=fom-fom_base

          ! Check that perturbing the parameter actually changes the figure of merit
          IF(df==0.) THEN
             WRITE(*,*) 'SET_PARAMETER_RANGES: Parameter:', i
             WRITE(*,*) 'SET_PARAMETER_RANGES: sigma:', p_rge(i)
             IF(p_log(i)) THEN
                WRITE(*,*) 'SET_PARAMETER_RANGES: Original value:', 10**p(i)
                WRITE(*,*) 'SET_PARAMETER_RANGES: Perturbed value:', 10**p2(i)
             ELSE
                WRITE(*,*) 'SET_PARAMETER_RANGES: Original value:', p(i)
                WRITE(*,*) 'SET_PARAMETER_RANGES: Perturbed value:', p2(i)
             END IF
             WRITE(*,*) 'SET_PARAMETER_RANGES: Original figure-of-merit:', fom_base
             WRITE(*,*) 'SET_PARAMETER_RANGES: Perturbed figure-of-merit:', fom
             WRITE(*,*) 'SET_PARAMETER_RANGES: Change in figure-of-merit:', df
             STOP 'SET_PARAMETER_RANGES: Error, changing parameter does not change power spectra'
          END IF

          IF(p_cov(i)) THEN

             ! Set sigma so that it gives a change of 'delta' in fom
             p_rge(i)=p_rge(i)*delta/df
             p2(i)=p(i)+p_rge(i)

          END IF

          ! Write parameters to screen
          IF(verbose) THEN
             IF(p_log(i)) THEN
                WRITE(*,fmt='(I10,A15,4F14.7)') i, trim(p_nme(i)), 10.**p(i), 10.**p2(i), 10.**p2(i)-10.**p(i), p_rge(i)/p(i)
             ELSE
                WRITE(*,fmt='(I10,A15,4F14.7)') i, trim(p_nme(i)), p(i), p2(i), p_rge(i), p_rge(i)/p(i)
             END IF
          END IF

       END IF

    END DO

    ! Write to screen
    IF(verbose) THEN
       WRITE(*,*) '================================================================================'
       WRITE(*,*) 'SET_PARAMETER_RANGES: Done initial setting'
       WRITE(*,*)
    END IF

    ! Write to screen
    IF(verbose) THEN
       WRITE(*,*) 'SET_PARAMETER_RANGES: Baseline figure of merit:', fom_base
       WRITE(*,*) '================================================================================'
       WRITE(*,*) 'Parameter           Name    Base value     New Value         Sigma         Ratio'
       WRITE(*,*) '================================================================================'
    END IF

    DO i=1,np

       IF(p_set(i)) THEN

          DO

             p2=p
             p2(i)=p(i)+p_rge(i)
             CALL fom_multiple(fields,nf,fom,p2,p_log,np,k,nk,z,nz,pow,pow_sim,weight,hmod,cosm,n)
             fom_diff=fom-fom_base
             
             ! Write parameters to screen
             IF(verbose) THEN
                IF(p_log(i)) THEN
                   WRITE(*,fmt='(I10,A15,4F14.7)') i, trim(p_nme(i)), 10**p(i), 10**p2(i), 10**p2(i)-10**p(i), fom_diff/delta
                ELSE
                   WRITE(*,fmt='(I10,A15,4F14.7)') i, trim(p_nme(i)), p(i), p2(i), p_rge(i), fom_diff/delta
                END IF
             END IF

             ! Make a new guess if the sigma was not correct
             IF(p_cov(i)) THEN
                IF(abs(fom_diff) > delta*eps) THEN
                   p_rge(i)=p_rge(i)/(fom_diff/delta)
                ELSE IF(abs(fom_diff) < delta/eps) THEN
                   p_rge(i)=p_rge(i)/(fom_diff/delta)
                ELSE
                   EXIT
                END IF
             ELSE
                EXIT
             END IF

          END DO

       END IF

    END DO

    ! Write to screen
    IF(verbose) THEN
       WRITE(*,*) '================================================================================'
       WRITE(*,*) 'Done check'
       WRITE(*,*)
    END IF

  END SUBROUTINE set_parameter_ranges

  SUBROUTINE read_simulation_power_spectrum(k,Pk,n,infile,kmin,kmax,cut_nyquist,subtract_shot,verbose)

    IMPLICIT NONE
    REAL, ALLOCATABLE, INTENT(OUT) :: k(:), Pk(:) ! Output simulation k and power
    INTEGER, INTENT(OUT) :: n ! Number of output k values
    CHARACTER(len=*), INTENT(IN) :: infile ! Input file location
    REAL, OPTIONAL, INTENT(IN) :: kmin, kmax ! Minimum and maximum k values to cut at
    LOGICAL, OPTIONAL, INTENT(IN) :: cut_nyquist ! Logical to cut Nyquist or not
    LOGICAL, OPTIONAL, INTENT(IN) :: subtract_shot ! Logical to subtract shot noise or not
    LOGICAL, OPTIONAL, INTENT(IN) :: verbose ! Logical verbose
    INTEGER :: i, i1, i2, m
    REAL :: shot, kbig
    LOGICAL :: lexist

    ! Deallocate arrays if they are already allocated
    IF(ALLOCATED(k))  DEALLOCATE(k)
    IF(ALLOCATED(Pk)) DEALLOCATE(Pk)

    ! Check file exists
    INQUIRE(file=infile,exist=lexist)
    IF(.NOT. lexist) THEN
       WRITE(*,*) 'READ_SIMULATION_POWER_SPECTRUM: File: ', trim(infile)
       STOP 'READ_SIMULATION_POWER_SPECTRUM: File does not exist'
    END IF

    ! Get file length and allocate arrays for output
    n=file_length(infile,verbose=.FALSE.)
    ALLOCATE(k(n),Pk(n))

    ! Write to screen
    IF(PRESENT(verbose)) THEN
       IF(verbose) THEN
          WRITE(*,*) 'READ_SIMULATION_POWER_SPECTRUM: Reading in data'
          WRITE(*,*) 'READ_SIMULATION_POWER_SPECTRUM: File: ', trim(infile)
          WRITE(*,*) 'READ_SIMULATION_POWER_SPECTRUM: nk:', nk
       END IF
    END IF

    ! Read in data from file
    OPEN(9,file=infile,status='old')
    DO i=1,n
       READ(9,*) k(i), Pk(i), shot
       IF(PRESENT(subtract_shot)) THEN
          IF(subtract_shot) Pk(i)=Pk(i)-shot
       END IF
    END DO
    CLOSE(9)

    IF(PRESENT(cut_nyquist)) THEN
       IF(cut_nyquist) THEN

          ! Find position in array of half-Nyquist
          kbig=k(n)
          DO i=1,n
             IF(k(i)>kbig/2.) EXIT
          END DO

          ! Cut arrays down to half-Nyquist
          CALL amputate(k,n,i)
          CALL amputate(Pk,n,i)
          n=i

          ! Write to screen
          IF(PRESENT(verbose)) THEN
             IF(verbose) THEN
                WRITE(*,*) 'READ_SIMULATION_POWER_SPECTRUM: Trimmed to Nyquist frequency'
             END IF
          END IF

       END IF
    END IF

    IF(PRESENT(kmin) .AND. PRESENT(kmax)) THEN

       i1=0
       i2=0
       DO i=1,n-1
          IF(k(i)<kmin .AND. k(i+1)>kmin) THEN
             i1=i
          END IF
          IF(k(i)<kmax .AND. k(i+1)>kmax) THEN
             i2=i+1
          END IF
       END DO

       IF(i1==0 .OR. i2==0) THEN
          STOP 'READ_SIMULATION_POWER_SPECTRUM: Error, something went wrong with kmin, kmax'
       END IF

       CALL amputate_general(k,n,m,i1,i2)
       CALL amputate_general(Pk,n,m,i1,i2)
       n=m

    END IF

    ! Write to screen
    IF(PRESENT(verbose)) THEN
       IF(verbose) THEN
          WRITE(*,*) 'READ_SIMULATION_POWER_SPECTRUM: Done'
          WRITE(*,*)
       END IF
    END IF

  END SUBROUTINE read_simulation_power_spectrum

  CHARACTER(len=256) FUNCTION BAHAMAS_power_file_name(model,z,ip)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: model
    REAL, INTENT(IN) :: z
    INTEGER, INTENT(IN) :: ip(2)
    CHARACTER(len=64) :: dir
    CHARACTER(len=32) :: snap, field(2), f1, f2
    LOGICAL :: lexist
    INTEGER :: j

    ! Directory containing everything
    IF(computer==1) dir='/Users/Mead/Physics/data/BAHAMAS/power/M1536'
    IF(computer==2) dir='/home/amead/data/BAHAMAS/power/M1536'

    ! Set the redshift
    IF(z==0.0) THEN
       snap='snap32'
    ELSE IF(z==0.5) THEN
       snap='snap28'
    ELSE IF(z==1.0) THEN
       snap='snap26'
    ELSE IF(z==2.0) THEN
       snap='snap22'
    ELSE
       STOP 'BAHAMAS_POWER_FILE_NAME: Error, redshift specified incorrectly'
    END IF

    ! Set the fields
    DO j=1,2
       IF(ip(j)==field_matter) THEN
          field(j)='all'
       ELSE IF(ip(j)==field_cdm) THEN
          field(j)='dm'
       ELSE IF(ip(j)==field_gas) THEN
          field(j)='gas'
       ELSE IF(ip(j)==field_star) THEN
          field(j)='stars'
       ELSE IF(ip(j)==field_electron_pressure) THEN
          field(j)='epressure'
       ELSE
          WRITE(*,*) 'BAHAMAS_POWER_FILE_NAME: Field number', j
          WRITE(*,*) 'BAHAMAS_POWER_FILE_NAME: Halo type', ip(j)
          STOP 'BAHAMAS_POWER_FILE_NAME: Error, ip specified incorrectly'
       END IF
    END DO

    DO j=1,2

       IF(j==1) THEN
          f1=field(1)
          f2=field(2)
       ELSE IF(j==2) THEN
          f1=field(2)
          f2=field(1)
       ELSE
          STOP 'BAHAMAS_POWER_FILE_NAME: Error, something fucked up'
       END IF

       ! File name
       BAHAMAS_power_file_name=trim(dir)//'/'//trim(model)//'_L400N1024_WMAP9_'//trim(snap)//'_'//trim(f1)//'_'//trim(f2)//'_power.dat'

       ! Check it exists
       INQUIRE(file=BAHAMAS_power_file_name,exist=lexist)

       IF(lexist) THEN
          EXIT
       ELSE IF(j==2) THEN
          WRITE(*,*) 'BAHAMAS_POWER_FILE_NAME: ', trim(BAHAMAS_power_file_name)
          STOP 'BAHAMAS_POWER_FILE_NAME: Error, file does not exist'
       END IF

    END DO

  END FUNCTION BAHAMAS_power_file_name

  SUBROUTINE read_BAHAMAS_power(k,Pk,nk,z,name,field,cosm,kmin,kmax)

    IMPLICIT NONE
    REAL, ALLOCATABLE, INTENT(OUT) :: k(:)
    REAL, ALLOCATABLE, INTENT(OUT) :: Pk(:)
    INTEGER, INTENT(OUT) :: nk
    REAL, INTENT(IN) :: z
    CHARACTER(len=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: field(2)
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL, OPTIONAL, INTENT(IN) :: kmin, kmax
    REAL, ALLOCATABLE :: Pk_DM(:), Pk_HMcode(:)
    CHARACTER(len=256) :: infile, dmonly

    INTEGER, PARAMETER :: field_all_matter(2)=field_matter
    REAL, PARAMETER :: mmin=1e7
    REAL, PARAMETER :: mmax=1e17
    LOGICAL, PARAMETER :: cut_nyquist=.FALSE.
    LOGICAL, PARAMETER :: subtract_shot=.TRUE.
    LOGICAL, PARAMETER :: verbose=.TRUE.

    dmonly=BAHAMAS_power_file_name('DMONLY_2fluid_nu0',z,field_all_matter)
    infile=BAHAMAS_power_file_name(name,z,field)

    CALL read_simulation_power_spectrum(k,Pk_DM,nk,dmonly,kmin,kmax,cut_nyquist,subtract_shot,verbose)
    CALL read_simulation_power_spectrum(k,Pk,   nk,infile,kmin,kmax,cut_nyquist,subtract_shot,verbose)
    Pk=Pk/Pk_DM

    ALLOCATE(Pk_HMcode(nk))
    CALL calculate_HMcode_a(k,scale_factor_z(z),Pk_HMcode,nk,cosm)
    Pk=Pk*Pk_HMcode

  END SUBROUTINE read_BAHAMAS_power

  SUBROUTINE write_fitting_power(base,k,pow_hm,pow_si,ncos,nf,nk,na)

    ! Write fitting data to disk
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: base
    REAL, INTENT(IN) :: k(nk)
    REAL, INTENT(IN) :: pow_hm(ncos,nf,nf,nk,na)
    REAL, INTENT(IN) :: pow_si(ncos,nf,nf,nk,na)
    INTEGER, INTENT(IN) :: ncos
    INTEGER, INTENT(IN) :: nf
    INTEGER, INTENT(IN) :: nk
    INTEGER, INTENT(IN) :: na
    CHARACTER(len=256) :: outfile, outbit
    CHARACTER(len=10) :: uscore, nothing, mid, ext
    INTEGER :: icos, ia, i1, i2, ik

    ! Bits for file name
    uscore='_'
    nothing=''
    mid='_z'
    ext='.dat'

    ! Loop over everything
    DO icos=1,ncos
       DO ia=1,na
          DO i1=1,nf
             DO i2=i1,nf
                outbit=number_file(base,icos,uscore)
                outbit=number_file2(outbit,i1,nothing,i2,mid)
                outfile=number_file(outbit,ia,ext)
                WRITE(*,*) 'WRITE_FITTING_POWER: Outfile: ', trim(outfile)
                OPEN(7,file=outfile)
                DO ik=1,nk
                   WRITE(7,*) k(ik), pow_hm(icos,i1,i2,ik,ia), pow_si(icos,i1,i2,ik,ia)
                END DO
                CLOSE(7)
                WRITE(*,*) 'WRITE_FITTING_POWER: Done'
                WRITE(*,*)
             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE write_fitting_power

  SUBROUTINE finesse_parameters(p,n,name)

    IMPLICIT NONE
    REAL, INTENT(INOUT) :: p(n)
    INTEGER, INTENT(IN) :: n
    CHARACTER(len=*), INTENT(IN) ::  name

    IF(name=='AGN_TUNED_nu0' .OR. name=='AGN_TUNED_nu0_v2' .OR. name=='AGN_TUNED_nu0_v3' .OR. name=='AGN_TUNED_nu0_v4') THEN
       p(param_alpha)=1.54074240    
       p(param_eps)=1.01597583    
       p(param_Gamma)=0.24264216    
       p(param_M0)=log10(6.51173707E+13)
       p(param_Astar)=3.45238000E-02
       p(param_Twhim)=log10(1389362.50)    
       p(param_cstar)=10.4484358    
       p(param_fcold)=6.32560020E-03
       p(param_Mstar)=log10(2.34850982E+12)
       p(param_sstar)=1.25916803    
       p(param_alphap)=-0.513900995    
       p(param_Gammap)=-5.66239981E-03
       p(param_cstarp)=-6.11630008E-02
       p(param_fhot)=1.26100000E-04
       p(param_alphaz)=0.340969592    
       p(param_Gammaz)=0.295596898    
       p(param_M0z)=-8.13719034E-02
       p(param_Astarz)=-0.545276821    
       p(param_Twhimz)=-0.122411400    
       p(param_eta)=-0.266318709
    ELSE IF(name=='AGN_7p6_nu0') THEN
       p(param_alpha)=1.52016437    
       p(param_eps)=1.06684244    
       p(param_Gamma)=0.24494147    
       p(param_M0)=log10(2.84777660E+13)
       p(param_Astar)=3.79242003E-02
       p(param_Twhim)=log10(987513.562)  
       p(param_cstar)=9.78754425    
       p(param_fcold)=1.83899999E-02
       p(param_Mstar)=log10(2.68670049E+12)
       p(param_sstar)=1.13123488    
       p(param_alphap)=-0.531507611    
       p(param_Gammap)=-6.00000005E-03
       p(param_cstarp)=-0.113370098    
       p(param_fhot)=1.08200002E-04
       p(param_alphaz)=0.346108794    
       p(param_Gammaz)=0.231210202    
       p(param_M0z)=-9.32227001E-02
       p(param_Astarz)=-0.536369383    
       p(param_Twhimz)=-0.139042795    
       p(param_eta)=-0.222550094 
    ELSE IF(name=='AGN_8p0_nu0') THEN
       p(param_alpha)=1.45703220    
       p(param_eps)=0.872408926    
       p(param_Gamma)=0.24960959    
       p(param_M0)=log10(1.15050950E+14)
       p(param_Astar)=3.86818014E-02
       p(param_Twhim)=log10(1619561.00) 
       p(param_cstar)=19.4119701    
       p(param_fcold)=1.10999999E-05
       p(param_Mstar)=log10(2.18769510E+12)
       p(param_sstar)=0.803893507    
       p(param_alphap)=-0.528370678    
       p(param_Gammap)=-3.31420009E-03
       p(param_cstarp)=-0.355121315    
       p(param_fhot)=3.90500005E-04
       p(param_alphaz)=0.740169585    
       p(param_Gammaz)=0.354409009    
       p(param_M0z)=-2.40819994E-02
       p(param_Astarz)=-0.425019890    
       p(param_Twhimz)=-8.60318989E-02
       p(param_eta)=-0.243649304   
    ELSE
       STOP 'FINESSE_PARAMETERS: Error, model name not specified correctly'
    END IF

  END SUBROUTINE finesse_parameters

  SUBROUTINE fix_star_parameters(p,n,name)

    IMPLICIT NONE
    REAL, INTENT(INOUT) :: p(n)
    INTEGER, INTENT(IN) :: n
    CHARACTER(len=*), INTENT(IN) ::  name

    IF(name=='AGN_TUNED_nu0' .OR. name=='AGN_TUNED_nu0_v2' .OR. name=='AGN_TUNED_nu0_v3' .OR. name=='AGN_TUNED_nu0_v4') THEN
       p_ori(param_Astar)=0.0486
       p_ori(param_cstar)=6.78
       p_ori(param_mstar)=12.38
       p_ori(param_sstar)=0.711
       p_ori(param_cstarp)=-0.351
       p_ori(param_Astarz)=-0.426
       p_ori(param_eta)=0.000
    ELSE IF(name=='AGN_7p6_nu0') THEN
       p_ori(param_Astar)=0.0466
       p_ori(param_cstar)=7.31
       p_ori(param_mstar)=12.40
       p_ori(param_sstar)=0.844
       p_ori(param_cstarp)=-0.319
       p_ori(param_Astarz)=-0.468
       p_ori(param_eta)=0.000
       p_ori(param_Gammaz)=0.2 ! Otherwise it got stuck
    ELSE IF(name=='AGN_8p0_nu0') THEN
       p_ori(param_Astar)=0.0471
       p_ori(param_cstar)=7.24
       p_ori(param_mstar)=12.32
       p_ori(param_sstar)=0.648
       p_ori(param_cstarp)=-0.396
       p_ori(param_Astarz)=-0.381
       p_ori(param_eta)=0.000
       p_ori(param_Twhimz)=-0.2 ! Otherwise it got stuck
    ELSE
       STOP 'FINESSE_PARAMETERS: Error, model name not specified correctly'
    END IF

  END SUBROUTINE fix_star_parameters

  SUBROUTINE set_parameters(p_nme,p_ori,p_lim,p_lam,p_min,p_max,p_log,n)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(OUT) :: p_nme(n) ! Name
    REAL, INTENT(OUT) :: p_ori(n)    ! Starting value
    REAL, INTENT(OUT) :: p_lim(n)    ! Sensible lower limit
    REAL, INTENT(OUT) :: p_lam(n)    ! Sensible upper limit
    REAL, INTENT(OUT) :: p_min(n)    ! Extreme lower limit
    REAL, INTENT(OUT) :: p_max(n)    ! Extreme upper limit
    LOGICAL, INTENT(OUT) :: p_log(n) ! Explore in log?
    INTEGER , INTENT(IN) :: n        ! Total number of fitting parameters

    p_nme(param_alpha)='alpha'
    p_ori(param_alpha)=1.60!*(1.+z(1))
    p_lim(param_alpha)=1.40
    p_lam(param_alpha)=1.80
    p_min(param_alpha)=1e-2
    p_max(param_alpha)=1e1
    p_log(param_alpha)=.TRUE.

    ! Can not be one because log(1)=0 and this messes things up in setting the ranges
    p_nme(param_eps)='epsilon'
    p_ori(param_eps)=1.1
    p_lim(param_eps)=0.9
    p_lam(param_eps)=1.3
    p_min(param_eps)=1e-2
    p_max(param_eps)=1e2
    p_log(param_eps)=.TRUE.

    p_nme(param_Gamma)='Gamma-1'
    p_ori(param_Gamma)=0.24
    p_lim(param_Gamma)=0.2
    p_lam(param_Gamma)=0.3
    p_min(param_Gamma)=0.01
    p_max(param_Gamma)=2.
    p_log(param_Gamma)=.FALSE.

    ! This must depend on z to ensure parameter has an effect at the higher z when 10^14 haloes are rare
    p_nme(param_M0)='log_M0'
    p_ori(param_M0)=log10(10.**13.77)!/(10.**z(1))) ! Okay with z dependence as long as z(1)=0
    p_lim(param_M0)=log10(10.**13)
    p_lam(param_M0)=log10(10.**14)
    p_min(param_M0)=log10(1e8)
    p_max(param_M0)=log10(1e16)
    p_log(param_M0)=.FALSE.

    p_nme(param_Astar)='A_*'
    !p_ori(param_Astar)=0.042
    p_ori(param_Astar)=0.05
    p_lim(param_Astar)=0.02
    p_lam(param_Astar)=0.06
    p_min(param_Astar)=1e-3
    p_max(param_Astar)=1e-1
    p_log(param_Astar)=.TRUE.

    ! For some reason when 1e6 is set a range is calculated that gives too large a change in figure-of-merit
    p_nme(param_Twhim)='log_T_whim'
    p_ori(param_Twhim)=log10(10.**6.11)
    p_lim(param_Twhim)=log10(10.**6)
    p_lam(param_Twhim)=log10(10.**7)
    p_min(param_Twhim)=log10(1e2)
    p_max(param_Twhim)=log10(1e8)
    p_log(param_Twhim)=.FALSE.

    p_nme(param_cstar)='c_*'
    p_ori(param_cstar)=7.
    p_lim(param_cstar)=5.
    p_lam(param_cstar)=50.
    p_min(param_cstar)=1e0
    p_max(param_cstar)=1e3
    p_log(param_cstar)=.TRUE.

    ! Needed to boost original here so that it has an effect
    p_nme(param_fcold)='f_cold'
    p_ori(param_fcold)=0.00126
    p_lim(param_fcold)=1e-4
    p_lam(param_fcold)=1e-1
    p_min(param_fcold)=1e-5
    p_max(param_fcold)=0.5
    p_log(param_fcold)=.TRUE.

    p_nme(param_Mstar)='log_M_*'
    p_ori(param_Mstar)=log10(10.**12.5)
    p_lim(param_Mstar)=log10(10.**12)
    p_lam(param_Mstar)=log10(10.**13)
    p_min(param_Mstar)=log10(1e8)
    p_max(param_Mstar)=log10(1e16)
    p_log(param_Mstar)=.FALSE.

    p_nme(param_sstar)='sigma_*'
    p_ori(param_sstar)=0.8
    p_lim(param_sstar)=0.5
    p_lam(param_sstar)=2.
    p_min(param_sstar)=0.1
    p_max(param_sstar)=10.
    p_log(param_sstar)=.FALSE.

    p_nme(param_fhot)='f_hot'
    p_ori(param_fhot)=0.01
    p_lim(param_fhot)=1e-3
    p_lam(param_fhot)=1e-1
    p_min(param_fhot)=1e-5
    p_max(param_fhot)=0.5
    p_log(param_fhot)=.TRUE.

    p_nme(param_eta)='eta'
    p_ori(param_eta)=-0.2
    p_lim(param_eta)=-0.6
    p_lam(param_eta)=-0.1
    p_min(param_eta)=-2.0
    p_max(param_eta)=0.0
    p_log(param_eta)=.FALSE.

    p_nme(param_alphap)='alpha_M_pow'
    p_ori(param_alphap)=-0.5
    p_lim(param_alphap)=-0.75
    p_lam(param_alphap)=-0.25
    p_min(param_alphap)=-1.
    p_max(param_alphap)=1.
    p_log(param_alphap)=.FALSE.

    p_nme(param_Gammap)='Gamma_M_pow'
    p_ori(param_Gammap)=-0.02
    p_lim(param_Gammap)=-0.05
    p_lam(param_Gammap)=0.05
    p_min(param_Gammap)=-0.2
    p_max(param_Gammap)=0.2
    p_log(param_Gammap)=.FALSE.

    p_nme(param_cstarp)='c_*_M_pow'
    p_ori(param_cstarp)=-0.2
    p_lim(param_cstarp)=-0.5
    p_lam(param_cstarp)=-0.1
    p_min(param_cstarp)=-1.
    p_max(param_cstarp)=1.
    p_log(param_cstarp)=.FALSE.

    p_nme(param_alphaz)='alpha_z_pow'
    p_ori(param_alphaz)=0.43
    p_lim(param_alphaz)=0.1
    p_lam(param_alphaz)=1.
    p_min(param_alphaz)=-3.
    p_max(param_alphaz)=3.
    p_log(param_alphaz)=.FALSE.

    p_nme(param_Gammaz)='Gamma_z_pow'
    p_ori(param_Gammaz)=0.3
    p_lim(param_Gammaz)=0.
    p_lam(param_Gammaz)=0.6
    p_min(param_Gammaz)=-1.
    p_max(param_Gammaz)=1.
    p_log(param_Gammaz)=.FALSE.

    p_nme(param_M0z)='M0_z_pow'
    p_ori(param_M0z)=-0.08
    p_lim(param_M0z)=-0.2
    p_lam(param_M0z)=-0.01
    p_min(param_M0z)=-1.
    p_max(param_M0z)=1.
    p_log(param_M0z)=.FALSE.

    p_nme(param_Astarz)='Astar_z_pow'
    p_ori(param_Astarz)=-0.45
    p_lim(param_Astarz)=-0.65
    p_lam(param_Astarz)=-0.15
    p_min(param_Astarz)=-1.
    p_max(param_Astarz)=1.
    p_log(param_Astarz)=.FALSE.

    p_nme(param_Twhimz)='Twhim_z_pow'
    p_ori(param_Twhimz)=-0.11
    p_lim(param_Twhimz)=-0.5
    p_lam(param_Twhimz)=-0.01
    p_min(param_Twhimz)=-1.
    p_max(param_Twhimz)=1.
    p_log(param_Twhimz)=.FALSE.

    p_nme(param_HMcode_Dv0)='Dv0'
    p_ori(param_HMcode_Dv0)=418.
    p_lim(param_HMcode_Dv0)=100.
    p_lam(param_HMcode_Dv0)=600.
    p_min(param_HMcode_Dv0)=50.
    p_max(param_HMcode_Dv0)=1000.
    p_log(param_HMcode_Dv0)=.FALSE.

    p_nme(param_HMcode_Dvp)='Dvp'
    p_ori(param_HMcode_Dvp)=-0.352    
    p_lim(param_HMcode_Dvp)=-0.5
    p_lam(param_HMcode_Dvp)=-0.2
    p_min(param_HMcode_Dvp)=-5.
    p_max(param_HMcode_Dvp)=5.
    p_log(param_HMcode_Dvp)=.FALSE.

    p_nme(param_HMcode_dc0)='dc0'
    p_ori(param_HMcode_dc0)=1.59
    p_lim(param_HMcode_dc0)=1.55
    p_lam(param_HMcode_dc0)=1.65
    p_min(param_HMcode_dc0)=1.
    p_max(param_HMcode_dc0)=2.
    p_log(param_HMcode_dc0)=.FALSE.

    p_nme(param_HMcode_dcp)='dcp'
    p_ori(param_HMcode_dcp)=0.0314
    p_lim(param_HMcode_dcp)=0.01
    p_lam(param_HMcode_dcp)=0.1
    p_min(param_HMcode_dcp)=-5.
    p_max(param_HMcode_dcp)=5.
    p_log(param_HMcode_dcp)=.FALSE.

    p_nme(param_HMcode_eta0)='eta0'
    p_ori(param_HMcode_eta0)=0.603
    p_lim(param_HMcode_eta0)=-1.
    p_lam(param_HMcode_eta0)=1.
    p_min(param_HMcode_eta0)=-5.
    p_max(param_HMcode_eta0)=5.
    p_log(param_HMcode_eta0)=.FALSE.

    p_nme(param_HMcode_eta1)='eta1'
    p_ori(param_HMcode_eta1)=0.300
    p_lim(param_HMcode_eta1)=0.1
    p_lam(param_HMcode_eta1)=0.7
    p_min(param_HMcode_eta1)=-5.
    p_max(param_HMcode_eta1)=5.
    p_log(param_HMcode_eta1)=.FALSE.

    p_nme(param_HMcode_f0)='f0'
    !p_ori(param_HMcode_f0)=0.0095 ! Mead (2016) 
    p_ori(param_HMcode_f0)=0.188 ! Mead (2015) damping
    p_lim(param_HMcode_f0)=-1.
    p_lam(param_HMcode_f0)=1.
    p_min(param_HMcode_f0)=-5.
    p_max(param_HMcode_f0)=5.
    p_log(param_HMcode_f0)=.FALSE.

    p_nme(param_HMcode_fp)='fp'
    !p_ori(param_HMcode_fp)=1.37 ! Mead (2016)
    p_ori(param_HMcode_fp)=4.29 ! Mead (2015)
    p_lim(param_HMcode_fp)=1.
    p_lam(param_HMcode_fp)=6.
    p_min(param_HMcode_fp)=-10.
    p_max(param_HMcode_fp)=10.
    p_log(param_HMcode_fp)=.FALSE.

    p_nme(param_HMcode_kstar)='kstar'
    p_ori(param_HMcode_kstar)=0.584
    p_lim(param_HMcode_kstar)=-1.
    p_lam(param_HMcode_kstar)=1.
    p_min(param_HMcode_kstar)=-5.
    p_max(param_HMcode_kstar)=5.
    p_log(param_HMcode_kstar)=.FALSE.

    p_nme(param_HMcode_As)='As'
    p_ori(param_HMcode_As)=3.13
    p_lim(param_HMcode_As)=2.
    p_lam(param_HMcode_As)=6.
    p_min(param_HMcode_As)=1.
    p_max(param_HMcode_As)=10.
    p_log(param_HMcode_As)=.FALSE.

    p_nme(param_HMcode_alpha0)='alpha0'
    p_ori(param_HMcode_alpha0)=3.24
    p_lim(param_HMcode_alpha0)=1.
    p_lam(param_HMcode_alpha0)=4.
    p_min(param_HMcode_alpha0)=-5.
    p_max(param_HMcode_alpha0)=5.
    p_log(param_HMcode_alpha0)=.FALSE.

    p_nme(param_HMcode_alpha1)='alpha1'
    p_ori(param_HMcode_alpha1)=1.85
    p_lim(param_HMcode_alpha1)=1.
    p_lam(param_HMcode_alpha1)=3.
    p_min(param_HMcode_alpha1)=-5.
    p_max(param_HMcode_alpha1)=5.
    p_log(param_HMcode_alpha1)=.FALSE.

  END SUBROUTINE set_parameters

END PROGRAM HMx_fitting
