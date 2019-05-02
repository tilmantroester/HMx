PROGRAM HMx_fitting

  USE HMx
  USE cosmology_functions
  USE cosmic_emu_stuff
  USE string_operations
  USE random_numbers
  USE special_functions
  USE owls
  USE owls_extras

  IMPLICIT NONE

  TYPE fitting
     REAL, ALLOCATABLE :: minimum(:)            ! Minimum value this parameter is allowed to take
     REAL, ALLOCATABLE :: maximum(:)            ! Maximum value this parameter is allowed to take
     REAL, ALLOCATABLE :: original(:)           ! Starting value for this parameter
     REAL, ALLOCATABLE :: sigma(:)              ! Sigma of jump size for this parameter
     REAL, ALLOCATABLE :: best(:)               ! Best-fitting value for this parameter
     LOGICAL, ALLOCATABLE :: set(:)             ! Should this parameter be fitted?
     LOGICAL, ALLOCATABLE :: log(:)             ! Should this parameter be fitted in log?     
     LOGICAL, ALLOCATABLE :: cov(:)             ! Should we determine the jump size for this parameter?     
     CHARACTER(len=256), ALLOCATABLE :: name(:) ! Name of this parameter
     INTEGER :: n                               ! Number of parameters
  END TYPE fitting

  INTEGER :: imode
  INTEGER :: ncos, nfields, nz, nk
  REAL, ALLOCATABLE :: k(:), z(:), pow(:,:,:,:,:), weight(:,:,:,:,:)
  INTEGER, ALLOCATABLE :: fields(:)
  TYPE(fitting) :: fit
  TYPE(halomod), ALLOCATABLE :: hmod(:)
  TYPE(cosmology), ALLOCATABLE :: cosm(:)
  CHARACTER(len=256) :: name, base, mode, zin, outbase, numchain, maxtime, accuracy, response
  REAL :: delta
  REAL :: tmax
  INTEGER :: nchain
  LOGICAL :: resp_BAHAMAS

  ! Halo model calculation parameters
  REAL, PARAMETER :: mmin=mmin_HMx ! Minimum halo mass for the calculation
  REAL, PARAMETER :: mmax=mmax_HMx ! Maximum halo mass for the calculation

  INTEGER, PARAMETER :: m=huge(m)            ! Re-evaluate range every 'm' points
  !INTEGER, PARAMETER :: m=100                ! Re-evaluate range every 'm' points
  INTEGER, PARAMETER :: seed=0               ! Random-number seed
  LOGICAL, PARAMETER :: random_start=.FALSE. ! Start from a random point within the prior range
  LOGICAL, PARAMETER :: mcmc=.TRUE.          ! Accept worse figure of merit with some probability
  INTEGER, PARAMETER :: computer=1           ! Which computer are you on?
  LOGICAL, PARAMETER :: set_ranges=.TRUE.    ! Set the parameter ranges or not

  ! k cuts for the BAHAMAS power spectra
  REAL, PARAMETER :: kmin_default=1e-4       ! Minimum wavenumber to read in [h/Mpc]
  REAL, PARAMETER :: kmax_default=1e2        ! Minimum wavenumber to read in [h/Mpc]
  REAL, PARAMETER :: kmin_BAHAMAS=0.03       ! Minimum wavenumber to use [h/Mpc]
  !REAL, PARAMETER :: kmin_BAHAMAS=0.15       ! Minimum wavenumber to use [h/Mpc]
  !REAL, PARAMETER :: kmax_BAHAMAS_z0p0=10.   ! Maximum wavenumber to use [h/Mpc] at z=0.0
  !REAL, PARAMETER :: kmax_BAHAMAS_z0p5=4.    ! Maximum wavenumber to use [h/Mpc] at z=0.5
  !REAL, PARAMETER :: kmax_BAHAMAS_z1p0=2.    ! Maximum wavenumber to use [h/Mpc] at z=1.0
  !REAL, PARAMETER :: kmax_BAHAMAS_z2p0=1.    ! Maximum wavenumber to use [h/Mpc] at z=2.0
  REAL, PARAMETER :: kmax_BAHAMAS=10.        ! Maximum wavenumber to us [h/Mpc]
  REAL, PARAMETER :: z_default=0.0           ! Default z to fit if not specified
  INTEGER, PARAMETER :: mesh=1024            ! BAHAMAS mesh size power to use
  LOGICAL, PARAMETER :: cut_nyquist=.TRUE.   ! Should the BAHAMAS measured P(k) be cut above the Nyquist frequency?
  LOGICAL, PARAMETER :: subtract_shot=.TRUE. ! Should the BAHAMAS measured P(k) have shot-noise subtracted?
  LOGICAL, PARAMETER :: response_default=.TRUE.                   ! Should I treat the BAHAMAS P(k) as HMcode response?
  CHARACTER(len=256), PARAMETER :: outbase_default='fitting/test' ! Default output file

  ! Read in starting option for which parameter fitting to run
  CALL get_command_argument(1,mode)
  IF(mode=='') THEN
     imode=-1
  ELSE
     READ(mode,*) imode
  END IF

  ! Decide what to do
  IF(imode==-1) THEN
     WRITE(*,*)
     WRITE(*,*) 'HMx_FITTING: Choose what to do'
     WRITE(*,*) '=============================='
     WRITE(*,*) ' 1 - Mira Titan nodes'
     WRITE(*,*) ' 2 - FrankenEmu nodes'
     WRITE(*,*) ' 3 - Random Mira Titan'
     WRITE(*,*) ' 4 - Random FrankenEmu'
     WRITE(*,*) '11 - Hydro: fixed z; final parameters; gas (isothermal)'
     WRITE(*,*) '12 - Hydro: fixed z; final parameters; pressure (isothermal)'
     WRITE(*,*) '13 - Hydro: fixed z; final parameters; matter and pressure (no pressure-pressure; isothermal)'
     WRITE(*,*) '14 - Hydro: fixed z; final parameters; matter (isothermal)'
     WRITE(*,*) '15 - Hydro: fixed z; final parameters; CDM, gas, stars'
     WRITE(*,*) '16 - Hydro: fixed z; final parameters; matter, CDM, gas, stars'
     WRITE(*,*) '17 - N/A'
     WRITE(*,*) '18 - N/A'
     WRITE(*,*) '19 - N/A'
     WRITE(*,*) '20 - Hydro: fixed z; final parameters; matter and pressure (no pressure-pressure)'
     WRITE(*,*) '21 - Hydro: extended z; final parameters; matter and pressure (no pressure-pressure)'
     WRITE(*,*) '22 - N/A'
     WRITE(*,*) '23 - N/A'
     WRITE(*,*) '24 - N/A'
     WRITE(*,*) '25 - N/A'
     WRITE(*,*) '26 - Hydro: fixed z; basic parameters; CDM'
     WRITE(*,*) '27 - Hydro: fixed z; basic parameters; matter, CDM, gas, stars'
     WRITE(*,*) '28 - Hydro: fixed z; basic parameters; matter, CDM, gas, stars, pressure (no pressure-pressure)'
     WRITE(*,*) '29 - Hydro: fixed z; basic parameters; matter'
     WRITE(*,*) '30 - Hydro: fixed z; basic parameters; pressure'
     WRITE(*,*) '31 - Hydro: fixed z; basic parameters; stars'
     WRITE(*,*) '32 - Hydro: fixed z; basic parameters; gas'
     WRITE(*,*) '33 - Hydro: fixed z; final parameters; CDM'
     WRITE(*,*) '34 - Hydro: fixed z; final parameters; gas'
     WRITE(*,*) '35 - Hydro: fixed z; final parameters; stars'
     WRITE(*,*) '36 - Hydro: fixed z; final parameters; pressure'
     WRITE(*,*) '37 - Hydro: fixed z; final parameters; CDM, gas, stars'
     WRITE(*,*) '38 - Hydro: fixed z; final parameters; matter, CDM, gas, stars'
     WRITE(*,*) '39 - Hydro: fixed z; final parameters; matter'
     WRITE(*,*) '40 - Hydro: fixed z; final parameters; matter, CDM, gas, stars, pressure (no pressure-pressure)'
     WRITE(*,*) '41 - Hydro: extended z; final parameters; CDM'
     WRITE(*,*) '42 - Hydro: extended z; final parameters; gas'
     WRITE(*,*) '43 - Hydro: extended z; final parameters; stars'
     WRITE(*,*) '44 - Hydro: extended z; final parameters; matter'
     WRITE(*,*) '45 - Hydro: extended z; final parameters; CDM, gas, stars'
     WRITE(*,*) '46 - Hydro: extended z; final parameters; matter, CDM, gas, stars'
     READ(*,*) imode 
     WRITE(*,*)
  END IF

  ! Read in chain length
  CALL get_command_argument(2,numchain)
  IF(numchain=='') THEN
     WRITE(*,*) 'HMx_FITTING: Specify points in fitting chain'
     READ(*,*) nchain
     WRITE(*,*)
  ELSE
     READ(numchain,*) nchain
  END IF

  ! Read in maximum time [mins]
  CALL get_command_argument(3,maxtime)
  IF(maxtime=='') THEN
     !tmax=tmax_default
     WRITE(*,*) 'HMx_FITTING: Specify maximum time [min]:'
     READ(*,*) tmax
     WRITE(*,*)
  ELSE
     READ(maxtime,*) tmax
  END IF

  ! Accuracy
  CALL get_command_argument(4,accuracy)
  IF(accuracy=='') THEN
     !delta=delta_default
     WRITE(*,*) 'HMx_FITTING: Specify fitting accuracy:'
     READ(*,*) delta
     WRITE(*,*)
  ELSE
     READ(accuracy,*) delta
  END IF

  ! Read in outfile
  CALL get_command_argument(5,outbase)
  IF(outbase=='') outbase=outbase_default

  ! Read in BAHAMAS simulation name
  CALL get_command_argument(6,name)

  ! Read in response logical that determines if we fit the BAHAMAS response or not
  CALL get_command_argument(7,response)
  IF(response=='') THEN
     resp_BAHAMAS=response_default
  ELSE IF(response=='TRUE') THEN
     resp_BAHAMAS=.TRUE.
  ELSE IF(response=='FALSE') THEN
     resp_BAHAMAS=.FALSE.
  ELSE
     STOP 'HMx_FITTING: Error, response logical specified incorrectly (TRUE or FALSE)'
  END IF

  ! Read in BAHAMAS simulation redshift if doing fixed z
  CALL get_command_argument(8,zin)

  ! Set the random-number generator
  CALL RNG_set(seed)

  ! Initial white space
  WRITE(*,*)

  ! Write useful info to screen
  WRITE(*,*) 'HMx_FITTING: Fitting routine'
  WRITE(*,*) 'HMx_FITTING: Mode:', imode
  WRITE(*,*) 'HMx_FITTING: Number of points in chain:', nchain
  WRITE(*,*) 'HMx_FITTING: Maximum run time [mins]:', tmax
  WRITE(*,*) 'HMx_FITTING: Maximum run time [hours]:', tmax/60.
  WRITE(*,*) 'HMx_FITTING: Accuracy:', delta
  WRITE(*,*) 'HMx_FITTING: Output file base: ', TRIM(outbase)
  WRITE(*,*) 'HMx_FITTING: Random number seed:', seed
  WRITE(*,*) 'HMx_FITTING: Random start:', random_start
  WRITE(*,*) 'HMx_FITTING: MCMC mode:', mcmc
  !WRITE(*,*) 'HMx_FITTING: Re-evaluate covariance:', m
  WRITE(*,*)

  ! Set the cosmological models
  CALL init_cosmologies(imode,cosm,ncos)

  ! Set the fields
  CALL init_fields(imode,fields,nfields)

  ! Set the redshifts
  CALL init_redshifts(imode,z,nz)

  ! Set the halo models
  CALL init_halomods(imode,hmod,ncos)

  ! Print one halo model out to check it looks sensible
  CALL init_halomod(mmin,mmax,scale_factor_z(z(1)),hmod(1),cosm(1),verbose=.TRUE.)
  CALL print_halomod(hmod(1),cosm(1),verbose=.TRUE.)

  ! Read in the simulation power spectra
  CALL read_simulation_power_spectra(imode,k,nk,pow,cosm,ncos,z,nz,fields,nfields)

  ! Set the weights
  CALL init_weights(imode,weight,ncos,nfields,nk,nz)

  ! Set the initial parameter values, ranges, names, etc.
  CALL init_parameters(fit)

  ! Set the parameters that will be varied
  CALL init_mode(imode,fit)

  ! Do the actual fitting algorithm
  CALL actual_fitting(fit%original,fit%n,delta,tmax,nchain,k,nk,z,nz,fields,nfields,weight,pow,fit,hmod,cosm,ncos)

CONTAINS

  SUBROUTINE init_parameters(fit)

    ! Assigns values to parameter name, original value, sigma, minimum and maximuim values and logexp
    IMPLICIT NONE
    TYPE(fitting), INTENT(INOUT) :: fit
    INTEGER, PARAMETER ::  n=param_n
    INTEGER :: i

    CALL allocate_arrays(fit,n)

    !! Hydro !!

    fit%name(param_alpha)='alpha'
    fit%original(param_alpha)=1. ! Standard value from ihm=3
    !fit%original(param_alpha)=0.90!*(1.+z(1))
    fit%sigma(param_alpha)=0.1
    fit%minimum(param_alpha)=1e-2
    fit%maximum(param_alpha)=1e1
    fit%log(param_alpha)=.TRUE.

    fit%name(param_beta)='beta'
    fit%original(param_beta)=1. ! Standard value from ihm=3
    !fit%original(param_beta)=0.90!*(1.+z(1))
    fit%sigma(param_beta)=0.1
    fit%minimum(param_beta)=1e-2
    fit%maximum(param_beta)=1e1
    fit%log(param_beta)=.TRUE.

    fit%name(param_eps)='epsilon'
    !fit%original(param_eps)=1. ! Standard value from ihm=3
    fit%original(param_eps)=1.1
    fit%sigma(param_eps)=0.05
    fit%minimum(param_eps)=1e-2
    fit%maximum(param_eps)=1e2
    !fit%log(param_eps)=.FALSE. ! Changed on 25/04/2019; may cause problems?
    fit%log(param_eps)=.TRUE.

    ! NOTE: This parameter is Gamma-1, not regular Gamma
    fit%name(param_Gamma)='Gamma-1'
    fit%original(param_Gamma)=0.17 ! Standard value from ihm=3
    !fit%original(param_Gamma)=0.24
    fit%sigma(param_Gamma)=0.005
    fit%minimum(param_Gamma)=0.01
    fit%maximum(param_Gamma)=2.
    fit%log(param_Gamma)=.FALSE.

    ! This should depend on z to ensure parameter has an effect at the higher z when 10^14 haloes are rare
    ! NOTE: This parameter is already log10
    fit%name(param_M0)='log_M0'
    !fit%original(param_M0)=log10(10.**13.5) ! Standard value from ihm=3
    !fit%original(param_M0)=log10((10.**13.5)/(10.**(z(1)/2.))) ! Okay with z dependence as long as z(1)=0
    fit%original(param_M0)=13.5-z(1)/2.
    fit%sigma(param_M0)=0.1
    fit%minimum(param_M0)=log10(1e8)
    fit%maximum(param_M0)=log10(1e16)
    fit%log(param_M0)=.FALSE.

    fit%name(param_Astar)='A_*'
    fit%original(param_Astar)=0.03 ! Standard value from ihm=3
    !fit%original(param_Astar)=0.042
     fit%sigma(param_Astar)=0.002
    fit%minimum(param_Astar)=1e-4
    fit%maximum(param_Astar)=0.2
    fit%log(param_Astar)=.TRUE.

    ! NOTE: This parameter is already log10
    fit%name(param_Twhim)='log_T_whim'
    !fit%original(param_Twhim)=log10(10.**6.47712) ! Standard value from ihm=3
    fit%original(param_Twhim)=log10(10.**6.5)
    fit%sigma(param_Twhim)=log10(10.**0.1)
    fit%minimum(param_Twhim)=log10(1e2)
    fit%maximum(param_Twhim)=log10(1e8)
    fit%log(param_Twhim)=.FALSE.

    fit%name(param_cstar)='c_*'
    fit%original(param_cstar)=10. ! Standard value from ihm=3
    !fit%original(param_cstar)=7.
    fit%sigma(param_cstar)=0.1
    fit%minimum(param_cstar)=1e-3
    fit%maximum(param_cstar)=1e3
    fit%log(param_cstar)=.TRUE.

    fit%name(param_eta)='eta'
    !fit%original(param_eta)=0.0 ! Standard value from ihm=3
    fit%original(param_eta)=-0.3 
    fit%sigma(param_eta)=0.05
    fit%minimum(param_eta)=-2.0
    fit%maximum(param_eta)=0.0
    fit%log(param_eta)=.FALSE.

    fit%name(param_sstar)='sigma_*'
    fit%original(param_sstar)=1.2 ! Standard value from ihm=3
    !fit%original(param_sstar)=0.8     
    fit%sigma(param_sstar)=0.02
    fit%minimum(param_sstar)=0.1
    fit%maximum(param_sstar)=10.
    fit%log(param_sstar)=.FALSE.

    ! NOTE: This parameter is already log10
    fit%name(param_Mstar)='log_M_*'
    !fit%original(param_Mstar)=log10(10.**12.6987) ! Standard value from ihm=3
    fit%original(param_Mstar)=log10(10.**12.5)    
    fit%sigma(param_Mstar)=0.1
    fit%minimum(param_Mstar)=log10(1e8)
    fit%maximum(param_Mstar)=log10(1e16)
    fit%log(param_Mstar)=.FALSE.

    fit%name(param_fcold)='f_cold'
    !fit%original(param_fcold)=0. ! Standard value from ihm=3
    fit%original(param_fcold)=1e-2   
    fit%sigma(param_fcold)=0.1
    fit%minimum(param_fcold)=1e-5
    fit%maximum(param_fcold)=0.5
    fit%log(param_fcold)=.TRUE.

    fit%name(param_fhot)='f_hot'
    !fit%original(param_fhot)=0. ! Standard value from ihm=3
    fit%original(param_fhot)=1e-2    
    fit%sigma(param_fhot)=0.1
    fit%minimum(param_fhot)=1e-5
    fit%maximum(param_fhot)=0.5
    fit%log(param_fhot)=.TRUE.

    fit%name(param_ibeta)='iso_beta'
    fit%original(param_ibeta)=2./3. ! Standard value from ihm=3
    fit%sigma(param_ibeta)=0.01
    fit%minimum(param_ibeta)=0.1
    fit%maximum(param_ibeta)=2.
    fit%log(param_ibeta)=.FALSE.

    !! !!

    !! Hydro - M indices !!

    fit%name(param_alphap)='alpha_M_pow'
    fit%original(param_alphap)=0.0 ! Standard value from ihm=3
    !fit%original(param_alphap)=-0.5
    fit%sigma(param_alphap)=0.01
    fit%minimum(param_alphap)=-1.
    fit%maximum(param_alphap)=1.
    fit%log(param_alphap)=.FALSE.

    fit%name(param_betap)='beta_M_pow'
    fit%original(param_betap)=0.0 ! Standard value from ihm=3
    !fit%original(param_betap)=-0.5
    fit%sigma(param_betap)=0.01
    fit%minimum(param_betap)=-1.
    fit%maximum(param_betap)=1.
    fit%log(param_betap)=.FALSE.

    fit%name(param_Gammap)='Gamma_M_pow'
    fit%original(param_Gammap)=0.0 ! Standard value from ihm=3
    !fit%original(param_Gammap)=-0.02
    fit%sigma(param_Gammap)=0.001
    fit%minimum(param_Gammap)=-0.2
    fit%maximum(param_Gammap)=0.2
    fit%log(param_Gammap)=.FALSE.

    fit%name(param_Astarp)='Astar_M_pow'
    fit%original(param_Astarp)=0.0 ! Standard value from ihm=3
    fit%sigma(param_Astarp)=0.01
    fit%minimum(param_Astarp)=-1.
    fit%maximum(param_Astarp)=1.
    fit%log(param_Astarp)=.FALSE.

    fit%name(param_cstarp)='c_*_M_pow'
    fit%original(param_cstarp)=0.0 ! Standard value from ihm=3
    !fit%original(param_cstarp)=-0.2
    fit%sigma(param_cstarp)=0.01
    fit%minimum(param_cstarp)=-1.
    fit%maximum(param_cstarp)=1.
    fit%log(param_cstarp)=.FALSE.

    fit%name(param_ibetap)='iso_beta_M_pow'
    fit%original(param_ibetap)=0.0 ! Standard value from ihm=3
    fit%sigma(param_ibetap)=0.01
    fit%minimum(param_ibetap)=-0.3 ! Problem if max set to higher than 0.5
    fit%maximum(param_ibetap)=0.3 ! Problem if set to higher than 0.5
    fit%log(param_ibetap)=.FALSE.

    !! !!

    !! Hydro - z indices !!

    fit%name(param_alphaz)='alpha_z_pow'
    fit%original(param_alphaz)=0.0 ! Standard value from ihm=3
    !fit%original(param_alphaz)=0.43
    fit%sigma(param_alphaz)=0.01
    fit%minimum(param_alphaz)=-3.
    fit%maximum(param_alphaz)=3.
    fit%log(param_alphaz)=.FALSE.

    fit%name(param_betaz)='beta_z_pow'
    fit%original(param_betaz)=0.0 ! Standard value from ihm=3
    !fit%original(param_betaz)=0.43
    fit%sigma(param_betaz)=0.01
    fit%minimum(param_betaz)=-3.
    fit%maximum(param_betaz)=3.
    fit%log(param_betaz)=.FALSE.

    fit%name(param_epsz)='eps_z_pow'
    fit%original(param_epsz)=0.0 ! Standard value from ihm=3
    !fit%original(param_epsz)=0.5
    fit%sigma(param_epsz)=0.01
    fit%minimum(param_epsz)=-3.
    fit%maximum(param_epsz)=3.
    fit%log(param_epsz)=.FALSE.

    fit%name(param_Gammaz)='Gamma_z_pow'
    fit%original(param_Gammaz)=0.0 ! Standard value from ihm=3
    !fit%original(param_Gammaz)=0.3
    fit%sigma(param_Gammaz)=0.01
    fit%minimum(param_Gammaz)=-1.
    fit%maximum(param_Gammaz)=1.
    fit%log(param_Gammaz)=.FALSE.

    fit%name(param_M0z)='M0_z_pow'
    fit%original(param_M0z)=0.0 ! Standard value from ihm=3
    !fit%original(param_M0z)=-0.08
    fit%sigma(param_M0z)=0.01
    fit%minimum(param_M0z)=-1.
    fit%maximum(param_M0z)=1.
    fit%log(param_M0z)=.FALSE.

    fit%name(param_Astarz)='Astar_z_pow'
    fit%original(param_Astarz)=0.0 ! Standard value from ihm=3
    !fit%original(param_Astarz)=-0.45
    fit%sigma(param_Astarz)=0.01
    fit%minimum(param_Astarz)=-1.
    fit%maximum(param_Astarz)=1.
    fit%log(param_Astarz)=.FALSE.

    fit%name(param_Twhimz)='Twhim_z_pow'
    fit%original(param_Twhimz)=0.0 ! Standard value from ihm=3
    !fit%original(param_Twhimz)=-0.11
    fit%sigma(param_Twhimz)=0.01
    fit%minimum(param_Twhimz)=-1.
    fit%maximum(param_Twhimz)=1.
    fit%log(param_Twhimz)=.FALSE.

    fit%name(param_cstarz)='c_*_z_pow'
    fit%original(param_cstarz)=0.0 ! Standard value from ihm=3
    fit%sigma(param_cstarz)=0.01
    fit%minimum(param_cstarz)=-3.
    fit%maximum(param_cstarz)=3.
    fit%log(param_cstarz)=.FALSE.

    fit%name(param_Mstarz)='M_*_z_pow'
    fit%original(param_Mstarz)=0.0 ! Standard value from ihm=3
    fit%sigma(param_Mstarz)=0.01
    fit%minimum(param_Mstarz)=-1.
    fit%maximum(param_Mstarz)=1.
    fit%log(param_Mstarz)=.FALSE.

    fit%name(param_ibetaz)='iso_beta_z_pow'
    fit%original(param_ibetaz)=0.0 ! Standard value from ihm=3
    fit%sigma(param_ibetaz)=0.01
    fit%minimum(param_ibetaz)=-1.
    fit%maximum(param_ibetaz)=1.
    fit%log(param_ibetaz)=.FALSE.

    !!

    !! HMcode !!

    fit%name(param_HMcode_Dv0)='Dv0'
    fit%original(param_HMcode_Dv0)=418.
    fit%sigma(param_HMcode_Dv0)=1.
    fit%minimum(param_HMcode_Dv0)=50.
    fit%maximum(param_HMcode_Dv0)=1000.
    fit%log(param_HMcode_Dv0)=.FALSE.

    fit%name(param_HMcode_Dvp)='Dvp'
    fit%original(param_HMcode_Dvp)=-0.352
    fit%sigma(param_HMcode_Dvp)=0.02
    fit%minimum(param_HMcode_Dvp)=-5.
    fit%maximum(param_HMcode_Dvp)=5.
    fit%log(param_HMcode_Dvp)=.FALSE.

    fit%name(param_HMcode_dc0)='dc0'
    fit%original(param_HMcode_dc0)=1.590
    fit%sigma(param_HMcode_dc0)=0.005
    fit%minimum(param_HMcode_dc0)=1.
    fit%maximum(param_HMcode_dc0)=2.
    fit%log(param_HMcode_dc0)=.FALSE.

    fit%name(param_HMcode_dcp)='dcp'
    fit%original(param_HMcode_dcp)=0.0314
    fit%sigma(param_HMcode_dcp)=0.002
    fit%minimum(param_HMcode_dcp)=-5.
    fit%maximum(param_HMcode_dcp)=5.
    fit%log(param_HMcode_dcp)=.FALSE.

    fit%name(param_HMcode_eta0)='eta0'
    fit%original(param_HMcode_eta0)=0.603
    fit%sigma(param_HMcode_eta0)=0.02
    fit%minimum(param_HMcode_eta0)=-5.
    fit%maximum(param_HMcode_eta0)=5.
    fit%log(param_HMcode_eta0)=.FALSE.

    fit%name(param_HMcode_eta1)='eta1'
    fit%original(param_HMcode_eta1)=0.300
    fit%sigma(param_HMcode_eta1)=0.02
    fit%minimum(param_HMcode_eta1)=-5.
    fit%maximum(param_HMcode_eta1)=5.
    fit%log(param_HMcode_eta1)=.FALSE.

    fit%name(param_HMcode_f0)='f0'
    !fit%original(param_HMcode_f0)=0.0095 ! Mead (2016) 
    fit%original(param_HMcode_f0)=0.188 ! Mead (2015) damping
    fit%sigma(param_HMcode_f0)=0.01
    fit%minimum(param_HMcode_f0)=-5.
    fit%maximum(param_HMcode_f0)=5.
    fit%log(param_HMcode_f0)=.FALSE.

    fit%name(param_HMcode_fp)='fp'
    !fit%original(param_HMcode_fp)=1.37 ! Mead (2016)
    fit%original(param_HMcode_fp)=4.29 ! Mead (2015)
    fit%sigma(param_HMcode_fp)=0.1
    fit%minimum(param_HMcode_fp)=-10.
    fit%maximum(param_HMcode_fp)=10.
    fit%log(param_HMcode_fp)=.FALSE.

    fit%name(param_HMcode_kstar)='kstar'
    fit%original(param_HMcode_kstar)=0.584
    fit%sigma(param_HMcode_kstar)=0.05
    fit%minimum(param_HMcode_kstar)=-5.
    fit%maximum(param_HMcode_kstar)=5.
    fit%log(param_HMcode_kstar)=.FALSE.

    fit%name(param_HMcode_As)='As'
    fit%original(param_HMcode_As)=3.13
    fit%sigma(param_HMcode_As)=0.1
    fit%minimum(param_HMcode_As)=1.
    fit%maximum(param_HMcode_As)=10.
    fit%log(param_HMcode_As)=.FALSE.

    fit%name(param_HMcode_alpha0)='alpha0'
    fit%original(param_HMcode_alpha0)=3.24
    fit%sigma(param_HMcode_alpha0)=0.1
    fit%minimum(param_HMcode_alpha0)=-5.
    fit%maximum(param_HMcode_alpha0)=5.
    fit%log(param_HMcode_alpha0)=.FALSE.

    fit%name(param_HMcode_alpha1)='alpha1'
    fit%original(param_HMcode_alpha1)=1.85
    fit%sigma(param_HMcode_alpha1)=0.1
    fit%minimum(param_HMcode_alpha1)=-5.
    fit%maximum(param_HMcode_alpha1)=5.
    fit%log(param_HMcode_alpha1)=.FALSE.

    !! !!

    ! Initially assume all parameters are not being varied
    fit%set=.FALSE.

    ! Initially assume we calculate the step size for each parameter that is being varied
    fit%cov=.TRUE.

    ! Check that starting value is not outside minimum and maximum values
    DO i=1,fit%n
       IF(fit%original(i)<fit%minimum(i) .OR. fit%original(i)>fit%maximum(i)) THEN
          WRITE(*,*) 'INIT_PARAMETERS: Parameter number:', i
          WRITE(*,*) 'INIT_PARAMETERS: Parameter name:', trim(fit%name(i))
          WRITE(*,*) 'INIT_PARAMETERS: Parameter original:', fit%original(i)
          WRITE(*,*) 'INIT_PARAMETERS: Parameter minimum:', fit%minimum(i)
          WRITE(*,*) 'INIT_PARAMETERS: Parameter maximum:', fit%maximum(i)
          STOP 'INIT_PARAMETERS: Error, parameter original value is below minimum or above maximum value'
       END IF
    END DO

    ! Loop over parameters and take logarithms of those that need it, the internal parameters are then set
    DO i=1,fit%n
       IF(fit%log(i)) THEN
          IF(fit%original(i)<=0. .OR. fit%minimum(i)<=0. .OR. fit%maximum(i)<=0.) THEN
             WRITE(*,*) 'INIT_PARAMETERS: Parameter number:', i
             WRITE(*,*) 'INIT_PARAMETERS: Parameter name:', trim(fit%name(i))
             WRITE(*,*) 'INIT_PARAMETERS: Parameter original:', fit%original(i)
             WRITE(*,*) 'INIT_PARAMETERS: Parameter minimum:', fit%minimum(i)
             WRITE(*,*) 'INIT_PARAMETERS: Parameter maximum:', fit%maximum(i)
             STOP 'INIT_PARAMETERS: Error, you are about to take the logarithm of zero or a negative number'
          END IF
          fit%original(i)=log10(fit%original(i))
          fit%minimum(i)=log10(fit%minimum(i))
          fit%maximum(i)=log10(fit%maximum(i))
       END IF
    END DO

    ! If we are doing a random start then pick parameter uniform-randomly between minimum and maximum allowed value
    IF(random_start) THEN
       DO i=1,fit%n
          fit%original(i)=random_uniform(fit%minimum(i),fit%maximum(i))
       END DO
    END IF

    ! Set the 'best' parameters to be equal to the 'original' parameters
    fit%best=fit%original
    
  END SUBROUTINE init_parameters

  SUBROUTINE allocate_arrays(fit,n)

    ! Allocate all arrays associated with the fitting
    IMPLICIT NONE
    TYPE(fitting), INTENT(INOUT) :: fit
    INTEGER, INTENT(IN) :: n

    fit%n=n

    ALLOCATE(fit%minimum(n))
    ALLOCATE(fit%maximum(n))
    ALLOCATE(fit%original(n))
    ALLOCATE(fit%sigma(n))
    ALLOCATE(fit%best(n))

    ALLOCATE(fit%set(n))
    ALLOCATE(fit%log(n))
    ALLOCATE(fit%cov(n))

    ALLOCATE(fit%name(n))

  END SUBROUTINE allocate_arrays

  SUBROUTINE init_cosmologies(im,cosm,ncos)

    ! Set the number of and types of cosmology depending on imode
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: im
    TYPE(cosmology), ALLOCATABLE, INTENT(OUT) :: cosm(:)
    INTEGER, INTENT(OUT) :: ncos
    INTEGER :: i, icosmo

    ! Set the number of cosmological models
    IF(im==1 .OR. im==2 .OR. im==3 .OR. im==4) THEN

       IF(im==1 .OR. im==3) THEN
          ncos=9 ! Number of Mita Titan nodes - only 10 have Omega_nu = 0. (ignore 1 because it is weird)
       ELSE IF(im==2 .OR. im==4) THEN
          ncos=37 ! Number of FrankenEmu nodes
       ELSE
          STOP 'HMx_FITTING: Error, something went wrong with setting fields'
       END IF

    ELSE IF(im>=11) THEN

       ! Set to the number of different cosmologies
       ncos=1 

    ELSE
       STOP 'HMx_FITTING: Error, something went wrong with setting fields'
    END IF

    ! Allocate arrays for cosmology and fields
    ALLOCATE(cosm(ncos))

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
          STOP 'HMx_FITTING: Error, im mode not specified correctly'
       END IF

       CALL assign_cosmology(icosmo,cosm(i),verbose=.TRUE.)
       CALL init_cosmology(cosm(i))
       CALL print_cosmology(cosm(i))

    END DO

  END SUBROUTINE init_cosmologies

  SUBROUTINE init_fields(im,fields,nf)

    ! Set the number and types of fields to fit to depending on imode
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: im
    INTEGER, ALLOCATABLE, INTENT(OUT) :: fields(:)
    INTEGER, INTENT(OUT) :: nf

    IF(im==1 .OR. im==2 .OR. im==3 .OR. im==4) THEN

       nf=1 ! Matter-matter only

    ELSE IF(im>=11) THEN

       ! Set the number of different fields
       IF(im==28 .OR. im==40) THEN
          nf=5   
       ELSE IF(im==26 .OR. im==29 .OR. im==30 .OR. im==31 .OR. im==32 .OR. &
            im==33 .OR. im==34 .OR. im==35 .OR. im==36 .OR. im==39 .OR. im==41 .OR. &
            im==42 .OR. im==43 .OR. im==44 .OR. im==11 .OR. im==12 .OR. im==14) THEN
          nf=1
       ELSE IF(im==20 .OR. im==21 .OR. im==13) THEN
          nf=2
       ELSE IF(im==16 .OR. im==27 .OR. im==38 .OR. im==46) THEN
          nf=4
       ELSE IF(im==15 .OR. im==37 .OR. im==45) THEN
          nf=3
       ELSE
          STOP 'INIT_FIELDS: Error, something went wrong with setting fields'
       END IF

    ELSE
       STOP 'INIT_FIELDS: Error, something went wrong with setting fields'
    END IF

    ALLOCATE(fields(nf))

    ! Set the fields
    IF(im==1 .OR. im==2 .OR. im==3 .OR. im==4) THEN
       fields(1)=field_matter    
    ELSE IF(im==28 .OR. im==40) THEN
       ! Matter, CDM, gas, stars, electron pressure
       fields(1)=field_matter
       fields(2)=field_cdm
       fields(3)=field_gas
       fields(4)=field_star
       fields(5)=field_electron_pressure     
    ELSE IF(im==16 .OR. im==27 .OR. im==38 .OR. im==46) THEN
       ! Matter, CDM, gas, stars
       fields(1)=field_matter
       fields(2)=field_cdm
       fields(3)=field_gas
       fields(4)=field_star       
    ELSE IF(im==29 .OR. im==39 .OR. im==44 .OR. im==14) THEN
       ! Matter
       fields(1)=field_matter
    ELSE IF(im==26 .OR. im==33 .OR. im==41) THEN
       ! CDM
       fields(1)=field_cdm
    ELSE IF(im==32 .OR. im==34 .OR. im==42 .OR. im==11) THEN
       ! Gas
       fields(1)=field_gas
    ELSE IF(im==31 .OR. im==35 .OR. im==43) THEN
       ! Stars
       fields(1)=field_star
    ELSE IF(im==30 .OR. im==36 .OR. im==12) THEN
       ! Electron pressure
       fields(1)=field_electron_pressure
    ELSE IF(im==15 .OR. im==37 .OR. im==45) THEN
       ! CDM, gas, stars
       fields(1)=field_cdm
       fields(2)=field_gas
       fields(3)=field_star
    ELSE IF(im==20 .OR. im==21 .OR. im==13) THEN
       ! Matter and pressure
       fields(1)=field_matter
       fields(2)=field_electron_pressure
    ELSE
       STOP 'INIT_FIELDS: Error, something went wrong'
    END IF

  END SUBROUTINE init_fields

  SUBROUTINE init_redshifts(im,z,nz)

    ! Set the redshifts to fit over depending on imode
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: im
    REAL, ALLOCATABLE, INTENT(OUT) :: z(:)
    INTEGER, INTENT(OUT) :: nz

    ! Set the number of redshifts
    IF(im==1 .OR. im==2 .OR. im==3 .OR. im==4) THEN

       ! Mira Titan or FrankenEmu
       nz=4

    ELSE IF(im>=11) THEN

       ! Hydro simulations
       IF(name=='') name='AGN_TUNED_nu0'

       IF(im==20 .OR. im==26 .OR. im==27 .OR. im==28 .OR. im==29 .OR. im==30 .OR. im==31 .OR. im==32 .OR. &
            im==33 .OR. im==34 .OR. im==35 .OR. im==36 .OR. im==37 .OR. im==38 .OR. im==39 .OR. &
            im==40 .OR. im==11 .OR. im==12 .OR. im==13 .OR. im==14 .OR. im==15 .OR. im==16) THEN
          nz=1
       ELSE IF(im==21 .OR. im==41 .OR. im==42 .OR. im==43 .OR. im==44 .OR. im==45 .OR. im==46) THEN
          nz=11
       ELSE
          STOP 'HMx_FITTING: Error, mode im not specified correctly for nz'
       END IF

    ELSE

       STOP 'HMx_FITTING: Error, mode im not specified correctly for redshifts'

    END IF

    ! Allocate arrays
    ALLOCATE(z(nz))

    ! Set the actual redshifts
    IF(im==1 .OR. im==2 .OR. im==3 .OR. im==4) THEN

       ! Mira Titan or FrankenEmu
       z(1)=0.0
       z(2)=0.5
       z(3)=1.0
       z(4)=2.0

    ELSE IF(im>=11) THEN      

       ! Set the redshifts
       IF(im==20 .OR. im==26 .OR. im==27 .OR. im==28 .OR. im==29 .OR. im==30 .OR. im==31 .OR. im==32 .OR. &
            im==33 .OR. im==34 .OR. im==35 .OR. im==36 .OR. im==37 .OR. im==38 .OR. im==39 .OR. &
            im==40 .OR. im==11 .OR. im==12 .OR. im==13 .OR. im==14 .OR. im==15 .OR. im==16) THEN        
          IF((zin)=='') THEN
             z(1)=z_default
          ELSE
             READ(zin,*) z(1)
          END IF
       ELSE IF(im==21 .OR. im==41 .OR. im==42 .OR. im==43 .OR. im==44 .OR. im==45 .OR. im==46) THEN
          ! All 11 snapshots between z=0 and z=2 from BAHAMAS
          z(1)=0.000
          z(2)=0.125
          z(3)=0.250
          z(4)=0.375
          z(5)=0.500
          z(6)=0.750
          z(7)=1.000
          z(8)=1.250
          z(9)=1.500
          z(10)=1.750
          z(11)=2.000
       ELSE
          STOP 'HMx_FITTING: Error, mode im not specified correctly for redshifts'
       END IF

    ELSE

       STOP 'HMx_FITTING: Error, mode im not specified correctly for redshifts'

    END IF

  END SUBROUTINE init_redshifts

  SUBROUTINE init_halomods(im,hmod,ncos)

    ! Set the halo models depending on the mode selected by the user
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: im
    TYPE(halomod), ALLOCATABLE, INTENT(OUT) :: hmod(:)
    INTEGER, INTENT(IN) :: ncos
    INTEGER :: i, ihm

    ! Choose halo model type
    IF(im==1 .OR. im==2 .OR. im==3 .OR. im==4) THEN
       ihm=15 ! 15 - HMcode 2018
    ELSE IF(im==20 .OR. im==21 .OR. im==26 .OR. im==27 .OR. im==28 .OR. im==29 .OR. im==30 .OR. im==31 .OR. &
         im==32 .OR. im==33 .OR. im==34 .OR. im==35 .OR. im==36 .OR. im==37 .OR. im==38 .OR. im==39 .OR. &
         im==40 .OR. im==41 .OR. im==42 .OR. im==43 .OR. im==44 .OR. im==45 .OR. im==46) THEN
       !ihm=3  ! 3 - Standard halo-model
       ihm=20 ! 20 - Standard halo model in response
    ELSE IF(im==11 .OR. im==12 .OR. im==13 .OR. im==14 .OR. im==15 .OR. im==16) THEN
       ihm=47 ! Isothermal beta model in response
    END IF

    ! Allocate array for halo models
    ALLOCATE(hmod(ncos))

    ! Set the halo model parameters
    DO i=1,ncos
       CALL assign_halomod(ihm,hmod(i),verbose=.FALSE.)
    END DO

  END SUBROUTINE init_halomods

  SUBROUTINE init_weights(im,weight,ncos,nf,nk,nz)

    ! Set the weights depending on the mode selected by the user
    ! This allows the user to exclude some scales, redshifts, fields combinations etc. etc.
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: im
    REAL, ALLOCATABLE, INTENT(OUT) :: weight(:,:,:,:,:)
    INTEGER, INTENT(IN) :: ncos
    INTEGER, INTENT(IN) :: nf
    INTEGER, INTENT(IN) :: nk
    INTEGER, INTENT(IN) :: nz
    REAL :: kmin, kmax
    INTEGER :: i, j

    ! Standard to give equal weight to everything
    ALLOCATE(weight(ncos,nf,nf,nk,nz))
    weight=1.

    ! No weight to pressure-pressure
    IF(im==20 .OR. im==21 .OR. im==13) THEN
       weight(:,2,2,:,:)=0.
    END IF

    ! No weight to pressure-pressure
    IF(im==28 .OR. im==40) THEN
       weight(:,5,5,:,:)=0.
    END IF

    ! k range for multi-z
    IF(im>=11) THEN

       DO j=1,nz

          kmin=kmin_BAHAMAS
          kmax=kmax_BAHAMAS
!!$        IF(z(j)==0.0) THEN
!!$           kmax=kmax_BAHAMAS_z0p0
!!$        ELSE IF(z(j)==0.5) THEN
!!$           kmax=kmax_BAHAMAS_z0p5
!!$        ELSE IF(z(j)==1.0) THEN
!!$           kmax=kmax_BAHAMAS_z1p0
!!$        ELSE IF(z(j)==2.0) THEN
!!$           kmax=kmax_BAHAMAS_z2p0
!!$        ELSE
!!$           STOP 'HMx_FITTING: Error, something went wrong setting kmax for this z'
!!$        END IF

          DO i=1,nk
             IF(k(i)<kmin .OR. k(i)>kmax) weight(:,:,:,i,j)=0.
          END DO

       END DO

    END IF

  END SUBROUTINE init_weights

  SUBROUTINE init_mode(im,fit)

    ! Chooses which parameters are to be varied based on imode selected by the user
    ! The sole purpose of this routine is to change some of the fit%set(i) from FALSE to TRUE
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: im
    TYPE(fitting), INTENT(INOUT) :: fit

    IF(im==1 .OR. im==2 .OR. im==3 .OR. im==4) THEN

       ! 1,2,3,4 - HMcode
       fit%set(param_HMcode_Dv0)=.TRUE.
       fit%set(param_HMcode_Dvp)=.TRUE.
       fit%set(param_HMcode_dc0)=.TRUE.
       fit%set(param_HMcode_dcp)=.TRUE.
       fit%set(param_HMcode_eta0)=.TRUE.
       fit%set(param_HMcode_eta1)=.TRUE.
       fit%set(param_HMcode_f0)=.TRUE.
       fit%set(param_HMcode_fp)=.TRUE.
       fit%set(param_HMcode_kstar)=.TRUE.
       fit%set(param_HMcode_As)=.TRUE.
       fit%set(param_HMcode_alpha0)=.TRUE.
       fit%set(param_HMcode_alpha1)=.TRUE.

    ELSE IF(im>=11) THEN

       !! CDM !!

       IF(im==26 .OR. im==33) THEN

          ! 26 - fixed z; basic parameterts; CDM
          ! 33 - fixed z; final parameters; CDM
          fit%set(param_eps)=.TRUE.
          !fit%set(param_M0)=.TRUE. ! Does not help

       ELSE IF(im==41) THEN

          ! 41 - extended z; final parameters; CDM
          fit%set(param_eps)=.TRUE.
          fit%set(param_epsz)=.TRUE.

          !! !!

          !! gas !!

       ELSE IF(im==32) THEN

          ! 32 - fixed z; basic parameters; gas
          fit%set(param_eps)=.TRUE.
          fit%set(param_Gamma)=.TRUE.
          fit%set(param_M0)=.TRUE.

       ELSE IF(im==34) THEN

          ! 34 - fixed z; final parameters; gas
          fit%set(param_Gamma)=.TRUE.
          fit%set(param_Gammap)=.TRUE. 
          fit%set(param_M0)=.TRUE.
          fit%set(param_Astar)=.TRUE. ! Needed because it governs how much overall gas there is

       ELSE IF(im==11) THEN

          ! 47 - fixed z; final parameters; gas
          fit%set(param_ibeta)=.TRUE.
          fit%set(param_ibetap)=.TRUE. 
          fit%set(param_M0)=.TRUE.
          fit%set(param_Astar)=.TRUE. ! Needed because it governs how much overall gas there is

       ELSE IF(im==42) THEN

          ! 42 - extended z; final parameters; gas
          fit%set(param_Gamma)=.TRUE.
          fit%set(param_Gammap)=.TRUE.
          fit%set(param_Gammaz)=.TRUE.
          fit%set(param_M0)=.TRUE.
          fit%set(param_M0z)=.TRUE.
          fit%set(param_Astar)=.TRUE. ! Needed because it governs how much overall gas there is
          fit%set(param_Astarz)=.TRUE.

          !! !!

          !! stars !!
          
       ELSE IF(im==31) THEN

          ! 31 - fixed z; basic parameters; stars
          fit%set(param_Astar)=.TRUE.       

       ELSE IF(im==35) THEN

          ! 35 - fixed z; final parameters; stars
          fit%set(param_Astar)=.TRUE.
          fit%set(param_cstar)=.TRUE.       
          fit%set(param_Mstar)=.TRUE.
          !fit%set(param_sstar)=.TRUE.  ! Does not provide a significant improvement
          !fit%set(param_Astarp)=.TRUE. ! Does not provide a significant improvement
          !fit%set(param_cstarp)=.TRUE. ! Does not provide a significant improvement
          !fit%set(param_eta)=.TRUE.    ! Does not provide a significant improvement

          ! None of the commented-out parameters provide a significant improvement
          ! Main problem seems to be model fitting poorly at transition region at high z
          ! In all cases, at high z the one-halo term looks like a power law

          ! I also tried removing the high-mass saturation of the f_* relation (to A_*/3), this did not help
          ! I also tried using the Schneider & Teyssier stellar-density profile. This did not help.

       ELSE IF(im==43) THEN

          ! 43 - extended z; final parameters; stars
          fit%set(param_Astar)=.TRUE.
          fit%set(param_Astarz)=.TRUE.
          fit%set(param_cstar)=.TRUE.
          !fit%set(param_cstarp)=.TRUE. ! Did not improve anything in the z-fixed case
          fit%set(param_cstarz)=.TRUE.
          fit%set(param_Mstar)=.TRUE.
          fit%set(param_Mstarz)=.TRUE.

          !! !!

          !! matter !!

       ELSE IF(im==29) THEN

          ! 29 - fixed z; basic parameters; matter
          fit%set(param_eps)=.TRUE.
          fit%set(param_Gamma)=.TRUE.
          fit%set(param_M0)=.TRUE.
          fit%set(param_Astar)=.TRUE.

       ELSE IF(im==39) THEN

          ! 39 - fixed z; final parameters; matter
          fit%set(param_eps)=.TRUE.
          fit%set(param_Gamma)=.TRUE.
          fit%set(param_Gammap)=.TRUE.
          fit%set(param_M0)=.TRUE.
          fit%set(param_Astar)=.TRUE.
          !fit%set(param_cstar)=.TRUE. ! Does not cause a large enough change in figure of merit      
          fit%set(param_Mstar)=.TRUE.

       ELSE IF(im==14) THEN

          ! 39 - fixed z; final parameters; matter
          fit%set(param_eps)=.TRUE.
          fit%set(param_ibeta)=.TRUE.
          fit%set(param_ibetap)=.TRUE.
          fit%set(param_M0)=.TRUE.
          fit%set(param_Astar)=.TRUE.
          !fit%set(param_cstar)=.TRUE. ! Does not cause a large enough change in figure of merit      
          fit%set(param_Mstar)=.TRUE.

       ELSE IF(im==44) THEN

          ! 44 - extended z; final parameters; matter
          fit%set(param_eps)=.TRUE.
          fit%set(param_epsz)=.TRUE.
          fit%set(param_Gamma)=.TRUE.
          fit%set(param_Gammap)=.TRUE.
          fit%set(param_Gammaz)=.TRUE.
          fit%set(param_M0)=.TRUE.   
          fit%set(param_M0z)=.TRUE.  
          fit%set(param_Astar)=.TRUE.
          fit%set(param_Astarz)=.TRUE.
          !fit%set(param_cstar)=.TRUE.  ! Did not cause a large enough change in figure-of-merit
          !fit%set(param_cstarp)=.TRUE. ! Did not cause a large enough change in figure-of-merit
          !fit%set(param_cstarz)=.TRUE. ! Did not cause a large enough change in figure-of-merit
          !fit%set(param_Mstar)=.TRUE.  ! Did not cause a large enough change in figure-of-merit
          !fit%set(param_Mstarz)=.TRUE. ! Did not cause a large enough change in figure-of-merit

          !! !!

          !! (matter), CDM, gas, stars !! 

       ELSE IF(im==27) THEN

          ! 27 - fixed z; basic parameters; matter, CDM, gas, stars
          fit%set(param_eps)=.TRUE.
          fit%set(param_Gamma)=.TRUE.
          fit%set(param_M0)=.TRUE.
          fit%set(param_Astar)=.TRUE.

       ELSE IF(im==37 .OR. im==38) THEN

          ! 37 - fixed z; final parameters; CDM, gas stars
          ! 38 - fixed z; final parameters; matter, CDM, gas, stars
          fit%set(param_eps)=.TRUE.
          fit%set(param_Gamma)=.TRUE.
          fit%set(param_Gammap)=.TRUE.
          fit%set(param_M0)=.TRUE.
          fit%set(param_Astar)=.TRUE.
          fit%set(param_cstar)=.TRUE.
          fit%set(param_Mstar)=.TRUE.

       ELSE IF(im==15 .OR. im==16) THEN

          ! 15 - fixed z; final parameters; CDM, gas stars (isothermal)
          ! 16 - fixed z; final parameters; matter, CDM, gas stars (isothermal)
          fit%set(param_eps)=.TRUE.
          fit%set(param_ibeta)=.TRUE.
          fit%set(param_ibetap)=.TRUE.
          fit%set(param_M0)=.TRUE.
          fit%set(param_Astar)=.TRUE.
          fit%set(param_cstar)=.TRUE.
          fit%set(param_Mstar)=.TRUE.

       ELSE IF(im==45 .OR. im==46) THEN

          ! 45 - extended z; final parameters; CDM, gas, stars
          ! 46 - extended z; final parameters; matter, CDM, gas, stars
          fit%set(param_eps)=.TRUE.
          fit%set(param_epsz)=.TRUE.
          fit%set(param_Gamma)=.TRUE.
          fit%set(param_Gammap)=.TRUE.
          fit%set(param_Gammaz)=.TRUE.
          fit%set(param_M0)=.TRUE.   
          fit%set(param_M0z)=.TRUE.  
          fit%set(param_Astar)=.TRUE.
          fit%set(param_Astarz)=.TRUE.
          fit%set(param_cstar)=.TRUE.
          fit%set(param_cstarp)=.TRUE.
          fit%set(param_cstarz)=.TRUE.
          fit%set(param_Mstar)=.TRUE.
          fit%set(param_Mstarz)=.TRUE.

          !! !!

          !! electron pressure !!

       ELSE IF(im==30) THEN

          ! 30 - fixed z; basic parameters; electron pressure
          fit%set(param_alpha)=.TRUE.
          fit%set(param_Gamma)=.TRUE.
          fit%set(param_M0)=.TRUE.
          fit%set(param_Twhim)=.TRUE.

       ELSE IF(im==36) THEN

          ! 36 - fixed z; final parameters; electron pressure
          fit%set(param_alpha)=.TRUE.
          !fit%set(param_alphap)=.TRUE. ! Does this make a difference?
          fit%set(param_Gamma)=.TRUE.
          fit%set(param_Gammap)=.TRUE. ! Does this make a difference?
          fit%set(param_M0)=.TRUE.
          fit%set(param_Twhim)=.TRUE.

       ELSE IF(im==12) THEN

          ! 12 - fixed z; final parameters; electron pressure
          fit%set(param_alpha)=.TRUE.
          fit%set(param_ibeta)=.TRUE.
          fit%set(param_ibetap)=.TRUE.
          fit%set(param_M0)=.TRUE.
          fit%set(param_Twhim)=.TRUE.

          !! !!

          !! matter, electron pressure !!

       ELSE IF(im==20) THEN

          ! 20 - fixed z; final parameters; matter, electron pressure
          fit%set(param_alpha)=.TRUE.
          fit%set(param_alphap)=.TRUE.
          fit%set(param_Twhim)=.TRUE.
          fit%set(param_eps)=.TRUE.
          fit%set(param_Gamma)=.TRUE.
          fit%set(param_Gammap)=.TRUE.
          fit%set(param_M0)=.TRUE.
          fit%set(param_Astar)=.TRUE.
          fit%set(param_Mstar)=.TRUE.

       ELSE IF(im==13) THEN

          ! 13 - fixed z; final parameters; matter, electron pressure (isothermal)
          fit%set(param_alpha)=.TRUE.
          fit%set(param_alphap)=.TRUE.
          fit%set(param_Twhim)=.TRUE.
          fit%set(param_eps)=.TRUE.
          fit%set(param_ibeta)=.TRUE.
          fit%set(param_ibetap)=.TRUE.
          fit%set(param_M0)=.TRUE.
          fit%set(param_Astar)=.TRUE.
          fit%set(param_Mstar)=.TRUE.
          
       ELSE IF(im==21) THEN

          ! 21 - extended z; final parameters; matter, electron pressure
          fit%set(param_alpha)=.TRUE.
          fit%set(param_alphap)=.TRUE.
          fit%set(param_alphaz)=.TRUE.
          fit%set(param_Twhim)=.TRUE.
          fit%set(param_Twhimz)=.TRUE.
          fit%set(param_eps)=.TRUE.
          fit%set(param_epsz)=.TRUE.
          fit%set(param_Gamma)=.TRUE.
          fit%set(param_Gammap)=.TRUE.
          fit%set(param_Gammaz)=.TRUE.
          fit%set(param_M0)=.TRUE.
          fit%set(param_M0z)=.TRUE.  
          fit%set(param_Astar)=.TRUE.
          fit%set(param_Astarz)=.TRUE.

          !! !!

          !! matter, CDM, gas, stars, electron pressure !!

       ELSE IF(im==28) THEN

          ! 28 - fixed z; basic parameters; matter, CDM, gas, stars, electron pressure
          fit%set(param_alpha)=.TRUE.
          fit%set(param_eps)=.TRUE.
          fit%set(param_Gamma)=.TRUE.
          fit%set(param_M0)=.TRUE.
          fit%set(param_Astar)=.TRUE.
          fit%set(param_Twhim)=.TRUE.

       ELSE IF(im==40) THEN

          ! 40 - fixed z; final parameters; matter, CDM, gas, stars, electron pressure
          fit%set(param_alpha)=.TRUE.
          fit%set(param_eps)=.TRUE.
          fit%set(param_Gamma)=.TRUE.
          fit%set(param_M0)=.TRUE.
          fit%set(param_Astar)=.TRUE.
          fit%set(param_cstar)=.TRUE.        
          fit%set(param_Mstar)=.TRUE.
          fit%set(param_Twhim)=.TRUE.

          !! !!

       ELSE

          STOP 'INIT_MODE: Something went wrong with setting parameters'

       END IF

    ELSE

       STOP 'INIT_MODE: Something went wrong with setting parameters'

    END IF

  END SUBROUTINE init_mode

  SUBROUTINE read_simulation_power_spectra(im,k,nk,pow,cosm,ncos,z,nz,fields,nf)

    ! Read in the simulation P(k) data that is to be fitted
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: im
    REAL, ALLOCATABLE, INTENT(OUT) :: k(:)
    INTEGER, INTENT(OUT) :: nk
    REAL, ALLOCATABLE, INTENT(OUT) :: pow(:,:,:,:,:)
    TYPE(cosmology), INTENT(INOUT) :: cosm(ncos)
    INTEGER, INTENT(IN) :: ncos
    REAL, INTENT(IN) :: z(nz)
    INTEGER, INTENT(IN) :: nz
    INTEGER, INTENT(IN) :: fields(nf)
    INTEGER, INTENT(IN) :: nf
    INTEGER :: i, j, j1, j2
    INTEGER :: ip(2)
    REAL, ALLOCATABLE :: k_sim(:), pow_sim(:)

    LOGICAL, PARAMETER :: rebin_emu=.TRUE.
    LOGICAL, PARAMETER :: verbose_BAHAMAS=.TRUE.

    ! Loop over cosmologies and redshifts
    DO i=1,ncos
       DO j=1,nz

          ! Loop over fields
          DO j1=1,nf
             DO j2=j1,nf

                ! Read in power spectra
                IF(im==1 .OR. im==3) THEN
                   CALL read_Mira_Titan_power(k_sim,pow_sim,nk,z(j),cosm(i),rebin=rebin_emu)
                ELSE IF(im==2 .OR. im==4) THEN
                   CALL read_FrankenEmu_power(k_sim,pow_sim,nk,z(j),cosm(i),rebin=rebin_emu)
                ELSE IF(im>=11) THEN
                   ip(1)=fields(j1)
                   ip(2)=fields(j2)
                   CALL read_BAHAMAS_power(k_sim,pow_sim,nk,z(j),name,mesh,ip,cosm(i),&
                        kmin_default,kmax_default,cut_nyquist,subtract_shot,response=resp_BAHAMAS,verbose=verbose_BAHAMAS)
                ELSE
                   STOP 'READ_SIMULATION_POWER_SPECTRA: Error, something went wrong reading data'
                END IF

                ! Allocate big arrays for P(k,z,cosm)
                IF(.NOT. ALLOCATED(k))   ALLOCATE(k(nk))
                IF(.NOT. ALLOCATED(pow)) ALLOCATE(pow(ncos,nf,nf,nk,nz))
                k=k_sim
                pow(i,j1,j2,:,j)=pow_sim

             END DO
          END DO

       END DO
    END DO

  END SUBROUTINE read_simulation_power_spectra

  SUBROUTINE set_HMx_parameters(p_in,n,fit,hmod)

    ! This routine is important.
    ! It translates the fitting parameters (p_in) to parameters in the halomod structure.
    ! It should be the ONLY place in the code that these translations are performed
    IMPLICIT NONE
    REAL, INTENT(IN) :: p_in(n)          ! Parameter set to map to HMx
    INTEGER, INTENT(IN) :: n             ! Total number of parameters
    TYPE(fitting), INTENT(INOUT) :: fit  ! Fitting parameters
    TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model
    REAL :: p_out(n)
    INTEGER :: i

    ! Exponentiate those parameters that need it
    DO i=1,n
       IF(fit%log(i)) THEN
          p_out(i)=10**p_in(i)
       ELSE
          p_out(i)=p_in(i)
       END IF
    END DO

    ! HMcode
    IF(fit%set(param_HMcode_Dv0))    hmod%Dv0=p_out(param_HMcode_Dv0)
    IF(fit%set(param_HMcode_Dvp))    hmod%Dv1=p_out(param_HMcode_Dvp)
    IF(fit%set(param_HMcode_dc0))    hmod%dc0=p_out(param_HMcode_dc0)
    IF(fit%set(param_HMcode_dcp))    hmod%dc1=p_out(param_HMcode_dcp)
    IF(fit%set(param_HMcode_eta0))   hmod%eta0=p_out(param_HMcode_eta0)
    IF(fit%set(param_HMcode_eta1))   hmod%eta1=p_out(param_HMcode_eta1)
    IF(fit%set(param_HMcode_f0))     hmod%f0=p_out(param_HMcode_f0)
    IF(fit%set(param_HMcode_fp))     hmod%f1=p_out(param_HMcode_fp)
    IF(fit%set(param_HMcode_kstar))  hmod%ks=p_out(param_HMcode_kstar)
    IF(fit%set(param_HMcode_As))     hmod%As=p_out(param_HMcode_As)
    IF(fit%set(param_HMcode_alpha0)) hmod%alp0=p_out(param_HMcode_alpha0)
    IF(fit%set(param_HMcode_alpha1)) hmod%alp1=p_out(param_HMcode_alpha1)

    ! Hydro
    IF(fit%set(param_alpha))  hmod%alpha=p_out(param_alpha)
    IF(fit%set(param_beta))   hmod%beta=p_out(param_beta)
    IF(fit%set(param_eps))    hmod%eps=p_out(param_eps)
    IF(fit%set(param_Gamma))  hmod%Gamma=1.+p_out(param_Gamma)   ! CARE: Funny, internal parameter is Gamma-1
    IF(fit%set(param_M0))     hmod%M0=10.**p_out(param_M0)       ! CARE: Log, internal parameter is log10(M0)
    IF(fit%set(param_Astar))  hmod%Astar=p_out(param_Astar)
    IF(fit%set(param_Twhim))  hmod%Twhim=10.**p_out(param_Twhim) ! CARE: Log, internal parameter is log10(Twhim)
    IF(fit%set(param_cstar))  hmod%cstar=p_out(param_cstar)
    IF(fit%set(param_fcold))  hmod%fcold=p_out(param_fcold)
    IF(fit%set(param_Mstar))  hmod%Mstar=10.**p_out(param_Mstar) ! CARE: Log, internal parameter is log10(M*)
    IF(fit%set(param_sstar))  hmod%sstar=p_out(param_sstar)
    IF(fit%set(param_fhot))   hmod%fhot=p_out(param_fhot)
    IF(fit%set(param_eta))    hmod%eta=p_out(param_eta)
    IF(fit%set(param_ibeta))  hmod%ibeta=p_out(param_ibeta)

    ! Hydro - mass indices
    IF(fit%set(param_alphap)) hmod%alphap=p_out(param_alphap)
    IF(fit%set(param_betap))  hmod%betap=p_out(param_betap)
    IF(fit%set(param_Gammap)) hmod%Gammap=p_out(param_Gammap)
    IF(fit%set(param_Astarp)) hmod%Astarp=p_out(param_Astarp)
    IF(fit%set(param_cstarp)) hmod%cstarp=p_out(param_cstarp)
    IF(fit%set(param_ibetap)) hmod%ibetap=p_out(param_ibetap)

    ! Hydro - z indices
    IF(fit%set(param_alphaz)) hmod%alphaz=p_out(param_alphaz)
    IF(fit%set(param_betaz))  hmod%betaz=p_out(param_betaz)
    IF(fit%set(param_epsz))   hmod%epsz=p_out(param_epsz)
    IF(fit%set(param_Gammaz)) hmod%Gammaz=p_out(param_Gammaz)
    IF(fit%set(param_M0z))    hmod%M0z=p_out(param_M0z)
    IF(fit%set(param_Astarz)) hmod%Astarz=p_out(param_Astarz)
    IF(fit%set(param_cstarz)) hmod%cstarz=p_out(param_cstarz)
    IF(fit%set(param_mstarz)) hmod%mstarz=p_out(param_mstarz)  
    IF(fit%set(param_Twhimz)) hmod%Twhimz=p_out(param_Twhimz)
    IF(fit%set(param_ibetaz)) hmod%ibetaz=p_out(param_ibetaz)

  END SUBROUTINE set_HMx_parameters

  SUBROUTINE fom_multiple(p,np,fields,nf,fom,k,nk,z,nz,pow_mod,pow_sim,weight,fit,hmod,cosm,ncos)

    USE special_functions

    ! Calculates the halo model power and then figures out the figure of merit between this and the simulated P(k)
    IMPLICIT NONE
    REAL, INTENT(IN) :: p(np)         ! Fitting parameters
    INTEGER, INTENT(IN) :: np         ! Total number of fitting parameters
    INTEGER, INTENT(IN) :: fields(nf) ! Field types
    INTEGER, INTENT(IN) :: nf         ! Number of fields
    REAL, INTENT(OUT) :: fom          ! Output figure of merit
    REAL, INTENT(IN) :: k(nk)         ! Array of k values for comparison data
    INTEGER, INTENT(IN) :: nk         ! Number of k values for comparison data
    REAL, INTENT(IN) :: z(nz)         ! Array of z values for comparison data
    INTEGER, INTENT(IN) :: nz         ! Number of z values
    REAL, INTENT(OUT) :: pow_mod(ncos,nf,nf,nk,nz)  ! Output array for power as a function of k, z, cosmology
    REAL, INTENT(IN) :: pow_sim(ncos,nf,nf,nk,nz)   ! Comparison power as a function of k, z, cosmology
    REAL, INTENT(IN) :: weight(ncos,nf,nf,nk,nz)    ! Weight array
    TYPE(fitting), INTENT(INOUT) :: fit          ! Fitting structure
    TYPE(halomod), INTENT(INOUT) :: hmod(ncos)   ! Array of halo models for each comparison
    TYPE(cosmology), INTENT(INOUT) :: cosm(ncos) ! Array of cosmological models for each comparison
    INTEGER, INTENT(IN) :: ncos                  ! Number of cosmological models being compared
    INTEGER :: icos, iz, if1, if2, ik
    REAL :: pow_li(ncos,nk,nz), pow_2h(ncos,nf,nf,nk,nz), pow_1h(ncos,nf,nf,nk,nz)
    REAL :: neff

    ! Set this counting output variable to zero
    fom=0.
    neff=0.

    ! Set this to zero too, for the banter
    pow_mod=0.

    ! Loop over cosmologies
    DO icos=1,ncos

       CALL set_HMx_parameters(p,np,fit,hmod(icos))

       ! Loop over redshifts
       DO iz=1,nz

          ! Initialise the halo-model calculation
          CALL init_halomod(mmin,mmax,scale_factor_z(z(iz)),hmod(icos),cosm(icos),verbose=.FALSE.)
          CALL print_halomod(hmod(icos),cosm(icos),verbose=.FALSE.)

          ! Calculate the halo-model power spectrum
          CALL calculate_HMx_a(fields,nf,k,nk,&
               pow_li(icos,:,iz),pow_2h(icos,:,:,:,iz),pow_1h(icos,:,:,:,iz),pow_mod(icos,:,:,:,iz),&
               hmod(icos),cosm(icos),verbose=.FALSE.,response=.FALSE.)

          ! Calculate figure of merit and add to total
          DO if1=1,nf
             DO if2=if1,nf
                DO ik=1,nk
                   fom=fom+weight(icos,if1,if2,ik,iz)*(pow_mod(icos,if1,if2,ik,iz)/pow_sim(icos,if1,if2,ik,iz)-1.)**2
                   neff=neff+weight(icos,if1,if2,ik,iz)
                END DO
             END DO
          END DO

       END DO

    END DO

    ! Calculate the final figure-of-merit by dividing by the effective number of data points and sqrt
    fom=sqrt(fom/neff)

  END SUBROUTINE fom_multiple

  SUBROUTINE set_parameter_sigma(p_original,np,delta,fields,nf,k,nk,z,nz,pow_sim,weight,fit,hmod,cosm,ncos,verbose)

    USE interpolate

    ! Calculates the sigma to jump parameters by
    IMPLICIT NONE
    REAL, INTENT(IN) :: p_original(np)            ! Fitting parameters
    INTEGER, INTENT(IN) :: np                     ! Total number of parameters
    REAL, INTENT(IN) :: delta                     ! Changed required in figure-of-merit by changing parameter
    INTEGER, INTENT(IN) :: fields(nf)             ! Array of all the fields
    INTEGER, INTENT(IN) :: nf                     ! Total number of fields
    REAL, INTENT(IN) :: k(nk)                     ! Array of wavenumbers [h/Mpc]
    INTEGER, INTENT(IN) :: nk                     ! Number of wavenumbers
    REAL, INTENT(IN) :: z(nz)                     ! Array of redshifts
    INTEGER, INTENT(IN) :: nz                     ! Number of redshifts
    REAL, INTENT(IN) :: pow_sim(ncos,nf,nf,nk,nz) ! Array of simulation power spectra
    REAL, INTENT(IN) :: weight(ncos,nf,nf,nk,nz)  ! Weight array for the simulation data
    TYPE(fitting), INTENT(INOUT) :: fit           ! Fitting
    TYPE(halomod), INTENT(INOUT) :: hmod(ncos)    ! Halo model to use
    TYPE(cosmology), INTENT(INOUT) :: cosm(ncos)  ! Cosmology to use
    INTEGER, INTENT(IN) :: ncos                   ! Number of halomodels/cosmologies
    LOGICAL, INTENT(IN) :: verbose                ! Verboseness
    REAL, ALLOCATABLE :: multis(:), sigmas(:), dfom(:)
    REAL, ALLOCATABLE :: sigmas_unique(:), dfom_unique(:)
    INTEGER :: i, j, nv, nnm
    REAL :: fom_base, fom_diff, fom, pow(np,nf,nf,nk,nz), dp
    REAL :: p_perturbed(np)
    REAL :: original, perturbed, sigma, ratio
    REAL, ALLOCATABLE :: perturbation(:,:)

    REAL, PARAMETER :: mult_min=1e-6    ! Multiple minimum (What fraction of the way to the max to we start trying?)
    REAL, PARAMETER :: mult_max=1.      ! Multiple maximum (really should be set to unity)
    INTEGER, PARAMETER :: nm=32         ! Number of parameter sigmas to check
    REAL, PARAMETER :: eps=2.0          ! Tolerated error in fom difference when setting range if doing finesse
    LOGICAL, PARAMETER :: check=.FALSE. ! Do we check eps?
    LOGICAL, PARAMETER :: debug=.FALSE. ! Debug?
    INTEGER, PARAMETER :: iorder=1      ! Order for interpolation to find correct dfom (cubic causes issues with non-monotonic)

    IF(verbose) THEN
       WRITE(*,*) 'SET_PARAMETER_SIGMA: Setting parameter step sizes'
    END IF

    ALLOCATE(multis(nm),sigmas(nm),dfom(nm))
    CALL fill_array(log(mult_min),log(mult_max),multis,nm)
    multis=exp(multis)
    
    ! Get the figure of merit for the base set of parameters
    CALL fom_multiple(p_original,np,fields,nf,fom_base,k,nk,z,nz,pow,pow_sim,weight,fit,hmod,cosm,ncos)

    ! Count the number of parameters being varied
    nv=0
    DO i=1,fit%n
       IF(fit%set(i)) nv=nv+1
    END DO

    ! Write to screen
    IF(verbose) THEN      
       WRITE(*,*) 'SET_PARAMETER_SIGMA: Number of varied parameters:', nv
       WRITE(*,*) 'SET_PARAMETER_SIGMA: Total number of parameters:', np
       WRITE(*,*) 'SET_PARAMETER_SIGMA: Number of cosmologies:', ncos
       WRITE(*,*) 'SET_PARAMETER_SIGMA: Number of fields:', nf
       WRITE(*,*) 'SET_PARAMETER_SIGMA: Number of wavenumbers:', nk
       WRITE(*,*) 'SET_PARAMETER_SIGMA: Number of redshifts:', nz
       WRITE(*,*) 'SET_PARAMETER_SIGMA: Fixing sigma to give change in fom:', delta
       WRITE(*,*) 'SET_PARAMETER_SIGMA: Baseline fom:', fom_base
       WRITE(*,*) '================================================================================'
       WRITE(*,*) 'Parameter           Name    Base value     New Value         Sigma         Ratio'
       WRITE(*,*) '================================================================================'
    END IF

    ! Allocate array to save all the perturbed parameter values
    ALLOCATE(perturbation(np,nm))
    perturbation=0.

    ! Loop over parameters
    DO i=1,np

       ! If we are setting the parameter range here then go do it
       IF(fit%set(i)) THEN

          ! Set the range of parameter 'p(i)' over which to take the derivative
          p_perturbed=p_original ! Reset all parameters to their starting values

          IF(fit%cov(i)) THEN

             sigmas=0.
             dfom=0.
             !sigmas=multis*fit%sigma(i)
             dp=fit%maximum(i)-p_original(i)
             sigmas=multis*dp

             ! Loop over number of attempts to find the correct jump
             DO j=1,nm

                ! Perturb parameter i
                p_perturbed(i)=p_original(i)+sigmas(j) 

                ! Check to see that we have not gone outside range
                IF(p_perturbed(i)<fit%minimum(i)) THEN
                   p_perturbed(i)=fit%minimum(i)
                ELSE IF(p_perturbed(i)>fit%maximum(i)) THEN
                   p_perturbed(i)=fit%maximum(i)
                END IF

                perturbation(i,j)=p_perturbed(i)

                ! Get the figure of merit for the updated parameter
                CALL fom_multiple(p_perturbed,np,fields,nf,fom,k,nk,z,nz,pow,pow_sim,weight,fit,hmod,cosm,ncos)

                ! Calculate the change in the figure of merit for this parameter
                dfom(j)=fom-fom_base

                IF(debug) THEN
                   IF(fit%log(i)) THEN
                      original=10.**p_original(i)
                      perturbed=10.**p_perturbed(i)
                   ELSE
                      original=p_original(i)
                      perturbed=p_perturbed(i)
                   END IF
                   WRITE(*,*) 'DEBUG:', j, original, perturbed, abs(dfom(j))
                END IF

             END DO

!!$             IF(all(dfom==0.)) STOP 'SET_PARAM_SIGMA: Error, changing this parameter does not change power spectra'

             IF(.NOT. within_array(delta,abs(dfom),nm)) THEN
                WRITE(*,*)
                WRITE(*,*) 'SET_PARAMETER_SIGMA: Parameter:', i
                WRITE(*,*) 'SET_PARAMETER_SIGMA: Name: ', trim(fit%name(i))
                WRITE(*,*) 'SET_PARAMETER_SIGMA: Desired change in figure of merit:', delta               
                WRITE(*,*) 'SET_PARAMETER_SIGMA: Original sigma:', fit%sigma(i)
                WRITE(*,*) '========================================================================================='
                WRITE(*,*) '    Attempt             Original                 Perturbed                 delta_fom'
                WRITE(*,*) '========================================================================================='
                DO j=1,nm
                   WRITE(*,*) j, p_original(i), perturbation(i,j), abs(dfom(j))
                END DO
                WRITE(*,*) '========================================================================================='
                STOP 'SET_PARAMETER_SIGMA: Error, sigma is outside table range'
             ELSE

                ! This is necessary to remove repeated entries from the dfom table
                ! If this is not done then the fitting routine can go haywire
                ! This happens if we butt up against the parameter boundaries
                ALLOCATE(dfom_unique(nm),sigmas_unique(nm))
                dfom_unique=dfom
                sigmas_unique=sigmas
                CALL remove_repeated_two_array_elements(dfom_unique,sigmas_unique,nm,nnm)

                ! Use find to interpolate to get the sigma that results in the correct dfom
                fit%sigma(i)=exp(find(log(delta),log(abs(dfom_unique)),log(abs(sigmas_unique)),nnm,&
                     iorder=iorder,ifind=3,imeth=2))

                ! Deallocate arrays for unique figure of merit and sigma values
                DEALLOCATE(dfom_unique,sigmas_unique)
                
             END IF

             ! Set the value for the perturbed parameter
             p_perturbed(i)=p_original(i)+fit%sigma(i)

          END IF

          ! Get the figure of merit for the value of perturbed parameter found to give the required change in fom
          ! This is just to check everything is looking sensible
          CALL fom_multiple(p_perturbed,np,fields,nf,fom,k,nk,z,nz,pow,pow_sim,weight,fit,hmod,cosm,ncos)
          fom_diff=fom-fom_base

          ! Write parameters to screen
          IF(verbose) THEN
             IF(fit%log(i)) THEN
                original=10.**p_original(i)
                perturbed=10.**p_perturbed(i)              
             ELSE
                original=p_original(i)
                perturbed=p_perturbed(i)
             END IF
             sigma=perturbed-original
             ratio=abs(fom_diff/delta)
             WRITE(*,fmt='(I10,A15,4F14.7)') i, trim(fit%name(i)), original, perturbed, sigma, ratio
             IF(check .AND. (ratio>eps .OR. ratio<1./eps)) STOP 'SET_PARAMETER_SIGMA: Error, ratio to desired delta not in range'
          END IF

       END IF

    END DO

    ! Write to screen
    IF(verbose) THEN
       WRITE(*,*) '================================================================================'
       WRITE(*,*) 'SET_PARAMETER_SIGMA: Done setting'
       WRITE(*,*)
    END IF

  END SUBROUTINE set_parameter_sigma

  SUBROUTINE write_fitting_power(base,k,z,pow_mod,pow_sim,ncos,nf,nk,nz)

    ! Write fitting data to disk
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: base
    REAL, INTENT(IN) :: k(nk)
    REAL, INTENT(IN) :: z(nz)
    REAL, INTENT(IN) :: pow_mod(ncos,nf,nf,nk,nz)
    REAL, INTENT(IN) :: pow_sim(ncos,nf,nf,nk,nz)
    INTEGER, INTENT(IN) :: ncos
    INTEGER, INTENT(IN) :: nf
    INTEGER, INTENT(IN) :: nk
    INTEGER, INTENT(IN) :: nz
    REAL :: pow_hmcode(nk), a(nz)
    CHARACTER(len=256) :: outfile, outbit
    CHARACTER(len=10) :: uscore, nothing, mid, ext
    INTEGER :: icos, iz, i1, i2, ik

    ! Bits for file name
    uscore='_'
    nothing=''
    mid='_z'
    ext='.dat'

    ! Need 'a' for HMcode
    DO iz=1,nz
       a(iz)=scale_factor_z(z(iz))
    END DO

    ! Loop over everything
    DO icos=1,ncos
       DO iz=1,nz
          CALL calculate_HMcode_a(k,a(iz),pow_hmcode,nk,cosm(icos))
          DO i1=1,nf
             DO i2=i1,nf
                outbit=number_file(base,icos,uscore)
                outbit=number_file2(outbit,i1,nothing,i2,mid)
                outfile=number_file(outbit,iz,ext)
                WRITE(*,*) 'WRITE_FITTING_POWER: Outfile: ', trim(outfile)
                OPEN(7,file=outfile)
                DO ik=1,nk
                   WRITE(7,*) k(ik), pow_mod(icos,i1,i2,ik,iz), pow_sim(icos,i1,i2,ik,iz), pow_hmcode(ik)
                END DO
                CLOSE(7)
                WRITE(*,*) 'WRITE_FITTING_POWER: Done'
                WRITE(*,*)
             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE write_fitting_power

  SUBROUTINE write_model(p,np,fit,hmod,cosm)

    ! Writes out the halo model with parameters p
    IMPLICIT NONE
    REAL, INTENT(IN) :: p(np)              ! Fitting parameters
    INTEGER, INTENT(IN) :: np              ! Total number of fitting parameters
    TYPE(fitting), INTENT(INOUT) :: fit    ! Fitting structure
    TYPE(halomod), INTENT(INOUT) :: hmod   ! Halo model
    TYPE(cosmology), INTENT(INOUT) :: cosm ! Array of cosmological models for each comparison

    CALL set_HMx_parameters(p,np,fit,hmod)

    CALL print_halomod(hmod,cosm,verbose=.TRUE.)

  END SUBROUTINE write_model

  SUBROUTINE actual_fitting(p_start,np,delta,tmax,nchain,k,nk,z,nz,fields,nf,weight,pow_sim,fit,hmod,cosm,ncos)

    ! Does the actual fitting of parameters
    ! TODO: Clean this up, this is the least clean of all the routines here
    IMPLICIT NONE
    REAL, INTENT(IN) :: p_start(np)
    INTEGER, INTENT(IN) :: np
    REAL, INTENT(IN) :: delta
    REAL, INTENT(IN) :: tmax
    INTEGER, INTENT(IN) :: nchain
    REAL, INTENT(IN) :: k(nk)
    INTEGER, INTENT(IN) :: nk
    REAL, INTENT(IN) :: z(nz)
    INTEGER, INTENT(IN) :: nz
    INTEGER, INTENT(IN) :: fields(nf)
    INTEGER, INTENT(IN) :: nf
    REAL, INTENT(IN) :: weight(ncos,nf,nf,nk,nz)
    REAL, INTENT(IN) :: pow_sim(ncos,nf,nf,nk,nz)
    TYPE(fitting), INTENT(INOUT) :: fit
    TYPE(halomod), INTENT(INOUT) :: hmod(ncos)
    TYPE(cosmology), INTENT(INOUT) :: cosm(ncos)
    INTEGER, INTENT(IN) :: ncos
    INTEGER :: i, j, l
    REAL :: p_new(fit%n), p_old(fit%n)
    REAL :: fom_old, fom_new, fom_bst, fom_ori
    REAL :: original, best, maximum, minimum
    INTEGER :: i_bet, i_wor, i_acc, i_fai, i_tot, i_bst
    INTEGER :: out
    REAL :: t1, t2
    LOGICAL :: accept
    REAL, ALLOCATABLE :: pow_mod(:,:,:,:,:)
    CHARACTER(len=256) :: outfile
    LOGICAL :: verbose

    ! Fix these
    p_new=p_start
    p_old=p_start

    ! Set the best figures-of-merit to some huge value
    fom_old=HUGE(fom_old)
    fom_new=HUGE(fom_new)
    fom_bst=HUGE(fom_bst)

    ! Loop over number of runs
    WRITE(*,*) 'HMx_FITTING: Starting fitting'
    WRITE(*,*) 'HMx_FITTING: Number of points in chain:', nchain
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

    ! Allocate arrays for halo-model power
    ALLOCATE(pow_mod(ncos,nf,nf,nk,nz))

    ! Do the chain
    DO l=1,nchain+1

       IF(l==1 .OR. mod(l,m)==0) THEN
          IF(l==1) THEN
             verbose=.TRUE.
          ELSE
             verbose=.TRUE.
          END IF
          IF(set_ranges) THEN
             CALL set_parameter_sigma(p_old,np,delta,fields,nf,k,nk,z,nz,pow_sim,weight,fit,hmod,cosm,ncos,verbose)
          END IF
          IF(l==1) THEN
             outfile=trim(outbase)//'_chain.dat'
             OPEN(10,file=outfile)
          END IF

       END IF

       IF(l==1) THEN
          ! Do nothing
       ELSE IF(l==nchain+1 .OR. (t2-t1)>60.*tmax) THEN
          ! Set to best-fitting parameters on last go
          p_new=fit%best           
       ELSE
          ! Randomly jump parameters
          DO i=1,fit%n
             IF(fit%set(i)) THEN
                p_new(i)=random_Gaussian(p_old(i),fit%sigma(i))
             ELSE
                p_new(i)=p_old(i)
             END IF
             IF(p_new(i)<fit%minimum(i)) p_new(i)=fit%minimum(i)
             IF(p_new(i)>fit%maximum(i)) p_new(i)=fit%maximum(i)
          END DO
       END IF

       ! Calculate the figure-of-merit
       CALL fom_multiple(p_new,np,fields,nf,fom_new,k,nk,z,nz,pow_mod,pow_sim,weight,fit,hmod,cosm,ncos)   

       ! Write original power spectra to disk
       IF(l==1) THEN

          ! Set the original figure-of-merit to the new figure-of-merit
          fom_old=fom_new
          fom_ori=fom_new
          fom_bst=fom_new

          ! Write out the original data
          base=trim(outbase)//'_orig_cos'
          CALL write_fitting_power(base,k,z,pow_mod,pow_sim,ncos,nf,nk,nz)

          accept=.TRUE.

          i_bst=1
          i_bet=i_bet+1
          i_tot=i_tot+1

       ELSE IF(l==nchain+1 .OR. (t2-t1)>60.*tmax) THEN

          WRITE(*,*)

          ! Output the best-fitting model
          base=trim(outbase)//'_best_cos'
          CALL write_fitting_power(base,k,z,pow_mod,pow_sim,ncos,nf,nk,nz)

          accept=.TRUE.
          EXIT

       ELSE

          ! Decide on acceptance and add to tallies
          i_tot=i_tot+1
          IF(fom_new < fom_bst) THEN
             ! If fit is the best then always accept...
             fit%best=p_new
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

       IF(l .NE. nchain+1) THEN
          WRITE(*,fmt='(I10,3F14.7,L3)') l, fom_bst, fom_old, fom_new, accept
          !WRITE(*,*) l, p_old(param_ibeta), p_new(param_ibeta)
       END IF

       IF(accept) THEN
          i_acc=i_acc+1
          p_old=p_new
          fom_old=fom_new
          WRITE(10,*) l, fom_old, (p_old(j), j=1,fit%n)
       END IF

       CALL cpu_time(t2)

    END DO
    CLOSE(10)
    WRITE(*,*) 'HMx_FITTING: Done'
    WRITE(*,*)

    ! Write useful information to screen and file
    DO j=1,2

       IF(j==1) THEN
          ! For writing to screen
          out=6
       ELSE IF(j==2) THEN
          ! For writing to file
          out=77
          outfile=trim(outbase)//'_params.dat'
          OPEN(out,file=outfile)
       ELSE
          STOP 'HMx_FITTING: Error, output fucked up badly'
       END IF

       WRITE(out,*) 'HMx_FITTING: Maximum time [mins]:', tmax
       WRITE(out,*) 'HMx_FITTING: Accuracy parameter:', delta
       WRITE(out,*) 'HMx_FITTING: Requested number of attempts:', nchain
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
       DO i=1,fit%n
          IF(fit%set(i)) THEN
             IF(fit%log(i)) THEN
                ! If the parameter is explored in log space then you need to do this
                original=10.**fit%original(i)
                best=10.**fit%best(i)
                minimum=10.**fit%minimum(i)
                maximum=10.**fit%maximum(i)
             ELSE
                ! Otherwise do nothing
                original=fit%original(i)
                best=fit%best(i)
                minimum=fit%minimum(i)
                maximum=fit%maximum(i)
             END IF
             WRITE(out,fmt='(I10,A15,4F14.7)') i, trim(fit%name(i)), original, best, minimum, maximum
          END IF
       END DO
       WRITE(out,*) '================================================================================'
       WRITE(out,*)

       IF(j==2) THEN
          CLOSE(out)
       END IF

    END DO

!!$    DO i=1,ncos
!!$       CALL write_model(fit%best,np,fit,hmod(i),cosm(i))
!!$    END DO

  END SUBROUTINE actual_fitting

END PROGRAM HMx_fitting
