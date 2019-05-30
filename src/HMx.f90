MODULE HMx

  ! Module usage statements
  USE constants
  USE array_operations
  USE solve_equations
  USE special_functions
  USE string_operations
  USE calculus_table
  USE cosmology_functions
  USE logical_operations
  USE interpolate
  USE file_info

  IMPLICIT NONE

  PRIVATE

  ! Main routines
  PUBLIC :: mmin_HMx
  PUBLIC :: mmax_HMx
  PUBLIC :: halomod ! Type
  PUBLIC :: assign_halomod
  PUBLIC :: init_halomod
  PUBLIC :: print_halomod
  PUBLIC :: calculate_HMx
  PUBLIC :: calculate_HMcode
  PUBLIC :: calculate_HMx_a
  PUBLIC :: calculate_HMcode_a
  PUBLIC :: set_halo_type
  PUBLIC :: halo_type
  PUBLIC :: M_nu
  PUBLIC :: b_nu
  PUBLIC :: g_nu
  PUBLIC :: nu_M
  PUBLIC :: mean_bias_number_weighted
  PUBLIC :: mean_nu_number_weighted
  PUBLIC :: mean_halo_mass_number_weighted
  PUBLIC :: mean_bias_mass_weighted
  PUBLIC :: mean_nu_mass_weighted
  PUBLIC :: mean_halo_mass_mass_weighted  
  PUBLIC :: mean_halo_number_density
  PUBLIC :: virial_radius
  PUBLIC :: convert_mass_definitions
  PUBLIC :: win_type
  PUBLIC :: UPP              ! TODO: Retire
  PUBLIC :: p_1void          ! TODO: Retire
  PUBLIC :: halo_HI_fraction ! TODO: Retire
  PUBLIC :: T_1h
  PUBLIC :: BNL

  ! Diagnostics
  PUBLIC :: halo_definitions
  PUBLIC :: halo_diagnostics
  PUBLIC :: halo_properties
  PUBLIC :: write_halo_profiles
  PUBLIC :: write_mass_fractions

  ! Winint functions
  PUBLIC :: winint_diagnostics
  PUBLIC :: winint_speed_tests

  ! Mass functions and bias
  PUBLIC :: mass_function
  PUBLIC :: multiplicity_function
  PUBLIC :: halo_bias
  
  ! HMx functions
  PUBLIC :: HMx_alpha
  PUBLIC :: HMx_beta
  PUBLIC :: HMx_eps
  PUBLIC :: HMx_Gamma
  PUBLIC :: HMx_M0
  PUBLIC :: HMx_Astar
  PUBLIC :: HMx_Twhim

  ! Fields
  PUBLIC :: field_dmonly
  PUBLIC :: field_matter
  PUBLIC :: field_cdm
  PUBLIC :: field_gas
  PUBLIC :: field_star
  PUBLIC :: field_bound_gas
  PUBLIC :: field_free_gas
  PUBLIC :: field_electron_pressure
  PUBLIC :: field_void
  PUBLIC :: field_compensated_void
  PUBLIC :: field_central_galaxies
  PUBLIC :: field_satellite_galaxies
  PUBLIC :: field_galaxies
  PUBLIC :: field_HI
  PUBLIC :: field_cold_gas
  PUBLIC :: field_hot_gas
  PUBLIC :: field_static_gas
  PUBLIC :: field_central_stars
  PUBLIC :: field_satellite_stars
  PUBLIC :: field_CIB_353
  PUBLIC :: field_CIB_545
  PUBLIC :: field_CIB_857
  PUBLIC :: field_halo_11p0_11p5
  PUBLIC :: field_halo_11p5_12p0
  PUBLIC :: field_halo_12p0_12p5
  PUBLIC :: field_halo_12p5_13p0
  PUBLIC :: field_halo_13p0_13p5
  PUBLIC :: field_halo_13p5_14p0
  PUBLIC :: field_halo_14p0_14p5
  PUBLIC :: field_halo_14p5_15p0
  PUBLIC :: i1_fields
  PUBLIC :: i2_fields

  ! Total number of fitting parameters
  PUBLIC :: param_n

  ! Fitting parameters - hydro
  PUBLIC :: param_alpha
  PUBLIC :: param_beta
  PUBLIC :: param_eps
  PUBLIC :: param_Gamma
  PUBLIC :: param_M0
  PUBLIC :: param_Astar
  PUBLIC :: param_Twhim
  PUBLIC :: param_cstar
  PUBLIC :: param_fcold
  PUBLIC :: param_mstar
  PUBLIC :: param_sstar
  PUBLIC :: param_fhot
  PUBLIC :: param_eta
  PUBLIC :: param_ibeta
  
  ! Fitting parameters - hydro mass indices  
  PUBLIC :: param_alphap
  PUBLIC :: param_betap
  PUBLIC :: param_Gammap
  PUBLIC :: param_Astarp
  PUBLIC :: param_cstarp
  PUBLIC :: param_ibetap

  ! Fitting parameters - hydro z indicies
  PUBLIC :: param_alphaz
  PUBLIC :: param_betaz
  PUBLIC :: param_Gammaz
  PUBLIC :: param_M0z
  PUBLIC :: param_Astarz
  PUBLIC :: param_Twhimz
  PUBLIC :: param_cstarz
  PUBLIC :: param_Mstarz
  PUBLIC :: param_epsz
  PUBLIC :: param_ibetaz

  ! Fitting parameters - HMcode
  PUBLIC :: param_HMcode_Dv0
  PUBLIC :: param_HMcode_Dvp
  PUBLIC :: param_HMcode_dc0
  PUBLIC :: param_HMcode_dcp
  PUBLIC :: param_HMcode_eta0
  PUBLIC :: param_HMcode_eta1
  PUBLIC :: param_HMcode_f0
  PUBLIC :: param_HMcode_fp
  PUBLIC :: param_HMcode_kstar
  PUBLIC :: param_HMcode_As
  PUBLIC :: param_HMcode_alpha0
  PUBLIC :: param_HMcode_alpha1
  !PUBLIC :: param_HMcode_Dvnu
  !PUBLIC :: param_HMcode_dcnu

  ! Halo-model stuff that needs to be recalculated for each new z
  TYPE halomod

     ! Redshift and scale factor
     REAL :: z, a

     ! Switches
     INTEGER :: ip2h, ibias, imf, iconc, iDolag, iAs, ip2h_corr
     INTEGER :: idc, iDv, ieta, ikstar, i2hdamp, i1hdamp, itrans

     ! Void stuff
     LOGICAL :: voids

     ! Spherical collapse parameters
     REAL :: dc, Dv

     ! HMx baryon parameters
     REAL :: alpha, beta, eps, Gamma, M0, Astar, Twhim, ibeta ! HMx baryon parameters
     REAL :: cstar, sstar, mstar, Theat, fcold, fhot, eta ! HMx baryon parameters
     REAL :: alphap, betap, Gammap, cstarp, Astarp, ibetap ! HMx mass-power parameters
     REAL :: alphaz, betaz, epsz, Gammaz, M0z, Astarz, Twhimz, cstarz, mstarz, ibetaz ! HMx z-power parameters
     REAL :: A_alpha, B_alpha, C_alpha, D_alpha, E_alpha
     REAL :: A_eps, B_eps, C_eps, D_eps
     REAL :: A_Gamma, B_Gamma, C_Gamma, D_Gamma, E_gamma
     REAL :: A_M0, B_M0, C_M0, D_M0, E_M0
     REAL :: A_Astar, B_Astar, C_Astar, D_Astar
     REAL :: A_Twhim, B_Twhim, C_Twhim, D_Twhim

     ! Look-up tables
     REAL, ALLOCATABLE :: c(:), rv(:), nu(:), sig(:), zc(:), m(:), rr(:), sigf(:), log_m(:)
     REAL, ALLOCATABLE :: r500(:), m500(:), c500(:), r200(:), m200(:), c200(:)
     REAL, ALLOCATABLE :: r500c(:), m500c(:), c500c(:), r200c(:), m200c(:), c200c(:)
     INTEGER :: n

     ! Window-function (not used?)
     REAL, ALLOCATABLE :: k(:), wk(:,:,:)
     INTEGER :: nk
     
     REAL :: sigv, sigv100, c3, knl, rnl, mnl, neff, sig8z, Rh, Mh, Mp

     ! HOD etc.
     REAL :: mhalo_min, mhalo_max, HImin, HImax, rcore, hmass
     REAL :: n_c, n_s, n_g, rho_HI, dlnc

     ! Mass function integrals
     REAL :: gmin, gmax, gbmin, gbmax, gnorm
     
     ! HMcode parameters
     REAL :: Dv0, Dv1, dc0, dc1, eta0, eta1, f0, f1, ks, As, alp0, alp1, Dvnu, dcnu     

     ! Halo types
     INTEGER :: halo_DMONLY, halo_CDM, halo_static_gas, halo_cold_gas, halo_hot_gas, halo_free_gas
     INTEGER :: halo_central_stars, halo_satellite_stars, halo_HI
     INTEGER :: halo_void, halo_compensated_void, electron_pressure

     ! Halo components
     INTEGER :: frac_central_stars, frac_stars, frac_HI
     INTEGER :: frac_bound_gas, frac_cold_bound_gas, frac_hot_bound_gas
     
     LOGICAL :: one_parameter_baryons
     
     LOGICAL :: has_HI, has_galaxies, has_mass_conversions, safe_negative, has_dewiggle
     
     LOGICAL :: simple_pivot
     INTEGER :: response
     REAL :: acc_HMx, small_nu, large_nu
     CHARACTER(len=256) :: name
     REAL, ALLOCATABLE :: log_k_pdamp(:), log_pdamp(:)
     INTEGER :: n_pdamp, HMx_mode

     ! Mass function and bias parameters
     REAL :: Tinker_alpha, Tinker_beta, Tinker_gamma, Tinker_phi, Tinker_eta
     REAL :: alpha_numu
     REAL :: ST_p, ST_q, ST_A
     LOGICAL :: has_mass_function

     ! Non-linear halo bias parameters
     REAL, ALLOCATABLE :: bnl(:,:,:), k_bnl(:), nu_bnl(:)
     INTEGER :: nk_bnl, nnu_bnl
     LOGICAL :: has_bnl=.FALSE.
     
  END TYPE halomod

  ! Window integration
  REAL, PARAMETER :: acc_win=1e-3           ! Window-function integration accuracy parameter
  INTEGER, PARAMETER :: imeth_win=12        ! Window-function integration method
  INTEGER, PARAMETER :: winint_order=3      ! Window-function integration order
  REAL, PARAMETER :: winint_test_seconds=1. ! Approximately how many seconds should each timing test take
  INTEGER, PARAMETER :: nmeth_win=13        ! Number of different winint methods
  INTEGER, PARAMETER :: nlim_bumps=2        ! Do the bumps approximation after this number of bumps
  LOGICAL, PARAMETER :: winint_exit=.FALSE. ! Exit when the contributions to the integral become small

  ! Halomodel
  REAL, PARAMETER :: mmin_HMx=1e7                    ! Standard minimum mass for integration
  REAL, PARAMETER :: mmax_HMx=1e17                   ! Standard maximum mass for integration
  LOGICAL, PARAMETER :: verbose_convert_mass=.FALSE. ! Verbosity when running the mass-definition-conversion routines
  REAL, PARAMETER :: eps_p2h=1e-3                    ! Tolerance for requal for lower limit of two-halo power

  ! Dolag
  ! I changed this from 10 -> 100 to make it more infinite, I used zinf=10 for Mead (2017)!
  REAL, PARAMETER :: zinf_Dolag=100. ! An approximate infinite z

  ! HMcode
  REAL, PARAMETER :: mmin_HMcode=1e7       ! Minimum mass to consider for one-halo integration
  REAL, PARAMETER :: mmax_HMcode=1e17      ! Maximum mass to consider for one-halo integration
  REAL, PARAMETER :: fdamp_HMcode_min=1e-3 ! Minimum value for f_damp parameter
  REAL, PARAMETER :: fdamp_HMcode_max=0.99 ! Maximum value for f_damp parameter
  REAL, PARAMETER :: alpha_HMcode_min=0.5  ! Minimum value for alpha transition parameter
  REAL, PARAMETER :: alpha_HMcode_max=2.0  ! Maximum value for alpha transition parameter

  ! HMx
  REAL, PARAMETER :: HMx_alpha_min=1e-2 ! Minimum alpha parameter; needs to be set at not zero
  REAL, PARAMETER :: HMx_beta_min=1e-2  ! Minimum alpha parameter; needs to be set at not zero
  REAL, PARAMETER :: HMx_Gamma_min=1.10 ! Minimum polytropic index
  REAL, PARAMETER :: HMx_Gamma_max=2.00 ! Maximum polytropic index
  REAL, PARAMETER :: HMx_Astar_min=1e-4 ! Minimum halo star fraction; needs to be set at not zero

  ! Minimum fraction before haloes are treated as delta functions (dangerous)
  ! TODO: Hydro tests passed if 1e-4, but I worry a bit, also it was a fairly minor speed up (25%)
  ! TODO: Also, not obvious that a constant frac_min is correct because some species have low abundance but we ususally care about perturbations to this abundance
  REAL, PARAMETER :: frac_min=0.

  ! Halo types
  LOGICAL, PARAMETER :: verbose_galaxies=.FALSE. ! Verbosity when doing the galaxies initialisation
  LOGICAL, PARAMETER :: verbose_HI=.TRUE.       ! Verbosity when doing the HI initialisation

  ! Field types
  ! TODO: Have this run from 1->n, rather than -1->n
  INTEGER, PARAMETER :: field_dmonly=-1
  INTEGER, PARAMETER :: field_matter=0
  INTEGER, PARAMETER :: field_cdm=1
  INTEGER, PARAMETER :: field_gas=2
  INTEGER, PARAMETER :: field_star=3
  INTEGER, PARAMETER :: field_bound_gas=4
  INTEGER, PARAMETER :: field_free_gas=5
  INTEGER, PARAMETER :: field_electron_pressure=6
  INTEGER, PARAMETER :: field_void=7
  INTEGER, PARAMETER :: field_compensated_void=8
  INTEGER, PARAMETER :: field_central_galaxies=9
  INTEGER, PARAMETER :: field_satellite_galaxies=10
  INTEGER, PARAMETER :: field_galaxies=11
  INTEGER, PARAMETER :: field_HI=12
  INTEGER, PARAMETER :: field_cold_gas=13
  INTEGER, PARAMETER :: field_hot_gas=14
  INTEGER, PARAMETER :: field_static_gas=15
  INTEGER, PARAMETER :: field_central_stars=16
  INTEGER, PARAMETER :: field_satellite_stars=17
  INTEGER, PARAMETER :: field_CIB_353=18
  INTEGER, PARAMETER :: field_CIB_545=19
  INTEGER, PARAMETER :: field_CIB_857=20
  INTEGER, PARAMETER :: field_halo_11p0_11p5=21
  INTEGER, PARAMETER :: field_halo_11p5_12p0=22
  INTEGER, PARAMETER :: field_halo_12p0_12p5=23
  INTEGER, PARAMETER :: field_halo_12p5_13p0=24
  INTEGER, PARAMETER :: field_halo_13p0_13p5=25
  INTEGER, PARAMETER :: field_halo_13p5_14p0=26
  INTEGER, PARAMETER :: field_halo_14p0_14p5=27
  INTEGER, PARAMETER :: field_halo_14p5_15p0=28
  INTEGER, PARAMETER :: i1_fields=-1
  INTEGER, PARAMETER :: i2_fields=28

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
  INTEGER, PARAMETER :: param_epsz=33
  INTEGER, PARAMETER :: param_beta=34
  INTEGER, PARAMETER :: param_betap=35
  INTEGER, PARAMETER :: param_betaz=36
  INTEGER, PARAMETER :: param_Astarp=37
  INTEGER, PARAMETER :: param_cstarz=38
  INTEGER, PARAMETER :: param_mstarz=39
  INTEGER, PARAMETER :: param_ibeta=40
  INTEGER, PARAMETER :: param_ibetap=41
  INTEGER, PARAMETER :: param_ibetaz=42
  INTEGER, PARAMETER :: param_n=42

CONTAINS

  SUBROUTINE assign_halomod(ihm,hmod,verbose)

    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: ihm
    TYPE(halomod), INTENT(OUT) :: hmod
    LOGICAL, INTENT(IN) :: verbose
    INTEGER :: i

    ! Names of pre-defined halo models
    INTEGER, PARAMETER :: nhalomod=50 ! Total number of pre-defined halo-model types (TODO: this is stupid)
    CHARACTER(len=256):: names(nhalomod)    
    names(1)='HMcode (Mead et al. 2016)'
    names(2)='Basic halo-model (Two-halo term is linear)'
    names(3)='Standard halo-model (Seljak 2000)'
    names(4)='Standard halo-model but with Mead et al. (2015) smoothed transition'
    names(5)='Standard halo-model but with Delta_v=200 and delta_c=1.686 and Bullock c(M)'
    names(6)='Half-accurate HMcode (Mead et al. 2015, 2016)'
    names(7)='HMcode (Mead et al. 2015)'
    names(8)='Including scatter in halo properties at fixed mass'
    names(9)='Parameters for CCL tests (high accuracy)'
    names(10)='Comparison of mass conversions with Wayne Hu code'
    names(11)='UPP for electron pressure'
    names(12)='Spherical collapse used for Mead (2017) results'
    names(13)='Experimental sigmoid transition'
    names(14)='Experimental scale-dependent halo bias'
    names(15)='HMcode (Mead et al. in prep)'
    names(16)='Halo-void model'
    names(17)='Tilman HMx - AGN 7.6'
    names(18)='Tilman HMx - AGN tuned'
    names(19)='Tilman HMx - AGN 8.0'
    names(20)='Standard halo-model (Seljak 2000) in response'
    names(21)='Cored profile model'
    names(22)='Delta function-NFW star profile model response'
    names(23)='Tinker mass function and bias; virial mass'
    names(24)='Non-linear halo bias for M200c haloes with Tinker'
    names(25)='Villaescusa-Navarro HI halo model'
    names(26)='Delta-function mass function'
    names(27)='Press & Schecter mass function'
    names(28)='One-parameter baryon test'
    names(29)='Adding in cold gas'
    names(30)='Adding in hot gas'
    names(31)='HMcode (2016) with damped BAO'
    names(32)='HMx: AGN 7.6; fixed z; simple pivot'
    names(33)='HMx: AGN tuned; fixed z; simple pivot'
    names(34)='HMx: AGN 8.0; fixed z; simple pivot'
    names(35)='HMx: AGN 7.6; fixed z; clever pivot; f_hot'
    names(36)='HMx: AGN tuned; fixed z; clever pivot; f_hot'
    names(37)='HMx: AGN 8.0; fixed z; clever pivot; f_hot'
    names(38)='HMx: AGN 7.6'
    names(39)='HMx: AGN tuned'
    names(40)='HMx: AGN 8.0'
    names(41)='Put some galaxy mass in the halo/satellites'
    names(42)='Tinker mass function and bias; M200c'
    names(43)='Standard halo-model (Seljak 2000) in matter response'
    names(44)='Tinker mass function and bias; M200'
    names(45)='No stars'
    names(46)='Isothermal beta model for gas'
    names(47)='Isothermal beta model for gas in response'
    names(48)='Non-linear halo bias for standard model'
    names(49)='Non-linear halo bias with Tinker but virial haloes'
    names(50)='HMcode (Mead et al. 2016) with Dolag pow=1 bug'

    IF(verbose) WRITE(*,*) 'ASSIGN_HALOMOD: Assigning halo model'

    ! Default options

    ! Number of points in integration (128 is okay, 1024 is better)
    hmod%n=128

    ! Accuracy for continuous integrals (1e-3 is okay, 1e-4 is better)
    hmod%acc_HMx=1e-3

    ! Small and large values for nu (6 is okay, corrections are suppressed by exp(-large_nu^2)
    hmod%small_nu=1e-6
    hmod%large_nu=6.

    ! Two-halo term
    ! 1 - Linear theory
    ! 2 - Standard from Seljak (2000)
    ! 3 - Linear theory with damped wiggles
    ! 4 - No two-halo term
    hmod%ip2h=2

    ! Method to correct the two-halo integral
    ! NB. This cannot be a parameter here because the value needs to be changed if doing cumulative distributions of power with mass
    ! 1 - Do nothing
    ! 2 - Add value of missing integral assuming that W(k)=1
    ! 3 - Put the missing part of the integrand as a delta function at lower mass limit
    hmod%ip2h_corr=3

    ! Order of halo bias to go to
    ! 1 - Linear order (standard)
    ! 2 - Second order
    ! 3 - Full non-linear halo bias
    ! 4 - Fudge from Fedeli (2014b)
    ! 5 - Fudge from me
    hmod%ibias=1

    ! One-halo term large-scale damping
    ! 1 - No damping
    ! 2 - Mead et al. (2015)
    ! 3 - k^4 at large scales
    hmod%i1hdamp=1

    ! Mass and halo bias function pair
    ! 1 - Press & Schecter (1974)
    ! 2 - Sheth & Tormen (1999)
    ! 3 - Tinker et al. (2010)
    ! 4 - Delta function in mass
    hmod%imf=2

    ! Concentration-mass relation
    ! 1 - Full Bullock et al. (2001; astro-ph/9909159)
    ! 2 - Simple Bullock et al. (2001; astro-ph/9909159)
    ! 3 - Duffy et al. (2008; astro-ph/0804.2486): full 200m
    ! 4 - Duffy et al. (2008; astro-ph/0804.2486): full virial
    ! 5 - Duffy et al. (2008; astro-ph/0804.2486): full 200c
    ! 6 - Duffy et al. (2008; astro-ph/0804.2486): relaxed 200m
    ! 7 - Duffy et al. (2008; astro-ph/0804.2486): relaxed virial
    ! 8 - Duffy et al. (2008; astro-ph/0804.2486): relaxed 200c
    hmod%iconc=4

    ! Linear collapse threshold delta_c
    ! 1 - Fixed 1.686
    ! 2 - Nakamura & Suto (1997) fitting function
    ! 3 - Mead et al. (2016)
    ! 4 - Mead (2017) fitting function
    ! 5 - Spherical-collapse calculation
    ! 6 - Mead et al. (2015)
    hmod%idc=2

    ! Virial density Delta_v
    ! 1 - Fixed 200
    ! 2 - Bryan & Norman (1998; arXiv:astro-ph/9710107) fitting function
    ! 3 - Mead et al. (2015)
    ! 4 - Mead (2017) fitting function
    ! 5 - Spherical-collapse calculation
    ! 6 - Fixed to unity to give Lagrangian radius
    hmod%iDv=2

    ! eta for halo window function
    ! 1 - No
    ! 2 - Mead et al. (2015)
    hmod%ieta=1

    ! k* for one-halo term large-scale damping
    ! 1 - No
    ! 2 - Mead et al. (2015)
    hmod%ikstar=1

    ! Concentration-mass rescaling
    ! 1 - No
    ! 2 - Mead et al. (2015, 2016)
    hmod%iAs=1

    ! fdamp for two-halo term damping
    ! 1 - No
    ! 2 - Mead et al. (2015)
    ! 3 - Mead et al. (2016)
    hmod%i2hdamp=1

    ! alpha for two- to one-halo transition region
    ! 1 - No
    ! 2 - Smoothed transition with alpha
    ! 3 -
    ! 4 - New HMx transition
    ! 5 - Tanh transition
    hmod%itrans=1

    ! Use the Dolag c(M) correction for dark energy?
    ! 1 - No
    ! 2 - Exactly as in Dolag et al. (2004)
    ! 3 - As in Dolag et al. (2004) but with a ^1.5 power
    ! 4 - My version of Dolag (2004) with some redshift dependence
    hmod%iDolag=2

    ! Halo gas fraction
    ! 1 - Fedeli (2014a) bound gas model
    ! 2 - Schneider & Teyssier (2015) bound gas
    ! 3 - Universal baryon fraction
    hmod%frac_bound_gas=2

    ! Halo cold gas fraction
    ! 1 - Constant fraction of bound halo gas
    hmod%frac_cold_bound_gas=1

    ! Halo hot gas fraction
    ! 1 - Constant fraction of bound halo gas
    hmod%frac_hot_bound_gas=1

    ! Halo central star fraction
    ! 1 - Fedeli (2014)
    ! 2 - Constant stellar fraction
    ! 3 - Fedeli (2014) but saturates at high halo mass
    ! 4 - No stars
    hmod%frac_stars=3

    ! Halo central star fraction
    ! 1 - All central
    ! 2 - Schneider et al. (2018)
    hmod%frac_central_stars=2

    ! Halo HI fraction
    ! 1 - Simple model
    ! 2 - Villaescusa-Navarro et al. (2018; 1804.09180)
    hmod%frac_HI=1

    ! DMONLY halo profile
    ! 1 - Analyical NFW
    ! 2 - Non-analytical NFW (good for testing W(k) functions)
    ! 3 - Tophat
    ! 4 - Delta function
    ! 5 - Cored NFW
    ! 6 - Isothermal
    ! 7 - Shell
    hmod%halo_DMONLY=1

    ! CDM halo profile
    ! 1 - NFW
    hmod%halo_CDM=1

    ! Bound static gas halo profile
    ! 1 - Simplified Komatsu & Seljak (2001) gas model
    ! 2 - Isothermal beta model
    ! 3 - Full Komatsu & Seljak (2001) gas model
    hmod%halo_static_gas=1

    ! Bound cold gas halo profile
    ! 1 - Delta function
    hmod%halo_cold_gas=1

    ! Bound hot gas halo profile
    ! 1 - Isothermal
    hmod%halo_hot_gas=1

    ! Free gas halo profile
    ! 1 - Isothermal model (out to 2rv)
    ! 2 - Ejected gas model from Schneider & Teyssier (2015)
    ! 3 - Isothermal shell that connects electron pressure and density to bound gas at rv
    ! 4 - Komatsu-Seljak continuation
    ! 5 - Power-law continuation
    ! 6 - Cubic profile
    ! 7 - Smoothly distributed
    ! 8 - Delta function
    hmod%halo_free_gas=7

    ! Cenrtal stars halo profile
    ! 1 - Fedeli (2014) stellar distribution
    ! 2 - Schneider & Teyssier (2015) stellar distribution
    ! 3 - Delta function
    ! 4 - Delta function at low mass and NFW at high mass
    hmod%halo_central_stars=1

    ! Satellite stars halo profile
    ! 1 - NFW
    ! 2 - Fedeli (2014) stellar distribution
    hmod%halo_satellite_stars=1

    ! HI halo profile
    ! 1 - NFW
    ! 2 - Delta function
    ! 3 - Polynomial with internal exponential hole
    ! 4 - NFW with internal exponential hole
    ! 5 - Modified NFW
    hmod%halo_HI=1

    ! Electron pressure
    ! 1 - UPP (Arnaud et al. 2010)
    ! 2 - From gas profiles
    hmod%electron_pressure=2

    ! Do voids?
    hmod%voids=.FALSE.

    ! Set the void model
    ! 1 - Top-hat void
    hmod%halo_void=1

    ! Set the compensated void model
    ! 1 - Top-hat void
    hmod%halo_compensated_void=1

    ! Safeguard against negative terms in cross correlations
    hmod%safe_negative=.FALSE.

    ! HMcode parameters
    hmod%Dv0=418.
    hmod%Dv1=-0.352
    hmod%dc0=1.59
    hmod%dc1=0.0314
    hmod%eta0=0.603
    hmod%eta1=0.300
    hmod%f0=0.0095
    hmod%f1=1.37
    hmod%ks=0.584
    hmod%As=3.13
    hmod%alp0=3.24
    hmod%alp1=1.85
    hmod%Dvnu=0.916
    hmod%dcnu=0.262

    ! Baryon model (relates eta0 and A)
    hmod%one_parameter_baryons=.FALSE.

    ! HMx parameters
    ! 1 - Fixed parameters
    ! 2 - Mass dependent parameters
    ! 3 - Mass and redshift dependent parameters
    ! 4 - Tilman model
    hmod%HMx_mode=3

    ! Pivot is 1e14 or M_h
    hmod%simple_pivot=.FALSE.

    ! Fixed parameters
    hmod%alpha=1.0     ! Non-virial temperature correction for static gas
    hmod%beta=1.0      ! Non-virial temperature correction for hot gas
    hmod%eps=1.        ! Concentration modification
    hmod%Gamma=1.17    ! Polytropic gas index
    hmod%M0=1e14       ! Halo mass that has lost half gas
    hmod%Astar=0.03    ! Maximum star-formation efficiency
    hmod%Twhim=10**6.5 ! WHIM temperature [K]
    hmod%cstar=10.     ! Stellar concentration r_* = rv/c
    hmod%sstar=1.2     ! sigma_* for f_* distribution
    hmod%Mstar=10**12.5! M* for most efficient halo mass for star formation
    hmod%fcold=0.0     ! Fraction of bound gas that is cold
    hmod%fhot=0.0      ! Fraction of bound gas that is hot
    hmod%eta=0.0       ! Power-law for central galaxy mass fraction
    hmod%ibeta=2./3.   ! Isothermal beta power index
    
    ! Mass indices
    hmod%alphap=0.0    ! Power-law index of alpha with halo mass
    hmod%betap=0.0     ! Power-law index of beta with halo mass
    hmod%Gammap=0.0    ! Power-law index of Gamma with halo mass
    hmod%cstarp=0.0    ! Power-law index of c* with halo mass
    hmod%Astarp=0.0    ! Power-law index of A* with halo mass
    hmod%ibetap=0.0    ! Power-law index of isothermal beta index with halo mass

    ! Redshift indices
    hmod%alphaz=0.0    ! Power-law index of alpha with redshift
    hmod%betaz=0.0     ! Power-law index of beta with redshift
    hmod%epsz=0.0      ! Power-law index of eps with redshift
    hmod%Gammaz=0.0    ! Power-law index of Gamma with redshift
    hmod%M0z=0.0       ! Power-law index of log(M0) with redshift
    hmod%Astarz=0.0    ! Power-law index of Astar with redshift
    hmod%Twhimz=0.0    ! Power-law index of log(Twhim) with redshift
    hmod%cstarz=0.0    ! Power-law index of c* with redshift
    hmod%mstarz=0.0    ! Power-law index of log(M*) with redshift
    hmod%ibetaz=0.0    ! Power-law index of isothermal beta index with redshift

    ! $\alpha$ z and Theat variation
    hmod%A_alpha=-0.005
    hmod%B_alpha=0.022
    hmod%C_alpha=0.865
    hmod%D_alpha=-13.565
    hmod%E_alpha=52.516

    ! $\log_{10} \epsilon_c$ z and Theat variation
    hmod%A_eps=-0.289
    hmod%B_eps=2.147
    hmod%C_eps=0.129
    hmod%D_eps=-0.867

    ! $\Gamma$ z and Theat variation
    hmod%A_Gamma=0.026
    hmod%B_Gamma=-0.064
    hmod%C_Gamma=1.150
    hmod%D_Gamma=-17.011
    hmod%E_Gamma=66.289

    ! $\log_{10}M_0$ z and Theat variation
    hmod%A_M0=-0.007
    hmod%B_M0=0.018
    hmod%C_M0=6.788
    hmod%D_M0=-103.453
    hmod%E_M0=406.705

    ! $A_*$ z and Theat variation
    hmod%A_Astar=0.004
    hmod%B_Astar=-0.034
    hmod%C_Astar=-0.017
    hmod%D_Astar=0.165

    ! $\log_{10}T_{WHIM}$ z and Theat variation
    hmod%A_Twhim=-0.024
    hmod%B_Twhim=0.077
    hmod%C_Twhim=0.454
    hmod%D_Twhim=1.717

    ! Do we treat the halomodel as a response model (multiply by HMcode) or not
    ! 0 - No
    ! 1 - Yes, to all spectra
    ! 2 - Yes, only to matter spectra
    hmod%response=0  

    ! Halo mass if the mass function is a delta function
    hmod%hmass=1e13

    ! Scatter in ln(c)
    hmod%dlnc=0.

    ! HOD parameters
    hmod%mhalo_min=1e12
    hmod%mhalo_max=1e16

    ! HI parameters
    hmod%HImin=1e9
    hmod%HImax=1e12

    ! NFW core radius
    hmod%rcore=0.1

    ! Sheth & Tormen (1999) mass function parameters
    hmod%ST_p=0.3
    hmod%ST_q=0.707

    ! Index to make the mass function integration easier
    ! Should be related to how the mass function diverges at low nu
    ! e.g., ST diverges like nu^-0.6, alpha=1/(1-0.6)=2.5
    hmod%alpha_numu=2.5

    ! Set flags to false
    hmod%has_galaxies=.FALSE.
    hmod%has_HI=.FALSE.
    hmod%has_mass_conversions=.FALSE.
    hmod%has_dewiggle=.FALSE.
    hmod%has_mass_function=.FALSE.
    hmod%has_bnl=.FALSE.

    IF(ihm==-1) THEN
       WRITE(*,*) 'ASSIGN_HALOMOD: Choose your halo model'
       DO i=1,nhalomod
          WRITE(*,*) i, TRIM(names(i))
       END DO
       READ(*,*) ihm
       WRITE(*,*)
    END IF

    IF(ihm==1 .OR. ihm==7 .OR. ihm==15 .OR. ihm==28 .OR. ihm==31 .OR. ihm==50) THEN
       !  1 - HMcode (Mead et al. 2016)
       !  7 - HMcode (Mead et al. 2015)
       ! 15 - HMcode (Mead et al. 2018)
       ! 28 - HMcode (Mead et al. 2016 w/ one parameter baryon model)
       ! 31 - HMcode (Mead et al. 2016 w/ additional BAO damping)
       ! 50 - HMcode (Mead et al. 2016 w/ pow=1 bug in Dolag)
       hmod%ip2h=1
       hmod%i1hdamp=2
       hmod%iconc=1
       hmod%idc=3
       hmod%iDv=3
       hmod%ieta=2
       hmod%ikstar=2
       hmod%iAs=2
       hmod%i2hdamp=3
       hmod%itrans=2
       hmod%iDolag=3
       IF(ihm==7) THEN
          ! Mead et al. (2015)
          hmod%idc=6
          hmod%i2hdamp=2
          hmod%itrans=2
          hmod%f0=0.188
          hmod%f1=4.29
          hmod%alp0=2.93
          hmod%alp1=1.77
          hmod%iDolag=2
          hmod%Dvnu=0.
          hmod%dcnu=0.
       ELSE IF(ihm==15) THEN
          ! Mead et al. (2018)
          hmod%i1hdamp=3 ! k^4 at large scales for one-halo term
          hmod%ip2h=3    ! Linear theory with damped wiggles
          hmod%i2hdamp=2 ! Change back to Mead (2015) model for two-halo damping
          hmod%Dv0=444.87
          hmod%Dv1=-0.3170
          hmod%dc0=1.6077
          hmod%dc1=0.01461
          hmod%eta0=0.5739
          hmod%eta1=0.2705
          hmod%f0=0.08974
          hmod%f1=3.698
          hmod%ks=0.6358
          hmod%As=3.0745
          hmod%alp0=3.129
          hmod%alp1=1.850         
       ELSE IF(ihm==28) THEN
          ! One-parameter baryon model
          hmod%one_parameter_baryons=.TRUE.
          hmod%As=3.13!*4.
       ELSE IF(ihm==31) THEN
          ! Damped BAO
          hmod%ip2h=3
       ELSE IF(ihm==50) THEN
          ! Bug present in CAMB when Dolag power is set to 1 rather than 1.5
          hmod%iDolag=2
       END IF
    ELSE IF(ihm==2) THEN
       ! Basic halo model with linear two halo term (Delta_v = 200, delta_c = 1.686))
       hmod%ip2h=1
       hmod%idc=1
       hmod%iDv=1
       hmod%iconc=1
    ELSE IF(ihm==3) THEN
       ! Standard halo-model calculation (Seljak 2000)
       ! This is the default, so do nothing here
    ELSE IF(ihm==4) THEN
       ! Standard halo-model calculation but with Mead et al. (2015) smoothed two- to one-halo transition and one-halo damping
       hmod%itrans=3
       hmod%ikstar=2
       hmod%i1hdamp=3
       hmod%safe_negative=.TRUE.
    ELSE IF(ihm==5) THEN
       ! Standard halo-model calculation but with Delta_v = 200 and delta_c = 1.686 fixed and Bullock c(M)
       hmod%idc=1
       hmod%iDv=1
       hmod%iconc=1
    ELSE IF(ihm==6) THEN
       ! Half-accurate halo-model calculation, inspired by (Mead et al. 2015, 2016)
       hmod%ikstar=2
       hmod%i1hdamp=3
       hmod%i2hdamp=3
       hmod%safe_negative=.TRUE.
       hmod%idc=3
       hmod%iDv=3
       hmod%ieta=2
       hmod%iAs=2
       hmod%itrans=2
       hmod%iconc=1
       hmod%iDolag=3
    ELSE IF(ihm==8) THEN
       ! Include scatter in halo properties
       hmod%dlnc=0.25
    ELSE IF(ihm==9) THEN
       ! For CCL comparison and benchmark generation
       hmod%n=2048 ! Increase accuracy for the CCL benchmarks
       hmod%acc_HMx=1e-5 ! Increase accuracy for the CCL benchmarks
       hmod%ip2h=2
       hmod%ip2h_corr=3
       hmod%ibias=1
       hmod%i1hdamp=1
       hmod%imf=2
       hmod%iconc=4 ! Virial Duffy relation for full sample
       hmod%idc=2   ! Virial dc
       hmod%iDv=2   ! Virial Dv
       hmod%ieta=1
       hmod%ikstar=1
       hmod%iAs=1
       hmod%i2hdamp=1
       hmod%itrans=1
       hmod%iDolag=1
    ELSE IF(ihm==10) THEN
       ! For mass conversions comparison with Wayne Hu's code
       hmod%iconc=2
       hmod%idc=1
       hmod%iDv=1
    ELSE IF(ihm==11) THEN
       ! UPP
       hmod%electron_pressure=1
    ELSE IF(ihm==12) THEN
       ! Spherical-collapse model to produce Mead (2017) results
       hmod%iconc=1
       hmod%idc=5
       hmod%iDv=5      
       hmod%imf=2    ! Sheth & Tormen mass function
       hmod%iconc=1  ! Bullock et al. c(M) relation
       hmod%iDolag=2 ! This is important for the accuracy of the z=0 results presented in Mead (2017)
    ELSE IF(ihm==13) THEN
       ! Experimental sigmoid transition
       hmod%itrans=5
    ELSE IF(ihm==14) THEN
       ! Experimental scale-dependent halo bias
       hmod%ibias=5
       hmod%ikstar=2
       hmod%i1hdamp=3
    ELSE IF(ihm==16) THEN
       ! Halo-void model
       hmod%voids=.TRUE.
    ELSE IF(ihm==17 .OR. ihm==18 .OR. ihm==19) THEN
       ! 17 - HMx AGN 7.6
       ! 18 - HMx AGN tuned
       ! 19 - HMx AGN 8.0
       hmod%HMx_mode=4
       hmod%itrans=4
       hmod%ikstar=2
       hmod%i1hdamp=3
       hmod%safe_negative=.TRUE.
       hmod%response=1
       IF(ihm==17) THEN
          ! AGN 7.6
          hmod%Theat=10**7.6
       ELSE IF(ihm==18) THEN
          ! AGN tuned
          hmod%Theat=10**7.8
       ELSE IF(ihm==19) THEN
          ! AGN 8.0
          hmod%Theat=10**8.0
       END IF
    ELSE IF(ihm==20) THEN
       ! Standard halo model but as response with HMcode
       hmod%response=1
    ELSE IF(ihm==21) THEN
       ! Cored NFW halo profile model
       hmod%halo_DMONLY=5 ! Cored profile
    ELSE IF(ihm==22) THEN
       ! Different stellar profile
       hmod%halo_central_stars=2 ! Schneider & Teyssier (2015)
       hmod%response=1
    ELSE IF(ihm==23) THEN
       ! Tinker (2010) mass function and bias; virial halo mass
       hmod%imf=3 ! Tinker mass function and bias
    ELSE IF(ihm==24) THEN
       ! Non-linear halo bias for M200c haloes
       hmod%ibias=3   ! Non-linear halo bias
       hmod%iDv=7     ! M200c
       hmod%imf=3     ! Tinker mass function and bias
       hmod%iconc=5   ! Duffy M200c concentrations for full sample
       hmod%idc=1     ! Fixed to 1.686
       !hmod%i1hdamp=3 ! One-halo damping like k^4
       !hmod%ikstar=2  ! One-halo damping via k* from Mead et al. (2015)
    ELSE IF(ihm==25) THEN
       ! Villaescusa-Navarro HI halo model
       hmod%imf=3     ! Tinker mass function
       hmod%frac_HI=3 ! HI mass fraction with z evolution (Villaescusa-Navarro et al. 1804.09180)
       !hmod%halo_HI=3 ! Exponentially cored polynomial (Villaescusa-Navarro et al. 1804.09180)
       !hmod%halo_HI=2 ! Delta function
       hmod%halo_HI=5 ! Modified NFW (Padmanabhan & Refregier 2017; 1607.01021)
       hmod%ibias=3
    ELSE IF(ihm==26) THEN
       ! Delta function mass function
       hmod%ip2h=2        ! Standard two-halo term
       hmod%iDv=6         ! Lagrangian radius haloes
       hmod%imf=4         ! Delta function mass function
       hmod%halo_DMONLY=3 ! Top-hat halo profile
    ELSE IF(ihm==27) THEN
       hmod%imf=1 ! Press & Schecter (1974) mass function
    ELSE IF(ihm==29) THEN
       ! Adding some cold gas
       hmod%fcold=0.1
    ELSE IF(ihm==30) THEN
       ! Adding some hot gas
       hmod%fhot=0.2
    ELSE IF(ihm==32 .OR. ihm==33 .OR. ihm==34 .OR. ihm==35 .OR. ihm==36 .OR. ihm==37) THEN
       ! 32 - Response model for AGN 7.6;   z = 0.0; simple pivot
       ! 33 - Response model for AGN tuned; z = 0.0; simple pivot
       ! 34 - Response model for AGN 8.0;   z = 0.0; simple pivot
       ! 35 - Response model for AGN 7.6;   z = 0.0; clever pivot; f_hot
       ! 36 - Response model for AGN tuned; z = 0.0; clever pivot; f_hot
       ! 37 - Response model for AGN 8.0;   z = 0.0; clever pivot; f_hot
       hmod%response=1
       IF(ihm==32) THEN
          ! AGN 7.6
          ! Mpiv = 1e14; z = 0.0
          hmod%simple_pivot=.TRUE.
          hmod%alpha=0.709
          hmod%eps=1.12
          hmod%Gamma=1.236
          hmod%M0=10.**13.6
          hmod%Astar=0.045
          hmod%Twhim=10.**6.07
          hmod%cstar=11.2
          hmod%fcold=0.021
          hmod%Mstar=10.**11.5
          hmod%sstar=0.85
          hmod%alphap=-0.479
          hmod%Gammap=-0.020
          hmod%cstarp=-0.459
       ELSE IF(ihm==33) THEN
          ! AGN tuned
          ! Mpiv = 1e14; z = 0.0
          hmod%simple_pivot=.TRUE.
          hmod%alpha=0.802
          hmod%eps=1.06
          hmod%Gamma=1.268
          hmod%M0=10.**13.9
          hmod%Astar=0.043
          hmod%Twhim=10.**6.09
          hmod%cstar=12.3
          hmod%fcold=0.019
          hmod%Mstar=10.**10.6
          hmod%sstar=1.00
          hmod%alphap=-0.518
          hmod%Gammap=-0.033
          hmod%cstarp=-0.523
       ELSE IF(ihm==34) THEN
          ! AGN 8.0
          ! Mpiv = 1e14; z = 0.0
          hmod%simple_pivot=.TRUE.
          hmod%alpha=0.769
          hmod%eps=0.91
          hmod%Gamma=1.274
          hmod%M0=10.**14.3
          hmod%Astar=0.039
          hmod%Twhim=10.**6.11
          hmod%cstar=12.7
          hmod%fcold=0.012
          hmod%Mstar=10.**10.
          hmod%sstar=1.16
          hmod%alphap=-0.409
          hmod%Gammap=-0.029
          hmod%cstarp=-0.510
       ELSE IF(ihm==35) THEN
          ! AGN 7.6
          ! Mpiv = Mh; z = 0.0; f_hot
          hmod%alpha=1.247          
          hmod%eps=1.087
          hmod%Gamma=1.253
          hmod%M0=10.**13.63
          hmod%Astar=0.042
          hmod%Twhim=10.**6.08
          hmod%cstar=10.80
          hmod%fcold=0.0153
          hmod%Mstar=10.**12.12
          hmod%sstar=0.916
          hmod%alphap=-0.488
          hmod%Gammap=-0.0164
          hmod%cstarp=-0.222
          hmod%fhot=0.0323
       ELSE IF(ihm==36) THEN
          ! AGN tuned
          ! Mpiv = Mh; z = 0.0; f_hot
          hmod%alpha=1.251
          hmod%eps=0.866
          hmod%Gamma=1.257
          hmod%M0=10.**13.856
          hmod%Astar=0.0413
          hmod%Twhim=10.**6.09
          hmod%cstar=11.77
          hmod%fcold=0.0066
          hmod%Mstar=10.**12.12
          hmod%sstar=0.856
          hmod%alphap=-0.481
          hmod%Gammap=-0.019
          hmod%cstarp=-0.303
          hmod%fhot=0.0282       
       ELSE IF(ihm==37) THEN
          ! AGN 8.0
          ! Mpiv = Mh; z = 0.0; f_hot
          hmod%alpha=1.024
          hmod%eps=0.821
          hmod%Gamma=1.264
          hmod%M0=10.**14.33
          hmod%Astar=0.0383
          hmod%Twhim=10.**6.10
          hmod%cstar=17.09
          hmod%fcold=0.00267
          hmod%Mstar=10.**11.52
          hmod%sstar=0.99
          hmod%alphap=-0.342
          hmod%Gammap=-0.018
          hmod%cstarp=-0.413
          hmod%fhot=0.0031
       ELSE
          STOP 'ASSIGN_HALOMOD: Error, ihm specified incorrectly'
       END IF
       hmod%beta=hmod%alpha
       hmod%betap=hmod%alphap
    ELSE IF(ihm==38 .OR. ihm==39 .OR. ihm==40) THEN
       ! 38 - AGN 7p6
       ! 39 - AGN tuned
       ! 49 - AGN 8p0
       hmod%response=1
       IF(ihm==38) THEN
          ! AGN 7p6
          hmod%alpha =   1.52016437    
          hmod%eps =   1.06684244    
          hmod%Gamma =   1.24494147    
          hmod%M0 =   2.84777660E+13
          hmod%Astar =   3.79242003E-02
          hmod%Twhim =   987513.562    
          hmod%cstar =   9.78754425    
          hmod%fcold =   1.83899999E-02
          hmod%mstar =   2.68670049E+12
          hmod%sstar =   1.13123488    
          hmod%alphap = -0.531507611    
          hmod%Gammap =  -6.00000005E-03
          hmod%cstarp = -0.113370098    
          hmod%fhot =   1.08200002E-04
          hmod%alphaz =  0.346108794    
          hmod%Gammaz =  0.231210202    
          hmod%M0z =  -9.32227001E-02
          hmod%Astarz = -0.536369383    
          hmod%Twhimz = -0.139042795    
          hmod%eta = -0.222550094
       ELSE IF(ihm==39) THEN
          ! AGN tuned
          hmod%alpha =   1.54074240    
          hmod%eps =   1.01597583    
          hmod%Gamma =   1.24264216    
          hmod%M0 =   6.51173707E+13
          hmod%Astar =   3.45238000E-02
          hmod%Twhim =   1389362.50    
          hmod%cstar =   10.4484358    
          hmod%fcold =   6.32560020E-03
          hmod%mstar =   2.34850982E+12
          hmod%sstar =   1.25916803    
          hmod%alphap = -0.513900995    
          hmod%Gammap =  -5.66239981E-03
          hmod%cstarp =  -6.11630008E-02
          hmod%fhot =   1.26100000E-04
          hmod%alphaz =  0.340969592    
          hmod%Gammaz =  0.295596898    
          hmod%M0z =  -8.13719034E-02
          hmod%Astarz = -0.545276821    
          hmod%Twhimz = -0.122411400    
          hmod%eta = -0.266318709
       ELSE IF(ihm==40) THEN
          ! AGN 8p0
          hmod%alpha =   1.45703220    
          hmod%eps =  0.872408926    
          hmod%Gamma =   1.24960959    
          hmod%M0 =   1.15050950E+14
          hmod%Astar =   3.86818014E-02
          hmod%Twhim =   1619561.00    
          hmod%cstar =   19.4119701    
          hmod%fcold =   1.10999999E-05
          hmod%mstar =   2.18769510E+12
          hmod%sstar =  0.803893507    
          hmod%alphap = -0.528370678    
          hmod%Gammap =  -3.31420009E-03
          hmod%cstarp = -0.355121315    
          hmod%fhot =   3.90500005E-04
          hmod%alphaz =  0.740169585    
          hmod%Gammaz =  0.354409009    
          hmod%M0z =  -2.40819994E-02
          hmod%Astarz = -0.425019890    
          hmod%Twhimz =  -8.60318989E-02
          hmod%eta = -0.243649304
       ELSE
          STOP 'ASSIGN_HALOMOD: Error, ihm specified incorrectly'
       END IF
       hmod%beta=hmod%alpha
       hmod%betap=hmod%alphap
       hmod%betaz=hmod%alphaz
    ELSE IF(ihm==41) THEN
       ! Some stellar mass in satellite galaxies
       hmod%eta=-0.3
    ELSE IF(ihm==42) THEN
       ! Tinker and stuff apprpriate for M200c
       hmod%imf=3   ! Tinker mass function and bias
       hmod%iDv=7   ! M200c
       hmod%iconc=5 ! Duffy for M200c for full sample
       hmod%idc=1   ! Fixed to 1.686
    ELSE IF(ihm==43) THEN
       ! Standard halo model but as response with HMcode but only for matter spectra
       hmod%response=2
    ELSE IF(ihm==44) THEN
       ! Tinker and stuff apprpriate for M200
       hmod%imf=3   ! Tinker mass function and bias
       hmod%iDv=1   ! M200
       hmod%iconc=3 ! Duffy for M200 for full sample
       hmod%idc=1   ! Fixed to 1.686
    ELSE IF(ihm==45) THEN
       ! No stars
       hmod%Astar=0.
    ELSE IF(ihm==46) THEN
       ! Isothermal beta model
       hmod%halo_static_gas=2
    ELSE IF(ihm==47) THEN
       ! Isothermal beta model, response
       hmod%halo_static_gas=2
       hmod%response=1
    ELSE IF(ihm==48) THEN
       ! Non-linear halo bias for standard halo model with virial haloes
       hmod%ibias=3   ! Non-linear halo bias
       !hmod%i1hdamp=3 ! One-halo damping like k^4
       !hmod%ikstar=2  ! One-halo damping via k* from Mead et al. (2015)
    ELSE IF(ihm==49) THEN
       ! Non-linear halo bias and Tinker mass function and virial mass haloes
       hmod%imf=3     ! Tinker mass function and bias
       hmod%ibias=3   ! Non-linear halo bias
       !hmod%i1hdamp=3 ! One-halo damping like k^4
       !hmod%ikstar=2  ! One-halo damping via k* from Mead et al. (2015)
    ELSE
       STOP 'ASSIGN_HALOMOD: Error, ihm specified incorrectly'
    END IF
    hmod%name=names(ihm)

    IF(verbose) THEN
       WRITE(*,*) 'ASSIGN_HALOMOD: ', TRIM(names(ihm))
       WRITE(*,*) 'ASSIGN_HALOMOD: Done'
       WRITE(*,*)
    END IF

  END SUBROUTINE assign_halomod

  SUBROUTINE init_halomod(mmin,mmax,a,hmod,cosm,verbose)

    ! Halo-model initialisation routine
    ! The computes other tables necessary for the one-halo integral
    IMPLICIT NONE
    REAL, INTENT(IN) :: mmin, mmax
    REAL, INTENT(IN) :: a
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    LOGICAL, INTENT(IN) :: verbose
    INTEGER :: i
    REAL :: Dv, dc, m, nu, R, sig, z, Om_stars
    REAL, PARAMETER :: glim=-1e-4
    REAL, PARAMETER :: gblim=-1e-4

    ! Set the redshift (this routine needs to be called anew for each z)
    hmod%a=a
    z=redshift_a(a)
    hmod%z=z

    IF(ALLOCATED(hmod%log_m)) CALL deallocate_HMOD(hmod)
    CALL allocate_HMOD(hmod)

    ! Reset flags to false
    hmod%has_galaxies=.FALSE.
    hmod%has_HI=.FALSE.
    hmod%has_mass_conversions=.FALSE.
    hmod%has_dewiggle=.FALSE.
    hmod%has_mass_function=.FALSE.
    hmod%has_bnl=.FALSE.

    ! Find value of sigma_V
    hmod%sigv=sigmaV(0.,a,cosm)
    IF(hmod%i2hdamp==3) hmod%sigv100=sigmaV(100.,a,cosm)
    hmod%sig8z=sigma(8.,a,cosm)

    IF(verbose) THEN
       WRITE(*,*) 'INIT_HALOMOD: Filling look-up tables'
       WRITE(*,*) 'INIT_HALOMOD: Number of entries:', hmod%n
       WRITE(*,*) 'INIT_HALOMOD: Tables being filled at redshift:', REAL(z)
       WRITE(*,*) 'INIT_HALOMOD: Tables being filled at scale-factor:', REAL(a)
       WRITE(*,*) 'INIT_HALOMOD: sigma_V [Mpc/h]:', REAL(hmod%sigv)
       IF(hmod%i2hdamp==3) WRITE(*,*) 'INIT_HALOMOD: sigmaV_100 [Mpc/h]:', REAL(hmod%sigv100)
       WRITE(*,*) 'INIT_HALOMOD: sigma_8(z):', REAL(hmod%sig8z)
    END IF

    ! Loop over halo masses and fill arrays
    DO i=1,hmod%n

       m=exp(progression(log(mmin),log(mmax),i,hmod%n))
       R=radius_m(m,cosm)
       sig=sigma(R,a,cosm)
       nu=nu_R(R,hmod,cosm)

       hmod%m(i)=m
       hmod%rr(i)=R
       hmod%sig(i)=sig
       hmod%nu(i)=nu

    END DO

    ! log-mass array
    hmod%log_m=log(hmod%m)

    IF(verbose) WRITE(*,*) 'INIT_HALOMOD: M, R, nu, sigma tables filled'

    ! Get delta_c
    dc=delta_c(hmod,cosm)
    hmod%dc=dc

    ! Fill virial radius table using real radius table
    Dv=Delta_v(hmod,cosm)
    hmod%Dv=Dv
    hmod%rv=hmod%rr/(Dv**(1./3.))

    ! Write some useful information to the screen
    IF(verbose) THEN
       WRITE(*,*) 'INIT_HALOMOD: virial radius tables filled'
       WRITE(*,*) 'INIT_HALOMOD: Delta_v:', REAL(Dv)
       WRITE(*,*) 'INIT_HALOMOD: delta_c:', REAL(dc)
       WRITE(*,*) 'INIT_HALOMOD: Minimum nu:', REAL(hmod%nu(1))
       WRITE(*,*) 'INIT_HALOMOD: Maximum nu:', REAL(hmod%nu(hmod%n))
       WRITE(*,*) 'INIT_HALOMOD: Minimum R_v [Mpc/h]:', REAL(hmod%rv(1))
       WRITE(*,*) 'INIT_HALOMOD: Maximum R_v [Mpc/h]:', REAL(hmod%rv(hmod%n))
       WRITE(*,*) 'INIT_HALOMOD: Minimum log10(M) [Msun/h]:', REAL(log10(hmod%m(1)))
       WRITE(*,*) 'INIT_HALOMOD: Maximum log10(M) [Msun/h]:', REAL(log10(hmod%m(hmod%n)))
    END IF

    IF(hmod%imf==4) THEN

       ! Do nothing

    ELSE

       ! Calculate missing mass things if necessary
       IF(hmod%ip2h_corr==2 .OR. hmod%ip2h_corr==3) THEN

          hmod%gmin=1.-integrate_g_nu(hmod%nu(1),hmod%large_nu,hmod)
          hmod%gmax=integrate_g_nu(hmod%nu(hmod%n),hmod%large_nu,hmod)
          hmod%gbmin=1.-integrate_gb_nu(hmod%nu(1),hmod%large_nu,hmod)
          hmod%gbmax=integrate_gb_nu(hmod%nu(hmod%n),hmod%large_nu,hmod)          
          !hmod%gnorm=integrate_g_nu(hmod%small_nu,hmod%large_nu,hmod)
          hmod%gnorm=integrate_g_mu(hmod%small_nu,hmod%large_nu,hmod)

          IF(verbose) THEN
             WRITE(*,*) 'INIT_HALOMOD: Missing g(nu) at low end:', REAL(hmod%gmin)
             WRITE(*,*) 'INIT_HALOMOD: Missing g(nu) at high end:', REAL(hmod%gmax)
             WRITE(*,*) 'INIT_HALOMOD: Missing g(nu)b(nu) at low end:', REAL(hmod%gbmin)
             WRITE(*,*) 'INIT_HALOMOD: Missing g(nu)b(nu) at high end:', REAL(hmod%gbmax)
             WRITE(*,*) 'INIT_HALOMOD: Total g(nu) integration:', REAL(hmod%gnorm) ! Could do better and actually integrate to zero
          END IF

          IF(hmod%gmin<0.)     STOP 'INIT_HALOMOD: Error, missing g(nu) at low end is less than zero'       
          IF(hmod%gmax<glim)   STOP 'INIT_HALOMOD: Error, missing g(nu) at high end is less than zero'
          IF(hmod%gbmin<0.)    STOP 'INIT_HALOMOD: Error, missing g(nu)b(nu) at low end is less than zero'
          IF(hmod%gbmax<gblim) STOP 'INIT_HALOMOD: Error, missing g(nu)b(nu) at high end is less than zero'   

       END IF

    END IF

    ! Find non-linear radius and scale
    ! This is defined as nu(M_star)=1 *not* sigma(M_star)=1, so depends on delta_c
    hmod%rnl=r_nl(hmod)
    hmod%mnl=mass_r(hmod%rnl,cosm)
    hmod%knl=1./hmod%rnl

    IF(verbose) THEN
       WRITE(*,*) 'INIT_HALOMOD: Non-linear mass [log10(M*) [Msun/h]]:', REAL(log10(hmod%mnl))
       WRITE(*,*) 'INIT_HALOMOD: Non-linear halo virial radius [Mpc/h]:', REAL(virial_radius(hmod%mnl,hmod,cosm))
       WRITE(*,*) 'INIT_HALOMOD: Non-linear Lagrangian radius [Mpc/h]:', REAL(hmod%rnl)
       WRITE(*,*) 'INIT_HALOMOD: Non-linear wavenumber [h/Mpc]:', REAL(hmod%knl)
    END IF
    
    !WRITE(*,*) 'INIT_HALOMOD: Cumulative halo number density above M* [(Mpc/h)^-3]:', REAL(cumulative_halo_density(hmod%mnl,hmod,cosm))
    !STOP

    hmod%neff=effective_index(hmod,cosm)

    IF(verbose) WRITE(*,*) 'INIT_HALOMOD: Collapse n_eff:', REAL(hmod%neff)

    CALL fill_halo_concentration(hmod,cosm)

    IF(verbose) THEN
       WRITE(*,*) 'INIT_HALOMOD: Halo concentration tables filled'
       WRITE(*,*) 'INIT_HALOMOD: Minimum concentration:', REAL(hmod%c(hmod%n))
       WRITE(*,*) 'INIT_HALOMOD: Maximum concentration:', REAL(hmod%c(1))
    END IF

    !IF(slow_hmod) THEN
    hmod%Rh=one_halo_amplitude(hmod,cosm)
    hmod%Mh=hmod%Rh*comoving_matter_density(cosm)
    IF(verbose) THEN
       WRITE(*,*) 'INIT_HALOMOD: One-halo amplitude [Mpc/h]^3:', hmod%Rh
       WRITE(*,*) 'INIT_HALOMOD: One-halo amplitude [log10(M) [Msun/h]]:', log10(hmod%Mh)
    END IF
    !END IF

    ! Calculate the total stellar mass fraction
    IF(verbose) THEN
       Om_stars=Omega_stars(hmod,cosm)
       WRITE(*,*) 'INIT_HALOMOD: Omega_*:', Om_stars
       WRITE(*,*) 'INIT_HALOMOD: Omega_* / Omega_m:', Om_stars/cosm%Om_m
    END IF

    IF(verbose) THEN
       WRITE(*,*) 'INIT_HALOMOD: Done'
       WRITE(*,*)
    END IF

  END SUBROUTINE init_halomod

  REAL FUNCTION integrate_g_mu(nu1,nu2,hmod)

    ! Integrate g(nu) between nu1 and nu2; this is the fraction of mass in the Universe in haloes between nu1 and nu2
    ! Integrating this over all nu should give unity
    ! Could do better but need to explicity remove divergence by making functions where alpha*mu^(alpha-1) is multiplied by g(nu)
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu1, nu2         ! Range in nu
    TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model
    INTEGER, PARAMETER :: iorder=3       ! Order for integration
    REAL :: mu1, mu2
    REAL :: alpha

    ! Relation between nu and mu: nu=mu^alpha
    alpha=hmod%alpha_numu
    mu1=nu1**(1./alpha)
    mu2=nu2**(1./alpha)

    integrate_g_mu=integrate_hmod(mu1,mu2,g_mu,hmod,hmod%acc_HMx,iorder)

  END FUNCTION integrate_g_mu

  REAL FUNCTION integrate_g_nu(nu1,nu2,hmod)

    ! Integrate g(nu) between nu1 and nu2; this is the fraction of mass in the Universe in haloes between nu1 and nu2
    ! Previously called 'mass_interval'
    ! Integrating this over all nu should give unity
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu1, nu2         ! Range in nu
    TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model
    INTEGER, PARAMETER :: iorder=3       ! Order for integration

    integrate_g_nu=integrate_hmod(nu1,nu2,g_nu,hmod,hmod%acc_HMx,iorder)

  END FUNCTION integrate_g_nu

  REAL FUNCTION integrate_g_nu_on_M(nu1,nu2,hmod)

    ! Integrate g(nu)/M between nu1 and nu2; this is the number density of haloes in the Universe between nu1 and nu2
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu1, nu2         ! Range in nu
    TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model
    INTEGER, PARAMETER :: iorder=3       ! Order for integration

    integrate_g_nu_on_M=integrate_hmod(nu1,nu2,g_nu_on_M,hmod,hmod%acc_HMx,iorder)

  END FUNCTION integrate_g_nu_on_M

  REAL FUNCTION integrate_gb_nu(nu1,nu2,hmod)

    ! Integrate b(nu)*g(nu) between nu1 and nu2
    ! Previously called 'bias_interval'
    ! This is the unnormalised mass-weighted mean bias in the range nu1 to nu2
    ! Integrating this over all nu should give unity
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu1, nu2         ! Range in nu
    TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model
    INTEGER, PARAMETER :: iorder=3       ! Order for integration

    integrate_gb_nu=integrate_hmod(nu1,nu2,gb_nu,hmod,hmod%acc_HMx,iorder)

  END FUNCTION integrate_gb_nu

  REAL FUNCTION integrate_gb_nu_on_M(nu1,nu2,hmod)

    ! Integrate b(nu)*g(nu) between nu1 and nu2
    ! Previously called 'bias_interval'
    ! This is the unnormalised mass-weighted mean bias in the range nu1 to nu2
    ! Integrating this over all nu should give unity
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu1, nu2         ! Range in nu
    TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model
    INTEGER, PARAMETER :: iorder=3       ! Order for integration

    integrate_gb_nu_on_M=integrate_hmod(nu1,nu2,gb_nu_on_M,hmod,hmod%acc_HMx,iorder)

  END FUNCTION integrate_gb_nu_on_M

  REAL FUNCTION integrate_nug_nu(nu1,nu2,hmod)

    ! Integrate nu*g(nu) between nu1 and nu2
    ! Previously called 'nu interval'
    ! This is the unnormalised mass-weighted nu in the range nu1 to nu2
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu1, nu2         ! Range in nu
    TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model
    INTEGER, PARAMETER :: iorder=3       ! Order for integration

    integrate_nug_nu=integrate_hmod(nu1,nu2,nug_nu,hmod,hmod%acc_HMx,iorder)

  END FUNCTION integrate_nug_nu

  REAL FUNCTION integrate_Mg_nu(nu1,nu2,hmod)

    ! Integrate nu*g(nu) between nu1 and nu2
    ! Previously called 'nu interval'
    ! This is the unnormalised mass-weighted nu in the range nu1 to nu2
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu1, nu2         ! Range in nu
    TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model
    INTEGER, PARAMETER :: iorder=3       ! Order for integration

    integrate_Mg_nu=integrate_hmod(nu1,nu2,Mg_nu,hmod,hmod%acc_HMx,iorder)

  END FUNCTION integrate_Mg_nu

  REAL FUNCTION integrate_nug_nu_on_M(nu1,nu2,hmod)

    ! Integrate nu*g(nu) between nu1 and nu2
    ! Previously called 'nu interval'
    ! This is the unnormalised mass-weighted nu in the range nu1 to nu2
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu1, nu2         ! Range in nu
    TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model
    INTEGER, PARAMETER :: iorder=3       ! Order for integration

    integrate_nug_nu_on_M=integrate_hmod(nu1,nu2,nug_nu_on_M,hmod,hmod%acc_HMx,iorder)

  END FUNCTION integrate_nug_nu_on_M

  REAL FUNCTION mean_bias_number_weighted(nu1,nu2,hmod)

    ! Calculate the mean number-weighted bias in the range nu1 to nu2
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu1, nu2         ! Range in nu
    TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model
    INTEGER, PARAMETER :: iorder=3       ! Order for integration

    mean_bias_number_weighted=integrate_gb_nu_on_M(nu1,nu2,hmod)/integrate_g_nu_on_M(nu1,nu2,hmod)

  END FUNCTION mean_bias_number_weighted

  REAL FUNCTION mean_nu_number_weighted(nu1,nu2,hmod)

    ! Calculate the mean number-weighted nu in the range nu1 to nu2
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu1, nu2         ! Range in nu
    TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model
    INTEGER, PARAMETER :: iorder=3       ! Order for integration

    mean_nu_number_weighted=integrate_nug_nu_on_M(nu1,nu2,hmod)/integrate_g_nu_on_M(nu1,nu2,hmod)

  END FUNCTION mean_nu_number_weighted

  REAL FUNCTION mean_halo_mass_number_weighted(nu1,nu2,hmod)

    ! Calculate the mean number-weighted halo mass in the range nu1 to nu2
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu1, nu2         ! Range in nu
    TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model
    INTEGER, PARAMETER :: iorder=3       ! Order for integration

    mean_halo_mass_number_weighted=integrate_g_nu(nu1,nu2,hmod)/integrate_g_nu_on_M(nu1,nu2,hmod)

  END FUNCTION mean_halo_mass_number_weighted

  REAL FUNCTION mean_bias_mass_weighted(nu1,nu2,hmod)

    ! Calculate the mean mass-weighted bias in the range nu1 to nu2
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu1, nu2         ! Range in nu
    TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model
    INTEGER, PARAMETER :: iorder=3       ! Order for integration

    mean_bias_mass_weighted=integrate_gb_nu(nu1,nu2,hmod)/integrate_g_nu(nu1,nu2,hmod)

  END FUNCTION mean_bias_mass_weighted

  REAL FUNCTION mean_nu_mass_weighted(nu1,nu2,hmod)

    ! Calculate the mean mass-weighted nu in the range nu1 to nu2
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu1, nu2         ! Range in nu
    TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model
    INTEGER, PARAMETER :: iorder=3       ! Order for integration

    mean_nu_mass_weighted=integrate_nug_nu(nu1,nu2,hmod)/integrate_g_nu(nu1,nu2,hmod)

  END FUNCTION mean_nu_mass_weighted

  REAL FUNCTION mean_halo_mass_mass_weighted(nu1,nu2,hmod)

    ! Calculate the mean mass-weighted halo mass in the range nu1 to nu2
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu1, nu2         ! Range in nu
    TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model
    INTEGER, PARAMETER :: iorder=3       ! Order for integration

    mean_halo_mass_mass_weighted=integrate_Mg_nu(nu1,nu2,hmod)/integrate_g_nu(nu1,nu2,hmod)

  END FUNCTION mean_halo_mass_mass_weighted

  REAL FUNCTION mean_halo_number_density(nu1,nu2,hmod,cosm)

    ! Calculate N(m) where N is the number density of haloes above mass m
    ! Obtained by integrating the mass function
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu1, nu2           ! Range in nu
    TYPE(halomod), INTENT(INOUT) :: hmod   ! Halo model
    TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology
    INTEGER, PARAMETER :: iorder=3         ! Order for integration

    mean_halo_number_density=integrate_hmod(nu1,nu2,g_nu_on_M,hmod,hmod%acc_HMx,iorder)
    mean_halo_number_density=mean_halo_number_density*comoving_matter_density(cosm)
    
  END FUNCTION mean_halo_number_density

!!$  REAL FUNCTION halo_number_density(m1,m2,hmod,cosm)
!!$
!!$    IMPLICIT NONE
!!$    REAL, INTENT(IN) :: m1
!!$    REAL, INTENT(IN) :: m2
!!$    TYPE(halomod), INTENT(INOUT) :: hmod
!!$    TYPE(cosmology), INTENT(INOUT) :: cosm
!!$    REAL :: nu1, nu2
!!$
!!$    nu1=nu_M(m1,hmod,cosm)
!!$    nu2=nu_M(m2,hmod,cosm)
!!$
!!$    halo_number_density=comoving_matter_density(cosm)*integrate_g_nu_on_M(nu1,nu2,hmod)
!!$    
!!$  END FUNCTION halo_number_density

  SUBROUTINE print_halomod(hmod,cosm,verbose)

    !This subroutine writes out the physical halo-model parameters at some redshift 
    !(e.g., Delta_v) rather than the model parameters
    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    LOGICAL, INTENT(IN) :: verbose

    IF(verbose) THEN

       WRITE(*,*) '============================================'
       WRITE(*,*) 'HALOMODEL: ', TRIM(hmod%name)
       WRITE(*,*) '============================================'
       WRITE(*,*) 'HALOMODEL: Accuracy parameters'
       WRITE(*,*) '============================================'
       WRITE(*,*) 'Number of points in look-up tables:', hmod%n
       WRITE(*,*) 'Halo-model accuracy parameter:', hmod%acc_HMx
       WRITE(*,*) 'Large value of nu:', hmod%large_nu
       WRITE(*,*) '============================================'

       ! Form of the two-halo term
       IF(hmod%ip2h==1) WRITE(*,*) 'HALOMODEL: Linear two-halo term'
       IF(hmod%ip2h==2) WRITE(*,*) 'HALOMODEL: Standard two-halo term (Seljak 2000)'
       IF(hmod%ip2h==3) WRITE(*,*) 'HALOMODEL: Linear two-halo term with damped wiggles'
       IF(hmod%ip2h==4) WRITE(*,*) 'HALOMODEL: No two-halo term'

       ! Order to go to in halo bias
       IF(hmod%ip2h==2) THEN
          IF(hmod%ibias==1) WRITE(*,*) 'HALOMODEL: Linear halo bias'
          IF(hmod%ibias==2) WRITE(*,*) 'HALOMODEL: Second-order halo bias'
          IF(hmod%ibias==3) WRITE(*,*) 'HALOMODEL: Full non-linear halo bias'
          IF(hmod%ibias==4) WRITE(*,*) 'HALOMODEL: Scale-dependent halo bias fudge from Fedeli (2014b)'
          IF(hmod%ibias==5) WRITE(*,*) 'HALOMODEL: Scale-dependent halo bias fudge from me'
       END IF

       ! Correction for missing low-mass haloes
       IF(hmod%ip2h==2) THEN
          IF(hmod%ip2h_corr==1) WRITE(*,*) 'HALOMODEL: No two-halo correction applied for missing low-mass haloes'
          IF(hmod%ip2h_corr==2) WRITE(*,*) 'HALOMODEL: Two-halo term corrected by adding missing g(nu)b(nu)' 
          IF(hmod%ip2h_corr==3) WRITE(*,*) 'HALOMODEL: Two-halo term corrected via delta function at low mass end'      
       END IF

       ! Halo mass function
       IF(hmod%imf==1) WRITE(*,*) 'HALOMODEL: Press & Schecter (1974) mass function'
       IF(hmod%imf==2) WRITE(*,*) 'HALOMODEL: Sheth & Tormen (1999) mass function'
       IF(hmod%imf==3) WRITE(*,*) 'HALOMODEL: Tinker et al. (2010) mass function'
       IF(hmod%imf==4) WRITE(*,*) 'HALOMODEL: Delta function mass function'

       ! Concentration-mass relation
       IF(hmod%iconc==1) WRITE(*,*) 'HALOMODEL: Full Bullock et al. (2001) concentration-mass relation'
       IF(hmod%iconc==2) WRITE(*,*) 'HALOMODEL: Simple Bullock et al. (2001) concentration-mass relation'
       IF(hmod%iconc==3) WRITE(*,*) 'HALOMODEL: Full sample for M200 Duffy et al. (2008) concentration-mass relation'
       IF(hmod%iconc==4) WRITE(*,*) 'HALOMODEL: Full sample for Mv Duffy et al. (2008) concentration-mass relation'
       IF(hmod%iconc==5) WRITE(*,*) 'HALOMODEL: Full sample for M200c Duffy et al. (2008) concentration-mass relation'
       IF(hmod%iconc==6) WRITE(*,*) 'HALOMODEL: Relaxed sample for M200 Duffy et al. (2008) concentration-mass relation'
       IF(hmod%iconc==7) WRITE(*,*) 'HALOMODEL: Relaxed sample for Mv Duffy et al. (2008) concentration-mass relation'
       IF(hmod%iconc==8) WRITE(*,*) 'HALOMODEL: Relaxed sample for M200c Duffy et al. (2008) concentration-mass relation'

       ! Concentration-mass relation correction
       IF(hmod%iDolag==1) WRITE(*,*) 'HALOMODEL: No concentration-mass correction for dark energy'
       IF(hmod%iDolag==2) WRITE(*,*) 'HALOMODEL: Dolag (2004) dark energy halo concentration correction'
       IF(hmod%iDolag==3) WRITE(*,*) 'HALOMODEL: Dolag (2004) dark energy halo concentration correction with 1.5 exponent'
       IF(hmod%iDolag==4) WRITE(*,*) 'HALOMODEL: Dolag (2004) dark energy halo concentration correction with redshift dependence'

       ! Bound gas fraction
       IF(hmod%frac_bound_gas==1) WRITE(*,*) 'HALOMODEL: Halo bound gas fraction: Fedeli (2014)'
       IF(hmod%frac_bound_gas==2) WRITE(*,*) 'HALOMODEL: Halo bound gas fraction: Schneider & Teyssier (2015)'
       IF(hmod%frac_bound_gas==3) WRITE(*,*) 'HALOMODEL: Halo bound gas fraction: Universal baryon fraction'

       ! Cold gas fraction
       IF(hmod%frac_cold_bound_gas==1) WRITE(*,*) 'HALOMODEL: Cold gas fraction: Constant fraction of bound halo gas'

       ! Hot gas fraction
       IF(hmod%frac_hot_bound_gas==1) WRITE(*,*) 'HALOMODEL: Hot gas fraction: Constant fraction of bound halo gas'

       ! Central star fraction
       IF(hmod%frac_stars==1) WRITE(*,*) 'HALOMODEL: Halo star fraction: Fedeli (2014)'
       IF(hmod%frac_stars==2) WRITE(*,*) 'HALOMODEL: Halo star fraction: Constant'
       IF(hmod%frac_stars==3) WRITE(*,*) 'HALOMODEL: Halo star fraction: Fedeli (2014) but saturated at high halo mass'
       IF(hmod%frac_stars==4) WRITE(*,*) 'HALOMODEL: Halo star fraction: No stars in haloes'

       ! Satellite star fraction
       IF(hmod%frac_central_stars==1) WRITE(*,*) 'HALOMODEL: Central star fraction: All central stars'
       IF(hmod%frac_central_stars==2) WRITE(*,*) 'HALOMODEL: Central star fraction: Schneider et al. (2018)'

       ! HI fraction
       IF(hmod%frac_HI==1) WRITE(*,*) 'HALOMODEL: Halo HI fraction: Hard truncation at high and low masses'
       IF(hmod%frac_HI==2) WRITE(*,*) 'HALOMODEL: Halo HI fraction: Villaescusa-Navarro et al. (2018; 1804.09180)'

       ! DMONLY halo model
       IF(hmod%halo_DMONLY==1) WRITE(*,*) 'HALOMODEL: DMONLY halo profile: NFW'
       IF(hmod%halo_DMONLY==2) WRITE(*,*) 'HALOMODEL: DMONLY halo profile: Non-analytical NFW (good for testing W(k) functions)'
       IF(hmod%halo_DMONLY==3) WRITE(*,*) 'HALOMODEL: DMONLY halo profile: Tophat'
       IF(hmod%halo_DMONLY==4) WRITE(*,*) 'HALOMODEL: DMONLY halo profile: Delta function'
       IF(hmod%halo_DMONLY==5) WRITE(*,*) 'HALOMODEL: DMONLY halo profile: Cored NFW'
       IF(hmod%halo_DMONLY==6) WRITE(*,*) 'HALOMODEL: DMONLY halo profile: Isothermal'
       IF(hmod%halo_DMONLY==7) WRITE(*,*) 'HALOMODEL: DMONLY halo profile: Shell'

       ! CDM halo profile
       IF(hmod%halo_CDM==1) WRITE(*,*) 'HALOMODEL: CDM halo profile: NFW'

       ! Bound gas halo profile
       IF(hmod%halo_static_gas==1) WRITE(*,*) 'HALOMODEL: Static gas profile: Simplified Komatsu & Seljak (2001)'
       IF(hmod%halo_static_gas==2) WRITE(*,*) 'HALOMODEL: Static gas profile: Isothermal beta profile'
       IF(hmod%halo_static_gas==3) WRITE(*,*) 'HALOMODEL: Static gas profile: Full Komatsu & Seljak (2001)'

       ! Cold gas profile
       IF(hmod%halo_cold_gas==1) WRITE(*,*) 'HALOMODEL: Cold gas profile: Delta function'

       ! Hot gas profile
       IF(hmod%halo_hot_gas==1) WRITE(*,*) 'HALOMODEL: Hot gas profile: Isothermal'

       ! Free gas halo profile
       IF(hmod%halo_free_gas==1) WRITE(*,*) 'HALOMODEL: Free gas profile: Isothermal model (out to 2rv)'
       IF(hmod%halo_free_gas==2) WRITE(*,*) 'HALOMODEL: Free gas profile: Ejected gas model from Schneider & Teyssier (2015)'
       IF(hmod%halo_free_gas==3) WRITE(*,*) 'HALOMODEL: Free gas profile: Isothermal shell that connects electron pressure and density to bound gas at rv'
       IF(hmod%halo_free_gas==4) WRITE(*,*) 'HALOMODEL: Free gas profile: Komatsu-Seljak continuation'
       IF(hmod%halo_free_gas==5) WRITE(*,*) 'HALOMODEL: Free gas profile: Power-law continuation'
       IF(hmod%halo_free_gas==6) WRITE(*,*) 'HALOMODEL: Free gas profile: Cubic profile'
       IF(hmod%halo_free_gas==7) WRITE(*,*) 'HALOMODEL: Free gas profile: Smoothly distributed'
       IF(hmod%halo_free_gas==8) WRITE(*,*) 'HALOMODEL: Free gas profile: Delta function'

       ! Central star halo profile
       IF(hmod%halo_central_stars==1) WRITE(*,*) 'HALOMODEL: Central star profile: Fedeli (2014)'
       IF(hmod%halo_central_stars==2) WRITE(*,*) 'HALOMODEL: Central star profile: Schneider & Teyssier (2015)'
       IF(hmod%halo_central_stars==3) WRITE(*,*) 'HALOMODEL: Central star profile: Delta function'
       IF(hmod%halo_central_stars==4) WRITE(*,*) 'HALOMODEL: Central star profile: Delta function at low mass and NFW at high mass'

       ! Satellite star halo profile
       IF(hmod%halo_satellite_stars==1) WRITE(*,*) 'HALOMODEL: Satellite star profile: NFW'
       IF(hmod%halo_satellite_stars==2) WRITE(*,*) 'HALOMODEL: Satellite star profile: Schneider & Teyssier (2015)'

       ! HI halo profile
       IF(hmod%halo_HI==1) WRITE(*,*) 'HALOMODEL: HI profile: NFW'
       IF(hmod%halo_HI==2) WRITE(*,*) 'HALOMODEL: HI profile: Delta function'
       IF(hmod%halo_HI==3) WRITE(*,*) 'HALOMODEL: HI profile: Polynomial with exponential hole (Villaescusa-Navarro et al. 2018)'
       IF(hmod%halo_HI==4) WRITE(*,*) 'HALOMODEL: HI profile: Modified NFW with exponential hole (Villaescusa-Navarro et al. 2018)'
       IF(hmod%halo_HI==5) WRITE(*,*) 'HALOMODEL: HI profile: Modified NFW (Padmanabhan & Refregier 2017)'

       ! Electron pressure profile
       IF(hmod%electron_pressure==1) WRITE(*,*) 'HALOMODEL: Electron pressure: Using UPP'
       IF(hmod%electron_pressure==2) WRITE(*,*) 'HALOMODEL: Electron pressure: Using gas profiles'

       ! Are voids being done?
       IF(hmod%voids) THEN

          WRITE(*,*) 'HALOMODEL: Considering voids'

          ! Void model
          IF(hmod%halo_void==1) WRITE(*,*) 'HALOMODEL: Void model: Top-hat'

          ! Compensated void model
          IF(hmod%halo_compensated_void==1) WRITE(*,*) 'HALOMODEL: Compensated void model: Top-hat'

       END IF

       ! delta_c
       IF(hmod%idc==1) WRITE(*,*) 'HALOMODEL: Fixed delta_c = 1.686'
       IF(hmod%idc==2) WRITE(*,*) 'HALOMODEL: delta_c from Nakamura & Suto (1997) fitting function'
       IF(hmod%idc==3) WRITE(*,*) 'HALOMODEL: delta_c from Mead et al. (2016) power spectrum fit'
       IF(hmod%idc==4) WRITE(*,*) 'HALOMODEL: delta_c from Mead (2017) fitting function'
       IF(hmod%idc==5) WRITE(*,*) 'HALOMODEL: delta_c from spherical-collapse calculation'
       IF(hmod%idc==6) WRITE(*,*) 'HALOMODEL: delta_c from Mead et al. (2015) power spectrum fit'

       ! Delta_v
       IF(hmod%iDv==1) WRITE(*,*) 'HALOMODEL: Fixed Delta_v = 200'
       IF(hmod%iDv==2) WRITE(*,*) 'HALOMODEL: Delta_v from Bryan & Norman (1998) fitting function'
       IF(hmod%iDv==3) WRITE(*,*) 'HALOMODEL: Delta_v from Mead et al. (2016) power spectrum fit'
       IF(hmod%iDv==4) WRITE(*,*) 'HALOMODEL: Delta_v from Mead (2017) fitting function'
       IF(hmod%iDv==5) WRITE(*,*) 'HALOMODEL: Delta_v from spherical-collapse calculation'
       IF(hmod%iDv==6) WRITE(*,*) 'HALOMODEL: Delta_v to give haloes Lagrangian radius'
       IF(hmod%iDv==7) WRITE(*,*) 'HALOMODEL: Fixed Delta_v for M200c'
       IF(hmod%iDv==8) WRITE(*,*) 'HALOMODEL: Delta_v from Mead et al. (2015) power spectrum fit'

       ! eta for halo window function
       IF(hmod%ieta==1) WRITE(*,*) 'HALOMODEL: eta = 0 fixed'
       IF(hmod%ieta==2) WRITE(*,*) 'HALOMODEL: eta from Mead et al. (2015,2016) power spectrum fit'

       ! Small-scale two-halo term damping coefficient
       IF(hmod%i2hdamp==1) WRITE(*,*) 'HALOMODEL: No two-halo term damping at small scales'
       IF(hmod%i2hdamp==2) WRITE(*,*) 'HALOMODEL: Two-halo term damping from Mead et al. (2015)'
       IF(hmod%i2hdamp==3) WRITE(*,*) 'HALOMODEL: Two-halo term damping from Mead et al. (2016)'

       ! Large-scale one-halo term damping function
       IF(hmod%i1hdamp==1) WRITE(*,*) 'HALOMODEL: No damping in one-halo term at large scales'
       IF(hmod%i1hdamp==2) WRITE(*,*) 'HALOMODEL: One-halo term large-scale damping via an exponential'
       IF(hmod%i1hdamp==3) WRITE(*,*) 'HALOMODEL: One-halo term large-scale damping like Delta^2 ~ k^7'

       ! Large-scale one-halo term damping coefficient
       IF(hmod%i1hdamp .NE. 1) THEN
          IF(hmod%ikstar==1) WRITE(*,*) 'HALOMODEL: No damping in one-halo term at large scales'
          IF(hmod%ikstar==2) WRITE(*,*) 'HALOMODEL: One-halo term damping function from Mead et al. (2015, 2016)'
       END IF

       ! Concentration-mass scaling
       IF(hmod%iAs==1) WRITE(*,*) 'HALOMODEL: No rescaling of concentration-mass relation'
       IF(hmod%iAs==2) WRITE(*,*) 'HALOMODEL: Concentration-mass relation rescaled mass independetly (Mead et al. 2015, 2016)'

       ! Two- to one-halo transition region
       IF(hmod%itrans==1) WRITE(*,*) 'HALOMODEL: Standard sum of two- and one-halo terms'
       IF(hmod%itrans==2) WRITE(*,*) 'HALOMODEL: Smoothed transition using alpha'
       IF(hmod%itrans==3) WRITE(*,*) 'HALOMODEL: Experimental smoothed transition'
       IF(hmod%itrans==4) WRITE(*,*) 'HALOMODEL: Experimental smoothed transition for HMx (2.5 exponent)'
       IF(hmod%itrans==5) WRITE(*,*) 'HALOMODEL: Tanh transition with k_nl'

       ! Response
       IF(hmod%response==1) WRITE(*,*) 'HALOMODEL: Power computed as response with HMcode multiplication'
       IF(hmod%response==2) WRITE(*,*) 'HALOMODEL: Matter power computed as response with HMcode multiplication'

       ! Numerical parameters
       WRITE(*,*) '======================================='
       WRITE(*,*) 'HALOMODEL: Standard parameters'
       WRITE(*,*) '======================================='
       WRITE(*,fmt='(A30,F10.5)') 'redshift:', hmod%z
       WRITE(*,fmt='(A30,F10.5)') 'scale factor:', hmod%a
       WRITE(*,fmt='(A30,F10.5)') 'Dv:', Delta_v(hmod,cosm)
       WRITE(*,fmt='(A30,F10.5)') 'dc:', delta_c(hmod,cosm)
       IF(hmod%imf==2 .OR. hmod%imf==3 .OR. hmod%imf==4) THEN
          WRITE(*,*) '======================================='
          WRITE(*,*) 'HALOMODEL: Mass function parameters'
          WRITE(*,*) '======================================='
          IF(hmod%imf==2) THEN
             WRITE(*,fmt='(A30,F10.5)') 'Sheth & Tormen p:', hmod%ST_p
             WRITE(*,fmt='(A30,F10.5)') 'Sheth & Tormen q:', hmod%ST_q
             WRITE(*,fmt='(A30,F10.5)') 'Sheth & Tormen A:', hmod%ST_A            
          ELSE IF(hmod%imf==3) THEN
             WRITE(*,fmt='(A30,F10.5)') 'Tinker alpha:', hmod%Tinker_alpha
             WRITE(*,fmt='(A30,F10.5)') 'Tinker beta:', hmod%Tinker_beta
             WRITE(*,fmt='(A30,F10.5)') 'Tinker gamma:', hmod%Tinker_gamma
             WRITE(*,fmt='(A30,F10.5)') 'Tinker eta:', hmod%Tinker_eta
             WRITE(*,fmt='(A30,F10.5)') 'Tinker phi:', hmod%Tinker_phi
          ELSE IF(hmod%imf==4) THEN
             WRITE(*,*) 'Halo mass log10(M) [Msun/h]:', log10(hmod%hmass)
          ELSE
             STOP 'HALOMODEL: Error, something went wrong with imf'
          END IF
       END IF
       WRITE(*,*) '======================================='
       WRITE(*,*) 'HALOMODEL: HMcode parameters'
       WRITE(*,*) '======================================='
       WRITE(*,fmt='(A30,F10.5)') 'Dv0:', hmod%Dv0
       WRITE(*,fmt='(A30,F10.5)') 'Dv1:', hmod%Dv1
       WRITE(*,fmt='(A30,F10.5)') 'dc0:', hmod%dc0
       WRITE(*,fmt='(A30,F10.5)') 'dc1:', hmod%dc1 
       WRITE(*,fmt='(A30,F10.5)') 'eta0:', hmod%eta0
       WRITE(*,fmt='(A30,F10.5)') 'eta1:', hmod%eta1
       WRITE(*,fmt='(A30,F10.5)') 'f0:', hmod%f0
       WRITE(*,fmt='(A30,F10.5)') 'f1:', hmod%f1
       WRITE(*,fmt='(A30,F10.5)') 'ks:', hmod%ks
       WRITE(*,fmt='(A30,F10.5)') 'As:', hmod%As 
       WRITE(*,fmt='(A30,F10.5)') 'alpha0:', hmod%alp0
       WRITE(*,fmt='(A30,F10.5)') 'alpha1:', hmod%alp1
       WRITE(*,fmt='(A30,F10.5)') 'Dvnu:', hmod%Dvnu
       WRITE(*,fmt='(A30,F10.5)') 'dcnu:', hmod%dcnu
       WRITE(*,*) '======================================='
       WRITE(*,*) 'HALOMODEL: HMcode variables'
       WRITE(*,*) '======================================='
       WRITE(*,fmt='(A30,F10.5)') 'eta:', eta_HMcode(hmod,cosm)
       WRITE(*,fmt='(A30,F10.5)') 'k*:', kstar_HMcode(hmod,cosm)
       WRITE(*,fmt='(A30,F10.5)') 'A:', A_HMcode(hmod,cosm)
       WRITE(*,fmt='(A30,F10.5)') 'fdamp:', fdamp_HMcode(hmod,cosm)
       WRITE(*,fmt='(A30,F10.5)') 'alpha:', alpha_HMcode(hmod,cosm)
       WRITE(*,*) '======================================='
       WRITE(*,*) 'HALOMODEL: HMx parameters'
       WRITE(*,*) '======================================='
       IF(hmod%HMx_mode==1 .OR. hmod%HMx_mode==2 .OR. hmod%HMx_mode==3) THEN
          WRITE(*,fmt='(A30,F10.5)') 'alpha:', hmod%alpha
          WRITE(*,fmt='(A30,F10.5)') 'beta:', hmod%beta
          WRITE(*,fmt='(A30,F10.5)') 'epsilon:', hmod%eps
          WRITE(*,fmt='(A30,F10.5)') 'Gamma:', hmod%Gamma
          WRITE(*,fmt='(A30,F10.5)') 'log10(M0) [Msun/h]:', log10(hmod%M0)          
          WRITE(*,fmt='(A30,F10.5)') 'A*:', hmod%Astar          
          WRITE(*,fmt='(A30,F10.5)') 'log10(T_WHIM) [K]:', log10(hmod%Twhim)          
          WRITE(*,fmt='(A30,F10.5)') 'c*:', hmod%cstar
          WRITE(*,fmt='(A30,F10.5)') 'eta:', hmod%eta
          WRITE(*,fmt='(A30,F10.5)') 'iso beta:', hmod%ibeta
       END IF
       WRITE(*,fmt='(A30,F10.5)') 'sigma*:', hmod%sstar
       WRITE(*,fmt='(A30,F10.5)') 'log10(M*) [Msun/h]:', log10(hmod%Mstar)
       WRITE(*,fmt='(A30,F10.5)') 'f_cold:', hmod%fcold
       WRITE(*,fmt='(A30,F10.5)') 'f_hot:', hmod%fhot
       IF(hmod%HMx_mode==2 .OR. hmod%HMx_mode==3) THEN
          WRITE(*,fmt='(A30,F10.5)') 'alpha mass index:', hmod%alphap
          WRITE(*,fmt='(A30,F10.5)') 'beta mass index:', hmod%betap
          WRITE(*,fmt='(A30,F10.5)') 'Gammma mass index:', hmod%Gammap
          WRITE(*,fmt='(A30,F10.5)') 'iso beta mass index:', hmod%ibetap
          WRITE(*,fmt='(A30,F10.5)') 'A* mass index:', hmod%Astarp
          WRITE(*,fmt='(A30,F10.5)') 'c* mass index:', hmod%cstarp
          WRITE(*,fmt='(A30,F10.5)') 'iso beta mass index:', hmod%ibetap
       END IF
       IF(hmod%HMx_mode==3) THEN
          WRITE(*,fmt='(A30,F10.5)') 'alpha z index:', hmod%alphaz
          WRITE(*,fmt='(A30,F10.5)') 'beta z index:', hmod%betaz
          WRITE(*,fmt='(A30,F10.5)') 'epsilon z index:', hmod%epsz
          WRITE(*,fmt='(A30,F10.5)') 'Gammma z index:', hmod%Gammaz
          WRITE(*,fmt='(A30,F10.5)') 'log10(M0) z index:', hmod%M0z
          WRITE(*,fmt='(A30,F10.5)') 'A* z index:', hmod%Astarz
          WRITE(*,fmt='(A30,F10.5)') 'T_WHIM z index:', hmod%Twhimz
          WRITE(*,fmt='(A30,F10.5)') 'c* z index:', hmod%cstarz
          WRITE(*,fmt='(A30,F10.5)') 'M* z index:', hmod%Mstarz
          WRITE(*,fmt='(A30,F10.5)') 'iso beta z index:', hmod%ibetaz
       END IF
       IF(hmod%HMx_mode==4) THEN
          WRITE(*,fmt='(A30,F10.5)') 'log10(T_heat) [K]:', log10(hmod%Theat)
          WRITE(*,fmt='(A30,F10.5)') 'alpha:', HMx_alpha(hmod%Mh,hmod)
          WRITE(*,fmt='(A30,F10.5)') 'beta:', HMx_beta(hmod%Mh,hmod)
          WRITE(*,fmt='(A30,F10.5)') 'epsilon:', HMx_eps(hmod)
          WRITE(*,fmt='(A30,F10.5)') 'Gamma:', HMx_Gamma(hmod%Mh,hmod)
          WRITE(*,fmt='(A30,F10.5)') 'log10(M0) [Msun/h]:', log10(HMx_M0(hmod))
          WRITE(*,fmt='(A30,F10.5)') 'A*:', HMx_Astar(hmod%Mh,hmod)
          WRITE(*,fmt='(A30,F10.5)') 'log10(T_WHIM) [K]:', log10(HMx_Twhim(hmod))
          WRITE(*,fmt='(A30,F10.5)') 'c*:', HMx_cstar(hmod%Mh,hmod)
       END IF
       WRITE(*,*) '======================================='
       WRITE(*,*) 'HALOMODEL: HOD parameters'
       WRITE(*,*) '======================================='
       WRITE(*,fmt='(A30,F10.5)') 'log10(M_halo_min) [Msun/h]:', log10(hmod%mhalo_min)
       WRITE(*,fmt='(A30,F10.5)') 'log10(M_halo_max) [Msun/h]:', log10(hmod%mhalo_max)
       WRITE(*,fmt='(A30,F10.5)') 'log10(M_HI_min) [Msun/h]:', log10(hmod%HImin)
       WRITE(*,fmt='(A30,F10.5)') 'log10(M_HI_max) [Msun/h]:', log10(hmod%HImax)
       WRITE(*,*) '======================================='
       IF(hmod%halo_DMONLY==5) WRITE(*,fmt='(A30,F10.5)') 'r_core [Mpc/h]:', hmod%rcore
       IF(hmod%dlnc .NE. 0.) WRITE(*,fmt='(A30,F10.5)') 'delta ln(c):', hmod%dlnc
       WRITE(*,*)

    END IF

  END SUBROUTINE print_halomod

  SUBROUTINE init_dewiggle(hmod,cosm)

    ! Initialise the dewiggled power spectrum
    USE fix_polynomial
    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: kv(4), pv(4), sigv, a
    REAL, ALLOCATABLE :: Pk(:), Pkraw(:), k(:)
    INTEGER :: i, nk

    IF(.NOT. ALLOCATED(cosm%log_k_plin)) STOP 'DEWIGGLE_INIT: Error, P(k) needs to be tabulated for this to work'

    nk=cosm%n_plin
    hmod%n_pdamp=nk
    sigv=hmod%sigv
    a=hmod%a

    ! Allocate a raw array so as to plot the processed and unprocessed results
    ALLOCATE(Pkraw(nk), Pk(nk), k(nk))

    ! Allocate the internal arrays from the cosmology arrays
    k=exp(cosm%log_k_plin)
    Pk=exp(cosm%log_plin)

    ! Fixed k values - CAMB!
    kv(1)=0.008
    kv(2)=0.01
    kv(3)=0.8
    kv(4)=1.0

    ! Fix p values
    DO i=1,4
       pv(i)=exp(find(log(kv(i)),log(k),log(Pk),nk,3,3,2))      
    END DO

    ! Create new 'raw' spectrum which has no wiggles
    DO i=1,nk
       IF(k(i)<=kv(1) .OR. k(i)>=kv(4)) THEN
          Pkraw(i)=Pk(i)
       ELSE
          Pkraw(i)=exp(Lagrange_polynomial(log(k(i)),3,log(kv),log(pv)))
       END IF
    END DO

    ! Isolate just the wiggles
    Pk=Pk-Pkraw

    ! Damp the wiggles
    Pk=Pk*exp(-(sigv*k)**2)

    ! Add the damped wiggles back in
    Pk=Pk+Pkraw

    ! Create the damped power array
    IF(ALLOCATED(hmod%log_k_pdamp)) DEALLOCATE(hmod%log_k_pdamp)
    IF(ALLOCATED(hmod%log_pdamp))   DEALLOCATE(hmod%log_pdamp)
    ALLOCATE(hmod%log_k_pdamp(nk),hmod%log_pdamp(nk))

    ! Fill the k array
    hmod%log_k_pdamp=log(k)

    ! Grow damped power to the correct redshift and fill array
    hmod%log_pdamp=log(Pk*grow(a,cosm)**2)

    ! Set the flag
    hmod%has_dewiggle=.TRUE.

  END SUBROUTINE init_dewiggle

  REAL FUNCTION p_dewiggle(k,hmod,cosm)

    ! Call the dewiggled power spectrum
    IMPLICIT NONE
    REAL, INTENT(IN) :: k
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(hmod%has_dewiggle .EQV. .FALSE.) CALL init_dewiggle(hmod,cosm)
    p_dewiggle=exp(find(log(k),hmod%log_k_pdamp,hmod%log_pdamp,hmod%n_pdamp,3,3,2))

  END FUNCTION p_dewiggle

  FUNCTION halo_type(i)

    ! Name function for halo types
    ! TODO: This must be able to be combined with set_halo_type
    ! TODO: Can this be in the header? Is this even necessary?
    IMPLICIT NONE
    CHARACTER(len=256) :: halo_type
    INTEGER :: i

    halo_type=''
    IF(i==field_dmonly)             halo_type='DMONLY'
    IF(i==field_matter)             halo_type='Matter'
    IF(i==field_cdm)                halo_type='CDM'
    IF(i==field_gas)                halo_type='Gas'
    IF(i==field_star)               halo_type='Star'
    IF(i==field_bound_gas)          halo_type='Bound gas'
    IF(i==field_free_gas)           halo_type='Free gas'
    IF(i==field_electron_pressure)  halo_type='Electron pressure'
    IF(i==field_void)               halo_type='Void'
    IF(i==field_compensated_void)   halo_type='Compensated void'
    IF(i==field_central_galaxies)   halo_type='Central galaxies/haloes'
    IF(i==field_satellite_galaxies) halo_type='Satellite galaxies'
    IF(i==field_galaxies)           halo_type='Galaxies'
    IF(i==field_HI)                 halo_type='HI'
    IF(i==field_cold_gas)           halo_type='Cold gas'
    IF(i==field_hot_gas)            halo_type='Hot gas'
    IF(i==field_static_gas)         halo_type='Static gas'
    IF(i==field_central_stars)      halo_type='Central stars'
    IF(i==field_satellite_stars)    halo_type='Satellite stars'
    IF(i==field_CIB_353)            halo_type='CIB 353 GHz'
    IF(i==field_CIB_545)            halo_type='CIB 545 GHz'
    IF(i==field_CIB_857)            halo_type='CIB 857 GHz'
    IF(i==field_halo_11p0_11p5)     halo_type='Haloes (10^11.0 -> 10^11.5)'
    IF(i==field_halo_11p5_12p0)     halo_type='Haloes (10^11.5 -> 10^12.0)'
    IF(i==field_halo_12p0_12p5)     halo_type='Haloes (10^12.0 -> 10^12.5)'
    IF(i==field_halo_12p5_13p0)     halo_type='Haloes (10^12.5 -> 10^13.0)'
    IF(i==field_halo_13p0_13p5)     halo_type='Haloes (10^13.0 -> 10^13.5)'
    IF(i==field_halo_13p5_14p0)     halo_type='Haloes (10^13.5 -> 10^14.0)'
    IF(i==field_halo_14p0_14p5)     halo_type='Haloes (10^14.0 -> 10^14.5)'
    IF(i==field_halo_14p5_15p0)     halo_type='Haloes (10^14.5 -> 10^15.0)'
    IF(halo_type=='') STOP 'HALO_TYPE: Error, i not specified correctly'

  END FUNCTION halo_type

  SUBROUTINE set_halo_type(ip)

    ! Set the halo types
    ! TODO: This must be able to be combined with halo_type
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: ip
    INTEGER :: i

    WRITE(*,*) 'SET_HALO_TYPE: Choose halo type'
    WRITE(*,*) '==============================='
    DO i=i1_fields,i2_fields
       WRITE(*,fmt='(I3,A3,A26)') i, '- ', TRIM(halo_type(i))
    END DO
    READ(*,*) ip
    WRITE(*,*) '==============================='
    WRITE(*,*)

    IF(ip<i1_fields .OR. ip>i2_fields) STOP 'SET_HALO_TYPE: Error, you have chosen a bad halo'

  END SUBROUTINE set_halo_type

  SUBROUTINE calculate_HMcode(k,a,Pk,nk,na,cosm)

    ! Get the HMcode prediction for a cosmology for a range of k and a
    IMPLICIT NONE
    REAL, INTENT(IN) :: k(nk)              ! Array of wavenumbers [h/Mpc]
    REAL, INTENT(IN) :: a(na)              ! Array of scale factors
    REAL, INTENT(OUT) :: Pk(nk,na)         ! Output power array, note that this is Delta^2(k), not P(k)
    INTEGER, INTENT(IN) :: nk              ! Number of wavenumbers
    INTEGER, INTENT(IN) :: na              ! Number of scale factors
    TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology
    REAL :: Pkk(nk)
    INTEGER :: j

    DO j=1,na
       CALL calculate_HMcode_a(k,a(j),Pkk,nk,cosm)
       Pk(:,j)=Pkk
    END DO
    
  END SUBROUTINE calculate_HMcode

  SUBROUTINE calculate_HMcode_a(k,a,Pk,nk,cosm)

    ! Get the HMcode prediction at this z for this cosmology
    IMPLICIT NONE
    REAL, INTENT(IN) :: k(nk)              ! Array of wavenumbers [h/Mpc]
    REAL, INTENT(IN) :: a                  ! Array of scale factors
    REAL, INTENT(OUT) :: Pk(nk)            ! Output power array, note that this is Delta^2(k), not P(k)
    INTEGER, INTENT(IN) :: nk              ! Number of wavenumbers
    TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology
    REAL :: pow_lin(nk), pow_2h(nk), pow_1h(nk)
    TYPE(halomod) :: hmod

    INTEGER :: ihm=1 ! Set HMcode, could/should be parameter
    INTEGER, PARAMETER :: dmonly(1)=field_dmonly ! DMONLY

    ! Do an HMcode run
    CALL assign_halomod(ihm,hmod,verbose=.FALSE.)
    CALL init_halomod(mmin_HMcode,mmax_HMcode,a,hmod,cosm,verbose=.FALSE.)
    CALL calculate_HMx_a(dmonly,1,k,nk,pow_lin,pow_2h,pow_1h,Pk,hmod,cosm,verbose=.FALSE.,response=.FALSE.)

  END SUBROUTINE calculate_HMcode_a

  SUBROUTINE calculate_HMx(ifield,nf,mmin,mmax,k,nk,a,na,pow_li,pow_2h,pow_1h,pow_hm,hmod,cosm,verbose,response)

    ! Public facing function, calculates the halo model power for k and a range
    ! TODO: Change (:,:,k,a) to (k,a,:,:) for speed or (a,k,:,:)?
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ifield(nf) ! Indices for different fields
    INTEGER, INTENT(IN) :: nf         ! Number of different fields
    REAL, INTENT(IN) :: mmin          ! Minimum halo mass [Msun/h]
    REAL, INTENT(IN) :: mmax          ! Maximum halo mass [Msun/h]
    REAL, INTENT(IN) :: k(nk)         ! k array [h/Mpc]
    INTEGER, INTENT(IN) :: nk         ! Number of k points
    REAL, INTENT(IN) :: a(na)         ! a array 
    INTEGER, INTENT(IN) :: na         ! Number of a points
    REAL, ALLOCATABLE, INTENT(OUT) :: pow_li(:,:)     ! Pow(k,a)
    REAL, ALLOCATABLE, INTENT(OUT) :: pow_2h(:,:,:,:) ! Pow(f1,f2,k,a)
    REAL, ALLOCATABLE, INTENT(OUT) :: pow_1h(:,:,:,:) ! Pow(f1,f2,k,a)
    REAL, ALLOCATABLE, INTENT(OUT) :: pow_hm(:,:,:,:) ! Pow(f1,f2,k,a)
    TYPE(halomod), INTENT(INOUT) :: hmod   ! Halo model
    TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology
    LOGICAL, INTENT(IN) :: verbose
    LOGICAL, INTENT(IN) :: response
    REAL :: z
    INTEGER :: i, j
    LOGICAL :: verbose2

    ! To avoid splurge of stuff printed to screen
    verbose2=verbose

    ! Deallocate arrays if they are already allocated
    IF(ALLOCATED(pow_li)) DEALLOCATE(pow_li)
    IF(ALLOCATED(pow_2h)) DEALLOCATE(pow_2h)
    IF(ALLOCATED(pow_1h)) DEALLOCATE(pow_1h)
    IF(ALLOCATED(pow_hm)) DEALLOCATE(pow_hm)

    ! Allocate power arrays
    ALLOCATE(pow_li(nk,na),pow_2h(nf,nf,nk,na),pow_1h(nf,nf,nk,na),pow_hm(nf,nf,nk,na))

    ! Do the halo-model calculation by looping over scale factor index
    DO i=na,1,-1
       
       z=redshift_a(a(i))
       CALL init_halomod(mmin,mmax,a(i),hmod,cosm,verbose2)
       CALL print_halomod(hmod,cosm,verbose2)
       CALL calculate_HMx_a(ifield,nf,k,nk,pow_li(:,i),pow_2h(:,:,:,i),pow_1h(:,:,:,i),pow_hm(:,:,:,i),hmod,cosm,verbose2,response)
       
       IF(i==na .and. verbose) THEN
          WRITE(*,*) 'CALCULATE_HMx: Doing calculation'
          DO j=1,nf
             WRITE(*,*) 'CALCULATE_HMx: Haloes:', ifield(j), TRIM(halo_type(ifield(j)))
          END DO
          WRITE(*,*) '======================================='
          WRITE(*,*) '                            a         z'
          WRITE(*,*) '======================================='
       END IF
       
       IF(verbose) WRITE(*,fmt='(A15,I5,2F10.3)') 'CALCULATE_HMx:', i, a(i), z
       verbose2=.FALSE.
       
    END DO
    
    IF(verbose) THEN
       WRITE(*,*) '======================================='
       WRITE(*,*) 'CALCULATE_HMx: Done'
       WRITE(*,*)
    END IF

  END SUBROUTINE calculate_HMx

  SUBROUTINE calculate_HMx_a(ifield,nf,k,nk,pow_li,pow_2h,pow_1h,pow_hm,hmod,cosm,verbose,response)

    ! Calculate halo model Pk(a) for a k range
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ifield(nf)
    INTEGER, INTENT(IN) :: nf
    REAL, INTENT(IN) :: k(nk)
    INTEGER, INTENT(IN) :: nk
    REAL, INTENT(OUT) :: pow_li(nk)
    REAL, INTENT(OUT) :: pow_2h(nf,nf,nk)
    REAL, INTENT(OUT) :: pow_1h(nf,nf,nk)
    REAL, INTENT(OUT) :: pow_hm(nf,nf,nk)
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    LOGICAL, INTENT(IN) :: verbose
    LOGICAL, INTENT(IN) :: response
    REAL :: plin
    REAL :: powg_2h(nk), powg_1h(nk), powg_hm(nk)
    REAL, ALLOCATABLE :: upow_2h(:,:,:), upow_1h(:,:,:), upow_hm(:,:,:)
    REAL :: hmcode_2h(nk), hmcode_1h(nk), hmcode_hm(nk)
    INTEGER :: i, j, ii, jj, ihmcode, match(nf), nnf
    INTEGER, ALLOCATABLE :: iifield(:)
    TYPE(halomod) :: hmcode    
    INTEGER, PARAMETER :: dmonly(1)=field_dmonly ! Needed because it needs to be an array(1)

    ! Make a new indexing scheme for only the unique arrays
    CALL unique_index(ifield,nf,iifield,nnf,match)
    ALLOCATE(upow_2h(nnf,nnf,nk),upow_1h(nnf,nnf,nk),upow_hm(nnf,nnf,nk))

    ! Write to screen
    IF(verbose) THEN
       DO i=1,nnf
          WRITE(*,*) 'CALCULATE_HMX_A: Halo type:', iifield(i), TRIM(halo_type(iifield(i)))
       END DO
       WRITE(*,*) 'CALCULATE_HMX_A: k min [h/Mpc]:', REAL(k(1))
       WRITE(*,*) 'CALCULATE_HMX_A: k max [h/Mpc]:', REAL(k(nk))
       WRITE(*,*) 'CALCULATE_HMX_A: number of k:', nk
       WRITE(*,*) 'CALCULATE_HMX_A: a:', REAL(hmod%a)
       WRITE(*,*) 'CALCULATE_HMX_A: z:', REAL(hmod%z)
       WRITE(*,*) 'CALCULATE_HMX_A: Calculating halo-model power spectrum'
       WRITE(*,*)
    END IF

    ! Do an HMcode calculation for multiplying the response
    IF(hmod%response==1 .OR. hmod%response==2) THEN
       ihmcode=1
       CALL assign_halomod(ihmcode,hmcode,verbose=.FALSE.)
       CALL init_halomod(mmin_HMx,mmax_HMx,hmod%a,hmcode,cosm,verbose=.FALSE.)
    END IF

    ! Loop over k values
    ! TODO: add OMP support properly. What is private and what is shared? CHECK THIS!
!!$OMP PARALLEL DO DEFAULT(SHARED)!, private(k,plin,pow_2h,pow_1h,pow,pow_lin)
!!$OMP PARALLEL DO DEFAULT(PRIVATE)
!!$OMP PARALLEL DO FIRSTPRIVATE(nk,cosm,compute_p_lin,k,a,pow_lin,plin,itype1,itype2,z,pow_2h,pow_1h,pow,hmod)
!!$OMP PARALLEL DO
    DO i=1,nk

       ! Get the linear power
       plin=p_lin(k(i),hmod%a,cosm)
       pow_li(i)=plin

       ! Do the halo model calculation
       ! TODO: slow array accessing
       CALL calculate_HMx_ka(iifield,nnf,k(i),plin,upow_2h(:,:,i),upow_1h(:,:,i),upow_hm(:,:,i),hmod,cosm)

       IF(response) THEN

          ! If doing a response then calculate a DMONLY prediction too
          CALL calculate_HMx_ka(dmonly,1,k(i),plin,powg_2h(i),powg_1h(i),powg_hm(i),hmod,cosm)
          pow_li(i)=1.                             ! This is just linear-over-linear, which is one
          upow_2h(:,:,i)=upow_2h(:,:,i)/powg_2h(i) ! Two-halo response (slow array accessing)
          upow_1h(:,:,i)=upow_1h(:,:,i)/powg_1h(i) ! One-halo response (slow array accessing)
          upow_hm(:,:,i)=upow_hm(:,:,i)/powg_hm(i) ! Full model response (slow array accessing)
          
       ELSE IF(hmod%response==1 .OR. hmod%response==2) THEN

          ! If doing a response then calculate a DMONLY prediction too
          CALL calculate_HMx_ka(dmonly,1,k(i),plin,powg_2h(i),powg_1h(i),powg_hm(i),hmod,cosm)
          CALL calculate_HMx_ka(dmonly,1,k(i),plin,hmcode_2h(i),hmcode_1h(i),hmcode_hm(i),hmcode,cosm)

          IF(hmod%response==1) THEN

             ! Do the response for everything
             upow_2h(:,:,i)=upow_2h(:,:,i)*hmcode_2h(i)/powg_2h(i) ! Two-halo response times HMcode (slow array accessing)
             upow_1h(:,:,i)=upow_1h(:,:,i)*hmcode_1h(i)/powg_1h(i) ! One-halo response times HMcode (slow array accessing)
             upow_hm(:,:,i)=upow_hm(:,:,i)*hmcode_hm(i)/powg_hm(i) ! Full model response times HMcode (slow array accessing)

          ELSE IF(hmod%response==2) THEN

             ! Exclude pressure from the response
             DO ii=1,nf
                DO jj=1,nf
                   IF((iifield(ii) .NE. field_electron_pressure) .AND. (iifield(jj) .NE. field_electron_pressure)) THEN 
                      upow_2h(ii,jj,i)=upow_2h(ii,jj,i)*hmcode_2h(i)/powg_2h(i) ! Two-halo response times HMcode (slow accessing)
                      upow_1h(ii,jj,i)=upow_1h(ii,jj,i)*hmcode_1h(i)/powg_1h(i) ! One-halo response times HMcode (slow accessing)
                      upow_hm(ii,jj,i)=upow_hm(ii,jj,i)*hmcode_hm(i)/powg_hm(i) ! Full model response times HMcode (slow accessing)
                   END IF
                END DO
             END DO
             
          END IF

       END IF

    END DO
!!$OMP END PARALLEL DO

    ! Now fill the full arrays from the unique arrays
    DO i=1,nf
       DO j=1,nf
          ii=match(i)
          jj=match(j)
          pow_1h(i,j,:)=upow_1h(ii,jj,:)
          pow_2h(i,j,:)=upow_2h(ii,jj,:)
          pow_hm(i,j,:)=upow_hm(ii,jj,:)
       END DO
    END DO

  END SUBROUTINE calculate_HMx_a

  SUBROUTINE calculate_HMx_ka(ifield,nf,k,plin,pow_2h,pow_1h,pow_hm,hmod,cosm)

    ! Gets the one- and two-halo terms and combines them
    ! TODO: Include scatter in two-halo term
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ifield(nf)
    INTEGER, INTENT(IN) :: nf
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: plin
    REAL, INTENT(OUT) :: pow_2h(nf,nf)
    REAL, INTENT(OUT) :: pow_1h(nf,nf)
    REAL, INTENT(OUT) :: pow_hm(nf,nf)
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: wk(hmod%n,nf), wk2(hmod%n,2), wk_product(hmod%n)
    INTEGER :: i, j, f1, f2, ih(2)

    ! Calls expressions for two- and one-halo terms and then combines to form the full power spectrum
    IF(k==0.) THEN

       ! This should really never be called for k=0
       pow_2h=0.
       pow_1h=0.

    ELSE

       ! Get the window functions
       CALL init_windows(k,ifield,nf,wk,hmod%n,hmod,cosm)

       ! Loop over fields and get the one-halo term for each pair
       DO f1=1,nf
          DO f2=f1,nf
             IF(hmod%dlnc==0.) THEN
                wk_product=wk(:,f1)*wk(:,f2)
             ELSE
                ih(1)=ifield(f1)
                ih(2)=ifield(f2)
                wk_product=wk_product_scatter(hmod%n,ih,k,hmod,cosm)
             END IF
             pow_1h(f1,f2)=p_1h(wk_product,hmod%n,k,hmod,cosm)            
          END DO
       END DO

       ! If linear theory is used for two-halo term we recalculate the window functions for the two-halo term with k=0 fixed
       IF(hmod%ip2h==1 .OR. hmod%ip2h==3) THEN
          CALL init_windows(0.,ifield,nf,wk,hmod%n,hmod,cosm)      
       ELSE IF(hmod%halo_free_gas==7) THEN
          ! If we are worrying about unbound gas then we need this
          CALL add_smooth_component_to_windows(ifield,nf,wk,hmod%n,hmod,cosm)
       END IF

       ! Get the two-halo term
       DO i=1,nf
          DO j=i,nf
             ih(1)=ifield(i)
             ih(2)=ifield(j)
             wk2(:,1)=wk(:,i)
             wk2(:,2)=wk(:,j)
             pow_2h(i,j)=p_2h(ih,wk2,hmod%n,k,plin,hmod,cosm)
          END DO
       END DO

    END IF

    ! Loop over fields and get the total halo-model power
    DO i=1,nf
       DO j=i,nf
          pow_hm(i,j)=p_hm(k,pow_2h(i,j),pow_1h(i,j),hmod,cosm)
       END DO
    END DO

    ! Construct symmetric parts using ij=ji symmetry of spectra
    DO i=1,nf
       DO j=i,nf
          pow_1h(j,i)=pow_1h(i,j)
          pow_2h(j,i)=pow_2h(i,j)
          pow_hm(j,i)=pow_hm(i,j)
       END DO
    END DO
    
  END SUBROUTINE calculate_HMx_ka

  SUBROUTINE init_windows(k,fields,nf,wk,nm,hmod,cosm)

    ! Fill the window functions for all the different fields
    IMPLICIT NONE
    REAL, INTENT(IN) :: k
    INTEGER, INTENT(IN) :: fields(nf)
    REAL, INTENT(OUT) :: wk(nm,nf)
    INTEGER, INTENT(IN) :: nf
    INTEGER, INTENT(IN) :: nm
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: i, j
    REAL :: m, rv, c, rs, nu, et
    INTEGER :: i_all, i_cdm, i_gas, i_sta
    LOGICAL :: quick_matter

    IF(repeated_entries(fields,nf)) STOP 'INIT_WINDOWS: Error, repeated fields'

    ! This should be set to false initially
    quick_matter=.FALSE.

    ! Get the array positions corresponding to all, cdm, gas, stars if they exist
    i_all=array_position(field_matter,fields,nf)
    i_cdm=array_position(field_cdm,fields,nf)
    i_gas=array_position(field_gas,fields,nf)
    i_sta=array_position(field_star,fields,nf)

    ! If all, cdm, gas and stars exist then activate the quick-matter mode
    IF((i_all .NE. 0) .AND. (i_cdm .NE. 0) .AND. (i_gas .NE. 0) .AND. (i_sta .NE. 0)) THEN
       quick_matter=.TRUE.
    END IF

    ! Get eta
    et=eta_HMcode(hmod,cosm)

    ! Calculate the halo window functions for each field
    DO j=1,nf

       IF(quick_matter .AND. j==i_all) CYCLE

       ! Loop over masses to fill window-function array
       DO i=1,nm

          m=hmod%m(i)
          rv=hmod%rv(i)
          c=hmod%c(i)
          rs=rv/c
          nu=hmod%nu(i)

          wk(i,j)=win_type(.FALSE.,fields(j),k*nu**et,m,rv,rs,hmod,cosm)

       END DO

    END DO

    ! If quick-matter mode is active then create the total matter window by summing contributions
    IF(quick_matter) THEN
       wk(:,i_all)=wk(:,i_cdm)+wk(:,i_gas)+wk(:,i_sta)
    END IF

  END SUBROUTINE init_windows

  SUBROUTINE add_smooth_component_to_windows(fields,nf,wk,nm,hmod,cosm)

    ! Refills the window functions for the two-halo term if this is necessary
    ! This is for contributions due to unbound gas, and the effect of this on electron pressure
    ! TODO: Have I inculded the halo bias corresponding to the free component correctly?
    ! BUG: Fixed conversion of thermal to electron pressure
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: fields(nf)
    REAL, INTENT(INOUT) :: wk(nm,nf)
    INTEGER, INTENT(IN) :: nf
    INTEGER, INTENT(IN) :: nm
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: m, rv, c, rs, nu, et, fc, pc
    REAL :: rho0, T0
    INTEGER :: i, i_all, i_gas, i_pre

    ! Check for repeated fields
    IF(repeated_entries(fields,nf)) STOP 'ADD_SMOOTH_COMPONENT_TO_WINDOWS: Error, repeated fields'

    ! Get the array positions corresponding to all, cdm, gas, stars if they existza
    i_all=array_position(field_matter,fields,nf)
    i_gas=array_position(field_gas,fields,nf)
    i_pre=array_position(field_electron_pressure,fields,nf)

    IF(i_all==0 .AND. i_gas==0 .AND. i_pre==0) THEN

       ! Do nothing

    ELSE

       ! Get eta
       et=eta_HMcode(hmod,cosm)

       ! Loop over mass and apply corrections
       DO i=1,hmod%n

          ! Halo variables
          m=hmod%m(i)
          rv=hmod%rv(i)
          c=hmod%c(i)
          rs=rv/c
          nu=hmod%nu(i)

          ! Correction factor for the gas density profiles
          fc=halo_free_gas_fraction(m,hmod,cosm)*m/comoving_matter_density(cosm)
          IF(i_all .NE. 0) wk(i,i_all)=wk(i,i_all)+fc
          IF(i_gas .NE. 0) wk(i,i_gas)=wk(i,i_gas)+fc
          
          IF(i_pre .NE. 0) THEN

             ! Calculate the value of the density profile prefactor [(Msun/h)/(Mpc/h)^3] and change units from cosmological to SI
             rho0=m*halo_free_gas_fraction(m,hmod,cosm) ! rho0 in [(Msun/h)/(Mpc/h)^3]
             rho0=rho0*msun/Mpc/Mpc/Mpc                 ! Overflow with REAL(4) if you use Mpc**3, this converts to SI units [h^2 kg/m^3]
             rho0=rho0*cosm%h**2                        ! Absorb factors of h, so now [kg/m^3]

             ! This is the total thermal pressure of the WHIM
             T0=HMx_Twhim(hmod) ! [K]

             ! Factors to convert from Temp x density -> electron pressure (Temp x n; n is all particle number density) 
             pc=(rho0/(mp*cosm%mup))*(kb*T0) ! Multiply window by *number density* (all particles) times temperature time k_B [J/m^3]
             pc=pc/(eV*(0.01)**(-3)) ! Change units to pressure in [eV/cm^3]
             pc=pc*cosm%mup/cosm%mue ! Convert from total thermal pressure to electron pressure
             
             ! Add correction to 'electron pressure' haloes
             wk(i,i_pre)=wk(i,i_pre)+pc

          END IF

       END DO

    END IF

  END SUBROUTINE add_smooth_component_to_windows

  FUNCTION wk_product_scatter(n,ifield,k,hmod,cosm)

    IMPLICIT NONE
    REAL :: wk_product_scatter(n)
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: ifield(2)
    REAL, INTENT(IN) :: k
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: i
    REAL :: m, rv, c

    DO i=1,n
       m=hmod%m(i)
       rv=hmod%rv(i)
       c=hmod%c(i)
       wk_product_scatter(i)=integrate_scatter(c,hmod%dlnc,ifield,k,m,rv,hmod,cosm,hmod%acc_HMx,3)
    END DO

  END FUNCTION wk_product_scatter

  REAL FUNCTION p_2h(ih,wk,n,k,plin,hmod,cosm)

    ! Produces the 'two-halo' power
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ih(2)
    REAL, INTENT(IN) :: wk(n,2)
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: plin
    TYPE(halomod), INTENT(INOUT) :: hmod    
    TYPE(cosmology), INTENT(INOUT) :: cosm    
    REAL :: sigv, frac, rhom
    REAL :: m0, b0, wk0(2), nu0
    REAL :: m1, m2, rv1, rv2, b1, b2, g1, g2, u1, u2, nu1, nu2
    REAL :: Inl, Inl_11, Inl_21, Inl_12, Inl_22, B_NL
    REAL :: I_11, I_12(n), I_21(n), I_22(n,n)
    REAL :: I2h, I2hs(2)
    INTEGER :: i, j
    LOGICAL, PARAMETER :: add_missing_bits=.TRUE.

    ! Necessary to prevent warning for some reason
    I_11=0.

    rhom=comoving_matter_density(cosm)

    IF(hmod%ip2h==1) THEN

       ! Simply linear theory
       p_2h=plin

    ELSE IF(hmod%ip2h==3) THEN

       ! Damped BAO linear theory
       p_2h=p_dewiggle(k,hmod,cosm)

    ELSE IF(hmod%ip2h==4) THEN

       ! No two-halo term
       p_2h=0.

    ELSE

       IF(hmod%imf==4) THEN

          ! In this case the mass function is a delta function...
          m0=hmod%hmass
          nu0=nu_M(m0,hmod,cosm)
          b0=b_nu(nu0,hmod)
          DO j=1,2
             wk0(j)=find(log(m0),hmod%log_m,wk(:,j),n,3,3,2)
             I2hs(j)=b0*rhom*wk0(j)/m0
          END DO

       ELSE

          ! ... otherwise we need to do an integral
          DO j=1,2
             CALL I_2h(ih(j),I2h,wk(:,j),n,hmod,cosm,ibias=1)
             I2hs(j)=I2h
          END DO
          
       END IF

       p_2h=plin*I2hs(1)*I2hs(2)

       ! Second-order bias correction
       ! This needs to have the property that \int f(nu)b2(nu) du = 0
       ! This means it is hard to check that the normalisation is correct
       ! e.g., how much do low mass haloes matter
       IF(hmod%ibias==2) THEN

          ! ... otherwise we need to do an integral
          DO j=1,2
             CALL I_2h(ih(j),I2h,wk(:,j),n,hmod,cosm,ibias=2)
             I2hs(j)=I2h
          END DO
          p_2h=p_2h+(plin**2)*I2hs(1)*I2hs(2)*rhom**2 ! TODO: Should rhom^2 really be here?

       ELSE IF(hmod%ibias==3) THEN

          ! Full non-linear bias correction
          DO j=1,n

             m2=hmod%m(j)
             nu2=hmod%nu(j)
             rv2=hmod%rv(j)
             b2=b_nu(nu2,hmod)
             g2=g_nu(nu2,hmod)
             u2=rhom*wk(j,2)/m2

             DO i=1,n

                m1=hmod%m(i)
                nu1=hmod%nu(i)
                rv1=hmod%rv(i)
                b1=b_nu(nu1,hmod)        
                g1=g_nu(nu1,hmod)               
                u1=rhom*wk(i,1)/m1

                !F(i,j)=B_NL(nu1,nu2,k,hmod%z)*b1*b2*g1*g2*u1*u2
                !F(i,j)=B_NL(k,nu1,nu2,rv1,rv2)*b1*b2*g1*g2*u1*u2
                B_NL=BNL(k,nu1,nu2,hmod)

                IF(i==1 .AND. j==1) THEN
                   I_11=B_NL*u1*u2
                END IF

                ! Integrand for upper-left quadrant. Only single integral over nu2. Only nu2 properties varying.
                IF(i==1) THEN
                   I_12(j)=B_NL*b2*g2*u2
                END IF

                ! Integrand for bottom-right quadrant. Only single integral over nu1. Only nu1 properties varying.
                IF(j==1) THEN
                   I_21(i)=B_NL*b1*g1*u1
                END IF

                ! Integrand for upper-right quadrant. Main double integral over nu1 and nu2
                I_22(i,j)=B_NL*b1*b2*g1*g2*u1*u2

             END DO
             
          END DO          
          
          !Inl_11=(wk(1,1)*rhom/hmod%m(1))*(wk(1,2)*rhom/hmod%m(1))*BNL(k,hmod%nu(1),hmod%nu(1),hmod%rv(1),hmod%rv(1),hmod)*hmod%gbmin**2
          Inl_11=I_11*hmod%gbmin**2
          Inl_12=integrate_table(hmod%nu,I_12,n,1,n,iorder=1)*hmod%gbmin*(wk(1,1)*rhom/hmod%m(1))
          Inl_21=integrate_table(hmod%nu,I_21,n,1,n,iorder=1)*hmod%gbmin*(wk(1,2)*rhom/hmod%m(1))
          Inl_22=integrate_table(hmod%nu,hmod%nu,I_22,n,n)

          !WRITE(*,*) 'B_NL:', k, Inl_11, Inl_12, Inl_21, Inl_22

          IF(add_missing_bits) THEN
             Inl=Inl_11+Inl_21+Inl_12+Inl_22
          ELSE
             Inl=Inl_22
          END IF

          p_2h=p_2h+plin*Inl
          
       ELSE IF(hmod%ibias==4) THEN

          ! Fedeli (2014b) 'scale-dependent halo bias'
          p_2h=p_2h*(1.+k)**0.54

       ELSE IF(hmod%ibias==5) THEN

          ! My own experimental model
          p_2h=p_2h*(1.+(k/hmod%knl))**0.5

       END IF

    END IF

    ! Get the normalisation correct for non-matter fields if you are using linear theory or damped BAO
    IF(hmod%ip2h==1 .OR. hmod%ip2h==3) THEN
       DO j=1,2
          IF(ih(j)==field_dmonly .OR. ih(j)==field_matter) THEN
             I2hs(j)=1.
          ELSE
             CALL I_2h(ih(j),I2h,wk(:,j),n,hmod,cosm,ibias=1)
             I2hs(j)=I2h
          END IF
       END DO
       p_2h=p_2h*I2hs(1)*I2hs(2)
    END IF

    ! Apply the damping to the two-halo term
    IF(hmod%i2hdamp .NE. 1) THEN
       ! Two-halo damping parameters
       sigv=hmod%sigv
       frac=fdamp_HMcode(hmod,cosm)
       IF(frac==0.) THEN
          p_2h=p_2h
       ELSE
          p_2h=p_2h*(1.-frac*(tanh(k*sigv/sqrt(ABS(frac))))**2)
       END IF
    END IF

!!$    ! For some extreme cosmologies frac>1. so this must be added to prevent p_2h<0
!!$    IF(p_2h<0.) THEN
!!$       WRITE(*,*) 'P_2H: Halo type 1:', ih(1)
!!$       WRITE(*,*) 'P_2H: Halo type 1:', ih(2)
!!$       WRITE(*,*) 'P_2H: k [h/Mpc]:', k
!!$       WRITE(*,*) 'P_2H: z:', hmod%z
!!$       WRITE(*,*) 'P_2H: Delta^2_{2H}:', p_2h
!!$       WRITE(*,*) 'P_2H: Caution! P_2h < 0, this was previously fixed by setting P_2h = 0 explicitly'
!!$       !p_2h=0.
!!$       STOP
!!$    END IF

  END FUNCTION p_2h

  SUBROUTINE I_2h(ih,int,wk,n,hmod,cosm,ibias)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ih
    REAL, INTENT(OUT) :: int
    REAL, INTENT(IN) :: wk(n)
    INTEGER, INTENT(IN) :: n
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER, INTENT(IN) :: ibias
    REAL :: rhom, b, m0, m, nu, integrand(n)
    INTEGER :: i

    rhom=comoving_matter_density(cosm)

    ! ...otherwise you need to do an integral
    
    DO i=1,n

       ! Some variables to make equations cleaner below
       m=hmod%m(i)
       nu=hmod%nu(i)

       ! Linear bias term, standard two-halo term integral
       IF(ibias==1) THEN
          b=b_nu(nu,hmod)
       ELSE IF(ibias==2) THEN
          b=b2_nu(nu,hmod)
       ELSE
          STOP 'I_2H: Error, ibias specified incorrectly'
       END IF

       integrand(i)=g_nu(nu,hmod)*b*(rhom*wk(i)/m)

    END DO

    ! Evaluate these integrals from the tabulated values
    int=integrate_table(hmod%nu,integrand,n,1,n,3)

    IF(ibias==1) THEN

       IF(hmod%ip2h_corr==1) THEN
          ! Do nothing in this case
          ! There will be large errors if any signal is from low-mass haloes
          ! e.g., for the matter power spectrum
       ELSE IF(hmod%ip2h_corr==2) THEN
          ! Add on the value of integral b(nu)*g(nu) assuming W(k)=1
          ! Advised by Yoo et al. (????) and Cacciato et al. (2012)
          STOP 'P_2H: This will not work for fields that do not have mass fractions defined'
          int=int+hmod%gbmin*halo_fraction(ih,m,hmod,cosm)
       ELSE IF(hmod%ip2h_corr==3) THEN
          ! Put the missing part of the integrand as a delta function at the low-mass limit of the integral
          ! I think this is the best thing to do
          !IF(requal(hmod%m(1),mmin_HMx,eps_p2h)) THEN
          m0=hmod%m(1)           
          int=int+hmod%gbmin*(rhom*wk(1)/m0)
          !END IF          
       ELSE
          STOP 'P_2h: Error, ip2h_corr not specified correctly'
       END IF

    END IF
    
  END SUBROUTINE I_2h

  REAL FUNCTION p_1h(wk_product,n,k,hmod,cosm)

    ! Calculates the one-halo term
    IMPLICIT NONE
    REAL, INTENT(IN) :: wk_product(n)
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(IN) :: k
    TYPE(halomod), INTENT(INOUT) :: hmod    
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: m, g, fac, ks, wk0_product, m0, rhom, integrand(n)
    INTEGER :: i

    INTEGER, PARAMETER :: iorder=1 ! Use basic trapezium rule because the integrand has rapid oscillations in W(k)

    ! Matter density
    rhom=comoving_matter_density(cosm)

    IF(hmod%imf==4) THEN

       ! In this case the mass function is a delta function...

       m0=hmod%hmass
       wk0_product=find(log(m0),hmod%log_m,wk_product,n,3,3,2)
       p_1h=rhom*wk0_product/m0

    ELSE

       ! ...otherwise you need to do an integral

       ! Calculates the value of the integrand at all nu values!
       DO i=1,n
          g=g_nu(hmod%nu(i),hmod)
          m=hmod%m(i)
          integrand(i)=g*wk_product(i)/m
       END DO

       ! Carries out the integration   
       p_1h=rhom*integrate_table(hmod%nu,integrand,n,1,n,iorder)

    END IF

    ! Convert from P(k) -> Delta^2(k)
    p_1h=p_1h*(4.*pi)*(k/twopi)**3

    ! Damping of the 1-halo term at very large scales
    ks=kstar_HMcode(hmod,cosm)       

    IF(ks>0.) THEN

       IF(hmod%i1hdamp==1) THEN
          ! Do nothing in this case
       ELSE IF(hmod%i1hdamp==2) THEN
          IF((k/ks)**2>7.) THEN
             ! Prevents problems if k/ks is very large
             fac=0.
          ELSE
             fac=exp(-((k/ks)**2))
          END IF
          p_1h=p_1h*(1.-fac)
       ELSE IF(hmod%i1hdamp==3) THEN
          ! Note that the power here should be 4 because it multiplies Delta^2(k) ~ k^3 at low k (NOT 7)
          ! Want f(k<<ks) ~ k^4; f(k>>ks) = 1
          fac=1./(1.+(ks/k)**4)
          p_1h=p_1h*fac
       ELSE
          STOP 'P_1H: Error, i1hdamp not specified correctly'          
       END IF

    END IF

  END FUNCTION p_1h

  REAL FUNCTION p_hm(k,pow_2h,pow_1h,hmod,cosm)

    IMPLICIT NONE
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: pow_2h
    REAL, INTENT(IN) :: pow_1h
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: alpha, pow, con

    ! alpha is set to one sometimes, which is just the standard halo-model sum of terms
    ! No need to have an IF statement around this
    IF(hmod%itrans==2 .OR. hmod%itrans==3 .OR. hmod%itrans==4) THEN

       ! If either term is less than zero then we need to be careful
       IF(pow_2h<0. .OR. pow_1h<0.) THEN

          ! Either the safe option and just sum the components
          IF(hmod%safe_negative) THEN
             p_hm=pow_2h+pow_1h
          ELSE
             WRITE(*,*) 'P_HM: Two-halo term:', pow_2h
             WRITE(*,*) 'P_1h: One-halo term:', pow_1h
             STOP 'P_HM: Error, either pow_2h or pow_1h is negative, which is a problem for smoothed transition'
          END IF

       ELSE

          ! Do the standard smoothed transition
          alpha=alpha_HMcode(hmod,cosm)
          p_hm=(pow_2h**alpha+pow_1h**alpha)**(1./alpha)

       END IF

    ELSE IF(hmod%itrans==5) THEN

       ! Sigmoid transition
       con=1. ! Constant for the transition
       pow=1. ! Power for the transition
       p_hm=pow_2h+sigmoid_log(con*(1.*k/hmod%knl),pow)*(pow_1h-pow_2h)
       !p_hm=pow_2h+0.5*(1.+tanh(con*(k-hmod%knl)))*(pow_1h-pow_2h)

    ELSE

       ! Standard sum
       p_hm=pow_2h+pow_1h

    END IF

    ! If we are adding in power from voids
    IF(hmod%voids) THEN
       p_hm=p_hm+p_1void(k,hmod)
    END IF

  END FUNCTION p_hm

  REAL FUNCTION p_1void(k,hmod)

    IMPLICIT NONE
    REAL, INTENT(IN) :: k
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: dc, wk, V, rvoid, rcomp, nu
    REAL :: integrand(hmod%n)
    INTEGER :: i, n

    REAL, PARAMETER :: dv=-1.
    REAL, PARAMETER :: fvoid=1.1
    REAL, PARAMETER :: rvoid_simple=10.
    LOGICAL, PARAMETER :: compensate=.TRUE.
    LOGICAL, PARAMETER :: simple=.FALSE.

    IF(simple) THEN
       n=1
    ELSE
       n=hmod%n
    END IF

    DO i=1,n

       !Get the void radius and compensation radius
       IF(simple) THEN
          rvoid=rvoid_simple
       ELSE         
          rvoid=hmod%rr(i)
          nu=hmod%nu(i)        
       END IF
       rcomp=fvoid*rvoid

       !Calculate the compensation over-density
       dc=-dv*rvoid**3/(rcomp**3-rvoid**3)

       !Calculate the void Fourier transform
       IF(compensate) THEN
          wk=(4.*pi/3.)*((dv-dc)*wk_tophat(k*rvoid)*rvoid**3+dc*wk_tophat(k*rcomp)*rcomp**3)
       ELSE
          wk=(4.*pi/3.)*dv*wk_tophat(k*rvoid)*rvoid**3
       END IF

       !Calculate the void volume
       IF(compensate) THEN
          V=rcomp**3
       ELSE
          V=rvoid**3
       END IF

       IF(simple .EQV. .FALSE.) THEN
          integrand(i)=g_nu(nu,hmod)*wk**2/V
       END IF

    END DO

    !Calculate the void one-halo term
    IF(simple) THEN
       p_1void=wk**2/V
    ELSE
       p_1void=integrate_table(hmod%nu,integrand,n,1,n,1)
    END IF

    p_1void=p_1void*(4.*pi)*(k/twopi)**3

  END FUNCTION p_1void

!!$  REAL FUNCTION B_NL(nu1,nu2,k,z)
!!$
!!$    IMPLICIT NONE
!!$    REAL, INTENT(IN) :: nu1, nu2
!!$    REAL, INTENT(IN) :: k, z
!!$    REAL :: A, k0
!!$    REAL :: A0, A1, k00, k01
!!$
!!$    ! Set the model
!!$    INTEGER, PARAMETER :: model=3
!!$
!!$    ! Redshifts model was calibrated at
!!$    REAL, PARAMETER :: zs(4)=[0.0,0.5,1.0,2.0]
!!$
!!$    ! Do we limit from below?
!!$    LOGICAL, PARAMETER :: impose_limit=.TRUE.
!!$    REAL, PARAMETER :: limit=-1.
!!$
!!$    ! Model 1 - ?
!!$    REAL, PARAMETER :: As(4)=[3.14,2.60,2.40,1.93]
!!$    REAL, PARAMETER :: k0s(4)=[1.79,1.61,1.66,1.66]
!!$
!!$    ! Model 2 - ?
!!$    REAL, PARAMETER :: A0s(4)=[9.83,39.85,106.42,153.41]
!!$    REAL, PARAMETER :: k00s(4)=[3.38,5.14,7.88,14.03]
!!$    
!!$    IF(model==1) THEN
!!$
!!$       ! Model 1 - ?
!!$       
!!$       A=Lagrange_Polynomial(z,3,zs,As)
!!$       k0=Lagrange_Polynomial(z,3,zs,k0s)
!!$
!!$       !B_NL=A*(k/k0)*(exp(-(k/(2.*k0**2))
!!$       B_NL=A*(k/k0)*(1.-k**2/(2.*k0**2))
!!$
!!$    ELSE IF(model==2) THEN
!!$
!!$       ! Model 2 - ?
!!$
!!$       A0=Lagrange_Polynomial(z,3,zs,A0s)
!!$       A1=-3
!!$       
!!$       k00=Lagrange_Polynomial(z,3,zs,k00s)
!!$       k01=-1.5
!!$
!!$       A=A0*(nu1+nu2)**A1
!!$       k0=k00*(nu1+nu2)**k01
!!$       
!!$       B_NL=A*(k/k0)*(1.-(k**2/(2.*k0**2)))
!!$
!!$    ELSE IF(model==3) THEN
!!$
!!$       !WRITE(*,*) 'Model 3!'
!!$       !STOP
!!$       
!!$       A=4.5
!!$       k0=5.
!!$
!!$       B_NL=A*(k/k0)*(1.-(k**2/(2.*k0**2)))
!!$
!!$    ELSE
!!$
!!$       STOP 'B_NL: Error, model not specified correctly'
!!$
!!$    END IF
!!$
!!$    IF(impose_limit .AND. B_NL<limit) B_NL=limit
!!$    
!!$  END FUNCTION B_NL

!!$  REAL FUNCTION B_NL(k,nu1,nu2,rv1,rv2)
!!$
!!$    IMPLICIT NONE
!!$    REAL, INTENT(IN) :: k
!!$    REAL, INTENT(IN) :: nu1, nu2
!!$    REAL, INTENT(IN) :: rv1, rv2
!!$    REAL :: A0, A1, k1, k2, k3, b0
!!$
!!$    !A=1.+1./(nu1*nu2)
!!$    !k0=1.
!!$
!!$    !B_NL=A*(k/k0)*(1.-rv1*rv2*k**2)
!!$
!!$    A0=0.63
!!$    A1=1.91
!!$    k1=1.23
!!$    k2=4.02
!!$    k3=99.99
!!$    b0=0.0
!!$
!!$    B_NL=(A0+A1/(nu1*nu2))*((k/k1)**1+(k/k2)**2+(k/k3)**4)*(1.-rv1*rv2*k**2)+b0
!!$
!!$    IF(B_NL<-1.) B_NL=-1.
!!$    
!!$  END FUNCTION B_NL

  !REAL FUNCTION BNL(k,nu1,nu2,rv1,rv2,hmod)
  REAL FUNCTION BNL(k,nu1,nu2,hmod)

    IMPLICIT NONE
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: nu1
    REAL, INTENT(IN) :: nu2
    !REAL, INTENT(IN) :: rv1
    !REAL, INTENT(IN) :: rv2
    TYPE(halomod), INTENT(INOUT) :: hmod
    INTEGER  :: nk, nnu
    REAL :: nuu1, nuu2, kk
    INTEGER, PARAMETER :: iorder=1 ! 1 - Linear interpolation
    INTEGER, PARAMETER :: ifind=3  ! 3 - Midpoint finding scheme
    INTEGER, PARAMETER :: imeth=1  ! 1 - Polynomial method    
    REAL, PARAMETER :: kmin=3e-2   ! Below this wavenumber set BNL to zero
    REAL, PARAMETER :: numin=0.    ! Below this halo mass set  BNL to zero
    REAL, PARAMETER :: numax=10.   ! Above this halo mass set  BNL to zero
    REAL, PARAMETER :: BNL_min=-1. ! Minimum value that BNL is allowed to be (could be below -1)

    IF(.NOT. hmod%has_bnl) CALL init_BNL(hmod)

    ! Ensure that nu1 is not outside array boundary
    nuu1=nu1
    CALL fix_min(nuu1,hmod%nu_bnl(1))
    CALL fix_max(nuu1,hmod%nu_bnl(hmod%nnu_bnl))

    ! Ensure that nu2 is not outside array boundary
    nuu2=nu2
    CALL fix_min(nuu2,hmod%nu_bnl(1))
    CALL fix_max(nuu2,hmod%nu_bnl(hmod%nnu_bnl))

    ! Ensure that k is not outside array boundary at high end
    kk=k
    CALL fix_max(kk,hmod%k_bnl(hmod%nk_bnl))

    IF(kk<kmin) THEN
       BNL=0.
    ELSE IF(nu1<numin .OR. nu2<numin) THEN
       BNL=0.
    ELSE IF(nu1>numax .OR. nu2>numax) THEN
       BNL=0.  
    ELSE
       nk=hmod%nk_bnl
       nnu=hmod%nnu_bnl
       BNL=find(log(kk),log(hmod%k_bnl),nuu1,hmod%nu_bnl,nuu2,hmod%nu_bnl,hmod%bnl,nk,nnu,nnu,iorder,ifind,imeth)
    END IF

    !IF(BNL<BNL_min) BNL=BNL_min
    CALL fix_min(BNL,BNL_min)
    
  END FUNCTION BNL

  SUBROUTINE init_BNL(hmod)

    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    INTEGER :: i, ibin, jbin, ik
    INTEGER :: nbin, nk
    REAL :: crap
    CHARACTER(len=256) :: infile, inbase, fbase, fmid, fext
    !CHARACTER(len=256) :: base='/Users/Mead/Physics/data/Bolshoi/power/M512'
    CHARACTER(len=256) :: base='/Users/Mead/Physics/data/Multidark/power/M512'
    REAL, PARAMETER :: eps=1e-2

    WRITE(*,*) 'INIT_BNL: Running'

    ! Read in the nu values from the binstats file
    IF(requal(hmod%z,0.00,eps)) THEN
       inbase=trim(base)//'/BNL_1.00109'
    ELSE IF(requal(hmod%z,0.47,eps)) THEN
       inbase=trim(base)//'/BNL_0.68215'
    ELSE IF(requal(hmod%z,0.69,eps)) THEN
       inbase=trim(base)//'/BNL_0.59103'
    ELSE IF(requal(hmod%z,1.00,eps)) THEN
       inbase=trim(base)//'/BNL_0.49990'
    ELSE IF(requal(hmod%z,2.89,eps)) THEN
       inbase=trim(base)//'/BNL_0.25690'
    ELSE
       STOP 'INIT_BNL: Error, your redshift is not supported'
    END IF

    infile=trim(inbase)//'_binstats.dat'
    nbin=file_length(infile)
    WRITE(*,*) 'INIT_BNL: Number of nu bins: ', nbin
    hmod%nnu_bnl=nbin
    ALLOCATE(hmod%nu_bnl(nbin))
    WRITE(*,*) 'INIT_BNL: Reading: ', trim(infile)
    OPEN(7,file=infile)
    DO i=1,nbin
       READ(7,*) crap, crap, crap, crap, crap, hmod%nu_bnl(i)
       WRITE(*,*) 'INIT_BNL: nu bin', i,  hmod%nu_bnl(i)
    END DO
    CLOSE(7)
    WRITE(*,*) 'INIT_BNL: Done with nu'

    ! Read in k and Bnl(k,nu1,nu2)
    !infile=trim(base)//'/BNL_1.00109_bin1_bin1_power.dat'
    infile=trim(inbase)//'_bin1_bin1_power.dat'
    nk=file_length(infile)
    hmod%nk_bnl=nk
    WRITE(*,*) 'INIT_BNL: Number of k values: ', nk
    ALLOCATE(hmod%k_bnl(nk),hmod%bnl(nk,nbin,nbin))
    DO ibin=1,nbin
       DO jbin=1,nbin
          !fbase=trim(base)//'/BNL_1.00109_bin'
          fbase=trim(inbase)//'_bin'
          fmid='_bin'
          fext='_power.dat'
          infile=number_file2(fbase,ibin,fmid,jbin,fext)
          WRITE(*,*) 'INIT_BNL: Reading: ', trim(infile)
          OPEN(7,file=infile)
          DO ik=1,nk
             READ(7,*) hmod%k_bnl(ik), crap, crap, crap, crap, hmod%bnl(ik,ibin,jbin)
          END DO
          CLOSE(7)
       END DO
    END DO

    ! Convert from Y to B_NL
    hmod%bnl=hmod%bnl-1.

    hmod%has_bnl=.TRUE.

    WRITE(*,*) 'INIT_BNL: Done'
    WRITE(*,*)
    
  END SUBROUTINE init_BNL

  REAL FUNCTION T_1h(k1,k2,ih,hmod,cosm)

    ! Halo model one-halo trispectrum
    IMPLICIT NONE
    REAL, INTENT(IN) :: k1, k2
    INTEGER, INTENT(IN) :: ih(2)
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm   
    REAL :: wk1(hmod%n,2), wk2(hmod%n,2)
    REAL :: uk1(hmod%n,2), uk2(hmod%n,2), uk_quad(hmod%n), uk0_quad
    REAL :: rhom, g, m, m0, integrand(hmod%n)
    INTEGER :: i
    INTEGER, PARAMETER :: iorder=1 ! Use basic trapezium rule because the integrand has rapid oscillations in W(k)

    ! Matter density
    rhom=comoving_matter_density(cosm)

    ! Get the window functions for the two different k values
    CALL init_windows(k1,ih,2,wk1,hmod%n,hmod,cosm)
    CALL init_windows(k2,ih,2,wk2,hmod%n,hmod,cosm)

    ! Create normalised wk
    DO i=1,hmod%n
       m=hmod%m(i)
       uk1(i,:)=wk1(i,:)*rhom/m
       uk2(i,:)=wk2(i,:)*rhom/m
    END DO
  
    ! Not sure if this is the correct k and field combinations
    uk_quad=uk1(:,1)*uk1(:,2)*uk2(:,1)*uk2(:,2)
    !uk_quad=uk1(:,1)*uk1(:,1)*uk2(:,2)*uk2(:,2) ! Could be this
    !uk_quad=uk1(:,1)*uk2(:,2)

    IF(hmod%imf==4) THEN

       ! In this case the mass function is a delta function...
       m0=hmod%hmass
       uk0_quad=find(log(m0),hmod%log_m,uk_quad,hmod%n,3,3,2)
       T_1h=rhom*uk0_quad/m0

    ELSE

       ! ...otherwise you need to do an integral

       ! Calculates the value of the integrand at all nu values!
       DO i=1,hmod%n
          g=g_nu(hmod%nu(i),hmod)
          m=hmod%m(i)
          integrand(i)=g*uk_quad(i)*rhom/m
       END DO

       ! Carries out the integration 
       T_1h=integrate_table(hmod%nu,integrand,hmod%n,1,hmod%n,iorder)

    END IF
    
  END FUNCTION T_1h

  SUBROUTINE halo_diagnostics(hmod,cosm,dir)

    ! Writes out to file a whole set of halo diagnostics
    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm    
    CHARACTER(len=*), INTENT(IN) :: dir
    REAL :: mass    
    CHARACTER(len=64) :: ext
    CHARACTER(len=512) :: base
    CHARACTER(len=1024) :: outfile
    INTEGER :: m, mi1, mi2

    ! Integer 10^m to produce haloes between
    REAL, PARAMETER :: m1=1e10
    REAL, PARAMETER :: m2=1e16

    WRITE(*,*) 'HALO_DIAGNOSTICS: Outputting diagnostics'

    outfile=TRIM(dir)//'/mass_fractions.dat'
    WRITE(*,*) 'HALO_DIAGNOSTICS: ', TRIM(outfile)
    CALL write_mass_fractions(hmod,cosm,outfile)

    IF(hmod%z==0.0) THEN
       ext='_z0.0.dat'
    ELSE IF(hmod%z==0.5) THEN
       ext='_z0.5.dat'
    ELSE IF(hmod%z==1.0) THEN
       ext='_z1.0.dat'
    ELSE IF(hmod%z==2.0) THEN
       ext='_z2.0.dat'
    ELSE
       STOP 'HALO_DIAGNOSTICS: Error, need to make this better with z'
    END IF

    mi1=NINT(log10(m1))
    mi2=NINT(log10(m2))

    DO m=mi1,mi2

       mass=10.**m

       base=TRIM(dir)//'/halo_profile_m'
       outfile=number_file(base,m,ext)
       WRITE(*,*) 'HALO_DIAGNOSTICS: ', TRIM(outfile)
       CALL write_halo_profiles(mass,hmod,cosm,outfile)

       base=TRIM(dir)//'/halo_window_m'
       outfile=number_file(base,m,ext)
       WRITE(*,*) 'HALO_DIAGNOSTICS: ', TRIM(outfile)
       CALL write_halo_transforms(mass,hmod,cosm,outfile)

    END DO

    WRITE(*,*) 'HALO_DIAGNOSTICS: Done'
    WRITE(*,*)

  END SUBROUTINE halo_diagnostics

  SUBROUTINE halo_definitions(hmod,cosm,dir)

    ! Writes out to files the different halo definitions
    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    CHARACTER(len=*), INTENT(IN) :: dir
    CHARACTER(len=256) :: fradius, fmass, fconc
    CHARACTER(len=64) :: ext
    INTEGER :: i

    WRITE(*,*) 'HALO_DEFINITIONS: Outputting definitions'

    IF(hmod%z==0.0) THEN
       ext='_z0.0.dat'
    ELSE IF(hmod%z==0.5) THEN
       ext='_z0.5.dat'
    ELSE IF(hmod%z==1.0) THEN
       ext='_z1.0.dat'
    ELSE IF(hmod%z==2.0) THEN
       ext='_z2.0.dat'
    ELSE
       STOP 'HALO_DEFINITIONS: Error, need to make this better with z'
    END IF

    fradius=TRIM(dir)//'/radius'//TRIM(ext)
    fmass=TRIM(dir)//'/mass'//TRIM(ext)
    fconc=TRIM(dir)//'/concentration'//TRIM(ext)

    WRITE(*,*) 'HALO_DEFINITIONS: ', TRIM(fradius)
    WRITE(*,*) 'HALO_DEFINITIONS: ', TRIM(fmass)
    WRITE(*,*) 'HALO_DEFINITIONS: ', TRIM(fconc)

    IF(hmod%has_mass_conversions .EQV. .FALSE.) CALL convert_mass_definitions(hmod,cosm)

    OPEN(7,file=fradius)
    OPEN(8,file=fmass)
    OPEN(9,file=fconc)
    DO i=1,hmod%n
       WRITE(7,*) hmod%rv(i), hmod%r200(i), hmod%r500(i), hmod%r200c(i), hmod%r500c(i)
       WRITE(8,*) hmod%m(i),  hmod%m200(i), hmod%m500(i), hmod%m200c(i), hmod%m500c(i)
       WRITE(9,*) hmod%c(i),  hmod%c200(i), hmod%c500(i), hmod%c200c(i), hmod%c500c(i)
    END DO
    CLOSE(7)
    CLOSE(8)
    CLOSE(9)

    WRITE(*,*) 'HALO_DEFINITIONS: Done'
    WRITE(*,*)

  END SUBROUTINE halo_definitions

  SUBROUTINE halo_properties(hmod,dir)

    !Writes out to files the different halo definitions
    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    CHARACTER(len=*), INTENT(IN) :: dir
    CHARACTER(len=256) :: output
    CHARACTER(len=64) :: ext
    INTEGER :: i

    WRITE(*,*) 'HALO_PROPERTIES: Outputting definitions'

    IF(hmod%z==0.0) THEN
       ext='_z0.0.dat'
    ELSE IF(hmod%z==0.5) THEN
       ext='_z0.5.dat'
    ELSE IF(hmod%z==1.0) THEN
       ext='_z1.0.dat'
    ELSE IF(hmod%z==2.0) THEN
       ext='_z2.0.dat'
    ELSE
       STOP 'HALO_PROPERTIES: Error, need to make this better with z'
    END IF

    output=TRIM(dir)//'/properies'//TRIM(ext)
    WRITE(*,*) 'HALO_PROPERTIES: ', TRIM(output)

    OPEN(7,file=output)
    DO i=1,hmod%n
       WRITE(7,*) hmod%m(i), hmod%rr(i), hmod%rv(i), hmod%nu(i), hmod%c(i), hmod%sig(i), hmod%sigf(i), hmod%zc(i)
    END DO
    CLOSE(7)

    WRITE(*,*) 'HALO_PROPERTIES: Done'
    WRITE(*,*)

  END SUBROUTINE halo_properties

  SUBROUTINE write_mass_fractions(hmod,cosm,outfile)

    ! Writes out the halo mass fractions
    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    CHARACTER(len=*), INTENT(IN) :: outfile   
    REAL :: m
    INTEGER :: i, j

    REAL, PARAMETER :: mmin=1e10
    REAL, PARAMETER :: mmax=1e16
    INTEGER, PARAMETER :: n=101

    OPEN(7,file=outfile)
    DO i=1,n
       m=exp(progression(log(mmin),log(mmax),i,n))
       WRITE(7,*) m, (halo_fraction(j,m,hmod,cosm), j=1,5)
    END DO
    CLOSE(7)

  END SUBROUTINE write_mass_fractions

  SUBROUTINE write_halo_profiles(m,hmod,cosm,outfile)

    ! Writes out the halo density profiles
    IMPLICIT NONE
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    CHARACTER(len=*), INTENT(IN) :: outfile    
    REAL :: r, rv, rs, c
    INTEGER :: i, j, nf
    INTEGER, ALLOCATABLE :: fields(:)

    REAL, PARAMETER :: rmin=1e-3     ! Mininum r/rv
    REAL, PARAMETER :: rmax=1.1e0    ! Maximum r/rv
    INTEGER, PARAMETER :: n=512      ! Number of points
    LOGICAL, PARAMETER :: real_space=.TRUE. ! Real profiles

    ! Calculate halo attributes
    rv=exp(find(log(m),hmod%log_m,log(hmod%rv),hmod%n,3,3,2))
    c=find(log(m),hmod%log_m,hmod%c,hmod%n,3,3,2)
    rs=rv/c

    ! Field types
    nf=5
    ALLOCATE(fields(nf))
    fields(1)=field_matter
    fields(2)=field_cdm
    fields(3)=field_gas
    fields(4)=field_star
    fields(5)=field_electron_pressure

    ! Write file
    OPEN(7,file=outfile)
    DO i=1,n
       r=exp(progression(log(rmin),log(rmax),i,n))
       r=r*rv
       WRITE(7,*) r/rv, (win_type(real_space,fields(j),r,m,rv,rs,hmod,cosm)*rv**3, j=1,nf) ! rv**3 here is from r^2 dr in integral
    END DO
    CLOSE(7)

  END SUBROUTINE write_halo_profiles

  SUBROUTINE write_halo_transforms(m,hmod,cosm,outfile)

    ! Writes out to file the Fourier transform of the halo density profiles
    IMPLICIT NONE
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    CHARACTER(len=*), INTENT(IN) :: outfile  
    REAL :: x, rv, c, rs, k, rhobar
    INTEGER :: i, j, nf
    INTEGER, ALLOCATABLE :: fields(:)

    REAL, PARAMETER :: xmin=1e-1      ! Minimum r/rv
    REAL, PARAMETER :: xmax=1e2       ! Maximum r/rv
    INTEGER, PARAMETER :: n=512       ! Number of points
    LOGICAL, PARAMETER :: rsp=.FALSE. ! Fourier profiles

    ! Calculate halo attributes
    rv=exp(find(log(m),hmod%log_m,log(hmod%rv),hmod%n,3,3,2))
    c=find(log(m),hmod%log_m,hmod%c,hmod%n,3,3,2)
    rs=rv/c

    ! Field types
    nf=5
    ALLOCATE(fields(nf))
    fields(1)=field_matter
    fields(2)=field_cdm
    fields(3)=field_gas
    fields(4)=field_star
    fields(5)=field_electron_pressure

    ! Need mean density
    rhobar=comoving_matter_density(cosm)

    ! Write file
    OPEN(7,file=outfile)
    DO i=1,n
       x=exp(progression(log(xmin),log(xmax),i,n))
       k=x/rv
       WRITE(7,*) x, (win_type(rsp,fields(j),k,m,rv,rs,hmod,cosm)*rhobar/m, j=1,nf)
    END DO
    CLOSE(7)

  END SUBROUTINE write_halo_transforms

  REAL FUNCTION delta_c(hmod,cosm)

    ! Linear collapse density
    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: a

    a=hmod%a

    IF(hmod%idc==1) THEN
       ! Fixed value
       delta_c=1.686
    ELSE IF(hmod%idc==2) THEN
       ! From Nakamura & Suto (1997) LCDM fitting function
       delta_c=dc_NakamuraSuto(a,cosm)
    ELSE IF(hmod%idc==3 .OR. hmod%idc==6) THEN
       ! From Mead et al. (2015)
       delta_c=hmod%dc0+hmod%dc1*log(sigma(8.,a,cosm))       
       IF(hmod%idc==3) THEN
          ! Mead et al. (2016) addition of small cosmology and neutrino dependence
          delta_c=delta_c*(1.+hmod%dcnu*cosm%f_nu)
          delta_c=delta_c*(dc_NakamuraSuto(a,cosm)/dc0)          
       END IF
    ELSE IF(hmod%idc==4) THEN
       ! From Mead (2017) fitting function
       delta_c=dc_Mead(a,cosm)
    ELSE IF(hmod%idc==5) THEN
       ! From spheircal-collapse calculation
       delta_c=dc_spherical(a,cosm)
    ELSE
       WRITE(*,*) 'DELTA_C: idc:', hmod%idc
       STOP 'DELTA_C: Error, idc defined incorrectly'
    END IF

  END FUNCTION delta_c

  REAL FUNCTION Delta_v(hmod,cosm)

    ! Virialised overdensity
    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: a

    a=hmod%a

    IF(hmod%iDv==1) THEN
       ! Fixed value; M200
       Delta_v=200.
    ELSE IF(hmod%iDv==2) THEN
       ! From Bryan & Norman (1998; arXiv:astro-ph/9710107) fitting functions
       Delta_v=Dv_BryanNorman(a,cosm)    
    ELSE IF(hmod%iDv==3 .OR. hmod%iDv==8) THEN
       ! From Mead et al. (2015, 2016)
       Delta_v=hmod%Dv0*Omega_m(a,cosm)**hmod%Dv1
       IF(hmod%iDv==3) THEN
          ! Mead et al. (2016) neutrino addition
          Delta_v=Delta_v*(1.+hmod%Dvnu*cosm%f_nu) 
       END IF
    ELSE IF(hmod%iDv==4) THEN
       ! From Mead (2017) fitting function
       Delta_v=Dv_Mead(a,cosm)
    ELSE IF(hmod%iDv==5) THEN
       ! From spheircal-collapse calculation
       Delta_v=Dv_spherical(a,cosm)
    ELSE IF(hmod%iDv==6) THEN
       ! Lagrangian radius
       Delta_v=1.
    ELSE IF(hmod%iDv==7) THEN
       ! M200c
       Delta_v=200./Omega_m(a,cosm)
    ELSE
       STOP 'DELTA_V: Error, iDv defined incorrectly'
    END IF

  END FUNCTION Delta_v

  REAL FUNCTION eta_HMcode(hmod,cosm)

    ! Calculates the eta that comes into the bastardised one-halo term
    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: eta0

    IF(hmod%ieta==1) THEN
       eta_HMcode=0.
    ELSE IF(hmod%ieta==2) THEN
       ! From Mead et al. (2015; arXiv 1505.07833, 2016)
       IF(hmod%one_parameter_baryons) THEN
          eta0=0.98-hmod%As*0.12
       ELSE
          eta0=hmod%eta0
       END IF
       eta_HMcode=eta0-hmod%eta1*(sigma(8.,hmod%a,cosm))
    ELSE
       STOP 'ETA_HMcode: Error, ieta defined incorrectly'
    END IF

  END FUNCTION eta_HMcode

  REAL FUNCTION kstar_HMcode(hmod,cosm)

    ! Calculates the one-halo damping wave number
    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm    
    REAL :: crap

    ! To prevent compile-time warnings
    crap=cosm%A

    IF(hmod%ikstar==1) THEN
       ! Set to zero for the standard Poisson one-halo term
       kstar_HMcode=0.
    ELSE IF(hmod%ikstar==2) THEN
       ! One-halo cut-off wavenumber from Mead et al. (2015, 2016)
       kstar_HMcode=hmod%ks/hmod%sigv
    ELSE
       STOP 'KSTAR_HMcode: Error, ikstar defined incorrectly'
    END IF

  END FUNCTION kstar_HMcode

  REAL FUNCTION A_HMcode(hmod,cosm)

    ! Halo concentration pre-factor
    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: crap

    ! To prevent compile-time warnings
    crap=cosm%A

    IF(hmod%iAs==1) THEN
       ! Set to 4 for the standard Bullock value
       A_HMcode=4.
    ELSE IF(hmod%iAs==2) THEN
       ! This is the 'A' halo-concentration parameter in Mead et al. (2015; arXiv 1505.07833, 2016)
       A_HMcode=hmod%As
    ELSE
       STOP 'A_HMcode: Error, iAs defined incorrectly'
    END IF

    ! Now this is divided by 4 so as to be relative to the Bullock base result
    A_HMcode=A_HMcode/4.

  END FUNCTION A_HMcode

  REAL FUNCTION fdamp_HMcode(hmod,cosm)

    ! Calculates the linear-theory damping factor
    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm    
    REAL :: crap

    ! To prevent compile-time warnings
    crap=cosm%A

    IF(hmod%i2hdamp==1) THEN
       ! Set to 0 for the standard linear theory two halo term
       fdamp_HMcode=0.
    ELSE IF(hmod%i2hdamp==2) THEN
       ! Mead et al. (2015)
       fdamp_HMcode=hmod%f0*hmod%sig8z**hmod%f1
    ELSE IF(hmod%i2hdamp==3) THEN
       ! Mead et al. (2016)
       fdamp_HMcode=hmod%f0*hmod%sigv100**hmod%f1
    ELSE
       STOP 'FDAMP_HMcode: Error, i2hdamp defined incorrectly'
    END IF

    ! Catches extreme values of fdamp that occur for ridiculous cosmologies
    IF(fdamp_HMcode<fdamp_HMcode_min) fdamp_HMcode=0.
    IF(fdamp_HMcode>fdamp_HMcode_max) fdamp_HMcode=fdamp_HMcode_max

  END FUNCTION fdamp_HMcode

  REAL FUNCTION alpha_HMcode(hmod,cosm)

    ! Calculates the alpha to smooth the two- to one-halo transition
    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: crap

    ! To prevent compile-time warnings
    crap=cosm%A

    IF(hmod%itrans==2) THEN
       ! From Mead et al. (2015, 2016)   
       alpha_HMcode=hmod%alp0*hmod%alp1**hmod%neff
    ELSE IF(hmod%itrans==3) THEN
       ! Specially for HMx, exponentiated Mead et al. (2016) result
       alpha_HMcode=(hmod%alp0*hmod%alp1**hmod%neff)**1.5
    ELSE IF(hmod%itrans==4) THEN
       ! Specially for HMx, exponentiated Mead et al. (2016) result
       alpha_HMcode=(hmod%alp0*hmod%alp1**hmod%neff)**2.5
    ELSE
       alpha_HMcode=1.
    END IF

    ! Catches values of alpha that are crazy
    IF(alpha_HMcode<alpha_HMcode_min) alpha_HMcode=alpha_HMcode_min
    IF(alpha_HMcode>alpha_HMcode_max) alpha_HMcode=alpha_HMcode_max 

  END FUNCTION alpha_HMcode

  REAL FUNCTION pivot_mass(hmod)

    ! Returns the 'pivot mass' [Msun/h]
    IMPLICIT NONE
    TYPE(halomod), INTENT(IN) :: hmod
    
    IF(hmod%simple_pivot) THEN
       pivot_mass=1e14
    ELSE
       pivot_mass=hmod%Mh
    END IF
    
  END FUNCTION pivot_mass

  REAL FUNCTION HMx_alpha(m,hmod)

    IMPLICIT NONE
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: z, T, A, B, C, D, E, Mp

    IF(hmod%HMx_mode==1) THEN

       HMx_alpha=hmod%alpha

    ELSE IF(hmod%HMx_mode==2) THEN

       Mp=pivot_mass(hmod)
       HMx_alpha=hmod%alpha*((m/Mp)**hmod%alphap)

    ELSE IF(hmod%HMx_mode==3) THEN

       Mp=pivot_mass(hmod)
       z=hmod%z
       HMx_alpha=hmod%alpha*((m/Mp)**hmod%alphap)*((1.+z)**hmod%alphaz)

    ELSE IF(hmod%HMx_mode==4) THEN

       A=hmod%A_alpha
       B=hmod%B_alpha
       C=hmod%C_alpha
       D=hmod%D_alpha
       E=hmod%E_alpha

       z=hmod%z
       T=log10(hmod%Theat)
       !HMx_alpha=A*(1.+z)*T+B*(1.+z)+C*T+D
       HMx_alpha=(A*(1.+z)**2+B*(1.+z)+C)*T**2+D*T+E
       HMx_alpha=HMx_alpha/(1.+z) ! NEW: (1+z) accounts for previous wrong definition of temperature

    ELSE

       STOP 'HMx_ALPHA: Error, HMx_mode not specified correctly'

    END IF

    IF(HMx_alpha<HMx_alpha_min) HMx_alpha=HMx_alpha_min

  END FUNCTION HMx_alpha

  REAL FUNCTION HMx_beta(m,hmod)

    IMPLICIT NONE
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: z, Mp

    IF(hmod%HMx_mode==1) THEN

       HMx_beta=hmod%beta

    ELSE IF(hmod%HMx_mode==2) THEN

       Mp=pivot_mass(hmod)
       HMx_beta=hmod%beta*((m/Mp)**hmod%betap)

    ELSE IF(hmod%HMx_mode==3) THEN

       Mp=pivot_mass(hmod)
       z=hmod%z
       HMx_beta=hmod%beta*((m/Mp)**hmod%betap)*((1.+z)**hmod%betaz)

    ELSE IF(hmod%HMx_mode==4) THEN

       HMx_beta=HMx_alpha(m,hmod)

    ELSE

       STOP 'HMx_BETA: Error, HMx_mode not specified correctly'

    END IF

    IF(HMx_beta<HMx_beta_min) HMx_beta=HMx_beta_min

  END FUNCTION HMx_beta

  REAL FUNCTION HMx_eps(hmod)

    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: z, T, A, B, C, D

    IF(hmod%HMx_mode==1 .OR. hmod%HMx_mode==2) THEN

       HMx_eps=hmod%eps

    ELSE IF(hmod%HMx_mode==3) THEN

       z=hmod%z
       HMx_eps=hmod%eps*((1.+z)**hmod%epsz)

    ELSE IF(hmod%HMx_mode==4) THEN

       A=hmod%A_eps
       B=hmod%B_eps
       C=hmod%C_eps
       D=hmod%D_eps

       z=hmod%z
       T=log10(hmod%Theat)
       HMx_eps=A*(1.+z)*T+B*(1.+z)+C*T+D
       HMx_eps=10**HMx_eps

    ELSE

       STOP 'HMx_EPS: Error, HMx_mode not specified correctly'

    END IF

  END FUNCTION HMx_eps

  REAL FUNCTION HMx_Gamma(m,hmod)

    IMPLICIT NONE
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: z, T, A, B, C, D, E, Mp

    IF(hmod%HMx_mode==1) THEN

       HMx_Gamma=hmod%Gamma

    ELSE IF(hmod%HMx_mode==2) THEN

       Mp=pivot_mass(hmod)
       HMx_Gamma=hmod%Gamma*((m/Mp)**hmod%Gammap)

    ELSE IF(hmod%HMx_mode==3) THEN

       Mp=pivot_mass(hmod)
       z=hmod%z
       HMx_Gamma=hmod%Gamma*((m/Mp)**hmod%Gammap)*((1.+z)**hmod%Gammaz)

    ELSE IF(hmod%HMx_mode==4) THEN

       A=hmod%A_Gamma
       B=hmod%B_Gamma
       C=hmod%C_Gamma
       D=hmod%D_Gamma
       E=hmod%E_Gamma

       z=hmod%z
       T=log10(hmod%Theat)
       !HMx_Gamma=A*(1.+z)*T+B*(1.+z)+C*T+D
       HMx_Gamma=(A*(1.+z)**2+B*(1.+z)+C)*T**2+D*T+E

    ELSE

       STOP 'HMx_GAMMA: Error, HMx_mode not specified correctly'

    END IF

    IF(HMx_Gamma<HMx_Gamma_min) HMx_Gamma=HMx_Gamma_min
    IF(HMx_Gamma>HMx_Gamma_max) HMx_Gamma=HMx_Gamma_max

  END FUNCTION HMx_Gamma

  REAL FUNCTION HMx_iso_beta(m,hmod)

    IMPLICIT NONE
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: z, Mp

    IF(hmod%HMx_mode==1) THEN

       HMx_iso_beta=hmod%ibeta

    ELSE IF(hmod%HMx_mode==2) THEN

       Mp=pivot_mass(hmod)
       HMx_iso_beta=hmod%ibeta*((m/Mp)**hmod%ibetap)

    ELSE IF(hmod%HMx_mode==3) THEN

       Mp=pivot_mass(hmod)
       z=hmod%z
       HMx_iso_beta=hmod%ibeta*((m/Mp)**hmod%ibetap)*((1.+z)**hmod%ibetaz)

    ELSE

       STOP 'HMx_ISO_BETA: Error, HMx_mode not specified correctly'

    END IF

  END FUNCTION HMx_iso_beta

  REAL FUNCTION HMx_M0(hmod)

    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: z, T, A, B, C, D, E

    IF(hmod%HMx_mode==1 .OR. hmod%HMx_mode==2) THEN

       HMx_M0=hmod%M0

    ELSE IF(hmod%HMx_mode==3) THEN

       ! Note, exponential here for z dependence because scaling applies to exponent
       z=hmod%z
       HMx_M0=hmod%M0**((1.+z)**hmod%M0z)

    ELSE IF(hmod%HMx_mode==4) THEN

       A=hmod%A_M0
       B=hmod%B_M0
       C=hmod%C_M0
       D=hmod%D_M0
       E=hmod%E_M0

       z=hmod%z
       T=log10(hmod%Theat)
       !HMx_M0=A*(1.+z)*T+B*(1.+z)+C*T+D
       HMx_M0=(A*(1.+z)**2+B*(1.+z)+C)*T**2+D*T+E
       HMx_M0=10**HMx_M0

    ELSE

       STOP 'HMx_M0: Error, HMx_mode not specified correctly'

    END IF

  END FUNCTION HMx_M0

  REAL FUNCTION HMx_Astar(m,hmod)

    IMPLICIT NONE
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: z, T, A, B, C, D, Mp

    IF(hmod%HMx_mode==1) THEN

       HMx_Astar=hmod%Astar

    ELSE IF(hmod%HMx_mode==2) THEN

       Mp=pivot_mass(hmod)
       HMx_Astar=hmod%Astar*((m/Mp)**hmod%Astarp)

    ELSE IF(hmod%HMx_mode==3) THEN

       Mp=pivot_mass(hmod)
       z=hmod%z
       HMx_Astar=hmod%Astar*((m/Mp)**hmod%Astarp)*((1.+z)**hmod%Astarz)

    ELSE IF(hmod%HMx_mode==4) THEN

       A=hmod%A_Astar
       B=hmod%B_Astar
       C=hmod%C_Astar
       D=hmod%D_Astar

       z=hmod%z
       T=log10(hmod%Theat)
       HMx_Astar=A*(1.+z)*T+B*(1.+z)+C*T+D

    ELSE

       STOP 'HMx_ASTAR: Error, HMx_mode not specified correctly'

    END IF

    IF(HMx_Astar<HMx_Astar_min) HMx_Astar=HMx_Astar_min

  END FUNCTION HMx_Astar

  REAL FUNCTION HMx_Twhim(hmod)

    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: z, T, A, B, C, D

    IF(hmod%HMx_mode==1 .OR. hmod%HMx_mode==2) THEN

       HMx_Twhim=hmod%Twhim

    ELSE IF(hmod%HMx_mode==3) THEN

       ! Note, exponential here for z dependence because scaling applies to exponent
       z=hmod%z
       HMx_Twhim=hmod%Twhim**((1.+z)**hmod%Twhimz)

    ELSE IF(hmod%HMx_mode==4) THEN

       A=hmod%A_Twhim
       B=hmod%B_Twhim
       C=hmod%C_Twhim
       D=hmod%D_Twhim

       z=hmod%z
       T=log10(hmod%Theat)
       !HMx_Twhim=A*(1.+z)*T+B*(1.+z)+C*T+D
       HMx_Twhim=(A*(1.+z)**2+B*(1.+z)+C)*T+D
       HMx_Twhim=10**HMx_Twhim

    ELSE

       STOP 'HMx_TWHIM: Error, HMx_mode not specified correctly'

    END IF

  END FUNCTION HMx_Twhim

  REAL FUNCTION HMx_cstar(m,hmod)

    IMPLICIT NONE
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: Mp, z

    IF(hmod%HMx_mode==1) THEN

       HMx_cstar=hmod%cstar

    ELSE IF(hmod%HMx_mode==2) THEN

       Mp=pivot_mass(hmod)
       HMx_cstar=hmod%cstar*((m/Mp)**hmod%cstarp)

    ELSE IF(hmod%HMx_mode==3) THEN

       Mp=pivot_mass(hmod)
       z=hmod%z
       HMx_cstar=hmod%cstar*((m/Mp)**hmod%cstarp)*((1.+z)**hmod%cstarz)

    ELSE IF(hmod%HMx_mode==4) THEN

       HMx_cstar=hmod%cstar

    ELSE

       STOP 'HMx_CSTAR: Error, HMx_mode not specified correctly'

    END IF

  END FUNCTION HMx_cstar

  REAL FUNCTION HMx_Mstar(hmod)

    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: z

    IF(hmod%HMx_mode==1 .OR. hmod%HMx_mode==2) THEN

       HMx_Mstar=hmod%mstar

    ELSE IF(hmod%HMx_mode==3) THEN

       ! Note, exponential here for z dependence because scaling applies to exponent
       z=hmod%z
       HMx_Mstar=hmod%mstar**((1.+z)**hmod%mstarz)

    ELSE IF(hmod%HMx_mode==4) THEN

       HMx_Mstar=hmod%mstar

    ELSE

       STOP 'HMx_MSTAR: Error, HMx_mode not specified correctly'

    END IF

  END FUNCTION HMx_Mstar

  FUNCTION r_nl(hmod)

    ! Calculates R_nl where nu(R_nl)=1.
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: r_nl  

    IF(hmod%nu(1)>1.) THEN
       ! This catches some very strange values
       r_nl=hmod%rr(1)
    ELSE
       r_nl=exp(find(log(1.),log(hmod%nu),log(hmod%rr),hmod%n,3,3,2))
    END IF

  END FUNCTION r_nl

  SUBROUTINE allocate_HMOD(hmod)

    ! Allocates memory for the look-up tables
    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    INTEGER :: n

    ! Number of entries in look-up table
    n=hmod%n

    ALLOCATE(hmod%log_m(n))
    ALLOCATE(hmod%zc(n),hmod%m(n),hmod%c(n),hmod%rv(n))
    ALLOCATE(hmod%nu(n),hmod%rr(n),hmod%sigf(n),hmod%sig(n))
    ALLOCATE(hmod%m500(n),hmod%r500(n),hmod%c500(n))
    ALLOCATE(hmod%m500c(n),hmod%r500c(n),hmod%c500c(n))
    ALLOCATE(hmod%m200(n),hmod%r200(n),hmod%c200(n))
    ALLOCATE(hmod%m200c(n),hmod%r200c(n),hmod%c200c(n))

    hmod%log_m=0.
    hmod%zc=0.
    hmod%m=0.
    hmod%c=0.
    hmod%rv=0.
    hmod%nu=0.
    hmod%rr=0.
    hmod%sigf=0.
    hmod%sig=0.

    hmod%m500=0.
    hmod%r500=0.
    hmod%c500=0.

    hmod%m500c=0.
    hmod%r500c=0.
    hmod%c500c=0.

    hmod%m200=0.
    hmod%r200=0.
    hmod%c200=0.

    hmod%m200c=0.
    hmod%r200c=0.
    hmod%c200c=0.

  END SUBROUTINE allocate_HMOD

  SUBROUTINE deallocate_HMOD(hmod)

    !Deallocates the look-up tables
    IMPLICIT NONE
    TYPE(halomod) :: hmod

    !Deallocates look-up tables
    DEALLOCATE(hmod%log_m)
    DEALLOCATE(hmod%zc,hmod%m,hmod%c,hmod%rv,hmod%nu,hmod%rr,hmod%sigf,hmod%sig)
    DEALLOCATE(hmod%m500,hmod%r500,hmod%c500,hmod%m500c,hmod%r500c,hmod%c500c)
    DEALLOCATE(hmod%m200,hmod%r200,hmod%c200,hmod%m200c,hmod%r200c,hmod%c200c)

  END SUBROUTINE deallocate_HMOD

  REAL FUNCTION Omega_stars(hmod,cosm)

    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm

    Omega_stars=rhobar(hmod%nu(1),hmod%large_nu,rhobar_star_integrand,hmod,cosm)
    Omega_stars=Omega_stars/comoving_critical_density(hmod%a,cosm)

  END FUNCTION Omega_stars

  SUBROUTINE init_galaxies(hmod,cosm)

    !Calcuate the number densities of galaxies
    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: nu_min, nu_max

    nu_min=nu_M(hmod%mhalo_min,hmod,cosm)
    nu_max=nu_M(hmod%mhalo_max,hmod,cosm)
    hmod%n_c=rhobar(nu_min,nu_max,rhobar_central_integrand,hmod,cosm)
    hmod%n_s=rhobar(nu_min,nu_max,rhobar_satellite_integrand,hmod,cosm)
    hmod%n_g=hmod%n_c+hmod%n_s
    IF(verbose_galaxies) THEN
       WRITE(*,*) 'INIT_GALAXIES: Comoving density of central galaxies [(Mpc/h)^-3]:', REAL(hmod%n_c)
       WRITE(*,*) 'INIT_GALAXIES: Comoving density of satellite galaxies [(Mpc/h)^-3]:', REAL(hmod%n_s)
       WRITE(*,*) 'INIT_GALAXIES: Comoving density of all galaxies [(Mpc/h)^-3]:', REAL(hmod%n_g)
       WRITE(*,*)
    END IF

    hmod%has_galaxies=.TRUE.

  END SUBROUTINE init_galaxies

  SUBROUTINE init_HI(hmod,cosm)

    ! Calculates the background HI density by integrating the HI mass function
    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: nu_min, nu_max

    nu_min=hmod%nu(1)
    nu_max=hmod%nu(hmod%n)
    hmod%rho_HI=rhobar(nu_min,hmod%large_nu,rhobar_HI_integrand,hmod,cosm)
    IF(verbose_HI) THEN
       WRITE(*,*) 'INIT_HI: z:', hmod%z
       WRITE(*,*) 'INIT_HI: HI density [log10(rho/(Msun/h)/(Mpc/h)^3)]:', REAL(log10(hmod%rho_HI))
       WRITE(*,*) 'INIT_HI: Omega_HI(z):', hmod%rho_HI/comoving_critical_density(hmod%a,cosm)
       WRITE(*,*) 'INIT_HI: Omega_HI(z) relative to z=0 critical density:', hmod%rho_HI/comoving_critical_density(1.,cosm)
       WRITE(*,*) 'INIT_HI: rho_HI / rho_matter:', hmod%rho_HI/comoving_matter_density(cosm)
       WRITE(*,*)
    END IF

    hmod%has_HI=.TRUE.

  END SUBROUTINE init_HI

  REAL FUNCTION nu_R(R,hmod,cosm)

    ! Calculates nu(R) where R is the comoving Lagrangian halo radius
    IMPLICIT NONE
    REAL, INTENT(IN) :: R
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: a

    a=hmod%a
    nu_R=delta_c(hmod,cosm)/sigma(R,a,cosm)

  END FUNCTION nu_R

  REAL FUNCTION nu_M(M,hmod,cosm)

    ! Calculates nu(M) where M is the halo mass
    IMPLICIT NONE
    REAL, INTENT(IN) :: M
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: R

    R=radius_m(M,cosm)
    nu_M=nu_R(R,hmod,cosm)

  END FUNCTION nu_M

  REAL FUNCTION M_nu(nu,hmod)

    !Calculates M(nu) where M is the halo mass and nu is the peak height
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu
    TYPE(halomod), INTENT(INOUT) :: hmod

    M_nu=exp(find(nu,hmod%nu,hmod%log_m,hmod%n,3,3,2))

  END FUNCTION M_nu

  REAL FUNCTION rhobar_central_integrand(nu,hmod,cosm)

    !Integrand for the number density of central galaxies
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm ! Could remove
    REAL :: M, crap

    crap=cosm%A

    M=M_nu(nu,hmod)
    rhobar_central_integrand=N_centrals(M,hmod)*g_nu(nu,hmod)/M

  END FUNCTION rhobar_central_integrand

  REAL FUNCTION rhobar_satellite_integrand(nu,hmod,cosm)

    !Integrand for the number density of satellite galaxies
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm ! Could remove
    REAL :: M, crap

    crap=cosm%A

    M=M_nu(nu,hmod)    
    rhobar_satellite_integrand=N_satellites(M,hmod)*g_nu(nu,hmod)/M

  END FUNCTION rhobar_satellite_integrand

  REAL FUNCTION rhobar_star_integrand(nu,hmod,cosm)

    !Integrand for the matter density of stars
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: M

    M=M_nu(nu,hmod)    
    rhobar_star_integrand=halo_star_fraction(M,hmod,cosm)*g_nu(nu,hmod)

  END FUNCTION rhobar_star_integrand

  REAL FUNCTION rhobar_HI_integrand(nu,hmod,cosm)

    !Integrand for the HI mass density
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: M

    M=M_nu(nu,hmod)    
    rhobar_HI_integrand=halo_HI_fraction(M,hmod,cosm)*g_nu(nu,hmod)

  END FUNCTION rhobar_HI_integrand

  REAL FUNCTION rhobar(nu_min,nu_max,integrand,hmod,cosm)

    ! Calculate the mean density of a tracer
    ! Integrand here is a function of mass, i.e. I(M); R = rho * Int I(M)dM
    ! TODO: This function is comparatively slow and should be accelerated somehow
    ! TODO: This uses integrate_hmod_cosm_exp, which is weird, surely can use some transformed integrand instead?
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu_min, nu_max
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm

    INTERFACE
       REAL FUNCTION integrand(M,hmod,cosm)
         IMPORT :: halomod, cosmology
         REAL, INTENT(IN) :: M
         TYPE(halomod), INTENT(INOUT) :: hmod
         TYPE(cosmology), INTENT(INOUT) :: cosm
       END FUNCTION integrand
    END INTERFACE

    rhobar=comoving_matter_density(cosm)*integrate_hmod_cosm_exp(log(nu_min),log(nu_max),integrand,hmod,cosm,hmod%acc_HMx,3)

  END FUNCTION rhobar

  FUNCTION one_halo_amplitude(hmod,cosm)

    !Calculates the amplitude of the shot-noise plateau of the one-halo term
    IMPLICIT NONE
    REAL :: one_halo_amplitude
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: integrand(hmod%n), g, m
    INTEGER :: i

    !Calculates the value of the integrand at all nu values!
    DO i=1,hmod%n
       g=g_nu(hmod%nu(i),hmod)
       m=hmod%m(i)
       integrand(i)=g*m
    END DO

    one_halo_amplitude=integrate_table(hmod%nu,integrand,hmod%n,1,hmod%n,1)/comoving_matter_density(cosm)

  END FUNCTION one_halo_amplitude

  REAL FUNCTION mass_function(m,hmod,cosm)

    ! Calculates the halo mass function, what I call n(M); some people call dn/dM [(Msun/h)^-1(Mpc/h)^-3]
    IMPLICIT NONE
    REAL, INTENT(IN) :: m                  ! Halo mass [Msun/h]
    TYPE(halomod), INTENT(INOUT) :: hmod   ! Halo model
    TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology

    mass_function=multiplicity_function(m,hmod,cosm)*comoving_matter_density(cosm)/m**2

  END FUNCTION mass_function

  REAL FUNCTION multiplicity_function(m,hmod,cosm)

    ! Returns the dimensionless multiplicity function: M^2n(M)/rho; n(M) = dn/dM sometimes
    IMPLICIT NONE
    REAL, INTENT(IN) :: m                  ! Halo mass [Msun/h]
    TYPE(halomod), INTENT(INOUT) :: hmod   ! Halo model
    TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology
    REAL :: nu, dnu_dlnm

    nu=nu_M(m,hmod,cosm)
    dnu_dlnm=derivative_table(log(m),log(hmod%m),hmod%nu,hmod%n,3,3)
    multiplicity_function=g_nu(nu,hmod)*dnu_dlnm
    
  END FUNCTION multiplicity_function

  REAL FUNCTION halo_bias(m,hmod,cosm)

    ! Calculate b(M)
    IMPLICIT NONE
    REAL, INTENT(IN) :: m                  ! Halo mass [Msun/h]
    TYPE(halomod), INTENT(INOUT) :: hmod   ! Halo model
    TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology
    REAL :: nu

    nu=nu_M(m,hmod,cosm)
    halo_bias=b_nu(nu,hmod)

  END FUNCTION halo_bias

  SUBROUTINE convert_mass_definitions(hmod,cosm)

    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: rhom, rhoc, Dv, a

    ! Scale factor
    a=hmod%a

    ! Get the densities
    rhom=comoving_matter_density(cosm)
    rhoc=comoving_critical_density(a,cosm)
    Dv=Delta_v(hmod,cosm)

    ! Calculate Delta = 200, 500 and Delta_c = 200, 500 quantities
    CALL convert_mass_definition(hmod%rv,hmod%c,hmod%m,Dv,1.,hmod%r500,hmod%c500,hmod%m500,500.,1.,hmod%n)
    CALL convert_mass_definition(hmod%rv,hmod%c,hmod%m,Dv,1.,hmod%r200,hmod%c200,hmod%m200,200.,1.,hmod%n)
    CALL convert_mass_definition(hmod%rv,hmod%c,hmod%m,Dv,rhom,hmod%r500c,hmod%c500c,hmod%m500c,500.,rhoc,hmod%n)
    CALL convert_mass_definition(hmod%rv,hmod%c,hmod%m,Dv,rhom,hmod%r200c,hmod%c200c,hmod%m200c,200.,rhoc,hmod%n)

    hmod%has_mass_conversions=.TRUE.

  END SUBROUTINE convert_mass_definitions

  SUBROUTINE convert_mass_definition(r1,c1,m1,D1,rho1,r2,c2,m2,D2,rho2,n)

    !Converts mass definition from Delta_1 rho_1 overdense to Delta_2 rho_2 overdense
    IMPLICIT NONE
    REAL, INTENT(IN) :: r1(n) ! Array of initial virial radii [Mpc/h]
    REAL, INTENT(IN) :: c1(n) ! Array of initial halo concentration
    REAL, INTENT(IN) :: m1(n) ! Array of initial halo mass [Msun/h]
    REAL, INTENT(IN) :: D1 ! Initial halo overdensity definition (e.g., 200, 500)
    REAL, INTENT(IN) :: rho1 ! Initial halo overdensity defintion (critical or mass)
    REAL, INTENT(OUT) :: r2(n) ! Output array of virial radii with new definition [Mpc/h]
    REAL, INTENT(OUT) :: c2(n) ! Output array of halo concentration with new definition
    REAL, INTENT(OUT) :: m2(n) ! Output array of halo mass with new definition [Msun/h]
    REAL, INTENT(IN) :: D2 ! Final halo overdensity definition (e.g., 200, 500)
    REAL, INTENT(IN) :: rho2 ! Final halo overdensity defintion (critical or mass)
    INTEGER, INTENT(IN) :: n  ! Number of entries in tables
    REAL :: f(n)
    REAL :: rmin, rmax, rs
    REAL, ALLOCATABLE :: r(:)
    INTEGER :: i, j

    IF(verbose_convert_mass) THEN
       WRITE(*,*) 'CONVERT_MASS_DEFINITION: Converting mass definitions:'
       WRITE(*,*) 'CONVERT_MASS_DEFINITION: Initial overdensity:', D1
       WRITE(*,*) 'CONVERT_MASS_DEFINITION: Final overdensity:', D2
    END IF

    !Make an array of general 'r' values for the solution later
    rmin=r1(1)/10. !Should be sufficient for reasonable cosmologies
    rmax=r1(n)*10. !Should be sufficient for reasonable cosmologies
    CALL fill_array(log(rmin),log(rmax),r,n) !Necessary to log space
    r=exp(r)

    IF(verbose_convert_mass) THEN
       WRITE(*,*) 'CONVERT_MASS_DEFINITION: rmin:', rmin
       WRITE(*,*) 'CONVERT_MASS_DEFINITION: rmax:', rmax
       WRITE(*,*) 'CONVERT_MASS_DEFINITION: nr:', n
    END IF

    !Now use the find algorithm to invert L(r_i)=R(r_j) so that r_j=R^{-1}[L(r_i)]
    IF(verbose_convert_mass) THEN
       WRITE(*,*) '========================================================================================================'
       WRITE(*,*) '         M_old         rv_old          c_old    M_new/M_old    r_new/r_old    c_new/c_old  M`_new/M`_old' 
       WRITE(*,*) '========================================================================================================'
    END IF
    DO i=1,n

       !Calculate the halo scale radius
       rs=r1(i)/c1(i)

       !Fill up the f(r) table which needs to be solved for R | f(R)=0 
       DO j=1,n
          f(j)=r1(i)**3*rho1*D1/NFW_factor(r1(i)/rs)-r(j)**3*rho2*D2/NFW_factor(r(j)/rs) !NFW
          !f(j)=r1(i)**2*rho1*D1-r(j)**2*rho2*D2 !Isothemral sphere
       END DO

       !First find the radius R | f(R)=0; I am fairly certain that I can use log on 'r' here
       r2(i)=exp(find_solve(0.,log(r),f,n))

       !Now do the concentration and mass conversions
       c2(i)=r2(i)/rs
       m2(i)=m1(i)*(rho2*D2*r2(i)**3)/(rho1*D1*r1(i)**3)
       !m2(i)=m1(i)*NFW_factor(c2(i))/NFW_factor(c1(i))

       IF(verbose_convert_mass) THEN
          WRITE(*,fmt='(ES15.7,6F15.7)') m1(i), r1(i), c1(i), m2(i)/m1(i), r2(i)/r1(i), c2(i)/c1(i), NFW_factor(c2(i))/NFW_factor(c1(i))
       END IF

    END DO

    IF(verbose_convert_mass) WRITE(*,*) '========================================================================================================'

    IF(verbose_convert_mass) THEN
       WRITE(*,*) 'CONVERT_MASS_DEFINITION: Done'
       WRITE(*,*)
       STOP
    END IF

  END SUBROUTINE convert_mass_definition

  REAL FUNCTION NFW_factor(x)

    ! The NFW 'mass' factor that crops up all the time
    ! This is X(c) in M(r) = M X(r/rs) / X(c)
    IMPLICIT NONE
    REAL, INTENT(IN) :: x

    NFW_factor=log(1.+x)-x/(1.+x)

  END FUNCTION NFW_factor

  FUNCTION radius_m(m,cosm)

    ! The comoving radius corresponding to mass M in a homogeneous universe
    IMPLICIT NONE
    REAL :: radius_m
    REAL, INTENT(IN) :: m
    TYPE(cosmology), INTENT(INOUT) :: cosm

    radius_m=(3.*m/(4.*pi*comoving_matter_density(cosm)))**(1./3.)

  END FUNCTION radius_m

  FUNCTION virial_radius(m,hmod,cosm)

    ! The comoving halo virial radius
    IMPLICIT NONE
    REAL :: virial_radius
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm

    virial_radius=(3.*m/(4.*pi*comoving_matter_density(cosm)*Delta_v(hmod,cosm)))**(1./3.)

  END FUNCTION virial_radius

  REAL FUNCTION effective_index(hmod,cosm)

    ! Power spectrum slope a the non-linear scale
    IMPLICIT NONE
    TYPE(cosmology), INTENT(INOUT) :: cosm
    TYPE(halomod), INTENT(INOUT) :: hmod

    ! Numerical differentiation to find effective index at collapse
    effective_index=-3.-derivative_table(log(hmod%rnl),log(hmod%rr),log(hmod%sig**2),hmod%n,3,3)

    ! For some bizarre cosmologies r_nl is very small, so almost no collapse has occured
    ! In this case the n_eff calculation goes mad and needs to be fixed using this fudge.
    IF(effective_index<cosm%n-4.) effective_index=cosm%n-4.
    IF(effective_index>cosm%n)      effective_index=cosm%n

  END FUNCTION effective_index

  SUBROUTINE fill_halo_concentration(hmod,cosm)

    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm    
    REAL :: mnl, m, zc, z
    INTEGER :: i

    ! iconc = 1: Full Bullock et al. (2001)
    ! iconc = 2: Simple Bullock et al. (2001)
    ! iconc = 3: Duffy et al. (2008): mean
    ! iconc = 4: Duffy et al. (2008): virial
    ! iconc = 5: Duffy et al. (2008): relaxed

    ! Get the redshift
    z=hmod%z

    ! Any initialisation for the c(M) relation goes here
    IF(hmod%iconc==1) THEN
       ! Fill the collapse z look-up table
       CALL zcoll_Bullock(z,hmod,cosm)
    ELSE IF(hmod%iconc==2) THEN
       mnl=hmod%mnl
    END IF

    ! Fill concentration-mass for all halo masses
    DO i=1,hmod%n

       ! Halo mass
       m=hmod%m(i)

       ! Choose concentration-mass relation
       IF(hmod%iconc==1) THEN
          zc=hmod%zc(i)
          hmod%c(i)=conc_Bullock(z,zc)
       ELSE IF(hmod%iconc==2) THEN         
          hmod%c(i)=conc_Bullock_simple(m,mnl)
       ELSE IF(hmod%iconc==3) THEN
          hmod%c(i)=conc_Duffy_full_M200(m,z)
       ELSE IF(hmod%iconc==4) THEN
          hmod%c(i)=conc_Duffy_full_virial(m,z)
       ELSE IF(hmod%iconc==5) THEN
          hmod%c(i)=conc_Duffy_full_M200c(m,z)
       ELSE IF(hmod%iconc==6) THEN
          hmod%c(i)=conc_Duffy_relaxed_M200(m,z)
       ELSE IF(hmod%iconc==7) THEN
          hmod%c(i)=conc_Duffy_relaxed_virial(m,z)
       ELSE IF(hmod%iconc==8) THEN
          hmod%c(i)=conc_Duffy_relaxed_M200c(m,z)
       ELSE
          STOP 'FILL_HALO_CONCENTRATION: Error, iconc specified incorrectly'
       END IF

       ! Rescale halo concentrations via the 'A' HMcode parameter
       hmod%c(i)=hmod%c(i)*A_HMcode(hmod,cosm)   

    END DO

    ! Dolag2004 prescription for adding DE dependence
    IF(hmod%iDolag .NE. 1) CALL Dolag_correction(hmod,cosm)

    ! Rescale the concentration-mass relation for gas the epsilon parameter
    ! This only rescales the concentrations of haloes that *contain* substantial amounts of gas
    DO i=1,hmod%n
       m=hmod%m(i)
       hmod%c(i)=hmod%c(i)*gas_correction(m,hmod,cosm)
       !WRITE(*,*) i, m, hmod%c(i), gas_correction(m,hmod,cosm)
    END DO
    !STOP
    
  END SUBROUTINE fill_halo_concentration

  REAL FUNCTION gas_correction(m,hmod,cosm)

    ! BUG: Fixed so that concentration modification applies to low-mass haloes
    IMPLICIT NONE
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: gas_supply ! Fraction of universal baryon ratio remaining in halo

    gas_supply=halo_bound_gas_fraction(m,hmod,cosm)/(cosm%Om_b/cosm%Om_m)
    gas_correction=1.+(HMx_eps(hmod)-1.)*(1.-gas_supply)
    
  END FUNCTION gas_correction

  SUBROUTINE Dolag_correction(hmod,cosm)

    ! Applies the Dolag et al. (2004) correction to concentration-mass relation
    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: g_LCDM, g_wCDM, g, a
    REAL :: ginf_LCDM, ginf_wCDM, ainf, f
    TYPE(cosmology) :: cosm_LCDM

    ! The 'infinite' scale factor
    ainf=scale_factor_z(zinf_Dolag)

    ! Save the growth function in the current cosmology
    ginf_wCDM=grow(ainf,cosm)
    
    ! Make a flat LCDM cosmology and calculate growth
    cosm_LCDM=cosm
    cosm_LCDM%iw=1
    cosm_LCDM%w=-1.
    cosm_LCDM%wa=0.
    cosm_LCDM%Om_w=0.
    cosm_LCDM%Om_v=1.-cosm%Om_m ! Added this so that 'making a LCDM cosmology' works for curved models.
    cosm_LCDM%verbose=.FALSE.
    CALL init_cosmology(cosm_LCDM) ! This is **essential**

    ! Growth factor in LCDM at 'infinity' calculated using Linder approximation
    ginf_LCDM=grow_Linder(ainf,cosm_LCDM)

    ! Fractional difference compared to LCDM
    f=ginf_wCDM/ginf_LCDM
    
    IF(hmod%iDolag==2) THEN
       ! Standard correction (Mead et al. 2015)
       hmod%c=hmod%c*f
    ELSE IF(hmod%iDolag==3) THEN
       ! Changed this to a power of 1.5 in Mead et al. 2016, produces more accurate results for extreme DE
       hmod%c=hmod%c*f**1.5
    ELSE IF(hmod%iDolag==4) THEN
       ! Correction with a sensible redshift dependence
       a=hmod%a
       g_wCDM=grow(a,cosm)
       g_LCDM=grow_Linder(a,cosm_LCDM)
       g=g_wCDM/g_LCDM       
       hmod%c=hmod%c*f/g
    ELSE
       STOP 'DOLAG_CORRECTION: Error, iDolag specified incorrectly'
    END IF

  END SUBROUTINE Dolag_correction

  REAL FUNCTION conc_Bullock(z,zc)

    IMPLICIT NONE
    REAL, INTENT(IN) :: z
    REAL, INTENT(IN) :: zc

    REAL, PARAMETER :: A=4. !Pre-factor for Bullock relation

    conc_Bullock=A*(1.+zc)/(1.+z)

  END FUNCTION conc_Bullock

  SUBROUTINE zcoll_Bullock(z,hmod,cosm)

    !This fills up the halo collapse redshift table as per Bullock relations   
    IMPLICIT NONE
    REAL, INTENT(IN) :: z
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm    
    REAL :: dc
    REAL :: af, zf, RHS, a, growz
    INTEGER :: i

    REAL, PARAMETER :: f=0.01**(1./3.) ! This is the f=0.01 parameter in the Bullock realtion sigma(fM,z)

    a=scale_factor_z(z)

    ! Fills up a table for sigma(fM) for Bullock c(m) relation    
    IF(hmod%iconc==1) THEN
       DO i=1,hmod%n
          hmod%sigf(i)=sigma(hmod%rr(i)*f,a,cosm)
       END DO
    END IF

    ! Do numerical inversion
    DO i=1,hmod%n

       ! I don't think this is really consistent with dc varying as a function of z
       ! but the change will *probably* be *very* small
       dc=delta_c(hmod,cosm)

       RHS=dc*grow(a,cosm)/hmod%sigf(i)

       ! growz=find(a,af_tab,grow_tab,cosm%ng,3,3,2)
       ! growz=exp(find(log(a),cosm%a_growth,cosm%growth,cosm%n_growth,3,3,2))
       growz=grow(a,cosm)

       IF(RHS>growz) THEN
          zf=z
       ELSE
          af=exp(find(log(RHS),cosm%log_growth,cosm%log_a_growth,cosm%n_growth,3,3,2))
          zf=redshift_a(af)
       END IF

       hmod%zc(i)=zf

    END DO

  END SUBROUTINE zcoll_Bullock

  FUNCTION conc_Bullock_simple(m,mstar)

    ! The simple concentration-mass relation from Bullock et al. (2001; astro-ph/9908159v3 equation 18)
    IMPLICIT NONE
    REAL :: conc_Bullock_simple
    REAL, INTENT(IN) :: m, mstar

    conc_Bullock_simple=9.*(m/mstar)**(-0.13)

  END FUNCTION conc_Bullock_simple

  REAL FUNCTION conc_Duffy_full_M200c(m,z)

    ! Duffy et al (2008; 0804.2486) c(M) relation for WMAP5, See Table 1
    IMPLICIT NONE
    REAL, INTENT(IN) :: m, z

    REAL, PARAMETER :: m_piv=2e12 ! Pivot mass in Msun/h
    REAL, PARAMETER :: A=5.71
    REAL, PARAMETER :: B=-0.084
    REAL, PARAMETER :: C=-0.47

    ! Equation (4) in 0804.2486, parameters from 4th row of Table 1
    conc_Duffy_full_M200c=A*(m/m_piv)**B*(1.+z)**C

  END FUNCTION conc_Duffy_full_M200c

  REAL FUNCTION conc_Duffy_relaxed_M200c(m,z)

    ! Duffy et al (2008; 0804.2486) c(M) relation for WMAP5, See Table 1
    IMPLICIT NONE
    REAL, INTENT(IN) :: m, z

    REAL, PARAMETER :: m_piv=2e12 ! Pivot mass in Msun/h
    REAL, PARAMETER :: A=6.71
    REAL, PARAMETER :: B=-0.091
    REAL, PARAMETER :: C=-0.44

    ! Equation (4) in 0804.2486, parameters from 4th row of Table 1
    conc_Duffy_relaxed_M200c=A*(m/m_piv)**B*(1.+z)**C

  END FUNCTION conc_Duffy_relaxed_M200c

  REAL FUNCTION conc_Duffy_full_virial(m,z)

    ! Duffy et al (2008; 0804.2486) c(M) relation for WMAP5, See Table 1
    IMPLICIT NONE
    REAL, INTENT(IN) :: m, z

    REAL, PARAMETER :: m_piv=2e12 ! Pivot mass in Msun/h
    REAL, PARAMETER :: A=7.85
    REAL, PARAMETER :: B=-0.081
    REAL, PARAMETER :: C=-0.71

    ! Equation (4) in 0804.2486, parameters from 6th row of Table 1
    conc_Duffy_full_virial=A*(m/m_piv)**B*(1.+z)**C

  END FUNCTION conc_Duffy_full_virial

  REAL FUNCTION conc_Duffy_relaxed_virial(m,z)

    ! Duffy et al (2008; 0804.2486) c(M) relation for WMAP5, See Table 1
    IMPLICIT NONE
    REAL, INTENT(IN) :: m, z

    REAL, PARAMETER :: m_piv=2e12 ! Pivot mass in Msun/h
    REAL, PARAMETER :: A=9.23
    REAL, PARAMETER :: B=-0.090
    REAL, PARAMETER :: C=-0.69

    ! Equation (4) in 0804.2486, parameters from 6th row of Table 1
    conc_Duffy_relaxed_virial=A*(m/m_piv)**B*(1.+z)**C

  END FUNCTION conc_Duffy_relaxed_virial  

  REAL FUNCTION conc_Duffy_full_M200(m,z)

    ! Duffy et al (2008; 0804.2486) c(M) relation for WMAP5, See Table 1
    IMPLICIT NONE
    REAL, INTENT(IN) :: m, z

    REAL, PARAMETER :: m_piv=2e12 ! Pivot mass in Msun/h
    REAL, PARAMETER :: A=10.14
    REAL, PARAMETER :: B=-0.081
    REAL, PARAMETER :: C=-1.01

    ! Equation (4) in 0804.2486, parameters from 10th row of Table 1
    conc_Duffy_full_M200=A*(m/m_piv)**B*(1.+z)**C

  END FUNCTION conc_Duffy_full_M200

  REAL FUNCTION conc_Duffy_relaxed_M200(m,z)

    ! Duffy et al (2008; 0804.2486) c(M) relation for WMAP5, See Table 1
    IMPLICIT NONE
    REAL, INTENT(IN) :: m, z

    REAL, PARAMETER :: m_piv=2e12 ! Pivot mass in Msun/h
    REAL, PARAMETER :: A=11.93
    REAL, PARAMETER :: B=-0.090
    REAL, PARAMETER :: C=-0.99

    ! Equation (4) in 0804.2486, parameters from 10th row of Table 1
    conc_Duffy_relaxed_M200=A*(m/m_piv)**B*(1.+z)**C

  END FUNCTION conc_Duffy_relaxed_M200

  REAL FUNCTION mass_r(r,cosm)

    ! Calcuates the mass contains in a sphere of comoving radius 'r' in a homogeneous universe
    IMPLICIT NONE
    REAL, INTENT(IN) :: r
    TYPE(cosmology), INTENT(INOUT) :: cosm

    ! Relation between mean cosmological mass and radius
    mass_r=(4.*pi/3.)*comoving_matter_density(cosm)*(r**3)

  END FUNCTION mass_r

  REAL FUNCTION win_type(real_space,ifield,k,m,rv,rs,hmod,cosm)

    ! Selects the halo profile type
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: real_space
    INTEGER, INTENT(IN) :: ifield
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    REAL, INTENT(IN) :: rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: nu, mmin, mmax

    IF(ifield==field_dmonly) THEN
       win_type=win_DMONLY(real_space,k,m,rv,rs,hmod,cosm)
    ELSE IF(ifield==field_matter) THEN
       win_type=win_total(real_space,k,m,rv,rs,hmod,cosm)
    ELSE IF(ifield==field_cdm) THEN
       win_type=win_CDM(real_space,k,m,rv,rs,hmod,cosm)
    ELSE IF(ifield==field_gas) THEN
       win_type=win_gas(real_space,ifield,k,m,rv,rs,hmod,cosm)
    ELSE IF(ifield==field_star) THEN
       win_type=win_stars(real_space,k,m,rv,rs,hmod,cosm)
    ELSE IF(ifield==field_bound_gas) THEN
       win_type=win_bound_gas(real_space,ifield,k,m,rv,rs,hmod,cosm)
    ELSE IF(ifield==field_free_gas) THEN
       win_type=win_free_gas(real_space,ifield,k,m,rv,rs,hmod,cosm)
    ELSE IF(ifield==field_electron_pressure) THEN
       win_type=win_electron_pressure(real_space,k,m,rv,rs,hmod,cosm)
    ELSE IF(ifield==field_void) THEN
       win_type=win_void(real_space,k,m,rv,rs,hmod,cosm)
    ELSE IF(ifield==field_compensated_void) THEN
       win_type=win_compensated_void(real_space,k,m,rv,rs,hmod,cosm)
    ELSE IF(ifield==field_central_galaxies) THEN
       win_type=win_centrals(real_space,k,m,rv,rs,hmod,cosm)
    ELSE IF(ifield==field_satellite_galaxies) THEN
       win_type=win_satellites(real_space,k,m,rv,rs,hmod,cosm)
    ELSE IF(ifield==field_galaxies) THEN
       win_type=win_galaxies(real_space,k,m,rv,rs,hmod,cosm)
    ELSE IF(ifield==field_HI) THEN
       win_type=win_HI(real_space,k,m,rv,rs,hmod,cosm)
    ELSE IF(ifield==field_cold_gas) THEN
       win_type=win_cold_gas(real_space,ifield,k,m,rv,rs,hmod,cosm)
    ELSE IF(ifield==field_hot_gas) THEN
       win_type=win_hot_gas(real_space,ifield,k,m,rv,rs,hmod,cosm)
    ELSE IF(ifield==field_static_gas) THEN
       win_type=win_static_gas(real_space,ifield,k,m,rv,rs,hmod,cosm)
    ELSE IF(ifield==field_central_stars) THEN
       win_type=win_central_stars(real_space,k,m,rv,rs,hmod,cosm)
    ELSE IF(ifield==field_satellite_stars) THEN
       win_type=win_satellite_stars(real_space,k,m,rv,rs,hmod,cosm)
    ELSE IF(ifield==field_CIB_353 .OR. ifield==field_CIB_545 .OR. ifield==field_CIB_857) THEN
       IF(ifield==field_CIB_353) THEN
          nu=353.e9 ! Frequency [Hz]
       ELSE IF(ifield==field_CIB_545) THEN
          nu=545.e9 ! Frequency [Hz]
       ELSE IF(ifield==field_CIB_857) THEN
          nu=857.e9 ! Frequency [Hz]
       ELSE
          STOP 'WIN_TYPE: Error, ifield specified incorrectly' 
       END IF
       win_type=win_CIB(real_space,nu,k,m,rv,rs,hmod,cosm)
    ELSE IF(ifield==field_halo_11p0_11p5 .OR. ifield==field_halo_11p5_12p0 .OR. &
         ifield==field_halo_12p0_12p5 .OR. ifield==field_halo_12p5_13p0 .OR. &
         ifield==field_halo_13p0_13p5 .OR. ifield==field_halo_13p5_14p0 .OR. &
         ifield==field_halo_14p0_14p5 .OR. ifield==field_halo_14p5_15p0) THEN        
       IF(ifield==field_halo_11p0_11p5) THEN
          mmin=10**11.0 ! Minimum halo mass [Msun/h]
          mmax=10**11.5 ! Maximum halo mass [Msun/h]
       ELSE IF(ifield==field_halo_11p5_12p0) THEN
          mmin=10**11.5 ! Minimum halo mass [Msun/h]
          mmax=10**12.0 ! Maximum halo mass [Msun/h]
       ELSE IF(ifield==field_halo_12p0_12p5) THEN
          mmin=10**12.0 ! Minimum halo mass [Msun/h]
          mmax=10**12.5 ! Maximum halo mass [Msun/h]
       ELSE IF(ifield==field_halo_12p5_13p0) THEN
          mmin=10**12.5 ! Minimum halo mass [Msun/h]
          mmax=10**13.0 ! Maximum halo mass [Msun/h]
       ELSE IF(ifield==field_halo_13p0_13p5) THEN
          mmin=10**13.0 ! Minimum halo mass [Msun/h]
          mmax=10**13.5 ! Maximum halo mass [Msun/h]
       ELSE IF(ifield==field_halo_13p5_14p0) THEN
          mmin=10**13.5 ! Minimum halo mass [Msun/h]
          mmax=10**14.0 ! Maximum halo mass [Msun/h]
       ELSE IF(ifield==field_halo_14p0_14p5) THEN
          mmin=10**14.0 ! Minimum halo mass [Msun/h]
          mmax=10**14.5 ! Maximum halo mass [Msun/h]
       ELSE IF(ifield==field_halo_14p5_15p0) THEN
          mmin=10**14.5 ! Minimum halo mass [Msun/h]
          mmax=10**15.0 ! Maximum halo mass [Msun/h]          
       ELSE
          STOP 'WIN_TYPE: Error, ifield specified incorrectly' 
       END IF
       win_type=win_haloes(real_space,mmin,mmax,k,m,rv,rs,hmod,cosm)
    ELSE
       WRITE(*,*) 'WIN_TYPE: ifield:', ifield
       STOP 'WIN_TYPE: Error, ifield specified incorreclty' 
    END IF

  END FUNCTION win_type

  REAL FUNCTION win_DMONLY(real_space,k,m,rv,rs,hmod,cosm)

    ! Halo profile for all matter under the assumption that it is all CDM
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: real_space
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    REAL, INTENT(IN) :: rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: irho
    REAL :: r, rmin, rmax, c, rss, p1, p2

    ! Set additional halo parameters to zero
    p1=0.
    p2=0.

    rmin=0.
    rmax=rv

    IF(hmod%halo_DMONLY==1) THEN
       ! Analytical NFW
       irho=5 
    ELSE IF(hmod%halo_DMONLY==2) THEN
       ! Non-analyical NFW
       irho=4 
    ELSE IF(hmod%halo_DMONLY==3) THEN
       ! Tophat
       irho=2 
    ELSE IF(hmod%halo_DMONLY==4) THEN
       ! Delta function
       irho=0 
    ELSE IF(hmod%halo_DMONLY==5) THEN
       ! Cored NFW
       irho=24
       p1=hmod%rcore
    ELSE IF(hmod%halo_DMONLY==6) THEN
       ! Isothermal
       irho=1
    ELSE IF(hmod%halo_DMONLY==7) THEN
       ! Shell
       irho=28
    ELSE
       STOP 'WIN_DMONLY: Error, halo_DMONLY specified incorrectly'
    END IF

    ! Force it to use the gravity-only concentration relation (unapply gas correction)
    ! TODO: This is super ugly and should be improved somehow; also is unncessary calculations
    ! TODO: Somehow should store modified and unmodified halo concentrations
    c=rv/rs
    c=c/gas_correction(m,hmod,cosm)
    rss=rv/c

    IF(real_space) THEN
       r=k
       win_DMONLY=rho(r,rmin,rmax,rv,rss,p1,p2,irho)
       win_DMONLY=win_DMONLY/normalisation(rmin,rmax,rv,rss,p1,p2,irho)
    ELSE
       !Properly normalise and convert to overdensity
       win_DMONLY=m*win_norm(k,rmin,rmax,rv,rss,p1,p2,irho)/comoving_matter_density(cosm)
    END IF

  END FUNCTION win_DMONLY

  REAL FUNCTION win_total(real_space,k,m,rv,rs,hmod,cosm)

    ! The halo profile of all the matter
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: real_space
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    REAL, INTENT(IN) :: rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: CDM, gas, stars

    CDM=win_CDM(real_space,k,m,rv,rs,hmod,cosm)
    gas=win_gas(real_space,field_gas,k,m,rv,rs,hmod,cosm)
    stars=win_stars(real_space,k,m,rv,rs,hmod,cosm)

    win_total=CDM+gas+stars

  END FUNCTION win_total

  REAL FUNCTION win_CDM(real_space,k,m,rv,rs,hmod,cosm)

    ! The halo profile for CDM
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: real_space
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    REAL, INTENT(IN) :: rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: irho
    REAL :: r, rmin, rmax, p1, p2, frac

    frac=halo_CDM_fraction(m,hmod,cosm)

    IF(frac==0.) THEN

       win_CDM=0.

    ELSE

       ! Default maximum and minimum radii
       rmin=0.
       rmax=rv

       ! Default additional halo parameters
       p1=0.
       p2=0.

       IF(hmod%halo_CDM==1) THEN      
          irho=5 ! Analytical NFW
       ELSE
          STOP 'WIN_CDM: Error, halo_CDM specified incorrectly'
       END IF

       IF(real_space) THEN
          r=k
          win_CDM=rho(r,rmin,rmax,rv,rs,p1,p2,irho)
          win_CDM=win_CDM/normalisation(rmin,rmax,rv,rs,p1,p2,irho)
       ELSE
          ! Properly normalise and convert to overdensity
          win_CDM=m*win_norm(k,rmin,rmax,rv,rs,p1,p2,irho)/comoving_matter_density(cosm)
       END IF

       win_CDM=frac*win_CDM

    END IF
    
  END FUNCTION win_CDM

  REAL FUNCTION win_gas(real_space,itype,k,m,rv,rs,hmod,cosm)

    ! Halo profile for gas density
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: real_space
    INTEGER, INTENT(IN) :: itype
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    REAL, INTENT(IN) :: rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: win_bound, win_free

    win_bound=win_bound_gas(real_space,itype,k,m,rv,rs,hmod,cosm)
    win_free=win_free_gas(real_space,itype,k,m,rv,rs,hmod,cosm)

    win_gas=win_bound+win_free

  END FUNCTION win_gas

  REAL FUNCTION win_bound_gas(real_space,itype,k,m,rv,rs,hmod,cosm)

    ! Halo profile for bound gas density
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: real_space
    INTEGER, INTENT(IN) :: itype
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    REAL, INTENT(IN) :: rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: win_static, win_cold, win_hot

    win_static=win_static_gas(real_space,itype,k,m,rv,rs,hmod,cosm)
    win_cold=win_cold_gas(real_space,itype,k,m,rv,rs,hmod,cosm)
    win_hot=win_hot_gas(real_space,itype,k,m,rv,rs,hmod,cosm)  

    win_bound_gas=win_static+win_cold+win_hot  

  END FUNCTION win_bound_gas

  REAL FUNCTION win_static_gas(real_space,itype,k,m,rv,rs,hmod,cosm)

    ! Halo profile for the bound gas component
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: real_space
    INTEGER, INTENT(IN) :: itype
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    REAL, INTENT(IN) :: rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: rho0, T0, r, a
    REAL :: rmin, rmax, p1, p2, frac
    INTEGER :: irho_density, irho_electron_pressure

    frac=halo_static_gas_fraction(m,hmod,cosm)

    IF(frac==0.) THEN

       win_static_gas=0.

    ELSE

       ! Initially set the halo parameters to zero
       p1=0.
       p2=0.

       ! Set maximum and minimum integration radius
       rmin=0.
       rmax=rv

       IF(frac<frac_min) THEN
          ! Treat as delta functions if there is not much abundance
          irho_density=0
          irho_electron_pressure=0       
       ELSE IF(hmod%halo_static_gas==1 .OR. hmod%halo_static_gas==3) THEN
          ! Komatsu & Seljak (2001) profile          
          IF(hmod%halo_static_gas==1) THEN
             ! Simplified KS model
             irho_density=11
             irho_electron_pressure=13
          ELSE IF(hmod%halo_static_gas==3) THEN
             ! Full KS model
             irho_density=21
             irho_electron_pressure=23
          END IF          
          p1=HMx_Gamma(m,hmod)         
       ELSE IF(hmod%halo_static_gas==2) THEN
          ! Set cored isothermal profile
          !irho_density=6 ! Isothermal beta model with beta=2/3
          irho_density=15 ! Isothermal beta model with general beta
          p1=HMx_iso_beta(m,hmod)
          irho_electron_pressure=irho_density ! Okay to use density for pressure because temperature is constant (isothermal)
       ELSE       
          STOP 'WIN_STATIC_GAS: Error, halo_static_gas not specified correctly'
       END IF

       IF(itype==field_matter .OR. itype==field_gas .OR. itype==field_bound_gas .OR. itype==field_static_gas) THEN

          ! Density profile of bound gas
          IF(real_space) THEN
             r=k
             win_static_gas=rho(r,rmin,rmax,rv,rs,p1,p2,irho_density)
             win_static_gas=win_static_gas/normalisation(rmin,rmax,rv,rs,p1,p2,irho_density)
          ELSE
             ! Properly normalise and convert to overdensity
             win_static_gas=m*win_norm(k,rmin,rmax,rv,rs,p1,p2,irho_density)/comoving_matter_density(cosm)
          END IF

          win_static_gas=frac*win_static_gas

       ELSE IF(itype==field_electron_pressure) THEN

          ! Electron pressure profile of bound gas
          IF(real_space) THEN
             r=k
             win_static_gas=rho(r,rmin,rmax,rv,rs,p1,p2,irho_electron_pressure)
          ELSE
             ! The electron pressure window is T(r) x rho_e(r), we want unnormalised, so multiply through by normalisation
             ! TODO: Can I make the code more efficient here by having an unnorm window function?
             win_static_gas=win_norm(k,rmin,rmax,rv,rs,p1,p2,irho_electron_pressure)*normalisation(rmin,rmax,rv,rs,p1,p2,irho_electron_pressure)
             !win_static_gas=winint(k,rmin,rmax,rv,rs,p1,p2,irho_electron_pressure)
          END IF

          ! Calculate the value of the density profile prefactor and change units from cosmological to SI
          rho0=m*frac/normalisation(rmin,rmax,rv,rs,p1,p2,irho_density)
          rho0=rho0*msun/mpc/mpc/mpc ! Overflow with REAL*4 if you use mpc**3
          rho0=rho0*cosm%h**2 ! Absorb factors of h, so now [kg/m^3]

          ! Calculate the value of the temperature prefactor [K]
          a=hmod%a
          T0=HMx_alpha(m,hmod)*virial_temperature(m,rv,hmod%a,cosm)

          ! Convert from Temp x density -> electron prsessure (Temp x n; n is all particle number density) 
          win_static_gas=win_static_gas*(rho0/(mp*cosm%mup))*(kb*T0) ! Multiply window by *number density* (all particles) times temperature time k_B [J/m^3]
          win_static_gas=win_static_gas/(eV*(0.01)**(-3)) ! Change units to pressure in [eV/cm^3]
          win_static_gas=win_static_gas*cosm%mup/cosm%mue ! Convert from total thermal pressure to electron pressure

       ELSE

          STOP 'WIN_STATIC_GAS: Error, itype not specified correctly'

       END IF

    END IF

  END FUNCTION win_static_gas

  REAL FUNCTION win_cold_gas(real_space,itype,k,m,rv,rs,hmod,cosm)

    ! Halo profile for the cold gas component
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: real_space
    INTEGER, INTENT(IN) :: itype
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    REAL, INTENT(IN) :: rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: r, rmin, rmax, p1, p2, frac
    INTEGER :: irho

    frac=halo_cold_gas_fraction(m,hmod,cosm)

    IF(frac==0.) THEN

       win_cold_gas=0.

    ELSE

       ! Initially set the halo parameters to zero
       p1=0.
       p2=0.

       ! Set maximum and minimum integration radius
       rmin=0.
       rmax=rv

       IF(frac<frac_min) THEN
          ! Treat as delta functions if there is not much abundance
          irho=0
       ELSE IF(hmod%halo_cold_gas==1) THEN
          ! Delta function
          irho=0       
       ELSE    
          STOP 'WIN_COLD_GAS: Error, halo_cold_gas not specified correctly'
       END IF

       IF(itype==field_matter .OR. itype==field_gas .OR. itype==field_bound_gas .OR. itype==field_cold_gas) THEN

          ! Density profile of cold gas
          IF(real_space) THEN
             r=k
             win_cold_gas=rho(r,rmin,rmax,rv,rs,p1,p2,irho)
             win_cold_gas=win_cold_gas/normalisation(rmin,rmax,rv,rs,p1,p2,irho)
          ELSE
             ! Properly normalise and convert to overdensity
             win_cold_gas=m*win_norm(k,rmin,rmax,rv,rs,p1,p2,irho)/comoving_matter_density(cosm)
          END IF

          win_cold_gas=frac*win_cold_gas

       ELSE IF(itype==field_electron_pressure) THEN

          ! No electron-pressure contribution from the cold gas
          win_cold_gas=0.

       ELSE

          STOP 'WIN_COLD_GAS: Error, itype not specified correctly'

       END IF

    END IF

  END FUNCTION win_cold_gas

  REAL FUNCTION win_hot_gas(real_space,itype,k,m,rv,rs,hmod,cosm)

    ! BUG: Fixed conversion of thermal to electron pressure
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: real_space
    INTEGER, INTENT(IN) :: itype
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    REAL, INTENT(IN) :: rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: r, rmin, rmax, p1, p2, frac
    INTEGER :: irho_density, irho_electron_pressure
    REAL :: rho0, a, T0

    frac=halo_hot_gas_fraction(m,hmod,cosm)

    IF(frac==0.) THEN

       win_hot_gas=0.

    ELSE

       ! Initially set the halo parameters to zero
       p1=0.
       p2=0.

       ! Set maximum and minimum integration radius
       rmin=0.
       rmax=rv

       IF(frac<frac_min) THEN
          ! Treat as delta functions if there is not much abundance
          irho_density=0
          irho_electron_pressure=0
       ELSE IF(hmod%halo_hot_gas==1) THEN
          ! Isothermal
          irho_density=1
          irho_electron_pressure=1
       ELSE    
          STOP 'WIN_HOT_GAS: Error, halo_hot_gas not specified correctly'
       END IF

       IF(itype==field_matter .OR. itype==field_gas .OR. itype==field_bound_gas .OR. itype==field_hot_gas) THEN

          ! Density profile of hot gas
          IF(real_space) THEN
             r=k
             win_hot_gas=rho(r,rmin,rmax,rv,rs,p1,p2,irho_density)
             win_hot_gas=win_hot_gas/normalisation(rmin,rmax,rv,rs,p1,p2,irho_density)
          ELSE
             ! Properly normalise and convert to overdensity
             win_hot_gas=m*win_norm(k,rmin,rmax,rv,rs,p1,p2,irho_density)/comoving_matter_density(cosm)
          END IF

          win_hot_gas=frac*win_hot_gas

       ELSE IF(itype==field_electron_pressure) THEN

          ! Electron-pressure profile of bound gas
          IF(real_space) THEN
             r=k
             win_hot_gas=rho(r,rmin,rmax,rv,rs,p1,p2,irho_electron_pressure)
          ELSE
             ! The electron pressure window is T(r) x rho_e(r), we want unnormalised, so multiply through by normalisation
             ! TODO: Can I make the code more efficient here by having an unnorm window function?
             win_hot_gas=win_norm(k,rmin,rmax,rv,rs,p1,p2,irho_electron_pressure)*normalisation(rmin,rmax,rv,rs,p1,p2,irho_electron_pressure)
          END IF

          ! Calculate the value of the density profile prefactor and change units from cosmological to SI
          rho0=m*frac/normalisation(rmin,rmax,rv,rs,p1,p2,irho_density)
          rho0=rho0*msun/mpc/mpc/mpc ! Overflow with REAL*4 if you use mpc**3
          rho0=rho0*cosm%h**2 ! Absorb factors of h, so now [kg/m^3]

          ! Calculate the value of the temperature prefactor [K]
          a=hmod%a
          T0=HMx_beta(m,hmod)*virial_temperature(m,rv,hmod%a,cosm)

          ! Convert from Temp x density -> electron pressure (Temp x n; n is all particle number density) 
          win_hot_gas=win_hot_gas*(rho0/(mp*cosm%mup))*(kb*T0) ! Multiply window by *number density* (all particles) times temperature time k_B [J/m^3]
          win_hot_gas=win_hot_gas/(eV*(0.01)**(-3)) ! Change units to pressure in [eV/cm^3]
          win_hot_gas=win_hot_gas*cosm%mup/cosm%mue ! Convert from total thermal pressure to electron pressure

       ELSE

          STOP 'WIN_HOT_GAS: Error, itype not specified correctly'

       END IF

    END IF

  END FUNCTION win_hot_gas

  REAL FUNCTION win_free_gas(real_space,itype,k,m,rv,rs,hmod,cosm)

    ! Halo profile for the free gas component
    ! BUG: Fixed conversion of thermal to electron pressure
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: real_space
    INTEGER, INTENT(IN) :: itype
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    REAL, INTENT(IN) :: rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: re, rmin, rmax, r, A, rho0, rhov, T0, p1, p2, beta, c, thing, m0, frac
    INTEGER :: irho_density, irho_electron_pressure

    ! Enable to force the electron pressure to be matched at the virial radius
    ! This is enabled by default for some halo gas/pressure models
    LOGICAL :: match_electron_pressure=.FALSE.

    frac=halo_free_gas_fraction(m,hmod,cosm)

    IF(frac==0.) THEN

       ! Sometimes the free-gas fraction will be zero, in which case this avoids problems
       win_free_gas=0.

    ELSE

       ! Set the halo 'parameter' variables to zero initially
       p1=0.
       p2=0.

       ! Set default min and max radii
       rmin=0.
       rmax=rv

       IF(frac<frac_min) THEN
          
          ! Treat as delta functions if there is not much abundance
          irho_density=0
          irho_electron_pressure=0

       ELSE IF(hmod%halo_free_gas==1) THEN

          ! Simple isothermal model, motivated by constant velocity and rate expulsion
          irho_density=1
          irho_electron_pressure=irho_density ! Okay because T is constant
          rmin=0.
          rmax=2.*rv

       ELSE IF(hmod%halo_free_gas==2) THEN

          ! Ejected gas model from Schneider & Teyssier (2015)
          irho_density=10
          irho_electron_pressure=irho_density ! Okay because T is constant
          rmin=rv
          re=rv
          p1=re
          rmax=15.*re ! Needs to be such that integral converges (15rf seems okay)

       ELSE IF(hmod%halo_free_gas==3) THEN

          ! Now do isothermal shell connected to the KS profile continuously
          irho_density=16
          irho_electron_pressure=irho_density ! Okay because T is constant

          ! Isothermal model with continuous link to KS
          rhov=win_static_gas(.TRUE.,2,rv,m,rv,rs,hmod,cosm) ! This is the value of the density at the halo boundary for the bound gas           
          A=rhov/rho(rv,0.,rv,rv,rs,p1,p2,irho_density) ! This is A, as in A/r^2

          rmin=rv
          rmax=rv+frac/(4.*pi*A) ! This ensures density continuity and mass conservation

          c=10. ! How many times larger than the virial radius can the gas cloud go?          
          IF(rmax>c*rv) rmax=c*rv ! This needs to be set otherwise get huge decrement in gas power at large scales
          match_electron_pressure=.TRUE. ! Match the electron pressure at the boundary

       ELSE IF(hmod%halo_free_gas==4) THEN

          ! Ejected gas is a continuation of the KS profile
          irho_density=11 ! KS
          irho_electron_pressure=13 ! KS
          rmin=rv
          rmax=2.*rv
          p1=HMx_Gamma(m,hmod)

       ELSE IF(hmod%halo_free_gas==5) THEN

          m0=1e14

          IF(m<m0) THEN

             irho_density=0
             irho_electron_pressure=irho_density
             rmin=0.
             rmax=rv

          ELSE

             ! Set the density profile to be the power-law profile
             irho_density=17
             irho_electron_pressure=irho_density ! Not okay

             ! Calculate the KS index at the virial radius
             c=rv/rs
             beta=(c-(1.+c)*log(1.+c))/((1.+c)*log(1.+c))
             beta=beta/(HMx_Gamma(m,hmod)-1.) ! This is the power-law index at the virial radius for the KS gas profile
             p1=beta
             !WRITE(*,*) 'Beta:', beta, log10(m)
             IF(beta<=-3.) beta=-2.9 ! If beta<-3 then there is only a finite amount of gas allowed in the free component

             ! Calculate the density at the boundary of the KS profile
             rhov=win_static_gas(.TRUE.,2,rv,m,rv,rs,hmod,cosm)
             !WRITE(*,*) 'rho_v:', rhov

             ! Calculate A as in rho(r)=A*r**beta
             A=rhov/rho(rv,0.,rv,rv,rs,p1,p2,irho_density)
             !WRITE(*,*) 'A:', A

             ! Set the minimum radius for the power-law to be the virial radius
             rmin=rv
             !WRITE(*,*) 'rmin:', rmin

             ! Set the maximum radius so that it joins to KS profile seamlessly
             thing=(beta+3.)*frac/(4.*pi*A)+(rhov*rv**3)/A
             !WRITE(*,*) 'thing:', thing
             IF(thing>0.) THEN
                ! This then fixes the condition of contiunity in amplitude and gradient
                rmax=thing**(1./(beta+3.))
             ELSE
                ! If there are no sohmodions then fix to 10rv and accept discontinuity
                ! There may be no sohmodion if there is a lot of free gas and if beta<-3
                rmax=10.*rv
             END IF
             !WRITE(*,*) 'rmax 2:', rmax

          END IF

       ELSE IF(hmod%halo_free_gas==6) THEN

          ! Cubic profile
          rmin=rv
          rmax=3.*rv
          irho_density=18
          irho_electron_pressure=irho_density

       ELSE IF(hmod%halo_free_gas==7) THEN

          ! Smooth profile (rho=0)
          rmin=0.
          rmax=rv
          irho_density=19
          irho_electron_pressure=irho_density

       ELSE IF(hmod%halo_free_gas==8) THEN

          ! Delta function
          rmin=0.
          rmax=rv
          irho_density=0
          irho_electron_pressure=irho_density

       ELSE
          STOP 'WIN_FREE_GAS: Error, halo_free_gas specified incorrectly'
       END IF

       ! Density profile
       IF(itype==field_matter .OR. itype==field_gas .OR. itype==field_free_gas) THEN

          ! Density profile of free gas
          IF(real_space) THEN
             r=k
             win_free_gas=rho(r,rmin,rmax,rv,rs,p1,p2,irho_density)
             win_free_gas=win_free_gas/normalisation(rmin,rmax,rv,rs,p1,p2,irho_density)
          ELSE
             ! Properly normalise and convert to overdensity
             win_free_gas=m*win_norm(k,rmin,rmax,rv,rs,p1,p2,irho_density)/comoving_matter_density(cosm)
          END IF

          win_free_gas=frac*win_free_gas

          ! Electron pressure profile
       ELSE IF(itype==field_electron_pressure) THEN

          ! If we are applying a pressure-matching condition
          IF(match_electron_pressure) THEN

             STOP 'WIN_FREE_GAS: Check the units and stuff here *very* carefully'

             r=k
             IF(r>rmin .AND. r<rmax) THEN
                ! Only works for isothermal profile
                win_free_gas=win_static_gas(.TRUE.,6,rv,m,rv,rs,hmod,cosm)*(r/rv)**(-2)
             ELSE
                win_free_gas=0.
             END IF

          ELSE

             ! Electron pressure profile of free gas
             IF(real_space) THEN
                r=k
                win_free_gas=rho(r,rmin,rmax,rv,rs,p1,p2,irho_electron_pressure)
             ELSE  
                win_free_gas=win_norm(k,rmin,rmax,rv,rs,p1,p2,irho_electron_pressure)*normalisation(rmin,rmax,rv,rs,p1,p2,irho_electron_pressure)
             END IF

             ! Calculate the value of the density profile prefactor [(Msun/h)/(Mpc/h)^3] and change units from cosmological to SI
             rho0=m*frac/normalisation(rmin,rmax,rv,rs,p1,p2,irho_density) ! rho0 in [(Msun/h)/(Mpc/h)^3]
             rho0=rho0*msun/Mpc/Mpc/Mpc ! Overflow with REAL(4) if you use Mpc**3, this converts to SI units [h^2 kg/m^3]
             rho0=rho0*cosm%h**2 ! Absorb factors of h, so now [kg/m^3]

             ! This is the total thermal pressure of the WHIM
             T0=HMx_Twhim(hmod) ! [K]

             ! Factors to convert from Temp x density -> electron pressure (Temp x n; n is all particle number density) 
             win_free_gas=win_free_gas*(rho0/(mp*cosm%mup))*(kb*T0) ! Multiply window by *number density* (all particles) times temperature time k_B [J/m^3]
             win_free_gas=win_free_gas/(eV*(0.01)**(-3)) ! Change units to pressure in [eV/cm^3]
             win_free_gas=win_free_gas*cosm%mup/cosm%mue ! Convert from total thermal pressure to electron pressure

          END IF

       ELSE

          STOP 'WIN_FREE_GAS: Error, itype not specified correctly'

       END IF

    END IF

  END FUNCTION win_free_gas

  REAL FUNCTION win_stars(real_space,k,m,rv,rs,hmod,cosm)

    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: real_space
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    REAL, INTENT(IN) :: rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: win_central, win_satellite

    win_central=win_central_stars(real_space,k,m,rv,rs,hmod,cosm)
    win_satellite=win_satellite_stars(real_space,k,m,rv,rs,hmod,cosm)

    win_stars=win_central+win_satellite

  END FUNCTION win_stars

  REAL FUNCTION win_central_stars(real_space,k,m,rv,rs,hmod,cosm)

    ! Halo profile for stars
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: real_space
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    REAL, INTENT(IN) :: rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: irho
    REAL :: rstar, r, rmin, rmax, p1, p2, frac

    frac=halo_central_star_fraction(m,hmod,cosm)

    IF(frac==0.) THEN

       win_central_stars=0.

    ELSE

       ! Initially set p1, p2
       p1=0.
       p2=0.

       ! Set default maximuim and minimum radii
       rmin=0.
       rmax=rv

       IF(frac<frac_min) THEN         
          ! Treat as delta functions if there is not much abundance
          irho=0
       ELSE IF(hmod%halo_central_stars==1) THEN
          ! Fedeli (2014)
          irho=7
          rstar=rv/HMx_cstar(m,hmod)
          p1=rstar
          rmax=rv ! Set so that not too much bigger than rstar, otherwise bumps integration goes tits
       ELSE IF(hmod%halo_central_stars==2) THEN
          ! Schneider & Teyssier (2015), following Mohammed (2014)
          irho=9
          rstar=rv/HMx_cstar(m,hmod)
          p1=rstar
          rmax=min(10.*rstar,rv) ! Set so that not too much bigger than rstar, otherwise bumps integration goes crazy
       ELSE IF(hmod%halo_central_stars==3) THEN
          ! Delta function
          irho=0
       ELSE IF(hmod%halo_central_stars==4) THEN
          ! Transition mass between NFW and delta function
          ! NOTE: mstar here is the same as in the stellar halo-mass fraction. It should probably not be this
          IF(m<HMx_Mstar(hmod)) THEN
             irho=0 ! Delta function
          ELSE
             irho=5 ! NFW
          END IF
       ELSE
          STOP 'WIN_CENTRAL_STARS: Error, halo_central_stars specified incorrectly'
       END IF

       IF(real_space) THEN
          r=k
          win_central_stars=rho(r,rmin,rmax,rv,rs,p1,p2,irho)
          win_central_stars=win_central_stars/normalisation(rmin,rmax,rv,rs,p1,p2,irho)
       ELSE
          ! Properly normalise and convert to overdensity
          win_central_stars=m*win_norm(k,rmin,rmax,rv,rs,p1,p2,irho)/comoving_matter_density(cosm)
       END IF

       win_central_stars=frac*win_central_stars

    END IF

  END FUNCTION win_central_stars

  REAL FUNCTION win_satellite_stars(real_space,k,m,rv,rs,hmod,cosm)

    ! Halo profile for stars
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: real_space
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    REAL, INTENT(IN) :: rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: irho
    REAL :: rstar, r, rmin, rmax, p1, p2, frac

    frac=halo_satellite_star_fraction(m,hmod,cosm)

    IF(frac==0.) THEN

       win_satellite_stars=0.

    ELSE

       ! Initially set p1, p2
       p1=0.
       p2=0.

       ! Set default maximuim and minimum radii
       rmin=0.
       rmax=rv

       IF(frac<frac_min) THEN         
          ! Treat as delta functions if there is not much abundance
          irho=0
       ELSE IF(hmod%halo_satellite_stars==1) THEN
          ! NFW
          irho=5
       ELSE IF(hmod%halo_satellite_stars==2) THEN
          ! Fedeli (2014)
          irho=7
          rstar=rv/HMx_cstar(m,hmod)
          p1=rstar
          rmax=rv ! Set so that not too much bigger than rstar, otherwise bumps integration goes tits
       ELSE
          STOP 'WIN_SATELLITE_STARS: Error, halo_satellite_stars specified incorrectly'
       END IF

       IF(real_space) THEN
          r=k
          win_satellite_stars=rho(r,rmin,rmax,rv,rs,p1,p2,irho)
          win_satellite_stars=win_satellite_stars/normalisation(rmin,rmax,rv,rs,p1,p2,irho)
       ELSE
          ! Properly normalise and convert to overdensity
          win_satellite_stars=m*win_norm(k,rmin,rmax,rv,rs,p1,p2,irho)/comoving_matter_density(cosm)
       END IF

       win_satellite_stars=frac*win_satellite_stars

    END IF

  END FUNCTION win_satellite_stars

  REAL FUNCTION win_electron_pressure(real_space,k,m,rv,rs,hmod,cosm)

    !Halo electron pressure profile function for the sum of bound + unbound electron gas
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: real_space
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    REAL, INTENT(IN) :: rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(hmod%electron_pressure==1) THEN
       ! This overrides everything and just uses the UPP
       win_electron_pressure=UPP(real_space,k,m,rv,rs,hmod,cosm)
    ELSE IF(hmod%electron_pressure==2) THEN
       ! Otherwise use...
       !win_electron_pressure=win_bound_gas(real_space,itype_pressure,k,m,rv,rs,hmod,cosm)+win_free_gas(real_space,itype_pressure,k,m,rv,rs,hmod,cosm)
       win_electron_pressure=win_gas(real_space,field_electron_pressure,k,m,rv,rs,hmod,cosm)
    ELSE
       STOP 'WIN_ELECTRON_PRESSURE: Error, electron_pressure specified incorrectly'
    END IF

  END FUNCTION win_electron_pressure

  REAL FUNCTION win_void(real_space,k,m,rv,rs,hmod,cosm)

    !Void profile
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: real_space
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    REAL, INTENT(IN) :: rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: irho
    REAL :: r, rmin, rmax, p1, p2

    ! Initially set p1, p2
    p1=0.
    p2=0.

    ! Set default maximuim and minimum radii
    rmin=0.
    rmax=rv

    IF(hmod%halo_void==1) THEN
       !Top-hat void
       irho=2
       rmin=0.
       rmax=10.*rv
    ELSE
       STOP 'WIN_VOID: Error, halo_void specified incorrectly'
    END IF

    IF(real_space) THEN
       r=k
       win_void=rho(r,rmin,rmax,rv,rs,p1,p2,irho)
       win_void=win_void/normalisation(rmin,rmax,rv,rs,p1,p2,irho)
    ELSE       
       win_void=m*win_norm(k,rmin,rmax,rv,rs,p1,p2,irho)/comoving_matter_density(cosm)
    END IF

  END FUNCTION win_void

  REAL FUNCTION win_compensated_void(real_space,k,m,rv,rs,hmod,cosm)

    !Profile for compensated voids
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: real_space
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    REAL, INTENT(IN) :: rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: irho
    REAL :: r, rmin, rmax, p1, p2

    ! Initially set p1, p2
    p1=0.
    p2=0.

    ! Set default maximuim and minimum radii
    rmin=0.
    rmax=rv

    IF(hmod%halo_compensated_void==1) THEN
       !Top-hat
       irho=2
       rmin=0.
       rmax=10.*rv
    ELSE
       STOP 'WIN_COMPENSATED_VOID: Error, halo_compensated_void specified incorrectly'
    END IF

    IF(real_space) THEN
       r=k
       win_compensated_void=rho(r,rmin,rmax,rv,rs,p1,p2,irho)
       win_compensated_void=win_compensated_void/normalisation(rmin,rmax,rv,rs,p1,p2,irho)
    ELSE       
       win_compensated_void=m*win_norm(k,rmin,rmax,rv,rs,p1,p2,irho)/comoving_matter_density(cosm)
    END IF

  END FUNCTION win_compensated_void

  REAL FUNCTION win_galaxies(real_space,k,m,rv,rs,hmod,cosm)

    ! Halo profile for all galaxies
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: real_space
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    REAL, INTENT(IN) :: rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm

    win_galaxies=win_centrals(real_space,k,m,rv,rs,hmod,cosm)+win_satellites(real_space,k,m,rv,rs,hmod,cosm)

  END FUNCTION win_galaxies

  REAL FUNCTION win_centrals(real_space,k,m,rv,rs,hmod,cosm)

    !Halo profile for central galaxies
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: real_space
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    REAL, INTENT(IN) :: rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: irho
    REAL :: r, rmin, rmax, p1, p2, N

    IF(hmod%has_galaxies .EQV. .FALSE.) CALL init_galaxies(hmod,cosm)

    N=N_centrals(m,hmod)

    IF(N==0.) THEN

       win_centrals=0.

    ELSE

       ! Default minimum and maximum radii
       rmin=0.
       rmax=rv

       ! Default additional halo parameters
       p1=0.
       p2=0.

       ! Delta functions
       irho=0

       IF(real_space) THEN
          r=k
          win_centrals=rho(r,rmin,rmax,rv,rs,p1,p2,irho)
          win_centrals=win_centrals/normalisation(rmin,rmax,rv,rs,p1,p2,irho)
       ELSE      
          win_centrals=win_norm(k,rmin,rmax,rv,rs,p1,p2,irho)/hmod%n_c
       END IF

       win_centrals=N*win_centrals

    END IF

  END FUNCTION win_centrals

  REAL FUNCTION win_haloes(real_space,mmin,mmax,k,m,rv,rs,hmod,cosm)

    ! Halo profile function for haloes
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: real_space
    REAL, INTENT(IN) :: mmin
    REAL, INTENT(IN) :: mmax
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    REAL, INTENT(IN) :: rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: irho
    REAL :: r, rmin, rmax, p1, p2, N, nhalo
    REAL :: nu1, nu2
    
    IF(m<mmin .OR. m>mmax) THEN
       N=0.
    ELSE
       N=1.
    END IF

    IF(N==0.) THEN

       win_haloes=0.

    ELSE

       ! Default minimum and maximum radii
       rmin=0.
       rmax=rv

       ! This shitty calculation really only needs to be done once
       ! This could slow down the calculation by a large amount
       ! TODO: Ensure this is only done once
       nu1=nu_M(mmin,hmod,cosm)
       nu2=nu_M(mmax,hmod,cosm)
       nhalo=mean_halo_number_density(nu1,nu2,hmod,cosm)

       ! Default additional halo parameters
       p1=0.
       p2=0.

       ! Delta function
       irho=0

       IF(real_space) THEN
          r=k
          win_haloes=rho(r,rmin,rmax,rv,rs,p1,p2,irho)
          win_haloes=win_haloes/normalisation(rmin,rmax,rv,rs,p1,p2,irho)
       ELSE      
          win_haloes=win_norm(k,rmin,rmax,rv,rs,p1,p2,irho)/nhalo
       END IF

    END IF

  END FUNCTION win_haloes

  REAL FUNCTION win_satellites(real_space,k,m,rv,rs,hmod,cosm)

    ! Halo profile for satellite galaxies
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: real_space
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    REAL, INTENT(IN) :: rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: irho
    REAL :: r, rmin, rmax, p1, p2, N

    IF(hmod%has_galaxies .EQV. .FALSE.) CALL init_galaxies(hmod,cosm)

    N=N_satellites(m,hmod)

    IF(N==0.) THEN

       win_satellites=0.

    ELSE

       ! Initially set p1, p2
       p1=0.
       p2=0.

       ! Default maximum and minimum radii
       rmin=0.
       rmax=rv

       ! NFW profile
       irho=5

       IF(real_space) THEN
          r=k
          win_satellites=rho(r,rmin,rmax,rv,rs,p1,p2,irho)
          win_satellites=win_satellites/normalisation(rmin,rmax,rv,rs,p1,p2,irho)
       ELSE
          win_satellites=win_norm(k,rmin,rmax,rv,rs,p1,p2,irho)/hmod%n_s
       END IF

       win_satellites=N*win_satellites

    END IF

  END FUNCTION win_satellites

  REAL FUNCTION win_CIB(real_space,nu,k,m,rv,rs,hmod,cosm)

    ! Halo profile for all matter under the assumption that it is all CDM
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: real_space
    REAL, INTENT(IN) :: nu
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    REAL, INTENT(IN) :: rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: irho
    REAL :: r, rmin, rmax, p1, p2, z
    REAL :: crap
    
    REAL, PARAMETER :: a=1. ! Dust blob size relative to halo virial radius
    REAL, PARAMETER :: T=15. ! Dust temperature [K]
    REAL, PARAMETER :: beta=1.6 ! Grey-body power-law index

    ! Set additional halo parameters to zero
    p1=0.
    p2=0.

    rmin=0.
    rmax=rv

    ! Prevent compile warnings
    crap=m
    crap=cosm%Om_m

    ! Halo type
    ! 0 - Delta function
    ! 1 - Isothermal
    ! 5 - NFW
    irho=1

    IF(real_space) THEN
       r=k
       win_CIB=rho(r,rmin,rmax,rv,rs,p1,p2,irho)
       win_CIB=win_CIB/normalisation(rmin,rmax,rv,rs,p1,p2,irho)
    ELSE
       !Properly normalise and convert to overdensity
       win_CIB=win_norm(k,rmin,rmax,rv,rs,p1,p2,irho)
    END IF

    z=hmod%z
    win_CIB=grey_body_nu((1.+z)*nu,T,beta) ! Get the black-body radiance [W m^-2 Sr^-1 Hz^-1]
    win_CIB=win_CIB/Jansky     ! Convert units to Jansky [Jy Sr^-1]
    win_CIB=win_CIB*(a*rv)**2  ! [Jy Sr^-1 (Mpc/h)^2]
    !win_CIB=win_CIB/(1e-3+luminosity_distance(a,cosm))**2 ! Bad idea because divide by zero when z=0

  END FUNCTION win_CIB

  REAL FUNCTION grey_body_nu(nu,T,beta)

    ! Grey body irradiance [W m^-2 Hz^-1 Sr^-1]
    USE physics
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu ! Observing frequency [Hz]
    REAL, INTENT(IN) :: T !  Grey body temperature [K]
    REAL, INTENT(IN) :: beta ! Power-law index
    REAL, PARAMETER :: tau=1. ! Emissivity (degenerate with R in CIB work)
    REAL, PARAMETER :: nu0=545.e9 ! Reference frequency [Hz]
    
    grey_body_nu=tau*((nu/nu0)**beta)*black_body_nu(nu,T)
    
  END FUNCTION grey_body_nu

  REAL FUNCTION N_centrals(m,hmod)

    ! The number of central galaxies as a function of halo mass
    IMPLICIT NONE
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod

    IF(m>hmod%mhalo_min .AND. m<hmod%mhalo_max) THEN
       N_centrals=1
    ELSE
       N_centrals=0
    END IF

  END FUNCTION N_centrals

  REAL FUNCTION N_satellites(m,hmod)

    ! The number of satellite galxies as a function of halo mass
    IMPLICIT NONE
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod

    IF(m<hmod%mhalo_min) THEN
       N_satellites=0
    ELSE
       N_satellites=ceiling(m/hmod%mhalo_min)-1
    END IF

  END FUNCTION N_satellites

  REAL FUNCTION N_galaxies(m,hmod)

    ! The number of central galaxies as a function of halo mass
    IMPLICIT NONE
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod

    N_galaxies=N_centrals(m,hmod)+N_satellites(m,hmod)

  END FUNCTION N_galaxies

  REAL FUNCTION win_HI(real_space,k,m,rv,rs,hmod,cosm)

    ! Returns the real or Fourier space HI halo profile
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: real_space ! Real space or Fourier space
    REAL, INTENT(IN) :: k  ! Comoving wave vector (or radius)
    REAL, INTENT(IN) :: m  ! Halo mass
    REAL, INTENT(IN) :: rv ! Halo virial radius
    REAL, INTENT(IN) :: rs ! Halo scale radius 
    TYPE(halomod), INTENT(INOUT) :: hmod ! Halomodel
    TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology
    INTEGER :: irho
    REAL :: r, rmin, rmax, p1, p2, frac
    REAL :: r0, alpha, c_HI, r_HI, z

    IF(hmod%has_HI .EQV. .FALSE.) CALL init_HI(hmod,cosm)

    frac=halo_HI_fraction(m,hmod,cosm)

    IF(frac==0.) THEN

       win_HI=0.

    ELSE

       ! Default minimum and maximum radii
       rmin=0.
       rmax=rv

       ! Default additional halo parameters
       p1=0.
       p2=0.

       IF(hmod%halo_HI==1) THEN
          ! NFW profile
          irho=5
       ELSE IF(hmod%halo_HI==2) THEN
          ! Delta function
          irho=0
       ELSE IF(hmod%halo_HI==3) THEN
          ! Polynomial with exponential cut off Villaescusa-Navarro et al. (1804.09180)
          irho=25
          r0=10**(-2.5) ! 0.003 Mpc/h (really small)
          alpha=3.00 ! Tending to homogeneity  
          p1=r0
          p2=alpha
       ELSE IF(hmod%halo_HI==4) THEN
          ! Modified NFW with exponential cut off Villaescusa-Navarro et al. (1804.09180)
          irho=26
          r0=10**(-2.5) ! 0.003 Mpc/h (really small)
          r_HI=10**(-3.0)    
          p1=r0
          p2=r_HI
       ELSE IF(hmod%halo_HI==5) THEN
          ! Modified NFW from Padmanabhan & Refreiger (1607.01021)
          irho=27
          z=hmod%z
          c_HI=130.
          c_HI=4.*c_HI*((M/1e11)**(-0.109))/(1.+z)
          r_HI=rv/c_HI
          p1=r_HI
       ELSE
          STOP 'win_HI: Error, halo_HI not specified correctly'
       END IF

       IF(real_space) THEN
          r=k
          win_HI=rho(r,rmin,rmax,rv,rs,p1,p2,irho)
          win_HI=win_HI/normalisation(rmin,rmax,rv,rs,p1,p2,irho)
       ELSE
          !win_HI=win_norm(k,rmin,rmax,rv,rs,p1,p2,irho) ! Wrong, but is what I first sent Richard and Kiyo
          win_HI=m*win_norm(k,rmin,rmax,rv,rs,p1,p2,irho)/hmod%rho_HI
       END IF

       win_HI=frac*win_HI

    END IF

  END FUNCTION win_HI

  REAL FUNCTION virial_temperature(M,rv,a,cosm)

    ! Halo physical virial temperature [K]
    ! Calculates the temperature as if pristine gas falls into the halo
    ! Energy is equally distributed between the particles
    IMPLICIT NONE
    REAL :: M ! virial mass
    REAL, INTENT(IN) :: rv ! virial radius
    REAL, INTENT(IN) :: a ! scale factor
    TYPE(cosmology), INTENT(INOUT) :: cosm ! cosmology

    REAL, PARAMETER :: modes=3. ! 1/2 k_BT per mode, 3 modes for 3 dimensions 

    virial_temperature=bigG*((M*msun)*mp*cosm%mup)/(a*rv*mpc) ! NEW: a in denominator converts comoving->physical halo radius
    virial_temperature=virial_temperature/(kb*modes/2.) ! Convert to temperature from energy

  END FUNCTION virial_temperature

  REAL FUNCTION UPP(real_space,k,m,rv,rs,hmod,cosm)

    ! Universal electron pressure profile (Arnaud et al. 2010; arxiv:0910.1234)
    ! Note *very* well that this is for *electron* pressure (see Arnaud 2010 and my notes on this)
    ! Note that is is also the physical pressure, and relates to the comoving pressure via (1+z)^3
    ! The units of the pressure profile are [eV/cm^3]
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: real_space
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    REAL, INTENT(IN) :: rs
    TYPE(halomod), INTENT(INOUT) :: hmod 
    TYPE(cosmology), INTENT(INOUT) :: cosm 
    REAL :: r500c, rmin, rmax, a, z, r, m500c, E

    REAL, PARAMETER :: alphap=0.12 ! Exponent correction
    REAL, PARAMETER :: b=0.        ! Hydrostatic mass bias
    INTEGER, PARAMETER :: irho=14  ! Set UPP profile

    IF(hmod%has_mass_conversions .EQV. .FALSE.) CALL convert_mass_definitions(hmod,cosm)

    ! Get r500 for UPP
    r500c=exp(find(log(m),hmod%log_m,log(hmod%r500c),hmod%n,3,3,2)) ! [Mpc/h]

    ! Set the radius range for the profile
    rmin=0.
    rmax=rv

    IF(real_space) THEN
       r=k
       UPP=rho(r,rmin,rmax,rv,rs,r500c,zero,irho)
    ELSE
       UPP=winint(k,rmin,rmax,rv,rs,r500c,zero,irho,imeth_win)
    END IF

    ! UPP, P(x), equation 4.1 in Ma et al. (2015)
    m500c=exp(find(log(m),hmod%log_m,log(hmod%m500c),hmod%n,3,3,2)) ![Msun/h]
    m500c=m500c*(1.-b) ![Msun/h]

    ! Dimensionless Hubble parameter
    a=hmod%a ! Scale-factor appears explicitly here
    E=sqrt(Hubble2(a,cosm))

    ! Pre-factors from equation 4.1 in Ma et al. (2015) [eV cm^-3 - no h factors]
    UPP=UPP*((m500c/2.1e14)**(alphap+2./3.))*(E**(8./3.))*1.65*(cosm%h/0.7)**2

    ! The standard UPP is written is in physical units
    ! scale to comoving using (1+z)^3 because pressure ~energy density
    z=hmod%z ! Redshift appears explicitly here
    UPP=UPP/(1.+z)**3

  END FUNCTION UPP

  FUNCTION rho(r,rmin,rmax,rv,rs,p1,p2,irho)

    ! This is an UNNORMALISED halo profile of any sort

    ! Types of profile
    ! ================
    !  0 - Delta function at r=0
    !  1 - Isothermal: r^-2
    !  2 - Top hat: r^0
    !  3 - Moore (1999)
    !  4 - NFW (1997)
    !  5 - Analytic NFW
    !  6 - Beta model with beta=2/3
    !  7 - Star profile
    !  8 - Komatsu & Seljak (2001) according to Schneider & Teyssier (2015)
    !  9 - Stellar profile from Schneider & Teyssier (2015)
    ! 10 - Ejected gas profile (Schneider & Teyssier 2015)
    ! 11 - Simplified Komatsu & Seljak (2001) density
    ! 12 - Simplified Komatsu & Seljak (2001) temperature
    ! 13 - Simplified Komatsu & Seljak (2001) pressure
    ! 14 - Universal pressure profile
    ! 15 - Isothermal beta model, beta=0.86 (Ma et al. 2015)
    ! 16 - Isothermal exterior
    ! 17 - Power-law profile
    ! 18 - Cubic profile: r^-3
    ! 19 - Smooth profile (rho = 0, not really physical)
    ! 20 - Exponential profile
    ! 21 - Full Komatsu & Seljak (2001) density
    ! 22 - Full Komatsu & Seljak (2001) temperature
    ! 23 - Full Komatsu & Seljak (2001) pressure
    ! 24 - Cored NFW profile (Copeland, Taylor & Hall 2018)
    ! 25 - Polynomial with central hole (Villaescusa-Navarro et al. 2018)
    ! 26 - Modified NFW with central hole (Villaescusa-Navarro et al. 2018)
    ! 27 - Modified NFW (Padmanabhan & Refregier 2018)
    ! 28 - Shell

    IMPLICIT NONE
    REAL :: rho
    REAL, INTENT(IN) :: r, rmin, rmax, rv, rs, p1, p2 ! Standard profile parameters
    INTEGER, INTENT(IN) :: irho
    REAL :: y, ct, t, c, beta, Gamma, r500c, rt, A, re, rstar, B, rb, r0, alpha, rh, eta0
    REAL :: f1, f2
    REAL :: crap

    ! UPP parameters
    REAL, PARAMETER :: P0=6.41
    REAL, PARAMETER :: c500=1.81
    REAL, PARAMETER :: alpha_UPP=1.33
    REAL, PARAMETER :: beta_UPP=4.13
    REAL, PARAMETER :: gamma_UPP=0.31

    ! To stop compile-time warnings
    crap=p2

    IF(r<rmin .OR. r>rmax) THEN
       ! The profile is considered to be zero outside this region
       rho=0.
    ELSE
       IF(irho==0) THEN
          ! Delta function
          ! Do not assign any value to rho as this gets handled properly elsewhere
          ! STOP 'RHO: You should not be here for a delta-function profile'
          rho=0. 
       ELSE IF(irho==1 .OR. irho==16) THEN
          ! Isothermal
          rho=1./r**2
       ELSE IF(irho==2) THEN
          ! Top hat
          rho=1.
       ELSE IF(irho==3) THEN
          ! Moore (1999)
          y=r/rs
          rho=1./((y**1.5)*(1.+y**1.5))
       ELSE IF(irho==4 .OR. irho==5) THEN
          ! NFW (1997)
          y=r/rs
          rho=1./(y*(1.+y)**2)
       ELSE IF(irho==6) THEN
          ! Isothermal beta model (X-ray gas; SZ profiles; beta=2/3 fixed)
          ! Also known as 'cored isothermal profile'
          ! In general this is (1+(r/rs)**2)**(-3*beta/2); beta=2/3 cancels the last power
          rho=1./(1.+(r/rs)**2)
       ELSE IF(irho==7) THEN
          ! Stellar profile from Fedeli (2014a)
          rstar=p1
          y=r/rstar
          rho=(1./y)*exp(-y)
       ELSE IF(irho==8) THEN
          ! Komatsu & Seljak (2001) profile with NFW transition radius
          ! VERY slow to calculate the W(k) for some reason
          ! Also creates a weird upturn in P(k) that I do not think can be correct
          STOP 'RHO: This is fucked'
          t=sqrt(5.)
          rt=rv/t
          y=r/rs
          c=rs/rv
          ct=c/t
          Gamma=(1.+3.*ct)*log(1.+ct)/((1.+ct)*log(1.+ct)-ct)
          IF(r<=rt) THEN
             ! Komatsu Seljak in the interior
             rho=(log(1.+y)/y)**Gamma
          ELSE
             ! NFW in the outskirts
             A=((rt/rs)*(1.+rt/rs)**2)*(log(1.+rt/rs)/(rt/rs))**Gamma
             rho=A/(y*(1.+y)**2)
          END IF
       ELSE IF(irho==9) THEN
          ! Stellar profile from Schneider & Teyssier (2015) via Mohammed (2014)
          rstar=p1
          rho=exp(-(r/(2.*rstar))**2)/r**2
          ! Converting to y caused the integration to crash for some reason !?!
          !y=r/rs
          !rho=exp(-(y/2.)**2.)/y**2.
       ELSE IF(irho==10) THEN
          ! Ejected gas profile from Schneider & Teyssier (2015)
          re=p1
          rho=exp(-0.5*(r/re)**2)
       ELSE IF(irho==11 .OR. irho==12 .OR. irho==13 .OR. irho==21 .OR. irho==22 .OR. irho==23) THEN
          ! Komatsu & Seljak (2001) profile
          !Gamma=1.18 ! Recommended by Rabold (2017)          
          IF(irho==11 .OR. irho==12 .OR. irho==13) THEN
             Gamma=p1
             y=r/rs
             rho=log(1.+y)/y
          ELSE IF(irho==21 .OR. irho==22 .OR. irho==23) THEN
             c=rv/rs
             Gamma=p1+0.01*(c-6.5)
             !WRITE(*,*) 'Gamma:', Gamma
             eta0=0.00676*(c-6.5)**2+0.206*(c-6.5)+2.48
             !WRITE(*,*) 'eta0:', eta0
             f1=(3./eta0)*(Gamma-1.)/Gamma
             !WRITE(*,*) 'f1:', f1
             f2=c/NFW_factor(c)
             !WRITE(*,*) 'f2:', f2
             B=f1*f2
             IF(B>1.) B=1.
             !WRITE(*,*) 'B:', B
             y=r/rs
             !WRITE(*,*) 'y:', y
             rho=1.-B*(1.-log(1.+y)/y)
             !WRITE(*,*) 'rho:', rho
          ELSE
             STOP 'RHO: Error, irho specified incorrectly'
          END IF
          IF(irho==11 .OR. irho==21) THEN
             ! KS density profile
             rho=rho**(1./(Gamma-1.))
          ELSE IF(irho==12 .OR. irho==22) THEN
             ! KS temperature profile
             rho=rho
          ELSE IF(irho==13 .OR. irho==23) THEN
             ! KS pressure profile
             rho=rho**(Gamma/(Gamma-1.))
          END IF
       ELSE IF(irho==14) THEN
          ! UPP is in terms of r500c, not rv
          r500c=p1
          ! UPP funny-P(x), equation 4.2 in Ma et al. (2015)
          f1=(c500*r/r500c)**gamma_UPP
          f2=(1.+(c500*r/r500c)**alpha_UPP)**((beta_UPP-gamma_UPP)/alpha_UPP)
          rho=P0/(f1*f2)
       ELSE IF(irho==15) THEN
          ! Isothermal beta model with general beta
          !beta=0.86 in Ma et al. (2015)
          beta=p1
          !WRITE(*,*) 'Beta:', beta
          rho=(1.+(r/rs)**2)**(-3.*beta/2.)
       ELSE IF(irho==16) THEN
          ! Isothermal exterior
          rho=1./r**2
       ELSE IF(irho==17) THEN
          ! Power-law profile
          beta=p1
          rho=r**beta
       ELSE IF(irho==18) THEN
          ! Cubic profile
          rho=r**(-3)
       ELSE IF(irho==19) THEN
          ! Smooth profile
          rho=0.
       ELSE IF(irho==20) THEN
          ! Exponential profile (HI from Padmanabhan et al. 2017)
          re=p1
          rho=exp(-r/re)
       ELSE IF(irho==24) THEN
          ! Cored NFW (Copeland, Taylor & Hall 2018)
          rb=p1
          f1=(r+rb)/rs
          f2=(1.+r/rs)**2
          rho=1./(f1*f2)
       ELSE IF(irho==25) THEN
          ! polynomial with central exponential hole
          r0=p1
          alpha=p2
          rho=(r**(-alpha))*exp(-r0/r)
       ELSE IF(irho==26) THEN
          ! modified NFW with central exponential hole
          r0=p1
          rh=p2
          rho=(1./((0.75+r/rh)*(1.+r/rh)**2))*exp(-r0/r)
       ELSE IF(irho==27) THEN
          ! modified NFW from Padmanabhan & Refregier (2017; 1607.01021)
          rh=p1
          rho=(1./((0.75+r/rh)*(1.+r/rh)**2))
       ELSE IF(irho==28) THEN
          ! Shell
          rho=0.
       ELSE
          STOP 'RHO: Error, irho not specified correctly'
       END IF

    END IF

  END FUNCTION rho

  REAL FUNCTION rhor2at0(irho)

    ! This is the value of rho(r)*r^2 at r=0
    ! For most profiles this is zero, BUT not if rho(r->0) -> r^-2
    ! Note if rho(r->0) -> r^n with n<-2 then the profile mass would diverge!

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: irho

    IF(irho==0) THEN
       STOP 'RHOR2AT0: You should not be here for a delta-function profile'
    ELSE IF(irho==1 .OR. irho==9) THEN
       !1 - Isothermal
       !9 - Stellar profile from Schneider & Teyssier (2015)
       rhor2at0=1.
    ELSE IF(irho==18) THEN
       STOP 'RHOR2AT0: Error, profile diverges at the origin'
    ELSE
       rhor2at0=0.
    END IF

  END FUNCTION rhor2at0

  FUNCTION normalisation(rmin,rmax,rv,rs,p1,p2,irho)

    ! This calculates the normalisation of a halo
    ! This is the integral of 4pir^2*rho(r)*dr between rmin and rmax

    ! Profile results
    !  0 - Delta function (M = 1)
    !  1 - Isothermal (M = 4pi*rv)
    !  2 - Top hat (M = (4pi/3)*rv^3)
    !  3 - Moore (M = (8pi/3)*rv^3*ln(1+c^1.5)/c^3)
    !  4 - NFW (M = 4pi*rs^3*[ln(1+c)-c/(1+c)])
    !  5 - NFW (M = 4pi*rs^3*[ln(1+c)-c/(1+c)])
    !  6 - Beta model with beta=2/3 (M = 4pi*rs^3*(rv/rs-atan(rv/rs)))
    !  7 - Fedeli stellar model (M = 4pi*rstar^2 * [1-exp(-rmax/rstar)*(1.+rmax/rstar)]
    !  8 - No
    !  9 - Stellar profile (Schneider & Teyssier 2015; M = 4(pi^1.5)*rstar; assumed rmax ~ infinity)
    ! 10 - Ejected gas profile (Schneider & Teyssier 2015; M = 4pi*sqrt(pi/2)*re^3; assumed rmax ~ infinity)
    ! 11 - No
    ! 12 - No
    ! 13 - No
    ! 14 - No
    ! 15 - No
    ! 16 - Isothermal shell (M = 4pi*(rmax-rmin))
    ! 17 - No
    ! 18 - Cubic profile (M = 4pi*log(rmax/rmin))
    ! 19 - Smooth profile (M = 1; prevents problems)
    ! 20 - No
    ! 21 - No
    ! 22 - No
    ! 23 - No
    ! 24 - No
    ! 25 - No
    ! 26 - No
    ! 27 - No
    ! 28 - Shell (M = 4pi*rv^3)

    IMPLICIT NONE
    REAL :: normalisation
    REAL, INTENT(IN) :: rmin, rmax, rv, rs, p1, p2
    INTEGER, INTENT(IN) :: irho
    REAL :: cmax, re, rstar, beta, rb, c, b

    IF(irho==0) THEN
       ! Delta function
       normalisation=1.
    ELSE IF(irho==1 .OR. irho==16) THEN
       ! Isothermal
       normalisation=4.*pi*(rmax-rmin)
    ELSE IF(irho==2) THEN
       ! Top hat
       normalisation=4.*pi*(rmax**3-rmin**3)/3.
    ELSE IF(irho==3) THEN
       ! Moore et al. (1999)
       IF(rmin .NE. 0.) THEN
          normalisation=winint(0.,rmin,rmax,rv,rs,p1,p2,irho,imeth_win)
       ELSE
          cmax=rmax/rs
          normalisation=(2./3.)*4.*pi*(rs**3)*log(1.+cmax**1.5)
       END IF
    ELSE IF(irho==4 .OR. irho==5) THEN
       ! NFW (1997)
       IF(rmin .NE. 0.) THEN
          normalisation=winint(0.,rmin,rmax,rv,rs,p1,p2,irho,imeth_win)
       ELSE
          cmax=rmax/rs
          normalisation=4.*pi*(rs**3)*NFW_factor(cmax)!(log(1.+cmax)-cmax/(1.+cmax))
       END IF
    ELSE IF(irho==6) THEN
       ! Isothermal beta model with beta=2/3
       IF(rmin .NE. 0.) THEN
          normalisation=winint(0.,rmin,rmax,rv,rs,p1,p2,irho,imeth_win)
       ELSE
          cmax=rmax/rs
          normalisation=4.*pi*(rs**3)*(cmax-atan(cmax))
       END IF
    ELSE IF(irho==7) THEN
       ! Fedeli (2014) stellar model
       IF(rmin .NE. 0) THEN
          ! I could actually derive an analytical expression here if this was ever necessary
          normalisation=winint(0.,rmin,rmax,rv,rs,p1,p2,irho,imeth_win)
       ELSE
          ! This would be even easier if rmax -> infinity (just 4*pi*rstar^2)
          rstar=p1
          normalisation=4.*pi*(rstar**3)*(1.-exp(-rmax/rstar)*(1.+rmax/rstar))
          !normalisation=4.*pi*rstar**3 ! rmax/rstar -> infinity limit (rmax >> rstar)
       END IF
    ELSE IF(irho==9) THEN
       ! Stellar profile from Schneider & Teyssier (2015)       
       IF(rmin .NE. 0.) THEN
          normalisation=winint(0.,rmin,rmax,rv,rs,p1,p2,irho,imeth_win)
       ELSE
          ! Assumed to go on to r -> infinity
          rstar=p1
          normalisation=4.*(pi**(3./2.))*rstar
       END IF
    ELSE IF(irho==10) THEN
       ! Ejected gas profile from Schneider & Teyssier (2015)
       ! Assumed to go on to r -> infinity
       IF(rmin .NE. 0.) THEN
          normalisation=winint(0.,rmin,rmax,rv,rs,p1,p2,irho,imeth_win)
       ELSE
          ! Assumed to go on to r -> infinity
          re=p1
          normalisation=4.*pi*sqrt(pi/2.)*re**3
       END IF
    ELSE IF(irho==17) THEN
       ! Power-law profile
       beta=p1
       normalisation=(4.*pi/(beta+3.))*(rmax**(beta+3.)-rmin**(beta+3.))
    ELSE IF(irho==18) THEN
       ! Cubic profile
       normalisation=4.*pi*log(rmax/rmin)
    ELSE IF(irho==19) THEN
       ! Smooth profile, needs a normalisation I think
       normalisation=1.
    ELSE IF(irho==24) THEN
       ! Cored NFW profile
       rb=p1
       IF(rb==0.) THEN
          ! This is then the standard NFW case
          c=rv/rs
          normalisation=4.*pi*(rs**3)*NFW_factor(c)
       ELSE
          ! Otherwise there is actually a core
          b=rv/rb
          c=rv/rs
          normalisation=(4.*pi*rs**3)/(b-c)**2
          normalisation=normalisation*(b*(b-2.*c)*NFW_factor(c)+(log(1.+b)-b/(1.+c))*c**2)
       END IF
    ELSE IF(irho==28) THEN
       ! Shell
       normalisation=4.*pi*rv**3
    ELSE
       ! Otherwise need to do the integral numerically
       ! k=0 gives normalisation
       normalisation=winint(0.,rmin,rmax,rv,rs,p1,p2,irho,imeth_win)
    END IF

  END FUNCTION normalisation

  REAL FUNCTION win_norm(k,rmin,rmax,rv,rs,p1,p2,irho)

    ! This is an UNNORMALISED halo profile of any sort

    ! Types of profile
    ! ================
    !  0 - Delta function at r=0
    !  1 - Isothermal: r^-2
    !  2 - Top hat: constant
    !  3 - No
    !  4 - No
    !  5 - NFW
    !  6 - No
    !  7 - Star profile
    !  8 - No
    !  9 - Stellar profile from Schneider & Teyssier (2015)
    ! 10 - Ejected gas profile (Schneider & Teyssier 2015)
    ! 11 - No
    ! 12 - No
    ! 13 - No
    ! 14 - No
    ! 15 - No
    ! 16 - Isothermal shell
    ! 17 - No
    ! 18 - No
    ! 19 - Smooth profile
    ! 20 - Exponential profile
    ! 21 - No
    ! 22 - No
    ! 23 - No
    ! 24 - Cored NFW profile (Copeland, Taylor & Hall 2018)
    ! 25 - No
    ! 26 - No
    ! 27 - No
    ! 28 - Shell

    ! Calculates the normalised spherical Fourier Transform of the density profile
    ! Note that this means win_norm(k->0)=1
    ! and that win must be between 0 and 1
    IMPLICIT NONE
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: rmin
    REAL, INTENT(IN) :: rmax
    REAL, INTENT(IN) :: rv
    REAL, INTENT(IN) :: rs
    REAL, INTENT(IN) :: p1
    REAL, INTENT(IN) :: p2
    INTEGER, INTENT(IN) :: irho
    REAL :: re, f1, f2, rstar, kstar, rb

    IF(k==0.) THEN

       ! If called for the zero mode (e.g. for the normalisation)
       win_norm=1.

    ELSE

       IF(irho==0) THEN
          ! Delta function profile is not localised in Fourier Space
          win_norm=1.
       ELSE IF(irho==1) THEN
          win_norm=wk_isothermal(k*rmax)
       ELSE IF(irho==2) THEN
          ! Analytic for top hat
          win_norm=wk_tophat(k*rmax)
       ELSE IF(irho==5) THEN
          ! Analytic for NFW
          win_norm=win_NFW(k,rmax,rs)
       ELSE IF(irho==7) THEN
          ! Analytic for Fedeli (2014) stellar profile
          rstar=p1
          kstar=k*rstar
          f1=kstar-exp(-rmax/rstar)*(sin(k*rmax)+kstar*cos(k*rmax))
          f2=kstar*(1.+kstar**2)
          win_norm=f1/f2
          !win_norm=1./(1.+kstar**2) !bigRstar -> infinity limit (rmax >> rstar)
       ELSE IF(irho==9) THEN
          ! Only valid if rmin=0 and rmax=inf
          rstar=p1
          win_norm=(sqrt(pi)/2.)*erf(k*rstar)/(k*rstar)
       ELSE IF(irho==10) THEN
          ! Ejected gas profile
          re=p1
          win_norm=exp(-1.5*(k*re)**2.)
       ELSE IF(irho==16) THEN
          ! Isothermal shells
          win_norm=wk_isothermal_2(k*rmax,k*rmin)
       ELSE IF(irho==19) THEN
          ! Smooth profile
          win_norm=0.
       ELSE IF(irho==20) THEN
          ! Exponential profile
          re=p1
          win_norm=1./(1.+(k*re)**2)**2
       ELSE IF(irho==24) THEN
          ! Cored NFW profile
          rb=p1
          IF(rb==0.) THEN
             ! In this case there is no core
             win_norm=win_NFW(k,rmax,rs)
          ELSE
             ! Otherwise there is a core
             win_norm=win_cored_NFW(k,rmax,rs,rb)
          END IF
       ELSE IF(irho==28) THEN
          ! Shell
          win_norm=sinc(k*rv)
       ELSE
          ! Numerical integral over the density profile (slower)
          win_norm=winint(k,rmin,rmax,rv,rs,p1,p2,irho,imeth_win)/normalisation(rmin,rmax,rv,rs,p1,p2,irho)
       END IF

    END IF

  END FUNCTION win_norm

  SUBROUTINE winint_speed_tests(k,nk,rmin,rmax,rv,rs,p1,p2,irho)

    IMPLICIT NONE
    REAL, INTENT(IN) :: k(nk), rmin, rmax, rv, rs, p1, p2
    INTEGER, INTENT(IN) :: nk, irho
    CHARACTER(len=256) :: base, ext, outfile
    INTEGER :: j, i, ii, imeth, n, ntime
    LOGICAL :: timing
    REAL :: t1, t2, w

    base='winint/results_'
    ext='.dat'

    ! j = 1: Time calculation
    ! j = 2: Write calculation out
    DO j=1,2

       timing=.FALSE.
       IF(j==2) THEN
          timing=.TRUE.
          WRITE(*,*) 'WININT_SPEED_TESTS: Doing this many evaluations for timing test:', ntime
          WRITE(*,*)
       END IF

       DO imeth=1,nmeth_win

          WRITE(*,*) 'WININT_SPEED_TESTS: Method:', imeth
          IF(imeth==1)  WRITE(*,*) 'WININT_SPEED_TESTS: Simple integration'
          IF(imeth==2)  WRITE(*,*) 'WININT_SPEED_TESTS: Simple integration over bumps'
          IF(imeth==3)  WRITE(*,*) 'WININT_SPEED_TESTS: Standard intergation'
          IF(imeth==4)  WRITE(*,*) 'WININT_SPEED_TESTS: Standard integration over bumps'
          IF(imeth==5)  WRITE(*,*) 'WININT_SPEED_TESTS: Constant approximation for bumps'
          IF(imeth==6)  WRITE(*,*) 'WININT_SPEED_TESTS: Linear approximation for bumps'
          IF(imeth==7)  WRITE(*,*) 'WININT_SPEED_TESTS: Quadratic approximation for bumps'
          IF(imeth==8)  WRITE(*,*) 'WININT_SPEED_TESTS: Cubic approximation for bumps'
          IF(imeth==9)  WRITE(*,*) 'WININT_SPEED_TESTS: Hybrid with constant approximation for bumps'
          IF(imeth==10) WRITE(*,*) 'WININT_SPEED_TESTS: Hybrid with linear approximation for bumps'
          IF(imeth==11) WRITE(*,*) 'WININT_SPEED_TESTS: Hybrid with quadratic approximation for bumps'
          IF(imeth==12) WRITE(*,*) 'WININT_SPEED_TESTS: Hybrid with cubic approximation for bumps'
          IF(imeth==13) WRITE(*,*) 'WININT_SPEED_TESTS: Full hybrid'

          IF(timing .EQV. .FALSE.) THEN
             outfile=number_file(base,imeth,ext)
             WRITE(*,*) 'WININT_SPEED_TESTS: Writing data: ', TRIM(outfile)
             OPEN(7,file=outfile)            
             n=1
             IF(imeth==1) CALL cpu_time(t1)
          ELSE
             CALL cpu_time(t1)
             n=ntime
          END IF

          ! Loop over number of iterations
          DO ii=1,n

             ! Loop over wave number and do integration
             DO i=1,nk
                w=winint(k(i),rmin,rmax,rv,rs,p1,p2,irho,imeth)/normalisation(rmin,rmax,rv,rs,p1,p2,irho)
                IF(.NOT. timing) THEN
                   WRITE(7,*) k(i), w
                END IF
             END DO

          END DO

          IF(.NOT. timing) THEN
             CLOSE(7)
             IF(imeth==1) THEN
                CALL cpu_time(t2)
                ntime=CEILING(winint_test_seconds/(t2-t1))
             END IF
          ELSE
             CALL cpu_time(t2)
             WRITE(*,fmt='(A30,F10.5)') 'WININT_SPEED_TESTS: Time [s]:', t2-t1
          END IF

          WRITE(*,*)

       END DO

    END DO

  END SUBROUTINE winint_speed_tests

  SUBROUTINE winint_diagnostics(rmin,rmax,rv,rs,p1,p2,irho,outfile)

    ! Write out the winint integrand as a function of k
    IMPLICIT NONE
    REAL, INTENT(IN) :: rmin, rmax, rv, rs, p1, p2
    INTEGER, INTENT(IN) :: irho
    CHARACTER(len=256), INTENT(IN) :: outfile
    INTEGER :: i, j
    REAL :: r, k
    REAL, ALLOCATABLE :: integrand(:)

    REAL, PARAMETER :: kmin=1e-1
    REAL, PARAMETER :: kmax=1e2
    INTEGER, PARAMETER :: nr=256 ! Number of points in r
    INTEGER, PARAMETER :: nk=16  ! Number of points in k

    WRITE(*,*) 'WININT_DIAGNOSTICS: Doing these'
    WRITE(*,*) 'WININT_DIAGNOSTICS: minimum r [Mpc/h]:', REAL(rmin)
    WRITE(*,*) 'WININT_DIAGNOSTICS: maximum r [Mpc/h]:', REAL(rmax)
    WRITE(*,*) 'WININT_DIAGNOSTICS: virial radius [Mpc/h]:', REAL(rv)
    WRITE(*,*) 'WININT_DIAGNOSTICS: scale radius [Mpc/h]:', REAL(rs)
    WRITE(*,*) 'WININT_DIAGNOSTICS: concentration:', REAL(rv/rs)
    WRITE(*,*) 'WININT_DIAGNOSTICS: halo parameter 1:', p1
    WRITE(*,*) 'WININT_DIAGNOSTICS: halo parameter 2:', p2
    WRITE(*,*) 'WININT_DIAGNOSTICS: profile number:', irho
    WRITE(*,*) 'WININT_DIAGNOSTICS: outfile: ', TRIM(outfile)

    ALLOCATE(integrand(nk))

    OPEN(7,file=outfile)
    DO i=1,nr
       r=progression(0.,rmax,i,nr)
       DO j=1,nk
          k=exp(progression(log(kmin),log(kmax),j,nk))
          integrand(j)=winint_integrand(r,rmin,rmax,rv,rs,p1,p2,irho)*sinc(r*k)
       END DO
       WRITE(7,*) r/rv, (integrand(j), j=1,nk)
    END DO
    CLOSE(7)

    WRITE(*,*) 'WININT_DIAGNOSTICS: Done'
    WRITE(*,*)

  END SUBROUTINE winint_diagnostics

  FUNCTION winint(k,rmin,rmax,rv,rs,p1,p2,irho,imeth)

    ! Calculates W(k,M)
    IMPLICIT NONE
    REAL :: winint
    REAL, INTENT(IN) :: k, rmin, rmax, rv, rs, p1, p2
    INTEGER, INTENT(IN) :: irho, imeth

    ! Integration method
    ! imeth =  1: Simple integration
    ! imeth =  2: Bumps with simple integration
    ! imeth =  3: Standard integration
    ! imeth =  4: Bumps with standard integration
    ! imeth =  5: Constant approximation for bumps
    ! imeth =  6: Linear approximation for bumps
    ! imeth =  7: Quadratic approximation for bumps
    ! imeth =  8: Cubic approximation for bumps
    ! imeth =  9: Hybrid with constant approximation for bumps
    ! imeth = 10: Hybrid with linear approximation for bumps
    ! imeth = 11: Hybrid with quadratic approximation for bumps
    ! imeth = 12: Hybrid with cubic approximation for bumps

    ! Bump methods go crazy with some star profiles (those that drop too fast)
    ! You need to make sure that the rmax for the integration does not extend too far out

    IF(imeth==1) THEN
       winint=winint_normal(rmin,rmax,k,rmin,rmax,rv,rs,p1,p2,irho,winint_order,acc_win)
    ELSE IF(imeth==3) THEN
       winint=winint_store(rmin,rmax,k,rmin,rmax,rv,rs,p1,p2,irho,winint_order,acc_win)
    ELSE
       winint=winint_bumps(k,rmin,rmax,rv,rs,p1,p2,irho,winint_order,acc_win,imeth)
    END IF

  END FUNCTION winint

  FUNCTION winint_normal(a,b,k,rmin,rmax,rv,rs,p1,p2,irho,iorder,acc)

    ! Integration routine using 'normal' method to calculate the normalised halo FT
    IMPLICIT NONE
    REAL :: winint_normal
    REAL, INTENT(IN) :: k, rmin, rmax, rv, rs, p1, p2
    INTEGER, INTENT(IN) :: irho
    INTEGER, INTENT(IN) :: iorder
    REAL, INTENT(IN) :: acc
    REAL, INTENT(IN) :: a, b
    DOUBLE PRECISION :: sum
    REAL :: r, dr, winold, weight
    INTEGER :: n, i, j

    INTEGER, PARAMETER :: jmin=5
    INTEGER, PARAMETER :: jmax=30
    INTEGER, PARAMETER :: ninit=2

    winold=0.

    IF(a==b) THEN

       winint_normal=0.

    ELSE

       !Integrates to required accuracy!
       DO j=1,jmax

          !Increase the number of integration points each go until convergence
          n=ninit*(2**(j-1))

          !Set the integration sum variable to zero
          sum=0.

          DO i=1,n

             !Get the weights
             IF(iorder==1) THEN
                !Composite trapezium weights
                IF(i==1 .OR. i==n) THEN
                   weight=0.5
                ELSE
                   weight=1.
                END IF
             ELSE IF(iorder==2) THEN
                !Composite extended formula weights
                IF(i==1 .OR. i==n) THEN
                   weight=0.416666666666
                ELSE IF(i==2 .OR. i==n-1) THEN
                   weight=1.083333333333
                ELSE
                   weight=1.
                END IF
             ELSE IF(iorder==3) THEN
                !Composite Simpson weights
                IF(i==1 .OR. i==n) THEN
                   weight=0.375
                ELSE IF(i==2 .OR. i==n-1) THEN
                   weight=1.166666666666
                ELSE IF(i==3 .OR. i==n-2) THEN
                   weight=0.958333333333
                ELSE
                   weight=1.
                END IF
             ELSE
                STOP 'WININT_NORMAL: Error, order specified incorrectly'
             END IF

             !Now get r and do the function evaluations
             r=progression(a,b,i,n)
             sum=sum+weight*winint_integrand(r,rmin,rmax,rv,rs,p1,p2,irho)*sinc(r*k)

          END DO

          !The dr are all equally spaced
          dr=(b-a)/REAL(n-1)

          winint_normal=REAL(sum)*dr

          IF((j>jmin) .AND. winint_normal==0.) THEN
             EXIT
          ELSE IF((j>jmin) .AND. (ABS(-1.+winint_normal/winold)<acc)) THEN
             EXIT
          ELSE
             winold=winint_normal
          END IF

       END DO

    END IF

  END FUNCTION winint_normal

  FUNCTION winint_store(a,b,k,rmin,rmax,rv,rs,p1,p2,irho,iorder,acc)

    !Integrates between a and b until desired accuracy is reached
    !Stores information to reduce function calls
    IMPLICIT NONE
    REAL :: winint_store
    REAL, INTENT(IN) :: k, rmin, rmax, rv, rs, p1, p2
    REAL, INTENT(IN) :: acc
    INTEGER, INTENT(IN) :: iorder, irho
    REAL, INTENT(IN) :: a, b
    INTEGER :: i, j, n
    REAL :: x, dx
    REAL :: f1, f2, fx
    DOUBLE PRECISION :: sum_n, sum_2n, sum_new, sum_old

    INTEGER, PARAMETER :: jmin=5
    INTEGER, PARAMETER :: jmax=30

    IF(a==b) THEN

       !Fix the answer to zero if the integration limits are identical
       winint_store=0.

    ELSE

       !Reset the sum variables for the integration
       sum_2n=0.d0
       sum_n=0.d0
       sum_old=0.d0
       sum_new=0.d0

       DO j=1,jmax

          !Note, you need this to be 1+2**n for some integer n
          !j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
          n=1+2**(j-1)

          !Calculate the dx interval for this value of 'n'
          dx=(b-a)/REAL(n-1)

          IF(j==1) THEN

             !The first go is just the trapezium of the end points
             f1=winint_integrand(a,rmin,rmax,rv,rs,p1,p2,irho)*sinc(a*k)
             f2=winint_integrand(b,rmin,rmax,rv,rs,p1,p2,irho)*sinc(b*k)
             sum_2n=0.5*(f1+f2)*dx
             sum_new=sum_2n

          ELSE

             !Loop over only new even points to add these to the integral
             DO i=2,n,2
                x=progression(a,b,i,n)
                fx=winint_integrand(x,rmin,rmax,rv,rs,p1,p2,irho)*sinc(x*k)
                sum_2n=sum_2n+fx
             END DO

             !Now create the total using the old and new parts
             sum_2n=sum_n/2.+sum_2n*dx

             !Now calculate the new sum depending on the integration order
             IF(iorder==1) THEN  
                sum_new=sum_2n
             ELSE IF(iorder==3) THEN         
                sum_new=(4.*sum_2n-sum_n)/3. !This is Simpson's rule and cancels error
             ELSE
                STOP 'WININT_STORE: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (ABS(-1.+sum_new/sum_old)<acc)) THEN
             !jmin avoids spurious early convergence
             winint_store=REAL(sum_new)
             EXIT
          ELSE IF(j==jmax) THEN
             winint_store=0.d0
             STOP 'WININT_STORE: Integration timed out'
          ELSE
             !Integral has not converged so store old sums and reset sum variables
             winint_store=0.d0
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.
          END IF

       END DO

    END IF

  END FUNCTION winint_store

  FUNCTION winint_bumps(k,rmin,rmax,rv,rs,p1,p2,irho,iorder,acc,imeth)

    ! Integration routine to calculate the normalised halo FT
    IMPLICIT NONE
    REAL :: winint_bumps
    REAL, INTENT(IN) :: k, rmin, rmax, rv, rs, p1, p2
    INTEGER, INTENT(IN) :: irho
    INTEGER, INTENT(IN) :: iorder, imeth
    REAL, INTENT(IN) :: acc
    REAL :: sum, w
    REAL :: r1, r2
    INTEGER :: i, n

    ! This MUST be set to zero for this routine
    IF(rmin .NE. 0.) STOP 'WININT_BUMPS: Error, rmin must be zero'

    ! Calculate the number of nodes of sinc(k*rmax) for 0<=r<=rmax
    n=FLOOR(k*rmax/pi)

    ! Set the sum variable to zero
    sum=0.

    ! Integrate over each chunk between nodes separately
    DO i=0,n

       ! Set the lower integration limit
       IF(k==0.) THEN
          ! Special case when k=0 to avoid division by zero
          r1=0.
       ELSE
          r1=i*pi/k
       END IF

       ! Set the upper integration limit
       IF(k==0. .OR. i==n) THEN
          ! Special case when on last section because end is rmax, not a node!
          r2=rmax
       ELSE
          r2=(i+1)*pi/k
       END IF

       ! Now do the integration along a section
       IF(imeth==2) THEN
          w=winint_normal(r1,r2,k,rmin,rmax,rv,rs,p1,p2,irho,iorder,acc)
       ELSE IF(i==0 .OR. i==n .OR. k==0. .OR. imeth==4) THEN
          w=winint_store(r1,r2,k,rmin,rmax,rv,rs,p1,p2,irho,iorder,acc)
       ELSE IF(imeth==13) THEN
          w=winint_hybrid(r1,r2,i,k,rmin,rmax,rv,rs,p1,p2,irho,iorder,acc)
       ELSE
          IF(imeth==5 .OR. (imeth==9 .AND. n>nlim_bumps)) THEN
             w=winint_approx(r1,r2,i,k,rmin,rmax,rv,rs,p1,p2,irho,iorder=0)
          ELSE IF(imeth==6 .OR. (imeth==10 .AND. n>nlim_bumps)) THEN
             w=winint_approx(r1,r2,i,k,rmin,rmax,rv,rs,p1,p2,irho,iorder=1)
          ELSE IF(imeth==7 .OR. (imeth==11 .AND. n>nlim_bumps)) THEN
             w=winint_approx(r1,r2,i,k,rmin,rmax,rv,rs,p1,p2,irho,iorder=2)
          ELSE IF(imeth==8 .OR. (imeth==12 .AND. n>nlim_bumps)) THEN
             w=winint_approx(r1,r2,i,k,rmin,rmax,rv,rs,p1,p2,irho,iorder=3)
          ELSE
             w=winint_store(r1,r2,k,rmin,rmax,rv,rs,p1,p2,irho,iorder,acc)
          END IF
       END IF

       sum=sum+w

       ! Exit if the contribution to the sum is very tiny, this seems to be necessary to prevent crashes
       IF(winint_exit .AND. ABS(w)<acc*ABS(sum)) EXIT

    END DO

    winint_bumps=sum

  END FUNCTION winint_bumps

  REAL FUNCTION winint_approx(rn,rm,i,k,rmin,rmax,rv,rs,p1,p2,irho,iorder)

    ! Approximate forms for the integral over the sine bump times a polynomial
    USE fix_polynomial
    IMPLICIT NONE
    REAL, INTENT(IN) :: rn, rm
    INTEGER, INTENT(IN) :: i
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: rmin, rmax
    REAL, INTENT(IN) :: rv, rs, p1, p2
    INTEGER, INTENT(IN) :: irho
    INTEGER, INTENT(IN) :: iorder
    REAL :: x0, x1, x2, x3
    REAL :: y0, y1, y2, y3
    REAL :: a0, a1, a2, a3
    REAL :: w

    IF(iorder==0) THEN
       x0=(rn+rm)/2. ! Middle
    ELSE IF(iorder==1) THEN
       x0=rn ! Beginning
       x1=rm ! End
    ELSE IF(iorder==2) THEN
       x0=rn ! Beginning
       x1=(rn+rm)/2. ! Middle
       x2=rm ! End
    ELSE IF(iorder==3) THEN
       x0=rn ! Beginning
       x1=rn+1.*(rm-rn)/3. ! Middle
       x2=rn+2.*(rm-rn)/3. ! Middle
       x3=rm ! End
    ELSE
       STOP 'WININT_APPROX: Error, iorder specified incorrectly'
    END IF

    y0=winint_integrand(x0,rmin,rmax,rv,rs,p1,p2,irho)/x0
    IF(iorder==1 .OR. iorder==2 .OR. iorder==3) y1=winint_integrand(x1,rmin,rmax,rv,rs,p1,p2,irho)/x1
    IF(iorder==2 .OR. iorder==3) y2=winint_integrand(x2,rmin,rmax,rv,rs,p1,p2,irho)/x2
    IF(iorder==3) y3=winint_integrand(x3,rmin,rmax,rv,rs,p1,p2,irho)/x3

    IF(iorder==0) THEN
       a0=y0
    ELSE IF(iorder==1) THEN
       CALL fix_line(a1,a0,x0,y0,x1,y1)
    ELSE IF(iorder==2) THEN
       CALL fix_quadratic(a2,a1,a0,x0,y0,x1,y1,x2,y2)
    ELSE IF(iorder==3) THEN
       CALL fix_cubic(a3,a2,a1,a0,x0,y0,x1,y1,x2,y2,x3,y3)
    ELSE
       STOP 'WININT_APPROX: Error, iorder specified incorrectly'
    END IF

    w=2.*a0
    IF(iorder==1 .OR. iorder==2 .OR. iorder==3) w=w+(rm+rn)*a1
    IF(iorder==2 .OR. iorder==3) w=w+(rm**2+rn**2-4./k**2)*a2
    IF(iorder==3) w=w+(rm**3+rn**3-6.*(rm+rn)/k**2)*a3

    w=((-1)**i)*w/k**2

    winint_approx=w

  END FUNCTION winint_approx

  REAL FUNCTION winint_hybrid(rn,rm,i,k,rmin,rmax,rv,rs,p1,p2,irho,iorder,acc)

    ! An attempt to do an automatic combination of approximation and proper integration
    ! It turned out to be very slow
    USE fix_polynomial
    IMPLICIT NONE
    REAL, INTENT(IN) :: rn, rm
    INTEGER, INTENT(IN) :: i
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: rmin, rmax
    REAL, INTENT(IN) :: rv, rs, p1, p2
    INTEGER, INTENT(IN) :: irho
    INTEGER, INTENT(IN) :: iorder
    REAL, INTENT(IN) :: acc
    REAL :: x0, x1, x2, x3
    REAL :: y0, y1, y2, y3
    REAL :: a0, a1, a2, a3
    REAL :: w, rmid, epsa0

    REAL, PARAMETER :: eps=1e-2

    rmid=(rn+rm)/2.

    x0=rn
    x1=rn+1.*(rm-rn)/3.
    x2=rn+2.*(rm-rn)/3.
    x3=rm

    y0=winint_integrand(x0,rmin,rmax,rv,rs,p1,p2,irho)/x0
    y1=winint_integrand(x1,rmin,rmax,rv,rs,p1,p2,irho)/x1
    y2=winint_integrand(x2,rmin,rmax,rv,rs,p1,p2,irho)/x2
    y3=winint_integrand(x3,rmin,rmax,rv,rs,p1,p2,irho)/x3

    CALL fix_cubic(a3,a2,a1,a0,x0,y0,x1,y1,x2,y2,x3,y3)

    epsa0=eps*a0
    IF(ABS(a3*rmid**3)<epsa0 .AND. ABS(a2*rmid**2)<epsa0 .AND. ABS(a1*rmid)<epsa0) THEN
       w=(rm**3+rn**3-6.*(rm+rn)/k**2)*a3+(rm**2+rn**2-4./k**2)*a2+(rm+rn)*a1+2.*a0
       w=w*(-1)**i
    ELSE
       w=winint_store(rn,rm,k,rmin,rmax,rv,rs,p1,p2,irho,iorder,acc)
    END IF

    winint_hybrid=w

  END FUNCTION winint_hybrid

  FUNCTION winint_integrand(r,rmin,rmax,rv,rs,p1,p2,irho)

    !The integrand for the W(k) integral
    !Note that the sinc function is *not* included
    IMPLICIT NONE
    REAL :: winint_integrand
    REAL, INTENT(IN) :: r, rmin, rmax, rv, rs, p1, p2
    INTEGER, INTENT(IN) :: irho

    IF(r==0.) THEN
       winint_integrand=4.*pi*rhor2at0(irho)
    ELSE
       winint_integrand=4.*pi*(r**2)*rho(r,rmin,rmax,rv,rs,p1,p2,irho)
    END IF

  END FUNCTION winint_integrand

  REAL FUNCTION win_NFW(k,rv,rs)

    ! The analytic normalised (W(k=0)=1) Fourier Transform of the NFW profile
    IMPLICIT NONE
    REAL, INTENT(IN) :: k, rv, rs
    REAL :: c, ks
    REAL :: p1, p2, p3
    REAL :: rmin, rmax

    c=rv/rs
    ks=k*rv/c

    p1=cos(ks)*F_NFW(k,rv,c)
    p2=sin(ks)*G_NFW(k,rv,c)
    p3=sin(ks*c)/(ks*(1.+c))

    rmin=0.
    rmax=rv
    win_NFW=4.*pi*(rs**3)*(p1+p2-p3)/normalisation(rmin,rmax,rv,rs,zero,zero,4)

  END FUNCTION win_NFW

  REAL FUNCTION win_cored_NFW(k,rv,rs,rb)

    ! The analytic normalised (W(k=0)=1) Fourier Transform of the cored NFW profile
    ! Appendix A of Copeland, Taylor & Hall (1712.07112)
    IMPLICIT NONE
    REAL, INTENT(IN) :: k, rv, rs, rb
    REAL :: b, c
    REAL :: rmin, rmax

    b=rv/rb
    c=rv/rs

    rmin=0.
    rmax=rv

    win_cored_NFW=NFW_factor(c)*win_NFW(k,rv,rs)
    win_cored_NFW=win_cored_NFW+(c/(b-c))*(1./(k*rs))*(G_NFW(k,rv,c)*cos(k*rs)-F_NFW(k,rv,c)*sin(k*rs))
    win_cored_NFW=win_cored_NFW-(c/(b-c))*(1./(k*rs))*(G_NFW(k,rv,b)*cos(k*rb)-F_NFW(k,rv,b)*sin(k*rb))
    win_cored_NFW=4.*pi*(rs**3)*(b/(b-c))*win_cored_NFW/normalisation(rmin,rmax,rv,rs,rb,zero,24)

  END FUNCTION win_cored_NFW

  REAL FUNCTION F_NFW(k,rv,c)

    ! Equation (A3; top) from Copeland, Taylor & Hall (1712.07112)
    IMPLICIT NONE
    REAL, INTENT(IN) :: k, rv, c
    REAL :: ks

    ks=k*rv/c

    F_NFW=Ci(ks*(1.+c))-Ci(ks)

  END FUNCTION F_NFW

  REAL FUNCTION G_NFW(k,rv,c)

    ! Equation (A3; bottom) from Copeland, Taylor & Hall (1712.07112)
    IMPLICIT NONE
    REAL, INTENT(IN) :: k, rv, c
    REAL :: ks

    ks=k*rv/c

    G_NFW=Si(ks*(1.+c))-Si(ks)

  END FUNCTION G_NFW

  REAL FUNCTION b_nu(nu,hmod)

    ! Bias function selection
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu
    TYPE(halomod), INTENT(INOUT) :: hmod

    IF(.NOT. hmod%has_mass_function) CALL init_mass_function(hmod)

    IF(hmod%imf==1) THEN
       b_nu=b_ps(nu,hmod)
    ELSE IF(hmod%imf==2) THEN
       b_nu=b_st(nu,hmod)
    ELSE IF(hmod%imf==3) THEN
       b_nu=b_Tinker(nu,hmod)
    ELSE IF(hmod%imf==4) THEN
       b_nu=1.
    ELSE
       STOP 'B_NU: Error, imf not specified correctly'
    END IF

  END FUNCTION b_nu

  REAL FUNCTION b_ps(nu,hmod)

    ! Press & Scheter (1974) halo bias
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: dc

    dc=hmod%dc

    b_ps=1.+(nu**2-1.)/dc

  END FUNCTION b_ps

  REAL FUNCTION b_st(nu,hmod)

    ! Sheth & Tormen (1999) halo bias (equation 12 in 9901122)
    ! Comes from peak-background split
    ! Haloes defined with SO relative to mean matter density with SC Delta_v relation
    ! A redshift dependent delta_c is used for barrier height, again from SC
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: dc, p, q

    p=hmod%ST_p
    q=hmod%ST_q
    dc=hmod%dc

    b_st=1.+(q*(nu**2)-1.+2.*p/(1.+(q*nu**2)**p))/dc

  END FUNCTION b_st

  REAL FUNCTION b_Tinker(nu,hmod)

    ! Tinker et al. (2010; 1001.3162) halo bias
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: dc
    REAL :: beta, gamma, phi, eta

    !IF(.NOT. hmod%has_Tinker) CALL init_Tinker(hmod)

    beta=hmod%Tinker_beta
    gamma=hmod%Tinker_gamma
    phi=hmod%Tinker_phi
    eta=hmod%Tinker_eta

    dc=hmod%dc

    b_Tinker=1.+(gamma*nu**2-(1.+2.*eta))/dc+(2.*phi/dc)/(1.+(beta*nu)**(2.*phi))

  END FUNCTION b_Tinker

  REAL FUNCTION b2_nu(nu,hmod)

    ! Bias function selection
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu
    TYPE(halomod), INTENT(INOUT) :: hmod

    IF(hmod%imf==1) THEN
       b2_nu=b2_ps(nu,hmod)
    ELSE IF(hmod%imf==2) THEN
       b2_nu=b2_st(nu,hmod)
    ELSE IF(hmod%imf==3) THEN
       STOP 'B2_NU: Error, second-order bias not specified for Tinker mass function'
    ELSE
       STOP 'B2_NU: Error, imf not specified correctly'
    END IF

  END FUNCTION b2_nu

  REAL FUNCTION b2_ps(nu,hmod)

    ! Press & Schechter (1974) second order bias
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: eps1, eps2, E1, E2, dc

    REAL, PARAMETER :: a2=-17./21.
    REAL, PARAMETER :: p=0.0
    REAL, PARAMETER :: q=1.0

    dc=hmod%dc

    STOP 'B2_PS: Check this very carefully'
    ! I just took the ST form and set p=0 and q=1

    eps1=(q*nu**2-1.)/dc
    eps2=(q*nu**2)*(q*nu**2-3.)/dc**2
    E1=(2.*p)/(dc*(1.+(q*nu**2)**p))
    E2=((1.+2.*p)/dc+2.*eps1)*E1

    b2_ps=2.*(1.+a2)*(eps1+E1)+eps2+E2

  END FUNCTION b2_ps

  REAL FUNCTION b2_st(nu,hmod)

    ! Sheth, Mo & Tormen (2001) second-order bias
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: eps1, eps2, E1, E2, dc

    ! Notation follows from Cooray & Sheth (2002) pp 25-26

    REAL, PARAMETER :: a2=-17./21.
    REAL, PARAMETER :: p=0.3
    REAL, PARAMETER :: q=0.707

    dc=hmod%dc

    eps1=(q*nu**2-1.)/dc
    eps2=(q*nu**2)*(q*nu**2-3.)/dc**2
    E1=(2.*p)/(dc*(1.+(q*nu**2)**p))
    E2=((1.+2.*p)/dc+2.*eps1)*E1

    b2_st=2.*(1.+a2)*(eps1+E1)+eps2+E2

  END FUNCTION b2_st

  REAL FUNCTION g_nu(nu,hmod)

    ! Mass function
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu
    TYPE(halomod), INTENT(INOUT) :: hmod

    IF(.NOT. hmod%has_mass_function) CALL init_mass_function(hmod)

    IF(hmod%imf==1) THEN
       g_nu=g_ps(nu,hmod)
    ELSE IF(hmod%imf==2) THEN
       g_nu=g_st(nu,hmod)
    ELSE IF(hmod%imf==3) THEN
       g_nu=g_Tinker(nu,hmod)
    ELSE
       STOP 'G_NU: Error, imf specified incorrectly'
    END IF

  END FUNCTION g_nu

  REAL FUNCTION g_mu(mu,hmod)

    ! Mass function g(nu) times dnu/dmu = alpha*mu^(alpha-1); nu=mu^alpha
    ! This transformation makes the integral easier at low nu
    ! TODO: Implement a Taylor expansion for low nu/mu
    IMPLICIT NONE
    REAL, INTENT(IN) :: mu
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: nu, alpha

    ! Relation between nu and mu: nu=mu**alpha
    alpha=hmod%alpha_numu
    nu=mu**alpha

    g_mu=g_nu(nu,hmod)*alpha*mu**(alpha-1.)

  END FUNCTION g_mu

  REAL FUNCTION g_ps(nu,hmod)

    ! Press & Scheter (1974) mass function!
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: crap
    REAL, PARAMETER :: A=sqrt(2./pi)

    ! Stop compile-time warnings
    crap=hmod%a

    g_ps=A*exp(-(nu**2)/2.)

  END FUNCTION g_ps

  REAL FUNCTION g_st(nu,hmod)

    ! Sheth & Tormen (1999) mass function, equation (10) in arXiv:9901122
    ! Note I use nu=dc/sigma(M) and this Sheth & Tormen (1999) use nu=(dc/sigma)^2, which accounts for some small differences
    ! Haloes defined with SO relative to mean matter density with SC Delta_v relation
    ! A redshift dependent delta_c is used for barrier height, again from SC
    ! TODO: Implement a Taylor expansion for low nu
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: p, q, A

    ! Mass-function parameters
    A=hmod%ST_A
    p=hmod%ST_p
    q=hmod%ST_q

    g_st=A*(1.+((q*nu**2)**(-p)))*exp(-q*nu**2/2.)

  END FUNCTION g_st

  REAL FUNCTION g_Tinker(nu,hmod)

    ! Tinker et al. (2010; 1001.3162) mass function (also 2008; xxxx.xxxx)
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: alpha, beta, gamma, phi, eta

    !IF(.NOT. hmod%has_Tinker) CALL init_Tinker(hmod)

    alpha=hmod%Tinker_alpha
    beta=hmod%Tinker_beta
    gamma=hmod%Tinker_gamma
    phi=hmod%Tinker_phi
    eta=hmod%Tinker_eta

    ! The actual mass function
    g_Tinker=alpha*(1.+(beta*nu)**(-2.*phi))*nu**(2.*eta)*exp(-0.5*gamma*nu**2)

  END FUNCTION g_Tinker

  SUBROUTINE init_mass_function(hmod)

    ! Initialise anything to do with the halo mass function
    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod

    IF(hmod%imf==1 .OR. hmod%imf==4) THEN
       hmod%has_mass_function=.TRUE.
    ELSE IF(hmod%imf==2) THEN
       CALL init_ST(hmod)
    ELSE IF(hmod%imf==3) THEN
       CALL init_Tinker(hmod)
    ELSE
       STOP 'INIT_MASS_FUNCTION: Error, something went wrong'
    END IF
    
  END SUBROUTINE init_mass_function

  SUBROUTINE init_ST(hmod)

    ! Normalises ST mass function
    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: p, q

    ! ST parameters
    p=hmod%ST_p
    q=hmod%ST_q

    ! Normalisation of ST mass function (involves Gamma function)
    hmod%ST_A=1./(sqrt(pi/(2.*q))+(1./sqrt(q))*(2.**(-p-0.5))*Gamma(0.5-p))

    ! Set the flag to true
    hmod%has_mass_function=.TRUE.
    
  END SUBROUTINE init_ST

  SUBROUTINE init_Tinker(hmod)

    ! Initialise the parameters of the Tinker et al. (2010) mass function and bias
    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: beta, Gamma, phi, eta
    REAL :: z, Dv

    ! Parameter arrays from Tinker (2010): Table 4
    INTEGER, PARAMETER :: n=9 ! Number of entries in parameter lists
    REAL, PARAMETER :: Delta_v(n)=[200.,300.,400.,600.,800.,1200.,1600.,2400.,3200.]
    !REAL, PARAMETER :: alpha0(n)=[0.368,0.363,0.385,0.389,0.393,0.365,0.379,0.355,0.327] ! Not needed if you normalise explicitly
    REAL, PARAMETER :: beta0(n)=[0.589,0.585,0.544,0.543,0.564,0.623,0.637,0.673,0.702]
    REAL, PARAMETER :: gamma0(n)=[0.864,0.922,0.987,1.09,1.20,1.34,1.50,1.68,1.81]
    REAL, PARAMETER :: phi0(n)=[-0.729,-0.789,-0.910,-1.05,-1.20,-1.26,-1.45,-1.50,-1.49]
    REAL, PARAMETER :: eta0(n)=[-0.243,-0.261,-0.261,-0.273,-0.278,-0.301,-0.301,-0.319,-0.336]
    REAL, PARAMETER :: beta_z_exp=0.20
    REAL, PARAMETER :: gamma_z_exp=-0.01
    REAL, PARAMETER :: phi_z_exp=-0.08
    REAL, PARAMETER :: eta_z_exp=0.27

    INTEGER, PARAMETER :: iorder=1 ! Order for interpolation
    INTEGER, PARAMETER :: ifind=3  ! Scheme for finding       (3 - Mid-point splitting)
    INTEGER, PARAMETER :: imeth=2  ! Method for interpolation (2 - Lagrange polynomial)

    LOGICAL, PARAMETER :: z_dependence=.TRUE. ! Do redshift dependence or not

    ! Get these from the halo-model structure
    z=hmod%z
    IF(z>3.) z=3. ! Recommendation from Tinker et al. (2010)
    Dv=hmod%Dv

    ! Delta_v dependence (changed to log Dv finding)
    !alpha=find(log(Dv),log(Delta_v),alpha0,n,iorder,ifind,imeth) ! Not needed if you normalise explicitly
    beta=find(log(Dv),log(Delta_v),beta0,n,iorder,ifind,imeth)
    gamma=find(log(Dv),log(Delta_v),gamma0,n,iorder,ifind,imeth)
    phi=find(log(Dv),log(Delta_v),phi0,n,iorder,ifind,imeth)
    eta=find(log(Dv),log(Delta_v),eta0,n,iorder,ifind,imeth)

    ! Redshift dependence
    IF(z_dependence) THEN
       beta=beta*(1.+z)**beta_z_exp     ! Equation (9)
       gamma=gamma**(1.+z)**gamma_z_exp ! Equation (12)
       phi=phi*(1.+z)**phi_z_exp        ! Equation (10)
       eta=eta*(1.+z)**eta_z_exp        ! Equation (11)
    END IF

    ! Set the Tinker parameters
    hmod%Tinker_alpha=1. ! Fix to unity before normalisation
    hmod%Tinker_beta=beta
    hmod%Tinker_gamma=gamma
    hmod%Tinker_phi=phi
    hmod%Tinker_eta=eta

    ! Set the flag
    hmod%has_mass_function=.TRUE.

    ! Explicitly normalise
    !hmod%Tinker_alpha=hmod%Tinker_alpha/integrate_g_nu(hmod%small_nu,hmod%large_nu,hmod)
    hmod%Tinker_alpha=hmod%Tinker_alpha/integrate_g_mu(hmod%small_nu,hmod%large_nu,hmod)
    
  END SUBROUTINE init_Tinker

  REAL FUNCTION gb_nu(nu,hmod)

    ! g(nu)*b(nu); useful as an integrand
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu               ! Proxy mass variable
    TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model

    gb_nu=g_nu(nu,hmod)*b_nu(nu,hmod)

  END FUNCTION gb_nu

  REAL FUNCTION gb_nu_on_M(nu,hmod)

    ! g(nu)*b(nu)/M(nu); useful as an integrand
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu               ! Proxy mass variable
    TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model

    gb_nu_on_M=g_nu(nu,hmod)*b_nu(nu,hmod)/M_nu(nu,hmod)

  END FUNCTION gb_nu_on_M

  REAL FUNCTION nug_nu(nu,hmod)

    ! nu*g(nu); useful as an integrand
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu               ! Proxy mass variable
    TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model

    nug_nu=nu*g_nu(nu,hmod)

  END FUNCTION nug_nu

  REAL FUNCTION Mg_nu(nu,hmod)

    ! M(nu)*g(nu); useful as an integrand
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu               ! Proxy mass variable
    TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model

    Mg_nu=M_nu(nu,hmod)*g_nu(nu,hmod)

  END FUNCTION Mg_nu

  REAL FUNCTION nug_nu_on_M(nu,hmod)

    ! nu*g(nu)/M(nu); useful as an integrand
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu               ! Proxy mass variable
    TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model

    nug_nu_on_M=nu*g_nu(nu,hmod)/M_nu(nu,hmod)

  END FUNCTION nug_nu_on_M

  REAL FUNCTION g_nu_on_M(nu,hmod)

    ! g(nu)/M(nu); useful as an integrand
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu               ! Proxy mass variable
    TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model

    g_nu_on_M=g_nu(nu,hmod)/M_nu(nu,hmod)
    
  END FUNCTION g_nu_on_M

  FUNCTION wk_isothermal(x)

    ! The normlaised Fourier Transform of an isothermal profile
    IMPLICIT NONE
    REAL :: wk_isothermal
    REAL, INTENT(IN) :: x

    REAL, PARAMETER :: dx=1e-3

    ! Taylor expansion used for low |x| to avoid cancellation problems

    IF(ABS(x)<ABS(dx)) THEN
       ! Taylor series at low x
       wk_isothermal=1.-(x**2)/18.
    ELSE
       wk_isothermal=Si(x)/x
    END IF

  END FUNCTION wk_isothermal

  FUNCTION wk_isothermal_2(x,y)

    ! The normlaised Fourier Transform of an isothemral profile from x -> y
    IMPLICIT NONE
    REAL :: wk_isothermal_2
    REAL, INTENT(IN) :: x, y

    wk_isothermal_2=(Si(x)-Si(y))/(x-y)

  END FUNCTION wk_isothermal_2

  REAL FUNCTION halo_fraction(itype,m,hmod,cosm)

    ! Mass fraction of a type within a halo
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: itype
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm

    If(itype==field_dmonly .OR. itype==field_matter) THEN
       halo_fraction=1.
    ELSE IF(itype==field_cdm) THEN
       halo_fraction=halo_CDM_fraction(m,hmod,cosm)
    ELSE IF(itype==field_gas) THEN
       halo_fraction=halo_gas_fraction(m,hmod,cosm)
    ELSE IF(itype==field_star) THEN
       halo_fraction=halo_star_fraction(m,hmod,cosm)
    ELSE IF(itype==field_bound_gas) THEN
       halo_fraction=halo_bound_gas_fraction(m,hmod,cosm)
    ELSE IF(itype==field_free_gas) THEN
       halo_fraction=halo_free_gas_fraction(m,hmod,cosm)
    ELSE IF(itype==field_cold_gas) THEN
       halo_fraction=halo_cold_gas_fraction(m,hmod,cosm)
    ELSE IF(itype==field_hot_gas) THEN
       halo_fraction=halo_hot_gas_fraction(m,hmod,cosm)
    ELSE IF(itype==field_static_gas) THEN
       halo_fraction=halo_static_gas_fraction(m,hmod,cosm)
    ELSE IF(itype==field_central_stars) THEN
       halo_fraction=halo_central_star_fraction(m,hmod,cosm)
    ELSE IF(itype==field_satellite_stars) THEN
       halo_fraction=halo_satellite_star_fraction(m,hmod,cosm)
    ELSE
       STOP 'HALO_FRACTION: Error, itype not specified correcntly'
    END IF

  END FUNCTION halo_fraction

  REAL FUNCTION halo_CDM_fraction(m,hmod,cosm)

    ! Mass fraction of a halo in CDM
    IMPLICIT NONE
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: crap

    ! To prevent compile-time warning
    crap=m
    crap=hmod%a

    ! Always the universal value
    halo_CDM_fraction=cosm%om_c/cosm%om_m

  END FUNCTION halo_CDM_fraction

  REAL FUNCTION halo_gas_fraction(m,hmod,cosm)

    ! Mass fraction of a halo in gas
    IMPLICIT NONE
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm

    halo_gas_fraction=halo_bound_gas_fraction(m,hmod,cosm)+halo_free_gas_fraction(m,hmod,cosm)

  END FUNCTION halo_gas_fraction

  REAL FUNCTION halo_bound_gas_fraction(m,hmod,cosm)

    ! Mass fraction of a halo in bound gas
    IMPLICIT NONE
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm    
    REAL :: m0, sigma, beta

    IF(hmod%frac_bound_gas==1) THEN
       ! From Fedeli (2014a)
       m0=1e12
       sigma=3.
       IF(m<m0) THEN
          halo_bound_gas_fraction=0.
       ELSE
          halo_bound_gas_fraction=erf(log10(m/m0)/sigma)*cosm%om_b/cosm%om_m
       END IF
    ELSE IF(hmod%frac_bound_gas==2) THEN
       ! From Schneider & Teyssier (2015)
       M0=HMx_M0(hmod)
       beta=0.6
       halo_bound_gas_fraction=(cosm%om_b/cosm%om_m)/(1.+(M0/m)**beta)
    ELSE IF(hmod%frac_bound_gas==3) THEN
       ! Universal baryon fraction model (account for stellar contribution)
       halo_bound_gas_fraction=cosm%om_b/cosm%om_m-halo_star_fraction(m,hmod,cosm)
    ELSE
       STOP 'HALO_BOUND_GAS_FRACTION: Error, frac_bound_gas not specified correctly'
    END IF

  END FUNCTION halo_bound_gas_fraction

  REAL FUNCTION halo_static_gas_fraction(m,hmod,cosm)

    ! Mass fraction of total gas that is static and bound in haloes
    IMPLICIT NONE
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: bound, cold, hot

    bound=halo_bound_gas_fraction(m,hmod,cosm)
    cold=halo_cold_gas_fraction(m,hmod,cosm)
    hot=halo_hot_gas_fraction(m,hmod,cosm)

    halo_static_gas_fraction=bound-cold-hot

  END FUNCTION halo_static_gas_fraction

  REAL FUNCTION halo_cold_gas_fraction(m,hmod,cosm)

    ! Mass fraction of total gas that is cold and bound in haloes
    IMPLICIT NONE
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(hmod%frac_cold_bound_gas==1) THEN
       ! Constant fraction of cold halo gas
       halo_cold_gas_fraction=hmod%fcold*halo_bound_gas_fraction(m,hmod,cosm)
    ELSE
       STOP 'HALO_COLD_GAS_FRACTION: Error, frac_cold_bound_gas not specified correctly'
    END IF

  END FUNCTION halo_cold_gas_fraction

  REAL FUNCTION halo_hot_gas_fraction(m,hmod,cosm)

    ! Mass fraction of total gas that is hot and bound in haloes
    IMPLICIT NONE
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(hmod%frac_hot_bound_gas==1) THEN
       ! Constant fraction of hot halo gas
       halo_hot_gas_fraction=hmod%fhot*halo_bound_gas_fraction(m,hmod,cosm)
    ELSE
       STOP 'HALO_HOT_GAS_FRACTION: Error, frac_hot_bound_gas not specified correctly'
    END IF

  END FUNCTION halo_hot_gas_fraction

  REAL FUNCTION halo_free_gas_fraction(m,hmod,cosm)

    ! Mass fraction of a halo in free gas
    IMPLICIT NONE
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm

    ! This is necessarily all the gas that is not bound or in stars
    halo_free_gas_fraction=cosm%om_b/cosm%om_m-halo_star_fraction(m,hmod,cosm)-halo_bound_gas_fraction(m,hmod,cosm)
    IF(halo_free_gas_fraction<0.) halo_free_gas_fraction=0.

  END FUNCTION halo_free_gas_fraction

  REAL FUNCTION halo_star_fraction(m,hmod,cosm)

    ! Mass fraction of a halo in stars
    IMPLICIT NONE
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: m0, sigma, A
    REAL :: crap

    crap=cosm%A

    IF(hmod%frac_stars==1 .OR. hmod%frac_stars==3) THEN
       ! Fedeli (2014)
       A=HMx_Astar(m,hmod)
       !m0=hmod%Mstar
       m0=HMx_Mstar(hmod)
       sigma=hmod%sstar
       halo_star_fraction=A*exp(-((log10(m/m0))**2)/(2.*sigma**2))
       IF(hmod%frac_stars==3) THEN
          ! Suggested by Ian, the relation I have is for the central stellar mass
          ! in reality this saturates for high-mass haloes (due to satellite contribution)
          IF(halo_star_fraction<A/3. .AND. m>m0) halo_star_fraction=A/3.
       END IF
    ELSE IF(hmod%frac_stars==2) THEN
       ! Constant star fraction
       A=0.005
       halo_star_fraction=A
    ELSE IF(hmod%frac_stars==4) THEN
       ! No stars (actually, things crash with exactly zero stars)
       halo_star_fraction=1e-4
    ELSE
       STOP 'HALO_STAR_FRACTION: Error, frac_stars specified incorrectly'
    END IF
    
  END FUNCTION halo_star_fraction

  REAL FUNCTION halo_central_star_fraction(m,hmod,cosm)

    ! Mass fraction of a halo in central stars
    IMPLICIT NONE
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: r

    IF(hmod%frac_central_stars==1) THEN
       ! All stellar mass in centrals
       halo_central_star_fraction=halo_star_fraction(m,hmod,cosm)
    ELSE IF(hmod%frac_central_stars==2) THEN
       ! Schnedier et al. (2018)
       IF(M<=HMx_Mstar(hmod)) THEN
          ! For M<M* all galaxy mass is in centrals
          r=1.
       ELSE
          ! Otherwise there is a power law, eta ~ -0.3, higher mass haloes have more mass in satellites
          r=(M/HMx_Mstar(hmod))**hmod%eta
       END IF
!!$       IF(r>1.) THEN
!!$          STOP 'HALO_CENTRAL_STAR_FRACTION: Error, r cannot be greater than one'
!!$       END IF
       halo_central_star_fraction=r*halo_star_fraction(m,hmod,cosm)
    ELSE
       STOP 'HALO_CENTRAL_STAR_FRACTION: Error, frac_central_stars specified incorrectly'
    END IF
    
  END FUNCTION halo_central_star_fraction

  REAL FUNCTION halo_satellite_star_fraction(m,hmod,cosm)

    ! Mass fraction of a halo in satellite stars
    IMPLICIT NONE
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm

    halo_satellite_star_fraction=halo_star_fraction(m,hmod,cosm)-halo_central_star_fraction(m,hmod,cosm)

  END FUNCTION halo_satellite_star_fraction

  REAL FUNCTION halo_HI_fraction(M,hmod,cosm)

    ! Dimensionless M_HI/M_halo
    ! TODO: This should probably be related to the cold gas fraction
    IMPLICIT NONE
    REAL, INTENT(IN) :: M ! Host-halo mass
    TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model
    TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology
    REAL :: alpha, M0, Mmin
    REAL :: crap

    ! Prevent compile-time warnings
    crap=cosm%A

    IF(hmod%frac_HI==1) THEN
       ! Simple model with hard truncation
       IF(m>=hmod%HImin .AND. m<=hmod%HImax) THEN
          halo_HI_fraction=1.
       ELSE
          halo_HI_fraction=0.
       END IF
    ELSE IF(hmod%frac_HI==2 .OR. hmod%frac_HI==3) THEN
       ! From Villaescusa-Navarro et al. (2018; 1804.09180)
       ! 2 - Just z=0 results
       ! 3 - Using z evolution
       IF(hmod%frac_HI==2 .OR. hmod%z==0.) THEN
          ! FoF values from Table 1
          !alpha=0.24
          !M0=4.3e10
          !Mmin=2e12
          ! FoF-SO values from Table 1
          alpha=0.16
          M0=4.1e10
          Mmin=2.4e12          
       ELSE IF(hmod%z==1.) THEN
          ! FoF values from Table 1
          !alpha=0.53
          !M0=1.5e10
          !Mmin=6e11
          ! FoF-SO values from Table 1
          alpha=0.43
          M0=1.8e10
          Mmin=8.6e11          
       ELSE IF(hmod%z==2.) THEN
          ! FoF values from Table 1
          !alpha=0.60
          !M0=1.3e10
          !Mmin=3.6e11
          ! FoF-SO values from Table 1
          alpha=0.51
          M0=1.5e10
          Mmin=4.6e11        
       ELSE IF(hmod%z==3.) THEN
          ! FoF values from Table 1
          !alpha=0.76
          !M0=2.9e9
          !Mmin=6.7e10
          ! FoF-SO values from Table 1
          alpha=0.69
          M0=3.7e9
          Mmin=9.6e10
       ELSE IF(hmod%z==4.) THEN
          ! FoF values from Table 1
          !alpha=0.79
          !M0=1.4e9
          !Mmin=2.1e10
          ! FoF-SO values from Table 1
          alpha=0.61
          M0=4.5e9
          Mmin=7.6e10         
       ELSE IF(hmod%z==5.) THEN
          ! FoF values from Table 1
          !alpha=0.74
          !M0=1.9e9
          !Mmin=2e10
          ! FoF-SO values from Table 1
          alpha=0.59
          M0=4.1e9
          Mmin=5.4e10          
       ELSE
          STOP 'HALO_HI_FRACTION: Error, redshift not supported here'
       END IF
       halo_HI_fraction=(M0/M)*((M/Mmin)**alpha)*exp(-(M/Mmin)**(-0.35))
    ELSE
       STOP 'HALO_HI_FRACTION: Error, frac_HI specified incorrectly'
    END IF

  END FUNCTION halo_HI_fraction

  FUNCTION integrate_hmod(a,b,f,hmod,acc,iorder)

    ! Integrates between a and b until desired accuracy is reached
    ! Stores information to reduce function calls
    IMPLICIT NONE
    REAL :: integrate_hmod
    REAL, INTENT(IN) :: a, b, acc
    INTEGER, INTENT(IN) :: iorder
    TYPE(halomod), INTENT(INOUT) :: hmod
    INTEGER :: i, j
    INTEGER :: n
    REAL :: x, dx
    REAL :: f1, f2, fx
    DOUBLE PRECISION :: sum_n, sum_2n, sum_new, sum_old

    INTEGER, PARAMETER :: jmin=5
    INTEGER, PARAMETER :: jmax=30

    INTERFACE
       REAL FUNCTION f(x,hmod)
         IMPORT :: halomod
         REAL, INTENT(IN) :: x
         TYPE(halomod), INTENT(INOUT) :: hmod
       END FUNCTION f
    END INTERFACE

    IF(a==b) THEN

       ! Fix the answer to zero if the integration limits are identical
       integrate_hmod=0.

    ELSE

       ! Set the sum variable for the integration
       sum_2n=0.d0
       sum_n=0.d0
       sum_old=0.d0
       sum_new=0.d0

       DO j=1,jmax

          ! Note, you need this to be 1+2**n for some integer n
          ! j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
          n=1+2**(j-1)

          ! Calculate the dx interval for this value of 'n'
          dx=(b-a)/REAL(n-1)

          IF(j==1) THEN

             ! The first go is just the trapezium of the end points
             f1=f(a,hmod)
             f2=f(b,hmod)
             sum_2n=0.5d0*(f1+f2)*dx
             sum_new=sum_2n

          ELSE

             ! Loop over only new even points to add these to the integral
             DO i=2,n,2
                x=a+(b-a)*REAL(i-1)/REAL(n-1)
                fx=f(x,hmod)
                sum_2n=sum_2n+fx
             END DO

             ! Now create the total using the old and new parts
             sum_2n=sum_n/2.d0+sum_2n*dx

             ! Now calculate the new sum depending on the integration order
             IF(iorder==1) THEN  
                sum_new=sum_2n
             ELSE IF(iorder==3) THEN         
                sum_new=(4.d0*sum_2n-sum_n)/3.d0 ! This is Simpson's rule and cancels error
             ELSE
                STOP 'INTEGRATE_HMOD: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (ABS(-1.d0+sum_new/sum_old)<acc)) THEN
             EXIT
          ELSE IF(j==jmax) THEN
             STOP 'INTEGRATE_HMOD: Integration timed out'
          ELSE
             ! Integral has not converged so store old sums and reset sum variables
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.d0
          END IF

       END DO

       integrate_hmod=REAL(sum_new)

    END IF

  END FUNCTION integrate_hmod

  FUNCTION integrate_hmod_cosm(a,b,f,hmod,cosm,acc,iorder)

    ! Integrates between a and b until desired accuracy is reached
    ! Stores information to reduce function calls
    IMPLICIT NONE
    REAL :: integrate_hmod_cosm
    REAL, INTENT(IN) :: a, b, acc
    INTEGER, INTENT(IN) :: iorder
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: i, j
    INTEGER :: n
    REAL :: x, dx
    REAL :: f1, f2, fx
    DOUBLE PRECISION :: sum_n, sum_2n, sum_new, sum_old

    INTEGER, PARAMETER :: jmin=5
    INTEGER, PARAMETER :: jmax=30

    INTERFACE
       REAL FUNCTION f(x,hmod,cosm)
         IMPORT :: halomod
         IMPORT :: cosmology
         REAL, INTENT(IN) :: x
         TYPE(halomod), INTENT(INOUT) :: hmod
         TYPE(cosmology), INTENT(INOUT) :: cosm
       END FUNCTION f
    END INTERFACE

    IF(a==b) THEN

       ! Fix the answer to zero if the integration limits are identical
       integrate_hmod_cosm=0.

    ELSE

       ! Set the sum variable for the integration
       sum_2n=0.d0
       sum_n=0.d0
       sum_old=0.d0
       sum_new=0.d0

       DO j=1,jmax

          ! Note, you need this to be 1+2**n for some integer n
          ! j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
          n=1+2**(j-1)

          ! Calculate the dx interval for this value of 'n'
          dx=(b-a)/REAL(n-1)

          IF(j==1) THEN

             ! The first go is just the trapezium of the end points
             f1=f(a,hmod,cosm)
             f2=f(b,hmod,cosm)
             sum_2n=0.5d0*(f1+f2)*dx
             sum_new=sum_2n

          ELSE

             ! Loop over only new even points to add these to the integral
             DO i=2,n,2
                x=a+(b-a)*REAL(i-1)/REAL(n-1)
                fx=f(x,hmod,cosm)
                sum_2n=sum_2n+fx
             END DO

             ! Now create the total using the old and new parts
             sum_2n=sum_n/2.d0+sum_2n*dx

             ! Now calculate the new sum depending on the integration order
             IF(iorder==1) THEN  
                sum_new=sum_2n
             ELSE IF(iorder==3) THEN         
                sum_new=(4.d0*sum_2n-sum_n)/3.d0 ! This is Simpson's rule and cancels error
             ELSE
                STOP 'INTEGRATE_HMOD: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (ABS(-1.d0+sum_new/sum_old)<acc)) THEN
             EXIT
          ELSE IF(j==jmax) THEN
             STOP 'INTEGRATE_HMOD: Integration timed out'
          ELSE
             ! Integral has not converged so store old sums and reset sum variables
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.d0
          END IF

       END DO

       integrate_hmod_cosm=REAL(sum_new)

    END IF

  END FUNCTION integrate_hmod_cosm

  FUNCTION integrate_hmod_cosm_exp(a,b,f,hmod,cosm,acc,iorder)

    ! Integrates between a and b until desired accuracy is reached
    ! Stores information to reduce function calls
    IMPLICIT NONE
    REAL :: integrate_hmod_cosm_exp
    REAL, INTENT(IN) :: a, b, acc
    INTEGER, INTENT(IN) :: iorder
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: i, j
    INTEGER :: n
    REAL :: x, dx
    REAL :: f1, f2, fx
    DOUBLE PRECISION :: sum_n, sum_2n, sum_new, sum_old

    INTEGER, PARAMETER :: jmin=5
    INTEGER, PARAMETER :: jmax=30

    INTERFACE
       REAL FUNCTION f(nu,hmod,cosm)
         IMPORT :: halomod
         IMPORT :: cosmology
         REAL, INTENT(IN) :: nu
         TYPE(halomod), INTENT(INOUT) :: hmod
         TYPE(cosmology), INTENT(INOUT) :: cosm
       END FUNCTION f
    END INTERFACE

    IF(a==b) THEN

       ! Fix the answer to zero if the integration limits are identical
       integrate_hmod_cosm_exp=0.

    ELSE

       ! Set the sum variable for the integration
       sum_2n=0.d0
       sum_n=0.d0
       sum_old=0.d0
       sum_new=0.d0

       DO j=1,jmax

          ! Note, you need this to be 1+2**n for some integer n
          ! j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
          n=1+2**(j-1)

          ! Calculate the dx interval for this value of 'n'
          dx=(b-a)/REAL(n-1)

          IF(j==1) THEN

             ! The first go is just the trapezium of the end points
             f1=f(exp(a),hmod,cosm)*exp(a)
             f2=f(exp(b),hmod,cosm)*exp(b)
             sum_2n=0.5d0*(f1+f2)*dx
             sum_new=sum_2n

          ELSE

             ! Loop over only new even points to add these to the integral
             DO i=2,n,2
                x=a+(b-a)*REAL(i-1)/REAL(n-1)
                fx=f(exp(x),hmod,cosm)*exp(x)
                sum_2n=sum_2n+fx
             END DO

             ! Now create the total using the old and new parts
             sum_2n=sum_n/2.d0+sum_2n*dx

             ! Now calculate the new sum depending on the integration order
             IF(iorder==1) THEN  
                sum_new=sum_2n
             ELSE IF(iorder==3) THEN         
                sum_new=(4.d0*sum_2n-sum_n)/3.d0 ! This is Simpson's rule and cancels error
             ELSE
                STOP 'INTEGRATE_HMOD: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (ABS(-1.d0+sum_new/sum_old)<acc)) THEN
             EXIT
          ELSE IF(j==jmax) THEN
             STOP 'INTEGRATE_HMOD: Integration timed out'
          ELSE
             ! Integral has not converged so store old sums and reset sum variables
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.d0
          END IF

       END DO

       integrate_hmod_cosm_exp=REAL(sum_new)

    END IF

  END FUNCTION integrate_hmod_cosm_exp

  REAL FUNCTION integrate_scatter(c,dc,ih,k,m,rv,hmod,cosm,acc,iorder)

    ! Integrates between a and b until desired accuracy is reached
    ! Stores information to reduce function calls
    IMPLICIT NONE
    REAL, INTENT(IN) :: c
    REAL, INTENT(IN) :: dc    
    INTEGER, INTENT(IN) :: ih(2)
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL, INTENT(IN) :: acc
    INTEGER, INTENT(IN) :: iorder
    REAL :: a, b
    INTEGER :: i, j
    INTEGER :: n
    REAL :: x, dx
    REAL :: f1, f2, fx
    DOUBLE PRECISION :: sum_n, sum_2n, sum_new, sum_old

    INTEGER, PARAMETER :: jmin=5
    INTEGER, PARAMETER :: jmax=30
    REAL, PARAMETER :: nsig=5

    a=c/(1.+nsig*dc)
    b=c*(1.+nsig*dc)

    IF(a==b) THEN

       ! Fix the answer to zero if the integration limits are identical
       integrate_scatter=0.

    ELSE

       ! Set the sum variable for the integration
       sum_2n=0.d0
       sum_n=0.d0
       sum_old=0.d0
       sum_new=0.d0

       DO j=1,jmax

          ! Note, you need this to be 1+2**n for some integer n
          ! j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
          n=1+2**(j-1)

          ! Calculate the dx interval for this value of 'n'
          dx=(b-a)/REAL(n-1)

          IF(j==1) THEN

             ! The first go is just the trapezium of the end points
             f1=scatter_integrand(a,c,dc,ih,k,m,rv,hmod,cosm)
             f2=scatter_integrand(b,c,dc,ih,k,m,rv,hmod,cosm)
             sum_2n=0.5d0*(f1+f2)*dx
             sum_new=sum_2n

          ELSE

             ! Loop over only new even points to add these to the integral
             DO i=2,n,2
                x=a+(b-a)*REAL(i-1)/REAL(n-1)
                fx=scatter_integrand(x,c,dc,ih,k,m,rv,hmod,cosm)
                sum_2n=sum_2n+fx
             END DO

             ! Now create the total using the old and new parts
             sum_2n=sum_n/2.d0+sum_2n*dx

             ! Now calculate the new sum depending on the integration order
             IF(iorder==1) THEN  
                sum_new=sum_2n
             ELSE IF(iorder==3) THEN         
                sum_new=(4.d0*sum_2n-sum_n)/3.d0 ! This is Simpson's rule and cancels error
             ELSE
                STOP 'INTEGRATE_SCATTER: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (ABS(-1.d0+sum_new/sum_old)<acc)) THEN
             EXIT
          ELSE IF(j==jmax) THEN
             STOP 'INTEGRATE_SCATTER: Integration timed out'
          ELSE
             ! Integral has not converged so store old sums and reset sum variables
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.d0
          END IF

       END DO

       integrate_scatter=REAL(sum_new)

    END IF

  END FUNCTION integrate_scatter

  REAL FUNCTION scatter_integrand(c,mean_c,sigma_lnc,ih,k,m,rv,hmod,cosm)

    ! Integrand for computing halo profiles with scatter
    IMPLICIT NONE
    REAL, INTENT(IN) :: c
    REAL, INTENT(IN) :: mean_c
    REAL, INTENT(IN) :: sigma_lnc
    INTEGER, INTENT(IN) :: ih(2)
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: wk(2), pc, rs
    INTEGER :: j

    !Halo profiles
    DO j=1,2
       rs=rv/c
       wk(j)=win_type(.FALSE.,ih(j),k,m,rv,rs,hmod,cosm)
    END DO

    !Probability distribution
    pc=lognormal(c,mean_c,sigma_lnc)

    !The full integrand
    scatter_integrand=wk(1)*wk(2)*pc

  END FUNCTION scatter_integrand

END MODULE HMx
