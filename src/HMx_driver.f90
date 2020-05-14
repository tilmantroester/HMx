PROGRAM HMx_driver

   USE constants
   USE HMx
   USE Limber
   USE file_info
   USE random_numbers
   USE cosmic_emu_stuff
   USE calculus_table
   USE string_operations
   USE special_functions
   USE basic_operations
   USE owls
   USE owls_extras
   USE interpolate
   USE cosmology_functions
   USE array_operations
   USE multidark_stuff

   IMPLICIT NONE

   ! Parameter definitions
   CHARACTER(len=256) :: mode, halomodel, cosmo
   INTEGER :: iimode, iicosmo, iihm

   CALL get_command_argument(1, mode)
   IF (mode == '') THEN
      iimode = -1
   ELSE
      READ (mode, *) iimode
   END IF

   CALL get_command_argument(2, cosmo)
   IF (cosmo == '') THEN
      iicosmo = -1
   ELSE
      READ (cosmo, *) iicosmo
   END IF

   CALL get_command_argument(3, halomodel)
   IF (halomodel == '') THEN
      iihm = -1
   ELSE
      READ (halomodel, *) iihm
   END IF

   ! Initial white space
   WRITE (*, *)

   ! Choose mode
   IF (iimode == -1) THEN
      WRITE (*, *) 'HMx_DRIVER: Choose what to do'
      WRITE (*, *) '============================='
      WRITE (*, *) '  0 - Gravity-only power spectrum at z=0'
      WRITE (*, *) '  1 - 3D Matter power spectrum over multiple z'
      WRITE (*, *) '  2 - Hydrodynamical halo model'
      WRITE (*, *) '  3 - Run diagnostics for haloes'
      WRITE (*, *) '  4 - Do random baryon parameters for bug testing'
      WRITE (*, *) '  5 - HMx2020'
      WRITE (*, *) '  6 - HMx2020 as a function of AGN temperature'
      WRITE (*, *) '  7 - Do general angular cross correlation'
      WRITE (*, *) '  8 - Angular cross correlation as a function of cosmology'
      WRITE (*, *) '  9 - Breakdown angular correlations in halo mass'
      WRITE (*, *) ' 10 - Breakdown angular correlations in redshift'
      WRITE (*, *) ' 11 - Do general angular cross correlation and correlation functions'
      WRITE (*, *) ' 12 - Triad'
      WRITE (*, *) ' 13 - Cross-correlation coefficient'
      WRITE (*, *) ' 14 - Effect of parameter variations on hydro power'
      WRITE (*, *) ' 15 - HMx2020 matter-matter suppression as a function of cosmology'
      WRITE (*, *) ' 16 - TEST: Equality of matter vs. DMONLY implementation'
      WRITE (*, *) ' 17 - 3D spectra for user choice of fields'
      WRITE (*, *) ' 18 - 3D bias'
      WRITE (*, *) ' 19 - Make CCL benchmark data'
      WRITE (*, *) ' 20 - Make data for Ma et al. (2015) Fig. 1'
      WRITE (*, *) ' 21 - PAPER: HMx2020 across AGN temperature range for matter'
      WRITE (*, *) ' 22 - PAPER: HMx2020 across AGN temperature range for matter, electron pressure'
      WRITE (*, *) ' 23 - Produce DE response results from Mead (2017)'
      WRITE (*, *) ' 24 - TEST: Projection'
      WRITE (*, *) ' 25 - Halo-void model'
      WRITE (*, *) ' 26 - TEST: DMONLY spectra; HMcode'
      WRITE (*, *) ' 27 - Comparison with Mira Titan nodes'
      WRITE (*, *) ' 28 - Comparison with Franken Emu nodes'
      WRITE (*, *) ' 29 - Comparison with random Mira Titan cosmology'
      WRITE (*, *) ' 30 - Comparison with random Franken Emu cosmology'
      WRITE (*, *) ' 31 - PAPER: Hydro power with varying upper-mass limits on halo mass integrals'
      WRITE (*, *) ' 32 - PAPER: Hydrodynamical halo model power for baseline model'
      WRITE (*, *) ' 33 - PAPER: Effect of parameter variations on baseline model hydro power'
      WRITE (*, *) ' 34 - Write Tilman HMx hydro parameters as a function of T_AGN and z'
      WRITE (*, *) ' 35 - Power spectra of cored halo profiles'
      WRITE (*, *) ' 36 - TEST: hydro spectra'
      WRITE (*, *) ' 37 - CFHTLenS correlation functions (Kilbinger et al. 2013)'
      WRITE (*, *) ' 38 - Tilman AGN model triad for all feedback models'
      WRITE (*, *) ' 39 - Tilman AGN model Triad for all feedback models (Triad 3 ell)'
      WRITE (*, *) ' 40 - Halo bias'
      WRITE (*, *) ' 41 - PAPER: Matter, pressure power as a function of sigma8'
      WRITE (*, *) ' 42 - PAPER: Contributions to k-k C(l) integral'
      WRITE (*, *) ' 43 - PAPER: Contributions to k-y C(l) integral'
      WRITE (*, *) ' 44 - Triad with Tilman model'
      WRITE (*, *) ' 45 - Comparison of Sheth-Tormen vs. Tinker mass function'
      WRITE (*, *) ' 46 - Mass function and bias plots'
      WRITE (*, *) ' 47 - Make CMB lensing to compare with CAMB'
      WRITE (*, *) ' 48 - HI bias'
      WRITE (*, *) ' 49 - HI mass fractions'
      WRITE (*, *) ' 50 - Mass function changes with Lbox'
      WRITE (*, *) ' 51 - Compare power with and without scatter'
      WRITE (*, *) ' 52 - Hydrodynamical halo model with BAHAMAS k range'
      WRITE (*, *) ' 53 - Breakdown 3D hydro power in halo mass bins'
      WRITE (*, *) ' 54 - Trispectrum test'
      WRITE (*, *) ' 55 - Triad for all feedback models'
      WRITE (*, *) ' 56 - Triad for all feedback models (Triad 3 ell)'
      WRITE (*, *) ' 57 - Direct integration of measured 3D BAHAMAS spectra: Triad 3'
      WRITE (*, *) ' 58 - Multiple same fields check'
      WRITE (*, *) ' 59 - Tinker (2010) bias plot check'
      WRITE (*, *) ' 60 - Make data for Limber comparison with CCL'
      WRITE (*, *) ' 61 - Direct integration of measured 3D BAHAMAS spectra: Triad 4'
      WRITE (*, *) ' 62 - Direct integration of measured 3D BAHAMAS spectra: Triad 5'
      WRITE (*, *) ' 63 - Contributions to y-y C(l) integral'
      WRITE (*, *) ' 64 - CIB 545 GHz'
      WRITE (*, *) ' 65 - Contributions to Limber integrals from BAHAMAS'
      WRITE (*, *) ' 66 - Direct integration of measured 3D BAHAMAS spectra: Triad 5 but larger ell range'
      WRITE (*, *) ' 67 - Compute matter power using non-linear halo bias'
      WRITE (*, *) ' 68 - Compute non-linear halo bias integrand'
      WRITE (*, *) ' 69 - Matter, halo power spectra with non-linear bias Multidark comparison'
      WRITE (*, *) ' 70 - Comparison with Cosmic Emu nodes'
      WRITE (*, *) ' 71 - Comparison with random Cosmic Emu cosmology'
      WRITE (*, *) ' 72 - Run halo model over lots of random cosmological parameters'
      WRITE (*, *) ' 73 - Check halo-mass function amplitude parameter'
      WRITE (*, *) ' 74 - Matter, halo power spectra with non-linear bias (low sigma_8) Multidark comparison'
      WRITE (*, *) ' 75 - Compare HMcode (2016) HMx vs CAMB'
      WRITE (*, *) ' 76 - Comparison with massless neutrino Mira Titan nodes'
      WRITE (*, *) ' 77 - DMONLY comparison at z = 0 with regular halo model and user choice cosmology' 
      WRITE (*, *) ' 78 - DMONLY comparison at multiple z with regular halo model and user choice cosmology'    
      WRITE (*, *) ' 79 - DMONLY comparison at z = 0 with HMcode (2016) and user choice cosmology'
      WRITE (*, *) ' 80 - DMONLY comparison at multiple z with HMcode (2016) and user choice cosmology'
      WRITE (*, *) ' 81 - DMONLY comparison at z = 0 with regular halo model and boring cosmology'
      WRITE (*, *) ' 82 - DMONLY comparison at multiple z with regular halo model and boring cosmology'
      WRITE (*, *) ' 83 - DMONLY comparison at z = 0 with HMcode (2016) and boring cosmology'
      WRITE (*, *) ' 84 - DMONLY comparison at multiple z with HMcode (2016) and boring cosmology'
      WRITE (*, *) ' 85 - Power comparison at multiple z for bump cosmologies'
      WRITE (*, *) ' 86 - PAPER: Run diagnostics for haloes'
      WRITE (*, *) ' 87 - Hydro power with varying upper-mass limits on halo mass integrals'
      WRITE (*, *) ' 88 - Matter, pressure power as a function of sigma8'
      WRITE (*, *) ' 89 - Compare to Mira Titan as a function of neutrino mass'
      WRITE (*, *) ' 90 - Compare HMcode (2015) HMx vs CAMB'
      WRITE (*, *) ' 91 - Compare HALOFIT (Smith et al.) HMx vs CAMB'
      WRITE (*, *) ' 92 - Compare HALOFIT (Bird et al.) HMx vs CAMB'
      WRITE (*, *) ' 93 - Compare HALOFIT (Takahashi et al.) HMx vs CAMB'
      WRITE (*, *) ' 94 - PAPER: Big Franken Emu node comparison'
      WRITE (*, *) ' 95 - PAPER: Big Mira Titan node comparison'
      WRITE (*, *) ' 96 - Big Cosmic Emu node comparison'
      WRITE (*, *) ' 97 - Compare HALOFIT (CAMB parameters) HMx vs CAMB'
      WRITE (*, *) ' 98 - Compare HMcode (2020) HMx vs CAMB'
      WRITE (*, *) ' 99 - PAPER: Emulator mean and variance'
      WRITE (*, *) '100 - Comparison with M000 of FrankenEmu'
      WRITE (*, *) '101 - HMcode speed tests'
      WRITE (*, *) '102 - Gas power spectra in bound and ejected'
      WRITE (*, *) '103 - Compare HMcode (2016) with neutrino fix HMx vs CAMB'
      READ (*, *) iimode
      WRITE (*, *) '============================'
      WRITE (*, *)
   END IF

   IF (iimode == 0) THEN
      CALL power_single(iicosmo, iihm)
   ELSE IF (iimode == 1) THEN
      CALL power_multiple(iicosmo, iihm)
   ELSE IF (iimode == 2 .OR. iimode == 5 .OR. iimode == 32 .OR. iimode == 52) THEN
      CALL hydro_stuff(iimode, iicosmo, iihm)
   ELSE IF (iimode == 6 .OR. iimode == 21 .OR. iimode == 22) THEN
      CALL hydro_across_temperature_range(iimode, iicosmo, iihm)
   ELSE IF (iimode == 15) THEN
      CALL hydro_suppression_as_a_function_of_cosmology(iicosmo, iihm)
   ELSE IF (iimode == 16) THEN
      CALL compare_matter_vs_DMONLY_implementation(iicosmo)
   ELSE IF (iimode == 3 .OR. iimode == 86) THEN
      CALL halo_stuff(iimode, iicosmo, iihm)
   ELSE IF (iimode == 4 .OR. iimode == 72) THEN
      CALL test_random_cosmologies(iimode, iicosmo, iihm)
   ELSE IF(is_in_array(iimode, [7, 8, 9, 10, 11, 37, 42, 43, 47, 63, 64])) THEN
      CALL general_projection(iimode, iicosmo, iihm)
   ELSE IF (is_in_array(iimode, [12, 38, 39, 44, 55, 56])) THEN
      CALL triad_stuff(iimode, iicosmo, iihm)
   ELSE IF (iimode == 13) THEN
      CALL cross_correlation_coefficient(iicosmo, iihm)
   ELSE IF (iimode == 14 .OR. iimode == 33) THEN
      CALL baryon_parameter_variations(iimode, iicosmo, iihm)
   ELSE IF (iimode == 17) THEN
      CALL spectra_fields_3D(iicosmo, iihm)
   ELSE IF (iimode == 18 .OR. iimode == 40 .OR. iimode == 48) THEN
      CALL bias_fields_3D(iimode, iicosmo, iihm)
   ELSE IF (iimode == 19) THEN
      CALL create_CCL_benchmark(iicosmo, iihm)
   ELSE IF (iimode == 20) THEN
      CALL Ma2015_Fig1(iicosmo, iihm)
   ELSE IF (iimode == 23) THEN
      CALL Mead2017(iicosmo, iihm)
   ELSE IF (iimode == 24) THEN
      CALL projection_tests(iicosmo, iihm)
   ELSE IF (iimode == 25) THEN
      CALL halo_void_model(iicosmo, iihm)
   ELSE IF (iimode == 26) THEN
      CALL halo_model_tests(iicosmo, iihm)
   ELSE IF (is_in_array(iimode, [27, 28, 29, 30, 70, 71, 76, 89, 100])) THEN
      CALL emulator_comparison(iimode, iicosmo, iihm)
   ELSE IF (iimode == 31 .OR. iimode == 53 .OR. iimode == 87) THEN
      CALL power_breakdown_halomass(iimode, iicosmo, iihm)
   ELSE IF (iimode == 34) THEN
      CALL write_HMx_variations(iicosmo, iihm)
   ELSE IF (iimode == 35) THEN
      CALL halo_cores(iicosmo, iihm)
   ELSE IF (iimode == 36) THEN
      CALL hydro_tests(iicosmo, iihm)
   ELSE IF (iimode == 41 .OR. iimode == 88) THEN
      CALL power_cosmology_dependence(iimode, iicosmo, iihm)
   ELSE IF (iimode == 45) THEN
      CALL power_different_halo_mass_functions(iicosmo, iihm)
   ELSE IF (iimode == 46) THEN
      CALL mass_function_plots(iicosmo, iihm)
   ELSE IF (iimode == 49) THEN
      CALL HI_mass_fractions(iicosmo, iihm)
   ELSE IF (iimode == 50) THEN
      CALL mass_function_Lbox(iicosmo, iihm)
   ELSE IF (iimode == 51) THEN
      CALL power_scatter(iicosmo, iihm)
   ELSE IF (iimode == 54) THEN
      CALL Trispectrum_test(iicosmo, iihm)
   ELSE IF (iimode == 57 .OR. iimode == 61 .OR. iimode == 62 .OR. iimode == 65 .OR. iimode == 66) THEN
      CALL Cl_direct_integration(iimode, iicosmo)
   ELSE IF (iimode == 58) THEN
      CALL check_fields_symmetry(iicosmo, iihm)
   ELSE IF (iimode == 59) THEN
      CALL Tinker2010_Fig1(iicosmo, iihm)
   ELSE IF (iimode == 60) THEN
      CALL Limber_CCL_comparison(iicosmo)
   ELSE IF (iimode == 67) THEN
      CALL Non_linear_halo_bias_model(iicosmo, iihm)
   ELSE IF (iimode == 68) THEN
      CALL non_linear_halo_bias_integrand(iicosmo, iihm)
   ELSE IF (iimode == 69 .OR. iimode == 74) THEN
      CALL halo_power_multidark(iimode, iicosmo, iihm)
   ELSE IF (iimode == 73) THEN
      CALL HMF_amplitude(iicosmo, iihm)
   ELSE IF (is_in_array(iimode, [75, 90, 91, 92, 93, 97, 98, 103])) THEN
      CALL compare_HMx_CAMB(iimode, iicosmo)
   ELSE IF (iimode == 77 .OR. iimode == 79 .OR. iimode == 81 .OR. iimode == 83) THEN
      CALL power_single_comparison(iimode, iicosmo, iihm)
   ELSE IF (iimode == 78 .OR. iimode == 80 .OR. iimode == 82 .OR. iimode == 84 .OR. iimode == 85) THEN
      CALL power_multiple_comparison(iimode, iicosmo, iihm)
   ELSE IF (is_in_array(iimode, [94, 95, 96])) THEN
      CALL big_emulator_comparison(iimode)
   ELSE IF (iimode == 99) THEN
      CALL emulator_mean_variance()
   ELSE IF (iimode == 101) THEN
      CALL halomod_speed_tests(iicosmo, iihm)
   ELSE IF (iimode == 102) THEN
      CALL gas_power_spectra(iicosmo, iihm)
   ELSE
      STOP 'HMx_DRIVER: Error, you have specified the mode incorrectly'
   END IF

CONTAINS

   SUBROUTINE gas_power_spectra(icos, ihm)

      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: icos
      INTEGER, INTENT(INOUT) :: ihm
      REAL, ALLOCATABLE :: k(:), a(:)
      REAL, ALLOCATABLE :: Pk(:, :, :, :)
      TYPE(cosmology) :: cosm
      INTEGER :: ik

      REAL, PARAMETER :: kmin = 1e-2
      REAL, PARAMETER :: kmax = 1e1
      INTEGER, PARAMETER :: nk = 128
      REAL, PARAMETER :: amin = 1.
      REAL, PARAMETER :: amax = 1.
      INTEGER, PARAMETER :: na = 1
      INTEGER, PARAMETER :: fields(3) = [field_gas, field_bound_gas, field_ejected_gas]
      INTEGER, PARAMETER :: nf = 3
      LOGICAL, PARAMETER :: verbose = .TRUE.
      INTEGER, PARAMETER :: unit = 7
      CHARACTER(len=256), PARAMETER :: outfile = 'data/power_gas.dat'

      CALL fill_array_log(kmin, kmax, k, nk)
      CALL fill_array(amin, amax, a, na)

      CALL assign_init_cosmology(icos, cosm, verbose)

      CALL calculate_HMx(fields, k, a, Pk, nf, nk, na, cosm, ihm)

      OPEN(unit, file=outfile)
      DO ik = 1, nk
         WRITE(unit, *) k(ik), Pk(1, 1, ik, 1), Pk(2, 2, ik, 1), Pk(3, 3, ik, 1), Pk(2, 3, ik, 1)
      END DO
      CLOSE(unit)

   END SUBROUTINE gas_power_spectra

   SUBROUTINE halomod_speed_tests(icos, ihm)

      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: icos
      INTEGER, INTENT(INOUT) :: ihm
      REAL, ALLOCATABLE :: k(:), a(:), Pk(:, :)
      REAL :: t1, t2, t
      INTEGER :: ia, na
      TYPE(cosmology) :: cosm
      REAL, PARAMETER :: kmin = 1e-3
      REAL, PARAMETER :: kmax = 1e2
      INTEGER, PARAMETER :: nk = 128
      REAL, PARAMETER :: amin = 0.1
      REAL, PARAMETER :: amax = 1.0
      INTEGER, PARAMETER :: naa = 7
      INTEGER, PARAMETER :: version = HMcode2016
      LOGICAL, PARAMETER :: verbose = .FALSE.
      INTEGER, PARAMETER :: unit = 7
      CHARACTER(len=256), PARAMETER :: outfile = 'data/halomod_timing.dat'

      ! Fill arrays for k
      CALL fill_array(kmin, kmax, k, nk)

      ! Loop over different number of a required
      OPEN(unit, file=outfile)
      DO ia = 0, naa

         ! Number of a at this point in the loop
         na = 2**ia

         ! Fill array for a
         IF (na == 1) THEN
            ALLOCATE (a(na))
            a(1) = 1.
         ELSE
            CALL fill_array(amin, amax, a, na)
         END IF

         ! Get the intial time
         CALL cpu_time(t1)

         ! Assign cosmology and do the HMcode calculation
         !CALL assign_cosmology(icos, cosm, verbose)
         !CALL init_cosmology(cosm)
         CALL assign_init_cosmology(icos, cosm, verbose)
         CALL calculate_halomod(k, a, Pk, nk, na, cosm, ihm)

         ! Get the final time
         CALL cpu_time(t2)

         ! Total time for calculation
         t = t2-t1

         ! Write to screen
         WRITE(*, *) 'HALOMOD_SPEED_TESTS: kmin [h/Mpc]:', kmin
         WRITE(*, *) 'HALOMOD_SPEED_TESTS: kmax [h/Mpc]:', kmax
         WRITE(*, *) 'HALOMOD_SPEED_TESTS: nk:', nk
         WRITE(*, *) 'HALOMOD_SPEED_TESTS: amin:', amin
         WRITE(*, *) 'HALOMOD_SPEED_TESTS: amax:', amax
         WRITE(*, *) 'HALOMOD_SPEED_TESTS: na:', na
         WRITE(*, *) 'HALOMOD_SPEED_TESTS: Total time [s]:', t
         WRITE(*, *)

         ! Write to file
         WRITE(unit, *) na, t

      END DO
      CLOSE(unit)

   END SUBROUTINE halomod_speed_tests

   SUBROUTINE compare_matter_vs_DMONLY_implementation(icosmo)

      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: icosmo
      REAL, ALLOCATABLE :: k(:), a(:)
      !REAL, ALLOCATABLE :: pow_li(:,:), pow_1h(:,:,:,:), pow_2h(:,:,:,:), pow_hm(:,:,:,:)
      REAL, ALLOCATABLE :: Pk(:, :, :, :)
      INTEGER :: ihm
      INTEGER:: ik, ia
      LOGICAL :: fail
      TYPE(cosmology) :: cosm
      !TYPE(halomod) :: hmod
      REAL, PARAMETER :: amin = 0.5
      REAL, PARAMETER :: amax = 1.0
      INTEGER, PARAMETER :: na = 2
      REAL, PARAMETER :: kmin = 1e-3
      REAL, PARAMETER :: kmax = 1e2
      INTEGER, PARAMETER :: nk = 101
      LOGICAL, PARAMETER :: verbose = .TRUE.
      INTEGER, PARAMETER :: fields(2) = [field_dmonly, field_matter]
      INTEGER, PARAMETER :: nf = 2
      !INTEGER, PARAMETER :: icosmo_default = 1 ! 60 - Boring but with 0.3 eV neutrinos
      INTEGER, PARAMETER :: ihm_default = 54 ! 54 - Matter masquerading as DMONLY
      REAL, PARAMETER :: eps = 1e-4

      CALL fill_array_log(kmin, kmax, k, nk)
      CALL fill_array(amin, amax, a, na)

      CALL assign_init_cosmology(icosmo, cosm, verbose)

      ihm = ihm_default
      !CALL assign_halomod(ihm, hmod, verbose)
      !CALL calculate_HMx(fields, nf, k, nk, a, na, pow_li, pow_2h, pow_1h, pow_hm, hmod, cosm, verbose)
      CALL calculate_HMx(fields, k, a, Pk, nk, na, nf, cosm, ihm)

      fail = .FALSE.
      DO ia = 1, na
         DO ik = 1, nk
            IF(.NOT. requal(Pk(1, 1, ik, ia), Pk(2, 2, ik, ia), eps)) THEN
               fail = .TRUE.
               WRITE(*, *) k(ik), a(ia), Pk(2, 2, ik, ia)/Pk(1, 1, ik, ia)
            END IF
         END DO
      END DO

      IF(fail) THEN
         WRITE(*, *)
         WRITE(*, *) 'Test failed'
      ELSE
         WRITE(*, *) 'Test passed'
      END IF
      WRITE(*, *)

   END SUBROUTINE compare_matter_vs_DMONLY_implementation

   SUBROUTINE compare_HMx_CAMB(imode, icosmo)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: imode
      INTEGER, INTENT(INOUT) :: icosmo
      REAL, ALLOCATABLE :: a(:), k(:), Pk_HMx(:, :), Pk_CAMB(:, :)
      REAL, ALLOCATABLE :: k_Tk(:), Tk(:,:)
      TYPE(cosmology) :: cosm
      REAL :: crap, max_error, error
      INTEGER :: i, j, ii, nfail
      INTEGER :: nk, nTk
      INTEGER :: HMx_version, HALOFIT_version, CAMB_nonlinear_version
      LOGICAL :: fail

      LOGICAL, PARAMETER :: verbose = .TRUE.
      REAL, PARAMETER :: amin = 1./3.
      REAL, PARAMETER :: amax = 1.0
      INTEGER, PARAMETER :: na = 16
      CHARACTER(len=256), PARAMETER :: outfile_HMx = 'data/power_HMcode_HMx.dat'
      CHARACTER(len=256), PARAMETER :: outfile_CAMB = 'data/power_HMcode_CAMB.dat'
      REAL, PARAMETER :: k_crap = 0.1      ! Wavenumber for first run to prevent function write [h/Mpc]
      REAL, PARAMETER :: a_crap = 1.0      ! Scale factor for first run to prevent function write
      REAL, PARAMETER :: eps = 0.01        ! Fractional error tolerance for the test
      REAL, PARAMETER :: kmin_test = 1e-2  ! Minimum wavenumber for the tests [h/Mpc]
      REAL, PARAMETER :: kmax_test = 1e1   ! Minimum wavenumber for the tests [h/Mpc]
      REAL, PARAMETER :: amin_test = 1./3. ! Minimum scale factor for test
      REAL, PARAMETER :: amax_test = 1.    ! Maximum scale factor for test
      INTEGER, PARAMETER :: iseed = 1      ! Seed for random number generator for tests
      INTEGER, PARAMETER :: ntest = 100    ! Number of tests to run
      LOGICAL, PARAMETER :: verbose_test = .FALSE.
      LOGICAL, PARAMETER :: stop_on_fail = .TRUE.

      IF (imode == 75) THEN
         ! HMcode 2016
         CAMB_nonlinear_version = CAMB_nonlinear_HMcode2016
         HMx_version = HMcode2016_CAMB
      ELSE IF (imode == 90) THEN
         ! HMcode 2015
         CAMB_nonlinear_version = CAMB_nonlinear_HMcode2015
         HMx_version = HMcode2015_CAMB
      ELSE IF (imode == 91) THEN
         ! HALOFIT Smith et al. (2003)
         CAMB_nonlinear_version = CAMB_nonlinear_HALOFIT_Smith
         HALOFIT_version = HALOFIT_Smith
      ELSE IF (imode == 92) THEN
         ! HALOFIT Bird et al. (2003)
         CAMB_nonlinear_version = CAMB_nonlinear_HALOFIT_Bird
         HALOFIT_version = HALOFIT_Bird
      ELSE IF (imode == 93) THEN
         ! HALOFIT Takahashi et al. (2012)
         CAMB_nonlinear_version = CAMB_nonlinear_HALOFIT_Takahashi
         HALOFIT_version = HALOFIT_Takahashi
      ELSE IF (imode == 97) THEN
         ! HALOFIT CAMB parameters
         CAMB_nonlinear_version = CAMB_nonlinear_HALOFIT_Takahashi
         HALOFIT_version = HALOFIT_CAMB
      ELSE IF (imode == 98) THEN
         ! HMcode 2020
         !CAMB_nonlinear_version = CAMB_nonlinear_HMcode2020
         HMx_version = HMcode2020
         STOP 'COMPARE_HMx_CAMB: Error, HMx 2020 needs to be implemented in CAMB'
      ELSE IF (imode == 103) THEN
         ! HMcode 2016 with neutrino fix
         CAMB_nonlinear_version = CAMB_nonlinear_HMcode2016_neutrinofix
         HMx_version = HMcode2016_neutrinofix_CAMB
      ELSE
         STOP 'COMPARE_HMx_CAMB: Error, something went wrong with imode'
      END IF

      ! Set the random number generator
      CALL RNG_set(iseed)  

      ! Loop over individual tests
      nfail = 0
      DO ii = 1, ntest

         ! Assigns the cosmological model
         CALL assign_cosmology(icosmo, cosm, verbose)
         CALL init_cosmology(cosm)
         CALL print_cosmology(cosm)

         ! Ensures power routines are called before calling CAMB below
         crap = plin(k_crap, a_crap, flag_power_total, cosm)

         ! Fill scale-factor arrays
         CALL fill_array(amin, amax, a, na)

         ! Calculate HMcode power via CAMB
         CALL get_CAMB_power(a, na, k, Pk_CAMB, nk, k_Tk, Tk, nTk, &
            non_linear=.TRUE., &
            halofit_version=CAMB_nonlinear_version, &
            cosm=cosm)

         ! Calculate HMcode power via HMx
         IF (ALLOCATED(Pk_HMx)) DEALLOCATE (Pk_HMx)
         ALLOCATE (Pk_HMx(nk, na))

         IF (is_in_array(imode, [75, 90, 98, 103])) THEN
            CALL calculate_HMcode(k, a, Pk_HMx, nk, na, cosm, HMx_version)
         ELSE IF (is_in_array(imode, [91, 92, 93, 97])) THEN
            CALL calculate_HALOFIT(k, a, Pk_HMx, nk, na, cosm, HALOFIT_version) 
         ELSE
            STOP 'COMPARE HMx CAMB: Error, imode not specified correctly'
         END IF

         ! Write data
         CALL write_power_a(k, a, Pk_HMx, nk, na, outfile_HMx, verbose=.FALSE.)
         CALL write_power_a(k, a, Pk_CAMB, nk, na, outfile_CAMB, verbose=.FALSE.)

         ! Now actually do tests
         fail = .FALSE.
         max_error = 0.
         DO j = 1, na
            DO i = 1, nk
               error = abs(-1.+Pk_HMx(i, j)/Pk_CAMB(i, j))
               IF ((k(i) >= kmin_test) .AND. (k(i) <= kmax_test) .AND. (a(j) >= amin_test) .AND. (a(j) <= amax_test)) THEN
                  IF (error > max_error) max_error = error
                  IF (error > eps) THEN
                     IF (verbose_test) THEN
                        WRITE (*, *) 'COMPARE_HMx_CAMB: Test failing'
                        WRITE (*, *) 'COMPARE_HMx_CAMB: k [h/Mpc]:', k(i)
                        WRITE (*, *) 'COMPARE_HMx_CAMB: a:', a(j)
                        WRITE (*, *) 'COMPARE_HMx_CAMB: Delta^2(k) [HMx]:', Pk_HMx(i, j)
                        WRITE (*, *) 'COMPARE_HMx_CAMB: Delta^2(k) [CAMB]:', Pk_CAMB(i, j)
                        WRITE (*, *) 'COMPARE_HMx_CAMB: Ratio:', Pk_HMx(i, j)/Pk_CAMB(i, j)
                        WRITE (*, *)
                     END IF
                     fail = .TRUE.
                  END IF
               END IF
            END DO
         END DO

         ! Write results to screen
         WRITE (*, *) 'COMPARE_HMx_CAMB: Seed:', iseed
         WRITE (*, *) 'COMPARE_HMx_CAMB: Test:', ii, 'of', ntest
         WRITE (*, *) 'COMPARE_HMx_CAMB: Minimum k considered in test [h/Mpc]:', kmin_test
         WRITE (*, *) 'COMPARE_HMx_CAMB: Maximum k considered in test [h/Mpc]:', kmax_test
         WRITE (*, *) 'COMPARE_HMx_CAMB: Minimum a considered in test:', amin_test
         WRITE (*, *) 'COMPARE_HMx_CAMB: Maximum a considered in test:', amax_test
         WRITE (*, *) 'COMPARE_HMx_CAMB: Maximum error:', max_error
         WRITE (*, *) 'COMPARE_HMx_CAMB: Tolerance:', eps
         IF (fail) THEN
            nfail = nfail+1
            WRITE (*, *) 'COMPARE_HMx_CAMB: Test failed'
            IF (stop_on_fail) STOP
         ELSE
            WRITE (*, *) 'COMPARE_HMx_CAMB: Test passed'
         END IF
         WRITE (*, *)

      END DO

      IF (stop_on_fail) THEN
         WRITE (*, *) 'COMPARE_HMx_CAMB: All tests passed'
      ELSE
         WRITE (*, *) 'COMPARE_HMx_CAMB: Number of failed tests:', nfail
         WRITE (*, *) 'COMPARE_HMx_CAMB: Total number of tests:', ntest
      END IF
      WRITE (*, *)

   END SUBROUTINE compare_HMx_CAMB

   SUBROUTINE non_linear_halo_bias_integrand(icosmo, ihm)

      ! Non-linear halo bias integrand
      USE table_integer
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: icosmo
      INTEGER, INTENT(INOUT) :: ihm
      REAL, ALLOCATABLE :: k(:), nu(:)
      REAL :: nu1, nu2, rv1, rv2, B_NL, I_NL
      INTEGER :: i, j, ik
      CHARACTER(len=256) :: outfile
      TYPE(cosmology) :: cosm
      TYPE(halomod) :: hmod
      LOGICAL, PARAMETER :: verbose = .TRUE.
      REAL, PARAMETER :: z = 0. ! Redshift
      REAL, PARAMETER :: kmin = 0.1
      REAL, PARAMETER :: kmax = 1.0
      INTEGER, PARAMETER :: nk = 10
      REAL, PARAMETER :: numin = 0.
      REAL, PARAMETER :: numax = 3.
      INTEGER, PARAMETER :: nnu = 100
      INTEGER, PARAMETER :: ifind = ifind_split
      INTEGER, PARAMETER :: iinterp = iinterp_Lagrange

      ! Set cosmology
      icosmo = 37 ! 37 - WMAP 5
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! Initiliasation for the halomodel calcualtion
      ihm = 49 ! 49 - Non-linear halo bias with Tinker
      CALL assign_halomod(ihm, hmod, verbose)
      CALL init_halomod(scale_factor_z(z), hmod, cosm, verbose)
      CALL print_halomod(hmod, cosm, verbose)

      ! Set k range and allocate
      CALL fill_array(kmin, kmax, k, nk)

      ! Set nu range and allocate
      CALL fill_pixels(numin, numax, nu, nnu)

      ! Needs to be called to prevent unit conflict below
      B_NL = BNL(k(1), nu(1), nu(1), 0., 0., hmod)

      ! Loop over k and two nu
      DO ik = 1, nk
         IF (ik == 1) outfile = 'data/bnl_k0.1.dat'
         IF (ik == 2) outfile = 'data/bnl_k0.2.dat'
         IF (ik == 3) outfile = 'data/bnl_k0.3.dat'
         IF (ik == 4) outfile = 'data/bnl_k0.4.dat'
         IF (ik == 5) outfile = 'data/bnl_k0.5.dat'
         IF (ik == 6) outfile = 'data/bnl_k0.6.dat'
         IF (ik == 7) outfile = 'data/bnl_k0.7.dat'
         IF (ik == 8) outfile = 'data/bnl_k0.8.dat'
         IF (ik == 9) outfile = 'data/bnl_k0.9.dat'
         IF (ik == 10) outfile = 'data/bnl_k1.0.dat'
         WRITE (*, *) 'HMx_DRIVER: Output file:', trim(outfile)
         WRITE (*, *) 'HMx_DRIVER: k [h/Mpc]:', k(ik)
         OPEN (7, file=outfile)
         DO i = 1, nnu
            DO j = 1, nnu
               nu1 = nu(i)
               nu2 = nu(j)
               rv1 = exp(find(nu1, hmod%nu, log(hmod%rv), hmod%n, iorder=3, ifind=ifind, iinterp=iinterp))
               rv2 = exp(find(nu2, hmod%nu, log(hmod%rv), hmod%n, iorder=3, ifind=ifind, iinterp=iinterp))
               B_NL = BNL(k(ik), nu1, nu2, rv1, rv2, hmod) ! Needed because otherwise function writes
               I_NL = B_NL*g_nu(nu1, hmod)*b_nu(nu1, hmod)*g_nu(nu2, hmod)*b_nu(nu2, hmod)
               WRITE (7, *) nu1, nu2, B_NL, I_NL
            END DO
         END DO
         CLOSE (7)
         WRITE (*, *) 'HMx_DRIVER: Done'
         WRITE (*, *)

      END DO

   END SUBROUTINE non_linear_halo_bias_integrand

   SUBROUTINE HMF_amplitude(icosmo, ihm)

      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: icosmo
      INTEGER, INTENT(INOUT) :: ihm
      TYPE(cosmology) :: cosm
      TYPE(halomod) :: hmod
      REAL, ALLOCATABLE :: k(:), pow_li(:), pow_1h(:), pow_2h(:), pow_hm(:)
      INTEGER :: i
      CHARACTER(len=256) :: outfile

      REAL, PARAMETER :: amp_min = 0.1
      REAL, PARAMETER :: amp_max = 2.
      INTEGER, PARAMETER :: n_amp = 16
      REAL, PARAMETER :: kmin = 1e-3
      REAL, PARAMETER :: kmax = 1e2
      INTEGER, PARAMETER :: nk = 128
      REAL, PARAMETER :: a = 1.
      REAL, PARAMETER :: mmin = 1e7
      REAL, PARAMETER :: mmax = 1e17
      INTEGER, PARAMETER :: field(1) = field_dmonly
      LOGICAL, PARAMETER :: verbose = .TRUE.
      CHARACTER(len=128), PARAMETER :: outfid = 'data/power_HMFamp_fid.dat'
      CHARACTER(len=128), PARAMETER :: fbase = 'data/power_HMFamp_'
      CHARACTER(len=128), PARAMETER :: fext = '.dat'

      ! 0 - Calculate halo model at a single z

      ! Set number of k points and k range (log spaced)
      CALL fill_array(log(kmin), log(kmax), k, nk)
      k = exp(k)

      ! Allocate arrays
      ALLOCATE (pow_li(nk), pow_2h(nk), pow_1h(nk), pow_hm(nk))

      ! Assigns the cosmological model
      icosmo = 1
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      DO i = 0, n_amp

         ! Initiliasation for the halomodel calcualtion
         ihm = 3
         CALL assign_halomod(ihm, hmod, verbose)

         IF (i == 0) THEN
            hmod%Amf = 1.
            outfile = outfid
         ELSE
            hmod%Amf = progression(amp_min, amp_max, i, n_amp)
            outfile = number_file(fbase, i, fext)
         END IF

         CALL init_halomod(a, hmod, cosm, verbose)
         CALL print_halomod(hmod, cosm, verbose)

         ! Do the halo-model calculation
         CALL calculate_HMx_a(field, 1, k, nk, pow_li, pow_2h, pow_1h, pow_hm, hmod, cosm, verbose)

         ! Write out the results
         CALL write_power(k, pow_li, pow_2h, pow_1h, pow_hm, nk, outfile, verbose)

      END DO

   END SUBROUTINE HMF_amplitude

   SUBROUTINE triad_tracers(version, i, ix, ix_name)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: version
      INTEGER, INTENT(IN) :: i
      INTEGER, INTENT(OUT) :: ix
      CHARACTER(len=256), INTENT(OUT) :: ix_name

      IF (version == 1) THEN
         CALL triad_1_tracers(i, ix, ix_name)
      ELSE IF (version == 2) THEN
         CALL triad_2_tracers(i, ix, ix_name)
      ELSE IF (version == 3) THEN
         CALL triad_3_tracers(i, ix, ix_name)
      ELSE IF (version == 4) THEN
         CALL triad_4_tracers(i, ix, ix_name)
      ELSE IF (version == 5) THEN
         CALL triad_5_tracers(i, ix, ix_name)
      ELSE IF (version == 6) THEN
         CALL triad_6_tracers(i, ix, ix_name)
      ELSE
         STOP 'TRIAD_TRACERS: Error, version of triad specified incorrectly'
      END IF

   END SUBROUTINE triad_tracers

   SUBROUTINE triad_1_tracers(i, ix, ix_name)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: i
      INTEGER, INTENT(OUT) :: ix
      CHARACTER(len=256), INTENT(OUT) :: ix_name

      IF (i == 1) THEN
         ix = tracer_KiDS
         ix_name = 'gal_z0.1-0.9'
      ELSE IF (i == 2) THEN
         ix = tracer_Compton_y
         ix_name = 'y'
      ELSE IF (i == 3) THEN
         ix = tracer_CMB_lensing
         ix_name = 'CMB'
      ELSE
         STOP 'TRIAD_1_TRACERS: Error, out of bounds'
      END IF

   END SUBROUTINE triad_1_tracers

   SUBROUTINE triad_2_tracers(i, ix, ix_name)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: i
      INTEGER, INTENT(OUT) :: ix
      CHARACTER(len=256), INTENT(OUT) :: ix_name

      IF (i == 1) THEN
         ix = tracer_KiDS_450
         ix_name = 'gal_z0.1-0.9'
      ELSE IF (i == 2) THEN
         ix = tracer_KiDS_450_fat_bin1
         ix_name = 'gal_z0.1-0.5'
      ELSE IF (i == 3) THEN
         ix = tracer_KiDS_450_fat_bin2
         ix_name = 'gal_z0.5-0.9'
      ELSE IF (i == 4) THEN
         ix = tracer_Compton_y
         ix_name = 'y'
      ELSE IF (i == 5) THEN
         ix = tracer_CMB_lensing
         ix_name = 'CMB'
      ELSE
         STOP 'TRIAD_2_TRACERS: Error, out of bounds'
      END IF

   END SUBROUTINE triad_2_tracers

   SUBROUTINE triad_3_tracers(i, ix, ix_name)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: i
      INTEGER, INTENT(OUT) :: ix
      CHARACTER(len=256), INTENT(OUT) :: ix_name

      IF (i == 1) THEN
         ix = tracer_KiDS_450
         ix_name = 'gal_z0.1-0.9'
      ELSE IF (i == 2) THEN
         ix = tracer_KiDS_450_fat_bin1
         ix_name = 'gal_z0.1-0.5'
      ELSE IF (i == 3) THEN
         ix = tracer_KiDS_450_fat_bin2
         ix_name = 'gal_z0.5-0.9'
      ELSE IF (i == 4) THEN
         ix = tracer_Compton_y
         ix_name = 'y'
      ELSE
         STOP 'TRIAD_3_TRACERS: Error, out of bounds'
      END IF

   END SUBROUTINE triad_3_tracers

   SUBROUTINE triad_4_tracers(i, ix, ix_name)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: i
      INTEGER, INTENT(OUT) :: ix
      CHARACTER(len=256), INTENT(OUT) :: ix_name

      IF (i == 1) THEN
         ix = tracer_lensing_z1p00
         ix_name = 'gal_z1.00'
      ELSE IF (i == 2) THEN
         ix = tracer_lensing_z0p75
         ix_name = 'gal_z0.75'
      ELSE IF (i == 3) THEN
         ix = tracer_lensing_z0p50
         ix_name = 'gal_z0.50'
      ELSE IF (i == 4) THEN
         ix = tracer_lensing_z0p25
         ix_name = 'gal_z0.25'
      ELSE IF (i == 5) THEN
         ix = tracer_Compton_y
         ix_name = 'y'
      ELSE
         STOP 'TRIAD_4_TRACERS: Error, out of bounds'
      END IF

   END SUBROUTINE triad_4_tracers

   SUBROUTINE triad_5_tracers(i, ix, ix_name)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: i
      INTEGER, INTENT(OUT) :: ix
      CHARACTER(len=256), INTENT(OUT) :: ix_name

      IF (i == 1) THEN
         ix = tracer_KiDS_450
         ix_name = 'gal_z0.1-0.9'
      ELSE IF (i == 2) THEN
         ix = tracer_KiDS_450_bin1
         ix_name = 'gal_z0.1-0.3'
      ELSE IF (i == 3) THEN
         ix = tracer_KiDS_450_bin2
         ix_name = 'gal_z0.3-0.5'
      ELSE IF (i == 4) THEN
         ix = tracer_KiDS_450_bin3
         ix_name = 'gal_z0.5-0.7'
      ELSE IF (i == 5) THEN
         ix = tracer_KiDS_450_bin4
         ix_name = 'gal_z0.7-0.9'
      ELSE IF (i == 6) THEN
         ix = tracer_Compton_y
         ix_name = 'y'
      ELSE
         STOP 'TRIAD_5_TRACERS: Error, out of bounds'
      END IF

   END SUBROUTINE triad_5_tracers

   SUBROUTINE triad_6_tracers(i, ix, ix_name)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: i
      INTEGER, INTENT(OUT) :: ix
      CHARACTER(len=256), INTENT(OUT) :: ix_name

      IF (i == 1) THEN
         ix = tracer_KiDS_450
         ix_name = 'gal_z0.1-0.9'
      ELSE IF (i == 2) THEN
         ix = tracer_KiDS_450_fat_bin1
         ix_name = 'gal_z0.1-0.5'
      ELSE IF (i == 3) THEN
         ix = tracer_KiDS_450_fat_bin2
         ix_name = 'gal_z0.5-0.9'
      ELSE IF (i == 4) THEN
         ix = tracer_KiDS_450_bin1
         ix_name = 'gal_z0.1-0.3'
      ELSE IF (i == 5) THEN
         ix = tracer_KiDS_450_bin2
         ix_name = 'gal_z0.3-0.5'
      ELSE IF (i == 6) THEN
         ix = tracer_KiDS_450_bin3
         ix_name = 'gal_z0.5-0.7'
      ELSE IF (i == 7) THEN
         ix = tracer_KiDS_450_bin4
         ix_name = 'gal_z0.7-0.9'
      ELSE IF (i == 8) THEN
         ix = tracer_lensing_z1p00
         ix_name = 'gal_z1.00'
      ELSE IF (i == 9) THEN
         ix = tracer_lensing_z0p75
         ix_name = 'gal_z0.75'
      ELSE IF (i == 10) THEN
         ix = tracer_lensing_z0p50
         ix_name = 'gal_z0.50'
      ELSE IF (i == 11) THEN
         ix = tracer_lensing_z0p25
         ix_name = 'gal_z0.25'
      ELSE IF (i == 12) THEN
         ix = tracer_Compton_y
         ix_name = 'y'
      ELSE IF (i == 13) THEN
         ix = tracer_CMB_lensing
         ix_name = 'CMB'
      ELSE
         STOP 'TRIAD_TRACERS: Error, out of bounds'
      END IF

   END SUBROUTINE triad_6_tracers

   SUBROUTINE read_k_values(infile, k, nk)

      ! Get k values from a data file, assumes they are the first column of a data file
      IMPLICIT NONE
      CHARACTER(len=256), INTENT(IN) :: infile ! Input file location
      REAL, ALLOCATABLE, INTENT(OUT) :: k(:) ! Output array of k values
      INTEGER, INTENT(OUT) :: nk ! Output number of k values
      INTEGER :: i

      ! Get the number of k values
      nk = file_length(infile, verbose=.FALSE.)

      ! Allocate the array in k
      ALLOCATE (k(nk))

      ! Read in the k values
      OPEN (7, file=infile)
      DO i = 1, nk
         READ (7, *) k(i)
      END DO
      CLOSE (7)

   END SUBROUTINE read_k_values

   SUBROUTINE random_baryon_parameters(hmod)

      IMPLICIT NONE
      TYPE(halomod), INTENT(INOUT) :: hmod

      REAL, PARAMETER :: alpha_min = 0.05
      REAL, PARAMETER :: alpha_max = 2.5

      REAL, PARAMETER :: eps_min = 0.5
      REAL, PARAMETER :: eps_max = 2.0

      REAL, PARAMETER :: Gamma_min = 1.05
      REAL, PARAMETER :: Gamma_max = 2.00

      REAL, PARAMETER :: M0_min = 10**(12.)
      REAL, PARAMETER :: M0_max = 10**(15.)

      REAL, PARAMETER :: Astar_min = 0.002
      REAL, PARAMETER :: Astar_max = 0.2

      REAL, PARAMETER :: Twhim_min = 10**(5.)
      REAL, PARAMETER :: Twhim_max = 10**(7.)

      REAL, PARAMETER :: cstar_min = 1.
      REAL, PARAMETER :: cstar_max = 1000.

      REAL, PARAMETER :: fcold_min = 1e-5
      REAL, PARAMETER :: fcold_max = 0.5

      REAL, PARAMETER :: mstar_min = 1e9
      REAL, PARAMETER :: mstar_max = 1e15

      REAL, PARAMETER :: sstar_min = 0.1
      REAL, PARAMETER :: sstar_max = 10.

      REAL, PARAMETER :: alphap_min = -1.
      REAL, PARAMETER :: alphap_max = 1.

      REAL, PARAMETER :: Gammap_min = -0.2
      REAL, PARAMETER :: Gammap_max = 0.2

      REAL, PARAMETER :: cstarp_min = -1.
      REAL, PARAMETER :: cstarp_max = 1.

      hmod%alpha = random_uniform(alpha_min, alpha_max)

      hmod%eps = random_uniform(log(eps_min), log(eps_max))
      hmod%eps = exp(hmod%eps)

      hmod%Gamma = random_uniform(Gamma_min, Gamma_max)

      hmod%M0 = random_uniform(log(M0_min), log(M0_max))
      hmod%M0 = exp(hmod%M0)

      hmod%Astar = random_uniform(Astar_min, Astar_max)

      hmod%Twhim = random_uniform(log(Twhim_min), log(Twhim_max))
      hmod%Twhim = exp(hmod%Twhim)

      hmod%cstar = random_uniform(log(cstar_min), log(cstar_max))
      hmod%cstar = exp(hmod%cstar)

      hmod%fcold = random_uniform(log(fcold_min), log(fcold_max))
      hmod%fcold = exp(hmod%fcold)

      hmod%mstar = random_uniform(log(mstar_min), log(mstar_max))
      hmod%mstar = exp(hmod%mstar)

      hmod%sstar = random_uniform(log(sstar_min), log(sstar_max))
      hmod%sstar = exp(hmod%sstar)

      hmod%alphap = random_uniform(alphap_min, alphap_max)

      hmod%Gammap = random_uniform(Gammap_min, Gammap_max)

      hmod%cstarp = random_uniform(cstarp_min, cstarp_max)

   END SUBROUTINE random_baryon_parameters

   SUBROUTINE triad_ell(ell, nl)

      IMPLICIT NONE
      REAL, ALLOCATABLE, INTENT(OUT) :: ell(:)
      INTEGER, INTENT(OUT) :: nl

      ! Triad 3 ell values
      nl = 10
      ALLOCATE (ell(nl))
      ell(1) = 1.018233764908628416e+02
      ell(2) = 1.553312629199899391e+02
      ell(3) = 2.323690199521311399e+02
      ell(4) = 3.400698045589385288e+02
      ell(5) = 4.735860422267907666e+02
      ell(6) = 6.606652213933027724e+02
      ell(7) = 9.337396960949569120e+02
      ell(8) = 1.314593569663487642e+03
      ell(9) = 1.844659111122426793e+03
      ell(10) = 2.590093228780283425e+03

   END SUBROUTINE triad_ell

   SUBROUTINE triad_ell_edges(ell, nl)

      IMPLICIT NONE
      REAL, ALLOCATABLE, INTENT(OUT) :: ell(:)
      INTEGER, INTENT(OUT) :: nl

      ! Triad 3 ell bin edge values
      nl = 11
      ALLOCATE (ell(nl))
      ell(1) = 1.000000000000000000e+02
      ell(2) = 1.405115826483646799e+02
      ell(3) = 1.974350485834819722e+02
      ell(4) = 2.774191114672182152e+02
      ell(5) = 3.898059840916188250e+02
      ell(6) = 5.477225575051662645e+02
      ell(7) = 7.696136340726083063e+02
      ell(8) = 1.081396297513015952e+03
      ell(9) = 1.519487052336355418e+03
      ell(10) = 2.135055305374795807e+03
      ell(11) = 3.000000000000001364e+03

   END SUBROUTINE triad_ell_edges

   SUBROUTINE bin_theory_ell(ell, Cl, nl, ell_bin, Cl_bin, n_bin)

      IMPLICIT NONE
      REAL, INTENT(IN) :: ell(nl)
      REAL, INTENT(IN) :: Cl(nl)
      INTEGER, INTENT(IN) :: nl
      REAL, INTENT(IN) :: ell_bin(n_bin+1)
      REAL, INTENT(OUT) :: Cl_bin(n_bin)
      INTEGER, INTENT(IN) :: n_bin
      INTEGER :: l
      INTEGER :: i
      REAL :: lmin, lmax, weight, weight_total

      ! Loop over new ell values
      DO i = 1, n_bin
         lmin = ell_bin(i)
         lmax = ell_bin(i+1)
         weight_total = 0.
         Cl_bin(i) = 0.
         DO l = 1, nl
            IF (ell(l) >= lmin .AND. ell(l) < lmax) THEN
               weight = ell(l)
               Cl_bin(i) = Cl_bin(i)+Cl(l)*weight
               weight_total = weight_total+weight
            END IF
         END DO
         Cl_bin(i) = Cl_bin(i)/weight_total
      END DO

   END SUBROUTINE bin_theory_ell

   SUBROUTINE add_highz_BAHAMAS(a, pow, nk, na, cosm)

      IMPLICIT NONE
      REAL, ALLOCATABLE, INTENT(INOUT) :: a(:)
      REAL, ALLOCATABLE, INTENT(INOUT) :: pow(:, :)
      INTEGER, INTENT(IN) :: nk
      INTEGER, INTENT(INOUT) :: na
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: a_save(na), pow_save(nk, na)
      INTEGER :: j
      REAL, PARAMETER :: z_new = 3.

      ! Save the original a and P(k,a) arrays
      a_save = a
      pow_save = pow

      ! Deallocate the original arrays
      DEALLOCATE (a, pow)

      ! Add an extra scale-factor at high-z
      na = na+1

      ! Reallocate original arrays with enough space for new z
      ALLOCATE (a(na), pow(nk, na))

      ! Add new scale factor
      a(1) = scale_factor_z(z_new)

      ! Add power as the previous highest z power times the squared growth ratio
      pow(:, 1) = pow_save(:, 1)*(grow(a(1), cosm)/grow(a_save(1), cosm))**2

      ! Fill in the rest of the arrays with the original values
      DO j = 1, na-1
         a(j+1) = a_save(j)
         pow(:, j+1) = pow_save(:, j)
      END DO

   END SUBROUTINE add_highz_BAHAMAS

   SUBROUTINE power_single(icosmo, ihm)

      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: icosmo
      INTEGER, INTENT(INOUT) :: ihm
      REAL, ALLOCATABLE :: k(:)
      REAL, ALLOCATABLE :: pow_li(:), pow_2h(:), pow_1h(:), pow_hm(:)
      CHARACTER(len=256) :: outfile
      TYPE(cosmology) :: cosm
      TYPE(halomod) :: hmod

      REAL, PARAMETER :: a = 1.
      REAL, PARAMETER :: kmin = 1e-3
      REAL, PARAMETER :: kmax = 1e2
      INTEGER, PARAMETER :: nk = 128
      INTEGER, PARAMETER :: field(1) = field_dmonly
      LOGICAL, PARAMETER :: verbose = .TRUE.

      ! Set number of k points and k range (log spaced)
      CALL fill_array(log(kmin), log(kmax), k, nk)
      k = exp(k)

      ! Assigns the cosmological model
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! Initiliasation for the halomodel calcualtion
      CALL assign_halomod(ihm, hmod, verbose)
      CALL init_halomod(a, hmod, cosm, verbose)
      CALL print_halomod(hmod, cosm, verbose)

      ! Allocate arrays
      ALLOCATE (pow_li(nk), pow_2h(nk), pow_1h(nk), pow_hm(nk))

      ! Do the halo-model calculation
      CALL calculate_HMx_a(field, 1, k, nk, pow_li, pow_2h, pow_1h, pow_hm, hmod, cosm, verbose)

      ! Write out the results
      outfile = 'data/power.dat'
      CALL write_power(k, pow_li, pow_2h, pow_1h, pow_hm, nk, outfile, verbose)

   END SUBROUTINE power_single

   SUBROUTINE power_single_comparison(imode, icosmo, ihm)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: imode
      INTEGER, INTENT(INOUT) :: icosmo
      INTEGER, INTENT(INOUT) :: ihm
      REAL, ALLOCATABLE :: k(:)
      REAL, ALLOCATABLE :: pow_li(:), pow_2h(:), pow_1h(:), pow_hm(:)
      CHARACTER(len=256) :: outfile
      INTEGER :: i
      INTEGER :: ihm_here, ihm_baseline
      INTEGER :: icosmo_here, icosmo_baseline
      TYPE(cosmology) :: cosm
      TYPE(halomod) :: hmod

      !INTEGER, PARAMETER :: icosmo_baseline = 1 ! 1 - Boring
      !INTEGER, PARAMETER :: ihm_baseline = 3    ! 3 - Seljak (2000)
      REAL, PARAMETER :: a = 1.
      REAL, PARAMETER :: kmin = 1e-3
      REAL, PARAMETER :: kmax = 1e2
      INTEGER, PARAMETER :: nk = 128
      INTEGER, PARAMETER :: field(1) = field_dmonly
      LOGICAL, PARAMETER :: verbose = .TRUE.

      ! Pick baseline halo model
      IF (imode == 77) THEN
         ihm_baseline = 3 ! 3 - Seljak (2000)
         icosmo_baseline = icosmo ! User choice
      ELSE IF (imode == 79) THEN
         ihm_baseline = 1 ! 1 - HMcode (2016)
         icosmo_baseline = icosmo ! User choice
      ELSE IF (imode == 81) THEN
         ihm_baseline = 3 ! 3 - Seljak (2000)
         icosmo_baseline = 1 ! 1 - Boring
      ELSE IF (imode == 83) THEN
         ihm_baseline = 1 ! 1 - HMcode (2016)
         icosmo_baseline = 1 ! 1 - Boring
      ELSE
         STOP 'POWER_SINGLE_COMPARISON: Error, imode not specified correctly'
      END IF

      ! Set number of k points and k range (log spaced)
      CALL fill_array(log(kmin), log(kmax), k, nk)
      k = exp(k)

      ! Allocate arrays
      ALLOCATE (pow_li(nk), pow_2h(nk), pow_1h(nk), pow_hm(nk))

      ! Loop over the baseline and the comparison model
      DO i = 1, 2

         ! Assigns the cosmological model
         IF (i == 1) THEN
            icosmo_here = icosmo
            ihm_here = ihm
            outfile = 'data/power.dat'
         ELSE IF (i == 2) THEN
            ! The comparison
            icosmo_here = icosmo_baseline
            ihm_here = ihm_baseline
            outfile = 'data/power_baseline.dat'
         ELSE
            STOP 'POWER_SINGLE_COMPARISON: A disaster occured'
         END IF
         !CALL assign_cosmology(icosmo_here, cosm, verbose)
         !CALL init_cosmology(cosm)
         !CALL print_cosmology(cosm)
         CALL assign_init_cosmology(icosmo_here, cosm, verbose)

         ! Initiliasation for the halomodel calcualtion
         CALL assign_halomod(ihm_here, hmod, verbose)
         CALL init_halomod(a, hmod, cosm, verbose)
         CALL print_halomod(hmod, cosm, verbose)

         ! Do the halo-model calculation
         CALL calculate_HMx_a(field, 1, k, nk, pow_li, pow_2h, pow_1h, pow_hm, hmod, cosm, verbose)

         ! Write out the results
         CALL write_power(k, pow_li, pow_2h, pow_1h, pow_hm, nk, outfile, verbose)

      END DO

   END SUBROUTINE power_single_comparison

   SUBROUTINE power_multiple(icosmo, ihm)

      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: icosmo
      INTEGER, INTENT(INOUT) :: ihm
      REAL, ALLOCATABLE :: k(:), a(:)
      REAL, ALLOCATABLE :: pow_li(:, :), pow_2h(:, :), pow_1h(:, :), pow_hm(:, :)
      CHARACTER(len=256) :: base
      TYPE(halomod) :: hmod
      TYPE(cosmology) :: cosm

      REAL, PARAMETER :: kmin = 1e-3
      REAL, PARAMETER :: kmax = 1e2
      INTEGER, PARAMETER :: nk = 128
      REAL, PARAMETER :: amin = 0.2
      REAL, PARAMETER :: amax = 1.
      INTEGER, PARAMETER :: na = 16
      LOGICAL, PARAMETER :: verbose = .TRUE.

      ! Set number of k points and k range (log spaced)
      CALL fill_array_log(kmin, kmax, k, nk)

      ! Set the number of scale factors and range
      CALL fill_array(amin, amax, a, na)

      ! Assigns the cosmological model
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! Assign the halo model
      CALL assign_halomod(ihm, hmod, verbose)
      CALL calculate_halomod_full(k, a, pow_li, pow_2h, pow_1h, pow_hm, nk, na, cosm, ihm)

      base = 'data/power'
      CALL write_power_a_multiple(k, a, pow_li, pow_2h, pow_1h, pow_hm, nk, na, base, verbose)

   END SUBROUTINE power_multiple

   SUBROUTINE power_multiple_comparison(imode, icosmo, ihm)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: imode
      INTEGER, INTENT(INOUT) :: icosmo
      INTEGER, INTENT(INOUT) :: ihm
      REAL, ALLOCATABLE :: k(:), a(:)
      !REAL, ALLOCATABLE :: pow_li(:, :), pow_2h(:, :, :, :), pow_1h(:, :, :, :), pow_hm(:, :, :, :)
      REAL, ALLOCATABLE :: pow_li(:, :), pow_2h(:, :), pow_1h(:, :), pow_hm(:, :)
      REAL :: zmin, zmax
      INTEGER :: icos, na, ncos
      INTEGER :: icosmo_here
      INTEGER :: ihm_here!, ihm_base, ihm_test
      CHARACTER(len=256) :: filebase!, filebase_test
      TYPE(halomod) :: hmod
      TYPE(cosmology) :: cosm

      !INTEGER, PARAMETER :: icosmo_baseline = 1 ! 1 - Boring
      !INTEGER, PARAMETER :: ihm_baseline = 3    ! 3 - Seljak (2000)
      REAL, PARAMETER :: kmin = 1e-3
      REAL, PARAMETER :: kmax = 1e2
      INTEGER, PARAMETER :: nk = 128
      LOGICAL, PARAMETER :: verbose = .TRUE.
      INTEGER, PARAMETER :: field(1) = field_dmonly

      ! Set number of k points and k range (log spaced)
      CALL fill_array(log(kmin), log(kmax), k, nk)
      k = exp(k)

      ! Set the number of redshifts and range (linearly spaced) and convert z -> a
      IF (imode == 85) THEN
         na = 6
         ALLOCATE(a(na))
         a(6) = 4.
         a(5) = 3.
         a(4) = 2.
         a(3) = 1.
         a(2) = 0.5
         a(1) = 0.
         a = scale_factor_z(a)
      ELSE
         zmin = 0.
         zmax = 4.
         na = 16
         CALL fill_array(zmin, zmax, a, na)
         a = scale_factor_z(a) ! Note that this is correct because 'a' here is actually 'z'
      END IF

      ! Fix the total number of cosmologies
      ncos = 1
      IF (imode == 85) ncos = 9

      ! Loop over cosmological models to compare
      DO icos = 0, ncos

         ! Assigns the cosmological model
         IF (imode == 85) THEN

            ihm_here = ihm ! 3 - Seljak (2000)             
            IF (icos == 0) THEN
               icosmo_here = 69 ! Planck-ish
               filebase = 'data/power_nobump'
            ELSE IF(icos == 1) THEN
               icosmo_here = 66 ! kbump = 0.05 h/Mpc; sigma = 1.0
               filebase = 'data/power_fatbump_k0p05'
            ELSE IF(icos == 2) THEN
               icosmo_here = 67 ! kbump = 0.10 h/Mpc; sigma = 1.0
               filebase = 'data/power_fatbump_k0p1'
            ELSE IF(icos == 3) THEN
               icosmo_here = 68 ! kbump = 1.00 h/Mpc; sigma = 1.0
               filebase = 'data/power_fatbump_k1'
            ELSE IF(icos == 4) THEN
               icosmo_here = 72 ! kbump = 0.05 h/Mpc; sigma = 0.1
               filebase = 'data/power_thinbump_k0p05'
            ELSE IF(icos == 5) THEN
               icosmo_here = 73 ! kbump = 0.10 h/Mpc; sigma = 0.1
               filebase = 'data/power_thinbump_k0p1'
            ELSE IF(icos == 6) THEN
               icosmo_here = 74 ! kbump = 1.00 h/Mpc; sigma = 0.1
               filebase = 'data/power_thinbump_k1'
            ELSE IF(icos == 7) THEN
               icosmo_here = 79 ! kbump = 0.05 h/Mpc; sigma = 0.3
               filebase = 'data/power_medbump_k0p05'
            ELSE IF(icos == 8) THEN
               icosmo_here = 80 ! kbump = 0.10 h/Mpc; sigma = 0.3
               filebase = 'data/power_medbump_k0p1'
            ELSE IF(icos == 9) THEN
               icosmo_here = 81 ! kbump = 1.00 h/Mpc; sigma = 0.3
               filebase = 'data/power_medbump_k1'
            ELSE
               STOP 'POWER_MULTIPLE_COMPARISON: Error, something went wrong with bumps'
            END IF

         ELSE

            IF (icos == 0) THEN

               ! The comparison
               IF (imode == 78) THEN
                  ihm_here = 3 ! 3 - Seljak (2000)        
                  icosmo_here = icosmo ! User choice        
               ELSE IF (imode == 80) THEN
                  ihm_here = 1 ! 1 - HMcode (2016)         
                  icosmo_here = icosmo ! User choice       
               ELSE IF (imode == 82) THEN
                  ihm_here = 3 ! 3 - Seljak (2000)
                  icosmo_here = 1 ! 1 - Boring        
               ELSE IF (imode == 84) THEN
                  ihm_here = 1 ! 1 - HMcode (2016)        
                  icosmo_here = 1 ! 1 - Boring        
               ELSE IF (imode == 85) THEN
                  ihm_here = ihm
                  icosmo_here = 69 ! Planck-ish
               END IF
               filebase = 'data/power_baseline'

            ELSE

               IF (imode == 78) THEN
                  ihm_here = ihm ! User choice
                  icosmo_here = icosmo ! User choice             
               ELSE IF (imode == 80) THEN
                  ihm_here = ihm ! User choice
                  icosmo_here = icosmo ! User choice
               ELSE IF (imode == 82) THEN
                  ihm_here = ihm ! User choice
                  icosmo_here = icosmo ! User choice
               ELSE IF (imode == 84) THEN
                  ihm_here = ihm ! User choice
                  icosmo_here = icosmo ! User choice
               END IF
               filebase = 'data/power'
         
            END IF

         END IF

         ! Assigns the cosmological model
         CALL assign_init_cosmology(icosmo_here, cosm, verbose)

         ! Assign the halo model
         CALL assign_halomod(ihm_here, hmod, verbose)

         !CALL calculate_HMx(field, 1, k, nk, a, na, pow_li, pow_2h, pow_1h, pow_hm, hmod, cosm, verbose)
         CALL calculate_halomod_full(k, a, pow_li, pow_2h, pow_1h, pow_hm, nk, na, cosm, ihm_here)

         !CALL write_power_a_multiple(k, a, pow_li, pow_2h(1, 1, :, :), pow_1h(1, 1, :, :), pow_hm(1, 1, :, :), nk, na, filebase, verbose)
         CALL write_power_a_multiple(k, a, pow_li, pow_2h, pow_1h, pow_hm, nk, na, filebase, verbose)

      END DO

   END SUBROUTINE power_multiple_comparison

   SUBROUTINE halo_stuff(imode, icosmo, ihm)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: imode
      INTEGER, INTENT(IN) :: icosmo
      INTEGER, INTENT(IN) :: ihm
      TYPE(halomod) :: hmod
      TYPE(cosmology) :: cosm
      REAL :: z
      INTEGER :: j, icosmo_here, ihm_here
      CHARACTER(len=256) :: dir

      LOGICAL, PARAMETER :: verbose = .TRUE.

      ! Halo diagnostics

      IF (imode == 86) THEN
         icosmo_here = 4
         ihm_here = 65
      ELSE
         icosmo_here = icosmo
         ihm_here = ihm
      END IF

      ! Assigns the cosmological model
      CALL assign_cosmology(icosmo_here, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! Loop over redshifts
      DO j = 1, 4

         ! Set the redshift
         IF (j == 1) z = 0.0
         IF (j == 2) z = 0.5
         IF (j == 3) z = 1.0
         IF (j == 4) z = 2.0

         ! Initiliasation for the halomodel calcualtion
         CALL assign_halomod(ihm_here, hmod, verbose)
         CALL init_halomod(scale_factor_z(z), hmod, cosm, verbose)
         CALL print_halomod(hmod, cosm, verbose)

         ! Runs the diagnostics
         dir = 'data'
         CALL halo_diagnostics(hmod, cosm, dir)
         CALL halo_definitions(hmod, cosm, dir)
         CALL halo_properties(hmod, dir)

      END DO

   END SUBROUTINE halo_stuff

   SUBROUTINE test_random_cosmologies(imode, icosmo, ihm)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: imode
      INTEGER, INTENT(INOUT) :: icosmo
      INTEGER, INTENT(INOUT) :: ihm
      REAL, ALLOCATABLE :: k(:)
      REAL :: z
      INTEGER :: i, ii, j1, j2
      INTEGER :: nf
      INTEGER, ALLOCATABLE :: fields(:)
      REAL, ALLOCATABLE :: pow_li(:), pow_2h(:, :, :), pow_1h(:, :, :), pow_hm(:, :, :)
      TYPE(cosmology) :: cosm
      TYPE(halomod) :: hmod

      REAL, PARAMETER :: zmin = 0.
      REAL, PARAMETER :: zmax = 2.
      REAL, PARAMETER :: kmin = 1e-3
      REAL, PARAMETER :: kmax = 1e2
      INTEGER, PARAMETER :: nk = 128
      INTEGER, PARAMETER :: n = 500
      LOGICAL, PARAMETER :: verbose = .TRUE.
      INTEGER, PARAMETER :: iseed_tests = 0

      ! Tests to ensure that the code does not crash
      !  4 - Random baryon parameters
      ! 72 - Random cosmological parameters

      ! Set the random number generator
      CALL RNG_set(iseed_tests)

      ! Set number of k points and k range (log spaced)
      !nk = 128
      !kmin = 1e-3
      !kmax = 1e2
      CALL fill_array(log(kmin), log(kmax), k, nk)
      k = exp(k)

      ! Set the number of fields
      IF (imode == 4) THEN
         nf = 5
      ELSE IF (imode == 72) THEN
         nf = 1
      ELSE
         STOP 'HMx_DRIVER: Error with random cosmology'
      END IF

      ! Allocate field array
      ALLOCATE (fields(nf))

      ! Set field types
      IF (imode == 4) THEN
         fields(1) = field_matter
         fields(2) = field_cdm
         fields(3) = field_gas
         fields(4) = field_stars
         fields(5) = field_electron_pressure
      ELSE IF (imode == 72) THEN
         fields = field_dmonly
      ELSE
         STOP 'HMx_DRIVER: Error with random cosmology'
      END IF

      ! Allocate arrays for power
      ALLOCATE (pow_li(nk), pow_2h(nf, nf, nk), pow_1h(nf, nf, nk), pow_hm(nf, nf, nk))

      ! Set the cosmological model
      IF (imode == 72) THEN
         CALL RNG_set(seed=0)
         !icosmo=39 ! 39 - Random cosmology
         icosmo = 40 ! Random CAMB cosmology
      END IF

      ! Loop forever
      DO ii = 1, n

         ! Set the cosmological model
         IF ((imode == 4 .AND. ii == 1) .OR. imode == 72) THEN
            CALL assign_cosmology(icosmo, cosm, verbose)
            CALL init_cosmology(cosm)
            CALL print_cosmology(cosm)
         END IF

         ! Random redshift in a range
         z = random_uniform(zmin, zmax)

         ! Initiliasation for the halomodel calcualtion
         CALL assign_halomod(ihm, hmod, verbose)
         IF (imode == 4) THEN
            CALL random_baryon_parameters(hmod)
         END IF
         CALL init_halomod(scale_factor_z(z), hmod, cosm, verbose)
         CALL print_halomod(hmod, cosm, verbose)

         CALL calculate_HMx_a(fields, nf, k, nk, pow_li, pow_2h, pow_1h, pow_hm, hmod, cosm, verbose=.FALSE.)

         ! Do the halo-model calculation for a range of halo types
         DO j1 = 1, nf
            DO j2 = 1, nf

               DO i = 1, nk
                  IF (ISNAN(pow_hm(j1, j2, i))) THEN
                     CALL print_halomod(hmod, cosm, verbose)
                     WRITE (*, *) 'HMx_DRIVER: Halo types:', fields(j1), fields(j2)
                     WRITE (*, *) 'HMx_DRIVER: k [h/Mpc]:', k(i)
                     STOP 'HMX_DRIVER: Error, NaN found in pow_full array'
                  END IF
               END DO

            END DO
         END DO

         WRITE (*, *) 'HMx_DRIVER: Done:', ii
         WRITE (*, *)

      END DO

   END SUBROUTINE test_random_cosmologies

   SUBROUTINE hydro_stuff(imode, icosmo, ihm)

      ! Make cross power spectra of all different components of haloes as well as pressure
      !  2 - Generic hydro
      !  5 - HMx 2020
      ! 32 - PAPER: Baseline hydro
      ! 52 - Generic hydro but with BAHAMAS k range
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: imode
      INTEGER, INTENT(INOUT) :: icosmo
      INTEGER, INTENT(INOUT) :: ihm
      REAL, ALLOCATABLE :: k(:), zs(:)
      REAL, ALLOCATABLE :: pow_li(:), pow_2h(:, :, :), pow_1h(:, :, :), pow_hm(:, :, :)
      REAL, ALLOCATABLE :: powd_li(:), powd_2h(:), powd_1h(:), powd_hm(:)
      INTEGER, ALLOCATABLE :: fields(:)
      REAL :: z
      INTEGER :: iz, j1, j2
      INTEGER :: iowl, field(1), nz, nowl
      CHARACTER(len=256) :: base, dir, ext, mid, outfile, name
      TYPE(cosmology) :: cosm
      TYPE(halomod) :: hmod

      REAL, PARAMETER :: kmin = 1e-3
      REAL, PARAMETER :: kmax = 1e2
      INTEGER :: nk = 128
      INTEGER, PARAMETER :: nf = 5
      LOGICAL, PARAMETER :: verbose = .TRUE.
      CHARACTER(len=256), PARAMETER :: infile_k = '/Users/Mead/Physics/BAHAMAS/power/M1024/DMONLY_nu0_L400N1024_WMAP9_snap32_all_all_power.dat'

      ! Set the redshift
      nz = 4
      ALLOCATE (zs(nz))
      zs(1) = 0.0
      zs(2) = 0.5
      zs(3) = 1.0
      zs(4) = 2.0

      IF (imode == 2 .OR. imode== 32 .OR. imode == 52) THEN
         ! Only do one 'model' here
         nowl = 1
      ELSE IF (imode == 5) THEN
         ! Do AGN_7p6, AGN_TUNED and AGN_8p0
         nowl = 3
      ELSE
         STOP 'HMx_DRIVER: Error, imode specified incorrectly'
      END IF

      ! Set number of k points and k range (log spaced)
      IF (imode == 2 .OR. imode == 32) THEN
         CALL fill_array(log(kmin), log(kmax), k, nk)
         k = exp(k)
      ELSE IF (imode == 5 .OR. imode == 52) THEN
         ! Get the k values from the simulation measured P(k)
         CALL read_k_values(infile_k, k, nk)
      ELSE
         STOP 'HMx_DRIVER: Error, imode specified incorrectly'
      END IF

      ! Field types
      ALLOCATE (fields(nf))
      fields(1) = field_matter
      fields(2) = field_cdm
      fields(3) = field_gas
      fields(4) = field_stars
      fields(5) = field_electron_pressure

      ! Allocate the arrays for P(k)
      ALLOCATE (powd_li(nk), powd_2h(nk), powd_1h(nk), powd_hm(nk))
      ALLOCATE (pow_li(nk), pow_2h(nf, nf, nk), pow_1h(nf, nf, nk), pow_hm(nf, nf, nk))

      DO iowl = 1, nowl

         DO iz = 1, nz

            z = zs(iz)

            !Assigns the cosmological model
            IF (imode == 2 .OR. imode == 32 .OR. imode == 52) THEN
               icosmo = 4
            END IF

            IF (imode == 32) THEN
               ihm = 65
            END IF

            IF (imode == 5) THEN
               ihm = 56 ! HMx2020 stars
               IF (iowl == 1) THEN
                  ! AGN 7.6
                  name = 'AGN_7p6_nu0'
                  icosmo = 4
               ELSE IF(iowl == 2) THEN
                  ! AGN TUNED
                  name = 'AGN_TUNED_nu0'
                  icosmo = 61
               ELSE IF (iowl == 3) THEN
                  ! AGN 8.0
                  name = 'AGN_8p0_nu0'
                  icosmo = 62
               END IF
            END IF

            ! Initialisation for cosmology
            CALL assign_cosmology(icosmo, cosm, verbose)
            CALL init_cosmology(cosm)
            CALL print_cosmology(cosm)
            
            ! Initiliasation for the halomodel calcualtion after variables changed
            CALL assign_halomod(ihm, hmod, verbose)
            CALL init_halomod(scale_factor_z(z), hmod, cosm, verbose)
            CALL print_halomod(hmod, cosm, verbose)

            ! Runs the diagnostics
            IF (imode == 2 .OR. imode == 52) THEN
               dir = 'data'
               CALL halo_diagnostics(hmod, cosm, dir)
               CALL halo_definitions(hmod, cosm, dir)
               CALL halo_properties(hmod, dir)
            END IF

            IF (imode == 2 .OR. imode == 32 .OR. imode == 52) THEN
               ! File base and extension
               IF (iz == 1) base = 'data/power_z0.0_'
               IF (iz == 2) base = 'data/power_z0.5_'
               IF (iz == 3) base = 'data/power_z1.0_'
               IF (iz == 4) base = 'data/power_z2.0_'
               mid = ''
               ext = '.dat'
            ! ELSE IF (imode == 15) THEN
            !    base = 'data/power_'//TRIM(fname)//'_'
            !    mid = ''
            !    ext = '.dat'
            ! ELSE IF (imode == 16) THEN
            !    IF (j == 1) base = 'data/power_'//TRIM(fname)//'_z0.0_'
            !    IF (j == 2) base = 'data/power_'//TRIM(fname)//'_z0.5_'
            !    IF (j == 3) base = 'data/power_'//TRIM(fname)//'_z1.0_'
            !    IF (j == 4) base = 'data/power_'//TRIM(fname)//'_z2.0_'
            !    mid = ''
            !    ext = '.dat'
            ELSE IF (imode == 5) THEN
               IF (iz == 1) base = 'data/power_'//trim(name)//'_z0.0_'
               IF (iz == 2) base = 'data/power_'//trim(name)//'_z0.5_'
               IF (iz == 3) base = 'data/power_'//trim(name)//'_z1.0_'
               IF (iz == 4) base = 'data/power_'//trim(name)//'_z2.0_'
               mid = ''
               ext = '.dat'
            END IF

            ! Dark-matter only
            IF (imode == 2 .OR. imode == 32 .OR. imode == 52) THEN
               IF (iz == 1) outfile = 'data/power_z0.0_11.dat'
               IF (iz == 2) outfile = 'data/power_z0.5_11.dat'
               IF (iz == 3) outfile = 'data/power_z1.0_11.dat'
               IF (iz == 4) outfile = 'data/power_z2.0_11.dat'
            ! ELSE IF (imode == 15) THEN
            !    outfile = 'data/power_DMONLY_11.dat'
            ! ELSE IF (imode == 16) THEN
            !    IF (j == 1) outfile = 'data/power_DMONLY_z0.0_11.dat'
            !    IF (j == 2) outfile = 'data/power_DMONLY_z0.5_11.dat'
            !    IF (j == 3) outfile = 'data/power_DMONLY_z1.0_11.dat'
            !    IF (j == 4) outfile = 'data/power_DMONLY_z2.0_11.dat'
            ELSE IF (imode == 5) THEN
               IF (iz == 1) outfile = 'data/power_DMONLY_z0.0_11.dat'
               IF (iz == 2) outfile = 'data/power_DMONLY_z0.5_11.dat'
               IF (iz == 3) outfile = 'data/power_DMONLY_z1.0_11.dat'
               IF (iz == 4) outfile = 'data/power_DMONLY_z2.0_11.dat'
            END IF

            ! Write some things to the screen
            field = field_dmonly
            WRITE (*, *) field(1), field(1), trim(outfile)

            ! Do the calculation for DMonly
            CALL calculate_HMx_a(field, 1, k, nk, powd_li, powd_2h, powd_1h, powd_hm, hmod, cosm, verbose)
            CALL write_power(k, powd_li, powd_2h, powd_1h, powd_hm, nk, outfile, verbose)

            ! Do the calculation for the rest of the fields
            CALL calculate_HMx_a(fields, nf, k, nk, pow_li, pow_2h, pow_1h, pow_hm, hmod, cosm, verbose)

            ! Loop over fields and write data
            DO j1 = 1, nf
               DO j2 = j1, nf

                  ! Fix output file and write to screen
                  outfile = number_file2(base, fields(j1), mid, fields(j2), ext)

                  ! Set the halo types and write to screen
                  WRITE (*, *) fields(j1), fields(j2), TRIM(outfile)

                  ! Write P(k)
                  CALL write_power(k, pow_li, pow_2h(j1, j2, :), pow_1h(j1, j2, :), pow_hm(j1, j2, :), nk, outfile, verbose=.FALSE.)

               END DO
            END DO

         END DO

      END DO

   END SUBROUTINE hydro_stuff

   SUBROUTINE hydro_across_temperature_range(imode, icosmo, ihm)

      ! Make cross power spectra of all different components of haloes as well as pressure
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: imode
      INTEGER, INTENT(INOUT) :: ihm
      INTEGER, INTENT(INOUT) :: icosmo
      REAL, ALLOCATABLE :: k(:), zs(:), logTs(:)
      REAL, ALLOCATABLE :: pow_li(:), pow_2h(:, :, :), pow_1h(:, :, :), pow_hm(:, :, :)
      REAL :: z, logT
      INTEGER :: iT, iz, j1, j2, nk
      CHARACTER(len=256) :: start, base, ext, mid, outfile
      TYPE(cosmology) :: cosm
      TYPE(halomod) :: hmod

      REAL, PARAMETER :: kmin = 1e-3
      REAL, PARAMETER :: kmax = 1e2
      REAL, PARAMETER :: logTmin = 7.5
      REAL, PARAMETER :: logTmax = 8.1
      !REAL, PARAMETER :: logTmin = 7.2
      !REAL, PARAMETER :: logTmax = 8.4
      INTEGER, PARAMETER :: nT = 7
      INTEGER, PARAMETER :: nf = 6
      INTEGER, PARAMETER :: fields(nf) = [field_dmonly, field_matter, field_cdm, field_gas, field_stars, field_electron_pressure]
      INTEGER, PARAMETER :: nz = 11
      LOGICAL, PARAMETER :: verbose = .TRUE.
      CHARACTER(len=256), PARAMETER :: infile_k = &
         '/Users/Mead/Physics/BAHAMAS/power/M1024/DMONLY_nu0_L400N1024_WMAP9_snap32_all_all_power.dat'

      ! Set models for paper plots
      IF (imode == 21) THEN
         ! matter-matter model
         icosmo = 4
         ihm = 60
         start = 'data/power_mm_T'
      ELSE IF(imode == 22) THEN
         ! matter-electron pressure model
         icosmo = 4
         ihm = 61
         start = 'data/power_me_T'
      ELSE
         start = 'data/power_T'
      END IF

      ! Set the redshifts
      ALLOCATE (zs(nz))
      zs(1) =  0.000
      zs(2) =  0.125
      zs(3) =  0.250
      zs(4) =  0.375
      zs(5) =  0.500
      zs(6) =  0.750
      zs(7) =  1.000
      zs(8) =  1.250
      zs(9) =  1.500
      zs(10) = 1.750
      zs(11) = 2.000

      CALL fill_array(logTmin, logTmax, logTs, nT)

      ! Get the k values from the simulation measured P(k)
      CALL read_k_values(infile_k, k, nk)

      ! Allocate the arrays for P(k)
      ALLOCATE (pow_li(nk), pow_2h(nf, nf, nk), pow_1h(nf, nf, nk), pow_hm(nf, nf, nk))

      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      CALL assign_halomod(ihm, hmod, verbose)

      DO iT = 1, nT

         logT = logTs(iT)
         cosm%Theat = 10**logT

         DO iz = 1, nz

            z = zs(iz)

            ! Initiliasation for the halomodel calculation
            CALL init_halomod(scale_factor_z(z), hmod, cosm, verbose)
            CALL print_halomod(hmod, cosm, verbose)

            base = trim(start)//trim(real_to_string(logT,1,1))//'_z'//trim(real_to_string(z,1,3))//'_'
            mid = ''
            ext = '.dat'

            ! Do the calculation for the rest of the fields
            CALL calculate_HMx_a(fields, nf, k, nk, pow_li, pow_2h, pow_1h, pow_hm, hmod, cosm, verbose)

            ! Loop over fields and write data
            DO j1 = 1, nf
               DO j2 = j1, nf

                  ! Fix output file and write to screen
                  outfile = number_file2(base, fields(j1), mid, fields(j2), ext)

                  ! Set the halo types and write to screen
                  WRITE (*, *) fields(j1), fields(j2), TRIM(outfile)

                  ! Write P(k)
                  CALL write_power(k, pow_li, pow_2h(j1, j2, :), pow_1h(j1, j2, :), pow_hm(j1, j2, :), nk, outfile, verbose=.FALSE.)

               END DO
            END DO

         END DO

      END DO

   END SUBROUTINE hydro_across_temperature_range

   SUBROUTINE hydro_suppression_as_a_function_of_cosmology(icosmo, ihm)

      ! Make cross power spectra of all different components of haloes as well as pressure
      IMPLICIT NONE
      !INTEGER, INTENT(IN) :: imode     
      INTEGER, INTENT(INOUT) :: icosmo
      INTEGER, INTENT(INOUT) :: ihm
      REAL, ALLOCATABLE :: k(:), zs(:)
      REAL, ALLOCATABLE :: pow_li(:), pow_2h(:, :, :), pow_1h(:, :, :), pow_hm(:, :, :)
      REAL, ALLOCATABLE :: suppression(:, :)
      REAL :: z
      INTEGER :: iz, nk, icos, ik
      TYPE(cosmology) :: cosm
      TYPE(halomod) :: hmod

      REAL, PARAMETER :: kmin = 1e-3
      REAL, PARAMETER :: kmax = 1e2
      REAL, PARAMETER :: logT = 7.8
      INTEGER, PARAMETER :: nf = 2
      INTEGER, PARAMETER :: fields(nf) = [field_dmonly, field_matter]
      INTEGER, PARAMETER :: nz = 1
      INTEGER, PARAMETER :: ncos = 2
      LOGICAL, PARAMETER :: verbose = .TRUE.
      CHARACTER(len=256), PARAMETER :: infile_k = &
         '/Users/Mead/Physics/BAHAMAS/power/M1024/DMONLY_nu0_L400N1024_WMAP9_snap32_all_all_power.dat'
      CHARACTER(len=256), PARAMETER :: outfile = 'data/hydro_suppression_cosmology.dat'

      ! Set the redshifts
      ALLOCATE (zs(nz))
      zs(1) = 0.

      ! Get the k values from the simulation measured P(k)
      CALL read_k_values(infile_k, k, nk)

      ! Allocate the arrays for P(k)
      ALLOCATE (pow_li(nk), pow_2h(nf, nf, nk), pow_1h(nf, nf, nk), pow_hm(nf, nf, nk))

      ALLOCATE(suppression(ncos, nk))

      DO icos = 1, ncos

         IF(icos == 1) THEN
            icosmo = 4 ! 4 - WMAP 9
         ELSE IF(icos == 2) THEN
            icosmo = 3 ! 3 - Planck 2013
         END IF

         CALL assign_cosmology(icosmo, cosm, verbose)
         CALL init_cosmology(cosm)
         CALL print_cosmology(cosm)

         CALL assign_halomod(ihm, hmod, verbose)
         hmod%Theat = 10**logT

         DO iz = 1, nz

            z = zs(iz)

            ! Initiliasation for the halomodel calculation
            CALL init_halomod(scale_factor_z(z), hmod, cosm, verbose)
            CALL print_halomod(hmod, cosm, verbose)

            ! Do the calculation for the rest of the fields
            CALL calculate_HMx_a(fields, nf, k, nk, pow_li, pow_2h, pow_1h, pow_hm, hmod, cosm, verbose)

            DO ik = 1, nk
               suppression(icos, ik) = pow_hm(2, 2, ik)/pow_hm(1, 1, ik)            
            END DO

         END DO

      END DO

      WRITE(*, *) 'HMx: Writing suppression data: ', trim(outfile)
      OPEN(7, file=outfile)    
      DO ik = 1, nk
         WRITE(7, *) k(ik), (suppression(icos, ik), icos = 1, ncos)
      END DO  
      CLOSE(7) 
      WRITE(*, *) 'HMx: Done'
      WRITE(*, *)

   END SUBROUTINE hydro_suppression_as_a_function_of_cosmology

   SUBROUTINE general_projection(imode, icosmo, ihm)

      ! TODO: Split this up into smaller chunks
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: imode
      INTEGER, INTENT(INOUT) :: icosmo
      INTEGER, INTENT(INOUT) :: ihm
      REAL, ALLOCATABLE :: k(:), a(:), pow_ka(:, :)
      REAL, ALLOCATABLE :: ell(:), Cl(:, :, :)
      REAL, ALLOCATABLE :: theta(:), xi(:, :)
      REAL, ALLOCATABLE :: pow_li(:, :), pow_2h(:, :, :, :), pow_1h(:, :, :, :), pow_hm(:, :, :, :)
      REAL :: spam
      REAL :: m1, m2, a1, a2, z1, z2, r1, r2
      REAL :: lmin, lmax, zmin, zmax
      INTEGER :: i, j
      INTEGER :: nl, ix(2), ip(2), nz
      INTEGER, ALLOCATABLE :: ibessel(:)
      CHARACTER(len=256) :: base, ext, mid, outfile
      TYPE(projection) :: proj(2)
      TYPE(cosmology) :: cosm
      TYPE(halomod) :: hmod

      REAL, PARAMETER :: kmin = 1e-3
      REAL, PARAMETER :: kmax = 1e1
      INTEGER, PARAMETER :: nk = 128
      REAL, PARAMETER :: amin = 0.2 ! Problems with one-halo term if amin is less than 0.1
      REAL, PARAMETER :: amax = 1.
      INTEGER, PARAMETER :: na = 16
      REAL, PARAMETER :: lmax_xi = 1e5
      REAL, PARAMETER :: sig8min = 0.7
      REAL, PARAMETER :: sig8max = 0.9
      REAL, PARAMETER :: thmin = 0.01
      REAL, PARAMETER :: thmax = 10.
      INTEGER, PARAMETER :: nth = 128
      INTEGER, PARAMETER :: ncos = 5 ! I may have changed this number inadvertantly
      INTEGER, PARAMETER :: nb = 3
      CHARACTER(len=256), PARAMETER :: dir = 'data/' ! Output directory
      LOGICAL, PARAMETER :: icumulative = .TRUE. ! Do cumlative distributions for breakdown
      LOGICAL, PARAMETER :: ifull = .FALSE.      ! Do only full halo model C(l), xi(theta) calculations (quicker, no breakdown ...)
      LOGICAL, PARAMETER :: verbose = .TRUE.

      ! General stuff for various 2D projections
      !  7 - Do general angular cross correlation
      !  8 - Angular cross correlation as a function of cosmology
      !  9 - Breakdown angular correlations in halo mass
      ! 10 - Breakdown angular correlations in redshift
      ! 11 - Breakdown angular correlations in halo radius
      ! 37 - CFHTLenS angular correlations
      ! 42 - PAPER: breakdown lensing-lensing per ell
      ! 43 - PAPER: breakdown lensing-y per ell
      ! 47 - Make CMB-lensing data to compare with CAMB
      ! 63 - Breakdown y-y per ell
      ! 64 - CIB

      ! Set the fields
      ix = -1
      !IF (imode == 37) ix = tracer_CFHTLenS_VanWaerbeke2013 ! CFHTLenS autospectrum
      IF (imode == 37) ix = tracer_CFHTLenS_Kilbinger2013 ! CFHTLenS autospectrum
      IF (imode == 42) ix = tracer_KiDS_450 ! KiDS-450 autospectrum
      IF (imode == 43) THEN
         ix(1) = tracer_KiDS_450  ! KiDS-450 z = 0.1->0.9
         ix(2) = tracer_Compton_y ! Compton y
      END IF
      IF (imode == 47) ix = tracer_CMB_lensing ! CMB lensing autospectrum
      IF (imode == 63) ix = tracer_Compton_y   ! tSZ Compton y parameter
      IF (imode == 64) ix = tracer_CIB_545     ! CIB at 545 GHz
      DO i = 1, 2
         CALL set_field_for_xpow(ix(i), ip(i))
      END DO
      IF (imode == 37 .OR. imode == 47) ip = field_dmonly ! Set DMONLY haloes

      ! Assign the cosmological model
      IF (imode == 42 .OR. imode == 43 .OR. imode == 63) icosmo = 4 ! 4 - WMAP9
      IF (imode == 47) icosmo = 26 ! 26 - Boring with CAMB linear spectrum
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)

      ! Set the k range
      CALL fill_array(log(kmin), log(kmax), k, nk)
      k = exp(k)

      ! Set the a range
      CALL fill_array(amin, amax, a, na)

      ! Need to call 'comoving_distance' at least once first so as to stop
      ! a write trying to happen while printing to screen
      spam = comoving_distance(1., cosm) ! CARE: This needs to be called before the write-to-screen below
      WRITE (*, *) 'HMx_DRIVER: lmax:', lmax
      WRITE (*, *) '======================================================='
      WRITE (*, *) '            a             z     r [Mpc/h]     k [h/MPc]'
      WRITE (*, *) '======================================================='
      DO i = 1, na
         WRITE (*, fmt='(4F14.4)') a(i), redshift_a(a(i)), comoving_distance(a(i), cosm), k_ell(lmax, a(i), cosm)
      END DO
      WRITE (*, *) '======================================================='
      WRITE (*, *)

      ! Set the ell range and fill l array
      IF (imode == 47) THEN
         lmin = 2.
         lmax = 10000.
         nl = nint(lmax)-nint(lmin)+1
         CALL fill_array(lmin, lmax, ell, nl)
      ELSE
         lmin = 1e0
         lmax = 1e4 ! Problems if this is pushed up to 10^5
         IF (imode == 11 .OR. imode == 37) lmax = 1e5 ! Need higher lmax for the correlation functions
         nl = 128
         CALL fill_array(log(lmin), log(lmax), ell, nl)
         ell = exp(ell)
      END IF

      ! Allocate arrays for l and C(l)
      ALLOCATE (Cl(2, 2, nl))

      ! Set the angular arrays in degrees and allocate arrays for theta and xi(theta)
      !thmin = 0.01 ! Minium theta [deg]
      !thmax = 10.  ! Maxium theta [deg]
      !INTEGER, PARAMETER :: nth = 128
      CALL fill_array(log(thmin), log(thmax), theta, nth)
      theta = exp(theta)
      !nb = 3 ! Number of correlation functions to compute
      ALLOCATE (xi(nb, nth), ibessel(nb))
      ibessel(1) = 0 ! J0 (xi plus for lensing)
      ibessel(2) = 2 ! J2 (standard angular correlation function)
      ibessel(3) = 4 ! J4 (xi minus for lensing)

      WRITE (*, *) 'HMx_DRIVER: Correlation function stuff'
      WRITE (*, *) 'HMx_DRIVER: Minium theta [deg]:', thmin
      WRITE (*, *) 'HMx_DRIVER: Minium theta [deg]:', thmax
      WRITE (*, *) 'HMx_DRIVER: Number of theta:', nth
      WRITE (*, *) 'HMx_DRIVER: lmax for xi:', lmax_xi
      WRITE (*, *)

      ! Write to screen
      WRITE (*, *) 'HMx_DRIVER: Cross-correlation information'
      WRITE (*, *) 'HMx_DRIVER: output directiory: ', TRIM(dir)
      WRITE (*, *) 'HMx_DRIVER: Profile type 1: ', TRIM(halo_type(ip(1)))
      WRITE (*, *) 'HMx_DRIVER: Profile type 2: ', TRIM(halo_type(ip(2)))
      WRITE (*, *) 'HMx_DRIVER: cross-correlation type 1: ', TRIM(xcorr_type(ix(1)))
      WRITE (*, *) 'HMx_DRIVER: cross-correlation type 2: ', TRIM(xcorr_type(ix(2)))
      WRITE (*, *) 'HMx_DRIVER: P(k) minimum k [h/Mpc]:', REAL(kmin)
      WRITE (*, *) 'HMx_DRIVER: P(k) maximum k [h/Mpc]:', REAL(kmax)
      WRITE (*, *) 'HMx_DRIVER: minimum a:', REAL(amin)
      WRITE (*, *) 'HMx_DRIVER: maximum a:', REAL(amax)
      WRITE (*, *) 'HMx_DRIVER: minimum ell:', REAL(lmin)
      WRITE (*, *) 'HMx_DRIVER: maximum ell:', REAL(lmax)
      WRITE (*, *)

      IF (imode == 7 .OR. &
          imode == 11 .OR. &
          imode == 37 .OR. &
          imode == 42 .OR. &
          imode == 43 .OR. &
          imode == 47 .OR. &
          imode == 63 .OR. &
          imode == 64) THEN

         ! Write cosmology to screen
         CALL print_cosmology(cosm)

         ! Write out diagnostics
         IF (imode == 47) ihm = 1 ! 1 - HMcode (2016)
         IF (imode == 63 .OR. imode == 64) ihm = 3  ! 3 - Standard halo model
         IF (imode == 42 .OR. imode == 43) ihm = 61 ! 55 - HMx2020
         !CALL assign_halomod(ihm, hmod, verbose)
         !CALL calculate_HMx(ip, 2, k, nk, a, na, pows_li, pows_2h, pows_1h, pows_hm, hmod, cosm, verbose)
         CALL calculate_HMx_full(ip, k, a, pow_li, pow_2h, pow_1h, pow_hm, 2, nk, na, cosm, ihm)

         base = TRIM(dir)//'power'
         CALL write_power_a_multiple(k, a, pow_li, pow_2h(1, 2, :, :), pow_1h(1, 2, :, :), pow_hm(1, 2, :, :), nk, na, base, verbose)

         ! Fill out the projection kernels
         CALL fill_projection_kernels(ix, proj, cosm)

         ! Set the distance range for the Limber integral
         !r1=100.
         !r2=proj%rs
         r1 = 0.
         r2 = maxdist(proj)

         ! Write information to screen
         WRITE (*, *) 'HMx_DRIVER: Computing C(l)'
         WRITE (*, *) 'HMx_DRIVER: r min [Mpc/h]:', r1
         WRITE (*, *) 'HMx_DRIVER: r max [Mpc/h]:', r2
         WRITE (*, *) 'HMx_DRIVER: ell min for C(l):', REAL(ell(1))
         WRITE (*, *) 'HMx_DRIVER: ell max for C(l):', REAL(ell(nl))
         WRITE (*, *) 'HMx_DRIVER: number of ell:', nl
         WRITE (*, *) 'HMx_DRIVER: lower limit of Limber integral [Mpc/h]:', REAL(r1)
         WRITE (*, *) 'HMx_DRIVER: upper limit of Limber integral [Mpc/h]:', REAL(r2)
         WRITE (*, *)

         ALLOCATE (pow_ka(nk, na))

         ! Loop over all types of C(l) to create
         DO j = 4, 1, -1

            IF (ifull .AND. (j .NE. 4)) CYCLE

            ! Write information to screen
            IF (j == 1) THEN
               WRITE (*, *) 'HMx_DRIVER: Doing linear'
               outfile = TRIM(dir)//'cl_linear.dat'
               IF (imode == 37) outfile = TRIM(dir)//'CFHTLenS_cl_linear.dat'
               IF (imode == 47) outfile = TRIM(dir)//'CMBlensing_cl_linear.dat'
               pow_ka = pow_li
            ELSE IF (j == 2) THEN
               WRITE (*, *) 'HMx_DRIVER: Doing 2-halo'
               outfile = TRIM(dir)//'cl_2h.dat'
               IF (imode == 37) outfile = TRIM(dir)//'CFHTLenS_cl_2h.dat'
               IF (imode == 47) outfile = TRIM(dir)//'CMBlensing_cl_2h.dat'
               pow_ka = pow_2h(1, 2, :, :)
            ELSE IF (j == 3) THEN
               WRITE (*, *) 'HMx_DRIVER: Doing 1-halo'
               outfile = TRIM(dir)//'cl_1h.dat'
               IF (imode == 37) outfile = TRIM(dir)//'CFHTLenS_cl_1h.dat'
               IF (imode == 47) outfile = TRIM(dir)//'CMBlensing_cl_1h.dat'
               pow_ka = pow_1h(1, 2, :, :)
            ELSE IF (j == 4) THEN
               WRITE (*, *) 'HMx_DRIVER: Doing full'
               outfile = TRIM(dir)//'cl_hm.dat'
               IF (imode == 37) outfile = TRIM(dir)//'CFHTLenS_cl_hm.dat'
               IF (imode == 47) outfile = TRIM(dir)//'CMBlensing_cl_hm.dat'
               pow_ka = pow_hm(1, 2, :, :)
            ELSE
               STOP 'HMx_DRIVER: Something went wrong'
            END IF

            WRITE (*, *) 'HMx_DRIVER: Output: ', TRIM(outfile)

            ! Actually calculate the C(ell)
            CALL calculate_Cl(r1, r2, ell, Cl(1, 2, :), nl, k, a, pow_ka, nk, na, proj, cosm)
            CALL write_Cl(ell, Cl(1, 2, :), nl, outfile, verbose)

            IF (j == 4 .AND. (imode == 7 .OR. imode == 11 .OR. imode == 42 .OR. imode == 43 .OR. imode == 63)) THEN
               base = 'data/Cl_contribution_ell_'
               ext = '.dat'
               CALL Cl_contribution_ell(r1, r2, k, a, pow_ka, nk, na, proj, cosm, base, ext)
               IF (imode == 42 .OR. imode == 43 .OR. imode == 63) EXIT
            END IF

            IF (imode == 11 .OR. imode == 37) THEN

               ! Set xi output files
               IF (j == 1) THEN
                  outfile = TRIM(dir)//'xi_linear.dat'
                  IF (imode == 37) outfile = TRIM(dir)//'CFHTLenS_xi_linear.dat'
               ELSE IF (j == 2) THEN
                  outfile = TRIM(dir)//'xi_2h.dat'
                  IF (imode == 37) outfile = TRIM(dir)//'CFHTLenS_xi_2h.dat'
               ELSE IF (j == 3) THEN
                  outfile = TRIM(dir)//'xi_1h.dat'
                  IF (imode == 37) outfile = TRIM(dir)//'CFHTLenS_xi_1h.dat'
               ELSE IF (j == 4) THEN
                  outfile = TRIM(dir)//'xi_hm.dat'
                  IF (imode == 37) outfile = TRIM(dir)//'CFHTLenS_xi_hm.dat'
               ELSE
                  STOP 'HMx_DRIVER: Something went wrong'
               END IF

               WRITE (*, *) 'HMx_DRIVER: Output: ', TRIM(outfile)

               ! Actually calculate the xi(theta)
               CALL calculate_angular_xi(ibessel, nb, theta, xi, nth, ell, Cl(1, 2, :), nl, NINT(lmax_xi))
               CALL write_angular_xi(theta, xi, nb, nth, outfile)

            END IF

         END DO
         WRITE (*, *) 'HMx_DRIVER: Done'
         WRITE (*, *)

      ELSE IF (imode == 8) THEN

         ! Allocate array for power
         ALLOCATE (pow_ka(nk, na))

         ! Set range in sigma_8
         !sig8min = 0.7
         !sig8max = 0.9
         !ncos = 5 ! I may have changed this number inadvertantly

         ! Loop over cosmology
         DO i = 1, ncos

            cosm%sig8 = progression(sig8min, sig8max, i, ncos)
            CALL init_cosmology(cosm)
            CALL print_cosmology(cosm)

            !CALL assign_halomod(ihm, hmod, verbose)
            !CALL calculate_HMx(ip, 2, k, nk, a, na, pows_li, pows_2h, pows_1h, pows_hm, hmod, cosm, verbose)
            CALL calculate_HMx_full(ip, k, a, pow_li, pow_2h, pow_1h, pow_hm, 2, nk, na, cosm, ihm)

            ! Fill out the projection kernels
            CALL fill_projection_kernels(ix, proj, cosm)

            ! Write to screen
            WRITE (*, *) 'HMx_DRIVER: Computing C(l)'
            WRITE (*, *) 'HMx_DRIVER: ell min:', REAL(ell(1))
            WRITE (*, *) 'HMx_DRIVER: ell max:', REAL(ell(nl))
            WRITE (*, *) 'HMx_DRIVER: number of ell:', nl
            WRITE (*, *)

            ! Loop over all types of C(l) to create
            !dir = 'data/'
            base = TRIM(dir)//'cosmology_'
            DO j = 1, 4

               IF (j == 1) THEN
                  WRITE (*, *) 'HMx_DRIVER: Doing C(l) linear'
                  ext = '_cl_linear.dat'
                  pow_ka = pow_li
               ELSE IF (j == 2) THEN
                  WRITE (*, *) 'HMx_DRIVER: Doing C(l) 2-halo'
                  ext = '_cl_2h.dat'
                  pow_ka = pow_2h(1, 2, :, :)
               ELSE IF (j == 3) THEN
                  WRITE (*, *) 'HMx_DRIVER: Doing C(l) 1-halo'
                  ext = '_cl_1h.dat'
                  pow_ka = pow_1h(1, 2, :, :)
               ELSE IF (j == 4) THEN
                  WRITE (*, *) 'HMx_DRIVER: Doing C(l) full'
                  ext = '_cl_hm.dat'
                  pow_ka = pow_hm(1, 2, :, :)
               END IF
               outfile = number_file(base, i, ext)
               WRITE (*, *) 'HMx_DRIVER: Output: ', TRIM(outfile)

               ! Actually calculate the C(l)
               CALL calculate_Cl(0., maxdist(proj), ell, Cl, nl, k, a, pow_ka, nk, na, proj, cosm)
               CALL write_Cl(ell, Cl, nl, outfile, verbose)

            END DO
            WRITE (*, *) 'HMx_DRIVER: Done'
            WRITE (*, *)

         END DO

      ELSE IF (imode == 9) THEN

         ! Breakdown cross-correlation in terms of mass

         ! Print to screen
         CALL print_cosmology(cosm)

         ! Fill out the projection kernels
         CALL fill_projection_kernels(ix, proj, cosm)

         ! Allocate arrays for power
         ALLOCATE (pow_li(nk, na), pow_2h(2, 2, nk, na), pow_1h(2, 2, nk, na), pow_hm(2, 2, nk, na))
         ALLOCATE (pow_ka(nk, na))

         DO i = 0, 6
            IF (icumulative .EQV. .FALSE.) THEN
               ! Set the mass intervals
               IF (i == 0) THEN
                  m1 = 1e7
                  m2 = 1e17
               ELSE IF (i == 1) THEN
                  m1 = 1e7
                  m2 = 1e11
               ELSE IF (i == 2) THEN
                  m1 = 1e11
                  m2 = 1e12
               ELSE IF (i == 3) THEN
                  m1 = 1e12
                  m2 = 1e13
               ELSE IF (i == 4) THEN
                  m1 = 1e13
                  m2 = 1e14
               ELSE IF (i == 5) THEN
                  m1 = 1e14
                  m2 = 1e15
               ELSE IF (i == 6) THEN
                  m1 = 1e15
                  m2 = 1e16
               END IF
            ELSE
               ! Set the mass intervals
               IF (i == 0) THEN
                  m1 = 1e7
                  m2 = 1e17
               ELSE IF (i == 1) THEN
                  m1 = 1e7
                  m2 = 1e11
               ELSE IF (i == 2) THEN
                  m1 = 1e7
                  m2 = 1e12
               ELSE IF (i == 3) THEN
                  m1 = 1e7
                  m2 = 1e13
               ELSE IF (i == 4) THEN
                  m1 = 1e7
                  m2 = 1e14
               ELSE IF (i == 5) THEN
                  m1 = 1e7
                  m2 = 1e15
               ELSE IF (i == 6) THEN
                  m1 = 1e7
                  m2 = 1e16
               END IF
            END IF

            ! Set the code to not 'correct' the two-halo power for missing
            ! mass when doing the calcultion binned in halo mass
            IF ((icumulative .EQV. .TRUE.) .AND. i > 0) hmod%ip2h = 0
            STOP 'ERROR: I am almost certain that hmod%ip2h=0 above needs to be changed to hmod%ip2h_corr=1'

            WRITE (*, fmt='(A16)') 'HMx_DRIVER: Mass range'
            WRITE (*, fmt='(A16,I5)') 'HMx_DRIVER: Iteration:', i
            WRITE (*, fmt='(A21,2ES15.7)') 'HMx_DRIVER: M_min [Msun/h]:', m1
            WRITE (*, fmt='(A21,2ES15.7)') 'HMx_DRIVER: M_max [Msun/h]:', m2
            WRITE (*, *)

            !Loop over redshifts
            DO j = 1, na

               !Initiliasation for the halomodel calcualtion
               CALL assign_halomod(ihm, hmod, verbose)
               hmod%mmin = m1
               hmod%mmax = m2
               CALL init_halomod(a(j), hmod, cosm, verbose)
               CALL print_halomod(hmod, cosm, verbose)
               CALL calculate_HMx_a(ip, 2, k, nk, pow_li(:, j), pow_2h(:, :, :, j), pow_1h(:, :, :, j), pow_hm(:, :, :, j), &
                                    hmod, cosm, verbose)

               !Write progress to screen
               IF (j == 1) THEN
                  WRITE (*, fmt='(A5,A7)') 'i', 'a'
                  WRITE (*, fmt='(A13)') '   ============'
               END IF
               WRITE (*, fmt='(I5,F8.3)') j, a(j)

            END DO
            WRITE (*, fmt='(A13)') '   ============'
            WRITE (*, *)

            !dir = 'data/'
            IF (i == 0) THEN
               outfile = TRIM(dir)//'power'
            ELSE
               base = TRIM(dir)//'mass_'
               mid = '_'
               ext = '_power'
               outfile = number_file2(base, NINT(log10(m1)), mid, NINT(log10(m2)), ext)
            END IF
            WRITE (*, *) 'HMx_DRIVER: File: ', TRIM(outfile)
            !CALL write_power_a(k,a,powa,nk,na,output)

            !Loop over all types of C(l) to create
            base = TRIM(dir)//'mass_'
            mid = '_'
            DO j = 1, 4

               !Skip the 1-halo C(l) because it takes ages (2017/02/06)
               IF (j == 3) CYCLE

               !Set output files
               IF (j == 1) THEN
                  outfile = TRIM(dir)//'cl_linear.dat'
                  ext = '_cl_linear.dat'
                  pow_ka = pow_li
               ELSE IF (j == 2) THEN
                  outfile = TRIM(dir)//'cl_2h.dat'
                  ext = '_cl_2h.dat'
                  pow_ka = pow_2h(1, 2, :, :)
               ELSE IF (j == 3) THEN
                  outfile = TRIM(dir)//'cl_1h.dat'
                  ext = '_cl_1h.dat'
                  pow_ka = pow_1h(1, 2, :, :)
               ELSE IF (j == 4) THEN
                  outfile = TRIM(dir)//'cl_hm.dat'
                  ext = '_cl_hm.dat'
                  pow_ka = pow_hm(1, 2, :, :)
               END IF

               IF (i > 0) outfile = number_file2(base, NINT(log10(m1)), mid, NINT(log10(m2)), ext)

               WRITE (*, *) 'HMx_DRIVER: File: ', TRIM(outfile)

               CALL calculate_Cl(0., maxdist(proj), ell, Cl, nl, k, a, pow_ka, nk, na, proj, cosm)
               CALL write_Cl(ell, Cl, nl, outfile, verbose)

            END DO
            WRITE (*, *) 'HMx_DRIVER: Done'
            WRITE (*, *)

         END DO

      ELSE IF (imode == 10) THEN

         ! Break down cross-correlation in terms of redshift

         ! Print to screen
         CALL print_cosmology(cosm)

         !CALL assign_halomod(ihm, hmod, verbose)
         !CALL calculate_HMx(ip, 2, k, nk, a, na, pow_li, pow_2h, pow_1h, pow_hm, hmod, cosm, verbose)
         CALL calculate_HMx_full(ip, k, a, pow_li, pow_2h, pow_1h, pow_hm, 2, nk, na, cosm, ihm)

         ! Allocate array for power
         ALLOCATE (pow_ka(nk, na))

         ! Fill out the projection kernels
         CALL fill_projection_kernels(ix, proj, cosm)

         ! Write to screen
         WRITE (*, *) 'HMx_DRIVER: Computing C(l)'
         WRITE (*, *) 'HMx_DRIVER: ell min:', REAL(ell(1))
         WRITE (*, *) 'HMx_DRIVER: ell max:', REAL(ell(nl))
         WRITE (*, *) 'HMx_DRIVER: number of ell:', nl
         WRITE (*, *)

         ! Set range of z
         zmin = 0.
         zmax = 1.
         nz = 8

         DO i = 0, nz

            IF (i == 0) THEN
               r1 = 0.
               r2 = maxdist(proj)
            ELSE
               IF (icumulative .EQV. .FALSE.) THEN
                  z1 = progression(zmin, zmax, i, nz)
               ELSE
                  z1 = zmin
               END IF
               z2 = zmin+(zmax-zmin)*float(i)/float(nz)
               a1 = scale_factor_z(z1)
               a2 = scale_factor_z(z2)
               r1 = comoving_distance(a1, cosm)
               r2 = comoving_distance(a2, cosm)
            END IF

            WRITE (*, *) 'HMx_DRIVER:', i
            IF (i > 0) THEN
               WRITE (*, *) 'HMx_DRIVER: z1:', REAL(z1)
               WRITE (*, *) 'HMx_DRIVER: z2:', REAL(z2)
            END IF
            WRITE (*, *) 'HMx_DRIVER: r1 [Mpc/h]:', REAL(r1)
            WRITE (*, *) 'HMx_DRIVER: r2 [Mpc/h]:', REAL(r2)

            ! Loop over all types of C(l) to create
            !dir = 'data/'
            base = TRIM(dir)//'redshift_'
            mid = '_'
            DO j = 1, 4

               !Set output files
               IF (j == 1) THEN
                  ext = '_cl_linear.dat'
                  outfile = TRIM(dir)//'cl_linear.dat'
                  pow_ka = pow_li
               ELSE IF (j == 2) THEN
                  ext = '_cl_2h.dat'
                  outfile = TRIM(dir)//'cl_2h.dat'
                  pow_ka = pow_2h(1, 2, :, :)
               ELSE IF (j == 3) THEN
                  ext = '_cl_1h.dat'
                  outfile = TRIM(dir)//'cl_1h.dat'
                  pow_ka = pow_1h(1, 2, :, :)
               ELSE IF (j == 4) THEN
                  ext = '_cl_hm.dat'
                  outfile = TRIM(dir)//'cl_hm.dat'
                  pow_ka = pow_hm(1, 2, :, :)
               END IF

               IF (i > 0 .AND. (icumulative .EQV. .FALSE.)) THEN
                  outfile = number_file2(base, i-1, mid, i, ext)
               ELSE IF (i > 0 .AND. icumulative) THEN
                  outfile = number_file2(base, 0, mid, i, ext)
               END IF
               WRITE (*, *) 'HMx_DRIVER: Output: ', TRIM(outfile)

               !This crashes for the low r2 values for some reason
               !Only a problem if lmax ~ 10^5
               !STOP 'This crashes for the low r2 values for high ell for some reason - should debug'
               CALL calculate_Cl(r1, r2, ell, Cl, nl, k, a, pow_ka, nk, na, proj, cosm)
               CALL write_Cl(ell, Cl, nl, outfile, verbose)

            END DO
            WRITE (*, *)

         END DO

         WRITE (*, *) 'HMx_DRIVER: Done'
         WRITE (*, *)

      ELSE

         STOP 'HMx_DRIVER: Error, you have specified the mode incorrectly'

      END IF

   END SUBROUTINE

   SUBROUTINE triad_stuff(imode, icosmo, ihm)

      ! Triad stuff
      ! 12 - Triad (T_5 cross correlations)
      ! 38 - AGN triad for all BAHAMAS feedback scenarios
      ! 39 - AGN triad for all BAHAMAS feedback scenarios and ell
      ! 44 - Triad for paper (fixed WMAP9 and feedback)
      ! 55 - Triad for all BAHAMAS feedback scenarios
      ! 56 - PAPER: Triad for all BAHAMAS feedback scenarios and ell
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: imode
      INTEGER, INTENT(INOUT) :: icosmo
      INTEGER, INTENT(INOUT) :: ihm
      CHARACTER(len=256), ALLOCATABLE :: ixx_names(:)
      CHARACTER(len=256) :: base, outfile
      INTEGER, ALLOCATABLE :: ixx(:)
      REAL, ALLOCATABLE :: ell(:), Cl(:, :, :)
      INTEGER :: i, j, ii, jj
      INTEGER :: nfeed
      TYPE(cosmology) :: cosm
      !TYPE(halomod) :: hmod

      REAL, PARAMETER :: lmin = 100.
      REAL, PARAMETER :: lmax = 4000.
      INTEGER :: nl = 64
      INTEGER, PARAMETER :: triad_version = 6
      INTEGER, PARAMETER :: nt = 13 !5
      CHARACTER(len=256), PARAMETER :: dir = 'data' ! Directory for data output
      LOGICAL, PARAMETER :: verbose = .TRUE.

      IF (imode == 12 .OR. imode == 44) THEN
         nfeed = 1
      ELSE IF (imode == 38 .OR. imode == 39 .OR. imode == 55 .OR. imode == 56) THEN
         nfeed = 3
      ELSE
         STOP 'HMx_DRIVER: Error, something fucked up'
      END IF

      ! Assign the cosmological model
      IF (imode == 38 .OR. imode == 39 .OR. imode == 44 .OR. imode == 55 .OR. imode == 56) icosmo = 4 ! WMAP 9 - BAHAMAS
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! Set the ell range
      IF (imode == 12 .OR. imode == 44 .OR. imode == 55) THEN
         ! Allocate arrays for l and C(l)
         CALL fill_array(log(lmin), log(lmax), ell, nl)
         ell = exp(ell)
      ELSE IF (imode == 39 .OR. imode == 56) THEN
         ! Triad 3 ell values
         CALL triad_ell(ell, nl)
      ELSE
         STOP 'HMx_DRIVER: Error, something fucked up, ell not set'
      END IF

      WRITE (*, *) 'HMx_DRIVER: Cross-correlation information'
      WRITE (*, *) 'HMx_DRIVER: output directiory: ', TRIM(dir)
      WRITE (*, *) 'HMx_DRIVER: minimum ell:', REAL(lmin)
      WRITE (*, *) 'HMx_DRIVER: maximum ell:', REAL(lmax)
      WRITE (*, *) 'HMx_DRIVER: number of ell:', nl
      WRITE (*, *)

      ! Set tracers
      ALLOCATE (ixx(nt), ixx_names(nt))
      DO i = 1, nt
         CALL triad_tracers(triad_version, i, ixx(i), ixx_names(i))
      END DO

      ! Allocate array for Cl
      ALLOCATE (Cl(nl, nt, nt))

      ! Loop over feedback scenarios
      DO j = 1, nfeed

         IF (imode == 44) ihm = 18 ! AGN tuned

         IF (imode == 38 .OR. imode == 39) THEN
            IF (j == 1) ihm = 18 ! AGN tuned
            IF (j == 2) ihm = 17 ! AGN low
            IF (j == 3) ihm = 19 ! AGN high
         END IF

         IF (imode == 55 .OR. imode == 56) THEN
            IF (j == 1) ihm = 39 ! AGN tuned
            IF (j == 2) ihm = 38 ! AGN low
            IF (j == 3) ihm = 40 ! AGN high
         END IF     

         IF (j == 1) base = 'triad_Cl_AGN_TUNED_nu0'
         IF (j == 2) base = 'triad_Cl_AGN_7p6_nu0'
         IF (j == 3) base = 'triad_Cl_AGN_8p0_nu0'

         ! Assign the halo
         !CALL assign_halomod(ihm, hmod, verbose)
         !CALL xpow_halomod(ixx, nt, ell, Cl, nl, hmod, cosm, verbose=.TRUE.)
         CALL xpow_halomod(ixx, nt, ell, Cl, nl, cosm, ihm, verbose=.TRUE.)

         ! Write data
         DO ii = 1, nt
            DO jj = 1, nt
               outfile = TRIM(dir)//'/'//TRIM(base)//'_'//TRIM(ixx_names(ii))//'-'//TRIM(ixx_names(jj))//'.dat'
               CALL write_Cl(ell, Cl(:, ii, jj), nl, outfile, verbose=.TRUE.)
            END DO
         END DO

      END DO

   END SUBROUTINE triad_stuff

   SUBROUTINE cross_correlation_coefficient(icosmo, ihm)

      ! Calculate a cross-correlation coefficient
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: icosmo
      INTEGER, INTENT(INOUT) :: ihm
      CHARACTER(len=256) :: outfile
      INTEGER :: ip(2), ix(2)
      REAL, ALLOCATABLE :: ell(:), Cl(:, :, :)
      INTEGER :: i, j
      TYPE(cosmology) :: cosm
      !TYPE(halomod) :: hmod

      CHARACTER(len=256), PARAMETER :: dir = 'data' ! Directory for outputs
      REAL, PARAMETER :: lmin = 1e0
      REAL, PARAMETER :: lmax = 1e4 ! Strange errors and crashes if this is increased to 10^5
      INTEGER, PARAMETER :: nl = 64
      LOGICAL, PARAMETER :: verbose = .TRUE.

      ! Assign the cosmology
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! Set the ell range and allocate arrays for l and C(l)
      CALL fill_array(log(lmin), log(lmax), ell, nl)
      ell = exp(ell)
      ALLOCATE (Cl(nl, 2, 2))

      ! Set the necessary fields
      ix = -1
      DO i = 1, 2
         CALL set_field_for_xpow(ix(i), ip(i))
      END DO

      ! Do the cross correlation
      ! Assign the halo model
      !CALL assign_halomod(ihm, hmod, verbose)
      !CALL xpow_halomod(ix, 2, ell, Cl, nl, hmod, cosm, verbose)
      CALL xpow_halomod(ix, 2, ell, Cl, nl, cosm, ihm, verbose)

      DO i = 1, 2
         DO j = i, 2

            ! Set output files
            IF (i == 1 .AND. j == 1) THEN
               outfile = TRIM(dir)//'/cl_first.dat'
            ELSE IF (i == 2 .AND. j == 2) THEN
               outfile = TRIM(dir)//'/cl_second.dat'
            ELSE IF (i == 1 .AND. j == 2) THEN
               outfile = TRIM(dir)//'/cl_hm.dat'
            END IF

            ! Write data
            CALL write_Cl(ell, Cl(:, i, j), nl, outfile, verbose)

         END DO
      END DO

   END SUBROUTINE cross_correlation_coefficient

   SUBROUTINE baryon_parameter_variations(imode, icosmo, ihm)

      ! Make power spectra variations as a function of baryon parameter variations
      ! 14 - General version
      ! 33 - Paper version
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: imode
      INTEGER, INTENT(INOUT) :: icosmo
      INTEGER, INTENT(INOUT) :: ihm
      REAL, ALLOCATABLE :: k(:)
      REAL :: param_min, param_max, param, param_neat
      LOGICAL :: ilog
      INTEGER :: ipa
      INTEGER :: i, j, j1, j2
      REAL :: mass
      REAL, ALLOCATABLE :: pow_li(:), pow_2h(:, :, :), pow_1h(:, :, :), pow_hm(:, :, :)
      REAL, ALLOCATABLE :: powd_li(:), powd_2h(:), powd_1h(:), powd_hm(:)
      INTEGER, ALLOCATABLE :: fields(:)
      CHARACTER(len=256) :: base, ext, mid, outfile
      TYPE(cosmology) :: cosm
      TYPE(halomod) :: hmod

      REAL, PARAMETER :: kmin = 1e-3
      REAL, PARAMETER :: kmax = 1e1
      INTEGER, PARAMETER :: nk = 128
      REAL, PARAMETER :: z = 0.0 ! Set the redshift
      INTEGER, PARAMETER :: m = 9 ! Number of values to try for each parameter
      INTEGER, PARAMETER :: nf = 5
      INTEGER, PARAMETER :: field(1) = field_dmonly
      !CHARACTER(len=256), PARAMETER :: outfile = 'data/DMONLY.dat'
      LOGICAL, PARAMETER :: verbose = .TRUE.

      ! Set number of k points and k range (log spaced)
      CALL fill_array(log(kmin), log(kmax), k, nk)
      k = exp(k)

      ! Set the fields
      ALLOCATE (fields(nf))
      fields(1) = field_matter
      fields(2) = field_cdm
      fields(3) = field_gas
      fields(4) = field_stars
      fields(5) = field_electron_pressure

      ! Allocate arrays for the power spectra
      ALLOCATE (powd_li(nk), powd_2h(nk), powd_1h(nk), powd_hm(nk))
      ALLOCATE (pow_li(nk), pow_2h(nf, nf, nk), pow_1h(nf, nf, nk), pow_hm(nf, nf, nk))

      ! Assigns the cosmological model
      IF (imode == 33) icosmo = 4
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! Initiliasation for the halo-model calcualtion
      IF (imode == 33) ihm = 55
      CALL assign_halomod(ihm, hmod, verbose)
      CALL init_halomod(scale_factor_z(z), hmod, cosm, verbose)
      CALL print_halomod(hmod, cosm, verbose)

      ! DMONLY
      outfile = 'data/DMONLY.dat'
      CALL calculate_HMx_a(field, 1, k, nk, powd_li, powd_2h, powd_1h, powd_hm, hmod, cosm, verbose)
      CALL write_power(k, powd_li, powd_2h, powd_1h, powd_hm, nk, outfile, verbose)

      ! Prevents warning
      ilog = .FALSE.

      ! Loop over parameters
      DO ipa = 1, param_n

         ! Set maximum and minimum parameter values and linear or log range
         IF (ipa == param_alpha) THEN
            ! 1 - alpha - virial temperature pre factor
            param_min = hmod%alpha-0.3
            param_max = hmod%alpha+0.3
            ilog = .FALSE.
         ELSE IF (ipa == param_eps) THEN
            ! 2 - epsilon - concentration change due to gas
            param_min = hmod%eps-0.5
            param_max = hmod%eps+0.5
            ilog = .FALSE.
         ELSE IF (ipa == param_gamma) THEN
            ! 3 - Gamma - KS polytropic index
            param_min = hmod%Gamma-0.05
            param_max = hmod%Gamma+0.05
            ilog = .FALSE.
         ELSE IF (ipa == param_M0) THEN
            ! 4 - M0 - bound gas transition in Msun/h
            param_min = hmod%M0/10.
            param_max = hmod%M0*10.
            ilog = .TRUE.
         ELSE IF (ipa == param_Astar) THEN
            ! 5 - A* - Stellar mass fraction
            param_min = hmod%Astar-0.01
            param_max = hmod%Astar+0.01
            ilog = .FALSE.
         ELSE IF (ipa == param_Twhim) THEN
            ! 6 - WHIM temperature in K
            param_min = hmod%Twhim/10.
            param_max = hmod%Twhim*10.
            ilog = .TRUE.
         ELSE IF (ipa == param_cstar) THEN
            ! 7 - c* - stellar concentation
            param_min = hmod%cstar/2.
            param_max = hmod%cstar*2.
            ilog = .TRUE.
         ELSE IF (ipa == param_fcold) THEN
            ! 8 - fcold - fraction of bound gas that is cold
            param_min = hmod%fcold
            param_max = hmod%fcold+0.25
            ilog = .FALSE.
         ELSE IF (ipa == param_mstar) THEN
            ! 9 - M* - peak of stellar halo mass filling
            param_min = hmod%mstar/sqrt(10.)
            param_max = hmod%mstar*sqrt(10.)
            ilog = .TRUE.
         ELSE IF (ipa == param_sstar) THEN
            ! 10 - sigma* - width of stellar halo mass filling distribution
            param_min = hmod%sstar-0.5
            param_max = hmod%sstar+0.5
            ilog = .FALSE.
         ELSE IF (ipa == param_alphap) THEN
            ! 11 - alphap - alpha mass index
            param_min = hmod%alphap-0.1
            param_max = hmod%alphap+0.1
            ilog = .FALSE.
         ELSE IF (ipa == param_Gammap) THEN
            ! 12 - Gammap - Gamma mass index
            param_min = hmod%Gammap-0.01
            param_max = hmod%Gammap+0.01
            ilog = .FALSE.
         ELSE IF (ipa == param_cstarp) THEN
            ! 13 - cstarp - mass power law for c*
            param_min = hmod%cstarp-0.1
            param_max = hmod%cstarp+0.1
            ilog = .FALSE.
         ELSE IF (ipa == param_fhot) THEN
            ! 14 - fhot - hot gas fraction
            param_min = hmod%fcold
            param_max = hmod%fcold+0.25
            ilog = .FALSE.
         ELSE IF (ipa == param_alphaz) THEN
            ! 15 - alphaz - z power-law for alpha
            param_min = hmod%alphaz-0.1
            param_max = hmod%alphaz+0.1
            ilog = .FALSE.
         ELSE IF (ipa == param_gammaz) THEN
            ! 16 - Gammaz - z power-law for Gamma
            param_min = hmod%gammaz-0.1
            param_max = hmod%gammaz+0.1
            ilog = .FALSE.
         ELSE IF (ipa == param_M0z) THEN
            ! 17 - M0z - z power-law for M0
            param_min = hmod%M0z-0.1
            param_max = hmod%M0z+0.1
            ilog = .FALSE.
         ELSE IF (ipa == param_Astarz) THEN
            ! 18 - Astarz - z power-law for A*
            param_min = hmod%Astarz-0.1
            param_max = hmod%Astarz+0.1
            ilog = .FALSE.
         ELSE IF (ipa == param_Twhimz) THEN
            ! 19 - Twhimz - z power-law for Twhim
            param_min = hmod%Twhimz-0.1
            param_max = hmod%Twhimz+0.1
            ilog = .FALSE.
         ELSE IF (ipa == param_eta) THEN
            ! 20 - eta - power-law for central galaxy mass fraction
            IF(imode == 33) THEN
               ! In this case the default eta is -0.3
               param_min = hmod%eta-0.2
               param_max = hmod%eta+0.2
            ELSE
               ! Otherwise the default value for eta is 0.
               param_min = hmod%eta-0.5
               param_max = hmod%eta
            END IF
            ilog = .FALSE.
         ELSE IF (ipa == param_epsz) THEN
            ! 33 - epsz - z power-law for epsilon
            param_min = hmod%epsz-0.1
            param_max = hmod%epsz+0.1
            ilog = .FALSE.
         ELSE IF (ipa == param_beta) THEN
            ! 34 - beta - temperature of hot gas
            param_min = hmod%beta-0.3
            param_max = hmod%beta+0.3
            ilog = .FALSE.
         ELSE IF (ipa == param_betap) THEN
            ! 35 - betap - mass power law for beta
            param_min = hmod%betap-0.1
            param_max = hmod%betap+0.1
            ilog = .FALSE.
         ELSE IF (ipa == param_betaz) THEN
            ! 36 - betaz - z power law for beta
            param_min = hmod%betaz-0.1
            param_max = hmod%betaz+0.1
            ilog = .FALSE.
         ELSE
            CYCLE
         END IF

         ! Loop over parameter values
         DO i = 1, m

            ! Set the parameter value that is being varied
            IF (ilog) THEN
               param = progression_log(param_min, param_max, i, m)
            ELSE
               param = progression(param_min, param_max, i, m)
            END IF

            CALL assign_halomod(ihm, hmod, verbose)

            !IF (hmod%HMx_mode == 1 .OR. hmod%HMx_mode == 2 .OR. hmod%HMx_mode == 3 .OR. hmod%HMx_mode == 5 .OR. hmod%HMx_mode == 6) THEN
            IF (is_in_array(hmod%HMx_mode, [1, 2, 3, 4, 5, 6])) THEN  
               IF (ipa == param_alpha)  hmod%alpha = param
               IF (ipa == param_eps)    hmod%eps = param
               IF (ipa == param_Gamma)  hmod%Gamma = param
               IF (ipa == param_M0)     hmod%M0 = param
               IF (ipa == param_Astar)  hmod%Astar = param
               IF (ipa == param_Twhim)  hmod%Twhim = param
               IF (ipa == param_cstar)  hmod%cstar = param
               IF (ipa == param_fcold)  hmod%fcold = param
               IF (ipa == param_mstar)  hmod%mstar = param
               IF (ipa == param_sstar)  hmod%sstar = param
               IF (ipa == param_alphap) hmod%alphap = param
               IF (ipa == param_Gammap) hmod%Gammap = param
               IF (ipa == param_cstarp) hmod%cstarp = param
               IF (ipa == param_fhot)   hmod%fhot = param
               IF (ipa == param_alphaz) hmod%alphaz = param
               IF (ipa == param_Gammaz) hmod%Gammaz = param
               IF (ipa == param_M0z)    hmod%M0z = param
               IF (ipa == param_Astarz) hmod%Astarz = param
               IF (ipa == param_Twhimz) hmod%Twhimz = param
               IF (ipa == param_eta)    hmod%eta = param
               IF (ipa == param_epsz)   hmod%epsz = param
               IF (ipa == param_beta)   hmod%beta = param
               IF (ipa == param_betap)  hmod%betap = param
               IF (ipa == param_betaz)  hmod%betaz = param
            ELSE IF (hmod%HMx_mode == 4) THEN
               IF (ipa == param_alpha) THEN
                  hmod%A_alpha = 0.
                  hmod%B_alpha = 0.
                  hmod%C_alpha = 0.
                  hmod%D_alpha = param
               ELSE IF (ipa == param_eps) THEN
                  hmod%A_eps = 0.
                  hmod%B_eps = 0.
                  hmod%C_eps = 0.
                  hmod%D_eps = log10(param)
               ELSE IF (ipa == param_Gamma) THEN
                  hmod%A_Gamma = 0.
                  hmod%B_Gamma = 0.
                  hmod%C_Gamma = 0.
                  hmod%D_Gamma = param
               ELSE IF (ipa == param_M0) THEN
                  hmod%A_M0 = 0.
                  hmod%B_M0 = 0.
                  hmod%C_M0 = 0.
                  hmod%D_M0 = log10(param)
               ELSE IF (ipa == param_Astar) THEN
                  hmod%A_Astar = 0.
                  hmod%B_Astar = 0.
                  hmod%C_Astar = 0.
                  hmod%D_Astar = param
               ELSE IF (ipa == param_Twhim) THEN
                  hmod%A_Twhim = 0.
                  hmod%B_Twhim = 0.
                  hmod%C_Twhim = 0.
                  hmod%D_Twhim = log10(param)
               ELSE
                  STOP 'HMx_DRIVER: Parameter variations not supported for this HMx_mode'
               END IF
            END IF

            CALL init_halomod(scale_factor_z(z), hmod, cosm, verbose)
            CALL print_halomod(hmod, cosm, verbose)

            ! DO NOT DELETE THIS
            ! It is only used to print values to the screen later
            ! For example, mass, which is inconvenient if written out in full
            IF (ilog) THEN
               param_neat = log10(param)
            ELSE
               param_neat = param
            END IF

            ! Write out halo matter and electron-pressure profile information
            ! All the string crap is in the loop for a reason
            DO j = 10, 16
               base = 'data/profile_mass_'
               ext = '_param_'
               base = number_file(base, j, ext)
               mid = '_value_'
               ext = '.dat'
               outfile = number_file2(base, ipa, mid, i, ext)
               mass = 10.**j
               CALL write_halo_profiles(mass, hmod, cosm, outfile)
            END DO

            ! Write out halo mass fraction information
            base = 'data/mass_fractions_param_'
            outfile = number_file2(base, ipa, mid, i, ext)
            CALL write_halo_fractions(hmod, cosm, outfile)

            CALL calculate_HMx_a(fields, nf, k, nk, pow_li, pow_2h, pow_1h, pow_hm, hmod, cosm, verbose)

            ! Loop over field combinations and do the calculation
            DO j1 = 1, nf
               DO j2 = j1, nf

                  ! Set output file
                  base = 'data/power_param_'
                  mid = '_value_'
                  ext = '_'
                  outfile = number_file2(base, ipa, mid, i, ext)

                  mid = ''
                  ext = '.dat'
                  outfile = number_file2(outfile, fields(j1), mid, fields(j2), ext)

                  ! Write progress to screen
                  WRITE (*, fmt='(4I5,F14.7,A50)') ipa, i, fields(j1), fields(j2), param_neat, TRIM(outfile)

                  ! Do the halo-model calculation and write to file
                  CALL write_power(k, pow_li, pow_2h(j1, j2, :), pow_1h(j1, j2, :), pow_hm(j1, j2, :), nk, outfile, verbose=.FALSE.)

               END DO
            END DO

         END DO

      END DO

   END SUBROUTINE baryon_parameter_variations

   SUBROUTINE spectra_fields_3D(icosmo, ihm)

      ! 3D spectra for a user choice of fields
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: icosmo
      INTEGER, INTENT(INOUT) :: ihm
      REAL, ALLOCATABLE :: k(:), a(:)
      REAL, ALLOCATABLE :: pow_li(:, :), pow_2h(:, :, :, :), pow_1h(:, :, :, :), pow_hm(:, :, :, :)
      INTEGER :: i
      INTEGER :: na, ip(2)
      CHARACTER(len=256) :: base
      TYPE(cosmology) :: cosm
      !TYPE(halomod) :: hmod

      REAL, PARAMETER :: kmin = 1e-3
      REAL, PARAMETER :: kmax = 1e2
      INTEGER, PARAMETER :: nk = 128
      REAL, PARAMETER :: zmin = 0.
      REAL, PARAMETER :: zmax = 4.
      INTEGER, PARAMETER :: nz = 16
      INTEGER, PARAMETER :: nf = 2
      LOGICAL, PARAMETER :: verbose = .TRUE.

      ! Assigns the cosmological model
      !icosmo = -1 ! 1 - Boring
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)  

      ! Set number of k points and k range (log spaced)
      ! The range kmin=1e-3 to kmax=1e4 is necessary to compare to HMcode
      CALL fill_array(log(kmin), log(kmax), k, nk)
      k = exp(k)

      !Set the number of redshifts and range (linearly spaced) and convert z -> a
      CALL fill_array(zmin, zmax, a, nz)
      a = 1./(1.+a)
      na = nz

      ! Choose the field types
      DO i = 1, nf
         WRITE (*, *) 'HMx_driver: Choose halo', i
         CALL set_halo_type(ip(i))
      END DO

      ! User chooses halo model
      !CALL assign_halomod(ihm, hmod, verbose)
      !CALL calculate_HMx(ip, nf, k, nk, a, na, pow_li, pow_2h, pow_1h, pow_hm, hmod, cosm, verbose)
      CALL calculate_HMx_full(ip, k, a, pow_li, pow_2h, pow_1h, pow_hm, nf, nk, na, cosm, ihm)

      ! Write out data
      base = 'data/power'
      CALL write_power_a_multiple(k, a, pow_li, pow_2h(1, 2, :, :), pow_1h(1, 2, :, :), &
         pow_hm(1, 2, :, :), nk, na, base, verbose)

   END SUBROUTINE spectra_fields_3D

   SUBROUTINE bias_fields_3D(imode, icosmo, ihm)

      ! Create 3D bias function
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: imode
      INTEGER, INTENT(INOUT) :: icosmo
      INTEGER, INTENT(INOUT) :: ihm
      REAL, ALLOCATABLE :: k(:), a(:)
      REAL, ALLOCATABLE :: pow_li(:, :), pow_2h(:, :, :, :), pow_1h(:, :, :, :), pow_hm(:, :, :, :)
      REAL :: zmin, zmax
      INTEGER :: j1, j2
      INTEGER :: nz, na
      INTEGER, ALLOCATABLE :: fields(:)
      CHARACTER(len=256) :: base
      TYPE(cosmology) :: cosm
      TYPE(halomod) :: hmod

      REAL, PARAMETER :: kmin = 1e-3
      REAL, PARAMETER :: kmax = 1e2
      INTEGER, PARAMETER :: nk = 128
      INTEGER, PARAMETER :: nf = 2
      LOGICAL, PARAMETER :: verbose = .TRUE.

      ! 18 - User choice
      ! 40 - Halo bias
      ! 48 - HI bias

      ! Assign the cosmological model
      IF (imode == 40) icosmo = 4  ! BAHAMAS
      IF (imode == 48) icosmo = 27 ! Illustris TNG 75
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! Assign the halo model
      !IF(imode==48) ihm=3 ! This is what I sent Kiyo & Richard
      IF (imode == 48) ihm = 25 ! Villaescua-Navarro halo model
      CALL assign_halomod(ihm, hmod, verbose)

      ! Set number of k points and k range (log spaced)
      CALL fill_array(log(kmin), log(kmax), k, nk)
      k = exp(k)

      ! Set the number of redshifts and range (linearly spaced) and convert z -> a
      IF (imode == 18 .OR. imode == 48) THEN
         IF (imode == 18) THEN
            nz = 16
            zmin = 0.
            zmax = 4.
         ELSE IF (imode == 48) THEN
            nz = 6
            zmin = 0.
            zmax = 5.
         END IF
         CALL fill_array(zmin, zmax, a, nz)
         a = 1./(1.+a)
         na = nz
      ELSE IF (imode == 40) THEN
         na = 4
         ALLOCATE (a(na))
         a(1) = 0.  ! This is actually z
         a(2) = 0.5 ! This is actually z
         a(3) = 1.  ! This is actually z
         a(4) = 2.0 ! This is actually z
         a = 1./(1.+a) ! Now convert z to a
      END IF

      ! Allocate arrays for fields
      ALLOCATE (fields(nf))
      fields(1) = field_dmonly

      ! Select field type for bias study
      IF (imode == 18) CALL set_halo_type(fields(2))
      IF (imode == 40) fields(2) = 9  ! Central galaxies/haloes
      IF (imode == 48) fields(2) = 12 ! HI

      !CALL calculate_HMx(fields, nf, k, nk, a, na, pows_li, pows_2h, pows_1h, pows_hm, hmod, cosm, verbose)
      CALL calculate_HMx_full(fields, k, a, pow_li, pow_2h, pow_1h, pow_hm, nf, nk, na, cosm, ihm)

      DO j1 = 1, nf
         DO j2 = j1, nf

            IF (j1 == 1 .AND. j2 == 1) THEN
               ! DMONLY-DMONLY
               base = 'data/power_mm'
            ELSE IF (j1 == 1 .AND. j2 == 2) THEN
               ! DMONLY-field
               base = 'data/power_mf'
            ELSE IF (j1 == 2 .AND. j2 == 2) THEN
               ! field-field
               base = 'data/power_ff'
            END IF

            CALL write_power_a_multiple(k, a, pow_li, pow_2h(j1, j2, :, :), pow_1h(j1, j2, :, :), pow_hm(j1, j2, :, :), nk, na, base, verbose)

         END DO
      END DO

   END SUBROUTINE bias_fields_3D

   SUBROUTINE create_CCL_benchmark(icosmo, ihm)

      ! Create CCL benchmark data
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: icosmo
      INTEGER, INTENT(INOUT) :: ihm
      REAL, ALLOCATABLE :: k(:)
      REAL :: z
      REAL :: kmin, kmax
      INTEGER :: i, j
      INTEGER :: nk
      REAL, ALLOCATABLE :: pow_li(:), pow_2h(:, :, :), pow_1h(:, :, :), pow_hm(:, :, :)
      CHARACTER(len=256) :: outfile
      TYPE(cosmology) :: cosm
      TYPE(halomod) :: hmod

      INTEGER, PARAMETER :: field(1) = field_dmonly
      LOGICAL, PARAMETER :: Alonso_k = .TRUE.
      LOGICAL, PARAMETER :: verbose = .TRUE.

      ! Set number of k points and k range (log spaced)
      IF (Alonso_k) THEN
         ! This is *almost* kmin=1e-4, kmax=1e2, but minus the end points
         nk = 256
         kmin = 1.027350768179302566e-04
         kmax = 9.733773809039202263e+01
      ELSE
         nk = 128
         kmin = 1e-3
         kmax = 1e1
      END IF
      CALL fill_array(log(kmin), log(kmax), k, nk)
      k = exp(k)
      ALLOCATE (pow_li(nk), pow_2h(1, 1, nk), pow_1h(1, 1, nk), pow_hm(1, 1, nk))

      ! Set the halo model to that for CCL tests
      ihm = 9

      ! Loop over tests
      DO j = 1, 3

         ! Assigns the cosmological model
         IF (j == 1) icosmo = 1
         IF (j == 2) icosmo = 2
         IF (j == 3) icosmo = 3
         CALL assign_cosmology(icosmo, cosm, verbose)
         CALL init_cosmology(cosm)
         CALL print_cosmology(cosm)

         ! Loop over redshifts
         DO i = 1, 2

            ! Sets the redshift and file names
            IF (i == 1) THEN
               z = 0.0
               IF (j == 1) outfile = '/Users/Mead/Physics/HMx/CCL/HMx_power_model1_z0.txt'
               IF (j == 2) outfile = '/Users/Mead/Physics/HMx/CCL/HMx_power_model2_z0.txt'
               IF (j == 3) outfile = '/Users/Mead/Physics/HMx/CCL/HMx_power_model3_z0.txt'
            ELSE IF (i == 2) THEN
               z = 1.0
               IF (j == 1) outfile = '/Users/Mead/Physics/HMx/CCL/HMx_power_model1_z1.txt'
               IF (j == 2) outfile = '/Users/Mead/Physics/HMx/CCL/HMx_power_model2_z1.txt'
               IF (j == 3) outfile = '/Users/Mead/Physics/HMx/CCL/HMx_power_model3_z1.txt'
            END IF

            ! Initialise the halo-model calculation
            CALL assign_halomod(ihm, hmod, verbose)
            CALL init_halomod(scale_factor_z(z), hmod, cosm, verbose)
            CALL print_halomod(hmod, cosm, verbose)

            ! Do the halo-model calculation
            CALL calculate_HMx_a(field, 1, k, nk, pow_li, pow_2h, pow_1h, pow_hm, hmod, cosm, verbose)

            ! Write out the results
            CALL write_power(k, pow_li, pow_2h(1, 1, :), pow_1h(1, 1, :), pow_hm(1, 1, :), nk, outfile, verbose)

         END DO

      END DO

   END SUBROUTINE

   SUBROUTINE Ma2015_Fig1(icosmo, ihm)

      ! Ma et al. Fig. 1
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: icosmo
      INTEGER, INTENT(INOUT) :: ihm
      TYPE(cosmology) :: cosm
      TYPE(halomod) :: hmod
      INTEGER :: i
      REAL :: r, rs, rv, c, Mh, rh, r500c, m500c

      REAL, PARAMETER :: z = 0.0
      REAL, PARAMETER :: M = 1e15    ! Halo virial? mass [Msun]
      REAL, PARAMETER :: rmin = 1e-3 ! Minimum radius [Mpc]
      REAL, PARAMETER :: rmax = 8    ! Maximum radius [Mpc]
      INTEGER, PARAMETER :: nr = 512 ! Number of points in radius
      LOGICAL, PARAMETER :: verbose = .TRUE.
      LOGICAL, PARAMETER :: real_space = .TRUE.
      INTEGER, PARAMETER :: itype = field_electron_pressure ! electron pressure

      ! Set the cosmology
      icosmo = 3
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! Set the halo model
      ihm = 4
      CALL assign_halomod(ihm, hmod, verbose)
      CALL init_halomod(scale_factor_z(z), hmod, cosm, verbose)
      CALL print_halomod(hmod, cosm, verbose)

      IF (hmod%has_mass_conversions .EQV. .FALSE.) CALL convert_mass_definitions(hmod, cosm)

      Mh = M*cosm%h ! This virial mass is now [Msun/h]

      rv = exp(find(log(Mh), hmod%log_m, log(hmod%rv), hmod%n, 3, 3, 2)) ! [Mpc/h]
      c = find(log(Mh), hmod%log_m, hmod%c, hmod%n, 3, 3, 2)
      rs = rv/c ! [Mpc/h]

      m500c = exp(find(log(Mh), hmod%log_m, log(hmod%m500c), hmod%n, 3, 3, 2)) ! [Mpc/h]
      r500c = exp(find(log(Mh), hmod%log_m, log(hmod%r500c), hmod%n, 3, 3, 2)) ! [Mpc/h]

      WRITE (*, *) 'MA2015_FIG1: Making data for this figure'
      WRITE (*, *) 'MA2015_FIG1: Redshift:', hmod%z
      WRITE (*, *) 'MA2015_FIG1: Virial radius [Mpc]:', rv/cosm%h
      WRITE (*, *) 'MA2015_FIG1: Virial radius [Mpc/h]:', rv
      WRITE (*, *) 'MA2015_FIG1: r_500,c [Mpc]:', r500c/cosm%h
      WRITE (*, *) 'MA2015_FIG1: r_500,c [Mpc/h]:', r500c
      WRITE (*, *) 'MA2015_FIG1: r_500,c / r_v:', r500c/rv
      WRITE (*, *) 'MA2015_FIG1: Virial halo mass [log10 Msun]:', log10(M)
      WRITE (*, *) 'MA2015_FIG1: Virial halo mass [log10 Msun/h]:', log10(Mh)
      WRITE (*, *) 'MA2015_FIG1: M_500,c [log10 Msun]:', log10(M500c/cosm%h)
      WRITE (*, *) 'MA2015_FIG1: M_500,c [log10 Msun/h]:', log10(M500c)
      WRITE (*, *) 'MA2015_FIG1: M_500,c / M_v:', M500c/Mh
      WRITE (*, *) 'MA2015_FIG1: Halo concentraiton:', c

      OPEN (7, file='data/YinZhe_Fig1.dat')
      DO i = 1, nr
         r = progression(rmin, rmax, i, nr) ! Radius [Mpc]
         rh = r*cosm%h ! Convert [Mpc/h]
         WRITE (7, *) r, UPP(real_space, rh, Mh, rv, rs, hmod, cosm)*r**2, win(real_space, itype, rh, Mh, rv, rs, hmod, cosm)*r**2
      END DO
      CLOSE (7)

      WRITE (*, *) 'MA2015_FIG1: Done'
      WRITE (*, *)

   END SUBROUTINE Ma2015_Fig1

   SUBROUTINE Mead2017(icosmo, ihm)

      ! Mead (2017) dark energy results
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: icosmo
      INTEGER, INTENT(INOUT) :: ihm
      REAL, ALLOCATABLE :: k(:)
      REAL :: z
      REAL, ALLOCATABLE :: pow_li(:), pow_2h(:, :, :), pow_1h(:, :, :), pow_hm(:, :, :)
      INTEGER :: i
      CHARACTER(len=256) :: outfile
      TYPE(cosmology) :: cosm
      TYPE(halomod) :: hmod

      REAL, PARAMETER :: kmin = 1e-3
      REAL, PARAMETER :: kmax = 1e1
      INTEGER, PARAMETER :: nk = 128
      INTEGER, PARAMETER :: ncos = 25
      INTEGER, PARAMETER :: field(1) = field_dmonly
      CHARACTER(len=256), PARAMETER :: dir = 'data' ! Directory for output
      LOGICAL, PARAMETER :: verbose = .TRUE.

      ! Set number of k points and k range (log spaced)

      CALL fill_array(log(kmin), log(kmax), k, nk)
      k = exp(k)
      ALLOCATE (pow_li(nk), pow_2h(1, 1, nk), pow_1h(1, 1, nk), pow_hm(1, 1, nk))

      ! Loop over cosmologies
      DO i = 1, ncos

         IF (i == 1) THEN
            ! LCDM; z=0
            icosmo = 1
            outfile = 'LCDM_z0'
            z = 0.0
         ELSE IF (i == 2) THEN
            ! OCDM; z=0
            icosmo = 5
            outfile = 'OCDM_z0'
            z = 0.0
         ELSE IF (i == 3) THEN
            ! EdS; z=0
            icosmo = 15
            outfile = 'SCDM_z0'
            z = 0.0
         ELSE IF (i == 4) THEN
            ! w = -0.7; z=0
            icosmo = 16
            outfile = 'w-0.7_z0'
            z = 0.0
         ELSE IF (i == 5) THEN
            ! w = -1.3; z=0
            icosmo = 17
            outfile = 'w-1.3_z0'
            z = 0.0
         ELSE IF (i == 6) THEN
            ! wa = 0.5; z=0
            icosmo = 18
            outfile = 'wa0.5_z0'
            z = 0.0
         ELSE IF (i == 7) THEN
            ! wa = -0.5; z=0
            icosmo = 19
            outfile = 'wa-0.5_z0'
            z = 0.0
         ELSE IF (i == 8) THEN
            ! w = -0.7; wa = -1.5; z=0
            icosmo = 20
            outfile = 'w0-0.7wa-1.5_z0'
            z = 0.0
         ELSE IF (i == 9) THEN
            ! w = -1.3; wa = 0.5; z=0
            icosmo = 21
            outfile = 'w0-1.3wa0.5_z0'
            z = 0.0
         ELSE IF (i == 10) THEN
            ! OCDM; z=1
            icosmo = 5
            outfile = 'OCDM_z1'
            z = 1.0
         ELSE IF (i == 11) THEN
            ! w = -0.7; z=1
            icosmo = 16
            outfile = 'w-0.7_z1'
            z = 1.0
         ELSE IF (i == 12) THEN
            ! w = -1.3; z=1
            icosmo = 17
            outfile = 'w-1.3_z1'
            z = 1.0
         ELSE IF (i == 13) THEN
            ! wa = 0.5; z=1
            icosmo = 18
            outfile = 'wa0.5_z1'
            z = 1.0
         ELSE IF (i == 14) THEN
            ! wa = -0.5; z=1
            icosmo = 19
            outfile = 'wa-0.5_z1'
            z = 1.0
         ELSE IF (i == 15) THEN
            ! w = -0.7; wa = -1.5; z=1
            icosmo = 20
            outfile = 'w0-0.7wa-1.5_z1'
            z = 1.0
         ELSE IF (i == 16) THEN
            ! w = -1.3; wa = 0.5; z=1
            icosmo = 21
            outfile = 'w0-1.3wa0.5_z1'
            z = 1.0
         ELSE IF (i == 17) THEN
            ! OCDM; z=1; LCDM
            icosmo = 29
            outfile = 'LCDM_OCDM_z1'
            z = 1.0
         ELSE IF (i == 18) THEN
            ! w = -0.7; z=1; LCDM
            icosmo = 30
            outfile = 'LCDM_w-0.7_z1'
            z = 1.0
         ELSE IF (i == 19) THEN
            ! w = -1.3; z=1; LCDM
            icosmo = 31
            outfile = 'LCDM_w-1.3_z1'
            z = 1.0
         ELSE IF (i == 20) THEN
            ! wa = 0.5; z=1; LCDM
            icosmo = 32
            outfile = 'LCDM_wa0.5_z1'
            z = 1.0
         ELSE IF (i == 21) THEN
            ! wa = -0.5; z=1; LCDM
            icosmo = 33
            outfile = 'LCDM_wa-0.5_z1'
            z = 1.0
         ELSE IF (i == 22) THEN
            ! w = -0.7; wa = -1.5; z=1; LCDM
            icosmo = 34
            outfile = 'LCDM_w0-0.7wa-1.5_z1'
            z = 1.0
         ELSE IF (i == 23) THEN
            ! w = -1.3; wa = 0.5; z=1; LCDM
            icosmo = 35
            outfile = 'LCDM_w0-1.3wa0.5_z1'
            z = 1.0
         ELSE IF (i == 24) THEN
            ! EdS; z=1
            icosmo = 15
            outfile = 'SCDM_z1'
            z = 1.0
         ELSE IF (i == 25) THEN
            ! EdS; z=1; LCDM
            icosmo = 36
            outfile = 'LCDM_SCDM_z1'
            z = 1.0
         END IF

         ! Assigns the cosmological model
         CALL assign_cosmology(icosmo, cosm, verbose)
         CALL init_cosmology(cosm)
         CALL print_cosmology(cosm)

         ! Initiliasation for the halomodel calcualtion
         ihm = 12
         CALL assign_halomod(ihm, hmod, verbose)
         CALL init_halomod(scale_factor_z(z), hmod, cosm, verbose)
         CALL print_halomod(hmod, cosm, verbose)

         ! Do the halo-model calculation
         !field = field_dmonly
         CALL calculate_HMx_a(field, 1, k, nk, pow_li, pow_2h, pow_1h, pow_hm, hmod, cosm, verbose)

         ! Write out the results
         outfile = TRIM(dir)//'/power_'//TRIM(outfile)//'.dat'
         CALL write_power(k, pow_li, pow_2h(1, 1, :), pow_1h(1, 1, :), pow_hm(1, 1, :), nk, outfile, verbose)

      END DO

   END SUBROUTINE Mead2017

   SUBROUTINE projection_tests(icosmo, ihm)

      ! Projection tests
      ! TODO: Fix why these tests are super slow if ihm=3 and Cl_yy looks really fucked up
      ! TODO: I think this is probably because of the pressure evolution in ihm=3
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: icosmo
      INTEGER, INTENT(INOUT) :: ihm
      REAL :: error
      CHARACTER(len=1) :: crap
      CHARACTER(len=256), ALLOCATABLE :: ixx_names(:)
      CHARACTER(len=256) :: infile, outfile
      INTEGER :: i, ii, jj
      INTEGER, ALLOCATABLE :: ixx(:)
      REAL, ALLOCATABLE :: ell(:), Cl(:, :, :), Cl_bm(:)
      TYPE(cosmology) :: cosm
      !TYPE(halomod) :: hmod

      REAL, PARAMETER :: lmin = 1e0
      REAL, PARAMETER :: lmax = 1e4
      INTEGER, PARAMETER :: nl = 128
      LOGICAL, PARAMETER :: verbose = .TRUE.
      INTEGER, PARAMETER :: nx = 3 ! Number of tests
      REAL, PARAMETER :: tolerance = 3e-3
      LOGICAL :: ifail = .FALSE. ! Initially assume tests pass

      ALLOCATE (ixx(nx))
      ixx(1) = tracer_RCSLenS
      ixx(2) = tracer_Compton_y
      ixx(3) = tracer_CMB_lensing

      ALLOCATE (ixx_names(nx))
      ixx_names(1) = 'RCSLenS'
      ixx_names(2) = 'y'
      ixx_names(3) = 'CMB'

      ! Set the ell range (should probably read this in from a file)
      CALL fill_array(log(lmin), log(lmax), ell, nl)
      ell = exp(ell)
      ALLOCATE (Cl_bm(nl), Cl(nl, nx, nx))

      ! Assigns the cosmological model
      icosmo = 1
      CALL assign_cosmology(icosmo, cosm, verbose=.FALSE.)
      CALL init_cosmology(cosm)

      ! Assign the halo model
      ihm = 18
      !CALL assign_halomod(ihm, hmod, verbose=.FALSE.)
      !CALL xpow_halomod(ixx, nx, ell, Cl, nl, hmod, cosm, verbose=.TRUE.)
      CALL xpow_halomod(ixx, nx, ell, Cl, nl, cosm, ihm, verbose=.TRUE.)

      ! Write data
      DO ii = 1, nx
         DO jj = 1, nx
            outfile = 'data/cl_'//TRIM(ixx_names(ii))//'_'//TRIM(ixx_names(jj))//'.dat'
            CALL write_Cl(ell, Cl(:, ii, jj), nl, outfile, verbose)
         END DO
      END DO

      ! Loop over and check values
      DO ii = 1, nx
         DO jj = 1, nx

            infile = 'benchmarks/cl_'//TRIM(ixx_names(ii))//'_'//TRIM(ixx_names(jj))//'.txt'
            IF (file_length(infile, verbose=.FALSE.) .NE. nl) STOP 'HMx_DRIVER: Error, benchmark file is not the correct lenght'

            OPEN (7, file=infile)
            DO i = 1, nl
               READ (7, *) crap, Cl_bm(i)
            END DO
            CLOSE (7)

            DO i = 1, nl
               error = ABS(-1.+Cl(i, ii, jj)/Cl_bm(i))
               IF (error > tolerance) THEN
                  WRITE (*, *) 'HMx_DRIVER: Tracer 1: ', TRIM(ixx_names(ii))
                  WRITE (*, *) 'HMx_DRIVER: Tracer 2: ', TRIM(ixx_names(jj))
                  WRITE (*, *) 'HMx_DRIVER: ell:', ell(i)
                  WRITE (*, *) 'HMx_DRIVER: Benchmark C(l):', Cl_bm(i)
                  WRITE (*, *) 'HMx_DRIVER: Calculated C(l):', Cl(i, ii, jj)
                  WRITE (*, *) 'HMx_DRIVER: Error:', error
                  WRITE (*, *)
                  ifail = .TRUE.
               END IF
            END DO

         END DO
      END DO

      WRITE (*, *) 'HMx_DRIVER: Limber tests should take ~18s to run'
      IF (ifail) THEN
         WRITE (*, *) 'HMx_DRIVER: Limber tests failed'
      ELSE
         WRITE (*, *) 'HMx_DRIVER: Limber tests passed'
      END IF
      WRITE (*, *)

   END SUBROUTINE projection_tests

   SUBROUTINE halo_void_model(icosmo, ihm)

      ! Halo-void model
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: icosmo
      INTEGER, INTENT(INOUT) :: ihm
      INTEGER :: i
      REAL, ALLOCATABLE :: k(:)
      REAL, ALLOCATABLE :: powd_li(:), powd_2h(:), powd_1h(:), powd_hm(:)
      TYPE(cosmology) :: cosm
      TYPE(halomod) :: hmod

      REAL, PARAMETER :: kmin = 1e-3
      REAL, PARAMETER :: kmax = 1e2
      INTEGER, PARAMETER :: nk = 128
      REAL, PARAMETER :: z = 0.0
      INTEGER, PARAMETER :: field(1) = field_dmonly
      LOGICAL, PARAMETER :: verbose = .TRUE.

      ! Set number of k points and k range (log spaced)
      CALL fill_array(log(kmin), log(kmax), k, nk)
      k = exp(k)

      ! Assigns the cosmological model
      icosmo = 1 ! 1 - Boring
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! Initiliasation for the halomodel calcualtion
      ihm = 16 ! 16 - ?
      CALL assign_halomod(ihm, hmod, verbose)
      CALL init_halomod(scale_factor_z(z), hmod, cosm, verbose)
      CALL print_halomod(hmod, cosm, verbose)

      ! Allocate arrays
      ALLOCATE (powd_li(nk), powd_2h(nk), powd_1h(nk), powd_hm(nk))

      ! Do the halo-model calculation
      CALL calculate_HMx_a(field, 1, k, nk, powd_li, powd_2h, powd_1h, powd_hm, hmod, cosm, verbose)

      ! Write the one-void term if necessary
      OPEN (7, file='data/power_halovoid.dat')
      DO i = 1, nk
         WRITE (7, *) k(i), powd_li(i), powd_2h(i), powd_1h(i), powd_hm(i), p_1void(k(i), hmod)
      END DO
      CLOSE (7)

   END SUBROUTINE halo_void_model

   SUBROUTINE halo_model_tests(icosmo, ihm)

      ! Automated testing
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: icosmo
      INTEGER, INTENT(INOUT) :: ihm
      REAL, ALLOCATABLE :: k(:), a(:), pow_ka(:, :)
      REAL, ALLOCATABLE :: pow_li(:, :), pow_2h(:, :), pow_1h(:, :), pow_hm(:, :)
      REAL :: error, error_max
      INTEGER :: i, j, itest
      INTEGER :: nk, na
      CHARACTER(len=1) :: crap
      CHARACTER(len=256) :: base, infile
      TYPE(cosmology) :: cosm
      !TYPE(halomod) :: hmod

      REAL, PARAMETER :: kmin_test = 1e-3
      REAL, PARAMETER :: kmax_test = 1e2
      REAL, PARAMETER :: tolerance = 3e-3
      INTEGER, PARAMETER :: field(1) = field_dmonly
      LOGICAL :: ifail = .FALSE.
      LOGICAL, PARAMETER :: verbose = .FALSE.

      ! Loop over different tests
      DO itest = 1, 3

         IF (itest == 1) THEN
            infile = 'benchmarks/power_HMcode_Mead.txt'
            ihm = 1
            base = 'data/power_HMx_Mead'
            icosmo = 1
         ELSE IF (itest == 2) THEN
            infile = 'benchmarks/power_HMcode_basic.txt'
            ihm = 2
            base = 'data/power_HMx_basic'
            icosmo = 1
         ELSE IF (itest == 3) THEN
            infile = 'benchmarks/power_HMcode_standard.txt'
            ihm = 5
            base = 'data/power_HMx_standard'
            icosmo = 1
         END IF

         ! Allocate arrays of k and a
         nk = 128
         na = 16
         ALLOCATE (k(nk), a(na), pow_ka(nk, na))

         ! Read-in test data
         OPEN (7, file=infile, status='old')
         DO i = 0, nk
            IF (i == 0) THEN
               READ (7, *) crap, (a(j), j=1, na)
            ELSE
               READ (7, *) k(i), (pow_ka(i, j), j=1, na)
            END IF
         END DO
         CLOSE (7)

         ! Convert read-in z to a; note that 'a' is the correct argument
         ! for the function here because actually read-in a is z.
         DO i = 1, na
            a(i) = scale_factor_z(a(i))
         END DO

         ! Assigns the cosmological model
         CALL assign_cosmology(icosmo, cosm, verbose)
         CALL init_cosmology(cosm)
         CALL print_cosmology(cosm)

         ! Assign the halo model
         !CALL assign_halomod(ihm, hmod, verbose)
         !field = field_dmonly
         !CALL calculate_HMx(field, 1, k, nk, a, na, pows_li, pows_2h, pows_1h, pows_hm, hmod, cosm, verbose)
         CALL calculate_halomod_full(k, a, pow_li, pow_2h, pow_1h, pow_hm, nk, na, cosm, ihm)

         ! Loop over k and a
         error_max = 0.
         DO j = 1, na
            DO i = 1, nk
               error = ABS(-1.+pow_hm(i, j)/pow_ka(i, j))
               IF (kmin_test <= k(i) .AND. k(i) <= kmax_test) THEN
                  IF (error > error_max) error_max = error
                  IF (error > tolerance) THEN
                     WRITE (*, *) 'HMx_DRIVER: Test:', itest
                     WRITE (*, *) 'HMx_DRIVER: Wavenumber [h/Mpc]:', k(i)
                     WRITE (*, *) 'HMx_DRIVER: Scale-factor:', a(j)
                     WRITE (*, *) 'HMx_DRIVER: Expected power:', pow_ka(i, j)
                     WRITE (*, *) 'HMx_DRIVER: Model power:', pow_hm(i, j)
                     WRITE (*, *) 'HMx_DRIVER: Tolerance:', tolerance
                     WRITE (*, *) 'HMx_DRIVER: Error:', error
                     WRITE (*, *)
                     ifail = .TRUE.
                  END IF
               END IF
            END DO
         END DO

         WRITE (*, *) 'HMx_DRIVER: Test:', itest
         WRITE (*, *) 'HMx_DRIVER: Tolerance:', tolerance
         WRITE (*, *) 'HMx_DRIVER: Max error:', error_max
         WRITE (*, *)

         ! Write data to file
         !base='data/power'
         !CALL write_power_a_multiple(k,a,powa_lin,powa_2h,powa_1h,powa_full,nk,na,base,verbose_tests)
         CALL write_power_a_multiple(k, a, pow_li, pow_2h, pow_1h, pow_hm, nk, na, base, verbose)

         DEALLOCATE (k, a, pow_ka)

      END DO

      IF (ifail) THEN
         STOP 'HMx_DRIVER: Error, tests failed'
      ELSE
         WRITE (*, *) 'HMx_DRIVER: Tests should take around 0.50 seconds to run'
         WRITE (*, *) 'HMx_DRIVER: Tests passed'
         WRITE (*, *)
      END IF

   END SUBROUTINE halo_model_tests

   SUBROUTINE emulator_comparison(imode, icosmo, ihm)

      ! Comparison with FrankenEmu or Mira Titan
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: imode
      INTEGER, INTENT(INOUT) :: icosmo
      INTEGER, INTENT(INOUT) :: ihm
      REAL, ALLOCATABLE :: z(:), a(:)
      REAL, ALLOCATABLE :: pow_ql(:), pow_oh(:), pow_hf(:)
      REAL, ALLOCATABLE :: k_sim(:), pow_sim(:)
      REAL, ALLOCATABLE :: pow_li(:), pow_2h(:, :, :), pow_1h(:, :, :), pow_hm(:, :, :)
      REAL :: HALOFIT_knl, HALOFIT_neff, HALOFIT_ncur
      INTEGER :: icos, iz, ik
      INTEGER :: field(1), nz, na, nk, ncos
      CHARACTER(len=256) :: outfile
      TYPE(cosmology) :: cosm
      TYPE(halomod) :: hmod

      REAL, PARAMETER :: m_nu_min = 0.00 ! Minimum neutrino mass [eV]
      REAL, PARAMETER :: m_nu_max = 0.94 ! Maximum neutrino mass [eV]
      LOGICAL, PARAMETER :: verbose = .TRUE.
      CHARACTER(len=256), PARAMETER :: base = 'data/cosmo'
      CHARACTER(len=256), PARAMETER :: mid = '_z'
      CHARACTER(len=256), PARAMETER :: ext = '.dat'
      !CHARACTER(len=256), PARAMETER :: cosmology_file = 'data/emulator_cosmologies.dat'
      INTEGER, PARAMETER :: ihf = halofit_CAMB

      ! Number of cosmological models
      IF (imode == 28 .OR. imode == 30) ncos = 37 ! Franken Emu
      IF (imode == 27 .OR. imode == 29) ncos = 36 ! Mira Titan
      IF (imode == 70 .OR. imode == 71) ncos = 37 ! Cosmic Emu
      IF (imode == 76) ncos = 10 ! Mira Titan with massless neutrinos
      IF (imode == 89) ncos = 10
      IF (imode == 100) ncos = 0

      ! Set number of redshifts
      IF (imode == 28 .OR. imode == 30 .OR. imode == 100) THEN
         nz = 6 ! Franken Emu
      ELSE IF (imode == 27 .OR. imode == 29 .OR. imode == 76 .OR. imode == 89) THEN
         nz = 4 ! Mira Titan
      ELSE IF (imode == 70 .OR. imode == 71) THEN
         nz = 3 ! Cosmic Emu
      ELSE
         STOP 'EMULATOR_COMPARISON:  Error, imode not specified correctly for EMU'
      END IF

      ! Set the random number generator for the random cosmologies
      IF (imode == 30 .OR. imode == 29 .OR. imode == 71) THEN
         CALL RNG_set(seed=0)
      END IF

      ! Allocate arrays
      na = nz
      ALLOCATE (z(nz), a(nz))

      ! Set redshifts
      IF (imode == 28 .OR. imode == 30 .OR. imode == 100) THEN
         ! Franken Emu (z up to 4)
         z(1) = 0.0
         z(2) = 0.5
         z(3) = 1.0
         z(4) = 2.0
         z(5) = 3.0
         z(6) = 4.0
      ELSE IF (imode == 27 .OR. imode == 29 .OR. imode == 76 .OR. imode == 89) THEN
         ! Mira Titan (z up to 2)
         z(1) = 0.0
         z(2) = 0.5
         z(3) = 1.0
         z(4) = 2.0
      ELSE IF (imode == 70 .OR. imode == 71) THEN
         ! Cosmic Emu (z up to 1)
         z(1) = 0.0
         z(2) = 0.5
         z(3) = 1.0
      ELSE
         STOP 'EMULATOR_COMPARISON:  Error, imode not specified correctly for EMU'
      END IF

      ! Set scale factors
      a = scale_factor_z(z)

      ! Initiliasation for the halomodel calcualtion
      CALL assign_halomod(ihm, hmod, verbose)

      ! Loop over cosmologies
      !OPEN(8, file=cosmology_file)
      DO icos = 0, ncos

         IF (imode == 27 .OR. imode == 76) THEN
            icosmo = 100+icos ! Mira Titan nodes
         ELSE IF (imode == 28 .OR. imode == 100) THEN
            icosmo = 200+icos ! Franken Emu nodes (same as cosmic emu nodes)
         ELSE IF (imode == 70) THEN
            icosmo = 300+icos ! Cosmic Emu nodes
         ELSE IF (imode == 29) THEN
            icosmo = 24 ! Random Mira Titan
         ELSE IF (imode == 30) THEN
            icosmo = 25 ! Random Franken Emu
         ELSE IF (imode == 71) THEN
            icosmo = 38 ! Random Cosmic Emu
         ELSE IF (imode == 89) THEN
            icosmo = 100 ! Mira Titan central node
         ELSE
            WRITE (*, *) 'EMULATOR_COMPARISON:', imode
            STOP 'EMULATOR_COMPARISON: Error, imode not specified correctly for EMU'
         END IF
         CALL assign_cosmology(icosmo, cosm, verbose)
         IF (imode == 89) THEN
            cosm%m_nu = progression(m_nu_min, m_nu_max, icos+1, ncos+1)
         END IF
         CALL init_cosmology(cosm)
         CALL print_cosmology(cosm)

         !WRITE(8, *) cosm%Om_m, cosm%Om_b, cosm%Om_nu, cosm%ns, cosm%h, cosm%sig8, cosm%w, cosm%wa

         ! Loop over redshift
         DO iz = 1, nz

            IF (imode == 27 .OR. imode == 29 .OR. imode == 76 .OR. imode == 89) THEN
               CALL get_MiraTitan_power_z(k_sim, pow_sim, nk, z(iz), cosm, rebin=.FALSE.)
            ELSE IF (imode == 28 .OR. imode == 30 .OR. imode == 100) THEN
               CALL get_FrankenEmu_power_z(k_sim, pow_sim, nk, z(iz), cosm, rebin=.FALSE.)
            ELSE IF (imode == 70 .OR. imode == 71) THEN
               CALL get_CosmicEmu_power_z(k_sim, pow_sim, nk, z(iz), cosm, rebin=.FALSE.)
            ELSE
               STOP 'EMULATOR_COMPARISON: Error, imode not specified correctly for EMU'
            END IF

            ALLOCATE (pow_li(nk), pow_2h(1, 1, nk), pow_1h(1, 1, nk), pow_hm(1, 1, nk))

            CALL init_halomod(a(iz), hmod, cosm, verbose=.TRUE.)
            CALL print_halomod(hmod, cosm, verbose=.TRUE.)
            field = field_dmonly
            CALL calculate_HMx_a(field, 1, k_sim, nk, pow_li, pow_2h, pow_1h, pow_hm, hmod, cosm, verbose=.FALSE.)

            ALLOCATE (pow_ql(nk), pow_oh(nk), pow_hf(nk))
            CALL calculate_halofit_a(k_sim, a(iz), pow_li, pow_ql, pow_oh, pow_hf, nk, cosm, verbose=.FALSE., ihf=ihf)

            CALL HALOFIT_init(HALOFIT_knl, HALOFIT_neff, HALOFIT_ncur, a(iz), cosm, verbose=.FALSE.)

            ! Write data
            outfile = number_file2(base, icos, mid, iz, ext)
            OPEN (7, file=outfile)
            DO ik = 1, nk
               WRITE (7, *) k_sim(ik), pow_li(ik), pow_sim(ik), &
                  pow_ql(ik), pow_oh(ik), pow_hf(ik), &
                  pow_2h(1, 1, ik), pow_1h(1, 1, ik), pow_hm(1, 1, ik), &
                  cosm%Om_m, cosm%Om_b, cosm%Om_nu, cosm%ns, cosm%h, cosm%sig8, cosm%w, cosm%wa, &
                  HALOFIT_knl, HALOFIT_neff, HALOFIT_ncur
            END DO
            CLOSE (7)

            DEALLOCATE (pow_li, pow_2h, pow_1h, pow_hm)
            DEALLOCATE (pow_ql, pow_oh, pow_hf)
            DEALLOCATE (k_sim, pow_sim)

         END DO

      END DO
      !CLOSE(8)

   END SUBROUTINE emulator_comparison

   SUBROUTINE emulator_mean_variance()

      USE statistics
      IMPLICIT NONE
      REAL, ALLOCATABLE :: k(:), a(:), Pk(:, :)
      REAL, ALLOCATABLE :: Pk_hm(:, :, :, :), Pk_emu(:, :, :)
      REAL, ALLOCATABLE :: residual(:)
      INTEGER, ALLOCATABLE :: ihms(:)
      INTEGER :: icos, ik, ia, icosmo, ihm
      INTEGER :: nk, na, ncos, nhm
      CHARACTER(len=256) :: outfile
      CHARACTER(len=32), ALLOCATABLE :: zlab(:)
      TYPE(cosmology) :: cosm

      INTEGER, PARAMETER :: emulator_version = emulator_FrankenEmu
      LOGICAL, PARAMETER :: rebin = .TRUE.
      LOGICAL, PARAMETER :: verbose = .FALSE.
      CHARACTER(len=256), PARAMETER :: outbase = 'data/emulator_mean_variance'

      ! Set redshifts
      na = 4
      ALLOCATE(a(na), zlab(na))
      a(1) = scale_factor_z(0.0)
      a(2) = scale_factor_z(0.5)
      a(3) = scale_factor_z(1.0)
      a(4) = scale_factor_z(2.0)
      zlab(1) = '0p0'
      zlab(2) = '0p5'
      zlab(3) = '1p0'
      zlab(4) = '2p0'

      ! Set cosmologies
      ncos = 37

      ! Set number of halo models
      nhm = 23
      ALLOCATE(ihms(nhm))
      ihms = [1, 3, 7, 15, 23, 27, 42, 44, 52, 68, 69, 70, 71, 72, 73, 74, 75, 76, 80, 87, 88, 90, 91]
      
      ! Loop over cosmologies
      DO icos = 1, ncos

         ! Set the cosmology
         icosmo = 200+icos
         CALL assign_cosmology(icosmo, cosm, verbose)
         CALL init_cosmology(cosm)

         ! Get emulator power
         CALL get_emulator_power(k, a, Pk, nk, na, cosm, rebin, emulator_version)
         IF(.NOT. ALLOCATED(Pk_emu)) ALLOCATE(Pk_emu(ncos, nk, na))
         Pk_emu(icos, :, :) = Pk

         ! Loop over halo models
         DO ihm = 1, nhm

            ! Calculate halo model power for each model
            CALL calculate_halomod(k, a, Pk, nk, na, cosm, ihms(ihm))
            IF(.NOT. ALLOCATED(Pk_hm)) THEN
               ALLOCATE(Pk_hm(ncos, nhm, nk, na))
               Pk_hm = 0.
            END IF
            Pk_hm(icos, ihm, :, :) = Pk

         END DO

      END DO

      ALLOCATE(residual(ncos))
      DO ihm = 1, nhm
         DO ia = 1, na
            outfile = trim(outbase)//'_hm'//trim(integer_to_string(ihms(ihm)))//'_z'//trim(zlab(ia))//'.dat'
            OPEN(7, file=outfile)
            DO ik = 1, nk
               residual = Pk_hm(:, ihm, ik, ia)/Pk_emu(:, ik, ia)
               WRITE(7, *) k(ik), mean(residual, ncos), standard_deviation(residual, ncos)
            END DO
            CLOSE(7)
         END DO
      END DO
   
   END SUBROUTINE

   SUBROUTINE big_emulator_comparison(imode)

      ! Comparison with FrankenEmu or Mira Titan
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: imode
      REAL, ALLOCATABLE :: z(:), a(:), k(:)
      REAL, ALLOCATABLE :: Pk_emu(:, :), Pk(:, :)
      INTEGER :: icosmo, icos, icomp, ik, iz
      INTEGER :: nz, na, nk, ncos
      INTEGER :: emulator_version
      CHARACTER(len=256) :: outfile, zlab, base, emulator, comp
      TYPE(cosmology) :: cosm

      REAL, PARAMETER :: m_nu_min = 0.00 ! Minimum neutrino mass [eV]
      REAL, PARAMETER :: m_nu_max = 0.94 ! Maximum neutrino mass [eV]
      LOGICAL, PARAMETER :: verbose = .TRUE.
      ! CHARACTER(len=256), PARAMETER :: base = 'data/cosmo'
      ! CHARACTER(len=256), PARAMETER :: mid = '_z'
      ! CHARACTER(len=256), PARAMETER :: ext = '.dat'
      LOGICAL, PARAMETER :: rebin = .FALSE.
      INTEGER, PARAMETER :: ihf = halofit_CAMB
      INTEGER, PARAMETER :: ncomp = 10

      ! Number of cosmological models
      IF (imode == 96) THEN
         ncos = 37 ! Cosmic Emu
         emulator_version = emulator_CosmicEmu
         emulator = 'CosmicEmu'
      ELSE IF (imode == 94) THEN
         ncos = 37 ! Franken Emu
         emulator_version = emulator_FrankenEmu
         emulator = 'FrankenEmu'
      ELSE IF (imode == 95) THEN
         ncos = 36 ! Mira Titan
         emulator_version = emulator_MiraTitan
         emulator = 'MiraTitan'
      ELSE 
         STOP 'BIG_EMULATOR_COMPARISON: Error, emulator_version not specified correctly'
      END IF

      ! Set number of redshifts
      IF (emulator_version == emulator_CosmicEmu) THEN
         nz = 3 ! Cosmic Emu
      ELSE IF (emulator_version == emulator_FrankenEmu) THEN
         nz = 6 ! Franken Emu
      ELSE IF (emulator_version == emulator_MiraTitan) THEN
         nz = 4 ! Mira Titan 
      ELSE
         STOP 'BIG_EMULATOR_COMPARISON:  Error, imode not specified correctly for EMU'
      END IF

      ! Allocate arrays
      na = nz
      ALLOCATE (z(nz), a(nz))

      ! Set redshifts
      IF (emulator_version == emulator_CosmicEmu) THEN
         ! Cosmic Emu (z up to 1)
         z(1) = 0.0
         z(2) = 0.5
         z(3) = 1.0
      ELSE IF (emulator_version == emulator_FrankenEmu) THEN
         ! Franken Emu (z up to 4)
         z(1) = 0.0
         z(2) = 0.5
         z(3) = 1.0
         z(4) = 2.0
         z(5) = 3.0
         z(6) = 4.0
      ELSE IF (emulator_version == emulator_MiraTitan) THEN
         ! Mira Titan (z up to 2)
         z(1) = 0.0
         z(2) = 0.5
         z(3) = 1.0
         z(4) = 2.0
      ELSE
         STOP 'BIG_EMULATOR_COMPARISON:  Error, imode not specified correctly for EMU'
      END IF

      ! Set scale factors
      a = scale_factor_z(z)

      ! Loop over cosmologies
      DO icos = 0, ncos
         
         IF (emulator_version == emulator_CosmicEmu) THEN
            icosmo = 300+icos ! Cosmic Emu nodes
         ELSE IF (emulator_version == emulator_FrankenEmu) THEN
            icosmo = 200+icos ! Franken Emu nodes
         ELSE IF (emulator_version == emulator_MiraTitan) THEN
            icosmo = 100+icos ! Mira Titan nodes 
         ELSE
            WRITE (*, *) 'BIG_EMULATOR_COMPARISON:', imode
            STOP 'BIG_EMULATOR_COMPARISON: Error, imode not specified correctly'
         END IF

         ! Assign and init the cosmology
         CALL assign_cosmology(icosmo, cosm, verbose)
         CALL init_cosmology(cosm)
         CALL print_cosmology(cosm)

         ! Get the emulator power
         CALL get_emulator_power(k, a, Pk_emu, nk, na, cosm, rebin, emulator_version)

         ! Loop over the different comparisons to be made
         DO icomp = 1, ncomp

            IF (icomp == 1) THEN
               CALL calculate_plin(k, a, Pk, nk, na, cosm)
               comp = 'linear'
            ELSE IF (icomp == 2) THEN
               CALL calculate_HALOFIT(k, a, Pk, nk, na, cosm, version=HALOFIT_Smith)
               comp = 'HALOFIT_Smith'
            ELSE IF (icomp == 3) THEN
               CALL calculate_HALOFIT(k, a, Pk, nk, na, cosm, version=HALOFIT_Bird)
               comp = 'HALOFIT_Bird'
            ELSE IF (icomp == 4) THEN
               CALL calculate_HALOFIT(k, a, Pk, nk, na, cosm, version=HALOFIT_Takahashi)
               comp = 'HALOFIT_Takahashi'
            ELSE IF (icomp == 5) THEN
               CALL calculate_HALOFIT(k, a, Pk, nk, na, cosm, version=HALOFIT_CAMB)
               comp = 'HALOFIT_CAMB'
            ELSE IF (icomp == 6) THEN
               CALL calculate_HMcode(k, a, Pk, nk, na, cosm, version=HMcode2015)
               comp = 'HMcode_2015'
            ELSE IF (icomp == 7) THEN
               CALL calculate_HMcode(k, a, Pk, nk, na, cosm, version=HMcode2016)
               comp = 'HMcode_2016'
            ELSE IF (icomp == 8) THEN
               CALL calculate_HMcode(k, a, Pk, nk, na, cosm, version=HMcode2016_neutrinofix)
               comp = 'HMcode_2016_neutrinofix'
            ELSE IF (icomp == 9) THEN
               CALL calculate_HMcode(k, a, Pk, nk, na, cosm, version=HMcode2019)
               comp = 'HMcode_2019'
            ELSE IF (icomp == 10) THEN
               CALL calculate_HMcode(k, a, Pk, nk, na, cosm, version=HMcode2020)
               comp = 'HMcode_2020'
            ELSE
               STOP 'BIG_EMULATOR_COMPARION: Error, icomp not supported'
            END IF

            base = 'data/power_'//trim(comp)//'_'//trim(emulator)//'_cos'//trim(integer_to_string(icos))

            ! Loop over redshift
            DO iz = 1, nz

               IF (z(iz) == 0.0) THEN
                  zlab = '0p0'
               ELSE IF (z(iz) == 0.5) THEN
                  zlab = '0p5'
               ELSE IF (z(iz) == 1.0) THEN
                  zlab = '1p0'
               ELSE IF (z(iz) == 2.0) THEN
                  zlab = '2p0'
               ELSE IF (z(iz) == 3.0) THEN
                  zlab = '3p0'
               ELSE IF (z(iz) == 4.0) THEN
                  zlab = '4p0'
               ELSE
                  STOP 'BIG_EMULATOR_COMPARISON: Error, z could not be converted into label'
               END IF

               ! Write data
               outfile = trim(base)//'_z'//trim(zlab)//'.dat'
               OPEN (7, file=outfile)
               DO ik = 1, nk
                  WRITE (7, *) k(ik), Pk(ik, iz), Pk_emu(ik, iz), cosm%m_nu
               END DO
               CLOSE (7)

            END DO

         END DO

      END DO

   END SUBROUTINE big_emulator_comparison

   SUBROUTINE power_breakdown_halomass(imode, icosmo, ihm)

      ! Power breakdown as a function of mass (paper)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: imode
      INTEGER, INTENT(INOUT) :: icosmo
      INTEGER, INTENT(INOUT) :: ihm
      REAL, ALLOCATABLE :: k(:)
      LOGICAL, PARAMETER :: verbose = .TRUE.
      REAL :: m1, m2
      REAL, ALLOCATABLE :: pow_li(:), pow_2h(:, :, :), pow_1h(:, :, :), pow_hm(:, :, :)
      INTEGER, ALLOCATABLE :: fields(:)
      INTEGER :: i, j1, j2
      CHARACTER(len=256) :: base, ext, mid, outfile
      TYPE(cosmology) :: cosm
      TYPE(halomod) :: hmod

      INTEGER, PARAMETER :: imin = 10
      INTEGER, PARAMETER :: imax = 16
      REAL, PARAMETER :: kmin = 1e-3
      REAL, PARAMETER :: kmax = 1e2
      INTEGER, PARAMETER :: nk = 128
      REAL, PARAMETER :: z = 0.0
      INTEGER, PARAMETER :: nf = 2

      ! Set number of k points and k range (log spaced)
      CALL fill_array(log(kmin), log(kmax), k, nk)
      k = exp(k)

      ! Set the fields
      ALLOCATE (fields(nf))
      fields(1) = field_matter
      fields(2) = field_electron_pressure

      ! Allocate arrays for power
      ALLOCATE (pow_li(nk), pow_2h(nf, nf, nk), pow_1h(nf, nf, nk), pow_hm(nf, nf, nk))

      ! Assign the cosmological model
      IF (imode == 31) icosmo = 4 ! 4 - BAHAMAS cosmology
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! Set the halo model
      IF (imode == 31) ihm = 65
      CALL assign_halomod(ihm, hmod, verbose)

      ! Loop over upper limit of mass integral
      DO i = imin, imax

         ! Set the upper limit for the mass integration
         ! Needs to be 10. to enforce that m2 is real
         IF (imode == 31 .OR. imode == 87) THEN
            m1 = 1e7
            m2 = 10.**i
         ELSE IF (imode == 53) THEN
            IF (i == imax) THEN
               m1 = 1e7
               m2 = 1e17
            ELSE
               m1 = 10.**i
               m2 = 10.**(i+1)
            END IF
         ELSE
            STOP 'HMx_DRIVER: Error, imode specified incorrectly'
         END IF
         hmod%mmin = m1
         hmod%mmax = m2

         ! Initiliasation for the halomodel calcualtion
         CALL init_halomod(scale_factor_z(z), hmod, cosm, verbose)
         CALL print_halomod(hmod, cosm, verbose)

         ! Do the halo-model calculation
         CALL calculate_HMx_a(fields, nf, k, nk, pow_li, pow_2h, pow_1h, pow_hm, hmod, cosm, verbose)

         ! Loop over fields
         DO j1 = 1, nf
            DO j2 = j1, nf

               ! Write out the results
               base = 'data/power_'
               IF (imode == 31 .OR. imode == 87) THEN
                  mid = ''
                  ext = '_m'
                  base = number_file2(base, fields(j1), mid, fields(j2), ext)
                  ext = '.dat'
                  outfile = number_file(base, i, ext)
               ELSE IF (imode == 53) THEN
                  IF (i == imax) THEN
                     mid = ''
                     ext = '.dat'
                     outfile = number_file2(base, fields(j1), mid, fields(j2), ext)
                  ELSE
                     mid = ''
                     ext = '_m'
                     base = number_file2(base, fields(j1), mid, fields(j2), ext)
                     mid = '_m'
                     ext = '.dat'
                     outfile = number_file2(base, i, mid, i+1, ext)
                  END IF
               ELSE
                  STOP 'HMx_DRIVER: Error, imode specified incorrectly'
               END IF

               CALL write_power(k, pow_li, pow_2h(j1, j2, :), pow_1h(j1, j2, :), pow_hm(j1, j2, :), nk, outfile, verbose)

            END DO
         END DO

      END DO

   END SUBROUTINE power_breakdown_halomass

   SUBROUTINE write_HMx_variations(icosmo, ihm)

      ! Write out the variation of HMx parameters with T_AGN and z
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: icosmo
      INTEGER, INTENT(INOUT) :: ihm
      INTEGER :: i, j
      CHARACTER(len=256) :: outfile
      TYPE(cosmology) :: cosm
      TYPE(halomod) :: hmod

      REAL, PARAMETER :: z = 0.0
      REAL, PARAMETER :: zmin = 0.
      REAL, PARAMETER :: zmax = 4.
      REAL, PARAMETER :: mass = 1e14 ! Mass to write out for [Msun/h]
      INTEGER, PARAMETER :: nz = 101
      LOGICAL, PARAMETER :: verbose = .TRUE.

      !Assign the cosmological model
      icosmo = 4
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      !Initiliasation for the halomodel calcualtion
      ihm = 18
      CALL assign_halomod(ihm, hmod, verbose)
      CALL init_halomod(scale_factor_z(z), hmod, cosm, verbose)
      CALL print_halomod(hmod, cosm, verbose)

      DO j = 1, 3

         IF (j == 1) THEN
            outfile = 'data/HMx_params_AGN_7p6_nu0.dat'
            hmod%Theat = 10**7.6
         ELSE IF (j == 2) THEN
            outfile = 'data/HMx_params_AGN_TUNED_nu0.dat'
            hmod%Theat = 10**7.8
         ELSE IF (j == 3) THEN
            outfile = 'data/HMx_params_AGN_8p0_nu0.dat'
            hmod%Theat = 10**8.0
         END IF

         OPEN (7, file=outfile)
         DO i = 1, nz
            hmod%z = progression(zmin, zmax, i, nz)
            WRITE (7, *) hmod%z, &
               HMx_alpha(mass, hmod, cosm), &
               HMx_eps(hmod, cosm), &
               HMx_Gamma(mass, hmod, cosm), &
               HMx_M0(hmod, cosm), &
               HMx_Astar(hmod, cosm), &
               HMx_Twhim(hmod, cosm)
         END DO
         CLOSE (7)

      END DO

   END SUBROUTINE write_HMx_variations

   SUBROUTINE halo_cores(icosmo, ihm)

      ! Look at the effect of cores on halo profiles
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: icosmo
      INTEGER, INTENT(INOUT) :: ihm
      REAL, ALLOCATABLE :: k(:)
      LOGICAL :: verbose2
      INTEGER :: field(1)
      INTEGER :: i
      REAL, ALLOCATABLE :: pow_li(:), pow_2h(:, :, :), pow_1h(:, :, :), pow_hm(:, :, :)
      CHARACTER(len=256) :: base, ext, outfile
      TYPE(cosmology) :: cosm
      TYPE(halomod) :: hmod

      REAL, PARAMETER :: kmin = 1e-3
      REAL, PARAMETER :: kmax = 1e2
      INTEGER, PARAMETER :: nk = 128
      REAL, PARAMETER :: z = 0.0
      REAL, PARAMETER :: rcore_min = 0.
      REAL, PARAMETER :: rcore_max = 0.1
      INTEGER, PARAMETER :: ncore = 16
      LOGICAL, PARAMETER :: verbose = .TRUE.

      ! Set number of k points and k range (log spaced)
      CALL fill_array(log(kmin), log(kmax), k, nk)
      k = exp(k)

      ! Allocate arrays for the power calculation
      ALLOCATE (pow_li(nk), pow_2h(1, 1, nk), pow_1h(1, 1, nk), pow_hm(1, 1, nk))

      ! Assign the cosmological model
      icosmo = 1 ! 1 - Boring
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! Initiliasation for the halomodel calcualtion
      ihm = 21
      CALL assign_halomod(ihm, hmod, verbose)

      ! Range of cores to explore
      DO i = 1, ncore

         ! Set the core radius
         hmod%rcore = progression(rcore_min, rcore_max, i, ncore)

         ! Initialise the halo-model calculation
         IF (i == 1) THEN
            verbose2 = verbose
         ELSE
            verbose2 = .FALSE.
         END IF
         CALL init_halomod(scale_factor_z(z), hmod, cosm, verbose2)
         CALL print_halomod(hmod, cosm, verbose2)

         ! Do the halo-model calculation
         field = field_dmonly
         CALL calculate_HMx_a(field, 1, k, nk, pow_li, pow_2h, pow_1h, pow_hm, hmod, cosm, verbose=.FALSE.)

         ! Write out the results
         base = 'data/power_cored_'
         ext = '.dat'
         outfile = number_file(base, i, ext)
         CALL write_power(k, pow_li, pow_2h(1, 1, :), pow_1h(1, 1, :), pow_hm(1, 1, :), nk, outfile, verbose)

      END DO

   END SUBROUTINE halo_cores

   SUBROUTINE hydro_tests(icosmo, ihm)

      ! Automated testing of hydro models
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: icosmo
      INTEGER, INTENT(INOUT) :: ihm
      REAL, ALLOCATABLE :: k(:), a(:)
      REAL :: error, error_max
      LOGICAL :: verbose_tests
      REAL :: spam
      REAL :: z
      INTEGER :: i, j, itest, jtest
      INTEGER :: nk
      REAL, ALLOCATABLE :: pow_li(:, :), pow_2h(:, :, :, :), pow_1h(:, :, :, :), pow_hm(:, :, :, :), powb_hm(:, :, :, :)
      INTEGER, ALLOCATABLE :: fields(:)
      CHARACTER(len=256) :: inbase, inext, infile, mid, outbase, outext, outfile
      TYPE(cosmology) :: cosm
      !TYPE(halomod) :: hmod

      REAL, PARAMETER :: kmin_test = 1e-3
      REAL, PARAMETER :: kmax_test = 10.
      REAL, PARAMETER :: tolerance = 3e-3
      INTEGER, PARAMETER :: nf = 5 ! Number of fields
      INTEGER, PARAMETER :: na = 4 ! Number of redshifts
      LOGICAL :: ifail = .FALSE.
      LOGICAL, PARAMETER :: verbose = .TRUE.

      ! Assigns the cosmological model
      icosmo = 4
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! Set the field combinations
      ALLOCATE (fields(nf))
      fields(1) = field_matter
      fields(2) = field_cdm
      fields(3) = field_gas
      fields(4) = field_stars
      fields(5) = field_electron_pressure

      ! Set the redshifts
      ALLOCATE (a(na))
      DO i = 1, na
         IF (i == 1) z = 0.0
         IF (i == 2) z = 0.5
         IF (i == 3) z = 1.0
         IF (i == 4) z = 2.0
         a(i) = scale_factor_z(z)
      END DO

      ! File naming things
      mid = ''
      inext = '.txt'

      !! Read benchmark data

      ! Loop over fields
      DO itest = 1, nf
         DO jtest = itest, nf

            ! Loop over redshifts
            DO j = 1, na

               ! Set the redshift and input and output file bases
               IF (j == 1) THEN
                  inbase = 'benchmarks/power_z0.0_'
               ELSE IF (j == 2) THEN
                  inbase = 'benchmarks/power_z0.5_'
               ELSE IF (j == 3) THEN
                  inbase = 'benchmarks/power_z1.0_'
               ELSE IF (j == 4) THEN
                  inbase = 'benchmarks/power_z2.0_'
               ELSE
                  STOP 'HMX_DRIVER: Error, iz specified incorrectly'
               END IF

               ! Input file name
               infile = number_file2(inbase, fields(itest), mid, fields(jtest), inext)

               ! Allocate arrays
               IF (itest == 1 .AND. jtest == 1 .AND. j == 1) THEN
                  nk = file_length(infile, verbose=.FALSE.)
                  ALLOCATE (k(nk), powb_hm(nf, nf, nk, na))
               END IF

               ! Loop over k and read in benchmark data
               ! TODO: Add tests for one- and two-halo terms individually
               OPEN (7, file=infile, status='old')
               DO i = 1, nk
                  READ (7, *) k(i), spam, spam, spam, powb_hm(itest, jtest, i, j)
               END DO
               CLOSE (7)

            END DO
         END DO
      END DO

      !!

      ! Assign the halo model
      ihm = 3
      !CALL assign_halomod(ihm, hmod, verbose)
      !CALL calculate_HMx(fields, nf, k, nk, a, na, pows_li, pows_2h, pows_1h, pows_hm, hmod, cosm, verbose)
      CALL calculate_HMx_full(fields, k, a, pow_li, pow_2h, pow_1h, pow_hm, nf, nk, na, cosm, ihm)

      ! File naming things
      outext = '.dat'

      ! Loop over fields
      DO itest = 1, nf
         DO jtest = itest, nf

            ! Loop over redshift
            DO j = 1, na

               ! Set the redshift and input and output file bases
               IF (j == 1) THEN
                  outbase = 'data/power_z0.0_'
               ELSE IF (j == 2) THEN
                  outbase = 'data/power_z0.5_'
               ELSE IF (j == 3) THEN
                  outbase = 'data/power_z1.0_'
               ELSE IF (j == 4) THEN
                  outbase = 'data/power_z2.0_'
               ELSE
                  STOP 'HMX_DRIVER: Error, iz specified incorrectly'
               END IF

               ! Write data to file
               outfile = number_file2(outbase, fields(itest), mid, fields(jtest), outext)
               CALL write_power(k, pow_li(:, j), pow_2h(itest, jtest, :, j), pow_1h(itest, jtest, :, j), pow_hm(itest, jtest, :, j), nk, outfile, verbose)

               ! Set the error max to be zero before each test
               error_max = 0.
               verbose_tests = .TRUE.

               ! Loop over k and check for accuracy and write out data
               DO i = 1, nk

                  ! This is what I use as the error
                  error = ABS(-1.+pow_hm(itest, jtest, i, j)/powb_hm(itest, jtest, i, j))

                  IF (kmin_test <= k(i) .AND. k(i) <= kmax_test) THEN

                     ! Update the maximium error if it is exceeded
                     IF (error > error_max) error_max = error

                     ! If the test has failed then write out this diagnostic stuff
                     IF (verbose_tests .AND. error > tolerance) THEN
                        WRITE (*, *) 'HMx_DRIVER: Test failing'
                        WRITE (*, *) 'HMx_DRIVER: Test fields:', fields(itest), fields(jtest)
                        WRITE (*, *) 'HMx_DRIVER: Wavenumber [h/Mpc]:', k(i)
                        WRITE (*, *) 'HMx_DRIVER: Redshift:', redshift_a(a(j))
                        WRITE (*, *) 'HMx_DRIVER: Benchmark power:', powb_hm(itest, jtest, i, j)
                        WRITE (*, *) 'HMx_DRIVER: HMx power:', pow_hm(itest, jtest, i, j)
                        WRITE (*, *) 'HMx_DRIVER: Tolerance:', tolerance
                        WRITE (*, *) 'HMx_DRIVER: Error:', error
                        WRITE (*, *)
                        verbose_tests = .FALSE.
                        ifail = .TRUE.
                     END IF
                     
                  END IF

               END DO

               ! Write test information to the screen
               WRITE (*, *) 'HMx_DRIVER: Test fields:', fields(itest), fields(jtest)
               WRITE (*, *) 'HMx_DRIVER: Redshift:', redshift_a(a(j))
               WRITE (*, *) 'HMx_DRIVER: Tolerance:', tolerance
               WRITE (*, *) 'HMx_DRIVER: Max error:', error_max
               WRITE (*, *)

            END DO

         END DO
      END DO

      ! Write pass/fail information to the screen
      IF (ifail) THEN
         WRITE (*, *) 'HMx_DRIVER: Hydro tests failed'
      ELSE
         WRITE (*, *) 'HMx_DRIVER: Hydro tests should take around 2.16 seconds to run'
         WRITE (*, *) 'HMx_DRIVER: Hydro tests passed'
      END IF
      WRITE (*, *)

   END SUBROUTINE hydro_tests

   SUBROUTINE power_cosmology_dependence(imode, icosmo, ihm)

      ! Assess 3D power as a function of cosmology
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: imode
      INTEGER, INTENT(INOUT) :: icosmo
      INTEGER, INTENT(INOUT) :: ihm
      REAL, ALLOCATABLE :: k(:)
      INTEGER, ALLOCATABLE :: fields(:)
      INTEGER :: i
      REAL, ALLOCATABLE :: pow_li(:), pow_2h(:, :, :), pow_1h(:, :, :), pow_hm(:, :, :)
      CHARACTER(len=256) :: dir, base, ext, outfile
      TYPE(cosmology) :: cosm
      TYPE(halomod) :: hmod

      INTEGER, PARAMETER :: nf = 2
      REAL, PARAMETER :: kmin = 1e-3
      REAL, PARAMETER :: kmax = 1e2
      INTEGER, PARAMETER :: nk = 128
      REAL, PARAMETER :: z = 0.0  ! Set the redshift
      REAL, PARAMETER :: sig8min = 0.7 ! Set range in sigma_8
      REAL, PARAMETER :: sig8max = 0.9 ! Set range in sigma_8
      INTEGER, PARAMETER :: ncos = 15
      LOGICAL, PARAMETER :: verbose = .TRUE.

      ! Assign the cosmology
      IF(imode == 41) icosmo = 44
      CALL assign_cosmology(icosmo, cosm, verbose)

      ! Assign the halo model
      IF(imode == 41) ihm = 65
      CALL assign_halomod(ihm, hmod, verbose=.FALSE.)

      ! Allocate arrays for the fields
      ALLOCATE (fields(nf))
      fields(1) = field_matter
      fields(2) = field_electron_pressure

      ! Set number of k points and k range (log spaced)

      CALL fill_array(log(kmin), log(kmax), k, nk)
      k = exp(k)

      ! Allocate arrays for the power
      ALLOCATE (pow_li(nk), pow_2h(nf, nf, nk), pow_1h(nf, nf, nk), pow_hm(nf, nf, nk))

      ! Loop over cosmology
      DO i = 1, ncos

         ! Change the sigma_8 value
         cosm%sig8 = progression(sig8min, sig8max, i, ncos)
         CALL init_cosmology(cosm)
         CALL print_cosmology(cosm)

         ! Initiliasation for the halomodel calcualtion
         CALL init_halomod(scale_factor_z(z), hmod, cosm, verbose=.FALSE.)
         CALL print_halomod(hmod, cosm, verbose=.FALSE.)

         ! Do the halo-model calculation
         CALL calculate_HMx_a(fields, 2, k, nk, pow_li, pow_2h, pow_1h, pow_hm, hmod, cosm, verbose=.FALSE.)

         ! Write data
         dir = 'data/'
         base = TRIM(dir)//'cosmology_'
         ext = '_matter-matter_power.dat'
         outfile = number_file(base, i, ext)
         CALL write_power(k, pow_li, pow_2h(1, 1, :), pow_1h(1, 1, :), pow_hm(1, 1, :), nk, outfile, verbose=.TRUE.)
         ext = '_matter-epressure_power.dat'
         outfile = number_file(base, i, ext)
         CALL write_power(k, pow_li, pow_2h(1, 2, :), pow_1h(1, 2, :), pow_hm(1, 2, :), nk, outfile, verbose=.TRUE.)

      END DO

   END SUBROUTINE power_cosmology_dependence

   SUBROUTINE power_different_halo_mass_functions(icosmo, ihm)

      ! Comparison of power spectra from Sheth & Tormen (1999) vs. Tinker et al. (2010) mass function
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: icosmo
      INTEGER, INTENT(INOUT) :: ihm
      REAL, ALLOCATABLE :: k(:)
      INTEGER :: field(1)
      INTEGER :: j
      LOGICAL :: verbose2
      REAL, ALLOCATABLE :: pow_li(:), pow_2h(:, :, :), pow_1h(:, :, :), pow_hm(:, :, :)
      CHARACTER(len=256) :: outfile
      TYPE(cosmology) :: cosm
      TYPE(halomod) :: hmod

      REAL, PARAMETER :: kmin = 1e-3
      REAL, PARAMETER :: kmax = 1e2
      INTEGER, PARAMETER :: nk = 128
      INTEGER, PARAMETER :: nhm = 4
      REAL, PARAMETER :: z = 0.0
      LOGICAL, PARAMETER :: verbose = .TRUE.

      ! Set number of k points and k range (log spaced)
      CALL fill_array(log(kmin), log(kmax), k, nk)
      k = exp(k)

      ! Allocate arrays for power
      ALLOCATE (pow_li(nk), pow_2h(1, 1, nk), pow_1h(1, 1, nk), pow_hm(1, 1, nk))

      ! Assigns the cosmological model
      icosmo = 1
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      DO j = 1, nhm

         IF (j == 1) THEN
            ! Sheth & Tormen (1999) mass function
            ihm = 3
            outfile = 'data/power_ShethTormen.dat'
         ELSE IF (j == 2) THEN
            ihm = 23 ! Tinker et al. (2010) mass function
            outfile = 'data/power_Tinker.dat'
         ELSE IF (j == 3) THEN
            ihm = 27 ! Press & Schecter (1974) mass function
            outfile = 'data/power_PressSchecter.dat'
         ELSE IF (j == 4) THEN
            ihm = 80 ! Jenkins mass function
            outfile = 'data/power_Jenkins.dat'
         ELSE
            STOP 'POWER_DIFFERENT_MASS_FUNCTIONS: Something went wrong'
         END IF

         IF (j == 1) THEN
            verbose2 = verbose
         ELSE
            verbose2 = .FALSE.
         END IF

         !Initiliasation for the halomodel calcualtion
         CALL assign_halomod(ihm, hmod, verbose2)
         CALL init_halomod(scale_factor_z(z), hmod, cosm, verbose2)
         CALL print_halomod(hmod, cosm, verbose2)

         ! Do the halo-model calculation
         field = field_dmonly
         CALL calculate_HMx_a(field, 1, k, nk, pow_li, pow_2h, pow_1h, pow_hm, hmod, cosm, verbose2)

         ! Write out data
         CALL write_power(k, pow_li, pow_2h(1, 1, :), pow_1h(1, 1, :), pow_hm(1, 1, :), nk, outfile, verbose)

      END DO

   END SUBROUTINE power_different_halo_mass_functions

   SUBROUTINE mass_function_plots(icosmo, ihm)

      ! Mass function plots for different mass functions
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: icosmo
      INTEGER, INTENT(INOUT) :: ihm
      REAL :: nu
      REAL, ALLOCATABLE :: m(:)
      INTEGER :: i, j
      TYPE(cosmology) :: cosm
      TYPE(halomod), ALLOCATABLE :: hmods(:)

      REAL, PARAMETER :: z = 0.0 ! Sets the redshift
      INTEGER, PARAMETER :: nhm = 4 ! Number of halo models
      LOGICAL, PARAMETER :: verbose = .TRUE.
      REAL, PARAMETER :: numin = 0.1  ! Minimum value of nu
      REAL, PARAMETER :: numax = 6.   ! Maximum value of nu
      INTEGER, PARAMETER :: nnu = 256 ! Number of nu values

      ! Assigns the cosmological model
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      !Initiliasation for the halomodel calcualtion
      ALLOCATE (hmods(nhm), m(nhm))
      DO i = 1, nhm
         IF (i == 1) THEN
            ihm = 3  ! 3 - Sheth & Tormen mass function, virial mass haloes
         ELSE IF (i == 2) THEN
            ihm = 27 ! 27 - Press & Schecter mass function, virial mass haloes
         ELSE IF (i == 3) THEN
            ihm = 23 ! 23 - Tinker mass function, virial mass haloes
         ELSE IF (i == 4) THEN
            ihm = 80 ! 80 - Jenkins (2001) FoF b = 0.2 haloes
         ELSE
            STOP 'HMx_DRIVER: Error, something went wrong setting halo model'
         END IF
         CALL assign_halomod(ihm, hmods(i), verbose)
         CALL init_halomod(scale_factor_z(z), hmods(i), cosm, verbose)
         CALL print_halomod(hmods(i), cosm, verbose)
      END DO

      ! Loop over nu and write out mass function
      OPEN (10, file='data/bnu_functions.dat')
      OPEN (11, file='data/gnu_functions.dat')
      OPEN (12, file='data/mass_functions.dat')
      OPEN (13, file='data/multiplicity_functions.dat')
      DO i = 1, nnu
         nu = progression(numin, numax, i, nnu)
         DO j = 1, nhm
            m(j) = M_nu(nu, hmods(j))
         END DO
         IF ((m(1) .NE. m(2)) .OR. (m(1) .NE. m(3))) THEN
            STOP 'A disaster has occured with halo mass'
         END IF
         WRITE (10, *) nu, m(1), (b_nu(nu, hmods(j)), j=1, nhm)
         WRITE (11, *) nu, m(1), (g_nu(nu, hmods(j)), j=1, nhm)
         WRITE (12, *) nu, m(1), (mass_function(m(j), hmods(j), cosm), j=1, nhm)
         WRITE (13, *) nu, m(1), (multiplicity_function(m(j), hmods(j), cosm), j=1, nhm)
      END DO
      CLOSE (10)
      CLOSE (11)
      CLOSE (12)
      CLOSE (13)

   END SUBROUTINE mass_function_plots

   SUBROUTINE HI_mass_fractions(icosmo, ihm)

      ! HI mass fractions
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: icosmo
      INTEGER, INTENT(INOUT) :: ihm
      REAL, ALLOCATABLE :: a(:), z(:)
      INTEGER :: i, j
      REAL :: m
      REAL, ALLOCATABLE :: HI_frac(:)
      TYPE(cosmology) :: cosm
      TYPE(halomod), ALLOCATABLE :: hmod(:)

      REAL, PARAMETER :: z1 = 0.
      REAL, PARAMETER :: z2 = 5.
      INTEGER, PARAMETER :: nz = 6
      REAL, PARAMETER :: m1 = 1e7
      REAL, PARAMETER :: m2 = 1e17
      INTEGER, PARAMETER :: nm = 256
      LOGICAL, PARAMETER :: verbose = .TRUE.

      ! Assigns the cosmological model
      icosmo = 1 ! 1 - Boring
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! Set the range of redshifts
      ALLOCATE (hmod(nz), HI_frac(nz), z(nz), a(nz))

      ! Initiliasation for the halomodel calcualtion at all z
      ihm = 25
      DO j = 1, nz
         z(j) = progression(z1, z2, j, nz)
         a(j) = scale_factor_z(z(j))
         CALL assign_halomod(ihm, hmod(j), verbose)
         CALL init_halomod(a(j), hmod(j), cosm, verbose)
         CALL print_halomod(hmod(j), cosm, verbose)
      END DO

      ! Loop over mass and z and do the calculation
      OPEN (7, file='data/HI_mass_fraction.dat')
      DO i = 1, nm
         m = exp(progression(log(m1), log(m2), i, nm))
         DO j = 1, nz
            HI_frac(j) = halo_HI_fraction(m, hmod(j), cosm)
         END DO
         WRITE (7, *) m, (HI_frac(j), j=1, nz)
      END DO
      CLOSE (7)
      CLOSE (8)

   END SUBROUTINE HI_mass_fractions

   SUBROUTINE mass_function_Lbox(icosmo, ihm)

      ! Mass function plots as Lbox is varied
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: icosmo
      INTEGER, INTENT(INOUT) :: ihm
      LOGICAL :: verbose2
      REAL :: nu
      REAL, ALLOCATABLE :: m(:)
      INTEGER :: i, j
      TYPE(cosmology) :: cosm
      TYPE(cosmology), ALLOCATABLE :: cosms(:)
      TYPE(halomod), ALLOCATABLE :: hmods(:)

      REAL, PARAMETER :: numin = 0.1
      REAL, PARAMETER :: numax = 6.
      INTEGER, PARAMETER :: nnu = 256
      REAL, PARAMETER :: Lmin = 32.
      REAL, PARAMETER :: Lmax = 2048.
      INTEGER, PARAMETER :: ncos = 8
      REAL, PARAMETER :: z = 0.0 ! Set the redshift
      LOGICAL, PARAMETER :: verbose = .TRUE.

      ! Assigns the cosmological model
      icosmo = 1
      CALL assign_cosmology(icosmo, cosm, verbose=.FALSE.)

      ! Minimum and maximum box sizes [Mpc/h] and number of cosmologies

      ! Set cosmological models
      ALLOCATE (cosms(ncos))
      cosms = cosm
      DO i = 1, ncos
         IF (i .NE. 1) THEN
            cosms(i)%box = .TRUE.
         END IF
         cosms(i)%Lbox = 2**(12-i)
         CALL init_cosmology(cosms(i))
      END DO
      CALL print_cosmology(cosms(1))

      !Initilisation for the halomodel calcualtion
      ihm = 3
      ALLOCATE (hmods(ncos))
      DO i = 1, ncos
         IF (i == 1) THEN
            verbose2 = verbose
         ELSE
            verbose2 = .FALSE.
         END IF
         CALL assign_halomod(ihm, hmods(i), verbose2)
         CALL init_halomod(scale_factor_z(z), hmods(i), cosms(i), verbose2)
      END DO
      CALL print_halomod(hmods(1), cosms(1), verbose)

      ALLOCATE (m(ncos))
      OPEN (7, file='data/nu_mass_Lbox.dat')
      DO i = 1, nnu
         DO j = 1, ncos
            nu = progression(numin, numax, i, nnu)
            m(j) = M_nu(nu, hmods(j))
         END DO
         WRITE (7, *) nu, (m(j), j=1, ncos)
      END DO
      CLOSE (7)

   END SUBROUTINE mass_function_Lbox

   SUBROUTINE power_scatter(icosmo, ihm)

      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: icosmo
      INTEGER, INTENT(INOUT) :: ihm
      REAL, ALLOCATABLE :: a(:), k(:)
      REAL, ALLOCATABLE :: pow_li(:, :), pow_2h(:, :, :, :), pow_1h(:, :, :, :), pow_hm(:, :, :, :)
      REAL :: c
      INTEGER :: i
      INTEGER :: field(1)
      CHARACTER(len=256) :: base
      TYPE(cosmology) :: cosm
      TYPE(halomod) :: hmod

      REAL, PARAMETER :: kmin = 1e-3
      REAL, PARAMETER :: kmax = 1e2
      INTEGER, PARAMETER :: nk = 128
      REAL, PARAMETER :: amin = 0.1
      REAL, PARAMETER :: amax = 1.
      INTEGER, PARAMETER :: na = 16
      REAL, PARAMETER :: cmin = 1.
      REAL, PARAMETER :: cmax = 10.
      INTEGER, PARAMETER :: n = 256
      REAL, PARAMETER :: cbar = 4.
      LOGICAL, PARAMETER :: verbose = .TRUE.
      CHARACTER(len=256), PARAMETER :: outfile = 'data/p_conc.dat'

      ! Assigns the cosmological model
      icosmo = 1
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! Set number of k points and k range (log spaced)
      CALL fill_array(log(kmin), log(kmax), k, nk)
      k = exp(k)

      ! Set the number of redshifts and range (linearly spaced) and convert z -> a

      CALL fill_array(amin, amax, a, na)

      ! Set the field
      field = field_dmonly

      DO i = 1, 2

         ! Assign halo model
         IF (i == 1) THEN
            ihm = 3
            base = 'data/power'
         ELSE IF (i == 2) THEN
            ihm = 8
            base = 'data/power_scatter'
         ELSE
            STOP 'HMX_DRIVER: Error, something went wrong'
         END IF

         !CALL assign_halomod(ihm, hmod, verbose)
         !CALL calculate_HMx(field, 1, k, nk, a, na, pows_li, pows_2h, pows_1h, pows_hm, hmod, cosm, verbose)
         CALL calculate_HMx_full(field, k, a, pow_li, pow_2h, pow_1h, pow_hm, 1, nk, na, cosm, ihm)

         CALL write_power_a_multiple(k, a, pow_li, pow_2h(1, 1, :, :), pow_1h(1, 1, :, :), pow_hm(1, 1, :, :), nk, na, base, verbose)

      END DO

      !outfile = 'data/p_conc.dat'
      CALL assign_halomod(ihm, hmod, verbose)
      OPEN (7, file=outfile)
      DO i = 1, n
         c = progression(cmin, cmax, i, n)
         WRITE (7, *) c, lognormal(c, cbar, hmod%dlnc)
      END DO
      CLOSE (7)

   END SUBROUTINE power_scatter

   SUBROUTINE Trispectrum_test(icosmo, ihm)

      ! Trispectrum test
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: icosmo
      INTEGER, INTENT(INOUT) :: ihm
      REAL, ALLOCATABLE :: k(:)
      INTEGER, ALLOCATABLE :: fields(:)
      INTEGER :: i, j
      TYPE(cosmology) :: cosm
      TYPE(halomod) :: hmod

      REAL, PARAMETER :: z = 0.       ! Set the redshift
      REAL, PARAMETER :: kmin = 1e-3
      REAL, PARAMETER :: kmax = 1e2
      INTEGER, PARAMETER :: nk = 64
      LOGICAL, PARAMETER :: verbose = .TRUE.

      ! Set number of k points and k range (log spaced)
      CALL fill_array(log(kmin), log(kmax), k, nk)
      k = k*real(nk)/real(nk+1) ! To stop plot spilling over border
      k = exp(k)

      ! Set the fields
      ALLOCATE (fields(2))
      fields = field_dmonly

      ! Assigns the cosmological model
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! Assign the halo model
      CALL assign_halomod(ihm, hmod, verbose)
      CALL init_halomod(scale_factor_z(z), hmod, cosm, verbose)
      CALL print_halomod(hmod, cosm, verbose)

      ! Write data
      OPEN (7, file='data/trispectrum.dat')
      DO j = 1, nk
         DO i = 1, nk
            WRITE (7, *) k(i), k(j), T_1h(k(i), k(j), fields, hmod, cosm)
         END DO
      END DO
      CLOSE (7)

   END SUBROUTINE Trispectrum_test

   SUBROUTINE Cl_direct_integration(imode, icosmo)

      ! Calculate C(l) by direct integration of the measured 3D BAHAMAS power spectra
      ! 57 - Triad 3
      ! 61 - Triad 4
      ! 62 - Triad 5
      ! 65 - Triad 5 with contribution to Limber integrals per ell
      ! 66 - Triad 5 but with extended ell range
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: imode
      INTEGER, INTENT(INOUT) :: icosmo
      REAL, ALLOCATABLE :: a(:), k(:), pow_ka(:, :)
      REAL :: lmin, lmax
      INTEGER :: i, j, ii, jj
      INTEGER :: nl, ix(2), ip(2), n_edge, n_all, nt, nk, na
      REAL, ALLOCATABLE :: pow_sim(:), err_sim(:)
      REAL, ALLOCATABLE :: ell(:), ell_edge(:), all_ell(:), Cl(:, :, :), all_Cl(:, :, :)
      LOGICAL :: verbose2
      CHARACTER(len=256) :: base, name, outbase, ext, outfile
      CHARACTER(len=256), ALLOCATABLE :: ixx_names(:)
      TYPE(cosmology) :: cosm

      LOGICAL, PARAMETER :: bin_theory = .FALSE.     ! Should the theory be binned in the same way as measurements?
      LOGICAL, PARAMETER :: cut_nyquist = .TRUE.     ! Should the BAHAMAS measured P(k) be cut above the Nyquist frequency?
      LOGICAL, PARAMETER :: subtract_shot = .TRUE.   ! Should the BAHAMAS measured P(k) have shot-noise subtracted?
      LOGICAL, PARAMETER :: response_triad = .FALSE. ! Should I treat the BAHAMAS P(k) as HMcode response?
      LOGICAL, PARAMETER :: add_highz = .FALSE.      ! Add in z=3 power
      INTEGER, PARAMETER :: mesh = 1024              ! Mesh size for BAHAMAS P(k) measurements
      INTEGER, PARAMETER :: nz_BAHAMAS = 15          ! Number of BAHAMAS redshift slices to use (4, 11 or 15)
      REAL, PARAMETER :: kmax_BAHAMAS = 1000.        ! Maximum k to trust the BAHAMAS P(k) [h/Mpc]
      INTEGER, PARAMETER :: nsim = 5
      LOGICAL, PARAMETER :: verbose = .TRUE.

      ! Assigns the cosmological model
      icosmo = 4 ! 4 - WMAP9
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      IF (imode == 57 .OR. imode == 61 .OR. imode == 62 .OR. imode == 65) THEN
         ! Triad 3 ell values
         CALL triad_ell(ell, nl)
      ELSE IF (imode == 66) THEN
         ! Another ell range
         lmin = 1e2
         lmax = 1e4
         nl = 128
         CALL fill_array(log(lmin), log(lmax), ell, nl)
         ell = exp(ell)
      ELSE
         STOP 'HMx_DRIVER: Error, something went wrong setting ell range'
      END IF
      ALLOCATE (Cl(nl, 2, 2))

      IF (bin_theory) THEN

         ! Triad 3,4,5 ell bin-edge values
         CALL triad_ell_edges(ell_edge, n_edge)

         ! Get range for all ell
         lmin = ell_edge(1)
         lmax = ell_edge(nl+1)
         n_all = NINT(lmax-lmin+1)
         CALL fill_array(lmin, lmax, all_ell, n_all)
         ALLOCATE (all_Cl(n_all, 2, 2))

      END IF

      ! Allocate array for cross-correlation names
      ALLOCATE (ixx_names(2))

      ! Number of tracers
      IF (imode == 57) THEN
         nt = 4
      ELSE IF (imode == 61) THEN
         nt = 5
      ELSE IF (imode == 62 .OR. imode == 65 .OR. imode == 66) THEN
         nt = 6
      ELSE
         STOP 'HMX_DRIVER: Error, imode specified incorrectly'
      END IF

      ! Loop over BAHAMAS models
      DO j = 1, nsim

         ! Set BAHAMAS models
         IF (j == 1) THEN
            name = 'AGN_TUNED_nu0'
            outbase = 'data/triad_Cl_direct_AGN_TUNED_nu0'
         ELSE IF (j == 2) THEN
            name = 'AGN_7p6_nu0'
            outbase = 'data/triad_Cl_direct_AGN_7p6_nu0'
         ELSE IF (j == 3) THEN
            name = 'AGN_8p0_nu0'
            outbase = 'data/triad_Cl_direct_AGN_8p0_nu0'
         ELSE IF (j == 4) THEN
            name = 'AGN_TUNED_nu0_v2'
            outbase = 'data/triad_Cl_direct_AGN_TUNED_nu0_v2'
         ELSE IF (j == 5) THEN
            name = 'AGN_TUNED_nu0_v3'
            outbase = 'data/triad_Cl_direct_AGN_TUNED_nu0_v3'
         ELSE
            STOP 'HMx_DRIVER: Error, feedback senario not supported'
         END IF

         ! Loop over first tracer
         verbose2 = verbose
         DO ii = 1, nt

            IF (imode == 57) THEN
               CALL triad_3_tracers(ii, ix(1), ixx_names(1))
            ELSE IF (imode == 61) THEN
               CALL triad_4_tracers(ii, ix(1), ixx_names(1))
            ELSE IF (imode == 62 .OR. imode == 65 .OR. imode == 66) THEN
               CALL triad_5_tracers(ii, ix(1), ixx_names(1))
            ELSE
               STOP 'HMX_DRIVER: Error, imode specified incorrectly'
            END IF

            ! Loop over second tracer
            DO jj = 1, nt

               IF (imode == 57) THEN
                  CALL triad_3_tracers(jj, ix(2), ixx_names(2))
               ELSE IF (imode == 61) THEN
                  CALL triad_4_tracers(jj, ix(2), ixx_names(2))
               ELSE IF (imode == 62 .OR. imode == 65 .OR. imode == 66) THEN
                  CALL triad_5_tracers(jj, ix(2), ixx_names(2))
               ELSE
                  STOP 'HMX_DRIVER: Error, imode specified incorrectly'
               END IF

               ! Use the xpowlation type to set the necessary halo profiles
               DO i = 1, 2
                  CALL set_field_for_xpow(ix(i), ip(i))
               END DO

               ! Choose the set of redshifts to use
               na = nz_BAHAMAS
               CALL BAHAMAS_scale_factors(a, na)

               ! Read in power
               DO i = 1, na
                  CALL BAHAMAS_read_power(k, pow_sim, err_sim, nk, redshift_a(a(i)), name, mesh, ip, cosm, &
                                          kmax=kmax_BAHAMAS, &
                                          cut_nyquist=cut_nyquist, &
                                          subtract_shot=subtract_shot, &
                                          response=response_triad, &
                                          verbose=verbose2)
                  verbose2 = .FALSE.
                  IF (i == 1) THEN
                     ALLOCATE (pow_ka(nk, na))
                  END IF
                  pow_ka(:, i) = pow_sim
               END DO

               IF (add_highz) CALL add_highz_BAHAMAS(a, pow_ka, nk, na, cosm)

               IF (imode == 65) THEN
                  IF (name == 'AGN_TUNED_nu0' .AND. ixx_names(1) == 'gal_z0.1-0.9' .AND. ixx_names(2) == 'gal_z0.1-0.9') THEN
                     base = 'data/Cl_kk_BAHAMAS_contribution_ell_'
                     ext = '.dat'
                     CALL Cl_contribution(ix, k, a, pow_ka, nk, na, cosm, base, ext)
                  ELSE IF (name == 'AGN_TUNED_nu0' .AND. ixx_names(1) == 'gal_z0.1-0.9' .AND. ixx_names(2) == 'y') THEN
                     base = 'data/Cl_ky_BAHAMAS_contribution_ell_'
                     ext = '.dat'
                     CALL Cl_contribution(ix, k, a, pow_ka, nk, na, cosm, base, ext)
                  ELSE IF (name == 'AGN_TUNED_nu0' .AND. ixx_names(1) == 'y' .AND. ixx_names(2) == 'y') THEN
                     base = 'data/Cl_yy_BAHAMAS_contribution_ell_'
                     ext = '.dat'
                     CALL Cl_contribution(ix, k, a, pow_ka, nk, na, cosm, base, ext)
                  ELSE
                     DEALLOCATE (pow_ka)
                     CYCLE
                  END IF
               END IF

               ! Do the cross correlation
               IF (bin_theory) THEN
                  CALL xpow_pka(ix, all_ell, all_Cl(:, 1, 2), n_all, k, a, pow_ka, nk, na, cosm)
                  CALL bin_theory_ell(all_ell, all_Cl(:, 1, 2), n_all, ell_edge, Cl(:, 1, 2), nl)
               ELSE
                  CALL xpow_pka(ix, ell, Cl(:, 1, 2), nl, k, a, pow_ka, nk, na, cosm)
               END IF

               ! Write data
               outfile = TRIM(outbase)//'_'//TRIM(ixx_names(1))//'-'//TRIM(ixx_names(2))//'.dat'
               OPEN (7, file=outfile)
               DO i = 1, nl
                  WRITE (7, *) ell(i), Cl(i, 1, 2), ell(i)*(ell(i)+1)*Cl(i, 1, 2)/twopi
               END DO
               CLOSE (7)

               ! Deallocate array
               DEALLOCATE (pow_ka)

            END DO
         END DO

      END DO

   END SUBROUTINE Cl_direct_integration

   SUBROUTINE check_fields_symmetry(icosmo, ihm)

      ! Check to see power is the same
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: icosmo
      INTEGER, INTENT(INOUT) :: ihm
      REAL, ALLOCATABLE :: a(:), k(:)
      INTEGER, ALLOCATABLE :: fields(:)
      INTEGER :: i
      REAL, ALLOCATABLE :: pow_li(:, :), pow_2h(:, :, :, :), pow_1h(:, :, :, :), pow_hm(:, :, :, :)
      TYPE(cosmology) :: cosm
      !TYPE(halomod) :: hmod

      REAL, PARAMETER :: kmin = 1e-3
      REAL, PARAMETER :: kmax = 1e2
      INTEGER, PARAMETER :: nk = 128
      REAL, PARAMETER :: amin = 0.1
      REAL, PARAMETER :: amax = 1.0
      INTEGER, PARAMETER :: na = 4
      INTEGER, PARAMETER :: nf = 2
      LOGICAL, PARAMETER :: verbose = .TRUE.

      ! Assigns the cosmological model
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! Set number of k points and k range (log spaced)
      CALL fill_array(log(kmin), log(kmax), k, nk)
      k = exp(k)

      ! Set the scale factor and range (linearly spaced)
      CALL fill_array(amin, amax, a, na)

      ALLOCATE (fields(nf))
      !fields=field_matter
      !fields=field_electron_pressure
      fields = field_gas

      ! Assign the halo model
      !CALL assign_halomod(ihm, hmod, verbose)
      !CALL calculate_HMx(fields, nf, k, nk, a, na, pows_li, pows_2h, pows_1h, pows_hm, hmod, cosm, verbose)
      CALL calculate_HMx_full(fields, k, a, pow_li, pow_2h, pow_1h, pow_hm, nf, nk, na, cosm, ihm)

      ! Write to screen to check they are the same
      WRITE (*, *) 'HMx_DRIVER: All these columns should be idential'
      DO i = 1, nk
         WRITE (*, *) pow_hm(1, 1, i, na), pow_hm(1, 2, i, na), pow_hm(2, 1, i, na), pow_hm(2, 2, i, na)
      END DO
      !IF (pow_hm(1, 1, : , :) .NE. pow_hm(1, 2, :, :)) STOP 'HMX_DRIVER: You have a problem'
      !IF (pow_hm(1, 1, : , :) .NE. pow_hm(2, 1, :, :)) STOP 'HMX_DRIVER: You have a problem'
      !IF (pow_hm(1, 1, : , :) .NE. pow_hm(2, 2, :, :)) STOP 'HMX_DRIVER: You have a problem'
      WRITE (*, *) 'HMx_DRIVER: Done'
      WRITE (*, *)

   END SUBROUTINE check_fields_symmetry

   SUBROUTINE Tinker2010_Fig1(icosmo, ihm)

      ! Produce data for Fig. 1 of Tinker et al. (2010)
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: icosmo
      INTEGER, INTENT(INOUT) :: ihm
      REAL :: nu
      INTEGER :: i
      TYPE(cosmology) :: cosm
      TYPE(halomod) :: hmod

      REAL, PARAMETER :: numin = 0.1
      REAL, PARAMETER :: numax = 10.
      INTEGER, PARAMETER :: n = 256
      CHARACTER(len=256), PARAMETER :: outfile = 'data/Tinker_bias.dat'
      LOGICAL, PARAMETER :: verbose = .TRUE.

      ! Assigns the cosmological model
      icosmo = 1
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ihm = 44
      CALL assign_halomod(ihm, hmod, verbose)
      CALL init_halomod(1., hmod, cosm, verbose)

      OPEN (7, file=outfile)
      DO i = 1, n
         nu = exp(progression(log(numin), log(numax), i, n))
         WRITE (7, *) nu, b_nu(nu, hmod)
      END DO
      CLOSE (7)

   END SUBROUTINE Tinker2010_Fig1

   SUBROUTINE Limber_CCL_comparison(icosmo)

      ! Make data for Limber comparison with CCL
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: icosmo
      REAL, ALLOCATABLE :: k(:), a(:), pow_ka(:, :)
      REAL, ALLOCATABLE :: ell(:), Cl_bm(:)
      INTEGER :: i, j
      INTEGER :: nl, ix(2)
      CHARACTER(len=256) :: outfile
      TYPE(cosmology) :: cosm

      LOGICAL, PARAMETER :: verbose = .TRUE.
      REAL, PARAMETER :: kmin = 1e-3 ! Minimum k to calculate P(k); it will be extrapolated below this by Limber
      REAL, PARAMETER :: kmax = 1e2  ! Maximum k to calculate P(k); it will be extrapolated above this by Limber
      INTEGER, PARAMETER :: nk = 256 ! Number of log-spaced k values (used to be 32)
      REAL, PARAMETER :: amin = 0.1  ! Minimum scale factor (problems with one-halo term if amin is less than 0.1)
      REAL, PARAMETER :: amax = 1.0  ! Maximum scale factor
      INTEGER, PARAMETER :: na = 16  ! Number of linearly-spaced scale factores
      REAL, PARAMETER :: lmin = 2.
      REAL, PARAMETER :: lmax = 10000.

      ! Assigns the cosmological model
      icosmo = 1 ! 1 - Boring
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! k range and array
      CALL fill_array(kmin, kmax, k, nk)

      ! a range and array
      CALL fill_array(amin, amax, a, na)

      ! Allocate and fill array with linear power; P(k,a)
      ALLOCATE (pow_ka(nk, na))
      DO j = 1, na
         DO i = 1, nk
            pow_ka(i, j) = plin(k(i), a(j), flag_power_total, cosm)
         END DO
      END DO

      ! l range
      nl = nint(lmax)-nint(lmin)+1 ! Get all ell
      CALL fill_array(lmin, lmax, ell, nl)
      ALLOCATE (Cl_bm(nl))

      ! Set lensing tracer and call C(l) routine
      !ix = tracer_CFHTLenS_vanWaerbeke2013
      ix = tracer_CFHTLenS_Kilbinger2013
      CALL xpow_pka(ix, ell, Cl_bm, nl, k, a, pow_ka, nk, na, cosm)

      ! Write data
      outfile = 'data/HMx_Cl_test.dat'
      OPEN (7, file=outfile)
      DO i = 1, nl
         WRITE (7, *) ell(i), Cl_bm(i)
      END DO
      CLOSE (7)

   END SUBROUTINE Limber_CCL_comparison

   SUBROUTINE non_linear_halo_bias_model(icosmo, ihm)

      ! Non-linear halo bias model
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: icosmo
      INTEGER, INTENT(INOUT) :: ihm
      REAL, ALLOCATABLE :: k(:), pow_sim(:)
      INTEGER :: field(1)
      INTEGER :: i, j
      INTEGER :: nk
      REAL, ALLOCATABLE :: pow_li(:), pow_2h(:), pow_1h(:), pow_hm(:)
      CHARACTER(len=256) :: outfile
      TYPE(cosmology) :: cosm
      TYPE(halomod) :: hmod

      REAL, PARAMETER :: z = 0.0 ! Set the redshift
      LOGICAL, PARAMETER :: verbose = .TRUE.

      ! Assigns the cosmological model
      icosmo = 200 ! 200 - Frankenemu M000
      !icosmo=5 ! 5 - WMAP5 (Multidark)
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! Get the emulator power
      CALL get_FrankenEmu_power_z(k, pow_sim, nk, z, cosm, rebin=.FALSE.)
      outfile = 'data/power_emu.dat'
      OPEN (7, file=outfile)
      DO i = 1, nk
         WRITE (7, *) k(i), pow_sim(i)
      END DO
      CLOSE (7)

      ! Loop over two different halo-model calculations
      DO j = 1, 2

         ! Set the halo model and output file
         IF (j == 1) THEN
            !ihm=42 ! 42 - Standard halo model with M200c and Tinker
            !ihm=3  ! 3 - Standard halo model with Mv and Sheth-Torman
            ihm = 23 ! 23 - Standard halo model with Mv and Tinker
            outfile = 'data/power.dat'
         ELSE IF (j == 2) THEN
            !ihm=24 ! 24 - Non-linear bias with M200c and Tinker
            !ihm=48 ! 48 - Non-linear bias with Mv and Sheth-Torman
            ihm = 49 ! 49 - Non-linear bias with Mv and Tinker
            outfile = 'data/power_bnl.dat'
         ELSE
            STOP 'HMX_DRIVER: Error, something went wrong in non-linear bias'
         END IF

         ! Initiliasation for the halomodel calcualtion
         CALL assign_halomod(ihm, hmod, verbose)
         CALL init_halomod(scale_factor_z(z), hmod, cosm, verbose)
         CALL print_halomod(hmod, cosm, verbose)

         ! Allocate arrays
         ALLOCATE (pow_li(nk), pow_2h(nk), pow_1h(nk), pow_hm(nk))

         ! Do the halo-model calculation
         field = field_dmonly
         CALL calculate_HMx_a(field, 1, k, nk, pow_li, pow_2h, pow_1h, pow_hm, hmod, cosm, verbose)

         ! Write out the results
         CALL write_power(k, pow_li, pow_2h, pow_1h, pow_hm, nk, outfile, verbose)

         ! Deallocate arrays
         DEALLOCATE (pow_li, pow_2h, pow_1h, pow_hm)

      END DO

   END SUBROUTINE non_linear_halo_bias_model

   SUBROUTINE halo_power_Multidark(imode, icosmo, ihm)

      ! Halo power to compare with Multidark
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: imode
      INTEGER, INTENT(INOUT) :: icosmo
      INTEGER, INTENT(INOUT) :: ihm
      REAL, ALLOCATABLE :: k(:), a(:)
      INTEGER, ALLOCATABLE :: bins(:), snaps(:), fields(:)
      REAL, ALLOCATABLE :: pow_li(:), pow_2h(:, :, :), pow_1h(:, :, :), pow_hm(:, :, :)
      INTEGER :: i, j, j1, j2
      INTEGER :: nk
      CHARACTER(len=256) :: base, ext, mid, outfile
      CHARACTER(len=256), ALLOCATABLE :: bases(:)
      TYPE(cosmology) :: cosm
      TYPE(halomod) :: hmod

      INTEGER, PARAMETER :: ia1 = 1 ! Actual range of a values to do
      INTEGER, PARAMETER :: ia2 = 5 ! Actual range of a values to do
      INTEGER, PARAMETER :: nf = 4  ! Number of fields, either 4 or 9
      INTEGER, PARAMETER :: na = 5 ! Number of redshifts
      LOGICAL, PARAMETER :: verbose = .TRUE.
      CHARACTER(len=256), PARAMETER :: infile = '/Users/Mead/Physics/Multidark/data/power/M512/BDMV_85_bin0_bin0_power.dat'

      IF (imode == 69) THEN
         icosmo = 37 ! 37 - WMAP 5
      ELSE IF (imode == 74) THEN
         icosmo = 43 ! 43 - WMAP 5 (low sigma_8)
      ELSE
         STOP 'HMX_DRIVER: Error, imode not set correctly'
      END IF
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! Set the redshift
      ALLOCATE (snaps(na), a(na))
      snaps(1) = 85
      snaps(2) = 62
      snaps(3) = 58
      snaps(4) = 52
      snaps(5) = 36
      DO i = 1, na
         a(i) = multidark_scale_factor(snaps(i))
         IF (a(i) > 1.) a(i) = 1. ! Prevent negative z
      END DO

      ! Output files
      ALLOCATE (bases(na))
      base = 'data/power_'
      ext = '_hh_'
      DO i = 1, na
         bases(i) = number_file(base, snaps(i), ext)
      END DO

      ! Allocate k range
      nk = file_length(infile)
      ALLOCATE (k(nk))
      OPEN (7, file=infile)
      DO i = 1, nk
         READ (7, *) k(i)
      END DO
      CLOSE (7)

      ! Set the fields
      ALLOCATE (fields(nf), bins(nf))
      IF (nf == 9) THEN
         fields(1) = field_dmonly
         fields(2) = field_halo_11p0_11p5
         fields(3) = field_halo_11p5_12p0
         fields(4) = field_halo_12p0_12p5
         fields(5) = field_halo_12p5_13p0
         fields(6) = field_halo_13p0_13p5
         fields(7) = field_halo_13p5_14p0
         fields(8) = field_halo_14p0_14p5
         fields(9) = field_halo_14p5_15p0
         bins(1) = 0
         bins(2) = 1
         bins(3) = 2
         bins(4) = 3
         bins(5) = 4
         bins(6) = 5
         bins(7) = 6
         bins(8) = 7
         bins(9) = 8
      ELSE IF (nf == 4) THEN
         fields(1) = field_dmonly
         fields(2) = field_halo_12p5_13p0
         fields(3) = field_halo_13p0_13p5
         fields(4) = field_halo_13p5_14p0
         bins(1) = 0
         bins(2) = 4
         bins(3) = 5
         bins(4) = 6
      ELSE
         STOP 'HMX_DRIVER: Error, number of fields specified incorrectly'
      END IF

      ! Allocate arrays
      ALLOCATE (pow_li(nk), pow_2h(nf, nf, nk), pow_1h(nf, nf, nk), pow_hm(nf, nf, nk))

      DO i = ia1, ia2

         !base='data/power_hh_'
         base = bases(i)
         mid = ''
         DO j = 1, 2

            IF (j == 1) THEN
               ihm = 23 ! 23 - Standard with Tinker
               ext = '_standard.dat'
            ELSE IF (j == 2) THEN
               !ihm=48 ! 48 - Non-linear halo bias with Sheth & Torman
               ihm = 49 ! 49 - Non-linear halo bias with Tinker
               ext = '_bnl.dat'
            ELSE
               STOP 'HMx_DRIVER: Error, something went wrong'
            END IF

            ! Initiliasation for the halomodel calcualtion
            CALL assign_halomod(ihm, hmod, verbose)
            hmod%mmin = 1e10 ! 1e7 is the default, be really careful here
            hmod%mmax = 1e16 ! 1e17 is the default, be really careful here
            hmod%n = 1024    ! 512 is at least necessary due to thin halo bins (lots of zero points)
            CALL init_halomod(a(i), hmod, cosm, verbose)
            CALL print_halomod(hmod, cosm, verbose)

            ! Do the halo-model calculation
            CALL calculate_HMx_a(fields, nf, k, nk, pow_li, pow_2h, pow_1h, pow_hm, hmod, cosm, verbose)

            ! Write out the results
            DO j1 = 1, nf
               DO j2 = 1, nf
                  outfile = number_file2(base, bins(j1), mid, bins(j2), ext)
                  WRITE (*, *) 'HMx_DRIVER: Output: ', trim(outfile)
                  CALL write_power(k, pow_li, pow_2h(j1, j2, :), pow_1h(j1, j2, :), pow_hm(j1, j2, :), nk, outfile, verbose)
               END DO
            END DO
            WRITE (*, *)

         END DO

      END DO

   END SUBROUTINE halo_power_Multidark

END PROGRAM HMx_driver
