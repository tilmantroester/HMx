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
   USE logical_operations
   USE owls
   USE owls_extras
   USE interpolate
   USE cosmology_functions
   USE array_operations

   ! TODO: Many pow(1,1,nk) could be pow(nk)
   ! TODO: Too many different pow(:,:), pow(:,:,:), pow(:,:,:,:), ...
   ! TODO: Maybe could use pow1(:), pow2(:,:), ...

   IMPLICIT NONE

   ! Parameter definitions
   REAL, ALLOCATABLE :: k(:), a(:)
   REAL, ALLOCATABLE :: k_sim(:), pow_sim(:), pow_ql(:), pow_oh(:), pow_hf(:)
   REAL, ALLOCATABLE :: pow_li(:), pow_2h(:, :, :), pow_1h(:, :, :), pow_hm(:, :, :)
   REAL, ALLOCATABLE :: pow_ka(:, :)
   REAL, ALLOCATABLE :: powd_li(:), powd_2h(:), powd_1h(:), powd_hm(:)
   REAL, ALLOCATABLE :: pows_li(:, :), pows_2h(:, :, :, :), pows_1h(:, :, :, :), pows_hm(:, :, :, :)
   REAL, ALLOCATABLE :: powb_hm(:, :, :, :)
   REAL, ALLOCATABLE :: ell(:), ell_edge(:), Cl(:, :, :), Cl_bm(:), theta(:), xi(:, :)
   REAL, ALLOCATABLE :: all_ell(:), all_Cl(:, :, :)
   REAL, ALLOCATABLE :: zs(:), masses(:)
   REAL, ALLOCATABLE :: z_tab(:), HI_frac(:), nu_tab(:)
   INTEGER :: i, j, ii, jj, nk, na, j1, j2, nf, nt, nx, n_all, l, n_edge, ik
   INTEGER :: n, nl, nz, nth, nnz, m, ipa, ncos, ncore, nfeed, nsim, nnu
   INTEGER :: ip(2), ix(2), field(1)
   INTEGER, ALLOCATABLE :: fields(:), ixx(:)
   REAL :: kmin, kmax, amin, amax, lmin, lmax, thmin, thmax, zmin, zmax
   REAL :: rcore_min, rcore_max, lmax_xi
   REAL :: z, z1, z2, r1, r2, a1, a2, nu1, nu2, B_NL, I_NL
   TYPE(cosmology) :: cosm
   TYPE(cosmology), ALLOCATABLE :: cosms(:)
   TYPE(halomod) :: hmod
   TYPE(halomod), ALLOCATABLE :: hmods(:)
   TYPE(projection) :: proj(2), pro
   CHARACTER(len=256) :: infile, outfile, base, mid, ext, dir, name, fname, inbase, outbase, inext, outext
   CHARACTER(len=256) :: mode, halomodel, cosmo
   CHARACTER(len=256), ALLOCATABLE :: ixx_names(:), bases(:)
   INTEGER :: imode, icosmo, iowl, ihm, irho, itest, jtest, nhm
   INTEGER :: imin, imax, mesh
   REAL :: sig8min, sig8max
   REAL :: mass, m1, m2, nu, numin, numax, mf
   REAL :: c, rmin, rmax, rv, rs, p1, p2, cmin, cmax, cbar
   REAL :: spam, mmin_bnl, mmax_bnl
   CHARACTER(len=1) :: crap
   LOGICAL :: verbose2
   INTEGER :: ia1, ia2
   INTEGER, ALLOCATABLE :: bins(:)

   ! Baryon stuff
   REAL :: param_min, param_max, param, param_neat
   LOGICAL :: ilog

   ! Tests
   REAL :: error, error_max

   ! Halo-model Parameters
   LOGICAL, PARAMETER :: verbose = .TRUE. ! Verbosity
   REAL, PARAMETER :: mmin = mmin_HMx     ! Minimum halo mass for the calculation
   REAL, PARAMETER :: mmax = mmax_HMx     ! Maximum halo mass for the calculation

   ! Test parameters
   INTEGER, PARAMETER :: iseed = 1
   REAL, PARAMETER :: tolerance = 3e-3
   LOGICAL :: verbose_tests = .FALSE.
   LOGICAl :: ifail = .FALSE.

   ! Benchmark parameters
   LOGICAL, PARAMETER :: Alonso_k = .TRUE.

   ! Output choices
   LOGICAL, PARAMETER :: icumulative = .TRUE. ! Do cumlative distributions for breakdown
   LOGICAL, PARAMETER :: ifull = .FALSE.      ! Do only full halo model C(l), xi(theta) calculations (quicker, no breakdown ...)

   ! Triad direct comparison
   LOGICAL, PARAMETER :: bin_theory = .FALSE.     ! Should the theory be binned in the same way as measurements?
   LOGICAL, PARAMETER :: cut_nyquist = .TRUE.     ! Should the BAHAMAS measured P(k) be cut above the Nyquist frequency?
   LOGICAL, PARAMETER :: subtract_shot = .TRUE.   ! Should the BAHAMAS measured P(k) have shot-noise subtracted?
   LOGICAL, PARAMETER :: response_triad = .FALSE. ! Should I treat the BAHAMAS P(k) as HMcode response?
   LOGICAL, PARAMETER :: add_highz = .FALSE.      ! Add in z=3 power
   INTEGER, PARAMETER :: mesh_BAHAMAS = 1024      ! Mesh size for BAHAMAS P(k) measurements
   INTEGER, PARAMETER :: nz_BAHAMAS = 15          ! Number of BAHAMAS redshift slices to use (4, 11 or 15)
   REAL, PARAMETER :: kmax_BAHAMAS = 1000.        ! Maximum k to trust the BAHAMAS P(k) [h/Mpc]

   ! Cross correlation
   REAL, PARAMETER :: kmin_xpow = 1e-3 ! Minimum k to calculate P(k); it will be extrapolated below this by Limber
   REAL, PARAMETER :: kmax_xpow = 1e2  ! Maximum k to calculate P(k); it will be extrapolated above this by Limber
   INTEGER, PARAMETER :: nk_xpow = 256 ! Number of log-spaced k values (used to be 32)
   REAL, PARAMETER :: amin_xpow = 0.1  ! Minimum scale factor (problems with one-halo term if amin is less than 0.1)
   REAL, PARAMETER :: amax_xpow = 1.0  ! Maximum scale factor
   INTEGER, PARAMETER :: na_xpow = 16  ! Number of linearly-spaced scale factores

   CALL get_command_argument(1, mode)
   IF (mode == '') THEN
      imode = -1
   ELSE
      READ (mode, *) imode
   END IF

   CALL get_command_argument(2, cosmo)
   IF (cosmo == '') THEN
      icosmo = -1
   ELSE
      READ (cosmo, *) icosmo
   END IF

   CALL get_command_argument(3, halomodel)
   IF (halomodel == '') THEN
      ihm = -1
   ELSE
      READ (halomodel, *) ihm
   END IF

   ! Initial white space
   WRITE (*, *)

   ! Choose mode
   IF (imode == -1) THEN
      WRITE (*, *) 'HMx_DRIVER: Choose what to do'
      WRITE (*, *) '============================='
      WRITE (*, *) ' 0 - Gravity-only power spectrum at z=0'
      WRITE (*, *) ' 1 - 3D Matter power spectrum over multiple z'
      WRITE (*, *) ' 2 - Hydrodynamical halo model'
      WRITE (*, *) ' 3 - Run diagnostics for haloes'
      WRITE (*, *) ' 4 - Do random baryon parameters for bug testing'
      WRITE (*, *) ' 5 - Lensing diagnostics'
      WRITE (*, *) ' 6 - n(z) check'
      WRITE (*, *) ' 7 - Do general angular cross correlation'
      WRITE (*, *) ' 8 - Angular cross correlation as a function of cosmology'
      WRITE (*, *) ' 9 - Breakdown angular correlations in halo mass'
      WRITE (*, *) '10 - Breakdown angular correlations in redshift'
      WRITE (*, *) '11 - Do general angular cross correlation and correlation functions'
      WRITE (*, *) '12 - Triad'
      WRITE (*, *) '13 - Cross-correlation coefficient'
      WRITE (*, *) '14 - 3D spectra for variations in baryon parameters'
      WRITE (*, *) '15 - 3D spectra for cosmo-OWLS models (Feedback models and k values) NOT SUPPORTED'
      WRITE (*, *) '16 - 3D spectra for BAHAMAS models (AGN models and k values)'
      WRITE (*, *) '17 - 3D spectra for user choice of fields'
      WRITE (*, *) '18 - 3D bias'
      WRITE (*, *) '19 - Make CCL benchmark data'
      WRITE (*, *) '20 - Make data for Ma et al. (2015) Fig. 1'
      WRITE (*, *) '21 - W(k) integrand diagnostics'
      WRITE (*, *) '22 - Time W(k) integration methods'
      WRITE (*, *) '23 - Produce DE response results from Mead (2017)'
      WRITE (*, *) '24 - TESTS: Projection'
      WRITE (*, *) '25 - Halo-void model'
      WRITE (*, *) '26 - TESTS: DMONLY spectra; HMcode'
      WRITE (*, *) '27 - Comparison with Mira Titan nodes'
      WRITE (*, *) '28 - Comparison with FrankenEmu nodes'
      WRITE (*, *) '29 - Comparison with random Mira Titan cosmology'
      WRITE (*, *) '30 - Comparison with random FrankenEmu cosmology'
      WRITE (*, *) '31 - PAPER: Breakdown 3D hydro power in halo mass'
      WRITE (*, *) '32 - PAPER: Hydro power for baseline model'
      WRITE (*, *) '33 - PAPER: Effect of parameter variations on baseline model hydro power'
      WRITE (*, *) '34 - Write HMx hydro parameter variations with T_AGN and z'
      WRITE (*, *) '35 - Power spectra of cored halo profiles'
      WRITE (*, *) '36 - TESTS: hydro spectra'
      WRITE (*, *) '37 - Produce CFHTLenS correlation functions'
      WRITE (*, *) '38 - Tilman AGN model triad for all feedback models'
      WRITE (*, *) '39 - Tilman AGN model Triad for all feedback models (Triad 3 ell)'
      WRITE (*, *) '40 - Halo bias'
      WRITE (*, *) '41 - 3D power as a function of cosmology'
      WRITE (*, *) '42 - PAPER: Contributions to k-k C(l) integral'
      WRITE (*, *) '43 - PAPER: Contributions to k-y C(l) integral'
      WRITE (*, *) '44 - Triad with Tilman model'
      WRITE (*, *) '45 - Comparison of Sheth-Tormen vs. Tinker mass function'
      WRITE (*, *) '46 - Mass function and bias plots'
      WRITE (*, *) '47 - Make CMB lensing to compare with CAMB'
      WRITE (*, *) '48 - HI bias'
      WRITE (*, *) '49 - HI mass fractions'
      WRITE (*, *) '50 - Mass function changes with Lbox'
      WRITE (*, *) '51 - Compare power with and without scatter'
      WRITE (*, *) '52 - Hydrodynamical halo model with BAHAMAS k range'
      WRITE (*, *) '53 - Breakdown 3D hydro power in halo mass II'
      WRITE (*, *) '54 - Trispectrum test'
      WRITE (*, *) '55 - Triad for all feedback models'
      WRITE (*, *) '56 - PAPER: Triad for all feedback models (Triad 3 ell)'
      WRITE (*, *) '57 - Direct integration of measured 3D BAHAMAS spectra: Triad 3'
      WRITE (*, *) '58 - Multiple same fields check'
      WRITE (*, *) '59 - Tinker (2010) bias plot check'
      WRITE (*, *) '60 - Make data for Limber comparison with CCL'
      WRITE (*, *) '61 - Direct integration of measured 3D BAHAMAS spectra: Triad 4'
      WRITE (*, *) '62 - Direct integration of measured 3D BAHAMAS spectra: Triad 5'
      WRITE (*, *) '63 - Contributions to y-y C(l) integral'
      WRITE (*, *) '64 - CIB 545 GHz'
      WRITE (*, *) '65 - Contributions to Limber integrals from BAHAMAS'
      WRITE (*, *) '66 - Direct integration of measured 3D BAHAMAS spectra: Triad 5 but larger ell range'
      WRITE (*, *) '67 - Compute matter power using non-linear halo bias'
      WRITE (*, *) '68 - Compute non-linear halo bias integrand'
      WRITE (*, *) '69 - Matter, halo power spectra with non-linear bias'
      WRITE (*, *) '70 - Comparison with Cosmic Emu nodes'
      WRITE (*, *) '71 - Comparison with random Cosmic Emu cosmology'
      WRITE (*, *) '72 - Run halo model over lots of random cosmological parameters'
      READ (*, *) imode
      WRITE (*, *) '============================'
      WRITE (*, *)
   END IF

   IF (imode == 0) THEN

      ! 0 - Calculate halo model at a single z

      ! Set number of k points and k range (log spaced)
      nk = 128
      kmin = 1e-3
      kmax = 1e2
      CALL fill_array(log(kmin), log(kmax), k, nk)
      k = exp(k)

      ! Assigns the cosmological model
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! Sets the redshift
      z = 1.0

      ! Initiliasation for the halomodel calcualtion
      CALL assign_halomod(ihm, hmod, verbose)
      CALL init_halomod(mmin, mmax, scale_factor_z(z), hmod, cosm, verbose)
      CALL print_halomod(hmod, cosm, verbose)

      ! Allocate arrays
      ALLOCATE (powd_li(nk), powd_2h(nk), powd_1h(nk), powd_hm(nk))

      ! Do the halo-model calculation
      field = field_dmonly
      CALL calculate_HMx_a(field, 1, k, nk, powd_li, powd_2h, powd_1h, powd_hm, hmod, cosm, verbose, response=.FALSE.)

      ! Write out the results
      outfile = 'data/power.dat'
      CALL write_power(k, powd_li, powd_2h, powd_1h, powd_hm, nk, outfile, verbose)

   ELSE IF (imode == 1) THEN

      ! Assigns the cosmological model
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! Assign the halo model
      CALL assign_halomod(ihm, hmod, verbose)

      ! Set number of k points and k range (log spaced)
      nk = 128
      kmin = 1e-3
      kmax = 1e2
      CALL fill_array(log(kmin), log(kmax), k, nk)
      k = exp(k)

      ! Set the scale factor and range (linearly spaced)
      !na=16
      !amin=0.2
      !amax=1.0
      !CALL fill_array(amin,amax,a,na)

      ! Set the number of redshifts and range (linearly spaced) and convert z -> a
      zmin = 0.
      zmax = 4.
      na = 16
      CALL fill_array(zmin, zmax, a, na)
      DO i = 1, na
         a(i) = scale_factor_z(a(i)) ! Note that this is correct because 'a' here is actually 'z'
      END DO

      field = field_dmonly
      CALL calculate_HMx(field, 1, mmin, mmax, k, nk, a, na, pows_li, pows_2h, pows_1h, pows_hm, hmod, cosm, verbose, response=.FALSE.)

      base = 'data/power'
      CALL write_power_a_multiple(k, a, pows_li, pows_2h(1, 1, :, :), pows_1h(1, 1, :, :), pows_hm(1, 1, :, :), nk, na, base, verbose)

   ELSE IF (imode == 2 .OR. imode == 15 .OR. imode == 16 .OR. imode == 32 .OR. imode == 52) THEN

      ! Make cross power spectra of all different components of haloes as well as pressure
      !  2 - Generic hydro
      ! 15 - cosmo-OWLS
      ! 16 - BAHAMAS
      ! 32 - PAPER: Baseline hydro

      IF (imode == 2 .OR. imode == 52) THEN

         ! Generic hydro

         ! Only do one 'model' here
         n = 1

         ! Set the redshift
         nz = 4
         ALLOCATE (z_tab(nz))
         z_tab(1) = 0.0
         z_tab(2) = 0.5
         z_tab(3) = 1.0
         z_tab(4) = 2.0

         ! Set number of k points and k range (log spaced)
         IF (imode == 2) THEN

            nk = 128
            kmin = 1e-3
            kmax = 1e2
            CALL fill_array(log(kmin), log(kmax), k, nk)
            k = exp(k)

         ELSE IF (imode == 52) THEN

            ! Get the k values from the simulation measured P(k)
            infile = '/Users/Mead/Physics/BAHAMAS/power/M1536/DMONLY_nu0_L400N1024_WMAP9_snap32_all_all_power.dat'
            CALL read_k_values(infile, k, nk)

         ELSE

            STOP 'HMx_DRIVER: imode specified incorrectly'

         END IF

      ELSE IF (imode == 32) THEN

         ! Generic hydro

         ! Only do one 'model' here
         n = 1

         ! Set the redshift
         nz = 4
         ALLOCATE (z_tab(nz))
         z_tab(1) = 0.0
         z_tab(2) = 0.5
         z_tab(3) = 1.0
         z_tab(4) = 2.0

         ! Set number of k points and k range (log spaced)
         nk = 128
         kmin = 1e-3
         kmax = 1e2
         CALL fill_array(log(kmin), log(kmax), k, nk)
         k = exp(k)

      ELSE IF (imode == 15) THEN

         ! cosmo-OWLS
         STOP 'HMx_DRIVER: not tested in ages, be very careful'

         ! Do from REF, NOCOOL, AGN, AGN 8.5, AGN 8.7
         n = 5

         ! Set the redshift
         nz = 1
         ALLOCATE (z_tab(nz))
         z_tab(1) = 0.

         ! Get the k values from the simulation measured P(k)
         infile = '/Users/Mead/Physics/cosmo-OWLS/power/N800/DMONLY_all_all_power.dat'
         CALL read_k_values(infile, k, nk)

      ELSE IF (imode == 16) THEN

         ! BAHAMAS

         ! Do AGN, AGN-lo and AGN-hi
         n = 3

         ! Set the redshift
         nz = 4
         ALLOCATE (z_tab(nz))
         z_tab(1) = 0.0
         z_tab(2) = 0.5
         z_tab(3) = 1.0
         z_tab(4) = 2.0

         ! Get the k values from the simulation measured P(k)
         infile = '/Users/Mead/Physics/BAHAMAS/power/M1536/DMONLY_nu0_L400N1024_WMAP9_snap32_all_all_power.dat'
         CALL read_k_values(infile, k, nk)

      END IF

      ! Field types
      nf = 5
      ALLOCATE (fields(nf))
      fields(1) = field_matter
      fields(2) = field_cdm
      fields(3) = field_gas
      fields(4) = field_star
      fields(5) = field_electron_pressure

      ! Allocate the arrays for P(k)
      ALLOCATE (powd_li(nk), powd_2h(nk), powd_1h(nk), powd_hm(nk))
      ALLOCATE (pow_li(nk), pow_2h(nf, nf, nk), pow_1h(nf, nf, nk), pow_hm(nf, nf, nk))

      DO iowl = 1, n

         DO j = 1, nz

            z = z_tab(j)

            !Assigns the cosmological model
            IF (imode == 2 .OR. imode == 52) THEN
               icosmo = 4
            ELSE IF (imode == 15) THEN
               icosmo = 2
            ELSE IF (imode == 16) THEN
               icosmo = 4
            ELSE IF (imode == 32) THEN
               icosmo = 4
               ihm = 3
            END IF

            CALL assign_cosmology(icosmo, cosm, verbose)
            CALL init_cosmology(cosm)
            CALL print_cosmology(cosm)
            CALL assign_halomod(ihm, hmod, verbose)

            IF (imode == 15) THEN

               ! cosmo-OWLS
               IF (imode == 15 .AND. iowl == 1) THEN
                  name = 'REF'
                  fname = name
                  ! From my fitting by eye
                  hmod%alpha = 2.
                  hmod%eps = 1.
                  hmod%Gamma = 1.24
                  hmod%M0 = 1e13
                  hmod%Astar = 0.055
               ELSE IF (imode == 15 .AND. iowl == 2) THEN
                  name = 'NOCOOL'
                  fname = name
                  ! From my fitting by eye
                  hmod%alpha = 2.
                  hmod%eps = 1.
                  hmod%Gamma = 1.1
                  hmod%M0 = 0.
                  hmod%Astar = 0.
               ELSE IF (imode == 15 .AND. iowl == 3) THEN
                  name = 'AGN'
                  fname = name
                  ! From Tilman's preliminary results
                  hmod%alpha = 0.52
                  hmod%eps = 1.
                  hmod%Gamma = 1.17
                  hmod%M0 = 1.047e14
                  hmod%Astar = 0.02
               ELSE IF (imode == 15 .AND. iowl == 4) THEN
                  name = 'AGN 8.5'
                  fname = 'AGN8p5'
                  ! From Tilman's preliminary results
                  hmod%alpha = 0.56
                  hmod%eps = 1.
                  hmod%Gamma = 1.19
                  hmod%M0 = 3.548e14
                  hmod%Astar = 0.01
               ELSE IF (imode == 15 .AND. iowl == 5) THEN
                  name = 'AGN 8.7'
                  fname = 'AGN8p7'
                  ! From Tilman's preliminary results
                  hmod%alpha = 0.53
                  hmod%eps = 1.
                  hmod%Gamma = 1.21
                  hmod%M0 = 7.586e14
                  hmod%Astar = 0.01
               END IF

            END IF

            !BAHAMAS
            IF (imode == 16) THEN

               IF (iowl == 1) THEN

                  ! Simulation name and file name
                  name = 'AGN'
                  fname = 'AGN'

                  ! Best z=0 fit on 21/06/2018
                  IF (ihm == 4) THEN

                     IF (z == 0.) THEN
                        hmod%alpha = 0.379
                        hmod%eps = 10**(-0.061)
                        hmod%Gamma = 1.205
                        hmod%M0 = 10**(13.823)
                        hmod%Astar = 0.029
                        hmod%Twhim = 10**(5.754)
                     ELSE IF (z == 0.5) THEN
                        hmod%alpha = 0.537
                        hmod%eps = 10**(-0.209)
                        hmod%Gamma = 1.163
                        hmod%M0 = 10**(13.964)
                        hmod%Astar = 0.024
                        hmod%Twhim = 10**(5.750)
                     ELSE IF (z == 1. .OR. z == 2.) THEN
                        hmod%alpha = 0.755
                        hmod%eps = 10**(-0.985)
                        hmod%Gamma = 1.162
                        hmod%M0 = 10**(13.673)
                        hmod%Astar = 0.019
                        hmod%Twhim = 10**(5.057)
                     END IF

                  ELSE IF (ihm == 6) THEN

                     IF (z == 0.) THEN
                        hmod%alpha = 0.195
                        hmod%eps = 10**(-0.484)
                        hmod%Gamma = 1.399
                        hmod%M0 = 10**(13.807)
                        hmod%Astar = 0.020
                        hmod%Twhim = 10**(6.049)
                     ELSE IF (z == 0.5) THEN
                        hmod%alpha = 0.619
                        hmod%eps = 10**(-0.309)
                        hmod%Gamma = 1.507
                        hmod%M0 = 10**(14.937)
                        hmod%Astar = 0.021
                        hmod%Twhim = 10**(5.987)
                     ELSE IF (z == 1. .OR. z == 2.) THEN
                        hmod%alpha = 0.384
                        hmod%eps = 10**(-0.475)
                        hmod%Gamma = 1.183
                        hmod%M0 = 10**(14.562)
                        hmod%Astar = 0.017
                        hmod%Twhim = 10**(5.796)
                     END IF

                  ELSE IF (ihm == 3 .OR. ihm == 14) THEN

                     IF (z == 0.) THEN
                        hmod%alpha = 0.428
                        hmod%eps = 10**(0.015)
                        hmod%Gamma = 1.287
                        hmod%M0 = 10**(13.233)
                        hmod%Astar = 0.030
                        hmod%Twhim = 10**(5.404)
                     ELSE IF (z == 0.5) THEN
                        hmod%alpha = 0.742
                        hmod%eps = 10**(0.148)
                        hmod%Gamma = 1.516
                        hmod%M0 = 10**(12.688)
                        hmod%Astar = 0.026
                        hmod%Twhim = 10**(5.531)
                     ELSE IF (z == 1. .OR. z == 2.) THEN
                        hmod%alpha = 0.830
                        hmod%eps = 10**(0.169)
                        hmod%Gamma = 1.487
                        hmod%M0 = 10**(12.004)
                        hmod%Astar = 0.024
                        hmod%Twhim = 10**(5.643)
                     END IF

                  ELSE IF (ihm == 17 .OR. ihm == 18 .OR. ihm == 19) THEN

                     hmod%Theat = 10**7.8

                  END IF

               ELSE IF (iowl == 3) THEN

                  ! Simulation name and file name
                  name = 'AGN high'
                  fname = 'AGN-hi'
                  hmod%Theat = 10**8.0

                  ! Best z=0 fit on 21/06/2018
                  IF (ihm == 4) THEN

                     IF (z == 0.) THEN
                        hmod%alpha = 0.421
                        hmod%eps = 10**(-0.154)
                        hmod%Gamma = 1.211
                        hmod%M0 = 10**(14.316)
                        hmod%Astar = 0.026
                        hmod%Twhim = 10**(5.801)
                     ELSE IF (z == 0.5) THEN
                        hmod%alpha = 0.599
                        hmod%eps = 10**(-0.686)
                        hmod%Gamma = 1.158
                        hmod%M0 = 10**(14.455)
                        hmod%Astar = 0.022
                        hmod%Twhim = 10**(5.849)
                     ELSE IF (z == 1. .OR. z == 2.) THEN
                        hmod%alpha = 0.659
                        hmod%eps = 10**(-1.993)
                        hmod%Gamma = 1.151
                        hmod%M0 = 10**(14.386)
                        hmod%Astar = 0.018
                        hmod%Twhim = 10**(5.829)
                     END IF

                  ELSE IF (ihm == 6) THEN

                     IF (z == 0.) THEN
                        hmod%alpha = 0.302
                        hmod%eps = 10**(-0.429)
                        hmod%Gamma = 1.339
                        hmod%M0 = 10**(14.662)
                        hmod%Astar = 0.026
                        hmod%Twhim = 10**(6.080)
                     ELSE IF (z == 0.5) THEN
                        hmod%alpha = 0.751
                        hmod%eps = 10**(-0.013)
                        hmod%Gamma = 1.502
                        hmod%M0 = 10**(14.958)
                        hmod%Astar = 0.014
                        hmod%Twhim = 10**(5.959)
                     ELSE IF (z == 1. .OR. z == 2.) THEN
                        hmod%alpha = 0.371
                        hmod%eps = 10**(-0.127)
                        hmod%Gamma = 1.101
                        hmod%M0 = 10**(14.966)
                        hmod%Astar = 0.018
                        hmod%Twhim = 10**(6.040)
                     END IF

                  ELSE IF (ihm == 3 .OR. ihm == 14) THEN

                     IF (z == 0.) THEN
                        hmod%alpha = 0.528
                        hmod%eps = 10**(0.038)
                        hmod%Gamma = 1.505
                        hmod%M0 = 10**(13.638)
                        hmod%Astar = 0.027
                        hmod%Twhim = 10**(5.078)
                     ELSE IF (z == 0.5) THEN
                        hmod%alpha = 0.742
                        hmod%eps = 10**(0.125)
                        hmod%Gamma = 1.547
                        hmod%M0 = 10**(13.481)
                        hmod%Astar = 0.024
                        hmod%Twhim = 10**(5.786)
                     ELSE IF (z == 1. .OR. z == 2.) THEN
                        hmod%alpha = 0.918
                        hmod%eps = 10**(0.289)
                        hmod%Gamma = 1.996
                        hmod%M0 = 10**(13.022)
                        hmod%Astar = 0.022
                        hmod%Twhim = 10**(5.849)
                     END IF

                  ELSE IF (ihm == 17 .OR. ihm == 18 .OR. ihm == 19) THEN

                     hmod%Theat = 10**8.0

                  END IF

               ELSE IF (iowl == 2) THEN

                  ! Simulation name and file name
                  name = 'AGN low'
                  fname = 'AGN-lo'
                  hmod%Theat = 10**7.6

                  ! Best z=0 21/06/2018
                  IF (ihm == 4) THEN

                     IF (z == 0.) THEN
                        hmod%alpha = 0.353
                        hmod%eps = 10**(-0.017)
                        hmod%Gamma = 1.199
                        hmod%M0 = 10**(13.517)
                        hmod%Astar = 0.031
                        hmod%Twhim = 10**(5.767)
                     ELSE IF (z == 0.5) THEN
                        hmod%alpha = 0.491
                        hmod%eps = 10**(-0.104)
                        hmod%Gamma = 1.158
                        hmod%M0 = 10**(13.663)
                        hmod%Astar = 0.025
                        hmod%Twhim = 10**(5.721)
                     ELSE IF (z == 1. .OR. z == 2.) THEN
                        hmod%alpha = 0.532
                        hmod%eps = 10**(-0.990)
                        hmod%Gamma = 1.079
                        hmod%M0 = 10**(13.815)
                        hmod%Astar = 0.021
                        hmod%Twhim = 10**(5.601)
                     END IF

                  ELSE IF (ihm == 6) THEN

                     IF (z == 0.) THEN
                        hmod%alpha = 0.181
                        hmod%eps = 10**(-0.489)
                        hmod%Gamma = 1.432
                        hmod%M0 = 10**(13.632)
                        hmod%Astar = 0.023
                        hmod%Twhim = 10**(6.099)
                     ELSE IF (z == 0.5) THEN
                        hmod%alpha = 0.420
                        hmod%eps = 10**(-0.413)
                        hmod%Gamma = 1.467
                        hmod%M0 = 10**(14.606)
                        hmod%Astar = 0.021
                        hmod%Twhim = 10**(5.874)
                     ELSE IF (z == 1. .OR. z == 2.) THEN
                        hmod%alpha = 0.715
                        hmod%eps = 10**(-0.405)
                        hmod%Gamma = 1.154
                        hmod%M0 = 10**(14.786)
                        hmod%Astar = 0.019
                        hmod%Twhim = 10**(5.434)
                     END IF

                  ELSE IF (ihm == 3 .OR. ihm == 14) THEN

                     IF (z == 0.) THEN
                        hmod%alpha = 0.409
                        hmod%eps = 10**(0.045)
                        hmod%Gamma = 1.275
                        hmod%M0 = 10**(12.737)
                        hmod%Astar = 0.032
                        hmod%Twhim = 10**(5.357)
                     ELSE IF (z == 0.5) THEN
                        hmod%alpha = 0.698
                        hmod%eps = 10**(0.159)
                        hmod%Gamma = 1.393
                        hmod%M0 = 10**(12.012)
                        hmod%Astar = 0.028
                        hmod%Twhim = 10**(5.491)
                     ELSE IF (z == 1. .OR. z == 2.) THEN
                        hmod%alpha = 0.832
                        hmod%eps = 10**(0.312)
                        hmod%Gamma = 1.344
                        hmod%M0 = 10**(12.020)
                        hmod%Astar = 0.025
                        hmod%Twhim = 10**(5.723)
                     END IF

                  ELSE IF (ihm == 17 .OR. ihm == 18 .OR. ihm == 19) THEN

                     hmod%Theat = 10**7.6

                  END IF

               END IF

            END IF

            IF (imode == 15) WRITE (*, *) 'Comparing to OWLS model: ', TRIM(name)
            IF (imode == 16) WRITE (*, *) 'Comparing to BAHAMAS model: ', TRIM(name)

            ! Initiliasation for the halomodel calcualtion after variables changed
            CALL init_halomod(mmin, mmax, scale_factor_z(z), hmod, cosm, verbose)
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
               IF (j == 1) base = 'data/power_z0.0_'
               IF (j == 2) base = 'data/power_z0.5_'
               IF (j == 3) base = 'data/power_z1.0_'
               IF (j == 4) base = 'data/power_z2.0_'
               mid = ''
               ext = '.dat'
            ELSE IF (imode == 15) THEN
               base = 'data/power_'//TRIM(fname)//'_'
               mid = ''
               ext = '.dat'
            ELSE IF (imode == 16) THEN
               IF (j == 1) base = 'data/power_'//TRIM(fname)//'_z0.0_'
               IF (j == 2) base = 'data/power_'//TRIM(fname)//'_z0.5_'
               IF (j == 3) base = 'data/power_'//TRIM(fname)//'_z1.0_'
               IF (j == 4) base = 'data/power_'//TRIM(fname)//'_z2.0_'
               mid = ''
               ext = '.dat'
            END IF

            ! Dark-matter only
            IF (imode == 2 .OR. imode == 32 .OR. imode == 52) THEN
               IF (j == 1) outfile = 'data/power_z0.0.dat'
               IF (j == 2) outfile = 'data/power_z0.5.dat'
               IF (j == 3) outfile = 'data/power_z1.0.dat'
               IF (j == 4) outfile = 'data/power_z2.0.dat'
            ELSE IF (imode == 15) THEN
               outfile = 'data/power_DMONLY_00.dat'
            ELSE IF (imode == 16) THEN
               IF (j == 1) outfile = 'data/power_DMONLY_z0.0_00.dat'
               IF (j == 2) outfile = 'data/power_DMONLY_z0.5_00.dat'
               IF (j == 3) outfile = 'data/power_DMONLY_z1.0_00.dat'
               IF (j == 4) outfile = 'data/power_DMONLY_z2.0_00.dat'
            END IF

            ! Write some things to the screen
            field = field_dmonly
            WRITE (*, *) field(1), field(1), TRIM(outfile)

            ! Do the calculation for DMonly
            CALL calculate_HMx_a(field, 1, k, nk, powd_li, powd_2h, powd_1h, powd_hm, hmod, cosm, verbose, response=.FALSE.)
            CALL write_power(k, powd_li, powd_2h, powd_1h, powd_hm, nk, outfile, verbose)

            ! Do the calculation for the rest of the fields
            CALL calculate_HMx_a(fields, nf, k, nk, pow_li, pow_2h, pow_1h, pow_hm, hmod, cosm, verbose, response=.FALSE.)

            ! Loop over fields and write data
            DO j1 = 1, nf
               DO j2 = j1, nf

                  ! Fix output file and write to screen
                  outfile = number_file2(base, fields(j1), mid, fields(j2), ext)

                  ! Set the halo types and write to screen
                  WRITE (*, *) fields(j1), fields(j2), TRIM(outfile)

                  ! Write P(k
                  CALL write_power(k, pow_li, pow_2h(j1, j2, :), pow_1h(j1, j2, :), pow_hm(j1, j2, :), nk, outfile, verbose=.FALSE.)

               END DO
            END DO

         END DO

      END DO

   ELSE IF (imode == 3) THEN

      ! Halo diagnostics

      ! Assigns the cosmological model
      CALL assign_cosmology(icosmo, cosm, verbose)
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
         CALL assign_halomod(ihm, hmod, verbose)
         CALL init_halomod(mmin, mmax, scale_factor_z(z), hmod, cosm, verbose)
         CALL print_halomod(hmod, cosm, verbose)

         ! Runs the diagnostics
         dir = 'data'
         CALL halo_diagnostics(hmod, cosm, dir)
         CALL halo_definitions(hmod, cosm, dir)
         CALL halo_properties(hmod, dir)

      END DO

   ELSE IF (imode == 4 .OR. imode == 72) THEN

      !  4 - Random baryon parameters
      ! 72 - Random cosmological parameters

      ! Set the random number generator
      CALL RNG_set(iseed)

      ! Set number of k points and k range (log spaced)
      nk = 128
      kmin = 1e-3
      kmax = 1e2
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
         fields(4) = field_star
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

      ! Set redshift range
      zmin = 0.0
      zmax = 2.0

      ! Loop forever
      n = 50
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
         CALL init_halomod(mmin, mmax, scale_factor_z(z), hmod, cosm, verbose)
         CALL print_halomod(hmod, cosm, verbose)

         CALL calculate_HMx_a(fields, nf, k, nk, pow_li, pow_2h, pow_1h, pow_hm, hmod, cosm, verbose=.FALSE., response=.FALSE.)

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

   ELSE IF (imode == 5) THEN

      ! Projection diagnostics

      ! Assigns the cosmological model
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! Set the field types
      ix = -1
      DO i = 1, 2
         CALL set_field_for_xpow(ix(i), ip(i))
      END DO

      ! Fill the projection kernels (plural)
      CALL fill_projection_kernels(ix, proj, cosm)

      DO j = 1, 2

         IF (j == 1) outfile = 'data/nz1.dat'
         IF (j == 2) outfile = 'data/nz2.dat'
         OPEN (7, file=outfile)
         IF (ALLOCATED(proj(j)%nz)) THEN
            CALL write_nz(proj(j), outfile)
         ELSE
            WRITE (7, *) 0., 0.
         END IF
         CLOSE (7)

         IF (j == 1) outfile = 'data/efficiency1.dat'
         IF (j == 2) outfile = 'data/efficiency2.dat'
         OPEN (7, file=outfile)
         IF (ALLOCATED(proj(j)%q)) THEN
            CALL write_efficiency(proj(j), cosm, outfile)
         ELSE
            WRITE (7, *) 0., 0.
         END IF
         CLOSE (7)

         IF (j == 1) outfile = 'data/kernel1.dat'
         IF (j == 2) outfile = 'data/kernel2.dat'
         CALL write_projection_kernel(proj(j), cosm, outfile)

      END DO

   ELSE IF (imode == 6) THEN

      ! n(z) normalisation check

      WRITE (*, *) 'HMx_DRIVER: Checking n(z) functions'
      WRITE (*, *)

      ! Number of n(z) to check
      nnz = 15
      DO i = 1, nnz
         IF (i == 1) nz = tracer_RCSLenS
         IF (i == 2) nz = tracer_CFHTLenS
         IF (i == 3) nz = tracer_KiDS
         IF (i == 4) nz = tracer_KiDS_bin1
         IF (i == 5) nz = tracer_KiDS_bin2
         IF (i == 6) nz = tracer_KiDS_bin3
         IF (i == 7) nz = tracer_KiDS_bin4
         IF (i == 8) nz = tracer_KiDS_450
         IF (i == 9) nz = tracer_KiDS_450_fat_bin1
         IF (i == 10) nz = tracer_KiDS_450_fat_bin2
         IF (i == 11) nz = tracer_KiDS_450_highz
         IF (i == 12) nz = tracer_KiDS_450_bin1
         IF (i == 13) nz = tracer_KiDS_450_bin2
         IF (i == 14) nz = tracer_KiDS_450_bin3
         IF (i == 15) nz = tracer_KiDS_450_bin4
         WRITE (*, *) 'HMx_DRIVER: n(z) number:', nz
         WRITE (*, *)
         CALL read_nz(nz, pro)
         WRITE (*, *) 'HMx_DRIVER: n(z) integral (linear):', integrate_table(pro%z_nz, pro%nz, pro%nnz, 1, pro%nnz, 1)
         WRITE (*, *) 'HMx_DRIVER: n(z) integral (quadratic):', integrate_table(pro%z_nz, pro%nz, pro%nnz, 1, pro%nnz, 2)
         WRITE (*, *) 'HMx_DRIVER: n(z) integral (cubic):', integrate_table(pro%z_nz, pro%nz, pro%nnz, 2, pro%nnz, 3)
         WRITE (*, *)
      END DO

   ELSE IF (imode == 7 .OR. &
            imode == 8 .OR. &
            imode == 9 .OR. &
            imode == 10 .OR. &
            imode == 11 .OR. &
            imode == 37 .OR. &
            imode == 42 .OR. &
            imode == 43 .OR. &
            imode == 47 .OR. &
            imode == 63 .OR. &
            imode == 64) THEN

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
      IF (imode == 37) ix = tracer_CFHTLenS    ! CFHTLenS autospectrum
      IF (imode == 42) ix = tracer_KiDS_450    ! KiDS-450 autospectrum
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
      IF (imode == 37) icosmo = 4 ! WMAP9 fits CFHTLenS okay
      IF (imode == 42 .OR. imode == 43 .OR. imode == 63) icosmo = 4 ! WMAP9
      IF (imode == 47) icosmo = 26 ! Boring with CAMB linear spectrum
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)

      ! Set the k range
      kmin = 1e-3
      kmax = 1e1
      nk = 128
      CALL fill_array(log(kmin), log(kmax), k, nk)
      k = exp(k)

      ! Set the a range
      amin = 0.1 ! Problems with one-halo term if amin is less than 0.1
      amax = 1.
      na = 16
      CALL fill_array(amin, amax, a, na)

      ! Need to call 'comoving_distance' at least once first so as to stop
      ! a write trying to happen while printing to screen
      spam = comoving_distance(1., cosm) ! CARE: This needs to be called before the write-to-screen below
      lmax = 5000.
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
      thmin = 0.01 ! Minium theta [deg]
      thmax = 10.  ! Maxium theta [deg]
      nth = 128
      CALL fill_array(log(thmin), log(thmax), theta, nth)
      theta = exp(theta)
      ALLOCATE (xi(3, nth))

      lmax_xi = 1e5

      WRITE (*, *) 'HMx_DRIVER: Correlation function stuff'
      WRITE (*, *) 'HMx_DRIVER: Minium theta [deg]:', thmin
      WRITE (*, *) 'HMx_DRIVER: Minium theta [deg]:', thmax
      WRITE (*, *) 'HMx_DRIVER: Number of theta:', nth
      WRITE (*, *) 'HMx_DRIVER: lmax for xi:', lmax_xi
      WRITE (*, *)

      ! Output directory
      dir = 'data/'

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
         IF (imode == 37 .OR. imode == 47) ihm = 1  ! HMcode (2016)
         IF (imode == 42 .OR. imode == 43 .OR. imode == 63) ihm = 3 ! 3 - Standard halo model
         IF (imode == 64) ihm = 3 ! Standard
         CALL assign_halomod(ihm, hmod, verbose)
         CALL calculate_HMx(ip, 2, mmin, mmax, k, nk, a, na, pows_li, pows_2h, pows_1h, pows_hm, hmod, cosm, verbose, response=.FALSE.)

         base = TRIM(dir)//'power'
         CALL write_power_a_multiple(k, a, pows_li, pows_2h(1, 2, :, :), pows_1h(1, 2, :, :), pows_hm(1, 2, :, :), nk, na, base, verbose)

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
               pow_ka = pows_li
            ELSE IF (j == 2) THEN
               WRITE (*, *) 'HMx_DRIVER: Doing 2-halo'
               outfile = TRIM(dir)//'cl_2h.dat'
               IF (imode == 37) outfile = TRIM(dir)//'CFHTLenS_cl_2h.dat'
               IF (imode == 47) outfile = TRIM(dir)//'CMBlensing_cl_2h.dat'
               pow_ka = pows_2h(1, 2, :, :)
            ELSE IF (j == 3) THEN
               WRITE (*, *) 'HMx_DRIVER: Doing 1-halo'
               outfile = TRIM(dir)//'cl_1h.dat'
               IF (imode == 37) outfile = TRIM(dir)//'CFHTLenS_cl_1h.dat'
               IF (imode == 47) outfile = TRIM(dir)//'CMBlensing_cl_1h.dat'
               pow_ka = pows_1h(1, 2, :, :)
            ELSE IF (j == 4) THEN
               WRITE (*, *) 'HMx_DRIVER: Doing full'
               outfile = TRIM(dir)//'cl_hm.dat'
               IF (imode == 37) outfile = TRIM(dir)//'CFHTLenS_cl_hm.dat'
               IF (imode == 47) outfile = TRIM(dir)//'CMBlensing_cl_hm.dat'
               pow_ka = pows_hm(1, 2, :, :)
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
               CALL calculate_angular_xi(theta, xi, nth, ell, Cl(1, 2, :), nl, NINT(lmax_xi))
               CALL write_angular_xi(theta, xi, nth, outfile)

            END IF

         END DO
         WRITE (*, *) 'HMx_DRIVER: Done'
         WRITE (*, *)

      ELSE IF (imode == 8) THEN

         ! Assess cross-correlation as a function of cosmology

         ! Allocate array for power
         ALLOCATE (pow_ka(nk, na))

         ! Set range in sigma_8
         sig8min = 0.7
         sig8max = 0.9
         ncos = 5 ! I may have changed this number inadvertantly

         ! Loop over cosmology
         DO i = 1, ncos

            cosm%sig8 = progression(sig8min, sig8max, i, ncos)
            CALL init_cosmology(cosm)
            CALL print_cosmology(cosm)

            CALL assign_halomod(ihm, hmod, verbose)
            CALL calculate_HMx(ip, 2, mmin, mmax, k, nk, a, na, pows_li, pows_2h, pows_1h, pows_hm, hmod, cosm, verbose, response=.FALSE.)

            ! Fill out the projection kernels
            CALL fill_projection_kernels(ix, proj, cosm)

            ! Write to screen
            WRITE (*, *) 'HMx_DRIVER: Computing C(l)'
            WRITE (*, *) 'HMx_DRIVER: ell min:', REAL(ell(1))
            WRITE (*, *) 'HMx_DRIVER: ell max:', REAL(ell(nl))
            WRITE (*, *) 'HMx_DRIVER: number of ell:', nl
            WRITE (*, *)

            ! Loop over all types of C(l) to create
            dir = 'data/'
            base = TRIM(dir)//'cosmology_'
            DO j = 1, 4

               IF (j == 1) THEN
                  WRITE (*, *) 'HMx_DRIVER: Doing C(l) linear'
                  ext = '_cl_linear.dat'
                  pow_ka = pows_li
               ELSE IF (j == 2) THEN
                  WRITE (*, *) 'HMx_DRIVER: Doing C(l) 2-halo'
                  ext = '_cl_2h.dat'
                  pow_ka = pows_2h(1, 2, :, :)
               ELSE IF (j == 3) THEN
                  WRITE (*, *) 'HMx_DRIVER: Doing C(l) 1-halo'
                  ext = '_cl_1h.dat'
                  pow_ka = pows_1h(1, 2, :, :)
               ELSE IF (j == 4) THEN
                  WRITE (*, *) 'HMx_DRIVER: Doing C(l) full'
                  ext = '_cl_hm.dat'
                  pow_ka = pows_hm(1, 2, :, :)
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
         ALLOCATE (pows_li(nk, na), pows_2h(2, 2, nk, na), pows_1h(2, 2, nk, na), pows_hm(2, 2, nk, na))
         ALLOCATE (pow_ka(nk, na))

         DO i = 0, 6
            IF (icumulative .EQV. .FALSE.) THEN
               ! Set the mass intervals
               IF (i == 0) THEN
                  m1 = mmin
                  m2 = mmax
               ELSE IF (i == 1) THEN
                  m1 = mmin
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
                  m1 = mmin
                  m2 = mmax
               ELSE IF (i == 1) THEN
                  m1 = mmin
                  m2 = 1e11
               ELSE IF (i == 2) THEN
                  m1 = mmin
                  m2 = 1e12
               ELSE IF (i == 3) THEN
                  m1 = mmin
                  m2 = 1e13
               ELSE IF (i == 4) THEN
                  m1 = mmin
                  m2 = 1e14
               ELSE IF (i == 5) THEN
                  m1 = mmin
                  m2 = 1e15
               ELSE IF (i == 6) THEN
                  m1 = mmin
                  m2 = 1e16
               END IF
            END IF

            ! Set the code to not 'correct' the two-halo power for missing
            ! mass when doing the calcultion binned in halo mass
            IF ((icumulative .EQV. .TRUE.) .AND. i > 0) hmod%ip2h = 0

            WRITE (*, fmt='(A16)') 'HMx_DRIVER: Mass range'
            WRITE (*, fmt='(A16,I5)') 'HMx_DRIVER: Iteration:', i
            WRITE (*, fmt='(A21,2ES15.7)') 'HMx_DRIVER: M_min [Msun/h]:', m1
            WRITE (*, fmt='(A21,2ES15.7)') 'HMx_DRIVER: M_max [Msun/h]:', m2
            WRITE (*, *)

            !Loop over redshifts
            DO j = 1, na

               !Initiliasation for the halomodel calcualtion
               CALL assign_halomod(ihm, hmod, verbose)
               CALL init_halomod(m1, m2, a(j), hmod, cosm, verbose)
               CALL print_halomod(hmod, cosm, verbose)
               CALL calculate_HMx_a(ip, 2, k, nk, pows_li(:, j), pows_2h(:, :, :, j), pows_1h(:, :, :, j), pows_hm(:, :, :, j), hmod, cosm, verbose, response=.FALSE.)

               !Write progress to screen
               IF (j == 1) THEN
                  WRITE (*, fmt='(A5,A7)') 'i', 'a'
                  WRITE (*, fmt='(A13)') '   ============'
               END IF
               WRITE (*, fmt='(I5,F8.3)') j, a(j)

            END DO
            WRITE (*, fmt='(A13)') '   ============'
            WRITE (*, *)

            dir = 'data/'
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
                  pow_ka = pows_li
               ELSE IF (j == 2) THEN
                  outfile = TRIM(dir)//'cl_2h.dat'
                  ext = '_cl_2h.dat'
                  pow_ka = pows_2h(1, 2, :, :)
               ELSE IF (j == 3) THEN
                  outfile = TRIM(dir)//'cl_1h.dat'
                  ext = '_cl_1h.dat'
                  pow_ka = pows_1h(1, 2, :, :)
               ELSE IF (j == 4) THEN
                  outfile = TRIM(dir)//'cl_hm.dat'
                  ext = '_cl_hm.dat'
                  pow_ka = pows_hm(1, 2, :, :)
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

         CALL assign_halomod(ihm, hmod, verbose)

         CALL calculate_HMx(ip, 2, mmin, mmax, k, nk, a, na, pows_li, pows_2h, pows_1h, pows_hm, hmod, cosm, verbose, response=.FALSE.)

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
            dir = 'data/'
            base = TRIM(dir)//'redshift_'
            mid = '_'
            DO j = 1, 4

               !Set output files
               IF (j == 1) THEN
                  ext = '_cl_linear.dat'
                  outfile = TRIM(dir)//'cl_linear.dat'
                  pow_ka = pows_li
               ELSE IF (j == 2) THEN
                  ext = '_cl_2h.dat'
                  outfile = TRIM(dir)//'cl_2h.dat'
                  pow_ka = pows_2h(1, 2, :, :)
               ELSE IF (j == 3) THEN
                  ext = '_cl_1h.dat'
                  outfile = TRIM(dir)//'cl_1h.dat'
                  pow_ka = pows_1h(1, 2, :, :)
               ELSE IF (j == 4) THEN
                  ext = '_cl_hm.dat'
                  outfile = TRIM(dir)//'cl_hm.dat'
                  pow_ka = pows_hm(1, 2, :, :)
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

   ELSE IF (imode == 12 .OR. imode == 44 .OR. imode == 38 .OR. imode == 39 .OR. imode == 55 .OR. imode == 56) THEN

      ! Triad stuff
      ! 12 - Triad (T_5 cross correlations)
      ! 38 - AGN triad for all BAHAMAS feedback scenarios
      ! 39 - AGN triad for all BAHAMAS feedback scenarios and ell
      ! 44 - Triad for paper (fixed WMAP9 and feedback)
      ! 55 - Triad for all BAHAMAS feedback scenarios
      ! 56 - PAPER: Triad for all BAHAMAS feedback scenarios and ell

      IF (imode == 12 .OR. imode == 44) THEN
         nfeed = 1
      ELSE IF (imode == 38 .OR. imode == 39 .OR. imode == 55 .OR. imode == 56) THEN
         nfeed = 3
      ELSE
         STOP 'HMx_DRIVER: Error, something fucked up'
      END IF

      ! Directory for data output
      dir = 'data'

      ! Assign the cosmological model
      IF (imode == 38 .OR. imode == 39 .OR. imode == 44 .OR. imode == 55 .OR. imode == 56) icosmo = 4 ! WMAP 9 - BAHAMAS
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! Set the ell range
      IF (imode == 12 .OR. imode == 44 .OR. imode == 55) THEN
         lmin = 100.
         lmax = 4000.
         nl = 64
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
      nt = 13
      !nt=5
      ALLOCATE (ixx(nt), ixx_names(nt))
      DO i = 1, nt
         CALL triad_tracers(i, ixx(i), ixx_names(i))
         !CALL triad_4_tracers(i,ixx(i),ixx_names(i))
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

         ! Assign the halo
         CALL assign_halomod(ihm, hmod, verbose)

         IF (j == 1) base = 'triad_Cl_AGN'
         IF (j == 2) base = 'triad_Cl_AGN-lo'
         IF (j == 3) base = 'triad_Cl_AGN-hi'

         CALL xpow(ixx, nt, mmin, mmax, ell, Cl, nl, hmod, cosm, verbose=.TRUE.)

         ! Write data
         DO ii = 1, nt
            DO jj = 1, nt
               outfile = TRIM(dir)//'/'//TRIM(base)//'_'//TRIM(ixx_names(ii))//'-'//TRIM(ixx_names(jj))//'.dat'
               CALL write_Cl(ell, Cl(:, ii, jj), nl, outfile, verbose=.TRUE.)
            END DO
         END DO

      END DO

   ELSE IF (imode == 13) THEN

      ! Calculate a cross-correlation coefficient

      ! Assign the cosmology
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! Assign the halo model
      CALL assign_halomod(ihm, hmod, verbose)

      ! Set the ell range and allocate arrays for l and C(l)
      lmin = 1e0
      lmax = 1e4 ! Strange errors and crashes if this is increased to 10^5
      nl = 64
      CALL fill_array(log(lmin), log(lmax), ell, nl)
      ell = exp(ell)
      ALLOCATE (Cl(nl, 2, 2))

      ! Directory for outputs
      dir = 'data'

      ! Set the necessary fields
      ix = -1
      DO i = 1, 2
         CALL set_field_for_xpow(ix(i), ip(i))
      END DO

      ! Do the cross correlation
      CALL xpow(ix, 2, mmin, mmax, ell, Cl, nl, hmod, cosm, verbose)

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

   ELSE IF (imode == 14 .OR. imode == 33) THEN

      ! Make power spectra variations as a function of baryon parameter variations
      ! 14 - General version
      ! 33 - Paper version

      ! Number of values to try for each parameter
      m = 9

      ! Set number of k points and k range (log spaced)
      kmin = 1e-3
      kmax = 1e1
      nk = 128
      CALL fill_array(log(kmin), log(kmax), k, nk)
      k = exp(k)

      ! Set the fields
      nf = 5
      ALLOCATE (fields(nf))
      fields(1) = field_matter
      fields(2) = field_cdm
      fields(3) = field_gas
      fields(4) = field_star
      fields(5) = field_electron_pressure

      ! Allocate arrays for the power spectra
      ALLOCATE (powd_li(nk), powd_2h(nk), powd_1h(nk), powd_hm(nk))
      ALLOCATE (pow_li(nk), pow_2h(nf, nf, nk), pow_1h(nf, nf, nk), pow_hm(nf, nf, nk))

      ! Set the redshift
      z = 0.0

      ! Assigns the cosmological model
      icosmo = 4
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! Initiliasation for the halo-model calcualtion
      IF (imode == 33) ihm = 3
      CALL assign_halomod(ihm, hmod, verbose)
      CALL init_halomod(mmin, mmax, scale_factor_z(z), hmod, cosm, verbose)
      CALL print_halomod(hmod, cosm, verbose)

      ! DMONLY
      field = field_dmonly
      outfile = 'data/DMONLY.dat'
      CALL calculate_HMx_a(field, 1, k, nk, powd_li, powd_2h, powd_1h, powd_hm, hmod, cosm, verbose, response=.FALSE.)
      CALL write_power(k, powd_li, powd_2h, powd_1h, powd_hm, nk, outfile, verbose)

      ! Prevents warning
      ilog = .FALSE.

      ! Loop over parameters
      DO ipa = 1, param_n

         ! Cycle HMcode parameters
         IF (ipa >= 21 .AND. ipa <= 32) CYCLE

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
            param_min = hmod%mstar/10.
            param_max = hmod%mstar*10.
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
            param_min = hmod%eta-0.5
            param_max = hmod%eta
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

            IF (hmod%HMx_mode == 1 .OR. hmod%HMx_mode == 2 .OR. hmod%HMx_mode == 3) THEN
               IF (ipa == param_alpha) hmod%alpha = param
               IF (ipa == param_eps) hmod%eps = param
               IF (ipa == param_Gamma) hmod%Gamma = param
               IF (ipa == param_M0) hmod%M0 = param
               IF (ipa == param_Astar) hmod%Astar = param
               IF (ipa == param_Twhim) hmod%Twhim = param
               IF (ipa == param_cstar) hmod%cstar = param
               IF (ipa == param_fcold) hmod%fcold = param
               IF (ipa == param_mstar) hmod%mstar = param
               IF (ipa == param_sstar) hmod%sstar = param
               IF (ipa == param_alphap) hmod%alphap = param
               IF (ipa == param_Gammap) hmod%Gammap = param
               IF (ipa == param_cstarp) hmod%cstarp = param
               IF (ipa == param_fhot) hmod%fhot = param
               IF (ipa == param_alphaz) hmod%alphaz = param
               IF (ipa == param_Gammaz) hmod%Gammaz = param
               IF (ipa == param_M0z) hmod%M0z = param
               IF (ipa == param_Astarz) hmod%Astarz = param
               IF (ipa == param_Twhimz) hmod%Twhimz = param
               IF (ipa == param_eta) hmod%eta = param
               IF (ipa == param_epsz) hmod%epsz = param
               IF (ipa == param_beta) hmod%beta = param
               IF (ipa == param_betap) hmod%betap = param
               IF (ipa == param_betaz) hmod%betaz = param
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

            CALL init_halomod(mmin, mmax, scale_factor_z(z), hmod, cosm, verbose)
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
            CALL write_mass_fractions(hmod, cosm, outfile)

            CALL calculate_HMx_a(fields, nf, k, nk, pow_li, pow_2h, pow_1h, pow_hm, hmod, cosm, verbose, response=.FALSE.)

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

   ELSE IF (imode == 17) THEN

      ! 3D spectra for a user choice of fields

      ! Assigns the cosmological model
      icosmo = 1
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      CALL assign_halomod(ihm, hmod, verbose)

      ! Set number of k points and k range (log spaced)
      ! The range kmin=1e-3 to kmax=1e4 is necessary to compare to HMcode
      nk = 128
      kmin = 1e-3
      kmax = 1e2
      CALL fill_array(log(kmin), log(kmax), k, nk)
      k = exp(k)

      !Set the number of redshifts and range (linearly spaced) and convert z -> a
      nz = 16
      zmin = 0.
      zmax = 4.
      CALL fill_array(zmin, zmax, a, nz)
      a = 1./(1.+a)
      na = nz

      ! Choose the field types
      nf = 2
      DO i = 1, nf
         WRITE (*, *) 'HMx_driver: Choose halo', i
         CALL set_halo_type(ip(i))
      END DO

      ! User chooses halo model
      CALL calculate_HMx(ip, nf, mmin, mmax, k, nk, a, na, pows_li, pows_2h, pows_1h, pows_hm, hmod, cosm, verbose, response=.FALSE.)

      base = 'data/power'
      CALL write_power_a_multiple(k, a, pows_li, pows_2h(ip(1), ip(2), :, :), pows_1h(ip(1), ip(2), :, :), pows_hm(ip(1), ip(2), :, :), nk, na, base, verbose)

   ELSE IF (imode == 18 .OR. imode == 40 .OR. imode == 48) THEN

      ! Create 3D bias function
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
      nk = 128
      kmin = 1e-3
      kmax = 1e2
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
      nf = 2
      ALLOCATE (fields(nf))
      fields(1) = field_dmonly

      ! Select field type for bias study
      IF (imode == 18) CALL set_halo_type(fields(2))
      IF (imode == 40) fields(2) = 9  ! Central galaxies/haloes
      IF (imode == 48) fields(2) = 12 ! HI

      CALL calculate_HMx(fields, nf, mmin, mmax, k, nk, a, na, pows_li, pows_2h, pows_1h, pows_hm, hmod, cosm, verbose, response=.FALSE.)

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

            CALL write_power_a_multiple(k, a, pows_li, pows_2h(j1, j2, :, :), pows_1h(j1, j2, :, :), pows_hm(j1, j2, :, :), nk, na, base, verbose)

         END DO
      END DO

   ELSE IF (imode == 19) THEN

      ! Create CCL benchmark data

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
            CALL init_halomod(mmin, mmax, scale_factor_z(z), hmod, cosm, verbose)
            CALL print_halomod(hmod, cosm, verbose)

            ! Do the halo-model calculation
            field = field_dmonly
            CALL calculate_HMx_a(field, 1, k, nk, pow_li, pow_2h, pow_1h, pow_hm, hmod, cosm, verbose, response=.FALSE.)

            ! Write out the results
            CALL write_power(k, pow_li, pow_2h(1, 1, :), pow_1h(1, 1, :), pow_hm(1, 1, :), nk, outfile, verbose)

         END DO

      END DO

   ELSE IF (imode == 20) THEN

      ! Ma et al. Fig. 1

      ! Set the cosmology
      icosmo = 3
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! Set the halo model
      z = 0.0
      ihm = 4
      CALL assign_halomod(ihm, hmod, verbose)
      CALL init_halomod(mmin, mmax, scale_factor_z(z), hmod, cosm, verbose)
      CALL print_halomod(hmod, cosm, verbose)

      ! Make the Figure
      CALL YinZhe_Fig1(hmod, cosm)

   ELSE IF (imode == 21) THEN

      ! Stuff for diagnosing problems with the window function integrand

      outfile = 'winint/integrand.dat'
      irho = 11
      rv = 1.
      c = 4.
      rs = rv/c
      p1 = 1.18
      p2 = 0.
      rmin = 0.
      rmax = rv
      CALL winint_diagnostics(rmin, rmax, rv, rs, p1, p2, irho, outfile)

   ELSE IF (imode == 22) THEN

      ! Speed tests for W(M,k) integrals

      ! k range
      kmin = 1e-2
      kmax = 1e3
      nk = 512
      CALL fill_array(log(kmin), log(kmax), k, nk)
      k = exp(k)

      ! Halo parameters
      rv = 1.    ! Virial radius
      c = 4.     ! Concentration
      rs = rv/c
      p1 = 1.2   ! Gamma
      p2 = 0.
      irho = 11  ! KS density profile
      rmin = 0.
      rmax = rv

      CALL winint_speed_tests(k, nk, rmin, rmax, rv, rs, p1, p2, irho)

   ELSE IF (imode == 23) THEN

      ! Mead (2017) dark energy results

      ! Set number of k points and k range (log spaced)
      nk = 128
      kmin = 1e-3
      kmax = 1e1
      CALL fill_array(log(kmin), log(kmax), k, nk)
      k = exp(k)
      ALLOCATE (pow_li(nk), pow_2h(1, 1, nk), pow_1h(1, 1, nk), pow_hm(1, 1, nk))

      ! Directory for output
      dir = 'data'

      ! Set the halo model
      ihm = 12

      ! Number of different cosmologies
      ncos = 25

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
         CALL assign_halomod(ihm, hmod, verbose)
         CALL init_halomod(mmin, mmax, scale_factor_z(z), hmod, cosm, verbose)
         CALL print_halomod(hmod, cosm, verbose)

         ! Do the halo-model calculation
         field = field_dmonly
         CALL calculate_HMx_a(field, 1, k, nk, pow_li, pow_2h, pow_1h, pow_hm, hmod, cosm, verbose, response=.FALSE.)

         ! Write out the results
         outfile = TRIM(dir)//'/power_'//TRIM(outfile)//'.dat'
         CALL write_power(k, pow_li, pow_2h(1, 1, :), pow_1h(1, 1, :), pow_hm(1, 1, :), nk, outfile, verbose)

      END DO

   ELSE IF (imode == 24) THEN

      ! Projection tests
      ! TODO: Fix why these tests are super slow if ihm=3 and Cl_yy looks really fucked up
      ! TODO: I think this is probably because of the pressure evolution in ihm=3

      ! Number of tests
      nx = 3

      ! Initially assume tests pass
      ifail = .FALSE.

      ALLOCATE (ixx(nx))
      ixx(1) = tracer_RCSLenS
      ixx(2) = tracer_Compton_y
      ixx(3) = tracer_CMB_lensing

      ALLOCATE (ixx_names(nx))
      ixx_names(1) = 'RCSLenS'
      ixx_names(2) = 'y'
      ixx_names(3) = 'CMB'

      ! Set the ell range (should probably read this in from a file)
      lmin = 1e0
      lmax = 1e4
      nl = 128
      CALL fill_array(log(lmin), log(lmax), ell, nl)
      ell = exp(ell)
      ALLOCATE (Cl_bm(nl), Cl(nl, nx, nx))

      ! Set the cosmology
      icosmo = 1

      ! Assigns the cosmological model
      CALL assign_cosmology(icosmo, cosm, verbose=.FALSE.)
      CALL init_cosmology(cosm)

      ! Set the halo model
      ihm = 18 ! Incredibly slow if I set ihm=3

      ! Assign the halo model
      CALL assign_halomod(ihm, hmod, verbose=.FALSE.)

      ! Do the cross correlation
      CALL xpow(ixx, nx, mmin, mmax, ell, Cl, nl, hmod, cosm, verbose=.TRUE.)

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

      WRITE (*, *) 'HMx_DRIVER: Limber tests should take ~22s to run'
      IF (ifail) THEN
         WRITE (*, *) 'HMx_DRIVER: Limber tests failed'
      ELSE
         WRITE (*, *) 'HMx_DRIVER: Limber tests passed'
      END IF
      WRITE (*, *)

   ELSE IF (imode == 25) THEN

      ! Halo-void model

      ! Set number of k points and k range (log spaced)
      nk = 128
      kmin = 1e-3
      kmax = 1e2
      CALL fill_array(log(kmin), log(kmax), k, nk)
      k = exp(k)

      ! Assigns the cosmological model
      icosmo = 1 ! 1 - Boring
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! Sets the redshift
      z = 0.0

      ! Initiliasation for the halomodel calcualtion
      ihm = 16 ! 16 - Halo-void model
      CALL assign_halomod(ihm, hmod, verbose)
      CALL init_halomod(mmin, mmax, scale_factor_z(z), hmod, cosm, verbose)
      CALL print_halomod(hmod, cosm, verbose)

      ! Allocate arrays
      ALLOCATE (powd_li(nk), powd_2h(nk), powd_1h(nk), powd_hm(nk))

      ! Do the halo-model calculation
      field = field_dmonly
      CALL calculate_HMx_a(field, 1, k, nk, powd_li, powd_2h, powd_1h, powd_hm, hmod, cosm, verbose, response=.FALSE.)

      ! Write the one-void term if necessary
      OPEN (7, file='data/power_halovoid.dat')
      DO i = 1, nk
         WRITE (7, *) k(i), powd_li(i), powd_2h(i), powd_1h(i), powd_hm(i), p_1void(k(i), hmod)
      END DO
      CLOSE (7)

   ELSE IF (imode == 26) THEN

      ! Automated testing

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
         CALL assign_cosmology(icosmo, cosm, verbose_tests)
         CALL init_cosmology(cosm)
         CALL print_cosmology(cosm)

         ! Assign the halo model
         CALL assign_halomod(ihm, hmod, verbose_tests)
         !CALL print_halomod(hmod)

         field = field_dmonly
         CALL calculate_HMx(field, 1, mmin, mmax, k, nk, a, na, pows_li, pows_2h, pows_1h, pows_hm, hmod, cosm, verbose_tests, response=.FALSE.)

         ! Loop over k and a
         error_max = 0.
         DO j = 1, na
            DO i = 1, nk
               error = ABS(-1.+pows_hm(1, 1, i, j)/pow_ka(i, j))
               IF (error > error_max) error_max = error
               IF (error > tolerance) THEN
                  WRITE (*, *) 'HMx_DRIVER: Test:', itest
                  WRITE (*, *) 'HMx_DRIVER: Wavenumber [h/Mpc]:', k(i)
                  WRITE (*, *) 'HMx_DRIVER: Scale-factor:', a(j)
                  WRITE (*, *) 'HMx_DRIVER: Expected power:', pow_ka(i, j)
                  WRITE (*, *) 'HMx_DRIVER: Model power:', pows_hm(1, 1, i, j)
                  WRITE (*, *) 'HMx_DRIVER: Tolerance:', tolerance
                  WRITE (*, *) 'HMx_DRIVER: Error:', error
                  WRITE (*, *)
                  ifail = .TRUE.
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
         CALL write_power_a_multiple(k, a, pows_li, pows_2h(1, 1, :, :), pows_1h(1, 1, :, :), pows_hm(1, 1, :, :), nk, na, base, verbose_tests)

         DEALLOCATE (k, a, pow_ka)

      END DO

      IF (ifail) THEN
         STOP 'HMx_DRIVER: Error, tests failed'
      ELSE
         WRITE (*, *) 'HMx_DRIVER: Tests should take around 0.45 seconds to run'
         WRITE (*, *) 'HMx_DRIVER: Tests passed'
         WRITE (*, *)
      END IF

   ELSE IF (imode == 27 .OR. imode == 28 .OR. imode == 29 .OR. imode == 30 .OR. imode == 70 .OR. imode == 71) THEN

      ! Comparison with FrankenEmu or Mira Titan

      ! Number of cosmological models (+1)
      IF (imode == 28 .OR. imode == 30) n = 37 ! Franken Emu
      IF (imode == 27 .OR. imode == 29) n = 36 ! Mira Titan
      IF (imode == 70 .OR. imode == 71) n = 37 ! Cosmic Emu

      ! Set number of redshifts
      IF (imode == 28 .OR. imode == 30) THEN
         nz = 6 ! Franken Emu
      ELSE IF (imode == 27 .OR. imode == 29) THEN
         nz = 4 ! Mira Titan
      ELSE IF (imode == 70 .OR. imode == 71) THEN
         nz = 3 ! Cosmic Emu
      ELSE
         STOP 'HMx_DRIVER: Error, imode not specified correctly for EMU'
      END IF

      ! Set the random number generator for the random cosmologies
      IF (imode == 30 .OR. imode == 29 .OR. imode == 71) THEN
         CALL RNG_set(seed=0)
      END IF

      ! Allocate arrays
      na = nz
      ALLOCATE (zs(nz), a(nz))

      ! Set redshifts
      IF (imode == 28 .OR. imode == 30) THEN
         ! Franken Emu (z up to 4)
         zs(1) = 0.0
         zs(2) = 0.5
         zs(3) = 1.0
         zs(4) = 2.0
         zs(5) = 3.0
         zs(6) = 4.0
      ELSE IF (imode == 27 .OR. imode == 29) THEN
         ! Mira Titan (z up to 2)
         zs(1) = 0.0
         zs(2) = 0.5
         zs(3) = 1.0
         zs(4) = 2.0
      ELSE IF (imode == 70 .OR. imode == 71) THEN
         ! Cosmic Emu (z up to 1)
         zs(1) = 0.0
         zs(2) = 0.5
         zs(3) = 1.0
      ELSE
         STOP 'HMx_DRIVER: Error, imode not specified correctly for EMU'
      END IF

      ! Set scale factors
      DO i = 1, na
         a(i) = scale_factor_z(zs(i))
      END DO

      ! Output files
      base = 'data/cosmo'
      mid = '_z'
      ext = '.dat'

      ! Initiliasation for the halomodel calcualtion
      CALL assign_halomod(ihm, hmod, verbose)

      ! Loop over cosmologies
      DO i = 0, n

         IF (imode == 27) THEN
            icosmo = 100+i ! Mira Titan nodes
         ELSE IF (imode == 28) THEN
            icosmo = 200+i ! Franken Emu nodes (same as cosmic emu nodes)
         ELSE IF (imode == 70) THEN
            icosmo = 300+i ! Cosmic Emu nodes
         ELSE IF (imode == 29) THEN
            icosmo = 24    ! Random Mira Titan
         ELSE IF (imode == 30) THEN
            icosmo = 25    ! Random Franken Emu
         ELSE IF (imode == 71) THEN
            icosmo = 38    ! Random Cosmic Emu
         ELSE
            STOP 'HMx_DRIVER: Error, imode not specified correctly for EMU'
         END IF
         CALL assign_cosmology(icosmo, cosm, verbose)
         CALL init_cosmology(cosm)
         CALL print_cosmology(cosm)

         ! Loop over redshift
         DO j = 1, nz

            IF (imode == 27 .OR. imode == 29) THEN
               CALL get_Mira_Titan_power(k_sim, pow_sim, nk, zs(j), cosm, rebin=.FALSE.)
            ELSE IF (imode == 28 .OR. imode == 30) THEN
               CALL get_Franken_Emu_power(k_sim, pow_sim, nk, zs(j), cosm, rebin=.FALSE.)
            ELSE IF (imode == 70 .OR. imode == 71) THEN
               CALL get_Cosmic_Emu_power(k_sim, pow_sim, nk, zs(j), cosm, rebin=.FALSE.)
            ELSE
               STOP 'HMx_DRIVER: Error, imode not specified correctly for EMU'
            END IF

            ALLOCATE (pow_li(nk), pow_2h(1, 1, nk), pow_1h(1, 1, nk), pow_hm(1, 1, nk))

            CALL init_halomod(mmin, mmax, a(j), hmod, cosm, verbose=.TRUE.)
            CALL print_halomod(hmod, cosm, verbose=.TRUE.)
            field = field_dmonly
            CALL calculate_HMx_a(field, 1, k_sim, nk, pow_li, pow_2h, pow_1h, pow_hm, hmod, cosm, verbose=.FALSE., response=.FALSE.)

            ALLOCATE (pow_ql(nk), pow_oh(nk), pow_hf(nk))
            CALL calculate_halofit_a(k_sim, a(j), pow_li, pow_ql, pow_oh, pow_hf, nk, cosm, verbose=.TRUE., ihf=4)

            ! Write data
            outfile = number_file2(base, i, mid, j, ext)
            OPEN (7, file=outfile)
            DO ii = 1, nk
               WRITE (7, *) k_sim(ii), pow_li(ii), pow_sim(ii), &
                  pow_ql(ii), pow_oh(ii), pow_hf(ii), &
                  pow_2h(1, 1, ii), pow_1h(1, 1, ii), pow_hm(1, 1, ii)
            END DO
            CLOSE (7)

            DEALLOCATE (pow_li, pow_2h, pow_1h, pow_hm)
            DEALLOCATE (pow_ql, pow_oh, pow_hf)
            DEALLOCATE (k_sim, pow_sim)

         END DO

      END DO

   ELSE IF (imode == 31 .OR. imode == 53) THEN

      ! Power breakdown as a function of mass (paper)

      ! Set number of k points and k range (log spaced)
      nk = 128
      kmin = 1e-3
      kmax = 1e2
      CALL fill_array(log(kmin), log(kmax), k, nk)
      k = exp(k)

      ! Set the fields
      nf = 2
      ALLOCATE (fields(2))
      fields(1) = field_matter
      fields(2) = field_electron_pressure

      ! Allocate arrays for power
      ALLOCATE (pow_li(nk), pow_2h(nf, nf, nk), pow_1h(nf, nf, nk), pow_hm(nf, nf, nk))

      ! Assign the cosmological model
      icosmo = 4 ! BAHAMAS cosmology
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! Set the halo model
      ihm = 3
      CALL assign_halomod(ihm, hmod, verbose)

      ! Set the redshift
      z = 0.0

      IF (imode == 31) THEN
         imin = 10
         imax = 16
      ELSE IF (imode == 53) THEN
         imin = 10
         imax = 16
      ELSE
         STOP 'HMx_DRIVER: Error, imode specified incorrectly'
      END IF

      ! Loop over upper limit of mass integral
      DO i = imin, imax

         ! Set the upper limit for the mass integration
         ! Needs to be 10. to enforce that m2 is real
         IF (imode == 31) THEN
            m1 = mmin
            m2 = 10.**i
         ELSE IF (imode == 53) THEN
            IF (i == imax) THEN
               m1 = mmin
               m2 = mmax
            ELSE
               m1 = 10.**i
               m2 = 10.**(i+1)
            END IF
         ELSE
            STOP 'HMx_DRIVER: Error, imode specified incorrectly'
         END IF

         ! Initiliasation for the halomodel calcualtion
         CALL init_halomod(m1, m2, scale_factor_z(z), hmod, cosm, verbose)
         CALL print_halomod(hmod, cosm, verbose)

         ! Do the halo-model calculation
         CALL calculate_HMx_a(fields, nf, k, nk, pow_li, pow_2h, pow_1h, pow_hm, hmod, cosm, verbose, response=.FALSE.)

         ! Loop over fields
         DO j1 = 1, nf
            DO j2 = j1, nf

               ! Write out the results
               base = 'data/power_'
               IF (imode == 31) THEN
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

   ELSE IF (imode == 34) THEN

      ! Write out the variation of HMx parameters with T_AGN and z

      !Assign the cosmological model
      icosmo = 4
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      !Sets the redshift
      z = 0.0

      !Initiliasation for the halomodel calcualtion
      ihm = 18
      CALL assign_halomod(ihm, hmod, verbose)
      CALL init_halomod(mmin, mmax, scale_factor_z(z), hmod, cosm, verbose)
      CALL print_halomod(hmod, cosm, verbose)

      zmin = 0.
      zmax = 4.
      nz = 101

      mass = 1e14 ! Mass to write out for [Msun/h]

      DO j = 1, 3

         IF (j == 1) THEN
            outfile = 'data/HMx_params_AGN-lo.dat'
            hmod%Theat = 10**7.6
         ELSE IF (j == 2) THEN
            outfile = 'data/HMx_params_AGN.dat'
            hmod%Theat = 10**7.8
         ELSE IF (j == 3) THEN
            outfile = 'data/HMx_params_AGN-hi.dat'
            hmod%Theat = 10**8.0
         END IF

         OPEN (7, file=outfile)
         DO i = 1, nz
            hmod%z = progression(zmin, zmax, i, nz)
            WRITE (7, *) hmod%z, HMx_alpha(mass, hmod), HMx_eps(hmod), HMx_Gamma(mass, hmod), HMx_M0(hmod), HMx_Astar(mass, hmod), HMx_Twhim(hmod)
         END DO
         CLOSE (7)

      END DO

   ELSE IF (imode == 35) THEN

      ! Look at the effect of cores on halo profiles

      ! Set number of k points and k range (log spaced)
      nk = 128
      kmin = 1e-3
      kmax = 1e2
      CALL fill_array(log(kmin), log(kmax), k, nk)
      k = exp(k)

      ! Allocate arrays for the power calculation
      ALLOCATE (pow_li(nk), pow_2h(1, 1, nk), pow_1h(1, 1, nk), pow_hm(1, 1, nk))

      ! Assign the cosmological model
      icosmo = 1 ! Boring
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! Set the redshift
      z = 0.0

      ! Initiliasation for the halomodel calcualtion
      ihm = 21
      CALL assign_halomod(ihm, hmod, verbose)

      ! Range of cores to explore
      rcore_min = 0.
      rcore_max = 0.1
      ncore = 16

      DO i = 1, ncore

         ! Set the core radius
         hmod%rcore = progression(rcore_min, rcore_max, i, ncore)

         ! Initialise the halo-model calculation
         IF (i == 1) THEN
            verbose2 = verbose
         ELSE
            verbose2 = .FALSE.
         END IF
         CALL init_halomod(mmin, mmax, scale_factor_z(z), hmod, cosm, verbose2)
         CALL print_halomod(hmod, cosm, verbose2)

         ! Do the halo-model calculation
         field = field_dmonly
         CALL calculate_HMx_a(field, 1, k, nk, pow_li, pow_2h, pow_1h, pow_hm, hmod, cosm, verbose=.FALSE., response=.FALSE.)

         ! Write out the results
         base = 'data/power_cored_'
         ext = '.dat'
         outfile = number_file(base, i, ext)
         CALL write_power(k, pow_li, pow_2h(1, 1, :), pow_1h(1, 1, :), pow_hm(1, 1, :), nk, outfile, verbose)

      END DO

   ELSE IF (imode == 36) THEN

      ! Automated testing of hydro models

      ! Assigns the cosmological model
      icosmo = 4
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! Set the field combinations
      nf = 5
      ALLOCATE (fields(nf))
      fields(1) = field_matter
      fields(2) = field_cdm
      fields(3) = field_gas
      fields(4) = field_star
      fields(5) = field_electron_pressure

      ! Set the redshifts
      na = 4
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
      CALL assign_halomod(ihm, hmod, verbose)
      CALL calculate_HMx(fields, nf, mmin, mmax, k, nk, a, na, pows_li, pows_2h, pows_1h, pows_hm, hmod, cosm, verbose, response=.FALSE.)

      ! Initially assume that all the tests will pass
      ifail = .FALSE.

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
               CALL write_power(k, pows_li(:, j), pows_2h(itest, jtest, :, j), pows_1h(itest, jtest, :, j), pows_hm(itest, jtest, :, j), nk, outfile, verbose)

               ! Set the error max to be zero before each test
               error_max = 0.
               verbose_tests = .TRUE.

               ! Loop over k and check for accuracy and write out data
               DO i = 1, nk

                  ! This is what I use as the error
                  error = ABS(-1.+pows_hm(itest, jtest, i, j)/powb_hm(itest, jtest, i, j))

                  ! Update the maximium error if it is exceeded
                  IF (error > error_max) error_max = error

                  ! If the test has failed then write out this diagnostic stuff
                  IF (verbose_tests .AND. error > tolerance) THEN
                     WRITE (*, *) 'HMx_DRIVER: Test failing'
                     WRITE (*, *) 'HMx_DRIVER: Test fields:', fields(itest), fields(jtest)
                     WRITE (*, *) 'HMx_DRIVER: Wavenumber [h/Mpc]:', k(i)
                     WRITE (*, *) 'HMx_DRIVER: Redshift:', redshift_a(a(j))
                     WRITE (*, *) 'HMx_DRIVER: Benchmark power:', powb_hm(itest, jtest, i, j)
                     WRITE (*, *) 'HMx_DRIVER: HMx power:', pows_hm(itest, jtest, i, j)
                     WRITE (*, *) 'HMx_DRIVER: Tolerance:', tolerance
                     WRITE (*, *) 'HMx_DRIVER: Error:', error
                     WRITE (*, *)
                     verbose_tests = .FALSE.
                     ifail = .TRUE.
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
         WRITE (*, *) 'HMx_DRIVER: Hydro tests should take around 2.00 seconds to run'
         WRITE (*, *) 'HMx_DRIVER: Hydro tests passed'
      END IF
      WRITE (*, *)

   ELSE IF (imode == 41) THEN

      ! Assess 3D power as a function of cosmology

      ! Assign the cosmology
      icosmo = 1
      CALL assign_cosmology(icosmo, cosm, verbose)

      ! Assign the halo model
      ihm = 3
      CALL assign_halomod(ihm, hmod, verbose=.FALSE.)

      ! Set the redshift
      z = 0.

      ! Set range in sigma_8
      sig8min = 0.7
      sig8max = 0.9
      ncos = 15

      ! Allocate arrays for the fields
      nf = 2
      ALLOCATE (fields(nf))
      fields(1) = field_matter
      fields(2) = field_electron_pressure

      ! Set number of k points and k range (log spaced)
      nk = 128
      kmin = 1e-3
      kmax = 1e2
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
         CALL init_halomod(mmin, mmax, scale_factor_z(z), hmod, cosm, verbose=.FALSE.)
         CALL print_halomod(hmod, cosm, verbose=.FALSE.)

         ! Do the halo-model calculation
         CALL calculate_HMx_a(fields, 2, k, nk, pow_li, pow_2h, pow_1h, pow_hm, hmod, cosm, verbose=.FALSE., response=.FALSE.)

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

   ELSE IF (imode == 45) THEN

      ! Comparison of power spectra from Sheth & Tormen (1999) vs. Tinker et al. (2010) mass function

      ! Set number of k points and k range (log spaced)
      nk = 128
      kmin = 1e-3
      kmax = 1e2
      CALL fill_array(log(kmin), log(kmax), k, nk)
      k = exp(k)

      ! Allocate arrays for power
      ALLOCATE (pow_li(nk), pow_2h(1, 1, nk), pow_1h(1, 1, nk), pow_hm(1, 1, nk))

      ! Assigns the cosmological model
      icosmo = 1
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! Sets the redshift
      z = 0.0

      DO j = 1, 3

         IF (j == 1) THEN
            ! Sheth & Tormen (1999) mass function and bias
            ihm = 3
            outfile = 'data/power_ShethTormen.dat'
         ELSE IF (j == 2) THEN
            ihm = 23 ! Tinker et al. (2010) mass function and bias
            outfile = 'data/power_Tinker.dat'
         ELSE IF (j == 3) THEN
            ihm = 27 ! Press & Schecter (1974) mass function and bias
            outfile = 'data/power_PressSchecter.dat'
         END IF

         IF (j == 1) THEN
            verbose2 = verbose
         ELSE
            verbose2 = .FALSE.
         END IF

         !Initiliasation for the halomodel calcualtion
         CALL assign_halomod(ihm, hmod, verbose2)
         CALL init_halomod(mmin, mmax, scale_factor_z(z), hmod, cosm, verbose2)
         CALL print_halomod(hmod, cosm, verbose2)

         ! Do the halo-model calculation
         field = field_dmonly
         CALL calculate_HMx_a(field, 1, k, nk, pow_li, pow_2h, pow_1h, pow_hm, hmod, cosm, verbose2, response=.FALSE.)

         ! Write out data
         CALL write_power(k, pow_li, pow_2h(1, 1, :), pow_1h(1, 1, :), pow_hm(1, 1, :), nk, outfile, verbose)

      END DO

   ELSE IF (imode == 46) THEN

      ! Mass function plots for different mass functions

      ! Assigns the cosmological model
      icosmo = 1
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! Sets the redshift
      z = 0.0

      !Initiliasation for the halomodel calcualtion
      nhm = 3
      ALLOCATE (hmods(nhm))
      DO i = 1, nhm
         IF (i == 1) THEN
            ihm = 3  ! 3 - Sheth & Tormen mass function, virial mass haloes
         ELSE IF (i == 2) THEN
            ihm = 27 ! 27 - Press & Schecter mass function, virial mass haloes
         ELSE IF (i == 3) THEN
            ihm = 23 ! 23 - Tinker mass function, virial mass haloes
         ELSE
            STOP 'HMx_DRIVER: Error, something went wrong setting halo model'
         END IF
         CALL assign_halomod(ihm, hmods(i), verbose)
         CALL init_halomod(mmin, mmax, scale_factor_z(z), hmods(i), cosm, verbose)
         CALL print_halomod(hmods(i), cosm, verbose)
      END DO

      ! Range in nu
      numin = 0.1
      numax = 6.
      n = 256

      ! Loop over nu and write out mass function
      ! TODO: Actually mass function, n(M), for different theory models, currenly only hmods(1)
      OPEN (7, file='data/bnu_functions.dat')
      OPEN (8, file='data/gnu_functions.dat')
      OPEN (9, file='data/mass_functions.dat')
      DO i = 1, n
         nu = progression(numin, numax, i, n)
         !mass=M_nu(nu,hmod)
         !mf=mass_function(mass,hmod,cosm)
         mass = M_nu(nu, hmods(1))
         mf = mass_function(mass, hmods(1), cosm)
         !WRITE(7,*) nu, mass, b_ps(nu,hmod), b_st(nu,hmod), b_Tinker(nu,hmod)
         !WRITE(8,*) nu, mass, g_ps(nu,hmod), g_st(nu,hmod), g_Tinker(nu,hmod)
         !WRITE(9,*) nu, mass, mf, mf*mass/comoving_matter_density(cosm), mf*mass**2/comoving_matter_density(cosm)
         WRITE (7, *) nu, mass, (b_nu(nu, hmods(j)), j=1, nhm)
         WRITE (8, *) nu, mass, (g_nu(nu, hmods(j)), j=1, nhm)
         WRITE (9, *) nu, mass, mf, mf*mass/comoving_matter_density(cosm), mf*mass**2/comoving_matter_density(cosm)
      END DO
      CLOSE (7)
      CLOSE (8)
      CLOSE (9)

   ELSE IF (imode == 49) THEN

      ! HI mass fractions

      ! Assigns the cosmological model
      icosmo = 1
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! Set the range of redshifts
      z1 = 0.
      z2 = 5.
      nz = 6
      ALLOCATE (hmods(nz), HI_frac(nz), z_tab(nz), a(nz))

      ! Initiliasation for the halomodel calcualtion at all zs
      ihm = 25
      DO j = 1, nz
         z_tab(j) = progression(z1, z2, j, nz)
         a(j) = scale_factor_z(z_tab(j))
         CALL assign_halomod(ihm, hmods(j), verbose)
         CALL init_halomod(mmin, mmax, a(j), hmods(j), cosm, verbose)
         CALL print_halomod(hmods(j), cosm, verbose)
      END DO

      ! Set the range of halo masses
      m1 = mmin
      m2 = mmax
      n = 256

      ! Loop over mass and z and do the calculation
      OPEN (7, file='data/HI_mass_fraction.dat')
      DO i = 1, n
         mass = exp(progression(log(m1), log(m2), i, n))
         DO j = 1, nz
            HI_frac(j) = halo_HI_fraction(mass, hmods(j), cosm)
         END DO
         WRITE (7, *) mass, (HI_frac(j), j=1, nz)
      END DO
      CLOSE (7)
      CLOSE (8)

   ELSE IF (imode == 50) THEN

      ! Mass function plots as Lbox is varied

      ! Assigns the cosmological model
      icosmo = 1
      CALL assign_cosmology(icosmo, cosm, verbose=.FALSE.)

      ! Minimum and maximum box sizes [Mpc/h] and number of cosmologies
      Lmin = 32.
      Lmax = 2048.
      ncos = 8

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

      ! Set the redshift
      z = 0.0

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
         CALL init_halomod(mmin, mmax, scale_factor_z(z), hmods(i), cosms(i), verbose2)
      END DO
      CALL print_halomod(hmods(1), cosms(1), verbose)

      ! Range in nu
      numin = 0.1
      numax = 6.
      n = 256

      ALLOCATE (masses(ncos))
      OPEN (7, file='data/nu_mass_Lbox.dat')
      DO i = 1, n
         DO j = 1, ncos
            nu = progression(numin, numax, i, n)
            masses(j) = M_nu(nu, hmods(j))
            !mf, mf*mass/comoving_matter_density(cosm), mf*mass**2/comoving_matter_density(cosm)
         END DO
         WRITE (7, *) nu, (masses(j), j=1, ncos)
      END DO
      CLOSE (7)

   ELSE IF (imode == 51) THEN

      ! Assigns the cosmological model
      icosmo = 1
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! Set number of k points and k range (log spaced)
      nk = 128
      kmin = 1e-3
      kmax = 1e2
      CALL fill_array(log(kmin), log(kmax), k, nk)
      k = exp(k)

      ! Set the number of redshifts and range (linearly spaced) and convert z -> a
      na = 16
      amin = 0.1
      amax = 1.
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

         CALL assign_halomod(ihm, hmod, verbose)

         CALL calculate_HMx(field, 1, mmin, mmax, k, nk, a, na, pows_li, pows_2h, pows_1h, pows_hm, hmod, cosm, verbose, response=.FALSE.)

         CALL write_power_a_multiple(k, a, pows_li, pows_2h(1, 1, :, :), pows_1h(1, 1, :, :), pows_hm(1, 1, :, :), nk, na, base, verbose)

      END DO

      cmin = 1.
      cmax = 10.
      n = 256
      outfile = 'data/p_conc.dat'
      cbar = 4.
      OPEN (7, file=outfile)
      DO i = 1, n
         c = progression(cmin, cmax, i, n)
         WRITE (7, *) c, lognormal(c, cbar, hmod%dlnc)
      END DO
      CLOSE (7)

   ELSE IF (imode == 54) THEN

      ! Set number of k points and k range (log spaced)
      nk = 64
      kmin = 1e-3
      kmax = 1e2
      CALL fill_array(log(kmin), log(kmax), k, nk)
      k = k*real(nk)/real(nk+1) ! To stop plot spilling over border
      k = exp(k)

      ! Set the redshift
      z = 0.

      ! Set the fields
      ALLOCATE (fields(2))
      fields = field_dmonly

      ! Assigns the cosmological model
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! Assign the halo model
      CALL assign_halomod(ihm, hmod, verbose)
      CALL init_halomod(mmin, mmax, scale_factor_z(z), hmod, cosm, verbose)
      CALL print_halomod(hmod, cosm, verbose)

      ! Write data
      OPEN (7, file='data/trispectrum.dat')
      DO j = 1, nk
         DO i = 1, nk
            WRITE (7, *) k(i), k(j), T_1h(k(i), k(j), fields, hmod, cosm)
         END DO
      END DO
      CLOSE (7)

   ELSE IF (imode == 57 .OR. imode == 61 .OR. imode == 62 .OR. imode == 65 .OR. imode == 66) THEN

      ! Calculate C(l) by direct integration of the measured 3D BAHAMAS power spectra
      ! 57 - Triad 3
      ! 61 - Triad 4
      ! 62 - Triad 5
      ! 65 - Triad 5 with contribution to Limber integrals per ell
      ! 66 - Triad 5 but with extended ell range

      ! Assigns the cosmological model
      icosmo = 4 ! 4 - WMAP9
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! Choose the set of redshifts to use
      !na=nz_BAHAMAS
      !CALL BAHAMAS_zs(a,na)

      ! Set the mesh size used for the measured P(k)
      mesh = mesh_BAHAMAS

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
      nsim = 5
      DO j = 1, nsim

         ! Set BAHAMAS models
         IF (j == 1) THEN
            name = 'AGN_TUNED_nu0'
            outbase = 'data/triad_Cl_direct_AGN'
         ELSE IF (j == 2) THEN
            name = 'AGN_7p6_nu0'
            outbase = 'data/triad_Cl_direct_AGN-lo'
         ELSE IF (j == 3) THEN
            name = 'AGN_8p0_nu0'
            outbase = 'data/triad_Cl_direct_AGN-hi'
         ELSE IF (j == 4) THEN
            name = 'AGN_TUNED_nu0_v2'
            outbase = 'data/triad_Cl_direct_AGN_v2'
         ELSE IF (j == 5) THEN
            name = 'AGN_TUNED_nu0_v3'
            outbase = 'data/triad_Cl_direct_AGN_v3'
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
                  CALL read_BAHAMAS_power(k, pow_sim, nk, redshift_a(a(i)), name, mesh, ip, cosm, &
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

   ELSE IF (imode == 58) THEN

      ! Check to see power is the same

      ! Assigns the cosmological model
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! Assign the halo model
      CALL assign_halomod(ihm, hmod, verbose)

      ! Set number of k points and k range (log spaced)
      nk = 32
      kmin = 1e-3
      kmax = 1e2
      CALL fill_array(log(kmin), log(kmax), k, nk)
      k = exp(k)

      ! Set the scale factor and range (linearly spaced)
      na = 4
      amin = 0.1
      amax = 1.0
      CALL fill_array(amin, amax, a, na)

      nf = 2
      ALLOCATE (fields(nf))
      !fields=field_matter
      !fields=field_electron_pressure
      fields = field_gas

      CALL calculate_HMx(fields, 2, mmin, mmax, k, nk, a, na, pows_li, pows_2h, pows_1h, pows_hm, hmod, cosm, verbose, response=.FALSE.)

      ! Write to screen to check they are the same
      WRITE (*, *) 'HMx_DRIVER: All these columns should be idential'
      DO i = 1, nk
         WRITE (*, *) pows_hm(1, 1, i, na), pows_hm(1, 2, i, na), pows_hm(2, 1, i, na), pows_hm(2, 2, i, na)
      END DO
      WRITE (*, *) 'HMx_DRIVER: Done'
      WRITE (*, *)

   ELSE IF (imode == 59) THEN

      ! Produce data for Fig. 1 of Tinker et al. (2010)

      ! Assigns the cosmological model
      icosmo = 1
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ihm = 44
      CALL assign_halomod(ihm, hmod, verbose)
      CALL init_halomod(mmin, mmax, 1., hmod, cosm, verbose)

      ! Range in nu
      numin = 0.1
      numax = 10.
      n = 256

      outfile = 'data/Tinker_bias.dat'
      OPEN (7, file=outfile)
      DO i = 1, n
         nu = exp(progression(log(numin), log(numax), i, n))
         WRITE (7, *) nu, b_nu(nu, hmod)
      END DO
      CLOSE (7)

   ELSE IF (imode == 60) THEN

      ! Make data for Limber comparison with CCL

      ! Assigns the cosmological model
      icosmo = 1 ! 1 - Boring
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! k range and array
      kmin = kmin_xpow
      kmax = kmax_xpow
      nk = nk_xpow
      CALL fill_array(kmin, kmax, k, nk)

      ! a range and array
      amin = amin_xpow
      amax = amax_xpow
      na = na_xpow
      CALL fill_array(amin, amax, a, na)

      ! Allocate and fill array with linear power; P(k,a)
      ALLOCATE (pow_ka(nk, na))
      DO j = 1, na
         DO i = 1, nk
            pow_ka(i, j) = p_lin(k(i), a(j), cosm)
         END DO
      END DO

      ! l range
      lmin = 2.
      lmax = 10000.
      nl = nint(lmax)-nint(lmin)+1 ! Get all ell
      CALL fill_array(lmin, lmax, ell, nl)
      ALLOCATE (Cl_bm(nl))

      ! Set lensing tracer and call C(l) routine
      ix = tracer_CFHTLenS
      CALL xpow_pka(ix, ell, Cl_bm, nl, k, a, pow_ka, nk, na, cosm)

      ! Write data
      outfile = 'data/HMx_Cl_test.dat'
      OPEN (7, file=outfile)
      DO i = 1, nl
         WRITE (7, *) ell(i), Cl_bm(i)
      END DO
      CLOSE (7)

   ELSE IF (imode == 67) THEN

      ! 67 - Non-linear halo bias model

      ! Assigns the cosmological model
      icosmo = 200 ! 200 - Frankenemu M000
      !icosmo=5 ! 5 - WMAP5 (Multidark)
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! Set the redshift
      z = 0.0

      ! Get the emulator power
      CALL get_Franken_Emu_power(k, pow_sim, nk, z, cosm, rebin=.FALSE.)
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
            STOP 'HMX_DRIVER: Error, something went wrong in 68'
         END IF

         ! Initiliasation for the halomodel calcualtion
         CALL assign_halomod(ihm, hmod, verbose)
         CALL init_halomod(mmin, mmax, scale_factor_z(z), hmod, cosm, verbose)
         CALL print_halomod(hmod, cosm, verbose)

         ! Allocate arrays
         ALLOCATE (powd_li(nk), powd_2h(nk), powd_1h(nk), powd_hm(nk))

         ! Do the halo-model calculation
         field = field_dmonly
         CALL calculate_HMx_a(field, 1, k, nk, powd_li, powd_2h, powd_1h, powd_hm, hmod, cosm, verbose, response=.FALSE.)

         ! Write out the results
         CALL write_power(k, powd_li, powd_2h, powd_1h, powd_hm, nk, outfile, verbose)

         ! Deallocate arrays
         DEALLOCATE (powd_li, powd_2h, powd_1h, powd_hm)

      END DO

   ELSE IF (imode == 68) THEN

      ! Non-linear halo bias integrand

      icosmo = 37 ! 37 - WMAP 5
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! Initiliasation for the halomodel calcualtion
      !ihm=48 ! 48 - Non-linear halo bias with Sheth & Torman
      ihm = 49 ! 49 - Non-linear halo bias with Tinker
      CALL assign_halomod(ihm, hmod, verbose)
      CALL init_halomod(mmin, mmax, scale_factor_z(z), hmod, cosm, verbose)
      CALL print_halomod(hmod, cosm, verbose)

      ! Set k range and allocate
      kmin = 0.1
      kmax = 1.0
      nk = 10
      CALL fill_array(kmin, kmax, k, nk)

      ! Set nu range and allocate
      numin = 0.
      numax = 3.
      nnu = 100
      CALL fill_pixels(numin, numax, nu_tab, nnu)

      ! Needs to be called to prevent unit conflict below
      B_NL = BNL(k(1), nu_tab(1), nu_tab(1), hmod)

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
               nu1 = nu_tab(i)
               nu2 = nu_tab(j)
               B_NL = BNL(k(ik), nu1, nu2, hmod) ! Needed because otherwise function writes
               I_NL = B_NL*g_nu(nu1, hmod)*b_nu(nu1, hmod)*g_nu(nu2, hmod)*b_nu(nu2, hmod)
               WRITE (7, *) nu1, nu2, B_NL, I_NL
            END DO
         END DO
         CLOSE (7)
         WRITE (*, *) 'HMx_DRIVER: Done'
         WRITE (*, *)

      END DO

   ELSE IF (imode == 69) THEN

      ! 69 - Halo power to compare with Multidark

      icosmo = 37 ! 37 - WMAP 5
      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

      ! Set the redshift
      na = 5
      ALLOCATE (a(na))
      a(1) = 1.00000 ! Should be 1.00109, but rounded to 1 exactly to prevent weird negative-z stuff
      a(2) = 0.68215
      a(3) = 0.59103
      a(4) = 0.49990
      a(5) = 0.25690

      ! Actual range of a values to do
      ia1 = 1
      ia2 = 5

      ALLOCATE (bases(na))
      bases(1) = 'data/power_a1.00109_hh_'
      bases(2) = 'data/power_a0.68215_hh_'
      bases(3) = 'data/power_a0.59103_hh_'
      bases(4) = 'data/power_a0.49990_hh_'
      bases(5) = 'data/power_a0.25690_hh_'

      ! Allocate k range
!!$     !kmin=1e-3
!!$     !kmax=1e1
!!$     !nk=101
!!$     !CALL fill_array(log(kmin),log(kmax),k,nk)
!!$     !k=exp(k)
      infile = '/Users/Mead/Physics/data/Multidark/power/M512/halo_1.00109_bin0_bin0_power.dat'
      nk = file_length(infile)
      ALLOCATE (k(nk))
      OPEN (7, file=infile)
      DO i = 1, nk
         READ (7, *) k(i)
      END DO
      CLOSE (7)

      ! Set the fields
      nf = 4
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
            hmod%n = 512 ! 512 is at least necessary due to thin halo bins (lots of zero points)
            mmin_bnl = 1e12 ! 1e7 is the default, be really careful here
            mmax_bnl = 1e16 ! 1e17 is the default, be really careful here
            CALL init_halomod(mmin_bnl, mmax_bnl, a(i), hmod, cosm, verbose)
            CALL print_halomod(hmod, cosm, verbose)

            ! Do the halo-model calculation
            CALL calculate_HMx_a(fields, nf, k, nk, pow_li, pow_2h, pow_1h, pow_hm, hmod, cosm, verbose, response=.FALSE.)

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

   ELSE

      STOP 'HMx_DRIVER: Error, you have specified the mode incorrectly'

   END IF

CONTAINS

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

   SUBROUTINE triad_tracers(i, ix, ix_name)

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

   END SUBROUTINE triad_tracers

   SUBROUTINE read_k_values(infile, k, nk)

      ! Get k values from a data file, assumes they are the first column of a data file
      IMPLICIT NONE
      CHARACTER(len=256), INTENT(IN) :: infile ! Input file location
      REAL, ALLOCATABLE, INTENT(OUT) :: k(:) ! Output array of k values
      INTEGER, INTENT(OUT) :: nk ! Output number of k values

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

   SUBROUTINE YinZhe_Fig1(hmod, cosm)

      IMPLICIT NONE
      TYPE(halomod), INTENT(INOUT) :: hmod   ! Halo model
      TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmological model
      INTEGER :: i
      REAL :: r, rs, rv, c, Mh, rh, r500c, m500c

      REAL, PARAMETER :: M = 1e15    ! Halo virial? mass [Msun]
      REAL, PARAMETER :: rmin = 1e-3 ! Minimum radius [Mpc]
      REAL, PARAMETER :: rmax = 8    ! Maximum radius [Mpc]
      INTEGER, PARAMETER :: nr = 512 ! Number of points in radius

      LOGICAL, PARAMETER :: real_space = .TRUE.
      INTEGER, PARAMETER :: itype = field_electron_pressure ! electron pressure

      IF (hmod%has_mass_conversions .EQV. .FALSE.) CALL convert_mass_definitions(hmod, cosm)

      Mh = M*cosm%h ! This virial mass is now [Msun/h]

      rv = exp(find(log(Mh), hmod%log_m, log(hmod%rv), hmod%n, 3, 3, 2)) ! [Mpc/h]
      c = find(log(Mh), hmod%log_m, hmod%c, hmod%n, 3, 3, 2)
      rs = rv/c ! [Mpc/h]

      m500c = exp(find(log(Mh), hmod%log_m, log(hmod%m500c), hmod%n, 3, 3, 2)) ! [Mpc/h]
      r500c = exp(find(log(Mh), hmod%log_m, log(hmod%r500c), hmod%n, 3, 3, 2)) ! [Mpc/h]

      WRITE (*, *) 'YINZHE_FIG1: Making data for this figure'
      WRITE (*, *) 'YINZHE_FIG1: Redshift:', z
      WRITE (*, *) 'YINZHE_FIG1: Virial radius [Mpc]:', rv/cosm%h
      WRITE (*, *) 'YINZHE_FIG1: Virial radius [Mpc/h]:', rv
      WRITE (*, *) 'YINZHE_FIG1: r_500,c [Mpc]:', r500c/cosm%h
      WRITE (*, *) 'YINZHE_FIG1: r_500,c [Mpc/h]:', r500c
      WRITE (*, *) 'YINZHE_FIG1: r_500,c / r_v:', r500c/rv
      WRITE (*, *) 'YINZHE_FIG1: Virial halo mass [log10 Msun]:', log10(M)
      WRITE (*, *) 'YINZHE_FIG1: Virial halo mass [log10 Msun/h]:', log10(Mh)
      WRITE (*, *) 'YINZHE_FIG1: M_500,c [log10 Msun]:', log10(M500c/cosm%h)
      WRITE (*, *) 'YINZHE_FIG1: M_500,c [log10 Msun/h]:', log10(M500c)
      WRITE (*, *) 'YINZHE_FIG1: M_500,c / M_v:', M500c/Mh
      WRITE (*, *) 'YINZHE_FIG1: Halo concentraiton:', c

      OPEN (7, file='data/YinZhe_Fig1.dat')
      DO i = 1, nr
         r = progression(rmin, rmax, i, nr) ! Radius [Mpc]
         rh = r*cosm%h ! Convert [Mpc/h]
         WRITE (7, *) r, UPP(real_space, rh, Mh, rv, rs, hmod, cosm)*r**2, win_type(real_space, itype, rh, Mh, rv, rs, hmod, cosm)*r**2
      END DO
      CLOSE (7)

      WRITE (*, *) 'YINZHE_FIG1: Done'
      WRITE (*, *)

   END SUBROUTINE YinZhe_Fig1

   SUBROUTINE write_power(k, pow_lin, pow_2h, pow_1h, pow, nk, output, verbose)

      IMPLICIT NONE
      CHARACTER(len=*), INTENT(IN) :: output
      INTEGER, INTENT(IN) :: nk
      REAL, INTENT(IN) :: k(nk), pow_lin(nk), pow_2h(nk), pow_1h(nk), pow(nk)
      LOGICAL, INTENT(IN) :: verbose
      INTEGER :: i

      IF (verbose) WRITE (*, *) 'WRITE_POWER: Writing power to ', TRIM(output)

      ! Loop over k values
      ! Fill the tables with one- and two-halo terms as well as total
      OPEN (7, file=output)
      DO i = 1, nk
         WRITE (7, fmt='(5ES20.10)') k(i), pow_lin(i), pow_2h(i), pow_1h(i), pow(i)
      END DO
      CLOSE (7)

      IF (verbose) THEN
         WRITE (*, *) 'WRITE_POWER: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE write_power

   SUBROUTINE write_power_a_multiple(k, a, pow_lin, pow_2h, pow_1h, pow_full, nk, na, base, verbose)

      IMPLICIT NONE
      CHARACTER(len=*), INTENT(IN) :: base
      INTEGER, INTENT(IN) :: nk, na
      REAL, INTENT(IN) :: k(nk), a(na), pow_lin(nk, na), pow_2h(nk, na), pow_1h(nk, na), pow_full(nk, na)
      LOGICAL, INTENT(IN) :: verbose
      REAL :: pow(nk, na)
      INTEGER :: i
      CHARACTER(len=512) :: output
      LOGICAL :: verbose2

      DO i = 1, 4
         IF (i == 1) THEN
            output = TRIM(base)//'_linear.dat'
            pow = pow_lin
         ELSE IF (i == 2) THEN
            output = TRIM(base)//'_2h.dat'
            pow = pow_2h
         ELSE IF (i == 3) THEN
            output = TRIM(base)//'_1h.dat'
            pow = pow_1h
         ELSE IF (i == 4) THEN
            output = TRIM(base)//'_hm.dat'
            pow = pow_full
         ELSE
            STOP 'WRITE_POWER_A_MULTIPLE: Error, something went FUBAR'
         END IF
         IF (i == 1) THEN
            verbose2 = verbose
         ELSE
            verbose2 = .FALSE.
         END IF
         CALL write_power_a(k, a, pow, nk, na, output, verbose2)
      END DO

   END SUBROUTINE write_power_a_multiple

   SUBROUTINE write_power_a(k, a, pow, nk, na, output, verbose)

      IMPLICIT NONE
      CHARACTER(len=*), INTENT(IN) :: output
      INTEGER, INTENT(IN) :: nk, na
      REAL, INTENT(IN) :: k(nk), a(na), pow(nk, na)
      LOGICAL, INTENT(IN) :: verbose
      INTEGER :: i, j

      ! Print to screen
      IF (verbose) THEN
         WRITE (*, *) 'WRITE_POWER_A: The first entry of the file is hashes - #####'
         WRITE (*, *) 'WRITE_POWER_A: The remainder of the first row are the scale factors - a'
         WRITE (*, *) 'WRITE_POWER_A: The remainder of the first column are the wave numbers - k'
         WRITE (*, *) 'WRITE_POWER_A: Each row then gives the power at that k and a'
         WRITE (*, *) 'WRITE_POWER_A: Output:', TRIM(output)
      END IF

      ! Write out data to files
      OPEN (7, file=output)
      DO i = 0, nk
         IF (i == 0) THEN
            WRITE (7, fmt='(A20,40F20.10)') '#####', (a(j), j=1, na)
         ELSE
            WRITE (7, fmt='(F20.10,40E20.10)') k(i), (pow(i, j), j=1, na)
         END IF
      END DO
      CLOSE (7)

      ! Print to screen
      IF (verbose) THEN
         WRITE (*, *) 'WRITE_POWER_A: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE write_power_a

!!$  SUBROUTINE write_distances(cosm)
!!$
!!$    ! Write file of z vs. r(z)
!!$    IMPLICIT NONE
!!$    TYPE(cosmology), INTENT(INOUT) :: cosm
!!$    CHARACTER(len=256) :: output
!!$    INTEGER :: i
!!$    REAL :: z
!!$
!!$    ! Now write the results of r(z) calculation
!!$    output='data/distance.dat'
!!$    WRITE(*,*) 'WRITE_DISTANCE: Writing r(a): ', TRIM(output)
!!$    OPEN(7,file=output)
!!$    DO i=1,cosm%n_r
!!$       z=redshift_a(cosm%a_r(i))
!!$       WRITE(7,*) z, cosm%r(i), f_k(cosm%r(i),cosm)
!!$    END DO
!!$    CLOSE(7)
!!$    WRITE(*,*) 'WRITE_DISTANCE: Done'
!!$    WRITE(*,*)
!!$
!!$  END SUBROUTINE write_distances

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

   SUBROUTINE xpow(ix, nx, mmin, mmax, ell, Cl, nl, hmod, cosm, verbose)

      ! Calculates the C(l) for the cross correlation of fields ix(1) and ix(2)
      ! TODO: Speed up if there are repeated fields in ix(n) (should this ever happen?)
      ! TODO: Add bin theory option
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: ix(nx)
      INTEGER, INTENT(IN) :: nx
      REAL, INTENT(IN) :: mmin, mmax
      REAL, INTENT(IN) :: ell(nl)
      REAL, INTENT(OUT) :: Cl(nl, nx, nx)
      INTEGER, INTENT(IN) :: nl
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      LOGICAL, INTENT(IN) :: verbose
      REAL, ALLOCATABLE :: a(:), k(:), pow_li(:, :), pow_2h(:, :, :, :), pow_1h(:, :, :, :), pow_hm(:, :, :, :)
      REAL :: lmin, lmax
      INTEGER :: ixx(2), ip(nx), nk, na, i, j, nnx, match(nx)
      INTEGER, ALLOCATABLE :: iix(:)
      REAL, ALLOCATABLE :: uCl(:, :, :)

      !IF(repeated_entries(ix,n)) STOP 'XPOW: Error, repeated tracers'
      CALL unique_index(ix, nx, iix, nnx, match)
      ALLOCATE (uCl(nl, nnx, nnx))

      ! Set the k range
      nk = nk_xpow
      CALL fill_array(log(kmin_xpow), log(kmax_xpow), k, nk)
      k = exp(k)

      ! Set the a range
      na = na_xpow
      CALL fill_array(amin_xpow, amax_xpow, a, na)

      ! Set the ell range
      lmin = ell(1)
      lmax = ell(nl)

      ! Write to screen
      IF (verbose) THEN
         WRITE (*, *) 'XPOW: Cross-correlation information'
         WRITE (*, *) 'XPOW: Number of tracers:', nx
         WRITE (*, *) 'XPOW: P(k) minimum k [h/Mpc]:', REAL(kmin_xpow)
         WRITE (*, *) 'XPOW: P(k) maximum k [h/Mpc]:', REAL(kmax_xpow)
         WRITE (*, *) 'XPOW: Number of k:', nk
         WRITE (*, *) 'XPOW: minimum a:', REAL(amin_xpow)
         WRITE (*, *) 'XPOW: maximum a:', REAL(amax_xpow)
         WRITE (*, *) 'XPOW: number of a:', na
         WRITE (*, *) 'XPOW: minimum ell:', REAL(lmin)
         WRITE (*, *) 'XPOW: maximum ell:', REAL(lmax)
         WRITE (*, *) 'XPOW: number of ell:', nl
         WRITE (*, *)
      END IF

      ! Use the xpowlation type to set the necessary halo profiles
      DO i = 1, nnx
         CALL set_field_for_xpow(ix(i), ip(i))
      END DO

      ! Do the halo model power spectrum calculation
      CALL calculate_HMx(ip, nnx, mmin, mmax, k, nk, a, na, pow_li, pow_2h, pow_1h, pow_hm, hmod, cosm, verbose, response=.FALSE.)

      ! Set the Cl to zero initially
      uCl = 0.

      ! Loop over triangle combinations
      DO i = 1, nnx
         ixx(1) = ix(i)
         DO j = i, nnx
            ixx(2) = ix(j)
            CALL xpow_pka(ixx, ell, uCl(:, i, j), nl, k, a, pow_hm(i, j, :, :), nk, na, cosm)
         END DO
      END DO

      ! Fill the symmetric cross terms
      DO i = 1, nnx
         DO j = i+1, nnx
            uCl(:, j, i) = uCl(:, i, j)
         END DO
      END DO

      ! Now fill the full arrays from the unique arrays
      DO i = 1, nx
         DO j = 1, nx
            ii = match(i)
            jj = match(j)
            Cl(:, i, j) = uCl(:, ii, jj)
         END DO
      END DO

      ! Write to screen
      IF (verbose) THEN
         WRITE (*, *) 'XPOW: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE xpow

   SUBROUTINE set_field_for_xpow(ix, ip)

      ! Set the cross-correlation type
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: ix ! Tracer for cross correlation
      INTEGER, INTENT(OUT) :: ip   ! Corresponding field for power spectrum
      INTEGER :: j

      IF (ix == -1) THEN
         WRITE (*, *) 'SET_FIELDS_FOR_XPOW: Choose field'
         WRITE (*, *) '================================='
         DO j = 1, n_tracers
            WRITE (*, fmt='(I3,A3,A30)') j, '- ', TRIM(xcorr_type(j))
         END DO
         READ (*, *) ix
         WRITE (*, *) '==================================='
         WRITE (*, *)
      END IF

      IF (ix == tracer_Compton_y) THEN
         ! Compton y
         ip = field_electron_pressure
      ELSE IF (ix == tracer_gravity_wave) THEN
         ! Gravitational waves
         ip = field_dmonly
      ELSE IF (ix == tracer_RCSLenS .OR. &
               ix == tracer_CFHTLenS .OR. &
               ix == tracer_CMB_lensing .OR. &
               ix == tracer_KiDS .OR. &
               ix == tracer_KiDS_bin1 .OR. &
               ix == tracer_KiDS_bin2 .OR. &
               ix == tracer_KiDS_bin3 .OR. &
               ix == tracer_KiDS_bin4 .OR. &
               ix == tracer_KiDS_450 .OR. &
               ix == tracer_KiDS_450_fat_bin1 .OR. &
               ix == tracer_KiDS_450_fat_bin2 .OR. &
               ix == tracer_KiDS_450_highz .OR. &
               ix == tracer_lensing_z1p00 .OR. &
               ix == tracer_lensing_z0p75 .OR. &
               ix == tracer_lensing_z0p50 .OR. &
               ix == tracer_lensing_z0p25 .OR. &
               ix == tracer_KiDS_450_bin1 .OR. &
               ix == tracer_KiDS_450_bin2 .OR. &
               ix == tracer_KiDS_450_bin3 .OR. &
               ix == tracer_KiDS_450_bin4) THEN
         ! Lensing
         ip = field_matter
      ELSE IF (ix == tracer_CIB_353) THEN
         ip = field_CIB_353
      ELSE IF (ix == tracer_CIB_545) THEN
         ip = field_CIB_545
      ELSE IF (ix == tracer_CIB_857) THEN
         ip = field_CIB_857
      ELSE IF (ix == tracer_galaxies) THEN
         ip = field_central_galaxies
      ELSE
         STOP 'SET_FIELD_FOR_XPOW: Error, tracer specified incorrectly'
      END IF

   END SUBROUTINE set_field_for_xpow

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

!!$  SUBROUTINE BAHAMAS_zs(a,n)
!!$
!!$    IMPLICIT NONE
!!$    REAL, ALLOCATABLE, INTENT(OUT) :: a(:)
!!$    INTEGER, INTENT(IN) :: n
!!$
!!$    IF(ALLOCATED(a)) DEALLOCATE(a)
!!$    ALLOCATE(a(n))
!!$
!!$    IF(n==4) THEN
!!$       a(1)=scale_factor_z(2.0)
!!$       a(2)=scale_factor_z(1.0)
!!$       a(3)=scale_factor_z(0.5)
!!$       a(4)=scale_factor_z(0.0)
!!$    ELSE IF(n==11) THEN
!!$       a(1)=scale_factor_z(2.0)
!!$       a(2)=scale_factor_z(1.75)
!!$       a(3)=scale_factor_z(1.5)
!!$       a(4)=scale_factor_z(1.25)
!!$       a(5)=scale_factor_z(1.0)
!!$       a(6)=scale_factor_z(0.75)
!!$       a(7)=scale_factor_z(0.5)
!!$       a(8)=scale_factor_z(0.375)
!!$       a(9)=scale_factor_z(0.25)
!!$       a(10)=scale_factor_z(0.125)
!!$       a(11)=scale_factor_z(0.0)
!!$    ELSE IF(n==15) THEN
!!$       a(1)=scale_factor_z(3.0)
!!$       a(2)=scale_factor_z(2.75)
!!$       a(3)=scale_factor_z(2.5)
!!$       a(4)=scale_factor_z(2.25)
!!$       a(5)=scale_factor_z(2.0)
!!$       a(6)=scale_factor_z(1.75)
!!$       a(7)=scale_factor_z(1.5)
!!$       a(8)=scale_factor_z(1.25)
!!$       a(9)=scale_factor_z(1.0)
!!$       a(10)=scale_factor_z(0.75)
!!$       a(11)=scale_factor_z(0.5)
!!$       a(12)=scale_factor_z(0.375)
!!$       a(13)=scale_factor_z(0.25)
!!$       a(14)=scale_factor_z(0.125)
!!$       a(15)=scale_factor_z(0.0)
!!$    ELSE
!!$       STOP 'BAHAMAS_ZS: Error, nz specified incorrectly'
!!$    END IF
!!$
!!$  END SUBROUTINE BAHAMAS_zs

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

      !DO j=1,na
      !   WRITE(*,*) j, a(j)
      !END DO
      !STOP

   END SUBROUTINE add_highz_BAHAMAS

END PROGRAM HMx_driver
