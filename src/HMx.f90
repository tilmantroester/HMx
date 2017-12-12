MODULE cosdef

  TYPE cosmology
     !Contains cosmological parameters that need only be calculated once
     REAL :: om_m, om_b, om_v, om_c, h, n, sig8, w, wa, om_nu
     REAL :: om, k, z_cmb, om_r, T_cmb
     REAL :: A
     REAL, ALLOCATABLE :: logsigma(:), logr_sigma(:)
     REAL, ALLOCATABLE :: growth(:), a_growth(:)
     REAL, ALLOCATABLE :: r(:), a_r(:)
     INTEGER :: nsig, ng, nr
     CHARACTER(len=256) :: name
     !Varying parameters
     INTEGER :: np=5
     REAL :: param(5), param_defaults(5), param_min(5), param_max(5)
     CHARACTER(len=256) :: param_names(5)
     LOGICAL :: param_log(5)
  END TYPE cosmology

  TYPE tables
     !Halo-model stuff that needs to be recalculated for each new z
     REAL, ALLOCATABLE :: c(:), rv(:), nu(:), sig(:), zc(:), m(:), rr(:), sigf(:)
     REAL, ALLOCATABLE :: r500(:), m500(:), c500(:), r200(:), m200(:), c200(:)
     REAL, ALLOCATABLE :: r500c(:), m500c(:), c500c(:), r200c(:), m200c(:), c200c(:)
     REAL, ALLOCATABLE :: log_m(:)
     REAL :: sigv, sigv100, c3, knl, rnl, neff, sig8z
     REAL :: gmin, gmax, gbmin, gbmax
     INTEGER :: n
  END TYPE tables

  TYPE projection
     !Projection quantities that need to be calculated only once
     !These relate to the Limber integrals
     REAL, ALLOCATABLE :: X(:), r_X(:)
     INTEGER :: nX
  END TYPE projection

  !Possibly this could usefully be merged with projection
  TYPE lensing
     !Quantities that are necessary for lensing specifically
     REAL, ALLOCATABLE :: q(:), r_q(:)
     REAL, ALLOCATABLE :: nz(:), z_nz(:)
     INTEGER :: nq, nnz
  END TYPE lensing

END MODULE cosdef

MODULE HMx

  !Module usage statements
  USE cosdef
  USE constants
  USE random_numbers
  USE array_operations
  USE file_info
  USE special_functions
  USE interpolate
  USE string_operations
  USE calculus
  USE calculus_table
  
  !Parameter definitions
  IMPLICIT NONE
  REAL :: plin
  REAL, ALLOCATABLE :: k(:), a(:)
  REAL, ALLOCATABLE :: pow_lin(:), pow_2h(:), pow_1h(:), pow_full(:)
  REAL, ALLOCATABLE :: powa(:,:), powa_lin(:,:), powa_2h(:,:), powa_1h(:,:), powa_full(:,:)
  REAL, ALLOCATABLE :: ell(:), Cell(:), theta(:), xi(:,:)
  INTEGER :: i, j, nk, na, j1, j2, n, nl, nz, nth, nnz, m, ipa
  INTEGER :: ip(2), ix(2), ixx(2)
  REAL :: kmin, kmax, amin, amax, lmin, lmax, thmin, thmax, zmin, zmax
  REAL :: z, z1, z2, r1, r2
  TYPE(cosmology) :: cosm
  TYPE(tables) :: lut
  TYPE(projection) :: proj(2)
  TYPE(lensing) :: lens
  CHARACTER(len=256) :: outfile, base, mid, ext, dir, name
  CHARACTER(len=256) :: mode
  INTEGER :: imode, icosmo, iowl
  REAL :: sig8min, sig8max
  INTEGER :: ncos
  REAL :: m1, m2, mass
  !REAL :: rv, rs, rmax

  !Varying parameters
  !INTEGER, PARAMETER :: np=5
  !REAL :: param(np), param_defaults(np), param_min(np), param_max(np)
  !REAL :: mass
  !CHARACTER(len=256) :: param_names(np)
  !INTEGER :: ipa
  !LOGICAL :: param_log(np)

  !Mathematical constants
  !REAL, PARAMETER :: pi=3.141592654
  !REAL, PARAMETER :: onethird=0.3333333333
  REAL, PARAMETER :: zero=0.

  !Halo-model Parameters
  INTEGER, PARAMETER :: imf=2 !Set mass function (1 - PS, 2 - ST)
  !INTEGER :: ihm=1 !Set verbosity
  LOGICAL :: verbose=.TRUE.
  INTEGER, PARAMETER :: imead=0 !Set to do Mead et al. (2015,2016) accurate calculation
  REAL, PARAMETER :: mmin=1e7 !Minimum halo mass for the calculation
  REAL, PARAMETER :: mmax=1e17 !Maximum halo mass for the calculation
  INTEGER :: ip2h=2 !Method to 'correct' the 2-halo integral
  INTEGER, PARAMETER :: ibias=1 !Bias order to go to
  INTEGER, PARAMETER :: ibox=0 !Consider the simulation volume
  REAL, PARAMETER :: Lbox=400. !Simulation box size
  INTEGER, PARAMETER :: icumulative=1 !Do cumlative distributions for breakdown
  LOGICAL, PARAMETER :: ixi=.FALSE. !Do correlation functions from C(l)
  LOGICAL, PARAMETER :: ifull=.FALSE. !Do only full halo model C(l), xi(theta) calculations
  REAL, PARAMETER :: acc=1e-4 !Global integration-accuracy parameter
  LOGICAL, PARAMETER :: void=.FALSE. !Do voids or not
  REAL, PARAMETER :: lcorr=0.5 !Set to zero for k=l/fk or 0.5 for k=(l+0.5)/fk

  !Physical constants
  REAL, PARAMETER :: yfac=8.125561e-16 !sigma_T/m_e*c^2 in SI
  REAL, PARAMETER :: kb=1.38065e-23 !Boltzmann constant in SI
  !REAL, PARAMETER :: epn=0.875 !1/mu_e electron per nucleon (~8/7 for x=0.25 He mass frac)
  !REAL, PARAMETER :: mue=1.14 !mu_e nucleons per electron (~7/8 for x=0.25 He mass frac)
  REAL, PARAMETER :: fh=0.76 !Hydrogen mass fraction
  REAL, PARAMETER :: mue=2./(1.+fh) !Nucleons per electron (~1.136 if fh=0.76)
  REAL, PARAMETER :: pfac=(5.*fh+3.)/(2.*(fh+1.)) !Pressure factor (Hill & Pajer 2013; I do not understand; ~1.932 if fh=0.76)
  REAL, PARAMETER :: conH0=2998. !(c/H0) in Mpc/h
  REAL, PARAMETER :: mp=1.6726219e-27 !Proton mass in kg
  REAL, PARAMETER :: msun=1.989e30 ! kg/Msun
  REAL, PARAMETER :: mpc=3.086e22 ! m/Mpc
  REAL, PARAMETER :: bigG=6.67408e-11 !Gravitational constant in kg^-1 m^3 s^-2
  REAL, PARAMETER :: eV=1.60218e-19 !Electronvolt in Joules
  REAL, PARAMETER :: cm=0.01 !Centimetre in metres
  REAL, PARAMETER :: rad2deg=180./pi !Radians-to-degrees conversion
  REAL, PARAMETER :: neff=3.046 !Effective number of neutrinos
  REAL, PARAMETER :: critical_density=2.775d11 !Critical density at z=0

  !Name parameters (cannot do PARAMETER with mixed length strings)
  CHARACTER(len=256) :: halo_type(-1:8), xcorr_type(10)!, kernel_type(3)

CONTAINS

  SUBROUTINE init_HMx()

    IMPLICIT NONE

    !Names of variable parameters
    cosm%param_names(1)='alpha' !alpha in virial temperature (turbulence?)
    cosm%param_names(2)='Dc' !Change in NFW concentration due to gas
    cosm%param_names(3)='Gamma' !Gamma in Komatsu-Seljak profile
    cosm%param_names(4)='M_B' !Halo mass at which free and unbound gas are equal
    cosm%param_names(5)='A_*' !Prefactor for stellar fraction

    !Set default values for variable parameters
    cosm%param(1)=1.
    cosm%param(2)=0.
    cosm%param(3)=1.18
    cosm%param(4)=1.2d14
    cosm%param(5)=0.02

    !Set a protected set of defaults
    cosm%param_defaults=cosm%param

    !Minimum parameter values in variation
    cosm%param_min(1)=0.4
    cosm%param_min(2)=0.
    cosm%param_min(3)=1.10
    cosm%param_min(4)=1e13
    cosm%param_min(5)=0.015

    !Maximum parameter values in variation
    cosm%param_max(1)=2.
    cosm%param_max(2)=2.
    cosm%param_max(3)=1.26
    cosm%param_max(4)=1e15
    cosm%param_max(5)=0.055

    !Should the range be explored in log?
    cosm%param_log(1)=.FALSE.
    cosm%param_log(2)=.FALSE.
    cosm%param_log(3)=.FALSE.
    cosm%param_log(4)=.TRUE.
    cosm%param_log(5)=.FALSE.

    !Name halo types
    halo_type(-1)='DMONLY'
    halo_type(0)='Matter'
    halo_type(1)='CDM'
    halo_type(2)='Gas'
    halo_type(3)='Star'
    halo_type(4)='Bound gas'
    halo_type(5)='Free gas'
    halo_type(6)='Pressure'
    halo_type(7)='Void'
    halo_type(8)='Compensated void'

    !Name kernel types
    !kernel_type(1)='Lensing'
    !kernel_type(2)='Compton-y'
    !kernel_type(3)='Gravitational waves'

    !Name cross-correlation field types
    xcorr_type(1)='RCSLenS lensing'
    xcorr_type(2)='Compton y'
    xcorr_type(3)='CMB lensing'
    xcorr_type(4)='CFHTLenS lensing'
    xcorr_type(5)='KiDS lensing (z = 0.1 -> 0.9)'
    xcorr_type(6)='KiDS lensing (z = 0.1 -> 0.3)'
    xcorr_type(7)='KiDS lensing (z = 0.3 -> 0.5)'
    xcorr_type(8)='KiDS lensing (z = 0.5 -> 0.7)'
    xcorr_type(9)='KiDS lensing (z = 0.7 -> 0.9)'
    xcorr_type(10)='Gravitational waves'

  END SUBROUTINE init_HMx

  SUBROUTINE calculate_HMx(itype,k,nk,a,na,powa_lin,powa_2h,powa_1h,powa_full,cosm)

    IMPLICIT NONE
    REAL, INTENT(IN) :: k(:), a(:)
    INTEGER, INTENT(IN) :: nk, na, itype(2)
    REAL, ALLOCATABLE, INTENT(OUT) :: powa_lin(:,:), powa_2h(:,:), powa_1h(:,:), powa_full(:,:)
    TYPE(cosmology), INTENT(IN) :: cosm
    INTEGER :: i
    REAL :: z
    TYPE(tables) :: lut

    IF(ALLOCATED(powa_lin))  DEALLOCATE(powa_lin)
    IF(ALLOCATED(powa_2h))   DEALLOCATE(powa_2h)
    IF(ALLOCATED(powa_1h))   DEALLOCATE(powa_1h)
    IF(ALLOCATED(powa_full)) DEALLOCATE(powa_full)

    !Allocate power arrays
    ALLOCATE(powa_lin(nk,na),powa_2h(nk,na),powa_1h(nk,na),powa_full(nk,na))

    !Do the halo-model calculation
    IF(verbose) WRITE(*,*) 'CALCULATE_HMx: Doing calculation'
    DO i=na,1,-1
       z=redshift_a(a(i))
       CALL halomod_init(mmin,mmax,z,lut,cosm)
       IF(verbose) WRITE(*,fmt='(A5,I5,F10.2)') 'HMx:', i, REAL(z)
       CALL calculate_halomod(itype(1),itype(2),k,nk,z,powa_lin(:,i),powa_2h(:,i),powa_1h(:,i),powa_full(:,i),lut,cosm)
    END DO
    IF(verbose) THEN
       WRITE(*,*) 'CALCULATE_HMx: Done'
       WRITE(*,*)
    END IF

  END SUBROUTINE calculate_HMx

  SUBROUTINE set_ix(ix,ip)

    !Set the cross-correlation type
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: ix(2)
    INTEGER, INTENT(OUT) :: ip(2)
    INTEGER :: i, j

    DO i=1,2

       IF(ix(i)==-1) THEN
          WRITE(*,fmt='(A20,I3)') 'SET_IX: Choose field: ', i
          WRITE(*,*) '========================='
          DO j=1,SIZE(xcorr_type)
             WRITE(*,fmt='(I3,A3,A30)') j, '- ', TRIM(xcorr_type(j))
          END DO
          READ(*,*) ix(i)
          WRITE(*,*) '========================='
          WRITE(*,*)
       END IF

       IF(ix(i)==2) THEN
          !Compton y
          ip(i)=6 !Profile type: 6 - Pressure
          !ik(i)=2 !Kernel type: 2 - y
          !inz(i)=-1 !n(z) distribution - NOT used here
       ELSE IF(ix(i)==10) THEN
          !Gravitational waves
          ip(i)=-1 !Profile type: -1 DMONLY
          !ik(i)=3 !Kernel type: 3 gravity waves
          !inz(i)=-1 !n(z) distribution - NOT used here
       ELSE
          !Gravitational lensing
          !ip(i)=-1 !Profile type: -1 - DMONLY
          ip(i)=0 !Profile type: 0 - Matter
          !ik(i)=1 !Kernel type: 1 - Lensing
          !IF(ix(i)==1) inz(i)=1 !n(z): 1 - RCSLenS
          !IF(ix(i)==3) inz(i)=0 !n(z): 0 - CMB
          !IF(ix(i)==4) inz(i)=7 !n(z): 7 - CFHTLenS
          !IF(ix(i)==5) inz(i)=2 !n(z): 2 - KiDS 0.1 -> 0.9
          !IF(ix(i)==6) inz(i)=3 !n(z): 2 - KiDS 0.1 -> 0.3
          !IF(ix(i)==7) inz(i)=4 !n(z): 2 - KiDS 0.3 -> 0.5
          !IF(ix(i)==8) inz(i)=5 !n(z): 2 - KiDS 0.5 -> 0.7
          !IF(ix(i)==9) inz(i)=6 !n(z): 2 - KiDS 0.7 -> 0.9
          !ELSE
          !    STOP 'SET_IX: inx specified incorrectly'
       END IF
    END DO

  END SUBROUTINE set_ix

  SUBROUTINE xcorr(ix,ell,Cell,nl,cosm,verbose)

    !Calculates the C(l) for the cross correlation of fields ix(1) and ix(2)
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: ix(2)
    INTEGER, INTENT(IN) :: nl
    REAL, INTENT(IN) :: ell(nl)
    REAL, INTENT(OUT) :: Cell(nl)
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL, ALLOCATABLE :: k(:), pow_lin(:), pow_2h(:), pow_1h(:), pow(:,:)
    !TYPE(tables) :: lut
    !TYPE(lensing) :: lens
    TYPE(projection) :: proj(2)
    LOGICAL, INTENT(IN) :: verbose
    REAL :: kmin, kmax, amin, amax
    INTEGER :: nk, na, ip(2)
    REAL :: r1, r2

    !Set the k range
    kmin=1e-3
    kmax=1e1
    nk=32
    CALL fill_array(log(kmin),log(kmax),k,nk)
    k=exp(k)   

    !Set the a range
    amin=0.1 !scale_factor(cosm%z_cmb) !Problems with one-halo term if amin is less than 0.1
    amax=1.
    na=16
    CALL fill_array(amin,amax,a,na)

    !Allocate power arrays
    ALLOCATE(pow(nk,na),pow_lin(nk),pow_2h(nk),pow_1h(nk))

    IF(verbose) THEN
       WRITE(*,*) 'XCORR: Cross-correlation information'
       WRITE(*,*) 'XCORR: P(k) minimum k [h/Mpc]:', REAL(kmin)
       WRITE(*,*) 'XCORR: P(k) maximum k [h/Mpc]:', REAL(kmax)
       WRITE(*,*) 'XCORR: Number of k:', nk
       WRITE(*,*) 'XCORR: minimum a:', REAL(amin)
       WRITE(*,*) 'XCORR: maximum a:', REAL(amax)
       WRITE(*,*) 'XCORR: number of a:', na
       WRITE(*,*) 'XCORR: minimum ell:', REAL(lmin)
       WRITE(*,*) 'XCORR: maximum ell:', REAL(lmax)
       WRITE(*,*) 'XCORR: number of ell:', nl
       WRITE(*,*)
    END IF

    CALL set_ix(ix,ip)

    !Loop over scale factors
    !DO j=na,1,-1
    !
    !   z=-1+1./a(j)
    !
    !   !Initiliasation for the halomodel calcualtion
    !   CALL halomod_init(mmin,mmax,z,lut,cosm)
    !   CALL calculate_halomod(ip(1),ip(2),k,nk,z,pow_lin,pow_2h,pow_1h,pow(:,j),lut,cosm)
    !
    !   !Write progress to screen
    !   IF(j==na .AND. verbose) THEN
    !      WRITE(*,fmt='(A5,A7)') 'i', 'a'
    !      WRITE(*,fmt='(A13)') '   ============'
    !   END IF
    !   WRITE(*,fmt='(I5,F8.3)') j, a(j)
    !
    !END DO
    !WRITE(*,fmt='(A13)') '   ============'
    !WRITE(*,*)

    CALL calculate_HMx(ip,k,nk,a,na,powa_lin,powa_2h,powa_1h,powa_full,cosm)

    !Fill out the projection kernels
    CALL fill_projection_kernels(ix,proj,cosm)
    IF(verbose) CALL write_projection_kernels(proj,cosm)

    !Set the distance range for the Limber integral
    r1=0.
    r2=maxdist(proj)

    !Actually calculate the C(ell), but only for the full halo model part
    CALL calculate_Cell(r1,r2,ell,Cell,nl,k,a,pow,nk,na,proj,cosm)

    IF(verbose) THEN
       WRITE(*,*) 'XCORR: Done'
       WRITE(*,*)
    END IF

  END SUBROUTINE xcorr

  SUBROUTINE write_nz(lens,output)

    IMPLICIT NONE
    TYPE(lensing), INTENT(IN) :: lens
    CHARACTER(len=256), INTENT(IN) :: output
    INTEGER :: i

    OPEN(7,file=output)
    DO i=1,lens%nnz
       WRITE(7,*) lens%z_nz(i), lens%nz(i)
    END DO
    CLOSE(7)

  END SUBROUTINE write_nz

  SUBROUTINE fill_projection_kernels(ix,proj,cosm)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ix(2)
    TYPE(cosmology), INTENT(IN) :: cosm
    !TYPE(lensing), INTENT(IN) :: lens
    TYPE(projection) :: proj(2)
    INTEGER :: nk, i

    IF(ix(1)==ix(2)) THEN
       nk=1
    ELSE
       nk=2
    END IF

    !Loop over the two kernels
    DO i=1,nk
       !Fill out the projection kernels
       !Repetition is a bit ugly, but probably cannot be avoided because the size
       !of the projection X and r_X arrays will be different in general
       !IF(i==1) THEN
       !   IF(ik(i)==1) CALL fill_lensing_kernel(inz(i),proj%r_x1,proj%x1,proj%nx1,lens,cosm)
       !   IF(ik(i)==2 .OR. ik(i)==3) CALL fill_kernel(ik(i),proj%r_x1,proj%x1,proj%nx1,cosm)
       !ELSE IF(i==2) THEN
       !   IF(ik(i)==1) CALL fill_lensing_kernel(inz(i),proj%r_x2,proj%x2,proj%nx2,lens,cosm)
       !   IF(ik(i)==2 .OR. ik(i)==3) CALL fill_kernel(ik(i),proj%r_x2,proj%x2,proj%nx2,cosm)
       !END IF
       CALL fill_projection_kernel(ix(i),proj(i),cosm)
    END DO

    !In case the autospectrum is being considered
    IF(nk==1) THEN
       !proj%nx2=proj%nx1
       !IF(ALLOCATED(proj%r_x2)) DEALLOCATE(proj%r_x2)
       !IF(ALLOCATED(proj%x2))   DEALLOCATE(proj%x2)
       !ALLOCATE(proj%r_x2(proj%nx2),proj%x2(proj%nx2))
       !proj%r_x2=proj%r_x1
       !proj%x2=proj%x1
       proj(2)=proj(1)
    END IF

    !Get the maximum distance to be considered for the Limber integral
    !CALL maxdist(proj,cosm)

  END SUBROUTINE fill_projection_kernels

  SUBROUTINE fill_projection_kernel(ix,proj,cosm)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ix
    TYPE(projection), INTENT(OUT) :: proj
    TYPE(cosmology), INTENT(IN) :: cosm
    TYPE(lensing) :: lens

    IF(ix==2 .OR. ix==3) THEN
       CALL fill_kernel(ix,proj,cosm)
    ELSE
       CALL fill_lensing_kernel(ix,proj,lens,cosm)
    END IF

  END SUBROUTINE fill_projection_kernel

  SUBROUTINE calculate_Cell(r1,r2,ell,Cell,nl,k,a,pow,nk,na,proj,cosm)

    !Calculates C(l) using the Limber approximation
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nl, nk, na
    REAL, INTENT(IN) :: ell(nl)
    REAL, INTENT(OUT) :: Cell(nl)
    REAL, INTENT(IN) :: k(nk), a(na), pow(nk,na)
    REAL, INTENT(IN) :: r1, r2
    TYPE(projection), INTENT(IN) :: proj(2)
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: logk(nk), loga(na), logpow(nk,na)
    INTEGER :: i, j
    !REAL, PARAMETER :: acc=1d-4

    !Note that using Limber and flat-sky for sensible results limits lmin to ~10

    !Create log tables to speed up 2D find routine in find_pkz
    !WRITE(*,*) 'Cocker'
    logk=log(k)
    loga=log(a)
    DO j=1,na
       logpow(:,j)=log((2.*pi**2)*pow(:,j)/k**3)
    END DO
    !WRITE(*,*) 'Cocker'
    !WRITE(*,*)

    !Write some useful things to the screen
    IF(verbose) THEN
       WRITE(*,*) 'CALCULATE_CELL: ell min:', REAL(ell(1))
       WRITE(*,*) 'CALCULATE_CELL: ell max:', REAL(ell(nl))
       WRITE(*,*) 'CALCULATE_CELL: number of ell:', nl
       WRITE(*,*) 'CALCULATE_CELL: number of k:', nk
       WRITE(*,*) 'CALCULATE_CELL: number of a:', na
       WRITE(*,*) 'CALCULATE_CELL: Minimum distance [Mpc/h]:', REAL(r1)
       WRITE(*,*) 'CALCULATE_CELL: Maximum distance [Mpc/h]:', REAL(r2)
    END IF

    !Finally do the integration
    IF(verbose) WRITE(*,*) 'CALCULATE CELL: Doing calculation'
    DO i=1,nl
       !WRITE(*,*) i, ell(i)
       Cell(i)=integrate_Limber(ell(i),r1,r2,logk,loga,logpow,nk,na,acc,3,proj,cosm)
       !WRITE(*,*) i, ell(i), Cell(i)
    END DO
    IF(verbose) THEN
       WRITE(*,*) 'CALCULATE_CELL: Done'
       WRITE(*,*)
    END IF

  END SUBROUTINE calculate_Cell

  SUBROUTINE Cell_contribution(r1,r2,k,a,pow,nk,na,proj,cosm)

    !Calculates C(l) using the Limber approximation
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nk, na
    REAL, INTENT(IN) :: k(nk), a(na), pow(nk,na)
    REAL, INTENT(IN) :: r1, r2
    TYPE(projection), INTENT(IN) :: proj(2)
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: logk(nk), loga(na), logpow(nk,na)
    INTEGER :: i, j, l
    CHARACTER(len=256) :: fbase, fext, outfile

    INTEGER, PARAMETER :: n=16 !Number of ell values to take, from ell=1 to ell=2**(n-1)

    !Note that using Limber and flat-sky for sensible results limits lmin to ~10

    !Create log tables to speed up 2D find routine in find_pkz
    logk=log(k)
    loga=log(a)
    DO j=1,na
       logpow(:,j)=log((2.*pi**2)*pow(:,j)/k**3)
    END DO

    !Now call the contribution subroutine
    fbase='Limber/Cell_contrib_ell_'
    fext='.dat'
    DO i=1,n
       l=2**(i-1)
       outfile=number_file(fbase,i,fext)
       CALL Limber_contribution(REAL(l),r1,r2,logk,loga,logpow,nk,na,proj,cosm,outfile)
    END DO

  END SUBROUTINE Cell_contribution

  SUBROUTINE write_Cell(ell,Cell,nl,output)

    !Write C(ell) to a file
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nl
    REAL, INTENT(IN) :: ell(nl), Cell(nl)
    CHARACTER(len=256) :: output
    INTEGER :: i

    OPEN(7,file=output)
    DO i=1,nl
       WRITE(7,*) ell(i), Cell(i), ell(i)*(1.+ell(i))*Cell(i)/(2.*pi)
    END DO
    CLOSE(7)

  END SUBROUTINE write_Cell

  SUBROUTINE calculate_xi(th_tab,xi_tab,nth,l_tab,cl_tab,nl,lmax)

    !Calcuate the correlation functions given a C(ell) table
    IMPLICIT NONE
    REAL, INTENT(IN) :: l_tab(nl), cl_tab(nl)
    REAL, INTENT(OUT) :: th_tab(nth), xi_tab(3,nth)
    INTEGER, INTENT(IN) :: nl, lmax, nth
    INTEGER :: i, j
    REAL :: logl(nl), logCl(nl)
    REAL :: theta, Cl, l, xi0, xi2, xi4

    !Speed up find routine by doing logarithms in advance
    logl=log(l_tab)
    logCl=log(cl_tab)

    !WRITE(*,*) 'CALCULATE_XI: Computing correlation functions via sum'
    DO i=1,nth

       !Get theta value and convert from degrees to radians
       theta=th_tab(i)/rad2deg

       !Set values to zero before summing
       xi0=0.
       xi2=0.
       xi4=0.

       !Do the conversion from Cl to xi as a summation over integer ell
       DO j=1,lmax

          l=REAL(j)
          !Cl=exp(find(log(l),log(l_tab),log(Cl_tab),nl,3,3,2))
          Cl=exp(find(log(l),logl,logCl,nl,3,3,2))

          xi0=xi0+(2.*l+1.)*Cl*Bessel(0,l*theta)
          xi2=xi2+(2.*l+1.)*Cl*Bessel(2,l*theta)
          xi4=xi4+(2.*l+1.)*Cl*Bessel(4,l*theta)

       END DO

       !Divide by correct pre-factor
       xi0=xi0/(4.*pi)
       xi2=xi2/(4.*pi)
       xi4=xi4/(4.*pi)

       !Convert theta from radians to degrees
       theta=theta*rad2deg

       !Populate tables
       th_tab(i)=theta
       xi_tab(1,i)=xi0
       xi_tab(2,i)=xi2
       xi_tab(3,i)=xi4

    END DO
    !WRITE(*,*) 'CALCULATE_XI: Done'
    !WRITE(*,*)

  END SUBROUTINE calculate_xi

  SUBROUTINE write_xi(th_tab,xi_tab,nth,output)

    IMPLICIT NONE
    REAL, INTENT(IN) :: th_tab(nth), xi_tab(3,nth)
    INTEGER, INTENT(IN) :: nth
    CHARACTER(len=256), INTENT(IN) :: output

    OPEN(7,file=output)
    DO i=1,nth
       WRITE(7,*) th_tab(i), xi_tab(1,i), xi_tab(2,i), xi_tab(3,i)
    END DO
    CLOSE(7)

  END SUBROUTINE write_xi

  SUBROUTINE calculate_halomod(itype1,itype2,k,nk,z,pow_lin,pow_2h,pow_1h,pow,lut,cosm)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: itype1, itype2
    INTEGER, INTENT(IN) :: nk
    REAL, INTENT(IN) :: k(nk), z
    REAL, INTENT(OUT) :: pow_lin(nk), pow_2h(nk), pow_1h(nk), pow(nk)
    TYPE(cosmology), INTENT(IN) :: cosm
    TYPE(tables), INTENT(IN) :: lut
    INTEGER :: i

    !Write to screen
    IF(verbose) THEN
       WRITE(*,*) 'CALCULATE_HALOMOD: k min:', k(1)
       WRITE(*,*) 'CALCULATE_HALOMOD: k max:', k(nk)
       WRITE(*,*) 'CALCULATE_HALOMOD: number of k:', nk
       WRITE(*,*) 'CALCULATE_HALOMOD: z:', z
       WRITE(*,*) 'CALCULATE_HALOMOD: Calculating halo-model power spectrum'
    END IF

    !Loop over k values
    !ADD OMP support properly. What is private and shared? CHECK THIS!
!!$OMP PARALLEL DO DEFAULT(SHARED), private(k,plin, pfull,p1h,p2h)
    DO i=1,nk

       !Get the linear power
       plin=p_lin(k(i),z,cosm)
       pow_lin(i)=plin

       !Do the halo model calculation
       CALL halomod(itype1,itype2,k(i),z,pow_2h(i),pow_1h(i),pow(i),plin,lut,cosm)

    END DO
!!$OMP END PARALLEL DO

    IF(verbose) THEN
       WRITE(*,*) 'CALCULATE_HALOMOD: Done'
       WRITE(*,*)
    END IF

  END SUBROUTINE calculate_halomod

  SUBROUTINE write_power(k,pow_lin,pow_2h,pow_1h,pow,nk,output)

    IMPLICIT NONE
    CHARACTER(len=256), INTENT(IN) :: output
    INTEGER, INTENT(IN) :: nk
    REAL, INTENT(IN) :: k(nk), pow_lin(nk), pow_2h(nk), pow_1h(nk), pow(nk)
    INTEGER :: i

    IF(verbose) WRITE(*,*) 'WRITE_POWER: Writing power to ', TRIM(output)

    !Loop over k values
    OPEN(7,file=output)
    DO i=1,nk
       !Fill the tables with one- and two-halo terms as well as total
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
    CHARACTER(len=256), INTENT(IN) :: base
    INTEGER, INTENT(IN) :: nk, na
    REAL, INTENT(IN) :: k(nk), a(na), pow_lin(nk,na), pow_2h(nk,na), pow_1h(nk,na), pow_full(nk,na)
    LOGICAL, INTENT(IN) :: verbose
    REAL :: pow(nk,na)
    INTEGER :: i
    CHARACTER(len=256) :: output
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
    CHARACTER(len=256), INTENT(IN) :: output
    INTEGER, INTENT(IN) :: nk, na
    REAL, INTENT(IN) :: k(nk), a(na), pow(nk,na)
    LOGICAL, INTENT(IN) :: verbose
    INTEGER :: i, j
    !CHARACTER(len=256) :: output!_2halo, output_1halo, output_full, output_lin, output

    !output_lin=TRIM(base)//'_linear.dat'
    !output_2halo=TRIM(base)//'_2halo.dat'
    !output_1halo=TRIM(base)//'_1halo.dat'
    !output_full=TRIM(base)//'_full.dat'

    !Write out data to files
    IF(verbose) THEN
       !   WRITE(*,*) 'WRITE_POWER_A: Writing 2-halo power to ', TRIM(output_2halo)
       !   WRITE(*,*) 'WRITE_POWER_A: Writing 1-halo power to ', TRIM(output_1halo)
       !   WRITE(*,*) 'WRITE_POWER_A: Writing full power to ',   TRIM(output_full)
       WRITE(*,*) 'WRITE_POWER_A: The top row of the file contains the redshifts '
       WRITE(*,*) 'WRITE_POWER_A: The first entry in the file is hashes - #####'
       WRITE(*,*) 'WRITE_POWER_A: Subsequent rows first contain ''k''...'
       WRITE(*,*) 'WRITE_POWER_A: ...and then the halo-model power for each scale factor'
       WRITE(*,*) 'WRITE_POWER_A: Output:', TRIM(output)
    END IF

    !DO o=1,4
    !   IF(o==1) output=output_lin
    !   IF(o==2) output=output_2halo
    !   IF(o==3) output=output_1halo
    !   IF(o==4) output=output_full
    OPEN(7,file=output)
    DO i=0,nk
       IF(i==0) THEN
          WRITE(7,fmt='(A20,40F20.10)') '#####', (a(j), j=1,na)
       ELSE
          WRITE(7,fmt='(F20.10,40E20.10)') k(i), (pow(i,j), j=1,na)
       END IF
    END DO
    CLOSE(7)
    !END DO

    IF(verbose) THEN
       WRITE(*,*) 'WRITE_POWER_A: Done'
       WRITE(*,*)
    END IF

  END SUBROUTINE write_power_a

  SUBROUTINE halo_diagnostics(z,lut,cosm,dir)

    IMPLICIT NONE
    TYPE(cosmology), INTENT(IN) :: cosm
    TYPE(tables), INTENT(IN) :: lut
    REAL, INTENT(IN) :: z
    CHARACTER(len=256), INTENT(IN) :: dir
    REAL :: mass    
    CHARACTER(len=256) :: base, ext, outfile
    INTEGER :: m, m1, m2

    !REAL, PARAMETER :: mass1=1d13
    !REAL, PARAMETER :: mass2=1d14
    !REAL, PARAMETER :: mass3=1d15

    WRITE(*,*) 'HALO_DIAGNOSTICS: Outputting diagnostics'

    outfile=TRIM(dir)//'/mass_fractions.dat'
    CALL write_mass_fractions(cosm,outfile)

    m1=10
    m2=16
    DO m=m1,m2

       mass=10.**m

       base=TRIM(dir)//'/halo_profile_m'
       ext='.dat'
       outfile=number_file(base,m,ext)
       CALL write_halo_profiles(mass,z,lut,cosm,outfile)

       !base=TRIM(dir)//'/halo_profile_m'
       !ext='noh.dat'
       !outfile=number_file(base,m,ext)
       !CALL write_halo_profiles(cosm%h*mass,z,lut,cosm,outfile)

       base=TRIM(dir)//'/halo_window_m'
       ext='.dat'
       outfile=number_file(base,m,ext)
       CALL write_halo_transforms(mass,z,lut,cosm,outfile)

    END DO

    WRITE(*,*) 'HALO_DIAGNOSTICS: Done'
    WRITE(*,*)

  END SUBROUTINE halo_diagnostics

  SUBROUTINE halo_definitions(lut,dir)

    IMPLICIT NONE
    TYPE(tables), INTENT(IN) :: lut
    CHARACTER(len=256), INTENT(IN) :: dir
    CHARACTER(len=256) :: fradius, fmass, fconc

    WRITE(*,*) 'HALO_DEFINITIONS: Outputting diagnostics'

    fradius=TRIM(dir)//'/radius.dat'
    fmass=TRIM(dir)//'/mass.dat'
    fconc=TRIM(dir)//'/concentration.dat'

    OPEN(7,file=fradius)
    OPEN(8,file=fmass)
    OPEN(9,file=fconc)
    DO i=1,lut%n
       WRITE(7,*) lut%rv(i), lut%r200(i), lut%r500(i), lut%r200c(i), lut%r500c(i)
       WRITE(8,*) lut%m(i),  lut%m200(i), lut%m500(i), lut%m200c(i), lut%m500c(i)
       WRITE(9,*) lut%c(i),  lut%c200(i), lut%c500(i), lut%c200c(i), lut%c500c(i)
    END DO
    CLOSE(7)
    CLOSE(8)
    CLOSE(9)

    WRITE(*,*) 'HALO_DEFINITIONS: Done'
    WRITE(*,*)

  END SUBROUTINE halo_definitions

  SUBROUTINE write_mass_fractions(cosm,outfile)

    IMPLICIT NONE
    CHARACTER(len=256), INTENT(IN) :: outfile
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: m, mmin, mmax
    INTEGER :: i, j, n

    mmin=1e10
    mmax=1e16
    n=101

    OPEN(7,file=outfile)
    DO i=1,n
       m=exp(progression(log(mmin),log(mmax),i,n))
       WRITE(7,*) m, (halo_fraction(j,m,cosm), j=1,5)
       !halo_fraction(1,m,cosm), halo_fraction(2,m,cosm), halo_fraction(3,m,cosm), halo_fraction(4,m,cosm), halo_fraction(5,m,cosm)
    END DO
    CLOSE(7)

  END SUBROUTINE write_mass_fractions

  SUBROUTINE write_halo_profiles(m,z,lut,cosm,outfile)

    IMPLICIT NONE
    CHARACTER(len=256), INTENT(IN) :: outfile
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL, INTENT(IN) :: m, z
    REAL :: r, rv, rs, c
    INTEGER :: i
    TYPE(tables), INTENT(IN) :: lut

    REAL, PARAMETER :: rmin=1e-3 !Mininum r/rv
    REAL, PARAMETER :: rmax=1e1 !Maximum r/rv
    INTEGER, PARAMETER :: n=512 !Number of points

    !Calculate halo attributes
    rv=exp(find(log(m),log(lut%m),log(lut%rv),lut%n,3,3,2))
    c=find(log(m),log(lut%m),lut%c,lut%n,3,3,2)
    !c=2.*c !To mimic baryonic contraction, or some such bullshit
    rs=rv/c

    !Max and min r/rv and number of points
    !xmin=1d-3
    !xmax=1d1
    !n=201

    OPEN(7,file=outfile)
    DO i=1,n
       !x=exp(log(xmin)+log(xmax/xmin)*float(i-1)/float(n-1))
       r=exp(progression(log(rmin),log(rmax),i,n))
       !r=x*rv
       WRITE(7,*) r, win_type(0,1,r,m,rv,rs,z,lut,cosm), win_type(0,2,r,m,rv,rs,z,lut,cosm), win_type(0,3,r,m,rv,rs,z,lut,cosm), win_type(0,4,r,m,rv,rs,z,lut,cosm), win_type(0,5,r,m,rv,rs,z,lut,cosm), win_type(0,6,r,m,rv,rs,z,lut,cosm)
    END DO
    CLOSE(7)

  END SUBROUTINE write_halo_profiles

  SUBROUTINE write_halo_transforms(m,z,lut,cosm,outfile)

    IMPLICIT NONE
    CHARACTER(len=256), INTENT(IN) :: outfile
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL, INTENT(IN) :: m, z
    TYPE(tables), INTENT(IN) :: lut
    REAL :: x, rv, c, rs, k, rhobar
    INTEGER :: i

    REAL, PARAMETER :: xmin=1e-1 !Mininum r/rv
    REAL, PARAMETER :: xmax=1e2 !Maximum r/rv
    INTEGER, PARAMETER :: n=512 !Number of points

    !Calculate halo attributes
    rv=exp(find(log(m),log(lut%m),log(lut%rv),lut%n,3,3,2))
    c=find(log(m),log(lut%m),lut%c,lut%n,3,3,2)
    !c=2.*c !To mimic baryonic contraction, or some such bullshit
    rs=rv/c

    !Max and min k*rv and number of points
    !xmin=1d-1
    !xmax=1d2
    !n=201

    rhobar=comoving_matter_density(cosm)

    OPEN(7,file=outfile)
    DO i=1,n
       !x=exp(log(xmin)+log(xmax/xmin)*float(i-1)/float(n-1))
       x=exp(progression(log(xmin),log(xmax),i,n))
       k=x/rv
       WRITE(7,*) x, win_type(1,1,k,m,rv,rs,z,lut,cosm)*rhobar/m, win_type(1,2,k,m,rv,rs,z,lut,cosm)*rhobar/m, win_type(1,3,k,m,rv,rs,z,lut,cosm)*rhobar/m, win_type(1,4,k,m,rv,rs,z,lut,cosm)*rhobar/m, win_type(1,5,k,m,rv,rs,z,lut,cosm)*rhobar/m
    END DO
    CLOSE(7)

  END SUBROUTINE write_halo_transforms

  FUNCTION Delta_v(z,cosm)

    IMPLICIT NONE
    REAL :: Delta_v
    REAL, INTENT(IN) :: z
    TYPE(cosmology), INTENT(IN) :: cosm

    !Virialised overdensity
    IF(imead==0 .OR. imead==-1) THEN
       !Delta_v=200.
       Delta_v=Dv_brynor(z,cosm)
    ELSE IF(imead==1) THEN
       Delta_v=418.*(omega_m(z,cosm)**(-0.352))
    ELSE
       STOP 'Error, imead defined incorrectly'
    END IF

  END FUNCTION Delta_v

  FUNCTION Dv_brynor(z,cosm)

    !Bryan & Norman (1998) spherical over-density calculation
    IMPLICIT NONE
    REAL :: Dv_brynor
    REAL :: x, om_m
    REAL, INTENT(IN) :: z
    TYPE(cosmology), INTENT(IN) :: cosm

    om_m=omega_m(z,cosm)
    x=om_m-1.

    IF(cosm%om_v==0.) THEN
       STOP 'Dv_BRYNOR: Should not be in here'
       !Open model results
       Dv_brynor=177.65+60.*x-32.*x**2
       Dv_brynor=dv_brynor/om_m
    ELSE
       !LCDM results
       Dv_brynor=177.65+82.*x-39.*x**2
       Dv_brynor=dv_brynor/om_m
    END IF

  END FUNCTION Dv_brynor

  FUNCTION delta_c(z,cosm)

    IMPLICIT NONE
    REAL :: delta_c
    REAL, INTENT(IN) :: z
    TYPE(cosmology), INTENT(IN) :: cosm

    !Linear collapse density
    IF(imead==0 .OR. imead==-1) THEN
       !Nakamura & Suto (1997) fitting formula for LCDM
       delta_c=1.686*(1.+0.0123*log10(omega_m(z,cosm)))
    ELSE IF(imead==1) THEN
       delta_c=1.59+0.0314*log(sigma_cb(8.,z,cosm))
       delta_c=delta_c*(1.+0.0123*log10(omega_m(z,cosm)))
    ELSE
       STOP 'Error, imead defined incorrectly'
    END IF

  END FUNCTION delta_c

  FUNCTION eta(z,cosm)

    IMPLICIT NONE
    REAL :: eta
    REAL, INTENT(IN) :: z
    TYPE(cosmology), INTENT(IN) :: cosm

    IF(imead==0 .OR. imead==-1) THEN
       eta=0.
    ELSE IF(imead==1) THEN
       !The first parameter here is 'eta_0' in Mead et al. (2015; arXiv 1505.07833)
       eta=0.603-0.3*(sigma_cb(8.,z,cosm))
    ELSE
       STOP 'Error, imead defined incorrectly'
    END IF

  END FUNCTION eta

  FUNCTION kstar(lut,cosm)

    IMPLICIT NONE
    REAL :: kstar
    TYPE(cosmology), INTENT(IN) :: cosm
    TYPE(tables), INTENT(IN) :: lut
    REAL :: crap

    !To prevent compile-time warnings
    crap=cosm%A

    IF(imead==0 .OR. imead==-1) THEN
       !Set to zero for the standard Poisson one-halo term
       kstar=0.
    ELSE IF(imead==1) THEN
       !One-halo cut-off wavenumber
       kstar=0.584*(lut%sigv)**(-1)
    ELSE
       STOP 'Error, imead defined incorrectly'
    END IF

  END FUNCTION kstar

  FUNCTION As(cosm)

    IMPLICIT NONE
    REAL :: As
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: crap

    !To prevent compile-time warnings
    crap=cosm%A

    !Halo concentration pre-factor
    IF(imead==0 .OR. imead==-1) THEN
       !Set to 4 for the standard Bullock value
       As=4.
    ELSE IF(imead==1) THEN
       !This is the 'A' halo-concentration parameter in Mead et al. (2015; arXiv 1505.07833)
       As=3.13
    ELSE
       STOP 'Error, imead defined incorrectly'
    END IF

  END FUNCTION As

  FUNCTION fdamp(z,cosm)

    IMPLICIT NONE
    REAL ::fdamp
    REAL, INTENT(IN) :: z
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: crap

    !To prevent compile-time warnings
    crap=cosm%A
    crap=z

    !Linear theory damping factor
    IF(imead==0 .OR. imead==-1) THEN
       !Set to 0 for the standard linear theory two halo term
       fdamp=0.
    ELSE IF(imead==1) THEN
       !fdamp=0.188*sigma_cb(8.,z,cosm)**4.29
       fdamp=0.0095*lut%sigv100**1.37
       !Catches extreme values of fdamp that occur for ridiculous cosmologies
       IF(fdamp<1.e-3) fdamp=0.
       IF(fdamp>0.99)  fdamp=0.99
    ELSE
       STOP 'Error, imead defined incorrectly'
    END IF

  END FUNCTION fdamp

  FUNCTION alpha_transition(lut,cosm)

    IMPLICIT NONE
    REAL :: alpha_transition
    TYPE(tables), INTENT(IN) :: lut
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: crap

    !To prevent compile-time warnings
    crap=cosm%A

    IF(imead==0 .OR. imead==-1) THEN
       !Set to 1 for the standard halo model addition of one- and two-halo terms
       alpha_transition=1.
    ELSE IF(imead==1) THEN
       !This uses the top-hat defined neff
       alpha_transition=3.24*1.85**lut%neff
    ELSE
       STOP 'Error, imead defined incorrectly'
    END IF

    !Catches values of alpha that are crazy
    IF(alpha_transition>2.)  alpha_transition=2.
    IF(alpha_transition<0.5) alpha_transition=0.5

  END FUNCTION alpha_transition

  SUBROUTINE print_halomodel_parameters(z,lut,cosm)

    IMPLICIT NONE
    REAL, INTENT(IN) :: z
    TYPE(cosmology), INTENT(IN) :: cosm
    TYPE(tables), INTENT(IN) :: lut

    !This subroutine writes out the physical parameters at some redshift 
    !(e.g. Delta_v) rather than the model parameters

    IF(verbose) THEN
       WRITE(*,*) 'PRINT_HALOMODEL_PARAMETERS: Writing out halo-model parameters'
       WRITE(*,*) 'PRINT_HALOMODEL_PARAMETERS: Halo-model parameters at your redshift'
       WRITE(*,*) '==========================='
       WRITE(*,fmt='(A10,F10.5)') 'z:', z
       WRITE(*,fmt='(A10,F10.5)') 'Dv:', Delta_v(z,cosm)
       WRITE(*,fmt='(A10,F10.5)') 'dc:', delta_c(z,cosm)
       WRITE(*,fmt='(A10,F10.5)') 'eta:', eta(z,cosm)
       WRITE(*,fmt='(A10,F10.5)') 'k*:', kstar(lut,cosm)
       WRITE(*,fmt='(A10,F10.5)') 'A:', As(cosm)
       WRITE(*,fmt='(A10,F10.5)') 'fdamp:', fdamp(z,cosm)
       WRITE(*,fmt='(A10,F10.5)') 'alpha:', alpha_transition(lut,cosm)
       WRITE(*,*) '==========================='
       WRITE(*,*) 'PRINT_HALOMODEL_PARAMETERS: Done'
       WRITE(*,*)
    END IF

  END SUBROUTINE print_halomodel_parameters

  FUNCTION r_nl(lut)

    !Calculates k_nl as 1/R where nu(R)=1.
    TYPE(tables), INTENT(IN) :: lut
    REAL :: r_nl  

    IF(lut%nu(1)>1.) THEN
       !This catches some very strange values
       r_nl=lut%rr(1)
    ELSE
       r_nl=exp(find(log(1.),log(lut%nu),log(lut%rr),lut%n,3,3,2))
    END IF

  END FUNCTION r_nl

  SUBROUTINE print_cosmology(cosm)

    IMPLICIT NONE
    TYPE(cosmology) :: cosm

    IF(verbose) THEN
       WRITE(*,*) 'COSMOLOGY: ', TRIM(cosm%name)
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'Omega_m:', cosm%om_m
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'Omega_b:', cosm%om_b
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'Omega_c:', cosm%om_c
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'Omega_v:', cosm%om_v
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'h:', cosm%h
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'w_0:', cosm%w
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'w_a:', cosm%wa
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'sig8:', cosm%sig8
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'n:', cosm%n
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'Omega:', cosm%om
       ! WRITE(*,fmt='(A11,A15,F10.5)') 'COSMOLOGY:', 'k / (Mpc/h)^-2:', cosm%k
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'z_CMB:', cosm%z_cmb
       WRITE(*,*)
    END IF

  END SUBROUTINE print_cosmology

  SUBROUTINE assign_cosmology(icosmo,cosm)

    IMPLICIT NONE
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER, INTENT(INOUT) :: icosmo
    CHARACTER(len=256) :: names(0:2)

    names(0)='Boring'
    names(1)='WMAP9 (OWLS version)'
    names(2)='Planck 2013 (OWLS version)'

    IF(icosmo==-1) THEN
       WRITE(*,*) 'ASSIGN_COSMOLOGY: Choose cosmological model'
       WRITE(*,*) '==========================================='
       !DO i=0,SIZE(names)-1
       !   WRITE(*,*) i, '- ', TRIM(names(i))
       !END DO
       WRITE(*,*) '0 - Boring'
       WRITE(*,*) '1 - WMAP9 (OWLS version)'
       WRITE(*,*) '2 - Planck 2013 (OWLS version)'
       READ(*,*) icosmo
       WRITE(*,*) '==========================================='
    END IF

    cosm%name=names(icosmo)

    !Boring defaults
    cosm%om_m=0.3
    cosm%om_b=0.05
    cosm%om_v=1.-cosm%om_m
    cosm%om_nu=0.
    cosm%h=0.7
    cosm%sig8=0.8
    cosm%n=0.96
    cosm%w=-1.
    cosm%wa=0.
    cosm%T_cmb=2.72
    cosm%z_cmb=1100.

    IF(icosmo==0) THEN
       !Boring - do nothing
    ELSE IF(icosmo==1) THEN
       !OWLS - WMAP9
       cosm%om_m=0.272
       cosm%om_b=0.0455
       cosm%om_v=1.-cosm%om_m
       cosm%om_nu=0.
       cosm%h=0.704
       cosm%sig8=0.81
       !cosm%sig8=0.797
       !cosm%sig8=0.823
       cosm%n=0.967
    ELSE IF(icosmo==2) THEN
       !OWLS - Planck 2013
       cosm%om_m=0.3175
       cosm%om_b=0.0490
       cosm%om_v=1.-cosm%om_m
       cosm%h=0.6711
       cosm%n=0.9624
       cosm%sig8=0.834
    ELSE
       STOP 'ASSIGN_COSMOLOGY: Error, icosmo not specified correctly'
    END IF

    WRITE(*,*) 'ASSIGN_COSMOLOGY: Cosmology assigned'
    WRITE(*,*)

  END SUBROUTINE assign_cosmology

  SUBROUTINE initialise_cosmology(cosm)

    IMPLICIT NONE
    REAL :: sigi
    TYPE(cosmology) :: cosm

    !Derived cosmological parameters
    cosm%om_r=2.5e-5*(1.+0.227*neff)/cosm%h**2
    cosm%om_m=cosm%om_m-cosm%om_r !Maintain flatness
    cosm%om_c=cosm%om_m-cosm%om_b-cosm%om_nu
    cosm%om=cosm%om_m+cosm%om_v
    cosm%k=(cosm%om-1.)/(conH0**2)

    IF(verbose) THEN
       WRITE(*,*) 'INITIALISE_COSMOLOGY: Omega_r:', REAL(cosm%om_r)
       WRITE(*,*) 'INITIALISE_COSMOLOGY: Omega_c:', REAL(cosm%om_c)
       WRITE(*,*) 'INITIALISE_COSMOLOGY: Omega:', REAL(cosm%om)
       WRITE(*,*) 'INITIALISE_COSMOLOGY: k:', REAL(cosm%k)
       WRITE(*,*) !Need this space because 'growth' is called in the middle
    END IF

    !Fill the tables of g(z)
    CALL fill_growtab(cosm)

    !Set the normalisation to 1 initially
    cosm%A=1.

    !Calculate the initial sigma_8 value (will not be correct)
    sigi=sigma(8.,zero,cosm)

    IF(verbose) WRITE(*,*) 'INITIALISE_COSMOLOGY: Initial sigma_8:', REAL(sigi)

    !Reset the normalisation to give the correct sigma8
    cosm%A=cosm%sig8/sigi
    !cosm%A=391.0112 !Appropriate for sig8=0.8 in the boring model (for tests)

    !Recalculate sigma8, should be correct this time
    sigi=sigma(8.,zero,cosm)

    !Write to screen
    IF(verbose) THEN
       WRITE(*,*) 'INITIALISE_COSMOLOGY: Normalisation factor:', REAL(cosm%A)
       WRITE(*,*) 'INITIALISE_COSMOLOGY: Target sigma_8:', REAL(cosm%sig8)
       WRITE(*,*) 'INITIALISE_COSMOLOGY: Final sigma_8 (calculated):', REAL(sigi)
       WRITE(*,*) 'INITIALISE_COSMOLOGY: Complete'
       WRITE(*,*)
    END IF

    !Fill tables of r vs. sigma(r)
    CALL fill_sigtab(cosm)

  END SUBROUTINE initialise_cosmology

  SUBROUTINE initialise_distances(cosm)

    !Fill up tables of a vs. r(a) (comoving distance)
    IMPLICIT NONE
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: zmin, zmax, amin, amax
    REAL :: rh
    INTEGER :: i

    INTEGER, PARAMETER :: nr=128

    zmin=0.
    zmax=cosm%z_cmb
    amin=scale_factor_z(zmax)
    amax=scale_factor_z(zmin)
    IF(verbose) THEN
       WRITE(*,*) 'INITIALISE_DISTANCE: Redshift range for r(z) tables'
       WRITE(*,*) 'INITIALISE_DISTANCE: minimum z:', REAL(zmin)
       WRITE(*,*) 'INITIALISE_DISTANCE: maximum z:', REAL(zmax)
       WRITE(*,*) 'INITIALISE_DISTANCE: minimum a:', REAL(amin)
       WRITE(*,*) 'INITIALISE_DISTANCE: maximum a:', REAL(amax)
    END IF
    cosm%nr=nr
    CALL fill_array(amin,amax,cosm%a_r,cosm%nr)
    IF(ALLOCATED(cosm%r)) DEALLOCATE(cosm%r)
    ALLOCATE(cosm%r(cosm%nr))

    !Now do the r(z) calculation
    DO i=1,cosm%nr
       cosm%r(i)=integrate_distance(cosm%a_r(i),1.,acc,3,cosm)
    END DO
    IF(verbose) THEN
       WRITE(*,*) 'INITIALISE_DISTANCE: minimum r [Mpc/h]:', REAL(cosm%r(cosm%nr))
       WRITE(*,*) 'INITIALISE_DISTANCE: maximum r [Mpc/h]:', REAL(cosm%r(1))
    END IF

    !Find the horizon distance in your cosmology
    rh=integrate_distance(0.,1.,acc,3,cosm)
    IF(verbose) THEN
       WRITE(*,*) 'INITIALISE_DISTANCE: Horizon distance [Mpc/h]:', REAL(rh)
       WRITE(*,*) 'INITIALISE_DISTANCE: Done'
       WRITE(*,*)
    END IF

  END SUBROUTINE initialise_distances

  SUBROUTINE write_distances(cosm)

    !Write file of z vs. r(z)
    IMPLICIT NONE
    TYPE(cosmology), INTENT(IN) :: cosm
    CHARACTER(len=256) :: output
    INTEGER :: i
    REAL :: z

    !Now do the r(z) calculation
    output='projection/distance.dat'
    WRITE(*,*) 'WRITE_DISTANCE: Writing r(a): ', TRIM(output)
    OPEN(7,file=output)
    DO i=1,cosm%nr
       z=redshift_a(cosm%a_r(i))
       WRITE(7,*) z, cosm%r(i), f_k(cosm%r(i),cosm)
    END DO
    CLOSE(7)
    WRITE(*,*) 'WRITE_DISTANCE: Done'
    WRITE(*,*)

  END SUBROUTINE write_distances

  FUNCTION redshift_a(a)

    IMPLICIT NONE
    REAL :: redshift_a
    REAL, INTENT(IN) :: a

    IF(a==0.) STOP 'REDSHIFT_A: Error, routine called with a=0'

    redshift_a=-1.+1./a

  END FUNCTION redshift_a

  FUNCTION scale_factor_z(z)

    IMPLICIT NONE
    REAL :: scale_factor_z
    REAL, INTENT(IN) :: z

    scale_factor_z=1./(1.+z)

  END FUNCTION scale_factor_z

  FUNCTION cosmic_distance(z,cosm)

    IMPLICIT NONE
    REAL :: cosmic_distance
    REAL, INTENT(IN) :: z
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: a

    a=scale_factor_z(z)

    cosmic_distance=find(a,cosm%a_r,cosm%r,cosm%nr,3,3,2)

  END FUNCTION cosmic_distance

  FUNCTION redshift_r(r,cosm)

    IMPLICIT NONE
    REAL :: redshift_r
    REAL, INTENT(IN) :: r
    TYPE(cosmology), INTENT(IN) :: cosm

    redshift_r=redshift_a(find(r,cosm%r,cosm%a_r,cosm%nr,3,3,2))

  END FUNCTION redshift_r

  FUNCTION maxdist(proj)!,cosm)

    !Calculates the maximum distance necessary for the lensing integration
    IMPLICIT NONE
    REAL :: maxdist
    TYPE(projection), INTENT(IN) :: proj(2)
    !TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: rmax1, rmax2

    REAL, PARAMETER :: dr=0.01

    !Fix the maximum redshift and distance (which may fixed at source plane)
    rmax1=MAXVAL(proj(1)%r_x)
    rmax2=MAXVAL(proj(2)%r_x)
    
    !Subtract a small distance here because of rounding errors in recalculating zmax
    maxdist=MIN(rmax1,rmax2)-dr
    
  END FUNCTION maxdist

  SUBROUTINE write_projection_kernels(proj,cosm)

    IMPLICIT NONE
    TYPE(projection), INTENT(IN) :: proj(2)
    TYPE(cosmology), INTENT(IN) :: cosm
    CHARACTER(len=256) :: output

    output=TRIM('projection/kernel1.dat')
    CALL write_projection_kernel(proj(1),cosm,output)

    output=TRIM('projection/kernel2.dat')
    CALL write_projection_kernel(proj(2),cosm,output)

  END SUBROUTINE write_projection_kernels

  SUBROUTINE write_projection_kernel(proj,cosm,output)

    IMPLICIT NONE
    TYPE(projection), INTENT(IN) :: proj
    TYPE(cosmology), INTENT(IN) :: cosm
    CHARACTER(len=256), INTENT(IN) :: output
    INTEGER :: i
    REAL :: r, z

    !Kernel 1
    IF(verbose) WRITE(*,*) 'WRITE_PROJECTION_KERNEL: Writing out kernel: ', TRIM(output)
    OPEN(7,file=output)
    DO i=1,proj%nX    
       r=proj%r_X(i)
       z=redshift_r(r,cosm)
       WRITE(7,*) r, z, proj%X(i)
    END DO
    CLOSE(7)
    IF(verbose) THEN
       WRITE(*,*) 'WRITE_PROJECTION_KERNEL: Writing done'
       WRITE(*,*)
    END IF

  END SUBROUTINE write_projection_kernel

  SUBROUTINE fill_lensing_kernel(ix,proj,lens,cosm)

    !Fill the lensing projection kernel
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ix
    TYPE(lensing) :: lens
    TYPE(cosmology), INTENT(IN) :: cosm
    TYPE(projection), INTENT(OUT) :: proj
    !REAL, ALLOCATABLE, INTENT(OUT) :: r_x(:), x(:)
    !INTEGER, INTENT(OUT) :: nx_out
    REAL :: zmin, zmax, rmax, r
    CHARACTER(len=256) :: output
    INTEGER :: i

    !Parameters
    REAL, PARAMETER :: rmin=0. !Minimum distance in integral    
    INTEGER, PARAMETER :: nx=128 !Number of entries in X(r) table

    !IF(inz==-1) THEN
    !   WRITE(*,*) 'FILL_LENSING_KERNEL: Choose n(z)'
    !   WRITE(*,*) '================================'
    !   WRITE(*,*) '0 - Fixed source plane'
    !   WRITE(*,*) '1 - Realistic n(z) distribution'
    !   READ(*,*) inz
    !   WRITE(*,*) '================================'
    !   WRITE(*,*)
    !END IF

    IF(ix==2 .OR. ix==10) STOP 'FILL_LENSING_KERNEL: Error, trying to do this for a non-lensing ix'

    !Choose either n(z) or fixed z_s
    IF(ix==3) THEN      
       zmin=0.
       zmax=cosm%z_cmb
       WRITE(*,*) 'FILL_LENSING_KERNEL: Source plane redshift:', REAL(zmax)
    ELSE
       CALL get_nz(ix,lens)
       zmin=lens%z_nz(1)
       zmax=lens%z_nz(lens%nnz)
       output='lensing/nz.dat'
       CALL write_nz(lens,output)
    END IF

    !Get the distance range for the lensing kernel
    rmax=cosmic_distance(zmax,cosm)
    WRITE(*,*) 'FILL_LENSING_KERNEL: minimum r [Mpc/h]:', REAL(rmin)
    WRITE(*,*) 'FILL_LENSING_KERNEL: maximum r [Mpc/h]:', REAL(rmax)
    WRITE(*,*) 'FILL_LENSING_KERNEL: minimum z:', REAL(zmin)
    WRITE(*,*) 'FILL_LENSING_KERNEL: maximum z:', REAL(zmax)
    WRITE(*,*)

    !Fill the q(r) table
    CALL fill_lensing_efficiency(ix,rmin,rmax,zmax,lens,cosm)

    !Write the q(r) table to file
    output='lensing/efficiency.dat'
    CALL write_lensing_efficiency(lens,cosm,output)

    !Assign arrays for the projection function
    proj%nx=nx
    CALL fill_array(rmin,rmax,proj%r_x,nx)
    WRITE(*,*) 'FILL_LENSING_KERNEL: number of points:', nx
    IF(ALLOCATED(proj%x)) DEALLOCATE(proj%x)
    ALLOCATE(proj%x(proj%nx))

    DO i=1,nx
       !Get r and fill X(r)
       r=proj%r_x(i)
       proj%x(i)=lensing_kernel(r,lens,cosm)  
       !Enforce that the kernel must not be negative (is this necessary?)
       IF(proj%x(i)<0.) proj%x(i)=0.
    END DO
    WRITE(*,*) 'FILL_LENSING_KERNEL: Done'
    WRITE(*,*)

  END SUBROUTINE fill_lensing_kernel

  SUBROUTINE fill_lensing_efficiency(ix,rmin,rmax,zmax,lens,cosm)

    !Fill the table for lensing efficiency
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ix
    REAL, INTENT(IN) :: rmin, rmax, zmax
    TYPE(lensing), INTENT(INOUT) :: lens
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: r, z
    INTEGER :: i

    INTEGER, PARAMETER :: nq=128 !Number of entries in q(r) table

    !Fill the r vs. q(r) tables
    lens%nq=nq
    WRITE(*,*) 'FILL_EFFICIENCY: number of points:', lens%nq
    CALL fill_array(rmin,rmax,lens%r_q,lens%nq)
    IF(ALLOCATED(lens%q)) DEALLOCATE(lens%q)
    ALLOCATE(lens%q(lens%nq))

    DO i=1,lens%nq
       r=lens%r_q(i)
       z=redshift_r(r,cosm)
       IF(r==0.) THEN
          !To avoid division by zero
          lens%q(i)=1.
       ELSE
          IF(ix==7) THEN
             !q(r) for a fixed source plane
             lens%q(i)=f_k(rmax-r,cosm)/f_k(rmax,cosm)
          ELSE
             !q(r) for a n(z) distribution 
             lens%q(i)=integrate_q(r,z,zmax,acc,3,lens,cosm)
          END IF
       END IF
    END DO
    WRITE(*,*) 'FILL_EFFICIENCY: Done writing'
    WRITE(*,*)

  END SUBROUTINE fill_lensing_efficiency

  FUNCTION q_r(r,lens)

    !Interpolation function for q(r)
    IMPLICIT NONE
    REAL :: q_r
    REAL, INTENT(IN) :: r
    TYPE(lensing), INTENT(IN) :: lens

    q_r=find(r,lens%r_q,lens%q,lens%nq,3,3,2)

  END FUNCTION q_r

  SUBROUTINE write_lensing_efficiency(lens,cosm,output)

    !Write lensing efficiency q(r) function to a file
    IMPLICIT NONE
    TYPE(lensing), INTENT(INOUT) :: lens
    TYPE(cosmology), INTENT(IN) :: cosm
    CHARACTER(len=256), INTENT(IN) :: output
    REAL :: r, z, q
    INTEGER :: i

    WRITE(*,*) 'WRITE_EFFICIENCY: Writing q(r): ', TRIM(output)
    OPEN(7,file=output)
    DO i=1,lens%nq
       r=lens%r_q(i)
       z=redshift_r(r,cosm)
       q=lens%q(i)
       WRITE(7,*) r, z, q
    END DO
    CLOSE(7)
    WRITE(*,*) 'WRITE_EFFICIENCY: Done'
    WRITE(*,*)

  END SUBROUTINE write_lensing_efficiency

  FUNCTION lensing_kernel(r,lens,cosm)

    !The lensing projection kernel
    IMPLICIT NONE
    REAL :: lensing_kernel
    REAL, INTENT(IN) :: r
    TYPE(lensing) :: lens
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: z, q

    !Get z(r)
    z=redshift_r(r,cosm)

    !Get the lensing efficiency
    q=q_r(r,lens)

    !This is then the projection kernel (X_kappa)
    lensing_kernel=(1.+z)*f_k(r,cosm)*q
    lensing_kernel=lensing_kernel*1.5*cosm%om_m/(conH0**2)

  END FUNCTION lensing_kernel

  SUBROUTINE fill_kernel(ix,proj,cosm)

    !Fill the Compton-y look-up table
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ix
    !REAL, ALLOCATABLE, INTENT(OUT) :: r_x(:), x(:)
    !INTEGER, INTENT(OUT) :: nx
    TYPE(projection), INTENT(OUT) :: proj
    TYPE(cosmology), INTENT(IN) :: cosm
    INTEGER :: i
    REAL :: rmax, r

    REAL, PARAMETER :: rmin=0. !Minimum r for table
    INTEGER, PARAMETER :: nx=128 !Entires in look-up table

    proj%nx=nx 

    !Get the distance range for the projection function
    !Use the same as that for the distance calculation
    !Assign arrays for the kernel function
    rmax=MAXVAL(cosm%r)
    CALL fill_array(rmin,rmax,proj%r_X,proj%nX)
    WRITE(*,*) 'FILL_KERNEL: minimum r [Mpc/h]:', REAL(rmin)
    WRITE(*,*) 'FILL_KERNEL: maximum r [Mpc/h]:', REAL(rmax)
    WRITE(*,*) 'FILL_KERNEL: number of points:', nx

    IF(ALLOCATED(proj%X)) DEALLOCATE(proj%X)
    ALLOCATE(proj%X(nX))

    !Now fill the y-kernel (which is simply proportional to 'a')
    DO i=1,nX
       r=proj%r_x(i)
       IF(ix==2) THEN
          proj%x(i)=y_kernel(r,cosm)
       ELSE IF(ix==10) THEN
          proj%x(i)=gravity_kernel(r,cosm)
       ELSE
          STOP 'FILL_KERNEL: Error, ix not specified correctly'
       END IF
    END DO

    WRITE(*,*) 'FILL_KERNEL: Done:'
    WRITE(*,*)

  END SUBROUTINE fill_kernel

  FUNCTION y_kernel(r,cosm)

    !The Compton-y projeciton kernel
    IMPLICIT NONE
    REAL :: y_kernel
    REAL, INTENT(IN) :: r
    TYPE(cosmology), INTENT(IN) :: cosm
    !REAL :: z, a
    REAL :: crap

    !z=redshift_r(r,cosm)
    !a=scale_factor_z(z)

    !To stop compile-time warnings
    crap=cosm%om_m
    crap=r 

    y_kernel=yfac*mpc !Convert some units; note that there is no factor of 'a'
    y_kernel=y_kernel*eV*cm**(-3) !Convert from eV cm^-3 to J m^-3

  END FUNCTION y_kernel

  FUNCTION gravity_kernel(r,cosm)

    !Projection kernel for gravitational waves (~ 1/r)
    IMPLICIT NONE
    REAL :: gravity_kernel
    REAL, INTENT(IN) :: r
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: crap

    REAL, PARAMETER :: A=1.
    REAL, PARAMETER :: rmin=10.

    !To stop compile-time warnings
    crap=cosm%om_m

    IF(r<rmin) THEN
       gravity_kernel=A/rmin
    ELSE
       gravity_kernel=A/r
    END IF

  END FUNCTION gravity_kernel

  FUNCTION f_k(r,cosm)

    !Curvature function
    IMPLICIT NONE
    REAL :: f_k
    REAL, INTENT(IN) :: r
    TYPE(cosmology), INTENT(IN) :: cosm

    IF(cosm%k==0.) THEN
       f_k=r
    ELSE IF(cosm%k<0.) THEN
       f_k=sinh(sqrt(-cosm%k)*r)/sqrt(-cosm%k)
    ELSE IF(cosm%k>0.) THEN
       f_k=sin(sqrt(cosm%k)*r)/sqrt(cosm%k)
    ELSE
       STOP 'F_K: Something went wrong'
    END IF

  END FUNCTION f_k

  FUNCTION fdash_k(r,cosm)

    !Derivative of curvature function
    IMPLICIT NONE
    REAL :: fdash_k
    REAL, INTENT(IN) :: r
    TYPE(cosmology), INTENT(IN) :: cosm

    IF(cosm%k==0.) THEN
       fdash_k=1.
    ELSE IF(cosm%k<0.) THEN
       fdash_k=cosh(sqrt(-cosm%k)*r)
    ELSE IF(cosm%k>0.) THEN
       fdash_k=cos(sqrt(cosm%k)*r)
    ELSE
       STOP 'F_K: Something went wrong'
    END IF

  END FUNCTION fdash_k

  SUBROUTINE random_cosmology(cosm)

    IMPLICIT NONE
    TYPE(cosmology) :: cosm
    REAL :: om_m_min, om_m_max, om_b_min, om_b_max, n_min, n_max
    REAL :: w_min, w_max, h_min, h_max, sig8_min, sig8_max, wa_min, wa_max

    !Needs to be set to normalise P_lin
    cosm%A=1.

    om_m_min=0.1
    om_m_max=1.
    cosm%om_m=uniform(om_m_min,om_m_max)

    cosm%om_v=1.-cosm%om_m

    om_b_min=0.005
    om_b_max=MIN(0.095,cosm%om_m)
    cosm%om_b=uniform(om_b_min,om_b_max)

    cosm%om_c=cosm%om_m-cosm%om_b

    n_min=0.5
    n_max=1.5
    cosm%n=uniform(n_min,n_max)

    h_min=0.4
    h_max=1.2
    cosm%h=uniform(h_min,h_max)

    w_min=-1.5
    w_max=-0.5
    cosm%w=uniform(w_min,w_max)

    wa_min=-1.
    wa_max=-cosm%w*0.8
    cosm%wa=uniform(wa_min,wa_max)

    sig8_min=0.2
    sig8_max=1.5
    cosm%sig8=uniform(sig8_min,sig8_max)

  END SUBROUTINE random_cosmology

  SUBROUTINE allocate_LUT(lut,n)

    IMPLICIT NONE
    TYPE(tables) :: lut
    INTEGER, INTENT(IN) :: n

    !Allocates memory for the look-up tables
    lut%n=n

    ALLOCATE(lut%zc(n),lut%m(n),lut%c(n),lut%rv(n))
    ALLOCATE(lut%nu(n),lut%rr(n),lut%sigf(n),lut%sig(n))
    ALLOCATE(lut%m500(n),lut%r500(n),lut%c500(n))
    ALLOCATE(lut%m500c(n),lut%r500c(n),lut%c500c(n))
    ALLOCATE(lut%m200(n),lut%r200(n),lut%c200(n))
    ALLOCATE(lut%m200c(n),lut%r200c(n),lut%c200c(n))

    !Experimental window look-up table
    !lut%nk=nk
    !ALLOCATE(lut%log_m(n),lut%log_k(nk),lut%log_win(n,nk))
    !lut%log_k=0.
    !lut%log_win=0.
    !lut%iwin=.FALSE.

    lut%zc=0.
    lut%m=0.
    lut%c=0.
    lut%rv=0.
    lut%nu=0.
    lut%rr=0.
    lut%sigf=0.
    lut%sig=0.

    lut%m500=0.
    lut%r500=0.
    lut%c500=0.

    lut%m500c=0.
    lut%r500c=0.
    lut%c500c=0.

    lut%m200=0.
    lut%r200=0.
    lut%c200=0.

    lut%m200c=0.
    lut%r200c=0.
    lut%c200c=0.

    !Experimental log tables
    ALLOCATE(lut%log_m(n))
    lut%log_m=0.

  END SUBROUTINE allocate_LUT

  SUBROUTINE deallocate_LUT(lut)

    IMPLICIT NONE
    TYPE(tables) :: lut

    !Deallocates look-up tables
    DEALLOCATE(lut%zc,lut%m,lut%c,lut%rv,lut%nu,lut%rr,lut%sigf,lut%sig)
    DEALLOCATE(lut%m500,lut%r500,lut%c500,lut%m500c,lut%r500c,lut%c500c)
    DEALLOCATE(lut%m200,lut%r200,lut%c200,lut%m200c,lut%r200c,lut%c200c)

    !Deallocate experimental window tables
    !DEALLOCATE(lut%log_win,lut%log_k)

    !Deallocate experimental log tables
    DEALLOCATE(lut%log_m)

  END SUBROUTINE deallocate_LUT

  SUBROUTINE halomod_init(mmin,mmax,z,lut,cosm)

    IMPLICIT NONE
    REAL, INTENT(IN) :: z
    REAL, INTENT(IN) :: mmin, mmax
    INTEGER :: i
    REAL :: Dv, dc, f, m, nu, r, sig, A0, rhom, rhoc
    TYPE(cosmology) :: cosm
    TYPE(tables) :: lut

    INTEGER, PARAMETER :: n=64 !Number of mass entries in look-up table

    !Halo-model initialisation routine
    !The computes other tables necessary for the one-halo integral

    !Find value of sigma_v
    lut%sigv=sqrt(dispint(0.,z,cosm)/3.)
    lut%sigv100=sqrt(dispint(100.,z,cosm)/3.)
    lut%sig8z=sigma(8.,z,cosm)

    IF(verbose) THEN
       WRITE(*,*) 'HALOMOD_INIT: Filling look-up tables'
       WRITE(*,*) 'HALOMOD_INIT: tables being filled at redshift:', REAL(z)
       WRITE(*,*) 'HALOMOD_INIT: sigv [Mpc/h]:', REAL(lut%sigv)
       WRITE(*,*) 'HALOMOD_INIT: sigv100 [Mpc/h]:', REAL(lut%sigv100)
       WRITE(*,*) 'HALOMOD_INIT: sig8(z):', REAL(lut%sig8z)
    END IF

    IF(ALLOCATED(lut%rr)) CALL deallocate_LUT(lut)

    CALL allocate_LUT(lut,n)

    dc=delta_c(z,cosm)

    DO i=1,n

       !m=exp(log(mmin)+log(mmax/mmin)*float(i-1)/float(n-1))
       m=exp(progression(log(mmin),log(mmax),i,n))
       r=radius_m(m,cosm)
       sig=sigma_cb(r,z,cosm)
       nu=dc/sig

       lut%m(i)=m
       lut%rr(i)=r
       lut%sig(i)=sig
       lut%nu(i)=nu

    END DO

    IF(verbose) WRITE(*,*) 'HALOMOD_INIT: m, r, nu, sig tables filled'

    !Fills up a table for sigma(fM) for Bullock c(m) relation
    !This is the f=0.01 parameter in the Bullock realtion sigma(fM,z)
    f=0.01**(1./3.)
    DO i=1,lut%n
       lut%sigf(i)=sigma_cb(lut%rr(i)*f,z,cosm)
    END DO
    IF(verbose) WRITE(*,*) 'HALOMOD_INIT: sigf tables filled'  

    !Fill virial radius table using real radius table
    Dv=Delta_v(z,cosm)
    lut%rv=lut%rr/(Dv**(1./3.))

    IF(verbose) THEN
       WRITE(*,*) 'HALOMOD_INIT: virial radius tables filled'
       WRITE(*,*) 'HALOMOD_INIT: Delta_v:', REAL(Dv)
       WRITE(*,*) 'HALOMOD_INIT: minimum nu:', REAL(lut%nu(1))
       WRITE(*,*) 'HALOMOD_INIT: maximum nu:', REAL(lut%nu(lut%n))
       WRITE(*,*) 'HALOMOD_INIT: minimum R_v [Mpc/h]:', REAL(lut%rv(1))
       WRITE(*,*) 'HALOMOD_INIT: maximum R_v [Mpc/h]:', REAL(lut%rv(lut%n))
       WRITE(*,*) 'HALOMOD_INIT: minimum M [Msun/h]:', REAL(lut%m(1))
       WRITE(*,*) 'HALOMOD_INIT: maximum M [Msun/h]:', REAL(lut%m(lut%n))
    END IF

    lut%gmin=1.-integrate(lut%nu(1),10.,gnu,acc,3)
    lut%gmax=integrate(lut%nu(lut%n),10.,gnu,acc,3)
    lut%gbmin=1.-integrate(lut%nu(1),10.,gnubnu,acc,3)
    lut%gbmax=integrate(lut%nu(lut%n),10.,gnubnu,acc,3)
    IF(verbose) THEN
       WRITE(*,*) 'HALOMOD_INIT: missing g(nu) at low end:', REAL(lut%gmin)
       WRITE(*,*) 'HALOMOD_INIT: missing g(nu) at high end:', REAL(lut%gmax)
       WRITE(*,*) 'HALOMOD_INIT: missing g(nu)b(nu) at low end:', REAL(lut%gbmin)
       WRITE(*,*) 'HALOMOD_INIT: missing g(nu)b(nu) at high end:', REAL(lut%gbmax)
    END IF

    !Find non-linear radius and scale
    lut%rnl=r_nl(lut)
    lut%knl=1./lut%rnl

    IF(verbose) THEN
       WRITE(*,*) 'HALOMOD_INIT: non-linear radius [Mpc/h]:', REAL(lut%rnl)
       WRITE(*,*) 'HALOMOD_INIT: non-linear wavenumber [h/Mpc]:', REAL(lut%knl)
    END IF

    lut%neff=effective_index(lut,cosm)

    IF(verbose) WRITE(*,*) 'HALOMOD_INIT: n_eff:', REAL(lut%neff)

    CALL conc_bull(z,cosm,lut)

    IF(verbose) THEN
       WRITE(*,*) 'HALOMOD_INIT: halo concentration tables filled'
       WRITE(*,*) 'HALOMOD_INIT: minimum concentration:', REAL(lut%c(lut%n))
       WRITE(*,*) 'HALOMOD_INIT: maximum concentration:', REAL(lut%c(1))
    END IF

    A0=one_halo_amplitude(z,lut,cosm)
    IF(verbose) THEN
       WRITE(*,*) 'HALOMOD_INIT: A0 [Mpc/h]^3:', REAL(A0)
       WRITE(*,*) 'HALOMOD_INIT: M* [Msun/h]:', REAL(A0*comoving_matter_density(cosm))
       WRITE(*,*) 'HALOMOD_INIT: Done'
       WRITE(*,*)
    END IF

    rhom=comoving_matter_density(cosm)
    rhoc=comoving_critical_density(z,cosm)

    !Calculate Delta = 200, 500 and Delta_c = 200, 500 quantities
    CALL convert_mass_definition(lut%rv,lut%c,lut%m,Dv,1.,lut%r500,lut%c500,lut%m500,500.,1.,lut%n)
    CALL convert_mass_definition(lut%rv,lut%c,lut%m,Dv,1.,lut%r200,lut%c200,lut%m200,200.,1.,lut%n)
    CALL convert_mass_definition(lut%rv,lut%c,lut%m,Dv,rhom,lut%r500c,lut%c500c,lut%m500c,500.,rhoc,lut%n)
    CALL convert_mass_definition(lut%rv,lut%c,lut%m,Dv,rhom,lut%r200c,lut%c200c,lut%m200c,200.,rhoc,lut%n)

    CALL print_halomodel_parameters(z,lut,cosm)

    IF(verbose) verbose=.FALSE.

  END SUBROUTINE halomod_init

  FUNCTION one_halo_amplitude(z,lut,cosm)

    IMPLICIT NONE
    REAL :: one_halo_amplitude
    REAL, INTENT(IN) :: z
    TYPE(tables), INTENT(IN) :: lut
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: wk(2,lut%n)
    INTEGER :: i

    REAL, PARAMETER :: ksmall=1e-3 !A 'small' wavenumber

    DO i=1,lut%n
       wk(1,i)=lut%m(i)
    END DO
    wk(1,:)=wk(1,:)/comoving_matter_density(cosm)
    wk(2,:)=wk(1,:)

    one_halo_amplitude=p_1h(wk,ksmall,z,lut,cosm)
    one_halo_amplitude=one_halo_amplitude/(4.*pi*(ksmall/(2.*pi))**3)

  END FUNCTIOn one_halo_amplitude

  SUBROUTINE convert_mass_definition(ri,ci,mi,Di,rhoi,rj,cj,mj,Dj,rhoj,n)

    !Converts mass definition from Delta_i rho_i overdense to Delta_j rho_j overdense
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(IN) :: ri(n), ci(n), mi(n)
    REAL, INTENT(OUT) :: rj(n), cj(n), mj(n)
    REAL, INTENT(IN) :: Di, Dj, rhoi, rhoj
    REAL, ALLOCATABLE :: LHS(:), RHS(:)
    REAL :: rs
    INTEGER :: i

    !Ensure these are all zero
    rj=0.
    cj=0.
    mj=0.

    !Allocate arrays for the LHS and RHS of the equation
    ALLOCATE(LHS(n),RHS(n))

    !Fill arrays for LHS and RHS of the equation - can use same r(i) table
    !The equation: (r_i^3 x rho_i x Delta_i / X_i(r_i/rs) = same for j)
    DO i=1,n
       rs=ri(i)/ci(i)
       LHS(i)=(ri(i)**3)*Di*rhoi/nfw_factor(ri(i)/rs)
       RHS(i)=(ri(i)**3)*Dj*rhoj/nfw_factor(ri(i)/rs)
       !RHS(i)=LHS(i)*Dj*rhoj/(Di*rhoi)
       !WRITE(*,fmt='(I5,3ES15.5)') i, ri(i), LHS(i), RHS(i)
    END DO

    !Now use the find algorithm to invert L(r_i)=R(r_j) so that
    !r_j=R^{-1}[L(r_i)]
    DO i=1,n

       !First find the radius
       rj(i)=exp(find_solve(log(LHS(i)),log(ri),log(RHS),n))

       !This is to check the solution is correct
       !LH=LHS(i)
       !RH=exp(find(log(rj(i)),log(ri),log(RHS),n,3,3,2))
       !WRITE(*,fmt='(I5,2F15.5)') i, LH, RH

       !NOTE VERY WELL - this does *NOT* mean that:
       !LHS(i)=(rj(i)**3)*Dj*rhoj/nfw_factor(rj(i)/rs)
       !Because the integer 'i' does not correspond to the solution
       !LH=LHS(i)
       !RH=(rj(i)**3)*Dj*rhoj/nfw_factor(rj(i)/rs)
       !WRITE(*,fmt='(I5,2F15.5)') i, LH, RH

       !Now do concentration and mass
       rs=ri(i)/ci(i)
       cj(i)=rj(i)/rs
       mj(i)=mi(i)*nfw_factor(cj(i))/nfw_factor(ci(i))

    END DO

    DEALLOCATE(LHS,RHS)

  END SUBROUTINE convert_mass_definition

  FUNCTION find_solve(a,xtab,ytab,n)

    !Solves y(x)=a for x
    IMPLICIT NONE
    REAL :: find_solve
    REAL, INTENT(IN) :: a, xtab(n), ytab(n)
    INTEGER, INTENT(IN) :: n

    find_solve=find(a,ytab,xtab,n,3,3,2)

  END FUNCTION find_solve

  FUNCTION nfw_factor(x)

    !The NFW 'mass' factor that crops up all the time
    IMPLICIT NONE
    REAL :: nfw_factor
    REAL, INTENT(IN) :: x

    nfw_factor=log(1.+x)-x/(1.+x)

  END FUNCTION nfw_factor

  FUNCTION radius_m(m,cosm)

    !The comoving radius corresponding to mass M in a homogeneous universe
    IMPLICIT NONE
    REAL :: radius_m
    REAL, INTENT(IN) :: m
    TYPE(cosmology), INTENT(IN) :: cosm

    radius_m=(3.*m/(4.*pi*comoving_matter_density(cosm)))**(1./3.)

  END FUNCTION radius_m

  FUNCTION effective_index(lut,cosm)

    !Power spectrum slope a the non-linear scale
    IMPLICIT NONE
    REAL :: effective_index
    TYPE(cosmology) :: cosm
    TYPE(tables) :: lut

    !Numerical differentiation to find effective index at collapse
    effective_index=-3.-derivative_table(log(lut%rnl),log(lut%rr),log(lut%sig**2),lut%n,3,3)

    !For some bizarre cosmologies r_nl is very small, so almost no collapse has occured
    !In this case the n_eff calculation goes mad and needs to be fixed using this fudge.
    IF(effective_index<cosm%n-4.) effective_index=cosm%n-4.
    IF(effective_index>cosm%n)      effective_index=cosm%n

  END FUNCTION effective_index

  SUBROUTINE conc_bull(z,cosm,lut)

    !Calculates the Bullock et al. (2001) mass-concentration relation
    IMPLICIT NONE
    REAL, INTENT(IN) :: z
    TYPE(cosmology) :: cosm, cos_lcdm
    TYPE(tables) :: lut
    REAL :: A, zinf, ainf, zf, g_lcdm, g_wcdm
    INTEGER :: i   

    A=As(cosm)

    !Fill the collapse z look-up table
    CALL zcoll_bull(z,cosm,lut)

    !Fill the concentration look-up table
    DO i=1,lut%n

       zf=lut%zc(i)
       lut%c(i)=A*(1.+zf)/(1.+z)

       !Dolag2004 prescription for adding DE dependence
       IF(imead==1) THEN

          !IF((cosm%w .NE. -1.) .OR. (cosm%wa .NE. 0)) THEN

          !The redshift considered to be infinite
          zinf=10.
          !ainf=1./(1.+zinf)
          ainf=scale_factor_z(zinf)

          !Save the growth function in the current cosmology
          g_wcdm=grow(zinf,cosm)

          !Make a LCDM cosmology
          cos_lcdm=cosm
          DEALLOCATE(cos_lcdm%growth)
          DEALLOCATE(cos_lcdm%a_growth)
          cos_lcdm%w=-1.
          cos_lcdm%wa=0.
          cos_lcdm%om_v=1.-cosm%om_m !Added this so that 'making a LCDM cosmology' works for curved models.

          !Needs to use grow_int explicitly in case tabulated values are stored
          g_lcdm=growint(ainf,cos_lcdm)

          !Changed this to a power of 1.5, which produces more accurate results for extreme DE
          lut%c(i)=lut%c(i)*((g_wcdm/g_lcdm)**1.5)

       END IF

    END DO

  END SUBROUTINE conc_bull

  SUBROUTINE zcoll_bull(z,cosm,lut)

    !This fills up the halo collapse redshift table as per Bullock relations   
    IMPLICIT NONE
    REAL, INTENT(IN) :: z
    TYPE(cosmology) :: cosm
    TYPE(tables) :: lut
    REAL :: dc
    REAL :: af, zf, RHS, a, growz
    REAL, ALLOCATABLE :: af_tab(:), grow_tab(:)
    INTEGER :: i, ntab  

    ntab=SIZE(cosm%growth)
    ALLOCATE(af_tab(ntab),grow_tab(ntab))

    af_tab=cosm%a_growth
    grow_tab=cosm%growth

    !Do numerical inversion
    DO i=1,lut%n

       !I don't think this is really consistent with dc varying as a function of z
       !But the change will be very small
       dc=delta_c(z,cosm)

       RHS=dc*grow(z,cosm)/lut%sigf(i)

       !a=1./(1.+z)
       a=scale_factor_z(z)
       growz=find(a,af_tab,grow_tab,cosm%ng,3,3,2)

       IF(RHS>growz) THEN
          zf=z
       ELSE
          af=find(RHS,grow_tab,af_tab,cosm%ng,3,3,2)
          !zf=-1.+1./af
          zf=redshift_a(af)
       END IF

       lut%zc(i)=zf

    END DO

    DEALLOCATE(af_tab,grow_tab)

  END SUBROUTINE zcoll_bull

  FUNCTION mass_r(r,cosm)

    !Calcuates the mass contains in a sphere of comoving radius 'r' in a homogeneous universe
    IMPLICIT NONE
    REAL :: mass_r, r
    TYPE(cosmology) :: cosm

    !Relation between mean cosmological mass and radius
    mass_r=(4.*pi/3.)*comoving_matter_density(cosm)*(r**3)

  END FUNCTION mass_r

  FUNCTION comoving_critical_density(z,cosm)

    !Comoving critical density in (Msun/h) / (Mpc/h)^3
    !This constant is (3/8pi) x H0^2/G
    IMPLICIT NONE
    REAL :: comoving_critical_density
    REAL, INTENT(IN) :: z
    TYPE(cosmology), INTENT(IN) :: cosm

    comoving_critical_density=physical_critical_density(z,cosm)/(1.+z)**3

  END FUNCTION comoving_critical_density

  FUNCTION physical_critical_density(z,cosm)

    !Physical critical density in (Msun/h) / (Mpc/h)^3
    !This constant is (3/8pi) x H0^2/G
    IMPLICIT NONE
    REAL :: physical_critical_density
    REAL, INTENT(IN) :: z
    TYPE(cosmology), INTENT(IN) :: cosm

    physical_critical_density=critical_density*Hubble2(z,cosm)

  END FUNCTION physical_critical_density

  FUNCTION comoving_matter_density(cosm)

    !Comoving matter density in (Msun/h) / (Mpc/h)^3
    !This constant is (3/8pi) x H0^2/G x Omega_m(z=0)
    !Not a function of redshift!
    IMPLICIT NONE
    REAL :: comoving_matter_density
    TYPE(cosmology), INTENT(IN) :: cosm

    comoving_matter_density=critical_density*cosm%om_m

  END FUNCTION comoving_matter_density

  FUNCTION physical_matter_density(z,cosm)

    !Physical matter density in (Msun/h) / (Mpc/h)^3
    !This constant is (3/8pi) x H0^2/G x Omega_m(z=0)
    IMPLICIT NONE
    REAL :: physical_matter_density
    REAL, INTENT(IN) :: z
    TYPE(cosmology), INTENT(IN) :: cosm

    physical_matter_density=comoving_matter_density(cosm)*(1.+z)**3

  END FUNCTION physical_matter_density

  FUNCTION Tk(k,cosm)

    !Transfer function selection
    IMPLICIT NONE
    REAL :: Tk, k
    TYPE(cosmology) :: cosm

    Tk=Tk_eh(k,cosm)

  END FUNCTION Tk

  FUNCTION Tk_eh(yy,cosm)

    !Eisenstein & Hu fitting function
    !JP - the astonishing D.J. Eisenstein & W. Hu fitting formula (ApJ 496 605 [1998])
    !JP - remember I use k/h, whereas they use pure k, om_m is cdm + baryons
    IMPLICIT NONE
    REAL :: Tk_eh
    REAL, INTENT(IN) :: yy
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: rk, e, thet, b1, b2, zd, ze, rd, re, rke, s, rks
    REAL :: q
    REAL :: y, g, ab
    REAL :: a1, a2, ac
    REAL :: bc
    REAL :: f, fac
    REAL :: c1, c2, tc
    REAL :: bb, bn, ss, tb
    REAL :: om_m, om_b, h

    om_m=cosm%om_m
    om_b=cosm%om_b
    h=cosm%h

    rk=yy*h

    e=exp(1.)

    thet=2.728/2.7
    b1=0.313*(om_m*h*h)**(-0.419)*(1+0.607*(om_m*h*h)**0.674)
    b2=0.238*(om_m*h*h)**0.223
    zd=1291.*(1+b1*(om_b*h*h)**b2)*(om_m*h*h)**0.251/(1.+0.659*(om_m*h*h)**0.828)
    ze=2.50e4*om_m*h*h/thet**4.
    rd=31500.*om_b*h*h/thet**4./zd !Should this be 1+zd (Steven Murray enquirey)?
    re=31500.*om_b*h*h/thet**4./ze
    rke=7.46e-2*om_m*h*h/thet**2.
    s=(2./3./rke)*sqrt(6./re)*log((sqrt(1.+rd)+sqrt(rd+re))/(1+sqrt(re)))
    rks=1.6*( (om_b*h*h)**0.52 ) * ( (om_m*h*h)**0.73 ) * (1.+(10.4*om_m*h*h)**(-0.95))

    q=rk/13.41/rke

    y=(1.+ze)/(1.+zd)
    g=y*(-6.*sqrt(1+y)+(2.+3.*y)*log((sqrt(1.+y)+1.)/(sqrt(1.+y)-1.)))
    ab=g*2.07*rke*s/(1.+rd)**(0.75)

    a1=(46.9*om_m*h*h)**0.670*(1+(32.1*om_m*h*h)**(-0.532))
    a2=(12.0*om_m*h*h)**0.424*(1+(45.0*om_m*h*h)**(-0.582))
    ac=(a1**(-om_b/om_m)) * (a2**(-(om_b/om_m)**3.))

    b1=0.944/(1+(458.*om_m*h*h)**(-0.708))
    b2=(0.395*om_m*h*h)**(-0.0266)
    bc=1./(1.+b1*((1.-om_b/om_m)**b2-1.))

    f=1./(1.+(rk*s/5.4)**4.)

    c1=14.2 + 386./(1.+69.9*q**1.08)
    c2=14.2/ac + 386./(1.+69.9*q**1.08)
    tc=f*log(e+1.8*bc*q)/(log(e+1.8*bc*q)+c1*q*q) +(1.-f)*log(e+1.8*bc*q)/(log(e+1.8*bc*q)+c2*q*q)

    bb=0.5+(om_b/om_m) + (3.-2.*om_b/om_m)*sqrt((17.2*om_m*h*h)**2.+1.)
    bn=8.41*(om_m*h*h)**0.435
    ss=s/(1.+(bn/rk/s)**3.)**(1./3.)
    tb=log(e+1.8*q)/(log(e+1.8*q)+c1*q*q)/(1+(rk*s/5.2)**2.)
    !IF((rk/rks**1.4)>7.) THEN
    !   fac=0.
    !ELSE
    !Removed this IF statement as it produced a discontinuity in P_lin(k) as cosmology
    !was varied - thanks David Copeland for pointing this out
    fac=exp(-(rk/rks)**1.4)
    !END IF
    tb=(tb+ab*fac/(1.+(bb/rk/s)**3.))*sin(rk*ss)/rk/ss

    tk_eh=real((om_b/om_m)*tb+(1-om_b/om_m)*tc)

  END FUNCTION TK_EH

  FUNCTION p_lin(k,z,cosm)

    !Linear matter power spectrum
    !P(k) should have been previously normalised so as to get the amplitude 'A' correct
    IMPLICIT NONE
    REAL :: p_lin
    REAL, INTENT (IN) :: k, z
    TYPE(cosmology), INTENT(IN) :: cosm 

    IF(k==0.) THEN
       !If p_lin happens to be foolishly called for 0 mode (which should never happen, but might in integrals)
       p_lin=0.
    ELSE IF(k>1.e8) THEN
       !Avoids some issues if p_lin is called for very (absurdly) high k values
       !For some reason crashes can occur if this is the case
       p_lin=0.
    ELSE IF(ibox==1 .AND. k<2.*pi/Lbox) THEN
       p_lin=0.
    ELSE
       !In this case look for the transfer function
       p_lin=(cosm%A**2)*(grow(z,cosm)**2)*(Tk(k,cosm)**2)*(k**(cosm%n+3.))
    END IF

  END FUNCTION p_lin

  SUBROUTINE halomod(ih1,ih2,k,z,p2h,p1h,pfull,plin,lut,cosm)

    !Gets the one- and two-halo terms and combines them
    IMPLICIT NONE
    REAL, INTENT(OUT) :: p1h, p2h, pfull
    REAL, INTENT(IN) :: plin, k, z
    INTEGER, INTENT(IN) :: ih1, ih2
    TYPE(cosmology), INTENT(IN) :: cosm
    TYPE(tables), INTENT(IN) :: lut
    REAL :: alp
    REAL :: wk(2,lut%n), m, rv, rs
    INTEGER :: i, j, ih(2)

    !Initially fill this small array 
    ih(1)=ih1
    ih(2)=ih2

    !For the i's
    !-1 - DMonly
    ! 0 - All matter
    ! 1 - CDM
    ! 2 - Gas
    ! 3 - Stars
    ! 4 - Bound gas
    ! 5 - Free gas
    ! 6 - Pressure

    !Calls expressions for one- and two-halo terms and then combines
    !to form the full power spectrum
    IF(k==0.) THEN

       !This should really never be called for k=0
       p1h=0.
       p2h=0.

    ELSE

       !Calculate the halo window functions
       DO j=1,2
          DO i=1,lut%n
             m=lut%m(i)
             rv=lut%rv(i)
             rs=rv/lut%c(i)
             wk(j,i)=win_type(1,ih(j),k,m,rv,rs,z,lut,cosm)
          END DO
          IF(ih(2)==ih(1)) THEN
             !Avoid having to call win_type twice if doing auto spectrum
             wk(2,:)=wk(1,:)
             EXIT
          END IF
       END DO

       !Get the one-halo term
       p1h=p_1h(wk,k,z,lut,cosm)

       !Only if imead=-1 do we need to recalcualte the window
       !functions for the two-halo term with k=0 fixed
       IF(imead==-1) THEN
          DO j=1,2
             DO i=1,lut%n
                m=lut%m(i)
                rv=lut%rv(i)
                rs=rv/lut%c(i)
                wk(j,i)=win_type(1,ih(j),0.,m,rv,rs,z,lut,cosm)
             END DO
             IF(ih(2)==ih(1)) THEN
                !Avoid having to call win_type twice if doing auto spectrum
                wk(2,:)=wk(1,:)
                EXIT
             END IF
          END DO
       END IF

       !Get the two-halo term
       p2h=p_2h(ih,wk,k,z,plin,lut,cosm)

    END IF

    !Construct the 'full' halo-model power spectrum
    IF(imead==0 .OR. imead==-1) THEN
       pfull=p2h+p1h
    ELSE IF(imead==1) THEN
       alp=alpha_transition(lut,cosm)
       pfull=(p2h**alp+p1h**alp)**(1./alp)
    END IF

    !If we are worrying about voids
    IF(void) THEN
       pfull=pfull+p_1v(k,lut)
    END IF

  END SUBROUTINE halomod

  FUNCTION p_2h(ih,wk,k,z,plin,lut,cosm)

    !Produces the 'two-halo' power
    IMPLICIT NONE
    REAL :: p_2h
    REAL, INTENT(IN) :: k, plin
    REAL, INTENT(IN) :: z
    TYPE(tables), INTENT(IN) :: lut
    REAL, INTENT(IN) :: wk(2,lut%n)
    TYPE(cosmology), INTENT(IN) :: cosm
    INTEGER, INTENT(IN) :: ih(2)
    REAL :: sigv, frac, rhom
    REAL, ALLOCATABLE :: integrand11(:), integrand12(:)
    REAL, ALLOCATABLE :: integrand21(:), integrand22(:)
    REAL :: nu, m, wk1, wk2, m0
    REAL :: sum11, sum12
    REAL :: sum21, sum22
    INTEGER :: i

    rhom=comoving_matter_density(cosm)

    IF(imead==0 .OR. imead==-1) THEN

       ALLOCATE(integrand11(lut%n),integrand12(lut%n))

       IF(ibias==2) THEN
          !Only necessary for second-order bias integral
          ALLOCATE(integrand21(lut%n),integrand22(lut%n))
       END IF

       DO i=1,lut%n

          m=lut%m(i)
          nu=lut%nu(i)

          !Linear bias term
          integrand11(i)=gnu(nu)*bnu(nu)*wk(1,i)/m
          integrand12(i)=gnu(nu)*bnu(nu)*wk(2,i)/m

          IF(ibias==2) THEN
             !Second-order bias term
             integrand21(i)=gnu(nu)*b2nu(nu)*wk(1,i)/m
             integrand22(i)=gnu(nu)*b2nu(nu)*wk(2,i)/m
          END IF

       END DO

       !Evaluate these integrals from the tabled values
       sum11=integrate_table(lut%nu,integrand11,lut%n,1,lut%n,3)
       sum12=integrate_table(lut%nu,integrand12,lut%n,1,lut%n,3)

       IF(ip2h==0) THEN
          !Do nothing in this case
       ELSE IF(ip2h==1) THEN
          !Add on the value of integral b(nu)*g(nu) assuming w=1
          sum11=sum11+lut%gbmin*halo_fraction(ih(1),m,cosm)/rhom
          sum12=sum12+lut%gbmin*halo_fraction(ih(2),m,cosm)/rhom
       ELSE IF(ip2h==2) THEN
          !Put the missing part of the integrand as a delta function at nu1
          m0=lut%m(1)
          wk1=wk(1,1)
          wk2=wk(2,1)
          sum11=sum11+lut%gbmin*wk1/m0
          sum12=sum12+lut%gbmin*wk2/m0
       ELSE
          STOP 'P_2h: Error, ip2h not specified correctly'
       END IF

       p_2h=plin*sum11*sum12*(rhom**2)

       IF(ibias==2) THEN
          !Second order bias correction
          !This needs to have the property that \int f(nu)b2(nu) du = 0
          !This means it is hard to check that the normalisation is correct
          !e.g., how much do low mass haloes matter
          !Varying mmin does make a difference to the values of the integrals
          sum21=integrate_table(lut%nu,integrand21,lut%n,1,lut%n,3)
          sum22=integrate_table(lut%nu,integrand22,lut%n,1,lut%n,3)
          p_2h=p_2h+(plin**2)*sum21*sum22*(rhom**2)
       END IF

    ELSE IF(imead==1) THEN

       sigv=lut%sigv
       frac=fdamp(z,cosm)

       IF(frac==0.) THEN
          p_2h=plin
       ELSE
          p_2h=plin*(1.-frac*(tanh(k*sigv/sqrt(ABS(frac))))**2)
       END IF

       !For some extreme cosmologies frac>1. so this must be added to prevent p_2h<0.
       IF(p_2h<0.) p_2h=0.

    END IF

  END FUNCTION p_2h

  FUNCTION p_1h(wk,k,z,lut,cosm)

    !Calculates the one-halo term
    IMPLICIT NONE
    REAL :: p_1h
    REAL, INTENT(IN) :: k, z
    TYPE(tables), INTENT(IN) :: lut
    REAL, INTENT(IN) :: wk(2,lut%n)
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: m, g, fac, et, ks
    REAL, ALLOCATABLE :: integrand(:)
    INTEGER :: i

    ALLOCATE(integrand(lut%n))
    integrand=0.

    !Only call eta once
    et=eta(z,cosm)

    !Calculates the value of the integrand at all nu values!
    DO i=1,lut%n
       g=gnu(lut%nu(i))
       m=lut%m(i)
       integrand(i)=g*wk(1,i)*wk(2,i)/m
    END DO

    !Carries out the integration
    !Important to use basic trapezium rule because the integrand is messy due to rapid oscillations in W(k)
    p_1h=comoving_matter_density(cosm)*integrate_table(lut%nu,integrand,lut%n,1,lut%n,1)*(4.*pi)*(k/(2.*pi))**3

    DEALLOCATE(integrand)

    IF(imead==1) THEN

       !Damping of the 1-halo term at very large scales
       ks=kstar(lut,cosm)

       !Prevents problems if k/ks is very large

       IF(ks>0.) THEN

          IF((k/ks)**2>7.) THEN
             fac=0.
          ELSE
             fac=exp(-((k/ks)**2))
          END IF

          p_1h=p_1h*(1.-fac)

       END IF

    END IF

  END FUNCTION p_1h

  FUNCTION p_1v(k,lut)!,cosm)

    IMPLICIT NONE
    REAL :: p_1v
    REAL, INTENT(IN) :: k
    TYPE(tables), INTENT(IN) :: lut
    !TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: dc, wk, V, rvoid, rcomp, nu
    REAL :: integrand(lut%n)
    INTEGER :: i, n

    !Parameters
    REAL, PARAMETER :: dv=-1.
    REAL, PARAMETER :: fvoid=1.1
    LOGICAL, PARAMETER :: compensate=.TRUE.
    LOGICAL, PARAMETER :: simple=.FALSE.

    IF(simple) THEN
       n=1
    ELSE
       n=lut%n
    END IF

    DO i=1,n

       !Get the void radius and compensation radius
       IF(simple) THEn
          rvoid=10.
       ELSE         
          rvoid=lut%rr(i)
          nu=lut%nu(i)        
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
          integrand(i)=gnu(nu)*wk**2/V
       END IF

    END DO

    !Calculate the void one-halo term
    IF(simple) THEN
       p_1v=wk**2/V
    ELSE
       p_1v=integrate_table(lut%nu,integrand,n,1,n,1)
    END IF

    p_1v=p_1v*(4.*pi)*(k/(2.*pi))**3

  END FUNCTION p_1v

  SUBROUTINE fill_sigtab(cosm)

    !This fills up tables of r vs. sigma(r) across a range in r!
    !It is used only in look-up for further calculations of sigma(r) and not otherwise!
    !and prevents a large number of calls to the sigint functions
    IMPLICIT NONE
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL, ALLOCATABLE :: rtab(:), sigtab(:)
    REAL :: r, sig
    INTEGER :: i

    !These values of 'r' work fine for any power spectrum of cosmological importance
    !Having nsig as a 2** number is most efficient for the look-up routines
    !rmin and rmax need to be decided in advance and are chosen such that
    !R vs. sigma(R) is a power-law below and above these values of R   
    INTEGER, PARAMETER :: nsig=64 !Number of entries for sigma(R) tables
    REAL, PARAMETER :: rmin=1e-4 !Minimum r value (NB. sigma(R) needs to be power-law below)
    REAL, PARAMETER :: rmax=1e3 !Maximum r value (NB. sigma(R) needs to be power-law above)

    !These must be not allocated before sigma calculations otherwise when sigma(r) is called
    !otherwise sigma(R) looks for the result in the tables
    !IF(ALLOCATED(cosm%r_sigma)) DEALLOCATE(cosm%r_sigma)
    !IF(ALLOCATED(cosm%sigma)) DEALLOCATE(cosm%sigma)
    IF(ALLOCATED(cosm%logr_sigma)) DEALLOCATE(cosm%logr_sigma)
    IF(ALLOCATED(cosm%logsigma)) DEALLOCATE(cosm%logsigma)

    cosm%nsig=nsig
    ALLOCATE(rtab(nsig),sigtab(nsig))

    IF(verbose) THEN
       WRITE(*,*) 'SIGTAB: Filling sigma interpolation table'
       WRITE(*,*) 'SIGTAB: R minimum [Mpc/h]:', rmin !If I put REAL() here I get an error with goslow for some reason!?
       WRITE(*,*) 'SIGTAB: R maximum [Mpc/h]:', REAL(rmax)
       WRITE(*,*) 'SIGTAB: number of points:', nsig
    END IF

    DO i=1,nsig

       !Equally spaced r in log
       !r=exp(log(rmin)+log(rmax/rmin)*float(i-1)/float(nsig-1))
       r=exp(progression(log(rmin),log(rmax),i,nsig))

       sig=sigma(r,0.,cosm)

       rtab(i)=r
       sigtab(i)=sig

    END DO

    !Must be allocated after the sigtab calulation above
    ALLOCATE(cosm%logr_sigma(nsig),cosm%logsigma(nsig))

    cosm%logr_sigma=log(rtab)
    cosm%logsigma=log(sigtab)

    DEALLOCATE(rtab,sigtab)

    IF(verbose) THEN
       WRITE(*,*) 'SIGTAB: Done'
       WRITE(*,*)
    END IF

  END SUBROUTINE fill_sigtab

  FUNCTION sigma(r,z,cosm)

    !Gets sigma(R)
    IMPLICIT NONE
    REAL :: sigma
    REAL, INTENT(IN) :: r, z
    TYPE(cosmology), INTENT(IN) :: cosm

    !REAL, PARAMETER :: acc=1d-4
    INTEGER, PARAMETER :: iorder=3
    REAL, PARAMETER :: rsplit=1e-2

    IF(r>=rsplit) THEN
       sigma=sqrt(sigint0(r,z,cosm,acc,iorder))
    ELSE IF(r<rsplit) THEN
       sigma=sqrt(sigint1(r,z,cosm,acc,iorder)+sigint2(r,z,cosm,acc,iorder))
    ELSE
       STOP 'SIGMA: Error, something went wrong'
    END IF

  END FUNCTION sigma

  FUNCTION sigma_integrand(k,R,z,cosm)

    !The integrand for the sigma(R) integrals
    IMPLICIT NONE
    REAL :: sigma_integrand
    REAL, INTENT(IN) :: k, R, z
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: y, w_hat

    IF(k==0.) THEN
       sigma_integrand=0.
    ELSE
       y=k*R
       w_hat=wk_tophat(y)
       sigma_integrand=p_lin(k,z,cosm)*(w_hat**2)/k
    END IF

  END FUNCTION sigma_integrand

  FUNCTION sigma_integrand_transformed(t,R,f,z,cosm)

    !The integrand for the sigma(R) integrals
    IMPLICIT NONE
    REAL :: sigma_integrand_transformed
    REAL, INTENT(IN) :: t, R, z
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: k, y, w_hat

    INTERFACE
       FUNCTION f(x)
         REAL :: f
         REAL, INTENT(IN) :: x
       END FUNCTION f
    END INTERFACE

    !Integrand to the sigma integral in terms of t. Defined by k=(1/t-1)/f(R) where f(R) is *any* function

    IF(t==0.) THEN
       !t=0 corresponds to k=infintiy when W(kR)=0.
       sigma_integrand_transformed=0.
    ELSE IF(t==1.) THEN
       !t=1 corresponds to k=0. when P(k)=0.
       sigma_integrand_transformed=0.
    ELSE
       !f(R) can be *any* function of R here to improve integration speed
       k=(-1.+1./t)/f(R)
       y=k*R
       w_hat=wk_tophat(y)
       sigma_integrand_transformed=p_lin(k,z,cosm)*(w_hat**2)/(t*(1.-t))
    END IF

  END FUNCTION sigma_integrand_transformed

  FUNCTION sigint0(r,z,cosm,acc,iorder)

    !Integrates between a and b until desired accuracy is reached
    !Stores information to reduce function calls
    IMPLICIT NONE
    REAL :: sigint0
    REAL, INTENT(IN) :: r, z
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL, INTENT(IN) :: acc
    INTEGER, INTENT(IN) :: iorder
    INTEGER :: i, j
    INTEGER :: n
    REAL :: x, dx
    REAL :: f1, f2, fx
    DOUBLE PRECISION :: sum_n, sum_2n, sum_new, sum_old

    INTEGER, PARAMETER :: jmin=5
    INTEGER, PARAMETER :: jmax=30
    REAL, PARAMETER :: a=0. !Integration lower limit (corresponts to k=inf)
    REAL, PARAMETER :: b=1. !Integration upper limit (corresponds to k=0)

    IF(a==b) THEN

       !Fix the answer to zero if the integration limits are identical
       sigint0=0.

    ELSE

       !Reset the sum variable for the integration
       sum_2n=0.
       sum_n=0.
       sum_old=0.
       sum_new=0.

       DO j=1,jmax

          !Note, you need this to be 1+2**n for some integer n
          !j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
          n=1+2**(j-1)

          !Calculate the dx interval for this value of 'n'
          dx=(b-a)/REAL(n-1)

          IF(j==1) THEN

             !The first go is just the trapezium of the end points
             f1=sigma_integrand_transformed(a,r,f0_rapid,z,cosm)
             f2=sigma_integrand_transformed(b,r,f0_rapid,z,cosm)
             sum_2n=0.5*(f1+f2)*dx
             sum_new=sum_2n

          ELSE

             !Loop over only new even points to add these to the integral
             DO i=2,n,2
                !x=a+(b-a)*REAL(i-1)/REAL(n-1)
                x=progression(a,b,i,n)
                fx=sigma_integrand_transformed(x,r,f0_rapid,z,cosm)
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
                STOP 'SIGINT0: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (ABS(-1.+sum_new/sum_old)<acc)) THEN
             !jmin avoids spurious early convergence
             sigint0=REAL(sum_new)
             EXIT
          ELSE IF(j==jmax) THEN
             STOP 'SIGINT0: Integration timed out'
          ELSE
             !Integral has not converged so store old sums and reset sum variables
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.
          END IF

       END DO

    END IF

  END FUNCTION sigint0

  FUNCTION f0_rapid(r)

    !This is the 'rapidising' function to increase integration speed
    !for sigma(R). Found by trial-and-error
    IMPLICIT NONE
    REAL :: f0_rapid
    REAL, INTENT(IN) :: r
    REAL :: alpha

    REAL, PARAMETER :: rsplit=1e-2

    IF(r>rsplit) THEN
       !alpha 0.3-0.5 works well
       alpha=0.5
    ELSE
       !If alpha=1 this goes tits up
       !alpha 0.7-0.9 works well
       alpha=0.8
    END IF

    f0_rapid=r**alpha

  END FUNCTION f0_rapid

  FUNCTION sigint1(r,z,cosm,acc,iorder)

    !Integrates between a and b until desired accuracy is reached
    !Stores information to reduce function calls
    IMPLICIT NONE
    REAL :: sigint1
    REAL, INTENT(IN) :: r, z
    TYPE(cosmology), INTENT(IN) :: cosm
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

    a=r/(r+r**.5)
    b=1.

    IF(a==b) THEN

       !Fix the answer to zero if the integration limits are identical
       sigint1=0.

    ELSE

       !Reset the sum variable for the integration
       sum_2n=0.
       sum_n=0.
       sum_old=0.
       sum_new=0.

       DO j=1,jmax

          !Note, you need this to be 1+2**n for some integer n
          !j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
          n=1+2**(j-1)

          !Calculate the dx interval for this value of 'n'
          dx=(b-a)/REAL(n-1)

          IF(j==1) THEN

             !The first go is just the trapezium of the end points
             f1=sigma_integrand_transformed(a,r,f1_rapid,z,cosm)
             f2=sigma_integrand_transformed(b,r,f1_rapid,z,cosm)
             sum_2n=0.5*(f1+f2)*dx
             sum_new=sum_2n

          ELSE

             !Loop over only new even points to add these to the integral
             DO i=2,n,2
                !x=a+(b-a)*REAL(i-1)/REAL(n-1)
                x=progression(a,b,i,n)
                fx=sigma_integrand_transformed(x,r,f1_rapid,z,cosm)
                sum_2n=sum_2n+fx
             END DO

             !Now create the total using the old and new parts
             sum_2n=sum_n/2.+sum_2n*dx

             !Now calculate the new sum depending on the integration order
             IF(iorder==1) THEN  
                sum_new=REAL(sum_2n)
             ELSE IF(iorder==3) THEN         
                sum_new=(4.*sum_2n-sum_n)/3. !This is Simpson's rule and cancels error
             ELSE
                STOP 'SIGINT1: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (ABS(-1.+sum_new/sum_old)<acc)) THEN
             !jmin avoids spurious early convergence
             sigint1=REAL(sum_new)
             EXIT
          ELSE IF(j==jmax) THEN
             STOP 'SIGINT1: Integration timed out'
          ELSE
             !Integral has not converged so store old sums and reset sum variables
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.
          END IF

       END DO

    END IF

  END FUNCTION sigint1

  FUNCTION f1_rapid(r)

    !This is the 'rapidising' function to increase integration speed
    !for sigma(R). Found by trial-and-error
    IMPLICIT NONE
    REAL :: f1_rapid
    REAL, INTENT(IN) :: r

    REAL, PARAMETER :: alpha=0.5

    f1_rapid=r**alpha

  END FUNCTION f1_rapid

  FUNCTION sigint2(r,z,cosm,acc,iorder)

    !Integrates between a and b until desired accuracy is reached
    !Stores information to reduce function calls
    IMPLICIT NONE
    REAL :: sigint2
    REAL, INTENT(IN) :: r, z
    TYPE(cosmology), INTENT(IN) :: cosm
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
    REAL, PARAMETER :: C=10. !How far to go out in 1/r units for integral

    a=1./r
    b=C/r

    IF(a==b) THEN

       !Fix the answer to zero if the integration limits are identical
       sigint2=0.

    ELSE

       !Reset the sum variable for the integration
       sum_2n=0.
       sum_n=0.
       sum_old=0.
       sum_new=0.

       DO j=1,jmax

          !Note, you need this to be 1+2**n for some integer n
          !j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
          n=1+2**(j-1)

          !Calculate the dx interval for this value of 'n'
          dx=(b-a)/REAL(n-1)

          IF(j==1) THEN

             !The first go is just the trapezium of the end points
             f1=sigma_integrand(a,r,z,cosm)
             f2=sigma_integrand(b,r,z,cosm)
             sum_2n=0.5*(f1+f2)*dx
             sum_new=sum_2n

          ELSE

             !Loop over only new even points to add these to the integral
             DO i=2,n,2
                !x=a+(b-a)*REAL(i-1)/REAL(n-1)
                x=progression(a,b,i,n)
                fx=sigma_integrand(x,r,z,cosm)
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
                STOP 'SIGINT2: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (ABS(-1.+sum_new/sum_old)<acc)) THEN
             !jmin avoids spurious early convergence
             sigint2=REAL(sum_new)
             !WRITE(*,*) 'INTEGRATE_STORE: Nint:', n
             EXIT
          ELSE IF(j==jmax) THEN
             STOP 'SIGINT2: Integration timed out'
          ELSE
             !Integral has not converged so store old sums and reset sum variables
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.
          END IF

       END DO

    END IF

  END FUNCTION sigint2

  FUNCTION sigma_cb(r,z,cosm)

    !USE cosdef
    IMPLICIT NONE
    REAL :: sigma_cb
    REAL, INTENT(IN) :: r, z
    TYPE(cosmology), INTENT(IN) :: cosm

    !Finds sigma_cold from look-up tables
    !In this version sigma_cold=sigma

    !sigma_cb=grow(z,cosm)*exp(find(log(r),log(cosm%r_sigma),log(cosm%sigma),cosm%nsig,3,3,2))
    sigma_cb=grow(z,cosm)*exp(find(log(r),cosm%logr_sigma,cosm%logsigma,cosm%nsig,3,3,2))

  END FUNCTION sigma_cb

!!$  FUNCTION inttab(x,y,n,iorder)
!!$
!!$    !Integrates tables y(x)dx
!!$    IMPLICIT NONE
!!$    REAL :: inttab
!!$    INTEGER, INTENT(IN) :: n
!!$    REAL, INTENT(IN) :: x(n), y(n)
!!$    REAL :: a, b, c, d, h
!!$    REAL :: q1, q2, q3, qi, qf
!!$    REAL :: x1, x2, x3, x4, y1, y2, y3, y4, xi, xf
!!$    REAL :: sum
!!$    INTEGER :: i, i1, i2, i3, i4
!!$    INTEGER, INTENT(IN) :: iorder
!!$
!!$    sum=0.
!!$
!!$    IF(iorder==1) THEN
!!$
!!$       !Sums over all Trapezia (a+b)*h/2
!!$       DO i=1,n-1
!!$          a=y(i+1)
!!$          b=y(i)
!!$          h=x(i+1)-x(i)
!!$          sum=sum+(a+b)*h/2.
!!$       END DO
!!$
!!$    ELSE IF(iorder==2) THEN
!!$
!!$       DO i=1,n-2
!!$
!!$          x1=x(i)
!!$          x2=x(i+1)
!!$          x3=x(i+2)
!!$
!!$          y1=y(i)
!!$          y2=y(i+1)
!!$          y3=y(i+2)
!!$
!!$          CALL fix_quadratic(a,b,c,x1,y1,x2,y2,x3,y3)
!!$
!!$          q1=a*(x1**3.)/3.+b*(x1**2.)/2.+c*x1
!!$          q2=a*(x2**3.)/3.+b*(x2**2.)/2.+c*x2
!!$          q3=a*(x3**3.)/3.+b*(x3**2.)/2.+c*x3
!!$
!!$          !Takes value for first and last sections but averages over sections where you
!!$          !have two independent estimates of the area
!!$          IF(n==3) THEN
!!$             sum=sum+q3-q1
!!$          ELSE IF(i==1) THEN
!!$             sum=sum+(q2-q1)+(q3-q2)/2.
!!$          ELSE IF(i==n-2) THEN
!!$             sum=sum+(q2-q1)/2.+(q3-q2)
!!$          ELSE
!!$             sum=sum+(q3-q1)/2.
!!$          END IF
!!$
!!$       END DO
!!$
!!$    ELSE IF(iorder==3) THEN
!!$
!!$       DO i=1,n-1
!!$
!!$          !First choose the integers used for defining cubics for each section
!!$          !First and last are different because the section does not lie in the *middle* of a cubic
!!$
!!$          IF(i==1) THEN
!!$
!!$             i1=1
!!$             i2=2
!!$             i3=3
!!$             i4=4
!!$
!!$          ELSE IF(i==n-1) THEN
!!$
!!$             i1=n-3
!!$             i2=n-2
!!$             i3=n-1
!!$             i4=n
!!$
!!$          ELSE
!!$
!!$             i1=i-1
!!$             i2=i
!!$             i3=i+1
!!$             i4=i+2
!!$
!!$          END IF
!!$
!!$          x1=x(i1)
!!$          x2=x(i2)
!!$          x3=x(i3)
!!$          x4=x(i4)
!!$
!!$          y1=y(i1)
!!$          y2=y(i2)
!!$          y3=y(i3)
!!$          y4=y(i4)
!!$
!!$          CALL fix_cubic(a,b,c,d,x1,y1,x2,y2,x3,y3,x4,y4)
!!$
!!$          !These are the limits of the particular section of integral
!!$          xi=x(i)
!!$          xf=x(i+1)
!!$
!!$          qi=a*(xi**4.)/4.+b*(xi**3.)/3.+c*(xi**2.)/2.+d*xi
!!$          qf=a*(xf**4.)/4.+b*(xf**3.)/3.+c*(xf**2.)/2.+d*xf
!!$
!!$          sum=sum+qf-qi
!!$
!!$       END DO
!!$
!!$    ELSE
!!$
!!$       STOP 'INTTAB: Error, order not specified correctly'
!!$
!!$    END IF
!!$
!!$    inttab=REAL(sum)
!!$
!!$  END FUNCTION inttab

  FUNCTION win_type(ik,itype,k,m,rv,rs,z,lut,cosm)

    IMPLICIT NONE
    REAL :: win_type
    REAL, INTENT(IN) :: k, m, rv, rs, z
    INTEGER, INTENT(IN) :: itype, ik
    TYPE(cosmology), INTENT(IN) :: cosm
    TYPE(tables), INTENT(IN) :: lut

    !IF(ik .NE. 0 .OR. ik .NE. 1) STOP 'WIN_TYPE: ik should be either 0 or 1'

    IF(itype==-1) THEN
       !Overdensity if all the matter were CDM
       win_type=win_DMONLY(ik,k,m,rv,rs,cosm)
    ELSE IF(itype==0) THEN
       !matter overdensity (sum of CDM, gas, stars)
       win_type=win_total(ik,k,m,rv,rs,cosm)
    ELSE IF(itype==1) THEN
       !CDM overdensity
       win_type=win_CDM(ik,k,m,rv,rs,cosm)
    ELSE IF(itype==2) THEN
       !All gas, both bound and free overdensity
       win_type=win_gas(ik,k,m,rv,rs,cosm)
    ELSE IF(itype==3) THEN
       !Stellar overdensity
       win_type=win_star(ik,k,m,rv,rs,cosm)
    ELSE IF(itype==4) THEN
       !Bound gas overdensity
       win_type=win_boundgas(ik,1,k,m,rv,rs,cosm)
    ELSE IF(itype==5) THEN
       !Free gas overdensity
       win_type=win_freegas(ik,1,k,m,rv,rs,cosm)
    ELSE IF(itype==6) THEN
       !Pressure
       win_type=win_pressure(ik,k,m,rv,rs,z,lut,cosm)
    ELSE IF(itype==7) THEN
       !Compensated void
       win_type=win_void(ik,k,m,rv,rs,cosm)
    ELSE IF(itype==8) THEN
       !Compensated void
       win_type=win_compensated_void(ik,k,m,rv,rs,cosm)
    ELSE
       STOP 'WIN_TYPE: Error, itype not specified correclty' 
    END IF

  END FUNCTION win_type

  FUNCTION win_total(ik,k,m,rv,rs,cosm)

    IMPLICIT NONE
    REAL :: win_total
    INTEGER, INTENT(IN) :: ik
    REAL, INTENT(IN) :: k, rv, rs, m
    TYPE(cosmology), INTENT(IN) :: cosm

    win_total=win_CDM(ik,k,m,rv,rs,cosm)+win_gas(ik,k,m,rv,rs,cosm)+win_star(ik,k,m,rv,rs,cosm)

  END FUNCTION win_total

  FUNCTION win_gas(ik,k,m,rv,rs,cosm)

    IMPLICIT NONE
    REAL :: win_gas
    INTEGER, INTENT(IN) :: ik
    REAL, INTENT(IN) :: k, m, rv, rs
    TYPE(cosmology), INTENT(IN) :: cosm

    win_gas=win_boundgas(ik,1,k,m,rv,rs,cosm)+win_freegas(ik,1,k,m,rv,rs,cosm)

  END FUNCTION win_gas

  FUNCTION win_void(ik,k,m,rv,rs,cosm)

    IMPLICIT NONE
    REAL :: win_void
    INTEGER, INTENT(IN) :: ik
    REAL, INTENT(IN) :: k, m, rv, rs
    TYPE(cosmology), INTENT(IN) :: cosm
    INTEGER :: irho
    REAL :: r, rmin, rmax

    !Set the void model
    !1 - Top-hat void
    INTEGER, PARAMETER :: imod=1

    IF(imod==1) THEN
       !Top-hat
       irho=2
    ELSE
       STOP 'WIN_VOID: Error, imod specified incorrectly'
    END IF

    rmin=0.
    rmax=10.*rv

    IF(ik==0) THEN
       r=k
       win_void=rho(r,rmin,rmax,rv,rs,irho)
       win_void=win_void/normalisation(rmin,rmax,rv,rs,irho)
    ELSE       
       win_void=m*win_norm(k,rmin,rmax,rv,rs,irho)/comoving_matter_density(cosm)
    END IF

  END FUNCTION win_void

  FUNCTION win_compensated_void(ik,k,m,rv,rs,cosm)

    IMPLICIT NONE
    REAL:: win_compensated_void
    INTEGER, INTENT(IN) :: ik
    REAL, INTENT(IN) :: k, m, rv, rs
    TYPE(cosmology), INTENT(IN) :: cosm
    INTEGER :: irho
    REAL :: r, rmin, rmax

    !Set the void model
    !1 - Top-hat void
    INTEGER, PARAMETER :: imod=1

    IF(imod==1) THEN
       !Top-hat
       irho=2
    ELSE
       STOP 'WIN_COMPENSATED_VOID: Error, imod specified incorrectly'
    END IF

    rmin=0.
    rmax=10.*rv

    IF(ik==0) THEN
       r=k
       win_compensated_void=rho(r,rmin,rmax,rv,rs,irho)
       win_compensated_void=win_compensated_void/normalisation(rmin,rmax,rv,rs,irho)
    ELSE       
       win_compensated_void=m*win_norm(k,rmin,rmax,rv,rs,irho)/comoving_matter_density(cosm)
    END IF

  END FUNCTION win_compensated_void

  FUNCTION win_DMONLY(ik,k,m,rv,rs,cosm)

    IMPLICIT NONE
    REAL :: win_DMONLY
    INTEGER, INTENT(IN) :: ik
    REAL, INTENT(IN) :: k, m, rv, rs
    TYPE(cosmology), INTENT(IN) :: cosm
    INTEGER :: irho
    REAL :: r, rmin, rmax

    !Set the DMONLY halo model
    !1 - Analyical NFW
    !2 - Non-analytical NFW (for testing W(k) functions)
    !3 - Tophat
    !4 - Delta function
    INTEGER, PARAMETER :: imod=1

    IF(imod==1) THEN
       irho=5 !Analytical NFW
    ELSE IF(imod==2) THEN
       irho=4 !Non-analyical NFW
    ELSE IF(imod==3) THEN
       irho=2 !Tophat
    ELSE IF(imod==4) THEN
       irho=0 !Delta function
    ELSE
       STOP 'WIN_DMONLY: Error, imod specified incorrectly'
    END IF

    rmin=0.
    rmax=rv

    IF(ik==0) THEN
       r=k
       win_DMONLY=rho(r,rmin,rmax,rv,rs,irho)
       win_DMONLY=win_DMONLY/normalisation(rmin,rmax,rv,rs,irho)
    ELSE IF(ik==1) THEN
       !Properly normalise and convert to overdensity
       win_DMONLY=m*win_norm(k,rmin,rmax,rv,rs,irho)/comoving_matter_density(cosm)
    ELSE
       STOP 'WIN_DMONLY: ik not specified correctly'
    END IF

  END FUNCTION win_DMONLY

  FUNCTION win_CDM(ik,k,m,rv,rs,cosm)

    IMPLICIT NONE
    REAL :: win_CDM
    INTEGER, INTENT(IN) :: ik
    REAL, INTENT(IN) :: k, m, rv, rs
    TYPE(cosmology), INTENT(IN) :: cosm
    INTEGER :: irho
    REAL :: rss, dc, r, rmin, rmax

    !Set the model
    !1 - NFW
    !2 - NFW with concentration change due to gas
    INTEGER, PARAMETER :: imod=2

    IF(imod==1) THEN
       !Analytical NFW
       irho=5
       rss=rs
    ELSE IF(imod==2) THEN
       !NFW with increase concentation
       irho=5
       !dc=1.
       dc=cosm%param(2)*halo_boundgas_fraction(m,cosm)*cosm%om_m/cosm%om_b
       rss=1./(1./rs+dc/rv)
    ELSE
       STOP 'WIN_CDM: Error, imod specified incorrectly'
    END IF

    rmin=0.
    rmax=rv

    IF(ik==0) THEN
       r=k
       win_CDM=rho(r,rmin,rmax,rv,rss,irho)
       win_CDM=win_CDM/normalisation(rmin,rmax,rv,rss,irho)
    ELSE IF(ik==1) THEN
       !Properly normalise and convert to overdensity
       win_CDM=m*win_norm(k,rmin,rmax,rv,rss,irho)/comoving_matter_density(cosm)
    ELSE
       STOP 'WIN_CDM: ik not specified correctly'
    END IF

    win_CDM=halo_CDM_fraction(m,cosm)*win_CDM

  END FUNCTION win_CDM

  FUNCTION win_star(ik,k,m,rv,rs,cosm)

    IMPLICIT NONE
    REAL :: win_star
    INTEGER, INTENT(IN) :: ik
    REAL, INTENT(IN) :: k, m, rv, rs
    TYPE(cosmology), INTENT(IN) :: cosm
    INTEGER :: irho
    REAL :: rstar, r, rmin, rmax
    REAL :: crap

    !Set the model
    !1 - Fedeli (2014) stellar distribution
    !2 - Schneider (2015) stellar distribution
    !3 - Delta function
    INTEGER, PARAMETER :: imod=1 !Set the model

    !To prevent compile-time warnings
    crap=rs

    IF(imod==1) THEN
       !Fedeli (2014)
       irho=7
       rstar=0.1*rv
       rmax=10.*rstar !Set so that not too much bigger than rstar, otherwise bumps integration goes tits
    ELSE IF(imod==2) THEN
       !Schneider (2015), following Mohammed (2014)
       irho=9
       rstar=0.01*rv
       rmax=10.*rstar !Set so that not too much bigger than rstar, otherwise bumps integration goes tits
    ELSE IF(imod==3) THEN
       !Delta function
       irho=0
       rmax=rv !Set this although it does not matter
       rstar=rv !Set this although it does not matter
    ELSE
       STOP 'WIN_STAR: Error, imod_star specified incorrectly'
    END IF

    rmin=0.

    IF(ik==0) THEN
       r=k
       win_star=rho(r,rmin,rmax,rv,rstar,irho)
       win_star=win_star/normalisation(rmin,rmax,rv,rstar,irho)
    ELSE IF(ik==1) THEN
       !Properly normalise and convert to overdensity
       win_star=m*win_norm(k,rmin,rmax,rv,rstar,irho)/comoving_matter_density(cosm)
    ELSE
       STOP 'WIN_STAR: ik not specified correctly'
    END IF

    win_star=halo_star_fraction(m,cosm)*win_star

  END FUNCTION win_star

  FUNCTION win_pressure(ik,k,m,rv,rs,z,lut,cosm)

    !Window function for bound + unbound pressures
    IMPLICIT NONE
    REAL :: win_pressure
    INTEGER, INTENT(IN) :: ik
    REAL, INTENT(IN) :: k, m, rv, rs, z
    TYPE(cosmology), INTENT(IN) :: cosm
    TYPE(tables), INTENT(IN) :: lut
    LOGICAL, PARAMETER :: use_UPP=.FALSE.

    IF(use_UPP) THEN
       !This overrides everything and just makes a UPP
       win_pressure=UPP(ik,k,m,rv,rs,z,lut,cosm)
    ELSE
       !win_pressure=win_pressure_bound(ik,k,m,rv,rs,z,lut,cosm)+win_pressure_free(ik,k,m,rv,rs,z,lut,cosm)
       win_pressure=win_boundgas(ik,2,k,m,rv,rs,cosm)+win_freegas(ik,2,k,m,rv,rs,cosm)
    END IF

  END FUNCTION win_pressure

  FUNCTION UPP(ik,k,m,rv,rs,z,lut,cosm)

    IMPLICIT NONE
    REAL :: UPP
    INTEGER, INTENT(IN) :: ik
    REAL, INTENT(IN) :: k, m, rv, rs, z
    TYPE(cosmology), INTENT(IN) :: cosm
    TYPE(tables), INTENT(IN) :: lut  
    INTEGER :: irho
    REAL :: r500c, rmin, rmax, a, r, alphap, b, m500c, E

    !Set UPP profile
    irho=14

    !Get r500 for UPP
    r500c=exp(find(log(m),log(lut%m),log(lut%r500c),lut%n,3,3,2))

    !Set the radius range for the integration
    rmin=0.
    rmax=1.*rv

    !UPP is written in terms of physical coordinates (?)
    a=scale_factor_z(z)
    IF(ik==0) THEN
       r=k
       !win_pressure_bound=rho(a*r,a*rmax,a*r500c,a*rs,irho)
       UPP=rho(r,rmin,rmax,r500c,rs,irho)
    ELSE IF(ik==1) THEN
       !win_pressure_bound=winint(k/a,a*rmax,a*r500c,a*rs,irho)
       UPP=winint(k,rmin,rmax,r500c,rs,irho)
    ELSE
       STOP 'WIN_PRESSURE_BOUND: Error, ik not specified correctly'
    END IF

    !UPP parameter
    alphap=0.12
    b=0. !How different is the inferred hydrostatic mass from true mass? (M_obs = (1-b) * M_true)

    !Upp, P(x), equation 4.1 in Ma et al. (2015)
    m500c=exp(find(log(m),log(lut%m),log(lut%m500c),lut%n,3,3,2))
    m500c=m500c*(1.-b)

    !Dimensionless Hubble
    E=sqrt(Hubble2(z,cosm))

    !WRITE(*,*) 'M [Msun/h]:', REAL(m)
    !WRITE(*,*) 'M500c [Msun/h]:', REAL(m500c)
    !WRITE(*,*) 'r500c [Mpc/h]:', r500c
    !WRITE(*,*)

    !Pre-factors from equation 4.1 in Ma et al. (2015)
    UPP=UPP*((m500c/2.1e14)**(alphap+2./3.))*(E**(8./3.))*3.37
    !win_pressure_bound=win_pressure_bound*eV*0.01**(-3) !Convert from eV cm^-3 to J m^-3

    !Is pressure comoving or not?
    !win_pressure_bound=win_pressure_bound/(1.+z)**2.

    !Convert thermal pressure to electron pressure
    UPP=UPP/pfac

  END FUNCTION UPP

  FUNCTION win_boundgas(ik,itype,k,m,rv,rs,cosm)

    !Window function for the pressure of the bound component
    IMPLICIT NONE
    REAL :: win_boundgas
    INTEGER, INTENT(IN) :: ik, itype
    REAL, INTENT(IN) :: k, m, rv, rs
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: rho0, T0, r, alpha
    REAL :: rmin, rmax, rb
    INTEGER :: irho_density, irho_pressure

    !Select model
    !1 - Komatsu-Seljak gas model
    !2 - Isothermal beta model
    INTEGER, PARAMETER :: imod=1

    IF(imod==1) THEN
       !Set KS profile
       irho_density=11
       irho_pressure=13
       rmin=0.
       rmax=rv
       rb=rs
    ELSE IF(imod==2) THEN
       irho_density=6 !Set cored isothermal profile with beta=2/3 
       irho_pressure=irho_density !okay to use density for pressure because temperature is constant
       rmin=0.
       rmax=rv
       rb=rs
    ELSE        
       STOP 'WIN_BOUNDGAD: Error, imod not specified correctly'
    END IF

    IF(itype==1) THEN

       !Density profile of bound gas
       IF(ik==0) THEN
          r=k
          win_boundgas=rho(r,rmin,rmax,rv,rb,irho_density)
          win_boundgas=win_boundgas/normalisation(rmin,rmax,rv,rb,irho_density)
       ELSE IF(ik==1) THEN
          !Properly normalise and convert to overdensity
          win_boundgas=m*win_norm(k,rmin,rmax,rv,rb,irho_density)/comoving_matter_density(cosm)
       ELSE
          STOP 'WIN_BOUNDGAS: ik not specified correctly'
       END IF

       win_boundgas=halo_boundgas_fraction(m,cosm)*win_boundgas

    ELSE IF(itype==2) THEN

       !Pressure profile of bound gas
       IF(ik==0) THEN
          r=k
          win_boundgas=rho(r,rmin,rmax,rv,rb,irho_pressure)
       ELSE IF(ik==1) THEN
          !The pressure window is T(r) x rho(r), we want unnormalised, so multiply by normalisation
          win_boundgas=win_norm(k,rmin,rmax,rv,rb,irho_pressure)*normalisation(rmin,rmax,rv,rb,irho_pressure) 
       ELSE
          STOP 'WIN_BOUNDGAS: Error, ik not specified correctly'
       END IF

       !Calculate the value of the density profile prefactor
       !also change units from cosmological to SI
       rho0=m*halo_boundgas_fraction(m,cosm)/normalisation(rmin,rmax,rv,rb,irho_density)
       rho0=rho0*msun/mpc**3 !OVERFLOW ERROR WITH REAL(4)

       !Calculate the value of the temperature prefactor
       !f=p=pac=1.
       !IF(variation) fac=param(1) !Fudge factor (turbulence?)
       alpha=cosm%param(1)
       T0=alpha*virial_temperature(m,rv)

       !Get the units correct
       win_boundgas=win_boundgas*rho0*T0*kb/(mp*mue)
       win_boundgas=win_boundgas/(eV*cm**(-3))
       win_boundgas=win_boundgas/pfac

    ELSE

       STOP 'WIN_BOUNDGAS: Error, itype not specified correctly'

    END IF

  END FUNCTION win_boundgas

  FUNCTION win_freegas(ik,itype,k,m,rv,rs,cosm)

    IMPLICIT NONE
    REAL :: win_freegas
    INTEGER, INTENT(IN) :: ik, itype
    REAL, INTENT(IN) :: k, m, rv, rs
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: rf, rmin, rmax, r, A, gamma, rho0, T0
    INTEGER :: irho_density, irho_pressure
    LOGICAL :: imatch_pressure

    !Set the model
    !1 - Isothermal model (out to 2rv)
    !2 - Ejected gas model from Schneider (2015)
    !3 - Isothermal shell that connects pressure and density to boundgas at rv
    INTEGER, PARAMETER :: imod=3

    imatch_pressure=.FALSE.

    IF(halo_freegas_fraction(m,cosm)==0.) THEN

       !Sometimes the freegas fraction will be zero, in which case this avoids problems
       win_freegas=0.

    ELSE

       IF(imod==1) THEN
          !Simple isothermal model, motivated by constant velocity and rate expulsion
          irho_density=1
          irho_pressure=irho_density !Okay because T is constant
          rmin=0.
          rmax=2.*rv
          rf=rs !Does not need to be set
       ELSE IF(imod==2) THEN
          !Ejected gas model from Schneider (2015)
          irho_density=10
          irho_pressure=irho_density !Okay because T is constant
          rmin=0.
          rf=rv
          !rf=param(3)*rv
          rmax=15.*rf !Needs to be such that integral converges (15rf seems okay)
       ELSE IF(imod==3) THEN
          !Isothermal model with continuous link to KS
          !!
          !Komatsu-Seljak density at the virial radius
          A=win_boundgas(0,1,rv,m,rv,rs,cosm) !This is A as in A/r^2
          rho0=A !This is the value of rho at the halo boundary for the gas
          A=A*rv**2 !This is A as in A/r^2
          !!
          !B=win_boundgas(0,2,rv,m,rv,rs,cosm)
          !!
          !!
          !Now do isothermal shell connected to the KS profile continuously
          irho_density=16
          irho_pressure=irho_density !Okay because T is constant
          rf=rs !Does not need to be set
          rmin=rv
          rmax=rv+halo_freegas_fraction(m,cosm)/(4.*pi*A) !This ensures density continuity and mass conservation
          gamma=5.
          IF(rmax>gamma*rv) rmax=gamma*rv !This needs to be set otherwise get huge decrement in gas power at large scales
          imatch_pressure=.TRUE. !Match the pressure at the boundary
          !!
       ELSE
          STOP 'WIN_FREEGAS: Error, imod_freegas specified incorrectly'
       END IF

       IF(itype==1) THEN

          !Density profile of free gas
          IF(ik==0) THEN
             r=k
             win_freegas=rho(r,rmin,rmax,rv,rf,irho_density)
             win_freegas=win_freegas/normalisation(rmin,rmax,rv,rf,irho_density)
          ELSE IF(ik==1) THEN
             !Properly normalise and convert to overdensity
             win_freegas=m*win_norm(k,rmin,rmax,rv,rf,irho_density)/comoving_matter_density(cosm)
          ELSE
             STOP 'WIN_FREEGAS: ik not specified correctly'
          END IF

          win_freegas=halo_freegas_fraction(m,cosm)*win_freegas

       ELSE IF(itype==2) THEN

          IF(imatch_pressure) THEN

             r=k
             IF(r>rmin .AND. r<rmax) THEN
                !Only works for isothermal profile
                win_freegas=win_boundgas(0,2,rv,m,rv,rs,cosm)*(r/rv)**(-2)
             ELSE
                win_freegas=0.
             END IF

          ELSE

             !Pressure profile of free gas
             IF(ik==0) THEN
                r=k
                win_freegas=rho(r,rmin,rmax,rv,rf,irho_pressure)
             ELSE IF(ik==1) THEN  
                win_freegas=win_norm(k,rmin,rmax,rv,rf,irho_pressure)*normalisation(rmin,rmax,rv,rf,irho_pressure)              
             ELSE
                STOP 'WIN_PRESSURE_FREE: Error, ik not specified correctly'
             END IF

             !Calculate the value of the density profile prefactor
             !and change units from cosmological to SI
             rho0=m*halo_freegas_fraction(m,cosm)/normalisation(rmin,rmax,rv,rf,irho_density)
             rho0=rho0*msun/mpc**3 !OVERFLOW ERROR WITH REAL(4)

             !Calculate the value of the temperature prefactor
             !beta=1.
             !beta=param(2)
             T0=virial_temperature(m,rv)

             !Pre factors to convert from Temp x density -> pressure (Temp x n_e)          
             win_freegas=win_freegas*rho0*T0*kb/(mp*mue)
             win_freegas=win_freegas/(eV*cm**(-3))
             win_freegas=win_freegas/pfac

          END IF

       ELSE

          STOP 'WIN_FREEGAS: Error, itype not specified correctly'

       END IF

    END IF

  END FUNCTION win_freegas

  FUNCTION virial_temperature(M,R)

    !Halo virial temperature in K
    IMPLICIT NONE
    REAL :: virial_temperature
    REAL :: M, R !Virial mass and radius

    REAL, PARAMETER :: fac=0.5 !Virial relation pre-factor (1/2, 3/2, ... ?)

    virial_temperature=fac*bigG*((m*msun)*mp*mue)/(r*mpc)
    virial_temperature=virial_temperature/kb !Convert to temperature from energy

  END FUNCTION virial_temperature

  FUNCTION win_norm(k,rmin,rmax,rv,rs,irho)

    !Calculates the normalised spherical Fourier Transform of the density profile
    !Note that this means win_norm(k->0)=1
    !and that win must be between 0 and 1
    !USE cosdef
    IMPLICIT NONE
    REAL :: win_norm
    REAL, INTENT(IN) :: rmin, rmax, k, rv, rs
    INTEGER, INTENT(IN) :: irho

    IF(k==0.) THEN

       !If called for the zero mode (e.g. for the normalisation)
       win_norm=1.

    ELSE

       IF(irho==0) THEN
          !Delta function profile is not localised in Fourier Space
          win_norm=1.
       ELSE IF(irho==1) THEN
          win_norm=wk_isothermal(k*rmax)
       ELSE IF(irho==2) THEN
          !Analytic for top hat
          win_norm=wk_tophat(k*rmax)
       ELSE IF(irho==5) THEN
          !Analytic for NFW
          win_norm=winnfw(k,rmax,rs)
       ELSE IF(irho==10) THEN
          !For ejected gas profile
          win_norm=exp(-1.5*(k*rs)**2.)
       ELSE IF(irho==16) THEN
          win_norm=wk_isothermal_2(k*rmax,k*rmin)
       ELSE
          !Numerical integral over the density profile (slower)
          win_norm=winint(k,rmin,rmax,rv,rs,irho)/normalisation(rmin,rmax,rv,rs,irho)
       END IF

    END IF

  END FUNCTION win_norm

  FUNCTION rhor2at0(irho)

    !This is the value of r^2/rho(r) at r=0. For most profiles this is zero

    IMPLICIT NONE
    REAL :: rhor2at0
    !REAL, INTENT(IN) :: rmax, rv, rs
    INTEGER, INTENT(IN) :: irho

    IF(irho==0) THEN
       STOP 'RHOR2AT0: You should not be here for a delta-function profile'
    ELSE IF(irho==1 .OR. irho==9) THEN
       !1 - Isothermal
       !9 - Stellar profile from Schneider (2015)
       rhor2at0=1.
    ELSE
       rhor2at0=0.
    END IF

  END FUNCTION rhor2at0

  FUNCTION rho(r,rmin,rmax,rv,rs,irho)

    !This is an UNNORMALISED halo profile of some sort (density, temperature, ...)

    !Types of profile
    !================
    ! 0 - Delta function at r=0
    ! 1 - Isothermal
    ! 2 - Top hat
    ! 3 - Moore (1999)
    ! 4 - NFW (1997)
    ! 5 - Analytic NFW
    ! 6 - Beta model with beta=2/3
    ! 7 - Star profile
    ! 8 - Komatsu & Seljak (2002) according to Schneider (2015)
    ! 9 - Stellar profile from Schneider (2015)
    !10 - Ejected gas profile (Schneider 2015)
    !11 - KS density
    !12 - KS temperature
    !13 - KS pressure
    !14 - Universal pressure profile
    !15 - Isothermal beta model, beta=0.86 (Ma et al. 2015)
    !16 - Isothermal shell

    IMPLICIT NONE
    REAL :: rho
    REAL, INTENT(IN) :: r, rmin, rmax, rv, rs
    INTEGER, INTENT(IN) :: irho
    REAL :: y, ct, t, c, gamma, rt, A
    REAL :: P0, c500, alpha, beta, r500
    REAL :: f1, f2

    IF(r<rmin .OR. r>rmax) THEN
       !The profile is considered to be zero outside this region
       rho=0.
    ELSE
       IF(irho==0) THEN
          !Delta function
          !STOP 'RHO: You should not be here for a delta-function profile'
          rho=0. 
       ELSE IF(irho==1 .OR. irho==16) THEN
          !Isothermal
          rho=1./r**2
       ELSE IF(irho==2) THEN
          !Top hat
          rho=1.
       ELSE IF(irho==3) THEN
          !Moore (1999)
          y=r/rs
          rho=1./((y**1.5)*(1.+y**1.5))
       ELSE IF(irho==4 .OR. irho==5) THEN
          !NFW (1997)
          y=r/rs
          rho=1./(y*(1.+y)**2)
       ELSE IF(irho==6) THEN
          !Isothermal beta model (X-ray gas; SZ profiles; beta=2/3 fixed)
          !AKA 'cored isothermal profile'
          y=r/rs
          beta=2./3.
          rho=1./((1.+y**2)**(3.*beta/2.))
       ELSE IF(irho==7) THEN
          !Stellar profile from Fedeli (2014a)
          y=r/rs
          rho=(1./y)*exp(-y)
       ELSE IF(irho==8) THEN
          !Komatsu & Seljak (2001) profile with NFW transition radius
          !VERY slow to calculate the W(k) for some reason
          !Also creates a weird upturn in P(k) that I do not think can be correct
          t=sqrt(5.)
          rt=rv/t
          y=r/rs
          c=rs/rv
          ct=c/t
          gamma=(1.+3.*ct)*log(1.+ct)/((1.+ct)*log(1.+ct)-ct)
          IF(r<=rt) THEN
             !Komatsu Seljak in the interior
             rho=(log(1.+y)/y)**gamma
          ELSE
             !NFW in the outskirts
             A=((rt/rs)*(1.+rt/rs)**2)*(log(1.+rt/rs)/(rt/rs))**gamma
             rho=A/(y*(1.+y)**2)
          END IF
       ELSE IF(irho==9) THEN
          !Stellar profile from Schneider (2015) via Mohammed (2014)    
          rho=exp(-(r/(2.*rs))**2)/r**2
          !Converting to y caused the integration to crash for some reason !?!
          !y=r/rs
          !rho=exp(-(y/2.)**2.)/y**2.
       ELSE IF(irho==10) THEN
          !Ejected gas profile from Schneider (2015)
          rho=exp(-0.5*(r/rs)**2.)
       ELSE IF(irho==11 .OR. irho==12 .OR. irho==13) THEN
          !Komatsu & Seljak (2001) profile
          !gamma=1.18 !Recommended by Rabold (2017)
          gamma=cosm%param(3)
          y=r/rs
          rho=(log(1.+y)/y)
          IF(irho==11) THEN
             !KS density profile
             rho=rho**(1./(gamma-1.))
          ELSE IF(irho==12) THEN
             !KS temperature profile
             rho=rho
          ELSE IF(irho==13) THEN
             !KS pressure profile
             rho=rho**(gamma/(gamma-1.))
          END IF
       ELSE IF(irho==14) THEN
          !UPP is in terms of r500c, not rv
          r500=rv       
          !UPP parameters from Planck V (2013) also in Ma et al. (2015)
          P0=6.41
          c500=1.81
          alpha=1.33
          beta=4.13
          gamma=0.31
          !UPP funny-P(x), equation 4.2 in Ma et al. (2015)
          f1=(c500*r/r500)**gamma
          f2=(1.+(c500*r/r500)**alpha)**((beta-gamma)/alpha)
          rho=P0/(f1*f2)
       ELSE IF(irho==15) THEN
          !Isothermal beta model
          !Parameter from Ma et al. (2015)
          beta=0.86
          rho=(1.+(r/rs)**2)**(-3.*beta/2.)
       ELSE
          STOP 'RHO: Error, irho not specified correctly'
       END IF

    END IF

  END FUNCTION rho

  FUNCTION winint(k,rmin,rmax,rv,rs,irho)

    !Calculates W(k,M)
    IMPLICIT NONE
    REAL :: winint
    REAL, INTENT(IN) :: k, rmin, rmax, rv, rs
    INTEGER, INTENT(IN) :: irho

    !REAL, PARAMETER :: acc=1d-4
    INTEGER, PARAMETER :: iorder=3 !Integration order
    INTEGER, PARAMETER :: imeth=3 !Integration method
    !imeth = 1 - normal integration
    !imeth = 2 - bumps with normal integration
    !imeth = 3 - storage integration
    !imeth = 4 - bumps with storage integration
    !imeth = 5 - linear bumps
    !imeth = 6 - cubic bumps
    !imeth = 7 - Hybrid with storage and cubic bumps

    !Bump methods go crazy with some star profiles (those that drop too fast)
    !You need to make sure that the rmax for the integration does not extend too far out

    !The hybrid method seems not to be faster for practical calculations here

    IF(imeth==1) THEN
       winint=winint_normal(rmin,rmax,k,rmin,rmax,rv,rs,irho,iorder,acc)
    ELSE IF(imeth==2 .OR. imeth==4 .OR. imeth==5 .OR. imeth==6 .OR. imeth==7) THEN
       IF(rmin .NE. 0.) STOP 'WININT: This cannot cope with rmin to rmax - probably could be fixed quickly'
       winint=winint_bumps(k,rmax,rv,rs,irho,iorder,acc,imeth)
    ELSE IF(imeth==3) THEN
       winint=winint_store(rmin,rmax,k,rmin,rmax,rv,rs,irho,iorder,acc)
    ELSE
       STOP 'WININT: Error, imeth not specified correctly'
    END IF

  END FUNCTION winint

!!$  SUBROUTINE winint_diagnostics(rmin,rmax,rv,rs,irho,outfile)
!!$
!!$    !Write out the winint integrand as a function of k
!!$    IMPLICIT NONE
!!$    REAL, INTENT(IN) :: rmin, rmax, rv, rs
!!$    INTEGER, INTENT(IN) :: irho
!!$    CHARACTER(len=256), INTENT(IN) :: outfile
!!$    INTEGER :: i, j
!!$    REAL :: r, k
!!$    REAL, ALLOCATABLE :: integrand(:)
!!$
!!$    REAL, PARAMETER :: kmin=1d-1
!!$    REAL, PARAMETER :: kmax=1d2
!!$    INTEGER, PARAMETER :: nr=256 !Number of points in r
!!$    INTEGER, PARAMETER :: nk=16 !Number of points in k
!!$
!!$    WRITE(*,*) 'WININT_DIAGNOSTICS: Doing these'
!!$    WRITE(*,*) 'WININT_DIAGNOSTICS: maximum r [Mpc/h]:', REAL(rmax)
!!$    WRITE(*,*) 'WININT_DIAGNOSTICS: virial radius [Mpc/h]:', REAL(rv)
!!$    WRITE(*,*) 'WININT_DIAGNOSTICS: scale radius [Mpc/h]:', REAL(rs)
!!$    WRITE(*,*) 'WININT_DIAGNOSTICS: concentration:', REAL(rv/rs)
!!$    WRITE(*,*) 'WININT_DIAGNOSTICS: profile number:', irho
!!$    WRITE(*,*) 'WININT_DIAGNOSTICS: outfile: ', TRIM(outfile)
!!$
!!$    ALLOCATE(integrand(nk))
!!$    
!!$    OPEN(7,file=outfile)
!!$    DO i=1,nr
!!$       r=progression(0.,rmax,i,nr)
!!$       DO j=1,nk
!!$          k=exp(progression(log(kmin),log(kmax),j,nk))
!!$          integrand(j)=winint_integrand(rmin,rmax,r,rv,rs,irho)*sinc(r*k)
!!$       END DO
!!$       WRITE(7,*) r, (integrand(j), j=1,nk)
!!$    END DO
!!$    CLOSE(7)
!!$
!!$    WRITE(*,*) 'WININT_DIAGNOSTICS: Done'
!!$    WRITE(*,*)
!!$    
!!$  END SUBROUTINE winint_diagnostics

  FUNCTION winint_normal(a,b,k,rmin,rmax,rv,rs,irho,iorder,acc)

    !Integration routine using 'normal' method to calculate the normalised halo FT
    IMPLICIT NONE
    REAL :: winint_normal
    REAL, INTENT(IN) :: k, rmin, rmax, rv, rs
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
             !r=a+(b-a)*DBLE(i-1)/DBLE(n-1)
             r=progression(a,b,i,n)
             sum=sum+weight*winint_integrand(r,rmin,rmax,rv,rs,irho)*sinc(r*k)

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

  FUNCTION winint_store(a,b,k,rmin,rmax,rv,rs,irho,iorder,acc)

    !Integrates between a and b until desired accuracy is reached
    !Stores information to reduce function calls
    IMPLICIT NONE
    REAL :: winint_store
    REAL, INTENT(IN) :: k, rmin, rmax, rv, rs, acc
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
       sum_2n=0.
       sum_n=0.
       sum_old=0.
       sum_new=0.

       DO j=1,jmax

          !Note, you need this to be 1+2**n for some integer n
          !j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
          n=1+2**(j-1)

          !Calculate the dx interval for this value of 'n'
          dx=(b-a)/REAL(n-1)

          IF(j==1) THEN

             !The first go is just the trapezium of the end points
             f1=winint_integrand(a,rmin,rmax,rv,rs,irho)*sinc(a*k)
             f2=winint_integrand(b,rmin,rmax,rv,rs,irho)*sinc(b*k)
             sum_2n=0.5*(f1+f2)*dx
             sum_new=sum_2n

          ELSE

             !Loop over only new even points to add these to the integral
             DO i=2,n,2
                x=progression(a,b,i,n)
                fx=winint_integrand(x,rmin,rmax,rv,rs,irho)*sinc(x*k)
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
             STOP 'WININT_STORE: Integration timed out'
          ELSE
             !Integral has not converged so store old sums and reset sum variables
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.
          END IF

       END DO

    END IF

  END FUNCTION winint_store

  FUNCTION winint_integrand(r,rmin,rmax,rv,rs,irho)

    !The integrand for the W(k) integral
    !Note that the sinc function is not included
    IMPLICIT NONE
    REAL :: winint_integrand
    REAL, INTENT(IN) :: r, rmin, rmax, rv, rs
    INTEGER, INTENT(IN) :: irho
    !REAL, PARAMETER :: rmin=0.
    !REAL, PARAMETER :: rmax=1d8 !Some large number

    IF(r==0.) THEN
       winint_integrand=4.*pi*rhor2at0(irho)
    ELSE
       winint_integrand=4.*pi*(r**2)*rho(r,rmin,rmax,rv,rs,irho)
    END IF

  END FUNCTION winint_integrand

  FUNCTION winint_bumps(k,rmax,rv,rs,irho,iorder,acc,imeth)

    !Integration routine to calculate the normalised halo FT
    IMPLICIT NONE
    REAL :: winint_bumps
    REAL, INTENT(IN) :: k, rmax, rv, rs
    INTEGER, INTENT(IN) :: irho
    INTEGER, INTENT(IN) :: iorder, imeth
    REAL, INTENT(IN) :: acc
    REAL :: sum, w, rn, rmin
    REAL :: r1, r2
    REAL :: a3, a2, a1, a0
    REAL :: x1, x2, x3, x4
    REAL :: y1, y2, y3, y4
    INTEGER :: i, n

    INTEGER, PARAMETER :: nlim=3 !Do the bumps approximation after this number of bumps

    !This MUST be set to zero for this routine
    rmin=0.

    !Calculate the number of nodes of sinc(k*rmax) for 0<=r<=rmax
    n=FLOOR(k*rmax/pi)

    !Set the sum variable to zero
    sum=0.

    !Integrate over each chunk between nodes separately
    DO i=0,n

       !Set the lower integration limit
       IF(k==0.) THEN
          !Special case when k=0 to avoid division by zero
          r1=0.
       ELSE
          r1=i*pi/k
       END IF

       !Set the upper integration limit
       IF(k==0. .OR. i==n) THEN
          !Special case when on last section because end is rmax, not a node!
          r2=rmax
       ELSE
          r2=(i+1)*pi/k
       END IF

       !WRITE(*,*) i, REAL(r1), REAL(r2)

       !Now do the integration along a section
       IF(imeth==2) THEN
          w=winint_normal(r1,r2,k,rmin,rmax,rv,rs,irho,iorder,acc)
       ELSE IF(k==0. .OR. imeth==4 .OR. (imeth==7 .AND. n<=nlim)) THEN
          w=winint_store(r1,r2,k,rmin,rmax,rv,rs,irho,iorder,acc)
       ELSE IF(imeth==5 .OR. imeth==6 .OR. imeth==7) THEN
          IF(i==0 .OR. i==n) THEN
             !First piece done 'normally' because otherwise /0 occurs in cubic
             !Last piece will not generally be over one full oscillation
             w=winint_store(r1,r2,k,rmin,rmax,rv,rs,irho,iorder,acc)
          ELSE
             IF(imeth==5) THEN
                !Linear approximation to integral between nodes - see notes
                !All results from the analytic integral of a linear polynomial vs. one sine oscillation
                rn=pi*(2*i+1)/(2.*k)
                w=(2./k**2)*winint_integrand(rn,rmin,rmax,rv,rs,irho)*((-1.)**i)/rn !Note there is no sinc here
             ELSE IF(imeth==6 .OR. (imeth==7 .AND. n>nlim)) THEN
                !Cubic approximation to integral between nodes - see notes
                !All results from the analytic integral of a cubic polynomial vs. one sine oscillation
                x1=r1 !Beginning
                x2=r1+1.*(r2-r1)/3. !Middle
                x3=r1+2.*(r2-r1)/3. !Middle
                x4=r2 !End
                y1=winint_integrand(x1,rmin,rmax,rv,rs,irho)/x1 !Note there is no sinc here
                y2=winint_integrand(x2,rmin,rmax,rv,rs,irho)/x2 !Note there is no sinc here
                y3=winint_integrand(x3,rmin,rmax,rv,rs,irho)/x3 !Note there is no sinc here
                y4=winint_integrand(x4,rmin,rmax,rv,rs,irho)/x4 !Note there is no sinc here
                CALL fix_cubic(a3,a2,a1,a0,x1,y1,x2,y2,x3,y3,x4,y4)
                w=-6.*a3*(r2+r1)-4.*a2
                w=w+(k**2)*(a3*(r2**3+r1**3)+a2*(r2**2+r1**2)+a1*(r2+r1)+2.*a0)
                w=w*((-1)**i)/k**4
             ELSE
                STOP 'BUMPS: Error, imeth specified incorrectly'
             END IF
          END IF
       ELSE
          STOP 'BUMPS: Error, imeth specified incorrectly'
       END IF

       !WRITE(*,*) i, REAL(r1), REAL(r2), REAL(w)

       sum=sum+w

       !Exit if the contribution to the sum is very tiny
       !This seems to be necessary to prevent crashes
       IF(ABS(w)<acc*ABS(sum)) EXIT

    END DO

    winint_bumps=sum

  END FUNCTION winint_bumps

  FUNCTION winnfw(k,rv,rs)

    !The analytic normalised (W(k=0)=1) Fourier Transform of the NFW profile
    IMPLICIT NONE
    REAL :: winnfw
    REAL, INTENT(IN) :: k, rv, rs
    REAL :: c, ks
    REAL :: si1, si2, ci1, ci2
    REAL :: p1, p2, p3
    REAL :: rmin, rmax

    c=rv/rs
    ks=k*rv/c

    si1=si(ks)
    si2=si((1.+c)*ks)
    ci1=ci(ks)
    ci2=ci((1.+c)*ks)

    p1=cos(ks)*(ci2-ci1)
    p2=sin(ks)*(si2-si1)
    p3=sin(ks*c)/(ks*(1.+c))

    winnfw=p1+p2-p3
    rmin=0.
    rmax=rv
    winnfw=4.*pi*winnfw*(rs**3.)/normalisation(rmin,rmax,rv,rs,4)

  END FUNCTION winnfw

  FUNCTION normalisation(rmin,rmax,rv,rs,irho)

    !This calculates the normalisation of a halo of concentration c
    !See your notes for details of what this means!

    !Factors of 4\pi have been *RESTORED*

    ! 0 - Delta function (M = 1)
    ! 1 - Isothermal (M = 4pi*rv)
    ! 2 - Top hat (M = (4pi/3)*rv^3)
    ! 3 - Moore (M = (8pi/3)*rv^3*ln(1+c^1.5)/c^3)
    ! 4,5 - NFW (M = 4pi*rs^3*[ln(1+c)-c/(1+c)])
    ! 6 - Beta model with beta=2/3 (M = 4*pi*rs^3*(rv/rs-atan(rv/rs)))
    ! 9 - Stellar profile (Schneider (2015)
    !10 - Ejected gas profile (Schneider 2015)
    !16 - Isothermal shell (M = 4pi*(rmax-rmin))

    IMPLICIT NONE
    REAL :: normalisation
    REAL, INTENT(IN) :: rmin, rmax, rv, rs
    INTEGER, INTENT(IN) :: irho
    REAL :: cmax, k

    IF(irho==0) THEN
       !Delta function
       normalisation=1.
    ELSE IF(irho==1 .OR. irho==16) THEN
       !Isothermal
       normalisation=4.*pi*(rmax-rmin)
    ELSE IF(irho==2) THEN
       !Top hat
       normalisation=4.*pi*(rmax**3-rmin**3)/3.
    ELSE IF(irho==3) THEN
       !Moore
       cmax=rmax/rs
       normalisation=(2./3.)*4.*pi*(rs**3)*log(1.+cmax**1.5)
    ELSE IF(irho==4 .OR. irho==5) THEN
       !NFW
       cmax=rmax/rs
       normalisation=4.*pi*(rs**3)*(log(1.+cmax)-cmax/(1.+cmax))
    ELSE IF(irho==6) THEN
       !Beta model with beta=2/3
       normalisation=4.*pi*(rs**3)*(rmax/rs-atan(rmax/rs))
    ELSE IF(irho==9) THEN
       !Stellar profile from Schneider (2015)
       !Assumed to go on to r -> infinity
       normalisation=4.*pi*rs*sqrt(pi)
    ELSE IF(irho==10) THEN
       !Ejected gas profile from Schneider (2015)
       !Assumed to go on to r -> infinity
       normalisation=4.*pi*sqrt(pi/2.)*rs**3
    ELSE
       !Otherwise need to do the integral numerically
       k=0. !k=0 gives normalisation
       normalisation=winint(k,rmin,rmax,rv,rs,irho)
    END IF

  END FUNCTION normalisation

  FUNCTION bnu(nu)

    !Bias function selection!
    IMPLICIT NONE
    REAL :: bnu
    REAL, INTENT(IN) :: nu    

    IF(imf==1) THEN
       bnu=bps(nu)
    ELSE IF(imf==2) THEN
       bnu=bst(nu)
    ELSE
       STOP 'BNU: Error, imf not specified correctly'
    END IF

  END FUNCTION bnu

  FUNCTION bps(nu)

    !Press Scheter bias
    IMPLICIT NONE
    REAL :: bps
    REAL, INTENT(IN) :: nu
    REAL :: dc

    dc=1.686

    bps=1.+(nu**2-1.)/dc

  END FUNCTION bps

  FUNCTION bst(nu)

    !Sheth Tormen bias
    IMPLICIT NONE
    REAL :: bst
    REAL, INTENT(IN) :: nu
    REAL :: p, q, dc

    p=0.3
    q=0.707
    dc=1.686

    bst=1.+(q*(nu**2)-1.+2.*p/(1.+(q*nu**2)**p))/dc

  END FUNCTION bst

  FUNCTION b2nu(nu)

    !Bias function selection!
    IMPLICIT NONE
    REAL :: b2nu
    REAL, INTENT(IN) :: nu    

    IF(imf==1) THEN
       b2nu=b2ps(nu)
    ELSE IF(imf==2) THEN
       b2nu=b2st(nu)
    ELSE
       STOP 'B2NU: Error, imf not specified correctly'
    END IF

  END FUNCTION b2nu

  FUNCTION b2ps(nu)

    !Press & Schechter second order bias
    IMPLICIT NONE
    REAL :: b2ps
    REAL, INTENT(IN) :: nu
    REAL :: p, q, dc
    REAL :: eps1, eps2, E1, E2, a2

    STOP 'B2PS: Check this very carefully'
    !I just took the ST form and set p=0 and q=1

    a2=-17./21.
    p=0.0
    q=1.0
    dc=1.686

    eps1=(q*nu**2-1.)/dc
    eps2=(q*nu**2)*(q*nu**2-3.)/dc**2
    E1=(2.*p)/(dc*(1.+(q*nu**2)**p))
    E2=((1.+2.*p)/dc+2.*eps1)*E1

    b2ps=2.*(1.+a2)*(eps1+E1)+eps2+E2

  END FUNCTION b2ps

  FUNCTION b2st(nu)

    !Sheth & Tormen second-order bias
    IMPLICIT NONE
    REAL :: b2st
    REAL, INTENT(IN) :: nu
    REAL :: p, q, dc
    REAL :: eps1, eps2, E1, E2, a2

    !Notation follows from Cooray & Sheth (2002) pp 25-26

    a2=-17./21.
    p=0.3
    q=0.707
    dc=1.686

    eps1=(q*nu**2-1.)/dc
    eps2=(q*nu**2)*(q*nu**2-3.)/dc**2
    E1=(2.*p)/(dc*(1.+(q*nu**2)**p))
    E2=((1.+2.*p)/dc+2.*eps1)*E1

    b2st=2.*(1.+a2)*(eps1+E1)+eps2+E2

  END FUNCTION b2st

  FUNCTION gnu(nu)

    !Mass function
    IMPLICIT NONE
    REAL :: gnu
    REAL, INTENT(IN) :: nu

    IF(imf==1) THEN
       gnu=gps(nu)
    ELSE IF(imf==2) THEN
       gnu=gst(nu)
    ELSE
       STOP 'GNU: Error, imf specified incorrectly'
    END IF

  END FUNCTION gnu

  FUNCTION gps(nu)

    !Press Scheter mass function!
    IMPLICIT NONE
    REAL :: gps
    REAL, INTENT(IN) :: nu

    REAL, PARAMETER :: A=0.7978846

    gps=A*exp(-(nu**2)/2.)

  END FUNCTION gps

  FUNCTION gst(nu)

    !Sheth Tormen mass function!
    !Note I use nu=dc/sigma(M) and this Sheth & Tormen (1999) use nu=(dc/sigma)^2
    !This accounts for some small differences
    IMPLICIT NONE
    REAL :: gst
    REAL, INTENT(IN) :: nu

    REAL, PARAMETER :: p=0.3
    REAL, PARAMETER :: q=0.707
    REAL, PARAMETER :: A=0.21616

    !p=0.3
    !q=0.707
    !A=0.21616

    gst=A*(1.+((q*nu*nu)**(-p)))*exp(-q*nu*nu/2.)

  END FUNCTION gst

  FUNCTION gnubnu(nu)

    !g(nu) times b(nu)
    IMPLICIT NONE
    REAL :: gnubnu
    REAL, INTENT(IN) :: nu

    gnubnu=gnu(nu)*bnu(nu)

  END FUNCTION gnubnu

  FUNCTION Hubble2(z,cosm)

    !Calculates Hubble^2 in units such that H^2(z=0)=1.
    IMPLICIT NONE
    REAL :: Hubble2
    REAL, INTENT(IN) :: z
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: a

    !a=1./(1.+z)
    a=scale_factor_z(z)

    Hubble2=cosm%om_m*(1.+z)**3+cosm%om_r*(1.+z)**4+cosm%om_v*X_de(a,cosm)+(1.-cosm%om)*(1.+z)**2

  END FUNCTION Hubble2

  FUNCTION Hubble2a4_highz(cosm)

    !Calculates Hubble^2a^4 in units such that H^2(z=0)=1.
    !This is only valid at high z, when only radiation is important
    !Makes some assumptions that DE is not important at high z
    !Need to worry if Omega_de is scaling anything like a^-4 (e.g., kinetic dominated a^-6)
    IMPLICIT NONE
    REAL :: Hubble2a4_highz
    TYPE(cosmology), INTENT(IN) :: cosm

    Hubble2a4_highz=cosm%om_r

  END FUNCTION Hubble2a4_highz

!!$  FUNCTION InverseHubble(z,cosm)
!!$
!!$    !1/H(z) in units of Mpc/h
!!$    !USE cosdef
!!$    IMPLICIT NONE
!!$    REAL :: InverseHubble
!!$    REAL, INTENT(IN) :: z
!!$    TYPE(cosmology) :: cosm
!!$
!!$    InverseHubble=conH0/sqrt(Hubble2(z,cosm))
!!$
!!$  END FUNCTION InverseHubble

  FUNCTION Omega_m(z,cosm)

    !This calculates Omega_m variations with z!
    IMPLICIT NONE
    REAL :: Omega_m
    REAL, INTENT(IN) :: z
    TYPE(cosmology), INTENT(IN) :: cosm

    Omega_m=(cosm%om_m*(1.+z)**3)/hubble2(z,cosm)

  END FUNCTION Omega_m

  FUNCTION grow(z,cosm)

    !Scale-independent growth function | g(z=0)=1
    !USE cosdef
    IMPLICIT NONE
    REAL :: grow
    REAL, INTENT(IN) :: z
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: a

    IF(z==0.) THEN
       grow=1.
    ELSE
       a=scale_factor_z(z)
       grow=find(a,cosm%a_growth,cosm%growth,cosm%ng,3,3,2)
    END IF

  END FUNCTION grow

  FUNCTION growint(a,cosm)

    !Integrates between a and b until desired accuracy is reached
    !Stores information to reduce function calls
    IMPLICIT NONE
    REAL :: growint
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: b
    INTEGER :: i, j
    INTEGER :: n
    REAL :: x, dx
    REAL :: f1, f2, fx
    DOUBLE PRECISION :: sum_n, sum_2n, sum_new, sum_old

    INTEGER, PARAMETER :: jmin=5
    INTEGER, PARAMETER :: jmax=30
    INTEGER, PARAMETER :: iorder=3   

    !Integration range for integration parameter
    !Note a -> 1
    b=1.

    IF(a==b) THEN

       !Fix the answer to zero if the integration limits are identical
       growint=exp(0.)

    ELSE

       !Reset the sum variable for the integration
       sum_2n=0.
       sum_n=0.
       sum_old=0.
       sum_new=0.

       DO j=1,jmax

          !Note, you need this to be 1+2**n for some integer n
          !j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
          n=1+2**(j-1)

          !Calculate the dx interval for this value of 'n'
          dx=(b-a)/REAL(n-1)

          IF(j==1) THEN

             !The first go is just the trapezium of the end points
             f1=growint_integrand(a,cosm)
             f2=growint_integrand(b,cosm)
             sum_2n=0.5*(f1+f2)*dx
             sum_new=sum_2n

          ELSE

             !Loop over only new even points to add these to the integral
             DO i=2,n,2
                !x=a+(b-a)*REAL(i-1)/REAL(n-1)
                x=progression(a,b,i,n)
                fx=growint_integrand(x,cosm)
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
                STOP 'GROWINT: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (ABS(-1.+sum_new/sum_old)<acc)) THEN
             !jmin avoids spurious early convergence
             growint=exp(REAL(sum_new))
             EXIT
          ELSE IF(j==jmax) THEN
             STOP 'GROWINT: Integration timed out'
          ELSE
             !Integral has not converged so store old sums and reset sum variables
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.
          END IF

       END DO

    END IF

  END FUNCTION growint

  FUNCTION growint_integrand(a,cosm)

    !Integrand for the approximate growth integral
    IMPLICIT NONE
    REAL :: growint_integrand
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: gam

    IF(cosm%w<-1.) THEN
       gam=0.55+0.02*(1.+cosm%w)
    ELSE IF(cosm%w>-1) THEN
       gam=0.55+0.05*(1.+cosm%w)
    ELSE
       gam=0.55
    END IF

    !Note the minus sign here
    growint_integrand=-(Omega_m(-1.+1./a,cosm)**gam)/a

  END FUNCTION growint_integrand

  FUNCTION dispint(R,z,cosm)

    !Integrates between a and b until desired accuracy is reached
    !Stores information to reduce function calls
    IMPLICIT NONE
    REAL :: dispint
    REAL, INTENT(IN) :: z, R
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: a, b
    INTEGER :: i, j
    INTEGER :: n
    REAL :: x, dx
    REAL :: f1, f2, fx
    DOUBLE PRECISION :: sum_n, sum_2n, sum_new, sum_old

    INTEGER, PARAMETER :: jmin=5
    INTEGER, PARAMETER :: jmax=30
    !REAL, PARAMETER :: acc=1d-4
    INTEGER, PARAMETER :: iorder=3   

    !Integration range for integration parameter
    !Note 0 -> infinity in k has changed to 0 -> 1 in x
    a=0.
    b=1.

    IF(a==b) THEN

       !Fix the answer to zero if the integration limits are identical
       dispint=0.

    ELSE

       !Reset the sum variable for the integration
       sum_2n=0.
       sum_n=0.
       sum_old=0.
       sum_new=0.

       DO j=1,jmax

          !Note, you need this to be 1+2**n for some integer n
          !j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
          n=1+2**(j-1)

          !Calculate the dx interval for this value of 'n'
          dx=(b-a)/REAL(n-1)

          IF(j==1) THEN

             !The first go is just the trapezium of the end points
             f1=dispint_integrand(a,R,z,cosm)
             f2=dispint_integrand(b,R,z,cosm)
             sum_2n=0.5*(f1+f2)*dx
             sum_new=sum_2n

          ELSE

             !Loop over only new even points to add these to the integral
             DO i=2,n,2
                !x=a+(b-a)*REAL(i-1)/REAL(n-1)
                x=progression(a,b,i,n)
                fx=dispint_integrand(x,R,z,cosm)
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
                STOP 'DISPINT: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (ABS(-1.+sum_new/sum_old)<acc)) THEN
             !jmin avoids spurious early convergence
             dispint=REAL(sum_new)
             EXIT
          ELSE IF(j==jmax) THEN
             STOP 'DISPINT: Integration timed out'
          ELSE
             !Integral has not converged so store old sums and reset sum variables
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.
          END IF

       END DO

    END IF

  END FUNCTION dispint

  FUNCTION dispint_integrand(theta,R,z,cosm)

    !This is the integrand for the velocity dispersion integral
    IMPLICIT NONE
    REAL :: dispint_integrand
    REAL, INTENT(IN) :: theta, z, R
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: k

    REAL, PARAMETER :: alpha=1.65 !Speeds up integral for large 'R'
    REAL, PARAMETER :: Rsplit=10. !Value to impliment speed up

    !Note that I have not included the speed up alpha and Rsplit
    !The choice of alpha=1.65 seemed to work well for R=100.
    !Rsplit=10 is thoughlessly chosen (only because 100.>10.)
    !Including this seems to make things slower (faster integration but slower IF statements?)

    IF(theta==0. .OR. theta==1.) THEN
       dispint_integrand=0.
    ELSE
       !IF(r>Rsplit) THEN
       !   k=(-1.+1./theta)/r**alpha
       !ELSE
       k=(-1.+1./theta)
       !END IF
       dispint_integrand=(p_lin(k,z,cosm)/k**2)*(wk_tophat(k*r)**2)/(theta*(1.-theta))
    END IF

  END FUNCTION dispint_integrand

  SUBROUTINE fill_growtab(cosm)

    !Fills a table of the growth function vs. a
    !USE cosdef
    IMPLICIT NONE
    TYPE(cosmology) :: cosm
    INTEGER :: i
    REAL :: a, norm
    REAL, ALLOCATABLE :: d_tab(:), v_tab(:), a_tab(:)
    REAL :: ainit, amax, dinit, vinit

    !REAL, PARAMETER :: acc=1d-4
    INTEGER, PARAMETER :: n=64 !Number of entries for growth tables

    !The calculation should start at a z when Om_m(z)=1., so that the assumption
    !of starting in the g\propto a growing mode is valid (this will not work for early DE)
    ainit=0.001
    !Final should be a=1. unless considering models in the future
    amax=1.

    !These set the initial conditions to be the Om_m=1. growing mode
    dinit=ainit
    vinit=1.

    IF(verbose) WRITE(*,*) 'GROWTH: Solving growth equation'
    CALL ode_growth(d_tab,v_tab,a_tab,0.,ainit,amax,dinit,vinit,acc,3,cosm)
    IF(verbose) WRITE(*,*) 'GROWTH: ODE done'

    !Normalise so that g(z=0)=1
    norm=find(1.,a_tab,d_tab,SIZE(a_tab),3,3,2)
    IF(verbose) WRITE(*,*) 'GROWTH: unnormalised g(a=1):', REAL(norm)
    d_tab=d_tab/norm
    IF(verbose) WRITE(*,*)

    !This downsamples the tables that come out of the ODE solver (which can be a bit long)
    !Could use some table-interpolation routine here to save time
    IF(ALLOCATED(cosm%a_growth)) DEALLOCATE(cosm%a_growth,cosm%growth)
    cosm%ng=n

    ALLOCATE(cosm%a_growth(n),cosm%growth(n))
    DO i=1,n
       !a=ainit+(amax-ainit)*float(i-1)/float(n-1)
       a=progression(ainit,amax,i,n)
       cosm%a_growth(i)=a
       cosm%growth(i)=find(a,a_tab,d_tab,SIZE(a_tab),3,3,2)
    END DO

  END SUBROUTINE fill_growtab

  SUBROUTINE ode_growth(x,v,t,kk,ti,tf,xi,vi,acc,imeth,cosm)

    !Solves 2nd order ODE x''(t) from ti to tf and writes out array of x, v, t values 
    IMPLICIT NONE
    REAL :: xi, ti, tf, dt, acc, vi, x4, v4, t4, kk
    REAL :: kx1, kx2, kx3, kx4, kv1, kv2, kv3, kv4
    REAL, ALLOCATABLE :: x8(:), t8(:), v8(:), xh(:), th(:), vh(:)
    REAL, ALLOCATABLE :: x(:), v(:), t(:)
    INTEGER :: i, j, k, n, np, ifail, kn, imeth
    TYPE(cosmology) :: cosm

    INTEGER, PARAMETER :: jmax=30
    INTEGER, PARAMETER :: ninit=100

    !xi and vi are the initial values of x and v (i.e. x(ti), v(ti))
    !fx is what x' is equal to
    !fv is what v' is equal to
    !acc is the desired accuracy across the entire solution
    !imeth selects method

    IF(ALLOCATED(x)) DEALLOCATE(x)
    IF(ALLOCATED(v)) DEALLOCATE(v)
    IF(ALLOCATED(t)) DEALLOCATE(t)

    DO j=1,jmax

       !Set the number of points for the forward integration
       n=ninit*(2**(j-1))
       n=n+1  

       !Allocate arrays
       ALLOCATE(x8(n),t8(n),v8(n))

       !Set the arrays to initialy be zeroes (is this neceseary?)
       x8=0.
       t8=0.
       v8=0.

       !Set the intial conditions at the intial time
       x8(1)=xi
       v8(1)=vi

       !Fill up a table for the time values
       CALL fill_array(ti,tf,t8,n)

       !Set the time interval
       dt=(tf-ti)/float(n-1)

       !Intially fix this to zero. It will change to 1 if method is a 'failure'
       ifail=0

       DO i=1,n-1

          x4=real(x8(i))
          v4=real(v8(i))
          t4=real(t8(i))

          IF(imeth==1) THEN

             !Crude method
             kx1=dt*fd(x4,v4,kk,t4,cosm)
             kv1=dt*fv(x4,v4,kk,t4,cosm)

             x8(i+1)=x8(i)+kx1
             v8(i+1)=v8(i)+kv1

          ELSE IF(imeth==2) THEN

             !Mid-point method
             !2017/06/18 - There was a bug in this part before. Luckily it was not used. Thanks Dipak Munshi.
             kx1=dt*fd(x4,v4,kk,t4,cosm)
             kv1=dt*fv(x4,v4,kk,t4,cosm)
             kx2=dt*fd(x4+kx1/2.,v4+kv1/2.,kk,t4+dt/2.,cosm)
             kv2=dt*fv(x4+kx1/2.,v4+kv1/2.,kk,t4+dt/2.,cosm)

             x8(i+1)=x8(i)+kx2
             v8(i+1)=v8(i)+kv2

          ELSE IF(imeth==3) THEN

             !4th order Runge-Kutta method (fast!)
             kx1=dt*fd(x4,v4,kk,t4,cosm)
             kv1=dt*fv(x4,v4,kk,t4,cosm)
             kx2=dt*fd(x4+kx1/2.,v4+kv1/2.,kk,t4+dt/2.,cosm)
             kv2=dt*fv(x4+kx1/2.,v4+kv1/2.,kk,t4+dt/2.,cosm)
             kx3=dt*fd(x4+kx2/2.,v4+kv2/2.,kk,t4+dt/2.,cosm)
             kv3=dt*fv(x4+kx2/2.,v4+kv2/2.,kk,t4+dt/2.,cosm)
             kx4=dt*fd(x4+kx3,v4+kv3,kk,t4+dt,cosm)
             kv4=dt*fv(x4+kx3,v4+kv3,kk,t4+dt,cosm)

             x8(i+1)=x8(i)+(kx1+(2.*kx2)+(2.*kx3)+kx4)/6.
             v8(i+1)=v8(i)+(kv1+(2.*kv2)+(2.*kv3)+kv4)/6.

          END IF

          !t8(i+1)=t8(i)+dt

       END DO

       IF(j==1) ifail=1

       IF(j .NE. 1) THEN

          np=1+(n-1)/2

          DO k=1,1+(n-1)/2

             kn=2*k-1

             IF(ifail==0) THEN

                IF(xh(k)>acc .AND. x8(kn)>acc .AND. (ABS(xh(k)/x8(kn))-1.)>acc) ifail=1
                IF(vh(k)>acc .AND. v8(kn)>acc .AND. (ABS(vh(k)/v8(kn))-1.)>acc) ifail=1

                IF(ifail==1) THEN
                   DEALLOCATE(xh,th,vh)
                   EXIT
                END IF

             END IF
          END DO

       END IF

       IF(ifail==0) THEN
          ALLOCATE(x(n),t(n),v(n))
          x=real(x8)
          v=real(v8)
          t=real(t8)
          EXIT
       END IF

       ALLOCATE(xh(n),th(n),vh(n))
       xh=x8
       vh=v8
       th=t8
       DEALLOCATE(x8,t8,v8)

    END DO

  END SUBROUTINE ode_growth

  FUNCTION fd(d,v,k,a,cosm)

    !USE cosdef
    IMPLICIT NONE
    REAL :: fd
    REAL, INTENT(IN) :: d, v, k, a
    REAL :: crap
    TYPE(cosmology), INTENT(IN) :: cosm

    !To prevent compile-time warnings
    crap=d
    crap=k
    crap=cosm%A
    crap=a

    !Needed for growth function solution
    !This is the fd in \dot{\delta}=fd

    fd=v

  END FUNCTION fd

  FUNCTION fv(d,v,k,a,cosm)

    !USE cosdef
    IMPLICIT NONE
    REAL :: fv
    REAL, INTENT(IN) :: d, v, k, a
    REAL :: f1, f2, z
    REAL :: crap
    TYPE(cosmology), INTENT(IN) :: cosm

    !To prevent compile-time warning
    crap=k

    !Needed for growth function solution
    !This is the fv in \ddot{\delta}=fv

    !z=-1.+(1./a)
    z=redshift_a(a)

    f1=3.*omega_m(z,cosm)*d/(2.*(a**2))
    f2=(2.+AH(z,cosm)/Hubble2(z,cosm))*(v/a)

    fv=f1-f2

  END FUNCTION fv

  FUNCTION AH(z,cosm)

    !USE cosdef
    IMPLICIT NONE
    REAL :: AH
    REAL, INTENT(IN) :: z
    REAL :: a
    TYPE(cosmology), INTENT(IN) :: cosm

    !\ddot{a}/a

    !a=1./(1.+z)
    a=scale_factor_z(z)

    AH=cosm%om_m*(a**(-3))+cosm%om_v*(1.+3.*w_de(a,cosm))*X_de(a,cosm)

    AH=-AH/2.

  END FUNCTION AH

  FUNCTION X_de(a,cosm)

    !USE cosdef
    IMPLICIT NONE
    REAL :: X_de
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(IN) :: cosm

    !The time evolution for Om_w for w(a) DE models
    X_de=(a**(-3.*(1.+cosm%w+cosm%wa)))*exp(-3.*cosm%wa*(1.-a))

  END FUNCTION X_de

  FUNCTION w_de(a,cosm)

    !w(a) for DE models
    !USE cosdef
    IMPLICIT NONE
    REAL :: w_de
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(IN) :: cosm

    w_de=cosm%w+(1.-a)*cosm%wa

  END FUNCTION w_de

  FUNCTION wk_isothermal(x)

    !The normlaised Fourier Transform of a top-hat
    IMPLICIT NONE
    REAL :: wk_isothermal
    REAL, INTENT(IN) :: x

    REAL, PARAMETER :: dx=1e-3

    !Taylor expansion used for low |x| to avoid cancellation problems

    IF(ABS(x)<ABS(dx)) THEN
       !Taylor series at low x
       wk_isothermal=1.-(x**2)/18.
    ELSE
       wk_isothermal=Si(x)/x
    END IF

  END FUNCTION wk_isothermal

  FUNCTION wk_isothermal_2(x,y)

    !The normlaised Fourier Transform of a top-hat
    IMPLICIT NONE
    REAL :: wk_isothermal_2
    REAL, INTENT(IN) :: x, y

    REAL, PARAMETER :: dx=1e-3

    !Taylor expansion used for low |x| to avoid cancellation problems

    !IF(ABS(x)<ABS(dx)) THEN
    !   !Taylor series at low x
    !   wk_isothermal_2=1.-(x**2)/18.
    !ELSE
    wk_isothermal_2=(Si(x)-Si(y))/(x-y)
    !END IF

  END FUNCTION wk_isothermal_2

  FUNCTION halo_fraction(itype,m,cosm)

    !Mass fraction of a type within a halo
    IMPLICIT NONE
    REAL :: halo_fraction
    INTEGER, INTENT(IN) :: itype
    REAL, INTENT(IN) :: m
    TYPE(cosmology) :: cosm

    If(itype==-1 .OR. itype==0) THEN
       halo_fraction=1.
    ELSE IF(itype==1) THEN
       halo_fraction=halo_CDM_fraction(m,cosm)
    ELSE IF(itype==2) THEN
       halo_fraction=halo_gas_fraction(m,cosm)
    ELSE IF(itype==3) THEN
       halo_fraction=halo_star_fraction(m,cosm)
    ELSE IF(itype==4) THEN
       halo_fraction=halo_boundgas_fraction(m,cosm)
    ELSE IF(itype==5) THEN
       halo_fraction=halo_freegas_fraction(m,cosm)
    ELSE
       STOP 'HALO_FRACTION: Error, itype not specified correcntly'
    END IF

  END FUNCTION halo_fraction

  FUNCTION halo_gas_fraction(m,cosm)

    !Mass fraction of a halo in gas
    IMPLICIT NONE
    REAL :: halo_gas_fraction
    REAL, INTENT(IN) :: m
    TYPE(cosmology), INTENT(IN) :: cosm

    halo_gas_fraction=halo_boundgas_fraction(m,cosm)+halo_freegas_fraction(m,cosm)

  END FUNCTION halo_gas_fraction

  FUNCTION halo_CDM_fraction(m,cosm)

    !Mass fraction of a halo in CDM
    IMPLICIT NONE
    REAL :: halo_CDM_fraction
    REAL, INTENT(IN) :: m
    REAL :: crap
    TYPE(cosmology), INTENT(IN) :: cosm

    !To prevent compile-time warning
    crap=m

    !Always the universal value
    halo_CDM_fraction=cosm%om_c/cosm%om_m

  END FUNCTION halo_CDM_fraction

  FUNCTION halo_freegas_fraction(m,cosm)

    !Mass fraction of a halo in free gas
    IMPLICIT NONE
    REAL :: halo_freegas_fraction
    REAL, INTENT(IN) :: m
    TYPE(cosmology), INTENT(IN) :: cosm

    !This is always all the gas that is not bound or in stars
    halo_freegas_fraction=cosm%om_b/cosm%om_m-halo_star_fraction(m,cosm)-halo_boundgas_fraction(m,cosm)
    IF(halo_freegas_fraction<0.) halo_freegas_fraction=0.

  END FUNCTION halo_freegas_fraction

  FUNCTION halo_boundgas_fraction(m,cosm)

    !Fraction of a halo in bound gas
    IMPLICIT NONE
    REAL :: halo_boundgas_fraction
    REAL, INTENT(IN) :: m
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: m0, sigma, beta

    !Set the model
    !1 - Fedeli (2014a) bound gas model
    !2 - Schneider (2015) bound gas
    !3 - Universal baryon fraction
    INTEGER, PARAMETER :: imod=2 !Set the model

    IF(imod==1) THEN
       !From Fedeli (2014a)
       m0=1.e12
       sigma=3.
       IF(m<m0) THEN
          halo_boundgas_fraction=0.
       ELSE
          halo_boundgas_fraction=erf(log10(m/m0)/sigma)*cosm%om_b/cosm%om_m
       END IF
    ELSE IF(imod==2) THEN
       !From Schneider (2015)
       !m0=1.2d14
       m0=cosm%param(4)
       beta=0.6
       halo_boundgas_fraction=(cosm%om_b/cosm%om_m)/(1.+(m0/m)**beta)
    ELSE IF(imod==3) THEN
       !Universal baryon fraction model (account for stellar contribution)
       halo_boundgas_fraction=cosm%om_b/cosm%om_m-halo_star_fraction(m,cosm)
    ELSE
       STOP 'HALO_BOUNDGAS_FRACTION: Error, imod_boundfrac not specified correctly'
    END IF

  END FUNCTION halo_boundgas_fraction

  FUNCTION halo_star_fraction(m,cosm)

    !Mass fraction of a halo in stars
    IMPLICIT NONE
    REAL :: halo_star_fraction
    REAL, INTENT(IN) :: m
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: m0, sigma, A, min
    REAL :: crap

    !Set the model
    !1 - Fedeli (2014)
    !2 - Constant stellar fraction
    !3 - Fedeli (2014) but saturates at high mass
    INTEGER, PARAMETER :: imod=3

    !To prevent compile-time warning
    crap=cosm%A

    IF(imod==1 .OR. imod==3) THEN
       !Fedeli (2014)
       !A=0.02
       !IF(variation) A=param(5)
       A=cosm%param(5)
       m0=5.e12
       sigma=1.2
       halo_star_fraction=A*exp(-((log10(m/m0))**2)/(2.*sigma**2))
       IF(imod==3) THEN
          !Suggested by Ian, the relation I have is for the central stellar mass
          !in reality this saturates for high-mass haloes (due to satellite contribution)
          min=0.01
          IF(halo_star_fraction<min .AND. m>m0) halo_star_fraction=min
       END IF
    ELSE IF(imod==2) THEN
       !Constant star fraction
       A=0.005
       halo_star_fraction=A
    ELSE
       STOP 'HALO_STAR_FRACTION: Error, imod_starfrac specified incorrectly'
    END IF

  END FUNCTION halo_star_fraction

  SUBROUTINE get_nz(ix,lens)

    !The the n(z) function for lensing
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ix
    TYPE(lensing), INTENT(INOUT) :: lens

    !CHARACTER(len=256) :: names(7)
    !names(1)='1 - RCSLenS'
    !names(2)='2 - KiDS (z = 0.1 -> 0.9)'
    !names(3)='3 - KiDS (z = 0.1 -> 0.3)'
    !names(4)='4 - KiDS (z = 0.3 -> 0.5)'
    !names(5)='5 - KiDS (z = 0.5 -> 0.7)'
    !names(6)='6 - KiDS (z = 0.7 -> 0.9)'
    !names(7)='7 - CFHTLenS (Van Waerbeke 2013)'

    !IF(inz==-1) THEN
    !   WRITE(*,*) 'GET_NZ: Choose n(z)'
    !   WRITE(*,*) '==================='
    !   DO i=1,SIZE(names)
    !      WRITE(*,*) TRIM(names(i))
    !   END DO
    !   READ(*,*) inz
    !   WRITE(*,*) '==================='
    !   WRITE(*,*)
    !END IF

    IF(ix==1 .OR. ix==4) THEN
       CALL fill_analytic_nz_table(ix,lens)
    ELSE
       CALL fill_nz_table(ix,lens)
    END IF

    !WRITE(*,*) 'GET_NZ: ', TRIM(names(inz))
    WRITE(*,*) 'GET_NZ: zmin:', lens%z_nz(1)
    WRITE(*,*) 'GET_NZ: zmax:', lens%z_nz(lens%nnz)
    WRITE(*,*) 'GET_NZ: nz:', lens%nnz
    WRITE(*,*)

  END SUBROUTINE get_nz

  SUBROUTINE fill_analytic_nz_table(ix,lens)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ix
    TYPE(lensing), INTENT(INOUT) :: lens
    INTEGER :: i

    REAL, PARAMETER :: zmin=0.
    REAL, PARAMETER :: zmax=2.5
    INTEGER, PARAMETER :: n=128

    !From analytical function
    lens%nnz=n
    IF(ALLOCATED(lens%z_nz)) DEALLOCATE(lens%z_nz)
    IF(ALLOCATED(lens%nz)) DEALLOCATE(lens%nz)
    ALLOCATE(lens%z_nz(lens%nnz),lens%nz(lens%nnz))

    !Fill the look-up tables
    CALL fill_array(zmin,zmax,lens%z_nz,lens%nnz)
    DO i=1,n
       lens%nz(i)=nz_lensing(lens%z_nz(i),ix)
    END DO

  END SUBROUTINE fill_analytic_nz_table

  SUBROUTINE fill_nz_table(ix,lens)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ix
    TYPE(lensing), INTENT(INOUT) :: lens
    INTEGER :: i
    CHARACTER(len=256) :: input

    !Get file name
    IF(ix==5) THEN
       input='/Users/Mead/Physics/KiDS/nz/KiDS_z0.1-0.9_MEAD.txt'
    ELSE IF(ix==6) THEN
       input='/Users/Mead/Physics/KiDS/nz/KiDS_z0.1-0.3.txt'
    ELSE IF(ix==7) THEN
       input='/Users/Mead/Physics/KiDS/nz/KiDS_z0.3-0.5.txt'
    ELSE IF(ix==8) THEN
       input='/Users/Mead/Physics/KiDS/nz/KiDS_z0.5-0.7.txt'
    ELSE IF(ix==9) THEN
       input='/Users/Mead/Physics/KiDS/nz/KiDS_z0.7-0.9.txt'
    ELSE
       STOP 'GET_NZ: ix not specified correctly'
    END IF
    WRITE(*,*) 'GET_NZ: Input file:', TRIM(input)

    !Allocate arrays
    lens%nnz=count_number_of_lines(input)
    IF(ALLOCATED(lens%z_nz)) DEALLOCATE(lens%z_nz)
    IF(ALLOCATED(lens%nz))   DEALLOCATE(lens%nz)
    ALLOCATE(lens%z_nz(lens%nnz),lens%nz(lens%nnz))

    !Read in n(z) table
    OPEN(7,file=input)
    DO i=1,lens%nnz
       READ(7,*) lens%z_nz(i), lens%nz(i)
    END DO
    CLOSE(7)

  END SUBROUTINE fill_nz_table

  FUNCTION nz_lensing(z,ix)

    IMPLICIT NONE
    REAL :: nz_lensing
    REAL, INTENT(IN) :: z
    INTEGER, INTENT(IN) :: ix
    REAL :: a, b, c, d, e, f, g, h, i
    REAL :: n1, n2, n3

    IF(ix==1) THEN
       !RCSLenS
       a=2.94
       b=-0.44
       c=1.03
       d=1.58
       e=0.40
       f=0.25
       g=0.38
       h=0.81
       i=0.12
       n1=a*z*exp(-(z-b)**2/c**2)
       n2=d*z*exp(-(z-e)**2/f**2)
       n3=g*z*exp(-(z-h)**2/i**2)
       nz_lensing=n1+n2+n3
    ELSE IF(ix==4) THEN
       !CFHTLenS
       z1=0.7 !Not a free parameter in Van Waerbeke 2013
       z2=1.2 !Not a free parameter in Van Waerbeke 2013
       a=1.50
       b=0.32
       c=0.20
       d=0.46
       nz_lensing=a*exp(-((z-z1)/b)**2)+c*exp(-((z-z2)/d)**2)
    ELSE
       STOP 'NZ_LENSING: ix specified incorrectly'
    END IF

  END FUNCTION nz_lensing

  FUNCTION integrate_distance(a,b,acc,iorder,cosm)

    !Integrates between a and b until desired accuracy is reached
    !Stores information to reduce function calls
    IMPLICIT NONE
    REAL :: integrate_distance
    REAL, INTENT(IN) :: a, b, acc
    INTEGER, INTENT(IN) :: iorder
    TYPE(cosmology), INTENT(IN) :: cosm
    INTEGER :: i, j
    INTEGER :: n
    REAL :: x, dx
    REAL :: f1, f2, fx
    DOUBLE PRECISION :: sum_n, sum_2n, sum_new, sum_old

    INTEGER, PARAMETER :: jmin=5
    INTEGER, PARAMETER :: jmax=30

    IF(a==b) THEN

       !Fix the answer to zero if the integration limits are identical
       integrate_distance=0.

    ELSE

       !Set the sum variable for the integration
       sum_2n=0.
       sum_n=0.
       sum_old=0.
       sum_new=0.

       DO j=1,jmax

          !Note, you need this to be 1+2**n for some integer n
          !j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
          n=1+2**(j-1)

          !Calculate the dx interval for this value of 'n'
          dx=(b-a)/REAL(n-1)

          IF(j==1) THEN

             !The first go is just the trapezium of the end points
             f1=distance_integrand(a,cosm)
             f2=distance_integrand(b,cosm)
             sum_2n=0.5*(f1+f2)*dx
             sum_new=sum_2n

          ELSE

             !Loop over only new even points to add these to the integral
             DO i=2,n,2
                x=progression(a,b,i,n)
                fx=distance_integrand(x,cosm)
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
                STOP 'INTEGRATE_DISTANCE: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (ABS(-1.+sum_new/sum_old)<acc)) THEN
             !jmin avoids spurious early convergence
             integrate_distance=REAL(sum_new)
             EXIT
          ELSE IF(j==jmax) THEN
             STOP 'INTEGRATE_DISTANCE: Integration timed out'
          ELSE
             !Integral has not converged so store old sums and reset sum variables
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.
          END IF

       END DO

    END IF

  END FUNCTION integrate_distance

  FUNCTION distance_integrand(a,cosm)

    IMPLICIT NONE
    REAL :: distance_integrand
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: z

    REAL, PARAMETER :: amin=1e-5

    IF(a<amin) THEN
       distance_integrand=conH0/sqrt(Hubble2a4_highz(cosm))
    ELSE
       z=redshift_a(a)
       distance_integrand=conH0/(sqrt(Hubble2(z,cosm))*a**2)
    END IF

  END FUNCTION distance_integrand

  FUNCTION integrate_q(r,a,b,acc,iorder,lens,cosm)

    !Integrates between a and b until desired accuracy is reached
    !Stores information to reduce function calls
    IMPLICIT NONE
    REAL :: integrate_q
    REAL, INTENT(IN) :: a, b, r, acc
    INTEGER, INTENT(IN) :: iorder
    TYPE(cosmology), INTENT(IN) :: cosm
    TYPE(lensing), INTENT(IN) :: lens
    INTEGER :: i, j
    INTEGER :: n
    REAL :: x, dx
    REAL :: f1, f2, fx
    DOUBLE PRECISION :: sum_n, sum_2n, sum_new, sum_old

    INTEGER, PARAMETER :: jmin=5
    INTEGER, PARAMETER :: jmax=30

    IF(a==b) THEN

       !Fix the answer to zero if the integration limits are identical
       integrate_q=0.

    ELSE

       !Set the sum variable for the integration
       sum_2n=0.
       sum_n=0.
       sum_old=0.
       sum_new=0.

       DO j=1,jmax

          !Note, you need this to be 1+2**n for some integer n
          !j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
          n=1+2**(j-1)

          !Calculate the dx interval for this value of 'n'
          dx=(b-a)/REAL(n-1)

          IF(j==1) THEN

             !The first go is just the trapezium of the end points
             f1=q_integrand(a,r,lens,cosm)
             f2=q_integrand(b,r,lens,cosm)
             sum_2n=0.5*(f1+f2)*dx
             sum_new=sum_2n

          ELSE

             !Loop over only new even points to add these to the integral
             DO i=2,n,2
                !x=a+(b-a)*DBLE(i-1)/DBLE(n-1)
                x=progression(a,b,i,n)
                fx=q_integrand(x,r,lens,cosm)
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
                STOP 'INTEGRATE_Q: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (ABS(-1.+sum_new/sum_old)<acc)) THEN
             !jmin avoids spurious early convergence
             integrate_q=REAL(sum_new)
             EXIT
          ELSE IF(j==jmax) THEN
             STOP 'INTEGRATE_Q: Integration timed out'
          ELSE
             !Integral has not converged so store old sums and reset sum variables
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.
          END IF

       END DO

    END IF

  END FUNCTION integrate_q

  FUNCTION q_integrand(z,r,lens,cosm)

    !The lensing efficiency integrand, which is a function of z
    !z is integrated over while r is just a parameter
    !This is only called for n(z)
    IMPLICIT NONE
    REAL :: q_integrand
    REAL, INTENT(IN) :: r, z
    TYPE(cosmology), INTENT(IN) :: cosm
    TYPE(lensing), INTENT(IN) :: lens
    REAL :: rdash, nz

    IF(z==0.) THEN
       q_integrand=0.
    ELSE
       !Find the r'(z) variable that is integrated over
       !rdash=find(z,cosm%z_r,cosm%r,cosm%nr,3,3,2)
       rdash=cosmic_distance(z,cosm)
       !Find the n(z)
       nz=find(z,lens%z_nz,lens%nz,lens%nnz,3,3,2)
       !This is then the integrand
       q_integrand=nz*f_k(rdash-r,cosm)/f_k(rdash,cosm)
    END IF

  END FUNCTION q_integrand

  FUNCTION integrate_Limber(l,a,b,logktab,logatab,logptab,nk,na,acc,iorder,proj,cosm)

    !Integrates between a and b until desired accuracy is reached
    !Stores information to reduce function calls
    IMPLICIT NONE
    REAL :: integrate_Limber
    REAL, INTENT(IN) :: a, b, acc
    INTEGER, INTENT(IN) :: iorder
    REAL, INTENT(IN) :: l
    REAL, INTENT(IN) :: logktab(nk), logatab(na), logptab(nk,na)
    INTEGER, INTENT(IN) :: nk, na
    TYPE(projection), INTENT(IN) :: proj(2)
    TYPE(cosmology), INTENT(IN) :: cosm
    INTEGER :: i, j
    INTEGER :: n
    REAL :: x, dx
    REAL :: f1, f2, fx
    DOUBLE PRECISION :: sum_n, sum_2n, sum_new, sum_old

    INTEGER, PARAMETER :: jmin=5
    INTEGER, PARAMETER :: jmax=25

    IF(a==b) THEN

       !Fix the answer to zero if the integration limits are identical
       integrate_Limber=0.

    ELSE

       !Set the sum variable for the integration
       sum_2n=0.
       sum_n=0.
       sum_old=0.
       sum_new=0.

       DO j=1,jmax

          !Note, you need this to be 1+2**n for some integer n
          !j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
          n=1+2**(j-1)

          !Calculate the dx interval for this value of 'n'
          dx=(b-a)/REAL(n-1)

          !WRITE(*,*) j, n

          IF(j==1) THEN

             !The first go is just the trapezium of the end points
             f1=Limber_integrand(a,l,logktab,logatab,logptab,nk,na,proj,cosm)
             f2=Limber_integrand(b,l,logktab,logatab,logptab,nk,na,proj,cosm)
             sum_2n=0.5*(f1+f2)*dx
             sum_new=sum_2n

          ELSE

             !Loop over only new even points to add these to the integral
             DO i=2,n,2
                !x=a+(b-a)*DBLE(i-1)/DBLE(n-1)
                x=progression(a,b,i,n)
                fx=Limber_integrand(x,l,logktab,logatab,logptab,nk,na,proj,cosm)
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
                STOP 'INTEGRATE_LIMBER: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (ABS(-1.+sum_new/sum_old)<acc)) THEN
             !jmin avoids spurious early convergence
             integrate_Limber=REAL(sum_new)
             EXIT
          ELSE IF(j==jmax) THEN
             STOP 'INTEGRATE_LIMBER: Integration timed out'
          ELSE
             !Integral has not converged so store old sums and reset sum variables
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.
          END IF

       END DO

    END IF

  END FUNCTION integrate_Limber

  FUNCTION Limber_integrand(r,l,logktab,logatab,logptab,nk,na,proj,cosm)

    !The integrand for the Limber integral
    IMPLICIT NONE
    REAL :: Limber_integrand
    REAL, INTENT(IN) :: r, l
    REAL, INTENT(IN) :: logktab(nk), logatab(na), logptab(nk,na)
    INTEGER, INTENT(IN) :: nk, na
    TYPE(cosmology), INTENT(IN) :: cosm
    TYPE(projection), INTENT(IN) :: proj(2)
    REAL :: z, a, k, X(2)
    INTEGER :: i

    IF(r==0.) THEN

       Limber_integrand=0.

    ELSE

       !Get the two kernels
       DO i=1,2
          X(i)=find(r,proj(i)%r_X,proj(i)%X,proj(i)%nX,3,3,2)
       END DO
       !x2=find(r,proj%r_x2,proj%x2,proj%nx2,3,3,2)
       !x1=exp(find(log(r),log(proj%r_x1),log(proj%x1),proj%nx1,3,3,2)) !Barfed with this
       !x2=exp(find(log(r),log(proj%r_x2),log(proj%x2),proj%nx2,3,3,2)) !Barfed with this

       !Get variables r, z(r) and k(r) for P(k,z)
       z=redshift_r(r,cosm)
       a=scale_factor_z(z)
       k=(l+lcorr)/f_k(r,cosm) !LoVerde et al. (2008) Limber correction

       !Construct the integrand
       Limber_integrand=X(1)*X(2)*find_pka(k,a,logktab,logatab,logptab,nk,na)/f_k(r,cosm)**2
       !Limber_integrand=x1*x2/f_k(r,cosm)**2 !Note that this is very hard to integrate for kappa-y or y-y
       !Limber_integrand=find_pkz(k,z,logktab,ztab,logptab,nk,nz)/f_k(r,cosm)

    END IF

  END FUNCTION Limber_integrand

  SUBROUTINE Limber_contribution(l,a,b,logktab,logatab,logptab,nk,na,proj,cosm,outfile)

    IMPLICIT NONE
    REAL, INTENT(IN) :: l, a, b
    REAL, INTENT(IN) :: logktab(nk), logatab(na), logptab(nk,na)
    INTEGER, INTENT(IN) :: nk, na
    TYPE(cosmology), INTENT(IN) :: cosm
    TYPE(projection), INTENT(IN) :: proj(2)
    CHARACTER(len=256), INTENT(IN) :: outfile
    REAL :: k, r, z, total, int
    INTEGER :: i

    INTEGER, PARAMETER  :: n=1024 !Number of samples to take in r

    !Calculate the integral for this value of ell
    total=integrate_Limber(l,a,b,logktab,logatab,logptab,nk,na,acc,3,proj,cosm)

    !Now split up the contributions
    !You need the Jacobian and to remember that the contribution is split in ln(k), ln(z) and ln(R)
    !This means factors of:
    !k - f_k(r)/f'_k(r) (=r if flat)
    !z - z/H(z)
    !r - r
    OPEN(7,file=TRIM(outfile))
    DO i=1,n
       !r=a+(b-a)*DBLE(i-1)/DBLE(n-1)
       r=progression(a,b,i,n)
       IF(r==0.) THEN
          CYCLE
       ELSE
          k=(l+lcorr)/f_k(r,cosm)
          !z=find(r,cosm%r,cosm%z_r,cosm%nr,3,3,2)
          z=redshift_r(r,cosm)
          int=Limber_integrand(r,l,logktab,logatab,logptab,nk,na,proj,cosm)
          WRITE(7,*) k, int*f_k(r,cosm)/(fdash_k(r,cosm)*total), z, int*z/(sqrt(Hubble2(z,cosm))*total), r, int*r/total
       END IF
    END DO
    CLOSE(7)

  END SUBROUTINE Limber_contribution

  FUNCTION find_pka(k,a,logktab,logatab,logptab,nk,na)

    !Looks up the power as a 2D function of k and a
    IMPLICIT NONE
    REAL :: find_pka
    INTEGER, INTENT(IN) :: nk, na
    REAL, INTENT(IN) :: k, a
    REAL, INTENT(IN) :: logktab(nk), logatab(na), logptab(nk,na)

    REAL, PARAMETER :: kmin=1e-3 !kmin value
    REAL, PARAMETER :: kmax=1e2 !kmax value

    IF(k<=kmin .OR. k>=kmax) THEN
       find_pka=0.
    ELSE
       find_pka=exp(find2d(log(k),logktab,log(a),logatab,logptab,nk,na,3,3,1))
    END IF

    !Convert from Delta^2 -> P(k) - with dimensions of (Mpc/h)^3
    !find_pkz=(2.*pi**2)*find_pkz/k**3

  END FUNCTION find_pka

  SUBROUTINE print_baryon_parameters(cosm)!param,param_names,n)

    IMPLICIT NONE
    !REAL, INTENT(IN) :: param(n)
    !CHARACTER(len=256), INTENT(IN) :: param_names(n)
    !INTEGER, INTENT(IN) :: n
    INTEGER :: i
    TYPE(cosmology), INTENT(IN) :: cosm

    WRITE(*,*) 'PRINT_BARYON_PARAMETERS: Writing to screen'
    WRITE(*,*) '=========================================='
    DO i=1,cosm%np
       WRITE(*,*) i, TRIM(cosm%param_names(i)), ':', cosm%param(i)
    END DO
    WRITE(*,*) '=========================================='
    WRITE(*,*) 'PRINT_BARYON_PARAMETERS: Done'
    WRITE(*,*)

  END SUBROUTINE print_baryon_parameters

END MODULE HMx

PROGRAM HMx_driver

  USE HMx
  USE cosdef
  IMPLICIT NONE

  CALL init_HMx()

  CALL get_command_argument(1,mode)
  IF(mode=='') THEN
     imode=-1
  ELSE
     READ(mode,*) imode
  END IF

  !HMx developed by Alexander Mead
  WRITE(*,*)
  WRITE(*,*) 'HMx: Welcome to HMx'
  IF(imead==-1) THEN
     WRITE(*,*) 'HMx: Doing basic halo-model calculation (Two-halo term is linear)'
  ELSE IF(imead==0) THEN
     WRITE(*,*) 'HMx: Doing standard halo-model calculation (Seljak 2000)'
  ELSE IF(imead==1) THEN
     WRITE(*,*) 'HMx: Doing accurate halo-model calculation (Mead et al. 2015)'
  ELSE
     STOP 'HMx: imead specified incorrectly'
  END IF
  WRITE(*,*)

  IF(imode==-1) THEN
     WRITE(*,*) 'HMx: Choose what to do'
     WRITE(*,*) '======================'
     WRITE(*,*) ' 0 - Matter power spectrum at z = 0'
     WRITE(*,*) ' 1 - Matter power spectrum over multiple z'
     WRITE(*,*) ' 2 - Comparison with cosmo-OWLS'
     WRITE(*,*) ' 3 - Run diagnostics'
     WRITE(*,*) ' 4 - Do random cosmologies for bug testing'
     WRITE(*,*) ' 5 - Pressure field comparison'
     WRITE(*,*) ' 6 - n(z) check'
     WRITE(*,*) ' 7 - Do cross correlation'
     WRITE(*,*) ' 8 - Cross correlation as a function of cosmology'
     WRITE(*,*) ' 9 - Breakdown correlations in halo mass'
     WRITE(*,*) '10 - Breakdown correlations in redshift'
     WRITE(*,*) '11 - Breakdown correlations in halo radius'
     WRITE(*,*) '12 - Project triad'
     WRITE(*,*) '13 - Cross-correlation coefficient'
     WRITE(*,*) '14 - 3D spectra as HMx parameters vary'
     !WRITE(*,*) '15 - matter spectra as HMx parameters vary'
     READ(*,*) imode
     WRITE(*,*) '======================'
     WRITE(*,*)
  END IF

  IF(imode==0) THEN

     !Set number of k points and k range (log spaced)
     nk=200
     kmin=1e-3
     kmax=1e2
     CALL fill_array(log(kmin),log(kmax),k,nk)
     k=exp(k)
     ALLOCATE(pow_lin(nk),pow_2h(nk),pow_1h(nk),pow_full(nk))

     !Assigns the cosmological model
     icosmo=0
     CALL assign_cosmology(icosmo,cosm)

     !Normalises power spectrum (via sigma_8) and fills sigma(R) look-up tables
     CALL initialise_cosmology(cosm)
     CALL print_cosmology(cosm)

     !Sets the redshift
     z=0.

     !Initiliasation for the halomodel calcualtion
     CALL halomod_init(mmin,mmax,z,lut,cosm)

     !Do the halo-model calculation
     CALL calculate_halomod(-1,-1,k,nk,z,pow_lin,pow_2h,pow_1h,pow_full,lut,cosm)

     !na=1
     !ALLOCATE(a(na))
     !a=scale_factor_z(z)
     !a=1.

     !ip=-1
     !CALL calculate_halomod(ip,k,nk,a,na,powa_lin,powa_2h,powa_1h,powa_full,cosm)

     !Write out the answer
     outfile='data/power.dat'
     !CALL write_power(k,powa_lin(:,1),powa_2h(:,1),powa_1h(:,1),powa_full(:,1),nk,outfile)
     CALL write_power(k,pow_lin,pow_2h,pow_1h,pow_full,nk,outfile)

     !Write the one-void term if necessary
     IF(void) THEN
        OPEN(8,file='data/power_1void.dat')
        DO i=1,nk     
           WRITE(8,*) k(i), p_1v(k(i),lut)!,cosm)
        END DO
        CLOSE(8)
     END IF

  ELSE IF(imode==1) THEN

     !Assigns the cosmological model
     icosmo=0
     CALL assign_cosmology(icosmo,cosm)

     !Normalises power spectrum (via sigma_8) and fills sigma(R) look-up tables
     CALL initialise_cosmology(cosm)
     CALL print_cosmology(cosm)

     !Set number of k points and k range (log spaced)
     nk=200
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

     !Instead could set an a-range
     !amin=scale_factor(cosm%z_cmb)
     !amax=1.
     !na=16
     !CALL fill_array(amin,amax,a,na)

     !Allocate power arrays
     !ALLOCATE(powa_lin(nk,na),powa_2h(nk,na),powa_1h(nk,na),powa_full(nk,na))

     !Do the halo-model calculation
     !WRITE(*,*) 'HMx: Doing calculation'
     !DO j=na,1,-1
     !   z=redshift_a(a(j))
     !   CALL halomod_init(mmin,mmax,z,lut,cosm)
     !   WRITE(*,fmt='(A5,I5,F10.2)') 'HMx:', j, REAL(z)
     !   CALL calculate_halomod(-1,-1,k,nk,z,powa_lin(:,j),powa_2h(:,j),powa_1h(:,j),powa_full(:,j),lut,cosm)
     !END DO
     !WRITE(*,*)

     ip=-1
     CALL calculate_HMx(ip,k,nk,a,na,powa_lin,powa_2h,powa_1h,powa_full,cosm)

     base='data/power'
     CALL write_power_a_multiple(k,a,powa_lin,powa_2h,powa_1h,powa_full,nk,na,base,.TRUE.)

  ELSE IF(imode==2) THEN

     !Compare to cosmo-OWLS models

     !Set OWLS model, or set to zero for defaults
     iowl=4
     IF(iowl==1) THEN
        name='DMONLY'
     ELSE IF(iowl==2) THEN
        name='REF'
        cosm%param(1)=2.
        cosm%param(2)=1.4
        cosm%param(3)=1.24
        cosm%param(4)=1e13
        cosm%param(5)=0.055
     ELSE IF(iowl==3) THEN
        name='NOCOOL'
        cosm%param(1)=2.
        cosm%param(2)=0.8
        cosm%param(3)=1.1
        cosm%param(4)=0.
        cosm%param(5)=0.
     ELSE IF(iowl==4) THEN
        name='AGN'
        cosm%param(1)=2.
        cosm%param(2)=0.5
        cosm%param(3)=1.18
        cosm%param(4)=8e13
        cosm%param(5)=0.0225
     ELSE IF(iowl==5) THEN
        name='AGN 8.5'
        cosm%param(1)=2.
        cosm%param(2)=-0.5
        cosm%param(3)=1.26
        cosm%param(4)=2d14
        cosm%param(5)=0.0175
     ELSE IF(iowl==6) THEN
        name='AGN 8.7'
        cosm%param(1)=2.
        cosm%param(2)=-2.
        cosm%param(3)=1.3
        cosm%param(4)=1e15
        cosm%param(5)=0.015
     END IF

     IF(iowl .NE. 0) WRITE(*,*) 'Comparing to OWLS model: ', TRIM(name)
     CALL print_baryon_parameters(cosm)

     !Set number of k points and k range (log spaced)
     nk=200
     kmin=1e-3
     kmax=1e2
     CALL fill_array(log(kmin),log(kmax),k,nk)
     k=exp(k)
     ALLOCATE(pow_lin(nk),pow_2h(nk),pow_1h(nk),pow_full(nk))

     !Set the redshift
     !na=1
     z=0.
     !ALLOCATE(a(na))
     !a(1)=scale_factor_z(z)

     !Assigns the cosmological model
     icosmo=1
     CALL assign_cosmology(icosmo,cosm)

     !Normalises power spectrum (via sigma_8) and fills sigma(R) look-up tables
     CALL initialise_cosmology(cosm)
     CALL print_cosmology(cosm)

     !Initiliasation for the halomodel calcualtion
     CALL halomod_init(mmin,mmax,z,lut,cosm)

     !Runs the diagnostics
     dir='diagnostics'
     CALL halo_diagnostics(z,lut,cosm,dir)
     CALL halo_definitions(lut,dir)

     !File base and extension
     base='cosmo-OWLS/data/power_'
     mid=''
     ext='.dat'

     !Dark-matter only
     outfile='cosmo-OWLS/data/DMONLY.dat'
     WRITE(*,fmt='(2I5,A30)') -1, -1, TRIM(outfile)
     CALL calculate_halomod(-1,-1,k,nk,z,pow_lin,pow_2h,pow_1h,pow_full,lut,cosm)
     CALL write_power(k,pow_lin,pow_2h,pow_1h,pow_full,nk,outfile)

     !ip=-1
     !CALL calculate_halomod(ip,k,nk,a,na,powa_lin,powa_2h,powa_1h,powa_full,cosm)
     !outfile='cosmo-OWLS/data/DMONLY.dat'
     !CALL write_power(k,powa_lin(:,1),powa_2h(:,1),powa_1h(:,1),powa_full(:,1),nk,outfile)

     !Loop over matter types and do auto and cross-spectra
     DO j1=0,6
        DO j2=j1,6

           !Fix output file and write to screen
           outfile=number_file2(base,j1,mid,j2,ext)
           WRITE(*,fmt='(2I5,A30)') j1, j2, TRIM(outfile)

           CALL calculate_halomod(j1,j2,k,nk,z,pow_lin,pow_2h,pow_1h,pow_full,lut,cosm)
           CALL write_power(k,pow_lin,pow_2h,pow_1h,pow_full,nk,outfile)

        END DO
     END DO

  ELSE IF(imode==3) THEN

     !Assigns the cosmological model
     icosmo=0
     CALL assign_cosmology(icosmo,cosm)

     !Normalises power spectrum (via sigma_8) and fills sigma(R) look-up tables
     CALL initialise_cosmology(cosm)
     CALL print_cosmology(cosm)

     !WRITE(*,*) 'Redshift:'
     !READ(*,*) z
     !WRITE(*,*)
     z=0.

     !Initiliasation for the halomodel calcualtion
     CALL halomod_init(mmin,mmax,z,lut,cosm)

     !Runs the diagnostics
     dir='diagnostics'
     CALL halo_diagnostics(z,lut,cosm,dir)
     CALL halo_definitions(lut,dir)

     !output='winint/integrand.dat'
     !irho=14
     !rv=1.
     !rs=0.25
     !rmax=rv
     !CALL winint_diagnostics(rmax,rv,rs,irho,output)

  ELSE IF(imode==4) THEN

     STOP 'Error, random mode not implemented yet'

     !Ignore this, only useful for bug tests
     CALL RNG_set(0)

     !Only not uncommented to suppress compile-time warnings
     DO
        CALL random_cosmology(cosm)   
     END DO
     !Ignore this, only useful for bug tests

  ELSE IF(imode==5) THEN

     !Compare to cosmo-OWLS models for pressure

     !Set number of k points and k range (log spaced)
     nk=200
     kmin=1.e-3
     kmax=1.e2
     CALL fill_array(log(kmin),log(kmax),k,nk)
     k=exp(k)
     ALLOCATE(pow_lin(nk),pow_2h(nk),pow_1h(nk),pow_full(nk))

     !Set the redshift
     z=0.

     !Assigns the cosmological model
     icosmo=1
     CALL assign_cosmology(icosmo,cosm)

     !Normalises power spectrum (via sigma_8) and fills sigma(R) look-up tables
     CALL initialise_cosmology(cosm)
     CALL print_cosmology(cosm)

     !Initiliasation for the halomodel calcualtion
     CALL halomod_init(mmin,mmax,z,lut,cosm)

     !Runs the diagnostics
     dir='diagnostics'
     CALL halo_diagnostics(z,lut,cosm,dir)
     CALL halo_definitions(lut,dir)

     !File base and extension
     base='pressure/power_'
     ext='.dat'

     !Number of different spectra
     n=3

     !Do the calculation
     DO j=0,n

        IF(j==0) THEN
           !DMONLY
           j1=-1
           j2=-1
           outfile='pressure/DMONLY.dat'
        ELSE IF(j==1) THEN
           !Matter x matter
           j1=0
           j2=0
           outfile='dd'
        ELSE IF(j==2) THEN
           !Matter x pressure
           j1=0
           j2=6
           outfile='dp'
        ELSE IF(j==3) THEN
           !Pressure x pressure
           j1=6
           j2=6
           outfile='pp'
        END IF

        IF(j .NE. 0) outfile=TRIM(base)//TRIM(outfile)//TRIM(ext)

        WRITE(*,fmt='(3I5,A30)') j, j1, j2, TRIM(outfile)

        CALL calculate_halomod(j1,j2,k,nk,z,pow_lin,pow_2h,pow_1h,pow_full,lut,cosm)
        CALL write_power(k,pow_lin,pow_2h,pow_1h,pow_full,nk,outfile)

     END DO

  ELSE IF(imode==6) THEN

     !n(z) normalisation check

     WRITE(*,*) 'HMx: Checking n(z) functions'

     !inz(1)=-1
     !inz(2)=0

     nnz=7
     DO i=1,nnz
        IF(i==1) nz=1
        IF(i==2) nz=4
        IF(i==3) nz=5
        IF(i==4) nz=6
        IF(i==5) nz=7
        IF(i==6) nz=8
        IF(i==7) nz=9
        WRITE(*,*) 'HMx: n(z) number:', nz
        CALL get_nz(nz,lens)
        WRITE(*,*) 'HMx: n(z) integral (linear):', integrate_table(lens%z_nz,lens%nz,lens%nnz,1,lens%nnz,1)
        WRITE(*,*) 'HMx: n(z) integral (quadratic):', integrate_table(lens%z_nz,lens%nz,lens%nnz,1,lens%nnz,2)
        WRITE(*,*) 'HMx: n(z) integral (cubic):', integrate_table(lens%z_nz,lens%nz,lens%nnz,2,lens%nnz,3)
        WRITE(*,*)
     END DO

  ELSE IF(imode==7 .OR. imode==8 .OR. imode==9 .OR. imode==10 .OR. imode==11) THEN

     !General stuff for all cross correlations

     !Set the fields
     ix=-1
     CALL set_ix(ix,ip)

     !Assign the cosmological model
     icosmo=-1
     CALL assign_cosmology(icosmo,cosm)

     !Set the k range
     kmin=1e-3
     kmax=1e2
     nk=200

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
     lmax=1e5
     nl=128

     !Allocate arrays for l and C(l)
     CALL fill_array(log(lmin),log(lmax),ell,nl)
     ell=exp(ell)
     ALLOCATE(Cell(nl))

     !Set the angular arrays in degrees
     thmin=0.01
     thmax=10.
     nth=128

     !Allocate arrays for theta and xi(theta)
     CALL fill_array(log(thmin),log(thmax),theta,nth)
     theta=exp(theta)
     ALLOCATE(xi(3,nth))

     WRITE(*,*) 'HMx: Cross-correlation information'
     WRITE(*,*) 'HMx: output directiory: ', TRIM(dir)
     WRITE(*,*) 'HMx: Profile type 1: ', TRIM(halo_type(ip(1)))
     WRITE(*,*) 'HMx: Profile type 2: ', TRIM(halo_type(ip(2)))
     !WRITE(*,*) 'HMx: Kernel type 1: ', TRIM(kernel_type(ik(1)))
     !WRITE(*,*) 'HMx: Kernel type 2: ', TRIM(kernel_type(ik(2)))
     WRITE(*,*) 'HMx: Xcorr type 1: ', TRIM(xcorr_type(ix(1)))
     WRITE(*,*) 'HMx: Xcorr type 2: ', TRIM(xcorr_type(ix(2)))
     WRITE(*,*) 'HMx: P(k) minimum k [h/Mpc]:', REAL(kmin)
     WRITE(*,*) 'HMx: P(k) maximum k [h/Mpc]:', REAL(kmax)
     WRITE(*,*) 'HMx: minimum a:', REAL(amin)
     WRITE(*,*) 'HMx: maximum a:', REAL(amax)
     WRITE(*,*) 'HMx: minimum ell:', REAL(lmin)
     WRITE(*,*) 'HMx: maximum ell:', REAL(lmax)
     WRITE(*,*)     

     IF(imode==7) THEN

        !Normalises power spectrum (via sigma_8) and fills sigma(R) look-up tables
        CALL initialise_cosmology(cosm)
        CALL print_cosmology(cosm)

        !Initialise the lensing part of the calculation
        CALL initialise_distances(cosm)
        CALL write_distances(cosm)

        !Write out diagnostics
        CALL halomod_init(mmin,mmax,z,lut,cosm)
        dir='diagnostics'
        CALL halo_diagnostics(z,lut,cosm,dir)
        CALL halo_definitions(lut,dir)

        CALL calculate_HMx(ip,k,nk,a,na,powa_lin,powa_2h,powa_1h,powa_full,cosm)

        !Loop over scale factors
        !DO j=na,1,-1
        !
        !   z=-1+1./a(j)
        !   
        !   !Initiliasation for the halomodel calcualtion
        !   CALL halomod_init(mmin,mmax,z,lut,cosm)
        !   CALL calculate_halomod(ip(1),ip(2),k,nk,z,powa_lin(:,j),powa_2h(:,j),powa_1h(:,j),powa_full(:,j),lut,cosm)
        !
        !   IF(z==0.) THEN
        !      dir='diagnostics'
        !      CALL halo_diagnostics(z,lut,cosm,dir)
        !      CALL halo_definitions(lut,dir)
        !   END IF

        !  !Write progress to screen
        !  IF(j==na) THEN
        !     WRITE(*,fmt='(A5,A7)') 'i', 'a'
        !     WRITE(*,fmt='(A13)') '   ============'
        !  END IF
        !  WRITE(*,fmt='(I5,F8.3)') j, a(j)

        !END DO
        !WRITE(*,fmt='(A13)') '   ============'
        !WRITE(*,*)

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
        CALL write_power_a_multiple(k,a,powa_lin,powa_2h,powa_1h,powa_full,nk,na,base,.TRUE.)

        !Fill out the projection kernels
        CALL fill_projection_kernels(ix,proj,cosm)
        CALL write_projection_kernels(proj,cosm)
        !output='projection/kernel1.dat'
        !CALL write_projection_kernel(proj(1),cosm,output)
        !output='projection/kernel2.dat'
        !CALL write_projection_kernel(proj(2),cosm,output)

        !Set the distance range for the Limber integral
        r1=0. !100.
        r2=maxdist(proj)!proj%rs

        !Write to screen
        WRITE(*,*) 'HMx: Computing C(l)'
        WRITE(*,*) 'HMx: ell min:', REAL(ell(1))
        WRITE(*,*) 'HMx: ell max:', REAL(ell(nl))
        WRITE(*,*) 'HMx: number of ell:', nl
        WRITE(*,*) 'HMx: lower limit of Limber integral [Mpc/h]:', REAL(r1)
        WRITE(*,*) 'HMx: upper limit of Limber integral [Mpc/h]:', REAL(r2)
        WRITE(*,*)

        !Loop over all types of C(l) to create
        DO j=1,4

           IF(ifull .AND. (j .NE. 4)) CYCLE
           !IF(j==3) CYCLE !Skip the fucking one-halo term

           !Write information to screen
           IF(j==1) THEN
              WRITE(*,*) 'HMx: Doing linear'
              outfile=TRIM(dir)//'cl_linear.dat'
              powa=powa_lin
           ELSE IF(j==2) THEN
              WRITE(*,*) 'HMx: Doing 2-halo'
              outfile=TRIM(dir)//'cl_2halo.dat'
              powa=powa_2h
           ELSE IF(j==3) THEN
              WRITE(*,*) 'HMx: Doing 1-halo'
              outfile=TRIM(dir)//'cl_1halo.dat'
              powa=powa_1h
           ELSE IF(j==4) THEN
              WRITE(*,*) 'HMx: Doing full'
              outfile=TRIM(dir)//'cl_full.dat'
              powa=powa_full
           END IF

           WRITE(*,*) 'HMx: Output: ', TRIM(outfile)

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
              WRITE(*,*) 'HMx: Output: ', TRIM(outfile)

              !Actually calculate the xi(theta)
              CALL calculate_xi(theta,xi,nth,ell,Cell,nl,NINT(lmax))
              CALL write_xi(theta,xi,nth,outfile)

           END IF

        END DO
        WRITE(*,*) 'HMx: Done'
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

           !Normalises power spectrum (via sigma_8) and fills sigma(R) look-up tables
           CALL initialise_cosmology(cosm)
           CALL print_cosmology(cosm)
           CALL initialise_distances(cosm)
           !CALL write_distances(cosm)

           !Loop over redshifts
           !DO j=1,na
           !
           !   z=redshift_a(a(j))
           !   
           !   !Initiliasation for the halomodel calcualtion
           !   CALL halomod_init(mmin,mmax,z,lut,cosm)
           !   CALL calculate_halomod(ip(1),ip(2),k,nk,z,powa_lin(:,j),powa_2h(:,j),powa_1h(:,j),powa_full(:,j),lut,cosm)
           !   !Write progress to screen
           !   IF(j==1) THEN
           !      WRITE(*,fmt='(A5,A7)') 'i', 'a'
           !      WRITE(*,fmt='(A13)') '   ============'
           !   END IF
           !   WRITE(*,fmt='(I5,F8.3)') j, a(j)
           !
           !END DO
           !WRITE(*,fmt='(A13)') '   ============'
           !WRITE(*,*)

           CALL calculate_HMx(ip,k,nk,a,na,powa_lin,powa_2h,powa_1h,powa_full,cosm)

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
           WRITE(*,*) 'HMx: Computing C(l)'
           WRITE(*,*) 'HMx: ell min:', REAL(ell(1))
           WRITE(*,*) 'HMx: ell max:', REAL(ell(nl))
           WRITE(*,*) 'HMx: number of ell:', nl
           WRITE(*,*)

           !Loop over all types of C(l) to create  
           base=TRIM(dir)//'cosmology_'
           DO j=1,4
              IF(j==1) THEN
                 WRITE(*,*) 'HMx: Doing C(l) linear'
                 ext='_cl_linear.dat'
                 powa=powa_lin
              ELSE IF(j==2) THEN
                 WRITE(*,*) 'HMx: Doing C(l) 2-halo'
                 ext='_cl_2halo.dat'
                 powa=powa_2h
              ELSE IF(j==3) THEN
                 WRITE(*,*) 'HMx: Doing C(l) 1-halo'
                 ext='_cl_1halo.dat'
                 powa=powa_1h
              ELSE IF(j==4) THEN
                 WRITE(*,*) 'HMx: Doing C(l) full'
                 ext='_cl_full.dat'
                 powa=powa_full
              END IF
              outfile=number_file(base,i,ext)
              !Actually calculate the C(ell)
              CALL calculate_Cell(0.,maxdist(proj),ell,Cell,nl,k,a,powa,nk,na,proj,cosm)
              CALL write_Cell(ell,Cell,nl,outfile)
           END DO
           WRITE(*,*) 'HMx: Done'
           WRITE(*,*)

        END DO

     ELSE IF(imode==9) THEN

        !Break down cross-correlation in terms of mass

        !Normalises power spectrum (via sigma_8) and fills sigma(R) look-up tables
        CALL initialise_cosmology(cosm)
        CALL print_cosmology(cosm)

        !Initialise the lensing part of the calculation
        CALL initialise_distances(cosm)
        CALL write_distances(cosm)

        !Fill out the projection kernels
        CALL fill_projection_kernels(ix,proj,cosm)
        CALL write_projection_kernels(proj,cosm)

        DO i=0,6
           IF(icumulative==0) THEN
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
           ELSE IF(icumulative==1) THEN
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
           ELSE
              STOP 'HMx: Error, icumulative not set correctly.'
           END IF

           !Set the code to not 'correct' the two-halo power for missing
           !mass when doing the calcultion binned in halo mass
           IF(icumulative==0 .AND. i>1) ip2h=0
           !IF(icumulative==1 .AND. i>0) ip2h=0

           WRITE(*,fmt='(A16)') 'HMx: Mass range'
           WRITE(*,fmt='(A16,I5)') 'HMx: Iteration:', i
           WRITE(*,fmt='(A21,2ES15.7)') 'HMx: M_min [Msun/h]:', m1
           WRITE(*,fmt='(A21,2ES15.7)') 'HMx: M_max [Msun/h]:', m2
           WRITE(*,*)

           !Loop over redshifts
           DO j=1,na

              z=redshift_a(a(j))

              !Initiliasation for the halomodel calcualtion
              CALL halomod_init(m1,m2,z,lut,cosm)
              CALL calculate_halomod(ip(1),ip(2),k,nk,z,powa_lin(:,j),powa_2h(:,j),powa_1h(:,j),powa_full(:,j),lut,cosm)

              !Write progress to screen
              IF(j==1) THEN
                 WRITE(*,fmt='(A5,A7)') 'i', 'a'
                 WRITE(*,fmt='(A13)') '   ============'
              END IF
              WRITE(*,fmt='(I5,F8.3)') j, a(j)

           END DO
           WRITE(*,fmt='(A13)') '   ============'
           WRITE(*,*)

           IF(i==0) THEN
              outfile=TRIM(dir)//'power'
           ELSE
              base=TRIM(dir)//'mass_'
              mid='_'
              ext='_power'
              outfile=number_file2(base,NINT(log10(m1)),mid,NINT(log10(m2)),ext)
           END IF
           WRITE(*,*) 'HMx: File: ', TRIM(outfile)
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

              WRITE(*,*) 'HMx: File: ', TRIM(outfile)

              CALL calculate_Cell(0.,maxdist(proj),ell,Cell,nl,k,a,powa,nk,na,proj,cosm)
              CALL write_Cell(ell,Cell,nl,outfile)

           END DO
           WRITE(*,*) 'HMx: Done'
           WRITE(*,*)

        END DO

     ELSE IF(imode==10) THEN

        !Break down cross-correlation in terms of redshift

        !Normalises power spectrum (via sigma_8) and fills sigma(R) look-up tables
        CALL initialise_cosmology(cosm)
        CALL print_cosmology(cosm)

        !Loop over redshifts
        !DO j=1,na
        !
        !   !z=-1.+1./a(j)
        !   z=redshift_a(a(j))
        !
        !   !Initiliasation for the halomodel calcualtion
        !   CALL halomod_init(mmin,mmax,z,lut,cosm)
        !   CALL calculate_halomod(ip(1),ip(2),k,nk,z,powa_lin(:,j),powa_2h(:,j),powa_1h(:,j),powa_full(:,j),lut,cosm)
        !
        !   !Write progress to screen
        !   IF(j==1) THEN
        !      WRITE(*,fmt='(A5,A7)') 'i', 'a'
        !      WRITE(*,fmt='(A13)') '   ============'
        !   END IF
        !   WRITE(*,fmt='(I5,F8.3)') j, a(j)
        !
        !END DO
        !WRITE(*,fmt='(A13)') '   ============'
        !WRITE(*,*)

        CALL calculate_HMx(ip,k,nk,a,na,powa_lin,powa_2h,powa_1h,powa_full,cosm)

        !output=TRIM(base)//'power'
        !CALL write_power_a(k,a,powa,nk,na,output)

        !Initialise the lensing part of the calculation
        CALL initialise_distances(cosm)
        CALL write_distances(cosm)

        !Fill out the projection kernels
        CALL fill_projection_kernels(ix,proj,cosm)
        CALL write_projection_kernels(proj,cosm)

        !Write to screen
        WRITE(*,*) 'HMx: Computing C(l)'
        WRITE(*,*) 'HMx: ell min:', REAL(ell(1))
        WRITE(*,*) 'HMx: ell max:', REAL(ell(nl))
        WRITE(*,*) 'HMx: number of ell:', nl
        WRITE(*,*)

        zmin=0.
        zmax=1.
        nz=8

        DO i=0,nz

           IF(i==0) THEN
              !z1=0.
              !z2=3.99 !Just less than z=4 to avoid rounding error
              r1=0.
              r2=maxdist(proj)
           ELSE
              IF(icumulative==0) THEN
                 !z1=zmin+(zmax-zmin)*float(i-1)/float(nz)
                 z1=progression(zmin,zmax,i,nz)
              ELSE IF(icumulative==1) THEN
                 z1=zmin
              END IF
              z2=zmin+(zmax-zmin)*float(i)/float(nz)
              r1=cosmic_distance(z1,cosm)
              r2=cosmic_distance(z2,cosm)
           END IF

           WRITE(*,*) 'HMx:', i
           IF(i>0) THEN
              WRITE(*,*) 'HMx: z1', REAL(z1)
              WRITE(*,*) 'HMx: z2', REAL(z2)
           END IF
           WRITE(*,*) 'HMx: r1 [Mpc/h]', REAL(r1)
           WRITE(*,*) 'HMx: r2 [Mpc/h]', REAL(r2)

           !Loop over all types of C(l) to create
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

              IF(i>0 .AND. icumulative==0) THEN
                 outfile=number_file2(base,i-1,mid,i,ext)
              ELSE IF(i>0 .AND. icumulative==1) THEN
                 outfile=number_file2(base,0,mid,i,ext)
              END IF
              WRITE(*,*) 'HMx: Output: ', TRIM(outfile)

              !This crashes for the low r2 values for some reason
              CALL calculate_Cell(r1,r2,ell,Cell,nl,k,a,powa,nk,na,proj,cosm)
              CALL write_Cell(ell,Cell,nl,outfile)

           END DO
           WRITE(*,*)

        END DO

        WRITE(*,*) 'HMx: Done'
        WRITE(*,*)

     ELSE IF(imode==11) THEN

        STOP 'HMx: Error, breakdown in radius is not supported yet'

     ELSE

        STOP 'HMx: Error, you have specified the mode incorrectly'

     END IF

  ELSE IF(imode==12) THEN

     !Project triad

     dir='data'

     !Assigns the cosmological model
     icosmo=0
     CALL assign_cosmology(icosmo,cosm)

     !Normalises power spectrum (via sigma_8) and fills sigma(R) look-up tables
     CALL initialise_cosmology(cosm)
     CALL print_cosmology(cosm)

     !Initialise the lensing part of the calculation
     CALL initialise_distances(cosm)
     CALL write_distances(cosm)

     !Set the ell range
     lmin=100.
     lmax=3000.
     nl=64

     !Allocate arrays for l and C(l)
     CALL fill_array(log(lmin),log(lmax),ell,nl)
     ell=exp(ell)
     ALLOCATE(Cell(nl))

     WRITE(*,*) 'HMx: Cross-correlation information'
     WRITE(*,*) 'HMx: output directiory: ', TRIM(dir)
     WRITE(*,*) 'HMx: minimum ell:', REAL(lmin)
     WRITE(*,*) 'HMx: maximum ell:', REAL(lmax)
     WRITE(*,*) 'HMx: number of ell:', nl
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

        CALL xcorr(ix,ell,Cell,nl,cosm,.TRUE.)
        CALL write_Cell(ell,Cell,nl,outfile)

        WRITE(*,*) 'HMx: Done'
        WRITE(*,*)

     END DO

  ELSE IF(imode==13) THEN

     !Calculate the cross-correlation coefficient

     !Assign the cosmology
     icosmo=-1
     CALL assign_cosmology(icosmo,cosm)

     !Normalises power spectrum (via sigma_8) and fills sigma(R) look-up tables
     CALL initialise_cosmology(cosm)
     CALL print_cosmology(cosm)

     !Initialise the lensing part of the calculation
     CALL initialise_distances(cosm)
     CALL write_distances(cosm)

     !Set the ell range and allocate arrays for l and C(l)
     lmin=1e0
     lmax=1e5
     nl=64 
     CALL fill_array(log(lmin),log(lmax),ell,nl)
     ell=exp(ell)
     ALLOCATE(Cell(nl))

     dir='data'

     ixx=-1
     CALL set_ix(ixx,ip)

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
        CALL xcorr(ix,ell,Cell,nl,cosm,.TRUE.)
        CALL write_Cell(ell,Cell,nl,outfile)
     END DO

  ELSE IF(imode==14) THEN

     !Make power spectra as a function of parameter variations

     !Number of values to try for each parameter
     m=9

     !Set number of k points and k range (log spaced)
     nk=128
     kmin=1.e-3
     kmax=1.e1
     CALL fill_array(log(kmin),log(kmax),k,nk)
     k=exp(k)
     ALLOCATE(pow_lin(nk),pow_2h(nk),pow_1h(nk),pow_full(nk))

     !Set the redshift
     z=0.

     !Assigns the cosmological model
     icosmo=1
     CALL assign_cosmology(icosmo,cosm)

     !Normalises power spectrum (via sigma_8) and fills sigma(R) look-up tables
     CALL initialise_cosmology(cosm)
     CALL print_cosmology(cosm)

     !Initiliasation for the halomodel calcualtion
     CALL halomod_init(mmin,mmax,z,lut,cosm)  

     !DMONLY
     j1=-1
     j2=-1
     outfile='variations/DMONLY.dat'
     CALL calculate_halomod(j1,j2,k,nk,z,pow_lin,pow_2h,pow_1h,pow_full,lut,cosm)
     CALL write_power(k,pow_lin,pow_2h,pow_1h,pow_full,nk,outfile)

     !Loop over parameters
     DO ipa=1,cosm%np
        !DO ipa=2,2

        cosm%param=cosm%param_defaults

        !Loop over parameter values
        DO i=1,m

           !Set the parameter value that is being varied
           IF(cosm%param_log(ipa)) THEN
              cosm%param(ipa)=progression(log(cosm%param_min(ipa)),log(cosm%param_max(ipa)),i,m)
              cosm%param(ipa)=exp(cosm%param(ipa))
           ELSE
              cosm%param(ipa)=progression(cosm%param_min(ipa),cosm%param_max(ipa),i,m)
           END IF

           !Write out halo matter and pressure profile information
           !All the string crap is in the loop for a reason
           DO j=10,16
              base='variations/profile_mass_'
              ext='_param_'
              base=number_file(base,j,ext)
              mid='_value_'
              ext='.dat'
              outfile=number_file2(base,ipa,mid,i,ext)
              mass=10.**j
              CALL write_halo_profiles(mass,z,lut,cosm,outfile)
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
                 !matter-matter
                 j1=0
                 j2=0
                 ext='_dd.dat'
              ELSE IF(j==2) THEN
                 !matter-pressure
                 j1=0
                 j2=6
                 ext='_dp.dat'
              ELSE IF(j==3) THEN
                 !pressure-pressure
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
              WRITE(*,fmt='(4I5,A50)') ipa, i, j1, j2, TRIM(outfile)

              !Do the halo-model calculation and write to file
              CALL calculate_halomod(j1,j2,k,nk,z,pow_lin,pow_2h,pow_1h,pow_full,lut,cosm)
              CALL write_power(k,pow_lin,pow_2h,pow_1h,pow_full,nk,outfile)

           END DO

        END DO

     END DO

  ELSE

     STOP 'HMx: Error, you have specified the mode incorrectly'

  END IF

END PROGRAM HMx_driver
