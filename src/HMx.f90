MODULE HMx

  !Module usage statements
  !USE cosdef
  USE constants
  USE array_operations
  USE file_info
  USE solve_equations
  USE special_functions
  USE interpolate
  USE string_operations
  USE calculus
  USE calculus_table
  USE cosmology_functions
  
  IMPLICIT NONE
  INTEGER, PARAMETER :: imf=2 !Set mass function (1 - PS, 2 - ST) !Move to 'tables' type eventually 
  INTEGER, PARAMETER :: imead=0 !Set to do Mead et al. (2015,2016) accurate calculation !Move to 'tables' type eventually 
  REAL, PARAMETER :: acc=1e-4 !Global integration-accuracy parameter
  REAL, PARAMETER :: null=0.d0

  !Halo-model stuff that needs to be recalculated for each new z
  TYPE tables     
     REAL, ALLOCATABLE :: c(:), rv(:), nu(:), sig(:), zc(:), m(:), rr(:), sigf(:)
     REAL, ALLOCATABLE :: r500(:), m500(:), c500(:), r200(:), m200(:), c200(:)
     REAL, ALLOCATABLE :: r500c(:), m500c(:), c500c(:), r200c(:), m200c(:), c200c(:)
     REAL, ALLOCATABLE :: log_m(:)
     REAL :: sigv, sigv100, c3, knl, rnl, neff, sig8z
     REAL :: gmin, gmax, gbmin, gbmax
     INTEGER :: ip2h, ibias
     INTEGER :: n
     LOGICAL :: void !Do voids or not
  END TYPE tables

CONTAINS

!!$  SUBROUTINE init_HMx(cosm)
!!$
!!$    !Sets values for the baryon parameters. Probably a poor choice of subroutine name.
!!$    IMPLICIT NONE
!!$    TYPE(cosmology), INTENT(INOUT) :: cosm
!!$
!!$    !Names of variable parameters
!!$    !cosm%param_names(1)='alpha' !alpha in virial temperature (turbulence?)
!!$    !cosm%param_names(2)='Dc' !Change in NFW concentration due to gas
!!$    !cosm%param_names(3)='Gamma' !Gamma in Komatsu-Seljak profile
!!$    !cosm%param_names(4)='M_B' !Halo mass at which free and unbound gas are equal
!!$    !cosm%param_names(5)='A_*' !Prefactor for stellar fraction
!!$
!!$    !Set default values for variable parameters
!!$    !cosm%param(1)=1.
!!$    !cosm%param(2)=0.
!!$    !cosm%param(3)=1.18
!!$    !cosm%param(4)=1.2d14
!!$    !cosm%param(5)=0.02
!!$
!!$    !Set default values for variable parameters
!!$    cosm%alpha=1.
!!$    cosm%Dc=0.
!!$    cosm%Gamma=1.18
!!$    cosm%M0=1.2e14
!!$    cosm%Astar=0.02
!!$
!!$    !Set some default parameters
!!$    !cosm%param_defaults=cosm%param
!!$
!!$    !Minimum parameter values in variation
!!$    !cosm%param_min(1)=0.4
!!$    !cosm%param_min(2)=0.
!!$    !cosm%param_min(3)=1.10
!!$    !cosm%param_min(4)=1e13
!!$    !cosm%param_min(5)=0.015
!!$
!!$    !Maximum parameter values in variation
!!$    !cosm%param_max(1)=2.
!!$    !cosm%param_max(2)=2.
!!$    !cosm%param_max(3)=1.26
!!$    !cosm%param_max(4)=1e15
!!$    !cosm%param_max(5)=0.055
!!$
!!$    !Should the range be explored in log?
!!$    !cosm%param_log(1)=.FALSE.
!!$    !cosm%param_log(2)=.FALSE.
!!$    !cosm%param_log(3)=.FALSE.
!!$    !cosm%param_log(4)=.TRUE.
!!$    !cosm%param_log(5)=.FALSE.
!!$
!!$  END SUBROUTINE init_HMx

  FUNCTION xcorr_type(ix)

    !Names for cross-correlation field types
    IMPLICIT NONE
    CHARACTER(len=256) :: xcorr_type
    INTEGER, INTENT(IN) :: ix

    xcorr_type=''
    IF(ix==1)  xcorr_type='RCSLenS lensing'
    IF(ix==2)  xcorr_type='Compton y'
    IF(ix==3)  xcorr_type='CMB lensing'
    IF(ix==4)  xcorr_type='CFHTLenS lensing'
    IF(ix==5)  xcorr_type='KiDS lensing (z = 0.1 -> 0.9)'
    IF(ix==6)  xcorr_type='KiDS lensing (z = 0.1 -> 0.3)'
    IF(ix==7)  xcorr_type='KiDS lensing (z = 0.3 -> 0.5)'
    IF(ix==8)  xcorr_type='KiDS lensing (z = 0.5 -> 0.7)'
    IF(ix==9)  xcorr_type='KiDS lensing (z = 0.7 -> 0.9)'
    IF(ix==10) xcorr_type='Gravitational waves'
    IF(xcorr_type=='') STOP 'XCORR_TYPE: Error, ix not specified correctly'
    
  END FUNCTION xcorr_type

  FUNCTION halo_type(i)

    !Name halo types
    IMPLICIT NONE
    CHARACTER(len=256) :: halo_type
    INTEGER :: i
    
    halo_type=''
    IF(i==-1) halo_type='DMONLY'
    IF(i==0)  halo_type='Matter'
    IF(i==1)  halo_type='CDM'
    IF(i==2)  halo_type='Gas'
    IF(i==3)  halo_type='Star'
    IF(i==4)  halo_type='Bound gas'
    IF(i==5)  halo_type='Free gas'
    IF(i==6)  halo_type='Pressure'
    IF(i==7)  halo_type='Void'
    IF(i==8)  halo_type='Compensated void'
    IF(halo_type=='') STOP 'HALO_TYPE: Error, i not specified correctly'
    
  END FUNCTION halo_type

  SUBROUTINE calculate_HMx(itype,mmin,mmax,k,nk,a,na,powa_lin,powa_2h,powa_1h,powa_full,cosm,verbose)

    IMPLICIT NONE
    REAL, INTENT(IN) :: k(:), a(:)
    INTEGER, INTENT(IN) :: nk, na, itype(2)
    !REAL, ALLOCATABLE, INTENT(INOUT) :: powa_lin(:,:) !Mead - commented out
    REAL, ALLOCATABLE, INTENT(OUT) :: powa_2h(:,:), powa_1h(:,:), powa_full(:,:), powa_lin(:,:) !Mead - added powa_lin here instead
    TYPE(cosmology), INTENT(IN) :: cosm
    LOGICAL, INTENT(IN) :: verbose
    REAL, INTENT(IN) :: mmin, mmax
    LOGICAL :: compute_p_lin
    INTEGER :: i
    REAL :: z
    TYPE(tables) :: lut
    LOGICAL :: verbose2

    !Mead - fixed this to always be false to avoid splurge of stuff printed to screen
    verbose2=.FALSE.

    !Tilman - added this
    IF(ALLOCATED(powa_lin)) THEN
       WRITE(*,*) 'Linear power spectrum provided.'
       compute_p_lin = .FALSE.
    ELSE
       !ALLOCATE(powa_lin(nk,na)) !Mead - commented out
       compute_p_lin = .TRUE.
    END IF
    !Tilman - End

    IF(ALLOCATED(powa_lin))  DEALLOCATE(powa_lin)
    IF(ALLOCATED(powa_2h))   DEALLOCATE(powa_2h)
    IF(ALLOCATED(powa_1h))   DEALLOCATE(powa_1h)
    IF(ALLOCATED(powa_full)) DEALLOCATE(powa_full)

    !Allocate power arrays
    !Mead - re-added powa_lin to be allocated here
    ALLOCATE(powa_lin(nk,na),powa_2h(nk,na),powa_1h(nk,na),powa_full(nk,na))

    !Do the halo-model calculation
    DO i=na,1,-1
       z=redshift_a(a(i))
       CALL halomod_init(mmin,mmax,z,lut,cosm,verbose2)
       !IF(verbose) WRITE(*,fmt='(A5,I5,F10.2)') 'HMx:', i, REAL(z)
       IF(i==na .and. verbose) WRITE(*,*) 'CALCULATE_HMx: Doing calculation'
       IF(verbose) WRITE(*,fmt='(A5,I5,F10.2)') 'HMx:', i, REAL(z) !Mead - re-added verbose dependence
       CALL calculate_halomod(itype(1),itype(2),k,nk,z,powa_lin(:,i),powa_2h(:,i),powa_1h(:,i),powa_full(:,i),lut,cosm,verbose2,compute_p_lin)
       verbose2=.FALSE.
    END DO
    IF(verbose) THEN
       WRITE(*,*) 'CALCULATE_HMx: Done'
       WRITE(*,*)
    END IF

  END SUBROUTINE calculate_HMx

  SUBROUTINE calculate_halomod(itype1,itype2,k,nk,z,pow_lin,pow_2h,pow_1h,pow,lut,cosm,verbose,compute_p_lin_arg)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: itype1, itype2
    INTEGER, INTENT(IN) :: nk
    REAL, INTENT(IN) :: k(nk), z
    REAL, INTENT(INOUT) :: pow_lin(nk)
    REAL, INTENT(OUT) ::pow_2h(nk), pow_1h(nk), pow(nk)
    TYPE(cosmology), INTENT(IN) :: cosm
    TYPE(tables), INTENT(IN) :: lut
    LOGICAL, OPTIONAL, INTENT(IN) :: compute_p_lin_arg
    LOGICAL, INTENT(IN) :: verbose
    INTEGER :: i
    REAL :: plin
    LOGICAL :: compute_p_lin

    !Tilman - added this
    IF(PRESENT(compute_p_lin_arg)) THEN
       compute_p_lin = compute_p_lin_arg
    ELSE
       compute_p_lin = .TRUE.
    END IF
    !compute_p_lin = .TRUE. !Mead - added this to make things work!!!
    !Tilman - done

    !Write to screen
    IF(verbose) THEN
       WRITE(*,*) 'CALCULATE_HALOMOD: k min:', k(1)
       WRITE(*,*) 'CALCULATE_HALOMOD: k max:', k(nk)
       WRITE(*,*) 'CALCULATE_HALOMOD: number of k:', nk
       WRITE(*,*) 'CALCULATE_HALOMOD: z:', z
       WRITE(*,*) 'CALCULATE_HALOMOD: Calculating halo-model power spectrum'
    END IF

    !Loop over k values
    !ADD OMP support properly. What is private and what is shared? CHECK THIS!
!!$OMP PARALLEL DO DEFAULT(SHARED), private(k,plin, pfull,p1h,p2h)
    DO i=1,nk

       !Tilman - added this
       IF(compute_p_lin) THEN
          !Get the linear power
          plin=p_lin(k(i),z,cosm)
          pow_lin(i)=plin
       END IF
       !Tilman - done

       !Do the halo model calculation
       CALL halomod(itype1,itype2,k(i),z,pow_2h(i),pow_1h(i),pow(i),plin,lut,cosm)

    END DO
!!$OMP END PARALLEL DO

    IF(verbose) THEN
       WRITE(*,*) 'CALCULATE_HALOMOD: Done'
       WRITE(*,*)
    END IF

  END SUBROUTINE calculate_halomod

  SUBROUTINE write_power(k,pow_lin,pow_2h,pow_1h,pow,nk,output,verbose)

    IMPLICIT NONE
    CHARACTER(len=256), INTENT(IN) :: output
    INTEGER, INTENT(IN) :: nk
    REAL, INTENT(IN) :: k(nk), pow_lin(nk), pow_2h(nk), pow_1h(nk), pow(nk)
    LOGICAL, INTENT(IN) :: verbose
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
       WRITE(*,*) 'WRITE_POWER_A: The first entry of the file is hashes - #####'
       WRITE(*,*) 'WRITE_POWER_A: The remainder of the first row are the scale factors - a'
       WRITE(*,*) 'WRITE_POWER_A: The remainder of the first column are the wave numbers - k'
       WRITE(*,*) 'WRITE_POWER_A: Each row then gives the power at that k and a'
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
    INTEGER :: m

    INTEGER, PARAMETER :: m1=10
    INTEGER, PARAMETER :: m2=16

    WRITE(*,*) 'HALO_DIAGNOSTICS: Outputting diagnostics'

    outfile=TRIM(dir)//'/mass_fractions.dat'
    WRITE(*,*) 'HALO_DIAGNOSTICS: ', TRIM(outfile)
    CALL write_mass_fractions(cosm,outfile)

    DO m=m1,m2

       mass=10.**m

       base=TRIM(dir)//'/halo_profile_m'
       ext='.dat'
       outfile=number_file(base,m,ext)
       WRITE(*,*) 'HALO_DIAGNOSTICS: ', TRIM(outfile)
       CALL write_halo_profiles(mass,z,lut,cosm,outfile)

       !base=TRIM(dir)//'/halo_profile_m'
       !ext='noh.dat'
       !outfile=number_file(base,m,ext)
       !CALL write_halo_profiles(cosm%h*mass,z,lut,cosm,outfile)

       base=TRIM(dir)//'/halo_window_m'
       ext='.dat'
       outfile=number_file(base,m,ext)
       WRITE(*,*) 'HALO_DIAGNOSTICS: ', TRIM(outfile)
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
    INTEGER :: i

    WRITE(*,*) 'HALO_DEFINITIONS: Outputting definitions'

    fradius=TRIM(dir)//'/radius.dat'
    fmass=TRIM(dir)//'/mass.dat'
    fconc=TRIM(dir)//'/concentration.dat'

    WRITE(*,*) 'HALO_DEFINITIONS: ', TRIM(fradius)
    WRITE(*,*) 'HALO_DEFINITIONS: ', TRIM(fmass)
    WRITE(*,*) 'HALO_DEFINITIONS: ', TRIM(fconc)

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

  FUNCTION fdamp(z,lut,cosm)

    IMPLICIT NONE
    REAL ::fdamp
    REAL, INTENT(IN) :: z
    TYPE(cosmology), INTENT(IN) :: cosm
    TYPE(tables), INTENT(IN) :: lut
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

    WRITE(*,*) 'PRINT_HALOMODEL_PARAMETERS: Writing out halo-model parameters'
    WRITE(*,*) 'PRINT_HALOMODEL_PARAMETERS: Halo-model parameters at your redshift'
    WRITE(*,*) '==========================='
    WRITE(*,fmt='(A10,F10.5)') 'z:', z
    WRITE(*,fmt='(A10,F10.5)') 'Dv:', Delta_v(z,cosm)
    WRITE(*,fmt='(A10,F10.5)') 'dc:', delta_c(z,cosm)
    WRITE(*,fmt='(A10,F10.5)') 'eta:', eta(z,cosm)
    WRITE(*,fmt='(A10,F10.5)') 'k*:', kstar(lut,cosm)
    WRITE(*,fmt='(A10,F10.5)') 'A:', As(cosm)
    WRITE(*,fmt='(A10,F10.5)') 'fdamp:', fdamp(z,lut,cosm)
    WRITE(*,fmt='(A10,F10.5)') 'alpha:', alpha_transition(lut,cosm)
    WRITE(*,*) '==========================='
    WRITE(*,*) 'PRINT_HALOMODEL_PARAMETERS: Done'
    WRITE(*,*)

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

  SUBROUTINE halomod_init(mmin,mmax,z,lut,cosm,verbose)

    !Halo-model initialisation routine
    !The computes other tables necessary for the one-halo integral
    IMPLICIT NONE
    REAL, INTENT(IN) :: z
    REAL, INTENT(IN) :: mmin, mmax
    LOGICAL, INTENT(IN) :: verbose
    TYPE(tables), INTENT(OUT) :: lut !Or is this just OUT?
    TYPE(cosmology), INTENT(IN) :: cosm
    INTEGER :: i
    REAL :: Dv, dc, f, m, nu, r, sig, A0, rhom, rhoc
    
    INTEGER, PARAMETER :: n=64 !Number of mass entries in look-up table
    REAL, PARAMETER :: large_nu=10. !Value for nu such that there are no haloes larger

    !Set method to correct for missing integrand in two-halo term
    !0 - Do nothing
    !1 - Add value of integral assuming that W(k)=1
    !2 - Put the missing part of the integrand as a delta function at nu1
    lut%ip2h=2

    !Order to go to in bias
    !1 - First order
    !2 - Second order
    lut%ibias=1

    !Do voids or not
    lut%void=.FALSE.

    !Find value of sigma_v
    lut%sigv=sqrt(dispint(0.,z,acc,cosm)/3.)
    lut%sigv100=sqrt(dispint(100.,z,acc,cosm)/3.)
    lut%sig8z=sigma(8.,z,cosm)

    IF(verbose) THEN
       WRITE(*,*) 'HALOMOD_INIT: Filling look-up tables'
       WRITE(*,*) 'HALOMOD_INIT: tables being filled at redshift:', REAL(z)
       WRITE(*,*) 'HALOMOD_INIT: sigv [Mpc/h]:', REAL(lut%sigv)
       WRITE(*,*) 'HALOMOD_INIT: sigv100 [Mpc/h]:', REAL(lut%sigv100)
       WRITE(*,*) 'HALOMOD_INIT: sig8(z):', REAL(lut%sig8z)
    END IF

    !Remove this if LUT is INTENT(OUT)
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
       WRITE(*,*) 'HALOMOD_INIT: minimum log10(M/[Msun/h]):', REAL(log10(lut%m(1)))
       WRITE(*,*) 'HALOMOD_INIT: maximum log10(M/[Msun/h]):', REAL(log10(lut%m(lut%n)))
    END IF

    lut%gmin=1.-integrate(lut%nu(1),large_nu,gnu,acc,3)
    lut%gmax=integrate(lut%nu(lut%n),large_nu,gnu,acc,3)
    lut%gbmin=1.-integrate(lut%nu(1),large_nu,gnubnu,acc,3)
    lut%gbmax=integrate(lut%nu(lut%n),large_nu,gnubnu,acc,3)
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
       WRITE(*,*) 'HALOMOD_INIT: one-halo amplitude [Mpc/h]^3:', REAL(A0)
       WRITE(*,*) 'HALOMOD_INIT: log10(M*/[Msun/h]):', REAL(log10(A0*comoving_matter_density(cosm)))
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

    IF(verbose) CALL print_halomodel_parameters(z,lut,cosm)

    !IF(verbose) verbose=.FALSE.

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
          g_lcdm=growint(ainf,acc,cos_lcdm)

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
    IF(lut%void) THEN
       STOP 'HALOMOD: Extreme caution, parameter void is defined twice'
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

       IF(lut%ibias==2) THEN
          !Only necessary for second-order bias integral
          ALLOCATE(integrand21(lut%n),integrand22(lut%n))
       END IF

       DO i=1,lut%n

          m=lut%m(i)
          nu=lut%nu(i)

          !Linear bias term
          integrand11(i)=gnu(nu)*bnu(nu)*wk(1,i)/m
          integrand12(i)=gnu(nu)*bnu(nu)*wk(2,i)/m

          IF(lut%ibias==2) THEN
             !Second-order bias term
             integrand21(i)=gnu(nu)*b2nu(nu)*wk(1,i)/m
             integrand22(i)=gnu(nu)*b2nu(nu)*wk(2,i)/m
          END IF

       END DO

       !Evaluate these integrals from the tabled values
       sum11=integrate_table(lut%nu,integrand11,lut%n,1,lut%n,3)
       sum12=integrate_table(lut%nu,integrand12,lut%n,1,lut%n,3)

       IF(lut%ip2h==0) THEN
          !Do nothing in this case
       ELSE IF(lut%ip2h==1) THEN
          !Add on the value of integral b(nu)*g(nu) assuming w=1
          sum11=sum11+lut%gbmin*halo_fraction(ih(1),m,cosm)/rhom
          sum12=sum12+lut%gbmin*halo_fraction(ih(2),m,cosm)/rhom
       ELSE IF(lut%ip2h==2) THEN
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

       IF(lut%ibias==2) THEN
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
       frac=fdamp(z,lut,cosm)

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
       win_void=rho(r,rmin,rmax,rv,rs,null,null,irho)
       win_void=win_void/normalisation(rmin,rmax,rv,rs,null,null,irho)
    ELSE       
       win_void=m*win_norm(k,rmin,rmax,rv,rs,null,null,irho)/comoving_matter_density(cosm)
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
       win_compensated_void=rho(r,rmin,rmax,rv,rs,null,null,irho)
       win_compensated_void=win_compensated_void/normalisation(rmin,rmax,rv,rs,null,null,irho)
    ELSE       
       win_compensated_void=m*win_norm(k,rmin,rmax,rv,rs,null,null,irho)/comoving_matter_density(cosm)
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
       win_DMONLY=rho(r,rmin,rmax,rv,rs,null,null,irho)
       win_DMONLY=win_DMONLY/normalisation(rmin,rmax,rv,rs,null,null,irho)
    ELSE IF(ik==1) THEN
       !Properly normalise and convert to overdensity
       win_DMONLY=m*win_norm(k,rmin,rmax,rv,rs,null,null,irho)/comoving_matter_density(cosm)
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
       dc=cosm%Dc*halo_boundgas_fraction(m,cosm)*cosm%om_m/cosm%om_b
       rss=1./(1./rs+dc/rv)
    ELSE
       STOP 'WIN_CDM: Error, imod specified incorrectly'
    END IF

    rmin=0.
    rmax=rv

    IF(ik==0) THEN
       r=k
       win_CDM=rho(r,rmin,rmax,rv,rss,null,null,irho)
       win_CDM=win_CDM/normalisation(rmin,rmax,rv,rss,null,null,irho)
    ELSE IF(ik==1) THEN
       !Properly normalise and convert to overdensity
       win_CDM=m*win_norm(k,rmin,rmax,rv,rss,null,null,irho)/comoving_matter_density(cosm)
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
       win_star=rho(r,rmin,rmax,rv,rstar,null,null,irho)
       win_star=win_star/normalisation(rmin,rmax,rv,rstar,null,null,irho)
    ELSE IF(ik==1) THEN
       !Properly normalise and convert to overdensity
       win_star=m*win_norm(k,rmin,rmax,rv,rstar,null,null,irho)/comoving_matter_density(cosm)
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
    
    LOGICAL, PARAMETER :: use_UPP=.FALSE. !Use UPP or not

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
    REAL :: r500c, rmin, rmax, a, r, alphap, b, m500c, E

    INTEGER, PARAMETER :: irho=14 !Set UPP profile 

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
       UPP=rho(r,rmin,rmax,r500c,rs,null,null,irho)
    ELSE IF(ik==1) THEN
       !win_pressure_bound=winint(k/a,a*rmax,a*r500c,a*rs,irho)
       UPP=winint(k,rmin,rmax,r500c,rs,null,null,irho)
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

    !Pre-factors from equation 4.1 in Ma et al. (2015) [eV cm^-3]
    UPP=UPP*((m500c/2.1e14)**(alphap+2./3.))*(E**(8./3.))*3.37

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
    REAL :: rho0, T0, r, alpha, gamma
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
       Gamma=cosm%Gamma
       !STOP 'WIN_BOUNDGAS: Caution, gamma not being varied in KS profile'
    ELSE IF(imod==2) THEN
       irho_density=6 !Set cored isothermal profile with beta=2/3 
       irho_pressure=irho_density !okay to use density for pressure because temperature is constant
       rmin=0.
       rmax=rv
       rb=rs
       Gamma=0. !Should probably set this to something
    ELSE        
       STOP 'WIN_BOUNDGAS: Error, imod not specified correctly'
    END IF

    IF(itype==1) THEN

       !Density profile of bound gas
       IF(ik==0) THEN
          r=k
          win_boundgas=rho(r,rmin,rmax,rv,rb,Gamma,null,irho_density)
          win_boundgas=win_boundgas/normalisation(rmin,rmax,rv,rb,Gamma,null,irho_density)
       ELSE IF(ik==1) THEN
          !Properly normalise and convert to overdensity
          win_boundgas=m*win_norm(k,rmin,rmax,rv,rb,Gamma,null,irho_density)/comoving_matter_density(cosm)
       ELSE
          STOP 'WIN_BOUNDGAS: ik not specified correctly'
       END IF

       win_boundgas=halo_boundgas_fraction(m,cosm)*win_boundgas

    ELSE IF(itype==2) THEN

       !Pressure profile of bound gas
       IF(ik==0) THEN
          r=k
          win_boundgas=rho(r,rmin,rmax,rv,rb,Gamma,null,irho_pressure)
       ELSE IF(ik==1) THEN
          !The pressure window is T(r) x rho(r), we want unnormalised, so multiply by normalisation
          win_boundgas=win_norm(k,rmin,rmax,rv,rb,Gamma,null,irho_pressure)*normalisation(rmin,rmax,rv,rb,Gamma,null,irho_pressure) 
       ELSE
          STOP 'WIN_BOUNDGAS: Error, ik not specified correctly'
       END IF

       !Calculate the value of the density profile prefactor
       !also change units from cosmological to SI
       rho0=m*halo_boundgas_fraction(m,cosm)/normalisation(rmin,rmax,rv,rb,Gamma,null,irho_density)
       rho0=rho0*msun/mpc/mpc/mpc !Overflow with REAL(4) if you use mpc**3

       !Calculate the value of the temperature prefactor
       !f=p=pac=1.
       !IF(variation) fac=param(1) !Fudge factor (turbulence?)
       !alpha=cosm%param(1)
       T0=cosm%alpha*virial_temperature(m,rv)

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
             win_freegas=rho(r,rmin,rmax,rv,rf,null,null,irho_density)
             win_freegas=win_freegas/normalisation(rmin,rmax,rv,rf,null,null,irho_density)
          ELSE IF(ik==1) THEN
             !Properly normalise and convert to overdensity
             win_freegas=m*win_norm(k,rmin,rmax,rv,rf,null,null,irho_density)/comoving_matter_density(cosm)
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
                win_freegas=rho(r,rmin,rmax,rv,rf,null,null,irho_pressure)
             ELSE IF(ik==1) THEN  
                win_freegas=win_norm(k,rmin,rmax,rv,rf,null,null,irho_pressure)*normalisation(rmin,rmax,rv,rf,null,null,irho_pressure)              
             ELSE
                STOP 'WIN_PRESSURE_FREE: Error, ik not specified correctly'
             END IF

             !Calculate the value of the density profile prefactor
             !and change units from cosmological to SI
             rho0=m*halo_freegas_fraction(m,cosm)/normalisation(rmin,rmax,rv,rf,null,null,irho_density)
             rho0=rho0*msun/mpc/mpc/mpc !Overflow with REAL(4) if you use mpc**3

             !Calculate the value of the temperature prefactor
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

  FUNCTION win_norm(k,rmin,rmax,rv,rs,p1,p2,irho)

    !Calculates the normalised spherical Fourier Transform of the density profile
    !Note that this means win_norm(k->0)=1
    !and that win must be between 0 and 1
    IMPLICIT NONE
    REAL :: win_norm
    REAL, INTENT(IN) :: rmin, rmax, k, rv, rs, p1, p2
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
          win_norm=winint(k,rmin,rmax,rv,rs,p1,p2,irho)/normalisation(rmin,rmax,rv,rs,p1,p2,irho)
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

  FUNCTION rho(r,rmin,rmax,rv,rs,p1,p2,irho)

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
    REAL, INTENT(IN) :: r, rmin, rmax, rv, rs, p1, p2 !Standard profile parameters
    INTEGER, INTENT(IN) :: irho
    REAL :: y, ct, t, c, Gamma, rt, A !Derived parameters
    REAL :: P0, c500, alpha, beta, r500 !UPP parameters
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
          !Gamma=1.18 !Recommended by Rabold (2017)
          Gamma=p1
          y=r/rs
          rho=(log(1.+y)/y)
          IF(irho==11) THEN
             !KS density profile
             rho=rho**(1./(Gamma-1.))
          ELSE IF(irho==12) THEN
             !KS temperature profile
             rho=rho
          ELSE IF(irho==13) THEN
             !KS pressure profile
             rho=rho**(Gamma/(Gamma-1.))
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

  FUNCTION winint(k,rmin,rmax,rv,rs,p1,p2,irho)

    !Calculates W(k,M)
    IMPLICIT NONE
    REAL :: winint
    REAL, INTENT(IN) :: k, rmin, rmax, rv, rs, p1, p2
    INTEGER, INTENT(IN) :: irho

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
       winint=winint_normal(rmin,rmax,k,rmin,rmax,rv,rs,p1,p2,irho,iorder,acc)
    ELSE IF(imeth==2 .OR. imeth==4 .OR. imeth==5 .OR. imeth==6 .OR. imeth==7) THEN
       IF(rmin .NE. 0.) STOP 'WININT: This cannot cope with rmin to rmax - probably could be fixed quickly'
       winint=winint_bumps(k,rmax,rv,rs,p1,p2,irho,iorder,acc,imeth)
    ELSE IF(imeth==3) THEN
       winint=winint_store(rmin,rmax,k,rmin,rmax,rv,rs,p1,p2,irho,iorder,acc)
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

  FUNCTION winint_normal(a,b,k,rmin,rmax,rv,rs,p1,p2,irho,iorder,acc)

    !Integration routine using 'normal' method to calculate the normalised halo FT
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

  FUNCTION winint_bumps(k,rmax,rv,rs,p1,p2,irho,iorder,acc,imeth)

    !Integration routine to calculate the normalised halo FT
    IMPLICIT NONE
    REAL :: winint_bumps
    REAL, INTENT(IN) :: k, rmax, rv, rs, p1, p2
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
          w=winint_normal(r1,r2,k,rmin,rmax,rv,rs,p1,p2,irho,iorder,acc)
       ELSE IF(k==0. .OR. imeth==4 .OR. (imeth==7 .AND. n<=nlim)) THEN
          w=winint_store(r1,r2,k,rmin,rmax,rv,rs,p1,p2,irho,iorder,acc)
       ELSE IF(imeth==5 .OR. imeth==6 .OR. imeth==7) THEN
          IF(i==0 .OR. i==n) THEN
             !First piece done 'normally' because otherwise /0 occurs in cubic
             !Last piece will not generally be over one full oscillation
             w=winint_store(r1,r2,k,rmin,rmax,rv,rs,p1,p2,irho,iorder,acc)
          ELSE
             IF(imeth==5) THEN
                !Linear approximation to integral between nodes - see notes
                !All results from the analytic integral of a linear polynomial vs. one sine oscillation
                rn=pi*(2*i+1)/(2.*k)
                w=(2./k**2)*winint_integrand(rn,rmin,rmax,rv,rs,p1,p2,irho)*((-1.)**i)/rn !Note there is no sinc here
             ELSE IF(imeth==6 .OR. (imeth==7 .AND. n>nlim)) THEN
                !Cubic approximation to integral between nodes - see notes
                !All results from the analytic integral of a cubic polynomial vs. one sine oscillation
                x1=r1 !Beginning
                x2=r1+1.*(r2-r1)/3. !Middle
                x3=r1+2.*(r2-r1)/3. !Middle
                x4=r2 !End
                y1=winint_integrand(x1,rmin,rmax,rv,rs,p1,p2,irho)/x1 !Note there is no sinc here
                y2=winint_integrand(x2,rmin,rmax,rv,rs,p1,p2,irho)/x2 !Note there is no sinc here
                y3=winint_integrand(x3,rmin,rmax,rv,rs,p1,p2,irho)/x3 !Note there is no sinc here
                y4=winint_integrand(x4,rmin,rmax,rv,rs,p1,p2,irho)/x4 !Note there is no sinc here
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

  FUNCTION winint_integrand(r,rmin,rmax,rv,rs,p1,p2,irho)

    !The integrand for the W(k) integral
    !Note that the sinc function is not included
    IMPLICIT NONE
    REAL :: winint_integrand
    REAL, INTENT(IN) :: r, rmin, rmax, rv, rs, p1, p2
    INTEGER, INTENT(IN) :: irho
    !REAL, PARAMETER :: rmin=0.
    !REAL, PARAMETER :: rmax=1d8 !Some large number

    IF(r==0.) THEN
       winint_integrand=4.*pi*rhor2at0(irho)
    ELSE
       winint_integrand=4.*pi*(r**2)*rho(r,rmin,rmax,rv,rs,p1,p2,irho)
    END IF

  END FUNCTION winint_integrand

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
    winnfw=4.*pi*winnfw*(rs**3.)/normalisation(rmin,rmax,rv,rs,null,null,4)

  END FUNCTION winnfw

  FUNCTION normalisation(rmin,rmax,rv,rs,p1,p2,irho)

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
    REAL, INTENT(IN) :: rmin, rmax, rv, rs, p1, p2
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
       normalisation=winint(k,rmin,rmax,rv,rs,p1,p2,irho)
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

    REAL, PARAMETER :: dc=1.686

    bps=1.+(nu**2-1.)/dc

  END FUNCTION bps

  FUNCTION bst(nu)

    !Sheth Tormen bias
    IMPLICIT NONE
    REAL :: bst
    REAL, INTENT(IN) :: nu

    REAL, PARAMETER :: p=0.3
    REAL, PARAMETER :: q=0.707
    REAL, PARAMETER :: dc=1.686

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
       M0=cosm%M0
       beta=0.6
       halo_boundgas_fraction=(cosm%om_b/cosm%om_m)/(1.+(M0/m)**beta)
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
    !crap=cosm%A

    IF(imod==1 .OR. imod==3) THEN
       !Fedeli (2014)
       !A=0.02
       !IF(variation) A=param(5)
       A=cosm%Astar
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

!!$  SUBROUTINE print_baryon_parameters(cosm)!param,param_names,n)
!!$
!!$    IMPLICIT NONE
!!$    !REAL, INTENT(IN) :: param(n)
!!$    !CHARACTER(len=256), INTENT(IN) :: param_names(n)
!!$    !INTEGER, INTENT(IN) :: n
!!$    INTEGER :: i
!!$    TYPE(cosmology), INTENT(IN) :: cosm
!!$
!!$    WRITE(*,*) 'PRINT_BARYON_PARAMETERS: Writing to screen'
!!$    WRITE(*,*) '=========================================='
!!$    DO i=1,cosm%np
!!$       WRITE(*,*) i, TRIM(cosm%param_names(i)), ':', cosm%param(i)
!!$    END DO
!!$    WRITE(*,*) '=========================================='
!!$    WRITE(*,*) 'PRINT_BARYON_PARAMETERS: Done'
!!$    WRITE(*,*)
!!$
!!$  END SUBROUTINE print_baryon_parameters

END MODULE HMx
