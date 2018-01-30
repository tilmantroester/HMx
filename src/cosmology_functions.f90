MODULE cosmology_functions

  !USE cosdef
  USE interpolate
  USE constants

  IMPLICIT NONE
  REAL, PARAMETER :: acc_cos=1e-4
  INTEGER, PARAMETER :: ibox=0 !Consider the simulation volume
  REAL, PARAMETER :: Lbox=400. !Simulation box size

  !Contains cosmological parameters that need only be calculated once
  TYPE cosmology     
     REAL :: om_m, om_b, om_v, om_c, h, n, sig8, w, wa, om_nu
     REAL :: om, k, z_cmb, om_r, T_cmb
     REAL :: A
     REAL, ALLOCATABLE :: logsigma(:), logr_logsigma(:)
     REAL, ALLOCATABLE :: growth(:), a_growth(:)
     REAL, ALLOCATABLE :: r(:), a_r(:)
     REAL, ALLOCATABLE :: logplin(:), logk_logplin(:) !Added for input linear Pk
     INTEGER :: nsig, ng, nr, nplin
     CHARACTER(len=256) :: name = ""
     LOGICAL :: external_plin
     !Varying baryon parameters
     !INTEGER :: np=5
     !REAL :: param(5), param_defaults(5), param_min(5), param_max(5)
     REAL :: alpha, Dc, Gamma, M0, Astar
     !CHARACTER(len=256) :: param_names(5)
     !LOGICAL :: param_log(5)
  END TYPE cosmology

CONTAINS

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

  FUNCTION AH(z,cosm)

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

  FUNCTION Omega_m(z,cosm)

    !This calculates Omega_m variations with z!
    IMPLICIT NONE
    REAL :: Omega_m
    REAL, INTENT(IN) :: z
    TYPE(cosmology), INTENT(IN) :: cosm

    Omega_m=(cosm%om_m*(1.+z)**3)/Hubble2(z,cosm)

  END FUNCTION Omega_m

  FUNCTION X_de(a,cosm)

    IMPLICIT NONE
    REAL :: X_de
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(IN) :: cosm

    !The time evolution for Om_w for w(a) DE models
    X_de=(a**(-3.*(1.+cosm%w+cosm%wa)))*exp(-3.*cosm%wa*(1.-a))

  END FUNCTION X_de

  FUNCTION w_de(a,cosm)

    !w(a) for DE models
    IMPLICIT NONE
    REAL :: w_de
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(IN) :: cosm

    w_de=cosm%w+(1.-a)*cosm%wa

  END FUNCTION w_de

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

  FUNCTION redshift_r(r,cosm)

    IMPLICIT NONE
    REAL :: redshift_r
    REAL, INTENT(IN) :: r
    TYPE(cosmology), INTENT(IN) :: cosm

    redshift_r=redshift_a(find(r,cosm%r,cosm%a_r,cosm%nr,3,3,2))

  END FUNCTION redshift_r

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

  FUNCTION cosmic_distance(z,cosm)

    IMPLICIT NONE
    REAL :: cosmic_distance
    REAL, INTENT(IN) :: z
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: a

    a=scale_factor_z(z)

    cosmic_distance=find(a,cosm%a_r,cosm%r,cosm%nr,3,3,2)

  END FUNCTION cosmic_distance

  SUBROUTINE print_cosmology(cosm)

    IMPLICIT NONE
    TYPE(cosmology) :: cosm

    WRITE(*,*) 'COSMOLOGY: ', TRIM(cosm%name)
    WRITE(*,*) '===================================='
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
    WRITE(*,*) '===================================='
    WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'alpha:', cosm%alpha
    WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'Dc:', cosm%Dc
    WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'Gamma:', cosm%Gamma
    WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'log10(M0):', log10(cosm%M0)
    WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'Astar:', cosm%Astar
    WRITE(*,*) '===================================='
    WRITE(*,*)

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

    !Boring default cosmology
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
    cosm%external_plin=.FALSE.

    !Default values of baryon parameters
    cosm%alpha=1.
    cosm%Dc=0.
    cosm%Gamma=1.18
    cosm%M0=1.2e14
    cosm%Astar=0.02

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

  SUBROUTINE initialise_cosmology(verbose,cosm)

    IMPLICIT NONE
    TYPE(cosmology), INTENT(INOUT) :: cosm
    LOGICAL, INTENT(IN) :: verbose
    REAL :: sigi    
    
    !cosm%name = ''
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
    CALL fill_growtab(verbose,cosm)

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
    CALL fill_sigtab(verbose,cosm)

  END SUBROUTINE initialise_cosmology

  SUBROUTINE initialise_distances(verbose,cosm)

    !Fill up tables of a vs. r(a) (comoving distance)
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: verbose
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
       cosm%r(i)=integrate_distance(cosm%a_r(i),1.,acc_cos,3,cosm)
    END DO
    IF(verbose) THEN
       WRITE(*,*) 'INITIALISE_DISTANCE: minimum r [Mpc/h]:', REAL(cosm%r(cosm%nr))
       WRITE(*,*) 'INITIALISE_DISTANCE: maximum r [Mpc/h]:', REAL(cosm%r(1))
    END IF

    !Find the horizon distance in your cosmology
    rh=integrate_distance(0.,1.,acc_cos,3,cosm)
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

  SUBROUTINE random_cosmology(cosm)

    USE random_numbers
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
    
    REAL, PARAMETER :: kmax=1e8

    IF(k==0.) THEN
       !If p_lin happens to be foolishly called for 0 mode
       !This call should never happen, but may in integrals
       p_lin=0.
    ELSE IF(k>kmax) THEN
       !Avoids some issues if p_lin is called for very (absurdly) high k values
       !For some reason crashes can occur if this is the case
       p_lin=0.
    ELSE IF(ibox==1 .AND. k<2.*pi/Lbox) THEN
       !If investigating effects caused by a finite box size
       p_lin=0.
    ELSE
       IF(cosm%external_plin) THEN
          p_lin=cosm%A**2*exp(find(log(k),cosm%logk_logplin,cosm%logplin,cosm%nplin,3,3,2))
       ELSE
          !In this case get the power from the transfer function
          p_lin=(cosm%A**2)*(Tk(k,cosm)**2)*(k**(cosm%n+3.))
       END IF
       !'Grow' the power from z=0 to the redshift of interest
       p_lin=p_lin*grow(z,cosm)**2
    END IF

  END FUNCTION p_lin

  SUBROUTINE fill_sigtab(verbose,cosm)

    !This fills up tables of r vs. sigma(r) across a range in r!
    !It is used only in look-up for further calculations of sigma(r) and not otherwise!
    !and prevents a large number of calls to the sigint functions
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: verbose
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
    IF(ALLOCATED(cosm%logr_logsigma)) DEALLOCATE(cosm%logr_logsigma)
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
    ALLOCATE(cosm%logr_logsigma(nsig),cosm%logsigma(nsig))

    cosm%logr_logsigma=log(rtab)
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

    INTEGER, PARAMETER :: iorder=3
    REAL, PARAMETER :: rsplit=1e-2

    IF(r>=rsplit) THEN
       sigma=sqrt(sigint0(r,z,cosm,acc_cos,iorder))
    ELSE IF(r<rsplit) THEN
       sigma=sqrt(sigint1(r,z,cosm,acc_cos,iorder)+sigint2(r,z,cosm,acc_cos,iorder))
    ELSE
       STOP 'SIGMA: Error, something went wrong'
    END IF

  END FUNCTION sigma

  FUNCTION sigma_integrand(k,R,z,cosm)

    !The integrand for the sigma(R) integrals
    USE special_functions
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
    USE special_functions
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

    IMPLICIT NONE
    REAL :: sigma_cb
    REAL, INTENT(IN) :: r, z
    TYPE(cosmology), INTENT(IN) :: cosm

    !Finds sigma_cold from look-up tables
    !In this version sigma_cold=sigma

    !sigma_cb=grow(z,cosm)*exp(find(log(r),log(cosm%r_sigma),log(cosm%sigma),cosm%nsig,3,3,2))
    sigma_cb=grow(z,cosm)*exp(find(log(r),cosm%logr_logsigma,cosm%logsigma,cosm%nsig,3,3,2))

  END FUNCTION sigma_cb

   FUNCTION grow(z,cosm)

    !Scale-independent growth function | g(z=0)=1
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

  FUNCTION growint(a,acc,cosm)

    !Integrates between a and b until desired accuracy is reached
    !Stores information to reduce function calls
    IMPLICIT NONE
    REAL :: growint
    REAL, INTENT(IN) :: a
    REAL, INTENT(IN) :: acc
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

  FUNCTION dispint(R,z,acc,cosm)

    !Integrates between a and b until desired accuracy is reached
    !Stores information to reduce function calls
    IMPLICIT NONE
    REAL :: dispint
    REAL, INTENT(IN) :: z, R, acc
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: a, b
    INTEGER :: i, j
    INTEGER :: n
    REAL :: x, dx
    REAL :: f1, f2, fx
    DOUBLE PRECISION :: sum_n, sum_2n, sum_new, sum_old

    INTEGER, PARAMETER :: jmin=5
    INTEGER, PARAMETER :: jmax=30
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
    USE special_functions
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

  SUBROUTINE fill_growtab(verbose,cosm)

    !Fills a table of the growth function vs. a
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: verbose
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: i
    REAL :: a, norm
    REAL, ALLOCATABLE :: d_tab(:), v_tab(:), a_tab(:)
    REAL :: ainit, amax, dinit, vinit
    
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
    CALL ode_growth(d_tab,v_tab,a_tab,0.,ainit,amax,dinit,vinit,acc_cos,3,cosm)
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

END MODULE cosmology_functions
