MODULE cosmology_functions

  USE cosdef
  USE interpolate

  IMPLICIT NONE

CONTAINS

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

  FUNCTION Omega_m(z,cosm)

    !This calculates Omega_m variations with z!
    IMPLICIT NONE
    REAL :: Omega_m
    REAL, INTENT(IN) :: z
    TYPE(cosmology), INTENT(IN) :: cosm

    Omega_m=(cosm%om_m*(1.+z)**3)/Hubble2(z,cosm)

  END FUNCTION Omega_m

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

  FUNCTION cosmic_distance(z,cosm)

    IMPLICIT NONE
    REAL :: cosmic_distance
    REAL, INTENT(IN) :: z
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: a

    a=scale_factor_z(z)

    cosmic_distance=find(a,cosm%a_r,cosm%r,cosm%nr,3,3,2)

  END FUNCTION cosmic_distance

END MODULE cosmology_functions
