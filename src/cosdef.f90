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
