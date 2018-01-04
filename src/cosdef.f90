MODULE cosdef

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
     CHARACTER(len=256) :: name
     LOGICAL :: external_plin
     !Varying baryon parameters
     INTEGER :: np=5
     REAL :: param(5), param_defaults(5), param_min(5), param_max(5)
     CHARACTER(len=256) :: param_names(5)
     LOGICAL :: param_log(5)
  END TYPE cosmology

  !Halo-model stuff that needs to be recalculated for each new z
  TYPE tables     
     REAL, ALLOCATABLE :: c(:), rv(:), nu(:), sig(:), zc(:), m(:), rr(:), sigf(:)
     REAL, ALLOCATABLE :: r500(:), m500(:), c500(:), r200(:), m200(:), c200(:)
     REAL, ALLOCATABLE :: r500c(:), m500c(:), c500c(:), r200c(:), m200c(:), c200c(:)
     REAL, ALLOCATABLE :: log_m(:)
     REAL :: sigv, sigv100, c3, knl, rnl, neff, sig8z
     REAL :: gmin, gmax, gbmin, gbmax
     INTEGER :: n
  END TYPE tables

  !Projection quantities that need to be calculated only once
  !These relate to the Limber integrals
  TYPE projection    
     REAL, ALLOCATABLE :: X(:), r_X(:)
     INTEGER :: nX
  END TYPE projection

  !Quantities that are necessary for lensing specifically
  !Possibly this could usefully be merged with the projection type
  TYPE lensing
     REAL, ALLOCATABLE :: q(:), r_q(:)
     REAL, ALLOCATABLE :: nz(:), z_nz(:)
     INTEGER :: nq, nnz
  END TYPE lensing

END MODULE cosdef
