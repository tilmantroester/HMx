module HMx_wrapper
    use iso_c_binding,  only : c_int, c_double, c_float, c_bool, c_loc, c_ptr, c_f_pointer
    use HMx, only: calculate_HMx, halomod, assign_halomod, init_halomod
    use cosmology_functions, only: cosmology, assign_cosmology, init_cosmology, norm_none, itk_external
    use constants

    implicit none

    contains
        function run_HMx(Om_m, Om_b, Om_w, m_nu, &
                           h, ns, sigma8, w, wa, &
                           halo_model_mode, Theat, eta0, As, &
                           nf, ifield, &
                           nk, k, &
                           na, a, &
                           nk_a_plin, k_plin, na_a_plin, a_plin, &
                           nk_plin, na_plin, pk_lin, &
                        !    nq, q, &
                        !    nk_Tcold, na_Tcold, Tcold, &
                           nf1_pk, nf2_pk, nk_pk, na_pk, pk_hmx, &
                           verbose) result(status) bind(c, name="run_HMx") 
            
            REAL(kind=c_double), INTENT(IN) :: Om_m, Om_b, Om_w, m_nu        ! Primary parameters
            REAL(kind=c_double), INTENT(IN) :: h, ns, sigma8, w, wa           ! Primary parameters
            
            INTEGER(kind=c_int), INTENT(IN) :: halo_model_mode              ! Halo model mode
            REAL(kind=c_double), INTENT(IN) :: Theat, eta0, As               ! Halo model parameters

            INTEGER(kind=c_int), INTENT(IN) :: nf                       ! Number of different fields
            INTEGER(kind=c_int), INTENT(IN) :: ifield(nf)               ! Indices for different fields
            INTEGER(kind=c_int), INTENT(IN) :: nk                       ! Number of requested k points
            REAL(kind=c_double), INTENT(IN) :: k(nk)                    ! Requested k array [h/Mpc]
            INTEGER(kind=c_int), INTENT(IN) :: na                       ! Number of requested a points
            REAL(kind=c_double), INTENT(IN) :: a(na)                    ! Requested a array
            INTEGER(kind=c_int), INTENT(IN) :: nk_a_plin, na_a_plin     ! Number of k, a arrays for p_lin
            REAL(kind=c_double), INTENT(IN) :: k_plin(nk_a_plin), a_plin(na_a_plin)  ! plin k and a arrays
            INTEGER(kind=c_int), INTENT(IN) :: nk_plin, na_plin         ! Number of p_lin points
            REAL(kind=c_double), INTENT(IN) :: pk_lin(nk_plin, na_plin) ! p_lin array
            ! INTEGER(kind=c_int), INTENT(IN) :: nq         ! Number of q points
            ! REAL(kind=c_double), INTENT(IN) :: q(nq)         ! q array [h/Mpc]
            ! INTEGER(kind=c_int), INTENT(IN) :: nk_Tcold, na_Tcold        ! Number of Tcold points
            ! REAL(kind=c_double), INTENT(IN) :: Tcold(nk_Tcold, na_Tcold)         ! Tcold array
            INTEGER(kind=c_int), INTENT(IN) :: nf1_pk, nf2_pk, nk_pk, na_pk
            REAL(kind=c_double), INTENT(OUT) :: pk_hmx(nf1_pk, nf2_pk, nk_pk, na_pk) ! Pow(f1,f2,k,a)

            
            LOGICAL(kind=c_bool), INTENT(IN) :: verbose

            INTEGER(kind=c_int) :: status
            
            INTEGER :: icosmo, ihm, i
            TYPE(halomod) :: hmod   ! Halo model
            TYPE(cosmology) :: cosm ! Cosmology


            ! Working arrays
            REAL, ALLOCATABLE :: pow_li(:, :)       ! Pow(k,a)
            REAL, ALLOCATABLE :: pow_2h(:, :, :, :) ! Pow(f1,f2,k,a)
            REAL, ALLOCATABLE :: pow_1h(:, :, :, :) ! Pow(f1,f2,k,a)
            REAL, ALLOCATABLE :: pow_hm(:, :, :, :) ! Pow(f1,f2,k,a)

            if(nk /= nk_pk) then
                write(*,*) "Inconsistent requested k array sizes:", nk, nk_pk
                status = 1
                return
            end if

            if(na /= na_pk) then
                write(*,*) "Inconsistent requested a array sizes:", na, na_pk
                status = 1
                return
            end if

            if(nk_a_plin /= nk_plin) then
                write(*,*) "Inconsistent plin k array sizes:", nk_a_plin, nk_plin
                status = 1
                return
            end if

            if(na_a_plin /= na_plin) then
                write(*,*) "Inconsistent plin a array sizes:", na_a_plin, na_plin
                status = 1
                return
            end if

            ! write(*,*) "Input:"
            ! write(*,*) "k:", k
            ! write(*,*) "a:", a
            ! write(*,*) "pk_lin:", pk_lin
            ! write(*,*) "k**3", k**3
            ! write(*,*) "pi:", pi, (2*pi**2)
            ! write(*,*) "Delta_lin(a=1):", pk_lin(:,na_plin)*k**3/(2*pi**2)
            ! STOP

            ! For cosmology
            icosmo = 1 ! 'boring' cosmology
            CALL assign_cosmology(icosmo, cosm, LOGICAL(verbose))

            cosm%Om_m = Om_m
            cosm%Om_b = Om_b
            cosm%Om_w = Om_w
            cosm%Om_v = 0.0

            cosm%m_nu = m_nu

            cosm%h = h
            cosm%ns = ns
            cosm%w = w
            cosm%wa = wa

            cosm%sig8 = sigma8
            cosm%norm_method = norm_none
            ! Copy in linear power
            cosm%itk = itk_external
            cosm%has_power = .true.
            cosm%nk_plin = nk_plin
            cosm%na_plin = na_plin

            if(.not. allocated(cosm%log_k_plin)) allocate(cosm%log_k_plin(nk_a_plin))
            if(.not. allocated(cosm%log_a_plin)) allocate(cosm%log_a_plin(na_a_plin))
            if(.not. allocated(cosm%log_plin)) allocate(cosm%log_plin(nk_plin))
            if(.not. allocated(cosm%log_plina)) allocate(cosm%log_plina(nk_plin, na_plin))
            cosm%log_k_plin = log(k_plin)
            cosm%log_plin = log(pk_lin(:,na_plin)*k_plin**3/(2*pi**2))
            cosm%log_a_plin = log(a_plin)
            forall (i=1:nk_plin) cosm%log_plina(i,:) = log(pk_lin(i,:)*k_plin(i)**3/(2*pi**2))

            ! if(m_nu > 0.0) then
            !     ! Need to deal with scale-depedent growth
            !     if(nq /= nk_Tcold) then
            !         write(*,*) "Inconsistent q (Tcold) array sizes:", nq, nk_Tcold
            !         status = 1
            !         return
            !     end if
    
            !     if(na /= na_Tcold) then
            !         write(*,*) "Inconsistent a (Tcold)  array sizes:", na, na_Tcold
            !         status = 1
            !         return
            !     end if
            !     if(.not. allocated(cosm%log_k_Tcold)) allocate(cosm%log_k_Tcold(nq))
            !     if(.not. allocated(cosm%Tcold)) allocate(cosm%Tcold(nk_Tcold, na_Tcold))
            !     cosm%log_k_Tcold = log(q)
            !     cosm%Tcold = Tcold
            ! end if

            ! write(*,*) "cosm structure:"
            ! write(*,*) "cosm%log_k_plin:", cosm%log_k_plin
            ! write(*,*) "cosm%log_plin:", cosm%log_plin
            ! write(*,*) "cosm%log_a_plin:", cosm%log_a_plin
            ! write(*,*) "cosm%log_plina:", cosm%log_plina

            cosm%Theat = Theat

            CALL init_cosmology(cosm)

            ! For halo model
            ihm = halo_model_mode
            CALL assign_halomod(ihm, hmod, LOGICAL(verbose))

            hmod%eta0 = eta0
            hmod%As = As

            write(*,*) "Internal sigma8:", cosm%sig8

            call calculate_HMx(ifield, nf, k, nk, a, na, pow_li, pow_2h, pow_1h, pow_hm, hmod, cosm, LOGICAL(verbose))
            write(*,*) "pow_hm shape:", shape(pow_hm)

            forall (i=1:nk_pk) pk_hmx(:,:,i,:) = pow_hm(:,:,i,:)/(k(i)**3/(2*pi**2))



            status = 0
        end function

end module