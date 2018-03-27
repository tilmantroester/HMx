module HMx_setup
    use cosmology_functions, only: cosmology
    implicit none

    type HMx_setup_config
        real(8) :: kmin, kmax
        integer :: nk

        real(8) :: zmin, zmax, amin, amax, mmin, mmax
        integer :: nz

        integer :: compute_p_lin

        integer, dimension(2) :: fields

        type(cosmology) :: cosm

        real(8), dimension(:), allocatable :: k, a
        logical :: verbose
    end type HMx_setup_config
end module HMx_setup

function setup(options) result(result)
	use HMx_setup
    use HMx
    use array_operations
	use cosmosis_modules
	implicit none

    !Arguments
	integer(cosmosis_block), value :: options
    !Return value
	type(c_ptr) :: result
    !Variables
	integer(cosmosis_status) :: status
    type(HMx_setup_config), pointer :: HMx_config
    integer :: verbose
	
	allocate(HMx_config)

	if(datablock_get(options, option_section, "nz", HMx_config%nz) /= 0) then
        write(*,*) "Could not load nz."
        stop
    end if
    if(datablock_get(options, option_section, "zmin", HMx_config%zmin) /= 0) then
        write(*,*) "Could not load zmin."
        stop
    end if
    if(datablock_get(options, option_section, "zmax", HMx_config%zmax) /= 0) then
        write(*,*) "Could not load zmax."
        stop
    end if
    HMx_config%amin = 1.0/(1+HMx_config%zmax)
    HMx_config%amax = 1.0/(1+HMx_config%zmin)

    if(datablock_get(options, option_section, "nk", HMx_config%nk) /= 0) then
        write(*,*) "Could not load nk."
        stop
    end if
    if(datablock_get(options, option_section, "kmin", HMx_config%kmin) /= 0) then
        write(*,*) "Could not load kmin."
        stop
    end if
    if(datablock_get(options, option_section, "kmax", HMx_config%kmax) /= 0) then
        write(*,*) "Could not load kmax."
        stop
    end if

    status = datablock_get(options, option_section, "field1", HMx_config%fields(1))
    if(status /= 0) then
        write(*,*) "Could not load field1:", status
        stop
    end if
    if(datablock_get(options, option_section, "field2", HMx_config%fields(2)) /= 0) then
        write(*,*) "Could not load field2."
        stop
    end if

    status = datablock_get_double_default(options, option_section, "mmin", 1e7, HMx_config%mmin)
    status = datablock_get_double_default(options, option_section, "mmax", 1e17, HMx_config%mmax)

    status = datablock_get_int_default(options, option_section, "compute_p_lin", 1, HMx_config%compute_p_lin)
    status = datablock_get_int_default(options, option_section, "verbose", 1, verbose)
    HMx_config%verbose = verbose > 0

    !Create k array (log spacing)
    call fill_array(log(HMx_config%kmin), log(HMx_config%kmax), HMx_config%k, HMx_config%nk)
    HMx_config%k = exp(HMx_config%k)
    !Create a arrays. The projection module wants increasing z, so a is decreasing
    call fill_array(HMx_config%amax, HMx_config%amin, HMx_config%a, HMx_config%nz)
	
    result = c_loc(HMx_config)

end function setup

function execute(block, config) result(status)
    use cosmosis_modules
    use HMx_setup
    use HMx, only : initialise_cosmology, print_cosmology, calculate_HMx
    use constants

    implicit none
    !Arguments
    integer(cosmosis_block), value :: block
    type(c_ptr), value :: config
    !Return value
    integer(cosmosis_status) :: status
    !Variables
    type(HMx_setup_config), pointer :: HMx_config
    integer :: i
    integer, dimension(:) :: fields(2)
    real(8) :: log10_M0
    real(8), dimension(:), allocatable :: k_plin, z_plin
    real(8), dimension(:,:), allocatable :: pk_lin, pk_1h, pk_2h, pk_full
    character(len=256) :: pk_section

    call c_f_pointer(config, HMx_config)

    HMx_config%cosm%wa=0.
    HMx_config%cosm%T_cmb=2.72
    HMx_config%cosm%z_cmb=1100.

    ! Cosmology parameters
    status = datablock_get_double_default(block, cosmological_parameters_section, "omega_m", 0.3, HMx_config%cosm%om_m)
    status = datablock_get_double_default(block, cosmological_parameters_section, "omega_lambda", 1.0-HMx_config%cosm%om_m, HMx_config%cosm%om_v)
    status = datablock_get_double_default(block, cosmological_parameters_section, "omega_b", 0.05, HMx_config%cosm%om_b)
    status = datablock_get_double_default(block, cosmological_parameters_section, "omega_nu", 0.0, HMx_config%cosm%om_nu)
    status = datablock_get_double_default(block, cosmological_parameters_section, "h0", 0.7, HMx_config%cosm%h)
    status = datablock_get_double_default(block, cosmological_parameters_section, "sigma_8", 0.8, HMx_config%cosm%sig8)
    status = datablock_get_double_default(block, cosmological_parameters_section, "n_s", 0.96, HMx_config%cosm%n)
    status = datablock_get_double_default(block, cosmological_parameters_section, "w", -1.0, HMx_config%cosm%w)

    ! Baryon parameters
    status = datablock_get_double_default(block, halo_model_parameters_section, "alpha", 1.0, HMx_config%cosm%alpha)
    status = datablock_get_double_default(block, halo_model_parameters_section, "Dc", 0.0, HMx_config%cosm%Dc)
    status = datablock_get_double_default(block, halo_model_parameters_section, "Gamma", 1.18, HMx_config%cosm%Gamma)
    status = datablock_get_double_default(block, halo_model_parameters_section, "log10_M0", 14.08, log10_M0)
    status = datablock_get_double_default(block, halo_model_parameters_section, "Astar", 0.02, HMx_config%cosm%Astar)
    status = datablock_get_double_default(block, halo_model_parameters_section, "log10_whim", 6.0, HMx_config%cosm%whim)

    HMx_config%cosm%M0 = 10**log10_M0
    HMx_config%cosm%whim = 10**log10_whim

    if(HMx_config%compute_p_lin == 0) then
        status = datablock_get_double_grid(block, matter_power_lin_section, &
                                        "k_h", k_plin, &
                                        "z", z_plin, &
                                        "p_k", pk_lin)
        if(status /= 0) then
            write(*,*) "Could not load load linear power spectrum."
            stop
        end if
        HMx_config%cosm%external_plin = .true.
        HMx_config%cosm%nplin = size(k_plin)
        allocate(HMx_config%cosm%logk_logplin, source=log(k_plin))
        allocate(HMx_config%cosm%logplin, source=log(pk_lin(:,1)*k_plin**3/(2*pi**2)))
    else
        HMx_config%cosm%external_plin = .false.
    end if

    call initialise_cosmology(HMx_config%verbose, HMx_config%cosm)
    if(HMx_config%verbose) then
        call print_cosmology(HMx_config%cosm)
    end if

    call calculate_HMx(HMx_config%fields, &
                       HMx_config%mmin, HMx_config%mmax, &
                       HMx_config%k, HMx_config%nk, &
                       HMx_config%a, HMx_config%nz, &
                       pk_lin, pk_2h, pk_1h, pk_full, &
                       HMx_config%cosm, HMx_config%verbose)

    ! Remove the k^3/2pi^2 factor
    forall (i=1:HMx_config%nk) pk_full(i,:) = pk_full(i,:)*2*pi**2/HMx_config%k(i)**3
    
    ! Write power spectra to relevant section
    if(all(HMx_config%fields == (0, 0))) then
        pk_section = "matter_matter_power_spectrum"
    else if(all(HMx_config%fields == (1, 1))) then
        pk_section = "dm_dm_power_spectrum"
    else if(all(HMx_config%fields == (2, 2))) then
        pk_section = "gas_gas_power_spectrum"
    else if(all(HMx_config%fields == (3, 3))) then
        pk_section = "stars_stars_power_spectrum"
    else if(all(HMx_config%fields == (6, 6))) then
        pk_section = "pressure_pressure_power_spectrum"
    else if(any(HMx_config%fields == 0) .and. any(HMx_config%fields == 1)) then
        pk_section = "matter_dm_power_spectrum"
    else if(any(HMx_config%fields == 0) .and. any(HMx_config%fields == 2)) then
        pk_section = "matter_gas_power_spectrum"
    else if(any(HMx_config%fields == 0) .and. any(HMx_config%fields == 3)) then
        pk_section = "matter_stars_power_spectrum"
    else if(any(HMx_config%fields == 0) .and. any(HMx_config%fields == 6)) then
        pk_section = "matter_pressure_power_spectrum"
    else
        write(*,*) "Unsupported combination of fields:", HMx_config%fields
        stop
    end if

    status = datablock_put_double_grid(block, pk_section, &
                                        "k_h", HMx_config%k, &
                                        "z", 1.0/HMx_config%a-1.0, &
                                        "p_k", pk_full)
    
    if(status /= 0) then
        write(*,*) "Failed to write NL power spectrum to datablock."
    end if

end function execute

function cleanup(config) result(status)
    use cosmosis_modules
    use HMx_setup
    
    !Arguments
    type(c_ptr), value :: config
    !Return value
    integer(cosmosis_status) :: status
    !Variables
    type(HMx_setup_config), pointer :: HMx_config  

    !Free memory allocated in the setup function
    call c_f_pointer(config, HMx_config)
    !deallocate(HMx_config)

    status = 0
end function cleanup