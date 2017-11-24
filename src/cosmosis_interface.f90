module HMx_setup
    implicit none

    type HMx_setup_config
        real(8) :: kmin, kmax
        integer :: nk

        real(8) :: zmin, zmax, amin, amax
        integer :: nz

        real(8), dimension(:), allocatable :: k, a
    end type HMx_setup_config
end module HMx_setup

function setup(options) result(result)
	use HMx_setup
    use HMx
	use cosmosis_modules
	implicit none

    !Arguments
	integer(cosmosis_block), value :: options
    !Return value
	type(c_ptr) :: result
    !Variables
	integer(cosmosis_status) :: status
	type(HMx_setup_config), pointer :: HMx_config
	
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

    call init_HMx()
    !Create k array (log spacing)
    call fill_table(log(HMx_config%kmin), log(HMx_config%kmin), HMx_config%k, HMx_config%nk)
    HMx_config%k = exp(HMx_config%k)
    !Create a arrays
    call fill_table(HMx_config%amin, HMx_config%amax, HMx_config%a, HMx_config%nz)
	
    result = c_loc(HMx_config)

end function setup

function execute(block, config) result(status)
    use cosmosis_modules
    use HMx_setup
    use HMx
    use cosdef

    implicit none
    !Arguments
    integer(cosmosis_block), value :: block
    type(c_ptr), value :: config
    !Return value
    integer(cosmosis_status) :: status
    !Variables
    type(HMx_setup_config), pointer :: HMx_config
    !type(cosmology) :: cosm
    integer, dimension(:) :: fields(2)
    real(8), dimension(:,:), allocatable :: pk_lin, pk_1h, pk_2h, pk_full

    call c_f_pointer(config, HMx_config)

    cosm%wa=0.
    cosm%T_cmb=2.72
    cosm%z_cmb=1100.

    if(datablock_get(block, option_section, "field1", fields(1)) /= 0) then
        write(*,*) "Could not load field1."
        stop
    end if
    if(datablock_get(block, option_section, "field2", fields(2)) /= 0) then
        write(*,*) "Could not load field2."
        stop
    end if

    status = datablock_get_double_default(block, cosmological_parameters_section, "omega_m", 0.3d0, cosm%om_m)
    status = datablock_get_double_default(block, cosmological_parameters_section, "omega_lambda", 1.d0-cosm%om_m, cosm%om_v)
    status = datablock_get_double_default(block, cosmological_parameters_section, "omega_b", 0.05d0, cosm%om_b)
    status = datablock_get_double_default(block, cosmological_parameters_section, "omega_nu", 0.0d0, cosm%om_nu)
    status = datablock_get_double_default(block, cosmological_parameters_section, "h0", 0.7d0, cosm%h)
    status = datablock_get_double_default(block, cosmological_parameters_section, "sigma_8", 0.8d0, cosm%sig8)
    status = datablock_get_double_default(block, cosmological_parameters_section, "n_s", 0.96d0, cosm%n)
    status = datablock_get_double_default(block, cosmological_parameters_section, "w", -1.0d0, cosm%w)

    call initialise_cosmology(cosm)
    call print_cosmology(cosm)

    call calculate_HMx(fields, &
                       HMx_config%k, HMx_config%nk, &
                       HMx_config%a, HMx_config%nz, &
                       pk_lin, pk_2h, pk_1h, pk_full, &
                       cosm)

    status = datablock_put_double_grid(block, matter_power_nl_section, &
                                       "k_h", HMx_config%k, &
                                       "z", 1.0-1.0/HMx_config%a, &
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