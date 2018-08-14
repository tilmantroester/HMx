PROGRAM test_gas_gas

    USE HMx
    USE cosmology_functions

    IMPLICIT NONE

    ! Parameter definitions
    REAL, ALLOCATABLE :: k(:)
    REAL, ALLOCATABLE :: pow_lin(:), pow_2h(:), pow_1h(:), pow_full(:)
    REAL :: kmin, kmax
    REAL :: z
    INTEGER :: ihm, icosmo, nk
    TYPE(cosmology) :: cosm
    TYPE(halomod) :: hmod
    CHARACTER(len=256) :: outfile


    ! Halo-model Parameters
    LOGICAL, PARAMETER :: verbose=.TRUE. ! Verbosity
    REAL, PARAMETER :: mmin=1e7 ! Minimum halo mass for the calculation
    REAL, PARAMETER :: mmax=1e17 ! Maximum halo mass for the calculation

  
    !Set number of k points and k range (log spaced)
    nk=128
    kmin=1e-3
    kmax=1e2
    CALL fill_array(log(kmin),log(kmax),k,nk)
    k=exp(k)
    ALLOCATE(pow_lin(nk),pow_2h(nk),pow_1h(nk),pow_full(nk))

    !Assigns the cosmological model
    icosmo = 1
    CALL assign_cosmology(icosmo,cosm,verbose)
    CALL init_cosmology(cosm)
    CALL print_cosmology(cosm)

    !Sets the redshift
    z=0.

    !Initiliasation for the halomodel calcualtion
    ihm = 6
    CALL assign_halomod(ihm,hmod,verbose)
    CALL init_halomod(mmin,mmax,z,hmod,cosm,verbose)

    !Do the halo-model calculation
    CALL calculate_halomod(2,2,k,nk,z,pow_lin,pow_2h,pow_1h,pow_full,hmod,cosm,verbose)

    !Write out the results
    outfile='tests/output/power_gas_gas.dat'
    CALL write_power(k,pow_lin,pow_2h,pow_1h,pow_full,nk,outfile,verbose)

    CONTAINS
        SUBROUTINE write_power(k,pow_lin,pow_2h,pow_1h,pow,nk,output,verbose)

            IMPLICIT NONE
            CHARACTER(len=*), INTENT(IN) :: output
            INTEGER, INTENT(IN) :: nk
            REAL, INTENT(IN) :: k(nk), pow_lin(nk), pow_2h(nk), pow_1h(nk), pow(nk)
            LOGICAL, INTENT(IN) :: verbose
            INTEGER :: i

            IF(verbose) WRITE(*,*) 'WRITE_POWER: Writing power to ', TRIM(output)

            !Loop over k values
            !Fill the tables with one- and two-halo terms as well as total
            OPEN(7,file=output)
            DO i=1,nk       
            WRITE(7,fmt='(5ES20.10)') k(i), pow_lin(i), pow_2h(i), pow_1h(i), pow(i)
            END DO
            CLOSE(7)

            IF(verbose) THEN
            WRITE(*,*) 'WRITE_POWER: Done'
            WRITE(*,*)
            END IF

        END SUBROUTINE write_power

END PROGRAM test_gas_gas