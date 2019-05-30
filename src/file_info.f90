MODULE file_info

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: file_length
  PUBLIC :: count_number_of_lines
  
CONTAINS

  INTEGER FUNCTION file_length(file_name,verbose)

    ! Get the number of lines in the file
    USE logical_operations
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: file_name
    LOGICAL, OPTIONAL, INTENT(IN) :: verbose
    INTEGER :: n
    LOGICAL :: lexist

    !IF(verbose) WRITE(*,*) 'FILE_LENGTH: File: ', trim(file_name)
    IF(present_and_correct(verbose)) WRITE(*,*) 'FILE_LENGTH: File: ', trim(file_name)
    INQUIRE(file=file_name,exist=lexist)
    IF(.NOT. lexist) THEN
       WRITE(*,*) 'FILE_LENGTH: File: ', trim(file_name)
       STOP 'FILE_LENGTH: Error, file does not exist'
    END IF
    OPEN(7,file=file_name,status='old')

    ! Newer version that lacks 'data' seems okay
    n=0
    DO
       n=n+1
       READ(7,*, end=301)
    END DO

    ! 301 is just the label to jump to when the end of the file is reached
301 CLOSE(7)

    file_length=n-1

    !IF(verbose) THEN
    IF(present_and_correct(verbose)) THEN
       WRITE(*,*) 'FILE_LENGTH: Length:', file_length
       WRITE(*,*)
    END IF

  END FUNCTION file_length

  FUNCTION count_number_of_lines(filename) result(n)

    ! Tilman's version of file_length (nice because no GOTO)
    IMPLICIT NONE
    CHARACTER(len=*), intent(in) :: filename
    INTEGER :: n, file_unit, iostat
    CHARACTER :: c

    OPEN(newunit=file_unit, file=filename, status='old', iostat=iostat)
    IF(iostat > 0) THEN
        PRINT *, "Failed to open file ", filename
        n = -1
        RETURN
    END IF

    n = 0
    DO
        READ(unit=file_unit, fmt=*, iostat=iostat) c
        IF(iostat == 0) THEN
            IF(c == '#') THEN
                print *, "Detected comment lines in file ", filename, ". Check your data file."
                CYCLE
            END IF

            n = n + 1
        ELSE IF(iostat < 0) THEN
            !print *, "Reached end of file."
            EXIT
        ELSE
            PRINT *, "Error reading file ", filename
            n = -1
            RETURN
        END IF
    END DO

    CLOSE(file_unit)
    
  END FUNCTION count_number_of_lines

END MODULE file_info

