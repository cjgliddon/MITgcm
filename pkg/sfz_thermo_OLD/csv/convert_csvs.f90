PROGRAM convert_csvs
    IMPLICIT NONE
    INTEGER :: ios
    INTEGER, PARAMETER :: NP = 126, NT = 26, NM = 46        ! dimension of the data arrays
    REAL, DIMENSION(NT*NM) :: dummy                         ! a dummy array for reading in data

    OPEN (UNIT=8, FILE='G.csv', STATUS='OLD', ACTION='READ', IOSTAT=ios)
    openif: IF (ios == 0) THEN
        ! successfully opened file
        readloop: DO
            READ (8,*, IOSTAT=ios) dummy
            IF (ios /= 0 ) EXIT
            WRITE (*,*) "dummy array first index : ", dummy(1)
        END DO readloop
    END IF openif
    CLOSE (UNIT=8)
END PROGRAM convert_csvs