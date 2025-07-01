PROGRAM subroutines_test
    ! test driver program for the seafreeze subroutines
    IMPLICIT NONE

END PROGRAM subroutines_test

SUBROUTINE get_spline_coeffs_sf(p, T, S, a)
!   Accepts a pressure, temperature, and salinity, and then calculates a
!   4 x 4 x 4 array of coefficients for the tricubic spline interpolation as 
!   described in Lekien & Marsden
    IMPLICIT NONE

    REAL, INTENT(IN) :: p, T, S
    REAL, DIMENSION(4, 4, 4), INTENT(OUT) :: a

    
END SUBROUTINE get_spline_coeffs_sf

SUBROUTINE rescale_pTS(p, T, S, pp, TT, SS)
!   Subroutine for "rescaling" the pressure, temperature, and salinity coordinates
!   into the form used by 
    IMPLICIT NONE
END SUBROUTINE rescale_pTS

SUBROUTINE find_rho_sf(p, T, S, rho)
!   Subroutine which calculates rho at the given pressure, temperature, and
!   salinity using a tricubic interpolation of the SeaFreeze data.
    IMPLICIT NONE

    REAL, INTENT(IN) :: p, T, S     ! pressure, temperature, salinity
    REAL, INTENT(OUT) :: rho        ! density

    INTEGER :: i, j, k              ! indices of the coefficient array
    REAL, DIMENSION(4, 4, 4) :: a   ! coefficient array
    REAL :: pp, TT, SS              ! transformed variables used for calculating rho
    REAL :: sum0                    ! 

    ! Get the coefficients of the data
    CALL get_spline_coeffs_sf(p, T, S, a)
    ! Transform coordinates
    CALL rescale_pTS(p, T, S, pp, TT, SS)
    ! finally, calculate rho. Iterate over coefficient array to calculate terms of rho
    DO i = 1, 4
        DO j = 1, 4
            DO k = 1, 4
                sum0 = sum0 + a(i, j, k) * i * (pp**(i - 1)) * (TT**j) * (SS**k)
            END DO
        END DO
    END DO
    rho = 1.0/sum0

    RETURN 
END SUBROUTINE find_rho_sf