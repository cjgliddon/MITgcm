PROGRAM test_driver
    USE sf_spline_data
    IMPLICIT NONE

    REAL :: p_Pa, T_C, S        ! inputs: pressure, temp (oC), salinity
    REAL :: p_MPa, T_K, M       ! outputs: pressure, temp (K), molality

    WRITE (*,*) "Enter pressure in pascals: "
    READ (*,*) p_Pa
    WRITE (*,*) "Enter temperature in degrees C: "
    READ (*,*) T_C
    WRITE (*,*) "Enter salinity in g/kg: "
    READ (*,*) S

    CALL convert_PTS(p_Pa, T_C, S, p_MPa, T_K, M)

    WRITE (*,*) "The pressure in MPa is: ", p_MPa
    WRITE (*,*) "The temperature in K is: ", T_K
    WRITE (*,*) "The molality in mol/kg is: ", M

END PROGRAM test_driver

SUBROUTINE get_spline_coeffs_sf(p, T, S, a)
    IMPLICIT NONE

    ! input/output variables
    REAL, INTENT(IN) :: p, T, S         ! pressure, temperature, salinity
    REAL, DIMENSION(4, 4, 4), INTENT(OUT) :: a          ! coefficient array

    ! internal variables
    REAL :: p_MPa, T_K, M                  ! pressure in MPa, temperature in K, molality in mol/kg

    CALL convert_PTS(p, T, S, p_MPa, T_K, M)
    CALL get_spline_indices(p_MPa, T_K, M, i, j, k)

END SUBROUTINE get_spline_coeffs_sf

SUBROUTINE convert_PTS(p, T, S, p_MPa, T_K, M)
!   Subroutine for converting the p, T, S coordinates of MITgcm into the pTM coordinates used by SeaFreeze
    IMPLICIT NONE
    REAL, INTENT(IN) :: p, T, S
    REAL, INTENT(OUT) :: p_MPa, T_K, M
    REAL, PARAMETER :: mm_salt = 58.44             ! molecular mass of salt, in g/mol

    p_MPa = p / 1.0E6
    T_K = T - 273.15
    M = S / mm_salt

    RETURN
END SUBROUTINE convert_PTS

SUBROUTINE get_spline_indices(p, T, M, i, j, k)
!   Given pressure, temperature, and molality (in MPa, K, mol/kg, respectively), returns
!   the indices of the p, T, M array grid corresponding to the origin of the parallelepiped
!   used to calculating the spline.
    USE sf_spline_data
    IMPLICIT NONE

    REAL, INTENT(IN) :: p, T, M
    INTEGER, INTENT(OUT) :: i, j, k

    
END SUBROUTINE get_spline_indices