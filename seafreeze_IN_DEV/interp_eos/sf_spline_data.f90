MODULE sf_spline_data
    IMPLICIT NONE
    SAVE

    INTEGER, PARAMETER :: NP = 126, NT = 26, NM = 46        ! dimension of the data arrays
    REAL, DIMENSION(NP) ::  P_grid                          ! pressure grid points, MPa
    REAL, DIMENSION(NT) ::  T_grid                          ! temperature grid points, K
    REAL, DIMENSION(NM) ::  M_grid                          ! molality grid points, mol/kg

    ! These are the arrays which store the Gibbs data and its derivatives
    REAL, DIMENSION(NP, NT, NM) :: G, dG_dP, dG_dT, dG_dM, d2G_dPdT, d2G_dPdM, d2G_dTdM, d3G_dPdTdM, d2G_dT2
    REAL, DIMENSION(NT*NM) :: dummy                         ! a dummy array for reading in data

    CONTAINS
        SUBROUTINE set_grid()
        !   Subroutine for filling in the coordinate arrays.
            IMPLICIT NONE
            INTEGER :: i_p, i_t, i_m                                ! indices for the do-loops

            DO i_p = 1, NP
                P_grid(i_p) = 0.1 + (i_p - 1) * 8.      ! 0.1 MPa to 1000.1 MPa
            END DO
            DO i_t = 1, NT
                T_grid(i_t) = 250. + (i_t - 1) * 4.     ! 250 K to 350 K
            END DO
            DO i_m = 1, NM
                M_grid(i_m) = 0.0 + (i_m - 1) * 0.08   ! 0 to 3.6 
            END DO
            RETURN
        END SUBROUTINE set_grid

END MODULE sf_spline_data