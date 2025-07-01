MODULE NaClaq_new_TD
    ! This module contains data used to compute the SeaFreeze tricubic interpolation.
    IMPLICIT none
    SAVE

    INTEGER :: i_p, i_t, i_m                                ! indices for the do-loops
    INTEGER :: rc
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

        SUBROUTINE read_sf_data()
        !   Subroutine for reading in the Gibbs data.
            ! Gibbs function
            OPEN (ACTION='read', FILE='csv/G.csv', IOSTAT=rc, NEWUNIT=20)
            outer: DO i_p = 1, NP
                read (20, *, IOSTAT=rc) dummy(:)        ! read data into dummy array
                inner: DO i_t = 1, NT
                    G(i_p, i_t, :) = dummy( (i_t - 1)*46 : i_t*46 )     ! transfer from dummy array to formal storage
                END DO inner
            END DO outer
            CLOSE (20)

            ! dG_dP
            OPEN (ACTION='read', FILE='csv/dGdP.csv', IOSTAT=rc, NEWUNIT=21)
            outer: DO i_p = 1, NP
                read (21, *, IOSTAT=rc) dummy(:)
                inner: DO i_t = 1, NT
                    dG_dP(i_p, i_t, :) = dummy( (i_t - 1)*46 : i_t*46 )
                END DO inner
            END DO outer
            CLOSE (21)

            ! dG_dT
            OPEN (ACTION='read', FILE='csv/dGdT.csv', IOSTAT=rc, NEWUNIT=22)
            outer: DO i_p = 1, NP
                read (22, *, IOSTAT=rc) dummy(:)
                inner: DO i_t = 1, NT
                    dG_dT(i_p, i_t, :) = dummy( (i_t - 1)*46 : i_t*46 )
                END DO inner
            END DO outer
            CLOSE (22)

            ! dG_dM
            OPEN (ACTION='read', FILE='csv/dGdM.csv', IOSTAT=rc, NEWUNIT=23)
            outer: DO i_p = 1, NP
                read (23, *, IOSTAT=rc) dummy(:)
                inner: DO i_t = 1, NT
                    dG_dM(i_p, i_t, :) = dummy( (i_t - 1)*46 : i_t*46 )
                END DO inner
            END DO outer
            CLOSE (23)

            ! now the mixed partials...
            OPEN (ACTION='read', FILE='csv/d2GdPT.csv', IOSTAT=rc, NEWUNIT=24)
            outer: DO i_p = 1, NP
                read (24, *, IOSTAT=rc) dummy(:)
                inner: DO i_t = 1, NT
                    d2G_dPdT(i_p, i_t, :) = dummy( (i_t - 1)*46 : i_t*46 )
                END DO inner
            END DO outer
            CLOSE (24)

            OPEN (ACTION='read', FILE='csv/d2GdPM.csv', IOSTAT=rc, NEWUNIT=25)
            outer: DO i_p = 1, NP
                read (25, *, IOSTAT=rc) dummy(:)
                inner: DO i_t = 1, NT
                    d2G_dPdM(i_p, i_t, :) = dummy( (i_t - 1)*46 : i_t*46 )
                END DO inner
            END DO outer
            CLOSE (25)

            OPEN (ACTION='read', FILE='csv/d2GdTM.csv', IOSTAT=rc, NEWUNIT=26)
            outer: DO i_p = 1, NP
                read (26, *, IOSTAT=rc) dummy(:)
                inner: DO i_t = 1, NT
                    d2G_dTdM(i_p, i_t, :) = dummy( (i_t - 1)*46 : i_t*46 )
                END DO inner
            END DO outer
            CLOSE (26)

            ! the tertiary mixed partial
            OPEN (ACTION='read', FILE='csv/d3GdPTM.csv', IOSTAT=rc, NEWUNIT=27)
            outer: DO i_p = 1, NP
                read (27, *, IOSTAT=rc) dummy(:)
                inner: DO i_t = 1, NT
                    d3G_dPdTdM(i_p, i_t, :) = dummy( (i_t - 1)*46 : i_t*46 )
                END DO inner
            END DO outer
            CLOSE (27)

            ! we also want to load d2G_dT2 for calculating theta/adiabatic temp. gradient
            OPEN (ACTION='read', FILE='csv/d2GdT2.csv', IOSTAT=rc, NEWUNIT=28)
            outer: DO i_p = 1, NP
                read (28, *, IOSTAT=rc) dummy(:)
                inner: DO i_t = 1, NT
                    d2G_dT2(i_p, i_t, :) = dummy( (i_t - 1)*46 : i_t*46 )
                END DO inner
            END DO outer
            CLOSE (28)
        END SUBROUTINE read_sf_data

END MODULE NaClaq_new_TD


!   subroutine 

SUBROUTINE get_spline_coeffs_sf(p, T, S, a)
!   Accepts a pressure, temperature, and salinity, and then calculates a
!   4 x 4 x 4 array of coefficients for the tricubic spline interpolation as 
!   described in Lekien & Marsden
!
    IMPLICIT NONE
END SUBROUTINE get_spline_coeffs_sf

SUBROUTINE rescale_pTS(p, T, S, pp, TT, SS)
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