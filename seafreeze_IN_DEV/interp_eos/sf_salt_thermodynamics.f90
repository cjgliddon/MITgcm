MODULE sf_salt_thermodynamics
!   This module contains routines.
    IMPLICIT NONE

    REAL, PARAMETER :: mm_salt = 58.44                  ! molecular mass of salt, in g/mol
    INTEGER, PARAMETER :: NP = 21, NT = 21, NM = 21        ! dimensions of the data arrays
    REAL(8), PARAMETER :: P0 = 0.0, T0 = 250., M0 = 0.0        ! "initial" values of pressure, temp., molality
    REAL(8), PARAMETER :: dP = 50., dT = 5., dM = 0.18          ! grid spacing of pressure, temperature, molality

    REAL(8), DIMENSION(NP) ::  P_grid                          ! pressure grid points, MPa
    REAL(8), DIMENSION(NT) ::  T_grid                          ! temperature grid points, K
    REAL(8), DIMENSION(NM) ::  M_grid                          ! molality grid points, mol/kg
    LOGICAL :: grid_is_full

    ! These are the arrays which store the Gibbs data and its derivatives
    REAL(8), DIMENSION(NP, NT, NM) :: G, dG_dP, dG_dT, dG_dM, d2G_dPdT, d2G_dPdM, d2G_dTdM, d3G_dPdTdM, d2G_dT2
    REAL(8), DIMENSION(NT*NM) :: dummy                         ! a dummy array for reading in data



    CONTAINS
        SUBROUTINE convert_PTS(p, T, S, p_MPa, T_K, M)
        !   Subroutine for converting the p, T, S coordinates of MITgcm into the pTM coordinates used by SeaFreeze
            IMPLICIT NONE
            REAL(8), INTENT(IN) :: p, T, S
            REAL(8), INTENT(OUT) :: p_MPa, T_K, M

            p_MPa = p / 1.0E6
            T_K = T + 273.15
            M = S / mm_salt

            RETURN
        END SUBROUTINE convert_PTS

        SUBROUTINE get_spline_indices(p, T, M, i, j, k)
        !   Subroutine which finds the indices used to calculate the interpolating spline at
        !   the specified (p, T, M) values.
            IMPLICIT NONE
            REAL(8), INTENT(IN) :: p, T, M
            INTEGER, INTENT(OUT) :: i, j, k

            i = AINT((p - P0)/dP, KIND(i)) + 1
            j = AINT((T - T0)/dT, KIND(j)) + 1
            k = AINT((M - M0)/dM, KIND(k)) + 1

            PRINT *, i, j, k

            RETURN
        END SUBROUTINE get_spline_indices

        SUBROUTINE load_gibbs_data()
        !   Loads the Gibbs data
            IMPLICIT NONE
            INTEGER :: i_p, i_t, i_m, rc

            OPEN (ACTION='read', FILE='csv/G.csv', IOSTAT=rc, UNIT=20)
            DO i_p = 1, NP
                read (20, *, IOSTAT=rc) dummy(:)        ! read data into dummy array
                DO i_t = 1, NT
                    G(i_p, i_t, :) = dummy( (i_t - 1)*NM + 1: i_t*NM )     ! transfer from dummy array to formal storage
                END DO
            END DO
            CLOSE (20)

            ! dG_dP
            OPEN (ACTION='read', FILE='csv/dGdP.csv', IOSTAT=rc, UNIT=21)
            DO i_p = 1, NP
                read (21, *, IOSTAT=rc) dummy(:)
                DO i_t = 1, NT
                    dG_dP(i_p, i_t, :) = dummy( (i_t - 1)*NM + 1: i_t*NM )
                END DO
            END DO
            CLOSE (21)

            ! dG_dT
            OPEN (ACTION='read', FILE='csv/dGdT.csv', IOSTAT=rc, UNIT=22)
            DO i_p = 1, NP
                read (22, *, IOSTAT=rc) dummy(:)
                DO i_t = 1, NT
                    dG_dT(i_p, i_t, :) = dummy( (i_t - 1)*NM + 1: i_t*NM )
                END DO
            END DO
            CLOSE (22)

            ! dG_dM
            OPEN (ACTION='read', FILE='csv/dGdM.csv', IOSTAT=rc, UNIT=23)
            DO i_p = 1, NP
                read (23, *, IOSTAT=rc) dummy(:)
                DO i_t = 1, NT
                    dG_dM(i_p, i_t, :) = dummy( (i_t - 1)*NM + 1: i_t*NM )
                END DO
            END DO
            CLOSE (23)

            ! now the mixed partials...
            OPEN (ACTION='read', FILE='csv/d2GdPT.csv', IOSTAT=rc, UNIT=24)
            DO i_p = 1, NP
                read (24, *, IOSTAT=rc) dummy(:)
                DO i_t = 1, NT
                    d2G_dPdT(i_p, i_t, :) = dummy( (i_t - 1)*NM + 1: i_t*NM )
                END DO
            END DO
            CLOSE (24)

            OPEN (ACTION='read', FILE='csv/d2GdPM.csv', IOSTAT=rc, UNIT=25)
            DO i_p = 1, NP
                read (25, *, IOSTAT=rc) dummy(:)
                DO i_t = 1, NT
                    d2G_dPdM(i_p, i_t, :) = dummy( (i_t - 1)*NM + 1: i_t*NM )
                END DO
            END DO
            CLOSE (25)

            OPEN (ACTION='read', FILE='csv/d2GdTM.csv', IOSTAT=rc, UNIT=26)
            DO i_p = 1, NP
                read (26, *, IOSTAT=rc) dummy(:)
                DO i_t = 1, NT
                    d2G_dTdM(i_p, i_t, :) = dummy( (i_t - 1)*NM + 1: i_t*NM )
                END DO
            END DO
            CLOSE (26)

            ! the tertiary mixed partial
            OPEN (ACTION='read', FILE='csv/d3GdPTM.csv', IOSTAT=rc, UNIT=27)
            DO i_p = 1, NP
                read (27, *, IOSTAT=rc) dummy(:)
                DO i_t = 1, NT
                    d3G_dPdTdM(i_p, i_t, :) = dummy( (i_t - 1)*NM + 1: i_t*NM )
                END DO
            END DO
            CLOSE (27)

            ! we also want to load d2G_dT2 for calculating theta/adiabatic temp. gradient
            OPEN (ACTION='read', FILE='csv/d2GdT2.csv', IOSTAT=rc, UNIT=28)
            DO i_p = 1, NP
                read (28, *, IOSTAT=rc) dummy(:)
                DO i_t = 1, NT
                    d2G_dT2(i_p, i_t, :) = dummy( (i_t - 1)*NM + 1: i_t*NM )
                END DO
            END DO
            CLOSE (28)
        END SUBROUTINE load_gibbs_data

        SUBROUTINE calc_spline_coeffs(p, T, M, a)
        ! Subroutine for calculating the 64 spline coefficients for the tricubic interpolator
            IMPLICIT NONE
            ! input/output variables 
            REAL(8), INTENT(IN) :: p, T, M
            REAL(8), DIMENSION(4, 4, 4), INTENT(OUT) :: a
            ! internal variables
            REAL(8) :: pp, TT, MM
            REAL(8), DIMENSION(4*4*4) :: a_1D, b_1D
            INTEGER, DIMENSION(64,64) :: B_inv
            INTEGER :: i, j, k, iter, iter_1, iter_2, rc, ii, jj, kk

            ! Load the Gibbs data
            CALL load_gibbs_data()

            ! Convert to proper units/dimensions and get indices
            CALL get_spline_indices(p, T, M, i, j, k)

            ! Fill out the vector b_1D with the derivative values
            b_1D(1:8) = (/ G(i, j, k), G(i+1, j, k), G(i, j+1, k), G(i+1, j+1, k), & 
                & G(i, j, k+1), G(i+1, j, k+1), G(i, j+1, k+1), G(i+1, j+1, k+1) /)
            b_1D(9:16) = (/ dG_dP(i, j, k), dG_dP(i+1, j, k), dG_dP(i, j+1, k), dG_dP(i+1, j+1, k), & 
                & dG_dP(i, j, k+1), dG_dP(i+1, j, k+1), dG_dP(i, j+1, k+1), dG_dP(i+1, j+1, k+1) /)
            b_1D(17:24) = (/ dG_dT(i, j, k), dG_dT(i+1, j, k), dG_dT(i, j+1, k), dG_dT(i+1, j+1, k), & 
                & dG_dT(i, j, k+1), dG_dT(i+1, j, k+1), dG_dT(i, j+1, k+1), dG_dT(i+1, j+1, k+1) /)
            b_1D(25:32) = (/ dG_dM(i, j, k), dG_dM(i+1, j, k), dG_dM(i, j+1, k), dG_dM(i+1, j+1, k), & 
                & dG_dM(i, j, k+1), dG_dM(i+1, j, k+1), dG_dM(i, j+1, k+1), dG_dM(i+1, j+1, k+1) /)
            b_1D(33:40) = (/ d2G_dPdT(i, j, k), d2G_dPdT(i+1, j, k), d2G_dPdT(i, j+1, k), d2G_dPdT(i+1, j+1, k), & 
                & d2G_dPdT(i, j, k+1), d2G_dPdT(i+1, j, k+1), d2G_dPdT(i, j+1, k+1), d2G_dPdT(i+1, j+1, k+1) /)
            b_1D(41:48) = (/ d2G_dPdM(i, j, k), d2G_dPdM(i+1, j, k), d2G_dPdM(i, j+1, k), d2G_dPdM(i+1, j+1, k), & 
                & d2G_dPdM(i, j, k+1), d2G_dPdM(i+1, j, k+1), d2G_dPdM(i, j+1, k+1), d2G_dPdM(i+1, j+1, k+1) /)
            b_1D(49:56) = (/ d2G_dTdM(i, j, k), d2G_dTdM(i+1, j, k), d2G_dTdM(i, j+1, k), d2G_dTdM(i+1, j+1, k), & 
                & d2G_dTdM(i, j, k+1), d2G_dTdM(i+1, j, k+1), d2G_dTdM(i, j+1, k+1), d2G_dTdM(i+1, j+1, k+1) /)
            b_1D(57:64) = (/ d3G_dPdTdM(i, j, k), d3G_dPdTdM(i+1, j, k), d3G_dPdTdM(i, j+1, k), d3G_dPdTdM(i+1, j+1, k), & 
                & d3G_dPdTdM(i, j, k+1), d3G_dPdTdM(i+1, j, k+1), d3G_dPdTdM(i, j+1, k+1), d3G_dPdTdM(i+1, j+1, k+1) /)
            
            ! load the B-matrix inverse found at https://github.com/nbigaouette/libtricubic/blob/master/tricubic-1.0/src/libtricubic/coeff.h
            OPEN (ACTION='read', FILE='csv/Binv.csv', IOSTAT=rc, UNIT=100)
            DO iter = 1, 64
                read (100, *, IOSTAT=rc) B_inv(iter,:)
            END DO
            CLOSE (100)

            ! now calculate the a_1D values
            a_1D(:) = 0
            DO iter_1 = 1, 64
                DO iter_2 = 1, 64
                    a_1D(iter_1) = a_1D(iter_1) + B_inv(iter_1, iter_2) * b_1D(iter_2)
                END DO
            END DO

            ! finally, reshape into the 4*4*4 array
            DO ii = 1, 4
                DO jj = 1, 4
                    DO kk = 1, 4
                        a(ii, jj, kk) = a_1D(16*(kk - 1) + 4*(jj - 1) + ii)
                    END DO
                END DO
            END DO
            RETURN 
        END SUBROUTINE calc_spline_coeffs

        SUBROUTINE find_rho_sf(p, T, S, rho)
        !   Calculates density as a function of p, T, S
            IMPLICIT NONE
            REAL(8), INTENT(IN) :: p, T, S
            REAL(8), INTENT(OUT) :: rho
            REAL(8) :: pp, TT, MM, denom
            REAL(8), DIMENSION(4, 4, 4) :: a
            INTEGER :: i, j, k

            CALL convert_PTS(p, T, S, pp, TT, MM)
            CALL calc_spline_coeffs(pp, TT, MM, a)

            denom = 0
            DO i = 0, 3
                DO j = 0, 3
                    DO k = 0, 3
                        denom = denom + a(i+1, j+1, k+1) * i * (pp**(i - 1)) * (TT**j) * (MM**k)
                    END DO
                END DO
            END DO

            rho = 1/(denom * 1.0E-6)        ! multiply by 1E-6 to convert from MPa to Pa
            RETURN

        END SUBROUTINE find_rho_sf

        SUBROUTINE find_alpha_sf(p, T, S, alpha) 
        !   Calculates thermal expansion coefficient as a function of p, T, S
            IMPLICIT NONE
            REAL(8), INTENT(IN) :: p, T, S
            REAL(8), INTENT(OUT) :: alpha
            REAL(8) :: pp, TT, MM, num, denom
            REAL(8), DIMENSION(4, 4, 4) :: a
            INTEGER :: i, j, k

            num = 0
            denom = 0
            DO i = 0, 3
                DO j = 0, 3
                    DO k = 0, 3
                        num = num + a(i+1, j+1, k+1) * i * j * (pp**(i - 1)) * (TT**(j - 1)) * (MM**k)
                        denom = denom + a(i+1, j+1, k+1) * i * (pp**(i - 1)) * (TT**j) * (MM**k)
                    END DO
                END DO
            END DO

        alpha = num/denom
        RETURN 
        END SUBROUTINE find_alpha_sf

END MODULE sf_salt_thermodynamics

! PROGRAM test_module
!     USE sf_salt_thermodynamics
!     IMPLICIT NONE
!     REAL(8) :: p, T, S, rho, alpha
! 
!     p = 1.0E7           ! Pa
!     T = -3.             ! oC
!     S = 15.             ! g/kg
! 
!     CALL find_rho_sf(p, T, S, rho)
!     CALL find_alpha_sf(p, T, S, alpha)
! 
!     PRINT *, rho
!     PRINT *, alpha
! 
! END PROGRAM test_module 

PROGRAM test_module
    USE sf_salt_thermodynamics
    IMPLICIT NONE
    REAL(8) :: p, T, S, pp, TT, MM, gg
    REAL(8), DIMENSION(4, 4, 4) :: a
    INTEGER :: i, j, k

    p = 1.0E7           ! Pa
    T = -3.             ! oC
    S = 15.             ! g/kg

    CALL convert_PTS(p, T, S, pp, TT, MM)
    CALL calc_spline_coeffs(pp, TT, MM, a)

    gg = 0
    DO i = 0, 3
        DO j = 0, 3
            DO k = 0, 3
                gg = gg + a(i+1, j+1, k+1) * (pp**i) * (TT**j) * (MM**k)
            END DO
        END DO
    END DO

    PRINT *, a(1, 1, 1)
    PRINT *, gg

END PROGRAM test_module