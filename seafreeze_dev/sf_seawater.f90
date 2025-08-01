MODULE sf_seawater
!   This module contains routines for calculating thermodynamic properties of seawater up to 1 GPa
    USE tricubic
    IMPLICIT NONE

    REAL, PARAMETER :: mm_salt = 58.44                         ! molecular mass of salt, in g/mol
    INTEGER, PARAMETER :: NP = 21, NT = 21, NM = 21        ! dimensions of the data arrays
    REAL(8), PARAMETER :: P0 = 0.0, T0 = 250., M0 = 0.0        ! "initial" values of pressure, temp., molality
    REAL(8), PARAMETER :: dP = 50., dT = 5., dM = 0.18          ! grid spacing of pressure, temperature, molality

    ! These are the arrays which store the Gibbs data and its derivatives
    REAL(8), DIMENSION(NP, NT, NM) :: G, dG_dP, dG_dT, dG_dM, d2G_dPdT, d2G_dPdM, d2G_dTdM, d3G_dPdTdM, d2G_dT2
    REAL(8), DIMENSION(NT*NM) :: dummy                         ! a dummy array for reading in data
    REAL(8), DIMENSION(NP-1, NT-1, NM-1, 4, 4, 4) :: spline_coeffs
    LOGICAL :: TD_data_is_loaded = .false.


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

        SUBROUTINE load_gibbs_data()
        !   Loads the Gibbs data
            IMPLICIT NONE
            INTEGER :: i_p, i_t, i_m, rc

            CALL set_grid_params(P0, T0, M0, dP, dT, dM, NP, NT, NM)    ! set the parameters of the grid

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
            ! finally, fill out the spline coefficient array
            CALL calc_all_coeffs(G, dG_dP, dG_dT, dG_dM, d2G_dPdT, d2G_dPdM, d2G_dTdM, d3G_dPdTdM, spline_coeffs)
            TD_data_is_loaded = .true.

        END SUBROUTINE load_gibbs_data

        SUBROUTINE find_gibbs_sf(p, T, S, gg)
        !   Returns the (interpolated) value of the Gibbs function at the specified pressure, temperature, and salinity.
            IMPLICIT NONE
            REAL(8), INTENT(IN) :: p, T, S
            REAL(8), INTENT(OUT) :: gg      ! the Gibbs function value
            REAL(8) :: pp, TT, MM

            IF (TD_data_is_loaded .neqv. .true. ) THEN
                CALL load_gibbs_data()
            END IF
            CALL convert_PTS(p, T, S, pp, TT, MM)
            CALL calc_spline(pp, TT, MM, 0, 0, 0, spline_coeffs, gg)
            RETURN
        END SUBROUTINE find_gibbs_sf

        SUBROUTINE find_rho_sf(p, T, S, rho)
        !   Calculates the density of seawater at the specified pressure, temperature, and salinity
            IMPLICIT NONE
            REAL(8), INTENT(IN) :: p, T, S
            REAL(8), INTENT(OUT) :: rho
            REAL(8) :: pp, TT, MM, dgdp_int

            IF (TD_data_is_loaded .neqv. .true. ) THEN
                CALL load_gibbs_data()
            END IF

            CALL convert_PTS(p, T, S, pp, TT, MM) 
            CALL calc_spline(pp, TT, MM, 1, 0, 0, spline_coeffs, dgdp_int)
            rho = 1 / (dgdp_int * 1.0E-6)        ! multiply by 1E-6 b/c derivative is in terms of MPa
            RETURN 
        END SUBROUTINE find_rho_sf

        SUBROUTINE find_alpha_sf(p, T, S, alpha) 
        !   Calculates thermal expansion coefficient as a function of p, T, S
            IMPLICIT NONE
            REAL(8), INTENT(IN) :: p, T, S
            REAL(8), INTENT(OUT) :: alpha
            REAL(8) :: pp, TT, MM, dgdp_int, d2gdpt_int

            IF (TD_data_is_loaded .neqv. .true. ) THEN
                CALL load_gibbs_data()
            END IF

            CALL convert_PTS(p, T, S, pp, TT, MM) 
            CALL calc_spline(pp, TT, MM, 1, 0, 0, spline_coeffs, dgdp_int)
            CALL calc_spline(pp, TT, MM, 1, 1, 0, spline_coeffs, d2gdpt_int)
            alpha = d2gdpt_int/dgdp_int
            RETURN 
        END SUBROUTINE find_alpha_sf

        FUNCTION ptmp_func(PP, TT, MM, theta)
        !   Function used in the subroutine calc_ptmp, whose root (to be found using the modified N-R method) is
        !   the potential temperature theta.
            IMPLICIT NONE
            REAL(8), INTENT(IN) :: PP, TT, MM, theta
            REAL(8) :: dgdp_T, dgdp_theta
            REAL(8), PARAMETER :: Pref = 0.1        ! reference pressure, in MPa
            REAL(8) :: ptmp_func

            CALL calc_spline(PP, TT, MM, 0, 1, 0, spline_coeffs, dgdp_T)
            CALL calc_spline(Pref, theta, MM, 0, 1, 0, spline_coeffs, dgdp_theta)
            ptmp_func = dgdp_theta - dgdp_T
        END FUNCTION ptmp_func

        FUNCTION ptmp_func_dv(MM, theta)
        !   Function used in the subroutine calc_ptmp, whose value is the partial derivative of ptmp_func with respect
        !   to the potential temperature theta.
            IMPLICIT NONE
            REAL(8), INTENT(IN) :: MM, theta
            REAL(8), PARAMETER :: Pref = 0.1
            REAL(8) :: ptmp_func_dv
            
            CALL calc_spline(Pref, theta, MM, 0, 2, 0, spline_coeffs, ptmp_func_dv)

        END FUNCTION ptmp_func_dv

        FUNCTION sw_ptmp(p, T, S)
        !   Function which calculates the potential temperature of seawater given the pressure, temperature,
        !   and salinity. The function uses an iterative method based on Newton-Raphson root-finding, but with
        !   modifications as described in McDougall & Wotherspoon (2014). 
            IMPLICIT NONE
            REAL(8), INTENT(IN) :: p, T, S
            REAL(8) :: sw_ptmp          ! returned value

            INTEGER :: i
            INTEGER, PARAMETER :: n_iter = 2                            ! number of iterations of the root-finder to perform
            REAL(8) :: PP, TT, MM, theta_kl, theta_ikl, theta_k, theta_ik, theta_kn
            
            IF (TD_data_is_loaded .neqv. .true. ) THEN
                CALL load_gibbs_data()
            END IF

            CALL convert_PTS(p, T, S, PP, TT, MM)
            ! First, we need to make an initial guess at the potential temperature theta. 
            ! Guessing theta = T is usually a pretty reasonable option...
            ! calculate theta_0, theta*_0, theta_1
            theta_kl = TT
            theta_ikl = TT
            theta_k = TT - ( ptmp_func(PP, TT, MM, TT)/ptmp_func_dv(MM, TT) )
            DO i = 1, n_iter
                theta_ik = theta_k - ptmp_func(PP, TT, MM, theta_k) / ptmp_func_dv(MM, 0.5 * (theta_kl + theta_ikl))
                theta_kn = theta_k - ptmp_func(PP, TT, MM, theta_k) / ptmp_func_dv(MM, 0.5 * (theta_k + theta_ik))
                ! update values
                theta_ikl = theta_ik
                theta_kl = theta_k
                theta_k = theta_kn
            END DO

            sw_ptmp = theta_k
        END FUNCTION sw_ptmp

END MODULE sf_seawater

PROGRAM test_module
    USE sf_seawater
    IMPLICIT NONE
    REAL(8) :: p, T, S, gg, rho, alpha, theta

    p = 1.0E7           ! Pa
    T = -3.             ! oC
    S = 15.             ! g/kg

    CALL load_gibbs_data()

    CALL find_gibbs_sf(p, T, S, gg)
    CALL find_rho_sf(p, T, S, rho)
    CALL find_alpha_sf(p, T, S, alpha)

    PRINT *, gg
    PRINT *, rho
    PRINT *, alpha

    theta = sw_ptmp(p, T, S)
    PRINT *, theta


END PROGRAM test_module 
