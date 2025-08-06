MODULE sfz_ice_mod
!   This module contains routines for calculating thermodynamic properties of seawater up to 1 GPa

    USE bspline_sub_mod
    IMPLICIT NONE

    REAL, PARAMETER :: mm_salt = 58.44                          ! molecular mass of salt, in g/mol

    ! coefficient array dimensions for each spline
    integer, parameter :: nP_Ih = 10, nT_Ih = 55
    integer, parameter :: nP_II = 10, nT_II = 34
    integer, parameter :: nP_III = 10, nT_III = 34
    integer, parameter :: nP_V = 30, nT_V = 77
    integer, parameter :: nP_VI = 20, nT_VI = 77
!    integer, parameter :: nP_salt = 25, nT_salt = 25, nm_salt = 25
    integer, parameter :: oP = 6, oT = 6
!    integer, parameter :: om = 5

    ! the knot sequences for each spline
    real(wp)    :: P_Ih(nP_Ih + oP), P_II(nP_II + oP), P_III(nP_III + oP), P_V(nP_V + oP), P_VI(nP_VI + oP), &
                   T_Ih(nT_Ih + oT), T_II(nT_II + oT), T_III(nT_III + oT), T_V(nT_V + oT), T_VI(nT_VI + oT)
!    real(wp)    :: P_salt(nP_salt + oP), T_salt(nT_salt + oT), m_salt(nm_salt + om)

    ! and now, the coefficients themselves
    real(wp), dimension(nP_Ih,  nT_Ih)           :: bcoefs_Ih
    real(wp), dimension(nP_II,  nT_II)           :: bcoefs_II
    real(wp), dimension(nP_III, nT_III)          :: bcoefs_III
    real(wp), dimension(nP_V,   nT_V)            :: bcoefs_V
    real(wp), dimension(nP_VI,  nT_VI)           :: bcoefs_VI
!    real(wp), dimension(nP_salt,nT_salt,nm_salt) :: bcoefs_salt

    ! internal variables
   integer :: inbvx = 1, inbvy = 1, iloy = 1
!    integer :: jnbvx = 1, jnbvy = 1, jnbvz = 1, jloy = 1, jloz = 1
!    real(wp), dimension(nT_salt*nm_salt)        :: bcoefs_dummy
!    real(wp), dimension(oT,om)           :: w2_salt  !! work array
!    real(wp), dimension(om)              :: w1_salt  !! work array
!    real(wp), dimension(3_ip*oP)         :: w0_salt  !! work array

!    character(38)       :: spline_dir = '/home/cgliddon/MITgcm/my_exp/sfz_csvs/'


    CONTAINS
        subroutine load_gibbs_data
        
            implicit none
            integer :: ios, i, j, k
            integer :: u = 20
 
            ! load in the knot data -----------------------------------------------------------------------
            OPEN (ACTION='read', FILE='/home/cgliddon/MITgcm/my_exp/sfz_csvs/P_knots_Ih.csv', IOSTAT=ios, UNIT=u)
            read (u, *, IOSTAT=ios) P_Ih(:)
            CLOSE (u)
            u = u + 1

            OPEN (ACTION='read', FILE='/home/cgliddon/MITgcm/my_exp/sfz_csvs/P_knots_II.csv', IOSTAT=ios, UNIT=u)
            read (u, *, IOSTAT=ios) P_II(:)
            CLOSE (u)
            u = u + 1

            OPEN (ACTION='read', FILE='/home/cgliddon/MITgcm/my_exp/sfz_csvs/P_knots_III.csv', IOSTAT=ios, UNIT=u)
            read (u, *, IOSTAT=ios) P_III(:)
            CLOSE (u)
            u = u + 1

            OPEN (ACTION='read', FILE='/home/cgliddon/MITgcm/my_exp/sfz_csvs/P_knots_V.csv', IOSTAT=ios, UNIT=u)
            read (u, *, IOSTAT=ios) P_V(:)
            CLOSE (u)
            u = u + 1

            OPEN (ACTION='read', FILE='/home/cgliddon/MITgcm/my_exp/sfz_csvs/P_knots_VI.csv', IOSTAT=ios, UNIT=u)
            read (u, *, IOSTAT=ios) P_VI(:)
            CLOSE (u)
            u = u + 1

!            OPEN (ACTION='read', FILE='/home/cgliddon/MITgcm/my_exp/sfz_csvs/P_knots_NaCl.csv', IOSTAT=ios, UNIT=u)
!            read (u, *, IOSTAT=ios) P_salt(:)
!            CLOSE (u)
!            u = u + 1

            OPEN (ACTION='read', FILE='/home/cgliddon/MITgcm/my_exp/sfz_csvs/T_knots_Ih.csv', IOSTAT=ios, UNIT=u)
            read (u, *, IOSTAT=ios) T_Ih(:)
            CLOSE (u)
            u = u + 1

            OPEN (ACTION='read', FILE='/home/cgliddon/MITgcm/my_exp/sfz_csvs/T_knots_II.csv', IOSTAT=ios, UNIT=u)
            read (u, *, IOSTAT=ios) T_II(:)
            CLOSE (u)
            u = u + 1

            OPEN (ACTION='read', FILE='/home/cgliddon/MITgcm/my_exp/sfz_csvs/T_knots_III.csv', IOSTAT=ios, UNIT=u)
            read (u, *, IOSTAT=ios) T_III(:)
            CLOSE (u)
            u = u + 1

            OPEN (ACTION='read', FILE='/home/cgliddon/MITgcm/my_exp/sfz_csvs/T_knots_V.csv', IOSTAT=ios, UNIT=u)
            read (u, *, IOSTAT=ios) T_V(:)
            CLOSE (u)
            u = u + 1

            OPEN (ACTION='read', FILE='/home/cgliddon/MITgcm/my_exp/sfz_csvs/T_knots_VI.csv', IOSTAT=ios, UNIT=u)
            read (u, *, IOSTAT=ios) T_VI(:)
            CLOSE (u)
            u = u + 1

!            OPEN (ACTION='read', FILE='/home/cgliddon/MITgcm/my_exp/sfz_csvs/T_knots_NaCl.csv', IOSTAT=ios, UNIT=u)
!            read (u, *, IOSTAT=ios) T_salt(:)
!            CLOSE (u)
!            u = u + 1

!            OPEN (ACTION='read', FILE='/home/cgliddon/MITgcm/my_exp/sfz_csvs/m_knots_NaCl.csv', IOSTAT=ios, UNIT=u)
!            read (u, *, IOSTAT=ios) m_salt(:)
!            CLOSE (u)
!            u = u + 1

            ! load in the spline coefficients -------------------------------------------------------------
            OPEN (ACTION='read', FILE='/home/cgliddon/MITgcm/my_exp/sfz_csvs/bcoefs_Ih.csv', IOSTAT=ios, UNIT=u)
            do i = 1, nP_Ih
                read (u, *, IOSTAT=ios) bcoefs_Ih(i,:)
            enddo 
            CLOSE (u)
            u = u + 1

            OPEN (ACTION='read', FILE='/home/cgliddon/MITgcm/my_exp/sfz_csvs/bcoefs_II.csv', IOSTAT=ios, UNIT=u)
            do i = 1, nP_II
                read (u, *, IOSTAT=ios) bcoefs_II(i,:)
            enddo 
            CLOSE (u)
            u = u + 1

            OPEN (ACTION='read', FILE='/home/cgliddon/MITgcm/my_exp/sfz_csvs/bcoefs_III.csv', IOSTAT=ios, UNIT=u)
            do i = 1, nP_III
                read (u, *, IOSTAT=ios) bcoefs_III(i,:)
            enddo 
            CLOSE (u)
            u = u + 1

            OPEN (ACTION='read', FILE='/home/cgliddon/MITgcm/my_exp/sfz_csvs/bcoefs_V.csv', IOSTAT=ios, UNIT=u)
            do i = 1, nP_V
                read (u, *, IOSTAT=ios) bcoefs_V(i,:)
            enddo 
            CLOSE (u)
            u = u + 1

            OPEN (ACTION='read', FILE='/home/cgliddon/MITgcm/my_exp/sfz_csvs/bcoefs_VI.csv', IOSTAT=ios, UNIT=u)
            do i = 1, nP_VI
                read (u, *, IOSTAT=ios) bcoefs_VI(i,:)
            enddo 
            CLOSE (u)
            u = u + 1

            OPEN (ACTION='read', FILE='/home/cgliddon/MITgcm/my_exp/sfz_csvs/bcoefs_VI.csv', IOSTAT=ios, UNIT=u)
            do i = 1, nP_VI
                read (u, *, IOSTAT=ios) bcoefs_VI(i,:)
            enddo 
            CLOSE (u)
            u = u + 1

!            OPEN (ACTION='read', FILE='/home/cgliddon/MITgcm/my_exp/sfz_csvs/bcoefs_NaCl.csv', IOSTAT=ios, UNIT=u)
!            do i = 1, nP_salt
!              read (u, *, IOSTAT=ios) bcoefs_dummy(:)
!              do j = 1, nm_salt
!                  bcoefs_salt(i,:,j) = bcoefs_dummy(25*j-24:25*j)
!              enddo
!            enddo 
!            CLOSE (u)
        end subroutine load_gibbs_data

        subroutine convert_PTS(p, T, S, p_MPa, T_K, M)
        !   Subroutine for converting the p, T, S coordinates of MITgcm into the pTM coordinates used by SeaFreeze
            IMPLICIT NONE
            REAL(8), INTENT(IN) :: p, T, S
            REAL(8), INTENT(OUT) :: p_MPa, T_K, M
            
            p_MPa = p / 1.0E6
            T_K = T + 273.15
            M = S / mm_salt
            RETURN
        end subroutine convert_PTS

!        SUBROUTINE find_rhosfz_scalar(p, T, S, rho)
!        !   Calculates the density of seawater at the specified pressure, temperature, and salinity
!            IMPLICIT NONE
!            REAL(8), INTENT(IN) :: p, T, S
!            REAL(8), INTENT(OUT) :: rho
!            INTEGER :: iflag
!            REAL(8) :: pp, TT, MM, spline_res
!           
!            CALL convert_PTS(p, T, S, pp, TT, MM)
!            call db3val(pp, TT, MM, 1, 0, 0, P_salt, T_salt, m_salt, nP_salt, nT_salt, nm_salt, oP, oT, om, &
!                bcoefs_salt, spline_res, iflag, jnbvx, jnbvy, jnbvz, jloy, jloz, w2_salt, w1_salt, w0_salt, .false.)            
!            rho = 1 / (spline_res * 1.0E-6)        ! multiply by 1E-6 b/c derivative is in terms of MPa
!            RETURN 
!        END SUBROUTINE find_rhosfz_scalar
!
!        SUBROUTINE find_alphasfz_scalar(p, T, S, alpha) 
!        !   Calculates thermal expansion coefficient as a function of p, T, S
!            IMPLICIT NONE
!            REAL(8), INTENT(IN) :: p, T, S
!            REAL(8), INTENT(OUT) :: alpha
!            INTEGER :: iflag
!            REAL(8) :: pp, TT, MM, spline_num, spline_den
!
!            CALL convert_PTS(p, T, S, pp, TT, MM) 
!            call db3val(pp, TT, MM, 1, 0, 0, P_salt, T_salt, m_salt, nP_salt, nT_salt, nm_salt, oP, oT, om, &
!                bcoefs_salt, spline_den, iflag, jnbvx, jnbvy, jnbvz, jloy, jloz, w2_salt, w1_salt, w0_salt, .false.)
!            call db3val(pp, TT, MM, 1, 1, 0, P_salt, T_salt, m_salt, nP_salt, nT_salt, nm_salt, oP, oT, om, &
!                bcoefs_salt, spline_num, iflag, jnbvx, jnbvy, jnbvz, jloy, jloz, w2_salt, w1_salt, w0_salt, .false.)  
!            alpha = spline_num/spline_den
!            RETURN 
!        END SUBROUTINE find_alphasfz_scalar
!
!        SUBROUTINE find_betasfz_scalar(p, T, S, beta)
!        !   Calculates saline contraction coefficient as a function of p, T, S
!            IMPLICIT NONE
!            REAL(8), INTENT(IN) :: p, T, S
!            REAL(8), INTENT(OUT) :: beta
!            INTEGER :: iflag
!            REAL(8) :: pp, TT, MM, spline_num, spline_den
!
!            CALL convert_PTS(p, T, S, pp, TT, MM) 
!            call db3val(pp, TT, MM, 1, 0, 0, P_salt, T_salt, m_salt, nP_salt, nT_salt, nm_salt, oP, oT, om, &
!                bcoefs_salt, spline_den, iflag, jnbvx, jnbvy, jnbvz, jloy, jloz, w2_salt, w1_salt, w0_salt, .false.)
!            call db3val(pp, TT, MM, 1, 1, 0, P_salt, T_salt, m_salt, nP_salt, nT_salt, nm_salt, oP, oT, om, &
!                bcoefs_salt, spline_num, iflag, jnbvx, jnbvy, jnbvz, jloy, jloz, w2_salt, w1_salt, w0_salt, .false.)  
!            beta = -(spline_num/spline_den)/mm_salt   ! convert back into dimensions of (g/kg)^-1
!            RETURN 
!
!        END SUBROUTINE find_betasfz_scalar

END MODULE sfz_ice_mod
!
!PROGRAM test_module
!    USE sfz_seawater_mod
!    IMPLICIT NONE
!    REAL(8) :: p, T, S, gg, rho, alpha, theta
!    REAL(8), DIMENSION(NP-1, NT-1, NM-1, 4, 4, 4) :: new_spline_coeffs
!
!    p = 147120.6        ! Pa
!    T = 1.0             ! oC
!    S = 20.             ! g/kg
!
!    CALL load_gibbs_data()
!
!    CALL find_gibbs_sf(p, T, S, gg)
!    CALL find_rho_sfz(p, T, S, rho)
!    CALL find_alpha_sf(p, T, S, alpha)
!
!    PRINT *, gg
!    PRINT *, rho
!    PRINT *, alpha
!
!    theta = sw_ptmp(p, T, S)
!    PRINT *, theta
!
!    ! write out coefficients as a binary file
!    OPEN(25, FILE="sf_coeffs.dat", FORM="unformatted")
!    WRITE(25) spline_coeffs
!    CLOSE(25)
!
!    OPEN(25, FILE="/home/cgliddon/MITgcm/sfz_G_coeffs.dat", FORM="unformatted", ACCESS="STREAM", STATUS="OLD")
!    READ(25) new_spline_coeffs
!
!    PRINT *, new_spline_coeffs(13,7,5,2,1,3)
!
!    CALL find_rhosfz_scalar(p, T, S, new_spline_coeffs, rho)
!    CALL find_alphasfz_scalar(p, T, S, new_spline_coeffs, alpha)
!
!    PRINT *, rho
!    PRINT *, alpha
!    
!
!END PROGRAM test_module
