MODULE tricubic_mod
!   Purpose:
!   Contains subroutines for performing tricubic spline interpolation on gridded data,
!   based on the method of Lekien & Marsden (Int. J. Numer. Meth. Engng., 2005).
    IMPLICIT NONE
    REAL(8) :: x0, y0, z0, dx, dy, dz
    INTEGER :: Nx, Ny, Nz
    ! B_inverse is a matrix used in calculating the interpolator coefficients
    INTEGER, DIMENSION(64, 64) :: B_inverse
    LOGICAL :: b_inverse_loaded = .false.           ! has the matrix been loaded into memory?

    PUBLIC  :: set_grid_params, load_b_inverse, calc_all_coeffs, calc_spline
    PRIVATE :: x0, y0, z0, dx, dy, dz, Nx, Ny, Nz

    CONTAINS
        SUBROUTINE load_b_inverse()
        !   Loads the array B_inverse from the csv file Binv.csv. If loading is successful,
        !   prints out a statement to that effect.
            IMPLICIT NONE
            INTEGER :: i, ios

            OPEN (ACTION='read', FILE='/home/cgliddon/MITgcm/Binv.csv', IOSTAT=ios, UNIT=100, STATUS='OLD')
            DO i = 1, 64
                read (100, *, IOSTAT=ios) B_inverse(i,:)
            END DO
            CLOSE (100)

            b_inverse_loaded = .true.       ! successfully loaded B_inverse into memory
        !    PRINT *, "Successfully loaded B_inverse into memory"
        END SUBROUTINE load_b_inverse

        SUBROUTINE set_grid_params(x0_, y0_, z0_, dx_, dy_, dz_, Nx_, Ny_, Nz_)
        !   Sets the values of grid parameters to be used in the model.
            IMPLICIT NONE
            REAL(8), INTENT(IN) :: x0_, y0_, z0_, dx_, dy_, dz_
            INTEGER, INTENT(IN) :: Nx_, Ny_, Nz_

            x0 = x0_
            y0 = y0_
            z0 = z0_
            dx = dx_
            dy = dy_
            dz = dz_
            Nx = Nx_
            Ny = Ny_
            Nz = Nz_

        END SUBROUTINE set_grid_params

        SUBROUTINE compute_coeffs(x, y, z,                                                        &
                                    f, dfdx, dfdy, dfdz, d2fdxdy, d2fdxdz, d2fdydz, d3fdxdydz,    &
                                    a)
        !   Computes the 4x4x4 coefficient tensor used in the tricubic interpolation. The coefficients
        !   are computed with respect to nondimensionalized, "normalized" coordinates 
            IMPLICIT NONE

            ! input/output variables
            REAL(8), INTENT(IN) :: x, y, z              ! coordinates at which to evaluate
            ! gridded values of f and its derivatives
            REAL(8), DIMENSION(Nx, Ny, Nz), INTENT(IN) :: f, dfdx, dfdy, dfdz, d2fdxdy, d2fdxdz, d2fdydz, d3fdxdydz
            REAL(8), DIMENSION(4, 4, 4), INTENT(OUT) :: a           ! coefficient matrix

            ! internal variables
            INTEGER :: i_x, i_y, i_z, i, j, k
            REAL(8), DIMENSION(64) :: a_vec, b_vec                  ! vector storing 
            REAL(8), DIMENSION(Nx, Ny, Nz) :: dfdx_nd, dfdy_nd, dfdz_nd, d2fdxdy_nd, d2fdxdz_nd, d2fdydz_nd, d3fdxdydz_nd

            ! load B_inverse matrix into memory if we haven't done it yet
            IF (b_inverse_loaded .neqv. .true.) THEN
                CALL load_b_inverse()
            END IF

            ! find the coordinates at which we need to access gridded values
            i_x = AINT((x - x0)/dx, KIND(i_x)) + 1
            i_y = AINT((y - y0)/dy, KIND(i_y)) + 1
            i_z = AINT((z - z0)/dz, KIND(i_z)) + 1

            ! convert the derivative arrays to nondimensional form
            dfdx_nd(:,:,:) = dfdx(:,:,:) * dx
            dfdy_nd(:,:,:) = dfdy(:,:,:) * dy
            dfdz_nd(:,:,:) = dfdz(:,:,:) * dz
            d2fdxdy_nd(:,:,:) = d2fdxdy(:,:,:) * dx * dy
            d2fdxdz_nd(:,:,:) = d2fdxdz(:,:,:) * dx * dz
            d2fdydz_nd(:,:,:) = d2fdydz(:,:,:) * dy * dz
            d3fdxdydz_nd(:,:,:) = d3fdxdydz(:,:,:) * dx * dy * dz
            
            ! load the function and derivative values into b_vec
            b_vec(1:8)   = (/ f(i_x, i_y, i_z),       f(i_x+1, i_y, i_z),     &
                              f(i_x, i_y+1, i_z),     f(i_x+1, i_y+1, i_z),   &
                              f(i_x, i_y, i_z+1),     f(i_x+1, i_y, i_z+1),   &
                              f(i_x, i_y+1, i_z+1),   f(i_x+1, i_y+1, i_z+1) /)
            b_vec(9:16)  = (/ dfdx_nd(i_x, i_y, i_z),       dfdx_nd(i_x+1, i_y, i_z),     &
                              dfdx_nd(i_x, i_y+1, i_z),     dfdx_nd(i_x+1, i_y+1, i_z),   &
                              dfdx_nd(i_x, i_y, i_z+1),     dfdx_nd(i_x+1, i_y, i_z+1),   &
                              dfdx_nd(i_x, i_y+1, i_z+1),   dfdx_nd(i_x+1, i_y+1, i_z+1) /)
            b_vec(17:24) = (/ dfdy_nd(i_x, i_y, i_z),       dfdy_nd(i_x+1, i_y, i_z),     &
                              dfdy_nd(i_x, i_y+1, i_z),     dfdy_nd(i_x+1, i_y+1, i_z),   &
                              dfdy_nd(i_x, i_y, i_z+1),     dfdy_nd(i_x+1, i_y, i_z+1),   &
                              dfdy_nd(i_x, i_y+1, i_z+1),   dfdy_nd(i_x+1, i_y+1, i_z+1) /)
            b_vec(25:32) = (/ dfdz_nd(i_x, i_y, i_z),       dfdz_nd(i_x+1, i_y, i_z),     &
                              dfdz_nd(i_x, i_y+1, i_z),     dfdz_nd(i_x+1, i_y+1, i_z),   &
                              dfdz_nd(i_x, i_y, i_z+1),     dfdz_nd(i_x+1, i_y, i_z+1),   &
                              dfdz_nd(i_x, i_y+1, i_z+1),   dfdz_nd(i_x+1, i_y+1, i_z+1) /)
            b_vec(33:40) = (/ d2fdxdy_nd(i_x, i_y, i_z),       d2fdxdy_nd(i_x+1, i_y, i_z),     &
                              d2fdxdy_nd(i_x, i_y+1, i_z),     d2fdxdy_nd(i_x+1, i_y+1, i_z),   &
                              d2fdxdy_nd(i_x, i_y, i_z+1),     d2fdxdy_nd(i_x+1, i_y, i_z+1),   &
                              d2fdxdy_nd(i_x, i_y+1, i_z+1),   d2fdxdy_nd(i_x+1, i_y+1, i_z+1) /)
            b_vec(41:48) = (/ d2fdxdz_nd(i_x, i_y, i_z),       d2fdxdz_nd(i_x+1, i_y, i_z),     &
                              d2fdxdz_nd(i_x, i_y+1, i_z),     d2fdxdz_nd(i_x+1, i_y+1, i_z),   &
                              d2fdxdz_nd(i_x, i_y, i_z+1),     d2fdxdz_nd(i_x+1, i_y, i_z+1),   &
                              d2fdxdz_nd(i_x, i_y+1, i_z+1),   d2fdxdz_nd(i_x+1, i_y+1, i_z+1) /)
            b_vec(49:56) = (/ d2fdydz_nd(i_x, i_y, i_z),       d2fdydz_nd(i_x+1, i_y, i_z),     &
                              d2fdydz_nd(i_x, i_y+1, i_z),     d2fdydz_nd(i_x+1, i_y+1, i_z),   &
                              d2fdydz_nd(i_x, i_y, i_z+1),     d2fdydz_nd(i_x+1, i_y, i_z+1),   &
                              d2fdydz_nd(i_x, i_y+1, i_z+1),   d2fdydz_nd(i_x+1, i_y+1, i_z+1) /)
            b_vec(57:64) = (/ d3fdxdydz_nd(i_x, i_y, i_z),       d3fdxdydz_nd(i_x+1, i_y, i_z),     &
                              d3fdxdydz_nd(i_x, i_y+1, i_z),     d3fdxdydz_nd(i_x+1, i_y+1, i_z),   &
                              d3fdxdydz_nd(i_x, i_y, i_z+1),     d3fdxdydz_nd(i_x+1, i_y, i_z+1),   &
                              d3fdxdydz_nd(i_x, i_y+1, i_z+1),   d3fdxdydz_nd(i_x+1, i_y+1, i_z+1) /)
            
            ! now calculate coefficients
            a_vec(:) = 0
            DO j = 1, 64
                DO k = 1, 64
                    a_vec(j) = a_vec(j) + B_inverse(j, k) * b_vec(k)
                END DO 
            END DO 

            ! finally, reshape into the 4x4x4 tensor
            DO i = 0, 3
                DO j = 0, 3
                    DO k = 0, 3
                        a(i+1, j+1, k+1) = a_vec(1 + i + 4*j + 16*k)
                    END DO
                END DO
            END DO
            RETURN
        END SUBROUTINE compute_coeffs

!         SUBROUTINE calc_spline(x, y, z, xder, yder, zder, f, dfdx, dfdy, dfdz, d2fdxdy, d2fdxdz, d2fdydz, d3fdxdydz, result)
!         !   Subroutine for calculating the value of the spline or one of its derivatives
!             IMPLICIT NONE
!             REAL(8), INTENT(IN) :: x, y, z
!             REAL(8), DIMENSION(Nx, Ny, Nz), INTENT(IN) :: f, dfdx, dfdy, dfdz, d2fdxdy, d2fdxdz, d2fdydz, d3fdxdydz
!             INTEGER, INTENT(IN) :: xder, yder, zder
!             REAL(8), INTENT(OUT) :: result
! 
!             REAL(8) :: x_nd, y_nd, z_nd
!             REAL(8), DIMENSION(4, 4, 4) :: a
!             INTEGER :: i, j, k
! 
!             ! calculate coefficients
!             CALL compute_coeffs(x, y, z, f, dfdx, dfdy, dfdz, d2fdxdy, d2fdxdz, d2fdydz, d3fdxdydz, a)
! 
!             ! convert x, y, z to nondimensional coordinates
!             x_nd = (x - x0)/dx - AINT((x - x0)/dx) 
!             y_nd = (y - y0)/dy - AINT((y - y0)/dy)
!             z_nd = (z - z0)/dz - AINT((z - z0)/dz) 
! 
!             result = 0
!             IF ((xder == 0) .AND. (yder == 0) .AND. (zder == 0)) THEN
!                 DO i = 0, 3
!                     DO j = 0, 3
!                         DO k = 0, 3
!                             result = result + a(i+1, j+1, k+1) * (x_nd ** i) * (y_nd ** j) * (z_nd ** k)
!                         END DO
!                     END DO
!                 END DO
!             ELSE IF ((xder == 1) .AND. (yder == 0) .AND. (zder == 0)) THEN
!                 DO i = 0, 3
!                     DO j = 0, 3
!                         DO k = 0, 3
!                             result = result + a(i+1, j+1, k+1) * i*(x_nd ** (i - 1)) * (y_nd ** j) * (z_nd ** k)
!                         END DO
!                     END DO
!                 END DO
!                 ! convert result back to dimensionful units
!                 result = result / dx
!             ELSE IF ((xder == 0) .AND. (yder == 1) .AND. (zder == 0)) THEN
!                 DO i = 0, 3
!                     DO j = 0, 3
!                         DO k = 0, 3
!                             result = result + a(i+1, j+1, k+1) * j*(x_nd ** i) * (y_nd ** (j - 1)) * (z_nd ** k)
!                         END DO
!                     END DO
!                 END DO
!                 ! convert result back to dimensionful units
!                 result = result / dy
!             ELSE IF ((xder == 1) .AND. (yder == 1) .AND. (zder == 0)) THEN
!                 DO i = 0, 3
!                     DO j = 0, 3
!                         DO k = 0, 3
!                             result = result + a(i+1, j+1, k+1) * i * j * (x_nd ** (i - 1)) * (y_nd ** (j - 1)) * (z_nd ** k)
!                         END DO
!                     END DO
!                 END DO
!                 ! convert result back to dimensionful units
!                 result = result / (dx * dy)
!             ELSE IF ((xder == 0) .AND. (yder == 2) .AND. (zder == 0)) THEN
!                 DO i = 0, 3
!                     DO j = 0, 3
!                         DO k = 0, 3
!                             result = result + a(i+1, j+1, k+1) * j * (j - 1) * (x_nd ** i) * (y_nd ** (j - 2)) * (z_nd ** k)
!                         END DO
!                     END DO
!                 END DO
!                 ! convert result back to dimensionful units
!                 result = result / (dy ** 2)
!             END IF
! 
!             RETURN
!         END SUBROUTINE calc_spline

        SUBROUTINE calc_all_coeffs(f, dfdx, dfdy, dfdz, d2fdxdy, d2fdxdz, d2fdydz, d3fdxdydz, outdat)
            IMPLICIT NONE
            REAL(8), DIMENSION(Nx, Ny, Nz), INTENT(IN) :: f, dfdx, dfdy, dfdz, d2fdxdy, d2fdxdz, d2fdydz, d3fdxdydz
            REAL(8), DIMENSION(Nx - 1, Ny - 1, Nz - 1, 4, 4, 4), INTENT(OUT) :: outdat     
            ! note that the coefficient matrix has a dimension one less than Nx/Ny/Nz 
            ! b/c there's one set of coefficients for each "cube" not each grid point
 
            INTEGER :: i_x, i_y, i_z, i, j, k
            REAL(8), DIMENSION(Nx, Ny, Nz) :: dfdx_nd, dfdy_nd, dfdz_nd, d2fdxdy_nd, d2fdxdz_nd, d2fdydz_nd, d3fdxdydz_nd
            REAL(8), DIMENSION(64) :: a_vec, b_vec                  ! vector storing
            REAL(8), DIMENSION(4, 4, 4) :: a 

            IF (b_inverse_loaded .neqv. .true.) THEN
                CALL load_b_inverse()
            END IF

            ! convert the derivative arrays to nondimensional form
            dfdx_nd(:,:,:) = dfdx(:,:,:) * dx
            dfdy_nd(:,:,:) = dfdy(:,:,:) * dy
            dfdz_nd(:,:,:) = dfdz(:,:,:) * dz
            d2fdxdy_nd(:,:,:) = d2fdxdy(:,:,:) * dx * dy
            d2fdxdz_nd(:,:,:) = d2fdxdz(:,:,:) * dx * dz
            d2fdydz_nd(:,:,:) = d2fdydz(:,:,:) * dy * dz
            d3fdxdydz_nd(:,:,:) = d3fdxdydz(:,:,:) * dx * dy * dz

            DO i_x = 1, Nx - 1
                DO i_y = 1, Ny - 1
                    DO i_z = 1, Nz - 1
                        b_vec(1:8)   = (/ f(i_x, i_y, i_z),     f(i_x+1, i_y, i_z),     &
                                        f(i_x, i_y+1, i_z),     f(i_x+1, i_y+1, i_z),   &
                                        f(i_x, i_y, i_z+1),     f(i_x+1, i_y, i_z+1),   &
                                        f(i_x, i_y+1, i_z+1),   f(i_x+1, i_y+1, i_z+1) /)
                        b_vec(9:16)  = (/ dfdx_nd(i_x, i_y, i_z),     dfdx_nd(i_x+1, i_y, i_z),     &
                                        dfdx_nd(i_x, i_y+1, i_z),     dfdx_nd(i_x+1, i_y+1, i_z),   &
                                        dfdx_nd(i_x, i_y, i_z+1),     dfdx_nd(i_x+1, i_y, i_z+1),   &
                                        dfdx_nd(i_x, i_y+1, i_z+1),   dfdx_nd(i_x+1, i_y+1, i_z+1) /)
                        b_vec(17:24) = (/ dfdy_nd(i_x, i_y, i_z),     dfdy_nd(i_x+1, i_y, i_z),     &
                                        dfdy_nd(i_x, i_y+1, i_z),     dfdy_nd(i_x+1, i_y+1, i_z),   &
                                        dfdy_nd(i_x, i_y, i_z+1),     dfdy_nd(i_x+1, i_y, i_z+1),   &
                                        dfdy_nd(i_x, i_y+1, i_z+1),   dfdy_nd(i_x+1, i_y+1, i_z+1) /)
                        b_vec(25:32) = (/ dfdz_nd(i_x, i_y, i_z),     dfdz_nd(i_x+1, i_y, i_z),     &
                                        dfdz_nd(i_x, i_y+1, i_z),     dfdz_nd(i_x+1, i_y+1, i_z),   &
                                        dfdz_nd(i_x, i_y, i_z+1),     dfdz_nd(i_x+1, i_y, i_z+1),   &
                                        dfdz_nd(i_x, i_y+1, i_z+1),   dfdz_nd(i_x+1, i_y+1, i_z+1) /)
                        b_vec(33:40) = (/ d2fdxdy_nd(i_x, i_y, i_z),     d2fdxdy_nd(i_x+1, i_y, i_z),     &
                                        d2fdxdy_nd(i_x, i_y+1, i_z),     d2fdxdy_nd(i_x+1, i_y+1, i_z),   &
                                        d2fdxdy_nd(i_x, i_y, i_z+1),     d2fdxdy_nd(i_x+1, i_y, i_z+1),   &
                                        d2fdxdy_nd(i_x, i_y+1, i_z+1),   d2fdxdy_nd(i_x+1, i_y+1, i_z+1) /)
                        b_vec(41:48) = (/ d2fdxdz_nd(i_x, i_y, i_z),     d2fdxdz_nd(i_x+1, i_y, i_z),     &
                                        d2fdxdz_nd(i_x, i_y+1, i_z),     d2fdxdz_nd(i_x+1, i_y+1, i_z),   &
                                        d2fdxdz_nd(i_x, i_y, i_z+1),     d2fdxdz_nd(i_x+1, i_y, i_z+1),   &
                                        d2fdxdz_nd(i_x, i_y+1, i_z+1),   d2fdxdz_nd(i_x+1, i_y+1, i_z+1) /)
                        b_vec(49:56) = (/ d2fdydz_nd(i_x, i_y, i_z),     d2fdydz_nd(i_x+1, i_y, i_z),     &
                                        d2fdydz_nd(i_x, i_y+1, i_z),     d2fdydz_nd(i_x+1, i_y+1, i_z),   &
                                        d2fdydz_nd(i_x, i_y, i_z+1),     d2fdydz_nd(i_x+1, i_y, i_z+1),   &
                                        d2fdydz_nd(i_x, i_y+1, i_z+1),   d2fdydz_nd(i_x+1, i_y+1, i_z+1) /)
                        b_vec(57:64) = (/ d3fdxdydz_nd(i_x, i_y, i_z),     d3fdxdydz_nd(i_x+1, i_y, i_z),     &
                                        d3fdxdydz_nd(i_x, i_y+1, i_z),     d3fdxdydz_nd(i_x+1, i_y+1, i_z),   &
                                        d3fdxdydz_nd(i_x, i_y, i_z+1),     d3fdxdydz_nd(i_x+1, i_y, i_z+1),   &
                                        d3fdxdydz_nd(i_x, i_y+1, i_z+1),   d3fdxdydz_nd(i_x+1, i_y+1, i_z+1) /)
                        ! now calculate coefficients
                        a_vec(:) = 0
                        DO j = 1, 64
                            DO k = 1, 64
                                a_vec(j) = a_vec(j) + B_inverse(j, k) * b_vec(k)
                            END DO 
                        END DO 

                        ! finally, reshape into the 4x4x4 tensor
                        DO i = 0, 3
                            DO j = 0, 3
                                DO k = 0, 3
                                    a(i+1, j+1, k+1) = a_vec(1 + i + 4*j + 16*k)
                                END DO
                            END DO
                        END DO
                        outdat(i_x, i_y, i_z, :, :, :) = a(:, :, :)

                    END DO
                END DO
            END DO
        !    PRINT *, "Gibbs function spline coefficients successfully calculated"
            RETURN

        END SUBROUTINE calc_all_coeffs 

        SUBROUTINE calc_spline(x, y, z, xder, yder, zder, coeffs, result)
        !   Subroutine for calculating the value of the spline or one of its derivatives
            IMPLICIT NONE
            REAL(8), INTENT(IN) :: x, y, z
            INTEGER, INTENT(IN) :: xder, yder, zder
            REAL(8), DIMENSION(Nx-1, Ny-1, Nz-1, 4, 4, 4), INTENT(IN) :: coeffs
            REAL(8), INTENT(OUT) :: result

            REAL(8) :: x_nd, y_nd, z_nd
            REAL(8) :: a000, a001, a002, a003, a010, a011, a012, a013, a020, a021, a022, a023, a030, a031, a032, a033, &
                       a100, a101, a102, a103, a110, a111, a112, a113, a120, a121, a122, a123, a130, a131, a132, a133, &
                       a200, a201, a202, a203, a210, a211, a212, a213, a220, a221, a222, a223, a230, a231, a232, a233, &
                       a300, a301, a302, a303, a310, a311, a312, a313, a320, a321, a322, a323, a330, a331, a332, a333
            INTEGER :: i, j, k, i_x, i_y, i_z
            
            ! load B_inverse matrix into memory if we haven't done it yet
            IF (b_inverse_loaded .neqv. .true.) THEN
                CALL load_b_inverse()
            END IF

!           PRINT *, 'Test coefficient : ', coeffs(13, 7, 5, 2, 1, 3)

            ! find the coordinates at which we need to access gridded values
            i_x = AINT((x - x0)/dx, KIND(i_x)) + 1
            i_y = AINT((y - y0)/dy, KIND(i_y)) + 1
            i_z = AINT((z - z0)/dz, KIND(i_z)) + 1

            ! and get the coefficients
            a000 = coeffs(i_x,i_y,i_z,0,0,0)
            a001 = coeffs(i_x,i_y,i_z,0,0,1)
            a002 = coeffs(i_x,i_y,i_z,0,0,2)
            a003 = coeffs(i_x,i_y,i_z,0,0,3)
            a010 = coeffs(i_x,i_y,i_z,0,1,0)
            a011 = coeffs(i_x,i_y,i_z,0,1,1)
            a012 = coeffs(i_x,i_y,i_z,0,1,2)
            a013 = coeffs(i_x,i_y,i_z,0,1,3)
            a020 = coeffs(i_x,i_y,i_z,0,2,0)
            a021 = coeffs(i_x,i_y,i_z,0,2,1)
            a022 = coeffs(i_x,i_y,i_z,0,2,2)
            a023 = coeffs(i_x,i_y,i_z,0,2,3)
            a030 = coeffs(i_x,i_y,i_z,0,3,0)
            a031 = coeffs(i_x,i_y,i_z,0,3,1)
            a032 = coeffs(i_x,i_y,i_z,0,3,2)
            a033 = coeffs(i_x,i_y,i_z,0,3,3)
            a100 = coeffs(i_x,i_y,i_z,1,0,0)
            a101 = coeffs(i_x,i_y,i_z,1,0,1)
            a102 = coeffs(i_x,i_y,i_z,1,0,2)
            a103 = coeffs(i_x,i_y,i_z,1,0,3)
            a110 = coeffs(i_x,i_y,i_z,1,1,0)
            a111 = coeffs(i_x,i_y,i_z,1,1,1)
            a112 = coeffs(i_x,i_y,i_z,1,1,2)
            a113 = coeffs(i_x,i_y,i_z,1,1,3)
            a120 = coeffs(i_x,i_y,i_z,1,2,0)
            a121 = coeffs(i_x,i_y,i_z,1,2,1)
            a122 = coeffs(i_x,i_y,i_z,1,2,2)
            a123 = coeffs(i_x,i_y,i_z,1,2,3)
            a130 = coeffs(i_x,i_y,i_z,1,3,0)
            a131 = coeffs(i_x,i_y,i_z,1,3,1)
            a132 = coeffs(i_x,i_y,i_z,1,3,2)
            a133 = coeffs(i_x,i_y,i_z,1,3,3)
            a200 = coeffs(i_x,i_y,i_z,2,0,0)
            a201 = coeffs(i_x,i_y,i_z,2,0,1)
            a202 = coeffs(i_x,i_y,i_z,2,0,2)
            a203 = coeffs(i_x,i_y,i_z,2,0,3)
            a210 = coeffs(i_x,i_y,i_z,2,1,0)
            a211 = coeffs(i_x,i_y,i_z,2,1,1)
            a212 = coeffs(i_x,i_y,i_z,2,1,2)
            a213 = coeffs(i_x,i_y,i_z,2,1,3)
            a220 = coeffs(i_x,i_y,i_z,2,2,0)
            a221 = coeffs(i_x,i_y,i_z,2,2,1)
            a222 = coeffs(i_x,i_y,i_z,2,2,2)
            a223 = coeffs(i_x,i_y,i_z,2,2,3)
            a230 = coeffs(i_x,i_y,i_z,2,3,0)
            a231 = coeffs(i_x,i_y,i_z,2,3,1)
            a232 = coeffs(i_x,i_y,i_z,2,3,2)
            a233 = coeffs(i_x,i_y,i_z,2,3,3)
            a300 = coeffs(i_x,i_y,i_z,3,0,0)
            a301 = coeffs(i_x,i_y,i_z,3,0,1)
            a302 = coeffs(i_x,i_y,i_z,3,0,2)
            a303 = coeffs(i_x,i_y,i_z,3,0,3)
            a310 = coeffs(i_x,i_y,i_z,3,1,0)
            a311 = coeffs(i_x,i_y,i_z,3,1,1)
            a312 = coeffs(i_x,i_y,i_z,3,1,2)
            a313 = coeffs(i_x,i_y,i_z,3,1,3)
            a320 = coeffs(i_x,i_y,i_z,3,2,0)
            a321 = coeffs(i_x,i_y,i_z,3,2,1)
            a322 = coeffs(i_x,i_y,i_z,3,2,2)
            a323 = coeffs(i_x,i_y,i_z,3,2,3)
            a330 = coeffs(i_x,i_y,i_z,3,3,0)
            a331 = coeffs(i_x,i_y,i_z,3,3,1)
            a332 = coeffs(i_x,i_y,i_z,3,3,2)
            a333 = coeffs(i_x,i_y,i_z,3,3,3)

            ! convert x, y, z to nondimensional coordinates
            x_nd = (x - x0)/dx - AINT((x - x0)/dx) 
            y_nd = (y - y0)/dy - AINT((y - y0)/dy)
            z_nd = (z - z0)/dz - AINT((z - z0)/dz) 

            result = 0
            IF ((xder == 0) .AND. (yder == 0) .AND. (zder == 0)) THEN
                result =   a000 + z_nd*(a001 + z_nd*(a002 + z_nd*a003))                                         &
                         + y_nd*(a010 + y_nd*(a020 + y_nd*a030) + z_nd*( a011 + z_nd*(a012 + z_nd*a013))        &
                                 + y_nd*(a021 + y_nd*a031 + z_nd*(a022 + z_nd*a023 + y_nd*(a032 + z_nd*a033)))) &
                         + x_nd*(a100 + x_nd*(a200 + x_nd*a300) + z_nd*(a101 + z_nd*(a102 + z_nd*a103) &
                                                                        + x_nd*( a201 + x_nd*a301 + z_nd*(a202 + z_nd*a203 + x_nd*(a302 + z_nd*a303)))) &
                                 + y_nd*(a110 + y_nd*(a120 + y_nd*a130) + x_nd*( a210 + x_nd*a310 + y_nd*(a220 + y_nd*a230 + x_nd*(a320 + y_nd*a330)))  &
                                        + z_nd*(a111 + z_nd*(a112 + z_nd*a113) + y_nd*(a121 + y_nd*a131 + z_nd*(a112 + z_nd*a123 + y_nd*(a132 + z_nd*a133))) &
                                                + a211 + x_nd*a311 + z_nd*(a121 + z_nd*a213 + x_nd*(a312 + z_nd*a313))           &
                                                                           + y_nd*(a221 + y_nd*a231 + x_nd*(a321 + y_nd*a331)    &
                                                                              + z_nd*(a222 + z_nd*a223 + y_nd*(a232 + z_nd*a233) &
                                                                               +x_nd*(a322 + z_nd*a323 + y_nd*(a332 + z_nd*a333)))))))
            ELSE IF ((xder == 1) .AND. (yder == 0) .AND. (zder == 0)) THEN
                result = a100 + 2*x_nd*(a200 + 3*x_nd*a300)                                                                                                 &
                         + z_nd*(a101 + x_nd*(2*a201 + 3*a301*x_nd) + z_nd*(a102 + a103*z_nd + x_nd*(2*a202 + 3*a302*x_nd + z_nd*(2*a203 + 3*a303*x_nd))))  &
                         + y_nd*(a110 + x_nd*(2*a210 + 3*a310*x_nd) + y_nd*(a120 + a130*y_nd + x_nd*(2*a220 + 3*a320*x_nd + y_nd*(2*a230 + 3*a330*x_nd)))   &
                          + z_nd*(a111 + x_nd*(2*a211 + 3*a311*x_nd) + z_nd*(a112 + a113*z_nd + x_nd*(2*a212 + 3*a312*x_nd + z_nd*(2*a213 + 3*a313*x_nd)))  &
                           + y_nd*(a121 + y_nd*a131 + x_nd*(2*a221 + 3*a321*x_nd + y_nd*(2*a231 + 3*a331*x_nd))                             &
                            + z_nd*(a122 + y_nd*a132 + z_nd*(a123 + y_nd*a133) + x_nd*(2*a222 + 3*a322*x_nd + y_nd*(2*a232 + 3*a332*x_nd)   &
                              + z_nd*(2*a223 + 2*y_nd*a233 + 3*x_nd*(a323 + a333*y_nd)))))))
                ! convert result back to dimensionful units
                result = result / dx
            ELSE IF ((xder == 0) .AND. (yder == 1) .AND. (zder == 0)) THEN
                result = a010 + 2*y_nd*(a020 + 3*y_nd*a030)                                                                                                 &
                         + z_nd*(a011 + y_nd*(2*a021 + 3*a031*y_nd) + z_nd*(a012 + a013*z_nd + y_nd*(2*a022 + 3*a032*y_nd + z_nd*(2*a023 + 3*a033*y_nd))))  &
                         + x_nd*(a110 + y_nd*(2*a120 + 3*a130*y_nd) + x_nd*(a210 + a310*x_nd + y_nd*(2*a220 + 3*a230*y_nd + x_nd*(2*a320 + 3*a330*y_nd)))   &
                          + z_nd*(a111 + y_nd*(2*a121 + 3*a131*y_nd) + z_nd*(a112 + a113*z_nd + y_nd*(2*a122 + 3*a132*y_nd + z_nd*(2*a123 + 3*a133*y_nd)))  &
                           + x_nd*(a211 + x_nd*a311 + y_nd*(2*a221 + 3*a231*y_nd + x_nd*(2*a321 + 3*a331*y_nd))                             &
                            + z_nd*(a212 + x_nd*a312 + z_nd*(a213 + x_nd*a313) + y_nd*(2*a222 + 3*a232*y_nd + x_nd*(2*a322 + 3*a332*y_nd)   &
                              + z_nd*(2*a223 + 2*x_nd*a323 + 3*y_nd*(a233 + a333*x_nd)))))))
                ! convert result back to dimensionful units
                result = result / dy
            ELSE IF ((xder == 1) .AND. (yder == 1) .AND. (zder == 0)) THEN
                result = a110 + 2*x_nd*a210 + 3*a310*x_nd**2 + 2*y_nd*a120 + 3*a130*y_nd**2 + 4*a220*x_nd*y_nd + 6*a320*y_nd*x_nd**2 + 6*a230*x_nd*y_nd**2  &
                            + 9*a330*(x_nd**2)*(y_nd**2)    &
                         + a111*z_nd + 2*x_nd*a211*z_nd + 3*a311*z_nd*x_nd**2 + 2*y_nd*a121*z_nd + 3*a131*z_nd*y_nd**2 + 4*a221*x_nd*y_nd*z_nd &
                         + 6*a321*z_nd*y_nd*x_nd**2 + 6*a231*z_nd*x_nd*y_nd**2 + 9*a331*z_nd*(x_nd**2)*(y_nd**2)    & 
                         + a112*z_nd**2 + 2*x_nd*a212*z_nd**2 + 3*a312*(z_nd**2)*x_nd**2 + 2*y_nd*a122*z_nd**2 + 3*a132*(z_nd**2)*y_nd**2 + 4*a222*x_nd*y_nd*z_nd**2 &
                         + 6*a322*(z_nd**2)*y_nd*x_nd**2 + 6*a232*(z_nd**2)*x_nd*y_nd**2 + 9*a332*(z_nd**2)*(x_nd**2)*(y_nd**2)    & 
                         + a113*z_nd**3 + 2*x_nd*a213*z_nd**3 + 3*a313*(z_nd**3)*x_nd**2 + 2*y_nd*a123*z_nd**3 + 3*a133*(z_nd**3)*y_nd**2 + 4*a223*x_nd*y_nd*z_nd**3 &
                         + 6*a323*(z_nd**3)*y_nd*x_nd**2 + 6*a233*(z_nd**3)*x_nd*y_nd**2 + 9*a333*(z_nd**3)*(x_nd**2)*(y_nd**2)
                ! convert result back to dimensionful units
                result = result / (dx * dy)
            ELSE IF ((xder == 1) .AND. (yder == 0) .AND. (zder == 1)) THEN
                result = a101 + 2*x_nd*a201 + 3*a301*x_nd**2 + 2*z_nd*a102 + 3*a103*z_nd**2 + 4*a202*x_nd*z_nd + 6*a302*z_nd*x_nd**2 + 6*a203*x_nd*z_nd**2  &
                            + 9*a303*(x_nd**2)*(z_nd**2)    &
                         + a111*y_nd + 2*x_nd*a211*y_nd + 3*a311*y_nd*x_nd**2 + 2*z_nd*a112*y_nd + 3*a113*y_nd*z_nd**2 + 4*a212*x_nd*z_nd*y_nd &
                         + 6*a312*y_nd*z_nd*x_nd**2 + 6*a213*y_nd*x_nd*z_nd**2 + 9*a313*y_nd*(x_nd**2)*(z_nd**2)    & 
                         + a121*y_nd**2 + 2*x_nd*a221*y_nd**2 + 3*a321*(y_nd**2)*x_nd**2 + 2*z_nd*a122*y_nd**2 + 3*a123*(y_nd**2)*z_nd**2 + 4*a222*x_nd*z_nd*y_nd**2 &
                         + 6*a322*(y_nd**2)*z_nd*x_nd**2 + 6*a223*(y_nd**2)*x_nd*z_nd**2 + 9*a323*(y_nd**2)*(x_nd**2)*(z_nd**2)    & 
                         + a131*y_nd**3 + 2*x_nd*a231*y_nd**3 + 3*a331*(y_nd**3)*x_nd**2 + 2*z_nd*a132*y_nd**3 + 3*a133*(y_nd**3)*z_nd**2 + 4*a232*x_nd*z_nd*y_nd**3 &
                         + 6*a332*(y_nd**3)*z_nd*x_nd**2 + 6*a233*(y_nd**3)*x_nd*z_nd**2 + 9*a333*(y_nd**3)*(x_nd**2)*(z_nd**2)
                ! convert result back to dimensionful units
                result = result / (dx * dz)
            ELSE IF ((xder == 0) .AND. (yder == 2) .AND. (zder == 0)) THEN
                ! TO DO: write up code for this
                result = 0
            END IF

            RETURN
        END SUBROUTINE calc_spline
    
END MODULE tricubic_mod
