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

    CONTAINS
        SUBROUTINE load_b_inverse()
        !   Loads the array B_inverse from the csv file Binv.csv. If loading is successful,
        !   prints out a statement to that effect.
            IMPLICIT NONE
            INTEGER :: i, ios

            OPEN (ACTION='read', FILE='Binv.csv', IOSTAT=ios, UNIT=100)
            DO i = 1, 64
                read (100, *, IOSTAT=ios) B_inverse(i,:)
            END DO
            CLOSE (100)

            b_inverse_loaded = .true.       ! successfully loaded B_inverse into memory
            PRINT *, "Successfully loaded B_inverse into memory"
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
            PRINT *, "Coefficients successfully calculated"
            PRINT *, SHAPE(outdat)
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
            REAL(8), DIMENSION(4, 4, 4) :: a
            INTEGER :: i, j, k, i_x, i_y, i_z

            ! load B_inverse matrix into memory if we haven't done it yet
            IF (b_inverse_loaded .neqv. .true.) THEN
                CALL load_b_inverse()
            END IF

            ! find the coordinates at which we need to access gridded values
            i_x = AINT((x - x0)/dx, KIND(i_x)) + 1
            i_y = AINT((y - y0)/dy, KIND(i_y)) + 1
            i_z = AINT((z - z0)/dz, KIND(i_z)) + 1

            ! and get the coefficients
            a(:, :, :) = coeffs(i_x, i_y, i_z, :, :, :)

            ! convert x, y, z to nondimensional coordinates
            x_nd = (x - x0)/dx - AINT((x - x0)/dx) 
            y_nd = (y - y0)/dy - AINT((y - y0)/dy)
            z_nd = (z - z0)/dz - AINT((z - z0)/dz) 

            result = 0
            IF ((xder == 0) .AND. (yder == 0) .AND. (zder == 0)) THEN
                DO i = 0, 3
                    DO j = 0, 3
                        DO k = 0, 3
                            result = result + a(i+1, j+1, k+1) * (x_nd ** i) * (y_nd ** j) * (z_nd ** k)
                        END DO
                    END DO
                END DO
            ELSE IF ((xder == 1) .AND. (yder == 0) .AND. (zder == 0)) THEN
                DO i = 0, 3
                    DO j = 0, 3
                        DO k = 0, 3
                            result = result + a(i+1, j+1, k+1) * i*(x_nd ** (i - 1)) * (y_nd ** j) * (z_nd ** k)
                        END DO
                    END DO
                END DO
                ! convert result back to dimensionful units
                result = result / dx
            ELSE IF ((xder == 0) .AND. (yder == 1) .AND. (zder == 0)) THEN
                DO i = 0, 3
                    DO j = 0, 3
                        DO k = 0, 3
                            result = result + a(i+1, j+1, k+1) * j*(x_nd ** i) * (y_nd ** (j - 1)) * (z_nd ** k)
                        END DO
                    END DO
                END DO
                ! convert result back to dimensionful units
                result = result / dy
            ELSE IF ((xder == 1) .AND. (yder == 1) .AND. (zder == 0)) THEN
                DO i = 0, 3
                    DO j = 0, 3
                        DO k = 0, 3
                            result = result + a(i+1, j+1, k+1) * i * j * (x_nd ** (i - 1)) * (y_nd ** (j - 1)) * (z_nd ** k)
                        END DO
                    END DO
                END DO
                ! convert result back to dimensionful units
                result = result / (dx * dy)
            ELSE IF ((xder == 0) .AND. (yder == 2) .AND. (zder == 0)) THEN
                DO i = 0, 3
                    DO j = 0, 3
                        DO k = 0, 3
                            result = result + a(i+1, j+1, k+1) * j * (j - 1) * (x_nd ** i) * (y_nd ** (j - 2)) * (z_nd ** k)
                        END DO
                    END DO
                END DO
                ! convert result back to dimensionful units
                result = result / (dy ** 2)
            END IF

            RETURN
        END SUBROUTINE calc_spline
    
END MODULE tricubic_mod