subroutine floorice_forcing_T(gT_arr, iMin,iMax,jMin,jMax, kLev, bi, bj, myTime, myIter, myThid)

    implicit none

    ! input/output parameters
    ! gT_arr is the temperature-tendency array. It contains contributions from *all* forcing sources
    !   (latent heat, geothermal, etc.) which are summed up in the apply_forcing subroutines. 
    real, intent(out), dimension(1-OLx:sNx+OLx,1-OLy:sNy+OLy) :: gT_arr
    integer, intent(in)     :: iMin,iMax,jMin,jMax, kLev, bi, bj
    real, intent(in)        :: myTime
    integer, intent(in)     :: myIter, myThid

    ! internal variables
    integer :: i, j, Kp1, Km1
    real    :: drLoc, gTloc

    ! other "external" variables
    ! _hFacC: the fraction of a cell which is "wet." Will be < 0 if we have bottom topography or ice.

    if ( FLOORICEboundaryLayer ) then
     do j=1,sNy
      do i=1,sNx
       ! If the lowermost wet cell, kLowC, has a small amount of water, we want to "redistribute" the
       ! forcing so we don't get huge temperature changes in this cell.
       if (kLev > 1 .AND. kLev == kLowC(i,j,bi,bj) ) then
        Km1 = MAX(kLev - 1, 1)
        drLoc = drF(kLev)*(1. _d 0 - _hFacC(i,j,kLev,bi,bj) )   ! will be zero if cell is ice-free!!
        drLoc = MIN( drLoc, drF(Km1) * _hFacC(i,j,Kp1,bi,bj) )  ! depth of water in cell above
        drLoc = MAX( drLoc, 0. _d 0)
        gTloc = flooriceForcingT(i,j,bi,bj)/( drF(kLev)*_hFacC(i,j,kLev,bi,bj) + drLoc )
        gT_arr(i,j) = gT_arr(i,j) + gTloc
       else if ( kLev < Nr .AND. kLev + 1 == kLowC(i,j,bi,bj) ) then
        Kp1 = MIN(kLev + 1, Nr)
        drLoc = drF(Kp1)*( 1. _d 0 - _hFacC(i, j, Kp1, bi, bj) )
        drLoc = MIN(drLoc, drF(kLev) * _hFacC(i,j,kLev,bi,bj) )
        drLoc = MAX( drLoc, 0. _d 0)
        gTloc = flooriceForcingT(i,j,bi,bj)/(drF(Kp1)*_hFacC(i,j,Kp1,bi,bj) + drLoc)
        gT_arr(i,j) = gT_arr(i,j) + gTloc*drLoc*recip_drF(kLev)* _recip_hFacC(i,j,kLev,bi,bj)
       endif
      enddo
     enddo
    endif

end subroutine floorice_forcing_T

subroutine floorice_forcing_S(gS_arr, iMin,iMax,jMin,jMax, kLev, bi, bj, myTime, myIter, myThid)

    implicit none

    ! input/output parameters
    ! gT_arr is the salinity-tendency array. It contains contributions from *all* forcing sources
    !   which are summed up in the apply_forcing subroutines. 
    real, intent(out), dimension(1-OLx:sNx+OLx,1-OLy:sNy+OLy) :: gS_arr
    integer, intent(in)     :: iMin,iMax,jMin,jMax, kLev, bi, bj
    real, intent(in)        :: myTime
    integer, intent(in)     :: myIter, myThid

    ! internal variables
    integer :: i, j, Kp1, Km1
    real    :: drLoc, gSloc

    ! other "external" variables
    ! _hFacC: the fraction of a cell which is "wet." Will be < 0 if we have bottom topography or ice.

    if ( FLOORICEboundaryLayer ) then
     do j=1,sNy
      do i=1,sNx
       ! If the lowermost wet cell, kLowC, has a small amount of water, we want to "redistribute" the
       ! forcing so we don't get huge temperature changes in this cell.
       if (kLev > 1 .AND. kLev == kLowC(i,j,bi,bj) ) then
        Km1 = MAX(kLev - 1, 1)
        drLoc = drF(kLev)*(1. _d 0 - _hFacC(i,j,kLev,bi,bj) )   ! will be zero if cell is ice-free!!
        drLoc = MIN( drLoc, drF(Km1) * _hFacC(i,j,Kp1,bi,bj) )  ! depth of water in cell above
        drLoc = MAX( drLoc, 0. _d 0)
        gSloc = flooriceForcingS(i,j,bi,bj)/( drF(kLev)*_hFacC(i,j,kLev,bi,bj) + drLoc )
        gS_arr(i,j) = gS_arr(i,j) + gSloc
       else if ( kLev < Nr .AND. kLev + 1 == kLowC(i,j,bi,bj) ) then
        Kp1 = MIN(kLev + 1, Nr)
        drLoc = drF(Kp1)*( 1. _d 0 - _hFacC(i, j, Kp1, bi, bj) )
        drLoc = MIN(drLoc, drF(kLev) * _hFacC(i,j,kLev,bi,bj) )
        drLoc = MAX( drLoc, 0. _d 0)
        gSloc = flooriceForcingS(i,j,bi,bj)/(drF(Kp1)*_hFacC(i,j,Kp1,bi,bj) + drLoc)
        gS_arr(i,j) = gS_arr(i,j) + gSloc*drLoc*recip_drF(kLev)* _recip_hFacC(i,j,kLev,bi,bj)
       endif
      enddo
     enddo
    endif

end subroutine floorice_forcing_S

