module floorice_thermodynamics_mod

#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "EOS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "FFIELDS.h"

    implicit none
    real, dimension(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy) :: flooriceForcingT, flooriceForcingS, fliTransCoeffT, fliTransCoeffS, Qtidal, ice_dmdt
    character, dimension(3) :: iceType      
    real :: aEOS, bEOS, cEOS              ! array of EOS coefficients
    real :: FLOORICElatentHeat            ! latent heat of melting/freezing
    real :: CpIce, rhoIce, sIce           ! TODO: change to be ice-type specific

    public :: floorice_init_fixed, floorice_forcing_T, floorice_forcing_S, floorice_thermodynamics, write_floorice_diagnostics

contains

! *************************************************************************************************
subroutine floorice_init_fixed( myThid )

    ! initializes variables for floorice calculations
    implicit none
    integer, intent(in) :: myThid

    if trim(iceType) == 'Ih' then
     print *, 'Ice polymporph Ih selected'
    else if trim(iceType) == 'III' then
     print *, 'Ice polymporph III selected'
    else if trim(iceType) == 'V' then
     print *, 'Ice polymporph V selected'
    else if trim(iceType) == 'VI' then
     print *, 'Ice polymporph VI selected'
    end if 

    sIce = 0.0
    CpIce = 2000.0
    rhoIce = 900.0

    do bj = myByLo(myThid), myByHi(myThid)
     do bi = myBxLo(myThid), myBxHi(myThid)
      fliTransCoeffT(:,:,bi,bj) = 1.0
      fliTransCoeffS(:,:,bi,bj) = 1.0
      ice_dmdt(:,:,bi,bj) = 0.0         ! TODO: figure out how to read this in from file
     enddo
    enddo

end subroutine floorice_init()

! *************************************************************************************************
subroutine floorice_forcing_T(gT_arr, iMin,iMax,jMin,jMax, kLev, bi, bj, myTime, myIter, myThid)

    ! Function which adds the contributions of melting/freezing at the lower ice-ocean interface to
    ! the total temperature-tendency array gT_arr.
    implicit none

    ! arguments
    real, intent(inout), dimension(1-OLx:sNx+OLx,1-OLy:sNy+OLy) :: gT_arr
    integer, intent(in)     :: iMin,iMax,jMin,jMax, kLev, bi, bj
    real, intent(in)        :: myTime
    integer, intent(in)     :: myIter, myThid

    ! internal variables
    integer :: i, j, Kp1, Km1
    real    :: drLoc, gTloc

    ! other "external" variables
    ! _hFacC: the fraction of a cell which is "wet." Will be < 0 if we have bottom topography or ice.

!    if ( FLOORICEboundaryLayer ) then
     do j=1,sNy         ! loop over all y-tiles handled by current processor
      do i=1,sNx        ! loop over all x-tiles handled by current processor
       ! If the lowermost wet cell, kLowC, has a small amount of water, we want to "redistribute"
       ! the forcing so we don't get huge temperature changes in this cell. See SHELFICE 
       ! documentation for a detailed discussion:
       ! https://mitgcm.readthedocs.io/en/latest/phys_pkgs/shelfice.html#shelfice-description
       if (kLev > 1 .AND. kLev == kLowC(i,j,bi,bj) ) then       ! if we are in the lowest (at least partially) "wet" cell...
        Km1 = MAX(kLev - 1, 1)
        drLoc = drF(kLev)*(1. _d 0 - _hFacC(i,j,kLev,bi,bj) )   ! will be zero if cell is ice-free!!
        drLoc = MIN( drLoc, drF(Km1) * _hFacC(i,j,Kp1,bi,bj) )  ! depth of water in cell above
        drLoc = MAX( drLoc, 0. _d 0)
        gTloc = flooriceForcingT(i,j,bi,bj)/( drF(kLev)*_hFacC(i,j,kLev,bi,bj) + drLoc )
        gT_arr(i,j) = gT_arr(i,j) + gTloc
       else if ( kLev < Nr .AND. kLev + 1 == kLowC(i,j,bi,bj) ) then    ! if we are in the second-lowest (at least partially) "wet" cell...
        Kp1 = MIN(kLev + 1, Nr)
        drLoc = drF(Kp1)*( 1. _d 0 - _hFacC(i, j, Kp1, bi, bj) )
        drLoc = MIN(drLoc, drF(kLev) * _hFacC(i,j,kLev,bi,bj) )
        drLoc = MAX( drLoc, 0. _d 0)
        gTloc = flooriceForcingT(i,j,bi,bj)/(drF(Kp1)*_hFacC(i,j,Kp1,bi,bj) + drLoc)
        gT_arr(i,j) = gT_arr(i,j) + gTloc*drLoc*recip_drF(kLev)* _recip_hFacC(i,j,kLev,bi,bj)
       endif
      enddo
     enddo
!    endif

end subroutine floorice_forcing_T

! *************************************************************************************************

subroutine floorice_forcing_S(gS_arr, iMin,iMax,jMin,jMax, kLev, bi, bj, myTime, myIter, myThid)

    ! Function which adds the contributions of melting/freezing at the lower ice-ocean interface to
    ! the total temperature-tendency array gT_arr.
    implicit none

    ! arguments
    real, intent(inout), dimension(1-OLx:sNx+OLx,1-OLy:sNy+OLy) :: gS_arr
    integer, intent(in)     :: iMin,iMax,jMin,jMax, kLev, bi, bj
    real, intent(in)        :: myTime
    integer, intent(in)     :: myIter, myThid

    ! internal variables
    integer :: i, j, Kp1, Km1
    real    :: drLoc, gSloc

    ! other "external" variables
    ! _hFacC: the fraction of a cell which is "wet." Will be < 0 if we have bottom topography or ice.

!    if ( FLOORICEboundaryLayer ) then
     do j=1,sNy         ! loop over all y-tiles handled by current processor
      do i=1,sNx        ! loop over all x-tiles handled by current processor
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
!    endif

end subroutine floorice_forcing_S

! *************************************************************************************************
subroutine floorice_thermodynamics(myTime, myIter, myThid)
    ! A subroutine for calculating the temperature and salinity forcings due to floor ice freezing/melting. 
    implicit none

    ! arguments
    real, intent(in)      :: myTime
    integer, intent(in)   :: myIter, myThid

    ! internal variables
    integer     :: I, J, bi, bj         ! loop counters
    integer     :: K
    real, dimension(1-OLx:sNx+OLx,1-OLy:sNy+OLy) :: sLoc, pLoc, rhoLoc, tLoc, tBl, TFluxOc


    do bj = myByLo(myThid), myByHi(myThid)
     do bi = myBxLo(myThid), myBxHi(myThid)
      do J = 1,sNy
       do I = 1, sNx
        K = kLowC(I,J,bi,bj)        ! lowest wet cell
        sLoc(I,J) = MAX(salt(I,J,K,bi,bj), zeroRL)  ! local salinity
        pLoc(I,J) = pRef4EOS(K)                     ! a lowest-order approx. to the local pressure. TODO: make more accurate
        rhoLoc(I,J) = rhoInSitu(I,J,K,bi,bj)        ! local density
        tLoc(I,J) = SW_TEMP(sLoc(I,J),theta(I,J,K,bi,bj),pLoc(I,J),zeroRL)  ! calculate temperature from potential temp.
       enddo
      enddo
      sBL(:,:) = (sLoc(:,:)*rhoLoc(:,:)*gammaOceanS - sIce*ice_dmdt(:,:))/(rhoLoc(:,:)*gammaOceanS - ice_dmdt(:,:))
      tBL(:,:) = aEOS*pLoc(:,:) + bEOS*sBL(:,:) + cEOS
      TFluxOc(:,:) = (rhoLoc(:,:)*gammaOceanT + ice_dmdt(:,:))*(tLoc(:,:) - tBL(:,:))
      Qtidal(:,:,bi,bj) = rhoIce*CpIce*kappaIce*TGradIce - FLOORICELatentHeat*ice_dmdt(:,:) - HeatCapacity_Cp*TFluxOc(:,:)
      flooriceForcingT(:,:,bi,bj) = TFluxOc(:,:)/rhoLoc(:,:)
      flooriceForcingS(:,:,bi,bj) = ice_dmdt(:,:)*(sBL(:,:) - sIce)/rhoLoc(:,:)
     enddo
    enddo

end subroutine floorice_thermodynamics
! *************************************************************************************************
subroutine write_floorice_diagnostics()

    implicit none
    ! TODO: figure out how to write quantities as diagnostics

end subroutine write_floorice_diagnostics

end module floorice_thermodynamics_mod