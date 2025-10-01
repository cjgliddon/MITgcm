subroutine floorice_thermodynamics(myTime, myIter, myThid)
  ! A subroutine for calculating the temperature and salinity forcings due to floor ice freezing/melting. 
  implicit none

#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "EOS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "FFIELDS.h"
! #include "FLOORICE.h"
  
  real, intent(in)      :: myTime
  integer, intent(in)   :: myIter, myThid

  ! Shared variables:
  ! floorIceHeatFlux, floorIceFreshWaterFlux, maskC, fliTransCoeffT
  ! sNy, sNx, kBotC, HeatCapacity_Cp, FLOORICElatentHeat
  ! phasenum (make it fully homogeneous for now)
  ! need a new shared variable iceMassFlux, but for now we make it local...
  ! also aEOS, bEOS, cEOS
  ! TODO: how do we find the proper bi, bj?
  ! ice_dmdt is a prescribed quantity
  
! #ifdef ALLOW_SEAFRZ
  real, dimension(1:sNx, sNy)   :: tLoc, sLoc, pLoc     ! in-situ temperature, salinity, and pressure
  real, dimension(sNx, sNy)     :: iceMassFlux          ! prescribed ice mass flux consistent with ice bathymetry
  real          :: xx_flifwflx_loc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)       ! freshwater flux, 0.0 by default
  real          :: tFreeze, recip_Cp, cFac
  real          :: TransCoeffTConst = 0 d_ 0
  real          :: gammaTurb = 1 d_ 1, gammaTmoleT = 109.0
  real, parameter :: rhoIce = 900.0, CpIce = 2000       ! TODO: get phase-dependent values
  real          :: TGradIce
  real          :: aEOS, bEOS, cEOS
  real          :: kappaIce, gammaOceanT, gammaOceanS, sIce
  integer       :: I, J, K, bi, bj
  real, dimension(1:sNx, sNy)   :: sBL, tBL, Qtidal, heatFluxOc

  cFac = 0.0                  ! TODO: make this some kind of shared/namelist variable
  recip_Cp = 1. _d 0 / HeatCapacity_Cp

  DO J = 1, sNy
    DO I = 1, sNx
      K = kLowC(I,J,bi,bj)
    ! TODO: add something for this
      K = kBotC(I,J,bi,bj)    ! global, set by topography (TODO: instantiate)
!     TODO: think about if this is the form we really want?
      sLoc(I,J) = MAX(salt(I,J,K,bi,bj), zeroRL)
      tLoc(I,J) = theta(I,J,K,bi,bj)
      pLoc(I,J) = pRef4EOS(K)
      rhoLoc(I,J) = rhoInSitu(I,J,K,bi,bj)
      ! convert tLoc to absolute rather than potential temperature
      tLoc(I,J) = SW_TEMP(sLoc(I,J),tLoc(I,J),pLoc(I,J),zeroRL)

! TODO: calculate parameters for linear fit for each phase
      tFreeze = aEOS*pLoc(I,J) + bEOS*sLoc(I,J) + cEOS



! TODO: modify this if we want to make a more sophisticated representation of the transfer coefficient
      TransCoeffTConst = 0.5/(gammaTurb + gammaTmoleT)
      fliTransCoeffT(:,:,:,:) = TransCoeffTConst

      floorIceHeatFlux(I,J,bi,bj) = maskC(I,J,K,bi,bj) &
                  * fliTransCoeffT(i,j,bi,bj)          &
                  * ( tLoc(I,J) - tFreeze )            &
                  * HeatCapacity_Cp*rUnit2mass         &
                  - xx_flifwflx_loc(I,J,bi,bj)*FLOORICElatentHeat     ! will we need to modify this based on ice phases?
      flooriceForcingT(i,j,bi,bj) = - floorIceHeatFlux(I,J,bi,bj)*recip_Cp*mass2rUnit &
                  - cFac * floorIceFreshWaterFlux(I,J,bi,bj) * mass2rUnit * ( tFreeze - tLoc(I,J) )
    ENDDO
  ENDDO
  ! solve three-equation system straightforwardly
  TFluxOc(:,:) = (rhoLoc(:,:)*gammaOceanT + ice_dmdt(:,:))*(tLoc(:,:) - tBL(:,:))
  sBL(:,:) = (sLoc(:,:)*rhoLoc(:,:)*gammaOceanS - sIce*ice_dmdt(:,:))/(rhoLoc(:,:)*gammaOceanS - ice_dmdt(:,:))
  tBL(:,:) = aEOS*pLoc(:,:) + bEOS*sBL(:,:) + cEOS
  Qtidal(:,:) = rhoIce*CpIce*kappaIce*TGradIce - SHELFICELatentHeat*ice_dmdt(:,:) - HeatCapacity_Cp*TFluxOc(:,:)
  flooriceForcingT(:,:) = TFluxOc(:,:)/rhoLoc(:,:)
  flooriceForcingS(:,:) = ice_dmdt(:,:)*(sBL(:,:) - sIce)/rhoLoc(:,:)
  ! TODO: understand the rFWinBL bit from original code. Maybe change from ice_dmdt to a FWFlux?

! #endif /* ALLOW_SEAFRZ */

end subroutine 