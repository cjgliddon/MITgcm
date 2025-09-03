subroutine floorice_thermodynamics(myTime, myIter, myThid)
  ! A subroutine for calculating the heat and mass fluxes associated with seafloor ice. 
  implicit none

#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "EOS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "FFIELDS.h"
#include "FLOORICE.h"
  
  real, intent(in)      :: myTime
  integer, intent(in)   :: myIter, myThid

  ! Shared variables:
  ! floorIceHeatFlux, floorIceFreshWaterFlux, maskC, fliTransCoeffT
  ! sNy, sNx, kBotC, HeatCapacity_Cp, FLOORICElatentHeat
  ! phasenum (make it fully homogeneous for now)

#ifdef ALLOW_FLOORICE
  real, dimension(1:sNx, sNy)   :: tLoc, sLoc, pLoc     ! in-situ temperature, salinity, and pressure
  real          :: xx_flifwflx_loc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)       ! freshwater flux, 0.0 by default
  real          :: tFreeze, recip_Cp, cFac
  real          :: TransCoeffTConst = 0 d_ 0
  real          :: gammaTurb = 1 d_ 1, gammaTmoleT = 109.0
  integer       :: I, J, K, bi, bj

  cFac = 0.0                  ! TODO: make this some kind of shared/namelist variable
  recip_Cp = 1. _d 0 / HeatCapacity_Cp

  DO J = 1, sNy
    DO I = 1, sNx
    ! TODO: add something for this
      K = kBotC(I,J,bi,bj)    ! global, set by topography (TODO: instantiate)

!     TODO: write pressure equation. Need to make locPres, the output of pressure_for_eos.F, a
!     global variable... either that or call it here. 
!     call PRESSURE_FOR_EOS(bi, bj, iMin, iMax, jMin, jMax,  k, dpRef, locPres, myThid ))
      sLoc(I,J) = MAX(salt(I,J,K,bi,bj), zeroRL)
      tLoc(I,J) = theta(I,J,K,bi,bj)
      ! convert tLoc to absolute rather than potential temperature
      tLoc(I,J) = SW_TEMP(sLoc(I,J),tLoc(I,J),pLoc(I,J),zeroRL)

! TODO: Calculate tFreeze using Gibbs data splines
      call getTFreezeSFZ(phasenum, pLoc(I,J), sLoc(I,J), tFreeze)

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
      

#endif /* ALLOW_FLOORICE */

end subroutine 