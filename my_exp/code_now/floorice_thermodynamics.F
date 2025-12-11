#include "SHELFICE_OPTIONS.h"
C     Floorice thermodynamics routines.

C     *=============================================================*
CBOP
C     *=============================================================*
      SUBROUTINE floorice_init_fixed(
     I                     myThid )
C     !uses:      
      IMPLICIT NONE
C     !Global variables:
C     !Input/output variables:
      INTEGER myThid
C     !Internal variables
CEOP

      
      END

C     *=============================================================*
CBOP
C     *=============================================================*
      SUBROUTINE floorice_init_varia(
     I                     myThid )
C     !uses:      
      IMPLICIT NONE
C     !Global variables:
C     !Input/output variables:
      INTEGER myThid
C     !Internal variables
CEOP
      END

C     *=============================================================*
CBOP
C     *=============================================================*
      SUBROUTINE floorice_forcing_T(
     U                     gT_arr,
     I                     iMin,iMax,jMin,jMax, kLev, bi, bj,
     I                     myTime, myIter, myThid )
C     !uses:      
      IMPLICIT NONE
C     !Global variables:
#include "SIZE.h"
#include "EEPARAMS.h"
#include "SHELFICE.h"
#include "GRID.h"

C     !Input/output variables:
      INTEGER myIter, myThid
      INTEGER iMin,iMax,jMin,jMax,kLev,bi,bj
      _RL     myTime
      _RL     gT_arr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
C     !Internal variables
      INTEGER i,j,Kp1,Km1
      _RL     drLoc, gTLoc
CEOP
     
      do j=1,sNy
       do i=1,sNx
        if (kLev>1 .AND. kLev == kLowC(i,j,bi,bj)) then
         Km1 = MAX(kLev-1, 1)
         drLoc = drF(kLev)*(1. _d 0 - _hFacC(i,j,kLev,bi,bj))
         drLoc = MIN(drLoc, drF(Km1) * _hFacC(i,j,Km1,bi,bj))
         drLoc = MAX(drLoc, 0. _d 0)
         gTLoc = flooriceForcingT(i,j,bi,bj)/
     &            ( drF(kLev)*_hFacC(i,j,kLev,bi,bj) + drLoc )
         gT_arr(i,j) = gT_arr(i,j) + gTLoc
        else if (kLev<Nr .AND. kLev+1 == kLowC(i,j,bi,bj)) then
         Kp1 = MIN(kLev + 1, Nr)
         drLoc = drF(Kp1)*( 1. _d 0 - _hFacC(i, j, Kp1, bi, bj) )
         drLoc = MIN(drLoc, drF(kLev) * _hFacC(i,j,kLev,bi,bj) )
         drLoc = MAX( drLoc, 0. _d 0)
         gTloc = flooriceForcingT(i,j,bi,bj)/
     &            ( drF(Kp1)*_hFacC(i,j,Kp1,bi,bj) + drLoc )
         gT_arr(i,j) = gT_arr(i,j) + gTloc*drLoc*recip_drF(kLev)* 
     &            _recip_hFacC(i,j,kLev,bi,bj)
        endif
       enddo
      enddo
      END

C     *=============================================================*
CBOP
C     *=============================================================*
      SUBROUTINE floorice_forcing_S(
     U                     gS_arr,
     I                     iMin,iMax,jMin,jMax, kLev, bi, bj,
     I                     myTime, myIter, myThid )
C     !uses:      
      IMPLICIT NONE
C     !Global variables:
#include "SIZE.h"
#include "EEPARAMS.h"
#include "SHELFICE.h"
#include "GRID.h"

C     !Input/output variables:
      INTEGER myIter, myThid
      INTEGER iMin,iMax,jMin,jMax,kLev,bi,bj
      _RL     myTime
      _RL     gS_arr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
C     !Internal variables
      INTEGER i,j,Kp1,Km1
      _RL     drLoc, gSLoc
CEOP
     
      do j=1,sNy
       do i=1,sNx
        if (kLev>1 .AND. kLev == kLowC(i,j,bi,bj)) then
         Km1 = MAX(kLev-1, 1)
         drLoc = drF(kLev)*(1. _d 0 - _hFacC(i,j,kLev,bi,bj))
         drLoc = MIN(drLoc, drF(Km1) * _hFacC(i,j,Km1,bi,bj))
         drLoc = MAX(drLoc, 0. _d 0)
         gSLoc = flooriceForcingS(i,j,bi,bj)/
     &            ( drF(kLev)*_hFacC(i,j,kLev,bi,bj) + drLoc )
         gS_arr(i,j) = gS_arr(i,j) + gSLoc
        else if (kLev<Nr .AND. kLev+1 == kLowC(i,j,bi,bj)) then
         Kp1 = MIN(kLev + 1, Nr)
         drLoc = drF(Kp1)*( 1. _d 0 - _hFacC(i, j, Kp1, bi, bj) )
         drLoc = MIN(drLoc, drF(kLev) * _hFacC(i,j,kLev,bi,bj) )
         drLoc = MAX( drLoc, 0. _d 0)
         gSloc = flooriceForcingS(i,j,bi,bj)/
     &            ( drF(Kp1)*_hFacC(i,j,Kp1,bi,bj) + drLoc )
         gS_arr(i,j) = gS_arr(i,j) + gSloc*drLoc*recip_drF(kLev)* 
     &            _recip_hFacC(i,j,kLev,bi,bj)
        endif
       enddo
      enddo
      END

C     *=============================================================*
CBOP
C     *=============================================================*
      SUBROUTINE floorice_thermodynamics(
     I                     myTime, myIter, myThid )
C     !uses:      
      IMPLICIT NONE
C     !Global variables:
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DYNVARS.h"
#include "SHELFICE.h"
#include "GRID.h"
#include "MNC_PARAMS.h"

C     !Input/output variables:
      _RL     myTime
      INTEGER myIter, myThid
C     !External functions
      _RL     SW_TEMP
      EXTERNAL SW_TEMP
C     !Internal variables
      INTEGER I,J,K,bi,bj
      _RL     sLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     pLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     rhoLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     tLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     tBl(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     sBl(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     TFluxOc(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
C     !Eventually we'll want to make the following variables part of
C     !a namelist, but we keep them "local" for now...
      _RL     Qtidal(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL     ice_dmdt(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     aEOS,bEOS,cEOS
      _RL     CpIce,rhoIce,sIce,FLOORICElatentHeat,kappaIce
      _RL     gammaOceanT, gammaOceanS, TGradIce
C     CHARACTER iceType(3)
CEOP

C     ! initialization
      aEOS = 2.  _d -4
      bEOS = 7.4 _d -4
      cEOS = 0.  _d 0
      sIce = 0.  _d 0
      CpIce = 2. _d 3
      rhoIce = 9. _d 2
      FLOORICElatentHeat = 3.3355 _d 5
      gammaOceanT = 1. _d 0
      gammaOceanS = 1. _d 0
      kappaIce = 1.54 _d -06
C     This is an order of magnitude estimate from Vance et al. 2014
C     TODO: make a namelist variable
      TGradIce = 1. _d -4

      do bj = myByLo(myThid), myByHi(myThid)
       do bi = myBxLo(myThid), myBxHi(myThid)
        do J = 1,sNy
         do I = 1, sNx
          flooriceforcingT(I,J,bi,bj) = 0. _d 0
          flooriceforcingS(I,J,bi,bj) = 0. _d 0
          K = kLowC(I,J,bi,bj)
          if (K > 0) then       
           sLoc(I,J) = MAX(salt(I,J,K,bi,bj), zeroRL)
C     TODO: make pressure more accurate
           pLoc(I,J) = pRef4EOS(K)*1. _d -4                     
C     TODO: Check if this bug fix is reasonable. Since rhoInSitu
C     doesn't seem to have been calculated at this point I wrote
C     this line           
           rhoLoc(I,J) = MAX(rhoInSitu(I,J,K,bi,bj), rhoRef(K))
C           PRINT *, "Diagnostics for tLoc at I,J,K:", I,J,K
C           PRINT *, "Salinity:", sloc(I,J)
C           PRINT *, "PTemp:", theta(I,J,K,bi,bj)
C           PRINT *, "Pressure:", pLoc(I,J)
           tLoc(I,J) = SW_TEMP(sLoc(I,J),theta(I,J,K,bi,bj),
     &                         pLoc(I,J),zeroRL)
C     Derived quantities 
C           PRINT *, "rhoLoc:", rhoLoc(I,J)
C           PRINT *, "ice_dmdt", ice_dmdt(I,J)
           sBL(I,J) = (sLoc(I,J)*rhoLoc(I,J)*gammaOceanS - 
     &      sIce*ice_dmdt(I,J))/(rhoLoc(I,J)*gammaOceanS - 
     &      ice_dmdt(I,J))
           tBL(I,J) = aEOS*pLoc(I,J) + bEOS*sBL(I,J) + cEOS
           TFluxOc(I,J) = (rhoLoc(I,J)*gammaOceanT + 
     &      ice_dmdt(I,J))*(tLoc(I,J) - tBL(I,J))
           Qtidal(I,J,bi,bj) = rhoIce*CpIce*kappaIce*TGradIce - 
     &      FLOORICELatentHeat*ice_dmdt(I,J) - 
     &      HeatCapacity_Cp*TFluxOc(I,J)
C           PRINT *, "Diagnostics at (I, J, bi, bj):", I,J,bi,bj
C           PRINT *, "Local temperature:", tLoc(I,J)
C           PRINT *, "BL temperature:", tBL(I,J)
C           PRINT *, "Local salinity:", sBL(I,J)
C           PRINT *, "Temp flux:", TFluxOc(I,J)
C           PRINT *, "Tidal heating: ", Qtidal(I,J,bi,bj)
           flooriceForcingT(I,J,bi,bj) = TFluxOc(I,J)/rhoLoc(I,J)
           flooriceForcingS(I,J,bi,bj) = ice_dmdt(I,J)*
     &      (sBL(I,J) - sIce)/rhoLoc(I,J)
          endif
         enddo
        enddo
C           TODO: add diagnostics to init_fixed
        CALL DIAGNOSTICS_FILL( TFluxOc, 'TFlux_LB', 
     &            0, 1, 2, bi, bj, myThid )
       enddo
      enddo
      CALL DIAGNOSTICS_FILL( Qtidal, 'QtidalLB', 
     &            0, 1, 0, 1, 1, myThid )
      END
