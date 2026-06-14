#include "HPIMM_OPTIONS.h"

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
#include "GRID.h"
#include "HPIMM.h"

C
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
         gTLoc = sfz_flooriceForcingT(i,j,bi,bj)/
     &            ( drF(kLev)*_hFacC(i,j,kLev,bi,bj) + drLoc )
         gT_arr(i,j) = gT_arr(i,j) + gTLoc
        else if (kLev<Nr .AND. kLev+1 == kLowC(i,j,bi,bj)) then
         Kp1 = MIN(kLev + 1, Nr)
         drLoc = drF(Kp1)*( 1. _d 0 - _hFacC(i, j, Kp1, bi, bj) )
         drLoc = MIN(drLoc, drF(kLev) * _hFacC(i,j,kLev,bi,bj) )
         drLoc = MAX( drLoc, 0. _d 0)
         gTloc = sfz_flooriceForcingT(i,j,bi,bj)/
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
#include "GRID.h"
#include "HPIMM.h"

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
         gSLoc = sfz_flooriceForcingS(i,j,bi,bj)/
     &            ( drF(kLev)*_hFacC(i,j,kLev,bi,bj) + drLoc )
         gS_arr(i,j) = gS_arr(i,j) + gSLoc
        else if (kLev<Nr .AND. kLev+1 == kLowC(i,j,bi,bj)) then
         Kp1 = MIN(kLev + 1, Nr)
         drLoc = drF(Kp1)*( 1. _d 0 - _hFacC(i, j, Kp1, bi, bj) )
         drLoc = MIN(drLoc, drF(kLev) * _hFacC(i,j,kLev,bi,bj) )
         drLoc = MAX( drLoc, 0. _d 0)
         gSloc = sfz_flooriceForcingS(i,j,bi,bj)/
     &            ( drF(Kp1)*_hFacC(i,j,Kp1,bi,bj) + drLoc )
         gS_arr(i,j) = gS_arr(i,j) + gSloc*drLoc*recip_drF(kLev)* 
     &            _recip_hFacC(i,j,kLev,bi,bj)
        endif
       enddo
      enddo
      END
