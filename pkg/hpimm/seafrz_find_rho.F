#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"
#include "HPIMM_OPTIONS.h"

      SUBROUTINE FIND_RHO_SFZ_SCALAR(
     I                      ct, sa, p,
     O                      rho )

#include "SIZE.h"
#include "EEPARAMS.h"
#include "HPIMM.h"

C     Arguments
      _RL ct, sa, p
      _RL rho
      
      _RL SFZ_TEMP
      EXTERNAL SFZ_TEMP

C     Local variables
      _RL PaTodb,PaToMPa
      _RL tReal,rhoRecip
      
      PaTodb  = 1.0     _d -4
      PaToMPa = 1.0     _d -6
C      dbToMPa = 1.0     _d -2

      CALL TRIC_EVAL(sfz_rhoCoefs, p, ct+273.15, sa,
     &               0, 0, 0,
     &               rho)
      RETURN

      END SUBROUTINE


      SUBROUTINE FIND_RHO_SFZ_2D(
     I              iMin, iMax, jMin, jMax,
     I              locPres, tFld, sFld,
     O              rhoLoc,
     I              myThid )

C     !USES:
      IMPLICIT NONE
C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "EOS.h"
#include "HPIMM.h"


C     == Routine arguments ==
      INTEGER iMin,iMax,jMin,jMax
      _RL locPres(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL tFld   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sFld   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL rhoLoc (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER myThid
C     == Functions
      _RL SFZ_TEMP
      EXTERNAL SFZ_TEMP
C     == Local variables ==
      INTEGER i,j
      _RL     ct, sa, p, tReal, rho
      _RL     Patodb,PaToMPa

      Patodb  = 1.0 _d -4
      PaToMPa = 1.0 _d -6

      DO j=jMin,jMax
       DO i=iMin,iMax
C       convert quantities to K, molal, MPa
        ct = tFld(i,j)
C       PRINT *, "Conservative temp: ", ct
        sa = sFld(i,j) 
        p  = locPres(i,j)
        CALL TRIC_EVAL(sfz_rhoCoefs, p, ct+273.15, sa,
     &                 0, 0, 0,
     &                 rho)

        rhoLoc(i,j) = rho
C        PRINT *, "p,rho =",p,rhoLoc(i,j)
       END DO
      END DO

      RETURN

      END SUBROUTINE
