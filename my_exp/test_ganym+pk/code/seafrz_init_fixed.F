#include "SEAFRZ_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: SEAFRZ_INIT_FIXED

C !INTERFACE:
      SUBROUTINE SEAFRZ_INIT_FIXED( myThid )

C     !DESCRIPTION:
C     Initialize SEAFRZ variables that are kept fixed during the run.

C     !USES:
      IMPLICIT NONE
C     == Global variables ===
#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#include "SEAFRZ.h"

C     !INPUT/OUTPUT PARAMETERS:
C     myThid ::  my Thread Id number
      INTEGER myThid
CEOP

C     !LOCAL VARIABLES:
      INTEGER iCoef

#ifdef ALLOW_MNC
      IF (useMNC) THEN
        CALL SEAFRZ_MNC_INIT( myThid )
      ENDIF
#endif /* ALLOW_MNC */

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN
        CALL SEAFRZ_DIAGNOSTICS_INIT( myThid )
      ENDIF
#endif

      IF ( sfz_flooriceType .eq. "III" ) THEN
        sfz_FPcoef_l(1) =  2.4577 _d  2
        sfz_FPcoef_l(2) =  3.3174 _d -2
        sfz_FPcoef_l(3) = -2.7079 _d -1
      ELSE IF ( sfz_flooriceType .eq. "V  ") THEN
        sfz_FPcoef_l(1) =  2.4091 _d  2
        sfz_FPcoef_l(2) =  5.5101 _d -2
        sfz_FPcoef_l(3) = -3.4015 _d -1
      ELSE IF ( sfz_flooriceType .eq. "VI ") THEN
        sfz_FPcoef_l(1) =  2.3398 _d  2
        sfz_FPcoef_l(2) =  6.8995 _d -2
        sfz_FPcoef_l(3) = -3.6907 _d -1
      END IF


C     dummy code for unset coefficients
      DO iCoef=1,20
       sfz_EOScoef_w(iCoef) = 0. _d 0
      ENDDO
      DO iCoef=1,5
       sfz_FPcoef_nl(iCoef) = 0. _d 0
      ENDDO      

      RETURN
      END
