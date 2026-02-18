#include "SEAFRZ_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: SEAFRZ_DIAGNOSTICS_STATE

C !INTERFACE:
      SUBROUTINE SEAFRZ_DIAGNOSTICS_STATE(
     I                     myTime, myIter, myThid )

C     !DESCRIPTION:
C     Fill-in diagnostics array for SEAFRZ state variables

C     !USES:
      IMPLICIT NONE
C     == Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
C  put following include within #ifdef in case this file is copied to local
C  code dir but pkg/seafrz is not compiled (missing file SEAFRZ.h).
#ifdef ALLOW_SEAFRZ
#include "SEAFRZ.h"
#endif

C     !INPUT PARAMETERS:
C     myTime   :: Current time of simulation ( s )
C     myIter   :: Current iteration number in simulation
C     myThid   :: my Thread Id number
      _RL     myTime
      INTEGER myIter, myThid
CEOP

#ifdef ALLOW_SEAFRZ
C     !LOCAL VARIABLES:

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN

#ifdef SEAFRZ_3D_STATE
         CALL DIAGNOSTICS_FILL( sfz_StatScal1,'SfzSta1',
     &                          0,Nr, 0, 1, 1, myThid )
         CALL DIAGNOSTICS_FILL( sfz_StatScal2,'SfzSta2',
     &                          0,Nr, 0, 1, 1, myThid )
         CALL DIAGNOSTICS_FILL( sfz_StatVelU, 'SfzStaU',
     &                          0,Nr, 0, 1, 1, myThid )
         CALL DIAGNOSTICS_FILL( sfz_StatVelV, 'SfzStaV',
     &                          0,Nr, 0, 1, 1, myThid )
#endif /* SEAFRZ_3D_STATE */

#ifdef SEAFRZ_2D_STATE
         CALL DIAGNOSTICS_FILL( sfz_Surf1 ,   'SfzSur1',
     &                          0, 1, 0, 1, 1, myThid )
         CALL DIAGNOSTICS_FILL( sfz_Surf2 ,   'SfzSur2',
     &                          0, 1, 0, 1, 1, myThid )
         CALL DIAGNOSTICS_FILL( sfz_floorice_Qflx,   'SfzLQflx',
     &                          0, 1, 0, 1, 1, myThid )
         CALL DIAGNOSTICS_FILL( sfz_floorice_Sflx,   'SfzLSflx',
     &                          0, 1, 0, 1, 1, myThid )
         CALL DIAGNOSTICS_FILL( sfz_floorice_dmdt,   'SfzLdmdt',
     &                          0, 1, 0, 1, 1, myThid )



#endif /* SEAFRZ_2D_STATE */

      ENDIF
#endif /* ALLOW_DIAGNOSTICS */

#endif /* ALLOW_SEAFRZ */

      RETURN
      END
