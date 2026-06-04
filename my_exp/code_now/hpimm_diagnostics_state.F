#include "HPIMM_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: HPIMM_DIAGNOSTICS_STATE

C !INTERFACE:
      SUBROUTINE HPIMM_DIAGNOSTICS_STATE(
     I                     myTime, myIter, myThid )

C     !DESCRIPTION:
C     Fill-in diagnostics array for HPIMM state variables

C     !USES:
      IMPLICIT NONE
C     == Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
C  put following include within #ifdef in case this file is copied to local
C  code dir but pkg/hpimm is not compiled (missing file HPIMM.h).
#ifdef ALLOW_HPIMM
#include "HPIMM.h"
#endif

C     !INPUT PARAMETERS:
C     myTime   :: Current time of simulation ( s )
C     myIter   :: Current iteration number in simulation
C     myThid   :: my Thread Id number
      _RL     myTime
      INTEGER myIter, myThid
CEOP

#ifdef ALLOW_HPIMM
C     !LOCAL VARIABLES:

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN

#ifdef HPIMM_3D_STATE
C         CALL DIAGNOSTICS_FILL( hpimm_StatScal1,'SfzSta1',
C     &                          0,Nr, 0, 1, 1, myThid )
#endif /* HPIMM_3D_STATE */

#ifdef HPIMM_2D_STATE
         CALL DIAGNOSTICS_FILL( hpimm_floorice_Qflx,   'HpiLQflx',
     &                          0, 1, 0, 1, 1, myThid )
         CALL DIAGNOSTICS_FILL( hpimm_floorice_Sflx,   'HpiLSflx',
     &                          0, 1, 0, 1, 1, myThid )
         CALL DIAGNOSTICS_FILL( hpimm_floorice_dmdt,   'HpiLdmdt',
     &                          0, 1, 0, 1, 1, myThid )
         CALL DIAGNOSTICS_FILL( hpimm_floorice_tBL,    'HpiLBtmp',
     &                          0, 1, 0, 1, 1, myThid )


#endif /* HPIMM_2D_STATE */

      ENDIF
#endif /* ALLOW_DIAGNOSTICS */

#endif /* ALLOW_HPIMM */

      RETURN
      END
