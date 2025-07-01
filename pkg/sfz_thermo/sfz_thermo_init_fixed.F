#include "SFZ_THERMO_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: SFZ_THERMO_INIT_FIXED

C !INTERFACE:
      SUBROUTINE SFZ_THERMO_INIT_FIXED( myThid )

C     !DESCRIPTION:
C     Initialize SFZ_THERMO variables that are kept fixed during the run.

C     !USES:
      IMPLICIT NONE
C     == Global variables ===
#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#include "SFZ_THERMO.h"

C     !INPUT/OUTPUT PARAMETERS:
C     myThid ::  my Thread Id number
      INTEGER myThid
CEOP

C     !LOCAL VARIABLES:

#ifdef ALLOW_MNC
      IF (useMNC) THEN
        CALL SFZ_THERMO_MNC_INIT( myThid )
      ENDIF
#endif /* ALLOW_MNC */

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN
        CALL SFZ_THERMO_DIAGNOSTICS_INIT( myThid )
      ENDIF
#endif

      RETURN
      END
