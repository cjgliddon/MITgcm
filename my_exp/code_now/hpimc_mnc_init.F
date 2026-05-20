#include "HPIMC_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: HPIMC_MNC_INIT

C     !INTERFACE:
      SUBROUTINE HPIMC_MNC_INIT( myThid )

C     !DESCRIPTION:
C     Initialize MNC lookup tables for HPIMC.

C     !USES:
      IMPLICIT NONE
C     == Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "HPIMC.h"

C     !INPUT PARAMETERS:
      INTEGER myThid
CEOP

#ifdef ALLOW_MNC

      IF (useMNC) THEN

         PRINT *, "Warning: no variables to initialize"

      ENDIF

#endif /* ALLOW_MNC */

      RETURN
      END
