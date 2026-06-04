#include "HPIMM_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: HPIMM_DIAGNOSTICS_INIT

C !INTERFACE:
      SUBROUTINE HPIMM_DIAGNOSTICS_INIT( myThid )

C     !DESCRIPTION:
C     Initialize list of all available diagnostics for the HPIMM package.

C     !USES:
      IMPLICIT NONE
C     == Global variables ===
#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#include "HPIMM.h"

C     !INPUT/OUTPUT PARAMETERS:
C     myThid ::  my Thread Id number
      INTEGER myThid
CEOP

#ifdef ALLOW_DIAGNOSTICS
C     !LOCAL VARIABLES:
C     diagNum   :: diagnostics number in the (long) list of available diag.
C     diagMate  :: diag. mate number in the (long) list of available diag.
C     diagName  :: local short name (8c) of a diagnostics
C     diagCode  :: local parser field with characteristics of the diagnostics
C              cf head of S/R DIAGNOSTICS_INIT_EARLY or DIAGNOSTICS_MAIN_INIT
C     diagUnits :: local string (16c): physical units of a diagnostic field
C     diagTitle :: local string (80c): description of field in diagnostic
      INTEGER       diagNum
      INTEGER       diagMate
      CHARACTER*8   diagName
      CHARACTER*16  diagCode
      CHARACTER*16  diagUnits
      CHARACTER*(80) diagTitle

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      IF ( useDiagnostics) THEN
        diagName  = 'HpiLQflx'
        diagTitle = 'Lower boundary heat flux'
        diagUnits = 'J m^-2 s^-1     '
        diagCode  = 'SM      L1      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'HpiLSflx'
        diagTitle = 'Lower boundary salt flux'
        diagUnits = 'g/(kg m^2 s)    '
        diagCode  = 'SM      L1      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'HpiLdmdt'
        diagTitle = 'Lower boundary freezing rate'
        diagUnits = 'kg m^-2 s^-1    '
        diagCode  = 'SM      L1      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'HpiLBtmp'
        diagTitle = 'Lower boundary in-situ temperature'
        diagUnits = 'degC            '
        diagCode  = 'SM      L1      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C--     an example of how to set "free levels option" diagnostics
C        diagName  = 'HpiLev3'
C        diagTitle = '3-levels field for package HPIMM'
C        diagUnits = 'unit            '
C        diagCode  = 'SM      MX      '
C        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
C     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
C        CALL DIAGNOSTICS_SETKLEV( diagName, 3, myThid )

      ENDIF

#endif /* ALLOW_DIAGNOSTICS */

      RETURN
      END
