#include "SEAFRZ_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: SEAFRZ_DIAGNOSTICS_INIT

C !INTERFACE:
      SUBROUTINE SEAFRZ_DIAGNOSTICS_INIT( myThid )

C     !DESCRIPTION:
C     Initialize list of all available diagnostics

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

      IF ( useDiagnostics ) THEN
        diagName  = 'SfzSta1'
        diagTitle = 'Scalar state variable 1 for package SEAFRZ'
        diagUnits = 'unit            '
        diagCode  = 'SM      MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'SfzSta2'
        diagTitle = 'Scalar state variable 2 for package SEAFRZ'
        diagUnits = 'unit            '
        diagCode  = 'SM      MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'SfzStaU'
        diagTitle = 'uVel-type state variable for package SEAFRZ'
        diagUnits = 'unit            '
        diagCode  = 'UM      MR      '
        diagMate  = diagNum + 2
        IF ( sfz_StaV_Cgrid ) diagCode(2:2) = 'U'
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I   diagName, diagCode, diagUnits, diagTitle, diagMate, myThid )

        diagName  = 'SfzStaV'
        diagTitle = 'vVel-type state variable for package SEAFRZ'
        diagUnits = 'unit            '
        diagCode  = 'VM      MR      '
        diagMate  = diagNum
        IF ( sfz_StaV_Cgrid ) diagCode(2:2) = 'V'
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I   diagName, diagCode, diagUnits, diagTitle, diagMate, myThid )

        diagName  = 'SfzSur1 '
        diagTitle = 'Surface field 1 for package SEAFRZ'
        diagUnits = 'unit            '
        diagCode  = 'SM      L1      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'SfzSur2 '
        diagTitle = 'Surface field 2 for package SEAFRZ'
        diagUnits = 'unit            '
        diagCode  = 'SM      L1      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'SfzLQflx'
        diagTitle = 'Lower boundary heat flux'
        diagUnits = 'J m^-2 s^-1     '
        diagCode  = 'SM      L1      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'SfzLSflx'
        diagTitle = 'Lower boundary salt flux'
        diagUnits = 'g/(kg m^2 s)    '
        diagCode  = 'SM      L1      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'SfzLdmdt'
        diagTitle = 'Lower boundary freezing rate'
        diagUnits = 'kg m^-2 s^-1    '
        diagCode  = 'SM      L1      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )



        diagName  = 'SfzdTdt'
        diagTitle = ' Pot.Temp. tendency from package SEAFRZ'
        diagUnits = 'K/s             '
        diagCode  = 'SMR     MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'SfzdSdt'
        diagTitle = ' Salinity  tendency from package SEAFRZ'
        diagUnits = 'g/kg/s          '
        diagCode  = 'SMR     MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'SfzdUdt'
        diagTitle = 'u-Velocity tendency from package SEAFRZ'
        diagUnits = 'm/s^2           '
        diagCode  = 'UMR     MR      '
        diagMate  = diagNum + 2
        IF ( sfz_Tend_Cgrid ) diagCode(2:2) = 'U'
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I   diagName, diagCode, diagUnits, diagTitle, diagMate, myThid )

        diagName  = 'SfzdVdt'
        diagTitle = 'v-Velocity tendency from package SEAFRZ'
        diagUnits = 'm/s^2           '
        diagCode  = 'VMR     MR      '
        diagMate  = diagNum
        IF ( sfz_Tend_Cgrid ) diagCode(2:2) = 'V'
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I   diagName, diagCode, diagUnits, diagTitle, diagMate, myThid )

C--     an example of how to set "free levels option" diagnostics
        diagName  = 'SfzLev3'
        diagTitle = '3-levels field for package SEAFRZ'
        diagUnits = 'unit            '
        diagCode  = 'SM      MX      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
        CALL DIAGNOSTICS_SETKLEV( diagName, 3, myThid )

      ENDIF

#endif /* ALLOW_DIAGNOSTICS */

      RETURN
      END
