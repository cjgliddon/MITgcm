#include "SEAFRZ_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: SEAFRZ_MNC_INIT

C     !INTERFACE:
      SUBROUTINE SEAFRZ_MNC_INIT( myThid )

C     !DESCRIPTION:
C     Initialize MNC lookup tables for Mypackage.

C     !USES:
      IMPLICIT NONE
C     == Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "SEAFRZ.h"

C     !INPUT PARAMETERS:
      INTEGER myThid
CEOP

#ifdef ALLOW_MNC

      IF (useMNC) THEN

        CALL MNC_CW_ADD_VNAME('sfz_StatScal1',
     &      'Cen_xy_Hn__L__t', 4,5,myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('sfz_StatScal1',
     &      'units','[unit]', myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('sfz_StatScal1','description',
     &      'scalar state variable sfz_StatScal1',myThid)

        CALL MNC_CW_ADD_VNAME('sfz_StatScal2',
     &      'Cen_xy_Hn__L__t', 4,5,myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('sfz_StatScal2',
     &      'units','[unit]', myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('sfz_StatScal2','description',
     &      'scalar state variable sfz_StatScal2',myThid)

        CALL MNC_CW_ADD_VNAME('sfz_StatVelU',
     &      'U_xy_Hn__C__t', 4,5,myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('sfz_StatVelU',
     &      'units','[unit]', myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('sfz_StatVelU','description',
     &       'uVel-type state variable sfz_StatVelU',myThid)

        CALL MNC_CW_ADD_VNAME('sfz_StatVelV',
     &      'V_xy_Hn__C__t', 4,5,myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('sfz_StatVelV',
     &      'units','[unit]', myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('sfz_StatVelV','description',
     &       'vVel-type state variable sfz_StatVelV',myThid)

        CALL MNC_CW_ADD_VNAME('sfz_Surf1',
     &      'Cen_xy_Hn__-__t', 3,4,myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('sfz_Surf1',
     &      'units','[unit]', myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('sfz_Surf1','description',
     &       'Surface scalar field sfz_Surf1',myThid)

        CALL MNC_CW_ADD_VNAME('sfz_Surf2',
     &      'Cen_xy_Hn__-__t', 3,4,myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('sfz_Surf2',
     &      'units','[unit]', myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('sfz_Surf2','description',
     &       'Surface scalar field sfz_Surf2',myThid)

      ENDIF

#endif /* ALLOW_MNC */

      RETURN
      END
