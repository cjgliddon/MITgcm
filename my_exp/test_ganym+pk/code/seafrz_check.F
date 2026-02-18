#include "SEAFRZ_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: SEAFRZ_CHECK

C !INTERFACE:
      SUBROUTINE SEAFRZ_CHECK( myThid )

C     !DESCRIPTION:
C     Check dependances with other packages

C     !USES:
      IMPLICIT NONE
C     == Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "SEAFRZ.h"

C     !INPUT PARAMETERS:
C     myThid   :: my Thread Id number
      INTEGER myThid
CEOP

C     !LOCAL VARIABLES:
C     msgBuf   :: Informational/error message buffer
      CHARACTER*(MAX_LEN_MBUF) msgBuf

      _BEGIN_MASTER(myThid)

       WRITE(msgBuf,'(A)') 'SEAFRZ_CHECK: #define SEAFRZ'
       CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )

C--  Print out some kee parameters :
       CALL WRITE_0D_L( sfz_StaV_Cgrid, INDEX_NONE,
     &  'sfz_StaV_Cgrid =', ' /* state vector on C-grid */')
       CALL WRITE_0D_L( sfz_Tend_Cgrid, INDEX_NONE,
     &  'sfz_Tend_Cgrid =', ' /* vector tendency on C-grid */')
       CALL WRITE_0D_L( sfz_applyTendT, INDEX_NONE,
     &  'sfz_applyTendT =',
     &  ' /* apply SFZ tendency to Temperature Eq.*/')
       CALL WRITE_0D_L( sfz_applyTendS, INDEX_NONE,
     &  'sfz_applyTendS =',
     &  ' /* apply SFZ tendency to Salinity Eq.*/')
       CALL WRITE_0D_L( sfz_applyTendU, INDEX_NONE,
     &  'sfz_applyTendU =',
     &  ' /* apply SFZ tendency to U momentum Eq.*/')
       CALL WRITE_0D_L( sfz_applyTendV, INDEX_NONE,
     &  'sfz_applyTendV =',
     &  ' /* apply SFZ tendency to V momentum Eq.*/')
C-   add some other as example:
       CALL WRITE_0D_I( sfz_index1, INDEX_NONE, ' sfz_index1 =',
     &  ' /* user defined parameter */')
       CALL WRITE_0D_RL( sfz_param1, INDEX_NONE, ' sfz_param1 =',
     &  ' /* user defined parameter */')
       CALL WRITE_0D_C( sfz_string1,-1,INDEX_NONE,'sfz_string1 =',
     &  ' /* user defined parameter */')

C--  Check parameters:

C-     Put stops here if some flags are wrongly chosen.
C-     For example, require this package to use tempStepping
       IF ( .NOT.tempStepping .AND. sfz_applyTendT ) THEN
        WRITE(msgBuf,'(A)') 'An example requirement is:'
        CALL PRINT_ERROR( msgBuf, myThid )
        WRITE(msgBuf,'(A)') 'sfz_applyTendT needs tempStepping=.true.'
        CALL PRINT_ERROR( msgBuf, myThid )
        CALL ALL_PROC_DIE( 0 )
        STOP 'ABNORMAL END: S/R SEAFRZ_CHECK'
       ENDIF

      _END_MASTER(myThid)

      RETURN
      END
