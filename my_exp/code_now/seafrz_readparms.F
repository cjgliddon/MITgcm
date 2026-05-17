#include "SEAFRZ_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: SEAFRZ_READPARMS

C !INTERFACE:
      SUBROUTINE SEAFRZ_READPARMS( myThid )

C     !DESCRIPTION:
C     Initialize SEAFRZ variables and constants.

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

C     !LOCAL VARIABLES:
C     msgBuf     :: Informational/error message buffer
C     iUnit      :: Work variable for IO unit number
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER iUnit

      NAMELIST /SEAFRZ_PARM01/
     &       sfz_MNC,
     &       sfz_StaV_Cgrid, sfz_Tend_Cgrid,
     &       sfz_applyTendT, sfz_applyTendS,
     &       sfz_applyTendU, sfz_applyTendV,
     &       sfz_doSwitch1,  sfz_doSwitch2,
     &       sfz_index1,     sfz_index2,
     &       sfz_QbotType,
     &       sfz_param1,     sfz_param2,
     &       sfz_string1,    sfz_string2,
     &       sfz_Scal1File,  sfz_Scal2File,
     &       sfz_VelUFile,   sfz_VelVFile,
     &       sfz_Surf1File,  sfz_Surf2File,
     &       sfz_QbotFile,
     &       sfz_flooriceType, dTdzBot,
     &       sfz_allowFloFrz, sfz_doAutoTini

      IF ( .NOT.useSEAFRZ ) THEN
C-    pkg SEAFRZ is not used
        _BEGIN_MASTER(myThid)
C-    Track pkg activation status:
C     print a (weak) warning if data.seafrz is found
         CALL PACKAGES_UNUSED_MSG( 'useSEAFRZ', ' ', ' ' )
        _END_MASTER(myThid)
        RETURN
      ENDIF

      _BEGIN_MASTER(myThid)

C--   Default values for SEAFRZ
      sfz_MNC         = useMNC
      sfz_StaV_Cgrid  = .TRUE.
      sfz_Tend_Cgrid  = .TRUE.
      sfz_applyTendT  = .FALSE.
      sfz_applyTendS  = .FALSE.
      sfz_applyTendU  = .FALSE.
      sfz_applyTendV  = .FALSE.
C-    additional parameters:
      sfz_doSwitch1   = .FALSE.
      sfz_doSwitch2   = .FALSE.
      sfz_index1      = 0
      sfz_index2      = 0
      sfz_QbotType    = 0
      sfz_param1      = 0. _d 0
      sfz_param2      = 0. _d 0
      sfz_simple       = .TRUE.
      sfz_addFloorice  = .TRUE.
      sfz_allowFloFrz  = .FALSE.
      sfz_doAutoTIni   = .TRUE.
C-    these values are identical to those in the shelfice pkg
      kappaT_hpi       = 1.54 _d -06
      gammaTo_hpi      = 1.0 _d -05
      gammaSo_hpi      = 1.0 _d -05
C-    this value is an estimate based on Vance et al. (2014)
      dTdzBot          = -1.0 _d -04
      sfz_fpParamType  = 'LINEAR'
      sfz_flooriceType = 'VI '
      sfz_string1     = ' '
      sfz_string2     = ' '
C-    file names for initial conditions:
      sfz_EOSFile     = 'sfz_simple_tccoeffs.dat'
      sfz_Scal1File   = ' '
      sfz_Scal2File   = ' '
      sfz_VelUFile    = ' '
      sfz_VelVFile    = ' '
      sfz_Surf1File   = ' '
      sfz_Surf2File   = ' '
      sfz_QbotFile    = ' '

      WRITE(msgBuf,'(A)') 'SEAFRZ_READPARMS: opening data.seafrz'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      CALL OPEN_COPY_DATA_FILE(
     I                     'data.seafrz', 'SEAFRZ_READPARMS',
     O                     iUnit,
     I                     myThid )

C     Read parameters from open data file
      READ(UNIT=iUnit,NML=SEAFRZ_PARM01)
      WRITE(msgBuf,'(A)')
     &    'SEAFRZ_READPARMS: finished reading data.seafrz'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
C     Close the open data file
#ifdef SINGLE_DISK_IO
      CLOSE(iUnit)
#else
      CLOSE(iUnit,STATUS='DELETE')
#endif /* SINGLE_DISK_IO */

C     Make sure that we locally honor the global MNC on/off flag
      sfz_MNC = sfz_MNC .AND. useMNC
#ifndef ALLOW_MNC
C     Fix to avoid running without getting any output:
      sfz_MNC = .FALSE.
#endif
      sfz_MDSIO = (.NOT. sfz_MNC) .OR. outputTypesInclusive

      _END_MASTER(myThid)

C--   Everyone else must wait for the parameters to be loaded
      _BARRIER

      RETURN
      END
