#include "HPIMM_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: HPIMM_READPARMS

C !INTERFACE:
      SUBROUTINE HPIMM_READPARMS( myThid )

C     !DESCRIPTION:
C     Initialize HPIMM variables and constants.

C     !USES:
      IMPLICIT NONE
C     == Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "HPIMM.h"

C     !INPUT PARAMETERS:
      INTEGER myThid
CEOP

C     !LOCAL VARIABLES:
C     msgBuf     :: Informational/error message buffer
C     iUnit      :: Work variable for IO unit number
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER iUnit

      NAMELIST /HPIMM_PARM01/
     &       hpimm_MNC,
     &       hpimm_QbotType,   hpimm_QbotFile,
     &       hpimm_TbtopFile,  hpimm_TbbotFile,
     &       hpimm_Scal1File,  hpimm_Scal2File,
     &       hpimm_VelUFile,   hpimm_VelVFile,
     &       hpimm_Surf1File,  hpimm_Surf2File,
     &       sfz_RhoFile, sfz_EtaFile,
     &       hpimm_addFloorice, hpimm_doAutoTini,
     &       hpimm_prescribeTb,
     &       hpimm_flooriceType, dTdz_hpi,
     &       kappaT_hpi, gammaT_hpi, gammaS_hpi,
     &       sfz_pRef

     
      IF ( .NOT.useHPIMM ) THEN
C-    pkg HPIMM is not used
        _BEGIN_MASTER(myThid)
C-    Track pkg activation status:
C     print a (weak) warning if data.HPIMM is found
         CALL PACKAGES_UNUSED_MSG( 'useHPIMM', ' ', ' ' )
        _END_MASTER(myThid)
        RETURN
      ENDIF

      _BEGIN_MASTER(myThid)

C--   Default values for HPIMM
      hpimm_MNC         = useMNC
      hpimm_QbotType    = 0
      hpimm_addFloorice  = .TRUE.
      hpimm_doAutoTIni   = .TRUE.
      hpimm_prescribeTb = .TRUE.
      sfz_pRef         = 1.0 _d +05
C-    these values are identical to those in the shelfice pkg
      kappaT_hpi       = 1.54 _d -06
      gammaT_hpi       = 1.0 _d -05
      gammaS_hpi       = 1.0 _d -05
C-    this value is an estimate based on Vance et al. (2014)
      dTdz_hpi         = -1.0 _d -04
      hpimm_flooriceType = 'VI '
C-    file names for initial conditions:
      hpimm_Scal1File   = ' '
      hpimm_Scal2File   = ' '
      hpimm_VelUFile    = ' '
      hpimm_VelVFile    = ' '
      hpimm_Surf1File   = ' '
      hpimm_Surf2File   = ' '
      hpimm_QbotFile    = ' '
      hpimm_TbtopFile   = ' '
      hpimm_TbbotFile   = ' '


      WRITE(msgBuf,'(A)') 'HPIMM_READPARMS: opening data.hpimm'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      CALL OPEN_COPY_DATA_FILE(
     I                     'data.hpimm', 'HPIMM_READPARMS',
     O                     iUnit,
     I                     myThid )

C     Read parameters from open data file
      READ(UNIT=iUnit,NML=HPIMM_PARM01)
      WRITE(msgBuf,'(A)')
     &    'HPIMM_READPARMS: finished reading data.hpimm'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
C     Close the open data file
#ifdef SINGLE_DISK_IO
      CLOSE(iUnit)
#else
      CLOSE(iUnit,STATUS='DELETE')
#endif /* SINGLE_DISK_IO */

C     Make sure that we locally honor the global MNC on/off flag
      hpimm_MNC = hpimm_MNC .AND. useMNC
#ifndef ALLOW_MNC
C     Fix to avoid running without getting any output:
      hpimm_MNC = .FALSE.
#endif
      hpimm_MDSIO = (.NOT. hpimm_MNC) .OR. outputTypesInclusive

      _END_MASTER(myThid)

C--   Everyone else must wait for the parameters to be loaded
      _BARRIER

      RETURN
      END
