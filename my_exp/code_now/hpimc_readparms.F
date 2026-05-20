#include "HPIMC_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: HPIMC_READPARMS

C !INTERFACE:
      SUBROUTINE HPIMC_READPARMS( myThid )

C     !DESCRIPTION:
C     Initialize HPIMC variables and constants.

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

C     !LOCAL VARIABLES:
C     msgBuf     :: Informational/error message buffer
C     iUnit      :: Work variable for IO unit number
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER iUnit

      NAMELIST /HPIMC_PARM01/
     &       hpimc_MNC,
     &       hpimc_QbotType,   hpimc_QbotFile,
     &       hpimc_Scal1File,  hpimc_Scal2File,
     &       hpimc_VelUFile,   hpimc_VelVFile,
     &       hpimc_Surf1File,  hpimc_Surf2File,
     &       hpimc_addFloorice, hpimc_doAutoTini,
     &       hpimc_flooriceType, dTdz_hpi,
     &       kappaT_hpi, gammaT_hpi, gammaS_hpi

     
      IF ( .NOT.useHPIMC ) THEN
C-    pkg HPIMC is not used
        _BEGIN_MASTER(myThid)
C-    Track pkg activation status:
C     print a (weak) warning if data.HPIMC is found
         CALL PACKAGES_UNUSED_MSG( 'useHPIMC', ' ', ' ' )
        _END_MASTER(myThid)
        RETURN
      ENDIF

      _BEGIN_MASTER(myThid)

C--   Default values for HPIMC
      hpimc_MNC         = useMNC
      hpimc_QbotType    = 0
      hpimc_addFloorice  = .TRUE.
      hpimc_doAutoTIni   = .TRUE.
C-    these values are identical to those in the shelfice pkg
      kappaT_hpi       = 1.54 _d -06
      gammaT_hpi       = 1.0 _d -05
      gammaS_hpi       = 1.0 _d -05
C-    this value is an estimate based on Vance et al. (2014)
      dTdz_hpi         = -1.0 _d -04
      hpimc_flooriceType = 'VI '
C-    file names for initial conditions:
      hpimc_Scal1File   = ' '
      hpimc_Scal2File   = ' '
      hpimc_VelUFile    = ' '
      hpimc_VelVFile    = ' '
      hpimc_Surf1File   = ' '
      hpimc_Surf2File   = ' '
      hpimc_QbotFile    = ' '

      WRITE(msgBuf,'(A)') 'HPIMC_READPARMS: opening data.hpimc'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      CALL OPEN_COPY_DATA_FILE(
     I                     'data.hpimc', 'HPIMC_READPARMS',
     O                     iUnit,
     I                     myThid )

C     Read parameters from open data file
      READ(UNIT=iUnit,NML=HPIMC_PARM01)
      WRITE(msgBuf,'(A)')
     &    'HPIMC_READPARMS: finished reading data.hpimc'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
C     Close the open data file
#ifdef SINGLE_DISK_IO
      CLOSE(iUnit)
#else
      CLOSE(iUnit,STATUS='DELETE')
#endif /* SINGLE_DISK_IO */

C     Make sure that we locally honor the global MNC on/off flag
      hpimc_MNC = hpimc_MNC .AND. useMNC
#ifndef ALLOW_MNC
C     Fix to avoid running without getting any output:
      hpimc_MNC = .FALSE.
#endif
      hpimc_MDSIO = (.NOT. hpimc_MNC) .OR. outputTypesInclusive

      _END_MASTER(myThid)

C--   Everyone else must wait for the parameters to be loaded
      _BARRIER

      RETURN
      END
