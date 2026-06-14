#include "HPIMM_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: HPIMM_INIT_VARIA

C !INTERFACE:
      SUBROUTINE HPIMM_INIT_VARIA( myThid )

C     !DESCRIPTION:
C     Initialize HPIMM variables

C     !USES:
      IMPLICIT NONE
C     == Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "HPIMM.h"

C     !INPUT/OUTPUT PARAMETERS:
C     myThid ::  my Thread Id number
      INTEGER myThid
CEOP

C     !LOCAL VARIABLES:
      INTEGER i,j,k,bi,bj

      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)

C     Initialize arrays in common blocks :
C     3-dim. = nothing to initialize; the commented-out loop is an example
C     of the structure
C        DO k=1,Nr
C         DO j=1-OLy,sNy+OLy
C          DO i=1-OLx,sNx+OLx
C#ifdef HPIMM_3D_STATE
C           hpimm_StatScal1(i,j,k,bi,bj) = 0. _d 0
C#endif /* HPIMM_3D_STATE */
C#ifdef HPIMM_TENDENCY
C           hpimm_TendScal1(i,j,k,bi,bj) = 0. _d 0
C#endif /* HPIMM_TENDENCY */
C          ENDDO
C         ENDDO
C        ENDDO

C     2-dim.
#ifdef HPIMM_2D_STATE
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          hpimm_floorice_tBL(i,j,bi,bj)  = 0. _d 0
          hpimm_floorice_sBL(i,j,bi,bj)  = 0. _d 0
          hpimm_floorice_Qflx(i,j,bi,bj) = 0. _d 0
          hpimm_floorice_dmdt(i,j,bi,bj) = 0. _d 0
          hpimm_floorice_Sflx(i,j,bi,bj) = 0. _d 0
          hpimm_flooriceForcingT(i,j,bi,bj) = 0. _d 0
          hpimm_flooriceForcingS(i,j,bi,bj) = 0. _d 0
         ENDDO
        ENDDO
#endif /* HPIMM_2D_STATE */

C- end bi,bj loops
       ENDDO
      ENDDO

      IF ( startTime.EQ.baseTime .AND. nIter0.EQ.0
     &                           .AND. pickupSuff.EQ.' ' ) THEN

C     Prevent loading of new data before everyone has finished initialisation
C     this is done within IO routines => no longer needed
c      _BARRIER

#ifdef HPIMM_3D_STATE
       IF ( hpimm_Scal1File .NE. ' ' ) THEN
         CALL READ_FLD_XYZ_RL( hpimm_Scal1File, ' ',
     &                         hpimm_StatScal1, 0, myThid )
         CALL EXCH_3D_RL( hpimm_StatScal1, Nr, myThid )
       ENDIF
       IF ( hpimm_Scal2File .NE. ' ' ) THEN
         CALL READ_FLD_XYZ_RL( hpimm_Scal2File, ' ',
     &                         hpimm_StatScal2, 0, myThid )
         CALL EXCH_3D_RL( hpimm_StatScal2, Nr, myThid )
       ENDIF

       IF ( hpimm_VelUFile .NE. ' ' ) THEN
         CALL READ_FLD_XYZ_RL( hpimm_VelUFile, ' ',
     &                         hpimm_StatVelU, 0, myThid )
       ENDIF
       IF ( hpimm_VelVFile .NE. ' ' ) THEN
         CALL READ_FLD_XYZ_RL( hpimm_VelVFile, ' ',
     &                         hpimm_StatVelV, 0, myThid )
       ENDIF
       IF ( hpimm_VelUFile .NE. ' ' .OR.
     &      hpimm_VelVFile .NE. ' ' ) THEN
        IF ( hpimm_StaV_Cgrid ) THEN
         CALL EXCH_UV_3D_RL( hpimm_StatVelU, hpimm_StatVelV,
     &                       .TRUE., Nr, myThid )
        ELSE
C-      Assume Agrid position:
         CALL EXCH_UV_AGRID_3D_RL( hpimm_StatVelU, hpimm_StatVelV,
     &                       .TRUE., Nr, myThid )
        ENDIF
       ENDIF
#endif /* HPIMM_3D_STATE */

#ifdef HPIMM_2D_STATE
       IF ( hpimm_Surf1File .NE. ' ' ) THEN
         CALL READ_FLD_XY_RL( hpimm_Surf1File,' ',
     &                        hpimm_Surf1, 0, myThid )
         CALL EXCH_XY_RL( hpimm_Surf1, myThid )
       ENDIF
       IF ( hpimm_Surf2File .NE. ' ' ) THEN
         CALL READ_FLD_XY_RL( hpimm_Surf2File,' ',
     &                        hpimm_Surf2, 0, myThid )
         CALL EXCH_XY_RL( hpimm_Surf2, myThid )
       ENDIF
#endif /* HPIMM_2D_STATE */

      ELSE
C--   restart from a pickup:
         CALL HPIMM_READ_PICKUP( nIter0, myThid )

C--   end start-from-iter-zero if/else block
      ENDIF

      RETURN
      END
