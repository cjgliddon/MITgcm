#include "SEAFRZ_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: SEAFRZ_INIT_VARIA

C !INTERFACE:
      SUBROUTINE SEAFRZ_INIT_VARIA( myThid )

C     !DESCRIPTION:
C     Initialize SEAFRZ variables

C     !USES:
      IMPLICIT NONE
C     == Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "SEAFRZ.h"

C     !INPUT/OUTPUT PARAMETERS:
C     myThid ::  my Thread Id number
      INTEGER myThid
CEOP

C     !LOCAL VARIABLES:
      INTEGER i,j,k,bi,bj

      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)

C     Initialize arrays in common blocks :
C     3-dim.
        DO k=1,Nr
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
#ifdef SEAFRZ_3D_STATE
           sfz_StatScal1(i,j,k,bi,bj) = 0. _d 0
           sfz_StatScal2(i,j,k,bi,bj) = 0. _d 0
           sfz_StatVelU (i,j,k,bi,bj) = 0. _d 0
           sfz_StatVelV (i,j,k,bi,bj) = 0. _d 0
#endif /* SEAFRZ_3D_STATE */
#ifdef SEAFRZ_TENDENCY
           sfz_TendScal1(i,j,k,bi,bj) = 0. _d 0
           sfz_TendScal2(i,j,k,bi,bj) = 0. _d 0
           sfz_TendVelU (i,j,k,bi,bj) = 0. _d 0
           sfz_TendVelV (i,j,k,bi,bj) = 0. _d 0
#endif /* SEAFRZ_TENDENCY */
          ENDDO
         ENDDO
        ENDDO

C     2-dim.
#ifdef SEAFRZ_2D_STATE
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          sfz_Surf1(i,j,bi,bj) = 0. _d 0
          sfz_Surf2(i,j,bi,bj) = 0. _d 0
          sfz_floorice_tBL(i,j,bi,bj)  = 0. _d 0
          sfz_floorice_sBL(i,j,bi,bj)  = 0. _d 0
          sfz_floorice_Qflx(i,j,bi,bj) = 0. _d 0
          sfz_floorice_dmdt(i,j,bi,bj) = 0. _d 0
          sfz_floorice_Sflx(i,j,bi,bj) = 0. _d 0
          sfz_flooriceForcingT(i,j,bi,bj) = 0. _d 0
          sfz_flooriceForcingS(i,j,bi,bj) = 0. _d 0
         ENDDO
        ENDDO
#endif /* SEAFRZ_2D_STATE */

C- end bi,bj loops
       ENDDO
      ENDDO

      IF ( startTime.EQ.baseTime .AND. nIter0.EQ.0
     &                           .AND. pickupSuff.EQ.' ' ) THEN

C     Prevent loading of new data before everyone has finished initialisation
C     this is done within IO routines => no longer needed
c      _BARRIER

#ifdef SEAFRZ_3D_STATE
       IF ( sfz_Scal1File .NE. ' ' ) THEN
         CALL READ_FLD_XYZ_RL( sfz_Scal1File, ' ',
     &                         sfz_StatScal1, 0, myThid )
         CALL EXCH_3D_RL( sfz_StatScal1, Nr, myThid )
       ENDIF
       IF ( sfz_Scal2File .NE. ' ' ) THEN
         CALL READ_FLD_XYZ_RL( sfz_Scal2File, ' ',
     &                         sfz_StatScal2, 0, myThid )
         CALL EXCH_3D_RL( sfz_StatScal2, Nr, myThid )
       ENDIF

       IF ( sfz_VelUFile .NE. ' ' ) THEN
         CALL READ_FLD_XYZ_RL( sfz_VelUFile, ' ',
     &                         sfz_StatVelU, 0, myThid )
       ENDIF
       IF ( sfz_VelVFile .NE. ' ' ) THEN
         CALL READ_FLD_XYZ_RL( sfz_VelVFile, ' ',
     &                         sfz_StatVelV, 0, myThid )
       ENDIF
       IF ( sfz_VelUFile .NE. ' ' .OR.
     &      sfz_VelVFile .NE. ' ' ) THEN
        IF ( sfz_StaV_Cgrid ) THEN
         CALL EXCH_UV_3D_RL( sfz_StatVelU, sfz_StatVelV,
     &                       .TRUE., Nr, myThid )
        ELSE
C-      Assume Agrid position:
         CALL EXCH_UV_AGRID_3D_RL( sfz_StatVelU, sfz_StatVelV,
     &                       .TRUE., Nr, myThid )
        ENDIF
       ENDIF
#endif /* SEAFRZ_3D_STATE */

#ifdef SEAFRZ_2D_STATE
       IF ( sfz_Surf1File .NE. ' ' ) THEN
         CALL READ_FLD_XY_RL( sfz_Surf1File,' ',
     &                        sfz_Surf1, 0, myThid )
         CALL EXCH_XY_RL( sfz_Surf1, myThid )
       ENDIF
       IF ( sfz_Surf2File .NE. ' ' ) THEN
         CALL READ_FLD_XY_RL( sfz_Surf2File,' ',
     &                        sfz_Surf2, 0, myThid )
         CALL EXCH_XY_RL( sfz_Surf2, myThid )
       ENDIF
#endif /* SEAFRZ_2D_STATE */

      ELSE
C--   restart from a pickup:
         CALL SEAFRZ_READ_PICKUP( nIter0, myThid )

C--   end start-from-iter-zero if/else block
      ENDIF

      RETURN
      END
