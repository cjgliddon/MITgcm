#include "HPIMC_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: HPIMC_INIT_VARIA

C !INTERFACE:
      SUBROUTINE HPIMC_INIT_VARIA( myThid )

C     !DESCRIPTION:
C     Initialize HPIMC variables

C     !USES:
      IMPLICIT NONE
C     == Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "HPIMC.h"

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
C#ifdef HPIMC_3D_STATE
C           hpimc_StatScal1(i,j,k,bi,bj) = 0. _d 0
C#endif /* HPIMC_3D_STATE */
C#ifdef HPIMC_TENDENCY
C           hpimc_TendScal1(i,j,k,bi,bj) = 0. _d 0
C#endif /* HPIMC_TENDENCY */
C          ENDDO
C         ENDDO
C        ENDDO

C     2-dim.
#ifdef HPIMC_2D_STATE
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          hpimc_floorice_tBL(i,j,bi,bj)  = 0. _d 0
          hpimc_floorice_sBL(i,j,bi,bj)  = 0. _d 0
          hpimc_floorice_Qflx(i,j,bi,bj) = 0. _d 0
          hpimc_floorice_dmdt(i,j,bi,bj) = 0. _d 0
          hpimc_floorice_Sflx(i,j,bi,bj) = 0. _d 0
          hpimc_flooriceForcingT(i,j,bi,bj) = 0. _d 0
          hpimc_flooriceForcingS(i,j,bi,bj) = 0. _d 0
         ENDDO
        ENDDO
#endif /* HPIMC_2D_STATE */

C- end bi,bj loops
       ENDDO
      ENDDO

      IF ( startTime.EQ.baseTime .AND. nIter0.EQ.0
     &                           .AND. pickupSuff.EQ.' ' ) THEN

C     Prevent loading of new data before everyone has finished initialisation
C     this is done within IO routines => no longer needed
c      _BARRIER

#ifdef HPIMC_3D_STATE
       IF ( hpimc_Scal1File .NE. ' ' ) THEN
         CALL READ_FLD_XYZ_RL( hpimc_Scal1File, ' ',
     &                         hpimc_StatScal1, 0, myThid )
         CALL EXCH_3D_RL( hpimc_StatScal1, Nr, myThid )
       ENDIF
       IF ( hpimc_Scal2File .NE. ' ' ) THEN
         CALL READ_FLD_XYZ_RL( hpimc_Scal2File, ' ',
     &                         hpimc_StatScal2, 0, myThid )
         CALL EXCH_3D_RL( hpimc_StatScal2, Nr, myThid )
       ENDIF

       IF ( hpimc_VelUFile .NE. ' ' ) THEN
         CALL READ_FLD_XYZ_RL( hpimc_VelUFile, ' ',
     &                         hpimc_StatVelU, 0, myThid )
       ENDIF
       IF ( hpimc_VelVFile .NE. ' ' ) THEN
         CALL READ_FLD_XYZ_RL( hpimc_VelVFile, ' ',
     &                         hpimc_StatVelV, 0, myThid )
       ENDIF
       IF ( hpimc_VelUFile .NE. ' ' .OR.
     &      hpimc_VelVFile .NE. ' ' ) THEN
        IF ( hpimc_StaV_Cgrid ) THEN
         CALL EXCH_UV_3D_RL( hpimc_StatVelU, hpimc_StatVelV,
     &                       .TRUE., Nr, myThid )
        ELSE
C-      Assume Agrid position:
         CALL EXCH_UV_AGRID_3D_RL( hpimc_StatVelU, hpimc_StatVelV,
     &                       .TRUE., Nr, myThid )
        ENDIF
       ENDIF
#endif /* HPIMC_3D_STATE */

#ifdef HPIMC_2D_STATE
       IF ( hpimc_Surf1File .NE. ' ' ) THEN
         CALL READ_FLD_XY_RL( hpimc_Surf1File,' ',
     &                        hpimc_Surf1, 0, myThid )
         CALL EXCH_XY_RL( hpimc_Surf1, myThid )
       ENDIF
       IF ( hpimc_Surf2File .NE. ' ' ) THEN
         CALL READ_FLD_XY_RL( hpimc_Surf2File,' ',
     &                        hpimc_Surf2, 0, myThid )
         CALL EXCH_XY_RL( hpimc_Surf2, myThid )
       ENDIF
#endif /* HPIMC_2D_STATE */

      ELSE
C--   restart from a pickup:
         CALL HPIMC_READ_PICKUP( nIter0, myThid )

C--   end start-from-iter-zero if/else block
      ENDIF

      RETURN
      END
