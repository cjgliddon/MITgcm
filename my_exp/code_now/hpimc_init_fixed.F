#include "HPIMC_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: HPIMC_INIT_FIXED

C !INTERFACE:
      SUBROUTINE HPIMC_INIT_FIXED( myThid )

C     !DESCRIPTION:
C     Initialize HPIMC variables that are kept fixed during the run.

C     !USES:
      IMPLICIT NONE
C     == Global variables ===
#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#include "HPIMC.h"

C     !INPUT/OUTPUT PARAMETERS:
C     myThid ::  my Thread Id number
      INTEGER myThid
CEOP

C     !LOCAL VARIABLES:
C      INTEGER iUnit
      INTEGER iCoef, i, j, bi, bj
      _RL MW_H2O

#ifdef ALLOW_MNC
      IF (useMNC) THEN
        CALL HPIMC_MNC_INIT( myThid )
      ENDIF
#endif /* ALLOW_MNC */

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN
        CALL HPIMC_DIAGNOSTICS_INIT( myThid )
      ENDIF
#endif
      
C     rho, cp, and Lfus values taken from Seafreeze.
C     This updates an earlier version of the code where the Lfus values were
C     estimated using formulas from Dunaeva et al. (2010); now Lfus is just
C     the enthalpy difference between the ice phase and pure water at the
C     specified P/T.
      IF ( hpimc_flooriceType .eq. "III" ) THEN
        sfz_FPcoef_l(1) =  3.4648 _d -2
        sfz_FPcoef_l(2) = -6.1843 _d -1
        sfz_FPcoef_l(3) =  2.4444 _d  2
C     assume p = 300 MPa, T = 255 K
        rho_hpi  = 1.165 _d 3
        cp_hpi   = 1.843 _d 3
        Lfus_hpi = 2.656 _d 5
      ELSE IF ( hpimc_flooriceType .eq. "V  ") THEN
        sfz_FPcoef_l(1) =  5.7906 _d -2
        sfz_FPcoef_l(2) = -8.0362 _d -2
        sfz_FPcoef_l(3) =  2.3781 _d  2
C     assume p = 450 MPa, T = 265 K ;
C     these are characteristic values for Callisto seafloor
        rho_hpi  = 1.252 _d 3
        cp_hpi   = 1.824 _d 3
        Lfus_hpi = 2.741 _d 5
      ELSE IF ( hpimc_flooriceType .eq. "VI ") THEN
        sfz_FPcoef_l(1) =  7.3012 _d -2
        sfz_FPcoef_l(2) = -8.9250 _d -1
        sfz_FPcoef_l(3) =  2.2884 _d  2
C     assume p = 700 MPa, T = 280 K ;
C     these are characteristic values for Ganymede seafloor
        rho_hpi  = 1.332 _d 3
        cp_hpi   = 2.188 _d 3
        Lfus_hpi = 2.836 _d 5
      END IF


C     dummy code for unset coefficients
C      DO iCoef=1,20
C       sfz_EOScoef_w(iCoef) = 0. _d 0
C      ENDDO

C     create bottom heat flux array
      IF ( hpimc_addFloorice ) THEN      
C     instantiate bottom forcing array
       IF (hpimc_QbotType .EQ. 0) THEN
        DO bj=myByLo(myThid), myByHi(myThid)
         DO bi=myBxLo(myThid), myBxHi(myThid)
          DO j = 1-OLy,sNy+OLy
           DO i = 1-OLx,sNx+OLx
            hpimc_Qbot(i,j,bi,bj) = -rho_hpi*cp_hpi*kappaT_hpi*dTdz_hpi
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDIF
      ELSE IF ((hpimc_QbotType .EQ. 1) .AND. (hpimc_QbotFile .NE. ' ')) THEN
       CALL READ_FLD_XY_RL( hpimc_QbotFile, ' ', hpimc_Qbot, 0, myThid)
      END IF

C     Read in boundary-layer temperature fields if we're prescribing them
      IF (hpimc_prescribeTb) THEN
       IF (hpimc_TbtopFile .NE. ' ' ) THEN
        CALL READ_FLD_XY_RL( hpimc_TbtopFile, ' ', hpimc_floorice_tBL, 0, myThid)
       ENDIF
       IF (hpimc_TbbotFile .NE. ' ' ) THEN
        CALL READ_FLD_XY_RL( hpimc_TbbotFile, ' ', hpimc_shelfice_tBL, 0, myThid)
       ENDIF
      ENDIF

      RETURN
      END
