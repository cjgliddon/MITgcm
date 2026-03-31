#include "SEAFRZ_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: SEAFRZ_INIT_FIXED

C !INTERFACE:
      SUBROUTINE SEAFRZ_INIT_FIXED( myThid )

C     !DESCRIPTION:
C     Initialize SEAFRZ variables that are kept fixed during the run.

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

C     !LOCAL VARIABLES:
C      INTEGER iUnit
      INTEGER iCoef, i, j, bi, bj
      _RL MW_H2O

#ifdef ALLOW_MNC
      IF (useMNC) THEN
        CALL SEAFRZ_MNC_INIT( myThid )
      ENDIF
#endif /* ALLOW_MNC */

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN
        CALL SEAFRZ_DIAGNOSTICS_INIT( myThid )
      ENDIF
#endif
      
C      iUnit=100
C      IF (eosType .eq. "SEAFRZ") THEN
C     Read in coefficients from the data file
C        OPEN(UNIT=iUnit, FILE=sfz_EOSFile, STATUS='OLD',
C     &       FORM='unformatted', ACCESS='sequential')
C        READ(iUnit) sfz_EOSnorm_w, sfz_EOScoef_w
C        CLOSE(iUnit)
C        print *, "Read in coefficients successfully"
C        print *, sfz_EOSnorm_w
C      END IF


      MW_H2O = 0.0180153
C     rho and cp values taken from Seafreeze
C     heats of fusion estimated using formulae in Dunaeva et al. (2010)
      IF ( sfz_flooriceType .eq. "III" ) THEN
        sfz_FPcoef_l(3) =  2.4577 _d  2
        sfz_FPcoef_l(1) =  3.3174 _d -2
        sfz_FPcoef_l(2) = -2.7079 _d -1
C     assume p = 300 MPa, T = 255 K
        rho_hpi  = 1.165 _d 3
        cp_hpi   = 1.843 _d 3
        Lfus_hpi = 4.541 _d 3 / MW_H2O
      ELSE IF ( sfz_flooriceType .eq. "V  ") THEN
        sfz_FPcoef_l(3) =  2.4091 _d  2
        sfz_FPcoef_l(1) =  5.5101 _d -2
        sfz_FPcoef_l(2) = -3.4015 _d -1
C     assume p = 450 MPa, T = 265 K ;
C     these are characteristic values for Callisto seafloor
        rho_hpi  = 1.252 _d 3
        cp_hpi   = 1.824 _d 3
        Lfus_hpi = 4.973 _d 3 / MW_H2O
      ELSE IF ( sfz_flooriceType .eq. "VI ") THEN
        sfz_FPcoef_l(3) =  2.3398 _d  2
        sfz_FPcoef_l(1) =  6.8995 _d -2
        sfz_FPcoef_l(2) = -3.6907 _d -1
C     assume p = 700 MPa, T = 280 K ;
C     these are characteristic values for Ganymede seafloor
        rho_hpi  = 1.332 _d 3
        cp_hpi   = 2.188 _d 3
        Lfus_hpi = 5.192 _d 3 / MW_H2O
      END IF


C     dummy code for unset coefficients
C      DO iCoef=1,20
C       sfz_EOScoef_w(iCoef) = 0. _d 0
C      ENDDO
      DO iCoef=1,5
       sfz_FPcoef_nl(iCoef) = 0. _d 0
      ENDDO
      
C     create bottom heat flux array
      IF ( sfz_addFloorice ) THEN      
C     instantiate bottom forcing array
       DO bj=myByLo(myThid), myByHi(myThid)
        DO bi=myBxLo(myThid), myBxHi(myThid)
         DO j = 1-OLy,sNy+OLy
          DO i = 1-OLx,sNx+OLx
           sfz_Qbot(i,j,bi,bj) = -rho_hpi*cp_hpi*kappaT_hpi*dTdzBot
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDIF

      RETURN
      END
