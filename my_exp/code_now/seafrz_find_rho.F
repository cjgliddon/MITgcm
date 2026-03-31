#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"
#include "SEAFRZ_OPTIONS.h"

C     TODO: clean up this file
CBOP
C     !ROUTINE: DEVALTC
C     ------------------------------------------------------------
      SUBROUTINE DEVALTC(
     I                  p, T, ml, Mp, MT, Mml,
     O                  RES )
C     !USES:
      IMPLICIT NONE
C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "EOS.h"
#include "SEAFRZ.h"

C     !DESCRIPTION:
C     Evaluate a partial derivative of the tricubic polynomial.
C
C     Arguments:
C       p, T, ml      (in) -- evaluation point (original coordinates)
C       Mp, MT, Mml   (in) -- derivative orders w.r.t. x, y, z
C                            (0 = function value, 1 = first deriv, etc.)
C                            Must satisfy 0 <= Mp,MT,Mml <= 3.
C                            Returns 0 if any order exceeds 3.
C
C     Note: LOADTC must have been called first to populate /TRICUB/.
C
C     The chain rule for normalisation gives a factor of (1/H)^M
C     for an M-th derivative along an axis with half-range H.
C     The derivative of x^P is P!/(P-M)! * x^(P-M) for P>=M, else 0.
C     Both factors are folded into the DPX/DPY/DPZ arrays so that
C     the summation loop is identical to the undifferentiated case.
C     ------------------------------------------------------------
      _RL p, T, ml
      INTEGER Mp, MT, Mml
      _RL RES
C     !Local variables
      _RL pC,pH,TC,TH,mlC,mlH
      _RL pN, TN, mlN
      _RL DPp(0:3), DPT(0:3), DPml(0:3)
      _RL SUM, SCLEp, SCLET, SCLEml
      _RL DFACT
      INTEGER i, j, k

C     Return zero immediately for out-of-range derivative orders
      IF (Mp.LT.0 .OR. Mp.GT.3 .OR.
     &    MT.LT.0 .OR. MT.GT.3 .OR.
     &    Mml.LT.0 .OR. Mml.GT.3) THEN
       RES = 0.0D0
       RETURN
      END IF

C     Normalise coordinates
      pN = (p - sfz_EOSnorm_w(1)) / sfz_EOSnorm_w(2)
      TN = (T - sfz_EOSnorm_w(3)) / sfz_EOSnorm_w(4)
      mlN = (ml - sfz_EOSnorm_w(5)) / sfz_EOSnorm_w(6)

C     Chain-rule scale factors: d^M/dx^M = (1/H)^M * d^M/dpN^M
      SCLEp = 1.0D0 / sfz_EOSnorm_w(2)**Mp
      SCLET = 1.0D0 / sfz_EOSnorm_w(4)**MT
      SCLEml = 1.0D0 / sfz_EOSnorm_w(6)**Mml

C     Build derivative power vectors.
C
C     For axis u with normalised coordinate UN and derivative order M:
C       DP(P) = 0                          if P < M
C       DP(P) = [P!/(P-M)!] * UN^(P-M)    if P >= M
C
C      The falling-factorial P!/(P-M)! for P=0..3, M=0..3:
C
C      M=0:  1,  1,   1,    1        (no differentiation)
C      M=1:  0,  1,   2,    3
C      M=2:  0,  0,   2,    6
C      M=3:  0,  0,   0,    6

C         --- DPp ---
          DO 10 i = 0, 3
              IF (i .LT. Mp) THEN
                  DPp(i) = 0.0D0
              ELSE
                  DPp(i) = DFACT(i, Mp) * pN**(i-Mp)
              END IF
   10     CONTINUE

C         --- DPT ---
          DO 20 j = 0, 3
              IF (j .LT. MT) THEN
                  DPT(j) = 0.0D0
              ELSE
                  DPT(j) = DFACT(j, MT) * TN**(j-MT)
              END IF
   20     CONTINUE

C         --- DPml ---
          DO 30 k = 0, 3
              IF (k .LT. Mml) THEN
                  DPml(k) = 0.0D0
              ELSE
                  DPml(k) = DFACT(k, Mml) * mlN**(k-Mml)
              END IF
   30     CONTINUE

C         Sum over all 64 monomials (identical structure to EVALTC)
          SUM = 0.0D0
          DO 60 k = 0, 3
           DO 50 j = 0, 3
            DO 40 i = 0, 3
             SUM = SUM + sfz_EOScoef_w(i,j,k)*DPp(i)*DPT(j)*DPml(k)
   40       CONTINUE
   50      CONTINUE
   60     CONTINUE

          RES = SCLEp * SCLET * SCLEmL * SUM
          RETURN
      END


C     ------------------------------------------------------------
      _RL FUNCTION DFACT(P, M)
C
C     Falling factorial: P! / (P-M)! = P * (P-1) * ... * (P-M+1)
C     Gives the coefficient from differentiating x^P exactly M times.
C     Assumes 0 <= M <= P <= 3.
C     ------------------------------------------------------------
          INTEGER P, M
          INTEGER I
          _RL RESULT

          RESULT = 1.0D0
          DO 10 I = 0, M-1
              RESULT = RESULT * DBLE(P - I)
   10     CONTINUE
          DFACT = RESULT
          RETURN
      END

      SUBROUTINE FIND_RHO_SFZ_2D(
     I              iMin, iMax, jMin, jMax,
     I              locPres, tFld, sFld,
     O              rhoLoc,
     I              myThid )

C     !USES:
      IMPLICIT NONE
C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "EOS.h"
#include "SEAFRZ.h"


C     == Routine arguments ==
      INTEGER iMin,iMax,jMin,jMax
      _RL locPres(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL tFld   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sFld   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL rhoLoc (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER myThid
C     == Functions
      _RL SW_TEMP
      EXTERNAL SW_TEMP
C     == Local variables ==
      INTEGER i,j
      _RL     ct, sa, p, tReal, rhoRecip
      _RL     Patodb,PaToMPa,mwNaCl 

      Patodb  = 1.0 _d -4
      PaToMPa = 1.0 _d -6
      mwNaCl = 1.80153 _d -2

      DO j=jMin,jMax
       DO i=iMin,iMax
C       convert quantities to K, molal, MPa
        ct = tFld(i,j)
C       PRINT *, "Conservative temp: ", ct
        sa = sFld(i,j) 
        p  = locPres(i,j)
        tReal = SW_TEMP(sa, ct, p*PaTodb, zeroRL)
C        PRINT *, "State:", ct,sa,p,tReal
        IF ((tReal > 4.0 _d 2) .OR. (tReal < -1.0 _d 2)) THEN
         PRINT *, "Warning: Temperature going out of range"
         PRINT *, "tReal, i, j, pres:"
         PRINT *, tReal, i, j, p
        ENDIF
        CALL DEVALTC(p, tReal+273.15, sa, 1, 0, 0, rhoRecip)
        rhoLoc(i,j) = 1.0 _d 0 / rhoRecip
       END DO
      END DO

      RETURN

      END SUBROUTINE

      SUBROUTINE FIND_RHO_SFZ_SCALAR(
     I                      ct, sa, p,
     O                      rho )

#include "EEPARAMS.h"

C     Arguments
      _RL ct, sa, p
      _RL rho
      
      _RL SW_TEMP
      EXTERNAL SW_TEMP

C     Local variables
      _RL mwNaCl
      _RL PaTodb,PaToMPa
      _RL tReal,rhoRecip
      
      PaTodb  = 1.0     _d -4
      PaToMPa = 1.0     _d -6
C      dbToMPa = 1.0     _d -2
      mwNaCl  = 1.80153 _d -2
      
      tReal = SW_TEMP(sa, ct, p*PaTodb, zeroRL)
      PRINT *, "FIND_RHO_SFZ_SCALAR data (sa, ct, p, tReal, rho):"
      PRINT *, sa, ct, p, tReal
      CALL DEVALTC(p, tReal+273.15, sa, 1, 0, 0, rhoRecip)
      rho = 1 _d 0 / rhoRecip
      PRINT *, rho

      RETURN

      END SUBROUTINE
