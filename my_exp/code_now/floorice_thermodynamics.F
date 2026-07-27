#include "HPIMM_OPTIONS.h"

CBOP
C
C
      SUBROUTINE FLOORICE_THERMODYNAMICS(
     I                        myTime, myIter, myThid )
C
C     !USES:
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "HPIMM.h"
C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     myIter :: iteration counter for this thread
C     myTime :: time counter for this thread
C     myThid :: thread number for this instance of the
C     routine.
      _RL  myTime
      INTEGER myIter
      INTEGER myThid

C     !FUNCTIONS

C     !LOCAL VARIABLES:
C     pLoc      :: local pressure (dummy variable)
C     sLoc      :: local salinity (dummy variable)
C     tLoc      :: local temperature (dummy variable)
C     rhoLoc    :: local density (dummy variable)
C     ocQflxLB  :: lower boundary heat flux into ocean
C     ocSflxLB  :: lower boundary salt flux into ocean
C     sBL       :: lower boundary-layer salinity
C     tBL       :: lower boundary-layer temperature
C     dbToMPa   :: conversion factor from dbar to MPa (x100)
C     PaTodb    :: conversion factor from Pa to dbar (x 1/10000)
C     i,j,bi,bj :: loop counters
C     kLB       :: grid cell location of lower boundary

      _RL  pLoc,sLoc,tLoc,rhoLoc
      _RL  ocQflxLB, ocSflxLB
      _RL  sBL, tBL
      _RL  dbToMPa,PaTodb,MPa2Pa
C-    FP coefficient names for more readable code
      _RL  a1,a2,a3,a4,a5
C-    dummy variables for solving the "three-equations" system
      _RL  e1,e2,e3,e4,e5,e6,e7,e8,e9
      _RL  aqe, bqe, cqe
      _RL  root1,root2,sigma
      INTEGER i,j,bi,bj
      INTEGER kLB
C     added just for debugging
      _RL  thetaLoc
      _RL  MW_H2O
      LOGICAL checkBalances

CEOP
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

C     conversion factors & EOS constants
      dbToMPa  = 1.0 _d -2
      PaTodb   = 1.0 _d -4
      MPa2Pa   = 1.0 _d -6
C     define FP coefficient names

      a1 = sfz_FPcoef_l(1)
      a2 = sfz_FPcoef_l(2)
      a3 = sfz_FPcoef_l(3)


C     flag for debugging
      checkBalances = .FALSE.


      IF ( hpimm_addFloorice ) THEN      
       DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
         DO j = 1,sNy
          DO i = 1,sNx
C     Initialize forcing arrays with zeros
           hpimm_flooriceForcingT(i,j,bi,bj) = 0. _d 0
           hpimm_flooriceForcingS(i,j,bi,bj) = 0. _d 0 
C     Find grid cell where boundary layer is located
C     For now, we use the lowest "wet cell" for a given column
           kLB = kLowC(i,j,bi,bj) 
C     Find local values of thermodynamic variables (p,T,S,rho)
C     rho & p approximated as reference values at k-level
           IF ((kLB > 0) .AND. (kLB <= Nr)) THEN
            rhoLoc = rhoRef(kLB)
C     the next few lines calculate a pressure adjustment due to
C     varying bottom topography
            pLoc = (phiRef(2*kLB-1) 
     &                + gravity*gravFacC(kLB)*drF(kLB)
     &                * hFacC(i,j,kLB,bi,bj))
     &                * rhoLoc*PaTodb*dbToMPa
            sLoc = salt(i,j,kLB,bi,bj)
            thetaLoc = theta(i,j,kLB,bi,bj)
            IF ( (pLoc > 1.0 _d 17) .OR. (sLoc > 1.0 _d 10)
     &                  .OR.(thetaLoc > 1.0 _d 4) ) THEN
             print *, "pLoc, sLoc, thetaLoc", pLoc, sLoc, thetaLoc
             print *, "at point:", i, j, kLB, bi, bj
            END IF
            CALL SFZ_THETA2T(pLoc*MPa2Pa, theta(i,j,kLB,bi,bj)+273.15,
     &                       sLoc, sfz_pRef, tLoc)
C           PRINT *, "Theta: ", theta(i,j,kLB,bi,bj) + 273.15
C           PRINT *, "T_abs: ", tLoc
C           PRINT *, "Salin: ", sLoc
C           PRINT *, "prsdb: ", pLoc/dbToMPa
C           PRINT *, "Local temperature at i,j=", i,j
C           PRINT *, tLoc, "Kelvin"
C           PRINT *, "Local theta: ", theta(i,j,kLB,bi,bj)
C           PRINT *, "sLoc: ", sLoc
C           PRINT *, "pLoc (MPa): ", pLoc
            IF (hpimm_prescribeTb) THEN
             ocQflxLB = rhoLoc*gammaT_hpi*HeatCapacity_Cp* 
     &        (hpimm_floorice_tBL(i,j,bi,bj) - thetaLoc) 
             sigma = (ocQflxLb - hpimm_Qbot(i,j,bi,bj))/Lfus_hpi
             sBL = rhoLoc*gammaS_hpi*sLoc/(rhoLoc*gammaS_hpi - sigma)
             ocSflxLB = sBL*sigma
C           Diagnostics
             hpimm_floorice_dmdt(i,j,bi,bj) = sigma
             hpimm_floorice_Sflx(i,j,bi,bj) = -sigma*sBL
             hpimm_floorice_Qflx(i,j,bi,bj) = ocQflxLB
            ELSE
C            three-equations system
             e1 = hpimm_Qbot(i,j,bi,bj)
             e2 = rhoLoc*gammaT_hpi*HeatCapacity_Cp
             e3 = -rhoLoc*gammaS_hpi*Lfus_hpi
             e4 = e2*(tLoc-a1*pLoc-a3)
             aqe = e2*a2
             bqe = e3 - e1 - e4
             cqe = -sLoc*e3
             IF (aqe .NE. 0. _d 0) THEN
C              sBL = (-1*bqe - SQRT(bqe*bqe -4*aqe*cqe))/(2*aqe)
              sBL = -2*cqe/(bqe - SQRT(bqe*bqe - 4*aqe*cqe))
             ELSE
              sBL = -1*cqe/bqe
             ENDIF
             tBL = a1*pLoc + a2*sBL + a3
             IF (sBL .NE. 0. _d 0) THEN
              sigma = gammaS_hpi*rhoLoc*(1 - sLoc/sBL)
             ELSE
              sigma = 0. _d 0
             ENDIF

C             hpimm_floorice_Qflx(i,j,bi,bj) = sigma
C     &            * (HeatCapacity_Cp*(tLoc-tBL) - Lfus_hpi ) 
             ocSflxLB = sigma*sBL
             ocQflxLB = e2*(tBL-tLoc)
             IF (checkBalances) THEN
              PRINT *, "Checking balances for floorice at:",i,j,bi,bj
              PRINT *, "Lfus, e1-e4, aqe,bqe,cqe,discriminant"
              PRINT *, Lfus_hpi,e1,e2,e3,e4,aqe,bqe,cqe
              PRINT *, bqe*bqe - 4*aqe*cqe
              PRINT *, "Conduction term: "
              PRINT *, e1
              PRINT *, "Diffusive term: "
              PRINT *, ocQflxLB
              PRINT *, "Latent heating of BL: "
              PRINT *, Lfus_hpi*sigma
              PRINT *, "Diffusive salt flux:"
              PRINT *, rhoLoc*gammaS_hpi*(sLoc - sBL)
              PRINT *, "Boundary layer state:"
              PRINT *, pLoc, tBL, sBL
              PRINT *, "Ocean state:"
              PRINT *, pLoc, sLoc, tLoc, rhoLoc
             ENDIF 
             hpimm_floorice_dmdt(i,j,bi,bj) = sigma
             hpimm_floorice_Sflx(i,j,bi,bj) = -sigma*sBL
             hpimm_floorice_Qflx(i,j,bi,bj) = ocQflxLB
            ENDIF
C            hpimm_floorice_dmdt(i,j,bi,bj) = 0. _d 0
C            IF (sfz_fpParamType .eq. 'LINEAR') THEN
C             tBL = sfz_FPcoef_l(1) + sfz_FPcoef_l(2)*pLoc 
C     &           + sfz_FPcoef_l(3)*sLoc
C             PRINT *, "Boundary layer temperature: "
C             PRINT *, tBL, "K"
C            ENDIF
C-      Heat flux into BL from ocean
C            ocQflxLB = gammaTo_hpi*rhoLoc*HeatCapacity_Cp*(tBL-tLoc)
C            hpimm_floorice_Qflx(i,j,bi,bj) =
C     &        - rho_hpi*cp_hpi*kappaT_hpi*dTdzBot
C     &        + ocQflxLB
C            ocSflxLB = 0. _d 0
C            hpimm_floorice_Sflx(i,j,bi,bj) = 0. _d 0
                    
C     calculate per-unit-depth T and S forcings; will be converted
C     to actual tendencies in the floorice_forcing routines
           hpimm_flooriceForcingT(i,j,bi,bj) =
     &                   ocQflxLB/(HeatCapacity_Cp*rhoLoc)
           hpimm_flooriceForcingS(i,j,bi,bj) =
     &                   ocSflxLB/(rhoLoc)
           ENDIF
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDIF

      END
