#include "SEAFRZ_OPTIONS.h"

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
#include "SEAFRZ.h"
C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     myIter :: iteration counter for this thread
C     myTime :: time counter for this thread
C     myThid :: thread number for this instance of the
C     routine.
      _RL  myTime
      INTEGER myIter
      INTEGER myThid

C     !LOCAL VARIABLES:
C     rho_hpi   :: high-pressure ice density
C     cp_hpi    :: high-pressure ice specific heat
C     Lfus_hpi  :: high-pressure ice enthalpy of fusion
C     pLoc      :: local pressure (dummy variable)
C     sLoc      :: local salinity (dummy variable)
C     tLoc      :: local temperature (dummy variable)
C     rhoLoc    :: local density (dummy variable)
C     ocQflxLB  :: lower boundary heat flux into ocean
C     ocSflxLB  :: lower boundary salt flux into ocean
C     sLB       :: lower boundary-layer salinity
C     tLB       :: lower boundary-layer temperature
C     dbToMPa   :: conversion factor from dbar to MPa (x100)
C     i,j,bi,bj :: loop counters
C     kLB       :: grid cell location of lower boundary

      _RL  rho_hpi
      _RL  cp_hpi
      _RL  Lfus_hpi
      _RL  pLoc,sLoc,tLoc,rhoLoc
      _RL  ocQflxLB, ocSflxLB
      _RL  sLB, tLB
      _RL  dbToMPa
      INTEGER i,j,bi,bj
      INTEGER kLB
      
CEOP
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|


C     rho and cp values taken from Seafreeze
C     heats of fusion estimated using formulae in Dunaeva et al. (2010)
      IF ( sfz_flooriceType .eq. "III" ) THEN
C     assume p = 300 MPa, T = 255 K
        rho_hpi  = 1.165 _d 3
        cp_hpi   = 1.843 _d 3
        Lfus_hpi = 4.541 _d 3
      ELSE IF ( sfz_flooriceType .eq. "V  ") THEN
C     assume p = 450 MPa, T = 265 K ;
C     these are characteristic values for Callisto seafloor
        rho_hpi  = 1.252 _d 3
        cp_hpi   = 1.824 _d 3
        Lfus_hpi = 4.973 _d 3
      ELSE IF ( sfz_flooriceType .eq. "VI ") THEN
C     assume p = 700 MPa, T = 280 K ;
C     these are characteristic values for Ganymede seafloor
        rho_hpi  = 1.332 _d 3
        cp_hpi   = 2.188 _d 3
        Lfus_hpi = 5.192 _d 3
      END IF
      
      IF ( sfz_addFloorice ) THEN      
      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO j = 1-OLy,sNy+OLy
         DO i = 1-OLx,sNx+OLx
C     Initialize forcing arrays with zeros
          sfz_flooriceForcingT(i,j,bi,bj) = 0. _d 0
          sfz_flooriceForcingS(i,j,bi,bj) = 0. _d 0 
C     Find grid cell where boundary layer is located
C     For now, we use the lowest "wet cell" for a given column
          kBL = kLowC(i,j,bi,bj)
C     Find local values of thermodynamic variables (p,T,S,rho)
C     rho & p approximated as reference values at k-level
          rhoLoc = rhoRef(kBL)
          pLoc = pRef4EOS(kBL)*dbToMPa
          sLoc = salt(i,j,kBL,bi,bj)
          tLoc = SW_TEMP(sLoc, theta(i,j,kBL,bi,bj),
     &                   pLoc/dbToMPa, zeroRL)
          IF (sfz_allowFloFrz) THEN
           IF (sfz_fpParamType .eq. 'LINEAR') THEN
            tBL = sfz_FPcoef_l(1) + sfz_FPcoef_l(2)*pLoc 
     &                            + sfz_FPcoef_l(3)*sBL
           ENDIF
          ELSE
           sfz_floorice_dmdt(i,j,bi,bj) = 0. _d 0
           IF (sfz_fpParamType .eq. 'LINEAR') THEN
            tBL = sfz_FPcoef_l(1) + sfz_FPcoef_l(2)*pLoc 
     &                            + sfz_FPcoef_l(3)*sLoc
           ENDIF
           ocQflxLB = gammaTo_hpi*rhoLoc*HeatCapacity_Cp*(tLoc-tBL)
           sfz_floorice_Qflx(i,j,bi,bj) =
     &         rho_hpi*cp_hpi*kappaT_hpi*dTdzBot
     &       + ocQflxLB
           ocSflxLB = 0. _d 0
           sfz_floorice_Sflx(i,j,bi,bj) = 0. _d 0
          ENDIF
C     calculate per-unit-depth T and S forcings; will be converted
C     to actual tendencies in the floorice_forcing routines
          sfz_flooriceForcingT(i,j,bi,bj) =
     &                  ocQflxLB/(HeatCapacity_Cp*rhoLoc)
          sfz_flooriceForcingS(i,j,bi,bj) =
     &                  ocSflxLB/(rhoLoc)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      ENDIF
