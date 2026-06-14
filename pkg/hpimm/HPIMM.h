CBOP
C     !ROUTINE: HPIMM.h
C     !INTERFACE:
C     #include "HPIMM.h"

C     !DESCRIPTION:
C     *================================================================*
C     | HPIMM.h
C     | o Header file defining "hpimm" parameters and variables
C     |   The hpimm (High-Pressure Icy Moon Configuration) package is 
C     |   used for modeling the ocean dynamics of high pressure icy
C     |   (Solar System) satellites.
C     *================================================================*

CEOP

C     Package flag
      LOGICAL hpimm_MNC
      LOGICAL hpimm_MDSIO
      COMMON /HPIMM_PACKAGE/
     &                     hpimm_MNC, hpimm_MDSIO


C-----------------------------------------------------------------------
C     Logical variables:
C     hpimm_addFloorice       :: flag for adding ice to lower ocean boundary
C     hpimm_doAutoTIni        :: flag turning on/off calculation of initial
C                                temperature profile in-model

      LOGICAL hpimm_addFloorice
      LOGICAL hpimm_doAutoTIni

C-----------------------------------------------------------------------
C     Numerical variables:
C     hpimm_QbotType    :: flag for how bottom heating is specified.
C                          Calculated in-model if 0; read from input if 1
C     rho_hpi           :: high-pressure ice density [kg/m^3]
C     Lfus_hpi          :: high-pressure ice enthalpy of fusion [J/kg]
C     gammaT_hpi        :: lower-boundary thermal diffusion coefficient [m/s]
C     gammaS_hpi        :: lower-boundary salt diffusion coefficient [m/s]
C     kappaT_hpi        :: thermal conductivity of high-pressure ice [W/m/K]
C                          (used only if hpimm_QbotType .eq. 0)
C     dTdz_hpi          :: temperature gradient in high-pressure ice [K/m]
C                          (used only if hpimm_QbotType .eq. 0)
C     sfz_EOScoef_w     :: coefficients for the SeaFreeze polynomial EoS
C     sfz_EOSnorm_w     :: norms for the SeaFreeze polynomial EOS
C                          (used in evaluating the fitting function)
C     sfz_FPcoef_l      :: coeffs. for linear freezing point EoS

      INTEGER hpimm_QbotType	
      _RL rho_hpi
      _RL Lfus_hpi
      _RL gammaT_hpi
      _RL gammaS_hpi
      _RL dTdz_hpi
      _RL sfz_EOScoef_w(0:3,0:3,0:3)
      _RL sfz_EOSnorm_w(6)
      _RL sfz_FPcoef_l(3)

C-----------------------------------------------------------------------
C     Char/string variables:
C     hpimm_flooriceType      :: ice phase at ocean bottom ["III"/"V"/"VI"]
C     hpimm_QbotFile          :: file containing bottom heat flux profile

      CHARACTER*(3) hpimm_flooriceType
      CHARACTER*(MAX_LEN_FNAM) hpimm_QbotFile

C-    file names for initial conditions:
      CHARACTER*(MAX_LEN_FNAM) hpimm_Scal1File
      CHARACTER*(MAX_LEN_FNAM) hpimm_Scal2File
      CHARACTER*(MAX_LEN_FNAM) hpimm_VelUFile
      CHARACTER*(MAX_LEN_FNAM) hpimm_VelVFile
      CHARACTER*(MAX_LEN_FNAM) hpimm_Surf1File
      CHARACTER*(MAX_LEN_FNAM) hpimm_Surf2File

C-----------------------------------------------------------------------

      COMMON /HPIMM_PARAMS_L/
     &       hpimm_addFloorice, hpimm_doAutoTIni
      COMMON /HPIMM_PARAMS_I/
     &       hpimm_QbotType
      COMMON /HPIMM_PARAMS_R/
     &	 rho_hpi, Lfus_hpi,
     &	 gammaT_hpi, gammaS_hpi, dTdz_hpi,
     &	 sfz_EOScoef_w, sfz_EOSnorm_w,
     &	 sfz_FPcoef_l
      COMMON /HPIMM_PARAMS_C/
     &       hpimm_Scal1File, hpimm_Scal2File,
     &       hpimm_VelUFile,  hpimm_VelVFile,
     &       hpimm_Surf1File, hpimm_Surf2File,
     &	 hpimm_flooriceType,
     &	 hpimm_QbotFile     

C-----------------------------------------------------------------------

#ifdef HPIMM_2D_STATE
C     HPIMM 2-dim. field variables
C     hpimm_floorice_tBL      :: lower boundary layer temperature [degC]
C     hpimm_floorice_sBL      :: lower boundary layer salinity    [g/kg]
C     hpimm_floorice_Qflx     :: net upward heat flux at lower boundary [W/m^2]
C     hpimm_floorice_dmdt     :: net freezing rate at lower boundary [kg/m^2/s]
C     hpimm_floorice_Sflx     :: net upward salt flux at l.b. [g/kg/m^2/s]
C     hpimm_flooriceForcingT  :: l.b. temperature forcing [degC/m/s]
C     hpimm_flooriceForcingS  :: l.b. salinity forcing [g/kg/m/s]
C     hpimm_Qbot              :: bottom (internal) heat flux [W/m^2]

      _RL hpimm_floorice_tBL(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL hpimm_floorice_sBL(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL hpimm_floorice_Qflx(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL hpimm_floorice_dmdt(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL hpimm_floorice_Sflx(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL hpimm_flooriceForcingT(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL hpimm_flooriceForcingS(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL hpimm_Qbot(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /SFZ_STATE_2D/
     &    hpimm_floorice_Qflx, hpimm_floorice_dmdt,
     &    hpimm_floorice_Sflx,
     &    hpimm_floorice_tBL, hpimm_floorice_sBL,
     &    hpimm_flooriceForcingT, hpimm_flooriceForcingS,
     &    hpimm_Qbot
#endif /* HPIMM_2D_STATE */

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
