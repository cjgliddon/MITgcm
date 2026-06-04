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
C     hpimm_prescribeTb       :: if True, boundary-layer temperatures are 
C                                prescribed in input files

      LOGICAL hpimm_addFloorice
      LOGICAL hpimm_doAutoTIni
      LOGICAL hpimm_prescribeTb

C-----------------------------------------------------------------------
C     Numerical variables:
C     hpimm_QbotType    :: flag for how bottom heating is specified.
C                          Calculated in-model if 0; read from input if 1
C     rho_hpi           :: high-pressure ice density [kg/m^3]
C     Lfus_hpi          :: high-pressure ice enthalpy of fusion [J/kg]
C     gammaT_hpi        :: lower-boundary thermal diffusion coefficient [m/s]
C     gammaS_hpi        :: lower-boundary salt diffusion coefficient [m/s]
C     cp_hpi            :: heat capacity of high-pressure ice [J/kg]
C                          (used only if hpimm_QbotType .eq. 0)
C     kappaT_hpi        :: thermal conductivity of high-pressure ice [W/m/K]
C                          (used only if hpimm_QbotType .eq. 0)
C     dTdz_hpi          :: temperature gradient in high-pressure ice [K/m]
C                          (used only if hpimm_QbotType .eq. 0)
C     sfz_EOScoef_w     :: coefficients for the SeaFreeze polynomial EoS
C     sfz_EOSnorm_w     :: norms for the SeaFreeze polynomial EOS
C                          (used in evaluating the fitting function)
C     sfz_FPcoef_l      :: coeffs. for linear freezing point EoS
C     sfz_pRef          :: reference pressure used in T--PT conversions [Pa]

      INTEGER hpimm_QbotType	
      _RL rho_hpi
      _RL Lfus_hpi
      _RL gammaT_hpi
      _RL gammaS_hpi
      _RL cp_hpi
      _RL kappaT_hpi
      _RL dTdz_hpi
      _RL sfz_EOScoef_w(0:3,0:3,0:3)
      _RL sfz_EOSnorm_w(6)
      _RL sfz_FPcoef_l(3)
      _RL sfz_pRef

C-----------------------------------------------------------------------
C     Char/string variables:
C     hpimm_flooriceType      :: ice phase at ocean bottom ["III"/"V"/"VI"]
C     hpimm_QbotFile          :: file containing bottom heat flux profile
C      The next two files are read if hpimm_prescribeTb is .TRUE.:
C     hpimm_TbtopFile         :: file containing top BL temperature 
C     hpimm_TbbotFile         :: file containing bottom BL temperature

      CHARACTER*(3) hpimm_flooriceType
      CHARACTER*(MAX_LEN_FNAM) hpimm_QbotFile
      CHARACTER*(MAX_LEN_FNAM) hpimm_TbtopFile
      CHARACTER*(MAX_LEN_FNAM) hpimm_TbbotFile

C-    file names for initial conditions:
      CHARACTER*(MAX_LEN_FNAM) hpimm_Scal1File
      CHARACTER*(MAX_LEN_FNAM) hpimm_Scal2File
      CHARACTER*(MAX_LEN_FNAM) hpimm_VelUFile
      CHARACTER*(MAX_LEN_FNAM) hpimm_VelVFile
      CHARACTER*(MAX_LEN_FNAM) hpimm_Surf1File
      CHARACTER*(MAX_LEN_FNAM) hpimm_Surf2File

C-----------------------------------------------------------------------

      COMMON /HPIMM_PARAMS_L/
     &       hpimm_addFloorice, hpimm_doAutoTIni,
     &       hpimm_prescribeTb
      COMMON /HPIMM_PARAMS_I/
     &       hpimm_QbotType
      COMMON /HPIMM_PARAMS_R/
     &	 rho_hpi, Lfus_hpi,
     &	 gammaT_hpi, gammaS_hpi, 
     &       cp_hpi, kappaT_hpi, dTdz_hpi,
     &	 sfz_EOScoef_w, sfz_EOSnorm_w,
     &	 sfz_FPcoef_l, sfz_pRef
      COMMON /HPIMM_PARAMS_C/
     &       hpimm_TbtopFile, hpimm_TbbotFile,
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
C     hpimm_shelfice_tBL      :: upper boundary layer temperature [degC]
C                                (used only if hpimm_prescribeTb = .TRUE.)

      _RL hpimm_floorice_tBL(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL hpimm_floorice_sBL(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL hpimm_floorice_Qflx(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL hpimm_floorice_dmdt(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL hpimm_floorice_Sflx(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL hpimm_flooriceForcingT(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL hpimm_flooriceForcingS(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL hpimm_Qbot(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL hpimm_shelfice_tBL(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      COMMON /SFZ_STATE_2D/
     &    hpimm_floorice_Qflx, hpimm_floorice_dmdt,
     &    hpimm_floorice_Sflx,
     &    hpimm_floorice_tBL, hpimm_floorice_sBL,
     &    hpimm_flooriceForcingT, hpimm_flooriceForcingS,
     &    hpimm_Qbot, hpimm_shelfice_tBL
#endif /* HPIMM_2D_STATE */

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
