CBOP
C     !ROUTINE: HPIMC.h
C     !INTERFACE:
C     #include "HPIMC.h"

C     !DESCRIPTION:
C     *================================================================*
C     | HPIMC.h
C     | o Header file defining "hpimc" parameters and variables
C     |   The hpimc (High-Pressure Icy Moon Configuration) package is 
C     |   used for modeling the ocean dynamics of high pressure icy
C     |   (Solar System) satellites.
C     *================================================================*

CEOP

C     Package flag
      LOGICAL hpimc_MNC
      LOGICAL hpimc_MDSIO
      COMMON /HPIMC_PACKAGE/
     &                     hpimc_MNC, hpimc_MDSIO


C-----------------------------------------------------------------------
C     Logical variables:
C     hpimc_addFloorice       :: flag for adding ice to lower ocean boundary
C     hpimc_doAutoTIni        :: flag turning on/off calculation of initial
C                                temperature profile in-model
C     hpimc_prescribeTb       :: if True, boundary-layer temperatures are 
C                                prescribed in input files

      LOGICAL hpimc_addFloorice
      LOGICAL hpimc_doAutoTIni
      LOGICAL hpimc_prescribeTb

C-----------------------------------------------------------------------
C     Numerical variables:
C     hpimc_QbotType    :: flag for how bottom heating is specified.
C                          Calculated in-model if 0; read from input if 1
C     rho_hpi           :: high-pressure ice density [kg/m^3]
C     Lfus_hpi          :: high-pressure ice enthalpy of fusion [J/kg]
C     gammaT_hpi        :: lower-boundary thermal diffusion coefficient [m/s]
C     gammaS_hpi        :: lower-boundary salt diffusion coefficient [m/s]
C     cp_hpi            :: heat capacity of high-pressure ice [J/kg]
C                          (used only if hpimc_QbotType .eq. 0)
C     kappaT_hpi        :: thermal conductivity of high-pressure ice [W/m/K]
C                          (used only if hpimc_QbotType .eq. 0)
C     dTdz_hpi          :: temperature gradient in high-pressure ice [K/m]
C                          (used only if hpimc_QbotType .eq. 0)
C     sfz_EOScoef_w     :: coefficients for the SeaFreeze polynomial EoS
C     sfz_EOSnorm_w     :: norms for the SeaFreeze polynomial EOS
C                          (used in evaluating the fitting function)
C     sfz_FPcoef_l      :: coeffs. for linear freezing point EoS
C     sfz_pRef          :: reference pressure used in T--PT conversions [Pa]

      INTEGER hpimc_QbotType	
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
C     hpimc_flooriceType      :: ice phase at ocean bottom ["III"/"V"/"VI"]
C     hpimc_QbotFile          :: file containing bottom heat flux profile
C      The next two files are read if hpimc_prescribeTb is .TRUE.:
C     hpimc_TbtopFile         :: file containing top BL temperature 
C     hpimc_TbbotFile         :: file containing bottom BL temperature

      CHARACTER*(3) hpimc_flooriceType
      CHARACTER*(MAX_LEN_FNAM) hpimc_QbotFile
      CHARACTER*(MAX_LEN_FNAM) hpimc_TbtopFile
      CHARACTER*(MAX_LEN_FNAM) hpimc_TbbotFile

C-    file names for initial conditions:
      CHARACTER*(MAX_LEN_FNAM) hpimc_Scal1File
      CHARACTER*(MAX_LEN_FNAM) hpimc_Scal2File
      CHARACTER*(MAX_LEN_FNAM) hpimc_VelUFile
      CHARACTER*(MAX_LEN_FNAM) hpimc_VelVFile
      CHARACTER*(MAX_LEN_FNAM) hpimc_Surf1File
      CHARACTER*(MAX_LEN_FNAM) hpimc_Surf2File

C-----------------------------------------------------------------------

      COMMON /HPIMC_PARAMS_L/
     &       hpimc_addFloorice, hpimc_doAutoTIni,
     &       hpimc_prescribeTb
      COMMON /HPIMC_PARAMS_I/
     &       hpimc_QbotType
      COMMON /HPIMC_PARAMS_R/
     &	 rho_hpi, Lfus_hpi,
     &	 gammaT_hpi, gammaS_hpi, 
     &       cp_hpi, kappaT_hpi, dTdz_hpi,
     &	 sfz_EOScoef_w, sfz_EOSnorm_w,
     &	 sfz_FPcoef_l, sfz_pRef
      COMMON /HPIMC_PARAMS_C/
     &       hpimc_TbtopFile, hpimc_TbbotFile,
     &       hpimc_Scal1File, hpimc_Scal2File,
     &       hpimc_VelUFile,  hpimc_VelVFile,
     &       hpimc_Surf1File, hpimc_Surf2File,
     &	 hpimc_flooriceType,
     &	 hpimc_QbotFile     

C-----------------------------------------------------------------------

#ifdef HPIMC_2D_STATE
C     HPIMC 2-dim. field variables
C     hpimc_floorice_tBL      :: lower boundary layer temperature [degC]
C     hpimc_floorice_sBL      :: lower boundary layer salinity    [g/kg]
C     hpimc_floorice_Qflx     :: net upward heat flux at lower boundary [W/m^2]
C     hpimc_floorice_dmdt     :: net freezing rate at lower boundary [kg/m^2/s]
C     hpimc_floorice_Sflx     :: net upward salt flux at l.b. [g/kg/m^2/s]
C     hpimc_flooriceForcingT  :: l.b. temperature forcing [degC/m/s]
C     hpimc_flooriceForcingS  :: l.b. salinity forcing [g/kg/m/s]
C     hpimc_Qbot              :: bottom (internal) heat flux [W/m^2]
C     hpimc_shelfice_tBL      :: upper boundary layer temperature [degC]
C                                (used only if hpimc_prescribeTb = .TRUE.)

      _RL hpimc_floorice_tBL(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL hpimc_floorice_sBL(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL hpimc_floorice_Qflx(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL hpimc_floorice_dmdt(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL hpimc_floorice_Sflx(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL hpimc_flooriceForcingT(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL hpimc_flooriceForcingS(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL hpimc_Qbot(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      COMMON /SFZ_STATE_2D/
     &    hpimc_floorice_Qflx, hpimc_floorice_dmdt,
     &    hpimc_floorice_Sflx,
     &    hpimc_floorice_tBL, hpimc_floorice_sBL,
     &    hpimc_flooriceForcingT, hpimc_flooriceForcingS,
     &    hpimc_Qbot, hpimc_shelfice_tBL
#endif /* HPIMC_2D_STATE */

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
