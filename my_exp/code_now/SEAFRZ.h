CBOP
C     !ROUTINE: SEAFRZ.h
C     !INTERFACE:
C     #include "SEAFRZ.h"

C     !DESCRIPTION:
C     *================================================================*
C     | SEAFRZ.h
C     | o Header file defining "seafrz" parameters and variables
C     *================================================================*
CEOP

C     Package flag
      LOGICAL sfz_MNC
      LOGICAL sfz_MDSIO
      COMMON /SFZ_PACKAGE/
     &                     sfz_MNC, sfz_MDSIO

C     SFZ parameters
      LOGICAL sfz_StaV_Cgrid
      LOGICAL sfz_Tend_Cgrid
      LOGICAL sfz_applyTendT
      LOGICAL sfz_applyTendS
      LOGICAL sfz_applyTendU
      LOGICAL sfz_applyTendV

C-    additional parameters:
      LOGICAL sfz_doSwitch1
      LOGICAL sfz_doSwitch2
C     specifies whether to use polynomial or spline EoS 
      LOGICAL sfz_simple
C     specifies whether there's ice at lower boundary
      LOGICAL sfz_addFloorice
      LOGICAL sfz_allowFloFrz
C     specifies whether to calculate an initial temperature profile in-model
      LOGICAL sfz_doAutoTIni

      INTEGER sfz_index1
      INTEGER sfz_index2
      _RL sfz_param1
      _RL sfz_param2
C     high-pressure ice parameters
      _RL rho_hpi
      _RL cp_hpi
      _RL kappaT_hpi
      _RL gammaTo_hpi
      _RL gammaSo_hpi
      _RL Lfus_hpi
      _RL dTdzBot
C     coefficients for polynomial equations of state & f.p. calc
      _RL sfz_EOScoef_w(0:3,0:3,0:3)
      _RL sfz_EOSnorm_w(6)
      _RL sfz_FPcoef_l(3)
      _RL sfz_FPcoef_nl(5)
      CHARACTER*(MAX_LEN_FNAM) sfz_string1
      CHARACTER*(MAX_LEN_FNAM) sfz_string2
C     type of ice at seafloor & type of parameterization used for f.p.
      CHARACTER*(3) sfz_flooriceType
      CHARACTER*(MAX_LEN_FNAM) sfz_fpParamType
C     filename for where the water coefficients are stored
      CHARACTER*(MAX_LEN_FNAM) sfz_EOSFile

C-    file names for initial conditions:
      CHARACTER*(MAX_LEN_FNAM) sfz_Scal1File
      CHARACTER*(MAX_LEN_FNAM) sfz_Scal2File
      CHARACTER*(MAX_LEN_FNAM) sfz_VelUFile
      CHARACTER*(MAX_LEN_FNAM) sfz_VelVFile
      CHARACTER*(MAX_LEN_FNAM) sfz_Surf1File
      CHARACTER*(MAX_LEN_FNAM) sfz_Surf2File

      COMMON /SFZ_PARAMS_L/
     &       sfz_StaV_Cgrid, sfz_Tend_Cgrid,
     &       sfz_applyTendT, sfz_applyTendS,
     &       sfz_applyTendU, sfz_applyTendV,
     &       sfz_doSwitch1, sfz_doSwitch2,
     &	     sfz_simple, sfz_addFloorice,
     &       sfz_allowFloFrz, sfz_doAutoTIni
      COMMON /SFZ_PARAMS_I/ sfz_index1, sfz_index2
      COMMON /SFZ_PARAMS_R/
     &	     sfz_param1, sfz_param2,
     &	     kappaT_hpi, gammaTo_hpi, gammaSo_hpi,
     &	     dTdzBot, 
     &	     rho_hpi, cp_hpi, Lfus_hpi,
C	EOS parameters
     &	     sfz_EOScoef_w, sfz_EOSnorm_w,
     &	     sfz_FPcoef_l,  sfz_FPcoef_nl
      COMMON /SFZ_PARAMS_C/ sfz_string1, sfz_string2,
     &       sfz_Scal1File, sfz_Scal2File,
     &       sfz_VelUFile,  sfz_VelVFile,
     &       sfz_Surf1File, sfz_Surf2File,
     &	     sfz_flooriceType, sfz_fpParamType,
     &	     sfz_EOSFile     

#ifdef SEAFRZ_3D_STATE
C     SFZ 3-dim. fields
      _RL sfz_StatScal1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL sfz_StatScal2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL sfz_StatVelU (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL sfz_StatVelV (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      COMMON /SFZ_STATE_3D/
     &    sfz_StatScal1, sfz_StatScal2,
     &    sfz_StatVelU,  sfz_StatVelV
#endif /* SEAFRZ_3D_STATE */
#ifdef SEAFRZ_2D_STATE
C     SFZ 2-dim. fields
      _RL sfz_Surf1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL sfz_Surf2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C     seafloor boundary-layer temp & salt, heat, mass, salt fluxes
      _RL sfz_floorice_tBL(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL sfz_floorice_sBL(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL sfz_floorice_Qflx(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL sfz_floorice_dmdt(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL sfz_floorice_Sflx(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C     *per unit depth* temperature and salt forcings
      _RL sfz_flooriceForcingT(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL sfz_flooriceForcingS(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C     Bottom heat flux; assumed FIXED.
      _RL sfz_Qbot(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /SFZ_STATE_2D/
     &    sfz_Surf1, sfz_Surf2,
     &    sfz_floorice_Qflx, sfz_floorice_dmdt,
     &    sfz_floorice_Sflx,
     &    sfz_floorice_tBL, sfz_floorice_sBL,
     &    sfz_flooriceForcingT, sfz_flooriceForcingS,
     &	  sfz_Qbot
#endif /* SEAFRZ_2D_STATE */

#ifdef SEAFRZ_TENDENCY
      _RL sfz_TendScal1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL sfz_TendScal2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL sfz_TendVelU (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL sfz_TendVelV (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      COMMON /SFZ_TENDENCY/
     &    sfz_TendScal1, mypa_TendScal2,
     &    sfz_TendVelU,  mypa_TendVelV
#endif /* SEAFRZ_TENDENCY */

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
