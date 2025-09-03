#ifdef ALLOW_FLOORICE

CBOP
C !ROUTINE: FLOORICE.h

C !DESCRIPTION: \bv
C     *==========================================================*
C     | FLOORICE.h
C     | o Basic header thermodnynamic floor ice package.
C     |   Contains all FLOORICE field declarations.
C     *==========================================================*

C-----------------------------------------------------------------------
C
C--   Constants that can be set in data.floorice
C     FLOORICEtopoFile         :: File containing the topography of the
C                                 floorice draught (unit=m)
C     FLOORICEmassFile         :: name of floorice Mass file
C     FLOORICEloadAnomalyFile  :: name of floorice load anomaly file
C     FLOORICEMassDynTendFile  :: file name for other mass tendency
C                                 (e.g. dynamics)
C     useISOMIPTD              :: use simple ISOMIP thermodynamics, def: F
C     FLOORICEconserve         :: use conservative form of H&O-thermodynamics
C                                 following Jenkins et al. (2001, JPO), def: F
C     FLOORICEMassStepping     :: flag to step forward ice floor mass/thickness
C                                 accounts for melting/freezing & dynamics
C                                 (from file or from coupling), def: F
C     FLOORICEDynMassOnly      :: step ice mass ONLY with Flooricemassdyntendency
C                                 (not melting/freezing) def: F
C     FLOORICEboundaryLayer    :: turn on vertical merging of cells to for a
C                                 boundary layer of drF thickness, def: F
C     FLI_withBL_realFWflux    :: with above BL, allow to use real-FW flux (and
C                                 adjust advective flux at boundary accordingly)
C                                 def: F
C     FLI_withBL_uStarTopDz    :: with FLOORICEboundaryLayer, compute uStar from
C                                 uVel,vVel avergaged over top Dz thickness;
C                                 def: F
C     FLOORICEadvDiffHeatFlux  :: use advective-diffusive heat flux into the
C                                 floor instead of default diffusive heat
C                                 flux, see Holland and Jenkins (1999),
C                                 eq.21,22,26,31; def: F
C     FLOORICEsaltToHeatRatio  :: constant ratio giving
C                                 FLOORICEsaltTransCoeff/FLOORICEheatTransCoeff
C                                 (def: 5.05e-3)
C     FLOORICEheatTransCoeff   :: constant heat transfer coefficient that
C                                 determines heat flux into floorice
C                                 (def: 1e-4 m/s)
C     FLOORICEsaltTransCoeff   :: constant salinity transfer coefficient that
C                                 determines salt flux into floorice
C                                 (def: FLOORICEsaltToHeatRatio * FLOORICEheatTransCoeff)
C     -----------------------------------------------------------------------
C     FLOORICEuseGammaFrict    :: use velocity dependent exchange coefficients,
C                                 see Holland and Jenkins (1999), eq.11-18,
C                                 with the following parameters (def: F):
C     FLOORICE_oldCalcUStar    :: use old uStar averaging expression
C     fliCdrag                 :: quadratic drag coefficient to compute uStar
C                                 (def: 0.0015)
C     fliZetaN                 :: ??? (def: 0.052)
C     fliRc                    :: ??? (not used, def: 0.2)
C     fliPrandtl, fliSchmidt   :: constant Prandtl (13.8) and Schmidt (2432.0)
C                                 numbers used to compute gammaTurb
C     fliKinVisc               :: constant kinetic viscosity used to compute
C                                 gammaTurb (def: 1.95e-5)
C     FLI_update_kBotC         :: update lateral extension (kTopC) if ice-floor
C                                 retreats from or expands to model top level
C                                 (requires to define ALLOW_FLOORICE_REMEFLING)
C     FLOORICEremeshFrequency  :: Frequency (in seconds) of call to
C                                 FLOORICE_REMEFLING (def: 0. --> no remefling)
C     FLOORICEsplitThreshold   :: Thickness fraction remefling threshold above
C                                  which top-cell splits (no unit)
C     FLOORICEmergeThreshold   :: Thickness fraction remefling threshold below
C                                  which top-cell merges with below (no unit)
C     -----------------------------------------------------------------------
C     FLOORICEDragLinear       :: linear drag at bottom floorice (1/s)
C     FLOORICEDragQuadratic    :: quadratic drag at bottom floorice (default
C                                 = fliCdrag or bottomDragQuadratic)
C     no_slip_floorice         :: set slip conditions for floorice separately,
C                                 (by default the same as no_slip_bottom, but
C                                 really should be false when there is linear
C                                 or quadratic drag)
C     FLOORICElatentHeat       :: latent heat of fusion (def: 334000 J/kg)
C     FLOORICEwriteState       :: enable output
C     FLOORICEHeatCapacity_Cp  :: heat capacity of ice floor (def: 2000 J/K/kg)
C     rhoFloorIce              :: density of ice floor (def: 917.0 kg/m^3)
C
C     FLOORICE_dump_mnc        :: use netcdf for snapshot output
C     FLOORICE_tave_mnc        :: use netcdf for time-averaged output
C     FLOORICE_dumpFreq        :: analoguous to dumpFreq (= default)
C     FLOORICE_taveFreq        :: analoguous to taveFreq (= default)
C
C--   Fields
C     kBotC                  :: index of the top "wet cell" (2D)
C     R_floorIce             :: floorice topography [m]
C     flooriceMassInit       :: ice-floor mass (per unit area) (kg/m^2)
C     flooriceMass           :: ice-floor mass (per unit area) (kg/m^2)
C     floorIceMassDynTendency :: other mass balance tendency  (kg/m^2/s)
C                            ::  (e.g., from dynamics)
C     flooriceLoadAnomaly    :: pressure load anomaly of floorice (Pa)
C     flooriceHeatFlux       :: upward heat flux (W/m^2)
C     flooriceFreshWaterFlux :: upward fresh water flux (virt. salt flux)
C                               (kg/m^2/s)
C     flooriceForcingT       :: analogue of surfaceForcingT
C                               units are  r_unit.Kelvin/s (=Kelvin.m/s if r=z)
C     flooriceForcingS       :: analogue of surfaceForcingS
C                               units are  r_unit.g/kg/s (=g/kg.m/s if r=z)
#ifdef ALLOW_DIAGNOSTICS
C     flooriceDragU          :: Ice-Floor stress (for diagnostics), Zonal comp.
C                               Units are N/m^2 ;   > 0 increase top uVel
C     flooriceDragV          :: Ice-Floor stress (for diagnostics), Merid. comp.
C                               Units are N/m^2 ;   > 0 increase top vVel
#endif /* ALLOW_DIAGNOSTICS */
#ifdef ALLOW_CTRL
C   maskFLI           ::  Mask=1 where ice floor is present on surface
C                           layer, showing full 2D ice floor extent.
C                           =maskC for rest of k values
C                           Used with ice floor fwflx
C                           or fliTransCoeffT/S ctrl.
#endif
C-----------------------------------------------------------------------
C \ev
CEOP

      COMMON /FLOORICE_PARMS_L/
     &     FLOORICEisOn,
     &     useISOMIPTD,
     &     FLOORICEconserve,
     &     FLOORICEboundaryLayer,
     &     FLI_withBL_realFWflux,
     &     FLI_withBL_uStarTopDz,
     &     no_slip_floorice,
     &     FLOORICEwriteState,
     &     FLOORICE_dump_mdsio,
     &     FLOORICE_tave_mdsio,
     &     FLOORICE_dump_mnc,
     &     FLOORICE_tave_mnc,
     &     FLOORICEadvDiffHeatFlux,
     &     FLOORICEuseGammaFrict,
     &     FLOORICE_oldCalcUStar,
     &     FLOORICEMassStepping,
     &     FLOORICEDynMassOnly,
     &     FLI_update_kTopC
      LOGICAL FLOORICEisOn
      LOGICAL useISOMIPTD
      LOGICAL FLOORICEconserve
      LOGICAL FLOORICEboundaryLayer
      LOGICAL FLI_withBL_realFWflux
      LOGICAL FLI_withBL_uStarTopDz
      LOGICAL no_slip_floorice
      LOGICAL FLOORICEwriteState
      LOGICAL FLOORICE_dump_mdsio
      LOGICAL FLOORICE_tave_mdsio
      LOGICAL FLOORICE_dump_mnc
      LOGICAL FLOORICE_tave_mnc
      LOGICAL FLOORICEadvDiffHeatFlux
      LOGICAL FLOORICEuseGammaFrict
      LOGICAL FLOORICE_oldCalcUStar
      LOGICAL FLOORICEMassStepping
      LOGICAL FLOORICEDynMassOnly
      LOGICAL FLI_update_kTopC

      COMMON /FLOORICE_PARMS_I/
     &     FLOORICEselectDragQuadr
      INTEGER FLOORICEselectDragQuadr

      COMMON /FLOORICE_PARMS_R/
     &     FLOORICE_dumpFreq, FLOORICE_taveFreq,
     &     FLOORICEsaltToHeatRatio,
     &     FLOORICEheatTransCoeff, FLOORICEsaltTransCoeff,
     &     rhoFloorice, FLOORICEkappa,
     &     FLOORICElatentHeat,
     &     FLOORICEheatCapacity_Cp,
     &     FLOORICEthetaSurface,
     &     FLOORICEsalinity,
     &     FLOORICEDragLinear, FLOORICEDragQuadratic,
     &     fliCdrag, fliZetaN, fliRc,
     &     fliPrandtl, fliSchmidt, fliKinVisc,
     &     FLOORICEremeshFrequency,
     &     FLOORICEsplitThreshold, FLOORICEmergeThreshold

      _RL FLOORICE_dumpFreq, FLOORICE_taveFreq
      _RL FLOORICEsaltToHeatRatio
      _RL FLOORICEheatTransCoeff
      _RL FLOORICEsaltTransCoeff
      _RL FLOORICElatentHeat
      _RL FLOORICEheatCapacity_Cp
      _RL rhoFloorice
      _RL FLOORICEkappa
      _RL FLOORICEDragLinear
      _RL FLOORICEDragQuadratic
      _RL FLOORICEthetaSurface
      _RL fliCdrag, fliZetaN, fliRc
      _RL fliPrandtl, fliSchmidt, fliKinVisc
      _RL FLOORICEremeshFrequency
      _RL FLOORICEsplitThreshold
      _RL FLOORICEmergeThreshold
      _RL FLOORICEsalinity

      COMMON /FLOORICE_PARM_C/
     &     FLOORICEloadAnomalyFile,
     &     FLOORICEmassFile,
     &     FLOORICEtopoFile,
     &     FLOORICEMassDynTendFile,
     &     FLOORICETransCoeffTFile
      CHARACTER*(MAX_LEN_FNAM) FLOORICEloadAnomalyFile
      CHARACTER*(MAX_LEN_FNAM) FLOORICEmassFile
      CHARACTER*(MAX_LEN_FNAM) FLOORICEtopoFile
      CHARACTER*(MAX_LEN_FNAM) FLOORICEMassDynTendFile
      CHARACTER*(MAX_LEN_FNAM) FLOORICETransCoeffTFile

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      COMMON /FLOORICE_FIELDS_I/ kTopC
      INTEGER kTopC (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      COMMON /FLOORICE_FIELDS_RL/
     &     flooriceMass, flooriceMassInit,
     &     flooriceLoadAnomaly,
     &     flooriceForcingT, flooriceForcingS,
     &     fliCDragFld, fliDragQuadFld

      _RL flooriceMass          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL flooriceMassInit      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL flooriceLoadAnomaly   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL flooriceForcingT      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL flooriceForcingS      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL fliCDragFld           (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL fliDragQuadFld        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      COMMON /FLOORICE_GAMMA_RL/
     &     fliTransCoeffT, fliTransCoeffS
       _RL fliTransCoeffT       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
       _RL fliTransCoeffS       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      COMMON /FLOORICE_FIELDS_RS/
     &     R_floorIce,
     &     flooriceHeatFlux,
     &     floorIceFreshWaterFlux,
     &     floorIceMassDynTendency
      _RS R_floorIce            (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS flooriceHeatFlux      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS flooriceFreshWaterFlux(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS
     &   floorIceMassDynTendency(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

#ifdef ALLOW_CTRL
      COMMON /FLOORICE_MASKS_CTRL/ maskFLI
      _RS maskFLI  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
#endif /* ALLOW_CTRL */

#ifdef ALLOW_DIAGNOSTICS
      COMMON /FLOORICE_DIAG_DRAG/ flooriceDragU, flooriceDragV
      _RS flooriceDragU(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS flooriceDragV(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
#endif /* ALLOW_DIAGNOSTICS */

#endif /* ALLOW_FLOORICE */
