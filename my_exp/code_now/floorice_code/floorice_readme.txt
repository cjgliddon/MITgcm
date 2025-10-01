NOTES ON FLOORICE CODE (Coleman Gliddon, 2025)

The new MITgcm package/module FLOORICE enables the modeling of ice layers at the lower boundary of ocean basins, such as those which exist on the floors of the subsurface oceans of certain icy Solar System moons (e.g. Ganymede, Titan). Its overall structure parallels that of the SHELFICE package, which models surface ice layers.

Basic effects

The basic effects of adding ice to the lower interface of an ocean model are the addition of temperature and salinity forcing terms to the dynamics which arise from the freezing and melting of ice. In MITgcm, these forcings are added in the subroutines APPLY_FORCING_T and APPLY_FORCING_S (defined in the file "apply_forcing.F"). In each of these subroutines, temperature and salinity tendency arrays are successively augmented by contributions from different forcers, which are applied in function calls.

Within APPLY_FORCING_T, the "floor ice tempreature forcing" from freeze/melt is applied by calling:

>            CALL FLOORICE_FORCING_T(
>     U                   gT_arr,
>     I                   iMin,iMax,jMin,jMax, k, bi,bj,
>     I                   myTime, myIter, myThid )

A similar call is made for the salinity forcing (FLOORICE_FORCING_S).

The subroutines FLOORICE_FORCING_T and FLOORICE_FORCING_S do not themselves calculate the forcing; this is done in a subroutine called FLOORICE_THERMODYNAMICS, which produces the arrays flooriceForcingT and flooriceForcingS (note the lack of underscores). These arrays are four-dimensional, with size (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy); put simply, the first two axes correspond to the lon-lat coordinates of the grid cell, while the last two axes identify the processor being used for the calculation (for parallel-processed model runs). The flooriceForcing[T/S] arrays have dimensions of "temperature tendency-distance" (K*m/sec); the FLOORICE_FORCING_[T/S] subroutines convert these values into time tendencies with dimensions of (K/sec). Additionally, if the lowermost wet grid cell contains only a small amount of water, these subroutines "redistribute" the forcing across the lowermost two wet cells to avoid creating unphysically large temperature tendencies (gradients) at the bottom cell. 

>>>     *QUESTION*: why do the do-loops in FLOORICE_FORCING_[S/T] iterate from 1 to sNx/sNy?

In order to properly calculate these quantities, the floorice_thermodynamics module must first be initialized with the floorice_init_fixed and floorice_init_var subroutines. The first of these is a fairly simple subroutine which sets the values of constant parameters; the second (which currently does not exist) initializes variable quantities. 