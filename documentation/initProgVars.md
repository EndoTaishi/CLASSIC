# Initialization of Prognostic Variables {#initProgVar}

# Initialization of Physics Prognostic Variables {#initPhysProgVar}

CLASSIC requires initial values of the land surface prognostic variables, either from the most recent atmospheric model integration or from field measurements. These are listed below, with guidelines for specifying values for each if measured values are not available.

- TBARROT Temperature of soil layers [K]
- THICROT Volumetric frozen water content of soil layers \f$[m^3 m^{-3} ]\f$
- THLQROT Volumetric liquid water content of soil layers \f$[m^3 m^{-3} ]\f$

TBARROT, THLQROT and THICROT are required for each of the modelled soil layers. Thin soil layers near the surface equilibrate quickly, but thicker, deeper layers respond more slowly, and long-term biases can be introduced into the simulation if their temperatures and moisture contents are not initialized to reasonable values. If measured values are not available, for the moisture contents, it may be better to err on the low side, since soil moisture recharge typically takes place on shorter time scales than soil moisture loss. Field capacity is commonly used as an initial value. If the soil layer temperature is above freezing, the liquid moisture content would be set to the field capacity and the frozen moisture content to zero; if the layer temperature is below zero, the liquid moisture content would be set to the minimum value and the frozen moisture content to the field capacity minus the minimum value (THLMIN; see @ref CLASSB.f). Very deep soil temperatures do not have a large effect on surface fluxes, but errors in their initial values can adversely affect hydrological simulations. For rock or ice layers, THLQROT and THICROT should both be set to zero.

- TBASROT Temperature of bedrock in the last permeable soil layer [K]

If CLASSIC is run with more than 3 soil layers, a separate calculation of TBASROT is deemed unnecessary (see @ref TMCALC.f). If the model is run with only 3 layers then TBASROT is used. Generally this can be initialized to the temperature of the third soil layer (as is presently done in @ref model_state_drivers.read_initialstate) and is not required to be read in from the initialization file.

- SNOROT Mass of snow pack \f$[kg m^{-2} ]\f$
- TSNOROT Snowpack temperature [K]
- WSNOROT Liquid water content of snow pack \f$[kg m^{-2} ]\f$
- ALBSROT Snow albedo [ ]
- RHOSROT Density of snow \f$[kg m^{-3} ]\f$

It is best to begin a simulation in snow-free conditions, so that the snow simulation can start from the simplest possible state where SNOROT, TSNOROT, ALBSROT, RHOSROT and WSNOROT are all initialized to zero. If erroneous values of the snow variables are specified as initial conditions, this can lead to a persistent bias in the land surface simulation. It can also lead to instability and resulting model crashes. These variables are, however, written to the restart file allowing one to restart the model from the same state (for example, after a spinup to begin a transient simulation). WSNOROT is not written or read-in from the restart file as its influence is small.

- TACROT Temperature of air within vegetation canopy [K]
- TCANROT Vegetation canopy temperature [K]
- RCANROT Intercepted liquid water stored on canopy \f$[kg m^{-2} ]\f$
- SCANROT Intercepted frozen water stored on canopy \f$[kg m^{-2} ]\f$
- QACROT Specific humidity of air within vegetation canopy space \f$[kg kg^{-1} ]\f$
- CMAIROT Aggregated mass of vegetation canopy \f$[kg m^{-2} ]\f$

The vegetation canopy has a relatively small heat capacity and water storage capacity compared with the soil, so its temperature and intercepted water stores equilibrate quite quickly. TCANROT and TACROT can be initialized to the air temperature and QACROT to the air specific humidity. RCANROT and SCANROT can be initialized to zero. CMAIROT, which is used only in the diagnostic energy balance check during the time step, can also be set to zero. At present the initialization file does not contain TACROT, QACROT or CMAIROT. All are initialized in @ref model_state_drivers.read_initialstate

- GROROT Vegetation growth index [ ]

GROROT should be initialized to 1 during the growing season and to 0 otherwise.

- TPNDROT Temperature of ponded water [K]
- TSFSROT Ground surface temperature over subarea [K]
- ZPNDROT Depth of ponded water on surface [m]

Surface ponded water is a small term and is ephemeral in nature, so ZPNDROT and TPNDROT can both be initialized to zero. TSFSROT is included simply to provide a first guess for the surface temperature iteration in the next time step, so it can be initialized to an arbitrary value (presently it is not read in but rather set in @ref model_state_drivers.read_initialstate for the snow-covered subareas of the surface to the freezing point of water; for the snow-free subareas to the temperature of the first soil layer).

# Initialization of Biogeochemical Prognostic Variables {#initCTEMProgVar}

CTEM's initialization variables are generally PFT dependent. Several model options require extra initialization variables which are described in @ref CTEMaddInputs.

- gleafmas Green leaf mass \f$[kg C m^{-2} ]\f$
- bleafmas Brown leaf mass \f$[kg C m^{-2} ]\f$
  - Brown leaf mass is only for grasses, all of other PFTs have a value of zero.
- rootmass Root mass \f$[kg C m^{-2} ]\f$
  - Root mass is assumed to include both the fine and coarse root masses of the plants. No distinction is made.
- stemmass Green Stem mass \f$[kg C m^{-2} ]\f$
  - Grasses and crops have no stem mass so are given values of zero.
- soilcmas Soil C mass per soil layer \f$[kg C m^{-2} ]\f$
- litrmass Litter mass per soil layer \f$[kg C m^{-2} ]\f$
  - Both litrmass and soilcmas are per PFT but also for bare ground (dimension is icc + 1, where icc is the number of CTEM PFTs)

Specify initialization amounts for green leaf biomass, brown leaf biomass, stem biomass, litter and soil carbon mass for the CTEM PFTs. When growing vegetation from bare ground set these to zero or if no values are known.

- lfstatus Leaf pheological state []
  - The phenology subroutine of CTEM tracks leaf status according to four plants states, 1) leaf onset or maximum growth, 2) normal growth, 3) leaf fall, and 4) no leaves state. When leaf status is set to 4 the model thinks there are no leaves. When initializing from bare ground with no vegetation set this to 4.
- pandays Days with positive net photosynthesis [Days]
  - The phenology subroutine of CTEM tracks days with positive net photosynthesis to initiate leaf onset. When net photosynthesis is positive for seven days, favourable weather is assumed to arrive and, leaf onset begins. If this variable is set to 7 then the model will think that favourable weather has arrived. When initializing from bare ground with no vegetation set this variable to 0.
