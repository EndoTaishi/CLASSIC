CLASSIC main page {#mainpage}
============

# The Canadian Land Surface Scheme Including biogeochemical Cycles (CLASSIC)

The Canadian Land Surface Scheme Including biogeochemical Cycles (CLASSIC) simulates the exchanges of energy, water, carbon, and momentum at the earth's surface. CLASSIC is formed by the coupling of the Canadian Land Surface Scheme (CLASS) and the Canadian Terrestrial Ecosystem Model (CTEM). CLASS handles the model physics including fluxes of energy, water and momentum. CTEM simulates biogeochemical cycles including fluxes of carbon.

1. @ref overviewCLASS
   - @subpage devHistory
2. @ref overviewCTEM
3. @ref compvsmosaic
4. There are a basic four types of data that are required to run CLASSIC:
  - @ref forcingData
  - @ref vegetationData
    - @subpage vegCLASSonly
    - @subpage vegCTEMtoo
  - @ref soilData
  - @ref initProgVar
    - @subpage initPhysProgVar
    - @subpage initCTEMProgVar
5. @ref CTEMaddInputs
  - @subpage initClimComp
  - @subpage initLightFire
  - @subpage initRice
  - @subpage initWetSlope
  - @subpage initPeat
6. [Creation and formatting of model input files] (@ref makeInputFile)
  - Meteorological inputs
  - The model initialization and restart files
  - Greenhouse gas inputs
  - Other inputs
7. Running CLASSIC
  - @subpage compilingMod
  - @subpage setupJobOpts
  - @subpage xmlSystem
  - [Running CLASSIC for a single location](@subpage runStandAloneMode)
  - [Running CLASSIC over a grid (regional,global)](@subpage runGrid)
