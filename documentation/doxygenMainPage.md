CLASSIC main page {#mainpage}
============

# The Canadian Land Surface Scheme Including biogeochemical Cycles (CLASSIC) {#main}

The Canadian Land Surface Scheme Including biogeochemical Cycles (CLASSIC) simulates the exchanges of energy, water, carbon, and momentum at the earth's surface. CLASSIC is formed by the coupling of the Canadian Land Surface Scheme (CLASS) and the Canadian Terrestrial Ecosystem Model (CTEM). CLASS handles the model physics including fluxes of energy, water and momentum. CTEM simulates biogeochemical cycles including fluxes of carbon.

1. @subpage overviewCLASS
   - @subpage devHistory
2. @subpage overviewCTEM
3. @subpage compvsmosaic
4. There are a basic four types of data that are required to run CLASSIC:
  - @subpage forcingData
  - @subpage vegetationData
    - @subpage vegCLASSonly
    - @subpage vegCTEMtoo
  - @subpage soilData
  - @subpage initProgVar
    - @subpage initPhysProgVar
    - @subpage initCTEMProgVar
5. @subpage CTEMaddInputs
  - Greenhouse gases
    - @subpage initCO2
    - @subpage initCH4
  - Disturbrance (fire) inputs
    - @subpage initLightFire
    - @subpage initPopd
  - @subpage initClimComp "Competition for space between PFTs"
  - Prognostic simulation of methane emissions
    - @subpage initRice "Rice agriculture"
    - @subpage initWetSlope "Dynamically-determined wetlands"
    - @subpage initWetArea
  - @subpage initPeat "Peatlands"
  - @subpage inputLUC
6. @subpage makeInputFiles
  - @subpage makeMet "Meteorological inputs"
  - @subpage makeInit "The model initialization and restart files"
  - Greenhouse gas inputs
  - @subpage makeOther "Other inputs"
7. @subpage runPrep "Preparing a CLASSIC run"
  - @subpage compilingMod
  - @subpage setupJobOpts
  - @subpage xmlSystem "Configuring the model outputs"
8. Running CLASSIC
  - @subpage runStandAloneMode "Running CLASSIC for a point location"
  - @subpage runGrid "Running CLASSIC over a grid (regional,global)"
