var index =
[
    [ "The Canadian Land Surface Scheme Including biogeochemical Cycles (CLASSIC)", "index.html#main", null ],
    [ "Overview of the Canadian Land Surface Scheme (CLASS)", "overviewCLASS.html", null ],
    [ "Overview of the Canadian Terrestrial Ecosystem Model (CTEM)", "overviewCTEM.html", [
      [ "Plant Functional Types (PFTs) in CLASSIC", "overviewCTEM.html#PFTsCLASSIC", null ],
      [ "Rate change equations for carbon pools", "overviewCTEM.html#CTEMRateChgEqns", null ]
    ] ],
    [ "Composite vs. mosaic representation", "compvsmosaic.html", null ],
    [ "Basic model inputs to run CLASSIC", "basicInputs.html", [
      [ "Atmospheric Forcing Data", "basicInputs.html#forcingData", null ],
      [ "Input Vegetation Data", "basicInputs.html#vegetationData", [
        [ "Development history of CLASS", "overviewCLASS.html#devHistory", null ],
        [ "Required vegetation data", "basicInputs.html#vegCLASSonly", null ],
        [ "Required vegetation data for a biogeochemical simulation (CLASS+CTEM)", "basicInputs.html#vegCTEMtoo", null ]
      ] ],
      [ "Input Soil Data", "basicInputs.html#soilData", null ],
      [ "Initialization of Prognostic Variables", "basicInputs.html#initProgVar", [
        [ "Initialization of Physics Prognostic Variables", "basicInputs.html#initPhysProgVar", null ],
        [ "Initialization of Biogeochemical Prognostic Variables", "basicInputs.html#initCTEMProgVar", null ]
      ] ],
      [ "Example model setups", "basicInputs.html#exModSets", null ]
    ] ],
    [ "Additional inputs depending on model configuration", "CTEMaddInputs.html", [
      [ "Atmospheric carbon dioxide concentration", "CTEMaddInputs.html#initCO2", null ],
      [ "Atmospheric methane concentration", "CTEMaddInputs.html#initCH4", null ],
      [ "Lightning frequency for fire ignition", "CTEMaddInputs.html#initLightFire", null ],
      [ "Population density for fire ignition/suppresion", "CTEMaddInputs.html#initPopd", null ],
      [ "Climatic variables for PFT competition simulations", "CTEMaddInputs.html#initClimComp", null ],
      [ "Rice agriculature", "CTEMaddInputs.html#initRice", null ],
      [ "Orographic information for dynamic wetland scheme", "CTEMaddInputs.html#initWetSlope", null ],
      [ "Prescribed wetland area", "CTEMaddInputs.html#initWetArea", null ],
      [ "Peatland variables", "CTEMaddInputs.html#initPeat", null ],
      [ "Land use change (LUC)", "CTEMaddInputs.html#inputLUC", null ]
    ] ],
    [ "Creation and formatting of model input files", "makeInputFiles.html", [
      [ "Preparation of files for meteorological inputs", "makeInputFiles.html#makeMet", null ],
      [ "Preparation of the model initialization file", "makeInputFiles.html#makeInit", null ],
      [ "Making input files for other input variables", "makeInputFiles.html#makeOther", null ],
      [ "Some notes on input file format", "makeInputFiles.html#inputFileForm", null ]
    ] ],
    [ "Preparing a CLASSIC run", "runPrep.html", [
      [ "Setting up the runtime environment", "runPrep.html#Environ", null ],
      [ "Running CLASSIC in a Singularity Container", "runPrep.html#Containers", null ],
      [ "Compiling CLASSIC for serial and parallel simulations", "runPrep.html#compilingMod", null ],
      [ "Setting up the joboptions file", "runPrep.html#setupJobOpts", null ]
    ] ],
    [ "Configuring the model outputs via the CLASSIC code and Output Variable Editor (OVE)", "xmlSystem.html", [
      [ "The XML Document Structure", "xmlSystem.html#xmlStruct", null ],
      [ "Editing the CLASSIC code to write the output variable", "xmlSystem.html#writeOutput", null ]
    ] ],
    [ "Running CLASSIC", "runCLASSIC.html", [
      [ "Running CLASSIC for a point location", "runCLASSIC.html#runStandAloneMode", null ],
      [ "Running CLASSIC over a grid (Global or Regional)", "runCLASSIC.html#runGrid", null ]
    ] ],
    [ "How do I ... / Something has gone wrong!", "howDoI.html", [
      [ "Change the number/depth/etc of soil layers?", "howDoI.html#chgSoil", null ],
      [ "My run starts and I get an immediate fail!", "howDoI.html#failStart", null ],
      [ "My run won't compile but I didn't do anything I can think of.", "howDoI.html#makeClean", null ],
      [ "initFileConverter tells me I need a soil colour index", "howDoI.html#soilColourIndex", null ],
      [ "Gotcha for Windows users", "howDoI.html#windowsEndings", null ]
    ] ]
];