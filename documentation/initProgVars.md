# Initialization of Prognostic Variables {#initProgVar}

## Initialization of Physics Prognostic Variables {#initPhysProgVar}

CLASSIC requires initial values of the land surface prognostic variables, either from the most recent atmospheric model integration or from field measurements. These are listed below, with guidelines for specifying values for each if measured values are not available.

- TBARROT Temperature of soil layers [K]
- THICROT Volumetric frozen water content of soil layers \f$[m^3 m^{-3} ]\f$
- THLQROT Volumetric liquid water content of soil layers \f$[m^3 m^{-3} ]\f$

TBARROT, THLQROT and THICROT are required for each of the modelled soil layers. Thin soil layers near the surface equilibrate quickly, but thicker, deeper layers respond more slowly, and long-term biases can be introduced into the simulation if their temperatures and moisture contents are not initialized to reasonable values. If measured values are not available, for the moisture contents, it may be better to err on the low side, since soil moisture recharge typically takes place on shorter time scales than soil moisture loss. Field capacity is commonly used as an initial value. If the soil layer temperature is above freezing, the liquid moisture content would be set to the field capacity and the frozen moisture content to zero; if the layer temperature is below zero, the liquid moisture content would be set to the minimum value and the frozen moisture content to the field capacity minus the minimum value (THLMIN; see src/CLASSB.f). Very deep soil temperatures do not have a large effect on surface fluxes, but errors in their initial values can adversely affect hydrological simulations. For rock or ice layers, THLQROT and THICROT should both be set to zero.

- TBASROT Temperature of bedrock in the last permeable soil layer [K]

TBASROT can be set to the third soil layer temperature

- SNOROT Mass of snow pack \f$[kg m^{-2} ]\f$
- TSNOROT Snowpack temperature [K]
- WSNOROT Liquid water content of snow pack \f$[kg m^{-2} ]\f$ FIXME: Not presently in initfile!
- ALBSROT Snow albedo [ ]
- RHOSROT Density of snow \f$[kg m^{-3} ]\f$

It is best to begin a simulation in snow-free conditions, so that the snow simulation can start from the simplest possible state where SNOROT, TSNOROT, ALBSROT, RHOSROT and WSNOROT are all initialized to zero. If erroneous values of the snow variables are specified as initial conditions, this can lead to a persistent bias in the land surface simulation. It can also lead to instability and resulting model crashes.

- TACROT Temperature of air within vegetation canopy [K] FIXME: Not presently in initfile!
- TCANROT Vegetation canopy temperature [K]
- RCANROT Intercepted liquid water stored on canopy \f$[kg m^{-2} ]\f$
- SCANROT Intercepted frozen water stored on canopy \f$[kg m^{-2} ]\f$
- QACROT Specific humidity of air within vegetation canopy space \f$[kg kg^{-1} ]\f$ FIXME: Not presently in initfile!
- CMAIROT Aggregated mass of vegetation canopy \f$[kg m^{-2} ]\f$ FIXME: Not presently in initfile!

The vegetation canopy has a relatively small heat capacity and water storage capacity compared with the soil, so its temperature and intercepted water stores equilibrate quite quickly. TCANROT and TACROT can be initialized to the air temperature and QACROT to the air specific humidity. RCANROT and SCANROT can be initialized to zero. CMAIROT, which is used only in the diagnostic energy balance check during the time step, can also be set to zero.

- GROROT Vegetation growth index [ ]

GROROT should be initialized to 1 during the growing season and to 0 otherwise.

- TPNDROT Temperature of ponded water [K]
- TSFSROT Ground surface temperature over subarea [K] FIXME: I think there is no need.
- ZPNDROT Depth of ponded water on surface [m]

Surface ponded water is a small term and is ephemeral in nature, so ZPNDROT and TPNDROT can both be initialized to zero. TSFSROT is included simply to provide a first guess for the surface temperature iteration in the next time step, so it can be initialized to an arbitrary value. For the snow-covered subareas of the surface it can be set to the freezing point of water; for the snow-free subareas it can be set to the temperature of the first soil layer.

## Initialization of Biogeochemical Prognostic Variables {#initCTEMProgVar}

float bleafmas(tile, icc, lat, lon) ;
  bleafmas:_FillValue = -999.f ;
  bleafmas:units = "kgC/m2" ;
  bleafmas:long_name = "Brown leaf mass" ;

float gleafmas(tile, icc, lat, lon) ;
  gleafmas:_FillValue = -999.f ;
  gleafmas:units = "kgC/m2" ;
  gleafmas:long_name = "Green leaf mass" ;

  float lfstatus(tile, icc, lat, lon) ;
    lfstatus:_FillValue = -999.f ;
    lfstatus:units = "-" ;
    lfstatus:long_name = "Leaf status, see Phenology" ;
  float litrmass(tile, iccp1, lat, lon) ;
    litrmass:_FillValue = -999.f ;
    litrmass:units = "kgC/m2" ;
    litrmass:long_name = "Litter mass per soil layer" ;

    float pandays(tile, icc, lat, lon) ;
      pandays:_FillValue = -999.f ;
      pandays:units = "-" ;
      pandays:long_name = "Days with +ve new photosynthesis, see Phenology" ;
    float rootmass(tile, icc, lat, lon) ;
      rootmass:_FillValue = -999.f ;
      rootmass:units = "kgC/m2" ;
      rootmass:long_name = "Root mass" ;

    float soilcmas(tile, iccp1, lat, lon) ;
      soilcmas:_FillValue = -999.f ;
      soilcmas:units = "kgC/m2" ;
      soilcmas:long_name = "Soil C mass per soil layer" ;
    float stemmass(tile, icc, lat, lon) ;
      stemmass:_FillValue = -999.f ;
      stemmass:units = "kgC/m2" ;
      stemmass:long_name = "Stem mass" ;

CTEM INITIALIZATION FILE. CTEM's 9 PFTs ARE
     NDL     NDL     BDL     BDL     BDL    CROP    CROP   GRASS   GRASS
     EVG     DCD     EVG DCD-CLD DCD-DRY      C3      C4      C3      C4 !Note 2 types of BDL DCD – Cold & Dry
Data for the 1st mosaic tile
   0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00 !Min. LAI for use with CTEM1 option only. Obsolete
   0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00 !Max. LAI for use with CTEM1 option only. Obsolete
   1.00    0.00    1.00    0.00    0.00    1.00    0.00    1.00    0.00 !Divide CLASS’ 4 PFTs -> 9
   0.00    0.00    0.182   0.00    0.00    0.00    0.00    0.00    0.00 !Green leaf mass for 9 PFTs
   0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00 !Brown leaf mass for 9 PFTs
   0.00    0.00    7.46    4.87    4.87    0.00    0.00    0.00    0.00 !Stem mass for 9 PFTs
   0.00    0.00    1.126   0.82    0.82    0.00    0.00    0.00    0.00 !Root mass for 9 PFTs
   0.00    0.00    0.464   0.00    0.00    0.00    0.00    0.00    0.00    0.00 !Litter mass for 9 PFTs and bare fraction
   0.00    0.00    5.644   0.00    0.00    0.00    0.00    0.00    0.00    0.00 !Soil C mass for 9 PFTs and bare fraction
      4       4       2       4       4       4       4       4       4 !Leaf status, see Phenology
      0       0       7       0       0       0       0       0       0 !Days with +ve new photosynthesis, see Phenology
Data for the 2nd mosaic tile (optional and needed only when running in the MOSAIC mode)
   0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00 !Min. LAI for use with CTEM1 option only. Obsolete
   0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00 !Max. LAI for use with CTEM1 option only. Obsolete
   1.00    0.00    1.00    0.00    0.00    1.00    0.00    1.00    0.00 !Divide CLASS’ 4 PFTs -> 9
   0.00    0.00    0.182   0.00    0.00    0.00    0.00    0.00    0.00 !Green leaf mass for 9 PFTs
   0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00 !Brown leaf mass for 9 PFTs
   0.00    0.00    7.46    4.87    4.87    0.00    0.00    0.00    0.00 !Stem mass for 9 PFTs
   0.00    0.00    1.126   0.82    0.82    0.00    0.00    0.00    0.00 !Root mass for 9 PFTs
   0.00    0.00    0.464   0.00    0.00    0.00    0.00    0.00    0.00    0.00 !Litter mass for 9 PFTs and bare fraction
   0.00    0.00    5.644   0.00    0.00    0.00    0.00    0.00    0.00    0.00 !Soil C mass for 9 PFTs and bare fraction
      4       4       2       4       4       4       4       4       4 !Leaf status, see Phenology
      0       0       7       0       0       0       0       0       0 !Days with +ve new photosynthesis, see Phenology
 Data for 2nd mosaic end here, the following  numbers correspond  to the whole grid cell
 6.191   8.996   5.230   5.005   4.705  15.311  !Jan-June, mean monthly lightning freq.
 10.437  33.310  30.767  33.250  19.165  13.850 !July-Dec, mean monthly lightning freq.
 1.00  !Fire extinguishing probability
 0.00  !Probability of fire due to human causes
 1     !Stand alone operation. 0 for GCM cell, 1 otherwise
 28.84   24.87 7931.07	1.21	3.97   ! Bioclimatic indices row 1 [Optional]
 8.03  717.44  461.50 1317.52 1573.46	6.65  ! Bioclimatic indices row 2  [Optional]
 0.00512  0.02891  0.13878  0.29462  0.45412  0.59206  0.69981  0.78309  !Slope-based fraction for dynamic wetlands



Table 4: Explanation of fields in the CTEM initialization file (.CTM)

Lines 4 and 5
The minimum and maximum LAIs in lines 4 and 5 can be set to zero, since CTEM calculates its own values. This is an old option that was designed to run with CLASS’ simulated LAI varying between these minimum and maximum values, but with CTEM simulated stomatal conductance. This option was not used very often and is no longer available. These lines will be removed in a future version.


Line 6
CLASS tracks water and energy balance of four PFTs while CTEM tracks nine PFTs, as mentioned above. The structural attributes of vegetation are lumped according to Table 4 before being passed onto CLASS. The information in line 6 is also used to divide the fractional coverage of CLASS’ four PFTs into CTEM’s nine PFTs. This functionality is essential to drive the code with CTEM switched on. Of course, when competition between PFTs is modelled then the fractional coverage of all seven natural non-crop PFTs are simulated dynamically by CTEM. As a result, when competition is switched on the fractions the fractional coverage of 4 PFTs in the .INI file and their split into 9 PFTs based on this line serves as initialization.

As an example, if there are broadleaf evergreen (25% cover) and broadleaf drought
deciduous (50% cover) trees in a grid cell then their combined fractional coverage is specified in the .INI file (0.75) and how this is broken into CTEM PFTs is specified in the .CTM file under the BDL category (0.00 0.33 0.67).


Lines 7 to 12
Specify initialization amounts for green leaf biomass, brown leaf biomass, stem biomass,
litter and soil carbon mass (Kg C/m2) for nine CTEM PFTs. When growing vegetation from
 bare ground set these to zero, of course. Note that litter and soil carbon mass are also specified over the bare fraction of the grid cell.


Line 13
The phenology subroutine of CTEM tracks leaf status according to four plants states,
1) leaf onset or maximum growth, 2) normal growth, 3) leaf fall, and 4) no leaves state.
When leaf status is set to 4 the model thinks there are no leaves. When initializing from bare ground with no vegetation set this to 4.


Line 14
The phenology subroutine of CTEM tracks days with positive net photosynthesis to
 initiate leaf onset. When net photosynthesis is positive for seven days, favourable
weather is assumed to arrive and, leaf onset begins. If this variable is set to 7 then the
model will think that favourable weather has arrived. When initializing from bare
ground with no vegetation set this variable to 0.

Note: The Previous 11 Lines are repeated for every mosaic tile in the grid cell

7th and 8th lines
from the last line
Mean monthly lightning frequency (flashes/km2.year) for the 12 months for use in the
disturbance subroutine for fire. The code at present is set to use mean monthly lightning
frequency and interpolates daily values using these monthly values. It is unlikely anyone will have a daily observation-based time series of lightning data for a grid cell. If some one does have this data they may want to modify the code to read in daily lightning frequency straight away. Of course, when both lightning frequency and probability of fire due to human causes  are set to zero then fire will not occur.


6th last line
Fire extinguishing probability varies between 0 and 1.0. A default value of 0.5 is suggested.
Setting this to 1.0 will lead to no area burned at all. Setting this to 0.0 is not allowed and setting this to a very small number will lead to the entire grid cell being burned when the fire occurs. When POPDON is set to TRUE (as explained later, which implies population density file is being read in) then population density is used to calculate fire extinguishing probability and this overwrites the value read from the .CTM file.


5th last line
Probability of fire due to human causes also varies between 0 and 1. Setting this to 0 implies that only lightning governs ignition. Setting this to 1 implies that an ignition source is always present regardless of lightning. Values between 1 and 0 represent the probability of ignition provided enough fuel is available and it is sufficiently dry to catch fire. When POPDON is set to TRUE (as explained later, which implies population density file is being read in) then population density is used to calculate probability of fire due to human causes and this overwrites the value read from the .CTM file.


4th last line
This is remnant of old code and was needed to determine if the code is being run for an
 equivalent climate model grid cell or a plot scale with a representative area of 100 km2.  The area burned by fire was then calculated based on the area of the grid cell/plot area. However, the code now does not calculate area burned but rather just outputs fraction burned. This part of the code needs to be cleaned up. The following text is the original description for this parameter.

While most CTEM quantities are calculated on a per m2 basis it is not possible to do this for fire. Fire is
spatial processes which involves some amount of area burned. The area burned due to fire in CTEM depends on area of the grid cell for which the simulation is run. When operated in the GCM, the area of the grid cell depends on the model resolution and the location of the grid cell on the surface of the Earth. When stand alone operation parameter is set to 0, the model uses the grid cell area for the GCM grid cell that encloses the latitude and the longitude provided in the .INI file. When this parameter is set to 1, the model assumes that the grid cell area is 100 km2 (which is essentially a point scale). This value may be changed in loop 91 of ctem.for.

Dividing grid cell area by annual average area burned yields an average fire return interval. Grid cell area is written to one of the CTEM output files as explained later in section 8.2.


3rd and 2nd last lines
When competition is on it needs long term values of 11 bioclimatic indices which are used to determine which PFTs can exist in a grid cell.

These 11 bioclimatic indices are 1) mean temperature of the warmest month, 2) mean temperature of the coldest month, 3) growing degree days above 5 oC, 4) annual aridity index (ratio of potential evaporation to precipitation), 5) number of months in a year with surplus water i.e. precipitation more than potential evaporation, 6)  number of months in a year with water deficit i.e. precipitation less than potential evaporation, 7) annual water deficit (mm) i.e. daily values of potential evaporation that exceed precipitation accumulated over a year, 8) annual water surplus (mm) i.e. daily values of precipitation that exceed potential evaporation accumulated over a year, 9) annual precipitation, 10) annual potential evaporation and 11) length of consecutive dry season in months, where a dry month is defined as the month in which potential evaporation exceeds precipitation.

These 11 bioclimatic indices are updated in an e-folding sense using a time scale of 25 years. However, the specification of these indices is optional (see switch settings for competition later). When the values of these indices are not provided but competition is on, then the model initializes these at zero and then as it runs it starts updating the values of these indices based on climate data in the .MET file.


  Last line
Eight slope based fractions for calculating dynamic wetland fractions. CTEM can now calculate dynamic fraction of wetland in each grid cell and its methane emissions. However, this calculation is purely diagnostic. As the soil moisture in a grid cell increases above specified thresholds then the really flat portions of the grid cell are assumed to gradually turn into wetlands. The eight slope based fractions correspond to the fraction of the grid cell that have slope less than 0.025%, 0.05%, 0.1%, 0.15%, 0.20%, 0.25%, 0.3% and 0.35% . The numbers used by CTEM are based on 1/60th degree (1 minute) resolution digital elevation data.

Just like fire, mortality is also a spatial process. For example, in CTEM approximately 1-2% of trees are killed to account for age related mortality. Clearly, trees in a small plot do not die 1% every year, but on a landscape scale on average 1% of trees may die. So when comparing model simulated biomass to point scale observations it is desirable to switch off mortality. For now, mortality is switched on only when competition is switched on (through the COMPETE switch in the job options file explained next).


### Climatic

float anndefct(lat, lon) ;
  anndefct:_FillValue = -999.f ;
  anndefct:units = "mm" ;
  anndefct:long_name = "Annual water deficit (PFTCompetition variable)" ;
float annpcp(lat, lon) ;
  annpcp:_FillValue = -999.f ;
  annpcp:units = "mm" ;
  annpcp:long_name = "Annual precipitation (PFTCompetition variable)" ;
float annsrpls(lat, lon) ;
  annsrpls:_FillValue = -999.f ;
  annsrpls:units = "mm" ;
  annsrpls:long_name = "Annual water surplus (PFTCompetition variable)" ;
float aridity(lat, lon) ;
  aridity:_FillValue = -999.f ;
  aridity:units = "-" ;
  aridity:long_name = "Aridity index, ratio of potential evaporation to precipitation (PFTCompetition variable)" ;
float defctmon(lat, lon) ;
  defctmon:_FillValue = -999.f ;
  defctmon:units = "months" ;
  defctmon:long_name = "Number of months in a year with water deficit i.e. precipitation less than potential evaporation (PFTCompetition variable)" ;
float dry_season_length(lat, lon) ;
  dry_season_length:_FillValue = -999.f ;
  dry_season_length:units = "months" ;
  dry_season_length:long_name = "Length of dry season (PFTCompetition variable)" ;
float gdd5(lat, lon) ;
  gdd5:_FillValue = -999.f ;
  gdd5:units = "days" ;
  gdd5:long_name = "Growing degree days above 5 C (PFTCompetition variable)" ;
float srplsmon(lat, lon) ;
  srplsmon:_FillValue = -999.f ;
  srplsmon:units = "months" ;
  srplsmon:long_name = "Number of months in a year with surplus water i.e. precipitation more than potential evaporation (PFTCompetition variable)" ;
float tcoldm(lat, lon) ;
  tcoldm:_FillValue = -999.f ;
  tcoldm:units = "C" ;
  tcoldm:long_name = "Temperature of the coldest month (PFTCompetition variable)" ;
float twarmm(lat, lon) ;
  twarmm:_FillValue = -999.f ;
  twarmm:units = "C" ;
  twarmm:long_name = "Temperature of the warmest month (PFTCompetition variable)" ;

  ### Peatland

  int nmtest(lat, lon) ;
    nmtest:_FillValue = -999 ;
    nmtest:long_name = "Number of tiles in each grid cell" ;

  float Cmossmas(tile, lat, lon) ;
		Cmossmas:_FillValue = -999.f ;
		Cmossmas:units = "kgC/m2" ;
		Cmossmas:long_name = "C in moss biomass" ;

    float litrmsmoss(tile, lat, lon) ;
      litrmsmoss:_FillValue = -999.f ;
      litrmsmoss:units = "kgC/m2" ;
      litrmsmoss:long_name = "Moss litter mass" ;

      float dmoss(tile, lat, lon) ;
        dmoss:_FillValue = -999.f ;
        dmoss:units = "m" ;
        dmoss:long_name = "Depth of living moss" ;
