# Creation and formatting of model input files {#makeInputFiles}

1. @ref makeMet
2. @ref makeInit
3. Greenhouse gas inputs
4. @ref makeOther

# Preparation of files for meteorological inputs {#makeMet}

CLASSIC requires the meteorological inputs as described here. The input files are netCDF format. A netcdf dump of the files will look like this:

        netcdf pres_v1.1.5_T63_chunked_1700_2017 {
        dimensions:
            lon = 128 ;
            lat = 64 ;
            time = UNLIMITED ; // (464280 currently)
        variables:
            float lon(lon) ;
                lon:standard_name = "longitude" ;
                lon:long_name = "longitude" ;
                lon:units = "degrees_east" ;
                lon:axis = "X" ;
                lon:_Storage = "contiguous" ;
                lon:_Endianness = "little" ;
            float lat(lat) ;
                lat:standard_name = "latitude" ;
                lat:long_name = "latitude" ;
                lat:units = "degrees_north" ;
                lat:axis = "Y" ;
                lat:_Storage = "contiguous" ;
                lat:_Endianness = "little" ;
            double time(time) ;
                time:standard_name = "time" ;
                time:units = "day as %Y%m%d.%f" ;
                time:calendar = "proleptic_gregorian" ;
                time:axis = "T" ;
                time:_Storage = "chunked" ;
                time:_ChunkSizes = 1 ;
                time:_Endianness = "little" ;
            float pres(time, lat, lon) ;
                pres:long_name = "Pressure" ;
                pres:units = "Pa" ;
                pres:_FillValue = 9.96921e+36f ;
                pres:missing_value = 9.96921e+36f ;
                pres:_Storage = "chunked" ;
                pres:_ChunkSizes = 464280, 8, 16 ;
                pres:_Endianness = "little" ;

        // global attributes:
                :CDI = "Climate Data Interface version 1.9.0 (http://mpimet.mpg.de/cdi)" ;
                :Conventions = "CF-1.4" ;
                :source = "Data is provided from the Japanese 55-year Reanalysis (JRA-55) project carried out by the Japan Meteorological Agency (JMA)" ;
                :institution = "Produced at the Climatic Research Unit, UEA, Norwich UK" ;
                :title = "CRUJRA Pressure: a forcing dataset based on JRA-55 on a half-degree grid" ;
                :comment = "A copy of JRA-55 year 1958 regridded to CRU 0.5 grid" ;
                :version_control = "V1.1: Regridding algorithm improvements and other corrections" ;
                :CDO = "Climate Data Operators version 1.9.0 (http://mpimet.mpg.de/cdo)" ;
                :NCO = "4.3.4" ;
                :_SuperblockVersion = 2 ;
                :_IsNetcdf4 = 1 ;
                :_Format = "netCDF-4" ;
        }
Importantly,

- Only one variable per file besides lon, lat, and time.
- Note the time units.
- The file is chunked (which depends on the grid being used. What you see here is optimal for T63 global runs). Chunking is not needed for site-level runs
- CLASSIC expects the first time step to be 0 hour 0 minute (in hh:mm:ss format - 00:00:00).

If you have ACSCII met files from a site, see [ASCII to NetCDF met file loader](@ref asciiMet) to use the provided tool to convert them to the appropriate netCDF format.

# Preparation of the model initialization file {#makeInit}

If you have the old format .INI (and CTEM's .CTM) initialization files, there is a [tool to convert them to netCDF format](@ref initTool) for use in CLASSIC. The tool itself is located in tools/initFileConverter.

If you have a netCDF format initialization/restart file that you wish to edit, you can use the [script created for that purpose](@ref modifyRS). It is located in tools/modifyRestartFile.

An ncdump of a properly formatted, global-scale, initialization file is included below. Note that for a physics-only run (CTEM off), many of the variables will not be read in. See the [manual's mainpage](@ref main) for links to sections describing the variables.

        netcdf CLASSIC_initFile_bulkdetrital_nosnow_ZoblerSDEP_competervars {
        dimensions:
            lat = 64 ;
            lon = 128 ;
            tile = 1 ;
            layer = 20 ;
            icp1 = 5 ;
            iccp1 = 10 ;
            ic = 4 ;
            icc = 9 ;
            months = 12 ;
            slope = 8 ;
        variables:
            double lat(lat) ;
                lat:units = "degrees_north" ;
                lat:long_name = "Latitude" ;
                lat:standard_name = "latitude" ;
                lat:actual_range = "-87.8637987363926, -87.8637987363926" ;
            double lon(lon) ;
                lon:units = "degrees_east" ;
                lon:long_name = "Longitude" ;
                lon:standard_name = "longitude" ;
                lon:actual_range = "0.0, 0.0" ;
            double tile(tile) ;
                tile:long_name = "tiles" ;
            double layer(layer) ;
                layer:long_name = "ground column layers" ;
            double icp1(icp1) ;
                icp1:long_name = "CLASS PFTs + bareground" ;
            double iccp1(iccp1) ;
                iccp1:long_name = "CTEM PFTs + bareground" ;
            double ic(ic) ;
                ic:long_name = "CLASS PFTs (needleleaved tree, broadleaved tree, crops, grass)" ;
            double icc(icc) ;
                icc:long_name = "CTEM PFTs (tree:NDL-EVG, NDL-DCD, BDL-EVG, BDL-COLD, BDL-DRY; crops C3,C4; grasses C3,C4)" ;
            double months(months) ;
                months:long_name = "Months" ;
            double slope(slope) ;
                slope:long_name = "wetland slope fractions for 0.025, 0.05, 0.1, 0.15, 0.20, 0.25, 0.3 and 0.35 percent slope threshold" ;
            int nmtest(lat, lon) ;
                nmtest:_FillValue = -999 ;
                nmtest:long_name = "Number of tiles in each grid cell" ;
            float FARE(tile, lat, lon) ;
                FARE:_FillValue = -999.f ;
                FARE:units = "fraction" ;
                FARE:long_name = "Tile fractional area of gridcell" ;
            float DELZ(tile, layer, lat, lon) ;
                DELZ:_FillValue = -999.f ;
                DELZ:units = "m" ;
                DELZ:long_name = "Ground layer thickness" ;
            float ZBOT(tile, layer, lat, lon) ;
                ZBOT:_FillValue = -999.f ;
                ZBOT:units = "m" ;
                ZBOT:long_name = "Depth of bottom of ground layer" ;
            float ALIC(tile, icp1, lat, lon) ;
                ALIC:_FillValue = -999.f ;
                ALIC:units = "-" ;
                ALIC:long_name = "Average near-IR albedo of vegetation category when fully-leafed" ;
            float ALVC(tile, icp1, lat, lon) ;
                ALVC:_FillValue = -999.f ;
                ALVC:units = "-" ;
                ALVC:long_name = "Average visible albedo of vegetation category when fully-leafed" ;
            float CMAS(tile, ic, lat, lon) ;
                CMAS:_FillValue = -999.f ;
                CMAS:units = "$[kg m^{-2} ]$" ;
                CMAS:long_name = " Annual maximum canopy mass for vegetation category" ;
            float LNZ0(tile, icp1, lat, lon) ;
                LNZ0:_FillValue = -999.f ;
                LNZ0:units = "-" ;
                LNZ0:long_name = "Natural logarithm of maximum vegetation roughness length" ;
            float PAMN(tile, ic, lat, lon) ;
                PAMN:_FillValue = -999.f ;
                PAMN:units = "-" ;
                PAMN:long_name = "Annual minimum plant area index of vegetation category" ;
            float PAMX(tile, ic, lat, lon) ;
                PAMX:_FillValue = -999.f ;
                PAMX:units = "-" ;
                PAMX:long_name = "Annual maximum plant area index of vegetation category" ;
            float ROOT(tile, ic, lat, lon) ;
                ROOT:_FillValue = -999.f ;
                ROOT:units = "m" ;
                ROOT:long_name = "Annual maximum rooting depth of vegetation category" ;
            int SOCI(tile, lat, lon) ;
                SOCI:_FillValue = -999 ;
                SOCI:units = "index" ;
                SOCI:long_name = "Soil colour index" ;
            float grclarea(lat, lon) ;
                grclarea:_FillValue = -999.f ;
                grclarea:units = "km2" ;
                grclarea:long_name = "Area of grid cell" ;
            float SDEP(tile, lat, lon) ;
                SDEP:_FillValue = -999.f ;
                SDEP:units = "m" ;
                SDEP:long_name = "Soil permeable depth" ;
            float fcancmx(tile, icc, lat, lon) ;
                fcancmx:_FillValue = -999.f ;
                fcancmx:units = "-" ;
                fcancmx:long_name = "PFT fractional coverage per grid cell" ;
            float SAND(tile, layer, lat, lon) ;
                SAND:_FillValue = -999.f ;
                SAND:units = "%" ;
                SAND:long_name = "Percentage sand content" ;
            float CLAY(tile, layer, lat, lon) ;
                CLAY:_FillValue = -999.f ;
                CLAY:units = "%" ;
                CLAY:long_name = "Percentage clay content" ;
            float ZRFM(lat, lon) ;
                ZRFM:_FillValue = -999.f ;
                ZRFM:units = "m" ;
                ZRFM:long_name = "Reference height associated with forcing wind speed" ;
            float ZRFH(lat, lon) ;
                ZRFH:_FillValue = -999.f ;
                ZRFH:units = "m" ;
                ZRFH:long_name = "Reference height associated with forcing air temperature and humidity" ;
            float ZBLD(lat, lon) ;
                ZBLD:_FillValue = -999.f ;
                ZBLD:units = "m" ;
                ZBLD:long_name = "Atmospheric blending height for surface roughness length averaging" ;
            int GC(lat, lon) ;
                GC:_FillValue = -999 ;
                GC:units = "-" ;
                GC:long_name = "GCM surface descriptor - land surfaces (inc. inland water) is -1" ;
            float FCAN(tile, icp1, lat, lon) ;
                FCAN:_FillValue = -999.f ;
                FCAN:units = "-" ;
                FCAN:long_name = "Annual maximum fractional coverage of modelled area (read in for CLASS only runs)" ;
            float RSMN(tile, ic, lat, lon) ;
                RSMN:_FillValue = -999.f ;
                RSMN:units = "s/m" ;
                RSMN:long_name = "Minimum stomatal resistance of vegetation category" ;
            float QA50(tile, ic, lat, lon) ;
                QA50:_FillValue = -999.f ;
                QA50:units = "W/m2" ;
                QA50:long_name = "Reference value of incoming shortwave radiation (used in stomatal resistance calculation)" ;
            float VPDA(tile, ic, lat, lon) ;
                VPDA:_FillValue = -999.f ;
                VPDA:units = "-" ;
                VPDA:long_name = "Vapour pressure deficit coefficient (used in stomatal resistance calculation)" ;
            float VPDB(tile, ic, lat, lon) ;
                VPDB:_FillValue = -999.f ;
                VPDB:units = "-" ;
                VPDB:long_name = "Vapour pressure deficit coefficient (used in stomatal resistance calculation)" ;
            float PSGA(tile, ic, lat, lon) ;
                PSGA:_FillValue = -999.f ;
                PSGA:units = "-" ;
                PSGA:long_name = " Soil moisture suction coefficient (used in stomatal resistance calculation)" ;
            float PSGB(tile, ic, lat, lon) ;
                PSGB:_FillValue = -999.f ;
                PSGB:units = "-" ;
                PSGB:long_name = " Soil moisture suction coefficient (used in stomatal resistance calculation)" ;
            float DRN(tile, lat, lon) ;
                DRN:_FillValue = -999.f ;
                DRN:units = "-" ;
                DRN:long_name = "Soil drainage index" ;
            float XSLP(tile, lat, lon) ;
                XSLP:_FillValue = -999.f ;
                XSLP:units = "-" ;
                XSLP:long_name = "Not in Use: parameters lateral movement of soil water" ;
            float GRKF(tile, lat, lon) ;
                GRKF:_FillValue = -999.f ;
                GRKF:units = "-" ;
                GRKF:long_name = "Not in Use: parameters lateral movement of soil water" ;
            float WFSF(tile, lat, lon) ;
                WFSF:_FillValue = -999.f ;
                WFSF:units = "-" ;
                WFSF:long_name = "Not in Use: parameters lateral movement of soil water" ;
            float WFCI(tile, lat, lon) ;
                WFCI:_FillValue = -999.f ;
                WFCI:units = "-" ;
                WFCI:long_name = "Not in Use: parameters lateral movement of soil water" ;
            int MID(tile, lat, lon) ;
                MID:_FillValue = -999 ;
                MID:units = "-" ;
                MID:long_name = "Mosaic tile type identifier (1 for land surface, 0 for inland lake)" ;
            float ORGM(tile, layer, lat, lon) ;
                ORGM:_FillValue = -999.f ;
                ORGM:units = "%" ;
                ORGM:long_name = "Percentage organic matter content" ;
            float TBAR(tile, layer, lat, lon) ;
                TBAR:_FillValue = -999.f ;
                TBAR:units = "C" ;
                TBAR:long_name = "Temperature of soil layers" ;
            float THLQ(tile, layer, lat, lon) ;
                THLQ:_FillValue = -999.f ;
                THLQ:units = "m3/m3" ;
                THLQ:long_name = "Volumetric liquid water content of soil layers" ;
            float THIC(tile, layer, lat, lon) ;
                THIC:_FillValue = -999.f ;
                THIC:units = "m3/m3" ;
                THIC:long_name = "Volumetric frozen water content of soil layers" ;
            float TCAN(tile, lat, lon) ;
                TCAN:_FillValue = -999.f ;
                TCAN:units = "C" ;
                TCAN:long_name = "Vegetation canopy temperature" ;
            float TSNO(tile, lat, lon) ;
                TSNO:_FillValue = -999.f ;
                TSNO:units = "C" ;
                TSNO:long_name = "Snowpack temperature" ;
            float TPND(tile, lat, lon) ;
                TPND:_FillValue = -999.f ;
                TPND:units = "C" ;
                TPND:long_name = "Temperature of ponded water" ;
            float ZPND(tile, lat, lon) ;
                ZPND:_FillValue = -999.f ;
                ZPND:units = "m" ;
                ZPND:long_name = "Depth of ponded water on surface" ;
            float RCAN(tile, lat, lon) ;
                RCAN:_FillValue = -999.f ;
                RCAN:units = "kg/m2" ;
                RCAN:long_name = "Intercepted liquid water stored on canopy" ;
            float SCAN(tile, lat, lon) ;
                SCAN:_FillValue = -999.f ;
                SCAN:units = "kg/m2" ;
                SCAN:long_name = "Intercepted frozen water stored on canopy" ;
            float SNO(tile, lat, lon) ;
                SNO:_FillValue = -999.f ;
                SNO:units = "kg/m2" ;
                SNO:long_name = "Mass of snow pack" ;
            float ALBS(tile, lat, lon) ;
                ALBS:_FillValue = -999.f ;
                ALBS:units = "-" ;
                ALBS:long_name = "Snow albedo" ;
            float RHOS(tile, lat, lon) ;
                RHOS:_FillValue = -999.f ;
                RHOS:units = "kg/m3" ;
                RHOS:long_name = "Density of snow" ;
            float GRO(tile, lat, lon) ;
                GRO:_FillValue = -999.f ;
                GRO:units = "-" ;
                GRO:long_name = "Vegetation growth index" ;
            float gleafmas(tile, icc, lat, lon) ;
                gleafmas:_FillValue = -999.f ;
                gleafmas:units = "kgC/m2" ;
                gleafmas:long_name = "Green leaf mass" ;
            float bleafmas(tile, icc, lat, lon) ;
                bleafmas:_FillValue = -999.f ;
                bleafmas:units = "kgC/m2" ;
                bleafmas:long_name = "Brown leaf mass" ;
            float stemmass(tile, icc, lat, lon) ;
                stemmass:_FillValue = -999.f ;
                stemmass:units = "kgC/m2" ;
                stemmass:long_name = "Stem mass" ;
            float rootmass(tile, icc, lat, lon) ;
                rootmass:_FillValue = -999.f ;
                rootmass:units = "kgC/m2" ;
                rootmass:long_name = "Root mass" ;
            float litrmass(tile, iccp1, lat, lon) ;
                litrmass:_FillValue = -999.f ;
                litrmass:units = "kgC/m2" ;
                litrmass:long_name = "Litter mass per soil layer" ;
            float soilcmas(tile, iccp1, lat, lon) ;
                soilcmas:_FillValue = -999.f ;
                soilcmas:units = "kgC/m2" ;
                soilcmas:long_name = "Soil C mass per soil layer" ;
            float lfstatus(tile, icc, lat, lon) ;
                lfstatus:_FillValue = -999.f ;
                lfstatus:units = "-" ;
                lfstatus:long_name = "Leaf status, see Phenology" ;
            float pandays(tile, icc, lat, lon) ;
                pandays:_FillValue = -999.f ;
                pandays:units = "-" ;
                pandays:long_name = "Days with +ve new photosynthesis, see Phenology" ;
            float slopefrac(slope, tile, lat, lon) ;
                slopefrac:_FillValue = -999.f ;
                slopefrac:units = "-" ;
                slopefrac:long_name = "Slope-based fraction for dynamic wetlands" ;
            float rice(months, lat, lon) ;
                rice:_FillValue = -999.f ;
                rice:units = "-" ;
                rice:long_name = "Monthly irrigated rice ag. gridcell fraction" ;
            float ipeatland(tile, lat, lon) ;
                ipeatland:_FillValue = -999.f ;
                ipeatland:units = "-" ;
                ipeatland:long_name = "Peatland flag: 0 = not a peatland, 1= bog, 2 = fen" ;
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
            float twarmm(lat, lon) ;
                twarmm:_FillValue = -999.f ;
                twarmm:units = "C" ;
                twarmm:long_name = "Temperature of the warmest month (PFTCompetition variable)" ;
            float tcoldm(lat, lon) ;
                tcoldm:_FillValue = -999.f ;
                tcoldm:units = "C" ;
                tcoldm:long_name = "Temperature of the coldest month (PFTCompetition variable)" ;
            float gdd5(lat, lon) ;
                gdd5:_FillValue = -999.f ;
                gdd5:units = "days" ;
                gdd5:long_name = "Growing degree days above 5 C (PFTCompetition variable)" ;
            float aridity(lat, lon) ;
                aridity:_FillValue = -999.f ;
                aridity:units = "-" ;
                aridity:long_name = "Aridity index, ratio of potential evaporation to precipitation (PFTCompetition variable)" ;
            float srplsmon(lat, lon) ;
                srplsmon:_FillValue = -999.f ;
                srplsmon:units = "months" ;
                srplsmon:long_name = "Number of months in a year with surplus water i.e. precipitation more than potential evaporation (PFTCompetition variable)" ;
            float defctmon(lat, lon) ;
                defctmon:_FillValue = -999.f ;
                defctmon:units = "months" ;
                defctmon:long_name = "Number of months in a year with water deficit i.e. precipitation less than potential evaporation (PFTCompetition variable)" ;
            float anndefct(lat, lon) ;
                anndefct:_FillValue = -999.f ;
                anndefct:units = "mm" ;
                anndefct:long_name = "Annual water deficit (PFTCompetition variable)" ;
            float annsrpls(lat, lon) ;
                annsrpls:_FillValue = -999.f ;
                annsrpls:units = "mm" ;
                annsrpls:long_name = "Annual water surplus (PFTCompetition variable)" ;
            float annpcp(lat, lon) ;
                annpcp:_FillValue = -999.f ;
                annpcp:units = "mm" ;
                annpcp:long_name = "Annual precipitation (PFTCompetition variable)" ;
            float dry_season_length(lat, lon) ;
                dry_season_length:_FillValue = -999.f ;
                dry_season_length:units = "months" ;
                dry_season_length:long_name = "Length of dry season (PFTCompetition variable)" ;

        // global attributes:
                :history = "Created by /home/rjm/Documents/CTEM/Model_restart_netcdf/make_model_init_netcdf.py" ;
                :creation_date = "Tue Aug 14 14:03:11 2018" ;
                :created_by = "Joe Melton, CCCma" ;
        }

# Making input files for other input variables {#makeOther}

Most other CLASSIC input files are not desired or needed for point runs of the model. If you require input files for regional or global simulations we can likely provide you with versions we use in our runs. Please contact joe.melton@canada.ca for access to the files we presently have available. The possible inputs are listed in Additional inputs depending on model configuration. Because these files are not required of most users we have not set up tools to help generate the files.
