
# Preparing a CLASSIC run {#runPrep}

1. @ref compilingMod
2. @ref setupJobOpts
3. @ref xmlSystem

----

# Compiling CLASSIC for serial and parallel simulations {#compilingMod}

CLASSIC's Makefile (yourpath/CLASSIC/Makefile) is setup to allow easy compilation for serial or parallel model running.

There are three options when compiling:

- For a serial run, use the command "make mode=serial" or just "make"
- For a parallel run, use either the command "make mode=parallel"
- For a supercomputer run, use the command "make mode=supercomputer"

The serial compilation is presently set to use the GNU compiler (gfortran). Parallel compilation uses the GNU MPI compiler mpif90. The supercomputer compilation is setup for the ECCC supercomputer environment and is unlikely to be useful for non-ECCC users.

Upon compilation all object (.o) and module (.mod) files are placed in the objectFiles folder. The executable in placed in the bin folder.

A useful command is 'make clean', which removes all *.o *.mod and the model binary. This can allow a fresh compilation which can be handy if some parameters are changed that aren't being refreshed on a make.

# Setting up the joboptions file {#setupJobOpts}

The joboptions file controls the model configuration, inputs files used, and model outputs. The template joboptions file is located in the configurationFiles folder. Use this as your starting point.

If, for example, we wanted to run CLASSIC at a site with observed meteorology, with leap years, from 1991 to 2017 then we can set up the meteorological options as follows:

        &joboptions

        ! Meteorological options:
            readMetStartYear = 1991  !< First year of meteorological forcing to read in from the met file
            readMetEndYear = 2017    !< Last year of meteorological forcing to read in from the met file
            metLoop = 1 ,            !< no. of times to cycle over the read-in meteorology
            leap = .true. ,         !< True if your meteorological forcing includes leap years

If we are interested in spinning up the C pools, we could set metLoop to run the model over the readMetStartYear to readMetEndYear a metLoop number of times. We should also point to our meteorological forcing files like,

        ! Meteorological forcing files:
            metFileFss = pathToFile/dswrf.nc',        !< location of the incoming shortwave radiation meteorology file
            metFileFdl = 'pathToFile/dlwrf.nc',        !< location of the incoming longwave radiation meteorology file
            metFilePre = 'pathToFile/pre.nc',        !< location of the precipitation meteorology file
            metFileTa = 'pathToFile/tmp.nc',         !< location of the air temperature meteorology file
            metFileQa = 'pathToFile/spfh.nc',         !< location of the specific humidity meteorology file
            metFileUv = 'pathToFile/wind.nc',         !< location of the wind speed meteorology file
            metFilePres = 'pathToFile/pres.nc',       !< location of the atmospheric pressure meteorology file

The model initialization and restart files need to be pointed to. Note the rs_file_to_overwrite is **overwritten**. The simplest thing to do at a start of a run is to make a copy of the init_file to be the rs_file_to_overwrite. CLASSIC only overwrites the prognostic variables leaving any others unchanged (see model_state_drivers.write_restart). The model parameters file also needs to be pointed to.

        ! Initialization and restart files
            init_file = 'inputFiles/CLASSCTEM_initialization_bulkdetrital_nosnow.nc' ,     !< location of the model initialization file
            rs_file_to_overwrite = 'rsFile.nc' ,       !< location of the existing netcdf file that will be **overwritten** for the restart file
                                                       !! typically here you will just copy the init_file and rename it so it can be overwritten.
        ! Namelist of model parameters
            runparams_file = 'configurationFiles/template_run_parameters.txt' ,    !< location of the namelist file containing the model parameters

The next series of switches relate to CTEM. If you wish to run a physics only run with prescribed vegetation parameters, ctem_on is set to .false. and the remainder of this section is ignored. The physics switches are farther down.

        ! CTEM (biogeochemistry) switches:
         ctem_on = .true. ,     !< set this to true for using ctem simulated dynamic lai and canopy mass, else class simulated specified
                                !< lai and canopy mass are used. with this switch on, all the main ctem subroutines are run.
            icc = 9 ,           !< Number of CTEM level PFTs. NOTE: The number specified here must match the data in your init netcdf file.
            l2max = 3 ,         !< Maximum number of level 2 CTEM PFTs. This is the maximum number of CTEM PFTs associated with a single CLASS PFT.
            spinfast = 3 ,      !< Set this to a higher number up to 10 to spin up soil carbon pool faster. Set to 1 for final round of spin up and transient runs.

The CO2 and CH4 switches behave similarly. If a constant [CO2] is desired, transientCO2 is set to .false. and the year of observed CO2 is specified in fixedYearCO2. The [CO2] of the corresponding year is selected from the CO2File and used for the run. Simlarly for CH4. If are interested in simulating methane related variables then simply copy your CO2File to be your CH4file (e.g. cp CO2file.nc fakeCH4file.nc, and then specify the fakeCH4file.nc for the CH4File since the model expects a file there). If for some reason you wish to run with non-historical [CO2], you could edit your CO2File (using a combination of the NCO tools ncdump and ncgen, for example) to include a future year with associated CO2 value (like 2050 and 500ppm for example)

        !CO2 switches:
            transientCO2 = .true. , !< Read in time varying CO2 concentration from CO2File or if set to false then use only the year of fixedYearCO2 value
            CO2File = 'inputFiles/mole_fraction_of_carbon_dioxide_in_air_input4MIPs_GHGConcentrations_CMIP_UoM-CMIP-1-2-0_gr1-GMNHSH_1850-2014_absTimeIndex.nc' ,
            fixedYearCO2 = 1850 ,   !< If transientCO2 = .true., this is ignored.

        !CH4 switches:
            ! If you don't care about running with methane, you can just copy your CO2 file (cp CO2file.nc fakeCH4file.nc) and point to it here.
            ! It won't affect the rest of your run. However you do need to have a file specified if CTEM is on.
            transientCH4 = .true. , !< Read in time varying CH4 concentration from CH4File or if set to false then use only the year of fixedYearCH4 value
            CH4File = 'inputFiles/mole_fraction_of_methane_in_air_input4MIPs_GHGConcentrations_CMIP_UoM-CMIP-1-2-0_gr1-GMNHSH_1850-2014_absTimeAxis.nc',
            fixedYearCH4 = 1850 ,   !If transientCH4 = .true., this is ignored.

Disturbance in the form of fire is optional in CLASSIC. If fire is turned on then the model requires a population density input file (see @ref initPopd) which is used in a manner similar to CO2 and CH4, i.e. it can use a fixed or transient value. Lightning strikes (see @ref initLightFire) are also required and can be specified similarly to population density.

        !Fire switches
            dofire = .true. ,               !< If true the fire disturbance parameterization is turned on.

                transientPOPD = .true. ,    !< Read in time varying population density from POPDFile or if set to false then use only the year of fixedYearPOPD.
                POPDFile = 'inputFiles/POPD_annual_1850_2016_T63_Sep_2017_chunked.nc' ,
                fixedYearPOPD = 1850 ,      !< If transientPOPD = .true., this is ignored.

                transientLGHT= .true.      !< use lightning strike time series, otherwise use fixedYearLGHT
                LGHTFile = 'inputFiles/lisotd_1995_2014_climtlgl_lghtng_as_ts_1850_2050_chunked.nc' , !< Location of the netcdf file containing lightning strike values
                fixedYearLGHT = 1850 ,    !< set the year to use for lightning strikes if transientLGHT is false.

Competition for space between plant functional types is parameterized in CLASSIC. If PFTCompetition is set to false, the PFT fractional coverage follows rules as will be outlined next. If PFTCompetition is true then CLASSIC has two more switches of interest. Competition uses bioclimatic indices to determine whether PFTs should be allowed to compete within a gridcell (see competition_scheme.existence and @ref initClimComp). The bioclimatic indices are either read-in from the init_file or are calculated anew for the run underway. If the values are in the init_file from a spinup run then set inibioclim to true, otherwise set inibioclim to false. It is also possible to start the model from bare ground (rather than the PFT configuration found in your init_file) by setting start_bare to true.

        ! Competition switches:
            PFTCompetition = .false. ,      !< If true, competition between PFTs for space on a grid cell is implimented
                inibioclim = .false. ,      !< set this to true if competition between pfts is to be implimented and you have the mean climate values
                                            !< in the init netcdf file.
                start_bare = .false.,       !< Set this to true if competition is true, and if you wish to start from bare ground. if this is set to false, the
                                            !< init netcdf file info will be used to set up the run. NOTE: This still keeps the crop fractions
                                            !< (while setting all pools to zero)

Land use change is possible via a LUCFile (see @ref inputLUC) that has the fractional coverage for each PFT annually. This switch has an additional option as described in the comment below. Pay attention here.

        ! Land Use switches:
            !** If you wish to use your own PFT fractional covers (specified in the init_file), set fixedYearLUC to -9999, otherwise set it
            ! to the year of land cover you want to use. If you wish to have transient land cover changes, set
            ! lnduseon to true, it will update the fractional coverages from LUCFile. When lnduseon is false it is
            ! not updated beyond the initial read in of landcover for fixedYearLUC, or if -9999 then the LUCFile is
            ! not used at all.**
            lnduseon = .true. ,
            LUCFile = 'inputFiles/LUH_HYDE_based_crop_area_adjusted_land_cover_CTEM_fractions_1850_2017_T63_chunked.nc' ,
            fixedYearLUC = 1901 ,

CLASSIC can determine dynamics wetland locations for wetland methane emissions. Alternatively CLASSIC can read in time evolving wetland fractions from an external file.

        ! Wetland switches:
            ! If you wish to read in and use observed wetland fractions, there are two options. If you wish time
            ! evolving wetland fractions set transientOBSWETF to true and give a OBSWETFFile. If you wish to use
            ! a single year of that file set transientOBSWETF to false, give a OBSWETFFile, and set fixedYearOBSWETF
            ! to some valid year. If you wish to use only dynamically determined wetland fractions set transientOBSWETF
            ! to false and set fixedYearOBSWETF to -9999. The slope fractions in the init_file will then be used to
            ! dynamically determine wetland extent.
            transientOBSWETF = .false. ,  !< use observed wetland fraction time series, otherwise use fixedYearOBSWETF
            OBSWETFFile = '',             !< Location of the netcdf file containing observed wetland fraction
            fixedYearOBSWETF = -9999 ,    !< set the year to use for observed wetland fraction if transientOBSWETF is false.

CLASS switches determine the configuration of the physics only as well as CLASS+CTEM (physics and biogeochemistry) runs. COMBAK

        ! Physics switches:

            ican = 4 ,     !< Number of PFTs considered by the physics subroutines. NOTE: The number specified here must match the data in your init netcdf file.
            IDISP = 0 ,    !< if idisp=0, vegetation displacement heights are ignored, because the atmospheric model considers these to be part
                            !< of the "terrain". if idisp=1, vegetation displacement heights are calculated.
            IZREF = 2 ,    !< if izref=1, the bottom of the atmospheric model is taken lie at the ground surface.
                            !< if izref=2, the bottom of the atmospheric model is taken to lie at the local roughness height.
            ISLFD = 0 ,    !< if islfd=0, drcoef is called for surface stability corrections and the original gcm set of screen-level diagnostic calculations
                            !< is done. if islfd=1, drcoef is called for surface stability corrections and sldiag is called for screen-level diagnostic calculations.
                            !< if islfd=2, flxsurfz is called for surface stability corrections and diasurf is called for screen-level diagnostic calculations.
            IPCP = 1 ,     !< if ipcp=1, the rainfall-snowfall cutoff is taken to lie at 0 C. if ipcp=2, a linear partitioning of precipitation between
                            !< rainfall and snowfall is done between 0 C and 2 C. if ipcp=3, rainfall and snowfall are partitioned according to
                            !< a polynomial curve between 0 C and 6 C.
            IWF = 0 ,      !< if iwf=0, only overland flow and baseflow are modelled, and the ground surface slope is not modelled. if iwf=n (0<n<4) ,
                           !< the watflood calculations of overland flow and interflow are performed; interflow is drawn from the top n soil layers.
            isnoalb = 0 ,  !< if isnoalb is set to 0, the original two-band snow albedo algorithms are used. if it is set to 1, the new four-band routines are used.

The iteration scheme for canopy or ground surface temperatures can be either a bisection or Newton-Raphson method. This is usually set to bisection

         ! Iteration scheme
            !< ITC, ITCG and ITG are switches to choose the iteration scheme to be used in calculating the canopy or ground surface temperature
            !< respectively.  if the switch is set to 1, a bisection method is used; if to 2, the newton-raphson method is used.
            ITC = 1 ,   !< Canopy
            ITCG = 1 ,  !< Ground under canopy
            ITG = 1 ,   !< Ground

User supplied values can be used for plant area index, vegetation height, and canopy, soil, or snow albedos. If any of these inputs are supplied the model_state_drivers.f90 needs to be adapted to read-in the user-supplied values.

         ! User-supplied values:
            !< if ipai, ihgt, ialc, ials and ialg are zero, the values of plant area index, vegetation height, canopy albedo, snow albedo
            !< and soil albedo respectively calculated by class are used. if any of these switches is set to 1, the value of the
            !< corresponding parameter calculated by class is overridden by a user-supplied input value.
            IPAI = 0 ,  !< Plant area index
            IHGT = 0 ,  !< Vegetation height
            IALC = 0 ,  !< Canopy albedo
            IALS = 0 ,  !< Snow albedo
            IALG = 0 ,  !< Soil albedo

Model outputs are in netcdf format. The outputs metadata is read in from an xmlFile (see @ref xmlSystem). Outputs can be grid-cell average (default), per PFT, or per tile. In all cases the output has to be specified in the xmlFile and also properly handled in prepareOutputs.f90. Temporal resolution of output files are half-hourly (for physics variables as well as photosynthesis and canopy conductance only), daily, monthly and annually. For half-hourly and daily outputs the start and end days as well as start and end years can be specified. Monthly file can specify the year to start writing the outputs.

        ! Output options:

            output_directory = 'outputFiles' ,        !< Directory where the output netcdfs will be placed
            xmlFile = 'configurationFiles/outputVariableDescriptors_v1.2.xml' ,  !< location of the xml file that outlines the possible netcdf output files

            doperpftoutput = .true. ,   !< Switch for making extra output files that are at the per PFT level
            dopertileoutput = .false. , !< Switch for making extra output files that are at the per tile level

            dohhoutput = .false. ,      !< Switch for making half hourly output files (annual are always outputted)
            JHHSTD = 166 ,                !< day of the year to start writing the half-hourly output
            JHHENDD = 185 ,             !< day of the year to stop writing the half-hourly output
            JHHSTY = 1901 ,             !< simulation year (iyear) to start writing the half-hourly output
            JHHENDY = 1901 ,            !< simulation year (iyear) to stop writing the half-hourly output

            dodayoutput = .false. ,     !< Switch for making daily output files (annual are always outputted)
            JDSTD = 20 ,                 !< day of the year to start writing the daily output
            JDENDD = 30 ,              !< day of the year to stop writing the daily output
            JDSTY = 1902 ,              !< simulation year (iyear) to start writing the daily output
            JDENDY = 1903 ,             !< simulation year (iyear) to stop writing the daily output

            domonthoutput = .true. ,    !< Switch for making monthly output files (annual are always outputted)
            JMOSTY = 1901 ,             !< Year to start writing out the monthly output files.

Comments can be added to output files using the Comment field below. Also comments can be left in the joboptions file after the backslash.

            Comment = ' test '          !< Comment about the run that will be written to the output netcdfs

         /

        This area can be used for comments about the run.
