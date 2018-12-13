!> Contains CTEM globally accessible parameters
module ctem_params

!>\ingroup ctem_params_main
!!@{

! J. Melton
! Jun 23 2013

! Change History:

! Sep 2018 - EC - Relocate variables defined in "readin_params" to the main module and allocate them
!                 in "allocateParamsCTEM". Cray compiler doesn't currently support namelists using 
!                 dynamic arrays dimensioned by a parameter specified in the main module.
! Jun 2017 - JM - Change to using namelist.
! Jun 11 2015 - JM - Add in new disturb params (extnmois_duff and extnmois_veg) replacing duff_dry and extnmois.
! Mar 12 2014 - JM - Allow for two sets of paramters (competition and prescribed PFT fractions).
!                    all ctem subroutines except PHTSYN3 keep their parameters here.
!
! Jan 17 2014 - JM - Add in more parameters from ctem.f, phenology.f, and allocate.f.

implicit none


public :: allocateParamsCTEM
public :: readin_params
public :: prepareGlobalParams


! Constants

real, parameter :: zero     = 1.0e-20
real, parameter :: abszero  = 1e-12     !<this one is for runclassctem.f and allocate.f
real, parameter :: pi       = 3.1415926535898d0
real, parameter :: earthrad = 6371.22   !<radius of earth, km
real, parameter :: km2tom2  = 1.0e+06   !<changes from \f$km^2\f$ to \f$m^2\f$
real, parameter :: deltat   = 1.0       !<CTEM's time step in days
real, parameter :: convertkgC = 1.201e-8 !< Converts from umolCO2/m2/s to kgC/m2/s

! These month arrays are possibly overwritten in runclassctem due to leap years.
integer, dimension(12) :: monthdays = [ 31,28,31,30,31,30,31,31,30,31,30,31 ] !< days in each month
integer, dimension(13) :: monthend  = [ 0,31,59,90,120,151,181,212,243,273,304,334,365 ] !< calender day at end of each month
integer, dimension(12) :: mmday     = [ 16,46,75,106,136,167,197,228,259,289,320,350 ] !<mid-month day
integer, parameter :: nmon = 12 !< Number of months in a year

! ----
! Read in from the netcdf initialization file:
integer :: nlat        !< Number of cells we are running, read in from the initialization file. Offline this always 1.
integer :: nmos        !< Number of mosaic tiles, read in from the initialization file
integer :: ilg         !< nlat x nmos
integer :: ignd        !< Number of soil layers, read in from the initialization file

! Read in from the job options file:
integer :: ican              !< Number of CLASS (physics) pfts, read in from the job options file.
integer :: icc               !< Number of CTEM (biogeochemical) pfts, read in from the job options file.
integer :: l2max             !< Maximum number of level 2 CTEM PFTs. This is the maximum number of CTEM PFTs
                             !! associated with a single CLASS PFT. Read in from the job options file.

! ----
! Plant-related parameters that are calculated based on job options
integer :: icp1              !< ican + 1
integer :: iccp1             !< iccp1
integer :: kk                !< product of class pfts and l2max
integer :: numcrops          !< number of crop pfts
integer :: numtreepfts       !< number of tree pfts
integer :: numgrass          !< number of grass pfts
integer :: numshrubs         !< number of shrubs pfts
integer, dimension(:), allocatable :: nol2pfts  !< Number of level 2 PFTs calculated in readin_params
logical, dimension(:), allocatable :: crop      !< simple crop matrix, define the number and position of the crops (NOTE: dimension icc)
logical, dimension(:), allocatable :: grass     !< simple grass matric, define the number and position of grass (NOTE: dimension icc)

real, parameter :: seed    = 0.001    !< seed pft fraction, same as in competition \nin mosaic mode, all tiles are given this as a minimum
real, parameter :: minbare = 1.0e-5   !< minimum bare fraction when running competition on to prevent numerical problems.
real, parameter :: c2dom   = 450.0    !< gc / kg dry organic matter \nconversion factor from carbon to dry organic matter value is from Li et al. 2012 biogeosci
real, parameter :: wtCH4   = 16.044   !< Molar mass of CH4 (\f$g mol^{-1}\f$)

integer, parameter :: nbs = 4         !<Number of modelled shortwave radiation wavelength bands COMBAK Move to namelist? Useful? FLAG

real :: gasc = 8.314    !< gas constant (\f$J mol^{-1} K^{-1}\f$)  !COMBAK FLAG is there no global one I can use?
real :: tolrance = 0.0001d0 !< our tolerance for balancing c budget in kg c/m2 in one day (differs when competition on or not)
                            ! YW May 12, 2015 in peatland the C balance gap reaches 0.00016.
real :: tolrnce1 = 0.5      !< kg c, tolerance of total c balance (FOR LUC) !IDEA FLAG would be good to make this consistent with global tolerance so only one value.

!> Logical switch for using constant allocation factors (default value is false)
logical :: consallo = .false.

! ============================================================
! Read in from the namelist: ---------------------------------

integer, dimension(:), allocatable :: modelpft      !<Separation of pfts into level 1 (for class) and level 2 (for ctem) pfts.
character(8), dimension(:), allocatable :: pftlist  !<List of PFTs
character(8), dimension(:), allocatable :: vegtype  !<Type of vegetation, options: Tree, Grass, Crop, Shrub
real, dimension(:), allocatable :: kn               !< Canopy light/nitrogen extinction coefficient

! Moved definitions in "readin_params" here and allocate via "allocateParamsCTEM".

real, dimension(:), allocatable :: omega_compete
real, dimension(:), allocatable :: epsilonl_compete
real, dimension(:), allocatable :: epsilons_compete
real, dimension(:), allocatable :: epsilonr_compete
real, dimension(:), allocatable :: bsrtstem_compete
real, dimension(:), allocatable :: bsrtroot_compete
real, dimension(:), allocatable :: mxmortge_compete
real, dimension(:), allocatable :: maxage_compete
real, dimension(:), allocatable :: drlsrtmx_compete

!allocate.f parameters: ----------

real, dimension(:), allocatable :: omega            !< omega, parameter used in allocation formulae (values differ if using prescribed vs. competition run)
real, dimension(:), allocatable :: epsilonl         !< Epsilon leaf, parameter used in allocation formulae (values differ if using prescribed vs. competition run)
real, dimension(:), allocatable :: epsilons         !< Epsilon stem, parameter used in allocation formulae (values differ if using prescribed vs. competition run)
real, dimension(:), allocatable :: epsilonr         !< Epsilon root, parameter used in allocation formulae (values differ if using prescribed vs. competition run)
real, dimension(:), allocatable :: rtsrmin          !< Minimum root:shoot ratio mostly for support and stability
real, dimension(:), allocatable :: aldrlfon         !< Allocation to leaves during leaf onset
real, dimension(:), allocatable :: caleaf           !< Constant allocation fractions to leaves if not using dynamic allocation.
                                        !!(NOT thoroughly tested, and using dynamic allocation is preferable)
real, dimension(:), allocatable :: castem           !< Constant allocation fractions to stem if not using dynamic allocation.
                                        !!(NOT thoroughly tested, and using dynamic allocation is preferable)
real, dimension(:), allocatable :: caroot           !< Constant allocation fractions to roots if not using dynamic allocation.
                                        !!(NOT thoroughly tested, and using dynamic allocation is preferable)

! bio2str.f parameters: ---------

real, dimension(:), allocatable :: abar             !< parameter determining average root profile
real, dimension(:), allocatable :: avertmas         !< average root biomass (kg c/m2) for ctem's 8 pfts used for estimating rooting profile
real, dimension(:), allocatable :: alpha            !< parameter determining how the roots grow
real, dimension(:), allocatable :: prcnslai         !< storage/imaginary lai is this percentage of maximum leaf area index that a given root+stem biomass can support
real, dimension(:), allocatable :: minslai          !< minimum storage lai. this is what the model uses as lai when growing
                                                    !!vegetation for scratch. consider these as model seeds.
real, dimension(:), allocatable :: mxrtdpth         !< maximum rooting depth. this is used so that the rooting depths simulated by ctem's variable rooting depth
                                                    !! parameterzation are constrained to realistic values
real, dimension(:), allocatable :: albvis           !< visible albedos of the ctem pfts
real, dimension(:), allocatable :: albnir           !< near IR albedos of the 9 ctem pfts

!! competition_mod.f90 parameters: ------

!! the model basically uses the temperature of the coldest month as
!! the major constraint for pft distribution. a range of the coldest
!! month temperature is prescribed for each pft within which pfts are
!! allowed to exist. in addition for tropical broadleaf drought
!! deciduous trees measure(s) of aridity (function of precipitation
!! and potential evaporation) are used.

! existence subroutine: ----------

real, dimension(:), allocatable :: tcoldmin         !< minimum coldest month temperature
real, dimension(:), allocatable :: tcoldmax         !< maximum coldest month temperature
real, dimension(:), allocatable :: twarmmax         !< maximum warmest month temperature
real, dimension(:), allocatable :: gdd5lmt          !< minimum gdd above 5 c required to exist
real, dimension(:), allocatable :: aridlmt          !< aridity index limit for broadleaf drought/dry deciduous trees
real, dimension(:), allocatable :: dryseasonlmt     !< minimum length of dry season for PFT to exist
real, dimension(:), allocatable :: bio2sap          !< multiplying factor for converting biomass density to sapling density
                                        !! smaller numbers give faster colonization rates.
real :: bioclimrt                       !< mortality rate (1/year) for pfts that no longer exist within their pre-defined bioclimatic range

! ctem.f parameters: ----------

real, dimension(:), allocatable :: grescoef         !< Growth respiration coefficient
real, dimension(:), allocatable :: humicfac         !< Humification factor - used for transferring carbon from litter into soil c pool
real, dimension(:), allocatable :: laimin           !< Minimum lai below which a pft doesn't expand
real, dimension(:), allocatable :: laimax           !< Maximum lai above which a pft always expands and lambdamax fraction of npp is used for expansion
real :: lambdamax                       !< Max. fraction of npp that is allocated for reproduction/colonization
real :: repro_fraction                  !< Fraction of NPP that is used to create reproductive tissues

! disturbance parameters: ------------

real, dimension(2) :: bmasthrs_fire !< min. and max. vegetation biomass thresholds to initiate fire, \f$kg c/m^2\f$
real :: extnmois_veg                !< extinction moisture content for estimating vegetation fire likeliness due to soil moisture
real :: extnmois_duff               !< extinction moisture content for estimating duff layer fire likeliness due to soil moisture
real :: lwrlthrs                    !< lower cloud-to-ground lightning threshold for fire likelihood flashes/km^2.year
                                    ! FireMIP value: 0.025
real :: hgrlthrs                    !< higher cloud-to-ground lightning threshold for fire likelihood flashes/km^2.year
                                    ! FireMIP value: 1.0

!>parameter m (mean) and b of logistic distribution used for \n
!>**Parmlght was increased to 0.8 to make it so areas with higher amounts of
!>lightning have higher lterm. The saturation is still the same, but the
!>increase is more gradual at low lightning density. JM
real :: parmlght
real :: parblght                    !< estimating fire likelihood due to lightning
real :: reparea                     !< typical area representing ctem's fire parameterization (km2)
real :: popdthrshld                 !< threshold of population density (people/km2) [Kloster et al., biogeosci. 2010]
real :: f0                          !< Fire spread rate in the absence of wind
real, dimension(:), allocatable :: maxsprd      !< max. fire spread rate, km/hr
real, dimension(:), allocatable :: frco2glf     !< fraction of green leaf biomass converted to gases due to combustion
real, dimension(:), allocatable :: frco2blf     !< fraction of brown leaf biomass converted to gases due to combustion
real, dimension(:), allocatable :: frltrglf     !< fraction of green leaf biomass becoming litter after combustion
real, dimension(:), allocatable :: frltrblf     !< fraction of brown leaf biomass becoming litter after combustion
real, dimension(:), allocatable :: frco2stm     !< fraction of stem biomass converted to gases due to combustion
real, dimension(:), allocatable :: frltrstm     !< fraction of stem biomass becoming litter after combustion
real, dimension(:), allocatable :: frco2rt      !< fraction of root biomass converted to gases due to combustion
real, dimension(:), allocatable :: frltrrt      !< fraction of root biomass becoming litter after combustion
real, dimension(:), allocatable :: frltrbrn     !< fraction of litter burned during fire and emitted as gases
real, dimension(:), allocatable :: standreplace !< pft prevalence for stand replacing fire events (based on resistance to fire damage, ie. cambial kill)(unitless)
real, dimension(:), allocatable :: emif_co2     !< pft-specific emission factors for CO2,g species / (kg DOM)
real, dimension(:), allocatable :: emif_co      !< pft-specific emission factors for CO,g species / (kg DOM)
real, dimension(:), allocatable :: emif_ch4     !< pft-specific emission factors for CH4,g species / (kg DOM)
real, dimension(:), allocatable :: emif_nmhc    !< pft-specific emission factors for non-methane hydrocarbons,g species / (kg DOM)
real, dimension(:), allocatable :: emif_h2      !< pft-specific emission factors for H2,g species / (kg DOM)
real, dimension(:), allocatable :: emif_nox     !< pft-specific emission factors for NOx,g species / (kg DOM)
real, dimension(:), allocatable :: emif_n2o     !< pft-specific emission factors for N2O,g species / (kg DOM)
real, dimension(:), allocatable :: emif_pm25    !< pft-specific emission factors for particles <2.5 micrometers in diameter,g species / (kg DOM)
real, dimension(:), allocatable :: emif_tpm     !< pft-specific emission factors for total particulate matter,g species / (kg DOM)
real, dimension(:), allocatable :: emif_tc      !< pft-specific emission factors for total carbon,g species / (kg DOM)
real, dimension(:), allocatable :: emif_oc      !< pft-specific emission factors for organic carbon,g species / (kg DOM)
real, dimension(:), allocatable :: emif_bc      !< pft-specific emission factors for black carbon,g species / (kg DOM)

! hetres parameters: ----------

real, dimension(:), allocatable :: bsratelt     !< litter respiration rates at 15 c in in kg c/kg c.year
real, dimension(:), allocatable :: bsratesc     !< soil carbon respiration rates at 15 c in kg c/kg c.year
real, dimension(4) :: tanhq10       !< Constants used in tanh formulation of respiration Q10 determination
real :: alpha_hetres                !< parameter for finding litter temperature as a weighted average of top soil layer temperature and root temperature
real :: bsratelt_g                  !< bare ground litter respiration rate at 15 c in kg c/kg c.year
real :: bsratesc_g                  !< bare ground soil c respiration rates at 15 c in kg c/kg c.year
real :: a                           !< parameter describing exponential soil carbon profile. used for estimating temperature of the carbon pool

! landuse_change_mod.f90 parameters: --------------

real, dimension(3) :: combust       !< how much deforested/chopped off biomass is combusted (these absolutely must add to 1.0!)
real, dimension(3) :: paper         !< how much deforested/chopped off biomass goes into short term storage such as paper
real, dimension(3) :: furniture     !< how much deforested/chopped off biomass goes into long term storage such as furniture
real, dimension(2) :: bmasthrs      !< biomass thresholds for determining if deforested area is a forest, a shrubland, or a bush kg c/m2

! mainres.f parameters: ----------

!    Base respiration rates for stem and root for ctem pfts in kg c/kg c.year at 15 degrees celcius. note that maintenance
!    respiration rates of root are higher because they contain both wood (coarse roots) and fine roots.
!    New parameter values introduced to produce carbon use efficiencies more in
!    line with literature (zhang et al. 2009, luyssaert et al. gcb 2007)
!    values changed for bsrtstem and bsrtroot. jm 06.2012

real, dimension(:), allocatable :: bsrtstem     !< Base respiration rates for stem in kg c/kg c.year at 15 degrees celcius (values differ if using prescribed vs. competition run)
real, dimension(:), allocatable :: bsrtroot     !< Base respiration rates for root in kg c/kg c.year at 15 degrees celcius (values differ if using prescribed vs. competition run)
real :: minlvfr                     !< Minimum live wood fraction

! peatlands_mod.f90 parameters:

! mosspht subroutine parameters

real :: rmlmoss25                   !< Base dark respiration rate of moss umol CO2/m2/s
real :: tau25m                      !< tau coefficient (rate at 25 Celcius) iscomputed on the basis of the specificity
                                    !! factor (102.33) times Kco2/Kh2o (28.38) to convert for value in solution to that
                                    !! based in air.
real :: ektau                       !< \f$J mol^{-1}\f$  !FLAG this is???
real :: kc25                        !< kinetic coef for CO2 at 25 C (Pa)
real :: ko25                        !< kinetic coef for O2 at 25C (rate)  (Pa)
real :: ec                          !< Activation energy for Kc of CO2 (\f$J mol^{-1}\f$) at 25C
real :: ej                          !< Activation energy for electron transport, (\f$J mol^{-1}\f$)
real :: eo                          !< Activation energy for Ko (\f$J mol^{-1}\f$) at 25C
real :: evc                         !< Activation energy for carboxylation (\f$J mol^{-1}\f$) at 25C
real :: sj                          !< Constant affecting J at low temperature (\f$J mol^{-1}\f$) at 25C
real :: hj                          !< Constant affecting J at high temperature (\f$J mol^{-1}\f$) at 25C
real :: alpha_moss                  !< Efficiency of conversion of incident photons into electrons(mol electron/mol photon),

! decp subroutine parameters:

real::     dctmin                   !< minimum temperature of soil respiration,K (peatland soils)
real::     dcbaset                  !< base temperature for Q10, K (peatland soils)
real::     bsrateltms               !< heterotrophic respiration base rate for peatlands (yr-1)

! peatland parameters used in other subroutines:

real :: zolnmoss                    !< natual logarithm of roughness length of moss surface.
real :: thpmoss                     !< pore volume of moss (\f$m^3/m^3\f$)
real :: thrmoss                     !< Liquid water retention capacity of moss (\f$m^3/m^3\f$)
real :: thmmoss                     !< residual liquid water content after freezing or evaporation of moss (\f$m^3/m^3\f$)
real :: bmoss                       !< Clapp and Hornberger empirical "b" parameter of moss
real :: psismoss                    !< Soil moisure suction at saturation of moss (m)
real :: grksmoss                    !< Saturated hydrualic conductivity of moss (m)
real :: hcpmoss                     !< Volumetric heat capacity of moss (\f$J m^{-3} K^{-1}\f$)

!real :: sphms = 2.70E3      !< Specific heat of moss layer (\f$J m^{-2} K^{-1}\f$) ! FLAG NOT USED also check units - written as J/kg/K. JM Sep 26 2016
                            !! same as that specific heat for vegetation
                            !! (Berlinger et al 2001)

!real :: rhoms = 40.0        !< Density of dry moss (\f$kg m^{-3}\f$)  ! FLAG NOT USED JM Oct 2016.
                            !! Based upon:
                            !! 40.0 Price et al. 2008, Price and
                            !! Whittington 2010 value for feather moss
                            !! is lower than sphagnmum (20 kg/m3 in
                            !! O'Donnell et al. 2009)

!real :: slams = 20.0        !< Specific leaf area of moss (\f$m^2 kg^{-1}\f$), ! FLAG NOT USED JM Oct 2016.
                            !! Based upon:
                            !! S vensson 1995 sphangum value ranges from
                            !! 13.5~47.3 m2/kg (Lamberty et al. 2007)

real::   grescoefmoss               !< Moss growth respiration coefficient
real::   rmortmoss                  !< Moss mortality rate constant (year-1)
real::   humicfacmoss               !< Humification ratio of moss litter

! phenology.f parameters: ----------

integer, dimension(:), allocatable :: dayschk       !< Number of days over which to check if net photosynthetic rate is positive before initiating leaf onset
real, dimension(:), allocatable :: drgta            !< Parameter determining how fast soil dryness causes leaves to fall
real, dimension(:), allocatable :: eta              !< eta and kappa, parameters for estimating min. stem+root biomass
real, dimension(:), allocatable :: kappa            !< required to support green leaf biomass. kappa is 1.6 for trees and crops, and 1.2 for grasses.
real, dimension(2) :: flhrspan          !< Harvest span (time in days over which crops are harvested,  15 days),
                                        !< and  fall span (time in days over which bdl cold dcd plants shed their leaves,  30 days)
real :: fracbofg                        !< Parameter used to estimate lai of brown leaves. We assume that SLA of brown leaves is this fraction of SLA of green leaves
real, dimension(:), allocatable :: harvthrs         !< LAI threshold for harvesting crops. values are zero for all pftsother than c3 and c4 crops.
real, dimension(:), allocatable :: specsla          !< CTEM can use user-specified specific leaf areas (SLA) if the following specified values are greater than zero
real, dimension(:), allocatable :: thrprcnt         !< Percentage of max. LAI that can be supported which is used as a threshold for determining leaf phenology status
real, dimension(:), allocatable :: lwrthrsh         !< Lower temperature threshold for ctem's 9 pfts. these are used to estimate cold stress related leaf loss rate (degree c)
real, dimension(:), allocatable :: cdlsrtmx         !< Max. loss rate for cold stress for all 9 pfts, (1/day)
real :: kmort1                                      !< kmort1, parameter used in growth efficiency mortality formulation

real, dimension(:), allocatable :: mxmortge !< Maximum mortality when growth efficiency is zero (1/year) (values differ if using prescribed vs. competition run)
real, dimension(:), allocatable :: maxage   !< Maximum plant age. used to calculate intrinsic mortality rate.
                                !< maximum age for crops is set to zero since they will be harvested
                                !< anyway. grasses are treated the same way since the turnover time
                                !< for grass leaves is 1 year and for roots is 2 years. (values differ if using prescribed vs. competition run)

real, dimension(:), allocatable :: lfespany !< Leaf life span (in years) for CTEM's pfts
real, dimension(:), allocatable :: drlsrtmx !< Max. loss rate for drought stress for all 9 pfts, (1/day) (values differ if using prescribed vs. competition run)
real, dimension(:), allocatable :: colda    !< Parameter determining how fast cold temperatures causes leaves to fall
integer, dimension(2) :: coldlmt!< No. of days for which some temperature has to remain below a given threshold for initiating a process, days
real, dimension(2) :: coldthrs  !<1. -5 c threshold for initiating "leaf fall" mode for ndl dcd trees \n
                                !!2.  8 c threshold for initiating "harvest" for crops, the array colddays tracks days corresponding to these thresholds
real :: roothrsh                !< Root temperature threshold for initiating leaf onset for cold broadleaf deciduous pft, degrees celcius

! turnover.f parameters: ---------------------------------

real, dimension(:), allocatable :: stemlife !< Stemlife, turnover time scale for stem for different pfts
real, dimension(:), allocatable :: rootlife !< Rootlife, turnover time scale for root for different pfts
real :: stmhrspn                !< Stem harvest span. same as crop harvest span. period in days over which crops are harvested.

! wetland_methane.f90 parameters: ------------

real :: ratioch4                !< methane to carbon dioxide flux scaling factor.

!>ratio of wetland to upland respiration\n
!>Use the heterotrophic respiration outputs for soil and litter as the ecosystem basis.  These were summed as "hetrores".
!>This respiration is for upland soils; we multiply by wtdryres as the ratio of wetland to upland respiration
!>based on literature measurements: Dalva et al. 1997 found 0.5 factor; Segers 1998 found a 0.4 factor. use 0.45 here (unitless)
real :: wtdryres
real :: lat_thrshld1   !< Northern zone for wetland determination (degrees North)
real :: lat_thrshld2   !< Boundary with southern zone for wetland determination (degrees North)
real :: soilw_thrshN   !< Soil wetness threshold in the North zone
real :: soilw_thrshE   !< Soil wetness threshold in the Equatorial zone
real :: soilw_thrshS   !< Soil wetness threshold in the South zone

! Photosynthesis parameters: --------------------------------------------------

logical , dimension(:), allocatable :: isc4     !< Array telling which vegetation type is c4
real, dimension(:), allocatable :: tlow         !< lower temperature limits for photosynthesis, (kelvin)
real, dimension(:), allocatable :: tup          !< upper temperature limits for photosynthesis, (kelvin)
real, dimension(:), allocatable :: alpha_phtsyn !< quantum efficiencies, values of 0.08 & 0.04 are used for c3 and c4 plants, respectively
real, dimension(:), allocatable :: omega_phtsyn !< leaf scattering coefficients, values of 0.15 & 0.17 are used for c3 and c4 plants, respectively
real, dimension(:), allocatable :: mm           !< parameter m used in photosynthesis-stomatal conductance coupling.
real, dimension(:), allocatable :: bb           !< parameter b used in photosynthesis-stomatal conductance coupling.
real, dimension(:), allocatable :: vpd0         !< parameter vpd0 used in leuning type photosynthesis - stomatal conductance coupling, (Pa)
integer, dimension(:), allocatable :: sn        !< exponent for soil moisture stress. for sn equal to 1, photosynthesis decreases
                                                !! linearly with soil moisture, and of course non-linearly for values higher than 1.
                                                !! when sn is about 10, photosynthesis does not start decreasing until soil moisture
                                                !! is about half way between wilting point and field capacity.
real, dimension(:), allocatable :: smscale      !< additional constrain of soil moisture stress on photosynthesis. this can be used
                                                !! to simulate the effect of irrigation for crops.
real, dimension(:), allocatable :: vmax         !< max. photosynthetic rate, (\f$(mol CO_2 m^{-2} s^{-1}\f$) values are mainly derived from
                                                !!\cite kattge20090c0 which doesn't include c4. also see \cite alton2017-pd
integer :: reqiter                              !< no. of iterations for calculating intercellular co2 concentration
real :: co2imax                                 !< max. intercellular co2 concentration (Pa)
real :: beta1, beta2                            !< photosynthesis coupling or curvature coefficients
real, dimension(:), allocatable :: inico2i      !< parameter to initialize intercellular co2 conc.
real, dimension(:), allocatable :: rmlcoeff     !< leaf maintenance respiration coefficients
real, dimension(:), allocatable :: chi          !< additional parameters for two-leaf model leaf angle distribution
real :: gamma_w                                 !< photosynthesis down regulation parameters equivalent co2 fertilization effect that we want model to yield
real :: gamma_m                                 !< equivalent co2 fertilization effect that model actually gives without any photosynthesis down-regulation

! Passed variables:

character(350)    :: runParamsFile
logical          :: PFTCompetitionSwitch


! --------------------------------------------------------------------------

contains

!>\ingroup ctem_params_prepareGlobalParams
!!@{
!> Initialize and/or read in all global model parameters

subroutine prepareGlobalParams

    implicit none

    !> Prepare CLASS parameters
    CALL CLASSD

    !> Assign iccp1 and icp1. icc and ican are read in from the job options file.
    iccp1 = icc + 1
    icp1 = ican + 1

    !> Allocate the arrays that store the CTEM parameter values
    call allocateParamsCTEM

    !> Initialize the CTEM parameters, this reads them in from a namelist file.
    call readin_params

end subroutine prepareGlobalParams
!!@}
! --------------------------------------------------------------------------

!>\ingroup ctem_params_allocateParamsCTEM
!!@{
!> Allocate the arrays for CTEM params that require it.

subroutine allocateParamsCTEM()

    implicit none

    kk = l2max * icc

    allocate(modelpft(kk),&
            pftlist(kk),&
            vegtype(kk),&
            kn(kk),&
            omega(kk),&
            epsilonl(kk),&
            epsilons(kk),&
            epsilonr(kk),&
            rtsrmin(kk),&
            aldrlfon(kk),&
            caleaf(kk),&
            castem(kk),&
            caroot(kk),&
            abar(kk),&
            avertmas(kk),&
            alpha(kk),&
            prcnslai(kk),&
            minslai(kk),&
            mxrtdpth(kk),&
            albvis(kk),&
            albnir(kk),&
            tcoldmin(kk),&
            tcoldmax(kk),&
            twarmmax(kk),&
            gdd5lmt(kk),&
            aridlmt(kk),&
            dryseasonlmt(kk),&
            bio2sap(kk),&
            grescoef(kk),&
            humicfac(kk),&
            laimin(kk),&
            laimax(kk),&
            maxsprd(kk),&
            frco2glf(kk),&
            frco2blf(kk),&
            frltrglf(kk),&
            frltrblf(kk),&
            frco2stm(kk),&
            frltrstm(kk),&
            frco2rt(kk),&
            frltrrt(kk),&
            frltrbrn(kk),&
            standreplace(kk),&
            emif_co2(kk),&
            emif_co(kk),&
            emif_ch4(kk),&
            emif_nmhc(kk),&
            emif_h2(kk),&
            emif_nox(kk),&
            emif_n2o(kk),&
            emif_pm25(kk),&
            emif_tpm(kk),&
            emif_tc(kk),&
            emif_oc(kk),&
            emif_bc(kk),&
            bsratelt(kk),&
            bsratesc(kk),&
            bsrtstem(kk),&
            bsrtroot(kk),&
            dayschk(kk),&
            drgta(kk),&
            eta(kk),&
            kappa(kk),&
            harvthrs(kk),&
            specsla(kk),&
            thrprcnt(kk),&
            lwrthrsh(kk),&
            cdlsrtmx(kk),&
            mxmortge(kk),&
            maxage(kk),&
            lfespany(kk),&
            drlsrtmx(kk),&
            colda(kk),&
            stemlife(kk),&
            rootlife(kk),&
            tup(kk),&
            tlow(kk),&
            alpha_phtsyn(kk),&
            omega_phtsyn(kk),&
            isc4(kk),&
            mm(kk),&
            bb(kk),&
            vpd0(kk),&
            sn(kk),&
            smscale(kk),&
            vmax(kk),&
            inico2i(kk),&
            chi(kk),&
            !rmlcoeff(kk))
            rmlcoeff(kk),&
            omega_compete(kk),&
            epsilonl_compete(kk),&
            epsilons_compete(kk),&
            epsilonr_compete(kk),&
            bsrtstem_compete(kk),&
            bsrtroot_compete(kk),&
            mxmortge_compete(kk),&
            maxage_compete(kk),&
            drlsrtmx_compete(kk))

end subroutine allocateParamsCTEM
!!@}
! --------------------------------------------------------------------------


!>\ingroup ctem_params_readin_params
!!@{
!> Read in the CTEM parameters from a namelist file. Populate a few parameters
!! based on what was read in.

subroutine readin_params

    implicit none

    !real, dimension(kk):: omega_compete
    !real, dimension(kk):: epsilonl_compete
    !real, dimension(kk):: epsilons_compete
    !real, dimension(kk):: epsilonr_compete
    !real, dimension(kk):: bsrtstem_compete
    !real, dimension(kk):: bsrtroot_compete
    !real, dimension(kk):: mxmortge_compete
    !real, dimension(kk):: maxage_compete
    !real, dimension(kk):: drlsrtmx_compete
    integer :: i,n
    character(8) :: pftkind
    integer :: isumc,k1c,k2c

    namelist /classicparams/ &
        modelpft,&
        vegtype, &
        pftlist,&
        kn, &
        omega,&
        omega_compete,&
        epsilonl,&
        epsilonl_compete,&
        epsilons,&
        epsilons_compete,&
        epsilonr,&
        epsilonr_compete,&
        rtsrmin,&
        aldrlfon,&
        caleaf,&
        castem,&
        caroot,&
        abar,&
        avertmas,&
        alpha,&
        prcnslai,&
        minslai,&
        mxrtdpth,&
        albvis,&
        albnir,&
        tcoldmin,&
        tcoldmax,&
        twarmmax,&
        gdd5lmt,&
        aridlmt,&
        dryseasonlmt,&
        bio2sap,&
        bioclimrt,&
        grescoef,&
        humicfac,&
        laimin,&
        laimax,&
        lambdamax,&
        repro_fraction,&
        bmasthrs_fire,&
        extnmois_veg,&
        extnmois_duff,&
        lwrlthrs,&
        hgrlthrs,&
        parmlght,&
        parblght,&
        reparea,&
        popdthrshld,&
        f0,&
        maxsprd,&
        frco2glf,&
        frco2blf,&
        frltrglf,&
        frltrblf,&
        frco2stm,&
        frltrstm,&
        frco2rt,&
        frltrrt,&
        frltrbrn,&
        standreplace,&
        emif_co2,&
        emif_co,&
        emif_ch4,&
        emif_nmhc,&
        emif_h2,&
        emif_nox,&
        emif_n2o,&
        emif_pm25,&
        emif_tpm,&
        emif_tc,&
        emif_oc,&
        emif_bc,&
        bsratelt,&
        bsratesc,&
        tanhq10,&
        alpha_hetres,&
        bsratelt_g,&
        bsratesc_g,&
        a,&
        combust,&
        paper,&
        furniture,&
        bmasthrs,&
        bsrtstem,&
        bsrtstem_compete,&
        bsrtroot,&
        bsrtroot_compete,&
        minlvfr,&
        rmlmoss25,&
        tau25m,&
        ektau,&
        kc25,&
        ko25,&
        ec,&
        ej,&
        eo,&
        evc,&
        sj,&
        hj,&
        alpha_moss,&
        dctmin,&
        dcbaset,&
        bsrateltms,&
        zolnmoss,&
        thpmoss,&
        thrmoss,&
        thmmoss,&
        bmoss,&
        psismoss,&
        grksmoss,&
        hcpmoss,&
        grescoefmoss,&
        rmortmoss,&
        humicfacmoss,&
        dayschk,&
        drgta,&
        eta,&
        kappa,&
        flhrspan,&
        fracbofg,&
        harvthrs,&
        specsla,&
        thrprcnt,&
        lwrthrsh,&
        cdlsrtmx,&
        kmort1,&
        mxmortge,&
        mxmortge_compete,&
        maxage,&
        maxage_compete,&
        lfespany,&
        drlsrtmx,&
        drlsrtmx_compete,&
        colda,&
        coldlmt,&
        coldthrs,&
        roothrsh,&
        stemlife,&
        rootlife,&
        stmhrspn,&
        ratioch4,&
        wtdryres,&
        lat_thrshld1,&
        lat_thrshld2,&
        soilw_thrshN,&
        soilw_thrshE,&
        soilw_thrshS,&
        kn,&
        tup,&
        tlow,&
        alpha_phtsyn,&
        omega_phtsyn,&
        isc4,&
        mm,&
        bb,&
        vpd0,&
        sn,&
        smscale,&
        vmax,&
        reqiter,&
        co2imax,&
        beta1,&
        beta2,&
        inico2i,&
        chi,&
        rmlcoeff,&
        gamma_w,&
        gamma_m

    ! ----------

    open(10,file=trim(runParamsFile),action='read',status='old')

    read(10,nml = classicparams)

    close(10)

    numcrops    = 0        !< number of crop pfts
    numtreepfts = 0        !< number of tree pfts
    numgrass    = 0        !< number of grass pfts
    numshrubs   = 0        !< number of shrubs pfts
    allocate(grass(icc))
    allocate(crop(icc))
    allocate(nol2pfts(ican))
    grass = .false.
    crop = .false.

    ! Calculate number of level 2 pfts using modelpft
    do i = 1, ican
        isumc = 0
        k1c = (i-1)*l2max + 1
        k2c = k1c + (l2max - 1)
        do n = k1c, k2c
            if(modelpft(n).eq.1) isumc = isumc + 1
        enddo
        nol2pfts(i)=isumc  ! number of level 2 pfts
    end do

    do i = 1, icc
        pftkind=vegtype(i)
        select case(pftkind)
            case('Tree')
                numtreepfts = numtreepfts + 1
            case('Shrub')
                numshrubs = numshrubs + 1
            case('Grass')
                numgrass = numgrass + 1
                grass(i) = .true.
            case('Crop')
                numcrops = numcrops + 1
                crop(i) = .true.
        end select
    end do

    !Overwrite the prescribed vars with the compete ones if competition is on.
    if (PFTCompetitionSwitch) then
        omega = omega_compete
        epsilonl = epsilonl_compete
        epsilons = epsilons_compete
        epsilonr = epsilonr_compete
        bsrtstem = bsrtstem_compete
        bsrtroot = bsrtroot_compete
        mxmortge = mxmortge_compete
        maxage = maxage_compete
        drlsrtmx = drlsrtmx_compete
    end if

end subroutine readin_params
!!@}

!>\file
!> This module holds CTEM globally accessible parameters
!!
!> These parameters are used in all CTEM subroutines
!> via use statements pointing to this module.
!>The structure of this subroutine is variables that are common to competition/prescribe PFT fractions
!>first, then the remaining variables are assigned different variables if competition is on, or not.
!>
!> PFT parameters
!!
!!Note the structure of vectors which clearly shows the CLASS
!!PFTs (along rows) and CTEM sub-PFTs (along columns)
!!
!!\f[
!!\begin{tabular} { | l | c | c | c | c | c | }
!!\hline
!!needle leaf &  evg &      dcd &      --- & ---& ---\\ \hline
!!broad leaf  &  evg &  dcd-cld &  dcd-dry & EVG-shrubs & DCD-shrubs\\ \hline
!!crops       &   c3 &       c4 &      --- & ---& ---\\ \hline
!!grasses     &   c3 &       c4 &      sedges & ---& ---\\ \hline
!!\end{tabular}
!!\f]
!!
!! We introduced three new
!!PFTs for peatlands: evergreen shrubs, deciduous shrubs, and sedges. Evergreen
!!shrubs, for example the ericaceous shrubs, are the common dominant vascular
!!plants in bogs and poor fens while deciduous shrubs, such as the betulaceous
!!shrubs, often dominate rich fens. Both shrubs are categorized as broadleaf
!!trees in CLASS morphologically, but their phenological and physiological
!!characteristics are more similar to those of needleleaf trees. The shrub
!!tundra ecosystem is situated adjacent to needleleaf forest in the Northern
!!Hemisphere (Kaplan et al., 2003) and they share similar responses to climate
!!in ESMs (e.g. Bonan et al., 2002). Table~2 lists the key parameters for the
!!peatland PFTs used in this model. (The photosynthesis and autotrophic
!!respiration of vascular PFTs are modelled the same as in the original CTEM.)


end module ctem_params
