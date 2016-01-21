module ctem_params

! This module holds CTEM globally accessible parameters
! These parameters are used in all CTEM subroutines
! via use statements pointing to this module EXCEPT PHTSYN3.f 
! which has the information passed in via arguments. This is a legacy thing.

! Remember that changes to this module will usually only take effect
! after you have done a 'make clean' then a 'make' (because it is a module).

! This module is structured with general parameters first then the PFT specific
! parameters later.

! J. Melton
! Jun 23 2013

! Change History:

! Jun 11 2015 - JM - Add in new disturb params (extnmois_duff and extnmois_veg) replacing duff_dry and extnmois.
! Mar 12 2014 - JM - Allow for two sets of paramters (competition and prescribed PFT fractions).
!                    all ctem subroutines except PHTSYN3 keep their parameters here.   
!
! Jan 17 2014 - JM - Add in more parameters from ctem.f, phenology.f, and allocate.f. 

implicit none

public :: initpftpars

! Constants

real, parameter :: zero = 1.0e-20
real, parameter :: abszero=1e-12 ! this one is for runclassctem.f and allocate.f 

real, parameter :: pi =3.1415926535898d0
real, parameter :: earthrad = 6371.22   ! radius of earth, km
real, parameter :: km2tom2 = 1.0e+06 ! changes from km2 to m2
real, parameter :: deltat =1.0  !CTEM's time step in days

integer, parameter, dimension(12) :: monthdays = [ 31,28,31,30,31,30,31,31,30,31,30,31 ] ! days in each month
integer, parameter, dimension(13) :: monthend = [ 0,31,59,90,120,151,181,212,243,273,304,334,365 ] ! calender day at end of each month
integer, parameter, dimension(12) :: mmday= [ 16,46,75,106,136,167,197,228,259,289,320,350 ] !mid-month day

integer, parameter :: lon = 128 ! specify gcm resolution for longitude
integer, parameter :: lat = 64  ! specify gcm resolution for latitude

! latitudes of the edges of the gcm grid cells for 128/x64 resolution
real, parameter, dimension(lat+1) :: edgelat = &
                                    [ -90.0,-86.4802,-83.7047,-80.9193,-78.1313,-75.3422,-72.5527, &
                                    -69.7628,-66.9727,-64.1825,-61.3922,-58.6018,-55.8114,-53.021, &
                                    -50.2305,-47.44,-44.6495,-41.8589,-39.0684,-36.2778,-33.4872, &
                                    -30.6967,-27.9061,-25.1155,-22.3249,-19.5343,-16.7437,-13.9531, &
                                    -11.1625,-8.3718,-5.5812,-2.7906,0.0,2.7906,5.5812,8.3718,11.16245, &
                                    13.9531,16.7437,19.5343,22.3249,25.1155,27.9061,30.69665,33.4872, &
                                    36.2778,39.06835,41.8589,44.64945,47.43995,50.23045,53.02095, &
                                    55.8114,58.6018,61.3922,64.1825,66.9727,69.7628,72.55265,75.3422, &
                                    78.13125,80.91925,83.7047,86.48015,90.0 ]

! ----
! Model state
integer, parameter :: nlat=1            !
integer, parameter :: nmos=10           ! Number of mosaic tiles
integer, parameter :: ilg = nlat*nmos   !
integer, parameter :: nmon=12           ! Number of months in a year

! ----
! Plant-related

integer, parameter :: ican=4            ! Number of CLASS pfts
integer, parameter :: ignd=26            ! Number of soil layers
integer, parameter :: icp1 = ican + 1   !
integer,parameter  :: icc=9             ! Number of CTEM pfts
integer,parameter  :: iccp1 = icc + 1   !
integer, parameter :: l2max = 3         !
integer, parameter :: kk = 12           ! product of class pfts and l2max
integer, parameter :: numcrops = 2      ! number of crop pfts
integer, parameter :: numtreepfts = 5   ! number of tree pfts
integer, parameter :: numgrass = 2      ! number of grass pfts

integer, parameter :: nbs = 4           !

! Separation of pfts into level 1 (for class) and level 2 (for ctem) pfts.
integer, parameter, dimension(kk) :: modelpft= [ 1,     1,     0,  &  ! CLASS PFT 1 NDL
                             !                  EVG    DCD
                                                 1,     1,     1,  &  ! CLASS PFT 2 BDL
                             !                  EVG  DCD-CLD DCD-DRY  ! NOTE 2 TYPES OF BDL DCD - COLD & DRY
                                                 1,     1,     0,  &  ! CLASS PFT 3 CROP
                             !                  C3      C4
                                                 1,     1,     0 ]   ! CLASS PFT 4 GRASS
                             !                  C3      C4


real, parameter :: seed = 0.001 ! seed pft fraction, same as in competition
                                ! in mosaic mode, all tiles are given this
                                ! as a minimum

real, parameter :: minbare = 1.0e-5 ! minimum bare fraction when running competition on to prevent numerical problems.

! conversion factor from carbon to dry organic matter value is from Li et al. 2012 biogeosci
real, parameter :: c2dom = 450.0 ! gc / kg dry organic matter

real, parameter :: wtCH4 = 16.044        ! Molar mass of CH4 ($g mol^{-1}$)

! simple crop matrix, define the number and position of the crops (NOTE: dimension icc)
logical, parameter, dimension(icc) :: crop = [ .false.,.false.,.false.,.false.,.false.,.true.,.true.,.false.,.false. ]   

! simple grass matric, define the number and position of grass 
logical, parameter, dimension(icc) :: grass = [ .false.,.false.,.false.,.false.,.false.,.false.,.false.,.true.,.true. ] 

integer, parameter, dimension(numgrass) :: grass_ind = [ 8, 9 ] ! index of the grass pfts (only 2 grass pfts at present)        

! ==========================================================================================================================

! PFT specific parameters:

! Parameters used in more than one subroutine:

real :: tolrance ! our tolerance for balancing c budget in kg c/m2 in one day (differs when competition on or not)

real, dimension(kk) :: kn               ! Canopy light/nitrogen extinction coefficient

! allocate.f parameters: ---------------------------------

real, dimension(kk) :: omega            ! omega, parameter used in allocation formulae
real, dimension(kk) :: epsilonl         ! Epsilon leaf, parameter used in allocation formulae
real, dimension(kk) :: epsilons         ! Epsilon stem, parameter used in allocation formulae
real, dimension(kk) :: epsilonr         ! Epsilon root, parameter used in allocation formulae
logical :: consallo                     ! Logical switch for using constant allocation factors (default value is false)
real, dimension(kk) :: rtsrmin          ! Minimum root:shoot ratio mostly for support and stability
real, dimension(kk) :: aldrlfon         ! Allocation to leaves during leaf onset

!               Constant allocation fractions if not using dynamic allocation. (NOT thoroughly tested, and using dynamic allocation is preferable)
real, dimension(kk) :: caleaf 
real, dimension(kk) :: castem 
real, dimension(kk) :: caroot

! bio2str.f parameters: ---------

real, dimension(kk) :: abar             ! parameter determining average root profile
real, dimension(kk) :: avertmas         ! average root biomass (kg c/m2) for ctem's 8 pfts used for estimating rooting profile
real, dimension(kk) :: alpha            ! parameter determining how the roots grow
real, dimension(kk) :: prcnslai         ! storage/imaginary lai is this percentage of maximum leaf area index that a given root+stem biomass can support 
real, dimension(kk) :: minslai          ! minimum storage lai. this is what the model uses as lai when growing vegetation for scratch. consider these as model seeds.
real, dimension(kk) :: mxrtdpth         ! maximum rooting depth. this is used so that the rooting depths simulated by ctem's variable rooting depth parameterzation are
                                        !      constrained to realistic values visible and near ir albedos of the 9 ctem pfts
real, dimension(kk) :: albvis           ! visible albedos of the 9 ctem pfts
real, dimension(kk) :: albnir           ! near IR albedos of the 9 ctem pfts

! competition_mod.f90 parameters: ------

! the model basically uses the temperature of the coldest month as
! the major constraint for pft distribution. a range of the coldest
! month temperature is prescribed for each pft within which pfts are
! allowed to exist. in addition for tropical broadleaf drought 
! deciduous trees measure(s) of aridity (function of precipitation
! and potential evaporation) are used.

   ! existence subroutine:

real, dimension(kk) :: tcoldmin         ! minimum coldest month temperature
real, dimension(kk) :: tcoldmax         ! maximum coldest month temperature
real, dimension(kk) :: twarmmax         ! maximum warmest month temperature
real, dimension(kk) :: gdd5lmt          ! minimum gdd above 5 c required to exist
real, dimension(kk) :: aridlmt          ! aridity index limit for broadleaf drought/dry deciduous trees
real, dimension(kk) :: dryseasonlmt     ! minimum length of dry season for PFT to exist

   ! competition subroutine:

real, dimension(kk) :: bio2sap          ! multiplying factor for converting biomass density to sapling density
real :: bioclimrt                       ! mortality rate (1/year) for pfts that no longer exist within their pre-defined bioclimatic range

! ctem.f parameters: ----------

real, dimension(kk) :: grescoef         ! Growth respiration coefficient 
real, dimension(kk) :: humicfac         ! Humification factor - used for transferring carbon from litter into soil c pool
real, dimension(kk) :: laimin           ! Minimum lai below which a pft doesn't expand
real, dimension(kk) :: laimax           ! Maximum lai above which a pft always expands and lambdamax fraction of npp is used for expansion
real :: lambdamax                       ! Max. fraction of npp that is allocated for reproduction/colonization
real :: repro_fraction                  ! Fraction of NPP that is used to create reproductive tissues

! disturb.f90 parameters: -------------

real, dimension(2) :: bmasthrs_fire     ! min. and max. vegetation biomass thresholds to initiate fire, kg c/m^2
real :: extnmois_veg                    ! extinction moisture content for estimating vegetation fire likeliness due to soil moisture 
real :: extnmois_duff                   ! extinction moisture content for estimating duff layer fire likeliness due to soil moisture 
real :: lwrlthrs                        ! lower cloud-to-ground lightning threshold for fire likelihood flashes/km^2.year
real :: hgrlthrs                        ! higher cloud-to-ground lightning threshold for fire likelihood flashes/km^2.year
real :: parmlght                        ! parameter m (mean) and b of logistic distribution used for 
real :: parblght                        !     estimating fire likelihood due to lightning
real :: reparea                         ! typical area representing ctem's fire parameterization (km2)
real :: popdthrshld                     ! threshold of population density (people/km2) [Kloster et al., biogeosci. 2010]

!     **These values were being calculated each time, they shouldn't be as they
!     are just parameters. JM
!      real, parameter :: ymin = 0.01798621     !=1.0/( 1.0+exp((parmlght-0.0)/parblght) )
!      real, parameter :: ymax = 0.9975273768   !=1.0/( 1.0+exp((parmlght-1.0)/parblght) )
!      real, parameter :: slope = 0.0204588331  !=abs(0.0-ymin)+abs(1.0-ymax)
!      real :: ymin, ymax, slope

real :: alpha_fire                      ! parameter alpha_fire and f0 used for estimating wind function for fire spread rate
real :: f0                              ! Fire spread rate in the absence of wind  -- Not used in CTEM v 2.0!
real, dimension(kk) :: standreplace     ! pft prevalence for stand replacing fire events (based on resistance to fire damage, ie. cambial kill)(unitless)
real, dimension(kk) :: maxsprd          ! max. fire spread rate, km/hr
real, dimension(kk) :: frco2glf         ! fraction of green leaf biomass converted to gases due to combustion
real, dimension(kk) :: frco2blf         ! fraction of brown leaf biomass converted to gases due to combustion
real, dimension(kk) :: frltrglf         ! fraction of green leaf biomass becoming litter after combustion
real, dimension(kk) :: frltrblf         ! fraction of brown leaf biomass becoming litter after combustion
real, dimension(kk) :: frco2stm         ! fraction of stem biomass converted to gases due to combustion
real, dimension(kk) :: frltrstm         ! fraction of stem biomass becoming litter after combustion
real, dimension(kk) :: frco2rt          ! fraction of root biomass converted to gases due to combustion
real, dimension(kk) :: frltrrt          ! fraction of root biomass becoming litter after combustion
real, dimension(kk) :: frltrbrn         ! fraction of litter burned during fire and emitted as gases

!     emissions factors by chemical species
real, dimension(kk) :: emif_co2         ! pft-specific emission factors for CO2,g species / (kg DOM)
real, dimension(kk) :: emif_co          ! pft-specific emission factors for CO,g species / (kg DOM)
real, dimension(kk) :: emif_ch4         ! pft-specific emission factors for CH4,g species / (kg DOM)
real, dimension(kk) :: emif_nmhc        ! pft-specific emission factors for non-methane hydrocarbons,g species / (kg DOM)
real, dimension(kk) :: emif_h2          ! pft-specific emission factors for H2,g species / (kg DOM)
real, dimension(kk) :: emif_nox         ! pft-specific emission factors for NOx,g species / (kg DOM)
real, dimension(kk) :: emif_n2o         ! pft-specific emission factors for N2O,g species / (kg DOM)

!     emission factors for aerosols
real, dimension(kk) :: emif_pm25        ! pft-specific emission factors for particles <2.5 micrometers in diameter,g species / (kg DOM)
real, dimension(kk) :: emif_tpm         ! pft-specific emission factors for total particulate matter,g species / (kg DOM)
real, dimension(kk) :: emif_tc          ! pft-specific emission factors for total carbon,g species / (kg DOM)
real, dimension(kk) :: emif_oc          ! pft-specific emission factors for organic carbon,g species / (kg DOM)
real, dimension(kk) :: emif_bc          ! pft-specific emission factors for black carbon,g species / (kg DOM)

! hetres parameters: -------------

real, dimension(kk) :: bsratelt         ! litter respiration rates at 15 c in in kg c/kg c.year 
real, dimension(kk) :: bsratesc         ! soil carbon respiration rates at 15 c in kg c/kg c.year
real, dimension(4) :: tanhq10           ! Constants used in tanh formulation of respiration Q10 determination
real :: alpha_hetres                    ! parameter for finding litter temperature as a weighted average of top soil layer temperature and root temperature
real :: bsratelt_g                      ! bare ground litter respiration rate at 15 c in kg c/kg c.year
real :: bsratesc_g                      ! bare ground soil c respiration rates at 15 c in kg c/kg c.year
real :: a                               ! parameter describing exponential soil carbon profile. used for estimating temperature of the carbon pool

! landuse_change_mod.f90: -------------

real, dimension(3) :: combust           ! how much deforested/chopped off biomass is combusted (these absolutely must add to 1.0!)
real, dimension(3) :: paper             ! how much deforested/chopped off biomass goes into short term storage such as paper
real, dimension(3) :: furniture         ! how much deforested/chopped off biomass goes into long term storage such as furniture
real, dimension(2) :: bmasthrs          ! biomass thresholds for determining if deforested area is a forest, a shrubland, or a bush kg c/m2

real :: tolrnce1                        ! kg c, tolerance of total c balance 

! mainres.f parameters: ----------

!    Base respiration rates for stem and root for ctem pfts in kg c/kg c.year at 15 degrees celcius. note that maintenance
!    respiration rates of root are higher because they contain both wood (coarse roots) and fine roots.
!    New parameter values introduced to produce carbon use efficiencies more in
!    line with literature (zhang et al. 2009, luyssaert et al. gcb 2007)
!    values changed for bsrtstem and bsrtroot. jm 06.2012

real, dimension(kk) :: bsrtstem 
real, dimension(kk) :: bsrtroot 
real :: minlvfr                         ! Minimum live wood fraction


! mortality.f parameters: ----------

real :: kmort1                          ! kmort1, parameter used in growth efficiency mortality formulation
real, dimension(kk) :: mxmortge         ! Maximum mortality when growth efficiency is zero (1/year)
real, dimension(kk) :: maxage           ! Maximum plant age. used to calculate intrinsic mortality rate.
                                        !     maximum age for crops is set to zero since they will be harvested
                                        !     anyway. grasses are treated the same way since the turnover time
                                        !     for grass leaves is ~1 year and for roots is ~2 year. 

! phenology.f parameters: ----------

real, dimension(kk) :: eta              ! eta and kappa, parameters for estimating min. stem+root biomass
real, dimension(kk) :: kappa            ! required to support green leaf biomass. kappa is 1.6 for trees and crops, and 1.2 for grasses.
real, dimension(kk) :: lfespany         ! Leaf life span (in years) for CTEM's 9 pfts
real, dimension(kk) :: specsla          ! CTEM can use user-specified specific leaf areas (SLA) if the following specified values are greater than zero
real :: fracbofg                        ! Parameter used to estimate lai of brown leaves. We assume that SLA of brown leaves is this fraction of SLA of green leaves
real, dimension(kk) :: cdlsrtmx         ! Max. loss rate for cold stress for all 9 pfts, (1/day)
real, dimension(kk) :: drlsrtmx         ! Max. loss rate for drought stress for all 9 pfts, (1/day)
real, dimension(kk) :: drgta            ! Parameter determining how fast soil dryness causes leaves to fall
real, dimension(kk) :: colda            ! Parameter determining how fast cold temperatures causes leaves to fall
real, dimension(kk) :: lwrthrsh         ! Lower temperature threshold for ctem's 9 pfts. these are used to estimate cold stress related leaf loss rate (degree c)
integer, dimension(kk) :: dayschk       ! Number of days over which to check if net photosynthetic rate is positive before initiating leaf onset
integer, dimension(2) :: coldlmt        ! No. of days for which some temperature has to remain below a given threshold for initiating a process
real, dimension(2) :: coldthrs          !     1. -5 c threshold for initiating "leaf fall" mode for ndl dcd trees
                                        !     2.  8 c threshold for initiating "harvest" mode for crops the array colddays(i,2)
                                        !     tracks days corresponding to these thresholds
real, dimension(kk) :: harvthrs         ! LAI threshold for harvesting crops. values are zero for all pftsother than c3 and c4 crops.
real, dimension(2) :: flhrspan          ! Harvest span (time in days over which crops are harvested, ~15 days), 
                                        !     and  fall span (time in days over which bdl cold dcd plants shed their leaves, ~30 days)
real, dimension(kk) :: thrprcnt         ! Percentage of max. LAI that can be supported which is used as a threshold for determining leaf phenology status
real :: roothrsh                        ! Root temperature threshold for initiating leaf onset for cold broadleaf deciduous pft, degrees celcius

! turnover.f parameters: ---------------------------------

real, dimension(kk) :: stemlife         ! Stemlife, turnover time scale for stem for different pfts
real, dimension(kk) :: rootlife         ! Rootlife, turnover time scale for root for different pfts
real :: stmhrspn                        ! Stem harvest span. same as crop harvest span. period in days over which crops are harvested. 

! wetland_methane.f90 parameters: ------------

real :: ratioch4                        ! methane to carbon dioxide flux scaling factor.
real :: wtdryres                        ! ratio of wetland to upland respiration 
real :: factor2                         ! constant value for secondary (ch4wet2) methane emissions calculation
real :: lat_thrshld1                    ! Northern zone for wetland determination (degrees North)
real :: lat_thrshld2                    ! Boundary with southern zone for wetland determination (degrees North)
real :: soilw_thrshN                    ! Soil wetness threshold in the North zone
real :: soilw_thrshE                    ! Soil wetness threshold in the Equatorial zone
real :: soilw_thrshS                    ! Soil wetness threshold in the South zone

! ----=====-------=========-----------========---------=========--------========-----------==========---------=======---========**
contains

subroutine initpftpars(compete)

! PFT parameters

!     Note the structure of vectors which clearly shows the CLASS
!     PFTs (along rows) and CTEM sub-PFTs (along columns)
!
!     needle leaf |  evg       dcd       ---
!     broad leaf  |  evg   dcd-cld   dcd-dry
!     crops       |   c3        c4       ---
!     grasses     |   c3        c4       ---


implicit none

logical, intent(in) :: compete   !true if the competition scheme is on.

! The structure of this subroutine is variables that are common to competition/prescribe PFT fractions
! first, then the remaining variables are assigned different variables if competition is on, or not.

!   ********************************************************************************************
!   =============                                                     ==========================
!   =============                   PARAMETERS COMMON TO BOTH         ==========================
!   ============   COMPETITION-ON AND PRESCRIBED PFT FRACTIONAL COVER MODEL MODES    ===========
!   =============                                                     ==========================
!   ********************************************************************************************

! Parameters used in more than one subroutine:

tolrance = 0.0001d0 

! (kn -> CAREFUL: Separate set defined in PHTSYN3.f!!)
kn= [ 0.50, 0.50, 0.00, &
      0.50, 0.50, 0.50, &
      0.40, 0.48, 0.00, &
      0.46, 0.44, 0.00 ]
      
lfespany  =   [ 5.00, 1.00, 0.00, &
                1.50, 1.00, 1.00, &  !PFT 3 was 1.75 (from IBIS), 2.00 follows LPJ. JM Mar 2014.
!lfespany  =   [ 5.00, 0.40, 0.00, & !These values are for the tests with Ray and Bakr. JM Mar 2015. HAVE NOT TESTED WITH COMPETITION!!!!
!               1.50, 1.00, 1.00, &
                1.75, 1.75, 0.00, &
                1.00, 1.00, 0.00 ]      

! allocate.f parameters: --------------

consallo = .false.

rtsrmin  = [ 0.16, 0.16, 0.00, &
             0.16, 0.16, 0.32, &
             0.16, 0.16, 0.00, & 
             0.50, 0.50, 0.00 ]

aldrlfon = [ 1.00, 1.00, 0.00, &
             1.00, 1.00, 1.00, &
             1.00, 1.00, 0.00, &
             1.00, 1.00, 0.00 ]

caleaf = [ 0.275, 0.300, 0.000, &
           0.200, 0.250, 0.250, &
           0.400, 0.400, 0.000, &
           0.450, 0.450, 0.000 ]

castem = [ 0.475, 0.450, 0.000, &
           0.370, 0.400, 0.400, &
           0.150, 0.150, 0.000, &
           0.000, 0.000, 0.000 ]

caroot =  [ 0.250, 0.250, 0.000, &
            0.430, 0.350, 0.350, &
            0.450, 0.450, 0.000, & 
            0.550, 0.550, 0.000 ]

! bio2str.f parameters: --------

abar    = [ 4.70, 5.86, 0.00, &
            3.87, 3.46, 3.97, &
            3.97, 3.97, 0.00, &
            5.86, 4.92, 0.00 ]

avertmas = [ 1.85, 1.45, 0.00, &
             2.45, 2.10, 2.10, &
             0.10, 0.10, 0.00, &
             0.70, 0.70, 0.00 ]

alpha   = [ 0.80, 0.80, 0.00, &
            0.80, 0.80, 0.80, &
            0.80, 0.80, 0.00, &
            0.80, 0.80, 0.00 ]

prcnslai = [ 7.5, 7.5, 0.0, &
             7.5, 7.5, 7.5, &
             7.5, 7.5, 0.0, &
             2.5, 2.5, 0.0 ]

minslai = [ 0.3, 0.3, 0.0, &
            0.3, 0.3, 0.3, &
            0.2, 0.2, 0.0, &
            0.2, 0.2, 0.0 ]

mxrtdpth = [ 3.00, 3.00, 0.00, &
             5.00, 5.00, 3.00, &
             2.00, 2.00, 0.00, &
             1.00, 1.00, 0.00 ] 

albvis = [ 3.00, 3.00, 0.00, &
           3.00, 5.00, 5.00, &
           5.50, 5.50, 0.00, &
           5.00, 6.00, 0.00 ]

albnir = [ 19.0, 19.0, 0.00, &
           23.0, 29.0, 29.0, &
           34.0, 34.0, 0.00, &
           30.0, 34.0, 0.00 ] 

! ctem.f parameters: ----------

humicfac = [ 0.42, 0.42, 0.00, &
             0.53, 0.48, 0.48, &
             0.10, 0.10, 0.00, &
             0.42, 0.42, 0.00 ]

grescoef = [ 0.15, 0.15, 0.00, &
             0.15, 0.15, 0.15, &
             0.15, 0.15, 0.00, &
             0.15, 0.15, 0.00 ]

lambdamax = 0.10 

repro_fraction = 0.10

! disturbance parameters: ------------

reparea = 300.0  !500.0 ;  FireMIP value: 300.0

popdthrshld = 300.

f0 = 0.05

bmasthrs_fire = [ 0.2, 1.0 ] 

extnmois_veg = 0.3 !0.37 !flag test was 0.3

extnmois_duff = 0.5 !0.55 !flag test was 0.5

lwrlthrs = 0.25 ! FireMIP value: 0.025

hgrlthrs = 10.0 ! FireMIP value: 1.0

!     **Parmlght was increased to 0.8 to make it so areas with higher amounts of
!     lightning have higher lterm. The saturation is still the same, but the 
!     increase is more gradual at low lightning density. JM
parmlght = 0.8  
parblght = 0.1

! Li et al vals (except grass)
maxsprd = [  0.54, 0.54, 0.00, &
             0.40, 0.40, 0.40, &
             0.00, 0.00, 0.00, &
             0.72, 0.72, 0.00 ]

frco2glf = [ 0.21, 0.21, 0.00, & 
             0.21, 0.21, 0.21, & 
             0.00, 0.00, 0.00, &
             0.24, 0.24, 0.00 ]   

frco2blf = [ 0.00, 0.00, 0.00, & 
             0.00, 0.00, 0.00, & 
             0.00, 0.00, 0.00, &
             0.27, 0.27, 0.00 ]

frltrglf = [ 0.06, 0.06, 0.00, &
             0.06, 0.06, 0.06, & 
             0.00, 0.00, 0.00, &
             0.03, 0.03, 0.00 ]  

frltrblf = [ 0.00, 0.00, 0.00, &
             0.00, 0.00, 0.00, & 
             0.00, 0.00, 0.00, &
             0.02, 0.02, 0.00 ]

frco2stm = [ 0.06, 0.06, 0.00, &
             0.06, 0.03, 0.03, &  
             0.00, 0.00, 0.00, &
             0.00, 0.00, 0.00 ]

frltrstm = [ 0.15, 0.15, 0.00, &
             0.15, 0.09, 0.09, & 
             0.00, 0.00, 0.00, &
             0.00, 0.00, 0.00 ]

frco2rt = [ 0.0, 0.0, 0.0, &
            0.0, 0.0, 0.0, &
            0.0, 0.0, 0.0, &
            0.0, 0.0, 0.0 ]

frltrrt = [ 0.03, 0.03, 0.00, &
            0.03, 0.03, 0.03, &
            0.00, 0.00, 0.00, &
            0.08, 0.08, 0.00 ]

frltrbrn = [ 0.15, 0.15, 0.00, &
             0.18, 0.18, 0.18, &
             0.00, 0.00, 0.00, &
             0.21, 0.21, 0.00 ]

!     emissions factors by chemical species
!     
!     Values are from Andreae 2011 as described in Li et al. 2012
!     Biogeosci. Units: g species / (kg DOM)
 
!     Andreae 2011 as described in Li et al. 2012
emif_co2 = [ 1576.0, 1576.0,   0.00, &
             1604.0, 1576.0, 1654.0, &
             1576.0, 1654.0,   0.00, &
             1576.0, 1654.0,   0.00 ]

!    values from Wiedinmyer et al. 2011
!emif_co2 = [ 1514.0, 1514.0,   0.00, &
!             1643.0, 1630.0, 1716.0, &
!             1537.0, 1537.0,   0.00, & 
!             1692.0, 1692.0,   0.00 ]
 
!     Andreae 2011 as described in Li et al. 2012
emif_co = [ 106.0, 106.0, 0.00, &
            103.0, 106.0, 64.0, &
            106.0,  64.0, 0.00, &
            106.0,  64.0, 0.00 ]

!    values from Wiedinmyer et al. 2011
!emif_co = [ 118.0, 118.0, 0.00, &
!             92.0, 102.0, 68.0, &
!            111.0, 111.0, 0.00, &
!             59.0,  59.0, 0.00 ]

!     Andreae 2011 as described in Li et al. 2012
emif_ch4 = [ 4.8, 4.8, 0.0, &
             5.8, 4.8, 2.4, &
             4.8, 2.4, 0.0, &
             4.8, 2.4, 0.0 ]

!    values from Wiedinmyer et al. 2011
!emif_ch4 = [ 6.0, 6.0, 0.0, &
!             5.1, 5.0, 2.6, &
!             6.0, 6.0, 0.0, &
!             1.5, 1.5, 0.0 ]

!     Andreae 2011 as described in Li et al. 2012
emif_nmhc = [ 5.7, 5.7, 0.0, &
              6.4, 5.7, 3.7, &
              5.7, 3.7, 0.0, &
              5.7, 3.7, 0.0 ] 

!    values from Wiedinmyer et al. 2011
!emif_nmhc = [ 5.7, 5.7, 0.0, &
!              1.7, 5.7, 3.4, &
!              7.0, 7.0, 0.0, &
!              3.4, 3.4, 0.0 ]

!     Andreae 2011 as described in Li et al. 2012
emif_h2 = [ 1.80, 1.80, 0.00, &
            2.54, 1.80, 0.98, &
            1.80, 0.98, 0.00, &
            1.80, 0.98, 0.00 ]

!    values from Wiedinmyer et al. 2011
!emif_h2 = [ 2.30, 2.30, 0.00, &
!            3.20, 1.80, 0.97, &
!            2.40, 2.40, 0.00, &
!            0.97, 0.97, 0.00 ]

!     Andreae 2011 as described in Li et al. 2012
emif_nox = [ 3.24, 3.24, 0.00, &
             2.90, 3.24, 2.49, &
             3.24, 2.49, 0.00, &
             3.24, 2.49, 0.00 ]

!    values from Wiedinmyer et al. 2011 (species: "NOx (as NO)" from Table 1)
!emif_nox = [ 1.80, 2.30, 0.00, &
!             2.60, 1.30, 3.90, &
!             3.50, 3.50, 0.00, &
!             2.80, 2.80, 0.00 ]

!     Andreae 2011 as described in Li et al. 2012 
emif_n2o = [ 0.26, 0.26, 0.00, &
             0.23, 0.26, 0.20, &
             0.26, 0.20, 0.00, &
             0.26, 0.20, 0.00 ]

!     emission factors for aerosols

!     Andreae 2011 as described in Li et al. 2012
emif_pm25 = [ 12.7, 12.7, 0.0, &
              10.5, 12.7, 5.2, &
              12.7,  5.2, 0.0, &
              12.7,  5.2, 0.0 ]

!    values from Wiedinmyer et al. 2011
!emif_pm25 = [ 13.0, 13.0, 0.0, &
!               9.7, 13.0, 9.3, &
!               5.8,  5.8, 0.0, &
!               5.4,  5.4, 0.0 ] 

!     Andreae 2011 as described in Li et al. 2012
emif_tpm = [ 17.6, 17.6, 0.0, &
             14.7, 17.6, 8.5, &
             17.6,  8.5, 0.0, &
             17.6,  8.5, 0.0 ]

!    values from Wiedinmyer et al. 2011
!emif_tpm = [ 18.0, 18.0, 0.0, &
!             13.0, 18.0,15.4, &
!             13.0, 13.0, 0.0, &
!              8.3,  8.3, 0.0 ]

!     Andreae 2011 as described in Li et al. 2012
emif_tc = [ 8.3, 8.3, 0.0, &
            7.2, 8.3, 3.4, &
            8.3, 3.4, 0.0, &
            8.3, 3.4, 0.0 ]

!    values from Wiedinmyer et al. 2011 (TPC in Table 1)
!emif_tc = [ 8.3, 8.3, 0.0, &
!            5.2, 9.7, 7.1, &
!            4.0, 4.0, 0.0, &
!            3.0, 3.0, 0.0 ]

!     Andreae 2011 as described in Li et al. 2012
emif_oc = [ 9.1, 9.1, 0.0, &
            6.7, 9.1, 3.2, & 
            9.1, 3.2, 0.0, &
            9.1, 3.2, 0.0 ] 

!    values from Wiedinmyer et al. 2011 
!emif_oc = [ 7.8, 7.8, 0.0, &
!            4.7, 9.2, 6.6, &
!            3.3, 3.3, 0.0, &
!            2.6, 2.6, 0.0 ]

!     Andreae 2011 as described in Li et al. 2012
emif_bc = [ 0.56, 0.56, 0.00, &
            0.56, 0.56, 0.47, &
            0.56, 0.47, 0.00, &
            0.56, 0.47, 0.00 ]

!    values from Wiedinmyer et al. 2011
!emif_bc = [ 0.20, 0.20, 0.00, &
!            0.52, 0.56, 0.50, &
!            0.69, 0.69, 0.00, &
!            0.37, 0.37, 0.00 ]

! hetres parameters: ----------


bsratelt = [ 0.4453, 0.5986, 0.0000, &
             0.6339, 0.7576, 0.6957, &
             0.6000, 0.6000, 0.0000, &
             0.5260, 0.5260, 0.0000 ] 

bsratesc = [ 0.0260, 0.0260, 0.0000, &
             0.0208, 0.0208, 0.0208, &
             0.0350, 0.0350, 0.0000, &
             0.0125, 0.0125, 0.0000 ]
! FLAG test new vals, JM Dec 29 2014.
!bsratesc = [ 0.0208, 0.0208, 0.0000, &
!             0.0166, 0.0166, 0.0166, &
!             0.0280, 0.0280, 0.0000, &
!             0.0100, 0.0100, 0.0000 ]

tanhq10  = [ 1.44, 0.56, 0.075, 46.0 ] 
!tanhq10  = [ 1.2, 0.2, 0.075, 46.0 ] ! was 1.44/0.56/0.075/46 TEST FLAG mar 3 JM 2015
           !   a     b      c     d
           ! q10 = a + b * tanh[ c (d-temperature) ]
           ! when a = 2, b = 0, we get the constant q10 of 2. if b is non
           ! zero then q10 becomes temperature dependent

alpha_hetres = 0.7

bsratelt_g = 0.5605

bsratesc_g = 0.02258
! FLAG test new vals, JM Dec 29 2014.
!bsratesc_g = 0.01806


a = 4.0

! landuse_change_mod.f90 parameters: --------------

! NOTE: combust, paper, furniture ABSOLUTELY MUST add to 1.0!
combust =   [ 0.15, 0.30, 0.45 ] 
paper =     [ 0.70, 0.70, 0.55 ]
furniture = [ 0.15, 0.00, 0.00 ]

bmasthrs = [ 4.0, 1.0 ]

tolrnce1 = 0.50  !FLAG would be good to make this consistent with global tolerance so only one value.

! mainres.f parameters: ---------

minlvfr = 0.05

! mortality.f parameters: ---------

kmort1 = 0.3

! phenology.f parameters: ---------

colda = [ 3.0, 3.0, 0.0, &
          3.0, 3.0, 3.0, &
          3.0, 3.0, 0.0, &
          3.0, 3.0, 0.0 ]

coldlmt = [ 7 , 5 ]   ! days

coldthrs = [ -5.0 , 8.0 ]

dayschk = [ 7, 7, 0, &
            7, 7, 7, &
            7, 7, 0, &
            7, 7, 0 ]

drgta = [ 3.0, 3.0, 0.0, &
          3.0, 3.0, 3.0, &
          3.0, 3.0, 0.0, &
          3.0, 3.0, 0.0 ] 

eta = [ 10.0, 30.8, 0.00, &
        31.0, 50.0, 30.0, &
        7.0,  7.0, 0.00,  &
        3.0,  3.0, 0.00 ]

kappa =[ 1.6, 1.6, 0.0, &
         1.6, 1.6, 1.6, &
         1.6, 1.6, 0.0, &
         1.2, 1.2, 0.0 ]
  
flhrspan = [ 17.0, 45.0 ]

fracbofg = 0.55

harvthrs = [ 0.0, 0.0, 0.0, &
             0.0, 0.0, 0.0, &
             4.5, 3.5, 0.0, &
             0.0, 0.0, 0.0 ]

specsla =[  0.0, 0.0, 0.0, &  ! Not used in general. Only will take effect if these are non-zero.
            0.0, 0.0, 0.0, &
            0.0, 0.0, 0.0, &
            0.0, 0.0, 0.0 ]

thrprcnt = [ 40.0, 40.0,  0.0, &
             40.0, 50.0, 50.0, &
             50.0, 50.0,  0.0, &
             40.0, 40.0,  0.0 ]  
             
lwrthrsh = [ -50.0, -5.0, 0.0, & 
               5.0,  8.0, 5.0, &  
               5.0,  5.0, 0.0, &
               0.1,  5.0, 0.0 ]             

cdlsrtmx = [ 0.10, 0.30, 0.00, &  
             0.30, 0.40, 0.15, &
             0.15, 0.15, 0.00, &
             0.15, 0.15, 0.00 ]

roothrsh = 8.0
             
! turnover.f parameters: --------------

stmhrspn = 17.0

stemlife = [ 86.3, 86.3, 0.00, &  
             80.5, 80.5, 75.8, &  
             20.0, 20.0, 0.00, &
              0.00, 0.00, 0.00 ]

rootlife = [ 13.8,13.2, 0.0, &
             12.7,10.9, 9.8, &    
              3.0, 3.0, 0.0, &
              3.0, 3.0, 0.0 ]

! wetland_methane.f90 parameters: -------

!	Rita Wania's thesis suggests about 0.25, but we get a better agreement 
!	to outputs from the Walter's model if we use 0.16.  Note that this
!	scaling factor likely is temperature dependent, and increases with
!	temperature, but it is difficult to know the function, so leave
!	constant for now ratio is mol ch4 to mol co2

ratioch4 = 0.16

!	Use the heterotrophic respiration outputs for soil and litter 
!	as the ecosystem basis.  These were summed as "hetrores".
!	This respiration is for upland soils; we multiply by 
!	wtdryres as the ratio of wetland to upland respiration 
!	based on literature measurements: Dalva et al. 1997 found 0.5 factor; 
!	Segers 1998 found a 0.4 factor. use 0.45 here (unitless)

wtdryres = 0.45

factor2 = 0.015

lat_thrshld1 = 40.0 ! degrees North
lat_thrshld2 = -20.0 ! degrees North
soilw_thrshN = 0.55   ! Soil wetness threshold in the North zone
soilw_thrshE = 0.80   ! Soil wetness threshold in the Equatorial zone
soilw_thrshS = 0.70   ! Soil wetness threshold in the South zone

!   ********************************************************************************************
!   =============                                                     ==========================
!   =============DYNAMIC PFT FRACTIONAL COVER PARAMETERS=(COMPETITION)==========================
!   =============                                                     ==========================
!   ********************************************************************************************

if (compete) then 

! These parameters are used when competition is on. If you are using
! prescribed PFT fractional cover, then the parameters after this section
! are used. Parameters that are the same in both are above this if loop.

! allocate.f parameters: --------------

! Parameterization values based on comparison mostly with LUYSSAERT, S.et al. CO2 balance
! of boreal, temperate, and tropical forests derived from a global database,
! Glob. Chang. Biol., 13(12), 2509–2537, 2007. and informed by LITTON, et al. Carbon
! allocation in forest ecosystems, Glob. Chang. Biol., 13(10), 2089–2109, 2007. Further
! tuning was performed on these basic values. JM Dec 20 2013.

omega = [ 0.80, 0.50, 0.00, & ! IN BOTH PRESCRIBED AND COMPETE (TWO VERSIONS)
          0.80, 0.45, 0.80, &  
          0.05, 0.05, 0.00, &
          1.00, 1.00, 0.00 ]

epsilonl = [ 0.19, 0.45, 0.00, &  ! IN BOTH PRESCRIBED AND COMPETE (TWO VERSIONS)
             0.39, 0.50, 0.30, &  
             0.80, 0.80, 0.00, &
             0.10, 0.10, 0.00 ]

epsilons = [ 0.40, 0.34, 0.00, & ! IN BOTH PRESCRIBED AND COMPETE (TWO VERSIONS)
             0.21, 0.35, 0.10, & 
             0.15, 0.15, 0.00, &
             0.00, 0.00, 0.00 ]

epsilonr = [ 0.41, 0.21, 0.00, &  ! IN BOTH PRESCRIBED AND COMPETE (TWO VERSIONS)
             0.40, 0.15, 0.60, &  
             0.05, 0.05, 0.00, &
             0.90, 0.90, 0.00 ]

! competition_mod.f90 parameters

  ! existence subroutine:

tcoldmin = [-999.9, -999.9,   0.0, & 
               2.5,  -35.0,   4.0, & 
            -999.9, -999.9,   0.0, &
            -999.9, -999.9,   0.0 ]  

tcoldmax = [ 18.0,  -28.0,   0.0, & 
            999.9,   16.0, 900.0, &       
            999.9,  999.9,   0.0, &
            999.9,  999.9,   0.0 ]  

twarmmax = [ 99.9,  25.0,  0.0, & 
             99.9,  99.9, 99.9, &       
             99.9,  99.9,  0.0, &
             99.9,  99.9,  0.0 ]

gdd5lmt = [ 375.0,  600.0,  0.0, &  
           1200.0,  300.0,  9.9, & 
              9.9,    9.9,  0.0, &
              9.9,    9.9,  0.0 ]

aridlmt = [ 9.9,  9.9,  0.0, &
            9.9,  9.9,  0.9, & 
            9.9,  9.9,  0.0, &
            9.9,  9.9,  0.0 ]

dryseasonlmt=[ 99.0,  99.9,    0.0, &
               99.9,  99.9,    5.5, & 
               99.9,  99.9,    0.0, &
               99.9,  99.9,    0.0 ]
 
  ! competition subroutine

  ! smaller numbers give faster colonization rates.
bio2sap = [ 0.32, 0.20, 0.00, & 
            0.08, 0.14, 0.13, & 
            0.00, 0.00, 0.00, &
            0.20, 0.20, 0.00 ] 

bioclimrt = 0.25                
                     
! ctem.f parameters: ----------

laimin = [ 1.0, 1.0, 0.0, &
           1.5, 1.0, 1.0, &  
           1.0, 1.0, 0.0, &
           0.01, 0.01, 0.0 ] 

laimax = [ 4.0, 3.0, 0.0, & 
           6.0, 5.0, 5.0, & 
           8.0, 8.0, 0.0, &
           4.0, 4.0, 0.0 ] 

! disturbance parameters: ------------

standreplace = [ 0.20, 0.20, 0.00, &
                 0.50, 0.20, 0.15, &  
                 0.00, 0.00, 0.00, &
                 0.25, 0.25, 0.00 ] 

! mainres.f parameters: ---------

bsrtstem = [ 0.0700, 0.0550, 0.0000, & ! IN BOTH PRESCRIBED AND COMPETE (TWO VERSIONS)
             0.0500, 0.0335, 0.0350, &  
             0.0365, 0.0365, 0.0000, & 
             0.0000, 0.0000, 0.0000 ] ! no stem component for grasses

bsrtroot = [ 0.5000, 0.2850, 0.0000, & ! IN BOTH PRESCRIBED AND COMPETE (TWO VERSIONS)
             0.4000, 0.2250, 0.1500, & 
             0.1600, 0.1600, 0.0000, & 
             0.1000, 0.1000, 0.0000 ]

! mortality.f parameters: ---------

maxage = [ 800.0, 500.0,   0.0, &  ! IN BOTH PRESCRIBED AND COMPETE (TWO VERSIONS)
           700.0, 450.0, 500.0, &  
             0.0,   0.0,   0.0, &
             0.0,   0.0,   0.0 ]

mxmortge = [ 0.005, 0.005, 0.00, & ! IN BOTH PRESCRIBED AND COMPETE (TWO VERSIONS)
             0.005, 0.005, 0.005, & 
             0.00, 0.00, 0.00, & 
             0.05, 0.10, 0.00 ] 

! phenology.f parameters: ---------

drlsrtmx = [ 0.006 , 0.005, 0.000, & ! IN BOTH PRESCRIBED AND COMPETE (TWO VERSIONS)
             0.010 , 0.025, 0.030, & 
             0.005 , 0.005, 0.000, &
             0.020 , 0.020, 0.000 ] 

!   ********************************************************************************************
!   =============                                                     ==========================
!   ============================== PRESCRIBED COVER PARAMETERS =================================
!   =============                                                     ==========================
!   ********************************************************************************************

else ! Prescribed PFT fractional cover

! These parameters are used when the PFT fractional cover is read in from the 
! CTM and INI files, or when LUC is on, the LUC file.

! allocate.f parameters: --------------

! ORIGINAL PRESCRIBED VALUES: FLAG JM Jan 11 2016.
omega = [ 0.80, 0.50, 0.00, &
          0.80, 0.80, 0.80, &
          0.05, 0.05, 0.00, &
          1.00, 1.00, 0.00 ]

epsilonl = [ 0.20, 0.06, 0.00, &
             0.35, 0.35, 0.25, &
             0.80, 0.80, 0.00, &
             0.01, 0.01, 0.00 ]

epsilons = [ 0.15, 0.05, 0.00, &
             0.05, 0.10, 0.10, &
             0.15, 0.15, 0.00, &
             0.00, 0.00, 0.00 ]

epsilonr = [ 0.65, 0.89, 0.00, &
             0.60, 0.55, 0.65, &
             0.05, 0.05, 0.00, &
             0.99, 0.99, 0.00 ]

! TESTING using the compete versions:
! omega = [ 0.80, 0.50, 0.00, & ! IN BOTH PRESCRIBED AND COMPETE (TWO VERSIONS)
!           0.80, 0.45, 0.80, &
!           0.05, 0.05, 0.00, &
!           1.00, 1.00, 0.00 ]
!
! epsilonl = [ 0.19, 0.45, 0.00, &  ! IN BOTH PRESCRIBED AND COMPETE (TWO VERSIONS)
!              0.39, 0.50, 0.30, &
!              0.80, 0.80, 0.00, &
!              0.10, 0.10, 0.00 ]
!
! epsilons = [ 0.40, 0.34, 0.00, & ! IN BOTH PRESCRIBED AND COMPETE (TWO VERSIONS)
!              0.21, 0.35, 0.10, &
!              0.15, 0.15, 0.00, &
!              0.00, 0.00, 0.00 ]
!
! epsilonr = [ 0.41, 0.21, 0.00, &  ! IN BOTH PRESCRIBED AND COMPETE (TWO VERSIONS)
!              0.40, 0.15, 0.60, &
!              0.05, 0.05, 0.00, &
!              0.90, 0.90, 0.00 ]

! --- to here.

! mainres.f parameters: ---------

bsrtstem = [ 0.0900, 0.0550, 0.0000, &
             0.0600, 0.0335, 0.0300, &
             0.0365, 0.0365, 0.0000, &
             0.0000, 0.0000, 0.0000 ] ! no stem component for grasses

bsrtroot = [ 0.5000, 0.2850, 0.0000, &
             0.6500, 0.2250, 0.0550, &
             0.1600, 0.1600, 0.0000, &
             0.1000, 0.1000, 0.0000 ]


! mortality.f parameters: ---------

maxage = [ 250.0, 400.0,   0.0, &    
           600.0, 250.0, 500.0, &  
             0.0,   0.0,   0.0, &
             0.0,   0.0,   0.0 ]

mxmortge = [ 0.005, 0.005, 0.00, &   ! Same as competition except for grasses.
             0.005, 0.005, 0.005, & 
             0.00, 0.00, 0.00, & 
             0.00, 0.00, 0.00 ] 


! phenology.f parameters: ---------

drlsrtmx = [ 0.0025, 0.005, 0.000, &
             0.005, 0.005, 0.025, &
             0.005, 0.005, 0.000, &
             0.050, 0.050, 0.000 ]    

end if


end subroutine initpftpars

end module ctem_params
