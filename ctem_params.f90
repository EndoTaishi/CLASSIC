!>\defgroup ctem_params_initpftpars

!> This module holds CTEM globally accessible parameters
!> These parameters are used in all CTEM subroutines
!> via use statements pointing to this module EXCEPT PHTSYN3.f
!> which has the information passed in via arguments. This is a legacy thing.

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
!>\file

module ctem_params

!>\ingroup ctem_params_main

!!@{

! J. Melton
! Jun 23 2013

! Change History:

! Jun 11 2015 - JM - Add in new disturb params (extnmois_duff and extnmois_veg) replacing duff_dry and extnmois.
! Mar 12 2014 - JM - Allow for two sets of paramters (competition and prescribed PFT fractions).
!                    all ctem subroutines except PHTSYN3 keep their parameters here.   
!
! Jan 17 2014 - JM - Add in more parameters from ctem.f, phenology.f, and allocate.f. 

! Remember that changes to this module will usually only take effect
! after you have done a 'make clean' then a 'make' (because it is a module).


implicit none

!public :: initpftpars
public :: readin_params

! Constants

real, parameter :: zero     = 1.0e-20
real, parameter :: abszero  = 1e-12    !<this one is for runclassctem.f and allocate.f

real, parameter :: pi       = 3.1415926535898d0
real, parameter :: earthrad = 6371.22 !<radius of earth, km
real, parameter :: km2tom2  = 1.0e+06  !<changes from \f$km^2\f$ to \f$m^2\f$
real, parameter :: deltat   = 1.0       !<CTEM's time step in days

! These month arrays are possibly overwritten in runclassctem due to leap years.
integer, dimension(12) :: monthdays = [ 31,28,31,30,31,30,31,31,30,31,30,31 ] !< days in each month
integer, dimension(13) :: monthend  = [ 0,31,59,90,120,151,181,212,243,273,304,334,365 ] !< calender day at end of each month
integer, dimension(12) :: mmday     = [ 16,46,75,106,136,167,197,228,259,289,320,350 ] !<mid-month day
integer, parameter :: nmon = 12 !< Number of months in a year

integer, parameter :: lon = 128 !< specify gcm resolution for longitude !FLAG - read in from netcdf!
integer, parameter :: lat = 64  !< specify gcm resolution for latitude !FLAG - read in from netcdf!

!> latitudes of the edges of the gcm grid cells for 128/x64 resolution !FLAG - read in from netcdf!
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
integer :: nlat = 1         !< Number of cells we are running, read in from the initialization file
integer :: nmos = 10        !< Number of mosaic tiles, read in from the initialization file
integer :: ilg = 10         !< nlat x nmos
integer :: ignd = 20        !< Number of soil layers, read in from the initialization file

! ----
! Plant-related
integer, parameter :: ican        = 4        !< Number of CLASS pfts, read in from the initialization file
integer, parameter :: icp1        = ican + 1 !
integer,parameter  :: icc=12                 !< Number of CTEM pfts (Peatlands add 3: EVG shrub,DCD shrubs, sedge)

integer,parameter  :: iccp1       = icc + 1  !

real, parameter :: seed    = 0.001    !< seed pft fraction, same as in competition \nin mosaic mode, all tiles are given this as a minimum
real, parameter :: minbare = 1.0e-5   !< minimum bare fraction when running competition on to prevent numerical problems.
real, parameter :: c2dom   = 450.0    !< gc / kg dry organic matter \nconversion factor from carbon to dry organic matter value is from Li et al. 2012 biogeosci
real, parameter :: wtCH4   = 16.044   !< Molar mass of CH4 ($g mol^{-1}$)


integer, parameter :: nbs         = 4        !


real :: tolrance = 0.0001d0 !< our tolerance for balancing c budget in kg c/m2 in one day (differs when competition on or not)
                            ! YW May 12, 2015 in peatland the C balance gap reaches 0.00016.

!> Logical switch for using constant allocation factors (default value is false)
logical :: consallo = .false.


! How to fit this stuff in?

integer, parameter :: l2max       = 5        !
integer, parameter :: kk          = 20       !< product of class pfts and l2max
integer, parameter :: numcrops    = 2        !< number of crop pfts
integer, parameter :: numtreepfts = 5        !< number of tree pfts
integer, parameter :: numgrass    = 3        !< number of grass pfts
integer, parameter :: numshrubs   = 2        !< number of shrubs pfts

!> simple crop matrix, define the number and position of the crops (NOTE: dimension icc)
logical, parameter, dimension(icc) :: crop = [ .false.,.false.,.false.,.false.,.false.,.false.,.false.,.true.,.true.,.false.,.false.,.false. ]

!> simple grass matric, define the number and position of grass
logical, parameter, dimension(icc) :: grass = [ .false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.true.,.true.,.true. ]

integer, parameter, dimension(numgrass) :: grass_ind = [ 10, 11,12 ]  !< index of the grass pfts (3 grass pfts at present)
integer, parameter, dimension(numshrubs) :: shrub_ind = [ 6, 7 ]  !not used for now
integer, parameter, dimension(numtreepfts) :: tree_ind = [ 1, 2,3,4,5 ]
integer, parameter, dimension(numcrops) :: crop_ind = [ 8,9 ]


! Read in from the namelist

integer, dimension(kk) :: modelpft      !<Separation of pfts into level 1 (for class) and level 2 (for ctem) pfts.
character(8), dimension(kk) :: pftlist  !<List of PFTs
real, dimension(kk) :: kn               !< Canopy light/nitrogen extinction coefficient; CAREFUL: Separate set defined in PHTSYN3.f!

!allocate.f parameters: ---------------------------------

real, dimension(kk) :: omega            !< omega, parameter used in allocation formulae (values differ if using prescribed vs. competition run)
real, dimension(kk) :: epsilonl         !< Epsilon leaf, parameter used in allocation formulae (values differ if using prescribed vs. competition run)
real, dimension(kk) :: epsilons         !< Epsilon stem, parameter used in allocation formulae (values differ if using prescribed vs. competition run)
real, dimension(kk) :: epsilonr         !< Epsilon root, parameter used in allocation formulae (values differ if using prescribed vs. competition run)
real, dimension(kk) :: rtsrmin          !< Minimum root:shoot ratio mostly for support and stability
real, dimension(kk) :: aldrlfon         !< Allocation to leaves during leaf onset
real, dimension(kk) :: caleaf           !< Constant allocation fractions to leaves if not using dynamic allocation.
                                        !!(NOT thoroughly tested, and using dynamic allocation is preferable)
real, dimension(kk) :: castem           !< Constant allocation fractions to stem if not using dynamic allocation.
                                        !!(NOT thoroughly tested, and using dynamic allocation is preferable)
real, dimension(kk) :: caroot           !< Constant allocation fractions to roots if not using dynamic allocation.
                                        !!(NOT thoroughly tested, and using dynamic allocation is preferable)

! bio2str.f parameters: ---------

real, dimension(kk) :: abar             !< parameter determining average root profile
real, dimension(kk) :: avertmas         !< average root biomass (kg c/m2) for ctem's 8 pfts used for estimating rooting profile
real, dimension(kk) :: alpha            !< parameter determining how the roots grow
real, dimension(kk) :: prcnslai         !< storage/imaginary lai is this percentage of maximum leaf area index that a given root+stem biomass can support
real, dimension(kk) :: minslai          !< minimum storage lai. this is what the model uses as lai when growing vegetation for scratch. consider these as model seeds.
real, dimension(kk) :: mxrtdpth         !< maximum rooting depth. this is used so that the rooting depths simulated by ctem's variable rooting depth parameterzation are
                                        !< constrained to realistic values
real, dimension(kk) :: albvis           !< visible albedos of the ctem pfts
real, dimension(kk) :: albnir           !< near IR albedos of the 9 ctem pfts

!! competition_mod.f90 parameters: ------

!! the model basically uses the temperature of the coldest month as
!! the major constraint for pft distribution. a range of the coldest
!! month temperature is prescribed for each pft within which pfts are
!! allowed to exist. in addition for tropical broadleaf drought
!! deciduous trees measure(s) of aridity (function of precipitation
!! and potential evaporation) are used.

! existence subroutine:

real, dimension(kk) :: tcoldmin         !< minimum coldest month temperature
real, dimension(kk) :: tcoldmax         !< maximum coldest month temperature
real, dimension(kk) :: twarmmax         !< maximum warmest month temperature
real, dimension(kk) :: gdd5lmt          !< minimum gdd above 5 c required to exist
real, dimension(kk) :: aridlmt          !< aridity index limit for broadleaf drought/dry deciduous trees
real, dimension(kk) :: dryseasonlmt     !< minimum length of dry season for PFT to exist
real, dimension(kk) :: bio2sap          !< multiplying factor for converting biomass density to sapling density
                                        !! smaller numbers give faster colonization rates.
real :: bioclimrt                       !< mortality rate (1/year) for pfts that no longer exist within their pre-defined bioclimatic range

! ctem.f parameters: ----------

real, dimension(kk) :: grescoef         !< Growth respiration coefficient
real, dimension(kk) :: humicfac         !< Humification factor - used for transferring carbon from litter into soil c pool
real, dimension(kk) :: laimin           !< Minimum lai below which a pft doesn't expand
real, dimension(kk) :: laimax           !< Maximum lai above which a pft always expands and lambdamax fraction of npp is used for expansion
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


!============================================================================
! !Separation of pfts into level 1 (for class) and level 2 (for ctem) pfts.
! integer, parameter, dimension(kk) :: modelpft= [ 1,     1,     0,     0,        0, &  ! CLASS PFT 1 NDL
!                              !                  EVG    DCD
!                                                  1,     1,     1,     1,        1, &  ! CLASS PFT 2 BDL
!                              !                  EVG  DCD-CLD DCD-DRY EVG-shrubs DCD-shrubs
!                                                  1,     1,     0,     0,        0,&  ! CLASS PFT 3 CROP
!                              !                  C3      C4
!                                                  1,     1,     1,     0,        0 ]   ! CLASS PFT 4 GRASS
!                              !                  C3      C4     sedge
!
! character(8), parameter, dimension(icc) :: pftlist = [ 'NdlEvgTr' , 'NdlDcdTr', 'BdlEvgTr','BdlDCoTr', &
!                                                      'BdlDDrTr','BdlEvgSh','BdlDcdSh','CropC3  ', &
!                                                      'CropC4  ','GrassC3 ','GrassC4 ','Sedge   ' ]
!Separation of pfts into level 1 (for class) and level 2 (for ctem) pfts.
!
! ! PFT specific parameters:
!
! ! Parameters used in more than one subroutine:
!
!
!> Canopy light/nitrogen extinction coefficient
!> (kn -> CAREFUL: Separate set defined in PHTSYN3.f!)
! real, dimension(kk) :: kn != [ 0.50, 0.50, 0.00, 0.00, 0.00, &
!                           !   0.50, 0.50, 0.50, 0.50, 0.50, &
!                           !   0.40, 0.48, 0.00, 0.00, 0.00, &
!                           !   0.46, 0.44, 0.46, 0.00, 0.00 ]

!allocate.f parameters: ---------------------------------

! real, dimension(kk) :: omega            !< omega, parameter used in allocation formulae (values differ if using prescribed vs. competition run)
! real, dimension(kk) :: epsilonl         !< Epsilon leaf, parameter used in allocation formulae (values differ if using prescribed vs. competition run)
! real, dimension(kk) :: epsilons         !< Epsilon stem, parameter used in allocation formulae (values differ if using prescribed vs. competition run)
! real, dimension(kk) :: epsilonr         !< Epsilon root, parameter used in allocation formulae (values differ if using prescribed vs. competition run)

! !> Logical switch for using constant allocation factors (default value is false)
! logical :: consallo = .false.

! !> Minimum root:shoot ratio mostly for support and stability
! real, dimension(kk) :: rtsrmin = [ 0.16, 0.16, 0.00, 0.00, 0.00, &
!                                    0.16, 0.16, 0.32, 0.16, 0.16, &
!                                    0.16, 0.16, 0.00, 0.00, 0.00, &
!                                    0.50, 0.50, 0.30, 0.00, 0.00 ]       !YW sedge has less roots fraction than the real grass
!
! !> Allocation to leaves during leaf onset
! real, dimension(kk) :: aldrlfon = [ 1.00, 1.00, 0.00, 0.00, 0.00, &
!                                     1.00, 1.00, 1.00, 1.00, 1.00, &
!                                     1.00, 1.00, 0.00, 0.00, 0.00, &
!                                     1.00, 1.00, 1.00, 0.00, 0.00 ]
!
! !> Constant allocation fractions to leaves if not using dynamic allocation. (NOT thoroughly tested, and using dynamic allocation is preferable)
! real, dimension(kk) :: caleaf = [ 0.275, 0.300, 0.000, 0.000, 0.000, &
!                                   0.200, 0.250, 0.250, 0.275, 0.300, &
!                                   0.400, 0.400, 0.000, 0.000, 0.000, &
!                                   0.450, 0.450, 0.450, 0.000, 0.000]
!
! !> Constant allocation fractions to stem if not using dynamic allocation. (NOT thoroughly tested, and using dynamic allocation is preferable)
! real, dimension(kk) :: castem = [ 0.475, 0.450, 0.000, 0.000, 0.000, &
!                                   0.370, 0.400, 0.400, 0.475, 0.450, &
!                                   0.150, 0.150, 0.000, 0.000, 0.000, &
!                                   0.000, 0.000, 0.000, 0.000, 0.000 ]
!
! !> Constant allocation fractions to roots if not using dynamic allocation. (NOT thoroughly tested, and using dynamic allocation is preferable)
! real, dimension(kk) :: caroot = [ 0.250, 0.250, 0.000, 0.000, 0.000, &
!                                   0.430, 0.350, 0.350, 0.250, 0.250, &
!                                   0.450, 0.450, 0.000, 0.000, 0.000, &
!                                   0.550, 0.550, 0.550, 0.000, 0.000 ]

! ! bio2str.f parameters: ---------
!
! !> parameter determining average root profile
! real, dimension(kk) :: abar    = [ 4.70, 5.86, 0.00, 0.00, 0.00, &
!                                    3.87, 3.46, 3.97, 8.50, 9.50, &  !shrubs in should have large abar because of Water restrictions
!                                    3.97, 3.97, 0.00, 0.00, 0.00, &  !so that rootdepth is shallower
!                                    5.86, 4.92, 9.50, 0.00, 0.00 ]
!
!
! !> average root biomass (kg c/m2) for ctem's 8 pfts used for estimating rooting profile
! real, dimension(kk) :: avertmas = [ 1.85, 1.45, 0.00, 0.00, 0.00, &
!                                     2.45, 2.10, 2.10, 1.50, 1.20, & !shrubs root biomass lower than trees YW
!                                     0.10, 0.10, 0.00, 0.00, 0.00, &
!                                     0.70, 0.70, 0.20, 0.00, 0.00 ] !Herbaceous above/below biomass ratio is 2.91 vs graminoids 0.63
!                                                   !in Tundra (Koerner & Renhardt 1987)
!
! !> parameter determining how the roots grow
! real, dimension(kk) :: alpha   = [ 0.80, 0.80, 0.00, 0.00, 0.00, &
!                                    0.80, 0.80, 0.80, 0.80, 0.80, &
!                                    0.80, 0.80, 0.00, 0.00, 0.00, &
!                                    0.80, 0.80, 0.80, 0.00, 0.00 ]
!
! !> storage/imaginary lai is this percentage of maximum leaf area index that a given root+stem biomass can support
! real, dimension(kk) :: prcnslai = [ 7.5, 7.5, 0.0, 0.0, 0.0, &
!                                     7.5, 7.5, 7.5, 7.5, 7.5, &
!                                     7.5, 7.5, 0.0, 0.0, 0.0 ,&
!                                     2.5, 2.5, 2.5, 0.0, 0.0 ]
!
! !> minimum storage lai. this is what the model uses as lai when growing vegetation for scratch. consider these as model seeds.
! real, dimension(kk) :: minslai = [ 0.3, 0.3, 0.0, 0.0, 0.0, &
!                                    0.3, 0.3, 0.3, 0.3, 0.3, &
!                                    0.2, 0.2, 0.0, 0.0, 0.0, &
!                                    0.2, 0.2, 0.2, 0.0, 0.0  ]
!
! !> maximum rooting depth. this is used so that the rooting depths simulated by ctem's variable rooting depth parameterzation are
! !> constrained to realistic values visible and near ir albedos of the 9 ctem pfts
! real, dimension(kk) :: mxrtdpth = [ 3.00, 3.00, 0.00, 0.00, 0.00, &
!                                     5.00, 5.00, 3.00, 1.00, 1.00, & !YW peatland shrubs lower values April 22, 2015
!                                     2.00, 2.00, 0.00, 0.00, 0.00, &
!                                     1.00, 1.00, 1.00, 0.00, 0.00 ]
!
! !> visible albedos of the 9 ctem pfts
! real, dimension(kk) :: albvis = [ 3.00, 3.00, 0.00, 0.00, 0.00, &
!                                   3.00, 5.00, 5.00, 3.00, 3.00, &
!                                   5.50, 5.50, 0.00, 0.00, 0.00, &
!                                   5.00, 6.00, 5.00, 0.00, 0.00 ]
!
! !> near IR albedos of the 9 ctem pfts
! real, dimension(kk) :: albnir = [ 19.0, 19.0, 0.00, 0.00, 0.00, &
!                                   23.0, 29.0, 29.0, 19.0, 19.0, &
!                                   34.0, 34.0, 0.00, 0.00, 0.00, &
!                                   30.0, 34.0, 30.0, 0.00, 0.00 ]

! !! competition_mod.f90 parameters: ------
!
! !! the model basically uses the temperature of the coldest month as
! !! the major constraint for pft distribution. a range of the coldest
! !! month temperature is prescribed for each pft within which pfts are
! !! allowed to exist. in addition for tropical broadleaf drought
! !! deciduous trees measure(s) of aridity (function of precipitation
! !! and potential evaporation) are used.
!
! ! existence subroutine:
!
! !> minimum coldest month temperature
! real, dimension(kk) :: tcoldmin = [-999.9, -999.9,   0.0,    0.0,    0.0, &
!                                       2.5,  -35.0,   4.0, -999.0, -999.0, &
!                                    -999.9, -999.9,   0.0,    0.0,    0.0, &
!                                    -999.9, -999.9,-999.9,    0.0,    0.0 ]
!
! !> maximum coldest month temperature
! real, dimension(kk) :: tcoldmax = [ 18.0,  -28.0,   0.0,    0.0,    0.0, &
!                                    999.9,   16.0, 900.0,   18.0,  -28.0, &
!                                    999.9,  999.9,   0.0,    0.0,    0.0, &
!                                    999.9,  999.9, 999.9,    0.0,    0.0  ]
!
! !> maximum warmest month temperature
! real, dimension(kk) :: twarmmax = [ 99.9,  25.0,  0.0,   0.0,   0.0,   &
!                                     99.9,  99.9, 99.9,  99.9,  25.0,   &
!                                     99.9,  99.9,  0.0,   0.0,   0.0,   &
!                                     99.9,  99.9, 99.9,   0.0,   0.0 ]
!
! !> minimum gdd above 5 c required to exist
! real, dimension(kk) :: gdd5lmt = [ 375.0,  600.0,  0.0,   0.0,   0.0, &
!                                   1200.0,  300.0,  9.9, 375.0, 600.0, &
!                                      9.9,    9.9,  0.0,   0.0,   0.0, &
!                                      9.9,    9.9,  9.9,   0.0,   0.0 ]
!
! !> aridity index limit for broadleaf drought/dry deciduous trees
! real, dimension(kk) :: aridlmt = [ 9.9,  9.9,  0.0,  0.0,  0.0, &
!                                    9.9,  9.9,  0.9,  9.9,  9.9, &
!                                    9.9,  9.9,  0.0,  0.0,  0.0, &
!                                    9.9,  9.9,  9.9,  0.0,  0.0 ]
!
! !> minimum length of dry season for PFT to exist
! real, dimension(kk) :: dryseasonlmt =[  9.0,  99.9,    0.0,   0.0,  0.0, &
!                                        99.9,  99.9,    5.5,   9.0, 99.9, &
!                                        99.9,  99.9,    0.0,   0.0,  0.0, &
!                                        99.9,  99.9,   99.9,   0.0,  0.0 ]
!
!
! !> multiplying factor for converting biomass density to sapling density
! ! smaller numbers give faster colonization rates.
! real, dimension(kk) :: bio2sap = [ 0.32, 0.20, 0.00, 0.00, 0.00, &
!                                    0.08, 0.14, 0.13, 0.32, 0.20, &
!                                    0.00, 0.00, 0.00, 0.00, 0.00, &
!                                    0.20, 0.20, 0.20, 0.00, 0.00 ]
!
! real :: bioclimrt = 0.25 !< mortality rate (1/year) for pfts that no longer exist within their pre-defined bioclimatic range

! ! ctem.f parameters: ----------
!
! !> Growth respiration coefficient
! real, dimension(kk) :: grescoef = [ 0.15, 0.15, 0.00, 0.00, 0.00, &
!                                     0.15, 0.15, 0.15, 0.15, 0.15, &
!                                     0.15, 0.15, 0.00, 0.00, 0.00, &
!                                     0.15, 0.15, 0.15, 0.00, 0.00 ]
!
!
! !> Humification factor - used for transferring carbon from litter into soil c pool
! real, dimension(kk) :: humicfac = [ 0.42, 0.42, 0.00, 0.00, 0.00, &
!                                     0.53, 0.48, 0.48, 0.42, 0.42, &
!                                     0.10, 0.10, 0.00, 0.00, 0.00, &
!                                     0.42, 0.42, 0.42, 0.00, 0.00 ]
!
! !> Minimum lai below which a pft doesn't expand
! real, dimension(kk) :: laimin = [ 1.0,   1.0,  0.0, 0.0, 0.0, &
!                                   1.5,   1.0,  1.0, 1.0, 1.0, &
!                                   1.0,   1.0,  0.0, 0.0, 0.0, &
!                                   0.01, 0.01, 0.01, 0.0, 0.0 ]
! !> Maximum lai above which a pft always expands and lambdamax fraction of npp is used for expansion
! real, dimension(kk) :: laimax = [ 4.0, 3.0, 0.0, 0.0, 0.0, &
!                                   6.0, 5.0, 5.0, 4.0, 3.0, &
!                                   8.0, 8.0, 0.0, 0.0, 0.0, &
!                                   4.0, 4.0, 4.0, 0.0, 0.0 ]
!
! real :: lambdamax      = 0.10 !< Max. fraction of npp that is allocated for reproduction/colonization
!
! real :: repro_fraction = 0.10 !< Fraction of NPP that is used to create reproductive tissues

! ! disturbance parameters: ------------
!
! real, dimension(2) :: bmasthrs_fire = [ 0.4, 1.2 ] !< min. and max. vegetation biomass thresholds to initiate fire, \f$kg c/m^2\f$
!
! real :: extnmois_veg  = 0.3  !< extinction moisture content for estimating vegetation fire likeliness due to soil moisture
!
! real :: extnmois_duff = 0.5  !< extinction moisture content for estimating duff layer fire likeliness due to soil moisture
!
! real :: lwrlthrs      = 0.25 ! FireMIP value: 0.025
!                              !< lower cloud-to-ground lightning threshold for fire likelihood flashes/km^2.year
!
! real :: hgrlthrs      = 10.0 ! FireMIP value: 1.0
!                              !< higher cloud-to-ground lightning threshold for fire likelihood flashes/km^2.year
!
! !>parameter m (mean) and b of logistic distribution used for \n
! !>**Parmlght was increased to 0.8 to make it so areas with higher amounts of
! !>lightning have higher lterm. The saturation is still the same, but the
! !>increase is more gradual at low lightning density. JM
! real :: parmlght    = 0.8
! real :: parblght    = 0.1    !< estimating fire likelihood due to lightning
!
! real :: reparea     = 500.0  !1000  ! FireMIP value: 300.0
!                              !< typical area representing ctem's fire parameterization (km2)
!
! real :: popdthrshld = 300.   !< threshold of population density (people/km2) [Kloster et al., biogeosci. 2010]
!
! !     **These values were being calculated each time, they shouldn't be as they
! !     are just parameters. JM
! !      real, parameter :: ymin = 0.01798621     !=1.0/( 1.0+exp((parmlght-0.0)/parblght) )
! !      real, parameter :: ymax = 0.9975273768   !=1.0/( 1.0+exp((parmlght-1.0)/parblght) )
! !      real, parameter :: slope = 0.0204588331  !=abs(0.0-ymin)+abs(1.0-ymax)
! !      real :: ymin, ymax, slope
!
! real :: alpha_fire       !< parameter alpha_fire and f0 used for estimating wind function for fire spread rate
! real :: f0 = 0.05        !< Fire spread rate in the absence of wind  -- Not used in CTEM v 2.0!

!> max. fire spread rate, km/hr
real, dimension(kk) :: maxsprd = [  0.38, 0.38, 0.00, 0.00, 0.00, &
                                    0.28, 0.28, 0.28, 0.38, 0.38, & !JM edit YW vals to be in line with new vals
                                    0.00, 0.00, 0.00, 0.00, 0.00, &
                                    0.51, 0.75, 0.51, 0.00, 0.00 ]

!> fraction of green leaf biomass converted to gases due to combustion
real, dimension(kk) :: frco2glf = [ 0.42, 0.42, 0.00, 0.00, 0.00,&
                                    0.42, 0.42, 0.42, 0.42, 0.42, &
                                    0.00, 0.00, 0.00, 0.00, 0.00, &
                                    0.48, 0.48, 0.48, 0.00, 0.00 ]

!> fraction of brown leaf biomass converted to gases due to combustion
real, dimension(kk) :: frco2blf = [ 0.00, 0.00, 0.00, 0.00, 0.00, &
                                    0.00, 0.00, 0.00, 0.00, 0.00, &
                                    0.00, 0.00, 0.00, 0.00, 0.00, &
                                    0.54, 0.54, 0.54, 0.00, 0.00 ]


!> fraction of green leaf biomass becoming litter after combustion
real, dimension(kk) :: frltrglf = [ 0.20, 0.20, 0.00, 0.00, 0.00, &
                                    0.20, 0.20, 0.20, 0.20, 0.20, &
                                    0.00, 0.00, 0.00, 0.00, 0.00, &
                                    0.10, 0.10, 0.10, 0.00, 0.00 ]

!> fraction of brown leaf biomass becoming litter after combustion
real, dimension(kk) :: frltrblf = [ 0.00, 0.00, 0.00, 0.00, 0.00, &
                                    0.00, 0.00, 0.00, 0.00, 0.00, &
                                    0.00, 0.00, 0.00, 0.00, 0.00, &
                                    0.06, 0.06, 0.06, 0.00, 0.00 ]

!> fraction of stem biomass converted to gases due to combustion
real, dimension(kk) :: frco2stm = [ 0.12, 0.12, 0.00, 0.00, 0.00, &
                                    0.12, 0.06, 0.06, 0.12, 0.12, &
                                    0.00, 0.00, 0.00, 0.00, 0.00, &
                                    0.00, 0.00, 0.00, 0.00, 0.00 ]

!> fraction of stem biomass becoming litter after combustion
real, dimension(kk) :: frltrstm = [ 0.60, 0.60, 0.00, 0.00, 0.00, &
                                    0.60, 0.40, 0.40, 0.60, 0.60, &
                                    0.00, 0.00, 0.00, 0.00, 0.00, &
                                    0.00, 0.00, 0.00, 0.00, 0.00  ]

!> fraction of root biomass converted to gases due to combustion
real, dimension(kk) :: frco2rt = [ 0.0, 0.0, 0.0, 0.0, 0.0, &
                                   0.0, 0.0, 0.0, 0.0, 0.0, &
                                   0.0, 0.0, 0.0, 0.0, 0.0, &
                                   0.0, 0.0, 0.0, 0.0, 0.0  ]

!> fraction of root biomass becoming litter after combustion
real, dimension(kk) :: frltrrt = [ 0.10, 0.10, 0.00, 0.00, 0.00, &
                                   0.10, 0.10, 0.10, 0.10, 0.10, &
                                   0.00, 0.00, 0.00, 0.00, 0.00, &
                                   0.25, 0.25, 0.25, 0.00, 0.00 ]

!> fraction of litter burned during fire and emitted as gases
real, dimension(kk) :: frltrbrn = [ 0.30, 0.30, 0.00, 0.00, 0.00, &
                                    0.36, 0.36, 0.36, 0.36, 0.36, &
                                    0.00, 0.00, 0.00, 0.00, 0.00, &
                                    0.42, 0.42, 0.42, 0.00, 0.00 ]


!> pft prevalence for stand replacing fire events (based on resistance to fire damage, ie. cambial kill)(unitless)
real, dimension(kk) :: standreplace = [ 0.20, 0.20, 0.00, 0.00, 0.00, &
                                        0.50, 0.20, 0.15, 0.20, 0.20, &
                                        0.00, 0.00, 0.00, 0.00, 0.00, &
                                        0.25, 0.25, 0.25, 0.00, 0.00 ]

!     emissions factors by chemical species
!
!     Values are from Andreae 2011 as described in Li et al. 2012
!     Biogeosci. Units: g species / (kg DOM)

!     Andreae 2011 as described in Li et al. 2012
!> pft-specific emission factors for CO2,g species / (kg DOM)
real, dimension(kk) :: emif_co2 = [ 1576.0, 1576.0,   0.00,   0.00,   0.00, &
                                    1604.0, 1576.0, 1654.0, 1576.0, 1576.0, &
                                    1576.0, 1654.0,   0.00,   0.00,   0.00, &
                                    1576.0, 1654.0, 1576.00,   0.00,   0.00 ]

!    values from Wiedinmyer et al. 2011
!real, dimension(kk) :: emif_co2 = [ 1514.0, 1514.0,   0.00, &
!             1643.0, 1630.0, 1716.0, &
!             1537.0, 1537.0,   0.00, &
!             1692.0, 1692.0,   0.00 ]

!     Andreae 2011 as described in Li et al. 2012
!> pft-specific emission factors for CO,g species / (kg DOM)
real, dimension(kk) :: emif_co = [ 106.0, 106.0, 0.00,  0.00,  0.00, &
                                    103.0, 106.0, 64.0, 106.0, 106.0, &
                                    106.0,  64.0, 0.00,  0.00,  0.00, &
                                    106.0,  64.0, 106.00,  0.00,  0.00 ]

!    values from Wiedinmyer et al. 2011
!real, dimension(kk) :: emif_co = [ 118.0, 118.0, 0.00, &
!             92.0, 102.0, 68.0, &
!            111.0, 111.0, 0.00, &
!             59.0,  59.0, 0.00 ]

!     Andreae 2011 as described in Li et al. 2012
!> pft-specific emission factors for CH4,g species / (kg DOM)
real, dimension(kk) :: emif_ch4 = [ 4.8, 4.8, 0.0, 0.0, 0.0, &
                                    5.8, 4.8, 2.4, 4.8, 4.8, &
                                    4.8, 2.4, 0.0, 0.0, 0.0, &
                                    4.8, 2.4, 4.8, 0.0, 0.0 ]

!    values from Wiedinmyer et al. 2011
!real, dimension(kk) :: emif_ch4 = [ 6.0, 6.0, 0.0, &
!             5.1, 5.0, 2.6, &
!             6.0, 6.0, 0.0, &
!             1.5, 1.5, 0.0 ]

!     Andreae 2011 as described in Li et al. 2012
!> pft-specific emission factors for non-methane hydrocarbons,g species / (kg DOM)
real, dimension(kk) :: emif_nmhc = [ 5.7, 5.7, 0.0, 0.0, 0.0, &
                                     6.4, 5.7, 3.7, 5.7, 5.7, &
                                     5.7, 3.7, 0.0, 0.0, 0.0, &
                                     5.7, 3.7, 5.7, 0.0, 0.0 ]

!    values from Wiedinmyer et al. 2011
!real, dimension(kk) :: emif_nmhc = [ 5.7, 5.7, 0.0, &
!              1.7, 5.7, 3.4, &
!              7.0, 7.0, 0.0, &
!              3.4, 3.4, 0.0 ]

!     Andreae 2011 as described in Li et al. 2012
!> pft-specific emission factors for H2,g species / (kg DOM)
real, dimension(kk) :: emif_h2 = [ 1.80, 1.80, 0.00, 0.00, 0.00, &
                                   2.54, 1.80, 0.98, 1.80, 1.80, &
                                   1.80, 0.98, 0.00, 0.00, 0.00, &
                                   1.80, 0.98, 1.80, 0.00, 0.00 ]


!     Andreae 2011 as described in Li et al. 2012
!> pft-specific emission factors for NOx,g species / (kg DOM)
real, dimension(kk) :: emif_nox = [ 3.24, 3.24, 0.00, 0.00, 0.00, &
                                    2.90, 3.24, 2.49, 3.24, 3.24, &
                                    3.24, 2.49, 0.00, 0.00, 0.00, &
                                    3.24, 2.49, 3.24, 0.00, 0.00 ]

!    values from Wiedinmyer et al. 2011 (species: "NOx (as NO)" from Table 1)
!real, dimension(kk) :: emif_nox = [ 1.80, 2.30, 0.00, &
!             2.60, 1.30, 3.90, &
!             3.50, 3.50, 0.00, &
!             2.80, 2.80, 0.00 ]

!     Andreae 2011 as described in Li et al. 2012
!> pft-specific emission factors for N2O,g species / (kg DOM)
real, dimension(kk) :: emif_n2o = [ 0.26, 0.26, 0.00, 0.00, 0.00, &
                                    0.23, 0.26, 0.20, 0.26, 0.26, &
                                    0.26, 0.20, 0.00, 0.00, 0.00, &
                                    0.26, 0.20, 0.26, 0.00, 0.00 ]

!     emission factors for aerosols

!     Andreae 2011 as described in Li et al. 2012
!> pft-specific emission factors for particles <2.5 micrometers in diameter,g species / (kg DOM)
real, dimension(kk) :: emif_pm25 = [ 12.7, 12.7,  0.0,  0.0,  0.0, &
                                     10.5, 12.7,  5.2, 12.7, 12.7, &
                                     12.7,  5.2,  0.0,  0.0,  0.0, &
                                     12.7,  5.2, 12.7,  0.0,  0.0 ]

!    values from Wiedinmyer et al. 2011
!real, dimension(kk) :: emif_pm25 = [ 13.0, 13.0, 0.0, &
!               9.7, 13.0, 9.3, &
!               5.8,  5.8, 0.0, &
!               5.4,  5.4, 0.0 ]

!     Andreae 2011 as described in Li et al. 2012
!> pft-specific emission factors for total particulate matter,g species / (kg DOM)
real, dimension(kk) :: emif_tpm = [ 17.6, 17.6,  0.0,  0.0,  0.0, &
                                    14.7, 17.6,  8.5, 17.6, 17.6, &
                                    17.6,  8.5,  0.0, 0.0,   0.0, &
                                    17.6,  8.5, 17.6, 0.0,   0.0 ]

!    values from Wiedinmyer et al. 2011
!real, dimension(kk) :: emif_tpm = [ 18.0, 18.0, 0.0, &
!             13.0, 18.0,15.4, &
!             13.0, 13.0, 0.0, &
!              8.3,  8.3, 0.0 ]


!     Andreae 2011 as described in Li et al. 2012
!> pft-specific emission factors for total carbon,g species / (kg DOM)
real, dimension(kk) :: emif_tc = [ 8.3, 8.3, 0.0, 0.0, 0.0, &
                                   7.2, 8.3, 3.4, 8.3, 8.3, &
                                   8.3, 3.4, 0.0, 0.0, 0.0, &
                                   8.3, 3.4, 8.3, 0.0, 0.0 ]

!    values from Wiedinmyer et al. 2011 (TPC in Table 1)
!real, dimension(kk) :: emif_tc = [ 8.3, 8.3, 0.0, &
!            5.2, 9.7, 7.1, &
!            4.0, 4.0, 0.0, &
!            3.0, 3.0, 0.0 ]

!     Andreae 2011 as described in Li et al. 2012
!> pft-specific emission factors for organic carbon,g species / (kg DOM)
real, dimension(kk) :: emif_oc = [ 9.1, 9.1, 0.0, 0.0, 0.0, &
                                   6.7, 9.1, 3.2, 9.1, 9.1, &
                                   9.1, 3.2, 0.0, 0.0, 0.0, &
                                   9.1, 3.2, 9.1, 0.0, 0.0 ]

!    values from Wiedinmyer et al. 2011
!real, dimension(kk) :: emif_oc = [ 7.8, 7.8, 0.0, &
!            4.7, 9.2, 6.6, &
!            3.3, 3.3, 0.0, &
!            2.6, 2.6, 0.0 ]



!     Andreae 2011 as described in Li et al. 2012
!> pft-specific emission factors for black carbon,g species / (kg DOM)
real, dimension(kk) :: emif_bc = [ 0.56, 0.56, 0.00, 0.00, 0.00, &
                                   0.56, 0.56, 0.47, 0.56, 0.56, &
                                   0.56, 0.47, 0.00, 0.00, 0.00, &
                                   0.56, 0.47, 0.56, 0.00, 0.00 ]

!    values from Wiedinmyer et al. 2011
!real, dimension(kk) :: emif_bc = [ 0.20, 0.20, 0.00, &
!            0.52, 0.56, 0.50, &
!            0.69, 0.69, 0.00, &
!            0.37, 0.37, 0.00 ]

! hetres parameters: ----------

!> litter respiration rates at 15 c in in kg c/kg c.year
real, dimension(kk) :: bsratelt = [ 0.4453, 0.5986, 0.0000, 0.0000, 0.0000, &
                                    0.6339, 0.7576, 0.6957, 0.4453, 0.5986, &
                                    0.6000, 0.6000, 0.0000, 0.0000, 0.0000, &
                                    0.5260, 0.5260, 0.5260, 0.0000, 0.0000 ]

!> soil carbon respiration rates at 15 c in kg c/kg c.year
real, dimension(kk) :: bsratesc = [ 0.0260, 0.0260, 0.0000, 0.0000, 0.0000, &
                                    0.0208, 0.0208, 0.0208, 0.0260, 0.0260, &
                                    0.0350, 0.0350, 0.0000, 0.0000, 0.0000, &
                                    0.0125, 0.0125, 0.0125, 0.0000, 0.0000 ]

!> Constants used in tanh formulation of respiration Q10 determination
real, dimension(4) :: tanhq10= [ 1.44, 0.56, 0.075, 46.0 ]
                               !   a     b      c     d
                               ! q10 = a + b * tanh[ c (d-temperature) ]
                               ! when a = 2, b = 0, we get the constant q10 of 2. if b is non
                               ! zero then q10 becomes temperature dependent

real :: alpha_hetres = 0.7     !< parameter for finding litter temperature as a weighted average of top soil layer temperature and root temperature

real :: bsratelt_g = 0.5605    !< bare ground litter respiration rate at 15 c in kg c/kg c.year

real :: bsratesc_g = 0.02258   !< bare ground soil c respiration rates at 15 c in kg c/kg c.year

real :: a = 4.0                !< parameter describing exponential soil carbon profile. used for estimating temperature of the carbon pool

! landuse_change_mod.f90 parameters: --------------
! NOTE: combust, paper, furniture ABSOLUTELY MUST add to 1.0!

real, dimension(3) :: combust   = [ 0.15, 0.30, 0.45 ] !< how much deforested/chopped off biomass is combusted (these absolutely must add to 1.0!)
real, dimension(3) :: paper     = [ 0.70, 0.70, 0.55 ] !< how much deforested/chopped off biomass goes into short term storage such as paper
real, dimension(3) :: furniture = [ 0.15, 0.00, 0.00 ] !< how much deforested/chopped off biomass goes into long term storage such as furniture
real, dimension(2) :: bmasthrs  = [ 4.0, 1.0 ]         !< biomass thresholds for determining if deforested area is a forest, a shrubland, or a bush kg c/m2

real :: tolrnce1 = 0.50  !FLAG would be good to make this consistent with global tolerance so only one value.
                         !< kg c, tolerance of total c balance

! mainres.f parameters: ----------

!    Base respiration rates for stem and root for ctem pfts in kg c/kg c.year at 15 degrees celcius. note that maintenance
!    respiration rates of root are higher because they contain both wood (coarse roots) and fine roots.
!    New parameter values introduced to produce carbon use efficiencies more in
!    line with literature (zhang et al. 2009, luyssaert et al. gcb 2007)
!    values changed for bsrtstem and bsrtroot. jm 06.2012

real, dimension(kk) :: bsrtstem !< Base respiration rates for stem in kg c/kg c.year at 15 degrees celcius (values differ if using prescribed vs. competition run)
real, dimension(kk) :: bsrtroot !< Base respiration rates for root in kg c/kg c.year at 15 degrees celcius (values differ if using prescribed vs. competition run)

real :: minlvfr = 0.05          !< Minimum live wood fraction



! peatlands_mod.f90 parameters:

! mosspht subroutine parameters

real :: rmlmoss25 = 1.1 !< Base dark respiration rate of moss umol CO2/m2/s (best estimate is from fig. 6b in
                        !! Williams and Flanagan (1998) seasonal variation not considered

real :: tau25m=2321.0   !< tau coefficient (rate at 25 Celcius) iscomputed on the basis of the specificity
                        !! factor (102.33) times Kco2/Kh2o (28.38) to convert for value in solution to that
                        !! based in air. The old value was 2321.1. New value (2904.12) Quercus robor (Balaguer
                        !! et al., 1996). Similar number from Dreyeret al. 2001, Tree Physiol, tau= 2710
                        !! 2600 from Brooks and Farquhar (1985) used AS A CONSTANT; 2321 Harley et al. (1992)

real :: gasc = 8.314    !< gas constant ($J mol^{-1} K^{-1}$)

real :: ektau =-29000.0 !< $J mol^{-1}$ (Jordan and Ogren, 1984)

real :: kc25= 30.2      !< kinetic coef for CO2 at 25 C (Pa)
                        !! Based upon:
                        !! 30.2 Arbutus unedo  (Harley et al., 1986)
                        !! 40.4 (von Caemmerer et al., 1994)
                        !! 27.46 deciduous forest understory (Wilson
                        !! et al. 2000); 30.0 (Collatz et al., 1991)
                        !! 30.0 Populus tremula,Coylus avelana,Tilia
                        !! cordata (Kull and Niinemets, 1998)

real :: ko25=24800.0    !< kinetic coef for O2 at 25C (rate)  (Pa)
                        !! Based upon:
                        !! 24800 Arbutus unedo (Harley et al.,1986)
                        !! 41980 deciduous forest understory (Wilson et al. 2000)
                        !! 30000 (Collatz et al., 1991)
                        !! 30000 (Kull and Niinemets, 1998) see kc25

real :: ec=63500.0      !< Activation energy for Kc of CO2 ($J mol^{-1}$) at 25C
                        !! Based upon:
                        !! 59430 Arbutus unedo (Harley et al., 1986)
                        !! 80500 deciduous forest understory (Wilson et al., 2000)
                        !! 63500 mosses and shrubs Tenhunen et al. (1994)

real :: ej=37000.0      !< Activation energy for electron transport, ($J mol^{-1}$)
                        !! Based upon:
                        !! 37000 at 25C J mol-1(Farquhar et al.(1980)
                        !! 55000 J mol-1 at 311K  (Harley & Baldocchi,
                        !! 1995,PCE) lichis tree

real :: eo=35000.0      !< Activation energy for Ko ($J mol^{-1}$) at 25C
                        !! Based upon:
                        !! 36000 Arbutus unedo (Harley et al.,1986)
                        !! 14500 for deciduous forest understory (Wilson et al., 2000)
                        !! 35000 mosses and shrubs (Tenhunen et al.,1994)

real :: evc=53000.0     !< Activation energy for carboxylation ($J mol^{-1}$) at 25C
                        !! Based upon:
                        !! 53000 (Kirschbaum & Farquhar, 1984)
                        !! 55000 (at 311K) (Harley & Baldocchi, 1995,PCE)

real :: sj=719.0        !< Constant affecting J at low temperature ($J mol^{-1}$)
                        !! at 25C; Based upon: 714.0 in Lloyd et al. (1995)
                        !! Macadamia integrifolia and Litchi chinesis

real :: hj=220300.0     !< Constant affecting J at high temperature ($J mol^{-1}$)
                        !! at 25C; Based upon: 220300 Macadamia integrifolia
                        !! and Litchi chinesis Lloyd et al. (1995)
                        !! 206083 Pinus radiata (Wilson et al. 2000)

real :: alpha_moss=0.21 !< Efficiency of conversion of incident photons
                        !! into electrons(mol electron/mol photon),
                        !! not sensitive, alpha = 0.21 Flanagan's
                        !! value for sphagnum and pleurozium
                        !! alpha = 0.3 (Seel et al.,1992) for
                        !! T. ruraliformis and d. palustris

! decp subroutine parameters:

real::     dctmin= 269.0        !< minimum temperature of soil respiration,K (peatland soils)

real::     dcbaset=283.0        !< base temperature for Q10, K (peatland soils)

real::     bsrateltms = 0.20    !< heterotrophic respiration base rate for peatlands (yr-1)
                                !! Based upon 0.05 (T.Moore unpubished data)
                                !! for shrubs 0.2 from Moore and around
                                !! 0,5 for ctem. TM's values is natural
                                !! condition, CTEM is opt at 15degrees,
                                !! So TM's values can be scaled up
                                !! more rates in Aerts 1997, Moore 2007

! peatland parameters used in other subroutines:

real :: zolnmoss = -6.57    !< natual logarithm of roughness length of moss surface.
                            !! Based upon:
                            !! roughness length as 0.0014 for momentum (Oke 1997, Raddatz et al, 2009)

real :: thpmoss = 0.98      !< pore volume of moss ($m^3/m^3$)
                            !! Based upon:
                            !! 0.89 Price and Whittington 2010
                            !! 98.2 (Sphagnum) and 98.9 (feather) in
                            !! O'Donnell et al. 2009; 0.90 in Berlinger
                            !! et al.2001

real :: thrmoss = 0.20      !< Liquid water retention capacity of moss ($m^3/m^3$)
                            !! Based upon: Price and Whittington 2010 for upper 5cm
                            !! moss 0.275/0.6 comparing to 0.39/0.62 in
                            !! fibric for the feather moss.For spahgnum
                            !! moss(O'Donnell et al.2009)

real :: thmmoss = 0.01      !< residual liquid water content after freezing or
                            !! evaporation of moss ($m^3/m^3$)
                            !! Based upon: Price and Whittington 2010

real :: bmoss = 2.3         !< Clapp and Hornberger empirical "b" parameter of moss
                            !! Based upon: Berlinger et al.2001(they set 4.0 for peat)

real :: psismoss =0.0103    !< Soil moisure suction at saturation of moss (m)
                            !! Based upon: moss and peat using same values in Berlinger
                            !! et al. 2001)

real :: grksmoss = 1.83E-3  !< Saturated hydrualic conductivity of moss (m)
                            !! The saturated hydrualic conductivity is
                            !! much higher in living moss than lower
                            !! layer (Price et al. 2008). 1.83E-3 for
                            !! 5cm moss, their value of 10cm moss is
                            !! 2.45E-4, similar to fibric here

real :: hcpmoss = 2.50E6    !< Volumetric heat capacity of moss ($J m^{-3} K^{-1}$), same as hcp of peat in Beringer et al. 2001  2016

!real :: sphms = 2.70E3      !< Specific heat of moss layer ($J m^{-2} K^{-1}$) ! FLAG NOT USED also check units - written as J/kg/K. JM Sep 26 2016
                            !! same as that specific heat for vegetation
                            !! (Berlinger et al 2001)

!real :: rhoms = 40.0        !< Density of dry moss ($kg m^{-3}$)  ! FLAG NOT USED JM Oct 2016.
                            !! Based upon:
                            !! 40.0 Price et al. 2008, Price and
                            !! Whittington 2010 value for feather moss
                            !! is lower than sphagnmum (20 kg/m3 in
                            !! O'Donnell et al. 2009)

!real :: slams = 20.0        !< Specific leaf area of moss ($m^2 kg^{-1}$), ! FLAG NOT USED JM Oct 2016.
                            !! Based upon:
                            !! S vensson 1995 sphangum value ranges from
                            !! 13.5~47.3 m2/kg (Lamberty et al. 2007)

real::   grescoefmoss = 0.15!< Moss growth respiration coefficient

real::   rmortmoss = 0.7    !< Moss mortality rate constant (year-1),assumed based on
                            !! literature (eg. 54g litter/75g biomass in Moore et al,2002)

real::   humicfacmoss=2.0   !< Humification ratio of moss litter, higher in mosses than
                            !!  vascular pfts because the decomposibility of moss litter is
                            !!  lower.It's calibrated over the training sites.

! phenology.f parameters: ----------

!> Number of days over which to check if net photosynthetic rate is positive before initiating leaf onset
integer, dimension(kk) :: dayschk  = [ 7, 7, 0, 0, 0, &
                                       7, 7, 7, 7, 7, &
                                       7, 7, 0, 0, 0, &
                                       7, 7, 7, 0, 0 ]

!> Parameter determining how fast soil dryness causes leaves to fall
real, dimension(kk) :: drgta = [ 3.0, 3.0, 0.0, 0.0, 0.0, &
                                 3.0, 3.0, 3.0, 0.0, 0.0, &
                                 3.0, 3.0, 0.0, 3.0, 3.0, &
                                 3.0, 3.0, 3.0, 0.0, 0.0 ]


!> eta and kappa, parameters for estimating min. stem+root biomass
real, dimension(kk) :: eta = [ 10.0, 30.8, 0.00, 0.00, 0.00, &
                               31.0, 50.0, 30.0, 10.0, 30.8, &
                                7.0,  7.0,   0.0,  0.0, 0.00, &
                                3.0,  3.0,   3.0,  0.0, 0.00 ]

!> required to support green leaf biomass. kappa is 1.6 for trees and crops, and 1.2 for grasses.
real, dimension(kk) :: kappa =[ 1.6, 1.6, 0.0, 0.0, 0.0, &
                                1.6, 1.6, 1.6, 1.6, 1.6, &
                                1.6, 1.6, 0.0, 0.0, 0.0, &
                                1.2, 1.2, 1.2, 0.0, 0.0 ]

real, dimension(2) :: flhrspan = [ 17.0, 45.0 ]  !< Harvest span (time in days over which crops are harvested,  15 days),
                                                 !< and  fall span (time in days over which bdl cold dcd plants shed their leaves,  30 days)

real :: fracbofg = 0.55 !< Parameter used to estimate lai of brown leaves. We assume that SLA of brown leaves is this fraction of SLA of green leaves

!> LAI threshold for harvesting crops. values are zero for all pftsother than c3 and c4 crops.
real, dimension(kk) :: harvthrs = [ 0.0, 0.0, 0.0, 0.0, 0.0, &
                                    0.0, 0.0, 0.0, 0.0, 0.0, &
                                    4.5, 3.5, 0.0, 0.0, 0.0, &
                                    0.0, 0.0, 0.0, 0.0, 0.0 ]

!> CTEM can use user-specified specific leaf areas (SLA) if the following specified values are greater than zero
real, dimension(kk) :: specsla =[  0.0, 0.0, 0.0, 0.0, 0.0, &  ! Not used, will be used if these are non-zero.
                                   0.0, 0.0, 0.0, 0.0, 0.0, &
                                   0.0, 0.0, 0.0, 0.0, 0.0, &
                                   0.0, 0.0, 0.0, 0.0, 0.0 ]

!> Percentage of max. LAI that can be supported which is used as a threshold for determining leaf phenology status
real, dimension(kk) :: thrprcnt = [ 40.0, 40.0,  0.0,  0.0,  0.0, &
                                    40.0, 50.0, 50.0, 40.0, 40.0, &
                                    50.0, 50.0,  0.0,  0.0,  0.0, &
                                    40.0, 40.0, 40.0,  0.0,  0.0 ]

!> Lower temperature threshold for ctem's 9 pfts. these are used to estimate cold stress related leaf loss rate (degree c)
real, dimension(kk) :: lwrthrsh = [ -50.0, -5.0, 0.0,   0.0,  0.0, &
                                      5.0,  8.0, 5.0, -50.0, -5.0, &
                                      5.0,  5.0, 0.0,   0.0,  0.0, &
                                      0.1,  5.0, 0.1,   0.0,  0.0 ]

!> Max. loss rate for cold stress for all 9 pfts, (1/day)
real, dimension(kk) :: cdlsrtmx = [ 0.10, 0.30, 0.00, 0.00, 0.00, &
                                    0.30, 0.40, 0.15, 0.10, 0.30, &
                                    0.15, 0.15, 0.00, 0.00, 0.00, &
                                    0.15, 0.15, 0.15, 0.00, 0.00 ]

real :: kmort1 = 0.3            !< kmort1, parameter used in growth efficiency mortality formulation

real, dimension(kk) :: mxmortge !< Maximum mortality when growth efficiency is zero (1/year) (values differ if using prescribed vs. competition run)
real, dimension(kk) :: maxage   !< Maximum plant age. used to calculate intrinsic mortality rate.
                                !< maximum age for crops is set to zero since they will be harvested
                                !< anyway. grasses are treated the same way since the turnover time
                                !< for grass leaves is 1 year and for roots is 2 years. (values differ if using prescribed vs. competition run)



!> Leaf life span (in years) for CTEM's 9 pfts
real, dimension(kk) :: lfespany = [ 5.00, 1.00, 0.00, 0.00, 0.00, &          !YW July 13, 2015  add 4 PFT slots and 2 with values for shrubs
                                    1.50, 1.00, 1.00, 5.00, 0.40, &  !PFT 3 was 1.75 (from IBIS), 2.00 follows LPJ. JM Mar 2014.
                                    1.75, 1.75, 0.00, 0.00, 0.00, &
                                    1.00, 1.00, 1.00, 0.00, 0.00 ]


real, dimension(kk) :: drlsrtmx !< Max. loss rate for drought stress for all 9 pfts, (1/day) (values differ if using prescribed vs. competition run)

!> Parameter determining how fast cold temperatures causes leaves to fall
real, dimension(kk) :: colda = [ 3.0, 3.0, 0.0, 0.0, 0.0, &
                                 3.0, 3.0, 3.0, 3.0, 3.0, &
                                 3.0, 3.0, 0.0, 0.0, 0.0, &
                                 3.0, 3.0, 3.0, 0.0, 0.0 ]


!> No. of days for which some temperature has to remain below a given threshold for initiating a process
integer, dimension(2) :: coldlmt = [ 7 , 5 ]   ! days

!>1. -5 c threshold for initiating "leaf fall" mode for ndl dcd trees \n
!>2.  8 c threshold for initiating "harvest" mode for crops the array colddays(i,2) tracks days corresponding to these thresholds
real, dimension(2) :: coldthrs = [ -5.0 , 8.0 ]

real :: roothrsh = 8.0 !< Root temperature threshold for initiating leaf onset for cold broadleaf deciduous pft, degrees celcius

! turnover.f parameters: ---------------------------------

!> Stemlife, turnover time scale for stem for different pfts
real, dimension(kk) :: stemlife = [ 86.3, 86.3, 0.00, 0.00, 0.00, &
                                    80.5, 80.5, 75.8, 65.0, 75.0, & !shrub values < trees in ctem2.0 was 45 for pft4 and 5 was used in v6.
                                    20.0, 20.0, 0.00, 0.00, 0.00, &
                                    0.00, 0.00, 0.00, 0.00, 0.00 ]

!> Rootlife, turnover time scale for root for different pfts
real, dimension(kk) :: rootlife = [ 13.8,13.2, 0.0,  0.0,  0.0, &
                                    12.7,10.9, 9.8, 11.5, 12.0, & !shrubs values<trees, in ctem 2.0 was 5.5 for pft 4 and 5 and was used and shrubs in v6.
                                     3.0, 3.0, 0.0,  0.0,  0.0, &
                                     3.0, 3.0, 3.0,  0.0,  0.0  ]

real :: stmhrspn = 17.0 !< Stem harvest span. same as crop harvest span. period in days over which crops are harvested.

! wetland_methane.f90 parameters: ------------

!>methane to carbon dioxide flux scaling factor.\n
!>Rita Wania's thesis suggests about 0.25, but we get a better agreement to outputs from the Walter's model if we use 0.16.
!>Note that this scaling factor likely is temperature dependent, and increases with temperature, but it is difficult to
!>know the function, so leave constant for now ratio is \f$mol ch_4\f$ to \f$mol co_2\f$
!real :: ratioch4     = 0.16  ! old value.
real :: ratioch4     = 0.135 ! New value based on Vivek's work of Aug 2016.

!>ratio of wetland to upland respiration\n
!>Use the heterotrophic respiration outputs for soil and litter as the ecosystem basis.  These were summed as "hetrores".
!>This respiration is for upland soils; we multiply by wtdryres as the ratio of wetland to upland respiration
!>based on literature measurements: Dalva et al. 1997 found 0.5 factor; Segers 1998 found a 0.4 factor. use 0.45 here (unitless)
real :: wtdryres     = 0.45
real :: factor2      = 0.015  !< constant value for secondary (ch4wet2) methane emissions calculation
real :: lat_thrshld1 = 40.0   !< Northern zone for wetland determination (degrees North)
real :: lat_thrshld2 = -35.0  !< Boundary with southern zone for wetland determination (degrees North)
real :: soilw_thrshN = 0.55   !< Soil wetness threshold in the North zone
real :: soilw_thrshE = 0.80   !< Soil wetness threshold in the Equatorial zone
real :: soilw_thrshS = 0.70   !< Soil wetness threshold in the South zone

! ----=====-------=========-----------========---------=========--------========-----------==========---------=======---========**

contains

subroutine initpftpars(compete)

!>\ingroup ctem_params_initpftpars
!!@{

implicit none

!>true if the competition scheme is on.
logical, intent(in) :: compete


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

omega = [ 0.80, 0.50, 0.00, 0.00, 0.00, &
          0.80, 0.45, 0.80, 0.80, 0.50, &
          0.05, 0.05, 0.00, 0.00, 0.00, &
          1.00, 1.00, 1.00, 0.00, 0.00 ]

epsilonl = [ 0.19, 0.45, 0.00, 0.00, 0.00, &
             0.39, 0.50, 0.30, 0.19, 0.45, &
             0.80, 0.80, 0.00, 0.00, 0.00, &
             0.10, 0.10, 0.60, 0.00, 0.00]        !Sedge YW less allocation to root than real grass


epsilons = [ 0.40, 0.34, 0.00, 0.00, 0.00, &
             0.21, 0.35, 0.10, 0.40, 0.34, &
             0.15, 0.15, 0.00, 0.00, 0.00, &
             0.00, 0.00, 0.00, 0.00, 0.00 ]

epsilonr = [ 0.41, 0.21, 0.00, 0.00, 0.00, &
             0.40, 0.15, 0.60, 0.41, 0.21, &
             0.05, 0.05, 0.00, 0.00, 0.00, &
             0.90, 0.90, 0.40, 0.00, 0.00 ]       !Sedge YW less allocation to root than real grass

! mainres.f parameters: ---------

bsrtstem = [ 0.0700, 0.0550, 0.0000, 0.0000, 0.0000, &
             0.0500, 0.0335, 0.0350, 0.0700, 0.0335, &
             0.0365, 0.0365, 0.0000, 0.0000, 0.0000, &
             0.0000, 0.0000, 0.0000, 0.0000, 0.0000 ] ! no stem component for grasses

bsrtroot = [ 0.5000, 0.2850, 0.0000, 0.0000, 0.0000, &
             0.4000, 0.2250, 0.1500, 0.5000, 0.2850, &
             0.1600, 0.1600, 0.0000, 0.0000, 0.0000, &
             0.1000, 0.1000, 0.1000, 0.0000, 0.0000 ]

! mortality.f parameters: ---------

maxage = [ 800., 500.,   0.,   0.,   0., &
           700., 450., 500., 800., 500., &
             0.,   0.,   0.,   0.,   0., &
             0.,   0.,   0.,   0.,   0. ]


mxmortge = [ 0.005, 0.005, 0.000, 0.000, 0.000, &
             0.005, 0.005, 0.005, 0.005, 0.005, &
             0.000, 0.000, 0.000, 0.000, 0.000, &
             0.050, 0.100, 0.050, 0.000, 0.000 ]

! phenology.f parameters: ---------


drlsrtmx = [ 0.006 , 0.005, 0.000, 0.000, 0.000, &
             0.010 , 0.025, 0.030, 0.006, 0.005, &
             0.005 , 0.005, 0.000, 0.000, 0.000, &
             0.020 , 0.020, 0.020, 0.000, 0.000 ]


else ! Prescribed PFT fractional cover

!   ********************************************************************************************
!   =============                                                     ==========================
!   ============================== PRESCRIBED COVER PARAMETERS =================================
!   =============                                                     ==========================
!   ********************************************************************************************

! These parameters are used when the PFT fractional cover is read in from the
! CTM and INI files, or when LUC is on, the LUC file.


! allocate.f parameters: --------------

omega = [ 0.80, 0.50, 0.00, 0.00, 0.00, &
          0.80, 0.80, 0.80, 0.80, 0.50, &
          0.05, 0.05, 0.00, 0.00, 0.00, &
          1.00, 1.00, 1.00, 0.00, 0.00 ]

epsilonl = [ 0.20, 0.06, 0.00, 0.00, 0.00, &
             0.35, 0.35, 0.25, 0.20, 0.06, &
             0.80, 0.80, 0.00, 0.80, 0.80, &
             0.01, 0.01, 0.01, 0.00, 0.00 ]

epsilons = [ 0.15, 0.05, 0.00, 0.00, 0.00, &
             0.05, 0.10, 0.10, 0.15, 0.05, &
             0.15, 0.15, 0.00, 0.00, 0.00, &
             0.00, 0.00, 0.00, 0.00, 0.00 ]

epsilonr = [ 0.65, 0.89, 0.00, 0.00, 0.00, &
             0.60, 0.55, 0.65, 0.65, 0.89, &
             0.05, 0.05, 0.00, 0.00, 0.00, &
             0.99, 0.99, 0.99, 0.00, 0.00 ]

! mainres.f parameters: ---------

bsrtstem = [ 0.0900, 0.0550, 0.0000, 0.0000, 0.0000, &
             0.0600, 0.0335, 0.0300, 0.0900, 0.0550, &
             0.0365, 0.0365, 0.0000, 0.0000, 0.0000, &
             0.0000, 0.0000, 0.0000, 0.0000, 0.0000 ] ! no stem component for grasses

bsrtroot = [ 0.5000, 0.2850, 0.0000, 0.0000, 0.0000, &
             0.6500, 0.2250, 0.0550, 0.5000, 0.2850, &
             0.1600, 0.1600, 0.0000, 0.0000, 0.0000, &
             0.1000, 0.1000, 0.1000, 0.0000, 0.0000 ]


! mortality.f parameters: ---------

maxage = [ 250.0, 400.0,   0.0,   0.0,   0.0,  &    !same as comp
           600.0, 250.0, 500.0, 250.0, 400.0,  &
             0.0,   0.0,   0.0,   0.0,   0.0,  &
             0.0,   0.0,   0.0,   0.0,   0.0 ]

mxmortge = [ 0.005, 0.005,  0.00,  0.00,  0.00, &   ! Same as competition except for grasses.
             0.005, 0.005, 0.005, 0.005, 0.005, &
              0.00,  0.00,  0.00,  0.00,  0.00, &
              0.00,  0.00,  0.00,  0.00,  0.00 ]


! phenology.f parameters: ---------



drlsrtmx = [ 0.0025, 0.005, 0.000, 0.000, 0.000, &
             0.005, 0.005, 0.025, 0.0025, 0.005, &
             0.005, 0.005, 0.000,  0.000, 0.000, &
             0.050, 0.050, 0.050,  0.000, 0.000 ]

end if

end subroutine initpftpars
 !----------------------------------------------------\

subroutine readin_params

!>\ingroup readin_params
!!@{

implicit none

namelist /classicparams/ &
    modelpft,&
    pftlist,&
    kn, &
    omega,&
    epsilonl,&
    epsilons,&
    epsilonr,&
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
    f0


open(10,file='template_run_parameters.txt',action='read',status='old')

read(10,nml = classicparams)

close(10)

print *,popdthrshld

! FLAG ASSIGN the compete values vs. prescribed here!!!

end subroutine readin_params

!>@}
end module ctem_params
