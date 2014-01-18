module ctem_params

! This module holds CTEM globally accessible parameters
! These parameters are used in all CTEM subroutines
! via use statements pointing to this module EXCEPT PHTSYN3.f 
! which has the information passed in via arguments. This is a legacy thing.

! Remember that changes to this module will usually only take effect
! after you have done a make clean then a make.

! This module is structured with general parameters first then the PFT specific
! parameters later.

! J. Melton
! Jun 23 2013

! Change History:

! Jan 17 2014 - JM - Add in more parameters from ctem.f, phenology.f, and allocate.f. We
!                    are going to carry different parameters for the model run in competition and
!                    prescribed fraction mode (for now at least).

implicit none

public :: initpftpars

! Constants

real, parameter :: zero = 1.0e-20
real, parameter :: abszero=1e-12 ! this one is for runclassctem.f and allocate.f 

real, parameter :: pi =3.1415926535898d0
real, parameter :: earthrad = 6371.22   ! radius of earth, km
real, parameter :: km2tom2 = 1.0e+06 ! changes from km2 to m2
real, parameter :: deltat =1.0  !CTEM's time step in days


real, parameter :: tolrance = 0.00001d0 ! our tolerance for balancing c budget in kg c/m2 in one day, this is 1/100th of a gram of c

integer, parameter, dimension(12) :: monthdays = [ 31,28,31,30,31,30,31,31,30,31,30,31 ] ! days in each month
integer, parameter, dimension(13) :: monthend = [ 0,31,59,90,120,151,181,212,243,273,304,334,365 ] ! calender day at end of each month
integer, parameter, dimension(12) :: mmday= [ 16,46,75,106,136,167,197,228,259,289,320,350 ] !mid-month day

integer, parameter :: lon = 96  ! specify gcm resolution for longitude
integer, parameter :: lat = 48  ! specify gcm resolution for latitude

! latitudes of the edges of the gcm grid cells for 96x48 resolution
real, parameter, dimension(lat+1) :: edgelat = [ -90.0000, -85.3190, -81.6280, -77.9236, -74.2159,  &
                                        -70.5068, -66.7970, -63.0868, -59.3763, -55.6657, -51.9549, &
                                        -48.2441, -44.5331, -40.8221, -37.1111, -33.4001, -29.6890, &
                                        -25.9779, -22.2668, -18.5557, -14.8446, -11.1335,  -7.4223, &
                                         -3.7112,   0.0000,   3.7112,   7.4223,  11.1335,  14.8446, &
                                         18.5557,  22.2668,  25.9779,  29.6890,  33.4001,  37.1111, &
                                         40.8221,  44.5331,  48.2441,  51.9549,  55.6657,  59.3763, &
                                         63.0868,  66.7970,  70.5068,  74.2159,  77.9236,  81.6280, &
                                         85.3190,  90.0000 ]

! ----
! Model state
integer, parameter :: nlat=1            !
integer, parameter :: nmos=10           ! Number of mosaic tiles
integer, parameter :: ilg = nlat*nmos   !
integer, parameter :: nmon=12           ! Number of months in a year

! ----
! Plant-related

integer, parameter :: ican=4            ! Number of CLASS pfts
integer, parameter :: ignd=3            ! Number of soil layers
integer, parameter :: icp1 = ican + 1   !
integer,parameter  :: icc=9             ! Number of CTEM pfts
integer,parameter  :: iccp1 = icc + 1   !
integer, parameter :: l2max = 3         !
integer, parameter :: kk = 12           ! product of class pfts and l2max
integer, parameter :: numcrops = 2      ! number of crop pfts
integer, parameter :: numtreepfts = 5   ! number of tree pfts
integer, parameter :: numgrass = 2      ! number of grass pfts

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

! conversion factor from carbon to dry organic matter value is from Li et al. 2012 biogeosci
real, parameter :: c2dom = 450.0 ! gc / kg dry organic matter

! ================================================

! PFT specific parameters:

! Parameters used in more than one subroutine:

real, dimension(kk) :: kn               ! Canopy light/nitrogen extinction coefficient

! allocate.f parameters: ---------------------------------

! Parameterization values based on comparison mostly with LUYSSAERT, S.et al. CO2 balance
! of boreal, temperate, and tropical forests derived from a global database,
! Glob. Chang. Biol., 13(12), 2509–2537, 2007. and informed by LITTON, et al. Carbon
! allocation in forest ecosystems, Glob. Chang. Biol., 13(10), 2089–2109, 2007. Further
! tuning was performed on these basic values. JM Dec 20 2013.

real, dimension(kk) :: omega            ! omega, parameter used in allocation formulae
real, dimension(kk) :: epsilonl         ! Epsilon leaf, parameter used in allocation formulae
real, dimension(kk) :: epsilons         ! Epsilon stem, parameter used in allocation formulae
real, dimension(kk) :: epsilonr         ! Epsilon root, parameter used in allocation formulae
logical :: consallo                     ! Logical switch for using constant allocation factors (default value is false)
real, dimension(kk) :: rtsrmin          ! Minimum root:shoot ratio mostly for support and stability
real, dimension(kk) :: aldrlfon         ! Allocation to leaves during leaf onset

! Constant allocation fractions if not using dynamic allocation. (NOT thoroughly tested, and using dynamic allocation is preferable)
real, dimension(kk) :: caleaf 
real, dimension(kk) :: castem 
real, dimension(kk) :: caroot

! ctem.f parameters: ----------

real, dimension(kk) :: grescoef         ! Growth respiration coefficient 
real, dimension(kk) :: humicfac         ! Humification factor - used for transferring carbon from litter into soil c pool
real, dimension(kk) :: laimin           ! Minimum lai below which a pft doesn't expand
real, dimension(kk) :: laimax           ! Maximum lai above which a pft always expands and lambdamax fraction of npp is used for expansion
real :: lambdamax                       ! Max. fraction of npp that is allocated for reproduction/colonization

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
real, dimension(2) :: coldlmt           ! No. of days for which some temperature has to remain below a given threshold for initiating a process
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

contains

! ----=====-------=========-----------========---------=========--------========-----------==========---------=======---========**

subroutine initpftpars()

! PFT parameters

!     Note the structure of vectors which clearly shows the CLASS
!     PFTs (along rows) and CTEM sub-PFTs (along columns)
!
!     needle leaf |  evg       dcd       ---
!     broad leaf  |  evg   dcd-cld   dcd-dry
!     crops       |   c3        c4       ---
!     grasses     |   c3        c4       ---


implicit none

! Parameters used in more than one subroutine:

kn= [ 0.50, 0.50, 0.00, &
      0.50, 0.50, 0.50, &
      0.40, 0.48, 0.00, &
      0.46, 0.44, 0.00 ]


! allocate.f parameters: --------------

omega = [ 0.80, 0.50, 0.00, & 
          0.80, 0.50, 0.80, &
          0.05, 0.05, 0.00, &
          1.00, 1.00, 0.00 ]

consallo = .false.

epsilonl = [ 0.19, 0.45, 0.00, &  !pft 2 was 0.06 JM Dec 17 2013
             0.42, 0.50, 0.30, &  !pft 5 was 0.25
             0.80, 0.80, 0.00, &
             0.10, 0.10, 0.00 ]

epsilons = [ 0.40, 0.34, 0.00, &
             0.18, 0.35, 0.10, & !pft 3 was 0.05
             0.15, 0.15, 0.00, &
             0.00, 0.00, 0.00 ]

epsilonr = [ 0.41, 0.21, 0.00, &  !pft 2 was 0.89 JM Dec 17 2013
             0.40, 0.15, 0.60, &  !pft 5 was 0.65
             0.05, 0.05, 0.00, &
             0.90, 0.90, 0.00 ]

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


! ctem.f parameters: ----------

grescoef = [ 0.15, 0.15, 0.00, &
             0.15, 0.15, 0.15, &
             0.15, 0.15, 0.00, &
             0.15, 0.15, 0.00 ]

humicfac = [ 0.42, 0.42, 0.00, &
             0.53, 0.48, 0.48, &
             0.42, 0.42, 0.00, &
             0.42, 0.42, 0.00 ]

laimin = [ 1.0, 1.0, 0.0, &
           1.5, 1.0, 0.7, &  ! flag test PFT 4 was 1.5
           1.0, 1.0, 0.0, &
           0.1, 0.1, 0.0 ]  ! flag test PFT 8&9 were 0.5

laimax = [ 6.0, 3.0, 0.0, &  ! flag test pft 1 was 3.0
           6.0, 5.0, 5.0, & 
           8.0, 8.0, 0.0, &
           4.0, 4.0, 0.0 ] 

lambdamax = 0.10

! mainres.f parameters: ---------

bsrtstem = [ 0.0900, 0.0550, 0.0000, &
             0.0600, 0.0335, 0.0300, & 
             0.0365, 0.0365, 0.0000, & 
             0.0000, 0.0000, 0.0000 ] ! no stem component for grasses

bsrtroot = [ 0.5000, 0.2850, 0.0000, &
             0.6500, 0.2250, 0.0550, & 
             0.1600, 0.1600, 0.0000, & 
             0.1000, 0.1000, 0.0000 ]

minlvfr = 0.05

! mortality.f parameters: ---------

kmort1 = 0.3

mxmortge = [ 0.01, 0.01, 0.00, &
             0.01, 0.01, 0.01, & 
             0.00, 0.00, 0.00, & 
             0.00, 0.00, 0.00 ]

maxage = [ 250.0, 250.0,   0.0, &
           250.0, 250.0, 250.0, &
             0.0,   0.0,   0.0, &
             0.0,   0.0,   0.0 ]

! phenology.f parameters: ---------

eta = [ 10.0, 30.8, 0.00, &
        31.0, 50.0, 30.0, &
        7.0,  7.0, 0.00,  &
        3.0,  3.0, 0.00 ]

kappa =[ 1.6, 1.6, 0.0, &
         1.6, 1.6, 1.6, &
         1.6, 1.6, 0.0, &
         1.2, 1.2, 0.0 ]

lfespany  =   [ 5.0, 1.00, 0.00, &
                1.75, 1.00, 1.00, &
                1.75, 1.75, 0.00, &
                1.00, 1.00, 0.00 ]

specsla =[ 11.0, 0.0, 0.0, &
            0.0, 0.0, 0.0, &
            0.0, 0.0, 0.0, &
            0.0, 0.0, 0.0 ]

fracbofg = 0.55

cdlsrtmx = [ 0.15, 0.30, 0.00, &
             0.30, 0.15, 0.15, &
             0.15, 0.15, 0.00, &
             0.15, 0.15, 0.00 ]

drlsrtmx = [ 0.0025, 0.005, 0.000, &
             0.005, 0.005, 0.025, & !pft 5 was 0.005
             0.005, 0.005, 0.000, &
             0.025, 0.025, 0.000 ]  !PFT 8 and 9 were 0.05 

drgta = [ 3.0, 3.0, 0.0, &
          3.0, 3.0, 3.0, &
          3.0, 3.0, 0.0, &
          3.0, 3.0, 0.0 ]   !pft 8 and 9 were 3

colda = [ 3.0, 3.0, 0.0, &
          3.0, 3.0, 3.0, &
          3.0, 3.0, 0.0, &
          3.0, 3.0, 0.0 ]

lwrthrsh = [ -45.0, -5.0, 0.0, &
               5.0,  5.0, 5.0, &
               5.0,  5.0, 0.0, &
               0.1,  5.0, 0.0 ]

dayschk = [ 7, 7, 0, &
            7, 7, 7, &
            7, 7, 0, &
            7, 7, 0 ]

coldlmt = [ 7 , 5 ]   ! days

coldthrs = [ -5.0 , 8.0 ]

harvthrs = [ 0.0, 0.0, 0.0, &
             0.0, 0.0, 0.0, &
             4.5, 3.5, 0.0, &
             0.0, 0.0, 0.0 ]

flhrspan = [ 17.0, 45.0 ]

thrprcnt = [ 40.0, 40.0,  0.0, &
             40.0, 50.0, 50.0, &
             50.0, 50.0,  0.0, &
             40.0, 40.0,  0.0 ]

roothrsh = 15.0


! turnover.f parameters: --------------

stemlife = [ 65.0, 75.0, 0.00, &
             45.0, 70.0, 65.0, &  !flag pft 4 was 40.0, pft 5 was 45.0
             20.0, 20.0, 0.00, &
              0.00, 0.00, 0.00 ]

rootlife = [ 10.0,11.5, 0.0, &
              8.5, 9.5, 8.5, &     ! flag pft 4 and 5 were 5.5
              3.0, 3.0, 0.0, &
              2.5, 2.5, 0.0 ]

stmhrspn = 17.0

end subroutine initpftpars

end module ctem_params
