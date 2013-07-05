module ctem_params

! This module holds ctem globally accessible parameters

! J. Melton
! Jun 23 2013

implicit none

! constants

real, parameter :: zero = 1.0e-20
real, parameter :: pi =3.1415926535898d0
real, parameter :: earthrad = 6371.22   ! radius of earth, km
real, parameter :: deltat =1.0  !CTEM's time step in days
real, parameter :: abszero=1e-12 ! this one is from RUNCLASSCTEM, replace with zero?? JM Jul 2 2013.

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

! conversion factor from carbon to dry organic matter value is from li et al. 2012 biogeosci
real, parameter :: c2dom = 450.0 ! gc / kg dry organic matter

! ----
! model state
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

! ----
! PFT parameters

! eta and kappa, parameters for estimating min. stem+root biomass
! required to support green leaf biomass. kappa is 1.6 for trees
! and crops, and 1.2 for grasses.
real, parameter, dimension(kk) :: eta =[ 10.0, 30.8, 0.00, &
                                         31.0, 50.0, 30.0, &
                                          7.0,  7.0, 0.00, &
                                          3.0,  3.0, 0.00 ]

real, parameter, dimension(kk) :: kappa =[ 1.6, 1.6, 0.0, &
                                           1.6, 1.6, 1.6, &
                                           1.6, 1.6, 0.0, &
                                           1.2, 1.2, 0.0 ]

! leaf life span (in years) for ctem's 9 pfts
real, parameter, dimension(kk) :: lfespany =[ 5.0, 1.00, 0.00, &
                                             1.75, 1.00, 1.00, &
                                             1.75, 1.75, 0.00, &
                                             1.00, 1.00, 0.00 ]

! CTEM can use user-specified specific leaf areas (sla) if the
! following specified values are greater than zero
real, parameter, dimension(kk) :: specsla =[ 11.0, 0.0, 0.0, &
                                              0.0, 0.0, 0.0, &
                                              0.0, 0.0, 0.0, &
                                              0.0, 0.0, 0.0 ]

! canopy light/nitrogen extinction coefficient 
real, parameter, dimension(kk) :: kn= [ 0.50, 0.50, 0.00, &
                                        0.50, 0.50, 0.50, &
                                        0.40, 0.48, 0.00, &
                                        0.46, 0.44, 0.00 ]

! fracbofg, parameter used to estimate lai of brown leaves. we
! assume that sla of brown leaves is this fraction of sla of
! green leaves
real, parameter :: fracbofg = 0.55


end module ctem_params
