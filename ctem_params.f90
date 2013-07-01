module ctem_params

! This module holds ctem globally accessible parameters

! J. Melton
! Jun 23 2013

implicit none

! constants

real, parameter :: zero = 1.0e-20
real, parameter :: pi =3.1415926535898d0
real, parameter :: earthrad = 6371.22   ! radius of earth, km

integer, parameter, dimension(12) :: monthdays = [ 31,28,31,30,31,30,31,31,30,31,30,31 ] ! days in each month
integer, parameter, dimension(13) :: monthend = [ 0,31,59,90,120,151,181,212,243,273,304,334,365 ] ! calender day at end of each month

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


! model state
integer, parameter :: nlat=1            !
integer, parameter :: nmos=10           !
integer, parameter :: ilg = nlat*nmos   !
integer, parameter :: nmon=12           !
integer, parameter :: ican=4            !
integer, parameter :: ignd=3            !
integer, parameter :: icp1 = ican + 1   !

integer,parameter :: icc=9              !
integer,parameter :: iccp1 = icc + 1    !


! plant-related

integer, parameter :: kk = 12  ! product of class pfts and l2max
integer, parameter :: numcrops = 2    ! number of crop pfts
integer, parameter :: numtreepfts = 5 ! number of tree pfts
integer, parameter :: numgrass = 2    ! number of grass pfts

real, parameter :: seed = 0.001 ! seed pft fraction, same as in competition
                                ! in mosaic mode, all tiles are given this
                                ! as a minimum



end module ctem_params
