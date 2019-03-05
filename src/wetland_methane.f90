

!>\file
!>Prognostic determination of wetland fractional coverage and emission of wetland methane subroutine 
!!@author B. Amiro, J. Melton, V. Arora
!>
subroutine wetland_methane (hetrores, il1, il2, ta, wetfrac, &
                        thliqg, currlat, &
                        sand,  slopefrac, &  
! inputs above this line, outputs below -------------
                       ch4WetSpec,  wetfdyn, ch4WetDyn)

!     20  Feb   2019 - Update code to allow namelist of params. Remove no CH4 flux if soils are
!     J. Melton        frozen condition.
!     31  Aug   2016 - Change how we find wetlands from discrete limits to
!     V. Arora         smooth function
!
!      4  July  2014 - Convert to f90 and bring in classic_params
!     J. Melton
!
!      9  June. 2010 - this subroutine calculates methane flux
!     Brian Amiro      from wetlands of a grid cell (i.e. land only)

use classic_params,        only : wtdryres, ratioch4,lat_thrshld1, &
                               lat_thrshld2, soilw_thrshN, soilw_thrshE, &
                               soilw_thrshS, ilg, ignd

implicit none

integer, intent(in) :: il1                      !< il1=1
integer, intent(in) :: il2                      !< il2=ilg
real, dimension(ilg), intent(in) :: hetrores    !< heterotrophic respiration from main ctem program calculated as sum of litres + socres
real, dimension(ilg), intent(in) :: ta          !< air temperature, k
real, dimension(ilg), intent(in) :: wetfrac     !< prescribed fraction of wetlands in a grid cell
real, dimension(ilg), intent(in) :: currlat     !< centre latitude of grid cells in degrees
real, dimension(ilg,8), intent(in) :: slopefrac !< Fraction of gridcell flatter than slope thresholds (0.025, 0.05, 0.1, 0.15, 0.20, 0.25, 0.3 and 0.35 percent slope thresholds)
real, dimension(ilg,ignd), intent(in) :: thliqg !< liquid soil moisture content (fraction)
real, dimension(ilg,ignd), intent(in) :: sand   !< percentage sand in soil layers
real, dimension(ilg), intent(out) :: ch4WetSpec    !< methane flux from wetlands calculated using hetrores in umol ch4/m2.s
real, dimension(ilg), intent(out) :: wetfdyn    !< dynamic gridcell wetland fraction determined using  slope and soil moisture
real, dimension(ilg), intent(out) :: ch4WetDyn    !< methane flux from wetlands calculated using hetrores and wetfdyn, in umol ch4/m2.s
    
! local variables  
real, dimension(ilg) :: wetresp !<heterotrophic wetland respiration
real :: upvslow    ! ratio of wetland to upland respiration for this location
real :: porosity
real :: soil_wetness
integer :: i

real, parameter :: alpha = 0.45  ! Determines shape of curve.
real :: low_mois_lim
real :: mid_mois_lim
real :: upp_mois_lim
real :: x1
real :: x2
real :: y1
real :: y2
real :: slope
real :: intercept

!>
!>---------------------------------------------------------------
!>Constants and parameters are located in classic_params.f90
!>-----------------------------------------------------------------
!>
!>initialize required arrays to zero
!>
wetresp(:)=0.0     
wetfdyn(:)=0.0  
ch4WetSpec(:)=0.0     
ch4WetDyn(:)=0.0
!>
!>--------------------------------------------------
!>
!>Estimate the methane flux from wetlands for each grid cell scaling by the wetland fraction in a grid cell

! Set up the latitude bounds based on the paramters read in from the namelist file.
! If soil wetness meets a latitude specific threshold then the slope based wetland
!fraction is wet and is an actual wetland, else not. As well the ratio of upland to
! wetland respiration is assumed to be latitidionally varying. This is intended to
! reflect the differnt wetland types in the tropics (floodplain types) vs. those in
! higher latitudes (peatlands). This is a very coarse approximation and should be 
! replaced once we have CH4 in our peatland module. The choice of wtdryres instead
! of ratioCH4 to try and mimic this is arbitrary, either parameter could be used since
! they are simply multiplicative.

     if (currlat(i) >= lat_thrshld1) then ! Northern high lats, all area north of lat_thrshld1
        low_mois_lim = soilw_thrshN(1) 
        mid_mois_lim = soilw_thrshN(2) 
        !upp_mois_lim = soilw_thrshN(3)
        upvslow = wtdryres(1)
     elseif (currlat(i) < lat_thrshld1 .and. currlat(i) >= lat_thrshld2) then ! Tropics
        low_mois_lim = soilw_thrshE(1) 
        mid_mois_lim = soilw_thrshE(2) 
        !upp_mois_lim = soilw_thrshE(3) 
        upvslow = wtdryres(2)
     else ! S. Hemi,  everything else below lat_thrshld2
        low_mois_lim = soilw_thrshS(1) 
        mid_mois_lim = soilw_thrshS(2) 
        !upp_mois_lim = soilw_thrshS(3) 
        upvslow = wtdryres(3)
      end if

!> First calculate for the specified wetland fractions read in from OBSWETFFile
   do 210 i = il1, il2 
      wetresp(i) = hetrores(i) * upvslow * wetfrac(i)
      ch4WetSpec(i) = ratioch4 * wetresp(i)
210 continue
!>
!>Next dynamically find the wetland locations and determine their methane emissions
!>
   do 310 i = il1, il2
     porosity = (-0.126 * sand(i,1) + 48.9) / 100.0 ! top soil layer porosity
     soil_wetness = (thliqg(i,1) / porosity)
     soil_wetness = max(0.0, min(soil_wetness, 1.0))


!    implement Vivek's new way of modelling WETFDYN
     
     x1 = low_mois_lim * (1. - alpha) + mid_mois_lim * alpha
     x2 = 1.
     y1 = 0.
     y2 = slopefrac(i,5)
     slope = (y2 - y1) / (x2 - x1)
     intercept = slope * x1 * (-1)

     wetfdyn(i) = min(1.0, max(0.0, slope * soil_wetness + intercept))

!    new dynamic calculation
!    same as ch4WetSpec & 2, but wetfrac replaced by wetfdyn

     wetresp(i) = hetrores(i) * upvslow * wetfdyn(i)
     ch4WetDyn(i) = ratioch4 * wetresp(i)

310 continue

return
end
