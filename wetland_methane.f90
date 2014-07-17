subroutine wetland_methane (hetrores, il1, il2, ilg, ta, wetfrac, &
                       ignd, npp, tbar, thliqg, currlat, &
                        sand,  wetfrac_s, &! obswetf, &
!
!    -------------- inputs above this line, outputs below -------------
                       ch4wet1,    ch4wet2,    wetfdyn,  &
                      ch4dyn1,    ch4dyn2)

!               Canadian Terrestrial Ecosystem Model (CTEM)
!           	Wetland and wetland methane subroutine 

!      4  July  2014 - Convert to f90 and bring in ctem_params
!     J. Melton
!
!      9  June. 2010 - this subroutine calculates methane flux
!     Brian Amiro      from wetlands of a grid cell (i.e. land only)

use ctem_params,        only : wtdryres, ratioch4, factor2,lat_thrshld1, &
                               lat_thrshld2, soilw_thrshN, soilw_thrshE, &
                               soilw_thrshS 

implicit none

!logical, intent(in) :: obswetf                  ! if true the observed wetland fraction will be used.
integer, intent(in) :: ilg                      ! no. of grid cells in latitude circle
integer, intent(in) :: il1                      ! il1=1
integer, intent(in) :: il2                      ! il2=ilg
integer, intent(in) :: ignd                     ! no. of soil layers 
real, dimension(ilg), intent(in) :: hetrores    ! heterotrophic respiration from main ctem program calculated as sum of litres + socres
real, dimension(ilg), intent(in) :: ta          ! air temperature, k
real, dimension(ilg), intent(in) :: wetfrac     ! prescribed fraction of wetlands in a grid cell
real, dimension(ilg), intent(in) :: npp         ! grid-averaged npp from ctem (u-mol co2/m2.s)
real, dimension(ilg), intent(in) :: currlat     ! centre latitude of grid cells in degrees
real, dimension(ilg,5), intent(in) :: wetfrac_s ! prescribed fraction of wetlands based on slope only (2, 2.5, 3, 3.5, 4 %)
real, dimension(ilg,ignd), intent(in) :: tbar   ! temperature of soil layers
real, dimension(ilg,ignd), intent(in) :: thliqg ! liquid soil moisture content (fraction)
real, dimension(ilg,ignd), intent(in) :: sand   ! percentage sand in soil layers
real, dimension(ilg), intent(out) :: ch4wet1    ! methane flux from wetlands calculated using hetrores in umol ch4-c/m2.s
real, dimension(ilg), intent(out) :: ch4wet2    ! methane flux from wetlands calculated using npp in umol ch4-c/m2.s
real, dimension(ilg), intent(out) :: wetfdyn    ! dynamic gridcell wetland fraction determined using  slope and soil moisture
real, dimension(ilg), intent(out) :: ch4dyn1    ! methane flux from wetlands calculated using hetrores and wetfdyn, in umol ch4-c/m2.s
real, dimension(ilg), intent(out) :: ch4dyn2    ! methane flux from wetlands calculated using npp  and wetfdyn, in umol ch4-c/m2.s
    
! local variables  
real, dimension(ilg) :: wetresp
real :: porosity
real :: soil_wetness
integer :: i

!     ---------------------------------------------------------------
!     Constants and parameters are located in ctem_params.f90
!     -----------------------------------------------------------------

!     initialize required arrays to zero

do 110 i = il1, il2
        wetresp(i)=0.0     ! heterotrophic wetland respiration
        ch4wet1(i)=0.0      ! methane flux from wetlands
        ch4wet2(i)=0.0
        ch4dyn1(i)=0.0
        ch4dyn2(i)=0.0
110   continue

!     initialization ends  
!	--------------------------------------------------
!
!	Estimate the methane flux from wetlands for each grid cell
!	scaling by the wetland fraction in a grid cell
! 	and set the methane flux to zero when screen temperature (ta) is below or at freezing
!	this is consistent with recent flux measurements by the university of manitoba at churchill, manitoba
!
!if (obswetf) then  ! Use the read-in wetland locations as the CH4 producing area

   do 210 i = il1, il2 
      wetresp(i)=hetrores(i)*wtdryres*wetfrac(i)
      ch4wet1(i)=ratioch4*wetresp(i)
      ch4wet2(i)=factor2*wetfrac(i)*max(0.0,npp(i))*(2**((tbar(i,1)-273.2)/10.0)) 
      if (ta(i).lt.273.2) then
         ch4wet1(i)=0.0
         ch4wet2(i)=0.0
      endif
210 continue

!else ! dynamically find the wetland locations

   do 310 i = il1, il2 
     porosity=(-0.126*sand(i,1)+48.9)/100.0 ! top soil layer porosity
     soil_wetness=(thliqg(i,1)/porosity)
     soil_wetness=max(0.0,min(soil_wetness,1.0))

!    if soil wetness meets a latitude specific threshold then the slope
!    based wetland fraction is wet and is an actual wetland else not
!
     wetfdyn(i)=0.0  ! initialize dynamic wetland fraction to zero
!

! Original way:
!     if (currlat(i).ge.lat_thrshld1) then ! all area north of 35 n
!        if(soil_wetness.gt.soilw_thrshN)then ! ckw 3041-3050
!           wetfdyn(i)=wetfrac_s(i)
!        endif
!     elseif (currlat(i).lt.lat_thrshld1.and.currlat(i).ge. lat_thrshld2) then ! between 10 s and 35 n
!        if(soil_wetness.gt.soilw_thrshE)then   ! ckw 3041-3050
!            wetfdyn(i)=wetfrac_s(i)
!        endif
!     else ! everything else below 10 s
!        if(soil_wetness.gt.soilw_thrshS)then ! ckw 3041-3050
!            wetfdyn(i)=wetfrac_s(i)
!        endif
!     endif    

! Testing:  

!     if (currlat(i).ge.lat_thrshld1) then ! all area north of 35 n
        if(soil_wetness .gt. 0.55 .and. soil_wetness .le. 0.75) then 
           wetfdyn(i)=wetfrac_s(i,1)  !<2% slope class
!        endif
!     elseif (currlat(i).lt.lat_thrshld1.and.currlat(i).ge. lat_thrshld2) then ! between 10 s and 35 n
        elseif (soil_wetness .gt. 0.75 .and. soil_wetness .lt. 00.1) then
            wetfdyn(i)=wetfrac_s(i,3) !<3 % slope class
!        endif
!     else ! everything else below 10 s
        elseif (soil_wetness .eq. 1.00) then 
            wetfdyn(i)=wetfrac_s(i,5) ! <4% slope class
        endif
!     endif    
     

! new dynamic calculation
! same as ch4wet1 & 2, but wetfrac replaced by wetfdyn
           wetresp(i)=hetrores(i)*wtdryres*wetfdyn(i)
           ch4dyn1(i)=ratioch4*wetresp(i)
           ch4dyn2(i)=factor2*wetfdyn(i)*max(0.0,npp(i))*(2**((tbar(i,1)-273.2)/10.0))
           if (ta(i).lt.273.2) then
             ch4dyn1(i)=0.0
             ch4dyn2(i)=0.0
           endif

310 continue

!end if !obswetf

return
end
 




