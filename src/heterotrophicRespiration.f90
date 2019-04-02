!>\file
!>Central module for all heterotrophic respiration-related operations
module heterotrophic_respiration

implicit none

! Subroutines contained in this module:
public  :: hetresg
public  :: hetresv

contains

  !>\ingroup heterotrophic_respiration_hetresg
  !!@{
  !>Heterotrophic respiration subroutine for bare ground fraction
  !> @author Vivek Arora and Joe Melton
  subroutine hetresg (litrmass,  soilcmas,   delzw,    thpor,    & !In
                           il1,       il2,     ilg,     tbar,    & !In
                         psisat,        b,   thliq,   & !In
                         thice,     frac,   isand,    & !In
                         litres,   socres) ! Out

  !               Canadian Terrestrial Ecosystem Model (CTEM)
  !           Heterotrophic Respiration Subroutine For Bare Fraction
  !
  !     11  Apr. 2003 - this subroutine calculates heterotrophic respiration
  !     V. Arora        over the bare subarea of a grid cell (i.e. ground only
  !                     and snow over ground subareas).
  !
  !     change history:

  !      7  Jun 2016  - Bring in step function for reduction in resp when soil freezes.
  !     J. Melton
  !
  !      6  Jun 2016  - Het resp reduces at depth in soil column.
  !     J. Melton
  !
  !      8  Feb 2016  - Adapted subroutine for multilayer soilc and litter (fast decaying)
  !     J. Melton       carbon pools
  !
  !     19  Jan 2016  - Implemented new LUC litter and soil C pools
  !     J. Melton
  !
  !     14  Jan 2016  - Converted to F90 and made it so it can handle > 3 soil layers
  !     J. Melton
  !
  !     30  Jul 2015  - Based on work by Yuanqiao Wu, respiration was found to
  !     J. Melton       behave incorrectly if the soil froze as it thought the water
  !                     was leaving the soil. This is now fixed.
  !
  !     17  Jan 2014  - Moved parameters to global file (classic_params.f90)
  !     J. Melton
  !
  !     23  Jul 2013  - add in module for parameters
  !     J. Melton
  !     J. Melton and V.Arora - changed tanhq10 parameters, they were switched
  !               25 Sep 2012
  !     J. Melton 31 Aug 2012 - remove isnow, it is not used.
  !     J. Melton 23 Aug 2012 - bring in isand, converting sand to
  !                             int was missing some gridcells assigned
  !                             to bedrock in classb

  use classic_params,        only : icc, ignd, zero, tanhq10, a_hetr, &
                                    bsratelt_g, bsratesc_g, r_depthredu, &
                                    tcrit,frozered

  implicit none

  integer, intent(in) :: ilg   !<
  integer, intent(in) :: il1   !<il1=1
  integer, intent(in) :: il2   !<il2=ilg
  integer, intent(in) :: isand(:,:) !<
  real, intent(in) :: litrmass(:,:)   !<litter mass for the 8 pfts + bare in \f$kg c/m^2\f$
  real, intent(in) :: soilcmas(:,:)   !<soil carbon mass for the 8 pfts + bare in \f$kg c/m^2\f$
  real, intent(in) :: tbar(:,:)    !<soil temperature, k
  real, intent(in) :: thliq(:,:)   !<liquid soil moisture content in 3 soil layers
  real, intent(in) :: frac(:)         !<fraction of bare ground (fg)
  real, intent(in) :: delzw(:,:)   !<
  real, intent(in) :: thice(:,:)   !<
  real, intent(in) :: psisat(:,:)  !<saturation matric potential
  real, intent(in) :: b(:,:)       !<parameter b of clapp and hornberger
  real, intent(in) :: thpor(:,:)   !<porosity
      
  real, intent(out) :: litres(ilg,ignd)      !<litter respiration over the given unvegetated sub-area in umol co2/m2.s
  real, intent(out) :: socres(ilg,ignd)      !<soil c respiration over the given unvegetated sub-area in umol co2/m2.s
  
  ! Local 
  real :: litrq10          !<
  real :: soilcq10         !<
  real :: q10func          !<
  real :: tempq10l         !<
  real :: socmoscl(ilg,ignd)    !<soil moisture scalar for soil carbon decomposition
  real :: ltrmoscl(ilg,ignd)    !<soil moisture scalar for litter decomposition
  real :: psi(ilg,ignd)    !<
  real :: tempq10s         !<
  real :: reduceatdepth
  integer :: i,j,k

  ! -------------------------------------------------------------------------
  
  ! initialize required arrays to zero

  socmoscl(:,:)=0.0   ! soil moisture scalar for soil carbon decomposition
  ltrmoscl(:,:)=0.0   ! soil moisture scalar for litter decomposition
  litres(:,:)=0.0     ! litter resp. rate
  socres(:,:)=0.0     ! soil c resp. rate

  !     initialization ends

  !> Find moisture scalar for soil c decomposition
  !! this is modelled as function of logarithm of matric potential.
  !! we find values for all soil layers, and then find an average value
  !! based on fraction of carbon present in each layer.

  do 260 j = 1, ignd
    do 270 i = il1, il2

      if (isand(i,j) == -3 .or. isand(i,j) == -4) then
        socmoscl(i,j)=0.2
        psi (i,j) = 10000.0 ! set to large number so that ltrmoscl becomes 0.2
      else ! i.e., sand.ne.-3 or -4

        ! We don't place a lower limit on psi as it is only used here and the
        ! value of psi <= psisat is just used to select a scomotrm value.
        if (thice(i,j) <= thpor(i,j)) then ! flag new limits
          psi(i,j) = psisat(i,j) * (thliq(i,j) / (thpor(i,j) - thice(i,j)))**(-b(i,j)) 
        else      
          ! if the whole pore space is ice then suction is assumed to be very high.   
          psi(i,j) = 10000.0
        end if 

        if (psi(i,j) >= 10000.0) then
          socmoscl(i,j) = 0.2
        else if (psi(i,j) < 10000.0 .and.  psi(i,j) > 6.0) then
          socmoscl(i,j) = 1.0 - 0.8 * ( (log10(psi(i,j)) - log10(6.0)) &
                                      / (log10(10000.0) - log10(6.0)) )
        else if (psi(i,j) <= 6.0 .and. psi(i,j) >= 4.0) then
          socmoscl(i,j) = 1.0
        else if (psi(i,j) < 4.0 .and. psi(i,j) > psisat(i,j)) then
          socmoscl(i,j) = 1.0 - 0.5 * ((log10(4.0) - log10(psi(i,j))) &
                                     / (log10(4.0) - log10(psisat(i,j))) )
        else if ( psi(i,j) <= psisat(i,j)) then
          socmoscl(i,j) = 0.5
        end if
      end if 

      socmoscl(i,j) = max(0.2, min(socmoscl(i,j), 1.0))

270   continue
260 continue

  !> Find moisture scalar for litter decomposition for top soil layer
  !> The difference between moisture scalar for litter and soil c
  !! is that the litter decomposition is not constrained by high
  !! soil moisture (assuming that litter is always exposed to air).
  !! In addition, we use moisture content of the top soil layer
  !! as a surrogate for litter moisture content. So we use only
  !! psi(i,1) calculated in loops 260 and 270 above.

  do 300 i = il1, il2
    if (psi(i,1) > 10000.0) then
      ltrmoscl(i,1) = 0.2
    else if (psi(i,1) <= 10000.0 .and.  psi(i,1) > 6.0 ) then
      ltrmoscl(i,1) = 1.0 - 0.8 * ((log10(psi(i,1)) - log10(6.0)) &
                                  /(log10(10000.0) - log10(6.0)) )
    else if (psi(i,1) <= 6.0) then
      ltrmoscl(i,1) = 1.0
    end if
    ltrmoscl(i,1) = max(0.2, min(ltrmoscl(i,1), 1.0))
300   continue

  !> Treat the lower litter layers like the soil C ones:
  ltrmoscl(:,2:ignd) = socmoscl(:,2:ignd)

!> Use temperature of the litter and soil c pools, and their soil
!! moisture scalars to find respiration rates from these pools

  do 330 i = il1, il2
    do 340 j = 1, ignd
      if(frac(i).gt.zero)then

        !! First find the q10 response function to scale base respiration
        !! rate from 15 c to current temperature, we do litter first
        tempq10l = tbar(i,j) - 273.16
        litrq10 = tanhq10(1) + tanhq10(2) * (tanh(tanhq10(3) * (tanhq10(4) - tempq10l)))
        
        !! Apply a step reduction in q10func when the soil layer freezes. This is to reflect the
        !! lost mobility of the population due to less liquid water
        if (tbar(i,j) - 273.16 > tcrit) then
          q10func = litrq10**(0.1* (tbar(i,j) - 273.16 - 15.0))
        else
          q10func = litrq10**(0.1 * (tbar(i,j) - 273.16 - 15.0)) * frozered
        end if

        !> Reduce the respiration at depth due to unresolved depth dependent processes including
        !! soil microbial population dynamics, pore-scale oxygen availability, mineral sorption,
        !! priming effects, and other unrepresented processes. This is following Lawrence et al.
        !! Enviro Res Lett 2015 \cite Lawrence2015-tj. We apply this for all soils.
        reduceatdepth = exp(-delzw(i,j) / r_depthredu)
        ! 2.64 converts bsratelt_g from kg c/kg c.year to u-mol co2/kg c.s
        litres(i,j) = ltrmoscl(i,j) * litrmass(i,j) * bsratelt_g * 2.64 &
                      * q10func * reduceatdepth 

        ! Respiration from soil c pool
        tempq10s = tbar(i,j) - 273.16
        soilcq10 = tanhq10(1) + tanhq10(2) * (tanh(tanhq10(3) * (tanhq10(4) - tempq10s)))
        
        if (tbar(i,j) - 273.16 > tcrit) then
            q10func = soilcq10**(0.1 * (tbar(i,j) - 273.16 - 15.0))
        else
            q10func = soilcq10**(0.1 * (tbar(i,j) - 273.16 - 15.0)) * frozered
        end if
        ! 2.64 converts bsratelt_g from kg c/kg c.year to u-mol co2/kg c.s
        socres(i,j) = socmoscl(i,j) * soilcmas(i,j) * bsratesc_g * 2.64 * q10func * reduceatdepth 
      endif

340  continue
330 continue

  return
  
end subroutine hetresg
!>@}
! ------------------------------------------------------------------------------------

!>\ingroup heterotrophic_respiration_hetresv
!!@{
!>Heterotrophic Respiration Subroutine For Vegetated Fraction
!> @author Vivek Arora, Joe Melton, Yuanqiao Wu
  subroutine hetresv ( fcan,      fct,   litrmass, soilcmas,  & !In
                      delzw,    thpor,        il1,      il2,  & !In
                        ilg,     tbar,     psisat,    thliq,  & !In
                       sort,        b,  & !In
                       isand,  thice,  ipeatland,            & !In
                    ltresveg, scresveg) !Out

  !               Canadian Terrestrial Ecosystem Model (CTEM)
  !           Heterotrophic Respiration Subtoutine For Vegetated Fraction

  !     16  oct. 2001 - this subroutine calculates heterotrophic respiration
  !     V. Arora        for a given sub-area, from litter and soil carbon
  !                     pools.

  !     change history:

  !
  !      7  Jun 2016  - Bring in step function for reduction in resp when soil freezes.
  !     J. Melton
  !
  !      6  Jun 2016  - Het resp reduces at depth in soil column.
  !     J. Melton
  !
  !      8  Feb 2016  - Adapted subroutine for multilayer soilc and litter (fast decaying)
  !     J. Melton       carbon pools
  !
  !     14  Jan 2016  - Converted to F90 and made it so it can handle > 3 soil layers
  !     J. Melton

  !     10  April 2015 -Bring in peatland scheme
  !     Y. Wu
  !
  !     30  Jul 2015  - Based on work by Yuanqiao Wu, respiration was found to
  !     J. Melton       behave incorrectly if the soil froze as it thought the water
  !                     was leaving the soil. This is now fixed.
  !     17  Jan 2014  - Moved parameters to global file (classic_params.f90)
  !     J. Melton

  !     22  Jul 2013  - Add in module for parameters
  !     J. Melton

  !     J. Melton and V.Arora - changed tanhq10 parameters, they were switched
  !               25 sep 2012
  !     J. Melton 23 aug 2012 - bring in isand, converting sand to
  !                             int was missing some gridcells assigned
  !                             to bedrock in classb
  !     ------

  use classic_params,        only : icc, ignd, kk, zero, bsratelt,&
                           bsratesc, abar, tanhq10,&
                           alpha_hetres, r_depthredu,tcrit,frozered

  implicit none

  ! Arguments:
  integer, intent(in) :: ilg                              !< il1=1
  integer, intent(in) :: il1                              !< il1=1
  integer, intent(in) :: il2                              !< il2=ilg
  real, dimension(ilg,icc,ignd), intent(in) :: litrmass   !< litter mass for the 9 pfts + bare [ \f$kg C/m^2\f$ ]
  real, dimension(ilg,icc,ignd), intent(in) :: soilcmas   !< soil carbon mass for the 9 pfts + bare [ \f$kg C/m^2\f$ ]
  real, dimension(ilg,ignd), intent(in) :: thpor          !< Soil total porosity [ \f$(cm^3 cm^{-3})\f$ ] - daily average
  real, dimension(ilg,ignd), intent(in) :: tbar           !< Soil temperature [ K ]
  real, dimension(ilg,ignd), intent(in) :: psisat         !< Saturated soil matric potential [ m ]
  real, dimension(ilg,ignd), intent(in) :: b              !< Clapp and Hornberger empirical “b” parameter [ ]
  real, dimension(ilg,ignd), intent(in) :: thliq          !< liquid soil moisture content in soil layers [ \f$(cm^3 cm^{-3})\f$ ]
  real, dimension(ilg,ignd), intent(in) :: thice         !< frozen soil moisture content in 3 soil layers in canopy covered subarea [ \f$(cm^3 cm^{-3})\f$ ]
  real, dimension(ilg), intent(in) :: fct                 !< Sum of all fcans [ ]
  real, dimension(ilg,icc), intent(in) :: fcan            !< fractional coverage of ctem's 9 pfts [ ]
  integer, dimension(ilg,ignd), intent(in) :: isand       !< flag for soil/bedrock/ice/glacier
  integer, dimension(icc), intent(in) :: sort             !< index for correspondence between CTEM pfts and all values in the parameters vectors
  real, dimension(ilg,ignd), intent(in) :: delzw          !< thickness of the permeable soil layers [ m ]
  integer, intent(in) ::  ipeatland(ilg)                  !<

  real, dimension(ilg,icc,ignd), intent(out) :: ltresveg  !< litter respiration for the given vegetated sub-area [ \f$u-mol co2/m2.sec\f$ ]
  real, dimension(ilg,icc,ignd), intent(out) :: scresveg  !< soil carbon respiration for the given vegetated sub-area [ \f$u-mol co2/m2.sec\f$ ]

  ! Local vars:
  integer i, j, k
  real :: litrq10
  real :: soilcq10
  real :: q10func
  real :: tempq10l
  real :: tempq10s
  real :: reduceatdepth
  real, dimension(ilg,ignd) :: socmoscl
  real, dimension(ilg,ignd) :: psi
  real, dimension(ilg,ignd) :: ltrmoscl

  ! Initialize required arrays to zero
  socmoscl(:,:)=0.0      ! soil carbon moisture term
  ltrmoscl(:,:)=0.0     ! soil moisture scalar for litter decomposition
  ltresveg(:,:,:)=0.0   ! litter resp. rate for each pft
  scresveg(:,:,:)=0.0   ! soil c resp. rate for each pft

  !! Find moisture scalar for soil c decomposition
  !! This is modelled as function of logarithm of matric potential.
  !! we find values for all soil layers, and then find an average value
  !! based on fraction of carbon present in each layer. this makes
  !! moisture scalar a function of vegetation type.

  do 260 j = 1, ignd
    do 270 i = il1, il2
      if(isand(i,j).eq.-3.or.isand(i,j).eq.-4)then
        socmoscl(i,j)=0.2
        psi (i,j) = 10000.0 ! set to large number so that ltrmoscl becomes 0.2
      else ! i.e., sand.ne.-3 or -4

        if (ipeatland(i) > 0) then
          if (thliq(i,j) + thice(i,j) + 0.01 < thpor(i,j) .and. tbar(i,j) < 273.16) then
            psi(i,j) = 0.001
          elseif (thice(i,j) > thpor(i,j)) then
            psi(i,j) = 0.001   !set to saturation
          else
            psi(i,j) = psisat(i,j) * (thliq(i,j) / (thpor(i,j) - thice(i,j)))**(-b(i,j))
          end if
        else
          if (thice(i,j) <= thpor(i,j)) then ! flag new limits
            psi(i,j) = psisat(i,j) * (thliq(i,j) / (thpor(i,j) - thice(i,j)))**(-b(i,j)) 
          else 
            ! if the whole pore space is ice then suction is assumed to be very high. 
            psi(i,j) = 10000.0
          end if 
        end if

        if(psi(i,j) >= 10000.0) then
          socmoscl(i,j) = 0.2
        else if (psi(i,j) < 10000.0 .and.  psi(i,j) > 6.0) then
          socmoscl(i,j) = 1.0 - 0.8 * ((log10(psi(i,j)) - log10(6.0)) &
                                     / (log10(10000.0) - log10(6.0)))
        else if (psi(i,j) <= 6.0 .and. psi(i,j) >= 4.0) then
          socmoscl(i,j) = 1.0
        else if (psi(i,j) < 4.0 .and. psi(i,j) > psisat(i,j)) then
          socmoscl(i,j) = 1.0 -0.5 * ((log10(4.0) - log10(psi(i,j)))&
                                   / (log10(4.0) - log10(psisat(i,j))))
        else if (psi(i,j) <= psisat(i,j)) then
          socmoscl(i,j) = 0.5
        end if
      end if ! sand.eq.-3 or -4

      socmoscl(i,j) = max(0.2, min(1.0, socmoscl(i,j)))

270     continue
260   continue

  !! Find moisture scalar for litter decomposition
  !! the difference between moisture scalar for litter and soil c
  !! is that the litter decomposition is not constrained by high
  !! soil moisture (assuming that litter is always exposed to air).
  !! in addition, we use moisture content of the top soil layer
  !! as a surrogate for litter moisture content. so we use only
  !! psi(i,1) calculated in loops 260 and 270 above.

  !  FLAG right now assume that the first soil layer litter behaves like this.
  !  the layers below are more impeded by soil moisture (same as soil C). JM Feb 8 2016.
  do 300 i = il1, il2
    if (ipeatland(i) == 0) then   !not peatland
      if (psi(i,1) > 10000.0) then
        ltrmoscl(i,1) = 0.2
      else if (psi(i,1) <= 10000.0 .and. psi(i,1) > 6.0) then
        ltrmoscl(i,1) = 1.0 - 0.8 *((log10(psi(i,1)) - log10(6.0)) &
                                  / (log10(10000.0) - log10(6.0)))
      else if (psi(i,1) <= 6.0) then
        ltrmoscl(i,1) = 1.0
      end if
      ltrmoscl(i,1) = max(0.2, min(1.0, ltrmoscl(i,1)))
    else  !is peatland
      !  FLAG per layer implementation untested!! JM Mar 14 2019.
      !    test psi optimal at psisat
      !    peatland microbals performs better towards wet environment,
      !    for b = 2.3, thpor = 0.98 as soil layer 1,
      !    thliq = 0.01  0.1   0.2    0.3    0.4    0.5   0.6   0.7    0.8     0.9
      !    psi   =  391  1.0   0.38  0.15   0.08   0.05   0.03  0.022  0.016  0.012
      !
      !    set the upper boundary at 500, optimal psi between 0.05 and 0.03
      !    (Mayono et al. 2013)
      !
      !    limit of ltrmoscalms at saturation
      if (psi(i,1) >= 10000.0) then
        ltrmoscl(i,1) = 0.2
      else if (psi(i,1) <= 10000.0 .and. psi(i,1) > 6.0) then
        ltrmoscl(i,1) =1.0 - 0.8 * ((log10(psi(i,1)) - log10(6.0)) &
                                  / (log10(10000.0) - log10(6.0)))**1.
      else if (psi(i,1) <= 6.0 .and. psi(i,1) > 4.0) then
        ltrmoscl(i,1) = 1.0
      else if (psi(i,1) <= 4.0 .and. psi(i,1) > psisat(i,1)) then
        ltrmoscl(i,1) = 1.0 - 0.99*((log10(4.0) - log10(psi(i,1))) &
                                  / (log10(4.0) - log10(psisat(i,1))))
      else if (psi(i,1) <= psisat(i,1)) then
        ltrmoscl(i,1) = 0.01
      end if
      ltrmoscl(i,1) = max(0.0, min(ltrmoscl(i,1), 1.0))
    endif  !peatland
300   continue

  !> Set the lower levels to have the same moisture sensitivity to soil C.
  ltrmoscl(:,2:ignd)=socmoscl(:,2:ignd)

  !< Use temperature of the litter and soil c pools, and their soil
  !! moisture scalars to find respiration rates from these pools

  do 320 j = 1, icc
    do 330 i = il1, il2
      do 340 k= 1, ignd
        if (fcan(i,j) .gt. 0.) then

          ! First find the q10 response function to scale base respiration
          ! rate from 15 c to current temperature, we do litter first

          tempq10l = tbar(i,k) - 273.16
          litrq10 = tanhq10(1) + tanhq10(2) * (tanh(tanhq10(3) * (tanhq10(4) - tempq10l)))

          !! Apply a step reduction in q10func when the soil layer freezes. This is to reflect the
          !! lost mobility of the population due to less liquid water
          if (tbar(i,k) - 273.16 > tcrit) then
              q10func = litrq10**(0.1 * (tbar(i,k) - 273.16 - 15.0))
          else
              q10func = litrq10**(0.1 * (tbar(i,k) - 273.16 - 15.0)) * frozered
          end if

          !! Reduce the respiration at depth due to unresolved depth dependent processes including
          !! soil microbial population dynamics, pore-scale oxygen availability, mineral sorption,
          !! priming effects, and other unrepresented processes. This is following Lawrence et al.
          !! Enviro Res Lett 2015. We apply this for all soils.
          reduceatdepth = exp(-delzw(i,k) / r_depthredu)

          ! 2.64 converts bsratelt from kg c/kg c.year to u-mol co2/kg c.s
          ltresveg(i,j,k) = ltrmoscl(i,k) * litrmass(i,j,k) * bsratelt(sort(j)) * 2.64 &
                                          * q10func * reduceatdepth

          ! Respiration from soil c pool
          tempq10s = tbar(i,k) - 273.16
          soilcq10 = tanhq10(1) + tanhq10(2) *(tanh(tanhq10(3) * (tanhq10(4) - tempq10s)))

          if (tbar(i,k) - 273.16 > tcrit) then
            q10func = soilcq10**(0.1 * (tbar(i,k) - 273.16 - 15.0))
          else
            q10func = soilcq10**(0.1 * (tbar(i,k) - 273.16 - 15.0)) * frozered
          end if
          ! 2.64 converts bsratesc from kg c/kg c.year to u-mol co2/kg c.s
          scresveg(i,j,k) = socmoscl(i,k) * soilcmas(i,j,k) * bsratesc(sort(j)) * 2.64 &
                                          * q10func * reduceatdepth
        endif
340   continue
330 continue
320 continue

  return

  end subroutine hetresv
!!@}

!>\defgroup hetresg Heterotrophic Respiration Bare Ground
!!Heterotrophic Respiration Subroutine For Bare Fraction
!!
!>\defgroup hetresv Heterotrophic Respiration Vegetated
!!Heterotrophic Respiration Subroutine For Vegetated Fraction
!!
!>\file
!!Heterotrophic Respiration Module (Vegetated and Bare Ground)
!!
!!Central module for all heterotrophic respiration-related operations
!!
!!Heterotrophic respiration, \f$R_\mathrm{h}\f$ (\f$mol\,CO_2\,m^{-2}\,s^{-1}\f$), in CTEM is
!!based on respiration from the litter (which includes contributions from the stem, leaf
!!and root components), \f$R_{h,D}\f$, and soil carbon, \f$R_{h,H}\f$, pools,
!!
!!\f[ \label{hetres_all} R_\mathrm{h}=R_{h,D}+R_{h,H}. \hspace{10pt}[Eqn 1] \f]
!!
!!Heterotrophic respiration is regulated by soil temperature and moisture and is
!!calculated on a daily time step. The original heterotrophic respiration scheme is
!!described in Arora (2003) \cite Arora2003-3b7 while the modified parametrization used in CTEM v. 2.0
!!is detailed in Melton and Arora (2014) \cite Melton2014-xy.
!!Respiration from the litter and soil carbon pools takes the following basic form
!!
!!\f[ R_{\mathrm{h},i} = 2.64 \times 10^{-6}\,\varsigma_i C_i f_{15}(Q_{10}) f(\Psi)_i f(z),
!!\nonumber \\ i = \mathrm{D}, \mathrm{H}. \hspace{10pt}[Eqn 2] \f]
!!
!!The soil carbon and litter respiration depends on the amount of carbon in these components
!!(\f$C_\mathrm{H}\f$ and \f$C_\mathrm{D}\f$; \f$kg\,C\,m^{-2}\f$) and on a PFT-dependent
!!respiration rate specified at \f$15\,{C}\f$ (\f$\varsigma_\mathrm{H}\f$ and
!!\f$\varsigma_\mathrm{D}\f$; \f$kg\,C\,(kg\,C)^{-1}\,yr^{-1}\f$; see also
!! classic_params.f90). The constant \f$2.64 \times 10^{-6}\f$ converts units from
!!\f$kg\,C\,m^{-2}\,yr^{-1}\f$ to \f$mol\,CO_2\,m^{-2}\,s^{-1}\f$.
!!
!!The effect of soil moisture is accounted for via dependence on soil matric
!!potential (\f$f(\Psi)\f$), described later. The temperature dependency of
!!microbial soil respiration rates has been estimated by several different
!! formulations, ranging from simple \f$Q_{10}\f$ (exponential) to Arrhenius-type
!!formulations (see review by Lloyd and Taylor (1994) \cite Lloyd1994-ct). In CTEM, soil temperature
!!influences heterotrophic respiration through a temperature-dependent
!!\f$Q_{10}\f$ function (\f$f_{15}(Q_{10})\f$). The value of \f$Q_{10}\f$
!!itself is assumed to be a function of temperature following a hyperbolic
!!tan function:
!!
!!\f[ Q_{10} = 1.44 + 0.56\,\tanh[0.075 (46.0 - T_i)], \nonumber\\ i
!!= \mathrm{D}, \mathrm{H}, \hspace{10pt}[Eqn 3]\f]
!!
!!where \f$T_{\{D,H\}}\f$ is the temperature of either the litter or soil
!!carbon pool (\f$C\f$), respectively. The parametrization is a compromise
!!between the temperature-independent \f$Q_{10}\f$ commonly found in many
!!terrestrial ecosystem models (Cox, 2011) \cite Cox2001-am and the temperature-dependent
!!\f$Q_{10}\f$ of Kirschbaum (1995) \cite Kirschbaum1995-db. While a constant \f$Q_{10}\f$ yields
!!an indefinitely increasing respiration rate with increasing temperature, the
!!formulation of Kirschbaum (1995) \cite Kirschbaum1995-db gives a continuously increasing
!!\f$Q_{10}\f$ under decreasing temperature, which leads to unreasonably high
!!soil and litter carbon pools at high latitudes in CTEM. The CTEM
!!parametrization avoids these issues with a \f$Q_{10}\f$ value of about 2.0
!!for temperatures less than \f$20\,C\f$, while a decreasing value of
!!\f$Q_{10}\f$ at temperatures above \f$20\,C\f$ ensures that the
!!respiration rate does not increase indefinitely. As the soil temperature decreases below
!! \f$ T_{crit} \f$ (typically \f$1\,C\f$) a step function is applied to the \f$f_{15}(Q_{10})\f$
!! function to reflect lost mobility of the microbial populations due to less liquid water as:
!!\f[ f_{15}(Q_{10,(T_i < T_{crit})}) = f_{15}(Q_{10}) * 0.1 \f]
!!
!!   \image html Q10_response_sm.png
!!   \image latex Q10_response.eps
!!
!!The soil detrital pools are explictly tracked per soil layer.
!!
!!The response of heterotrophic respiration to soil moisture is formulated through
!!soil matric potential (\f$\Psi\f$; \f$MPa\f$). While soil matric potential values
!!are usually negative, the formulation uses absolute values to allow its logarithm
!!to be taken. Absolute values of soil matric potential are high when soil is dry
!!and low when it is wet. The primary premise of soil moisture control on heterotrophic
!!respiration is that heterotrophic respiration is constrained both when the soils
!!are dry (due to reduced microbial activity) and when they are wet (due to impeded
!!oxygen supply to microbes) with optimum conditions in-between. The exception is the
!!respiration from the litter component of the first soil layer, which is assumed to be continually exposed
!!to air, and thus never oxygen deprived, even when soil moisture content is high
!!(\f$0.04 > \vert \Psi \vert \geq \vert \Psi_{sat} \vert\f$, where \f$\Psi_{sat}\f$
!!is the soil matric potential at saturation). The soil moisture dependence for each
!!soil layer thus varies between 0 and 1 with matric potential as follows:
!!
!!for \f$0.04 > \vert\Psi_i\vert \geq \vert\Psi_{sat,i}\vert\f$
!!
!!\f[ f(\Psi_i)_\mathrm{H,(D, i>1)} = 1 - 0.5  \frac{\log(0.04)-\log\vert\Psi_i\vert}
!!{\log(0.04)-\log\vert\Psi_{sat,i}\vert} \hspace{10pt}[Eqn 4]\f]
!!
!!\f[f(\Psi_i)_{D,i=1} = 1\nonumber; \hspace{10pt}[Eqn 5]\f]
!!
!!for \f$0.06 \geq \vert\Psi_i\vert \geq 0.04\f$
!!\f[ f(\Psi_i)_{D,H} = 1; \hspace{10pt}[Eqn 6]\f]
!!
!!for \f$100.0 \geq \vert\Psi_i\vert > 0.06\f$
!!\f[ f(\Psi_i)_{D,H} = 1 - 0.8\frac{\log\vert\Psi_i\vert-\log(0.06)}{\log(100)-\log(0.06)}; \hspace{10pt}[Eqn 7]\f]
!!
!!for \f$\vert\Psi_i\vert > 100.0\f$
!!\f[ f(\Psi_i)_{D,H}=0.2. \hspace{10pt}[Eqn 8]\f]
!!
!!Respiration also is reduced at depth in soil (\f$f(z)\f$) following Lawrence et al. (2015)
!! \cite Lawrence2015-tj. This term is meant to represent unresolved depth dependent soil
!!processes (such as oxygen availability, microbial community changes, etc.). The reduction
!! in respiration with depth per soil layer is dependent upon the layer depth and
!!a term, \f$z_t\f$, which is given a value of 10.0 (see ctem_params.f90) as,
!!
!!\f[ f(z_i) =\exp (-z_i / z_t) \hspace{10pt}[Eqn 9]\f]
!!
!!   \image html decr_resp_wit_depth_sm.png
!!   \image latex decr_resp_wit_depth.eps
!!
!!Heterotrophic respiration for bare ground is treated separately in CTEM. The carbon
!!contributions to the bare ground litter and soil carbon pools come via processes
!!such as creation of bare ground due to fire, competition between PFTs and land use
!!change. The heterotrophic respiration is sensitive to temperature and moisture in
!!the same manner as vegetated areas using Eqs. (2)--(8). The
!!base respiration rates of \f$\varsigma_{D,bare}\f$ and \f$\varsigma_{H,bare}\f$ are
!!set to \f$0.5605\f$ and \f$0.02258\,kg\,C\,(kg\,C)^{-1}\,yr^{-1}\f$, respectively.
!!
!!The amount of humidified litter, which is transferred from the litter to the soil
!!carbon pool (\f$C_{\mathrm{D} \rightarrow \mathrm{H}}\f$) is modelled as a fraction
!!of litter respiration (\f$R_{h,D}\f$) as
!!
!!\f[ \label{cdtoh} C_{\mathrm{D} \rightarrow \mathrm{H}} = \chi\,R_{h,D} \hspace{10pt}[Eqn 10] \f]
!!
!!where \f$\chi\f$ (see also ctem_params.f90) is the PFT-dependent humification factor
!!and varies between 0.4 and 0.5. For crops, \f$\chi\f$ is set to 0.1 to account for
!!reduced transfer of humidified litter to the soil carbon pool which leads to loss in
!!soil carbon when natural vegetation is converted to croplands. Over the bare ground
!!fraction \f$\chi\f$ is set to 0.45.
!!
!!With heterotrophic respiration known, net ecosystem productivity (\f$NEP\f$) is
!!calculated as
!!\f[ NEP = G_{canopy} - R_\mathrm{m} - R_\mathrm{g} - R_\mathrm{h}. \hspace{10pt}[Eqn 11]\f]
!!

!>\file
end module heterotrophic_respiration
