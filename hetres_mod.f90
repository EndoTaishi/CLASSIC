module heterotrophic_respiration

! Central module for all heterotrophic respiration-related operations

! J. Melton. Jan 14 2016

implicit none

! Subroutines contained in this module:
public  :: hetresg
public  :: hetresv

contains

! ------------------------------------------------------------------

subroutine hetresg (litrmass, soilcmas, delzw,  thpor, &
                    il1,      il2,     tbar,  psisat, b, &
                    thliq,    zbotw,   thiceg, &
                        frac,    isnow,      isand, &
!    -------------- inputs above this line, outputs below -------------
                      litres,   socres)

!               Canadian Terrestrial Ecosystem Model (CTEM)
!           Heterotrophic Respiration Subroutine For Bare Fraction
!
!     11  Apr. 2003 - this subroutine calculates heterotrophic respiration
!     V. Arora        over the bare subarea of a grid cell (i.e. ground only
!                     and snow over ground subareas).
!
!     change history:
!
!     14  Jan 2016  - Converted to F90 and made it so it can handle > 3 soil layers
!     J. Melton
!
!     30  Jul 2015  - Based on work by Yuanqiao Wu, respiration was found to
!     J. Melton       behave incorrectly if the soil froze as it thought the water
!                     was leaving the soil. This is now fixed.
!
!     17  Jan 2014  - Moved parameters to global file (ctem_params.f90)
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

!     ------
!     inputs

!     litrmass  - litter mass for the 8 pfts + bare in kg c/m2
!     soilcmas  - soil carbon mass for the 8 pfts + bare in kg c/m2
!     icc       - no. of vegetation types (currently 8)
!     ignd        - no. of soil layers (currently 3)
!     ilg       - no. of grid cells in latitude circle
!     il1,il2   - il1=1, il2=ilg
!     tbar      - soil temperature, k
!     thliq     - liquid soil moisture content in 3 soil layers
!     sand      - percentage sand
!     zbotw     - bottom of soil layers
!     frac      - fraction of ground (fg) or snow over ground (fgs)
!     isnow     - integer telling if bare fraction is fg (0) or fgs (1)
!     delzw     - permeable thickness of the layer (m)
!     thiceg    - frozen soil moisture content over bare ground fraction

!     outputs

!     litres    - litter respiration over the given unvegetated sub-area
!                 in umol co2/m2.s
!     socres    - soil c respiration over the given unvegetated sub-area
!                 in umol co2/m2.s

use ctem_params,        only : icc, ilg, ignd, zero, tanhq10, a, &
                                     bsratelt_g, bsratesc_g

implicit none

integer il1,il2,i,j,k,isnow,isand(ilg,ignd)

real                  litrmass(ilg,icc+1),  soilcmas(ilg,icc+1),&
           tbar(ilg,ignd),    thliq(ilg,ignd),   psisat(ilg,ignd),  &
           zbotw(ilg,ignd),      litres(ilg),          socres(ilg), &
           frac(ilg),       delzw(ilg,ignd), thpor(ilg,ignd),  b(ilg,ignd), &
           thiceg(ilg,ignd)

real              litrq10,           soilcq10, &
        litrtemp(ilg),      solctemp(ilg),             q10func, &
           grksat(ilg,ignd),       &
              beta, &
      fracarb(ilg,ignd),         zcarbon, &
        tempq10l(ilg),      socmoscl(ilg),     scmotrm(ilg,ignd), &
        ltrmoscl(ilg),        psi(ilg,ignd),       tempq10s(ilg), &
               fcoeff,          zcarb_g

!     ------------------------------------------------------------------
!     Constants and parameters are located in ctem_params.f90
!     ---------------------------------------------------------------

!     initialize required arrays to zero

      do 100 k = 1, ignd
        do 100 i = il1, il2
          fracarb(i,k)=0.0  ! fraction of carbon in each soil layer
100   continue

      do 110 i = il1, il2
        litrtemp(i)=0.0     ! litter temperature
        solctemp(i)=0.0     ! soil carbon pool temperature
        socmoscl(i)=0.0     ! soil moisture scalar for soil carbon decomposition
        ltrmoscl(i)=0.0     ! soil moisture scalar for litter decomposition
        litres(i)=0.0       ! litter resp. rate
        tempq10l(i)=0.0
        socres(i)=0.0       ! soil c resp. rate
        tempq10s(i)=0.0
110   continue

!     initialization ends

!     ------------------------------------------------------------------

!     Estimate temperature of the litter and soil carbon pools.

!     Over the bare fraction there is no live root. So we make the
!     simplest assumption that litter temperature is same as temperature of the top soil layer.

      do 210 i = il1, il2
        litrtemp(i)=tbar(i,1)
210   continue

!     We estimate the temperature of the soil c pool assuming that soil carbon over the bare fraction is distributed exponentially. note
!     that bare fraction may contain dead roots from different pfts all of which may be distributed differently. For simplicity we do not
!     track each pft's dead root biomass and assume that distribution of soil carbon over the bare fraction can be described by a single
!     parameter.

      do 240 i = il1, il2

!         zcarbon=3.0/a                 ! 95% depth
!
!         if(zcarbon.le.zbotw(i,1)) then
!             fracarb(i,1)=1.0             ! fraction of carbon in soil layers
!         else
!             fcoeff=exp(-a*zcarbon)
!             fracarb(i,1)=1.0-(exp(-a*zbotw(i,1))-fcoeff)/(1.0-fcoeff)
!             if(zcarbon.le.zbotw(i,2)) then
!                 fracarb(i,2)=1.0-fracarb(i,1)
!                 fracarb(i,3)=0.0
!             else
!                 fracarb(i,3)=(exp(-a*zbotw(i,2))-fcoeff)/(1.0-fcoeff)
!                 fracarb(i,2)=1.0-fracarb(i,1)-fracarb(i,3)
!             endif
!         endif
! New, allows for many soil layers:
        zcarbon = 3.0 / a
        zcarb_g = 0.0
        do j=1,ignd
          zcarb_g = zcarb_g + delzw(i,j)
        end do
        zcarbon = min(zcarbon,zcarb_g)
        fcoeff=exp(-a*zcarbon)
        fracarb(i,1)=1.0-(exp(-a*zbotw(i,1))-fcoeff)/(1.0-fcoeff)
        do 245 j = 2,ignd
          if (zcarbon <= zbotw(i,j) - delzw(i,j)+0.0001) then
             fracarb(i,j) = 0.0
          elseif (zcarbon <= zbotw(i,j)) then
             fracarb(i,j) = (exp(-a * zbotw(i,j)-delzw(i,j))-exp(-a*zcarbon))/(1. - exp(-a*zcarbon))
          else
             fracarb(i,j) = (exp(-a * zbotw(i,j)-delzw(i,j))-exp(-a*zbotw(i,j)))/(1. - exp(-a*zcarbon))
          end if
245       continue

        ! -------------
        !solctemp(i)=tbar(i,1)*fracarb(i,1) +tbar(i,2)*fracarb(i,2) +tbar(i,3)*fracarb(i,3)
        !solctemp(i)=solctemp(i) /(fracarb(i,1)+fracarb(i,2)+fracarb(i,3))
        solctemp(i) = sum(tbar(i,:)*fracarb(i,:))/ sum(fracarb(i,:))


!       make sure we don't use temperatures of 2nd and 3rd soil layers
!       if they are specified bedrock via sand -3 flag
!        JM NOTE: This is not needed. zbotw already accounts for bedrock.
!         if(isand(i,3).eq.-3)then ! third layer bed rock
!           solctemp(i)=tbar(i,1)*fracarb(i,1) +tbar(i,2)*fracarb(i,2)
!           solctemp(i)=solctemp(i) /(fracarb(i,1)+fracarb(i,2))
!         endif
!         if(isand(i,2).eq.-3)then ! second layer bed rock
!           solctemp(i)=tbar(i,1)
!         endif

240   continue

!     find moisture scalar for soil c decomposition

!     this is modelled as function of logarithm of matric potential.
!     we find values for all soil layers, and then find an average value
!     based on fraction of carbon present in each layer.

      do 260 j = 1, ignd
        do 270 i = il1, il2

          if(isand(i,j).eq.-3.or.isand(i,j).eq.-4)then
            scmotrm (i,j)=0.2
            psi (i,j) = 10000.0 ! set to large number so that
                               ! ltrmoscl becomes 0.2
          else ! i.e., sand.ne.-3 or -4

!            psisat(i,j)= (10.0**(-0.0131*sand(i,j)+1.88))/100.0
!            b(i,j)     = 0.159*clay(i,j)+2.91
!            thpor(i,j) = (-0.126*sand(i,j)+48.9)/100.0
!            psi(i,j)   = psisat(i,j)*(thliq(i,j)/thpor(i,j))**(-b(i,j))

            psi(i,j)   = psisat(i,j)*(thliq(i,j)/(thpor(i,j)+0.005 -thiceg(i,j)))**(-b(i,j))

            if(psi(i,j).ge.10000.0) then
              scmotrm(i,j)=0.2
            else if( psi(i,j).lt.10000.0 .and.  psi(i,j).gt.6.0 ) then
              scmotrm(i,j)=1.0 - 0.8*( (log10(psi(i,j)) - log10(6.0))/(log10(10000.0)-log10(6.0)) )
            else if( psi(i,j).le.6.0 .and.  psi(i,j).ge.4.0 ) then
              scmotrm(i,j)=1.0
            else if( psi(i,j).lt.4.0.and.psi(i,j).gt.psisat(i,j) )then
              scmotrm(i,j)=1.0 -0.5*( (log10(4.0) - log10(psi(i,j))) /(log10(4.0)-log10(psisat(i,j))) )
            else if( psi(i,j).le.psisat(i,j) ) then
              scmotrm(i,j)=0.5
            endif
          endif ! if sand.eq.-3 or -4

          scmotrm(i,j)=max(0.0,min(scmotrm(i,j),1.0))
270     continue
260   continue

      do 290 i = il1, il2

        socmoscl(i) = sum(scmotrm(i,:)*fracarb(i,:)) / sum(fracarb(i,:))
!         socmoscl(i) = scmotrm(i,1)*fracarb(i,1) +scmotrm(i,2)*fracarb(i,2) +scmotrm(i,3)*fracarb(i,3)
!         socmoscl(i) = socmoscl(i) /(fracarb(i,1)+fracarb(i,2)+fracarb(i,3))

!       make sure we don't use scmotrm of 2nd and 3rd soil layers
!       if they are specified bedrock via sand -3 flag

!         if(isand(i,3).eq.-3)then ! third layer bed rock
!           socmoscl(i) = scmotrm(i,1)*fracarb(i,1) +scmotrm(i,2)*fracarb(i,2)
!           socmoscl(i) = socmoscl(i) /(fracarb(i,1)+fracarb(i,2))
!         endif
!         if(isand(i,2).eq.-3)then ! second layer bed rock
!           socmoscl(i) = scmotrm(i,1)
!         endif

        socmoscl(i)=max(0.2,min(socmoscl(i),1.0))

290   continue

!     find moisture scalar for litter decomposition

!     the difference between moisture scalar for litter and soil c
!     is that the litter decomposition is not constrained by high
!     soil moisture (assuming that litter is always exposed to air).
!     in addition, we use moisture content of the top soil layer
!     as a surrogate for litter moisture content. so we use only
!     psi(i,1) calculated in loops 260 and 270 above.

      do 300 i = il1, il2
        if(psi(i,1).gt.10000.0) then
          ltrmoscl(i)=0.2
        else if( psi(i,1).le.10000.0 .and.  psi(i,1).gt.6.0 ) then
          ltrmoscl(i)=1.0 - 0.8*( (log10(psi(i,1)) - log10(6.0))/(log10(10000.0)-log10(6.0)) )
        else if( psi(i,1).le.6.0 ) then
          ltrmoscl(i)=1.0
        endif
        ltrmoscl(i)=max(0.2,min(ltrmoscl(i),1.0))
300   continue

!     use temperature of the litter and soil c pools, and their soil
!     moisture scalars to find respiration rates from these pools

      do 330 i = il1, il2
      if(frac(i).gt.zero)then

!       first find the q10 response function to scale base respiration
!       rate from 15 c to current temperature, we do litter first

        tempq10l(i)=litrtemp(i)-273.16
        litrq10 = tanhq10(1) + tanhq10(2)*( tanh( tanhq10(3)*(tanhq10(4)-tempq10l(i))  ) )

        q10func = litrq10**(0.1*(litrtemp(i)-273.16-15.0))
        litres(i)= ltrmoscl(i) * litrmass(i,icc+1)*bsratelt_g*2.64*q10func ! 2.64 converts bsratelt_g from kg c/kg c.year
                                                                         ! to u-mol co2/kg c.s

!       respiration from soil c pool

        tempq10s(i)=solctemp(i)-273.16
        soilcq10= tanhq10(1) + tanhq10(2)*( tanh( tanhq10(3)*(tanhq10(4)-tempq10s(i))  ) )

        q10func = soilcq10**(0.1*(solctemp(i)-273.16-15.0))
        socres(i)= socmoscl(i)* soilcmas(i,icc+1)*bsratesc_g*2.64*q10func ! 2.64 converts bsratesc_g from kg c/kg c.year
                                                                          ! to u-mol co2/kg c.s

      endif
330   continue

      return
end subroutine hetresg

! ------------------------------------------------------------------------------------

subroutine hetresv ( fcan,      fct, litrmass, soilcmas, &
                      delzw,  thpor, il1, &
                      il2,     tbar,   psisat, b, thliq,  &
                     roottemp,    zbotw,     sort, &
                     isand, thicec, &
!    -------------- inputs above this line, outputs below -------------
                    ltresveg, scresveg)

!               Canadian Terrestrial Ecosystem Model (CTEM)
!           Heterotrophic Respiration Subtoutine For Vegetated Fraction

!     16  oct. 2001 - this subroutine calculates heterotrophic respiration
!     V. Arora        for a given sub-area, from litter and soil carbon
!                     pools.

!     change history:

!     30  Jul 2015  - Based on work by Yuanqiao Wu, respiration was found to
!     J. Melton       behave incorrectly if the soil froze as it thought the water
!                     was leaving the soil. This is now fixed.
!     17  Jan 2014  - Moved parameters to global file (ctem_params.f90)
!     J. Melton

!     22  Jul 2013  - Add in module for parameters
!     J. Melton

!     J. Melton and V.Arora - changed tanhq10 parameters, they were switched
!               25 sep 2012
!     J. Melton 23 aug 2012 - bring in isand, converting sand to
!                             int was missing some gridcells assigned
!                             to bedrock in classb
!     ------
!     inputs

!     fcan      - fractional coverage of ctem's 9 pfts
!     fct       - sum of all fcan
!                 fcan & fct are not used at this time but could
!                 be used at some later stage
!     litrmass  - litter mass for the 9 pfts + bare in kg c/m2
!     soilcmas  - soil carbon mass for the 9 pfts + bare in kg c/m2
!     icc       - no. of vegetation types (currently 9)
!     ignd        - no. of soil layers (currently 3)
!     ilg       - no. of grid cells in latitude circle
!     il1,il2   - il1=1, il2=ilg
!     tbar      - soil temperature, k
!     thliq     - liquid soil moisture content in 3 soil layers
!     thicec    - frozen soil moisture content in 3 soil layers in canopy
!                 covered subarea
!     sand      - percentage sand
!     clay      - percentage clay
!     roottemp  - root temperature as estimated in mainres subroutine
!     zbotw     - bottom of soil layers
!     sort      - index for correspondence between 9 pfts and 12 values
!                 in the parameters vectors
!     delzw     - permeable thickness of the layer (m)

!     outputs

!     ltresveg  - litter respiration for the given sub-area in
!                 umol co2/m2.s, for ctem's 9 pfts
!     scresveg  - soil carbon respiration for the given sub-area in
!                 umol co2/m2.s, for ctem's 9 pfts

      use ctem_params,        only : icc, ilg, ignd, kk, zero, bsratelt,&
                               bsratesc, abar, tanhq10,&
                               alpha_hetres

      implicit none

      integer  il1, il2, i, j, k, sort(icc), isand(ilg,ignd)

      real    fcan(ilg,icc),           fct(ilg),  litrmass(ilg,icc+1),&
         tbar(ilg,ignd),soilcmas(ilg,icc+1),      thliq(ilg,ignd),&
        roottemp(ilg,icc), psisat(ilg,ignd),  b(ilg,ignd), &
        zbotw(ilg,ignd),  ltresveg(ilg,icc),    scresveg(ilg,icc),&
        thicec(ilg,ignd),  delzw(ilg,ignd),  thpor(ilg,ignd)


      real           litrq10,           soilcq10,&
        litrtemp(ilg,icc),  solctemp(ilg,icc),             q10func,&
       grksat(ilg,ignd),       &
       fracarb(ilg,icc,ignd),                       zcarbon,&
         tempq10l(ilg,icc),  socmoscl(ilg,icc),     scmotrm(ilg,ignd),&
             ltrmoscl(ilg),        psi(ilg,ignd),   tempq10s(ilg,icc),&
                    fcoeff, zcarb_g

!     ------------------------------------------------------------------
!     Constants and parameters are located in ctem_params.f90
!     ---------------------------------------------------------------

!     initialize required arrays to zero

do 100 j = 1, icc
  do 110 i = il1, il2
          litrtemp(i,j)=0.0       ! litter temperature
          tempq10l(i,j)=0.0
          solctemp(i,j)=0.0       ! soil carbon pool temperature
          tempq10s(i,j)=0.0
          socmoscl(i,j)=0.0       ! soil moisture scalar for soil carbon decomposition
          ltresveg(i,j)=0.0       ! litter resp. rate for each pft
          scresveg(i,j)=0.0       ! soil c resp. rate for each pft
          do 120 k = 1, ignd
            scmotrm(i,k)=0.0        ! soil carbon moisture term
            fracarb(i,j,k)=0.0    ! fraction of carbon in each soil layer for each vegetation
          120   continue
110     continue
100   continue

do 130 i = il1, il2
    ltrmoscl(i)=0.0           ! soil moisture scalar for litter decomposition
130   continue

!     initialization ends

!     ------------------------------------------------------------------

!     estimate temperature of the litter and soil carbon pools. litter
!     temperature is weighted average of temperatue of top soil layer
!     (where the stem and leaf litter sits) and root temperature, because
!     litter pool is made of leaf, stem, and root litter.

      do 200 j = 1,icc
        do 210 i = il1, il2
         if (fcan(i,j) .gt. 0.) then
          litrtemp(i,j)=alpha_hetres*tbar(i,1)+roottemp(i,j)*(1.0-alpha_hetres)
         endif
210     continue
200   continue

!     estimation of soil carbon pool temperature is not straight forward.
!     ideally soil c pool temperature should be set same as root temperature,
!     since soil c profiles are similar to root profiles. but in the event
!     when the roots die then we may run into trouble. so we find the
!     temperature of the soil c pool assuming that soil carbon is
!     exponentially distributed, just like roots. but rather than using
!     the parameter of this exponential profile from our variable root
!     distribution we use fixed vegetation-dependent parameters.

      do 230 j = 1, icc
        do 240 i = il1, il2
         if (fcan(i,j) .gt. 0.) then

!           zcarbon=3.0/abar(sort(j))                ! 95% depth
!           if(zcarbon.le.zbotw(i,1)) then
!               fracarb(i,j,1)=1.0             ! fraction of carbon in
!               fracarb(i,j,2)=0.0             ! soil layers
!               fracarb(i,j,3)=0.0
!           else
!               fcoeff=exp(-abar(sort(j))*zcarbon)
!               fracarb(i,j,1)=1.0-(exp(-abar(sort(j))*zbotw(i,1))-fcoeff)/(1.0-fcoeff)
!               if(zcarbon.le.zbotw(i,2)) then
!                   fracarb(i,j,2)=1.0-fracarb(i,j,1)
!                   fracarb(i,j,3)=0.0
!               else
!                   fracarb(i,j,3)=(exp(-abar(sort(j))*zbotw(i,2))-fcoeff)/(1.0-fcoeff)
!                   fracarb(i,j,2)=1.0-fracarb(i,j,1)-fracarb(i,j,3)
!               endif
!           endif
        ! New, allows for many soil layers:
        zcarbon = 3.0 / abar(sort(j))
        zcarb_g = 0.0
        do k=1,ignd
          zcarb_g = zcarb_g + delzw(i,k)
        end do
        zcarbon = min(zcarbon,zcarb_g)
        fcoeff=exp(-abar(sort(j))*zcarbon)
        fracarb(i,j,1)=1.0-(exp(-abar(sort(j))*zbotw(i,1))-fcoeff)/(1.0-fcoeff)
        do 245 k = 2,ignd
          if (zcarbon <= zbotw(i,k) - delzw(i,k)+0.0001) then
             fracarb(i,j,k) = 0.0
          elseif (zcarbon <= zbotw(i,k)) then
             fracarb(i,j,k) = (exp(-abar(sort(j)) * zbotw(i,k)-delzw(i,k))- &
                               exp(-abar(sort(j))*zcarbon))/(1. - exp(-abar(sort(j))*zcarbon))
          else
             fracarb(i,j,k) = (exp(-abar(sort(j)) * zbotw(i,k)-delzw(i,k))- &
                               exp(-abar(sort(j))*zbotw(i,k)))/(1. - exp(-abar(sort(j))*zcarbon))
          end if
245       continue

!           solctemp(i,j)=tbar(i,1)*fracarb(i,j,1) +tbar(i,2)*fracarb(i,j,2) +tbar(i,3)*fracarb(i,j,3)
!           solctemp(i,j)=solctemp(i,j) /(fracarb(i,j,1)+fracarb(i,j,2)+fracarb(i,j,3))
          solctemp(i,j) = sum(tbar(i,:)*fracarb(i,j,:))/ sum(fracarb(i,j,:))

!         make sure we don't use temperatures of 2nd and 3rd soil layers
!         if they are specified bedrock via sand -3 flag

!           if(isand(i,3).eq.-3)then ! third layer bed rock
!             solctemp(i,j)=tbar(i,1)*fracarb(i,j,1) +tbar(i,2)*fracarb(i,j,2)
!             solctemp(i,j)=solctemp(i,j) /(fracarb(i,j,1)+fracarb(i,j,2))
!           endif
!           if(isand(i,2).eq.-3)then ! second layer bed rock
!             solctemp(i,j)=tbar(i,1)
!           endif
        endif
240     continue
230   continue

!     find moisture scalar for soil c decomposition

!     this is modelled as function of logarithm of matric potential.
!     we find values for all soil layers, and then find an average value
!     based on fraction of carbon present in each layer. this makes
!     moisture scalar a function of vegetation type.

      do 260 j = 1, ignd
        do 270 i = il1, il2

          if(isand(i,j).eq.-3.or.isand(i,j).eq.-4)then
            scmotrm (i,j)=0.2
            psi (i,j) = 10000.0 ! set to large number so that
                               ! ltrmoscl becomes 0.2
          else ! i.e., sand.ne.-3 or -4

!             psisat(i,j)= (10.0**(-0.0131*sand(i,j)+1.88))/100.0
!             b(i,j)     = 0.159*clay(i,j)+2.91
!             thpor(i,j) = (-0.126*sand(i,j)+48.9)/100.0

            !the 0.005 below prevents a divide by 0 situation.
            psi(i,j)   = psisat(i,j)*(thliq(i,j)/(thpor(i,j)+0.005 -thicec(i,j)))**(-b(i,j))

            if(psi(i,j).ge.10000.0) then
              scmotrm(i,j)=0.2
            else if( psi(i,j).lt.10000.0 .and.  psi(i,j).gt.6.0 ) then
              scmotrm(i,j)=1.0 - 0.8*( (log10(psi(i,j)) - log10(6.0))/(log10(10000.0)-log10(6.0)) )
            else if( psi(i,j).le.6.0 .and. psi(i,j).ge.4.0 ) then
              scmotrm(i,j)=1.0
            else if( psi(i,j).lt.4.0 .and. psi(i,j).gt.psisat(i,j) )then
              scmotrm(i,j)=1.0 -0.5*( (log10(4.0) - log10(psi(i,j))) /(log10(4.0)-log10(psisat(i,j))) )
            else if( psi(i,j).le.psisat(i,j) ) then
              scmotrm(i,j)=0.5
            endif
            scmotrm(i,j)=max(0.2,min(1.0,scmotrm(i,j)))
          endif ! sand.eq.-3 or -4

270     continue
260   continue

      do 280 j = 1, icc
        do 290 i = il1, il2
         if (fcan(i,j) .gt. 0.) then

          socmoscl(i,j) = sum(scmotrm(i,:)*fracarb(i,j,:)) / sum(fracarb(i,j,:))
!           socmoscl(i,j) = scmotrm(i,1)*fracarb(i,j,1) +scmotrm(i,2)*fracarb(i,j,2) +scmotrm(i,3)*fracarb(i,j,3)
!           socmoscl(i,j) = socmoscl(i,j) /(fracarb(i,j,1)+fracarb(i,j,2)+fracarb(i,j,3))

!         make sure we don't use scmotrm of 2nd and 3rd soil layers
!         if they are specified bedrock via sand -3 flag

!           if(isand(i,3).eq.-3)then ! third layer bed rock
!             socmoscl(i,j) = scmotrm(i,1)*fracarb(i,j,1) +scmotrm(i,2)*fracarb(i,j,2)
!             socmoscl(i,j) = socmoscl(i,j) /(fracarb(i,j,1)+fracarb(i,j,2))
!           endif
!           if(isand(i,2).eq.-3)then ! second layer bed rock
!             socmoscl(i,j) = scmotrm(i,1)
!           endif

          socmoscl(i,j)=max(0.2,min(1.0,socmoscl(i,j)))

         endif
290     continue
280   continue

!     find moisture scalar for litter decomposition

!     the difference between moisture scalar for litter and soil c
!     is that the litter decomposition is not constrained by high
!     soil moisture (assuming that litter is always exposed to air).
!     in addition, we use moisture content of the top soil layer
!     as a surrogate for litter moisture content. so we use only
!     psi(i,1) calculated in loops 260 and 270 above.

      do 300 i = il1, il2
        if(psi(i,1).gt.10000.0) then
          ltrmoscl(i)=0.2
        else if( psi(i,1).le.10000.0 .and. psi(i,1).gt.6.0 ) then
          ltrmoscl(i)=1.0 -  0.8*( (log10(psi(i,1)) - log10(6.0))/(log10(10000.0)-log10(6.0)) )
        else if( psi(i,1).le.6.0 ) then
          ltrmoscl(i)=1.0
        endif
        ltrmoscl(i)=max(0.2,min(1.0,ltrmoscl(i)))
300   continue

!     use temperature of the litter and soil c pools, and their soil
!     moisture scalars to find respiration rates from these pools

      do 320 j = 1, icc
        do 330 i = il1, il2
         if (fcan(i,j) .gt. 0.) then

!         first find the q10 response function to scale base respiration
!         rate from 15 c to current temperature, we do litter first

          tempq10l(i,j)=litrtemp(i,j)-273.16
          litrq10 = tanhq10(1) + tanhq10(2)*( tanh( tanhq10(3)*(tanhq10(4)-tempq10l(i,j))  ) )

          q10func = litrq10**(0.1*(litrtemp(i,j)-273.16-15.0))
          ltresveg(i,j)= ltrmoscl(i) * litrmass(i,j)*bsratelt(sort(j))*2.64*q10func ! 2.64 converts bsratelt from kg c/kg c.year
                                          ! to u-mol co2/kg c.s
!         respiration from soil c pool

          tempq10s(i,j)=solctemp(i,j)-273.16
          soilcq10= tanhq10(1) + tanhq10(2)*( tanh( tanhq10(3)*(tanhq10(4)-tempq10s(i,j))  ) )

          q10func = soilcq10**(0.1*(solctemp(i,j)-273.16-15.0))
          scresveg(i,j)= socmoscl(i,j)* soilcmas(i,j)*bsratesc(sort(j))*2.64*q10func  ! 2.64 converts bsratesc from kg c/kg c.year
                                                                                     ! to u-mol co2/kg c.s
         endif
330     continue
320   continue

      return

end subroutine hetresv

end module heterotrophic_respiration