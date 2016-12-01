!>\defgroup hetresg
!!
!!Canadian Terrestrial Ecosystem Model (CTEM)
!!Heterotrophic respiration Subroutine
!!
!!Heterotrophic Respiration Subroutine For Bare Fraction
!!
!! FLAG This needs to be edited for the peatland changes!
!!Heterotrophic respiration, \f$R_\mathrm{h}\f$ (\f$mol\,CO_2\,m^{-2}\,s^{-1}\f$), in CTEM is
!!based on respiration from the litter (which includes contributions from the stem, leaf
!!and root components), \f$R_{h,D}\f$, and soil carbon, \f$R_{h,H}\f$, pools,
!!\f[ \label{hetres_all} R_\mathrm{h}=R_{h,D}+R_{h,H}. \f]
!!
!!Heterotrophic respiration is regulated by soil temperature and moisture and is
!!calculated on a daily time step. The original heterotrophic respiration scheme is
!!described in \cite Arora2003-3b7 while the modified parametrization used in CTEM v. 2.0
!!is detailed in \cite Melton2014-xy and is briefly described here. Respiration from the
!!litter and soil carbon pools takes the following basic form
!!\f[ R_{\mathrm{h},i} = 2.64 \times 10^{-6}\,\varsigma_i C_i f_{15}(Q_{10}) f(\Psi)_i,
!!\nonumber \\ i = \mathrm{D}, \mathrm{H}.\label{lithet} \f]
!!
!!The soil carbon and litter respiration depends on the amount of carbon in these components
!!(\f$C_\mathrm{H}\f$ and \f$C_\mathrm{D}\f$; \f$kg\,C\,m^{-2}\f$) and on a PFT-dependent
!!respiration rate specified at \f$15\,{C}\f$ (\f$\varsigma_\mathrm{H}\f$ and
!!\f$\varsigma_\mathrm{D}\f$; \f$kg\,C\,(kg\,C)^{-1}\,yr^{-1}\f$; see also
!! ctem_params.f90). The constant \f$2.64 \times 10^{-6}\f$ converts units from
!!\f$kg\,C\,m^{-2}\,yr^{-1}\f$ to \f$mol\,CO_2\,m^{-2}\,s^{-1}\f$.
!!
!!The effect of soil moisture is accounted for via dependence on soil matric
!!potential (\f$f(\Psi)\f$), described later. The temperature dependency of
!!microbial soil respiration rates has been estimated by several different
!! formulations, ranging from simple \f$Q_{10}\f$ (exponential) to Arrhenius-type
!!formulations (see review by \cite Lloyd1994-ct). In CTEM, soil temperature
!!influences heterotrophic respiration through a temperature-dependent
!!\f$Q_{10}\f$ function (\f$f_{15}(Q_{10})\f$). The value of \f$Q_{10}\f$
!!itself is assumed to be a function of temperature following a hyperbolic
!!tan function:
!!\f[ Q_{10} = 1.44 + 0.56\,\tanh[0.075 (46.0 - T_i)], \nonumber\\ i
!!= \mathrm{D}, \mathrm{H}\label{hyper}, \f]
!!where \f$T_{\{D,H\}}\f$ is the temperature of either the litter or soil
!!carbon pool (\f$C\f$), respectively. The parametrization is a compromise
!!between the temperature-independent \f$Q_{10}\f$ commonly found in many
!!terrestrial ecosystem models \cite Cox2001-am and the temperature-dependent
!!\f$Q_{10}\f$ of \cite Kirschbaum1995-db. While a constant \f$Q_{10}\f$ yields
!!an indefinitely increasing respiration rate with increasing temperature, the
!!formulation of \cite Kirschbaum1995-db gives a continuously increasing
!!\f$Q_{10}\f$ under decreasing temperature, which leads to unreasonably high
!!soil and litter carbon pools at high latitudes in CTEM. The CTEM
!!parametrization avoids these issues with a $Q_{10}$ value of about 2.0
!!for temperatures less than \f$20\,C\f$, while a decreasing value of
!!\f$Q_{10}\f$ at temperatures above \f$20\,C\f$ ensures that the
!!respiration rate does not increase indefinitely.
!!
!!The temperature of the litter pool is a weighted average of the
!!temperature of the top soil layer (\f$T_1\f$) and the root temperature
!!(\f$T_\mathrm{R}\f$) as litter consists of leaf, stem, and root litter
!!(\f$T_\mathrm{D} = 0.7 T_1 + 0.3T_\mathrm{R}\f$). The temperature of the
!!soil carbon pool is calculated as the mean soil temperature in the rooting
!!zone based upon the fraction of roots in each soil layer and their temperature.
!!The carbon in each soil layer is not explicitly tracked but assumed to adopt an
!!exponential distribution \cite Jobbagy2000-pa.
!!
!!The response of heterotrophic respiration to soil moisture is formulated through
!!soil matric potential (\f$\Psi\f$; \f$MPa\f$). While soil matric potential values
!!are usually negative, the formulation uses absolute values to allow its logarithm
!!to be taken. Absolute values of soil matric potential are high when soil is dry
!!and low when it is wet. The primary premise of soil moisture control on heterotrophic
!!respiration is that heterotrophic respiration is constrained both when the soils
!!are dry (due to reduced microbial activity) and when they are wet (due to impeded
!!oxygen supply to microbes) with optimum conditions in-between. The exception is the
!!respiration from the litter component, which is assumed to be continually exposed
!!to air, and thus never oxygen deprived, even when soil moisture content is high
!!(\f$0.04 > \vert \Psi \vert \geq \vert \Psi_{sat} \vert\f$, where \f$\Psi_{sat}\f$
!!is the soil matric potential at saturation). The soil moisture dependence thus
!!varies between 0 and 1 with matric potential as follows:
!!
!!for \f$0.04 > \vert\Psi\vert \geq \vert\Psi_{sat}\vert\f$
!!\f[ f(\Psi)_\mathrm{H} = 1 - 0.5  \frac{\log(0.04)-\log\vert\Psi\vert}
!!{\log(0.04)-\log\vert\Psi_{sat}\vert}\\ f(\Psi)_D = 1\nonumber; \f]
!!for \f$0.06 \geq \vert\Psi\vert \geq 0.04\f$
!!\f[ f(\Psi)_\{D,H\} = 1; \f]
!!for \f$100.0 \geq \vert\Psi\vert > 0.06\f$
!!\f[ f(\Psi)_\{D,H\} = 1 - 0.8\frac{\log\vert\Psi\vert-\log(0.06)}{\log(100)-\log(0.06)}; \f]
!!for \f$\vert\Psi\vert > 100.0\f$
!!\f[ \label{lastpsi} f(\Psi)_\{D,H\}=0.2. \f]
!!
!!Heterotrophic respiration for bare ground is treated separately in CTEM. The carbon
!!contributions to the bare ground litter and soil carbon pools come via processes
!!such as creation of bare ground due to fire, competition between PFTs and land use
!!change. The heterotrophic respiration is sensitive to temperature and moisture in
!!the same manner as vegetated areas using Eqs. (\ref{lithet})--(\ref{lastpsi}). The
!!base respiration rates of \f$\varsigma_{D,bare}\f$ and \f$\varsigma_{H,bare}\f$ are
!!set to \f$0.5605\f$ and \f$0.02258\,kg\,C\,(kg\,C)^{-1}\,yr^{-1}\f$, respectively.
!!
!!The amount of humidified litter, which is transferred from the litter to the soil
!!carbon pool (\f$C_{\mathrm{D} \rightarrow \mathrm{H}}\f$) is modelled as a fraction
!!of litter respiration (\f$R_{h,D}\f$) as
!!\f[ \label{cdtoh} C_{\mathrm{D} \rightarrow \mathrm{H}} = \chi\,R_{h,D} \f]
!!where \f$\chi\f$ (see also ctem_params.f90) is the PFT-dependent humification factor
!!and varies between 0.4 and 0.5. For crops, \f$\chi\f$ is set to 0.1 to account for
!!reduced transfer of humidified litter to the soil carbon pool which leads to loss in
!!soil carbon when natural vegetation is converted to croplands. Over the bare ground
!!fraction \f$\chi\f$ is set to 0.45.
!!
!!With heterotrophic respiration known, net ecosystem productivity (\f$NEP\f$) is
!!calculated as
!!\f[ NEP = G_{canopy} - R_\mathrm{m} - R_\mathrm{g} - R_\mathrm{h}. \f]
!!

!>\defgroup hetresv

!>
!!Heterotrophic Respiration Subroutine For Vegetated Fraction

!>\file
!!Central module for all heterotrophic respiration-related operations
module heterotrophic_respiration

! J. Melton. Jan 14 2016

implicit none

! Subroutines contained in this module:
public  :: hetresg
public  :: hetresv

contains

! ------------------------------------------------------------------
!>\ingroup hetresg
!!@{

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

use ctem_params,        only : icc, ilg, ignd, zero, tanhq10, a, &
                                     bsratelt_g, bsratesc_g

implicit none

      integer il1   !<il1=1
      integer il2   !<il2=ilg
      integer i,j,k
      integer isnow !<integer telling if bare fraction is fg (0) or fgs (1), isnow
                    !<is changed to isnow(ilg) in classt of class version higher than 3.4 for coupling with ctem
      integer isand(ilg,ignd) !<

      real litrmass(ilg,icc+1)!<litter mass for the 8 pfts + bare in \f$kg c/m^2\f$
      real soilcmas(ilg,icc+1)!<soil carbon mass for the 8 pfts + bare in \f$kg c/m^2\f$
      real tbar(ilg,ignd)     !<soil temperature, k
      real thliq(ilg,ignd)    !<liquid soil moisture content in 3 soil layers
      real zbotw(ilg,ignd)    !<bottom of soil layers
      real litres(ilg)        !<litter respiration over the given unvegetated sub-area in umol co2/m2.s
      real socres(ilg)        !<soil c respiration over the given unvegetated sub-area in umol co2/m2.s
      real frac(ilg)          !<fraction of ground (fg) or snow over ground (fgs)

      real delzw(ilg,ignd)  !<
      real thiceg(ilg,ignd) !<
      real zcarb_g          !<

      real litrq10          !<
      real soilcq10         !<
      real litrtemp(ilg)    !<litter temperature
      real solctemp(ilg)    !<soil carbon pool temperature
      real q10func          !<
      real psisat(ilg,ignd) !<saturation matric potential
      real grksat(ilg,ignd) !<saturation hyd. conductivity
      real b(ilg,ignd)      !<parameter b of clapp and hornberger
      real thpor(ilg,ignd)  !<porosity
      real beta             !<
      real fracarb(ilg,ignd)!<fraction of carbon in each soil layer
      real zcarbon          !<
      real tempq10l(ilg)    !<
      real socmoscl(ilg)    !<soil moisture scalar for soil carbon decomposition
      real scmotrm(ilg,ignd)!<
      real ltrmoscl(ilg)    !<soil moisture scalar for litter decomposition
      real psi(ilg,ignd)    !<
      real tempq10s(ilg)    !<
      real fcoeff           !<
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

!>     Estimate temperature of the litter and soil carbon pools.

!!     Over the bare fraction there is no live root. So we make the
!!     simplest assumption that litter temperature is same as temperature of the top soil layer.

      do 210 i = il1, il2
        litrtemp(i)=tbar(i,1)
210   continue

!>     We estimate the temperature of the soil c pool assuming that soil carbon over the bare fraction is distributed exponentially. note
!!     that bare fraction may contain dead roots from different pfts all of which may be distributed differently. For simplicity we do not
!!     track each pft's dead root biomass and assume that distribution of soil carbon over the bare fraction can be described by a single
!!     parameter.

      do 240 i = il1, il2

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
        solctemp(i) = sum(tbar(i,:)*fracarb(i,:))/ sum(fracarb(i,:))

240   continue

!>     find moisture scalar for soil c decomposition

!!     this is modelled as function of logarithm of matric potential.
!!     we find values for all soil layers, and then find an average value
!!     based on fraction of carbon present in each layer.

      do 260 j = 1, ignd
        do 270 i = il1, il2

          if(isand(i,j).eq.-3.or.isand(i,j).eq.-4)then
            scmotrm (i,j)=0.2
            psi (i,j) = 10000.0 ! set to large number so that
                               ! ltrmoscl becomes 0.2

          else ! i.e., sand.ne.-3 or -4

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

!>     find moisture scalar for litter decomposition

!!     the difference between moisture scalar for litter and soil c
!!     is that the litter decomposition is not constrained by high
!!     soil moisture (assuming that litter is always exposed to air).
!!     in addition, we use moisture content of the top soil layer
!!     as a surrogate for litter moisture content. so we use only
!!     psi(i,1) calculated in loops 260 and 270 above.

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

!!     use temperature of the litter and soil c pools, and their soil
!!     moisture scalars to find respiration rates from these pools

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
!>@}
! ------------------------------------------------------------------------------------

!>\ingroup hetresv
!!@{
subroutine hetresv ( fcan,      fct, litrmass, soilcmas, &
                      delzw,  thpor, il1, &
                      il2,     tbar,   psisat, b, thliq,  &
                     roottemp,    zbotw,     sort, &
                     isand, thicec, ipeatland, &
!    -------------- inputs above this line, outputs below -------------
                    ltresveg, scresveg)

!               Canadian Terrestrial Ecosystem Model (CTEM)
!           Heterotrophic Respiration Subtoutine For Vegetated Fraction

!     16  oct. 2001 - this subroutine calculates heterotrophic respiration
!     V. Arora        for a given sub-area, from litter and soil carbon
!                     pools.

!     change history:

!     10  April 2015 -Bring in peatland scheme
!     Y. Wu
!
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

      use ctem_params,        only : icc, ilg, ignd, kk, zero, bsratelt,&
                               bsratesc, abar, tanhq10,&
                               alpha_hetres

      implicit none

      integer il1       !<il1=1
      integer il2       !<il2=ilg
      integer i, j, k
      integer sort(icc) !<index for correspondence between 9 pfts and 12 values in the parameters vectors
      integer isand(ilg,ignd) !<
      integer  ipeatland(ilg) !<

      real fcan(ilg,icc)      !<fractional coverage of ctem's 9 pfts
      real fct(ilg)           !<sum of all fcan, fcan & fct are not used at this time but could be used at some later stage
      real litrmass(ilg,icc+1)!<litter mass for the 9 pfts + bare in \f$kg c/m^2\f$
      real tbar(ilg,ignd)     !<soil temperature, k
      real soilcmas(ilg,icc+1)!<soil carbon mass for the 9 pfts + bare in \f$kg c/m^2\f$
      real thliq(ilg,ignd)    !<liquid soil moisture content in 3 soil layers
      real roottemp(ilg,icc)  !<root temperature as estimated in mainres subroutine
      real zbotw(ilg,ignd)    !<bottom of soil layers
      real ltresveg(ilg,icc)  !<litter respiration for the given sub-area in umol co2/m2.s, for ctem's 9 pfts
      real scresveg(ilg,icc)  !<soil carbon respiration for the given sub-area in umol co2/m2.s, for ctem's 9 pfts
      real thicec(ilg,ignd)   !<liquid soil moisture content in 3 soil layers in canopy covered subarea

      real delzw(ilg,ignd)    !<
      real zcarb_g            !<

      real litrq10            !<
      real soilcq10           !<
      real litrtemp(ilg,icc)  !<litter temperature
      real solctemp(ilg,icc)  !<soil carbon pool temperature
      real q10func            !<
      real psisat(ilg,ignd)   !<saturation matric potential
      real grksat(ilg,ignd)   !<saturation hyd. conductivity
      real b(ilg,ignd)        !<parameter b of clapp and hornberger
      real thpor(ilg,ignd)    !<porosity
      real fracarb(ilg,icc,ignd) !<fraction of carbon in each soil layer for each vegetation
      real zcarbon            !<
      real tempq10l(ilg,icc)  !<
      real socmoscl(ilg,icc)  !<soil moisture scalar for soil carbon decomposition
      real scmotrm(ilg,ignd)  !<soil carbon moisture term
      real ltrmoscl(ilg)      !<soil moisture scalar for litter decomposition
      real psi(ilg,ignd)      !<
      real tempq10s(ilg,icc)  !<
      real fcoeff             !<


!     ------------------------------------------------------------------
!!     Constants and parameters are located in ctem_params.f90
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

!<     estimate temperature of the litter and soil carbon pools. litter
!!     temperature is weighted average of temperatue of top soil layer
!!     (where the stem and leaf litter sits) and root temperature, because
!!     litter pool is made of leaf, stem, and root litter.

      do 200 j = 1,icc
        do 210 i = il1, il2
         if (fcan(i,j) .gt. 0.) then
          litrtemp(i,j)=alpha_hetres*tbar(i,1)+roottemp(i,j)*(1.0-alpha_hetres)
         endif
210     continue
200   continue

!<     estimation of soil carbon pool temperature is not straight forward.
!!     ideally soil c pool temperature should be set same as root temperature,
!!     since soil c profiles are similar to root profiles. but in the event
!!     when the roots die then we may run into trouble. so we find the
!!     temperature of the soil c pool assuming that soil carbon is
!!     exponentially distributed, just like roots. but rather than using
!!     the parameter of this exponential profile from our variable root
!!     distribution we use fixed vegetation-dependent parameters.

      do 230 j = 1, icc
        do 240 i = il1, il2
         if (fcan(i,j) .gt. 0.) then

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

          solctemp(i,j) = sum(tbar(i,:)*fracarb(i,j,:))/ sum(fracarb(i,j,:))

        endif
240     continue
230   continue

!!     find moisture scalar for soil c decomposition

!!     this is modelled as function of logarithm of matric potential.
!!     we find values for all soil layers, and then find an average value
!!     based on fraction of carbon present in each layer. this makes
!!     moisture scalar a function of vegetation type.

      do 260 j = 1, ignd
        do 270 i = il1, il2

          if(isand(i,j).eq.-3.or.isand(i,j).eq.-4)then
            scmotrm (i,j)=0.2
            psi (i,j) = 10000.0 ! set to large number so that
                               ! ltrmoscl becomes 0.2
          else ! i.e., sand.ne.-3 or -4

            !the 0.005 below prevents a divide by 0 situation.
            psi(i,j)   = psisat(i,j)*(thliq(i,j)/(thpor(i,j)+0.005 -thicec(i,j)))**(-b(i,j))

!         FLAG- check on this as I had to change a fair amount what YW had, JM. Sep 21 2016.
!           Also not sure if it is needed?
!               JM- Turn off for now, we'll see how testing looks. Nov 2016.
!             if (ipeatland(i) >0) then
!                 if (thliq(i,j)+ thicec(i,j)+0.01 < thpor(i,j) &
!                    .and.  tbar(i,j) <273.16)                   then
!                   psi(i,j) = 0.001
!                 elseif (thicec(i,j) > thpor(i,j))    then
!                   psi(i,j) = 0.001   !set to saturation
!                 !else
!                     ! leave as-is.
!                 endif
!             endif

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
          socmoscl(i,j)=max(0.2,min(1.0,socmoscl(i,j)))

         endif
290     continue
280   continue

!!     find moisture scalar for litter decomposition

!!     the difference between moisture scalar for litter and soil c
!!     is that the litter decomposition is not constrained by high
!!     soil moisture (assuming that litter is always exposed to air).
!!     in addition, we use moisture content of the top soil layer
!!     as a surrogate for litter moisture content. so we use only
!!     psi(i,1) calculated in loops 260 and 270 above.

      do 300 i = il1, il2
      if (ipeatland(i) == 0)        then   !not peatland
        if(psi(i,1).gt.10000.0) then
          ltrmoscl(i)=0.2
        else if( psi(i,1).le.10000.0 .and. psi(i,1).gt.6.0 ) then
          ltrmoscl(i)=1.0 -  0.8*( (log10(psi(i,1)) - log10(6.0))/(log10(10000.0)-log10(6.0)) )
        else if( psi(i,1).le.6.0 ) then
          ltrmoscl(i)=1.0
        endif
        ltrmoscl(i)=max(0.2,min(1.0,ltrmoscl(i)))
      else  !is peatland
!
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
          if (psi(i,1).ge. 10000.0) then
               ltrmoscl(i) = 0.2
          elseif (psi(i,1).le. 10000.0 .and.psi(i,1).gt. 6.0) then
               ltrmoscl(i)=1.0 - 0.8*((log10(psi(i,1))-log10(6.0)) &
                        /(log10(10000.0)-log10(6.0)))**1.
          elseif (psi(i,1).le. 6.0 .and. psi(i,1) .gt. 4.0) then
               ltrmoscl(i)=1.0
          elseif (psi(i,1).le. 4.0 .and. psi(i,1).gt.psisat(i,1))  then
               ltrmoscl(i)=1.0-0.99*((log10(4.0)-log10(psi(i,1)))/ &
                        (log10(4.0)-log10(psisat(i,1))))
          elseif (psi(i,1) .le. psisat(i,1))                     then
               ltrmoscl(i)=0.01
          endif
          ltrmoscl(i)=max(0.0,min(ltrmoscl(i),1.0))
        endif  !peatland
300   continue

!!    use temperature of the litter and soil c pools, and their soil
!!     moisture scalars to find respiration rates from these pools

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
!!@}
end module heterotrophic_respiration