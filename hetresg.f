      subroutine hetresg (litrmass, soilcmas,         
     1                         il1,      il2,     tbar,    
     2                       thliq,     sand,      clay,   zbotw,   
     3                        frac,    isnow,      isand,   ipeatland, !YW
c    -------------- inputs above this line, outputs below -------------
     4                      litres,   socres)  
c
C               Canadian Terrestrial Ecosystem Model (CTEM) 
C           Heterotrophic Respiration Subroutine For Bare Fraction
c
c     11  Apr. 2003 - this subroutine calculates heterotrophic respiration
c     V. Arora        over the bare subarea of a grid cell (i.e. ground only
c                     and snow over ground subareas).
c
c     change history:
c
c     17  Jan 2014  - Moved parameters to global file (ctem_params.f90)
c     J. Melton
c
c     23  Jul 2013  - add in module for parameters
c     J. Melton
c     J. Melton and V.Arora - changed tanhq10 parameters, they were switched
c               25 Sep 2012
c     J. Melton 31 Aug 2012 - remove isnow, it is not used.
c     J. Melton 23 Aug 2012 - bring in isand, converting sand to
c                             int was missing some gridcells assigned
c                             to bedrock in classb

c     ------
c     inputs 
c
c     litrmass  - litter mass for the 8 pfts + bare in kg c/m2
c     soilcmas  - soil carbon mass for the 8 pfts + bare in kg c/m2
c     icc       - no. of vegetation types (currently 8)
c     ignd        - no. of soil layers (currently 3)
c     ilg       - no. of grid cells in latitude circle
c     il1,il2   - il1=1, il2=ilg
c     tbar      - soil temperature, k
c     thliq     - liquid soil moisture content in 3 soil layers
c     sand      - percentage sand
c     clay      - percentage clay
c     zbotw     - bottom of soil layers
c     frac      - fraction of ground (fg) or snow over ground (fgs)
c     isnow     - integer telling if bare fraction is fg (0) or fgs (1)
c
c     outputs
c
c     litres    - litter respiration over the given unvegetated sub-area
c                 in umol co2/m2.s
c     socres    - soil c respiration over the given unvegetated sub-area
c                 in umol co2/m2.s
c
      use ctem_params,        only : icc, ilg, ignd, zero, tanhq10, a,
     1                               bsratelt_g, bsratesc_g

      implicit none
c
c     isnow is changed to isnow(ilg) in classt of class version higher 
c     than 3.4 for coupling with ctem
c
      integer il1,il2,i,j,k,isnow,isand(ilg,ignd),ipeatland(ilg)

      real                  litrmass(ilg,icc+1),  soilcmas(ilg,icc+1), 
     1           tbar(ilg,ignd),    thliq(ilg,ignd),    sand(ilg,ignd), 
     2          zbotw(ilg,ignd),      litres(ilg),          socres(ilg),
     3           clay(ilg,ignd),        frac(ilg)
      
      real              litrq10,           soilcq10,               
     2        litrtemp(ilg),      solctemp(ilg),             q10func,
     3       psisat(ilg,ignd),     grksat(ilg,ignd),        b(ilg,ignd),
     4        thpor(ilg,ignd),               beta,
     5      fracarb(ilg,ignd),         zcarbon,
     6        tempq10l(ilg),      socmoscl(ilg),     scmotrm(ilg,ignd),
     7        ltrmoscl(ilg),        psi(ilg,ignd),       tempq10s(ilg),
     8               fcoeff

c     ------------------------------------------------------------------
c     Constants and parameters are located in ctem_params.f90
c     ---------------------------------------------------------------
c
c     initialize required arrays to zero
c
      do 100 k = 1, ignd
        do 100 i = il1, il2
          fracarb(i,k)=0.0  ! fraction of carbon in each soil layer
100   continue
c
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
c
      do 120 j = 1, ignd
        do 130 i = il1, il2
          psisat(i,j) = 0.0       ! saturation matric potential
          grksat(i,j) = 0.0       ! saturation hyd. conductivity
          thpor(i,j) = 0.0        ! porosity
          b(i,j) = 0.0            ! parameter b of clapp and hornberger
130     continue
120   continue
c
c     initialization ends    
c
c     ------------------------------------------------------------------
c
c     estimate temperature of the litter and soil carbon pools. 
c
c     over the bare fraction there is no live root. so we make the
c     simplest assumption that litter temperature is same as temperature
c     of the top soil layer.
c     
      do 210 i = il1, il2
        litrtemp(i)=tbar(i,1)
210   continue
c
c     we estimate the temperature of the soil c pool assuming that soil 
c     carbon over the bare fraction is distributed exponentially. note
c     that bare fraction may contain dead roots from different pfts all of
c     which may be distributed differently. for simplicity we do not
c     track each pft's dead root biomass and assume that distribution of
c     soil carbon over the bare fraction can be described by a single
c     parameter.
c
      do 240 i = il1, il2
c
        zcarbon=3.0/a                 ! 95% depth
        if(zcarbon.le.zbotw(i,1)) then
            fracarb(i,1)=1.0             ! fraction of carbon in
            fracarb(i,2)=0.0             ! soil layers
            fracarb(i,3)=0.0
        else
            fcoeff=exp(-a*zcarbon)
            fracarb(i,1)=
     &        1.0-(exp(-a*zbotw(i,1))-fcoeff)/(1.0-fcoeff)
            if(zcarbon.le.zbotw(i,2)) then
                fracarb(i,2)=1.0-fracarb(i,1)
                fracarb(i,3)=0.0
            else
                fracarb(i,3)=
     &            (exp(-a*zbotw(i,2))-fcoeff)/(1.0-fcoeff)
                fracarb(i,2)=1.0-fracarb(i,1)-fracarb(i,3)
            endif
        endif
c
        solctemp(i)=tbar(i,1)*fracarb(i,1) +
     &     tbar(i,2)*fracarb(i,2) +
     &     tbar(i,3)*fracarb(i,3)
        solctemp(i)=solctemp(i) /
     &     (fracarb(i,1)+fracarb(i,2)+fracarb(i,3))
c
c
c       make sure we don't use temperatures of 2nd and 3rd soil layers
c       if they are specified bedrock via sand -3 flag
c
        if(isand(i,3).eq.-3)then ! third layer bed rock
          solctemp(i)=tbar(i,1)*fracarb(i,1) +
     &       tbar(i,2)*fracarb(i,2) 
          solctemp(i)=solctemp(i) /
     &       (fracarb(i,1)+fracarb(i,2))
        endif
        if(isand(i,2).eq.-3)then ! second layer bed rock
          solctemp(i)=tbar(i,1)
        endif
c
240   continue     
c
c     find moisture scalar for soil c decomposition
c
c     this is modelled as function of logarithm of matric potential. 
c     we find values for all soil layers, and then find an average value 
c     based on fraction of carbon present in each layer. 
c
      do 260 j = 1, ignd
        do 270 i = il1, il2
c
          if(isand(i,j).eq.-3.or.isand(i,j).eq.-4)then
            scmotrm (i,j)=0.2
            psi (i,j) = 10000.0 ! set to large number so that
c                               ! ltrmoscl becomes 0.2
          else ! i.e., sand.ne.-3 or -4
            psisat(i,j)= (10.0**(-0.0131*sand(i,j)+1.88))/100.0
            b(i,j)     = 0.159*clay(i,j)+2.91
            thpor(i,j) = (-0.126*sand(i,j)+48.9)/100.0
            psi(i,j)   = psisat(i,j)*(thliq(i,j)/thpor(i,j))**(-b(i,j)) 
c   
            if(psi(i,j).ge.10000.0) then
              scmotrm(i,j)=0.2
            else if( psi(i,j).lt.10000.0 .and.  psi(i,j).gt.6.0 ) then
              scmotrm(i,j)=1.0 - 0.8*
     &   ( (log10(psi(i,j)) - log10(6.0))/(log10(10000.0)-log10(6.0)) )         
            else if( psi(i,j).le.6.0 .and.  psi(i,j).ge.4.0 ) then
              scmotrm(i,j)=1.0
            else if( psi(i,j).lt.4.0.and.psi(i,j).gt.psisat(i,j) )then 
              scmotrm(i,j)=1.0 - 
     &          0.5*( (log10(4.0) - log10(psi(i,j))) / 
     &         (log10(4.0)-log10(psisat(i,j))) )
            else if( psi(i,j).le.psisat(i,j) ) then
              scmotrm(i,j)=0.5
            endif
          endif ! if sand.eq.-3 or -4
c
          scmotrm(i,j)=max(0.0,min(scmotrm(i,j),1.0))
270     continue     
260   continue     
c
      do 290 i = il1, il2
        socmoscl(i) = scmotrm(i,1)*fracarb(i,1) + 
     &     scmotrm(i,2)*fracarb(i,2) +
     &     scmotrm(i,3)*fracarb(i,3)
        socmoscl(i) = socmoscl(i) /
     &     (fracarb(i,1)+fracarb(i,2)+fracarb(i,3))    
c
c       make sure we don't use scmotrm of 2nd and 3rd soil layers
c       if they are specified bedrock via sand -3 flag
c
        if(isand(i,3).eq.-3)then ! third layer bed rock
          socmoscl(i) = scmotrm(i,1)*fracarb(i,1) +
     &                    scmotrm(i,2)*fracarb(i,2)
          socmoscl(i) = socmoscl(i) /
     &     (fracarb(i,1)+fracarb(i,2))
        endif
        if(isand(i,2).eq.-3)then ! second layer bed rock
          socmoscl(i) = scmotrm(i,1)
        endif
c
        socmoscl(i)=max(0.2,min(socmoscl(i),1.0))
290   continue     
c
c     find moisture scalar for litter decomposition
c
c     the difference between moisture scalar for litter and soil c
c     is that the litter decomposition is not constrained by high
c     soil moisture (assuming that litter is always exposed to air).
c     in addition, we use moisture content of the top soil layer
c     as a surrogate for litter moisture content. so we use only 
c     psi(i,1) calculated in loops 260 and 270 above.
c
      do 300 i = il1, il2
c     peatland branch where litter respiration can be limtied by water\

        if (ipeatland(i) == 0)        then           
          if(psi(i,1).gt.10000.0) then
            ltrmoscl(i)=0.2
          else if( psi(i,1).le.10000.0 .and.  psi(i,1).gt.6.0 ) then
            ltrmoscl(i)=1.0 - 0.8*
     &    ( (log10(psi(i,1)) - log10(6.0))/(log10(10000.0)-log10(6.0)) )
          else if( psi(i,1).le.6.0 ) then
            ltrmoscl(i)=1.0 
          endif
          ltrmoscl(i)=max(0.2,min(ltrmoscl(i),1.0))
        else                       !is peatland
          if (psi(i,1).ge. 10000.0) then
	    		ltrmoscl(i) = 0.2
        	elseif (psi(i,1).le.10000.0 .and.psi(i,1).gt.6.0) then
	          ltrmoscl(i)=1.0 - 0.8*((log10(psi(i,1))-log10(6.0))
	1				/(log10(10000.0)-log10(6.0)) )
         	elseif (psi(i,1).le.6.0 .and. psi(i,1) .ge. 4.0) then
              	ltrmoscl(i)=1.0
        	elseif (psi(i,1).lt.4.0 .and. psi(i,1).gt.psisat(i,1))  then 
              	ltrmoscl(i)=1.0-0.2*((log10(4.0)-log10(psi(i,1)))/   
     1         		(log10(4.0)-log10(psisat(i,1))))
         	elseif (psi(i,1) .le. psisat(i,1)) 				then
             	ltrmoscl(i)=0.2		
        	endif
         	ltrmoscl(i)=max(0.0,min(ltrmoscl(i),1.0))
        endif  !peatland  
c    -------------------YW March 30, 2015 -----------------------------/
          
300   continue
c
c     use temperature of the litter and soil c pools, and their soil
c     moisture scalars to find respiration rates from these pools
c
      do 330 i = il1, il2
      if(frac(i).gt.zero)then
c
c       first find the q10 response function to scale base respiration
c       rate from 15 c to current temperature, we do litter first
c
        tempq10l(i)=litrtemp(i)-273.16
        litrq10 = tanhq10(1) + tanhq10(2)*
     &            ( tanh( tanhq10(3)*(tanhq10(4)-tempq10l(i))  ) )
c
        q10func = litrq10**(0.1*(litrtemp(i)-273.16-15.0))
        litres(i)= ltrmoscl(i) * litrmass(i,icc+1)*
     &    bsratelt_g*2.64*q10func ! 2.64 converts bsratelt_g from kg c/kg c.year
c                                  ! to u-mol co2/kg c.s
c
c       respiration from soil c pool
c
        tempq10s(i)=solctemp(i)-273.16
        soilcq10= tanhq10(1) + tanhq10(2)*
     &            ( tanh( tanhq10(3)*(tanhq10(4)-tempq10s(i))  ) )
c
        q10func = soilcq10**(0.1*(solctemp(i)-273.16-15.0))
        socres(i)= socmoscl(i)* soilcmas(i,icc+1)*
     &    bsratesc_g*2.64*q10func ! 2.64 converts bsratesc_g from kg c/kg c.year
c                                  ! to u-mol co2/kg c.s
c
      endif
330   continue
c
      return
      end

