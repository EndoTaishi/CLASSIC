      subroutine mainres (  fcan,      fct,     stemmass,   rootmass, 
     1                       il1,
     2                       il2,     tcan,         tbar,   rmatctem,
     3                      sort, nol2pfts,        isand,
c    -------------- inputs above this line, outputs below ----------
     4                      rmsveg, rmrveg,     roottemp)
c
c
C               Canadian Terrestrial Ecosystem Model (CTEM)
C                    Maintenance Respiration Subtoutine
c
c     20  sep. 2001 - this subroutine calculates maintenance respiration,
c     V. Arora        over a given sub-area, for stem and root components.
c                     leaf respiration is estimated within the phtsyn
c                     subroutine.

c     change history:

c     J. Melton 22  Jul 2013 - Add in module for parameters
C     
c     J. Melton 20 sep 2012 - made it so does not do calcs for pfts with
c                             fcan = 0.
c     J. Melton 23 aug 2012 - change sand to isand, converting sand to
c                             int was missing some gridcells assigned
c                             to bedrock in classb. isand is now passed
c                             in.
c
c     inputs 
c
c     fcan      - fractional coverage of ctem's 9 pfts over the given
c                 sub-area
c     fct       - sum of all fcan
c                 fcan & fct are not used at this time but could
c                 be used at some later stage
c     stemmass  - stem biomass for the 9 pfts in kg c/m2
c     rootmass  - root biomass for the 9 pfts in kg c/m2
c     icc       - no. of ctem pfts (currently 9)
c     ignd        - no. of soil layers (currently 3)
c     ilg       - no. of grid cells in latitude circle
c     il1,il2   - il1=1, il2=ilg
c     tcan      - canopy temperature, k
c     tbar      - soil temperature, k
c     rmatctem  - fraction of roots in each layer for each pft
c     sort      - index for correspondence between 9 pfts and 12 values
c                 in the parameter vectors
c     nol2pfts  - number of level 2 ctem pfts
c     ican        - number of class pfts, currently 4
c     isand      - flag for bedrock or ice in a soil layer
c
c     outputs 
c
c     rmsveg    - maintenance respiration for stem for the 9 pfts
c     rmrveg    - maintenance respiration for root for the 9 pfts
c                 both in u mol co2/m2. sec
c     roottemp  - root temperature (k)
c
      use ctem_params,        only : icc, ilg, ignd, ican, kk, zero

      implicit none
c
      integer il1, il2, i, j, k, sort(icc), 
     1        nol2pfts(ican),   k1,   k2,  m,  isand(ilg,ignd)
c
      real  fcan(ilg,icc),         fct(ilg),      stemmass(ilg,icc), 
     1          tcan(ilg),     tbar(ilg,ignd),      rootmass(ilg,icc),   
     2    rmsveg(ilg,icc),  rmrveg(ilg,icc),   rmatctem(ilg,icc,ignd)   
c
      real  bsrtstem(kk),      bsrtroot(kk),      tempq10r(ilg,icc), 
     3     tempq10s(ilg), roottemp(ilg,icc),                    q10, 
     4           q10func, livstmfr(ilg,icc),
     5 livrotfr(ilg,icc),           minlvfr
c
      logical consq10
c
c     constants and parameters 
c
c     note the structure of vectors which clearly shows the class
c     pfts (along rows) and ctem sub-pfts (along columns)
c
c     needle leaf |  evg       dcd       ---
c     broad leaf  |  evg   dcd-cld   dcd-dry
c     crops       |   c3        c4       ---
c     grasses     |   c3        c4       ---
c
c     ---------------------------------------------------
c
c     base respiration rates for stem and root for ctem pfts in
c     kg c/kg c.year at 15 degrees celcius. note that maintenance
c     respiration rates of root are higher because they contain
c     both wood (coarse roots) and fine roots.

c    new parameter values to produce carbon use efficiencies more in
c    line with literature (zhang et al. 2009, luyssaert et al. gcb 2007)
c    values changed for bsrtstem and bsrtroot. jm 06.2012

      data bsrtstem/0.0900, 0.0550, 0.0000,
     &              0.0600, 0.0335, 0.0300,
     &              0.0365, 0.0365, 0.0000,
     &              0.0000, 0.0000, 0.0000/ ! no stem component for grasses

      data bsrtroot/0.5000, 0.2850, 0.0000,
     &              0.6500, 0.2250, 0.0550,
     &              0.1600, 0.1600, 0.0000,
     &              0.1000, 0.1000, 0.0000/
c
c     set the following switch to .true. for using constant temperature
c     indepedent q10 specified below
      data consq10 /.false./
c
c     q10 - if using a constant temperature independent value, i.e.
c     if consq10 is set to true
      data q10/2.00/
c
c     minimum live wood fraction
      data minlvfr/0.05/
c
c     ---------------------------------------------------
c
c     initialize required arrays to zero
c
      do 100 j = 1, icc
        do 110 i = il1, il2
          roottemp(i,j) = 0.0        ! root temperature
          rmsveg(i,j) = 0.0          ! stem maintenance respiration
          rmrveg(i,j) = 0.0          ! root maintenance respiration
          livstmfr(i,j)= 0.0         ! live stem fraction
          livrotfr(i,j)= 0.0         ! live root fraction
110     continue 
100   continue 
c
c     initialization ends
c
c     based on root and stem biomass, find fraction which is live.
c     for stem this would be the sapwood to total wood ratio.
c
      k1=0
      do 120 j = 1, ican
        if(j.eq.1) then
          k1 = k1 + 1
        else
          k1 = k1 + nol2pfts(j-1)
        endif
        k2 = k1 + nol2pfts(j) - 1
        do 125 m = k1, k2
         do 130 i = il1, il2
          if(j.le.2)then     ! trees
            livstmfr(i,m) = exp(-0.2835*stemmass(i,m))  !following century model              
            livstmfr(i,m) = max(minlvfr,min(livstmfr(i,m),1.0))
            livrotfr(i,m) = exp(-0.2835*rootmass(i,m))               
            livrotfr(i,m) = max(minlvfr,min(livrotfr(i,m),1.0))
          else                 ! crop and grass are all live
            livstmfr(i,m) = 1.0
            livrotfr(i,m) = 1.0
          endif
130     continue 
125    continue 
120   continue 
c
c     fraction of roots for each vegetation type, for each soil layer, 
c     in each grid cell is given by rmatctem (grid cell, veg type, soil layer) 
c     which bio2str subroutine calculates. rmatctem can thus be used 
c     to find average root temperature for each plant functional type 
c
      do 180 j = 1, icc
        do 190 i = il1, il2
         if (fcan(i,j) .gt. 0.) then
          roottemp(i,j)=tbar(i,1)*rmatctem(i,j,1) + 
     &       tbar(i,2)*rmatctem(i,j,2) +  
     &       tbar(i,3)*rmatctem(i,j,3)
          roottemp(i,j)=roottemp(i,j) /
     &       (rmatctem(i,j,1)+rmatctem(i,j,2)+rmatctem(i,j,3))

c
c        make sure that i do not use temperatures from 2nd and 3rd layers
c        if they are bed rock
c
          if(isand(i,3).eq.-3)then ! third layer bed rock
            roottemp(i,j)=tbar(i,1)*rmatctem(i,j,1) + 
     &        tbar(i,2)*rmatctem(i,j,2)   
            roottemp(i,j)=roottemp(i,j) /
     &        (rmatctem(i,j,1)+rmatctem(i,j,2))
          endif         
          if(isand(i,2).eq.-3)then ! second layer bed rock
            roottemp(i,j)=tbar(i,1)
          endif    
         endif !fcan check.     
190     continue 
180   continue 
c
c     we assume that stem temperature is same as canopy temperature tcan.
c     using stem and root temperatures we can find their maintenance 
c     respirations rates
c
      do 200 i = il1, il2
c
c         first find the q10 response function to scale base respiration
c         rate from 15 c to current temperature, we do the stem first.
c
          if (.not.consq10) then
c           when finding temperature dependent q10, use temperature which
c           is close to average of actual temperature and the temperature
c           at which base rate is specified
            tempq10s(i)=(15.0+273.16+tcan(i))/1.9
            q10 = 3.22 - 0.046*(tempq10s(i)-273.16)       
            q10 = min(4.0, max(1.5, q10))
          endif
c
          q10func = q10**(0.1*(tcan(i)-288.16))
c
        do 210 j = 1, icc
         if (fcan(i,j) .gt. 0.) then

c         this q10 value is then used with the base rate of respiration
c         (commonly taken at some reference temperature (15 deg c), see tjoelker et
c         al. 2009 new phytologist or atkin et al. 2000 new phyto for 
c         an example.). long-term acclimation to temperature could be occuring 
c         see king et al. 2006 nature som for a possible approach. jm.

          rmsveg(i,j)=stemmass(i,j)* livstmfr(i,j)* q10func*
     &     (bsrtstem(sort(j))/365.0)
c
c         convert kg c/m2.day -> u mol co2/m2.sec
          rmsveg(i,j)= rmsveg(i,j) * 963.62
c
c         root respiration
c    
          if (.not.consq10) then
            tempq10r(i,j)=(15.0+273.16+roottemp(i,j))/1.9
            q10 = 3.22 - 0.046*(tempq10r(i,j)-273.16)       
            q10 = min(4.0, max(1.5, q10))
          endif
c
          q10func = q10**(0.1*(roottemp(i,j)-288.16))
          rmrveg(i,j)=rootmass(i,j)* livrotfr(i,j)* q10func*
     &     (bsrtroot(sort(j))/365.0)
c
c         convert kg c/m2.day -> u mol co2/m2.sec
          rmrveg(i,j)= rmrveg(i,j) * 963.62 
c
         endif !fcan check.   
210     continue 
200   continue 
c     
      return
      end

