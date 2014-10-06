      subroutine mortalty (stemmass, rootmass,    ailcg, gleafmas,
     1                     bleafmas,     il1, 
     2                          il2,     iday,   do_age_mort, sort,
     3                      fcancmx,
c    + ------------------ inputs above this line ----------------------   
     4                     lystmmas, lyrotmas, tymaxlai, grwtheff,
c    + -------------- inputs updated above this line ------------------
     5                     stemltrm, rootltrm, glealtrm, geremort,
     6                     intrmort)
c    + ------------------outputs above this line ----------------------
c
c               Canadian Terrestrial Ecosystem Model (CTEM) 
C                             Mortality Subroutine
c
c     17  Jan 2014  - Moved parameters to global file (ctem_params.f90)
c     J. Melton
c
c     22  Jul 2013  - Add in module for parameters
C     J. Melton
c
c     24  sep 2012  - add in checks to prevent calculation of non-present
c     j. melton       pfts
c
c     07  may 2003  - this subroutine calculates the litter generated
c     v. arora        from leaves, stem, and root components after
c                     vegetation dies due to reduced growth efficiency
c                     or due to aging (the intrinsic mortality)  

c     inputs 
c
c     stemmass  - stem mass for each of the 9 ctem pfts, kg c/m2
c     rootmass  - root mass for each of the 9 ctem pfts, kg c/m2
c     ailcg     - green or live lai
c     gleafmas  - green leaf mass for each of the 9 ctem pfts, kg c/m2
c     bleafmas  - brown leaf mass for each of the 9 ctem pfts, kg c/m2
c     lystmmas  - stem mass at the end of last year
c     lyrotmas  - root mass at the end of last year
c     tymaxlai  - this year's maximum lai
c     grwtheff  - growth efficiency. change in biomass per year per
c                 unit max. lai (g c/m2)/(m2/m2)
c     icc       - no. of ctem plant function types, currently 8
c     ilg       - no. of grid cells in latitude circle
c     il1,il2   - il1=1, il2=ilg
c     iday      - day of the year
c     do_age_mort - switch to control calc of age mortality (false=not done)
c     sort      - index for correspondence between ctem 9 pfts and size
c                 12 of parameters vectors
c
c     outputs
c
c     stemltrm  - stem litter generated due to mortality (kg c/m2)
c     rootltrm  - root litter generated due to mortality (kg c/m2)
c     glealtrm  - green leaf litter generated due to mortality (kg c/m2)
c     geremort  - growth efficiency related mortality (1/day)
c     intrmort  - intrinsic mortality (1/day)

      use ctem_params,        only : icc, ilg, kk, zero, mxmortge, 
     1                               kmort1, maxage
c
      implicit none
c
      integer il1, il2, i, j, k, iday, n
c
      integer       sort(icc)
c
      logical      do_age_mort
c
      real  stemmass(ilg,icc), rootmass(ilg,icc), gleafmas(ilg,icc),
     1         ailcg(ilg,icc), grwtheff(ilg,icc), lystmmas(ilg,icc),
     2      lyrotmas(ilg,icc), tymaxlai(ilg,icc), bleafmas(ilg,icc)
c
      real  stemltrm(ilg,icc), rootltrm(ilg,icc), glealtrm(ilg,icc),
     1      geremort(ilg,icc), intrmort(ilg,icc), fcancmx(ilg,icc)
c
c     ------------------------------------------------------------------
c     Constants and parameters are located in ctem_params.f90
c     ---------------------------------------------------------------
c
c     initialize required arrays to zero
c
      do 140 j = 1,icc
        do 150 i = il1, il2
          stemltrm(i,j)=0.0     !stem litter due to mortality
          rootltrm(i,j)=0.0     !root litter due to mortality
          glealtrm(i,j)=0.0     !green leaf litter due to mortality
          geremort(i,j)=0.0     !growth efficiency related mortality rate 
          intrmort(i,j)=0.0     !intrinsic mortality rate 
150     continue                  
140   continue
c
c     initialization ends    
c
c     ------------------------------------------------------------------
c
c     at the end of every year, i.e. when iday equals 365, we calculate
c     growth related mortality. rather than using this number to kill
c     plants at the end of every year, this mortality rate is applied
c     gradually over the next year.
c
      do 200 j = 1, icc
        n = sort(j)
        do 210 i = il1, il2
         if (fcancmx(i,j).gt.0.0) then
          if(iday.eq.1)then
            tymaxlai(i,j) =0.0
          endif
c
          if(ailcg(i,j).gt.tymaxlai(i,j))then
            tymaxlai(i,j)=ailcg(i,j)
          endif
c
          if(iday.eq.365)then
            if(tymaxlai(i,j).gt.zero)then
              grwtheff(i,j)= ( (stemmass(i,j)+rootmass(i,j))-
     &         (lystmmas(i,j)+lyrotmas(i,j)) )/tymaxlai(i,j) 
            else
              grwtheff(i,j)= 0.0
            endif
            grwtheff(i,j)=max(0.0,grwtheff(i,j))*1000.0
            lystmmas(i,j)=stemmass(i,j)
            lyrotmas(i,j)=rootmass(i,j)
          endif
c
c         calculate growth related mortality using last year's growth
c         efficiency or the new growth efficiency if day is 365 and
c         growth efficiency estimate has been updated above.
c
          geremort(i,j)=mxmortge(n)/(1.0+kmort1*grwtheff(i,j))
c
c         convert (1/year) rate into (1/day) rate   
          geremort(i,j)=geremort(i,j)/365.0
         endif
210     continue
200   continue
c
c     calculate intrinsic mortality rate due to aging which implicity
c     includes effects of frost, hail, wind throw etc. it is assumed 
c     that only 1% of the plants exceed maximum age (which is a pft-
c     dependent parameter). to achieve this some fraction of the plants
c     need to be killed every year. 
c
      do 250 j = 1, icc
        n = sort(j)
        do 260 i = il1, il2
         if (fcancmx(i,j).gt.0.0) then

!          if(do_age_mort)then

!             if(maxage(n).gt.zero)then
               intrmort(i,j)=1.0-exp(-4.605/maxage(n))
!             else
!               intrmort(i,j)=0.0
!             endif

!          else
!            geremort(i,j)=0.0
!            intrmort(i,j)=0.0
!          end if

c         convert (1/year) rate into (1/day) rate   
          intrmort(i,j)=intrmort(i,j)/365.0
         endif
260     continue
250   continue 
c
c     now that we have both growth related and intrinsic mortality rates,
c     lets combine these rates for every pft and estimate litter generated
c
      do 300 j = 1, icc
        do 310 i = il1, il2
         if (fcancmx(i,j).gt.0.0) then
          stemltrm(i,j)=stemmass(i,j)*
     &    ( 1.0-exp(-1.0*(geremort(i,j)+intrmort(i,j))) )
          rootltrm(i,j)=rootmass(i,j)*
     &    ( 1.0-exp(-1.0*(geremort(i,j)+intrmort(i,j))) )
          glealtrm(i,j)=gleafmas(i,j)*
     &    ( 1.0-exp(-1.0*(geremort(i,j)+intrmort(i,j))) )
         endif
310     continue
300   continue
c
      return
      end

