      subroutine allocate(lfstatus,    thliq,    ailcg,     ailcb, 
     1                         il1,     il2,     sand,     clay,  
     2                    rmatctem, gleafmas, stemmass, rootmass,
     4                        sort, nol2pfts, fcancmx,
c    5 ------------------ inputs above this line ----------------------   
     6                     afrleaf,  afrstem,  afrroot,  wiltsm,
     7                     fieldsm, wtstatus, ltstatus)
c    8 ------------------outputs  above this line ---------------------
c
C               Canadian Terrestrial Ecosystem Model (CTEM)
C                            Allocation Subroutine
C
c     17  Jan 2014  - Moved parameters to global file (ctem_params.f90)
c     J. Melton
c   
c     5   Jul 2013  - Fixed bug with initializing the variables. Brought in
c     J. Melton       the modules for global parameters
c
c     22  Nov 2012  - Calling this version 1.1 since a fair bit of ctem
c     V. Arora        subroutines were changed for compatibility with class
c                     version 3.6 including the capability to run ctem in
c                     mosaic/tile version along with class.
c
c     24  Sep 2012  - Add in checks to prevent calculation of non-present
c     J. Melton       pfts
c
c     05  May 2003  - This subroutine calculates the allocation fractions
c     V. Arora        for leaf, stem, and root components for ctem's pfts 
c
c     inputs 
c
c     lfstatus  - leaf status. an integer indicating if leaves are  
c                 in "max. growth", "normal growth", "fall/harvest",
c                 or "no leaves" mode. see phenolgy subroutine for 
c                 more details.
c     thliq     - liquid soil moisture content in 3 soil layers
c     ailcg     - green or live leaf area index
c     ailcb     - brown or dead leaf area index
c     icc       - no. of ctem plant function types, currently 9
c     ignd        - no. of soil layers (currently 3)
c     ilg       - no. of grid cells in latitude circle
c     il1,il2   - il1=1, il2=ilg
c     sand      - percentage sand
c     clay      - percentage clay
c     rmatctem  - fraction of roots in each soil layer for each pft
c     gleafmas  - green or live leaf mass in kg c/m2, for the 9 pfts
c     stemmass  - stem mass for each of the 9 ctem pfts, kg c/m2
c     rootmass  - root mass for each of the 9 ctem pfts, kg c/m2
c     sort      - index for correspondence between 9 pfts and the
c                 12 values in parameters vectors
c     nol2pfts  - number of level 2 ctem pfts
c     ican        - number of class pfts
c     fcancmx   - max. fractional coverage of ctem's 9 pfts, but this can be
c                modified by land-use change, and competition between pfts
c
c     outputs
c
c     afrleaf   - allocation fraction for leaves
c     afrstem   - allocation fraction for stem
c     afrroot   - allocation fraction for root
c     wiltsm    - wilting point soil moisture content
c     fieldsm   - field capacity soil moisture content
c     wtstatus  - soil water status (0 dry -> 1 wet)
c     ltstatus  - light status

      use ctem_params,        only : eta, kappa, kn, abszero, icc, ilg,
     1                               ignd, kk, ican, omega, epsilonl,
     2                               epsilons, epsilonr, caleaf, castem,
     3                               caroot, consallo, rtsrmin, aldrlfon
c
      implicit none
c
      integer il1, il2, i, j, k,  lfstatus(ilg,icc), 
     1         n, k1,  k2,   m
c
      integer       sort(icc),      nol2pfts(ican)
c
      real     ailcg(ilg,icc),    ailcb(ilg,icc),       thliq(ilg,ignd), 
     1         wiltsm(ilg,ignd),   fieldsm(ilg,ignd), rootmass(ilg,icc),
     2   rmatctem(ilg,icc,ignd), gleafmas(ilg,icc),   stemmass(ilg,icc),
     3           sand(ilg,ignd),      clay(ilg,ignd),  thpor(ilg,ignd),
     4         psisat(ilg,ignd),         b(ilg,ignd),  grksat(ilg,ignd)
c
      real   afrleaf(ilg,icc),  afrstem(ilg,icc),     afrroot(ilg,icc),
     1       fcancmx(ilg,icc)
c
      real  avwiltsm(ilg,icc),  afieldsm(ilg,icc),    avthliq(ilg,icc),
     1      wtstatus(ilg,icc),  ltstatus(ilg,icc),    nstatus(ilg,icc),
     2      wnstatus(ilg,icc),              denom,   mnstrtms(ilg,icc),
     3                   diff,              term1,               term2,
     4         aleaf(ilg,icc),     astem(ilg,icc),      aroot(ilg,icc)
c
c
c     ------------------------------------------------------------------
c     Constants and parameters are located in ctem_params.f90
c     ---------------------------------------------------------------
c
c     initialize required arrays to 0
c
      do 140 j = 1,icc
        do 150 i = il1, il2
          afrleaf(i,j)=0.0    !allocation fraction for leaves
          afrstem(i,j)=0.0    !allocation fraction for stem
          afrroot(i,j)=0.0    !allocation fraction for root
c
            aleaf(i,j)=0.0    !temporary variable
            astem(i,j)=0.0    !temporary variable
            aroot(i,j)=0.0    !temporary variable
c
c                                 !averaged over the root zone
          avwiltsm(i,j)=0.0   !wilting point soil moisture
          afieldsm(i,j)=0.0   !field capacity soil moisture
           avthliq(i,j)=0.0   !liquid soil moisture content
c
          wtstatus(i,j)=0.0   !water status
          ltstatus(i,j)=0.0   !light status
           nstatus(i,j)=0.0   !nitrogen status, if and when we
c                                 !will have n cycle in the model
          wnstatus(i,j)=0.0   !min. of water & n status
c
          mnstrtms(i,j)=0.0   !min. (stem+root) biomass needed to
c                                 !support leaves
150     continue                  
140   continue
c
c     initialization ends    
c
c     ------------------------------------------------------------------
c     Estimate field capacity and wilting point soil moisture contents
c
c     Wilting point corresponds to matric potential of 150 m
c     field capacity corresponds to hydarulic conductivity of
c     0.10 mm/day -> 1.157x1e-09 m/s
c
      do 160 j = 1, ignd
        do 170 i = il1, il2
c
          psisat(i,j)= (10.0**(-0.0131*sand(i,j)+1.88))/100.0
          grksat(i,j)= (10.0**(0.0153*sand(i,j)-0.884))*7.0556e-6
          thpor(i,j) = (-0.126*sand(i,j)+48.9)/100.0
          b(i,j)     = 0.159*clay(i,j)+2.91
c
          wiltsm(i,j) = (150./psisat(i,j))**(-1.0/b(i,j))
          wiltsm(i,j) = thpor(i,j) * wiltsm(i,j)
c
          fieldsm(i,j) = (1.157e-09/grksat(i,j))**
     &      (1./(2.*b(i,j)+3.))
          fieldsm(i,j) = thpor(i,j) *  fieldsm(i,j)
c
170     continue
160   continue
c
c
c     Calculate liquid soil moisture content, and wilting and field capacity 
c     soil moisture contents averaged over the root zone. note that while
c     the soil moisture content is same under the entire gcm grid cell,
c     soil moisture averaged over the rooting depth is different for each
c     pft because of different fraction of roots present in each soil layer.
c
      do 200 j = 1, icc
        do 210 i = il1, il2
         if (fcancmx(i,j).gt.0.0) then 
         avwiltsm(i,j) =  wiltsm(i,1)*rmatctem(i,j,1) +
     &                    wiltsm(i,2)*rmatctem(i,j,2) +
     &                    wiltsm(i,3)*rmatctem(i,j,3)
         avwiltsm(i,j) = avwiltsm(i,j) /
     &    (rmatctem(i,j,1)+rmatctem(i,j,2)+rmatctem(i,j,3))
c
         afieldsm(i,j) =  fieldsm(i,1)*rmatctem(i,j,1) +
     &                    fieldsm(i,2)*rmatctem(i,j,2) +
     &                    fieldsm(i,3)*rmatctem(i,j,3)
         afieldsm(i,j) = afieldsm(i,j) /
     &    (rmatctem(i,j,1)+rmatctem(i,j,2)+rmatctem(i,j,3))
c
         avthliq(i,j)  =  thliq(i,1)*rmatctem(i,j,1) +
     &                    thliq(i,2)*rmatctem(i,j,2) +
     &                    thliq(i,3)*rmatctem(i,j,3)
         avthliq(i,j)  = avthliq(i,j) /
     &    (rmatctem(i,j,1)+rmatctem(i,j,2)+rmatctem(i,j,3))
         endif
210     continue
200   continue
c
c     Using liquid soil moisture content together with wilting and field 
c     capacity soil moisture contents averaged over the root zone, find
c     soil water status.
c
      do 230 j = 1, icc
        do 240 i = il1, il2
         if (fcancmx(i,j).gt.0.0) then 
          if(avthliq(i,j).le.avwiltsm(i,j))then
            wtstatus(i,j)=0.0
          else if(avthliq(i,j).gt.avwiltsm(i,j).and.
     &    avthliq(i,j).lt.afieldsm(i,j))then
            wtstatus(i,j)=(avthliq(i,j)-avwiltsm(i,j))/
     &      (afieldsm(i,j)-avwiltsm(i,j))
          else
            wtstatus(i,j)=1.0
          endif
         endif
240     continue
230   continue
c
c     Calculate light status as a function of lai and light extinction
c     parameter. for now set nitrogen status equal to 1, which means 
c     nitrogen is non-limiting.
c
      k1=0
      do 250 j = 1, ican
       if(j.eq.1) then
         k1 = k1 + 1
       else
         k1 = k1 + nol2pfts(j-1)
       endif
       k2 = k1 + nol2pfts(j) - 1
       do 255 m = k1, k2
        do 260 i = il1, il2
          if(j.eq.4) then  ! grasses
            ltstatus(i,m)=max(0.0, (1.0-(ailcg(i,m)/4.0)) )
          else             ! trees and crops
            ltstatus(i,m)=exp(-kn(sort(m))*ailcg(i,m))
          endif
          nstatus(i,m) =1.0
260     continue 
255    continue
250   continue
c
c     allocation to roots is determined by min. of water and nitrogen
c     status
c
      do 380 j = 1,icc
        do 390 i = il1, il2
         if (fcancmx(i,j).gt.0.0) then 
          wnstatus(i,j)=min(nstatus(i,j),wtstatus(i,j))
         endif
390     continue
380   continue
c
c     now that we know water, light, and nitrogen status we can find
c     allocation fractions for leaves, stem, and root components. note
c     that allocation formulae for grasses are different from those
c     for trees and crops, since there is no stem component in grasses. 
c
      k1=0
      do 400 j = 1, ican
       if(j.eq.1) then
         k1 = k1 + 1
       else
         k1 = k1 + nol2pfts(j-1)
       endif
       k2 = k1 + nol2pfts(j) - 1
       do 405 m = k1, k2
        do 410 i = il1, il2
          n = sort(m)
          if(j.le.3)then           !trees and crops
            denom = 1.0 + (omega(n)*( 2.0-ltstatus(i,m)-wnstatus(i,m) ))    
            afrstem(i,m)=( epsilons(n)+omega(n)*(1.0-ltstatus(i,m)) )/
     &                     denom  
            afrroot(i,m)=( epsilonr(n)+omega(n)*(1.0-wnstatus(i,m)) )/
     &                     denom  
            afrleaf(i,m)=  epsilonl(n)/denom 
          else if (j.eq.4) then     !grasses
            denom = 1.0 + (omega(n)*( 1.0+ltstatus(i,m)-wnstatus(i,m) ))
            afrleaf(i,m)=( epsilonl(n) + omega(n)*ltstatus(i,m) ) /denom  
            afrroot(i,m)=( epsilonr(n)+omega(n)*(1.0-wnstatus(i,m)) )/
     &                     denom  
            afrstem(i,m)= 0.0
          endif
410     continue
405    continue
400   continue
c
c     if using constant allocation factors then replace the dynamically
c     calculated allocation fractions.
c
      if(consallo)then
        do 420 j = 1, icc
          do 421 i = il1, il2
           if (fcancmx(i,j).gt.0.0) then 
            afrleaf(i,j)=caleaf(sort(j))
            afrstem(i,j)=castem(sort(j))
            afrroot(i,j)=caroot(sort(j))
           endif
421       continue
420     continue
      endif
c
c     make sure allocation fractions add to one
c
      do 430 j = 1, icc
        do 440 i = il1, il2 
         if (fcancmx(i,j).gt.0.0) then 
          if(abs(afrstem(i,j)+afrroot(i,j)+afrleaf(i,j)-1.0).gt.abszero) 
     &    then  
           write(6,2000) i,j,(afrstem(i,j)+afrroot(i,j)+afrleaf(i,j))
2000       format(' at (i) = (',i3,'), pft=',i2,'  allocation fractions
     &not adding to one. sum  = ',e12.7)
       write(*,*)abs(afrstem(i,j)+afrroot(i,j)+afrleaf(i,j)-1.0)-abszero
          call xit('allocate',-2)
          endif
         endif
440     continue
430   continue
c
c     the allocation fractions calculated above are overridden by two
c     rules. 
c
c     rule 1 which states that at the time of leaf onset which corresponds 
c     to leaf status equal to 1, more c is allocated to leaves so 
c     that they can grow asap. in addition when leaf status is 
c     "fall/harvest" then nothing is allocated to leaves.
c
      k1=0
      do 500 j = 1, ican
       if(j.eq.1) then
         k1 = k1 + 1
       else
         k1 = k1 + nol2pfts(j-1)
       endif
       k2 = k1 + nol2pfts(j) - 1
       do 505 m = k1, k2
        do 510 i = il1, il2
         if (fcancmx(i,m).gt.0.0) then 
          if(lfstatus(i,m).eq.1) then
            aleaf(i,m)=aldrlfon(sort(m))
c
c           for grasses we use the usual allocation even at leaf onset
c
            if(j.eq.4)then
              aleaf(i,m)=afrleaf(i,m)
            endif
c
            diff  = afrleaf(i,m)-aleaf(i,m)
            if((afrstem(i,m)+afrroot(i,m)).gt.abszero)then 
              term1 = afrstem(i,m)/(afrstem(i,m)+afrroot(i,m))
              term2 = afrroot(i,m)/(afrstem(i,m)+afrroot(i,m))
            else
              term1 = 0.0
              term2 = 0.0
            endif 
            astem(i,m) = afrstem(i,m) + diff*term1
            aroot(i,m) = afrroot(i,m) + diff*term2
            afrleaf(i,m)=aleaf(i,m)
            afrstem(i,m)=max(0.0,astem(i,m))
            afrroot(i,m)=max(0.0,aroot(i,m))
          else if(lfstatus(i,m).eq.3)then
            aleaf(i,m)=0.0
            diff  = afrleaf(i,m)-aleaf(i,m)
            if((afrstem(i,m)+afrroot(i,m)).gt.abszero)then 
              term1 = afrstem(i,m)/(afrstem(i,m)+afrroot(i,m))
              term2 = afrroot(i,m)/(afrstem(i,m)+afrroot(i,m))
            else
              term1 = 0.0
              term2 = 0.0
            endif 
            astem(i,m) = afrstem(i,m) + diff*term1
            aroot(i,m) = afrroot(i,m) + diff*term2
            afrleaf(i,m)=aleaf(i,m)
            afrstem(i,m)=astem(i,m)
            afrroot(i,m)=aroot(i,m)
          endif
         endif
510     continue
505    continue
500   continue
c
c
c     rule 2 overrides rule 1 above and makes sure that we do not allow the 
c     amount of leaves on trees and crops (i.e. pfts 1 to 7) to exceed 
c     an amount such that the remaining woody biomass cannot support. 
c     if this happens, allocation to leaves is reduced and most npp 
c     is allocated to stem and roots, in a proportion based on calculated 
c     afrstem and afrroot. for grasses this rule essentially constrains 
c     the root:shoot ratio, meaning that the model grasses can't have 
c     lots of leaves without having a reasonable amount of roots.
c
      do 530 j = 1, icc
        n=sort(j)
        do 540 i = il1, il2
         if (fcancmx(i,j).gt.0.0) then 
c         find min. stem+root biomass needed to support the green leaf 
c         biomass.
          mnstrtms(i,j)=eta(n)*(gleafmas(i,j)**kappa(n))
c
          if( (stemmass(i,j)+rootmass(i,j)).lt.mnstrtms(i,j)) then   
            if( (afrstem(i,j)+afrroot(i,j)).gt.abszero ) then
              aleaf(i,j)=min(0.05,afrleaf(i,j))
              diff  = afrleaf(i,j)-aleaf(i,j)
              term1 = afrstem(i,j)/(afrstem(i,j)+afrroot(i,j))
              term2 = afrroot(i,j)/(afrstem(i,j)+afrroot(i,j))
              astem(i,j) = afrstem(i,j) + diff*term1
              aroot(i,j) = afrroot(i,j) + diff*term2
              afrleaf(i,j)=aleaf(i,j)
              afrstem(i,j)=astem(i,j)
              afrroot(i,j)=aroot(i,j)
            else
              aleaf(i,j)=min(0.05,afrleaf(i,j))
              diff  = afrleaf(i,j)-aleaf(i,j)
              afrleaf(i,j)=aleaf(i,j)
              afrstem(i,j)=diff*0.5 + afrstem(i,j)
              afrroot(i,j)=diff*0.5 + afrroot(i,j)
            endif
          endif
         endif
540     continue
530   continue
c
c     make sure that root:shoot ratio is at least equal to rtsrmin. if not
c     allocate more to root and decrease allocation to stem.
c
      do 541 j = 1, icc
        n=sort(j)
        do 542 i = il1, il2
         if (fcancmx(i,j).gt.0.0) then 
          if( (stemmass(i,j)+gleafmas(i,j)).gt.0.05)then
            if( (rootmass(i,j)/(stemmass(i,j)+gleafmas(i,j))).
     &      lt.rtsrmin(n) ) then  
              astem(i,j)=min(0.05,afrstem(i,j))
              diff = afrstem(i,j)-astem(i,j)
              afrstem(i,j)=afrstem(i,j)-diff
              afrroot(i,j)=afrroot(i,j)+diff
            endif
          endif
         endif
542     continue
541   continue
c
c     finally check if all allocation fractions are positive and check
c     again they all add to one.
c
      do 550 j = 1, icc
        do 560 i = il1, il2
         if (fcancmx(i,j).gt.0.0) then 
          if( (afrleaf(i,j).lt.0.0).or.(afrstem(i,j).lt.0.0).or.
     &    (afrroot(i,j).lt.0.0))then
           write(6,2200) i,j
2200       format(' at (i) = (',i3,'), pft=',i2,'  allocation fractions 
     & negative') 
           write(6,2100)afrleaf(i,j),afrstem(i,j),afrroot(i,j)
2100       format(' aleaf = ',f12.9,' astem = ',f12.9,' aroot = ',f12.9)
           call xit('allocate',-3)
          endif
         endif
560     continue
550   continue
c
      do 580 j = 1, icc
        do 590 i = il1, il2
         if (fcancmx(i,j).gt.0.0) then 
          if(abs(afrstem(i,j)+afrroot(i,j)+afrleaf(i,j)-1.0).gt.abszero) 
     &    then  
           write(6,2300) i,j,(afrstem(i,j)+afrroot(i,j)+afrleaf(i,j))
2300       format(' at (i) = (',i3,'), pft=',i2,'  allocation fractions
     &not adding to one. sum  = ',f12.7)
           call xit('allocate',-4)
          endif
         endif
590     continue
580   continue
c
      return
      end

