      subroutine turnover (stemmass, rootmass,  lfstatus,    ailcg,
     1                          il1,      il2,   
     2                         sort, nol2pfts,  fcancmx,
c    3 ------------------ inputs above this line ----------------------   
     4                     stmhrlos, rothrlos,
c    5 ----------- inputs which are updated above this line -----------
     6                     stemlitr, rootlitr)
c    7 ------------------outputs above this line ----------------------
c
c               Canadian Terrestrial Ecosystem Model (CTEM) 
C                       Stem And Root Turnover Subroutine
c
c     22  Jul 2013  - Add in module for parameters
C     J. Melton
c
c     24  Sep 2012  - add in checks to prevent calculation of non-present
c     J. Melton       pfts
c
c     07  May 2003  - this subroutine calculates the litter generated
c     V. Arora        from stem and root turnover
c
c     inputs 
c
c     stemmass  - stem mass for each of the 9 ctem pfts, kg c/m2
c     rootmass  - root mass for each of the 9 ctem pfts, kg c/m2
c     lfstatus  - leaf status. an integer indicating if leaves are  
c                 in "max. growth", "normal growth", "fall/harvest", 
c                 or "no leaves" mode. see phenolgy subroutine for 
c                 more details.
c     ailcg     - green or live lai
c     icc       - no. of ctem plant function types, currently 9
c     ilg       - no. of grid cells in latitude circle
c     il1,il2   - il1=1, il2=ilg
c     sort      - index for correspondence between 9 ctem pfts and
c                 size 12 of parameter vectors
c     nol2pfts  - number of level 2 ctem pfts
c     ican        - number of class pfts
c     fcancmx   - max. fractional coverage of ctem's 9 pfts, but this can be
c                modified by land-use change, and competition between pfts
c
c     updates
c
c     stmhrlos  - stem harvest loss for crops. when in "harvest" 
c                 mode for crops, stem is also assumed to be
c                 harvested and this generates litter.
c     rothrlos  - root death for crops. when in "harvest" 
c                 mode for crops, root is assumed to die in a
c                 similar way as stem is harvested.
c
c     outputs
c
c     stemlitr  - stem litter (kg c/m2)
c     rootlitr  - root litter (kg c/m2)

      use ctem_params,        only : icc, ilg, ican, kk, zero
c
      implicit none
c
      integer il1, il2, i, j, k, lfstatus(ilg,icc), n,
     1          m,  k1,  k2
c
      integer       sort(icc),      nol2pfts(ican)
c
      real  stemmass(ilg,icc), rootmass(ilg,icc),    ailcg(ilg,icc),
     1       fcancmx(ilg,icc)
c
      real  stemlitr(ilg,icc), rootlitr(ilg,icc), nrmlsmlr(ilg,icc),
     1      nrmlrtlr(ilg,icc), rothrlos(ilg,icc), stmhrlos(ilg,icc)
c
      real       stemlife(kk),      rootlife(kk),   stmhrspn
c
c     ------------------------------------------------------------------
c
c     parameters used. also note the structure of parameter vectors
c     which clearly shows the class pfts (along rows) and ctem
c     sub-pfts (along columns)
c
c     needle leaf |  evg       dcd       ---
c     broad leaf  |  evg   dcd-cld   dcd-dry
c     crops       |   c3        c4       ---
c     grasses     |   c3        c4       ---
c
c     stemlife, turnover time scale for stem for different pfts
c       flag testing new parameter values needle evg jm 07.06.2012
c      data  stemlife/60.0, 75.0, 0.00,
      data  stemlife/65.0, 75.0, 0.00,
     &               45.0, 70.0, 65.0,  !flag pft 4 was 40.0, pft 5 was 45.0
     &               20.0, 20.0, 0.00,
     &               0.00, 0.00, 0.00/
c
c     rootlife, turnover time scale for root for different pfts
c       flag testing new parameter values needle evg jm 07.06.2012
c      data  rootlife/8.5,11.5, 0.0,
      data  rootlife/10.0,11.5, 0.0,
     &               8.5, 9.5, 8.5,     ! flag pft 4 and 5 were 5.5
     &               3.0, 3.0, 0.0,
     &               2.5, 2.5, 0.0/
c
c     stem harvest span. same as crop harvest span. period in days
c     over which crops are harvested. 
      data stmhrspn/17.0/
c
c     ---------------------------------------------------------------
c
c     initialize required arrays to zero
c
      do 140 j = 1,icc
        do 150 i = il1, il2
          stemlitr(i,j)=0.0          ! total stem litter
          rootlitr(i,j)=0.0          ! total root litter
          nrmlsmlr(i,j)=0.0          ! stem litter from normal turnover
          nrmlrtlr(i,j)=0.0          ! root litter from normal turnover
150     continue                  
140   continue
c
c     initialization ends    
c
c     ------------------------------------------------------------------
c
c     calculate normal stem and root litter using the amount of stem and
c     root biomass and their turnover time scales.
c
      do 200 j = 1, icc
       n = sort(j)
       do 210 i = il1, il2
       if (fcancmx(i,j).gt.0.0) then
        if(stemlife(n).gt.zero)then
         nrmlsmlr(i,j)=stemmass(i,j)*(1.0-exp(-1.0/(365.0*stemlife(n))))  
        endif
        if(rootlife(n).gt.zero)then
         nrmlrtlr(i,j)=rootmass(i,j)*(1.0-exp(-1.0/(365.0*rootlife(n))))  
        endif
       endif
210    continue
200   continue
c
c     if crops are in harvest mode then we start harvesting stem as well.
c     if stem has already been harvested then we set the stem harvest
c     loss equal to zero. the roots of the crop die in a similar way.
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
         if (fcancmx(i,m).gt.0.0) then
          if(j.eq.3)then     !stem/root harvest/death for crops
c
            if(lfstatus(i,m).eq.3.and.stmhrlos(i,m).le.zero.and.
     &      stemmass(i,m).gt.zero)then          
              stmhrlos(i,m)=stemmass(i,m)*(1.0/stmhrspn)
            endif
c
            if(lfstatus(i,m).eq.3.and.rothrlos(i,m).le.zero.and.
     &      rootmass(i,m).gt.zero)then          
              rothrlos(i,m)=rootmass(i,m)*(1.0/stmhrspn)   
            endif
c
            if(stemmass(i,m).le.zero.or.lfstatus(i,m).eq.1.or.
     &      lfstatus(i,m).eq.2)then
              stmhrlos(i,m)=0.0
            endif
c
            if(rootmass(i,m).le.zero.or.lfstatus(i,m).eq.1.or.
     &      lfstatus(i,m).eq.2)then
              rothrlos(i,m)=0.0
            endif
c
          else
            stmhrlos(i,m)=0.0
            rothrlos(i,m)=0.0
          endif
         endif
260     continue
255    continue   
250   continue   
c
c     add stem and root litter from all sources
c
      do 350 j = 1, icc
        do 360 i = il1, il2
         if (fcancmx(i,j).gt.0.0) then
          stemlitr(i,j)=nrmlsmlr(i,j)+stmhrlos(i,j)
          rootlitr(i,j)=nrmlrtlr(i,j)+rothrlos(i,j)
         endif
360     continue
350   continue
c
      return
      end

