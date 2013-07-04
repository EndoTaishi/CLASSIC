      subroutine    bio2str( gleafmas, bleafmas, stemmass, rootmass,                        
     1                            il1,      il2,  fcancmx,    zbotw,
     2                          delzw, nol2pfts,  soildpth,
c    4--------------- inputs above this line, outputs below --------
     5                          ailcg,    ailcb,     ailc,    zolnc,
     6                          rmatc, rmatctem,     slai,  bmasveg,
     7                       cmasvegc,  veghght, rootdpth,   alvisc,
     8                         alnirc,     paic,    slaic )
c
c     ----------------------------------------------------------------
c
C           Canadian Terrestrial Ecosystem Model (CTEM)
C        Biomass To Structural Attributes Conversion Subroutine 
c

c     2   Jul 2013  - Integreated ctem_params module
c     J. Melton       

c     22  Nov 2012  - calling this version 1.1 since a fair bit of ctem
c     V. Arora        subroutines were changed for compatibility with class
c                     version 3.6 including the capability to run ctem in
c                     mosaic/tile version along with class.
c
c     24  Sep 2012  - add in checks to prevent calculation of non-present
c     J. Melton       pfts. cleaned up initialization
c
c     18  May  2012 - fix bug for soils with the depth close to the
c     J. Melton       boundary between layers, was resulting in roots 
c                     being placed incorrectly.
c
c     28  Nov. 2011 - make changes for coupling with class 3.5
c     Y. Peng
c
c     31  Aug. 2009 - make changes for coupling with class 3.4. include
c     V. Arora        new output variable called plant area index
c                     (paic) and storage leaf area index (slaic) 
c
c     14  Mar. 2003 - this subroutine converts leaf, stem, and root biomass 
c     V. Arora        into lai, vegetation height, and fraction of roots
c                     in each soil layer. storage lai is also calculated.
c                     
c                     note that while ctem keeps track of 9 pfts, class 2.7
c                     keeps track of 4 pfts, so all these vegetation 
c                     structural attributes are calculated for 9 pfts
c                     sepatarely but then lumped into 4 pfts for use in
c                     energy & water balance calculations by class
c               
c                     also, this subroutine does not estimate zolnc(i,5)
c                     the log of roughness length over the bare fraction
c                     of the grid cell. only roughness lengths over the
c                     vegetated fraction are updated
c 
c              ctem 1.0                 class 2.7
c         -----------------            -------------
c        1. needle leaf evg
c        2. needle leaf dcd           1. needle leaf
c        3. broad leaf evg
c        4. broad leaf cold dcd       2. broad leaf
c        5. broad leaf dry dcd            
c        6. c3 crop
c        7. c4 crop                   3. crop
c        8. c3 grass
c        9. c4 grass                  4. grass
c
c     ----------------------------------------------------------------
c     inputs
c
c     gleafmas  - green or live leaf mass in kg c/m2, for the 9 pfts
c     bleafmas  - brown or dead leaf mass in kg c/m2, for the 9 pfts
c     stemmass  - stem biomass in kg c/m2, for the 9 pfts
c     rootmass  - root biomass in kg c/m2, for the 9 pfts
c     fcancmx   - max. fractional coverages of ctem's 9 pfts. this is
c                 different from fcanc and fcancs (which may vary with
c                 snow depth). fcancmx doesn't change, unless of course
c                 its changed by land use change or dynamic vegetation.
c     delzw     - thicknesses of the 3 soil layers
c     zbotw     - bottom of soil layers
c     il1, il2  - il1=1, il2=ilg
c     nol2pfts  - number of level 2 pfts
c     soildpth  - soil depth (m)
c    
c     outputs
c
c     ailcg     - green lai for ctem's 9 pfts
c     ailcb     - brown lai for ctem's 9 pfts. for now we assume only
c                 grasses can have brown leaves.
c     ailc      - lai to be used by class
c     zolnc     - log of roughness length to be used by class
c     rmatctem  - fraction of live roots in each soil layer for each of 
c                 ctem's 9 pfts
c     rmatc     - fraction of live roots in each soil layer for each of the
c                 class' 4 pfts
c     slai      - storage or imaginary lai for phenology purposes
c     bmasveg   - total (gleaf + stem + root) biomass for each ctem pft, kg c/m2
c     cmasvegc  - total canopy mass for each of the 4 class pfts. recall that 
c                 class requires canopy mass as an input, and this is now 
c                 provided by ctem. kg/m2.
c     veghght   - vegetation height (meters)
c     rootdpth  - 99% soil rooting depth (meters)
c                 both veghght & rootdpth can be used as diagnostics to see
c                 how vegetation grows above and below ground, respectively.
c     alvisc    - albedo for 4 class pfts simulated by ctem, visible 
c     alnirc    - albedo for 4 class pfts simulated by ctem, near ir
c     paic      - plant area index for class' 4 pfts. this is the sum of
c                 leaf area index and stem area index.
c     slaic     - storage lai. this will be used as min. lai that class
c                 sees so that it doesn't blow up in its stomatal
c                 conductance calculations.
c 
c     ---------------------------------------------------------------- 

      use ctem_params,        only : ignd, icc, ilg, ican, zero, l2max,
     1                               kk, eta, kappa, kn, lfespany, 
     2                               fracbofg, specsla     
c
      implicit none

      integer il1, il2, i, j, k, m, n, k1c, k2c
      integer nol2pfts(ican), icount, sort(icc),kend

      logical deeproots
        
      real  gleafmas(ilg,icc),  bleafmas(ilg,icc),  stemmass(ilg,icc),
     1      rootmass(ilg,icc),     ailcg(ilg,icc),     ailcb(ilg,icc),   
     2           ailc(ilg,ican),      zolnc(ilg,ican), paic(ilg, ican),
     3       rmatc(ilg,ican,ignd),   fcancmx(ilg,icc), delzw(ilg,ignd),
     4          zbotw(ilg,ignd),               rmatctem(ilg,icc,ignd),
     5          slai(ilg,icc),   bmasveg(ilg,icc),   cmasvegc(ilg,ican),
     6           sai(ilg,icc),       saic(ilg,ican), sfcancmx(ilg,ican),
     7         alvisc(ilg,ican),     alnirc(ilg,ican),  pai(ilg,icc),
     8          slaic(ilg,ican)

      real           sla(icc),   veghght(ilg,icc),
     1      lnrghlth(ilg,icc),   averough(ilg,ican),           abar(kk),
     2           avertmas(kk),          alpha(kk),             b(icc),
     3      rootdpth(ilg,icc),  usealpha(ilg,icc),         a(ilg,icc),
     4          useb(ilg,icc),              zroot,      soildpth(ilg),
     5           prcnslai(kk),         minslai(kk),        fcoeff,            
     7           mxrtdpth(kk),         albvis(kk),         albnir(kk),
     8       etmp(ilg,icc,ignd),    totala(ilg,icc),       rmat_sum
c
c     ---------------------------------------------------------------
c     parameters used. also note the structure of parameter vectors
c     which clearly shows the class pfts (along rows) and ctem sub-pfts
c     (along columns) 
c
c     needle leaf |  evg       dcd       ---    
c     broad leaf  |  evg   dcd-cld   dcd-dry   
c     crops       |   c3        c4       ---  
c     grasses     |   c3        c4       --- 
c
c     parameter determining average root profile
      data abar/4.70, 5.86, 0.00,
     &          3.87, 3.46, 3.97,
     &          3.97, 3.97, 0.00,
     &          5.86, 4.92, 0.00/
c
c     average root biomass (kg c/m2) for ctem's 8 pfts used for
c     estimating rooting profile
      data avertmas/1.85, 1.45, 0.00,
     &              2.45, 2.10, 2.10,
     &              0.10, 0.10, 0.00,
     &              0.70, 0.70, 0.00/
c
c     alpha, parameter determining how the roots grow
      data alpha/0.80, 0.80, 0.00,
     &           0.80, 0.80, 0.80,
     &           0.80, 0.80, 0.00,
     &           0.80, 0.80, 0.00/
c
c     prcnslai - storage/imaginary lai is this percentage of maximum 
c     leaf area index that a given root+stem biomass can supoort 
      data prcnslai/7.5, 7.5, 0.0,
     &              7.5, 7.5, 7.5, 
     &              7.5, 7.5, 0.0,
     &              2.5, 2.5, 0.0/
c
c     minslai - minimum storage lai. this is what the model uses as
c     lai when growing vegetation for scratch. consider these as model
c     seeds.
      data minslai/0.3, 0.3, 0.0,
     &             0.3, 0.3, 0.3,
     &             0.2, 0.2, 0.0,
     &             0.2, 0.2, 0.0/
c
c      maximum rooting depth. this is used so that the rooting depths
c      simulated by ctem's variable rooting depth parameterzation are
c      constrained to realistic values
c
       data mxrtdpth/3.00, 3.00, 0.00,
     &               5.00, 5.00, 3.00,
     &               2.00, 2.00, 0.00,
     &               1.00, 1.00, 0.00/
c
c     visible and near ir albedos of the 9 ctem pfts
c
       data albvis/3.00, 3.00, 0.00,
     &             3.00, 5.00, 5.00,
     &             5.50, 5.50, 0.00,
     &             5.00, 6.00, 0.00/
c
       data albnir/19.0, 19.0, 0.00,
     &             23.0, 29.0, 29.0,
     &             34.0, 34.0, 0.00,
     &             30.0, 34.0, 0.00/
c
c     class' original root parameterization has deeper roots than ctem's
c     default values based on literature. in the coupled model this leads
c     to lower evapotranspiration (et) values. an option is provided here to
c     deepen roots, but this will also increase photosynthesis and vegetation
c     biomass slightly, due to more access to soil water. so while use of
c     deeper roots is desirable in the coupled global model, one may decide
c     to use ctem's default parameterizarion for stand alone simulations, and
c     set deeproots to .false. 
c
      data deeproots/.false./
c     ---------------------------------------------------------------
c
c     initialization
c
      do 30 j = 1,icc
        do 40 k = 1,ignd
          do 50 i = il1,il2
            rmatctem(i,j,k)=0.0
            etmp(i,j,k)    =0.0
50        continue
40      continue
30    continue
c
      do 31 j = 1,ican
        do 41 k = 1,ignd
          do 51 i = il1,il2
            rmatc(i,j,k)=0.0
51        continue
41      continue
31    continue
c
      icount=0
      do 52 j = 1, ican
        do 53 m = 1, nol2pfts(j)
          n = (j-1)*l2max + m
          icount = icount + 1
          sort(icount)=n
53      continue
52    continue
c
      do 60 j = 1,ican
        do 70 i =  il1, il2
          ailc(i,j)=0.0
          saic(i,j)=0.0
          paic(i,j)=0.0
          slaic(i,j)=0.0
          zolnc(i,j)=0.0
          averough(i,j)=0.0
          alvisc(i,j)=0.0
          alnirc(i,j)=0.0
          cmasvegc(i,j)=0.0
          sfcancmx(i,j)=0.0    ! sum of fcancmxs

70      continue
60    continue
c
      do 80 j = 1,icc
        sla(j)=0.0
        do 90 i =  il1, il2
          usealpha(i,j)=alpha(sort(j))
          useb(i,j)=0.0
          ailcg(i,j)=0.0
          ailcb(i,j)=0.0
          veghght(i,j)=0.0
          lnrghlth(i,j)=0.0
          a(i,j)=0.0
          slai(i,j)=0.0
          sai(i,j)=0.0
          bmasveg(i,j)=0.0
          pai(i,j)=0.0
90      continue
80    continue
c
c
c     ------ 1. conversion of leaf biomass into leaf area index -------
c
c     find specific leaf area (sla, m2/kg) using leaf life span
c
      icount=0
      do 100 j = 1, ican
        do 101 m = 1, nol2pfts(j)
          n = (j-1)*l2max + m
          icount = icount + 1
          sla(icount) = 25.0*(lfespany(n)**(-0.50))
          if(specsla(n).gt.zero) sla(icount)=specsla(n)  
101     continue
100   continue
c
c     convert leaf biomass into lai. brown leaves could have less
c     lai than the green leaves for the same leaf mass. for now we
c     assume sla of brown leaves is fracbofg times that of green 
c     leaves. 
c
c     also find stem area index as a function of stem biomass
c
      do 150 j = 1,icc
        do 160 i = il1,il2
         if (fcancmx(i,j).gt.0.0) then
          ailcg(i,j)=sla(j)*gleafmas(i,j)
          ailcb(i,j)=sla(j)*bleafmas(i,j)*fracbofg
          sai(i,j)=0.55*(1.0-exp(-0.175*stemmass(i,j))) !stem area index
c         plant area index is sum of green and brown leaf area indices
c         and stem area index
          pai(i,j)=ailcg(i,j)+ailcb(i,j)+sai(i,j)

c         make class see some minimum pai, otherwise it runs into numerical
c         problems
            pai(i,j)=max(0.3,pai(i,j))

         endif
160     continue
150   continue
c
c     get fcancmx weighted leaf area index for use by class
c       needle leaf evg + dcd = total needle leaf    
c       broad leaf evg + dcd cld + dcd dry = total broad leaf    
c       crop c3 + c4 = total crop
c       grass c3 + c4 = total grass
c     also add brown lai. note that although green + brown
c     lai is to be used by class for energy and water balance
c     calculations, stomatal conductance estimated by the 
c     photosynthesis subroutine is only based on green lai.
c     that is although both green+brown leaves intercept
c     water and light, only the green portion photosynthesizes.
c     also lump stem and plant area indices for class' 4 pfts
c
      k1c=0
      do 200 j = 1, ican
        if(j.eq.1) then
          k1c = k1c + 1
        else
          k1c = k1c + nol2pfts(j-1)
        endif
        k2c = k1c + nol2pfts(j) - 1
        do 210 m = k1c, k2c
          do 220 i = il1, il2
            sfcancmx(i,j)=sfcancmx(i,j)+fcancmx(i,m)
            ailc(i,j)=ailc(i,j)+(fcancmx(i,m)*(ailcg(i,m)+ailcb(i,m)))
            saic(i,j)=saic(i,j)+(fcancmx(i,m)*sai(i,m))
            paic(i,j)=paic(i,j)+(fcancmx(i,m)*pai(i,m))
            slaic(i,j)=slaic(i,j)+(fcancmx(i,m)*slai(i,m))
220       continue
210     continue
200   continue
c
      do 230 j = 1, ican
        do 240 i = il1, il2
c
          if(sfcancmx(i,j).gt.zero)then
             ailc(i,j)=ailc(i,j)/sfcancmx(i,j)
             saic(i,j)=saic(i,j)/sfcancmx(i,j)
             paic(i,j)=paic(i,j)/sfcancmx(i,j)
             slaic(i,j)=slaic(i,j)/sfcancmx(i,j)
          else
             ailc(i,j)=0.0
             saic(i,j)=0.0
             paic(i,j)=0.0
             slaic(i,j)=0.0
          endif
c         for crops and grasses set the minimum lai to a small number, other
c         wise class will never run tsolvc and thus phtsyn and ctem will not
c         be able to grow crops or grasses.
          if(j.eq.3.or.j.eq.4) ailc(i,j)=max(ailc(i,j),0.1)
c
240     continue
230   continue
c
c     ------ 2. conversion of stem biomass into roughness length -------
c 
c     class uses log of roughness length (zoln) as an input parameter. when 
c     vegetation grows and dies as per ctem, then zoln is provided by ctem.
c
c     1. convert stem biomass into vegetation height for trees and crops,
c        and convert leaf biomass into vegetation height for grass
c     2. convert vegetation height into roughness length & take its log
c     3. lump this for ctem's 9 pfts into class' 4 pfts
c
      k1c=0
      do 250 j = 1, ican
        if(j.eq.1) then
          k1c = k1c + 1
        else
          k1c = k1c + nol2pfts(j-1)
        endif
        k2c = k1c + nol2pfts(j) - 1
        do 260 m = k1c, k2c
          do 270 i = il1, il2
c          
          if (j.le.2) then                            ! trees
           veghght(i,m)=10.0*stemmass(i,m)**0.385
           veghght(i,m)=min(veghght(i,m),45.0)
          else if (j.eq.3) then                       ! crops
           veghght(i,m)=1.0*(stemmass(i,m)+gleafmas(i,m))**0.385
          else if (j.eq.4) then                       ! grasses
           veghght(i,m)=3.5*(gleafmas(i,m)+fracbofg*bleafmas(i,m))**0.50   
          endif
          lnrghlth(i,m)= log(0.10 * max(veghght(i,m),0.10))
c
270       continue
260     continue
250   continue
c

      k1c=0
      do 300 j = 1, ican
        if(j.eq.1) then
          k1c = k1c + 1
        else
          k1c = k1c + nol2pfts(j-1)
        endif
        k2c = k1c + nol2pfts(j) - 1
        do 310 m = k1c, k2c
          do 320 i = il1, il2
            averough(i,j)=averough(i,j)+(fcancmx(i,m)*lnrghlth(i,m))
320       continue
310     continue
300   continue
c
      do 330 j = 1, ican
        do 340 i = il1, il2
c
          if(sfcancmx(i,j).gt.zero)then
             averough(i,j)=averough(i,j)/sfcancmx(i,j)
          else
            averough(i,j)=-4.605
          endif
          zolnc(i,j)=averough(i,j)
c
340     continue
330   continue
c
c
c     ------ 3. estimating fraction of roots in each soil layer for -----
c     ------      ctem's each vegetation type, using root biomass   -----
c
c     
c     estimate parameter b of variable root profile parameterization
c
      icount=0
      do 350 j = 1, ican
        do 360 m = 1, nol2pfts(j)
          n = (j-1)*l2max + m
          icount = icount + 1
c         use decreased abar for all pfts to deepen roots to account
c         for hydraulic redistrubution
          if(deeproots) then
              b(icount) = (abar(n)-1.5) * (avertmas(n)**alpha(n))
          else
              b(icount) = abar(n) * (avertmas(n)**alpha(n))
          endif
360     continue
350   continue
c
c     use b to estimate 99% rooting depth
c
      k1c=0
      do 370 j = 1,ican
        if(j.eq.1) then
          k1c = k1c + 1
        else
          k1c = k1c + nol2pfts(j-1)
        endif
        k2c = k1c + nol2pfts(j) - 1
        do 380 m = k1c, k2c
          do 390 i = il1, il2
c
            useb(i,m)=b(m)
            usealpha(i,m)=alpha(sort(m))
            rootdpth(i,m) = (4.605*(rootmass(i,m)**alpha(sort(m))))/b(m)

c
c           if estimated rooting depth is greater than soil depth, or
c           the maximum rooting depth then adjust rooting depth and
c           parameter alpha
c
c           also find "a" (parameter determining root profile). this is 
c           the "a" which depends on time varying root biomass 
c
c
            if(rootdpth(i,m).gt.min(soildpth(i),zbotw(i,ignd),
     1                                mxrtdpth(sort(m))))then
              rootdpth(i,m) = min(soildpth(i),zbotw(i,ignd),
     1                                  mxrtdpth(sort(m)))
              if(rootdpth(i,m).le.zero)then
                a(i,m)=100.0
              else
                a(i,m)=4.605/rootdpth(i,m)
              endif
            else
              if(rootmass(i,m).le.zero)then
                a(i,m)=100.0
              else
                a(i,m)=useb(i,m)/(rootmass(i,m)**usealpha(i,m))
              endif
            endif
c
390       continue
380     continue
370   continue
c
      do 400 j = 1,icc
        do 410 i = il1, il2

           kend=9999  ! initialize with a dummy value
c
c          using parameter "a" we can find fraction of roots in each soil layer
c          just like class
c
           zroot=rootdpth(i,j)

           totala(i,j) = 1.0-exp(-a(i,j)*zroot)

           if(zroot.le.zbotw(i,1))then
c           if rootdepth is shallower than the bottom of the first layer
            rmatctem(i,j,1)=1.0
            do 414 k=2,ignd
             rmatctem(i,j,k)=0.0
414         continue
            kend=1
           else
c
            do 415 k=2,ignd
             if(zroot.le.zbotw(i,k).and.zroot.gt.zbotw(i,k-1))then
c             if rootdepth is shallower than the bottom of current layer and
c             is deeper than bottom of the previous top layer
              kend=k ! kend = soil layer number in which the roots end  
             endif
415         continue
c
            if (kend .eq. 9999) then
              write(6,2100) i,j,k,kend
2100          format(' at (i) = (',i3,'), pft=',i2,', depth=',i2,' kend 
     & is not assigned. kend  = ',i5)
              call xit('bio2str',-3)
            end if

            etmp(i,j,1)=exp(-a(i,j)*zbotw(i,1))
            rmatctem(i,j,1)=(1.0-etmp(i,j,1))/totala(i,j)
            if (kend .eq. 2) then
c             if rootdepth is shallower than the bottom of 2nd layer
               etmp(i,j,kend)=exp(-a(i,j)*zroot)
               rmatctem(i,j,kend)=(etmp(i,j,kend-1)-etmp(i,j,kend))
     1                          /totala(i,j)
            elseif (kend .gt. 2) then
c             if rootdepth is shallower than the bottom of 3rd layer 
c             or even the deeper layer (ignd>3)

              do 416 k=2,kend-1
                etmp(i,j,k)=exp(-a(i,j)*zbotw(i,k))
                rmatctem(i,j,k)=(etmp(i,j,k-1)-etmp(i,j,k))/totala(i,j)
416           continue

              etmp(i,j,kend)=exp(-a(i,j)*zroot)
              rmatctem(i,j,kend)=(etmp(i,j,kend-1)-etmp(i,j,kend))
     1                          /totala(i,j)
            endif
           endif
c
410     continue
400   continue
c
c    make sure all fractions (of roots in each layer) add to one.
c
      do 411 j = 1, icc
        do 412 i = il1, il2
         if (fcancmx(i,j).gt.0.0) then
          rmat_sum = 0.0
          do 413 k = 1, ignd
            rmat_sum = rmat_sum + rmatctem(i,j,k)
413       continue
c
          if( abs(rmat_sum-1.0).gt.1e-10) then
           write(6,2300) i,j,rmat_sum
2300       format(' at (i) = (',i3,'), pft=',i2,' fractions of roots
     &not adding to one. sum  = ',f12.7)
           call xit('bio2str',-3)
          endif
         endif
412     continue
411   continue
c
c     lump rmatctem(i,9,ignd)  into rmatc(i,4,ignd) for use by class
c
      k1c=0
      do 420 j = 1, ican
        if(j.eq.1) then
          k1c = k1c + 1
        else
          k1c = k1c + nol2pfts(j-1)
        endif
        k2c = k1c + nol2pfts(j) - 1
        do 430 m = k1c, k2c
          do 440 i = il1, il2
c   
            do 441 k = 1, ignd
              rmatc(i,j,k)=rmatc(i,j,k)+(fcancmx(i,m)*rmatctem(i,m,k))  
441         continue
c
440       continue
430     continue
420   continue
c
      do 450 j = 1, ican
        do 460 i = il1, il2
c
          if(sfcancmx(i,j).gt.zero)then
             do 461 k = 1, ignd
               rmatc(i,j,k)=rmatc(i,j,k)/sfcancmx(i,j)
461          continue
          else
             rmatc(i,j,1)=1.0
             do 462 k = 2, ignd
               rmatc(i,j,k)=0.0
462          continue
          endif
c
460     continue
450   continue
c
c
c     -------------------  4. calculate storage lai  --------------------
c
      do 500 j = 1, icc
        do 510 i = il1, il2
         if (fcancmx(i,j).gt.0.0) then
          slai(i,j)=((stemmass(i,j)+rootmass(i,j))/eta(sort(j)))
     &     **(1./kappa(sort(j)))
          slai(i,j)=(prcnslai(sort(j))/100.0)*sla(j)*slai(i,j)
c
c         need a minimum slai to be able to grow from scratch. consider
c         this as model seeds.
          slai(i,j)=max(slai(i,j),minslai(sort(j)))
         endif
510     continue
500   continue
c
c
c     --- 5. calculate total vegetation biomass for each ctem pft, and --
c     ---------------- canopy mass for each class pft ------------------
c
      do 550 j = 1, icc
        do 560 i = il1, il2
         if (fcancmx(i,j).gt.0.0) then
          bmasveg(i,j)=gleafmas(i,j)+stemmass(i,j)+rootmass(i,j)
         endif
560     continue
550   continue
c
c     since class uses canopy mass and not total vegetation biomass as an
c     input, we find canopy mass as a sum of stem and leaf mass, for each
c     class pft, i.e. only above ground biomass. 
c
      k1c=0
      do 600 j = 1, ican
        if(j.eq.1) then
          k1c = k1c + 1
        else
          k1c = k1c + nol2pfts(j-1)
        endif
        k2c = k1c + nol2pfts(j) - 1
        do 610 m = k1c, k2c
          do 620 i = il1, il2
            cmasvegc(i,j)= cmasvegc(i,j) +
     &      (fcancmx(i,m)*(bleafmas(i,m)+gleafmas(i,m)+stemmass(i,m)))  
620       continue
610     continue
600   continue
c
      do 630 j = 1, ican
        do 640 i = il1, il2
c
          if(sfcancmx(i,j).gt.zero)then
            cmasvegc(i,j)=cmasvegc(i,j)/sfcancmx(i,j)
            cmasvegc(i,j)=cmasvegc(i,j)*(1.0/0.50) !assuming biomass is 50% c
          else
            cmasvegc(i,j)=0.0
          endif
c       
c         if there is no vegetation canopy mass will be zero. this should 
c         essentially mean more bare ground, but since we are not changing
c         fractional coverages at present, we pass a minimum canopy mass
c         to class so that it doesn't run into numerical problems.
c
          cmasvegc(i,j)=max(cmasvegc(i,j),3.0)
c
640     continue
630   continue
c
c     --- 6. calculate albedo for class' 4 pfts based on specified ----
c     ------ albedos of ctem 9 pfts and their fractional coveraes -----
c
      k1c=0
      do 700 j = 1, ican
        if(j.eq.1) then
          k1c = k1c + 1
        else
          k1c = k1c + nol2pfts(j-1)
        endif
        k2c = k1c + nol2pfts(j) - 1
        do 710 m = k1c, k2c
          do 720 i = il1, il2
            alvisc(i,j)= alvisc(i,j) + (fcancmx(i,m)*albvis(sort(m)))  
            alnirc(i,j)= alnirc(i,j) + (fcancmx(i,m)*albnir(sort(m)))  
720       continue
710     continue
700   continue
c
      do 730 j = 1, ican
        do 740 i = il1, il2

          if(sfcancmx(i,j).gt.zero)then
            alvisc(i,j)=(alvisc(i,j)/sfcancmx(i,j))/100.0
            alnirc(i,j)=(alnirc(i,j)/sfcancmx(i,j))/100.0
          else
            alvisc(i,j)=0.0
            alnirc(i,j)=0.0
          endif

740     continue
730   continue

      return
      end


