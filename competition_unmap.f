      subroutine competition_unmap( nml, ilmos, jlmos, nol2pfts,
     2                             fare_cmp,  nppveg_cmp,geremort_cmp,
     3                         intrmort_cmp,gleafmas_cmp,bleafmas_cmp,
     4                         stemmass_cmp,rootmass_cmp,litrmass_cmp,
     5                         soilcmas_cmp,pftexist_cmp,  lambda_cmp,
     6                       bmasveg_cmp,burnvegf_cmp,add2allo_cmp,
     7                               cc_cmp,      mm_cmp,  fcanmx_cmp,
     8                         vgbiomas_cmp, grclarea_cmp,
     9                         gavgltms_cmp,gavgscms_cmp,
     1                               ta_cmp,  precip_cmp,  netrad_cmp, 
     2                            tcurm_cmp,srpcuryr_cmp,dftcuryr_cmp,
     3                           tmonth_cmp,anpcpcur_cmp, anpecur_cmp, 
     4                          gdd5cur_cmp,surmncur_cmp,defmncur_cmp,
     5                         srplscur_cmp,defctcur_cmp,  twarmm_cmp, 
     6                           tcoldm_cmp,    gdd5_cmp, aridity_cmp,
     7                         srplsmon_cmp,defctmon_cmp,anndefct_cmp,
     8                         annsrpls_cmp,  annpcp_cmp,
     &                         dry_season_length_cmp,
     9                         lucemcom_cmp, lucltrin_cmp, lucsocin_cmp,
     1                         pfcancmx_cmp, nfcancmx_cmp,pstemmass_cmp,
     2                         pgleafmass_cmp,
c    ------------------- inputs above this line ----------------------
     a                            netradrow,
c    ------------------- saved for intermediate above this line ------
     a                              faregat, fcancmx,  nppveg,geremort,  
     b                             intrmort,gleafmas,bleafmas,stemmass,
     c                             rootmass,litrmass,soilcmas,grclarea,
     d                             pftexist,lambda,bmasveg,burnvegf,
     e                             add2allo,      cc,      mm,  fcanmx,
     f                             vgbiomas,gavgltms,gavgscms,
     g                             ta,precip, netrad,   tcurm,srpcuryr,
     h                             dftcuryr,  tmonth,anpcpcur, anpecur, 
     i                              gdd5cur,surmncur,defmncur,srplscur,  
     j                             defctcur,  twarmm,  tcoldm,    gdd5, 
     k                              aridity,srplsmon,defctmon,anndefct,
     l                             annsrpls,  annpcp,
     &                         dry_season_length,
     m                           lucemcom,  lucltrin, lucsocin,
     n                            pfcancmx, nfcancmx, pstemmass, 
     o                            pgleafmass )
c    ------------------- updates above this line ---------------------
c
c
C               Canadian Terrestrial Ecosystem Model (CTEM)
C                       Unmapping For Competition 
c
c     23  Jul 2013  - Add in module for parameters
C     J. Melton

c    28 Aug 2012 - this subroutine unmap the pft fractions in 
c    Y. Peng       each grid cell back to the only pft in each 
c                  mosaic after the competition calculations are done
c
c
c                  input:        array_cmp(nlat,icc)
c                    | 
c                    | unmapping      
c                    v
c                  intermediate: arrayrow(nlat,nmos,icc)
c                    | 
c                    | gathering
c                    v
c                  output(updates): array(ilg,icc)
c
c                  array_cmp: mapped array used for competition
c
c    -----------------------------------------------------------------
c
c    indices
c
c     nlat    - max. number of grid cells in the latitude circle, which
c               is prescribed in runclassXXctem.f
c     nmos    - max. number of mosaic tiles in each latitudinal grid cell, 
c               which is prescribed in runclassXXctem.f
c     ilg     - ilg=nlat*nmos
c     nml     - total number of mosaic tiles with pft fractions larger than 1,
c               see gatprep.f
c     ilmos   - indices for scattering, see gatprep.f
c     jlmos   - indices for scattering, see gatprep.f
c     icc     - number of pfts for use by ctem, currently 10
c     ican      - number of pfts for use by class, currently 4
c
c    inputs
c
c     fare_cmp    - fractional coverage of ctem's 10 pfts in each 
c                   latitudinal grid cell
c     nppveg_cmp  - npp for each pft type of vegetated area in each 
c                   latitudinal grid cell
c     geremort_cmp- growth related mortality in each latitudinal grid cell
c     intrmort_cmp- intrinsic (age related) mortality in each latitudinal grid cell
c     gleafmas_cmp- green leaf mass for each of the 10 ctem pfts in each 
c                   latitudinal grid cell
c     bleafmas_cmp- brown leaf mass for each of the 10 ctem pfts in each 
c                   latitudinal grid cell
c     stemmass_cmp- stem mass for each of the 10 ctem pfts in each 
c                   latitudinal grid cell
c     rootmass_cmp- root mass for each of the 10 ctem pfts in each 
c                   latitudinal grid cell
c     litrmass_cmp- litter mass for each of the 10 ctem pfts + bare, in each 
c                   latitudinal grid cell
c     soilcmas_cmp- soil carbon mass for each of the 10 ctem pfts + bare, in each 
c                   latitudinal grid cell 
c     pftexist_cmp- logical array indicating pfts exist (t) or not (f) in each 
c                   latitudinal grid cell 
c     lambda_cmp  - fraction of npp that is used for horizontal expansion in each 
c                   latitudinal grid cell 
c     bmasveg_cmp - total (gleaf + stem + root) biomass for each ctem pft, kg c/m2 in each 
c                   latitudinal grid cell 
c     burnvegf_cmp - fractional areas burned, for ctem pfts in each 
c                   latitudinal grid cell 
c     add2allo_cmp- npp kg c/m2.day in each latitudinal grid cell that is used 
c                   for expansion and subsequently allocated to leaves, stem,  
c                   and root via the allocation part of the model.
c     cc,mm_cmp   - colonization rate & mortality rate in each 
c                   latitudinal grid cell  
c     fcanmx_cmp  - fractional coverage of class' 4 pfts in each 
c                   latitudinal grid cell 
c     vgbiomas_cmp- grid averaged vegetation biomass, kg c/m2
c     gavgltms_cmp- grid averaged litter mass, kg c/m2
c     gavgscms_cmp- grid averaged soil c mass, kg c/m2
c     grclarea_cmp- area of the grid cell, km^2
c
c     ta_cmp      - mean daily temperature (k) in each latitudinal grid cell 
c     precip_cmp  - daily precipitation (mm/day) in each latitudinal grid cell 
c     netrad_cmp  - daily net radiation (w/m2) in each latitudinal grid cell 
c     tcurm_cmp   - temperature of the current month (c) in each latitudinal grid cell 
c     srpcuryr_cmp- water surplus for the current year in each latitudinal grid cell 
c     dftcuryr_cmp- water deficit for the current year in each latitudinal grid cell 
c     tmonth_cmp  - monthly temperatures in each latitudinal grid cell 
c     anpcpcur_cmp- annual precipitation for current year (mm) in each latitudinal grid cell  
c     anpecur_cmp - annual potential evaporation for current year (mm) in each 
c                   latitudinal grid cell 
c     gdd5cur_cmp - growing degree days above 5 c for current year in each 
c                   latitudinal grid cell 
c     surmncur_cmp- number of months with surplus water for current year in each 
c                   latitudinal grid cell 
c     defmncur_cmp- number of months with water deficit for current year in each 
c                   latitudinal grid cell 
c     srplscur_cmp- water surplus for the current month in each latitudinal grid cell  
c     defctcur_cmp- water deficit for the current month in each latitudinal grid cell 
c     twarmm_cmp  - temperature of the warmest month (c) in each latitudinal grid cell 
c     tcoldm_cmp  - temperature of the coldest month (c) in each latitudinal grid cell 
c     gdd5_cmp    - growing degree days above 5 c in each latitudinal grid cell 
c     aridity_cmp - aridity index, ratio of potential evaporation to
c                   precipitation in each latitudinal grid cell 
c     srplsmon_cmp- number of months in a year with surplus water i.e.
c                   precipitation more than potential evaporation in each 
c                   latitudinal grid cell 
c     defctmon_cmp- number of months in a year with water deficit i.e.
c                   precipitation less than potential evaporation in each 
c                   latitudinal grid cell 
c     anndefct_cmp- annual water deficit (mm) in each latitudinal grid cell   
c     annsrpls_cmp- annual water surplus (mm) in each latitudinal grid cell 
c     annpcp_cmp  - annual precipitation (mm) in each latitudinal grid cell 
c     dry_season_length_cmp - length of dry season (months) in each latitudional grid cell
c     lucemcom_cmp- land use change (luc) related combustion emission losses
c                   in each latitudional grid cell, u-mol co2/m2.sec 
c     lucltrin_cmp- luc related inputs to litter pool, in each latitudional 
c                   grid cell, u-mol co2/m2.sec
c     lucsocin_cmp- luc related inputs to soil c pool, in each latitudional 
c                   grid cell, u-mol co2/m2.sec
c     todfrac_cmp - max. fractional coverage of ctem's 9 pfts by the end
c                   of the dayin each latitudinal grid cell, for use by land use subroutine
c     pfcancmx_cmp- previous year's fractional coverages of pfts in each latitudinal grid cell
c     nfcancmx_cmp- next year's fractional coverages of pfts in each latitudinal grid cell
c     pstemmass_cmp - stem mass from previous timestep, is value before fire. used by burntobare subroutine
c     pgleafmass_cmp - root mass from previous timestep, is value before fire. used by burntobare subroutine


c
c    updates
c 
c     faregat - fractional coverage of each ctem pft in each mosaic tile   
c     fcancmx - fractional coverage of ctem's 10 pfts in each mosaic 
c     nppveg  - npp for each pft type /m2 of vegetated area [u-mol co2-c/m2.sec]
c     geremort- growth related mortality (1/day)
c     intrmort- intrinsic (age related) mortality (1/day)
c     gleafmas- green leaf mass for each of the 10 ctem pfts, kg c/m2
c     bleafmas- brown leaf mass for each of the 10 ctem pfts, kg c/m2
c     stemmass- stem mass for each of the 10 ctem pfts, kg c/m2
c     rootmass- root mass for each of the 10 ctem pfts, kg c/m2
c     litrmass- litter mass for each of the 10 ctem pfts + bare, kg c/m2
c     soilcmas- soil carbon mass for each of the 10 ctem pfts + bare, kg c/m2
c     pftexist- logical array indicating pfts exist (t) or not (f)
c     lambda  - fraction of npp that is used for horizontal expansion
c     bmasveg - total (gleaf + stem + root) biomass for each ctem pft, kg c/m2
c     burnvegf - fractional areas burned,for ctem pfts
c     add2allo- npp kg c/m2.day that is used for expansion and
c               subsequently allocated to leaves, stem, and root via 
c               the allocation part of the model.
c     cc,mm   - colonization rate & mortality rate 
c     fcanmx  - fractional coverage of class' 4 pfts
c     vgbiomas- grid averaged vegetation biomass, kg c/m2
c     gavgltms- grid averaged litter mass, kg c/m2
c     gavgscms- grid averaged soil c mass, kg c/m2
c     lucemcom - land use change (luc) related combustion emission losses,
c                u-mol co2/m2.sec 
c     lucltrin - luc related inputs to litter pool, u-mol co2/m2.sec
c     lucsocin - luc related inputs to soil c pool, u-mol co2/m2.sec
c     todfrac  - max. fractional coverage of ctem's 9 pfts by the end
c                of the day, for use by land use subroutine
c     pfcancmx - previous year's fractional coverages of pfts
c     nfcancmx - next year's fractional coverages of pfts
c     grclarea - area of the grid cell, km^2
c     ta      - mean daily temperature, k
c     precip  - daily precipitation (mm/day)
c     netrad  - daily net radiation (w/m2)
c     tcurm   - temperature of the current month (c)
c     srpcuryr- water surplus for the current year
c     dftcuryr- water deficit for the current year
c     tmonth  - monthly temperatures
c     anpcpcur- annual precipitation for current year (mm)
c     anpecur - annual potential evaporation for current year (mm)
c     gdd5cur - growing degree days above 5 c for current year
c     surmncur- number of months with surplus water for current year
c     defmncur- number of months with water deficit for current year
c     srplscur- water surplus for the current month 
c     defctcur- water deficit for the current month
c     twarmm  - temperature of the warmest month (c)
c     tcoldm  - temperature of the coldest month (c)
c     gdd5    - growing degree days above 5 c
c     aridity - aridity index, ratio of potential evaporation to
c               precipitation
c     srplsmon- number of months in a year with surplus water i.e.
c               precipitation more than potential evaporation
c     defctmon- number of months in a year with water deficit i.e.
c               precipitation less than potential evaporation
c     anndefct- annual water deficit (mm) 
c     annsrpls- annual water surplus (mm)
c     annpcp  - annual precipitation (mm)
c     dry_season_length - length of dry season (months)
c
      use ctem_params,        only : nlat,nmos,icc, ilg, ican,iccp1

      implicit none
c
      integer i,m,l,k,mn,
     1        nml,ilmos(ilg),jlmos(ilg),nol2pfts(ican)
c
c--------input arrays for unmapping------------------------------------------
c
      real  fare_cmp(nlat,icc), nppveg_cmp(nlat,icc),
     1      geremort_cmp(nlat,icc),intrmort_cmp(nlat,icc),
     2      gleafmas_cmp(nlat,icc),bleafmas_cmp(nlat,icc),
     3      stemmass_cmp(nlat,icc),rootmass_cmp(nlat,icc),
     4      litrmass_cmp(nlat,iccp1),soilcmas_cmp(nlat,iccp1),
     5      lambda_cmp(nlat,icc),    todfrac_cmp(nlat,icc),
     6      bmasveg_cmp(nlat,icc),   burnvegf_cmp(nlat,icc),
     7      add2allo_cmp(nlat,icc),  cc_cmp(nlat,icc),mm_cmp(nlat,icc),
     8      fcanmx_cmp(nlat,ican),   grclarea_cmp(nlat),  
     9      vgbiomas_cmp(nlat),    
     1      gavgltms_cmp(nlat),      gavgscms_cmp(nlat),
     2      lucemcom_cmp(nlat),  lucltrin_cmp(nlat), lucsocin_cmp(nlat),
     3      pfcancmx_cmp(nlat,icc), nfcancmx_cmp(nlat,icc),
     3      pstemmass_cmp(nlat,icc), pgleafmass_cmp(nlat,icc)
c
      logical pftexist_cmp(nlat,icc)
c
      real  ta_cmp(nlat),       precip_cmp(nlat),  netrad_cmp(nlat), 
     1      tcurm_cmp(nlat),    srpcuryr_cmp(nlat),dftcuryr_cmp(nlat),
     2      tmonth_cmp(12,nlat),anpcpcur_cmp(nlat),anpecur_cmp(nlat), 
     3      gdd5cur_cmp(nlat),  surmncur_cmp(nlat),defmncur_cmp(nlat),
     4      srplscur_cmp(nlat), defctcur_cmp(nlat),twarmm_cmp(nlat), 
     5      tcoldm_cmp(nlat),   gdd5_cmp(nlat),    aridity_cmp(nlat),
     6      srplsmon_cmp(nlat), defctmon_cmp(nlat),anndefct_cmp(nlat),
     7      annsrpls_cmp(nlat), annpcp_cmp(nlat),
     8      dry_season_length_cmp(nlat)
c
c--------intermediate arrays for unmapping-----------------------------------
c
      real  fcancmxrow(nlat,nmos,icc), nppvegrow(nlat,nmos,icc),
     1      geremortrow(nlat,nmos,icc),intrmortrow(nlat,nmos,icc),
     2      gleafmasrow(nlat,nmos,icc),bleafmasrow(nlat,nmos,icc),
     3      stemmassrow(nlat,nmos,icc),rootmassrow(nlat,nmos,icc),
     4      litrmassrow(nlat,nmos,iccp1),soilcmasrow(nlat,nmos,iccp1),
     5      lambdarow(nlat,nmos,icc),    todfracrow(nlat,nmos,icc),
     6      bmasvegrow(nlat,nmos,icc),   burnvegfrow(nlat,nmos,icc),
     7      add2allorow(nlat,nmos,icc),  ccrow(nlat,nmos,icc),
     8      mmrow(nlat,nmos,icc),        fcanmxrow(nlat,nmos,ican),
     9      farerow(nlat,nmos),          grclarearow(nlat,nmos),
     1      vgbiomasrow(nlat,nmos),    
     2      gavgltmsrow(nlat,nmos),      gavgscmsrow(nlat,nmos),
     3      lucemcomrow(nlat,nmos),      lucltrinrow(nlat,nmos),
     4      lucsocinrow(nlat,nmos),
     5      pfcancmxrow(nlat,nmos,icc), nfcancmxrow(nlat,nmos,icc),
     5      pstemmassrow(nlat,nmos,icc), pgleafmassrow(nlat,nmos,icc)
c
      logical pftexistrow(nlat,nmos,icc)      
c
c--------these intermediate arrays were transferred from mapping------------\\
c
      real  netradrow(nlat,nmos)
c
c---------------------------------------------------------------------------//
c
      real  tarow(nlat,nmos),       preciprow(nlat,nmos),
     1      tcurmrow(nlat,nmos),
     2      srpcuryrrow(nlat,nmos), dftcuryrrow(nlat,nmos),
     3      tmonthrow(12,nlat,nmos),anpcpcurrow(nlat,nmos), 
     4      anpecurrow(nlat,nmos),  gdd5currow(nlat,nmos), 
     5      surmncurrow(nlat,nmos), defmncurrow(nlat,nmos),
     6      srplscurrow(nlat,nmos), defctcurrow(nlat,nmos), 
     7      twarmmrow(nlat,nmos),   tcoldmrow(nlat,nmos),  
     8      gdd5row(nlat,nmos),     aridityrow(nlat,nmos),
     9      srplsmonrow(nlat,nmos), defctmonrow(nlat,nmos), 
     1      anndefctrow(nlat,nmos), annsrplsrow(nlat,nmos),
     2      annpcprow(nlat,nmos),
     3      dry_season_lengthrow(nlat,nmos)
c
c--------updated arrays after unmapping--------------------------------------
c
      real  fcancmx(ilg,icc), faregat(ilg),
     1      nppveg(ilg,icc),  geremort(ilg,icc),  
     2      intrmort(ilg,icc),gleafmas(ilg,icc),
     3      bleafmas(ilg,icc),stemmass(ilg,icc),
     4      rootmass(ilg,icc),litrmass(ilg,iccp1),
     5      soilcmas(ilg,iccp1),
     6      lambda(ilg,icc),  todfrac(ilg,icc),
     7      bmasveg(ilg,icc), burnvegf(ilg,icc),
     8      add2allo(ilg,icc),cc(ilg,icc),mm(ilg,icc),
     9      fcanmx(ilg,ican),   
     1      vgbiomas(ilg),    gavgltms(ilg),
     2      gavgscms(ilg),  grclarea(ilg),
     3      lucemcom(ilg),  lucltrin(ilg), lucsocin(ilg),
     4      pfcancmx(ilg,icc), nfcancmx(ilg,icc),
     4      pstemmass(ilg,icc), pgleafmass(ilg,icc)
c
      logical pftexist(ilg,icc)
c
      real  ta(ilg),        precip(ilg),   netrad(ilg),
     1      tcurm(ilg),     srpcuryr(ilg), dftcuryr(ilg),
     2      tmonth(12,ilg), anpcpcur(ilg), anpecur(ilg),
     3      gdd5cur(ilg),   surmncur(ilg), defmncur(ilg),
     4      srplscur(ilg),  defctcur(ilg), twarmm(ilg),
     5      tcoldm(ilg),    gdd5(ilg),     aridity(ilg),
     6      srplsmon(ilg),  defctmon(ilg), anndefct(ilg),
     7      annsrpls(ilg),  annpcp(ilg),
     8      dry_season_length(ilg) 
c
c--------internal arrays-----------------------------------------------------
c
      integer mcount, mcount1,mcount2
      real tpftfrac(nlat)
      logical bareexist(nlat)
c
c     ------------------------------------------------------------------
c                           parameters used 
c
c     note the structure of parameter vectors which clearly shows the
c     class pfts (along rows) and ctem sub-pfts (along columns)
c
c     needle leaf |  evg1      evg2      dcd
c     broad leaf  |  evg   dcd-cld   dcd-dry
c     crops       |   c3        c4       ---
c     grasses     |   c3        c4       ---
c
c     ---------------------------------------------------------------
c
      if(icc.ne.9)                     call xit('compete_map',-1)
      if(ican.ne.4)                       call xit('compete_map',-2)
c
c     initialization
c
      do 30 i=1,nlat
       tpftfrac(i)=0.0
       do 31 m=1,nmos
c
        do l=1,icc
         fcancmxrow(i,m,l) = 0.0
         nppvegrow(i,m,l)  = 0.0
         geremortrow(i,m,l)= 0.0
         intrmortrow(i,m,l)= 0.0
         gleafmasrow(i,m,l)= 0.0
         bleafmasrow(i,m,l)= 0.0
         stemmassrow(i,m,l)= 0.0
         rootmassrow(i,m,l)= 0.0
         pftexistrow(i,m,l)= .false.
         lambdarow(i,m,l)  = 0.0
         todfracrow(i,m,l) = 0.0
         pfcancmxrow(i,m,l)= 0.0
         nfcancmxrow(i,m,l)= 0.0
         pstemmassrow(i,m,l)= 0.0
         pgleafmassrow(i,m,l)= 0.0
         bmasvegrow(i,m,l) = 0.0
         burnvegfrow(i,m,l) = 0.0
         add2allorow(i,m,l)= 0.0
         ccrow(i,m,l)      = 0.0
         mmrow(i,m,l)      = 0.0
        enddo 
c
        do l=1,iccp1
         litrmassrow(i,m,l)= 0.0
         soilcmasrow(i,m,l)= 0.0
        enddo
c
        do l=1,ican
         fcanmxrow(i,m,l)  = 0.0            
        enddo
c
         farerow(i,m)      = 0.0
         grclarearow(i,m)  = 0.0
         vgbiomasrow(i,m)  = 0.0
         gavgltmsrow(i,m)  = 0.0
         gavgscmsrow(i,m)  = 0.0  
         lucemcomrow(i,m)  = 0.0
         lucltrinrow(i,m)  = 0.0
         lucsocinrow(i,m)  = 0.0 
         tarow(i,m)        = 0.0  
         preciprow(i,m)    = 0.0  
c        note that netradrow is not initialized here 
c        because it is updated
         tcurmrow(i,m)     = 0.0  
         srpcuryrrow(i,m)  = 0.0  
         dftcuryrrow(i,m)  = 0.0  
c
         do mn=1,12
         tmonthrow(mn,i,m) = 0.0
         enddo
c
         anpcpcurrow(i,m)  = 0.0  
         anpecurrow(i,m)   = 0.0  
         gdd5currow(i,m)   = 0.0  
         surmncurrow(i,m)  = 0.0  
         defmncurrow(i,m)  = 0.0  
         srplscurrow(i,m)  = 0.0    
         defctcurrow(i,m)  = 0.0  
         twarmmrow(i,m)    = 0.0  
         tcoldmrow(i,m)    = 0.0  
         gdd5row(i,m)      = 0.0   
         aridityrow(i,m)   = 0.0  
         srplsmonrow(i,m)  = 0.0  
         defctmonrow(i,m)  = 0.0  
         anndefctrow(i,m)  = 0.0  
         annsrplsrow(i,m)  = 0.0  
         annpcprow(i,m)    = 0.0  
         dry_season_lengthrow(i,m) = 0.0
c
31     continue
30    continue
c
      do 40 i=1,ilg
c
       do l=1,icc
        fcancmx(i,l) = 0.0
        nppveg(i,l)  = 0.0
        geremort(i,l)= 0.0
        intrmort(i,l)= 0.0
        gleafmas(i,l)= 0.0
        bleafmas(i,l)= 0.0
        stemmass(i,l)= 0.0
        rootmass(i,l)= 0.0
        pftexist(i,l)= .false.
        lambda(i,l)  = 0.0
        todfrac(i,l) = 0.0
        pfcancmx(i,l)= 0.0
        nfcancmx(i,l)= 0.0
        pstemmass(i,l)= 0.0
        pgleafmass(i,l)= 0.0
        bmasveg(i,l) = 0.0
        burnvegf(i,l) = 0.0
        add2allo(i,l)= 0.0
        cc(i,l)      = 0.0
        mm(i,l)      = 0.0
       enddo
c
       do l=1,iccp1
        litrmass(i,l)= 0.0
        soilcmas(i,l)= 0.0
       enddo
c
       do l=1,ican
        fcanmx(i,l)  = 0.0  
       enddo
c
       faregat(i)   = 0.0
       grclarea(i)  = 0.0
       vgbiomas(i)  = 0.0
       gavgltms(i)  = 0.0
       gavgscms(i)  = 0.0 
       lucemcom(i)  = 0.0 
       lucltrin(i)  = 0.0
       lucsocin(i)  = 0.0
       ta(i)        = 0.0  
       precip(i)    = 0.0  
       netrad(i)    = 0.0  
       tcurm(i)     = 0.0  
       srpcuryr(i)  = 0.0  
       dftcuryr(i)  = 0.0  
c
       do mn=1,12
        tmonth(mn,i) = 0.0
       enddo
c
       anpcpcur(i)  = 0.0  
       anpecur(i)   = 0.0  
       gdd5cur(i)   = 0.0  
       surmncur(i)  = 0.0  
       defmncur(i)  = 0.0  
       srplscur(i)  = 0.0    
       defctcur(i)  = 0.0  
       twarmm(i)    = 0.0  
       tcoldm(i)    = 0.0  
       gdd5(i)      = 0.0   
       aridity(i)   = 0.0  
       srplsmon(i)  = 0.0  
       defctmon(i)  = 0.0  
       anndefct(i)  = 0.0  
       annsrpls(i)  = 0.0  
       annpcp(i)    = 0.0  
       dry_season_length(i) = 0.0
40    continue
c
c     unmapping the pft fraction in each grid cell (fare_cmp) 
c     back to the pft index in each mosaic of each grid cell (fcancmxrow) 
c     and update the pft areal fraction in each mosaic of each grid cell (farerow)
c
c     check if bare fraction is existed 
c 
      do 50 i=1,nlat
c
       do 55 l=1,icc
        tpftfrac(i)=tpftfrac(i)+fare_cmp(i,l)
55     continue
c
       if (tpftfrac(i) .lt. 1.) then 
        bareexist(i) = .true.
       else   
        bareexist(i) = .false.  
         write(*,*)'minimal bare fraction has been eliminated'
         stop
       endif
c
50    continue
c
c     unmapping back to fcancmxrow and
c     update farerow
c
      do 100 i=1,nlat
c
         if (bareexist(i)) then 
          farerow(i,nmos)=1.000-tpftfrac(i)
         endif
c
         do l=1,icc
          if(fare_cmp(i,l) .gt. 0.0) then
           mcount=l
           fcancmxrow(i,mcount,l) = 1.
           farerow(i,mcount)      = fare_cmp(i,l)
           nppvegrow(i,mcount,l)  = nppveg_cmp(i,l)
           geremortrow(i,mcount,l)= geremort_cmp(i,l)
           intrmortrow(i,mcount,l)= intrmort_cmp(i,l)
           gleafmasrow(i,mcount,l)= gleafmas_cmp(i,l)
           bleafmasrow(i,mcount,l)= bleafmas_cmp(i,l)
           stemmassrow(i,mcount,l)= stemmass_cmp(i,l)
           rootmassrow(i,mcount,l)= rootmass_cmp(i,l)
           pftexistrow(i,mcount,l)= pftexist_cmp(i,l)
           lambdarow(i,mcount,l)  = lambda_cmp(i,l)
           todfracrow(i,mcount,l) = todfrac_cmp(i,l)
           pfcancmxrow(i,mcount,l)= pfcancmx_cmp(i,l)
           nfcancmxrow(i,mcount,l)= nfcancmx_cmp(i,l)
           pstemmassrow(i,mcount,l)= pstemmass_cmp(i,l)
           pgleafmassrow(i,mcount,l)= pgleafmass_cmp(i,l)
           bmasvegrow(i,mcount,l) = bmasveg_cmp(i,l)
           burnvegfrow(i,mcount,l) = burnvegf_cmp(i,l)
           add2allorow(i,mcount,l)= add2allo_cmp(i,l)
           ccrow(i,mcount,l)      = cc_cmp(i,l)
           mmrow(i,mcount,l)      = mm_cmp(i,l)
           vgbiomasrow(i,mcount)  = gleafmas_cmp(i,l)+bleafmas_cmp(i,l)+
     &                              stemmass_cmp(i,l)+rootmass_cmp(i,l)
           lucemcomrow(i,mcount) = lucemcom_cmp(i) 
           lucltrinrow(i,mcount) = lucltrin_cmp(i) 
           lucsocinrow(i,mcount) = lucsocin_cmp(i) 

          endif 
         enddo
c
         do m=1,nmos
           tarow(i,m)        = ta_cmp(i)
           preciprow(i,m)    = precip_cmp(i)
           tcurmrow(i,m)     = tcurm_cmp(i)
           srpcuryrrow(i,m ) = srpcuryr_cmp(i)
           dftcuryrrow(i,m)  = dftcuryr_cmp(i)
           grclarearow(i,m)  = grclarea_cmp(i)
c
           do mn=1,12
            tmonthrow(mn,i,m) = tmonth_cmp(mn,i)
           enddo
c
           anpcpcurrow(i,m)  = anpcpcur_cmp(i)
           anpecurrow(i,m)   = anpecur_cmp(i)  
           gdd5currow(i,m)   = gdd5cur_cmp(i)  
           surmncurrow(i,m)  = surmncur_cmp(i)  
           defmncurrow(i,m)  = defmncur_cmp(i)  
           srplscurrow(i,m)  = srplscur_cmp(i)    
           defctcurrow(i,m)  = defctcur_cmp(i)  
           twarmmrow(i,m)    = twarmm_cmp(i) 
           tcoldmrow(i,m)    = tcoldm_cmp(i)  
           gdd5row(i,m)      = gdd5_cmp(i)  
           aridityrow(i,m)   = aridity_cmp(i)  
           srplsmonrow(i,m)  = srplsmon_cmp(i)  
           defctmonrow(i,m)  = defctmon_cmp(i)  
           anndefctrow(i,m)  = anndefct_cmp(i) 
           annsrplsrow(i,m)  = annsrpls_cmp(i)   
           annpcprow(i,m)    = annpcp_cmp(i) 
           dry_season_lengthrow(i,m) = dry_season_length_cmp(i)

         enddo
c
         do l=1,iccp1
          if(litrmass_cmp(i,l) .gt. 0.0) then
           mcount1 = l
           litrmassrow(i,mcount1,l)= litrmass_cmp(i,l)
           soilcmasrow(i,mcount1,l)= soilcmas_cmp(i,l)
           gavgltmsrow(i,mcount1)  = litrmass_cmp(i,l) 
           gavgscmsrow(i,mcount1)  = soilcmas_cmp(i,l)
          endif
         enddo
c
         mcount2=1
         do l=1,ican
          if(fcanmx_cmp(i,l) .gt. 0.0) then
c
           do m=mcount2,mcount2+nol2pfts(l)-1                          
            fcanmxrow(i,m,l)=1.

           enddo
c
           mcount2=mcount2+nol2pfts(l)
          endif
         enddo

 100  continue
c
c     gathering the pft index in each mosaic of each grid cell (fcancmxrow) 
c     to the pft index in each mosaic (fcancmx) 
c     nml, ilmos and jlmos are referring to gatprep.f
c           
      do 200 l=1,icc
       do 200 k=1,nml
         fcancmx(k,l)  = fcancmxrow(ilmos(k),jlmos(k),l)
         nppveg(k,l)   = nppvegrow(ilmos(k),jlmos(k),l)
         geremort(k,l) = geremortrow(ilmos(k),jlmos(k),l)
         intrmort(k,l) = intrmortrow(ilmos(k),jlmos(k),l)
         gleafmas(k,l) = gleafmasrow(ilmos(k),jlmos(k),l)
         bleafmas(k,l) = bleafmasrow(ilmos(k),jlmos(k),l)
         stemmass(k,l) = stemmassrow(ilmos(k),jlmos(k),l)
         rootmass(k,l) = rootmassrow(ilmos(k),jlmos(k),l)
         pftexist(k,l) = pftexistrow(ilmos(k),jlmos(k),l)
         lambda(k,l)   = lambdarow(ilmos(k),jlmos(k),l)
         todfrac(k,l)  = todfracrow(ilmos(k),jlmos(k),l)
         pfcancmx(k,l) = pfcancmxrow(ilmos(k),jlmos(k),l)
         nfcancmx(k,l) = nfcancmxrow(ilmos(k),jlmos(k),l)
         pstemmass(k,l) = pstemmassrow(ilmos(k),jlmos(k),l)
         pgleafmass(k,l) = pgleafmassrow(ilmos(k),jlmos(k),l)
         bmasveg(k,l)  = bmasvegrow(ilmos(k),jlmos(k),l)
         burnvegf(k,l)  = burnvegfrow(ilmos(k),jlmos(k),l)
         add2allo(k,l) = add2allorow(ilmos(k),jlmos(k),l)
         cc(k,l)       = ccrow(ilmos(k),jlmos(k),l)
         mm(k,l)       = mmrow(ilmos(k),jlmos(k),l)
 200  continue
c
      do 210 l=1,iccp1
       do 210 k=1,nml
         litrmass(k,l) = litrmassrow(ilmos(k),jlmos(k),l)
         soilcmas(k,l) = soilcmasrow(ilmos(k),jlmos(k),l)
 210  continue
c
      do 220 l=1,ican
       do 220 k=1,nml
         fcanmx(k,l) = fcanmxrow(ilmos(k),jlmos(k),l)
 220  continue
c
      do 230 k=1,nml
        faregat(k)  = farerow(ilmos(k),jlmos(k))
        grclarea(k) = grclarearow(ilmos(k),jlmos(k))
        vgbiomas(k) = vgbiomasrow(ilmos(k),jlmos(k))
        gavgltms(k) = gavgltmsrow(ilmos(k),jlmos(k))
        gavgscms(k) = gavgscmsrow(ilmos(k),jlmos(k))
        lucemcom(k) = lucemcomrow(ilmos(k),jlmos(k))
        lucltrin(k) = lucltrinrow(ilmos(k),jlmos(k))
        lucsocin(k) = lucsocinrow(ilmos(k),jlmos(k))
        ta(k)       = tarow(ilmos(k),jlmos(k))
        precip(k)   = preciprow(ilmos(k),jlmos(k))
        netrad(k)   = netradrow(ilmos(k),jlmos(k))
        tcurm(k)    = tcurmrow(ilmos(k),jlmos(k))
        srpcuryr(k) = srpcuryrrow(ilmos(k),jlmos(k))
        dftcuryr(k) = dftcuryrrow(ilmos(k),jlmos(k))
c
        do mn=1,12
          tmonth(mn,k) = tmonthrow(mn,ilmos(k),jlmos(k))
        enddo
c
        anpcpcur(k) = anpcpcurrow(ilmos(k),jlmos(k))
        anpecur(k)  = anpecurrow(ilmos(k),jlmos(k))  
        gdd5cur(k)  = gdd5currow(ilmos(k),jlmos(k))  
        surmncur(k) = surmncurrow(ilmos(k),jlmos(k))  
        defmncur(k) = defmncurrow(ilmos(k),jlmos(k))  
        srplscur(k) = srplscurrow(ilmos(k),jlmos(k))    
        defctcur(k) = defctcurrow(ilmos(k),jlmos(k))  
        twarmm(k)   = twarmmrow(ilmos(k),jlmos(k)) 
        tcoldm(k)   = tcoldmrow(ilmos(k),jlmos(k))  
        gdd5(k)     = gdd5row(ilmos(k),jlmos(k))  
        aridity(k)  = aridityrow(ilmos(k),jlmos(k))  
        srplsmon(k) = srplsmonrow(ilmos(k),jlmos(k))  
        defctmon(k) = defctmonrow(ilmos(k),jlmos(k))  
        anndefct(k) = anndefctrow(ilmos(k),jlmos(k)) 
        annsrpls(k) = annsrplsrow(ilmos(k),jlmos(k))   
        annpcp(k)   = annpcprow(ilmos(k),jlmos(k)) 
        dry_season_length(k) = dry_season_lengthrow(ilmos(k),jlmos(k))
 230  continue
c
      return
      end


