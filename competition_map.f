      subroutine competition_map(nml, ilmos, jlmos, grclarea,
     2                           faregat,  fcancmx,  nppveg,geremort,
     3                           intrmort,gleafmas,bleafmas,stemmass,
     4                           rootmass,litrmass,soilcmas,
     5                           pftexist,  lambda, bmasveg,burnvegf,
     6                           add2allo,      cc,      mm,  fcanmx,
     7                           vgbiomas,gavgltms,gavgscms,
     8                           ta,precip, netrad,   tcurm,srpcuryr,
     9                           dftcuryr,  tmonth,anpcpcur, anpecur, 
     1                            gdd5cur,surmncur,defmncur,srplscur,  
     2                           defctcur,  twarmm,  tcoldm,    gdd5, 
     3                            aridity,srplsmon,defctmon,anndefct,
     4                           annsrpls,  annpcp,anpotevp,
     5                           lucemcom,  lucltrin, lucsocin,
     6                           pfcancmx,  nfcancmx,
c    ------------------- inputs above this line ---------------------
     a                             netradrow,
c    ------------------- intermediate and to save above this line ---
     a                           fare_cmp,  nppveg_cmp,geremort_cmp,
     b                          intrmort_cmp,gleafmas_cmp,bleafmas_cmp,
     c                          stemmass_cmp,rootmass_cmp,litrmass_cmp,
     d                          soilcmas_cmp,pftexist_cmp,  lambda_cmp,
     e                         bmasveg_cmp,burnvegf_cmp,add2allo_cmp,
     f                                cc_cmp,      mm_cmp,  fcanmx_cmp,
     g                          vgbiomas_cmp, grclarea_cmp,
     h                          gavgltms_cmp,gavgscms_cmp,
     i                                ta_cmp,  precip_cmp,  netrad_cmp, 
     j                             tcurm_cmp,srpcuryr_cmp,dftcuryr_cmp,
     k                            tmonth_cmp,anpcpcur_cmp, anpecur_cmp, 
     l                           gdd5cur_cmp,surmncur_cmp,defmncur_cmp,
     m                          srplscur_cmp,defctcur_cmp,  twarmm_cmp, 
     n                            tcoldm_cmp,    gdd5_cmp, aridity_cmp,
     o                          srplsmon_cmp,defctmon_cmp,anndefct_cmp,
     p                          annsrpls_cmp,  annpcp_cmp,anpotevp_cmp,
     q                         lucemcom_cmp, lucltrin_cmp, lucsocin_cmp,
     r                           pfcancmx_cmp,  nfcancmx_cmp)      
c    ------------------- outputs above this line----------------------   
c
c
C               Canadian Terrestrial Ecosystem Model (CTEM) 
C                    Mapping For Competition Subroutine
c
c    23  Jul 2013 - Add in module for parameters
C    J. Melton

c    28 Aug 2012 - this subroutine prepares for the competition 
c    Y. Peng       calculation by mapping the only pft in each 
c                  mosaic to pft fractions in each grid cell
c
c                  input:        array(ilg,icc)
c                    | 
c                    | scattering      
c                    v
c                  intermediate: arrayrow(nlat,nmos,icc)
c                    | 
c                    | mapping
c                    v
c                  output:       array_cmp(nlat,icc)
c
c                  array_cmp: mapped array prepared for competition
c
c    -----------------------------------------------------------------
c
c    indices
c
c     nlat    - max. number of grid cells in the latitude circle, which
c               is prescribed in runclass36ctem.f
c     nmos    - max. number of mosaic tiles in each latitudinal grid cell, 
c               which is prescribed in runclass35ctem.f
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
c     burnvegf - fractional areas burned, for ctem pfts
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
c
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
c     anpotevp- annual potential evaporation (mm)
c
c    outputs
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
c     grclarea_cmp - area of the grid cell, km^2
c     gavgltms_cmp- grid averaged litter mass, kg c/m2
c     gavgscms_cmp- grid averaged soil c mass, kg c/m2
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
c     anpotevp_cmp- annual potential evaporation (mm) in each latitudinal grid cell 
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
c
      use ctem_params,        only : nlat,nmos,icc, ilg, ican,iccp1

      implicit none
c
      integer i,m,l,k,mn,
     1        nml,ilmos(ilg),jlmos(ilg)
c
c--------input arrays for mapping------------------------------------------
c
      real  fcancmx(ilg,icc), faregat(ilg), grclarea(ilg),
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
     2      gavgscms(ilg),
     3      lucemcom(ilg),  lucltrin(ilg), lucsocin(ilg),
     4      pfcancmx(ilg,icc), nfcancmx(ilg,icc)
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
     7      annsrpls(ilg),  annpcp(ilg),   anpotevp(ilg)
 
c
c--------intermediate arrays for mapping-----------------------------------
c
      real  fcancmxrow(nlat,nmos,icc),   nppvegrow(nlat,nmos,icc),
     1      geremortrow(nlat,nmos,icc),  intrmortrow(nlat,nmos,icc),
     2      gleafmasrow(nlat,nmos,icc),  bleafmasrow(nlat,nmos,icc),
     3      stemmassrow(nlat,nmos,icc),  rootmassrow(nlat,nmos,icc),
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
     5      pfcancmxrow(nlat,nmos,icc), nfcancmxrow(nlat,nmos,icc)
c
      logical pftexistrow(nlat,nmos,icc)
c
c--------these intermediate arrays will be saved for unmapping------------\\
c
      real  netradrow(nlat,nmos)
c
c-------------------------------------------------------------------------//
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
     2      annpcprow(nlat,nmos),   anpotevprow(nlat,nmos) 
c
c--------output arrays after mapping---------------------------------------
c
      real  fare_cmp(nlat,icc),   nppveg_cmp(nlat,icc),
     1      geremort_cmp(nlat,icc),  intrmort_cmp(nlat,icc),
     2      gleafmas_cmp(nlat,icc),  bleafmas_cmp(nlat,icc),
     3      stemmass_cmp(nlat,icc),  rootmass_cmp(nlat,icc),
     4      litrmass_cmp(nlat,iccp1),soilcmas_cmp(nlat,iccp1),
     5      lambda_cmp(nlat,icc),    todfrac_cmp(nlat,icc),
     6      bmasveg_cmp(nlat,icc),   burnvegf_cmp(nlat,icc),
     7      add2allo_cmp(nlat,icc),  cc_cmp(nlat,icc),mm_cmp(nlat,icc),
     8      fcanmx_cmp(nlat,ican),   grclarea_cmp(nlat),  
     9      vgbiomas_cmp(nlat),    
     1      gavgltms_cmp(nlat),      gavgscms_cmp(nlat),
     2      lucemcom_cmp(nlat),  lucltrin_cmp(nlat), lucsocin_cmp(nlat),
     3      pfcancmx_cmp(nlat,icc), nfcancmx_cmp(nlat,icc)
      logical pftexist_cmp(nlat,icc)
      real  ta_cmp(nlat),       precip_cmp(nlat),  netrad_cmp(nlat), 
     1      tcurm_cmp(nlat),    srpcuryr_cmp(nlat),dftcuryr_cmp(nlat),
     2      tmonth_cmp(12,nlat),anpcpcur_cmp(nlat),anpecur_cmp(nlat), 
     3      gdd5cur_cmp(nlat),  surmncur_cmp(nlat),defmncur_cmp(nlat),
     4      srplscur_cmp(nlat), defctcur_cmp(nlat),twarmm_cmp(nlat), 
     5      tcoldm_cmp(nlat),   gdd5_cmp(nlat),    aridity_cmp(nlat),
     6      srplsmon_cmp(nlat), defctmon_cmp(nlat),anndefct_cmp(nlat),
     7      annsrpls_cmp(nlat), annpcp_cmp(nlat),  anpotevp_cmp(nlat)
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
      if(icc.ne.9)                    call xit('compete_unmap',-1)
      if(ican.ne.4)                     call xit('compete_unmap',-2)
c
c     initialization
c
      do 90 i = 1, nlat
       do 91 m = 1, nmos
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
         vgbiomasrow(i,m)  = 0.0
         gavgltmsrow(i,m)  = 0.0
         gavgscmsrow(i,m)  = 0.0 
         lucemcomrow(i,m)  = 0.0
         lucltrinrow(i,m)  = 0.0
         lucsocinrow(i,m)  = 0.0 
         tarow(i,m)        = 0.0  
         preciprow(i,m)    = 0.0  
         netradrow(i,m)    = 0.0  
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
         anpotevprow(i,m)  = 0.0  
91     continue
c
       do l=1,icc
        fare_cmp(i,l) = 0.0
        nppveg_cmp(i,l)  = 0.0
        geremort_cmp(i,l)= 0.0
        intrmort_cmp(i,l)= 0.0
        gleafmas_cmp(i,l)= 0.0
        bleafmas_cmp(i,l)= 0.0
        stemmass_cmp(i,l)= 0.0
        rootmass_cmp(i,l)= 0.0
        pftexist_cmp(i,l)= .false.
        lambda_cmp(i,l)  = 0.0
        todfrac_cmp(i,l) = 0.0
        pfcancmx_cmp(i,l)= 0.0
        nfcancmx_cmp(i,l)= 0.0
        bmasveg_cmp(i,l) = 0.0
        burnvegf_cmp(i,l) = 0.0
        add2allo_cmp(i,l)= 0.0
        cc_cmp(i,l)      = 0.0
        mm_cmp(i,l)      = 0.0
       enddo
c
       do l=1,iccp1
        litrmass_cmp(i,l)= 0.0
        soilcmas_cmp(i,l)= 0.0
       enddo
c
       do l=1,ican
        fcanmx_cmp(i,l)  = 0.0  
       enddo
c
        vgbiomas_cmp(i)  = 0.0
        gavgltms_cmp(i)  = 0.0
        gavgscms_cmp(i)  = 0.0 
        lucemcom_cmp(i)  = 0.0 
        lucltrin_cmp(i)  = 0.0
        lucsocin_cmp(i)  = 0.0
        ta_cmp(i)        = 0.0  
        precip_cmp(i)    = 0.0  
        netrad_cmp(i)    = 0.0  
        tcurm_cmp(i)     = 0.0  
        srpcuryr_cmp(i)  = 0.0  
        dftcuryr_cmp(i)  = 0.0  
c
        do mn=1,12
         tmonth_cmp(mn,i) = 0.0
        enddo
c
        anpcpcur_cmp(i)  = 0.0  
        anpecur_cmp(i)   = 0.0  
        gdd5cur_cmp(i)   = 0.0  
        surmncur_cmp(i)  = 0.0  
        defmncur_cmp(i)  = 0.0  
        srplscur_cmp(i)  = 0.0    
        defctcur_cmp(i)  = 0.0  
        twarmm_cmp(i)    = 0.0  
        tcoldm_cmp(i)    = 0.0  
        gdd5_cmp(i)      = 0.0   
        aridity_cmp(i)   = 0.0  
        srplsmon_cmp(i)  = 0.0  
        defctmon_cmp(i)  = 0.0  
        anndefct_cmp(i)  = 0.0  
        annsrpls_cmp(i)  = 0.0  
        annpcp_cmp(i)    = 0.0  
        anpotevp_cmp(i)  = 0.0  
        grclarea_cmp(i)  = 0.0
90    continue 
c
c     scattering the pft index in each mosaic (fcancmx) to 
c     pft index in each mosaic of each grid cell (fcancmxrow)
c     nml, ilmos and jlmos are referring to gatprep.f
c
      do 100 l=1,icc
       do 100 k=1,nml
         fcancmxrow(ilmos(k),jlmos(k),l)  = fcancmx(k,l)
         nppvegrow(ilmos(k),jlmos(k),l)   = nppveg(k,l)
         geremortrow(ilmos(k),jlmos(k),l) = geremort(k,l)
         intrmortrow(ilmos(k),jlmos(k),l) = intrmort(k,l)
         gleafmasrow(ilmos(k),jlmos(k),l) = gleafmas(k,l)
         bleafmasrow(ilmos(k),jlmos(k),l) = bleafmas(k,l)
         stemmassrow(ilmos(k),jlmos(k),l) = stemmass(k,l)
         rootmassrow(ilmos(k),jlmos(k),l) = rootmass(k,l)
         pftexistrow(ilmos(k),jlmos(k),l) = pftexist(k,l)
         lambdarow(ilmos(k),jlmos(k),l)   = lambda(k,l)
         todfracrow(ilmos(k),jlmos(k),l)  = todfrac(k,l)
         pfcancmxrow(ilmos(k),jlmos(k),l) = pfcancmx(k,l)
         nfcancmxrow(ilmos(k),jlmos(k),l) = nfcancmx(k,l)
         bmasvegrow(ilmos(k),jlmos(k),l)  = bmasveg(k,l)
         burnvegfrow(ilmos(k),jlmos(k),l)  = burnvegf(k,l)
         add2allorow(ilmos(k),jlmos(k),l) = add2allo(k,l)
         ccrow(ilmos(k),jlmos(k),l)       = cc(k,l)
         mmrow(ilmos(k),jlmos(k),l)       = mm(k,l)   
 100  continue
c
      do 110 l=1,iccp1
       do 110 k=1,nml
         litrmassrow(ilmos(k),jlmos(k),l) = litrmass(k,l)
         soilcmasrow(ilmos(k),jlmos(k),l) = soilcmas(k,l)
 110  continue
c
      do 120 l=1,ican
       do 120 k=1,nml
         fcanmxrow(ilmos(k),jlmos(k),l) = fcanmx(k,l)
 120  continue
c
      do 130 k=1,nml
         farerow(ilmos(k),jlmos(k))     = faregat(k)
         vgbiomasrow(ilmos(k),jlmos(k)) = vgbiomas(k)
         gavgltmsrow(ilmos(k),jlmos(k)) = gavgltms(k)
         gavgscmsrow(ilmos(k),jlmos(k)) = gavgscms(k) 
         lucemcomrow(ilmos(k),jlmos(k)) = lucemcom(k)
         lucltrinrow(ilmos(k),jlmos(k)) = lucltrin(k) 
         lucsocinrow(ilmos(k),jlmos(k)) = lucsocin(k)  
         tarow(ilmos(k),jlmos(k))       = ta(k)
         preciprow(ilmos(k),jlmos(k))   = precip(k)
         netradrow(ilmos(k),jlmos(k))   = netrad(k)
         tcurmrow(ilmos(k),jlmos(k))    = tcurm(k)
         srpcuryrrow(ilmos(k),jlmos(k)) = srpcuryr(k)
         dftcuryrrow(ilmos(k),jlmos(k)) = dftcuryr(k)
         grclarearow(ilmos(k),jlmos(k)) = grclarea(k)

         do mn=1,12
          tmonthrow(mn,ilmos(k),jlmos(k)) = tmonth(mn,k)
         enddo
c
         anpcpcurrow(ilmos(k),jlmos(k)) = anpcpcur(k)
         anpecurrow(ilmos(k),jlmos(k))  = anpecur(k)  
         gdd5currow(ilmos(k),jlmos(k))  = gdd5cur(k)  
         surmncurrow(ilmos(k),jlmos(k)) = surmncur(k)  
         defmncurrow(ilmos(k),jlmos(k)) = defmncur(k)  
         srplscurrow(ilmos(k),jlmos(k)) = srplscur(k)    
         defctcurrow(ilmos(k),jlmos(k)) = defctcur(k)  
         twarmmrow(ilmos(k),jlmos(k))   = twarmm(k) 
         tcoldmrow(ilmos(k),jlmos(k))   = tcoldm(k)  
         gdd5row(ilmos(k),jlmos(k))     = gdd5(k)  
         aridityrow(ilmos(k),jlmos(k))  = aridity(k)  
         srplsmonrow(ilmos(k),jlmos(k)) = srplsmon(k)  
         defctmonrow(ilmos(k),jlmos(k)) = defctmon(k)  
         anndefctrow(ilmos(k),jlmos(k)) = anndefct(k) 
         annsrplsrow(ilmos(k),jlmos(k)) = annsrpls(k)   
         annpcprow(ilmos(k),jlmos(k))   = annpcp(k) 
         anpotevprow(ilmos(k),jlmos(k)) = anpotevp(k)  
 130  continue
c
c     mapping the pft areal fraction in each mosaic of each
c     grid cell (farerow) to pft fraction in each grid cell (fare_cmp)
c
      do 200 i=1,nlat
      do 210 m=1,nmos
c
       do l=1,icc
        if (fcancmxrow(i,m,l) .eq. 1.0) then
         fare_cmp(i,l)    = farerow(i,m)
         nppveg_cmp(i,l)  = nppvegrow(i,m,l)
         geremort_cmp(i,l)= geremortrow(i,m,l)
         intrmort_cmp(i,l)= intrmortrow(i,m,l)
         gleafmas_cmp(i,l)= gleafmasrow(i,m,l)
         bleafmas_cmp(i,l)= bleafmasrow(i,m,l)
         stemmass_cmp(i,l)= stemmassrow(i,m,l)
         rootmass_cmp(i,l)= rootmassrow(i,m,l)
         pftexist_cmp(i,l)= pftexistrow(i,m,l)
         lambda_cmp(i,l)  = lambdarow(i,m,l)
         todfrac_cmp(i,l) = todfracrow(i,m,l)
         pfcancmx_cmp(i,l)= pfcancmxrow(i,m,l)
         nfcancmx_cmp(i,l)= nfcancmxrow(i,m,l)
         bmasveg_cmp(i,l) = bmasvegrow(i,m,l)
         burnvegf_cmp(i,l) = burnvegfrow(i,m,l)
         add2allo_cmp(i,l)= add2allorow(i,m,l)
         cc_cmp(i,l)      = ccrow(i,m,l)
         mm_cmp(i,l)      = mmrow(i,m,l)

        endif
       enddo
c
       do l=1,iccp1
        if (litrmassrow(i,m,l) .gt. 0.) then
         litrmass_cmp(i,l) = litrmassrow(i,m,l)
         soilcmas_cmp(i,l) = soilcmasrow(i,m,l)
        endif
       enddo 
c
       do l=1,ican
        if (fcanmxrow(i,m,l) .eq. 1.0) then
         fcanmx_cmp(i,l)  = fcanmxrow(i,m,l)
        endif
       enddo     
c
210   continue
c
       do m=1,nmos                                    
        vgbiomas_cmp(i) = vgbiomas_cmp(i)+vgbiomasrow(i,m)*farerow(i,m)
        gavgltms_cmp(i) = gavgltms_cmp(i)+gavgltmsrow(i,m)*farerow(i,m)
        gavgscms_cmp(i) = gavgscms_cmp(i)+gavgscmsrow(i,m)*farerow(i,m) 
        lucemcom_cmp(i) = lucemcom_cmp(i)+lucemcomrow(i,m)*farerow(i,m)
        lucltrin_cmp(i) = lucltrin_cmp(i)+lucltrinrow(i,m)*farerow(i,m)
        lucsocin_cmp(i) = lucsocin_cmp(i)+lucsocinrow(i,m)*farerow(i,m)

        ta_cmp(i)       = tarow(i,m)
        precip_cmp(i)   = preciprow(i,m)
        netrad_cmp(i)   = netrad_cmp(i)+netradrow(i,m)*farerow(i,m)
        tcurm_cmp(i)    = tcurmrow(i,m)
        srpcuryr_cmp(i) = srpcuryrrow(i,m)
        dftcuryr_cmp(i) = dftcuryrrow(i,m)
c
        do mn=1,12
         tmonth_cmp(mn,i) = tmonthrow(mn,i,m)
        enddo
c
        anpcpcur_cmp(i) = anpcpcurrow(i,m)
        anpecur_cmp(i)  = anpecurrow(i,m)  
        gdd5cur_cmp(i)  = gdd5currow(i,m)  
        surmncur_cmp(i) = surmncurrow(i,m)  
        defmncur_cmp(i) = defmncurrow(i,m)  
        srplscur_cmp(i) = srplscurrow(i,m)    
        defctcur_cmp(i) = defctcurrow(i,m)  
        twarmm_cmp(i)   = twarmmrow(i,m) 
        tcoldm_cmp(i)   = tcoldmrow(i,m)  
        gdd5_cmp(i)     = gdd5row(i,m)  
        aridity_cmp(i)  = aridityrow(i,m)  
        srplsmon_cmp(i) = srplsmonrow(i,m)  
        defctmon_cmp(i) = defctmonrow(i,m)  
        anndefct_cmp(i) = anndefctrow(i,m) 
        annsrpls_cmp(i) = annsrplsrow(i,m)   
        annpcp_cmp(i)   = annpcprow(i,m) 
        anpotevp_cmp(i) = anpotevprow(i,m)
        grclarea_cmp(i) = grclarearow(i,m)
       enddo 
c
200   continue
c
      return
      end

