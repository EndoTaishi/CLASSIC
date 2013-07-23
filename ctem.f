      subroutine     ctem( fcancmx,    fsnow,     sand,      clay,  
     2                      il1,       il2,      iday,      radj, 
     4                       tcano,    tcans,    tbarc,    tbarcs,    
     5                       tbarg,   tbargs,       ta,     delzw,
     6                     ancsveg,  ancgveg, rmlcsveg,  rmlcgveg,    
     7                       zbotw,   thliqc,   thliqg,    deltat,
     8                       uwind,    vwind,  lightng,  prbfrhuc, 
     9                    extnprob,   stdaln,     tbar,    
     a                    nol2pfts, pfcancmx, nfcancmx,  lnduseon,
     b                      thicec, soildpth, spinfast,   todfrac,
     &                     compete,   netrad,   precip,   
     &                      popdin,   dofire,     isand,  faregat,
     &                      mosaic,
c
c    -------------- inputs used by ctem are above this line ---------
c
     c                    stemmass, rootmass, litrmass,  gleafmas,
     d                    bleafmas, soilcmas,    ailcg,      ailc,
     e                       zolnc, rmatctem,    rmatc,     ailcb,
     f                    flhrloss,  pandays, lfstatus,  grwtheff,
     g                    lystmmas, lyrotmas, tymaxlai,  vgbiomas,
     h                    gavgltms, gavgscms, stmhrlos,      slai, 
     i                     bmasveg, cmasvegc, colddays,  rothrlos,
     j                      fcanmx,   alvisc,   alnirc,   gavglai,
c    ------- following 5 lines are competition related variables ----
     k                       tcurm, srpcuryr, dftcuryr,inibioclim,
     l                      tmonth, anpcpcur,  anpecur,   gdd5cur,
     m                    surmncur, defmncur, srplscur,  defctcur,
     n                    geremort, intrmort,   lambda,  lyglfmas,
     o                    pftexist, twarmm,    tcoldm,       gdd5,
     1                     aridity, srplsmon, defctmon,  anndefct,
     2                    annsrpls,  annpcp,  anpotevp,   burnveg,
c
c    -------------- inputs updated by ctem are above this line ------
c
     p                        npp,       nep, hetrores,   autores,
     q                   soilresp,        rm,       rg,       nbp,
     r                     litres,    socres,      gpp, dstcemls1,
     s                   litrfall,  humiftrs,  veghght,  rootdpth,
     t                        rml,       rms,      rmr,  tltrleaf,
     u                   tltrstem,  tltrroot, leaflitr,  roottemp,
     v                    afrleaf,   afrstem,  afrroot,  wtstatus,
     w                   ltstatus,  burnfrac, probfire,  lucemcom,
     x                   lucltrin,  lucsocin,   nppveg,  grclarea,
     y                   dstcemls3, paicgat,  slaicgat,    
     z                    emit_co2, emit_co,  emit_ch4, emit_nmhc,
     1                    emit_h2,  emit_nox, emit_n2o, emit_pm25,
     2                    emit_tpm, emit_tc,  emit_oc,  emit_bc,
     3                         cc,       mm,
     4                      rmlveg,  rmsveg,   rmrveg,    rgveg,
     5                vgbiomas_veg,  gppveg,   nepveg,   nbpveg,
     6                  hetrsveg,autoresveg, ltresveg, scresveg,
     7                 nml,    ilmos, jlmos)
c
c    ---------------- outputs are listed above this line ------------ 
c
c
c             Canadian Terrestrial Ecosystem Model (CTEM) 
C             Main Ctem Subroutine Compatible With CLASS 
C
c     Dec 6   2012   Make it so competition and luc can function in both
c     J. Melton      composite and mosaic modes.
c
c     sep 25  2012   Add competition_map and competition_unmap
c     Y. Peng   
c
c     Sep 12  2012
c     J. Melton     Add in bottom limit to bleafmas to ensure it does not
c                   slowly decay to infintesimaly small number.
c
c     Aug 23  2012
c     J. Melton     Pass in isand to ensure soil levels are properly
c                   labelled as bedrock if assigned so in classb
c
c     Jan 10  2012 
c     Yiran         Re-test bioclim and existence for competition
c
c     19  Sep. 2001 - This is the main terrestrial carbon model subroutine
c     V. Arora        
c                     all primary ctem subroutines are called from here,
c                     except phtsyn which is called from tsolvc
c
c     08 June  2001 - Add calls to three new subroutines (bioclim,
c     V. Arora        existence, and competition) to add dynamic
c                     competition between pfts. changes are also made 
c                     to land use change (luc) subroutine and the manner
c                     in which distrubance is handles. fire now creates
c                     bare ground which is subsequently available for
c                     colonization. other changes are also made to keep
c                     everything consistent with changing vegetation
c                     fractions. 
c    -----------------------------------------------------------------

      use ctem_params,        only : kk, lon, lat, pi, earthrad, zero,
     1                               edgelat, kn,iccp1, ican, ilg, nlat,
     2                               ignd, icc, nmos, l2max
      use landuse_change,     only : luc
      use competition_scheme, only : bioclim, existence, competition

c
c     inputs
c
c     fcancmx  - max. fractional coverage of ctem's 9 pfts, but this can be
c                modified by land-use change, and competition between pfts
c     fsnow    - fraction of snow simulated by class
c     sand     - percentage sand
c     clay     - percentage clay
c     icc      - no of pfts for use by ctem, currently 9
c     ican       - no of pfts for use by class, currently 4
c     ig       - no. of soil layers, 3
c     ilg      - no. of grid cells in latitude circle
c     il1,il2  - il1=1, il2=ilg
c     iday     - day of year
c     mosaic   - true if the simulation is a mosaic, otherwise it is composite
c     radj     - latitude in radians
c     tcano    - canopy temperature for canopy over ground subarea, k
c     tcans    - canopy temperature for canopy over snow subarea
c     tbarc    - soil temperature for canopy over ground subarea, k
c     tbarcs   - soil temperature for canopy over snow subarea
c     tbarg    - soil temperature for ground subarea
c     tbargs   - soil temperature for snow over ground subarea
c     ta       - air temperatuure, k
c     ancsveg  - net photosynthetic rate for ctems 9 pfts for
c                canopy over snow subarea
c     ancgveg  - net photosynthetic rate for ctems 9 pfts for
c                canopy over ground subarea
c     rmlcsveg - leaf respiration rate for ctems 9 pfts for
c                canopy over snow subarea
c     rmlcgveg - leaf respiration rate for ctems 9 pfts for
c                canopy over ground subarea
c     delzw    - thicknesses of the 3 soil layers
c     zbotw    - bottom of soil layers
c     thliqc   - liquid mois. content of 3 soil layers, for canopy
c                over snow and canopy over ground subareas
c     thliqg   - liquid mois. content of 3 soil layers, for ground
c                and snow over ground subareas
c     deltat   - ctem time step in days
c     uwind    - u wind speed, m/s
c     vwind    - v wind speed, m/s
c     lightng  - total lightning frequency, flashes/km2.year
c     prbfrhuc - probability of fire due to human causes
c     extnprob - fire extingusinging probability
c     stdaln   - an integer telling if ctem is operated within gcm (=0)
c                or in stand alone mode (=1). this is used for fire
c                purposes. see comments just above where disturb 
c                subroutine is called.
c     tbar     - soil temperature, k
c     l2max    - max. number of level 2 ctem pfts
c     nol2pfts - number of level 2 ctem pfts
c     pfcancmx - previous year's fractional coverages of pfts
c     nfcancmx - next year's fractional coverages of pfts
c     lnduseon - logical switch to run the land use change subroutine
c                or not. 
c     thicec   - frozen mois. content of 3 soil layers, for canopy
c                over snow and canopy over ground subareas
c     soildpth - soil depth (m)
c     spinfast - spinup factor for soil carbon whose default value is
c                1. as this factor increases the soil c pool will come
c                into equilibrium faster. reasonable value for spinfast
c                is between 5 and 10. when spinfast.ne.1 then the 
c                balcar subroutine is not run.
c     todfrac  - max. fractional coverage of ctem's 9 pfts by the end
c                of the day, for use by land use subroutine
c     compete  - logical boolean telling if competition between pfts is
c                on or not
c     netrad   - daily net radiation (w/m2)
c     precip   - daily precipitation (mm/day)
c     popdin   - population density (people / km^2)
c
c     updates
c
c     stemmass - stem mass for each of the 9 ctem pfts, kg c/m2
c     rootmass - root mass for each of the 9 ctem pfts, kg c/m2
c     gleafmas - green leaf mass for each of the 9 ctem pfts, kg c/m2
c     bleafmas - brown leaf mass for each of the 9 ctem pfts, kg c/m2
c     litrmass - litter mass for each of the 9 ctem pfts + bare, kg c/m2
c     soilcmas - soil carbon mass for each of the 9 ctem pfts 
c                 + bare, kg c/m2
c     ailcg    - green lai for ctem's 9 pfts
c     ailc     - lumped lai for class' 4 pfts
c     zolnc    - lumped log of roughness length for class' 4 pfts
c     rmatctem - fraction of roots for each of ctem's 9 pfts in each
c                soil layer
c     rmatc    - fraction of roots for each of class' 4 pfts in each
c                soil layer
c     ailcb    - brown lai for ctem's 9 pfts. for now we assume only
c                grasses can have brown lai
c     flhrloss - fall or harvest loss for deciduous trees and crops,
c                respectively, kg c/m2
c     pandays  - days with positive net photosynthesis (an) for use in
c                the phenology subroutine
c     lfstatus - leaf phenology status
c     grwtheff - growth efficiency. change in biomass per year per
c                unit max. lai (kg c/m2)/(m2/m2), for use in mortality
c                subroutine
c     lystmmas - stem mass at the end of last year
c     lyrotmas - root mass at the end of last year
c     tymaxlai - this year's maximum lai
c     vgbiomas - grid averaged vegetation biomass, kg c/m2
c     gavgltms - grid averaged litter mass, kg c/m2
c     gavgscms - grid averaged soil c mass, kg c/m2
c     stmhrlos - stem harvest loss for crops, kg c/m2
c     slai     - storage/imaginary lai for phenology purposes
c     bmasveg  - total (gleaf + stem + root) biomass for each ctem pft, kg c/m2
c     cmasvegc - total canopy mass for each of the 4 class pfts. recall that 
c                class requires canopy mass as an input, and this is now 
c                provided by ctem. kg/m2.
c     colddays - cold days counter for tracking days below a certain
c                temperature threshold for ndl dcd and crop pfts.
c     rothrlos - root death as crops are harvested, kg c/m2
c     fcanmx   - fractional coverage of class' 4 pfts
c     alvisc   - visible albedo for class' 4 pfts
c     alnirc   - near ir albedo for class' 4 pfts
c     gavglai  - grid averaged green leaf area index
c
c     competition related variables
c
c     tcurm     - temperature of the current month (c)
c     srpcuryr  - water surplus for the current year
c     dftcuryr  - water deficit for the current year
c     inibioclim- switch telling if bioclimatic parameters are being
c                 initialized from scratch (false) or being initialized
c                 from some spun up values(true).
c     tmonth    - monthly temperatures
c     anpcpcur  - annual precipitation for current year (mm)
c     anpecur   - annual potential evaporation for current year (mm)
c     gdd5cur   - growing degree days above 5 c for current year
c     surmncur  - number of months with surplus water for current year
c     defmncur  - number of months with water deficit for current year
c     srplscur  - water surplus for the current month
c     defctcur  - water deficit for the current month
c
c     outputs
c                grid-averaged fluxes in u-mol co2/m2.sec
c
c     npp      - net primary productivity
c     nep      - net ecosystem productivity
c     hetrores - heterotrophic respiration
c     autores  - autotrophic respiration
c     soilresp - soil respiration. this includes root respiration
c                and respiration from litter and soil carbon pools.
c                note that soilresp is different from socres, which is
c                respiration from the soil c pool.
c     rm       - maintenance respiration
c     rg       - growth respiration
c     nbp      - net biome productivity
c     gpp      - gross primary productivity
c     litres   - litter respiration
c     socres   - soil carbon respiration
c     dstcemls1- carbon emission losses due to disturbance (fire at present)
c                from vegetation  
c     litrfall - total litter fall (from leaves, stem, and root) due
c                to all causes (mortality, turnover, and disturbance)
c     humiftrs - transfer of humidified litter from litter to soil c
c                pool
c     lucemcom - land use change (luc) related combustion emission losses,
c                u-mol co2/m2.sec 
c     lucltrin - luc related inputs to litter pool, u-mol co2/m2.sec
c     lucsocin - luc related inputs to soil c pool, u-mol co2/m2.sec
c     nppveg   - npp for individual pfts,  u-mol co2/m2.sec
c     grclarea - area of the grid cell, km^2
c     dstcemls3- carbon emission losses due to disturbance (fire at present)
c                from litter pool

c     twarmm    - temperature of the warmest month (c)
c     tcoldm    - temperature of the coldest month (c)
c     gdd5      - growing degree days above 5 c
c     aridity   - aridity index, ratio of potential evaporation to precipitation
c     srplsmon  - number of months in a year with surplus water i.e. precipitation more than potential evaporation
c     defctmon  - number of months in a year with water deficit i.e.precipitation less than potential evaporation
c     anndefct  - annual water deficit (mm) 
c     annsrpls  - annual water surplus (mm)
c     annpcp    - annual precipitation (mm)
c     anpotevp  - annual potential evaporation (mm)
c
c                other quantities
c
c     veghght  - vegetation height (meters)
c     rootdpth - 99% soil rooting depth (meters)
c                both veghght & rootdpth can be used as diagnostics to see
c                how vegetation grows above and below ground, respectively
c     rml      - leaf maintenance respiration (u-mol co2/m2.sec)
c     rms      - stem maintenance respiration (u-mol co2/m2.sec)
c     rmr      - root maintenance respiration (u-mol co2/m2.sec)
c     tltrleaf - total leaf litter fall rate (u-mol co2/m2.sec)
c     tltrstem - total stem litter fall rate (u-mol co2/m2.sec)
c     tltrroot - total root litter fall rate (u-mol co2/m2.sec)
c     leaflitr - leaf litter fall rate (u-mol co2/m2.sec). this leaf litter
c                does not include litter generated due to mortality/fire
c     roottemp - root temperature, k
c     afrleaf  - allocation fraction for leaves
c     afrstem  - allocation fraction for stem
c     afrroot  - allocation fraction for root
c     wtstatus - soil water status used for calculating allocation fractions
c     ltstatus - light status used for calculating allocation fractions
c     burnfrac - areal fraction burned due to fire for every grid cell (%)
c     probfire - probability of fire for every grid cell
c
c     emitted compounds from biomass burning in g of compound
c
c      emit_co2  - carbon dioxide
c      emit_co   - carbon monoxide
c      emit_ch4  - methane
c      emit_nmhc - non-methane hydrocarbons
c      emit_h2   - hydrogen gas
c      emit_nox  - nitrogen oxides
c      emit_n2o  - nitrous oxide
c      emit_pm25 - particulate matter less than 2.5 um in diameter
c      emit_tpm  - total particulate matter
c      emit_tc   - total carbon
c      emit_oc   - organic carbon
c      emit_bc   - black carbon
c      dofire    - boolean, if true allow fire, if false no fire.
c
c     competition related variables
c
c     twarmm    - temperature of the warmest month (c)
c     tcoldm    - temperature of the coldest month (c)
c     gdd5      - growing degree days above 5 c
c     aridity   - aridity index, ratio of potential evaporation to
c                 precipitation
c     srplsmon  - number of months in a year with surplus water i.e.
c                 precipitation more than potential evaporation
c     defctmon  - number of months in a year with water deficit i.e.
c                 precipitation less than potential evaporation
c     anndefct  - annual water deficit (mm)
c     annsrpls  - annual water surplus (mm)
c     annpcp    - annual precipitation (mm)
c     anpotevp  - annual potential evaporation (mm)
c     burnveg   - areas burned, km^2, for 9 ctem pfts
c
      implicit none
c
      logical   lnduseon,  dofire, do_mortality, mosaic

      integer      il1,       il2,
     1           iday,        i,        j,        k,    stdaln,    lath,
     2         icount,        n,        m,  sort(icc),
     3   nol2pfts(ican),       k1,       k2,            spinfast,
     4           nml,    ilmos(ilg), jlmos(ilg)
c
      integer       pandays(ilg,icc), curlatno(ilg),    colddays(ilg,2),
     1             lfstatus(ilg,icc), isand(ilg,ignd)                 
c
      real fsnow(ilg),  sand(ilg,ignd), clay(ilg,ignd),thliqc(ilg,ignd),
     1     tcano(ilg), tcans(ilg),tbarc(ilg,ignd),rmatc(ilg,ican,ignd),
     2  zbotw(ilg,ignd),      rml(ilg),   gpp(ilg),   fcancmx(ilg,icc),
     3 tbarcs(ilg,ignd),tbarg(ilg,ignd),tbargs(ilg,ignd),
     4   radj(ilg),   ta(ilg), deltat, delzw(ilg,ignd),thliqg(ilg,ignd),
     5   tbar(ilg,ignd),thicec(ilg,ignd), soildpth(ilg),todfrac(ilg,icc)

c
      real fare_cmp(nlat,icc), nppveg_cmp(nlat,icc),
     1     geremort_cmp(nlat,icc),intrmort_cmp(nlat,icc),
     2     gleafmas_cmp(nlat,icc),bleafmas_cmp(nlat,icc),
     3     stemmass_cmp(nlat,icc),rootmass_cmp(nlat,icc),
     4     litrmass_cmp(nlat,iccp1),soilcmas_cmp(nlat,iccp1),
     5     lambda_cmp(nlat,icc),
     6     bmasveg_cmp(nlat,icc),   burnveg_cmp(nlat,icc),
     7     add2allo_cmp(nlat,icc),  cc_cmp(nlat,icc),mm_cmp(nlat,icc),
     8     fcanmx_cmp(nlat,ican),     
     9     vgbiomas_cmp(nlat),    
     1     gavgltms_cmp(nlat),      gavgscms_cmp(nlat),
     2     yesfrac_mos(nlat,icc),   todfrac_cmp(nlat),
     3     pfcancmx_cmp(nlat,icc),  nfcancmx_cmp(nlat,icc)
c
      integer surmncur_cmp(nlat), defmncur_cmp(nlat)

      logical pftexist_cmp(nlat,icc)
c
      real vgbiomasrow(nlat,nmos),      netradrow(nlat,nmos),
     1     gavgltmsrow(nlat,nmos),      gavgscmsrow(nlat,nmos)
c
      real ta_cmp(nlat),       precip_cmp(nlat),  netrad_cmp(nlat), 
     1     tcurm_cmp(nlat),    srpcuryr_cmp(nlat),dftcuryr_cmp(nlat),
     2     tmonth_cmp(12,nlat),anpcpcur_cmp(nlat),anpecur_cmp(nlat), 
     3     gdd5cur_cmp(nlat),  
     4     srplscur_cmp(nlat), defctcur_cmp(nlat),twarmm_cmp(nlat), 
     5     tcoldm_cmp(nlat),   gdd5_cmp(nlat),    aridity_cmp(nlat),
     6     srplsmon_cmp(nlat), defctmon_cmp(nlat),anndefct_cmp(nlat),
     7     annsrpls_cmp(nlat), annpcp_cmp(nlat),  anpotevp_cmp(nlat),
     8     lucemcom_cmp(nlat),  lucltrin_cmp(nlat), lucsocin_cmp(nlat)
c
      real  stemmass(ilg,icc),   rootmass(ilg,icc), litrmass(ilg,iccp1),
     1      gleafmas(ilg,icc),   bleafmas(ilg,icc), soilcmas(ilg,iccp1),
     2       ancsveg(ilg,icc),    ancgveg(ilg,icc),   rmlcsveg(ilg,icc),
     3      rmlcgveg(ilg,icc),      ailcg(ilg,icc),     ailc(ilg,ican),
     4   rmatctem(ilg,icc,ignd),       zolnc(ilg,ican),  ailcb(ilg,icc),
     5          vgbiomas(ilg),       gavgltms(ilg),       gavgscms(ilg),
     6          slai(ilg,icc),    bmasveg(ilg,icc),  cmasvegc(ilg,ican),
     7       veghght(ilg,icc),   rootdpth(ilg,icc),   gppcsveg(ilg,icc),
     8      gppcgveg(ilg,icc),   pfcancmx(ilg,icc),    fcanmx(ilg,ican),
     9      nfcancmx(ilg,icc),      alvisc(ilg,ican),  alnirc(ilg,ican),
     a           gavglai(ilg),    yesfrac_comp(ilg,icc)
c
      real   npp(ilg),      nep(ilg),  hetrores(ilg),      autores(ilg),
     1  soilresp(ilg),       rm(ilg),        rg(ilg),          nbp(ilg),
     2 dstcemls1(ilg), litrfall(ilg),  humiftrs(ilg),     galtcels(ilg),
     3 dstcemls2(ilg), lucemcom(ilg),  lucltrin(ilg),     lucsocin(ilg),
     4 dstcemls3(ilg)
c
      real    fc(ilg),       fg(ilg),       fcs(ilg),         fgs(ilg),
     1  fcans(ilg,ican),  fcan(ilg,ican),         rms(ilg),
     2       rmr(ilg),  grescoef(kk),    litres(ilg),      socres(ilg),
     3   humicfac(kk),        term,       laimin(kk),       laimax(kk)
c
      real pglfmass(ilg,icc),   pblfmass(ilg,icc),   pstemass(ilg,icc),
     1     protmass(ilg,icc), plitmass(ilg,iccp1), psocmass(ilg,iccp1),
     2         pvgbioms(ilg),       pgavltms(ilg),       pgavscms(ilg)
c
      real   fcancs(ilg,icc),      fcanc(ilg,icc),   rmscgveg(ilg,icc),
     1     rmscsveg(ilg,icc),   rmrcgveg(ilg,icc),   rmrcsveg(ilg,icc),
     2       rmsveg(ilg,icc),     rmrveg(ilg,icc),      anveg(ilg,icc),
     3       rmlveg(ilg,icc),     gppveg(ilg,icc),     nppveg(ilg,icc),
     4        rgveg(ilg,icc),      rmveg(ilg,icc),     nepveg(ilg,icc),
     5     rttempcs(ilg,icc),   rttempcg(ilg,icc),     nbpveg(ilg,icc),
     6     pheanveg(ilg,icc),   pancsveg(ilg,icc),   pancgveg(ilg,icc)
c
      real ltrsvgcs(ilg,icc),   ltrsvgcg(ilg,icc),   scrsvgcs(ilg,icc),
     1     scrsvgcg(ilg,icc), ltresveg(ilg,iccp1), scresveg(ilg,iccp1),
     2          ltrsbrg(ilg),        scrsbrg(ilg),       ltrsbrgs(ilg),
     3         scrsbrgs(ilg), hetrsveg(ilg,iccp1), humtrsvg(ilg,iccp1),
     4   soilrsvg(ilg,iccp1), autoresveg(ilg,icc)             
c
      real ltrestep(ilg,iccp1),screstep(ilg,iccp1), hutrstep(ilg,iccp1) 
c
      real roottemp(ilg,icc),     tbarccs(ilg,ignd),leaflitr(ilg,icc),
     1       fieldsm(ilg,ignd),   flhrloss(ilg,icc),   wiltsm(ilg,ignd)
c
      real rootlitr(ilg,icc),   stemlitr(ilg,icc),   stmhrlos(ilg,icc),
     1     rothrlos(ilg,icc)
c
      real  afrleaf(ilg,icc),    afrstem(ilg,icc),    afrroot(ilg,icc),
     1     wtstatus(ilg,icc),   ltstatus(ilg,icc)
c
      real nppvgstp(ilg,icc),   rmlvgstp(ilg,icc),   rmsvgstp(ilg,icc),
     1     rmrvgstp(ilg,icc),   gppvgstp(ilg,icc),   ntchlveg(ilg,icc),
     2     ntchsveg(ilg,icc),   ntchrveg(ilg,icc)
c
      real grwtheff(ilg,icc),   lystmmas(ilg,icc),   lyrotmas(ilg,icc), 
     1     tymaxlai(ilg,icc),   stemltrm(ilg,icc),   rootltrm(ilg,icc), 
     2     glealtrm(ilg,icc),   geremort(ilg,icc),   intrmort(ilg,icc)
c
      real      currlat(ilg),            wl(lat),
     1             radl(lat),          wossl(lat),             sl(lat),
     2               cl(lat),             ml(ilg),       grclarea
c
      real        uwind(ilg),          vwind(ilg),        lightng(ilg),
     1         prbfrhuc(ilg),       extnprob(ilg),   pftareab(ilg,icc),
     2     stemltdt(ilg,icc),   rootltdt(ilg,icc),   glfltrdt(ilg,icc),
     3     blfltrdt(ilg,icc),   glcaemls(ilg,icc),   blcaemls(ilg,icc),
     4     rtcaemls(ilg,icc),   stcaemls(ilg,icc),   ltrcemls(ilg,icc),
     5     pftareaa(ilg,icc),       burnfrac(ilg),   dscemlv1(ilg,icc),
     6     dscemlv2(ilg,icc),       probfire(ilg),    burnveg(ilg,icc)
c
      real     emit_co2(ilg),        emit_co(ilg),       emit_ch4(ilg),
     1        emit_nmhc(ilg),        emit_h2(ilg),       emit_nox(ilg),
     2         emit_n2o(ilg),      emit_pm25(ilg),       emit_tpm(ilg),
     3          emit_tc(ilg),        emit_oc(ilg),        emit_bc(ilg)
c
      real tltrleaf(ilg,icc),   tltrstem(ilg,icc),   tltrroot(ilg,icc),
     1                popdin
c
      real  faregat(ilg), paicgat(ilg,ican),slaicgat(ilg,ican)  
c
      real  vgbiomas_veg(ilg,icc)
c  
      real       precip(ilg),         netrad(ilg),         tcurm(ilg),
     1           annpcp(ilg),       anpotevp(ilg),        twarmm(ilg),
     2           tcoldm(ilg),           gdd5(ilg),       aridity(ilg),
     3         srplsmon(ilg),       defctmon(ilg),     tmonth(12,ilg),
     4         anpcpcur(ilg),        anpecur(ilg),       gdd5cur(ilg),
     5         srplscur(ilg), 
     6         defctcur(ilg),       srpcuryr(ilg),      dftcuryr(ilg),
     7         anndefct(ilg),       annsrpls(ilg)
c     
      real     barefrac(ilg),       pbarefrc(ilg),           tolrance,
     1       lambda(ilg,icc),   add2allo(ilg,icc),  lyglfmas(ilg,icc),
     2     expbalvg(ilg,icc),       expnbaln(ilg), ltrflcom(ilg,iccp1),
     3             lambdamax,         cc(ilg,icc),         mm(ilg,icc)
c
      integer   surmncur(ilg),       defmncur(ilg)
c
      logical compete, inibioclim, pftexist(ilg,icc)
c
c     ---------------------------------------------------------------
c                     constants and parameters 
c
c     note the structure of vectors which clearly shows the class
c     pfts (along rows) and ctem sub-pfts (along columns)
c
c     needle leaf |  evg       dcd       ---
c     broad leaf  |  evg   dcd-cld   dcd-dry
c     crops       |   c3        c4       ---
c     grasses     |   c3        c4       ---
c
c     growth respiration coefficient
      data grescoef/0.15, 0.15, 0.00,
     &              0.15, 0.15, 0.15,
     &              0.15, 0.15, 0.00,
     &              0.15, 0.15, 0.00/
c
c     humification factor - used for transferring carbon from litter into 
c     soil c pool
      data humicfac/0.42, 0.42, 0.00,
     &              0.53, 0.48, 0.48,
     &              0.42, 0.42, 0.00,
     &              0.42, 0.42, 0.00/
c
c     minimum lai below which a pft doesn't expand
      data laimin/1.0, 1.0, 0.0,
     &            2.0, 1.5, 1.0,
     &            1.0, 1.0, 0.0,
     &            0.5, 0.5, 0.0/
c
c     maximum lai above which a pft always expands and lambdamax fraction
c     of npp is used for expansion
      data laimax/3.0, 3.0, 0.0,
     &            8.0, 8.0, 5.0,
     &            8.0, 8.0, 0.0,
     &            4.0, 4.0, 0.0/
c
c     max. fraction of npp that is allocated for reproduction/colonization
      data lambdamax/0.10/
c
c     -----------------------------------------------------------------
c     find area of the gcm grid cells. this is needed for land use change
c     and disturbance subroutines
c
      if(stdaln.eq.0)then         ! i.e. when operated in a gcm mode 
c
        do 50 i = il1, il2
          currlat(i)=radj(i)*180.0/pi                             
          curlatno(i)=0
50     continue
c
c       find current latitude number
        do 60 k = 1, lat
          do 61 i = il1, il2
            if(currlat(i).ge.edgelat(k).and.
     &      currlat(i).lt.edgelat(k+1))then   
              curlatno(i)=k
            endif
61        continue
60      continue
c
        do 70 i = il1, il2
          if(curlatno(i).eq.0)then
            write(6,2000)i
2000        format('cannot find current latitude no. for i = ',i3)  
            call xit ('ctem',-5)
          endif
70      continue
c
        lath = lat/2
        call gaussg(lath,sl,wl,cl,radl,wossl)
        call trigl(lath,sl,wl,cl,radl,wossl)
c
c       wl contains zonal weights, lets find meridional weights
c
        do 80 i = il1,il2
          ml(i) = 1.0/real(lon)
80      continue 
c
C        do 81 i = il1, il2
          grclarea = 4.0*pi*(earthrad**2)*wl(curlatno(i))*ml(i)
     &                   *faregat(i)/2.0  !km^2, faregat is areal fraction of each mosaic
c         dividing by 2.0 because wl(1 to lat) add to 2.0 not 1.0
C81      continue  
c
      else if(stdaln.eq.1)then    ! i.e. when operated at point scale
c
        do 90 j = 1, icc
C          do 91 i = il1, il2
            grclarea=1.0
C91        continue
90      continue
c
      endif
c
c     ---------------------------------------------------------------
c

c     generate the sort index for correspondence between 9 pfts and the
c     12 values in the parameter vectors
c
      icount=0
      do 95 j = 1, ican
        do 96 m = 1, nol2pfts(j)
          n = (j-1)*l2max + m
          icount = icount + 1
          sort(icount)=n
 96    continue
 95   continue
c
c     ---------------------------------------------------------------
c
      if(compete .or. lnduseon)then

c      land use change and competition for mosaics needs mapping and
c      unmapping of the pfts. composite does not require these extra steps.

       if (mosaic) then
c
c       check if number of mosaics is equal to the number of pfts plus one
c       bare, e.g., nmos=iccp1
c
        if (nmos.ne.iccp1) then 
         write(*,2050) 'number of mosaics, nmos= ',nmos,
     &                 ' is not equal to the number of pfts plus',
     &                 ' one bare, iccp1= ',iccp1
         write(*,2051) 'competition works properly only when all pfts',
     &                 ' and bare are considered.                    '
         call xit ('ctem',-11)
        endif 
2050    format(a25,i2,a40,a18,i2,a1)
2051    format(a45,a40)
c
c       check for fcancmx(i,j). this should be either 0 or 1 for competition
c       to work.
c
        do j=1, icc
         do i=il1, il2
          if(fcancmx(i,j).ne.1.0 .and. fcancmx(i,j).ne.0.0) then
           write(*,2100) 
     &                'mosaic ',i,' has pft fraction: ', 
     &                'fcancmx(',i,',',j,')=',fcancmx(i,j)
           write(*,2101) 
     &                'mosaic competition and luc work only when ',
     &                'each mosaic is 100% occupied with one pft'
           call xit ('ctem',-12)
          endif
         enddo
        enddo
2100    format(a7,i2,a19,a8,i2,a1,i2,a2,f8.3)
2101    format(a40,a40)
c   
c       competition_map scatters and maps the array with indices 
c       of (ilg,icc) to (nlat,icc) for preparation for competition
c  
          call competition_map(    nml,    ilmos,   jlmos,  
     b                         faregat,   fcancmx,  nppveg,  geremort,
     c                         intrmort, gleafmas, bleafmas, stemmass,
     d                         rootmass, litrmass, soilcmas,
     e                         pftexist,   lambda,  bmasveg,  burnveg,
     f                         add2allo,       cc,       mm,   fcanmx,
     g                         vgbiomas, gavgltms, gavgscms,
     h                               ta,   precip,   netrad,    tcurm,
     i                         srpcuryr, dftcuryr,   tmonth, anpcpcur, 
     j                          anpecur,  gdd5cur, surmncur, defmncur,
     k                         srplscur, defctcur,   twarmm,   tcoldm,
     l                             gdd5,  aridity, srplsmon, defctmon,
     m                         anndefct, annsrpls,   annpcp, anpotevp,
     &                         lucemcom, lucltrin, lucsocin, pfcancmx,
     &                         nfcancmx,
c    ------------------- inputs above this line ---------------------
     n                        netradrow,
c    ------------------- intermediate and saved above this line -----
     o                         fare_cmp,    nppveg_cmp,  geremort_cmp,
     p                     intrmort_cmp,  gleafmas_cmp,  bleafmas_cmp,
     q                     stemmass_cmp,  rootmass_cmp,  litrmass_cmp,
     r                     soilcmas_cmp,  pftexist_cmp,    lambda_cmp,
     s                      bmasveg_cmp,   burnveg_cmp,  add2allo_cmp,
     t                           cc_cmp,        mm_cmp,    fcanmx_cmp,
     u                     vgbiomas_cmp,
     v                     gavgltms_cmp,  gavgscms_cmp,
     w                           ta_cmp,    precip_cmp,    netrad_cmp, 
     x                        tcurm_cmp,  srpcuryr_cmp,  dftcuryr_cmp,
     y                       tmonth_cmp,  anpcpcur_cmp,   anpecur_cmp, 
     z                      gdd5cur_cmp,  surmncur_cmp,  defmncur_cmp,
     1                     srplscur_cmp,  defctcur_cmp,    twarmm_cmp, 
     2                       tcoldm_cmp,      gdd5_cmp,   aridity_cmp,
     3                     srplsmon_cmp,  defctmon_cmp,  anndefct_cmp,
     4                     annsrpls_cmp,    annpcp_cmp,  anpotevp_cmp,
     5                     lucemcom_cmp,  lucltrin_cmp,  lucsocin_cmp,
     6                     pfcancmx_cmp,   nfcancmx_cmp )
c    ------------------- outputs above this line --------------------
c   
        if (compete) then

c        calculate bioclimatic parameters for estimating pfts existence
c
         call  bioclim (iday,       ta_cmp,    precip_cmp,  netrad_cmp,
     1                    1,         nlat,          nlat,
     2            tcurm_cmp, srpcuryr_cmp,  dftcuryr_cmp,  inibioclim,
     3           tmonth_cmp, anpcpcur_cmp,   anpecur_cmp, gdd5cur_cmp,
     4         surmncur_cmp, defmncur_cmp,  srplscur_cmp,defctcur_cmp,
     5           twarmm_cmp,   tcoldm_cmp,      gdd5_cmp, aridity_cmp,
     6         srplsmon_cmp, defctmon_cmp, anndefct_cmp, annsrpls_cmp,
     7           annpcp_cmp, anpotevp_cmp)

       if (inibioclim) then
c
c        if first day of year then based on updated bioclimatic parameters
c        find if pfts can exist or not
c
          call existence(iday,            1,         nlat,         nlat,
     1                   sort,     nol2pfts,      
     2             twarmm_cmp,   tcoldm_cmp,     gdd5_cmp,  aridity_cmp,
     3           srplsmon_cmp, defctmon_cmp, anndefct_cmp, annsrpls_cmp,
     4             annpcp_cmp, anpotevp_cmp, pftexist_cmp)
c     
c
c        call competition subroutine which on the basis of previous day's
c        npp estimates changes in fractional coverage of pfts
c
         call competition (iday,          1,          nlat,        nlat,
     1               nol2pfts,   nppveg_cmp,
     2           pftexist_cmp, geremort_cmp, intrmort_cmp,
     3           gleafmas_cmp, bleafmas_cmp, stemmass_cmp, rootmass_cmp,
     4           litrmass_cmp, soilcmas_cmp, grclarea,   lambda_cmp,
     5            bmasveg_cmp,  burnveg_cmp,         sort,
c
c    ------------------- inputs above this line -------------------
c
     6               fare_cmp,   fcanmx_cmp, vgbiomas_cmp, gavgltms_cmp,
     7           gavgscms_cmp,
c
c    ------------------- updates above this line ------------------
c
     8           add2allo_cmp,      cc_cmp,      mm_cmp)


         else  ! It is the first year of a run that you do not have the 
              ! climatological means already in the CTM file. After one
              ! year inibioclim is set to true and the climatological means
              ! are used from the first year.

           do i = il1, nlat
             do j = 1, icc
              add2allo_cmp(i,j)=0.0
             enddo
           enddo

         end if !inibioclim
c
c    ------------------- outputs above this line ------------------
c      
        else !compete = false

           do i = il1, nlat
             do j = 1, icc
              add2allo_cmp(i,j)=0.0
             enddo
           enddo

        endif !compete check
c     -----------------------------------------------------------------

        if(lnduseon)then

         do i = il1, nlat
           do j = 1, icc  
            yesfrac_mos(i,j)=fare_cmp(i,j)
           enddo
         enddo

         call luc(il1,      nlat,     nlat,     nol2pfts, 
     2           grclarea,    pfcancmx_cmp, nfcancmx_cmp,     iday,
     3           todfrac_cmp,  yesfrac_mos,   .true.,
     4           gleafmas_cmp, bleafmas_cmp, stemmass_cmp, rootmass_cmp,
     5           litrmass_cmp, soilcmas_cmp, vgbiomas_cmp, gavgltms_cmp,
     6           gavgscms_cmp,     fare_cmp,   fcanmx_cmp,
     7           lucemcom_cmp, lucltrin_cmp, lucsocin_cmp)

        endif !lnduseon check

c     -----------------------------------------------------------------
c       competition_unmap unmaps and gathers the array with  
c       indices (nlat,icc) back to (ilg,icc) after competition is done 
c
        call competition_unmap( nml,      ilmos,    jlmos,   nol2pfts,
     b                           fare_cmp,   nppveg_cmp, geremort_cmp,
     c                       intrmort_cmp, gleafmas_cmp, bleafmas_cmp,
     d                       stemmass_cmp, rootmass_cmp, litrmass_cmp,
     e                       soilcmas_cmp, pftexist_cmp,   lambda_cmp,
     f                        bmasveg_cmp,  burnveg_cmp, add2allo_cmp,
     g                             cc_cmp,       mm_cmp,   fcanmx_cmp,
     h                       vgbiomas_cmp,
     i                       gavgltms_cmp, gavgscms_cmp,
     j                             ta_cmp,   precip_cmp,   netrad_cmp, 
     k                          tcurm_cmp, srpcuryr_cmp, dftcuryr_cmp,
     l                         tmonth_cmp, anpcpcur_cmp,  anpecur_cmp, 
     m                        gdd5cur_cmp, surmncur_cmp, defmncur_cmp,
     n                       srplscur_cmp, defctcur_cmp,   twarmm_cmp, 
     o                         tcoldm_cmp,     gdd5_cmp,  aridity_cmp,
     p                       srplsmon_cmp, defctmon_cmp, anndefct_cmp,
     q                       annsrpls_cmp,   annpcp_cmp, anpotevp_cmp,
     &                     lucemcom_cmp,  lucltrin_cmp,  lucsocin_cmp,
     &                     pfcancmx_cmp,   nfcancmx_cmp,
c
c    ------------------- inputs above this line ---------------------
c
     r                            netradrow,
c
c    ------------------- saved for intermediate above this line -----
c
     s                        faregat,  fcancmx,    nppveg, geremort,  
     t                       intrmort, gleafmas,  bleafmas, stemmass,
     u                       rootmass, litrmass,  soilcmas,
     v                        pftexist,  lambda,   bmasveg,  burnveg,
     w                        add2allo,      cc,        mm,   fcanmx,
     x                        vgbiomas, gavgltms, gavgscms,
     y                     ta,  precip,   netrad,    tcurm, srpcuryr,
     z                        dftcuryr ,  tmonth, anpcpcur,  anpecur, 
     1                         gdd5cur, surmncur, defmncur, srplscur,  
     2                        defctcur,   twarmm,   tcoldm,     gdd5, 
     3                         aridity, srplsmon, defctmon, anndefct,
     4                        annsrpls,   annpcp, anpotevp,
     5                         lucemcom, lucltrin, lucsocin, pfcancmx,
     6                         nfcancmx )
c    ------------------- updates above this line --------------------
c
      else !composite
c
       if (compete) then

c       calculate bioclimatic parameters for estimating pfts existence
c
        call  bioclim (iday,       ta,    precip,  netrad,
     1                    1,     il2,    ilg,
     2                 tcurm, srpcuryr,  dftcuryr,  inibioclim,
     3                 tmonth, anpcpcur,   anpecur, gdd5cur,
     4                 surmncur, defmncur,  srplscur,defctcur,
     5                 twarmm,   tcoldm,      gdd5, aridity,
     6                 srplsmon, defctmon, anndefct, annsrpls,
     7                 annpcp, anpotevp )
c
        if (inibioclim) then

c       if first day of year then based on updated bioclimatic parameters
c       find if pfts can exist or not
c
        call existence(iday,            1,         il2,         ilg,
     1                     sort,     nol2pfts,        
     2                   twarmm,   tcoldm,     gdd5,  aridity,
     3                   srplsmon, defctmon, anndefct, annsrpls,
     4                   annpcp, anpotevp, pftexist )
c     
c       call competition subroutine which on the basis of previous day's
c       npp estimates changes in fractional coverage of pfts
c
        call competition (iday,          1,          il2,         ilg,
     1                    nol2pfts,       nppveg,
     2                    pftexist, geremort, intrmort,
     3                    gleafmas, bleafmas, stemmass, rootmass,
     4                    litrmass, soilcmas, grclarea,   lambda,
     5                    bmasveg,  burnveg,        sort,
c
c    ------------------- inputs above this line -------------------
c
     6                    fcancmx,   fcanmx, vgbiomas, gavgltms,
     7                    gavgscms,
c
c    ------------------- updates above this line ------------------
c
     8                    add2allo,      cc,      mm)
c
c    ------------------- outputs above this line ------------------
c
        else  ! It is the first year of a run that you do not have the 
              ! climatological means already in the CTM file. After one
              ! year inibioclim is set to true and the climatological means
              ! are used from the first year.
          do j = 1, icc
            do i = il1, il2
              add2allo(i,j)=0.0
            enddo
          enddo
        end if

       else  ! competition is not on

        do j = 1, icc
          do i = il1, il2
            add2allo(i,j)=0.0
          enddo
        enddo

       endif  ! if (compete)
c
c     -----------------------------------------------------------------
c
c      if landuse is on, then implelement luc, change fractional coverages, 
c      move biomasses around, and estimate luc related combustion emission 
c      losses.
c
        if(lnduseon)then
         
         do j = 1, icc
           do i = il1, il2  
             yesfrac_comp(i,j)=fcancmx(i,j)
           enddo
         enddo

         call luc(    il1,      il2,   ilg,  nol2pfts, 
     2                  grclarea, pfcancmx, nfcancmx,     iday,
     3                   todfrac,yesfrac_comp,.true.,
     4                  gleafmas, bleafmas, stemmass, rootmass,
     5                  litrmass, soilcmas, vgbiomas, gavgltms,
     6                  gavgscms,  fcancmx,   fcanmx,
     7                  lucemcom, lucltrin, lucsocin)

        endif !lnduseon
c
       endif ! mosaic vs. composite

       else  !not compete/lnduseon

        do j = 1, icc
          do i = il1, il2
            add2allo(i,j)=0.0
          enddo
        enddo

      endif !compete/lnduseon

c     ---------------------------------------------------------------
c
c     find pft areas before (these are required for disturb subroutine)
c
      if(stdaln.eq.0)then         ! i.e. when operated in a gcm mode 
        do 82 j = 1, icc
          do  83 i = il1, il2
            pftareab(i,j)=grclarea*fcancmx(i,j)   ! area in km^2
83        continue
82      continue
      else if(stdaln.eq.1)then    ! i.e. when operated at point scale
        do 92 j = 1, icc
          do 93 i = il1, il2
            grclarea=0.01
            pftareab(i,j)=0.01*fcancmx(i,j)      ! 0.01 km2 = 1 hectare
93        continue
92      continue
      endif
c
c     initialize required arrays to zero
c
      do 100 i = il1, il2
        rms(i) = 0.0         !grid ave. stem maintenance respiration
        rmr(i) = 0.0         !grid ave. root maintenance respiration
        rml(i) = 0.0         !grid ave. leaf maintenance respiration
        rm(i) = 0.0          !grid ave. total maintenance respiration
        rg(i) = 0.0          !grid ave. growth respiration
        npp(i) = 0.0         !grid ave. net primary productivity
        gpp(i) = 0.0         !grid ave. gross primary productivity
        nep(i)=0.0           !grid ave. net ecosystem productivity
        nbp(i)=0.0           !grid ave. net biome productivity
c
        litres(i)=0.0        !grid ave. litter respiration
        socres(i)=0.0        !grid ave. soil carbon respiration
c
        hetrores(i)=0.0      !grid ave. heterotrophic respiration
        autores(i)=0.0       !grid ave. autotrophic respiration
        soilresp(i)=0.0      !grid ave. soil respiration
        humiftrs(i)=0.0      !grid ave. humification rate
        dstcemls1(i)=0.0     !grid ave. carbon emission losses due to disturbance, vegetation
        dstcemls2(i)=0.0     !grid ave. carbon emission losses due to disturbance, total
        dstcemls3(i)=0.0     !grid ave. carbon emission losses due to disturbance, litter
        galtcels(i)=0.0      !grid ave. litter fire emission losses (redundant, same as dstcemls3)
c
        fc(i)=0.0            !fraction of canopy over ground subarea 
        fcs(i)=0.0           !fraction of canopy over snow subarea
        fg(i)=0.0            !fraction of bare ground subarea 
        fgs(i)=0.0           !fraction of snow over ground subarea
c
        tbarccs(i,1)=0.0     !avg. soil temperature over canopy over snow
        tbarccs(i,2)=0.0     !and canopy over ground subareas.
        tbarccs(i,3)=0.0     
c
c                              over bare fraction of the grid cell
        screstep(i,iccp1)=0.0  !soil c respiration in kg c/m2 over the time step
        ltrestep(i,iccp1)=0.0  !litter c respiration in kg c/m2 over the time step
        soilrsvg(i,iccp1)=0.0  !soil respiration over the bare fraction
        humtrsvg(i,iccp1)=0.0  !humified rate the bare fraction
c
        ltresveg(i,iccp1)=0.0  !litter respiration rate over bare fraction
        scresveg(i,iccp1)=0.0  !soil c respiration rate over bare fraction
        hetrsveg(i,iccp1)=0.0  !heterotrophic resp. rate over bare fraction
c       expnbaln  is used for competition
        expnbaln(i)=0.0        !amount of c related to spatial expansion

100   continue 
c
      do 110 j = i,icc
        do 120 i = il1, il2
          fcanc(i,j) =0.0
          fcancs(i,j)=0.0
c
          rmsveg(i,j)=0.0    !stem maintenance resp. rate for each pft
          rmrveg(i,j)=0.0    !root maintenance resp. rate for each pft
          rmlveg(i,j)=0.0    !leaf maintenance resp. rate for each pft
           rmveg(i,j)=0.0    !total maintenance resp. rate for each pft
           rgveg(i,j)=0.0    !growth resp. rate for each pft
           anveg(i,j)=0.0    !net photosynthesis rate for each pft
          pheanveg(i,j)=0.0  !net photosynthesis rate, for phenology purposes
          pancsveg(i,j)=0.0  !net photosynthesis rate, canopy over snow subarea, for phenology purposes
          pancgveg(i,j)=0.0  !net photosynthesis rate, canopy over ground subarea, for phenology purposes
c
          gppveg(i,j)=0.0    !gross primary productity for each pft
          nppveg(i,j)=0.0    !net primary productity for each pft
          nbpveg(i,j)=0.0    !net biome productity for each pft
          nepveg(i,j)=0.0    !net ecosystem productity for each pft
c
          ltresveg(i,j)=0.0  !litter respiration rate for each pft
          scresveg(i,j)=0.0  !soil c respiration rate for each pft
          hetrsveg(i,j)=0.0  !heterotrophic resp. rate for each pft
          soilrsvg(i,j)=0.0  !soil respiration rate for each pft
          humtrsvg(i,j)=0.0  !humification rate for each pft
          screstep(i,j)=0.0  !soil c respiration in kg c/m2 over the time step
          ltrestep(i,j)=0.0  !litter c respiration in kg c/m2 over the time step
          hutrstep(i,j)=0.0  !humification rate in kg c/m2 over the time step
c
          roottemp(i,j)=0.0  !root temperature
          nppvgstp(i,j)=0.0  !npp (kg c/m2) sequestered over the model time step
          gppvgstp(i,j)=0.0  !gpp (kg c/m2) sequestered over the model time step
          rmlvgstp(i,j)=0.0  !leaf maintenance resp. (kg c/m2) respired over the model time step
          rmsvgstp(i,j)=0.0  !stem maintenance resp. (kg c/m2) respired over the model time step
          rmrvgstp(i,j)=0.0  !root maintenance resp. (kg c/m2) respired over the model time step
c 
          ntchlveg(i,j)=0.0  !net change in gleaf biomass after auto. resp. & allocation
          ntchsveg(i,j)=0.0  !net change in stem biomass after auto. resp. & allocation
          ntchrveg(i,j)=0.0  !net change in root biomass after auto. resp. & allocation
c
          dscemlv1(i,j)=0.0  !total carbon emission losses (kg c/m2), mainly due to fire
          dscemlv2(i,j)=0.0  !total carbon emission losses (kg c/m2), mainly due to fire
c
          tltrleaf(i,j)=0.0  !total leaf litter
          tltrstem(i,j)=0.0  !total stem litter
          tltrroot(i,j)=0.0  !total root litter
c
          vgbiomas_veg(i,j)=0.0 !vegetation biomass for each pft
c
c         following are competition related
          lambda(i,j)=0.0    !fraction of npp that is used for horizontal expansion
          expbalvg(i,j)=0.0  !amount of c related to spatial expansion
          add2allo(i,j)=0.0
 
c
120     continue
110   continue
c
c     store green and brown leaf, stem, and root biomass, and litter and 
c     soil c pool mass in arrays. knowing initial sizes of all pools and
c     final sizes at the end of this subroutine, we check for conservation
c     of mass.
c
      do 130 j = 1, icc
        do 140 i = il1, il2
          pglfmass(i,j)=gleafmas(i,j)    !green leaf mass from last time step
          pblfmass(i,j)=bleafmas(i,j)    !brown leaf mass from last time step
          pstemass(i,j)=stemmass(i,j)    !stem mass from last time step
          protmass(i,j)=rootmass(i,j)    !root mass from last time step
          plitmass(i,j)=litrmass(i,j)    !litter mass from last time step
          psocmass(i,j)=soilcmas(i,j)    !soil c mass from last time step
140     continue
130   continue
c
      do 145 i = il1, il2
        pvgbioms(i)=vgbiomas(i)          !vegetation biomass from last time step
        vgbiomas(i)= 0.0
        pgavltms(i)=gavgltms(i)          !litter mass from last time step
        gavgltms(i)=0.0
        pgavscms(i)=gavgscms(i)          !soil c mass from last time step
        gavgscms(i)=0.0
        litrfall(i)=0.0                  !combined total litter fall rate
        gavglai (i)=0.0                  !grid averaged green lai
c
        plitmass(i,iccp1)=litrmass(i,iccp1)  !litter mass over bare fraction
        psocmass(i,iccp1)=soilcmas(i,iccp1)  !soil c mass over bare fraction
145   continue
c
c     initialization ends
c
c     find fc and fcs based on fcancmx
c
      do 150 j = 1, icc
        do 160 i = il1, il2
          fcancs(i,j) = fcancmx(i,j)*fsnow(i)
          fcanc(i,j)  = fcancmx(i,j)*(1.-fsnow(i))
          fcs(i) = fcs(i) + fcancs(i,j)
          fc(i)  = fc(i)  + fcanc(i,j)
160     continue 
150   continue 
c
      do 170 i = il1, il2
        fgs(i)=(1.0-fcs(i)-fc(i))*fsnow(i) 
        fg(i)=(1.0-fcs(i)-fc(i))*(1.0-fsnow(i)) 
170   continue
c
c     ------------------------------------------------------------------
c
c     Autotrophic respiration part starts
c
c     Leaf respiration is calculated in phtsyn subroutine, while stem
c     and root maintenance respiration are calculated here.
c
c     We treat canopy over ground and canopy over snow subareas
c     separately because stem temperature (for which we use canopy
c     temperature as a surrogate) can be different for these two
c     subareas.
c
c     Find maintenance respiration for canopy over snow sub-area
c     in umol co2/m2/sec
c
      call   mainres (fcancs,      fcs,     stemmass,   rootmass,        
     1                  il1,
     2                   il2,       ta,       tbarcs,   rmatctem,
     3                  sort, nol2pfts,        isand,
     4              rmscsveg, rmrcsveg,     rttempcs)
c
c     Find maintenance respiration for canopy over ground sub-area
c
      call   mainres ( fcanc,       fc,     stemmass,   rootmass,        
     1                   il1,
     2                   il2,       ta,        tbarc,   rmatctem,
     3                  sort, nol2pfts,        isand,
     4              rmscgveg, rmrcgveg,     rttempcg)
c
c
c     If ailcg/gleafmas is zero, i.e. real leaves are not on, then
c     make maintenance respiration and gpp from storage/imaginary lai 
c     equal to zero so that we don't use these numbers in carbon budget.
c
      do 180 j = 1, icc
        do 190 i = il1, il2
          gppcsveg(i,j)=ancsveg(i,j)+rmlcsveg(i,j)
          gppcgveg(i,j)=ancgveg(i,j)+rmlcgveg(i,j)
c
          if (lfstatus(i,j).eq.4) then
            rmlcgveg(i,j)=0.0
            rmlcsveg(i,j)=0.0
            pancsveg(i,j)=ancsveg(i,j)   ! to be used for phenology
            pancgveg(i,j)=ancgveg(i,j)   ! purposes
            ancsveg(i,j)=0.0
            ancgveg(i,j)=0.0
          else
            pancsveg(i,j)=ancsveg(i,j)   ! to be used for phenology
            pancgveg(i,j)=ancgveg(i,j)   ! purposes
            if(slai(i,j).gt.ailcg(i,j))then
             term=((1.0/kn(sort(j)))*(1.0-exp(-kn(sort(j))*ailcg(i,j))) 
     &          /(1.0/kn(sort(j)))*(1.0-exp(-kn(sort(j))* slai(i,j))))
             rmlcgveg(i,j)=rmlcgveg(i,j)*term
             rmlcsveg(i,j)=rmlcsveg(i,j)*term
            endif
          endif
190     continue
180   continue
c
c     find vegetation averaged leaf, stem, and root respiration, and
c     gpp using values from canopy over ground and canopy over snow
c     subareas
c
      do 270 j = 1, icc
        do 280 i = il1, il2
          if( (fcanc(i,j)+fcancs(i,j)).gt.zero) then
            rmsveg(i,j)= (fcanc(i,j)*rmscgveg(i,j) + 
     &        fcancs(i,j)*rmscsveg(i,j)) / ( fcanc(i,j) + fcancs(i,j))     
            rmrveg(i,j)= (fcanc(i,j)*rmrcgveg(i,j) + 
     &        fcancs(i,j)*rmrcsveg(i,j)) / ( fcanc(i,j) + fcancs(i,j))
            rmlveg(i,j)= (fcanc(i,j)*rmlcgveg(i,j) + 
     &        fcancs(i,j)*rmlcsveg(i,j)) / ( fcanc(i,j) + fcancs(i,j))
            anveg(i,j)= (fcanc(i,j)*ancgveg(i,j) + 
     &        fcancs(i,j)*ancsveg(i,j)) / ( fcanc(i,j) + fcancs(i,j))
            gppveg(i,j)= (fcanc(i,j)*gppcgveg(i,j) + 
     &        fcancs(i,j)*gppcsveg(i,j)) / ( fcanc(i,j) + fcancs(i,j))
            pheanveg(i,j)= (fcanc(i,j)*pancgveg(i,j) + 
     &        fcancs(i,j)*pancsveg(i,j)) / ( fcanc(i,j) + fcancs(i,j))
          else
            rmsveg(i,j)= 0.0
            rmrveg(i,j)= 0.0
            rmlveg(i,j)= 0.0
            anveg(i,j)= 0.0
            gppveg(i,j)= 0.0
            pheanveg(i,j)= 0.0
          endif

c         
          if(lfstatus(i,j).eq.4)then
            gppveg(i,j) = anveg(i,j) + rmlveg(i,j)
          endif
c
          rmveg(i,j)  = rmlveg(i,j) + rmrveg(i,j) + rmsveg(i,j)
          nppveg(i,j) = gppveg(i,j) - rmveg(i,j)
280     continue 
270   continue 
c
c     now that we know maintenance respiration from leaf, stem, and root,
c     and gpp, we can find growth respiration for each vegetation 
c
      do 300 j = 1, icc
        do 310 i = il1, il2
          if( nppveg(i,j).gt.zero ) then
            rgveg(i,j)=grescoef(sort(j))*nppveg(i,j)
          else
            rgveg(i,j)=0.0
          endif
          nppveg(i,j) = nppveg(i,j) - rgveg(i,j)
310     continue
300   continue
c
c     calculate grid-averaged rates of rm, rg, npp, and gpp
c
      do 320 j = 1,icc
        do 330 i = il1, il2
          rml(i)=rml(i)+fcancmx(i,j)*rmlveg(i,j)
          rms(i)=rms(i)+fcancmx(i,j)*rmsveg(i,j)
          rmr(i)=rmr(i)+fcancmx(i,j)*rmrveg(i,j)
          rm(i) =rm(i)+fcancmx(i,j)*rmveg(i,j)
          rg(i) =rg(i)+fcancmx(i,j)*rgveg(i,j)
          npp(i)=npp(i)+fcancmx(i,j)*nppveg(i,j)
          gpp(i)=gpp(i)+fcancmx(i,j)*gppveg(i,j)
          autores(i)=rg(i)+rm(i)
          autoresveg(i,j)=rmveg(i,j) + rgveg(i,j)
330     continue
320   continue
c
c
c     autotrophic respiration part ends
c
c     ------------------------------------------------------------------
c
c     heterotrophic respiration part starts
c
c     find heterotrophic respiration rates (umol co2/m2/sec) for canopy
c     over snow subarea
c
       call    hetresv ( fcancs,      fcs, litrmass, soilcmas,
     1                      il1,
     2                      il2,   tbarcs,   thliqc,     sand,
     3                     clay, rttempcs,    zbotw,     sort,
     4                     isand,
     5                 ltrsvgcs, scrsvgcs) 
c
c     find heterotrophic respiration rates for canopy over ground 
c     subarea
c
       call    hetresv (  fcanc,       fc, litrmass, soilcmas,
     1                      il1,
     2                      il2,    tbarc,   thliqc,     sand,
     3                     clay, rttempcg,    zbotw,     sort,
     4                     isand,
     5                 ltrsvgcg, scrsvgcg) 

c
c     find heterotrophic respiration rates from bare ground subarea
c
       call  hetresg  (litrmass, soilcmas,            
     1                      il1,      il2,     tbarg,   
     2                   thliqg,     sand,      clay,   zbotw,   
     3                       fg,        0,
     4                     isand,
     5                   ltrsbrg,  scrsbrg)
c
c     find heterotrophic respiration rates from snow over ground 
c     subarea
c
       call  hetresg  (litrmass, soilcmas,            
     1                      il1,      il2,    tbargs,   
     2                   thliqg,     sand,      clay,   zbotw,   
     3                      fgs,        1,
     4                     isand,
     5                   ltrsbrgs, scrsbrgs)
c
c
c     find vegetation averaged litter and soil c respiration rates
c     using values from canopy over ground and canopy over snow subareas
c
      do 340 j = 1, icc
        do 350 i = il1, il2
          if( (fcanc(i,j)+fcancs(i,j)).gt.zero) then
            ltresveg(i,j)= (fcanc(i,j)*ltrsvgcg(i,j) + 
     &        fcancs(i,j)*ltrsvgcs(i,j)) / ( fcanc(i,j) + fcancs(i,j))     
            scresveg(i,j)= (fcanc(i,j)*scrsvgcg(i,j) + 
     &        fcancs(i,j)*scrsvgcs(i,j)) / ( fcanc(i,j) + fcancs(i,j))
            hetrsveg(i,j) =  ltresveg(i,j) + scresveg(i,j)
          else
            ltresveg(i,j)= 0.0
            scresveg(i,j)= 0.0
            hetrsveg(i,j)= 0.0
          endif
          nepveg(i,j)=nppveg(i,j)-hetrsveg(i,j)
350     continue 
340   continue 
c
c     find litter and soil c respiration rates averaged over the bare 
c     fraction of the grid cell using values from ground and snow over
c     ground sub-areas.
c
      do 355 i = il1, il2
        if( (fg(i)+fgs(i)).gt.zero) then
          ltresveg(i,iccp1)= (fg(i)*ltrsbrg(i) + 
     &      fgs(i)*ltrsbrgs(i)) / ( fg(i) + fgs(i) )     
          scresveg(i,iccp1)= (fg(i)*scrsbrg(i) + 
     &      fgs(i)*scrsbrgs(i)) / ( fg(i) + fgs(i) )     
          hetrsveg(i,iccp1) =  ltresveg(i,iccp1) + scresveg(i,iccp1)
        else
          ltresveg(i,iccp1)= 0.0
          scresveg(i,iccp1)= 0.0
          hetrsveg(i,iccp1)= 0.0
        endif
355   continue
c
c     find grid averaged litter and soil c respiration rates
c
      do 360 j = 1,icc
        do 370 i = il1, il2
          litres(i)=litres(i)+fcancmx(i,j)*ltresveg(i,j)
          socres(i)=socres(i)+fcancmx(i,j)*scresveg(i,j)
370     continue
360   continue
c
      do 380 i = il1, il2
        litres(i)=litres(i)+( (fg(i)+fgs(i))*ltresveg(i,iccp1))
        socres(i)=socres(i)+( (fg(i)+fgs(i))*scresveg(i,iccp1))
        hetrores(i)= litres(i)+socres(i)
        nep(i)=npp(i)-hetrores(i)
380   continue
c
c     ---------------------------------------------------------------
c
c     update the litter and soil c pools based on litter and soil c
c     respiration rates found above. also transfer humidified litter 
c     to the soil c pool.
c
      do 420 j = 1, iccp1
        do 430 i = il1, il2
c         convert u mol co2/m2.sec -> kg c/m2 respired over the model
c         time step
          ltrestep(i,j)=ltresveg(i,j)*(1.0/963.62)*deltat
          screstep(i,j)=scresveg(i,j)*(1.0/963.62)*deltat
c
c         update litter and soil c pools
          if (j .ne. iccp1) then
           litrmass(i,j)=litrmass(i,j)-(ltrestep(i,j)*
     &                   (1.0+humicfac(sort(j))))
           hutrstep(i,j)=(humicfac(sort(j))* ltrestep(i,j))
          else
           litrmass(i,j)=litrmass(i,j)-(ltrestep(i,j)*(1.0+0.45))
           hutrstep(i,j)=(0.45 * ltrestep(i,j))
          endif
c
          humtrsvg(i,j)=hutrstep(i,j)*(963.62/deltat) ! u-mol co2/m2.sec
          soilcmas(i,j)=soilcmas(i,j) + 
     &          real(spinfast) * (hutrstep(i,j) -  screstep(i,j)) 
c
          if(litrmass(i,j).lt.zero) litrmass(i,j)=0.0
          if(soilcmas(i,j).lt.zero) soilcmas(i,j)=0.0
430     continue
420   continue
c
c     estimate soil respiration. this is sum of heterotrophic respiration
c     and root maintenance respiration.
c
      do 440 j = 1, icc
        do 450 i = il1, il2
          soilrsvg(i,j)=ltresveg(i,j)+scresveg(i,j)+rmrveg(i,j)
450     continue
440   continue
c
c     but over the bare fraction there is no live root.
c
      do 460 i = il1, il2
        soilrsvg(i,iccp1)=ltresveg(i,iccp1)+scresveg(i,iccp1)
460   continue
c
c     find grid averaged humification and soil respiration rates
c
      do 470 j = 1,icc
        do 480 i = il1, il2
          soilresp(i)=soilresp(i)+fcancmx(i,j)*soilrsvg(i,j)
          humiftrs(i)=humiftrs(i)+fcancmx(i,j)*humtrsvg(i,j)
480     continue
470   continue
c
      do 490 i = il1, il2
        soilresp(i)=soilresp(i)+( (fg(i)+fgs(i))*soilrsvg(i,iccp1))
        humiftrs(i)=humiftrs(i)+( (fg(i)+fgs(i))*humtrsvg(i,iccp1))
490   continue
c
c     heterotrophic respiration part ends
c
c     ------------------------------------------------------------------
c
c     estimate allocation fractions for leaf, stem, and root components.
c
           call allocate (lfstatus,   thliqc,    ailcg,     ailcb,
     1                     il1, il2,     sand,     clay,  
     3                    rmatctem,   gleafmas, stemmass, rootmass,      
     4                       sort,    nol2pfts,  fcancmx,
     5                     afrleaf,  afrstem,  afrroot,    wiltsm,
     6                     fieldsm, wtstatus, ltstatus)
c  
c     estimate fraction of npp that is to be used for horizontal
c     expansion (lambda) during the next day. the remaining fraction (1
c     - lambda) is what we will use for vertical expansion. 
c

      if (compete) then
       do 500 j = 1, icc
        if(j.ne.6.and.j.ne.7) then   ! not for crops
         do 501 i = il1, il2
c
           n = sort(j)
           if(ailcg(i,j).le.laimin(n))then
              lambda(i,j)=0.0
           else if(ailcg(i,j).ge.laimax(n))then
              lambda(i,j)=lambdamax
           else
              lambda(i,j)=((ailcg(i,j)-laimin(n))*lambdamax)/
     &                    (laimax(n)-laimin(n))
           endif
           lambda(i,j)=max(0.0, min(lambdamax, lambda(i,j)))
c
c          if tree and leaves still coming out, or if npp is negative, then
c          do not expand
           if((j.le.5.and.lfstatus(i,j).eq.1).or.nppveg(i,j).lt.0.0
     &     .or..not.pftexist(i,j))then
             lambda(i,j)=0.0
           endif
c
501      continue
        endif
500    continue
      endif !compete       
c
c    ------------------------------------------------------------------
c
c     maintenance respiration also reduces leaf, stem, and root biomass.
c     when npp for a given pft is positive then this is taken care by
c     allocating +ve npp amongst the leaves, stem, and root component.
c     when npp for a given pft is negative then maintenance respiration
c     loss is explicitly deducted from each component.
c
      do 600 j = 1, icc
        do 610 i = il1, il2
c
c         convert npp and maintenance respiration from different components
c         from units of u mol co2/m2.sec -> kg c/m2 sequestered or respired
c         over the model time step          
          gppvgstp(i,j)=gppveg(i,j)*(1.0/963.62)*deltat + add2allo(i,j)
          nppvgstp(i,j)=nppveg(i,j)*(1.0/963.62)*deltat*(1.-lambda(i,j))
     &                  + add2allo(i,j)
c
c         amount of c related to horizontal expansion
c
          expbalvg(i,j)=-1.0*nppveg(i,j)*deltat*lambda(i,j)
     &                  + add2allo(i,j)*(963.62/1.0)
c
          rmlvgstp(i,j)=rmlveg(i,j)*(1.0/963.62)*deltat
          rmsvgstp(i,j)=rmsveg(i,j)*(1.0/963.62)*deltat
          rmrvgstp(i,j)=rmrveg(i,j)*(1.0/963.62)*deltat
c
          if(lfstatus(i,j).ne.4)then
            if(nppvgstp(i,j).gt.0.0) then
              ntchlveg(i,j)=afrleaf(i,j)*nppvgstp(i,j)
              ntchsveg(i,j)=afrstem(i,j)*nppvgstp(i,j)
              ntchrveg(i,j)=afrroot(i,j)*nppvgstp(i,j)
            else
              ntchlveg(i,j)=-rmlvgstp(i,j)+afrleaf(i,j)*gppvgstp(i,j)
              ntchsveg(i,j)=-rmsvgstp(i,j)+afrstem(i,j)*gppvgstp(i,j)
              ntchrveg(i,j)=-rmrvgstp(i,j)+afrroot(i,j)*gppvgstp(i,j)
            endif
          else  ! i.e. if lfstatus.eq.4
c           and since we do not have any real leaves on then we do not take
c           into account co2 uptake by imaginary leaves in carbon budget.
c           rmlvgstp(i,j) should be zero because we set maintenance
c           respiration from storage/imaginary leaves equal to zero. 
c           in loop 180 
c
            ntchlveg(i,j)=-rmlvgstp(i,j) 
            ntchsveg(i,j)=-rmsvgstp(i,j)
            ntchrveg(i,j)=-rmrvgstp(i,j)
c
c           since no real leaves are on, make allocation fractions equal to
c           zero.
c
            afrleaf(i,j)=0.0
            afrstem(i,j)=0.0
            afrroot(i,j)=0.0
          endif
c
          gleafmas(i,j)=gleafmas(i,j)+ntchlveg(i,j)
          stemmass(i,j)=stemmass(i,j)+ntchsveg(i,j)
          rootmass(i,j)=rootmass(i,j)+ntchrveg(i,j)
c
c
          if(gleafmas(i,j).lt.0.0)then
            write(6,1900)'gleafmas lt zero at i=',i,' for pft=',j,''   
            write(6,1901)'gleafmas = ',gleafmas(i,j)
            write(6,1901)'ntchlveg = ',ntchlveg(i,j)
            write(6,1902)'lfstatus = ',lfstatus(i,j)
            write(6,1901)'ailcg    = ',ailcg(i,j)
            write(6,1901)'slai     = ',slai(i,j)
1900        format(a23,i4,a10,i2,a1)
1902        format(a11,i4)
            call xit ('ctem',-6)
          endif
c
          if(stemmass(i,j).lt.0.0)then
            write(6,1900)'stemmass lt zero at i=(',i,') for pft=',j,')'   
            write(6,1901)'stemmass = ',stemmass(i,j)
            write(6,1901)'ntchsveg = ',ntchsveg(i,j)
            write(6,1902)'lfstatus = ',lfstatus(i,j)
            write(6,1901)'rmsvgstp = ',rmsvgstp(i,j)
            write(6,1901)'afrstem  = ',afrstem(i,j)
            write(6,1901)'gppvgstp = ',gppvgstp(i,j)
            write(6,1901)'rmscsveg = ',rmscsveg(i,j)
            write(6,1901)'rmscgveg = ',rmscgveg(i,j)
1901        format(a11,f12.8)
            call xit ('ctem',-7)
          endif
c
          if(rootmass(i,j).lt.0.0)then
            write(6,1900)'rootmass lt zero at i=(',i,') for pft=',j,')'   
            write(6,1901)'rootmass = ',rootmass(i,j)
            call xit ('ctem',-8)
          endif
c
c         convert net change in leaf, stem, and root biomass into 
c         u-mol co2/m2.sec for use in balcar subroutine
c          
          ntchlveg(i,j)=ntchlveg(i,j)*(963.62/deltat)         
          ntchsveg(i,j)=ntchsveg(i,j)*(963.62/deltat)         
          ntchrveg(i,j)=ntchrveg(i,j)*(963.62/deltat)         
c
c         to avoid over/underflow problems set gleafmas, stemmass, and
c         rootmass to zero if they get too small
c
          if(bleafmas(i,j).lt.zero) bleafmas(i,j)=0.0
          if(gleafmas(i,j).lt.zero) gleafmas(i,j)=0.0
          if(stemmass(i,j).lt.zero) stemmass(i,j)=0.0
          if(rootmass(i,j).lt.zero) rootmass(i,j)=0.0
c
610     continue
600   continue
c  
c     calculate grid averaged value of c related to spatial expansion
c
      do 620 j = 1,icc
        do 621 i = il1, il2
c      if (compete .or. lnduseon) then
          expnbaln(i)=expnbaln(i)+fcancmx(i,j)*expbalvg(i,j)
c      else
c          expnbaln(i)=0.0
c      endif
621     continue
620   continue
c
c
c    ------------------------------------------------------------------
c
c     phenology part starts
c
c     the phenology subroutine determines leaf status for each pft and 
c     calculates leaf litter. the phenology subroutine uses soil 
c     temperature (tbar) and root temperature. however, since ctem
c     doesn't make the distinction between canopy over ground, and
c     canopy over snow sub-areas for phenology purposes (for  example,
c     leaf onset is not assumed to occur at different times over these
c     sub-areas) we use average soil and root temperature in the phenology
c     subroutine.
c
c     calculate average soil temperature and root temperature using
c     values for canopy over ground and canopy over snow sub-areas, for
c     each vegetation type.
c
      do 650 j = 1, icc
        do 660 i = il1, il2
          if( (fcanc(i,j)+fcancs(i,j)).gt.zero) then
            roottemp(i,j)= (fcanc(i,j)*rttempcg(i,j) + 
     &        fcancs(i,j)*rttempcs(i,j)) / ( fcanc(i,j) + fcancs(i,j))     
          else
            roottemp(i,j)= rttempcg(i,j)
          endif
660     continue
650   continue
c
      do 680 j = 1, ignd
        do 690 i = il1, il2
          if( (fc(i)+fcs(i)).gt.zero) then
            tbarccs(i,j)= (fc(i)*tbarc(i,j) + 
     &        fcs(i)*tbarcs(i,j)) / ( fc(i) + fcs(i))     
          else
            tbarccs(i,j)= tbar(i,j)
          endif
690     continue
680   continue
c
c    -------------------------------------------------------------------
c
c     call the phenology subroutine, which determines the leaf growth
c     status, calculates leaf litter, and converts green grass into
c     brown.
c
            call phenolgy(gleafmas, bleafmas, 
     1                         il1,      il2,  tbarccs,
     2                      thliqc,   wiltsm,  fieldsm,       ta,
     3                    pheanveg,     iday,     radj, roottemp,
     4                    rmatctem, stemmass, rootmass,     sort,
     5                    nol2pfts,  fcancmx,
     6                    flhrloss, leaflitr, lfstatus,  pandays,
     7                    colddays)
c
c    -------------------------------------------------------------------
c

c     while leaf litter is calculated in the phenology subroutine, stem
c     and root turnover is calculated in the turnover subroutine.
c
            call turnover (stemmass, rootmass,  lfstatus,    ailcg,
     1                          il1,      il2,
     2                         sort, nol2pfts,  fcancmx,
     3                     stmhrlos, rothrlos,
     4                     stemlitr, rootlitr)
c
c    -------------------------------------------------------------------
c
c     update green leaf biomass for trees and crops and brown leaf biomass 
c     for grasses
c
      k1=0
      do 700 j = 1, ican 
       if(j.eq.1) then
         k1 = k1 + 1
       else
         k1 = k1 + nol2pfts(j-1)
       endif
       k2 = k1 + nol2pfts(j) - 1
       do 705 m = k1, k2
        do 710 i = il1, il2

          if(j.le.3)then    ! trees and crops
            gleafmas(i,m)=gleafmas(i,m)-leaflitr(i,m)
            if( gleafmas(i,m).lt.0.0) then
              leaflitr(i,m)=leaflitr(i,m)+gleafmas(i,m)
              gleafmas(i,m)=0.0
            endif
          else              ! grasses
            bleafmas(i,m)=bleafmas(i,m)-leaflitr(i,m)
            if( bleafmas(i,m).lt.0.0) then
              leaflitr(i,m)=leaflitr(i,m)+bleafmas(i,m)
              bleafmas(i,m)=0.0
            endif
          endif

710     continue
705    continue
700   continue
c
c     update stem and root biomass for litter deductions
c
      do 780 j = 1, icc
        do 790 i = il1, il2
          stemmass(i,j)=stemmass(i,j)-stemlitr(i,j)
          if( stemmass(i,j).lt.0.0) then
            stemlitr(i,j)=stemlitr(i,j)+stemmass(i,j)
            stemmass(i,j)=0.0
          endif
c
          rootmass(i,j)=rootmass(i,j)-rootlitr(i,j)
          if( rootmass(i,j).lt.0.0) then
            rootlitr(i,j)=rootlitr(i,j)+rootmass(i,j)
            rootmass(i,j)=0.0
          endif
790     continue
780   continue
c
c     update litter pool with leaf litter calculated in the phenology 
c     subroutine and stem and root litter calculated in the turnover
c     subroutine.
c
      do 800 j = 1, icc
        do 810 i = il1, il2
          litrmass(i,j)=litrmass(i,j)+leaflitr(i,j)+stemlitr(i,j)+
     &                  rootlitr(i,j)
810     continue
800   continue
c
c
c    ------------------------------------------------------------------
c
c     call the mortaliy subroutine which calculates mortality due to 
c     reduced growth and aging. exogenous mortality due to fire and other 
c     disturbances and the subsequent litter that is generated is 
c     calculated in the disturb subroutine.
c    
c     set do_mortality=.false. to switch off mortality due to age and 
c     reduced growth. mortality is linked to the competition parameterization 
c     and generates bare fraction.
c
      if (compete) then
        do_mortality=.true.
      else
        do_mortality=.false.
      end if  
c
      call       mortalty (stemmass, rootmass,        ailcg, gleafmas,
     1                     bleafmas,      il1, 
     2                          il2,     iday, do_mortality,     sort,
     3                      fcancmx, lystmmas,     lyrotmas, tymaxlai,
     4                     grwtheff, stemltrm,     rootltrm, glealtrm,
     5                     geremort, intrmort)
c
c    ------------------------------------------------------------------
c
c     update leaf, stem, and root biomass pools to take into loss
c     due to mortality, and put the litter into the litter pool. the 
c     mortality for green grasses doesn't generate litter, instead
c     they turn brown.
c
      k1=0
      do 830 j = 1, ican 
       if(j.eq.1) then
         k1 = k1 + 1
       else
         k1 = k1 + nol2pfts(j-1)
       endif
       k2 = k1 + nol2pfts(j) - 1
       do 835 m = k1, k2
        do 840 i = il1, il2
c
          stemmass(i,m)=stemmass(i,m)-stemltrm(i,m)
          rootmass(i,m)=rootmass(i,m)-rootltrm(i,m)
          litrmass(i,m)=litrmass(i,m)+stemltrm(i,m)+rootltrm(i,m)  
c
          if(j.eq.4)then    ! grasses
            gleafmas(i,m)=gleafmas(i,m)-glealtrm(i,m)
            bleafmas(i,m)=bleafmas(i,m)+glealtrm(i,m)
            glealtrm(i,m)=0.0
          else              ! trees and crops
            gleafmas(i,m)=gleafmas(i,m)-glealtrm(i,m)
          endif
          litrmass(i,m)=litrmass(i,m)+glealtrm(i,m)
c
840     continue
835    continue 
830   continue 
c
c    ------------------------------------------------------------------
c
c     call the disturbance subroutine which calculates mortality due to 
c     fire and other disturbances. the primary output from from 
c     disturbance subroutine is litter generated, c emissions due to 
c     fire and area burned, which may be used to estimate change in 
c     fractional coverages.
c
c     disturbance is spatial and requires area of gcm grid cell and areas
c     of different pfts present in a given grid cell. however, when ctem is
c     operated at a point scale then it is assumed that the spatial scale 
c     is 1 hectare = 10,000 m2. the disturbance subroutine may be stopped 
c     from simulating any fire by specifying fire extingushing probability
c     equal to 1.
c
c
            call disturb (stemmass, rootmass, gleafmas, bleafmas,
     1                      thliqc,   wiltsm,  fieldsm,    uwind,
     2                       vwind,  lightng,  fcancmx, litrmass,
     3                    prbfrhuc, rmatctem, extnprob, pftareab,
     4                         il1,      il2,     sort, nol2pfts,
     6                    grclarea,   thicec,   popdin, lucemcom,
     7                      dofire,
c    in above, out below 
     8                    stemltdt, rootltdt, glfltrdt, blfltrdt,
     9                    pftareaa, glcaemls, rtcaemls, stcaemls,
     a                    blcaemls, ltrcemls, burnfrac, probfire,
     b                    emit_co2, emit_co,  emit_ch4, emit_nmhc,
     c                    emit_h2,  emit_nox, emit_n2o, emit_pm25,
     d                    emit_tpm, emit_tc,  emit_oc,  emit_bc,
     e                    burnveg )
c
c    ------------------------------------------------------------------
c
c     update litter pool and reduce leaf, stem, and root biomass to
c     take into account losses from disturbance. calculate nbp (net biome
c     production) for each pft by taking into account c emission losses.
c
      do 1000 j = 1, icc
        do 1010 i = il1, il2
          gleafmas(i,j)=gleafmas(i,j) - glfltrdt(i,j) - glcaemls(i,j)
          bleafmas(i,j)=bleafmas(i,j) - blfltrdt(i,j) - blcaemls(i,j)
          stemmass(i,j)=stemmass(i,j) - stemltdt(i,j) - stcaemls(i,j)
          rootmass(i,j)=rootmass(i,j) - rootltdt(i,j) - rtcaemls(i,j)
          litrmass(i,j)=litrmass(i,j) + glfltrdt(i,j) + blfltrdt(i,j) +  
     &                  stemltdt(i,j) + rootltdt(i,j) - ltrcemls(i,j)
          dscemlv1(i,j)=glcaemls(i,j) + blcaemls(i,j) + stcaemls(i,j) +
     &                  rtcaemls(i,j) 
          dscemlv2(i,j)=glcaemls(i,j) + blcaemls(i,j) + stcaemls(i,j) +
     &                  rtcaemls(i,j) + ltrcemls(i,j)
c         convert kg c/m2 emitted in one day into u mol co2/m2.sec before
c         subtracting emission losses from nep
          nbpveg(i,j)  =nepveg(i,j)   - dscemlv2(i,j)*(963.62/deltat)
1010    continue
1000  continue
c
c     calculate grid. averaged rate of carbon emissions due to fire in
c     u-mol co2/m2.sec. convert all emission losses from kg c/m2
c     emitted in 1 day to u-mol co2/m2.sec. calculate grid averaged
c     carbon emission losses from litter.
c
      do 1030 j = 1,icc
        do 1040 i = il1, il2
          dstcemls1(i)=dstcemls1(i) +
     &     fcancmx(i,j)*dscemlv1(i,j)*(963.62/deltat)
          dstcemls2(i)=dstcemls2(i) +
     &     fcancmx(i,j)*dscemlv2(i,j)*(963.62/deltat)
          galtcels(i)=galtcels(i) +
     &     fcancmx(i,j)*ltrcemls(i,j)*(963.62/deltat)
          glcaemls(i,j)=glcaemls(i,j)*(963.62/deltat)
          blcaemls(i,j)=blcaemls(i,j)*(963.62/deltat)
          stcaemls(i,j)=stcaemls(i,j)*(963.62/deltat)
          rtcaemls(i,j)=rtcaemls(i,j)*(963.62/deltat)
          ltrcemls(i,j)=ltrcemls(i,j)*(963.62/deltat)
1040    continue
1030  continue
c

      do 1041 i = il1, il2
        nbp(i)=nep(i)-dstcemls2(i)
        dstcemls3(i)=dstcemls2(i)-dstcemls1(i)
1041  continue
c
c     calculate total litter fall from each component (leaves, stem, and
c     root) from all causes (normal turnover, drought and cold stress for
c     leaves, mortality, and disturbance) for use in balcar subroutine
c
      do 1050 j = 1,icc
        do 1060 i = il1, il2
c     
c        units here are kg c/m2.day
         tltrleaf(i,j)=leaflitr(i,j)+glealtrm(i,j)+glfltrdt(i,j)+
     &                 blfltrdt(i,j)
         tltrstem(i,j)=stemlitr(i,j)+stemltrm(i,j)+stemltdt(i,j)
         tltrroot(i,j)=rootlitr(i,j)+rootltrm(i,j)+rootltdt(i,j)
c          
c        convert units to u-mol co2/m2.sec
         leaflitr(i,j)=leaflitr(i,j)*(963.62/deltat)
         tltrleaf(i,j)=tltrleaf(i,j)*(963.62/deltat)
         tltrstem(i,j)=tltrstem(i,j)*(963.62/deltat)
         tltrroot(i,j)=tltrroot(i,j)*(963.62/deltat)
1060    continue
1050  continue
c
c     calculate grid-average vegetation biomass, litter mass, and soil
c     carbon mass, and litter fall rate
c
      do 1100 j = 1, icc
        do 1110 i = il1, il2
          vgbiomas(i)=vgbiomas(i)+fcancmx(i,j)*(gleafmas(i,j)+
     &     bleafmas(i,j)+stemmass(i,j)+rootmass(i,j))
          litrfall(i)=litrfall(i)+fcancmx(i,j)*(tltrleaf(i,j)+
     &     tltrstem(i,j)+tltrroot(i,j))
          gavgltms(i)=gavgltms(i)+fcancmx(i,j)*litrmass(i,j)
          gavgscms(i)=gavgscms(i)+fcancmx(i,j)*soilcmas(i,j)
c          gavglai (i)=gavglai (i)+fcancmx(i,j)*ailcg(i,j)
          vgbiomas_veg(i,j)=gleafmas(i,j)+
     &     bleafmas(i,j)+stemmass(i,j)+rootmass(i,j) !vegetation biomass for each pft
1110    continue
1100  continue
c
      do 1020 i = il1, il2
        gavgltms(i)=gavgltms(i)+( (fg(i)+fgs(i))*litrmass(i,iccp1))
        gavgscms(i)=gavgscms(i)+( (fg(i)+fgs(i))*soilcmas(i,iccp1))
1020  continue
c
c     -----------------------------------------------------------------
c
c     at this stage we have all required fluxes in u-mol co2/m2.sec and
c     initial (loop 140 and 145) and updated sizes of all pools 
c     (in kg c/m2). now we call the balcar subroutine and make sure that
c     c in leaves, stem, root, litter and soil c pool balances within a
c     certain tolerance.
c
      if(spinfast.eq.1)then
             call  balcar(gleafmas, stemmass, rootmass,  bleafmas,
     1                    litrmass, soilcmas, ntchlveg,  ntchsveg,
     2                    ntchrveg, tltrleaf, tltrstem,  tltrroot,
     3                    glcaemls, blcaemls, stcaemls,  rtcaemls,
     4                    ltrcemls, ltresveg, scresveg,  humtrsvg,
     5                    pglfmass, pblfmass, pstemass,  protmass,
     6                    plitmass, psocmass, vgbiomas,
     7                    pvgbioms, gavgltms, pgavltms,  gavgscms,
     8                    pgavscms, galtcels, expnbaln,
     9                         npp,  autores, hetrores,       gpp,
     a                         nep,   litres,   socres, dstcemls1,
     b                         nbp, litrfall, humiftrs,
     c                         il1,       il2)
      endif
c
c     -----------------------------------------------------------------
c
c     the luc subroutine also estimates carbon emissions due to
c     combustion associated with luc, and this flux is spread out over 
c     the whole year and is therefore subtracted to get nbp of each pft
c     as well as the grid averaged value of nbp.
c
      if(lnduseon)then
          do 1150 j = 1, icc
            do 1151 i = il1, il2
c             luc related combustion flux is assumed to be spread
c             uniformly over the grid cell and thus reduces nbp of each
c             pft
              nbpveg(i,j)=nbpveg(i,j)-lucemcom(i) 
1151        continue
1150      continue
c
          do 1160 i = il1, il2
            nbp(i)=nbp(i)-lucemcom(i)
1160      continue
      endif        
c
c     -----------------------------------------------------------------
c
c     finally find vegetation structural attributes which can be passed
c     to the land surface scheme using leaf, stem, and root biomass. 
c
              call bio2str( gleafmas, bleafmas, stemmass, rootmass,
     1                            il1,      il2, fcancmx,    zbotw,
     3                          delzw, nol2pfts,  soildpth,
     4                          ailcg,    ailcb,     ailc,    zolnc,
     5                          rmatc, rmatctem,     slai,  bmasveg,
     6                       cmasvegc,  veghght, rootdpth,   alvisc,
     8                         alnirc,  paicgat,  slaicgat  )

c
c    calculation of gavglai is moved from loop 1100 to here 
c    since ailcg is updated by bio2str
c
      do j = 1, icc
        do i = il1, il2
          gavglai (i)=gavglai (i)+fcancmx(i,j)*ailcg(i,j)
        enddo
      enddo
c
c
c     -----------------------------------------------------------------
c
      return
      end

