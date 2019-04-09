!>\file
!!Principle driver for CTEM
!!@author V. Arora, J. Melton, Y. Peng, R. Shrestha
!!
module ctemDriver

implicit none

! subroutines contained in this module:

public  :: ctem   
public  :: calcNBP 
!public  :: calcNPP

contains

  subroutine ctem( fcancmx,    fsnow,     sand,      clay,  & ! In
                 ilg,   il1,      il2,     iday,      radj,  & ! In
                         ta,     delzw, ancgveg,   rmlcgveg, & ! In
                      zbotw,  & ! In
                      uwind,    vwind,  lightng,      tbar,  & ! In
                  pfcancmx, nfcancmx,             & ! In
                   soildpth, spinfast,   todfrac,& ! In
                     netrad,   precip,    psisat,            & ! In
                   grclarea,   popdin,    isand,             & ! In
                   wetfrac, slopefrac,       bi,             & ! In
                     thpor,   currlat,  ch4conc,  & ! In
                      THFC,      THLW,     thliq,  thice,    & ! In
                 ipeatland,    anmoss,   rmlmoss,  gppmoss,  & ! In
                  wtable,   maxAnnualActLyr,         & ! In
                  PFTCompetition,  dofire,  lnduseon,  inibioclim,  & ! In
                     leapnow,                                   & ! In
                   stemmass, rootmass, litrmass,  gleafmas,& ! In/ Out
                   bleafmas, soilcmas,    ailcg,      ailc,& ! In/ Out
                      zolnc, rmatctem,    rmatc,     ailcb,& ! In/ Out
                   flhrloss,  pandays, lfstatus,  grwtheff,& ! In/ Out
                   lystmmas, lyrotmas, tymaxlai,  vgbiomas,& ! In/ Out
                   gavgltms, gavgscms, stmhrlos,      slai, & ! In/ Out
                    bmasveg, cmasvegc, colddays,  rothrlos,& ! In/ Out
                     fcanmx,   alvisc,   alnirc,   gavglai,& ! In/ Out
                   Cmossmas, litrmsmoss,     peatdep,      &! In/ Out
                   geremort,   intrmort,   pstemmass,    pgleafmass, &! In/ Out
                      tcurm,   srpcuryr,    dftcuryr,      lambda,   &! In/ Out
                     tmonth, anpcpcur,  anpecur,   gdd5cur,&! In/ Out
                   surmncur, defmncur, srplscur,  defctcur,&! In/ Out
                    aridity, srplsmon, defctmon,  anndefct,&! In/ Out
                   annsrpls,  annpcp,dry_season_length,&! In/ Out
                   pftexist,   twarmm,       tcoldm,         gdd5,   &! In/ Out
                       npp,       nep, hetrores,   autores,& ! Out (Primary)
                  soilresp,        rm,       rg,       nbp,& ! Out (Primary)
                    litres,    socres,      gpp, dstcemls1,& ! Out (Primary)
                  litrfall,  humiftrs,  veghght,  rootdpth,& ! Out (Primary)
                       rml,       rms,      rmr,  tltrleaf,& ! Out (Primary)
                  tltrstem,  tltrroot, leaflitr,  roottemp,& ! Out (Primary)
                  burnfrac,                 lucemcom,    lucltrin,   & ! Out (Primary)
               lucsocin,   dstcemls3,                             & ! Out (Primary)
                ch4WetSpec,  ch4WetDyn,      wetfdyn,   ch4soills,   & ! Out (Primary)
                                paicgat,    slaicgat,                & ! Out (Primary)
                 emit_co2,   emit_ch4, reprocost, blfltrdt, glfltrdt, &  ! Out (Primary)
                 glcaemls, blcaemls, rtcaemls, stcaemls,  ltrcemls, &  ! Out (Primary)
                 ntchlveg, ntchsveg, ntchrveg,  mortLeafGtoB,       &  ! Out (Primary)
                 phenLeafGtoB,  turbLitter, turbSoilC,             &  ! Out (Primary)
                 gLeafLandCompChg, bLeafLandCompChg, stemLandCompChg, &! Out (Primary)
                 rootLandCompChg, litterLandCompChg, soilCLandCompChg, &! Out (Primary)
                  emit_co,   emit_nmhc,  smfunc_veg,                & ! Out (Secondary)
                   emit_h2,  emit_nox, emit_n2o, emit_pm25,& ! Out (Secondary)
                   emit_tpm, emit_tc,  emit_oc,    emit_bc,& ! Out (Secondary)
                 bterm_veg,      lterm,    mterm_veg,  burnvegf,     & ! Out (Secondary)
               litrfallveg,    humtrsvg,    ltstatus,      nppveg,   & ! Out (Secondary)
                   afrleaf,     afrstem,     afrroot,    wtstatus,   & ! Out (Secondary)
                     rmlveg,  rmsveg,   rmrveg,    rgveg,& ! Out (Secondary)
               vgbiomas_veg,  gppveg,   nepveg,   nbpveg,& ! Out (Secondary)
                 hetrsveg,autoresveg, ltresveg, scresveg,& ! Out (Secondary)
                   nppmoss,     armoss,                  & ! Out (Secondary)
                        cc,         mm) ! Out (Secondary)

!
!             Canadian Terrestrial Ecosystem Model (CTEM)
!             Main Ctem Subroutine Compatible With CLASS

!     29 Mar 2019 - Rewrote for clarity, removed sub-area TBAR, THLIQ, THICE, ANVEG, rmlveg,
!     J. Melton     was needlessly confusing and didn't have impact. Added subroutines to make 
!                   this a driver, not a main calculation subroutine.
!     28  Nov 2018  - Clean up argument list before implementation in AGCM
!     V. Arora        
!
!     14  Mar 2016  - Remove grid cell area calculation to the driver. This will
!     J. Melton       harmonize this subroutine with the coupled code.
!
!      8  Feb 2016  - Adapted subroutine for multilayer soilc and litter (fast decaying)
!     J. Melton       carbon pools
!
!     3   Feb 2016  - Bring in onetile_perPFT switch so now in mosaic mode the
!     J. Melton       tiles in a grid cell are treated independently by competition
!                     and LUC.
!     19  Jan 2016  - Implemented new LUC litter and soil C pools
!     J. Melton
!
!     3   Jul 2014  - Bring in wetland and wetland methane code
!     R. Shrestha
!
!     12  Jun 2014  - Bring in a constant reproductive cost, remove expnbaln,
!     J. Melton       add a smoothing function for lambda calculation for competition,
!                     made it so NEP and NBP work with competition on.
!
!     17  Jan 2014  - Moved parameters to global file (classic_params.f90)
!     J. Melton
!
!     Dec 6   2012   Make it so competition and luc can function in both
!     J. Melton      composite and mosaic modes.
!
!     sep 25  2012   Add competition_map and competition_unmap
!     Y. Peng
!
!     Sep 12  2012
!     J. Melton     Add in bottom limit to bleafmas to ensure it does no
!                   slowly decay to infintesimaly small number.
!
!     Aug 23  2012
!     J. Melton     Pass in isand to ensure soil levels are properly
!                   labelled as bedrock if assigned so in classb
!
!     Jan 10  2012
!     Yiran         Re-test bioclim and existence for competition
!
!     19  Sep. 2001 - This is the main terrestrial carbon model subrouti
!     V. Arora
!                     all primary ctem subroutines are called from here,
!                     except phtsyn which is called from tsolvc
!
!     08 June  2001 - Add calls to three new subroutines (bioclim,
!     V. Arora        existence, and competition) to add dynamic
!                     competition between pfts. changes are also made
!                     to land use change (luc) subroutine and the manner
!                     in which distrubance is handles. fire now creates
!                     bare ground which is subsequently available for
!                     colonization. other changes are also made to keep
!                     everything consistent with changing vegetation
!                     fractions.
!    -----------------------------------------------------------------

  use classic_params,        only : kk, pi, zero, icp1, &
                                kn,iccp1, ican, nlat,&
                                ignd, icc, nmos, l2max, grescoef,&
                                humicfac,laimin,laimax,lambdamax,&
                                crop,repro_fraction,grescoefmoss,&
                                rmortmoss,humicfacmoss,GRAV,RHOW,RHOICE,&
                               classpfts,ctempfts,iccp2,humicfac_bg,&
                               nol2pfts,deltat
     

  use landuse_change,     only : luc
  use competition_scheme, only : bioclim, existence, competition, expansion 
  use disturbance_scheme, only : disturb
  use heterotrophic_respiration, only : hetresg, hetresv
  use peatlands_mod, only : hetres_peat,peatDayEnd,peatDepth
  use ctemUtilities, only : genSortIndex
  use autotrophicRespiration, only : mainres 
  use balanceCarbon, only : balcar, prepBalanceC
  use mortality, only : mortalty, updatePoolsMortality
  use turnover, only : turnoverStemRoot, updatePoolsTurnover
  use soilC_processes, only : turbation
  use applyAllometry, only : allometry

  implicit none

  ! inputs

  logical, intent(in) :: lnduseon                         !<logical switch to run the land use change subroutine or not.
  logical, intent(in) :: PFTCompetition                   !<logical boolean telling if competition between pfts is on or not
  logical, intent(in) :: dofire                           !<boolean, if true allow fire, if false no fire.
  logical, intent(in) :: leapnow                          !< true if this year is a leap year. Only used if the switch 'leap' is true.
  integer, intent(in) :: iday                             !<day of year
  integer, intent(in) ::  spinfast                        !<spinup factor for soil carbon whose default value is 1. as this factor increases the
                                                          !<soil c pool will come into equilibrium faster. reasonable value for spinfast is
                                                          !<between 5 and 10. when spinfast.ne.1 then the balcar subroutine is not run.
  integer, intent(in) :: ilg                              !<ilg=no. of grid cells in latitude circle
  integer, intent(in) :: il1                              !<il1=1
  integer, intent(in) :: il2                              !<il2=ilg (no. of grid cells in latitude circle)
  integer, dimension(ilg,ignd), intent(in) :: isand       !<
  integer, dimension(ilg), intent(in) :: ipeatland        !< Peatland flag: 0 = not a peatland, 1= bog, 2 = fen
  real, dimension(ilg), intent(in) :: fsnow               !< fraction of snow simulated by class
  real, dimension(ilg,ignd), intent(in) :: sand           !< percentage sand
  real, dimension(ilg,ignd), intent(in) :: clay           !< percentage clay
  real, dimension(ilg), intent(in) :: radj                !< latitude in radians
  real, dimension(ilg,ignd), intent(in) ::  tbar          !<soil temperature, k
  real, dimension(ilg,ignd), intent(in) :: psisat         !< Saturated soil matric potential (m)
  real, dimension(ilg,ignd), intent(in) :: bi             !< Brooks and Corey/Clapp and Hornberger b term
  real, dimension(ilg,ignd), intent(in) :: thpor          !< Soil total porosity \f$(cm^3 cm^{-3})\f$ - daily average
  real, dimension(ilg), intent(in) :: ta                  !< air temp, K
  real, dimension(ilg,ignd), intent(in) :: delzw          !< thicknesses of the soil layers
  real, dimension(ilg,ignd), intent(in) :: zbotw          !< bottom of soil layers
  real, dimension(ilg), intent(in) :: soildpth            !<soil depth (m)
  real, dimension(ilg,ignd), intent(in) :: thliq          !< liquid mois. content of soil layers
  real, dimension(ilg,ignd), intent(in) :: thice          !< Frozen soil moisture content
  real, dimension(ilg), intent(in) ::  grclarea           !< area of the grid cell, \f$km^2\f$
  real, dimension(ilg), intent(in) ::  currlat            !< centre latitude of grid cells in degrees
  real, dimension(ilg), intent(in) :: uwind               !< u wind speed, m/s
  real, dimension(ilg), intent(in) :: vwind               !< v wind speed, m/s
  real, dimension(ilg), intent(in) ::  precip             !<daily precipitation (mm/day)
  real, dimension(ilg), intent(in) ::  netrad             !<daily net radiation (w/m2)
  real, dimension(ilg), intent(in) :: lightng             !< total lightning frequency, flashes/km2.year
  real, dimension(ilg,icc), intent(in) :: pfcancmx        !<previous year's fractional coverages of pfts
  real, dimension(ilg,icc), intent(in) :: nfcancmx        !<next year's fractional coverages of pfts
  real, dimension(ilg,icc), intent(in) :: todfrac         !<max. fractional coverage of ctem's 9 pfts by the end of the day, for use by land use subroutine
  real, dimension(ilg), intent(in) :: ch4conc             !< Atmospheric \f$CH_4\f$ concentration at the soil surface (ppmv)
  real, dimension(ilg), intent(in) :: wetfrac             !< Prescribed fraction of wetlands in a grid cell
  real, dimension(ilg,8), intent(in) :: slopefrac         !<
  real, dimension(ilg), intent(in) :: anmoss              !< moss net photoysnthesis -daily averaged C fluxes rates (umol/m2/s)
  real, dimension(ilg), intent(in) :: rmlmoss             !< moss maintainance respiration -daily averaged C fluxes rates (umol/m2/s)
  real, dimension(ilg), intent(in) :: gppmoss             !< moss GPP -daily averaged C fluxes rates (umol/m2/s)
  real, dimension(ilg), intent(in) :: wtable              !< water table (m)
  real, dimension(ilg,ignd), intent(in) :: THFC           !<
  real, dimension(ilg,ignd), intent(in) :: THLW           !<
  real, dimension(ilg), intent(in) :: maxAnnualActLyr     !< Active layer depth maximum over the e-folding period specified by parameter eftime (m).

  !
  !     updates
  !
  logical, intent(inout) :: pftexist(ilg,icc)             !<
  logical, intent(inout) :: inibioclim                    !<switch telling if bioclimatic parameters are being initialized from scratch (false)
                                                          !<or being initialized from some spun up values(true).
  integer, dimension(ilg,icc), intent(inout) :: pandays   !<days with positive net photosynthesis (an) for use in the phenology subroutine
  integer, dimension(ilg,2), intent(inout) :: colddays    !<cold days counter for tracking days below a certain temperature threshold for ndl dcd and crop pfts.
  integer, dimension(ilg,icc), intent(inout) :: lfstatus  !<leaf phenology status
  real, dimension(ilg,icc), intent(inout) :: ancgveg      !< net photosynthetic rate for CTEM's pfts 
  real, dimension(ilg,icc), intent(inout) :: rmlcgveg     !< leaf respiration rate for CTEM's pfts
  real, dimension(ilg,icc), intent(inout) :: fcancmx      !< max. fractional coverage of CTEM's pfts, but this can be
                                                          !< modified by land-use change, and competition between pfts
  real, dimension(ilg,ican,ignd), intent(inout) :: rmatc  !<fraction of roots for each of class' 4 pfts in each soil layer
  real, dimension(ilg), intent(inout) :: surmncur         !<number of months with surplus water for current year
  real, dimension(ilg), intent(inout) :: defmncur         !<number of months with water deficit for current year
  real, dimension(ilg), intent(inout) :: tcurm            !<temperature of the current month (c)
  real, dimension(ilg), intent(inout) :: annpcp           !<annual precipitation (mm)
  real, dimension(ilg), intent(inout) :: dry_season_length!<length of the dry season (months)
  real, dimension(ilg), intent(inout) :: twarmm           !<temperature of the warmest month (c)
  real, dimension(ilg), intent(inout) :: tcoldm           !<temperature of the coldest month (c)
  real, dimension(ilg), intent(inout) :: gdd5             !<growing degree days above 5 c
  real, dimension(ilg), intent(inout) :: aridity          !<aridity index, ratio of potential evaporation to precipitation
  real, dimension(ilg), intent(inout) :: srplsmon         !<number of months in a year with surplus water i.e. precipitation more than potential evaporation
  real, dimension(ilg), intent(inout) :: defctmon         !<number of months in a year with water deficit i.e. precipitation less than potential evaporation
  real, dimension(12,ilg), intent(inout) :: tmonth        !<monthly temperatures
  real, dimension(ilg), intent(inout) :: anpcpcur         !<annual precipitation for current year (mm)
  real, dimension(ilg), intent(inout) :: anpecur          !<annual potential evaporation for current year (mm)
  real, dimension(ilg), intent(inout) :: gdd5cur          !<growing degree days above 5 c for current year
  real, dimension(ilg), intent(inout) :: srplscur         !<water surplus for the current month
  real, dimension(ilg), intent(inout) :: defctcur         !<water deficit for the current month
  real, dimension(ilg), intent(inout) :: srpcuryr         !<water surplus for the current year
  real, dimension(ilg), intent(inout) :: dftcuryr         !<water deficit for the current year
  real, dimension(ilg), intent(inout) :: anndefct         !<annual water deficit (mm)
  real, dimension(ilg), intent(inout) :: annsrpls         !<annual water surplus (mm)
  real, dimension(ilg,icc), intent(inout) :: stemmass     !<stem mass for each of the ctem pfts, \f$(kg C/m^2)\f$
  real, dimension(ilg,icc), intent(inout) :: rootmass     !<root mass for each of the ctem pfts, \f$(kg C/m^2)\f$
  real, dimension(ilg,iccp2,ignd), intent(inout) :: litrmass   !<litter mass for each of the ctem pfts + bare + LUC product pools, \f$(kg C/m^2)\f$
  real, dimension(ilg,icc), intent(inout) :: gleafmas     !<green leaf mass for each of the ctem pfts, \f$(kg C/m^2)\f$
  real, dimension(ilg,icc), intent(inout) :: bleafmas     !<brown leaf mass for each of the ctem pfts, \f$(kg C/m^2)\f$
  real, dimension(ilg,iccp2,ignd), intent(inout) :: soilcmas   !<soil carbon mass for each of the ctem pfts + bare + LUC product pools, \f$(kg C/m^2)\f$
  real, dimension(ilg,icc), intent(inout) :: ailcg        !<green lai for ctem's 9 pfts
  real, dimension(ilg,ican), intent(inout) :: ailc        !<lumped lai for class' 4 pfts
  real, dimension(ilg,icc,ignd), intent(inout) :: rmatctem!<fraction of roots for each of ctem's 9 pfts in each soil layer
  real, dimension(ilg,ican), intent(inout) :: zolnc       !<lumped log of roughness length for class' 4 pfts
  real, dimension(ilg,icc), intent(inout) :: ailcb        !<brown lai for ctem's 9 pfts. for now we assume only grasses can have brown lai
  real, dimension(ilg), intent(inout) :: vgbiomas         !<grid averaged vegetation biomass, \f$(kg C/m^2)\f$
  real, dimension(ilg), intent(inout) :: gavgltms         !<grid averaged litter mass, \f$(kg C/m^2)\f$
  real, dimension(ilg), intent(inout) :: gavgscms         !<grid averaged soil c mass, \f$(kg C/m^2)\f$
  real, dimension(ilg), intent(inout) :: gavglai          !<grid averaged green leaf area index
  real, dimension(ilg,icc), intent(inout) :: bmasveg      !<total (gleaf + stem + root) biomass for each ctem pft, \f$(kg C/m^2)\f$
  real, dimension(ilg,ican), intent(inout) :: cmasvegc    !<total canopy mass for each of the 4 class pfts. recall that class requires canopy
                                                          !<mass as an input, and this is now provided by ctem. \f$kg/m^2\f$.
  real, dimension(ilg,icp1), intent(inout) :: fcanmx      !<fractional coverage of class' 4 pfts
  real, dimension(ilg,ican), intent(inout) :: alvisc      !<visible albedo for class' 4 pfts
  real, dimension(ilg,ican), intent(inout) :: alnirc      !<near ir albedo for class' 4 pfts
  real, dimension(ilg,icc), intent(inout) :: pstemmass    !<stem mass from previous timestep, is value before fire. used by burntobare subroutine
  real, dimension(ilg,icc), intent(inout) :: pgleafmass   !<root mass from previous timestep, is value before fire. used by burntobare subroutine
  real, dimension(ilg,icc), intent(inout) :: flhrloss     !<fall or harvest loss for deciduous trees and crops, respectively, \f$(kg C/m^2)\f$
  real, dimension(ilg,icc), intent(inout) :: stmhrlos     !<stem harvest loss for crops, \f$(kg C/m^2)\f$
  real, dimension(ilg,icc), intent(inout) :: rothrlos     !<root death as crops are harvested, \f$(kg C/m^2)\f$
  real, dimension(ilg,icc), intent(inout) :: grwtheff     !<growth efficiency. change in biomass per year per unit max. lai (\f$(kg C/m^2)\f$)/(m2/m2),
                                                        !<for use in mortality subroutine
  real, dimension(ilg,icc), intent(inout) :: lystmmas     !<stem mass at the end of last year
  real, dimension(ilg,icc), intent(inout) :: lyrotmas     !<root mass at the end of last year
  real, dimension(ilg,icc), intent(inout) :: tymaxlai     !<this year's maximum lai
  real, dimension(ilg,icc), intent(inout) ::  geremort    !<
  real, dimension(ilg,icc), intent(inout) :: intrmort     !<
  real, dimension(ilg,icc), intent(inout) ::  burnvegf    !<per PFT fraction burned of that PFT's area
  real, dimension(ilg), intent(inout) :: popdin           !<population density \f$(people / km^2)\f$
  real, dimension(ilg), intent(inout) :: Cmossmas         !< moss biomass C pool (kgC/m2)
  real, dimension(ilg), intent(inout) :: litrmsmoss       !< moss litter C pool (kgC/m2)
  real, dimension(ilg), intent(inout) :: peatdep          !< peat depth (m)
  !
  !   outputs
  !
  real, dimension(ilg), intent(out) :: ch4soills          !< Methane uptake into the soil column \f$(mg CH_4 m^{-2} s^{-1})\f$
  real, dimension(ilg), intent(out) :: rml                !<leaf maintenance respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
  real, dimension(ilg), intent(out) :: gpp                !<gross primary productivity
  real, dimension(ilg,icc), intent(out) :: slai           !<storage/imaginary lai for phenology purposes
  real, dimension(ilg,icc), intent(out) :: veghght        !<vegetation height (meters)
  real, dimension(ilg,icc), intent(out) :: rootdpth       !<99% soil rooting depth (meters) both veghght & rootdpth can be used as diagnostics
                                                          !<to see how vegetation grows above and below ground, respectively
  real, dimension(ilg), intent(out) :: npp                !<net primary productivity
  real, dimension(ilg), intent(out) :: nep                !<net ecosystem productivity
  real, dimension(ilg), intent(out) :: hetrores           !<heterotrophic respiration
  real, dimension(ilg), intent(out) :: autores            !<autotrophic respiration
  real, dimension(ilg), intent(out) :: soilresp           !<soil respiration. this includes root respiration and respiration from litter and soil
                                                          !<carbon pools. note that soilresp is different from socres, which is respiration from the soil c pool.
  real, dimension(ilg), intent(out) :: rm                 !<maintenance respiration
  real, dimension(ilg), intent(out) :: rg                 !<growth respiration
  real, dimension(ilg), intent(out) :: nbp                !<net biome productivity
  real, dimension(ilg), intent(out) :: dstcemls1          !<carbon emission losses due to disturbance (fire at present) from vegetation
  real, dimension(ilg), intent(out) :: litrfall           !<total litter fall (from leaves, stem, and root) due to all causes (mortality, turnover, and disturbance)
  real, dimension(ilg), intent(out) :: humiftrs           !<transfer of humidified litter from litter to soil c pool
  real, dimension(ilg), intent(out) :: lucemcom           !<land use change (luc) related combustion emission losses, u-mol co2/m2.sec
  real, dimension(ilg), intent(out) :: lucltrin           !<luc related inputs to litter pool, u-mol co2/m2.sec
  real, dimension(ilg), intent(out) :: lucsocin           !<luc related inputs to soil c pool, u-mol co2/m2.sec
  real, dimension(ilg), intent(out) :: dstcemls3          !<carbon emission losses due to disturbance (fire at present) from litter pool
  real, dimension(ilg), intent(out) :: rms                !<stem maintenance respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
  real, dimension(ilg), intent(out) :: rmr                !<root maintenance respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
  real, dimension(ilg), intent(out) :: litres             !<litter respiration
  real, dimension(ilg), intent(out) :: socres             !<soil carbon respiration
  real, dimension(ilg,icc), intent(out) :: rmsveg         !<
  real, dimension(ilg,icc), intent(out) :: rmrveg         !<
  real, dimension(ilg,icc), intent(out) :: rmlveg         !<
  real, dimension(ilg,icc), intent(out) :: gppveg         !<
  real, dimension(ilg,icc), intent(out) :: nppveg         !<npp for individual pfts,  u-mol co2/m2.sec
  real, dimension(ilg,icc), intent(out) :: rgveg          !<
  real, dimension(ilg,iccp1), intent(out) :: nepveg       !<
  real, dimension(ilg,iccp1), intent(out) :: nbpveg       !<
  real, dimension(ilg,iccp2,ignd), intent(out) :: ltresveg     !<fluxes for each pft: litter respiration for each pft + bare fraction
  real, dimension(ilg,iccp2,ignd), intent(out) :: scresveg     !<soil carbon respiration for the given sub-area in umol co2/m2.s, for ctem's pfts
  real, dimension(ilg,iccp1), intent(out) :: hetrsveg     !<
  real, dimension(ilg,iccp2,ignd), intent(out) :: humtrsvg     !<transfer of humidified litter from litter to soil c pool per PFT.
  real, dimension(ilg,icc), intent(out) :: autoresveg     !<
  real, dimension(ilg,icc), intent(out) :: litrfallveg    !<
  real, dimension(ilg,icc), intent(out) :: roottemp       !<root temperature, k
  real, dimension(ilg,icc), intent(out) :: emit_co2       !<carbon dioxide emitted from biomass burning in g of compound
  real, dimension(ilg,icc), intent(out) :: emit_co        !<carbon monoxide emitted from biomass burning in g of compound
  real, dimension(ilg,icc), intent(out) :: emit_ch4       !<methane emitted from biomass burning in g of compound
  real, dimension(ilg,icc), intent(out) :: emit_nmhc      !<non-methane hydrocarbons emitted from biomass burning in g of compound
  real, dimension(ilg,icc), intent(out) :: emit_h2        !<hydrogen gas emitted from biomass burning in g of compound
  real, dimension(ilg,icc), intent(out) :: emit_nox       !<nitrogen oxides emitted from biomass burning in g of compound
  real, dimension(ilg,icc), intent(out) :: emit_n2o       !<nitrous oxide emitted from biomass burning in g of compound
  real, dimension(ilg,icc), intent(out) :: emit_pm25      !<particulate matter less than 2.5 um in diameter emitted from biomass burning in g of compound
  real, dimension(ilg,icc), intent(out) :: emit_tpm       !<total particulate matter emitted from biomass burning in g of compound
  real, dimension(ilg,icc), intent(out) :: emit_tc        !<total carbon emitted from biomass burning in g of compound
  real, dimension(ilg,icc), intent(out) :: emit_oc        !<organic carbon emitted from biomass burning in g of compound
  real, dimension(ilg,icc), intent(out) :: emit_bc        !<black carbon emitted from biomass burning in g of compound
  real, dimension(ilg,icc), intent(out) :: bterm_veg      !<biomass term for fire probabilty calc
  real, dimension(ilg), intent(out) :: lterm              !<lightning term for fire probabilty calc
  real, dimension(ilg,icc), intent(out) :: mterm_veg      !<moisture term for fire probabilty calc
  real, dimension(ilg), intent(out) :: ch4WetSpec            !<
  real, dimension(ilg), intent(out) :: wetfdyn            !<
  real, dimension(ilg), intent(out) :: ch4WetDyn            !<
  real, dimension(ilg,icc), intent(out) :: cc             !<
  real, dimension(ilg,icc), intent(out) :: mm             !<
  real, dimension(ilg,icc), intent(out) :: lambda         !<
  real, dimension(ilg,icc), intent(out) :: afrleaf        !<allocation fraction for leaves
  real, dimension(ilg,icc), intent(out) :: afrstem        !<allocation fraction for stem
  real, dimension(ilg,icc), intent(out) :: afrroot        !<allocation fraction for root
  real, dimension(ilg,icc), intent(out) :: wtstatus       !<soil water status used for calculating allocation fractions
  real, dimension(ilg,icc), intent(out) :: ltstatus       !<light status used for calculating allocation fractions
  real, dimension(ilg), intent(out) :: burnfrac           !<areal fraction burned due to fire for every grid cell (%)
  real, dimension(ilg,icc), intent(out) :: leaflitr       !<leaf litter fall rate (\f$\mu mol CO_2 m^{-2} s^{-1}\f$). this leaf litter does not
                                                        !<include litter generated due to mortality/fire
  real, dimension(ilg,icc), intent(out) :: smfunc_veg     !<soil moisture dependence on fire spread rate
  real, dimension(ilg,icc), intent(out) :: tltrleaf       !<total leaf litter fall rate (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
  real, dimension(ilg,icc), intent(out) :: tltrstem       !<total stem litter fall rate (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
  real, dimension(ilg,icc), intent(out) :: tltrroot       !<total root litter fall rate (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
  real, dimension(ilg,ican), intent(out) :: paicgat       !<
  real, dimension(ilg,ican), intent(out) :: slaicgat      !<
  real, dimension(ilg,icc), intent(out) :: vgbiomas_veg   !<
  real, dimension(ilg), intent(out) :: armoss             !<autotrophic respiration of moss (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
  real, dimension(ilg), intent(out) :: nppmoss            !<net primary production of moss (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
  real, intent(out) :: reprocost(ilg,icc) !< Cost of making reproductive tissues, only non-zero when NPP is positive (\f$\mu mol CO_2 m^{-2} s^{-1}\f$) 
  real, intent(out) :: blfltrdt(ilg,icc)  !<brown leaf litter generated due to disturbance \f$(kg c/m^2)\f$
  real, intent(out) :: glfltrdt(ilg,icc)  !<green leaf litter generated due to disturbance \f$(kg c/m^2)\f$
  real, intent(out) :: glcaemls(ilg,icc)  !<green leaf carbon emission disturbance losses, \f$kg c/m^2\f$
  real, intent(out) :: blcaemls(ilg,icc)  !<brown leaf carbon emission disturbance losses, \f$kg c/m^2\f$
  real, intent(out) :: rtcaemls(ilg,icc)  !<root carbon emission disturbance losses, \f$kg c/m^2\f$
  real, intent(out) :: stcaemls(ilg,icc)  !<stem carbon emission disturbance losses, \f$kg c/m^2\f$
  real, intent(out) :: ltrcemls(ilg,icc)  !<litter carbon emission disturbance losses, \f$kg c/m^2\f$
  real, intent(out) :: ntchlveg(ilg,icc)  !<fluxes for each pft: Net change in leaf biomass, u-mol CO2/m2.sec
  real, intent(out) :: ntchsveg(ilg,icc)  !<fluxes for each pft: Net change in stem biomass, u-mol CO2/m2.sec
  real, intent(out) :: ntchrveg(ilg,icc)  !<fluxes for each pft: Net change in root biomass, 
                                          !! the net change is the difference between allocation and
                                          !! autotrophic respiratory fluxes, u-mol CO2/m2.sec
  real, intent(out) :: mortLeafGtoB(ilg,icc)  !< Green leaf mass converted to brown due to mortality \f$(kg C/m^2)\f$
  real, intent(out) :: phenLeafGtoB(ilg,icc)  !< Green leaf mass converted to brown due to phenology \f$(kg C/m^2)\f$
  real, intent(out) :: turbLitter(ilg,iccp2,ignd)  !< Litter gains/losses due to turbation [ \f$kg C/m^2\f$ ], negative is a gain.
  real, intent(out) :: turbSoilC(ilg,iccp2,ignd)   !< Soil C gains/losses due to turbation [ \f$kg C/m^2\f$ ], negative is a gain.
  real, intent(out) :: gLeafLandCompChg(ilg,icc)    !< Tracker variable for C movement due to competition and LUC in the green leaf pool  [ \f$kg C/m^2\f$ ], negative is a gain.
  real, intent(out) :: bLeafLandCompChg(ilg,icc)    !< Tracker variable for C movement due to competition and LUC in the brown leaf pool  [ \f$kg C/m^2\f$ ], negative is a gain.
  real, intent(out) :: stemLandCompChg(ilg,icc)   !< Tracker variable for C movement due to competition and LUC in the stem pool  [ \f$kg C/m^2\f$ ], negative is a gain.
  real, intent(out) :: rootLandCompChg(ilg,icc)   !< Tracker variable for C movement due to competition and LUC in the root pool  [ \f$kg C/m^2\f$ ], negative is a gain.
  real, intent(out) :: litterLandCompChg(ilg,iccp2,ignd) !< Tracker variable for C movement due to competition and LUC in the litter pool  [ \f$kg C/m^2\f$ ], negative is a gain.
  real, intent(out) :: soilCLandCompChg(ilg,iccp2,ignd) !< Tracker variable for C movement due to competition and LUC in the soil C pool  [ \f$kg C/m^2\f$ ], negative is a gain.

  ! ---------------------------------------------
  ! Local variables:

  integer i, j, k, n, m
  integer sort(icc)
  real yesfrac_comp(ilg,icc) !<
  real fc(ilg)  !<
  real fg(ilg)  !<
  real term        !<
  real pglfmass(ilg,icc)  !<
  real pblfmass(ilg,icc)  !<
  real pstemass(ilg,icc)  !<
  real protmass(ilg,icc)  !<
  real plitmass(ilg,iccp2)!<
  real psocmass(ilg,iccp2)!<
  real pvgbioms(ilg)      !<
  real pgavltms(ilg)      !<
  real pgavscms(ilg)      !<
  real anveg(ilg,icc)    !<
  real rmveg(ilg,icc)    !<
  real pheanveg(ilg,icc) !<
  real soilrsvg(ilg,iccp2) !<
  real ltrestep(ilg,iccp2,ignd) !<
  real screstep(ilg,iccp2,ignd) !<
  real hutrstep(ilg,iccp2,ignd) !<
  real plitmasspl(ilg,iccp2,ignd) !<
  real psocmasspl(ilg,iccp2,ignd) !<
  real rootlitr(ilg,icc) !<
  real stemlitr(ilg,icc) !<
  real nppvgstp(ilg,icc) !<
  real rmlvgstp(ilg,icc) !<
  real rmsvgstp(ilg,icc) !<
  real rmrvgstp(ilg,icc) !<
  real gppvgstp(ilg,icc) !<
  real stemltrm(ilg,icc) !<
  real rootltrm(ilg,icc) !<
  real glealtrm(ilg,icc) !<
  real stemltdt(ilg,icc)  !<
  real rootltdt(ilg,icc)  !<
  real dscemlv1(ilg,icc)  !<
  real dscemlv2(ilg,icc)  !<
  real add2allo(ilg,icc)  !<
  real repro_cost_g(ilg)  !< Tile-level cost of making reproductive tissues, only non-zero when NPP is positive (\f$\mu mol CO_2 m^{-2} s^{-1}\f$) 
  real rgmoss(ilg)        !< moss growth respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
  real litresmoss(ilg)    !< moss litter respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
  real socres_peat(ilg)   !< heterotrophic repsiration from peat soil (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
  real resoxic(ilg)       !< oxic respiration from peat soil (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
  real resanoxic(ilg)     !< anoxic respiration from peat soil (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
  real litrfallmoss(ilg)  !< moss litter fall (kgC/m2/timestep)
  real ltrestepmoss(ilg)  !< litter respiration from moss (kgC/m2/timestep)
  real humstepmoss(ilg)   !< moss humification (kgC/m2/timestep)
  real pCmossmas(ilg)     !< moss biomass C at the previous time step (kgC/m2)
  real plitrmsmoss(ilg)   !< moss litter C at the previous time step (kgC/m2)
  real nppmosstep(ilg)    !< moss npp (kgC/m2/timestep)
  real socrestep(ilg)     !< heterotrophic respiration from soil (kgC/m2/timestep)
  real hutrstep_g(ilg)    !< grid sum of humification from vascualr litter (kgC/m2/timestep)

  !> Initialize the tracer pool trackers for competition and land use change to zero.
  gLeafLandCompChg = 0.0
  bLeafLandCompChg = 0.0
  stemLandCompChg = 0.0
  rootLandCompChg = 0.0
  litterLandCompChg = 0.0
  soilCLandCompChg = 0.0

  !> Begin calculations

  !> Generate the sort index for correspondence between CTEM pfts and the
  !>  values in the parameter vectors
  sort = genSortIndex()

  if (PFTCompetition .or. lnduseon) then
    !>Store green and brown leaf, stem, and root biomass, and litter and
    !!soil c pool mass in arrays. We use these to track C movements due to
    !! LUC and competition for tracer.
    pglfmass = gleafmas    !<green leaf mass from last time step
    pblfmass = bleafmas    !<brown leaf mass from last time step
    pstemass = stemmass    !<stem mass from last time step
    protmass = rootmass    !<root mass from last time step
    plitmasspl = litrmass  !<litter mass from last time step
    psocmasspl = soilcmas  !< soil C mass from last time step
  end if 
  
  if (PFTCompetition) then
          
    !>Calculate bioclimatic parameters for estimating pfts existence
    call  bioclim (iday,       ta,    precip,  netrad,&
                         1,     il2,    ilg, leapnow, &
                      tcurm, srpcuryr,  dftcuryr,  inibioclim,&
                      tmonth, anpcpcur,   anpecur, gdd5cur,&
                      surmncur, defmncur,  srplscur,defctcur,&
                      twarmm,   tcoldm,      gdd5, aridity,&
                      srplsmon, defctmon, anndefct, annsrpls,&
                      annpcp, dry_season_length )

    if (inibioclim) then
      !>
      !>If first day of year then based on updated bioclimatic parameters
      !!find if pfts can exist or not.
      !!If .not. inibioclim then it is the first year of a run that you do not have the
      !!climatological means already in the CTM file. After one
      !!year inibioclim is set to true and the climatological means
      !!are used from the first year.
      !!
      call existence(iday,            1,         il2, ilg,&
                      sort,       twarmm,     tcoldm,  &    
                      gdd5,      aridity,  srplsmon,   defctmon,  &
                  anndefct,     annsrpls,    annpcp,   pftexist,  &
                  dry_season_length )
      !>
      !!Call competition subroutine which on the basis of previous day's
      !!npp estimates changes in fractional coverage of pfts
      !!
      call competition (iday,     1,        il2,      ilg,& !In
                          nppveg,   dofire, leapnow,& !In
                         pftexist, geremort, intrmort,& !In
                         gleafmas, bleafmas, stemmass, rootmass,& !In
                         litrmass, soilcmas, grclarea,   lambda,& !In
                         burnvegf, sort,  pstemmass, & !In
                         pgleafmass, rmatctem,& !In
                         fcancmx,   fcanmx, vgbiomas, gavgltms,& !In/Out
                         gavgscms, bmasveg,  & !In/Out
                         add2allo,      cc,      mm) !Out 

    end if ! inibioclim
  endif  ! if (PFTCompetition)
  !>
  !!If landuse is on, then implelement luc, change fractional coverages,
  !!move biomasses around, and estimate luc related combustion emission losses.
  !!
  if (lnduseon) then

    do j = 1, icc
     do i = il1, il2
       yesfrac_comp(i,j)=fcancmx(i,j)
     enddo
    enddo

    call luc(    il1,      il2,   ilg,  &
                   grclarea, pfcancmx, nfcancmx,     iday,&
               todfrac,  yesfrac_comp,     .true.,  PFTCompetition,  &
                   leapnow,                               &
                   gleafmas, bleafmas, stemmass, rootmass,&
                   litrmass, soilcmas, vgbiomas, gavgltms,&
                   gavgscms,  fcancmx,   fcanmx,&
                   lucemcom, lucltrin, lucsocin)
  else
    lucemcom = 0.
    lucltrin = 0.
    lucsocin = 0.
  endif !lnduseon
  
  if (PFTCompetition .or. lnduseon) then
    !>Save the change in the C pools for the tracer subroutines.
    gLeafLandCompChg = pglfmass - gleafmas    
    bLeafLandCompChg = pblfmass - bleafmas    
    stemLandCompChg = pstemass - stemmass    
    rootLandCompChg = protmass - rootmass    
    litterLandCompChg = plitmasspl - litrmass 
    soilCLandCompChg = psocmasspl - soilcmas  
  end if 

  !>
  !!Store green and brown leaf, stem, and root biomass, and litter and
  !!soil c pool mass in arrays. knowing initial sizes of all pools and
  !!final sizes at the end of this subroutine, we check for conservation of mass.
  !!
  plitmass = 0.0
  psocmass = 0.0
  do 130 j = 1, icc
    do 140 i = il1, il2
      pglfmass(i,j)=gleafmas(i,j)    !<green leaf mass from last time step
      pblfmass(i,j)=bleafmas(i,j)    !<brown leaf mass from last time step
      pstemass(i,j)=stemmass(i,j)    !<stem mass from last time step
      protmass(i,j)=rootmass(i,j)    !<root mass from last time step
      do k = 1, ignd !FLAG at this stage keep as per pft and per tile. JM Feb8 2016.
        plitmass(i,j)=plitmass(i,j) + litrmass(i,j,k)    !litter mass from last time step
        psocmass(i,j)=psocmass(i,j) + soilcmas(i,j,k)    !soil c mass from last time step
      end do
140     continue
130   continue
!
  do 145 i = il1, il2
    pvgbioms(i)=vgbiomas(i)          !<vegetation biomass from last time step    
    pgavltms(i)=gavgltms(i)          !<litter mass from last time step
    pgavscms(i)=gavgscms(i)          !<soil c mass from last time step
    do j = iccp1, iccp2 ! do over the bare fraction and the LUC pool
      do k = 1, ignd !FLAG at this stage keep as per pft and per tile. JM Feb8 2016.
        plitmass(i,j)=plitmass(i,j) + litrmass(i,j,k)  !litter mass over bare fraction
        psocmass(i,j)=psocmass(i,j) + soilcmas(i,j,k)  !soil c mass over bare fraction
      end do
    end do

    pCmossmas(i)  = Cmossmas(i)
    plitrmsmoss(i)= litrmsmoss(i)
145   continue

  ! Find the canopy covered fraction and the bare fraction of the tiles:
  do i = il1, il2
    fc(i) = sum(fcancmx(i,:))
    fg(i)= 1.0 - fc(i)
  end do

  !     ------------------------------------------------------------------
  !>Initialization ends

  !>Autotrophic respiration
  !!
  !!Leaf respiration is calculated in phtsyn subroutine, while stem
  !!and root maintenance respiration are calculated here. We use air 
  !!temperature as a surrogate for  stem temperature
  !!
  !!Find stem and root maintenance respiration in umol co2/m2/sec

  call mainres (fcancmx, fc, stemmass, rootmass, & !In
                il1, il2, ilg, leapnow, & !In
               ta, tbar, rmatctem,& !In
               sort, isand,& !In
               rmsveg, rmrveg, roottemp) ! Out

  ! NOTE: This next bit is a little tricky. Remember ancgveg is the net photosynthesis 
  ! so it is ancgveg = gpp - rmlcgveg. Here we assign the daily mean net photosynthesis
  ! (ancgveg) and mean rml (rmlveg) to the anveg and rmlveg variables 
  ! and their sum to gpp veg since gpp = anveg + rml if we both have some coverage of the 
  ! PFT (fcancmx > 1) and the leaves are not imaginary (lfstatus /= 4).
  ! If no real leaves are in existence we leave rml, anveg and gpp set to 0.
  ! However, if we have real leaves, but they are quite small we are going to
  ! reduce rml, but we need to use the original rml to find the gppveg otherwise 
  ! the gpp will not be corrected properly. We store the original ancgveg for 
  ! phenology to test if the leaves should be coming out.

  pheanveg = 0. 
  rmlveg=0.
  anveg=0.
  gppveg=0.
  do 180 j = 1, icc
    do 190 i = il1, il2
      
      if (fcancmx(i,j) > zero) then
        
        pheanveg(i,j) = ancgveg(i,j) ! to be used for phenology purposes
        
        if (lfstatus(i,j) /= 4) then ! real leaves so use values

          anveg(i,j) = ancgveg(i,j)
          rmlveg(i,j) = rmlcgveg(i,j)
          gppveg(i,j) = anveg(i,j) + rmlveg(i,j)

          if(slai(i,j) > ailcg(i,j))then
            term = ((1.0 / kn(sort(j))) * (1.0 - exp(-kn(sort(j)) * ailcg(i,j))) &
                  / (1.0 / kn(sort(j))) * (1.0 - exp(-kn(sort(j)) * slai(i,j))))
            rmlveg(i,j) = rmlveg(i,j) * term
          end if
        !else 
          ! the leaves were imaginary so leave variables set to the initialized 0 value.
        endif
      end if
190     continue
180   continue

  !! Find total maintenance respiration values and net primary productivity.
  rml(:) = 0.
  rms(:) = 0.
  rmr(:) = 0.
  rm(:) = 0.
  rg(:) = 0.
  npp(:) = 0.
  gpp(:) = 0.
  do 270 j = 1, icc
    do 280 i = il1, il2
      rmveg(i,j)  = rmlveg(i,j) + rmrveg(i,j) + rmsveg(i,j)
      nppveg(i,j) = gppveg(i,j) - rmveg(i,j)

      !>Now that we know maintenance respiration from leaf, stem, and root
      !>and gpp, we can find growth respiration for each vegetation type

      if( nppveg(i,j).gt.zero ) then
          rgveg(i,j)=grescoef(sort(j))*nppveg(i,j)
      else
          rgveg(i,j)=0.0
      endif
      nppveg(i,j) = nppveg(i,j) - rgveg(i,j)

      !>Calculate grid/tile-averaged rates of rm, rg, npp, and gpp
      rml(i)=rml(i)+fcancmx(i,j)*rmlveg(i,j)
      rms(i)=rms(i)+fcancmx(i,j)*rmsveg(i,j)
      rmr(i)=rmr(i)+fcancmx(i,j)*rmrveg(i,j)
      rm(i) =rm(i)+fcancmx(i,j)*rmveg(i,j)
      rg(i) =rg(i)+fcancmx(i,j)*rgveg(i,j)
      npp(i)=npp(i)+fcancmx(i,j)*nppveg(i,j)
      gpp(i)=gpp(i)+fcancmx(i,j)*gppveg(i,j)
      autores(i)=rg(i)+rm(i)
      autoresveg(i,j)=rmveg(i,j) + rgveg(i,j)

  280     continue 
270   continue 
  !
  !>    Add moss GPP and rml to the grid/tile average C fluxes
  !>    for grid cells which have peatlands
  !
  do 335 i = il1, il2
    if (ipeatland(i) > 0) then
      rgmoss(i) = anmoss(i)*grescoefmoss
      rml(i)= rml(i) + rmlmoss(i)
      rm(i) = rm(i)  + rmlmoss(i)
      rg(i) = rg(i)  + rgmoss (i)
      armoss(i) = rmlmoss(i) + rgmoss(i)
      nppmoss(i) = anmoss(i) - rgmoss(i)
      npp(i)= npp(i) + nppmoss (i)
      gpp(i)= gpp(i) + gppmoss(i)
      autores(i) = autores(i) + armoss(i)
    endif
335   continue
  !
  !     autotrophic respiration part ends
  !
  !     ------------------------------------------------------------------
  !>
  !!     Heterotrophic respiration part starts

  !! Find heterotrophic respiration rates (umol co2/m2/sec) 

  ! CAUTION says Vivek
  ! Note that ipeatland is passed to hetresv to calculate psi (matric potential)
  ! in a different way if ipeatland(i)=1. psi is used to calculated soil moisture
  ! dependence scalars for litter and soil cabron respiration in hetresv. Yet, 
  ! soil carbon respiration for gridcells/tiles which are peatlands is overwritten
  ! below in loop 380 by socres_peat(i) making passing of ipeatland to hetresv and
  ! calculation of psi in a different way useless.
  call    hetresv ( fcancmx, fc, litrmass(:,1:icc,:), soilcmas(:,1:icc,:),& ! In
                    delzw, thpor, il1, il2, & ! In
                    ilg,   tbar, psisat, thliq,& ! In
                    sort,bi,  & ! In
                    isand, thice, ipeatland, & ! In
                    ltresveg(:,1:icc,:), scresveg(:,1:icc,:)) ! Out
  
  !! Find heterotrophic respiration rates from bare ground subarea
  call  hetresg  (litrmass(:,iccp1,:),soilcmas(:,iccp1,:), delzw,thpor,& ! In
                    il1,       il2,       ilg,   tbar,    & ! In
                 psisat,        bi,    thliq,      & ! In  
                 thice,        fg,     isand,    & ! In
                 ltresveg(:,iccp1,:),  scresveg(:,iccp1,:)) ! Out

  !>  Find heterotrophic respiration rates for the LUC litter and soilc pools
  !!  The LUC litter and soil C respiration rates are assumed to
  !!  be applied over the entire tile but kept in layer 1  
  call  hetresg  (litrmass(:,iccp2,:),soilcmas(:,iccp2,:), delzw,thpor,& ! In
                      il1,       il2,       ilg,   tbar,    & ! In
                   psisat,        bi,    thliq,   & ! In
                   thice,        fg,     isand,    & ! In
                   ltresveg(:,iccp2,:),  scresveg(:,iccp2,:)) ! Out

  !! When peatlands are simulated find the peat heterotrophic respiration.
  !! This is called even when there are no peatlands as it just sets the outputs 
  !! to zero and exits.
  call hetres_peat       (    il1,          il2,          ilg,   ipeatland,   & ! In
                          isand,   litrmsmoss,      peatdep,      wtable,   & ! In
                           tbar,        thliq,        thice,       thpor,   & ! In
                             bi,        zbotw,        delzw,      psisat,   & ! In
                     litresmoss,  socres_peat,      resoxic,   resanoxic) ! Out

  !>
  !!Find vegetation averaged litter and soil c respiration rates
  !!using values from canopy over ground and canopy over snow subareas
  !!
  hetrsveg(:,:)= 0.0
  do 340 j = 1, icc
    do 350 i = il1, il2     
      if(fcancmx(i,j) > zero) then
        do k = 1,ignd
          ! hetrsveg is kept per PFT and tile (not per layer) at the moment.
          hetrsveg(i,j) =  hetrsveg(i,j) + ltresveg(i,j,k) + scresveg(i,j,k)
        end do
      endif
      nepveg(i,j)=nppveg(i,j)-hetrsveg(i,j)

350     continue
340   continue
  !>
  !!Find litter and soil c respiration rates averaged over the bare
  !!fraction of the grid cell using values from ground and snow over ground sub-areas.
  !!
  do 355 i = il1, il2
    if(fg(i) > zero) then
      do k = 1,ignd
	      ! hetrsveg is kept per PFT and tile (not per layer) at the moment.
        hetrsveg(i,iccp1) = hetrsveg(i,iccp1) + ltresveg(i,iccp1,k) + scresveg(i,iccp1,k)
      end do
      nepveg(i,iccp1)=0.-hetrsveg(i,iccp1)
    endif
355   continue
  !>
  !!Find grid averaged litter and soil c respiration rates
  !!
  litres(:) = 0.
  socres(:) = 0.
  do 360 j = 1,icc
    do 370 i = il1, il2
      do k =1, ignd
        litres(i)=litres(i)+fcancmx(i,j)*ltresveg(i,j,k)
        socres(i)=socres(i)+fcancmx(i,j)*scresveg(i,j,k)
      end do
370     continue
360   continue
  !>
  !! Add moss and peat soil respiration to the grid
  !! In loop 360/370 we have aggregated across all pfts for litres and socres, In loop 380 we add
  !! bareground (iccp1) and LUC pool (iccp2) values to the grid sum if it's not peatland. If it is a peatland, we add litresmoss to
  !! the grid sum but no bareground values as we assume peatlands have no bareground.

  nep(:) = 0.
  do 380 i = il1, il2
    if (ipeatland(i) == 0) then
      do k = 1, ignd
        litres(i) = litres(i) + fg(i) * ltresveg(i,iccp1,k)
        socres(i) = socres(i) + fg(i) * scresveg(i,iccp1,k)
      end do
    else  ! peatlands
      litres(i) = litres(i)+ litresmoss(i) !add the moss litter, which is assumed to cover whole tile.
      !
      ! CAUTION says Vivek
      ! Note that the following line overwrites socres(i) calculated above in loop 370, although
      ! socres(i) is based on scresveg(i,j) which isn't modified. Similarly, the implementation
      ! of peatlands also means hetrores(i), calculated below, is now inconsistent with hetrsveg(i,j).
      ! The implementation of peatlands is based on overwriting variables without making things
      ! consistent. Of course, overwriting variables is not advised because it makes things confusing. 
      !
      socres(i) = socres_peat(i) ! since this is only peatland on this tile, use just the peat value.
    endif

    hetrores(i)= litres(i)+socres(i)
    nep(i)=npp(i)-hetrores(i)

380   continue
  
  !> Update the litter and soil c pools based on litter and soil c respiration rates
  !! found above. also transfer humidified litter to the soil c pool.
  !!
  do 420 j = 1, iccp2
    do 430 i = il1, il2
      do 435 k = 1, ignd

        !> Convert u mol co2/m2.sec -> \f$(kg C/m^2)\f$ respired over the model time step
        ltrestep(i,j,k) = ltresveg(i,j,k) * deltat / 963.62
        screstep(i,j,k) = scresveg(i,j,k) * deltat / 963.62

        !> Update litter and soil c pools
        if (j < iccp1) then
          litrmass(i,j,k) = litrmass(i,j,k) - (ltrestep(i,j,k) * (1.0+humicfac(sort(j))))
          hutrstep(i,j,k) = humicfac(sort(j)) * ltrestep(i,j,k)
        else
          !> Next we add bareground and LUC pool litter mass and humification for non-peatlands.
          if (ipeatland(i) == 0) then
              litrmass(i,j,k)=litrmass(i,j,k)-(ltrestep(i,j,k)*(1.0+humicfac_bg))
              hutrstep(i,j,k)= humicfac_bg * ltrestep(i,j,k)
          ! else for peatlands:
            ! In peatlands there is no bareground litter mass since it is the moss layer.
          endif
        endif

        humtrsvg(i,j,k) = hutrstep(i,j,k)*(963.62/deltat) ! u-mol co2/m2.sec

        soilcmas(i,j,k)=soilcmas(i,j,k) + real(spinfast) * (hutrstep(i,j,k) -  screstep(i,j,k))

        if(litrmass(i,j,k) < zero) litrmass(i,j,k)=0.0
        if(soilcmas(i,j,k) < zero) soilcmas(i,j,k)=0.0
435   continue
430 continue
420 continue

  !>Estimate soil respiration. this is sum of heterotrophic respiration and root maintenance respiration.
  soilrsvg(:,:) = 0.
  do 440 j = 1, icc
    do 445 i = il1, il2
      do 450 k = 1, ignd
        ! soilrsvg kept as per pft/per tile for now (not per layer)
        soilrsvg(i,j) = soilrsvg(i,j) + ltresveg(i,j,k) + scresveg(i,j,k)
450   continue
      soilrsvg(i,j) = soilrsvg(i,j) + rmrveg(i,j)
445 continue
440 continue

  !> But over the bare fraction and LUC product pool there is no live root.

  do 460 i = il1, il2
    do 465 k = 1, ignd
      soilrsvg(i,iccp1) = soilrsvg(i,iccp1) + ltresveg(i,iccp1,k)+scresveg(i,iccp1,k)
465 continue
460 continue

  !> Find grid averaged humification and soil respiration rates
  soilresp(:) = 0.0
  humiftrs(:) = 0.0
  hutrstep_g(:) = 0.0
  do 470 i = il1, il2
    do 480 j = 1,icc
      soilresp(i)=soilresp(i)+fcancmx(i,j)*soilrsvg(i,j)      
    do k = 1,ignd
      hutrstep_g(i) = hutrstep_g(i)+fcancmx(i,j)*hutrstep(i,j,k) 
      humiftrs(i)=humiftrs(i)+fcancmx(i,j)*humtrsvg(i,j,k)
    end do
480     continue

    !> After aggregation of humification and soil respiration rates for non-peatlands
    !! to the grid/tile level, the same must be done for bareground (iccp1).
    !! For peatlands, we additionally add moss values to the grid (litter respiration
    !! and moss root respiration). Note in loop 430 iccp1 is passed for peatlands

    if (ipeatland(i) ==0 ) then !non peatland
      
      soilresp(i) = soilresp(i) + fg(i) * soilrsvg(i,iccp1)
      do k = 1,ignd
        humiftrs(i) = humiftrs(i) + fg(i) * humtrsvg(i,iccp1,k)
      end do

      !Set all peatland vars to 0.
      litrfallmoss(i)= 0.
      ltrestepmoss(i)= 0.
      nppmosstep(i)= 0.
      humstepmoss(i)= 0.
      socrestep(i) = 0. 
      
    else !peatland
   
      ! CAUTION says Vivek
      ! Here again soilresp(i) is overwritten with socres(i)=socres_peat(i) as calculated
      ! above in loop 380. This makes soilresp(i) inconsistent with soilrsvg(i,j) for
      ! peatland gridcells/tile.      
      
      soilresp(i) = socres(i)+litres(i)+rmr(i) !moss root and litter respiration. No bareground!

      ! Calculate moss time step C fluxes, '/365*deltat' convert year-1
      ! to time step-1, 'deltat/963.62' convert umol CO2/m2/s to kgC/m2/deltat.
      ! note that hutrstep_g aggregation for icc was done in loop 480

      litrfallmoss(i)= Cmossmas(i)*rmortmoss/365*deltat !kgC/m2/day(dt)
      ltrestepmoss(i)= litresmoss(i)*(1.0/963.62)*deltat   !kgC/m2/dt
      nppmosstep(i)= nppmoss(i)*(1.0/963.62)*deltat    !kgC/m2/dt
      socrestep(i) = socres(i)*(1.0/963.62)*deltat     !kgC/m2/dt
      soilresp(i)  = soilresp(i)*(1.0/963.62)*deltat   !kgC/m2/dt
      humstepmoss(i)= humicfacmoss * ltrestepmoss(i)        !kgC/m2/dt
      hutrstep_g(i)= hutrstep_g(i) + humstepmoss(i)     !kgC/m2/dt
      humiftrs(i)  = humiftrs(i)+humstepmoss(i)*(963.62/deltat)!umol/m2/s
      
    end if          
470 continue
  !
  !     heterotrophic respiration part ends
  !
  !     ------------------------------------------------------------------

  !> Find CH4 wetland area (if not prescribed) and emissions:
  call  wetland_methane (hetrores,       il1,       il2,      ilg,  & !In
                           wetfrac,    thliq,   currlat,     sand,  &   !In
                         slopefrac,        ta,                       & !In
                        ch4WetSpec,   wetfdyn, ch4WetDyn) ! Out

  !> Calculate the methane that is oxidized by the soil sink
  call soil_ch4uptake(        il1,       il2,       ilg,     tbar,  & !In
                                bi,    thliq,     thice,  psisat,  & !In
                            fcanmx,    wetfdyn, wetfrac,  & !In
                             isand,   ch4conc,  thpor,& !In
                         ch4soills) ! Out

  !> Estimate allocation fractions for leaf, stem, and root components.
  call allocate (lfstatus,   thliq,    ailcg,     ailcb,& !In
                          il1,        il2,       ilg,       sand,  & !In
                         clay,   rmatctem,  gleafmas,   stemmass,  & !In
                     rootmass,       sort,    fcancmx,  & !In
                        isand,       THFC,       THLW,  & !In
                          afrleaf,  afrstem,  afrroot, & !Out
                         wtstatus, ltstatus) !Out

  !>
  !! Estimate fraction of npp that is to be used for horizontal
  !! expansion (lambda) during the next day (i.e. this will be determining
  !! the colonization rate in competition).
  !!
  if (PFTCompetition)  lambda = expansion(il1,il2,ilg,sort,ailcg,lfstatus,nppveg,pftexist)

  !> Maintenance respiration also reduces leaf, stem, and root biomass.
  !! when npp for a given pft is positive then this is taken care by
  !! allocating +ve npp amongst the leaves, stem, and root component.
  !! when npp for a given pft is negative then maintenance respiration
  !! loss is explicitly deducted from each component.

  do 600 j = 1, icc
    do 610 i = il1, il2

        !! Convert npp and maintenance respiration from different components
        !! from units of u mol co2/m2.sec -> \f$(kg C/m^2)\f$ sequestered 
        !! or respired over the model time step (deltat)

        gppvgstp(i,j)=gppveg(i,j)*(1.0/963.62)*deltat !+ add2allo(i,j)

        !> Remove the cost of making reproductive tissues. This cost can
        !! only be removed when NPP is positive.
        reprocost(i,j) =max(0.,nppveg(i,j)*repro_fraction)

        ! Not in use. We now use a constant reproductive cost as the prior formulation
        ! produces perturbations that do not allow closing of the C balance. JM Jun 2014.
        ! nppvgstp(i,j)=nppveg(i,j)*(1.0/963.62)*deltat*(1.-lambda(i,j))
        !     &                  + add2allo(i,j)
        nppvgstp(i,j)=(nppveg(i,j)-reprocost(i,j))*(1.0/963.62)*deltat

        ! Amount of c related to horizontal expansion
        ! Not in use. JM Jun 2014
        ! expbalvg(i,j)=-1.0*nppveg(i,j)*deltat*lambda(i,j)+ add2allo(i,j)*(963.62/1.0)
        
        rmlvgstp(i,j)=rmlveg(i,j)*(1.0/963.62)*deltat
        rmsvgstp(i,j)=rmsveg(i,j)*(1.0/963.62)*deltat
        rmrvgstp(i,j)=rmrveg(i,j)*(1.0/963.62)*deltat
        
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
        else  !>i.e. if lfstatus.eq.4
          
          !> And since we do not have any real leaves on then we do not take into
          !! account co2 uptake by imaginary leaves in carbon budget. rmlvgstp(i,j)
          !! should be zero because we set maintenance respiration from storage/imaginary 
          !! leaves equal to zero. in loop 180
          !!
          ntchlveg(i,j)=-rmlvgstp(i,j)
          ntchsveg(i,j)=-rmsvgstp(i,j)
          ntchrveg(i,j)=-rmrvgstp(i,j)
          
          !> Since no real leaves are on, make allocation fractions equal to zero.
          
          afrleaf(i,j)=0.0
          afrstem(i,j)=0.0
          afrroot(i,j)=0.0
        endif

        gleafmas(i,j)=gleafmas(i,j)+ntchlveg(i,j)
        stemmass(i,j)=stemmass(i,j)+ntchsveg(i,j)
        rootmass(i,j)=rootmass(i,j)+ntchrveg(i,j)

        if(gleafmas(i,j).lt.0.0)then
          write(6,1900)'gleafmas < zero at i=',i,' for pft=',j,''
          write(6,1901)'gleafmas = ',gleafmas(i,j)
          write(6,1901)'ntchlveg = ',ntchlveg(i,j)
          write(6,1902)'lfstatus = ',lfstatus(i,j)
          write(6,1901)'ailcg    = ',ailcg(i,j)
          write(6,1901)'slai     = ',slai(i,j)
  1900    format(a23,i4,a10,i2,a1)
  1902    format(a11,i4)
          call xit ('ctem',-2)
        endif

        if(stemmass(i,j).lt.0.0)then
          write(6,1900)'stemmass < zero at i=(',i,') for pft=',j,')'
          write(6,1901)'stemmass = ',stemmass(i,j)
          write(6,1901)'ntchsveg = ',ntchsveg(i,j)
          write(6,1902)'lfstatus = ',lfstatus(i,j)
          write(6,1901)'rmsvgstp = ',rmsvgstp(i,j)
          write(6,1901)'afrstem  = ',afrstem(i,j)
          write(6,1901)'gppvgstp = ',gppvgstp(i,j)
          write(6,1901)'rmsveg = ',rmsveg(i,j)
  1901    format(a11,f12.8)
          call xit ('ctem',-3)
        endif

        if(rootmass(i,j).lt.0.0)then
          write(6,1900)'rootmass < zero at i=(',i,') for pft=',j,')'
          write(6,1901)'rootmass = ',rootmass(i,j)
          call xit ('ctem',-4)
        endif

        !! Convert net change in leaf, stem, and root biomass into
        !! u-mol co2/m2.sec for use in balcar subroutine
        !!
        ntchlveg(i,j) = ntchlveg(i,j) * (963.62 / deltat)
        ntchsveg(i,j) = ntchsveg(i,j) * (963.62 / deltat)
        ntchrveg(i,j) = ntchrveg(i,j) * (963.62 / deltat)

        !! To avoid over/underflow problems set gleafmas, stemmass, and
        !! rootmass to zero if they get too small
        !!
        if(bleafmas(i,j) < zero) bleafmas(i,j) = 0.0
        if(gleafmas(i,j) < zero) gleafmas(i,j) = 0.0
        if(stemmass(i,j) < zero) stemmass(i,j) = 0.0
        if(rootmass(i,j) < zero) rootmass(i,j) = 0.0
      
610 continue
600 continue
  
  !> Calculate grid averaged value of C related to spatial expansion
  
  repro_cost_g(:)=0.0    !< amount of C for production of reproductive tissues
  do 620 j = 1,icc
    do 621 i = il1, il2
      repro_cost_g(i) = repro_cost_g(i) + fcancmx(i,j) * reprocost(i,j)
621 continue
620 continue

  !    ------------------------------------------------------------------

  !! The phenology subroutine determines leaf status for each pft and calculates leaf litter.
  !! the phenology subroutine uses soil temperature (tbar) and root temperature. however,
  !! since CTEM doesn't make the distinction between canopy over ground, and canopy over
  !! snow sub-areas for phenology purposes (for example, leaf onset is not assumed to occur
  !! at different times over these sub-areas) we use average soil and root temperature in 
  !! the phenology subroutine.
  call phenolgy(il1,  il2,  ilg,    leapnow,  tbar,  thice,    &!In 
                   thliq,     THLW,     THFC,       ta,&!In 
                 pheanveg,     iday,     radj, roottemp,&!In 
                 rmatctem, stemmass, rootmass,     sort,&!In 
                 fcancmx,  isand, &!In 
                 lfstatus,  pandays, colddays, gleafmas, bleafmas, & !In/Out
                 flhrloss, leaflitr, phenLeafGtoB ) ! Out
                         
  !> While leaf litter is calculated in the phenology subroutine, stem
  !! and root turnover is calculated in the turnoverStemRoot subroutine.
  call turnoverStemRoot (stemmass, rootmass,  lfstatus,    ailcg,& !In 
                   il1,       il2,       ilg,  leapnow,  &!In
                              sort, fcancmx,&!In
                          stmhrlos, rothrlos,& !In/Out
                          stemlitr, rootlitr) !Out

  !> Update green leaf biomass for trees and crops, brown leaf biomass for grasses,
  !! stem and root biomass for litter deductions, and update litter pool with leaf
  !! litter calculated in the phenology subroutine and stem and root litter 
  !! calculated in the turnoverStemRoot subroutine. Also add the reproduction
  !!  carbon directly to the litter pool. We only add to non-perennially frozen soil
  !! layers so first check which layers are unfrozen and then do the allotment 
  !! appropriately. For defining which layers are frozen, we use the active layer depth.
  call updatePoolsTurnover(il1, il2, ilg, reprocost, rmatctem, & !In
                          stemmass, rootmass, litrmass, rootlitr,& !In/Out
                         gleafmas, bleafmas, leaflitr, stemlitr) !In/Out

  !>
  !> Call the mortaliy subroutine which calculates mortality due to reduced growth and aging. 
  !! Exogenous mortality due to fire and other disturbances and the subsequent litter 
  !! that is generated is calculated in the disturb subroutine.
  !!
  !! Set maxage >0 in classic_params.f90 to switch on mortality due to age and
  !! reduced growth. Mortality is linked to the competition parameterization and generates bare fraction.
  call mortalty (stemmass,   rootmass,    ailcg,   gleafmas,  & !In
                    bleafmas,        il1,      il2,        ilg,  & !In
                     leapnow,       iday,     sort,    fcancmx,  & !In
                    lystmmas,   lyrotmas, tymaxlai,   grwtheff,  & !In/Out
                    stemltrm,   rootltrm, glealtrm,   geremort,  & !Out
                    intrmort) !Out

  !> Update leaf, stem, and root biomass pools to take into loss due to mortality, and put the
  !! litter into the litter pool. the mortality for green grasses doesn't generate litter, instead they turn brown.
  call updatePoolsMortality(il1, il2, ilg, stemltrm, rootltrm, & ! In 
                            rmatctem, & !In
                            stemmass, rootmass, litrmass, & !In/Out
                            glealtrm, gleafmas, bleafmas,& !In/Out
                            mortLeafGtoB) ! Out
  !    ------------------------------------------------------------------
  !>
  !> Call the disturbance subroutine which calculates mortality due to fire and other disturbances.
  !> the primary output from from disturbance subroutine is litter generated, c emissions due to fire
  !! and area burned, which may be used to estimate change in fractional coverages.
  !!
  !! Disturbance is spatial and requires area of gcm grid cell and areas of different pfts present in
  !! a given grid cell. however, when ctem is operated at a point scale then it is assumed that the
  !! spatial scale is 1 hectare = 10,000 m2. the disturbance subroutine may be stopped from simulating
  !! any fire by specifying fire extingushing probability equal to 1.
  call disturb (stemmass, rootmass, gleafmas, bleafmas,& !In
                           thliq,    THLW,      THFC,    uwind,& !In
                            vwind,  lightng,  fcancmx, litrmass,& !In
                         rmatctem,     ilg,  il1,      il2,     sort, & !In
                         grclarea,   thice,   popdin, lucemcom,& !In
                           dofire,  currlat,     iday, fsnow,& !In
                            isand,  & !In
                         stemltdt, rootltdt, glfltrdt, blfltrdt,& !Out (Primary)
                         glcaemls, rtcaemls, stcaemls,& !Out (Primary)
                         blcaemls,   ltrcemls,   burnfrac,               & !Out (Primary)
                         pstemmass,  pgleafmass, emit_co2,    emit_ch4,   & !Out (Primary)
                         emit_co,  emit_nmhc,                           & ! Out (Secondary)
                         emit_h2,  emit_nox, emit_n2o, emit_pm25,& ! Out (Secondary)
                         emit_tpm, emit_tc,  emit_oc,  emit_bc,& ! Out (Secondary)
                         burnvegf, bterm_veg,mterm_veg,  lterm,& ! Out (Secondary)
                        smfunc_veg) ! Out (Secondary) 

  !> Calculate NBP (net biome production) for each pft by taking into account
  !! C emission losses. The disturbance routine produces emissions due to fire
  !! and while the land use change subroutine calculates emissions due to LUC. 
  !! The LUC related combustion flux is assumed to be spread uniformly over the
  !! tile as it is no longer associated with any one PFT. To calculate the NBP 
  !! we do not subtract LUC emissions from the PFT-level NBP but we do subtract 
  !! it from the per tile NBP. 
  call calcNBP(il1, il2, ilg,  deltat, nepveg, fcancmx, & !In
                   lucemcom, ltresveg, scresveg, nep, & !In
                   glcaemls, blcaemls, stcaemls, rtcaemls, ltrcemls, & ! In/Out
                   nbpveg, dstcemls1, dstcemls3, nbp) ! Out 

  !> Allow cryoturbation and bioturbation to move the soil C between
  !! layers. Since this is neither consuming nor adding C, this does not
  !! affect our C balance in balcar. There is also an internal C balance check.
  call turbation(il1,il2,delzw,zbotw,isand,maxAnnualActLyr,spinfast, &!In
                litrmass,soilcmas, &! In/Out
                turbLitter, turbSoilC) ! Out

  !> Prepare for the carbon balance check. Calculate total litter fall from each 
  !! component (leaves, stem, and root) from all causes (normal turnover, drought
  !! and cold stress for leaves, mortality, and disturbance), calculate grid-average
  !! vegetation biomass, litter mass, and soil carbon mass, and litter fall rate.
  !! Also add the bare ground values to the grid-average. If a peatland, we assume no bareground and
  !! add the moss values instead. Note: peatland soil C is not aggregated from plants but updated
  !! by humification and respiration from the previous stored value
  call prepBalanceC(il1, il2, ilg, fcancmx, glealtrm, glfltrdt, &  ! In 
                        blfltrdt, stemltrm, stemltdt, rootltrm, rootltdt, & ! In
                        ipeatland, nppmosstep, pgavscms, humstepmoss, & ! In
                        ltrestepmoss, stemlitr, rootlitr, rootmass, & ! In
                        litrmass, soilCmas, hutrstep_g, stemmass, bleafmas, & ! In
                        gleafmas, socrestep, fg, litrfallmoss, & ! In
                        leaflitr, Cmossmas, litrmsmoss,  & ! In/Out
                        tltrleaf, tltrstem, tltrroot, & ! Out 
                        vgbiomas, litrfall, gavgltms, litrfallveg, &! Out 
                        gavgscms, vgbiomas_veg) ! Out 
  
  !>At this stage we have all required fluxes in u-mol co2/m2.sec and initial (loop 140 and 145)
  !!and updated sizes of all pools (in \f$(kg C/m^2)\f$). Now we call the balcar subroutine and make sure
  !!that C in leaves, stem, root, litter and soil C pool balances within a certain tolerance.
  if (spinfast == 1) call  balcar(gleafmas, stemmass, rootmass,  bleafmas,&
                                  litrmass, soilcmas, ntchlveg,  ntchsveg,&
                                  ntchrveg, tltrleaf, tltrstem,  tltrroot,&
                                  glcaemls, blcaemls, stcaemls,  rtcaemls,&
                                  ltrcemls, ltresveg, scresveg,  humtrsvg,&
                                  pglfmass, pblfmass, pstemass,  protmass,&
                                  plitmass, psocmass, vgbiomas,  reprocost,&
                                  pvgbioms, gavgltms, pgavltms,  gavgscms,&
                                  pgavscms, dstcemls3, repro_cost_g,       &
                                   autores, hetrores,      gpp,    &
                                    litres,   socres, dstcemls1,    &
                                  litrfall, humiftrs,                 &
                                       il1,      il2,      ilg,     &
                                 ipeatland, Cmossmas, pCmossmas,                &
                                nppmosstep, litrfallmoss, litrmsmoss,&
                               plitrmsmoss, ltrestepmoss, humstepmoss)
  
  !>
  !> Finally find vegetation structural attributes which can be passed 
  !! to the land surface scheme using leaf, stem, and root biomass.
  !>
  call allometry( gleafmas, bleafmas, stemmass, rootmass,& !In
                   il1,        il2,        ilg,      zbotw,  & !In
                 delzw,   soildpth,    fcancmx,  & !In
             ipeatland,  maxAnnualActLyr,              & !In
                 ailcg,    ailcb,     ailc,    zolnc,& !Out
                 rmatc, rmatctem,     slai,  bmasveg,& !Out
             cmasvegc,  veghght, rootdpth,   alvisc,& !Out
                alnirc,    paicgat,   slaicgat) !Out

  !> Calculation of gavglai is moved from loop 1100 to here since ailcg is updated by allometry
  gavglai (:) = 0.0
  do j = 1, icc
    do i = il1, il2
      gavglai(i) = gavglai(i) + fcancmx(i,j) * ailcg(i,j)
    end do
  end do
  
  !> At the end of the day, find the depth of the peat, update the degree days for moss photosynthesis and the peat bottom layer depth
  do i = il1, il2
    if (ipeatland(i) > 0) peatdep(i) = peatDepth(gavgscms(i))
  end do  

  call peatDayEnd(il2)

  return
    
end subroutine ctem 
!>@}
    
!>\ingroup ctem_calcnbp
!!@{
!> Calculate NBP (net biome production) for each pft by taking into account
!! C emission losses. The disturbance routine produces emissions due to fire
!! and while the land use change subroutine calculates emissions due to LUC. 
!! The LUC related combustion flux is assumed to be spread uniformly over the
!! tile as it is no longer associated with any one PFT. To calculate the NBP 
!! we do not subtract LUC emissions from the PFT-level NBP but we do subtract 
!! it from the per tile NBP. 
!!@author V. Arora, J. Melton

subroutine calcNBP(il1, il2, ilg,  deltat, nepveg, fcancmx, & !In
                  lucemcom, ltresveg, scresveg, nep, & !In
                  glcaemls, blcaemls, stcaemls, rtcaemls, ltrcemls, & ! In/Out
                  nbpveg, dstcemls1, dstcemls3, nbp) ! Out 
  
  use classic_params, only : icc, iccp1, iccp2

  implicit none

  ! arguments
  integer, intent(in) :: il1             !< il1=1
  integer, intent(in) :: il2             !< il2=ilg (no. of grid cells in latitude circle)
  integer, intent(in) :: ilg
  real, intent(in)    :: fcancmx(:,:)    !< max. fractional coverage of ctem's 9 pfts, but this can be
                                         !! modified by land-use change, and competition between pfts
  real, intent(in) :: deltat             !< CTEM (biogeochemical) time step (days)
  real, intent(in) :: nepveg(:,:)        !< Net ecosystem productivity,  \f$\mu mol CO_2 m^{-2} s^{-1}\f$
  real, intent(in) :: lucemcom(:)        !< Land use change (LUC) related combustion emission losses,  \f$\mu mol CO_2 m^{-2} s^{-1}\f$
  real, intent(in) :: ltresveg(:,:,:)    !< Litter respiration for each pft, bare fraction, and LUC product pool, \f$\mu mol CO_2 m^{-2} s^{-1}\f$
  real, intent(in) :: scresveg(:,:,:)    !< Soil carbon respiration for each pft, bare fraction, and LUC product pool, \f$\mu mol CO_2 m^{-2} s^{-1}\f$
  real, intent(in) :: nep(:)             !< Net ecosystem productivity, tile average,  \f$\mu mol CO_2 m^{-2} s^{-1}\f$

  real, intent(inout) ::  glcaemls(:,:)  !< Green leaf carbon emission losses, \f$(kg C/m^2)\f$
  real, intent(inout) ::  blcaemls(:,:)  !< Brown leaf carbon emission losses, \f$(kg C/m^2)\f$
  real, intent(inout) ::  rtcaemls(:,:)  !< Root carbon emission losses, \f$(kg C/m^2)\f$
  real, intent(inout) ::  stcaemls(:,:)  !< Stem carbon emission losses, \f$(kg C/m^2)\f$
  real, intent(inout) ::  ltrcemls(:,:)  !< Litter carbon emission losses, \f$(kg C/m^2)\f$
  
  real, intent(out) :: nbpveg(ilg,iccp1) !< Net biome productivity,  \f$\mu mol CO_2 m^{-2} s^{-1}\f$
  real, intent(out) :: dstcemls1(ilg)    !< grid ave. carbon emission losses due to disturbance, vegetation,  \f$\mu mol CO_2 m^{-2} s^{-1}\f$
  real, intent(out) :: dstcemls3(ilg)    !< grid ave. carbon emission losses due to disturbance, litter,  \f$\mu mol CO_2 m^{-2} s^{-1}\f$
  real, intent(out) :: nbp(ilg)          !< Net biome productivity, tile average,  \f$\mu mol CO_2 m^{-2} s^{-1}\f$
  
  ! Local vars 
  integer :: i, j 
  real :: dscemlv1(ilg,icc)  !< Disturbance emission losses from plants, \f$(kg C/m^2)\f$
  real :: dscemlv2(ilg,icc)  !< Disturbance emission losses from plants and litter, \f$(kg C/m^2)\f$
  real :: dstcemls2(ilg)     !<grid ave. carbon emission losses due to disturbance, total,  \f$\mu mol CO_2 m^{-2} s^{-1}\f$

  !--
  
  do 100 i = il1, il2
    do 101 j = 1, icc
      dscemlv1(i,j) = glcaemls(i,j) + blcaemls(i,j) + stcaemls(i,j) + rtcaemls(i,j)
      dscemlv2(i,j) = dscemlv1(i,j) + ltrcemls(i,j)

      ! Convert \f$(kg C/m^2)\f$ emitted in one day into u mol co2/m2.sec before
      ! subtracting emission losses from nep.
      nbpveg(i,j) = nepveg(i,j) - dscemlv2(i,j)*(963.62/deltat)

101    continue

    ! For accounting purposes, we also need to account for the bare fraction 
    ! NBP. Since there is no fire on the bare, we use 0.
    nbpveg(i,iccp1) = nepveg(i,iccp1)   - 0.

100 continue
  !>
  !! Calculate grid. averaged rate of carbon emissions due to fire in u-mol co2/m2.sec. 
  !! Convert all emission losses from \f$(kg C/m^2)\f$ emitted in 1 day to u-mol co2/m2.sec.
  !! Calculate grid averaged carbon emission losses from litter.
  !!
  dstcemls1 = 0.0
  dstcemls2 = 0.0
    do 103 j = 1,icc
      do 104 i = il1, il2
        dstcemls1(i) = dstcemls1(i) + fcancmx(i,j) * dscemlv1(i,j) * (963.62 / deltat)
        dstcemls2(i) = dstcemls2(i) + fcancmx(i,j) * dscemlv2(i,j) * (963.62 / deltat)
        glcaemls(i,j)= glcaemls(i,j)*(963.62/deltat)
        blcaemls(i,j)=blcaemls(i,j)*(963.62/deltat)
        stcaemls(i,j)=stcaemls(i,j)*(963.62/deltat)
        rtcaemls(i,j)=rtcaemls(i,j)*(963.62/deltat)
        ltrcemls(i,j)=ltrcemls(i,j)*(963.62/deltat)
104    continue
103  continue

  ! For the tile-level NBP, we include the disturbance emissions as well as
  ! respiration from the paper (litter) and furniture (soil carbon) pools (LUC
  ! product pools). Also include here the instantaneous emissions due to LUC.
  nbp(:) = 0.
  dstcemls3 = 0.0
  do 105 i = il1, il2
    nbp(i) = nep(i) - dstcemls2(i) - (ltresveg(i,iccp2,1) + scresveg(i,iccp2,1)) - lucemcom(i)
    dstcemls3(i) = dstcemls2(i) - dstcemls1(i)  !litter is total - vegetation.
105  continue

end subroutine calcNBP
!>@}
!>\namespace ctem

!>Central module that contains the ctem driver and associated subroutines.
!!
!!The basic model structure of CTEM includes three live vegetation components
!!(leaf (L), stem (S) and root (R)) and two dead carbon pools (litter or
!!detritus (D) and soil carbon (H)). The amount of carbon in these pools
!!(\f$C_\mathrm{L}\f$, \f$C_\mathrm{S}\f$, \f$C_\mathrm{R}\f$, \f$C_\mathrm{D}\f$,
!!\f$C_\mathrm{H}\f$, \f$kgC m^{-2}\f$) is tracked prognostically through the
!!fluxes in and out of them. The rate change equations for carbon in these
!!pools are summarized in Sect. \ref{rate_change_eqns} after the processes
!!leading to the calculation of fluxes in and out of these pools are introduced
!!in the following sections.
!!
!>\file
end module ctemDriver