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

      subroutine     ctem( fcancmx,    fsnow,     sand,      clay,  &
     &                  ilg,   il1,      il2,     iday,      radj,  &
     &                          ta,     delzw, ancgveg,   rmlcgveg, &
     &                       zbotw,  &
     &                       uwind,    vwind,  lightng,      tbar,  &
     &                   pfcancmx, nfcancmx,             &
     &                    soildpth, spinfast,   todfrac,&
     &                      netrad,   precip,    psisat,            &
     &                    grclarea,   popdin,    isand,             &
     &                    wetfrac, slopefrac,       bi,             &
     &                      thpor,   currlat,  ch4conc,  &
     &                       THFC,      THLW,     thliq,  thice,    &
     &                  ipeatland,    anmoss,   rmlmoss,  gppmoss,  &
     &                   wtable,                                    &
!
!    ------------- logical switches determining model behaviour 
!
     &               PFTCompetition,  dofire,  lnduseon,  inibioclim,  &
     &                      leapnow,                                   &
!
!    -------------- all inputs used by ctem are above this line ---------
!
     &                    stemmass, rootmass, litrmass,  gleafmas,&
     &                    bleafmas, soilcmas,    ailcg,      ailc,&
     &                       zolnc, rmatctem,    rmatc,     ailcb,&
     &                    flhrloss,  pandays, lfstatus,  grwtheff,&
     &                    lystmmas, lyrotmas, tymaxlai,  vgbiomas,&
     &                    gavgltms, gavgscms, stmhrlos,      slai, &
     &                     bmasveg, cmasvegc, colddays,  rothrlos,&
     &                      fcanmx,   alvisc,   alnirc,   gavglai,&
     &                    Cmossmas, litrmsmoss,     peatdep,               &
!
!    ------------- the following are all competition related variables ---
!
     &                    geremort,   intrmort,   pstemmass,    pgleafmass, &
     &                       tcurm,   srpcuryr,    dftcuryr,      lambda,   &
     &                      tmonth, anpcpcur,  anpecur,   gdd5cur,&
     &                    surmncur, defmncur, srplscur,  defctcur,&
     &                     aridity, srplsmon, defctmon,  anndefct,&
     &                    annsrpls,  annpcp,dry_season_length,&
     &                    pftexist,   twarmm,       tcoldm,         gdd5,   &
!
!    -------------- inputs updated by ctem are above this line ------
!    ------------- these include all prognostic variables -----------
!
!    
!    --------- and finally all output is below this line ------------------
!
!    ---- OUTPUT COMMON TO AGCM AND OFFLINE RUNS ----\
     &                        npp,       nep, hetrores,   autores,&
     &                   soilresp,        rm,       rg,       nbp,&
     &                     litres,    socres,      gpp, dstcemls1,&
     &                   litrfall,  humiftrs,  veghght,  rootdpth,&
     &                        rml,       rms,      rmr,  tltrleaf,&
     &                   tltrstem,  tltrroot, leaflitr,  roottemp,&
     &                   burnfrac,                 lucemcom,    lucltrin,   &
     &                   lucsocin,   dstcemls3,                             &
     &                 ch4WetSpec,  ch4WetDyn,      wetfdyn,   ch4soills,   &
     &                                 paicgat,    slaicgat,                &
   &                    emit_co2,   emit_ch4,                              & 
!    ---- OUTPUT COMMON TO AGCM AND OFFLINE RUNS ----/

!    ---- OUTPUT EXCLUSIVE TO OFFLINE RUNS ----\
     &                   emit_co,   emit_nmhc,  smfunc_veg,                &
     &                    emit_h2,  emit_nox, emit_n2o, emit_pm25,&
     &                    emit_tpm, emit_tc,  emit_oc,    emit_bc,&
     &                  bterm_veg,      lterm,    mterm_veg,  burnvegf,     &
     &                litrfallveg,    humtrsvg,    ltstatus,      nppveg,   &
     &                    afrleaf,     afrstem,     afrroot,    wtstatus,   &
     &                      rmlveg,  rmsveg,   rmrveg,    rgveg,&
     &                vgbiomas_veg,  gppveg,   nepveg,   nbpveg,&
     &                  hetrsveg,autoresveg, ltresveg, scresveg,&
     &                    nppmoss,     armoss,                              &
     &                         cc,         mm                               &
!    ---- OUTPUT EXCLUSIVE TO OFFLINE RUNS ----/
     &                  )

!
!             Canadian Terrestrial Ecosystem Model (CTEM)
!             Main Ctem Subroutine Compatible With CLASS

!     28  Nov 2018  - Clean up argument list before implementation in AGCM
!     V. Arora        
!
!     14  Mar 2016  - Remove grid cell area calculation to the driver. This will
!     J. Melton       harmonize this subroutine with the coupled code.
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

implicit none
!
!     inputs
!
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
!
!     updates
!
logical, intent(inout) :: pftexist(ilg,icc)             !<
logical, intent(inout) :: inibioclim                    !<switch telling if bioclimatic parameters are being initialized from scratch (false)
                                                        !<or being initialized from some spun up values(true).
integer, dimension(ilg,icc), intent(inout) :: pandays   !<days with positive net photosynthesis (an) for use in the phenology subroutine
integer, dimension(ilg,2), intent(inout) :: colddays    !<cold days counter for tracking days below a certain temperature threshold for ndl dcd and crop pfts.
integer, dimension(ilg,icc), intent(inout) :: lfstatus  !<leaf phenology status
! real, dimension(ilg,icc), intent(inout) :: ancsveg      !< net photosynthetic rate for ctems 9 pfts for canopy over snow subarea
real, dimension(ilg,icc), intent(inout) :: ancgveg      !< net photosynthetic rate for ctems 9 pfts for canopy over ground subarea
! real, dimension(ilg,icc), intent(inout) :: rmlcsveg     !< leaf respiration rate for ctems 9 pfts forcanopy over snow subarea
real, dimension(ilg,icc), intent(inout) :: rmlcgveg     !< leaf respiration rate for ctems 9 pfts forcanopy over ground subarea
real, dimension(ilg,icc), intent(inout) :: fcancmx      !< max. fractional coverage of ctem's 9 pfts, but this can be
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
  real, dimension(ilg,iccp2), intent(inout) :: litrmass   !<litter mass for each of the ctem pfts + bare + LUC product pools, \f$(kg C/m^2)\f$
  real, dimension(ilg,icc), intent(inout) :: gleafmas     !<green leaf mass for each of the ctem pfts, \f$(kg C/m^2)\f$
  real, dimension(ilg,icc), intent(inout) :: bleafmas     !<brown leaf mass for each of the ctem pfts, \f$(kg C/m^2)\f$
  real, dimension(ilg,iccp2), intent(inout) :: soilcmas   !<soil carbon mass for each of the ctem pfts + bare + LUC product pools, \f$(kg C/m^2)\f$
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
real, dimension(ilg,iccp2), intent(out) :: ltresveg     !<fluxes for each pft: litter respiration for each pft + bare fraction
real, dimension(ilg,iccp2), intent(out) :: scresveg     !<soil carbon respiration for the given sub-area in umol co2/m2.s, for ctem's pfts
real, dimension(ilg,iccp1), intent(out) :: hetrsveg     !<
real, dimension(ilg,iccp2), intent(out) :: humtrsvg     !<transfer of humidified litter from litter to soil c pool per PFT.
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

! ---------------------------------------------
! Local variables:

integer i
integer j
integer k
integer icount
integer n
integer m
integer sort(icc)
integer k1
integer k2
integer nml
integer ilmos(ilg)
integer jlmos(ilg)
real gppcgveg(ilg,icc)   !<
real yesfrac_comp(ilg,icc) !<
real galtcels(ilg) !<
real dstcemls2(ilg)!<
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
real rmrcgveg(ilg,icc) !<
real anveg(ilg,icc)    !<
real rmveg(ilg,icc)    !<
real pheanveg(ilg,icc) !<
real pancgveg(ilg,icc) !<
real soilrsvg(ilg,iccp2) !<
real ltrestep(ilg,iccp2) !<
real screstep(ilg,iccp2) !<
real hutrstep(ilg,iccp2) !<
real, dimension(ignd) :: unfrzrt        !< root distribution only over unfrozen layers
integer :: botlyr                       !< bottom layer of the unfrozen soil column
real :: frznrtlit                       !< fraction of root distribution in frozen layers
!
real rootlitr(ilg,icc) !<
real stemlitr(ilg,icc) !<
real nppvgstp(ilg,icc) !<
real rmlvgstp(ilg,icc) !<
real rmsvgstp(ilg,icc) !<
real rmrvgstp(ilg,icc) !<
real gppvgstp(ilg,icc) !<
real ntchlveg(ilg,icc) !<
real ntchsveg(ilg,icc) !<
real ntchrveg(ilg,icc) !<
real stemltrm(ilg,icc) !<
real rootltrm(ilg,icc) !<
real glealtrm(ilg,icc) !<
real stemltdt(ilg,icc)  !<
real rootltdt(ilg,icc)  !<
real glfltrdt(ilg,icc)  !<
real blfltrdt(ilg,icc)  !<
real glcaemls(ilg,icc)  !<
real blcaemls(ilg,icc)  !<
real rtcaemls(ilg,icc)  !<
real stcaemls(ilg,icc)  !<
real ltrcemls(ilg,icc)  !<
real dscemlv1(ilg,icc)  !<
real dscemlv2(ilg,icc)  !<
real add2allo(ilg,icc)  !<
real reprocost(ilg,icc) !<
real repro_cost_g(ilg)  !<
real lambdaalt !<
real rgmoss(ilg)        !< moss growth respiration ($\mu mol CO2 m^{-2} s^{-1}$)
real litresmoss(ilg)    !< moss litter respiration ($\mu mol CO2 m^{-2} s^{-1}$)
real socres_peat(ilg)   !< heterotrophic repsiration from peat soil ($\mu mol CO2 m^{-2} s^{-1}$)
real resoxic(ilg)       !< oxic respiration from peat soil ($\mu mol CO2 m^{-2} s^{-1}$)
real resanoxic(ilg)     !< anoxic respiration from peat soil ($\mu mol CO2 m^{-2} s^{-1}$)
real litrfallmoss(ilg)  !< moss litter fall (kgC/m2/timestep)
real ltrestepmoss(ilg)  !< litter respiration from moss (kgC/m2/timestep)
real humstepmoss(ilg)   !< moss humification (kgC/m2/timestep)
real pCmossmas(ilg)     !< moss biomass C at the previous time step (kgC/m2)
real plitrmsmoss(ilg)   !< moss litter C at the previous time step (kgC/m2)
real nppmosstep(ilg)    !< moss npp (kgC/m2/timestep)
real socrestep(ilg)     !< heterotrophic respiration from soil (kgC/m2/timestep)
real hutrstep_g(ilg)    !< grid sum of humification from vascualr litter (kgC/m2/timestep)

!> Begin calculations

!> Generate the sort index for correspondence between 9 pfts and the
!> 12 values in the parameter vectors
  sort = genSortIndex()

if (PFTCompetition) then

!>Calculate bioclimatic parameters for estimating pfts existence

        call  bioclim (iday,       ta,    precip,  netrad,&
     &                    1,     il2,    ilg, leapnow, &
     &                 tcurm, srpcuryr,  dftcuryr,  inibioclim,&
     &                 tmonth, anpcpcur,   anpecur, gdd5cur,&
     &                 surmncur, defmncur,  srplscur,defctcur,&
     &                 twarmm,   tcoldm,      gdd5, aridity,&
     &                 srplsmon, defctmon, anndefct, annsrpls,&
     &                 annpcp, dry_season_length )

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
     &                 sort,       twarmm,     tcoldm,  &    
     &                 gdd5,      aridity,  srplsmon,   defctmon,  &
     &             anndefct,     annsrpls,    annpcp,   pftexist,  &
     &     dry_season_length )
!>
!!Call competition subroutine which on the basis of previous day's
!!npp estimates changes in fractional coverage of pfts
!!
       call competition (iday,     1,        il2,      ilg,&
     &                      nppveg,   dofire, leapnow,&
     &                    pftexist, geremort, intrmort,&
     &                    gleafmas, bleafmas, stemmass, rootmass,&
     &                    litrmass, soilcmas, grclarea,   lambda,&
     &                    burnvegf, sort,  pstemmass, &
     &                    pgleafmass,&
!    ------------------- inputs above this line -------------------
     &                    fcancmx,   fcanmx, vgbiomas, gavgltms,&
     &                    gavgscms, bmasveg,  &
!    ------------------- updates above this line ------------------
     &                    add2allo,      cc,      mm)
!    ------------------- outputs above this line ------------------
!
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

!     ---------------------------------------------------------------
!>
!>initialize required arrays to zero
!>
do 100 i = il1, il2
    rms(i) = 0.0         !<grid ave. stem maintenance respiration
    rmr(i) = 0.0         !<grid ave. root maintenance respiration
    rml(i) = 0.0         !<grid ave. leaf maintenance respiration
    rm(i) = 0.0          !<grid ave. total maintenance respiration
    rg(i) = 0.0          !<grid ave. growth respiration
    npp(i) = 0.0         !<grid ave. net primary productivity
    gpp(i) = 0.0         !<grid ave. gross primary productivity
    nep(i)=0.0           !<grid ave. net ecosystem productivity
    nbp(i)=0.0           !<grid ave. net biome productivity
    litres(i)=0.0        !<grid ave. litter respiration
    socres(i)=0.0        !<grid ave. soil carbon respiration
    hetrores(i)=0.0      !<grid ave. heterotrophic respiration
    autores(i)=0.0       !<grid ave. autotrophic respiration
    soilresp(i)=0.0      !<grid ave. soil respiration
    humiftrs(i)=0.0      !<grid ave. humification rate
    dstcemls1(i)=0.0     !<grid ave. carbon emission losses due to disturbance, vegetation
    dstcemls2(i)=0.0     !<grid ave. carbon emission losses due to disturbance, total
    dstcemls3(i)=0.0     !<grid ave. carbon emission losses due to disturbance, litter
    galtcels(i)=0.0      !<grid ave. litter fire emission losses (redundant, same as dstcemls3) FLAG
    fc(i)=0.0            !<fraction of canopy over ground subarea
    ! fcs(i)=0.0           !<fraction of canopy over snow subarea
    fg(i)=0.0            !<fraction of bare ground subarea
    ! fgs(i)=0.0           !<fraction of snow over ground subarea
    screstep(i,iccp1:iccp2)=0.0  !<soil c respiration in \f$kg c/m^2\f$ over the time step
    ltrestep(i,iccp1:iccp2)=0.0  !<litter c respiration in \f$kg c/m^2\f$ over the time step
    soilrsvg(i,iccp1:iccp2)=0.0  !<soil respiration over the bare fraction
    humtrsvg(i,iccp1:iccp2)=0.0  !<humified rate the bare fraction
    ltresveg(i,iccp1:iccp2)=0.0  !<litter respiration rate over bare fraction
    scresveg(i,iccp1:iccp2)=0.0  !<soil c respiration rate over bare fraction
    hetrsveg(i,iccp1)=0.0  !<heterotrophic resp. rate over bare fraction
    nbpveg(i,iccp1) = 0.0  !<net biome productity for bare fraction
    nepveg(i,iccp1 ) = 0.0  !<net ecosystem productity for bare fraction
    !        expnbaln(i)=0.0        !amount of c related to spatial expansion !Not used JM Jun 2014
    repro_cost_g(i)=0.0    !<amount of C for production of reproductive tissues
100   continue
!
do 110 j = 1,icc
    do 120 i = il1, il2
        ! fcanc(i,j) =0.0
        ! fcancs(i,j)=0.0
        rmsveg(i,j)=0.0    !<stem maintenance resp. rate for each pft
        rmrveg(i,j)=0.0    !<root maintenance resp. rate for each pft
        rmlveg(i,j)=0.0    !<leaf maintenance resp. rate for each pft
        rmveg(i,j)=0.0    !<total maintenance resp. rate for each pft
        rgveg(i,j)=0.0    !<growth resp. rate for each pft
        anveg(i,j)=0.0    !<net photosynthesis rate for each pft
        pheanveg(i,j)=0.0  !<net photosynthesis rate, for phenology purposes
        ! pancsveg(i,j)=0.0  !<net photosynthesis rate, canopy over snow subarea, for phenology purposes
        pancgveg(i,j)=0.0  !<net photosynthesis rate, canopy over ground subarea, for phenology purposes
        gppveg(i,j)=0.0    !<gross primary productity for each pft
        nppveg(i,j)=0.0    !<net primary productity for each pft
        nbpveg(i,j)=0.0    !<net biome productity for each pft
        nepveg(i,j)=0.0    !<net ecosystem productity for each pft
        ltresveg(i,j)=0.0  !<litter respiration rate for each pft
        scresveg(i,j)=0.0  !<soil c respiration rate for each pft
        hetrsveg(i,j)=0.0  !<heterotrophic resp. rate for each pft
        soilrsvg(i,j)=0.0  !<soil respiration rate for each pft
        humtrsvg(i,j)=0.0  !<humification rate for each pft
        litrfallveg(i,j)=0.0 !<litter fall in \f$kg c/m^2\f$ for each pft
        screstep(i,j)=0.0  !<soil c respiration in \f$kg c/m^2\f$ over the tim
        ltrestep(i,j)=0.0  !<litter c respiration in \f$kg c/m^2\f$ over the t
        hutrstep(i,j)=0.0  !<humification rate in \f$kg c/m^2\f$ over the time
        roottemp(i,j)=0.0  !<root temperature
        nppvgstp(i,j)=0.0  !<npp (\f$kg c/m^2\f$) sequestered over the model time step
        gppvgstp(i,j)=0.0  !<gpp (\f$kg c/m^2\f$) sequestered over the model time step
        rmlvgstp(i,j)=0.0  !<leaf maintenance resp. (\f$kg c/m^2\f$) respired over the model time step
        rmsvgstp(i,j)=0.0  !<stem maintenance resp. (\f$kg c/m^2\f$) respired  over the model time step
        rmrvgstp(i,j)=0.0  !<root maintenance resp. (\f$kg c/m^2\f$) respired over the model time step
        ntchlveg(i,j)=0.0  !<net change in gleaf biomass after auto. resp. & allocation
        ntchsveg(i,j)=0.0  !<net change in stem biomass after auto. resp. & allocation
        ntchrveg(i,j)=0.0  !<net change in root biomass after auto. resp. & allocation
        dscemlv1(i,j)=0.0  !<total carbon emission losses (\f$kg c/m^2\f$), mainly due to fire
        dscemlv2(i,j)=0.0  !<total carbon emission losses (\f$kg c/m^2\f$), mainly due to fire
        tltrleaf(i,j)=0.0  !<total leaf litter
        tltrstem(i,j)=0.0  !<total stem litter
        tltrroot(i,j)=0.0  !<total root litter
        vgbiomas_veg(i,j)=0.0 !<vegetation biomass for each pft
        lambda(i,j)=0.0    !< Used to determine the colonization rate
        reprocost(i,j) = 0.0 !< cost of producing reproductive tissues
    !          expbalvg(i,j)=0.0  !amount of c related to spatial expansion !Not used JM Jun 2014
120     continue
110   continue
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
        !do k = 1, ignd !FLAG at this stage keep as per pft and per tile. JM Feb8 2016.
            plitmass(i,j)=plitmass(i,j) + litrmass(i,j)!,k)    !litter mass from last time step
            psocmass(i,j)=psocmass(i,j) + soilcmas(i,j)!,k)    !soil c mass from last time step
        !end do
140     continue
130   continue
!
do 145 i = il1, il2
    pvgbioms(i)=vgbiomas(i)          !<vegetation biomass from last time step
    vgbiomas(i)= 0.0
    pgavltms(i)=gavgltms(i)          !<litter mass from last time step
    gavgltms(i)=0.0
    pgavscms(i)=gavgscms(i)          !<soil c mass from last time step
    gavgscms(i)=0.0
    litrfall(i)=0.0                  !<combined total litter fall rate
    gavglai (i)=0.0                  !<grid averaged green lai
    do j = iccp1, iccp2 ! do over the bare fraction and the LUC pool
        !do k = 1, ignd !FLAG at this stage keep as per pft and per tile. JM Feb8 2016.
            plitmass(i,j)=plitmass(i,j) + litrmass(i,j)!,k)  !litter mass over bare fraction
            psocmass(i,j)=psocmass(i,j) + soilcmas(i,j)!,k)  !soil c mass over bare fraction
        !end do
    end do

    pCmossmas(i)  = Cmossmas(i)
    plitrmsmoss(i)= litrmsmoss(i)
    litrfallmoss(i) = 0.0
    litresmoss(i)   = 0.0
    socres_peat(i)    = 0.0
    resoxic(i)    = 0.0
    resanoxic(i)  = 0.0
    ltrestepmoss(i) = 0.0
    humstepmoss(i) = 0.0
    nppmosstep(i) = 0.0
    rgmoss(i)     = 0.0
    hutrstep_g(i) = 0.0
    nppmoss(i)    = 0.0
    armoss(i)     = 0.0

145   continue

  ! Find the canopy covered fraction and the bare fraction of the tiles:
  do i = il1, il2
    fc(i) = sum(fcancmx(i,:))
    fg(i)= 1.0 - fc(i)
  end do

!     ------------------------------------------------------------------
!>Initialization ends
!>
!>Autotrophic respiration
!!
!!Leaf respiration is calculated in phtsyn subroutine, while stem
  !!and root maintenance respiration are calculated here. We use air 
  !!temperature as a surrogate for  stem temperature
!!
  !!Find stem and root maintenance respiration in umol co2/m2/sec
!!
  call mainres (fcancmx, fc, stemmass, rootmass, & !In
                il1, il2, ilg, leapnow, & !In
               ta, tbar, rmatctem,& !In
               sort, isand,& !In
               rmsveg, rmrveg, roottemp) ! Out

  !call calcNPP()

!>If ailcg/gleafmas is zero, i.e. real leaves are not on, then
!>make maintenance respiration and gpp from storage/imaginary lai
!>equal to zero so that we don't use these numbers in carbon budget.

do 180 j = 1, icc
  do 190 i = il1, il2

    !gppcsveg(i,j)=ancsveg(i,j)+rmlcsveg(i,j)
    gppcgveg(i,j)=ancgveg(i,j)+rmlcgveg(i,j)
!
    if (lfstatus(i,j).eq.4) then
        rmlcgveg(i,j)=0.0
        ! rmlcsveg(i,j)=0.0
        ! pancsveg(i,j)=ancsveg(i,j)   !< to be used for phenology
        pancgveg(i,j)=ancgveg(i,j)   !< purposes
        ! ancsveg(i,j)=0.0
        ancgveg(i,j)=0.0
    else
        ! pancsveg(i,j)=ancsveg(i,j)   !< to be used for phenology
        pancgveg(i,j)=ancgveg(i,j)   !< purposes
        if(slai(i,j).gt.ailcg(i,j))then
            term=((1.0/kn(sort(j)))*(1.0-exp(-kn(sort(j))*ailcg(i,j))) &
    &          /(1.0/kn(sort(j)))*(1.0-exp(-kn(sort(j))* slai(i,j))))
            rmlcgveg(i,j)=rmlcgveg(i,j)*term
            ! rmlcsveg(i,j)=rmlcsveg(i,j)*term
        endif
    endif
190     continue
180   continue
!>
!!Find vegetation averaged leaf, stem, and root respiration, and
!!gpp using values from canopy over ground and canopy over snow subareas
!!
do 270 j = 1, icc
    do 280 i = il1, il2
        ! if( (fcanc(i,j)+fcancs(i,j)).gt.zero) then
        if(fcancmx(i,j) > zero) then  
        !     rmsveg(i,j)= (fcanc(i,j)*rmscgveg(i,j) + &
        ! &        fcancs(i,j)*rmscsveg(i,j)) / ( fcanc(i,j) + fcancs(i,j))
        !     rmrveg(i,j)= (fcanc(i,j)*rmrcgveg(i,j) + &
        ! &        fcancs(i,j)*rmrcsveg(i,j)) / ( fcanc(i,j) + fcancs(i,j))
        !     rmlveg(i,j)= (fcanc(i,j)*rmlcgveg(i,j) + &
        ! &        fcancs(i,j)*rmlcsveg(i,j)) / ( fcanc(i,j) + fcancs(i,j))
            rmlveg(i,j) = rmlcgveg(i,j)
        !     anveg(i,j)= (fcanc(i,j)*ancgveg(i,j) + &
        ! &        fcancs(i,j)*ancsveg(i,j)) / ( fcanc(i,j) + fcancs(i,j))
            anveg(i,j)= ancgveg(i,j)
            gppveg(i,j)= gppcgveg(i,j)          
        !     gppveg(i,j)= (fcanc(i,j)*gppcgveg(i,j) + &
        ! &        fcancs(i,j)*gppcsveg(i,j)) / ( fcanc(i,j) + fcancs(i,j))
        !     pheanveg(i,j)= (fcanc(i,j)*pancgveg(i,j) + &
        ! &        fcancs(i,j)*pancsveg(i,j)) / ( fcanc(i,j) + fcancs(i,j))
            pheanveg(i,j)= pancgveg(i,j) 
        else
            rmsveg(i,j)= 0.0
            rmrveg(i,j)= 0.0
            rmlveg(i,j)= 0.0
            anveg(i,j)= 0.0
            gppveg(i,j)= 0.0
            pheanveg(i,j)= 0.0
        endif

        if(lfstatus(i,j).eq.4)then
            gppveg(i,j) = anveg(i,j) + rmlveg(i,j)
        endif

        rmveg(i,j)  = rmlveg(i,j) + rmrveg(i,j) + rmsveg(i,j)
        nppveg(i,j) = gppveg(i,j) - rmveg(i,j)

280     continue
270   continue
!>
!>Now that we know maintenance respiration from leaf, stem, and root
!>and gpp, we can find growth respiration for each vegetation
!>
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

!>Calculate grid/tile-averaged rates of rm, rg, npp, and gpp

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
  call    hetresv ( fcancmx, fc, litrmass(:,1:icc), soilcmas(:,1:icc),& ! In
                    delzw, thpor, il1, il2, & ! In
                    ilg,   tbar, psisat, thliq,& ! In
                    roottemp,    zbotw, sort,bi,  & ! In
                    isand, thice, ipeatland, & ! In
                    ltresveg, scresveg) ! Out
!>
!! Find heterotrophic respiration rates from bare ground subarea

call  hetresg  (litrmass(:,iccp1),soilcmas(:,iccp1), delzw,thpor,&
     &               il1,       il2,       ilg,   tbar,    &
     &            psisat,        bi,    thliq,   zbotw,    &
     &            thice,        fg,     isand,    &
     &            ltresveg(:,iccp1),  scresveg(:,iccp1))


  !>  Find heterotrophic respiration rates for the LUC litter and soilc pools
  !!  The LUC litter and soil C respiration rates are assumed to
  !!  be applied over the entire tile but kept in layer 1  

  call  hetresg  (litrmass(:,iccp2),soilcmas(:,iccp2), delzw,thpor,&
       &               il1,       il2,       ilg,   tbar,    &
       &            psisat,        bi,    thliq,   zbotw,    &
       &            thice,        fg,     isand,    &
       &            ltresveg(:,iccp2),  scresveg(:,iccp2))

!>
!! When peatlands are simulated find the peat heterotrophic respiration

call hetres_peat       (    il1,          il2,          ilg,   ipeatland,   &
                          isand,   litrmsmoss,      peatdep,      wtable,   &
                           tbar,        thliq,        thice,       thpor,   &
                             bi,        zbotw,        delzw,      psisat,   &
                     litresmoss,  socres_peat,      resoxic,   resanoxic)

!>
!!Find vegetation averaged litter and soil c respiration rates
!!using values from canopy over ground and canopy over snow subareas
!!
do 340 j = 1, icc
  do 350 i = il1, il2
    if(fcancmx(i,j) > zero) then
        hetrsveg(i,j) =  ltresveg(i,j) + scresveg(i,j)
    else
        hetrsveg(i,j)= 0.0
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
            hetrsveg(i,iccp1) =  ltresveg(i,iccp1) + scresveg(i,iccp1)
            nepveg(i,iccp1)=0.-hetrsveg(i,iccp1)
    else
            hetrsveg(i,iccp1)= 0.0
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
    litres(i)=litres(i)+fcancmx(i,j)*ltresveg(i,j)
    socres(i)=socres(i)+fcancmx(i,j)*scresveg(i,j)
370     continue
360   continue
!>
!!    add moss and peat soil respiration to the grid
!!    In loop 360/370 we have aggregated across all pfts for litres and socres, In loop 380 we add
!!    bareground (iccp1) and LUC pool (iccp2) values to the grid sum if it's not peatland. If it is a peatland, we add litresmoss to
!!    the grid sum but no bareground values as we assume peatlands have no bareground.
!
  nep(:) = 0.
      do 380 i = il1, il2
          if (ipeatland(i) == 0) then
      litres(i) = litres(i) + fg(i) * ltresveg(i,iccp1)
      socres(i) = socres(i) + fg(i) * scresveg(i,iccp1)
          else  ! peatlands
              litres(i) = litres(i)+ litresmoss(i) !add the moss litter, which is assumed to cover whole tile.
!
!              CAUTION says Vivek
!              Note that the following line overwrites socres(i) calculated above in loop 370, although
!              socres(i) is based on scresveg(i,j) which isn't modified. Similarly, the implementation
!              of peatlands also means hetrores(i), calculated below, is now inconsistent with hetrsveg(i,j).
!              The implementation of peatlands is based on overwriting variables without making things
!              consistent. Of course, overwriting variables is not advised because it makes things confusing. 
!
              socres(i) = socres_peat(i) ! since this is only peatland on this tile, use just the peat value.
          endif

              hetrores(i)= litres(i)+socres(i)
              nep(i)=npp(i)-hetrores(i)

380   continue
!>
!!Update the litter and soil c pools based on litter and soil c respiration rates
!!found above. also transfer humidified litter to the soil c pool.
!!
do 420 j = 1, iccp2
  do 430 i = il1, il2

!>Convert u mol co2/m2.sec -> \f$(kg C/m^2)\f$ respired over the model time step
    ltrestep(i,j)=ltresveg(i,j)*(1.0/963.62)*deltat
    screstep(i,j)=scresveg(i,j)*(1.0/963.62)*deltat

!>Update litter and soil c pools
    if (j < iccp1) then
        litrmass(i,j) = litrmass(i,j) - ltrestep(i,j) * (1.0+humicfac(sort(j)))
        hutrstep(i,j) = humicfac(sort(j)) * ltrestep(i,j)
    else
!>         Next we add bareground and LUC pool litter mass and humification for non-peatlands.
        if (ipeatland(i) == 0) then
            litrmass(i,j) = litrmass(i,j) -  ltrestep(i,j)*(1.0+humicfac_bg)
            hutrstep(i,j) = humicfac_bg * ltrestep(i,j)
        !else for peatlands:
        ! In peatlands there is no bareground litter mass since it is the moss layer.
        endif
    endif
!
    humtrsvg(i,j)=hutrstep(i,j)*(963.62/deltat) ! u-mol co2/m2.sec
    soilcmas(i,j)=soilcmas(i,j) + real(spinfast) * (hutrstep(i,j) -  screstep(i,j))

    if(litrmass(i,j).lt.zero) litrmass(i,j)=0.0
    if(soilcmas(i,j).lt.zero) soilcmas(i,j)=0.0
430     continue
420   continue

!>Estimate soil respiration. this is sum of heterotrophic respiration and root maintenance respiration.

do 440 j = 1, icc
  do 450 i = il1, il2
    soilrsvg(i,j)=ltresveg(i,j)+scresveg(i,j)+rmrveg(i,j)
450     continue
440   continue

!> But over the bare fraction and LUC product pool there is no live root.

do 460 i = il1, il2
    soilrsvg(i,iccp1)=ltresveg(i,iccp1)+scresveg(i,iccp1)
    soilrsvg(i,iccp2) = ltresveg(i,iccp2) + scresveg(i,iccp2)
460   continue

!>Find grid averaged humification and soil respiration rates
  soilresp(:) = 0.0
  humiftrs(:) = 0.0
  hutrstep_g(:) = 0.0
do 470 i = il1, il2
  do 480 j = 1,icc
    soilresp(i)=soilresp(i)+fcancmx(i,j)*soilrsvg(i,j)
    humiftrs(i)=humiftrs(i)+fcancmx(i,j)*humtrsvg(i,j)
      hutrstep_g(i) = hutrstep_g(i) + fcancmx(i,j) * hutrstep(i,j) 
480     continue

!
!>    After aggregation of humification and soil respiration rates for non-peatlands
!!   to the grid/tile level, the same must be done for bareground (iccp1).
!!    For peatlands, we additionally add moss values to the grid (litter respiration
!!   and moss root respiration). Note in loop 430 iccp1 is passed for peatlands
!
          if (ipeatland(i) ==0 ) then !non peatland
    
    soilresp(i) = soilresp(i) + fg(i) * soilrsvg(i,iccp1)
    humiftrs(i) = humiftrs(i) + fg(i) * humtrsvg(i,iccp1)
            
          else !peatland
! 
!             CAUTION says Vivek
!             Here again soilresp(i) is overwritten with socres(i)=socres_peat(i) as calculated
!             above in loop 380. This makes soilresp(i) inconsistent with soilrsvg(i,j) for
!             peatland gridcells/tile.      
 
            soilresp(i) = socres(i)+litres(i)+rmr(i) !moss root and litter respiration. No bareground!

  ! Calculate moss time step C fluxes, '/365*deltat' convert year-1
!    to time step-1, 'deltat/963.62' convert umol CO2/m2/s to kgC/m2/deltat.
!    note that hutrstep_g aggregation for icc was done in loop 480
!
              litrfallmoss(i)= Cmossmas(i)*rmortmoss/365*deltat !kgC/m2/day(dt)
              ltrestepmoss(i)= litresmoss(i)*(1.0/963.62)*deltat   !kgC/m2/dt
              nppmosstep(i)= nppmoss(i)*(1.0/963.62)*deltat    !kgC/m2/dt
              socrestep(i) = socres(i)*(1.0/963.62)*deltat     !kgC/m2/dt
              soilresp(i)  = soilresp(i)*(1.0/963.62)*deltat   !kgC/m2/dt
              humstepmoss(i)= humicfacmoss * ltrestepmoss(i)        !kgC/m2/dt
              hutrstep_g(i)= hutrstep_g(i) + humstepmoss(i)     !kgC/m2/dt
              humiftrs(i)  = humiftrs(i)+humstepmoss(i)*(963.62/deltat)!umol/m2/s
          endif
          
470   continue
!     ------------------------------------------------------------------
!
!     heterotrophic respiration part ends
!
!     ------------------------------------------------------------------

!>Find CH4 wetland area (if not prescribed) and emissions:
     call  wetland_methane (hetrores,       il1,       il2,      ilg,  &
     &                      wetfrac,    thliq,   currlat,     sand,  &  
     &                    slopefrac,        ta,                       &
     &                   ch4WetSpec,   wetfdyn, ch4WetDyn)


!> Calculate the methane that is oxidized by the soil sink
    call soil_ch4uptake(        il1,       il2,       ilg,     tbar,  &
     &                           bi,    thliq,     thice,  psisat,  &
     &                       fcanmx,    wetfdyn, wetfrac,  &
     &                        isand,   ch4conc,  &
     &                        thpor,                                  &
! ------------------------- inputs above this line, outputs below ----
     &                    ch4soills)

!    -------------------------------------------------------------------

!>Estimate allocation fractions for leaf, stem, and root components.

       call allocate (lfstatus,   thliq,    ailcg,     ailcb,&
     &                     il1,        il2,       ilg,       sand,  &
     &                    clay,   rmatctem,  gleafmas,   stemmass,  &
     &                rootmass,       sort,    fcancmx,  &
     &                   isand,       THFC,       THLW,             & 
! --------------------- inputs above this line, outputs below ------------
     &                     afrleaf,  afrstem,  afrroot, &
     &                    wtstatus, ltstatus) 


!>
!!
!!Estimate fraction of npp that is to be used for horizontal
!!expansion (lambda) during the next day (i.e. this will be determining
!!the colonization rate in competition).
!!
  if (PFTCompetition)  lambda = expansion(il1,il2,ilg,sort,ailcg,lfstatus,nppveg,pftexist)

!    ------------------------------------------------------------------
!>
!>Maintenance respiration also reduces leaf, stem, and root biomass.
!!when npp for a given pft is positive then this is taken care by
!!allocating +ve npp amongst the leaves, stem, and root component.
!!when npp for a given pft is negative then maintenance respiration
!!loss is explicitly deducted from each component.

      do 600 j = 1, icc
        do 610 i = il1, il2
!>
!!Convert npp and maintenance respiration from different components
  !!from units of u mol co2/m2.sec -> \f$(kg C/m^2)\f$ sequestered or respired over the model time step (deltat)

          gppvgstp(i,j)=gppveg(i,j)*(1.0/963.62)*deltat !+ add2allo(i,j)

!>Remove the cost of making reproductive tissues. This cost can only be removed when NPP is positive.
            reprocost(i,j) =max(0.,nppveg(i,j)*repro_fraction)

!         Not in use. We now use a constant reproductive cost as the prior formulation
!         produces perturbations that do not allow closing of the C balance. JM Jun 2014.
!          nppvgstp(i,j)=nppveg(i,j)*(1.0/963.62)*deltat*(1.-lambda(i,j))
!     &                  + add2allo(i,j)
          nppvgstp(i,j)=(nppveg(i,j)-reprocost(i,j))*(1.0/963.62)*deltat
!
!         Amount of c related to horizontal expansion
!         Not in use. JM Jun 2014
!         expbalvg(i,j)=-1.0*nppveg(i,j)*deltat*lambda(i,j)+ add2allo(i,j)*(963.62/1.0)
!
          rmlvgstp(i,j)=rmlveg(i,j)*(1.0/963.62)*deltat
          rmsvgstp(i,j)=rmsveg(i,j)*(1.0/963.62)*deltat
          rmrvgstp(i,j)=rmrveg(i,j)*(1.0/963.62)*deltat
!
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
!>and since we do not have any real leaves on then we do not take into account co2 uptake by imaginary leaves in carbon budget.
!!rmlvgstp(i,j) should be zero because we set maintenance respiration from storage/imaginary leaves equal to zero. in loop 180
!!
            ntchlveg(i,j)=-rmlvgstp(i,j)
            ntchsveg(i,j)=-rmsvgstp(i,j)
            ntchrveg(i,j)=-rmrvgstp(i,j)
!>
!>since no real leaves are on, make allocation fractions equal to zero.
!>
            afrleaf(i,j)=0.0
            afrstem(i,j)=0.0
            afrroot(i,j)=0.0
          endif
!
          gleafmas(i,j)=gleafmas(i,j)+ntchlveg(i,j)
          stemmass(i,j)=stemmass(i,j)+ntchsveg(i,j)
          rootmass(i,j)=rootmass(i,j)+ntchrveg(i,j)
!
          if(gleafmas(i,j).lt.0.0)then
            write(6,1900)'gleafmas < zero at i=',i,' for pft=',j,''
            write(6,1901)'gleafmas = ',gleafmas(i,j)
            write(6,1901)'ntchlveg = ',ntchlveg(i,j)
            write(6,1902)'lfstatus = ',lfstatus(i,j)
            write(6,1901)'ailcg    = ',ailcg(i,j)
            write(6,1901)'slai     = ',slai(i,j)
1900        format(a23,i4,a10,i2,a1)
1902        format(a11,i4)
            call xit ('ctem',-2)
          endif
!
          if(stemmass(i,j).lt.0.0)then
            write(6,1900)'stemmass < zero at i=(',i,') for pft=',j,')'
            write(6,1901)'stemmass = ',stemmass(i,j)
            write(6,1901)'ntchsveg = ',ntchsveg(i,j)
            write(6,1902)'lfstatus = ',lfstatus(i,j)
            write(6,1901)'rmsvgstp = ',rmsvgstp(i,j)
            write(6,1901)'afrstem  = ',afrstem(i,j)
            write(6,1901)'gppvgstp = ',gppvgstp(i,j)
            write(6,1901)'rmsveg = ',rmsveg(i,j)
1901        format(a11,f12.8)
            call xit ('ctem',-3)
          endif
!
          if(rootmass(i,j).lt.0.0)then
            write(6,1900)'rootmass < zero at i=(',i,') for pft=',j,')'
            write(6,1901)'rootmass = ',rootmass(i,j)
            call xit ('ctem',-4)
          endif
!>
!!convert net change in leaf, stem, and root biomass into
!!u-mol co2/m2.sec for use in balcar subroutine
!!
          ntchlveg(i,j)=ntchlveg(i,j)*(963.62/deltat)
          ntchsveg(i,j)=ntchsveg(i,j)*(963.62/deltat)
          ntchrveg(i,j)=ntchrveg(i,j)*(963.62/deltat)
!>
!!to avoid over/underflow problems set gleafmas, stemmass, and
!!rootmass to zero if they get too small
!!
          if(bleafmas(i,j).lt.zero) bleafmas(i,j)=0.0
          if(gleafmas(i,j).lt.zero) gleafmas(i,j)=0.0
          if(stemmass(i,j).lt.zero) stemmass(i,j)=0.0
          if(rootmass(i,j).lt.zero) rootmass(i,j)=0.0
!
610     continue
600   continue
!>
!>calculate grid averaged value of C related to spatial expansion
!>
  repro_cost_g(:)=0.0    !<amount of C for production of reproductive tissues
      do 620 j = 1,icc
        do 621 i = il1, il2
         !if (PFTCompetition .or. lnduseon) then
!           Not in use. We now use the constant reproductive cost below. JM Jun 2014
!           expnbaln(i)=expnbaln(i)+fcancmx(i,j)*expbalvg(i,j)
            repro_cost_g(i)=repro_cost_g(i)+fcancmx(i,j)*reprocost(i,j)
         !endif
621     continue
620   continue
!
!    ------------------------------------------------------------------

!>Phenology part starts
!!
!!the phenology subroutine determines leaf status for each pft and calculates leaf litter.
!!the phenology subroutine uses soil temperature (tbar) and root temperature. however,
!!since ctem doesn't make the distinction between canopy over ground, and canopy over
!!snow sub-areas for phenology purposes (for  example, leaf onset is not assumed to occur
!!at different times over these sub-areas) we use average soil and root temperature in the phenology subroutine.
!!
!>
!!Call the phenology subroutine, which determines the leaf growth
!!status, calculates leaf litter, and converts green grass into brown.
!!

call phenolgy(gleafmas,   bleafmas,        il1,        il2,    &
     &             ilg,    leapnow,    tbar,     thice,    &
     &                      thliq,     THLW,     THFC,       ta,&
     &                    pheanveg,     iday,     radj, roottemp,&
     &                    rmatctem, stemmass, rootmass,     sort,&
     &                    fcancmx,  isand, &
!     ------------------ inputs above this line ----------------------
     &                    flhrloss, leaflitr, lfstatus,  pandays,&
     &                    colddays)
!     --- variables which are updated and outputs above this line ----
!

  !> While leaf litter is calculated in the phenology subroutine, stem
  !! and root turnover is calculated in the turnoverStemRoot subroutine.
!!
call turnoverStemRoot (stemmass, rootmass,  lfstatus,    ailcg,&
     &              il1,       il2,       ilg,  leapnow,  &
     &                         sort, fcancmx,&
!  ------------------ inputs above this line ----------------------
     &                     stmhrlos, rothrlos,&
! ----------- inputs which are updated above this line -----------
     &                     stemlitr, rootlitr)
! ------------------outputs above this line ----------------------


  !> Update green leaf biomass for trees and crops, brown leaf biomass for grasses,
  !! stem and root biomass for litter deductions, and update litter pool with leaf
  !! litter calculated in the phenology subroutine and stem and root litter 
  !! calculated in the turnoverStemRoot subroutine. Also add the reproduction
  !!  carbon directly to the litter pool
  !!
  call updatePoolsTurnover(il1, il2, reprocost, & !In
                          stemmass, rootmass, litrmass, rootlitr,& !In/Out
                         gleafmas, bleafmas, leaflitr, stemlitr) !In/Out

!    ------------------------------------------------------------------
!>
!>Call the mortaliy subroutine which calculates mortality due to reduced growth and aging. exogenous mortality due to fire and other
!!disturbances and the subsequent litter that is generated is calculated in the disturb subroutine.
!!
!!set maxage >0 in classic_params.f90 to switch on mortality due to age and
!!reduced growth. Mortality is linked to the competition parameterization and generates bare fraction.
!!

call       mortalty (stemmass,   rootmass,    ailcg,   gleafmas,  &
     &               bleafmas,        il1,      il2,        ilg,  &
     &                leapnow,       iday,     sort,    fcancmx,  &
! ----------------- inputs above this line ----------------------
     &               lystmmas,   lyrotmas, tymaxlai,   grwtheff,  &
! -------------- inputs updated above this line ------------------
     &               stemltrm,   rootltrm, glealtrm,   geremort,  &
     &               intrmort)
! ------------------outputs above this line ----------------------

!>Update leaf, stem, and root biomass pools to take into loss due to mortality, and put the
!!litter into the litter pool. the mortality for green grasses doesn't generate litter, instead they turn brown.
!!
!  call updatePoolsMortality(il1, il2, stemltrm, rootltrm, & !In
!                            stemmass, rootmass, litrmass, & !In/Out
!                            glealtrm, gleafmas, bleafmas) !In/Out
! to:
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
          stemmass(i,m)=stemmass(i,m)-stemltrm(i,m)
          rootmass(i,m)=rootmass(i,m)-rootltrm(i,m)
          litrmass(i,m)=litrmass(i,m)+stemltrm(i,m)+rootltrm(i,m)  
          select case(classpfts(j))
            case ('NdlTr' , 'BdlTr', 'Crops', 'BdlSh')
              gleafmas(i,m)=gleafmas(i,m)-glealtrm(i,m)
            case('Grass')    ! grasses
            gleafmas(i,m)=gleafmas(i,m)-glealtrm(i,m)
            bleafmas(i,m)=bleafmas(i,m)+glealtrm(i,m)
            glealtrm(i,m)=0.0
            case default
              print*,'Unknown CLASS PFT in ctem ',classpfts(j)
              call XIT('ctem',-6)                                                                       
          end select
          litrmass(i,m)=litrmass(i,m)+glealtrm(i,m)
840     continue
835    continue
830   continue
! here

!    ------------------------------------------------------------------
!>
!>call the disturbance subroutine which calculates mortality due to fire and other disturbances.
!>the primary output from from disturbance subroutine is litter generated, c emissions due to fire
!!and area burned, which may be used to estimate change in fractional coverages.
!!
!!disturbance is spatial and requires area of gcm grid cell and areas of different pfts present in
!!a given grid cell. however, when ctem is operated at a point scale then it is assumed that the
!!spatial scale is 1 hectare = 10,000 m2. the disturbance subroutine may be stopped from simulating
!!any fire by specifying fire extingushing probability equal to 1.
!!
call disturb (stemmass, rootmass, gleafmas, bleafmas,&
     &                      thliq,    THLW,      THFC,    uwind,&
     &                       vwind,  lightng,  fcancmx, litrmass,&
     &                    rmatctem,     ilg,           &
     &                         il1,      il2,     sort, &
     &                    grclarea,   thice,   popdin, lucemcom,&
     &                      dofire,  currlat,     iday, fsnow,&
     &                       isand,  &

!    ------------------- inputs above this line ---------------------
!    ------------ outputs below are the primary outputs -------------

     &                    stemltdt, rootltdt, glfltrdt, blfltrdt,&
     &                    glcaemls, rtcaemls, stcaemls,&
     &                    blcaemls,   ltrcemls,   burnfrac,               &
     &                    pstemmass,  pgleafmass, emit_co2,    emit_ch4,   &

!    ------------ outputs below are the secondary outputs -------------
!                         which may be omitted in AGCM

     &                    emit_co,  emit_nmhc,                           &
     &                    emit_h2,  emit_nox, emit_n2o, emit_pm25,&
     &                    emit_tpm, emit_tc,  emit_oc,  emit_bc,&
     &                    burnvegf, bterm_veg,mterm_veg,  lterm,&
     &                   smfunc_veg                                       & 
     &                   )

!    ------------------------------------------------------------------
!>
!> Calculate NBP (net biome production) for each pft by taking into account
!! C emission losses. The disturbance routine produces emissions due to fire
!! and while the land use change subroutine calculates emissions due to LUC. 
!! The LUC related combustion flux is assumed to be spread uniformly over the
!! tile as it is no longer associated with any one PFT. To calculate the NBP 
!! we do not subtract LUC emissions from the PFT-level NBP but we do subtract 
!! it from the per tile NBP. 
!!
!  call calcNBP(il1, il2, ilg,  deltat, nepveg, fcancmx, & !In
!                    lucemcom, ltresveg, scresveg, nep, & !In
!                    glcaemls, blcaemls, stcaemls, rtcaemls, ltrcemls, & ! In/Out
!                    nbpveg, dstcemls1, dstcemls3, nbp) ! Out 
!to:
      do 1000 i = il1, il2
        do 1010 j = 1, icc
          dscemlv1(i,j) = glcaemls(i,j) + blcaemls(i,j) + stcaemls(i,j) + rtcaemls(i,j)
          dscemlv2(i,j) = dscemlv1(i,j) + ltrcemls(i,j)

!         convert \f$kg c/m^2\f$ emitted in one day into u mol co2/m2.sec before
!         subtracting emission losses from nep.
          nbpveg(i,j)  =nepveg(i,j) - dscemlv2(i,j)*(963.62/deltat)

1010    continue

!       For accounting purposes, we also need to account for the bare fraction
!       NBP. Since there is no fire on the bare, we use 0.
        nbpveg(i,iccp1)  =nepveg(i,iccp1)   - 0.

1000  continue
!>
!! Calculate grid. averaged rate of carbon emissions due to fire in u-mol co2/m2.sec. 
!! Convert all emission losses from \f$kg c/m^2\f$ emitted in 1 day to u-mol co2/m2.sec.
!! Calculate grid averaged carbon emission losses from litter.
!!
      do 1030 j = 1,icc
        do 1040 i = il1, il2
          dstcemls1(i) = dstcemls1(i) + fcancmx(i,j) * dscemlv1(i,j) * (963.62 / deltat)
          dstcemls2(i) = dstcemls2(i) + fcancmx(i,j) * dscemlv2(i,j) * (963.62 / deltat)
          galtcels(i)  = galtcels(i) + fcancmx(i,j) * ltrcemls(i,j) * (963.62 / deltat)
          glcaemls(i,j)= glcaemls(i,j)*(963.62/deltat)
          blcaemls(i,j)=blcaemls(i,j)*(963.62/deltat)
          stcaemls(i,j)=stcaemls(i,j)*(963.62/deltat)
          rtcaemls(i,j)=rtcaemls(i,j)*(963.62/deltat)
          ltrcemls(i,j)=ltrcemls(i,j)*(963.62/deltat)
1040    continue
1030  continue
!
!       For the tile-level NBP, we include the disturbance emissions as well as
!       respiration from the paper (litter) and furniture (soil carbon) pools (LUC
!       product pools). Also include here the instantaneous emissions due to LUC.
      do 1041 i = il1, il2
        nbp(i) = nep(i) - dstcemls2(i) - (ltresveg(i,iccp2) + scresveg(i,iccp2)) - lucemcom(i)
        dstcemls3(i) = dstcemls2(i) - dstcemls1(i)  !litter is total - vegetation.
1041  continue

! here.

!> Prepare for the carbon balance check. Calculate total litter fall from each 
  !! component (leaves, stem, and root) from all causes (normal turnover, drought
  !! and cold stress for leaves, mortality, and disturbance), calculate grid-average
  !! vegetation biomass, litter mass, and soil carbon mass, and litter fall rate.
  !! Also add the bare ground values to the grid-average. If a peatland, we assume no bareground and
  !! add the moss values instead. Note: peatland soil C is not aggregated from plants but updated
  !! by humification and respiration from the previous stored value
  !!
  !call prepBalanceC(il1, il2, ilg, fcancmx, glealtrm, glfltrdt, &  ! In 
  !                       blfltrdt, stemltrm, stemltdt, rootltrm, rootltdt, &
  !                       ipeatland, nppmosstep, pgavscms, humstepmoss, &
  !                       ltrestepmoss, stemlitr, rootlitr, rootmass, &
  !                       litrmass, soilCmas, hutrstep_g, stemmass, bleafmas, &
  !                       gleafmas, socrestep, fg, litrfallmoss, &
  !                       leaflitr, Cmossmas, litrmsmoss,  & ! In/Out
  !                       tltrleaf, tltrstem, tltrroot, & ! Out 
  !                       vgbiomas, litrfall, gavgltms, litrfallveg, &
  !                       gavgscms, vgbiomas_veg)
   ! to:            

!>
!! Calculate total litter fall from each component (leaves, stem, and root)
!!  from all causes (normal turnover, drought and cold stress for leaves, mortality, 
!! and disturbance) for use in balcar subroutine
!!
      do 1050 j = 1,icc
        do 1060 i = il1, il2
          
        ! units here are \f$kg c/m^2 .day\f$
         tltrleaf(i,j)=leaflitr(i,j)+glealtrm(i,j)+glfltrdt(i,j)+&
     &                 blfltrdt(i,j)
         tltrstem(i,j)=stemlitr(i,j)+stemltrm(i,j)+stemltdt(i,j)
         tltrroot(i,j)=rootlitr(i,j)+rootltrm(i,j)+rootltdt(i,j)
!>
!>convert units to u-mol co2/m2.sec
         leaflitr(i,j)=leaflitr(i,j)*(963.62/deltat)
         tltrleaf(i,j)=tltrleaf(i,j)*(963.62/deltat)
         tltrstem(i,j)=tltrstem(i,j)*(963.62/deltat)
         tltrroot(i,j)=tltrroot(i,j)*(963.62/deltat)
1060    continue
1050  continue
!>
!>calculate grid-average vegetation biomass, litter mass, and soil carbon mass, and litter fall rate
!>
      do 1100 j = 1, icc
        do 1110 i = il1, il2
          vgbiomas(i)=vgbiomas(i)+fcancmx(i,j)*(gleafmas(i,j)+&
     &     bleafmas(i,j)+stemmass(i,j)+rootmass(i,j))
          litrfall(i)=litrfall(i)+fcancmx(i,j)*(tltrleaf(i,j)+&
     &     tltrstem(i,j)+tltrroot(i,j))
          ! store the per PFT litterfall for outputting.
          litrfallveg(i,j)=(tltrleaf(i,j)+tltrstem(i,j)+tltrroot(i,j))
          gavgltms(i)=gavgltms(i)+fcancmx(i,j)*litrmass(i,j)

          if (ipeatland(i)==0) then ! Non-peatlands
               gavgscms(i)=gavgscms(i)+fcancmx(i,j)*soilcmas(i,j)
          !else
             !Peatland soil C is calculated from peat depth (peatdep) in the peatland
          endif
          vgbiomas_veg(i,j)=gleafmas(i,j)+&
     &     bleafmas(i,j)+stemmass(i,j)+rootmass(i,j) !vegetation biomass for each pft
1110    continue
1100  continue
!
!>    Add the bare ground values to the grid-average. If a peatland, we assume no bareground and
!!    add the moss values instead.
!!    Note: peatland soil C is not aggregated from plants but updated
!!    by humification and respiration from the previous stored value
!
      do 1020 i = il1, il2
          if (ipeatland(i)==0) then
          ! Add the bare fraction dead C
             ! gavgltms(i)=gavgltms(i)+( (fg(i)+fgs(i))*litrmass(i,iccp1))
             ! gavgscms(i)=gavgscms(i)+( (fg(i)+fgs(i))*soilcmas(i,iccp1))
             gavgltms(i)=gavgltms(i)+(fg(i)*litrmass(i,iccp1))
             gavgscms(i)=gavgscms(i)+(fg(i)*soilcmas(i,iccp1))
             
          else
             litrmsmoss(i)= litrmsmoss(i)+litrfallmoss(i)-&
      &                     ltrestepmoss(i)-humstepmoss(i)     !kg/m2
             Cmossmas(i)= Cmossmas(i)+nppmosstep(i)-litrfallmoss(i)
             vgbiomas(i) = vgbiomas(i) + Cmossmas(i)
             litrfall(i) = litrfall(i) + litrfallmoss(i)*(963.62/deltat)!umolCO2/m2/s
             gavgltms(i) = gavgltms(i) + litrmsmoss(i)
             gavgscms(i) = pgavscms(i) + hutrstep_g(i)- socrestep(i)
             ! Calculate the peat depth based on equation 18 in Wu, Verseghy, Melton 2016 GMD.
             peatdep(i)=(-72067.0+sqrt((72067.0**2.0)-(4.0*4056.6*&
     &              (-gavgscms(i)*1000/0.487))))/(2*4056.6)
          endif
1020  continue
! here.
!     -----------------------------------------------------------------
!>
!>At this stage we have all required fluxes in u-mol co2/m2.sec and initial (loop 140 and 145)
  !!and updated sizes of all pools (in \f$(kg C/m^2)\f$). Now we call the balcar subroutine and make sure
!!that C in leaves, stem, root, litter and soil C pool balances within a certain tolerance.
!!
if(spinfast.eq.1)then
        call  balcar(gleafmas, stemmass, rootmass,  bleafmas,&
&                    litrmass, soilcmas, ntchlveg,  ntchsveg,&
&                    ntchrveg, tltrleaf, tltrstem,  tltrroot,&
&                    glcaemls, blcaemls, stcaemls,  rtcaemls,&
&                    ltrcemls, ltresveg, scresveg,  humtrsvg,&
&                    pglfmass, pblfmass, pstemass,  protmass,&
&                    plitmass, psocmass, vgbiomas,  reprocost,&
&                    pvgbioms, gavgltms, pgavltms,  gavgscms,&
  &                    pgavscms, dstcemls3, repro_cost_g,       &
&                     autores, hetrores,      gpp,    &
&                      litres,   socres, dstcemls1,    &
&                    litrfall, humiftrs,                 &
&                         il1,      il2,      ilg,     &
&                   ipeatland, Cmossmas, pCmossmas,                &
&                  nppmosstep, litrfallmoss, litrmsmoss,&
&                 plitrmsmoss, ltrestepmoss, humstepmoss)

endif

!     -----------------------------------------------------------------
!>
!>Finally find vegetation structural attributes which can be passed to the land surface scheme using leaf, stem, and root biomass.
!>
call bio2str( gleafmas, bleafmas, stemmass, rootmass,&
     &              il1,        il2,        ilg,      zbotw,  &
     &            delzw,   soildpth,    fcancmx,  &
     &        ipeatland,                                      &
! --------------- inputs above this line, outputs below --------
     &                          ailcg,    ailcb,     ailc,    zolnc,&
     &                          rmatc, rmatctem,     slai,  bmasveg,&
     &                       cmasvegc,  veghght, rootdpth,   alvisc,&
     &           alnirc,    paicgat,   slaicgat)

!>
!>Calculation of gavglai is moved from loop 1100 to here since ailcg is updated by bio2str
!>
do j = 1, icc
    do i = il1, il2
        gavglai (i)=gavglai (i)+fcancmx(i,j)*ailcg(i,j)
    enddo
enddo

!>
  !> At the end of the day, find the depth of the peat, update the degree days for moss photosynthesis and the peat bottom layer depth
!>
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
  real, intent(in) :: ltresveg(:,:)      !< Litter respiration for each pft, bare fraction, and LUC product pool, \f$\mu mol CO_2 m^{-2} s^{-1}\f$
  real, intent(in) :: scresveg(:,:)      !< Soil carbon respiration for each pft, bare fraction, and LUC product pool, \f$\mu mol CO_2 m^{-2} s^{-1}\f$
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

!     Convert \f$(kg C/m^2)\f$ emitted in one day into u mol co2/m2.sec before
!     subtracting emission losses from nep.
      nbpveg(i,j) = nepveg(i,j) - dscemlv2(i,j)*(963.62/deltat)

101    continue

!   For accounting purposes, we also need to account for the bare fraction
!   NBP. Since there is no fire on the bare, we use 0.
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
!
!       For the tile-level NBP, we include the disturbance emissions as well as
!       respiration from the paper (litter) and furniture (soil carbon) pools (LUC
!       product pools). Also include here the instantaneous emissions due to LUC.
  nbp(:) = 0.
  dstcemls3 = 0.0
  do 105 i = il1, il2
    nbp(i) = nep(i) - dstcemls2(i) - (ltresveg(i,iccp2) + scresveg(i,iccp2)) - lucemcom(i)
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
