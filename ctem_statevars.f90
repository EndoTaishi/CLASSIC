!>\defgroup ctem_statevars

!>
!!this module contains the variable type structures:
!! 1. c_switch - switches for running CTEM, read from the joboptions file
!! 2. vrot - CTEM's 'rot' vars
!! 3. vgat - CTEM's 'gat' vars
!! 4. class_out - CLASS's monthly outputs
!! 5. ctem_grd - CTEM's grid average variables
!! 6. ctem_tile - CTEM's variables per tile
!! 7. ctem_mo - CTEM's variables monthly averaged (per pft)
!! 8. ctem_grd_mo - CTEM's grid average monthly values
!! 9. ctem_tile_mo - CTEM's variables per tile monthly values
!! 10. ctem_yr - CTEM's average annual values (per PFT)
!! 11. ctem_grd_yr - CTEM's grid average annual values
!! 12. ctem_tile_yr - CTEM's variables per tile annual values
!
!>\file
module ctem_statevars

! J. Melton Apr 2015

use ctem_params,  only : initpftpars, nlat, nmos, ilg, nmon,ican, ignd,icp1, icc, iccp1, &
                    monthend, mmday,modelpft, l2max,deltat, abszero, monthdays,seed, crop, NBS

implicit none

public :: alloc_ctem_vars
public :: initrowvars
public :: resetclassmon
public :: resetclassyr
public :: resetdaily
public :: resetmonthend
public :: resetyearend
public :: resetclassaccum
public :: resetgridavg
public :: finddaylength

!=================================================================================
!>switches for running the model, read from the joboptions file

type ctem_switches

    logical :: ctem_on     !<
    logical :: parallelrun !<set this to be true if model is run in parallel mode for 
                           !<multiple grid cells, output is limited to monthly & yearly 
                    	   !<grid-mean only. else the run is in stand alone mode, in which 
                     	   !<output includes half-hourly and daily and mosaic-mean as well.
    logical :: cyclemet    !<to cycle over only a fixed number of years 
                 	   !<(nummetcylyrs) starting at a certain year (metcylyrst)
                 	   !<if cyclemet, then put co2on = false and set an appopriate setco2conc, also
                 	   !<if popdon is true, it will choose the popn and luc data for year
                 	   !<metcylyrst and cycle on that.
    logical :: dofire      !<boolean, if true allow fire, if false no fire.
    logical :: run_model   !<
    logical :: met_rewound !<
    logical :: reach_eof   !<
    logical :: compete     !<logical boolean telling if competition between pfts is on or not
    logical :: start_bare  !<set this to true if competition is true, and if you wish to start from bare ground.
                           !<if this is set to false, the ini and ctm file info will be used to set up the run.
                           !<NOTE: This still keeps the crop fractions (while setting all pools to zero)
    logical :: rsfile      !<set this to true if restart files (.ini_rs and .ctm_rs) are written at the end of each
                           !<year. these files are necessary for checking whether the model reaches equilibrium after
                           !<running for a certain years. set this to false if restart files are not needed
                           !<(known how many years the model will run)
    logical :: lnduseon    !<logical switch to run the land use change subroutine or not.
    logical :: co2on       !<use \f$co_2\f$ time series, set to false if cyclemet is true
    logical :: ch4on       !<use \f$CH_4\f$ time series, set to false if cyclemet is true the \f$CO_2\f$ timeseries is in the 
                           !<same input file as the \f$CO_2\f$ one.
    logical :: popdon      !<if set true use population density data to calculate fire extinguishing probability and
                           !<probability of fire due to human causes, or if false, read directly from .ctm file
    logical :: inibioclim  !<switch telling if bioclimatic parameters are being initialized
                           !<from scratch (false) or being initialized from some spun up
                           !<values(true).
    logical :: start_from_rs!<if true, this option copies the _RS INI and CTM files to be the .INI and .CTM files and
                           !<then starts the run as per normal. it is handy when spinning up so you don't have to do a
                           !<complicated copying of the RS files to restart from them. NOTE! This will not work on
                           !<hadar or spica, instead you have to manually move the files and set this to .false.    
    logical :: dowetlands   !<if true allow wetland methane emission
    logical :: obswetf      !<observed wetland fraction
    logical :: transient_run!<
    logical :: use_netcdf        !< If true all model inputs and outputs are handled via netcdfs
    character(180) :: met_file   !< location of the netcdf meteorological dataset
    character(180) :: init_file  !< location of the netcdf initialization file

    ! CLASS switches:
!     
!     integer :: idisp    !< if idisp=0, vegetation displacement heights are ignored,
!                                  !< because the atmospheric model considers these to be part
!                                  !< of the "terrain".
!                                  !< if idisp=1, vegetation displacement heights are calculated.
! 
!     integer :: izref    !< if izref=1, the bottom of the atmospheric model is taken
!                                  !< to lie at the ground surface.
!                                  !< if izref=2, the bottom of the atmospheric model is taken
!                                  !< to lie at the local roughness height.
! 
!     integer :: islfd    !< if islfd=0, drcoef is called for surface stability corrections
!                                  !< and the original gcm set of screen-level diagnostic calculations 
!                                  !< is done.
!                                  !< if islfd=1, drcoef is called for surface stability corrections
!                                  !< and sldiag is called for screen-level diagnostic calculations. 
!                                  !< if islfd=2, flxsurfz is called for surface stability corrections
!                                  !< and diasurf is called for screen-level diagnostic calculations. 
! 
!     integer :: ipcp     !< if ipcp=1, the rainfall-snowfall cutoff is taken to lie at 0 c.
!                                  !< if ipcp=2, a linear partitioning of precipitation betweeen 
!                                  !< rainfall and snowfall is done between 0 c and 2 c.
!                                  !< if ipcp=3, rainfall and snowfall are partitioned according to
!                                  !< a polynomial curve between 0 c and 6 c.
! 
!     integer :: iwf     !< if iwf=0, only overland flow and baseflow are modelled, and
!                                 !< the ground surface slope is not modelled.
!                                 !< if iwf=n (0<n<4), the watflood calculations of overland flow 
!                                 !< and interflow are performed; interflow is drawn from the top 
!                                 !< n soil layers.
! 
!     integer :: ITC!< itc, itcg and itg are switches to choose the iteration scheme to
!                            !< be used in calculating the canopy or ground surface temperature
!                            !< respectively.  if the switch is set to 1, a bisection method is
!                            !< used; if to 2, the newton-raphson method is used.
!     integer :: ITCG!< itc, itcg and itg are switches to choose the iteration scheme to
!                            !< be used in calculating the canopy or ground surface temperature
!                            !< respectively.  if the switch is set to 1, a bisection method is
!                            !< used; if to 2, the newton-raphson method is used.
!     integer :: ITG!< itc, itcg and itg are switches to choose the iteration scheme to
!                            !< be used in calculating the canopy or ground surface temperature
!                            !< respectively.  if the switch is set to 1, a bisection method is
!                            !< used; if to 2, the newton-raphson method is used.
!    
!     integer :: IPAI !< if ipai, ihgt, ialc, ials and ialg are zero, the values of 
!                            !< plant area index, vegetation height, canopy albedo, snow albedo
!                            !< and soil albedo respectively calculated by class are used.
!                            !< if any of these switches is set to 1, the value of the
!                            !< corresponding parameter calculated by class is overridden by
!                            !< a user-supplied input value.  
!     integer :: IHGT !< if ipai, ihgt, ialc, ials and ialg are zero, the values of 
!                            !< plant area index, vegetation height, canopy albedo, snow albedo
!                            !< and soil albedo respectively calculated by class are used.
!                            !< if any of these switches is set to 1, the value of the
!                            !< corresponding parameter calculated by class is overridden by
!                            !< a user-supplied input value. 
!     integer :: IALC !< if ipai, ihgt, ialc, ials and ialg are zero, the values of 
!                            !< plant area index, vegetation height, canopy albedo, snow albedo
!                            !< and soil albedo respectively calculated by class are used.
!                            !< if any of these switches is set to 1, the value of the
!                            !< corresponding parameter calculated by class is overridden by
!                            !< a user-supplied input value. 
!     integer :: IALS !< if ipai, ihgt, ialc, ials and ialg are zero, the values of 
!                            !< plant area index, vegetation height, canopy albedo, snow albedo
!                            !< and soil albedo respectively calculated by class are used.
!                            !< if any of these switches is set to 1, the value of the
!                            !< corresponding parameter calculated by class is overridden by
!                            !< a user-supplied input value. 
!     integer :: IALG !< if ipai, ihgt, ialc, ials and ialg are zero, the values of 
!                            !< plant area index, vegetation height, canopy albedo, snow albedo
!                            !< and soil albedo respectively calculated by class are used.
!                            !< if any of these switches is set to 1, the value of the
!                            !< corresponding parameter calculated by class is overridden by
!                            !< a user-supplied input value. 
! 
!     ! -------------
!     ! >>>> note: if you wish to use the values in the .ini file, set all to -9999 in the job options file
!     !            and the .ini file will be used.
! 
!     integer :: jhhstd  !< day of the year to start writing the half-hourly output
!     integer :: jhhendd !< day of the year to stop writing the half-hourly output
!     integer :: jdstd   !< day of the year to start writing the daily output
!     integer :: jdendd  !< day of the year to stop writing the daily output
!     integer :: jhhsty  !< simulation year (iyear) to start writing the half-hourly output
!     integer :: jhhendy !< simulation year (iyear) to stop writing the half-hourly output
!     integer :: jdsty   !< simulation year (iyear) to start writing the daily output
!     integer :: jdendy  !< simulation year (iyear) to stop writing the daily output
! 
! 
!     integer :: isnoalb !< if isnoalb is set to 0, the original two-band snow albedo algorithms are used.  
!                                 !< if it is set to 1, the new four-band routines are used.
                                
    character(80) :: titlec1!<
    character(80) :: titlec2!<
    character(80) :: titlec3!<

end type ctem_switches

type (ctem_switches), save, target :: c_switch

!=================================================================================
!>CTEM's 'rot' vars
type veg_rot

    logical, allocatable, dimension(:,:,:) :: pftexist  !<logical array indicating pfts exist (t) or not (f)
    integer, allocatable, dimension(:,:,:) :: lfstatus  !<leaf phenology status
    integer, allocatable, dimension(:,:,:) :: pandays   !<days with positive net photosynthesis (an) for use in
                                                   !<the phenology subroutine
    real, allocatable, dimension(:,:,:) :: ailcmin      !<
    real, allocatable, dimension(:,:,:) :: ailcmax      !<
    real, allocatable, dimension(:,:,:) :: dvdfcan      !<
    real, allocatable, dimension(:,:,:) :: gleafmas     !<green leaf mass for each of the 9 ctem pfts, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:,:) :: bleafmas     !<brown leaf mass for each of the 9 ctem pfts, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:,:) :: stemmass     !<stem mass for each of the 9 ctem pfts, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:,:) :: rootmass     !<root mass for each of the 9 ctem pfts, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:,:) :: pstemmass    !<stem mass from previous timestep, is value before fire. used by burntobare subroutine
    real, allocatable, dimension(:,:,:) :: pgleafmass   !<root mass from previous timestep, is value before fire. used by burntobare subroutine
    real, allocatable, dimension(:,:,:) :: fcancmx      !<max. fractional coverage of ctem's 9 pfts, but this can be
                                                   !<modified by land-use change, and competition between pfts
    real, allocatable, dimension(:,:,:) :: ailcg        !<green lai for ctem's 9 pfts
    real, allocatable, dimension(:,:,:) :: ailcgs       !<GREEN LAI FOR CANOPY OVER SNOW SUB-AREA
    real, allocatable, dimension(:,:,:) :: fcancs       !<FRACTION OF CANOPY OVER SNOW FOR CTEM's 9 PFTs
    real, allocatable, dimension(:,:,:) :: fcanc        !<FRACTIONAL COVERAGE OF 8 CARBON PFTs, CANOPY OVER SNOW
    real, allocatable, dimension(:,:,:) :: co2i1cg      !<INTERCELLULAR CO2 CONC FOR 8 PFTs FOR CANOPY OVER GROUND SUBAREA (Pa) - FOR SINGLE/SUNLIT LEAF
    real, allocatable, dimension(:,:,:) :: co2i1cs      !<SAME AS ABOVE BUT FOR SHADED LEAF (above being co2i1cg)
    real, allocatable, dimension(:,:,:) :: co2i2cg      !<INTERCELLULAR CO2 CONC FOR 8 PFTs FOR CANOPY OVER SNOWSUBAREA (Pa) - FOR SINGLE/SUNLIT LEAF
    real, allocatable, dimension(:,:,:) :: co2i2cs      !<SAME AS ABOVE BUT FOR SHADED LEAF (above being co2i2cg)
    real, allocatable, dimension(:,:,:) :: ancsveg      !<net photosynthetic rate for ctems 9 pfts for canopy over snow subarea
    real, allocatable, dimension(:,:,:) :: ancgveg      !<net photosynthetic rate for ctems 9 pfts for canopy over ground subarea
    real, allocatable, dimension(:,:,:) :: rmlcsveg     !<leaf respiration rate for ctems 9 pfts forcanopy over snow subarea
    real, allocatable, dimension(:,:,:) :: rmlcgveg     !<leaf respiration rate for ctems 9 pfts forcanopy over ground subarea
    real, allocatable, dimension(:,:,:) :: slai         !<storage/imaginary lai for phenology purposes
    real, allocatable, dimension(:,:,:) :: ailcb        !<brown lai for ctem's 9 pfts. for now we assume only grasses can have brown lai
    real, allocatable, dimension(:,:,:) :: flhrloss     !<fall or harvest loss for deciduous trees and crops, respectively, \f$kg c/m^2\f$il1
    real, allocatable, dimension(:,:,:) :: grwtheff     !<growth efficiency. change in biomass per year per unit max.
                                                   !<lai (\f$kg c/m^2\f$)/(m2/m2), for use in mortality subroutine
    real, allocatable, dimension(:,:,:) :: lystmmas     !<stem mass at the end of last year
    real, allocatable, dimension(:,:,:) :: lyrotmas     !<root mass at the end of last year
    real, allocatable, dimension(:,:,:) :: tymaxlai     !<this year's maximum lai
    real, allocatable, dimension(:,:,:) :: stmhrlos     !<stem harvest loss for crops, \f$kg c/m^2\f$                                                   
    real, allocatable, dimension(:,:,:) :: vgbiomas_veg !<vegetation biomass for each pft
    real, allocatable, dimension(:,:,:) :: emit_co2     !<carbon dioxide
    real, allocatable, dimension(:,:,:) :: emit_co      !<carbon monoxide
    real, allocatable, dimension(:,:,:) :: emit_ch4     !<methane
    real, allocatable, dimension(:,:,:) :: emit_nmhc    !<non-methane hydrocarbons
    real, allocatable, dimension(:,:,:) :: emit_h2      !<hydrogen gas
    real, allocatable, dimension(:,:,:) :: emit_nox     !<nitrogen oxides
    real, allocatable, dimension(:,:,:) :: emit_n2o     !<nitrous oxide
    real, allocatable, dimension(:,:,:) :: emit_pm25    !<particulate matter less than 2.5 um in diameter
    real, allocatable, dimension(:,:,:) :: emit_tpm     !<total particulate matter
    real, allocatable, dimension(:,:,:) :: emit_tc      !<total carbon
    real, allocatable, dimension(:,:,:) :: emit_oc      !<organic carbon
    real, allocatable, dimension(:,:,:) :: emit_bc      !<black carbon
    real, allocatable, dimension(:,:,:) :: burnvegf     !<per PFT fraction burned of that PFT's area
    real, allocatable, dimension(:,:,:) :: smfuncveg    !<
    real, allocatable, dimension(:,:,:) :: bterm        !<biomass term for fire probabilty calc                                
    real, allocatable, dimension(:,:,:) :: mterm        !<moisture term for fire probabilty calc
    real, allocatable, dimension(:,:,:) :: bmasveg      !<total (gleaf + stem + root) biomass for each ctem pft, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:,:) :: veghght      !<vegetation height (meters)
    real, allocatable, dimension(:,:,:) :: rootdpth     !<99% soil rooting depth (meters)
                                                   !<both veghght & rootdpth can be used as diagnostics to see
                                                   !<how vegetation grows above and below ground, respectively
    real, allocatable, dimension(:,:,:) :: tltrleaf     !<total leaf litter fall rate (u-mol co2/m2.sec)
    real, allocatable, dimension(:,:,:) :: tltrstem     !<total stem litter fall rate (u-mol co2/m2.sec)
    real, allocatable, dimension(:,:,:) :: tltrroot     !<total root litter fall rate (u-mol co2/m2.sec)
    real, allocatable, dimension(:,:,:) :: leaflitr     !<leaf litter fall rate (u-mol co2/m2.sec). this leaf litter 
                                                   !<does not include litter generated due to mortality/fire
    real, allocatable, dimension(:,:,:) :: roottemp     !<root temperature, k
    real, allocatable, dimension(:,:,:) :: afrleaf      !<allocation fraction for leaves
    real, allocatable, dimension(:,:,:) :: afrstem      !<allocation fraction for stem
    real, allocatable, dimension(:,:,:) :: afrroot      !<allocation fraction for root
    real, allocatable, dimension(:,:,:) :: wtstatus     !<soil water status used for calculating allocation fractions
    real, allocatable, dimension(:,:,:) :: ltstatus     !<light status used for calculating allocation fractions
    real, allocatable, dimension(:,:,:) :: gppveg       !<!gross primary productity for each pft
    real, allocatable, dimension(:,:,:) :: nppveg       !<npp for individual pfts,  u-mol co2/m2.sec
    real, allocatable, dimension(:,:,:) :: autoresveg   !<
    real, allocatable, dimension(:,:,:) :: rmlvegacc    !<
    real, allocatable, dimension(:,:,:) :: rmsveg       !<stem maintenance resp. rate for each pft
    real, allocatable, dimension(:,:,:) :: rmrveg       !<root maintenance resp. rate for each pft
    real, allocatable, dimension(:,:,:) :: rgveg        !<growth resp. rate for each pft
    real, allocatable, dimension(:,:,:) :: litrfallveg  !<litter fall in \f$kg c/m^2\f$ for each pft
    real, allocatable, dimension(:,:,:) :: rothrlos     !<root death as crops are harvested, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:,:) :: pfcancmx     !<previous year's fractional coverages of pfts
    real, allocatable, dimension(:,:,:) :: nfcancmx     !<next year's fractional coverages of pfts
    real, allocatable, dimension(:,:,:) :: anveg        !<net photosynthesis rate for each pft
    real, allocatable, dimension(:,:,:) :: rmlveg       !<leaf maintenance resp. rate for each pft
    
! allocated with nlat,nmos:   
    integer, allocatable, dimension(:,:)     :: icount         !<
    integer, allocatable, dimension(:,:)     :: stdaln         !<an integer telling if ctem is operated within gcm (=0) or in stand
                                                        !<alone mode (=1). this is used for fire purposes. see comments just
                                                        !<above where disturb subroutine is called.
    real, allocatable, dimension(:,:) :: gavglai               !<grid averaged green leaf area index
    real, allocatable, dimension(:,:) :: co2conc               !<ATMOS. CO2 CONC. IN PPM
    real, allocatable, dimension(:,:) :: ch4conc               !<
    real, allocatable, dimension(:,:) :: canres                !<
    real, allocatable, dimension(:,:) :: vgbiomas              !<grid averaged vegetation biomass, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:) :: gavgltms              !<grid averaged litter mass, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:) :: gavgscms              !<grid averaged soil c mass, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:) :: burnfrac              !<areal fraction burned due to fire for every grid cell (%)
    real, allocatable, dimension(:,:) :: popdin                !<population density \f$(people / km^2)\f$
    real, allocatable, dimension(:,:) :: lterm                 !<lightning term for fire probabilty calc
    real, allocatable, dimension(:,:) :: extnprob              !<fire extingusinging probability
    real, allocatable, dimension(:,:) :: prbfrhuc              !<probability of fire due to human causes
    real, allocatable, dimension(:,:) :: rml                   !<leaf maintenance respiration (u-mol co2/m2.sec)
    real, allocatable, dimension(:,:) :: rms                   !<stem maintenance respiration (u-mol co2/m2.sec)
    real, allocatable, dimension(:,:) :: rmr                   !<root maintenance respiration (u-mol co2/m2.sec)
    real, allocatable, dimension(:,:) :: ch4wet1               !<methane flux from wetlands calculated using hetrores in umol ch4/m2.s
    real, allocatable, dimension(:,:) :: ch4wet2               !<methane flux from wetlands calculated using npp in umol ch4/m2.s
    real, allocatable, dimension(:,:) :: wetfdyn               !<dynamic wetland fraction
    real, allocatable, dimension(:,:) :: ch4dyn1               !<methane flux from wetlands calculated using hetrores 
                                                        !<and wetfdyn, in umol ch4/m2.s
    real, allocatable, dimension(:,:) :: ch4dyn2               !<methane flux from wetlands calculated using npp and wetfdyn, 
                                                        !<in umol ch4/m2.s
    real, allocatable, dimension(:,:) :: ch4_soills            !<Methane uptake into the soil column ($mg CH_4 m^{-2} s^{-1}$)
    real, allocatable, dimension(:,:) :: lucemcom              !<land use change (luc) related combustion emission losses, u-mol co2/m2.sec
    real, allocatable, dimension(:,:) :: lucltrin              !<luc related inputs to litter pool, u-mol co2/m2.sec
    real, allocatable, dimension(:,:) :: lucsocin              !<luc related inputs to soil c pool, u-mol co2/m2.sec
    real, allocatable, dimension(:,:) :: npp                   !<net primary productivity
    real, allocatable, dimension(:,:) :: nep                   !<net ecosystem productivity
    real, allocatable, dimension(:,:) :: nbp                   !<net biome productivity
    real, allocatable, dimension(:,:) :: gpp                   !<gross primary productivity
    real, allocatable, dimension(:,:) :: hetrores              !<heterotrophic respiration
    real, allocatable, dimension(:,:) :: autores               !<autotrophic respiration
    real, allocatable, dimension(:,:) :: soilcresp             !<
    real, allocatable, dimension(:,:) :: rm                    !<maintenance respiration
    real, allocatable, dimension(:,:) :: rg                    !<growth respiration
    real, allocatable, dimension(:,:) :: litres                !<litter respiration
    real, allocatable, dimension(:,:) :: socres                !<soil carbon respiration
    real, allocatable, dimension(:,:) :: dstcemls              !<carbon emission losses due to disturbance, mainly fire
    real, allocatable, dimension(:,:) :: litrfall              !<total litter fall (from leaves, stem, and root) due to 
                                                        !<all causes (mortality, turnover, and disturbance)
    real, allocatable, dimension(:,:) :: humiftrs              !<transfer of humidified litter from litter to soil c pool
    real, allocatable, dimension(:,:) :: cfluxcg               !<
    real, allocatable, dimension(:,:) :: cfluxcs               !<
    real, allocatable, dimension(:,:) :: dstcemls3             !<carbon emission losses due to disturbance (fire at present) from litter pool
    real, allocatable, dimension(:,:) :: PREACC_M              !<
    real, allocatable, dimension(:,:) :: GTACC_M               !<
    real, allocatable, dimension(:,:) :: QEVPACC_M             !<
    real, allocatable, dimension(:,:) :: HFSACC_M              !<
    real, allocatable, dimension(:,:) :: HMFNACC_M             !<
    real, allocatable, dimension(:,:) :: ROFACC_M              !<
    real, allocatable, dimension(:,:) :: SNOACC_M              !<
    real, allocatable, dimension(:,:) :: OVRACC_M              !<
    real, allocatable, dimension(:,:) :: WTBLACC_M             !<
    real, allocatable, dimension(:,:) :: ALVSACC_M             !<
    real, allocatable, dimension(:,:) :: ALIRACC_M             !<
    real, allocatable, dimension(:,:) :: RHOSACC_M             !<
    real, allocatable, dimension(:,:) :: TSNOACC_M             !<
    real, allocatable, dimension(:,:) :: WSNOACC_M             !<
    real, allocatable, dimension(:,:) :: SNOARE_M              !<
    real, allocatable, dimension(:,:) :: TCANACC_M             !<
    real, allocatable, dimension(:,:) :: RCANACC_M             !<
    real, allocatable, dimension(:,:) :: SCANACC_M             !<
    real, allocatable, dimension(:,:) :: ALTOTACC_M            !<Daily broadband albedo
    real, allocatable, dimension(:,:) :: GROACC_M              !<
    real, allocatable, dimension(:,:) :: FSINACC_M             !<
    real, allocatable, dimension(:,:) :: FLINACC_M             !<
    real, allocatable, dimension(:,:) :: TAACC_M               !<
    real, allocatable, dimension(:,:) :: UVACC_M               !<
    real, allocatable, dimension(:,:) :: PRESACC_M             !<
    real, allocatable, dimension(:,:) :: QAACC_M               !<
    real, allocatable, dimension(:,:) :: EVAPACC_M             !<
    real, allocatable, dimension(:,:) :: FLUTACC_M             !<
    real, allocatable, dimension(:,:) :: tcanrs                !<
    real, allocatable, dimension(:,:) :: tsnors                !<
    real, allocatable, dimension(:,:) :: tpndrs                !<
    real, allocatable, dimension(:,:) :: tcanoaccrow_m         !<
    real, allocatable, dimension(:,:) :: uvaccrow_m            !<
    real, allocatable, dimension(:,:) :: vvaccrow_m            !<
    real, allocatable, dimension(:,:) :: tcanoaccrow_out       !<
    real, allocatable, dimension(:,:) :: qevpacc_m_save        !<
    real, allocatable, dimension(:,:) :: twarmm                !< temperature of the warmest month (c)
    real, allocatable, dimension(:,:) :: tcoldm                !< temperature of the coldest month (c)
    real, allocatable, dimension(:,:) :: gdd5                  !< growing degree days above 5 c
    real, allocatable, dimension(:,:) :: aridity               !< aridity index, ratio of potential evaporation to precipitation
    real, allocatable, dimension(:,:) :: srplsmon              !< number of months in a year with surplus water i.e. precipitation more than potential evaporation
    real, allocatable, dimension(:,:) :: defctmon              !< number of months in a year with water deficit i.e. precipitation less than potential evaporation
    real, allocatable, dimension(:,:) :: anndefct              !< annual water deficit (mm)
    real, allocatable, dimension(:,:) :: annsrpls              !< annual water surplus (mm)
    real, allocatable, dimension(:,:) :: annpcp                !< annual precipitation (mm)
    real, allocatable, dimension(:,:) :: dry_season_length     !< length of dry season (months)

    
! allocated with nlat,nmos,ican:     
    real, allocatable, dimension(:,:,:) :: zolnc            !<lumped log of roughness length for class' 4 pfts
    real, allocatable, dimension(:,:,:) :: ailc             !<lumped lai for class' 4 pfts
    real, allocatable, dimension(:,:,:) :: cmasvegc         !<total canopy mass for each of the 4 class pfts. recall that
                                                        !<class requires canopy mass as an input, and this is now provided by ctem. \f$kg/m^2\f$.
    real, allocatable, dimension(:,:,:) :: alvsctm          !<
    real, allocatable, dimension(:,:,:) :: paic             !<plant area index for class' 4 pfts. this is the sum of leaf
                                                        !<area index and stem area index.
    real, allocatable, dimension(:,:,:) :: slaic            !<storage lai. this will be used as min. lai that class sees
                                                        !<so that it doesn't blow up in its stomatal conductance calculations.
    real, allocatable, dimension(:,:,:) :: alirctm          !<
    real, allocatable, dimension(:,:,:) :: csum             !<
    
! allocated with nlat,nmos,ignd:        
    real, allocatable, dimension(:,:,:) :: TBARACC_M        !<
    real, allocatable, dimension(:,:,:) :: THLQACC_M        !<
    real, allocatable, dimension(:,:,:) :: THICACC_M        !<
    real, allocatable, dimension(:,:,:) :: THALACC_M        !<
    real, allocatable, dimension(:,:,:) :: tbaraccrow_m     !<
    
! allocated with nlat,nmos,ican,ignd:       
    real, allocatable, dimension(:,:,:,:) :: rmatc       !<fraction of roots for each of class' 4 pfts in each soil layer
 
 ! allocated with nlat,nmos,icc,ignd: 
    real, allocatable, dimension(:,:,:,:) :: rmatctem     !<fraction of roots for each of ctem's 9 pfts in each soil layer
    
! allocated with nlat,nmos,iccp1:     
    real, allocatable, dimension(:,:,:) :: litrmass    !<litter mass for each of the 9 ctem pfts + bare, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:,:) :: soilcmas    !<soil carbon mass for each of the 9 ctem pfts + bare, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:,:) :: nepveg      !<net ecosystem productity for bare fraction expnbaln(i)=0.0 amount
                                                    !<of c related to spatial expansion Not used JM Jun 2014 
                                                    !<OR net ecosystem productity for each pft
    real, allocatable, dimension(:,:,:) :: nbpveg      !<net biome productity for bare fraction OR net biome productity for each pft
    real, allocatable, dimension(:,:,:) :: hetroresveg !<
    real, allocatable, dimension(:,:,:) :: litresveg   !<
    real, allocatable, dimension(:,:,:) :: soilcresveg !<
    real, allocatable, dimension(:,:,:) :: humiftrsveg !<
    
! allocated with nlat,nmos,{some number}: 
    integer, allocatable, dimension(:,:,:):: colddays          !<cold days counter for tracking days below a certain
                                                        !<temperature threshold for ndl dcd and crop pfts.
    real, allocatable, dimension(:,:,:)  :: slopefrac          !<prescribed fraction of wetlands based on slope
                                                        !<only(0.025, 0.05, 0.1, 0.15, 0.20, 0.25, 0.3 and 0.35 percent slope thresholds)
    real, allocatable, dimension(:,:,:) :: mlightng           !<
    real, allocatable, dimension(:,:,:) :: wetfrac_mon        !<
    
! allocated with nlat:     
    real, allocatable, dimension(:)    :: dayl_max        !< maximum daylength for that location (hours)
    real, allocatable, dimension(:)    :: dayl            !< daylength for that location (hours)
    integer, allocatable, dimension(:) :: altotcntr_d     !<Used to count the number of time steps with the sun above the horizon

end type veg_rot

type (veg_rot), allocatable, save, target :: vrot

!=================================================================================
!>CTEM's 'gat' vars
type veg_gat

    ! This is the basic data structure that contains the state variables
    ! for the Plant Functional Type (PFT). The dimensions are ilg,{icc,iccp1}

    real, allocatable,dimension(:,:) :: ailcmin     !<
    real, allocatable, dimension(:,:) :: ailcmax    !<
    real, allocatable, dimension(:,:) :: dvdfcan    !<
    real, allocatable, dimension(:,:) :: gleafmas   !<green leaf mass for each of the 9 ctem pfts, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:) :: bleafmas   !<brown leaf mass for each of the 9 ctem pfts, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:) :: stemmass   !<stem mass for each of the 9 ctem pfts, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:) :: rootmass   !<root mass for each of the 9 ctem pfts, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:) :: pstemmass  !<stem mass from previous timestep, is value before fire. used by burntobare subroutine
    real, allocatable, dimension(:,:) :: pgleafmass !<root mass from previous timestep, is value before fire. used by burntobare subroutine
    real, allocatable, dimension(:,:) :: fcancmx    !<max. fractional coverage of ctem's 9 pfts, but this can be
                                                    !<modified by land-use change, and competition between pfts
    real, allocatable, dimension(:) :: gavglai        !<grid averaged green leaf area index

    real, allocatable, dimension(:) :: lightng        !<total lightning frequency, flashes/km2.year
    real, allocatable, dimension(:) :: tcanoaccgat_out!<

    real, allocatable, dimension(:,:) :: zolnc     !<lumped log of roughness length for class' 4 pfts
    real, allocatable, dimension(:,:) :: ailc      !<lumped lai for class' 4 pfts

    real, allocatable, dimension(:,:) :: ailcg      !<green lai for ctem's 9 pfts
    real, allocatable, dimension(:,:) :: ailcgs     !<GREEN LAI FOR CANOPY OVER SNOW SUB-AREA
    real, allocatable, dimension(:,:) :: fcancs     !<FRACTION OF CANOPY OVER SNOW FOR CTEM's 9 PFTs
    real, allocatable, dimension(:,:) :: fcanc      !<FRACTIONAL COVERAGE OF 8 CARBON PFTs, CANOPY OVER GROUND

    real, allocatable, dimension(:)     :: co2conc    !<ATMOS. CO2 CONC. IN PPM
    real, allocatable, dimension(:)     :: ch4conc    !<

    real, allocatable, dimension(:,:) :: co2i1cg    !<INTERCELLULAR CO2 CONC FOR 8 PFTs FOR CANOPY OVER GROUND SUBAREA (Pa) - FOR SINGLE/SUNLIT LEAF
    real, allocatable, dimension(:,:) :: co2i1cs    !<SAME AS ABOVE BUT FOR SHADED LEAF (above being co2i1cg)
    real, allocatable, dimension(:,:) :: co2i2cg    !<INTERCELLULAR CO2 CONC FOR 8 PFTs FOR CANOPY OVER SNOWSUBAREA (Pa) - FOR SINGLE/SUNLIT LEAF
    real, allocatable, dimension(:,:) :: co2i2cs    !<SAME AS ABOVE BUT FOR SHADED LEAF (above being co2i2cg)
    real, allocatable, dimension(:,:) :: ancsveg    !<net photosynthetic rate for ctems 9 pfts for canopy over snow subarea
    real, allocatable, dimension(:,:) :: ancgveg    !<net photosynthetic rate for ctems 9 pfts for canopy over ground subarea
    real, allocatable, dimension(:,:) :: rmlcsveg   !<leaf respiration rate for ctems 9 pfts forcanopy over snow subarea
    real, allocatable, dimension(:,:) :: rmlcgveg   !<leaf respiration rate for ctems 9 pfts forcanopy over ground subarea
    real, allocatable, dimension(:,:) :: slai       !<storage/imaginary lai for phenology purposes
    real, allocatable, dimension(:,:) :: ailcb      !<brown lai for ctem's 9 pfts. for now we assume only grasses can have brown lai
    real, allocatable, dimension(:)   :: canres     !<
    real, allocatable, dimension(:,:) :: flhrloss   !<fall or harvest loss for deciduous trees and crops, respectively, \f$kg c/m^2\f$il1

    real, allocatable, dimension(:,:) :: grwtheff   !<growth efficiency. change in biomass per year per unit max.
                                                    !<lai (\f$kg c/m^2\f$)/(m2/m2), for use in mortality subroutine
    real, allocatable, dimension(:,:) :: lystmmas   !<stem mass at the end of last year
    real, allocatable, dimension(:,:) :: lyrotmas   !<root mass at the end of last year
    real, allocatable, dimension(:,:) :: tymaxlai   !<this year's maximum lai
    real, allocatable, dimension(:)   :: vgbiomas   !<grid averaged vegetation biomass, \f$kg c/m^2\f$
    real, allocatable, dimension(:)   :: gavgltms   !<grid averaged litter mass, \f$kg c/m^2\f$
    real, allocatable, dimension(:)   :: gavgscms   !<grid averaged soil c mass, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:) :: stmhrlos   !<stem harvest loss for crops, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:,:) :: rmatc    !<fraction of roots for each of class' 4 pfts in each soil layer
    real, allocatable, dimension(:,:,:) :: rmatctem !<fraction of roots for each of ctem's 9 pfts in each soil layer
    real, allocatable, dimension(:,:) :: litrmass   !<litter mass for each of the 9 ctem pfts + bare, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:) :: soilcmas   !<soil carbon mass for each of the 9 ctem pfts + bare, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:) :: vgbiomas_veg !<vegetation biomass for each pft

    real, allocatable, dimension(:,:) :: emit_co2   !<carbon dioxide
    real, allocatable, dimension(:,:) :: emit_co    !<carbon monoxide
    real, allocatable, dimension(:,:) :: emit_ch4   !<methane
    real, allocatable, dimension(:,:) :: emit_nmhc  !<non-methane hydrocarbons
    real, allocatable, dimension(:,:) :: emit_h2    !<hydrogen gas
    real, allocatable, dimension(:,:) :: emit_nox   !<nitrogen oxides
    real, allocatable, dimension(:,:) :: emit_n2o   !<nitrous oxide
    real, allocatable, dimension(:,:) :: emit_pm25  !<particulate matter less than 2.5 um in diameter
    real, allocatable, dimension(:,:) :: emit_tpm   !<total particulate matter
    real, allocatable, dimension(:,:) :: emit_tc    !<total carbon
    real, allocatable, dimension(:,:) :: emit_oc    !<organic carbon
    real, allocatable, dimension(:,:) :: emit_bc    !<black carbon
    real, allocatable, dimension(:)   :: burnfrac   !<areal fraction burned due to fire for every grid cell (%)
    real, allocatable, dimension(:,:) :: burnvegf   !<per PFT fraction burned of that PFT's area
    real, allocatable, dimension(:,:) :: smfuncveg  !<
    real, allocatable, dimension(:)   :: popdin     !<population density (people / \f$km^2\f$)
    real, allocatable, dimension(:,:) :: bterm      !<biomass term for fire probabilty calc
    real, allocatable, dimension(:)   :: lterm      !<lightning term for fire probabilty calc
    real, allocatable, dimension(:,:) :: mterm      !<moisture term for fire probabilty calc

    real, allocatable, dimension(:)     :: extnprob   !<fire extingusinging probability
    real, allocatable, dimension(:)     :: prbfrhuc   !<probability of fire due to human causes
    real, allocatable, dimension(:,:)   :: mlightng   !<
    real, allocatable, dimension(:)     :: dayl_max   !< maximum daylength for that location (hours)
    real, allocatable, dimension(:)     :: dayl       !< daylength for that location (hours)

    real, allocatable, dimension(:,:) :: bmasveg    !<total (gleaf + stem + root) biomass for each ctem pft, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:) :: cmasvegc   !<total canopy mass for each of the 4 class pfts. recall that
                                                    !<class requires canopy mass as an input, and this is now provided by ctem. \f$kg/m^2\f$.
    real, allocatable, dimension(:,:) :: veghght    !<vegetation height (meters)
    real, allocatable, dimension(:,:) :: rootdpth   !<99% soil rooting depth (meters)
                                                    !<both veghght & rootdpth can be used as diagnostics to see
                                                    !<how vegetation grows above and below ground, respectively
    real, allocatable, dimension(:)   :: rml        !<leaf maintenance respiration (u-mol co2/m2.sec)
    real, allocatable, dimension(:)   :: rms        !<stem maintenance respiration (u-mol co2/m2.sec)
    real, allocatable, dimension(:,:) :: tltrleaf   !<total leaf litter fall rate (u-mol co2/m2.sec)
    real, allocatable, dimension(:,:) :: tltrstem   !<total stem litter fall rate (u-mol co2/m2.sec)
    real, allocatable, dimension(:,:) :: tltrroot   !<total root litter fall rate (u-mol co2/m2.sec)
    real, allocatable, dimension(:,:) :: leaflitr   !<leaf litter fall rate (u-mol co2/m2.sec). this leaf litter
                                                    !<does not include litter generated due to mortality/fire
    real, allocatable, dimension(:,:) :: roottemp   !<root temperature, k
    real, allocatable, dimension(:,:) :: afrleaf    !<allocation fraction for leaves
    real, allocatable, dimension(:,:) :: afrstem    !<allocation fraction for stem
    real, allocatable, dimension(:,:) :: afrroot    !<allocation fraction for root
    real, allocatable, dimension(:,:) :: wtstatus   !<soil water status used for calculating allocation fractions
    real, allocatable, dimension(:,:) :: ltstatus   !<light status used for calculating allocation fractions
    real, allocatable, dimension(:)   :: rmr        !<root maintenance respiration (u-mol co2/m2.sec)

    real, allocatable, dimension(:,:) :: slopefrac      !<prescribed fraction of wetlands based on slope
                                                        !<only(0.025, 0.05, 0.1, 0.15, 0.20, 0.25, 0.3 and 0.35 percent slope thresholds)
    real, allocatable, dimension(:)    :: wetfrac_pres  !<
    real, allocatable, dimension(:,:)  :: wetfrac_mon   !<
    real, allocatable, dimension(:)    :: ch4wet1       !<methane flux from wetlands calculated using hetrores in umol ch4/m2.s
    real, allocatable, dimension(:)    :: ch4wet2       !<methane flux from wetlands calculated using npp in umol ch4/m2.s
    real, allocatable, dimension(:)    :: wetfdyn       !<dynamic wetland fraction
    real, allocatable, dimension(:)    :: ch4dyn1       !<methane flux from wetlands calculated using hetrores
                                                        !<and wetfdyn, in umol ch4/m2.s
    real, allocatable, dimension(:)    :: ch4dyn2       !<methane flux from wetlands calculated using npp and wetfdyn,
                                                        !<in umol ch4/m2.s
    real, allocatable, dimension(:)    :: ch4_soills    !<Methane uptake into the soil column ($mg CH_4 m^{-2} s^{-1}$)

    real, allocatable, dimension(:) :: lucemcom   !<land use change (luc) related combustion emission losses, u-mol co2/m2.sec
    real, allocatable, dimension(:) :: lucltrin   !<luc related inputs to litter pool, u-mol co2/m2.sec
    real, allocatable, dimension(:) :: lucsocin   !<luc related inputs to soil c pool, u-mol co2/m2.sec

    real, allocatable, dimension(:) :: npp        !<net primary productivity
    real, allocatable, dimension(:) :: nep        !<net ecosystem productivity
    real, allocatable, dimension(:) :: nbp        !<net biome productivity
    real, allocatable, dimension(:) :: gpp        !<gross primary productivity
    real, allocatable, dimension(:) :: hetrores   !<heterotrophic respiration
    real, allocatable, dimension(:) :: autores    !<autotrophic respiration
    real, allocatable, dimension(:) :: soilcresp  !<
    real, allocatable, dimension(:) :: rm         !<maintenance respiration
    real, allocatable, dimension(:) :: rg         !<growth respiration
    real, allocatable, dimension(:) :: litres     !<litter respiration
    real, allocatable, dimension(:) :: socres     !<soil carbon respiration
    real, allocatable, dimension(:) :: dstcemls   !<carbon emission losses due to disturbance, mainly fire
    real, allocatable, dimension(:) :: litrfall   !<total litter fall (from leaves, stem, and root) due to
                                                  !<all causes (mortality, turnover, and disturbance)
    real, allocatable, dimension(:) :: humiftrs   !<transfer of humidified litter from litter to soil c pool

    real, allocatable, dimension(:,:) :: gppveg     !<gross primary productity for each pft
    real, allocatable, dimension(:,:) :: nepveg     !<net ecosystem productity for bare fraction expnbaln(i)=0.0 amount
                                                    !<of c related to spatial expansion Not used JM Jun 2014
                                                    !<OR net ecosystem productity for each pft
    real, allocatable, dimension(:,:) :: nbpveg     !<net biome productity for bare fraction OR net biome productity for each pft
    real, allocatable, dimension(:,:) :: nppveg     !<npp for individual pfts,  u-mol co2/m2.sec
    real, allocatable, dimension(:,:) :: hetroresveg!<
    real, allocatable, dimension(:,:) :: autoresveg !<
    real, allocatable, dimension(:,:) :: litresveg  !<
    real, allocatable, dimension(:,:) :: soilcresveg!<
    real, allocatable, dimension(:,:) :: rmlvegacc  !<
    real, allocatable, dimension(:,:) :: rmsveg     !<stem maintenance resp. rate for each pft
    real, allocatable, dimension(:,:) :: rmrveg     !<root maintenance resp. rate for each pft
    real, allocatable, dimension(:,:) :: rgveg      !<growth resp. rate for each pft
    real, allocatable, dimension(:,:) :: litrfallveg!<litter fall in \f$kg c/m^2\f$ for each pft
    real, allocatable, dimension(:,:) :: humiftrsveg!<

    real, allocatable, dimension(:,:) :: rothrlos !<root death as crops are harvested, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:) :: pfcancmx !<previous year's fractional coverages of pfts
    real, allocatable, dimension(:,:) :: nfcancmx !<next year's fractional coverages of pfts
    real, allocatable, dimension(:,:) :: alvsctm  !<
    real, allocatable, dimension(:,:) :: paic     !<plant area index for class' 4 pfts. this is the sum of leaf
                                                  !<area index and stem area index.
    real, allocatable, dimension(:,:) :: slaic    !<storage lai. this will be used as min. lai that class sees
                                                  !<so that it doesn't blow up in its stomatal conductance calculations.
    real, allocatable, dimension(:,:) :: alirctm  !<
    real, allocatable, dimension(:)   :: cfluxcg  !<
    real, allocatable, dimension(:)   :: cfluxcs  !<
    real, allocatable, dimension(:)   :: dstcemls3!<carbon emission losses due to disturbance (fire at present) from litter pool
    real, allocatable, dimension(:,:) :: anveg    !<net photosynthesis rate for each pft
    real, allocatable, dimension(:,:) :: rmlveg   !<leaf maintenance resp. rate for each pft

    real, allocatable, dimension(:) :: twarmm            !< temperature of the warmest month (c)
    real, allocatable, dimension(:) :: tcoldm            !< temperature of the coldest month (c)
    real, allocatable, dimension(:) :: gdd5              !< growing degree days above 5 c
    real, allocatable, dimension(:) :: aridity           !< aridity index, ratio of potential evaporation to precipitation
    real, allocatable, dimension(:) :: srplsmon          !< number of months in a year with surplus water i.e. precipitation more than potential evaporation
    real, allocatable, dimension(:) :: defctmon          !< number of months in a year with water deficit i.e. precipitation less than potential evaporation
    real, allocatable, dimension(:) :: anndefct          !< annual water deficit (mm)
    real, allocatable, dimension(:) :: annsrpls          !< annual water surplus (mm)
    real, allocatable, dimension(:) :: annpcp            !< annual precipitation (mm)
    real, allocatable, dimension(:) :: dry_season_length !< length of dry season (months)

    ! These go into CTEM and are used to keep track of the bioclim limits.
    real, allocatable, dimension(:) :: tcurm     !<temperature of the current month (c)
    real, allocatable, dimension(:) :: srpcuryr  !<water surplus for the current year
    real, allocatable, dimension(:) :: dftcuryr  !<water deficit for the current year
    real, allocatable, dimension(:,:) :: tmonth  !<monthly temperatures
    real, allocatable, dimension(:) :: anpcpcur  !<annual precipitation for current year (mm)
    real, allocatable, dimension(:) :: anpecur   !<annual potential evaporation for current year (mm)
    real, allocatable, dimension(:) :: gdd5cur   !<growing degree days above 5 c for current year
    real, allocatable, dimension(:) :: surmncur  !<number of months with surplus water for current year
    real, allocatable, dimension(:) :: defmncur  !<number of months with water deficit for current year
    real, allocatable, dimension(:) :: srplscur  !<water surplus for the current month
    real, allocatable, dimension(:) :: defctcur  !<water deficit for the current month

    real, allocatable, dimension(:,:) :: geremort !<growth efficiency related mortality (1/day)
    real, allocatable, dimension(:,:) :: intrmort !<intrinsic (age related) mortality (1/day)
    real, allocatable, dimension(:,:) :: lambda   !<Used to determine the colonization rate
    real, allocatable, dimension(:,:) :: cc       !<colonization rate & mortality rate
    real, allocatable, dimension(:,:) :: mm       !<colonization rate & mortality rate

    logical, allocatable, dimension(:,:) :: pftexist !<logical array indicating pfts exist (t) or not (f)
    integer, allocatable, dimension(:,:) :: colddays !<cold days counter for tracking days below a certain
                                                     !<temperature threshold for ndl dcd and crop pfts.
    integer, allocatable, dimension(:)   :: icount   !<
    integer, allocatable, dimension(:,:) :: lfstatus !<leaf phenology status
    integer, allocatable, dimension(:,:) :: pandays  !<days with positive net photosynthesis (an) for use in
                                                     !<the phenology subroutine
    integer, allocatable, dimension(:)   :: stdaln   !<an integer telling if ctem is operated within gcm (=0) or in stand
                                                     !<alone mode (=1). this is used for fire purposes. see comments just
                                                     !<above where disturb subroutine is called.

end type veg_gat

type (veg_gat), allocatable, save, target :: vgat
!=================================================================================
!>CLASS's monthly outputs
type class_moyr_output

!   MONTHLY OUTPUT FOR CLASS GRID-MEAN

! allocated with nlat:
    real, allocatable, dimension(:) :: ALVSACC_MO   !<
    real, allocatable, dimension(:) :: ALIRACC_MO   !<
    real, allocatable, dimension(:) :: FLUTACC_MO   !<
    real, allocatable, dimension(:) :: FSINACC_MO   !<
    real, allocatable, dimension(:) :: FLINACC_MO   !<
    real, allocatable, dimension(:) :: HFSACC_MO    !<
    real, allocatable, dimension(:) :: QEVPACC_MO   !<
    real, allocatable, dimension(:) :: SNOACC_MO    !<
    real, allocatable, dimension(:) :: WSNOACC_MO   !<
    real, allocatable, dimension(:) :: ROFACC_MO    !<
    real, allocatable, dimension(:) :: PREACC_MO    !<
    real, allocatable, dimension(:) :: EVAPACC_MO   !<
    real, allocatable, dimension(:) :: TRANSPACC_MO !<
    real, allocatable, dimension(:) :: TAACC_MO     !<
    real, allocatable, dimension(:) :: ACTLYR_MO
    real, allocatable, dimension(:) :: FTABLE_MO
    real, allocatable, dimension(:) :: ACTLYR_MIN_MO
    real, allocatable, dimension(:) :: ACTLYR_MAX_MO
    real, allocatable, dimension(:) :: FTABLE_MIN_MO
    real, allocatable, dimension(:) :: FTABLE_MAX_MO
    real, allocatable, dimension(:) :: ALTOTACC_MO  !< Broadband albedo
    real, allocatable, dimension(:) :: GROUNDEVAP   !< evaporation and sublimation from the ground surface (formed from QFG and QFN), kg /m/mon
    real, allocatable, dimension(:) :: CANOPYEVAP   !< evaporation and sublimation from the canopy (formed from QFCL and QFCF), kg /m/mon
    integer, allocatable, dimension(:) :: altotcntr_m!< Used to count the number of time steps with the sun above the horizon

! allocated with nlat,ignd:
    real, allocatable, dimension(:,:) :: TBARACC_MO !<
    real, allocatable, dimension(:,:) :: THLQACC_MO !<
    real, allocatable, dimension(:,:) :: THICACC_MO !<

!   YEARLY OUTPUT FOR CLASS GRID-MEAN

    real, allocatable, dimension(:) :: ALVSACC_YR  !<
    real, allocatable, dimension(:) :: ALIRACC_YR  !<
    real, allocatable, dimension(:) :: FLUTACC_YR  !<
    real, allocatable, dimension(:) :: FSINACC_YR  !<
    real, allocatable, dimension(:) :: FLINACC_YR  !<
    real, allocatable, dimension(:) :: HFSACC_YR   !<
    real, allocatable, dimension(:) :: QEVPACC_YR  !<
    real, allocatable, dimension(:) :: ROFACC_YR   !<
    real, allocatable, dimension(:) :: PREACC_YR   !<
    real, allocatable, dimension(:) :: EVAPACC_YR  !<
    real, allocatable, dimension(:) :: TRANSPACC_YR!<
    real, allocatable, dimension(:) :: TAACC_YR    !<
    real, allocatable, dimension(:) :: ACTLYR_YR
    real, allocatable, dimension(:) :: ACTLYR_MIN_YR
    real, allocatable, dimension(:) :: ACTLYR_MAX_YR
    real, allocatable, dimension(:) :: FTABLE_YR
    real, allocatable, dimension(:) :: FTABLE_MIN_YR
    real, allocatable, dimension(:) :: FTABLE_MAX_YR
    real, allocatable, dimension(:) :: ALTOTACC_YR !< Broadband albedo
    integer, allocatable, dimension(:) :: altotcntr_yr !<Used to count the number of time steps with the sun above the horizon


end type class_moyr_output

type (class_moyr_output), allocatable, save, target :: class_out

!=================================================================================
!>CTEM's grid average variables
type ctem_gridavg

! Grid-averaged variables (denoted with an ending of "_g")

! allocated with nlat:
      real, allocatable, dimension(:) :: WSNOROT_g !<
      real, allocatable, dimension(:) :: ROFSROT_g !<
      real, allocatable, dimension(:) :: SNOROT_g  !<
      real, allocatable, dimension(:) :: RHOSROT_g !<
      real, allocatable, dimension(:) :: ROFROT_g  !<
      real, allocatable, dimension(:) :: ZPNDROT_g !<
      real, allocatable, dimension(:) :: RCANROT_g !<
      real, allocatable, dimension(:) :: SCANROT_g !<
      real, allocatable, dimension(:) :: TROFROT_g !<
      real, allocatable, dimension(:) :: TROOROT_g !<
      real, allocatable, dimension(:) :: TROBROT_g !<
      real, allocatable, dimension(:) :: ROFOROT_g !<
      real, allocatable, dimension(:) :: ROFBROT_g !<
      real, allocatable, dimension(:) :: TROSROT_g !<
      real, allocatable, dimension(:) :: FSGVROT_g !<
      real, allocatable, dimension(:) :: FSGSROT_g !<
      real, allocatable, dimension(:) :: FLGVROT_g !<
      real, allocatable, dimension(:) :: FLGSROT_g !<
      real, allocatable, dimension(:) :: HFSCROT_g !<
      real, allocatable, dimension(:) :: HFSSROT_g !<
      real, allocatable, dimension(:) :: HEVCROT_g !<
      real, allocatable, dimension(:) :: HEVSROT_g !<
      real, allocatable, dimension(:) :: HMFCROT_g !<
      real, allocatable, dimension(:) :: HMFNROT_g !<
      real, allocatable, dimension(:) :: HTCSROT_g !<
      real, allocatable, dimension(:) :: HTCCROT_g !<
      real, allocatable, dimension(:) :: FSGGROT_g !<
      real, allocatable, dimension(:) :: FLGGROT_g !<
      real, allocatable, dimension(:) :: HFSGROT_g !<
      real, allocatable, dimension(:) :: HEVGROT_g !<
      real, allocatable, dimension(:) :: CDHROT_g  !<
      real, allocatable, dimension(:) :: CDMROT_g  !<
      real, allocatable, dimension(:) :: SFCUROT_g !<
      real, allocatable, dimension(:) :: SFCVROT_g !<
      real, allocatable, dimension(:) :: ACTLYR_g
      real, allocatable, dimension(:) :: FTABLE_g
      real, allocatable, dimension(:) :: fc_g      !<
      real, allocatable, dimension(:) :: fg_g      !<
      real, allocatable, dimension(:) :: fcs_g     !<
      real, allocatable, dimension(:) :: fgs_g     !<
      real, allocatable, dimension(:) :: PCFCROT_g !<
      real, allocatable, dimension(:) :: PCLCROT_g !<
      real, allocatable, dimension(:) :: PCPGROT_g !<
      real, allocatable, dimension(:) :: QFCFROT_g !<
      real, allocatable, dimension(:) :: QFGROT_g  !<
      real, allocatable, dimension(:,:) :: QFCROT_g !<
      real, allocatable, dimension(:) :: ROFCROT_g  !<
      real, allocatable, dimension(:) :: ROFNROT_g  !<
      real, allocatable, dimension(:) :: WTRSROT_g  !<
      real, allocatable, dimension(:) :: WTRGROT_g  !<
      real, allocatable, dimension(:) :: PCPNROT_g  !<
      real, allocatable, dimension(:) :: QFCLROT_g  !<
      real, allocatable, dimension(:) :: QFNROT_g   !<
      real, allocatable, dimension(:) :: WTRCROT_g  !<
      real, allocatable, dimension(:) :: gpp_g      !<
      real, allocatable, dimension(:) :: npp_g      !<
      real, allocatable, dimension(:) :: nbp_g      !<
      real, allocatable, dimension(:) :: socres_g   !<
      real, allocatable, dimension(:) :: autores_g  !<
      real, allocatable, dimension(:) :: litres_g   !<
      real, allocatable, dimension(:) :: dstcemls3_g!<
      real, allocatable, dimension(:) :: litrfall_g !<
      real, allocatable, dimension(:) :: rml_g      !<
      real, allocatable, dimension(:) :: rms_g      !<
      real, allocatable, dimension(:) :: rg_g       !<
      real, allocatable, dimension(:) :: leaflitr_g !<
      real, allocatable, dimension(:) :: tltrstem_g !<
      real, allocatable, dimension(:) :: tltrroot_g !<
      real, allocatable, dimension(:) :: nep_g      !<
      real, allocatable, dimension(:) :: hetrores_g !<
      real, allocatable, dimension(:) :: dstcemls_g !<
      real, allocatable, dimension(:) :: humiftrs_g !<
      real, allocatable, dimension(:) :: rmr_g      !<
      real, allocatable, dimension(:) :: tltrleaf_g !<
      real, allocatable, dimension(:) :: gavgltms_g !<

      real, allocatable, dimension(:) :: vgbiomas_g !<
      real, allocatable, dimension(:) :: gavglai_g  !<
      real, allocatable, dimension(:) :: gavgscms_g !<
      real, allocatable, dimension(:) :: gleafmas_g !<
      real, allocatable, dimension(:) :: bleafmas_g !<
      real, allocatable, dimension(:) :: stemmass_g !<
      real, allocatable, dimension(:) :: rootmass_g !<
      real, allocatable, dimension(:) :: litrmass_g !<
      real, allocatable, dimension(:) :: soilcmas_g !<
      real, allocatable, dimension(:) :: slai_g     !<
      real, allocatable, dimension(:) :: ailcg_g    !<
      real, allocatable, dimension(:) :: ailcb_g    !<
      real, allocatable, dimension(:) :: veghght_g  !<
      real, allocatable, dimension(:) :: rootdpth_g !<
      real, allocatable, dimension(:) :: roottemp_g !<
      real, allocatable, dimension(:) :: totcmass_g !<
      real, allocatable, dimension(:) :: tcanoacc_out_g!<
      real, allocatable, dimension(:) :: burnfrac_g !<
      real, allocatable, dimension(:) :: smfuncveg_g!<
      real, allocatable, dimension(:) :: lucemcom_g !<
      real, allocatable, dimension(:) :: lucltrin_g !<
      real, allocatable, dimension(:) :: lucsocin_g !<
      real, allocatable, dimension(:) :: emit_co2_g !<
      real, allocatable, dimension(:) :: emit_co_g  !<
      real, allocatable, dimension(:) :: emit_ch4_g !<
      real, allocatable, dimension(:) :: emit_nmhc_g!<
      real, allocatable, dimension(:) :: emit_h2_g  !<
      real, allocatable, dimension(:) :: emit_nox_g !<
      real, allocatable, dimension(:) :: emit_n2o_g !<
      real, allocatable, dimension(:) :: emit_pm25_g!<
      real, allocatable, dimension(:) :: emit_tpm_g !<
      real, allocatable, dimension(:) :: emit_tc_g  !<
      real, allocatable, dimension(:) :: emit_oc_g  !<
      real, allocatable, dimension(:) :: emit_bc_g  !<
      real, allocatable, dimension(:) :: bterm_g    !<
      real, allocatable, dimension(:) :: lterm_g    !<
      real, allocatable, dimension(:) :: mterm_g    !<
      real, allocatable, dimension(:) :: ch4wet1_g  !<
      real, allocatable, dimension(:) :: ch4wet2_g  !<
      real, allocatable, dimension(:) :: wetfdyn_g  !<
      real, allocatable, dimension(:) :: ch4dyn1_g  !<
      real, allocatable, dimension(:) :: ch4dyn2_g  !<
      real, allocatable, dimension(:) :: ch4_soills_g   !<
      real, allocatable, dimension(:,:) :: afrleaf_g  !<
      real, allocatable, dimension(:,:) :: afrstem_g  !<
      real, allocatable, dimension(:,:) :: afrroot_g  !<
      real, allocatable, dimension(:,:) :: lfstatus_g !<
      real, allocatable, dimension(:,:) :: rmlvegrow_g!<
      real, allocatable, dimension(:,:) :: anvegrow_g !<
      real, allocatable, dimension(:,:) :: rmatctem_g!<
      real, allocatable, dimension(:,:) :: HMFGROT_g !<
      real, allocatable, dimension(:,:) :: HTCROT_g  !<
      real, allocatable, dimension(:,:) :: TBARROT_g !<
      real, allocatable, dimension(:,:) :: THLQROT_g !<
      real, allocatable, dimension(:,:) :: THICROT_g !<
      real, allocatable, dimension(:,:) :: GFLXROT_g !<

      real, allocatable, dimension(:) :: fsstar_g !<
      real, allocatable, dimension(:) :: flstar_g !<
      real, allocatable, dimension(:) :: qh_g     !<
      real, allocatable, dimension(:) :: qe_g     !<
      real, allocatable, dimension(:) :: snomlt_g !<
      real, allocatable, dimension(:) :: beg_g    !<
      real, allocatable, dimension(:) :: gtout_g  !<
      real, allocatable, dimension(:) :: tpn_g    !<
      real, allocatable, dimension(:) :: altot_g  !<
      real, allocatable, dimension(:) :: tcn_g    !<
      real, allocatable, dimension(:) :: tsn_g    !<
      real, allocatable, dimension(:) :: zsn_g    !<


end type ctem_gridavg

type (ctem_gridavg), allocatable, save, target :: ctem_grd

!=================================================================================
!>CTEM's variables per tile
type ctem_tile_level

!   Tile-level variables (denoted by an ending of "_t")

! allocated with nlat,nmos:
      real, allocatable, dimension(:,:) :: leaflitr_t !<
      real, allocatable, dimension(:,:) :: tltrleaf_t !<
      real, allocatable, dimension(:,:) :: tltrstem_t !<
      real, allocatable, dimension(:,:) :: tltrroot_t !<
      real, allocatable, dimension(:,:) :: ailcg_t    !<
      real, allocatable, dimension(:,:) :: ailcb_t    !<
      real, allocatable, dimension(:,:,:) :: rmatctem_t !< nlat,nmos,ignd
      real, allocatable, dimension(:,:) :: veghght_t  !<
      real, allocatable, dimension(:,:) :: rootdpth_t !<
      real, allocatable, dimension(:,:) :: roottemp_t !<
      real, allocatable, dimension(:,:) :: slai_t     !<
      real, allocatable, dimension(:,:) :: afrroot_t  !<
      real, allocatable, dimension(:,:) :: afrleaf_t  !<
      real, allocatable, dimension(:,:) :: afrstem_t  !<
      real, allocatable, dimension(:,:) :: laimaxg_t  !<
      real, allocatable, dimension(:,:) :: stemmass_t !<
      real, allocatable, dimension(:,:) :: rootmass_t !<
      real, allocatable, dimension(:,:) :: litrmass_t !<
      real, allocatable, dimension(:,:) :: gleafmas_t !<
      real, allocatable, dimension(:,:) :: bleafmas_t !<
      real, allocatable, dimension(:,:) :: soilcmas_t !<
      real, allocatable, dimension(:,:) :: emit_co2_t !<
      real, allocatable, dimension(:,:) :: emit_co_t  !<
      real, allocatable, dimension(:,:) :: emit_ch4_t !<
      real, allocatable, dimension(:,:) :: emit_nmhc_t!<
      real, allocatable, dimension(:,:) :: emit_h2_t  !<
      real, allocatable, dimension(:,:) :: emit_nox_t !<
      real, allocatable, dimension(:,:) :: emit_n2o_t !<
      real, allocatable, dimension(:,:) :: emit_pm25_t!<
      real, allocatable, dimension(:,:) :: emit_tpm_t !<
      real, allocatable, dimension(:,:) :: emit_tc_t  !<
      real, allocatable, dimension(:,:) :: emit_oc_t  !<
      real, allocatable, dimension(:,:) :: emit_bc_t  !<
      real, allocatable, dimension(:,:) :: bterm_t    !<
      real, allocatable, dimension(:,:) :: mterm_t    !<
      real, allocatable, dimension(:,:) :: smfuncveg_t!<

      real, allocatable, dimension(:) :: fsnowacc_t       !<
      real, allocatable, dimension(:) :: tcansacc_t       !<
      real, allocatable, dimension(:) :: tcanoaccgat_t    !<
      real, allocatable, dimension(:) :: taaccgat_t       !<
      real, allocatable, dimension(:) :: uvaccgat_t       !<
      real, allocatable, dimension(:) :: vvaccgat_t       !<

! allocated with ilg,ignd:
      real, allocatable, dimension(:,:) :: tbaraccgat_t!<
      real, allocatable, dimension(:,:) :: tbarcacc_t  !<
      real, allocatable, dimension(:,:) :: tbarcsacc_t !<
      real, allocatable, dimension(:,:) :: tbargacc_t  !<
      real, allocatable, dimension(:,:) :: tbargsacc_t !<
      real, allocatable, dimension(:,:) :: thliqcacc_t !<
      real, allocatable, dimension(:,:) :: thliqgacc_t !<
      real, allocatable, dimension(:,:) :: thliqacc_t  !<
      real, allocatable, dimension(:,:) :: thicecacc_t !<
      real, allocatable, dimension(:,:) :: thicegacc_t !<

! allocated with ilg,icc:
      real, allocatable, dimension(:,:) :: ancsvgac_t  !<
      real, allocatable, dimension(:,:) :: ancgvgac_t  !<
      real, allocatable, dimension(:,:) :: rmlcsvga_t  !<
      real, allocatable, dimension(:,:) :: rmlcgvga_t  !<

end type ctem_tile_level

type (ctem_tile_level), allocatable, save, target :: ctem_tile

!=================================================================================
!>CTEM's variables monthly averaged (per pft)
type ctem_monthly

!     Tile-level monthly variables (denoted by name ending in "_mo_t")

! allocated with nlat,nmos,icc/iccp1:

      real, allocatable, dimension(:,:,:)   :: laimaxg_mo    !<
      real, allocatable, dimension(:,:,:)   :: stemmass_mo   !<
      real, allocatable, dimension(:,:,:)   :: rootmass_mo   !<
      real, allocatable, dimension(:,:,:)   :: litrfallveg_mo!<
      real, allocatable, dimension(:,:,:) :: humiftrsveg_mo!<
      real, allocatable, dimension(:,:,:)   :: npp_mo        !<
      real, allocatable, dimension(:,:,:)   :: gpp_mo        !<
      real, allocatable, dimension(:,:,:)   :: vgbiomas_mo   !<
      real, allocatable, dimension(:,:,:)   :: autores_mo    !<
      real, allocatable, dimension(:,:,:) :: totcmass_mo   !<
      real, allocatable, dimension(:,:,:) :: litrmass_mo   !<
      real, allocatable, dimension(:,:,:) :: soilcmas_mo   !<
      real, allocatable, dimension(:,:,:) :: nep_mo        !<
      real, allocatable, dimension(:,:,:) :: litres_mo     !<
      real, allocatable, dimension(:,:,:) :: soilcres_mo   !<
      real, allocatable, dimension(:,:,:) :: hetrores_mo   !<
      real, allocatable, dimension(:,:,:) :: nbp_mo        !<
      real, allocatable, dimension(:,:,:) :: emit_co2_mo  !<
      real, allocatable, dimension(:,:,:) :: emit_co_mo   !<
      real, allocatable, dimension(:,:,:) :: emit_ch4_mo  !<
      real, allocatable, dimension(:,:,:) :: emit_nmhc_mo !<
      real, allocatable, dimension(:,:,:) :: emit_h2_mo   !<
      real, allocatable, dimension(:,:,:) :: emit_nox_mo  !<
      real, allocatable, dimension(:,:,:) :: emit_n2o_mo  !<
      real, allocatable, dimension(:,:,:) :: emit_pm25_mo !<
      real, allocatable, dimension(:,:,:) :: emit_tpm_mo  !<
      real, allocatable, dimension(:,:,:) :: emit_tc_mo   !<
      real, allocatable, dimension(:,:,:) :: emit_oc_mo   !<
      real, allocatable, dimension(:,:,:) :: emit_bc_mo   !<
      real, allocatable, dimension(:,:,:) :: burnfrac_mo  !<
      real, allocatable, dimension(:,:,:) :: bterm_mo     !<
      real, allocatable, dimension(:,:,:) :: mterm_mo     !<
      real, allocatable, dimension(:,:,:) :: smfuncveg_mo !<

end type ctem_monthly

type (ctem_monthly), allocatable, save, target :: ctem_mo

!=================================================================================
!>CTEM's grid average monthly values
type ctem_gridavg_monthly

!  Grid averaged monthly variables (denoted by name ending in "_mo_g")

! allocated with nlat:
    real, allocatable, dimension(:) :: laimaxg_mo_g  !<
    real, allocatable, dimension(:) :: stemmass_mo_g !<
    real, allocatable, dimension(:) :: rootmass_mo_g !<
    real, allocatable, dimension(:) :: litrmass_mo_g !<
    real, allocatable, dimension(:) :: soilcmas_mo_g !<
    real, allocatable, dimension(:) :: litrfall_mo_g !<
    real, allocatable, dimension(:) :: humiftrs_mo_g !<
    real, allocatable, dimension(:) :: npp_mo_g      !<
    real, allocatable, dimension(:) :: gpp_mo_g      !<
    real, allocatable, dimension(:) :: nep_mo_g      !<
    real, allocatable, dimension(:) :: nbp_mo_g      !<
    real, allocatable, dimension(:) :: hetrores_mo_g !<
    real, allocatable, dimension(:) :: autores_mo_g  !<
    real, allocatable, dimension(:) :: litres_mo_g   !<
    real, allocatable, dimension(:) :: soilcres_mo_g !<
    real, allocatable, dimension(:) :: vgbiomas_mo_g !<
    real, allocatable, dimension(:) :: totcmass_mo_g !<
    real, allocatable, dimension(:) :: emit_co2_mo_g !<
    real, allocatable, dimension(:) :: emit_co_mo_g  !<
    real, allocatable, dimension(:) :: emit_ch4_mo_g !<
    real, allocatable, dimension(:) :: emit_nmhc_mo_g!<
    real, allocatable, dimension(:) :: emit_h2_mo_g  !<
    real, allocatable, dimension(:) :: emit_nox_mo_g !<
    real, allocatable, dimension(:) :: emit_n2o_mo_g !<
    real, allocatable, dimension(:) :: emit_pm25_mo_g!<
    real, allocatable, dimension(:) :: emit_tpm_mo_g !<
    real, allocatable, dimension(:) :: emit_tc_mo_g  !<
    real, allocatable, dimension(:) :: emit_oc_mo_g  !<
    real, allocatable, dimension(:) :: emit_bc_mo_g  !<
    real, allocatable, dimension(:) :: smfuncveg_mo_g!<
    real, allocatable, dimension(:) :: luc_emc_mo_g  !<
    real, allocatable, dimension(:) :: lucltrin_mo_g !<
    real, allocatable, dimension(:) :: lucsocin_mo_g !<
    real, allocatable, dimension(:) :: burnfrac_mo_g !<
    real, allocatable, dimension(:) :: bterm_mo_g    !<
    real, allocatable, dimension(:) :: lterm_mo_g    !<
    real, allocatable, dimension(:) :: mterm_mo_g    !<
    real, allocatable, dimension(:) :: ch4wet1_mo_g  !<
    real, allocatable, dimension(:) :: ch4wet2_mo_g  !<
    real, allocatable, dimension(:) :: wetfdyn_mo_g  !<
    real, allocatable, dimension(:) :: ch4dyn1_mo_g  !<
    real, allocatable, dimension(:) :: ch4dyn2_mo_g  !<
    real, allocatable, dimension(:) :: ch4soills_mo_g!<

end type ctem_gridavg_monthly

type (ctem_gridavg_monthly), allocatable, save, target :: ctem_grd_mo

!=================================================================================
!>CTEM's variables per tile monthly values
type ctem_tileavg_monthly

!     Tile-level monthly variables (denoted by name ending in "_mo_t")

! allocated with nlat,nmos:

      real, allocatable, dimension(:,:) :: laimaxg_mo_t  !<
      real, allocatable, dimension(:,:) :: stemmass_mo_t !<
      real, allocatable, dimension(:,:) :: rootmass_mo_t !<
      real, allocatable, dimension(:,:) :: litrfall_mo_t !<
      real, allocatable, dimension(:,:) :: humiftrs_mo_t !<
      real, allocatable, dimension(:,:) :: npp_mo_t      !<
      real, allocatable, dimension(:,:) :: gpp_mo_t      !<
      real, allocatable, dimension(:,:) :: vgbiomas_mo_t !<
      real, allocatable, dimension(:,:) :: autores_mo_t  !<
      real, allocatable, dimension(:,:) :: totcmass_mo_t !<
      real, allocatable, dimension(:,:) :: litrmass_mo_t !<
      real, allocatable, dimension(:,:) :: soilcmas_mo_t !<
      real, allocatable, dimension(:,:) :: nep_mo_t      !<
      real, allocatable, dimension(:,:) :: litres_mo_t   !<
      real, allocatable, dimension(:,:) :: soilcres_mo_t !<
      real, allocatable, dimension(:,:) :: hetrores_mo_t !<
      real, allocatable, dimension(:,:) :: nbp_mo_t      !<
      real, allocatable, dimension(:,:) :: emit_co2_mo_t !<
      real, allocatable, dimension(:,:) :: emit_co_mo_t  !<
      real, allocatable, dimension(:,:) :: emit_ch4_mo_t !<
      real, allocatable, dimension(:,:) :: emit_nmhc_mo_t!<
      real, allocatable, dimension(:,:) :: emit_h2_mo_t  !<
      real, allocatable, dimension(:,:) :: emit_nox_mo_t !<
      real, allocatable, dimension(:,:) :: emit_n2o_mo_t !<
      real, allocatable, dimension(:,:) :: emit_pm25_mo_t!<
      real, allocatable, dimension(:,:) :: emit_tpm_mo_t !<
      real, allocatable, dimension(:,:) :: emit_tc_mo_t  !<
      real, allocatable, dimension(:,:) :: emit_oc_mo_t  !<
      real, allocatable, dimension(:,:) :: emit_bc_mo_t  !<
      real, allocatable, dimension(:,:) :: burnfrac_mo_t !<
      real, allocatable, dimension(:,:) :: smfuncveg_mo_t!<
      real, allocatable, dimension(:,:) :: bterm_mo_t    !<
      real, allocatable, dimension(:,:) :: luc_emc_mo_t  !<
      real, allocatable, dimension(:,:) :: lterm_mo_t    !<
      real, allocatable, dimension(:,:) :: lucsocin_mo_t !<
      real, allocatable, dimension(:,:) :: mterm_mo_t    !<
      real, allocatable, dimension(:,:) :: lucltrin_mo_t !<
      real, allocatable, dimension(:,:) :: ch4wet1_mo_t  !<
      real, allocatable, dimension(:,:) :: ch4wet2_mo_t  !<
      real, allocatable, dimension(:,:) :: wetfdyn_mo_t  !<
      real, allocatable, dimension(:,:) :: ch4dyn1_mo_t  !<
      real, allocatable, dimension(:,:) :: ch4dyn2_mo_t  !<
      real, allocatable, dimension(:,:) :: ch4soills_mo_t!<
      real, allocatable, dimension(:,:) :: wind_mo_t     !<

end type ctem_tileavg_monthly

type (ctem_tileavg_monthly), allocatable, save, target :: ctem_tile_mo


!=================================================================================

!>CTEM's average annual values (per PFT)
type ctem_annual

! c      Annual output for CTEM mosaic variables:
! c      (denoted by name ending in "_yr_m")
!
! allocated with nlat,nmos,icc/iccp1:

      real, allocatable, dimension(:,:,:) :: laimaxg_yr   !<
      real, allocatable, dimension(:,:,:) :: stemmass_yr  !<
      real, allocatable, dimension(:,:,:) :: rootmass_yr  !<
      real, allocatable, dimension(:,:,:) :: npp_yr       !<
      real, allocatable, dimension(:,:,:) :: gpp_yr       !<
      real, allocatable, dimension(:,:,:) :: vgbiomas_yr  !<
      real, allocatable, dimension(:,:,:) :: autores_yr   !<
      real, allocatable, dimension(:,:,:) :: totcmass_yr!<
      real, allocatable, dimension(:,:,:) :: litrmass_yr!<
      real, allocatable, dimension(:,:,:) :: soilcmas_yr!<
      real, allocatable, dimension(:,:,:) :: nep_yr     !<
      real, allocatable, dimension(:,:,:) :: litres_yr  !<
      real, allocatable, dimension(:,:,:) :: soilcres_yr!<
      real, allocatable, dimension(:,:,:) :: hetrores_yr!<
      real, allocatable, dimension(:,:,:) :: nbp_yr     !<
      real, allocatable, dimension(:,:,:) :: emit_co2_yr  !<
      real, allocatable, dimension(:,:,:) :: emit_co_yr   !<
      real, allocatable, dimension(:,:,:) :: emit_ch4_yr  !<
      real, allocatable, dimension(:,:,:) :: emit_nmhc_yr !<
      real, allocatable, dimension(:,:,:) :: emit_h2_yr   !<
      real, allocatable, dimension(:,:,:) :: emit_nox_yr  !<
      real, allocatable, dimension(:,:,:) :: emit_n2o_yr  !<
      real, allocatable, dimension(:,:,:) :: emit_pm25_yr !<
      real, allocatable, dimension(:,:,:) :: emit_tpm_yr  !<
      real, allocatable, dimension(:,:,:) :: emit_tc_yr   !<
      real, allocatable, dimension(:,:,:) :: emit_oc_yr   !<
      real, allocatable, dimension(:,:,:) :: emit_bc_yr   !<
      real, allocatable, dimension(:,:,:) :: bterm_yr     !<
      real, allocatable, dimension(:,:,:) :: mterm_yr     !<
      real, allocatable, dimension(:,:,:) :: burnfrac_yr  !<
      real, allocatable, dimension(:,:,:) :: smfuncveg_yr !<
      real, allocatable, dimension(:,:,:) :: veghght_yr   !<

end type ctem_annual

type (ctem_annual), allocatable, save, target :: ctem_yr

!=================================================================================

!>CTEM's grid average annual values
type ctem_gridavg_annual

! Annual output for CTEM grid-averaged variables:
! (denoted by name ending in "_yr_g")

! allocated with nlat:
    real, allocatable, dimension(:) :: laimaxg_yr_g  !<
    real, allocatable, dimension(:) :: stemmass_yr_g !<
    real, allocatable, dimension(:) :: rootmass_yr_g !<
    real, allocatable, dimension(:) :: litrmass_yr_g !<
    real, allocatable, dimension(:) :: soilcmas_yr_g !<
    real, allocatable, dimension(:) :: npp_yr_g      !<
    real, allocatable, dimension(:) :: gpp_yr_g      !<
    real, allocatable, dimension(:) :: nep_yr_g      !<
    real, allocatable, dimension(:) :: nbp_yr_g      !<
    real, allocatable, dimension(:) :: hetrores_yr_g !<
    real, allocatable, dimension(:) :: autores_yr_g  !<
    real, allocatable, dimension(:) :: litres_yr_g   !<
    real, allocatable, dimension(:) :: soilcres_yr_g !<
    real, allocatable, dimension(:) :: vgbiomas_yr_g !<
    real, allocatable, dimension(:) :: totcmass_yr_g !<
    real, allocatable, dimension(:) :: emit_co2_yr_g !<
    real, allocatable, dimension(:) :: emit_co_yr_g  !<
    real, allocatable, dimension(:) :: emit_ch4_yr_g !<
    real, allocatable, dimension(:) :: emit_nmhc_yr_g!<
    real, allocatable, dimension(:) :: emit_h2_yr_g  !<
    real, allocatable, dimension(:) :: emit_nox_yr_g !<
    real, allocatable, dimension(:) :: emit_n2o_yr_g !<
    real, allocatable, dimension(:) :: emit_pm25_yr_g!<
    real, allocatable, dimension(:) :: emit_tpm_yr_g !<
    real, allocatable, dimension(:) :: emit_tc_yr_g  !<
    real, allocatable, dimension(:) :: emit_oc_yr_g  !<
    real, allocatable, dimension(:) :: emit_bc_yr_g  !<
    real, allocatable, dimension(:) :: smfuncveg_yr_g!<
    real, allocatable, dimension(:) :: luc_emc_yr_g  !<
    real, allocatable, dimension(:) :: lucltrin_yr_g !<
    real, allocatable, dimension(:) :: lucsocin_yr_g !<
    real, allocatable, dimension(:) :: burnfrac_yr_g !<
    real, allocatable, dimension(:) :: bterm_yr_g    !<
    real, allocatable, dimension(:) :: lterm_yr_g    !<
    real, allocatable, dimension(:) :: mterm_yr_g    !<
    real, allocatable, dimension(:) :: ch4wet1_yr_g  !<
    real, allocatable, dimension(:) :: ch4wet2_yr_g  !<
    real, allocatable, dimension(:) :: wetfdyn_yr_g  !<
    real, allocatable, dimension(:) :: ch4dyn1_yr_g  !<
    real, allocatable, dimension(:) :: ch4dyn2_yr_g  !<
    real, allocatable, dimension(:) :: ch4soills_yr_g!<
    real, allocatable, dimension(:) :: veghght_yr_g  !<

end type ctem_gridavg_annual

type (ctem_gridavg_annual), allocatable, save, target :: ctem_grd_yr

!=================================================================================
!>CTEM's variables per tile annual values
type ctem_tileavg_annual

! c      Annual output for CTEM mosaic variables:
! c      (denoted by name ending in "_yr_m")
!
! allocated with nlat,nmos:
      real, allocatable, dimension(:,:) :: laimaxg_yr_t  !<
      real, allocatable, dimension(:,:) :: stemmass_yr_t !<
      real, allocatable, dimension(:,:) :: rootmass_yr_t !<
      real, allocatable, dimension(:,:) :: npp_yr_t      !<
      real, allocatable, dimension(:,:) :: gpp_yr_t      !<
      real, allocatable, dimension(:,:) :: vgbiomas_yr_t !<
      real, allocatable, dimension(:,:) :: autores_yr_t  !<
      real, allocatable, dimension(:,:) :: totcmass_yr_t !<
      real, allocatable, dimension(:,:) :: litrmass_yr_t !<
      real, allocatable, dimension(:,:) :: soilcmas_yr_t !<
      real, allocatable, dimension(:,:) :: nep_yr_t      !<
      real, allocatable, dimension(:,:) :: litres_yr_t   !<
      real, allocatable, dimension(:,:) :: soilcres_yr_t !<
      real, allocatable, dimension(:,:) :: hetrores_yr_t !<
      real, allocatable, dimension(:,:) :: nbp_yr_t      !<
      real, allocatable, dimension(:,:) :: emit_co2_yr_t !<
      real, allocatable, dimension(:,:) :: emit_co_yr_t  !<
      real, allocatable, dimension(:,:) :: emit_ch4_yr_t !<
      real, allocatable, dimension(:,:) :: emit_nmhc_yr_t!<
      real, allocatable, dimension(:,:) :: emit_h2_yr_t  !<
      real, allocatable, dimension(:,:) :: emit_nox_yr_t !<
      real, allocatable, dimension(:,:) :: emit_n2o_yr_t !<
      real, allocatable, dimension(:,:) :: emit_pm25_yr_t!<
      real, allocatable, dimension(:,:) :: emit_tpm_yr_t !<
      real, allocatable, dimension(:,:) :: emit_tc_yr_t  !<
      real, allocatable, dimension(:,:) :: emit_oc_yr_t  !<
      real, allocatable, dimension(:,:) :: emit_bc_yr_t  !<
      real, allocatable, dimension(:,:) :: burnfrac_yr_t !<
      real, allocatable, dimension(:,:) :: smfuncveg_yr_t!<
      real, allocatable, dimension(:,:) :: bterm_yr_t    !<
      real, allocatable, dimension(:,:) :: luc_emc_yr_t  !<
      real, allocatable, dimension(:,:) :: lterm_yr_t    !<
      real, allocatable, dimension(:,:) :: lucsocin_yr_t !<
      real, allocatable, dimension(:,:) :: mterm_yr_t    !<
      real, allocatable, dimension(:,:) :: lucltrin_yr_t !<
      real, allocatable, dimension(:,:) :: ch4wet1_yr_t  !<
      real, allocatable, dimension(:,:) :: ch4wet2_yr_t  !<
      real, allocatable, dimension(:,:) :: wetfdyn_yr_t  !<
      real, allocatable, dimension(:,:) :: ch4dyn1_yr_t  !<
      real, allocatable, dimension(:,:) :: ch4dyn2_yr_t  !<
      real, allocatable, dimension(:,:) :: ch4soills_yr_t!<
      real, allocatable, dimension(:,:) :: veghght_yr_t  !<

end type ctem_tileavg_annual

type (ctem_tileavg_annual), allocatable, save, target :: ctem_tile_yr


contains

!=================================================================================

subroutine alloc_ctem_vars()

use ctem_params, only : ican, icc, iccp1,ilg,nlat,nmos,ignd

implicit none

!-----------

! allocated with nlat,nmos, icc:

allocate(vrot%ailcmin(nlat,nmos,icc))
allocate(vrot%pftexist(nlat,nmos,icc))
allocate(vrot%lfstatus(nlat,nmos,icc))
allocate(vrot%pandays (nlat,nmos,icc))
allocate(vrot%ailcmin (nlat,nmos,icc))
allocate(vrot%ailcmax (nlat,nmos,icc))
allocate(vrot%dvdfcan (nlat,nmos,icc))
allocate(vrot%gleafmas(nlat,nmos,icc))
allocate(vrot%bleafmas(nlat,nmos,icc))
allocate(vrot%stemmass(nlat,nmos,icc))
allocate(vrot%rootmass(nlat,nmos,icc))
allocate(vrot%pstemmass (nlat,nmos,icc))
allocate(vrot%pgleafmass (nlat,nmos,icc))
allocate(vrot%fcancmx (nlat,nmos,icc))
allocate(vrot%ailcg   (nlat,nmos,icc))
allocate(vrot%ailcgs  (nlat,nmos,icc))
allocate(vrot%fcancs  (nlat,nmos,icc))
allocate(vrot%fcanc   (nlat,nmos,icc))
allocate(vrot%co2i1cg (nlat,nmos,icc))
allocate(vrot%co2i1cs (nlat,nmos,icc))
allocate(vrot%co2i2cg (nlat,nmos,icc))
allocate(vrot%co2i2cs (nlat,nmos,icc))
allocate(vrot%ancsveg (nlat,nmos,icc))
allocate(vrot%ancgveg (nlat,nmos,icc))
allocate(vrot%rmlcsveg(nlat,nmos,icc))
allocate(vrot%rmlcgveg(nlat,nmos,icc))
allocate(vrot%slai    (nlat,nmos,icc))
allocate(vrot%ailcb   (nlat,nmos,icc))
allocate(vrot%flhrloss(nlat,nmos,icc))
allocate(vrot%grwtheff(nlat,nmos,icc))
allocate(vrot%lystmmas(nlat,nmos,icc))
allocate(vrot%lyrotmas(nlat,nmos,icc))
allocate(vrot%tymaxlai(nlat,nmos,icc))
allocate(vrot%stmhrlos(nlat,nmos,icc))
allocate(vrot%vgbiomas_veg(nlat,nmos,icc))
allocate(vrot%emit_co2(nlat,nmos,icc))
allocate(vrot%emit_co (nlat,nmos,icc))
allocate(vrot%emit_ch4(nlat,nmos,icc))
allocate(vrot%emit_nmhc(nlat,nmos,icc))
allocate(vrot%emit_h2 (nlat,nmos,icc))
allocate(vrot%emit_nox(nlat,nmos,icc))
allocate(vrot%emit_n2o(nlat,nmos,icc))
allocate(vrot%emit_pm25(nlat,nmos,icc))
allocate(vrot%emit_tpm(nlat,nmos,icc))
allocate(vrot%emit_tc (nlat,nmos,icc))
allocate(vrot%emit_oc (nlat,nmos,icc))
allocate(vrot%emit_bc (nlat,nmos,icc))
allocate(vrot%burnvegf(nlat,nmos,icc))
allocate(vrot%smfuncveg(nlat,nmos,icc))
allocate(vrot%bterm   (nlat,nmos,icc))
allocate(vrot%mterm   (nlat,nmos,icc))
allocate(vrot%bmasveg (nlat,nmos,icc))
allocate(vrot%veghght (nlat,nmos,icc))
allocate(vrot%rootdpth(nlat,nmos,icc))
allocate(vrot%tltrleaf(nlat,nmos,icc))
allocate(vrot%tltrstem(nlat,nmos,icc))
allocate(vrot%tltrroot(nlat,nmos,icc))
allocate(vrot%leaflitr(nlat,nmos,icc))
allocate(vrot%roottemp(nlat,nmos,icc))
allocate(vrot%afrleaf (nlat,nmos,icc))
allocate(vrot%afrstem (nlat,nmos,icc))
allocate(vrot%afrroot (nlat,nmos,icc))
allocate(vrot%wtstatus(nlat,nmos,icc))
allocate(vrot%ltstatus(nlat,nmos,icc))
allocate(vrot%gppveg  (nlat,nmos,icc))
allocate(vrot%nppveg  (nlat,nmos,icc))
allocate(vrot%autoresveg(nlat,nmos,icc))
allocate(vrot%rmlvegacc (nlat,nmos,icc))
allocate(vrot%rmsveg  (nlat,nmos,icc))
allocate(vrot%rmrveg  (nlat,nmos,icc))
allocate(vrot%rgveg   (nlat,nmos,icc))
allocate(vrot%litrfallveg(nlat,nmos,icc))
allocate(vrot%rothrlos(nlat,nmos,icc))
allocate(vrot%pfcancmx(nlat,nmos,icc))
allocate(vrot%nfcancmx(nlat,nmos,icc))
allocate(vrot%anveg   (nlat,nmos,icc))
allocate(vrot%rmlveg  (nlat,nmos,icc))
    
! allocated with nlat,nmos:   
allocate(vrot%icount(nlat,nmos))
allocate(vrot%stdaln(nlat,nmos))
allocate(vrot%gavglai (nlat,nmos))
allocate(vrot%co2conc (nlat,nmos))
allocate(vrot%ch4conc (nlat,nmos))
allocate(vrot%canres  (nlat,nmos))
allocate(vrot%vgbiomas(nlat,nmos))
allocate(vrot%gavgltms(nlat,nmos))
allocate(vrot%gavgscms(nlat,nmos))
allocate(vrot%burnfrac(nlat,nmos))
allocate(vrot%popdin  (nlat,nmos))
allocate(vrot%lterm   (nlat,nmos))
allocate(vrot%extnprob(nlat,nmos))
allocate(vrot%prbfrhuc(nlat,nmos))
allocate(vrot%rml     (nlat,nmos))
allocate(vrot%rms     (nlat,nmos))
allocate(vrot%rmr     (nlat,nmos))
allocate(vrot%ch4wet1 (nlat,nmos))
allocate(vrot%ch4wet2 (nlat,nmos))
allocate(vrot%wetfdyn (nlat,nmos))
allocate(vrot%ch4dyn1 (nlat,nmos))
allocate(vrot%ch4dyn2 (nlat,nmos))
allocate(vrot%ch4_soills(nlat,nmos))
allocate(vrot%lucemcom(nlat,nmos))
allocate(vrot%lucltrin(nlat,nmos))
allocate(vrot%lucsocin(nlat,nmos))
allocate(vrot%npp     (nlat,nmos))
allocate(vrot%nep     (nlat,nmos))
allocate(vrot%nbp     (nlat,nmos))
allocate(vrot%gpp     (nlat,nmos))
allocate(vrot%hetrores(nlat,nmos))
allocate(vrot%autores (nlat,nmos))
allocate(vrot%soilcresp(nlat,nmos))
allocate(vrot%rm      (nlat,nmos))
allocate(vrot%rg      (nlat,nmos))
allocate(vrot%litres  (nlat,nmos))
allocate(vrot%socres  (nlat,nmos))
allocate(vrot%dstcemls(nlat,nmos))
allocate(vrot%litrfall(nlat,nmos))
allocate(vrot%humiftrs(nlat,nmos))
allocate(vrot%cfluxcg (nlat,nmos))
allocate(vrot%cfluxcs (nlat,nmos))
allocate(vrot%dstcemls3 (nlat,nmos))
allocate(vrot%PREACC_M(nlat,nmos))
allocate(vrot%GTACC_M (nlat,nmos))
allocate(vrot%QEVPACC_M (nlat,nmos))
allocate(vrot%HFSACC_M(nlat,nmos))
allocate(vrot%HMFNACC_M (nlat,nmos))
allocate(vrot%ROFACC_M(nlat,nmos))
allocate(vrot%SNOACC_M(nlat,nmos))
allocate(vrot%OVRACC_M(nlat,nmos))
allocate(vrot%WTBLACC_M(nlat,nmos))
allocate(vrot%ALVSACC_M(nlat,nmos))
allocate(vrot%ALIRACC_M(nlat,nmos))
allocate(vrot%RHOSACC_M(nlat,nmos))
allocate(vrot%TSNOACC_M(nlat,nmos))
allocate(vrot%WSNOACC_M(nlat,nmos))
allocate(vrot%SNOARE_M(nlat,nmos))
allocate(vrot%TCANACC_M(nlat,nmos))
allocate(vrot%RCANACC_M(nlat,nmos))
allocate(vrot%SCANACC_M(nlat,nmos))
allocate(vrot%ALTOTACC_M(nlat,nmos))
allocate(vrot%GROACC_M(nlat,nmos))
allocate(vrot%FSINACC_M (nlat,nmos))
allocate(vrot%FLINACC_M(nlat,nmos))
allocate(vrot%TAACC_M (nlat,nmos))
allocate(vrot%UVACC_M (nlat,nmos))
allocate(vrot%PRESACC_M (nlat,nmos))
allocate(vrot%QAACC_M (nlat,nmos))
allocate(vrot%EVAPACC_M (nlat,nmos))
allocate(vrot%FLUTACC_M(nlat,nmos))
allocate(vrot%tcanrs  (nlat,nmos))
allocate(vrot%tsnors  (nlat,nmos))
allocate(vrot%tpndrs  (nlat,nmos))
allocate(vrot%tcanoaccrow_m(nlat,nmos))
allocate(vrot%uvaccrow_m(nlat,nmos))
allocate(vrot%vvaccrow_m (nlat,nmos))
allocate(vrot%tcanoaccrow_out (nlat,nmos))
allocate(vrot%qevpacc_m_save(nlat,nmos))
allocate(vrot%twarmm  (nlat,nmos))
allocate(vrot%tcoldm  (nlat,nmos))
allocate(vrot%gdd5    (nlat,nmos))
allocate(vrot%aridity (nlat,nmos))
allocate(vrot%srplsmon(nlat,nmos))
allocate(vrot%defctmon(nlat,nmos))
allocate(vrot%anndefct(nlat,nmos))
allocate(vrot%annsrpls(nlat,nmos))
allocate(vrot%annpcp  (nlat,nmos))
allocate(vrot%dry_season_length(nlat,nmos))

! allocated with nlat,nmos,ican:     
allocate(vrot%zolnc(nlat,nmos,ican))
allocate(vrot%ailc(nlat,nmos,ican))
allocate(vrot%cmasvegc(nlat,nmos,ican))
allocate(vrot%alvsctm(nlat,nmos,ican))
allocate(vrot%paic(nlat,nmos,ican))
allocate(vrot%slaic(nlat,nmos,ican))
allocate(vrot%alirctm(nlat,nmos,ican))
allocate(vrot%csum(nlat,nmos,ican))
    
! allocated with nlat,nmos,ignd:        
allocate(vrot%TBARACC_M(nlat,nmos,ignd))
allocate(vrot%THLQACC_M(nlat,nmos,ignd))
allocate(vrot%THICACC_M(nlat,nmos,ignd))
allocate(vrot%THALACC_M(nlat,nmos,ignd))
allocate(vrot%tbaraccrow_m(nlat,nmos,ignd))
    
! allocated with nlat,nmos,ican,ignd:       
allocate(vrot%rmatc(nlat,nmos,ican,ignd))
 
 ! allocated with nlat,nmos,icc,ignd: 
allocate(vrot%rmatctem(nlat,nmos,icc,ignd))
    
! allocated with nlat,nmos,iccp1:     
allocate(vrot%litrmass(nlat,nmos,iccp1))
allocate(vrot%soilcmas(nlat,nmos,iccp1))
allocate(vrot%nepveg(nlat,nmos,iccp1))
allocate(vrot%nbpveg(nlat,nmos,iccp1))
allocate(vrot%hetroresveg(nlat,nmos,iccp1))
allocate(vrot%litresveg(nlat,nmos,iccp1))
allocate(vrot%soilcresveg(nlat,nmos,iccp1))
allocate(vrot%humiftrsveg(nlat,nmos,iccp1))
    
! allocated with nlat,nmos,{some number}: 
allocate(vrot%colddays(nlat,nmos,2))
allocate(vrot%slopefrac(nlat,nmos,8))
allocate(vrot%mlightng(nlat,nmos,12))
allocate(vrot%wetfrac_mon(nlat,nmos,12))
    
! allocated with nlat:     
allocate(vrot%dayl_max(nlat))
allocate(vrot%dayl(nlat))
allocate(vrot%altotcntr_d(nlat))

! Now on to the veg_gat vars

ilg = nlat * nmos

! allocated with ilg

allocate(vgat%icount(ilg))
allocate(vgat%gavglai (ilg))
allocate(vgat%lightng (ilg))
allocate(vgat%tcanoaccgat_out (ilg))
allocate(vgat%co2conc (ilg))
allocate(vgat%ch4conc (ilg))
allocate(vgat%canres (ilg))
allocate(vgat%vgbiomas (ilg))
allocate(vgat%gavgltms (ilg))
allocate(vgat%gavgscms (ilg))
allocate(vgat%lterm (ilg))
allocate(vgat%ch4wet1 (ilg))
allocate(vgat%ch4wet2 (ilg))
allocate(vgat%wetfdyn (ilg))
allocate(vgat%ch4dyn1 (ilg))
allocate(vgat%ch4dyn2 (ilg))
allocate(vgat%ch4_soills (ilg))
allocate(vgat%lucemcom (ilg))
allocate(vgat%lucltrin (ilg))
allocate(vgat%lucsocin (ilg))
allocate(vgat%npp (ilg))
allocate(vgat%nep (ilg))
allocate(vgat%nbp (ilg))
allocate(vgat%gpp (ilg))
allocate(vgat%hetrores (ilg))
allocate(vgat%autores (ilg))
allocate(vgat%soilcresp (ilg))
allocate(vgat%rm (ilg))
allocate(vgat%rg (ilg))
allocate(vgat%litres (ilg))
allocate(vgat%socres (ilg))
allocate(vgat%dstcemls (ilg))
allocate(vgat%litrfall (ilg))
allocate(vgat%humiftrs (ilg))
allocate(vgat%rml (ilg))
allocate(vgat%rms (ilg))
allocate(vgat%rmr (ilg))
allocate(vgat%burnfrac (ilg))
allocate(vgat%popdin (ilg))
allocate(vgat%extnprob (ilg))
allocate(vgat%prbfrhuc (ilg))
allocate(vgat%dayl_max (ilg))
allocate(vgat%dayl (ilg))
allocate(vgat%wetfrac_pres (ilg))
allocate(vgat%twarmm (ilg))
allocate(vgat%tcoldm (ilg))
allocate(vgat%gdd5 (ilg))
allocate(vgat%aridity (ilg))
allocate(vgat%srplsmon (ilg))
allocate(vgat%defctmon (ilg))
allocate(vgat%anndefct (ilg))
allocate(vgat%annsrpls (ilg))
allocate(vgat%annpcp (ilg))
allocate(vgat%dry_season_length (ilg))
allocate(vgat%cfluxcg (ilg))
allocate(vgat%cfluxcs (ilg))
allocate(vgat%dstcemls3 (ilg))
allocate(vgat%tcurm (ilg))
allocate(vgat%srpcuryr (ilg))
allocate(vgat%dftcuryr (ilg))
allocate(vgat%anpcpcur (ilg))
allocate(vgat%anpecur (ilg))
allocate(vgat%gdd5cur (ilg))
allocate(vgat%surmncur (ilg))
allocate(vgat%defmncur (ilg))
allocate(vgat%srplscur (ilg))
allocate(vgat%defctcur (ilg))
allocate(vgat%stdaln (ilg))

! allocated with ilg, icc
allocate(vgat%ailcmin (ilg,icc))
allocate(vgat%ailcmax (ilg,icc))
allocate(vgat%dvdfcan (ilg,icc))
allocate(vgat%gleafmas (ilg,icc))
allocate(vgat%bleafmas (ilg,icc))
allocate(vgat%stemmass (ilg,icc))
allocate(vgat%rootmass (ilg,icc))
allocate(vgat%pstemmass (ilg,icc))
allocate(vgat%pgleafmass (ilg,icc))
allocate(vgat%fcancmx (ilg,icc))
allocate(vgat%ailcg (ilg,icc))
allocate(vgat%ailcgs (ilg,icc))
allocate(vgat%fcancs (ilg,icc))
allocate(vgat%fcanc (ilg,icc))
allocate(vgat%co2i1cg (ilg,icc))
allocate(vgat%co2i1cs (ilg,icc))
allocate(vgat%co2i2cg (ilg,icc))
allocate(vgat%co2i2cs (ilg,icc))
allocate(vgat%ancsveg (ilg,icc))
allocate(vgat%ancgveg (ilg,icc))
allocate(vgat%rmlcsveg (ilg,icc))
allocate(vgat%rmlcgveg (ilg,icc))
allocate(vgat%slai (ilg,icc))
allocate(vgat%ailcb (ilg,icc))
allocate(vgat%flhrloss (ilg,icc))
allocate(vgat%grwtheff (ilg,icc))
allocate(vgat%lystmmas (ilg,icc))
allocate(vgat%lyrotmas (ilg,icc))
allocate(vgat%tymaxlai (ilg,icc))
allocate(vgat%stmhrlos (ilg,icc))
allocate(vgat%vgbiomas_veg (ilg,icc))
allocate(vgat%emit_co2 (ilg,icc))
allocate(vgat%emit_co (ilg,icc))
allocate(vgat%emit_ch4 (ilg,icc))
allocate(vgat%emit_nmhc (ilg,icc))
allocate(vgat%emit_h2 (ilg,icc))
allocate(vgat%emit_nox (ilg,icc))
allocate(vgat%emit_n2o (ilg,icc))
allocate(vgat%emit_pm25 (ilg,icc))
allocate(vgat%emit_tpm (ilg,icc))
allocate(vgat%emit_tc (ilg,icc))
allocate(vgat%emit_oc (ilg,icc))
allocate(vgat%emit_bc (ilg,icc))
allocate(vgat%burnvegf (ilg,icc))
allocate(vgat%smfuncveg (ilg,icc))
allocate(vgat%bterm (ilg,icc))
allocate(vgat%mterm (ilg,icc))
allocate(vgat%bmasveg (ilg,icc))
allocate(vgat%veghght (ilg,icc))
allocate(vgat%rootdpth (ilg,icc))
allocate(vgat%tltrleaf (ilg,icc))
allocate(vgat%tltrstem (ilg,icc))
allocate(vgat%tltrroot (ilg,icc))
allocate(vgat%leaflitr (ilg,icc))
allocate(vgat%roottemp (ilg,icc))
allocate(vgat%afrleaf (ilg,icc))
allocate(vgat%afrstem (ilg,icc))
allocate(vgat%afrroot (ilg,icc))
allocate(vgat%wtstatus (ilg,icc))
allocate(vgat%ltstatus (ilg,icc))
allocate(vgat%gppveg (ilg,icc))
allocate(vgat%nppveg (ilg,icc))
allocate(vgat%autoresveg (ilg,icc))
allocate(vgat%rmlvegacc (ilg,icc))
allocate(vgat%rmsveg (ilg,icc))
allocate(vgat%rmrveg (ilg,icc))
allocate(vgat%rgveg (ilg,icc))
allocate(vgat%litrfallveg (ilg,icc))
allocate(vgat%anveg (ilg,icc))
allocate(vgat%rmlveg (ilg,icc))
allocate(vgat%geremort (ilg,icc))
allocate(vgat%intrmort (ilg,icc))
allocate(vgat%lambda (ilg,icc))
allocate(vgat%cc (ilg,icc))
allocate(vgat%mm (ilg,icc))
allocate(vgat%pftexist (ilg,icc))
allocate(vgat%rothrlos (ilg,icc))
allocate(vgat%pfcancmx (ilg,icc))
allocate(vgat%nfcancmx (ilg,icc))
allocate(vgat%lfstatus (ilg,icc))
allocate(vgat%pandays (ilg,icc))

! allocated with ilg, ican:
allocate(vgat%zolnc (ilg,ican))
allocate(vgat%ailc (ilg,ican))
allocate(vgat%cmasvegc (ilg,ican))
allocate(vgat%alvsctm (ilg,ican))
allocate(vgat%paic (ilg,ican))
allocate(vgat%slaic (ilg,ican))
allocate(vgat%alirctm (ilg,ican))

! allocated with ilg, iccp1:
allocate(vgat%litrmass (ilg,iccp1))
allocate(vgat%soilcmas (ilg,iccp1))
allocate(vgat%nepveg (ilg,iccp1))
allocate(vgat%nbpveg (ilg,iccp1))
allocate(vgat%hetroresveg (ilg,iccp1))
allocate(vgat%litresveg (ilg,iccp1))
allocate(vgat%soilcresveg (ilg,iccp1))
allocate(vgat%humiftrsveg (ilg,iccp1))

! allocated with ilg, ican, ignd:
allocate(vgat%rmatc (ilg,ican,ignd))

! allocated with ilg, icc, ignd:
allocate(vgat%rmatctem (ilg,icc,ignd))

! allocated with ilg, icc, {some number}:
allocate(vgat%colddays (ilg,2))
allocate(vgat%slopefrac (ilg,8))
allocate(vgat%mlightng (ilg,12))
allocate(vgat%wetfrac_mon (ilg,12))
allocate(vgat%tmonth (12,ilg))

! allocated with nlat:
allocate(class_out%ALVSACC_MO(nlat))
allocate(class_out%ALIRACC_MO (nlat))
allocate(class_out%FLUTACC_MO (nlat))
allocate(class_out%FSINACC_MO (nlat))
allocate(class_out%FLINACC_MO (nlat))
allocate(class_out%HFSACC_MO (nlat))
allocate(class_out%QEVPACC_MO (nlat))
allocate(class_out%SNOACC_MO (nlat))
allocate(class_out%WSNOACC_MO (nlat))
allocate(class_out%ROFACC_MO (nlat))
allocate(class_out%PREACC_MO (nlat))
allocate(class_out%EVAPACC_MO (nlat))
allocate(class_out%TRANSPACC_MO (nlat))
allocate(class_out%TAACC_MO (nlat))
allocate(class_out%ACTLYR_MO (nlat))
allocate(class_out%FTABLE_MO (nlat))
allocate(class_out%ACTLYR_MIN_MO (nlat))
allocate(class_out%ACTLYR_MAX_MO (nlat))
allocate(class_out%FTABLE_MIN_MO (nlat))
allocate(class_out%FTABLE_MAX_MO (nlat))
allocate(class_out%ALTOTACC_MO (nlat))
allocate(class_out%GROUNDEVAP (nlat))
allocate(class_out%CANOPYEVAP (nlat))
allocate(class_out%altotcntr_m (nlat))

! allocated with nlat,ignd:
allocate(class_out%TBARACC_MO (nlat,ignd))
allocate(class_out%THLQACC_MO (nlat,ignd))
allocate(class_out%THICACC_MO (nlat,ignd))
allocate(class_out%ALVSACC_YR (nlat))
allocate(class_out%ALIRACC_YR (nlat))
allocate(class_out%FLUTACC_YR (nlat))
allocate(class_out%FSINACC_YR (nlat))
allocate(class_out%FLINACC_YR (nlat))
allocate(class_out%HFSACC_YR (nlat))
allocate(class_out%QEVPACC_YR (nlat))
allocate(class_out%ROFACC_YR (nlat))
allocate(class_out%PREACC_YR (nlat))
allocate(class_out%EVAPACC_YR (nlat))
allocate(class_out%TRANSPACC_YR (nlat))
allocate(class_out%TAACC_YR (nlat))
allocate(class_out%ACTLYR_YR (nlat))
allocate(class_out%ACTLYR_MIN_YR (nlat))
allocate(class_out%ACTLYR_MAX_YR (nlat))
allocate(class_out%FTABLE_YR (nlat))
allocate(class_out%FTABLE_MIN_YR (nlat))
allocate(class_out%FTABLE_MAX_YR (nlat))
allocate(class_out%ALTOTACC_YR (nlat))
allocate(class_out%altotcntr_yr (nlat))

allocate(ctem_grd%WSNOROT_g (nlat))
allocate(ctem_grd%ROFSROT_g (nlat))
allocate(ctem_grd%SNOROT_g (nlat))
allocate(ctem_grd%RHOSROT_g (nlat))
allocate(ctem_grd%ROFROT_g (nlat))
allocate(ctem_grd%ZPNDROT_g (nlat))
allocate(ctem_grd%RCANROT_g (nlat))
allocate(ctem_grd%SCANROT_g (nlat))
allocate(ctem_grd%TROFROT_g (nlat))
allocate(ctem_grd%TROOROT_g (nlat))
allocate(ctem_grd%TROBROT_g (nlat))
allocate(ctem_grd%ROFOROT_g (nlat))
allocate(ctem_grd%ROFBROT_g (nlat))
allocate(ctem_grd%TROSROT_g (nlat))
allocate(ctem_grd%FSGVROT_g (nlat))
allocate(ctem_grd%FSGSROT_g (nlat))
allocate(ctem_grd%FLGVROT_g (nlat))
allocate(ctem_grd%FLGSROT_g (nlat))
allocate(ctem_grd%HFSCROT_g (nlat))
allocate(ctem_grd%HFSSROT_g (nlat))
allocate(ctem_grd%HEVCROT_g (nlat))
allocate(ctem_grd%HEVSROT_g (nlat))
allocate(ctem_grd%HMFCROT_g (nlat))
allocate(ctem_grd%HMFNROT_g (nlat))
allocate(ctem_grd%HTCSROT_g (nlat))
allocate(ctem_grd%HTCCROT_g (nlat))
allocate(ctem_grd%FSGGROT_g (nlat))
allocate(ctem_grd%FLGGROT_g (nlat))
allocate(ctem_grd%HFSGROT_g (nlat))
allocate(ctem_grd%HEVGROT_g (nlat))
allocate(ctem_grd%CDHROT_g (nlat))
allocate(ctem_grd%CDMROT_g (nlat))
allocate(ctem_grd%SFCUROT_g (nlat))
allocate(ctem_grd%SFCVROT_g (nlat))
allocate(ctem_grd%ACTLYR_g (nlat))
allocate(ctem_grd%FTABLE_g (nlat))
allocate(ctem_grd%fc_g (nlat))
allocate(ctem_grd%fg_g (nlat))
allocate(ctem_grd%fcs_g (nlat))
allocate(ctem_grd%fgs_g (nlat))
allocate(ctem_grd%PCFCROT_g (nlat))
allocate(ctem_grd%PCLCROT_g (nlat))
allocate(ctem_grd%PCPGROT_g (nlat))
allocate(ctem_grd%QFCFROT_g (nlat))
allocate(ctem_grd%QFGROT_g (nlat))
allocate(ctem_grd%ROFCROT_g (nlat))
allocate(ctem_grd%QFCROT_g (nlat,ignd))
allocate(ctem_grd%ROFNROT_g (nlat))
allocate(ctem_grd%WTRSROT_g (nlat))
allocate(ctem_grd%WTRGROT_g (nlat))
allocate(ctem_grd%PCPNROT_g (nlat))
allocate(ctem_grd%QFCLROT_g (nlat))
allocate(ctem_grd%QFNROT_g (nlat))
allocate(ctem_grd%WTRCROT_g (nlat))
allocate(ctem_grd%gpp_g (nlat))
allocate(ctem_grd%npp_g (nlat))
allocate(ctem_grd%nbp_g (nlat))
allocate(ctem_grd%socres_g (nlat))
allocate(ctem_grd%autores_g (nlat))
allocate(ctem_grd%litres_g (nlat))
allocate(ctem_grd%dstcemls3_g (nlat))
allocate(ctem_grd%litrfall_g (nlat))
allocate(ctem_grd%rml_g (nlat))
allocate(ctem_grd%rms_g (nlat))
allocate(ctem_grd%rg_g (nlat))
allocate(ctem_grd%leaflitr_g (nlat))
allocate(ctem_grd%tltrstem_g (nlat))
allocate(ctem_grd%tltrroot_g (nlat))
allocate(ctem_grd%nep_g (nlat))
allocate(ctem_grd%hetrores_g (nlat))
allocate(ctem_grd%dstcemls_g (nlat))
allocate(ctem_grd%humiftrs_g (nlat))
allocate(ctem_grd%rmr_g (nlat))
allocate(ctem_grd%tltrleaf_g (nlat))
allocate(ctem_grd%gavgltms_g (nlat))
allocate(ctem_grd%vgbiomas_g (nlat))
allocate(ctem_grd%gavglai_g (nlat))
allocate(ctem_grd%gavgscms_g (nlat))
allocate(ctem_grd%gleafmas_g (nlat))
allocate(ctem_grd%bleafmas_g (nlat))
allocate(ctem_grd%stemmass_g (nlat))
allocate(ctem_grd%rootmass_g (nlat))
allocate(ctem_grd%litrmass_g (nlat))
allocate(ctem_grd%soilcmas_g (nlat))
allocate(ctem_grd%slai_g (nlat))
allocate(ctem_grd%ailcg_g (nlat))
allocate(ctem_grd%ailcb_g (nlat))
allocate(ctem_grd%veghght_g (nlat))
allocate(ctem_grd%rootdpth_g (nlat))
allocate(ctem_grd%roottemp_g (nlat))
allocate(ctem_grd%totcmass_g (nlat))
allocate(ctem_grd%tcanoacc_out_g (nlat))
allocate(ctem_grd%burnfrac_g (nlat))
allocate(ctem_grd%smfuncveg_g (nlat))
allocate(ctem_grd%lucemcom_g (nlat))
allocate(ctem_grd%lucltrin_g (nlat))
allocate(ctem_grd%lucsocin_g (nlat))
allocate(ctem_grd%emit_co2_g (nlat))
allocate(ctem_grd%emit_co_g (nlat))
allocate(ctem_grd%emit_ch4_g (nlat))
allocate(ctem_grd%emit_nmhc_g (nlat))
allocate(ctem_grd%emit_h2_g (nlat))
allocate(ctem_grd%emit_nox_g (nlat))
allocate(ctem_grd%emit_n2o_g (nlat))
allocate(ctem_grd%emit_pm25_g (nlat))
allocate(ctem_grd%emit_tpm_g (nlat))
allocate(ctem_grd%emit_tc_g (nlat))
allocate(ctem_grd%emit_oc_g (nlat))
allocate(ctem_grd%emit_bc_g (nlat))
allocate(ctem_grd%bterm_g (nlat))
allocate(ctem_grd%lterm_g (nlat))
allocate(ctem_grd%mterm_g (nlat))
allocate(ctem_grd%ch4wet1_g (nlat))
allocate(ctem_grd%ch4wet2_g (nlat))
allocate(ctem_grd%wetfdyn_g (nlat))
allocate(ctem_grd%ch4dyn1_g (nlat))
allocate(ctem_grd%ch4dyn2_g (nlat))
allocate(ctem_grd%ch4_soills_g (nlat))
allocate(ctem_grd%afrleaf_g (nlat,icc))
allocate(ctem_grd%afrstem_g (nlat,icc))
allocate(ctem_grd%afrroot_g (nlat,icc))
allocate(ctem_grd%lfstatus_g (nlat,icc))
allocate(ctem_grd%rmlvegrow_g (nlat,icc))
allocate(ctem_grd%anvegrow_g(nlat,icc))
allocate(ctem_grd%rmatctem_g (nlat,ignd))
allocate(ctem_grd%HMFGROT_g (nlat,ignd))
allocate(ctem_grd%HTCROT_g(nlat,ignd))
allocate(ctem_grd%TBARROT_g (nlat,ignd))
allocate(ctem_grd%THLQROT_g (nlat,ignd))
allocate(ctem_grd%THICROT_g (nlat,ignd))
allocate(ctem_grd%GFLXROT_g (nlat,ignd))
allocate(ctem_grd%fsstar_g (nlat))
allocate(ctem_grd%flstar_g (nlat))
allocate(ctem_grd%qh_g(nlat))
allocate(ctem_grd%qe_g (nlat))
allocate(ctem_grd%snomlt_g (nlat))
allocate(ctem_grd%beg_g (nlat))
allocate(ctem_grd%gtout_g (nlat))
allocate(ctem_grd%tpn_g (nlat))
allocate(ctem_grd%altot_g(nlat))
allocate(ctem_grd%tcn_g (nlat))
allocate(ctem_grd%tsn_g (nlat))
allocate(ctem_grd%zsn_g (nlat))

allocate(ctem_tile%leaflitr_t (nlat,nmos))
allocate(ctem_tile%tltrleaf_t (nlat,nmos))
allocate(ctem_tile%tltrstem_t (nlat,nmos))
allocate(ctem_tile%tltrroot_t (nlat,nmos))
allocate(ctem_tile%ailcg_t (nlat,nmos))
allocate(ctem_tile%ailcb_t (nlat,nmos))
allocate(ctem_tile%rmatctem_t (nlat,nmos,ignd))
allocate(ctem_tile%veghght_t (nlat,nmos))
allocate(ctem_tile%rootdpth_t (nlat,nmos))
allocate(ctem_tile%roottemp_t (nlat,nmos))
allocate(ctem_tile%slai_t (nlat,nmos))
allocate(ctem_tile%afrroot_t (nlat,nmos))
allocate(ctem_tile%afrleaf_t (nlat,nmos))
allocate(ctem_tile%afrstem_t (nlat,nmos))
allocate(ctem_tile%laimaxg_t (nlat,nmos))
allocate(ctem_tile%stemmass_t (nlat,nmos))
allocate(ctem_tile%rootmass_t (nlat,nmos))
allocate(ctem_tile%litrmass_t (nlat,nmos))
allocate(ctem_tile%gleafmas_t (nlat,nmos))
allocate(ctem_tile%bleafmas_t(nlat,nmos))
allocate(ctem_tile%soilcmas_t (nlat,nmos))
allocate(ctem_tile%emit_co2_t (nlat,nmos))
allocate(ctem_tile%emit_co_t (nlat,nmos))
allocate(ctem_tile%emit_ch4_t (nlat,nmos))
allocate(ctem_tile%emit_nmhc_t (nlat,nmos))
allocate(ctem_tile%emit_h2_t (nlat,nmos))
allocate(ctem_tile%emit_nox_t (nlat,nmos))
allocate(ctem_tile%emit_n2o_t (nlat,nmos))
allocate(ctem_tile%emit_pm25_t (nlat,nmos))
allocate(ctem_tile%emit_tpm_t (nlat,nmos))
allocate(ctem_tile%emit_tc_t (nlat,nmos))
allocate(ctem_tile%emit_oc_t (nlat,nmos))
allocate(ctem_tile%emit_bc_t (nlat,nmos))
allocate(ctem_tile%bterm_t (nlat,nmos))
allocate(ctem_tile%mterm_t (nlat,nmos))
allocate(ctem_tile%smfuncveg_t (nlat,nmos))
allocate(ctem_tile%fsnowacc_t (ilg))
allocate(ctem_tile%tcansacc_t (ilg))
allocate(ctem_tile%tcanoaccgat_t (ilg))
allocate(ctem_tile%taaccgat_t (ilg))
allocate(ctem_tile%uvaccgat_t (ilg))
allocate(ctem_tile%vvaccgat_t (ilg))
allocate(ctem_tile%tbaraccgat_t (ilg,ignd))
allocate(ctem_tile%tbarcacc_t (ilg,ignd))
allocate(ctem_tile%tbarcsacc_t (ilg,ignd))
allocate(ctem_tile%tbargacc_t (ilg,ignd))
allocate(ctem_tile%tbargsacc_t (ilg,ignd))
allocate(ctem_tile%thliqcacc_t (ilg,ignd))
allocate(ctem_tile%thliqgacc_t (ilg,ignd))
allocate(ctem_tile%thliqacc_t (ilg,ignd))
allocate(ctem_tile%thicecacc_t (ilg,ignd))
allocate(ctem_tile%thicegacc_t (ilg,ignd))
allocate(ctem_tile%ancsvgac_t (ilg,icc))
allocate(ctem_tile%ancgvgac_t (ilg,icc))
allocate(ctem_tile%rmlcsvga_t (ilg,icc))
allocate(ctem_tile%rmlcgvga_t (ilg,icc))

allocate(ctem_mo%laimaxg_mo (nlat,nmos,icc))
allocate(ctem_mo%stemmass_mo (nlat,nmos,icc))
allocate(ctem_mo%rootmass_mo (nlat,nmos,icc))
allocate(ctem_mo%litrfallveg_mo (nlat,nmos,icc))
allocate(ctem_mo%npp_mo (nlat,nmos,icc))
allocate(ctem_mo%gpp_mo (nlat,nmos,icc))
allocate(ctem_mo%vgbiomas_mo (nlat,nmos,icc))
allocate(ctem_mo%autores_mo (nlat,nmos,icc))
allocate(ctem_mo%humiftrsveg_mo (nlat,nmos,iccp1))
allocate(ctem_mo%totcmass_mo (nlat,nmos,iccp1))
allocate(ctem_mo%litrmass_mo (nlat,nmos,iccp1))
allocate(ctem_mo%soilcmas_mo (nlat,nmos,iccp1))
allocate(ctem_mo%nep_mo (nlat,nmos,iccp1))
allocate(ctem_mo%litres_mo (nlat,nmos,iccp1))
allocate(ctem_mo%soilcres_mo (nlat,nmos,iccp1))
allocate(ctem_mo%hetrores_mo (nlat,nmos,iccp1))
allocate(ctem_mo%nbp_mo (nlat,nmos,iccp1))
allocate(ctem_mo%emit_co2_mo (nlat,nmos,icc))
allocate(ctem_mo%emit_co_mo (nlat,nmos,icc))
allocate(ctem_mo%emit_ch4_mo (nlat,nmos,icc))
allocate(ctem_mo%emit_nmhc_mo (nlat,nmos,icc))
allocate(ctem_mo%emit_h2_mo (nlat,nmos,icc))
allocate(ctem_mo%emit_nox_mo (nlat,nmos,icc))
allocate(ctem_mo%emit_n2o_mo (nlat,nmos,icc))
allocate(ctem_mo%emit_pm25_mo (nlat,nmos,icc))
allocate(ctem_mo%emit_tpm_mo (nlat,nmos,icc))
allocate(ctem_mo%emit_tc_mo (nlat,nmos,icc))
allocate(ctem_mo%emit_oc_mo (nlat,nmos,icc))
allocate(ctem_mo%emit_bc_mo (nlat,nmos,icc))
allocate(ctem_mo%burnfrac_mo (nlat,nmos,icc))
allocate(ctem_mo%bterm_mo (nlat,nmos,icc))
allocate(ctem_mo%mterm_mo (nlat,nmos,icc))
allocate(ctem_mo%smfuncveg_mo (nlat,nmos,icc))

allocate(ctem_grd_mo%laimaxg_mo_g (nlat))
allocate(ctem_grd_mo%stemmass_mo_g (nlat))
allocate(ctem_grd_mo%rootmass_mo_g (nlat))
allocate(ctem_grd_mo%litrmass_mo_g (nlat))
allocate(ctem_grd_mo%soilcmas_mo_g (nlat))
allocate(ctem_grd_mo%litrfall_mo_g (nlat))
allocate(ctem_grd_mo%humiftrs_mo_g (nlat))
allocate(ctem_grd_mo%npp_mo_g (nlat))
allocate(ctem_grd_mo%gpp_mo_g (nlat))
allocate(ctem_grd_mo%nep_mo_g (nlat))
allocate(ctem_grd_mo%nbp_mo_g (nlat))
allocate(ctem_grd_mo%hetrores_mo_g (nlat))
allocate(ctem_grd_mo%autores_mo_g (nlat))
allocate(ctem_grd_mo%litres_mo_g (nlat))
allocate(ctem_grd_mo%soilcres_mo_g (nlat))
allocate(ctem_grd_mo%vgbiomas_mo_g (nlat))
allocate(ctem_grd_mo%totcmass_mo_g (nlat))
allocate(ctem_grd_mo%emit_co2_mo_g (nlat))
allocate(ctem_grd_mo%emit_co_mo_g (nlat))
allocate(ctem_grd_mo%emit_ch4_mo_g (nlat))
allocate(ctem_grd_mo%emit_nmhc_mo_g (nlat))
allocate(ctem_grd_mo%emit_h2_mo_g (nlat))
allocate(ctem_grd_mo%emit_nox_mo_g (nlat))
allocate(ctem_grd_mo%emit_n2o_mo_g (nlat))
allocate(ctem_grd_mo%emit_pm25_mo_g (nlat))
allocate(ctem_grd_mo%emit_tpm_mo_g (nlat))
allocate(ctem_grd_mo%emit_tc_mo_g (nlat))
allocate(ctem_grd_mo%emit_oc_mo_g (nlat))
allocate(ctem_grd_mo%emit_bc_mo_g (nlat))
allocate(ctem_grd_mo%smfuncveg_mo_g (nlat))
allocate(ctem_grd_mo%luc_emc_mo_g (nlat))
allocate(ctem_grd_mo%lucltrin_mo_g (nlat))
allocate(ctem_grd_mo%lucsocin_mo_g (nlat))
allocate(ctem_grd_mo%burnfrac_mo_g (nlat))
allocate(ctem_grd_mo%bterm_mo_g (nlat))
allocate(ctem_grd_mo%lterm_mo_g(nlat))
allocate(ctem_grd_mo%mterm_mo_g (nlat))
allocate(ctem_grd_mo%ch4wet1_mo_g (nlat))
allocate(ctem_grd_mo%ch4wet2_mo_g (nlat))
allocate(ctem_grd_mo%wetfdyn_mo_g (nlat))
allocate(ctem_grd_mo%ch4dyn1_mo_g (nlat))
allocate(ctem_grd_mo%ch4dyn2_mo_g (nlat))
allocate(ctem_grd_mo%ch4soills_mo_g (nlat))

allocate(ctem_tile_mo%laimaxg_mo_t (nlat,nmos))
allocate(ctem_tile_mo%stemmass_mo_t (nlat,nmos))
allocate(ctem_tile_mo%rootmass_mo_t (nlat,nmos))
allocate(ctem_tile_mo%litrfall_mo_t (nlat,nmos))
allocate(ctem_tile_mo%humiftrs_mo_t (nlat,nmos))
allocate(ctem_tile_mo%npp_mo_t (nlat,nmos))
allocate(ctem_tile_mo%gpp_mo_t (nlat,nmos))
allocate(ctem_tile_mo%vgbiomas_mo_t (nlat,nmos))
allocate(ctem_tile_mo%autores_mo_t (nlat,nmos))
allocate(ctem_tile_mo%totcmass_mo_t (nlat,nmos))
allocate(ctem_tile_mo%litrmass_mo_t (nlat,nmos))
allocate(ctem_tile_mo%soilcmas_mo_t (nlat,nmos))
allocate(ctem_tile_mo%nep_mo_t (nlat,nmos))
allocate(ctem_tile_mo%litres_mo_t (nlat,nmos))
allocate(ctem_tile_mo%soilcres_mo_t (nlat,nmos))
allocate(ctem_tile_mo%hetrores_mo_t (nlat,nmos))
allocate(ctem_tile_mo%nbp_mo_t (nlat,nmos))
allocate(ctem_tile_mo%emit_co2_mo_t (nlat,nmos))
allocate(ctem_tile_mo%emit_co_mo_t (nlat,nmos))
allocate(ctem_tile_mo%emit_ch4_mo_t (nlat,nmos))
allocate(ctem_tile_mo%emit_nmhc_mo_t (nlat,nmos))
allocate(ctem_tile_mo%emit_h2_mo_t (nlat,nmos))
allocate(ctem_tile_mo%emit_nox_mo_t (nlat,nmos))
allocate(ctem_tile_mo%emit_n2o_mo_t (nlat,nmos))
allocate(ctem_tile_mo%emit_pm25_mo_t (nlat,nmos))
allocate(ctem_tile_mo%emit_tpm_mo_t (nlat,nmos))
allocate(ctem_tile_mo%emit_tc_mo_t (nlat,nmos))
allocate(ctem_tile_mo%emit_oc_mo_t (nlat,nmos))
allocate(ctem_tile_mo%emit_bc_mo_t (nlat,nmos))
allocate(ctem_tile_mo%burnfrac_mo_t (nlat,nmos))
allocate(ctem_tile_mo%smfuncveg_mo_t (nlat,nmos))
allocate(ctem_tile_mo%bterm_mo_t (nlat,nmos))
allocate(ctem_tile_mo%luc_emc_mo_t (nlat,nmos))
allocate(ctem_tile_mo%lterm_mo_t (nlat,nmos))
allocate(ctem_tile_mo%lucsocin_mo_t (nlat,nmos))
allocate(ctem_tile_mo%mterm_mo_t (nlat,nmos))
allocate(ctem_tile_mo%lucltrin_mo_t (nlat,nmos))
allocate(ctem_tile_mo%ch4wet1_mo_t (nlat,nmos))
allocate(ctem_tile_mo%ch4wet2_mo_t (nlat,nmos))
allocate(ctem_tile_mo%wetfdyn_mo_t (nlat,nmos))
allocate(ctem_tile_mo%ch4dyn1_mo_t (nlat,nmos))
allocate(ctem_tile_mo%ch4dyn2_mo_t (nlat,nmos))
allocate(ctem_tile_mo%ch4soills_mo_t (nlat,nmos))
allocate(ctem_tile_mo%wind_mo_t (nlat,nmos))

allocate(ctem_yr%laimaxg_yr (nlat,nmos,icc))
allocate(ctem_yr%stemmass_yr (nlat,nmos,icc))
allocate(ctem_yr%rootmass_yr (nlat,nmos,icc))
allocate(ctem_yr%npp_yr (nlat,nmos,icc))
allocate(ctem_yr%gpp_yr (nlat,nmos,icc))
allocate(ctem_yr%vgbiomas_yr (nlat,nmos,icc))
allocate(ctem_yr%autores_yr (nlat,nmos,icc))
allocate(ctem_yr%totcmass_yr (nlat,nmos,iccp1))
allocate(ctem_yr%litrmass_yr (nlat,nmos,iccp1))
allocate(ctem_yr%soilcmas_yr (nlat,nmos,iccp1))
allocate(ctem_yr%nep_yr (nlat,nmos,iccp1))
allocate(ctem_yr%litres_yr (nlat,nmos,iccp1))
allocate(ctem_yr%soilcres_yr (nlat,nmos,iccp1))
allocate(ctem_yr%hetrores_yr (nlat,nmos,iccp1))
allocate(ctem_yr%nbp_yr (nlat,nmos,iccp1))
allocate(ctem_yr%emit_co2_yr (nlat,nmos,icc))
allocate(ctem_yr%emit_co_yr (nlat,nmos,icc))
allocate(ctem_yr%emit_ch4_yr (nlat,nmos,icc))
allocate(ctem_yr%emit_nmhc_yr (nlat,nmos,icc))
allocate(ctem_yr%emit_h2_yr (nlat,nmos,icc))
allocate(ctem_yr%emit_nox_yr (nlat,nmos,icc))
allocate(ctem_yr%emit_n2o_yr (nlat,nmos,icc))
allocate(ctem_yr%emit_pm25_yr (nlat,nmos,icc))
allocate(ctem_yr%emit_tpm_yr (nlat,nmos,icc))
allocate(ctem_yr%emit_tc_yr (nlat,nmos,icc))
allocate(ctem_yr%emit_oc_yr (nlat,nmos,icc))
allocate(ctem_yr%emit_bc_yr (nlat,nmos,icc))
allocate(ctem_yr%bterm_yr (nlat,nmos,icc))
allocate(ctem_yr%mterm_yr (nlat,nmos,icc))
allocate(ctem_yr%burnfrac_yr (nlat,nmos,icc))
allocate(ctem_yr%smfuncveg_yr (nlat,nmos,icc))
allocate(ctem_yr%veghght_yr (nlat,nmos,icc))

allocate(ctem_grd_yr%laimaxg_yr_g (nlat))
allocate(ctem_grd_yr%stemmass_yr_g (nlat))
allocate(ctem_grd_yr%rootmass_yr_g (nlat))
allocate(ctem_grd_yr%litrmass_yr_g (nlat))
allocate(ctem_grd_yr%soilcmas_yr_g (nlat))
allocate(ctem_grd_yr%npp_yr_g (nlat))
allocate(ctem_grd_yr%gpp_yr_g (nlat))
allocate(ctem_grd_yr%nep_yr_g (nlat))
allocate(ctem_grd_yr%nbp_yr_g (nlat))
allocate(ctem_grd_yr%hetrores_yr_g (nlat))
allocate(ctem_grd_yr%autores_yr_g (nlat))
allocate(ctem_grd_yr%litres_yr_g (nlat))
allocate(ctem_grd_yr%soilcres_yr_g (nlat))
allocate(ctem_grd_yr%vgbiomas_yr_g (nlat))
allocate(ctem_grd_yr%totcmass_yr_g (nlat))
allocate(ctem_grd_yr%emit_co2_yr_g (nlat))
allocate(ctem_grd_yr%emit_co_yr_g (nlat))
allocate(ctem_grd_yr%emit_ch4_yr_g (nlat))
allocate(ctem_grd_yr%emit_nmhc_yr_g (nlat))
allocate(ctem_grd_yr%emit_h2_yr_g (nlat))
allocate(ctem_grd_yr%emit_nox_yr_g (nlat))
allocate(ctem_grd_yr%emit_n2o_yr_g (nlat))
allocate(ctem_grd_yr%emit_pm25_yr_g (nlat))
allocate(ctem_grd_yr%emit_tpm_yr_g (nlat))
allocate(ctem_grd_yr%emit_tc_yr_g (nlat))
allocate(ctem_grd_yr%emit_oc_yr_g (nlat))
allocate(ctem_grd_yr%emit_bc_yr_g (nlat))
allocate(ctem_grd_yr%smfuncveg_yr_g (nlat))
allocate(ctem_grd_yr%luc_emc_yr_g (nlat))
allocate(ctem_grd_yr%lucltrin_yr_g (nlat))
allocate(ctem_grd_yr%lucsocin_yr_g (nlat))
allocate(ctem_grd_yr%burnfrac_yr_g (nlat))
allocate(ctem_grd_yr%bterm_yr_g (nlat))
allocate(ctem_grd_yr%lterm_yr_g (nlat))
allocate(ctem_grd_yr%mterm_yr_g (nlat))
allocate(ctem_grd_yr%ch4wet1_yr_g (nlat))
allocate(ctem_grd_yr%ch4wet2_yr_g (nlat))
allocate(ctem_grd_yr%wetfdyn_yr_g (nlat))
allocate(ctem_grd_yr%ch4dyn1_yr_g (nlat))
allocate(ctem_grd_yr%ch4dyn2_yr_g (nlat))
allocate(ctem_grd_yr%ch4soills_yr_g (nlat))
allocate(ctem_grd_yr%veghght_yr_g (nlat))

allocate(ctem_tile_yr%laimaxg_yr_t (nlat,nmos))
allocate(ctem_tile_yr%stemmass_yr_t (nlat,nmos))
allocate(ctem_tile_yr%rootmass_yr_t (nlat,nmos))
allocate(ctem_tile_yr%npp_yr_t (nlat,nmos))
allocate(ctem_tile_yr%gpp_yr_t (nlat,nmos))
allocate(ctem_tile_yr%vgbiomas_yr_t (nlat,nmos))
allocate(ctem_tile_yr%autores_yr_t (nlat,nmos))
allocate(ctem_tile_yr%totcmass_yr_t (nlat,nmos))
allocate(ctem_tile_yr%litrmass_yr_t (nlat,nmos))
allocate(ctem_tile_yr%soilcmas_yr_t (nlat,nmos))
allocate(ctem_tile_yr%nep_yr_t (nlat,nmos))
allocate(ctem_tile_yr%litres_yr_t (nlat,nmos))
allocate(ctem_tile_yr%soilcres_yr_t (nlat,nmos))
allocate(ctem_tile_yr%hetrores_yr_t (nlat,nmos))
allocate(ctem_tile_yr%nbp_yr_t (nlat,nmos))
allocate(ctem_tile_yr%emit_co2_yr_t (nlat,nmos))
allocate(ctem_tile_yr%emit_co_yr_t (nlat,nmos))
allocate(ctem_tile_yr%emit_ch4_yr_t (nlat,nmos))
allocate(ctem_tile_yr%emit_nmhc_yr_t (nlat,nmos))
allocate(ctem_tile_yr%emit_h2_yr_t (nlat,nmos))
allocate(ctem_tile_yr%emit_nox_yr_t (nlat,nmos))
allocate(ctem_tile_yr%emit_n2o_yr_t (nlat,nmos))
allocate(ctem_tile_yr%emit_pm25_yr_t (nlat,nmos))
allocate(ctem_tile_yr%emit_tpm_yr_t (nlat,nmos))
allocate(ctem_tile_yr%emit_tc_yr_t (nlat,nmos))
allocate(ctem_tile_yr%emit_oc_yr_t (nlat,nmos))
allocate(ctem_tile_yr%emit_bc_yr_t (nlat,nmos))
allocate(ctem_tile_yr%burnfrac_yr_t (nlat,nmos))
allocate(ctem_tile_yr%smfuncveg_yr_t (nlat,nmos))
allocate(ctem_tile_yr%bterm_yr_t (nlat,nmos))
allocate(ctem_tile_yr%luc_emc_yr_t (nlat,nmos))
allocate(ctem_tile_yr%lterm_yr_t (nlat,nmos))
allocate(ctem_tile_yr%lucsocin_yr_t (nlat,nmos))
allocate(ctem_tile_yr%mterm_yr_t (nlat,nmos))
allocate(ctem_tile_yr%lucltrin_yr_t (nlat,nmos))
allocate(ctem_tile_yr%ch4wet1_yr_t (nlat,nmos))
allocate(ctem_tile_yr%ch4wet2_yr_t (nlat,nmos))
allocate(ctem_tile_yr%wetfdyn_yr_t (nlat,nmos))
allocate(ctem_tile_yr%ch4dyn1_yr_t (nlat,nmos))
allocate(ctem_tile_yr%ch4dyn2_yr_t (nlat,nmos))
allocate(ctem_tile_yr%ch4soills_yr_t (nlat,nmos))
allocate(ctem_tile_yr%veghght_yr_t (nlat,nmos))

end subroutine alloc_ctem_vars

! -----------------------------------------------------


subroutine initrowvars()

use ctem_params, only : nlat, nmos, ican, ignd ,icc, iccp1

implicit none

integer :: j,k,l,m

 do j = 1,nlat

   do k = 1,nmos

!         vrot%PREACC_M(j,k) = 0.
!         vrot%GTACC_M(j,k) = 0.
!         vrot%QEVPACC_M(j,k) = 0.
!         vrot%HFSACC_M(j,k) = 0.
!         vrot%HMFNACC_M(j,k) = 0.
!         vrot%ROFACC_M(j,k) = 0.
!         vrot%SNOACC_M(j,k) = 0.
!         vrot%OVRACC_M(j,k) = 0.
!         vrot%WTBLACC_M(j,k) = 0.
!
!         vrot%ALVSACC_M(j,k) = 0.
!         vrot%ALIRACC_M(j,k) = 0.
!         vrot%RHOSACC_M(j,k) = 0.
!         vrot%TSNOACC_M(j,k) = 0.
!         vrot%WSNOACC_M(j,k) = 0.
!         vrot%TCANACC_M(j,k) = 0.
!         vrot%RCANACC_M(j,k) = 0.
!         vrot%SCANACC_M(j,k) = 0.
!         vrot%GROACC_M(j,k) = 0.
!         vrot%FSINACC_M(j,k) = 0.
!         vrot%FLINACC_M(j,k) = 0.
!         vrot%TAACC_M(j,k) = 0.
!         vrot%UVACC_M(j,k) = 0.
!         vrot%PRESACC_M(j,k) = 0.
!         vrot%QAACC_M(j,k) = 0.
!         vrot%EVAPACC_M(j,k) = 0.
!         vrot%FLUTACC_M(j,k) = 0.
        vrot%icount(j,k)           = 0
        vrot%co2conc(j,k)          = 0.0
        vrot%npp(j,k)              = 0.0
        vrot%nep(j,k)              = 0.0
        vrot%hetrores(j,k)         = 0.0
        vrot%autores(j,k)          = 0.0
        vrot%soilcresp(j,k)        = 0.0
        vrot%rm(j,k)               = 0.0
        vrot%rg(j,k)               = 0.0
        vrot%nbp(j,k)              = 0.0
        vrot%litres(j,k)           = 0.0
        vrot%socres(j,k)           = 0.0
        vrot%gpp(j,k)              = 0.0
        vrot%dstcemls(j,k)         = 0.0
        vrot%dstcemls3(j,k)        = 0.0
        vrot%litrfall(j,k)         = 0.0
        vrot%humiftrs(j,k)         = 0.0
        vrot%canres(j,k)           = 0.0
        vrot%rml(j,k)              = 0.0
        vrot%rms(j,k)              = 0.0
        vrot%rmr(j,k)              = 0.0
        vrot%lucemcom(j,k)         = 0.0
        vrot%lucltrin(j,k)         = 0.0
        vrot%lucsocin(j,k)         = 0.0
        vrot%burnfrac(j,k)         = 0.0
        vrot%lterm(j,k)            = 0.0
        vrot%cfluxcg(j,k)          = 0.0
        vrot%cfluxcs(j,k)          = 0.0
        !vrot%TCANOACC_M(j,k)       = 0.0
        !vrot%UVACC_M(j,k)          = 0.0
        !vrot%VVACC_M(j,k)          = 0.0
        !vrot%TCANOACC_OUT(j,k)     = 0.0
        vrot%ch4wet1(j,k)          = 0.0
        vrot%ch4wet2(j,k)          = 0.0
        vrot%wetfdyn(j,k)          = 0.0
        vrot%ch4dyn1(j,k)          = 0.0
        vrot%ch4dyn2(j,k)          = 0.0
        vrot%ch4_soills(j,k)       = 0.0


        do l=1,ignd
            vrot%tbaraccrow_m(j,k,l)  = 0.0
!             vrot%TBARACC_M(j,k,l) = 0.
!             vrot%THLQACC_M(j,k,l) = 0.
!             vrot%THICACC_M(j,k,l) = 0.
!             vrot%THALACC_M(j,k,l) = 0.
        end do

        do l=1,ican
            vrot%ZOLNC(j,k,l)        = 0.0
            vrot%AILC(j,k,l)         = 0.0
            vrot%CMASVEGC(j,k,l)     = 0.0
            vrot%ALVSCTM(j,k,l)      = 0.0
            vrot%ALIRCTM(j,k,l)      = 0.0
            vrot%CSUM(j,k,l)            = 0.0
            vrot%PAIC(j,k,l)         = 0.0
            vrot%SLAIC(j,k,l)        = 0.0

            do m = 1, 3
                vrot%RMATC(j,k,l,m)    = 0.0
            end do

        end do

        do l = 1,icc

            vrot%smfuncveg(j,k,l)         = 0.0
            vrot%ailcmin(j,k,l) = 0.
            vrot%ailcmax(j,k,l) = 0.
            vrot%dvdfcan(j,k,l) = 0.
            vrot%gleafmas(j,k,l) = 0.
            vrot%bleafmas(j,k,l) = 0.
            vrot%stemmass(j,k,l) = 0.
            vrot%rootmass(j,k,l) = 0.
            vrot%pstemmass(j,k,l) = 0.
            vrot%pgleafmass(j,k,l) = 0.
            vrot%litrfallveg(j,k,l)=0.
            vrot%bterm(j,k,l)        = 0.0
            vrot%mterm(j,k,l)        = 0.0
            vrot%ailcg(j,k,l)        = 0.0
            vrot%ailcgs(j,k,l)       = 0.0
            vrot%fcancs(j,k,l)       = 0.0
            vrot%fcanc(j,k,l)        = 0.0
            vrot%fcancmx(j,k,l)      = 0.0
            vrot%co2i1cg(j,k,l)      = 0.0
            vrot%co2i1cs(j,k,l)      = 0.0
            vrot%co2i2cg(j,k,l)      = 0.0
            vrot%co2i2cs(j,k,l)      = 0.0
            vrot%ancsveg(j,k,l)      = 0.0
            vrot%ancgveg(j,k,l)      = 0.0
            vrot%rmlcsveg(j,k,l)     = 0.0
            vrot%rmlcgveg(j,k,l)     = 0.0
            vrot%stemmass(j,k,l)     = 0.0
            vrot%rootmass(j,k,l)     = 0.0
            vrot%ailcb(j,k,l)        = 0.0
            vrot%grwtheff(j,k,l)     = 0.0
            vrot%dvdfcan(j,k,l)      = 0.0
            vrot%bmasveg(j,k,l)      = 0.0
            vrot%tltrleaf(j,k,l)     = 0.0
            vrot%tltrstem(j,k,l)     = 0.0
            vrot%tltrroot(j,k,l)     = 0.0
            vrot%leaflitr(j,k,l)     = 0.0
            vrot%roottemp(j,k,l)     = 0.0
            vrot%afrleaf(j,k,l)      = 0.0
            vrot%afrstem(j,k,l)      = 0.0
            vrot%afrroot(j,k,l)      = 0.0
            vrot%wtstatus(j,k,l)     = 0.0
            vrot%ltstatus(j,k,l)     = 0.0
            vrot%ailcmin(j,k,l)      = 0.0
            vrot%ailcmax(j,k,l)      = 0.0
            vrot%pfcancmx(j,k,l)     = 0.0
            vrot%nfcancmx(j,k,l)     = 0.0
            vrot%nppveg(j,k,l)       = 0.0
            vrot%veghght(j,k,l)      = 0.0
            vrot%rootdpth(j,k,l)     = 0.0
            vrot%gleafmas(j,k,l)     = 0.0
            vrot%bleafmas(j,k,l)     = 0.0
            vrot%anveg(j,k,l)        = 0.0
            vrot%rmlveg(j,k,l)       = 0.0
            vrot%rmlvegacc(j,k,l)    = 0.0
            vrot%rmsveg(j,k,l)       = 0.0
            vrot%rmrveg(j,k,l)       = 0.0
            vrot%rgveg(j,k,l)        = 0.0
            vrot%vgbiomas_veg(j,k,l) = 0.0
            vrot%gppveg(j,k,l) = 0.0
            vrot%autoresveg(j,k,l) = 0.0
            vrot%emit_co2(j,k,l)         =0.0
            vrot%emit_co(j,k,l)          =0.0
            vrot%emit_ch4(j,k,l)         =0.0
            vrot%emit_nmhc(j,k,l)        =0.0
            vrot%emit_h2(j,k,l)          =0.0
            vrot%emit_nox(j,k,l)         =0.0
            vrot%emit_n2o(j,k,l)         =0.0
            vrot%emit_pm25(j,k,l)        =0.0
            vrot%emit_tpm(j,k,l)         =0.0
            vrot%emit_tc(j,k,l)          =0.0
            vrot%emit_oc(j,k,l)          =0.0
            vrot%emit_bc(j,k,l)          =0.0
            vrot%burnvegf(j,k,l)         =0.0

                do m = 1, ignd
                    vrot%rmatctem(j,k,l,m) = 0.0
                end do

            end do !icc
            !
            do l = 1, iccp1
                vrot%litrmass(j,k,l)    = 0.0
                vrot%soilcmas(j,k,l)    = 0.0
                vrot%hetroresveg(j,k,l) = 0.0
                vrot%litresveg(j,k,l) = 0.0
                vrot%soilcresveg(j,k,l) = 0.0
                vrot%nepveg(j,k,l) = 0.0
                vrot%nbpveg(j,k,l) = 0.0
                vrot%humiftrsveg(j,k,l)=0.
            end do !iccp1

   end do !nmos
 end do !nlat

end subroutine initrowvars

!==================================================

subroutine resetclassmon(nltest)

use ctem_params, only : ignd

implicit none

integer, intent(in) :: nltest

integer :: i,j

do i=1,nltest
    class_out%ALVSACC_MO(I)=0.
    class_out%ALIRACC_MO(I)=0.
    class_out%FLUTACC_MO(I)=0.
    class_out%FSINACC_MO(I)=0.
    class_out%FLINACC_MO(I)=0.
    class_out%HFSACC_MO(I) =0.
    class_out%QEVPACC_MO(I)=0.
    class_out%TRANSPACC_MO(I)=0.
    class_out%SNOACC_MO(I) =0.
    class_out%WSNOACC_MO(I)=0.
    class_out%ROFACC_MO(I) =0.
    class_out%PREACC_MO(I) =0.
    class_out%EVAPACC_MO(I)=0.
    class_out%TAACC_MO(I)=0.
    class_out%ACTLYR_MO(I)=0.
    class_out%FTABLE_MO(I)=0.
    class_out%ACTLYR_MIN_MO(I)=100000.
    class_out%FTABLE_MIN_MO(I)=100000.
    class_out%ACTLYR_MAX_MO(I)=0.
    class_out%FTABLE_MAX_MO(I)=0.
    class_out%CANOPYEVAP(I)=0.
    class_out%GROUNDEVAP(I)=0.
    class_out%ALTOTACC_MO(I)=0.
    class_out%altotcntr_m(i)=0


    DO J=1,IGND
        class_out%TBARACC_MO(I,J)=0.
        class_out%THLQACC_MO(I,J)=0.
        class_out%THICACC_MO(I,J)=0.
    end do
end do

end subroutine resetclassmon

!==================================================

subroutine resetclassyr(nltest)

implicit none

integer, intent(in) :: nltest

integer :: i

do i=1,nltest
          class_out%ALVSACC_YR(I)=0.
          class_out%ALIRACC_YR(I)=0.
          class_out%FLUTACC_YR(I)=0.
          class_out%FSINACC_YR(I)=0.
          class_out%FLINACC_YR(I)=0.
          class_out%HFSACC_YR(I) =0.
          class_out%QEVPACC_YR(I)=0.
          class_out%ROFACC_YR(I) =0.
          class_out%PREACC_YR(I) =0.
          class_out%EVAPACC_YR(I)=0.
          class_out%TRANSPACC_YR(I)=0.
          class_out%TAACC_YR(I)=0.
          class_out%ALTOTACC_YR(I)=0.
          class_out%altotcntr_yr(i)=0
end do

end subroutine resetclassyr

!==================================================

subroutine resetdaily(nltest,nmtest)

use ctem_params, only : ignd,icc

implicit none

integer, intent(in) :: nltest
integer, intent(in) :: nmtest

integer :: i,j,k,m

! First reset the grid average
do i=1,nltest
    ctem_grd%gpp_g(i) =0.0
    ctem_grd%npp_g(i) =0.0
    ctem_grd%nep_g(i) =0.0
    ctem_grd%nbp_g(i) =0.0
    ctem_grd%autores_g(i) =0.0
    ctem_grd%hetrores_g(i)=0.0
    ctem_grd%litres_g(i) =0.0
    ctem_grd%socres_g(i) =0.0
    ctem_grd%dstcemls_g(i)=0.0
    ctem_grd%dstcemls3_g(i)=0.0
    ctem_grd%litrfall_g(i)=0.0
    ctem_grd%humiftrs_g(i)=0.0
    ctem_grd%rml_g(i) =0.0
    ctem_grd%rms_g(i) =0.0
    ctem_grd%rmr_g(i) =0.0
    ctem_grd%rg_g(i) =0.0
    ctem_grd%vgbiomas_g(i) =0.0
    ctem_grd%totcmass_g(i) =0.0
    ctem_grd%gavglai_g(i) =0.0
    ctem_grd%gavgltms_g(i) =0.0
    ctem_grd%gavgscms_g(i) =0.0
    ctem_grd%ailcg_g(i)=0.0
    ctem_grd%ailcb_g(i)=0.0
    ctem_grd%tcanoacc_out_g(i) =0.0
    ctem_grd%burnfrac_g(i) =0.0
    ctem_grd%smfuncveg_g(i) =0.0
    ctem_grd%lucemcom_g(i) =0.0
    ctem_grd%lucltrin_g(i) =0.0
    ctem_grd%lucsocin_g(i) =0.0
    ctem_grd%emit_co2_g(i) =0.0
    ctem_grd%emit_co_g(i)  =0.0
    ctem_grd%emit_ch4_g(i) =0.0
    ctem_grd%emit_nmhc_g(i) =0.0
    ctem_grd%emit_h2_g(i) =0.0
    ctem_grd%emit_nox_g(i) =0.0
    ctem_grd%emit_n2o_g(i) =0.0
    ctem_grd%emit_pm25_g(i) =0.0
    ctem_grd%emit_tpm_g(i) =0.0
    ctem_grd%emit_tc_g(i) =0.0
    ctem_grd%emit_oc_g(i) =0.0
    ctem_grd%emit_bc_g(i) =0.0
    ctem_grd%bterm_g(i)   =0.0
    ctem_grd%lterm_g(i)   =0.0
    ctem_grd%mterm_g(i)   =0.0
    ctem_grd%leaflitr_g(i)=0.0
    ctem_grd%tltrleaf_g(i)=0.0
    ctem_grd%tltrstem_g(i)=0.0
    ctem_grd%tltrroot_g(i)=0.0
    ctem_grd%gleafmas_g(i)=0.0
    ctem_grd%bleafmas_g(i)=0.0
    ctem_grd%stemmass_g(i)=0.0
    ctem_grd%rootmass_g(i)=0.0
    ctem_grd%litrmass_g(i)=0.0
    ctem_grd%soilcmas_g(i)=0.0
    ctem_grd%veghght_g(i)=0.0
    ctem_grd%rootdpth_g(i)=0.0
    ctem_grd%roottemp_g(i)=0.0
    ctem_grd%slai_g(i)=0.0
    ctem_grd%CH4WET1_G(i) = 0.0
    ctem_grd%CH4WET2_G(i) = 0.0
    ctem_grd%WETFDYN_G(i) = 0.0
    ctem_grd%CH4DYN1_G(i) = 0.0
    ctem_grd%CH4DYN2_G(i) = 0.0
    ctem_grd%ch4_soills_g(i) = 0.0

    do k=1,ignd
      ctem_grd%rmatctem_g(i,k)=0.0
    enddo

    do j=1,icc
      ctem_grd%afrleaf_g(i,j)=0.0
      ctem_grd%afrstem_g(i,j)=0.0
      ctem_grd%afrroot_g(i,j)=0.0
    enddo

    do m = 1, nmtest

        ctem_tile%leaflitr_t(i,m)=0.0
        ctem_tile%tltrleaf_t(i,m)=0.0
        ctem_tile%tltrstem_t(i,m)=0.0
        ctem_tile%tltrroot_t(i,m)=0.0
        ctem_tile%ailcg_t(i,m)=0.0
        ctem_tile%ailcb_t(i,m)=0.0
        ctem_tile%afrleaf_t(i,m)=0.0
        ctem_tile%afrstem_t(i,m)=0.0
        ctem_tile%afrroot_t(i,m)=0.0
        ctem_tile%veghght_t(i,m)=0.0
        ctem_tile%rootdpth_t(i,m)=0.0
        ctem_tile%roottemp_t(i,m)=0.0
        ctem_tile%slai_t(i,m)=0.0
        ctem_tile%gleafmas_t(i,m) = 0.0
        ctem_tile%bleafmas_t(i,m) = 0.0
        ctem_tile%stemmass_t(i,m) = 0.0
        ctem_tile%rootmass_t(i,m) = 0.0
        ctem_tile%litrmass_t(i,m) = 0.0
        ctem_tile%soilcmas_t(i,m) = 0.0
        ctem_tile%emit_co2_t(i,m) = 0.0
        ctem_tile%emit_co_t(i,m) = 0.0
        ctem_tile%emit_ch4_t(i,m) = 0.0
        ctem_tile%emit_nmhc_t(i,m) = 0.0
        ctem_tile%emit_h2_t(i,m) = 0.0
        ctem_tile%emit_nox_t(i,m) = 0.0
        ctem_tile%emit_n2o_t(i,m) = 0.0
        ctem_tile%emit_pm25_t(i,m) = 0.0
        ctem_tile%emit_tpm_t(i,m) = 0.0
        ctem_tile%emit_tc_t(i,m) = 0.0
        ctem_tile%emit_oc_t(i,m) = 0.0
        ctem_tile%emit_bc_t(i,m) = 0.0

        do k=1,ignd
            ctem_tile%rmatctem_t(i,m,k)=0.0
        enddo

    end do !nmtest
end do !nltest

end subroutine resetdaily

!==================================================
subroutine resetmonthend(nltest,nmtest)

use ctem_params, only : iccp1,icc

implicit none

integer, intent(in) :: nltest
integer, intent(in) :: nmtest

integer :: i,m,j

! These are assigned to mid-month, but are not accumulated so can be
! zeroed out at the same time as the other month-end vars.
do i=1,nltest
    ctem_grd_mo%stemmass_mo_g(i)=0.0
    ctem_grd_mo%rootmass_mo_g(i)=0.0
    ctem_grd_mo%litrmass_mo_g(i)=0.0
    ctem_grd_mo%soilcmas_mo_g(i)=0.0
    ctem_grd_mo%vgbiomas_mo_g(i)=0.0
    ctem_grd_mo%totcmass_mo_g(i)=0.0
  do m = 1, nmtest
        ctem_tile_mo%stemmass_mo_t(i,m)=0.0
        ctem_tile_mo%rootmass_mo_t(i,m)=0.0
        ctem_tile_mo%litrmass_mo_t(i,m)=0.0
        ctem_tile_mo%soilcmas_mo_t(i,m)=0.0
        ctem_tile_mo%vgbiomas_mo_t(i,m)=0.0
        ctem_tile_mo%totcmass_mo_t(i,m)=0.0
        do j = 1,icc
            ctem_mo%stemmass_mo(i,m,j)=0.0
            ctem_mo%rootmass_mo(i,m,j)=0.0
            ctem_mo%litrmass_mo(i,m,j)=0.0
            ctem_mo%soilcmas_mo(i,m,j)=0.0
            ctem_mo%vgbiomas_mo(i,m,j)=0.0
            ctem_mo%totcmass_mo(i,m,j)=0.0
        end do
            ctem_mo%litrmass_mo(i,m,iccp1)=0.0
            ctem_mo%soilcmas_mo(i,m,iccp1)=0.0
            ctem_mo%totcmass_mo(i,m,iccp1)=0.0
  end do
end do

! Now zero out the month end vars.
do i=1,nltest
    ! Grid avg
    ctem_grd_mo%laimaxg_mo_g(i)=0.0
    ctem_grd_mo%npp_mo_g(i)=0.0
    ctem_grd_mo%gpp_mo_g(i)=0.0
    ctem_grd_mo%nep_mo_g(i)=0.0
    ctem_grd_mo%nbp_mo_g(i)=0.0
    ctem_grd_mo%hetrores_mo_g(i)=0.0
    ctem_grd_mo%autores_mo_g(i)=0.0
    ctem_grd_mo%litres_mo_g(i)=0.0
    ctem_grd_mo%soilcres_mo_g(i)=0.0
    ctem_grd_mo%litrfall_mo_g(i)=0.0
    ctem_grd_mo%humiftrs_mo_g(i)=0.0
    ctem_grd_mo%emit_co2_mo_g(i)=0.0
    ctem_grd_mo%emit_co_mo_g(i) =0.0
    ctem_grd_mo%emit_ch4_mo_g(i) =0.0
    ctem_grd_mo%emit_nmhc_mo_g(i) =0.0
    ctem_grd_mo%emit_h2_mo_g(i) =0.0
    ctem_grd_mo%emit_nox_mo_g(i) =0.0
    ctem_grd_mo%emit_n2o_mo_g(i) =0.0
    ctem_grd_mo%emit_pm25_mo_g(i) =0.0
    ctem_grd_mo%emit_tpm_mo_g(i) =0.0
    ctem_grd_mo%emit_tc_mo_g(i) =0.0
    ctem_grd_mo%emit_oc_mo_g(i) =0.0
    ctem_grd_mo%emit_bc_mo_g(i) =0.0
    ctem_grd_mo%smfuncveg_mo_g(i) =0.0
    ctem_grd_mo%luc_emc_mo_g(i) =0.0
    ctem_grd_mo%lucsocin_mo_g(i) =0.0
    ctem_grd_mo%lucltrin_mo_g(i) =0.0
    ctem_grd_mo%burnfrac_mo_g(i) =0.0
    ctem_grd_mo%bterm_mo_g(i)    =0.0
    ctem_grd_mo%lterm_mo_g(i)    =0.0
    ctem_grd_mo%mterm_mo_g(i)    =0.0
    ctem_grd_mo%ch4wet1_mo_g(i)  =0.0
    ctem_grd_mo%ch4wet2_mo_g(i)  =0.0
    ctem_grd_mo%wetfdyn_mo_g(i)  =0.0
    ctem_grd_mo%ch4dyn1_mo_g(i)  =0.0
    ctem_grd_mo%ch4dyn2_mo_g(i)  =0.0
    ctem_grd_mo%ch4soills_mo_g(i)  =0.0

    do m = 1,nmtest
        ! Tile avg
        ctem_tile_mo%laimaxg_mo_t(i,m)=0.0
        ctem_tile_mo%npp_mo_t(i,m)=0.0
        ctem_tile_mo%gpp_mo_t(i,m)=0.0
        ctem_tile_mo%nep_mo_t(i,m)=0.0
        ctem_tile_mo%nbp_mo_t(i,m)=0.0
        ctem_tile_mo%hetrores_mo_t(i,m)=0.0
        ctem_tile_mo%autores_mo_t(i,m)=0.0
        ctem_tile_mo%litres_mo_t(i,m)=0.0
        ctem_tile_mo%soilcres_mo_t(i,m)=0.0
        ctem_tile_mo%litrfall_mo_t(i,m)=0.0
        ctem_tile_mo%humiftrs_mo_t(i,m)=0.0
        ctem_tile_mo%emit_co2_mo_t(i,m)=0.0
        ctem_tile_mo%emit_co_mo_t(i,m) =0.0
        ctem_tile_mo%emit_ch4_mo_t(i,m) =0.0
        ctem_tile_mo%emit_nmhc_mo_t(i,m) =0.0
        ctem_tile_mo%emit_h2_mo_t(i,m) =0.0
        ctem_tile_mo%emit_nox_mo_t(i,m) =0.0
        ctem_tile_mo%emit_n2o_mo_t(i,m) =0.0
        ctem_tile_mo%emit_pm25_mo_t(i,m) =0.0
        ctem_tile_mo%emit_tpm_mo_t(i,m) =0.0
        ctem_tile_mo%emit_tc_mo_t(i,m) =0.0
        ctem_tile_mo%emit_oc_mo_t(i,m) =0.0
        ctem_tile_mo%emit_bc_mo_t(i,m) =0.0
        ctem_tile_mo%smfuncveg_mo_t(i,m) =0.0
        ctem_tile_mo%luc_emc_mo_t(i,m) =0.0
        ctem_tile_mo%lucsocin_mo_t(i,m) =0.0
        ctem_tile_mo%lucltrin_mo_t(i,m) =0.0
        ctem_tile_mo%burnfrac_mo_t(i,m) =0.0
        ctem_tile_mo%bterm_mo_t(i,m)    =0.0
        ctem_tile_mo%lterm_mo_t(i,m)    =0.0
        ctem_tile_mo%mterm_mo_t(i,m)    =0.0
        ctem_tile_mo%ch4wet1_mo_t(i,m)  =0.0
        ctem_tile_mo%ch4wet2_mo_t(i,m)  =0.0
        ctem_tile_mo%wetfdyn_mo_t(i,m)  =0.0
        ctem_tile_mo%ch4dyn1_mo_t(i,m)  =0.0
        ctem_tile_mo%ch4dyn2_mo_t(i,m)  =0.0
        ctem_tile_mo%ch4soills_mo_t(i,m)  =0.0
        ctem_tile_mo%wind_mo_t(i,m) = 0.0

        do j=1,icc
            ! per pft
            ctem_mo%laimaxg_mo(i,m,j)=0.0
            ctem_mo%npp_mo(i,m,j)=0.0
            ctem_mo%gpp_mo(i,m,j)=0.0
            ctem_mo%nep_mo(i,m,j)=0.0
            ctem_mo%nbp_mo(i,m,j)=0.0
            ctem_mo%hetrores_mo(i,m,j)=0.0
            ctem_mo%autores_mo(i,m,j)=0.0
            ctem_mo%litres_mo(i,m,j)=0.0
            ctem_mo%soilcres_mo(i,m,j)=0.0
            ctem_mo%litrfallveg_mo(i,m,j)=0.0
            ctem_mo%humiftrsveg_mo(i,m,j)=0.0
            ctem_mo%emit_co2_mo(i,m,j)=0.0
            ctem_mo%emit_co_mo(i,m,j) =0.0
            ctem_mo%emit_ch4_mo(i,m,j) =0.0
            ctem_mo%emit_nmhc_mo(i,m,j) =0.0
            ctem_mo%emit_h2_mo(i,m,j) =0.0
            ctem_mo%emit_nox_mo(i,m,j) =0.0
            ctem_mo%emit_n2o_mo(i,m,j) =0.0
            ctem_mo%emit_pm25_mo(i,m,j) =0.0
            ctem_mo%emit_tpm_mo(i,m,j) =0.0
            ctem_mo%emit_tc_mo(i,m,j) =0.0
            ctem_mo%emit_oc_mo(i,m,j) =0.0
            ctem_mo%emit_bc_mo(i,m,j) =0.0
            ctem_mo%burnfrac_mo(i,m,j) =0.0
            ctem_mo%bterm_mo(i,m,j) =0.0
            ctem_mo%mterm_mo(i,m,j) =0.0
            ctem_mo%smfuncveg_mo(i,m,j) =0.0
        end do

        ctem_mo%nep_mo(i,m,iccp1)=0.0
        ctem_mo%nbp_mo(i,m,iccp1)=0.0
        ctem_mo%hetrores_mo(i,m,iccp1)=0.0
        ctem_mo%litres_mo(i,m,iccp1)=0.0
        ctem_mo%soilcres_mo(i,m,iccp1)=0.0
        ctem_mo%humiftrsveg_mo(i,m,iccp1)=0.0

    end do !nmtest
end do ! nltest

end subroutine resetmonthend

!==================================================

subroutine resetyearend(nltest,nmtest)

use ctem_params, only : iccp1,icc

implicit none

integer, intent(in) :: nltest
integer, intent(in) :: nmtest

integer :: i,m,j

do i=1,nltest
    ! Grid avg
    ctem_grd_yr%laimaxg_yr_g(i)=0.0
    ctem_grd_yr%stemmass_yr_g(i)=0.0
    ctem_grd_yr%rootmass_yr_g(i)=0.0
    ctem_grd_yr%litrmass_yr_g(i)=0.0
    ctem_grd_yr%soilcmas_yr_g(i)=0.0
    ctem_grd_yr%vgbiomas_yr_g(i)=0.0
    ctem_grd_yr%totcmass_yr_g(i)=0.0
    ctem_grd_yr%veghght_yr_g(i)=0.0
    ctem_grd_yr%npp_yr_g(i)=0.0
    ctem_grd_yr%gpp_yr_g(i)=0.0
    ctem_grd_yr%nep_yr_g(i)=0.0
    ctem_grd_yr%nbp_yr_g(i)=0.0
    ctem_grd_yr%hetrores_yr_g(i)=0.0
    ctem_grd_yr%autores_yr_g(i)=0.0
    ctem_grd_yr%litres_yr_g(i)=0.0
    ctem_grd_yr%soilcres_yr_g(i)=0.0
    ctem_grd_yr%emit_co2_yr_g(i)=0.0
    ctem_grd_yr%emit_co_yr_g(i)=0.0
    ctem_grd_yr%emit_ch4_yr_g(i)=0.0
    ctem_grd_yr%emit_nmhc_yr_g(i)=0.0
    ctem_grd_yr%emit_h2_yr_g(i)=0.0
    ctem_grd_yr%emit_nox_yr_g(i)=0.0
    ctem_grd_yr%emit_n2o_yr_g(i)=0.0
    ctem_grd_yr%emit_pm25_yr_g(i)=0.0
    ctem_grd_yr%emit_tpm_yr_g(i)=0.0
    ctem_grd_yr%emit_tc_yr_g(i)=0.0
    ctem_grd_yr%emit_oc_yr_g(i)=0.0
    ctem_grd_yr%emit_bc_yr_g(i)=0.0
    ctem_grd_yr%smfuncveg_yr_g(i)=0.0
    ctem_grd_yr%luc_emc_yr_g(i)=0.0
    ctem_grd_yr%lucsocin_yr_g(i)=0.0
    ctem_grd_yr%lucltrin_yr_g(i)=0.0
    ctem_grd_yr%burnfrac_yr_g(i)=0.0
    ctem_grd_yr%bterm_yr_g(i)=0.0
    ctem_grd_yr%lterm_yr_g(i)=0.0
    ctem_grd_yr%mterm_yr_g(i)=0.0
    ctem_grd_yr%ch4wet1_yr_g(i)  =0.0
    ctem_grd_yr%ch4wet2_yr_g(i)  =0.0
    ctem_grd_yr%wetfdyn_yr_g(i)  =0.0
    ctem_grd_yr%ch4dyn1_yr_g(i)  =0.0
    ctem_grd_yr%ch4dyn2_yr_g(i)  =0.0
    ctem_grd_yr%ch4soills_yr_g(i)  =0.0

    do m = 1,nmtest
        ! Tile avg
        ctem_tile_yr%laimaxg_yr_t(i,m)=0.0
        ctem_tile_yr%stemmass_yr_t(i,m)=0.0
        ctem_tile_yr%rootmass_yr_t(i,m)=0.0
        ctem_tile_yr%litrmass_yr_t(i,m)=0.0
        ctem_tile_yr%soilcmas_yr_t(i,m)=0.0
        ctem_tile_yr%vgbiomas_yr_t(i,m)=0.0
        ctem_tile_yr%totcmass_yr_t(i,m)=0.0
        ctem_tile_yr%veghght_yr_t(i,m)=0.0
        ctem_tile_yr%npp_yr_t(i,m)=0.0
        ctem_tile_yr%gpp_yr_t(i,m)=0.0
        ctem_tile_yr%nep_yr_t(i,m)=0.0
        ctem_tile_yr%nbp_yr_t(i,m)=0.0
        ctem_tile_yr%hetrores_yr_t(i,m)=0.0
        ctem_tile_yr%autores_yr_t(i,m)=0.0
        ctem_tile_yr%litres_yr_t(i,m)=0.0
        ctem_tile_yr%soilcres_yr_t(i,m)=0.0
        ctem_tile_yr%emit_co2_yr_t(i,m)=0.0
        ctem_tile_yr%emit_co_yr_t(i,m)=0.0
        ctem_tile_yr%emit_ch4_yr_t(i,m)=0.0
        ctem_tile_yr%emit_nmhc_yr_t(i,m)=0.0
        ctem_tile_yr%emit_h2_yr_t(i,m)=0.0
        ctem_tile_yr%emit_nox_yr_t(i,m)=0.0
        ctem_tile_yr%emit_n2o_yr_t(i,m)=0.0
        ctem_tile_yr%emit_pm25_yr_t(i,m)=0.0
        ctem_tile_yr%emit_tpm_yr_t(i,m)=0.0
        ctem_tile_yr%emit_tc_yr_t(i,m)=0.0
        ctem_tile_yr%emit_oc_yr_t(i,m)=0.0
        ctem_tile_yr%emit_bc_yr_t(i,m)=0.0
        ctem_tile_yr%smfuncveg_yr_t(i,m)=0.0
        ctem_tile_yr%luc_emc_yr_t(i,m)=0.0
        ctem_tile_yr%lucsocin_yr_t(i,m)=0.0
        ctem_tile_yr%lucltrin_yr_t(i,m)=0.0
        ctem_tile_yr%burnfrac_yr_t(i,m)=0.0
        ctem_tile_yr%bterm_yr_t(i,m)=0.0
        ctem_tile_yr%lterm_yr_t(i,m)=0.0
        ctem_tile_yr%mterm_yr_t(i,m)=0.0
        ctem_tile_yr%ch4wet1_yr_t(i,m)  =0.0
        ctem_tile_yr%ch4wet2_yr_t(i,m)  =0.0
        ctem_tile_yr%wetfdyn_yr_t(i,m)  =0.0
        ctem_tile_yr%ch4dyn1_yr_t(i,m)  =0.0
        ctem_tile_yr%ch4dyn2_yr_t(i,m)  =0.0
        ctem_tile_yr%ch4soills_yr_t(i,m)  =0.0

        do j=1,icc
            ! per pft
            ctem_yr%laimaxg_yr(i,m,j)=0.0
            ctem_yr%stemmass_yr(i,m,j)=0.0
            ctem_yr%rootmass_yr(i,m,j)=0.0
            ctem_yr%litrmass_yr(i,m,j)=0.0
            ctem_yr%soilcmas_yr(i,m,j)=0.0
            ctem_yr%vgbiomas_yr(i,m,j)=0.0
            ctem_yr%totcmass_yr(i,m,j)=0.0
            ctem_yr%veghght_yr(i,m,j)=0.0
            ctem_yr%npp_yr(i,m,j)=0.0
            ctem_yr%gpp_yr(i,m,j)=0.0
            ctem_yr%nep_yr(i,m,j)=0.0
            ctem_yr%nbp_yr(i,m,j)=0.0
            ctem_yr%hetrores_yr(i,m,j)=0.0
            ctem_yr%autores_yr(i,m,j)=0.0
            ctem_yr%litres_yr(i,m,j)=0.0
            ctem_yr%soilcres_yr(i,m,j)=0.0
            ctem_yr%emit_co2_yr(i,m,j)=0.0
            ctem_yr%emit_co_yr(i,m,j)=0.0
            ctem_yr%emit_ch4_yr(i,m,j)=0.0
            ctem_yr%emit_nmhc_yr(i,m,j)=0.0
            ctem_yr%emit_h2_yr(i,m,j)=0.0
            ctem_yr%emit_nox_yr(i,m,j)=0.0
            ctem_yr%emit_n2o_yr(i,m,j)=0.0
            ctem_yr%emit_pm25_yr(i,m,j)=0.0
            ctem_yr%emit_tpm_yr(i,m,j)=0.0
            ctem_yr%emit_tc_yr(i,m,j)=0.0
            ctem_yr%emit_oc_yr(i,m,j)=0.0
            ctem_yr%emit_bc_yr(i,m,j)=0.0
            ctem_yr%bterm_yr(i,m,j)=0.0
            ctem_yr%mterm_yr(i,m,j)=0.0
            ctem_yr%burnfrac_yr(i,m,j)=0.0
            ctem_yr%smfuncveg_yr(i,m,j)=0.0
        end do

        ctem_yr%hetrores_yr(i,m,iccp1)=0.0
        ctem_yr%litres_yr(i,m,iccp1)=0.0
        ctem_yr%soilcres_yr(i,m,iccp1)=0.0
        ctem_yr%nep_yr(i,m,iccp1)=0.0
        ctem_yr%nbp_yr(i,m,iccp1)=0.0
        ctem_yr%litrmass_yr(i,m,iccp1)=0.0
        ctem_yr%soilcmas_yr(i,m,iccp1)=0.0
        ctem_yr%totcmass_yr(i,m,iccp1)=0.0

    end do !nmtest
end do ! nltest

end subroutine resetyearend

!==================================================
subroutine resetclassaccum(nltest,nmtest)

use ctem_params, only : ignd

implicit none

integer, intent(in) :: nltest
integer, intent(in) :: nmtest

integer :: i,m,j


DO I=1,NLTEST

    vrot%altotcntr_d(i) = 0

  DO M=1,NMTEST

        vrot%PREACC_M(i,m) = 0.
        vrot%GTACC_M(i,m) = 0.
        vrot%QEVPACC_M(i,m) = 0.
        vrot%HFSACC_M(i,m) = 0.
        vrot%HMFNACC_M(i,m) = 0.
        vrot%ROFACC_M(i,m) = 0.
        vrot%SNOACC_M(i,m) = 0.
        vrot%OVRACC_M(i,m) = 0.
        vrot%WTBLACC_M(i,m) = 0.
        vrot%ALVSACC_M(i,m) = 0.
        vrot%ALIRACC_M(i,m) = 0.
        vrot%RHOSACC_M(i,m) = 0.
        vrot%TSNOACC_M(i,m) = 0.
        vrot%WSNOACC_M(i,m) = 0.
        vrot%SNOARE_M(i,m) = 0.
        vrot%TCANACC_M(i,m) = 0.
        vrot%RCANACC_M(i,m) = 0.
        vrot%SCANACC_M(i,m) = 0.
        vrot%GROACC_M(i,m) = 0.
        vrot%FSINACC_M(i,m) = 0.
        vrot%FLINACC_M(i,m) = 0.
        vrot%TAACC_M(i,m) = 0.
        vrot%UVACC_M(i,m) = 0.
        vrot%PRESACC_M(i,m) = 0.
        vrot%QAACC_M(i,m) = 0.
        vrot%ALTOTACC_M(i,m) = 0.
        vrot%EVAPACC_M(i,m) = 0.
        vrot%FLUTACC_M(i,m) = 0.

    DO J=1,IGND
        vrot%TBARACC_M(I,M,J)=0.
        vrot%THLQACC_M(I,M,J)=0.
        vrot%THICACC_M(I,M,J)=0.
        vrot%THALACC_M(I,M,J)=0.
    end do
  end do
end do

end subroutine resetclassaccum


!==================================================

subroutine resetgridavg(nltest)

use ctem_params, only : ignd,icc

implicit none

integer, intent(in) :: nltest

integer :: i,j

        do i = 1, nltest

        ctem_grd%fsstar_g(i) =0.0
        ctem_grd%flstar_g(i) =0.0
        ctem_grd%qh_g(i)     =0.0
        ctem_grd%qe_g(i)     =0.0
        ctem_grd%snomlt_g(i) =0.0
        ctem_grd%beg_g(i)    =0.0
        ctem_grd%gtout_g(i)  =0.0
        ctem_grd%SNOROT_g(i) =0.0
        ctem_grd%RHOSROT_g(i)=0.0
        ctem_grd%WSNOROT_g(i)=0.0
        ctem_grd%altot_g(i)  =0.0
        ctem_grd%ROFROT_g(i) =0.0
        ctem_grd%tpn_g(i)    =0.0
        ctem_grd%ZPNDROT_g(i)=0.0
        ctem_grd%tcn_g(i)=0.0
        ctem_grd%tsn_g(i)=0.0
        ctem_grd%zsn_g(i)=0.0

        do j=1,ignd
         ctem_grd%TBARROT_g(i,j)=0.0
         ctem_grd%THLQROT_g(i,j)=0.0
         ctem_grd%THICROT_g(i,j)=0.0
         ctem_grd%GFLXROT_g(i,j)=0.0
         ctem_grd%HMFGROT_g(i,j)=0.0
         ctem_grd%HTCROT_g(i,j)=0.0
         ctem_grd%QFCROT_g(i,j)=0.0
        end do

        
        ctem_grd%RCANROT_g(i) =0.0
        ctem_grd%SCANROT_g(i) =0.0
        ctem_grd%TROFROT_g(i)=0.0
        ctem_grd%TROOROT_g(i)=0.0
        ctem_grd%TROBROT_g(i)=0.0
        ctem_grd%TROSROT_g(i)=0.0
        ctem_grd%ROFOROT_g(i)=0.0
        ctem_grd%ROFSROT_g(i)=0.0
        ctem_grd%ROFBROT_g(i)=0.0
        ctem_grd%FSGVROT_g(i)=0.0
        ctem_grd%FSGSROT_g(i)=0.0
        ctem_grd%FSGGROT_g(i)=0.0
        ctem_grd%FLGVROT_g(i)=0.0
        ctem_grd%FLGSROT_g(i)=0.0
        ctem_grd%FLGGROT_g(i)=0.0
        ctem_grd%HFSCROT_g(i)=0.0
        ctem_grd%HFSSROT_g(i)=0.0
        ctem_grd%HFSGROT_g(i)=0.0
        ctem_grd%HEVCROT_g(i)=0.0
        ctem_grd%HEVSROT_g(i)=0.0
        ctem_grd%HEVGROT_g(i)=0.0
        ctem_grd%HMFCROT_g(i)=0.0
        ctem_grd%HMFNROT_g(i)=0.0
        ctem_grd%HTCCROT_g(i)=0.0
        ctem_grd%HTCSROT_g(i)=0.0
        ctem_grd%PCFCROT_g(i)=0.0
        ctem_grd%PCLCROT_g(i)=0.0
        ctem_grd%PCPNROT_g(i)=0.0
        ctem_grd%PCPGROT_g(i)=0.0
        ctem_grd%QFCFROT_g(i)=0.0
        ctem_grd%QFCLROT_g(i)=0.0
        ctem_grd%QFNROT_g(i)=0.0
        ctem_grd%QFGROT_g(i)=0.0
        ctem_grd%ROFCROT_g(i)=0.0
        ctem_grd%ROFNROT_g(i)=0.0
        ctem_grd%WTRCROT_g(i)=0.0
        ctem_grd%WTRSROT_g(i)=0.0
        ctem_grd%WTRGROT_g(i)=0.0
        ctem_grd%CDHROT_g(i)=0.0
        ctem_grd%CDMROT_g(i)=0.0
        ctem_grd%SFCUROT_g(i)=0.0
        ctem_grd%SFCVROT_g(i)=0.0
        ctem_grd%ACTLYR_g(i)=0.0
        ctem_grd%FTABLE_g(i)=0.0

       if (c_switch%ctem_on) then
          do j=1,icc
            ctem_grd%anvegrow_g(i,j)=0.0
            ctem_grd%rmlvegrow_g(i,j)=0.0
          end do
       end if

       end do

end subroutine resetgridavg



!==================================================

subroutine finddaylength(solday,radl,daylength)

! Calculate the daylength based on the latitude and day of year

! Joe Melton Dec 18 2015 (taken from phenlogy.f)

use ctem_params, only : pi

implicit none

real, intent(in) :: solday  !day of year
real, intent(in) :: radl    ! latitude
real, intent(out) :: daylength  ! calculated daylength
real :: theta               ! temp var
real :: decli               ! temp var
real :: term                ! temp var

    theta=0.2163108 + 2.0*atan(0.9671396*tan(0.0086*(solday-186.0)))
    decli=asin(0.39795*cos(theta))    !declination !note I see that CLASS does this also but with different formula...
    term=(sin(radl)*sin(decli))/(cos(radl)*cos(decli))
    term=max(-1.0,min(term,1.0))
    daylength=24.0-(24.0/pi)*acos(term)

end subroutine finddaylength

!==================================================


! separate one:
!
! c     reset mosaic accumulator arrays.
! c
!       do 655 i=1,nml
!          uvaccgat_t(i)=0.0
! 655   continue
! c
!       if (ctem_on) then
!         do 705 i = 1, nml
! c
! c         competitition related variables added by y. peng \\
!           fsinacc_gat(i)=0.
!           flinacc_gat(i)=0.
!           flutacc_gat(i)=0.
!           alswacc_gat(i)=0.
!           allwacc_gat(i)=0.
!           pregacc_gat(i)=0.
! c         competitition related variables added by y. peng //
! c
!           fsnowacc_t(i)=0.0
!           tcanoaccgat_out(i)=tcanoaccgat_t(i)
!           tcanoaccgat_t(i)=0.0
! c
!           tcansacc_t(i)=0.0
!           taaccgat_t(i)=0.0
!           vvaccgat_t(i)=0.0
! c
!           do 715 j=1,ignd
!              tbaraccgat_t(i,j)=0.0
!              tbarcacc_t(i,j)=0.0
!              tbarcsacc_t(i,j)=0.0
!              tbargacc_t(i,j)=0.0
!              tbargsacc_t(i,j)=0.0
!              thliqcacc_t(i,j)=0.0
!              thliqgacc_t(i,j)=0.0
!              thicecacc_t(i,j)=0.0
! 715       continue
! c
!           do 716 j = 1, icc
!             ancsvgac_t(i,j)=0.0
!             ancgvgac_t(i,j)=0.0
!             rmlcsvga_t(i,j)=0.0
!             rmlcgvga_t(i,j)=0.0
! 716       continue
! c
! 705     continue
!       endif  ! if(ctem_on)
!       endif  ! if(ncount.eq.nday)
!=================================================================================
!>@}
end module ctem_statevars
