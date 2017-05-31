!>\file
!! Principle driver program to run CLASSIC in stand-alone mode using specified boundary
!! conditions and atmospheric forcing.
!!
!! # Overview
!!
!! This driver program initializes the run, reads in CLASS input files, manages the run
!! and the coupling between CLASS and CTEM, writes the CLASS sub-monthly outputs, and
!! closes the run.

module main_driver

implicit none

public :: CLASSIC_driver

contains

subroutine CLASSIC_driver()
!
!>
!!------------------------------------------------------------------
!! ## Dimension statements.

!!     ### first set of definitions:
!!     background variables, and prognostic and diagnostic
!!     variables normally provided by and/or used by the gcm.
!!      the suffix "rot" refers to variables existing on the
!!      mosaic grid on the current latitude circle.  the suffix
!!      "gat" refers to the same variables after they have undergone
!!      a "gather" operation in which the two mosaic dimensions
!!      are collapsed into one.  the suffix "row" refers both to
!!      grid-constant input variables. and to grid-averaged
!!      diagnostic variables.
!!
!!      the first dimension element of the "rot" variables
!!      refers to the number of grid cells on the current
!!      latitude circle.  in this stand-alone version, this
!!      number is arbitrarily set to three, to allow up to three
!!      simultaneous tests to be run.  the second dimension
!!      element of the "rot" variables refers to the maximum
!!      number of tiles in the mosaic.  in this stand-alone
!!      version, this number is set to eight.  the first
!!      dimension element in the "gat" variables is given by
!!      the product of the first two dimension elements in the
!!      "rot" variables.
!!
!!     The majority of CTEM parameters are stored in ctem_params.f90.
!!     Also the CTEM variables are stored in modules that we point to
!!     in this driver. We access the variables and parameters
!!     through use statements for modules:

      use ctem_params,        only : initpftpars,nlat,nmos,ilg,nmon,&
     &                               ican, ignd,icp1, icc, iccp1,&
     &                               monthend, mmday,modelpft, l2max,&
     &                                deltat, abszero, monthdays,seed,&
     &                                crop, NBS, lat, edgelat,earthrad,&
     &                                lon

      use landuse_change,     only : initialize_luc, readin_luc

      use ctem_statevars,     only : vrot,vgat,c_switch,initrowvars,&
     &                               class_out,resetclassmon,&
     &                               resetclassyr,&
     &                               resetmonthend,resetyearend,&
     &                               resetclassaccum,ctem_grd,&
     &                               ctem_tile,resetgridavg,&
     &                               finddaylength,alloc_ctem_vars,&
     &                               ctem_mo,ctem_grd_mo,ctem_tile_mo, &
     &                               ctem_yr,ctem_grd_yr,ctem_tile_yr

      use class_statevars,    only : alloc_class_vars,class_gat,class_rot

      use io_driver,          only : read_from_ctm, create_outfiles,&
     &                               write_ctm_rs, class_monthly_aw,&
     &                               ctem_annual_aw,ctem_monthly_aw,&
     &                               close_outfiles,ctem_daily_aw,&
     &                               class_annual_aw,bounds

      use model_state_drivers, only : read_initialstate

      implicit none

      ! Flag test
      real :: longitude, latitude

!
!     * INTEGER CONSTANTS.
!
      INTEGER IDISP  !<Flag governing treatment of vegetation displacement height
      INTEGER IZREF  !<Flag governing treatment of surface roughness length
      INTEGER ISLFD  !<Flag governing options for surface stability functions and diagnostic calculations
      INTEGER IPCP   !<Flag selecting algorithm for dividing precipitation between rainfall and snowfall
      INTEGER IWF    !<Flag governing lateral soil water flow calculations
      INTEGER IPAI   !<Flag to enable use of user-specified plant area index
      INTEGER IHGT   !<Flag to enable use of user-specified vegetation height
      INTEGER IALC   !<Flag to enable use of user-specified canopy albedo
      INTEGER IALS   !<Flag to enable use of user-specified snow albedo
      INTEGER IALG   !<Flag to enable use of user-specified ground albedo
      INTEGER ITG    !<Flag to select iteration scheme for ground or snow surface
      INTEGER ITC    !<Flag to select iteration scheme for canopy temperature
      INTEGER ITCG   !<Flag to select iteration scheme for surface under canopy
      INTEGER isnoalb!<
      INTEGER igralb !<

      INTEGER NLTEST  !<Number of grid cells being modelled for this run
      INTEGER NMTEST  !<Number of mosaic tiles per grid cell being modelled for this run
      INTEGER NCOUNT  !<Counter for daily averaging
      INTEGER NDAY    !<
      INTEGER IMONTH  !<
      INTEGER NDMONTH !<
      INTEGER NT      !<
      INTEGER IHOUR   !<Hour of day
      INTEGER IMIN    !<Minutes elapsed in current hour
      INTEGER IDAY    !<Julian day of the year
      INTEGER IYEAR   !<Year of run
      INTEGER NML     !<Counter representing number of mosaic tiles on modelled domain that are land
      INTEGER NMW     !<Counter representing number of mosaic tiles on modelled domain that are lakes
      INTEGER JLAT    !<Integer index corresponding to latitude of grid cell
      INTEGER NLANDCS !<Number of modelled areas that contain subareas of canopy over snow
      INTEGER NLANDGS !<Number of modelled areas that contain subareas of snow over bare ground
      INTEGER NLANDC  !<Number of modelled areas that contain subareas of canopy over bare ground
      INTEGER NLANDG  !<Number of modelled areas that contain subareas of bare ground
      INTEGER NLANDI  !<Number of modelled areas that are ice sheets
      INTEGER I,J,K,L,M,N
      INTEGER NTLD    !<
!
      INTEGER K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11
      INTEGER ITA        !<
      INTEGER ITCAN      !<
      INTEGER ITD        !<
      INTEGER ITAC       !<
      INTEGER ITS        !<
      INTEGER ITSCR      !<
      INTEGER ITD2       !<
      INTEGER ITD3       !<
      INTEGER ITD4       !<
      INTEGER NFS        !<
      INTEGER NDRY       !<

      INTEGER*4 TODAY(3), NOW(3)

! The following are stored in the data structure: class_gat
! they are allocatted in alloc_class_vars in the class_statevars
! module and pointed to here.

! These will be allocated the dimension: 'ilg'

    integer, pointer, dimension(:) :: ILMOS     !<Index of grid cell corresponding to current element of gathered vector of land surface variables [ ]
    integer, pointer, dimension(:) :: JLMOS     !<Index of mosaic tile corresponding to current element of gathered vector of land surface variables [ ]
    integer, pointer, dimension(:) :: IWMOS     !<Index of grid cell corresponding to current element of gathered vector of inland water body variables [ ]
    integer, pointer, dimension(:) :: JWMOS     !<Index of mosaic tile corresponding to current element of gathered vector of inland water body variables [ ]
    integer, pointer, dimension(:) :: IGDRGAT   !<Index of soil layer in which bedrock is encountered

    real, pointer, dimension(:) :: DELZ    !<
    real, pointer, dimension(:) :: ZBOT    !<
    real, pointer, dimension(:) :: ALBSGAT !<Snow albedo [ ]
    real, pointer, dimension(:) :: CMAIGAT !<Aggregated mass of vegetation canopy \f$[kg m^{-2} ]\f$
    real, pointer, dimension(:) :: GROGAT  !<Vegetation growth index [ ]
    real, pointer, dimension(:) :: QACGAT  !<Specific humidity of air within vegetation canopy space \f$[kg kg^{-1} ]\f$
    real, pointer, dimension(:) :: RCANGAT !<Intercepted liquid water stored on canopy \f$[kg m^{-2} ]\f$
    real, pointer, dimension(:) :: RHOSGAT !<Density of snow \f$[kg m^{-3} ]\f$
    real, pointer, dimension(:) :: SCANGAT !<Intercepted frozen water stored on canopy \f$[kg m^{-2} ]\f$
    real, pointer, dimension(:) :: SNOGAT  !<Mass of snow pack [kg m^{-2} ]\f$
    real, pointer, dimension(:) :: TACGAT  !<Temperature of air within vegetation canopy [K]
    real, pointer, dimension(:) :: TBASGAT !<Temperature of bedrock in third soil layer [K]
    real, pointer, dimension(:) :: TCANGAT !<Vegetation canopy temperature [K]
    real, pointer, dimension(:) :: TPNDGAT !<Temperature of ponded water [K]
    real, pointer, dimension(:) :: TSNOGAT !<Snowpack temperature [K]
    real, pointer, dimension(:) :: WSNOGAT !<Liquid water content of snow pack \f$[kg m^{-2} ]\f$
    real, pointer, dimension(:) :: ZPNDGAT !<Depth of ponded water on surface [m]
    real, pointer, dimension(:) :: REFGAT  !<
    real, pointer, dimension(:) :: BCSNGAT !<
    real, pointer, dimension(:) :: AGIDGAT !<Optional user-specified value of ground near-infrared albedo to override CLASS-calculated value [ ]
    real, pointer, dimension(:) :: AGVDGAT !<Optional user-specified value of ground visible albedo to override CLASS-calculated value [ ]
    real, pointer, dimension(:) :: ALGDGAT !<Reference albedo for dry soil [ ]
    real, pointer, dimension(:) :: ALGWGAT !<Reference albedo for saturated soil [ ]
    real, pointer, dimension(:) :: ASIDGAT !<Optional user-specified value of snow near-infrared albedo to override CLASS-calculated value [ ]
    real, pointer, dimension(:) :: ASVDGAT !<Optional user-specified value of snow visible albedo to override CLASS-calculated value [ ]
    real, pointer, dimension(:) :: DRNGAT  !<Drainage index at bottom of soil profile [ ]
    real, pointer, dimension(:) :: GRKFGAT !<WATROF parameter used when running MESH code [ ]
    real, pointer, dimension(:) :: WFCIGAT !<WATROF parameter used when running MESH code [ ]
    real, pointer, dimension(:) :: WFSFGAT !<WATROF parameter used when running MESH code [ ]
    real, pointer, dimension(:) :: XSLPGAT !<Surface slope (used when running MESH code) [degrees]
    real, pointer, dimension(:) :: ZPLGGAT !<Maximum water ponding depth for snow-free subareas (user-specified when running MESH code) [m]
    real, pointer, dimension(:) :: ZPLSGAT !<Maximum water ponding depth for snow-covered subareas (user-specified when running MESH code) [m]
    real, pointer, dimension(:) :: ZSNLGAT !<Limiting snow depth below which coverage is < 100% [m]
    real, pointer, dimension(:) :: ALGWVGAT !<
    real, pointer, dimension(:) :: ALGWNGAT !<
    real, pointer, dimension(:) :: ALGDVGAT !<
    real, pointer, dimension(:) :: ALGDNGAT !<
    real, pointer, dimension(:) :: EMISGAT  !<
    real, pointer, dimension(:) :: CSZGAT  !<Cosine of solar zenith angle [ ]
    real, pointer, dimension(:) :: DLONGAT !<Longitude of grid cell (east of Greenwich) [degrees]
    real, pointer, dimension(:) :: DLATGAT !< Latitude of grid cell [degrees]
    real, pointer, dimension(:) :: FCLOGAT !<Fractional cloud cover [ ]
    real, pointer, dimension(:) :: FDLGAT  !<Downwelling longwave radiation at bottom of atmosphere (i.e. incident on modelled land surface elements \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: FSIHGAT !<Near-infrared radiation incident on horizontal surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: FSVHGAT !<Visible radiation incident on horizontal surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: GGEOGAT !<Geothermal heat flux at bottom of soil profile \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: PADRGAT !<Partial pressure of dry air [Pa]
    real, pointer, dimension(:) :: PREGAT  !<Surface precipitation rate \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: PRESGAT !<Surface air pressure [Pa]
    real, pointer, dimension(:) :: QAGAT   !<Specific humidity at reference height \f$[kg kg^{-1} ]\f$
    real, pointer, dimension(:) :: RADJGAT !<Latitude of grid cell (positive north of equator) [rad]
    real, pointer, dimension(:) :: RHOAGAT !<Density of air \f$[kg m^{-3} ]\f$
    real, pointer, dimension(:) :: RHSIGAT !<Density of fresh snow \f$[kg m^{-3} ]\f$
    real, pointer, dimension(:) :: RPCPGAT !<Rainfall rate over modelled area \f$[m s^{-1} ]\f$
    real, pointer, dimension(:) :: SPCPGAT !<Snowfall rate over modelled area \f$[m s^{-1} ]\f$
    real, pointer, dimension(:) :: TAGAT   !<Air temperature at reference height [K]
    real, pointer, dimension(:) :: TADPGAT !<Dew point temperature of air [K]
    real, pointer, dimension(:) :: TRPCGAT !<Rainfall temperature [K]
    real, pointer, dimension(:) :: TSPCGAT !<Snowfall temperature [K]
    real, pointer, dimension(:) :: ULGAT   !<Zonal component of wind velocity \f$[m s^{-1} ]\f$
    real, pointer, dimension(:) :: VLGAT   !<Meridional component of wind velocity \f$[m s^{-1} ]\f$
    real, pointer, dimension(:) :: VMODGAT !<Wind speed at reference height \f$[m s^{-1} ]\f$
    real, pointer, dimension(:) :: VPDGAT  !<Vapour pressure deficit [mb]
    real, pointer, dimension(:) :: Z0ORGAT !<Orographic roughness length [m]
    real, pointer, dimension(:) :: ZBLDGAT !<Atmospheric blending height for surface roughness length averaging [m]
    real, pointer, dimension(:) :: ZDHGAT  !<User-specified height associated with diagnosed screen-level variables [m]
    real, pointer, dimension(:) :: ZDMGAT  !<User-specified height associated with diagnosed anemometer-level wind speed [m]
    real, pointer, dimension(:) :: ZRFHGAT !<Reference height associated with forcing air temperature and humidity [m]
    real, pointer, dimension(:) :: ZRFMGAT !<Reference height associated with forcing wind speed [m]
    real, pointer, dimension(:) :: FSGGAT  !<
    real, pointer, dimension(:) :: FLGGAT  !<
    real, pointer, dimension(:) :: GUSTGAT !<
    real, pointer, dimension(:) :: DEPBGAT !<
    real, pointer, dimension(:) :: GTBS    !<
    real, pointer, dimension(:) :: SFCUBS  !<
    real, pointer, dimension(:) :: SFCVBS  !<
    real, pointer, dimension(:) :: USTARBS !<
    real, pointer, dimension(:) :: TCSNOW  !<
    real, pointer, dimension(:) :: GSNOW   !<
    real, pointer, dimension(:) :: ALIRGAT !<Diagnosed total near-infrared albedo of land surface [ ]
    real, pointer, dimension(:) :: ALVSGAT !<Diagnosed total visible albedo of land surface [ ]
    real, pointer, dimension(:) :: CDHGAT  !<Surface drag coefficient for heat [ ]
    real, pointer, dimension(:) :: CDMGAT  !<Surface drag coefficient for momentum [ ]
    real, pointer, dimension(:) :: DRGAT   !<Surface drag coefficient under neutral stability [ ]
    real, pointer, dimension(:) :: EFGAT   !<Evaporation efficiency at ground surface [ ]
    real, pointer, dimension(:) :: FLGGGAT !<Diagnosed net longwave radiation at soil surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: FLGSGAT !<Diagnosed net longwave radiation at snow surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: FLGVGAT !<Diagnosed net longwave radiation on vegetation canopy \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: FSGGGAT !<Diagnosed net shortwave radiation at soil surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: FSGSGAT !<Diagnosed net shortwave radiation at snow surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: FSGVGAT !<Diagnosed net shortwave radiation on vegetation canopy \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: FSNOGAT !<Diagnosed fractional snow coverage [ ]
    real, pointer, dimension(:) :: GAGAT   !<Diagnosed product of drag coefficient and wind speed over modelled area \f$[m s^{-1} ]\f$
    real, pointer, dimension(:) :: GTGAT   !<Diagnosed effective surface black-body temperature [K]
    real, pointer, dimension(:) :: HBLGAT  !<Height of the atmospheric boundary layer [m]
    real, pointer, dimension(:) :: HEVCGAT !<Diagnosed latent heat flux on vegetation canopy \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: HEVGGAT !<Diagnosed latent heat flux at soil surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: HEVSGAT !<Diagnosed latent heat flux at snow surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: HFSGAT  !<Diagnosed total surface sensible heat flux over modelled area \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: HFSCGAT !<Diagnosed sensible heat flux on vegetation canopy \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: HFSGGAT !<Diagnosed sensible heat flux at soil surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: HFSSGAT !<Diagnosed sensible heat flux at snow surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: HMFCGAT !<Diagnosed energy associated with phase change of water on vegetation \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: HMFNGAT !<Diagnosed energy associated with phase change of water in snow pack \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: HTCCGAT !<Diagnosed internal energy change of vegetation canopy due to conduction and/or change in mass \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: HTCSGAT !<Diagnosed internal energy change of snow pack due to conduction and/or change in mass \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: ILMOGAT !<Inverse of Monin-Obukhov roughness length \f$(m^{-1} ]\f$
    real, pointer, dimension(:) :: PCFCGAT !<Diagnosed frozen precipitation intercepted by vegetation \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: PCLCGAT !<Diagnosed liquid precipitation intercepted by vegetation \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: PCPGGAT !<Diagnosed precipitation incident on ground \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: PCPNGAT !<Diagnosed precipitation incident on snow pack \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: PETGAT  !<Diagnosed potential evapotranspiration \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: QEVPGAT !<Diagnosed total surface latent heat flux over modelled area \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: QFCFGAT !<Diagnosed vapour flux from frozen water on vegetation \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: QFCLGAT !<Diagnosed vapour flux from liquid water on vegetation \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: QFGGAT  !<Diagnosed water vapour flux from ground \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: QFNGAT  !<Diagnosed water vapour flux from snow pack \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: QFSGAT  !<Diagnosed total surface water vapour flux over modelled area \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: QFXGAT  !<Product of surface drag coefficient, wind speed and surface-air specific humidity difference \f$[m s^{-1} ]\f$
    real, pointer, dimension(:) :: QGGAT   !<Diagnosed surface specific humidity \f$[kg kg^{-1} ]\f$
    real, pointer, dimension(:) :: ROFGAT  !<Total runoff from soil \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: ROFBGAT !<Base flow from bottom of soil column \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: ROFCGAT !<Liquid/frozen water runoff from vegetation \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: ROFNGAT !<Liquid water runoff from snow pack \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: ROFOGAT !<Overland flow from top of soil column \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: ROFSGAT !<Interflow from sides of soil column \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: ROVGGAT !<Diagnosed liquid/frozen water runoff from vegetation to ground surface \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: SFCQGAT !<Diagnosed screen-level specific humidity \f$[kg kg^{-1} ]\f$
    real, pointer, dimension(:) :: SFCTGAT !<Diagnosed screen-level air temperature [K]
    real, pointer, dimension(:) :: SFCUGAT !<Diagnosed anemometer-level zonal wind \f$[m s^{-1} ]\f$
    real, pointer, dimension(:) :: SFCVGAT !<Diagnosed anemometer-level meridional wind \f$[m s^{-1} ]\f$
    real, pointer, dimension(:) :: TFXGAT  !<Product of surface drag coefficient, wind speed and surface-air temperature difference \f$[K m s^{-1} ]\f$
    real, pointer, dimension(:) :: TROBGAT !<Temperature of base flow from bottom of soil column [K]
    real, pointer, dimension(:) :: TROFGAT !<Temperature of total runoff [K]
    real, pointer, dimension(:) :: TROOGAT !<Temperature of overland flow from top of soil column [K]
    real, pointer, dimension(:) :: TROSGAT !<Temperature of interflow from sides of soil column [K]
    real, pointer, dimension(:) :: UEGAT   !<Friction velocity of air \f$[m s^{-1} ]\f$
    real, pointer, dimension(:) :: WTABGAT !<Depth of water table in soil [m]
    real, pointer, dimension(:) :: WTRCGAT !<Diagnosed residual water transferred off the vegetation canopy \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: WTRGGAT !<Diagnosed residual water transferred into or out of the soil \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: WTRSGAT !<Diagnosed residual water transferred into or out of the snow pack \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: QLWOGAT !<
    real, pointer, dimension(:) :: SFRHGAT !<
    real, pointer, dimension(:) :: FTEMP   !<
    real, pointer, dimension(:) :: FVAP    !<
    real, pointer, dimension(:) :: RIB     !<
    real, pointer, dimension(:) :: FC      !<
    real, pointer, dimension(:) :: FG      !<
    real, pointer, dimension(:) :: FCS     !<
    real, pointer, dimension(:) :: FGS     !<
    real, pointer, dimension(:) :: RBCOEF  !<
    real, pointer, dimension(:) :: ZSNOW   !<
    real, pointer, dimension(:) :: FSVF    !<
    real, pointer, dimension(:) :: FSVFS   !<
    real, pointer, dimension(:) :: ALVSCN  !<
    real, pointer, dimension(:) :: ALIRCN  !<
    real, pointer, dimension(:) :: ALVSG   !<
    real, pointer, dimension(:) :: ALIRG   !<
    real, pointer, dimension(:) :: ALVSCS  !<
    real, pointer, dimension(:) :: ALIRCS  !<
    real, pointer, dimension(:) :: ALVSSN  !<
    real, pointer, dimension(:) :: ALIRSN  !<
    real, pointer, dimension(:) :: ALVSGC  !<
    real, pointer, dimension(:) :: ALIRGC  !<
    real, pointer, dimension(:) :: ALVSSC  !<
    real, pointer, dimension(:) :: ALIRSC  !<
    real, pointer, dimension(:) :: TRVSCN  !<
    real, pointer, dimension(:) :: TRIRCN  !<
    real, pointer, dimension(:) :: TRVSCS  !<
    real, pointer, dimension(:) :: TRIRCS  !<
    real, pointer, dimension(:) :: RC      !<
    real, pointer, dimension(:) :: RCS     !<
    real, pointer, dimension(:) :: FRAINC  !<
    real, pointer, dimension(:) :: FSNOWC  !<
    real, pointer, dimension(:) :: FRAICS  !<
    real, pointer, dimension(:) :: FSNOCS  !<
    real, pointer, dimension(:) :: CMASSC  !<
    real, pointer, dimension(:) :: CMASCS  !<
    real, pointer, dimension(:) :: DISP    !<
    real, pointer, dimension(:) :: DISPS   !<
    real, pointer, dimension(:) :: ZOMLNC  !<
    real, pointer, dimension(:) :: ZOELNC  !<
    real, pointer, dimension(:) :: ZOMLNG  !<
    real, pointer, dimension(:) :: ZOELNG  !<
    real, pointer, dimension(:) :: ZOMLCS  !<
    real, pointer, dimension(:) :: ZOELCS  !<
    real, pointer, dimension(:) :: ZOMLNS  !<
    real, pointer, dimension(:) :: ZOELNS  !<
    real, pointer, dimension(:) :: TRSNOWC !<
    real, pointer, dimension(:) :: CHCAP   !<
    real, pointer, dimension(:) :: CHCAPS  !<
    real, pointer, dimension(:) :: GZEROC  !<
    real, pointer, dimension(:) :: GZEROG  !<
    real, pointer, dimension(:) :: GZROCS  !<
    real, pointer, dimension(:) :: GZROGS  !<
    real, pointer, dimension(:) :: G12C    !<
    real, pointer, dimension(:) :: G12G    !<
    real, pointer, dimension(:) :: G12CS   !<
    real, pointer, dimension(:) :: G12GS   !<
    real, pointer, dimension(:) :: G23C    !<
    real, pointer, dimension(:) :: G23G    !<
    real, pointer, dimension(:) :: G23CS   !<
    real, pointer, dimension(:) :: G23GS   !<
    real, pointer, dimension(:) :: QFREZC  !<
    real, pointer, dimension(:) :: QFREZG  !<
    real, pointer, dimension(:) :: QMELTC  !<
    real, pointer, dimension(:) :: QMELTG  !<
    real, pointer, dimension(:) :: EVAPC   !<
    real, pointer, dimension(:) :: EVAPCG  !<
    real, pointer, dimension(:) :: EVAPG   !<
    real, pointer, dimension(:) :: EVAPCS  !<
    real, pointer, dimension(:) :: EVPCSG  !<
    real, pointer, dimension(:) :: EVAPGS  !<
    real, pointer, dimension(:) :: TCANO   !<
    real, pointer, dimension(:) :: TCANS   !<
    real, pointer, dimension(:) :: RAICAN  !<
    real, pointer, dimension(:) :: SNOCAN  !<
    real, pointer, dimension(:) :: RAICNS  !<
    real, pointer, dimension(:) :: SNOCNS  !<
    real, pointer, dimension(:) :: CWLCAP  !<
    real, pointer, dimension(:) :: CWFCAP  !<
    real, pointer, dimension(:) :: CWLCPS  !<
    real, pointer, dimension(:) :: CWFCPS  !<
    real, pointer, dimension(:) :: TSNOCS  !<
    real, pointer, dimension(:) :: TSNOGS  !<
    real, pointer, dimension(:) :: RHOSCS  !<
    real, pointer, dimension(:) :: RHOSGS  !<
    real, pointer, dimension(:) :: WSNOCS  !<
    real, pointer, dimension(:) :: WSNOGS  !<
    real, pointer, dimension(:) :: TPONDC  !<
    real, pointer, dimension(:) :: TPONDG  !<
    real, pointer, dimension(:) :: TPNDCS  !<
    real, pointer, dimension(:) :: TPNDGS  !<
    real, pointer, dimension(:) :: ZPLMCS  !<
    real, pointer, dimension(:) :: ZPLMGS  !<
    real, pointer, dimension(:) :: ZPLIMC  !<
    real, pointer, dimension(:) :: ZPLIMG  !<
    !
    !     * DIAGNOSTIC ARRAYS USED FOR CHECKING ENERGY AND WATER
    !     * BALANCES.
    !
    real, pointer, dimension(:) :: CTVSTP !<
    real, pointer, dimension(:) :: CTSSTP !<
    real, pointer, dimension(:) :: CT1STP !<
    real, pointer, dimension(:) :: CT2STP !<
    real, pointer, dimension(:) :: CT3STP !<
    real, pointer, dimension(:) :: WTVSTP !<
    real, pointer, dimension(:) :: WTSSTP !<
    real, pointer, dimension(:) :: WTGSTP !<

! These will be allocated the dimension: 'ilg, ignd'
    integer, pointer, dimension(:,:) :: ISNDGAT !<Integer identifier associated with sand content
    real, pointer, dimension(:,:) :: TBARGAT !<Temperature of soil layers [K]
    real, pointer, dimension(:,:) :: THICGAT !<Volumetric frozen water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, pointer, dimension(:,:) :: THLQGAT !<Volumetric liquid water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, pointer, dimension(:,:) :: BIGAT   !<Clapp and Hornberger empirical “b” parameter [ ]
    real, pointer, dimension(:,:) :: DLZWGAT !<Permeable thickness of soil layer [m]
    real, pointer, dimension(:,:) :: GRKSGAT !<Saturated hydraulic conductivity of soil layers \f$[m s^{-1} ]\f$
    real, pointer, dimension(:,:) :: HCPSGAT !<Volumetric heat capacity of soil particles \f$[J m^{-3} ]\f$
    real, pointer, dimension(:,:) :: PSISGAT !<Soil moisture suction at saturation [m]
    real, pointer, dimension(:,:) :: PSIWGAT !<Soil moisture suction at wilting point [m]
    real, pointer, dimension(:,:) :: TCSGAT  !<Thermal conductivity of soil particles \f$[W m^{-1} K^{-1} ]\f$\
    real, pointer, dimension(:,:) :: THFCGAT !<Field capacity \f$[m^3 m^{-3} ]\f$
    real, pointer, dimension(:,:) :: THMGAT  !<Residual soil liquid water content remaining after freezing or evaporation \f$[m^3 m^{-3} ]\f$
    real, pointer, dimension(:,:) :: THPGAT  !<Pore volume in soil layer \f$[m^3 m^{-3} ]\f$
    real, pointer, dimension(:,:) :: THRGAT  !<Liquid water retention capacity for organic soil \f$[m^3 m^{-3} ]\f$
    real, pointer, dimension(:,:) :: THRAGAT !<Fractional saturation of soil behind the wetting front [ ]
    real, pointer, dimension(:,:) :: ZBTWGAT !<Depth to permeable bottom of soil layer [m]
    real, pointer, dimension(:,:) :: THLWGAT !<
    real, pointer, dimension(:,:) :: GFLXGAT !<Heat conduction between soil layers \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:) :: HMFGGAT !<Diagnosed energy associated with phase change of water in soil layers \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:) :: HTCGAT  !<Diagnosed internal energy change of soil layer due to conduction and/or change in mass \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:) :: QFCGAT  !<Diagnosed vapour flux from transpiration over modelled area \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:) :: TBARC  !<
    real, pointer, dimension(:,:) :: TBARG  !<
    real, pointer, dimension(:,:) :: TBARCS !<
    real, pointer, dimension(:,:) :: TBARGS !<
    real, pointer, dimension(:,:) :: THLIQC !<
    real, pointer, dimension(:,:) :: THLIQG !<
    real, pointer, dimension(:,:) :: THICEC !<
    real, pointer, dimension(:,:) :: THICEG !<
    real, pointer, dimension(:,:) :: FROOT  !<
    real, pointer, dimension(:,:) :: HCPC   !<
    real, pointer, dimension(:,:) :: HCPG   !<
    real, pointer, dimension(:,:) :: FROOTS !<
    real, pointer, dimension(:,:) :: TCTOPC !<
    real, pointer, dimension(:,:) :: TCBOTC !<
    real, pointer, dimension(:,:) :: TCTOPG !<
    real, pointer, dimension(:,:) :: TCBOTG !<

! These will be allocated the dimension: 'ilg, ican'
    real, pointer, dimension(:,:) :: ACIDGAT !<Optional user-specified value of canopy near-infrared albedo to override CLASS-calculated value [ ]
    real, pointer, dimension(:,:) :: ACVDGAT !<Optional user-specified value of canopy visible albedo to override CLASS-calculated value [ ]
    real, pointer, dimension(:,:) :: CMASGAT !<Maximum canopy mass for vegetation category \f$[kg m^{-2} ]\f$
    real, pointer, dimension(:,:) :: HGTDGAT !<Optional user-specified values of height of vegetation categories to override CLASS-calculated values [m]
    real, pointer, dimension(:,:) :: PAIDGAT !<Optional user-specified value of plant area indices of vegetation categories to override CLASS-calculated values [ ]
    real, pointer, dimension(:,:) :: PAMNGAT !<Minimum plant area index of vegetation category [ ]
    real, pointer, dimension(:,:) :: PAMXGAT !<Minimum plant area index of vegetation category [ ]
    real, pointer, dimension(:,:) :: PSGAGAT !<Soil moisture suction coefficient for vegetation category (used in stomatal resistance calculation) [ ]
    real, pointer, dimension(:,:) :: PSGBGAT !<Soil moisture suction coefficient for vegetation category (used in stomatal resistance calculation) [ ]
    real, pointer, dimension(:,:) :: QA50GAT !<Reference value of incoming shortwave radiation for vegetation category (used in stomatal resistance calculation) \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:) :: ROOTGAT !<Maximum rooting depth of vegetation category [m]
    real, pointer, dimension(:,:) :: RSMNGAT !<Minimum stomatal resistance of vegetation category \f$[s m^{-1} ]\f$
    real, pointer, dimension(:,:) :: VPDAGAT !<Vapour pressure deficit coefficient for vegetation category (used in stomatal resistance calculation) [ ]
    real, pointer, dimension(:,:) :: VPDBGAT !<Vapour pressure deficit coefficient for vegetation category (used in stomatal resistance calculation) [ ]


! These will be allocated the dimension: 'ilg, icp1'
    real, pointer, dimension(:,:) :: ALICGAT !<Background average near-infrared albedo of vegetation category [ ]
    real, pointer, dimension(:,:) :: ALVCGAT !<Background average visible albedo of vegetation category [ ]
    real, pointer, dimension(:,:) :: FCANGAT !<Maximum fractional coverage of modelled area by vegetation category [ ]
    real, pointer, dimension(:,:) :: LNZ0GAT !<Natural logarithm of maximum roughness length of vegetation category [ ]

! These will be allocated the dimension: 'ilg, nbs'
    real, pointer, dimension(:,:) :: FSDBGAT !<
    real, pointer, dimension(:,:) :: FSFBGAT !<
    real, pointer, dimension(:,:) :: FSSBGAT !<
    real, pointer, dimension(:,:) :: SALBGAT !<
    real, pointer, dimension(:,:) :: CSALGAT !<
    real, pointer, dimension(:,:) :: ALTG    !<
    real, pointer, dimension(:,:) :: ALSNO   !<
    real, pointer, dimension(:,:) :: TRSNOWG !<

! These will be allocated the dimension: 'ilg, 4'
    real, pointer, dimension(:,:) :: TSFSGAT !<Ground surface temperature over subarea [K]

! These will be allocated the dimension: 'ilg, 6, 50'
    integer, pointer, dimension(:,:,:) :: ITCTGAT !<Counter of number of iterations required to solve surface energy balance for the elements of the four subareas

! The following are stored in the data structure: class_rot
! they are allocatted in alloc_class_vars in the class_statevars
! module and pointed to here.

! These will be allocated the dimension: 'nlat'

    real, pointer, dimension(:) :: ALIRACC !<Diagnosed total near-infrared albedo of land surface [ ]
    real, pointer, dimension(:) :: ALVSACC !<Diagnosed total visible albedo of land surface [ ]
    real, pointer, dimension(:) :: EVAPACC !<Diagnosed total surface water vapour flux over modelled area \f$[kg m^{-2} ]\f$
    real, pointer, dimension(:) :: FLINACC !<Downwelling longwave radiation above surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: FLUTACC !<Upwelling longwave radiation from surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: FSINACC !<Downwelling shortwave radiation above surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: GROACC  !<Vegetation growth index [ ]
    real, pointer, dimension(:) :: GTACC   !<Diagnosed effective surface black-body temperature [K]
    real, pointer, dimension(:) :: HFSACC  !<Diagnosed total surface sensible heat flux over modelled area \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: HMFNACC !<Diagnosed energy associated with phase change of water in snow pack \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: OVRACC  !<Overland flow from top of soil column \f$[kg m^{-2} ]\f$
    real, pointer, dimension(:) :: PREACC  !<Surface precipitation rate \f$[kg m^{-2} ]\f$
    real, pointer, dimension(:) :: PRESACC !<Surface air pressure [Pa]
    real, pointer, dimension(:) :: QAACC   !<Specific humidity at reference height \f$[kg kg^{-1} ]\f$
    real, pointer, dimension(:) :: QEVPACC !<Diagnosed total surface latent heat flux over modelled area \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: RCANACC !<Intercepted liquid water stored on canopy \f$[kg m^{-2} ]\f$
    real, pointer, dimension(:) :: RHOSACC !<Density of snow \f$[kg m^{-3} ]\f$
    real, pointer, dimension(:) :: ROFACC  !<Total runoff from soil \f$[kg m^{-2} ]\f$
    real, pointer, dimension(:) :: SCANACC !<Intercepted frozen water stored on canopy \f$[kg m^{-2} ]\f$
    real, pointer, dimension(:) :: SNOACC  !<Mass of snow pack \f$[kg m^{-2} ]\f$
    real, pointer, dimension(:) :: TAACC   !<Air temperature at reference height [K]
    real, pointer, dimension(:) :: TCANACC !<Vegetation canopy temperature [K]
    real, pointer, dimension(:) :: TSNOACC !<Snowpack temperature [K]
    real, pointer, dimension(:) :: UVACC   !<Wind speed \f$[m s^{-1} ]\f$
    real, pointer, dimension(:) :: WSNOACC !<Liquid water content of snow pack \f$[kg m^{-2} ]\f$
    real, pointer, dimension(:) :: WTBLACC !<Depth of water table in soil [m]
    real, pointer, dimension(:) :: ALTOTACC!<Broadband albedo [-]
    real, pointer, dimension(:) :: CANARE  !<
    real, pointer, dimension(:) :: SNOARE  !<
    real, pointer, dimension(:) :: CSZROW  !<
    real, pointer, dimension(:) :: DLONROW !<
    real, pointer, dimension(:) :: DLATROW !<
    real, pointer, dimension(:) :: FCLOROW !<
    real, pointer, dimension(:) :: FDLROW  !<
    real, pointer, dimension(:) :: FSIHROW !<
    real, pointer, dimension(:) :: FSVHROW !<
    real, pointer, dimension(:) :: GCROW   !<Type identifier for grid cell (1 = sea ice, 0 = ocean, -1 = land)
    real, pointer, dimension(:) :: GGEOROW !<
    real, pointer, dimension(:) :: PADRROW !<
    real, pointer, dimension(:) :: PREROW  !<
    real, pointer, dimension(:) :: PRESROW !<
    real, pointer, dimension(:) :: QAROW   !<
    real, pointer, dimension(:) :: RADJROW !<
    real, pointer, dimension(:) :: RHOAROW !<
    real, pointer, dimension(:) :: RHSIROW !<
    real, pointer, dimension(:) :: RPCPROW !<
    real, pointer, dimension(:) :: RPREROW !<Rainfall rate over modelled area \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: SPCPROW !<
    real, pointer, dimension(:) :: SPREROW !<Snowfall rate over modelled area \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: TAROW   !<
    real, pointer, dimension(:) :: TADPROW !<
    real, pointer, dimension(:) :: TRPCROW !<
    real, pointer, dimension(:) :: TSPCROW !<
    real, pointer, dimension(:) :: ULROW   !<
    real, pointer, dimension(:) :: VLROW   !<
    real, pointer, dimension(:) :: VMODROW !<
    real, pointer, dimension(:) :: VPDROW  !<
    real, pointer, dimension(:) :: ZBLDROW !<
    real, pointer, dimension(:) :: ZDHROW  !<
    real, pointer, dimension(:) :: ZDMROW  !<
    real, pointer, dimension(:) :: ZRFHROW !<
    real, pointer, dimension(:) :: ZRFMROW !<
    real, pointer, dimension(:) :: UVROW   !<
    real, pointer, dimension(:) :: XDIFFUS !<
    real, pointer, dimension(:) :: Z0ORROW !<
    real, pointer, dimension(:) :: FSSROW  !< Shortwave radiation \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: PRENROW !<
    real, pointer, dimension(:) :: CLDTROW !<
    real, pointer, dimension(:) :: FSGROL  !<
    real, pointer, dimension(:) :: FLGROL  !<
    real, pointer, dimension(:) :: GUSTROL !<
    real, pointer, dimension(:) :: DEPBROW !<
    real, pointer, dimension(:) :: ALIRROW !<
    real, pointer, dimension(:) :: ALVSROW !<
    real, pointer, dimension(:) :: CDHROW  !<
    real, pointer, dimension(:) :: CDMROW  !<
    real, pointer, dimension(:) :: DRROW   !<
    real, pointer, dimension(:) :: EFROW   !<
    real, pointer, dimension(:) :: FLGGROW !<
    real, pointer, dimension(:) :: FLGSROW !<
    real, pointer, dimension(:) :: FLGVROW !<
    real, pointer, dimension(:) :: FSGGROW !<
    real, pointer, dimension(:) :: FSGSROW !<
    real, pointer, dimension(:) :: FSGVROW !<
    real, pointer, dimension(:) :: FSNOROW !<
    real, pointer, dimension(:) :: GAROW   !<
    real, pointer, dimension(:) :: GTROW   !<
    real, pointer, dimension(:) :: HBLROW  !<
    real, pointer, dimension(:) :: HEVCROW !<
    real, pointer, dimension(:) :: HEVGROW !<
    real, pointer, dimension(:) :: HEVSROW !<
    real, pointer, dimension(:) :: HFSROW  !<
    real, pointer, dimension(:) :: HFSCROW !<
    real, pointer, dimension(:) :: HFSGROW !<
    real, pointer, dimension(:) :: HFSSROW !<
    real, pointer, dimension(:) :: HMFCROW !<
    real, pointer, dimension(:) :: HMFNROW !<
    real, pointer, dimension(:) :: HTCCROW !<
    real, pointer, dimension(:) :: HTCSROW !<
    real, pointer, dimension(:) :: ILMOROW !<
    real, pointer, dimension(:) :: PCFCROW !<
    real, pointer, dimension(:) :: PCLCROW !<
    real, pointer, dimension(:) :: PCPGROW !<
    real, pointer, dimension(:) :: PCPNROW !<
    real, pointer, dimension(:) :: PETROW  !<
    real, pointer, dimension(:) :: QEVPROW !<
    real, pointer, dimension(:) :: QFCFROW !<
    real, pointer, dimension(:) :: QFCLROW !<
    real, pointer, dimension(:) :: QFGROW  !<
    real, pointer, dimension(:) :: QFNROW  !<
    real, pointer, dimension(:) :: QFSROW  !<
    real, pointer, dimension(:) :: QFXROW  !<
    real, pointer, dimension(:) :: QGROW   !<
    real, pointer, dimension(:) :: ROFROW  !<
    real, pointer, dimension(:) :: ROFBROW !<
    real, pointer, dimension(:) :: ROFCROW !<
    real, pointer, dimension(:) :: ROFNROW !<
    real, pointer, dimension(:) :: ROFOROW !<
    real, pointer, dimension(:) :: ROFSROW !<
    real, pointer, dimension(:) :: ROVGROW !<
    real, pointer, dimension(:) :: SFCQROW !<
    real, pointer, dimension(:) :: SFCTROW !<
    real, pointer, dimension(:) :: SFCUROW !<
    real, pointer, dimension(:) :: SFCVROW !<
    real, pointer, dimension(:) :: TFXROW  !<
    real, pointer, dimension(:) :: UEROW   !<
    real, pointer, dimension(:) :: WTABROW !<
    real, pointer, dimension(:) :: WTRCROW !<
    real, pointer, dimension(:) :: WTRGROW !<
    real, pointer, dimension(:) :: WTRSROW !<
    real, pointer, dimension(:) :: SFRHROW !<

    ! These will be allocated the dimension: 'nlat,nmos'

    integer, pointer, dimension(:,:) :: IGDRROT !<
    integer, pointer, dimension(:,:) :: MIDROT  !<Mosaic tile type identifier (1 for land surface, 0 for inland lake)
    real, pointer, DIMENSION(:,:) :: ALBSROT !<
    real, pointer, dimension(:,:) :: CMAIROT !<
    real, pointer, dimension(:,:) :: GROROT  !<
    real, pointer, dimension(:,:) :: QACROT  !<
    real, pointer, dimension(:,:) :: RCANROT !<
    real, pointer, dimension(:,:) :: RHOSROT !<
    real, pointer, dimension(:,:) :: SCANROT !<
    real, pointer, dimension(:,:) :: SNOROT  !<
    real, pointer, dimension(:,:) :: TACROT  !<
    real, pointer, dimension(:,:) :: TBASROT !<
    real, pointer, dimension(:,:) :: TCANROT !<
    real, pointer, dimension(:,:) :: TPNDROT !<
    real, pointer, dimension(:,:) :: TSNOROT !<
    real, pointer, dimension(:,:) :: WSNOROT !<
    real, pointer, dimension(:,:) :: ZPNDROT !<
    real, pointer, dimension(:,:) :: REFROT  !<
    real, pointer, dimension(:,:) :: BCSNROT !<
    real, pointer, dimension(:,:) :: AGIDROT !<
    real, pointer, dimension(:,:) :: AGVDROT !<
    real, pointer, dimension(:,:) :: ALGDROT !<
    real, pointer, dimension(:,:) :: ALGWROT !<
    real, pointer, dimension(:,:) :: ASIDROT !<
    real, pointer, dimension(:,:) :: ASVDROT !<
    real, pointer, dimension(:,:) :: DRNROT  !<
    real, pointer, dimension(:,:) :: FAREROT !<Fractional coverage of mosaic tile on modelled area
    real, pointer, dimension(:,:) :: GRKFROT !<
    real, pointer, dimension(:,:) :: WFCIROT !<
    real, pointer, dimension(:,:) :: WFSFROT !<
    real, pointer, dimension(:,:) :: XSLPROT !<
    real, pointer, dimension(:,:) :: ZPLGROT !<
    real, pointer, dimension(:,:) :: ZPLSROT !<
    real, pointer, dimension(:,:) :: ZSNLROT !<
    real, pointer, dimension(:,:) :: ZSNOROT  !<
    real, pointer, dimension(:,:) :: ALGWVROT !<
    real, pointer, dimension(:,:) :: ALGWNROT !<
    real, pointer, dimension(:,:) :: ALGDVROT !<
    real, pointer, dimension(:,:) :: ALGDNROT !<
    real, pointer, dimension(:,:) :: EMISROT  !<
    real, pointer, dimension(:,:) :: ALIRROT !<
    real, pointer, dimension(:,:) :: ALVSROT !<
    real, pointer, dimension(:,:) :: CDHROT  !<
    real, pointer, dimension(:,:) :: CDMROT  !<
    real, pointer, dimension(:,:) :: DRROT   !<
    real, pointer, dimension(:,:) :: EFROT   !<
    real, pointer, dimension(:,:) :: FLGGROT !<
    real, pointer, dimension(:,:) :: FLGSROT !<
    real, pointer, dimension(:,:) :: FLGVROT !<
    real, pointer, dimension(:,:) :: FSGGROT !<
    real, pointer, dimension(:,:) :: FSGSROT !<
    real, pointer, dimension(:,:) :: FSGVROT !<
    real, pointer, dimension(:,:) :: FSNOROT !<
    real, pointer, dimension(:,:) :: GAROT   !<
    real, pointer, dimension(:,:) :: GTROT   !<
    real, pointer, dimension(:,:) :: HBLROT  !<
    real, pointer, dimension(:,:) :: HEVCROT !<
    real, pointer, dimension(:,:) :: HEVGROT !<
    real, pointer, dimension(:,:) :: HEVSROT !<
    real, pointer, dimension(:,:) :: HFSROT  !<
    real, pointer, dimension(:,:) :: HFSCROT !<
    real, pointer, dimension(:,:) :: HFSGROT !<
    real, pointer, dimension(:,:) :: HFSSROT !<
    real, pointer, dimension(:,:) :: HMFCROT !<
    real, pointer, dimension(:,:) :: HMFNROT !<
    real, pointer, dimension(:,:) :: HTCCROT !<
    real, pointer, dimension(:,:) :: SDEPROT !<Depth to bedrock in the soil profile
    real, pointer, dimension(:,:) :: SOCIROT  !<
    real, pointer, dimension(:,:) :: HTCSROT !<
    real, pointer, dimension(:,:) :: ILMOROT !<
    real, pointer, dimension(:,:) :: PCFCROT !<
    real, pointer, dimension(:,:) :: PCLCROT !<
    real, pointer, dimension(:,:) :: PCPGROT !<
    real, pointer, dimension(:,:) :: PCPNROT !<
    real, pointer, dimension(:,:) :: PETROT  !<
    real, pointer, dimension(:,:) :: QEVPROT !<
    real, pointer, dimension(:,:) :: QFCFROT !<
    real, pointer, dimension(:,:) :: QFCLROT !<
    real, pointer, dimension(:,:) :: QFGROT  !<
    real, pointer, dimension(:,:) :: QFNROT  !<
    real, pointer, dimension(:,:) :: QFSROT  !<
    real, pointer, dimension(:,:) :: QFXROT  !<
    real, pointer, dimension(:,:) :: QGROT   !<
    real, pointer, dimension(:,:) :: ROFROT  !<
    real, pointer, dimension(:,:) :: ROFBROT !<
    real, pointer, dimension(:,:) :: ROFCROT !<
    real, pointer, dimension(:,:) :: ROFNROT !<
    real, pointer, dimension(:,:) :: ROFOROT !<
    real, pointer, dimension(:,:) :: ROFSROT !<
    real, pointer, dimension(:,:) :: ROVGROT !<
    real, pointer, dimension(:,:) :: SFCQROT !<
    real, pointer, dimension(:,:) :: SFCTROT !<
    real, pointer, dimension(:,:) :: SFCUROT !<
    real, pointer, dimension(:,:) :: SFCVROT !<
    real, pointer, dimension(:,:) :: TFXROT  !<
    real, pointer, dimension(:,:) :: TROBROT !<
    real, pointer, dimension(:,:) :: TROFROT !<
    real, pointer, dimension(:,:) :: TROOROT !<
    real, pointer, dimension(:,:) :: TROSROT !<
    real, pointer, dimension(:,:) :: UEROT   !<
    real, pointer, dimension(:,:) :: WTABROT !<
    real, pointer, dimension(:,:) :: WTRCROT !<
    real, pointer, dimension(:,:) :: WTRGROT !<
    real, pointer, dimension(:,:) :: WTRSROT !<
    real, pointer, dimension(:,:) :: SFRHROT !<

    ! There will be allocated the dimension: 'nlat,nmos,ignd'
    integer, pointer, dimension(:,:,:) :: ISNDROT !<
    real, pointer, dimension(:,:,:) :: TBARROT !<
    real, pointer, dimension(:,:,:) :: THICROT !<
    real, pointer, dimension(:,:,:) :: THLQROT !<
    real, pointer, dimension(:,:,:) :: BIROT   !<
    real, pointer, dimension(:,:,:) :: DLZWROT !<
    real, pointer, dimension(:,:,:) :: GRKSROT !<
    real, pointer, dimension(:,:,:) :: HCPSROT !<
    real, pointer, dimension(:,:,:) :: SANDROT !<Percentage sand content of soil
    real, pointer, dimension(:,:,:) :: CLAYROT !<Percentage clay content of soil
    real, pointer, dimension(:,:,:) :: ORGMROT !<Percentage organic matter content of soil
    real, pointer, dimension(:,:,:) :: PSISROT !<
    real, pointer, dimension(:,:,:) :: PSIWROT !<
    real, pointer, dimension(:,:,:) :: TCSROT  !<
    real, pointer, dimension(:,:,:) :: THFCROT !<
    real, pointer, dimension(:,:,:) :: THMROT  !<
    real, pointer, dimension(:,:,:) :: THPROT  !<
    real, pointer, dimension(:,:,:) :: THRROT  !<
    real, pointer, dimension(:,:,:) :: THRAROT !<
    real, pointer, dimension(:,:,:) :: ZBTWROT !<
    real, pointer, dimension(:,:,:) :: THLWROT !<
    real, pointer, dimension(:,:,:) :: GFLXROT !<
    real, pointer, dimension(:,:,:) :: HMFGROT !<
    real, pointer, dimension(:,:,:) :: HTCROT  !<
    real, pointer, dimension(:,:,:) :: QFCROT  !<

    ! These will be allocated the dimension: 'nlat,nmos,ican'
    real, pointer, dimension(:,:,:) :: ACIDROT !<
    real, pointer, dimension(:,:,:) :: ACVDROT !<
    real, pointer, dimension(:,:,:) :: CMASROT !<
    real, pointer, dimension(:,:,:) :: HGTDROT !<
    real, pointer, dimension(:,:,:) :: PAIDROT !<
    real, pointer, dimension(:,:,:) :: PAMNROT !<
    real, pointer, dimension(:,:,:) :: PAMXROT !<
    real, pointer, dimension(:,:,:) :: PSGAROT !<
    real, pointer, dimension(:,:,:) :: PSGBROT !<
    real, pointer, dimension(:,:,:) :: QA50ROT !<
    real, pointer, dimension(:,:,:) :: ROOTROT !<
    real, pointer, dimension(:,:,:) :: RSMNROT !<
    real, pointer, dimension(:,:,:) :: VPDAROT !<
    real, pointer, dimension(:,:,:) :: VPDBROT !<

    ! These will be allocated the dimension: 'nlat,nmos,icp1'
    real, pointer, dimension(:,:,:) :: ALICROT !<
    real, pointer, dimension(:,:,:) :: ALVCROT !<
    real, pointer, dimension(:,:,:) :: FCANROT !<
    real, pointer, dimension(:,:,:) :: LNZ0ROT !<

    ! These will be allocated the dimension: 'nlat,nmos,nbs'
    real, pointer, dimension(:,:,:)  :: SALBROT  !<
    real, pointer, dimension(:,:,:)  :: CSALROT  !<

    ! These will be allocated the dimension: 'nlat,nbs'
    real, pointer, dimension(:,:) :: FSDBROL  !<
    real, pointer, dimension(:,:) :: FSFBROL  !<
    real, pointer, dimension(:,:) :: FSSBROL  !<

    ! These will be allocated the dimension: 'nlat,ignd'

    real, pointer, dimension(:,:) :: TBARACC !<Temperature of soil layers [K]
    real, pointer, dimension(:,:) :: THALACC !<Total volumetric water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, pointer, dimension(:,:) :: THICACC !<Volumetric frozen water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, pointer, dimension(:,:) :: THLQACC !<Volumetric liquid water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, pointer, dimension(:,:) :: GFLXROW !<
    real, pointer, dimension(:,:) :: HMFGROW !<
    real, pointer, dimension(:,:) :: HTCROW  !<
    real, pointer, dimension(:,:) :: QFCROW  !<

    ! These will be allocated the dimension: 'nlat,nmos,6,50'
    integer, pointer, dimension(:,:,:,:) :: ITCTROT !<

    ! These will be allocated the dimension: 'nlat,nmos,4'
    real, pointer, dimension(:,:,:)  :: TSFSROT !<

!       REAL,DIMENSION(ILG)            :: ALBSGAT !<Snow albedo [ ]
!       REAL,DIMENSION(NLAT,NMOS)      :: ALBSROT !<
!       REAL,DIMENSION(ILG)            :: CMAIGAT !<Aggregated mass of vegetation canopy \f$[kg m^{-2} ]\f$
!       REAL,DIMENSION(NLAT,NMOS)      :: CMAIROT !<
!       REAL,DIMENSION(ILG)            :: GROGAT  !<Vegetation growth index [ ]
!       REAL,DIMENSION(NLAT,NMOS)      :: GROROT  !<
!       REAL,DIMENSION(ILG)            :: QACGAT  !<Specific humidity of air within vegetation canopy space \f$[kg kg^{-1} ]\f$
!       REAL,DIMENSION(NLAT,NMOS)      :: QACROT  !<
!       REAL,DIMENSION(ILG)            :: RCANGAT !<Intercepted liquid water stored on canopy \f$[kg m^{-2} ]\f$
!       REAL,DIMENSION(NLAT,NMOS)      :: RCANROT !<
!       REAL,DIMENSION(ILG)            :: RHOSGAT !<Density of snow \f$[kg m^{-3} ]\f$
!       REAL,DIMENSION(NLAT,NMOS)      :: RHOSROT !<
!       REAL,DIMENSION(ILG)            :: SCANGAT !<Intercepted frozen water stored on canopy \f$[kg m^{-2} ]\f$
!       REAL,DIMENSION(NLAT,NMOS)      :: SCANROT !<
!       REAL,DIMENSION(ILG)            :: SNOGAT  !<Mass of snow pack [kg m^{-2} ]\f$
!       REAL,DIMENSION(NLAT,NMOS)      :: SNOROT  !<
!       REAL,DIMENSION(ILG)            :: TACGAT  !<Temperature of air within vegetation canopy [K]
!       REAL,DIMENSION(NLAT,NMOS)      :: TACROT  !<
!       REAL,DIMENSION(ILG,IGND)       :: TBARGAT !<Temperature of soil layers [K]
!       REAL,DIMENSION(NLAT,NMOS,IGND) :: TBARROT !<
!       REAL,DIMENSION(ILG)            :: TBASGAT !<Temperature of bedrock in third soil layer [K]
!       REAL,DIMENSION(NLAT,NMOS)      :: TBASROT !<
!       REAL,DIMENSION(ILG)            :: TCANGAT !<Vegetation canopy temperature [K]
!       REAL,DIMENSION(NLAT,NMOS)      :: TCANROT !<
!       REAL,DIMENSION(ILG,IGND)       :: THICGAT !<Volumetric frozen water content of soil layers \f$[m^3 m^{-3} ]\f$
!       REAL,DIMENSION(NLAT,NMOS,IGND) :: THICROT !<
!       REAL,DIMENSION(ILG,IGND)       :: THLQGAT !<Volumetric liquid water content of soil layers \f$[m^3 m^{-3} ]\f$
!       REAL,DIMENSION(NLAT,NMOS,IGND) :: THLQROT !<
!       REAL,DIMENSION(ILG)            :: TPNDGAT !<Temperature of ponded water [K]
!       REAL,DIMENSION(NLAT,NMOS)      :: TPNDROT !<
!       REAL                       TSFSGAT(ILG,4) !<Ground surface temperature over subarea [K]
!       REAL                 TSFSROT(NLAT,NMOS,4) !<
!       REAL,DIMENSION(ILG)            :: TSNOGAT !<Snowpack temperature [K]
!       REAL,DIMENSION(NLAT,NMOS)      :: TSNOROT !<
!       REAL,DIMENSION(ILG)            :: WSNOGAT !<Liquid water content of snow pack \f$[kg m^{-2} ]\f$
!       REAL,DIMENSION(NLAT,NMOS)      :: WSNOROT !<
!       REAL,DIMENSION(ILG)            :: ZPNDGAT !<Depth of ponded water on surface [m]
!       REAL,DIMENSION(NLAT,NMOS)      :: ZPNDROT !<
!
! !
! !     * LAND SURFACE PROGNOSTIC VARIABLES.
! !
! !
!       REAL,DIMENSION(NLAT,NMOS) :: REFROT  !<
!       REAL,DIMENSION(NLAT,NMOS) :: BCSNROT !<
!       REAL,DIMENSION(ILG)       :: REFGAT  !<
!       REAL,DIMENSION(ILG)       :: BCSNGAT !<
! !
! !     * GATHER-SCATTER INDEX ARRAYS.
! !
! !       INTEGER ILMOS (ILG) !<Index of grid cell corresponding to current element of gathered vector of land surface variables [ ]
! !       INTEGER JLMOS (ILG) !<Index of mosaic tile corresponding to current element of gathered vector of land surface variables [ ]
! !       INTEGER IWMOS (ILG) !<Index of grid cell corresponding to current element of gathered vector of inland water body variables [ ]
! !       INTEGER JWMOS (ILG) !<Index of mosaic tile corresponding to current element of gathered vector of inland water body variables [ ]
! !
! !     * CANOPY AND SOIL INFORMATION ARRAYS.
! !     * (THE LAST DIMENSION OF MOST OF THESE ARRAYS IS GIVEN BY
! !     * THE NUMBER OF SOIL LAYERS (IGND), THE NUMBER OF BROAD
! !     * VEGETATION CATEGORIES (ICAN), OR ICAN+1.
! !
!
!       REAL,DIMENSION(ILG,ICAN)       :: ACIDGAT !<Optional user-specified value of canopy near-infrared albedo to override CLASS-calculated value [ ]
!       REAL,DIMENSION(NLAT,NMOS,ICAN) :: ACIDROT !<
!       REAL,DIMENSION(ILG,ICAN)       :: ACVDGAT !<Optional user-specified value of canopy visible albedo to override CLASS-calculated value [ ]
!       REAL,DIMENSION(NLAT,NMOS,ICAN) :: ACVDROT !<
!       REAL,DIMENSION(ILG)            :: AGIDGAT !<Optional user-specified value of ground near-infrared albedo to override CLASS-calculated value [ ]
!       REAL,DIMENSION(NLAT,NMOS)      :: AGIDROT !<
!       REAL,DIMENSION(ILG)            :: AGVDGAT !<Optional user-specified value of ground visible albedo to override CLASS-calculated value [ ]
!       REAL,DIMENSION(NLAT,NMOS)      :: AGVDROT !<
!       REAL,DIMENSION(ILG)            :: ALGDGAT !<Reference albedo for dry soil [ ]
!       REAL,DIMENSION(NLAT,NMOS)      :: ALGDROT !<
!       REAL,DIMENSION(ILG)            :: ALGWGAT !<Reference albedo for saturated soil [ ]
!       REAL,DIMENSION(NLAT,NMOS)      :: ALGWROT !<
!       REAL,DIMENSION(ILG,ICP1)       :: ALICGAT !<Background average near-infrared albedo of vegetation category [ ]
!       REAL,DIMENSION(NLAT,NMOS,ICP1) :: ALICROT !<
!       REAL,DIMENSION(ILG,ICP1)       :: ALVCGAT !<Background average visible albedo of vegetation category [ ]
!       REAL,DIMENSION(NLAT,NMOS,ICP1) :: ALVCROT !<
!       REAL,DIMENSION(ILG)            :: ASIDGAT !<Optional user-specified value of snow near-infrared albedo to override CLASS-calculated value [ ]
!       REAL,DIMENSION(NLAT,NMOS)      :: ASIDROT !<
!       REAL,DIMENSION(ILG)            :: ASVDGAT !<Optional user-specified value of snow visible albedo to override CLASS-calculated value [ ]
!       REAL,DIMENSION(NLAT,NMOS)      :: ASVDROT !<
!       REAL,DIMENSION(ILG,IGND)       :: BIGAT   !<Clapp and Hornberger empirical “b” parameter [ ]
!       REAL,DIMENSION(NLAT,NMOS,IGND) :: BIROT   !<
!       REAL CLAYROT(NLAT,NMOS,IGND)              !<Percentage clay content of soil
!       REAL,DIMENSION(ILG,ICAN)       :: CMASGAT !<Maximum canopy mass for vegetation category \f$[kg m^{-2} ]\f$
!       REAL,DIMENSION(NLAT,NMOS,ICAN) :: CMASROT !<
!       REAL,DIMENSION(ILG,IGND)       :: DLZWGAT !<Permeable thickness of soil layer [m]
!       REAL,DIMENSION(NLAT,NMOS,IGND) :: DLZWROT !<
!       REAL,DIMENSION(ILG)            :: DRNGAT  !<Drainage index at bottom of soil profile [ ]
!       REAL,DIMENSION(NLAT,NMOS)      :: DRNROT  !<
!       REAL FAREROT(NLAT,NMOS)                   !<Fractional coverage of mosaic tile on modelled area
!       REAL,DIMENSION(ILG,ICP1)       :: FCANGAT !<Maximum fractional coverage of modelled area by vegetation category [ ]
!       REAL,DIMENSION(NLAT,NMOS,ICP1) :: FCANROT !<
!       REAL,DIMENSION(ILG)            :: GRKFGAT !<WATROF parameter used when running MESH code [ ]
!       REAL,DIMENSION(NLAT,NMOS)      :: GRKFROT !<
!       REAL,DIMENSION(ILG,IGND)       :: GRKSGAT !<Saturated hydraulic conductivity of soil layers \f$[m s^{-1} ]\f$
!       REAL,DIMENSION(NLAT,NMOS,IGND) :: GRKSROT !<
!       REAL,DIMENSION(ILG,IGND)       :: HCPSGAT !<Volumetric heat capacity of soil particles \f$[J m^{-3} ]\f$
!       REAL,DIMENSION(NLAT,NMOS,IGND) :: HCPSROT !<
!       REAL,DIMENSION(ILG,ICAN)       :: HGTDGAT !<Optional user-specified values of height of vegetation categories to override CLASS-calculated values [m]
!       REAL,DIMENSION(NLAT,NMOS,ICAN) :: HGTDROT !<
!       INTEGER IGDRGAT(ILG)                      !<Index of soil layer in which bedrock is encountered
!       INTEGER IGDRROT(NLAT,NMOS)                !<
!       INTEGER ISNDGAT(ILG,IGND)                 !<Integer identifier associated with sand content
!       INTEGER ISNDROT(NLAT,NMOS,IGND)           !<
!       REAL,DIMENSION(ILG,ICP1)       :: LNZ0GAT !<Natural logarithm of maximum roughness length of vegetation category [ ]
!       REAL,DIMENSION(NLAT,NMOS,ICP1) :: LNZ0ROT !<
!       INTEGER MIDROT (NLAT,NMOS)                !<Mosaic tile type identifier (1 for land surface, 0 for inland lake)
!       REAL ORGMROT(NLAT,NMOS,IGND)              !<Percentage organic matter content of soil
!       REAL,DIMENSION(ILG,ICAN)       :: PAIDGAT !<Optional user-specified value of plant area indices of vegetation categories to override CLASS-calculated values [ ]
!       REAL,DIMENSION(NLAT,NMOS,ICAN) :: PAIDROT !<
!       REAL,DIMENSION(ILG,ICAN)       :: PAMNGAT !<Minimum plant area index of vegetation category [ ]
!       REAL,DIMENSION(NLAT,NMOS,ICAN) :: PAMNROT !<
!       REAL,DIMENSION(ILG,ICAN)       :: PAMXGAT !<Minimum plant area index of vegetation category [ ]
!       REAL,DIMENSION(NLAT,NMOS,ICAN) :: PAMXROT !<
!       REAL,DIMENSION(ILG,ICAN)       :: PSGAGAT !<Soil moisture suction coefficient for vegetation category (used in stomatal resistance calculation) [ ]
!       REAL,DIMENSION(NLAT,NMOS,ICAN) :: PSGAROT !<
!       REAL,DIMENSION(ILG,ICAN)       :: PSGBGAT !<Soil moisture suction coefficient for vegetation category (used in stomatal resistance calculation) [ ]
!       REAL,DIMENSION(NLAT,NMOS,ICAN) :: PSGBROT !<
!       REAL,DIMENSION(ILG,IGND)       :: PSISGAT !<Soil moisture suction at saturation [m]
!       REAL,DIMENSION(NLAT,NMOS,IGND) :: PSISROT !<
!       REAL,DIMENSION(ILG,IGND)       :: PSIWGAT !<Soil moisture suction at wilting point [m]
!       REAL,DIMENSION(NLAT,NMOS,IGND) :: PSIWROT !<
!       REAL,DIMENSION(ILG,ICAN)       :: QA50GAT !<Reference value of incoming shortwave radiation for vegetation category (used in stomatal resistance calculation) \f$[W m^{-2} ]\f$
!       REAL,DIMENSION(NLAT,NMOS,ICAN) :: QA50ROT !<
!       REAL,DIMENSION(ILG,ICAN)       :: ROOTGAT !<Maximum rooting depth of vegetation category [m]
!       REAL,DIMENSION(NLAT,NMOS,ICAN) :: ROOTROT !<
!       REAL,DIMENSION(ILG,ICAN)       :: RSMNGAT !<Minimum stomatal resistance of vegetation category \f$[s m^{-1} ]\f$
!       REAL,DIMENSION(NLAT,NMOS,ICAN) :: RSMNROT !<
!       REAL SANDROT(NLAT,NMOS,IGND)              !<Percentage sand content of soil
!       REAL SDEPROT(NLAT,NMOS)                   !<Depth to bedrock in the soil profile
!       REAL,DIMENSION(ILG,IGND)       :: TCSGAT  !<Thermal conductivity of soil particles \f$[W m^{-1} K^{-1} ]\f$
!       REAL,DIMENSION(NLAT,NMOS,IGND) :: TCSROT  !<
!       REAL,DIMENSION(ILG,IGND)       :: THFCGAT !<Field capacity \f$[m^3 m^{-3} ]\f$
!       REAL,DIMENSION(NLAT,NMOS,IGND) :: THFCROT !<
!       REAL,DIMENSION(ILG,IGND)       :: THMGAT  !<Residual soil liquid water content remaining after freezing or evaporation \f$[m^3 m^{-3} ]\f$
!       REAL,DIMENSION(NLAT,NMOS,IGND) :: THMROT  !<
!       REAL,DIMENSION(ILG,IGND)       :: THPGAT  !<Pore volume in soil layer \f$[m^3 m^{-3} ]\f$
!       REAL,DIMENSION(NLAT,NMOS,IGND) :: THPROT  !<
!       REAL,DIMENSION(ILG,IGND)       :: THRGAT  !<Liquid water retention capacity for organic soil \f$[m^3 m^{-3} ]\f$
!       REAL,DIMENSION(NLAT,NMOS,IGND) :: THRROT  !<
!       REAL,DIMENSION(ILG,IGND)       :: THRAGAT !<Fractional saturation of soil behind the wetting front [ ]
!       REAL,DIMENSION(NLAT,NMOS,IGND) :: THRAROT !<
!       REAL,DIMENSION(ILG,ICAN)       :: VPDAGAT !<Vapour pressure deficit coefficient for vegetation category (used in stomatal resistance calculation) [ ]
!       REAL,DIMENSION(NLAT,NMOS,ICAN) :: VPDAROT !<
!       REAL,DIMENSION(ILG,ICAN)       :: VPDBGAT !<Vapour pressure deficit coefficient for vegetation category (used in stomatal resistance calculation) [ ]
!       REAL,DIMENSION(NLAT,NMOS,ICAN) :: VPDBROT !<
!       REAL,DIMENSION(ILG)            :: WFCIGAT !<WATROF parameter used when running MESH code [ ]
!       REAL,DIMENSION(NLAT,NMOS)      :: WFCIROT !<
!       REAL,DIMENSION(ILG)            :: WFSFGAT !<WATROF parameter used when running MESH code [ ]
!       REAL,DIMENSION(NLAT,NMOS)      :: WFSFROT !<
!       REAL,DIMENSION(ILG)            :: XSLPGAT !<Surface slope (used when running MESH code) [degrees]
!       REAL,DIMENSION(NLAT,NMOS)      :: XSLPROT !<
!       REAL,DIMENSION(ILG,IGND)       :: ZBTWGAT !<Depth to permeable bottom of soil layer [m]
!       REAL,DIMENSION(NLAT,NMOS,IGND) :: ZBTWROT !<
!       REAL,DIMENSION(ILG)            :: ZPLGGAT !<Maximum water ponding depth for snow-free subareas (user-specified when running MESH code) [m]
!       REAL,DIMENSION(NLAT,NMOS)      :: ZPLGROT !<
!       REAL,DIMENSION(ILG)            :: ZPLSGAT !<Maximum water ponding depth for snow-covered subareas (user-specified when running MESH code) [m]
!       REAL,DIMENSION(NLAT,NMOS)      :: ZPLSROT !<
!       REAL,DIMENSION(ILG)            :: ZSNLGAT !<Limiting snow depth below which coverage is < 100% [m]
!       REAL,DIMENSION(NLAT,NMOS)      :: ZSNLROT !<
!
!
!
!       REAL,DIMENSION(NLAT,NMOS,IGND) :: THLWROT  !<
!       REAL,DIMENSION(NLAT,NMOS)      :: ZSNOROT  !<
!       REAL,DIMENSION(NLAT,NMOS)      :: ALGWVROT !<
!       REAL,DIMENSION(NLAT,NMOS)      :: ALGWNROT !<
!       REAL,DIMENSION(NLAT,NMOS)      :: ALGDVROT !<
!       REAL,DIMENSION(NLAT,NMOS)      :: ALGDNROT !<
!       REAL,DIMENSION(NLAT,NMOS)      :: EMISROT  !<
!       REAL,DIMENSION(NLAT,NMOS,NBS)  :: SALBROT  !<
!       REAL,DIMENSION(NLAT,NMOS,NBS)  :: CSALROT  !<
!       REAL,DIMENSION(NLAT,NBS)       :: FSDBROL  !<
!       REAL,DIMENSION(NLAT,NBS)       :: FSFBROL  !<
!       REAL,DIMENSION(NLAT,NBS)       :: FSSBROL  !<
!
!       REAL,DIMENSION(ILG,IGND)       :: THLWGAT  !<
!       REAL,DIMENSION(ILG)            :: ALGWVGAT !<
!       REAL,DIMENSION(ILG)            :: ALGWNGAT !<
!       REAL,DIMENSION(ILG)            :: ALGDVGAT !<
!       REAL,DIMENSION(ILG)            :: ALGDNGAT !<
!       REAL,DIMENSION(ILG)            :: EMISGAT  !<
!       REAL SOCIROT(NLAT,NMOS)                    !<
! !
!       REAL,DIMENSION(ILG,NBS) :: FSDBGAT !<
!       REAL,DIMENSION(ILG,NBS) :: FSFBGAT !<
!       REAL,DIMENSION(ILG,NBS) :: FSSBGAT !<
!       REAL,DIMENSION(ILG,NBS) :: SALBGAT !<
!       REAL,DIMENSION(ILG,NBS) :: CSALGAT !<
!
!     * ARRAYS ASSOCIATED WITH COMMON BLOCKS.
!FLAG! >>> Not in the new structure
      REAL THPORG (  3) !<
      REAL THRORG (  3) !<
      REAL THMORG (  3) !<
      REAL BORG   (  3) !<
      REAL PSISORG(  3) !<
      REAL GRKSORG(  3) !<

      REAL CANEXT(ICAN) !<
      REAL XLEAF (ICAN) !<
      REAL ZORAT (ICAN) !<
!FLAG! <<< Not in the new structure
!
!       REAL DELZ  (IGND) !<
!       REAL ZBOT  (IGND) !<
!FLAG! >>> Not in the new structure
      REAL GROWYR (  18,4,2) !< !
!FLAG! <<< Not in the new structure
!
! !     * ATMOSPHERIC AND GRID-CONSTANT INPUT VARIABLES.
! !
!       REAL,DIMENSION(ILG)  :: CSZGAT  !<Cosine of solar zenith angle [ ]
!       REAL,DIMENSION(NLAT) :: CSZROW  !<
!       REAL,DIMENSION(ILG)  :: DLONGAT !<Longitude of grid cell (east of Greenwich) [degrees]
!       REAL,DIMENSION(NLAT) :: DLONROW !<
!       REAL,DIMENSION(NLAT) :: DLATROW !< Latitude of grid cell [degrees]
!       REAL,DIMENSION(ILG)  :: DLATGAT !<
!
!       REAL,DIMENSION(ILG)  :: FCLOGAT !<Fractional cloud cover [ ]
!       REAL,DIMENSION(NLAT) :: FCLOROW !<
!       REAL,DIMENSION(ILG)  :: FDLGAT  !<Downwelling longwave radiation at bottom of atmosphere (i.e. incident on modelled land surface elements \f$[W m^{-2} ]\f$
!       REAL,DIMENSION(NLAT) :: FDLROW  !<
!       REAL,DIMENSION(ILG)  :: FSIHGAT !<Near-infrared radiation incident on horizontal surface \f$[W m^{-2} ]\f$
!       REAL,DIMENSION(NLAT) :: FSIHROW !<
!       REAL,DIMENSION(ILG)  :: FSVHGAT !<Visible radiation incident on horizontal surface \f$[W m^{-2} ]\f$
!       REAL,DIMENSION(NLAT) :: FSVHROW !<
!       REAL,DIMENSION(NLAT) :: GCROW   !<Type identifier for grid cell (1 = sea ice, 0 = ocean, -1 = land)
!       REAL,DIMENSION(ILG)  :: GGEOGAT !<Geothermal heat flux at bottom of soil profile \f$[W m^{-2} ]\f$
!       REAL,DIMENSION(NLAT) :: GGEOROW !<
!       REAL,DIMENSION(ILG)  :: PADRGAT !<Partial pressure of dry air [Pa]
!       REAL,DIMENSION(NLAT) :: PADRROW !<
!       REAL,DIMENSION(ILG)  :: PREGAT  !<Surface precipitation rate \f$[kg m^{-2} s^{-1} ]\f$
!       REAL,DIMENSION(NLAT) :: PREROW  !<
!       REAL,DIMENSION(ILG)  :: PRESGAT !<Surface air pressure [Pa]
!       REAL,DIMENSION(NLAT) :: PRESROW !<
!       REAL,DIMENSION(ILG)  :: QAGAT   !<Specific humidity at reference height \f$[kg kg^{-1} ]\f$
!       REAL,DIMENSION(NLAT) :: QAROW   !<
!       REAL,DIMENSION(ILG)  :: RADJGAT !<Latitude of grid cell (positive north of equator) [rad]
!       REAL,DIMENSION(NLAT) :: RADJROW !<
!       REAL,DIMENSION(ILG)  :: RHOAGAT !<Density of air \f$[kg m^{-3} ]\f$
!       REAL,DIMENSION(NLAT) :: RHOAROW !<
!       REAL,DIMENSION(ILG)  :: RHSIGAT !<Density of fresh snow \f$[kg m^{-3} ]\f$
!       REAL,DIMENSION(NLAT) :: RHSIROW !<
!       REAL,DIMENSION(ILG)  :: RPCPGAT !<Rainfall rate over modelled area \f$[m s^{-1} ]\f$
!       REAL,DIMENSION(NLAT) :: RPCPROW !<
!       REAL,DIMENSION(NLAT) :: RPREROW !<Rainfall rate over modelled area \f$[kg m^{-2} s^{-1} ]\f$
!       REAL,DIMENSION(ILG)  :: SPCPGAT !<Snowfall rate over modelled area \f$[m s^{-1} ]\f$
!       REAL,DIMENSION(NLAT) :: SPCPROW !<
!       REAL,DIMENSION(NLAT) :: SPREROW !<Snowfall rate over modelled area \f$[kg m^{-2} s^{-1} ]\f$
!       REAL,DIMENSION(ILG)  :: TAGAT   !<Air temperature at reference height [K]
!       REAL,DIMENSION(NLAT) :: TAROW   !<
!       REAL,DIMENSION(ILG)  :: TADPGAT !<Dew point temperature of air [K]
!       REAL,DIMENSION(NLAT) :: TADPROW !<
!       REAL,DIMENSION(ILG)  :: TRPCGAT !<Rainfall temperature [K]
!       REAL,DIMENSION(NLAT) :: TRPCROW !<
!       REAL,DIMENSION(ILG)  :: TSPCGAT !<Snowfall temperature [K]
!       REAL,DIMENSION(NLAT) :: TSPCROW !<
!       REAL,DIMENSION(ILG)  :: ULGAT   !<Zonal component of wind velocity \f$[m s^{-1} ]\f$
!       REAL,DIMENSION(NLAT) :: ULROW   !<
!       REAL,DIMENSION(ILG)  :: VLGAT   !<Meridional component of wind velocity \f$[m s^{-1} ]\f$
!       REAL,DIMENSION(NLAT) :: VLROW   !<
!       REAL,DIMENSION(ILG)  :: VMODGAT !<Wind speed at reference height \f$[m s^{-1} ]\f$
!       REAL,DIMENSION(NLAT) :: VMODROW !<
!       REAL,DIMENSION(ILG)  :: VPDGAT  !<Vapour pressure deficit [mb]
!       REAL,DIMENSION(NLAT) :: VPDROW  !<
!       REAL,DIMENSION(ILG)  :: Z0ORGAT !<Orographic roughness length [m]
!       REAL,DIMENSION(ILG)  :: ZBLDGAT !<Atmospheric blending height for surface roughness length averaging [m]
!       REAL,DIMENSION(NLAT) :: ZBLDROW !<
!       REAL,DIMENSION(ILG)  :: ZDHGAT  !<User-specified height associated with diagnosed screen-level variables [m]
!       REAL,DIMENSION(NLAT) :: ZDHROW  !<
!       REAL,DIMENSION(ILG)  :: ZDMGAT  !<User-specified height associated with diagnosed anemometer-level wind speed [m]
!       REAL,DIMENSION(NLAT) :: ZDMROW  !<
!       REAL,DIMENSION(ILG)  :: ZRFHGAT !<Reference height associated with forcing air temperature and humidity [m]
!       REAL,DIMENSION(NLAT) :: ZRFHROW !<
!       REAL,DIMENSION(ILG)  :: ZRFMGAT !<Reference height associated with forcing wind speed [m]
!       REAL,DIMENSION(NLAT) :: ZRFMROW !<
!
!
!
!       REAL,DIMENSION(NLAT) :: UVROW   !<
!       REAL,DIMENSION(NLAT) :: XDIFFUS !<
!       REAL,DIMENSION(NLAT) :: Z0ORROW !<
!
!
!       REAL,DIMENSION(NLAT) :: FSSROW  !< Shortwave radiation \f$[W m^{-2} ]\f$
!       REAL,DIMENSION(NLAT) :: PRENROW !<
!       REAL,DIMENSION(NLAT) :: CLDTROW !<
!       REAL,DIMENSION(NLAT) :: FSGROL  !<
!       REAL,DIMENSION(NLAT) :: FLGROL  !<
!       REAL,DIMENSION(NLAT) :: GUSTROL !<
!       REAL,DIMENSION(NLAT) :: DEPBROW !<
! !
!       REAL,DIMENSION(ILG)  :: FSGGAT  !<
!       REAL,DIMENSION(ILG)  :: FLGGAT  !<
!       REAL,DIMENSION(ILG)  :: GUSTGAT !<
!       REAL,DIMENSION(ILG)  :: DEPBGAT !<
!       REAL,DIMENSION(ILG)  :: GTBS    !<
!       REAL,DIMENSION(ILG)  :: SFCUBS  !<
!       REAL,DIMENSION(ILG)  :: SFCVBS  !<
!       REAL,DIMENSION(ILG)  :: USTARBS !<
!       REAL,DIMENSION(ILG)  :: TCSNOW  !<
!       REAL,DIMENSION(ILG)  :: GSNOW   !<
!
! !
! !     * LAND SURFACE DIAGNOSTIC VARIABLES.
! !
!
!       REAL,DIMENSION(ILG)       :: ALIRGAT !<Diagnosed total near-infrared albedo of land surface [ ]
!       REAL,DIMENSION(NLAT,NMOS) :: ALIRROT !<
!       REAL,DIMENSION(NLAT)      :: ALIRROW !<
!       REAL,DIMENSION(ILG)       :: ALVSGAT !<Diagnosed total visible albedo of land surface [ ]
!       REAL,DIMENSION(NLAT,NMOS) :: ALVSROT !<
!       REAL,DIMENSION(NLAT)      :: ALVSROW !<
!       REAL,DIMENSION(ILG)       :: CDHGAT  !<Surface drag coefficient for heat [ ]
!       REAL,DIMENSION(NLAT,NMOS) :: CDHROT  !<
!       REAL,DIMENSION(NLAT)      :: CDHROW  !<
!       REAL,DIMENSION(ILG)       :: CDMGAT  !<Surface drag coefficient for momentum [ ]
!       REAL,DIMENSION(NLAT,NMOS) :: CDMROT  !<
!       REAL,DIMENSION(NLAT)      :: CDMROW  !<
!       REAL,DIMENSION(ILG)       :: DRGAT   !<Surface drag coefficient under neutral stability [ ]
!       REAL,DIMENSION(NLAT,NMOS) :: DRROT   !<
!       REAL,DIMENSION(NLAT)      :: DRROW   !<
!       REAL,DIMENSION(ILG)       :: EFGAT   !<Evaporation efficiency at ground surface [ ]
!       REAL,DIMENSION(NLAT,NMOS) :: EFROT   !<
!       REAL,DIMENSION(NLAT)      :: EFROW   !<
!       REAL,DIMENSION(ILG)       :: FLGGGAT !<Diagnosed net longwave radiation at soil surface \f$[W m^{-2} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: FLGGROT !<
!       REAL,DIMENSION(NLAT)      :: FLGGROW !<
!       REAL,DIMENSION(ILG)       :: FLGSGAT !<Diagnosed net longwave radiation at snow surface \f$[W m^{-2} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: FLGSROT !<
!       REAL,DIMENSION(NLAT)      :: FLGSROW !<
!       REAL,DIMENSION(ILG)       :: FLGVGAT !<Diagnosed net longwave radiation on vegetation canopy \f$[W m^{-2} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: FLGVROT !<
!       REAL,DIMENSION(NLAT)      :: FLGVROW !<
!       REAL,DIMENSION(ILG)       :: FSGGGAT !<Diagnosed net shortwave radiation at soil surface \f$[W m^{-2} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: FSGGROT !<
!       REAL,DIMENSION(NLAT)      :: FSGGROW !<
!       REAL,DIMENSION(ILG)       :: FSGSGAT !<Diagnosed net shortwave radiation at snow surface \f$[W m^{-2} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: FSGSROT !<
!       REAL,DIMENSION(NLAT)      :: FSGSROW !<
!       REAL,DIMENSION(ILG)       :: FSGVGAT !<Diagnosed net shortwave radiation on vegetation canopy \f$[W m^{-2} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: FSGVROT !<
!       REAL,DIMENSION(NLAT)      :: FSGVROW !<
!       REAL,DIMENSION(ILG)       :: FSNOGAT !<Diagnosed fractional snow coverage [ ]
!       REAL,DIMENSION(NLAT,NMOS) :: FSNOROT !<
!       REAL,DIMENSION(NLAT)      :: FSNOROW !<
!       REAL,DIMENSION(ILG)       :: GAGAT   !<Diagnosed product of drag coefficient and wind speed over modelled area \f$[m s^{-1} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: GAROT   !<
!       REAL,DIMENSION(NLAT)      :: GAROW   !<
!       REAL GFLXGAT(ILG,IGND)               !<Heat conduction between soil layers \f$[W m^{-2} ]\f$
!       REAL GFLXROT(NLAT,NMOS,IGND)         !<
!       REAL GFLXROW(NLAT,IGND)              !<
!       REAL,DIMENSION(ILG)       :: GTGAT   !<Diagnosed effective surface black-body temperature [K]
!       REAL,DIMENSION(NLAT,NMOS) :: GTROT   !<
!       REAL,DIMENSION(NLAT)      :: GTROW   !<
!       REAL,DIMENSION(ILG)       :: HBLGAT  !<Height of the atmospheric boundary layer [m]
!       REAL,DIMENSION(NLAT,NMOS) :: HBLROT  !<
!       REAL,DIMENSION(NLAT)      :: HBLROW  !<
!       REAL,DIMENSION(ILG)       :: HEVCGAT !<Diagnosed latent heat flux on vegetation canopy \f$[W m^{-2} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: HEVCROT !<
!       REAL,DIMENSION(NLAT)      :: HEVCROW !<
!       REAL,DIMENSION(ILG)       :: HEVGGAT !<Diagnosed latent heat flux at soil surface \f$[W m^{-2} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: HEVGROT !<
!       REAL,DIMENSION(NLAT)      :: HEVGROW !<
!       REAL,DIMENSION(ILG)       :: HEVSGAT !<Diagnosed latent heat flux at snow surface \f$[W m^{-2} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: HEVSROT !<
!       REAL,DIMENSION(NLAT)      :: HEVSROW !<
!       REAL,DIMENSION(ILG)       :: HFSGAT  !<Diagnosed total surface sensible heat flux over modelled area \f$[W m^{-2} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: HFSROT  !<
!       REAL,DIMENSION(NLAT)      :: HFSROW  !<
!       REAL,DIMENSION(ILG)       :: HFSCGAT !<Diagnosed sensible heat flux on vegetation canopy \f$[W m^{-2} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: HFSCROT !<
!       REAL,DIMENSION(NLAT)      :: HFSCROW !<
!       REAL,DIMENSION(ILG)       :: HFSGGAT !<Diagnosed sensible heat flux at soil surface \f$[W m^{-2} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: HFSGROT !<
!       REAL,DIMENSION(NLAT)      :: HFSGROW !<
!       REAL,DIMENSION(ILG)       :: HFSSGAT !<Diagnosed sensible heat flux at snow surface \f$[W m^{-2} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: HFSSROT !<
!       REAL,DIMENSION(NLAT)      :: HFSSROW !<
!       REAL,DIMENSION(ILG)       :: HMFCGAT !<Diagnosed energy associated with phase change of water on vegetation \f$[W m^{-2} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: HMFCROT !<
!       REAL,DIMENSION(NLAT)      :: HMFCROW !<
!       REAL HMFGGAT(ILG,IGND)               !<Diagnosed energy associated with phase change of water in soil layers \f$[W m^{-2} ]\f$
!       REAL HMFGROT(NLAT,NMOS,IGND)         !<
!       REAL HMFGROW(NLAT,IGND)              !<
!       REAL,DIMENSION(ILG)       :: HMFNGAT !<Diagnosed energy associated with phase change of water in snow pack \f$[W m^{-2} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: HMFNROT !<
!       REAL,DIMENSION(NLAT)      :: HMFNROW !<
!       REAL HTCGAT (ILG,IGND)               !<Diagnosed internal energy change of soil layer due to conduction and/or change in mass \f$[W m^{-2} ]\f$
!       REAL HTCROT (NLAT,NMOS,IGND)         !<
!       REAL HTCROW (NLAT,IGND)              !<
!       REAL,DIMENSION(ILG)       :: HTCCGAT !<Diagnosed internal energy change of vegetation canopy due to conduction and/or change in mass \f$[W m^{-2} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: HTCCROT !<
!       REAL,DIMENSION(NLAT)      :: HTCCROW !<
!       REAL,DIMENSION(ILG)       :: HTCSGAT !<Diagnosed internal energy change of snow pack due to conduction and/or change in mass \f$[W m^{-2} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: HTCSROT !<
!       REAL,DIMENSION(NLAT)      :: HTCSROW !<
!       REAL,DIMENSION(ILG)       :: ILMOGAT !<Inverse of Monin-Obukhov roughness length \f$(m^{-1} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: ILMOROT !<
!       REAL,DIMENSION(NLAT)      :: ILMOROW !<
!FLAG! >>> Not in the new structure
      INTEGER ISUM(6)                      !<Total number of iterations required to solve surface energy balance for the elements of the four subareas for the current run
!FLAG! <<< Not in the new structure
!       INTEGER ITCTGAT(ILG,6,50)            !<Counter of number of iterations required to solve surface energy balance for the elements of the four subareas
!       INTEGER ITCTROT(NLAT,NMOS,6,50)      !<
!       REAL,DIMENSION(ILG)       :: PCFCGAT !<Diagnosed frozen precipitation intercepted by vegetation \f$[kg m^{-2} s^{-1} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: PCFCROT !<
!       REAL,DIMENSION(NLAT)      :: PCFCROW !<
!       REAL,DIMENSION(ILG)       :: PCLCGAT !<Diagnosed liquid precipitation intercepted by vegetation \f$[kg m^{-2} s^{-1} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: PCLCROT !<
!       REAL,DIMENSION(NLAT)      :: PCLCROW !<
!       REAL,DIMENSION(ILG)       :: PCPGGAT !<Diagnosed precipitation incident on ground \f$[kg m^{-2} s^{-1} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: PCPGROT !<
!       REAL,DIMENSION(NLAT)      :: PCPGROW !<
!       REAL,DIMENSION(ILG)       :: PCPNGAT !<Diagnosed precipitation incident on snow pack \f$[kg m^{-2} s^{-1} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: PCPNROT !<
!       REAL,DIMENSION(NLAT)      :: PCPNROW !<
!       REAL,DIMENSION(ILG)       :: PETGAT  !<Diagnosed potential evapotranspiration \f$[kg m^{-2} s^{-1} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: PETROT  !<
!       REAL,DIMENSION(NLAT)      :: PETROW  !<
!       REAL,DIMENSION(ILG)       :: QEVPGAT !<Diagnosed total surface latent heat flux over modelled area \f$[W m^{-2} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: QEVPROT !<
!       REAL,DIMENSION(NLAT)      :: QEVPROW !<
!       REAL QFCGAT (ILG,IGND)               !<Diagnosed vapour flux from transpiration over modelled area \f$[W m^{-2} ]\f$
!       REAL QFCROT (NLAT,NMOS,IGND)         !<
!       REAL QFCROW (NLAT,IGND)              !<
!       REAL,DIMENSION(ILG)       :: QFCFGAT !<Diagnosed vapour flux from frozen water on vegetation \f$[kg m^{-2} s^{-1} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: QFCFROT !<
!       REAL,DIMENSION(NLAT)      :: QFCFROW !<
!       REAL,DIMENSION(ILG)       :: QFCLGAT !<Diagnosed vapour flux from liquid water on vegetation \f$[kg m^{-2} s^{-1} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: QFCLROT !<
!       REAL,DIMENSION(NLAT)      :: QFCLROW !<
!       REAL,DIMENSION(ILG)       :: QFGGAT  !<Diagnosed water vapour flux from ground \f$[kg m^{-2} s^{-1} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: QFGROT  !<
!       REAL,DIMENSION(NLAT)      :: QFGROW  !<
!       REAL,DIMENSION(ILG)       :: QFNGAT  !<Diagnosed water vapour flux from snow pack \f$[kg m^{-2} s^{-1} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: QFNROT  !<
!       REAL,DIMENSION(NLAT)      :: QFNROW  !<
!       REAL,DIMENSION(ILG)       :: QFSGAT  !<Diagnosed total surface water vapour flux over modelled area \f$[kg m^{-2} s^{-1} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: QFSROT  !<
!       REAL,DIMENSION(NLAT)      :: QFSROW  !<
!       REAL,DIMENSION(ILG)       :: QFXGAT  !<Product of surface drag coefficient, wind speed and surface-air specific humidity difference \f$[m s^{-1} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: QFXROT  !<
!       REAL,DIMENSION(NLAT)      :: QFXROW  !<
!       REAL,DIMENSION(ILG)       :: QGGAT   !<Diagnosed surface specific humidity \f$[kg kg^{-1} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: QGROT   !<
!       REAL,DIMENSION(NLAT)      :: QGROW   !<
!       REAL,DIMENSION(ILG)       :: ROFGAT  !<Total runoff from soil \f$[kg m^{-2} s^{-1} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: ROFROT  !<
!       REAL,DIMENSION(NLAT)      :: ROFROW  !<
!       REAL,DIMENSION(ILG)       :: ROFBGAT !<Base flow from bottom of soil column \f$[kg m^{-2} s^{-1} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: ROFBROT !<
!       REAL,DIMENSION(NLAT)      :: ROFBROW !<
!       REAL,DIMENSION(ILG)       :: ROFCGAT !<Liquid/frozen water runoff from vegetation \f$[kg m^{-2} s^{-1} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: ROFCROT !<
!       REAL,DIMENSION(NLAT)      :: ROFCROW !<
!       REAL,DIMENSION(ILG)       :: ROFNGAT !<Liquid water runoff from snow pack \f$[kg m^{-2} s^{-1} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: ROFNROT !<
!       REAL,DIMENSION(NLAT)      :: ROFNROW !<
!       REAL,DIMENSION(ILG)       :: ROFOGAT !<Overland flow from top of soil column \f$[kg m^{-2} s^{-1} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: ROFOROT !<
!       REAL,DIMENSION(NLAT)      :: ROFOROW !<
!       REAL,DIMENSION(ILG)       :: ROFSGAT !<Interflow from sides of soil column \f$[kg m^{-2} s^{-1} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: ROFSROT !<
!       REAL,DIMENSION(NLAT)      :: ROFSROW !<
!       REAL,DIMENSION(ILG)       :: ROVGGAT !<Diagnosed liquid/frozen water runoff from vegetation to ground surface \f$[kg m^{-2} s^{-1} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: ROVGROT !<
!       REAL,DIMENSION(NLAT)      :: ROVGROW !<
!       REAL,DIMENSION(ILG)       :: SFCQGAT !<Diagnosed screen-level specific humidity \f$[kg kg^{-1} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: SFCQROT !<
!       REAL,DIMENSION(NLAT)      :: SFCQROW !<
!       REAL,DIMENSION(ILG)       :: SFCTGAT !<Diagnosed screen-level air temperature [K]
!       REAL,DIMENSION(NLAT,NMOS) :: SFCTROT !<
!       REAL,DIMENSION(NLAT)      :: SFCTROW !<
!       REAL,DIMENSION(ILG)       :: SFCUGAT !<Diagnosed anemometer-level zonal wind \f$[m s^{-1} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: SFCUROT !<
!       REAL,DIMENSION(NLAT)      :: SFCUROW !<
!       REAL,DIMENSION(ILG)       :: SFCVGAT !<Diagnosed anemometer-level meridional wind \f$[m s^{-1} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: SFCVROT !<
!       REAL,DIMENSION(NLAT)      :: SFCVROW !<
!       REAL,DIMENSION(ILG)       :: TFXGAT  !<Product of surface drag coefficient, wind speed and surface-air temperature difference \f$[K m s^{-1} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: TFXROT  !<
!       REAL,DIMENSION(NLAT)      :: TFXROW  !<
!       REAL,DIMENSION(ILG)       :: TROBGAT !<Temperature of base flow from bottom of soil column [K]
!       REAL,DIMENSION(NLAT,NMOS) :: TROBROT !<
!       REAL,DIMENSION(ILG)       :: TROFGAT !<Temperature of total runoff [K]
!       REAL,DIMENSION(NLAT,NMOS) :: TROFROT !<
!       REAL,DIMENSION(ILG)       :: TROOGAT !<Temperature of overland flow from top of soil column [K]
!       REAL,DIMENSION(NLAT,NMOS) :: TROOROT !<
!       REAL,DIMENSION(ILG)       :: TROSGAT !<Temperature of interflow from sides of soil column [K]
!       REAL,DIMENSION(NLAT,NMOS) :: TROSROT !<
!       REAL,DIMENSION(ILG)       :: UEGAT   !<Friction velocity of air \f$[m s^{-1} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: UEROT   !<
!       REAL,DIMENSION(NLAT)      :: UEROW   !<
!       REAL,DIMENSION(ILG)       :: WTABGAT !<Depth of water table in soil [m]
!       REAL,DIMENSION(NLAT,NMOS) :: WTABROT !<
!       REAL,DIMENSION(NLAT)      :: WTABROW !<
!       REAL,DIMENSION(ILG)       :: WTRCGAT !<Diagnosed residual water transferred off the vegetation canopy \f$[kg m^{-2} s^{-1} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: WTRCROT !<
!       REAL,DIMENSION(NLAT)      :: WTRCROW !<
!       REAL,DIMENSION(ILG)       :: WTRGGAT !<Diagnosed residual water transferred into or out of the soil \f$[kg m^{-2} s^{-1} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: WTRGROT !<
!       REAL,DIMENSION(NLAT)      :: WTRGROW !<
!       REAL,DIMENSION(ILG)       :: WTRSGAT !<Diagnosed residual water transferred into or out of the snow pack \f$[kg m^{-2} s^{-1} ]\f$
!       REAL,DIMENSION(NLAT,NMOS) :: WTRSROT !<
!       REAL,DIMENSION(NLAT)      :: WTRSROW !<
!
!       REAL,DIMENSION(ILG)       :: QLWOGAT !<
!       REAL,DIMENSION(ILG)       :: SFRHGAT !<
!       REAL,DIMENSION(NLAT,NMOS) :: SFRHROT !<
!       REAL,DIMENSION(NLAT)      :: SFRHROW !<
!
!       REAL,DIMENSION(ILG)       :: FTEMP   !<
!       REAL,DIMENSION(ILG)       :: FVAP    !<
!       REAL,DIMENSION(ILG)       :: RIB     !<
!
!    * ARRAYS USED FOR OUTPUT AND DISPLAY PURPOSES.
!    * (THE SUFFIX "ACC" REFERS TO ACCUMULATOR ARRAYS USED IN
!    * CALCULATING TIME AVERAGES.)
!FLAG! >>> Not in the new structure
      CHARACTER     TITLE1*4,     TITLE2*4,     TITLE3*4,&
     &              TITLE4*4,     TITLE5*4,     TITLE6*4
      CHARACTER     NAME1*4,      NAME2*4,      NAME3*4,&
     &              NAME4*4,      NAME5*4,      NAME6*4
      CHARACTER     PLACE1*4,     PLACE2*4,     PLACE3*4,&
     &              PLACE4*4,     PLACE5*4,     PLACE6*4
!FLAG! <<< Not in the new structure
!       REAL,DIMENSION(NLAT) :: ALIRACC !<Diagnosed total near-infrared albedo of land surface [ ]
!       REAL,DIMENSION(NLAT) :: ALVSACC !<Diagnosed total visible albedo of land surface [ ]
!       REAL,DIMENSION(NLAT) :: EVAPACC !<Diagnosed total surface water vapour flux over modelled area \f$[kg m^{-2} ]\f$
!       REAL,DIMENSION(NLAT) :: FLINACC !<Downwelling longwave radiation above surface \f$[W m^{-2} ]\f$
!       REAL,DIMENSION(NLAT) :: FLUTACC !<Upwelling longwave radiation from surface \f$[W m^{-2} ]\f$
!       REAL,DIMENSION(NLAT) :: FSINACC !<Downwelling shortwave radiation above surface \f$[W m^{-2} ]\f$
!       REAL,DIMENSION(NLAT) :: GROACC  !<Vegetation growth index [ ]
!       REAL,DIMENSION(NLAT) :: GTACC   !<Diagnosed effective surface black-body temperature [K]
!       REAL,DIMENSION(NLAT) :: HFSACC  !<Diagnosed total surface sensible heat flux over modelled area \f$[W m^{-2} ]\f$
!       REAL,DIMENSION(NLAT) :: HMFNACC !<Diagnosed energy associated with phase change of water in snow pack \f$[W m^{-2} ]\f$
!       REAL,DIMENSION(NLAT) :: OVRACC  !<Overland flow from top of soil column \f$[kg m^{-2} ]\f$
!       REAL,DIMENSION(NLAT) :: PREACC  !<Surface precipitation rate \f$[kg m^{-2} ]\f$
!       REAL,DIMENSION(NLAT) :: PRESACC !<Surface air pressure [Pa]
!       REAL,DIMENSION(NLAT) :: QAACC   !<Specific humidity at reference height \f$[kg kg^{-1} ]\f$
!       REAL,DIMENSION(NLAT) :: QEVPACC !<Diagnosed total surface latent heat flux over modelled area \f$[W m^{-2} ]\f$
!       REAL,DIMENSION(NLAT) :: RCANACC !<Intercepted liquid water stored on canopy \f$[kg m^{-2} ]\f$
!       REAL,DIMENSION(NLAT) :: RHOSACC !<Density of snow \f$[kg m^{-3} ]\f$
!       REAL,DIMENSION(NLAT) :: ROFACC  !<Total runoff from soil \f$[kg m^{-2} ]\f$
!       REAL,DIMENSION(NLAT) :: SCANACC !<Intercepted frozen water stored on canopy \f$[kg m^{-2} ]\f$
!       REAL,DIMENSION(NLAT) :: SNOACC  !<Mass of snow pack \f$[kg m^{-2} ]\f$
!       REAL,DIMENSION(NLAT) :: TAACC   !<Air temperature at reference height [K]
!       REAL TBARACC(NLAT,IGND)         !<Temperature of soil layers [K]
!       REAL THALACC(NLAT,IGND)         !<Total volumetric water content of soil layers \f$[m^3 m^{-3} ]\f$
!       REAL THICACC(NLAT,IGND)         !<Volumetric frozen water content of soil layers \f$[m^3 m^{-3} ]\f$
!       REAL THLQACC(NLAT,IGND)         !<Volumetric liquid water content of soil layers \f$[m^3 m^{-3} ]\f$
!       REAL,DIMENSION(NLAT) :: TCANACC !<Vegetation canopy temperature [K]
!       REAL,DIMENSION(NLAT) :: TSNOACC !<Snowpack temperature [K]
!       REAL,DIMENSION(NLAT) :: UVACC   !<Wind speed \f$[m s^{-1} ]\f$
!       REAL,DIMENSION(NLAT) :: WSNOACC !<Liquid water content of snow pack \f$[kg m^{-2} ]\f$
!       REAL,DIMENSION(NLAT) :: WTBLACC !<Depth of water table in soil [m]
!       REAL,DIMENSION(NLAT) :: ALTOTACC!<Broadband albedo [-]
!
!       REAL,DIMENSION(NLAT) :: CANARE  !<
!       REAL,DIMENSION(NLAT) :: SNOARE  !<
!
! !
! !     * ARRAYS DEFINED TO PASS INFORMATION BETWEEN THE THREE MAJOR
! !     * SUBSECTIONS OF CLASS ("CLASSA", "CLASST" AND "CLASSW").
!  !<
!       REAL,DIMENSION(ILG,IGND) :: TBARC  !<
!       REAL,DIMENSION(ILG,IGND) :: TBARG  !<
!       REAL,DIMENSION(ILG,IGND) :: TBARCS !<
!       REAL,DIMENSION(ILG,IGND) :: TBARGS !<
!       REAL,DIMENSION(ILG,IGND) :: THLIQC !<
!       REAL,DIMENSION(ILG,IGND) :: THLIQG !<
!       REAL,DIMENSION(ILG,IGND) :: THICEC !<
!       REAL,DIMENSION(ILG,IGND) :: THICEG !<
!       REAL,DIMENSION(ILG,IGND) :: FROOT  !<
!       REAL,DIMENSION(ILG,IGND) :: HCPC   !<
!       REAL,DIMENSION(ILG,IGND) :: HCPG   !<
!       REAL,DIMENSION(ILG,IGND) :: FROOTS !<
!       REAL,DIMENSION(ILG,IGND) :: TCTOPC !<
!       REAL,DIMENSION(ILG,IGND) :: TCBOTC !<
!       REAL,DIMENSION(ILG,IGND) :: TCTOPG !<
!       REAL,DIMENSION(ILG,IGND) :: TCBOTG !<
! !
!       REAL FC     (ILG)  !<
!       REAL FG     (ILG)  !<
!       REAL FCS    (ILG)  !<
!       REAL FGS    (ILG)  !<
!       REAL RBCOEF (ILG)  !<
!       REAL ZSNOW  (ILG)  !<
!       REAL FSVF   (ILG)  !<
!       REAL FSVFS  (ILG)  !<
!       REAL ALVSCN (ILG)  !<
!       REAL ALIRCN (ILG)  !<
!       REAL ALVSG  (ILG)  !<
!       REAL ALIRG  (ILG)  !<
!       REAL ALVSCS (ILG)  !<
!       REAL ALIRCS (ILG)  !<
!       REAL ALVSSN (ILG)  !<
!       REAL ALIRSN (ILG)  !<
!       REAL ALVSGC (ILG)  !<
!       REAL ALIRGC (ILG)  !<
!       REAL ALVSSC (ILG)  !<
!       REAL ALIRSC (ILG)  !<
!       REAL TRVSCN (ILG)  !<
!       REAL TRIRCN (ILG)  !<
!       REAL TRVSCS (ILG)  !<
!       REAL TRIRCS (ILG)  !<
!       REAL RC     (ILG)  !<
!       REAL RCS    (ILG)  !<
!       REAL FRAINC (ILG)  !<
!       REAL FSNOWC (ILG)  !<
!       REAL FRAICS (ILG)  !<
!       REAL FSNOCS (ILG)  !<
!       REAL CMASSC (ILG)  !<
!       REAL CMASCS (ILG)  !<
!       REAL DISP   (ILG)  !<
!       REAL DISPS  (ILG)  !<
!       REAL ZOMLNC (ILG)  !<
!       REAL ZOELNC (ILG)  !<
!       REAL ZOMLNG (ILG)  !<
!       REAL ZOELNG (ILG)  !<
!       REAL ZOMLCS (ILG)  !<
!       REAL ZOELCS (ILG)  !<
!       REAL ZOMLNS (ILG)  !<
!       REAL ZOELNS (ILG)  !<
!       REAL TRSNOWC (ILG) !<
!       REAL CHCAP  (ILG)  !<
!       REAL CHCAPS (ILG)  !<
!       REAL GZEROC (ILG)  !<
!       REAL GZEROG (ILG)  !<
!       REAL GZROCS (ILG)  !<
!       REAL GZROGS (ILG)  !<
!       REAL G12C   (ILG)  !<
!       REAL G12G   (ILG)  !<
!       REAL G12CS  (ILG)  !<
!       REAL G12GS  (ILG)  !<
!       REAL G23C   (ILG)  !<
!       REAL G23G   (ILG)  !<
!       REAL G23CS  (ILG)  !<
!       REAL G23GS  (ILG)  !<
!       REAL QFREZC (ILG)  !<
!       REAL QFREZG (ILG)  !<
!       REAL QMELTC (ILG)  !<
!       REAL QMELTG (ILG)  !<
!       REAL EVAPC  (ILG)  !<
!       REAL EVAPCG (ILG)  !<
!       REAL EVAPG  (ILG)  !<
!       REAL EVAPCS (ILG)  !<
!       REAL EVPCSG (ILG)  !<
!       REAL EVAPGS (ILG)  !<
!       REAL TCANO  (ILG)  !<
!       REAL TCANS  (ILG)  !<
!       REAL RAICAN (ILG)  !<
!       REAL SNOCAN (ILG)  !<
!       REAL RAICNS (ILG)  !<
!       REAL SNOCNS (ILG)  !<
!       REAL CWLCAP (ILG)  !<
!       REAL CWFCAP (ILG)  !<
!       REAL CWLCPS (ILG)  !<
!       REAL CWFCPS (ILG)  !<
!       REAL TSNOCS (ILG)  !<
!       REAL TSNOGS (ILG)  !<
!       REAL RHOSCS (ILG)  !<
!       REAL RHOSGS (ILG)  !<
!       REAL WSNOCS (ILG)  !<
!       REAL WSNOGS (ILG)  !<
!       REAL TPONDC (ILG)  !<
!       REAL TPONDG (ILG)  !<
!       REAL TPNDCS (ILG)  !<
!       REAL TPNDGS (ILG)  !<
!       REAL ZPLMCS (ILG)  !<
!       REAL ZPLMGS (ILG)  !<
!       REAL ZPLIMC (ILG)  !<
!       REAL ZPLIMG (ILG)  !<
! !
!       REAL ALTG(ILG,NBS)    !<
!       REAL ALSNO(ILG,NBS)   !<
!       REAL TRSNOWG(ILG,NBS) !<
!
! !
! !     * DIAGNOSTIC ARRAYS USED FOR CHECKING ENERGY AND WATER
! !     * BALANCES.
! !
!       REAL CTVSTP(ILG) !<
!       REAL CTSSTP(ILG) !<
!       REAL CT1STP(ILG) !<
!       REAL CT2STP(ILG) !<
!       REAL CT3STP(ILG) !<
!       REAL WTVSTP(ILG) !<
!       REAL WTSSTP(ILG) !<
!       REAL WTGSTP(ILG) !<

!     * CONSTANTS AND TEMPORARY VARIABLES.
!
      REAL DEGLON,DAY,DECL,HOUR,COSZ,CUMSNO,EVAPSUM,&
     &     QSUMV,QSUMS,QSUM1,QSUM2,QSUM3,WSUMV,WSUMS,WSUMG,ALTOT,&
     &     FSSTAR,FLSTAR,QH,QE,BEG,SNOMLT,ZSN,TCN,TSN,TPN,GTOUT,TAC,&
     &     TSURF
!
!     * COMMON BLOCK PARAMETERS.
!
      REAL X1,X2,X3,X4,G,GAS,X5,X6,CPRES,GASV,X7,CPI,X8,CELZRO,X9,&
     &     X10,X11,X12,X13,X14,X15,SIGMA,X16,DELTIM,DELT,TFREZ,&
     &     RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN,TCW,TCICE,TCSAND,TCCLAY,&
     &     TCOM,TCDRYS,RHOSOL,RHOOM,HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,&
     &     HCPCLY,SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,TCGLAC,CLHMLT,&
     &     CLHVAP,PI,ZOLNG,ZOLNS,ZOLNI,ZORATG,ALVSI,ALIRI,ALVSO,ALIRO,&
     &     ALBRCK,DELTA,CGRAV,CKARM,CPD,AS,ASX,CI,BS,BETA,FACTN,HMIN,&
     &     ANGMAX,A,B


!
!================= CTEM array declaration ===============================\
!
!     Local variables for coupling CLASS and CTEM
!
      integer strlen
      character*80   titlec1
      character*80   argbuff
      character*160  command

       integer   lopcount,  isumc,     k1c,       k2c,&
     &           jhhstd,    jhhendd,   jdstd,   jdendd,&
     &           jhhsty,     jhhendy,   jdsty,   jdendy,&
     &           month1,     month2,      xday,  ctemloop,&
     &           nummetcylyrs, ncyear,  co2yr,   spinfast,&
     &           nol2pfts(4),  popyr, metcylyrst, metcycendyr,&
     &           climiyear,   popcycleyr,    cypopyr, lucyr,&
     &           cylucyr, endyr,bigpftc(1), obswetyr,&
     &           cywetldyr, trans_startyr, jmosty, obslghtyr,&
     &           lath, testyr

      real      co2concin,    setco2conc, sumfare,&
     &           temp_var, barefrac,ch4concin, setch4conc

      integer, allocatable, dimension(:) :: curlatno !ilg
      integer, allocatable, dimension(:) :: altotcount_ctm !nlat
      real, allocatable, dimension(:,:)  :: todfrac  !(ilg,icc)
      real, allocatable, dimension(:,:)  :: barf  !(nlat,nmos)
      real, allocatable, dimension(:)    :: currlat !(ilg)
      real, allocatable, dimension(:)    :: wl !(lat)
      real, allocatable, dimension(:)    :: grclarea !(ilg)
      real, allocatable, dimension(:)    :: wossl !(lat)
      real, allocatable, dimension(:)    :: sl !(lat)
      real, allocatable, dimension(:)    :: radl !(lat)
      real, allocatable, dimension(:)    :: cl !(lat)
      real, allocatable, dimension(:)    :: ml !(ilg)
      real, allocatable, dimension(:)    :: fsinacc_gat !(ilg)
      real, allocatable, dimension(:)    :: flutacc_gat !(ilg)
      real, allocatable, dimension(:)    :: flinacc_gat !(ilg)
      real, allocatable, dimension(:)    :: alswacc_gat !(ilg)
      real, allocatable, dimension(:)    :: allwacc_gat !(ilg)
      real, allocatable, dimension(:)    :: pregacc_gat !(ilg)
      real, allocatable, dimension(:)    :: altotacc_gat !(ilg)
      real, allocatable, dimension(:)    :: netrad_gat !(ilg)
      real, allocatable, dimension(:)    :: preacc_gat !(ilg)
      real, allocatable, dimension(:)    :: sdepgat !(ilg)
      real, allocatable, dimension(:,:)  :: rgmgat !(ilg,ignd)
      real, allocatable, dimension(:,:)  :: sandgat !(ilg,ignd)
      real, allocatable, dimension(:,:)  :: claygat !(ilg,ignd)
      real, allocatable, dimension(:,:)  :: orgmgat !(ilg,ignd)
      real, allocatable, dimension(:)    :: xdiffusgat !(ilg) ! the corresponding ROW is CLASS's XDIFFUS
      real, allocatable, dimension(:)    :: faregat !(ilg)   ! the ROT is FAREROT
      real, allocatable, dimension(:,:)  :: FTABLE !(NLAT,NMOS) !,ALAVG,ALMAX,FTAVG,FTMAX
      real, allocatable, dimension(:,:)  :: ACTLYR !(NLAT,NMOS)

       real fsstar_gat, flstar_gat !FLAG should these have more dimensions? JM Feb 2016.

      ! Model switches:
      logical, pointer :: ctem_on
      logical, pointer :: parallelrun
      logical, pointer :: cyclemet
      logical, pointer :: dofire
      logical, pointer :: run_model
      logical, pointer :: met_rewound
      logical, pointer :: reach_eof
      logical, pointer :: compete
      logical, pointer :: start_bare
      logical, pointer :: rsfile
      logical, pointer :: lnduseon
      logical, pointer :: co2on
      logical, pointer :: ch4on
      logical, pointer :: popdon
      logical, pointer :: inibioclim
      logical, pointer :: dowetlands
      logical, pointer :: obswetf
      logical, pointer :: transient_run
      logical, pointer :: use_netcdf
      character(180), pointer :: met_file
      character(180), pointer :: init_file

      ! ROW vars:
      logical, pointer, dimension(:,:,:) :: pftexistrow
      integer, pointer, dimension(:,:,:) :: colddaysrow
      integer, pointer, dimension(:,:) :: icountrow
      integer, pointer, dimension(:,:,:) :: lfstatusrow
      integer, pointer, dimension(:,:,:) :: pandaysrow
      integer, pointer, dimension(:,:) :: stdalnrow
      real, pointer, dimension(:,:) :: tcanrs
      real, pointer, dimension(:,:) :: tsnors
      real, pointer, dimension(:,:) :: tpndrs
      real, pointer, dimension(:,:,:) :: csum
      real, pointer, dimension(:,:,:) :: tbaraccrow_m
      real, pointer, dimension(:,:) :: tcanoaccrow_m
      real, pointer, dimension(:,:) :: uvaccrow_m
      real, pointer, dimension(:,:) :: vvaccrow_m

      real, pointer, dimension(:,:,:) :: ailcminrow         !
      real, pointer, dimension(:,:,:) :: ailcmaxrow         !
      real, pointer, dimension(:,:,:) :: dvdfcanrow         !
      real, pointer, dimension(:,:,:) :: gleafmasrow        !
      real, pointer, dimension(:,:,:) :: bleafmasrow        !
      real, pointer, dimension(:,:,:) :: stemmassrow        !
      real, pointer, dimension(:,:,:) :: rootmassrow        !
      real, pointer, dimension(:,:,:) :: pstemmassrow       !
      real, pointer, dimension(:,:,:) :: pgleafmassrow      !
      real, pointer, dimension(:,:,:) :: fcancmxrow
      real, pointer, dimension(:,:) :: gavglairow
      real, pointer, dimension(:,:,:) :: zolncrow
      real, pointer, dimension(:,:,:) :: ailcrow
      real, pointer, dimension(:,:,:) :: ailcgrow
      real, pointer, dimension(:,:,:) :: ailcgsrow
      real, pointer, dimension(:,:,:) :: fcancsrow
      real, pointer, dimension(:,:,:) :: fcancrow
      real, pointer, dimension(:,:) :: co2concrow
      real, pointer, dimension(:,:) :: ch4concrow
      real, pointer, dimension(:,:,:) :: co2i1cgrow
      real, pointer, dimension(:,:,:) :: co2i1csrow
      real, pointer, dimension(:,:,:) :: co2i2cgrow
      real, pointer, dimension(:,:,:) :: co2i2csrow
      real, pointer, dimension(:,:,:) :: ancsvegrow
      real, pointer, dimension(:,:,:) :: ancgvegrow
      real, pointer, dimension(:,:,:) :: rmlcsvegrow
      real, pointer, dimension(:,:,:) :: rmlcgvegrow
      real, pointer, dimension(:,:,:) :: slairow
      real, pointer, dimension(:,:,:) :: ailcbrow
      real, pointer, dimension(:,:) :: canresrow
      real, pointer, dimension(:,:,:) :: flhrlossrow

      real, pointer, dimension(:,:,:) :: grwtheffrow
      real, pointer, dimension(:,:,:) :: lystmmasrow
      real, pointer, dimension(:,:,:) :: lyrotmasrow
      real, pointer, dimension(:,:,:) :: tymaxlairow
      real, pointer, dimension(:,:) :: vgbiomasrow
      real, pointer, dimension(:,:) :: gavgltmsrow
      real, pointer, dimension(:,:) :: gavgscmsrow
      real, pointer, dimension(:,:,:) :: stmhrlosrow
      real, pointer, dimension(:,:,:,:) :: rmatcrow
      real, pointer, dimension(:,:,:,:) :: rmatctemrow
      real, pointer, dimension(:,:,:) :: litrmassrow
      real, pointer, dimension(:,:,:) :: soilcmasrow
      real, pointer, dimension(:,:,:) :: vgbiomas_vegrow

      real, pointer, dimension(:,:,:) :: emit_co2row
      real, pointer, dimension(:,:,:) :: emit_corow
      real, pointer, dimension(:,:,:) :: emit_ch4row
      real, pointer, dimension(:,:,:) :: emit_nmhcrow
      real, pointer, dimension(:,:,:) :: emit_h2row
      real, pointer, dimension(:,:,:) :: emit_noxrow
      real, pointer, dimension(:,:,:) :: emit_n2orow
      real, pointer, dimension(:,:,:) :: emit_pm25row
      real, pointer, dimension(:,:,:) :: emit_tpmrow
      real, pointer, dimension(:,:,:) :: emit_tcrow
      real, pointer, dimension(:,:,:) :: emit_ocrow
      real, pointer, dimension(:,:,:) :: emit_bcrow
      real, pointer, dimension(:,:) :: burnfracrow
      real, pointer, dimension(:,:,:) :: burnvegfrow
      real, pointer, dimension(:,:,:) :: smfuncvegrow
      real, pointer, dimension(:,:) :: popdinrow
      real, pointer, dimension(:,:,:) :: btermrow
      real, pointer, dimension(:,:) :: ltermrow
      real, pointer, dimension(:,:,:) :: mtermrow

      real, pointer, dimension(:,:) :: extnprobrow
      real, pointer, dimension(:,:) :: prbfrhucrow
      real, pointer, dimension(:,:,:) :: mlightngrow
      real, pointer, dimension(:) :: dayl_maxrow
      real, pointer, dimension(:) :: daylrow


      real, pointer, dimension(:,:,:) :: bmasvegrow
      real, pointer, dimension(:,:,:) :: cmasvegcrow
      real, pointer, dimension(:,:,:) :: veghghtrow
      real, pointer, dimension(:,:,:) :: rootdpthrow
      real, pointer, dimension(:,:) :: rmlrow
      real, pointer, dimension(:,:) :: rmsrow
      real, pointer, dimension(:,:,:) :: tltrleafrow
      real, pointer, dimension(:,:,:) :: tltrstemrow
      real, pointer, dimension(:,:,:) :: tltrrootrow
      real, pointer, dimension(:,:,:) :: leaflitrrow
      real, pointer, dimension(:,:,:) :: roottemprow
      real, pointer, dimension(:,:,:) :: afrleafrow
      real, pointer, dimension(:,:,:) :: afrstemrow
      real, pointer, dimension(:,:,:) :: afrrootrow
      real, pointer, dimension(:,:,:) :: wtstatusrow
      real, pointer, dimension(:,:,:) :: ltstatusrow
      real, pointer, dimension(:,:) :: rmrrow

      real, pointer, dimension(:,:,:) :: slopefracrow
      real, pointer, dimension(:,:) :: ch4wet1row
      real, pointer, dimension(:,:) :: ch4wet2row
      real, pointer, dimension(:,:) :: wetfdynrow
      real, pointer, dimension(:,:) :: ch4dyn1row
      real, pointer, dimension(:,:) :: ch4dyn2row
      real, pointer, dimension(:,:,:) :: wetfrac_monrow
      real, pointer, dimension(:,:) :: ch4soillsrow

      real, pointer, dimension(:,:) :: lucemcomrow
      real, pointer, dimension(:,:) :: lucltrinrow
      real, pointer, dimension(:,:) :: lucsocinrow

      real, pointer, dimension(:,:) :: npprow
      real, pointer, dimension(:,:) :: neprow
      real, pointer, dimension(:,:) :: nbprow
      real, pointer, dimension(:,:) :: gpprow
      real, pointer, dimension(:,:) :: hetroresrow
      real, pointer, dimension(:,:) :: autoresrow
      real, pointer, dimension(:,:) :: soilcresprow
      real, pointer, dimension(:,:) :: rmrow
      real, pointer, dimension(:,:) :: rgrow
      real, pointer, dimension(:,:) :: litresrow
      real, pointer, dimension(:,:) :: socresrow
      real, pointer, dimension(:,:) :: dstcemlsrow
      real, pointer, dimension(:,:) :: litrfallrow
      real, pointer, dimension(:,:) :: humiftrsrow

      real, pointer, dimension(:,:,:) :: gppvegrow
      real, pointer, dimension(:,:,:) :: nepvegrow
      real, pointer, dimension(:,:,:) :: nbpvegrow
      real, pointer, dimension(:,:,:) :: nppvegrow
      real, pointer, dimension(:,:,:) :: hetroresvegrow
      real, pointer, dimension(:,:,:) :: autoresvegrow
      real, pointer, dimension(:,:,:) :: litresvegrow
      real, pointer, dimension(:,:,:) :: soilcresvegrow
      real, pointer, dimension(:,:,:) :: rmlvegaccrow
      real, pointer, dimension(:,:,:) :: rmsvegrow
      real, pointer, dimension(:,:,:) :: rmrvegrow
      real, pointer, dimension(:,:,:) :: rgvegrow
      real, pointer, dimension(:,:,:) :: litrfallvegrow
      real, pointer, dimension(:,:,:) :: humiftrsvegrow

      real, pointer, dimension(:,:,:) :: rothrlosrow
      real, pointer, dimension(:,:,:) :: pfcancmxrow
      real, pointer, dimension(:,:,:) :: nfcancmxrow
      real, pointer, dimension(:,:,:) :: alvsctmrow
      real, pointer, dimension(:,:,:) :: paicrow
      real, pointer, dimension(:,:,:) :: slaicrow
      real, pointer, dimension(:,:,:) :: alirctmrow
      real, pointer, dimension(:,:) :: cfluxcgrow
      real, pointer, dimension(:,:) :: cfluxcsrow
      real, pointer, dimension(:,:) :: dstcemls3row
      real, pointer, dimension(:,:,:) :: anvegrow
      real, pointer, dimension(:,:,:) :: rmlvegrow

      real, pointer, dimension(:,:) :: twarmmrow
      real, pointer, dimension(:,:) :: tcoldmrow
      real, pointer, dimension(:,:) :: gdd5row
      real, pointer, dimension(:,:) :: aridityrow
      real, pointer, dimension(:,:) :: srplsmonrow
      real, pointer, dimension(:,:) :: defctmonrow
      real, pointer, dimension(:,:) :: anndefctrow
      real, pointer, dimension(:,:) :: annsrplsrow
      real, pointer, dimension(:,:) :: annpcprow
      real, pointer, dimension(:,:) :: dry_season_lengthrow


      ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.
      ! GAT version:

      logical, pointer, dimension(:,:) :: pftexistgat
      integer, pointer, dimension(:,:) :: colddaysgat
      integer, pointer, dimension(:) :: icountgat
      integer, pointer, dimension(:,:) :: lfstatusgat
      integer, pointer, dimension(:,:) :: pandaysgat
      integer, pointer, dimension(:) :: stdalngat
      real, pointer, dimension(:) :: lightng

      real, pointer, dimension(:,:) :: ailcmingat         !
      real, pointer, dimension(:,:) :: ailcmaxgat         !
      real, pointer, dimension(:,:) :: dvdfcangat         !
      real, pointer, dimension(:,:) :: gleafmasgat        !
      real, pointer, dimension(:,:) :: bleafmasgat        !
      real, pointer, dimension(:,:) :: stemmassgat        !
      real, pointer, dimension(:,:) :: rootmassgat        !
      real, pointer, dimension(:,:) :: pstemmassgat       !
      real, pointer, dimension(:,:) :: pgleafmassgat      !
      real, pointer, dimension(:,:) :: fcancmxgat
      real, pointer, dimension(:) :: gavglaigat
      real, pointer, dimension(:,:) :: zolncgat
      real, pointer, dimension(:,:) :: ailcgat
      real, pointer, dimension(:,:) :: ailcggat
      real, pointer, dimension(:,:) :: ailcgsgat
      real, pointer, dimension(:,:) :: fcancsgat
      real, pointer, dimension(:,:) :: fcancgat
      real, pointer, dimension(:) :: co2concgat
      real, pointer, dimension(:) :: ch4concgat
      real, pointer, dimension(:,:) :: co2i1cggat
      real, pointer, dimension(:,:) :: co2i1csgat
      real, pointer, dimension(:,:) :: co2i2cggat
      real, pointer, dimension(:,:) :: co2i2csgat
      real, pointer, dimension(:,:) :: ancsveggat
      real, pointer, dimension(:,:) :: ancgveggat
      real, pointer, dimension(:,:) :: rmlcsveggat
      real, pointer, dimension(:,:) :: rmlcgveggat
      real, pointer, dimension(:,:) :: slaigat
      real, pointer, dimension(:,:) :: ailcbgat
      real, pointer, dimension(:) :: canresgat
      real, pointer, dimension(:,:) :: flhrlossgat

      real, pointer, dimension(:,:) :: grwtheffgat
      real, pointer, dimension(:,:) :: lystmmasgat
      real, pointer, dimension(:,:) :: lyrotmasgat
      real, pointer, dimension(:,:) :: tymaxlaigat
      real, pointer, dimension(:) :: vgbiomasgat
      real, pointer, dimension(:) :: gavgltmsgat
      real, pointer, dimension(:) :: gavgscmsgat
      real, pointer, dimension(:,:) :: stmhrlosgat
      real, pointer, dimension(:,:,:) :: rmatcgat
      real, pointer, dimension(:,:,:) :: rmatctemgat
      real, pointer, dimension(:,:) :: litrmassgat
      real, pointer, dimension(:,:) :: soilcmasgat
      real, pointer, dimension(:,:) :: vgbiomas_veggat

      real, pointer, dimension(:,:) :: emit_co2gat
      real, pointer, dimension(:,:) :: emit_cogat
      real, pointer, dimension(:,:) :: emit_ch4gat
      real, pointer, dimension(:,:) :: emit_nmhcgat
      real, pointer, dimension(:,:) :: emit_h2gat
      real, pointer, dimension(:,:) :: emit_noxgat
      real, pointer, dimension(:,:) :: emit_n2ogat
      real, pointer, dimension(:,:) :: emit_pm25gat
      real, pointer, dimension(:,:) :: emit_tpmgat
      real, pointer, dimension(:,:) :: emit_tcgat
      real, pointer, dimension(:,:) :: emit_ocgat
      real, pointer, dimension(:,:) :: emit_bcgat
      real, pointer, dimension(:) :: burnfracgat
      real, pointer, dimension(:,:) :: burnvegfgat
      real, pointer, dimension(:,:) :: smfuncveggat
      real, pointer, dimension(:) :: popdingat
      real, pointer, dimension(:,:) :: btermgat
      real, pointer, dimension(:) :: ltermgat
      real, pointer, dimension(:,:) :: mtermgat

      real, pointer, dimension(:) :: extnprobgat
      real, pointer, dimension(:) :: prbfrhucgat
      real, pointer, dimension(:,:) :: mlightnggat
      real, pointer, dimension(:) :: dayl_maxgat
      real, pointer, dimension(:) :: daylgat

      real, pointer, dimension(:,:) :: bmasveggat
      real, pointer, dimension(:,:) :: cmasvegcgat
      real, pointer, dimension(:,:) :: veghghtgat
      real, pointer, dimension(:,:) :: rootdpthgat
      real, pointer, dimension(:) :: rmlgat
      real, pointer, dimension(:) :: rmsgat
      real, pointer, dimension(:,:) :: tltrleafgat
      real, pointer, dimension(:,:) :: tltrstemgat
      real, pointer, dimension(:,:) :: tltrrootgat
      real, pointer, dimension(:,:) :: leaflitrgat
      real, pointer, dimension(:,:) :: roottempgat
      real, pointer, dimension(:,:) :: afrleafgat
      real, pointer, dimension(:,:) :: afrstemgat
      real, pointer, dimension(:,:) :: afrrootgat
      real, pointer, dimension(:,:) :: wtstatusgat
      real, pointer, dimension(:,:) :: ltstatusgat
      real, pointer, dimension(:) :: rmrgat

      real, pointer, dimension(:,:) :: slopefracgat
      real, pointer, dimension(:) :: wetfrac_presgat
      real, pointer, dimension(:,:) :: wetfrac_mongat
      real, pointer, dimension(:) :: ch4wet1gat
      real, pointer, dimension(:) :: ch4wet2gat
      real, pointer, dimension(:) :: wetfdyngat
      real, pointer, dimension(:) :: ch4dyn1gat
      real, pointer, dimension(:) :: ch4dyn2gat
      real, pointer, dimension(:) :: ch4soillsgat

      real, pointer, dimension(:) :: lucemcomgat
      real, pointer, dimension(:) :: lucltringat
      real, pointer, dimension(:) :: lucsocingat

      real, pointer, dimension(:) :: nppgat
      real, pointer, dimension(:) :: nepgat
      real, pointer, dimension(:) :: nbpgat
      real, pointer, dimension(:) :: gppgat
      real, pointer, dimension(:) :: hetroresgat
      real, pointer, dimension(:) :: autoresgat
      real, pointer, dimension(:) :: soilcrespgat
      real, pointer, dimension(:) :: rmgat
      real, pointer, dimension(:) :: rggat
      real, pointer, dimension(:) :: litresgat
      real, pointer, dimension(:) :: socresgat
      real, pointer, dimension(:) :: dstcemlsgat
      real, pointer, dimension(:) :: litrfallgat
      real, pointer, dimension(:) :: humiftrsgat

      real, pointer, dimension(:,:) :: gppveggat
      real, pointer, dimension(:,:) :: nepveggat
      real, pointer, dimension(:,:) :: nbpveggat
      real, pointer, dimension(:,:) :: nppveggat
      real, pointer, dimension(:,:) :: hetroresveggat
      real, pointer, dimension(:,:) :: autoresveggat
      real, pointer, dimension(:,:) :: litresveggat
      real, pointer, dimension(:,:) :: soilcresveggat
      real, pointer, dimension(:,:) :: rmlvegaccgat
      real, pointer, dimension(:,:) :: rmsveggat
      real, pointer, dimension(:,:) :: rmrveggat
      real, pointer, dimension(:,:) :: rgveggat
      real, pointer, dimension(:,:) :: litrfallveggat
      real, pointer, dimension(:,:) :: humiftrsveggat

      real, pointer, dimension(:,:) :: rothrlosgat
      real, pointer, dimension(:,:) :: pfcancmxgat
      real, pointer, dimension(:,:) :: nfcancmxgat
      real, pointer, dimension(:,:) :: alvsctmgat
      real, pointer, dimension(:,:) :: paicgat
      real, pointer, dimension(:,:) :: slaicgat
      real, pointer, dimension(:,:) :: alirctmgat
      real, pointer, dimension(:) :: cfluxcggat
      real, pointer, dimension(:) :: cfluxcsgat
      real, pointer, dimension(:) :: dstcemls3gat
      real, pointer, dimension(:,:) :: anveggat
      real, pointer, dimension(:,:) :: rmlveggat

      real, pointer, dimension(:) :: twarmmgat
      real, pointer, dimension(:) :: tcoldmgat
      real, pointer, dimension(:) :: gdd5gat
      real, pointer, dimension(:) :: ariditygat
      real, pointer, dimension(:) :: srplsmongat
      real, pointer, dimension(:) :: defctmongat
      real, pointer, dimension(:) :: anndefctgat
      real, pointer, dimension(:) :: annsrplsgat
      real, pointer, dimension(:) :: annpcpgat
      real, pointer, dimension(:) :: dry_season_lengthgat

      real, pointer, dimension(:) :: tcurm
      real, pointer, dimension(:) :: srpcuryr
      real, pointer, dimension(:) :: dftcuryr
      real, pointer, dimension(:,:) :: tmonth
      real, pointer, dimension(:) :: anpcpcur
      real, pointer, dimension(:) :: anpecur
      real, pointer, dimension(:) :: gdd5cur
      real, pointer, dimension(:) :: surmncur
      real, pointer, dimension(:) :: defmncur
      real, pointer, dimension(:) :: srplscur
      real, pointer, dimension(:) :: defctcur

      real, pointer, dimension(:,:) :: geremortgat
      real, pointer, dimension(:,:) :: intrmortgat
      real, pointer, dimension(:,:) :: lambdagat
      real, pointer, dimension(:,:) :: ccgat
      real, pointer, dimension(:,:) :: mmgat
      integer, pointer, dimension(:) :: altotcntr_d

      ! Mosaic level:

      real, pointer, dimension(:,:) :: PREACC_M
      real, pointer, dimension(:,:) :: GTACC_M
      real, pointer, dimension(:,:) :: QEVPACC_M
      real, pointer, dimension(:,:) :: HFSACC_M
      real, pointer, dimension(:,:) :: HMFNACC_M
      real, pointer, dimension(:,:) :: ROFACC_M
      real, pointer, dimension(:,:) :: SNOACC_M
      real, pointer, dimension(:,:) :: OVRACC_M
      real, pointer, dimension(:,:) :: WTBLACC_M
      real, pointer, dimension(:,:,:) :: TBARACC_M
      real, pointer, dimension(:,:,:) :: THLQACC_M
      real, pointer, dimension(:,:,:) :: THICACC_M
      real, pointer, dimension(:,:,:) :: THALACC_M
      real, pointer, dimension(:,:) :: ALVSACC_M
      real, pointer, dimension(:,:) :: ALIRACC_M
      real, pointer, dimension(:,:) :: RHOSACC_M
      real, pointer, dimension(:,:) :: TSNOACC_M
      real, pointer, dimension(:,:) :: WSNOACC_M
      real, pointer, dimension(:,:) :: SNOARE_M
      real, pointer, dimension(:,:) :: TCANACC_M
      real, pointer, dimension(:,:) :: RCANACC_M
      real, pointer, dimension(:,:) :: SCANACC_M
      real, pointer, dimension(:,:) :: GROACC_M
      real, pointer, dimension(:,:) :: FSINACC_M
      real, pointer, dimension(:,:) :: FLINACC_M
      real, pointer, dimension(:,:) :: TAACC_M
      real, pointer, dimension(:,:) :: UVACC_M
      real, pointer, dimension(:,:) :: PRESACC_M
      real, pointer, dimension(:,:) :: QAACC_M
      real, pointer, dimension(:,:) :: ALTOTACC_M
      real, pointer, dimension(:,:) :: EVAPACC_M
      real, pointer, dimension(:,:) :: FLUTACC_M

!      Outputs

       real, pointer, dimension(:,:) :: tcanoaccrow_out
       real, pointer, dimension(:) :: tcanoaccgat_out
       real, pointer, dimension(:,:) :: qevpacc_m_save

!     -----------------------
!      Tile-level variables (denoted by an ending of "_t")

      real, pointer, dimension(:) :: fsnowacc_t
      real, pointer, dimension(:) :: tcansacc_t
      real, pointer, dimension(:) :: tcanoaccgat_t
      real, pointer, dimension(:) :: taaccgat_t
      real, pointer, dimension(:) :: uvaccgat_t
      real, pointer, dimension(:) :: vvaccgat_t
      real, pointer, dimension(:,:) :: tbaraccgat_t
      real, pointer, dimension(:,:) :: tbarcacc_t
      real, pointer, dimension(:,:) :: tbarcsacc_t
      real, pointer, dimension(:,:) :: tbargacc_t
      real, pointer, dimension(:,:) :: tbargsacc_t
      real, pointer, dimension(:,:) :: thliqcacc_t
      real, pointer, dimension(:,:) :: thliqgacc_t
      real, pointer, dimension(:,:) :: thliqacc_t
      real, pointer, dimension(:,:) :: thicecacc_t
      real, pointer, dimension(:,:) :: thicegacc_t
      real, pointer, dimension(:,:) :: ancsvgac_t
      real, pointer, dimension(:,:) :: ancgvgac_t
      real, pointer, dimension(:,:) :: rmlcsvga_t
      real, pointer, dimension(:,:) :: rmlcgvga_t

!     -----------------------
!     Grid-averaged variables (denoted with an ending of "_g")

      real, pointer, dimension(:) ::  fsstar_g
      real, pointer, dimension(:) ::  flstar_g
      real, pointer, dimension(:) ::  qh_g
      real, pointer, dimension(:) ::  qe_g
      real, pointer, dimension(:) ::  snomlt_g
      real, pointer, dimension(:) ::  beg_g
      real, pointer, dimension(:) ::  gtout_g
      real, pointer, dimension(:) ::  tpn_g
      real, pointer, dimension(:) ::  altot_g
      real, pointer, dimension(:) ::  tcn_g
      real, pointer, dimension(:) ::  tsn_g
      real, pointer, dimension(:) ::  zsn_g

      real, pointer, dimension(:) :: WSNOROT_g
      real, pointer, dimension(:) :: ROFSROT_g
      real, pointer, dimension(:) :: SNOROT_g
      real, pointer, dimension(:) :: RHOSROT_g
      real, pointer, dimension(:) :: ROFROT_g
      real, pointer, dimension(:) :: ZPNDROT_g
      real, pointer, dimension(:) :: RCANROT_g
      real, pointer, dimension(:) :: SCANROT_g
      real, pointer, dimension(:) :: TROFROT_g
      real, pointer, dimension(:) :: TROOROT_g
      real, pointer, dimension(:) :: TROBROT_g
      real, pointer, dimension(:) :: ROFOROT_g
      real, pointer, dimension(:) :: ROFBROT_g
      real, pointer, dimension(:) :: TROSROT_g
      real, pointer, dimension(:) :: FSGVROT_g
      real, pointer, dimension(:) :: FSGSROT_g
      real, pointer, dimension(:) :: FLGVROT_g
      real, pointer, dimension(:) :: FLGSROT_g
      real, pointer, dimension(:) :: HFSCROT_g
      real, pointer, dimension(:) :: HFSSROT_g
      real, pointer, dimension(:) :: HEVCROT_g
      real, pointer, dimension(:) :: HEVSROT_g
      real, pointer, dimension(:) :: HMFCROT_g
      real, pointer, dimension(:) :: HMFNROT_g
      real, pointer, dimension(:) :: HTCSROT_g
      real, pointer, dimension(:) :: HTCCROT_g
      real, pointer, dimension(:) :: FSGGROT_g
      real, pointer, dimension(:) :: FLGGROT_g
      real, pointer, dimension(:) :: HFSGROT_g
      real, pointer, dimension(:) :: HEVGROT_g
      real, pointer, dimension(:) :: CDHROT_g
      real, pointer, dimension(:) :: CDMROT_g
      real, pointer, dimension(:) :: SFCUROT_g
      real, pointer, dimension(:) :: SFCVROT_g
      real, pointer, dimension(:) :: ACTLYR_g
      real, pointer, dimension(:) :: FTABLE_g
      real, pointer, dimension(:) :: fc_g
      real, pointer, dimension(:) :: fg_g
      real, pointer, dimension(:) :: fcs_g
      real, pointer, dimension(:) :: fgs_g
      real, pointer, dimension(:) :: PCFCROT_g
      real, pointer, dimension(:) :: PCLCROT_g
      real, pointer, dimension(:) :: PCPGROT_g
      real, pointer, dimension(:) :: QFCFROT_g
      real, pointer, dimension(:) :: QFGROT_g
      real, pointer, dimension(:,:) :: QFCROT_g
      real, pointer, dimension(:) :: ROFCROT_g
      real, pointer, dimension(:) :: ROFNROT_g
      real, pointer, dimension(:) :: WTRSROT_g
      real, pointer, dimension(:) :: WTRGROT_g
      real, pointer, dimension(:) :: PCPNROT_g
      real, pointer, dimension(:) :: QFCLROT_g
      real, pointer, dimension(:) :: QFNROT_g
      real, pointer, dimension(:) :: WTRCROT_g
      real, pointer, dimension(:,:) :: rmlvegrow_g
      real, pointer, dimension(:,:) :: anvegrow_g
      real, pointer, dimension(:,:) :: HMFGROT_g
      real, pointer, dimension(:,:) :: HTCROT_g
      real, pointer, dimension(:,:) :: TBARROT_g
      real, pointer, dimension(:,:) :: THLQROT_g
      real, pointer, dimension(:,:) :: THICROT_g
      real, pointer, dimension(:,:) :: GFLXROT_g

    ! Model Switches (rarely changed ones only! The rest are in joboptions file):

      logical, parameter :: obslght = .false.  ! if true the observed lightning will be used. False means you will use the
                                             ! lightning climatology from the CTM file. This was brought in for FireMIP runs.

    ! If you intend to have LUC BETWEEN tiles then set this to true:
      logical, parameter ::  onetile_perPFT = .False. ! NOTE: This is usually not the behaviour desired unless you are
                                                   ! running with one PFT on each tile and want them to compete for space
                                                   ! across tiles. In general keep this as False. JM Feb 2016.
!
!============= CTEM array declaration done =============================/
!
!=======================================================================
!     * PHYSICAL CONSTANTS.
!     * PARAMETERS IN THE FOLLOWING COMMON BLOCKS ARE NORMALLY DEFINED
!     * WITHIN THE GCM.

      COMMON /PARAMS/ X1,    X2,    X3,    X4,   G,GAS,   X5,&
     &                X6,    CPRES, GASV,  X7
      COMMON /PARAM1/ CPI,   X8,    CELZRO,X9,    X10,    X11
      COMMON /PARAM3/ X12,   X13,   X14,   X15,   SIGMA,  X16
      COMMON  /TIMES/ DELTIM,K1,    K2,    K3,    K4,     K5,&
     &                K6,    K7,    K8,    K9,    K10,    K11
!
!     * THE FOLLOWING COMMON BLOCKS ARE DEFINED SPECIFICALLY FOR USE
!     * IN CLASS, VIA BLOCK DATA AND THE SUBROUTINE "CLASSD".
!
      COMMON /CLASS1/ DELT,TFREZ
      COMMON /CLASS2/ RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN
      COMMON /CLASS3/ TCW,TCICE,TCSAND,TCCLAY,TCOM,TCDRYS,&
     &                RHOSOL,RHOOM
      COMMON /CLASS4/ HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,&
     &                SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,&
     &                TCGLAC,CLHMLT,CLHVAP
      COMMON /CLASS5/ THPORG,THRORG,THMORG,BORG,PSISORG,GRKSORG
      COMMON /CLASS6/ PI,GROWYR,ZOLNG,ZOLNS,ZOLNI,ZORAT,ZORATG
      COMMON /CLASS7/ CANEXT,XLEAF
      COMMON /CLASS8/ ALVSI,ALIRI,ALVSO,ALIRO,ALBRCK
      COMMON /PHYCON/ DELTA,CGRAV,CKARM,CPD
      COMMON /CLASSD2/ AS,ASX,CI,BS,BETA,FACTN,HMIN,ANGMAX
!
!===================== CTEM ==============================================\
    allocate(class_gat)
    allocate(class_rot)


    ! Point the CLASS pointers

    ILMOS   => class_gat%ILMOS
    JLMOS   => class_gat%JLMOS
    IWMOS   => class_gat%IWMOS
    JWMOS   => class_gat%JWMOS
    IGDRGAT => class_gat%IGDRGAT
    DELZ    => class_gat% DELZ
    ZBOT    => class_gat% ZBOT
    ALBSGAT => class_gat% ALBSGAT
    CMAIGAT => class_gat% CMAIGAT
    GROGAT  => class_gat% GROGAT
    QACGAT  => class_gat% QACGAT
    RCANGAT => class_gat%RCANGAT
    RHOSGAT=> class_gat%RHOSGAT
    SCANGAT=> class_gat%SCANGAT
    SNOGAT=> class_gat%SNOGAT
    TACGAT=> class_gat%TACGAT
    TBASGAT=> class_gat%TBASGAT
    TCANGAT=> class_gat%TCANGAT
    TPNDGAT=> class_gat%TPNDGAT
    TSNOGAT=> class_gat%TSNOGAT
    WSNOGAT=> class_gat%WSNOGAT
    ZPNDGAT=> class_gat%ZPNDGAT
    REFGAT=> class_gat%REFGAT
    BCSNGAT => class_gat%BCSNGAT
    AGIDGAT => class_gat%AGIDGAT
    AGVDGAT => class_gat%AGVDGAT
    ALGDGAT => class_gat%ALGDGAT
    ALGWGAT => class_gat%ALGWGAT
    ASIDGAT => class_gat%ASIDGAT
    ASVDGAT => class_gat%ASVDGAT
    DRNGAT => class_gat%DRNGAT
    GRKFGAT => class_gat%GRKFGAT
    WFCIGAT => class_gat%WFCIGAT
    WFSFGAT => class_gat%WFSFGAT
    XSLPGAT => class_gat%XSLPGAT
    ZPLGGAT => class_gat%ZPLGGAT
    ZPLSGAT => class_gat%ZPLSGAT
    ZSNLGAT => class_gat%ZSNLGAT
    ALGWVGAT => class_gat%ALGWVGAT
    ALGWNGAT => class_gat%ALGWNGAT
    ALGDVGAT => class_gat%ALGDVGAT
    ALGDNGAT => class_gat%ALGDNGAT
    EMISGAT => class_gat%EMISGAT
    CSZGAT => class_gat%CSZGAT
    DLONGAT => class_gat%DLONGAT
    DLATGAT => class_gat%DLATGAT
    FCLOGAT => class_gat%FCLOGAT
    FDLGAT => class_gat%FDLGAT
    FSIHGAT => class_gat%FSIHGAT
    FSVHGAT => class_gat%FSVHGAT
    GGEOGAT => class_gat%GGEOGAT
    PADRGAT => class_gat%PADRGAT
    PREGAT => class_gat%PREGAT
    PRESGAT => class_gat%PRESGAT
    QAGAT => class_gat%QAGAT
    RADJGAT => class_gat%RADJGAT
    RHOAGAT => class_gat%RHOAGAT
    RHSIGAT => class_gat%RHSIGAT
    RPCPGAT => class_gat%RPCPGAT
    SPCPGAT => class_gat%SPCPGAT
    TAGAT => class_gat%TAGAT
    TADPGAT => class_gat%TADPGAT
    TRPCGAT => class_gat%TRPCGAT
    TSPCGAT => class_gat%TSPCGAT
    ULGAT => class_gat%ULGAT
    VLGAT => class_gat%VLGAT
    VMODGAT => class_gat%VMODGAT
    VPDGAT => class_gat%VPDGAT
    Z0ORGAT => class_gat%Z0ORGAT
    ZBLDGAT => class_gat%ZBLDGAT
    ZDHGAT => class_gat%ZDHGAT
    ZDMGAT => class_gat%ZDMGAT
    ZRFHGAT => class_gat%ZRFHGAT
    ZRFMGAT => class_gat%ZRFMGAT
    FSGGAT => class_gat%FSGGAT
    FLGGAT => class_gat%FLGGAT
    GUSTGAT => class_gat%GUSTGAT
    DEPBGAT => class_gat%DEPBGAT
    GTBS => class_gat%GTBS
    SFCUBS => class_gat%SFCUBS
    SFCVBS => class_gat%SFCVBS
    USTARBS => class_gat%USTARBS
    TCSNOW => class_gat%TCSNOW
    GSNOW => class_gat%GSNOW
    ALIRGAT => class_gat%ALIRGAT
    ALVSGAT => class_gat%ALVSGAT
    CDHGAT => class_gat%CDHGAT
    CDMGAT => class_gat%CDMGAT
    DRGAT => class_gat%DRGAT
    EFGAT => class_gat%EFGAT
    FLGGGAT => class_gat%FLGGGAT
    FLGSGAT => class_gat%FLGSGAT
    FLGVGAT => class_gat%FLGVGAT
    FSGGGAT => class_gat%FSGGGAT
    FSGSGAT => class_gat% FSGSGAT
    FSGVGAT => class_gat% FSGVGAT
    FSNOGAT => class_gat% FSNOGAT
    GAGAT => class_gat% GAGAT
    GTGAT => class_gat% GTGAT
    HBLGAT => class_gat% HBLGAT
    HEVCGAT => class_gat% HEVCGAT
    HEVGGAT => class_gat% HEVGGAT
    HEVSGAT => class_gat% HEVSGAT
    HFSGAT => class_gat% HFSGAT
    HFSCGAT => class_gat% HFSCGAT
    HFSGGAT => class_gat% HFSGGAT
    HFSSGAT => class_gat% HFSSGAT
    HMFCGAT => class_gat% HMFCGAT
    HMFNGAT => class_gat% HMFNGAT
    HTCCGAT => class_gat% HTCCGAT
    HTCSGAT => class_gat% HTCSGAT
    ILMOGAT => class_gat% ILMOGAT
    PCFCGAT => class_gat% PCFCGAT
    PCLCGAT => class_gat% PCLCGAT
    PCPGGAT => class_gat% PCPGGAT
    PCPNGAT => class_gat% PCPNGAT
    PETGAT => class_gat% PETGAT
    QEVPGAT => class_gat% QEVPGAT
    QFCFGAT => class_gat% QFCFGAT
    QFCLGAT => class_gat% QFCLGAT
    QFGGAT => class_gat% QFGGAT
    QFNGAT => class_gat% QFNGAT
    QFSGAT => class_gat% QFSGAT
    QFXGAT => class_gat% QFXGAT
    QGGAT => class_gat% QGGAT
    ROFGAT => class_gat% ROFGAT
    ROFBGAT => class_gat% ROFBGAT
    ROFCGAT => class_gat% ROFCGAT
    ROFNGAT => class_gat% ROFNGAT
    ROFOGAT => class_gat% ROFOGAT
    ROFSGAT => class_gat% ROFSGAT
    ROVGGAT => class_gat% ROVGGAT
    SFCQGAT => class_gat% SFCQGAT
    SFCTGAT => class_gat% SFCTGAT
    SFCUGAT => class_gat% SFCUGAT
    SFCVGAT => class_gat% SFCVGAT
    TFXGAT => class_gat% TFXGAT
    TROBGAT => class_gat% TROBGAT
    TROFGAT => class_gat% TROFGAT
    TROOGAT => class_gat% TROOGAT
    TROSGAT => class_gat% TROSGAT
    UEGAT => class_gat% UEGAT
    WTABGAT => class_gat% WTABGAT
    WTRCGAT => class_gat% WTRCGAT
    WTRGGAT => class_gat% WTRGGAT
    WTRSGAT => class_gat% WTRSGAT
    QLWOGAT => class_gat% QLWOGAT
    SFRHGAT => class_gat% SFRHGAT
    FTEMP => class_gat% FTEMP
    FVAP => class_gat% FVAP
    RIB => class_gat% RIB
    FC => class_gat% FC
    FG => class_gat% FG
    FCS => class_gat% FCS
    FGS => class_gat% FGS
    RBCOEF => class_gat% RBCOEF
    ZSNOW => class_gat% ZSNOW
    FSVF => class_gat% FSVF
    FSVFS => class_gat% FSVFS
    ALVSCN => class_gat% ALVSCN
    ALIRCN => class_gat% ALIRCN
    ALVSG => class_gat% ALVSG
    ALIRG => class_gat% ALIRG
    ALVSCS => class_gat% ALVSCS
    ALIRCS => class_gat% ALIRCS
    ALVSSN => class_gat% ALVSSN
    ALIRSN => class_gat% ALIRSN
    ALVSGC => class_gat% ALVSGC
    ALIRGC => class_gat% ALIRGC
    ALVSSC => class_gat% ALVSSC
    ALIRSC => class_gat% ALIRSC
    TRVSCN => class_gat% TRVSCN
    TRIRCN => class_gat% TRIRCN
    TRVSCS => class_gat%TRVSCS
    TRIRCS => class_gat%TRIRCS
    RC => class_gat%RC
    RCS => class_gat%RCS
    FRAINC => class_gat%FRAINC
    FSNOWC => class_gat%FSNOWC
    FRAICS => class_gat%FRAICS
    FSNOCS => class_gat%FSNOCS
    CMASSC => class_gat%CMASSC
    CMASCS => class_gat%CMASCS
    DISP => class_gat%DISP
    DISPS => class_gat%DISPS
    ZOMLNC => class_gat%ZOMLNC
    ZOELNC => class_gat%ZOELNC
    ZOMLNG => class_gat%ZOMLNG
    ZOELNG => class_gat%ZOELNG
    ZOMLCS => class_gat%ZOMLCS
    ZOELCS => class_gat%ZOELCS
    ZOMLNS => class_gat%ZOMLNS
    ZOELNS => class_gat%ZOELNS
    TRSNOWC => class_gat%TRSNOWC
    CHCAP => class_gat%CHCAP
    CHCAPS => class_gat%CHCAPS
    GZEROC => class_gat%GZEROC
    GZEROG => class_gat%GZEROG
    GZROCS => class_gat%GZROCS
    GZROGS => class_gat%GZROGS
    G12C => class_gat%G12C
    G12G => class_gat%G12G
    G12CS => class_gat%G12CS
    G12GS => class_gat%G12GS
    G23C => class_gat%G23C
    G23G => class_gat%G23G
    G23CS => class_gat%G23CS
    G23GS => class_gat%G23GS
    QFREZC => class_gat%QFREZC
    QFREZG => class_gat%QFREZG
    QMELTC => class_gat%QMELTC
    QMELTG => class_gat%QMELTG
    EVAPC => class_gat%EVAPC
    EVAPCG => class_gat%EVAPCG
    EVAPG => class_gat%EVAPG
    EVAPCS => class_gat%EVAPCS
    EVPCSG => class_gat%EVPCSG
    EVAPGS => class_gat%EVAPGS
    TCANO => class_gat%TCANO
    TCANS => class_gat%TCANS
    RAICAN => class_gat%RAICAN
    SNOCAN => class_gat%SNOCAN
    RAICNS => class_gat%RAICNS
    SNOCNS => class_gat%SNOCNS
    CWLCAP => class_gat%CWLCAP
    CWFCAP => class_gat%CWFCAP
    CWLCPS => class_gat%CWLCPS
    CWFCPS => class_gat%CWFCPS
    TSNOCS => class_gat%TSNOCS
    TSNOGS => class_gat%TSNOGS
    RHOSCS => class_gat%RHOSCS
    RHOSGS => class_gat%RHOSGS
    WSNOCS => class_gat%WSNOCS
    WSNOGS => class_gat%WSNOGS
    TPONDC => class_gat%TPONDC
    TPONDG => class_gat%TPONDG
    TPNDCS => class_gat%TPNDCS
    TPNDGS => class_gat%TPNDGS
    ZPLMCS => class_gat%ZPLMCS
    ZPLMGS => class_gat%ZPLMGS
    ZPLIMC => class_gat%ZPLIMC
    ZPLIMG => class_gat%ZPLIMG
    CTVSTP => class_gat%CTVSTP
    CTSSTP => class_gat%CTSSTP
    CT1STP => class_gat%CT1STP
    CT2STP => class_gat%CT2STP
    CT3STP => class_gat%CT3STP
    WTVSTP => class_gat%WTVSTP
    WTSSTP => class_gat%WTSSTP
    WTGSTP => class_gat%WTGSTP
    ISNDGAT => class_gat% ISNDGAT
    TBARGAT => class_gat% TBARGAT
    THICGAT => class_gat% THICGAT
    THLQGAT => class_gat% THLQGAT
    BIGAT => class_gat% BIGAT
    DLZWGAT => class_gat% DLZWGAT
    GRKSGAT => class_gat% GRKSGAT
    HCPSGAT => class_gat% HCPSGAT
    PSISGAT => class_gat% PSISGAT
    PSIWGAT => class_gat% PSIWGAT
    TCSGAT => class_gat% TCSGAT
    THFCGAT => class_gat% THFCGAT
    THMGAT => class_gat% THMGAT
    THPGAT => class_gat% THPGAT
    THRGAT => class_gat% THRGAT
    THRAGAT => class_gat% THRAGAT
    ZBTWGAT => class_gat% ZBTWGAT
    THLWGAT => class_gat% THLWGAT
    GFLXGAT => class_gat% GFLXGAT
    HMFGGAT => class_gat% HMFGGAT
    HTCGAT => class_gat% HTCGAT
    QFCGAT => class_gat% QFCGAT
    TBARC => class_gat% TBARC
    TBARG => class_gat% TBARG
    TBARCS => class_gat% TBARCS
    TBARGS => class_gat% TBARGS
    THLIQC => class_gat% THLIQC
    THLIQG => class_gat% THLIQG
    THICEC => class_gat% THICEC
    THICEG => class_gat% THICEG
    FROOT => class_gat% FROOT
    HCPC => class_gat% HCPC
    HCPG => class_gat% HCPG
    FROOTS => class_gat% FROOTS
    TCTOPC => class_gat% TCTOPC
    TCBOTC => class_gat% TCBOTC
    TCTOPG => class_gat% TCTOPG
    TCBOTG => class_gat% TCBOTG
    ACIDGAT => class_gat% ACIDGAT
    ACVDGAT => class_gat% ACVDGAT
    CMASGAT => class_gat% CMASGAT
    HGTDGAT => class_gat% HGTDGAT
    PAIDGAT => class_gat% PAIDGAT
    PAMNGAT => class_gat% PAMNGAT
    PAMXGAT => class_gat% PAMXGAT
    PSGAGAT => class_gat% PSGAGAT
    PSGBGAT => class_gat% PSGBGAT
    QA50GAT => class_gat% QA50GAT
    ROOTGAT => class_gat% ROOTGAT
    RSMNGAT => class_gat% RSMNGAT
    VPDAGAT => class_gat% VPDAGAT
    VPDBGAT => class_gat% VPDBGAT
    ALICGAT => class_gat% ALICGAT
    ALVCGAT => class_gat% ALVCGAT
    FCANGAT => class_gat% FCANGAT
    LNZ0GAT => class_gat% LNZ0GAT
    FSDBGAT => class_gat% FSDBGAT
    FSFBGAT => class_gat% FSFBGAT
    FSSBGAT => class_gat% FSSBGAT
    SALBGAT => class_gat% SALBGAT
    CSALGAT => class_gat% CSALGAT
    ALTG => class_gat% ALTG
    ALSNO => class_gat% ALSNO
    TRSNOWG => class_gat% TRSNOWG
    TSFSGAT => class_gat% TSFSGAT
    ITCTGAT => class_gat% ITCTGAT

    ALIRACC => class_rot% ALIRACC
    ALVSACC => class_rot% ALVSACC
    EVAPACC => class_rot% EVAPACC
    FLINACC => class_rot% FLINACC
    FLUTACC => class_rot% FLUTACC
    FSINACC => class_rot% FSINACC
    GROACC => class_rot% GROACC
    GTACC => class_rot% GTACC
    HFSACC => class_rot% HFSACC
    HMFNACC => class_rot% HMFNACC
    OVRACC => class_rot% OVRACC
    PREACC => class_rot% PREACC
    PRESACC => class_rot% PRESACC
    QAACC => class_rot% QAACC
    QEVPACC => class_rot% QEVPACC
    RCANACC => class_rot% RCANACC
    RHOSACC => class_rot% RHOSACC
    ROFACC => class_rot% ROFACC
    SCANACC => class_rot% SCANACC
    SNOACC => class_rot% SNOACC
    TAACC => class_rot% TAACC
    TCANACC => class_rot% TCANACC
    TSNOACC => class_rot% TSNOACC
    UVACC => class_rot% UVACC
    WSNOACC => class_rot% WSNOACC
    WTBLACC => class_rot% WTBLACC
    ALTOTACC => class_rot% ALTOTACC
    CANARE => class_rot% CANARE
    SNOARE => class_rot% SNOARE
    CSZROW => class_rot% CSZROW
    DLONROW => class_rot% DLONROW
    DLATROW => class_rot% DLATROW
    FCLOROW => class_rot% FCLOROW
    FDLROW => class_rot% FDLROW
    FSIHROW => class_rot% FSIHROW
    FSVHROW => class_rot% FSVHROW
    GCROW => class_rot% GCROW
    GGEOROW => class_rot% GGEOROW
    PADRROW => class_rot% PADRROW
    PREROW => class_rot% PREROW
    PRESROW => class_rot% PRESROW
    QAROW => class_rot% QAROW
    RADJROW => class_rot% RADJROW
    RHOAROW => class_rot% RHOAROW
    RHSIROW => class_rot% RHSIROW
    RPCPROW => class_rot% RPCPROW
    RPREROW => class_rot% RPREROW
    SPCPROW => class_rot% SPCPROW
    SPREROW => class_rot% SPREROW
    TAROW => class_rot% TAROW
    TADPROW => class_rot% TADPROW
    TRPCROW => class_rot% TRPCROW
    TSPCROW => class_rot% TSPCROW
    ULROW => class_rot% ULROW
    VLROW => class_rot% VLROW
    VMODROW => class_rot% VMODROW
    VPDROW => class_rot% VPDROW
    ZBLDROW => class_rot% ZBLDROW
    ZDHROW => class_rot% ZDHROW
    ZDMROW => class_rot% ZDMROW
    ZRFHROW => class_rot% ZRFHROW
    ZRFMROW => class_rot% ZRFMROW
    UVROW => class_rot% UVROW
    XDIFFUS => class_rot% XDIFFUS
    Z0ORROW => class_rot% Z0ORROW
    FSSROW => class_rot% FSSROW
    PRENROW => class_rot% PRENROW
    CLDTROW => class_rot% CLDTROW
    FSGROL => class_rot% FSGROL
    FLGROL => class_rot% FLGROL
    GUSTROL => class_rot% GUSTROL
    DEPBROW => class_rot% DEPBROW
    ALIRROW => class_rot% ALIRROW
    ALVSROW => class_rot% ALVSROW
    CDHROW => class_rot% CDHROW
    CDMROW => class_rot% CDMROW
    DRROW => class_rot% DRROW
    EFROW => class_rot% EFROW
    FLGGROW => class_rot% FLGGROW
    FLGSROW => class_rot% FLGSROW
    FLGVROW => class_rot% FLGVROW
    FSGGROW => class_rot% FSGGROW
    FSGSROW => class_rot% FSGSROW
    FSGVROW => class_rot% FSGVROW
    FSNOROW => class_rot% FSNOROW
    GAROW => class_rot% GAROW
    GTROW => class_rot% GTROW
    HBLROW => class_rot% HBLROW
    HEVCROW => class_rot% HEVCROW
    HEVGROW => class_rot% HEVGROW
    HEVSROW => class_rot% HEVSROW
    HFSROW => class_rot% HFSROW
    HFSCROW => class_rot% HFSCROW
    HFSGROW => class_rot% HFSGROW
    HFSSROW => class_rot% HFSSROW
    HMFCROW => class_rot% HMFCROW
    HMFNROW => class_rot% HMFNROW
    HTCCROW => class_rot% HTCCROW
    HTCSROW => class_rot% HTCSROW
    ILMOROW => class_rot% ILMOROW
    PCFCROW => class_rot% PCFCROW
    PCLCROW => class_rot% PCLCROW
    PCPGROW => class_rot% PCPGROW
    PCPNROW => class_rot% PCPNROW
    PETROW => class_rot% PETROW
    QEVPROW => class_rot% QEVPROW
    QFCFROW => class_rot% QFCFROW
    QFCLROW => class_rot% QFCLROW
    QFGROW => class_rot% QFGROW
    QFNROW => class_rot% QFNROW
    QFSROW => class_rot% QFSROW
    QFXROW => class_rot% QFXROW
    QGROW => class_rot% QGROW
    ROFROW => class_rot% ROFROW
    ROFBROW => class_rot% ROFBROW
    ROFCROW => class_rot% ROFCROW
    ROFNROW => class_rot% ROFNROW
    ROFOROW => class_rot% ROFOROW
    ROFSROW => class_rot% ROFSROW
    ROVGROW => class_rot% ROVGROW
    SFCQROW => class_rot% SFCQROW
    SFCTROW => class_rot% SFCTROW
    SFCUROW => class_rot% SFCUROW
    SFCVROW => class_rot% SFCVROW
    TFXROW => class_rot% TFXROW
    UEROW => class_rot% UEROW
    WTABROW => class_rot% WTABROW
    WTRCROW => class_rot% WTRCROW
    WTRGROW => class_rot% WTRGROW
    WTRSROW => class_rot% WTRSROW
    SFRHROW => class_rot% SFRHROW
    IGDRROT => class_rot% IGDRROT
    MIDROT => class_rot% MIDROT
    ALBSROT => class_rot% ALBSROT
    CMAIROT => class_rot% CMAIROT
    GROROT => class_rot% GROROT
    QACROT => class_rot% QACROT
    RCANROT => class_rot% RCANROT
    RHOSROT => class_rot% RHOSROT
    SCANROT => class_rot% SCANROT
    SNOROT => class_rot% SNOROT
    TACROT => class_rot% TACROT
    TBASROT => class_rot% TBASROT
    TCANROT => class_rot% TCANROT
    TPNDROT => class_rot% TPNDROT
    TSNOROT => class_rot% TSNOROT
    WSNOROT => class_rot% WSNOROT
    ZPNDROT => class_rot% ZPNDROT
    REFROT => class_rot% REFROT
    BCSNROT => class_rot% BCSNROT
    AGIDROT => class_rot% AGIDROT
    AGVDROT => class_rot% AGVDROT
    ALGDROT => class_rot% ALGDROT
    ALGWROT => class_rot% ALGWROT
    ASIDROT => class_rot% ASIDROT
    ASVDROT => class_rot% ASVDROT
    DRNROT => class_rot% DRNROT
    FAREROT => class_rot% FAREROT
    GRKFROT => class_rot% GRKFROT
    WFCIROT => class_rot% WFCIROT
    WFSFROT => class_rot% WFSFROT
    XSLPROT => class_rot% XSLPROT
    ZPLGROT => class_rot% ZPLGROT
    ZPLSROT => class_rot% ZPLSROT
    ZSNLROT => class_rot% ZSNLROT
    ZSNOROT => class_rot% ZSNOROT
    ALGWVROT => class_rot% ALGWVROT
    ALGWNROT => class_rot% ALGWNROT
    ALGDVROT => class_rot% ALGDVROT
    ALGDNROT => class_rot% ALGDNROT
    EMISROT => class_rot% EMISROT
    ALIRROT => class_rot% ALIRROT
    ALVSROT => class_rot% ALVSROT
    CDHROT => class_rot% CDHROT
    CDMROT => class_rot% CDMROT
    DRROT => class_rot% DRROT
    EFROT => class_rot% EFROT
    FLGGROT => class_rot% FLGGROT
    FLGSROT => class_rot% FLGSROT
    FLGVROT => class_rot% FLGVROT
    FSGGROT => class_rot% FSGGROT
    FSGSROT => class_rot% FSGSROT
    FSGVROT => class_rot% FSGVROT
    FSNOROT => class_rot% FSNOROT
    GAROT => class_rot% GAROT
    GTROT => class_rot% GTROT
    HBLROT => class_rot% HBLROT
    HEVCROT => class_rot% HEVCROT
    HEVGROT => class_rot% HEVGROT
    HEVSROT => class_rot% HEVSROT
    HFSROT => class_rot% HFSROT
    HFSCROT => class_rot% HFSCROT
    HFSGROT => class_rot% HFSGROT
    HFSSROT => class_rot% HFSSROT
    HMFCROT => class_rot% HMFCROT
    HMFNROT => class_rot% HMFNROT
    HTCCROT => class_rot% HTCCROT
    SDEPROT => class_rot% SDEPROT
    SOCIROT => class_rot% SOCIROT
    HTCSROT => class_rot% HTCSROT
    ILMOROT => class_rot% ILMOROT
    PCFCROT => class_rot% PCFCROT
    PCLCROT => class_rot% PCLCROT
    PCPGROT => class_rot% PCPGROT
    PCPNROT => class_rot% PCPNROT
    PETROT => class_rot% PETROT
    QEVPROT => class_rot% QEVPROT
    QFCFROT => class_rot% QFCFROT
    QFCLROT => class_rot% QFCLROT
    QFGROT => class_rot% QFGROT
    QFNROT => class_rot% QFNROT
    QFSROT => class_rot% QFSROT
    QFXROT => class_rot% QFXROT
    QGROT => class_rot% QGROT
    ROFROT => class_rot% ROFROT
    ROFBROT => class_rot% ROFBROT
    ROFCROT => class_rot% ROFCROT
    ROFNROT => class_rot% ROFNROT
    ROFOROT => class_rot% ROFOROT
    ROFSROT => class_rot% ROFSROT
    ROVGROT => class_rot% ROVGROT
    SFCQROT => class_rot% SFCQROT
    SFCTROT => class_rot% SFCTROT
    SFCUROT => class_rot% SFCUROT
    SFCVROT => class_rot% SFCVROT
    TFXROT => class_rot% TFXROT
    TROBROT => class_rot% TROBROT
    TROFROT => class_rot% TROFROT
    TROOROT => class_rot% TROOROT
    TROSROT => class_rot% TROSROT
    UEROT => class_rot% UEROT
    WTABROT => class_rot% WTABROT
    WTRCROT => class_rot% WTRCROT
    WTRGROT => class_rot% WTRGROT
    WTRSROT => class_rot% WTRSROT
    SFRHROT => class_rot% SFRHROT
    ISNDROT => class_rot% ISNDROT
    TBARROT => class_rot% TBARROT
    THICROT => class_rot% THICROT
    THLQROT => class_rot% THLQROT
    BIROT => class_rot% BIROT
    DLZWROT => class_rot% DLZWROT
    GRKSROT => class_rot% GRKSROT
    HCPSROT => class_rot% HCPSROT
    SANDROT => class_rot% SANDROT
    CLAYROT => class_rot% CLAYROT
    ORGMROT => class_rot% ORGMROT
    PSISROT => class_rot% PSISROT
    PSIWROT => class_rot% PSIWROT
    TCSROT => class_rot% TCSROT
    THFCROT => class_rot% THFCROT
    THMROT => class_rot% THMROT
    THPROT => class_rot% THPROT
    THRROT => class_rot% THRROT
    THRAROT => class_rot% THRAROT
    ZBTWROT => class_rot% ZBTWROT
    THLWROT => class_rot% THLWROT
    GFLXROT => class_rot% GFLXROT
    HMFGROT => class_rot% HMFGROT
    HTCROT => class_rot% HTCROT
    QFCROT => class_rot% QFCROT
    ACIDROT => class_rot% ACIDROT
    ACVDROT => class_rot% ACVDROT
    CMASROT => class_rot% CMASROT
    HGTDROT => class_rot% HGTDROT
    PAIDROT => class_rot% PAIDROT
    PAMNROT => class_rot% PAMNROT
    PAMXROT => class_rot% PAMXROT
    PSGAROT => class_rot% PSGAROT
    PSGBROT => class_rot% PSGBROT
    QA50ROT => class_rot% QA50ROT
    ROOTROT => class_rot% ROOTROT
    RSMNROT => class_rot% RSMNROT
    VPDAROT => class_rot% VPDAROT
    VPDBROT => class_rot% VPDBROT
    ALICROT => class_rot% ALICROT
    ALVCROT => class_rot% ALVCROT
    FCANROT => class_rot% FCANROT
    LNZ0ROT => class_rot% LNZ0ROT
    SALBROT => class_rot% SALBROT
    CSALROT => class_rot% CSALROT
    FSDBROL => class_rot% FSDBROL
    FSFBROL => class_rot% FSFBROL
    FSSBROL => class_rot% FSSBROL
    TBARACC => class_rot% TBARACC
    THALACC => class_rot% THALACC
    THICACC => class_rot% THICACC
    THLQACC => class_rot% THLQACC
    GFLXROW => class_rot% GFLXROW
    HMFGROW => class_rot% HMFGROW
    HTCROW => class_rot% HTCROW
    QFCROW => class_rot% QFCROW
    ITCTROT => class_rot% ITCTROT
    TSFSROT => class_rot% TSFSROT

    ! Point CTEM pointers

      ctem_on           => c_switch%ctem_on
      parallelrun       => c_switch%parallelrun
      cyclemet          => c_switch%cyclemet
      dofire            => c_switch%dofire
      run_model         => c_switch%run_model
      met_rewound       => c_switch%met_rewound
      reach_eof         => c_switch%reach_eof
      compete           => c_switch%compete
      start_bare        => c_switch%start_bare
      rsfile            => c_switch%rsfile
      lnduseon          => c_switch%lnduseon
      co2on             => c_switch%co2on
      ch4on             => c_switch%ch4on
      popdon            => c_switch%popdon
      inibioclim        => c_switch%inibioclim
      dowetlands        => c_switch%dowetlands
      obswetf           => c_switch%obswetf
      transient_run     => c_switch%transient_run
      use_netcdf        => c_switch%use_netcdf
      met_file          => c_switch%met_file
      init_file         => c_switch%init_file

      ! Allocate the ctem structure vrot
      allocate(vrot)

      tcanrs            => vrot%tcanrs
      tsnors            => vrot%tsnors
      tpndrs            => vrot%tpndrs
      csum              => vrot%csum
      tbaraccrow_m      => vrot%tbaraccrow_m
      tcanoaccrow_m     => vrot%tcanoaccrow_m
      uvaccrow_m        => vrot%uvaccrow_m
      vvaccrow_m        => vrot%vvaccrow_m

      ! ROW:
      ailcminrow        => vrot%ailcmin
      ailcmaxrow        => vrot%ailcmax
      dvdfcanrow        => vrot%dvdfcan
      gleafmasrow       => vrot%gleafmas
      bleafmasrow       => vrot%bleafmas
      stemmassrow       => vrot%stemmass
      rootmassrow       => vrot%rootmass
      pstemmassrow      => vrot%pstemmass
      pgleafmassrow     => vrot%pgleafmass
      fcancmxrow        => vrot%fcancmx
      gavglairow        => vrot%gavglai
      zolncrow          => vrot%zolnc
      ailcrow           => vrot%ailc
      ailcgrow          => vrot%ailcg
      ailcgsrow         => vrot%ailcgs
      fcancsrow         => vrot%fcancs
      fcancrow          => vrot%fcanc
      co2concrow        => vrot%co2conc
      ch4concrow        => vrot%ch4conc
      co2i1cgrow        => vrot%co2i1cg
      co2i1csrow        => vrot%co2i1cs
      co2i2cgrow        => vrot%co2i2cg
      co2i2csrow        => vrot%co2i2cs
      ancsvegrow        => vrot%ancsveg
      ancgvegrow        => vrot%ancgveg
      rmlcsvegrow       => vrot%rmlcsveg
      rmlcgvegrow       => vrot%rmlcgveg
      slairow           => vrot%slai
      ailcbrow          => vrot%ailcb
      canresrow         => vrot%canres
      flhrlossrow       => vrot%flhrloss

      tcanoaccrow_out   => vrot%tcanoaccrow_out
      qevpacc_m_save    => vrot%qevpacc_m_save

      grwtheffrow       => vrot%grwtheff
      lystmmasrow       => vrot%lystmmas
      lyrotmasrow       => vrot%lyrotmas
      tymaxlairow       => vrot%tymaxlai
      vgbiomasrow       => vrot%vgbiomas
      gavgltmsrow       => vrot%gavgltms
      gavgscmsrow       => vrot%gavgscms
      stmhrlosrow       => vrot%stmhrlos
      rmatcrow          => vrot%rmatc
      rmatctemrow       => vrot%rmatctem
      litrmassrow       => vrot%litrmass
      soilcmasrow       => vrot%soilcmas
      vgbiomas_vegrow   => vrot%vgbiomas_veg

      emit_co2row       => vrot%emit_co2
      emit_corow        => vrot%emit_co
      emit_ch4row       => vrot%emit_ch4
      emit_nmhcrow      => vrot%emit_nmhc
      emit_h2row        => vrot%emit_h2
      emit_noxrow       => vrot%emit_nox
      emit_n2orow       => vrot%emit_n2o
      emit_pm25row      => vrot%emit_pm25
      emit_tpmrow       => vrot%emit_tpm
      emit_tcrow        => vrot%emit_tc
      emit_ocrow        => vrot%emit_oc
      emit_bcrow        => vrot%emit_bc
      burnfracrow       => vrot%burnfrac
      burnvegfrow       => vrot%burnvegf
      smfuncvegrow      => vrot%smfuncveg
      popdinrow         => vrot%popdin
      btermrow          => vrot%bterm
      ltermrow          => vrot%lterm
      mtermrow          => vrot%mterm

      extnprobrow       => vrot%extnprob
      prbfrhucrow       => vrot%prbfrhuc
      mlightngrow       => vrot%mlightng
      daylrow           => vrot%dayl
      dayl_maxrow       => vrot%dayl_max

      bmasvegrow        => vrot%bmasveg
      cmasvegcrow       => vrot%cmasvegc
      veghghtrow        => vrot%veghght
      rootdpthrow       => vrot%rootdpth
      rmlrow            => vrot%rml
      rmsrow            => vrot%rms
      tltrleafrow       => vrot%tltrleaf
      tltrstemrow       => vrot%tltrstem
      tltrrootrow       => vrot%tltrroot
      leaflitrrow       => vrot%leaflitr
      roottemprow       => vrot%roottemp
      afrleafrow        => vrot%afrleaf
      afrstemrow        => vrot%afrstem
      afrrootrow        => vrot%afrroot
      wtstatusrow       => vrot%wtstatus
      ltstatusrow       => vrot%ltstatus
      rmrrow            => vrot%rmr

      slopefracrow      => vrot%slopefrac
      ch4wet1row        => vrot%ch4wet1
      ch4wet2row        => vrot%ch4wet2
      wetfdynrow        => vrot%wetfdyn
      ch4dyn1row        => vrot%ch4dyn1
      ch4dyn2row        => vrot%ch4dyn2
      wetfrac_monrow    => vrot%wetfrac_mon
      ch4soillsrow      => vrot%ch4_soills

      lucemcomrow       => vrot%lucemcom
      lucltrinrow       => vrot%lucltrin
      lucsocinrow       => vrot%lucsocin

      npprow            => vrot%npp
      neprow            => vrot%nep
      nbprow            => vrot%nbp
      gpprow            => vrot%gpp
      hetroresrow       => vrot%hetrores
      autoresrow        => vrot%autores
      soilcresprow      => vrot%soilcresp
      rmrow             => vrot%rm
      rgrow             => vrot%rg
      litresrow         => vrot%litres
      socresrow         => vrot%socres
      dstcemlsrow       => vrot%dstcemls
      litrfallrow       => vrot%litrfall
      humiftrsrow       => vrot%humiftrs

      gppvegrow         => vrot%gppveg
      nepvegrow         => vrot%nepveg
      nbpvegrow         => vrot%nbpveg
      nppvegrow         => vrot%nppveg
      hetroresvegrow    => vrot%hetroresveg
      autoresvegrow     => vrot%autoresveg
      litresvegrow      => vrot%litresveg
      soilcresvegrow    => vrot%soilcresveg
      rmlvegaccrow      => vrot%rmlvegacc
      rmsvegrow         => vrot%rmsveg
      rmrvegrow         => vrot%rmrveg
      rgvegrow          => vrot%rgveg
      litrfallvegrow    => vrot%litrfallveg
      humiftrsvegrow    => vrot%humiftrsveg

      rothrlosrow       => vrot%rothrlos
      pfcancmxrow       => vrot%pfcancmx
      nfcancmxrow       => vrot%nfcancmx
      alvsctmrow        => vrot%alvsctm
      paicrow           => vrot%paic
      slaicrow          => vrot%slaic
      alirctmrow        => vrot%alirctm
      cfluxcgrow        => vrot%cfluxcg
      cfluxcsrow        => vrot%cfluxcs
      dstcemls3row      => vrot%dstcemls3
      anvegrow          => vrot%anveg
      rmlvegrow         => vrot%rmlveg

      pftexistrow       => vrot%pftexist
      colddaysrow       => vrot%colddays
      icountrow         => vrot%icount
      lfstatusrow       => vrot%lfstatus
      pandaysrow        => vrot%pandays
      stdalnrow         => vrot%stdaln

      twarmmrow            => vrot%twarmm
      tcoldmrow            => vrot%tcoldm
      gdd5row              => vrot%gdd5
      aridityrow           => vrot%aridity
      srplsmonrow          => vrot%srplsmon
      defctmonrow          => vrot%defctmon
      anndefctrow          => vrot%anndefct
      annsrplsrow          => vrot%annsrpls
      annpcprow            => vrot%annpcp
      dry_season_lengthrow => vrot%dry_season_length

      altotcntr_d       => vrot%altotcntr_d

      ! >>>>>>>>>>>>>>>>>>>>>>>>>>
      ! GAT:

      ! Allocate the ctem structure vgat
      allocate(vgat)

      lightng           => vgat%lightng
      tcanoaccgat_out   => vgat%tcanoaccgat_out

      ailcmingat        => vgat%ailcmin
      ailcmaxgat        => vgat%ailcmax
      dvdfcangat        => vgat%dvdfcan
      gleafmasgat       => vgat%gleafmas
      bleafmasgat       => vgat%bleafmas
      stemmassgat       => vgat%stemmass
      rootmassgat       => vgat%rootmass
      pstemmassgat      => vgat%pstemmass
      pgleafmassgat     => vgat%pgleafmass
      fcancmxgat        => vgat%fcancmx
      gavglaigat        => vgat%gavglai
      zolncgat          => vgat%zolnc
      ailcgat           => vgat%ailc
      ailcggat          => vgat%ailcg
      ailcgsgat         => vgat%ailcgs
      fcancsgat         => vgat%fcancs
      fcancgat          => vgat%fcanc
      co2concgat        => vgat%co2conc
      ch4concgat        => vgat%ch4conc
      co2i1cggat        => vgat%co2i1cg
      co2i1csgat        => vgat%co2i1cs
      co2i2cggat        => vgat%co2i2cg
      co2i2csgat        => vgat%co2i2cs
      ancsveggat        => vgat%ancsveg
      ancgveggat        => vgat%ancgveg
      rmlcsveggat       => vgat%rmlcsveg
      rmlcgveggat       => vgat%rmlcgveg
      slaigat           => vgat%slai
      ailcbgat          => vgat%ailcb
      canresgat         => vgat%canres
      flhrlossgat       => vgat%flhrloss

      grwtheffgat       => vgat%grwtheff
      lystmmasgat       => vgat%lystmmas
      lyrotmasgat       => vgat%lyrotmas
      tymaxlaigat       => vgat%tymaxlai
      vgbiomasgat       => vgat%vgbiomas
      gavgltmsgat       => vgat%gavgltms
      gavgscmsgat       => vgat%gavgscms
      stmhrlosgat       => vgat%stmhrlos
      rmatcgat          => vgat%rmatc
      rmatctemgat       => vgat%rmatctem
      litrmassgat       => vgat%litrmass
      soilcmasgat       => vgat%soilcmas
      vgbiomas_veggat   => vgat%vgbiomas_veg
      litrfallveggat    => vgat%litrfallveg
      humiftrsveggat    => vgat%humiftrsveg

      emit_co2gat       => vgat%emit_co2
      emit_cogat        => vgat%emit_co
      emit_ch4gat       => vgat%emit_ch4
      emit_nmhcgat      => vgat%emit_nmhc
      emit_h2gat        => vgat%emit_h2
      emit_noxgat       => vgat%emit_nox
      emit_n2ogat       => vgat%emit_n2o
      emit_pm25gat      => vgat%emit_pm25
      emit_tpmgat       => vgat%emit_tpm
      emit_tcgat        => vgat%emit_tc
      emit_ocgat        => vgat%emit_oc
      emit_bcgat        => vgat%emit_bc
      burnfracgat       => vgat%burnfrac
      burnvegfgat       => vgat%burnvegf
      popdingat         => vgat%popdin
      smfuncveggat      => vgat%smfuncveg
      btermgat          => vgat%bterm
      ltermgat          => vgat%lterm
      mtermgat          => vgat%mterm

      extnprobgat       => vgat%extnprob
      prbfrhucgat       => vgat%prbfrhuc
      mlightnggat       => vgat%mlightng
      daylgat           => vgat%dayl
      dayl_maxgat       => vgat%dayl_max

      bmasveggat        => vgat%bmasveg
      cmasvegcgat       => vgat%cmasvegc
      veghghtgat        => vgat%veghght
      rootdpthgat       => vgat%rootdpth
      rmlgat            => vgat%rml
      rmsgat            => vgat%rms
      tltrleafgat       => vgat%tltrleaf
      tltrstemgat       => vgat%tltrstem
      tltrrootgat       => vgat%tltrroot
      leaflitrgat       => vgat%leaflitr
      roottempgat       => vgat%roottemp
      afrleafgat        => vgat%afrleaf
      afrstemgat        => vgat%afrstem
      afrrootgat        => vgat%afrroot
      wtstatusgat       => vgat%wtstatus
      ltstatusgat       => vgat%ltstatus
      rmrgat            => vgat%rmr

      slopefracgat      => vgat%slopefrac
      wetfrac_presgat   => vgat%wetfrac_pres
      wetfrac_mongat    => vgat%wetfrac_mon
      ch4wet1gat        => vgat%ch4wet1
      ch4wet2gat        => vgat%ch4wet2
      wetfdyngat        => vgat%wetfdyn
      ch4dyn1gat        => vgat%ch4dyn1
      ch4dyn2gat        => vgat%ch4dyn2
      ch4soillsgat      => vgat%ch4_soills

      lucemcomgat       => vgat%lucemcom
      lucltringat       => vgat%lucltrin
      lucsocingat       => vgat%lucsocin

      nppgat            => vgat%npp
      nepgat            => vgat%nep
      nbpgat            => vgat%nbp
      gppgat            => vgat%gpp
      hetroresgat       => vgat%hetrores
      autoresgat        => vgat%autores
      soilcrespgat      => vgat%soilcresp
      rmgat             => vgat%rm
      rggat             => vgat%rg
      litresgat         => vgat%litres
      socresgat         => vgat%socres
      dstcemlsgat       => vgat%dstcemls
      litrfallgat       => vgat%litrfall
      humiftrsgat       => vgat%humiftrs

      gppveggat         => vgat%gppveg
      nepveggat         => vgat%nepveg
      nbpveggat         => vgat%nbpveg
      nppveggat         => vgat%nppveg
      hetroresveggat    => vgat%hetroresveg
      autoresveggat     => vgat%autoresveg
      litresveggat      => vgat%litresveg
      soilcresveggat    => vgat%soilcresveg
      rmlvegaccgat      => vgat%rmlvegacc
      rmsveggat         => vgat%rmsveg
      rmrveggat         => vgat%rmrveg
      rgveggat          => vgat%rgveg

      rothrlosgat       => vgat%rothrlos
      pfcancmxgat       => vgat%pfcancmx
      nfcancmxgat       => vgat%nfcancmx
      alvsctmgat        => vgat%alvsctm
      paicgat           => vgat%paic
      slaicgat          => vgat%slaic
      alirctmgat        => vgat%alirctm
      cfluxcggat        => vgat%cfluxcg
      cfluxcsgat        => vgat%cfluxcs
      dstcemls3gat      => vgat%dstcemls3
      anveggat          => vgat%anveg
      rmlveggat         => vgat%rmlveg

      twarmmgat            => vgat%twarmm
      tcoldmgat            => vgat%tcoldm
      gdd5gat              => vgat%gdd5
      ariditygat           => vgat%aridity
      srplsmongat          => vgat%srplsmon
      defctmongat          => vgat%defctmon
      anndefctgat          => vgat%anndefct
      annsrplsgat          => vgat%annsrpls
      annpcpgat            => vgat%annpcp
      dry_season_lengthgat => vgat%dry_season_length

      tcurm             => vgat%tcurm
      srpcuryr          => vgat%srpcuryr
      dftcuryr          => vgat%dftcuryr
      tmonth            => vgat%tmonth
      anpcpcur          => vgat%anpcpcur
      anpecur           => vgat%anpecur
      gdd5cur           => vgat%gdd5cur
      surmncur          => vgat%surmncur
      defmncur          => vgat%defmncur
      srplscur          => vgat%srplscur
      defctcur          => vgat%defctcur

      geremortgat       => vgat%geremort
      intrmortgat       => vgat%intrmort
      lambdagat         => vgat%lambda
      ccgat             => vgat%cc
      mmgat             => vgat%mm

      pftexistgat       => vgat%pftexist
      colddaysgat       => vgat%colddays
      icountgat         => vgat%icount
      lfstatusgat       => vgat%lfstatus
      pandaysgat        => vgat%pandays
      stdalngat         => vgat%stdaln

      ! Mosaic-level (CLASS vars):

      PREACC_M          => vrot%PREACC_M
      GTACC_M           => vrot%GTACC_M
      QEVPACC_M         => vrot%QEVPACC_M
      HFSACC_M          => vrot%HFSACC_M
      HMFNACC_M         => vrot%HMFNACC_M
      ROFACC_M          => vrot%ROFACC_M
      SNOACC_M          => vrot%SNOACC_M
      OVRACC_M          => vrot%OVRACC_M
      WTBLACC_M         => vrot%WTBLACC_M
      TBARACC_M         => vrot%TBARACC_M
      THLQACC_M         => vrot%THLQACC_M
      THICACC_M         => vrot%THICACC_M
      THALACC_M         => vrot%THALACC_M
      ALVSACC_M         => vrot%ALVSACC_M
      ALIRACC_M         => vrot%ALIRACC_M
      RHOSACC_M         => vrot%RHOSACC_M
      TSNOACC_M         => vrot%TSNOACC_M
      WSNOACC_M         => vrot%WSNOACC_M
      SNOARE_M          => vrot%SNOARE_M
      TCANACC_M         => vrot%TCANACC_M
      RCANACC_M         => vrot%RCANACC_M
      SCANACC_M         => vrot%SCANACC_M
      GROACC_M          => vrot%GROACC_M
      FSINACC_M         => vrot%FSINACC_M
      FLINACC_M         => vrot%FLINACC_M
      TAACC_M           => vrot%TAACC_M
      UVACC_M           => vrot%UVACC_M
      PRESACC_M         => vrot%PRESACC_M
      QAACC_M           => vrot%QAACC_M
      ALTOTACC_M        => vrot%ALTOTACC_M
      EVAPACC_M         => vrot%EVAPACC_M
      FLUTACC_M         => vrot%FLUTACC_M

      ! grid-averaged (CLASS vars)

      ! Allocate the ctem structure ctem_grd
      allocate(ctem_grd)

      WSNOROT_g         => ctem_grd%WSNOROT_g
      ROFSROT_g         => ctem_grd%ROFSROT_g
      SNOROT_g          => ctem_grd%SNOROT_g
      RHOSROT_g         => ctem_grd%RHOSROT_g
      ROFROT_g          => ctem_grd%ROFROT_g
      ZPNDROT_g         => ctem_grd%ZPNDROT_g
      RCANROT_g         => ctem_grd%RCANROT_g
      SCANROT_g         => ctem_grd%SCANROT_g
      TROFROT_g         => ctem_grd%TROFROT_g
      TROOROT_g         => ctem_grd%TROOROT_g
      TROBROT_g         => ctem_grd%TROBROT_g
      ROFOROT_g         => ctem_grd%ROFOROT_g
      ROFBROT_g         => ctem_grd%ROFBROT_g
      TROSROT_g         => ctem_grd%TROSROT_g
      FSGVROT_g         => ctem_grd%FSGVROT_g
      FSGSROT_g         => ctem_grd%FSGSROT_g
      FLGVROT_g         => ctem_grd%FLGVROT_g
      FLGSROT_g         => ctem_grd%FLGSROT_g
      HFSCROT_g         => ctem_grd%HFSCROT_g
      HFSSROT_g         => ctem_grd%HFSSROT_g
      HEVCROT_g         => ctem_grd%HEVCROT_g
      HEVSROT_g         => ctem_grd%HEVSROT_g
      HMFCROT_g         => ctem_grd%HMFCROT_g
      HMFNROT_g         => ctem_grd%HMFNROT_g
      HTCSROT_g         => ctem_grd%HTCSROT_g
      HTCCROT_g         => ctem_grd%HTCCROT_g
      FSGGROT_g         => ctem_grd%FSGGROT_g
      FLGGROT_g         => ctem_grd%FLGGROT_g
      HFSGROT_g         => ctem_grd%HFSGROT_g
      HEVGROT_g         => ctem_grd%HEVGROT_g
      CDHROT_g          => ctem_grd%CDHROT_g
      CDMROT_g          => ctem_grd%CDMROT_g
      SFCUROT_g         => ctem_grd%SFCUROT_g
      SFCVROT_g         => ctem_grd%SFCVROT_g
      ACTLYR_g          => ctem_grd%ACTLYR_g
      FTABLE_g          => ctem_grd%FTABLE_g
      fc_g              => ctem_grd%fc_g
      fg_g              => ctem_grd%fg_g
      fcs_g             => ctem_grd%fcs_g
      fgs_g             => ctem_grd%fgs_g
      PCFCROT_g         => ctem_grd%PCFCROT_g
      PCLCROT_g         => ctem_grd%PCLCROT_g
      PCPGROT_g         => ctem_grd%PCPGROT_g
      QFCFROT_g         => ctem_grd%QFCFROT_g
      QFGROT_g          => ctem_grd%QFGROT_g
      QFCROT_g          => ctem_grd%QFCROT_g
      ROFCROT_g         => ctem_grd%ROFCROT_g
      ROFNROT_g         => ctem_grd%ROFNROT_g
      WTRSROT_g         => ctem_grd%WTRSROT_g
      WTRGROT_g         => ctem_grd%WTRGROT_g
      PCPNROT_g         => ctem_grd%PCPNROT_g
      QFCLROT_g         => ctem_grd%QFCLROT_g
      QFNROT_g          => ctem_grd%QFNROT_g
      WTRCROT_g         => ctem_grd%WTRCROT_g
      rmlvegrow_g       => ctem_grd%rmlvegrow_g
      anvegrow_g        => ctem_grd%anvegrow_g
      HMFGROT_g         => ctem_grd%HMFGROT_g
      HTCROT_g          => ctem_grd%HTCROT_g
      TBARROT_g         => ctem_grd%TBARROT_g
      THLQROT_g         => ctem_grd%THLQROT_g
      THICROT_g         => ctem_grd%THICROT_g
      GFLXROT_g         => ctem_grd%GFLXROT_g

       fsstar_g         => ctem_grd%fsstar_g
       flstar_g         => ctem_grd%flstar_g
       qh_g             => ctem_grd%qh_g
       qe_g             => ctem_grd%qe_g
       snomlt_g         => ctem_grd%snomlt_g
       beg_g            => ctem_grd%beg_g
       gtout_g          => ctem_grd%gtout_g
       tpn_g            => ctem_grd%tpn_g
       altot_g          => ctem_grd%altot_g
       tcn_g            => ctem_grd%tcn_g
       tsn_g            => ctem_grd%tsn_g
       zsn_g            => ctem_grd%zsn_g

      ! mosaic level variables (CLASS):

      ! Allocate the ctem structure ctem_tile
      allocate(ctem_tile)

      fsnowacc_t        => ctem_tile%fsnowacc_t
      tcansacc_t        => ctem_tile%tcansacc_t
      tcanoaccgat_t     => ctem_tile%tcanoaccgat_t
      taaccgat_t        => ctem_tile%taaccgat_t
      uvaccgat_t        => ctem_tile%uvaccgat_t
      vvaccgat_t        => ctem_tile%vvaccgat_t
      tbaraccgat_t      => ctem_tile%tbaraccgat_t
      tbarcacc_t        => ctem_tile%tbarcacc_t
      tbarcsacc_t       => ctem_tile%tbarcsacc_t
      tbargacc_t        => ctem_tile%tbargacc_t
      tbargsacc_t       => ctem_tile%tbargsacc_t
      thliqcacc_t       => ctem_tile%thliqcacc_t
      thliqgacc_t       => ctem_tile%thliqgacc_t
      thliqacc_t        => ctem_tile%thliqacc_t
      thicecacc_t       => ctem_tile%thicecacc_t
      thicegacc_t       => ctem_tile%thicegacc_t
      ancsvgac_t        => ctem_tile%ancsvgac_t
      ancgvgac_t        => ctem_tile%ancgvgac_t
      rmlcsvga_t        => ctem_tile%rmlcsvga_t
      rmlcgvga_t        => ctem_tile%rmlcgvga_t

!    =================================================================================
!    =================================================================================

!    Allocate some other data structures:
    allocate(class_out,ctem_mo,ctem_grd_mo,ctem_tile_mo, &
      &      ctem_yr,ctem_grd_yr,ctem_tile_yr)

!    Declarations are complete, run preparations begin
      CALL CLASSD

!     Initialize the CTEM parameters
      call initpftpars(compete)

!> Since we know the nlat, nmos, ignd, and ilg we can allocate the CLASS and
!! CTEM variable structures.
      call alloc_class_vars()
      call alloc_ctem_vars()

      ! Allocate the local variables that rely on nlat, ilg, etc.
      allocate(curlatno(ilg),&
               altotcount_ctm(nlat),&
               todfrac(ilg,icc),&
               barf(nlat,nmos),&
               currlat(ilg),&
               wl(lat),&
               grclarea(ilg),&
               wossl(lat),&
               sl(lat),&
               radl(lat),&
               cl(lat),&
               ml(ilg),&
               fsinacc_gat(ilg),&
               flutacc_gat(ilg),&
               flinacc_gat(ilg),&
               alswacc_gat(ilg),&
               allwacc_gat(ilg),&
               pregacc_gat(ilg),&
               altotacc_gat(ilg),&
               netrad_gat(ilg),&
               preacc_gat(ilg),&
               sdepgat(ilg),&
               rgmgat(ilg,ignd),&
               sandgat(ilg,ignd),&
               claygat(ilg,ignd),&
               orgmgat(ilg,ignd),&
               xdiffusgat(ilg),& ! the corresponding ROW is CLASS's XDIFFUS
               faregat(ilg),&    ! the ROT is FAREROT
               FTABLE(nlat,nmos),&
               ACTLYR(nlat,nmos))

!>    This reads in from the restart file (replacing the INI and CTM
!!    files).
      call read_initialstate()
      write(*,*)zrfmrow
        print *,'after readinitialstate'

! #ED - ignore after this.
!>     Open the met netcdf file. This also sets up the run boundaries
!!     based on the metadata in the netcdf. It is important to ensure the
!!     netcdf is of the same dimensions as the intended output files.
!!     Based upon the bounds used to call the model, this will figure out how
!!     big the NLAT vector is.
!      call openmet()
!      write(*,*)'done openmet'

      ! #ED - later on, we will read in the MET forcing data using netcdf like this
      ! at present we have to still rely on the ASCII text files so this is commented
      ! out (and also not really coded up).

      !call readin_met(1,dlatgat,dlongat)
!
!     checking the time spent for running model
!
!      call idate(today)
!      call itime(now)
!      write(*,1000)   today(2), today(1), 2000+today(3), now
! 1000 format( 'start date: ', i2.2, '/', i2.2, '/', i4.4,
!     &      '; start time: ', i2.2, ':', i2.2, ':', i2.2 )
!
!     INITIALIZATION FOR COUPLING CLASS AND CTEM
!
       call initrowvars()
       call resetclassaccum(nltest,nmtest)

      end subroutine CLASSIC_driver

      end module main_driver

