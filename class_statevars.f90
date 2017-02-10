!>\defgroup class_statevars

!>
!!this module contains the variable type structures:
!! 1. class_rot - CLASS's 'rot' and 'row' vars
!! 2. class_gat - CLASS's 'gat' vars
!
!>\file
module class_statevars

! J. Melton Nov 2016

use ctem_params,  only : ican, icp1, NBS

implicit none

public :: alloc_class_vars

!=================================================================================
!>
type class_gather

! These will be allocated the dimension: 'ilg'

    integer, allocatable, dimension(:) :: ILMOS     !<Index of grid cell corresponding to current element of gathered vector of land surface variables [ ]
    integer, allocatable, dimension(:) :: JLMOS     !<Index of mosaic tile corresponding to current element of gathered vector of land surface variables [ ]
    integer, allocatable, dimension(:) :: IWMOS     !<Index of grid cell corresponding to current element of gathered vector of inland water body variables [ ]
    integer, allocatable, dimension(:) :: JWMOS     !<Index of mosaic tile corresponding to current element of gathered vector of inland water body variables [ ]
    integer, allocatable, dimension(:) :: IGDRGAT   !<Index of soil layer in which bedrock is encountered

    real, allocatable, dimension(:) :: DELZ    !<
    real, allocatable, dimension(:) :: ZBOT    !<
    real, allocatable, dimension(:) :: ALBSGAT !<Snow albedo [ ]
    real, allocatable, dimension(:) :: CMAIGAT !<Aggregated mass of vegetation canopy \f$[kg m^{-2} ]\f$
    real, allocatable, dimension(:) :: GROGAT  !<Vegetation growth index [ ]
    real, allocatable, dimension(:) :: QACGAT  !<Specific humidity of air within vegetation canopy space \f$[kg kg^{-1} ]\f$
    real, allocatable, dimension(:) :: RCANGAT !<Intercepted liquid water stored on canopy \f$[kg m^{-2} ]\f$
    real, allocatable, dimension(:) :: RHOSGAT !<Density of snow \f$[kg m^{-3} ]\f$
    real, allocatable, dimension(:) :: SCANGAT !<Intercepted frozen water stored on canopy \f$[kg m^{-2} ]\f$
    real, allocatable, dimension(:) :: SNOGAT  !<Mass of snow pack [kg m^{-2} ]\f$
    real, allocatable, dimension(:) :: TACGAT  !<Temperature of air within vegetation canopy [K]
    real, allocatable, dimension(:) :: TBASGAT !<Temperature of bedrock in third soil layer [K]
    real, allocatable, dimension(:) :: TCANGAT !<Vegetation canopy temperature [K]
    real, allocatable, dimension(:) :: TPNDGAT !<Temperature of ponded water [K]
    real, allocatable, dimension(:) :: TSNOGAT !<Snowpack temperature [K]
    real, allocatable, dimension(:) :: WSNOGAT !<Liquid water content of snow pack \f$[kg m^{-2} ]\f$
    real, allocatable, dimension(:) :: ZPNDGAT !<Depth of ponded water on surface [m]
    real, allocatable, dimension(:) :: REFGAT  !<
    real, allocatable, dimension(:) :: BCSNGAT !<
    real, allocatable, dimension(:) :: AGIDGAT !<Optional user-specified value of ground near-infrared albedo to override CLASS-calculated value [ ]
    real, allocatable, dimension(:) :: AGVDGAT !<Optional user-specified value of ground visible albedo to override CLASS-calculated value [ ]
    real, allocatable, dimension(:) :: ALGDGAT !<Reference albedo for dry soil [ ]
    real, allocatable, dimension(:) :: ALGWGAT !<Reference albedo for saturated soil [ ]
    real, allocatable, dimension(:) :: ASIDGAT !<Optional user-specified value of snow near-infrared albedo to override CLASS-calculated value [ ]
    real, allocatable, dimension(:) :: ASVDGAT !<Optional user-specified value of snow visible albedo to override CLASS-calculated value [ ]
    real, allocatable, dimension(:) :: DRNGAT  !<Drainage index at bottom of soil profile [ ]
    real, allocatable, dimension(:) :: GRKFGAT !<WATROF parameter used when running MESH code [ ]
    real, allocatable, dimension(:) :: WFCIGAT !<WATROF parameter used when running MESH code [ ]
    real, allocatable, dimension(:) :: WFSFGAT !<WATROF parameter used when running MESH code [ ]
    real, allocatable, dimension(:) :: XSLPGAT !<Surface slope (used when running MESH code) [degrees]
    real, allocatable, dimension(:) :: ZPLGGAT !<Maximum water ponding depth for snow-free subareas (user-specified when running MESH code) [m]
    real, allocatable, dimension(:) :: ZPLSGAT !<Maximum water ponding depth for snow-covered subareas (user-specified when running MESH code) [m]
    real, allocatable, dimension(:) :: ZSNLGAT !<Limiting snow depth below which coverage is < 100% [m]
    real, allocatable, dimension(:) :: ALGWVGAT !<
    real, allocatable, dimension(:) :: ALGWNGAT !<
    real, allocatable, dimension(:) :: ALGDVGAT !<
    real, allocatable, dimension(:) :: ALGDNGAT !<
    real, allocatable, dimension(:) :: EMISGAT  !<
    real, allocatable, dimension(:) :: CSZGAT  !<Cosine of solar zenith angle [ ]
    real, allocatable, dimension(:) :: DLONGAT !<Longitude of grid cell (east of Greenwich) [degrees]
    real, allocatable, dimension(:) :: DLATGAT !< Latitude of grid cell [degrees]
    real, allocatable, dimension(:) :: FCLOGAT !<Fractional cloud cover [ ]
    real, allocatable, dimension(:) :: FDLGAT  !<Downwelling longwave radiation at bottom of atmosphere (i.e. incident on modelled land surface elements \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: FSIHGAT !<Near-infrared radiation incident on horizontal surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: FSVHGAT !<Visible radiation incident on horizontal surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: GGEOGAT !<Geothermal heat flux at bottom of soil profile \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: PADRGAT !<Partial pressure of dry air [Pa]
    real, allocatable, dimension(:) :: PREGAT  !<Surface precipitation rate \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: PRESGAT !<Surface air pressure [Pa]
    real, allocatable, dimension(:) :: QAGAT   !<Specific humidity at reference height \f$[kg kg^{-1} ]\f$
    real, allocatable, dimension(:) :: RADJGAT !<Latitude of grid cell (positive north of equator) [rad]
    real, allocatable, dimension(:) :: RHOAGAT !<Density of air \f$[kg m^{-3} ]\f$
    real, allocatable, dimension(:) :: RHSIGAT !<Density of fresh snow \f$[kg m^{-3} ]\f$
    real, allocatable, dimension(:) :: RPCPGAT !<Rainfall rate over modelled area \f$[m s^{-1} ]\f$
    real, allocatable, dimension(:) :: SPCPGAT !<Snowfall rate over modelled area \f$[m s^{-1} ]\f$
    real, allocatable, dimension(:) :: TAGAT   !<Air temperature at reference height [K]
    real, allocatable, dimension(:) :: TADPGAT !<Dew point temperature of air [K]
    real, allocatable, dimension(:) :: TRPCGAT !<Rainfall temperature [K]
    real, allocatable, dimension(:) :: TSPCGAT !<Snowfall temperature [K]
    real, allocatable, dimension(:) :: ULGAT   !<Zonal component of wind velocity \f$[m s^{-1} ]\f$
    real, allocatable, dimension(:) :: VLGAT   !<Meridional component of wind velocity \f$[m s^{-1} ]\f$
    real, allocatable, dimension(:) :: VMODGAT !<Wind speed at reference height \f$[m s^{-1} ]\f$
    real, allocatable, dimension(:) :: VPDGAT  !<Vapour pressure deficit [mb]
    real, allocatable, dimension(:) :: Z0ORGAT !<Orographic roughness length [m]
    real, allocatable, dimension(:) :: ZBLDGAT !<Atmospheric blending height for surface roughness length averaging [m]
    real, allocatable, dimension(:) :: ZDHGAT  !<User-specified height associated with diagnosed screen-level variables [m]
    real, allocatable, dimension(:) :: ZDMGAT  !<User-specified height associated with diagnosed anemometer-level wind speed [m]
    real, allocatable, dimension(:) :: ZRFHGAT !<Reference height associated with forcing air temperature and humidity [m]
    real, allocatable, dimension(:) :: ZRFMGAT !<Reference height associated with forcing wind speed [m]
    real, allocatable, dimension(:) :: FSGGAT  !<
    real, allocatable, dimension(:) :: FLGGAT  !<
    real, allocatable, dimension(:) :: GUSTGAT !<
    real, allocatable, dimension(:) :: DEPBGAT !<
    real, allocatable, dimension(:) :: GTBS    !<
    real, allocatable, dimension(:) :: SFCUBS  !<
    real, allocatable, dimension(:) :: SFCVBS  !<
    real, allocatable, dimension(:) :: USTARBS !<
    real, allocatable, dimension(:) :: TCSNOW  !<
    real, allocatable, dimension(:) :: GSNOW   !<
    real, allocatable, dimension(:) :: ALIRGAT !<Diagnosed total near-infrared albedo of land surface [ ]
    real, allocatable, dimension(:) :: ALVSGAT !<Diagnosed total visible albedo of land surface [ ]
    real, allocatable, dimension(:) :: CDHGAT  !<Surface drag coefficient for heat [ ]
    real, allocatable, dimension(:) :: CDMGAT  !<Surface drag coefficient for momentum [ ]
    real, allocatable, dimension(:) :: DRGAT   !<Surface drag coefficient under neutral stability [ ]
    real, allocatable, dimension(:) :: EFGAT   !<Evaporation efficiency at ground surface [ ]
    real, allocatable, dimension(:) :: FLGGGAT !<Diagnosed net longwave radiation at soil surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: FLGSGAT !<Diagnosed net longwave radiation at snow surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: FLGVGAT !<Diagnosed net longwave radiation on vegetation canopy \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: FSGGGAT !<Diagnosed net shortwave radiation at soil surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: FSGSGAT !<Diagnosed net shortwave radiation at snow surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: FSGVGAT !<Diagnosed net shortwave radiation on vegetation canopy \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: FSNOGAT !<Diagnosed fractional snow coverage [ ]
    real, allocatable, dimension(:) :: GAGAT   !<Diagnosed product of drag coefficient and wind speed over modelled area \f$[m s^{-1} ]\f$
    real, allocatable, dimension(:) :: GTGAT   !<Diagnosed effective surface black-body temperature [K]
    real, allocatable, dimension(:) :: HBLGAT  !<Height of the atmospheric boundary layer [m]
    real, allocatable, dimension(:) :: HEVCGAT !<Diagnosed latent heat flux on vegetation canopy \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: HEVGGAT !<Diagnosed latent heat flux at soil surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: HEVSGAT !<Diagnosed latent heat flux at snow surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: HFSGAT  !<Diagnosed total surface sensible heat flux over modelled area \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: HFSCGAT !<Diagnosed sensible heat flux on vegetation canopy \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: HFSGGAT !<Diagnosed sensible heat flux at soil surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: HFSSGAT !<Diagnosed sensible heat flux at snow surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: HMFCGAT !<Diagnosed energy associated with phase change of water on vegetation \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: HMFNGAT !<Diagnosed energy associated with phase change of water in snow pack \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: HTCCGAT !<Diagnosed internal energy change of vegetation canopy due to conduction and/or change in mass \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: HTCSGAT !<Diagnosed internal energy change of snow pack due to conduction and/or change in mass \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: ILMOGAT !<Inverse of Monin-Obukhov roughness length \f$(m^{-1} ]\f$
    real, allocatable, dimension(:) :: PCFCGAT !<Diagnosed frozen precipitation intercepted by vegetation \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: PCLCGAT !<Diagnosed liquid precipitation intercepted by vegetation \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: PCPGGAT !<Diagnosed precipitation incident on ground \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: PCPNGAT !<Diagnosed precipitation incident on snow pack \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: PETGAT  !<Diagnosed potential evapotranspiration \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: QEVPGAT !<Diagnosed total surface latent heat flux over modelled area \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: QFCFGAT !<Diagnosed vapour flux from frozen water on vegetation \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: QFCLGAT !<Diagnosed vapour flux from liquid water on vegetation \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: QFGGAT  !<Diagnosed water vapour flux from ground \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: QFNGAT  !<Diagnosed water vapour flux from snow pack \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: QFSGAT  !<Diagnosed total surface water vapour flux over modelled area \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: QFXGAT  !<Product of surface drag coefficient, wind speed and surface-air specific humidity difference \f$[m s^{-1} ]\f$
    real, allocatable, dimension(:) :: QGGAT   !<Diagnosed surface specific humidity \f$[kg kg^{-1} ]\f$
    real, allocatable, dimension(:) :: ROFGAT  !<Total runoff from soil \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: ROFBGAT !<Base flow from bottom of soil column \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: ROFCGAT !<Liquid/frozen water runoff from vegetation \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: ROFNGAT !<Liquid water runoff from snow pack \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: ROFOGAT !<Overland flow from top of soil column \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: ROFSGAT !<Interflow from sides of soil column \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: ROVGGAT !<Diagnosed liquid/frozen water runoff from vegetation to ground surface \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: SFCQGAT !<Diagnosed screen-level specific humidity \f$[kg kg^{-1} ]\f$
    real, allocatable, dimension(:) :: SFCTGAT !<Diagnosed screen-level air temperature [K]
    real, allocatable, dimension(:) :: SFCUGAT !<Diagnosed anemometer-level zonal wind \f$[m s^{-1} ]\f$
    real, allocatable, dimension(:) :: SFCVGAT !<Diagnosed anemometer-level meridional wind \f$[m s^{-1} ]\f$
    real, allocatable, dimension(:) :: TFXGAT  !<Product of surface drag coefficient, wind speed and surface-air temperature difference \f$[K m s^{-1} ]\f$
    real, allocatable, dimension(:) :: TROBGAT !<Temperature of base flow from bottom of soil column [K]
    real, allocatable, dimension(:) :: TROFGAT !<Temperature of total runoff [K]
    real, allocatable, dimension(:) :: TROOGAT !<Temperature of overland flow from top of soil column [K]
    real, allocatable, dimension(:) :: TROSGAT !<Temperature of interflow from sides of soil column [K]
    real, allocatable, dimension(:) :: UEGAT   !<Friction velocity of air \f$[m s^{-1} ]\f$
    real, allocatable, dimension(:) :: WTABGAT !<Depth of water table in soil [m]
    real, allocatable, dimension(:) :: WTRCGAT !<Diagnosed residual water transferred off the vegetation canopy \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: WTRGGAT !<Diagnosed residual water transferred into or out of the soil \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: WTRSGAT !<Diagnosed residual water transferred into or out of the snow pack \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: QLWOGAT !<
    real, allocatable, dimension(:) :: SFRHGAT !<
    real, allocatable, dimension(:) :: FTEMP   !<
    real, allocatable, dimension(:) :: FVAP    !<
    real, allocatable, dimension(:) :: RIB     !<
    real, allocatable, dimension(:) :: FC      !<
    real, allocatable, dimension(:) :: FG      !<
    real, allocatable, dimension(:) :: FCS     !<
    real, allocatable, dimension(:) :: FGS     !<
    real, allocatable, dimension(:) :: RBCOEF  !<
    real, allocatable, dimension(:) :: ZSNOW   !<
    real, allocatable, dimension(:) :: FSVF    !<
    real, allocatable, dimension(:) :: FSVFS   !<
    real, allocatable, dimension(:) :: ALVSCN  !<
    real, allocatable, dimension(:) :: ALIRCN  !<
    real, allocatable, dimension(:) :: ALVSG   !<
    real, allocatable, dimension(:) :: ALIRG   !<
    real, allocatable, dimension(:) :: ALVSCS  !<
    real, allocatable, dimension(:) :: ALIRCS  !<
    real, allocatable, dimension(:) :: ALVSSN  !<
    real, allocatable, dimension(:) :: ALIRSN  !<
    real, allocatable, dimension(:) :: ALVSGC  !<
    real, allocatable, dimension(:) :: ALIRGC  !<
    real, allocatable, dimension(:) :: ALVSSC  !<
    real, allocatable, dimension(:) :: ALIRSC  !<
    real, allocatable, dimension(:) :: TRVSCN  !<
    real, allocatable, dimension(:) :: TRIRCN  !<
    real, allocatable, dimension(:) :: TRVSCS  !<
    real, allocatable, dimension(:) :: TRIRCS  !<
    real, allocatable, dimension(:) :: RC      !<
    real, allocatable, dimension(:) :: RCS     !<
    real, allocatable, dimension(:) :: FRAINC  !<
    real, allocatable, dimension(:) :: FSNOWC  !<
    real, allocatable, dimension(:) :: FRAICS  !<
    real, allocatable, dimension(:) :: FSNOCS  !<
    real, allocatable, dimension(:) :: CMASSC  !<
    real, allocatable, dimension(:) :: CMASCS  !<
    real, allocatable, dimension(:) :: DISP    !<
    real, allocatable, dimension(:) :: DISPS   !<
    real, allocatable, dimension(:) :: ZOMLNC  !<
    real, allocatable, dimension(:) :: ZOELNC  !<
    real, allocatable, dimension(:) :: ZOMLNG  !<
    real, allocatable, dimension(:) :: ZOELNG  !<
    real, allocatable, dimension(:) :: ZOMLCS  !<
    real, allocatable, dimension(:) :: ZOELCS  !<
    real, allocatable, dimension(:) :: ZOMLNS  !<
    real, allocatable, dimension(:) :: ZOELNS  !<
    real, allocatable, dimension(:) :: TRSNOWC !<
    real, allocatable, dimension(:) :: CHCAP   !<
    real, allocatable, dimension(:) :: CHCAPS  !<
    real, allocatable, dimension(:) :: GZEROC  !<
    real, allocatable, dimension(:) :: GZEROG  !<
    real, allocatable, dimension(:) :: GZROCS  !<
    real, allocatable, dimension(:) :: GZROGS  !<
    real, allocatable, dimension(:) :: G12C    !<
    real, allocatable, dimension(:) :: G12G    !<
    real, allocatable, dimension(:) :: G12CS   !<
    real, allocatable, dimension(:) :: G12GS   !<
    real, allocatable, dimension(:) :: G23C    !<
    real, allocatable, dimension(:) :: G23G    !<
    real, allocatable, dimension(:) :: G23CS   !<
    real, allocatable, dimension(:) :: G23GS   !<
    real, allocatable, dimension(:) :: QFREZC  !<
    real, allocatable, dimension(:) :: QFREZG  !<
    real, allocatable, dimension(:) :: QMELTC  !<
    real, allocatable, dimension(:) :: QMELTG  !<
    real, allocatable, dimension(:) :: EVAPC   !<
    real, allocatable, dimension(:) :: EVAPCG  !<
    real, allocatable, dimension(:) :: EVAPG   !<
    real, allocatable, dimension(:) :: EVAPCS  !<
    real, allocatable, dimension(:) :: EVPCSG  !<
    real, allocatable, dimension(:) :: EVAPGS  !<
    real, allocatable, dimension(:) :: TCANO   !<
    real, allocatable, dimension(:) :: TCANS   !<
    real, allocatable, dimension(:) :: RAICAN  !<
    real, allocatable, dimension(:) :: SNOCAN  !<
    real, allocatable, dimension(:) :: RAICNS  !<
    real, allocatable, dimension(:) :: SNOCNS  !<
    real, allocatable, dimension(:) :: CWLCAP  !<
    real, allocatable, dimension(:) :: CWFCAP  !<
    real, allocatable, dimension(:) :: CWLCPS  !<
    real, allocatable, dimension(:) :: CWFCPS  !<
    real, allocatable, dimension(:) :: TSNOCS  !<
    real, allocatable, dimension(:) :: TSNOGS  !<
    real, allocatable, dimension(:) :: RHOSCS  !<
    real, allocatable, dimension(:) :: RHOSGS  !<
    real, allocatable, dimension(:) :: WSNOCS  !<
    real, allocatable, dimension(:) :: WSNOGS  !<
    real, allocatable, dimension(:) :: TPONDC  !<
    real, allocatable, dimension(:) :: TPONDG  !<
    real, allocatable, dimension(:) :: TPNDCS  !<
    real, allocatable, dimension(:) :: TPNDGS  !<
    real, allocatable, dimension(:) :: ZPLMCS  !<
    real, allocatable, dimension(:) :: ZPLMGS  !<
    real, allocatable, dimension(:) :: ZPLIMC  !<
    real, allocatable, dimension(:) :: ZPLIMG  !<
    !
    !     * DIAGNOSTIC ARRAYS USED FOR CHECKING ENERGY AND WATER
    !     * BALANCES.
    !
    real, allocatable, dimension(:) :: CTVSTP !<
    real, allocatable, dimension(:) :: CTSSTP !<
    real, allocatable, dimension(:) :: CT1STP !<
    real, allocatable, dimension(:) :: CT2STP !<
    real, allocatable, dimension(:) :: CT3STP !<
    real, allocatable, dimension(:) :: WTVSTP !<
    real, allocatable, dimension(:) :: WTSSTP !<
    real, allocatable, dimension(:) :: WTGSTP !<

! These will be allocated the dimension: 'ilg, ignd'
    integer, allocatable, dimension(:,:) :: ISNDGAT !<Integer identifier associated with sand content
    real, allocatable, dimension(:,:) :: TBARGAT !<Temperature of soil layers [K]
    real, allocatable, dimension(:,:) :: THICGAT !<Volumetric frozen water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, allocatable, dimension(:,:) :: THLQGAT !<Volumetric liquid water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, allocatable, dimension(:,:) :: BIGAT   !<Clapp and Hornberger empirical “b” parameter [ ]
    real, allocatable, dimension(:,:) :: DLZWGAT !<Permeable thickness of soil layer [m]
    real, allocatable, dimension(:,:) :: GRKSGAT !<Saturated hydraulic conductivity of soil layers \f$[m s^{-1} ]\f$
    real, allocatable, dimension(:,:) :: HCPSGAT !<Volumetric heat capacity of soil particles \f$[J m^{-3} ]\f$
    real, allocatable, dimension(:,:) :: PSISGAT !<Soil moisture suction at saturation [m]
    real, allocatable, dimension(:,:) :: PSIWGAT !<Soil moisture suction at wilting point [m]
    real, allocatable, dimension(:,:) :: TCSGAT  !<Thermal conductivity of soil particles \f$[W m^{-1} K^{-1} ]\f$\
    real, allocatable, dimension(:,:) :: THFCGAT !<Field capacity \f$[m^3 m^{-3} ]\f$
    real, allocatable, dimension(:,:) :: THMGAT  !<Residual soil liquid water content remaining after freezing or evaporation \f$[m^3 m^{-3} ]\f$
    real, allocatable, dimension(:,:) :: THPGAT  !<Pore volume in soil layer \f$[m^3 m^{-3} ]\f$
    real, allocatable, dimension(:,:) :: THRGAT  !<Liquid water retention capacity for organic soil \f$[m^3 m^{-3} ]\f$
    real, allocatable, dimension(:,:) :: THRAGAT !<Fractional saturation of soil behind the wetting front [ ]
    real, allocatable, dimension(:,:) :: ZBTWGAT !<Depth to permeable bottom of soil layer [m]
    real, allocatable, dimension(:,:) :: THLWGAT !<
    real, allocatable, dimension(:,:) :: GFLXGAT !<Heat conduction between soil layers \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:,:) :: HMFGGAT !<Diagnosed energy associated with phase change of water in soil layers \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:,:) :: HTCGAT  !<Diagnosed internal energy change of soil layer due to conduction and/or change in mass \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:,:) :: QFCGAT  !<Diagnosed vapour flux from transpiration over modelled area \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:,:) :: TBARC  !<
    real, allocatable, dimension(:,:) :: TBARG  !<
    real, allocatable, dimension(:,:) :: TBARCS !<
    real, allocatable, dimension(:,:) :: TBARGS !<
    real, allocatable, dimension(:,:) :: THLIQC !<
    real, allocatable, dimension(:,:) :: THLIQG !<
    real, allocatable, dimension(:,:) :: THICEC !<
    real, allocatable, dimension(:,:) :: THICEG !<
    real, allocatable, dimension(:,:) :: FROOT  !<
    real, allocatable, dimension(:,:) :: HCPC   !<
    real, allocatable, dimension(:,:) :: HCPG   !<
    real, allocatable, dimension(:,:) :: FROOTS !<
    real, allocatable, dimension(:,:) :: TCTOPC !<
    real, allocatable, dimension(:,:) :: TCBOTC !<
    real, allocatable, dimension(:,:) :: TCTOPG !<
    real, allocatable, dimension(:,:) :: TCBOTG !<

! These will be allocated the dimension: 'ilg, ican'
    real, allocatable, dimension(:,:) :: ACIDGAT !<Optional user-specified value of canopy near-infrared albedo to override CLASS-calculated value [ ]
    real, allocatable, dimension(:,:) :: ACVDGAT !<Optional user-specified value of canopy visible albedo to override CLASS-calculated value [ ]
    real, allocatable, dimension(:,:) :: CMASGAT !<Maximum canopy mass for vegetation category \f$[kg m^{-2} ]\f$
    real, allocatable, dimension(:,:) :: HGTDGAT !<Optional user-specified values of height of vegetation categories to override CLASS-calculated values [m]
    real, allocatable, dimension(:,:) :: PAIDGAT !<Optional user-specified value of plant area indices of vegetation categories to override CLASS-calculated values [ ]
    real, allocatable, dimension(:,:) :: PAMNGAT !<Minimum plant area index of vegetation category [ ]
    real, allocatable, dimension(:,:) :: PAMXGAT !<Minimum plant area index of vegetation category [ ]
    real, allocatable, dimension(:,:) :: PSGAGAT !<Soil moisture suction coefficient for vegetation category (used in stomatal resistance calculation) [ ]
    real, allocatable, dimension(:,:) :: PSGBGAT !<Soil moisture suction coefficient for vegetation category (used in stomatal resistance calculation) [ ]
    real, allocatable, dimension(:,:) :: QA50GAT !<Reference value of incoming shortwave radiation for vegetation category (used in stomatal resistance calculation) \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:,:) :: ROOTGAT !<Maximum rooting depth of vegetation category [m]
    real, allocatable, dimension(:,:) :: RSMNGAT !<Minimum stomatal resistance of vegetation category \f$[s m^{-1} ]\f$
    real, allocatable, dimension(:,:) :: VPDAGAT !<Vapour pressure deficit coefficient for vegetation category (used in stomatal resistance calculation) [ ]
    real, allocatable, dimension(:,:) :: VPDBGAT !<Vapour pressure deficit coefficient for vegetation category (used in stomatal resistance calculation) [ ]


! These will be allocated the dimension: 'ilg, icp1'
    real, allocatable, dimension(:,:) :: ALICGAT !<Background average near-infrared albedo of vegetation category [ ]
    real, allocatable, dimension(:,:) :: ALVCGAT !<Background average visible albedo of vegetation category [ ]
    real, allocatable, dimension(:,:) :: FCANGAT !<Maximum fractional coverage of modelled area by vegetation category [ ]
    real, allocatable, dimension(:,:) :: LNZ0GAT !<Natural logarithm of maximum roughness length of vegetation category [ ]

! These will be allocated the dimension: 'ilg, nbs'
    real, allocatable, dimension(:,:) :: FSDBGAT !<
    real, allocatable, dimension(:,:) :: FSFBGAT !<
    real, allocatable, dimension(:,:) :: FSSBGAT !<
    real, allocatable, dimension(:,:) :: SALBGAT !<
    real, allocatable, dimension(:,:) :: CSALGAT !<
    real, allocatable, dimension(:,:) :: ALTG    !<
    real, allocatable, dimension(:,:) :: ALSNO   !<
    real, allocatable, dimension(:,:) :: TRSNOWG !<

! These will be allocated the dimension: 'ilg, 4'
    real, allocatable, dimension(:,:) :: TSFSGAT !<Ground surface temperature over subarea [K]

! These will be allocated the dimension: 'ilg, 6, 50'
    integer, allocatable, dimension(:,:,:) :: ITCTGAT !<Counter of number of iterations required to solve surface energy balance for the elements of the four subareas

end type class_gather

type (class_gather), allocatable, save, target :: class_gat

! ================================================================================

type class_rotated

! These will be allocated the dimension: 'nlat'

    real, allocatable, dimension(:) :: ALIRACC !<Diagnosed total near-infrared albedo of land surface [ ]
    real, allocatable, dimension(:) :: ALVSACC !<Diagnosed total visible albedo of land surface [ ]
    real, allocatable, dimension(:) :: EVAPACC !<Diagnosed total surface water vapour flux over modelled area \f$[kg m^{-2} ]\f$
    real, allocatable, dimension(:) :: FLINACC !<Downwelling longwave radiation above surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: FLUTACC !<Upwelling longwave radiation from surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: FSINACC !<Downwelling shortwave radiation above surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: GROACC  !<Vegetation growth index [ ]
    real, allocatable, dimension(:) :: GTACC   !<Diagnosed effective surface black-body temperature [K]
    real, allocatable, dimension(:) :: HFSACC  !<Diagnosed total surface sensible heat flux over modelled area \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: HMFNACC !<Diagnosed energy associated with phase change of water in snow pack \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: OVRACC  !<Overland flow from top of soil column \f$[kg m^{-2} ]\f$
    real, allocatable, dimension(:) :: PREACC  !<Surface precipitation rate \f$[kg m^{-2} ]\f$
    real, allocatable, dimension(:) :: PRESACC !<Surface air pressure [Pa]
    real, allocatable, dimension(:) :: QAACC   !<Specific humidity at reference height \f$[kg kg^{-1} ]\f$
    real, allocatable, dimension(:) :: QEVPACC !<Diagnosed total surface latent heat flux over modelled area \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: RCANACC !<Intercepted liquid water stored on canopy \f$[kg m^{-2} ]\f$
    real, allocatable, dimension(:) :: RHOSACC !<Density of snow \f$[kg m^{-3} ]\f$
    real, allocatable, dimension(:) :: ROFACC  !<Total runoff from soil \f$[kg m^{-2} ]\f$
    real, allocatable, dimension(:) :: SCANACC !<Intercepted frozen water stored on canopy \f$[kg m^{-2} ]\f$
    real, allocatable, dimension(:) :: SNOACC  !<Mass of snow pack \f$[kg m^{-2} ]\f$
    real, allocatable, dimension(:) :: TAACC   !<Air temperature at reference height [K]
    real, allocatable, dimension(:) :: TCANACC !<Vegetation canopy temperature [K]
    real, allocatable, dimension(:) :: TSNOACC !<Snowpack temperature [K]
    real, allocatable, dimension(:) :: UVACC   !<Wind speed \f$[m s^{-1} ]\f$
    real, allocatable, dimension(:) :: WSNOACC !<Liquid water content of snow pack \f$[kg m^{-2} ]\f$
    real, allocatable, dimension(:) :: WTBLACC !<Depth of water table in soil [m]
    real, allocatable, dimension(:) :: ALTOTACC!<Broadband albedo [-]
    real, allocatable, dimension(:) :: CANARE  !<
    real, allocatable, dimension(:) :: SNOARE  !<
    real, allocatable, dimension(:) :: CSZROW  !<
    real, allocatable, dimension(:) :: DLONROW !<
    real, allocatable, dimension(:) :: DLATROW !<
    real, allocatable, dimension(:) :: FCLOROW !<
    real, allocatable, dimension(:) :: FDLROW  !<
    real, allocatable, dimension(:) :: FSIHROW !<
    real, allocatable, dimension(:) :: FSVHROW !<
    real, allocatable, dimension(:) :: GCROW   !<Type identifier for grid cell (1 = sea ice, 0 = ocean, -1 = land)
    real, allocatable, dimension(:) :: GGEOROW !<
    real, allocatable, dimension(:) :: PADRROW !<
    real, allocatable, dimension(:) :: PREROW  !<
    real, allocatable, dimension(:) :: PRESROW !<
    real, allocatable, dimension(:) :: QAROW   !<
    real, allocatable, dimension(:) :: RADJROW !<
    real, allocatable, dimension(:) :: RHOAROW !<
    real, allocatable, dimension(:) :: RHSIROW !<
    real, allocatable, dimension(:) :: RPCPROW !<
    real, allocatable, dimension(:) :: RPREROW !<Rainfall rate over modelled area \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: SPCPROW !<
    real, allocatable, dimension(:) :: SPREROW !<Snowfall rate over modelled area \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: TAROW   !<
    real, allocatable, dimension(:) :: TADPROW !<
    real, allocatable, dimension(:) :: TRPCROW !<
    real, allocatable, dimension(:) :: TSPCROW !<
    real, allocatable, dimension(:) :: ULROW   !<
    real, allocatable, dimension(:) :: VLROW   !<
    real, allocatable, dimension(:) :: VMODROW !<
    real, allocatable, dimension(:) :: VPDROW  !<
    real, allocatable, dimension(:) :: ZBLDROW !<
    real, allocatable, dimension(:) :: ZDHROW  !<
    real, allocatable, dimension(:) :: ZDMROW  !<
    real, allocatable, dimension(:) :: ZRFHROW !<
    real, allocatable, dimension(:) :: ZRFMROW !<
    real, allocatable, dimension(:) :: UVROW   !<
    real, allocatable, dimension(:) :: XDIFFUS !<
    real, allocatable, dimension(:) :: Z0ORROW !<
    real, allocatable, dimension(:) :: FSSROW  !< Shortwave radiation \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: PRENROW !<
    real, allocatable, dimension(:) :: CLDTROW !<
    real, allocatable, dimension(:) :: FSGROL  !<
    real, allocatable, dimension(:) :: FLGROL  !<
    real, allocatable, dimension(:) :: GUSTROL !<
    real, allocatable, dimension(:) :: DEPBROW !<
    real, allocatable, dimension(:) :: ALIRROW !<
    real, allocatable, dimension(:) :: ALVSROW !<
    real, allocatable, dimension(:) :: CDHROW  !<
    real, allocatable, dimension(:) :: CDMROW  !<
    real, allocatable, dimension(:) :: DRROW   !<
    real, allocatable, dimension(:) :: EFROW   !<
    real, allocatable, dimension(:) :: FLGGROW !<
    real, allocatable, dimension(:) :: FLGSROW !<
    real, allocatable, dimension(:) :: FLGVROW !<
    real, allocatable, dimension(:) :: FSGGROW !<
    real, allocatable, dimension(:) :: FSGSROW !<
    real, allocatable, dimension(:) :: FSGVROW !<
    real, allocatable, dimension(:) :: FSNOROW !<
    real, allocatable, dimension(:) :: GAROW   !<
    real, allocatable, dimension(:) :: GTROW   !<
    real, allocatable, dimension(:) :: HBLROW  !<
    real, allocatable, dimension(:) :: HEVCROW !<
    real, allocatable, dimension(:) :: HEVGROW !<
    real, allocatable, dimension(:) :: HEVSROW !<
    real, allocatable, dimension(:) :: HFSROW  !<
    real, allocatable, dimension(:) :: HFSCROW !<
    real, allocatable, dimension(:) :: HFSGROW !<
    real, allocatable, dimension(:) :: HFSSROW !<
    real, allocatable, dimension(:) :: HMFCROW !<
    real, allocatable, dimension(:) :: HMFNROW !<
    real, allocatable, dimension(:) :: HTCCROW !<
    real, allocatable, dimension(:) :: HTCSROW !<
    real, allocatable, dimension(:) :: ILMOROW !<
    real, allocatable, dimension(:) :: PCFCROW !<
    real, allocatable, dimension(:) :: PCLCROW !<
    real, allocatable, dimension(:) :: PCPGROW !<
    real, allocatable, dimension(:) :: PCPNROW !<
    real, allocatable, dimension(:) :: PETROW  !<
    real, allocatable, dimension(:) :: QEVPROW !<
    real, allocatable, dimension(:) :: QFCFROW !<
    real, allocatable, dimension(:) :: QFCLROW !<
    real, allocatable, dimension(:) :: QFGROW  !<
    real, allocatable, dimension(:) :: QFNROW  !<
    real, allocatable, dimension(:) :: QFSROW  !<
    real, allocatable, dimension(:) :: QFXROW  !<
    real, allocatable, dimension(:) :: QGROW   !<
    real, allocatable, dimension(:) :: ROFROW  !<
    real, allocatable, dimension(:) :: ROFBROW !<
    real, allocatable, dimension(:) :: ROFCROW !<
    real, allocatable, dimension(:) :: ROFNROW !<
    real, allocatable, dimension(:) :: ROFOROW !<
    real, allocatable, dimension(:) :: ROFSROW !<
    real, allocatable, dimension(:) :: ROVGROW !<
    real, allocatable, dimension(:) :: SFCQROW !<
    real, allocatable, dimension(:) :: SFCTROW !<
    real, allocatable, dimension(:) :: SFCUROW !<
    real, allocatable, dimension(:) :: SFCVROW !<
    real, allocatable, dimension(:) :: TFXROW  !<
    real, allocatable, dimension(:) :: UEROW   !<
    real, allocatable, dimension(:) :: WTABROW !<
    real, allocatable, dimension(:) :: WTRCROW !<
    real, allocatable, dimension(:) :: WTRGROW !<
    real, allocatable, dimension(:) :: WTRSROW !<
    real, allocatable, dimension(:) :: SFRHROW !<

    ! These will be allocated the dimension: 'nlat,nmos'

    integer, allocatable, dimension(:,:) :: IGDRROT !<
    integer, allocatable, dimension(:,:) :: MIDROT  !<Mosaic tile type identifier (1 for land surface, 0 for inland lake)
    real, allocatable, DIMENSION(:,:) :: ALBSROT !<
    real, allocatable, dimension(:,:) :: CMAIROT !<
    real, allocatable, dimension(:,:) :: GROROT  !<
    real, allocatable, dimension(:,:) :: QACROT  !<
    real, allocatable, dimension(:,:) :: RCANROT !<
    real, allocatable, dimension(:,:) :: RHOSROT !<
    real, allocatable, dimension(:,:) :: SCANROT !<
    real, allocatable, dimension(:,:) :: SNOROT  !<
    real, allocatable, dimension(:,:) :: TACROT  !<
    real, allocatable, dimension(:,:) :: TBASROT !<
    real, allocatable, dimension(:,:) :: TCANROT !<
    real, allocatable, dimension(:,:) :: TPNDROT !<
    real, allocatable, dimension(:,:) :: TSNOROT !<
    real, allocatable, dimension(:,:) :: WSNOROT !<
    real, allocatable, dimension(:,:) :: ZPNDROT !<
    real, allocatable, dimension(:,:) :: REFROT  !<
    real, allocatable, dimension(:,:) :: BCSNROT !<
    real, allocatable, dimension(:,:) :: AGIDROT !<
    real, allocatable, dimension(:,:) :: AGVDROT !<
    real, allocatable, dimension(:,:) :: ALGDROT !<
    real, allocatable, dimension(:,:) :: ALGWROT !<
    real, allocatable, dimension(:,:) :: ASIDROT !<
    real, allocatable, dimension(:,:) :: ASVDROT !<
    real, allocatable, dimension(:,:) :: DRNROT  !<
    real, allocatable, dimension(:,:) :: FAREROT !<Fractional coverage of mosaic tile on modelled area
    real, allocatable, dimension(:,:) :: GRKFROT !<
    real, allocatable, dimension(:,:) :: WFCIROT !<
    real, allocatable, dimension(:,:) :: WFSFROT !<
    real, allocatable, dimension(:,:) :: XSLPROT !<
    real, allocatable, dimension(:,:) :: ZPLGROT !<
    real, allocatable, dimension(:,:) :: ZPLSROT !<
    real, allocatable, dimension(:,:) :: ZSNLROT !<
    real, allocatable, dimension(:,:) :: ZSNOROT  !<
    real, allocatable, dimension(:,:) :: ALGWVROT !<
    real, allocatable, dimension(:,:) :: ALGWNROT !<
    real, allocatable, dimension(:,:) :: ALGDVROT !<
    real, allocatable, dimension(:,:) :: ALGDNROT !<
    real, allocatable, dimension(:,:) :: EMISROT  !<
    real, allocatable, dimension(:,:) :: ALIRROT !<
    real, allocatable, dimension(:,:) :: ALVSROT !<
    real, allocatable, dimension(:,:) :: CDHROT  !<
    real, allocatable, dimension(:,:) :: CDMROT  !<
    real, allocatable, dimension(:,:) :: DRROT   !<
    real, allocatable, dimension(:,:) :: EFROT   !<
    real, allocatable, dimension(:,:) :: FLGGROT !<
    real, allocatable, dimension(:,:) :: FLGSROT !<
    real, allocatable, dimension(:,:) :: FLGVROT !<
    real, allocatable, dimension(:,:) :: FSGGROT !<
    real, allocatable, dimension(:,:) :: FSGSROT !<
    real, allocatable, dimension(:,:) :: FSGVROT !<
    real, allocatable, dimension(:,:) :: FSNOROT !<
    real, allocatable, dimension(:,:) :: GAROT   !<
    real, allocatable, dimension(:,:) :: GTROT   !<
    real, allocatable, dimension(:,:) :: HBLROT  !<
    real, allocatable, dimension(:,:) :: HEVCROT !<
    real, allocatable, dimension(:,:) :: HEVGROT !<
    real, allocatable, dimension(:,:) :: HEVSROT !<
    real, allocatable, dimension(:,:) :: HFSROT  !<
    real, allocatable, dimension(:,:) :: HFSCROT !<
    real, allocatable, dimension(:,:) :: HFSGROT !<
    real, allocatable, dimension(:,:) :: HFSSROT !<
    real, allocatable, dimension(:,:) :: HMFCROT !<
    real, allocatable, dimension(:,:) :: HMFNROT !<
    real, allocatable, dimension(:,:) :: HTCCROT !<
    real, allocatable, dimension(:,:) :: SDEPROT !<Depth to bedrock in the soil profile
    real, allocatable, dimension(:,:) :: SOCIROT  !<
    real, allocatable, dimension(:,:) :: HTCSROT !<
    real, allocatable, dimension(:,:) :: ILMOROT !<
    real, allocatable, dimension(:,:) :: PCFCROT !<
    real, allocatable, dimension(:,:) :: PCLCROT !<
    real, allocatable, dimension(:,:) :: PCPGROT !<
    real, allocatable, dimension(:,:) :: PCPNROT !<
    real, allocatable, dimension(:,:) :: PETROT  !<
    real, allocatable, dimension(:,:) :: QEVPROT !<
    real, allocatable, dimension(:,:) :: QFCFROT !<
    real, allocatable, dimension(:,:) :: QFCLROT !<
    real, allocatable, dimension(:,:) :: QFGROT  !<
    real, allocatable, dimension(:,:) :: QFNROT  !<
    real, allocatable, dimension(:,:) :: QFSROT  !<
    real, allocatable, dimension(:,:) :: QFXROT  !<
    real, allocatable, dimension(:,:) :: QGROT   !<
    real, allocatable, dimension(:,:) :: ROFROT  !<
    real, allocatable, dimension(:,:) :: ROFBROT !<
    real, allocatable, dimension(:,:) :: ROFCROT !<
    real, allocatable, dimension(:,:) :: ROFNROT !<
    real, allocatable, dimension(:,:) :: ROFOROT !<
    real, allocatable, dimension(:,:) :: ROFSROT !<
    real, allocatable, dimension(:,:) :: ROVGROT !<
    real, allocatable, dimension(:,:) :: SFCQROT !<
    real, allocatable, dimension(:,:) :: SFCTROT !<
    real, allocatable, dimension(:,:) :: SFCUROT !<
    real, allocatable, dimension(:,:) :: SFCVROT !<
    real, allocatable, dimension(:,:) :: TFXROT  !<
    real, allocatable, dimension(:,:) :: TROBROT !<
    real, allocatable, dimension(:,:) :: TROFROT !<
    real, allocatable, dimension(:,:) :: TROOROT !<
    real, allocatable, dimension(:,:) :: TROSROT !<
    real, allocatable, dimension(:,:) :: UEROT   !<
    real, allocatable, dimension(:,:) :: WTABROT !<
    real, allocatable, dimension(:,:) :: WTRCROT !<
    real, allocatable, dimension(:,:) :: WTRGROT !<
    real, allocatable, dimension(:,:) :: WTRSROT !<
    real, allocatable, dimension(:,:) :: SFRHROT !<

    ! There will be allocated the dimension: 'nlat,nmos,ignd'
    integer, allocatable, dimension(:,:,:) :: ISNDROT !<
    real, allocatable, dimension(:,:,:) :: TBARROT !<
    real, allocatable, dimension(:,:,:) :: THICROT !<
    real, allocatable, dimension(:,:,:) :: THLQROT !<
    real, allocatable, dimension(:,:,:) :: BIROT   !<
    real, allocatable, dimension(:,:,:) :: DLZWROT !<
    real, allocatable, dimension(:,:,:) :: GRKSROT !<
    real, allocatable, dimension(:,:,:) :: HCPSROT !<
    real, allocatable, dimension(:,:,:) :: SANDROT !<Percentage sand content of soil
    real, allocatable, dimension(:,:,:) :: CLAYROT !<Percentage clay content of soil
    real, allocatable, dimension(:,:,:) :: ORGMROT !<Percentage organic matter content of soil
    real, allocatable, dimension(:,:,:) :: PSISROT !<
    real, allocatable, dimension(:,:,:) :: PSIWROT !<
    real, allocatable, dimension(:,:,:) :: TCSROT  !<
    real, allocatable, dimension(:,:,:) :: THFCROT !<
    real, allocatable, dimension(:,:,:) :: THMROT  !<
    real, allocatable, dimension(:,:,:) :: THPROT  !<
    real, allocatable, dimension(:,:,:) :: THRROT  !<
    real, allocatable, dimension(:,:,:) :: THRAROT !<
    real, allocatable, dimension(:,:,:) :: ZBTWROT !<
    real, allocatable, dimension(:,:,:) :: THLWROT !<
    real, allocatable, dimension(:,:,:) :: GFLXROT !<
    real, allocatable, dimension(:,:,:) :: HMFGROT !<
    real, allocatable, dimension(:,:,:) :: HTCROT  !<
    real, allocatable, dimension(:,:,:) :: QFCROT  !<

    ! These will be allocated the dimension: 'nlat,nmos,ican'
    real, allocatable, dimension(:,:,:) :: ACIDROT !<
    real, allocatable, dimension(:,:,:) :: ACVDROT !<
    real, allocatable, dimension(:,:,:) :: CMASROT !<
    real, allocatable, dimension(:,:,:) :: HGTDROT !<
    real, allocatable, dimension(:,:,:) :: PAIDROT !<
    real, allocatable, dimension(:,:,:) :: PAMNROT !<
    real, allocatable, dimension(:,:,:) :: PAMXROT !<
    real, allocatable, dimension(:,:,:) :: PSGAROT !<
    real, allocatable, dimension(:,:,:) :: PSGBROT !<
    real, allocatable, dimension(:,:,:) :: QA50ROT !<
    real, allocatable, dimension(:,:,:) :: ROOTROT !<
    real, allocatable, dimension(:,:,:) :: RSMNROT !<
    real, allocatable, dimension(:,:,:) :: VPDAROT !<
    real, allocatable, dimension(:,:,:) :: VPDBROT !<

    ! These will be allocated the dimension: 'nlat,nmos,icp1'
    real, allocatable, dimension(:,:,:) :: ALICROT !<
    real, allocatable, dimension(:,:,:) :: ALVCROT !<
    real, allocatable, dimension(:,:,:) :: FCANROT !<
    real, allocatable, dimension(:,:,:) :: LNZ0ROT !<

    ! These will be allocated the dimension: 'nlat,nmos,nbs'
    real, allocatable, dimension(:,:,:)  :: SALBROT  !<
    real, allocatable, dimension(:,:,:)  :: CSALROT  !<

    ! These will be allocated the dimension: 'nlat,nbs'
    real, allocatable, dimension(:,:) :: FSDBROL  !<
    real, allocatable, dimension(:,:) :: FSFBROL  !<
    real, allocatable, dimension(:,:) :: FSSBROL  !<

    ! These will be allocated the dimension: 'nlat,ignd'

    real, allocatable, dimension(:,:) :: TBARACC !<Temperature of soil layers [K]
    real, allocatable, dimension(:,:) :: THALACC !<Total volumetric water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, allocatable, dimension(:,:) :: THICACC !<Volumetric frozen water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, allocatable, dimension(:,:) :: THLQACC !<Volumetric liquid water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, allocatable, dimension(:,:) :: GFLXROW !<
    real, allocatable, dimension(:,:) :: HMFGROW !<
    real, allocatable, dimension(:,:) :: HTCROW  !<
    real, allocatable, dimension(:,:) :: QFCROW  !<

    ! These will be allocated the dimension: 'nlat,nmos,6,50'
    integer, allocatable, dimension(:,:,:,:) :: ITCTROT !<

    ! These will be allocated the dimension: 'nlat,nmos,4'
    real, allocatable, dimension(:,:,:)  :: TSFSROT !<

end type class_rotated

type (class_rotated), allocatable, save, target :: class_rot

contains

!=================================================================================

subroutine alloc_class_vars(nlat,nmos,ignd)

use ctem_params, only : ican, nbs, icp1

implicit none

integer, intent(in) :: nlat
integer, intent(in) :: nmos
integer, intent(in) :: ignd

integer :: ilg

!-----------

ilg = nlat * nmos

!ilg
allocate(class_gat% ILMOS   (ilg))
allocate(class_gat% JLMOS   (ilg))
allocate(class_gat% IWMOS   (ilg))
allocate(class_gat% JWMOS   (ilg))
allocate(class_gat%IGDRGAT  (ilg))
allocate(class_gat% DELZ    (ilg))
allocate(class_gat% ZBOT    (ilg))
allocate(class_gat% ALBSGAT (ilg))
allocate(class_gat% CMAIGAT (ilg))
allocate(class_gat% GROGAT  (ilg))
allocate(class_gat% QACGAT  (ilg))
allocate(class_gat%RCANGAT(ilg))
allocate(class_gat%RHOSGAT(ilg))
allocate(class_gat%SCANGAT(ilg))
allocate(class_gat%SNOGAT(ilg))
allocate(class_gat%TACGAT(ilg))
allocate(class_gat%TBASGAT(ilg))
allocate(class_gat%TCANGAT(ilg))
allocate(class_gat%TPNDGAT(ilg))
allocate(class_gat%TSNOGAT(ilg))
allocate(class_gat%WSNOGAT(ilg))
allocate(class_gat%ZPNDGAT(ilg))
allocate(class_gat%REFGAT(ilg))
allocate(class_gat%BCSNGAT(ilg))
allocate(class_gat%AGIDGAT(ilg))
allocate(class_gat%AGVDGAT(ilg))
allocate(class_gat%ALGDGAT(ilg))
allocate(class_gat%ALGWGAT(ilg))
allocate(class_gat%ASIDGAT(ilg))
allocate(class_gat%ASVDGAT(ilg))
allocate(class_gat%DRNGAT(ilg))
allocate(class_gat%GRKFGAT(ilg))
allocate(class_gat%WFCIGAT(ilg))
allocate(class_gat%WFSFGAT(ilg))
allocate(class_gat%XSLPGAT(ilg))
allocate(class_gat%ZPLGGAT(ilg))
allocate(class_gat%ZPLSGAT(ilg))
allocate(class_gat%ZSNLGAT(ilg))
allocate(class_gat%ALGWVGAT(ilg))
allocate(class_gat%ALGWNGAT(ilg))
allocate(class_gat%ALGDVGAT(ilg))
allocate(class_gat%ALGDNGAT(ilg))
allocate(class_gat%EMISGAT(ilg))
allocate(class_gat%CSZGAT(ilg))
allocate(class_gat%DLONGAT(ilg))
allocate(class_gat%DLATGAT(ilg))
allocate(class_gat%FCLOGAT(ilg))
allocate(class_gat%FDLGAT(ilg))
allocate(class_gat%FSIHGAT(ilg))
allocate(class_gat%FSVHGAT(ilg))
allocate(class_gat%GGEOGAT(ilg))
allocate(class_gat%PADRGAT(ilg))
allocate(class_gat%PREGAT(ilg))
allocate(class_gat%PRESGAT(ilg))
allocate(class_gat%QAGAT(ilg))
allocate(class_gat%RADJGAT(ilg))
allocate(class_gat%RHOAGAT(ilg))
allocate(class_gat%RHSIGAT(ilg))
allocate(class_gat%RPCPGAT(ilg))
allocate(class_gat%SPCPGAT(ilg))
allocate(class_gat%TAGAT(ilg))
allocate(class_gat%TADPGAT(ilg))
allocate(class_gat%TRPCGAT(ilg))
allocate(class_gat%TSPCGAT(ilg))
allocate(class_gat%ULGAT(ilg))
allocate(class_gat%VLGAT(ilg))
allocate(class_gat%VMODGAT(ilg))
allocate(class_gat%VPDGAT(ilg))
allocate(class_gat%Z0ORGAT(ilg))
allocate(class_gat%ZBLDGAT(ilg))
allocate(class_gat%ZDHGAT(ilg))
allocate(class_gat%ZDMGAT(ilg))
allocate(class_gat%ZRFHGAT(ilg))
allocate(class_gat%ZRFMGAT(ilg))
allocate(class_gat%FSGGAT(ilg))
allocate(class_gat%FLGGAT(ilg))
allocate(class_gat%GUSTGAT(ilg))
allocate(class_gat%DEPBGAT(ilg))
allocate(class_gat%GTBS(ilg))
allocate(class_gat%SFCUBS(ilg))
allocate(class_gat%SFCVBS(ilg))
allocate(class_gat%USTARBS(ilg))
allocate(class_gat%TCSNOW(ilg))
allocate(class_gat%GSNOW(ilg))
allocate(class_gat%ALIRGAT(ilg))
allocate(class_gat%ALVSGAT(ilg))
allocate(class_gat%CDHGAT(ilg))
allocate(class_gat%CDMGAT(ilg))
allocate(class_gat%DRGAT(ilg))
allocate(class_gat%EFGAT(ilg))
allocate(class_gat%FLGGGAT(ilg))
allocate(class_gat%FLGSGAT(ilg))
allocate(class_gat%FLGVGAT(ilg))
allocate(class_gat%FSGGGAT(ilg))
allocate(class_gat% FSGSGAT (ilg))
allocate(class_gat% FSGVGAT (ilg))
allocate(class_gat% FSNOGAT (ilg))
allocate(class_gat% GAGAT   (ilg))
allocate(class_gat% GTGAT   (ilg))
allocate(class_gat% HBLGAT  (ilg))
allocate(class_gat% HEVCGAT (ilg))
allocate(class_gat% HEVGGAT (ilg))
allocate(class_gat% HEVSGAT (ilg))
allocate(class_gat% HFSGAT  (ilg))
allocate(class_gat% HFSCGAT (ilg))
allocate(class_gat% HFSGGAT (ilg))
allocate(class_gat% HFSSGAT (ilg))
allocate(class_gat% HMFCGAT (ilg))
allocate(class_gat% HMFNGAT (ilg))
allocate(class_gat% HTCCGAT (ilg))
allocate(class_gat% HTCSGAT (ilg))
allocate(class_gat% ILMOGAT (ilg))
allocate(class_gat% PCFCGAT (ilg))
allocate(class_gat% PCLCGAT (ilg))
allocate(class_gat% PCPGGAT (ilg))
allocate(class_gat% PCPNGAT (ilg))
allocate(class_gat% PETGAT  (ilg))
allocate(class_gat% QEVPGAT (ilg))
allocate(class_gat% QFCFGAT (ilg))
allocate(class_gat% QFCLGAT (ilg))
allocate(class_gat% QFGGAT  (ilg))
allocate(class_gat% QFNGAT  (ilg))
allocate(class_gat% QFSGAT  (ilg))
allocate(class_gat% QFXGAT  (ilg))
allocate(class_gat% QGGAT   (ilg))
allocate(class_gat% ROFGAT  (ilg))
allocate(class_gat% ROFBGAT (ilg))
allocate(class_gat% ROFCGAT (ilg))
allocate(class_gat% ROFNGAT (ilg))
allocate(class_gat% ROFOGAT (ilg))
allocate(class_gat% ROFSGAT (ilg))
allocate(class_gat% ROVGGAT (ilg))
allocate(class_gat% SFCQGAT (ilg))
allocate(class_gat% SFCTGAT (ilg))
allocate(class_gat% SFCUGAT (ilg))
allocate(class_gat% SFCVGAT (ilg))
allocate(class_gat% TFXGAT  (ilg))
allocate(class_gat% TROBGAT (ilg))
allocate(class_gat% TROFGAT (ilg))
allocate(class_gat% TROOGAT (ilg))
allocate(class_gat% TROSGAT (ilg))
allocate(class_gat% UEGAT   (ilg))
allocate(class_gat% WTABGAT (ilg))
allocate(class_gat% WTRCGAT (ilg))
allocate(class_gat% WTRGGAT (ilg))
allocate(class_gat% WTRSGAT (ilg))
allocate(class_gat% QLWOGAT (ilg))
allocate(class_gat% SFRHGAT (ilg))
allocate(class_gat% FTEMP   (ilg))
allocate(class_gat% FVAP    (ilg))
allocate(class_gat% RIB     (ilg))
allocate(class_gat% FC      (ilg))
allocate(class_gat% FG      (ilg))
allocate(class_gat% FCS     (ilg))
allocate(class_gat% FGS     (ilg))
allocate(class_gat% RBCOEF  (ilg))
allocate(class_gat% ZSNOW   (ilg))
allocate(class_gat% FSVF    (ilg))
allocate(class_gat% FSVFS   (ilg))
allocate(class_gat% ALVSCN  (ilg))
allocate(class_gat% ALIRCN  (ilg))
allocate(class_gat% ALVSG   (ilg))
allocate(class_gat% ALIRG   (ilg))
allocate(class_gat% ALVSCS  (ilg))
allocate(class_gat% ALIRCS  (ilg))
allocate(class_gat% ALVSSN  (ilg))
allocate(class_gat% ALIRSN  (ilg))
allocate(class_gat% ALVSGC  (ilg))
allocate(class_gat% ALIRGC  (ilg))
allocate(class_gat% ALVSSC  (ilg))
allocate(class_gat% ALIRSC  (ilg))
allocate(class_gat% TRVSCN  (ilg))
allocate(class_gat% TRIRCN  (ilg))
allocate(class_gat%TRVSCS  (ilg))
allocate(class_gat%TRIRCS  (ilg))
allocate(class_gat%RC      (ilg))
allocate(class_gat%RCS     (ilg))
allocate(class_gat%FRAINC  (ilg))
allocate(class_gat%FSNOWC  (ilg))
allocate(class_gat%FRAICS  (ilg))
allocate(class_gat%FSNOCS  (ilg))
allocate(class_gat%CMASSC  (ilg))
allocate(class_gat%CMASCS  (ilg))
allocate(class_gat%DISP    (ilg))
allocate(class_gat%DISPS   (ilg))
allocate(class_gat%ZOMLNC  (ilg))
allocate(class_gat%ZOELNC  (ilg))
allocate(class_gat%ZOMLNG  (ilg))
allocate(class_gat%ZOELNG  (ilg))
allocate(class_gat%ZOMLCS  (ilg))
allocate(class_gat%ZOELCS  (ilg))
allocate(class_gat%ZOMLNS  (ilg))
allocate(class_gat%ZOELNS  (ilg))
allocate(class_gat%TRSNOWC (ilg))
allocate(class_gat%CHCAP   (ilg))
allocate(class_gat%CHCAPS  (ilg))
allocate(class_gat%GZEROC  (ilg))
allocate(class_gat%GZEROG  (ilg))
allocate(class_gat%GZROCS  (ilg))
allocate(class_gat%GZROGS  (ilg))
allocate(class_gat%G12C    (ilg))
allocate(class_gat%G12G    (ilg))
allocate(class_gat%G12CS   (ilg))
allocate(class_gat%G12GS   (ilg))
allocate(class_gat%G23C    (ilg))
allocate(class_gat%G23G    (ilg))
allocate(class_gat%G23CS   (ilg))
allocate(class_gat%G23GS   (ilg))
allocate(class_gat%QFREZC  (ilg))
allocate(class_gat%QFREZG  (ilg))
allocate(class_gat%QMELTC  (ilg))
allocate(class_gat%QMELTG  (ilg))
allocate(class_gat%EVAPC   (ilg))
allocate(class_gat%EVAPCG  (ilg))
allocate(class_gat%EVAPG   (ilg))
allocate(class_gat%EVAPCS  (ilg))
allocate(class_gat%EVPCSG  (ilg))
allocate(class_gat%EVAPGS  (ilg))
allocate(class_gat%TCANO   (ilg))
allocate(class_gat%TCANS   (ilg))
allocate(class_gat%RAICAN  (ilg))
allocate(class_gat%SNOCAN  (ilg))
allocate(class_gat%RAICNS  (ilg))
allocate(class_gat%SNOCNS  (ilg))
allocate(class_gat%CWLCAP  (ilg))
allocate(class_gat%CWFCAP  (ilg))
allocate(class_gat%CWLCPS  (ilg))
allocate(class_gat%CWFCPS  (ilg))
allocate(class_gat%TSNOCS  (ilg))
allocate(class_gat%TSNOGS  (ilg))
allocate(class_gat%RHOSCS  (ilg))
allocate(class_gat%RHOSGS  (ilg))
allocate(class_gat%WSNOCS  (ilg))
allocate(class_gat%WSNOGS  (ilg))
allocate(class_gat%TPONDC  (ilg))
allocate(class_gat%TPONDG  (ilg))
allocate(class_gat%TPNDCS  (ilg))
allocate(class_gat%TPNDGS  (ilg))
allocate(class_gat%ZPLMCS  (ilg))
allocate(class_gat%ZPLMGS  (ilg))
allocate(class_gat%ZPLIMC  (ilg))
allocate(class_gat%ZPLIMG  (ilg))
allocate(class_gat%CTVSTP (ilg))
allocate(class_gat%CTSSTP (ilg))
allocate(class_gat%CT1STP (ilg))
allocate(class_gat%CT2STP (ilg))
allocate(class_gat%CT3STP (ilg))
allocate(class_gat%WTVSTP (ilg))
allocate(class_gat%WTSSTP (ilg))
allocate(class_gat%WTGSTP (ilg))

! These will be allocated the dimension: 'ilg, ignd'
allocate(class_gat% ISNDGAT (ilg,ignd))
allocate(class_gat% TBARGAT (ilg,ignd))
allocate(class_gat% THICGAT (ilg,ignd))
allocate(class_gat% THLQGAT (ilg,ignd))
allocate(class_gat% BIGAT   (ilg,ignd))
allocate(class_gat% DLZWGAT (ilg,ignd))
allocate(class_gat% GRKSGAT (ilg,ignd))
allocate(class_gat% HCPSGAT (ilg,ignd))
allocate(class_gat% PSISGAT (ilg,ignd))
allocate(class_gat% PSIWGAT (ilg,ignd))
allocate(class_gat% TCSGAT  (ilg,ignd))
allocate(class_gat% THFCGAT (ilg,ignd))
allocate(class_gat% THMGAT  (ilg,ignd))
allocate(class_gat% THPGAT  (ilg,ignd))
allocate(class_gat% THRGAT  (ilg,ignd))
allocate(class_gat% THRAGAT (ilg,ignd))
allocate(class_gat% ZBTWGAT (ilg,ignd))
allocate(class_gat% THLWGAT (ilg,ignd))
allocate(class_gat% GFLXGAT (ilg,ignd))
allocate(class_gat% HMFGGAT (ilg,ignd))
allocate(class_gat% HTCGAT  (ilg,ignd))
allocate(class_gat% QFCGAT  (ilg,ignd))
allocate(class_gat% TBARC  (ilg,ignd))
allocate(class_gat% TBARG  (ilg,ignd))
allocate(class_gat% TBARCS (ilg,ignd))
allocate(class_gat% TBARGS (ilg,ignd))
allocate(class_gat% THLIQC (ilg,ignd))
allocate(class_gat% THLIQG (ilg,ignd))
allocate(class_gat% THICEC (ilg,ignd))
allocate(class_gat% THICEG (ilg,ignd))
allocate(class_gat% FROOT  (ilg,ignd))
allocate(class_gat% HCPC   (ilg,ignd))
allocate(class_gat% HCPG   (ilg,ignd))
allocate(class_gat% FROOTS (ilg,ignd))
allocate(class_gat% TCTOPC (ilg,ignd))
allocate(class_gat% TCBOTC (ilg,ignd))
allocate(class_gat% TCTOPG (ilg,ignd))
allocate(class_gat% TCBOTG (ilg,ignd))

! These will be allocated the dimension: 'ilg, ican'
allocate(class_gat% ACIDGAT (ilg,ican))
allocate(class_gat% ACVDGAT (ilg,ican))
allocate(class_gat% CMASGAT (ilg,ican))
allocate(class_gat% HGTDGAT (ilg,ican))
allocate(class_gat% PAIDGAT (ilg,ican))
allocate(class_gat% PAMNGAT (ilg,ican))
allocate(class_gat% PAMXGAT (ilg,ican))
allocate(class_gat% PSGAGAT (ilg,ican))
allocate(class_gat% PSGBGAT (ilg,ican))
allocate(class_gat% QA50GAT (ilg,ican))
allocate(class_gat% ROOTGAT (ilg,ican))
allocate(class_gat% RSMNGAT (ilg,ican))
allocate(class_gat% VPDAGAT (ilg,ican))
allocate(class_gat% VPDBGAT (ilg,ican))


! These will be allocated the dimension: 'ilg, icp1'
allocate(class_gat% ALICGAT (ilg,icp1))
allocate(class_gat% ALVCGAT (ilg,icp1))
allocate(class_gat% FCANGAT (ilg,icp1))
allocate(class_gat% LNZ0GAT (ilg,icp1))

! These will be allocated the dimension: 'ilg, nbs'
allocate(class_gat% FSDBGAT (ilg,nbs))
allocate(class_gat% FSFBGAT (ilg,nbs))
allocate(class_gat% FSSBGAT (ilg,nbs))
allocate(class_gat% SALBGAT (ilg,nbs))
allocate(class_gat% CSALGAT (ilg,nbs))
allocate(class_gat% ALTG    (ilg,nbs))
allocate(class_gat% ALSNO   (ilg,nbs))
allocate(class_gat% TRSNOWG (ilg,nbs))

! These will be allocated the dimension: 'ilg, 4'
allocate(class_gat% TSFSGAT (ilg,4))

! These will be allocated the dimension: 'ilg, 6, 50'
allocate(class_gat% ITCTGAT (ilg,6,50))

! -----------------------------------------------------------
! Now allocate the class_rot structure:

! These will be allocated the dimension: 'nlat'

allocate(class_rot% ALIRACC (nlat))
allocate(class_rot% ALVSACC (nlat))
allocate(class_rot% EVAPACC (nlat))
allocate(class_rot% FLINACC (nlat))
allocate(class_rot% FLUTACC (nlat))
allocate(class_rot% FSINACC (nlat))
allocate(class_rot% GROACC  (nlat))
allocate(class_rot% GTACC   (nlat))
allocate(class_rot% HFSACC  (nlat))
allocate(class_rot% HMFNACC (nlat))
allocate(class_rot% OVRACC  (nlat))
allocate(class_rot% PREACC  (nlat))
allocate(class_rot% PRESACC (nlat))
allocate(class_rot% QAACC   (nlat))
allocate(class_rot% QEVPACC (nlat))
allocate(class_rot% RCANACC (nlat))
allocate(class_rot% RHOSACC (nlat))
allocate(class_rot% ROFACC  (nlat))
allocate(class_rot% SCANACC (nlat))
allocate(class_rot% SNOACC  (nlat))
allocate(class_rot% TAACC   (nlat))
allocate(class_rot% TCANACC (nlat))
allocate(class_rot% TSNOACC (nlat))
allocate(class_rot% UVACC   (nlat))
allocate(class_rot% WSNOACC (nlat))
allocate(class_rot% WTBLACC (nlat))
allocate(class_rot% ALTOTACC(nlat))
allocate(class_rot% CANARE  (nlat))
allocate(class_rot% SNOARE  (nlat))
allocate(class_rot% CSZROW  (nlat))
allocate(class_rot% DLONROW (nlat))
allocate(class_rot% DLATROW (nlat))
allocate(class_rot% FCLOROW (nlat))
allocate(class_rot% FDLROW  (nlat))
allocate(class_rot% FSIHROW (nlat))
allocate(class_rot% FSVHROW (nlat))
allocate(class_rot% GCROW   (nlat))
allocate(class_rot% GGEOROW (nlat))
allocate(class_rot% PADRROW (nlat))
allocate(class_rot% PREROW  (nlat))
allocate(class_rot% PRESROW (nlat))
allocate(class_rot% QAROW   (nlat))
allocate(class_rot% RADJROW (nlat))
allocate(class_rot% RHOAROW (nlat))
allocate(class_rot% RHSIROW (nlat))
allocate(class_rot% RPCPROW (nlat))
allocate(class_rot% RPREROW (nlat))
allocate(class_rot% SPCPROW (nlat))
allocate(class_rot% SPREROW (nlat))
allocate(class_rot% TAROW   (nlat))
allocate(class_rot% TADPROW (nlat))
allocate(class_rot% TRPCROW (nlat))
allocate(class_rot% TSPCROW (nlat))
allocate(class_rot% ULROW   (nlat))
allocate(class_rot% VLROW   (nlat))
allocate(class_rot% VMODROW (nlat))
allocate(class_rot% VPDROW  (nlat))
allocate(class_rot% ZBLDROW (nlat))
allocate(class_rot% ZDHROW  (nlat))
allocate(class_rot% ZDMROW  (nlat))
allocate(class_rot% ZRFHROW (nlat))
allocate(class_rot% ZRFMROW (nlat))
allocate(class_rot% UVROW   (nlat))
allocate(class_rot% XDIFFUS (nlat))
allocate(class_rot% Z0ORROW (nlat))
allocate(class_rot% FSSROW  (nlat))
allocate(class_rot% PRENROW (nlat))
allocate(class_rot% CLDTROW (nlat))
allocate(class_rot% FSGROL  (nlat))
allocate(class_rot% FLGROL  (nlat))
allocate(class_rot% GUSTROL (nlat))
allocate(class_rot% DEPBROW (nlat))
allocate(class_rot% ALIRROW (nlat))
allocate(class_rot% ALVSROW (nlat))
allocate(class_rot% CDHROW  (nlat))
allocate(class_rot% CDMROW  (nlat))
allocate(class_rot% DRROW   (nlat))
allocate(class_rot% EFROW   (nlat))
allocate(class_rot% FLGGROW (nlat))
allocate(class_rot% FLGSROW (nlat))
allocate(class_rot% FLGVROW (nlat))
allocate(class_rot% FSGGROW (nlat))
allocate(class_rot% FSGSROW (nlat))
allocate(class_rot% FSGVROW (nlat))
allocate(class_rot% FSNOROW (nlat))
allocate(class_rot% GAROW   (nlat))
allocate(class_rot% GTROW   (nlat))
allocate(class_rot% HBLROW  (nlat))
allocate(class_rot% HEVCROW (nlat))
allocate(class_rot% HEVGROW (nlat))
allocate(class_rot% HEVSROW (nlat))
allocate(class_rot% HFSROW  (nlat))
allocate(class_rot% HFSCROW (nlat))
allocate(class_rot% HFSGROW (nlat))
allocate(class_rot% HFSSROW (nlat))
allocate(class_rot% HMFCROW (nlat))
allocate(class_rot% HMFNROW (nlat))
allocate(class_rot% HTCCROW (nlat))
allocate(class_rot% HTCSROW (nlat))
allocate(class_rot% ILMOROW (nlat))
allocate(class_rot% PCFCROW (nlat))
allocate(class_rot% PCLCROW (nlat))
allocate(class_rot% PCPGROW (nlat))
allocate(class_rot% PCPNROW (nlat))
allocate(class_rot% PETROW  (nlat))
allocate(class_rot% QEVPROW (nlat))
allocate(class_rot% QFCFROW (nlat))
allocate(class_rot% QFCLROW (nlat))
allocate(class_rot% QFGROW  (nlat))
allocate(class_rot% QFNROW  (nlat))
allocate(class_rot% QFSROW  (nlat))
allocate(class_rot% QFXROW  (nlat))
allocate(class_rot% QGROW   (nlat))
allocate(class_rot% ROFROW  (nlat))
allocate(class_rot% ROFBROW (nlat))
allocate(class_rot% ROFCROW (nlat))
allocate(class_rot% ROFNROW (nlat))
allocate(class_rot% ROFOROW (nlat))
allocate(class_rot% ROFSROW (nlat))
allocate(class_rot% ROVGROW (nlat))
allocate(class_rot% SFCQROW (nlat))
allocate(class_rot% SFCTROW (nlat))
allocate(class_rot% SFCUROW (nlat))
allocate(class_rot% SFCVROW (nlat))
allocate(class_rot% TFXROW  (nlat))
allocate(class_rot% UEROW   (nlat))
allocate(class_rot% WTABROW (nlat))
allocate(class_rot% WTRCROW (nlat))
allocate(class_rot% WTRGROW (nlat))
allocate(class_rot% WTRSROW (nlat))
allocate(class_rot% SFRHROW (nlat))

    ! These will be allocated the dimension: 'nlat,nmos'

allocate(class_rot% IGDRROT (nlat,nmos))
allocate(class_rot% MIDROT  (nlat,nmos))
allocate(class_rot% ALBSROT (nlat,nmos))
allocate(class_rot% CMAIROT (nlat,nmos))
allocate(class_rot% GROROT  (nlat,nmos))
allocate(class_rot% QACROT  (nlat,nmos))
allocate(class_rot% RCANROT (nlat,nmos))
allocate(class_rot% RHOSROT (nlat,nmos))
allocate(class_rot% SCANROT (nlat,nmos))
allocate(class_rot% SNOROT  (nlat,nmos))
allocate(class_rot% TACROT  (nlat,nmos))
allocate(class_rot% TBASROT (nlat,nmos))
allocate(class_rot% TCANROT (nlat,nmos))
allocate(class_rot% TPNDROT (nlat,nmos))
allocate(class_rot% TSNOROT (nlat,nmos))
allocate(class_rot% WSNOROT (nlat,nmos))
allocate(class_rot% ZPNDROT (nlat,nmos))
allocate(class_rot% REFROT  (nlat,nmos))
allocate(class_rot% BCSNROT (nlat,nmos))
allocate(class_rot% AGIDROT (nlat,nmos))
allocate(class_rot% AGVDROT (nlat,nmos))
allocate(class_rot% ALGDROT (nlat,nmos))
allocate(class_rot% ALGWROT (nlat,nmos))
allocate(class_rot% ASIDROT (nlat,nmos))
allocate(class_rot% ASVDROT (nlat,nmos))
allocate(class_rot% DRNROT  (nlat,nmos))
allocate(class_rot% FAREROT (nlat,nmos))
allocate(class_rot% GRKFROT (nlat,nmos))
allocate(class_rot% WFCIROT (nlat,nmos))
allocate(class_rot% WFSFROT (nlat,nmos))
allocate(class_rot% XSLPROT (nlat,nmos))
allocate(class_rot% ZPLGROT (nlat,nmos))
allocate(class_rot% ZPLSROT (nlat,nmos))
allocate(class_rot% ZSNLROT (nlat,nmos))
allocate(class_rot% ZSNOROT  (nlat,nmos))
allocate(class_rot% ALGWVROT (nlat,nmos))
allocate(class_rot% ALGWNROT (nlat,nmos))
allocate(class_rot% ALGDVROT (nlat,nmos))
allocate(class_rot% ALGDNROT (nlat,nmos))
allocate(class_rot% EMISROT  (nlat,nmos))
allocate(class_rot% ALIRROT (nlat,nmos))
allocate(class_rot% ALVSROT (nlat,nmos))
allocate(class_rot% CDHROT  (nlat,nmos))
allocate(class_rot% CDMROT  (nlat,nmos))
allocate(class_rot% DRROT   (nlat,nmos))
allocate(class_rot% EFROT   (nlat,nmos))
allocate(class_rot% FLGGROT (nlat,nmos))
allocate(class_rot% FLGSROT (nlat,nmos))
allocate(class_rot% FLGVROT (nlat,nmos))
allocate(class_rot% FSGGROT (nlat,nmos))
allocate(class_rot% FSGSROT (nlat,nmos))
allocate(class_rot% FSGVROT (nlat,nmos))
allocate(class_rot% FSNOROT (nlat,nmos))
allocate(class_rot% GAROT   (nlat,nmos))
allocate(class_rot% GTROT   (nlat,nmos))
allocate(class_rot% HBLROT  (nlat,nmos))
allocate(class_rot% HEVCROT (nlat,nmos))
allocate(class_rot% HEVGROT (nlat,nmos))
allocate(class_rot% HEVSROT (nlat,nmos))
allocate(class_rot% HFSROT  (nlat,nmos))
allocate(class_rot% HFSCROT (nlat,nmos))
allocate(class_rot% HFSGROT (nlat,nmos))
allocate(class_rot% HFSSROT (nlat,nmos))
allocate(class_rot% HMFCROT (nlat,nmos))
allocate(class_rot% HMFNROT (nlat,nmos))
allocate(class_rot% HTCCROT (nlat,nmos))
allocate(class_rot% SDEPROT (nlat,nmos))
allocate(class_rot% SOCIROT  (nlat,nmos))
allocate(class_rot% HTCSROT (nlat,nmos))
allocate(class_rot% ILMOROT (nlat,nmos))
allocate(class_rot% PCFCROT (nlat,nmos))
allocate(class_rot% PCLCROT (nlat,nmos))
allocate(class_rot% PCPGROT (nlat,nmos))
allocate(class_rot% PCPNROT (nlat,nmos))
allocate(class_rot% PETROT  (nlat,nmos))
allocate(class_rot% QEVPROT (nlat,nmos))
allocate(class_rot% QFCFROT (nlat,nmos))
allocate(class_rot% QFCLROT (nlat,nmos))
allocate(class_rot% QFGROT  (nlat,nmos))
allocate(class_rot% QFNROT  (nlat,nmos))
allocate(class_rot% QFSROT  (nlat,nmos))
allocate(class_rot% QFXROT  (nlat,nmos))
allocate(class_rot% QGROT   (nlat,nmos))
allocate(class_rot% ROFROT  (nlat,nmos))
allocate(class_rot% ROFBROT (nlat,nmos))
allocate(class_rot% ROFCROT (nlat,nmos))
allocate(class_rot% ROFNROT (nlat,nmos))
allocate(class_rot% ROFOROT (nlat,nmos))
allocate(class_rot% ROFSROT (nlat,nmos))
allocate(class_rot% ROVGROT (nlat,nmos))
allocate(class_rot% SFCQROT (nlat,nmos))
allocate(class_rot% SFCTROT (nlat,nmos))
allocate(class_rot% SFCUROT (nlat,nmos))
allocate(class_rot% SFCVROT (nlat,nmos))
allocate(class_rot% TFXROT  (nlat,nmos))
allocate(class_rot% TROBROT (nlat,nmos))
allocate(class_rot% TROFROT (nlat,nmos))
allocate(class_rot% TROOROT (nlat,nmos))
allocate(class_rot% TROSROT (nlat,nmos))
allocate(class_rot% UEROT   (nlat,nmos))
allocate(class_rot% WTABROT (nlat,nmos))
allocate(class_rot% WTRCROT (nlat,nmos))
allocate(class_rot% WTRGROT (nlat,nmos))
allocate(class_rot% WTRSROT (nlat,nmos))
allocate(class_rot% SFRHROT (nlat,nmos))

    ! There will be allocated the dimension: 'nlat,nmos,ignd'
allocate(class_rot% ISNDROT (nlat,nmos,ignd))
allocate(class_rot% TBARROT (nlat,nmos,ignd))
allocate(class_rot% THICROT (nlat,nmos,ignd))
allocate(class_rot% THLQROT (nlat,nmos,ignd))
allocate(class_rot% BIROT   (nlat,nmos,ignd))
allocate(class_rot% DLZWROT (nlat,nmos,ignd))
allocate(class_rot% GRKSROT (nlat,nmos,ignd))
allocate(class_rot% HCPSROT (nlat,nmos,ignd))
allocate(class_rot% SANDROT (nlat,nmos,ignd))
allocate(class_rot% CLAYROT (nlat,nmos,ignd))
allocate(class_rot% ORGMROT (nlat,nmos,ignd))
allocate(class_rot% PSISROT (nlat,nmos,ignd))
allocate(class_rot% PSIWROT (nlat,nmos,ignd))
allocate(class_rot% TCSROT  (nlat,nmos,ignd))
allocate(class_rot% THFCROT (nlat,nmos,ignd))
allocate(class_rot% THMROT  (nlat,nmos,ignd))
allocate(class_rot% THPROT  (nlat,nmos,ignd))
allocate(class_rot% THRROT  (nlat,nmos,ignd))
allocate(class_rot% THRAROT (nlat,nmos,ignd))
allocate(class_rot% ZBTWROT (nlat,nmos,ignd))
allocate(class_rot% THLWROT (nlat,nmos,ignd))
allocate(class_rot% GFLXROT (nlat,nmos,ignd))
allocate(class_rot% HMFGROT (nlat,nmos,ignd))
allocate(class_rot% HTCROT  (nlat,nmos,ignd))
allocate(class_rot% QFCROT  (nlat,nmos,ignd))

    ! These will be allocated the dimension: 'nlat,nmos,ican'
allocate(class_rot% ACIDROT (nlat,nmos,ican))
allocate(class_rot% ACVDROT (nlat,nmos,ican))
allocate(class_rot% CMASROT (nlat,nmos,ican))
allocate(class_rot% HGTDROT (nlat,nmos,ican))
allocate(class_rot% PAIDROT (nlat,nmos,ican))
allocate(class_rot% PAMNROT (nlat,nmos,ican))
allocate(class_rot% PAMXROT (nlat,nmos,ican))
allocate(class_rot% PSGAROT (nlat,nmos,ican))
allocate(class_rot% PSGBROT (nlat,nmos,ican))
allocate(class_rot% QA50ROT (nlat,nmos,ican))
allocate(class_rot% ROOTROT (nlat,nmos,ican))
allocate(class_rot% RSMNROT (nlat,nmos,ican))
allocate(class_rot% VPDAROT (nlat,nmos,ican))
allocate(class_rot% VPDBROT (nlat,nmos,ican))

    ! These will be allocated the dimension: 'nlat,nmos,icp1'
allocate(class_rot% ALICROT (nlat,nmos,icp1))
allocate(class_rot% ALVCROT (nlat,nmos,icp1))
allocate(class_rot% FCANROT (nlat,nmos,icp1))
allocate(class_rot% LNZ0ROT (nlat,nmos,icp1))

    ! These will be allocated the dimension: 'nlat,nmos,nbs'
allocate(class_rot% SALBROT  (nlat,nmos,nbs))
allocate(class_rot% CSALROT  (nlat,nmos,nbs))

    ! These will be allocated the dimension: 'nlat,nbs'
allocate(class_rot% FSDBROL  (nlat,nbs))
allocate(class_rot% FSFBROL  (nlat,nbs))
allocate(class_rot% FSSBROL  (nlat,nbs))

    ! These will be allocated the dimension: 'nlat,ignd'

allocate(class_rot% TBARACC (nlat,ignd))
allocate(class_rot% THALACC (nlat,ignd))
allocate(class_rot% THICACC (nlat,ignd))
allocate(class_rot% THLQACC (nlat,ignd))
allocate(class_rot% GFLXROW (nlat,ignd))
allocate(class_rot% HMFGROW (nlat,ignd))
allocate(class_rot% HTCROW  (nlat,ignd))
allocate(class_rot% QFCROW  (nlat,ignd))

    ! These will be allocated the dimension: 'nlat,nmos,6,50'
allocate(class_rot% ITCTROT (nlat,nmos,6,50))

    ! These will be allocated the dimension: 'nlat,nmos,4'
allocate(class_rot% TSFSROT (nlat,nmos,4))

end subroutine alloc_class_vars


end module class_statevars