
module creator_module_bf

implicit none

! In the creator_module you define all of the variables that will be created
! in your netcdf files. Also contains global parameters that are referenced from
! other parts of the program. 

! Joe Melton
! June 19, 2013

! Parameters
integer, parameter :: ntile = 10  	! maximum number of mosiac tiles
integer, parameter :: ctemnpft = 10      ! number of CTEM pfts  ! this includes the bare!
integer, parameter :: classnpft = 4     ! number of CLASS pfts
!integer, parameter :: cntx = 96		! number of longitudes T47
!integer, parameter :: cnty = 48		! number of latitudes T47
integer, parameter :: cntx = 128		! number of longitudes T47
integer, parameter :: cnty = 64	    	! number of latitudes T47

integer, parameter :: num_land_cells = 1955 !1958 !Remove one because 110_58, 111_57, 96_61 keeps making a glacier!!

real :: fill_value = 1.e38         !value given for empty fields in the NetCDF output files

logical, parameter :: net4 = .false. ! Set to true if you wish to make a netcdf4 output file.
                                     ! NetCDF4 allow unlimited file size and groups. However some
                                     ! things (like cview) are not compatible. Also when running in
                                     ! NetCDF3, you must have different names for each variable as
                                     ! they can not be put into groups. The present NetCDF3 setup 
                                     ! does also allow unlimited file size.

! NOTE: these coordinates are specific. If you have different coordinates,
! you need to make new matrices for them here.

real, parameter, dimension(cnty) :: valslats= &
! T47:        [ -87.16, -83.48, -79.78, -76.07, -72.36, -68.65, &
!	  -64.94, -61.23, -57.52, -53.81, -50.10, -46.39, &
!	  -42.68, -38.97, -35.26, -31.54, -27.83, -24.12, &
!	  -20.41, -16.70, -12.99,  -9.28,  -5.57,  -1.86, &
!	    1.86,   5.57,   9.28,  12.99,  16.70,  20.41, &
!	   24.12,  27.83,  31.54,  35.26,  38.97,  42.68, &
!	   46.39,  50.10,  53.81,  57.52,  61.23,  64.94, &
!	   68.65,  72.36,  76.07,  79.78,  83.48,  87.16 ] 
! T63:
    [ -87.863, -85.096, -82.312, -79.525, -76.736, -73.947, -71.157, -68.367, &
    -65.577, -62.787, -59.997, -57.206, -54.416, -51.625, -48.835, -46.044, &
    -43.254, -40.463, -37.673, -34.882, -32.091, -29.301, -26.510, -23.720, &
    -20.929, -18.138, -15.348, -12.557,  -9.767,  -6.976,  -4.185,  -1.395, &
     1.395,  4.185,  6.976,  9.767, 12.557, 15.348, 18.138, 20.929, &
    23.720, 26.510, 29.301, 32.091, 34.882, 37.673, 40.463, 43.254, &
    46.044, 48.835, 51.625, 54.416, 57.206, 59.997, 62.787, 65.577, &
    68.367, 71.157, 73.947, 76.736, 79.525, 82.312, 85.096, 87.863 ]

real, parameter, dimension(cntx) :: valslons= &
! T47:        [ 1.875, 5.625, 9.375, 13.125, 16.875, 20.625, 24.375, 28.125, 31.875, 35.625, 39.375,&
!          43.125, 46.875, 50.625, 54.375, 58.125, 61.875, 65.625, 69.375, 73.125, 76.875,     &
!          80.625, 84.375, 88.125, 91.875, 95.625, 99.375, 103.125, 106.875, 110.625, 114.375, &
!          118.125, 121.875, 125.625, 129.375, 133.125, 136.875, 140.625, 144.375, 148.125,    &
!          151.875, 155.625, 159.375, 163.125, 166.875, 170.625, 174.375, 178.125, 181.875,    &
!          185.625, 189.375, 193.125, 196.875, 200.625, 204.375, 208.125, 211.875, 215.625,    &
!          219.375, 223.125, 226.875, 230.625, 234.375, 238.125, 241.875, 245.625, 249.375,    &
!          253.125, 256.875, 260.625, 264.375, 268.125, 271.875, 275.625, 279.375, 283.125,    &
!          286.875, 290.625, 294.375, 298.125, 301.875, 305.625, 309.375, 313.125, 316.875,    &
!          320.625, 324.375, 328.125, 331.875, 335.625, 339.375, 343.125, 346.875, 350.625,    &
!          354.375, 358.125 ]
! T63:
    [ 0.0, 2.8125, 5.625, 8.4375, 11.25, 14.0625, 16.875, 19.6875, 22.5, &
    25.3125, 28.125, 30.9375, 33.75, 36.5625, 39.375, 42.1875, 45.0, 47.8125, & 
    50.625, 53.4375, 56.25, 59.0625, 61.875, 64.6875, 67.5, 70.3125, 73.125, &
    75.9375, 78.75, 81.5625, 84.375, 87.1875, 90.0, 92.8125, 95.625, 98.4375, &
    101.25, 104.0625, 106.875, 109.6875, 112.5, 115.3125, 118.125, 120.9375, &
    123.75, 126.5625, 129.375, 132.1875, 135.0, 137.8125, 140.625, 143.4375, &
    146.25, 149.0625, 151.875, 154.6875, 157.5, 160.3125, 163.125, 165.9375, &
    168.75, 171.5625, 174.375, 177.1875, 180.0, 182.8125, 185.625, 188.4375, &
    191.25, 194.0625, 196.875, 199.6875, 202.5, 205.3125, 208.125, 210.9375, &
    213.75, 216.5625, 219.375, 222.1875, 225.0, 227.8125, 230.625, 233.4375, &
    236.25, 239.0625, 241.875, 244.6875, 247.5, 250.3125, 253.125, 255.9375, &
    258.75, 261.5625, 264.375, 267.1875, 270.0, 272.8125, 275.625, 278.4375, &
    281.25, 284.0625, 286.875, 289.6875, 292.5, 295.3125, 298.125, 300.9375, &
    303.75, 306.5625, 309.375, 312.1875, 315.0, 317.8125, 320.625, 323.4375, &
    326.25, 329.0625, 331.875, 334.6875, 337.5, 340.3125, 343.125, 345.9375, &
    348.75, 351.5625, 354.375, 357.1875 ]

integer, parameter, dimension(13) :: monthend = [ 0,31,59,90,120,151,181,212,243,273,304,334,365 ] ! calender day at end of each month

! ---------------- VARIABLES ---------------- 

integer :: lon,lat,tile,pft,month,time,layer,lat_bnds,lon_bnds,bnds
integer :: varid,ncid,ncid_m
integer :: status

real, dimension(2) :: xrange
real, dimension(2) :: yrange
integer, dimension(4) :: bounds

real, allocatable, dimension(:) :: lonvect
real, allocatable, dimension(:) :: latvect

real, allocatable, dimension(:) :: latboundsvect
real, allocatable, dimension(:) :: lonboundsvect

real, allocatable, dimension(:,:) :: latbound
real, allocatable, dimension(:,:) :: lonbound

! When the arrays are declared below, they share commonalities between the
! composite and tiled formats. The distinctions between the two occur in the
! netcdf_creator and writer subroutines. We do distinguish between disturbance
! and competition here.

! The order that the variables are written below **MUST** correspond to the same order 
! they appear in the CLASSCTEM output file (ignoring the month and year).

! NOTE! Make sure that the first variable in an array has excess spacing as the xlf compiler truncates
! the arrays of characters at the number of characters corresponding to the first in the array. Sigh.

!====================== Declare arrays for monthly CLASS =============

! .OF1M_G ***********************
! MONTH YEAR  SW     LW      QH      QE    SNOACC    WSNOACC    ROFACC      PCP      EVAP      TAir
!             W/m2    W/m2    W/m2    W/m2    kg/m2   kg/m2      mm.mon    mm.mon    mm.mon    degC

integer, parameter :: numclasvars_m = 14   !number of monthly CLASS vars to write

character(100), parameter, dimension(numclasvars_m) :: CLASS_M_VAR=['SW           ','LW','QH','QE','SNOACC','WSNOACC',                 &
                                                                    'ROFACC','PCP','EVAP','TA','TRANSP','TE','GROUNDEVAP','CANOPYEVAP']

character(100), parameter, dimension(numclasvars_m) :: CLASS_M_NAME=['Shortwave radiation                      ','Longwave radiation','Sensible heat flux','Latent heat flux',                     &
                                                                     'Mass of snow pack','Liquid water content of snow pack','Total runoff from soil',                       &
                                                                     'Precipitation','Evaporation','Air temperature','Transpiration','Transpiration / Evaporation',&
                                                                     'Ground evaporation+sublimation','Canopy evaporation+sublimation']

character(100), parameter, dimension(numclasvars_m) :: CLASS_M_UNIT=['W/$m^2$        ','W/$m^2$','W/$m^2$','W/$m^2$','kg/$m^2$','kg/$m^2$','mm/month',              &
                                                                     'mm/month','mm/month','$\circ$C','mm/month','ratio','kg/$m^2$','kg/$m^2$']

! CLASS Soil vars are stored in a separate file.

! .OF2M_G **************************
!  MONTH  YEAR  TG1  THL1  THI1     TG2  THL2  THI2     TG3  THL3  THI3  ACTLYR_MO ACTLYR_MAX_MO ACTLYR_MIN_MO FTABLE_MO FTABLE_MAX_MO FTABLE_MIN_MO
!              deg  m3/m3  m3/m3   deg  m3/m3  m3/m3   deg  m3/m3  m3/m3 m m m m m m

integer, parameter :: nclassoilvars_m = 9   !number of monthly CLASS vars to write

character(100), parameter, dimension(nclassoilvars_m) :: CLASS_M_S_VAR=['TG       ','THL','THI','ACTLYR_mean','ACTLYR_max','ACTLYR_min','FTABLE_mean','FTABLE_max','FTABLE_min']

character(100), parameter, dimension(nclassoilvars_m) :: CLASS_M_S_NAME=['Ground temperature in celsius of soil layer                    ', &
                                                                     'Volumetric liquid water content of soil layer','Volumetric frozen water content of soil layer', &
                                                                     'Mean active layer depth','Maximum active layer depth','Minimum active layer depth', &
                                                                     'Mean depth to frozen water table','Maximum depth to frozen water table','Minimum depth to frozen water table']

character(100), parameter, dimension(nclassoilvars_m) :: CLASS_M_S_UNIT=['$\circ$C           ','$m^3$/$m^3$','$m^3$/$m^3$','m','m','m','m','m','m']

!======================== Declare arrays for annual CLASS =============

integer, parameter :: numclasvars_a = 9  !number of annual CLASS vars to write

! .OF1Y_G
! YEAR   SW     LW      QH      QE     ROFACC    PCP     EVAP  
!        W/m2   W/m2    W/m2    W/m2    mm.yr    mm.yr    mm.yr

character(100), parameter, dimension(numclasvars_a) :: CLASS_Y_VAR=['ANN_SW        ','ANN_LW','ANN_QH','ANN_QE','ANN_ROFACC','ANN_PCP','ANN_EVAP','ANN_TRANSP','ANN_TE']

character(100), parameter, dimension(numclasvars_a) :: CLASS_Y_NAME=['Shortwave radiation          ','Longwave radiation','Sensible heat flux',     &
                                                                     'Latent heat flux','Total runoff from soil','Precipitation',         &
                                                                     'Evaporation','Transpiration', 'Transpiration / Evaporation']

character(100), parameter, dimension(numclasvars_a) :: CLASS_Y_UNIT=['W/$m^2$          ','W/$m^2$','W/$m^2$','W/$m^2$',&
                                                                     'mm/year','mm/year','mm/year','mm/year','ratio']

!====================== Declare arrays for monthly CTEM COMPOSITE AND tiled=============

! NOTE: These do not include disturbance or competition variables. 

integer, parameter :: numctemvars_m = 14  !number of annual CTEM vars to write

! .CT01M_G
!  MONTH  YEAR  LAIMAXG  VGBIOMAS  LITTER    SOIL C    NPP       GPP        NEP       NBP    HETRES   AUTORES    LITRES   SOILRES
!                 m2/m2  Kg C/m2  Kg C/m2   Kg C/m2  gC/m2.mon  gC/m2.mon  gC/m2.mon  g/m2.mon   g/m2.mon gC/m2.mon  gC/m2.mon  gC/m2.mon

character(100), parameter, dimension(numctemvars_m) :: CTEM_M_VAR=['LAIMAXG           ','VGBIOMAS','LITTER','SOIL_C','NPP','GPP','NEP',   &
                                                                   'NBP','HETRESP','AUTORESP','LITRESP','SOILCRESP','LITRFALL','HUMIFTRS' ]

! We respecify them here for the Grid-Averaged values
character(100), parameter, dimension(numctemvars_m) :: CTEM_M_VAR_GA=['LAIMAXG_GA          ','VGBIOMAS_GA','LITTER_GA','SOIL_C_GA','NPP_GA','GPP_GA','NEP_GA',   &
                                                                   'NBP_GA','HETRESP_GA','AUTORESP_GA','LITRESP_GA','SOILCRESP_GA','LITRFALL_GA','HUMIFTRS_GA' ]

! We respecify them here for the Grid-Averaged values
character(100), parameter, dimension(numctemvars_m) :: CTEM_M_VAR_TA=['LAIMAXG_TA          ','VGBIOMAS_TA','LITTER_TA','SOIL_C_TA','NPP_TA','GPP_TA','NEP_TA',   &
                                                                   'NBP_TA','HETRESP_TA','AUTORESP_TA','LITRESP_TA','SOILCRESP_TA','LITRFALL_TA','HUMIFTRS_TA' ]

character(100), parameter, dimension(numctemvars_m) :: CTEM_M_NAME=['Maximum LAI from month                  ','Grid averaged vegetation biomass','Litter mass',                              &
                                                                    'Soil C mass','Grid averaged monthly net primary productivity',                                         &
                                                                    'Monthly gross primary productivity','Monthly net ecosystem productivity',                              &
                                                                    'Monthly net biome productivity','Monthly heterotrophic respiration',                                   &
                                                                    'Monthly autotrophic respiration','Monthly litter respiration','Monthly soil carbon respiration',&
                                                                    'Monthly transfers from veg to litter','Monthly transfers from litter to soil C' ]

character(100), parameter, dimension(numctemvars_m) :: CTEM_M_UNIT=['$m^2$/$m^2$               ','kg C/$m^2$','kg C/$m^2$','kg C/$m^2$','gC $m^{-2}$ month$^{-1}$','gC $m^{-2}$ month$^{-1}$','gC $m^{-2}$ month$^{-1}$', &
                                                                    'gC $m^{-2}$ month$^{-1}$','gC $m^{-2}$ month$^{-1}$','gC $m^{-2}$ month$^{-1}$','gC $m^{-2}$ month$^{-1}$','gC $m^{-2}$ month$^{-1}$','gC $m^{-2}$ month$^{-1}$','gC $m^{-2}$ month$^{-1}$' ]

!====DISTURBANCE================== Declare arrays for monthly CTEM COMPOSITE AND tiled =============

integer, parameter :: nctemdistvars_m = 21  !number of annual CTEM disturbance vars to write

! .CT06M_G
!  MONTH  YEAR  CO2        CO        CH4      NMHC       H2       NOX       N2O       PM25       TPM        TC        OC        BC   PROBFIRE  LUC_CO2_E  LUC_LTRIN  LUC_SOCIN   BURNFRAC
!            g/m2.mon  g/m2.mon  g/m2.mon  g/m2.mon  g/m2.mon  g/m2.mon  g/m2.mon  g/m2.mon  g/m2.mon  g/m2.mon  g/m2.mon  g/m2.mon  prod/mon    g C/m2    g C/m2    g C/m2         %

character(100), parameter, dimension(nctemdistvars_m) :: CTEM_M_D_VAR=['CO2                 ','CO','CH4','NMHC','H2','NOX',                         &
                                                                   'N2O','PM25','TPM','TC','OC','BC','PROBFIRE','LUC_CO2',&
                                                                   'LUC_LITR','LUC_SOC','BURNFRAC','BTERM','LTERM','MTERM','WIND' ]

character(100), parameter, dimension(nctemdistvars_m) :: CTEM_M_D_VAR_GA=['CO2_GA                  ','CO_GA','CH4_GA','NMHC_GA','H2_GA','NOX_GA',                         &
                                                                   'N2O_GA','PM25_GA','TPM_GA','TC_GA','OC_GA','BC_GA','PROBFIRE_GA','LUC_CO2_GA',&
                                                                   'LUC_LITR_GA','LUC_SOC_GA','BURNFRAC_GA','BTERM_GA','LTERM_GA','MTERM_GA','WIND_GA' ]

character(100), parameter, dimension(nctemdistvars_m) :: CTEM_M_D_VAR_TA=['CO2_TA                  ','CO_TA','CH4_TA','NMHC_TA','H2_TA','NOX_TA',                         &
                                                                   'N2O_TA','PM25_TA','TPM_TA','TC_TA','OC_TA','BC_TA','PROBFIRE_TA','LUC_CO2_TA',&
                                                                   'LUC_LITR_TA','LUC_SOC_TA','BURNFRAC_TA','BTERM_TA','LTERM_TA','MTERM_TA','WIND_TA' ]

character(100), parameter, dimension(nctemdistvars_m) :: CTEM_M_D_NAME=['Monthly disturbance CO2 emissions                   ','Monthly disturbance CO emissions',                                 &
                                                                    'Monthly disturbance CH4 emissions','Monthly disturbance non-CH4 hydrocarbons emissions',               &
                                                                    'Monthly disturbance H2 gas emissions','Monthly disturbance nitrogen oxides emissions',                 &
                                                                    'Monthly disturbance N2O emissions','Monthly disturbance particulate matter less than 2.5um emissions', &
                                                                    'Monthly disturbance total particulate matter emissions','Monthly disturbance total C emissions',       &
                                                                    'Monthly disturbance organic C emissions','Monthly disturbance black C emissions',                      &
                                                                    'Monthly fire probability','Monthly land use CO2 emissions','Monthly land use litter additions',        &
                                                                    'Monthly land use soil C additions','Monthly burned fraction','Monthly fire probability-Biomass',       &
                                                                    'Monthly fire probability-Lightning','Monthly fire probability-Moisture','Monthly wind' ]

character(100), parameter, dimension(nctemdistvars_m) :: CTEM_M_D_UNIT=['g CO$_2$ $m^{-2}$ month$^{-1}$                     ','g CO $m^{-2}$ month$^{-1}$','g $CH_4$ $m^{-2}$ month$^{-1}$',&
                                                                    'g NMHC $m^{-2}$ month$^{-1}$',                   &
                                                                    'g $H_2$ $m^{-2}$ month$^{-1}$','g NOx $m^{-2}$ month$^{-1}$','g $N_2O$ $m^{-2}$ month$^{-1}$',&
                                                                    'g PM2.5 $m^{-2}$ month$^{-1}$',                                 &
                                                                    'g TPM $m^{-2}$ month$^{-1}$','g total C $m^{-2}$ month$^{-1}$','g OC $m^{-2}$ month$^{-1}$',&
                                                                    'g BC $m^{-2}$ month$^{-1}$',                                 &
                                                                    'kg $CO_2$/m$^2$','% ','kg C $m^{-2}$ month$^{-1}$','kg C $m^{-2}$ month$^{-1}$','fraction   ','% ','% ','% ','m/s' ]

!====COMPETITION/LUC================== Declare arrays for monthly CTEM COMPOSITE AND tiled =============

integer, parameter :: nctemcompvars_m = 3  !number of annual CTEM COMPETITION/LUC vars to write

! .CT01M_GM
! MONTH YEAR  FRAC #1   FRAC #2   FRAC #3   FRAC #4   FRAC #5   FRAC #6   FRAC #7   FRAC #8   FRAC #9   FRAC #10   SUMCHECK
!             %         %         %         %         %         %         %         %         %         %

character(100), parameter, dimension(nctemcompvars_m) :: CTEM_M_C_VAR=['TOT_PLANT_FRAC                     ','PFT_FRAC','PFT_EXIST' ]

character(100), parameter, dimension(nctemcompvars_m) :: CTEM_M_C_NAME=[ 'Monthly Total Plant Cover                 ','Monthly PFT Percent Cover', 'Monthly PFT existence based on climate' ]

character(100), parameter, dimension(nctemcompvars_m) :: CTEM_M_C_UNIT=['Percent                 ','Percent','Boolean' ]

!======================== Declare arrays for annual CTEM COMPOSITE AND tiled=============

! NOTE: These do not include disturbance or competition variables. 

integer, parameter :: numctemvars_a = 16  !number of annual CTEM vars to write

! .CT01Y_G
! YEAR   LAIMAXG  VGBIOMAS  STEMMASS  ROOTMASS  LITRMASS  SOILCMAS  TOTCMASS  ANNUALNPP ANNUALGPP ANNUALNEP ANNUALNBP ANNHETRSP ANAUTORSP ANNLITRES ANSOILRES
!          m2/m2   Kg C/m2   Kg C/m2   Kg C/m2    Kg C/m2  Kg C/m2   Kg C/m2   gC/m2.yr  gC/m2.yr  gC/m2.yr  

character(100), parameter, dimension(numctemvars_a) :: CTEM_Y_VAR=['LAIMAXG_A               ','VGBIOMAS_A','STEMMASS_A','ROOTMASS_A',   &
                                                                   'LITRMASS_A','SOILCMAS_A','TOTCMASS','ANNUALNPP',                &
                                                                   'ANNUALGPP','ANNUALNEP','ANNUALNBP','ANNUALHETRESP',         &
                                                                   'ANNUALAUTORESP','ANNUALLITRESP','ANNUALSOILCRESP','VEGHGHT' ]

! We respecify them here for the Grid-Averaged values
character(100), parameter, dimension(numctemvars_a) :: CTEM_Y_VAR_GA=['LAIMAXG_A_GA                  ','VGBIOMAS_A_GA','STEMMASS_A_GA','ROOTMASS_A_GA',   &
                                                                   'LITRMASS_A_GA','SOILCMAS_A_GA','TOTCMASS_GA','ANNUALNPP_GA',                &
                                                                   'ANNUALGPP_GA','ANNUALNEP_GA','ANNUALNBP_GA','ANNUALHETRESP_GA',         &
                                                                   'ANNUALAUTORESP_GA','ANNUALLITRESP_GA','ANNUALSOILCRESP_GA','VEGHGHT_GA' ]

character(100), parameter, dimension(numctemvars_a) :: CTEM_Y_VAR_TA=['LAIMAXG_A_TA                  ','VGBIOMAS_A_TA','STEMMASS_A_TA','ROOTMASS_A_TA',   &
                                                                   'LITRMASS_A_TA','SOILCMAS_A_TA','TOTCMASS_TA','ANNUALNPP_TA',                &
                                                                   'ANNUALGPP_TA','ANNUALNEP_TA','ANNUALNBP_TA','ANNUALHETRESP_TA',         &
                                                                   'ANNUALAUTORESP_TA','ANNUALLITRESP_TA','ANNUALSOILCRESP_TA','VEGHGHT_TA' ]

character(100), parameter, dimension(numctemvars_a) :: CTEM_Y_NAME=['Maximum LAI from month               ','Grid averaged vegetation biomass',        &
                                                                    'Grid averaged stem biomass', 'Grid averaged root biomass',         &
                                                                    'Litter mass','Soil C mass','Total C mass',               &
                                                                    'Annual net primary productivity',                   &
                                                                    'Annual gross primary productivity',                                &
                                                                    'Annual net ecosystem productivity',                                &
                                                                    'Annual net biome productivity','Annual heterotrophic respiration', &
                                                                    'Annual autotrophic respiration','Annual litter respiration',       &
                                                                    'Annual soil carbon respiration','Vegetation height' ]

character(100), parameter, dimension(numctemvars_a) :: CTEM_Y_UNIT = ['m$^2$/m$^2$           ','kg C/m$^2$','kg C/m$^2$','kg C/m$^2$','kg C/m$^2$',            &
                                                                      'kg C/m$^2$','kg C/m$^2$','g C m$^{-2}$ year$^{-1}$','g C m$^{-2}$ year$^{-1}$',&
                                                                      'g C m$^{-2}$ year$^{-1}$',  &
                                                                      'g C m$^{-2}$ year$^{-1}$','g C m$^{-2}$ year$^{-1}$','g C m$^{-2}$ year$^{-1}$',&
                                                                      'g C m$^{-2}$ year$^{-1}$',          &
                                                                      'g C m$^{-2}$ year$^{-1}$','m' ]

!====DISTURBANCE==================== Declare arrays for annual CTEM COMPOSITE AND tiled=============

integer, parameter :: nctemdistvars_a = 20  !number of annual CTEM vars to write

! .CT06Y_G
! YEAR   ANNUALCO2  ANNUALCO  ANNUALCH4  ANN_NMHC ANNUAL_H2 ANNUALNOX ANNUALN2O  ANN_PM25  ANNUALTPM ANNUAL_TC ANNUAL_OC ANNUAL_BC APROBFIRE ANNLUCCO2  ANNLUCLTR ANNLUCSOC ABURNFRAC
!        gC/m2.yr  gC/m2.yr  gC/m2.yr  gC/m2.yr  gC/m2.yr  gC/m2.yr  gC/m2.yr  gC/m2.yr  gC/m2.yr  gC/m2.yr  gC/m2.yr  gC/m2.yr  avgprob/d   gC/m2.yr  gC/m2.yr  gC/m2.yr    %

character(100), parameter, dimension(nctemdistvars_a) :: CTEM_Y_D_VAR=['ANNUALCO2             ','ANNUALCO',                                      &
                                                                   'ANNUALCH4','ANN_NMHC','ANNUAL_H2','ANNUALNOX',              &
                                                                   'ANNUALN2O','ANN_PM25','ANNUALTPM','ANNUAL_TC',              &
                                                                   'ANNUAL_OC','ANNUAL_BC','ANNUAL_PROBFIRE','ANNUAL_LUC_CO2',  &
                                                                   'ANNUAL_LUC_LITTER','ANNUAL_LUC_SOILC','ANNUAL_BURNFRAC',    &
                                                                   'ANN_BTERM','ANN_LTERM','ANN_MTERM' ]

character(100), parameter, dimension(nctemdistvars_a) :: CTEM_Y_D_VAR_GA=['ANNUALCO2_GA             ','ANNUALCO_GA',                                      &
                                                                   'ANNUALCH4_GA','ANN_NMHC_GA','ANNUAL_H2_GA','ANNUALNOX_GA',              &
                                                                   'ANNUALN2O_GA','ANN_PM25_GA','ANNUALTPM_GA','ANNUAL_TC_GA',              &
                                                                   'ANNUAL_OC_GA','ANNUAL_BC_GA','ANNUAL_PROBFIRE_GA','ANNUAL_LUC_CO2_GA',  &
                                                                   'ANNUAL_LUC_LITTER_GA','ANNUAL_LUC_SOILC_GA','ANNUAL_BURNFRAC_GA',       &
                                                                   'ANNUAL_BTERM_GA','ANNUAL_LTERM_GA','ANNUAL_MTERM_GA'  ]

character(100), parameter, dimension(nctemdistvars_a) :: CTEM_Y_D_VAR_TA=['ANNUALCO2_TA             ','ANNUALCO_TA',                                      &
                                                                   'ANNUALCH4_TA','ANN_NMHC_TA','ANNUAL_H2_TA','ANNUALNOX_TA',              &
                                                                   'ANNUALN2O_TA','ANN_PM25_TA','ANNUALTPM_TA','ANNUAL_TC_TA',              &
                                                                   'ANNUAL_OC_TA','ANNUAL_BC_TA','ANNUAL_PROBFIRE_TA','ANNUAL_LUC_CO2_TA',  &
                                                                   'ANNUAL_LUC_LITTER_TA','ANNUAL_LUC_SOILC_TA','ANNUAL_BURNFRAC_TA',       &
                                                                   'ANNUAL_BTERM_TA','ANNUAL_LTERM_TA','ANNUAL_MTERM_TA'  ]

 
character(100), parameter, dimension(nctemdistvars_a) :: CTEM_Y_D_NAME=['Annual disturbance CO2 emissions                  ',                                 &
                                                                    'Annual disturbance CO emissions',                                  &
                                                                    'Annual disturbance CH4 emissions',                                 &
                                                                    'Annual disturbance non-CH4 hydrocarbons emissions',                &
                                                                    'Annual disturbance H2 gas emissions',                              &
                                                                    'Annual disturbance nitrogen oxides emissions',                     & 
                                                                    'Annual disturbance N2O emissions',                                 &
                                                                    'Annual disturbance particulate matter less than 2.5um emissions',  &
                                                                    'Annual disturbance total particulate matter emissions',            &
                                                                    'Annual disturbance total C emissions',                             &
                                                                    'Annual disturbance organic C emissions',                           &
                                                                    'Annual disturbance black C emissions',                             &
                                                                    'Annual fire probability','Annual land use CO2 emissions',          &
                                                                    'Annual land use litter additions',                                 &
                                                                    'Annual land use soil C additions','Annual burned frac',            &
                                                                    'Annual fire probability-Biomass','Annual fire probability-Lightning',&
                                                                    'Annual fire probability-Moisture' ]

character(100), parameter, dimension(nctemdistvars_a) :: CTEM_Y_D_UNIT = ['g/m2.yr                 ','g/m2.yr','g/m2.yr','g/m2.yr',& 
                                                                      'g/m2.yr','g/m2.yr','g/m2.yr','g/m2.yr','g/m2.yr','g/m2.yr',      &
                                                                      'g/m2.yr','g/m2.yr','Kg CO2/m2.yr','%.yr','Kg C/m2.yr','Kg C/m2.yr',' ',&
                                                                      '%.yr','%.yr','%.yr' ] 

!====COMPETITION================== Declare arrays for annual CTEM COMPOSITE AND tiled =============

integer, parameter :: nctemcompvars_a = 3  !number of annual CTEM disturbance vars to write

! .CT01Y_GM
!  YEAR   FRAC #1   FRAC #2   FRAC #3   FRAC #4   FRAC #5   FRAC #6   FRAC #7   FRAC #8   FRAC #9   FRAC #10  TOT-PLANT SUMCHECK
!         %         %         %         %         %         %         %         %         %         %           % 

character(100), parameter, dimension(nctemcompvars_a) :: CTEM_Y_C_VAR=['ANN_TOT_BARE_FRAC                     ', 'ANN_PFT_FRAC','ANN_PFT_EXIST' ]

character(100), parameter, dimension(nctemcompvars_a) :: CTEM_Y_C_NAME=[ 'Annual Total Bare Ground                  ','Annual PFT Fractional Cover', 'Annual PFT existence based on climate' ]

character(100), parameter, dimension(nctemcompvars_a) :: CTEM_Y_C_UNIT=['percent              ','percent','boolean' ]

!====== DOWETLANDS ================== Declare arryas for annual CTEM COMPOSITE

integer, parameter :: nctemwetvars_a = 6  !number of annual CTEM dowetlands vars to write
!.CT08Y

character(100), parameter, dimension(nctemwetvars_a) :: CTEM_Y_W_VAR=['CH4WET1_A                 ','CH4WET2_A','WETFDYN_A', 'CH4DYN1_A','CH4DYN2_A','SOILCH4UP_A' ]

! per tile
character(100), parameter, dimension(nctemwetvars_a) :: CTEM_Y_W_T_VAR=['CH4WET1_A_T                 ','CH4WET2_A_T','WETFDYN_A_T', 'CH4DYN1_A_T','CH4DYN2_A_T','SOILCH4UP_A_T' ]

character(100), parameter, dimension(nctemwetvars_a) :: CTEM_Y_W_NAME=['Annual Methane flux from obswetf using Hetres                          ','Annual Methane flux from obswetf using NPP',&
                                                                             'Dynamic Wetland Fraction','Annual Methane Flux using Hetres', 'Annual Methane Flux using NPP', &
                                                                             'Annual Soil uptake of Methane']

character(100), parameter, dimension(nctemwetvars_a) :: CTEM_Y_W_UNIT=['CH4/m2.yr                ','CH4/m2.yr','Fraction','CH4/m2.yr','CH4/m2.yr','CH4/m2.yr' ]

!===========Dowetlands===========Declare arrays for monthly CTEM composite 

integer, parameter :: nctemwetvars_m = 6  !number of annual CTEM dowetlands vars to write

character(100), parameter, dimension(nctemwetvars_m) :: CTEM_M_W_VAR=['CH4WET1_M                    ','CH4WET2_M','WETFDYN_M', 'CH4DYN1_M','CH4DYN2_M','SOILCH4UP_M' ]

!per tile
character(100), parameter, dimension(nctemwetvars_m) :: CTEM_M_W_T_VAR=['CH4WET1_M_T                    ','CH4WET2_M_T','WETFDYN_M_T', 'CH4DYN1_M_T','CH4DYN2_M_T','SOILCH4UP_M_T' ]

character(100), parameter, dimension(nctemwetvars_m) :: CTEM_M_W_NAME=['Monthly Methane flux from obswetf using Hetres                          ','Monthly Methane flux from obswetf using NPP',&
                                                                             'Dynamic Wetland Fraction','Monthly Methane Flux using Hetres', 'Monthly Methane Flux using NPP', &
                                                                             'Monthly Soil uptake of Methane']

character(100), parameter, dimension(nctemwetvars_m) :: CTEM_M_W_UNIT=['CH4/m2.mo                ','CH4/m2.mo','Fraction','CH4/m2.mo','CH4/m2.mo','CH4/m2.mo' ]

end module creator_module_bf
