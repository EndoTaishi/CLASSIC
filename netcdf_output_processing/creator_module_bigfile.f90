
module creator_module_bf

implicit none

! In the creator_module you define all of the variables that will be created
! in your netcdf files. Also contains global parameters that are referenced from
! other parts of the program. 

! Joe Melton
! June 19, 2013

! Parameters
integer, parameter :: ntile = 10  	! maximum number of mosiac tiles
integer, parameter :: ctemnpft = 9      ! number of CTEM pfts 
integer, parameter :: classnpft = 4     ! number of CLASS pfts
integer, parameter :: nl = 3		! number of soil layers
integer, parameter :: cntx = 96		! number of longitudes 
integer, parameter :: cnty = 48		! number of latitudes

real :: fill_value = 1.e38         !value given for empty fields in the NetCDF output files

logical, parameter :: net4 = .false. ! Set to true if you wish to make a netcdf4 output file.
                                     ! NetCDF4 allow unlimited file size and groups. However some
                                     ! things (like cview) are not compatible. Also when running in
                                     ! NetCDF3, you must have different names for each variable as
                                     ! they can not be put into groups. The present NetCDF3 setup 
                                     ! does also allow unlimited file size.

! NOTE: these coordinates are specific for a T47 run. If you have different coordinates,
! you need to make new matrices for them here.

real, parameter, dimension(cnty) :: valslats= &
        [ -87.16, -83.48, -79.78, -76.07, -72.36, -68.65, &
	  -64.94, -61.23, -57.52, -53.81, -50.10, -46.39, &
	  -42.68, -38.97, -35.26, -31.54, -27.83, -24.12, &
	  -20.41, -16.70, -12.99,  -9.28,  -5.57,  -1.86, &
	    1.86,   5.57,   9.28,  12.99,  16.70,  20.41, &
	   24.12,  27.83,  31.54,  35.26,  38.97,  42.68, &
	   46.39,  50.10,  53.81,  57.52,  61.23,  64.94, &
	   68.65,  72.36,  76.07,  79.78,  83.48,  87.16 ] 

real, parameter, dimension(cntx) :: valslons= &
        [ 1.875, 5.625, 9.375, 13.125, 16.875, 20.625, 24.375, 28.125, 31.875, 35.625, 39.375,&
          43.125, 46.875, 50.625, 54.375, 58.125, 61.875, 65.625, 69.375, 73.125, 76.875,     &
          80.625, 84.375, 88.125, 91.875, 95.625, 99.375, 103.125, 106.875, 110.625, 114.375, &
          118.125, 121.875, 125.625, 129.375, 133.125, 136.875, 140.625, 144.375, 148.125,    &
          151.875, 155.625, 159.375, 163.125, 166.875, 170.625, 174.375, 178.125, 181.875,    &
          185.625, 189.375, 193.125, 196.875, 200.625, 204.375, 208.125, 211.875, 215.625,    &
          219.375, 223.125, 226.875, 230.625, 234.375, 238.125, 241.875, 245.625, 249.375,    &
          253.125, 256.875, 260.625, 264.375, 268.125, 271.875, 275.625, 279.375, 283.125,    &
          286.875, 290.625, 294.375, 298.125, 301.875, 305.625, 309.375, 313.125, 316.875,    &
          320.625, 324.375, 328.125, 331.875, 335.625, 339.375, 343.125, 346.875, 350.625,    &
          354.375, 358.125 ]

! ---------------- VARIABLES ---------------- 

integer :: lon,lat,tile,pft,month,time,layer
integer :: varid,ncid
integer :: status

real, dimension(2) :: xrange
real, dimension(2) :: yrange
real, dimension(4) :: bounds

real, allocatable, dimension(:) :: lonvect
real, allocatable, dimension(:) :: latvect

! When the arrays are declared below, they share commonalities between the
! composite and mosaic formats. The distinctions between the two occur in the 
! netcdf_creator and writer subroutines. We do distinguish between disturbance
! and competition here.

! The order that the variables are written below **must** correspond to the same order 
! they appear in the CLASSCTEM output file (ignoring the month and year).

!====================== Declare arrays for monthly CLASS =============

! .OF1M_G ***********************
! MONTH YEAR  SW     LW      QH      QE    SNOACC    WSNOACC    ROFACC      PCP      EVAP      TAir
!             W/m2    W/m2    W/m2    W/m2    kg/m2   kg/m2      mm.mon    mm.mon    mm.mon    degC

integer, parameter :: numclasvars_m = 10   !number of monthly CLASS vars to write

character(100), parameter, dimension(numclasvars_m) :: CLASS_M_VAR=['SW','LW','QH','QE','SNOACC','WSNOACC',                 &
                                                                    'ROFACC','PCP','EVAP','TA']

character(100), parameter, dimension(numclasvars_m) :: CLASS_M_NAME=['Shortwave radiation','Longwave radiation','Sensible heat flux','Latent heat flux',                     &
                                                                     'Mass of snow pack','Liquid water content of snow pack','Total runoff from soil',                       &
                                                                     'Precipitation','Evaporation','Air temperature']

character(100), parameter, dimension(numclasvars_m) :: CLASS_M_UNIT=['W/m^2','W/m^2','W/m^2','W/m^2','kg/m^2','kg/m^2','mm.mon',              &
                                                                     'mm.mon','mm.mon','degC']

! CLASS Soil vars are stored in a separate file.

! .OF2M_G **************************
!  MONTH  YEAR  TG1  THL1  THI1     TG2  THL2  THI2     TG3  THL3  THI3
!              deg  m3/m3  m3/m3   deg  m3/m3  m3/m3   deg  m3/m3  m3/m3

integer, parameter :: nclassoilvars_m = 3   !number of monthly CLASS vars to write

character(100), parameter, dimension(nclassoilvars_m) :: CLASS_M_S_VAR=['TG','THL','THI']

character(100), parameter, dimension(nclassoilvars_m) :: CLASS_M_S_NAME=['Ground temperature in celsius of soil layer',                                                        &
                                                                     'Volumetric liquid water content of soil layer','Volumetric frozen water content of soil layer']

character(100), parameter, dimension(nclassoilvars_m) :: CLASS_M_S_UNIT=['degC','m^3/m^3','m^3/m^3']

!======================== Declare arrays for annual CLASS =============

integer, parameter :: numclasvars_a = 7  !number of annual CLASS vars to write

! .OF1Y_G
! YEAR   SW     LW      QH      QE     ROFACC    PCP     EVAP  
!        W/m2   W/m2    W/m2    W/m2    mm.yr    mm.yr    mm.yr

character(100), parameter, dimension(numclasvars_a) :: CLASS_Y_VAR=['ANN_SW','ANN_LW','ANN_QH','ANN_QE','ANN_ROFACC','ANN_PCP','ANN_EVAP']

character(100), parameter, dimension(numclasvars_a) :: CLASS_Y_NAME=['Shortwave radiation','Longwave radiation','Sensible heat flux',     &
                                                                     'Latent heat flux','Total runoff from soil','Precipitation',         &
                                                                     'Evaporation']

character(100), parameter, dimension(numclasvars_a) :: CLASS_Y_UNIT=['W/m^2','W/m^2','W/m^2','W/m^2',&
                                                                     'mm.year','mm.year','mm.year']

!====================== Declare arrays for monthly CTEM COMPOSITE AND MOSAIC=============

! NOTE: These do not include disturbance or competition variables. 

integer, parameter :: numctemvars_m = 12  !number of annual CTEM vars to write

! .CT01M_G
!  MONTH  YEAR  LAIMAXG  VGBIOMAS  LITTER    SOIL C    NPP       GPP        NEP       NBP    HETRES   AUTORES    LITRES   SOILRES
!                 m2/m2  Kg C/m2  Kg C/m2   Kg C/m2  gC/m2.mon  gC/m2.mon  gC/m2.mon  g/m2.mon   g/m2.mon gC/m2.mon  gC/m2.mon  gC/m2.mon

character(100), parameter, dimension(numctemvars_m) :: CTEM_M_VAR=['LAIMAXG','VGBIOMAS','LITTER','SOIL_C','NPP','GPP','NEP',   &
                                                                   'NBP','HETRESP','AUTORESP','LITRESP','SOILCRESP' ]

character(100), parameter, dimension(numctemvars_m) :: CTEM_M_NAME=['Maximum LAI from month','Grid averaged vegetation biomass','Litter mass',                              &
                                                                    'Soil C mass','Grid averaged monthly net primary productivity',                                         &
                                                                    'Monthly gross primary productivity','Monthly net ecosystem productivity',                              &
                                                                    'Monthly net biome productivity','Monthly heterotrophic respiration',                                   &
                                                                    'Monthly autotrophic respiration','Monthly litter respiration','Monthly soil carbon respiration' ]

character(100), parameter, dimension(numctemvars_m) :: CTEM_M_UNIT=['m^2/m^2','Kg C/m^2','Kg C/m^2','Kg C/m^2','gC/m^2.month','gC/m^2.month','gC/m^2.month', &
                                                                    'g/m^2.month','g/m^2.month','g/m^2.month','g/m^2.month','g/m^2.month' ]

!====DISTURBANCE================== Declare arrays for monthly CTEM COMPOSITE AND MOSAIC =============

integer, parameter :: nctemdistvars_m = 17  !number of annual CTEM disturbance vars to write

! .CT06M_G
!  MONTH  YEAR  CO2        CO        CH4      NMHC       H2       NOX       N2O       PM25       TPM        TC        OC        BC   PROBFIRE  LUC_CO2_E  LUC_LTRIN  LUC_SOCIN   BURNFRAC
!            g/m2.mon  g/m2.mon  g/m2.mon  g/m2.mon  g/m2.mon  g/m2.mon  g/m2.mon  g/m2.mon  g/m2.mon  g/m2.mon  g/m2.mon  g/m2.mon  prod/mon    g C/m2    g C/m2    g C/m2         %

character(100), parameter, dimension(nctemdistvars_m) :: CTEM_M_D_VAR=['CO2','CO','CH4','NMHC','H2','NOX',                         &
                                                                   'N2O','PM25','TPM','TC','OC','BC','PROBFIRE','LUC_CO2',&
                                                                   'LUC_LITR','LUC_SOC','BURNFRAC' ]

character(100), parameter, dimension(nctemdistvars_m) :: CTEM_M_D_NAME=['Monthly disturbance CO2 emissions','Monthly disturbance CO emissions',                                 &
                                                                    'Monthly disturbance CH4 emissions','Monthly disturbance non-CH4 hydrocarbons emissions',               &
                                                                    'Monthly disturbance H2 gas emissions','Monthly disturbance nitrogen oxides emissions',                 &
                                                                    'Monthly disturbance N2O emissions','Monthly disturbance particulate matter less than 2.5um emissions', &
                                                                    'Monthly disturbance total particulate matter emissions','Monthly disturbance total C emissions',       &
                                                                    'Monthly disturbance organic C emissions','Monthly disturbance black C emissions',                      &
                                                                    'Monthly fire probability','Monthly land use CO2 emissions','Monthly land use litter additions',        &
                                                                    'Monthly land use soil C additions','Monthly burned fraction' ]

character(100), parameter, dimension(nctemdistvars_m) :: CTEM_M_D_UNIT=['g/m^2.month','g/m^2.month','g/m^2.month','g/m^2.month',                   &
                                                                    'g/m^2.month','g/m^2.month','g/m^2.month','g/m^2.month',                                 &
                                                                    'g/m^2.month','g/m^2.month','g/m^2.month','g/m^2.month',                                 &
                                                                    'Kg CO2/m^2','% ','Kg C/m^2.mo','Kg C/m^2.mo','           ' ]

!====COMPETITION/LUC================== Declare arrays for monthly CTEM COMPOSITE AND MOSAIC =============

integer, parameter :: nctemcompvars_m = 3  !number of annual CTEM COMPETITION/LUC vars to write

! .CT01M_GM
! MONTH YEAR  FRAC #1   FRAC #2   FRAC #3   FRAC #4   FRAC #5   FRAC #6   FRAC #7   FRAC #8   FRAC #9   FRAC #10   SUMCHECK
!             %         %         %         %         %         %         %         %         %         %

character(100), parameter, dimension(nctemcompvars_m) :: CTEM_M_C_VAR=['TOT_PLANT_FRAC','PFT_FRAC','PFT_EXIST' ]

character(100), parameter, dimension(nctemcompvars_m) :: CTEM_M_C_NAME=[ 'Monthly Total Plant Cover','Monthly PFT Percent Cover', 'Monthly PFT existence based on climate' ]

character(100), parameter, dimension(nctemcompvars_m) :: CTEM_M_C_UNIT=['Percent','Percent','Boolean' ]

!======================== Declare arrays for annual CTEM COMPOSITE AND MOSAIC=============

! NOTE: These do not include disturbance or competition variables. 

integer, parameter :: numctemvars_a = 15  !number of annual CTEM vars to write

! .CT01Y_G
! YEAR   LAIMAXG  VGBIOMAS  STEMMASS  ROOTMASS  LITRMASS  SOILCMAS  TOTCMASS  ANNUALNPP ANNUALGPP ANNUALNEP ANNUALNBP ANNHETRSP ANAUTORSP ANNLITRES ANSOILRES
!          m2/m2   Kg C/m2   Kg C/m2   Kg C/m2    Kg C/m2  Kg C/m2   Kg C/m2   gC/m2.yr  gC/m2.yr  gC/m2.yr  

character(100), parameter, dimension(numctemvars_a) :: CTEM_Y_VAR=['LAIMAXG_A','VGBIOMAS_A','STEMMASS_A','ROOTMASS_A',   &
                                                                   'LITRMASS_A','SOILCMAS_A','TOTCMASS','ANNUALNPP',                &
                                                                   'ANNUALGPP','ANNUALNEP','ANNUALNBP','ANNUALHETRESP',         &
                                                                   'ANNUALAUTORESP','ANNUALLITRESP','ANNUALSOILCRESP' ]
 
character(100), parameter, dimension(numctemvars_a) :: CTEM_Y_NAME=['Maximum LAI from month','Grid averaged vegetation biomass',        &
                                                                    'Grid averaged stem biomass', 'Grid averaged root biomass',         &
                                                                    'Litter mass','Soil C mass','Total C mass',               &
                                                                    'Annual net primary productivity',                   &
                                                                    'Annual gross primary productivity',                                &
                                                                    'Annual net ecosystem productivity',                                &
                                                                    'Annual net biome productivity','Annual heterotrophic respiration', &
                                                                    'Annual autotrophic respiration','Annual litter respiration',       &
                                                                    'Annual soil carbon respiration' ]

character(100), parameter, dimension(numctemvars_a) :: CTEM_Y_UNIT = ['m^2/m^2','Kg C/m^2','Kg C/m^2','Kg C/m^2','Kg C/m^2',            &
                                                                      'Kg C/m^2','Kg C/m^2','gC/m^2.year','gC/m^2.year','gC/m^2.year',  &
                                                                      'gC/m^2.year','gC/m^2.year','gC/m^2.year','gC/m^2.year',          &
                                                                      'gC/m^2.year' ] 

!====DISTURBANCE==================== Declare arrays for annual CTEM COMPOSITE AND MOSAIC=============

integer, parameter :: nctemdistvars_a = 17  !number of annual CTEM vars to write

! .CT06Y_G
! YEAR   ANNUALCO2  ANNUALCO  ANNUALCH4  ANN_NMHC ANNUAL_H2 ANNUALNOX ANNUALN2O  ANN_PM25  ANNUALTPM ANNUAL_TC ANNUAL_OC ANNUAL_BC APROBFIRE ANNLUCCO2  ANNLUCLTR ANNLUCSOC ABURNFRAC
!        gC/m2.yr  gC/m2.yr  gC/m2.yr  gC/m2.yr  gC/m2.yr  gC/m2.yr  gC/m2.yr  gC/m2.yr  gC/m2.yr  gC/m2.yr  gC/m2.yr  gC/m2.yr  avgprob/d   gC/m2.yr  gC/m2.yr  gC/m2.yr    %

character(100), parameter, dimension(nctemdistvars_a) :: CTEM_Y_D_VAR=['ANNUALCO2','ANNUALCO',                                      &
                                                                   'ANNUALCH4','ANN_NMHC','ANNUAL_H2','ANNUALNOX',              &
                                                                   'ANNUALN2O','ANN_PM25','ANNUALTPM','ANNUAL_TC',              &
                                                                   'ANNUAL_OC','ANNUAL_BC','ANNUAL_PROBFIRE','ANNUAL_LUC_CO2',  &
                                                                   'ANNUAL_LUC_LITTER','ANNUAL_LUC_SOILC','ANNUAL_BURNFRAC' ]
 
character(100), parameter, dimension(nctemdistvars_a) :: CTEM_Y_D_NAME=['Annual disturbance CO2 emissions',                                 &
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
                                                                    'Annual land use soil C additions','Annual burned frac' ]

character(100), parameter, dimension(nctemdistvars_a) :: CTEM_Y_D_UNIT = ['g/m2.yr','g/m2.yr','g/m2.yr','g/m2.yr',& 
                                                                      'g/m2.yr','g/m2.yr','g/m2.yr','g/m2.yr','g/m2.yr','g/m2.yr',      &
                                                                      'g/m2.yr','g/m2.yr','Kg CO2/m2.yr','%.yr','Kg C/m2.yr','Kg C/m2.yr',' ' ] 

!====COMPETITION================== Declare arrays for annual CTEM COMPOSITE AND MOSAIC =============

integer, parameter :: nctemcompvars_a = 3  !number of annual CTEM disturbance vars to write

! .CT01Y_GM
!  YEAR   FRAC #1   FRAC #2   FRAC #3   FRAC #4   FRAC #5   FRAC #6   FRAC #7   FRAC #8   FRAC #9   FRAC #10  TOT-PLANT SUMCHECK
!         %         %         %         %         %         %         %         %         %         %           % 

character(100), parameter, dimension(nctemcompvars_a) :: CTEM_Y_C_VAR=['ANN_TOT_BARE_FRAC', 'ANN_PFT_FRAC','ANN_PFT_EXIST' ]

character(100), parameter, dimension(nctemcompvars_a) :: CTEM_Y_C_NAME=[ 'Annual Total Bare Ground','Annual PFT Fractional Cover', 'Annual PFT existence based on climate' ]

character(100), parameter, dimension(nctemcompvars_a) :: CTEM_Y_C_UNIT=['percent','percent','boolean' ]



end module creator_module_bf
