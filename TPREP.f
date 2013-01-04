      SUBROUTINE TPREP(THLIQC, THLIQG, THICEC, THICEG, TBARC,  TBARG,             
     1                 TBARCS, TBARGS, HCPC,   HCPG,   TCTOPC, TCBOTC,
     2                 TCTOPG, TCBOTG, HCPSCS, HCPSGS, TCSNOW, TSNOCS, 
     3                 TSNOGS, WSNOCS, WSNOGS, RHOSCS, RHOSGS, TCANO,  
     4                 TCANS,  CEVAP,  IEVAP,  TBAR1P, WTABLE, ZERO,
     5                 EVAPC,  EVAPCG, EVAPG,  EVAPCS, EVPCSG, EVAPGS,            
     6                 GSNOWC, GSNOWG, GZEROC, GZEROG, GZROCS, GZROGS,
     7                 QMELTC, QMELTG, EVAP,
     8                 TPONDC, TPONDG, TPNDCS, TPNDGS, QSENSC, QSENSG, 
     9                 QEVAPC, QEVAPG, TACCO,  QACCO,  TACCS,  QACCS,  
     A                 ILMOX,  UEX,    HBLX,
     B                 ILMO,   UE,     HBL,    
     C                 ST,     SU,     SV,     SQ,     CDH,    CDM,              
     D                 QSENS,  QEVAP,  QLWAVG, 
     E                 FSGV,   FSGS,   FSGG,   FLGV,   FLGS,   FLGG,   
     F                 HFSC,   HFSS,   HFSG,   HEVC,   HEVS,   HEVG,   
     G                 HMFC,   HMFN,   QFCF,   QFCL,   EVPPOT, ACOND,  
     H                 DRAG,   THLIQ,  THICE,  TBAR,   ZPOND,  TPOND,  
     I                 THPOR,  THLMIN, THLRET, THFC,   HCPS,   TCS,    
     J                 TA,     RHOSNO, TSNOW,  ZSNOW,  WSNOW,  TCAN,
     K                 FC,     FCS,    DELZ,   DELZW,  ZBOTW,
     L                 ISAND,  ILG,    IL1,    IL2,    JL,     IG,  
     M                 FVEG,   TCSATU, TCSATF, FTEMP,  FTEMPX, FVAP,
     N                 FVAPX,  RIB,    RIBX  )           
C
C     Purpose: Initialize subarea variables and calculate various 
C     parameters for surface energy budget calculations.
C
C     * NOV 24/11 - R.HARVEY.   NEW SNOW THERMAL CONDUCTIVITY FROM
C     *                         STURM ET AL. (1997).
C     * OCT 12/11 - M.LAZARE.   REMOVED TSURF.     
C     * AUG   /08 - J.P.PAQUIN. ADD CALCULATION FOR FTEMP, FVAP AND
C     *                         RIB FOR OUTPUT IN GEM (IMPLEMENTED BY
C     *                         L. DUARTE ON OCT. 28/08).
C     * MAR 20/08 - D.VERSEGHY. REMOVE TBAR3, TCTOP3, TCBOT3.
C     * DEC 12/07 - D.VERSEGHY. MAJOR REVISIONS TO CALCULATION OF
C     *                         SOIL THERMAL CONDUCTIVITY.
C     * MAY 18/06 - D.VERSEGHY. ADJUST CALCULATION OF TBAR1P FOR ROCK
C     *                         SOILS; LIMIT CALCULATION OF TBAR3(I,3)
C     *                         TO UPPER 4.1 M OF SOIL; CORRECT WTABLE
C     *                         TO ACCOUNT FOR PRESENCE OF ICE.
C     * MAR 23/06 - D.VERSEGHY. MODIFY CALCULATION OF HCPSNO TO ACCOUNT
C     *                         FOR PRESENCE OF WATER IN SNOWPACK.
C     * MAR 21/06 - P.BARTLETT. INITIALIZE ADDITIONAL VARIABLES TO ZERO.
C     * OCT 04/05 - D.VERSEGHY. NEW VARIABLES TBAR3,TCTOP3,TCBOT3.
C     * APR 08/05 - Y.DELAGE. TCTOP VARIES GRADUALLY WITH ZPOND TO TCW.
C     * MAR 16/05 - D.VERSEGHY. TREAT FROZEN SOIL WATER AS ICE
C     *                         VOLUME RATHER THAN AS EQUIVALENT 
C     *                         LIQUID WATER VOLUME; REVERSE ORDER
C     *                         IN LOOP 500.
C     * NOV 03/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
C     * AUG 05/04 - Y.DELAGE/D.VERSEGHY. INITIALIZE NEW DIAGNOSTIC
C     *                         VARIABLES ILMO, UE AND HBL.
C     * JUL 30/02 - D.VERSEGHY. MOVE CALCULATION OF VEGETATION 
C     *                         STOMATAL RESISTANCE INTO APREP
C     *                         AND CANALB; SHORTENED CLASS3
C     *                         COMMON BLOCK.
C     * JUN 17/02 - D.VERSEGHY. NEW THERMAL ARRAYS FOR SURFACE 
C     *                         TEMPERATURE ITERATION, WITH PONDED
C     *                         WATER ROLLED INTO SOIL UPPER LAYER;
C     *                         SHORTENED CLASS4 COMMON BLOCK.
C     * MAR 20/02 - D.VERSEGHY. MOVE CALCULATION OF BACKGROUND SOIL 
C     *                         PROPERTIES INTO "CLASSB"; UPDATES 
C     *                         TO MAKE ZPOND A PROGNOSTIC VARIABLE.
C     * FEB 27/02 - D.VERSEGHY. RECALCULATE WILTING POINT BASED ON
C     *                         FIELD CAPACITY.
C     * JAN 18/02 - D.VERSEGHY. INTRODUCTION OF CALCULATION OF FIELD
C     *                         CAPACITY AND NEW BARE SOIL EVAPORATION
C     *                         PARAMETERS.
C     * APR 11/01 - M.LAZARE.   SHORTENED "CLASS2" COMMON BLOCK.     
C     * NOV 01/00 - A.WU/D.VERSEGHY. EXTEND MINERAL SOIL CALCULATION
C     *                              OF SOIL EVAPORATION "BETA" TO
C     *                              ORGANIC SOILS.
C     * SEP 19/00 - A.WU/D.VERSEGHY. CHANGE CALCULATION OF THERMAL
C     *                              CONDUCTIVITY FOR ORGANIC SOILS,
C     *                              USING METHOD OF FAROUKI (1981).
C     *                              ALSO, CALCULATE STOMATAL RESISTANCE
C     *                              USING VEGETATION-VARYING 
C     *                              COEFFICIENTS FOR ENVIRONMENTAL
C     *                              VARIABLES.
C     * FEB 14/00 - D.VERSEGHY. INSERT CALCULATION OF WATER TABLE DEPTH
C     *                         FOR ORGANIC SOILS.
C     * DEC 07/99 - A.WU/D.VERSEGHY.  INCORPORATE CALCULATION OF "BETA"
C     *                               PARAMETER FOR NEW SOIL EVAPORATION
C     *                               FORMULATION.
C     * JUN 20/97 - D.VERSEGHY. CLASS - VERSION 2.7.
C     *                         CHANGES RELATED TO VARIABLE SOIL DEPTH
C     *                         (MOISTURE HOLDING CAPACITY) AND DEPTH-
C     *                         VARYING SOIL PROPERTIES.
C     * JAN 24/97 - D.VERSEGHY. CLASS - VERSION 2.6.
C     *                         SET RC AND RCS TO ZERO FOR GRID CELLS
C     *                         WITH NO VEGETATION.
C     * JAN 02/96 - D.VERSEGHY. CLASS - VERSION 2.5.
C     *                         COMPLETION OF ENERGY BALANCE 
C     *                         DIAGNOSTICS.
C     * AUG 30/95 - D.VERSEGHY. CLASS - VERSION 2.4.
C     *                         REMOVE SUBTRACTION OF RESIDUAL SOIL
C     *                         MOISTURE CONTENT IN CALCULATIONS OF
C     *                         "PSIZRO" AND "PSII".
C     * AUG 18/95 - D.VERSEGHY. REVISIONS TO ALLOW FOR INHOMOGENEITY
C     *                         BETWEEN SOIL LAYERS AND FRACTIONAL
C     *                         ORGANIC MATTER CONTENT.
C     * DEC 16/94 - D.VERSEGHY. CLASS - VERSION 2.3.
C     *                         INITIALIZE THREE NEW DIAGNOSTIC FIELDS.
C     * NOV 12/94 - D.VERSEGHY. SET INITIAL TEMPERATURE OF EMERGING
C     *                         CANOPY TO TA INSTEAD OF TO ZERO.
C     * JAN 31/94 - D.VERSEGHY. CLASS - VERSION 2.2.
C     *                         INTRODUCE LIMITING VALUES INTO
C     *                         CALCULATION OF "PSIZRO" TO AVOID
C     *                         OVERFLOWS.
C     * JUL 27/93 - D.VERSEGHY/M.LAZARE. INITIALIZE NEW DIAGNOSTIC 
C     *                                  FIELDS FSGV,FSGG,FLGV,FLGG,
C     *                                  HFSC,HFSG,HMFC.
C     * MAY 06/93 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.1.
C     *                                  MODIFICATIONS TO CANOPY
C     *                                  RESISTANCE TO ADD "RCS"
C     *                                  FIELD FOR SNOW-COVERED
C     *                                  CANOPY. 
C     * JUL 04/92 - D.VERSEGHY/M.LAZARE. REVISED AND VECTORIZED CODE
C     *                                  FOR MODEL VERSION GCM7.                            
C     * AUG 12/91 - D.VERSEGHY. CODE FOR MODEL VERSION GCM7U -
C     *                         CLASS VERSION 2.0 (WITH CANOPY).            
C     * APR 11/89 - D.VERSEGHY. PREPARATION AND INITIALIZATION FOR
C     *                         LAND SURFACE ENERGY BUDGET 
C     *                         CALCULATIONS.
C
      IMPLICIT NONE
C                                                                                 
C     * INTEGER CONSTANTS.
C
      INTEGER ILG,IL1,IL2,JL,IG,I,J  
C                                                                                 
C     * OUTPUT ARRAYS.                                                            
C
C     (Suffix CS = vegetation over snow cover; GS = bare snow cover; C 
C     or CO = vegetation over ground; G or GO = bare ground.)
C
      REAL TBARC (ILG,IG)   !Subarea temperatures of soil layers [C]
      REAL TBARG (ILG,IG)   !Subarea temperatures of soil layers [C]
      REAL TBARCS(ILG,IG)   !Subarea temperatures of soil layers [C]
      REAL TBARGS(ILG,IG)   !Subarea temperatures of soil layers [C]
      REAL THLIQC(ILG,IG)   !Liquid water content of soil layers under 
                            !vegetation [m3 m-3]
      REAL THLIQG(ILG,IG)   !Liquid water content of soil layers in bare 
                            !areas [m3 m-3]
      REAL THICEC(ILG,IG)   !Frozen water content of soil layers under 
                            !vegetation [m3 m-3]
      REAL THICEG(ILG,IG)   !Frozen water content of soil layers in bare 
                            !areas [m3 m-3]
      REAL HCPC  (ILG,IG)   !Heat capacity of soil layers under 
                            !vegetation [J m-3 K-1] (Cg)
      REAL HCPG  (ILG,IG)   !Heat capacity of soil layers in bare areas 
                            ![J m-3 K1] (Cg)
      REAL TCTOPC(ILG,IG)   !Thermal conductivity of soil at top of 
                            !layer in canopy-covered subareas 
                            ![W m-1 K-1] (lambda)
      REAL TCBOTC(ILG,IG)   !Thermal conductivity of soil at bottom of 
                            !layer in canopy-covered subareas
                            ![W m-1 K-1] (lambda)
      REAL TCTOPG(ILG,IG)   !Thermal conductivity of soil at top of 
                            !layer in bare ground subareas [W m-1 K-1] 
                            !(lambda)
      REAL TCBOTG(ILG,IG)   !Thermal conductivity of soil at bottom of 
                            !layer in bare ground subareas [W m-1 K-1] 
                            !(lambda)
C                                                                                 
      REAL HCPSCS(ILG)  !Heat capacity of snow pack under vegetation 
                        !canopy [J m -3 K1] (Cs)
      REAL HCPSGS(ILG)  !Heat capacity of snow pack in bare areas 
                        ![J m-3 K1] (Cs)
      REAL TCSNOW(ILG)  !Thermal conductivity of snow [W m-1 K-1] 
      REAL TSNOGS(ILG)  !Temperature of snow pack in bare areas [K] 
      REAL TSNOCS(ILG)  !Temperature of snow pack under vegetation 
                        !canopy [K]
      REAL WSNOCS(ILG)  !Liquid water content of snow pack under 
                        !vegetation [kg m-2]
      REAL WSNOGS(ILG)  !Liquid water content of snow pack in bare areas 
                        ![kg m-2]
      REAL RHOSCS(ILG)  !Density of snow pack under vegetation canopy 
                        ![kg m-3]
      REAL RHOSGS(ILG)  !Density of snow pack in bare areas [kg m-3] 
      REAL TCANO (ILG)  !Temperature of canopy over ground [K] 
      REAL TCANS (ILG)  !Temperature of canopy over snow [K] 
      REAL CEVAP (ILG)  !Soil evaporation efficiency coefficient [ ] 
                        !(beta)
      REAL TBAR1P(ILG)  !Lumped temperature of ponded water and first 
                        !soil layer [K]
      REAL WTABLE(ILG)  !Depth of water table in soil [m] (zwt)
      REAL ZERO  (ILG)  !Dummy vector containing all zeros
      REAL TPONDC(ILG)  !Subarea temperature of surface ponded water [C] 
      REAL TPONDG(ILG)  !Subarea temperature of surface ponded water [C] 
      REAL TPNDCS(ILG)  !Subarea temperature of surface ponded water [C]
      REAL TPNDGS(ILG)  !Subarea temperature of surface ponded water [C]
C
      INTEGER             IEVAP (ILG)   !Flag indicating whether soil 
                                        !evaporation is occurring or not
C                                                                                 
C     * OUTPUT ARRAYS WHICH ARE INTERNAL WORK ARRAYS FOR CLASST                   
C     * AND ARE INITIALIZED TO ZERO HERE.                                               
C                                                                                 
      REAL EVAPC (ILG)  !Evaporation from vegetation over ground [m s-1]   
      REAL EVAPCG(ILG)  !Evaporation from ground under vegetation 
                        ![m s-1]
      REAL EVAPG (ILG)  !Evaporation from bare ground [m s-1] 
      REAL EVAPCS(ILG)  !Evaporation from vegetation over snow [m s-1]            
      REAL EVPCSG(ILG)  !Evaporation from snow under vegetation [m s-1] 
      REAL EVAPGS(ILG)  !Evaporation from snow on bare ground [m s-1] 
      REAL GSNOWC(ILG)  !Heat flux at top of snow pack under canopy 
                        ![W m-2] 
      REAL GSNOWG(ILG)  !Heat flux at top of snow pack over bare ground 
                        ![W m-2]
      REAL GZEROC(ILG)  !Subarea heat flux at soil surface [W m-2] 
      REAL GZEROG(ILG)  !Subarea heat flux at soil surface [W m-2]
      REAL GZROCS(ILG)  !Subarea heat flux at soil surface [W m-2]
      REAL GZROGS(ILG)  !Subarea heat flux at soil surface [W m-2]
      REAL QMELTC(ILG)  !Heat to be used for melting snow under canopy 
                        ![W m-2] 
      REAL QMELTG(ILG)  !Heat to be used for melting snow on bare ground 
                        ![W m-2]
      REAL EVAP  (ILG)  !Diagnosed total surface water vapour flux over 
                        !modelled area [kg m -2 s-1]
      REAL QSENSC(ILG)  !Sensible heat flux from vegetation canopy over 
                        !subarea [W m-2]
      REAL QSENSG(ILG)  !Sensible heat flux from ground over subarea 
                        ![W m-2]
      REAL QEVAPC(ILG)  !Latent heat flux from vegetation canopy over 
                        !subarea [W m-2]
      REAL QEVAPG(ILG)  !Latent heat flux from ground over subarea 
                        ![W m-2]
      REAL TACCO (ILG)  !Temperature of air within vegetation canopy 
                        !space over bare ground [K]
      REAL QACCO (ILG)  !Specific humidity of air within vegetation 
                        !canopy space over bare ground [kg kg-1]
      REAL TACCS (ILG)  !Temperature of air within vegetation canopy 
                        !space over snow [K]
      REAL QACCS (ILG)  !Specific humidity of air within vegetation 
                        !canopy space over snow [kg kg-1]
 
C
C     * DIAGNOSTIC ARRAYS.
C
      REAL ST    (ILG)  !Diagnosed screen-level air temperature [K]   
      REAL SU    (ILG)  !Diagnosed anemometer-level zonal wind [m s-1] 
      REAL SV    (ILG)  !Diagnosed anemometer-level meridional wind 
                        ![m s-1] 
      REAL SQ    (ILG)  !Diagnosed screen-level specific humidity 
                        ![kg kg-1]
      REAL CDH   (ILG)  !Surface drag coefficient for heat [ ] 
      REAL CDM   (ILG)  !Surface drag coefficient for momentum [ ]
      REAL QSENS (ILG)  !Diagnosed total surface sensible heat flux over 
                        !modelled area [W m-2]
      REAL QEVAP (ILG)  !Diagnosed total surface latent heat flux over 
                        !modelled area [W m-2]
      REAL QLWAVG(ILG)  !Upwelling longwave radiation over modelled area 
                        ![W m-2]
      REAL FSGV  (ILG)  !Diagnosed net shortwave radiation on vegetation 
                        !canopy [W m -2]
      REAL FSGS  (ILG)  !Diagnosed net shortwave radiation at snow 
                        !surface [W m-2]
      REAL FSGG  (ILG)  !Diagnosed net shortwave radiation at soil 
                        !surface [W m-2]
      REAL FLGV  (ILG)  !Diagnosed net longwave radiation on vegetation 
                        !canopy [W m-2]
      REAL FLGS  (ILG)  !Diagnosed net longwave radiation at snow 
                        !surface [W m-2]
      REAL FLGG  (ILG)  !Diagnosed net longwave radiation at soil 
                        !surface [W m-2] 
      REAL HFSC  (ILG)  !Diagnosed sensible heat flux on vegetation 
                        !canopy [W m-2]
      REAL HFSS  (ILG)  !Diagnosed sensible heat flux at snow surface 
                        ![W m-2] 
      REAL HFSG  (ILG)  !Diagnosed sensible heat flux at soil surface 
                        ![W m-2]
      REAL HEVC  (ILG)  !Diagnosed latent heat flux on vegetation canopy 
                        ![W m-2]
      REAL HEVS  (ILG)  !Diagnosed latent heat flux at snow surface 
                        ![W m-2]
      REAL HEVG  (ILG)  !Diagnosed latent heat flux at soil surface 
                        ![W m-2]
      REAL HMFC  (ILG)  !Diagnosed energy associated with phase change 
                        !of water on vegetation [W m-2]
      REAL HMFN  (ILG)  !Diagnosed energy associated with phase change 
                        !of water in snow pack [W m-2]
      REAL EVPPOT(ILG)  !Diagnosed potential evapotranspiration 
                        ![kg m-2 s-1]   
      REAL ACOND (ILG)  !Diagnosed product of drag coefficient and wind 
                        !speed over modelled area [m s-1]
      REAL DRAG  (ILG)  !Surface drag coefficient under neutral 
                        !stability [ ]
      REAL ILMO  (ILG)  !Surface drag coefficient under neutral 
                        !stability [ ]
      REAL UE    (ILG)  !Friction velocity of air [m s-1]   
      REAL HBL   (ILG)  !Height of the atmospheric boundary layer [m] 
      REAL ILMOX (ILG)  !Inverse of Monin-Obukhov roughness length over 
                        !each subarea [m-1]
      REAL UEX   (ILG)  !Friction velocity of air over each subarea 
                        ![m s-1]
      REAL HBLX  (ILG)  !Height of the atmospheric boundary layer over 
                        !each subarea [m]   

      REAL QFCF  (ILG),   QFCL  (ILG)
      REAL FTEMP (ILG),   FTEMPX(ILG),   FVAP  (ILG),
     1     FVAPX (ILG),   RIB   (ILG),   RIBX  (ILG)

C
C     * INPUT ARRAYS.                                                             
C                                                                                 
      REAL THLIQ (ILG,IG)   !Volumetric liquid water content of soil 
                            !layers [m3 m-3] (theta_l)
      REAL THICE (ILG,IG)   !Volumetric frozen water content of soil 
                            !layers [m3 m-3] (theta_i)
      REAL TBAR  (ILG,IG)   !Temperature of soil layers [K]
      REAL ZPOND (ILG)      !Depth of ponded water on surface [m]
      REAL TPOND (ILG)      !Temperature of ponded water [K]
C                                                                                 
      REAL TA    (ILG)  !Air temperature at reference height [K]   
      REAL RHOSNO(ILG)  !Density of snow [kg m-3] (rho_s) 
      REAL TSNOW (ILG)  !Snowpack temperature [K] 
      REAL ZSNOW (ILG)  !Depth of snow pack [m] (zs) 
      REAL WSNOW (ILG)  !Liquid water content of snow pack [kg m-2] (ws) 
      REAL TCAN  (ILG)  !Vegetation canopy temperature [K] 
      REAL FC    (ILG)  !Fractional coverage of canopy over bare ground 
                        !for modelled area [ ]
      REAL FCS   (ILG)  !Fractional coverage of canopy over snow for 
                        !modelled area [ ]
C
C     * SOIL PROPERTY ARRAYS.                                     
C                                                                                 
      REAL THPOR(ILG,IG)    !Pore volume in soil layer [m3 m-3] 
                            !(theta_p)
      REAL THLMIN(ILG,IG)   !Residual soil liquid water content 
                            !remaining after freezing or evaporation 
                            ![m3 m-3]
      REAL THLRET(ILG,IG)   !Liquid water retention capacity for organic 
                            !soil [m3 m-3 ] (theta_ret)
      REAL THFC  (ILG,IG)   !Field capacity [m3 m-3] (theta_fc)
      REAL HCPS  (ILG,IG)   !Heat capacity of soil material [J m-3 K-1] 
                            !(Cm)
      REAL TCS   (ILG,IG)   !Thermal conductivity of soil particles 
                            ![W m-1 K-1] (theta_s)
      REAL DELZW(ILG,IG)    !Permeable thickness of soil layer [m] 
                            !(delta_zw) 
      REAL ZBOTW(ILG,IG)    !Depth to permeable bottom of soil layer [m] 
                            !(zb,w)
      REAL DELZ(IG)         !Overall thickness of soil layer [m]
C                                                                                 
      INTEGER       ISAND (ILG,IG)  !Sand content flag
C                                                                                 
C     * INTERNAL WORK FIELDS FOR THIS ROUTINE.                                    
C                                                                                 
      REAL FVEG  (ILG),   TCSATU(ILG),   TCSATF(ILG)
C
C     * TEMPORARY VARIABLES.
C
      REAL SATRAT,THLSAT,THISAT,TCDRY,TCKAPU,TCKAPF,TCRATU,TCRATF,
     1     TCSOLU,TCSOLF,TCSOIL,TSUM1,TSUM2,ZSUM
C
C     * COMMON BLOCK PARAMETERS.
C
      REAL DELT     !Time step [s]
      REAL TFREZ    !Freezing point of water [K]
      REAL RGAS     !Gas Constant [J kg-1 K-1]
      REAL RGASV    !Gas constant for water vapour [J kg-1 K-1]
      REAL GRAV     !Acceleration due to gravity [m s-1]
      REAL SBC      !Stefan-Boltzmann constant [W m-2 K-4]
      REAL VKC      !Von Karman constant (0.40)
      REAL CT       !Drag coefficient for water (1.15*10^-3)
      REAL VMIN     !Minimum wind speed (0.1) [m s-1]
      REAL TCW      !Thermal conductivity of water (0.57) [W m-1 K-1]
      REAL TCICE    !Thermal conductivity of ice (2.24) [W m-1 K-1]
      REAL TCSAND   !Thermal conductivity of sand particles (2.5) 
                    ![W m-1 K-1]
      REAL TCCLAY   !Thermal conductivity of fine mineral particles 
                    !(2.5) [W m-1 K-1]
      REAL TCOM     !Thermal conductivity of organic matter (0.25) 
                    ![W m-1 K-1]
      REAL TCDRYS   !Thermal conductivity of dry mineral soil (0.275) 
                    ![W m-1 K-1]
      REAL RHOSOL   !Density of soil mineral matter (2.65*10^3) [kg m-3]
      REAL RHOOM    !Density of soil organic matter (1.30*10^3) [kg m-3]
      REAL HCPW     !Volumetric heat capacity of water (4.187*10^6) 
                    ![J m-3 K-1]
      REAL HCPICE   !Volumetric heat capacity of ice (1.9257*10^6) 
                    ![J m-3 K-1]
      REAL HCPSOL   !Volumetric heat capacity of mineral matter 
                    !(2.25*10^6) [J m-3 K-1]
      REAL HCPOM    !Volumetric heat capacity of organic matter 
                    !(2.50*10^6) [J m-3 K-1]
      REAL HCPSND   !Volumetric heat capacity of sand particles 
                    !(2.13*10^6) [J m-3 K-1]
      REAL HCPCLY   !Volumetric heat capacity of fine mineral particles 
                    !(2.38*10^6) [J m-3 K-1]
      REAL SPHW     !Specific heat of water (4.186*10^3) [J kg-1 K-1]
      REAL SPHICE   !Specific heat of ice (2.10*10^3) [J kg-1 K-1]
      REAL SPHVEG   !Specific heat of vegetation matter (2.70*10^3) 
                    ![J kg-1 K-1]
      REAL SPHAIR   !Specific heat of air [J kg-1 K-1]
      REAL RHOW     !Density of water (1.0*10^3) [kg m-3]
      REAL RHOICE   !Density of ice (0.917*10^3) [kg m-3]
      REAL TCGLAC   !Thermal conductivity of ice sheets (2.24) 
                    ![W m-1 K-1]
      REAL CLHMLT   !Latent heat of freezing of water (0.334*10^6) 
                    ![J kg-1]
      REAL CLHVAP   !Latent heat of vaporization of water (2.501*10^6) 
                    ![J kg-1]
C                                                                                 
      COMMON /CLASS1/ DELT,TFREZ
      COMMON /CLASS2/ RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN
      COMMON /CLASS3/ TCW,TCICE,TCSAND,TCCLAY,TCOM,TCDRYS,
     1                RHOSOL,RHOOM
      COMMON /CLASS4/ HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,
     1                SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,
     2                TCGLAC,CLHMLT,CLHVAP
C                                                                                 
C----------------------------------------------------------------------        
C
      !
      !In the first two loops, various subarea arrays and internal 
      !CLASST variables are initialized. The initial temperatures of the 
      !vegetation canopy above snow and above bare ground (TCANS and 
      !TCANO) are set to the temperature of the vegetation over the 
      !whole modelled area (TCAN) if TCAN is not effectively 0 K (the 
      !value it is assigned if vegetation is not present). Otherwise, 
      !the canopy temperatures are initialized to the air temperature 
      !TA.
      !
C     * INITIALIZE 2-D AND 3-D ARRAYS.                                                    
C                                                                                 
      DO 50 J=1,IG                                                                
      DO 50 I=IL1,IL2                                                             
          THLIQG(I,J)=THLIQ(I,J)                                                  
          THICEG(I,J)=THICE(I,J)                                                  
          THLIQC(I,J)=THLIQ(I,J)                                                  
          THICEC(I,J)=THICE(I,J)                                                  
          TBARCS(I,J)=0.0                                                         
          TBARGS(I,J)=0.0                                                         
          TBARC (I,J)=0.0                                                         
          TBARG (I,J)=0.0
          TCTOPC(I,J)=0.0
          TCBOTC(I,J)=0.0
          TCTOPG(I,J)=0.0
          TCBOTG(I,J)=0.0
   50 CONTINUE                                                                    
C                                                                                 
C     * INITIALIZE 1-D INTERNAL WORK FIELDS AND DIAGNOSTIC ARRAYS.
C                                                                                 
      DO 100 I=IL1,IL2                                                            
          FVEG  (I)=FC(I)+FCS(I)                                                
          IF(TCAN(I).GT.5.0) THEN
              TCANS (I)=TCAN(I)  
              TCANO (I)=TCAN(I) 
          ELSE
              TCANS (I)=TA(I)  
              TCANO (I)=TA(I) 
          ENDIF
          EVAPC (I)=0.                                                            
          EVAPCG(I)=0.                                                            
          EVAPG (I)=0.                                                            
          EVAPCS(I)=0.                                                            
          EVPCSG(I)=0.                                                            
          EVAPGS(I)=0.                                                            
          GSNOWC(I)=0.                                                            
          GSNOWG(I)=0.                                                            
          GZEROC(I)=0.                                                            
          GZEROG(I)=0.                                                            
          GZROCS(I)=0.
          GZROGS(I)=0.
          QMELTC(I)=0.                                                            
          QMELTG(I)=0.
          QSENSC(I)=0.
          QSENSG(I)=0.
          QEVAPC(I)=0.
          QEVAPG(I)=0.
          TPONDC(I)=0.
          TPONDG(I)=0.
          TPNDCS(I)=0.
          TPNDGS(I)=0.
          TACCS (I)=0.
          QACCS (I)=0.
          TACCO (I)=0.
          QACCO (I)=0.
          ST    (I)=0.
          SU    (I)=0.
          SV    (I)=0.
          SQ    (I)=0.
          CDH   (I)=0.
          CDM   (I)=0.
          QSENS (I)=0.
          QEVAP (I)=0.
          EVAP  (I)=0.
          QLWAVG(I)=0.
          FSGV  (I)=0.
          FSGS  (I)=0.
          FSGG  (I)=0.
          FLGV  (I)=0. 
          FLGS  (I)=0. 
          FLGG  (I)=0.
          HFSC  (I)=0.
          HFSS  (I)=0.
          HFSG  (I)=0.
          HEVC  (I)=0.
          HEVS  (I)=0.
          HEVG  (I)=0.
          HMFC  (I)=0.
          HMFN  (I)=0.
          QFCF  (I)=0.
          QFCL  (I)=0.
          EVPPOT(I)=0.
          ACOND (I)=0.
          DRAG  (I)=0.
          ILMO  (I)=0.
          UE    (I)=0.
          HBL   (I)=0.
          ILMOX (I)=0.
          UEX   (I)=0.
          HBLX  (I)=0.
          ZERO  (I)=0.
          FTEMP (I)=0.
          FVAP  (I)=0.
          RIB   (I)=0.
          FTEMPX(I)=0.
          FVAPX (I)=0.
          RIBX  (I)=0.
          WTABLE(I)=9999.
  100 CONTINUE                                                                    
C                                                                                 
C     * SURFACE EVAPORATION EFFICIENCY FOR BARE SOIL ENERGY BALANCE
C     * CALCULATIONS.                                 
C
      !
      !In loop 200 the soil surface evaporation flag IEVAP and the 
      !evaporation efficiency coefficient CEVAP are assigned. If the 
      !liquid water content of the first soil layer is effectively equal 
      !to the minimum water content THLMIN, IEVAP and CEVAP are set to 
      !zero. If the liquid water content of the first soil layer is 
      !greater than the field capacity THFC, IEVAP and CEVAP are set to 
      !unity. Otherwise, IEVAP is set to 1 and CEVAP (or beta as it is 
      !typically symbolized in the literature) is calculated using a 
      !relation presented by Lee and Pielke (1992):
      !
      !CEVAP = 0.25*[1 – cos(THLIQ*pi/THFC)]^2
      !
      !where THLIQ is the liquid water content of the first soil layer 
      !and THFC is its field capacity.
      !
      DO 200 I=IL1,IL2    
          IF(THLIQG(I,1).LT.(THLMIN(I,1)+0.001)) THEN    
              IEVAP(I)=0  
              CEVAP(I)=0.0
          ELSEIF(THLIQG(I,1).GT.THFC(I,1)) THEN
              IEVAP(I)=1   
              CEVAP(I)=1.0
          ELSE
              IEVAP(I)=1
              CEVAP(I)=0.25*(1.0-COS(3.14159*THLIQG(I,1)/THFC(I,1)))**2
          ENDIF
  200 CONTINUE  
C                                                                                 
C     * VOLUMETRIC HEAT CAPACITIES OF SOIL LAYERS.
C                   
      !
      !In loop 300 the volumetric heat capacities Cg of the soil layers 
      !under a bare surface (HCPG) and under vegetation (HCPC) are 
      !calculated, from their respective liquid and frozen water 
      !contents THLIQ and THICE:
      !
      !Cg = HCPW*THLIQ + HCPICE*THICE + HCPS(1 - THPOR)
      !
      !where HCPS is the heat capacity of the soil matter and THPOR is 
      !the pore volume. (The heat capacity of air is neglected.)
      !                                                              
      DO 300 J=1,IG                                                               
      DO 300 I=IL1,IL2                                                            
          IF(ISAND(I,1).GT.-4)                                     THEN          
              HCPG(I,J)=HCPW*THLIQG(I,J)+HCPICE*THICEG(I,J)+
     1            HCPS(I,J)*(1.0-THPOR(I,J))
              HCPC(I,J)=HCPW*THLIQC(I,J)+HCPICE*THICEC(I,J)+
     1            HCPS(I,J)*(1.0-THPOR(I,J))
          ELSE                                                                    
              HCPC(I,J)=HCPICE                                                    
              HCPG(I,J)=HCPICE                                                    
          ENDIF                                                                   
  300 CONTINUE                                                                    
C                                                                                 
C     * THERMAL PROPERTIES OF SNOW.
C                            
      !
      !In loop 400, the thermal properties of the snow pack under the 
      !vegetation canopy and over bare soil are assigned on the basis of 
      !the properties of the snow pack over the whole modelled area. The 
      !heat capacity of the snow pack Cs is calculated from the volume 
      !fractions of snow particles and liquid water in the pack. The 
      !former is obtained from the ratio of the densities of the snow 
      !pack and ice, and the latter from the ratio of the liquid water 
      !content, normalized by the snow depth, and the density of water:
      !
      !Cs = HCPICE*[RHOSNO/RHOICE] + HCPW*WSNOW/[RHOW*ZSNOW]
      !
      !The thermal conductivity of snow TCSNOW is obtained from the 
      !snow density using an empirical relationship derived by Sturm et 
      !al. (1997):
      !
      !TCSNOW = (3.233*10^-6)*RHOSNO^2 – (1.01*10^-3)*RHOSNO + 0.138 
      !RHOSNO >= 156.0
      !TCSNOW = (0.234*10-3)*RHOSNO + 0.023      RHOSNO < 156.0
      !                                                     
      DO 400 I=IL1,IL2                                                            
          IF(ZSNOW(I).GT.0.)                                        THEN          
              HCPSCS(I)=HCPICE*RHOSNO(I)/RHOICE+HCPW*WSNOW(I)/
     1            (RHOW*ZSNOW(I)) 
              HCPSGS(I)=HCPSCS(I)
C             TCSNOW(I)=2.576E-6*RHOSNO(I)*RHOSNO(I)+0.074                        
              IF(RHOSNO(I).LT.156.0) THEN
                  TCSNOW(I)=0.234E-3*RHOSNO(I)+0.023
              ELSE
                  TCSNOW(I)=3.233E-6*RHOSNO(I)*RHOSNO(I)-1.01E-3*
     1                RHOSNO(I)+0.138
              ENDIF
              IF(FVEG(I).LT.1.)                                 THEN              
                  TSNOGS(I)=TSNOW(I)                                              
                  WSNOGS(I)=WSNOW(I)
                  RHOSGS(I)=RHOSNO(I)
              ELSE                                                                
                  TSNOGS(I)=0.0                                                   
                  WSNOGS(I)=0.0
                  RHOSGS(I)=0.0
              ENDIF                                                               
              IF(FVEG(I).GT.0.)                                 THEN              
                  TSNOCS(I)=TSNOW(I)                                              
                  WSNOCS(I)=WSNOW(I)
                  RHOSCS(I)=RHOSNO(I)
              ELSE                                                                
                  TSNOCS(I)=0.0                                                   
                  WSNOCS(I)=0.0
                  RHOSCS(I)=0.0
              ENDIF                                                               
          ELSE                                                                    
              TSNOGS(I)=0.0                                                       
              WSNOGS(I)=0.0
              RHOSGS(I)=0.0
              TSNOCS(I)=0.0                                                       
              WSNOCS(I)=0.0
              RHOSCS(I)=0.0
              TCSNOW(I)=0.0
          ENDIF                                                                   
  400 CONTINUE                                                                    
C                
C     * THERMAL CONDUCTIVITIES OF SOIL LAYERS AND DEPTH OF WATER
C     * TABLE IN ORGANIC SOILS.                                         
C        
      !In loop 500, the thermal conductivities of the soil layers are 
      !assigned. If the ISAND flag for the first soil layer is -4 
      !(indicating glacier or ice sheet), or if the ISAND flag is -3 
      !(indicating rock), then literature values for glacier ice or sand 
      !particles respectively are assigned. If the ISAND flag is equal 
      !to -2, indicating organic soil, the depth of the water table 
      !WTABLE is first calculated. This is taken to lie within the first 
      !layer, counting from the bottom of the soil profile, in which the 
      !soil water content is larger than the retention capacity THLRET. 
      !The water table depth is deduced by assuming that the soil is 
      !saturated below the water table, and that the water content is at 
      !the retention capacity above it. Thus, if THLIQ + THICE = THPOR 
      !for the soil layer, the water table is located at the top of the 
      !soil layer; if THLIQ + THICE = THLRET, it is located at the 
      !permeable bottom of the soil layer; and if THLIQ + THICE is 
      !between these two values, its location is given by:
      !
      !WTABLE = ZBOTW - DELZW*[(THLIQ + THICE - THLRET)/(THPOR-THLRET)]
      !
      !where DELZW is the permeable thickness of the soil layer.
      !
      !The thermal conductivities of organic and mineral soils are 
      !calculated following the analysis of Côté and Konrad (2005). 
      !They model the soil thermal conductivity using the concept of a 
      !relative thermal conductivity lambda_r which has a value of 0 for 
      !dry soils and 1 at saturation:
      !
      !lambda = [ lambda_sat – lambda_dry ]*lambda_r + lambda_dry
      !
      !The relative thermal conductivity is obtained from the degree of 
      !saturation (the water content divided by the pore volume) Sr, 
      !using the following generalized relationship:
      !
      !lambda_r = kappa*Sr/[1 + (kappa-1)*Sr ]
      !
      !The empirical coefficient kappa takes the following values:
      !
      !Unfrozen coarse mineral soils:   kappa = 4.0
      !Frozen coarse mineral soils:     kappa = 1.2
      !Unfrozen fine mineral soils:     kappa = 1.9
      !Frozen fine mineral soils:       kappa = 0.85
      !Unfrozen organic soils:          kappa = 0.6
      !Frozen organic soils:            kappa = 0.25
      !
      !The dry thermal conductivity lambda_dry is calculated using an 
      !empirical relationship based on the pore volume THPOR, with 
      !different coefficients for mineral and organic soils:
      !
      !lambda_dry = 0.75*exp(-2.76*THPOR)   (mineral)
      !lambda_dry = 0.30*exp(-2.0*THPOR)    (organic)
      !
      !The saturated thermal conductivity lambda_sat is calculated by 
      !Cote and Konrad as a geometric mean of the conductivities of the 
      !soil components. However, other researchers (e.g. Zhang et al., 
      !2008) have found the linear averaging used by de Vries (1963) to 
      !be more generally accurate:
      !
      !lambda_sat = lambda_w*THPOR + lambda_s*(1 - THPOR)   (unfrozen)
      !lambda_sat = lambda_i*THPOR + lambda_s*(1 - THPOR)   (frozen)
      !
      !where lambda_w is the thermal conductivity of water, lambda_i is 
      !that of ice and lambda_s is that of the soil particles.
      !
      !In the 500 loop, thermal conductivities are calculated for the 
      !top and bottom of each soil layer. The degree of saturation 
      !SATRAT is calculated as the sum of liquid and frozen water 
      !contents, theta_w and THICE, divided by the pore volume. In 
      !organic soils, if the liquid water content of the soil layer is 
      !above the retention capacity THLRET, theta_w at the top of the 
      !soil layer is assumed to be equal to theta_re and Sr at the 
      !bottom of the layer is assumed to be 1. The relative liquid and 
      !frozen water contents, THLSAT and THISAT, are calculated from 
      !theta_w and THICE normalized by theta_w + THICE. The dry thermal 
      !conductivity, and the saturated thermal conductivity for unfrozen 
      !and frozen conditions, are evaluated using the equations above. 
      !The unfrozen and frozen relative thermal conductivity, TCRATU and 
      !TCRATF, are obtained from SATRAT and the appropriate values of 
      !the empirical coefficient kappa. For mineral soils, kappa is 
      !obtained as a weighted average over the percent sand content 
      !(ISAND converted to a real value) and the percentage of fine 
      !material (assumed to be 100-ISAND). The unfrozen and frozen soil 
      !thermal conductivities, TCSOLU and TCSOLF, are then calculated 
      !from TCRATU, TCRATF, and the dry and saturated thermal 
      !conductivities; and the actual thermal conductivity of the soil, 
      !TCSOIL, is determined as the average of TCSOLU and TCSOLF, 
      !weighted according to the relative liquid and frozen water 
      !contents THLSAT and THISAT. If the permeable thickness of the 
      !layer, DELZW, is greater than zero, the thermal conductivity at 
      !the top of the layer is set to TCSOIL; otherwise it is set to the 
      !rock value, TCSAND. If DELZW is less than the thermal thickness 
      !of the layer DELZ, the thermal conductivity at the bottom of the 
      !layer is set to TCSAND; otherwise it is set to TCSOIL. (In the 
      !case of organic soils in the latter branch, if theta_w was 
      !greater than THLRET, the thermal conductivity at the bottom of 
      !the layer is set to the average of the saturated unfrozen and 
      !frozen values, weighted by THLSAT and THISAT.) Finally, if there 
      !is ponded water present on the soil surface, the thermal 
      !conductivity at the top of the first soil layer is treated as
      ! varying linearly from the calculated soil thermal conductivity 
      !if the pond depth ZPOND is zero, to the thermal conductivity of 
      !water if ZPOND >= 10-2 m.
      !                                                                        
      DO 500 J=IG,1,-1
      DO 500 I=IL1,IL2                                                            
          IF    (ISAND(I,1).EQ.-4)                              THEN          
              TCTOPG(I,J)=TCGLAC                                                     
              TCBOTG(I,J)=TCGLAC
          ELSEIF(ISAND(I,J).EQ.-3)                              THEN
              TCTOPC(I,J)=TCSAND
              TCTOPG(I,J)=TCSAND
              TCBOTC(I,J)=TCSAND
              TCBOTG(I,J)=TCSAND
          ELSEIF(ISAND(I,J).EQ.-2)                          THEN 
              IF ((THLIQG(I,J)+THICEG(I,J)).GT.(THLRET(I,J)+0.0001)) 
     1                                                 THEN
                  WTABLE(I)=ZBOTW(I,J)-DELZW(I,J)*MIN(1.0,
     1                      (THLIQG(I,J)+THICEG(I,J)-THLRET(I,J))/
     2                      (THPOR(I,J)-THLRET(I,J)))
              ENDIF
              IF (THLIQG(I,J).GT.(THLRET(I,J)+0.0001)) THEN
                  SATRAT=MIN((THLRET(I,J)+THICEG(I,J))/
     1                   THPOR(I,J), 1.0)              
                  THLSAT=THLIQG(I,J)/(THLIQG(I,J)+THICEG(I,J))          
                  THISAT=THICEG(I,J)/(THLIQG(I,J)+THICEG(I,J))          
                  TCDRY=0.30*EXP(-2.0*THPOR(I,J))
                  TCSATU(I)=TCW*THPOR(I,J)+TCS(I,J)*(1.0-THPOR(I,J))
                  TCSATF(I)=TCICE*THPOR(I,J)+TCS(I,J)*(1.0-THPOR(I,J))
                  TCRATU=0.6*SATRAT/(1.0-0.4*SATRAT)
                  TCRATF=0.25*SATRAT/(1.0-0.75*SATRAT)
                  TCSOLU=(TCSATU(I)-TCDRY)*TCRATU+TCDRY                              
                  TCSOLF=(TCSATF(I)-TCDRY)*TCRATF+TCDRY
                  TCSOIL=TCSOLU*THLSAT+TCSOLF*THISAT
                  IF(DELZW(I,J).GT.0.0) THEN
                      TCTOPC(I,J)=TCSOIL
                      TCTOPG(I,J)=TCSOIL
                  ELSE
                      TCTOPC(I,J)=TCSAND
                      TCTOPG(I,J)=TCSAND
                  ENDIF
                  IF(DELZW(I,J).LT.DELZ(J)) THEN
                      TCBOTC(I,J)=TCSAND
                      TCBOTG(I,J)=TCSAND
                  ELSE
                      TCBOTC(I,J)=TCSATU(I)*THLSAT+TCSATF(I)*THISAT
                      TCBOTG(I,J)=TCSATU(I)*THLSAT+TCSATF(I)*THISAT
                  ENDIF
                  IF(J.EQ.1) THEN
                      TCTOPC(I,J)=TCTOPC(I,J)+(TCW-TCTOPC(I,J))*
     1                            MIN(ZPOND(I),1.0E-2)*100.0
                      TCTOPG(I,J)=TCTOPC(I,J)
                  ENDIF
              ELSE
                  SATRAT=MIN((THLIQG(I,J)+THICEG(I,J))/
     1                   THPOR(I,J), 1.0)              
                  THLSAT=THLIQG(I,J)/(THLIQG(I,J)+THICEG(I,J))          
                  THISAT=THICEG(I,J)/(THLIQG(I,J)+THICEG(I,J))          
                  TCDRY=0.30*EXP(-2.0*THPOR(I,J))
                  TCSATU(I)=TCW*THPOR(I,J)+TCS(I,J)*(1.0-THPOR(I,J))
                  TCSATF(I)=TCICE*THPOR(I,J)+TCS(I,J)*(1.0-THPOR(I,J))
                  TCRATU=0.6*SATRAT/(1.0-0.4*SATRAT)
                  TCRATF=0.25*SATRAT/(1.0-0.75*SATRAT)
                  TCSOLU=(TCSATU(I)-TCDRY)*TCRATU+TCDRY                              
                  TCSOLF=(TCSATF(I)-TCDRY)*TCRATF+TCDRY
                  TCSOIL=TCSOLU*THLSAT+TCSOLF*THISAT
                  IF(DELZW(I,J).GT.0.0) THEN
                      TCTOPC(I,J)=TCSOIL
                      TCTOPG(I,J)=TCSOIL
                  ELSE
                      TCTOPC(I,J)=TCSAND
                      TCTOPG(I,J)=TCSAND
                  ENDIF
                  IF(DELZW(I,J).LT.DELZ(J)) THEN
                      TCBOTC(I,J)=TCSAND
                      TCBOTG(I,J)=TCSAND
                  ELSE
                      TCBOTC(I,J)=TCSOIL
                      TCBOTG(I,J)=TCSOIL
                  ENDIF
                  IF(J.EQ.1) THEN
                      TCTOPC(I,J)=TCTOPC(I,J)+(TCW-TCTOPC(I,J))*
     1                            MIN(ZPOND(I),1.0E-2)*100.0
                      TCTOPG(I,J)=TCTOPC(I,J)
                  ENDIF
              ENDIF    
          ELSE
              SATRAT=MIN((THLIQG(I,J)+THICEG(I,J))/
     1               THPOR(I,J), 1.0)              
              THLSAT=THLIQG(I,J)/(THLIQG(I,J)+THICEG(I,J))          
              THISAT=THICEG(I,J)/(THLIQG(I,J)+THICEG(I,J))          
              TCDRY=0.75*EXP(-2.76*THPOR(I,J))
              TCSATU(I)=TCW*THPOR(I,J)+TCS(I,J)*(1.0-THPOR(I,J))
              TCSATF(I)=TCICE*THPOR(I,J)+TCS(I,J)*(1.0-THPOR(I,J))
              TCKAPU=(4.0*REAL(ISAND(I,J))+1.9*REAL(100-ISAND(I,J)))/
     1               100.0
              TCKAPF=(1.2*REAL(ISAND(I,J))+0.85*REAL(100-ISAND(I,J)))/
     1               100.0
              TCRATU=TCKAPU*SATRAT/(1.0+(TCKAPU-1.0)*SATRAT)
              TCRATF=TCKAPF*SATRAT/(1.0+(TCKAPF-1.0)*SATRAT)
              TCSOLU=(TCSATU(I)-TCDRY)*TCRATU+TCDRY                              
              TCSOLF=(TCSATF(I)-TCDRY)*TCRATF+TCDRY
              TCSOIL=TCSOLU*THLSAT+TCSOLF*THISAT
              IF(DELZW(I,J).GT.0.0) THEN
                  TCTOPC(I,J)=TCSOIL
                  TCTOPG(I,J)=TCSOIL
C                  IF(J.EQ.1) TCTOPC(I,J)=TCTOPC(I,J)*0.1
              ELSE
                  TCTOPC(I,J)=TCSAND
                  TCTOPG(I,J)=TCSAND
              ENDIF
              IF(DELZW(I,J).LT.DELZ(J)) THEN
                  TCBOTC(I,J)=TCSAND
                  TCBOTG(I,J)=TCSAND
              ELSE
                  TCBOTC(I,J)=TCSOIL
                  TCBOTG(I,J)=TCSOIL
              ENDIF
              IF(J.EQ.1) THEN
                  TCTOPC(I,J)=TCTOPC(I,J)+(TCW-TCTOPC(I,J))*
     1                        MIN(ZPOND(I),1.0E-2)*100.0
                  TCTOPG(I,J)=TCTOPC(I,J)
              ENDIF
          ENDIF                                                                   
  500 CONTINUE                                                                    
C                                                                           
C     * ADD PONDED WATER TEMPERATURE TO FIRST SOIL LAYER FOR USE
C     * IN GROUND HEAT FLUX CALCULATIONS.
C
      !
      !Finally, in loop 600, a variable TBAR1P is evaluated, 
      !representing the weighted average value of the first layer soil 
      !temperature and the ponded water, if any. (The heat capacity of 
      !the soil is determined as the weighted average of HCPG over the 
      !permeable thickness DELZW, and the heat capacity of rock, HCPSND, 
      !over the impermeable thickness, DELZ-DELZW.)
      !
      DO 600 I=IL1,IL2
          IF(ZPOND(I).GT.0.)                          THEN 
              TBAR1P(I)=(TPOND(I)*HCPW*ZPOND(I) + 
     1                  TBAR(I,1)*(HCPG(I,1)*DELZW(I,1)+
     2                  HCPSND*(DELZ(1)-DELZW(I,1))))/
     3                  (HCPW*ZPOND(I)+HCPG(I,1)*DELZW(I,1)+
     4                  HCPSND*(DELZ(1)-DELZW(I,1)))
          ELSE
              TBAR1P(I)=TBAR(I,1)
          ENDIF
  600 CONTINUE
C
      RETURN                                                                      
      END 
