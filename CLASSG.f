      SUBROUTINE CLASSG(TBARGAT,THLQGAT,THICGAT,TPNDGAT,ZPNDGAT,
     1                  TBASGAT,ALBSGAT,TSNOGAT,RHOSGAT,SNOGAT, 
     2                  TCANGAT,RCANGAT,SCANGAT,GROGAT, CMAIGAT,
     3                  FCANGAT,LNZ0GAT,ALVCGAT,ALICGAT,PAMXGAT,
     4                  PAMNGAT,CMASGAT,ROOTGAT,RSMNGAT,QA50GAT,
     5                  VPDAGAT,VPDBGAT,PSGAGAT,PSGBGAT,PAIDGAT,
     6                  HGTDGAT,ACVDGAT,ACIDGAT,TSFSGAT,WSNOGAT,
     7                  THPGAT, THRGAT, THMGAT, BIGAT,  PSISGAT,
     8                  GRKSGAT,THRAGAT,HCPSGAT,TCSGAT, IGDRGAT,
     9                  THFCGAT,PSIWGAT,DLZWGAT,ZBTWGAT,VMODGAT,
     A                  ZSNLGAT,ZPLGGAT,ZPLSGAT,TACGAT, QACGAT,
     B                  DRNGAT, XSLPGAT,GRKFGAT,WFSFGAT,WFCIGAT,
     C                  ALGWGAT,ALGDGAT,ASVDGAT,ASIDGAT,AGVDGAT,
     D                  AGIDGAT,ISNDGAT,RADJGAT,ZBLDGAT,Z0ORGAT,
     E                  ZRFMGAT,ZRFHGAT,ZDMGAT, ZDHGAT, FSVHGAT,
     F                  FSIHGAT,CSZGAT, FDLGAT, ULGAT,  VLGAT,  
     G                  TAGAT,  QAGAT,  PRESGAT,PREGAT, PADRGAT,
     H                  VPDGAT, TADPGAT,RHOAGAT,RPCPGAT,TRPCGAT,
     I                  SPCPGAT,TSPCGAT,RHSIGAT,FCLOGAT,DLONGAT,
     J                  GGEOGAT,
     K                  ILMOS,JLMOS,IWMOS,JWMOS,
     L                  NML,NL,NM,ILG,IG,IC,ICP1,
     M                  TBARROT,THLQROT,THICROT,TPNDROT,ZPNDROT,
     N                  TBASROT,ALBSROT,TSNOROT,RHOSROT,SNOROT, 
     O                  TCANROT,RCANROT,SCANROT,GROROT, CMAIROT,
     P                  FCANROT,LNZ0ROT,ALVCROT,ALICROT,PAMXROT,
     Q                  PAMNROT,CMASROT,ROOTROT,RSMNROT,QA50ROT,
     R                  VPDAROT,VPDBROT,PSGAROT,PSGBROT,PAIDROT,
     S                  HGTDROT,ACVDROT,ACIDROT,TSFSROT,WSNOROT,
     T                  THPROT, THRROT, THMROT, BIROT,  PSISROT,
     U                  GRKSROT,THRAROT,HCPSROT,TCSROT, IGDRROT,
     V                  THFCROT,PSIWROT,DLZWROT,ZBTWROT,VMODL,
     W                  ZSNLROT,ZPLGROT,ZPLSROT,TACROT, QACROT,
     X                  DRNROT, XSLPROT,GRKFROT,WFSFROT,WFCIROT,
     Y                  ALGWROT,ALGDROT,ASVDROT,ASIDROT,AGVDROT,
     Z                  AGIDROT,ISNDROT,RADJ   ,ZBLDROW,Z0ORROW,
     +                  ZRFMROW,ZRFHROW,ZDMROW, ZDHROW, FSVHROW,
     +                  FSIHROW,CSZROW, FDLROW, ULROW,  VLROW,  
     +                  TAROW,  QAROW,  PRESROW,PREROW, PADRROW,
     +                  VPDROW, TADPROW,RHOAROW,RPCPROW,TRPCROW,
     +                  SPCPROW,TSPCROW,RHSIROW,FCLOROW,DLONROW,
     +                  GGEOROW  )
C
C     Purpose: Gather variables from two-dimensional arrays (latitude 
C     circle x mosaic tiles) onto long vectors for optimum processing 
C     efficiency on vector supercomputers.
C
C     * OCT 18/11 - M.LAZARE.  ADD IGDR.
C     * OCT 07/11 - M.LAZARE.  ADD VMODL->VMODGAT.
C     * OCT 05/11 - M.LAZARE.  PUT BACK IN PRESGROW->PRESGAT
C     *                        REQUIRED FOR ADDED SURFACE RH 
C     *                        CALCULATION.
C     * OCT 03/11 - M.LAZARE.  REMOVE ALL INITIALIZATION TO
C     *                        ZERO OF GAT ARRAYS (NOW DONE
C     *                        IN CLASS DRIVER).
C     * SEP 16/11 - M.LAZARE.  - ROW->ROT AND GRD->ROW.
C     *                        - REMOVE INITIALIZATION OF
C     *                          {ALVS,ALIR} TO ZERO.
C     *                        - REMOVE PRESGROW->PRESGAT 
C     *                          (OCEAN-ONLY NOW).
C     *                        - RADJROW (64-BIT) NOW RADJ
C     *                          (32-BIT).
C     * MAR 23/06 - D.VERSEGHY. ADD WSNO,FSNO,GGEO.
C     * MAR 18/05 - D.VERSEGHY. ADDITIONAL VARIABLES.
C     * FEB 18/05 - D.VERSEGHY. ADD "TSFS" VARIABLES.
C     * NOV 03/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
C     * AUG 15/02 - D.VERSEGHY. GATHER OPERATION ON CLASS 
C     *                         VARIABLES.
C 
      IMPLICIT NONE
C
C     (Suffix GAT refers to variables on gathered long vectors; suffix 
C     ROT refers to variables on original two-dimensional arrays.)
C
C     * INTEGER CONSTANTS.
C
      INTEGER  NML,NL,NM,ILG,IG,IC,ICP1,K,L,M
C
C     * LAND SURFACE PROGNOSTIC VARIABLES.
C
      REAL TBARROT(NL,NM,IG)    !Temperature of soil layers [K]
      REAL THLQROT(NL,NM,IG)    !Volumetric liquid water content of soil 
                                !layers [m3 m-3]
      REAL THICROT(NL,NM,IG)    !Frozen water content of soil layers 
                                !under vegetation [m3 m-3]
      REAL TPNDROT(NL,NM)   !Temperature of ponded water [K]    
      REAL ZPNDROT(NL,NM)   !Depth of ponded water on surface [m]
      REAL TBASROT(NL,NM)   !Temperature of bedrock in third soil layer 
                            ![K]
      REAL ALBSROT(NL,NM)   !Snow albedo [ ]
      REAL TSNOROT(NL,NM)   !Snowpack temperature [K]
      REAL RHOSROT(NL,NM)   !Density of snow [kg m-3]
      REAL SNOROT (NL,NM)   !Mass of snow pack [kg m-2] 
      REAL TCANROT(NL,NM)   !Vegetation canopy temperature [K] 
      REAL RCANROT(NL,NM)   !Intercepted liquid water stored on canopy 
                            ![kg m-2]
      REAL SCANROT(NL,NM)   !Intercepted frozen water stored on canopy 
                            ![kg m-2]
      REAL GROROT (NL,NM)   !Vegetation growth index [ ] 
      REAL CMAIROT(NL,NM)   !Aggregated mass of vegetation canopy 
                            ![kg m-2]
      REAL TSFSROT(NL,NM,4) !Ground surface temperature over subarea [K] 
      REAL TACROT (NL,NM)   !Temperature of air within vegetation canopy 
                            ![K]
      REAL QACROT (NL,NM)   !Specific humidity of air within vegetation 
                            !canopy [kg kg-1]
      REAL WSNOROT(NL,NM)   !Liquid water content of snow pack [kg m-2]
C
      REAL    TBARGAT(ILG,IG),   THLQGAT(ILG,IG),   THICGAT(ILG,IG), 
     1        TPNDGAT(ILG),      ZPNDGAT(ILG),      TBASGAT(ILG),   
     2        ALBSGAT(ILG),      TSNOGAT(ILG),      RHOSGAT(ILG),   
     3        SNOGAT (ILG),      TCANGAT(ILG),      RCANGAT(ILG),   
     4        SCANGAT(ILG),      GROGAT (ILG),      CMAIGAT(ILG),
     5        TSFSGAT(ILG,4),    TACGAT (ILG),      QACGAT (ILG),
     6        WSNOGAT(ILG)
C
C     * GATHER-SCATTER INDEX ARRAYS.
C
      INTEGER  ILMOS (ILG)  !Index of latitude grid cell corresponding 
                            !to current element of gathered vector of 
                            !land surface variables [ ]
      INTEGER  JLMOS (ILG)  !Index of mosaic tile corresponding to 
                            !current element of gathered vector of land 
                            !surface variables [ ]
      INTEGER  IWMOS (ILG)  !Index of latitude grid cell corresponding 
                            !to current element of gathered vector of 
                            !inland water body variables [ ]
      INTEGER  JWMOS (ILG)  !Index of mosaic tile corresponding to 
                            !current element of gathered vector of 
                            !inland water body variables [ ]

C
C     * CANOPY AND SOIL INFORMATION ARRAYS.
C     * (THE LENGTH OF THESE ARRAYS IS DETERMINED BY THE NUMBER
C     * OF SOIL LAYERS (3) AND THE NUMBER OF BROAD VEGETATION
C     * CATEGORIES (4, OR 5 INCLUDING URBAN AREAS).)
C
      REAL FCANROT(NL,NM,ICP1)  !Maximum fractional coverage of modelled 
                                !area by vegetation category [ ]
      REAL LNZ0ROT(NL,NM,ICP1)  !Natural logarithm of maximum roughness 
                                !length of vegetation category [ ]
      REAL ALVCROT(NL,NM,ICP1)  !Background average visible albedo of 
                                !vegetation category [ ]
      REAL ALICROT(NL,NM,ICP1)  !Background average near-infrared albedo 
                                !of vegetation category [ ]
      REAL PAMXROT(NL,NM,IC)    !Maximum plant area index of vegetation 
                                !category [ ]
      REAL PAMNROT(NL,NM,IC)    !Minimum plant area index of vegetation 
                                !category [ ]
      REAL CMASROT(NL,NM,IC)    !Maximum canopy mass for vegetation 
                                !category [kg m-2]
      REAL ROOTROT(NL,NM,IC)    !Maximum rooting depth of vegetation 
                                !category [m]
      REAL RSMNROT(NL,NM,IC)    !Minimum stomatal resistance of 
                                !vegetation category [s m-1]
      REAL QA50ROT(NL,NM,IC)    !Reference value of incoming shortwave 
                                !radiation for vegetation category (used 
                                !in stomatal resistance calculation) 
                                ![W m-2]
      REAL VPDAROT(NL,NM,IC)    !Vapour pressure deficit coefficient for 
                                !vegetation category (used in stomatal 
                                !resistance calculation) [ ]
      REAL VPDBROT(NL,NM,IC)    !Vapour pressure deficit coefficient for 
                                !vegetation category (used in stomatal 
                                !resistance calculation) [ ]
      REAL PSGAROT(NL,NM,IC)    !Soil moisture suction coefficient for 
                                !vegetation category (used in stomatal 
                                !resistance calculation) [ ]
      REAL PSGBROT(NL,NM,IC)    !Soil moisture suction coefficient for 
                                !vegetation category (used in stomatal 
                                !resistance calculation) [ ]
      REAL PAIDROT(NL,NM,IC)    !Optional user-specified value of plant 
                                !area indices of vegetation categories 
                                !to override CLASS-calculated values [ ]
      REAL HGTDROT(NL,NM,IC)    !Optional user-specified values of 
                                !height of vegetation categories to 
                                !override CLASS-calculated values [m]
      REAL ACVDROT(NL,NM,IC)    !Optional user-specified value of canopy 
                                !visible albedo to override CLASS-
                                !calculated value [ ]
      REAL ACIDROT(NL,NM,IC)    !Optional user-specified value of canopy 
                                !near-infrared albedo to override CLASS-
                                !calculated value [ ]
C
      REAL          FCANGAT(ILG,ICP1),   LNZ0GAT(ILG,ICP1),
     1              ALVCGAT(ILG,ICP1),   ALICGAT(ILG,ICP1),
     2              PAMXGAT(ILG,IC),     PAMNGAT(ILG,IC),
     3              CMASGAT(ILG,IC),     ROOTGAT(ILG,IC),
     4              RSMNGAT(ILG,IC),     QA50GAT(ILG,IC),
     5              VPDAGAT(ILG,IC),     VPDBGAT(ILG,IC),
     6              PSGAGAT(ILG,IC),     PSGBGAT(ILG,IC),
     7              PAIDGAT(ILG,IC),     HGTDGAT(ILG,IC),
     8              ACVDGAT(ILG,IC),     ACIDGAT(ILG,IC)
C
      REAL THPROT (NL,NM,IG)    !Pore volume in soil layer [m3 m-3] 
      REAL THRROT (NL,NM,IG)    !Liquid water retention capacity for 
                                !organic soil [m3 m-3 ]
      REAL THMROT (NL,NM,IG)    !Residual soil liquid water content 
                                !remaining after freezing or evaporation 
                                ![m3 m-3]
      REAL BIROT  (NL,NM,IG)    !Clapp and Hornberger empirical “b” 
                                !parameter [ ]
      REAL PSISROT(NL,NM,IG)    !Soil moisture suction at saturation [m]
      REAL GRKSROT(NL,NM,IG)    !Saturated hydraulic conductivity of 
                                !soil layers [m s-1]
      REAL THRAROT(NL,NM,IG)    !Fractional saturation of soil behind 
                                !the wetting front [ ]
      REAL HCPSROT(NL,NM,IG)    !Volumetric heat capacity of soil 
                                !particles [J m-3]
      REAL TCSROT (NL,NM,IG)    !Thermal conductivity of soil particles 
                                ![W m-1 K-1]
      REAL THFCROT(NL,NM,IG)    !Field capacity [m3 m-3]
      REAL PSIWROT(NL,NM,IG)    !Soil moisture suction at wilting point 
                                ![m]
      REAL DLZWROT(NL,NM,IG)    !Permeable thickness of soil layer [m]
      REAL ZBTWROT(NL,NM,IG)    !Depth to permeable bottom of soil layer 
                                ![m]
      REAL DRNROT (NL,NM)       !Drainage index at bottom of soil 
                                !profile [ ]
      REAL XSLPROT(NL,NM)       !Surface slope (used when running MESH 
                                !code) [degrees]
      REAL GRKFROT(NL,NM)       !WATROF parameter used when running MESH 
                                !code [ ]
      REAL WFSFROT(NL,NM)       !WATROF parameter used when running MESH 
                                !code [ ]
      REAL WFCIROT(NL,NM)       !WATROF parameter used when running MESH 
                                !code [ ]
      REAL ALGWROT(NL,NM)       !Reference albedo for saturated soil [ ]
      REAL ALGDROT(NL,NM)       !Reference albedo for dry soil [ ]
      REAL ASVDROT(NL,NM)       !Optional user-specified value of snow 
                                !visible albedo to override CLASS-
                                !calculated value [ ]
      REAL ASIDROT(NL,NM)       !Optional user-specified value of snow 
                                !near-infrared albedo to override CLASS-
                                !calculated value [ ]
      REAL AGVDROT(NL,NM)       !Optional user-specified value of ground 
                                !visible albedo to override CLASS-
                                !calculated value [ ]
      REAL AGIDROT(NL,NM)       !Optional user-specified value of ground 
                                !near-infrared albedo to override CLASS-
                                !calculated value [ ]
      REAL ZSNLROT(NL,NM)       !Limiting snow depth below which 
                                !coverage is < 100% [m]
      REAL ZPLGROT(NL,NM)       !Maximum water ponding depth for snow-
                                !free subareas (user-specified when 
                                !running MESH code) [m]
      REAL ZPLSROT(NL,NM)       !Maximum water ponding depth for snow-
                                !covered subareas (user-specified when 
                                !running MESH code) [m]
C

      REAL    THPGAT (ILG,IG),   THRGAT (ILG,IG),   THMGAT (ILG,IG),
     1        BIGAT  (ILG,IG),   PSISGAT(ILG,IG),   GRKSGAT(ILG,IG),   
     2        THRAGAT(ILG,IG),   HCPSGAT(ILG,IG), 
     3        TCSGAT (ILG,IG),   THFCGAT(ILG,IG),   PSIWGAT(ILG,IG),  
     4        DLZWGAT(ILG,IG),   ZBTWGAT(ILG,IG),   
     5        DRNGAT (ILG),      XSLPGAT(ILG),      GRKFGAT(ILG),
     6        WFSFGAT(ILG),      WFCIGAT(ILG),      ALGWGAT(ILG),     
     7        ALGDGAT(ILG),      ASVDGAT(ILG),      ASIDGAT(ILG),     
     8        AGVDGAT(ILG),      AGIDGAT(ILG),      ZSNLGAT(ILG),
     9        ZPLGGAT(ILG),      ZPLSGAT(ILG)
C
      INTEGER ISNDROT(NL,NM,IG), ISNDGAT(ILG,IG)    !Sand content flag
      INTEGER IGDRROT(NL,NM),    IGDRGAT(ILG)   !Index of soil layer in 
                                                !which bedrock is 
                                                !encountered
C
C     * ATMOSPHERIC AND GRID-CONSTANT INPUT VARIABLES.
C
      REAL ZRFMROW( NL) !Reference height associated with forcing wind 
                        !speed [m]
      REAL ZRFHROW( NL) !Reference height associated with forcing air 
                        !temperature and humidity [m]
      REAL ZDMROW ( NL) !User-specified height associated with diagnosed 
                        !anemometer-level wind speed [m]
      REAL ZDHROW ( NL) !User-specified height associated with diagnosed 
                        !screen-level variables [m]
      REAL FSVHROW( NL) !Visible radiation incident on horizontal 
                        !surface [W m-2]
      REAL FSIHROW( NL) !Near-infrared radiation incident on horizontal 
                        !surface [W m -2]
      REAL CSZROW ( NL) !Cosine of solar zenith angle [ ]
      REAL FDLROW ( NL) !Downwelling longwave radiation at bottom of 
                        !atmosphere [W m-2]
      REAL ULROW  ( NL) !Zonal component of wind speed [m s-1]
      REAL VLROW  ( NL) !Meridional component of wind speed [m s-1]
      REAL TAROW  ( NL) !Air temperature at reference height [K]
      REAL QAROW  ( NL) !Specific humidity at reference height [kg kg-1]
      REAL PRESROW( NL) !Surface air pressure [Pa]
      REAL PREROW ( NL) !Surface precipitation rate [kg m-2 s-1]
      REAL PADRROW( NL) !Partial pressure of dry air [Pa]
      REAL VPDROW ( NL) !Vapour pressure deficit [mb]
      REAL TADPROW( NL) !Dew point temperature of air [K]
      REAL RHOAROW( NL) !Density of air [kg m-3]
      REAL ZBLDROW( NL) !Atmospheric blending height for surface 
                        !roughness length averaging [m]
      REAL Z0ORROW( NL) !Orographic roughness length [m]
      REAL RPCPROW( NL) !Rainfall rate over modelled area [m s-1]
      REAL TRPCROW( NL) !Rainfall temperature [K]
      REAL SPCPROW( NL) !Snowfall rate over modelled area [m s-1]
      REAL TSPCROW( NL) !Snowfall temperature [K]
      REAL RHSIROW( NL) !Density of fresh snow [kg m-3]
      REAL FCLOROW( NL) !Fractional cloud cover [ ]
      REAL DLONROW( NL) !Longitude of grid cell (east of Greenwich) 
                        ![degrees]
      REAL GGEOROW( NL) !Geothermal heat flux at bottom of soil profile 
                        ![W m-2]
      REAL RADJ   ( NL) !Latitude of grid cell (positive north of 
                        !equator) [rad]
      REAL VMODL  ( NL) !Wind speed at reference height [m s-1]
C
      REAL  ZRFMGAT(ILG), ZRFHGAT(ILG), ZDMGAT (ILG), ZDHGAT (ILG),
     1      FSVHGAT(ILG), FSIHGAT(ILG), CSZGAT (ILG), FDLGAT (ILG), 
     2      ULGAT  (ILG), VLGAT  (ILG), TAGAT  (ILG), QAGAT  (ILG), 
     3      PRESGAT(ILG), PREGAT (ILG), PADRGAT(ILG), VPDGAT (ILG), 
     4      TADPGAT(ILG), RHOAGAT(ILG), ZBLDGAT(ILG), Z0ORGAT(ILG),
     5      RPCPGAT(ILG), TRPCGAT(ILG), SPCPGAT(ILG), TSPCGAT(ILG),
     6      RHSIGAT(ILG), FCLOGAT(ILG), DLONGAT(ILG), GGEOGAT(ILG),
     7      RADJGAT(ILG), VMODGAT(ILG)
C----------------------------------------------------------------------
      !
      !The prognostic, background and input variables are gathered into 
      !long arrays (collapsing the latitude and mosaic dimensions into 
      !one, but retaining the soil level and canopy category dimensions) 
      !using the pointer vectors generated in GATPREP.
      !
      DO 100 K=1,NML
          TPNDGAT(K)=TPNDROT(ILMOS(K),JLMOS(K))  
          ZPNDGAT(K)=ZPNDROT(ILMOS(K),JLMOS(K))  
          TBASGAT(K)=TBASROT(ILMOS(K),JLMOS(K))  
          ALBSGAT(K)=ALBSROT(ILMOS(K),JLMOS(K))  
          TSNOGAT(K)=TSNOROT(ILMOS(K),JLMOS(K))  
          RHOSGAT(K)=RHOSROT(ILMOS(K),JLMOS(K))  
          SNOGAT (K)=SNOROT (ILMOS(K),JLMOS(K))  
          WSNOGAT(K)=WSNOROT(ILMOS(K),JLMOS(K))  
          TCANGAT(K)=TCANROT(ILMOS(K),JLMOS(K))  
          RCANGAT(K)=RCANROT(ILMOS(K),JLMOS(K))  
          SCANGAT(K)=SCANROT(ILMOS(K),JLMOS(K))  
          GROGAT (K)=GROROT (ILMOS(K),JLMOS(K))  
          CMAIGAT(K)=CMAIROT(ILMOS(K),JLMOS(K))  
          DRNGAT (K)=DRNROT (ILMOS(K),JLMOS(K))  
c         XSLPGAT(K)=XSLPROT(ILMOS(K),JLMOS(K))  
c         GRKFGAT(K)=GRKFROT(ILMOS(K),JLMOS(K))  
c         WFSFGAT(K)=WFSFROT(ILMOS(K),JLMOS(K))  
c         WFCIGAT(K)=WFCIROT(ILMOS(K),JLMOS(K))  
          ALGWGAT(K)=ALGWROT(ILMOS(K),JLMOS(K))  
          ALGDGAT(K)=ALGDROT(ILMOS(K),JLMOS(K))  
c         ASVDGAT(K)=ASVDROT(ILMOS(K),JLMOS(K))  
c         ASIDGAT(K)=ASIDROT(ILMOS(K),JLMOS(K))  
c         AGVDGAT(K)=AGVDROT(ILMOS(K),JLMOS(K))  
c         AGIDGAT(K)=AGIDROT(ILMOS(K),JLMOS(K))  
          ZSNLGAT(K)=ZSNLROT(ILMOS(K),JLMOS(K))  
c         ZPLGGAT(K)=ZPLGROT(ILMOS(K),JLMOS(K))  
c         ZPLSGAT(K)=ZPLSROT(ILMOS(K),JLMOS(K))  
          TACGAT (K)=TACROT (ILMOS(K),JLMOS(K))  
          QACGAT (K)=QACROT (ILMOS(K),JLMOS(K))  
          IGDRGAT(K)=IGDRROT(ILMOS(K),JLMOS(K))
          ZBLDGAT(K)=ZBLDROW(ILMOS(K))
          Z0ORGAT(K)=Z0ORROW(ILMOS(K))
          ZRFMGAT(K)=ZRFMROW(ILMOS(K))
          ZRFHGAT(K)=ZRFHROW(ILMOS(K))
          ZDMGAT (K)=ZDMROW(ILMOS(K))
          ZDHGAT (K)=ZDHROW(ILMOS(K))
          FSVHGAT(K)=FSVHROW(ILMOS(K))
          FSIHGAT(K)=FSIHROW(ILMOS(K))
          CSZGAT (K)=CSZROW (ILMOS(K))
          FDLGAT (K)=FDLROW (ILMOS(K))
          ULGAT  (K)=ULROW  (ILMOS(K))
          VLGAT  (K)=VLROW  (ILMOS(K))
          TAGAT  (K)=TAROW  (ILMOS(K))
          QAGAT  (K)=QAROW  (ILMOS(K))
          PRESGAT(K)=PRESROW(ILMOS(K))
          PREGAT (K)=PREROW (ILMOS(K))
          PADRGAT(K)=PADRROW(ILMOS(K))
          VPDGAT (K)=VPDROW (ILMOS(K))
          TADPGAT(K)=TADPROW(ILMOS(K))
          RHOAGAT(K)=RHOAROW(ILMOS(K))
          RPCPGAT(K)=RPCPROW(ILMOS(K))
          TRPCGAT(K)=TRPCROW(ILMOS(K))
          SPCPGAT(K)=SPCPROW(ILMOS(K))
          TSPCGAT(K)=TSPCROW(ILMOS(K))
          RHSIGAT(K)=RHSIROW(ILMOS(K))
          FCLOGAT(K)=FCLOROW(ILMOS(K))
          DLONGAT(K)=DLONROW(ILMOS(K))
          GGEOGAT(K)=GGEOROW(ILMOS(K))
          RADJGAT(K)=RADJ   (ILMOS(K))
          VMODGAT(K)=VMODL  (ILMOS(K))
  100 CONTINUE
C
      DO 250 L=1,IG
      DO 200 K=1,NML
          TBARGAT(K,L)=TBARROT(ILMOS(K),JLMOS(K),L)
          THLQGAT(K,L)=THLQROT(ILMOS(K),JLMOS(K),L)
          THICGAT(K,L)=THICROT(ILMOS(K),JLMOS(K),L)
          THPGAT (K,L)=THPROT (ILMOS(K),JLMOS(K),L)
          THRGAT (K,L)=THRROT (ILMOS(K),JLMOS(K),L)
          THMGAT (K,L)=THMROT (ILMOS(K),JLMOS(K),L)
          BIGAT  (K,L)=BIROT  (ILMOS(K),JLMOS(K),L)
          PSISGAT(K,L)=PSISROT(ILMOS(K),JLMOS(K),L)
          GRKSGAT(K,L)=GRKSROT(ILMOS(K),JLMOS(K),L)
          THRAGAT(K,L)=THRAROT(ILMOS(K),JLMOS(K),L)
          HCPSGAT(K,L)=HCPSROT(ILMOS(K),JLMOS(K),L)
          TCSGAT (K,L)=TCSROT (ILMOS(K),JLMOS(K),L)
          THFCGAT(K,L)=THFCROT(ILMOS(K),JLMOS(K),L)
          PSIWGAT(K,L)=PSIWROT(ILMOS(K),JLMOS(K),L)
          DLZWGAT(K,L)=DLZWROT(ILMOS(K),JLMOS(K),L)
          ZBTWGAT(K,L)=ZBTWROT(ILMOS(K),JLMOS(K),L)
          ISNDGAT(K,L)=ISNDROT(ILMOS(K),JLMOS(K),L)
  200 CONTINUE
  250 CONTINUE
C
      DO 300 L=1,ICP1
      DO 300 K=1,NML
          FCANGAT(K,L)=FCANROT(ILMOS(K),JLMOS(K),L)
          LNZ0GAT(K,L)=LNZ0ROT(ILMOS(K),JLMOS(K),L)
          ALVCGAT(K,L)=ALVCROT(ILMOS(K),JLMOS(K),L)
          ALICGAT(K,L)=ALICROT(ILMOS(K),JLMOS(K),L)
  300 CONTINUE
C
      DO 400 L=1,IC
      DO 400 K=1,NML
          PAMXGAT(K,L)=PAMXROT(ILMOS(K),JLMOS(K),L)
          PAMNGAT(K,L)=PAMNROT(ILMOS(K),JLMOS(K),L)
          CMASGAT(K,L)=CMASROT(ILMOS(K),JLMOS(K),L)
          ROOTGAT(K,L)=ROOTROT(ILMOS(K),JLMOS(K),L)
          RSMNGAT(K,L)=RSMNROT(ILMOS(K),JLMOS(K),L)
          QA50GAT(K,L)=QA50ROT(ILMOS(K),JLMOS(K),L)
          VPDAGAT(K,L)=VPDAROT(ILMOS(K),JLMOS(K),L)
          VPDBGAT(K,L)=VPDBROT(ILMOS(K),JLMOS(K),L)
          PSGAGAT(K,L)=PSGAROT(ILMOS(K),JLMOS(K),L)
          PSGBGAT(K,L)=PSGBROT(ILMOS(K),JLMOS(K),L)
c         PAIDGAT(K,L)=PAIDROT(ILMOS(K),JLMOS(K),L)
c         HGTDGAT(K,L)=HGTDROT(ILMOS(K),JLMOS(K),L)
c         ACVDGAT(K,L)=ACVDROT(ILMOS(K),JLMOS(K),L)
c         ACIDGAT(K,L)=ACIDROT(ILMOS(K),JLMOS(K),L)
          TSFSGAT(K,L)=TSFSROT(ILMOS(K),JLMOS(K),L)
400   CONTINUE

      RETURN
      END
