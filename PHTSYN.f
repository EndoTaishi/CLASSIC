      SUBROUTINE PHTSYN(   AILCG,  FCANC, TCAN, CO2CONC,  PRESSG,   FC,
     1                     CFLUX,     QA, QSWV,      IC,   THLIQ, SAND,
     2                        TA ,  CLAY, RMAT,   COSZS, XDIFFUS,  ILG,
     3                       IL1  ,  IL2,   IG,     ICC,   ISNOW, SLAI,
     4                   FCANCMX,  CTEM1, CTEM2,  L2MAX, NOL2PFTS,
C    ---------------------- INPUTS ABOVE, OUTPUTS BELOW ---------------
     5                        RC,  CO2I1, CO2I2, AN_VEG, RML_VEG)
C
C
C               CANADIAN TERRESTRIAL ECOSYSTEM MODEL (CTEM) V1.0
C                       PHOTOSYNTHESIS SUBROUTINE
C
C     * 7 SEP.  2001 - PROGRAM TO CALCULATE STOMATAL CONDUCTANCE, TO BE
C     * V. ARORA       USED BY CLASS, BY ESTIMATING PHOTOSYNTHESIS AND THEN
C     *                RELATING PHOTOSYNTHESIS TO STOMATAL CONDUCTANCE.
C
C     1.  SINGLE-LEAF & TWO-LEAF COMBINED VERSION, CAN USE EITHER APPROACH
C     2.  CAN USE EITHER BWB OR LEUNING TYPE STOMATAL CONDUCTANCE FORMULATION
C     3.  ALSO, CAN USE SMOOTHED AVERAGE OF THE 3 LIMITING RATES, MIN. OF
C         THE 3 LIMITING RATES, OR MIN. OF LIGHT AND RUBSICO RATES.
C
C     CLASS' 4 MAJOR VEGETATION TYPES ARE
C
C     1. NEEDLE LEAF OR TALL CONIFEROUS (C3, DECIDUOUS AND EVERGREEN)
C     2. BROAD LEAF (C3, DECD. AND EVRG.)
C     3. ARABLE & CROPS - (BOTH C3 AND C4)
C     4. GRASSES, TUNDRA, ETC. (BOTH C3 AND C4)
C
C     BUT FOR PHOTOSYNTHESIS WE NEED TO MAKE DISTINCTION BETWEEN C3 AND
C     C4, AND DECIDUOUS AND EVERGREEN. SO THESE 4 VEGETATION TYPES GET
C     CONVERTED INTO THE FOLLOWING 9
C
C     1. NEEDLE LEAF EVERGREEN, C3
C     2. NEEDLE LEAF DECIDUOUS, C3
C     3. BROAD LEAF EVERGREEN, C3
C     4. BROAD LEAF COLD DECIDUOUS, C3
C     5. BROAD LEAF DRY DECIDUOUS, C3
C     6. C3 CROP
C     7. C4 CROP
C     8. C3 GREEN GRASS
C     9. C4 GREEN GRASS
C
C     INPUTS
C
C     AILCG     - GREEN LEAF AREA INDEX FOR USE BY PHOTOSYNTHESIS, M^2/M^2
C     FCANC     - FRACTIONAL COVERAGE OF CTEM's 9 PFTs
C     TCAN      - CANOPTY TEMPERATURE, KELVIN
C     CO2CONC   - ATMOS. CO2 IN PPM, AND THEN CONVERT IT TO
C                 PARTIAL PRESSURE, PASCALS, CO2A, FOR USE IN THIS SUBROUTINE
C     PRESSG    - ATMOS. PRESSURE, PASCALS
C     FC        - SUM OF ALL FCANC OVER A GIVEN SUB-AREA
C     CFLUX     - AERODYNAMIC CONDUCTANCE, M/S
C     QA        - SCREEN LEVEL HUMIDITY IN KG/KG
C     THLIQ     - LIQUID MOIS. CONTENT OF 3 SOIL LAYERS
C     SAND      - PERCENTAGE SAND
C     CLAY      - PERCENTAGE CLAY
C     TA        - AIR TEMPERATURE IN KELVINS
C     RMAT      - FRACTION OF ROOTS IN EACH LAYER (grid cell, vegetation, layer)
C     COSZS     - COS OF ZENITH ANGLE
C     XDIFFUS   - FRACTION OF DIFFUSED PAR
C     QSWV      - ABSORBED VISIBLE PART OF SHORTWAVE RADIATION, W/M^2
C     IC        - NO. OF CLASS VEGETATION TYPES, 4
C     ICC       - NO. OF CTEM's PFTs, CURRENTLY 9
C     IG        - NO. OF SOIL LAYERS, 3
C     ILG       - NO. OF GRID CELLS IN LATITUDE CIRCLE
C     IL1, IL2  - IL1=1, IL2=ILG
C     ISNOW     - INTEGER (0 or 1) TELLING IF PHTSYN IS TO BE RUN OVER
C                 CANOPY OVER SNOW OR CANOPY OVER GROUND SUBAREA
C                 (RLi: ISNOW=0 -->GROUND; ISNOW=1 -->SNOW)
C     SLAI      - STORAGE LAI. THIS LAI IS USED FOR PHTSYN EVEN IF ACTUAL
C                 LAI IS ZERO. ESTIMATE OF NET PHOTOSYNTHESIS BASED ON
C                 SLAI IS USED FOR INITIATING LEAF ONSET. SEE PHENOLGY
C                 SUBROUTINE FOR MORE DETAILS.
C     FCANCMX   - MAX. FRACTIONAL COVERAGES OF CTEM's 8 PFTs. THIS IS
C                 DIFFERENT FROM FCANC AND FCANCS (WHICH MAY VARY WITH
C                 SNOW DEPTH). FCANCMX DOESN'T CHANGE, UNLESS OF
C                 COURSE ITS CHANGED BY LAND USE CHANGE OR DYNAMIC VEGETATION.
C     CTEM1     - LOGICAL BOOLEAN FOR USING CTEM's STOMATAL RESISTANCE OR NOT
C     CTEM2     - LOGICAL BOOLEAN FOR USING CTEM's STRUCTURAL ATTRIBUTES
C                 OR NOT.
C     L2MAX     - MAX. NUMBER OF LEVEL 2 PFTs
C     NOL2MAX   - NUMBER OF LEVEL 2 CTEM PFTs
C
C     OUTPUTS
C
C     RC        - GRID-AVERAGED STOMATAL RESISTANCE, S/M
C     CO2I1     - INTERCELLULAR CO2 CONCENTRATION FROM THE PREVIOUS TIME
C                 STEP WHICH GETS UPDATED FOR THE SINGLE LEAF OR THE SUNLIT
C                 PART OF THE TWO LEAF MODEL
C     CO2I2     - INTERCELLULAR CO2 CONCENTRATION FOR THE SHADED PART OF
C                 THE TWO LEAF MODEL FROM THE PREVIOUS TIME STEP
C     AN_VEG    - NET PHOTOSYNTHESIS RATE, u MOL CO2/M2.S FOR EACH PFT
C     RML_VEG   - LEAF RESPIRATION RATE, u MOL CO2/M2.S FOR EACH PFT
C
C     ------------------------------------------------------------------
C
      IMPLICIT NONE
C
      INTEGER KK
      PARAMETER (KK=12)  ! PRODUCT OF CLASS PFTs AND L2MAX (4 x 3 = 12)
C
      INTEGER       ILG,        IC,          I,         J,        IL1,
     1              IL2,  IT_COUNT,    REQITER,        IG,        ICC,
     2          LEAFOPT,   PS_COUP,   ISC4(KK),     ISNOW(ILG),   K1C,   
     3 USESLAI(ILG,ICC),        K2C,     ICOUNT, SORT(ICC), USEBB(ICC),
     4     NOL2PFTS(IC),         M,      L2MAX,        N, ISAND(ILG,IG)
C
      REAL FCANC(ILG,ICC),  AILCG(ILG,ICC),    TCAN(ILG),      FC(ILG),
     1         CFLUX(ILG),    SAND(ILG,IG), CLAY(ILG,IG), FC_TEST(ILG),
     2  USEAILCG(ILG,ICC),   SLAI(ILG,ICC),      QA(ILG),  INICO2I(KK),
     3         TTEMP(ILG),         RH(ILG)
C
      REAL   CO2CONC(ILG),       ALPHA(KK),    OMEGA(KK),  SMSCALE(KK),
     1              TFREZ,       STD_PRESS,
     2    SM_FUNC(ILG,IG),         AVE_SM_FUNC(ILG,ICC)
C
      REAL       VMAX(KK),          KN(KK),     TUP(KK),     Q10_FUNCN,
     1           TLOW(KK),  VMAXC(ILG,ICC),    Q10_FUNC,  JE3(ILG,ICC),
     2    VMUNS1(ILG,ICC), VMUNS2(ILG,ICC),            VMUNS3(ILG,ICC),
     3     VMUNS(ILG,ICC),     VM(ILG,ICC),  SIGMA(ILG),   TGAMMA(ILG),
     4            KC(ILG),        KO (ILG),              CO2I(ILG,ICC),
     5 PREV_CO2I(ILG,ICC),   FPAR(ILG,ICC), JC(ILG,ICC),  JC1(ILG,ICC),
     6       JC2(ILG,ICC),    JC3(ILG,ICC), JE(ILG,ICC),  JE1(ILG,ICC),
     7       JE2(ILG,ICC),       IPAR(ILG), JS(ILG,ICC),   PRESSG(ILG),
     8     A_VEG(ILG,ICC),RML_VEG(ILG,ICC),            AN_VEG(ILG,ICC),
     9            GB(ILG), RC_VEG(ILG,ICC),   QSWV(ILG), GCTU(ILG,ICC),
     A            TA(ILG),  GCMIN(ILG,ICC),           RMAT(ILG,ICC,IG),
     B     WILTSM(ILG,IG), FIELDSM(ILG,IG),             GCMAX(ILG,ICC),
     C     CO2I1(ILG,ICC),  CO2I2(ILG,ICC),  EASAT(ILG),       EA(ILG),
     D           VPD(ILG),              CA,          CB, THLIQ(ILG,IG),
     E  VPD_TERM(ILG,ICC), FCANCMX(ILG,ICC),  Q10_FUNCD
C
      REAL         MM(KK),          BB(KK),    VPD0(KK),           Q10,
     1     CO2LS(ILG,ICC),     GC(ILG,ICC),     RC(ILG),        SN(KK),
     3       O2_CONC(ILG),       CO2A(ILG),     CO2IMAX,    COSZS(ILG),
     4       XDIFFUS(ILG),            ZERO,               RMLCOEFF(KK)
C
      REAL         TEMP_B,          TEMP_C,      TEMP_R,       TEMP_Q1,
     1            TEMP_Q2,         TEMP_JP,       BETA1,         BETA2,
     2            TEMP_AN
C
C     FOR LIMITING CO2 UPTAKE
C
      REAL DELTA_CO2(ILG),   N_EFFECT(ILG),      GAMMA_W,       GAMMA_M
C
      LOGICAL SMOOTH, MIN2, MIN3, CTEM1, CTEM2
C
C     SOIL PROPERTY ARRAYS.
C
      REAL  THPOR(ILG,IG),  PSISAT(ILG,IG),   B(ILG,IG),  GRKSAT(ILG,IG)
C
C     -----------------------------------------------------------------
C                 VARIABLES USED ONLY FOR THE TWO LEAF MODEL
C
      REAL           CHI(KK),       GDIR(ILG,ICC),          KB(ILG,ICC),
     2     FPAR_SUN(ILG,ICC),   FPAR_SHA(ILG,ICC),   VMAXC_SUN(ILG,ICC),
     3    VMAXC_SHA(ILG,ICC), VMUNS1_SUN(ILG,ICC),  VMUNS1_SHA(ILG,ICC),
     4    VMUNS_SUN(ILG,ICC),  VMUNS_SHA(ILG,ICC),      VM_SUN(ILG,ICC),
     5       VM_SHA(ILG,ICC), CO2I_SUN(ILG,ICC), PREV_CO2I_SUN(ILG,ICC),
     6     CO2I_SHA(ILG,ICC),                    PREV_CO2I_SHA(ILG,ICC),
     7      JC1_SUN(ILG,ICC),    JC1_SHA(ILG,ICC),     JC3_SUN(ILG,ICC),
     8      JC3_SHA(ILG,ICC),     JC_SUN(ILG,ICC),      JC_SHA(ILG,ICC),
     9         IPAR_SUN(ILG),       IPAR_SHA(ILG),     JE1_SUN(ILG,ICC),
     A      JE1_SHA(ILG,ICC),    JE2_SUN(ILG,ICC),     JE2_SHA(ILG,ICC),
     B       JE_SUN(ILG,ICC),     JE_SHA(ILG,ICC),      JS_SUN(ILG,ICC),
     C       JS_SHA(ILG,ICC),  A_VEG_SUN(ILG,ICC),   A_VEG_SHA(ILG,ICC),
     D      RML_SUN(ILG,ICC),    RML_SHA(ILG,ICC),      AN_SUN(ILG,ICC),
     E       AN_SHA(ILG,ICC),  CO2LS_SUN(ILG,ICC),   CO2LS_SHA(ILG,ICC),
     F    AILCG_SUN(ILG,ICC),  AILCG_SHA(ILG,ICC),      GC_SUN(ILG,ICC),
     G       GC_SHA(ILG,ICC),   GCTU_SUN(ILG,ICC),    GCTU_SHA(ILG,ICC)
C
      REAL TEMP_PHI1, TEMP_PHI2
C     -----------------------------------------------------------------
C
C     CONSTANTS AND PARAMETERS USED IN THE PHOTOSYNTHESIS MODEL. ALSO
C     NOTE THE STRUCTURE OF VECTORS WHICH CLEARLY SHOWS THE CLASS
C     PFTs (ALONG ROWS) AND CTEM SUB-PFTs (ALONG COLUMNS)
C
C     NEEDLE LEAF |  EVG       DCD       ---
C     BROAD LEAF  |  EVG   DCD-CLD   DCD-DRY
C     CROPS       |   C3        C4       ---
C     GRASSES     |   C3        C4       ---
C
C     CANOPY LIGHT/NITROGEN EXTINCTION COEFFICIENT - THIS BASICALLY
C     ASSUMES THAT MEAN PROFILE OF NITROGEN IS SAME AS THAT FOR
C     TIME MEAN PROFILE OF RADIATION - THE ASSUMPTION MADE BY SINGLE
C     BIG-LEAF MODELS
      DATA   KN/0.50, 0.50, 0.00,
     &          0.50, 0.50, 0.50,
     &          0.40, 0.48, 0.00,
     &          0.46, 0.44, 0.00/
C
C     LOWER AND UPPER TEMPERATURE LIMITS FOR PHOTOSYNTHESIS, KELVIN
C     LOWER LIMIT IN CELCIUS /-5, -5, --,
C                              0,  0,  0,
C                             -3,  5, --,
C                             -1, 10, --/
      DATA TLOW/268.1, 268.1, 0.000,
     &          273.1, 273.1, 273.1,
     &          270.1, 278.1, 0.000,
     &          272.1, 283.1, 0.000/
C
C     UPPER LIMIT IN CELCIUS /34, 34, --,
C                             37, 37, 37,
C                             42, 42, --,
C                             40, 50, --/
      DATA  TUP/307.1, 307.1, 0.000,
     &          310.1, 310.1, 310.1,
     &          315.1, 315.1, 0.000,
     &          313.1, 323.1, 0.000/
C
C     ARRAY TELLING WHICH VEGETATION TYPE IS C4
      DATA  ISC4/0, 0, 0,
     &           0, 0, 0,
     &           0, 1, 0,
     &           0, 1, 0/
C
C     QUANTUM EFFICIENCIES, VALUES OF 0.08 & 0.04 ARE USED FOR C3 AND
C     C4 PLANTS, RESPECTIVELY
      DATA  ALPHA/0.08, 0.08, 0.00,
     &            0.08, 0.08, 0.08,
     &            0.08, 0.04, 0.00,
     &            0.08, 0.04, 0.00/
C
C     LEAF SCATERRING COEFFICIENTS, VALUES OF 0.15 & 0.17 ARE USED
C     FOR C3 AND C4 PLANTS, RESPECTIVELY
      DATA  OMEGA/0.15, 0.15, 0.00,
     &            0.15, 0.15, 0.15,
     &            0.15, 0.17, 0.00,
     &            0.15, 0.17, 0.00/
C
C     PARAMETER M USED IN BWB PHOTOSYNTHESIS-STOMATAL CONDUCTANCE
C     COUPLING. ASSUMING THAT ALL NEEDLE LEAF ARE CONIFERS ESSENTIALLY
C     AND USING A VALUE OF 6.0 FOR NEEDLE LEAF.
C
      DATA  MM/9.0, 9.0, 0.0,
     &        12.0,12.0,12.0,
     &        12.0, 6.0, 0.0,
     &        12.0, 6.0, 0.0/
C
C     PARAMETER B USED IN BWB PHOTOSYNTHESIS-STOMATAL CONDUCTANCE
C     COUPLING.
      DATA  BB/0.01, 0.01, 0.00,
     &         0.01, 0.01, 0.01,
     &         0.01, 0.04, 0.00,
     &         0.01, 0.04, 0.00/
C
C     PARAMETER VPD0 USED IN LEUNING TYPE PHOTOSYNTHESIS - STOMATAL
C     CONDUCTANCE COUPLING, IN PASCALS
      DATA VPD0/2000., 2000., 0.000,
     &          2000., 2000., 2000.,
     &          1500., 1500., 0.000,
     &          1500., 1500., 0.000/
C
C     EXPONENT FOR SOIL MOISTURE STRESS. FOR SN EQUAL TO 1, PHOTOSYNTHESIS
C     DECREASES LINEARLY WITH SOIL MOISTURE, AND OF COURSE NON-LINEARLY
C     FOR VALUES HIGHER THAN 1. WHEN SN IS ABOUT 10, PHOTOSYNTHESIS DOES
C     NOT START DECREASING UNTIL ABOUT SOIL MOISTURE IS HALF WAY BETWEEN
C     WILTING POINT AND FIELD CAPACITY.
C
      DATA SN/2.0, 2.0, 0.0,
     &        2.0, 2.0, 2.0,
     &        2.0, 2.0, 0.0,
     &        2.0, 2.0, 0.0/
C
C     ADDITIONAL CONSTRAIN OF SOIL MOISTURE STRESS ON PHOTOSYNTHESIS.
C     THIS CAN BE USED TO SIMULATE THE EFFECT OF IRRIGATION FOR CROPS.
C
      DATA SMSCALE/0.0, 0.0, 0.0,
     &             0.0, 0.0, 0.0,
     &             0.1, 0.1, 0.0,
     &             0.0, 0.0, 0.0/
C
C     MAX. PHOTOSYNTHETIC RATE, MOL CO2 M^-2 S^-1
      DATA VMAX/42.0E-06,  70.0E-06, 0.00E-06,
c      DATA VMAX/50.0E-06,  70.0E-06, 0.00E-06,
     &          35.0E-06,  65.0E-06, 65.0E-06,
     &          80.0E-06,  40.0E-06, 0.00E-06,
     &          70.0E-06,  30.0E-06, 0.00E-06/
C
C     NO. OF ITERATIONS FOR CALCULATING INTERCELLULAR CO2 CONCENTRATION
      DATA  REQITER/4/
C
C     MAX. INTERCELLULAR CO2 CONCENTRATION, PASCALS
      DATA CO2IMAX/2000.00/
C
C     PHOTOSYNTHESIS COUPLING OR CURVATURE COEFFICIENTS
      DATA BETA1/0.950/
      DATA BETA2/0.990/
C
C     PARAMETER TO INITIALIZE INTERCELLULAR CO2 CONC.
      DATA  INICO2I/0.65, 0.65, 0.00,
     &              0.65, 0.65, 0.65,
     &              0.65, 0.37, 0.00,
     &              0.65, 0.37, 0.00/
C
C     LEAF MAINTENANCE RESPIRATION COEFFICIENTS
      DATA  RMLCOEFF/0.015, 0.015, 0.000,
     &               0.015, 0.015, 0.015,
     &               0.015, 0.025, 0.000,
     &               0.015, 0.025, 0.000/
C
C     FREEZING TEMPERATURE
      DATA TFREZ/273.16/
C
C     STANDARD ATMOS. PRESSURE
      DATA STD_PRESS/101325.0/
C
C     ZERO
      DATA ZERO/1E-20/
C
C     ADDITIONAL PARAMETERS FOR TWO-LEAF MODEL
C     LEAF ANGLE DISTRIBUTION
      DATA  CHI/0.01,  0.01, 0.00,
     &          0.17,  0.17, 0.17,
     &         -0.30, -0.30, 0.00,
     &         -0.30, -0.30, 0.00/
C
C     PHOTOSYNTHESIS DOWN REGULATION PARAMETERS
C     EQUIVALENT CO2 FERTILIZATION EFFECT THAT WE WANT MODEL TO YIELD
      DATA GAMMA_W/0.45/
C     EQUIVALENT CO2 FERTILIZATION EFFECT THAT MODEL ACTUALLY GIVES
C     WITHOUT ANY PHOTOSYNTHESIS DOWN-REGULATION
      DATA GAMMA_M/0.95/

C     --------------------------------------------------------------
C     DECIDE HERE IF WE WANT TO USE SINGLE LEAF OR TWO-LEAF MODEL
C     CHOOSE 1 FOR SINGLE-LEAF MODEL, AND 2 FOR TWO-LEAF MODEL
C
      LEAFOPT=1
C
C     DECIDE IF WE WANT TO USE BWB (1) OR LEUNING TYPE (2) PHOTOSYNTHESIS
C     STOMATAL CONDUCTANCE COUPLING
C
      PS_COUP=2
C
C     DECIDE IF WE WANT TO ESTIMATE PHOTOSYNTHETIC RATE AS A SMOOTHED
C     AVERAGE OF THE 3 LIMITING RATES, OR AS A MINIMUM OF THE 3 LIMITING
C     RATES, OR AS A MINIMUM OF THE RUBSICO AND LIGHT LIMITING RATES
C
      SMOOTH=.TRUE.
      MIN2=.FALSE.
      MIN3=.FALSE.
C
C     --------------------------------------------------------------
c     test by replacing TCAN with TA
c
      DO I=IL1,IL2
       TTEMP(I)=TCAN(I)
       TCAN(I)=TA(I)
      ENDDO
c
C
C     INITIALIZATION
C
      IF(LEAFOPT.EQ.1)THEN
C       INITIALIZE REQUIRED ARRAYS TO ZERO FOR SINGLE LEAF MODEL
        DO 100 I = IL1, IL2
          CO2A(I) = 0.0
          IPAR(I) = 0.0
          SIGMA(I) = 0.0
          TGAMMA(I) = 0.0
          KC(I) = 0.0
          KO(I) = 0.0
          GB(I)=0.0
          RC(I)=0.0
          FC_TEST(I)=0.0
100     CONTINUE
C
        DO 200 J = 1, ICC
          USEBB(J)=0.0
          DO 201 I = IL1, IL2
            FPAR(I,J) = 0.0
            VMAXC(I,J) = 0.0
            VMUNS1(I,J) = 0.0
            VMUNS2(I,J) = 0.0
            VMUNS3(I,J) = 0.0
            VMUNS(I,J) = 0.0
            AVE_SM_FUNC(I,J) = 0.0
            VM(I,J) = 0.0
            JC1(I,J) = 0.0
            JC2(I,J) = 0.0
            JC3(I,J) = 0.0
            JC(I,J) = 0.0
            JE(I,J) = 0.0
            JE1(I,J) = 0.0
            JE2(I,J) = 0.0
            JS(I,J) = 0.0
            A_VEG(I,J) = 0.0
            RML_VEG(I,J)=0.0
            AN_VEG(I,J)=0.0
            CO2LS(I,J)=0.0
            GC(I,J)=0.0
            GCTU(I,J)=0.0
            RC_VEG(I,J)=5000.0
            USESLAI(I,J)=0
            USEAILCG(I,J)=0.0
201       CONTINUE
200     CONTINUE
C
      ELSE IF(LEAFOPT.EQ.2)THEN
C       INITIALIZE ARRAYS FOR THE TWO LEAF MODEL
        DO 210 I = IL1, IL2
          CO2A(I) = 0.0
          IPAR(I) = 0.0
          SIGMA(I) = 0.0
          TGAMMA(I) = 0.0
          KC(I) = 0.0
          KO(I) = 0.0
          GB(I)=0.0
          RC(I)=0.0
          IPAR_SUN(I) = 0.0
          IPAR_SHA(I) = 0.0
          FC_TEST(I)=0.0
210     CONTINUE
C
        DO 230 J = 1, ICC
          USEBB(J)=0.0
          DO 220 I = IL1, IL2
            GDIR(I,J) = 0.0
            KB(I,J) = 0.0
            AILCG_SUN(I,J) = 0.0
            AILCG_SHA(I,J) = 0.0
            FPAR_SUN(I,J) = 0.0
            FPAR_SHA(I,J) = 0.0
            VMAXC_SUN(I,J) = 0.0
            VMAXC_SHA(I,J) = 0.0
            VMUNS1_SUN(I,J) = 0.0
            VMUNS1_SHA(I,J) = 0.0
            VMUNS_SUN(I,J) = 0.0
            VMUNS_SHA(I,J) = 0.0
            AVE_SM_FUNC(I,J) = 0.0
            VM_SUN(I,J) = 0.0
            VM_SHA(I,J) = 0.0
            JC1_SUN(I,J) = 0.0
            JC1_SHA(I,J) = 0.0
            JC3_SUN(I,J) = 0.0
            JC3_SHA(I,J) = 0.0
            JC_SUN(I,J) = 0.0
            JC_SHA(I,J) = 0.0
            JE_SUN(I,J) = 0.0
            JE_SHA(I,J) = 0.0
            JE1_SUN(I,J) = 0.0
            JE1_SHA(I,J) = 0.0
            JE2_SUN(I,J) = 0.0
            JE2_SHA(I,J) = 0.0
            JS_SUN(I,J) = 0.0
            JS_SHA(I,J) = 0.0
            A_VEG_SUN(I,J) = 0.0
            A_VEG_SHA(I,J) = 0.0
            RML_SUN(I,J) = 0.0
            RML_SHA(I,J) = 0.0
            AN_SUN(I,J) = 0.0
            AN_SHA(I,J) = 0.0
            CO2LS_SUN(I,J) = 0.0
            CO2LS_SHA(I,J) = 0.0
            GC_SUN(I,J) = 0.0
            GC_SHA(I,J) = 0.0
            GCTU_SUN(I,J) = 0.0
            GCTU_SHA(I,J) = 0.0
            RC_VEG(I,J)=5000.0
            AN_VEG(I,J)=0.0
            RML_VEG(I,J)=0.0
            USESLAI(I,J)=0
            USEAILCG(I,J)=0.0
220       CONTINUE
230     CONTINUE
      ENDIF
C
C     FOLLOWING VARIABLES AND CONSTANTS ARE COMMON TO BOTH SINGLE AND TWO-LEAF MODEL
      DO 240 J = 1, IG
        DO 250 I = IL1, IL2
          PSISAT(I,J) = 0.0
          GRKSAT(I,J) = 0.0
          B(I,J) = 0.0
          WILTSM(I,J)=0.0
          FIELDSM(I,J)=0.0
          SM_FUNC(I,J) = 0.0
          ISAND(I,J)=NINT(SAND(I,J))
250     CONTINUE
240   CONTINUE
C
      DO 260 J = 1, ICC
        DO 270 I = IL1, IL2
          GCMIN(I,J) = 0.0
          GCMAX(I,J) = 0.0
270     CONTINUE
260   CONTINUE
C
C     GENERATE THE SORT INDEX FOR CORRESPONDENCE BETWEEN 9 PFTs AND THE
C     12 VALUES IN THE PARAMETER VECTORS
C
      ICOUNT=0
      DO 280 J = 1, IC
        DO 281 M = 1, NOL2PFTS(J)
          N = (J-1)*L2MAX + M
          ICOUNT = ICOUNT + 1
          SORT(ICOUNT)=N
281     CONTINUE
280   CONTINUE
C
C     INITIALIZATION ENDS
C
C     -------------------------------------------------------------------
C
C     IF LAI IS LESS THAN SLAI THAN WE USE STORAGE LAI TO PHOTOSYNTHESIZE.
C     HOWEVER, WE DO NOT USE THE STOMATAL RESISTANCE ESTIMATED IN THIS CASE,
C     BECAUSE STORAGE LAI IS AN IMAGINARY LAI, AND WE SET STOMATAL RESISTANCE
C     TO ITS MAX. VALUE. IF, HOWEVER, CTEM1 IS TRUE AND CTEM2 IS FALSE, WHICH
C     MEANS WE ARE ONLY USING STOMATAL RESISTANCE FROM CTEM, THEN THE CONCEPT
C     OF SLAI DOES NOT APPLY AND WE MAKE THIS ZERO. NOTE THAT THE CONCEPT OF
C     STORAGE/IMAGINARY LAI IS USED FOR PHENOLOGY PURPOSES AND THIS
C     IMAGINARY LAI ACTS AS MODEL SEEDS.
C
      DO 320 J = 1, ICC
        DO 330 I = IL1, IL2
          IF( CTEM1 .AND. (.NOT.CTEM2) ) THEN
            SLAI(I,J)=0.0
          ENDIF
330     CONTINUE
320   CONTINUE
C
      DO 340 J = 1, ICC
        DO 350 I = IL1, IL2
          IF(AILCG(I,J).LT.SLAI(I,J))THEN
            USESLAI(I,J)=1
            USEAILCG(I,J)=SLAI(I,J)
          ELSE
            USEAILCG(I,J)=AILCG(I,J)
          ENDIF
350     CONTINUE
340   CONTINUE
C
C     SET MIN. AND MAX. VALUES FOR STOMATAL CONDUCTANCE. WE MAKE SURE
C     THAT MAX. STOMATAL RESISTANCE IS AROUND 5000 S/M AND MIN. STOMATAL
C     RESISTANCE IS 51 S/M.
C
      DO 360 J = 1, ICC
        DO 370 I = IL1, IL2
          GCMIN(I,J)=0.0002 * (TFREZ/TCAN(I)) * (1./0.0224) *
     &      (PRESSG(I)/STD_PRESS)
C
          IF(LEAFOPT.EQ.1)THEN
C           GCMAX(I,J)=0.0196 * (TFREZ/TCAN(I)) * (1./0.0224) *
            GCMAX(I,J)=0.1    * (TFREZ/TCAN(I)) * (1./0.0224) *
     &        (PRESSG(I)/STD_PRESS)
          ELSEIF(LEAFOPT.EQ.2)THEN
C           GCMAX(I,J)=0.0196 * (TFREZ/TCAN(I)) * (1./0.0224) *
            GCMAX(I,J)=0.1    * (TFREZ/TCAN(I)) * (1./0.0224) *
     &        (PRESSG(I)/STD_PRESS) * 0.5
          ENDIF
C
370     CONTINUE
360   CONTINUE
C
C     IF WE ARE USING LEUNING TYPE PHOTOSYNTHESIS-STOMATAL CONDUCTANCE
C     COUPLING WE NEED VAPOR PRESSURE DEFICIT AS WELL. CALCULATE THIS
C     FROM THE RH AND AIR TEMPEATURE WE HAVE. WE FIND E_SAT, E, AND VPD
C     IN PASCALS.
C
      IF(PS_COUP.EQ.2)THEN
        DO 390 J = 1, ICC
          DO 400 I = IL1, IL2
            VPD_TERM(I,J)=0.0
400       CONTINUE
390     CONTINUE
C
        DO 420 I = IL1, IL2
          VPD(I)=0.0
          IF(TA(I).GE.TFREZ) THEN
              CA=17.269
              CB=35.86
          ELSE
              CA=21.874
              CB=7.66
          ENDIF
          EA(I)  = QA(I)*PRESSG(I)/(0.622+0.378*QA(I))
          EASAT(I)=611.0*EXP(CA*(TA(I)-TFREZ)/(TA(I)-CB))
          RH(I)=EA(I)/EASAT(I)
          VPD(I)=EASAT(I)-EA(I)
          VPD(I)=MAX(0.0,VPD(I))
420     CONTINUE
C
        K1C = 0
        DO 440 J = 1, IC
          IF(J.EQ.1) THEN
            K1C = K1C + 1
          ELSE
            K1C = K1C + NOL2PFTS(J-1)
          ENDIF
          K2C = K1C + NOL2PFTS(J) - 1
          DO 445 M = K1C, K2C
            DO 450 I = IL1, IL2
              VPD_TERM(I,M)=1.0/( 1.0 +( VPD(I)/VPD0(SORT(M)) ) )
450         CONTINUE
445       CONTINUE
440     CONTINUE
      ENDIF
C
C     ESTIMATE PARTIAL PRESSURE OF CO2 AND IPAR
C
      DO 460 I = IL1, IL2
C       CONVERT CO2CONC FROM PPM TO PASCALS
        CO2A(I)=CO2CONC(I)*PRESSG(I)/1E+06
C       CHANGE PAR FROM W/M^2 TO MOL/M^2.S
        IPAR(I) = QSWV(I)*4.6*1E-06
C
C       SUNLIT PART GETS BOTH DIRECT AND DIFFUSED, WHILE
C       THE SHADED PART GETS ONLY DIFFUSED
C
        IPAR_SUN(I) = QSWV(I)*4.6*1E-06
        IPAR_SHA(I) = QSWV(I)*4.6*1E-06* XDIFFUS(I)
460   CONTINUE
C
      K1C = 0
      DO 480 J = 1,IC
       IF(J.EQ.1) THEN
         K1C = K1C + 1
       ELSE
         K1C = K1C + NOL2PFTS(J-1)
       ENDIF
       K2C = K1C + NOL2PFTS(J) - 1
       DO 485 M = K1C, K2C
        DO 490 I = IL1, IL2
        IF(FCANC(I,M).GT.ZERO)THEN
C
C         FOR TWO-LEAF MODEL FIND Kb AS A FUNCTION OF COSZS AND
C         LEAF ANGLE DISTRIBUTION (VEGETATION DEPENDENT)
C
          IF(LEAFOPT.EQ.2)THEN
           IF(COSZS(I).GT.0.0)THEN
C
C           MAKE SURE -0.4 < CHI < 0.6
            CHI(SORT(M))=MIN (MAX (CHI(SORT(M)), -0.4), 0.6)
C           MAKE VALUES CLOSE TO ZERO EQUAL TO 0.01
            IF( ABS(CHI(SORT(M))).LE.0.01 ) CHI(SORT(M))=0.01
C
            TEMP_PHI1 =
     &       0.5-0.633*CHI(SORT(M))-0.33*CHI(SORT(M))*CHI(SORT(M))
            TEMP_PHI2 = 0.877*(1.-2.*TEMP_PHI1)
            GDIR(I,M) = TEMP_PHI1 + TEMP_PHI2*COSZS(I)
            KB(I,M) = (GDIR(I,M)/COSZS(I))
            KB(I,M) = KB(I,M) * (SQRT(1.-OMEGA(SORT(M)) ))
C
C           ALSO FIND SUNLIT AND SHADED LAI
            AILCG_SUN(I,M)=(1/KB(I,M)) *
     &        ( 1.0-EXP( -1.0*KB(I,M)*USEAILCG(I,M) ) )
            AILCG_SHA(I,M)=USEAILCG(I,M) - AILCG_SUN(I,M)
C
C           FOLLOWING FEW LINES TO MAKE SURE THAT ALL LEAVES ARE
C           SHADED WHEN XDIFFUS EQUALS 1. NOT DOING SO GIVES ERRATIC
C           RESULTS WHEN TWO LEAF OPTION IS USED
C
            IF(XDIFFUS(I).GT.0.99)THEN
              AILCG_SUN(I,M)=0.0
              AILCG_SHA(I,M)=USEAILCG(I,M)
            ENDIF
C
           ENDIF
          ENDIF
C
C         FIND FPAR - FACTOR FOR SCALING PHOTOSYNTHESIS TO CANOPY
C         BASED ON ASSUMPTION THAT NITROGEN IS OPTIMALLY DISTRIBUTED.
C         THE TWO-LEAF MODEL IS NOT THAT DIFFERENT FROM THE SINGLE-LEAF
C         MODEL. ALL WE DO IS USE TWO SCALING FACTORS (I.E. SCALING FROM
C         LEAF TO CANOPY) INSTEAD OF ONE, AND THUS PERFORM CALCULATIONS
C         TWICE, AND IN THE END ADD CONDUCTANCE AND NET PHOTOSYNTHESIS
C         FROM THE TWO LEAVES TO GET THE TOTAL.
C
          FPAR(I,M)=(1/KN(SORT(M)))*(1-EXP(-KN(SORT(M))*USEAILCG(I,M)))
          IF(LEAFOPT.EQ.2)THEN
            FPAR_SUN(I,M) = ( 1/(KN(SORT(M))+KB(I,M)) )*
     &        ( 1-EXP( -1.*(KN(SORT(M))+KB(I,M))*USEAILCG(I,M) ) )
            FPAR_SHA(I,M) = FPAR(I,M) - FPAR_SUN(I,M)
C
C           IF ALL RADIATION IS DIFFUSED, THEN ALL LEAVES ARE SHADED,
C           AND WE ADJUST FPARs ACCORDINGLY. WITHOUT THIS THE TWO LEAF
C           MODELS MAY BEHAVE ERRATICALLY
C
            IF(XDIFFUS(I).GT.0.99)THEN
              FPAR_SHA(I,M) = FPAR(I,M)
              FPAR_SUN(I,M) = 0.0
            ENDIF
          ENDIF
C
C         FIND Vmax,canopy, THAT IS Vmax SCALED BY LAI FOR THE SINGLE
C         LEAF MODEL
C
          VMAXC(I,M)=VMAX(SORT(M)) * FPAR(I,M)
C
C         FIND Vm,unstressed (DUE TO WATER) BUT STRESSED DUE TO TEMPERATURE
C
          Q10 = 2.00
          Q10_FUNC = Q10**(0.1*(TCAN(I)-298.16))
C
          IF(COSZS(I).GT.0.0)THEN
            IF(LEAFOPT.EQ.1)THEN
              VMUNS1(I,M) = VMAXC(I,M) * Q10_FUNC
            ELSE IF(LEAFOPT.EQ.2)THEN
              VMUNS1_SUN(I,M) = VMAXC_SUN(I,M) * Q10_FUNC
              VMUNS1_SHA(I,M) = VMAXC_SHA(I,M) * Q10_FUNC
            ENDIF
          ENDIF
C
C         ASSUMING THAT SUNLIT AND SHADED TEMEPERATURES ARE SAME
C
          VMUNS2(I,M) = (1 + EXP(0.3*(TCAN(I)-TUP(SORT(M)))))
          VMUNS3(I,M) = (1 + EXP(0.3*(TLOW(SORT(M))-TCAN(I))))
C
          IF(LEAFOPT.EQ.1)THEN
            VMUNS(I,M)=VMUNS1(I,M)/
     &        (VMUNS2(I,M)*VMUNS3(I,M))
          ELSE IF(LEAFOPT.EQ.2)THEN
            VMUNS_SUN(I,M)=VMUNS1_SUN(I,M)/
     &        (VMUNS2(I,M)*VMUNS3(I,M))
            VMUNS_SHA(I,M)=VMUNS1_SHA(I,M)/
     &        (VMUNS2(I,M)*VMUNS3(I,M))
          ENDIF
C
        ENDIF
490     CONTINUE
485    CONTINUE
480   CONTINUE
C
C     CALCULATE SOIL MOIS STRESS TO ACCOUNT FOR REDUCTION IN PHOTOSYN
C     DUE TO LOW SOIL MOISTURE, THREE STEPS HERE -> 1. FIND WILTING
C     POINT AND FIELD CAPACITY SOIL MOIS. CONTENT FOR ALL THREE LAYERS.
C     2. USING THESE FIND THE SOIL MOISTURE STRESS TERM FOR ALL
C     THREE LAYERS, AND 3. AVERAGE THIS SOIL MOISTURE STRESS TERM
C     OVER THE 3 LAYERS USING FRACTION OF ROOTS PRESENT IN EACH LAYER
C     FOR EACH PFT. NOTE THAT WHILE SOIL MOISTURE IS UNIFORM OVER
C     AN ENTIRE GCM GRID CELL, THE SOIL MOISTURE STRESS FOR EACH
C     PFT IS NOT BECAUSE OF DIFFERENCES IN ROOT DISTRIBUTION.
C
C     WILTING POINT CORRESPONDS TO MATRIC POTENTIAL OF 150 M
C     FIELD CAPACITY CORRESPONDS TO HYDARULIC CONDUCTIVITY OF
C     0.10 MM/DAY -> 1.157x1E-09 M/S
C
      DO 500 J = 1, IG
        DO 510 I = IL1, IL2
C
          IF(NINT(SAND(I,J)).EQ.-3.OR.NINT(SAND(I,J)).EQ.-4)THEN
            SM_FUNC(I,J)=0.01
          ELSE ! I.E., SAND.NE.-3 OR -4
C
           PSISAT(I,J)= (10.0**(-0.0131*SAND(I,J)+1.88))/100.0
           GRKSAT(I,J)= (10.0**(0.0153*SAND(I,J)-0.884))*7.0556E-6
           THPOR(I,J) = (-0.126*SAND(I,J)+48.9)/100.0
           B(I,J)     = 0.159*CLAY(I,J)+2.91
C
           WILTSM(I,J) = (150./PSISAT(I,J))**(-1/B(I,J))
           WILTSM(I,J) = THPOR(I,J) * WILTSM(I,J)
C
           FIELDSM(I,J) = (1.157E-09/GRKSAT(I,J))**
     &      (1/(2*B(I,J)+3))
           FIELDSM(I,J) = THPOR(I,J) *  FIELDSM(I,J)
C
           IF(THLIQ(I,J).LE.WILTSM(I,J)) THEN
            SM_FUNC(I,J)=0.01
           ELSE IF(THLIQ(I,J).GT.WILTSM(I,J) .AND.
     &      THLIQ(I,J).LT.FIELDSM(I,J)) THEN
            SM_FUNC(I,J)=(THLIQ(I,J)-WILTSM(I,J))/
     &          (FIELDSM(I,J)-WILTSM(I,J))
           ELSE IF(THLIQ(I,J).GE.FIELDSM(I,J)) THEN
            SM_FUNC(I,J)=1.0
           ENDIF
C
          ENDIF ! SAND.EQ.-3 OR -4
C
510     CONTINUE
500   CONTINUE
C
      K1C=0
      DO 520 J = 1, IC
       IF(J.EQ.1) THEN
         K1C = K1C + 1
       ELSE
         K1C = K1C + NOL2PFTS(J-1)
       ENDIF
       K2C = K1C + NOL2PFTS(J) - 1
       DO 525 M = K1C, K2C
        DO 530 I = IL1, IL2
         SM_FUNC(I,1)=( 1.0 - (1.0-SM_FUNC(I,1))**SN(SORT(M)) )
         SM_FUNC(I,2)=( 1.0 - (1.0-SM_FUNC(I,2))**SN(SORT(M)) )
         SM_FUNC(I,3)=( 1.0 - (1.0-SM_FUNC(I,3))**SN(SORT(M)) )

         SM_FUNC(I,1)=SM_FUNC(I,1)+(1.0-SM_FUNC(I,1))*SMSCALE(SORT(M))
         SM_FUNC(I,2)=SM_FUNC(I,2)+(1.0-SM_FUNC(I,2))*SMSCALE(SORT(M))
         SM_FUNC(I,3)=SM_FUNC(I,3)+(1.0-SM_FUNC(I,3))*SMSCALE(SORT(M))

         AVE_SM_FUNC(I,M)=SM_FUNC(I,1)*RMAT(I,M,1) +
     &                    SM_FUNC(I,2)*RMAT(I,M,2) +
     &                    SM_FUNC(I,3)*RMAT(I,M,3)
         AVE_SM_FUNC(I,M)= AVE_SM_FUNC(I,M) /
     &    (RMAT(I,M,1)+RMAT(I,M,2)+RMAT(I,M,3))
C
C        MAKE SURE THAT I DO NOT USE SM_FUNC FROM 2nd AND 3rd LAYERS
C        IF THEY ARE BED ROCK
C
         IF(ISAND(I,3).EQ.-3)THEN ! THIRD LAYER BED ROCK
           AVE_SM_FUNC(I,M)=SM_FUNC(I,1)*RMAT(I,M,1) +
     &                      SM_FUNC(I,2)*RMAT(I,M,2)
C
           AVE_SM_FUNC(I,M)= AVE_SM_FUNC(I,M) /
     &      (RMAT(I,M,1)+RMAT(I,M,2))
         ENDIF
         IF(ISAND(I,2).EQ.-3)THEN ! SECOND LAYER BED ROCK
           AVE_SM_FUNC(I,M)=SM_FUNC(I,1)
         ENDIF
C
         IF( (RMAT(I,M,1)+RMAT(I,M,2)+RMAT(I,M,3)).LT.0.9) THEN
           WRITE(6,*)'PFT = ',M,' I =',I
           WRITE(6,*)'RMAT ADD =',(RMAT(I,M,1)+RMAT(I,M,2)+RMAT(I,M,3))
           CALL XIT('PHTSYN', -99)
         ENDIF
C
530     CONTINUE
525    CONTINUE
520   CONTINUE
C
C     USE SOIL MOISTURE FUNCTION TO MAKE Vm,unstressed -> Vm STRESSED
C
      DO 540 J = 1, ICC
        DO 550 I = IL1, IL2
        IF(FCANC(I,J).GT.ZERO)THEN
C
          IF(COSZS(I).GT.0.0)THEN
            IF(LEAFOPT.EQ.1)THEN
              VM(I,J) = VMUNS(I,J) * AVE_SM_FUNC(I,J)
            ELSE IF(LEAFOPT.EQ.2)THEN
              VM_SUN(I,J) = VMUNS_SUN(I,J) * AVE_SM_FUNC(I,J)
              VM_SHA(I,J) = VMUNS_SHA(I,J) * AVE_SM_FUNC(I,J)
            ENDIF
          ENDIF
C
        ENDIF
550     CONTINUE
540   CONTINUE
C
C     FIND TEMPERATURE DEPENDENT PARAMETER VALUES
C
      DO 570 I = IL1, IL2
C
C       FIND RUBISCO SPECIFICITY FOR CO2 RELATIVE TO O2 - SIGMA
C
        Q10 = 0.57
        Q10_FUNC = Q10**(0.1*(TCAN(I)-298.16))
        SIGMA(I) = 2600.0 * Q10_FUNC
C
C       FIND CO2 COMPENSATION POINT USING RUBISCO SPECIFICITY - TGAMMA.
C       KEEP IN MIND THAT CO2 COMPENSATION POINT FOR C4 PLANTS IS ZERO,
C       SO THE FOLLOWING VALUE IS RELEVANT FOR C3 PLANTS ONLY
C
        O2_CONC(I) = .2095 * PRESSG(I)
        TGAMMA(I) = O2_CONC(I) / (2*SIGMA(I))
C
C       ESTIMATE MICHELIS-MENTON CONSTANTS FOR CO2 (Kc) and O2 (Ko) TO
C       BE USED LATER FOR ESTIMATING RUBISCO LIMITED PHOTOSYNTHETIC RATE
C
        Q10 = 2.10
        Q10_FUNC = Q10**(0.1*(TCAN(I)-298.16))
        KC(I) = 30.0 * Q10_FUNC
C
        Q10 = 1.20
        Q10_FUNC = Q10**(0.1*(TCAN(I)-298.16))
        KO(I) = 30000.0 * Q10_FUNC
C
570   CONTINUE
C
C     CHOOSE A VALUE OF INTERCELLULAR CO2 CONCENTRATION (CO2i) IF STARTING
C     FOR THE FIRST TIME, OR USE VALUE FROM THE PREVIOUS TIME STEP
C
      IT_COUNT = 0
999   CONTINUE
C
      DO 580 J = 1,ICC
        DO 590 I = IL1, IL2
C
          IF(LEAFOPT.EQ.1)THEN
            CO2I(I,J) = CO2I1(I,J)
            IF(CO2I(I,J).LE.ZERO) THEN
              CO2I(I,J) = INICO2I(SORT(J)) * CO2A(I)
            ENDIF
          ELSE IF(LEAFOPT.EQ.2)THEN
            CO2I_SUN(I,J) = CO2I1(I,J)
            IF(CO2I_SUN(I,J).LE.ZERO) THEN
              CO2I_SUN(I,J) = INICO2I(SORT(J)) * CO2A(I)
            ENDIF
            CO2I_SHA(I,J) = CO2I2(I,J)
            IF(CO2I_SHA(I,J).LE.ZERO) THEN
              CO2I_SHA(I,J) = INICO2I(SORT(J)) * CO2A(I)
            ENDIF
          ENDIF
C
590     CONTINUE
580   CONTINUE
C
C     ESTIMATE RUBISCO LIMITED PHOTOSYNTHETIC RATE
C
      DO 600 J = 1,ICC
        DO 610 I = IL1, IL2
        IF(FCANC(I,J).GT.ZERO)THEN
C
          IF(COSZS(I).GT.0.0)THEN
           IF(LEAFOPT.EQ.1)THEN
            JC1(I,J) = CO2I(I,J) - TGAMMA(I)
            JC2(I,J) = KC(I) * (1 + (O2_CONC(I)/KO(I)) )
            JC3(I,J) = CO2I(I,J) + JC2(I,J)
C
            IF (ISC4(SORT(J)).EQ.1) THEN
              JC(I,J)  = VM(I,J)
            ELSE IF (ISC4(SORT(J)).EQ.0) THEN
              JC(I,J)  = VM(I,J) * (JC1(I,J)/JC3(I,J))
            ENDIF
           ELSE IF(LEAFOPT.EQ.2)THEN
            JC1_SUN(I,J) = CO2I_SUN(I,J) - TGAMMA(I)
            JC1_SHA(I,J) = CO2I_SHA(I,J) - TGAMMA(I)
            JC2(I,J) = KC(I) * (1 + (O2_CONC(I)/KO(I)) )
            JC3_SUN(I,J) = CO2I_SUN(I,J) + JC2(I,J)
            JC3_SHA(I,J) = CO2I_SHA(I,J) + JC2(I,J)
C
            IF (ISC4(SORT(J)).EQ.1) THEN
              JC_SUN(I,J)=VM_SUN(I,J)
              JC_SHA(I,J)=VM_SHA(I,J)
            ELSE IF (ISC4(SORT(J)).EQ.0) THEN
              JC_SUN(I,J)=VM_SUN(I,J)*(JC1_SUN(I,J)/JC3_SUN(I,J))
              JC_SHA(I,J)=VM_SHA(I,J)*(JC1_SHA(I,J)/JC3_SHA(I,J))
            ENDIF
           ENDIF
          ENDIF
C
        ENDIF
610     CONTINUE
600   CONTINUE
C
C     ESTIMATE PHOTOSYNTHETIC RATE LIMITED BY AVAILABLE LIGHT
C
      DO 620 J = 1,ICC
        DO 630 I = IL1, IL2
        IF(FCANC(I,J).GT.ZERO)THEN
C
          IF(COSZS(I).GT.0.0)THEN
           IF(LEAFOPT.EQ.1)THEN
            JE1(I,J)= FPAR(I,J)*ALPHA(SORT(J))*(1-OMEGA(SORT(J)))
            JE2(I,J)=( CO2I(I,J)-TGAMMA(I) ) /
     &        ( CO2I(I,J)+ (2*TGAMMA(I)) )
C
            IF (ISC4(SORT(J)).EQ.1) THEN
              JE(I,J) = IPAR(I) * JE1(I,J)
            ELSE IF (ISC4(SORT(J)).EQ.0) THEN
              JE(I,J) = IPAR(I) * JE1(I,J) * JE2(I,J)
            ENDIF
           ELSE IF(LEAFOPT.EQ.2)THEN
            JE1_SUN(I,J)= FPAR_SUN(I,J)*ALPHA(SORT(J))
            JE1_SHA(I,J)= FPAR_SHA(I,J)*ALPHA(SORT(J))
            JE2_SUN(I,J)=( CO2I_SUN(I,J)-TGAMMA(I) ) /
     &        ( CO2I_SUN(I,J)+ (2*TGAMMA(I)) )
            JE2_SHA(I,J)=( CO2I_SHA(I,J)-TGAMMA(I) ) /
     &        ( CO2I_SHA(I,J)+ (2*TGAMMA(I)) )

            IF (ISC4(SORT(J)).EQ.1) THEN
              JE_SUN(I,J)=(IPAR_SUN(I)+IPAR_SHA(I))*JE1_SUN(I,J)
              JE_SHA(I,J)=IPAR_SHA(I)*JE1_SHA(I,J)
            ELSE IF (ISC4(SORT(J)).EQ.0) THEN
              JE_SUN(I,J)=(IPAR_SUN(I)+IPAR_SHA(I))*JE1_SUN(I,J)*
     &                    JE2_SUN(I,J)
              JE_SHA(I,J)=IPAR_SHA(I)*JE1_SHA(I,J)*JE2_SHA(I,J)
            ENDIF
           ENDIF
          ENDIF
C
        ENDIF
630     CONTINUE
620   CONTINUE
C
C     ESTIMATE PHOTOSYNTHETIC RATE LIMITED BY TRANSPORT CAPACITY
C
      DO 640 J = 1,ICC
        DO 650 I = IL1, IL2
        IF(FCANC(I,J).GT.ZERO)THEN
C
          IF(COSZS(I).GT.0.0)THEN
           IF(LEAFOPT.EQ.1)THEN
            IF (ISC4(SORT(J)).EQ.1) THEN
              JS(I,J) = 20000.0 * VM(I,J)*(CO2I(I,J)/PRESSG(I))
            ELSE IF (ISC4(SORT(J)).EQ.0) THEN
              JS(I,J) = 0.5 * VM(I,J)
            ENDIF
           ELSE IF(LEAFOPT.EQ.2)THEN
            IF (ISC4(SORT(J)).EQ.1) THEN
              JS_SUN(I,J)=20000.0*VM_SUN(I,J)*(CO2I_SUN(I,J)/PRESSG(I))
              JS_SHA(I,J)=20000.0*VM_SHA(I,J)*(CO2I_SHA(I,J)/PRESSG(I))
            ELSE IF (ISC4(SORT(J)).EQ.0) THEN
              JS_SUN(I,J) = 0.5 * VM_SUN(I,J)
              JS_SHA(I,J) = 0.5 * VM_SHA(I,J)
            ENDIF
           ENDIF
          ENDIF
C
        ENDIF
650     CONTINUE
640   CONTINUE
C
C     INCLUDE NUTRIENT LIMITATION EFFECT BY DOWN-REGULATING PHOTOSYNTHESIS
C     N_EFFECT DECREASES FROM 1.0 AS CO2 INCREASES ABOVE 288 PPM.
C
      DO 655 I = IL1, IL2

        DELTA_CO2(I) = CO2CONC(I) - 288.0
        TEMP_R = LOG(1.0 + (DELTA_CO2(I)/288.0))
        TEMP_B = 1.0 + TEMP_R * (GAMMA_W)
        TEMP_C = 1.0 + TEMP_R * (GAMMA_M)
        N_EFFECT(I) = TEMP_B / TEMP_C

655   CONTINUE
C
C     FIND THE SMOOTHED AVERAGE OF THREE PHOTOSYNTHETIC RATES JC, JE,
C     AND JS USING COLLATZ'S TWO QUADRATIC EQUATIONS, OR FIND THE MIN.
C     OF THIS TWO RATES OR FIND MIN. OF JC AND JE.
C
      DO 660 J = 1,ICC
        DO 670 I = IL1, IL2
        IF(FCANC(I,J).GT.ZERO)THEN
C
          IF(COSZS(I).GT.0.0)THEN
           IF(LEAFOPT.EQ.1)THEN
            IF(SMOOTH)THEN
              TEMP_B  = 0.0
              TEMP_C  = 0.0
              TEMP_R  = 0.0
              TEMP_Q1 = 0.0
              TEMP_Q2 = 0.0
              TEMP_JP = 0.0
C
              TEMP_B  = JC(I,J) + JE(I,J)
              TEMP_C  = JC(I,J) * JE(I,J)
              TEMP_R  = MAX( (TEMP_B**2 - 4. * BETA1 * TEMP_C), 0.0)
              TEMP_Q1 = (TEMP_B + SQRT(TEMP_R)) / (2 * BETA1)
              TEMP_Q2 = (TEMP_B - SQRT(TEMP_R)) / (2 * BETA1)
              TEMP_JP = MIN(TEMP_Q1, TEMP_Q2)
C
              TEMP_B  = TEMP_JP + JS(I,J)
              TEMP_C  = TEMP_JP * JS(I,J)
              TEMP_R  = MAX( (TEMP_B**2 - 4. * BETA2 * TEMP_C), 0.0)
              TEMP_Q1 = (TEMP_B + SQRT(TEMP_R)) / (2 * BETA2)
              TEMP_Q2 = (TEMP_B - SQRT(TEMP_R)) / (2 * BETA2)
              A_VEG(I,J) = MIN(TEMP_Q1, TEMP_Q2)
            ELSEIF(MIN2) THEN
              A_VEG(I,J)=MIN(JC(I,J), JE(I,J))
            ELSEIF(MIN3) THEN
              A_VEG(I,J)=MIN(JC(I,J), JE(I,J), JS(I,J))
            ELSE
              CALL XIT('PHTSYN',-1)
            ENDIF
C           DOWN-REGULATE PHOTOSYNTHESIS FOR C3 PLANTS
            IF(ISC4(SORT(J)).EQ.0)THEN
              A_VEG(I,J) = A_VEG(I,J) * N_EFFECT(I)
            ENDIF
            A_VEG(I,J) = MAX(0.0, A_VEG(I,J))
           ELSE IF(LEAFOPT.EQ.2)THEN
             IF(SMOOTH)THEN
               TEMP_B  = 0.0
               TEMP_C  = 0.0
               TEMP_R  = 0.0
               TEMP_Q1 = 0.0
               TEMP_Q2 = 0.0
               TEMP_JP = 0.0
C
               TEMP_B  = JC_SUN(I,J) + JE_SUN(I,J)
               TEMP_C  = JC_SUN(I,J) * JE_SUN(I,J)
               TEMP_R  = MAX( (TEMP_B**2 - 4. * BETA1 * TEMP_C), 0.0)
               TEMP_Q1 = (TEMP_B + SQRT(TEMP_R)) / (2 * BETA1)
               TEMP_Q2 = (TEMP_B - SQRT(TEMP_R)) / (2 * BETA1)
               TEMP_JP = MIN(TEMP_Q1, TEMP_Q2)

               TEMP_B  = TEMP_JP + JS_SUN(I,J)
               TEMP_C  = TEMP_JP * JS_SUN(I,J)
               TEMP_R  = MAX( (TEMP_B**2 - 4. * BETA2 * TEMP_C), 0.0)
               TEMP_Q1 = (TEMP_B + SQRT(TEMP_R)) / (2 * BETA2)
               TEMP_Q2 = (TEMP_B - SQRT(TEMP_R)) / (2 * BETA2)
               A_VEG_SUN(I,J) = MIN(TEMP_Q1, TEMP_Q2)
             ELSEIF(MIN2) THEN
               A_VEG_SUN(I,J)=MIN(JC_SUN(I,J), JE_SUN(I,J))
             ELSEIF(MIN3) THEN
               A_VEG_SUN(I,J)=MIN(JC_SUN(I,J),JE_SUN(I,J),JS_SUN(I,J))
             ELSE
               CALL XIT('PHTSYN',-2)
             ENDIF
             IF(ISC4(SORT(J)).EQ.0)THEN
               A_VEG_SUN(I,J) = A_VEG_SUN(I,J) * N_EFFECT(I)
             ENDIF
             A_VEG_SUN(I,J) = MAX(0.0, A_VEG_SUN(I,J))
C
             IF(SMOOTH)THEN
               TEMP_B  = 0.0
               TEMP_C  = 0.0
               TEMP_R  = 0.0
               TEMP_Q1 = 0.0
               TEMP_Q2 = 0.0
               TEMP_JP = 0.0
C
               TEMP_B  = JC_SHA(I,J) + JE_SHA(I,J)
               TEMP_C  = JC_SHA(I,J) * JE_SHA(I,J)
               TEMP_R  = MAX( (TEMP_B**2 - 4. * BETA1 * TEMP_C), 0.0)
               TEMP_Q1 = (TEMP_B + SQRT(TEMP_R)) / (2 * BETA1)
               TEMP_Q2 = (TEMP_B - SQRT(TEMP_R)) / (2 * BETA1)
               TEMP_JP = MIN(TEMP_Q1, TEMP_Q2)
C
               TEMP_B  = TEMP_JP + JS_SHA(I,J)
               TEMP_C  = TEMP_JP * JS_SHA(I,J)
               TEMP_R  = MAX( (TEMP_B**2 - 4. * BETA2 * TEMP_C), 0.0)
               TEMP_Q1 = (TEMP_B + SQRT(TEMP_R)) / (2 * BETA2)
               TEMP_Q2 = (TEMP_B - SQRT(TEMP_R)) / (2 * BETA2)
               A_VEG_SHA(I,J) = MIN(TEMP_Q1, TEMP_Q2)
             ELSEIF(MIN2) THEN
               A_VEG_SHA(I,J)=MIN(JC_SHA(I,J), JE_SHA(I,J))
             ELSEIF(MIN3) THEN
               A_VEG_SHA(I,J)=MIN(JC_SHA(I,J),JE_SHA(I,J),JS_SHA(I,J))
             ENDIF
             IF(ISC4(SORT(J)).EQ.0)THEN
               A_VEG_SHA(I,J) = A_VEG_SHA(I,J) * N_EFFECT(I)
             ENDIF
             A_VEG_SHA(I,J) = MAX(0.0, A_VEG_SHA(I,J))
           ENDIF
          ENDIF
C
        ENDIF
670     CONTINUE
660   CONTINUE
C
C     ESTIMATE LEAF MAINTENANCE RESPIRATION RATES AND NET PHOTOSYNTHETIC
C     RATE. THIS NET PHOSYNTHETIC RATE IS /M^2 OF VEGETATED LAND.
C
      DO 680 J = 1,ICC
        DO 690 I = IL1, IL2
        IF(FCANC(I,J).GT.ZERO)THEN
C
C         RECENT STUDIES SHOW RmL IS LESS TEMPERATURE SENSITIVE THAN
C         PHOTOSYNTHESIS DURING DAY, THAT'S WHY A SMALL Q10 VALUE IS
C         USED DURING DAY.
C
          Q10_FUNCN = 2.00**(0.1*(TCAN(I)-298.16))
          Q10_FUNCD = 1.30**(0.1*(TCAN(I)-298.16))
C
          IF(LEAFOPT.EQ.1)THEN
           IF(COSZS(I).GT.0.0)THEN
            RML_VEG(I,J) = RMLCOEFF(SORT(J))*VMAXC(I,J)*Q10_FUNCD
           ELSE
            RML_VEG(I,J) = RMLCOEFF(SORT(J))*VMAXC(I,J)*Q10_FUNCN* 0.5
           ENDIF
           AN_VEG(I,J) = A_VEG(I,J) - RML_VEG(I,J)
          ELSE IF(LEAFOPT.EQ.2)THEN
           IF(COSZS(I).GT.0.0)THEN
            RML_SUN(I,J) = RMLCOEFF(SORT(J))*VMAXC_SUN(I,J)*Q10_FUNCD
            RML_SHA(I,J) = RMLCOEFF(SORT(J))*VMAXC_SHA(I,J)*Q10_FUNCD
           ELSE
            RML_SUN(I,J)=RMLCOEFF(SORT(J))*VMAXC_SUN(I,J)*Q10_FUNCN*0.5
            RML_SHA(I,J)=RMLCOEFF(SORT(J))*VMAXC_SHA(I,J)*Q10_FUNCN*0.5
           ENDIF
           AN_SUN(I,J) = A_VEG_SUN(I,J) - RML_SUN(I,J)
           AN_SHA(I,J) = A_VEG_SHA(I,J) - RML_SHA(I,J)
          ENDIF
C
        ENDIF
690     CONTINUE
680   CONTINUE
C
C     FIND CO2 CONCENTRATION AT LEAF SURFACE FOR ALL VEGETATION TYPES.
C     ALTHOUGH WE ARE FINDING CO2 CONC AT THE LEAF SURFACE SEPARATELY
C     FOR ALL VEGETATION TYPES, THE BIG ASSUMPTION HERE IS THAT THE
C     AERODYNAMIC CONDUCTANCE IS SAME OVER ALL VEGETATION TYPES. CLASS
C     FINDS AERODYNAMIC RESISTANCE OVER ALL THE 4 SUB-AREAS, BUT NOT
C     FOR DIFFERENT VEGETATION TYPES WITHIN A SUB-AREA.
C     ALSO CHANGE AERODYNAMIC CONDUCTANCE, CFLUX, FROM M/S TO MOL/M^2/S
C
      DO 700 J = 1,ICC
        DO 710 I = IL1, IL2
        IF(FCANC(I,J).GT.ZERO)THEN
C
          GB(I)=CFLUX(I)*(TFREZ/TA(I))*(PRESSG(I)/STD_PRESS)*(1/0.0224)
C         GB(I)= MIN(0.224, MAX(2.24E-03, GB(I)) )
          GB(I)= MIN(10.0, MAX(0.1, GB(I)) )
C
          IF(LEAFOPT.EQ.1)THEN
            TEMP_AN=AN_VEG(I,J)
            CO2LS(I,J) = 0.5*(CO2LS(I,J)+
     &        (CO2A(I)-( (TEMP_AN*1.37*PRESSG(I)) / GB(I)))  )
            CO2LS(I,J) = MAX (1.05*TGAMMA(I) , CO2LS(I,J))
          ELSE IF(LEAFOPT.EQ.2)THEN
            TEMP_AN=AN_SUN(I,J)
            CO2LS_SUN(I,J)= 0.5*(CO2LS_SUN(I,J)+
     &        (CO2A(I)-( (TEMP_AN*1.37*PRESSG(I))/GB(I))) )
            CO2LS_SUN(I,J) = MAX (1.05*TGAMMA(I), CO2LS_SUN(I,J))
C
            TEMP_AN=AN_SHA(I,J)
            CO2LS_SHA(I,J)= 0.5*(CO2LS_SHA(I,J)+
     &        (CO2A(I)-( (TEMP_AN*1.37*PRESSG(I))/GB(I))) )
            CO2LS_SHA(I,J) = MAX (1.05*TGAMMA(I), CO2LS_SHA(I,J))
          ENDIF
C
        ENDIF
710     CONTINUE
700   CONTINUE
C
C     FIND STOMATAL CONDUCTANCE AS PER BALL-WOODROW-BERRY FORMULATION
C     USED BY COLLATZ ET AL. OR USE THE LEUNING TYPE FORMULATION WHICH
C     USES VPD INSTEAD OF RH
C
      DO 720 J = 1,ICC
        DO 730 I = IL1, IL2
        IF(FCANC(I,J).GT.ZERO)THEN
C
C         IF LIGHT IS TOO LESS MAKE PARAMETER BB VERY SMALL
          IF(QSWV(I).LT.2.0) THEN
            USEBB(J)=0.001
          ELSE
            USEBB(J)=BB(SORT(J))
          ENDIF
          RH(I)=MAX(0.3, RH(I))          !FOLLOWING IBIS
C
          IF(LEAFOPT.EQ.1)THEN
            TEMP_AN=AN_VEG(I,J)
            IF(PS_COUP.EQ.1)THEN
              GC(I,J)=( (MM(SORT(J))*RH(I)*PRESSG(I)*TEMP_AN )/
     &          CO2LS(I,J) ) + USEBB(J)*USEAILCG(I,J)*AVE_SM_FUNC(I,J)
            ELSE IF(PS_COUP.EQ.2) THEN
              GC(I,J)=( (MM(SORT(J))*VPD_TERM(I,J)*PRESSG(I)*TEMP_AN )/
     &          (CO2LS(I,J)-TGAMMA(I)) ) +
     &          USEBB(J)*USEAILCG(I,J)*AVE_SM_FUNC(I,J)
            ENDIF
C
            GC(I,J)=MAX(GCMIN(I,J),
     &          USEBB(J)*USEAILCG(I,J)*AVE_SM_FUNC(I,J), GC(I,J))
            GC(I,J)=MIN(GCMAX(I,J), GC(I,J))
          ELSE IF(LEAFOPT.EQ.2)THEN
            TEMP_AN=AN_SUN(I,J)
            IF(PS_COUP.EQ.1)THEN
             GC_SUN(I,J) = ( (MM(SORT(J))*RH(I)*PRESSG(I)*TEMP_AN )/
     &          CO2LS_SUN(I,J) ) +
     &          USEBB(J)*AILCG_SUN(I,J)*AVE_SM_FUNC(I,J)
            ELSE IF(PS_COUP.EQ.2) THEN
             GC_SUN(I,J)=((MM(SORT(J))*VPD_TERM(I,J)*PRESSG(I)*TEMP_AN)/
     &          ( CO2LS_SUN(I,J)-TGAMMA(I) ) ) +
     &          USEBB(J)*AILCG_SUN(I,J)*AVE_SM_FUNC(I,J)
            ENDIF
C
            GC_SUN(I,J)=MAX(GCMIN(I,J),
     &        USEBB(J)*AILCG_SUN(I,J)*AVE_SM_FUNC(I,J), GC_SUN(I,J))
            GC_SUN(I,J)=MIN(GCMAX(I,J), GC_SUN(I,J))
C
            TEMP_AN=AN_SHA(I,J)
            IF(PS_COUP.EQ.1)THEN
             GC_SHA(I,J) = ( (MM(SORT(J))*RH(I)*PRESSG(I)*TEMP_AN )/
     &          CO2LS_SHA(I,J) ) +
     &          USEBB(J)*AILCG_SHA(I,J)*AVE_SM_FUNC(I,J)
            ELSE IF(PS_COUP.EQ.2) THEN
             GC_SHA(I,J)=((MM(SORT(J))*VPD_TERM(I,J)*PRESSG(I)*TEMP_AN)/
     &          (CO2LS_SHA(I,J)-TGAMMA(I)) ) +
     &          USEBB(J)*AILCG_SHA(I,J)*AVE_SM_FUNC(I,J)
            ENDIF
C
            GC_SHA(I,J)=MAX(GCMIN(I,J),
     &        USEBB(J)*AILCG_SHA(I,J)*AVE_SM_FUNC(I,J), GC_SHA(I,J))
            GC_SHA(I,J)=MIN(GCMAX(I,J), GC_SHA(I,J))
          ENDIF
C
        ENDIF
730     CONTINUE
720   CONTINUE
C
C     FIND THE INTERCELLULAR CO2 CONCENTRATION BASED ON ESTIMATED
C     VALUE OF GC
C
      DO 740 J = 1,ICC
        DO 750 I = IL1, IL2
        IF(FCANC(I,J).GT.ZERO)THEN
C
          IF(LEAFOPT.EQ.1)THEN
            TEMP_AN=AN_VEG(I,J)
            CO2I(I,J)= 0.5*(CO2I(I,J) + (CO2LS(I,J) -
     &        ( (TEMP_AN*1.65*PRESSG(I))/GC(I,J) ) ) )
            CO2I(I,J)=MAX(1.05*TGAMMA(I), MIN(CO2IMAX, CO2I(I,J)))
            PREV_CO2I(I,J)=CO2I(I,J)
          ELSE IF(LEAFOPT.EQ.2)THEN
            TEMP_AN=AN_SUN(I,J)
            CO2I_SUN(I,J)=0.5*(CO2I_SUN(I,J)+(CO2LS_SUN(I,J)-
     &        ( (TEMP_AN*1.65*PRESSG(I))/GC_SUN(I,J))) )
            CO2I_SUN(I,J)=MAX(1.05*TGAMMA(I), MIN(CO2IMAX,
     &         CO2I_SUN(I,J)))
            PREV_CO2I_SUN(I,J)=CO2I_SUN(I,J)
C
            TEMP_AN=AN_SHA(I,J)
            CO2I_SHA(I,J)=0.5*(CO2I_SHA(I,J)+(CO2LS_SHA(I,J)-
     &        ( (TEMP_AN*1.65*PRESSG(I))/GC_SHA(I,J))) )
            CO2I_SHA(I,J)=MAX(1.05*TGAMMA(I), MIN(CO2IMAX,
     &         CO2I_SHA(I,J)))
            PREV_CO2I_SHA(I,J)=CO2I_SHA(I,J)
          ENDIF
C
        ENDIF
750     CONTINUE
740   CONTINUE
C
      DO 760 J = 1,ICC
        DO 770 I = IL1, IL2
        IF(FCANC(I,J).GT.ZERO)THEN

          IF(LEAFOPT.EQ.1)THEN
            CO2I1(I,J) = PREV_CO2I(I,J)
            CO2I2(I,J) = 0.0
          ELSE IF(LEAFOPT.EQ.2)THEN
            CO2I1(I,J) = PREV_CO2I_SUN(I,J)
            CO2I2(I,J) = PREV_CO2I_SHA(I,J)
          ENDIF

        ENDIF
770     CONTINUE
760   CONTINUE
C
      IT_COUNT = IT_COUNT + 1
C
C     SEE IF WE HAVE PERFORMED THE REQUIRED NO. OF ITERATIONS, IF NOT
C     THEN WE GO BACK AND DO ANOTHER ITERATION
C
      IF(IT_COUNT.LT.REQITER)THEN
        GO TO 999
      ENDIF
C
C     WHEN REQUIRED NO. OF ITERATIONS HAVE BEEN PERFORMED THEN FIND
C     STOMATAL CONDUCTANCES FOR ALL VEGETATION TYPES IN M/S AND THEN
C     USE CONDUCTANCES TO FIND RESISTANCES. GCTU IMPLIES GC IN TRADITIONAL
C     UNITS OF M/S
C
      DO 780 J = 1,ICC
        DO 790 I = IL1, IL2
        IF(FCANC(I,J).GT.ZERO)THEN
C
          IF(LEAFOPT.EQ.1)THEN
            GCTU(I,J)=GC(I,J)*(TCAN(I)/TFREZ)*
     &        (STD_PRESS/PRESSG(I))*0.0224
            RC_VEG(I,J) = 1/GCTU(I,J)
          ELSE IF(LEAFOPT.EQ.2)THEN
            GCTU_SUN(I,J)=GC_SUN(I,J)*(TCAN(I)/TFREZ)*
     &        (STD_PRESS/PRESSG(I))*0.0224
            GCTU_SHA(I,J)=GC_SHA(I,J)*(TCAN(I)/TFREZ)*
     &        (STD_PRESS/PRESSG(I))*0.0224
C
            IF(COSZS(I).LT.0.0.OR.QSWV(I).LT.2.0)THEN
C             DON'T WANT TO REDUCE RESISTANCE AT NIGHT TO LESS THAN
C             OUR MAX. VALUE OF AROUND 5000 S/M
              GCTU(I,J)=0.5*(GCTU_SUN(I,J)+GCTU_SHA(I,J))
            ELSE
              GCTU(I,J)=GCTU_SUN(I,J)+GCTU_SHA(I,J)
            ENDIF
C
            RC_VEG(I,J) = 1/GCTU(I,J)
            AN_VEG(I,J)=AN_SUN(I,J)+AN_SHA(I,J)
            RML_VEG(I,J)=RML_SUN(I,J)+RML_SHA(I,J)
          ENDIF
C
        ENDIF
790     CONTINUE
780   CONTINUE
C
C     IF USING STORAGE LAI THEN WE SET STOMATAL RESISTANCE TO ITS
C     MAXIMUM VALUE.
C
      DO 800 J = 1, ICC
        DO 810 I = IL1, IL2
          IF(USESLAI(I,J).EQ.1.AND.AILCG(I,J).LT.0.2)THEN
            RC_VEG(I,J)=5000.0
          ENDIF
810     CONTINUE
800   CONTINUE
C
C     AND FINALLY TAKE WEIGHTED AVERAGE OF RC_VEG BASED ON FRACTIONAL
C     COVERAGE OF OUR 4 VEGETATION TYPES
C
      DO 820 J = 1,ICC
        DO 830 I = IL1, IL2
          RC(I)=RC(I)+FCANC(I,J)*RC_VEG(I,J)
830     CONTINUE
820   CONTINUE
C
      DO 840 I = IL1, IL2
        FC_TEST(I)=FCANC(I,1)+FCANC(I,2)+FCANC(I,3)+FCANC(I,4)+
     &             FCANC(I,5)+FCANC(I,6)+FCANC(I,7)+FCANC(I,8)+
     &             FCANC(I,9)
        IF(FC_TEST(I).GT.ZERO)THEN
          RC(I)=RC(I)/FC_TEST(I)
        ELSE
          RC(I) = 5000.0
        ENDIF
840   CONTINUE
C
C     FOR CROPS AND GRASS, WHENEVER LAI DROPS BELOW 1.0, CLASS 2.7 FIXES
C     LAI TO 1.0 AND REDUCES THEIR FRACTIONAL COVERAGE. TO TAKE THIS
C     INTO ACCOUNT WE REDUCE THE PER M2 OF An AND RmL WELL. IF THIS IS
C     NOT DONE THEN PHOTOSYNTHETIC UPTAKE IS OVERESTIMATED. NOTE THAT
C     THIS HAS BEEN SWITCHED OFF IN APREP FOR NOW SO THE FOLLOWING IS
C     BASICALLY A PRECAUTION.
C
      K1C = NOL2PFTS(1) + NOL2PFTS(2)
      DO 850 J = 3, 4   ! FOR CROPS AND GRASSES
       IF(J.EQ.3) THEN
         K1C = K1C + 1
       ELSE
         K1C = K1C + NOL2PFTS(J-1)
       ENDIF
       K2C = K1C + NOL2PFTS(J) - 1
       DO 855 M = K1C, K2C
        DO 860 I = IL1, IL2
          IF(FCANCMX(I,M).GT.ZERO) THEN
            AN_VEG(I,M)=AN_VEG(I,M)*(FCANC(I,M)/FCANCMX(I,M))
            RML_VEG(I,M)=RML_VEG(I,M)*(FCANC(I,M)/FCANCMX(I,M))
          ELSE
            AN_VEG(I,M)=0.0
            RML_VEG(I,M)=0.0
          ENDIF
860     CONTINUE
855    CONTINUE
850   CONTINUE
C
C     CONVERT AN_VEG AND RML_VEG TO u-MOL CO2/M2.SEC
C
      DO 870 J = 1, ICC
        DO 880 I = IL1, IL2
          AN_VEG(I,J)=AN_VEG(I,J)*1.0E+06
          RML_VEG(I,J)=RML_VEG(I,J)*1.0E+06
880     CONTINUE
870   CONTINUE
C
c
c     test by replacing TCAN with TA (restore TCAN)
c     
      DO I=IL1,IL2
       TCAN(I)=TTEMP(I)
      ENDDO
c
      RETURN
      END


