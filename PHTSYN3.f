      SUBROUTINE PHTSYN3(  AILCG,  FCANC, TCAN, CO2CONC,  PRESSG,   FC,  
     1                     CFLUX,     QA, QSWV,      IC,   THLIQ,ISAND, 
     2                        TA,   RMAT,   COSZS, XDIFFUS,  ILG,
     3                       IL1,    IL2,   IG,     ICC,   ISNOW, SLAI,
     4                       THFC, THLW, FCANCMX,  L2MAX, NOL2PFTS,
C    ---------------------- INPUTS ABOVE, OUTPUTS BELOW ---------------
     5                        RC,  CO2I1, CO2I2, AN_VEG, RML_VEG,
     6                        LFSTATUS)
C     
C               CANADIAN TERRESTRIAL ECOSYSTEM MODEL (CTEM) 
C                       PHOTOSYNTHESIS SUBROUTINE

C     HISTORY:
C
C     * FEB 27/15 - J. Melton      Increase REQUITR to 10, it was not converging correctly at 4. Also rename
C                                  WILTSM and FIELDSM to THLW and THFC, respectively, for consistency with CLASS.
C     * SEP 15/14 - J. Melton      Since SN is converted to INT, I made the assignment explicitly INT, rather
C     *                            than a real that gets cast as INT later.
C     * APR 11/14 - J. Melton      When at wilting point SM_FUNC should be 0. Was incorrectly set to 0.01
C     * JAN 15/14 - J. Melton      Fixed bug in SM_FUNC calculation that prevented
C     *                            higher number PFTs from feeling water stress.
C     *                            Also have new parameter values in Vmax, SN, and RMLcoef        
C     * AUG 16/12 - J. MELTON      Fixed GB constraint, removed incorrect
C     *                            scaling of RML, made reals being treated
C     *                            as ints more explicit, changed declaration 
C     *                            of usebb to real. Also changed use of TA to
C     *                            TCAN for RH and GB calculations to ensure 
C     *                            consistency with rest of subroutine.
C     * NOV 15, 2011 - M.LAZARE.   New version for gcm15j+ to support
C     *                            class_v3.5:
C     *                            - SAND and CLAY no longer passed in.
C     *                            - ISAND passed in so automatic array
C     *                              removed.
C     *                            - FIELDSM and WILTSM now passed in
C     *                              (calculated elsewhere in model)
C     *                              so no need for other soil propery
C     *                              automatic arrays.
C     * NOV 12, 2010 - M.LAZARE.   Revised (cosmetic) to put if condition
c     *                            on execution of loop 460 based on "fc"
c     *                            (known in tsolvc6 as "fct") since
c     *                            "qswv" (known in tsolvc6 as "qswnvc")
c     *                            is only defined for fct>0. This modification
c     *                            has been tested to give bit-for-bit the same
c     *                            answer. The model runs ok without this for
c     *                            the usual "noxlfqinitauto=on" but otherwise
C     *                            crashes without the change.  
C     * APR 23, 2010 - V.ARORA/    PREVIOUS VERSION PHTSYN2 FOR GCM15I:
C     *                M.LAZARE.   - GAMMA_W NOW 0.25 INSTEAD OF 0.425.
C     *                            - QA IS NOW INPUT WITH RH AS AN
C     *                              ALLOCATABLE ARRAY, INSTEAD OF RH
C     *                              AS INPUT AND EA,EASAT AS ALLOCATABLE
C     *                              ARRAYS. THE LATTER TWO DON'T HAVE
C     *                              TO BE ARRAYS TO CALCULATE VPD AND
C     *                              DOING IT THIS WAY IS MORE CONSISTENT
C     *                              WITH THE REST OF THE MODEL.
C     * 7 SEP.  2001 - PROGRAM TO CALCULATE STOMATAL CONDUCTANCE, TO BE
C     * V. ARORA       USED BY CLASS, BY ESTIMATING PHOTOSYNTHESIS AND THEN
C     *                RELATING PHOTOSYNTHESIS TO STOMATAL CONDUCTANCE
C     *                PREVIOUS VERSION PHTSYN).
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
C     TCAN      - CANOPY TEMPERATURE, KELVIN
C     CO2CONC   - ATMOS. CO2 IN PPM, AND THEN CONVERT IT TO
C                 PARTIAL PRESSURE, PASCALS, CO2A, FOR USE IN THIS SUBROUTINE
C     PRESSG    - ATMOS. PRESSURE, PASCALS
C     FC        - SUM OF ALL FCANC OVER A GIVEN SUB-AREA
C     CFLUX     - AERODYNAMIC CONDUCTANCE, M/S
C     QA        - SCREEN LEVEL HUMIDITY IN KG/KG
C     THLIQ     - LIQUID MOIS. CONTENT OF 3 SOIL LAYERS
C     ISAND     - SAND INDEX.
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
C     SLAI      - STORAGE LAI. THIS LAI IS USED FOR PHTSYN EVEN IF ACTUAL
C                 LAI IS ZERO. ESTIMATE OF NET PHOTOSYNTHESIS BASED ON
C                 SLAI IS USED FOR INITIATING LEAF ONSET. SEE PHENOLGY
C                 SUBROUTINE FOR MORE DETAILS.
C     THFC   - SOIL FIELD CAPACITY.
C     THLW    - SOIL WILT CAPACITY.
C     FCANCMX   - MAX. FRACTIONAL COVERAGES OF CTEM's 8 PFTs. THIS IS
C                 DIFFERENT FROM FCANC AND FCANCS (WHICH MAY VARY WITH
C                 SNOW DEPTH). FCANCMX DOESN'T CHANGE, UNLESS OF
C                 COURSE ITS CHANGED BY LAND USE CHANGE OR DYNAMIC VEGETATION.
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

      INTEGER, DIMENSION(:,:), ALLOCATABLE  :: USESLAI
      INTEGER, DIMENSION(:), ALLOCATABLE    :: SORT
      
      REAL, DIMENSION(:), ALLOCATABLE       :: FC_TEST, SIGMA,  TGAMMA,
     1     KC, KO, IPAR, GB, RH, VPD, O2_CONC, CO2A, USEBB

      REAL, DIMENSION(:,:), ALLOCATABLE     :: USEAILCG, SM_FUNC,
     1     AVE_SM_FUNC, VMAXC, JE3,SM_FUNC2,
     2     VMUNS1, VMUNS2, VMUNS3, VMUNS, VM, CO2I, PREV_CO2I, 
     3     FPAR, JC,  JC1, JC2, JC3, JE, JE1, JE2, JS, A_VEG,
     4     RC_VEG, GCTU, GCMIN, GCMAX, VPD_TERM, CO2LS, GC
C    -----------------------------------------------------------------
C                 VARIABLES USED ONLY FOR THE TWO LEAF MODEL
C
      REAL, DIMENSION(:), ALLOCATABLE :: IPAR_SUN, IPAR_SHA

      REAL, DIMENSION(:,:), ALLOCATABLE     :: GDIR, KB, FPAR_SUN,
     1    FPAR_SHA, VMAXC_SUN, VMAXC_SHA, VMUNS1_SUN, VMUNS1_SHA,
     2    VMUNS_SUN, VMUNS_SHA, VM_SUN, VM_SHA, CO2I_SUN, PREV_CO2I_SUN,
     3    CO2I_SHA, PREV_CO2I_SHA, JC1_SUN, JC1_SHA, JC3_SUN, JC3_SHA,
     4    JC_SUN, JC_SHA, JE1_SUN, JE1_SHA,
     5    JE2_SUN, JE2_SHA, JE_SUN, JE_SHA, JS_SUN, JS_SHA, 
     6    A_VEG_SUN, A_VEG_SHA, RML_SUN, RML_SHA, AN_SUN, AN_SHA,
     7    CO2LS_SUN, CO2LS_SHA, AILCG_SUN, AILCG_SHA, GC_SUN, GC_SHA,
     8    GCTU_SUN, GCTU_SHA



      INTEGER       ILG,        IC,          I,         J,        IL1, 
     1              IL2,  IT_COUNT,    REQITER,        IG,        ICC,
     2          LEAFOPT,   PS_COUP,   ISC4(KK),     ISNOW,         K1,
     3               K2,     ICOUNT, 
     4     NOL2PFTS(IC),         M,      L2MAX,         N
C
      REAL FCANC(ILG,ICC),  AILCG(ILG,ICC),    TCAN(ILG),      FC(ILG), 
     1         CFLUX(ILG),   SLAI(ILG,ICC),      QA(ILG),  INICO2I(KK) 
C
      REAL   CO2CONC(ILG),       ALPHA(KK),    OMEGA(KK),  SMSCALE(KK),
     1              TFREZ,       STD_PRESS
C
      REAL       VMAX(KK),          KN(KK),     TUP(KK),     Q10_FUNCN,
     1           TLOW(KK),        Q10_FUNC,  
     2       PRESSG(ILG), RML_VEG(ILG,ICC),            AN_VEG(ILG,ICC),
     3            QSWV(ILG),       TA(ILG),           RMAT(ILG,ICC,IG), 
     4     CO2I1(ILG,ICC),  CO2I2(ILG,ICC), 
     5                 CA,              CB,              THLIQ(ILG,IG),          
     6     THFC(ILG,IG), THLW(ILG,IG), 
     6   FCANCMX(ILG,ICC),  Q10_FUNCD
C
      REAL         MM(KK),          BB(KK),    VPD0(KK),           Q10, 
     1            RC(ILG),         CO2IMAX,    COSZS(ILG),
     2       XDIFFUS(ILG),            ZERO,               RMLCOEFF(KK)
C        
      REAL         TEMP_B,          TEMP_C,      TEMP_R,       TEMP_Q1,
     1            TEMP_Q2,         TEMP_JP,       BETA1,         BETA2,
     2            TEMP_AN
C    
      INTEGER  ISAND(ILG,IG),  SN(KK), LFSTATUS(ILG,ICC) !FLAG test LFSTATUS DEC 4 2014. JM.

      REAL use_vmax !FLAG test LFSTATUS DEC 4 2014. JM.
C
C     FOR LIMITING CO2 UPTAKE
C
      REAL CHI(KK), DELTA_CO2(ILG),  N_EFFECT(ILG), GAMMA_W, GAMMA_M
      LOGICAL SMOOTH, MIN2, MIN3
C
C
      REAL TEMP_PHI1, TEMP_PHI2, EA, EASAT, T_TEMP(ILG)
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
C                             45, 37, 37,
C                             42, 42, --,
C                             40, 50, --/
C    JM CHANGED PFT 3 TO 45 DEG FOLLOWING
C    ITO AND OIKAWA 2000. 
      DATA  TUP/307.1, 307.1, 0.000,
     &          318.1, 310.1, 310.1,
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
C     LEAF SCATTERING COEFFICIENTS, VALUES OF 0.15 & 0.17 ARE USED
C     FOR C3 AND C4 PLANTS, RESPECTIVELY
      DATA  OMEGA/0.15, 0.15, 0.00,
     &            0.15, 0.15, 0.15,
     &            0.15, 0.17, 0.00,
     &            0.15, 0.17, 0.00/
C
C     PARAMETER M USED IN PHOTOSYNTHESIS-STOMATAL CONDUCTANCE
C     COUPLING. 
C
      DATA  MM/9.0, 9.0, 0.0,
     &        12.0,12.0,12.0,
     &        12.0, 6.0, 0.0,
     &        12.0, 6.0, 0.0/
C
C     PARAMETER B USED IN PHOTOSYNTHESIS-STOMATAL CONDUCTANCE
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
C     NOT START DECREASING UNTIL SOIL MOISTURE IS ABOUT HALF WAY BETWEEN
C     WILTING POINT AND FIELD CAPACITY.
C
      DATA SN/2, 2, 0,
     &        4, 2, 2,
     &        2, 2, 0,
     &        2, 2, 0/

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
C     VALUES ARE MAINLY DERIVED FROM KATTGE ET AL. 2009 WHICH 
C     DOESN'T INCLUDE C4
      DATA VMAX/62.0E-06, 47.0E-06, 0.00E-06, 
     &          35.0E-06, 57.0E-06, 40.0E-06, !FLAG test Fri Feb27th JM, was 48 for PFT3
     &          55.0E-06, 40.0E-06, 0.00E-06,
     &          75.0E-06, 15.0E-06, 0.00E-06/

C     NEEDLE LEAF |  EVG       DCD       ---
C     BROAD LEAF  |  EVG   DCD-CLD   DCD-DRY
C     CROPS       |   C3        C4       ---
C     GRASSES     |   C3        C4       ---

C
C     NO. OF ITERATIONS FOR CALCULATING INTERCELLULAR CO2 CONCENTRATION
      DATA  REQITER/10/ 
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
      DATA  RMLCOEFF/0.015, 0.021, 0.000, 
     &               0.025, 0.015, 0.015,  
     &               0.015, 0.025, 0.000,
     &               0.013, 0.025, 0.000/
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
C      DATA GAMMA_W/0.45/  
      DATA GAMMA_W/0.17/  ! Tests occurring in Sept 2014. JM.  
C
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
C
C
      ALLOCATE(USESLAI(ILG,ICC))
      ALLOCATE(SORT(ICC))
      ALLOCATE(USEBB(ICC))

      ALLOCATE(FC_TEST(ILG))
      ALLOCATE(USEAILCG(ILG,ICC))
      ALLOCATE(SM_FUNC(ILG,IG))
      ALLOCATE(SM_FUNC2(ILG,IG))
      ALLOCATE(AVE_SM_FUNC(ILG,ICC))

      ALLOCATE(VMAXC(ILG,ICC))
      ALLOCATE(JE3(ILG,ICC))
      ALLOCATE(VMUNS1(ILG,ICC), VMUNS2(ILG,ICC), VMUNS3(ILG,ICC))
      ALLOCATE(VMUNS(ILG,ICC))
      ALLOCATE(VM(ILG,ICC))
      ALLOCATE(SIGMA(ILG))
      ALLOCATE(TGAMMA(ILG))
      ALLOCATE(KC(ILG))
      ALLOCATE(KO(ILG))
      ALLOCATE(CO2I(ILG,ICC), PREV_CO2I(ILG,ICC))
      ALLOCATE(FPAR(ILG,ICC))
      ALLOCATE(JC(ILG,ICC))
      ALLOCATE(JC1(ILG,ICC), JC2(ILG,ICC), JC3(ILG,ICC))
      ALLOCATE(JE(ILG,ICC), JE1(ILG,ICC), JE2(ILG,ICC))
      ALLOCATE(IPAR(ILG))
      ALLOCATE(JS(ILG,ICC))
      ALLOCATE(A_VEG(ILG,ICC))
      ALLOCATE(GB(ILG))
      ALLOCATE(RC_VEG(ILG,ICC))
      ALLOCATE(GCTU(ILG,ICC))
      ALLOCATE(GCMIN(ILG,ICC))
      ALLOCATE(GCMAX(ILG,ICC))
      ALLOCATE(RH(ILG))
      ALLOCATE(VPD(ILG), VPD_TERM(ILG,ICC))

      ALLOCATE(CO2LS(ILG,ICC))
      ALLOCATE(GC(ILG,ICC))
      ALLOCATE(O2_CONC(ILG))
      ALLOCATE(CO2A(ILG))

      ALLOCATE(GDIR(ILG,ICC),   KB(ILG,ICC))
      ALLOCATE(FPAR_SUN(ILG,ICC),   FPAR_SHA(ILG,ICC))
      ALLOCATE(VMAXC_SUN(ILG,ICC), VMAXC_SHA(ILG,ICC))
      ALLOCATE(VMUNS1_SUN(ILG,ICC),  VMUNS1_SHA(ILG,ICC))
      ALLOCATE(VMUNS_SUN(ILG,ICC),  VMUNS_SHA(ILG,ICC))
      ALLOCATE(VM_SUN(ILG,ICC), VM_SHA(ILG,ICC))
      ALLOCATE(CO2I_SUN(ILG,ICC), PREV_CO2I_SUN(ILG,ICC))
      ALLOCATE(CO2I_SHA(ILG,ICC), PREV_CO2I_SHA(ILG,ICC))
      ALLOCATE(JC1_SUN(ILG,ICC),    JC1_SHA(ILG,ICC))
      ALLOCATE(JC3_SUN(ILG,ICC),    JC3_SHA(ILG,ICC))
      ALLOCATE(JC_SUN(ILG,ICC),      JC_SHA(ILG,ICC))
      ALLOCATE(IPAR_SUN(ILG),       IPAR_SHA(ILG))
      ALLOCATE(JE1_SUN(ILG,ICC), JE1_SHA(ILG,ICC))
      ALLOCATE(JE2_SUN(ILG,ICC), JE2_SHA(ILG,ICC))
      ALLOCATE(JE_SUN(ILG,ICC),  JE_SHA(ILG,ICC))
      ALLOCATE(JS_SUN(ILG,ICC), JS_SHA(ILG,ICC))
      ALLOCATE(A_VEG_SUN(ILG,ICC),   A_VEG_SHA(ILG,ICC))
      ALLOCATE(RML_SUN(ILG,ICC),    RML_SHA(ILG,ICC))
      ALLOCATE(AN_SUN(ILG,ICC), AN_SHA(ILG,ICC))
      ALLOCATE(CO2LS_SUN(ILG,ICC),   CO2LS_SHA(ILG,ICC))
      ALLOCATE(AILCG_SUN(ILG,ICC),  AILCG_SHA(ILG,ICC))
      ALLOCATE(GC_SUN(ILG,ICC), GC_SHA(ILG,ICC))
      ALLOCATE(GCTU_SUN(ILG,ICC),    GCTU_SHA(ILG,ICC))

C
C     --------------------------------------------------------------
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
          SM_FUNC(I,J) = 0.0
          SM_FUNC2(I,J) = 0.0
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
C     SINCE WE DON'T YET USE TCAN, STORE THE TCAN VALUE IN A TEMP VAR AND
C     SET TCAN TO TA 
C
      T_TEMP = TCAN
      TCAN   = TA
C
C     IF LAI IS LESS THAN SLAI THAN WE USE STORAGE LAI TO PHOTOSYNTHESIZE.
C     HOWEVER, WE DO NOT USE THE STOMATAL RESISTANCE ESTIMATED IN THIS CASE,
C     BECAUSE STORAGE LAI IS AN IMAGINARY LAI, AND WE SET STOMATAL RESISTANCE
C     TO ITS MAX. NOTE THAT THE CONCEPT OF
C     STORAGE/IMAGINARY LAI IS USED FOR PHENOLOGY PURPOSES AND THIS
C     IMAGINARY LAI ACTS AS MODEL SEEDS.
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
C     FROM THE RH AND AIR TEMPERATURE WE HAVE. WE FIND E_SAT, E, AND VPD
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
          IF(TCAN(I).GE.TFREZ) THEN
              CA=17.269
              CB=35.86
          ELSE
              CA=21.874
              CB=7.66
          ENDIF
          EA     = QA(I)*PRESSG(I)/(0.622+0.378*QA(I))  
          EASAT  = 611.0*EXP(CA*(TCAN(I)-TFREZ)/(TCAN(I)-CB))
          RH(I)  = EA/EASAT
          VPD(I) = EASAT-EA
          VPD(I)=MAX(0.0,VPD(I))

420     CONTINUE
C
        K1 = 0
        DO 440 J = 1, IC
          IF(J.EQ.1) THEN
            K1 = K1 + 1
          ELSE
            K1 = K1 + NOL2PFTS(J-1)
          ENDIF
          K2 = K1 + NOL2PFTS(J) - 1
          DO 445 M = K1, K2
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
       IF(FC(I).GT.0.)                                          THEN 
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
       ENDIF
460   CONTINUE  
C
      K1 = 0
      DO 480 J = 1,IC
       IF(J.EQ.1) THEN
         K1 = K1 + 1
       ELSE
         K1 = K1 + NOL2PFTS(J-1)
       ENDIF
       K2 = K1 + NOL2PFTS(J) - 1
       DO 485 M = K1, K2
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
            AILCG_SUN(I,M)=(1.0/KB(I,M)) *
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
          FPAR(I,M)=(1.0/KN(SORT(M)))*(1.0-EXP(-KN(SORT(M))
     &               *USEAILCG(I,M)))
          IF(LEAFOPT.EQ.2)THEN
            FPAR_SUN(I,M) = ( 1.0/(KN(SORT(M))+KB(I,M)) )*
     &        ( 1.0-EXP( -1.*(KN(SORT(M))+KB(I,M))*USEAILCG(I,M) ) )
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
!         ------------- Changing Vcmax seasonally -----------------------\\\
!         Based on Wilson, Baldocchi, and Hanson, Plant, Cell and Environment 2001, 24, 571-583.
!         there is good evidence for the Vcmax varying throughout the season for deciduous tree
!         species. We are adopting a parameterization based upon their paper with some differences.
!         First, we cannot use an initial Vcmax of 0 at the start of leaf out. This is because we
!         need the leaves to gain the carbon they need for growth and we have no non-structural
!         carbohydrates to push the leaves out. So we thus start with high Vcmax. We otherwise
!         generally follow their model of spring Vcmax maximums, summer declines, then rapid
!         declines once senescence is initiated.

          ! FLAG: test only done for one-leaf model!! JM Dec 4 2014./ Jan 21 2015.
!          if (lfstatus(i,m).eq.1 .and. (m .eq. 2 .or. m .eq. 4)) then
!            use_vmax = vmax(sort(m)) * 2.0
!          else if ((lfstatus(i,m).eq.3).and.(m.eq. 2 .or. m .eq. 4))then  !leaf fall
!            use_vmax = vmax(sort(m)) * 0.5
!          else
!            use_vmax = vmax(sort(m)) !normal growth
!          end if
!          vmaxc(i,m)=use_vmax * fpar(i,m)
          VMAXC(I,M)=VMAX(SORT(M)) * FPAR(I,M)
          IF(LEAFOPT.EQ.2)THEN
             VMAXC_SUN(I,M) = VMAX(SORT(M)) * FPAR_SUN(I,M)
             VMAXC_SHA(I,M) = VMAX(SORT(M)) * FPAR_SHA(I,M)
          ENDIF
C
!         ------------- Changing Vcmax seasonally -----------------------///


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
C         ASSUMING THAT SUNLIT AND SHADED TEMPERATURES ARE SAME
C
          VMUNS2(I,M) = (1. + EXP(0.3*(TCAN(I)-TUP(SORT(M)))))
          VMUNS3(I,M) = (1. + EXP(0.3*(TLOW(SORT(M))-TCAN(I))))
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
          IF(ISAND(I,J) .EQ. -3 .OR. ISAND(I,J) .EQ. -4)THEN
            SM_FUNC(I,J)=0.0 
          ELSE ! I.E., ISAND.NE.-3 OR -4
           IF(THLIQ(I,J).LE.THLW(I,J)) THEN
            SM_FUNC(I,J)=0.0   
           ELSE IF(THLIQ(I,J).GT.THLW(I,J) .AND.
     &      THLIQ(I,J).LT.THFC(I,J)) THEN
            SM_FUNC(I,J)=(THLIQ(I,J)-THLW(I,J))/
     &          (THFC(I,J)-THLW(I,J))
           ELSE IF(THLIQ(I,J).GE.THFC(I,J)) THEN
            SM_FUNC(I,J)=1.0
           ENDIF
C
          ENDIF ! ISAND.EQ.-3 OR -4
510     CONTINUE 
500   CONTINUE 
C
      K1=0
      DO 520 J = 1, IC
       IF(J.EQ.1) THEN
         K1 = K1 + 1
       ELSE
         K1 = K1 + NOL2PFTS(J-1)
       ENDIF
       K2 = K1 + NOL2PFTS(J) - 1
       DO 525 M = K1, K2
        DO 530 I = IL1, IL2
C
         SM_FUNC2(I,1)=( 1.0 - (1.0-SM_FUNC(I,1))**SN(SORT(M)) )
         SM_FUNC2(I,2)=( 1.0 - (1.0-SM_FUNC(I,2))**SN(SORT(M)) )
         SM_FUNC2(I,3)=( 1.0 - (1.0-SM_FUNC(I,3))**SN(SORT(M)) )

         SM_FUNC2(I,1)=SM_FUNC2(I,1)+(1.0-SM_FUNC2(I,1))
     &                 *SMSCALE(SORT(M))
         SM_FUNC2(I,2)=SM_FUNC2(I,2)+(1.0-SM_FUNC2(I,2))
     &                 *SMSCALE(SORT(M))
         SM_FUNC2(I,3)=SM_FUNC2(I,3)+(1.0-SM_FUNC2(I,3))
     &                 *SMSCALE(SORT(M))

         AVE_SM_FUNC(I,M)=SM_FUNC2(I,1)*RMAT(I,M,1) +
     &                    SM_FUNC2(I,2)*RMAT(I,M,2) +
     &                    SM_FUNC2(I,3)*RMAT(I,M,3)

         AVE_SM_FUNC(I,M)= AVE_SM_FUNC(I,M) /
     &    (RMAT(I,M,1)+RMAT(I,M,2)+RMAT(I,M,3))
C
C        MAKE SURE THAT I DO NOT USE SM_FUNC FROM 2nd AND 3rd LAYERS
C        IF THEY ARE BED ROCK
C
         IF(ISAND(I,3) .EQ. -3)THEN ! THIRD LAYER BED ROCK
           AVE_SM_FUNC(I,M)=SM_FUNC2(I,1)*RMAT(I,M,1) +
     &                      SM_FUNC2(I,2)*RMAT(I,M,2)
C
           AVE_SM_FUNC(I,M)= AVE_SM_FUNC(I,M) /
     &      (RMAT(I,M,1)+RMAT(I,M,2))
         ENDIF
C
         IF(ISAND(I,2) .EQ. -3)THEN ! SECOND LAYER BED ROCK
           AVE_SM_FUNC(I,M)=SM_FUNC2(I,1)
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
        TGAMMA(I) = O2_CONC(I) / (2.0*SIGMA(I))
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
            JC2(I,J) = KC(I) * (1.0 + (O2_CONC(I)/KO(I)) )
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
            JC2(I,J) = KC(I) * (1.0 + (O2_CONC(I)/KO(I)) )
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
            JE1(I,J)= FPAR(I,J)*ALPHA(SORT(J))*(1.0-OMEGA(SORT(J)))
            JE2(I,J)=( CO2I(I,J)-TGAMMA(I) ) /
     &        ( CO2I(I,J)+ (2.0*TGAMMA(I)) )
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
     &        ( CO2I_SUN(I,J)+ (2.0*TGAMMA(I)) )
            JE2_SHA(I,J)=( CO2I_SHA(I,J)-TGAMMA(I) ) /
     &        ( CO2I_SHA(I,J)+ (2.0*TGAMMA(I)) )

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
      DO 641 I = IL1, IL2
        DELTA_CO2(I) = CO2CONC(I) - 288.0
        TEMP_R = LOG(1.0 + (DELTA_CO2(I)/288.0))
        TEMP_B = 1.0 + TEMP_R * (GAMMA_W)
        TEMP_C = 1.0 + TEMP_R * (GAMMA_M)
        N_EFFECT(I) = TEMP_B / TEMP_C
C       LIMIT N_EFFECT TO MAX OF 1.0 SO THAT NO UP-REGULATION OCCURS
        N_EFFECT(I) = MIN(1.0, N_EFFECT(I)) 
641   CONTINUE
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
              TEMP_Q1 = (TEMP_B + SQRT(TEMP_R)) / (2. * BETA1)
              TEMP_Q2 = (TEMP_B - SQRT(TEMP_R)) / (2. * BETA1)
              TEMP_JP = MIN(TEMP_Q1, TEMP_Q2)
C
              TEMP_B  = TEMP_JP + JS(I,J)
              TEMP_C  = TEMP_JP * JS(I,J)
              TEMP_R  = MAX( (TEMP_B**2 - 4. * BETA2 * TEMP_C), 0.0)
              TEMP_Q1 = (TEMP_B + SQRT(TEMP_R)) / (2. * BETA2)
              TEMP_Q2 = (TEMP_B - SQRT(TEMP_R)) / (2. * BETA2)
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
               TEMP_Q1 = (TEMP_B + SQRT(TEMP_R)) / (2. * BETA1)
               TEMP_Q2 = (TEMP_B - SQRT(TEMP_R)) / (2. * BETA1)
               TEMP_JP = MIN(TEMP_Q1, TEMP_Q2)

               TEMP_B  = TEMP_JP + JS_SUN(I,J)
               TEMP_C  = TEMP_JP * JS_SUN(I,J)
               TEMP_R  = MAX( (TEMP_B**2 - 4. * BETA2 * TEMP_C), 0.0)
               TEMP_Q1 = (TEMP_B + SQRT(TEMP_R)) / (2. * BETA2)
               TEMP_Q2 = (TEMP_B - SQRT(TEMP_R)) / (2. * BETA2)
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
               TEMP_Q1 = (TEMP_B + SQRT(TEMP_R)) / (2. * BETA1)
               TEMP_Q2 = (TEMP_B - SQRT(TEMP_R)) / (2. * BETA1)
               TEMP_JP = MIN(TEMP_Q1, TEMP_Q2)
C
               TEMP_B  = TEMP_JP + JS_SHA(I,J)
               TEMP_C  = TEMP_JP * JS_SHA(I,J)
               TEMP_R  = MAX( (TEMP_B**2 - 4. * BETA2 * TEMP_C), 0.0)
               TEMP_Q1 = (TEMP_B + SQRT(TEMP_R)) / (2. * BETA2)
               TEMP_Q2 = (TEMP_B - SQRT(TEMP_R)) / (2. * BETA2)
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
            RML_VEG(I,J) = RMLCOEFF(SORT(J))*VMAXC(I,J)*Q10_FUNCN 
           ENDIF
           AN_VEG(I,J) = A_VEG(I,J) - RML_VEG(I,J)
          ELSE IF(LEAFOPT.EQ.2)THEN
           IF(COSZS(I).GT.0.0)THEN
            RML_SUN(I,J) = RMLCOEFF(SORT(J))*VMAXC_SUN(I,J)*Q10_FUNCD
            RML_SHA(I,J) = RMLCOEFF(SORT(J))*VMAXC_SHA(I,J)*Q10_FUNCD
           ELSE
            RML_SUN(I,J)=RMLCOEFF(SORT(J))*VMAXC_SUN(I,J)*Q10_FUNCN 
            RML_SHA(I,J)=RMLCOEFF(SORT(J))*VMAXC_SHA(I,J)*Q10_FUNCN 
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
          GB(I)=CFLUX(I)*(TFREZ/TCAN(I))*(PRESSG(I)/STD_PRESS)
     &          *(1./0.0224)
C
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
            RC_VEG(I,J) = 1./GCTU(I,J)
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
            RC_VEG(I,J) = 1./GCTU(I,J)
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
C     CONVERT AN_VEG AND RML_VEG TO u-MOL CO2/M2.SEC
C
      DO 870 J = 1, ICC
        DO 880 I = IL1, IL2
          AN_VEG(I,J)=AN_VEG(I,J)*1.0E+06
          RML_VEG(I,J)=RML_VEG(I,J)*1.0E+06
880     CONTINUE
870   CONTINUE
C
C     PUT TCAN BACK TO ITS ORIGINAL VALUE FROM THE TEMP VAR 
      TCAN = T_TEMP


      DEALLOCATE(USESLAI)
      DEALLOCATE(SORT)
      DEALLOCATE(USEBB)

      DEALLOCATE(FC_TEST)
      DEALLOCATE(USEAILCG)
      DEALLOCATE(SM_FUNC)
      DEALLOCATE(SM_FUNC2)
      DEALLOCATE(AVE_SM_FUNC)

      DEALLOCATE(VMAXC)
      DEALLOCATE(JE3)
      DEALLOCATE(VMUNS1, VMUNS2, VMUNS3)
      DEALLOCATE(VMUNS)
      DEALLOCATE(VM)
      DEALLOCATE(SIGMA)
      DEALLOCATE(TGAMMA)
      DEALLOCATE(KC)
      DEALLOCATE(KO)
      DEALLOCATE(CO2I, PREV_CO2I)
      DEALLOCATE(FPAR)
      DEALLOCATE(JC)
      DEALLOCATE(JC1, JC2, JC3)
      DEALLOCATE(JE, JE1, JE2)
      DEALLOCATE(IPAR)
      DEALLOCATE(JS)
      DEALLOCATE(A_VEG)
      DEALLOCATE(GB)
      DEALLOCATE(RC_VEG)
      DEALLOCATE(GCTU)
      DEALLOCATE(GCMIN)
      DEALLOCATE(GCMAX)
      DEALLOCATE(RH)

      DEALLOCATE(VPD, VPD_TERM)

      DEALLOCATE(CO2LS)
      DEALLOCATE(GC)
      DEALLOCATE(O2_CONC)
      DEALLOCATE(CO2A)

      DEALLOCATE(GDIR,   KB)
      DEALLOCATE(FPAR_SUN,   FPAR_SHA)
      DEALLOCATE(VMAXC_SUN, VMAXC_SHA)
      DEALLOCATE(VMUNS1_SUN,  VMUNS1_SHA)
      DEALLOCATE(VMUNS_SUN,  VMUNS_SHA)
      DEALLOCATE(VM_SUN, VM_SHA)
      DEALLOCATE(CO2I_SUN, PREV_CO2I_SUN)
      DEALLOCATE(CO2I_SHA, PREV_CO2I_SHA)
      DEALLOCATE(JC1_SUN,    JC1_SHA)
      DEALLOCATE(JC3_SUN,    JC3_SHA)
      DEALLOCATE(JC_SUN,      JC_SHA)
      DEALLOCATE(IPAR_SUN,       IPAR_SHA)
      DEALLOCATE(JE1_SUN, JE1_SHA)
      DEALLOCATE(JE2_SUN, JE2_SHA)
      DEALLOCATE(JE_SUN,  JE_SHA)
      DEALLOCATE(JS_SUN, JS_SHA)
      DEALLOCATE(A_VEG_SUN,   A_VEG_SHA)
      DEALLOCATE(RML_SUN,    RML_SHA)
      DEALLOCATE(AN_SUN, AN_SHA)
      DEALLOCATE(CO2LS_SUN,   CO2LS_SHA)
      DEALLOCATE(AILCG_SUN,  AILCG_SHA)
      DEALLOCATE(GC_SUN, GC_SHA)
      DEALLOCATE(GCTU_SUN, GCTU_SHA)


      RETURN
      END


