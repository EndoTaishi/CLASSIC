      SUBROUTINE PHENOLGY(GLEAFMAS, BLEAFMAS,      ICC,       IG,      
     1                         ILG,      IL1,      IL2,     TBAR,    
     2                       THLIQ,   WILTSM,  FIELDSM,       TA,  
     3                       ANVEG,     IDAY,     RADL, ROOTTEMP,
     4                    RMATCTEM, STEMMASS, ROOTMASS,     SORT,
     5                       L2MAX, NOL2PFTS,       IC,  
C    6 ------------------ INPUTS ABOVE THIS LINE ----------------------   
     7                    FLHRLOSS, LEAFLITR, LFSTATUS,  PANDAYS,
     8                    COLDDAYS)  
C    9 --- VARIABLES WHICH ARE UPDATED AND OUTPUTS ABOVE THIS LINE ----
C
C               CANADIAN TERRESTRIAL ECOSYSTEM MODEL (CTEM) V1.0
C               PHENOLOGY, LEAF TURNOVER & MORTALITY SUBROUTINE
C
C     15  APR. 2003 - THIS SUBROUTINE CALCULATES THE LEAF STATUS 
C     V. ARORA        FOR CTEM's PFTs AND LEAF LITTER GENERATED BY
C                     NORMAL TURNOVER OF LEAVES AND COLD AND DROUGHT 
C                     STRESS. CROP HARVEST IS ALSO MODELLED IN THIS
C                     SUBROUTINE, AND FOR GRASSES GREEN LEAVES ARE
C                     CONVERTED INTO BROWN.
C
C     INPUTS 
C
C     GLEAFMAS  - GREEN OR LIVE LEAF MASS IN KG C/M2, FOR THE 9 PFTs
C     BLEAFMAS  - BROWN OR DEAD LEAF MASS IN KG C/M2, FOR THE 9 PFTs
C     ICC       - NO. OF CTEM PLANT FUNCTION TYPES, CURRENTLY 9
C     IG        - NO. OF SOIL LAYERS (CURRENTLY 3)
C     ILG       - NO. OF GRID CELLS IN LATITUDE CIRCLE
C     IL1,IL2   - IL1=1, IL2=ILG
C     TBAR      - SOIL TEMPERATURE, K
C     THLIQ     - LIQUID SOIL MOISTURE CONTENT IN 3 SOIL LAYERS
C     WILTSM    - WILTING POINT SOIL MOISTURE CONTENT
C     FIELDSM   - FIELD CAPACITY SOIL MOISTURE CONTENT
C                 BOTH CALCULATED IN ALLOCATE SUBROUTINE
C     TA        - AIR TEMPERATURE, K
C     ANVEG     - NET PHOTOSYNTHESIS RATE OF CTEM's PFTs, uMOL CO2/M2.S
C     IDAY      - DAY OF YEAR
C     RADL      - LATITUDE IN RADIANS
C     ROOTTEMP  - ROOT TEMPERATURE, WHICH IS A FUNCTION OF SOIL TEMPERATURE
C                 OF COURSE, K.
C     RMATCTEM  - FRACTION OF ROOTS IN EACH SOIL LAYER FOR EACH PFT
C     STEMMASS  - STEM MASS FOR EACH OF THE 9 CTEM PFTs, Kg C/M2
C     ROOTMASS  - ROOT MASS FOR EACH OF THE 9 CTEM PFTs, Kg C/M2
C     SORT      - INDEX FOR CORRESPONDENCE BETWEEN 9 PFTs AND THE
C                 12 VALUES IN PARAMETERS VECTORS
C     L2MAX     - MAXIMUM NUMBER OF LEVEL 2 CTEM PFTs
C     NOL2PFTS  - NUMBER OF LEVEL 2 CTEM PFTs
C     IC        - NUMBER OF CLASS PFTs
C
C     UPDATES
C
C     FLHRLOSS  - FALL & HARVEST LOSS FOR BDL DCD PLANTS AND CROPS, 
C                 RESPECTIVELY, Kg C/M2.
C     PANDAYS   - COUNTER FOR POSITIVE NET PHOTOSYNTHESIS (An) DAYS FOR 
C                 INITIATING LEAF ONSET
C     LFSTATUS  - INTEGER INDICATING LEAF STATUS OR MODE
C                 1 - MAX. GROWTH OR ONSET, WHEN ALL NPP IS ALLOCATED TO
C                     LEAVES
C                 2 - NORMAL GROWTH, WHEN NPP IS ALLOCATED TO LEAVES, STEM,
C                     AND ROOT
C                 3 - FALL FOR DCD TREES/HARVEST FOR CROPS, WHEN ALLOCATION
C                     TO LEAVES IS ZERO.
C                 4 - NO LEAVES
C     COLDDAYS  - COLD DAYS COUNTER FOR TRACKING DAYS BELOW A CERTAIN 
C                 TEMPERATURE THRESHOLD FOR NDL DCD AND CROP PFTs.
C
C     OUTPUTS
C
C     LEAFLITR  - LEAF LITTER GENERATED BY NORMAL TURNOVER, COLD AND
C                 DROUGHT STRESS, AND LEAF FALL/HARVEST, Kg C/M2
C
C     ------------------------------------------------------------------    
C
      IMPLICIT NONE
C
      INTEGER ILG, ICC, IG, IL1, IL2, I, J, K, IDAY, KK, N, M, IC, K1,
     1         K2
C
      PARAMETER (KK=12)  ! PRODUCT OF CLASS PFTs AND L2MAX
C
      INTEGER        SORT(ICC),     NOL2PFTS(IC),                L2MAX
C
      REAL  GLEAFMAS(ILG,ICC), BLEAFMAS(ILG,ICC),              TA(ILG),
     1           TBAR(ILG,IG),     THLIQ(ILG,IG),         SAND(ILG,IG), 
     2           CLAY(ILG,IG),    ANVEG(ILG,ICC),    LEAFLITR(ILG,ICC),
     3      ROOTTEMP(ILG,ICC),                    RMATCTEM(ILG,ICC,IG),
     4      STEMMASS(ILG,ICC), ROOTMASS(ILG,ICC)
C
C
      REAL     LFESPANY(KK),       DRLSRTMX(KK),          CDLSRTMX(KK), 
     1            DRGTA(KK),          COLDA(KK),                  ZERO, 
     2         LWRTHRSH(KK),       HARVTHRS(KK),           COLDTHRS(2), 
     3          FLHRSPAN(2),       THRPRCNT(KK),               ETA(KK),
     4            KAPPA(KK),             KN(KK),           SPECSLA(KK),
     5             ROOTHRSH

      INTEGER  DAYSCHK(KK),    PANDAYS(ILG,ICC),     LFSTATUS(ILG,ICC),
     1    CHKMODE(ILG,ICC),     COLDDAYS(ILG,2),            COLDLMT(2)
C
      REAL        SLA(ICC),      AILCG(ILG,ICC),        AILCB(ILG,ICC)
     
C
      REAL                      FIELDSM(ILG,IG),        WILTSM(ILG,IG),
     1                DAY,            RADL(ILG),                 THETA,
     2              DECLI,                 TERM,         DAYLNGTH(ILG),
     3  NRMLLOSS(ILG,ICC),     BETADRGT(ILG,IG),     DRGTSTRS(ILG,ICC),
     4  DRGTLSRT(ILG,ICC),    DRGTLOSS(ILG,ICC),     COLDLOSS(ILG,ICC),
     5  COLDSTRS(ILG,ICC),    COLDLSRT(ILG,ICC),     FLHRLOSS(ILG,ICC),
     6    LFTHRS(ILG,ICC),             FRACBOFG
C
      COMMON /CTEM1/ ETA, KAPPA, KN
      COMMON /CTEM2/ LFESPANY, FRACBOFG, SPECSLA
C
C     ------------------------------------------------------------------
C                     CONSTANTS AND PARAMETERS USED 
C
C     ALSO NOTE THE STRUCTURE OF PARAMETER VECTORS WHICH CLEARLY 
C     SHOWS THE CLASS PFTs (ALONG ROWS) AND CTEM SUB-PFTs (ALONG 
C     COLUMNS)
C
C     NEEDLE LEAF |  EVG       DCD       ---
C     BROAD LEAF  |  EVG   DCD-CLD   DCD-DRY
C     CROPS       |   C3        C4       ---
C     GRASSES     |   C3        C4       ---
C
C     MAX. LOSS RATE FOR COLD STRESS FOR ALL 8 PFTs, (1/DAY)
      DATA CDLSRTMX/0.30, 0.30, 0.00,
c      DATA CDLSRTMX/0.15, 0.30, 0.00,
     &              0.30, 0.15, 0.15,
     &              0.15, 0.15, 0.00,
     &              0.15, 0.15, 0.00/
C
C     MAX. LOSS RATE FOR DROUGHT STRESS FOR ALL 9 PFTs, (1/DAY)
      DATA DRLSRTMX/0.005, 0.005, 0.000,
     &              0.005, 0.005, 0.025,
     &              0.005, 0.005, 0.000,
     &              0.050, 0.050, 0.000/    
C
C     PARAMETER DETERMINING HOW FAST SOIL DRYNESS CAUSES LEAVES TO FALL
      DATA DRGTA/3.0, 3.0, 0.0,
     &           3.0, 3.0, 3.0,
     &           3.0, 3.0, 0.0,
     &           3.0, 3.0, 0.0/
C
C     PARAMETER DETERMINING HOW FAST COLD TEMPERATURES CAUSES LEAVES TO FALL
      DATA COLDA/3.0, 3.0, 0.0,
c      DATA COLDA/1.5, 3.0, 0.0,
     &           3.0, 3.0, 3.0,
     &           3.0, 3.0, 0.0,
     &           3.0, 3.0, 0.0/
C
C     LOWER TEMPERATURE THRESHOLD FOR CTEM's 9 PFTs. THESE ARE USED
C     TO ESTIMATE COLD STRESS RELATED LEAF LOSS RATE (DEGREE C)
      DATA LWRTHRSH/-45.0, -5.0, 0.0,
c      DATA LWRTHRSH/-25.0, -5.0, 0.0,
     &                5.0,  5.0, 5.0,
     &                5.0,  5.0, 0.0,
     &                0.1,  5.0, 0.0/
C
C     NUMBER OF DAYS OVER WHICH TO CHECK IF NET PHOTOSYNTHETIC RATE IS
C     POSITIVE BEFORE INITIATING LEAF ONSET
      DATA DAYSCHK/7, 7, 0,
     &             7, 7, 7,
     &             7, 7, 0,
     &             7, 7, 0/
C
C     NO. OF DAYS FOR WHICH SOME TEMPERATURE HAS TO REMAIN BELOW 
C     A GIVEN THRESHOLD FOR INITIATING A PROCESS
C     1. -5 C THRESHOLD FOR INITIATING "LEAF FALL" MODE FOR NDL DCD TREES
C     2.  8 C THRESHOLD FOR INITIATING "HARVEST" MODE FOR CROPS
C     THE ARRAY COLDDAYS(I,2) TRACKS DAYS CORRESPONDING TO THESE THRESHOLDS
      DATA COLDLMT /7,5/   ! DAYS
      DATA COLDTHRS/-5.0,8.0/
C
C     LAI THRESHOLD FOR HARVESTING CROPS. VALUES ARE ZERO FOR ALL PFTs
C     OTHER THAN C3 AND C4 CROPS.
      DATA HARVTHRS/0.0, 0.0, 0.0,
     &              0.0, 0.0, 0.0,
     &              4.5, 3.5, 0.0,
     &              0.0, 0.0, 0.0/
C
C     HARVEST SPAN (TIME IN DAYS OVER WHICH CROPS ARE HARVESTED, ~15 DAYS), 
C     AND  FALL SPAN (TIME IN DAYS OVER WHICH BDL COLD DCD PLANTS SHED 
C     THEIR LEAVES, ~30 DAYS)
      DATA FLHRSPAN/17.0, 45.0/
C
C     PERCENTAGE OF MAX. LAI THAT CAN BE SUPPORTED WHICH IS USED AS A
C     THRESHOLD FOR DETERMINING LEAF PHENOLOGY STATUS
      DATA THRPRCNT/40.0, 40.0,  0.0,
     &              40.0, 50.0, 50.0,
     &              50.0, 50.0,  0.0,
     &              40.0, 40.0,  0.0/
C
C     ROOT TEMPERATURE THRESHOLD FOR INITIATING LEAF ONSET FOR COLD
C     BROADLEAF DECIDUOUS PFT, DEGREES CELCIUS
c      DATA ROOTHRSH/2.0/
      DATA ROOTHRSH/15.0/
C
C     ZERO
      DATA ZERO/1E-20/

      REAL, PARAMETER :: PI=3.1415926535898d0
C
C     ---------------------------------------------------------------
C
      IF(ICC.NE.9)                            CALL XIT('PHENOLGY',-1)
C
C     INITIALIZE REQUIRED ARRAYS TO ZERO

      DO 110 J = 1, ICC
        SLA(J)=0.0                ! SPECIFIC LEAF AREA
110   CONTINUE
C
      DO 120 J = 1, IG
        DO 130 I = IL1, IL2
          BETADRGT(I,J)=0.0       ! (1 - DROUGHT STRESS)
130     CONTINUE
120   CONTINUE
C
      DO 140 J = 1,ICC
        DO 150 I = IL1, IL2
          AILCG(I,J)=0.0          ! GREEN LAI
          AILCB(I,J)=0.0          ! BROWN LAI
          CHKMODE(I,J)=0          ! INDICATOR FOR MAKING SURE THAT LEAF STATUS 
C                                 ! IS UPDATED
          LEAFLITR(I,J)=0.0       ! LEAF LITTER
          NRMLLOSS(I,J)=0.0       ! LEAF LOSS DUE TO NORMAL TURNOVER
          DRGTSTRS(I,J)=0.0       ! DROUGHT STRESS TERM
          DRGTLSRT(I,J)=0.0       ! DROUGHT LOSS RATE
          DRGTLOSS(I,J)=0.0       ! LEAF LOSS DUE TO DROUGHT STRESS
          COLDSTRS(I,J)=0.0       ! COLD STRESS TERM
          COLDLSRT(I,J)=0.0       ! COLD LOSS RATE
          COLDLOSS(I,J)=0.0       ! LEAF LOSS DUE TO COLD STRESS
          LFTHRS(I,J)=0.0         ! THRESHOLD LAI FOR FINDING LEAF STATUS
150     CONTINUE                  
140   CONTINUE
C
C     INITIALIZATION ENDS    
C
C     ------------------------------------------------------------------
C
C     CONVERT GREEN LEAF MASS INTO LEAF AREA INDEX USING SPECIFIC LEAF
C     AREA (SLA, M2/KG C) ESTIMATED USING LEAF LIFE SPAN. SEE BIO2STR
C     SUBROUTINE FOR MORE DETAILS. 
C
      DO 160 J = 1,ICC
        SLA(J) = 25.0*(LFESPANY(SORT(J))**(-0.50))
        IF(SPECSLA(SORT(J)).GT.ZERO) SLA(J)=SPECSLA(SORT(J))
160   CONTINUE
C
      DO 170 J = 1,ICC
        N = SORT(J)
        DO 180 I = IL1,IL2
          AILCG(I,J)=SLA(J)*GLEAFMAS(I,J)
          AILCB(I,J)=SLA(J)*BLEAFMAS(I,J)*FRACBOFG
C         ALSO FIND THRESHOLD LAI AS A FUNCTION OF STEM+ROOT BIOMASS
C         WHICH IS USED TO DETERMINE LEAF STATUS
          LFTHRS(I,J)=((STEMMASS(I,J)+ROOTMASS(I,J))/ETA(N))
     &     **(1/KAPPA(N))   
          LFTHRS(I,J)=(THRPRCNT(N)/100.0)*SLA(J)*LFTHRS(I,J)
180     CONTINUE
170   CONTINUE
C
C     USING GREEN LEAF AREA INDEX (AILCG) DETERMINE THE LEAF STATUS FOR
C     EACH PFT. LOOPS 190 AND 200 THUS INITIALIZE LFSTATUS, IF THIS
C     THIS INFORMATION IS NOT PASSED SPECIFICALLY AS AN INITIALIZATION
C     QUANTITY. 
C
      DO 190 J = 1, ICC
        DO 200 I = IL1, IL2
          IF(LFSTATUS(I,J).EQ.0)THEN
            IF(AILCG(I,J).LE.ZERO)THEN            
              LFSTATUS(I,J)=4                      !NO LEAVES
            ELSE IF (AILCG(I,J).GT.LFTHRS(I,J))THEN 
              LFSTATUS(I,J)=2                      !NORMAL GROWTH
            ELSE                                  
              LFSTATUS(I,J)=4                      !TREAT THIS AS NO LEAVES
            ENDIF                                  !SO THAT WE START GROWING 
          ENDIF                                    !IF POSSIBLE
200     CONTINUE
190   CONTINUE
C
C     KNOWING LFSTATUS (AFTER INITIALIZATION ABOVE OR USING VALUE FROM
C     FROM PREVIOUS TIME STEP) WE DECIDE IF WE STAY IN A GIVEN LEAF
C     MODE OR WE MOVE TO SOME OTHER MODE.
C
C     WE START WITH THE "NO LEAVES" MODE
C     ----------------------------------
C
C     ADD ONE TO PANDAYS(I,J) IF DAILY An IS POSITIVE, OTHERWISE SET IT TO
C     ZERO.
C
      DO 220 J = 1, ICC
        N = SORT(J)
        DO 230 I = IL1, IL2
          IF(ANVEG(I,J).GT.ZERO) THEN
            PANDAYS(I,J)=PANDAYS(I,J)+1
            IF(PANDAYS(I,J).GT.DAYSCHK(N))THEN
              PANDAYS(I,J)=DAYSCHK(N)
            ENDIF
          ELSE
            PANDAYS(I,J)=0
          ENDIF
230     CONTINUE
220   CONTINUE
C
C     IF IN "NO LEAVES" MODE CHECK IF An HAS BEEN POSITIVE OVER LAST
C     DAYSCHK(J) DAYS TO MOVE INTO "MAX. GROWTH" MODE. IF NOT WE STAY
C     IN "NO LEAVES" MODE. ALSO SET THE CHKMODE(I,J) SWITCH TO 1.
C
      DO 240 J = 1, ICC
        N = SORT(J)
        DO 250 I = IL1, IL2
          IF(CHKMODE(I,J).EQ.0.AND.LFSTATUS(I,J).EQ.4)THEN
            IF(PANDAYS(I,J).GE.DAYSCHK(N))THEN
              LFSTATUS(I,J)=1        ! SWITCH TO "MAX. GROWTH" MODE
              CHKMODE(I,J)=1         ! MODE CHECKED, NO MORE CHECKS FURTHER DOWN
            ELSE
              LFSTATUS(I,J)=4        ! STAY IN "NO LEAVES" MODE
              CHKMODE(I,J)=1         ! MODE CHECKED, NO MORE CHECKS FURTHER DOWN
            ENDIF
          ENDIF
250     CONTINUE
240   CONTINUE
C
C     FIND DAY LENGTH USING DAY OF YEAR AND LATITUDE. THIS IS TO BE USED FOR
C     INITIATING LEAF OFFSET FOR BROAD LEAF DCD TREES.
C
CRL L      DAY=FLOAT(IDAY)
      DAY=REAL(IDAY)
      DO 260 I = IL1, IL2
       THETA=0.2163108 + 2*ATAN(0.9671396*TAN(0.0086*(DAY-186.0)))
       DECLI=ASIN(0.39795*COS(THETA))    !DECLINATION
       TERM=(SIN(RADL(I))*SIN(DECLI))/(COS(RADL(I))*COS(DECLI))
       TERM=MAX(-1.0,MIN(TERM,1.0))
       DAYLNGTH(I)=24.0-(24.0/PI)*ACOS(TERM)
260   CONTINUE
C
C     EVEN IF PANDAYS CRITERIA HAS BEEN SATISFIED DO NOT GO INTO MAX. GROWTH
C     MODE IF ENVIRONMENTAL CONDITIONS ARE SUCH THAT THEY WILL FORCE LEAF FALL
C     OR HARVEST MODE.
C
      DO 270 I = IL1, IL2
C       NEEDLE LEAF DCD
        IF(LFSTATUS(I,2).EQ.1.AND.CHKMODE(I,2).EQ.1)THEN
          IF(TA(I).LT.(COLDTHRS(1)+273.16))THEN
            LFSTATUS(I,2)=4
          ENDIF
        ENDIF
C
C       BROAD LEAF DCD CLD & DRY
        IF(LFSTATUS(I,4).EQ.1.AND.CHKMODE(I,4).EQ.1)THEN
          IF(ROOTTEMP(I,4).LT.(ROOTHRSH+273.16).OR.
     &    (DAYLNGTH(I).LT.11.0.AND.ROOTTEMP(I,4).LT.(11.15+273.16)))THEN
            LFSTATUS(I,4)=4
          ENDIF
        ENDIF

        IF(LFSTATUS(I,5).EQ.1.AND.CHKMODE(I,5).EQ.1)THEN
          IF(ROOTTEMP(I,5).LT.(ROOTHRSH+273.16).OR.
     &    (DAYLNGTH(I).LT.11.0.AND.ROOTTEMP(I,5).LT.(11.15+273.16)))THEN
            LFSTATUS(I,5)=4
          ENDIF
        ENDIF
C
C       CROPS
        IF(LFSTATUS(I,6).EQ.1.AND.CHKMODE(I,6).EQ.1)THEN
          IF(TA(I).LT.(COLDTHRS(2)+273.16))THEN
            LFSTATUS(I,6)=4
          ENDIF
        ENDIF
        IF(LFSTATUS(I,7).EQ.1.AND.CHKMODE(I,7).EQ.1)THEN
          IF(TA(I).LT.(COLDTHRS(2)+273.16))THEN
            LFSTATUS(I,7)=4
          ENDIF
        ENDIF
270   CONTINUE
C
C     SIMILAR TO THE WAY WE COUNT NO. OF DAYS WHEN An IS POSITIVE, WE
C     FIND NO. OF DAYS WHEN TEMPERATURE IS BELOW -5 C. WE NEED THIS TO
C     DETERMINE IF WE GO INTO "LEAF FALL" MODE FOR NEEDLE LEAF DCD TREES.
C
C     ALSO ESTIMATE NO. OF DAYS BELOW 8 C. WE USE THESE DAYS TO DECIDE 
C     IF ITS COLD ENOUGH TO HARVEST CROPS.
C
      DO 280 K = 1, 2
        DO 290 I = IL1, IL2
          IF(TA(I).LT.(COLDTHRS(K)+273.16)) THEN
            COLDDAYS(I,K)=COLDDAYS(I,K)+1
            IF(COLDDAYS(I,K).GT.COLDLMT(K))THEN
              COLDDAYS(I,K)=COLDLMT(K)
            ENDIF
          ELSE
            COLDDAYS(I,K)=0
          ENDIF
290     CONTINUE
280   CONTINUE
C
C     IF IN "MAX GROWTH" MODE
C     ------------------------
C
C     IF MODE HASN'T BEEN CHECKED AND WE ARE IN "MAX. GROWTH" MODE, THEN 
C     CHECK IF WE ARE ABOVE PFT-DEPENDENT LAI THRESHOLD. IF LAI IS MORE
C     THEN THIS THRESHOLD WE MOVE INTO "NORMAL GROWTH" MODE, OTHERWISE
C     WE STAY IN "MAX GROWTH" MODE SO THAT LEAVES CAN GROW AT THEIR
C     MAX. CLIMATE-DEPENDENT RATE
C
      DO 300 J = 1, ICC
        DO 310 I = IL1, IL2
          IF(CHKMODE(I,J).EQ.0.AND.LFSTATUS(I,J).EQ.1)THEN
            IF(AILCG(I,J).GE.LFTHRS(I,J))THEN
              LFSTATUS(I,J)=2        ! SWITCH TO "NORMAL GROWTH" MODE
              CHKMODE(I,J)=1         
            ELSE IF(AILCG(I,J).LE.ZERO) THEN
              LFSTATUS(I,J)=4        ! SWITCH TO "NO LEAVES" MODE
              CHKMODE(I,J)=1         
              PANDAYS(I,J)=0
            ELSE 
              LFSTATUS(I,J)=1        ! STAY IN "MAX. GROWTH" MODE
              CHKMODE(I,J)=1         
            ENDIF
C
C           FOR DCD TREES WE ALSO NEED TO GO INTO "LEAF FALL" MODE
C           DIRECTLY FROM "MAX. GROWTH" MODE.
C
C           NDL DCD
            IF(J.EQ.2)THEN
              IF(AILCG(I,J).LT.LFTHRS(I,J).AND.
     &        COLDDAYS(I,1).GE.COLDLMT(1).AND.
     &        AILCG(I,J).GT.ZERO)THEN
                LFSTATUS(I,J)=3        ! GO INTO "LEAF FALL" MODE
                CHKMODE(I,J)=1         
              ENDIF
            ENDIF
C
C           BDL DCD COLD 
            IF(J.EQ.4)THEN
              IF( AILCG(I,J).GT.ZERO.AND. 
     &        ((DAYLNGTH(I).LT.11.0.AND.ROOTTEMP(I,J).LT.(11.15+273.16)) 
     &        .OR. ROOTTEMP(I,J).LT.(ROOTHRSH+273.16)) )THEN
                LFSTATUS(I,J)=3        ! GO INTO "LEAF FALL" MODE
                CHKMODE(I,J)=1         
                FLHRLOSS(I,J)=GLEAFMAS(I,J)*(1.0/FLHRSPAN(2))
              ENDIF
            ENDIF

C           BDL DCD DRY
            IF(J.EQ.5)THEN
              IF( AILCG(I,J).GT.ZERO.AND. 
     &        ((DAYLNGTH(I).LT.11.0.AND.ROOTTEMP(I,J).LT.(11.15+273.16)) 
     &        .OR. ROOTTEMP(I,J).LT.(ROOTHRSH+273.16)) )THEN
                LFSTATUS(I,J)=3        ! GO INTO "LEAF FALL" MODE
                CHKMODE(I,J)=1         
              ENDIF
            ENDIF

          ENDIF
310     CONTINUE
300   CONTINUE
C
C     IF IN "NORMAL GROWTH" MODE
C     --------------------------
C
C     IF IN "NORMAL GROWTH" MODE THEN GO THROUGH EVERY PFT INDIVIDUALLY
C     AND FOLLOW SET OF RULES TO DETERMINE IF WE GO INTO "FALL/HARVEST"
C     MODE
C
      DO 320 I =  IL1, IL2
C
C       NEEDLE LEAF EVG
        IF(CHKMODE(I,1).EQ.0.AND.LFSTATUS(I,1).EQ.2)THEN
          IF(AILCG(I,1).LT.LFTHRS(I,1).AND.AILCG(I,1).GT.ZERO)THEN  
            LFSTATUS(I,1)=1         ! GO BACK TO "MAX. GROWTH" MODE
            CHKMODE(I,1)=1         
          ELSE IF(AILCG(I,1).LE.ZERO) THEN
            LFSTATUS(I,1)=4         ! SWITCH TO "NO LEAVES" MODE
            CHKMODE(I,1)=1         
            PANDAYS(I,1)=0
          ELSE
            LFSTATUS(I,1)=2         ! STAY IN "NORMAL GROWTH" MODE
            CHKMODE(I,1)=1         
          ENDIF
        ENDIF 
C
C       NEEDLE LEAF DCD
        IF(CHKMODE(I,2).EQ.0.AND.LFSTATUS(I,2).EQ.2)THEN
          IF(AILCG(I,2).LT.LFTHRS(I,2).AND.
     &      COLDDAYS(I,1).GE.COLDLMT(1).AND.
     &      AILCG(I,2).GT.ZERO)THEN
            LFSTATUS(I,2)=3         ! GO INTO "LEAF FALL" MODE
            CHKMODE(I,2)=1         
          ELSE IF(AILCG(I,2).LE.ZERO) THEN
            LFSTATUS(I,2)=4         ! SWITCH TO "NO LEAVES" MODE
            CHKMODE(I,2)=1         
            PANDAYS(I,2)=0
          ELSE
            LFSTATUS(I,2)=2         ! STAY IN "NORMAL GROWTH" MODE
            CHKMODE(I,2)=1         
          ENDIF
        ENDIF
C
C       BROAD LEAF EVG
        IF(CHKMODE(I,3).EQ.0.AND.LFSTATUS(I,3).EQ.2)THEN
          IF(AILCG(I,3).LT.LFTHRS(I,3).AND.AILCG(I,3).GT.ZERO)THEN  
            LFSTATUS(I,3)=1         ! GO BACK TO "MAX. GROWTH" MODE
            CHKMODE(I,3)=1         
          ELSE IF(AILCG(I,3).LE.ZERO) THEN
            LFSTATUS(I,3)=4         ! SWITCH TO "NO LEAVES" MODE
            CHKMODE(I,3)=1         
            PANDAYS(I,3)=0
          ELSE
            LFSTATUS(I,3)=2         ! STAY IN "NORMAL GROWTH" MODE
            CHKMODE(I,3)=1         
          ENDIF
        ENDIF 
C
C       BROAD LEAF DCD COLD
C       WE USE DAYLENGTH AND ROOTTEMP TO INITIATE LEAF OFFSET
        IF(CHKMODE(I,4).EQ.0.AND.LFSTATUS(I,4).EQ.2)THEN
          IF( AILCG(I,4).GT.ZERO.AND. 
     &    ((DAYLNGTH(I).LT.11.0.AND.ROOTTEMP(I,4).LT.(11.15+273.16))  
     &    .OR. ROOTTEMP(I,4).LT.(ROOTHRSH+273.16)) )THEN
            LFSTATUS(I,4)=3         ! GO INTO "LEAF FALL" MODE
            CHKMODE(I,4)=1         
            FLHRLOSS(I,4)=GLEAFMAS(I,4)*(1.0/FLHRSPAN(2))
          ELSE IF(AILCG(I,4).GT.ZERO.AND.AILCG(I,4).LT.LFTHRS(I,4))THEN   
            LFSTATUS(I,4)=1         ! SWITCH TO "MAX. GROWTH" MODE
            CHKMODE(I,4)=1         
          ELSE IF(AILCG(I,4).LE.ZERO) THEN
            LFSTATUS(I,4)=4         ! SWITCH TO "NO LEAVES" MODE
            CHKMODE(I,4)=1         
            PANDAYS(I,4)=0
            FLHRLOSS(I,4)=0.0
          ELSE
            LFSTATUS(I,4)=2         ! STAY IN "NORMAL GROWTH" MODE
            CHKMODE(I,4)=1         
          ENDIF
        ENDIF 
C
C       BROAD LEAF DCD DRY
C       WE STILL USE DAYLENGTH AND ROOTTEMP TO INITIATE LEAF OFFSET,
C       FOR THE PATHOLOGICAL CASES OF DRY DCD TREES BEING FURTHER
C       AWAY FROM THE EQUATOR THEN WE CAN IMAGINE. OTHER WISE LEAF
C       LOSS WILL OCCUR DUE TO DROUGHT ANYWAY.
        IF(CHKMODE(I,5).EQ.0.AND.LFSTATUS(I,5).EQ.2)THEN
          IF( AILCG(I,5).GT.ZERO.AND. 
     &    ((DAYLNGTH(I).LT.11.0.AND.ROOTTEMP(I,5).LT.(11.15+273.16))  
     &    .OR. ROOTTEMP(I,5).LT.(ROOTHRSH+273.16)) )THEN
            LFSTATUS(I,5)=3         ! GO INTO "LEAF FALL" MODE
            CHKMODE(I,5)=1         
          ELSE IF(AILCG(I,5).GT.ZERO.AND.AILCG(I,5).LT.LFTHRS(I,5))THEN   
            LFSTATUS(I,5)=1         ! SWITCH TO "MAX. GROWTH" MODE
            CHKMODE(I,5)=1         
          ELSE IF(AILCG(I,5).LE.ZERO) THEN
            LFSTATUS(I,5)=4         ! SWITCH TO "NO LEAVES" MODE
            CHKMODE(I,5)=1         
            PANDAYS(I,5)=0
          ELSE
            LFSTATUS(I,5)=2         ! STAY IN "NORMAL GROWTH" MODE
            CHKMODE(I,5)=1         
          ENDIF
        ENDIF 
C
320   CONTINUE
C
C     "NORMAL GROWTH" TO "FALL/HARVEST" TRANSITION FOR CROPS IS BASED ON
C     SPECIFIED LAI. WE HARVEST IF LAI OF CROPS REACHES A THRESHOLD.
C     IF LAI DOESN'T REACH THIS THRESHOLD (SAY DUE TO A BAD YEAR)
C     WE HARVEST ANYWAY IF IT STARTS GETTING COLD, OTHERWISE WE
C     DON'T HARVEST.
C
      DO 340 J = 6,7
       N = SORT(J)
        DO 350 I = IL1, IL2
          IF(CHKMODE(I,J).EQ.0.AND.LFSTATUS(I,J).EQ.2)THEN
            IF(AILCG(I,J).GE.HARVTHRS(N))THEN
              LFSTATUS(I,J)=3        ! GO INTO "HARVEST" MODE
              CHKMODE(I,J)=1
              FLHRLOSS(I,J)=GLEAFMAS(I,J)*(1.0/FLHRSPAN(1))
            ELSE IF( AILCG(I,J).GT.ZERO.AND.
     &      COLDDAYS(I,2).GE.COLDLMT(2) ) THEN
              LFSTATUS(I,J)=3        ! GO INTO "HARVEST" MODE
              CHKMODE(I,J)=1         ! REGARDLESS OF LAI
              FLHRLOSS(I,J)=GLEAFMAS(I,J)*(1.0/FLHRSPAN(1))
            ELSE IF(AILCG(I,J).LE.ZERO) THEN
              LFSTATUS(I,J)=4        ! SWITCH TO "NO LEAVES" MODE
              CHKMODE(I,J)=1         
              PANDAYS(I,J)=0
              FLHRLOSS(I,J)=0.0
            ELSE
              LFSTATUS(I,J)=2        ! STAY IN "NORMAL GROWTH" MODE
              CHKMODE(I,J)=1         
            ENDIF
          ENDIF
350     CONTINUE
340   CONTINUE
C
C     "NORMAL GROWTH" TO "MAX. GROWTH" TRANSITION FOR GRASSES 
C
      DO 370 J = 8,9
        DO 380 I = IL1, IL2
          IF(CHKMODE(I,J).EQ.0.AND.LFSTATUS(I,J).EQ.2)THEN
            IF(AILCG(I,J).LT.LFTHRS(I,J).AND.AILCG(I,J).GT.ZERO)THEN  
              LFSTATUS(I,J)=1        ! SWITCH BACK TO "MAX. GROWTH" MODE
              CHKMODE(I,J)=1
            ELSE IF(AILCG(I,J).LE.ZERO) THEN
              LFSTATUS(I,J)=4        ! SWITCH TO "NO LEAVES" MODE
              CHKMODE(I,J)=1         
              PANDAYS(I,J)=0
            ELSE
              LFSTATUS(I,J)=2        ! STAY IN "NORMAL GROWTH" MODE
              CHKMODE(I,J)=1
            ENDIF
          ENDIF
380     CONTINUE
370   CONTINUE 
C
C     IF IN "FALL/HARVEST" MODE
C     --------------------------
C
C     GRASSES AND EVG TREES DO NOT COME INTO THIS MODE, BECAUSE THEY WANT
C     TO STAY GREEN IF POSSIBLE. THIS MODE IS ACTIVATED FOR DCD PLANTS AND
C     CROPS. ONCE IN THIS MODE DCD TREES LOOSE THEIR LEAVES AND CROPS ARE
C     HARVESTED. NDL DCD TREES KEEP LOOSING THEIR LEAVES AT RATE DETERMINED
C     BY COLD STRESS, BDL DCD TREES LOOSE THEIR LEAVES AT A SPECIFIED
C     RATE, AND CROPS ARE HARVESTED OVER A PERIOD OF ~15 DAYS. DCD TREES
C     AND CROPS STAY IN "LEAF FALL/HARVEST" MODEL UNTIL ALL GREEN LEAVES
C     ARE GONE AT WHICH TIME THEY SWITCH INTO "NO LEAVES" MODE, AND THEN 
C     WAIT FOR THE CLIMATE TO BECOME FAVOURABLE TO GO INTO "MAX. GROWTH"
C     MODE
C

      DO 400 J = 1, ICC
        IF(J.EQ.2.OR.J.EQ.4.OR.J.EQ.5.OR.J.EQ.6.OR.J.EQ.7)THEN  !ONLY DCD TREES AND CROPS
          DO 410 I = IL1, IL2
            IF(CHKMODE(I,J).EQ.0.AND.LFSTATUS(I,J).EQ.3)THEN
              IF(AILCG(I,J).LE.0.01)THEN
                LFSTATUS(I,J)=4            ! GO INTO "NO LEAVES" MODE
                CHKMODE(I,J)=1
                PANDAYS(I,J)=0
                FLHRLOSS(I,J)=0.0
              ELSE
                IF(J.EQ.2)THEN             ! NDL DCD TREES
                  IF(PANDAYS(I,J).GE.DAYSCHK(J).AND.
     &            TA(I).GT.(COLDTHRS(1)+273.16)) THEN
                    IF(AILCG(I,J).LT.LFTHRS(I,J))THEN
                      LFSTATUS(I,J)=1      ! GO INTO "MAX. GROWTH" MODE
                      CHKMODE(I,J)=1
                    ELSE
                      LFSTATUS(I,J)=2      ! GO INTO "NORMAL GROWTH" MODE
                      CHKMODE(I,J)=1
                    ENDIF
                  ELSE  
                    LFSTATUS(I,J)=3        ! STAY IN "FALL/HARVEST" MODE 
                    CHKMODE(I,J)=1
                  ENDIF
                ELSE IF(J.EQ.4.OR.J.EQ.5)THEN        ! BDL DCD TREES
                  IF( (PANDAYS(I,J).GE.DAYSCHK(J)).AND.
     &            ((ROOTTEMP(I,4).GT.(ROOTHRSH+273.16)).AND.
     &             (DAYLNGTH(I).GT.11.0) ) )THEN 
                    IF(AILCG(I,J).LT.LFTHRS(I,J))THEN
                      LFSTATUS(I,J)=1      ! GO INTO "MAX. GROWTH" MODE
                      CHKMODE(I,J)=1
                    ELSE
                      LFSTATUS(I,J)=2      ! GO INTO "NORMAL GROWTH" MODE
                      CHKMODE(I,J)=1
                    ENDIF
                  ELSE  
                    LFSTATUS(I,J)=3        ! STAY IN "FALL/HARVEST" MODE 
                    CHKMODE(I,J)=1
                  ENDIF
                ELSE                       ! CROPS
                  LFSTATUS(I,J)=3          ! STAY IN "FALL/HARVEST" MODE 
                  CHKMODE(I,J)=1
                ENDIF
              ENDIF
            ENDIF
410       CONTINUE
        ENDIF
400   CONTINUE
C
C     CHECK THAT LEAF STATUS OF ALL VEGETATION TYPES IN ALL GRID CELLS HAS
C     BEEN UPDATED
C
      DO 411 J = 1, ICC
        DO 412 I = IL1, IL2
          IF(CHKMODE(I,J).EQ.0)THEN
           WRITE(6,2000) I,J
2000       FORMAT(' AT (I) = (',I3,'), PFT=',I2,' LFSTATUS NOT UPDATED')   
           CALL XIT('PHENOLGY',-2)
          ENDIF
412     CONTINUE
411   CONTINUE
C
C     ------------------------------------------------------------------     
C
C     HAVING DECIDED LEAF STATUS FOR EVERY PFT, WE NOW CALCULATE NORMAL
C     LEAF TURNOVER, COLD AND DROUGHT STRESS MORTALITY, AND FOR BDL DCD 
C     PLANTS WE ALSO CALCULATE SPECIFIED LOSS RATE IF THEY ARE IN "LEAF FALL"
C     MODE, AND FOR CROPS WE CALCULATE HARVEST LOSS, IF THEY ARE IN
C     "HARVEST" MODE.
C
C     ALL THESE LOSS CALCULATIONS WILL YIELD LEAF LITTER IN Kg C/M2 FOR
C     THE GIVEN DAY FOR ALL PFTs 
C
C     NORMAL LEAF TURN OVER
C
      DO 420 J = 1, ICC
        N = SORT(J)
        DO 430 I = IL1, IL2
          NRMLLOSS(I,J)=GLEAFMAS(I,J)*(1.0-EXP(-1.0/(365*LFESPANY(N))))     
430     CONTINUE
420   CONTINUE
C
C     FOR DROUGHT STRESS RELATED MORTALITY WE NEED FIELD CAPACITY 
C     AND WILTING POINT SOIL MOISTURE CONTENTS, WHICH WE CALCULATED
C     IN ALLOCATE SUBROUTINE
C
C
      DO 450 J = 1, IG
        DO 460 I = IL1, IL2
C
C         ESTIMATE (1-DROUGHT STRESS) 
C
          IF(THLIQ(I,J).LE.WILTSM(I,J)) THEN
            BETADRGT(I,J)=0.0
          ELSE IF(THLIQ(I,J).GT.WILTSM(I,J).AND.
     &      THLIQ(I,J).LT.FIELDSM(I,J))THEN
            BETADRGT(I,J)=(THLIQ(I,J)-WILTSM(I,J))
            BETADRGT(I,J)=BETADRGT(I,J)/(FIELDSM(I,J)-WILTSM(I,J))
          ELSE 
            BETADRGT(I,J)=1.0
          ENDIF          
          BETADRGT(I,J)=MAX(0.0, MIN(1.0,BETADRGT(I,J)))
C
460     CONTINUE
450   CONTINUE
C
C     ESTIMATE DROUGHT STRESS TERM AVERAGED OVER THE ROOTING DEPTH
C     FOR EACH PFT
C
      DO 480 J = 1, ICC
        DO 490 I = IL1, IL2
          DRGTSTRS(I,J) =  (1.0-BETADRGT(I,1))*RMATCTEM(I,J,1) +  
     &                     (1.0-BETADRGT(I,2))*RMATCTEM(I,J,2) +  
     &                     (1.0-BETADRGT(I,3))*RMATCTEM(I,J,3)   
          DRGTSTRS(I,J) = DRGTSTRS(I,J) /
     &     (RMATCTEM(I,J,1)+RMATCTEM(I,J,2)+RMATCTEM(I,J,3))  
          DRGTSTRS(I,J)=MAX(0.0, MIN(1.0,DRGTSTRS(I,J)))
490     CONTINUE
480   CONTINUE
C
C     USING THIS DROUGHT STRESS TERM AND OUR TWO VEGETATION-DEPENDENT
C     PARAMETERS WE FIND LEAF LOSS RATE ASSOCIATED WITH DROUGHT
C
      DO 500 J = 1, ICC
        N = SORT(J)
        DO 510 I = IL1, IL2
C         DROUGHT RELATED LEAF LOSS RATE
          DRGTLSRT(I,J)=DRLSRTMX(N)*(DRGTSTRS(I,J)**DRGTA(N))
510     CONTINUE
500   CONTINUE
C
C     ESTIMATE LEAF LOSS IN Kg C/M2 DUE TO DROUGHT STRESS
C
      DO 520 J = 1, ICC
        DO 530 I = IL1, IL2
          DRGTLOSS(I,J)=GLEAFMAS(I,J)*( 1.0-EXP(-DRGTLSRT(I,J)) )
530     CONTINUE
520   CONTINUE
C
C     SIMILAR TO DRGTSTRS WE FIND COLDSTRS FOR EACH PFT. WE ASSUME THAT
C     MAX. COLD STRESS RELATED LEAF LOSS OCCURS WHEN TEMPERATURE IS 5 C
C     OR MORE BELOW PFT's THRESHOLD
C
      DO 550 J = 1, ICC
        N = SORT(J)
        DO 560 I = IL1, IL2
          IF(TA(I).LE.(LWRTHRSH(N)-5.0+273.16))THEN
            COLDSTRS(I,J)=1.0
          ELSE IF(TA(I).GT.(LWRTHRSH(N)-5.0+273.16).AND.
     &    TA(I).LT.(LWRTHRSH(N)+273.16))THEN
            COLDSTRS(I,J)=1.0-((TA(I)-(LWRTHRSH(N)-5.0+273.16))/(5.0))   
          ELSE 
            COLDSTRS(I,J)=0.0
          ENDIF
          COLDSTRS(I,J)=MAX(0.0, MIN(1.0,COLDSTRS(I,J)))
560     CONTINUE
550   CONTINUE
C
C     USING THIS COLD STRESS TERM AND OUR TWO VEGETATION-DEPENDENT
C     PARAMETERS WE FIND LEAF LOSS RATE ASSOCIATED WITH COLD
C
      DO 570 J = 1, ICC
        N = SORT(J)
        DO 580 I = IL1, IL2
C         COLD RELATED LEAF LOSS RATE
          COLDLSRT(I,J)=CDLSRTMX(N)*(COLDSTRS(I,J)**COLDA(N))
580     CONTINUE
570   CONTINUE
c      write(*,*) TA(1),LWRTHRSH(1)
c      write(*,1200) (COLDSTRS(1,J),J=1,ICC)
c      write(*,1200) (COLDLSRT(1,J),J=1,ICC)
1200  FORMAT(9F10.2)
C
C     ESTIMATE LEAF LOSS IN Kg C/M2 DUE TO COLD STRESS
C
      DO 600 J = 1, ICC
        DO 610 I = IL1, IL2
          COLDLOSS(I,J)=GLEAFMAS(I,J)*( 1.0-EXP(-COLDLSRT(I,J)) )
610     CONTINUE
600   CONTINUE
C
C     NOW THAT WE HAVE ALL TYPES OF LEAF LOSSES (DUE TO NORMAL TURNOVER,
C     COLD AND DROUGHT STRESS, AND FALL/HARVEST) WE TAKE THE LOSSES
C     FOR GRASSES AND USE THOSE TO TURN LIVE GREEN GRASS INTO DEAD
C     BROWN GRASS. WE THEN FIND THE LEAF LITTER FROM THE BROWN GRASS
C     WHICH WILL THEN GO INTO THE LITTER POOL.
C
      DO 620 J = 8,9
        N = SORT(J)
        DO 630 I = IL1, IL2
          GLEAFMAS(I,J)   = GLEAFMAS(I,J)-NRMLLOSS(I,J)-DRGTLOSS(I,J)-  
     &                      COLDLOSS(I,J)
          IF( GLEAFMAS(I,J).LT.0.0) THEN
            BLEAFMAS(I,J) = BLEAFMAS(I,J)+NRMLLOSS(I,J)+DRGTLOSS(I,J)+  
     &                      COLDLOSS(I,J)+GLEAFMAS(I,J)
            GLEAFMAS(I,J)=0.0
          ELSE
            BLEAFMAS(I,J) = BLEAFMAS(I,J)+NRMLLOSS(I,J)+DRGTLOSS(I,J)+  
     &                      COLDLOSS(I,J)
          ENDIF
          NRMLLOSS(I,J) = 0.0
          DRGTLOSS(I,J) = 0.0
          COLDLOSS(I,J) = 0.0
C         WE ASSUME LIFE SPAN OF BROWN GRASS IS 10% THAT OF GREEN GRASS
C         BUT THIS IS AN ADJUSTABLE PARAMETER.
          NRMLLOSS(I,J) = BLEAFMAS(I,J)*
     &      (1.0-EXP(-1.0/(0.10*365*LFESPANY(N))))       
630     CONTINUE
620   CONTINUE 
C
C     COMBINE NRMLLOSS, DRGTLOSS, AND COLDLOSS TOGETHER TO GET TOTAL
C     LEAF LITTER GENERATED, WHICH IS WHAT WE USE TO UPDATE GLEAFMASS,
C     EXCEPT FOR GRASSES, FOR WHICH WE HAVE ALERADY UPDATED GLEAFMASS.
C
      DO 650 J = 1, ICC
        DO 660 I = IL1, IL2
          LEAFLITR(I,J) = NRMLLOSS(I,J)+DRGTLOSS(I,J)+COLDLOSS(I,J)+
     &                    FLHRLOSS(I,J)
660     CONTINUE
650   CONTINUE
c      write(*,'(I4,4F14.6)') IDAY,NRMLLOSS(1,1),DRGTLOSS(1,1),
c     1                       COLDLOSS(1,1),FLHRLOSS(1,1)
C
      RETURN
      END

