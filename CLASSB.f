      SUBROUTINE CLASSB(THPOR,THLRET,THLMIN,BI,PSISAT,GRKSAT,
     1                  THLRAT,HCPS,TCS,THFC,THLW,PSIWLT,
     2                  DELZW,ZBOTW,ALGWET,ALGDRY,
     +                  ALGWV,ALGWN,ALGDV,ALGDN,            
     3                  SAND,CLAY,ORGM,SOCI,DELZ,ZBOT,SDEPTH,           
     4                  ISAND,IGDR,NL,NM,IL1,IL2,IM,IG,ICTEMMOD)                 
C
C     * JAN 15/15 - D.VERSEGHY. CHANGE PSIWLT FOR MINERAL SOILS
C     *                         TO A CONSTANT VALUE OF 150 M.
C     *                         AND ADD NEW VARIABLE THLW.
C     * AUG 25/14 - M.LAZARE.   PASS IN NEW WET AND DRY SOIL
C     *                         BRIGHTNESS FIELDS FROM CLM.
C     * NOV 16/13 - M.LAZARE.   FINAL VERSION FOR GCM17:                
C     *                         - REVERT BACK TO CLASS2.7               
C     *                           SPECIFICATION FOR "ALGWET".           
C     * NOV 11/11 - M.LAZARE.   - IMPLEMENT CTEM CHOICE OF
C     *                           ALGDRY DETERMINED BY ADDED
C     *                           PASSED SWITCH "ICTEMMOD". 
C     * OCT 18/11 - M.LAZARE.   - REMOVE UNUSED "IORG".
C     *                         - CHANGE "THSAND", "THORG"
C     *                           AND "THFINE" FROM ARRAYS
C     *                           (INTERNAL ONLY) TO SCALAR. 
C     *                         - IGDR NOW PASSED OUT TO BE
C     *                           USED IN GRINFL/GRDRAN/WEND.
C     *                         - PASS IN IL1 AND IL2 TO
C     *                           DEFINE LOOPS.
C     * OCT 08/11 - M.LAZARE.   ALGDRY CHANGED BACK TO FORMULA
C     *                         USED IN GCM15I (.0056->.0046).
C     * SEP 27/11 - D.VERSEGHY. CONSTRAIN DELZW TO BE >= 5 CMS
C     *                         TO AVOID PROBLEMATIC UNDERSHOOTS.
C     * AUG 25/11 - D.VERSEGHY. USE THFC FORMULATION FOR BOTTOM
C     *                         LAYER AT BOTTOM OF SOIL PERMEABLE
C     *                         DEPTH.
C     * DEC 23/09 - V.FORTIN.   REVISE CALCULATION OF THFC FOR
C     *                         BOTTOM LAYER IN MINERAL SOILS 
C     *                         ACCORDING TO SOULIS ET AL. (2009).
C     * JAN 06/09 - D.VERSEGHY. REVERSE ORDER OF 200 AND 300 LOOPS.
C     * DEC 11/07 - D.VERSEGHY. CHANGE CALCULATION OF TCS FROM
C     *                         GEOMETRIC MEAN TO LINEAR MEAN.
C     * FEB 07/07 - D.VERSEGHY. SET THFC TO THLRET FOR ORGANIC SOILS;
C     *                         STREAMLINE SOME CALCULATIONS.
C     * SEP 15/05 - D.VERSEGHY. REMOVE HARD CODING OF IG=3 IN 300 LOOP.
C     * APR 06/05 - D.VERSEGHY. MOVE CALCULATION OF GRKTLD
C     *                         INTO GRINFL; REVISED CALCULATION
C     *                         OF ALGDRY (WITH M.LAZARE).
C     * NOV 03/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
C     * SEP 04/03 - D.VERSEGHY. PERMEABLE THICKNESS OF SOIL
C     *                         LAYERS CONSTRAINED TO >= 1 MM;
C     *                         PROTECT SENSITIVE CALCULATIONS
C     *                         AGAINST ROUNDOFF ERROR.
C     * JUN 28/02 - D.VERSEGHY. ASSIGN SOIL HYDROLOGICAL AND
C     *                         THERMAL PROPERTIES BASED ON
C     *                         SAND, CLAY AND ORGANIC MATTER
C     *                         CONTENT.
C
      IMPLICIT NONE
C
C     * INTEGER CONSTANTS.
C
      INTEGER NL,NM,IL1,IL2,IM,IG,ICTEMMOD,I,J,M
C
C     * OUTPUT ARRAYS.
C
      REAL THPOR (NL,NM,IG),  THLRET(NL,NM,IG),  THLMIN(NL,NM,IG),
     1     BI    (NL,NM,IG),  PSISAT(NL,NM,IG),  GRKSAT(NL,NM,IG),  
     2     THLRAT(NL,NM,IG),  HCPS  (NL,NM,IG),  
     3     TCS   (NL,NM,IG),  THFC  (NL,NM,IG),  PSIWLT(NL,NM,IG),
     4     DELZW (NL,NM,IG),  ZBOTW (NL,NM,IG),
     4     ALGWET(NL,NM),     ALGDRY(NL,NM), THLW(NL,NM,IG)
C
      INTEGER                 ISAND (NL,NM,IG),  IGDR  (NL,NM)
C
C     * INPUT ARRAYS.
C
      REAL SAND  (NL,NM,IG),  CLAY  (NL,NM,IG),  ORGM  (NL,NM,IG),
     1     DELZ  (IG),        ZBOT  (IG),        SDEPTH(NL,NM) 
C
      REAL SOCI  (NL,NM)                                                
C                                                                       
      REAL THPORG (3),      THRORG (3),      THMORG (3),
     1     BORG   (3),      PSISORG(3),      GRKSORG(3)
C
C     * TEMPORARY VARIABLES.
C
      REAL ALWV(20), ALWN(20), ALDV(20), ALDN(20)                       
C                                                                       
      REAL VSAND,VORG,VFINE,VTOT,AEXP,ABC,THSAND,THFINE,THORG
C
C     * COMMON BLOCK PARAMETERS.
C
      REAL TCW,TCICE,TCSAND,TCFINE,TCOM,TCDRYS,RHOSOL,RHOOM,
     1     HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPFIN,SPHW,SPHICE,SPHVEG,
     2     SPHAIR,RHOW,RHOICE,TCGLAC,CLHMLT,CLHVAP 
C
      COMMON /CLASS3/ TCW,TCICE,TCSAND,TCFINE,TCOM,TCDRYS,
     1                RHOSOL,RHOOM
      COMMON /CLASS4/ HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPFIN,
     1                SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,
     2                TCGLAC,CLHMLT,CLHVAP
      COMMON /CLASS5/ THPORG,THRORG,THMORG,BORG,PSISORG,GRKSORG
C                                                                       
      DATA ALWV /0.25,0.23,0.21,0.20,0.19,0.18,0.17,0.16,0.15,0.14,0.13,
     1           0.12,0.11,0.10,0.09,0.08,0.07,0.06,0.05,0.04/          
      DATA ALWN /0.50,0.46,0.42,0.40,0.38,0.36,0.34,0.32,0.30,0.28,0.26,
     1           0.24,0.22,0.20,0.18,0.16,0.14,0.12,0.10,0.08/          
      DATA ALDV /0.36,0.34,0.32,0.31,0.30,0.29,0.28,0.27,0.26,0.25,0.24,
     1           0.23,0.22,0.20,0.18,0.16,0.14,0.12,0.10,0.08/          
      DATA ALDN /0.61,0.57,0.53,0.51,0.49,0.48,0.45,0.43,0.41,0.39,0.37,
     1           0.35,0.33,0.31,0.29,0.27,0.25,0.23,0.21,0.16/          
C---------------------------------------------------------------------
C
      DO 50 M=1,IM
      DO 50 I=IL1,IL2
          IGDR(I,M)=1
50    CONTINUE
C
      DO 100 J=1,IG
      DO 100 M=1,IM
      DO 100 I=IL1,IL2
          ISAND (I,M,J)=NINT(SAND(I,M,J))                                               
          IF(ISAND(I,M,J).GT.-3) IGDR(I,M)=J
100   CONTINUE
C
      DO 200 M=1,IM
      DO 200 I=IL1,IL2
          DO 150 J=1,IG
              IF(ISAND(I,M,1).EQ.-4) THEN
                  DELZW(I,M,J)=DELZ(J)
                  ISAND(I,M,J)=-4
              ELSEIF(ISAND(I,M,J).EQ.-3) THEN
                  DELZW(I,M,J)=0.0
              ELSEIF(SDEPTH(I,M).GE.ZBOT(J)) THEN
                  DELZW(I,M,J)=DELZ(J)
              ELSEIF(SDEPTH(I,M).LT.(ZBOT(J)-DELZ(J)+0.025)) THEN
                  DELZW(I,M,J)=0.0
                  ISAND(I,M,J)=-3
              ELSE
                  DELZW(I,M,J)=MAX(0.05,(SDEPTH(I,M)-(ZBOT(J)-DELZ(J))))
              ENDIF
              ZBOTW(I,M,J)=MAX(0.0,ZBOT(J)-DELZ(J))+DELZW(I,M,J)
150       CONTINUE
          IF(SAND(I,M,1).GE.0.0) THEN
              ALGWET(I,M)=0.08+0.0022*SAND(I,M,1)
              IF(ICTEMMOD.EQ.0) THEN
                  ALGDRY(I,M)=MIN(0.14+0.0046*SAND(I,M,1), 0.45) ! FOR GLC2000
              ELSE
                  ALGDRY(I,M)=MIN(0.14+0.0027*SAND(I,M,1), 0.41) ! FOR CTEM   
              ENDIF
              ALGWV(I,M)=ALWV(NINT(SOCI(I,M)))                          
              ALGWN(I,M)=ALWN(NINT(SOCI(I,M)))                          
              ALGDV(I,M)=ALDV(NINT(SOCI(I,M)))                          
              ALGDN(I,M)=ALDN(NINT(SOCI(I,M)))
          ELSE
              ALGWET(I,M)=0.0
              ALGDRY(I,M)=0.0
              ALGWV(I,M)=0.0                                            
              ALGWN(I,M)=0.0                                            
              ALGDV(I,M)=0.0                                            
              ALGDN(I,M)=0.0                                            
          ENDIF
200   CONTINUE
C
      DO 300 J=1,IG
      DO 300 M=1,IM
      DO 300 I=IL1,IL2
          IF(ISAND(I,M,J).EQ.-4) THEN
              THPOR (I,M,J)=0.0
              THLRET(I,M,J)=0.0
              THLMIN(I,M,J)=0.0
              BI    (I,M,J)=0.0
              PSISAT(I,M,J)=0.0
              GRKSAT(I,M,J)=0.0
              THLRAT(I,M,J)=0.0
              HCPS(I,M,J)=HCPICE
              TCS(I,M,J)=TCICE
              THFC(I,M,J)=0.0
              THLW(I,M,J)=0.0
              PSIWLT(I,M,J)=0.0                                         
          ELSEIF(ISAND(I,M,J).EQ.-3) THEN
              THPOR (I,M,J)=0.0
              THLRET(I,M,J)=0.0
              THLMIN(I,M,J)=0.0
              BI    (I,M,J)=0.0
              PSISAT(I,M,J)=0.0
              GRKSAT(I,M,J)=0.0
              THLRAT(I,M,J)=0.0
              HCPS(I,M,J)=HCPSND
              TCS(I,M,J)=TCSAND
              THFC(I,M,J)=0.0
              THLW(I,M,J)=0.0              
              PSIWLT(I,M,J)=0.0                                         
          ELSEIF(ISAND(I,M,J).EQ.-2) THEN
              THPOR (I,M,J)=THPORG(MIN(J,3))
              THLRET(I,M,J)=THRORG(MIN(J,3))
              THLMIN(I,M,J)=THMORG(MIN(J,3))
              BI    (I,M,J)=BORG(MIN(J,3))
              PSISAT(I,M,J)=PSISORG(MIN(J,3))
              GRKSAT(I,M,J)=GRKSORG(MIN(J,3))
              THLRAT(I,M,J)=0.5**(1.0/(2.0*BI(I,M,J)+3.0))
              HCPS(I,M,J)=HCPOM
              TCS(I,M,J)=TCOM
              THFC(I,M,J)=THLRET(I,M,J)           
              PSIWLT(I,M,J)=PSISAT(I,M,J)*(THLMIN(I,M,J)/
     1            THPOR(I,M,J))**(-BI(I,M,J))
              THLW(I,M,J)=THLMIN(I,M,J)
          ELSEIF(SAND(I,M,J).GT.0.) THEN
              THPOR (I,M,J)=(-0.126*SAND(I,M,J)+48.9)/100.0
              THLRET(I,M,J)=0.04
              THLMIN(I,M,J)=0.04
              BI    (I,M,J)=0.159*CLAY(I,M,J)+2.91
              PSISAT(I,M,J)=0.01*EXP(-0.0302*SAND(I,M,J)+4.33)
              GRKSAT(I,M,J)=7.0556E-6*(EXP(0.0352*SAND(I,M,J)-2.035))
              THLRAT(I,M,J)=0.5**(1.0/(2.0*BI(I,M,J)+3.0))
              VSAND=SAND(I,M,J)/(RHOSOL*100.0)
              VORG=ORGM(I,M,J)/(RHOOM*100.0)
              VFINE=(100.0-SAND(I,M,J)-ORGM(I,M,J))/(RHOSOL*100.0)
              VTOT=VSAND+VFINE+VORG
              THSAND=(1.0-THPOR(I,M,J))*VSAND/VTOT
              THORG=(1.0-THPOR(I,M,J))*VORG/VTOT
              THFINE=1.0-THPOR(I,M,J)-THSAND-THORG
              HCPS(I,M,J)=(HCPSND*THSAND+HCPFIN*THFINE+
     1            HCPOM*THORG)/(1.0-THPOR(I,M,J))
              TCS(I,M,J)=(TCSAND*THSAND+TCOM*THORG+
     1            TCFINE*THFINE)/(1.0-THPOR(I,M,J))
              IF(J.NE.IGDR(I,M))                       THEN
                  THFC(I,M,J)=THPOR(I,M,J)*(1.157E-9/GRKSAT(I,M,J))**
     1                (1.0/(2.0*BI(I,M,J)+3.0))
              ELSE
                  AEXP=(BI(I,M,J)-1.0)/BI(I,M,J)
                  ABC=(3.0*BI(I,M,J)+2.0)**AEXP-
     1                (2.0*BI(I,M,J)+2.0)**AEXP
                  THFC(I,M,J)=(ABC*THPOR(I,M,J)/(BI(I,M,J)-1.0))*
     1                (PSISAT(I,M,J)*BI(I,M,J)/SDEPTH(I,M))**
     2                (1.0/BI(I,M,J))
              ENDIF

              PSIWLT(I,M,J)=150.0
              THLW(I,M,J)=THPOR(I,M,J)*(PSIWLT(I,M,J)/PSISAT(I,M,J))**
     1                    (-1.0/BI(I,M,J))
     
          ENDIF
300   CONTINUE
C
      RETURN
      END
