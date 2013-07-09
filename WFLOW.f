      SUBROUTINE WFLOW(WMOVE,TMOVE,LZF,NINF,TRMDR,TPOND,ZPOND,
     1                 R,TR,EVAP,PSIF,GRKINF,THLINF,THLIQX,TBARWX,
     2                 DELZX,ZBOTX,FMAX,ZF,DZF,DTFLOW,THLNLZ,
     3                 THLQLZ,DZDISP,WDISP,WABS,ITER,NEND,ISIMP,
     4                 IGRN,IG,IGP1,IGP2,ILG,IL1,IL2,JL,N )
     
C     * Evaluates infiltration of water into soil under saturated 
C     * conditions.
                                                                         
C     * FEB 09/09 - J.P.BLANCHETTE. INCREASE LIMITING VALUE OF NEND.
C     * JAN 06/09 - D.VERSEGHY. CHECKS ON WETTING FRONT LOCATION
C     *                         IN 500 LOOP.
C     * AUG 07/07 - D.VERSEGHY. INCREASE ITERATION COUNTER NEND;
C     *                         MOVE CALCULATION OF FMAX FROM
C     *                         GRINFL TO THIS ROUTINE.
C     * MAY 17/06 - D.VERSEGHY. PROTECT AGAINST DIVISIONS BY ZERO.
C     * SEP 13/05 - D.VERSEGHY. REPLACE HARD-CODED 4 WITH IGP1.
C     * SEP 23/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
C     * JUL 30/04 - D.VERSEGHY/Y.DELAGE. PROTECT SENSITIVE 
C     *                         CALCULATIONS AGAINST ROUNDOFF ERRORS.
C     * JUN 21/02 - D.VERSEGHY. UPDATE SUBROUTINE CALL.
C     * MAR 04/02 - D.VERSEGHY. DEFINE "NEND" FOR ALL CASES.
C     * DEC 16/94 - D.VERSEGHY. BUG FIX - SPECIFY TMOVE BEHIND
C     *                         WETTING FRONT AFTER ANY FULL-LAYER 
C     *                         JUMP DOWNWARD.
C     * APR 24/92 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.2.
C     *                                  REVISED AND VECTORIZED CODE 
C     *                                  FOR MODEL VERSION GCM7.
C     * AUG 12/91 - D.VERSEGHY. CODE FOR MODEL VERSION GCM7U -
C     *                         CLASS VERSION 2.0 (WITH CANOPY).
C     * APR 11/89 - D.VERSEGHY. SATURATED FLOW OF WATER THROUGH SOIL.
C
      IMPLICIT NONE
C
C     * INTEGER CONSTANTS.
C
      INTEGER IG,IGP1,IGP2,ILG,IL1,IL2,JL,I,J,NPNTS,NIT,N
C
C     * INPUT/OUTPUT FIELDS.
C                      
      REAL WMOVE (ILG,IGP2)  !Water movement matrix [m3 m-2]
      REAL TMOVE (ILG,IGP2)  !Temperature matrix associated with ground water movement [C]
C
      INTEGER  LZF   (ILG)  !Index of soil layer in which wetting front is located
      INTEGER  NINF  (ILG)  !Number of levels involved in water movement
C
      REAL TRMDR (ILG)  !Time remaining in current time step [s]
      REAL TPOND (ILG)  !Temperature of ponded water [C]
      REAL ZPOND (ILG)  !Depth of ponded water [m] (zp)
C
C     * INPUT FIELDS.
C
      REAL R     (ILG)       !Rainfall rate at ground surface [m s-1]
      REAL TR    (ILG)       !Temperature of rainfall [C]
      REAL EVAP  (ILG)       !Surface evaporation rate [m s-1]
      REAL PSIF  (ILG,IGP1)  !Soil water suction across the wetting front [m] (psi_f)
      REAL GRKINF(ILG,IGP1)  !Hydraulic conductivity of soil behind the wetting front [m s-1] (K)
      REAL THLINF(ILG,IGP1)  !Volumetric liquid water content behind the wetting front [m3 m-3]
      REAL THLIQX(ILG,IGP1)  !Volumetric liquid water content of soil layer [m3 m-3]
      REAL TBARWX(ILG,IGP1)  !Temperature of water in soil layer [C]
      REAL DELZX (ILG,IGP1)  !Permeable depth of soil layer [m] (delta_z)
      REAL ZBOTX (ILG,IGP1)  !Depth of bottom of soil layer [m]
      REAL FMAX  (ILG)       !Maximum infiltration rate, defined as minimum value of GRKINF
      REAL ZF    (ILG)       !Depth of the wetting front [m] (zf)
C
      INTEGER   IGRN  (ILG)  !Flag to indicate whether infiltration is occurring
C
C     * INTERNAL WORK FIELDS.
C
      REAL DZF   (ILG),        DTFLOW(ILG),        THLNLZ(ILG),        
     1     THLQLZ(ILG),        DZDISP(ILG),        WDISP (ILG),        
     2     WABS  (ILG)
C
      INTEGER                  ITER  (ILG),        NEND  (ILG),
     1                         ISIMP (ILG)    
C 
C     * TEMPORARY VARIABLES.
C
      REAL RESID,FINF,ZPTEST,WINF
C 
C     * COMMON BLOCK PARAMETERS.
C
      REAL DELT     !Time step [s]
      REAL TFREZ    !Freezing point of water [K]
C
      COMMON /CLASS1/ DELT,TFREZ
C-----------------------------------------------------------------------
C
C     General calculations performed:
C     *The infiltration rate Finf under saturated conditions is 
C     *calculated as
C     *
C     *Finf = GRKINF [(PSIF + ZF + ZPOND)/ZF ]
C     *
C     *where GRKINF is the hydraulic conductivity of the soil behind the
C     *wetting front, PSIF is the soil C*moisture suction across the 
C     *wetting front, ZF is the depth of the wetting front, and ZP is 
C     *the depth of water ponded on the surface. Since PSIF will vary 
C     *by soil layer, and zf and zp will vary with time, the period of 
C     *infiltration is divided into two-minute segments. A maximum
C     *iteration flag NEND is defined as the number of iteration 
C     *segments plus 100, as a safeguard against runaway iterations.
C     *Since the surface infiltration rate will be limited by the
C     *lowest infiltration rate encountered, a maximum infiltration rate
C     *FMAX is defined as the minimum value of GRKINF, the hydraulic
C     *conductivity behind the wetting front, in all of the soil layers
C     *above and including that containing the wetting front.
C
C     * CALCULATE ITERATION ENDPOINT "NEND", AND SET SWITCH "ITER" TO 1
C     * FOR POINTS OVER WHICH THIS SUBROUTINE IS TO BE PERFORMED.
C
      DO 50 I=IL1,IL2
          IF(IGRN(I).GT.0 .AND. TRMDR(I).GT.0.0)                    THEN
              RESID=MOD(TRMDR(I),120.)                                                      
              IF(RESID.GT.0.)                       THEN  
                  NEND(I)=NINT(TRMDR(I)/120.+0.5)+100
              ELSE                                                                        
                  NEND(I)=NINT(TRMDR(I)/120.)+100
              ENDIF
              ITER(I)=1
              FMAX(I)=999999. 
          ELSE
              NEND(I)=0  
              ITER(I)=0  
          ENDIF
   50 CONTINUE
      NIT=1
C                     
      DO 75 J=1,IGP1
      DO 75 I=IL1,IL2
          IF(ITER(I).GT.0 .AND. LZF(I).GT.1 .AND. J.LT.LZF(I))     THEN   
              FMAX(I)=MIN(GRKINF(I,J),FMAX(I))                                      
          ENDIF                                                                   
   75 CONTINUE
C
  100 CONTINUE

C     * BEGINNING OF ITERATION SEQUENCE.
C     * SET OR RESET NUMBER OF POINTS TO BE PROCESSED ON THE CURRENT 
C     * LATITUDE CIRCLE(S).
C
      NPNTS=0   
C      
C     * In loop 200, a check is performed to determine whether the 
C     * liquid water content of the current soil layer, THLIQX, equals 
C     * or exceeds that behind the wetting front, THLINF. If so, the 
C     * depth of the wetting front is reset to the depth of the soil 
C     * layer. The water content of the soil layer is assigned to the
C     * level of the water movement matrix WMOVE corresponding to the 
C     * soil layer index plus 1, and the temperature of the water in the
C     * layer is assigned to the same level of the matrix TMOVE. FMAX is
C     * recalculated, and the flags LZF and NINF are each incremented
C     * by 1. The flag ISIMP is set to 1, indicating that the following
C     * calculations are to be bypassed for this iteration.
C
C     * If ISIMP is not 1, the time period DTFLOW applying to the
C     * current iteration loop is set to two minutes or to the
C     * remainder of the time step, whichever is less. If GRKINF for
C     * the current soil layer is not vanishingly small, the flag ISIMP
C     * is set to -2. Otherwise, it is deemed that infiltration is
C     * suppressed, and the temperature and depth of water ponded on
C     * the surface are simply updated. The temperature of the ponded
C     * water is calculated as the weighted average of the current pond
C     * temperature and the rainfall added to it. The ponded water
C     * remaining after rainfall and evaporation have taken place is
C     * calculated as ZPTEST. If ZPTEST is less than zero, it is
C     * deduced that evaporation must have removed the ponded water
C     * before the end of the current iteration loop. If this is the
C     * case, the time period of the current iteration loop is
C     * recalculated as the amount of time required for the difference
C     * between the evaporation
C     * and the rainfall rates to consume the ponded water, and the
C     * pond depth ZPOND is set to zero. Otherwise, ZPOND is set to
C     * ZPTEST. Finally, ISIMP is set to -1.
                               
C
C     * IF THE WATER CONTENT OF THE CURRENT SOIL LAYER EQUALS OR EXCEEDS
C     * THE WATER CONTENT BEHIND THE WETTING FRONT, INSTANTANEOUSLY 
C     * RELOCATE WETTING FRONT TO BOTTOM OF CURRENT SOIL LAYER; 
C     * RE-EVALUATE INFILTRATION PARAMETERS, UPDATE WATER MOVEMENT 
C     * MATRIX, SET SWITCH "ISIMP" TO 1 AND DROP TO END OF ITERATION 
C     * LOOP.
C     * (SOME ARRAYS ARE GATHERED (ON LZF) TO AVOID MULTIPLE INDIRECT-
C     * ADDRESSING REFERENCES IN THE ENSUING LOOPS.)
C
      DO 200 I=IL1,IL2
          IF(ITER(I).EQ.1)                                          THEN
              THLNLZ(I)=THLINF(I,LZF(I))
              THLQLZ(I)=THLIQX(I,LZF(I))
              IF(THLQLZ(I).GT.(THLNLZ(I)-1.0E-6) .AND. LZF(I).LT.IGP1
     1                .AND. THLQLZ(I).GT.0.0001)               THEN             
                  ZF(I)=ZBOTX(I,LZF(I))                                                   
                  WMOVE(I,NINF(I))=THLQLZ(I)*DELZX(I,LZF(I))                              
                  TMOVE(I,NINF(I))=TBARWX(I,LZF(I))                              
                  FINF=GRKINF(I,LZF(I))*(ZF(I)+ZPOND(I)+PSIF(I,LZF(I)))/
     1                 ZF(I)                        
                  FMAX(I)=MIN(FMAX(I),FINF)                                      
                  LZF(I)=LZF(I)+1                                                       
                  NINF(I)=NINF(I)+1
                  ISIMP(I)=1                                                     
              ELSE                                            
                  ISIMP(I)=0
              ENDIF
          ELSE
              ISIMP(I)=0
          ENDIF
  200 CONTINUE
C
C     * INFILTRATION CALCULATIONS TAKING FINITE TIME. SET TIMESTEP OF
C     * CURRENT ITERATION PASS AND CHECK HYDRAULIC CONDUCTIVITY OF
C     * CURRENT SOIL LAYER. IF ZERO, RECALCULATE POND DEPTH AND POND
C     * TEMPERATURE AND SET "ISIMP" TO -1; ELSE SET "ISIMP" TO -2.
C
      DO 300 I=IL1,IL2
          IF(ITER(I).EQ.1 .AND. ISIMP(I).NE.1)                      THEN
              DTFLOW(I)=MIN(TRMDR(I),120.)
              IF(GRKINF(I,LZF(I)).GT.1.0E-12)                   THEN
                  ISIMP(I)=-2
              ELSE
                  TPOND(I)=(ZPOND(I)*TPOND(I)+R(I)*DTFLOW(I)*TR(I))/
     1                     (ZPOND(I)+R(I)*DTFLOW(I))
                  IF(ABS(R(I)-EVAP(I)).GT.1.0E-11)            THEN                                 
                      ZPTEST=ZPOND(I)+(R(I)-EVAP(I))*DTFLOW(I)                               
                  ELSE                                                        
                      ZPTEST=ZPOND(I)                                            
                  ENDIF                                                       
                  IF(ZPTEST.LT.0.)                        THEN                                      
                      DTFLOW(I)=ZPOND(I)/(EVAP(I)-R(I))                                       
                      ZPOND(I)=0.0                                               
                  ELSE                                                        
                      ZPOND(I)=ZPTEST                                            
                  ENDIF                                                       
                  ISIMP(I)=-1
              ENDIF
          ENDIF
  300 CONTINUE
C 
C     The 400 loop addresses the infiltration process under saturated 
C     conditions. Such infiltration is modelled as “piston flow”. First 
C     the current infiltration rate FINF is calculated using the 
C     equation given above. If the wetting front has passed the bottom 
C     of the lowest soil layer, PSIF Cis neglected. If FINF is greater 
C     than the rainfall rate R, FINF is set equal to R; if FINF is 
C     greater than CFMAX, FINF is set equal to FMAX. If the wetting 
C     front has not passed the bottom of the lowest soil layer, the 
C     change in depth of the wetting front DZF over the current time 
C     interval is calculated from the amount of infiltrating water WINF 
C     and the volumetric water content behind the wetting front, THLINF. 
C     The amount of soil water WDISP displaced by this movement of the 
C     wetting front is calculated from DZF and THLIQX, and the depth to 
C     which this water penetrates, DZDISP, is calculated as 
C     WDISP/(THLINF – THLIQX). The amount of soil water WABS entrained 
C     by the movement of WDISP itself is calculated as the product of 
C     DZDISP and THLIQX. The change in depth of the wetting front, 
C     behind which the liquid water content of the   soil layer is 
C     THLINF, is now the sum of the depth represented by infiltration of 
C     water at the surface, DZF, and the depth represented by 
C     displacement of the pre-existing soil water, DZDISP. If this 
C     change in depth causes the wetting front to move beyond the bottom 
C     of the current soil layer, DTFLOW is recalculated as the amount of 
C     time required for the composite wetting front to reach the bottom 
C     of the soil layer, and DZF, WDISP, DZDISP and WABS are likewise 
C     recalculated. As in the case for ISIMP = -1, the temperature of 
C     the ponded water on the surface is calculated as the weighted 
C     average of the current pond temperature and the rainfall added to 
C     it. The ponded water remaining after rainfall, infiltration and 
C     evaporation have taken place is calculated as ZPTEST. If ZPTEST is 
C     less than zero, it is deduced that infiltration and evaporation 
C     must have removed the ponded water before the end of the current 
C     iteration loop. If this is the case, the time period of the 
C     current iteration loop is recalculated as the amount of time 
C     required for the infiltration and evaporation minus the rainfall 
C     to consume the ponded water. The pond depth ZPOND is set to zero; 
C     if the wetting front has not passed the bottom of the lowest soil 
C     layer, DZF, WDISP, DZDISP and WABS are recalculated. Otherwise, 
C     ZPOND is set to ZPTEST. Finally, the first layer of the water 
C     movement matrix WMOVE is incremented with the amount of the 
C     infiltrated water, and the first layer of the matrix TMOVE with 
C     the temperature of the infiltrated water. The last layer of WMOVE 
C     is incremented by WDISP+WABS, and the depth of the wetting front 
C     by DZF+DZDISP.
C
C     * "ISIMP"=-2: NORMAL SATURATED INFILTRATION UNDER PISTON-FLOW
C     * CONDITIONS. CALCULATE CURRENT INFILTRATION RATE (FINF); WATER
C     * INFILTRATING DURING CURRENT ITERATION PASS (WINF); SOIL WATER
C     * DISPLACED INTO EMPTY PORES AHEAD OF WETTING FRONT (WDISP);
C     * AND SOIL WATER OVERTAKEN BY DISPLACED SOIL WATER (WABS).
C     * RE-EVALUATE POND TEMPERATURE AND POND DEPTH; UPDATE WATER
C     * MOVEMENT MATRIX; ADJUST CURRENT POSITION OF WETTING FRONT.
C
      DO 400 I=IL1,IL2
          IF(ITER(I).EQ.1 .AND. ISIMP(I).EQ.-2)                     THEN
              IF(LZF(I).LT.IGP1)                            THEN  
                  IF(ZF(I).GT.1.0E-7)                THEN
                      FINF=GRKINF(I,LZF(I))*(ZF(I)+ZPOND(I)+
     1                     PSIF(I,LZF(I)))/ZF(I)    
                  ELSE                                                    
                      FINF=GRKINF(I,1)                                      
                  ENDIF                                                   
              ELSE
                  FINF=GRKINF(I,LZF(I))*(ZF(I)+ZPOND(I))/ZF(I)                 
              ENDIF
              IF(ZPOND(I).LT.1.0E-8 .AND. FINF.GT.R(I))      FINF=R(I)
              IF(FINF.GT.FMAX(I)) FINF=FMAX(I)                              
              WINF=FINF*DTFLOW(I)                                            
              IF(LZF(I).LT.IGP1)                            THEN
                  DZF(I)=WINF/THLNLZ(I)                                    
                  WDISP(I)=DZF(I)*THLQLZ(I)                                   
                  DZDISP(I)=WDISP(I)/(THLNLZ(I)-THLQLZ(I))                  
                  WABS(I)=DZDISP(I)*THLQLZ(I)                                 
                  IF((ZF(I)+DZF(I)+DZDISP(I)).GT.ZBOTX(I,LZF(I))) THEN
                     DTFLOW(I)=(ZBOTX(I,LZF(I))-ZF(I))/
     1                         (FINF/THLNLZ(I)+
     2                         (FINF*THLQLZ(I))/
     3                         (THLNLZ(I)*               
     4                         (THLNLZ(I)-THLQLZ(I))))                     
                     WINF=FINF*DTFLOW(I)          
                     DZF(I)=WINF/THLNLZ(I)                                
                     WDISP(I)=DZF(I)*THLQLZ(I)                               
                     DZDISP(I)=WDISP(I)/(THLNLZ(I)-THLQLZ(I))              
                     WABS(I)=DZDISP(I)*THLQLZ(I)                             
                  ENDIF                                                   
              ENDIF
              IF(ZPOND(I)+R(I)*DTFLOW(I).GT.1.0E-8)          THEN
                  TPOND(I)=(ZPOND(I)*TPOND(I)+R(I)*DTFLOW(I)*TR(I))/
     1                     (ZPOND(I)+R(I)*DTFLOW(I))
              ENDIF
              IF(ABS(R(I)-FINF-EVAP(I)).GT.1.0E-11)          THEN 
                  ZPTEST=ZPOND(I)+R(I)*DTFLOW(I)-WINF-EVAP(I)*
     1                   DTFLOW(I)                      
              ELSE                                                    
                  ZPTEST=ZPOND(I)                                        
              ENDIF                                                   
              IF(ZPTEST.LT.0.)                               THEN  
                  DTFLOW(I)=ZPOND(I)/(FINF+EVAP(I)-R(I))      
                  WINF=FINF*DTFLOW(I)                                        
                  ZPOND(I)=0.0                                           
                  IF(LZF(I).LT.IGP1)                    THEN
                     DZF(I)=WINF/THLNLZ(I)                                
                     WDISP(I)=DZF(I)*THLQLZ(I)                               
                     DZDISP(I)=WDISP(I)/(THLNLZ(I)-THLQLZ(I))              
                     WABS(I)=DZDISP(I)*THLQLZ(I)                             
                  ENDIF
              ELSE                                                    
                  ZPOND(I)=ZPTEST                                        
              ENDIF                                                   
              IF((WMOVE(I,1)+WINF).GT.0.)                    THEN 
                 TMOVE(I,1)=(WMOVE(I,1)*TMOVE(I,1)+WINF*TPOND(I))/            
     1                      (WMOVE(I,1)+WINF)                                 
              ENDIF                                                   
              WMOVE(I,1)=WMOVE(I,1)+WINF                                  
          ENDIF
  400 CONTINUE
C
C     * (THIS PORTION OF THE ABOVE DO-LOOP WAS SPLIT OFF ON THE CRAY
C     * BECAUSE IT WOULD NOT VECTORIZE. ONE MIGHT TRY AND RE-COMBINE 
C     * IT ON THE SX-3 (GOES IN FIRST PART OF IF BLOCK)).
C
      DO 450 I=IL1,IL2
          IF(ITER(I).EQ.1 .AND. ISIMP(I).EQ.-2 .AND.
     1       LZF(I).LT.IGP1)                                        THEN   
              WMOVE(I,NINF(I))=WMOVE(I,NINF(I))+WDISP(I)+WABS(I)                      
              ZF(I)=ZF(I)+DZF(I)+DZDISP(I)
          ENDIF
  450 CONTINUE  
C
C     In the 500 loop, if ISIMP < 0 (i.e. water movement has occurred), 
C     the time remaining in the current time step is recalculated. If 
C     the wetting front is at a soil layer boundary, if the time 
C     remaining is non-zero, if the wetting front is still within the 
C     modelled soil column, and if the calculated water content behind 
C     the wetting front is greater than zero, FMAX is updated, and LZF 
C     and NINF are incremented by 1. The NINF level of the TMOVE matrix 
C     is set to the water temperature of the new soil layer.
C
C     * CALCULATE REMAINING ITERATION TIME; RE-EVALUATE INFILTRATION
C     * PARAMETERS.
C
      DO 500 I=IL1,IL2
          IF(ITER(I).EQ.1 .AND. ISIMP(I).NE.1)                      THEN
              TRMDR(I)=TRMDR(I)-DTFLOW(I)                                                  
              IF(ABS(ZF(I)-ZBOTX(I,LZF(I))).LT.1.0E-6 .AND.
     1                TRMDR(I).GT.0. .AND. LZF(I).LT.IGP1 .AND.
     1                THLQLZ(I).GT.0.0001)                       THEN             
                  FINF=GRKINF(I,LZF(I))*(ZBOTX(I,LZF(I))+ZPOND(I)+
     1                 PSIF(I,LZF(I)))/ZBOTX(I,LZF(I))
                  FMAX(I)=MIN(FMAX(I),FINF)                                              
                  LZF(I)=LZF(I)+1                                                   
                  NINF(I)=NINF(I)+1                                                 
                  TMOVE(I,NINF(I))=TBARWX(I,LZF(I))                                     
              ENDIF                                                           
          ENDIF                                                                   
  500 CONTINUE
C
C     * INCREMENT ITERATION COUNTER ("NIT") AND SEE IF ANY POINTS STILL
C     * REMAIN TO BE DONE (USING "NPNTS"). IF SO, RETURN TO BEGINNING 
C     * TO COMPLETE THESE REMAINING POINTS.
C
      NIT=NIT+1
C
C   At the end of the iteration pass, checks are done in the 600 loop to 
C   ascertain whether the number of iteration passes is still less than 
C   NEND, whether either ponded water still exists or rain is still 
C   falling, and whether there is still time remaining in the current 
C   time step. If these conditions are all fulfilled, the counter NPNTS 
C   representing the number of points in the current vector of mosaic 
C   tiles for which infiltration is still occurring is incremented by 1. 
C   Otherwise, the iteration flag for the current tile is changed from 1 
C   to 0, signaling the end of the saturated infiltration calculations 
C   for that tile.
C
      DO 600 I=IL1,IL2
          IF(IGRN(I).GT.0)                                          THEN 
              IF(NIT.LE.NEND(I) .AND. ITER(I).EQ.1 .AND.
     1           (ZPOND(I).GT.1.0E-8 .OR. R(I).GT.0.)  .AND.
     2           TRMDR(I).GT.0.)                           THEN
                  NPNTS=NPNTS+1
              ELSE
                  ITER(I)=0
              ENDIF
          ENDIF
  600 CONTINUE
C
      IF(NPNTS.GT.0)                                         GO TO 100                       
C                                                                                  
      RETURN                                                                      
      END        
