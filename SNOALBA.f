      SUBROUTINE SNOALBA(ALVSSN,ALIRSN,ALVSSC,ALIRSC,ALBSNO,
     1                   TRSNOW,ZSNOW,FSNOW,ASVDAT,ASIDAT,
     2                   ILG,IG,IL1,IL2,JL,IALS)
C
C     Purpose: Diagnose snowpack visible and near-IR albedos given the 
C     all-wave albedo at the current time step. Calculate snowpack 
C     transmissivity for shortwave radiation.
C
C     * FEB 05/07 - D.VERSEGHY. STREAMLINE CALCULATIONS OF
C     *                         ALVSSN AND ALIRSN.
C     * APR 13/06 - D.VERSEGHY. SEPARATE ALBEDOS FOR OPEN AND 
C     *                         CANOPY-COVERED SNOW.
C     * NOV 03/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
C     * MAR 18/02 - D.VERSEGHY. UPDATES TO ALLOW ASSIGNMENT OF
C     *                         USER-SPECIFIED VALUES TO SNOW
C     *                         ALBEDO.
C     * JUN 05/97 - D.VERSEGHY. CLASS - VERSION 2.7.
C     *                         SPECIFY LOCATION OF ICE SHEETS
C     *                         BY SOIL TEXTURE ARRAY RATHER
C     *                         THAN BY SOIL COLOUR INDEX.
C     * NOV 29/94 - M.LAZARE.   CLASS - VERSION 2.3.
C     *                         CALL ABORT CHANGED TO CALL XIT TO 
C     *                         ENABLE RUNNING ON PC'S.
C     * MAR 13/92 - M.LAZARE.   CODE FOR MODEL VERSION GCM7 -
C     *                         DIVIDE PREVIOUS SUBROUTINE 
C     *                         "SNOALB" INTO "SNOALBA" AND
C     *                         "SNOALBW" AND VECTORIZE.
C     * AUG 12/91 - D.VERSEGHY. CODE FOR MODEL VERSION GCM7U -
C     *                         CLASS VERSION 2.0 (WITH CANOPY).
C     * APR 11/89 - D.VERSEGHY. DISAGGREGATE SNOW ALBEDO INTO
C     *                         VISIBLE AND NEAR-IR PORTIONS;
C     *                         CALCULATE TRANSMISSIVITY TO
C     *                         SHORTWAVE RADIATION.
C
      IMPLICIT NONE
C    
C     * INTEGER CONSTANTS.
C
      INTEGER ILG,IG,IL1,IL2,JL,IALS,IPTBAD,I
C
C     * OUTPUT ARRAYS.
C
      REAL ALVSSN(ILG)  ! Visible albedo of snow pack on bare ground [ ] 
                        !(alpha_s,VIS)
      REAL ALIRSN(ILG)  !Near-IR albedo of snow pack on bare ground [ ] 
                        !(alpha_s,NIR)
      REAL ALVSSC(ILG)  !Visible albedo of snow on ground under 
                        !vegetation canopy [ ]
      REAL ALIRSC(ILG)  !Near-IR albedo of snow on ground under 
                        !vegetation canopy [ ]
      REAL TRSNOW(ILG)  !Transmissivity of snow to shortwave radiation 
                        ![ ] (tau_s)

C
C     * INPUT ARRAYS.
C
      REAL ALBSNO(ILG)  !All-wave albedo of snow pack [ ] (alpha_s,T)  
      REAL ZSNOW (ILG)  !Depth of snow [m] (zs )
      REAL FSNOW (ILG)  !Fractional coverage of snow on grid cell [ ]
      REAL ASVDAT(ILG)  !Assigned value of visible albedo of snow pack – 
                        !optional [ ]
      REAL ASIDAT(ILG)  !Assigned value of near-IR albedo of snow pack – 
                        !optional [ ]
C
C------------------------------------------------------------------
      !
      !In subroutine SNOALBW, called at the end of CLASSW, the change of 
      !total snow albedo over the current time step is calculated using 
      !an empirical exponential decay function, which has different 
      !coefficients depending on whether the snow is dry or melting. In 
      !this subroutine, the visible and near-IR components of the snow 
      !albedo are diagnosed from the total albedo. According to the 
      !literature (Aguado, 1985; Robinson and Kukla, 1984; Dirmhirn and 
      !Eaton, 1975), the following represent typical snow albedos for 
      !fresh snow, old dry snow and melting snow:
      !
      !----------------------------------------------------------------|
      !              | Total albedo  | Visible albedo | Near-IR albedo |
      !----------------------------------------------------------------|
      ! Fresh snow   |      0.84     |      0.95      |      0.73      |
      !----------------------------------------------------------------|
      ! Old dry snow |      0.70     |      0.84      |      0.56      |
      !----------------------------------------------------------------|
      ! Melting snow |      0.50     |      0.62      |      0.38      |
      !----------------------------------------------------------------|
      !
      !The same decay function is assumed to apply to all three albedo 
      !ranges, so the relative location of the visible and near-IR 
      !albedos, ALVSSN and ALIRSN, on the decay curve will be analogous 
      !to that of the total albedo, ALBSNO. Thus, for dry snow:
      !
      ![ALVSSN - 0.84]/[0.95-0.84] = [ALBSNO - 0.70]/[0.84-0.70]
      ![ALIRSN - 0.56]/[0.73-0.56] = [ALBSNO - 0.70]/[0.84-0.70]
      !
      !or, simplifying:
      !
      !ALVSSN = 0.79[ALBSNO - 0.70] + 0.84
      !ALIRSN = 1.21[ALBSNO - 0.70] + 0.56
      !
      !For melting snow:
      !
      ![ALVSSN - 0.62]/[0.95-0.62] = [ALBSNO - 0.50]/[0.84-0.50]
      ![ALIRSN - 0.38]/[0.73-0.38] = [ALBSNO - 0.50]/[0.84-0.50]
      !
      !or, simplifying:
      !
      !ALVSSN = 0.97[ALBSNO - 0.50] + 0.62
      !ALIRSN = 1.03[ALBSNO - 0.50] + 0.38
      !
      !The above calculations are performed if the flag IALS is set to 
      !zero. If IALS is set to one, indicating that assigned snow 
      !albedos are to be used instead of calculated values, ALVSSN and 
      !ALIRSN are set to the assigned values ASVDAT and ASIDAT 
      !respectively. The sub-canopy values of visible and near-IR albedo 
      !are currently set equal to the open snowpack values (this is 
      !expected to change if a canopy litterfall parametrization is 
      !developed).
      !
      !The transmissivity of snow τs is calculated from the snow depth 
      !ZSNOW using Beer’s law, with an empirical extinction coefficient 
      !of 25.0 m-1 derived from the literature (Grenfell and Maykut, 
      !1977; Thomas, 1963):
      !
      !TRSNOW = exp[-25.0*ZSNOW]
      !
      IPTBAD=0
      DO 100 I=IL1,IL2                                           
         IF(ALBSNO(I).LT.0.50.AND.ALBSNO(I).GT.0.499) ALBSNO(I)=0.50                      
         IF(FSNOW(I).GT.0.0 .AND. IALS.EQ.0)              THEN  
             IF(ALBSNO(I).GT.0.70)                    THEN
                 ALVSSN(I)=0.79*(ALBSNO(I)-0.70)+0.84                                         
                 ALIRSN(I)=1.21*(ALBSNO(I)-0.70)+0.56                                         
             ELSE
                 ALVSSN(I)=0.97*(ALBSNO(I)-0.50)+0.62                                         
                 ALIRSN(I)=1.03*(ALBSNO(I)-0.50)+0.38                                         
             ENDIF
             IF(ALVSSN(I).GT.0.999.OR.ALVSSN(I).LT.0.001) IPTBAD=I
             IF(ALIRSN(I).GT.0.999.OR.ALIRSN(I).LT.0.001) IPTBAD=I
         ELSE IF(FSNOW(I).GT.0.0 .AND. IALS.EQ.1)         THEN  
             ALVSSN(I)=ASVDAT(I)
             ALIRSN(I)=ASIDAT(I)
         ENDIF                                                                   
         ALVSSC(I)=ALVSSN(I)
         ALIRSC(I)=ALIRSN(I)
         TRSNOW(I)=EXP(-25.0*ZSNOW(I))                                                 
  100 CONTINUE
C
      IF(IPTBAD.NE.0) THEN
         WRITE(6,6100) IPTBAD,JL,ALVSSN(IPTBAD),ALIRSN(IPTBAD)
 6100    FORMAT('0AT (I,J)= (',I3,',',I3,'), ALVSSN,ALIRSN = ',2F10.5)
         CALL XIT('SNOALBA',-1)
      ENDIF
C
      RETURN
      END
