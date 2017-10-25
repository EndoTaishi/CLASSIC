!>\file
C!Diagnoses snowpack visible and near-IR albedos given the
C!all-wave albedo at the current time step. Calculates snowpack
C!transmissivity for shortwave radiation.
C!
      SUBROUTINE SNOALBA(ALVSSN,ALIRSN,ALVSSC,ALIRSC,ALBSNO,
     1                   TRSNOWC, ALSNO, TRSNOWG, FSDB, FSFB, RHOSNO,   
     2                   REFSN,BCSN,SNO,CSZ,ZSNOW,FSNOW,ASVDAT,ASIDAT,  
     3                   ALVSG, ALIRG,                                  
     4                   ILG,IG,IL1,IL2,JL,IALS,NBS,ISNOALB)            
C
C     * JAN 27/16 - D.VERSEGHY. REFINE CALCULATIONS OF ALVSSN AND ALIRSN.
C     * NOV 16/13 - J.COLE.     Final version for gcm17:                
C     *                         - Fixes to get the proper BC mixing ratio in 
C     *                           snow, which required passing in and using  
C     *                           the snow density RHON.                
C     * JUN 22/13 - J.COLE.     ADD CODE FOR "ISNOALB" OPTION,          
C     *                         WHICH IS BASED ON 4-BAND SOLAR.      
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
      INTEGER ILG,IG,IL1,IL2,JL,IALS,IPTBAD,I,IB

      INTEGER NBS       !<Number of modelled shortwave radiation wavelength bands
      INTEGER ISNOALB   !< Switch to model snow albedo in two or more wavelength bands
C
C     * OUTPUT ARRAYS.
C
      REAL ALVSSN(ILG)  !<Visible albedo of snow pack on bare ground \f$[ ] 
                        !<(alpha_{s,VIS})\f$
      REAL ALIRSN(ILG)  !<Near-IR albedo of snow pack on bare ground \f$[ ] 
                        !<(alpha_{s,NIR})\f$
      REAL ALVSSC(ILG)  !<Visible albedo of snow on ground under 
                        !<vegetation canopy [ ]
      REAL ALIRSC(ILG)  !<Near-IR albedo of snow on ground under 
                        !<vegetation canopy [ ]
      REAL TRSNOWC(ILG) !<Transmissivity of snow under vegetation to shortwave radiation
                        !<\f$[ ] (\tau_{s,c})\f$
      REAL ALSNO (ILG,NBS) !<Albedo of snow in each modelled wavelength band  [  ]
      REAL TRSNOWG(ILG,NBS)!<Transmissivity of snow in bare areas to shortwave radiation \f$[ ] (\tau_{s,g})\f$
      REAL ALVSG  (ILG)    !<Near-IR albedo f bare ground  [  ]
      REAL ALIRG  (ILG)    !<Visible albedo of bare ground  [  ]
      

C
C     * INPUT ARRAYS.
C
      REAL FSDB(ILG,NBS)!<Direct solar radiation in each modelled wavelength band \f$[W m^{-2}]\f$
      REAL FSFB(ILG,NBS)!<Diffuse solar radiation in each modelled wavelength band \f$[W m^{-2}]\f$
      REAL ALBSNO(ILG)  !<All-wave albedo of snow pack \f$[ ] (\alpha_{s,T})\f$  
      REAL ZSNOW (ILG)  !<Depth of snow \f$[m] (z_s)\f$
      REAL FSNOW (ILG)  !<Fractional coverage of snow on grid cell [ ]
      REAL ASVDAT(ILG)  !<Assigned value of visible albedo of snow pack – 
                        !<optional [ ]
      REAL ASIDAT(ILG)  !<Assigned value of near-IR albedo of snow pack – 
                        !<optional [ ]
      REAL REFSN (ILG)  !<Snow grain size  [m]
      REAL BCSN  (ILG)  !<Black carbon mixing ratio \f$[kg m^{-3}]\f$
      REAL CSZ   (ILG)  !<Cosine of solar zenith angle  [  ]
      REAL SNO   (ILG)  !<Mass of snow pack \f$[kg m^{-2}]\f$
      REAL RHOSNO(ILG)  !<Density of snow pack  \f$[kg m^{-3}]\f$

C
C     * LOCAL ARRAYS                                                    
C                                                                       
      REAL SALBG(ILG,NBS), ALDIR(ILG,NBS), ALDIF(ILG,NBS),              
     +                     TRDIR(ILG,NBS), TRDIF(ILG,NBS)               
      REAL REFSNO(ILG), BCSNO(ILG)                                      
      INTEGER C_FLAG(ILG)                                               
C                                                                       
C     * CONSTANTS.                                                      
C                                                                       
      REAL WDIRCT, WDIFF                                                
      INTEGER SUM_C_FLAG                                                
C------------------------------------------------------------------
      !>
      !!In subroutine SNOALBW, called at the end of CLASSW, the change of 
      !!total snow albedo over the current time step is calculated using 
      !!an empirical exponential decay function, which has different 
      !!coefficients depending on whether the snow is dry or melting. In 
      !!this subroutine, if the ISNOALB switch is set to 0, the visible and
      !!near-IR components of the snow albedo are diagnosed from the total albedo.
      !!According to the literature (Aguado, 1985 \cite Aguado1985-fv ; Robinson and Kukla, 1984; Dirmhirn and
      !!Eaton, 1975 \cite Dirmhirn1975-vx), the following represent typical snow albedos for 
      !!fresh snow, old dry snow and melting snow:
      !!
      !!\f[
      !!\begin{tabular} { | l | c | c | c | }
      !!\hline
      !!              & Total albedo  & Visible albedo & Near-IR albedo \\ \hline
      !! Fresh snow   &      0.84     &      0.95      &      0.73      \\ \hline
      !! Old dry snow &      0.70     &      0.84      &      0.56      \\ \hline
      !! Melting snow &      0.50     &      0.62      &      0.38      \\ \hline
      !!\end{tabular}
      !!\f]
      !!
      !!The same decay function is assumed to apply to all three albedo 
      !!ranges, so the relative location of the visible and near-IR 
      !!albedos, \f$\alpha_{s,VIS}\f$ and \f$\alpha_{s,NIR}\f$, on the decay curve will be analogous 
      !!to that of the total albedo, \f$\alpha_{s,T}\f$. Thus, for dry snow:
      !!
      !!\f$[\alpha_{s,VIS} - 0.84]/[0.95-0.84] = [\alpha_{s,T} - 0.70]/[0.84-0.70]\f$
      !!\f$[\alpha_{s,NIR} - 0.56]/[0.73-0.56] = [\alpha_{s,T} - 0.70]/[0.84-0.70]\f$
      !!
      !!or, simplifying:
      !!
      !!\f$\alpha_{s,VIS} = 0.7857 \alpha_{s,T} + 0.2900\f$
      !!\f$\alpha_{s,NIR} = 1.2142 \alpha_{s,T} - 0.2900\f$
      !!
      !!For melting snow:
      !!
      !![\f$\alpha_{s,VIS} - 0.62]/[0.95-0.62] = [\alpha_{s,T} - 0.50]/[0.84-0.50]\f$
      !![\f$\alpha_{s,NIR} - 0.38]/[0.73-0.38] = [\alpha_{s,T} - 0.50]/[0.84-0.50]\f$
      !!
      !!or, simplifying:
      !!
      !!\f$\alpha_{s,VIS} = 0.9706 \alpha_{s,T} + 0.1347\f$
      !!\f$\alpha_{s,NIR} = 1.0294 \alpha_{s,T} - 0.1347\f$
      !!
      !!The above calculations are performed if the flag IALS is set to 
      !!zero. If IALS is set to one, indicating that assigned snow 
      !!albedos are to be used instead of calculated values, \f$\alpha_{s,VIS}\f$ and 
      !!\f$\alpha_{s,NIR}\f$ are set to the assigned values ASVDAT and ASIDAT 
      !!respectively. The sub-canopy values of visible and near-IR albedo 
      !!are currently set equal to the open snowpack values (this is 
      !!expected to change if a canopy litterfall parametrization is 
      !!developed).
      !!
      !!The transmissivity of snow under vegetation \f$\tau_{s,c}\f$ is
      !! then calculated from the snow depth
      !!ZSNOW using Beer’s law, with an empirical extinction coefficient 
      !!of \f$25.0 m^{-1}\f$ derived from the literature (Grenfell and Maykut, 
      !!1977 \cite Grenfell1977-pi ; Thomas, 1963):
      !!
      !!\f$\tau_{s,c} = exp[-25.0 z_s]\f$
      !!
      !!If the ISNOALB switch is set to zero, the value of ALSNO in the first
      !!wavelength band is set to the previously calculated value of \f$\alpha_{s,VIS}\f$
      !!and the values for the remaining bands are set to the previously calculated value
      !! of \f$\alpha_{s,NIR}\f$; the value of \f$\tau_{s,g}\f$  is set to \f$\tau_{s,c}\f$.
      !!If the ISNOALB switch is set to 1, a new parameterization for the snow albedo and
      !!transmissivity in four shortwave radiation bands (one visible and three near-IR)
      !!is used, according to Cole et al. (2017).  This parameterization incorporates
      !!the effects of snow grain size and black carbon content, and makes use of lookup
      !!tables contained in the CCCma subroutines SNOW_ALBVAL and SNOW_TRANVAL.
      IPTBAD=0
      DO 100 I=IL1,IL2                                           
         IF(ALBSNO(I).LT.0.50.AND.ALBSNO(I).GT.0.499) ALBSNO(I)=0.50                      
         IF(FSNOW(I).GT.0.0 .AND. IALS.EQ.0)              THEN  
             IF(ALBSNO(I).GT.0.70)                    THEN
                 ALVSSN(I)=0.7857*ALBSNO(I)+0.2900
                 ALIRSN(I)=1.2142*ALBSNO(I)-0.2900
             ELSE
                 ALVSSN(I)=0.9706*ALBSNO(I)+0.1347
                 ALIRSN(I)=1.0294*ALBSNO(I)-0.1347
             ENDIF
             IF(ALVSSN(I).GT.0.999.OR.ALVSSN(I).LT.0.001) IPTBAD=I
             IF(ALIRSN(I).GT.0.999.OR.ALIRSN(I).LT.0.001) IPTBAD=I
         ELSE IF(FSNOW(I).GT.0.0 .AND. IALS.EQ.1)         THEN  
             ALVSSN(I)=ASVDAT(I)
             ALIRSN(I)=ASIDAT(I)
         ENDIF                                                                   
         ALVSSC(I)=ALVSSN(I)
         ALIRSC(I)=ALIRSN(I)
         TRSNOWC(I)=EXP(-25.0*ZSNOW(I))                                                 
  100 CONTINUE
C
      IF(IPTBAD.NE.0) THEN
         WRITE(6,6100) IPTBAD,JL,ALVSSN(IPTBAD),ALIRSN(IPTBAD)
 6100    FORMAT('0AT (I,J)= (',I3,',',I3,'), ALVSSN,ALIRSN = ',2F10.5)
         CALL XIT('SNOALBA',-1)
      ENDIF
C
      IF (ISNOALB .EQ. 0) THEN
         DO I = IL1, IL2                                                
            ALSNO(I,1) = ALVSSN(I)                                      
            ALSNO(I,2) = ALIRSN(I)                                      
            ALSNO(I,3) = ALIRSN(I)                                      
            ALSNO(I,4) = ALIRSN(I)

            TRSNOWG(I,1:NBS) = TRSNOWC(I)                               
         END DO ! I                                                     
      ELSE IF (ISNOALB .EQ. 1) THEN                                     
         DO IB = 1, NBS                                                 
            DO I = IL1, IL2                                             
               IF (IB .EQ. 1) THEN                                      
                  SALBG(I,IB) = ALVSG(I)                                
                  ALSNO(I,IB) = ALVSSN(I)                               
               ELSE                                                     
                  SALBG(I,IB) = ALIRG(I)                                
                  ALSNO(I,IB) = ALIRSN(I)                               
               END IF                                                   
            END DO ! I                                                  
         END DO ! IB                                                    
         SUM_C_FLAG = 0                                                 
         DO I = IL1, IL2                                                
            IF (ZSNOW(I) .GT. 0.0) THEN                                 
               C_FLAG(I) = 1                                            
            ELSE                                                        
               C_FLAG(I) = 0                                            
            END IF                                                      
            SUM_C_FLAG = SUM_C_FLAG + C_FLAG(I)                         
         END DO ! I                                                     
                                                                        
         IF (IALS .EQ. 0) THEN                                          
            IF (SUM_C_FLAG .GT. 0) THEN                                 
!>
!! Convert the units of the snow grain size and BC mixing ratio          
!! Snow grain size from meters to microns and BC from \f$kg BC/m^3\f$ to ng BC/kg SNOW 
!!
               DO I = IL1,IL2                                           
                 IF (C_FLAG(I) .EQ. 1) THEN                             
                  REFSNO(I) = REFSN(I)*1.0E6                            
                  BCSNO(I)  = (BCSN(I)/RHOSNO(I))*1.0E12                
                 END IF                                                 
               END DO ! I                                               
                                                                        
               CALL SNOW_ALBVAL(ALDIF, ! OUTPUT                         
     +                          ALDIR,                                  
     +                          CSZ,   ! INPUT                          
     +                          SALBG,                                  
     +                          BCSNO,                                  
     +                          REFSNO,                                 
     +                          SNO,                                    
     +                          C_FLAG,                                 
     +                          IL1,                                    
     +                          IL2,                                    
     +                          ILG,                                    
     +                          NBS)                                    
                                                                        
               CALL SNOW_TRANVAL(TRDIF, ! OUTPUT                        
     +                           TRDIR,                                 
     +                           CSZ,   ! INPUT                         
     +                           SALBG,                                 
     +                           BCSNO,                                 
     +                           REFSNO,                                
     +                           SNO,                                   
     +                           C_FLAG,                                
     +                           IL1,                                   
     +                           IL2,                                   
     +                           ILG,                                   
     +                           NBS)                                   
                                                                        
               DO IB = 1, NBS                                           
                  DO I = IL1, IL2                                       
                     IF (C_FLAG(I) .EQ. 1) THEN                         
                        WDIRCT = FSDB(I,IB)                             
     +                         /(FSDB(I,IB)+FSFB(I,IB)+1.E-10)          
                        WDIFF  = 1.0-WDIRCT                             
                        ALSNO(I,IB) = ALDIF(I,IB)*WDIFF                 
     +                              + ALDIR(I,IB)*WDIRCT                
                        TRSNOWG(I,IB) = TRDIF(I,IB)*WDIFF               
     +                                + TRDIR(I,IB)*WDIRCT              
                     END IF ! C_FLAG                                    
                  END DO ! I                                            
               END DO ! IB                                              
            ELSE ! SUM_C_FLAG .EQ. 0                                    
               DO I = IL1, IL2                                          
                  ALSNO(I,1)     = ALVSSN(I)                            
                  ALSNO(I,2:NBS) = ALIRSN(I)                            
                  TRSNOWG(I,1:NBS) = TRSNOWC(I)                         
               END DO ! I                                               
            ENDIF ! SUM_C_FLAG                                          
         ELSE IF (IALS .EQ. 1) THEN                                     
            DO I = IL1, IL2                                             
               ALSNO(I,1) = ASVDAT(I)                                   
               ALSNO(I,2) = ASIDAT(I)                                   
               ALSNO(I,3) = ASIDAT(I)                                   
               ALSNO(I,4) = ASIDAT(I)                                   
                                                                        
               TRSNOWG(I,1:NBS) = TRSNOWC(I)                            
            END DO ! I                                                  
         END IF ! IALS                                                  
      END IF ! ISNOALB                                                  
      RETURN
      END
