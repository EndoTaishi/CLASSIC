!>\file
!!Calculates visible and near-IR ground albedos.
!!@author D. Verseghy, M. Lazare
!
      SUBROUTINE GRALB(ALVSG,ALIRG,ALVSGC,ALIRGC,
     1                 ALGWV,ALGWN,ALGDV,ALGDN,
     2                 THLIQ,FSNOW,ALVSU,ALIRU,FCMXU,                   
     3                 AGVDAT,AGIDAT,FG,ISAND,  
     4                 ILG,IG,IL1,IL2,JL,IALG)        
C
C     * DEC 15/16 - D.VERSEGHY. ASSIGN ROCK ALBEDO USING SOIL COLOUR
C     *                         INDEX INSTEAD OF VIA LOOKUP TABLE.
C     * JAN 16/15 - D.VERSEGHY. CORRECT ACCOUNTING FOR URBAN ALBEDO.
C     * FEB 09/15 - D.VERSEGHY. New version for gcm18 and class 3.6:
C     *                         - Wet and dry albedoes for EACH of
C     *                           visible and near-ir are passed in
C     *                           instead of ALGWET and ALGDRY. These
C     *                           are used to calculate ALISG and ALIRG.
C     * NOV 16/13 - M.LAZARE.   FINAL VERSION FOR GCM17:                
C     *                         - REMOVE UNNECESSARY LOWER BOUND        
C     *                           OF 1.E-5 ON "FURB".                   
C     * APR 13/06 - D.VERSEGHY. SEPARATE ALBEDOS FOR OPEN AND
C     *                         CANOPY-COVERED GROUND.
C     * NOV 03/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
C     * SEP 04/03 - D.VERSEGHY. RATIONALIZE CALCULATION OF URBAN
C     *                         ALBEDO.
C     * MAR 18/02 - D.VERSEGHY. UPDATES TO ALLOW ASSIGNMENT OF USER-
C     *                         SPECIFIED VALUES TO GROUND ALBEDO.
C     *                         PASS IN ICE AND ORGANIC ALBEDOS
C     *                         VIA NEW COMMON BLOCK "CLASS8".
C     * FEB 07/02 - D.VERSEGHY. REVISIONS TO BARE SOIL ALBEDO
C     *                         CALCULATIONS; REMOVAL OF SOIL
C     *                         MOISTURE EXTRAPOLATION TO SURFACE.
C     * JUN 05/97 - D.VERSEGHY. CLASS - VERSION 2.7.
C     *                         CALCULATE SOIL ALBEDO FROM PERCENT
C     *                         SAND CONTENT RATHER THAN FROM COLOUR 
C     *                         INDEX.
C     * SEP 27/96 - D.VERSEGHY. CLASS - VERSION 2.5.
C     *                         FIX BUG TO CALCULATE GROUND ALBEDO
C     *                         UNDER CANOPIES AS WELL AS OVER BARE
C     *                         SOIL.
C     * NOV 29/94 - M.LAZARE.   CLASS - VERSION 2.3.
C     *                         "CALL ABORT" CHANGED TO "CALL XIT",
C     *                         TO ENABLE RUNNING RUN ON PC'S.
C     * FEB 12/93 - D.VERSEGHY/M.LAZARE. INCREASE DRY SOIL ALBEDO TO 
C     *                                  0.45 FROM 0.35. 
C     * MAR 13/92 - D.VERSEGHY/M.LAZARE. REVISED AND VECTORIZED CODE
C     *                                  FOR MODEL VERSION GCM7.
C     * AUG 12/91 - D.VERSEGHY. CLASS - VERSION 2.0.
C     *                         CODE FOR MODEL VERSION GCM7U (WITH
C     *                         CANOPY). 
C     * APR 11/89 - D.VERSEGHY. CALCULATE VISIBLE AND NEAR-IR SOIL
C     *                         ALBEDOS BASED ON TEXTURE AND SURFACE
C     *                         WETNESS. (SET TO ICE ALBEDOS OVER
C     *                         CONTINENTAL ICE SHEETS.)
C
      use classic_params, only : ALVSI,ALIRI,ALVSO,ALIRO
      
      IMPLICIT NONE
C                
C     * INTEGER CONSTANTS.
C
      INTEGER  ILG,IG,IL1,IL2,JL,IALG,IPTBAD,I
C
C     * OUTPUT ARRAYS.
C
      REAL ALVSG (ILG)  !<Visible albedo of bare ground \f$[ ] (\alpha_{g,VIS})\f$   
      REAL ALIRG (ILG)  !<Near-IR albedo of bare ground \f$[ ] (\alpha_{g,NIR})\f$
      REAL ALVSGC (ILG) !<Visible albedo of ground under vegetation canopy [ ]
      REAL ALIRGC (ILG) !<Near-IR albedo of ground under vegetation canopy [ ]
C
C     * INPUT ARRAYS.
C
      REAL ALGWV (ILG)    !<Visible albedo of wet soil for modelled area \f$[ ] (alpha_{g,NIR,wet})\f$
      REAL ALGWN (ILG)    !<Near-infrared albedo of wet soil for modelled area  \f$[ ] (alpha_{g,NIR,wet})\f$
      REAL ALGDV (ILG)    !<Visible albedo of dry soil for modelled area \f$[ ] (alpha_{g,VIS,dry})\f$
      REAL ALGDN (ILG)    !<Near-infrared albedo of dry soil for modelled area  \f$[ ] (alpha_{g,NIR,dry})\f$
      REAL THLIQ (ILG,IG) !<Volumetric liquid water content of soil layers \f$[m^3 m^{-3}]\f$
      REAL ALVSU (ILG)    !<Visible albedo of urban part of modelled area \f$[ ] (alpha_{u,VIS})\f$
      REAL ALIRU (ILG)    !<Near-IR albedo of urban part of modelled area \f$[ ] (alpha_{u,NIR})\f$
      REAL FCMXU (ILG)    !<Fractional coverage of urban part of modelled area \f$[ ] (X_u)\f$
      REAL AGVDAT(ILG)    !<Assigned value of visible albedo of ground – optional [ ]
      REAL AGIDAT(ILG)    !<Assigned value of near-IR albedo of ground – optional [ ]
      REAL FG    (ILG)    !<Fractional coverage of bare soil on modelled area \f$[ ] (X_g)\f$
      REAL FSNOW (ILG)    !<Fractional coverage of snow on modelled area  [  ]
C
      INTEGER ISAND (ILG,IG)!<Soil type flag based on sand content, assigned in subroutine CLASSB
C
C     * TEMPORARY VARIABLES.
C
      REAL FURB,ALBSOL
C
C---------------------------------------------------------------------
      !>
      !!If the ISAND flag for the surface soil layer is greater than or
      !!equal to zero (indicating mineral soil), first the urban area not
      !!covered by snow is evaluated.  Next the visible and near-IR open ground
      !!albedos, \f$\alpha_{g,VIS}\f$ and \f$\alpha_{g,NIR}\f$, for each wavelength range, are calculated on the basis of the
      !!wet and dry ground albedos \f$\alpha_{g,wet}\f$ and \f$\alpha_{g,dry}\f$ which were assigned 
      !!for the modelled area in CLASSB. Idso et al. (1975) found a 
      !!correlation between the soil liquid moisture content in the top 
      !!10 cm of soil (represented in CLASS by that of the first soil 
      !!layer, \f$\theta_{1,1}\f$) and the total surface albedo \f$\alpha_{g,T}\f$: for
      !!water contents less than 0.22 \f$m^3 m^{-3}\f$, \f$\alpha_{g,T}\f$ took on the value of 
      !!\f$\alpha_{g,dry}\f$; for liquid water contents greater than 0.26 \f$m^3 m^{-3}\f$,
      !!\f$\alpha_{g,T}\f$ took on the value of \f$\alpha_{g,wet}\f$. For values of \f$\theta_{1,1}\f$ between 
      !!these two limits, a linear relationship is assumed:
      !!
      !!\f$[\alpha_{g,T} - \alpha_{g,dry} ] / [\theta_{l,1} - 0.22] = [\alpha_{g,wet} - \alpha_{g,dry} ]/[0.26 - 0.22]\f$
      !!
      !!Thus, in GRALB for each of the two wavelength ranges \f$\alpha_{g}\f$ is calculated as follows:
      !!
      !!\f$\alpha_{g,T} = \alpha_{g,dry}\f$            \f$\theta_{l,1} \leq 0.22 \f$ \n
      !!\f$\alpha_{g,T} = \theta_{l,1} [\alpha_{g,wet} - \alpha_{g,dry} ]/0.04 - 5.50 [\alpha_{g,wet} - \alpha_{g,dry} ] + \alpha_{g,dry} \f$
      !!\f$0.22 < \theta_{l,1} < 0.26 \f$ \n
      !!\f$\alpha_{g,T} = \alpha{g,wet} \f$             \f$  0.26 \leq \theta_{l,1}\f$
      !!
      !!Afterwards, a correction is applied to \f$\alpha_{g,VIS}\f$ and \f$\alpha_{g,NIR}\f$ in order to
      !!account for the possible presence of urban surfaces over the
      !!modelled area. Visible and near-IR albedos are assigned for local 
      !!urban areas, \f$\alpha_{g,VIS}\f$ and \f$\alpha_{g,NIR}\f$, as part of the background data 
      !!(see the section on “Data Requirements”). A weighted average over the bare soil area \f$X_g\f$ is
      !!calculated from the fractional snow-free urban area \f$X_u\f$ as:
      !!
      !!\f$\alpha_{g,VIS} = [X_u \alpha_{u,VIS} + (X_g-X_u ) \alpha_{g,VIS}] / X_g \f$\n
      !!\f$\alpha_{g,NIR} = [X_u \alpha_{u,NIR} + (1.0-X_u ) \alpha_{g,NIR}] / X_g \f$
      !!
      !!If the soil on the modelled area is not mineral, i.e. if the 
      !!ISAND flag is less than zero, \f$\alpha_{g,VIS}\f$ and \f$\alpha_{g,NIR}\f$ are determined as 
      !!follows:
      !!
      !!If ISAND = -2, indicating organic soil, \f$\alpha_{g,VIS}\f$ and \f$\alpha_{g,NIR}\f$ are 
      !!assigned values of 0.05 and 0.30 respectively from the lookup 
      !!tables in the block data subroutine CLASSBD, corresponding to 
      !!average measured values reported in Comer et al. (2000) \cite Comer2000-mz.
      !!
      !!If ISAND = -3, indicating rock at the surface, \f$\alpha_{g,VIS}\f$ and
      !!\f$\alpha_{g,NIR}\f$ are assigned the dry ground values from CLASSB.
      !!
      !!If ISAND = -4, indicating continental ice sheet or glacier, \f$\alpha_{g,VIS}\f$ 
      !!and \f$\alpha_{g,NIR}\f$ are assigned values of 0.95 and 073 from CLASSBD, 
      !!reflecting values reported for Antarctica (e.g. Sellers, 1974).
      !!
      !!The above calculations are all performed if the flag IALG is set 
      !!to zero. If IALG is set to one, indicating that assigned ground 
      !!albedos are to be used instead of calculated values, \f$\alpha_{g,VIS}\f$ and 
      !!\f$\alpha_{g,NIR}\f$ are set to the assigned values AGVDAT and AGIDAT 
      !!respectively.
      !!
      !!Lastly, the ground values of visible and near-IR albedo under the 
      !!vegetation canopy are currently set equal to the open values 
      !!(this approach is under review).
      !!
      IPTBAD=0
      DO 100 I=IL1,IL2
         IF(IALG.EQ.0)                                          THEN
            IF(ISAND(I,1).GE.0)                          THEN
                FURB=FCMXU(I)*(1.0-FSNOW(I))                        

                   IF(THLIQ(I,1).GE.0.26) THEN                      
                      ALVSG(I)=ALGWV(I)                            
                      ALIRG(I)=ALGWN(I)                           
                   ELSEIF(THLIQ(I,1).LE.0.22) THEN                   
                      ALVSG(I)=ALGDV(I)                              
                      ALIRG(I)=ALGDN(I)                              
                   ELSE                                              
                      ALVSG(I)=THLIQ(I,1)*(ALGWV(I)-ALGDV(I))/0.04+  
     1                       ALGDV(I)-5.50*(ALGWV(I)-ALGDV(I))       
                      ALIRG(I)=THLIQ(I,1)*(ALGWN(I)-ALGDN(I))/0.04+  
     1                       ALGDN(I)-5.50*(ALGWN(I)-ALGDN(I))       
                   ENDIF                                             
C                                                                       
                IF(FG(I).GT.0.001)                          THEN
                    ALVSG(I)=((FG(I)-FURB)*ALVSG(I)+FURB*ALVSU(I))/FG(I)
                    ALIRG(I)=((FG(I)-FURB)*ALIRG(I)+FURB*ALIRU(I))/FG(I)
                ENDIF
                IF(ALVSG(I).GT.1.0.OR.ALVSG(I).LT.0.0) IPTBAD=I
                IF(ALIRG(I).GT.1.0.OR.ALIRG(I).LT.0.0) IPTBAD=I

            ELSE IF(ISAND(I,1).EQ.-4)                    THEN
                ALVSG(I)=ALVSI
                ALIRG(I)=ALIRI
            ELSE IF(ISAND(I,1).EQ.-3)                    THEN
                    ALVSG(I)=ALGDV(I)
                    ALIRG(I)=ALGDN(I)
            ELSE IF(ISAND(I,1).EQ.-2)                    THEN
                ALVSG(I)=ALVSO
                ALIRG(I)=ALIRO
            ENDIF
         ELSEIF(IALG.EQ.1)                                     THEN
            ALVSG(I)=AGVDAT(I)
            ALIRG(I)=AGIDAT(I)
         ENDIF     
         ALVSGC(I)=ALVSG(I)
         ALIRGC(I)=ALIRG(I)
  100 CONTINUE
C
      IF(IPTBAD.NE.0) THEN
         WRITE(6,6100) IPTBAD,JL,ALVSG(IPTBAD),ALIRG(IPTBAD)
 6100    FORMAT('0AT (I,J)= (',I3,',',I3,'), ALVSG,ALIRG = ',2F10.5)
         CALL XIT('GRALB',-1)
      ENDIF

      RETURN                                                                      
      END
