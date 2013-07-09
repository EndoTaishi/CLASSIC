      SUBROUTINE TWCALC(TBAR,THLIQ,THICE,HCP,TBARW,HMFG,HTC,
     1                  FI,EVAP,THPOR,THLMIN,HCPS,DELZW,
     2                  DELZZ,ISAND,IG,ILG,IL1,IL2,JL)
C
C   Check for freezing or thawing of liquid or frozen water in the soil 
C   layers, and adjust layer temperatures and water stores accordingly.
C
C     * SEP 23/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
C     * MAY 16/03 - Y.DELAGE/D.VERSEGHY. BUGFIX IN FREEZING/
C     *                                  THAWING CALCULATIONS
C     *                                  (PRESENT SINCE V.2.7)
C     * JUL 26/02 - D.VERSEGHY. SHORTENED CLASS4 COMMON BLOCK.
C     * JUN 20/97 - D.VERSEGHY. COSMETIC REARRANGEMENT OF 
C     *                         SUBROUTINE CALL.
C     * JUN 20/97 - D.VERSEGHY. CLASS - VERSION 2.7.
C     *                         MODIFICATIONS TO ALLOW FOR VARIABLE
C     *                         SOIL PERMEABLE DEPTH.
C     * JAN 02/96 - D.VERSEGHY. CLASS - VERSION 2.5.
C     *                         COMPLETION OF ENERGY BALANCE
C     *                         DIAGNOSTICS.
C     * AUG 18/95 - D.VERSEGHY. CLASS - VERSION 2.4.
C     *                         REVISIONS TO ALLOW FOR INHOMOGENEITY
C     *                         BETWEEN SOIL LAYERS AND FRACTIONAL
C     *                         ORGANIC MATTER CONTENT. 
C     * DEC 22/94 - D.VERSEGHY. CLASS - VERSION 2.3.
C     *                         REVISE CALCULATION OF HTC.
C     * JUL 30/93 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.2.
C     *                                  NEW DIAGNOSTIC FIELDS.
C     * APR 24/92 - D.VERSEGHY/M.LAZARE. REVISED AND VECTORIZED CODE
C     *                                  FOR MODEL VERSION GCM7.
C     * AUG 12/91 - D.VERSEGHY. CODE FOR MODEL VERSION GCM7U - 
C     *                         CLASS VERSION 2.0 (WITH CANOPY).
C     * APR 11/89 - D.VERSEGHY. ADJUST SOIL LAYER TEMPERATURES
C     *                         AND LIQUID/FROZEN MOISTURE CONTENTS
C     *                         FOR FREEZING/THAWING.
C
      IMPLICIT NONE
C
C     * INTEGER CONSTANTS.
C
      INTEGER IG,ILG,IL1,IL2,JL,I,J
C
C     * INPUT/OUTPUT ARRAYS.
C
      REAL TBAR  (ILG,IG)   !Temperature of soil layer [C] (Tg)
      REAL THLIQ (ILG,IG)   !Volumetric liquid water content of soil layer [m3 m-3] (theta_l)
      REAL THICE (ILG,IG)   !Volumetric frozen water content of soil layer [m3 m-3] (theta_i)
      REAL HCP   (ILG,IG)   !Heat capacity of soil layer [J m-3 K-1] (Cg)
      REAL TBARW (ILG,IG)   !Temperature of water in soil layer [C]
      REAL HMFG  (ILG,IG)   !Energy associated with freezing or thawing of water in soil layer [W m-2]
      REAL HTC   (ILG,IG)   !Internal energy change of soil layer due to conduction and/or change in mass [W m-2] (Ig)
C
C     * INPUT ARRAYS.
C
      REAL FI    (ILG)  !Fractional coverage of subarea in question on modelled area [ ] (Xi)
      REAL EVAP  (ILG)  !Calculated evaporation rate from soil surface [m s-1]
C
      REAL THPOR (ILG,IG)   !Pore volume in soil layer [mm] (theta_p)
      REAL THLMIN(ILG,IG)   !Residual soil liquid water content remaining after freezing or evaporation [m3 m-3] (theta_r)
      REAL HCPS  (ILG,IG)   !Heat capacity of soil material [J m-3 K-1] (Cm)
      REAL DELZW (ILG,IG)   !Permeable thickness of soil layer [m] (delta_zg,w)
      REAL DELZZ (ILG,IG)   !Soil layer thicknesses to bottom of permeable depth for standard three-layer configuration, or to bottom of thermal depth for multiple layers [m] (delta_zg,z)
C
      INTEGER              ISAND (ILG,IG)   !Sand content flag
C
C     * TEMPORARY VARIABLES.
C
      REAL THFREZ,THEVAP,HADD,THMELT
C
C     * COMMON BLOCK PARAMETERS.
C
      REAL DELT     !Time step [s]
      REAL TFREZ    !Freezing point of water [K]
      REAL HCPW     !Volumetric heat capacity of water (4.187*10^6) [J m-3 K-1]
      REAL HCPICE   !Volumetric heat capacity of ice (1.9257*10^6) [J m-3 K-1]
      REAL HCPSOL   !Volumetric heat capacity of mineral matter (2.25*10^6) [J m-3 K-1]
      REAL HCPOM    !Volumetric heat capacity of organic matter (2.50*10^6) [J m-3 K-1]
      REAL HCPSND   !Volumetric heat capacity of sand particles (2.13*10^6) [J m-3 K-1]
      REAL HCPCLY   !Volumetric heat capacity of fine mineral particles (2.38*10^6) [J m-3 K-1]
      REAL SPHW     !Specific heat of water (4.186*10^3) [J kg-1 K-1]
      REAL SPHICE   !Specific heat of ice (2.10*10^3) [J kg-1 K-1]
      REAL SPHVEG   !Specific heat of vegetation matter (2.70*10^3) [J kg-1 K-1]
      REAL SPHAIR   !Specific heat of air [J kg-1 K-1]
      REAL RHOW     !Density of water (1.0*10^3) [kg m-3]
      REAL RHOICE   !Density of ice (0.917*10^3) [kg m-3]
      REAL TCGLAC   !Thermal conductivity of ice sheets (2.24) [W m-1 K-1]
      REAL CLHMLT   !Latent heat of freezing of water (0.334*10^6) [J kg-1]
      REAL CLHVAP   !Latent heat of vaporization of water (2.501*10^6) [J kg-1]
C
      COMMON /CLASS1/ DELT,TFREZ
      COMMON /CLASS4/ HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,
     1                SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,
     2                TCGLAC,CLHMLT,CLHVAP
C-----------------------------------------------------------------------
C
C   The adjustments of soil layer temperature and water content in this 
C   routine are done over the whole soil profile in the case of multiple 
C   soil layers (see the section on assignment of background data), but 
C   only to the bottom of the permeable depth in the case of the 
C   standard three-layer configuration (0.10, 0.25 and 3.75 m). This is 
C   because if the permeable depth lies within the thick third soil 
C   layer, it is recognized as desirable to apply the temperature 
C   changes only to that upper part of the layer in which the phase 
C   change is occurring, in order to avoid systematic damping of the 
C   temperature response of the layer. Thus the local array DELZZ (set 
C   in subroutine WPREP) is used here instead of DELZ when referring to 
C   the total thermal thickness of the soil layer, where DELZZ=DELZW for 
C   the third soil layer when the three-layer configuration is being 
C   used, but DELZZ=DELZ for all other cases.
C
      DO 100 J=1,IG                                                             
      DO 100 I=IL1,IL2
          IF(FI(I).GT.0. .AND. DELZW(I,J).GT.0. .AND. ISAND(I,1).GT.-4)
     1                                                              THEN
          !The heat capacity HCP of the permeable part DELZW of the soil 
          !layer under consideration is calculated here and in various 
          !other places in the subroutine as the weighted average of the 
          !heat capacities of the liquid water content THLIQ, the frozen 
          !water content THICE, and the soil material (taken to apply to 
          !the volume fraction not occupied by the pore volume THPOR). 
          !The heat capacity of air in the pores is neglected:
          !
          !HCP = Cw*THLIQ + Ci*THICE + Cm(1 - THPOR)
          !
          !Over the impermeable portion of the layer, the heat capacity 
          !of rock Cr is assumed to apply. Thus an effective heat 
          !capacity Cg,e (in units of J m-2 K-1) over the soil layer in 
          !question, DELZZ, can be calculated as:
          !
          !Cg,e = HCP*DELZW + Cr(DELZZ - DELZW)
          !
          !The change of internal energy I in the soil layers as a 
          !result of freezing and thawing is calculated as the 
          !difference in I between the beginning and end of the 
          !subroutine:
          !
          !delta_Ij = Xi*delta(Cg,e*TBAR)/DELT
          !
          !where TBAR is the temperature of the layer, DELT the 
          !length of the time step, and Xi the fractional coverage of
          !the subarea under consideration relative to the modelled 
          !area.
          !
              HCP  (I,J)=HCPW*THLIQ(I,J)+HCPICE*THICE(I,J)+
     1                   HCPS(I,J)*(1.-THPOR(I,J))                                   
              HTC  (I,J)=HTC(I,J)-FI(I)*(HCP(I,J)*DELZW(I,J)+
     1                   HCPSND*(DELZZ(I,J)-DELZW(I,J)))*
     2                   (TBAR(I,J)+TFREZ)/DELT
              !
              !If the soil layer temperature is less than 0 C and the 
              !volumetric liquid water content THLIQ of the layer is 
              !greater than the residual water content THLMIN, the water 
              !content THFREZ that can be frozen by the available energy 
              !sink is calculated from Ce and TBAR. The volumetric water 
              !content THEVAP of the first layer that is required to 
              !satisfy the surface evaporative flux is determined. For 
              !each layer, if THLIQ is found to exceed THLMIN + THEVAP, 
              !THFREZ is compared to the available water. If THFREZ ≤ 
              !THLIQ – THLMIN – THEVAP, all of the available energy sink 
              !is used to freeze part of the liquid water content in 
              !the permeable part of the soil layer, the amount of 
              !energy involved is subtracted from HTC and added to HMFG, 
              !HCP is recalculated and the layer temperature is set to 
              !0 C. Otherwise, all of the liquid water content of the 
              !layer above THLMIN + THEVAP is converted to frozen water, 
              !and HMFG and HTC are recalculated to reflect this. Then 
              !HCP is recomputed, and the remaining energy sink is 
              !applied to decreasing the temperature of the soil layer 
              !(both the permeable and impermeable portions) using Ce.
              !
              IF(TBAR(I,J).LT.0. .AND. THLIQ(I,J).GT.THLMIN(I,J)) THEN                        
                  THFREZ=-(HCP(I,J)*DELZW(I,J)+HCPSND*(DELZZ(I,J)-
     1                    DELZW(I,J)))*TBAR(I,J)/(CLHMLT*RHOW*
     2                    DELZW(I,J))                              
                  IF(J.EQ.1)                               THEN 
                      THEVAP=EVAP(I)*DELT/DELZW(I,J)                                          
                  ELSE                                                                
                      THEVAP=0.0                                                      
                  ENDIF                                                               
                  IF((THLIQ(I,J)-THLMIN(I,J)-THEVAP).GT.0.0)      THEN 
                    IF(THFREZ.LE.(THLIQ(I,J)-THLMIN(I,J)-THEVAP)) THEN                         
                      HMFG(I,J)=HMFG(I,J)-FI(I)*THFREZ*CLHMLT*
     1                          RHOW*DELZW(I,J)/DELT
                      HTC(I,J)=HTC(I,J)-FI(I)*THFREZ*CLHMLT*
     1                          RHOW*DELZW(I,J)/DELT
                      THLIQ(I,J)=THLIQ(I,J)-THFREZ                                        
                      THICE(I,J)=THICE(I,J)+THFREZ*RHOW/RHOICE                            
                      HCP  (I,J)=HCPW*THLIQ(I,J)+HCPICE*THICE(I,J)+
     1                           HCPS(I,J)*(1.-THPOR(I,J))
                      TBAR (I,J)=0.0                                                   
                    ELSE                                                                
                      HMFG(I,J)=HMFG(I,J)-FI(I)*(THLIQ(I,J)-
     1                    THLMIN(I,J)-THEVAP)*CLHMLT*RHOW*DELZW(I,J)/
     2                    DELT
                      HTC(I,J)=HTC(I,J)-FI(I)*(THLIQ(I,J)-THLMIN(I,J)-
     1                          THEVAP)*CLHMLT*RHOW*DELZW(I,J)/DELT
                      HADD=(THFREZ-(THLIQ(I,J)-THLMIN(I,J)-THEVAP))*
     1                     CLHMLT*RHOW*DELZW(I,J)
                      THICE(I,J)=THICE(I,J)+(THLIQ(I,J)-THLMIN(I,J)-
     1                           THEVAP)*RHOW/RHOICE          
                      THLIQ(I,J)=THLMIN(I,J)+THEVAP                                          
                      HCP  (I,J)=HCPW*THLIQ(I,J)+HCPICE*THICE(I,J)+
     1                           HCPS(I,J)*(1.-THPOR(I,J)) 
                      TBAR (I,J)=-HADD/(HCP(I,J)*DELZW(I,J)+HCPSND*
     1                           (DELZZ(I,J)-DELZW(I,J)))
                    ENDIF                                                               
                  ENDIF                                                               
              ENDIF
C               
              !
              !If the soil layer temperature is greater than 0 C and the 
              !volumetric ice content THICE of the layer is greater than 
              !zero, the ice content THMELT that can be melted by the 
              !available energy is calculated from Ce and TBAR. For 
              !each layer, if THMELT <= THICE, all of the available 
              !energy is used to melt part of the frozen water content 
              !of the permeable part of the layer, the amount of energy 
              !involved is subtracted from HTC and added to HMFG, HCP is 
              !recalculated and the layer temperature is set to 0 C. 
              !Otherwise, all of the frozen water content of the layer 
              !is converted to liquid water, and HMFG and HTC are 
              !recalculated to reflect this. Then HCP is recomputed, 
              !and the remaining energy is applied to increasing the 
              !temperature of the soil layer (both the permeable and 
              !impermeable portions) using Ce.
              !                                                    
              IF(TBAR(I,J).GT.0. .AND. THICE(I,J).GT.0.)        THEN                           
                  THMELT=(HCP(I,J)*DELZW(I,J)+HCPSND*(DELZZ(I,J)-
     1                   DELZW(I,J)))*TBAR(I,J)/(CLHMLT*RHOICE*
     2                   DELZW(I,J))                             
                  IF(THMELT.LE.THICE(I,J))                 THEN 
                      HMFG(I,J)=HMFG(I,J)+FI(I)*THMELT*CLHMLT*
     1                          RHOICE*DELZW(I,J)/DELT
                      HTC(I,J)=HTC(I,J)+FI(I)*THMELT*CLHMLT*
     1                          RHOICE*DELZW(I,J)/DELT
                      THICE(I,J)=THICE(I,J)-THMELT                                        
                      THLIQ(I,J)=THLIQ(I,J)+THMELT*RHOICE/RHOW                            
                      HCP  (I,J)=HCPW*THLIQ(I,J)+HCPICE*THICE(I,J)+
     1                           HCPS(I,J)*(1.-THPOR(I,J)) 
                      TBAR (I,J)=0.0                                                   
                  ELSE                                                                
                      HMFG(I,J)=HMFG(I,J)+FI(I)*THICE(I,J)*CLHMLT*
     1                          RHOICE*DELZW(I,J)/DELT
                      HTC(I,J)=HTC(I,J)+FI(I)*THICE(I,J)*CLHMLT*
     1                          RHOICE*DELZW(I,J)/DELT
                      HADD=(THMELT-THICE(I,J))*CLHMLT*RHOICE*
     1                     DELZW(I,J)
                      THLIQ(I,J)=THLIQ(I,J)+THICE(I,J)*RHOICE/RHOW                          
                      THICE(I,J)=0.0                                                    
                      HCP  (I,J)=HCPW*THLIQ(I,J)+HCPICE*THICE(I,J)+
     1                           HCPS(I,J)*(1.-THPOR(I,J))
                      TBAR (I,J)=HADD/(HCP(I,J)*DELZW(I,J)+HCPSND*
     1                           (DELZZ(I,J)-DELZW(I,J)))
                  ENDIF                                                               
              ENDIF
              !
              !In the final cleanup, the internal energy calculations 
              !for this subroutine are completed, and the first half of
              !a new set of internal energy calculations is done to span 
              !the subroutines treating ground water movement, which 
              !will be completed in subroutine TMCALC. Lastly, TBARW, 
              !the liquid water temperature of each soil layer, is 
              !assigned using TBAR.
              !
              HTC  (I,J)=HTC(I,J)+FI(I)*(HCP(I,J)*DELZW(I,J)+
     1                   HCPSND*(DELZZ(I,J)-DELZW(I,J)))*
     2                   (TBAR(I,J)+TFREZ)/DELT
              HTC(I,J)=HTC(I,J)-FI(I)*(TBAR(I,J)+TFREZ)*
     1                 (HCPW*THLIQ(I,J)+HCPICE*THICE(I,J))*
     2                 DELZW(I,J)/DELT
          ENDIF                                                      
          TBARW(I,J)=TBAR(I,J)
  100 CONTINUE                                                                    
C                                                                                  
      RETURN                                                                      
      END 
