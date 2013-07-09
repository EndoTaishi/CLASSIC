      SUBROUTINE TNPOST(TBARPR,G12,G23,TPOND,GZERO,QFREZG,GCONST,
     1                  GCOEFF,TBAR,TCTOP,TCBOT,HCP,ZPOND,TSURF,
     2                  TBASE,TBAR1P,A1,A2,B1,B2,C2,FI,IWATER,
     3                  ISAND,DELZ,DELZW,ILG,IL1,IL2,JL,IG       )
C
C     Purpose: Soil heat flux calculations and cleanup after surface 
C     energy budget calculations.
C
C     * NOV 01/06 - D.VERSEGHY. ALLOW PONDING ON ICE SHEETS.
C     * OCT 04/05 - D.VERSEGHY. MODIFY 300 LOOP FOR CASES WHERE IG>3.
C     * NOV 04/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
C     * JUN 17/02 - D.VERSEGHY. RESET PONDED WATER TEMPERATURE
C     *                         USING CALCULATED GROUND HEAT FLUX;
C     *                         SHORTENED CLASS4 COMMON BLOCK.
C     * JUN 20/97 - D.VERSEGHY. CLASS - VERSION 2.7.
C     *                         INCORPORATE EXPLICITLY CALCULATED
C     *                         THERMAL CONDUCTIVITIES AT TOPS AND
C     *                         BOTTOMS OF SOIL LAYERS, AND 
C     *                         MODIFICATIONS TO ALLOW FOR VARIABLE
C     *                         SOIL PERMEABLE DEPTH.
C     * SEP 27/96 - D.VERSEGHY. CLASS - VERSION 2.5.
C     *                         FIX BUG IN CALCULATION OF FLUXES
C     *                         BETWEEN SOIL LAYERS (PRESENT SINCE 
C     *                         FIRST RELEASE OF VERSION 2.5).
C     * DEC 22/94 - D.VERSEGHY. CLASS - VERSION 2.3.
C     *                         REVISE CALCULATION OF TBARPR(I,1).
C     * APR 10/92 - M.LAZARE.   CLASS - VERSION 2.2.
C     *                         DIVIDE PREVIOUS SUBROUTINE "T3LAYR" INTO
C     *                         "TNPREP" AND "TNPOST" AND VECTORIZE.
C     * APR 11/89 - D.VERSEGHY. CALCULATE HEAT FLUXES BETWEEN SOIL 
C     *                         LAYERS; DISAGGREGATE FIRST SOIL LAYER
C     *                         TEMPERATURE INTO PONDED WATER AND
C     *                         SOIL TEMPERATURES; CONSISTENCY CHECK 
C     *                         ON CALCULATED SURFACE LATENT HEAT OF
C     *                         MELTING/FREEZING; CONVERT SOIL LAYER
C     *                         TEMPERATURES TO DEGREES C.
C
      IMPLICIT NONE
C                                                                                 
C     * INTEGER CONSTANTS.
C
      INTEGER ILG,IL1,IL2,JL,IG,I,J
C
C     * OUTPUT ARRAYS.
C
      REAL TBARPR(ILG,IG)   !Temperatures of soil layers for subarea [C]
C
      REAL G12   (ILG)  !Heat conduction between first and second soil 
                        !layers [W m-2 ] (G(delta_z1p))
      REAL G23   (ILG)  !Heat conduction between second and third soil 
                        !layers [W m-2 ]
      REAL TPOND (ILG)  !Temperature of ponded water [C] (Tp)
C
C     * INPUT/OUTPUT ARRAYS.
C
      REAL GZERO (ILG)  !Heat conduction into soil surface [W m-2 ] 
                        !(G(0))
      REAL QFREZG(ILG)  !Energy sink to be applied to freezing of ponded 
                        !water [W m-2 ]
C
C     * INPUT ARRAYS.
C
      REAL TBAR  (ILG,IG)   !Temperatures of soil layers, averaged over 
                            !modelled area [K]
      REAL TCTOP (ILG,IG)   !Thermal conductivity of soil at top of 
                            !layer [W m-1 K-1] (lambda)
      REAL TCBOT (ILG,IG)   !Thermal conductivity of soil at bottom of 
                            !layer [W m-1 K-1] (lambda)
      REAL HCP   (ILG,IG)   !Heat capacity of soil layer [J m-3 K-1]
      REAL DELZW (ILG,IG)   !Permeable thickness of soil layer [m]
C
      REAL ZPOND (ILG)  !Depth of ponded water on surface [m] (delta_zp)    
      REAL TSURF (ILG)  !Ground surface temperature [K]  
      REAL TBASE (ILG)  !Temperature of bedrock in third soil layer [K]  
      REAL TBAR1P(ILG)  !Lumped temperature of ponded water and first 
                        !soil layer [K] (T1p)
      REAL A1    (ILG)  !Work array used in calculation of GCONST and 
                        !GCOEFF
      REAL A2    (ILG)  !Work array used in calculation of GCONST and 
                        !GCOEFF
      REAL B1    (ILG)  !Work array used in calculation of GCONST and 
                        !GCOEFF
      REAL B2    (ILG)  !Work array used in calculation of GCONST and 
                        !GCOEFF
      REAL C2    (ILG)  !Work array used in calculation of GCONST and 
                        !GCOEFF
      REAL FI    (ILG)  !Fractional coverage of subarea in question on 
                        !modelled area [ ]
      REAL GCONST(ILG)  !Intercept used in equation relating ground 
                        !surface heat flux to surface temperature 
                        ![W m-2 ]
      REAL GCOEFF(ILG)  !Multiplier used in equation relating ground 
                        !surface heat flux to surface temperature 
                        ![W m-2 K-1]
C
      INTEGER              IWATER(ILG)      !Flag indicating condition 
                                            !of surface (dry, water-
                                            !covered or snow-covered)
      INTEGER              ISAND (ILG,IG)   !Sand content flag
C
      REAL DELZ  (IG)   !Overall thickness of soil layer [m] (delta_z)
C
C     * TEMPORARY VARIABLES.
C
      REAL GZROLD,DELZ1
C
C     * COMMON BLOCK PARAMETERS.
C
      REAL DELT     !Time step [s]
      REAL TFREZ    !Freezing point of water [K]
      REAL HCPW     !Volumetric heat capacity of water (4.187*10^6) 
                    ![J m-3 K-1]
      REAL HCPICE   !Volumetric heat capacity of ice (1.9257*10^6) 
                    ![J m-3 K-1]
      REAL HCPSOL   !Volumetric heat capacity of mineral matter 
                    !(2.25*10^6) [J m-3 K-1]
      REAL HCPOM    !Volumetric heat capacity of organic matter 
                    !(2.50*10^6) [J m-3 K-1]
      REAL HCPSND   !Volumetric heat capacity of sand particles 
                    !(2.13*10^6) [J m-3 K-1]
      REAL HCPCLY   !Volumetric heat capacity of fine mineral particles 
                    !(2.38*10^6) [J m-3 K-1]
      REAL SPHW     !Specific heat of water (4.186*10^3) [J kg-1 K-1]
      REAL SPHICE   !Specific heat of ice (2.10*10^3) [J kg-1 K-1]
      REAL SPHVEG   !Specific heat of vegetation matter (2.70*10^3) 
                    ![J kg-1 K-1]
      REAL SPHAIR   !Specific heat of air [J kg-1 K-1]
      REAL RHOW     !Density of water (1.0*10^3) [kg m-3]
      REAL RHOICE   !Density of ice (0.917*10^3) [kg m-3]
      REAL TCGLAC   !Thermal conductivity of ice sheets (2.24) 
                    ![W m-1 K-1]
      REAL CLHMLT   !Latent heat of freezing of water (0.334*10^6) 
                    ![J kg-1]
      REAL CLHVAP   !Latent heat of vaporization of water (2.501*10^6) 
                    ![J kg-1]
C
      COMMON /CLASS1/ DELT,TFREZ                                                  
      COMMON /CLASS4/ HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,
     1                SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,
     2                TCGLAC,CLHMLT,CLHVAP
C-----------------------------------------------------------------------
C
      !In the 100 loop, the heat flux into the ground (without 
      !adjustments that may have been applied owing to phase changes of 
      !water at the surface or partitioning the residual of the surface 
      !energy balance among the surface flux terms) is calculated from 
      !the ground surface temperature TSURF, using the GCOEFF and GCONST 
      !terms (see documentation of subroutine TNPREP). This flux is then 
      !used to back-calculate the heat conduction between the first and 
      !second, and the second and third soil layers.
      !
      !If the depth of water ponded on the surface is greater than zero, 
      !the total thickness of the lumped layer consisting of the first 
      !soil layer and the ponded water, delta_z1p, is calculated as the 
      !sum of the two depths. The temperature of the ponded water over 
      !the subarea in question is disaggregated from the temperature of 
      !the lumped layer, T1p, by making use of the calculated heat 
      !fluxes at the top and bottom of the layer, and the assumption 
      !(discussed in the documentation of subroutine TNPREP) that that 
      !the variation of temperature T with depth z within a soil layer 
      !can be modelled using a quadratic equation:
      !
      !T(z) = 1⁄2 T′′(0)z2 + T′(0)z +T(0)
      !
      !Integrating this equation separately over the ponded water depth 
      !ZPOND and over the thickness of the lumped layer produces 
      !expressions for the ponded water temperature TPOND and the 
      !temperature of the lumped layer. Making use of the fact that
      !
      !T′′(0) = [T′(DELZ) - T′(0)]/DELZ
      !
      !where DELZ is a depth interval, and
      !
      !G(z) = - lambda(z)T′(z)
      !where lambda represents the thermal conductivity, an expression 
      !for TPOND can be derived:
      !
      !TPOND = [G(0)/lambda(0) – G(delta_z1p)/lambda(delta_z1p)]*
      ![ZPOND^2- delta_z1p^2]/6*delta_z1p – 
      !G(0)*[ZPOND- delta_z1p]/2*lambda(0) – TBAR1P
      !
      !The temperature TBARPR of the first soil layer over the subarea 
      !in question can then be obtained by disaggregating the ponded 
      !water temperature from the lumped layer temperature using the 
      !respective heat capacities of the ponded water and the soil 
      !layer. (The heat capacity of the soil is determined as the 
      !weighted average of HCP over the permeable thickness DELZW, and 
      !the heat capacity of rock, HCPSND, over the impermeable 
      !thickness, DELZ-DELZW.) Both the ponded water temperature and the 
      !soil layer temperature are converted to C. Lastly, if there is 
      !ponded water on the surface (IWATER = 1) and QFREZG, the surface 
      !energy available for phase change of water, is positive 
      !(indicating an energy source), or if there is snow on the ground 
      !(IWATER = 2) and QFREZG is negative (indicating an energy sink), 
      !or if there is no liquid or frozen water on the surface 
      !(IWATER = 0), QFREZG is added to the heat flux into the ground, 
      !GZERO, and then reset to zero.
      !
      DO 100 I=IL1,IL2
          IF(FI(I).GT.0.)                                          THEN
              GZROLD=GCOEFF(I)*TSURF(I)+GCONST(I)
              G12(I)=(TSURF(I)-TBAR1P(I)-A1(I)*GZROLD)/B1(I)
              G23(I)=(TSURF(I)-TBAR(I,2)-A2(I)*GZROLD-B2(I)*G12(I))/
     1               C2(I)                                    
              IF(ZPOND(I).GT.0.)                                THEN 
                  DELZ1=DELZ(1)+ZPOND(I)
                  TPOND(I)=(GZROLD/TCTOP(I,1)-G12(I)/TCBOT(I,1))*
     1                     (ZPOND(I)*ZPOND(I)-DELZ1*DELZ1)/(6.0*DELZ1)-
     2                     GZROLD*(ZPOND(I)-DELZ1)/(2.0*TCTOP(I,1))+
     3                     TBAR1P(I)-TFREZ
                  TBARPR(I,1)=((HCP(I,1)*DELZW(I,1)+HCPSND*(DELZ(1)-
     1                        DELZW(I,1))+HCPW*ZPOND(I))*TBAR1P(I)-
     2                        HCPW*ZPOND(I)*(TPOND(I)+TFREZ))/
     3                        (HCP(I,1)*DELZW(I,1)+HCPSND*(DELZ(1)-
     4                        DELZW(I,1)))-TFREZ
              ELSE                                                                        
                  TPOND(I)=0.                                                               
                  TBARPR(I,1)=TBAR(I,1)-TFREZ                                             
              ENDIF           
C
              IF((IWATER(I).EQ.1 .AND. QFREZG(I).GT.0.) .OR.
     1           (IWATER(I).EQ.2 .AND. QFREZG(I).LT.0.) .OR.
     2            IWATER(I).EQ.0)                               THEN              
                  GZERO(I)=GZERO(I)+QFREZG(I)                                                      
                  QFREZG(I)=0.                                                              
              ENDIF
          ENDIF
  100 CONTINUE
C 
      !In loop 200, the subarea soil layer temperatures TBARPR are set 
      !for the remaining soil layers. In all cases the temperature of 
      !the layer is set to that for the modelled area, TBAR, converted 
      !to C, except in the case of the third soil layer if the standard 
      !three-layer configuration is being modelled (with a very thick 
      !third soil layer of 3.75 m). In this case TBARPR and the layer 
      !heat capacity HCP are considered to apply to the permeable depth 
      !DELZW of the layer, and the bedrock temperature TBASE and the 
      !rock heat capacity HCPSND apply to the remainder, DELZ-DELZW. The 
      !disaggregation of TBARPR from TBAR and TBASE is carried out on 
      !this basis. TBARPR for this layer is also converted to C.
      !
      DO 200 I=IL1,IL2
          IF(FI(I).GT.0.)                                          THEN
              TBARPR(I,2)=TBAR(I,2)-TFREZ
              IF(DELZW(I,3).GT.0.0 .AND. DELZW(I,3).LT.DELZ(3)
     1                             .AND. IG.EQ.3)              THEN
                  TBARPR(I,3)=(TBAR(I,3)*(HCP(I,3)*DELZW(I,3)+
     1                         HCPSND*(DELZ(3)-DELZW(I,3)))-TBASE(I)*
     2                         HCPSND*(DELZ(3)-DELZW(I,3)))/(HCP(I,3)*
     3                         DELZW(I,3))-TFREZ
              ELSE
                  DO 150 J=3,IG
                      TBARPR(I,J)=TBAR(I,J)-TFREZ
  150             CONTINUE
              ENDIF
          ENDIF
  200 CONTINUE                                                                    
C
      RETURN
      END 
