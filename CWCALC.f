      SUBROUTINE CWCALC(TCAN,RAICAN,SNOCAN,FRAINC,FSNOWC,CHCAP,
     1                  HMFC,HTCC,FI,CMASS,ILG,IL1,IL2,JL) 
C      
C     Purpose: Check for freezing or thawing of liquid or frozen water 
C     on the vegetation canopy, and adjust canopy temperature and 
C     intercepted water stores accordingly.
C
C                                                           
C     * MAR 25/08 - D.VERSEGHY. UPDATE FRAINC AND FSNOWC.
C     * SEP 24/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
C     * JUN 20/02 - D.VERSEGHY. COSMETIC REARRANGEMENT OF
C     *                         SUBROUTINE CALL; SHORTENED
C     *                         CLASS4 COMMON BLOCK.
C     * JUN 20/97 - D.VERSEGHY. CLASS - VERSION 2.7.
C     *                         PASS IN NEW "CLASS4" COMMON BLOCK.
C     * JAN 02/96 - D.VERSEGHY. CLASS - VERSION 2.5.
C     *                         COMPLETION OF ENERGY BALANCE
C     *                         DIAGNOSTICS.
C     * JUL 30/93 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.2.
C     *                                  NEW DIAGNOSTIC FIELDS.
C     * MAR 17/92 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.1.
C     *                                  REVISED AND VECTORIZED CODE
C     *                                  FOR MODEL VERSION GCM7.
C     * AUG 13/91 - D.VERSEGHY. ADJUST CANOPY TEMPERATURE AND
C     *                         INTERCEPTED LIQUID/FROZEN 
C     *                         MOISTURE STORES FOR FREEZING/
C     *                         THAWING.
C
      IMPLICIT NONE
C                                       
C     * INTEGER CONSTANTS.
C
      INTEGER ILG,IL1,IL2,JL,I
C
C     * I/O ARRAYS.
C
      REAL TCAN  (ILG)  !Temperature of vegetation canopy [K] (Tc)
      REAL RAICAN(ILG)  !Intercepted liquid water stored on the canopy [kg m-2] 
      REAL SNOCAN(ILG)  !Intercepted frozen water stored on the canopy [kg m-2]
      REAL FRAINC(ILG)  !Fractional coverage of canopy by liquid water [ ]  
      REAL FSNOWC(ILG)  !Fractional coverage of canopy by frozen water [ ]  
      REAL CHCAP (ILG)  !Heat capacity of vegetation canopy [J m-2 K-1] (Cc)  
      REAL HMFC  (ILG)  !Energy associated with freezing or thawing of water in canopy interception stores [W m-2]  
      REAL HTCC  (ILG)  !Internal energy change of canopy due to changes in temperature and/or mass [W m-2] (Ic)
C
C     * INPUT ARRAYS.
C
      REAL FI    (ILG)  !Fractional coverage of subarea in question on modelled area [ ] (Xi)
      REAL CMASS (ILG)  !Mass of vegetation canopy [kg m-2]

C
C     * TEMPORARY VARIABLES.
C
      REAL HFREZ,HCONV,RCONV,HCOOL,HMELT,SCONV,HWARM
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
C---------------------------------------------------------------------
C
C     The change of internal energy HTCC of the vegetation canopy as a 
C     result of the phase change processes treated here is calculated as 
C     the difference in HTCC between the beginning and end of the subroutine:
C
C     delta_HTCC = FI*delta[CHCAP*TCAN]/delta_t
C
C     where Cc represents the volumetric heat capacity of the canopy, Tc 
C     its temperature, delta_t the length of the time step, and FI the 
C     fractional coverage of the subarea under consideration relative to 
C     the modelled area.
C
      DO 100 I=IL1,IL2
          IF(FI(I).GT.0.)                                        THEN
              HTCC  (I)=HTCC(I)-FI(I)*TCAN(I)*CHCAP(I)/DELT
              !
              !If there is liquid water stored on the canopy and the
              !canopy temperature is less than 0 C, the available
              !energy sink HFREZ is calculated from CHCAP and the
              !difference between TCAN and 0 C, and
              !compared with HCONV, calculated as the energy sink 
              !required to freeze all of the liquid water on the
              !canopy. If HFREZ <= HCONV, the amount of water that 
              !can be frozen is calculated using the latent heat
              !of melting. The fractional coverages of frozen and 
              !liquid water FSNOWC and FRAINC and their masses
              !SNOCAN and RAICAN are adjusted accordingly, TCAN is 
              !set to 0 C, and the amount of energy
              !involved is subtracted from the internal energy HTCC
              !and added to HMFC. Otherwise all of the
              !intercepted liquid water is converted to frozen 
              !water, and the energy available for cooling the canopy is
              !calculated as HCOOL = HFREZ – HCONV. This available 
              !energy is applied to decreasing the
              !temperature of the canopy, using the specific heat of
              !the canopy elements, and the amount of energy that
              !was involved in the phase change is subtracted from 
              !HTCC and added to HMFC.
              !                                                 
              IF(RAICAN(I).GT.0. .AND. TCAN(I).LT.TFREZ)      THEN                                    
                  HFREZ=CHCAP(I)*(TFREZ-TCAN(I))                                                
                  HCONV=RAICAN(I)*CLHMLT  
                  IF(HFREZ.LE.HCONV)                       THEN 
                     RCONV=HFREZ/CLHMLT                                                  
                     FSNOWC(I)=FSNOWC(I)+FRAINC(I)*RCONV/RAICAN(I)
                     FRAINC(I)=FRAINC(I)-FRAINC(I)*RCONV/RAICAN(I)
                     SNOCAN(I)=SNOCAN(I)+RCONV                                                 
                     RAICAN(I)=RAICAN(I)-RCONV                                                 
                     TCAN  (I)=TFREZ                                                          
                     HMFC  (I)=HMFC(I)-FI(I)*CLHMLT*RCONV/DELT
                     HTCC  (I)=HTCC(I)-FI(I)*CLHMLT*RCONV/DELT
                  ELSE                                                                    
                     HCOOL=HFREZ-HCONV                                                   
                     SNOCAN(I)=SNOCAN(I)+RAICAN(I)                                                
                     FSNOWC(I)=FSNOWC(I)+FRAINC(I)
                     FRAINC(I)=0.0
                     TCAN  (I)=-HCOOL/(SPHVEG*CMASS(I)+SPHICE*
     1                         SNOCAN(I))+TFREZ  
                     HMFC  (I)=HMFC(I)-FI(I)*CLHMLT*RAICAN(I)/DELT
                     HTCC  (I)=HTCC(I)-FI(I)*CLHMLT*RAICAN(I)/DELT
                     RAICAN(I)=0.0                                                          
                  ENDIF                                                                   
              ENDIF
              ! 
              !If there is frozen water stored on the canopy and the 
              !canopy temperature is greater than 0 C, the available
              !energy for melting, HMELT, is calculated from CHCAP and 
              !the difference between TCAN and 0 C, and
              !compared with HCONV, calculated as the energy required to
              !melt all of the frozen water on the canopy.
              !If HMELT <= HCONV, the amount of frozen water that can be 
              !melted is calculated using the latent heat
              !of melting. The fractional coverages of frozen and liquid
              !water FSNOWC and FRAINC and their masses
              !SNOCAN and RAICAN are adjusted accordingly, TCAN is set 
              !to 0 C, and the amount of energy
              !involved is subtracted from HTCC and added to HMFC. 
              !Otherwise, all of the intercepted frozen water is
              !converted to liquid water, and the energy available for 
              !warming the canopy is calculated as HWARM =
              !HMELT – HCONV. This available energy is applied to 
              !increasing the temperature of the canopy, using
              !the specific heats of the canopy elements, and the amount
              !of energy that was involved in the phase
              !change is subtracted from HTCC and added to HMFC.
              !                                                        
              IF(SNOCAN(I).GT.0. .AND. TCAN(I).GT.TFREZ)        THEN 
                  HMELT=CHCAP(I)*(TCAN(I)-TFREZ)                                                
                  HCONV=SNOCAN(I)*CLHMLT                                                     
                  IF(HMELT.LE.HCONV)                       THEN 
                     SCONV=HMELT/CLHMLT                                                  
                     FRAINC(I)=FRAINC(I)+FSNOWC(I)*SCONV/SNOCAN(I)
                     FSNOWC(I)=FSNOWC(I)-FSNOWC(I)*SCONV/SNOCAN(I)
                     SNOCAN(I)=SNOCAN(I)-SCONV                                                 
                     RAICAN(I)=RAICAN(I)+SCONV                                                 
                     TCAN(I)=TFREZ                                                          
                     HMFC  (I)=HMFC(I)+FI(I)*CLHMLT*SCONV/DELT
                     HTCC  (I)=HTCC(I)+FI(I)*CLHMLT*SCONV/DELT
                  ELSE                                                                    
                     HWARM=HMELT-HCONV                                                   
                     RAICAN(I)=RAICAN(I)+SNOCAN(I)                                                
                     FRAINC(I)=FRAINC(I)+FSNOWC(I)
                     FSNOWC(I)=0.0
                     TCAN(I)=HWARM/(SPHVEG*CMASS(I)+SPHW*RAICAN(I))+
     1                       TFREZ                         
                     HMFC  (I)=HMFC(I)+FI(I)*CLHMLT*SNOCAN(I)/DELT
                     HTCC  (I)=HTCC(I)+FI(I)*CLHMLT*SNOCAN(I)/DELT
                     SNOCAN(I)=0.0                                                          
                  ENDIF                                                                   
              ENDIF 
              !
              !In the final cleanup, the canopy heat capacity is 
              !recomputed and the remaining internal energy calculations
              !are completed.
              !
              CHCAP(I)=SPHVEG*CMASS(I)+SPHW*RAICAN(I)+SPHICE*SNOCAN(I)
              HTCC (I)=HTCC(I)+FI(I)*TCAN(I)*CHCAP(I)/DELT
          ENDIF                                
  100 CONTINUE
C                                                                                  
      RETURN                                                                      
      END 
