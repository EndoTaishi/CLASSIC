      SUBROUTINE SNOVAP(RHOSNO,ZSNOW,HCPSNO,TSNOW,EVAP,QFN,QFG,HTCS,
     1                  WLOST,TRUNOF,RUNOFF,TOVRFL,OVRFLW,
     2                  FI,R,S,RHOSNI,WSNOW,ILG,IL1,IL2,JL)
C
C     * AUG 25/11 - D.VERSEGHY. CORRECT CALCULATION OF TRUNOF
C     *                         AND TOVRFL.
C     * FEB 22/07 - D.VERSEGHY. NEW ACCURACY LIMITS FOR R AND S.
C     * MAR 24/06 - D.VERSEGHY. ALLOW FOR PRESENCE OF WATER IN SNOW.
C     * SEP 23/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
C     * JUL 26/02 - D.VERSEGHY. CHANGE RHOSNI FROM CONSTANT TO
C     *                         VARIABLE.
C     * APR 11/01 - M.LAZARE.   CHECK FOR EXISTENCE OF SNOW BEFORE
C     *                         PERFORMING CALCULATIONS.
C     * JUN 20/97 - D.VERSEGHY. CLASS - VERSION 2.7.
C     *                         PASS IN NEW "CLASS4" COMMON BLOCK.
C     * JAN 02/96 - D.VERSEGHY. CLASS - VERSION 2.5.
C     *                         COMPLETION OF ENERGY BALANCE
C     *                         DIAGNOSTICS.
C     * AUG 16/95 - D.VERSEGHY. CLASS - VERSION 2.4.
C     *                         INCORPORATE DIAGNOSTIC ARRAY "WLOST". 
C     * DEC 22/94 - D.VERSEGHY. CLASS - VERSION 2.3.
C     *                         ADDITIONAL DIAGNOSTIC CALCULATION -
C     *                         UPDATE HTCS.
C     * JUL 30/93 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.2.
C     *                                  NEW DIAGNOSTIC FIELDS.
C     * APR 24/92 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.1.
C     *                                  REVISED AND VECTORIZED CODE
C     *                                  FOR MODEL VERSION GCM7.
C     * AUG 12/91 - D.VERSEGHY. CODE FOR MODEL VERSION GCM7U -
C     *                         CLASS VERSION 2.0 (WITH CANOPY).
C     * APR 11/89 - D.VERSEGHY. SUBLIMATION FROM SNOWPACK.
C                                          
      IMPLICIT NONE
C
C     * INTEGER CONSTANTS.
C
      INTEGER ILG,IL1,IL2,JL,I
C
C     * INPUT/OUTPUT ARRAYS.
C
      REAL RHOSNO(ILG)  !Density of snow pack [kg m-3] (rho_s )   
      REAL ZSNOW (ILG)  !Depth of snow pack [m] (zg) 
      REAL HCPSNO(ILG)  !Heat capacity of snow pack [J m-3 K-1] (Cs) 
      REAL TSNOW (ILG)  !Temperature of the snow pack [C] (Ts)
      REAL EVAP  (ILG)  !Sublimation rate from snow surface at start of 
                        !subroutine [m s-1]
      REAL QFN   (ILG)  !Sublimation from snow pack [kg m-2 s-1] 
      REAL QFG   (ILG)  !Evaporation from ground [kg m-2 s-1] 
      REAL HTCS  (ILG)  !Internal energy change of snow pack due to 
                        !conduction and/or change in mass [W m-2] (Is)
      REAL WLOST (ILG)  !Residual amount of water that cannot be 
                        !supplied by surface stores [kg m-2]
      REAL TRUNOF(ILG)  !Temperature of total runoff [K] 
      REAL RUNOFF(ILG)  !Total runoff [m s-1] 
      REAL TOVRFL(ILG)  !Temperature of overland flow [K]
      REAL OVRFLW(ILG)  !Overland flow from top of soil column [m s-1]
C
C     * INPUT ARRAYS.
C
      REAL FI    (ILG)  !Fractional coverage of subarea in question on 
                        !modelled area [ ] (Xi)
      REAL R     (ILG)  !Rainfall rate incident on snow pack [m -1s ] 
      REAL S     (ILG)  !Snowfall rate incident on snow pack 
                        ![kg m-2 s-1]
      REAL RHOSNI(ILG)  !Density of fresh snow [kg m-3]
      REAL WSNOW (ILG)  !Liquid water content of snow pack [kg m-2] (ws) 
C
C     * TEMPORARY VARIABLES.
C
      REAL ZADD,ZLOST,ZREM
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
      !
      !These calculations are done if a snowpack is present and there is 
      !no rainfall or snowfall occurring. The change of internal energy 
      !HTCS of the snow pack as a result of the change in its mass is 
      !calculated as the difference in Is between the beginning and end 
      !of the subroutine:
      !
      !HTCS = FI*delta[HCPSNO*ZSNOW*TSNOW]/DELT
      !
      !where HCPSNO represents the volumetric heat capacity of the snow 
      !pack, TSNOW its temperature, DELT the length of the time step, 
      !and FI the fractional coverage of the subarea under consideration 
      !relative to the modelled area.
      !
      !If the sublimation rate EVAP over the snow pack is negative 
      !(downward), the deposited depth of snow ZADD is calculated from 
      !EVAP by converting it from a liquid water flux to a fresh snow 
      !depth using RHOSNI, the fresh snow density. The snowpack density 
      !is updated as a weighted average of the original snow density 
      !RHOSNO and RHOSNI. The new snow depth is calculated as the sum of 
      !the old snow depth ZSNOW and ZADD. The new volumetric heat 
      !capacity of the snow pack is obtained from the heat capacities of 
      !ice and water HCPICE and HCPW, the snow, ice and water densities 
      !RHOSNO RHOICE, and RHOW, and the water content and depth of the 
      !snow pack WSNOW and ZSNOW, as:
      !
      !HCPSNO = HCPICE*[RHOSNO/RHOICE] + HCPW*WSNOW/[RHOW*ZSNOW]
      !
      !If the sublimation rate is positive, the depth of the snow pack 
      !ZLOST that is sublimated over the time step is calculated from 
      !EVAP using RHOSNO. If ZLOST <= ZSNOW, the snow depth is reduced 
      !and HCPSNO is recalculated. Otherwise the deficit amount ZREM is 
      !calculated from ZLOST – ZSNOW and converted to a depth of water. 
      !This amount is further converted to an evaporation rate by 
      !applying a correction factor of (Lm+ Lv)/Lv, where Lm is the 
      !latent heat of melting and Lv is the latent heat of vaporization 
      !(to account for the fact that the energy is now being used to 
      !evaporate water instead of sublimate snow). This necessarily 
      !leads to a small discrepancy between the overall vapour flux for 
      !the subarea that was originally calculated in CLASST, and the 
      !actual change of water storage in the subarea, and therefore this 
      !discrepancy is added to the housekeeping variable WLOST for use 
      !in the water balance checks done later in CHKWAT. If there was 
      !liquid water in the snow pack, WSNOW, it is assigned to overall 
      !runoff RUNOFF, and to overland flow OVRFLW. The resulting 
      !temperatures of the runoff and overland flow, TRUNOF and TOVRFL, 
      !are recalculated as weighted averages using the original runoff 
      !amounts and temperatures, and the original snow temperature TSNOW 
      !for WSNOW. The snow depth, heat capacity, temperature and water 
      !content are all set to zero. Finally, since ZREM now becomes soil 
      !evaporation rather than snow sublimation, the diagnostic 
      !variables QFN and QFG, representing the vapour flux from snow and 
      !soil respectively, are adjusted to reflect this.
      !
      DO 100 I=IL1,IL2
          IF(FI(I).GT.0. .AND. (S(I).LT.1.0E-11 .OR. R(I).LT.1.0E-11)
     1                .AND. ZSNOW(I).GT.0.)                       THEN
              HTCS(I)=HTCS(I)-FI(I)*HCPSNO(I)*(TSNOW(I)+TFREZ)*
     1                ZSNOW(I)/DELT
              IF(EVAP(I).LT.0.)                             THEN 
                  ZADD=-EVAP(I)*DELT*RHOW/RHOSNI(I)
                  RHOSNO(I)=(ZSNOW(I)*RHOSNO(I)+ZADD*RHOSNI(I))/
     1                      (ZSNOW(I)+ZADD)                          
                  ZSNOW (I)=ZSNOW(I)+ZADD                                                        
                  HCPSNO(I)=HCPICE*RHOSNO(I)/RHOICE+HCPW*WSNOW(I)/
     1                (RHOW*ZSNOW(I))
                  EVAP  (I)=0.0                                                                
              ELSE                                                                        
                  ZLOST=EVAP(I)*DELT*RHOW/RHOSNO(I)
                  IF(ZLOST.LE.ZSNOW(I))                     THEN 
                      ZSNOW(I)=ZSNOW(I)-ZLOST                                                   
                      EVAP (I)=0.0                                                            
                      HCPSNO(I)=HCPICE*RHOSNO(I)/RHOICE+HCPW*WSNOW(I)/
     1                    (RHOW*ZSNOW(I))
                  ELSE                                                                    
                      ZREM=(ZLOST-ZSNOW(I))*RHOSNO(I)/RHOW
                      ZSNOW(I)=0.0                                                           
                      HCPSNO(I)=0.0
                      EVAP(I)=ZREM*(CLHMLT+CLHVAP)/(CLHVAP*DELT)
                      WLOST(I)=WLOST(I)-ZREM*RHOW*CLHMLT/CLHVAP
                      IF(RUNOFF(I).GT.0. .OR. WSNOW(I).GT.0.)
     1                 TRUNOF(I)=(TRUNOF(I)*RUNOFF(I)+(TSNOW(I)+TFREZ)*
     1                      WSNOW(I)/RHOW)/(RUNOFF(I)+WSNOW(I)/RHOW)
                      RUNOFF(I)=RUNOFF(I)+WSNOW(I)/RHOW
                      IF(OVRFLW(I).GT.0. .OR. WSNOW(I).GT.0.)
     1                 TOVRFL(I)=(TOVRFL(I)*OVRFLW(I)+(TSNOW(I)+TFREZ)*
     1                      FI(I)*WSNOW(I)/RHOW)/(OVRFLW(I)+FI(I)*
     2                      WSNOW(I)/RHOW)
                      OVRFLW(I)=OVRFLW(I)+FI(I)*WSNOW(I)/RHOW
                      TSNOW(I)=0.0 
                      WSNOW(I)=0.0
                      QFN(I)=QFN(I)-FI(I)*ZREM*RHOW/DELT
                      QFG(I)=QFG(I)+FI(I)*EVAP(I)*RHOW
                  ENDIF                                                                   
              ENDIF 
              HTCS(I)=HTCS(I)+FI(I)*HCPSNO(I)*(TSNOW(I)+TFREZ)*
     1                ZSNOW(I)/DELT
          ENDIF                                                                      
  100 CONTINUE
C                                                                                  
      RETURN                                                                      
      END        
