      SUBROUTINE CLASSI(VPD,TADP,PADRY,RHOAIR,RHOSNI,
     1                  RPCP,TRPCP,SPCP,TSPCP,
     2                  TA,QA,PCPR,RRATE,SRATE,PRESSG,
     3                  IPCP,NL,IL1,IL2)
C
C     Purpose: Evaluate atmospheric variables and rainfall/snowfall 
C     rates over modelled area.
C
C     * NOV 17/11 - M.LAZARE.   REMOVE CALCULATION OF PCPR
C     *                         FOR IPCP=4 (REDUNDANT SINCE
C     *                         PCPR MUST BE PASSED IN FOR
C     *                         IF CONDITION ON LINE 100).
C     * NOV 22/06 - P.BARTLETT. CALCULATE PCPR IF IPCP=4.
C     * JUN 06/06 - V.FORTIN.   ADD OPTION FOR PASSING IN
C     *                         RAINFALL AND SNOWFALL RATES
C     *                         CALCULATED BY ATMOSPHERIC MODEL.
C     * NOV 03/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND;
C     *                         MOVE CLOUD FRACTION CALCULATION
C     *                         BACK INTO DRIVER.
C     * SEP 04/03 - D.VERSEGHY. NEW LOWER LIMIT ON PRECIPITATION
C     *                         RATE.
C     * AUG 09/02 - D.VERSEGHY. MOVE CALCULATION OF SOME
C     *                         ATMOSPHERIC VARIABLES HERE
C     *                         PRIOR TO GATHERING.
C     * JUL 26/02 - R.BROWN/S.FASSNACHT/D.VERSEGHY. PROVIDE 
C     *                         ALTERNATE METHODS OF ESTIMATING 
C     *                         RAINFALL/SNOWFALL PARTITIONING.
C     * JUN 27/02 - D.VERSEGHY. ESTIMATE FRACTIONAL CLOUD COVER
C     *                         AND RAINFALL/SNOWFALL RATES
C     *                         IF NECESSARY.
C
      IMPLICIT NONE
C
C     * INTEGER CONSTANTS.
C
      INTEGER IPCP,NL,IL1,IL2,I
C
C     * OUTPUT ARRAYS.
C
      REAL VPD   (NL)   !Vapour pressure deficit of air [mb] (ed)  
      REAL TADP  (NL)   !Dew point temperature of air [K]
      REAL PADRY (NL)   !Partial pressure of dry air [Pa] (pdry)
      REAL RHOAIR(NL)   !Density of air [kg m-3] (rho_a)
      REAL RHOSNI(NL)   !Density of fresh snow [kg m-3] (rho_s,i)
      REAL RPCP  (NL)   !Calculated rainfall rate over modelled area 
                        ![m s-1]
      REAL TRPCP (NL)   !Rainfall temperature over modelled area [C]
      REAL SPCP  (NL)   !Calculated snowfall rate over modelled area 
                        ![m s-1]
      REAL TSPCP (NL)   !Snowfall temperature over modelled area [C]
C
C     * INPUT ARRAYS.
C
      REAL TA    (NL)   !Air temperature at reference height [K] (Ta)  
      REAL QA    (NL)   !Specific humidity at reference height [kg kg-1] 
                        !(qa)
      REAL PRESSG(NL)   !Surface atmospheric pressure [Pa] (p)
      REAL PCPR  (NL)   !Precipitation rate over modelled area 
                        ![kg m-2 s-1]
      REAL RRATE (NL)   !Input rainfall rate over modelled area 
                        ![kg m-2 s-1]
      REAL SRATE (NL)   !Input snowfall rate over modelled area 
                        ![kg m-2 s-1]
C
C     * WORK ARRAYS.
C
      REAL PHASE (NL)
C
C     * TEMPORARY VARIABLES.
C
      REAL EA,CA,CB,EASAT,CONST
C
C     * COMMON BLOCK PARAMETERS.
C
      REAL DELT     !Time step [s]
      REAL TFREZ    !Freezing point of water [K]
      REAL RGAS     !Gas constant [J kg-1 K-1]
      REAL RGASV    !Gas constant for water vapour [J kg-1 K-1]
      REAL GRAV     !Acceleration due to gravity [m s-1]
      REAL SBC      !Stefan-Boltzmann constant [W m-2 K-4]
      REAL VKC      !Von Karman constant (0.40)
      REAL CT       !Drag coefficient for water (1.15*10^-3)
      REAL VMIN     !Minimum wind speed (0.1) [m s-1]
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
 
      COMMON /CLASS1/ DELT,TFREZ
      COMMON /CLASS2/ RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN
      COMMON /CLASS4/ HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,
     1                SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,
     2                TCGLAC,CLHMLT,CLHVAP
C----------------------------------------------------------------
C
C     * CALCULATION OF ATMOSPHERIC INPUT VARIABLES.
C
      !
      !In the first section, the air vapour pressure deficit, dry air 
      !pressure, air density and dew point temperature are calculated. 
      !The vapour pressure deficit VPD (in units of mb) is obtained from 
      !the saturated and actual vapour pressures of the air, ea and 
      !ea,sat respectively (in units of Pa), as
      !VPD = [ea,sat - ea] /100.0
      !
      !The air vapour pressure is obtained from the specific humidity QA 
      !using the formula
      !
      !ea = QA*PRESSG/[0.622 + 0.378*QA]
      !
      !where PRESSG is the surface atmospheric pressure. For the 
      !saturated vapour pressure, a standard empirical equation is 
      !utilized relating ea,sat to the temperature TA and the freezing 
      !point Tf: 
      !
      !ea,sat = 611.0*exp[17.269*(TA – Tf)/(TA – 35.86)]    TA >= Tf
      !ea,sat = 611.0*exp[21.874*(TA – Tf)/(TA – 7.66)]     TA < Tf
      !
      !The partial pressure of dry air, pdry, is obtained by subtracting 
      !ea from p, and the density of the air is calculated as the sum of 
      !the densities of the dry air and the water vapour:
      !
      !RHOAIR = pdry/RGAS*TA + ea/RGASV*TA
      !
      !where RGAS and RGASV are the gas constants for dry air and water 
      !vapour respectively. The dew point temperature of the air is 
      !evaluated by substituting ea for ea,sat on the left-hand side of 
      !the appropriate equation above, and solving for TA.
      !
      !Finally, if IPCP = 4, this indicates that the partitioning 
      !between rainfall and snowfall has been done outside of CLASS. 
      !The rainfall and snowfall rates RRATE and SRATE that have been 
      !passed into the subroutine are therefore simply assigned to RPCP 
      !and SPCP.
      !
      DO 100 I=IL1,IL2
          EA=QA(I)*PRESSG(I)/(0.622+0.378*QA(I))                              
          IF(TA(I).GE.TFREZ) THEN                                             
              CA=17.269                                                       
              CB=35.86                                                        
          ELSE                                                                
              CA=21.874                                                       
              CB=7.66                                                         
          ENDIF                                                               
          EASAT=611.0*EXP(CA*(TA(I)-TFREZ)/(TA(I)-CB))                        
          VPD(I)=MAX(0.0,(EASAT-EA)/100.0)                                     
          PADRY(I)=PRESSG(I)-EA                                                  
          RHOAIR(I)=PADRY(I)/(RGAS*TA(I))+EA/(RGASV*TA(I))                          
          CONST=LOG(EA/611.0)                                          
          TADP(I)=(CB*CONST-CA*TFREZ)/(CONST-CA)

C
C     * DENSITY OF FRESH SNOW.
C
          !
          !In the next section, the density of fresh snow RHOSNI is 
          !determined as an empirical function of the air temperature. 
          !For temperatures below 0 C, an equation presented by Hedstrom 
          !and Pomeroy (1998) is used. For temperatures >= 0 C, a 
          !relation following Pomeroy and Gray (1995) is used, with an 
          !upper limit of 200 kg m-3:
          !
          !RHOSNI = 67.92 + 51.25*exp[(TA – Tf)/2.59] TA < Tf
          !RHOSNI = 119.17 + 20.0*(TA – Tf)        TA >= Tf
          !
          IF(TA(I).LE.TFREZ) THEN
              RHOSNI(I)=67.92+51.25*EXP((TA(I)-TFREZ)/2.59)
          ELSE
              RHOSNI(I)=MIN((119.17+20.0*(TA(I)-TFREZ)),200.0)
          ENDIF

C
C     * PRECIPITATION PARTITIONING BETWEEN RAIN AND SNOW.
C
          !
          !In the last section, the partitioning of precipitation 
          !between rainfall RPCP and snowfall SPCP is addressed. Four 
          !options for doing so are provided; the user’s selection is 
          !indicated by the flag IPCP. In each case the rainfall and 
          !snowfall rates are converted to units of m s-1, by dividing 
          !by the density of water in the case of rain and by RHOSNI in 
          !the case of snow. The rainfall temperature is set to the 
          !maximum of 0 C and TA, and the snowfall temperature to the 
          !minimum of 0 C and TA. If IPCP = 1, the precipitation is 
          !simply diagnosed as rain if the air temperature is greater 
          !than 0 C, and as snow otherwise. If IPCP = 2, an empirical 
          !relation developed by Brown (2001) is used, where the 
          !precipitation is entirely snowfall when TA <= 0 C, and 
          !entirely rainfall when TA >= 2.0 C, and varies linearly 
          !between the two, with an equal mix of rain and snow at 
          !TA = 1.0 C. If IPCP = 3, the precipitation is assumed to be 
          !entirely snowfall when TA <= 0 C, and entirely rainfall when 
          !Ta ≥ 6.0 C, and between the two a polynomial function 
          !presented by Auer (1974) is used, relating the fraction of 
          !the precipitation that is snowfall, Xsf, to Ta:
          !
          !Xsf = [0.0202*TA^6 – 0.3660*TA^5 + 2.0399*TA^4 – 1.5089*TA^3
          !      – 15.038*TA^2 + 4.6664*TA + 100.0]/100.0
          !
          RPCP (I)=0.0  
          TRPCP(I)=0.0 
          SPCP (I)=0.0
          TSPCP(I)=0.0   
          IF(PCPR(I).GT.1.0E-8)                              THEN 
              IF(IPCP.EQ.1)                           THEN
                  IF(TA(I).GT.TFREZ) THEN
                      RPCP (I)=PCPR(I)/RHOW                   
                      TRPCP(I)=MAX((TA(I)-TFREZ),0.0)                   
                  ELSE
                      SPCP (I)=PCPR(I)/RHOSNI(I)
                      TSPCP(I)=MIN((TA(I)-TFREZ),0.0)               
                  ENDIF
              ELSEIF(IPCP.EQ.2)                       THEN
                  IF(TA(I).LE.TFREZ) THEN
                      PHASE(I)=1.0
                  ELSEIF(TA(I).GE.(TFREZ+2.0)) THEN
                      PHASE(I)=0.0 
                  ELSE
                      PHASE(I)=1.0-0.5*(TA(I)-TFREZ)
                  ENDIF
                  RPCP(I)=(1.0-PHASE(I))*PCPR(I)/RHOW
                  IF(RPCP(I).GT.0.0) TRPCP(I)=MAX((TA(I)-TFREZ),0.0) 
                  SPCP(I)=PHASE(I)*PCPR(I)/RHOSNI(I)
                  IF(SPCP(I).GT.0.0) TSPCP(I)=MIN((TA(I)-TFREZ),0.0)
              ELSEIF(IPCP.EQ.3)                       THEN
                  IF(TA(I).LE.TFREZ) THEN
                      PHASE(I)=1.0
                  ELSEIF(TA(I).GE.(TFREZ+6.0)) THEN
                      PHASE(I)=0.0
                  ELSE
                      PHASE(I)=(0.0202*(TA(I)-TFREZ)**6-0.3660*
     1                    (TA(I)-TFREZ)**5+2.0399*(TA(I)-TFREZ)**4-
     2                    1.5089*(TA(I)-TFREZ)**3-15.038*
     3                    (TA(I)-TFREZ)**2+4.6664*(TA(I)-TFREZ)+100.0)/
     4                    100.0
                      PHASE(I)=MAX(0.0,MIN(1.0,PHASE(I)))
                  ENDIF
                  RPCP(I)=(1.0-PHASE(I))*PCPR(I)/RHOW
                  IF(RPCP(I).GT.0.0) TRPCP(I)=MAX((TA(I)-TFREZ),0.0) 
                  SPCP(I)=PHASE(I)*PCPR(I)/RHOSNI(I)
                  IF(SPCP(I).GT.0.0) TSPCP(I)=MIN((TA(I)-TFREZ),0.0)
              ELSEIF(IPCP.EQ.4)                       THEN
                  RPCP(I)=RRATE(I)/RHOW
                  IF(RPCP(I).GT.0.0) TRPCP(I)=MAX((TA(I)-TFREZ),0.0) 
                  SPCP(I)=SRATE(I)/RHOSNI(I)
                  IF(SPCP(I).GT.0.0) TSPCP(I)=MIN((TA(I)-TFREZ),0.0)
              ENDIF
          ENDIF
100   CONTINUE
C
      RETURN
      END
