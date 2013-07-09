      SUBROUTINE TSPREP(GCOEFFS,GCONSTS,CPHCHG,IWATER,
     1                  FI,ZSNOW,TSNOW,TCSNOW,
     2                  ILG,IL1,IL2,JL      )
C
C     Calculate coefficients for solution of snow pack heat conduction.
C
C     * AUG 16/06 - D.VERSEGHY. MAJOR REVISION TO IMPLEMENT THERMAL
C     *                         SEPARATION OF SNOW AND SOIL.
C     * MAY 24/06 - D.VERSEGHY. LIMIT DELZ3 TO <= 4.1 M.
C     * OCT 04/05 - D.VERSEGHY. USE THREE-LAYER TBAR,TCTOP,TCBOT.
C     * NOV 04/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
C     * AUG 06/02 - D.VERSEGHY. SHORTENED CLASS3 COMMON BLOCK.
C     * JUN 17/02 - D.VERSEGHY. USE NEW LUMPED SOIL AND PONDED WATER
C     *                         TEMPERATURE FOR FIRST LAYER; SHORTENED
C     *                         CLASS4 COMMON BLOCK.
C     * JUN 20/97 - D.VERSEGHY. CLASS - VERSION 2.7.
C     *                         INCORPORATE EXPLICITLY CALCULATED
C     *                         THERMAL CONDUCTIVITIES AT TOPS AND 
C     *                         BOTTOMS OF SOIL LAYERS.
C     * AUG 24/95 - D.VERSEGHY. CLASS - VERSION 2.4.
C     *                         REVISIONS TO ALLOW FOR INHOMOGENEITY
C     *                         BETWEEN SOIL LAYERS AND FRACTIONAL
C     *                         ORGANIC MATTER CONTENT.
C     * NOV 28/94 - M. LAZARE.  CLASS - VERSION 2.3.
C     *                         TCSATW,TCSATI DECLARED REAL(16).
C     * APR 10/92 - M. LAZARE.  CLASS - VERSION 2.1.
C     *                         DIVIDE PREVIOUS SUBROUTINE "T4LAYR" 
C     *                         INTO "TSPREP" AND "TSPOST" AND
C     *                         VECTORIZE.
C     * APR 11/89 - D.VERSEGHY. CALCULATE COEFFICIENTS FOR GROUND HEAT
C     *                         FLUX, EXPRESSED AS A LINEAR FUNCTION OF
C     *                         SURFACE TEMPERATURE. COEFFICIENTS ARE
C     *                         CALCULATED FROM LAYER TEMPERATURES, 
C     *                         THICKNESSES AND THERMAL CONDUCTIVITIES,
C     *                         ASSUMING A QUADRATIC VARIATION OF
C     *                         TEMPERATURE WITH DEPTH WITHIN EACH
C     *                         SOIL/SNOW LAYER. SET THE SURFACE 
C     *                         LATENT HEAT OF VAPORIZATION OF WATER
C     *                         AND THE STARTING TEMPERATURE FOR THE
C     *                         ITERATION IN "TSOLVC"/"TSOLVE".
C
      IMPLICIT NONE
C                                                                                 
C     * INTEGER CONSTANTS.
C
      INTEGER ILG,IL1,IL2,JL,I,J
C
C     * OUTPUT ARRAYS.
C
      REAL GCOEFFS(ILG)     !Multiplier used in equation relating snow surface heat flux to snow surface
                            !temperature [W m-2 K-1]
      REAL GCONSTS(ILG)     !Intercept used in equation relating snow surface heat flux to snow surface
                            !temperature [W m-2 ]
      REAL CPHCHG(ILG)      !Latent heat of sublimation [J kg-1]

C
      INTEGER              IWATER(ILG)     !Flag indicating condition of surface (dry, water-covered or snow-covered)
C
C     * INPUT ARRAYS.
C
      REAL FI    (ILG)  !Fractional coverage of subarea in question on modelled area [ ]
      REAL ZSNOW (ILG)  !Depth of snow pack [m] (delta_zs)
      REAL TSNOW (ILG)  !Snowpack temperature [K] (Ts(delta_zs))
      REAL TCSNOW(ILG)  !Thermal conductivity of snow [W m-1 K-1]
C
C     * COMMON BLOCK PARAMETERS.
C
      REAL TCW      !Thermal conductivity of water (0.57) [W m-1 K-1]
      REAL TCICE    !Thermal conductivity of ice (2.24) [W m-1 K-1]
      REAL TCSAND   !Thermal conductivity of sand particles (2.5) [W m-1 K-1]
      REAL TCCLAY   !Thermal conductivity of fine mineral particles (2.5) [W m-1 K-1]
      REAL TCOM     !Thermal conductivity of organic matter (0.25) [W m-1 K-1]
      REAL TCDRYS   !Thermal conductivity of dry mineral soil (0.25) [W m-1 K-1]
      REAL RHOSOL   !Density of soil mineral matter (2.65*10^3) [kg m-3]
      REAL RHOOM    !Density of soil organic matter (1.30*10^3) [kg m-3]
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
      COMMON /CLASS3/ TCW,TCICE,TCSAND,TCCLAY,TCOM,TCDRYS,
     1                RHOSOL,RHOOM
      COMMON /CLASS4/ HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,
     1                SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,
     2                TCGLAC,CLHMLT,CLHVAP
C-----------------------------------------------------------------------
C
C   In this subroutine, coefficients are derived for an equation 
C   relating the heat flux at the snow surface to the snow surface 
C   temperature. It is assumed that the variation of temperature T with 
C   depth z in the snow pack can be modelled by using a quadratic 
C   equation:
C
C   T(z) = (1/2)*a*z^2 + b*z + c
C
C   By substituting 0 for z in the above equation and in the expressions 
C   for its first and second derivatives, it can be shown that a = 
C   T″(0), b = T′(0), and c = T(0). The term T″(0) can be evaluated from 
C   the expression for the first derivative evaluated at the bottom of 
C   the snow pack, T(delta_zs):
C
C   T″(0) = [T′(ZSNOW) - T′(0)]/ZSNOW
C
C   The temperature gradient T′(0) at the snow surface is related to the 
C   surface heat flux G(0) by the snow
C   thermal conductivity lambda_s:
C
C   G(0) = -lambda_s*T′(0)
C
C   The average snow temperature, TSNOW, can be obtained by integrating 
C   the resulting equation for T(z) between 0 and ZSNOW. Making use of 
C   all of the above expressions, and assuming as a first approximation 
C   that the heat flux at the bottom of the snow pack is zero, a linear 
C   equation can be derived relating G(0) to T(0):
C
C   G(0) = 3*lambda_s/ZSNOW [T(0) - TSNOW]
C
C   Just four calculations are performed in this subroutine. The slope 
C   and intercept of the G(0) vs. T(0) relation, GCOEFF and GCONST, are 
C   evaluated as 3*lambda_s/ZSNOW and -3*lambda_s*TSNOW/ZSNOW 
C   respectively; the flag IWATER is set to 2, indicating a snow 
C   surface; and the latent heat of vaporization at the surface, 
C   CPHCHG, is set to the value for sublimation (by adding the latent 
C   heat of melting to the latent heat of vaporization).
C
C     * CALCULATE COEFFICIENTS.
C
      DO 100 I=IL1,IL2
          IF(FI(I).GT.0.)                                          THEN
              GCOEFFS(I)=3.0*TCSNOW(I)/ZSNOW(I)
              GCONSTS(I)=-3.0*TCSNOW(I)*TSNOW(I)/ZSNOW(I)
              IWATER(I)=2                                                                    
              CPHCHG(I)=CLHVAP+CLHMLT
          ENDIF
  100 CONTINUE
C                                                                                  
      RETURN                                                                      
      END
