!> \file
!! Calculates near surface output variables
!! @author Y. Delage, D. Verseghy, M. Lazare, J. Melton
!
subroutine SLDIAG(SUT,SVT,STT,SQT,CDM,CDH,UA,VA,TA,QA,T0,Q0, &
                   Z0M,Z0E,F,ZA,ZU,ZT,ILG,IL1,IL2,JL)

  !     * JUN 23/14 - M.LAZARE.   New version for gcm18+:
  !     *                         - Bugfix to calculation of
  !     *                           screen temperature and
  !     *                           screen specific humidity.
  !     *                         - Accumulation removed (now
  !     *                           done in classt/oiflux11) so
  !     *                           that a screen relative humidity
  !     *                           can be calculated. Therefore,
  !     *                           "instantaneous" fields are
  !     *                           calculated and passed out
  !     *                           instead.
  !     * SEP 05/12 - J.MELTON.   Made some numbers explicitly reals
  !     * OCT 17/11 - D.VERSEGHY. ADD CODE TO CIRCUMVENT SPECIAL
  !     *                         CASE WHERE TA~T0 OR QA~QO, THUS
  !     *                         AVOIDING A DIVIDE BY ZERO.
  !     * NOV 04/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
  !     * JUL 19/96 - Y. DELAGE.
  !     ------------------------------------
  use classic_params, only : GRAV,VKC,VMIN

  implicit none
  !
  !     * INTEGER CONSTANTS
  !
  integer :: ILG   !< NUMBER OF POINTS TO BE TREATED
  integer :: IL1,IL2,JL,I
  !
  !     * INPUT FIELDS
  !
  real :: SUT(ILG)  !< U COMPONENT OF THE WIND AT ZU
  real :: SVT(ILG)  !< V COMPONENT OF THE WIND AT ZU
  real :: STT(ILG)  !< TEMPERATURE AT ZT
  real :: SQT(ILG)  !< SPECIFIC HUMIDITY AT ZT
  !
  !     * OUTPUT FIELDS
  !
  real :: CDM(ILG) !< DRAG COEFFICIENT
  real :: CDH(ILG) !< TRASFER COEFFICIENT FOR HEAT AND MOISTURE
  real :: UA(ILG)  !< U COMPONENT OF THE WIND AT ZA
  real :: VA(ILG)  !< V COMPONENT OF THE WIND AT ZA
  real :: TA(ILG)  !< POTENTIAL TEMPERATURE AT ZA
  real :: QA(ILG)  !< SPECIFIC HUMIDITY AT REFERENCE HEIGHT
  real :: Z0M(ILG) !< ROUGHNESS LENGTH FOR MOMENTUM
  real :: Z0E(ILG) !< ROUGHNESS LENGTH FOR HEAT AND MOISTURE
  real :: F(ILG)   !< FRACTION OF GRID POINT AFFECTED BY PROCESS
  real :: T0(ILG)  !< TEMPERATURE AT BOTTOM OF SURFACE LAYER
  real :: Q0(ILG)  !< SPECIFIC HUMIDITY AT BOTTOM OF SURFACE LAYER
  real :: ZA(ILG)  !< TOP OF SURFACE LAYER
  real :: ZU(ILG)  !< HEIGHT OF OUTPUT WIND
  real :: ZT(ILG)  !< HEIGHT OF OUTPUT TEMPERATURE AND HUMIDITY
  !
  !     TEMPORARY VARIABLES
  !
  real :: PR,WSPD,CM,US,TS,QS,L,UVA,RATIO,UVU,TTA,CE
  !
  real :: PSM,PSE,Y,PIM,PIE,X

  !     * STABILITY FUNCTIONS FOR THE STABLE CASE

  PSM(X) = - X - .667 * (X - 5.0 / .35) * EXP( - .35 * X)
  PSE(X) = - (1.0 + .667 * X) ** 1.5 - .667 * (X - 5.0 / .35) * EXP( - .35 * X)

  !     * STABILITY FUNCTIONS FOR THE UNSTABLE CASE

  Y(X) = (1.0 - 16.0 * X) ** .25
  PIM(X) = LOG((1.0 + X) ** 2 * (1.0 + X ** 2)) - 2.0 * ATAN(X)
  PIE(X) = 2.0 * LOG(1.0 + X ** 2)

  PR = 1.0
  do I = IL1,IL2 ! loop 100
    if (F(I) > 0.) then

      !     * CALCULATION OF SURFACE FLUXES AND MONIN-OBUKHOV LENGTH

      WSPD = MAX(VMIN,SQRT(UA(I) ** 2 + VA(I) ** 2))
      CM = SQRT(CDM(I))
      US = CM * WSPD

      if (ABS(TA(I) - T0(I)) < 0.01) then
        TS = - 0.01 * CDH(I) / CM
      else
        TS = CDH(I) * (TA(I) - T0(I)) / CM
      end if
      if (ABS(QA(I) - Q0(I)) < 1.0E-7) then
        QS = - 1.0E-7 * CDH(I) / CM
      else
        QS = CDH(I) * (QA(I) - Q0(I)) / CM
      end if

      L = TA(I) * US ** 2 / (VKC * GRAV * (TS * (1.0 + .61 * QA(I)) + .61 * TA(I) * QS))

      !     * CALCULATE CORRECTION FACTORS TO TAKE INTO ACCOUNT THE APPROXIMATIONS
      !     * IN DRCOEF

      if (L > 0.) then

        !     * STABLE CASE

        UVA = US / VKC * (LOG(ZA(I) / Z0M(I)) - PSM(ZA(I) / L) + PSM(Z0M(I) / L))
        RATIO = WSPD / UVA
        UVU = US / VKC * (LOG((ZU(I) + Z0M(I)) / Z0M(I)) - PSM((ZU(I) + Z0M(I)) / L) &
        + PSM(Z0M(I) / L)) * RATIO
        TTA = T0(I) + TS / VKC * PR * (LOG(ZA(I) / Z0E(I)) - PSE(ZA(I) / L) + &
           PSE(Z0E(I) / L))
        RATIO = (TA(I) - T0(I)) / SIGN(MAX(ABS(TTA - T0(I)),1.E-4),TTA - T0(I))
        CE = (LOG((ZT(I) + Z0M(I)) / Z0E(I)) - PSE((ZT(I) + Z0M(I)) / L) &
       + PSE(Z0E(I) / L)) * RATIO * PR / VKC

      else

        !     * UNSTABLE CASE

        UVA = US / VKC * (LOG(ZA(I) / Z0M(I)) - PIM(Y(ZA(I) / L)) + PIM(Y(Z0M(I) / L)))
        RATIO = WSPD / UVA
        UVU = US / VKC * (LOG((ZU(I) + Z0M(I)) / Z0M(I)) - PIM(Y((ZU(I) + Z0M(I)) / L)) &
         + PIM(Y(Z0M(I) / L))) * RATIO

        TTA = T0(I) + TS / VKC * PR * (LOG(ZA(I) / Z0E(I)) - PIE(Y(ZA(I) / L)) + &
           PIE(Y(Z0E(I) / L)))
        RATIO = (TA(I) - T0(I)) / SIGN(MAX(ABS(TTA - T0(I)),1.E-4),TTA - T0(I))
        CE = (LOG((ZT(I) + Z0M(I)) / Z0E(I)) - PIE(Y((ZT(I) + Z0M(I)) / L)) &
       + PIE(Y(Z0E(I) / L))) * RATIO * PR / VKC

      end if
      !
      SUT(I) = UVU * UA(I) / WSPD
      SVT(I) = UVU * VA(I) / WSPD
      STT(I) = T0(I) + TS * CE
      SQT(I) = Q0(I) + QS * CE
    end if
  end do ! loop 100

  return
end
