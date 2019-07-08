!> \file
!>
!! Scatters variables from long, gathered vectors back onto original two-dimensional arrays (latitude
!! circle x mosaic tiles). The suffix ROT refers to variables on original two-dimensional arrays.
!! The suffix GAT refers to variables on gathered long vectors.
!! @author D. Verseghy, M. Lazare
!!
subroutine classScatter (TBARROT,THLQROT,THICROT,TSFSROT,TPNDROT, & ! Formerly CLASSS
                    ZPNDROT,TBASROT,ALBSROT,TSNOROT,RHOSROT, &
                    SNOROT, GTROT,  TCANROT,RCANROT,SCANROT, &
                   GROROT, CMAIROT,TACROT, QACROT, WSNOROT, &
                    REFROT, BCSNROT,EMISROT,SALBROT,CSALROT, &
                    ILMOS,JLMOS, &
                    NML,NL,NT,NM,ILG,IG,IC,ICP1,NBS, &
                    TBARGAT,THLQGAT,THICGAT,TSFSGAT,TPNDGAT, &
                    ZPNDGAT,TBASGAT,ALBSGAT,TSNOGAT,RHOSGAT, &
                    SNOGAT, GTGAT,  TCANGAT,RCANGAT,SCANGAT, &
                    GROGAT, CMAIGAT,TACGAT, QACGAT, WSNOGAT, &
                    REFGAT, BCSNGAT,EMISGAT,SALBGAT,CSALGAT)
  !
  !     * DEC 23/16 - M.LAZARE.  PROMOTE DIMENSIONS OF WSNOROT TO
  !     *                        NLAT,NMOS (FOR LAKE MODEL)
  !     * Jun 20, 2014 - M.Lazare. New version for gcm18, called
  !     *                          by new "sfcproc2":
  !     *                          - Adds SALBGAT/SALBROT and
  !     *                            CSALGAT,CSALROT (need to pass
  !     *                            NBS as well).
  !     *                          - Adds EMISGAT/EMISROT.
  !     *                          - Adds GTGAT/GTROT.
  !     *                          - Adds NT (NTLD in sfcproc2
  !     *                            call) to dimension land-only
  !     *                            ROT fields, consistent with
  !     *                            new comrow12.
  !     *                          - Unused IWMOS,JWMOS removed.
  !     * Jun 12, 2013 - M.Lazare. Previous version for gcm17,
  !     *                          called by "sfcproc".
  !     *                          CLASS scatter routine called by
  !     *                          "sfcproc" in new version gcm17.
  !     * NOTE: This contains the following changes compared to the
  !     *       working temporary version used in conjunction with
  !     *       updates to gcm16 (ie not official):
  !     *         1) {REF,BCSN} added for Maryam's new code.
  !     *         2) GFLX removed.
  !
  !     * OCT 25/11 - M.LAZARE.   REMOVE OPERATIONS ON INTERNAL
  !     *                         ROT ARRAYS (NOW DONE DIRECTLY
  !     *                         GAT->ROW IN SFCPROC).
  !     * OCT 07/11 - M.LAZARE.   REMOVE TSF.
  !     * OCT 05/11 - M.LAZARE.   ADD SFCH.
  !     * OCT 04/11 - M.LAZARE.   REMOVE ITCT.
  !     * MAR 23/06 - D.VERSEGHY. ADD WSNO,FSNO.
  !     * MAR 18/05 - D.VERSEGHY. ADDITIONAL VARIABLES.
  !     * FEB 18/05 - D.VERSEGHY. ADD "TSFS" VARIABLES.
  !     * AUG 05/04 - D.VERSEGHY. ADD NEW DIAGNOSTIC VARIABLES
  !     *                         ILMO, UE AND HBL.
  !     * AUG 15/02 - D.VERSEGHY. SCATTER OPERATION ON CLASS
  !     *                         VARIABLES.
  !
  implicit none
  !
  !     * INTEGER CONSTANTS.
  !
  integer, intent(in) :: NML,NL,NT,NM,ILG,IG,IC,ICP1,NBS
  integer             :: K,L,M
  !
  !     * LAND SURFACE PROGNOSTIC VARIABLES.
  !

  !! Suffix ROT refers to variables on original two-dimensional arrays.
  real, intent(out)    :: SALBROT(NL,NM,NBS) !< All-sky albedo  [  ]
  real, intent(out)    :: CSALROT(NL,NM,NBS) !< Clear-sky albedo  [  ]
  real, intent(out)    :: TBARROT(NL,NT,IG) !< Temperature of soil layers [K]
  real, intent(out)    :: THLQROT(NL,NT,IG) !< Volumetric liquid water content of soil layers \f$[m^3 m^{-3} ]\f$
  real, intent(out)    :: THICROT(NL,NT,IG) !< Volumetric frozen water content of soil layers \f$[m^3 m^{-3} ]\f$
  real, intent(out)    :: TSFSROT(NL,NT,4) !< Ground surface temperature over subarea [K]
  real, intent(out)    :: TPNDROT(NL,NT) !< Temperature of ponded water [K]
  real, intent(out)    :: ZPNDROT(NL,NT) !< Depth of ponded water on surface [m]
  real, intent(out)    :: TBASROT(NL,NT) !< Temperature of bedrock in third soil layer [K]
  real, intent(out)    :: ALBSROT(NL,NM) !< Snow albedo [ ]
  real, intent(out)    :: TSNOROT(NL,NM) !< Snowpack temperature [K]
  real, intent(out)    :: RHOSROT(NL,NM) !< Density of snow \f$[kg m^{-3} ]\f$
  real, intent(out)    :: SNOROT (NL,NM) !< Mass of snow pack \f$[kg m^{-2} ]\f$
  real, intent(out)    :: GTROT  (NL,NM) !< Effective surface black-body temperature  [K]
  real, intent(out)    :: TCANROT(NL,NT) !< Vegetation canopy temperature [K]
  real, intent(out)    :: RCANROT(NL,NT) !< Intercepted liquid water stored on canopy \f$[kg m^{-2} ]\f$
  real, intent(out)    :: SCANROT(NL,NT) !< Intercepted frozen water stored on canopy \f$[kg m^{-2} ]\f$
  real, intent(out)    :: GROROT (NL,NT) !< Vegetation growth index [ ]
  real, intent(out)    :: TACROT (NL,NT) !< Temperature of air within vegetation canopy [K]
  real, intent(out)    :: QACROT (NL,NT) !< Specific humidity of air within vegetation canopy space \f$[kg kg^{-1} ]\f$
  real, intent(out)    :: WSNOROT(NL,NT) !< Liquid water content of snow pack \f$[kg m^{-2} ]\f$
  real, intent(out)    :: CMAIROT(NL,NT) !< Aggregated mass of vegetation canopy \f$[kg m^{-2} ]\f$
  real, intent(out)    :: REFROT (NL,NM) !< Snow grain size  [m]
  real, intent(out)    :: BCSNROT(NL,NM) !< Black carbon mixing ratio \f$[kg m^{-3} ]\f$
  real, intent(out)    :: EMISROT(NL,NM) !< Surface emissivity  [  ]
  !
  real, intent(in)    :: SALBGAT(ILG,NBS) !< All-sky albedo  [  ]
  real, intent(in)    :: CSALGAT(ILG,NBS) !< Clear-sky albedo  [  ]
  real, intent(in)    :: TBARGAT(ILG,IG) !< Temperature of soil layers [K]
  real, intent(in)    :: THLQGAT(ILG,IG) !< Volumetric liquid water content of soil layers \f$[m^3 m^{-3} ]\f$
  real, intent(in)    :: THICGAT(ILG,IG) !< Volumetric frozen water content of soil layers \f$[m^3 m^{-3} ]\f$
  real, intent(in)    :: TSFSGAT(ILG,4) !< Ground surface temperature over subarea [K]
  real, intent(in)    :: TPNDGAT(ILG) !< Temperature of ponded water [K]
  real, intent(in)    :: ZPNDGAT(ILG) !< Depth of ponded water on surface [m]
  real, intent(in)    :: TBASGAT(ILG) !< Temperature of bedrock in third soil layer [K]
  real, intent(in)    :: ALBSGAT(ILG) !< Snow albedo [ ]
  real, intent(in)    :: TSNOGAT(ILG) !< Snowpack temperature [K]
  real, intent(in)    :: RHOSGAT(ILG) !< Density of snow \f$[kg m^{-3} ]\f$
  real, intent(in)    :: SNOGAT (ILG) !< Mass of snow pack \f$[kg m^{-2} ]\f$
  real, intent(in)    :: GTGAT  (ILG) !< Effective surface black-body temperature  [K]
  real, intent(in)    :: TCANGAT(ILG) !< Vegetation canopy temperature [K]
  real, intent(in)    :: RCANGAT(ILG) !< Intercepted liquid water stored on canopy \f$[kg m^{-2} ]\f$
  real, intent(in)    :: SCANGAT(ILG) !< Intercepted frozen water stored on canopy \f$[kg m^{-2} ]\f$
  real, intent(in)    :: GROGAT (ILG) !< Vegetation growth index [ ]
  real, intent(in)    :: TACGAT (ILG) !< Temperature of air within vegetation canopy [K]
  real, intent(in)    :: QACGAT (ILG) !< Specific humidity of air within vegetation canopy space \f$[kg kg^{-1} ]\f$
  real, intent(in)    :: WSNOGAT(ILG) !< Liquid water content of snow pack \f$[kg m^{-2} ]\f$
  real, intent(in)    :: CMAIGAT(ILG) !< Aggregated mass of vegetation canopy \f$[kg m^{-2} ]\f$
  real, intent(in)    :: REFGAT (ILG) !< Snow grain size  [m]
  real, intent(in)    :: BCSNGAT(ILG) !< Black carbon mixing ratio \f$[kg m^{-3} ]\f$
  real, intent(in)    :: EMISGAT(ILG) !< Surface emissivity  [  ]

  !
  !     * GATHER-SCATTER INDEX ARRAYS.
  !
  integer, intent(in) :: ILMOS (ILG),  JLMOS  (ILG)
  !----------------------------------------------------------------------
  do K = 1,NML ! loop 100
    TPNDROT(ILMOS(K),JLMOS(K)) = TPNDGAT(K)
    ZPNDROT(ILMOS(K),JLMOS(K)) = ZPNDGAT(K)
    TBASROT(ILMOS(K),JLMOS(K)) = TBASGAT(K)
    ALBSROT(ILMOS(K),JLMOS(K)) = ALBSGAT(K)
    TSNOROT(ILMOS(K),JLMOS(K)) = TSNOGAT(K)
    RHOSROT(ILMOS(K),JLMOS(K)) = RHOSGAT(K)
    SNOROT (ILMOS(K),JLMOS(K)) = SNOGAT (K)
    GTROT  (ILMOS(K),JLMOS(K)) = GTGAT  (K)
    WSNOROT(ILMOS(K),JLMOS(K)) = WSNOGAT(K)
    TCANROT(ILMOS(K),JLMOS(K)) = TCANGAT(K)
    RCANROT(ILMOS(K),JLMOS(K)) = RCANGAT(K)
    SCANROT(ILMOS(K),JLMOS(K)) = SCANGAT(K)
    GROROT (ILMOS(K),JLMOS(K)) = GROGAT (K)
    TACROT (ILMOS(K),JLMOS(K)) = TACGAT (K)
    QACROT (ILMOS(K),JLMOS(K)) = QACGAT (K)
    CMAIROT(ILMOS(K),JLMOS(K)) = CMAIGAT(K)
    REFROT (ILMOS(K),JLMOS(K)) = REFGAT (K)
    BCSNROT(ILMOS(K),JLMOS(K)) = BCSNGAT(K)
    EMISROT(ILMOS(K),JLMOS(K)) = EMISGAT(K)

    !>
    !! The prognostic variables are scattered from the long, gathered arrays (collapsing the latitude and mosaic
    !! dimensions into one) back onto the original arrays using the pointer vectors generated in classGatherPrep.
    !!

  end do ! loop 100

  do L = 1,NBS ! loop 200
    do K = 1,NML
      SALBROT(ILMOS(K),JLMOS(K),L) = SALBGAT(K,L)
      CSALROT(ILMOS(K),JLMOS(K),L) = CSALGAT(K,L)
    end do
  end do ! loop 200

  do L = 1,IG ! loop 300
    do K = 1,NML
      TBARROT(ILMOS(K),JLMOS(K),L) = TBARGAT(K,L)
      THLQROT(ILMOS(K),JLMOS(K),L) = THLQGAT(K,L)
      THICROT(ILMOS(K),JLMOS(K),L) = THICGAT(K,L)
    end do
  end do ! loop 300

  do L = 1,4 ! loop 400
    do K = 1,NML
      TSFSROT(ILMOS(K),JLMOS(K),L) = TSFSGAT(K,L)
    end do
  end do ! loop 400

  return
end
