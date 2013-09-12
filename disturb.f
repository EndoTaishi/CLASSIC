      subroutine disturb (stemmass, rootmass, gleafmas, bleafmas,
     1                       thliq,   wiltsm,  fieldsm,    uwind,
     2                       vwind,  lightng,  fcancmx, litrmass,     
     3                    prbfrhuc, rmatctem, extnprob, pftareab,
     4                         il1,      il2,     sort, nol2pfts,
     6                    grclarea,    thice,   popdin, lucemcom,
     7                      dofire,   currlat,
c    8 ------------------ inputs above this line ----------------------   
     9                    stemltdt, rootltdt, glfltrdt, blfltrdt,
     a                    pftareaa, glcaemls, rtcaemls, stcaemls,
     b                    blcaemls, ltrcemls, burnfrac, probfire,
     c                    emit_co2, emit_co,  emit_ch4, emit_nmhc,
     d                    emit_h2,  emit_nox, emit_n2o, emit_pm25,
     e                    emit_tpm, emit_tc,  emit_oc,  emit_bc,
     f                    burnvegf, bterm,    mterm,    lterm )
c    g ------------------outputs above this line ----------------------
c
C               Canadian Terrestrial Ecosystem Model (CTEM)
C                           Disturbance Subroutine
c
c
c     25  Jul 2013  - Add in module for common params, cleaned up code 
c     J. Melton       
c
c     24  Sep 2012  - add in checks to prevent calculation of non-present
c     J. Melton       pfts
c
c     09  May 2012  - addition of emission factors and revising of the
c     J. Melton       fire scheme
c
c     15  May 2003  - this subroutine calculates the litter generated
c     V. Arora        and c emissions from leaves, stem, and root 
c                     components due to fire. c emissions from burned
c                     litter are also estimated. at present no other 
c                     form of disturbance is modelled.
c     inputs 
c
c     stemmass  - stem mass for each of the 9 ctem pfts, kg c/m2
c     rootmass  - root mass for each of the 9 ctem pfts, kg c/m2
c     gleafmas  - green leaf mass for each of the 9 ctem pfts, kg c/m2
c     bleafmas  - brown leaf mass 
c     thliq     - liquid soil moisture content
c     wiltsm    - wilting point soil moisture content
c     fieldsm   - field capacity soil moisture content
c     uwind     - wind speed, m/s
c     vwind     - wind speed, m/s
c     lightng   - total lightning, flashes/(km^2.year)
c                 it is assumed that cloud to ground lightning is
c                 some fixed fraction of total lightning.
c     fcancmx   - fractional coverages of ctem's 9 pfts
c     litrmass  - litter mass for each of the 9 pfts
c     prbfrhuc  - probability of fire due to human causes
c     rmatctem  - fraction of roots in each soil layer for each pft
c     extnprob  - fire extinguishing probability
c     pftareab  - areas of different pfts in a grid cell, before fire, km^2
c     il1,il2   - il1=1, il2=ilg
c     sort      - index for correspondence between 9 pfts and size 12 of
c                 parameters vectors
c     nol2pfts  - number of level 2 ctem pfts
c     grclarea  - gcm grid cell area, km^2
c     thice     - frozen soil moisture content over canopy fraction
c     popdin    - population density (people / km^2)
c     lucemcom  - land use change (luc) related combustion emission losses,
c                 u-mol co2/m2.sec 
c     dofire    - boolean, if true allow fire, if false no fire.
c
c     outputs
c
c     stemltdt  - stem litter generated due to disturbance (kg c/m2)
c     rootltdt  - root litter generated due to disturbance (kg c/m2)
c     glfltrdt  - green leaf litter generated due to disturbance (kg c/m2)
c     blfltrdt  - brown leaf litter generated due to disturbance (kg c/m2)
c     burnarea  - total area burned, km^2
c     burnfrac  - total areal fraction burned, (%)
c     probfire  - probability of fire
c     pftareaa  - areas of different pfts in a grid cell, after fire, km^2
c
c     note the following c burned will be converted to a trace gas 
c     emission or aerosol on the basis of emission factors.
c
c     glcaemls  - green leaf carbon emission losses, kg c/m2
c     blcaemls  - brown leaf carbon emission losses, kg c/m2
c     rtcaemls  - root carbon emission losses, kg c/m2
c     stcaemls  - stem carbon emission losses, kg c/m2
c     ltrcemls  - litter carbon emission losses, kg c/m2

c     emission factors for trace gases and aerosols. units are
c     g of compound emitted per kg of dry organic matter.
c     values are taken from li et al. 2012 biogeosci 
c     emif_co2  - carbon dioxide
c     emif_co   - carbon monoxide
c     emif_ch4  - methane
c     emif_nmhc - non-methane hydrocarbons
c     emif_h2   - hydrogen gas
c     emif_nox  - nitrogen oxides
c     emif_n2o  - nitrous oxide
c     emif_pm25 - particulate matter less than 2.5 um in diameter
c     emif_tpm  - total particulate matter
c     emif_tc   - total carbon
c     emif_oc   - organic carbon
c     emif_bc   - black carbon

c     emitted compounds from biomass burning in g of compound
c     emit_co2  - carbon dioxide
c     emit_co   - carbon monoxide
c     emit_ch4  - methane
c     emit_nmhc - non-methane hydrocarbons
c     emit_h2   - hydrogen gas
c     emit_nox  - nitrogen oxides
c     emit_n2o  - nitrous oxide
c     emit_pm25 - particulate matter less than 2.5 um in diameter
c     emit_tpm  - total particulate matter
c     emit_tc   - total carbon
c     emit_oc   - organic carbon
c     emit_bc   - black carbon

c     tot_emit  - sum of all pools to be converted to emissions/aerosols (g c/m2)
c     tot_emit_dom - tot_emit converted to kg dom / m2

c     hb_interm - interm calculation
c     hbratio   - head to back ratio of ellipse
c     burnvegf - total per PFT areal fraction burned, (%)

      use ctem_params,        only : ignd, icc, ilg, ican, zero,
     1                               kk, pi, c2dom
c
      implicit none
c
      integer    il1,        il2,      i,      j,    
     1             k,      iseed,   m,
     2            k1,      k2,        n
c
      integer       sort(icc),      nol2pfts(ican)

      logical dofire, fire(ilg)
c
      real  stemmass(ilg,icc), rootmass(ilg,icc), gleafmas(ilg,icc),
     1      bleafmas(ilg,icc),     thliq(ilg,ignd),    wiltsm(ilg,ignd),
     2        fieldsm(ilg,ignd),        uwind(ilg),        vwind(ilg),
     3       fcancmx(ilg,icc),      lightng(ilg),litrmass(ilg,icc+1),
     4          prbfrhuc(ilg),     extnprob(ilg), pftareab(ilg,icc),
     5   rmatctem(ilg,icc,ignd),     thice(ilg,ignd),            popdin,
     6          lucemcom(ilg),    tmpprob(ilg),  currlat(ilg) 
c
      real  stemltdt(ilg,icc), rootltdt(ilg,icc), glfltrdt(ilg,icc),
     1          burnarea(ilg), pftareaa(ilg,icc), glcaemls(ilg,icc),
     2      rtcaemls(ilg,icc), stcaemls(ilg,icc), ltrcemls(ilg,icc),
     3      blfltrdt(ilg,icc), blcaemls(ilg,icc),     burnfrac(ilg),
     4     emit_co2(ilg,icc),   emit_co(ilg,icc), emit_ch4(ilg,icc),
     5    emit_nmhc(ilg,icc),   emit_h2(ilg,icc), emit_nox(ilg,icc),
     6     emit_n2o(ilg,icc), emit_pm25(ilg,icc), emit_tpm(ilg,icc),
     7      emit_tc(ilg,icc),   emit_oc(ilg,icc),  emit_bc(ilg,icc)
c
      real        bmasthrs(2),                              
     1               extnmois,          lwrlthrs,         hgrlthrs,
     2               parmlght,          parblght,            alpha,
     3                     f0,       maxsprd(kk),      frco2lf(kk),
     4            frltrlf(kk),      frco2stm(kk),     frltrstm(kk),
     5            frco2rt(kk),       frltrrt(kk),     frltrbrn(kk),
     7           emif_co2(kk),       emif_co(kk),     emif_ch4(kk),
     8          emif_nmhc(kk),       emif_h2(kk),     emif_nox(kk),
     9           emif_n2o(kk),     emif_pm25(kk),     emif_tpm(kk),
     a            emif_tc(kk),       emif_oc(kk),      emif_bc(kk)
c
      real   biomass(ilg,icc),        bterm(ilg), drgtstrs(ilg,icc),
     1       betadrgt(ilg,ignd),     avgdryns(ilg),        fcsum(ilg),
     2          avgbmass(ilg),        mterm(ilg),     c2glgtng(ilg),
     3          betalght(ilg),            y(ilg),              ymin,
     4                   ymax,             slope,        lterm(ilg),
     5          probfire(ilg),             ctime,            random,
     6                   temp,     betmsprd(ilg),       smfunc(ilg),
     7              wind(ilg),      wndfunc(ilg),     sprdrate(ilg),
     8           lbratio(ilg),     arbn1day(ilg),     areamult(ilg),
     9       burnveg(ilg,icc),      vegarea(ilg),     grclarea(ilg),
     a                reparea,          tot_emit,      tot_emit_dom,
     b       burnvegf(ilg,icc)   

      real          hb_interm,      hbratio(ilg),       popdthrshld,
     1                 fden_m
c
c     ------------------------------------------------------------------
c                     Constants used in the model
c     Note the structure of vectors which shows the CLASS
c     PFTs (along rows) and CTEM sub-PFTs (along columns)
c
c     needle leaf |  evg       dcd       ---
c     broad leaf  |  evg   dcd-cld   dcd-dry
c     crops       |   c3        c4       ---
c     grasses     |   c3        c4       ---
c
c     min. and max. vegetation biomass thresholds to initiate fire, kg c/m^2
c      data bmasthrs/0.2, 1.0/ !ORIG
      data bmasthrs/0.4, 1.2/ ! testing
!     **Lowering the bmasthrs gets too much burning outside of the tropics, 
!     however maybe this higher value prevents the savannahs from burning more? -JM

c     extinction moisture content for estimating fire likeliness due
c     to soil moisture
      data extnmois/0.30/ !ORIG
c
c     lower cloud-to-ground lightning threshold for fire likelihood
c     flashes/km^2.year
      data lwrlthrs/0.25/
c
c     higher cloud-to-ground lightning threshold for fire likelihood
c     flashes/km^2.year
      data hgrlthrs/10.0/ !ORIG
c
c     parameter m (mean) and b of logistic distribution used for 
c     estimating fire likelihood due to lightning
      data parmlght/0.8/  !ORIG 0.4
      data parblght/0.1/
!     **Parmlght was increased to 0.8 to make it so areas with higher amounts of
!     lightning have higher lterm. The saturation is still the same, but the 
!     increase is more gradual at low lightning density. JM

!     **These values were being calculated each time, they shouldn't be as they
!     are just parameters. JM
      data ymin/0.01798621/     !=1.0/( 1.0+exp((parmlght-0.0)/parblght) )
      data ymax/0.9975273768/   !=1.0/( 1.0+exp((parmlght-1.0)/parblght) )
      data slope/0.0204588331/  !=abs(0.0-ymin)+abs(1.0-ymax)
c
c     parameter alpha and f0 used for estimating wind function for
c     fire spread rate
      data alpha/8.16326e-04/

c     Fire spread rate in the absence of wind
c      data f0/0.1/ !ORIG
c     new: 
c     **Now a derived quantity from the formulation of the wind speed fire spread
c     rate scalar
      data f0/0.05/
c
c     max. fire spread rate, km/hr
c     old:
c     data maxsprd/0.45/
c     **Li et al. values:
      data maxsprd/0.54, 0.54, 0.00,
     &             0.40, 0.40, 0.40,
     &             0.00, 0.00, 0.00,
     &             0.72, 0.72, 0.00/

c     fraction of leaf biomass converted to gases due to combustion
c      data frco2lf/0.70, 0.70, 0.00, !ORIG
c     &             0.70, 0.70, 0.70,
c     &             0.00, 0.00, 0.00,
c     &             0.80, 0.80, 0.00/

      data frco2lf/0.60, 0.60, 0.00, !TEST
     &             0.60, 0.60, 0.60,
     &             0.00, 0.00, 0.00,
     &             0.70, 0.70, 0.00/

c
c     fraction of leaf biomass becoming litter after combustion
      data frltrlf/0.20, 0.20, 0.00,
     &             0.20, 0.20, 0.20,
     &             0.00, 0.00, 0.00,
     &             0.10, 0.10, 0.00/
c
c     fraction of stem biomass converted to gases due to combustion
c      data frco2stm/0.20, 0.20, 0.00, !ORIG
c     &              0.20, 0.10, 0.10,
c     &              0.00, 0.00, 0.00,
c     &              0.00, 0.00, 0.00/

      data frco2stm/0.15, 0.15, 0.00, !TEST
     &              0.15, 0.05, 0.05,
     &              0.00, 0.00, 0.00,
     &              0.00, 0.00, 0.00/

c
c     fraction of stem biomass becoming litter after combustion
c      data frltrstm/0.60, 0.60, 0.00,  !ORIG
c     &              0.60, 0.40, 0.40,
c     &              0.00, 0.00, 0.00,
c     &              0.00, 0.00, 0.00/

      data frltrstm/0.50, 0.50, 0.00,  !TEST
     &              0.50, 0.30, 0.30,
     &              0.00, 0.00, 0.00,
     &              0.00, 0.00, 0.00/

c
c     fraction of root biomass converted to gases due to combustion
      data frco2rt/0.0, 0.0, 0.0,
     &             0.0, 0.0, 0.0, 
     &             0.0, 0.0, 0.0, 
     &             0.0, 0.0, 0.0/
c
c     fraction of root biomass becoming litter after combustion
      data frltrrt/0.10, 0.10, 0.00,
     &             0.10, 0.10, 0.10,
     &             0.00, 0.00, 0.00,
     &             0.25, 0.25, 0.00/
c
c     fraction of litter burned during fire and emitted as gases
      data frltrbrn/0.50, 0.50, 0.00,
     &              0.60, 0.60, 0.60,
     &              0.00, 0.00, 0.00,
     &              0.70, 0.70, 0.00/
c
c     ========================

c     emissions factors by chemical species
c     
c     Values are from Andreae 2011 as described in Li et al. 2012
c     Biogeosci. Units: g species / (kg DOM)

c     pft-specific emission factors for CO2 
c     Andreae 2011 as described in Li et al. 2012
c      data emif_co2/1576.0, 1576.0,   0.00,
c     &              1604.0, 1576.0, 1654.0,
c     &              1576.0, 1654.0,   0.00,
c     &              1576.0, 1654.0,   0.00/

c    values from Wiedinmyer et al. 2011
      data emif_co2/1514.0, 1514.0,   0.00,
     &              1643.0, 1630.0, 1716.0,
     &              1537.0, 1537.0,   0.00,
     &              1692.0, 1692.0,   0.00/

c     pft-specific emission factors for CO 
c     Andreae 2011 as described in Li et al. 2012
c      data emif_co /106.0, 106.0, 0.00,
c     &              103.0, 106.0, 64.0,
c     &              106.0,  64.0, 0.00,
c     &              106.0,  64.0, 0.00/

c    values from Wiedinmyer et al. 2011
      data emif_co /118.0, 118.0, 0.00,
     &               92.0, 102.0, 68.0,
     &              111.0, 111.0, 0.00,
     &               59.0,  59.0, 0.00/

c     pft-specific emission factors for CH4 
c     Andreae 2011 as described in Li et al. 2012
c      data emif_ch4/ 4.8, 4.8, 0.0,
c     &               5.8, 4.8, 2.4,
c     &               4.8, 2.4, 0.0,
c     &               4.8, 2.4, 0.0/

c    values from Wiedinmyer et al. 2011
      data emif_ch4/ 6.0, 6.0, 0.0,
     &               5.1, 5.0, 2.6,
     &               6.0, 6.0, 0.0,
     &               1.5, 1.5, 0.0/

c     pft-specific emission factors for NMHC
c     Andreae 2011 as described in Li et al. 2012
c      data emif_nmhc/ 5.7, 5.7, 0.0,
c     &                6.4, 5.7, 3.7,
c     &                5.7, 3.7, 0.0,
c     &                5.7, 3.7, 0.0/

c    values from Wiedinmyer et al. 2011
      data emif_nmhc/ 5.7, 5.7, 0.0,
     &                1.7, 5.7, 3.4,
     &                7.0, 7.0, 0.0,
     &                3.4, 3.4, 0.0/

c     pft-specific emission factors for H2
c     Andreae 2011 as described in Li et al. 2012
c      data emif_h2/ 1.80, 1.80, 0.00,
c     &              2.54, 1.80, 0.98,
c     &              1.80, 0.98, 0.00,
c     &              1.80, 0.98, 0.00/

c    values from Wiedinmyer et al. 2011
      data emif_h2/ 2.30, 2.30, 0.00,
     &              3.20, 1.80, 0.97,
     &              2.40, 2.40, 0.00,
     &              0.97, 0.97, 0.00/

c     pft-specific emission factors for NOx
c     Andreae 2011 as described in Li et al. 2012
c      data emif_nox/3.24, 3.24, 0.00,
c     &              2.90, 3.24, 2.49,
c     &              3.24, 2.49, 0.00,
c     &              3.24, 2.49, 0.00/

c    values from Wiedinmyer et al. 2011 (species: "NOx (as NO)" from Table 1)
      data emif_nox/1.80, 2.30, 0.00,
     &              2.60, 1.30, 3.90,
     &              3.50, 3.50, 0.00,
     &              2.80, 2.80, 0.00/

c     pft-specific emission factors for N2O
c     Andreae 2011 as described in Li et al. 2012 
      data emif_n2o/0.26, 0.26, 0.00,
     &              0.23, 0.26, 0.20,
     &              0.26, 0.20, 0.00,
     &              0.26, 0.20, 0.00/

c     emission factors for aerosols

c     pft-specific emission factors for PM2.5
c     (particles less than 2.5 micrometers in 
c     diameter)
c     Andreae 2011 as described in Li et al. 2012
c      data emif_pm25/12.7, 12.7, 0.0,
c     &               10.5, 12.7, 5.2,
c     &               12.7,  5.2, 0.0,
c     &               12.7,  5.2, 0.0/

c    values from Wiedinmyer et al. 2011
      data emif_pm25/13.0, 13.0, 0.0,
     &                9.7, 13.0, 9.3,
     &                5.8,  5.8, 0.0,
     &                5.4,  5.4, 0.0/

c     pft-specific emission factors for TPM 
c     (total particulate matter)
c     Andreae 2011 as described in Li et al. 2012
c      data emif_tpm/17.6, 17.6, 0.0,
c     &              14.7, 17.6, 8.5,
c     &              17.6,  8.5, 0.0,
c     &              17.6,  8.5, 0.0/

c    values from Wiedinmyer et al. 2011
      data emif_tpm/18.0, 18.0, 0.0,
     &              13.0, 18.0,15.4,
     &              13.0, 13.0, 0.0,
     &               8.3,  8.3, 0.0/

c     pft-specific emission factors for TC
c     (total carbon)
c     Andreae 2011 as described in Li et al. 2012
c      data emif_tc/ 8.3, 8.3, 0.0,
c     &              7.2, 8.3, 3.4,
c     &              8.3, 3.4, 0.0,
c     &              8.3, 3.4, 0.0/

c    values from Wiedinmyer et al. 2011 (TPC in Table 1)
      data emif_tc/ 8.3, 8.3, 0.0,
     &              5.2, 9.7, 7.1,
     &              4.0, 4.0, 0.0,
     &              3.0, 3.0, 0.0/

c     pft-specific emission factors for OC 
c     (organic carbon)
c     Andreae 2011 as described in Li et al. 2012
c      data emif_oc/ 9.1, 9.1, 0.0,
c     &              6.7, 9.1, 3.2,
c     &              9.1, 3.2, 0.0,
c     &              9.1, 3.2, 0.0/

c    values from Wiedinmyer et al. 2011 
      data emif_oc/ 7.8, 7.8, 0.0,
     &              4.7, 9.2, 6.6,
     &              3.3, 3.3, 0.0,
     &              2.6, 2.6, 0.0/

c     pft-specific emission factors for BC
c     (black carbon)
c     Andreae 2011 as described in Li et al. 2012
c      data emif_bc/ 0.56, 0.56, 0.00,
c     &              0.56, 0.56, 0.47,
c     &              0.56, 0.47, 0.00,
c     &              0.56, 0.47, 0.00/

c    values from Wiedinmyer et al. 2011
      data emif_bc/ 0.20, 0.20, 0.00,
     &              0.52, 0.56, 0.50,
     &              0.69, 0.69, 0.00,
     &              0.37, 0.37, 0.00/


c     typical area representing ctem's fire parameterization
      data reparea/1000.0/ ! km^2
c
c     threshold of population density (people/km2) [kloster et al., biogeosci. 2010]
      data popdthrshld/300./

c     ---------------------------------------------------------------
c
C     * if icc /= 9 or ignd /= 3 this subroutine will need changes.
      IF(ICC.NE.9)                            CALL XIT('DISTURB',-1)
      IF(IGND.NE.3)                            CALL XIT('DISTURB',-2)

c     initialize required arrays to zero
c
      do 140 j = 1,icc
        do 150 i = il1, il2
          stemltdt(i,j)=0.0     !stem litter due to disturbance
          rootltdt(i,j)=0.0     !root litter due to disturbance
          glfltrdt(i,j)=0.0     !green leaf litter due to disturbance
          blfltrdt(i,j)=0.0     !brown leaf litter due to disturbance
          biomass(i,j)=0.0      !total biomass for fire purposes
          drgtstrs(i,j)=0.0     !soil dryness factor for pfts
          burnveg(i,j)=0.0      !burn area for each pft
          pftareaa(i,j)=0.0     !pft area after fire
          glcaemls(i,j)=0.0     !green leaf carbon fire emissions
          blcaemls(i,j)=0.0     !brown leaf carbon fire emissions
          stcaemls(i,j)=0.0     !stem carbon fire emissions
          rtcaemls(i,j)=0.0     !root carbon fire emissions
          ltrcemls(i,j)=0.0     !litter carbon fire emissions

          emit_co2(i,j) = 0.0
          emit_co(i,j) = 0.0
          emit_ch4(i,j) = 0.0
          emit_nmhc(i,j) = 0.0
          emit_h2(i,j) = 0.0
          emit_nox(i,j) = 0.0
          emit_n2o(i,j) = 0.0
          emit_pm25(i,j) = 0.0
          emit_tpm(i,j) = 0.0
          emit_tc(i,j) = 0.0
          emit_oc(i,j) = 0.0
          emit_bc(i,j) = 0.0
          burnvegf(i,j)=0.0

150     continue                  
140   continue
c
      do 160 k = 1,ignd
        do 170 i = il1, il2
          betadrgt(i,k)=0.0     !dryness term for soil layers
170     continue                  
160   continue
c
      do 180 i = il1, il2
        avgbmass(i)=0.0         !avg. veg. biomass over the veg. fraction 
c                               !of grid cell
        avgdryns(i)=0.0         !avg. dryness over the vegetated fraction
        fcsum(i)=0.0            !total vegetated fraction
        bterm(i)=0.0            !biomass fire probability term
        mterm(i)=0.0            !moisture fire probability term
        c2glgtng(i)=0.0         !cloud-to-ground lightning
        betalght(i)=0.0         !0-1 lightning term
        y(i)=0.0                !logistic dist. for fire prob. due to lightning
        lterm(i)=0.0            !lightning fire probability term
        probfire(i)=0.0         !probability of fire
        fire(i)=.false.         !fire occuring 
        burnarea(i)=0.0         !total area burned due to fire
        burnfrac(i)=0.0         !total areal fraction burned due to fire
        betmsprd(i)=0.0         !beta moisture for calculating fire spread rate
        smfunc(i)=0.0           !soil moisture function used for fire spread rate
        wind(i)=0.0             !wind speed in km/hr
        wndfunc(i)=0.0          !wind function for fire spread rate
        sprdrate(i)=0.0         !fire spread rate
        lbratio(i)=0.0          !length to breadth ratio of fire
        arbn1day(i)=0.0         !area burned in 1 day, km^2
        areamult(i)=0.0         !multiplier to find area burned
        vegarea(i)=0.0          !total vegetated area in a grid cell

180   continue

c     if not simulating fire, leave the subroutine now.
      if (.not. dofire) goto 600

c
      do 190 i = il1, il2
        if(extnprob(i).le.zero) then
          write(6,*)'fire extinguishing prob. (',i,'= ',extnprob(i)
          write(6,*)'please use an appropriate value of this paramater'
          write(6,*)'else the whole grid cell will burn down leading to'
          write(6,*)'numerical problems.'
          call xit('disturb',-2)
        endif
190   continue
c
c     initialization ends    
c
c     ------------------------------------------------------------------
c
c     Find the probability of fire as a product of three functions
c     with dependence on total biomass, soil moisture, and lightning 
c
c     1. Dependence on total biomass
c
      do 200 j = 1, icc
        do 210 i = il1, il2

         if (fcancmx(i,j) .gt. 0.0) then
c         Root biomass is not used to initiate fire. for example if
c         the last fire burned all grass leaves, and some of the roots
c         were left, its unlikely these roots could catch fire. 
          biomass(i,j)=gleafmas(i,j)+bleafmas(i,j)+stemmass(i,j)+
     &                 litrmass(i,j)

c        Find average biomass over the vegetated fraction
          avgbmass(i) = avgbmass(i)+biomass(i,j)*fcancmx(i,j)
         endif

210     continue
200   continue
c
      do 250 i = il1, il2
        do j = 1, icc
         fcsum(i)=fcsum(i) + fcancmx(i,j)
        end do        

        if(fcsum(i) .gt. zero)then
          avgbmass(i)=avgbmass(i)/fcsum(i)
        else
          avgbmass(i)=0.0
        endif
c
        if(avgbmass(i) .ge. bmasthrs(2))then
          bterm(i)=1.0
        else if(avgbmass(i) .lt. bmasthrs(2) .and.
     &          avgbmass(i) .gt. bmasthrs(1)) then
            bterm(i)=(avgbmass(i)-bmasthrs(1))/(bmasthrs(2)-bmasthrs(1))     
        else if(avgbmass(i) .le. bmasthrs(1)) then
            bterm(i)=0.0  !no fire if biomass below the lower threshold
        endif
        bterm(i)=max(0.0, min(bterm(i),1.0))
250   continue 

c     2. dependence on soil moisture
c
c     This is calculated in a way such that the more dry the root zone
c     of a pft type is, and more fractional area is covered with that
c     pft, the more likely it is that fire will get started. that is
c     the dryness factor is weighted by fraction of roots in soil
c     layers, as well as according to the fractional coverage of 
c     different pfts. the assumption here is that if there is less 
c     moisture in root zone, then it is more likely the vegetation 
c     will be dry and thus the likeliness of fire is more.
c
c     First find the dryness factor for each soil layer.
c
      do 300 j = 1, ignd
        do 310 i = il1, il2
c
          if((thliq(i,j)+thice(i,j)) .le. wiltsm(i,j)) then
            betadrgt(i,j)=0.0
          else if((thliq(i,j)+thice(i,j)) .gt. wiltsm(i,j) .and.
     &            (thliq(i,j)+thice(i,j)) .lt. fieldsm(i,j)) then
            betadrgt(i,j)=(thliq(i,j)+thice(i,j)-wiltsm(i,j))
            betadrgt(i,j)=betadrgt(i,j)/(fieldsm(i,j)-wiltsm(i,j))   
          else
            betadrgt(i,j)=1.0
          endif
          betadrgt(i,j)=max(0.0, min(betadrgt(i,j),1.0))

310     continue
300   continue
c
c     Now find weighted value of this dryness factor averaged over 
c     the rooting depth, for each pft
c
      do 320 j = 1, icc
        do 330 i = il1, il2
         if (fcancmx(i,j) .gt. 0.0) then
          drgtstrs(i,j) =  (betadrgt(i,1))*rmatctem(i,j,1) +
     &                    (betadrgt(i,2))*rmatctem(i,j,2) +
     &                    (betadrgt(i,3))*rmatctem(i,j,3)
          drgtstrs(i,j) = drgtstrs(i,j) /
     &     (rmatctem(i,j,1)+rmatctem(i,j,2)+rmatctem(i,j,3))
          drgtstrs(i,j)=max(0.0, min(drgtstrs(i,j),1.0))
c        Next find this dryness factor averaged over the vegetated fraction
          avgdryns(i) = avgdryns(i)+drgtstrs(i,j)*fcancmx(i,j)
         endif
330     continue
320   continue

c     Use average root zone vegetation dryness to find likelihood of
c     fire due to moisture. 
c
      do 380 i = il1, il2
        if(fcsum(i) .gt. zero)then

          avgdryns(i)=avgdryns(i)/fcsum(i)
          mterm(i)=1.0-tanh((1.75*avgdryns(i)/extnmois)**2)

        else
          mterm(i)=0.0   !no fire likelihood due to moisture if no vegetation
          avgdryns(i)=0.0
        endif
        mterm(i)=max(0.0, min(mterm(i),1.0))
380   continue
c
c     3. dependence on lightning
c
c     Dependence on lightning is modelled in a simple way which implies that
c     a large no. of lightning flashes are more likely to cause fire than
c     few lightning flashes.
c
      do 400 i = il1, il2

        !c2glgtng(i)=0.25*lightng(i)   !ORIG
        !c2glgtng(i)=(1./(4.16+2.16*cos(currlat(i)*pi/180)))*lightng(i) !From Prentice and Mackerras
!       ** New version from Price and Rind. More complete dataset than Prentice and Mackerras.          
        c2glgtng(i)=0.219913*exp(0.0058899*abs(currlat(i)))*lightng(i)  ! approximation of Price and Rind equation.

        if( c2glgtng(i) .le. lwrlthrs) then
           betalght(i)=0.0
        else if ( c2glgtng(i) .gt. lwrlthrs .and.
     &  c2glgtng(i) .lt. hgrlthrs) then 
           betalght(i)=(c2glgtng(i)-lwrlthrs)/(hgrlthrs-lwrlthrs)
        else if (c2glgtng(i) .ge. hgrlthrs) then
           betalght(i)=1.0
        endif

        y(i)=1.0/( 1.0+exp((parmlght-betalght(i))/parblght) )
!        No need to calculate each time, moved to parameters. JM. Sep 4 2013.
!        ymin=1.0/( 1.0+exp((parmlght-0.0)/parblght) )
!        ymax=1.0/( 1.0+exp((parmlght-1.0)/parblght) )
!        slope=abs(0.0-ymin)+abs(1.0-ymax)
        temp=y(i)+(0.0-ymin)+betalght(i)*slope

c     Determine the probability of fire due to human causes
c     this is based upon the population density from the .popd
c     read-in file

!       FLAG testing. See old way below.

!       ** Attempt to raise the probability of ignition for the tropics.
        if (currlat(i) < 30. .and. currlat(i) > -30.) then

         tmpprob(i)=0.
!        ** First try - increase prbfrhuc to 1 always for grasslands.
         do j = 1, icc 
          if (j .eq. 8 .or. j .eq. 9) then !grass
            tmpprob(i)=tmpprob(i) + 1.0  * fcancmx(i,j)
          else  ! not grass - treat as per normal.
            tmpprob(i)=tmpprob(i) + min(1.0,(popdin/popdthrshld)**0.43)
     &                * fcancmx(i,j)
          end if !grass/tree
         end do !j loop

          prbfrhuc(i) = tmpprob(i)

          lterm(i)=temp+(1.0-temp)*prbfrhuc(i)
          lterm(i)=max(0.0, min(lterm(i),1.0))

!         ** Now we try to also raise lterm for grasslands (>50% grass cover)
!         and during the dry season.  
          if (fcancmx(i,8)+fcancmx(i,9) .gt. 0.50 .and. mterm(i) 
     &                                   .gt. 0.8) then
            lterm(i)=max(0.8,min(lterm(i),1.0))
          end if 

        else !extratropics. Treat as per normal.

          prbfrhuc(i)=min(1.0,(popdin/popdthrshld)**0.43) !new way

          lterm(i)=temp+(1.0-temp)*prbfrhuc(i)
          lterm(i)=max(0.0, min(lterm(i),1.0))

        end if

!       OLD WAY:
!        prbfrhuc(i)=min(1.0,(popdin/popdthrshld)**0.43) !Kloster way
!        prbfrhuc(i)=0.5 !ORIG (more or less) 

!        lterm(i)=temp+(1.0-temp)*prbfrhuc(i)
!        lterm(i)=max(0.0, min(lterm(i),1.0))


400   continue
c
c     Multiply the bterm, mterm, and the lterm to find probability of
c     fire.

      do 420 i = il1, il2
        probfire(i)=bterm(i)*mterm(i)*lterm(i)
        if (probfire(i) .gt. zero) fire(i)=.true.

420   continue
c
c     if fire is to be started then estimate burn area and litter generated
c     by the fire, else do nothing.

      do 430 i = il1, il2
        if(fire(i))then
c
c         find spread rate as a function of wind speed and soil moisture in the
c         root zone (as found above) which we use as a surrogate for moisture
c         content of vegetation.
          if( avgdryns(i) .gt. extnmois )then
            betmsprd(i)= 1.0   
          else
            betmsprd(i)= avgdryns(i)/extnmois   
          endif
          smfunc(i)=(1.0-betmsprd(i))**2.0
          wind(i)=sqrt(uwind(i)**2.0 + vwind(i)**2.0)
          wind(i)=wind(i)*3.60     ! change m/s to km/hr

c         length to breadth ratio of fire
c         orig a&b paper value:
c          lbratio(i)=1.0+10.0*(1.0-exp(-0.017*wind(i))) !ORIG
c         Li et al. value, derived quantity:
          lbratio(i)=1.0+10.0*(1.0-exp(-0.06*wind(i)))

c         flag testing:
c         **New. Calculate the head to back ratio of the fire
          hb_interm = (lbratio(i)**2 - 1.0)**0.5
          hbratio(i) = (lbratio(i) + hb_interm)/(lbratio(i) - hb_interm)

c         flag testing:
c         Following Li et al. 2012 this function has been derived
c         from the fire rate spread perpendicular to the wind 
c         direction, in the downwind direction, the head to back
c         ratio and the length to breadth ratio. f0 is also now
c         a derived quantity (0.05)

c         old:
c          wndfunc(i)=1.0 - ( (1.0-f0)*exp(-1.0*alpha*wind(i)**2) ) !ORIG

c         new:
          wndfunc(i)= (2.0 * lbratio(i)) / (1.0 + 1.0 
     &                     / hbratio(i)) * f0 

c         flag: try a pft-dependent spread rate
c         old:
c          sprdrate(i)=maxsprd* smfunc(i)* wndfunc(i)
c         new:
      do 435 j = 1, icc
          n = sort(j)  
          sprdrate(i)=sprdrate(i) + maxsprd(n)* fcancmx(i,j) * smfunc(i)
     &                            * wndfunc(i) 
435     continue

c
c         Area burned in 1 day, km^2

c         The assumed ratio of the head to back ratio was 5.0 in 
c         Arora and Boer 2005 JGR, this value can be calculated
c         as was pointed out in Li et al. 2012 Biogeosci. We 
c         adopt the calculated version below
c         flag testing:
c         old:
c          arbn1day(i)=(pi*0.36*24*24*sprdrate(i)**2)/lbratio(i) !ORIG
c         new:
          arbn1day(i)=(pi*24.0*24.0*sprdrate(i)**2)/(4.0 * lbratio(i)) 
     &                  *(1.0 + 1.0 / hbratio(i))**2

c
c         based on fire extinguishing probability we estimate the number 
c         which needs to be multiplied with arbn1day to estimate average 
c         area burned

c         flag testing: Kloster way
c         fire extinction is based upon population density
           extnprob(i)=max(0.0,0.9-exp(-0.025*popdin))
           extnprob(i)=0.5+extnprob(i)/2.0

c          extnprob(i)=0.5 !ORIG

          areamult(i)=((1.0-extnprob(i))*(2.0-extnprob(i)))/
     &      extnprob(i)**2                              
c
c         area burned, km^2
          burnarea(i)=arbn1day(i)*areamult(i)*(grclarea(i)/reparea)
     &                                                 *probfire(i)
          burnfrac(i)=100.*burnarea(i)/grclarea(i)
c
        endif
430   continue
c
c     make sure area burned is not greater than the vegetated area. 
c     distribute burned area equally amongst pfts present in the grid cell.
c 
      do 460 i = il1, il2

        do j = 1,icc
        vegarea(i)= vegarea(i) + pftareab(i,j)
        end do  
 
        if(burnarea(i) .gt. vegarea(i)) then
          burnarea(i)=vegarea(i)
          burnfrac(i)=100.*burnarea(i)/grclarea(i)
        endif
460   continue
c
      k1=0
      do 470 j = 1, ican
       if(j.eq.1) then
         k1 = k1 + 1
       else
         k1 = k1 + nol2pfts(j-1)
       endif
       k2 = k1 + nol2pfts(j) - 1
       do 475 m = k1, k2
        do 480 i = il1, il2
         if (fcancmx(i,m) .gt. 0.0) then
          if(vegarea(i) .gt. zero)then
            burnveg(i,m)= (burnarea(i)*pftareab(i,m)/vegarea(i))
            if(j .eq. 3)then  !crops not allowed to burn
              burnveg(i,m)= 0.0
            endif
          else
            burnveg(i,m)= 0.0
          endif
          pftareaa(i,m)=pftareab(i,m)-burnveg(i,m)
         endif
480     continue
475    continue 
470   continue 
c
      do 490 i = il1, il2
       burnarea(i) = 0.0
       do j = 1,icc
         burnarea(i)= burnarea(i) + burnveg(i,j)
       end do
       burnfrac(i)=100.*burnarea(i)/grclarea(i)


490   continue
c
c     Check that the sum of fraction of leaves, stem, and root 
c     that needs to be burned and converted into CO2, and fraction that 
c     needs to become litter doesn't exceed one.
c
      do 500 j = 1, icc
        n = sort(j)
        if( (frco2lf(n)+frltrlf(n)) .gt. 1.0 )    call xit('disturb',-3)
        if( (frco2stm(n)+frltrstm(n)) .gt. 1.0 )  call xit('disturb',-4)
        if( (frco2rt(n)+frltrrt(n)) .gt. 1.0 )    call xit('disturb',-5)
        if(  frltrbrn(n) .gt. 1.0 )               call xit('disturb',-6)
500   continue
c
c     Finally estimate amount of litter generated from each pft, and
c     each vegetation component (leaves, stem, and root) based on their
c     resistance to combustion. we also estimate co2 emissions from each
c     of these components. note that the litter which is generated due 
c     to disturbance is uniformly distributed over the entire area of 
c     a given pft, and this essentially thins the vegetation biomass. 
c     at this stage we do not make the burn area bare, and therefore
c     fire doesn't change the fractional coverage of pfts.  
c
      do 520 j = 1, icc
        n = sort(j)
        do 530 i = il1, il2
         if (fcancmx(i,j) .gt. 0.0) then
c
          if(pftareab(i,j) .gt. zero)then

            glfltrdt(i,j)=frltrlf(n) *gleafmas(i,j)*
     &       (burnveg(i,j)/pftareab(i,j))
c
            blfltrdt(i,j)=frltrlf(n) *bleafmas(i,j)*
     &       (burnveg(i,j)/pftareab(i,j))
c
            stemltdt(i,j)=frltrstm(n)*stemmass(i,j)*
     &       (burnveg(i,j)/pftareab(i,j))
c
            rootltdt(i,j)=frltrrt(n) *rootmass(i,j)*
     &       (burnveg(i,j)/pftareab(i,j))
c
            glcaemls(i,j)=frco2lf(n) *gleafmas(i,j)*
     &       (burnveg(i,j)/pftareab(i,j))
c
            blcaemls(i,j)=frco2lf(n) *bleafmas(i,j)*
     &       (burnveg(i,j)/pftareab(i,j))
c
            stcaemls(i,j)=frco2stm(n)*stemmass(i,j)*
     &       (burnveg(i,j)/pftareab(i,j))
c
            rtcaemls(i,j)=frco2rt(n) *rootmass(i,j)*
     &       (burnveg(i,j)/pftareab(i,j))
c
            ltrcemls(i,j)=frltrbrn(n)*litrmass(i,j)*
     &       (burnveg(i,j)/pftareab(i,j))

c          calculate the emissions of trace gases and aerosols based upon how
c          much plant matter was burnt

c          sum all pools that will be converted to emissions/aerosols (g c/m2)
           tot_emit = (glcaemls(i,j) + blcaemls(i,j) + rtcaemls(i,j)
     &            + stcaemls(i,j) + ltrcemls(i,j)) * 1000.0

c          add in the emissions due to luc fires (deforestation)
c          the luc emissions are converted from umol co2 m-2 s-1
c          to g c m-2 (day-1) before adding to tot_emit         
           tot_emit = tot_emit + (lucemcom(i) / 963.62 * 1000.0)

c          convert burnt plant matter from carbon to dry organic matter using 
c          a conversion factor, assume all parts of the plant has the same
c          ratio of carbon to dry organic matter. units: kg dom / m2
           tot_emit_dom = tot_emit / c2dom

c          convert the dom to emissions/aerosols using emissions factors
c          units: g compound / m2

           emit_co2(i,j)  = emif_co2(n) * tot_emit_dom
           emit_co(i,j)   = emif_co(n)  * tot_emit_dom
           emit_ch4(i,j)  = emif_ch4(n) * tot_emit_dom
           emit_nmhc(i,j) = emif_nmhc(n) * tot_emit_dom
           emit_h2(i,j)   = emif_h2(n) * tot_emit_dom
           emit_nox(i,j)  = emif_nox(n) * tot_emit_dom
           emit_n2o(i,j)  = emif_n2o(n) * tot_emit_dom
           emit_pm25(i,j) = emif_pm25(n) * tot_emit_dom
           emit_tpm(i,j)  = emif_tpm(n) * tot_emit_dom
           emit_tc(i,j)   = emif_tc(n) * tot_emit_dom
           emit_oc(i,j)   = emif_oc(n) * tot_emit_dom
           emit_bc(i,j)   = emif_bc(n) * tot_emit_dom

!         Output the burned area per PFT
          burnvegf(i,j)=100.*burnveg(i,j)/grclarea(i)

          endif
c
         endif
530     continue
520   continue
c
600   continue

      return
      end

