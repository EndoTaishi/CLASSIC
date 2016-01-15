module disturbance_scheme

! Central module for all disturbance scheme-related operations

! J. Melton. Mar 26, 2014

implicit none

! Subroutines contained in this module:
public  :: disturb
public  :: burntobare

contains

! ------------------------------------------------------------------

subroutine disturb (stemmass, rootmass, gleafmas, bleafmas, &
                            thliq,   wiltsm,  fieldsm,    uwind, &
                            vwind,  lightng,  fcancmx, litrmass, &    
                         prbfrhuc, rmatctem, extnprob, popdon,   &
                              il1,      il2,     sort, nol2pfts, &
                         grclarea,    thice,   popdin, lucemcom, &
                           dofire,   currlat,   iday,  fsnow,    &
!     ------------------ inputs above this line ----------------------          
                         stemltdt, rootltdt, glfltrdt, blfltrdt, &
                         glcaemls, rtcaemls, stcaemls, & 
                         blcaemls, ltrcemls, burnfrac, probfire, &
                         emit_co2, emit_co,  emit_ch4, emit_nmhc, &
                         emit_h2,  emit_nox, emit_n2o, emit_pm25, &
                         emit_tpm, emit_tc,  emit_oc,  emit_bc, &
                         burnvegf, bterm,    mterm,    lterm, &
                         pstemmass, pgleafmass )  
!     ------------------outputs above this line ----------------------
!
!               Canadian Terrestrial Ecosystem Model (CTEM)
!                           Disturbance Subroutine
!
!     11  Jun 2015  - Clean up to make all calculations PFT dependent and add in a
!     V. Arora        capability to do a fire count method for comparison with
!                     our usual parameterization
!
!     26  Mar 2014  - Split subroutine into two and create module. Move all fcancmx
!     J. Melton       adjustments into subroutine burntobare and call from competition
!
!     20  Feb 2014  - Adapt to deal with competition on. Bring in code that makes
!     J. Melton       bare fractions from competition module. Moved parameters to
!                     ctem_params.f90
!
!     4   Jan 2014  - Convert to f90 and include a saturation effect for lightning
!     J. Melton       strikes

!     25  Jul 2013  - Add in module for common params, cleaned up code 
!     J. Melton

!     24  Sep 2012  - add in checks to prevent calculation of non-present
!     J. Melton       pfts

!     09  May 2012  - addition of emission factors and revising of the
!     J. Melton       fire scheme

!     15  May 2003  - this subroutine calculates the litter generated
!     V. Arora        and c emissions from leaves, stem, and root 
!                     components due to fire. c emissions from burned
!                     litter are also estimated. at present no other 
!                     form of disturbance is modelled.
!     inputs 

!     stemmass  - stem mass for each of the 9 ctem pfts, kg c/m2
!     rootmass  - root mass for each of the 9 ctem pfts, kg c/m2
!     gleafmas  - green leaf mass for each of the 9 ctem pfts, kg c/m2
!     bleafmas  - brown leaf mass 
!     thliq     - liquid soil moisture content
!     wiltsm    - wilting point soil moisture content
!     fieldsm   - field capacity soil moisture content
!     uwind     - wind speed, m/s
!     vwind     - wind speed, m/s
!     lightng   - total lightning, flashes/(km^2.year)
!                 it is assumed that cloud to ground lightning is
!                 some fixed fraction of total lightning.
!     fcancmx   - fractional coverages of ctem's 9 pfts
!     litrmass  - litter mass for each of the 9 pfts
!     prbfrhuc  - probability of fire due to human causes
!     rmatctem  - fraction of roots in each soil layer for each pft
!     extnprob  - fire extinguishing probability
!     pftareab  - areas of different pfts in a grid cell, before fire, km^2
!     il1,il2   - il1=1, il2=ilg
!     sort      - index for correspondence between 9 pfts and size 12 of
!                 parameters vectors
!     nol2pfts  - number of level 2 ctem pfts
!     grclarea  - gcm grid cell area, km^2
!     thice     - frozen soil moisture content over canopy fraction
!     popdin    - population density (people / km^2)
!     lucemcom  - land use change (luc) related combustion emission losses,
!                 u-mol co2/m2.sec 
!     dofire    - boolean, if true allow fire, if false no fire.
!     fsnow    - fraction of snow simulated by class

!     Outputs

!     stemltdt  - stem litter generated due to disturbance (kg c/m2)
!     rootltdt  - root litter generated due to disturbance (kg c/m2)
!     glfltrdt  - green leaf litter generated due to disturbance (kg c/m2)
!     blfltrdt  - brown leaf litter generated due to disturbance (kg c/m2)
!     burnarea  - total area burned, km^2
!     burnfrac  - total areal fraction burned, (%)
!     probfire  - probability of fire

!     note the following c burned will be converted to a trace gas 
!     emission or aerosol on the basis of emission factors.

!     glcaemls  - green leaf carbon emission losses, kg c/m2
!     blcaemls  - brown leaf carbon emission losses, kg c/m2
!     rtcaemls  - root carbon emission losses, kg c/m2
!     stcaemls  - stem carbon emission losses, kg c/m2
!     ltrcemls  - litter carbon emission losses, kg c/m2

!     emission factors for trace gases and aerosols. units are
!     g of compound emitted per kg of dry organic matter.
!     values are taken from li et al. 2012 biogeosci 
!     emif_co2  - carbon dioxide
!     emif_co   - carbon monoxide
!     emif_ch4  - methane
!     emif_nmhc - non-methane hydrocarbons
!     emif_h2   - hydrogen gas
!     emif_nox  - nitrogen oxides
!     emif_n2o  - nitrous oxide
!     emif_pm25 - particulate matter less than 2.5 um in diameter
!     emif_tpm  - total particulate matter
!     emif_tc   - total carbon
!     emif_oc   - organic carbon
!     emif_bc   - black carbon

!     emitted compounds from biomass burning in g of compound
!     emit_co2  - carbon dioxide
!     emit_co   - carbon monoxide
!     emit_ch4  - methane
!     emit_nmhc - non-methane hydrocarbons
!     emit_h2   - hydrogen gas
!     emit_nox  - nitrogen oxides
!     emit_n2o  - nitrous oxide
!     emit_pm25 - particulate matter less than 2.5 um in diameter
!     emit_tpm  - total particulate matter
!     emit_tc   - total carbon
!     emit_oc   - organic carbon
!     emit_bc   - black carbon

!     tot_emit  - sum of all pools to be converted to emissions/aerosols (g c/m2)
!     tot_emit_dom - tot_emit converted to kg dom / m2

!     hb_interm - interm calculation
!     hbratio   - head to back ratio of ellipse
!     burnvegf - total per PFT areal fraction burned

use ctem_params, only : ignd, icc, ilg, ican, zero,kk, pi, c2dom, crop, &
                        iccp1, standreplace, tolrance, bmasthrs_fire, &
                        lwrlthrs, hgrlthrs, parmlght, parblght, reparea, popdthrshld, & 
                        f0, maxsprd, frco2glf, frco2blf, &
                        frltrglf, frltrblf, frco2stm, frltrstm, frco2rt, frltrrt, &
                        frltrbrn, emif_co2, emif_co, emif_ch4, emif_nmhc, emif_h2, &
                        emif_nox, emif_n2o, emif_pm25, emif_tpm, emif_tc, emif_oc, emif_bc, &
                        grass, extnmois_veg, extnmois_duff

implicit none


real, dimension(ilg,icc), intent(out) :: pstemmass 
real, dimension(ilg,icc), intent(out) :: pgleafmass

logical, intent(in) :: popdon   ! if set true use population density data to calculate fire extinguishing 
                                            ! probability and probability of fire due to human causes, 
                                            ! or if false, read directly from .ctm file

integer :: il1,il2,i,j,k,m,k1,k2,n

integer :: sort(icc), nol2pfts(ican), iday

logical :: dofire, fire(ilg)

real :: stemmass(ilg,icc), rootmass(ilg,icc), gleafmas(ilg,icc),      &
           bleafmas(ilg,icc),     thliq(ilg,ignd),    wiltsm(ilg,ignd),  &
             fieldsm(ilg,ignd),        uwind(ilg),        vwind(ilg),    &
            fcancmx(ilg,icc),      lightng(ilg),litrmass(ilg,icc+1),     &
               prbfrhuc(ilg),     extnprob(ilg), &
        rmatctem(ilg,icc,ignd),     thice(ilg,ignd),    popdin(ilg),  &
               lucemcom(ilg),    tmpprob(ilg),  currlat(ilg), fsnow(ilg) 

real ::  stemltdt(ilg,icc), rootltdt(ilg,icc), glfltrdt(ilg,icc), &
               burnarea(ilg), glcaemls(ilg,icc),    &
           rtcaemls(ilg,icc), stcaemls(ilg,icc), ltrcemls(ilg,icc),    &
           blfltrdt(ilg,icc), blcaemls(ilg,icc),     burnfrac(ilg),    &
          emit_co2(ilg,icc),   emit_co(ilg,icc), emit_ch4(ilg,icc),    &
         emit_nmhc(ilg,icc),   emit_h2(ilg,icc), emit_nox(ilg,icc),    &
          emit_n2o(ilg,icc), emit_pm25(ilg,icc), emit_tpm(ilg,icc),    &
           emit_tc(ilg,icc),   emit_oc(ilg,icc),  emit_bc(ilg,icc)    

real ::  biomass(ilg,icc),        bterm(ilg), drgtstrs(ilg,icc), &
            betadrgt(ilg,ignd),     avgdryns(ilg),        fcsum(ilg),  &
               avgbmass(ilg),        mterm(ilg),     c2glgtng(ilg),    &
               betalght(ilg),            y(ilg),        lterm(ilg),    &
               probfire(ilg),             ctime,            random,    &
                        temp,     betmsprd(ilg),       smfunc(ilg),    &
                   wind(ilg),      wndfunc(ilg),     sprdrate(ilg),    &
                lbratio(ilg),     arbn1day(ilg),     areamult(ilg),    &
                 vegarea(ilg),     grclarea(ilg),    &
                     tot_emit,      tot_emit_dom,     burnvegf(ilg,icc)   

real :: natural_ignitions(ilg),  anthro_ignitions(ilg), fire_supp(ilg), &
         num_ignitions(ilg, icc), num_fires(ilg, icc)

real :: hb_interm, hbratio(ilg)
real, dimension(ilg) :: surface_duff_f  ! fraction of biomass that is in the surface duff (grass brown leaves + litter) 
real, dimension(ilg,icc) :: pftareab    ! pft area before fire (km2)

real :: ymin, ymax, slope
real :: soilterm, duffterm              ! temporary variables


real, dimension(ilg,icc) :: bterm_veg     
real, dimension(ilg,icc) :: duff_frac_veg
real, dimension(ilg,icc) :: mterm_veg
real, dimension(ilg,icc) :: probfire_veg
real, dimension(ilg,icc) :: smfunc_veg
real, dimension(ilg,icc) :: sprdrate_veg
real, dimension(ilg,icc) :: arbn1day_veg
real, dimension(ilg,icc) :: burnarea_veg
logical, dimension(ilg,icc) :: fire_veg
real :: soilterm_veg, duffterm_veg, betmsprd_veg, betmsprd_duff      ! temporary variables

!     ------------------------------------------------------------------
!     Constants and parameters are located in ctem_params.f90
!     -----------------------------------------------------------------

!     * if icc /= 9 or ignd /= 3 this subroutine will need changes.
      IF(ICC.NE.12)      CALL XIT('DISTURB',-1)    !YW January 15, 2016 
      IF(IGND.NE.10)     CALL XIT('DISTURB',-2)   !YW April 14, 2015 

!     initialize required arrays to zero, or assign value

      do 140 j = 1,icc
        do 150 i = il1, il2
          stemltdt(i,j)=0.0     !stem litter due to disturbance
          rootltdt(i,j)=0.0     !root litter due to disturbance
          glfltrdt(i,j)=0.0     !green leaf litter due to disturbance
          blfltrdt(i,j)=0.0     !brown leaf litter due to disturbance
          biomass(i,j)=0.0      !total biomass for fire purposes
          drgtstrs(i,j)=0.0     !soil dryness factor for pfts
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

          bterm_veg(i,j)=0.0            !biomass fire probability term, Vivek
          mterm_veg(i,j)=0.0            !moisture fire probability term, Vivek
          probfire_veg(i,j)=0.0         !PFT fire probability term, Vivek
          fire_veg(i,j)=.false.         !fire occuring logical, Vivek
          duff_frac_veg(i,j)=0.0        !duff fraction for each PFT, Vivek
          smfunc_veg(i,j)=0.0           !soil moisture dependence on fire spread rate, Vivek
          sprdrate_veg(i,j)=0.0         !per PFT fire spread rate
          arbn1day_veg(i,j)=0.0         !per PFT area burned in 1 day
          burnarea_veg(i,j)=0.0         !per PFT area burned over fire duration
150     continue                  
140   continue

      do 160 k = 1,ignd
        do 170 i = il1, il2
          betadrgt(i,k)=1.0     !dryness term for soil layers
170     continue                  
160   continue

      do 180 i = il1, il2
        
        avgbmass(i)=0.0         !avg. veg. biomass over the veg. fraction of grid cell
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
        surface_duff_f(i)=0.0   

180   continue

!     if not simulating fire, leave the subroutine now.
      if (.not. dofire) goto 600


      do 190 i = il1, il2
        if(extnprob(i).le.zero) then
          write(6,*)'fire extinguishing prob. (',i,'= ',extnprob(i)
          write(6,*)'please use an appropriate value of this paramater'
          write(6,*)'else the whole grid cell will burn down leading to'
          write(6,*)'numerical problems.'
          call xit('disturb',-2)
        endif
190   continue

!     initialization ends    

!     Find pft areas before
        do 82 j = 1, icc
          do  83 i = il1, il2
            pftareab(i,j)=fcancmx(i,j)*grclarea(i)  ! area in km^2
83        continue
82      continue

!     ------------------------------------------------------------------

!     Find the probability of fire as a product of three functions
!     with dependence on total biomass, soil moisture, and lightning 

!     1. Dependence on total biomass

      do 200 j = 1, icc
        do 210 i = il1, il2

         if (.not. crop(j)) then !don't allow it to bring in crops since they are not allowed to burn. 

!          Root biomass is not used to initiate fire. For example if
!          the last fire burned all grass leaves, and some of the roots
!          were left, its unlikely these roots could catch fire. 
           biomass(i,j)=gleafmas(i,j)+bleafmas(i,j)+stemmass(i,j)+ &
                      litrmass(i,j)

!          Find average biomass over the vegetated fraction
           avgbmass(i) = avgbmass(i)+biomass(i,j)*fcancmx(i,j)

!          Sum up the vegetated area
           fcsum(i)=fcsum(i) + fcancmx(i,j)

         endif

210     continue
200   continue

!     calculate bterm for individual PFTs as well

      do 251 j=1, icc
       do 252 i = il1, il2
          bterm_veg(i,j)=min(1.0,max(0.0,(biomass(i,j)-bmasthrs_fire(1))/(bmasthrs_fire(2)-bmasthrs_fire(1))))     
252    continue
251   continue


!     2. Dependence on soil moisture

!     This is calculated in a way such that the more dry the root zone
!     of a pft type is, and more fractional area is covered with that
!     pft, the more likely it is that fire will get started. that is
!     the dryness factor is weighted by fraction of roots in soil
!     layers, as well as according to the fractional coverage of 
!     different pfts. the assumption here is that if there is less 
!     moisture in root zone, then it is more likely the vegetation 
!     will be dry and thus the likeliness of fire is more.

!     First find the dryness factor for each soil layer.

!     If there is snow on the ground, do not allow fire so set betadrgt to
!     0 for all soil layers otherwise calculate as per normal.  
      do i = il1, il2
        if (fsnow(i) .eq. 0.) then
          do j = 1, ignd
           betadrgt(i,j)=min(1.0,max(0.0,(thliq(i,j)+thice(i,j)-wiltsm(i,j))/(fieldsm(i,j)-wiltsm(i,j))))   
          end do        
        end if
      end do

!     Now find weighted value of this dryness factor averaged over 
!     the rooting depth, for each pft

      do 320 j = 1, icc
        do 330 i = il1, il2
         if (.not. crop(j)) then
     
          drgtstrs(i,j) =  (betadrgt(i,1))*rmatctem(i,j,1) + (betadrgt(i,2))*rmatctem(i,j,2) + &
                         (betadrgt(i,3))*rmatctem(i,j,3)

          drgtstrs(i,j) = min(1.0,max(0.0,drgtstrs(i,j)/(rmatctem(i,j,1)+rmatctem(i,j,2)+rmatctem(i,j,3))))

! !         Next find this dryness factor averaged over the vegetated fraction
!           avgdryns(i) = avgdryns(i) + drgtstrs(i,j)*fcancmx(i,j)
! 
          ! The litter and brown leaves are not affected by the soil water potential
          ! therefore they will react only to the moisture conditions (our proxy here
          ! is the upper soil moisture). If it is dry they increase the probability of 
          ! fire corresponding to the proportion of total C they contribute. Only allow
          ! if there is no snow. 
!           if (biomass(i,j) .gt. 0. .and. fsnow(i) .eq. 0.) then
!             ! The surface duff calculation ignores the litter on the bare fraction.
!             surface_duff_f(i) = surface_duff_f(i) + (bleafmas(i,j)+litrmass(i,j)) &
!                                                      /biomass(i,j) * fcancmx(i,j)
! 
!           end if

         endif
       write(6,6990) j, drgtstrs(1,j), avgdryns, surface_duff_f
6990   format(i3, 10f7.3)
330     continue
320   continue


!     Use average root zone vegetation dryness to find likelihood of
!     fire due to moisture. 

!       calculate mterm for each PFT

      do 381 j = 1, icc
        do 382 i = il1, il2
           ! duff fraction for each PFT, Vivek
            if (biomass(i,j) .gt. 0. .and. fsnow(i) .eq. 0.) then
               duff_frac_veg(i,j) = (bleafmas(i,j)+litrmass(i,j)) / biomass(i,j)
            end if

           ! drgtstrs(i,j) is phi_root in Melton and Arora GMDD (2015) paper
           soilterm_veg = 1.0-tanh((1.75*drgtstrs(i,j)/extnmois_veg)**2) 
           duffterm_veg = 1.0-tanh((1.75*betadrgt(i,1)/extnmois_duff)**2)

           if(fcancmx(i,j) .gt. zero)then
             mterm_veg(i,j)=soilterm_veg*(1.-duff_frac_veg(i,j)) + duffterm_veg*duff_frac_veg(i,j)
           else
             mterm_veg(i,j)=0.0   !no fire likelihood due to moisture if no vegetation
           endif

382     continue
381   continue


!     3. dependence on lightning

!     Dependence on lightning is modelled in a simple way which implies that
!     a large no. of lightning flashes are more likely to cause fire than
!     few lightning flashes.

      do 400 i = il1, il2

!       New approximation of Price and Rind equation. It was developed from a more complete dataset than Prentice and Mackerras. 
!       Lightning comes in in units of flashes/km2/yr so divide by 12 to make per month.
        !c2glgtng(i)=lightng(i) !FLAG FireMIP lightning is already C2G! If using 'normal' lightning climatology use line below.
        c2glgtng(i)=0.219913*exp(0.0058899*abs(currlat(i)))*lightng(i)

        betalght(i)=min(1.0,max(0.0,(c2glgtng(i)-lwrlthrs)/(hgrlthrs-lwrlthrs)))
        y(i)=1.0/( 1.0+exp((parmlght-betalght(i))/parblght) )

!       No need to calculate each time, once settled on parameters, precalc and moved into a parameter. JM. Feb 19 2014.
        ymin=1.0/( 1.0+exp((parmlght-0.0)/parblght) )
        ymax=1.0/( 1.0+exp((parmlght-1.0)/parblght) )
        slope=abs(0.0-ymin)+abs(1.0-ymax)

        temp=y(i)+(0.0-ymin)+betalght(i)*slope

!       Determine the probability of fire due to human causes
!       this is based upon the population density from the .popd read-in file
        if (popdon) then
            prbfrhuc(i)=min(1.0,(popdin(i)/popdthrshld)**0.43) !From Kloster et al. (2010)
        end if

        lterm(i)=max(0.0, min(1.0, temp+(1.0-temp)*prbfrhuc(i) ))

!       ----------------------- Number of fire calculations ----------------------\\

!       This is not used in CTEM in general but we keep the code here for testing purposes

!       calculate natural and anthorpogenic ignitions/km2.day
!       the constant 0.25 assumes not all c2g lightning hits causes ignitions, only 0.25 of them do
!       the constant (1/30.4) converts c2g lightning from flashes/km2.month to flashes/km2.day
!       MAKE SURE LIGHTNING IS IN UNITS OF FLASHES/KM2.MONTH

!       Eqs. (4) and (5) of Li et al. 2012 doi:10.5194/bg-9-2761-2012 + also see corrigendum

        !natural_ignitions(i)=c2glgtng(i) * 0.25 * (1/30.4)
        !anthro_ignitions(i)=9.72E-4 * (1/30.4) * popdin(i) * (6.8 * (popdin(i)**(-0.6)) )

!       calculate fire suppression also as a function of population density.
!       Li et al. (2012) formulation is quite similar to what we already use based
!       on Kloster et al. (2010, I think) but we use Kloster's formulation together
!       with our fire extingishing probability. Here, the fire suppression formulation
!       is just by itself

        !fire_supp(i)=0.99 - ( 0.98*exp(-0.025*popdin(i)*(i)) )

!       ----------------------- Number of fire calculations ----------------------//
400   continue

!     calculate fire probability for each PFT. Recall that lightning term is still grid
!     averaged

      do 421 j = 1, icc
        do 422 i = il1, il2
           probfire_veg(i,j)=bterm_veg(i,j)*mterm_veg(i,j)*lterm(i)
           if (probfire_veg(i,j) .gt. zero) fire_veg(i,j)=.true.

!          ----------------------- Number of fire calculations ----------------------\\

!       This is not used in CTEM in general but we keep the code here for testing purposes
!
!          calculate total number of ignitions based natural and anthorpogenic ignitions
!          for the whole grid cell
    
           !num_ignitions(i,j) = ( natural_ignitions(i) + anthro_ignitions(i) ) * fcancmx(i,j)*grclarea(i)

!          finally calculate number of fire, noting that not all ignitions turn into fire
!          because moisture and biomass may now allow that to happen, and some of those
!          will be suppressed due to fire fighting efforts

           !num_fires(i,j) = num_ignitions(i,j)*(1-fire_supp(i))*bterm_veg(i,j)*mterm_veg(i,j)
!
!          ----------------------- Number of fire calculations ----------------------//

422     continue
421   continue

!     Calculate area burned for each PFT, make sure it's not greater than the
!     area available, then find area and fraction burned for the whole gridcell

      do 500 j = 1, icc
        do 501 i = il1, il2
          if (fire_veg(i,j) ) then

!         soil moisture dependence on fire spread rate
          betmsprd_veg = (1. - min(1., (drgtstrs(i,j)/extnmois_veg) ))**2
          betmsprd_duff = (1. - min(1., (betadrgt(i,1)/extnmois_duff) ))**2
          smfunc_veg(i,j)= betmsprd_veg*(1-duff_frac_veg(i,j)) + betmsprd_duff*duff_frac_veg(i,j)

!         wind speed, which is gridcell specific
          wind(i)=sqrt(uwind(i)**2.0 + vwind(i)**2.0)
          wind(i)=wind(i)*3.60     ! change m/s to km/hr

!         Length to breadth ratio from Li et al. (2012)
          lbratio(i)=1.0+10.0*(1.0-exp(-0.06*wind(i)))

!         head to back ratio from Li et al. (2012). 
          hb_interm = (lbratio(i)**2 - 1.0)**0.5
          hbratio(i) = (lbratio(i) + hb_interm)/(lbratio(i) - hb_interm)

!         dependence of spread rate on wind
          wndfunc(i)= (2.0 * lbratio(i)) / (1.0 + 1.0 / hbratio(i)) * f0 

!         fire spread rate per PFT
          n = sort(j)  
          sprdrate_veg(i,j)= maxsprd(n) * smfunc_veg(i,j) * wndfunc(i) 

!         area burned in 1 day for that PFT
          arbn1day_veg(i,j)=(pi*24.0*24.0*sprdrate_veg(i,j)**2)/(4.0 * lbratio(i))*(1.0 + 1.0 / hbratio(i))**2

!         fire extinguishing probability as a function of grid-cell averaged population density
          if (popdon) then
            extnprob(i)=max(0.0,0.9-exp(-0.025*popdin(i)))  !FLAG changed from -0.015 to keep in line with Kloster. JM Jun 24 2015
            extnprob(i)=0.5+extnprob(i)/2.0
          end if

!         area multipler to calculate area burned over the duration of the fire
          areamult(i)=((1.0-extnprob(i))*(2.0-extnprob(i)))/ extnprob(i)**2                              

!         per PFT area burned, km^2
          burnarea_veg(i,j)=arbn1day_veg(i,j)*areamult(i)*(grclarea(i)*fcancmx(i,j)*probfire_veg(i,j))/reparea
          
!          ------- Area burned based on number of fire calculations ----------------------\\

!       This is not used in CTEM in general but we keep the code here for testing purposes

!         the constant 4 is suppose to address the fact that Li et al. (2012) suggested to
!         double the fire spread rates. However, if we do that than our usual calculations
!         based on CTEM's original parameterization will be affected. Rather than do that
!         we just use a multiplier of 4, since doubling fire spread rates means 4 times the
!         area burned
!
!          burnarea_veg(i,j)=arbn1day_veg(i,j)*num_fires(i,j)*2.0  !flag test was 4!
 !
!          ------- Area burned based on number of fire calculations ----------------------//

!         if area burned greater than area of PFT, set it to area of PFT

          if ( burnarea_veg(i,j) .gt. grclarea(i)*fcancmx(i,j) ) then
             burnarea_veg(i,j)=grclarea(i)*fcancmx(i,j)
          endif

         endif
501     continue
500   continue

!     Calculate gridcell area burned and fraction

      do 510 j = 1, icc
        do 511 i = il1, il2
          burnarea(i)=burnarea(i)+burnarea_veg(i,j)
511     continue       
510   continue

      do 512 i = il1, il2
        burnfrac(i)=burnarea(i)/grclarea(i)
512   continue       

!     Finally estimate amount of litter generated from each pft, and
!     each vegetation component (leaves, stem, and root) based on their
!     resistance to combustion. Update the veg pools due to combustion.

      do 520 j = 1, icc
       n = sort(j)
        do 530 i = il1, il2
          if(pftareab(i,j) .gt. zero)then

            !Set aside these pre-disturbance stem and root masses for use
            !in burntobare subroutine.
            pstemmass(i,j)=stemmass(i,j)
            pgleafmass(i,j)=gleafmas(i,j)

            glfltrdt(i,j)= frltrglf(n) *gleafmas(i,j) *(burnarea_veg(i,j) /pftareab(i,j))
            blfltrdt(i,j)= frltrblf(n) *bleafmas(i,j) *(burnarea_veg(i,j) /pftareab(i,j))
            stemltdt(i,j)= frltrstm(n) *stemmass(i,j) *(burnarea_veg(i,j) /pftareab(i,j))
            rootltdt(i,j)= frltrrt(n)  *rootmass(i,j) *(burnarea_veg(i,j) /pftareab(i,j))
            glcaemls(i,j)= frco2glf(n) *gleafmas(i,j) *(burnarea_veg(i,j) /pftareab(i,j))
            blcaemls(i,j)= frco2blf(n) *bleafmas(i,j) *(burnarea_veg(i,j) /pftareab(i,j))
            stcaemls(i,j)= frco2stm(n) *stemmass(i,j) *(burnarea_veg(i,j) /pftareab(i,j))
            rtcaemls(i,j)= frco2rt(n)  *rootmass(i,j) *(burnarea_veg(i,j) /pftareab(i,j))
            ltrcemls(i,j)= frltrbrn(n) *litrmass(i,j) *(burnarea_veg(i,j) /pftareab(i,j))

!           Update the pools:
            gleafmas(i,j)=gleafmas(i,j) - glfltrdt(i,j) - glcaemls(i,j)
            bleafmas(i,j)=bleafmas(i,j) - blfltrdt(i,j) - blcaemls(i,j)
            stemmass(i,j)=stemmass(i,j) - stemltdt(i,j) - stcaemls(i,j)
            rootmass(i,j)=rootmass(i,j) - rootltdt(i,j) - rtcaemls(i,j)
            litrmass(i,j)=litrmass(i,j) + glfltrdt(i,j) + blfltrdt(i,j) + stemltdt(i,j) + rootltdt(i,j) - ltrcemls(i,j)

!           Output the burned area per PFT
            burnvegf(i,j)=burnarea_veg(i,j)/grclarea(i)

          endif
530     continue
520   continue


600   continue ! If .not. dofire then we enter here and perform the calculations for the emissions
               ! since we might have some from luc. 

!     We also estimate CO2 emissions from each
!     of these components. Note that the litter which is generated due 
!     to disturbance is uniformly distributed over the entire area of 
!     a given pft, and this essentially thins the vegetation biomass. 
!     If compete is not on, this does not change the vegetation fractions,
!     if competition is on a fraction will become bare. That is handled in
!     burntobare subroutine called from competition subroutine.

      do 620 j = 1, icc
       n = sort(j)
        do 630 i = il1, il2

!          Calculate the emissions of trace gases and aerosols based upon how
!          much plant matter was burnt

!          Sum all pools that will be converted to emissions/aerosols (g c/m2)
           tot_emit = (glcaemls(i,j) + blcaemls(i,j) + rtcaemls(i,j)+ stcaemls(i,j) + ltrcemls(i,j)) * 1000.0

!          Add in the emissions due to luc fires (deforestation)
!          the luc emissions are converted from umol co2 m-2 s-1 to g c m-2 (day-1) before adding to tot_emit         
           tot_emit = tot_emit + (lucemcom(i) / 963.62 * 1000.0)

!          Convert burnt plant matter from carbon to dry organic matter using 
!          a conversion factor, assume all parts of the plant has the same
!          ratio of carbon to dry organic matter. units: kg dom / m2
           tot_emit_dom = tot_emit / c2dom

!          Convert the dom to emissions/aerosols using emissions factors
!          units: g species / m2

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

630     continue
620   continue

end subroutine disturb

! ------------------------------------------------------------------------------------

subroutine burntobare(il1, il2, nilg, sort,pvgbioms,pgavltms,pgavscms,fcancmx, burnvegf, stemmass, &
                      rootmass, gleafmas, bleafmas, litrmass, soilcmas, pstemmass, pgleafmass,&
                      nppveg)

!     Update fractional coverages of pfts to take into account the area
!     burnt by fire. Adjust all pools with new densities in their new
!     areas and increase bare fraction.

!     And while we are doing this also run a small check to make sure
!     grid averaged quantities do not get messed up.

!     J. Melton. Mar 26 2014  - Create subroutine


use ctem_params, only : crop, icc, seed, standreplace, grass, zero, &
                        iccp1, tolrance, numcrops

implicit none

integer, intent(in) :: il1
integer, intent(in) :: il2
integer, intent(in) :: nilg       ! no. of grid cells in latitude circle (this is passed in as either ilg or nlat depending on mos/comp)
integer, dimension(icc), intent(in) :: sort             ! index for correspondence between 9 ctem pfts and
                                                        ! size 12 of parameter vectors
real, dimension(nilg), intent(in) :: pvgbioms           ! initial veg biomass
real, dimension(nilg), intent(in) :: pgavltms           ! initial litter mass
real, dimension(nilg), intent(in) :: pgavscms           ! initial soil c mass
real, dimension(nilg,icc), intent(inout) :: fcancmx     ! initial fractions of the ctem pfts
real, dimension(nilg,icc), intent(in) :: burnvegf       ! total per PFT areal fraction burned
real, dimension(nilg,icc), intent(inout) :: gleafmas    ! green leaf carbon mass for each of the 9 ctem pfts, kg c/m2
real, dimension(nilg,icc), intent(inout) :: bleafmas    ! brown leaf carbon mass for each of the 9 ctem pfts, kg c/m2
real, dimension(nilg,icc), intent(inout) :: stemmass    ! stem carbon mass for each of the 9 ctem pfts, kg c/m2
real, dimension(nilg,icc), intent(inout) :: rootmass    ! roots carbon mass for each of the 9 ctem pfts, kg c/m2
real, dimension(nilg,icc), intent(inout) :: nppveg      ! npp for individual pfts,  u-mol co2/m2.sec
real, dimension(nilg,iccp1), intent(inout) :: soilcmas  ! soil carbon mass for each of the 9 ctem pfts + bare, kg c/m2
real, dimension(nilg,iccp1), intent(inout) :: litrmass  ! litter carbon mass for each of the 9 ctem pfts + bare, kg c/m2
real, dimension(nilg,icc), intent(in)    :: pstemmass   ! grid averaged stemmass prior to disturbance, kg c/m2
real, dimension(nilg,icc), intent(in)    :: pgleafmass  ! grid averaged rootmass prior to disturbance, kg c/m2

logical, dimension(nilg) :: shifts_occur      ! true if any fractions changed
integer :: i, j, n, k
real :: pftfraca_old
real :: term                                  ! temp variable for change in fraction due to fire
real, dimension(nilg) :: pbarefra             ! bare fraction prior to fire              
real, dimension(nilg) :: barefrac             ! bare fraction of grid cell
real, dimension(nilg) :: litr_lost            ! litter that is transferred to bare 
real, dimension(nilg) :: soilc_lost           ! soilc that is transferred to bare
real, dimension(nilg) :: vgbiomas_temp        ! grid averaged vegetation biomass for internal checks, kg c/m2
real, dimension(nilg) :: gavgltms_temp        ! grid averaged litter mass for internal checks, kg c/m2
real, dimension(nilg) :: gavgscms_temp        ! grid averaged soil c mass for internal checks, kg c/m2
real, dimension(nilg,icc) :: pftfracb         ! pft fractions before accounting for creation of bare ground
real, dimension(nilg,icc) :: pftfraca         ! pft fractions after accounting for creation of bare ground
real :: frac_chang                            ! pftfracb - pftfraca

! -----------------------------------------

! Do some initializations
do 10 i = il1, il2
        shifts_occur(i) = .false.
        pbarefra(i)=1.0
        barefrac(i)=1.0
        vgbiomas_temp(i)=0.0
        gavgltms_temp(i)=0.0
        gavgscms_temp(i)=0.0
        litr_lost(i)=0.0
        soilc_lost(i)=0.0
10  continue

!  Account for disturbance creation of bare ground. This occurs with relatively low
!  frequency and is PFT dependent. We must adjust the amount of bare ground created
!  to ensure that we do not increase the density of the remaining vegetation. 

       do 20 i = il1, il2
          do 25 j = 1, icc
            if(.not. crop(j))then  

              n = sort(j)

              pbarefra(i)=pbarefra(i)-fcancmx(i,j)
              pftfracb(i,j)=fcancmx(i,j)

              pftfraca(i,j) = max(seed,fcancmx(i,j) - burnvegf(i,j) * standreplace(n))
        
              fcancmx(i,j) = pftfraca(i,j)
   
              barefrac(i)=barefrac(i)-fcancmx(i,j)

            else  !crops

              pbarefra(i)=pbarefra(i)-fcancmx(i,j)
              barefrac(i)=barefrac(i)-fcancmx(i,j)

            endif

25      continue
20   continue

      do 40 i = il1, il2
       do 50 j = 1, icc
        if(.not. crop(j))then 
 
          ! Test the pftfraca to ensure it does not cause densification of the exisiting biomass  
          ! Trees compare the stemmass while grass compares the root mass. 
          if (pftfraca(i,j) .ne. pftfracb(i,j)) then

            shifts_occur(i) = .true.

            term = pftfracb(i,j)/pftfraca(i,j)

            if (.not. grass(j)) then
               if (stemmass(i,j)*term .gt. pstemmass(i,j) .and. pstemmass(i,j) .gt. 0.) then  !the pstemmass is from before the fire occurred, i.e. no thinning!

                 pftfraca_old = pftfraca(i,j)
                 pftfraca(i,j) = max(seed,stemmass(i,j) * pftfracb(i,j) / pstemmass(i,j))
                 fcancmx(i,j) = pftfraca(i,j)

                 ! adjust the bare frac to accomodate for the changes 
                 barefrac(i) = barefrac(i) + pftfraca_old - pftfraca(i,j)

               end if
            else !grasses

               if (gleafmas(i,j)*term .gt. pgleafmass(i,j) .and. pgleafmass(i,j) .gt. 0.) then !the pgleafmass is from before the fire occurred, i.e. no thinning!

                 pftfraca_old = pftfraca(i,j)
                 pftfraca(i,j) = max(seed,gleafmas(i,j) * pftfracb(i,j) / pgleafmass(i,j))
                 fcancmx(i,j) = pftfraca(i,j)

                 ! adjust the bare frac to accomodate for the changes 
                 barefrac(i) = barefrac(i) + pftfraca_old - pftfraca(i,j)

               end if
            end if

            term = pftfracb(i,j)/pftfraca(i,j)
            gleafmas(i,j)=gleafmas(i,j)*term
            bleafmas(i,j)=bleafmas(i,j)*term
            stemmass(i,j)=stemmass(i,j)*term
            rootmass(i,j)=rootmass(i,j)*term
            nppveg(i,j)  =nppveg(i,j)*term

!           Soil and litter carbon are treated such that we actually transfer the carbon
!           to the bare fraction since it would remain in place as a location was devegetated
!           In doing so we do not adjust the litter or soilc density on the 
!           remaining vegetated fraction. But we do adjust it on the bare fraction to ensure
!           our carbon balance works out.
            frac_chang = pftfracb(i,j) - pftfraca(i,j)
            litr_lost(i)= litr_lost(i) + litrmass(i,j) * frac_chang
            soilc_lost(i)= soilc_lost(i) + soilcmas(i,j) * frac_chang

        ! else  

        !    no changes

          end if

        endif  !crop
50     continue
40  continue

      do 100 i = il1, il2
       if (shifts_occur(i)) then !only do checks if we actually shifted fractions here. 

        if(barefrac(i).ge.zero .and. barefrac(i) .gt. pbarefra(i))then
          litrmass(i,iccp1) = (litrmass(i,iccp1)*pbarefra(i) + litr_lost(i)) / barefrac(i) 
          soilcmas(i,iccp1) = (soilcmas(i,iccp1)*pbarefra(i) + soilc_lost(i)) / barefrac(i) 
        else if (barefrac(i) .lt. 0.) then  
          
          write(6,*)' In burntobare you have negative bare area, which should be impossible...'
          write(6,*)' bare is',barefrac(i),' original was',pbarefra(i)
          call xit('disturb-burntobare',-6)

        endif
       end if
100  continue

!     check if total grid average biomass density is 
!           same before and after adjusting fractions

    do 200 i = il1, il2

      if (shifts_occur(i)) then !only do checks if we actually shifted fractions here. 

       do 250 j = 1, icc

          vgbiomas_temp(i)=vgbiomas_temp(i)+fcancmx(i,j)*(gleafmas(i,j)+&
          bleafmas(i,j)+stemmass(i,j)+rootmass(i,j))
          gavgltms_temp(i)=gavgltms_temp(i)+fcancmx(i,j)*litrmass(i,j)
          gavgscms_temp(i)=gavgscms_temp(i)+fcancmx(i,j)*soilcmas(i,j)

250    continue

        !then add the bare ground in.
        gavgltms_temp(i)=gavgltms_temp(i)+ barefrac(i)*litrmass(i,iccp1)
        gavgscms_temp(i)=gavgscms_temp(i)+ barefrac(i)*soilcmas(i,iccp1)

        if(abs(vgbiomas_temp(i)-pvgbioms(i)).gt.tolrance)then
          write(6,*)'grid averaged biomass densities do not balance'
          write(6,*)'after fractional coverages are changed to take'
          write(6,*)'into account burn area'
          write(6,*)'vgbiomas_temp(',i,')=',vgbiomas_temp(i)
          write(6,*)'pvgbioms(',i,')=',pvgbioms(i)
          call xit('disturb',-7)
        endif
        if(abs(gavgltms_temp(i)-pgavltms(i)).gt.tolrance)then
          write(6,*)'grid averaged litter densities do not balance'
          write(6,*)'after fractional coverages are changed to take'
          write(6,*)'into account burn area'
          write(6,*)'gavgltms_temp(',i,')=',gavgltms_temp(i)
          write(6,*)'pgavltms(',i,')=',pgavltms(i)
          call xit('disturb',-8)
        endif
        if(abs(gavgscms_temp(i)-pgavscms(i)).gt.tolrance)then
          write(6,*)'grid averaged soilc densities do not balance'
          write(6,*)'after fractional coverages are changed to take'
          write(6,*)'into account burn area'
          write(6,*)'gavgscms_temp(',i,')=',gavgscms_temp(i)
          write(6,*)'pgavscms(',i,')=',pgavscms(i)
          call xit('disturb',-9)
        endif

      end if

200  continue

      return

end subroutine burntobare

end module

