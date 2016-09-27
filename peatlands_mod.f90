!>\defgroup moss_photosynthesis

!! Moss photosynthesis subroutine

!>\defgroup peat_soil_het_resp


!>\file

!>Central module for all peatland-related operations

module peatlands_mod


! J. Melton. Sep 26, 2016

implicit none

! Subroutines contained in this module:
public  :: mosspht
public  :: decp

contains

! ------------------------------------------------------------------

!>\ingroup moss_photosynthesis
!!@{

!! Moss photosynthesis subroutine

subroutine mosspht(ilg,ignd,isand,iday,qswnv,thliq,tbar,thpor, &
            &    co2conc,tsurfk,zsnow,delzw,pres,qg,coszs,Cmossmas,dmoss, &
!         output below
            &    anmoss,rmlmoss,cevapms,ievapms,ipeatland &
!    testing
            &    ,iyear, ihour,imin,daylength,pdd)


! History

! J. Melton Sep 26 2016
!   - Bring into rest of model and model formatting, convert to doxygen compatible code

! Created by Yuanqiao Wu, 2015


! ----------

use ctem_params, only : rmlmoss25, tau25m,ektau,gasc,kc25,ko25,ec,ej,eo,evc,sj, &
                        hj,alpha_moss,thpms,thmms

implicit none

! arguments:

integer i               !>
integer j               !>
integer ilg             !>
integer ignd            !>
integer isand(ilg,ignd) !>
integer iday            !>
integer ihour           !>
integer iyear           !>
integer imin            !>
integer ievapms(ilg)    !>
integer ipeatland(ilg)  !>

real:: qswnv(ilg)  !< visible short wave radiation = qswnv in TSOLVE
                    !! and qswnvg in TSOLVC (W/m2)

!    environment related
real     thliq(ilg,ignd),    tbar(ilg,ignd),     thpor(ilg,ignd), &
        bi(ignd),           co2conc(ilg),       &
        zsnow(ilg),         delzw(ilg,ignd),    pres(ilg), &
        fcs(ilg),      fgs(ilg),      fc(ilg),       fg(ilg), &
        qg(ilg),       coszs(ilg)
!    moss related
real:: Cmossmas(ilg)    !< unit kg moss C updated in ctem
real:: mmoss(ilg)       !< unit kg moss mass in kg, update later
real:: dmoss(ilg)       !< unit m, depth of living moss. assume = 2 cm
                        !! can be related to mmoss as a variable
!    ---------------input above this line-------------------------------

real:: anmoss(ilg)      !< net photosynthesis (umol/m2/s)
real:: rmlmoss(ilg)     !< moss autotrophic respiration (umol/m2/s)

!    ---------------output above this line, parameters below------------

! Local parameters:

real, parameter :: tref = 298.16    !< unit K

! Remainder of parameters stored in ctem_params.f90

! -----------------------

integer:: pheno(ilg)    !<phenology flag of mosses, 1 = photosynthesis,
                        !!0 = does not photosynthesis
real:: parm(ilg)        !<par at the ground (moss layer) surface umol/m2/s
real:: tsurf(ilg)       !<grid average ground surface temperature in C
real:: tsurfk(ilg)      !<grid average ground surface temprature in K
real:: wmoss(ilg)       !<water content extraporated from the surface
                        !!humidity qg and thliq of the first soil layer
                        !!unit kg water/ kg dw
real:: wmosmin(ilg)     !<residual water content kg water /kg moss
real:: wmosmax(ilg)     !<maximum water content kg water /kg moss
real:: fwmoss(ilg)      !<relative water content of mosses in g fw /g dw
real:: dsmoss(ilg)      !<degree of moss saturation = relative water
                        !!content/maximum relative water content
real:: cevapms(ilg)     !<evaporation coefficent for moss surface
real:: g_moss(ilg)      !<moss conductance umol/m2/s (based on
                        !!Williams and Flanagan,1998 for Sphagnum)
real:: mwce (ilg)       !<moisture function of dark respiration of moss
real:: tmoss(ilg)       !<moss temperature extraporated from the tbar 1
                        !!and grid averaged ground surface temperature
                        !!tsurf
real:: tmossk(ilg)      !<moss temperature in K
real:: q10rmlmos(ilg)   !<temperature function of the moss dark respiration
real:: gamma(ilg)       !<compensation point for gross photosynthesis (Pa)
real:: o2(ilg)          !<partial presure of oxygen (Pa)
real:: co2a(ilg)        !<partical pressure of co2 (pa) same as in PHTSYN
real:: tau(ilg)         !<arrhenius funciton of temperature
real:: kc(ilg)          !<kinetic coeffficient of CO2 for photosynthesis
real:: ko(ilg)          !<kinetic coeffficient of O2 for photosynthesis
real:: bc(ilg)          !<coefficient used for wc
real:: vcmax25(ilg)     !<seasonal varied maximum carboylation at 25
                        !!sphagnum (fig. 6, Williams and Flanagan, 1998)
real:: vcmax(ilg)       !<max carboxylation rate (umol/m2/s)
real:: jmax25(ilg)      !<maximum electorn transport rate at 25 degrees (umol/m2/s)
real:: jmax(ilg)        !<maximum electorn transport rate (umol/m2/s)
real:: wj(ilg)          !<net co2 assimilation rate limited by
                        !!electron transport (umol/m2/s)=jE in PHTSYN
real:: wc(ilg)          !<net co2 assimilation rate limited by
real:: ws(ilg)          !<net co2 assimilation rate limited by
                        !!sucrose availability (umol/m2/s)=JE IN PHTSYN
real:: photon(ilg)      !<electron transport rate (umol/m2/s)

!    --------------internal variables above this line------------------

real:: term1(ilg), term2(ilg), term3(ilg)   !temporal terms for
                        !photosynthesis calculations
real:: psna(ilg), psnb(ilg), psne(ilg) !coefficients for quadratric
                        !solution of net photosynthesis
real:: mI(ilg), mII(ilg)!coefficients of the solutions for net psn
!real:: beta1, beta2
real:: temp_b, temp_c, temp_r, temp_q1, temp_q2, temp_jp

!    -----------testing------------YW May 06, 2015 --------------------
real::  dr2, SWin_ex, pdd(ilg), ta(ilg),daylength(ilg)
!     -----------temporal terms above this line------------------------

real     DELT,TFREZ, HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,SPHW,&
        SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,TCGLAC,CLHMLT,CLHVAP

!    -------------common block parameters above this line--------------

COMMON /CLASS1/ DELT,TFREZ
COMMON /CLASS4/ HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY, &
                SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE, &
                TCGLAC,CLHMLT,CLHVAP


! Do the
do   i = 1, ilg
    if (iday == 2)    then
        pdd(i) = 0.
    elseif (tsurfk(i)>tfrez)           then
        pdd(i)=pdd(i)+(tsurfk(i)-tfrez)*DELT/86400.
    endif
end do

!     PHOTOSYNTHESIS COUPLING OR CURVATURE COEFFICIENTS
!real, parameter :: BETA1 = 0.950
!real, parameter :: BETA2 = 0.990

!!    find the light level (parm) at the ground surface for moss photosynthesis,
!!    and a scaling factor degree of saturation of the moss layer phenology
!!    parm is in umol/m2/s and converted from qswnv in W/m2

do 100   i = 1, ilg
    parm(i)= qswnv(i)*4.6
    o2(i)  = 20.9/100.0 * pres(i)
    co2a(i) = co2conc(i)/1000000.0 * pres(i)
    wj(i)     = 0.0
    ws(i)     = 0.0
    wc(i)     = 0.0
    anmoss(i) = 0.0
    rmlmoss(i)= 0.0
    tsurf(i)= tsurfk(i)-tfrez
    tmossk(i) = tsurfk(i)
    tmoss(i)= tsurf(i)
100   continue

!!    phenology  water factor on mosses ,grow when temperature > -4.0
!!    and snowpack < 0.15 m
do 200 i = 1, ilg
    if (zsnow(i) .gt. 0.05 .or. tsurf (i) .lt. 0.5)  then
        pheno(i) = 0.0
    else
        pheno(i) = 1.0
    endif
200   continue

!!    ** water content used for the living moss (depth dmoss)
!!    dmoss is an input and site specific. Preferably make dmoss a function
!!    of Cmoss and ipeatland (different species in fens and bogs)
!!    observed range of wmoss: 5 to 40 in Robrek (2007, 2009), 5 to 25
!!    (Flanagen and Williams 1998)
!!    dmoss is between 2.5 to 5cm based on the species (Lamberty et al. 2006)

do 250   i = 1, ilg
    mmoss(i) = Cmossmas(i)/0.46
    wmoss(i)= thliq(i,1)*rhow/(mmoss(i)/dmoss(i))
    wmosmax(i) = min(45.0, thpms*dmoss(1)*rhow/mmoss(i))
    wmosmin(i) = max(5.0, thmms*dmoss(1)*rhow/mmoss(i))
    wmoss(i)= min(wmosmax(i),max(wmosmin(i),wmoss(i)))
    fwmoss(i)=wmoss(i)+1     !g fresh weight /g dry weight
250   continue

!!    ** moss conductance g_moss umol/m2/s
!!   (Williams and Flanagan, 1998 for Sphagnum). follow MWM, fwmoss is
!!    the mosswat_fd in MWM. Empirical equation is only valid up to
!!   fwmoss=13, above 13 apply a linear extension to the equation.

do 300 i = 1, ilg
    if (fwmoss(i) .le. 13.0)      then
    g_moss(i)=-0.195+0.134*fwmoss(i) - 0.0256*(fwmoss(i)) &
        **2 + 0.00228*(fwmoss(i))**3 - 0.0000984* &
        (fwmoss(i))**4 + 0.00000168*(fwmoss(i))**5
    else
    g_moss (i) = -0.000447 * fwmoss(i) + 0.0489
    endif
    g_moss (i) = g_moss(i) * 1000000.0
    g_moss(i) = max(0.0, g_moss(i))

!!    ** moss surface evaporation coefficient
!!    controled by the degree of saturation in moss, pass to TSOLVC and TSOLVE
!!    apply a similar equation of soil surface cevap in TPREP
!!   CEVAP = 0.25*[1 – cos(THLIQ*pi/THFC)]^2

    dsmoss(i) = (wmoss(i)-wmosmin(i))/(wmosmax(i)-wmosmin(i))
    if (dsmoss(i) .lt.   0.001)       then
        ievapms(i) = 0
        cevapms(i) = 0.
    elseif (dsmoss(i) .ge. 1.0)       then
        ievapms(i) = 1
        cevapms(i) = 1.0
    else
        ievapms(i) = 1
        cevapms(i) = 0.25*(1.0-cos(3.14159*dsmoss(i)))**2
    endif


!!    ** moss water content effect on dark respiration
!!   in MWM and PDM an optimal wmoss is at 5.8 gw/gdw(fig. 2e, Frolking et al.,1996)
!!   Recent studies show weak but significant increases of sphagnum dark respiration
!!   with moss water content above 5.8 gw/gdw (Adkinson and Humphreys, 2011 and ref.)
!!   this change has improved the ER simulation greatly
    if (wmoss(i) .lt. 0.4)   then
        mwce (i) = 0.0
    elseif(wmoss(i) .lt. 5.8 .and. wmoss(i) .gt. 0.4) then
        mwce(i) = 0.35*wmoss(i)**(2.0/3.0)-0.14
    else
!              mwce(i) =-0.04*wmoss(i)+1.232 ! in MWM and PDM
        mwce(i)= 0.01*wmoss(i)+0.942
    endif
300   continue

!!    ** moss dark respiration
!!    observed range of rmlmoss 0.60 to 1.60 umol/m2/s (e.g. Adkinson 2006)
do 350 i = 1, ilg
    q10rmlmos(i)=(3.22-(0.046*tmoss(i)))**((tmoss(i)-25.0)/10.0)
    rmlmoss(i) = rmlmoss25*mwce(i)*q10rmlmos(i)

!!   ** moss photosynthesis
!!   calculate bc (coefficient used for Wc, limited by Rubisco)

    tau(i) = tau25m*exp((tmossk(i)-tref)*ektau/(tref*gasc*tmossk(i)))
    gamma(i) = .5 * o2(i)/ tau(i)
    kc(i) = kc25*exp((tmossk(i)-tref)*ec/(tref*gasc*tmossk(i)))
    ko(i) = ko25*exp((tmossk(i)-tref)*eo/(tref*gasc*tmossk(i)))
    bc(i)  = kc(i) * (1.0 + (o2(i)/ ko(i)))

!    seasonal change of Vcmax, sphagnum (fig. 6, Williams and Flanagan, 1998)
!    from May 1st to september 1st Vmax is maximum
!         if(iday .lt. 121)        then
!              vcmax25(i) = 6.0
!         else if(iday .gt. 245)   then
!              vcmax25(i) = 7.0
!         else
!              vcmax25(i) = 13.5
!         endif

!!     use a function that is also valid for the southern hemisphere
!!      dr2= 1.0+0.033*cos(2.0*pi*day/365.0)    ! dr inverse of solar sun distance
!!      swin_ex=1367.0*dr2*coszs
    if (daylength(i)>14.0 .and. pdd(i)>200. .and. pdd(i)<2000. ) then
        vcmax25(i) = 14.0
    else
        vcmax25(i) = 6.5
    endif
!
!      IF (IYEAR ==2008)          THEN
!      write(90,6991) iday,imin,ihour,pdd,daylength,vcmax25,tmoss
!6991  format(3I5,5f9.3)
!      ENDIF

    vcmax(i) = vcmax25(i)*exp((tmossk(i)-tref)*evc/(tref*gasc*tmossk(i)))

!!     calculate ws (phototysnthesis rate limited by transport capacity)
!!    = js in PHTSYN3

    if (coszs(i).gt.0.0)     THEN
        ws(i) = 0.5*vcmax(i)
    endif


!!    calculate the maximum electorn transport rate Jmax (umol/m2/s)
!!   1.67 = vcmax25m/jmax25m ratio

    jmax25(i) = 1.67 * vcmax25(1)
    term1(i)=exp(((tmossk(i)/tref)-1.)*ej/(gasc * tmossk(i)))
    term2(i)=1+exp(((tref*sj)-hj)/(tref*gasc))
    term3(i)=1+exp(((sj*tmossk(i))-hj)/(gasc*tmossk(i)))
    jmax(i)=jmax25(i)*term1(i)*term2(i)*term3(i)
    if (jmax(i) .gt. 0.0)              then
        photon(i)=alpha_moss*parm(i)/sqrt(1.0+(alpha_moss**2*parm(i)**2/(jmax(i)**2)))
    else
        photon(i) = 0.0     !electron trasport rate in mosses
    end if

!!    calculate Wj, Wc (Farquhar and Caemmerer 1982)
!!   wj = light limited, = je in PHTSYN

    wj(i) = photon(i)*(co2a(i)-gamma(i))/(4.*co2a(i)+(8.*gamma(i)))

!!    carboxylase(rubisco) limitation = jc in PHTSYN
    wc(i) = vcmax(i)*(co2a(i)-gamma(i))/(co2a(i)+bc(i))

!!    Choose the minimum of Wj and Wj both having the form:
!!   W = (a Ci - ad) / (e Ci + b)
!!   Then set a, b, d and e for the quadratic solution for net photosynthesis.

    if(wj(i) < wc(i))then
        psnb(i) = 8. * gamma(i)
        psna(i) = photon(i)
        psne(i) = 4.0
    elseif(wc(i) < wj(i))then
        psnb(i) = bc(i)
        psna(i) = vcmax(i)
        psne(i) = 1.0
    end if
350    continue

!!    Calculate net and gross photosynthesis by solve the quadratic equation
!!   first root of solution is net photosynthesis An= min(Wj,Wc) - Rd
!!   gross photosynthesis GPP = min(Wc,Wj) = An + Rd

do 400   i = 1, ilg
    if (psna(i) .gt.0.0)                              then
        mI(i)=rmlmoss(i)-(psnb(i)*g_moss(i)/pres(i)/psne(i)) &
            - (co2a(i)*g_moss(i)/pres(i))-(psna(i)/psne(i))

        mII(i)=(psna(i)*co2a(i)*g_moss(i)/pres(i)/psne(i))- &
        (rmlmoss(i)*co2a(i)*g_moss(i)/pres(i))-(rmlmoss(i)* &
        psnb(i)*g_moss(i)/pres(i)/psne(i))-(psna(i)*gamma(i)* &
        g_moss(i)/pres(i)/psne(i))
    else
        mI(i)= 0.0
        mII(i)= 0.0
    endif
    anmoss(i) = (-mI(i)-(mI(i)*mI(i)-4*mII(i))**0.5)/2
!         anmoss(i) = min(anmoss(i), ws(i))
    anmoss(i) = pheno(i)*min(anmoss(i), ws(i))   !YW April 22, 2015  add phenology of mosses
400       continue

!    write to output file midday             CT16D_G
if (iyear .eq. 2004 .and. ihour==12)        then
    write(98,6998) iday,tmoss,cevapms,fwmoss,thliq(1,1),&
    dsmoss,g_moss, wmoss,rmlmoss,mwce,q10rmlmos,wmosmax,wmosmin

6998      format(I6,20f12.4)
endif

return
end subroutine mosspht
!>@}

! ---------------------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------------------

!>\ingroup peat_soil_het_resp
!!@{
!! grid average peat soil heterotrophic respiration subroutine

subroutine  decp(il1,il2,iyear,iday,ihour,imin,ipeatland,isand, &
                litrmassms,hpd, wtable,tbar, thliq, thice,thpor,bi,zbotw, &
                delzw,psisat,tfrez, jdsty,jdstd,jdendy,jdendd, &
!    -------------- inputs above this line, outputs below -------------
                litresms, socresp, resoxic, resanoxic)


!   History:

! J. Melton Sep 26 2016
!   - Bring into rest of model and model formatting, convert to doxygen compatible code

!   Created by Yuanqiao Wu, March 20, 2015

!   ----------------------------------

use ctem_params,      only :icc, ilg,ignd,zero,tanhq10,dctmin,dcbaset,bsrateltms

implicit none

!     inputs-----------------------------------------------------------
integer  iyear, iday, ihour, imin ,i,j, il1, il2, isand(ilg,ignd)
integer  jdsty,jdstd,jdendy,jdendd
integer    ipeatland(ilg) !0 = not peatland, 1 = bog, 2 = fen
integer:: lewtable(ilg) !layer index of the water table layer
real     hpd(ilg) , wtable(ilg), tbar(ilg,ignd), &  !(K)
            thliq(ilg,ignd),     thice(ilg,ignd),    thpor(ilg,ignd), &
            bi(ilg,ignd),       zbotw(ilg,ignd),    delzw(ilg,ignd), &
            frac(ilg),          tfrez,              litrmassms(ilg)

!     --------------outputs C fluxes in umolCO2/m2/s--------------------

real:: litresms(ilg)        !< moss litter respiration ($\mu mol CO_2 m^{-2} s^{-1}$)
real:: socresp(ilg)         !< soil C respiration ($\mu mol CO_2 m^{-2} s^{-1}$)
real:: resoxic(ilg)         !< respiration rate of the oxic compartment ($\mu mol CO_2 m^{-2} s^{-1}$)
real:: resanoxic(ilg)       !< respiration rate of the oxic compartment ($\mu mol CO_2 m^{-2} s^{-1}$)

!     internal variables------------------------------------------------

real:: Cso (ilg)            !< carbon mass in the oxic compartment ($kg C m^{-2}$)
real:: Csa (ilg)            !< carbon mass in the anxic compartment ($kg C m^{-2}$)
real:: fto(ilg)             !< temperature factor of the oxic soil respiration
real:: fta(ilg)             !< temperature factor of the anoxic soil respiration
real:: tsoilo(ilg)          !< average temperature in the oxic compartment(C)
real:: tsoila(ilg)          !< average temperature of the anoxic compartment (C)
real:: ewtable(ilg)         !< effective water table depth (m)
real:: ratescpo(ilg)        !< oxic respiration rate constant($\mu mol CO_2 [kg C]^{-2} s^{-1}$)
real:: ratescpa(ilg)        !< anoxic respiration rate constant($\mu mol CO_2 [kg C]^{-2} s^{-1}$)
real:: psi(ilg,ignd)        !< matrix potential of soil layers (Pa)
real:: psisat(ilg,ignd)     !< saturated matrix potential in soil (m)
real:: q10funcms(ilg)       !< q10 fuction for moss litter respiration
real:: ltrmosclms(ilg)      !< moisture scale for litter respiration
real:: litrtempms(ilg)      !< temperature of the moss litter
real:: litpsims(ilg)        !< matrix potential of litter
real:: litrq10ms(ilg)       !< q10 coefficient as a functon of T (as in CTEM)
real:: soilq10o(ilg)        !< q10 coefficient of oxic soil
real:: soilq10a(ilg)        !< q10 coefficient of anoxic soil

!    ------------------------------------------------------------------
!
!    initialization
do 10 i = il1, il2
    litresms(i) = 0.0
    socresp(i) = 0.0
    resoxic(i) = 0.0
    resanoxic(i) = 0.0
10    continue

!    ** calculate soil respiration from peat

!>    find the effective water table depth and the layer index to devide
!!    the peat soil into two compartment

do 20      i= il1, il2
    ewtable(i)= wtable(i)
    if (ewtable(i) .le. 0.0)                then
        lewtable(i) = 0
    elseif (ewtable(i) .le. 0.1)            then
        lewtable(i) = 1
    else
        do j = 1, ignd
            if (ewtable(i) .gt. zbotw(i,j))        then
                lewtable(i)=j+1
            endif
        enddo
    endif

!>    find the temperature in litter, oxic soil and anoxic soil in kalvin
!!    lewtable is the layer index of the water table layer, lewtable = 0
!!    indicates WTD is above the ground surface.
!!    Set the oxic layer temperature to dctmin (minimum soil respiration
!!    temperture) when the entire soil is in the anoxic zone.

    tsoilo(i) = 0.0
    tsoila(i) = 0.0
    if (lewtable(i) .eq. 0)   then  !WT is at or above the surface
        do j = 1, ignd
            tsoila(i) = tsoila(i)+tbar(i,j)*delzw(i,j)
        enddo
        tsoila(i)=tsoila(i)/zbotw(i,ignd)
        tsoilo(i) = dctmin
    elseif (lewtable(i) .eq. 1) then     !WT is at the first layer
        tsoilo(i)=tbar(i,1)
        do j= lewtable(i)+1, ignd
            tsoila(i)=tsoila(i)+ tbar(i,j)*delzw(i,j)
        enddo
        tsoila(i)=(tsoila(i)+tbar(i,1)*(zbotw(i,lewtable(i)) &
                -ewtable(i)))/(zbotw(i,ignd)-ewtable(i))
    else                                !WT is below layer 1
        do j = 1,lewtable(i)-1
            tsoilo(i) = tsoilo(i)+ tbar(i,j)*delzw(i,j)
        enddo
        tsoilo(i)= (tsoilo(i)+ tbar(i,lewtable(i))* &
                (ewtable(i)-zbotw(i,lewtable(i)-1)))/ewtable(i)
        do j = ignd,lewtable(i)+1,-1
            tsoila(i) = tsoila(i)+ tbar(i,j)*delzw(i,j)
        enddo
        tsoila(i)= (tsoila(i)+tbar(i,lewtable(i))* &
                (zbotw(i,lewtable(i))-ewtable(i)) )/ &
                (zbotw(i,ignd)-ewtable(i))
    endif
20    continue

!>    calculate the temperature multiplier (ftsocres) for oxic and anoxic
!!   soil compartments

do 30      i = il1, il2
        soilq10o(i) = tanhq10(1) + tanhq10(2)* &
&            ( tanh( tanhq10(3)*(tanhq10(4)-(tsoilo(i)-tfrez))))
    fto(i)= soilq10o(i)**(0.1*(tsoilo(i)-tfrez-15.0))
        soilq10a(i) = tanhq10(1) + tanhq10(2)* &
&            ( tanh( tanhq10(3)*(tanhq10(4)-(tsoila(i)-tfrez))))
    fta(i)= soilq10a(i)**(0.1*(tsoila(i)-tfrez-15.0))
30    continue

!>    find the heterotrophic respiration rate constant in tje oxic and
!!    anoxic (unit in yr-1), based on Fig.2b in Frolking 2001

do 40 i = il1, il2
    if (ipeatland(i) == 1)               then      !bogs
        if (ewtable(i) .lt. 0.0)                             then
            ratescpo(i)=0.0
            ratescpa(i)=-0.183*exp(-18.0*hpd(i))+0.03*hpd(i)+0.0134
        elseif (ewtable(i).lt.0.30 .and.ewtable(i).ge.0.0) then
            ratescpo(i)=0.009*(1-exp(-20.*ewtable(i)))+0.015*ewtable(i)
            ratescpa(i)=0.009*exp(-20.*ewtable(i))-0.183*exp(-18.*hpd(i))-0.015*ewtable(i)+0.0044
        elseif (ewtable(i) .ge. 0.30)                        then
            ratescpo(i)=0.0134-0.183*exp(-18.*ewtable(i))+0.003*ewtable(i)
            ratescpa(i)=-0.183*exp(-18.*hpd(i))+0.003*(hpd(i)-wtable(i))+0.183*exp(-18.*ewtable(i))
!                 ratescpa(i)=-0.183*exp(-18*hpd(i))+0.003*(hpd(i)
!    1                 -wtable(i))+0.183*exp(-18*ewtable(i))-0.004504   !for continuity
        endif
    elseif (ipeatland(i) == 2)           then      !fens
        if (ewtable(i) .lt. 0.0)                             then
            ratescpo(i) = 0.0
            ratescpa(i) = 0.01512 -1.12*exp(-25.*hpd(i))
        elseif(ewtable(i).lt. 0.30 .and.ewtable(i).ge. 0.0)then
            ratescpo(i) = -0.01*exp(-40.*ewtable(i))+0.015*ewtable(i)+0.01
            ratescpa(i) = abs(-0.01*exp(-40.*ewtable(i))-1.12*exp( &
                            -25.*hpd(i))+0.015*ewtable(i)+0.005119)
        elseif(ewtable(i) .ge. 0.30)                         then
        ratescpo(i) = 0.01512-1.12*exp(-25*ewtable(i))
        ratescpa(i) = -1.12*(exp(-25.*hpd(i))-exp(-25.*ewtable(i)))
        endif
    endif
!!    converts respiration rates from kg c/kg c.year to u-mol co2/kgC/s
    ratescpo(i) = 2.64 * ratescpo(i)
    ratescpa(i) = 2.64 * ratescpa(i)

!>    find the carbon storage in oxic and anoxic compartments (Cso. Csa)
!!    The water table depth delineates the oxic and anoxic compartments.
!!    functions (R**2 = 0.9999) determines the carbon content of each
!!    compartment from a peat bulk density profile based on unpulished
!!    data from P.J.H. Richard (described in fig. 1, Frokling et al.(2001)
!!    conversion of peat into carbon with 48.7% (Mer Bleue unpublished data,
!!    Moore)

    Cso(i) = (4056.6*ewtable(i)**2+72067.0*ewtable(i))*0.487/1000.0
    Csa(i) = ((4056.6*hpd(i)**2+72067.0*hpd(i))*0.487/1000.0)-Cso(i)

40    continue

!>    find the soil respiration rate in Cso and Csa umol/m2/s.
!!    Moisture multiplier (0.025) indicates rate reduction in decomposition due
!!    to anoxia (Frolking et al. 2001), only applied to anoxic layer
do 50 i = il1, il2
    resoxic(i)   = ratescpo(i)*Cso(i)*fto(i)
    resanoxic(i) = ratescpa(i)*Csa(i)*fta(i) *0.025
    socresp(i)   = resoxic(i) + resanoxic(i)
50    continue

!!    **calcualte litter respiration of moss

!!    first find the matrix potential of the soil layers
do 60 j = 1, ignd
    do 60 i = il1, il2
        if(isand(i,j).eq.-3.or.isand(i,j).eq.-4) then  !ice or rock
        psi(i,j) = 10000.0 ! a large number so that ltrmoscl = 0.2
        else
        if (thliq(i,j)+ thice(i,j)+0.01 < thpor(i,j).and.  tbar(i,j) <273.16) then
            psi(i,j) = 0.001
        elseif (thice(i,j) > thpor(i,j))    then
            psi(i,j) = 0.001   !set to saturation
        else
            psi(i,j)=psisat(i,j)*(thliq(i,j)/(thpor(i,j)-thice(i,j)))**(-bi(i,j))
            endif
    endif
60    continue

!!    litter in peatlands can be saturated so we limit the rate by high
!!    moisuture level similar to soil in CTEM, but less effectively (the
!!    min moisture factor is at 0.5 for moss litter but at 0.2 for soil).
do 70 i = il1, il2
    litpsims(i) = psi(i,1)
!    limit of ltrmoscalms at saturation YW April 10, 2015
    if (litpsims(i) .gt. 10000.0) then
            ltrmosclms(i) = 0.2
        elseif (litpsims(i).le. 10000.0 .and.litpsims(i).gt. 6.0) then
        ltrmosclms(i)=1.0-0.8*((log10(litpsims(i))-log10(6.0))/(log10(10000.0)-log10(6.0)) )**1
        elseif (litpsims(i).le. 6.0 .and.litpsims(i).gt. 4.0) then
            ltrmosclms(i)=1.0
        elseif (litpsims(i).le. 4.0 .and. litpsims(i).gt.psisat(i,1))then
            ltrmosclms(i)=1.0-0.99*((log10(4.0)-log10(litpsims(i)))/(log10(4.0)-log10(psisat(i,1))))
        elseif (litpsims(i) .le. psisat(i,1))                     then
        ltrmosclms(i)=0.01
        endif
!    -----------------------------------------------------------------

        ltrmosclms(i)=max(0.0,min(ltrmosclms(i),1.0))

!!    find the temperature factor for moss litter respiration
    litrtempms(i)=tbar(i,1)-tfrez
        litrq10ms(i) = tanhq10(1) + tanhq10(2)* &
&            ( tanh( tanhq10(3)*(tanhq10(4)-litrtempms(i))  ) )
    q10funcms(i)= litrq10ms(i)**(0.1*(litrtempms(i)-15.0))

!!    calculate the litter respiration rate in mosses and converts it
!!     from kg c/kg c.year to u-mol co2/kg c.s using 2.64
    litresms(i)=ltrmosclms(i)*litrmassms(i)*bsrateltms*2.64*q10funcms(i)
70        continue

!     write peat soil respiration details .CT15D_G

if ((iyear .ge. jdsty) .and. (iyear .le. jdendy)) then
    if ((iday .ge. jdstd) .and. (iday .le. jdendd)) then

write(97, 6997) litresms, litpsims(1), psisat(1,1),ltrmosclms, &
        litrmassms, tbar(1,1), q10funcms ,litrtempms, &
        ratescpo, ratescpa, Cso,Csa, fto,fta, &
        resoxic, resanoxic, frac, &
        tsoila-tfrez, tsoilo-tfrez,ewtable,real(lewtable), &
        tbar(1,1)-tfrez, tbar(1,2)-tfrez, tbar(1,3)-tfrez, &
        thliq(1,1)
6997  format(50f10.3)
    endif
endif

return

end subroutine decp
!>@}
end module