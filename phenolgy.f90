!>\file
!!Canadian Terrestrial Ecosystem Model (CTEM)
!!Phenology, Leaf Turnover & Mortality Subroutine
!!
!!
!!
!!The leaf phenology parametrization used in CTEM v. 1.0 is described in detail by \cite Arora2005-6b1. Changes between
!!version 1.0 and 2.0 are limited to parameter values and the parametrization is briefly described here. There are four
!!different leaf phenological states in which vegetation can be at a given instant: (i) no leaves or dormant,
!!(ii) maximum growth, (iii) normal growth and (iv) leaf fall or harvest. PFTs may go through only some, or all,
!!of these phenological states depending on their deciduousness. A broadleaf cold deciduous tree, for example,
!!transitions through all these four states in a year. In winter, the broadleaf cold deciduous trees are in the
!!no leaves/dormant state; favourable climatic conditions in spring trigger leaf growth and the tree enters the
!!maximum leaf growth state when all the NPP is allocated to leaves to accelerate leaf out; when the LAI reaches a
!!threshold (described below) the tree enters the normal leaf growth state and NPP is also allocated to stem and root
!!components; finally the arrival of autumn triggers leaf fall and the trees go into the leaf fall mode where no carbon
!!is allocated to leaves (but it continues for roots and stems). When all the leaves have been shed, the trees go into
!!the no leaves or dormant state again and the cycle is repeated the next year. The evergreen tree PFTs and the grass
!!PFTs do not enter the leaf fall state and maintain a leaf canopy as long as environmental conditions are favourable.
!!Although drought and cold stress cause accelerated leaf loss compared to the normal leaf turnover from these PFTs,
!!they do not explicitly go into the leaf fall mode where the intent is to lose all leaves in a specified amount of time.
!!
!!The leaf phenological state transitions are dependent upon environmental conditions. In particular, the transition
!!from no leaves/dormant state to the maximum growth state is based on the carbon-gain approach. CTEM uses
!!\f$\textit{virtual}\f$ leaves to assess favourable meteorological conditions for leaf out. The virtual leaves
!!photosynthesize and respire in a manner similar to normal leaves except the carbon gain or loss is not taken into
!!account in vegetation's carbon balance. A positive net leaf photosynthesis rate (\f$G_{canopy,net}\f$, Eq. \ref{Gnet})
!!for the virtual leaves over seven consecutive days indicates the arrival of favourable growth conditions and triggers
!!leaf onset and the associated transition from the no leaves/dormant state to the maximum leaf growth state, when the
!!entire positive NPP is allocated to leaves (\f$a_{fL} = 1\f$, \f$a_{fS} = a_{fR} = 0\f$). When LAI reaches \f${LAI}_{thrs}\f$
!!then the vegetation switches to the normal growth mode and positive NPP is allocated to all three vegetation
!!components -- leaves, stem and roots (\f$a_{fL}, a_{fS}, a_{fR} > 0\f$). \f${LAI}_{thrs}\f$ is calculated as
!!\f[{LAI}_{thrs} = L_f \left[ {SLA}\left(\frac{C_S + C_R}{\eta}\right)^{1/\kappa} \right].\f]
!!
!!The PFT-specific \f$L_f\f$ term (see also ctem_params.f90) calculates \f${LAI}_{thrs}\f$ to be typically between
!!40 and 50\,{\%} of the maximum LAI that a given amount of stem and root biomass can support (based on the terms in the square
!!brackets and Eq. (\ref{propwoody}). \f$SLA\f$ is the specific leaf area (Eq. \ref{sla}).  This rule for transition
!!from a maximum to a normal growth state is also used for evergreen tree PFTs and grass PFTs. Similar to
!!\f${LAI}_{thrs}\f$, the LAI of virtual leaves is \f$7.5\,{\%}\f$ of the maximum LAI a given amount of root
!!and stem biomass can support for tree and crop PFTsand \f$2.5\,{\%}\f$ for grass PFTs. In
!!addition, the LAI of virtual leaves is constrained to be, at least, \f$0.3\,m^2\,m^{-2}\f$
!!for tree PFTs and \f$0.2\,m^2\,m^{-2}\f$ for crop and grass PFTs.
!!
!!The transition from the normal growth state to the leaf fall state is triggered by unfavourable environmental conditions
!!and shorter day length. Broadleaf deciduous trees transition to the leaf fall state when either: (i) day length is less
!!than \f$11\,h\f$ and the rooting zone temperature drops below \f$11.15\,C\f$ or (ii) when the rooting zone temperature
!!drops below \f$8\,C\f$ regardless of the day length. Needleleaf deciduous tress begin leaf fall after seven consecutive
!!days with daily mean air temperature below \f$-5\,C\f$. Leaf fall occurs over a period of 15 days. In the leaf fall
!!state, the vegetation continues carbon allocation to its root and stem components, but not to leaves (\f$a_{fL} = 0\f$,
!!\f$a_{fS} + a_{fR} = 1\f$). Evergreen trees and grasses do not enter the leaf fall state and neither do the
!!broadleaf drought deciduous trees. The implication for the latter PFT is that if the climate changes and the
!!dry season becomes shorter, then the trees will keep their leaves on for a longer period of time since broadleaf
!!drought deciduous trees lose leaves due to soil moisture stress (described below).
!!
!!The model vegetation is able to transition between the different leaf phenological states in response to changing
!!conditions. For example, a leaf out in spring for broadleaf cold deciduous trees can be interrupted by a cold event
!!when the vegetation goes into a leaf fall state until the return of more favourable conditions.
!!
!!Leaf litter generation is caused by normal turnover of leaves (\f$\Omega_N\f$, \f$day^{-1}\f$) and also
!!due to cold (\f$\Omega_C\f$, \f$day^{-1}\f$) and drought (\f$\Omega_D\f$, \f$day^{-1}\f$) stress,
!!both of which contribute to seasonality of LAI. For example, the leaf loss associated with drought and reduced
!!photosynthesis during the dry season are the principal causes of the seasonality of LAI for the broadleaf drought deciduous tree PFT.
!!
!!The conversion of leaf carbon to leaf litter (\f$D_L\f$, \f$kg\,C\,m^{-2}\,day^{-1}\f$) is expressed as
!!
!!\f[ \label{cltod} D_L = C_L[1 - \exp(-\Omega_N - \Omega_C - \Omega_{D})],\f]
!!
!!where (\f$\Omega_{N,C,D}\f$, \f$day^{-1}\f$) are the leaf loss rates associated with normal turnover of leaves and
!!the cold and drought stress. The rate of normal turnover of leaves is governed by PFT-specific leaf lifespan
!!(\f$\tau_L\f$, \f$yr\f$) as \f$\Omega_N= 1/365 \tau_L\f$ (see also ctem_params.f90}
!!for PFT specific values of \f$\tau_L\f$).  The leaf loss rate associated with cold stress
!!(\f$\Omega_C\f$) is calculated as
!!
!!\f[ \label{gamma_cold} \Omega_C = \Omega_{C,max}L_{cold}^3, \f]
!!
!!where \f$\Omega_{C,max}\f$ (\f$day^{-1}\f$, see also ctem_params.f90) is the maximum cold stress loss rate.
!!\f$L_{cold}\f$ is a scalar that varies between 0 and 1 as
!!
!!\f[ \label{cldls} L_{cold} = \begin{cases} 1, \quad T_a < \left(T_{cold}^{leaf} -
!!5\right) \\ 1 - \frac{T_a - \left(T_{cold}^{leaf} - 5\right)}{5}, \\ \quad T_{cold}^{leaf}
!!> T_a > (T_{cold}^{leaf} - 5) \\ 0, \quad T_a > T_{cold}^{leaf} ,\\ \end{cases} \f]
!!
!!where \f$T_{cold}^{leaf}\f$ is a PFT-specific temperature threshold below which a PFT experiences damage
!!to its leaves promoting leaf loss (see also ctem_params.f90) and \f$T_a\f$ is the daily mean air
!!temperature (\f$C\f$).  The leaf loss rate due to drought stress is calculated in a similar manner
!!
!!\f[ \label{gamma_dry} \Omega_{D} = \Omega_{D,max}\,(1-\phi_{root})^3, \f]
!!
!!where \f$\Omega_{D,max}\f$ (\f$day^{-1}\f$, see also ctem_params.f90) is the maximum drought
!!stress loss rate and \f$\phi_{root}\f$ (Eq. \ref{degsoilsat}) is the degree of soil saturation in the rooting zone.
!!
!!
!!
subroutine phenolgy(gleafmas, bleafmas,  &
             &           il1,      il2,     tbar, &
             &          thliq,   THLW,  THFC,       ta,&
             &           anveg,     iday,     radl, roottemp,&
             &       rmatctem, stemmass, rootmass,     sort,&
             &       nol2pfts,  fcancmx, isand,&
!     ------------------ inputs above this line ----------------------
             &       flhrloss, leaflitr, lfstatus,  pandays,&
             &       colddays,thice)
!     --- variables which are updated and outputs above this line ----
!
!               Canadian Terrestrial Ecosystem Model (CTEM)
!               Phenology, Leaf Turnover & Mortality Subroutine
!
!     7   Nov 2016  - Allow adaptions due to peatland scheme of Yuanqiao Wu.

!     14  Jan 2016  - There was a bit of hardwired code for 3 soil layers, that has been
!                     fixed to allow >3.

!     17  Jan 2014  - Moved parameters to global file (ctem_params.f90)
!     J. Melton
!
!     22  Jul 2013  - Add in module for parameters
!     J. Melton
!
!     24  Sep 2012  - add in checks to prevent calculation of non-present
!     J. Melton       pfts
!
!     15  Apr. 2003 - this subroutine calculates the leaf status 
!     V. Arora        for ctem's pfts and leaf litter generated by
!                     normal turnover of leaves and cold and drought 
!                     stress. crop harvest is also modelled in this
!                     subroutine, and for grasses green leaves are
!                     converted into brown.
!
!     inputs 
!
!     icc       - no. of ctem plant function types, currently 9
!     ignd      - no. of soil layers (currently 3)
!     ilg       - no. of grid cells in latitude circle
!     l2max     - maximum number of level 2 ctem pfts
!     ican      - number of class pfts
!
      use ctem_params,        only : kn, pi, zero, kappa, eta, lfespany,&
     &                               fracbofg, specsla,ilg,ignd,icc,kk,&
     &                               ican, cdlsrtmx, drlsrtmx, drgta,&
     &                               colda, lwrthrsh, dayschk, coldlmt,&
     &                               coldthrs, harvthrs, flhrspan,&
     &                               thrprcnt, roothrsh, pftlist

      implicit none
!
      integer il1                 !<il1=1
      integer il2                 !<il2=ilg
      integer i, j, k
      integer iday                !<day of year
      integer isand(ilg,ignd) !<
      integer n, m, k1, k2
!
      integer sort(icc)           !<index for correspondence between 9 pfts and the 12 values in parameters vectors
      integer nol2pfts(ican)      !<number of level 2 ctem pfts
!
      real gleafmas(ilg,icc)      !<green or live leaf mass in \f$kg c/m
      real bleafmas(ilg,icc)      !<brown or dead leaf mass in \f$kg c/m
      real ta(ilg)                !<air temperature, k
      real tbar(ilg,ignd)         !<soil temperature, k
      real thliq(ilg,ignd)        !<liquid soil moisture content in 3 so
      real thice(ilg,ignd)    !YW May 04, 2015 
      real sand(ilg,ignd)         !<
      real clay(ilg,ignd)         !<
      real anveg(ilg,icc)         !<net photosynthesis rate of ctem's pfts, umol co2/m2.s
      real leaflitr(ilg,icc)      !<leaf litter generated by normal turnover, cold and drought stress, 
                                  !<and leaf fall/harvest, \f$kg c/m^2\f$
      real roottemp(ilg,icc)      !<root temperature, which is a function of soil temperature of course, k.
      real rmatctem(ilg,icc,ignd) !<fraction of roots in each soil layer for each pft
      real stemmass(ilg,icc)      !<stem mass for each of the 9 ctem pfts, \f$kg c/m^2\f$
      real rootmass(ilg,icc)      !<root mass for each of the 9 ctem pfts, \f$kg c/m^2\f$
      real fcancmx(ilg,icc)       !<max. fractional coverage of ctem's 9 pfts, but this can be modified 
                                  !<by land-use change, and competition between pfts
!
      integer pandays(ilg,icc)    !<counter for positive net photosynthesis (an) days for initiating leaf onset
      integer lfstatus(ilg,icc)   !<integer indicating leaf status or mode
                                  !<1 - max. growth or onset, when all npp is allocated to leaves
                                  !<2 - normal growth, when npp is allocated to leaves, stem, and root
                                  !<3 - fall for dcd trees/harvest for crops, when allocation to leaves is zero.
                                  !<4 - no leaves
      integer chkmode(ilg,icc)    !<indicator for making sure that leaf status is updated
      integer colddays(ilg,2)     !<cold days counter for tracking days below a certain temperature threshold 
                                  !<for ndl dcd and crop pfts.
!
      real sla(icc)               !<specific leaf area
      real ailcg(ilg,icc)         !<green lai
      real ailcb(ilg,icc)         !<brown lai

      real THFC(ilg,ignd)      !<field capacity soil moisture content both calculated in allocate subroutine
      real THLW(ilg,ignd)       !<wilting point soil moisture content
      real day                    !<
      real radl(ilg)              !<latitude in radians
      real theta                  !<
      real decli                  !<
      real term                   !<
      real daylngth(ilg)          !<
      real nrmlloss(ilg,icc)      !<leaf loss due to normal turnover
      real betadrgt(ilg,ignd)     !<(1 - drought stress)
      real drgtstrs(ilg,icc)      !<drought stress term
      real drgtlsrt(ilg,icc)      !<drought loss rate
      real drgtloss(ilg,icc)      !<leaf loss due to drought stress
      real coldloss(ilg,icc)      !<leaf loss due to cold stress
      real coldstrs(ilg,icc)      !<cold stress term
      real coldlsrt(ilg,icc)      !<cold loss rate
      real flhrloss(ilg,icc)      !<fall & harvest loss for bdl dcd plants and crops, respectively, \f$kg c/m^2\f$.
      real lfthrs(ilg,icc)        !<threshold lai for finding leaf status

      character(8) :: pftkind
!>
!!------------------------------------------------------------------
!!Constants and parameters are located in ctem_params.f90
!!
!!---------------------------------------------------------------
!!

      if(icc.ne.12)                            call xit('phenolgy',-1)
!>
!!initialize required arrays to zero
!!
      do 120 j = 1, ignd
        do 130 i = il1, il2
          betadrgt(i,j)=0.0       
130     continue
120   continue
!
      do 140 j = 1,icc
        sla(j)=0.0                
        do 150 i = il1, il2
          ailcg(i,j)=0.0           
          ailcb(i,j)=0.0           
          chkmode(i,j)=0           
          leaflitr(i,j)=0.0       
          nrmlloss(i,j)=0.0        
          drgtstrs(i,j)=0.0        
          drgtlsrt(i,j)=0.0        
          drgtloss(i,j)=0.0        
          coldstrs(i,j)=0.0        
          coldlsrt(i,j)=0.0        
          coldloss(i,j)=0.0        
          lfthrs(i,j)=0.0          
150     continue                  
140   continue
!>
!!initialization ends    
!!
!!------------------------------------------------------------------
!!
!!convert green leaf mass into leaf area index using specific leaf
!!area \f$(sla, m^2 /kg c)\f$ estimated using leaf life span. see bio2str
!!subroutine for more details. 
!!
      do 170 j = 1,icc

        sla(j) = 25.0*(lfespany(sort(j))**(-0.50))
        if(specsla(sort(j)).gt.zero) sla(j)=specsla(sort(j))

        n = sort(j)

        do 180 i = il1,il2
         if (fcancmx(i,j).gt.0.0) then 

          ailcg(i,j)=sla(j)*gleafmas(i,j)
          ailcb(i,j)=sla(j)*bleafmas(i,j)*fracbofg
!>also find threshold lai as a function of stem+root biomass
!!which is used to determine leaf status
          lfthrs(i,j)=((stemmass(i,j)+rootmass(i,j))/eta(n))&
     &     **(1.0/kappa(n))   
          lfthrs(i,j)=(thrprcnt(n)/100.0)*sla(j)*lfthrs(i,j)
!>
!!using green leaf area index (ailcg) determine the leaf status for
!!each pft. loops 190 and 200 thus initialize lfstatus, if this
!!this information is not passed specifically as an initialization quantity. 
!!
          if(lfstatus(i,j).eq.0)then
            if(ailcg(i,j).le.zero)then            
              lfstatus(i,j)=4                      !no leaves
            else if (ailcg(i,j).gt.lfthrs(i,j))then 
              lfstatus(i,j)=2                      !normal growth
            else                                  
              lfstatus(i,j)=4                      !treat this as no leaves
            endif                                  !so that we start growing 
           endif                                   !if possible
!
         endif !fcancmx
180     continue
170   continue
!>
!!knowing lfstatus (after initialization above or using value from
!!from previous time step) we decide if we stay in a given leaf
!!mode or we move to some other mode.
!!
!!we start with the "no leaves" mode
!!----------------------------------
!!
!!add one to pandays(i,j) if daily an is positive, otherwise set it to zero.
!!
      do 220 j = 1, icc
        n = sort(j)
        do 230 i = il1, il2
         if (fcancmx(i,j).gt.0.0) then 
          if(anveg(i,j).gt.zero) then
            pandays(i,j)=pandays(i,j)+1
            if(pandays(i,j).gt.dayschk(n))then
              pandays(i,j)=dayschk(n)
            endif
          else
            pandays(i,j)=0
          endif
         endif
230     continue
220   continue
!>
!!if in "no leaves" mode check if an has been positive over last
!!dayschk(j) days to move into "max. growth" mode. if not we stay
!!in "no leaves" mode. also set the chkmode(i,j) switch to 1.
!!
      do 240 j = 1, icc
        n = sort(j)
        do 250 i = il1, il2
         if (fcancmx(i,j).gt.0.0) then 
          if(chkmode(i,j).eq.0.and.lfstatus(i,j).eq.4)then
            if(pandays(i,j).ge.dayschk(n))then
              lfstatus(i,j)=1        ! switch to "max. growth" mode
              chkmode(i,j)=1         ! mode checked, no more checks further down
            else
              lfstatus(i,j)=4        ! stay in "no leaves" mode
              chkmode(i,j)=1         ! mode checked, no more checks further down
            endif
          endif
         endif
250     continue
240   continue
!>
!!find day length using day of year and latitude. this is to be used for
!!initiating leaf offset for broad leaf dcd trees.
!!
      day=real(iday)
      do 260 i = il1, il2
       theta=0.2163108 + 2.0*atan(0.9671396*tan(0.0086*(day-186.0)))
       decli=asin(0.39795*cos(theta))    !declination
       term=(sin(radl(i))*sin(decli))/(cos(radl(i))*cos(decli))
       term=max(-1.0,min(term,1.0))
       daylngth(i)=24.0-(24.0/pi)*acos(term)
260   continue

!!similar to the way we count no. of days when an is positive, we
!!find no. of days when temperature is below -5 c. we need this to
!!determine if we go into "leaf fall" mode for needle leaf dcd trees.
!!
!!also estimate no. of days below 8 c. we use these days to decide
!!if its cold enough to harvest crops.
!!
      do 280 k = 1, 2
        do 290 i = il1, il2
          if(ta(i).lt.(coldthrs(k)+273.16)) then
            colddays(i,k)=colddays(i,k)+1
            if(colddays(i,k).gt.coldlmt(k))then
              colddays(i,k)=coldlmt(k)
            endif
          else
            colddays(i,k)=0
          endif
290     continue
280   continue

!>
!!even if pandays criteria has been satisfied do not go into max. growth
!!mode if environmental conditions are such that they will force leaf fall or harvest mode.
!!
do 270 i = il1, il2
    do j = 1, icc
        if (fcancmx(i,j).gt.0.0) then
            if(lfstatus(i,j).eq.1.and.chkmode(i,j).eq.1)then
                pftkind = pftlist(j)

                select case (pftkind)
                case ('NdlDcdTr') ! Needle leaf deciduous tree
                    if(ta(i).lt.(coldthrs(1)+273.16))then
                        lfstatus(i,j)=4
                    endif
                case ('BdlDCoTr') ! Broadleaf deciduous cold tree
                    if(roottemp(i,j).lt.(roothrsh+273.16).or. &
                &    (daylngth(i).lt.11.0.and. &
                &     roottemp(i,j).lt.(11.15+273.16)))then
                        lfstatus(i,j)=4
                    endif
                case ('BdlDDrTr') ! Broadleaf deciduous drought dry tree
                    if(roottemp(i,j).lt.(roothrsh+273.16).or. &
                    &    (daylngth(i).lt.11.0.and. &
                    &     roottemp(i,j).lt.(11.15+273.16)))then
                        lfstatus(i,j)=4
                    endif
                case ('CropC3  ')
                    if(ta(i).lt.(coldthrs(2)+273.16))then
                        lfstatus(i,j)=4
                    endif
                case ('CropC4  ')
                    if(ta(i).lt.(coldthrs(2)+273.16))then
                        lfstatus(i,j)=4
                    endif
                case ('BdlDcdSh') ! Broadleaf deciduous shrub
                    if(ta(i).lt.(coldthrs(1)+273.16))then
                        lfstatus(i,j)=4
                    endif
                case ('GrassC3 ')
                    if(ta(i).lt.(coldthrs(2)+273.16))then
                        lfstatus(i,j)=4
                    endif
                case ('GrassC4 ')
                    if(ta(i).lt.(coldthrs(2)+273.16))then
                        lfstatus(i,j)=4
                    endif
                end select
            endif
        endif
    end do
270   continue
!       do 270 i = il1, il2
! !>needle leaf dcd
!        if (fcancmx(i,2).gt.0.0) then
!         if(lfstatus(i,2).eq.1.and.chkmode(i,2).eq.1)then
!           if(ta(i).lt.(coldthrs(1)+273.16))then
!             lfstatus(i,2)=4
!           endif
!         endif
!        endif
! !>
! !!broad leaf dcd cld & dry
!        if (fcancmx(i,4).gt.0.0) then
!         if(lfstatus(i,4).eq.1.and.chkmode(i,4).eq.1)then
!           if(roottemp(i,4).lt.(roothrsh+273.16).or.
!      &    (daylngth(i).lt.11.0.and.roottemp(i,4).lt.(11.15+273.16)))th
!             lfstatus(i,4)=4
! !       write(*,'(a5,2i4,3f10.3)')'lf1=',lfstatus(i,4),iday,roottemp(i
! !     & ,daylngth(i),ailcg(i,4)
!           endif
!         endif
!        endif
!
!        if (fcancmx(i,5).gt.0.0) then
!         if(lfstatus(i,5).eq.1.and.chkmode(i,5).eq.1)then
!           if(roottemp(i,5).lt.(roothrsh+273.16).or.
!      &    (daylngth(i).lt.11.0.and.roottemp(i,5).lt.(11.15+273.16)))th
!             lfstatus(i,5)=4
!           endif
!         endif
!        endif
! !>
!>crops
!       if (icc==9)        then           !YW April 15, 2015
!        if (fcancmx(i,6).gt.0.0) then
!         if(lfstatus(i,6).eq.1.and.chkmode(i,6).eq.1)then
!           if(ta(i).lt.(coldthrs(2)+273.16))then
!             lfstatus(i,6)=4
!           endif
!         endif
!        endif
!        if (fcancmx(i,7).gt.0.0) then
!         if(lfstatus(i,7).eq.1.and.chkmode(i,7).eq.1)then
!           if(ta(i).lt.(coldthrs(2)+273.16))then
!             lfstatus(i,7)=4
!           endif
!         endif
!        endif
!      
!       else               ! when icc !=9, here icc=11 for peatland
! c       dcd shrubs
!        if (fcancmx(i,7).gt.0.0) then
!         if(lfstatus(i,7).eq.1.and.chkmode(i,7).eq.1)then
!           if(ta(i).lt.(coldthrs(1)+273.16))then
!             lfstatus(i,7)=4
!           endif
!         endif
!        endif
!
! c    crops j = 8, 9
!       do j = 8, 9
!        if (fcancmx(i,j).gt.0.0) then
!         if(lfstatus(i,j).eq.1.and.chkmode(i,j).eq.1)then
!           if(ta(i).lt.(coldthrs(2)+273.16))then
!             lfstatus(i,j)=4
!           endif
!         endif
!        endif
!       enddo
!       endif         !if icc==9
! 270   continue
!>
! !!similar to the way we count no. of days when an is positive, we
! !!find no. of days when temperature is below -5 c. we need this to
! !!determine if we go into "leaf fall" mode for needle leaf dcd trees.
! !!
! !!also estimate no. of days below 8 c. we use these days to decide
! !!if its cold enough to harvest crops.
! !!
!       do 280 k = 1, 2
!         do 290 i = il1, il2
!           if(ta(i).lt.(coldthrs(k)+273.16)) then
!             colddays(i,k)=colddays(i,k)+1
!             if(colddays(i,k).gt.coldlmt(k))then
!               colddays(i,k)=coldlmt(k)
!             endif
!           else
!             colddays(i,k)=0
!           endif
! 290     continue
! 280   continue
!>
!!if in "max growth" mode
!!------------------------
!!
!!if mode hasn't been checked and we are in "max. growth" mode, then 
!!check if we are above pft-dependent lai threshold. if lai is more
!!then this threshold we move into "normal growth" mode, otherwise
!!we stay in "max growth" mode so that leaves can grow at their
!!max. climate-dependent rate
!!
      do 300 j = 1, icc
        do 310 i = il1, il2
         if (fcancmx(i,j).gt.0.0) then 
          if(chkmode(i,j).eq.0.and.lfstatus(i,j).eq.1)then
            if(ailcg(i,j).ge.lfthrs(i,j))then
              lfstatus(i,j)=2        ! switch to "normal growth" mode
              chkmode(i,j)=1         
            else if(ailcg(i,j).le.zero) then
              lfstatus(i,j)=4        ! switch to "no leaves" mode
              chkmode(i,j)=1         
              pandays(i,j)=0
            else 
              lfstatus(i,j)=1        ! stay in "max. growth" mode
              chkmode(i,j)=1         
            endif
!>
!!for dcd trees we also need to go into "leaf fall" mode directly from "max. growth" mode.
!!
!!ndl dcd
            pftkind = pftlist(j)

            select case (pftkind)
            case ('NdlDcdTr') ! Needle leaf deciduous tree
    !            if(j.eq.2)then
                if(ailcg(i,j).lt.lfthrs(i,j).and. &
        &        colddays(i,1).ge.coldlmt(1).and. &
        &        ailcg(i,j).gt.zero)then
                    lfstatus(i,j)=3        ! go into "leaf fall" mode
                    chkmode(i,j)=1
                endif
!            endif
!>
!>bdl dcd cold 
            case ('BdlDCoTr') ! Broadleaf deciduous cold tree
                !if(j.eq.4)then
                if( ailcg(i,j).gt.zero.and. &
        &        ((daylngth(i).lt.11.0.and.roottemp(i,j).lt.(11.15+273.16)) &
        &        .or. roottemp(i,j).lt.(roothrsh+273.16)) )then
                    lfstatus(i,j)=3        ! go into "leaf fall" mode
                    chkmode(i,j)=1
                    flhrloss(i,j)=gleafmas(i,j)*(1.0/flhrspan(2))
                endif
           ! endif

!>bdl dcd dry
          !  if(j.eq.5)then
            case ('BdlDDrTr') ! Broadleaf deciduous drought dry tree
                if( ailcg(i,j).gt.zero.and. &
        &        ((daylngth(i).lt.11.0.and.roottemp(i,j).lt.(11.15+273.16)) &
        &        .or. roottemp(i,j).lt.(roothrsh+273.16)) )then
                    lfstatus(i,j)=3        ! go into "leaf fall" mode
                    chkmode(i,j)=1
                endif
            !  endif

            end select

          endif
         endif
310     continue
300   continue
!>
!!if in "normal growth" mode
!!--------------------------
!!
!!if in "normal growth" mode then go through every pft individually
!!and follow set of rules to determine if we go into "fall/harvest" mode
!!
      do 320 i =  il1, il2
        do j = 1, icc
         if (fcancmx(i,j).gt.0.0) then
          if(chkmode(i,j).eq.0.and.lfstatus(i,j).eq.2)then

            pftkind = pftlist(j)

            select case (pftkind)
            case ('NdlEvgTr') ! needle leaf evg tree
                if(ailcg(i,j).lt.lfthrs(i,j).and.ailcg(i,j).gt.zero)then
                    lfstatus(i,j)=1         ! go back to "max. growth" mode
                    chkmode(i,j)=1
                else if(ailcg(i,j).le.zero) then
                    lfstatus(i,j)=4         ! switch to "no leaves" mode
                    chkmode(i,j)=1
                    pandays(i,j)=0
                else
                    lfstatus(i,j)=2         ! stay in "normal growth" mode
                    chkmode(i,j)=1
                endif

!       if (fcancmx(i,1).gt.0.0) then
!        if(chkmode(i,1).eq.0.and.lfstatus(i,1).eq.2)then
!           if(ailcg(i,1).lt.lfthrs(i,1).and.ailcg(i,1).gt.zero)then
!             lfstatus(i,1)=1         ! go back to "max. growth" mode
!             chkmode(i,1)=1
!           else if(ailcg(i,1).le.zero) then
!             lfstatus(i,1)=4         ! switch to "no leaves" mode
!             chkmode(i,1)=1
!             pandays(i,1)=0
!           else
!             lfstatus(i,1)=2         ! stay in "normal growth" mode
!             chkmode(i,1)=1
!           endif
!        endif
!       endif
!>
!!needle leaf dcd
            case ('NdlDcdTr') ! Needle leaf deciduous tree
                if(ailcg(i,j).lt.lfthrs(i,j).and. &
            &      colddays(i,1).ge.coldlmt(1).and. &
            &      ailcg(i,j).gt.zero)then
                    lfstatus(i,j)=3         ! go into "leaf fall" mode
                    chkmode(i,j)=1
                else if(ailcg(i,j).le.zero) then
                    lfstatus(i,j)=4         ! switch to "no leaves" mode
                    chkmode(i,j)=1
                    pandays(i,j)=0
                else
                    lfstatus(i,j)=2         ! stay in "normal growth" mode
                    chkmode(i,j)=1
                endif

!        if (fcancmx(i,2).gt.0.0) then
!         if(chkmode(i,2).eq.0.and.lfstatus(i,2).eq.2)then
!           if(ailcg(i,2).lt.lfthrs(i,2).and.
!      &      colddays(i,1).ge.coldlmt(1).and.
!      &      ailcg(i,2).gt.zero)then
!             lfstatus(i,2)=3         ! go into "leaf fall" mode
!             chkmode(i,2)=1
!           else if(ailcg(i,2).le.zero) then
!             lfstatus(i,2)=4         ! switch to "no leaves" mode
!             chkmode(i,2)=1
!             pandays(i,2)=0
!           else
!             lfstatus(i,2)=2         ! stay in "normal growth" mode
!             chkmode(i,2)=1
!           endif
!         endif
!       endif
!>
!!broad leaf evg

            case ('BdlEvgTr') ! broad leaf evg tree
                if(ailcg(i,j).lt.lfthrs(i,j).and.ailcg(i,j).gt.zero)then
                    lfstatus(i,j)=1         ! go back to "max. growth" mode
                    chkmode(i,j)=1
                else if(ailcg(i,j).le.zero) then
                    lfstatus(i,j)=4         ! switch to "no leaves" mode
                    chkmode(i,j)=1
                    pandays(i,j)=0
                else
                    lfstatus(i,j)=2         ! stay in "normal growth" mode
                    chkmode(i,j)=1
                endif

!        if (fcancmx(i,3).gt.0.0) then
!         if(chkmode(i,3).eq.0.and.lfstatus(i,3).eq.2)then
!           if(ailcg(i,3).lt.lfthrs(i,3).and.ailcg(i,3).gt.zero)then
!             lfstatus(i,3)=1         ! go back to "max. growth" mode
!             chkmode(i,3)=1
!           else if(ailcg(i,3).le.zero) then
!             lfstatus(i,3)=4         ! switch to "no leaves" mode
!             chkmode(i,3)=1
!             pandays(i,3)=0
!           else
!             lfstatus(i,3)=2         ! stay in "normal growth" mode
!             chkmode(i,3)=1
!           endif
!         endif
!        endif
!>
!!broad leaf dcd cold
!!we use daylength and roottemp to initiate leaf offset
        case ('BdlDCoTr') ! Broadleaf deciduous cold tree
          if( ailcg(i,j).gt.zero.and.&
     &    ((daylngth(i).lt.11.0.and.roottemp(i,j).lt.(11.15+273.16))&
     &    .or. roottemp(i,4).lt.(roothrsh+273.16)) )then
            lfstatus(i,j)=3         ! go into "leaf fall" mode
            chkmode(i,j)=1
            flhrloss(i,j)=gleafmas(i,j)*(1.0/flhrspan(2))
          else if(ailcg(i,j).gt.zero.and.ailcg(i,j).lt.lfthrs(i,j))then
            lfstatus(i,j)=1         ! switch to "max. growth" mode
            chkmode(i,j)=1
          else if(ailcg(i,j).le.zero) then
            lfstatus(i,j)=4         ! switch to "no leaves" mode
            chkmode(i,j)=1
            pandays(i,j)=0
            flhrloss(i,j)=0.0
          else
            lfstatus(i,j)=2         ! stay in "normal growth" mode
            chkmode(i,j)=1
          endif

!        if (fcancmx(i,4).gt.0.0) then
!         if(chkmode(i,4).eq.0.and.lfstatus(i,4).eq.2)then
!           if( ailcg(i,4).gt.zero.and.
!      &    ((daylngth(i).lt.11.0.and.roottemp(i,4).lt.(11.15+273.16))
!      &    .or. roottemp(i,4).lt.(roothrsh+273.16)) )then
!             lfstatus(i,4)=3         ! go into "leaf fall" mode
!             chkmode(i,4)=1
!             flhrloss(i,4)=gleafmas(i,4)*(1.0/flhrspan(2))
!           else if(ailcg(i,4).gt.zero.and.ailcg(i,4).lt.lfthrs(i,4))then
!             lfstatus(i,4)=1         ! switch to "max. growth" mode
!             chkmode(i,4)=1
!           else if(ailcg(i,4).le.zero) then
!             lfstatus(i,4)=4         ! switch to "no leaves" mode
!             chkmode(i,4)=1
!             pandays(i,4)=0
!             flhrloss(i,4)=0.0
!           else
!             lfstatus(i,4)=2         ! stay in "normal growth" mode
!             chkmode(i,4)=1
!           endif
!         endif
!        endif
!>
!!broad leaf dcd dry
!!we still use daylength and roottemp to initiate leaf offset,
!!for the pathological cases of dry dcd trees being further
!!away from the equator then we can imagine. other wise leaf
!!loss will occur due to drought anyway.
            case ('BdlDDrTr') ! Broadleaf deciduous drought dry tree
                if( ailcg(i,j).gt.zero.and. &
            &    ((daylngth(i).lt.11.0.and.roottemp(i,j).lt.(11.15+273.16)) &
            &    .or. roottemp(i,j).lt.(roothrsh+273.16)) )then
                    lfstatus(i,j)=3         ! go into "leaf fall" mode
                    chkmode(i,j)=1
                else if(ailcg(i,j).gt.zero.and.ailcg(i,j).lt.lfthrs(i,j))then
                    lfstatus(i,j)=1         ! switch to "max. growth" mode
                    chkmode(i,j)=1
                else if(ailcg(i,j).le.zero) then
                    lfstatus(i,j)=4         ! switch to "no leaves" mode
                    chkmode(i,j)=1
                    pandays(i,j)=0
                else
                    lfstatus(i,j)=2         ! stay in "normal growth" mode
                    chkmode(i,j)=1
                endif

!        if (fcancmx(i,5).gt.0.0) then
!         if(chkmode(i,5).eq.0.and.lfstatus(i,5).eq.2)then
!           if( ailcg(i,5).gt.zero.and.
!      &    ((daylngth(i).lt.11.0.and.roottemp(i,5).lt.(11.15+273.16))
!      &    .or. roottemp(i,5).lt.(roothrsh+273.16)) )then
!             lfstatus(i,5)=3         ! go into "leaf fall" mode
!             chkmode(i,5)=1
!           else if(ailcg(i,5).gt.zero.and.ailcg(i,5).lt.lfthrs(i,5))then
!             lfstatus(i,5)=1         ! switch to "max. growth" mode
!             chkmode(i,5)=1
!           else if(ailcg(i,5).le.zero) then
!             lfstatus(i,5)=4         ! switch to "no leaves" mode
!             chkmode(i,5)=1
!             pandays(i,5)=0
!           else
!             lfstatus(i,5)=2         ! stay in "normal growth" mode
!             chkmode(i,5)=1
!           endif
!         endif
!        endif
!
!>
!!"normal growth" to "fall/harvest" transition for crops is based on
!!specified lai. we harvest if lai of crops reaches a threshold.
!!if lai doesn't reach this threshold (say due to a bad year)
!!we harvest anyway if it starts getting cold, otherwise we don't harvest.
!!
!!
!!    add a branch for icc =12, j6, 7 = shrubs,j 8,9 = crops, j10,11,12
!!    = grasses, YW April 15, 2015 -----------------------------------\
!
!      if (icc==9)                           then

    case ('CropC3  ':'CropC4  ') ! Crops
!      do 340 j = 6,7
       n = sort(j)
!        do 350 i = il1, il2
!         if (fcancmx(i,j).gt.0.0) then
!          if(chkmode(i,j).eq.0.and.lfstatus(i,j).eq.2)then
            if(ailcg(i,j).ge.harvthrs(n))then
              lfstatus(i,j)=3        ! go into "harvest" mode
              chkmode(i,j)=1
              flhrloss(i,j)=gleafmas(i,j)*(1.0/flhrspan(1))
            else if( ailcg(i,j).gt.zero.and.&
     &      colddays(i,2).ge.coldlmt(2) ) then
              lfstatus(i,j)=3        ! go into "harvest" mode
              chkmode(i,j)=1         ! regardless of lai
              flhrloss(i,j)=gleafmas(i,j)*(1.0/flhrspan(1))
            else if(ailcg(i,j).le.zero) then
              lfstatus(i,j)=4        ! switch to "no leaves" mode
              chkmode(i,j)=1         
              pandays(i,j)=0
              flhrloss(i,j)=0.0
            else
              lfstatus(i,j)=2        ! stay in "normal growth" mode
              chkmode(i,j)=1         
            endif
!          endif
!         endif
!350     continue
!340   continue
!>
!!"normal growth" to "max. growth" transition for grasses 
!!
    case ('GrassC3 ':'GrassC4 ') ! Grasses
!      do 370 j = 8,9
!        do 380 i = il1, il2
!         if (fcancmx(i,j).gt.0.0) then
!          if(chkmode(i,j).eq.0.and.lfstatus(i,j).eq.2)then
            if(ailcg(i,j).lt.lfthrs(i,j).and.ailcg(i,j).gt.zero)then  
              lfstatus(i,j)=1        ! switch back to "max. growth" mode
              chkmode(i,j)=1
            else if(ailcg(i,j).le.zero) then
              lfstatus(i,j)=4        ! switch to "no leaves" mode
              chkmode(i,j)=1         
              pandays(i,j)=0
            else
              lfstatus(i,j)=2        ! stay in "normal growth" mode
              chkmode(i,j)=1
            endif
!          endif
!         endif
!380     continue
!370   continue
! FLAG!
!

!      else          !icc != 9 here icc == 12, 6,7 =shrubs, 8,9 = crops,
                    ! 10, 11,12 = grass 
!    EVG-SHRUB      treated the same as needel leaf EVG                 

    case ('BdlEvgSh') ! EVG-SHRUB, treated the same as needel leaf EVG
!      do 385 i = il1, il2
!       if (fcancmx(i,6).gt.0.0) then
!        if(chkmode(i,6).eq.0.and.lfstatus(i,6).eq.2)then

          if(ailcg(i,j).lt.lfthrs(i,j).and.ailcg(i,j).gt.zero)then
            lfstatus(i,j)=1         ! go back to "max. growth" mode
            chkmode(i,j)=1
          else if(ailcg(i,j).le.zero) then
            lfstatus(i,j)=4         ! switch to "no leaves" mode
            chkmode(i,j)=1
            pandays(i,j)=0
          else
            lfstatus(i,j)=2         ! stay in "normal growth" mode
            chkmode(i,j)=1
          endif
!        endif
!       endif
!
!     DCD-SHRUB     treated the same as needle leaf dcd
    case ('BdlDcdSh') ! DCD-SHRUB     treated the same as needle leaf dcd
!       if (fcancmx(i,7).gt.0.0) then
!        if(chkmode(i,7).eq.0.and.lfstatus(i,7).eq.2)then
          if(ailcg(i,j).lt.lfthrs(i,j).and.&
     &      colddays(i,1).ge.coldlmt(1).and.&
     &      ailcg(i,j).gt.zero)then
            lfstatus(i,j)=3         ! go into "leaf fall" mode
            chkmode(i,j)=1
          else if(ailcg(i,j).le.zero) then
            lfstatus(i,j)=4         ! switch to "no leaves" mode
            chkmode(i,j)=1
            pandays(i,j)=0
          else
            lfstatus(i,j)=2         ! stay in "normal growth" mode
            chkmode(i,j)=1
          endif
!        endif
!      endif

    case ('Sedge   ') ! Sedges
        if(ailcg(i,j).lt.lfthrs(i,j).and.ailcg(i,j).gt.zero)then
            lfstatus(i,j)=1        ! switch back to "max. growth" mode
            chkmode(i,j)=1
        else if(ailcg(i,j).le.zero) then
            lfstatus(i,j)=4        ! switch to "no leaves" mode
            chkmode(i,j)=1
            pandays(i,j)=0
        else
            lfstatus(i,j)=2        ! stay in "normal growth" mode
            chkmode(i,j)=1
        endif

!385    continue               !shrubs
      

!       do 390 j = 8,9
!        n = sort(j)
!         do 391 i = il1, il2
!          if (fcancmx(i,j).gt.0.0) then
!           if(chkmode(i,j).eq.0.and.lfstatus(i,j).eq.2)then
!             if(ailcg(i,j).ge.harvthrs(n))then
!               lfstatus(i,j)=3        ! go into "harvest" mode
!               chkmode(i,j)=1
!               flhrloss(i,j)=gleafmas(i,j)*(1.0/flhrspan(1))
!             else if( ailcg(i,j).gt.zero.and.&
!      &      colddays(i,2).ge.coldlmt(2) ) then
!               lfstatus(i,j)=3        ! go into "harvest" mode
!               chkmode(i,j)=1         ! regardless of lai
!               flhrloss(i,j)=gleafmas(i,j)*(1.0/flhrspan(1))
!             else if(ailcg(i,j).le.zero) then
!               lfstatus(i,j)=4        ! switch to "no leaves" mode
!               chkmode(i,j)=1
!               pandays(i,j)=0
!               flhrloss(i,j)=0.0
!             else
!               lfstatus(i,j)=2        ! stay in "normal growth" mode
!               chkmode(i,j)=1
!             endif
!           endif
!          endif
! 391     continue
! 390   continue
! !
! !     "normal growth" to "max. growth" transition for grasses
! !
!       do 395 j = 10,12
!         do 396 i = il1, il2
!          if (fcancmx(i,j).gt.0.0) then
!           if(chkmode(i,j).eq.0.and.lfstatus(i,j).eq.2)then
!             if(ailcg(i,j).lt.lfthrs(i,j).and.ailcg(i,j).gt.zero)then
!               lfstatus(i,j)=1        ! switch back to "max. growth" mode
!               chkmode(i,j)=1
!             else if(ailcg(i,j).le.zero) then
!               lfstatus(i,j)=4        ! switch to "no leaves" mode
!               chkmode(i,j)=1
!               pandays(i,j)=0
!             else
!               lfstatus(i,j)=2        ! stay in "normal growth" mode
!               chkmode(i,j)=1
!             endif
!           endif
!          endif
! 396     continue
! 395   continue
!
!       endif         !if icc==9
      
!          write(90,6990) '680', real(lfstatus)

            end select
        end if
        end if
        end do
320   continue


!>
!!if in "fall/harvest" mode
!!--------------------------
!!
!!grasses and evg trees do not come into this mode, because they want
!!to stay green if possible. this mode is activated for dcd plants and
!!crops. once in this mode dcd trees loose their leaves and crops are
!!harvested. ndl dcd trees keep loosing their leaves at rate determined
!!by cold stress, bdl dcd trees loose their leaves at a specified
!!rate, and crops are harvested over a period of  15 days. dcd trees
!!and crops stay in "leaf fall/harvest" model until all green leaves
!!are gone at which time they switch into "no leaves" mode, and then 
!!wait for the climate to become favourable to go into "max. growth" mode
!!

do 410 i = il1, il2
      do 400 j = 1, icc
        pftkind = pftlist(j)

        if  (pftkind == 'NdlDcdTr' .or. &  ! All deciduous PFTs
          &  pftkind == 'BdlDCoTr' .or. &
          &  pftkind == 'BdlDDrTr' .or. &
          &  pftkind == 'BdlDcdSh' .or. &
          &  pftkind == 'CropC3  ' .or. &
          &  pftkind == 'CropC4  ' ) then

!        if(j.eq.2.or.j.eq.4.or.j.eq.5.or.j.eq.6.or.j.eq.7) !YW for icc=12
!        if(j.eq.2.or.j.eq.4.or.j.eq.5.or.j.eq.7.or.j.eq.8.or. j.eq.9) then  !only dcd trees and crops
!          do 410 i = il1, il2
           if (fcancmx(i,j).gt.0.0 .and. chkmode(i,j).eq.0 .and. lfstatus(i,j).eq.3)then
!            if(chkmode(i,j).eq.0.and.lfstatus(i,j).eq.3)then
              if(ailcg(i,j).le.0.01)then
                lfstatus(i,j)=4            ! go into "no leaves" mode
                chkmode(i,j)=1
                pandays(i,j)=0
                flhrloss(i,j)=0.0
              else
!                if(j.eq.2)then             ! ndl dcd trees
!                if(j.eq. 2 .or. j.eq. 7)then !and dcd shrubs YW April 15 2015
                if (pftkind == 'NdlDcdTr' .or. pftkind == 'BdlDcdSh') then
!         
                  if(pandays(i,j).ge.dayschk(j).and. ta(i).gt.(coldthrs(1)+273.16)) then
                    if(ailcg(i,j).lt.lfthrs(i,j))then
                      lfstatus(i,j)=1      ! go into "max. growth" mode
                      chkmode(i,j)=1
                    else
                      lfstatus(i,j)=2      ! go into "normal growth" mode
                      chkmode(i,j)=1
                    endif
                  else  
                    lfstatus(i,j)=3        ! stay in "fall/harvest" mode
                    chkmode(i,j)=1
                  endif
!                else if(j.eq.4.or.j.eq.5)then        ! bdl dcd trees
                else if (pftkind == 'BdlDCoTr' .or. pftkind == 'BdlDDrTr') then
                  if( (pandays(i,j).ge.dayschk(j)).and.&
     &            ((roottemp(i,4).gt.(roothrsh+273.16)).and.&
     &             (daylngth(i).gt.11.0) ) )then 
                    if(ailcg(i,j).lt.lfthrs(i,j))then
                      lfstatus(i,j)=1      ! go into "max. growth" mode
                      chkmode(i,j)=1
                    else
                      lfstatus(i,j)=2      ! go into "normal growth" mode
                      chkmode(i,j)=1
                    endif
                  else                                 !crops 
                    lfstatus(i,j)=3        ! stay in "fall/harvest" mode
                    chkmode(i,j)=1
                  endif
                else  ! crops
                  lfstatus(i,j)=3          ! stay in "fall/harvest" mode
                  chkmode(i,j)=1
                endif
              endif !ailcg
            endif !fcancmax etc.
           endif !pftlist
400       continue
        !endif
410   continue
!
!         YW July 13, 2015  comment out the testing to reverse to 2.0.4
!       FLAG test done to see impact of no alloc to leaves after 20 days
!       do i = il1, il2
!          j = 2 !needle dcd
!            if (ailcg(i,j).gt.0.0) then
!              if (iday > 192 .and. radl(i) > 0. .and. lfstatus(i,j).ne.
!      &             4) then ! north hemi past summer solstice
!                  lfstatus(i,j) = 3 ! no allocation to leaves permitted
!              else if ((iday < 172 .and. iday > 10) .and. radl(i) < 0.
!      &               .and. lfstatus(i,j).ne. 4)then  ! southern hemi a
!                  lfstatus(i,j) = 3 ! no allocation to leaves permitted
!              end if
!            endif
!          j = 4  ! broad dcd
!            if (ailcg(i,j).gt.0.0) then
!               if (iday > 192  .and.  radl(i) >0. .and. lfstatus(i,j).n
!      &             4) then ! north hemi past summer solstice
!                  lfstatus(i,j) = 3 ! no allocation to leaves permitted
!               else if ((iday < 172 .and. iday > 10) .and. radl(i) < 0.
!      &               .and. lfstatus(i,j).ne. 4) then  ! southern hemi
!                  lfstatus(i,j) = 3 ! no allocation to leaves permitted
!               end if
!            endif
!       end do

!>
!!check that leaf status of all vegetation types in all grid cells has been updated
!!
      do 411 j = 1, icc
        do 412 i = il1, il2
         if (fcancmx(i,j).gt.0.0) then 
          if(chkmode(i,j).eq.0)then
           write(6,2000) i,j
2000       format(' at (i) = (',i3,'), pft=',i2,' lfstatus not updated')
           call xit('phenolgy',-2)
          endif
         endif
412     continue
411   continue
!>
!!------------------------------------------------------------------    
!!
!!having decided leaf status for every pft, we now calculate normal
!!leaf turnover, cold and drought stress mortality, and for bdl dcd 
!!plants we also calculate specified loss rate if they are in "leaf fall"
!!mode, and for crops we calculate harvest loss, if they are in "harvest" mode.
!!
!!all these loss calculations will yield leaf litter in \f$kg c/m^2\f$ for the given day for all pfts 
!!
!!normal leaf turn over
!!
      do 420 j = 1, icc
        n = sort(j)
        do 430 i = il1, il2
         if (fcancmx(i,j).gt.0.0) then 
         nrmlloss(i,j)=gleafmas(i,j)*(1.0-exp(-1.0/(365.0*lfespany(n))))
!          ! FLAG! TEST Dec 102014 JM. Testing the influence of only al
!          ! leaf aging turnover when the lfstatus is >1 (so normal allo
!          ! no alloc to leaves). When lfstatus is 1, it is not applied.
!            if (j == 2 .or. j == 4) then !only deciduous PFTs
!                 if (lfstatus(i,j) .ne. 1) then
!                     nrmlloss(i,j)=gleafmas(i,j)*(1.0-exp(-1.0/
!      &                          (365.0*lfespany(n))))
!                 else
!                     nrmlloss(i,j)=0. ! no loss during leaf out.
!                 end if
!             else ! pfts other than deciduous
!                 nrmlloss(i,j)=gleafmas(i,j)*(1.0-exp(-1.0/
!      &                (365.0*lfespany(n))))
!             end if  !decid/non
         endif   !fcancmx
    
430     continue
420   continue
!>
!!for drought stress related mortality we need field capacity and wilting point soil 
!!moisture contents, which we calculated in allocate subroutine
!!
      do 450 j = 1, ignd
        do 460 i = il1, il2
         if (isand(i,j) .ne. -3) then !for non-bedrock (bedrock keeps the initialization value of 0)
!         estimate (1-drought stress) 
!
          if(thliq(i,j).le.THLW(i,j)) then
            betadrgt(i,j)=0.0
          else if (thliq(i,j).gt.(THLW(i,j)-thice(i,j)) .and.&
     &      (thliq(i,j).lt.(THFC(i,j)-thice(i,j)))) then
            betadrgt(i,j)=(thliq(i,j)-(THLW(i,j)-thice(i,j)))
            betadrgt(i,j)=betadrgt(i,j)/(THFC(i,j)-thice(i,j)&
     &                   -THLW(i,j))
          else 
            betadrgt(i,j)=1.0
          endif          
          betadrgt(i,j)=max(0.0, min(1.0,betadrgt(i,j)))
          end if
460     continue
450   continue
!>
!!estimate drought stress term averaged over the rooting depth for each pft
!!
      do 480 j = 1, icc
        n = sort(j)

        do 490 i = il1, il2
         if (fcancmx(i,j).gt.0.0) then 
          drgtstrs(i,j) =  sum((1.0-betadrgt(i,:))*rmatctem(i,j,:))
          drgtstrs(i,j) = drgtstrs(i,j) / sum(rmatctem(i,j,:))
!    ---------YW May 04, 2015 ---------------------------------------/
!          drgtstrs(i,j)=max(0.0, min(1.0,drgtstrs(i,j))) !FLAG why is this commented out? JM Nov 2016.
       
!>         using this drought stress term and our two vegetation-dependent
!!         parameters we find leaf loss rate associated with drought

!!         drought related leaf loss rate
          drgtlsrt(i,j)=drlsrtmx(n)*(drgtstrs(i,j)**drgta(n))

!>estimate leaf loss in \f$kg c/m^2\f$ due to drought stress
          drgtloss(i,j)=gleafmas(i,j)*( 1.0-exp(-drgtlsrt(i,j)) )

!>
!>similar to drgtstrs we find coldstrs for each pft. we assume that
!!max. cold stress related leaf loss occurs when temperature is 5 c
!!or more below pft's threshold
!!
          if(ta(i).le.(lwrthrsh(n)-5.0+273.16))then
            coldstrs(i,j)=1.0
          else if(ta(i).gt.(lwrthrsh(n)-5.0+273.16).and.&
     &    ta(i).lt.(lwrthrsh(n)+273.16))then
            coldstrs(i,j)=1.0-((ta(i)-(lwrthrsh(n)-5.0+273.16))/(5.0))  
          else 
            coldstrs(i,j)=0.0
          endif
          coldstrs(i,j)=max(0.0, min(1.0,coldstrs(i,j)))
!>
!>using this cold stress term and our two vegetation-dependent
!>parameters we find leaf loss rate associated with cold
!>cold related leaf loss rate
          coldlsrt(i,j)=cdlsrtmx(n)*(coldstrs(i,j)**colda(n))
!>
!!estimate leaf loss in \f$kg c/m^2\f$ due to cold stress
          coldloss(i,j)=gleafmas(i,j)*( 1.0-exp(-coldlsrt(i,j)) )

         endif
490     continue 
480   continue 
!>
!!now that we have all types of leaf losses (due to normal turnover,
!!cold and drought stress, and fall/harvest) we take the losses
!!for grasses and use those to turn live green grass into dead
!!brown grass. we then find the leaf litter from the brown grass
!!which will then go into the litter pool.
!!
!     YW May 11, 2015 add ICC = 12  
!      if (icc==9) then
!      do 620 j = 8,9

    do 620 j = 1,icc
     pftkind = pftlist(j)
     if (pftkind == 'GrassC3 ' .or. &
         pftkind == 'GrassC4 ' .or. &
         pftkind == 'Sedge   ') then
        n = sort(j)
         do 630 i = il1, il2
         if (fcancmx(i,j).gt.0.0) then 
          gleafmas(i,j)   = gleafmas(i,j)-nrmlloss(i,j)-drgtloss(i,j)-  &
     &                      coldloss(i,j)
          if( gleafmas(i,j).lt.0.0) then
            bleafmas(i,j) = bleafmas(i,j)+nrmlloss(i,j)+drgtloss(i,j)+  &
     &                      coldloss(i,j)+gleafmas(i,j)
            gleafmas(i,j)=0.0
          else
            bleafmas(i,j) = bleafmas(i,j)+nrmlloss(i,j)+drgtloss(i,j)+  &
     &                      coldloss(i,j)
          endif
          nrmlloss(i,j) = 0.0
          drgtloss(i,j) = 0.0
          coldloss(i,j) = 0.0
!>we assume life span of brown grass is 10% that of green grass
!!but this is an adjustable parameter.
          nrmlloss(i,j) = bleafmas(i,j)*&
     &      (1.0-exp(-1.0/(0.10*365.0*lfespany(n)))) 
         endif
630     continue
    end if
620   continue

!       else          !if icc==12
!       do 640 j = 10,12
!         n = sort(j)
!         do 645 i = il1, il2
!          if (fcancmx(i,j).gt.0.0) then
!           gleafmas(i,j)   = gleafmas(i,j)-nrmlloss(i,j)-drgtloss(i,j)-  &
!      &                      coldloss(i,j)
!           if( gleafmas(i,j).lt.0.0) then
!             bleafmas(i,j) = bleafmas(i,j)+nrmlloss(i,j)+drgtloss(i,j)+  &
!      &                      coldloss(i,j)+gleafmas(i,j)
!             gleafmas(i,j)=0.0
!           else
!             bleafmas(i,j) = bleafmas(i,j)+nrmlloss(i,j)+drgtloss(i,j)+  &
!      &                      coldloss(i,j)
!           endif
!           nrmlloss(i,j) = 0.0
!           drgtloss(i,j) = 0.0
!           coldloss(i,j) = 0.0
! !         we assume life span of brown grass is 10% that of green grass
! !         but this is an adjustable parameter.
!           nrmlloss(i,j) = bleafmas(i,j)*&
!      &      (1.0-exp(-1.0/(0.10*365.0*lfespany(n))))
!          endif
! 645     continue
! 640   continue
!
!       endif         ! icc ==9

!>
!!combine nrmlloss, drgtloss, and coldloss together to get total
!!leaf litter generated, which is what we use to update gleafmass,
!!except for grasses, for which we have already updated gleafmass.
!!

      do 650 j = 1, icc
        do 660 i = il1, il2
         if (fcancmx(i,j).gt.0.0) then
          leaflitr(i,j) = nrmlloss(i,j)+drgtloss(i,j)+coldloss(i,j)+&
     &                    flhrloss(i,j)
         endif
660     continue
650   continue
!
!      write(90,6991)  iday, nrmlloss(1,1),drgtloss(1,1),
!     1    coldloss(1,1),flhrloss(1,1),drgtstrs(1,1),betadrgt(1,1),
!     2    thliq(1,1),thice(1,1),THLW(1,1),THFC(1,1)
!6991      format(I5, 12E10.3)
      return
      end

