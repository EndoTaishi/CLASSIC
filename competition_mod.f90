module competition_scheme

! Central module for all competition scheme-related operations

! J. Melton. Jun 22, 2013

implicit none

! Subroutines contained in this module:
public  :: bioclim
public  :: existence
public  :: competition

contains

!-------------------------------------------------------------------------------------------------------------

subroutine  bioclim (   iday,        ta,   precip,   netrad, &
                              il1,       il2,      nilg, &
                            tcurm,  srpcuryr, dftcuryr,  inibioclim, &
                           tmonth,  anpcpcur,  anpecur,   gdd5cur, &
                         surmncur,  defmncur, srplscur,  defctcur, &
                           twarmm,    tcoldm,     gdd5,  aridity, &
                         srplsmon,  defctmon, anndefct, annsrpls, &
                           annpcp,  anpotevp, dry_season_length)    

!               Canadian Terrestrial Ecosystem Model (CTEM)
!                Bioclimatic Parameters Estimation Subroutine 
!
!     25  Jun 2013  - Convert to f90.
!     J. Melton
!
!     22  Nov 2012  - Calling this version 1.1 since a fair bit of ctem
!     V. Arora        subroutines were changed for compatibility with class
!                     version 3.6 including the capability to run ctem in
!                     mosaic/tile version along with class.
!
!     25  May 2004  - This subroutine calculates the bioclimatic
!     V. Arora        parameters that are required for determining
!                     existence of pfts. the bioclimatic parameters 
!                     are the mean monthly temperature of the warmest
!                     and the coldest months, growing degree days 
!                     above 5 c, annual precipitation and potential
!                     evaporation and some aridity parameters that are
!                     function of potential evaporation and precipitation.
!
!                     In addition, all these parameters are updated
!                     in an e-folding sense at some specified time
!                     scale.
!
!                     Note that this subroutine is only necessary when
!                     competition is switched on.

use ctem_params, only : zero, monthdays, monthend

implicit none

! arguments

integer, intent(in) :: iday      ! day of the year
integer, intent(in) :: nilg       ! no. of grid cells in latitude circle (this is passed in as either ilg or nlat depending on mos/comp)
integer, intent(in) :: il1       ! il1=1
integer, intent(in) :: il2       ! il2=nilg
real, dimension(nilg), intent(in)    :: ta        ! mean daily temperature, k
real, dimension(nilg), intent(in)    :: precip    ! daily precipitation (mm/day)
real, dimension(nilg), intent(in)    :: netrad    ! daily net radiation (w/m2)

logical, intent(inout) :: inibioclim  ! switch telling if bioclimatic parameters are being
                                    ! initialized from scratch (false) or being initialized
                                    ! from some spun up values(true).
real, dimension(nilg), intent(inout) :: tcurm     ! temperature of the current month (c)
real, dimension(nilg), intent(inout) :: srpcuryr  ! water surplus for the current year
real, dimension(nilg), intent(inout) :: dftcuryr  ! water deficit for the current year
real, dimension(12,nilg), intent(inout) :: tmonth    ! monthly temperatures
real, dimension(nilg), intent(inout) :: anpcpcur  ! annual precipitation for current year (mm)
real, dimension(nilg), intent(inout) :: anpecur   ! annual potential evaporation for current year (mm)
real, dimension(nilg), intent(inout) :: gdd5cur   ! growing degree days above 5 c for current year
integer, dimension(nilg), intent(inout) :: surmncur  ! number of months with surplus water for current year
integer, dimension(nilg), intent(inout) :: defmncur  ! number of months with water deficit for current year
real, dimension(nilg), intent(inout) :: srplscur  ! water surplus for the current month
real, dimension(nilg), intent(inout) :: defctcur  ! water deficit for the current month

! the following are running averages in an e-folding sense

real, dimension(nilg), intent(inout) :: twarmm    ! temperature of the warmest month (c)
real, dimension(nilg), intent(inout) :: tcoldm    ! temperature of the coldest month (c)
real, dimension(nilg), intent(inout) :: gdd5      ! growing degree days above 5 c
real, dimension(nilg), intent(inout) :: aridity   ! aridity index, ratio of potential evaporation to precipitation
real, dimension(nilg), intent(inout) :: srplsmon  ! number of months in a year with surplus water i.e.
                                                  !  precipitation more than potential evaporation
real, dimension(nilg), intent(inout) :: defctmon  ! number of months in a year with water deficit i.e.
                                                  ! precipitation less than potential evaporation
real, dimension(nilg), intent(inout) :: anndefct  ! annual water deficit (mm) 
real, dimension(nilg), intent(inout) :: annsrpls  ! annual water surplus (mm)
real, dimension(nilg), intent(inout) :: annpcp    ! annual precipitation (mm)
real, dimension(nilg), intent(inout) :: anpotevp  ! annual potential evaporation (mm) 
real, dimension(nilg), intent(inout) :: dry_season_length ! annual maximum dry month length   !Rudra
 
! local variables
real, dimension(nilg) :: tccuryr
real, dimension(nilg) :: twcuryr
real, dimension(nilg) :: aridcur
real :: wtrbal
integer :: month, atmonthend, temp, nmax, i, j, k, curmonth, m, n, l
integer, save, dimension(:,:), allocatable :: wet_dry_mon_index   !Rudra
integer, save, dimension(:,:), allocatable :: wet_dry_mon_index2   !Rudra
real, dimension(:), allocatable, save :: dry_season_length_curyr !current year's maximum dry month length 

! local parameters
real, parameter :: eftime = 25.00 ! e-folding time scale for updating bioclimatic parameters (years)
real, parameter :: factor=exp(-1.0/eftime) !faster to calculate this only at compile time.

!     ---------------------------------------------------------------

!     initializations

      if(iday.eq.1)then

        ! Allocate the arrays to find the length of dry season
          allocate(wet_dry_mon_index(nilg,12))
          allocate(wet_dry_mon_index2(nilg,24))
          allocate(dry_season_length_curyr(nilg))

        do 100 i = il1, il2
          gdd5cur(i)=0.0    ! gdd5 for the current year
          anpcpcur(i)=0.0   ! annual precip. for the current year
          anpecur(i)=0.0    ! annual potential evap for the current year
          aridcur(i)=100.0  ! aridity index for the current year
          surmncur(i)=0     ! months with surplus water for current year
          defmncur(i)=0     ! months with water deficit for current year
          srpcuryr(i)=0.0   ! current year's water surplus
          dftcuryr(i)=0.0   ! current year's water deficit
          tcurm(i)=0.0      ! temperature of current month
          srplscur(i)=0.0   ! current month's water surplus
          defctcur(i)=0.0   ! current month's water deficit
          dry_season_length_curyr(i) = 0.   !current year's maximum dry month length    !Rudra

         do month = 1,12
          tmonth(month,i)=0.0
         end do
100     continue      
      endif

!     Find current month

      curmonth=0
      do 220 k = 2, 13
        if(iday.ge.(monthend(k-1)+1).and.iday.le.monthend(k))then
          curmonth=k-1
        endif
220   continue
        
      if(curmonth.eq.0)then
        call xit('bioclim',-1)
      endif

!     Find if we are at end of month or not
      atmonthend=0
      if (iday.eq.monthend(curmonth+1)) then
        atmonthend=1
      endif

!     Update monthly temperature for the current month, and other
!     variables. at the end of the month we will have average of 
!     all daily temperatures for the current month.

      do 240 i = il1, il2
          tcurm(i)=tcurm(i)+(ta(i)-273.16)*(1.0/real(monthdays(curmonth)))
          gdd5cur(i)=gdd5cur(i)+max(0.0, (ta(i)-273.16-5.0))
          anpcpcur(i)=anpcpcur(i) + precip(i)
!         net radiation (W/m2) x 12.87 = potential evap (mm)
          anpecur(i)=anpecur(i) + netrad(i)*12.87*(1.0/365.0)
          wtrbal=precip(i)-(netrad(i)*12.87*(1.0/365.0))
          if(wtrbal.ge.0.0)then
            srplscur(i)=srplscur(i)+wtrbal
          else if(wtrbal.lt.0.0)then
            defctcur(i)=defctcur(i)+abs(wtrbal)
          endif
240   continue

!     If its the end of the month then store the monthly temperature 
!     and set tcurm equal to zero. also check if this month had water
!     deficit or surplus

      do 250 i = il1, il2
        if(atmonthend.eq.1)then
          tmonth(curmonth,i)=tcurm(i)
          if( srplscur(i).ge.defctcur(i) )then
            surmncur(i) = surmncur(i) + 1
            wet_dry_mon_index(i,curmonth) = 1    !Rudra
          else if(srplscur(i).lt.defctcur(i) )then
            defmncur(i) = defmncur(i) + 1
            wet_dry_mon_index(i,curmonth) = -1   !Rudra
          endif
          srpcuryr(i)=srpcuryr(i)+srplscur(i)
          dftcuryr(i)=dftcuryr(i)+defctcur(i)

          tcurm(i)=0.0    ! temperature of current month
          srplscur(i)=0.0 ! current month's water surplus
          defctcur(i)=0.0 ! current month's water deficit

        endif

        if(iday.eq.365)then
          twcuryr(i)=-9000.0
          tccuryr(i)=9000.0
          if(anpcpcur(i).gt.zero)then
            aridcur(i)=anpecur(i)/anpcpcur(i)
          else
            aridcur(i)=100.0
          endif
!................................../Rudra
                             
            do j = 1,12          !this loop double up the size of the "wet_dry_mon_index" matrix 
                do k = 1,2
                   m = (k-1)*12 + j
                   wet_dry_mon_index2(i,m) = wet_dry_mon_index(i,j)
                end do
            end do
!...................MAXIMUM LENGTH OF DRY MONTH.................................

            n = 0       !number of dry month
            nmax = 0    !maximum length of dry month
               do l = 1, 24
                   temp = wet_dry_mon_index2(i,l) 
                   if(temp.eq.-1)then
                      n = n+1
                      nmax = max(nmax, n) 
                   else 
                      n = 0
                   end if 
                end do 
            nmax = min(nmax, 12)
            dry_season_length_curyr(i) = real(nmax)
!..................................\Rudra
        endif     !iday=365

250   continue

!     If its the end of year, then find the temperature of the warmest
!     and the coldest month
      if(iday.eq.365)then
          do 270 i = il1, il2
              twcuryr(i)=maxval(tmonth(:,i))
              tccuryr(i)=minval(tmonth(:,i))
270       continue

!       Update long term moving average of bioclimatic parameters in an 
!       e-folding sense
        do 280 i = il1, il2          
          if(.not. inibioclim)then
            twarmm(i)=twcuryr(i)
            tcoldm(i)=tccuryr(i)
            gdd5(i)=gdd5cur(i)
            aridity(i)=aridcur(i)
            srplsmon(i)=real(surmncur(i))
            defctmon(i)=real(defmncur(i))
            annsrpls(i)=srpcuryr(i)
            anndefct(i)=dftcuryr(i)
            annpcp(i)=anpcpcur(i)
            anpotevp(i)=anpecur(i)
            dry_season_length(i)=dry_season_length_curyr(i)    !Rudra
            inibioclim=.true.
          else
            twarmm(i)=twarmm(i)*factor + twcuryr(i)*(1.0-factor)
            tcoldm(i)=tcoldm(i)*factor + tccuryr(i)*(1.0-factor)
            gdd5(i)  =gdd5(i)*factor + gdd5cur(i)*(1.0-factor)
            aridity(i)=aridity(i)*factor + aridcur(i)*(1.0-factor)
            srplsmon(i)=srplsmon(i)*factor + real(surmncur(i))*(1.0-factor)
            defctmon(i)=defctmon(i)*factor + real(defmncur(i))*(1.0-factor)
            annsrpls(i)=annsrpls(i)*factor + srpcuryr(i)*(1.0-factor)
            anndefct(i)=anndefct(i)*factor + dftcuryr(i)*(1.0-factor)
            annpcp(i)=annpcp(i)*factor + anpcpcur(i)*(1.0-factor)
            anpotevp(i)=anpotevp(i)*factor + anpecur(i)*(1.0-factor)
            dry_season_length(i)=dry_season_length(i)*factor + dry_season_length_curyr(i)*(1.0-factor)    !Rudra
          endif
280     continue

         ! Deallocate the arrays for dry season length
          deallocate(wet_dry_mon_index)
          deallocate(wet_dry_mon_index2)
          deallocate(dry_season_length_curyr)

      endif

      return

end subroutine bioclim

!-------------------------------------------------------------------------------------------------------------

subroutine  existence(  iday,       il1,      il2,      nilg, &
                             sort,  nol2pfts,                 &
                           twarmm,    tcoldm,     gdd5,  aridity, &
                         srplsmon,  defctmon, anndefct, annsrpls, &
                           annpcp,  anpotevp,pftexist,dry_season_length) 

!               Canadian Terrestrial Ecosystem Model (CTEM)
!                          PFT Existence Subroutine 
!
!     27  Jan 2014  - Moved parameters to global file (ctem_params.f90)
!     J. Melton
!
!     26  Nov 2013  - Update parameters for global off-line runs
!     J. Melton
!
!     25  Jun 2013  - Convert to f90, incorporate modules, and into larger module.
!     J. Melton         

!     27  May 2004  - This subroutine calculates the existence of
!     V. Arora        pfts in grid cells based on a set of bioclimatic
!                     parameters that are estimated in a running average
!                     sense using some specified timescale. 
!
!                     If long term averaged bioclimatic parameters
!                     indicate non-existence of a pft then an additional
!                     mortality term kicks in the competition eqns. Also
!                     while a pft may be able to exist in a grid cell it
!                     may be excluded by competition from other pfts.
!
!                -->  Note that since the fractional coverage of c3 and
!                -->  c4 crops is going to be prescribed, the model
!                -->  assumes that these pfts can always exist. But, of
!                     course the prescribed fractional coverage of these 
!                     pfts will decide if they are present in a grid cell or
!                     not.

use ctem_params, only : zero, kk, icc, ican, tcoldmin, tcoldmax, twarmmax, &
                        gdd5lmt, aridlmt, dryseasonlmt

implicit none

! arguments

integer, intent(in) :: iday      ! day of the year
integer, intent(in) :: nilg       ! no. of grid cells in latitude circle (this is passed in as either ilg or nlat depending on mos/comp)
integer, intent(in) :: il1       ! il1=1
integer, intent(in) :: il2       ! il2=nilg

integer, dimension(icc), intent(in) :: sort ! index for correspondence between 9 ctem pfts and
                                            ! size 12 of parameter vectors
integer, dimension(ican), intent(in) :: nol2pfts ! number of level 2 ctem pfts
real, dimension(nilg), intent(in) :: twarmm    ! temperature of the warmest month (c)
real, dimension(nilg), intent(in) :: tcoldm    ! temperature of the coldest month (c)
real, dimension(nilg), intent(in) :: gdd5      ! growing degree days above 5 c
real, dimension(nilg), intent(in) :: aridity   ! aridity index, ratio of potential evaporation to precipitation
real, dimension(nilg), intent(in) :: srplsmon  ! number of months in a year with surplus water i.e.
                                                  !  precipitation more than potential evaporation
real, dimension(nilg), intent(in) :: defctmon  ! number of months in a year with water deficit i.e.
                                                  ! precipitation less than potential evaporation
real, dimension(nilg), intent(in) :: anndefct  ! annual water deficit (mm) 
real, dimension(nilg), intent(in) :: annsrpls  ! annual water surplus (mm)
real, dimension(nilg), intent(in) :: annpcp    ! annual precipitation (mm)
real, dimension(nilg), intent(in) :: anpotevp  ! annual potential evaporation (mm)
real, dimension(nilg), intent(in) :: dry_season_length ! length of dry season (months)

logical, dimension(nilg,icc), intent(out) :: pftexist(nilg,icc) !binary array indicating pfts exist (=1) or not (=0)

! local variables
integer :: i,j

! ----------------------------------------------------------------------
!     Constants and parameters are located in ctem_params.f90
! ----------------------------------------------------------------------

!     go through all grid cells and based on bioclimatic parameters
!     decide if a given pft should exist or not. 

      do 100 i = il1, il2

!       needleleaf evergreen
        j=1
        if(tcoldm(i).le.tcoldmax(sort(j)).and.gdd5(i).ge.gdd5lmt(sort(j)) &
           .and.dry_season_length(i).le.dryseasonlmt(sort(j)))then
           pftexist(i,j)=.true.
        else
           pftexist(i,j)=.false.
        endif

!       needleleaf deciduous
        j=2
        if(tcoldm(i).le.tcoldmax(sort(j)).and.twarmm(i).le.twarmmax(sort(j)).and. &
           gdd5(i).ge.gdd5lmt(sort(j)))then
           pftexist(i,j)=.true.
        else
           pftexist(i,j)=.false.
        endif

!       broadleaf evergreen
        j=3
        if(tcoldm(i).ge.tcoldmin(sort(j)).and. &
           gdd5(i).ge.gdd5lmt(sort(j)))then
           pftexist(i,j)=.true.
        else
           pftexist(i,j)=.false.
        endif

!       broadleaf deciduous cold  (see note below pft 5 too)
        j=4
        if(tcoldm(i).le.tcoldmax(sort(j)).and. &
           gdd5(i).ge.gdd5lmt(sort(j)).and. tcoldm(i).ge.tcoldmin &  !TEST FLAG added in this new constraint. JM Dec 17 2013.
           (sort(j)))then
           pftexist(i,j)=.true.
        else
           pftexist(i,j)=.false.
        endif

!       broadleaf deciduous dry
        j=5
        if(tcoldm(i).ge.tcoldmin(sort(j)).and. aridity(i).ge.aridlmt(sort(j)) &
           .and.dry_season_length(i).ge.dryseasonlmt(sort(j)))then
           pftexist(i,j)=.true.
!       We don't want both broadleaf species co-existing so if it has PFT 5
!       remove PFT 4. TEST JM flag Dec 6 2013.
           pftexist(i,j-1)=.false.
        else
           pftexist(i,j)=.false.
        endif

!       c3 and c4 crops
        pftexist(i,6)=.true.
        pftexist(i,7)=.true.

!       c3 grass
        j=8
        if(tcoldm(i).le.tcoldmax(sort(j)))then
           pftexist(i,j)=.true.
        else
           pftexist(i,j)=.false.
        endif

!       c4 grass
        j=9
        if(tcoldm(i).ge.tcoldmin(sort(j)))then
           pftexist(i,j)=.true.
        else
           pftexist(i,j)=.false.
        endif

100   continue

      return

end subroutine existence

!-------------------------------------------------------------------------------------------------------------

subroutine competition(  iday,      il1,       il2,      nilg, &
                          nol2pfts,   nppveg,   dofire, &
                          pftexist,  geremort, intrmort, &
                          gleafmas, bleafmas,  stemmass, rootmass, &
                          litrmass, soilcmas,  grclarea,   lambda, &
                           burnvegf,     sort, pstemmass, &
                           pgleafmass, &
                           fcancmx,   fcanmx,  vgbiomas, gavgltms, &
                          gavgscms,  bmasveg,   &
                          add2allo,        colrate,        mortrate) 

!               Canadian Terrestrial Ecosystem Model (CTEM) V1.1
!                          PFT Competition Subroutine 

!     26  Mar 2014  - Move disturbance adjustments back via a subroutine called
!     J. Melton       burntobare
!
!     20  Feb 2014  - Move adjustments due to disturbance out of here and into
!     J. Melton       disturbance subroutine.
!
!     27  Jan 2014  - Moved parameters to global file (ctem_params.f90)
!     J. Melton

!     25  Jun 2013  - Convert to f90, incorporate modules, and into larger module.
!     J. Melton
 
!     17  Oct 2012  - Adapt subroutine to any number of crops or grass 
!     J. Melton       pfts
!
!     1   Oct 2012  - Update subroutine and implement for running with
!     Y. Peng         mosaic version of class 3.6
!
!     27  May 2004  - This subroutine calculates the competition between
!     V. Arora        pFTs based on Lotka-Volterra eqns. ot its modified
!                     forms. either option may be used.
!
!                     PFTs that may exist are allowed to compete for
!                     available space in a grid cell. pfts that can't
!                     exist based on long term bioclimatic parameters
!                     are slowly killed by increasing their mortality.               
!

use ctem_params, only : zero, kk, numcrops, numgrass, numtreepfts, &
                        icc, ican, deltat, iccp1, seed, bio2sap, bioclimrt, &
                        tolrance, crop, grass, grass_ind

use disturbance_scheme, only : burntobare


implicit none

! arguments

integer, intent(in) :: iday      ! day of the year
integer, intent(in) :: nilg       ! no. of grid cells in latitude circle (this is passed in as either ilg or nlat depending on mos/comp)
integer, intent(in) :: il1       ! il1=1
integer, intent(in) :: il2       ! il2=nilg
logical ,intent(in) :: dofire    ! if true then we have disturbance on.
real,  dimension(nilg),  intent(in) :: grclarea  ! grid cell area, km^2
integer, dimension(icc), intent(in) :: sort ! index for correspondence between 9 ctem pfts and
                                            ! size 12 of parameter vectors
integer, dimension(ican), intent(in) :: nol2pfts ! number of level 2 ctem pfts
logical, dimension(nilg,icc), intent(in) :: pftexist(nilg,icc) !indicating pfts exist (T) or not (F)
real, dimension(nilg,icc), intent(inout) :: nppveg ! npp for each pft type /m2 of vegetated area u-mol co2-c/m2.sec
real, dimension(nilg,icc), intent(in) :: geremort  ! growth related mortality (1/day)
real, dimension(nilg,icc), intent(in) :: intrmort  ! intrinsic (age related) mortality (1/day)
real, dimension(nilg,icc), intent(in) :: lambda    ! fraction of npp that is used for spatial expansion
real, dimension(nilg,icc), intent(inout) :: bmasveg   ! total (gleaf + stem + root) biomass for each ctem pft, kg c/m2
real, dimension(nilg,icc), intent(in) :: burnvegf   ! fractional areas burned, for 9 ctem pfts

real, dimension(nilg,icc), intent(inout) :: gleafmas  ! green leaf mass for each of the 9 ctem pfts, kg c/m2
real, dimension(nilg,icc), intent(inout) :: bleafmas  ! brown leaf mass for each of the 9 ctem pfts, kg c/m2
real, dimension(nilg,icc), intent(inout) :: stemmass  ! stem mass for each of the 9 ctem pfts, kg c/m2
real, dimension(nilg,icc), intent(inout) :: rootmass  ! root mass for each of the 9 ctem pfts, kg c/m2
real, dimension(nilg,iccp1), intent(inout) :: litrmass  ! litter mass for each of the 9 ctem pfts + bare, kg c/m2
real, dimension(nilg,iccp1), intent(inout) :: soilcmas  ! soil carbon mass for each of the 9 ctem pfts + bare, kg c/m2
real, dimension(nilg,icc), intent(inout) :: fcancmx  ! fractional coverage of ctem's 9 pfts
real, dimension(nilg,ican), intent(inout)  :: fcanmx   ! fractional coverage of class' 4 pfts
real, dimension(nilg),     intent(inout) :: vgbiomas ! grid averaged vegetation biomass, kg c/m2
real, dimension(nilg),     intent(inout) :: gavgltms ! grid averaged litter mass, kg c/m2
real, dimension(nilg),     intent(inout) :: gavgscms ! grid averaged soil c mass, kg c/m2

real, dimension(nilg,icc), intent(out) :: add2allo   ! npp kg c/m2.day that is used for expansion and
                                                    ! subsequently allocated to leaves, stem, and root via 
                                                    ! the allocation part of the model.
real, dimension(nilg,icc), intent(out) :: colrate    ! colonization rate (1/day)    
real, dimension(nilg,icc), intent(out) :: mortrate   ! mortality rate
real, dimension(nilg,icc), intent(in) :: pstemmass   ! stem mass from previous timestep, is value before fire. used by burntobare subroutine
real, dimension(nilg,icc), intent(in) :: pgleafmass   ! root mass from previous timestep, is value before fire. used by burntobare subroutine

! local variables

integer :: i, j
integer :: n, k, k1, k2, l, a, b, g
integer :: sdfracin
integer, dimension(nilg) :: t1
integer, dimension(icc-numcrops) :: inirank
integer, dimension(nilg,icc-numcrops) :: rank
integer, dimension(nilg,icc-numcrops) :: exist1
integer, dimension(nilg,icc-numcrops) :: useexist
integer, dimension(nilg,icc) :: fraciord
integer, dimension(nilg) :: bareiord  

real :: befrmass, aftrmass
real :: sum1, sum2, sum3,term,sum4 
real :: colmult
real, dimension(nilg,icc) :: mrtboclm
real, dimension(nilg,icc) :: usenppvg
real, dimension(nilg) :: temp
real, dimension(nilg,icc-numcrops) :: usefrac, usec, usem
real, dimension(nilg,icc-numcrops) :: frac
real, dimension(nilg,icc-numcrops) :: c1
real, dimension(nilg,icc-numcrops) :: m1
real, dimension(nilg,icc-numcrops) :: term2, term3, term4, colterm, deathterm
real, dimension(nilg,icc-numcrops) :: delfrac
real, dimension(nilg) :: cropfrac, vegfrac    
real, dimension(nilg,icc) :: chngfrac
real, dimension(nilg,icc) :: expnterm, mortterm  
real, dimension(nilg,icc) :: pglfmass  
real, dimension(nilg,icc) :: pblfmass
real, dimension(nilg,icc) :: protmass
real, dimension(nilg,icc) :: pstmmass 
real, dimension(nilg,icc) :: pfcancmx
real, dimension(nilg) :: mincfrac  
real, dimension(nilg,icc) :: pbiomasvg, biomasvg
real, dimension(nilg,icc) :: putaside
real, dimension(nilg,icc) :: nppvegar 
real, dimension(nilg,iccp1) :: pltrmass
real, dimension(nilg,iccp1) :: psocmass
real, dimension(nilg,iccp1) :: deadmass 
real, dimension(nilg,iccp1) :: pdeadmas
real, dimension(nilg) :: barefrac   
real, dimension(nilg,icc) :: usebmsvg  
real, dimension(nilg,iccp1) ::ownsolc, ownlitr
real, dimension(nilg,icc) :: baresolc
real, dimension(nilg,icc) :: barelitr, baresoilc                
real, dimension(nilg,iccp1) :: incrlitr, incrsolc
real, dimension(nilg) :: pvgbioms
real, dimension(nilg) :: pgavltms
real, dimension(nilg) :: pgavscms
real, dimension(nilg) :: add2dead
real, dimension(nilg) :: gavgputa
real, dimension(nilg) :: gavgnpp       
real, dimension(nilg) :: pbarefra
real, dimension(nilg) :: grsumlit, grsumsoc

!     ---------------------------------------------------------------
!     Constants and parameters are located in ctem_params.f90
!     ---------------------------------------------------------------

! Model switches:

! set desired model to be used to .true. and all other to .false.
!                ** ONLY ONE can be true. **
logical, parameter :: lotvol=.false. ! original lotka-volterra eqns.
logical, parameter :: arora =.true.  ! modified form of lv eqns with f missing
logical, parameter :: boer  =.false. ! modified form of lv eqns with f missing and a modified self-thinning term

!     ---------------------------------------------------------------

      if(icc.ne.9)                      call xit('competition',-1)
      if(ican.ne.4)                       call xit('competition',-2)

!    set competition parameters according to the model chosen
!
      if (lotvol) then
        a=1 ! alpha. this is the b in the arora & boer (2006) paper
        b=1 ! beta
        g=0 ! gamma
        colmult=4.00  ! multiplier for colonization rate
       ! sdfracin=1   ! seed fraction index ! All seed values set to 'seed'
      else if (arora) then
        a=0 ! alpha
        b=1 ! beta
        g=0 ! gamma
        colmult=1.00  ! multiplier for colonization rate
      !  sdfracin=2   ! seed fraction index ! All seed values set to 'seed'
      else if (boer) then
        a=0 ! alpha
        b=1 ! beta
        g=1 ! gamma
        colmult=1.00 ! multiplier for colonization rate
        !sdfracin=2   ! seed fraction index ! All seed values set to 'seed'
      else
        write(6,*)'fool! choose competition model properly'
        call xit('competition',-4)
      endif

!     ---------------------------------------------------------------

!   First, let's adjust the fractions if fire is turned on.

    if (dofire) then

        call burntobare(il1, il2, vgbiomas, gavgltms, gavgscms,fcancmx, burnvegf, stemmass, &
                      rootmass, gleafmas, bleafmas, litrmass, soilcmas, pstemmass, pgleafmass, &
                      nppveg)

      ! Since the biomass pools could have changed, update bmasveg.
      do 190 i = il1, il2
        do 195 j = 1, icc
         if (fcancmx(i,j).gt.0.0) then
          bmasveg(i,j)=gleafmas(i,j)+stemmass(i,j)+rootmass(i,j)
         endif
195     continue
190   continue

    end if

!     Do our usual initialization

      do 150 j = 1, icc
        do 160 i = il1, il2
          colrate(i,j)=0.0        ! colonization rate
          mortrate(i,j)=0.0        ! mortality rate
          mrtboclm(i,j)=0.0 ! mortality rate if long-term bioclimatic 
!                           ! conditions become unfavourable
          usenppvg(i,j)=0.0
          chngfrac(i,j)=0.0
          expnterm(i,j)=0.0
          mortterm(i,j)=0.0
          add2allo(i,j)=0.0
          pglfmass(i,j)=gleafmas(i,j) ! save all biomasses before making
          pblfmass(i,j)=bleafmas(i,j) ! changes so that we can make sure
          protmass(i,j)=rootmass(i,j) ! mass balance is preserved.
          pstmmass(i,j)=stemmass(i,j)
          pltrmass(i,j)=litrmass(i,j)
          psocmass(i,j)=soilcmas(i,j)
          pfcancmx(i,j)=fcancmx(i,j)
          biomasvg(i,j)=0.0
          pbiomasvg(i,j)=0.0
          nppvegar(i,j)=0.0
          deadmass(i,j)=0.0
          pdeadmas(i,j)=0.0
          barelitr(i,j)=0.0    ! kg c of litter added to bare fraction
          baresolc(i,j)=0.0    ! and same for soil c
          fraciord(i,j)=0
          incrlitr(i,j)=0.0
          incrsolc(i,j)=0.0
          ownlitr(i,j)=0.0
          ownsolc(i,j)=0.0
160     continue
150   continue

      do 170 i = il1, il2
        cropfrac(i)=0.0
        vegfrac(i)=0.0
        temp(i)=0.0
        t1(i)=0
        mincfrac(i)=0.0
        barefrac(i)=1.0
        pbarefra(i)=1.0
        pvgbioms(i)=vgbiomas(i)  ! store grid average quantities in
        pgavltms(i)=gavgltms(i)  ! temporary arrays
        pgavscms(i)=gavgscms(i)
        vgbiomas(i)=0.0
        gavgltms(i)=0.0
        gavgscms(i)=0.0
        pltrmass(i,iccp1)=litrmass(i,iccp1)
        psocmass(i,iccp1)=soilcmas(i,iccp1)
        deadmass(i,iccp1)=0.0
        pdeadmas(i,iccp1)=0.0
        add2dead(i)=0.0
        gavgputa(i)=0.0 ! grid averaged value of c put aside for allocation
        gavgnpp(i)=0.0  ! grid averaged npp kg c/m2 for balance purposes
        bareiord(i)=0
        grsumlit(i)=0.0
        grsumsoc(i)=0.0
        ownlitr(i,iccp1)=0.0
        ownsolc(i,iccp1)=0.0
        incrlitr(i,iccp1)=0.0
        incrsolc(i,iccp1)=0.0
170   continue

      do 180 j = 1, icc-numcrops
        do 181 i = il1, il2
          rank(i,j)=j
          frac(i,j)=0.0
          c1(i,j)=0.0
          m1(i,j)=0.0
          usefrac(i,j)=0.0
          usec(i,j)=0.0
          usem(i,j)=0.0
          colterm(i,j)=0.0
          deathterm(i,j)=0.0
          term2(i,j)=0.0
          term3(i,j)=0.0
          term4(i,j)=0.0
          delfrac(i,j)=0.0
          exist1(i,j)=0
          useexist(i,j)=0
181     continue
180   continue

!     initial rank/superiority order for simulating competition. since
!     crops are not in competition their rank doesn't matter and
!     therefore we only have icc-2 ranks corresponding to the remaining 
!     pfts. the first icc-4 are tree pfts and the last two are the c3 and c4
!     grasses.
      do j = 1,icc-numcrops
         inirank(j)=j
      end do 

!    Estimate colonization and mortality rate for each pft, except for
!    crops whose fractional coverage is prescribed.

      do 200 j = 1, icc
       if(.not. crop(j))then  ! do not run for crops
        do 210 i = il1, il2

!         colonization rate (1/day). the factor (deltat/963.62) converts
!         npp from u-mol co2-c/m2.sec -> kg c/m2.day

          usebmsvg(i,j)= min(5.0, max(0.25, bmasveg(i,j)))

          colrate(i,j)=lambda(i,j)*max(0.0,nppveg(i,j))*(deltat/963.62)* &
                colmult*(1.0/(bio2sap(sort(j))*usebmsvg(i,j)))

!         mortality rate is the sum of growth related mortality,
!         intrinsic mortality, and an additional mortality that kicks in
!         when long term averaged bioclimatic conditions become
!         unfavourable for a pft. this last term is based on the
!         binary array pftexist.

          if(.not. pftexist(i,j))then
            mrtboclm(i,j)=bioclimrt/365.0
          else if(pftexist(i,j))then
            mrtboclm(i,j)=0.0
          endif

          mortrate(i,j)=geremort(i,j)+intrmort(i,j)+mrtboclm(i,j)

210     continue
       endif
200   continue

!    ---> from here on we assume that we only have icc-numcrops pfts <----
!              since crops are not part of the competition.

!     based on npp for each pft find the competition ranks / superiority 
!     order for simulating competition. note that crops
!     are not in competition, so the competition is between the
!     remaining pfts. in addition pfts which shouldn't exist in the
!     grid cell because of unfavourable values of long-term climatic
!     conditions are considered inferior.

!     find crop fraction

      cropfrac=0.0
      do 220 j = 1, icc
        if (crop(j)) then
         do 221 i = il1, il2
          cropfrac(i)=cropfrac(i)+fcancmx(i,j)
221      continue
        endif
220   continue

!     rank the tree pfts according to their colonization rates 
       
      do 250 j = 1, icc
        do 251 i = il1, il2
          usenppvg(i,j)=0.0 !assign to zero then check
          if (pftexist(i,j)) then
          usenppvg(i,j)=colrate(i,j)
          end if
251     continue
250   continue

      do 260 j = 1, icc-numcrops
        do 261 i = il1, il2
          rank(i,j)=inirank(j)
261     continue
260   continue

!     bubble sort according to colonization rates

      do 270 j = 1, numtreepfts
        do 280 n = 1, numtreepfts
          do 290 i = il1, il2
            if(usenppvg(i,n).lt.usenppvg(i,j))then
              temp(i)=usenppvg(i,n)
              usenppvg(i,n)=usenppvg(i,j)
              usenppvg(i,j)=temp(i)
              t1(i)=rank(i,n)
              rank(i,n)=rank(i,j)
              rank(i,j)=t1(i)
            endif
290       continue
280     continue
270   continue

!     the rank of c3 and c4 grass is also determined on the basis of
!     their npp but grasses are always assumed to be inferior to tree
!     pfts

      do 310 i = il1, il2 
        if(usenppvg(i,grass_ind(1)).ge.usenppvg(i,grass_ind(2)))then 
          rank(i,grass_ind(1)-numcrops)=grass_ind(1)-numcrops
          rank(i,grass_ind(2)-numcrops)=grass_ind(2)-numcrops
        elseif(usenppvg(i,grass_ind(1)).lt.usenppvg(i,grass_ind(2)))then
          rank(i,grass_ind(1)-numcrops)=grass_ind(2)-numcrops
          rank(i,grass_ind(2)-numcrops)=grass_ind(1)-numcrops
        endif
310   continue

!     with the ranks of all pfts in all grid cells we can now simulate
!     competition between them. for lotka-volterra eqns we need a
!     minimum seeding fraction otherwise the pfts will not expand at
!     all.

      do 330 j = 1, icc-numcrops   ! j now goes from 1 to icc-numcrops
        if(j.le.numtreepfts)then
          n=j
        else
          n=j+numcrops
        endif 
        do 340 i = il1, il2
          frac(i,j)=max(seed,fcancmx(i,n)) 
          if (pftexist(i,n)) then
           exist1(i,j)=1
           c1(i,j)=colrate(i,n)
          else
           exist1(i,j)=0
           c1(i,j)=0.0
          end if
          m1(i,j)=mortrate(i,n)
340     continue
330   continue

!     arrange colonization and mortality rates, and fractions, according
!     to superiority ranks

      do 350 n = 1, icc-numcrops   ! n now goes from 1 to icc-numcrops
        do 360 i = il1, il2
          usefrac(i,n)=frac(i,rank(i,n))
          usec(i,n)=c1(i,rank(i,n))
          usem(i,n)=m1(i,rank(i,n))
          useexist(i,n)=exist1(i,rank(i,n))
360     continue
350   continue

      do 400 n = 1, icc-numcrops   ! n now goes from 1 to icc-numcrops
        do 410 i = il1, il2

          colterm(i,n)=usec(i,n)*(usefrac(i,n)**a) ! colonization term

          sum1 = cropfrac(i)+seed !minbare
          do 420 k = 1, n-1, 1
            sum1 = sum1 + usefrac(i,k)
420       continue          

          term2(i,n)=usec(i,n)*(usefrac(i,n)**a)*(sum1+(usefrac(i,n)**b)) ! self & expansion thinning
          term3(i,n)=usem(i,n)*usefrac(i,n) ! mortality term

          sum2 = 0.0
          do 430 j = 1, n-1, 1
            sum3 = cropfrac(i)
            do 440 k = 1, j-1, 1
              sum3 = sum3 + usefrac(i,k)
440         continue
            sum4 = cropfrac(i)
            do 450 k = 1, j, 1
              sum4 = sum4 + usefrac(i,k)
450         continue
            sum2 = sum2 + ( &
            ( ((1.-sum3)**g)*usec(i,j)*(usefrac(i,j)**a)*usefrac(i,n) )/ &
            ( (1.-sum4)**g )  )
430       continue
          term4(i,n)=sum2  ! invasion
          deathterm(i,n) = term2(i,n) + term3(i,n) + term4(i,n)
          delfrac(i,n)=colterm(i,n)-deathterm(i,n) ! delta fraction

410     continue
400   continue

!     update fractions and check if all fractions are +ve 

      do 500 n = 1, icc-numcrops
        do 510 i = il1, il2
          usefrac(i,n)=usefrac(i,n)+delfrac(i,n)
          if(usefrac(i,n).lt.0.0)then
            write(6,*)'fractional coverage -ve for cell ',i,' and pft',n
            call xit('competition',-5)
          endif
          usefrac(i,n)=max(seed,usefrac(i,n))
510     continue
500   continue

!     with the minimum seeding fraction prescription, especially for
!     lotka volterra eqns the total veg fraction may exceed 1. to
!     prevent this we need to adjust fractional coverage of all non-crop
!     pfts that do not have the minimum fraction.

      do 530 i = il1, il2
        vegfrac(i)=cropfrac(i)+ seed  !total vegetation fraction
        mincfrac(i)=cropfrac(i)+ seed !sum of mininum prescribed & crop fractions
530   continue

      do 540 n = 1, icc-numcrops
        do 541 i = il1, il2
          vegfrac(i)=vegfrac(i)+usefrac(i,n)
          if(abs(usefrac(i,n)-seed).le.zero) then
            mincfrac(i)= mincfrac(i)+ usefrac(i,n)
          endif 
541     continue
540   continue

      do 550 n = 1, icc-numcrops
        do 551 i = il1, il2
          if(vegfrac(i).gt.1.0.and. &
          abs(usefrac(i,n)-seed).gt.zero) then
            term =(1.-mincfrac(i))/(vegfrac(i)-mincfrac(i)) 
            usefrac(i,n)=usefrac(i,n)*term
          endif
551     continue
550   continue

!     check again that total veg frac doesn't exceed 1.

      do 560 i = il1, il2
        vegfrac(i)=cropfrac(i)+ seed  !total vegetation fraction
560   continue

      do 570 n = 1, icc-numcrops
        do 571 i = il1, il2
          vegfrac(i)=vegfrac(i)+usefrac(i,n)
571     continue
570   continue

      do 580 i = il1, il2
        if(vegfrac(i).gt.1.0+1e-5)then
          write(6,*)'vegetation fraction in cell ',i,' greater than'
          write(6,*)'1.0 and equal to ',vegfrac(i) 
          call xit('competition',-6)
        endif
580   continue

!     map delfrac to chngfrac so that we get change in fraction
!     corresponding to the actual number of pfts

      do 590 j = 1, icc-numcrops   ! j now goes from 1 to icc-numcrops
        do 591 i = il1, il2

          if(rank(i,j).le.numtreepfts)then
            k=rank(i,j)
          else
            k=rank(i,j)+2
          endif 
          expnterm(i,k)=colterm(i,j)
          mortterm(i,k)=deathterm(i,j)
          fcancmx(i,k)=usefrac(i,j)
          chngfrac(i,k)=fcancmx(i,k)-pfcancmx(i,k)

591     continue
590   continue

!       ---> from here on we get back to our usual icc pfts <----

!     get bare fraction

      do 600 j = 1, icc
        do 601 i = il1, il2
          barefrac(i)=barefrac(i)-fcancmx(i,j)
          pbarefra(i)=pbarefra(i)-pfcancmx(i,j)
601     continue
600   continue

!     check if a pft's fractional cover is increasing or decreasing

      do 620 j = 1, icc
        do 621 i = il1, il2
          if( ( fcancmx(i,j).gt.pfcancmx(i,j)) .and. &
             (abs(pfcancmx(i,j)-fcancmx(i,j)).gt.zero) ) then
              fraciord(i,j)=1
          else if( ( fcancmx(i,j).lt.pfcancmx(i,j)) .and. &
                  (abs(pfcancmx(i,j)-fcancmx(i,j)).gt.zero) ) then
              fraciord(i,j)=-1
          endif
621     continue
620   continue

!     check if bare fraction increases or decreases

      do 640 i = il1, il2
        if( ( barefrac(i).gt.pbarefra(i)) .and. &
           (abs(pbarefra(i)-barefrac(i)).gt.zero) ) then
              bareiord(i)=1
        else if ( ( barefrac(i).lt.pbarefra(i)) .and. &
                 (abs(pbarefra(i)-barefrac(i)).gt.zero) ) then
              bareiord(i)=-1
        endif
640   continue

!     now that we know the change in fraction for every pft we use its
!     npp for spatial expansion and litter generation. we also spread
!     vegetation biomass uniformly over the new fractions, and generate
!     additional litter from mortality if the fractions decrease.
!
!     three things can happen here
!
!     1. fraciord = 0, which means all npp that was used for expansion 
!        becomes litter, due to self/expansion thinning and mortality.
!
!     2. fraciord = 1, which means a part of or full npp is used for
!        expansion but some litter may also be generated. the part of 
!        npp that is used for expansion needs to be allocated to leaves,
!        stem, and root. rather than doing this here we will let the
!        allocation part handle this. so allocation module will allocate
!        not only the npp that is used for pure vertical expansion but 
!        also this npp. but we will do our part here and spread the
!        vegetation biomass over the new increased fraction.
!
!     3. fraciord = -1, which means all of the npp is to be used for
!        litter generation but in addition some more litter will be
!        generated from mortality of the standing biomass.

      do 660 j = 1, icc
       if(.not. crop(j))then  ! do not run for crops
        do 661 i = il1, il2

          if(fraciord(i,j).eq.1)then

!           reduce biomass density by spreading over larger fraction

            term = (pfcancmx(i,j)/fcancmx(i,j))
            gleafmas(i,j) = gleafmas(i,j)*term
            bleafmas(i,j) = bleafmas(i,j)*term
            stemmass(i,j) = stemmass(i,j)*term
            rootmass(i,j) = rootmass(i,j)*term
            litrmass(i,j) = litrmass(i,j)*term
            soilcmas(i,j) = soilcmas(i,j)*term

!           only a fraction of npp becomes litter which for simplicity
!           and for now we spread over the whole grid cell

            if(expnterm(i,j).le.zero.and.mortterm(i,j).gt.zero)then
              write(6,*)'expansion term<= zero when fractional coverage'
              write(6,*)'is increasing for pft',j,' in grid cell',i
              write(*,*)'pfcancmx(',i,',',j,')=',pfcancmx(i,j)
              write(*,*)'fcancmx(',i,',',j,')=',fcancmx(i,j)
              write(*,*)'expnterm(',i,',',j,')=',expnterm(i,j)
              call xit('competition',-7)
            else if(expnterm(i,j).le.zero.and.mortterm(i,j).le.zero)then
              term = 1.0
            else
              term = (mortterm(i,j)/expnterm(i,j))
            endif

            add2allo(i,j)=(1.-term)

!           the factor (deltat/963.62) converts npp from u-mol co2-c/m2.sec 
!           -> kg c/m2.deltat

            term = term*max(0.0,nppveg(i,j))*(deltat/963.62)*lambda(i,j)*pfcancmx(i,j)

            incrlitr(i,j)=term
            grsumlit(i)=grsumlit(i)+incrlitr(i,j)

!           rest put aside for allocation

            add2allo(i,j) = add2allo(i,j)* max(0.0,nppveg(i,j))*(deltat/963.62)*lambda(i,j)*&
                          (pfcancmx(i,j)/fcancmx(i,j))


          else if(fraciord(i,j).eq.-1)then

!           all npp used for expansion becomes litter plus there is
!           additional mortality of the standing biomass. the npp that 
!           becomes litter is now spread over the whole grid cell.
!           all biomass from fraction that dies due to mortality is 
!           also distributed over the litter pool of whole grid cell.

            term = abs(chngfrac(i,j))*grclarea(i)*1.0e06*(gleafmas(i,j)+ &
               bleafmas(i,j)+stemmass(i,j)+rootmass(i,j)+litrmass(i,j))
            term = term/(grclarea(i)*1.0e06)

            incrlitr(i,j)=term

            term = max(0.0,nppveg(i,j))*(deltat/963.62)*lambda(i,j)*pfcancmx(i,j)

            incrlitr(i,j) = incrlitr(i,j)+term
            grsumlit(i)=grsumlit(i)+incrlitr(i,j)

!           chop off soil c from the fraction that goes down and
!           spread it uniformly over the soil c pool of entire grid cell

            term = abs(chngfrac(i,j))*grclarea(i)*1.0e06*soilcmas(i,j)
            term = term/(grclarea(i)*1.0e06)
            incrsolc(i,j)=term
            grsumsoc(i)=grsumsoc(i)+incrsolc(i,j)

          else if(fraciord(i,j).eq.0)then

!           all npp used for expansion becomes litter

            incrlitr(i,j) =max(0.0,nppveg(i,j))*(deltat/963.62)*lambda(i,j)* pfcancmx(i,j)
            grsumlit(i)=grsumlit(i)+incrlitr(i,j)

          endif

661     continue
       endif
660   continue

!     if bare fraction decreases then chop off the litter and soil c
!     from the decreased fraction and add it to grsumlit & grsumsoc
!     for spreading over the whole grid cell. if bare fraction increases
!     then spread its litter and soil c uniformly over the increased 
!     fraction.

      do 680 i = il1, il2
        if(bareiord(i).eq.-1)then

          term =(pbarefra(i)-barefrac(i))*litrmass(i,iccp1)
          incrlitr(i,iccp1)=term
          grsumlit(i)=grsumlit(i)+term

          term =(pbarefra(i)-barefrac(i))*soilcmas(i,iccp1)
          incrsolc(i,iccp1)=term
          grsumsoc(i)=grsumsoc(i)+term

        else if(bareiord(i).eq.1)then

          term = pbarefra(i)/barefrac(i)
          litrmass(i,iccp1)=litrmass(i,iccp1)*term
          soilcmas(i,iccp1)=soilcmas(i,iccp1)*term

        endif
680   continue

!     if a pft is not suppose to exist as indicated by pftexist and its 
!     fractional coverage is really small then get rid of the pft all
!     together and spread its live and dead biomass over the grid cell.

      do 690 j = 1, icc
        do 691 i = il1, il2
          if(.not. pftexist(i,j).and.fcancmx(i,j).lt.1.0e-05)then

            term = fcancmx(i,j)*(gleafmas(i,j)+bleafmas(i,j) &
              +stemmass(i,j)+rootmass(i,j)+litrmass(i,j))
            incrlitr(i,j)=incrlitr(i,j) + term
            grsumlit(i)=grsumlit(i)+incrlitr(i,j)

            term = fcancmx(i,j)*soilcmas(i,j)
            incrsolc(i,j)=incrsolc(i,j) + term
            grsumsoc(i)=grsumsoc(i)+incrsolc(i,j)

            barefrac(i)=barefrac(i)+fcancmx(i,j)

!           adjust litter and soil c mass densities for increase in
!           barefrac over the bare fraction.

            term = (barefrac(i)-fcancmx(i,j))/barefrac(i)
            litrmass(i,iccp1) = litrmass(i,iccp1)*term
            soilcmas(i,iccp1) = soilcmas(i,iccp1)*term

            fcancmx(i,j)=0.0
          endif
691     continue
690   continue

!     spread litter and soil c over all pfts and the barefrac

      do 700 j = 1, icc
        do 701 i = il1, il2
          if(fcancmx(i,j).gt.zero)then
            litrmass(i,j)=litrmass(i,j)+grsumlit(i)
            soilcmas(i,j)=soilcmas(i,j)+grsumsoc(i)
          else
            gleafmas(i,j)=0.0
            bleafmas(i,j)=0.0
            stemmass(i,j)=0.0
            rootmass(i,j)=0.0
            litrmass(i,j)=0.0
            soilcmas(i,j)=0.0
          endif
701     continue
700   continue

      do 720 i = il1, il2
        if(barefrac(i).gt.zero)then
          litrmass(i,iccp1)=litrmass(i,iccp1)+grsumlit(i)
          soilcmas(i,iccp1)=soilcmas(i,iccp1)+grsumsoc(i)
        else
          litrmass(i,iccp1)=0.0
          soilcmas(i,iccp1)=0.0
        endif
720   continue

!     get fcanmxs for use by class based on the new fcancmxs

      do 740 j = 1, ican
        do 741 i = il1, il2
           fcanmx(i,j)=0.0 ! fractional coverage of class' pfts
741     continue
740   continue

      k1=0
      do 750 j = 1, ican
        if(j.eq.1) then
          k1 = k1 + 1
        else
          k1 = k1 + nol2pfts(j-1)
        endif
        k2 = k1 + nol2pfts(j) - 1
        do 751 l = k1, k2
          do 752 i = il1, il2
            fcanmx(i,j)=fcanmx(i,j)+fcancmx(i,l)
752       continue
751     continue
750   continue

!     update grid averaged vegetation biomass, and litter and soil c densities

      do 800 j = 1, icc
        do 801 i = il1, il2
          vgbiomas(i)=vgbiomas(i)+fcancmx(i,j)*(gleafmas(i,j)+ &
                     bleafmas(i,j)+stemmass(i,j)+rootmass(i,j))
          gavgltms(i)=gavgltms(i)+fcancmx(i,j)*litrmass(i,j)
          gavgscms(i)=gavgscms(i)+fcancmx(i,j)*soilcmas(i,j)
801     continue
800   continue

      do 810 i = il1, il2
        gavgltms(i)=gavgltms(i)+( barefrac(i)*litrmass(i,iccp1) )
        gavgscms(i)=gavgscms(i)+( barefrac(i)*soilcmas(i,iccp1) )
810   continue

!     and finally we check the c balance. we were suppose to use a
!     fraction of npp for competition. some of it is used for expansion
!     (this is what we save for allocation), and the rest becomes litter. 
!     so for each pft the total c mass in vegetation and litter pools
!     must all add up to the same value as before competition.

!     JM: I changed this to deal with density, not kg/gridcell.

      do 830 j = 1, icc
       if (.not. crop(j)) then  
        do 831 i = il1, il2

          biomasvg(i,j)=fcancmx(i,j)* & !grclarea(i)*1.0e06* &
           (gleafmas(i,j)+bleafmas(i,j)+stemmass(i,j)+rootmass(i,j)) 
          pbiomasvg(i,j)=pfcancmx(i,j)* & !grclarea(i)*1.0e06* &
           (pglfmass(i,j)+pblfmass(i,j)+protmass(i,j)+pstmmass(i,j)) 

!         part of npp that we will use later for allocation
          putaside(i,j)=add2allo(i,j)*fcancmx(i,j) !*grclarea(i)*1.0e06

          gavgputa(i) = gavgputa(i) + putaside(i,j)

!         litter added to bare
          barelitr(i,j)=grsumlit(i)*fcancmx(i,j) !*grclarea(i)*1.0e06
          ownlitr(i,j)=incrlitr(i,j) !*grclarea(i)*1.0e06

!         soil c added to bare
          baresolc(i,j)=grsumsoc(i)*fcancmx(i,j) !*grclarea(i)*1.0e06
          ownsolc(i,j)=incrsolc(i,j) ! *grclarea(i)*1.0e06

          add2dead(i) = add2dead(i) + barelitr(i,j) + baresolc(i,j)

!         npp we had in first place to expand
          nppvegar(i,j)=max(0.0,nppveg(i,j))*(deltat/963.62)*lambda(i,j)*pfcancmx(i,j) !*grclarea(i)*1.0e06

          gavgnpp(i) = gavgnpp(i) + nppvegar(i,j)

          deadmass(i,j)=fcancmx(i,j)*(litrmass(i,j)+soilcmas(i,j)) !*grclarea(i)*1.0e06
          pdeadmas(i,j)=pfcancmx(i,j)*(pltrmass(i,j)+psocmass(i,j)) !*grclarea(i)*1.0e06

!         total mass before competition
          befrmass=pbiomasvg(i,j)+nppvegar(i,j)+pdeadmas(i,j)

!         total mass after competition
          aftrmass=biomasvg(i,j)+putaside(i,j)+deadmass(i,j)- &
                  barelitr(i,j)-baresolc(i,j)+ownlitr(i,j)+ownsolc(i,j)

          if(abs(befrmass-aftrmass).gt.tolrance)then
            write(6,*)'total biomass for pft',j,', and grid cell =',i 
            write(6,*)'does not balance before and after competition'
            write(6,*)' '
            write(6,*)'chngfrac(',i,',',j,')=',chngfrac(i,j)
            write(6,*)'fraciord(',i,',',j,')=',fraciord(i,j)
            write(6,*)' '
            write(6,*)'pbiomasvg(',i,',',j,')=',pbiomasvg(i,j)
            write(6,*)'pdeadmas(',i,',',j,')=',pdeadmas(i,j)
            write(6,*)'nppvegar(',i,',',j,')=',nppvegar(i,j)
            write(6,*)' '
            write(6,*)' biomasvg(',i,',',j,')=',biomasvg(i,j)
            write(6,*)'deadmass(',i,',',j,')=',deadmass(i,j)
            write(6,*)'putaside(',i,',',j,')=',putaside(i,j)
            write(6,*)' '
            write(6,*)'before biomass density = ',gleafmas(i,j)+ &
                  bleafmas(i,j)+stemmass(i,j)+rootmass(i,j)
            write(6,*)'after  biomass density = ',pglfmass(i,j)+ &
                 pblfmass(i,j)+protmass(i,j)+pstmmass(i,j)
            write(6,*)' '
            write(6,*)'barelitr(',i,',',j,')=',barelitr(i,j)
            write(6,*)'baresolc(',i,',',j,')=',baresolc(i,j)
            write(6,*)'ownlitr(',i,',',j,')=',ownlitr(i,j)
            write(6,*)'ownsolc(',i,',',j,')=',ownsolc(i,j)
            write(6,*)' '
           ! write(6,*)'grclarea(',i,')=',grclarea(i)
            write(6,*)'fcancmx(',i,',',j,')=',fcancmx(i,j)
            write(6,*)'pfcancmx(',i,',',j,')=',pfcancmx(i,j)
            write(6,*)' '
            write(6,*)'abs(befrmass-aftrmass)=',abs(befrmass-aftrmass)
            call xit('competition',-8)
          endif

831     continue
       endif
830   continue

!     check balance over the bare fraction

      j = iccp1
        do 851 i = il1, il2
          deadmass(i,j)=barefrac(i)*(litrmass(i,j)+soilcmas(i,j)) !*grclarea(i)*1.0e06
          pdeadmas(i,j)=pbarefra(i)*(pltrmass(i,j)+psocmass(i,j)) !*grclarea(i)*1.0e06

          add2dead(i)=(grsumlit(i)+grsumsoc(i))*barefrac(i) !*grclarea(i)*1.0e06

          ownlitr(i,j)=incrlitr(i,j) !*grclarea(i)*1.0e06
          ownsolc(i,j)=incrsolc(i,j) !*grclarea(i)*1.0e06

          befrmass=pdeadmas(i,j)+add2dead(i)
          aftrmass=deadmass(i,j)+ownlitr(i,j)+ownsolc(i,j)

          if(abs(befrmass-aftrmass).gt.tolrance)then
            write(6,*)'total dead mass for grid cell =',i,'does not balance over bare' 
            write(6,*)'pdeadmas(',i,',',j,')=',pdeadmas(i,j)
            write(6,*)'add2dead(',i,') term=',add2dead(i)
            write(6,*)'deadmass(',i,',',j,')=',deadmass(i,j)
            write(6,*)' '
            write(6,*)'ownlitr(',i,',',j,')=',ownlitr(i,j)
            write(6,*)'ownsolc(',i,',',j,')=',ownsolc(i,j)
            write(6,*)' '
            write(6,*)'pbarefra(',i,')=',pbarefra(i)
            write(6,*)'bareiord(',i,') =',bareiord(i)
            write(6,*)'barefrac(',i,')=',barefrac(i)
            write(6,*)' '
            write(6,*)'abs(befrmass-aftrmass)=',abs(befrmass-aftrmass)
            call xit('competition',-9)
          endif
851     continue


!     grid averaged densities must also balance

!     JM: Now done above.

!      do 870 i = il1, il2
!        befrmass=(pvgbioms(i)+pgavltms(i)+pgavscms(i))+gavgnpp(i) /(grclarea(i)*1.0e06)
!        aftrmass=(vgbiomas(i)+gavgltms(i)+gavgscms(i))+gavgputa(i)/(grclarea(i)*1.0e06)
!        if(abs(befrmass-aftrmass).gt.tolrance)then
!          write(6,*)'total (live+dead) mass for grid cell =',i,'does not balance' 
!          write(6,*)'abs(befrmass-aftrmass)',abs(befrmass-aftrmass),'is gt our tolerance of',tolrance
!          write(6,*)'pvgbioms(',i,')=',pvgbioms(i)*grclarea(i)*1.0e06
!          write(6,*)'pgavltms(',i,')=',pgavltms(i)*grclarea(i)*1.0e06
!          write(6,*)'pgavscms(',i,')=',pgavscms(i)*grclarea(i)*1.0e06
!          write(6,*)'gavgnpp(',i,')=',gavgnpp(i)
!          write(6,*)'befrmass*1.0e06 =',befrmass*1.0e06
!          write(6,*)' '
!          write(6,*)'vgbiomas(',i,')=',vgbiomas(i)*grclarea(i)*1.0e06
!          write(6,*)'gavgltms(',i,')=',gavgltms(i)*grclarea(i)*1.0e06
!          write(6,*)'gavgscms(',i,')=',gavgscms(i)*grclarea(i)*1.0e06
!          write(6,*)'gavgputa(',i,')=',gavgputa(i)
!          write(6,*)'aftrmass*1.0e06=',aftrmass*1.0e06
!          call xit('competition',-10)
!        endif
!870   continue

      return

end subroutine competition

end module


