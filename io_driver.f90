module io_driver

! Central module that handles all CTEM reading and writing of external files

! J. Melton Mar 30 2015

implicit none

! subroutines contained in this module:
public  :: read_from_ctm
public  :: write_ctm_rs
public  :: create_outfiles
public  :: class_monthly_aw
public  :: ctem_daily_aw
public  :: ctem_monthly_aw
public  :: ctem_annual_aw
public  :: close_outfiles

contains

!-------------------------------------------------------------------------------------------------------------

subroutine read_from_ctm(nltest,nmtest,FCANROT,FAREROT,RSMNROT,QA50ROT, &
                         VPDAROT,VPDBROT,PSGAROT,PSGBROT,DRNROT,SDEPROT, &
                         XSLPROT,GRKFROT,WFSFROT,WFCIROT,MIDROT,SANDROT, &
                         CLAYROT,ORGMROT,TBARROT,THLQROT,THICROT,TCANROT, &
                         TSNOROT,TPNDROT,ZPNDROT,RCANROT,SCANROT,SNOROT, &
                         ALBSROT,RHOSROT,GROROT,argbuff)

!   This subroutine reads in the restart/starting conditions from
!   the .CTM file. The input values are checked and possibly adjusted
!   if the run is intended to be with competition on and if the start_bare
!   flag is true. 

use ctem_params,        only : icc,iccp1,nmos,seed,ignd,ilg,icp1,nlat,ican,abszero
use ctem_statevars,     only : c_switch,vrot,vgat

implicit none


! arguments:
character(80), intent(in) :: argbuff
integer, intent(in) :: nltest
integer, intent(inout) :: nmtest
real, dimension(nlat,nmos,icp1), intent(inout) :: FCANROT
real, dimension(nlat,nmos), intent(inout) :: FAREROT
real, dimension(nlat,nmos,ican), intent(inout) :: RSMNROT
real, dimension(nlat,nmos,ican), intent(inout) :: QA50ROT
real, dimension(nlat,nmos,ican), intent(inout) :: VPDAROT
real, dimension(nlat,nmos,ican), intent(inout) :: VPDBROT
real, dimension(nlat,nmos,ican), intent(inout) :: PSGAROT
real, dimension(nlat,nmos,ican), intent(inout) :: PSGBROT
real, dimension(nlat,nmos), intent(inout) :: DRNROT
real, dimension(nlat,nmos), intent(inout) :: SDEPROT
real, dimension(nlat,nmos), intent(inout) :: XSLPROT
real, dimension(nlat,nmos), intent(inout) :: GRKFROT
real, dimension(nlat,nmos), intent(inout) :: WFSFROT
real, dimension(nlat,nmos), intent(inout) :: WFCIROT
integer, dimension(nlat,nmos), intent(inout) :: MIDROT
real, dimension(nlat,nmos,ignd), intent(inout) :: SANDROT
real, dimension(nlat,nmos,ignd), intent(inout) :: CLAYROT
real, dimension(nlat,nmos,ignd), intent(inout) :: ORGMROT
real, dimension(nlat,nmos,ignd), intent(inout) :: TBARROT
real, dimension(nlat,nmos,ignd), intent(inout) :: THLQROT
real, dimension(nlat,nmos,ignd), intent(inout) :: THICROT
real, dimension(nlat,nmos), intent(inout) :: TCANROT
real, dimension(nlat,nmos), intent(inout) :: TSNOROT
real, dimension(nlat,nmos), intent(inout) :: TPNDROT
real, dimension(nlat,nmos), intent(inout) :: ZPNDROT
real, dimension(nlat,nmos), intent(inout) :: RCANROT
real, dimension(nlat,nmos), intent(inout) :: SCANROT
real, dimension(nlat,nmos), intent(inout) :: SNOROT
real, dimension(nlat,nmos), intent(inout) :: ALBSROT
real, dimension(nlat,nmos), intent(inout) :: RHOSROT
real, dimension(nlat,nmos), intent(inout) :: GROROT

! pointers:
logical, pointer :: mosaic
logical, pointer :: dofire
logical, pointer :: compete
logical, pointer :: inibioclim
logical, pointer :: dowetlands
logical, pointer :: start_bare
logical, pointer :: lnduseon
logical, pointer :: obswetf
character(80), pointer :: titlec1
character(80), pointer :: titlec2
character(80), pointer :: titlec3
real, pointer, dimension(:,:,:) :: ailcminrow           !
real, pointer, dimension(:,:,:) :: ailcmaxrow           !
real, pointer, dimension(:,:,:) :: dvdfcanrow           !
real, pointer, dimension(:,:,:) :: gleafmasrow          !
real, pointer, dimension(:,:,:) :: bleafmasrow          !
real, pointer, dimension(:,:,:) :: stemmassrow          !
real, pointer, dimension(:,:,:) :: rootmassrow          !
real, pointer, dimension(:,:,:) :: pstemmassrow         !
real, pointer, dimension(:,:,:) :: pgleafmassrow        !     
real, pointer, dimension(:) :: twarmm                   ! temperature of the warmest month (c)
real, pointer, dimension(:) :: tcoldm                   ! temperature of the coldest month (c)
real, pointer, dimension(:) :: gdd5                     ! growing degree days above 5 c
real, pointer, dimension(:) :: aridity                  ! aridity index, ratio of potential evaporation to precipitation
real, pointer, dimension(:) :: srplsmon                 ! number of months in a year with surplus water i.e.precipitation more than potential evaporation
real, pointer, dimension(:) :: defctmon                 ! number of months in a year with water deficit i.e.precipitation less than potential evaporation
real, pointer, dimension(:) :: anndefct                 ! annual water deficit (mm) 
real, pointer, dimension(:) :: annsrpls                 ! annual water surplus (mm)
real, pointer, dimension(:) :: annpcp                   ! annual precipitation (mm)
real, pointer, dimension(:) :: dry_season_length        ! length of dry season (months)
real, pointer, dimension(:,:,:) :: litrmassrow
real, pointer, dimension(:,:,:) :: soilcmasrow
real, pointer, dimension(:) :: extnprobgrd
real, pointer, dimension(:) :: prbfrhucgrd
real, pointer, dimension(:,:) :: mlightnggrd  
integer, pointer, dimension(:,:,:) :: lfstatusrow
integer, pointer, dimension(:,:,:) :: pandaysrow
integer, pointer, dimension(:) :: stdalngrd
real, pointer, dimension(:,:) :: wetfrac_sgrd

! local variables

integer :: i,m,j,strlen
real, dimension(ilg,2) :: crop_temp_frac

! point pointers:
mosaic            => c_switch%mosaic
dofire            => c_switch%dofire
compete           => c_switch%compete
inibioclim        => c_switch%inibioclim
dowetlands        => c_switch%dowetlands
start_bare        => c_switch%start_bare    
lnduseon          => c_switch%lnduseon
obswetf           => c_switch%obswetf
titlec1           => c_switch%titlec1
titlec2           => c_switch%titlec2
titlec3           => c_switch%titlec3
ailcminrow        => vrot%ailcmin
ailcmaxrow        => vrot%ailcmax
dvdfcanrow        => vrot%dvdfcan
gleafmasrow       => vrot%gleafmas
bleafmasrow       => vrot%bleafmas
stemmassrow       => vrot%stemmass
rootmassrow       => vrot%rootmass
pstemmassrow      => vrot%pstemmass
pgleafmassrow     => vrot%pgleafmass
twarmm            => vgat%twarmm
tcoldm            => vgat%tcoldm
gdd5              => vgat%gdd5
aridity           => vgat%aridity
srplsmon          => vgat%srplsmon
defctmon          => vgat%defctmon
anndefct          => vgat%anndefct
annsrpls          => vgat%annsrpls
annpcp            => vgat%annpcp
dry_season_length => vgat%dry_season_length
litrmassrow       => vrot%litrmass
soilcmasrow       => vrot%soilcmas
extnprobgrd       => vrot%extnprob
prbfrhucgrd       => vrot%prbfrhuc
mlightnggrd       => vrot%mlightng 
wetfrac_sgrd      => vgat%wetfrac_s
stdalngrd         => vrot%stdaln
lfstatusrow       => vrot%lfstatus
pandaysrow        => vrot%pandays
      
! -----------------      
! Begin      

open(unit=11,file=argbuff(1:strlen(argbuff))//'.CTM', status='old')

read (11,7010) titlec1
read (11,7010) titlec2
read (11,7010) titlec3

7010  FORMAT(A80)

!  read from ctem initialization file (.CTM)

  do 71 i=1,nltest
    do 72 m=1,nmtest

!       The following three variables are needed to run ctem. 
!       min & max leaf area index are needed to break
!       class lai into dcd and evg for trees (for crops and grasses it
!       doesn't matter much).

!       dvdfcanrow is needed to divide needle & broad leaf into dcd and evg,
!       and crops & grasses into c3 and c4 fractions.

        read(11,*) (ailcminrow(i,m,j),j=1,icc)
        read(11,*) (ailcmaxrow(i,m,j),j=1,icc)
        read(11,*) (dvdfcanrow(i,m,j),j=1,icc)

!       rest of the initialization variables are needed to run ctem.
!       if starting from bare ground initialize all live and dead c pools from zero. suitable values
!       of extnprobgrd and prbfrhucgrd would still be required. set stdalngrd to
!       1 for operation in non-gcm stand alone mode, in the ctem
!       initialization file.

        read(11,*) (gleafmasrow(i,m,j),j=1,icc)
        read(11,*) (bleafmasrow(i,m,j),j=1,icc)
        read(11,*) (stemmassrow(i,m,j),j=1,icc)
        read(11,*) (rootmassrow(i,m,j),j=1,icc)

!       If fire and competition are on, save the stemmass and rootmass for
!       use in burntobare subroutine on the first timestep.
        if (dofire .and. compete) then
            do j =1,icc
            pstemmassrow(i,m,j)=stemmassrow(i,m,j)
            pgleafmassrow(i,m,j)=rootmassrow(i,m,j)    
            end do           
        end if

        read(11,*) (litrmassrow(i,m,j),j=1,iccp1)
        read(11,*) (soilcmasrow(i,m,j),j=1,iccp1)
        read(11,*) (lfstatusrow(i,m,j),j=1,icc)
        read(11,*) (pandaysrow(i,m,j),j=1,icc)

72      continue

        read(11,*) (mlightnggrd(i,j),j=1,6)  !mean monthly lightning frequency
        read(11,*) (mlightnggrd(i,j),j=7,12) !flashes/km2.year
        read(11,*) extnprobgrd(i)
        read(11,*) prbfrhucgrd(i)
        read(11,*) stdalngrd(i)

        if (compete .and. inibioclim) then  !read in the bioclimatic parameters
        read(11,*) twarmm(i), tcoldm(i), gdd5(i), aridity(i),srplsmon(i)
        read(11,*) defctmon(i), anndefct(i), annsrpls(i), annpcp(i), dry_season_length(i)
        else if (compete .and. .not. inibioclim) then ! set them to zero
        twarmm(i)=0.0
        tcoldm(i)=0.0
        gdd5(i)=0.0
        aridity(i)=0.0
        srplsmon(i)=0.0
        defctmon(i)=0.0
        anndefct(i)=0.0
        annsrpls(i)=0.0
        annpcp(i)=0.0
        dry_season_length(i) = 0.0
        endif

        if (dowetlands) then      ! Rudra !if true then read wetland fractions
            read(11,*) (wetfrac_sgrd(i,j),j=1,8)
        endif   
71    continue

close(11)
    

!     check that a competition or luc run has the correct number of mosaics
!     if it is not a start_bare run, then nmtest should equal nmos
      if (mosaic .and. (compete .or. lnduseon) .and. .not. start_bare) then
        if (nmtest .ne. nmos) then
           write(6,*)'compete or luc runs that do not start from bare'
           write(6,*)'ground need the number of mosaics to equal icc+1'
           write(6,*)'nmtest = ',nmtest,' nmos = ',nmos
            call xit('runclass36ctem', -2)
        endif
      endif

!     if this run uses the competition or lnduseon parameterization and
!     starts from bare ground, set up the model state here. this 
!     overwrites what was read in from the .ini and .ctm files. 
!     for composite runs (the composite set up is after this one for mosaics)
      if ((compete .or. lnduseon) .and. start_bare) then

       if (mosaic) then 

!       store the read-in crop fractions as we keep them even when we start bare. 
!       FLAG: this is setup assuming that crops are in mosaics 6 and 7. JM Apr 9 2014.
         do i=1,nltest
          crop_temp_frac(i,1)=FAREROT(i,6)
          crop_temp_frac(i,2)=FAREROT(i,7)
         end do

!       check the number of mosaics that came from the .ini file
        if (nmtest .ne. nmos) then

!        we need to transfer some initial parameterization info to all
!        mosaics, so set all values to that of the first mosaic.
         do i=1,nltest
          do m=nmtest+1,nmos

           do j=1,ican
             RSMNROT(i,m,j)=RSMNROT(i,1,j)
             QA50ROT(i,m,j)=QA50ROT(i,1,j)
             VPDAROT(i,m,j)=VPDAROT(i,1,j)
             VPDBROT(i,m,j)=VPDBROT(i,1,j)
             PSGAROT(i,m,j)=PSGAROT(i,1,j)
             PSGBROT(i,m,j)=PSGBROT(i,1,j)
           enddo

           DRNROT(i,m)=DRNROT(i,1)
           SDEPROT(i,m)=SDEPROT(i,1)
           FAREROT(i,m)=FAREROT(i,1)
           XSLPROT(i,m)=XSLPROT(i,1)
           GRKFROT(i,m)=GRKFROT(i,1)
           WFSFROT(i,m)=WFSFROT(i,1)
           WFCIROT(i,m)=WFCIROT(i,1)
           MIDROT(i,m)=MIDROT(i,1)

           do j=1,3
            SANDROT(i,m,j)=SANDROT(i,1,j)
            CLAYROT(i,m,j)=CLAYROT(i,1,j)
            ORGMROT(i,m,j)=ORGMROT(i,1,j)
            TBARROT(i,m,j)=TBARROT(i,1,j)
            THLQROT(i,m,j)=THLQROT(i,1,j)
            THICROT(i,m,j)=THICROT(i,1,j)
           enddo

           TCANROT(i,m)=TCANROT(i,1)
           TSNOROT(i,m)=TSNOROT(i,1)
           TPNDROT(i,m)=TPNDROT(i,1)
           ZPNDROT(i,m)=ZPNDROT(i,1)
           RCANROT(i,m)=RCANROT(i,1)
           SCANROT(i,m)=SCANROT(i,1)
           SNOROT(i,m)=SNOROT(i,1)
           ALBSROT(i,m)=ALBSROT(i,1)
           RHOSROT(i,m)=RHOSROT(i,1)
           GROROT(i,m)=GROROT(i,1)
           do j=1,icc
             lfstatusrow(i,m,j) = 4
           enddo !j

          enddo !m
         enddo !i

!       set the number of mosaics to icc+1        
        nmtest=nmos

        endif  !if (nmtest .ne. nmos)

!       set the initial conditions for the pfts
!       (bah, this is such an inelegant way to do this, but oh well...)

!       initalize to zero
        FCANROT=0.0
        dvdfcanrow=0.0
        FAREROT=0.0

        do i=1,nltest
         do m=1,nmtest

!         set the seed amount for each pft in its mosaic
          if (compete .or. lnduseon) then
            if (m .lt. icc+1) then
             FAREROT(i,m)=seed
            else
             FAREROT(i,m)=1.0 - (real(icc) * seed)
            endif
          endif

          do j = 1,icc
            ailcminrow(i,m,j)=0.0
            ailcmaxrow(i,m,j)=0.0
            gleafmasrow(i,m,j)=0.0
            bleafmasrow(i,m,j)=0.0
            stemmassrow(i,m,j)=0.0
            rootmassrow(i,m,j)=0.0
            lfstatusrow(i,m,j)=4
            pandaysrow(i,m,j)=0
          enddo
  
          lfstatusrow(i,m,1)=2

          do j = 1,iccp1
            litrmassrow(i,m,j)=0. 
            soilcmasrow(i,m,j)=0. 
          enddo

!         initial conditions always required
          dvdfcanrow(i,m,1)=1.0  !ndl
          dvdfcanrow(i,m,3)=1.0  !bdl
          dvdfcanrow(i,m,6)=1.0  !crop
          dvdfcanrow(i,m,8)=1.0  !grasses

!         then adjusted below for the actual mosaic makeup
          if (m .le. 2) then                     !ndl
           FCANROT(i,m,1)=1.0
           if (m .eq. 2) then
             dvdfcanrow(i,m,1)=0.0
             dvdfcanrow(i,m,2)=1.0        
           endif
          elseif (m .ge. 3 .and. m .le. 5) then  !bdl
           FCANROT(i,m,2)=1.0
           if (m .eq. 4) then
             dvdfcanrow(i,m,3)=0.0
             dvdfcanrow(i,m,4)=1.0        
           endif
           if (m .eq. 5) then
             dvdfcanrow(i,m,3)=0.0
             dvdfcanrow(i,m,5)=1.0        
           endif
          elseif (m .eq. 6 .or. m .eq. 7) then  !crop
           FCANROT(i,m,3)=1.0
           if (m .eq. 7) then
             dvdfcanrow(i,m,6)=0.0
             dvdfcanrow(i,m,7)=1.0        
           endif
          elseif (m .eq. 8 .or. m .eq. 9) then  !grasses
           FCANROT(i,m,4)=1.0
           if (m .eq. 9) then
             dvdfcanrow(i,m,8)=0.0
             dvdfcanrow(i,m,9)=1.0        
           endif
          else                                  !bare/urban? 
           FCANROT(i,m,5)=1.0
           endif !mosaic adjustments
         enddo  !m
        enddo  !i


         do i=1,nltest
          FAREROT(i,6)=crop_temp_frac(i,1)
          FAREROT(i,7)=crop_temp_frac(i,2)
         end do

      else if (.not. mosaic) then  !set up for composite runs when start_bare is on and compete or landuseon

!       store the read-in crop fractions as we keep them even when we start bare. 
!       FLAG: this is setup assuming that crops are in pft number 6 and 7. JM Apr 9 2014.
         do i=1,nltest
          crop_temp_frac(i,1)=FCANROT(i,1,3)*dvdfcanrow(i,1,6)
          crop_temp_frac(i,2)=FCANROT(i,1,3)*dvdfcanrow(i,1,7)
         end do

!      initalize to zero, these will be filled in by the luc or 
!      competition subroutines.
       FCANROT=0.0
       dvdfcanrow=0.0

       ! Added this as start_bare runs were not properly assigning 
       ! a TCAN on the very first day since the FCANROT was 0. JM Jan 14 2014. 
       do i=1,nltest
        do j=1,iccp1       
           if (j .lt. icc+1) then
            FCANROT(i,1,j)=seed
           else
            FCANROT(i,1,j)=1.0 - (real(icc) * seed)
           endif
        end do
       end do

       do i=1,nltest

!      initial conditions always required
         dvdfcanrow(i,1,1)=1.0  !ndl
         dvdfcanrow(i,1,3)=1.0  !bdl
         dvdfcanrow(i,1,6)=1.0  !crop
         dvdfcanrow(i,1,8)=1.0  !grasses

         do j = 1,icc
           ailcminrow(i,1,j)=0.0
           ailcmaxrow(i,1,j)=0.0
           gleafmasrow(i,1,j)=0.0
           bleafmasrow(i,1,j)=0.0
           stemmassrow(i,1,j)=0.0
           rootmassrow(i,1,j)=0.0
           lfstatusrow(i,1,j)=4
           pandaysrow(i,1,j)=0
         enddo

         lfstatusrow(i,1,1)=2

         do j = 1,iccp1
           litrmassrow(i,1,j)=0.0 
           soilcmasrow(i,1,j)=0.0 
         enddo
       enddo !nltest

         do i=1,nltest
          FCANROT(i,1,3) = crop_temp_frac(i,1) + crop_temp_frac(i,2)
          if (FCANROT(i,1,3) .gt. abszero) then
           dvdfcanrow(i,1,6) = crop_temp_frac(i,1) / FCANROT(i,1,3)
           dvdfcanrow(i,1,7) = crop_temp_frac(i,2) / FCANROT(i,1,3)
          else
           dvdfcanrow(i,1,6) = 1.0
           dvdfcanrow(i,1,7) = 0.0
          end if
         end do

      end if ! mosaic / composite
      end if !if (compete/landuseon .and. start_bare) 

end subroutine read_from_ctm

!==============================================================================================================

subroutine write_ctm_rs(nltest,nmtest,FCANROT,argbuff)

!   After a set period is complete the restart file for CTEM (.CTM_RS) is written
!   this restart file contains all of the CTEM level information needed to 
!   to restart the model to the same state.

use ctem_params,        only : ican,l2max,modelpft,icc,nmos,nlat,icp1,iccp1
use ctem_statevars,     only : c_switch,vrot,vgat

implicit none

! arguments:
character(80), intent(in) :: argbuff
integer, intent(in) :: nltest
integer, intent(inout) :: nmtest
real, dimension(nlat,nmos,icp1), intent(inout) :: FCANROT

! pointers:
logical, pointer :: lnduseon
logical, pointer :: compete
logical, pointer :: dowetlands
character(80), pointer :: titlec1
character(80), pointer :: titlec2
character(80), pointer :: titlec3
integer, pointer, dimension(:,:) :: icountrow
real, pointer, dimension(:,:,:) :: dvdfcanrow           !
real, pointer, dimension(:,:,:) :: fcancmxrow
real, pointer, dimension(:,:,:) :: ailcminrow           !
real, pointer, dimension(:,:,:) :: ailcmaxrow           !
real, pointer, dimension(:,:,:) :: gleafmasrow          !
real, pointer, dimension(:,:,:) :: bleafmasrow          !
real, pointer, dimension(:,:,:) :: stemmassrow          !
real, pointer, dimension(:,:,:) :: rootmassrow          !
real, pointer, dimension(:,:,:) :: litrmassrow
real, pointer, dimension(:,:,:) :: soilcmasrow
integer, pointer, dimension(:,:,:) :: lfstatusrow
integer, pointer, dimension(:,:,:) :: pandaysrow
real, pointer, dimension(:) :: extnprobgrd
real, pointer, dimension(:) :: prbfrhucgrd
real, pointer, dimension(:,:) :: mlightnggrd  
integer, pointer, dimension(:) :: stdalngrd
real, pointer, dimension(:) :: twarmm                   ! temperature of the warmest month (c)
real, pointer, dimension(:) :: tcoldm                   ! temperature of the coldest month (c)
real, pointer, dimension(:) :: gdd5                     ! growing degree days above 5 c
real, pointer, dimension(:) :: aridity                  ! aridity index, ratio of potential evaporation to precipitation
real, pointer, dimension(:) :: srplsmon                 ! number of months in a year with surplus water i.e.precipitation more than potential evaporation
real, pointer, dimension(:) :: defctmon                 ! number of months in a year with water deficit i.e.precipitation less than potential evaporation
real, pointer, dimension(:) :: anndefct                 ! annual water deficit (mm) 
real, pointer, dimension(:) :: annsrpls                 ! annual water surplus (mm)
real, pointer, dimension(:) :: annpcp                   ! annual precipitation (mm)
real, pointer, dimension(:) :: dry_season_length        ! length of dry season (months)
real, pointer, dimension(:,:) :: wetfrac_sgrd

! local variables

integer :: i,m,j,strlen
integer :: k1c,k2c,n
real, dimension(icc) :: rnded_pft

! point pointers:
lnduseon          => c_switch%lnduseon
compete           => c_switch%compete
dowetlands        => c_switch%dowetlands
titlec1           => c_switch%titlec1
titlec2           => c_switch%titlec2
titlec3           => c_switch%titlec3
icountrow         => vrot%icount
dvdfcanrow        => vrot%dvdfcan
fcancmxrow        => vrot%fcancmx
ailcminrow        => vrot%ailcmin
ailcmaxrow        => vrot%ailcmax
gleafmasrow       => vrot%gleafmas
bleafmasrow       => vrot%bleafmas
stemmassrow       => vrot%stemmass
rootmassrow       => vrot%rootmass
litrmassrow       => vrot%litrmass
soilcmasrow       => vrot%soilcmas
lfstatusrow       => vrot%lfstatus
pandaysrow        => vrot%pandays
extnprobgrd       => vrot%extnprob
prbfrhucgrd       => vrot%prbfrhuc
mlightnggrd       => vrot%mlightng 
stdalngrd         => vrot%stdaln
twarmm            => vgat%twarmm
tcoldm            => vgat%tcoldm
gdd5              => vgat%gdd5
aridity           => vgat%aridity
srplsmon          => vgat%srplsmon
defctmon          => vgat%defctmon
anndefct          => vgat%anndefct
annsrpls          => vgat%annsrpls
annpcp            => vgat%annpcp
dry_season_length => vgat%dry_season_length
wetfrac_sgrd      => vgat%wetfrac_s
      
! -----------------      
! Begin

open(unit=101,file=argbuff(1:strlen(argbuff))//'.CTM_RS')

write(101,7010) titlec1
write(101,7010) titlec2
write(101,7010) titlec3

7010  FORMAT(A80)

! if landuseon or competition, then we need to recreate the dvdfcanrow so do so now
if (lnduseon .or. compete ) then
  icountrow=0
  do j = 1, ican
    do i = 1,nltest
        do m = 1,nmtest 
        k1c = (j-1)*l2max + 1
        k2c = k1c + (l2max - 1)
    
        do n = k1c, k2c
          if (modelpft(n) .eq. 1) then
            icountrow(i,m) = icountrow(i,m) + 1
            if (FCANROT(i,m,j) .gt. 0.) then
              dvdfcanrow(i,m,icountrow(i,m)) = fcancmxrow(i,m,icountrow(i,m))/FCANROT(i,m,j)
            else
              dvdfcanrow(i,m,icountrow(i,m)) = 0.
            end if
          end if !modelpft
        end do !n
    
        ! check to ensure that the dvdfcanrow's add up to 1 across a class-level pft
        if (dvdfcanrow(i,m,1) .eq. 0. .and. dvdfcanrow(i,m,2) .eq. 0.) then
            dvdfcanrow(i,m,1)=1.0
        else if (dvdfcanrow(i,m,3) .eq. 0. .and. dvdfcanrow(i,m,4) .eq. 0. .and. dvdfcanrow(i,m,5) .eq. 0.) then
            dvdfcanrow(i,m,3)=1.0
        else if (dvdfcanrow(i,m,6) .eq. 0. .and. dvdfcanrow(i,m,7) .eq. 0.) then
            dvdfcanrow(i,m,6)=1.0
        else if (dvdfcanrow(i,m,8) .eq. 0. .and. dvdfcanrow(i,m,9) .eq. 0.) then
            dvdfcanrow(i,m,8)=1.0
        end if

        end do !m
    enddo !i 
  enddo !j

  do i=1,nltest
   do m=1,nmtest 
    do j = 1, icc
!            lastly check if the different pfts accidently add up > 1.0
!            after rounding to the number of sig figs used in the output
!            this rounds to 3 decimal places. if you are found to be over
!            or under, arbitrarily reduce one of the pfts. the amount of
!            the change will be inconsequential. 
       rnded_pft(j) =real(int(dvdfcanrow(i,m,j) * 1000.0))/ 1000.0
       dvdfcanrow(i,m,j) = rnded_pft(j)
    end do

    if (dvdfcanrow(i,m,1) + dvdfcanrow(i,m,2) .ne. 1.0) then
        dvdfcanrow(i,m,1) = 1.0 - rnded_pft(2)
        dvdfcanrow(i,m,2) = rnded_pft(2)
    end if 
    if (dvdfcanrow(i,m,3) + dvdfcanrow(i,m,4) +  dvdfcanrow(i,m,5) .ne. 1.0) then
        dvdfcanrow(i,m,3) = 1.0 - rnded_pft(4) - rnded_pft(5)
        dvdfcanrow(i,m,4) = rnded_pft(4)
        dvdfcanrow(i,m,5) = rnded_pft(5)
    end if 
    if (dvdfcanrow(i,m,6) + dvdfcanrow(i,m,7) .ne. 1.0) then
        dvdfcanrow(i,m,6) = 1.0 - rnded_pft(7)
        dvdfcanrow(i,m,7) = rnded_pft(7)
    end if 
    if (dvdfcanrow(i,m,8) + dvdfcanrow(i,m,9) .ne. 1.0) then
        dvdfcanrow(i,m,8) = 1.0 - rnded_pft(9)
        dvdfcanrow(i,m,9) = rnded_pft(9)
    end if
   enddo
  enddo

end if !lnuse/compete

do i=1,nltest
    do m=1,nmtest
        write(101,7011) (ailcminrow(i,m,j),j=1,icc)
        write(101,7011) (ailcmaxrow(i,m,j),j=1,icc)
        write(101,'(9f8.3)') (dvdfcanrow(i,m,j),j=1,icc)
        write(101,7011) (gleafmasrow(i,m,j),j=1,icc)
        write(101,7011) (bleafmasrow(i,m,j),j=1,icc)
        write(101,7011) (stemmassrow(i,m,j),j=1,icc)
        write(101,7011) (rootmassrow(i,m,j),j=1,icc)
        write(101,7013) (litrmassrow(i,m,j),j=1,iccp1)
        write(101,7013) (soilcmasrow(i,m,j),j=1,iccp1)
        write(101,7012) (lfstatusrow(i,m,j),j=1,icc)
        write(101,7012) (pandaysrow(i,m,j),j=1,icc)
    end do !nmtest

    write(101,"(6f8.3)") (mlightnggrd(i,j),j=1,6)  !mean monthly lightning frequency
    write(101,"(6f8.3)") (mlightnggrd(i,j),j=7,12) !flashes/km2.year
    write(101,"(f8.2)") extnprobgrd(i)
    write(101,"(f8.2)") prbfrhucgrd(i)
    write(101,"(i4)") stdalngrd(i)

    if (compete) then
        write(101,"(5f8.2)")twarmm(i),tcoldm(i),gdd5(i),aridity(i),srplsmon(i)
        write(101,"(5f8.2)")defctmon(i),anndefct(i),annsrpls(i), annpcp(i),dry_season_length(i)
    end if

    if (dowetlands) then     
        write(101,"(8f9.5)")(wetfrac_sgrd(i,j),j=1,8)
    end if   

end do !nltest

close(101)

7011  format(9f8.2)
7012  format(9i8)
7013  format(10f8.2)

end subroutine write_ctm_rs        

!==============================================================================================================

subroutine create_outfiles(argbuff,title1, title2, title3, title4, title5, title6, name1, name2, name3, &
                           name4, name5, name6, place1 ,place2, place3, place4, place5, place6)
                           
!   All output files are initialized in this subroutine

use ctem_statevars,     only : c_switch,vrot,vgat

implicit none

! arguments:
character(80), intent(in) :: argbuff
character(4), intent(in) :: title1, title2, title3, title4, &
                            title5, title6, name1, name2, name3, &
                            name4, name5, name6, place1 ,place2, &
                            place3, place4, place5, place6

! pointers:
logical, pointer :: mosaic
logical, pointer :: dofire
logical, pointer :: ctem_on
logical, pointer :: compete
logical, pointer :: dowetlands
logical, pointer :: lnduseon
logical, pointer :: obswetf
logical, pointer :: parallelrun

! local variables:
integer :: strlen
character(80) :: titlec1, titlec2, titlec3

! point pointers:
mosaic            => c_switch%mosaic
dofire            => c_switch%dofire
ctem_on           => c_switch%ctem_on
compete           => c_switch%compete
dowetlands        => c_switch%dowetlands    
lnduseon          => c_switch%lnduseon
obswetf           => c_switch%obswetf
parallelrun       => c_switch%parallelrun

!-----     
! begin:

!      the ctem output file suffix naming convention is as follows:
!                       ".CT##{time}_{mosaic/grid}"
!      where the ## is a numerical identifier, {time} is any of H, D, M,
!      or Y for half hourly, daily, monthly, or yearly, respectively. 
!      after the underscore M or G is used to denote mosaic or grid 
!      -averaged values, respectively. also possible is GM for competition
!      outputs since they are the same format in either composite or 
!      mosaic modes. 
     
6001  FORMAT('CLASS-CTEM TEST RUN:     ',6A4)
6002  FORMAT('RESEARCHER:         ',6A4)
6003  FORMAT('INSTITUTION:        ',6A4)

if (.not. parallelrun .and. ctem_on) then ! stand alone mode, includes half-hourly and daily output

    ! ctem half hourly output files
    open(unit=71, file=argbuff(1:strlen(argbuff))//'.CT01H_M')  
    open(unit=711,file=argbuff(1:strlen(argbuff))//'.CT01H_G')

    if (mosaic) then
    
    !   ctem daily output files (mosaic)
        open(unit=72,file=argbuff(1:strlen(argbuff))//'.CT01D_M')  
        open(unit=73,file=argbuff(1:strlen(argbuff))//'.CT02D_M')
        open(unit=74,file=argbuff(1:strlen(argbuff))//'.CT03D_M')
        open(unit=75,file=argbuff(1:strlen(argbuff))//'.CT04D_M')
        open(unit=76,file=argbuff(1:strlen(argbuff))//'.CT05D_M')
        
        if (dofire .or. lnduseon) then
        open(unit=78,file=argbuff(1:strlen(argbuff))//'.CT06D_M') ! disturbance vars
        endif

    end if ! mosaic

    ! ctem daily output files (grid-average)
    open(unit=721,file=argbuff(1:strlen(argbuff))//'.CT01D_G') 
    open(unit=731,file=argbuff(1:strlen(argbuff))//'.CT02D_G')
    open(unit=741,file=argbuff(1:strlen(argbuff))//'.CT03D_G')
    open(unit=751,file=argbuff(1:strlen(argbuff))//'.CT04D_G')

    if (dofire .or. lnduseon) then
        open(unit=781,file=argbuff(1:strlen(argbuff))//'.CT06D_G') ! disturbance vars
    endif

    if (compete .or. lnduseon) then
        open(unit=761,file=argbuff(1:strlen(argbuff))//'.CT07D_G') ! competition
    end if

    if (dowetlands .or. obswetf) then
        open(unit=762,file=argbuff(1:strlen(argbuff))//'.CT08D_G') ! Methane(Wetland)
    endif 

endif ! parallelrun & ctem_on

! monthly & yearly output for both parallel mode and stand alone mode

! CLASS MONTHLY OUTPUT FILES 
OPEN(UNIT=81,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF1M_G') 
OPEN(UNIT=82,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF2M_G')

! CLASS YEARLY OUTPUT FILES
OPEN(UNIT=83,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF1Y_G')

if (ctem_on) then

    if (mosaic) then

        ! Mosaic file names:

        ! CTEM monthly output files
        open(unit=84,file=argbuff(1:strlen(argbuff))//'.CT01M_M')

        ! CTEM yearly output files
        open(unit=86,file=argbuff(1:strlen(argbuff))//'.CT01Y_M')

        if (dofire .or. lnduseon) then
            open(unit=85,file=argbuff(1:strlen(argbuff))//'.CT06M_M') ! Monthly disturbance
            open(unit=87,file=argbuff(1:strlen(argbuff))//'.CT06Y_M') ! Annual disturbance
        endif 

    else

        ! Composite file names:

        ! CTEM monthly output files
        open(unit=84,file=argbuff(1:strlen(argbuff))//'.CT01M_G')

        ! CTEM yearly output files
        open(unit=86,file=argbuff(1:strlen(argbuff))//'.CT01Y_G')

        if (dofire .or. lnduseon) then
            open(unit=85,file=argbuff(1:strlen(argbuff))//'.CT06M_G') ! Monthly disturbance
            open(unit=87,file=argbuff(1:strlen(argbuff))//'.CT06Y_G') ! Annual disturbance
        endif 

    end if !mosiac/composite 
    
    if (compete .or. lnduseon) then
        open(unit=88,file=argbuff(1:strlen(argbuff))//'.CT07M_GM')! ctem pft fractions MONTHLY

        open(unit=89,file=argbuff(1:strlen(argbuff))//'.CT07Y_GM')! ctem pft fractions YEARLY
    endif
    
    if (dowetlands .or. obswetf) then
        open(unit=91,file=argbuff(1:strlen(argbuff))//'.CT08M_G')  !Methane(wetland) MONTHLY
        
        open(unit=92,file=argbuff(1:strlen(argbuff))//'.CT08Y_G')  !Methane(wetland) YEARLY
    endif !dowetlands
    
end if !ctem_on

!===========================
!
!     CTEM FILE TITLES
!
if (ctem_on .and. .not. parallelrun) then
    write(71,6001) title1,title2,title3,title4,title5,title6
    write(71,6002) name1,name2,name3,name4,name5,name6
    write(71,6003) place1,place2,place3,place4,place5,place6
    write(71,7020)
    write(71,7030)
    
    if (mosaic) then
        write(72,6001) title1,title2,title3,title4,title5,title6
        write(72,6002) name1,name2,name3,name4,name5,name6
        write(72,6003) place1,place2,place3,place4,place5,place6
        write(72,7020)
        write(72,7040)
        
        write(73,6001) title1,title2,title3,title4,title5,title6
        write(73,6002) name1,name2,name3,name4,name5,name6
        write(73,6003) place1,place2,place3,place4,place5,place6
        write(73,7020)
        write(73,7050)
        
        write(74,6001) title1,title2,title3,title4,title5,title6
        write(74,6002) name1,name2,name3,name4,name5,name6
        write(74,6003) place1,place2,place3,place4,place5,place6
        write(74,7020)
        write(74,7061)
        
        write(75,6001) title1,title2,title3,title4,title5,title6
        write(75,6002) name1,name2,name3,name4,name5,name6
        write(75,6003) place1,place2,place3,place4,place5,place6
        write(75,7020)
        write(75,7070)
        
        write(76,6001) title1,title2,title3,title4,title5,title6
        write(76,6002) name1,name2,name3,name4,name5,name6
        write(76,6003) place1,place2,place3,place4,place5,place6
        write(76,7020)
        write(76,7080)
    
        if (dofire .or. lnduseon) then
            write(78,6001) title1,title2,title3,title4,title5,title6
            write(78,6002) name1,name2,name3,name4,name5,name6
            write(78,6003) place1,place2,place3,place4,place5,place6
            write(78,7021)
            write(78,7110)
            write(78,7111)
        end if
        
    end if !mosaic

7010  FORMAT(A80)
7020  FORMAT('CANADIAN TERRESTRIAL ECOSYSTEM MODEL (CTEM) DAILY RESULTS')
7021  FORMAT('CANADIAN TERRESTRIAL ECOSYSTEM MODEL (CTEM) DAILY ',' DISTURBANCE RESULTS')
7030  FORMAT('HOUR MIN  DAY YEAR, An FOR 9 PFTs, RmL FOR 9 PFTs')
7040  FORMAT('  DAY YEAR       GPP       NPP       NEP       NBP', '   AUTORES  HETRORES    LITRES    SOCRES  DSTCEMLS  LITRFALL', &
     '  HUMIFTRS')
7050  FORMAT('  DAY YEAR       RML       RMS       RMR        RG','  LEAFLITR  TLTRLEAF  TLTRSTEM  TLTRROOT ')
7060  FORMAT('  DAY YEAR  VGBIOMAS   GAVGLAI  GAVGLTMS  GAVGSCMS  ','TOTCMASS  GLEAFMAS   BLEAFMAS STEMMASS   ROOTMASS  LITRMASS ',' SOILCMAS')
7061  FORMAT('  DAY YEAR  VGBIOMAS   GAVGLAI  GLEAFMAS   BLEAFMAS ', 'STEMMASS   ROOTMASS  LITRMASS SOILCMAS')
7070  FORMAT('  DAY YEAR     AILCG     AILCB    RMATCTEM ','LAYER 1,2, & 3     VEGHGHT  ROOTDPTH  ROOTTEMP      SLAI')
7075  FORMAT('  DAY YEAR   FRAC #1   FRAC #2   FRAC #3   FRAC #4   ','FRAC #5   FRAC #6   FRAC #7   FRAC #8   FRAC #9  ','FRAC #10[%] SUMCHECK')
7080  FORMAT('  DAY YEAR   AFRLEAF   AFRSTEM   AFRROOT  TCANOACC','  LFSTATUS')
7110  FORMAT('  DAY YEAR   EMIT_CO2','    EMIT_CO   EMIT_CH4  EMIT_NMHC    EMIT_H2   EMIT_NOX', &
            '   EMIT_N2O  EMIT_PM25   EMIT_TPM    EMIT_TC    EMIT_OC','    EMIT_BC   BURNFRAC   PROBFIRE   LUCEMCOM   LUCLTRIN',&
            '   LUCSOCIN   GRCLAREA   BTERM   LTERM   MTERM')
7111  FORMAT('               g/m2.D     g/m2.d','     g/m2.d     g/m2.d     g/m2.d     g/m2.d     g/m2.d',&
            '     g/m2.d     g/m2.d     g/m2.d     g/m2.d     g/m2.d   ','       %  avgprob/d uMOL-CO2/M2.S KgC/M2.D','   KgC/M2.D      KM^2    prob/d       prob/d       prob/d')
7112  FORMAT(' DAY  YEAR   CH4WET1    CH4WET2    WETFDYN   CH4DYN1  CH4DYN2 ')
7113  FORMAT('          umolCH4/M2.S    umolCH4/M2.S          umolCH4/M2.S  umolCH4/M2.S')

    write(711,6001) title1,title2,title3,title4,title5,title6
    write(711,6002) name1,name2,name3,name4,name5,name6
    write(711,6003) place1,place2,place3,place4,place5,place6
    write(711,7020)
    write(711,7030)

    write(721,6001) title1,title2,title3,title4,title5,title6
    write(721,6002) name1,name2,name3,name4,name5,name6
    write(721,6003) place1,place2,place3,place4,place5,place6
    write(721,7020)
    write(721,7040)

    write(731,6001) title1,title2,title3,title4,title5,title6
    write(731,6002) name1,name2,name3,name4,name5,name6
    write(731,6003) place1,place2,place3,place4,place5,place6
    write(731,7020)
    write(731,7050)

    write(741,6001) title1,title2,title3,title4,title5,title6
    write(741,6002) name1,name2,name3,name4,name5,name6
    write(741,6003) place1,place2,place3,place4,place5,place6
    write(741,7020)
    write(741,7060)

    write(751,6001) title1,title2,title3,title4,title5,title6
    write(751,6002) name1,name2,name3,name4,name5,name6
    write(751,6003) place1,place2,place3,place4,place5,place6
    write(751,7020)
    write(751,7070)

    if (compete .or. lnduseon) then
        write(761,6001) title1,title2,title3,title4,title5,title6
        write(761,6002) name1,name2,name3,name4,name5,name6
        write(761,6003) place1,place2,place3,place4,place5,place6
        write(761,7020)
        write(761,7075)
    end if
    
    if (dofire .or. lnduseon) then
        write(781,6001) title1,title2,title3,title4,title5,title6
        write(781,6002) name1,name2,name3,name4,name5,name6
        write(781,6003) place1,place2,place3,place4,place5,place6
        write(781,7020)
        write(781,7110)
        write(781,7111)
    end if
    ! methane(wetland) variables 
    if (dowetlands .or. obswetf) then     
        write(762,6001) title1,title2,title3,title4,title5,title6
        write(762,6002) name1,name2,name3,name4,name5,name6
        write(762,6003) place1,place2,place3,place4,place5,place6
        write(762,7020)
        write(762,7112)
        write(762,7113)
    end if 

end if !ctem_on & not parallelrun
 
!     CLASS MONTHLY & YEARLY OUTPUT FOR BOTH PARALLEL MODE AND STAND ALONE MODE
WRITE(81,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
WRITE(81,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
WRITE(81,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
WRITE(81,6021) 
WRITE(81,6121)

WRITE(82,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
WRITE(82,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
WRITE(82,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
WRITE(82,6022)
WRITE(82,6122)

WRITE(83,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
WRITE(83,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
WRITE(83,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
WRITE(83,6023)
WRITE(83,6123)

if (ctem_on) then

    write(84,6001) title1,title2,title3,title4,title5,title6
    write(84,6002) name1,name2,name3,name4,name5,name6
    write(84,6003) place1,place2,place3,place4,place5,place6
    write(84,6024)
    write(84,6124)
    write(84,6224)
    
    if (dofire .or. lnduseon) then
        write(85,6001) title1,title2,title3,title4,title5,title6
        write(85,6002) name1,name2,name3,name4,name5,name6
        write(85,6003) place1,place2,place3,place4,place5,place6
        write(85,6025)
        write(85,6125)
        write(85,6225)
    end if

    write(86,6001) title1,title2,title3,title4,title5,title6
    write(86,6002) name1,name2,name3,name4,name5,name6
    write(86,6003) place1,place2,place3,place4,place5,place6
    write(86,6026)
    write(86,6126)
    write(86,6226)

    if (dofire .or. lnduseon) then        
        write(87,6001) title1,title2,title3,title4,title5,title6
        write(87,6002) name1,name2,name3,name4,name5,name6
        write(87,6003) place1,place2,place3,place4,place5,place6
        write(87,6027)
        write(87,6127)
        write(87,6227)
    end if

    if (compete .or. lnduseon) then
        write(88,6001) title1,title2,title3,title4,title5,title6
        write(88,6002) name1,name2,name3,name4,name5,name6
        write(88,6003) place1,place2,place3,place4,place5,place6
        write(88,6028)
        write(88,6128)
        write(88,6228)

        write(89,6001) title1,title2,title3,title4,title5,title6
        write(89,6002) name1,name2,name3,name4,name5,name6
        write(89,6003) place1,place2,place3,place4,place5,place6
        write(89,6029)
        write(89,6129)
        write(89,6229)
    end if !compete

    if (dowetlands .or. obswetf) then
        write(91,6001) title1,title2,title3,title4,title5,title6
        write(91,6002) name1,name2,name3,name4,name5,name6
        write(91,6003) place1,place2,place3,place4,place5,place6
        write(91,6028)
        write(91,6230)
        write(91,6231)

        write(92,6001) title1,title2,title3,title4,title5,title6
        write(92,6002) name1,name2,name3,name4,name5,name6
        write(92,6003) place1,place2,place3,place4,place5,place6
        write(92,6028)
        write(92,6232)
        write(92,6233)
    end if 
    
end if !ctem_on & parallelrun

6021  FORMAT(2X,'MONTH YEAR  SW     LW      QH      QE    SNOACC    ','WSNOACC    ROFACC      PCP      EVAP       TAIR')
6121  FORMAT(2X,'           W/m2    W/m2    W/m2    W/m2    kg/m2   ','kg/m2      mm.mon    mm.mon    mm.mon      degC') 
6022  FORMAT(2X,'MONTH  YEAR  TG1  THL1  THI1     TG2  THL2  THI2','     TG3  THL3  THI3')
6122  FORMAT(2X,'             deg  m3/m3  m3/m3   deg  m3/m3  ','m3/m3   deg  m3/m3  m3/m3')
6023  FORMAT(2X,'YEAR   SW     LW      QH      QE     ROFACC   ',' PCP     EVAP  ')
6123  FORMAT(2X,'      W/m2   W/m2    W/m2    W/m2    mm.yr    ','mm.yr    mm.yr')
6024  FORMAT('CANADIAN TERRESTRIAL ECOSYSTEM MODEL (CTEM) MONTHLY ','RESULTS')
6124  FORMAT('  MONTH  YEAR  LAIMAXG  VGBIOMAS  LITTER    SOIL_C  ', '  NPP       GPP        NEP       NBP    HETRES','   AUTORES    LITRES   SOILCRES')
6224  FORMAT('                 m2/m2  Kg C/m2  Kg C/m2   Kg C/m2  ','gC/m2.mon  gC/m2.mon  gC/m2.mon  g/m2.mon   g/m2.mon ',&
            'gC/m2.mon  gC/m2.mon  gC/m2.mon')   
6025  FORMAT('CANADIAN TERRESTRIAL ECOSYSTEM MODEL (CTEM) MONTHLY ','RESULTS FOR DISTURBANCES')
6125  FORMAT('  MONTH  YEAR  CO2','        CO        CH4      NMHC       H2       NOX       N2O       PM25       TPM        TC        OC        BC  ',&
            ' PROBFIRE  LUC_CO2_E  LUC_LTRIN  LUC_SOCIN   BURNFRAC    BTERM',' LTERM   MTERM')
6225  FORMAT('            g/m2.mon  g/m2.mon','  g/m2.mon  g/m2.mon  g/m2.mon  g/m2.mon  g/m2.mon', &
            '  g/m2.mon  g/m2.mon  g/m2.mon  g/m2.mon  g/m2.mon','  prob/mon    g C/m2    g C/m2    g C/m2         %  prob/mon','  prob/mon  prob/mon')  
6026  FORMAT('CANADIAN TERRESTRIAL ECOSYSTEM MODEL (CTEM) YEARLY ','RESULTS')
6126  FORMAT('  YEAR   LAIMAXG  VGBIOMAS  STEMMASS  ROOTMASS  LITRMASS', '  SOILCMAS  TOTCMASS  ANNUALNPP ANNUALGPP ANNUALNEP ANNUALNBP',&
     ' ANNHETRSP ANAUTORSP ANNLITRES ANSOILCRES')
6226  FORMAT('          m2/m2   Kg C/m2   Kg C/m2   Kg C/m2    Kg C/m2','  Kg C/m2   Kg C/m2   gC/m2.yr  gC/m2.yr  gC/m2.yr  gC/m2.yr',&
     '  gC/m2.yr  gC/m2.yr  gC/m2.yr  gC/m2.yr')
6027  FORMAT('CANADIAN TERRESTRIAL ECOSYSTEM MODEL (CTEM) YEARLY ','RESULTS FOR DISTURBANCES')
6127  FORMAT('  YEAR   ANNUALCO2','  ANNUALCO  ANNUALCH4  ANN_NMHC ANNUAL_H2 ANNUALNOX ANNUALN2O','  ANN_PM25  ANNUALTPM ANNUAL_TC ANNUAL_OC ANNUAL_BC APROBFIRE',&
     ' ANNLUCCO2  ANNLUCLTR ANNLUCSOC ABURNFRAC ANNBTERM ANNLTERM',' ANNMTERM')
6227  FORMAT('         g/m2.yr','  g/m2.yr  g/m2.yr  g/m2.yr  g/m2.yr  g/m2.yr  g/m2.yr','  g/m2.yr  g/m2.yr  g/m2.yr  g/m2.yr  g/m2.yr  prob/yr ',&
            '  g/m2.yr  g/m2.yr  g/m2.yr    %     prob/yr  prob/yr','  prob/yr')
6028  FORMAT('CANADIAN TERRESTRIAL ECOSYSTEM MODEL (CTEM) MONTHLY ','RESULTS')
6128  FORMAT(' MONTH YEAR  FRAC #1   FRAC #2   FRAC #3   FRAC #4   ','FRAC #5   FRAC #6   FRAC #7   FRAC #8   FRAC #9   FRAC #10   ',&
            'SUMCHECK, PFT existence for each of the 9 pfts') 
6228  FORMAT('             %         %         %         %         ','%         %         %         %         %         %          ','     ')   
6029  FORMAT('CANADIAN TERRESTRIAL ECOSYSTEM MODEL (CTEM) YEARLY ','RESULTS')
6129  FORMAT('  YEAR   FRAC #1   FRAC #2   FRAC #3   FRAC #4   ','FRAC #5   FRAC #6   FRAC #7   FRAC #8   FRAC #9   FRAC #10   ',&
            'SUMCHECK, PFT existence for each of the 9 pfts')
6229  FORMAT('         %         %         %         %         ','%         %         %         %         %         %          ','     ')
     
6230  FORMAT('MONTH  YEAR   CH4WET1    CH4WET2    WETFDYN   CH4DYN1  CH4DYN2 ')
6231  FORMAT('       gCH4/M2.MON     gCH4/M2.MON        gCH4/M2.MON  gCH4/M2.MON')
6232  FORMAT('  YEAR   CH4WET1    CH4WET2    WETFDYN   CH4DYN1  CH4DYN2 ')
6233  FORMAT('      gCH4/M2.YR      gCH4/M2.YR         gCH4/M2.YR   gCH4/M2.YR ')
 
end subroutine create_outfiles

!==============================================================================================================


!==============================================================================================================

subroutine class_monthly_aw(IDAY,IYEAR,NCOUNT,NDAY,SBC,DELT,nltest,nmtest,&
                                  ALVSROT,FAREROT,FSVHROW,ALIRROT,FSIHROW,GTROT,FSSROW, &
                                  FDLROW,HFSROT,ROFROT,PREROW,QFSROT,QEVPROT,SNOROT, &
                                  TAROW,WSNOROT,TBARROT,THLQROT,THICROT,TFREZ)
                           
use ctem_statevars,     only : class_out,resetclassmon
use ctem_params, only : nmon, monthend, nlat, nmos, ignd

implicit none

! arguments
integer, intent(in) :: IDAY
integer, intent(in) :: IYEAR
integer, intent(in) :: NCOUNT
integer, intent(in) :: NDAY
integer, intent(in) :: nltest
integer, intent(in) :: nmtest
real, intent(in) :: SBC
real, intent(in) :: DELT
real, intent(in) :: TFREZ
real, dimension(nlat), intent(in) :: FSSROW
real, dimension(nlat), intent(in) :: FDLROW
real, dimension(nlat), intent(in) :: FSVHROW
real, dimension(nlat), intent(in) :: FSIHROW
real, dimension(nlat), intent(in) :: TAROW
real, dimension(nlat), intent(in) :: PREROW
real, dimension(nlat,nmos), intent(in) :: ALVSROT
real, dimension(nlat,nmos), intent(in) :: FAREROT
real, dimension(nlat,nmos), intent(in) :: ALIRROT
real, dimension(nlat,nmos), intent(in) :: GTROT
real, dimension(nlat,nmos), intent(in) :: HFSROT
real, dimension(nlat,nmos), intent(in) :: QEVPROT
real, dimension(nlat,nmos), intent(in) :: SNOROT
real, dimension(nlat,nmos), intent(in) :: WSNOROT
real, dimension(nlat,nmos), intent(in) :: ROFROT
real, dimension(nlat,nmos), intent(in) :: QFSROT
real, dimension(nlat,nmos,ignd), intent(in) :: TBARROT
real, dimension(nlat,nmos,ignd), intent(in) :: THLQROT
real, dimension(nlat,nmos,ignd), intent(in) :: THICROT

! pointers
real, pointer, dimension(:) :: ALVSACC_MO
real, pointer, dimension(:) :: ALIRACC_MO
real, pointer, dimension(:) :: FLUTACC_MO
real, pointer, dimension(:) :: FSINACC_MO
real, pointer, dimension(:) :: FLINACC_MO
real, pointer, dimension(:) :: HFSACC_MO
real, pointer, dimension(:) :: QEVPACC_MO
real, pointer, dimension(:) :: SNOACC_MO
real, pointer, dimension(:) :: WSNOACC_MO
real, pointer, dimension(:) :: ROFACC_MO
real, pointer, dimension(:) :: PREACC_MO
real, pointer, dimension(:) :: EVAPACC_MO
real, pointer, dimension(:) :: TAACC_MO
real, pointer :: FSSTAR_MO
real, pointer :: FLSTAR_MO
real, pointer :: QH_MO
real, pointer :: QE_MO
real, pointer, dimension(:,:) :: TBARACC_MO
real, pointer, dimension(:,:) :: THLQACC_MO
real, pointer, dimension(:,:) :: THICACC_MO
   
! local
real :: ALTOT_MO
integer :: NT
integer :: NDMONTH
integer :: i,m,j
integer :: IMONTH

! point pointers
ALVSACC_MO        => class_out%ALVSACC_MO  
ALIRACC_MO        => class_out%ALIRACC_MO
FLUTACC_MO        => class_out%FLUTACC_MO
FSINACC_MO        => class_out%FSINACC_MO
FLINACC_MO        => class_out%FLINACC_MO
HFSACC_MO         => class_out%HFSACC_MO
QEVPACC_MO        => class_out%QEVPACC_MO
SNOACC_MO         => class_out%SNOACC_MO
WSNOACC_MO        => class_out%WSNOACC_MO
ROFACC_MO         => class_out%ROFACC_MO
PREACC_MO         => class_out%PREACC_MO
EVAPACC_MO        => class_out%EVAPACC_MO
TAACC_MO          => class_out%TAACC_MO
FSSTAR_MO         => class_out%FSSTAR_MO
FLSTAR_MO         => class_out%FLSTAR_MO
QH_MO             => class_out%QH_MO
QE_MO             => class_out%QE_MO  
TBARACC_MO        => class_out%TBARACC_MO
THLQACC_MO        => class_out%THLQACC_MO    
THICACC_MO        => class_out%THICACC_MO

! ------------

! Accumulate output data for monthly averaged fields for class grid-mean.
! for both parallel mode and stand alone mode

FSSTAR_MO   =0.0
FLSTAR_MO   =0.0
QH_MO       =0.0
QE_MO       =0.0
ALTOT_MO    =0.0

DO 820 I=1,NLTEST
DO 821 M=1,NMTEST
    ALVSACC_MO(I)=ALVSACC_MO(I)+ALVSROT(I,M)*FAREROT(I,M)*FSVHROW(I)
    ALIRACC_MO(I)=ALIRACC_MO(I)+ALIRROT(I,M)*FAREROT(I,M)*FSIHROW(I) 
    FLUTACC_MO(I)=FLUTACC_MO(I)+SBC*GTROT(I,M)**4*FAREROT(I,M)
    FSINACC_MO(I)=FSINACC_MO(I)+FSSROW(I)*FAREROT(I,M)
    FLINACC_MO(I)=FLINACC_MO(I)+FDLROW(I)*FAREROT(I,M)
    HFSACC_MO(I) =HFSACC_MO(I)+HFSROT(I,M)*FAREROT(I,M)
    QEVPACC_MO(I)=QEVPACC_MO(I)+QEVPROT(I,M)*FAREROT(I,M)
    SNOACC_MO(I) =SNOACC_MO(I)+SNOROT(I,M)*FAREROT(I,M)
    TAACC_MO(I)=TAACC_MO(I)+TAROW(I)*FAREROT(I,M)

    IF(SNOROT(I,M).GT.0.0) THEN
        WSNOACC_MO(I)=WSNOACC_MO(I)+WSNOROT(I,M)*FAREROT(I,M)
    ENDIF

    ROFACC_MO(I) =ROFACC_MO(I)+ROFROT(I,M)*FAREROT(I,M)*DELT
    PREACC_MO(I) =PREACC_MO(I)+PREROW(I)*FAREROT(I,M)*DELT
    EVAPACC_MO(I)=EVAPACC_MO(I)+QFSROT(I,M)*FAREROT(I,M)*DELT

    DO 823 J=1,IGND
        TBARACC_MO(I,J)=TBARACC_MO(I,J)+TBARROT(I,M,J)*FAREROT(I,M)
        THLQACC_MO(I,J)=THLQACC_MO(I,J)+THLQROT(I,M,J)*FAREROT(I,M)
        THICACC_MO(I,J)=THICACC_MO(I,J)+THICROT(I,M,J)*FAREROT(I,M)
823 CONTINUE

821    CONTINUE
820   CONTINUE 

DO NT=1,NMON
    IF(IDAY.EQ.monthend(NT+1).AND.NCOUNT.EQ.NDAY)THEN
    
        IMONTH=NT
        NDMONTH=(monthend(NT+1)-monthend(NT))*NDAY

        DO 824 I=1,NLTEST
            IF(FSINACC_MO(I).GT.0.0) THEN
                ALVSACC_MO(I)=ALVSACC_MO(I)/(FSINACC_MO(I)*0.5)
                ALIRACC_MO(I)=ALIRACC_MO(I)/(FSINACC_MO(I)*0.5)
            ELSE
                ALVSACC_MO(I)=0.0
                ALIRACC_MO(I)=0.0
            ENDIF
            
            FLUTACC_MO(I)=FLUTACC_MO(I)/REAL(NDMONTH)
            FSINACC_MO(I)=FSINACC_MO(I)/REAL(NDMONTH)
            FLINACC_MO(I)=FLINACC_MO(I)/REAL(NDMONTH)
            HFSACC_MO(I) =HFSACC_MO(I)/REAL(NDMONTH)
            QEVPACC_MO(I)=QEVPACC_MO(I)/REAL(NDMONTH)
            SNOACC_MO(I) =SNOACC_MO(I)/REAL(NDMONTH)
            WSNOACC_MO(I)=WSNOACC_MO(I)/REAL(NDMONTH)
            ROFACC_MO(I) =ROFACC_MO(I)
            PREACC_MO(I) =PREACC_MO(I)
            EVAPACC_MO(I)=EVAPACC_MO(I)
            TAACC_MO(I)=TAACC_MO(I)/REAL(NDMONTH)
            DO J=1,IGND
                TBARACC_MO(I,J)=TBARACC_MO(I,J)/REAL(NDMONTH)
                THLQACC_MO(I,J)=THLQACC_MO(I,J)/REAL(NDMONTH)
                THICACC_MO(I,J)=THICACC_MO(I,J)/REAL(NDMONTH)
            ENDDO

            ALTOT_MO=(ALVSACC_MO(I)+ALIRACC_MO(I))/2.0
            FSSTAR_MO=FSINACC_MO(I)*(1.-ALTOT_MO)
            FLSTAR_MO=FLINACC_MO(I)-FLUTACC_MO(I)
            QH_MO=HFSACC_MO(I)
            QE_MO=QEVPACC_MO(I)

            WRITE(81,8100)IMONTH,IYEAR,FSSTAR_MO,FLSTAR_MO,QH_MO, &
                         QE_MO,SNOACC_MO(I),WSNOACC_MO(I), &
                         ROFACC_MO(I),PREACC_MO(I),EVAPACC_MO(I), &
                         TAACC_MO(I)-TFREZ
            IF (IGND.GT.3) THEN
            WRITE(82,8101)IMONTH,IYEAR,(TBARACC_MO(I,J)-TFREZ, &
                          THLQACC_MO(I,J),THICACC_MO(I,J),J=1,5)
            WRITE(82,8101)IMONTH,IYEAR,(TBARACC_MO(I,J)-TFREZ, &
                          THLQACC_MO(I,J),THICACC_MO(I,J),J=6,10)
            ELSE
            WRITE(82,8102)IMONTH,IYEAR,(TBARACC_MO(I,J)-TFREZ, &
                          THLQACC_MO(I,J),THICACC_MO(I,J),J=1,3)
            ENDIF   

          call resetclassmon(nltest)
          
 826      CONTINUE   

824     CONTINUE ! I
               
       END IF ! IF(IDAY.EQ.monthend(NT+1).AND.NCOUNT.EQ.NDAY)
      END DO ! NMON

8100  FORMAT(1X,I4,I5,5(F8.2,1X),F8.3,F12.4,3(E12.3,1X),2(A6,I2))
8101  FORMAT(1X,I4,I5,5(F7.2,1X,2F6.3,1X),2(A6,I2))
8102  FORMAT(1X,I4,I5,3(F8.2,1X,2F6.3,1X),2(A6,I2))


end subroutine class_monthly_aw

!==============================================================================================================

subroutine ctem_daily_aw(nltest,nmtest,iday,FAREROT,iyear,jdstd,jdsty,jdendd,jdendy,grclarea)
                           
use ctem_statevars,     only : ctem_tile, vrot, c_switch, &
                               resetctem_g, ctem_grd
use ctem_params, only : icc,ignd,nmos,iccp1

implicit none

! arguments
integer, intent(in) :: nltest
integer, intent(in) :: nmtest
integer, intent(in) :: iday
real, intent(in), dimension(:,:) :: FAREROT
integer, intent(in) :: iyear
integer, intent(in) :: jdstd
integer, intent(in) :: jdsty
integer, intent(in) :: jdendd
integer, intent(in) :: jdendy
real, intent(in), dimension(:) :: grclarea

! pointers

logical, pointer :: dofire
logical, pointer :: lnduseon
logical, pointer :: compete
logical, pointer :: dowetlands
logical, pointer :: obswetf
logical, pointer :: mosaic
real, pointer, dimension(:,:,:) :: fcancmxrow
real, pointer, dimension(:,:,:) :: gppvegrow
real, pointer, dimension(:,:,:) :: nepvegrow
real, pointer, dimension(:,:,:) :: nbpvegrow
real, pointer, dimension(:,:,:) :: nppvegrow
real, pointer, dimension(:,:,:) :: hetroresvegrow
real, pointer, dimension(:,:,:) :: autoresvegrow
real, pointer, dimension(:,:,:) :: litresvegrow
real, pointer, dimension(:,:,:) :: soilcresvegrow
real, pointer, dimension(:,:,:) :: rmlvegaccrow
real, pointer, dimension(:,:,:) :: rmsvegrow
real, pointer, dimension(:,:,:) :: rmrvegrow
real, pointer, dimension(:,:,:) :: rgvegrow
real, pointer, dimension(:,:,:) :: ailcgrow
real, pointer, dimension(:,:,:) :: emit_co2row
real, pointer, dimension(:,:,:) :: emit_corow
real, pointer, dimension(:,:,:) :: emit_ch4row
real, pointer, dimension(:,:,:) :: emit_nmhcrow
real, pointer, dimension(:,:,:) :: emit_h2row
real, pointer, dimension(:,:,:) :: emit_noxrow
real, pointer, dimension(:,:,:) :: emit_n2orow
real, pointer, dimension(:,:,:) :: emit_pm25row
real, pointer, dimension(:,:,:) :: emit_tpmrow
real, pointer, dimension(:,:,:) :: emit_tcrow
real, pointer, dimension(:,:,:) :: emit_ocrow
real, pointer, dimension(:,:,:) :: emit_bcrow
real, pointer, dimension(:,:) :: burnfracrow
real, pointer, dimension(:,:,:) :: burnvegfrow
real, pointer, dimension(:,:) :: probfirerow
real, pointer, dimension(:,:) :: btermrow
real, pointer, dimension(:,:) :: ltermrow
real, pointer, dimension(:,:) :: mtermrow
real, pointer, dimension(:,:) :: lucemcomrow
real, pointer, dimension(:,:) :: lucltrinrow
real, pointer, dimension(:,:) :: lucsocinrow
real, pointer, dimension(:,:) :: ch4wet1row
real, pointer, dimension(:,:) :: ch4wet2row
real, pointer, dimension(:,:) :: wetfdynrow
real, pointer, dimension(:,:) :: ch4dyn1row
real, pointer, dimension(:,:) :: ch4dyn2row
real, pointer, dimension(:,:,:) :: litrmassrow
real, pointer, dimension(:,:,:) :: soilcmasrow
real, pointer, dimension(:,:,:) :: vgbiomas_vegrow
real, pointer, dimension(:,:,:) :: stemmassrow        
real, pointer, dimension(:,:,:) :: rootmassrow        
real, pointer, dimension(:,:,:) :: gleafmasrow        !
real, pointer, dimension(:,:,:) :: bleafmasrow        !
real, pointer, dimension(:,:) :: gavglairow
real, pointer, dimension(:,:,:) :: slairow
real, pointer, dimension(:,:,:) :: ailcbrow
real, pointer, dimension(:,:,:) :: flhrlossrow
real, pointer, dimension(:,:) :: dstcemls3row
integer, pointer, dimension(:,:,:) :: lfstatusrow
real, pointer, dimension(:,:) :: vgbiomasrow
real, pointer, dimension(:,:) :: gavgltmsrow
real, pointer, dimension(:,:) :: gavgscmsrow
      
real, pointer, dimension(:,:) :: leaflitr_m
real, pointer, dimension(:,:) :: tltrleaf_m
real, pointer, dimension(:,:) :: tltrstem_m
real, pointer, dimension(:,:) :: tltrroot_m
real, pointer, dimension(:,:) :: ailcg_m
real, pointer, dimension(:,:) :: ailcb_m
real, pointer, dimension(:,:,:) :: rmatctem_m
real, pointer, dimension(:,:) :: veghght_m
real, pointer, dimension(:,:) :: rootdpth_m
real, pointer, dimension(:,:) :: roottemp_m
real, pointer, dimension(:,:) :: slai_m      
real, pointer, dimension(:,:) :: afrroot_m
real, pointer, dimension(:,:) :: afrleaf_m
real, pointer, dimension(:,:) :: afrstem_m
real, pointer, dimension(:,:) :: laimaxg_m
real, pointer, dimension(:,:) :: stemmass_m
real, pointer, dimension(:,:) :: rootmass_m
real, pointer, dimension(:,:) :: litrmass_m
real, pointer, dimension(:,:) :: gleafmas_m
real, pointer, dimension(:,:) :: bleafmas_m
real, pointer, dimension(:,:) :: soilcmas_m
integer, pointer, dimension(:,:) :: ifcancmx_m
real, pointer, dimension(:,:) :: tcanoaccrow_out

real, pointer, dimension(:,:) :: npprow
real, pointer, dimension(:,:) :: neprow
real, pointer, dimension(:,:) :: nbprow
real, pointer, dimension(:,:) :: gpprow
real, pointer, dimension(:,:) :: hetroresrow
real, pointer, dimension(:,:) :: autoresrow
real, pointer, dimension(:,:) :: soilcresprow
real, pointer, dimension(:,:) :: rgrow
real, pointer, dimension(:,:) :: litresrow
real, pointer, dimension(:,:) :: socresrow

real, pointer, dimension(:) :: gpp_g
real, pointer, dimension(:) :: npp_g
real, pointer, dimension(:) :: nbp_g
real, pointer, dimension(:) :: autores_g
real, pointer, dimension(:) :: socres_g
real, pointer, dimension(:) :: litres_g
real, pointer, dimension(:) :: dstcemls3_g
real, pointer, dimension(:) :: litrfall_g
real, pointer, dimension(:) :: rml_g
real, pointer, dimension(:) :: rms_g      
real, pointer, dimension(:) :: rg_g
real, pointer, dimension(:) :: leaflitr_g
real, pointer, dimension(:) :: tltrstem_g
real, pointer, dimension(:) :: tltrroot_g
real, pointer, dimension(:) :: nep_g
real, pointer, dimension(:) :: hetrores_g
real, pointer, dimension(:) :: dstcemls_g
real, pointer, dimension(:) :: humiftrs_g
real, pointer, dimension(:) :: rmr_g
real, pointer, dimension(:) :: tltrleaf_g
real, pointer, dimension(:) :: gavgltms_g
real, pointer, dimension(:) :: vgbiomas_g
real, pointer, dimension(:) :: gavglai_g
real, pointer, dimension(:) :: gavgscms_g
real, pointer, dimension(:) :: gleafmas_g
real, pointer, dimension(:) :: bleafmas_g
real, pointer, dimension(:) :: stemmass_g
real, pointer, dimension(:) :: rootmass_g
real, pointer, dimension(:) :: litrmass_g
real, pointer, dimension(:) :: soilcmas_g
real, pointer, dimension(:) :: slai_g
real, pointer, dimension(:) :: ailcg_g
real, pointer, dimension(:) :: ailcb_g
real, pointer, dimension(:) :: veghght_g
real, pointer, dimension(:) :: rootdpth_g
real, pointer, dimension(:) :: roottemp_g
real, pointer, dimension(:) :: totcmass_g
real, pointer, dimension(:) :: tcanoacc_out_g
real, pointer, dimension(:) :: burnfrac_g
real, pointer, dimension(:) :: probfire_g
real, pointer, dimension(:) :: lucemcom_g
real, pointer, dimension(:) :: lucltrin_g
real, pointer, dimension(:) :: lucsocin_g
real, pointer, dimension(:) :: emit_co2_g
real, pointer, dimension(:) :: emit_co_g
real, pointer, dimension(:) :: emit_ch4_g
real, pointer, dimension(:) :: emit_nmhc_g
real, pointer, dimension(:) :: emit_h2_g
real, pointer, dimension(:) :: emit_nox_g
real, pointer, dimension(:) :: emit_n2o_g
real, pointer, dimension(:) :: emit_pm25_g
real, pointer, dimension(:) :: emit_tpm_g
real, pointer, dimension(:) :: emit_tc_g
real, pointer, dimension(:) :: emit_oc_g
real, pointer, dimension(:) :: emit_bc_g
real, pointer, dimension(:) :: bterm_g
real, pointer, dimension(:) :: lterm_g
real, pointer, dimension(:) :: mterm_g
real, pointer, dimension(:) :: ch4wet1_g
real, pointer, dimension(:) :: ch4wet2_g
real, pointer, dimension(:) :: wetfdyn_g
real, pointer, dimension(:) :: ch4dyn1_g
real, pointer, dimension(:) :: ch4dyn2_g
real, pointer, dimension(:,:) :: afrleaf_g
real, pointer, dimension(:,:) :: afrstem_g     
real, pointer, dimension(:,:) :: afrroot_g
real, pointer, dimension(:,:) :: lfstatus_g
real, pointer, dimension(:,:) :: rmlvegrow_g
real, pointer, dimension(:,:) :: anvegrow_g
real, pointer, dimension(:,:) :: rmatctem_g

real, pointer, dimension(:,:,:) :: bmasvegrow
real, pointer, dimension(:,:,:) :: cmasvegcrow
real, pointer, dimension(:,:,:) :: veghghtrow
real, pointer, dimension(:,:,:) :: rootdpthrow
real, pointer, dimension(:,:) :: rmlrow
real, pointer, dimension(:,:) :: rmsrow
real, pointer, dimension(:,:,:) :: tltrleafrow
real, pointer, dimension(:,:,:) :: tltrstemrow
real, pointer, dimension(:,:,:) :: tltrrootrow
real, pointer, dimension(:,:,:) :: leaflitrrow
real, pointer, dimension(:,:,:) :: roottemprow
real, pointer, dimension(:,:,:) :: afrleafrow
real, pointer, dimension(:,:,:) :: afrstemrow
real, pointer, dimension(:,:,:) :: afrrootrow
real, pointer, dimension(:,:,:) :: wtstatusrow
real, pointer, dimension(:,:,:) :: ltstatusrow
real, pointer, dimension(:,:) :: rmrrow

real, pointer, dimension(:,:,:,:) :: rmatctemrow
real, pointer, dimension(:,:) :: dstcemlsrow
real, pointer, dimension(:,:) :: litrfallrow
real, pointer, dimension(:,:) :: humiftrsrow
   
! local
integer :: i,m,j,nt,k
real :: barefrac
real :: sumfare

! point pointers

dofire                => c_switch%dofire
lnduseon              => c_switch%lnduseon
compete               => c_switch%compete
dowetlands            => c_switch%dowetlands
obswetf               => c_switch%obswetf
mosaic                => c_switch%mosaic
fcancmxrow            => vrot%fcancmx
gppvegrow         => vrot%gppveg
nepvegrow         => vrot%nepveg
nbpvegrow         => vrot%nbpveg
nppvegrow         => vrot%nppveg
hetroresvegrow    => vrot%hetroresveg
autoresvegrow     => vrot%autoresveg   
litresvegrow      => vrot%litresveg
soilcresvegrow    => vrot%soilcresveg
rmlvegaccrow      => vrot%rmlvegacc
rmsvegrow         => vrot%rmsveg
rmrvegrow         => vrot%rmrveg
rgvegrow          => vrot%rgveg
ailcgrow          => vrot%ailcg
emit_co2row       => vrot%emit_co2
emit_corow        => vrot%emit_co
emit_ch4row       => vrot%emit_ch4
emit_nmhcrow      => vrot%emit_nmhc
emit_h2row        => vrot%emit_h2
emit_noxrow       => vrot%emit_nox
emit_n2orow       => vrot%emit_n2o
emit_pm25row      => vrot%emit_pm25
emit_tpmrow       => vrot%emit_tpm
emit_tcrow        => vrot%emit_tc
emit_ocrow        => vrot%emit_oc
emit_bcrow        => vrot%emit_bc
burnfracrow       => vrot%burnfrac
burnvegfrow       => vrot%burnvegf
probfirerow       => vrot%probfire
btermrow          => vrot%bterm
ltermrow          => vrot%lterm
mtermrow          => vrot%mterm
lucemcomrow       => vrot%lucemcom
lucltrinrow       => vrot%lucltrin
lucsocinrow       => vrot%lucsocin
ch4wet1row        => vrot%ch4wet1
ch4wet2row        => vrot%ch4wet2
wetfdynrow        => vrot%wetfdyn
ch4dyn1row        => vrot%ch4dyn1
ch4dyn2row        => vrot%ch4dyn2
litrmassrow       => vrot%litrmass
soilcmasrow       => vrot%soilcmas
vgbiomas_vegrow   => vrot%vgbiomas_veg
stemmassrow       => vrot%stemmass
rootmassrow       => vrot%rootmass
flhrlossrow       => vrot%flhrloss
dstcemls3row      => vrot%dstcemls3
lfstatusrow       => vrot%lfstatus

leaflitr_m        => ctem_tile%leaflitr_m
tltrleaf_m        => ctem_tile%tltrleaf_m
tltrstem_m        => ctem_tile%tltrstem_m
tltrroot_m        => ctem_tile%tltrroot_m
ailcg_m           => ctem_tile%ailcg_m
ailcb_m           => ctem_tile%ailcb_m
rmatctem_m        => ctem_tile%rmatctem_m
veghght_m         => ctem_tile%veghght_m
rootdpth_m        => ctem_tile%rootdpth_m
roottemp_m        => ctem_tile%roottemp_m
slai_m            => ctem_tile%slai_m
afrroot_m         => ctem_tile%afrroot_m
afrleaf_m         => ctem_tile%afrleaf_m
afrstem_m         => ctem_tile%afrstem_m
laimaxg_m         => ctem_tile%laimaxg_m
stemmass_m        => ctem_tile%stemmass_m    
rootmass_m        => ctem_tile%rootmass_m
litrmass_m        => ctem_tile%litrmass_m
gleafmas_m        => ctem_tile%gleafmas_m       
bleafmas_m        => ctem_tile%bleafmas_m
soilcmas_m        => ctem_tile%soilcmas_m      
ifcancmx_m        => ctem_tile%ifcancmx_m
 
tcanoaccrow_out   => vrot%tcanoaccrow_out

npprow            => vrot%npp
neprow            => vrot%nep
nbprow            => vrot%nbp
gpprow            => vrot%gpp
hetroresrow       => vrot%hetrores
autoresrow        => vrot%autores
soilcresprow      => vrot%soilcresp
rgrow             => vrot%rg
litresrow         => vrot%litres
socresrow         => vrot%socres
vgbiomasrow       => vrot%vgbiomas
gavgltmsrow       => vrot%gavgltms
gavgscmsrow       => vrot%gavgscms
      
gpp_g             => ctem_grd%gpp_g
npp_g             => ctem_grd%npp_g
nbp_g             => ctem_grd%nbp_g
autores_g         => ctem_grd%autores_g
socres_g          => ctem_grd%socres_g
litres_g          => ctem_grd%litres_g
dstcemls3_g       => ctem_grd%dstcemls3_g 
litrfall_g        => ctem_grd%litrfall_g
rml_g             => ctem_grd%rml_g
rms_g             => ctem_grd%rms_g    
rg_g              => ctem_grd%rg_g    
leaflitr_g        => ctem_grd%leaflitr_g
tltrstem_g        => ctem_grd%tltrstem_g
tltrroot_g        => ctem_grd%tltrroot_g
nep_g             => ctem_grd%nep_g
hetrores_g        => ctem_grd%hetrores_g
dstcemls_g        => ctem_grd%dstcemls_g
humiftrs_g        => ctem_grd%humiftrs_g
rmr_g             => ctem_grd%rmr_g
tltrleaf_g        => ctem_grd%tltrleaf_g
gavgltms_g        => ctem_grd%gavgltms_g
vgbiomas_g        => ctem_grd%vgbiomas_g
gavglai_g         => ctem_grd%gavglai_g
gavgscms_g        => ctem_grd%gavgscms_g
gleafmas_g        => ctem_grd%gleafmas_g
bleafmas_g        => ctem_grd%bleafmas_g
stemmass_g        => ctem_grd%stemmass_g
rootmass_g        => ctem_grd%rootmass_g
litrmass_g        => ctem_grd%litrmass_g
soilcmas_g        => ctem_grd%soilcmas_g
slai_g            => ctem_grd%slai_g
ailcg_g           => ctem_grd%ailcg_g
ailcb_g           => ctem_grd%ailcb_g
veghght_g         => ctem_grd%veghght_g
rootdpth_g        => ctem_grd%rootdpth_g  
roottemp_g        => ctem_grd%roottemp_g
totcmass_g        => ctem_grd%totcmass_g
tcanoacc_out_g    => ctem_grd%tcanoacc_out_g
burnfrac_g        => ctem_grd%burnfrac_g
probfire_g        => ctem_grd%probfire_g
lucemcom_g        => ctem_grd%lucemcom_g
lucltrin_g        => ctem_grd%lucltrin_g
lucsocin_g        => ctem_grd%lucsocin_g
emit_co2_g        => ctem_grd%emit_co2_g
emit_co_g         => ctem_grd%emit_co_g
emit_ch4_g        => ctem_grd%emit_ch4_g
emit_nmhc_g       => ctem_grd%emit_nmhc_g
emit_h2_g         => ctem_grd%emit_h2_g
emit_nox_g        => ctem_grd%emit_nox_g
emit_n2o_g        => ctem_grd%emit_n2o_g
emit_pm25_g       => ctem_grd%emit_pm25_g
emit_tpm_g        => ctem_grd%emit_tpm_g
emit_tc_g         => ctem_grd%emit_tc_g
emit_oc_g         => ctem_grd%emit_oc_g
emit_bc_g         => ctem_grd%emit_bc_g
bterm_g           => ctem_grd%bterm_g
lterm_g           => ctem_grd%lterm_g
mterm_g           => ctem_grd%mterm_g
ch4wet1_g         => ctem_grd%ch4wet1_g
ch4wet2_g         => ctem_grd%ch4wet2_g
wetfdyn_g         => ctem_grd%wetfdyn_g
ch4dyn1_g         => ctem_grd%ch4dyn1_g
ch4dyn2_g         => ctem_grd%ch4dyn2_g
afrleaf_g         => ctem_grd%afrleaf_g
afrstem_g         => ctem_grd%afrstem_g   
afrroot_g         => ctem_grd%afrroot_g
lfstatus_g        => ctem_grd%lfstatus_g
rmlvegrow_g       => ctem_grd%rmlvegrow_g
anvegrow_g        => ctem_grd%anvegrow_g
rmatctem_g        => ctem_grd%rmatctem_g    

bmasvegrow        => vrot%bmasveg
cmasvegcrow       => vrot%cmasvegc
veghghtrow        => vrot%veghght
rootdpthrow       => vrot%rootdpth
rmlrow            => vrot%rml
rmsrow            => vrot%rms
tltrleafrow       => vrot%tltrleaf
tltrstemrow       => vrot%tltrstem
tltrrootrow       => vrot%tltrroot
leaflitrrow       => vrot%leaflitr
roottemprow       => vrot%roottemp
afrleafrow        => vrot%afrleaf
afrstemrow        => vrot%afrstem
afrrootrow        => vrot%afrroot    
wtstatusrow       => vrot%wtstatus
ltstatusrow       => vrot%ltstatus
rmrrow            => vrot%rmr

gleafmasrow       => vrot%gleafmas
bleafmasrow       => vrot%bleafmas
gavglairow        => vrot%gavglai
slairow           => vrot%slai
ailcbrow          => vrot%ailcb

flhrlossrow       => vrot%flhrloss

rmatctemrow       => vrot%rmatctem
dstcemlsrow       => vrot%dstcemls
litrfallrow       => vrot%litrfall
humiftrsrow       => vrot%humiftrs
   
!--------

!do i=1,nltest
!    do j=1,icc
!        ifcancmx_g(i,j)=0   !pointless var? JM.
!    enddo
!enddo

do i=1,nltest
    do m=1,nmtest
        ifcancmx_m(i,m)=0   !0=bare soil tile; 1=tile with vegetation
        leaflitr_m(i,m)=0.0
        tltrleaf_m(i,m)=0.0
        tltrstem_m(i,m)=0.0
        tltrroot_m(i,m)=0.0
        ailcg_m(i,m)=0.0
        ailcb_m(i,m)=0.0
        afrleaf_m(i,m)=0.0
        afrstem_m(i,m)=0.0
        afrroot_m(i,m)=0.0
        veghght_m(i,m)=0.0
        rootdpth_m(i,m)=0.0
        roottemp_m(i,m)=0.0
        slai_m(i,m)=0.0
        gleafmas_m(i,m) = 0.0
        bleafmas_m(i,m) = 0.0
        stemmass_m(i,m) = 0.0
        rootmass_m(i,m) = 0.0
        litrmass_m(i,m) = 0.0
        soilcmas_m(i,m) = 0.0

        do j=1,icc
        if (fcancmxrow(i,m,j) .gt.0.0) then
        !ifcancmx_g(i,j)=1
        ifcancmx_m(i,m)=1
        endif 
        enddo

        do k=1,ignd
        rmatctem_m(i,m,k)=0.0
        enddo

    enddo ! m
enddo  ! i

!       ---------------------------------------------------------

call resetctem_g(nltest)

do 851 i=1,nltest
  do 852 m=1,nmtest
    do j=1,icc
        leaflitr_m(i,m)=leaflitr_m(i,m)+leaflitrrow(i,m,j)*fcancmxrow(i,m,j)
        tltrleaf_m(i,m)=tltrleaf_m(i,m)+tltrleafrow(i,m,j)*fcancmxrow(i,m,j)
        tltrstem_m(i,m)=tltrstem_m(i,m)+tltrstemrow(i,m,j)*fcancmxrow(i,m,j)
        tltrroot_m(i,m)=tltrroot_m(i,m)+tltrrootrow(i,m,j)*fcancmxrow(i,m,j)
        veghght_m(i,m)=veghght_m(i,m)+veghghtrow(i,m,j)*fcancmxrow(i,m,j)
        rootdpth_m(i,m)=rootdpth_m(i,m)+rootdpthrow(i,m,j)*fcancmxrow(i,m,j)
        roottemp_m(i,m)=roottemp_m(i,m)+roottemprow(i,m,j)*fcancmxrow(i,m,j)
        slai_m(i,m)=slai_m(i,m)+slairow(i,m,j)*fcancmxrow(i,m,j)

        afrleaf_m(i,m)=afrleaf_m(i,m)+afrleafrow(i,m,j)*fcancmxrow(i,m,j)
        afrstem_m(i,m)=afrstem_m(i,m)+afrstemrow(i,m,j)*fcancmxrow(i,m,j)
        afrroot_m(i,m)=afrroot_m(i,m)+afrrootrow(i,m,j)*fcancmxrow(i,m,j)

        ailcg_m(i,m)=ailcg_m(i,m)+ailcgrow(i,m,j)*fcancmxrow(i,m,j)
        ailcb_m(i,m)=ailcb_m(i,m)+ailcbrow(i,m,j)*fcancmxrow(i,m,j)

        gleafmas_m(i,m) = gleafmas_m(i,m) + gleafmasrow(i,m,j)*fcancmxrow(i,m,j)
        bleafmas_m(i,m) = bleafmas_m(i,m) + bleafmasrow(i,m,j)*fcancmxrow(i,m,j)
        stemmass_m(i,m) = stemmass_m(i,m) + stemmassrow(i,m,j)*fcancmxrow(i,m,j)
        rootmass_m(i,m) = rootmass_m(i,m) + rootmassrow(i,m,j)*fcancmxrow(i,m,j)
        litrmass_m(i,m) = litrmass_m(i,m) + litrmassrow(i,m,j)*fcancmxrow(i,m,j)
        soilcmas_m(i,m) = soilcmas_m(i,m) + soilcmasrow(i,m,j)*fcancmxrow(i,m,j)

        do k=1,ignd
        rmatctem_m(i,m,k)=rmatctem_m(i,m,k)+rmatctemrow(i,m,j,k)*fcancmxrow(i,m,j)
        enddo
    enddo

    npprow(i,m)     =npprow(i,m)*1.0377 ! convert to gc/m2.day
    gpprow(i,m)     =gpprow(i,m)*1.0377 ! convert to gc/m2.day
    neprow(i,m)     =neprow(i,m)*1.0377 ! convert to gc/m2.day
    nbprow(i,m)     =nbprow(i,m)*1.0377 ! convert to gc/m2.day
    lucemcomrow(i,m)=lucemcomrow(i,m)*1.0377 ! convert to gc/m2.day
    lucltrinrow(i,m)=lucltrinrow(i,m)*1.0377 ! convert to gc/m2.day
    lucsocinrow(i,m)=lucsocinrow(i,m)*1.0377 ! convert to gc/m2.day

    hetroresrow(i,m)=hetroresrow(i,m)*1.0377 ! convert to gc/m2.day
    autoresrow(i,m) =autoresrow(i,m)*1.0377  ! convert to gc/m2.day
    litresrow(i,m)  =litresrow(i,m)*1.0377   ! convert to gc/m2.day
    socresrow(i,m)  =socresrow(i,m)*1.0377   ! convert to gc/m2.day

    CH4WET1ROW(i,m) = CH4WET1ROW(i,m)*1.0377 * 16.044 / 12. ! convert from umolCH4/m2/s to gCH4/m2.day 
    CH4WET2ROW(i,m) = CH4WET2ROW(i,m)*1.0377 * 16.044 / 12. ! convert from umolCH4/m2/s to gCH4/m2.day
    CH4DYN1ROW(i,m) = CH4DYN1ROW(i,m)*1.0377 * 16.044 / 12. ! convert from umolCH4/m2/s to gCH4/m2.day
    CH4DYN2ROW(i,m) = CH4DYN2ROW(i,m)*1.0377 * 16.044 / 12. ! convert from umolCH4/m2/s to gCH4/m2.day 

!          write daily ctem results
    if ((iyear .ge. jdsty).and.(iyear.le.jdendy))then
     if ((iday .ge. jdstd).and.(iday .le.jdendd))then

!             write grid-averaged fluxes of basic quantities to 
!             file *.CT01D_M

        if (mosaic) then
        write(72,8200)iday,iyear,gpprow(i,m),npprow(i,m), &
                neprow(i,m),nbprow(i,m),autoresrow(i,m), &
                hetroresrow(i,m),litresrow(i,m),socresrow(i,m), &
                (dstcemlsrow(i,m)+dstcemls3row(i,m)), &
               litrfallrow(i,m),humiftrsrow(i,m),' TILE ',m,'AVGE'
        end if

!             write breakdown of some of basic fluxes to file *.CT3 
!             and selected litter fluxes for selected pft

!       First for the bare fraction of the grid cell.
        hetroresvegrow(i,m,iccp1)=hetroresvegrow(i,m,iccp1)*1.0377 ! convert to gc/m2.day
        litresvegrow(i,m,iccp1)=litresvegrow(i,m,iccp1)*1.0377 ! convert to gc/m2.day
        soilcresvegrow(i,m,iccp1)=soilcresvegrow(i,m,iccp1)*1.0377 ! convert to gc/m2.day

        do 853 j=1,icc
        
        if (fcancmxrow(i,m,j) .gt.0.0) then
        
            gppvegrow(i,m,j)=gppvegrow(i,m,j)*1.0377 ! convert to gc/m2.day
            nppvegrow(i,m,j)=nppvegrow(i,m,j)*1.0377 ! convert to gc/m2.day
            nepvegrow(i,m,j)=nepvegrow(i,m,j)*1.0377 ! convert to gc/m2.day
            nbpvegrow(i,m,j)=nbpvegrow(i,m,j)*1.0377 ! convert to gc/m2.day
            hetroresvegrow(i,m,j)=hetroresvegrow(i,m,j)*1.0377 ! convert to gc/m2.day
            autoresvegrow(i,m,j)=autoresvegrow(i,m,j)*1.0377 ! convert to gc/m2.day
            litresvegrow(i,m,j)=litresvegrow(i,m,j)*1.0377 ! convert to gc/m2.day
            soilcresvegrow(i,m,j)=soilcresvegrow(i,m,j)*1.0377 ! convert to gc/m2.day

!                write to file .CT01D_M 
            if (mosaic) then
            write(72,8201)iday,iyear,gppvegrow(i,m,j), &
            nppvegrow(i,m,j),nepvegrow(i,m,j), &
            ' TILE ',m,'PFT',j
            

!                write to file .CT02D_M 
            write(73,8300)iday,iyear,rmlvegaccrow(i,m,j), & 
           rmsvegrow(i,m,j),rmrvegrow(i,m,j),rgvegrow(i,m,j), &
           leaflitrrow(i,m,j),tltrleafrow(i,m,j), &
          tltrstemrow(i,m,j),tltrrootrow(i,m,j),' TILE ',m,'PFT',j


!                write grid-averaged pool sizes and component sizes for
!                seleced ctem pft to file *.CT03D_M

            write(74,8401)iday,iyear,vgbiomas_vegrow(i,m,j), &
               ailcgrow(i,m,j),gleafmasrow(i,m,j), &
               bleafmasrow(i,m,j), stemmassrow(i,m,j), &
               rootmassrow(i,m,j), litrmassrow(i,m,j),  &
               soilcmasrow(i,m,j),' TILE ',m,'PFT',j

!                write lai, rmatctem, & structural attributes for selected 
!                pft to file *.CT04D_M

            write(75,8500)iday,iyear, ailcgrow(i,m,j),  &
                ailcbrow(i,m,j),(rmatctemrow(i,m,j,k),k=1,3), &
                veghghtrow(i,m,j),rootdpthrow(i,m,j), &
              roottemprow(i,m,j),slairow(i,m,j),' TILE ',m,'PFT',j

!                write allocation fractions for selected pft to 
!                file *.CT05D_M

            write(76,8600)iday,iyear, afrleafrow(i,m,j),  &
                afrstemrow(i,m,j),afrrootrow(i,m,j),  &
                tcanoaccrow_out(i,m), lfstatusrow(i,m,j), &
                ' TILE ',m,'PFT',j

!             write fire and luc results to file *.CT06D_M

        if (dofire .or. lnduseon) then   !FLAG FIX THIS
        write(78,8800)iday,iyear, &
         emit_co2row(i,m,j),emit_corow(i,m,j),emit_ch4row(i,m,j), &
         emit_nmhcrow(i,m,j),emit_h2row(i,m,j),emit_noxrow(i,m,j), &
         emit_n2orow(i,m,j),emit_pm25row(i,m,j), &
         emit_tpmrow(i,m,j),emit_tcrow(i,m,j),emit_ocrow(i,m,j), &
         emit_bcrow(i,m,j),burnvegfrow(i,m,j),probfirerow(i,m), & 
         lucemcomrow(i,m),lucltrinrow(i,m), lucsocinrow(i,m), &
         grclarea(i),btermrow(i,m),ltermrow(i,m),mtermrow(i,m), &
         ' TILE ',m,'PFT',j
        endif

        end if !mosaic

        endif  !if (fcancmxrow(i,m,j) .gt.0.0) then

853           continue

        if (mosaic) then
        if (ifcancmx_m(i,m) .gt. 0) then


!               write to file .CT02D_M 
        write(73,8300)iday,iyear,rmlrow(i,m),rmsrow(i,m), &
          rmrrow(i,m),rgrow(i,m),leaflitr_m(i,m),tltrleaf_m(i,m), &
          tltrstem_m(i,m),tltrroot_m(i,m),' TILE ',m,'AVGE'

!               write to file .CT03D_M 
        write(74,8402)iday,iyear,vgbiomasrow(i,m), &
               gavglairow(i,m),gavgltmsrow(i,m), &
               gavgscmsrow(i,m),gleafmasrow(i,m,j), &
               bleafmasrow(i,m,j), stemmassrow(i,m,j), &
               rootmassrow(i,m,j), litrmassrow(i,m,j), & 
               soilcmasrow(i,m,j),' TILE ',m, 'AVGE'

!               write to file .CT04D_M
        write(75,8500)iday,iyear,ailcg_m(i,m), &
                ailcb_m(i,m),(rmatctem_m(i,m,k),k=1,3), &
                veghght_m(i,m),rootdpth_m(i,m), &
                roottemp_m(i,m),slai_m(i,m),' TILE ',m, 'AVGE'

!               write to file .CT05D_M
        write(76,8601)iday,iyear, afrleaf_m(i,m), &
                afrstem_m(i,m),afrroot_m(i,m),  &
                tcanoaccrow_out(i,m),  &
                ' TILE ',m,'AVGE'

        end if !if (ifcancmx_m(i,m) .gt.0.0) then
        endif !mosaic
!
    end if
    endif ! if write daily

8200       format(1x,i4,i5,11f10.5,2(a6,i2))
8201       format(1x,i4,i5,3f10.5,80x,2(a6,i2))
8300       format(1x,i4,i5,8f10.5,2(a6,i2))
8301       format(1x,i4,i5,4f10.5,40x,2(a6,i2)) 
8400       format(1x,i4,i5,11f10.5,2(a6,i2))
8401       format(1x,i4,i5,2f10.5,6f10.5,2(a6,i2))
8402       format(1x,i4,i5,10f10.5,2(a6,i2))
!                   8402       format(1x,i4,i5,2f10.5,40x,2f10.5,2(a6,i2))
8500       format(1x,i4,i5,9f10.5,2(a6,i2))
8600       format(1x,i4,i5,4f10.5,i8,2(a6,i2))
8601       format(1x,i4,i5,4f10.5,8x,2(a6,i2))   
8800       format(1x,i4,i5,20f11.4,2x,f9.2,2(a6,i2))
8810       format(1x,i4,i5,5f11.4,2(a6,i2))

!          Calculation of grid averaged variables

    gpp_g(i) =gpp_g(i) + gpprow(i,m)*FAREROT(i,m)
    npp_g(i) =npp_g(i) + npprow(i,m)*FAREROT(i,m)
    nep_g(i) =nep_g(i) + neprow(i,m)*FAREROT(i,m)
    nbp_g(i) =nbp_g(i) + nbprow(i,m)*FAREROT(i,m)
    autores_g(i) =autores_g(i) +autoresrow(i,m)*FAREROT(i,m)
    hetrores_g(i)=hetrores_g(i)+hetroresrow(i,m)*FAREROT(i,m)
    litres_g(i) =litres_g(i) + litresrow(i,m)*FAREROT(i,m)
    socres_g(i) =socres_g(i) + socresrow(i,m)*FAREROT(i,m)
    dstcemls_g(i)=dstcemls_g(i)+dstcemlsrow(i,m)*FAREROT(i,m)
    dstcemls3_g(i)=dstcemls3_g(i)+dstcemls3row(i,m)*FAREROT(i,m)

    litrfall_g(i)=litrfall_g(i)+litrfallrow(i,m)*FAREROT(i,m)
    humiftrs_g(i)=humiftrs_g(i)+humiftrsrow(i,m)*FAREROT(i,m)
    rml_g(i) =rml_g(i) + rmlrow(i,m)*FAREROT(i,m)
    rms_g(i) =rms_g(i) + rmsrow(i,m)*FAREROT(i,m)
    rmr_g(i) =rmr_g(i) + rmrrow(i,m)*FAREROT(i,m)
    rg_g(i) =rg_g(i) + rgrow(i,m)*FAREROT(i,m)
    leaflitr_g(i) = leaflitr_g(i) + leaflitr_m(i,m)*FAREROT(i,m)
    tltrleaf_g(i) = tltrleaf_g(i) + tltrleaf_m(i,m)*FAREROT(i,m)
    tltrstem_g(i) = tltrstem_g(i) + tltrstem_m(i,m)*FAREROT(i,m)
    tltrroot_g(i) = tltrroot_g(i) + tltrroot_m(i,m)*FAREROT(i,m)
    vgbiomas_g(i) =vgbiomas_g(i) + vgbiomasrow(i,m)*FAREROT(i,m)
    gavglai_g(i) =gavglai_g(i) + gavglairow(i,m)*FAREROT(i,m)
    gavgltms_g(i) =gavgltms_g(i) + gavgltmsrow(i,m)*FAREROT(i,m)
    gavgscms_g(i) =gavgscms_g(i) + gavgscmsrow(i,m)*FAREROT(i,m)
    tcanoacc_out_g(i) =tcanoacc_out_g(i)+tcanoaccrow_out(i,m)*FAREROT(i,m)
    totcmass_g(i) =vgbiomas_g(i) + gavgltms_g(i) + gavgscms_g(i)
    gleafmas_g(i) = gleafmas_g(i) + gleafmas_m(i,m)*FAREROT(i,m)
    bleafmas_g(i) = bleafmas_g(i) + bleafmas_m(i,m)*FAREROT(i,m)
    stemmass_g(i) = stemmass_g(i) + stemmass_m(i,m)*FAREROT(i,m)
    rootmass_g(i) = rootmass_g(i) + rootmass_m(i,m)*FAREROT(i,m)
    litrmass_g(i) = litrmass_g(i) + litrmass_m(i,m)*FAREROT(i,m)
    soilcmas_g(i) = soilcmas_g(i) + soilcmas_m(i,m)*FAREROT(i,m)

    burnfrac_g(i) =burnfrac_g(i)+ burnfracrow(i,m)*FAREROT(i,m) 

    probfire_g(i) =probfire_g(i)+probfirerow(i,m)*FAREROT(i,m)
    lucemcom_g(i) =lucemcom_g(i)+lucemcomrow(i,m)*FAREROT(i,m)
    lucltrin_g(i) =lucltrin_g(i)+lucltrinrow(i,m)*FAREROT(i,m)
    lucsocin_g(i) =lucsocin_g(i)+lucsocinrow(i,m)*FAREROT(i,m) 
    bterm_g(i)    =bterm_g(i)   +btermrow(i,m)*FAREROT(i,m) 
    lterm_g(i)    =lterm_g(i)   +ltermrow(i,m)*FAREROT(i,m) 
    mterm_g(i)    =mterm_g(i)   +mtermrow(i,m)*FAREROT(i,m)

    CH4WET1_G(i) = CH4WET1_G(i) + CH4WET1ROW(i,m)*FAREROT(i,m)
    CH4WET2_G(i) = CH4WET2_G(i) + CH4WET2ROW(i,m)*FAREROT(i,m)
    WETFDYN_G(i) = WETFDYN_G(i) + WETFDYNROW(i,m)*FAREROT(i,m)
    CH4DYN1_G(i) = CH4DYN1_G(i) + CH4DYN1ROW(i,m)*FAREROT(i,m)
    CH4DYN2_G(i) = CH4DYN2_G(i) + CH4DYN2ROW(i,m)*FAREROT(i,m)

    do j=1,icc  

    do k=1,ignd
        rmatctem_g(i,k)=rmatctem_g(i,k)+rmatctemrow(i,m,j,k)*FAREROT(i,m)*fcancmxrow(i,m,j)
    end do

    veghght_g(i) = veghght_g(i) + veghghtrow(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
    rootdpth_g(i) = rootdpth_g(i) + rootdpthrow(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
    roottemp_g(i) = roottemp_g(i) + roottemprow(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
    slai_g(i) = slai_g(i) + slairow(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)

    emit_co2_g(i) =emit_co2_g(i)+ emit_co2row(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
    emit_co_g(i)  =emit_co_g(i) + emit_corow(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
    emit_ch4_g(i) =emit_ch4_g(i)+ emit_ch4row(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
    emit_nmhc_g(i)=emit_nmhc_g(i)+emit_nmhcrow(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
    emit_h2_g(i)  =emit_h2_g(i) + emit_h2row(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
    emit_nox_g(i) =emit_nox_g(i)+ emit_noxrow(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
    emit_n2o_g(i) =emit_n2o_g(i)+ emit_n2orow(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
    emit_pm25_g(i)=emit_pm25_g(i)+emit_pm25row(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
    emit_tpm_g(i) =emit_tpm_g(i)+ emit_tpmrow(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
    emit_tc_g(i)  =emit_tc_g(i) + emit_tcrow(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
    emit_oc_g(i)  =emit_oc_g(i) + emit_ocrow(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
    emit_bc_g(i)  =emit_bc_g(i) + emit_bcrow(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)

    ailcg_g(i)=ailcg_g(i)+ailcgrow(i,m,j)*fcancmxrow(i,m,j)*FAREROT(i,m)
    ailcb_g(i)=ailcb_g(i)+ailcbrow(i,m,j)*fcancmxrow(i,m,j)*FAREROT(i,m)

    enddo 

852       continue

    if ((iyear .ge. jdsty).and.(iyear.le.jdendy))then
     if ((iday .ge. jdstd).and.(iday .le.jdendd))then

!           write to file .CT01D_G               
    write(721,8200)iday,iyear,gpp_g(i),npp_g(i), &
                nep_g(i),nbp_g(i),autores_g(i), &
                hetrores_g(i),litres_g(i),socres_g(i), &
                (dstcemls_g(i)+dstcemls3_g(i)), &
                litrfall_g(i),humiftrs_g(i)

!           write breakdown of some of basic fluxes to file  
!           *.CT02D_G and selected litter fluxes for selected pft
    write(731,8300)iday,iyear,rml_g(i),rms_g(i), &
          rmr_g(i),rg_g(i),leaflitr_g(i),tltrleaf_g(i), &
          tltrstem_g(i),tltrroot_g(i)

!           write to file .CT03D_G
    write(741,8400)iday,iyear,vgbiomas_g(i), &
               gavglai_g(i),gavgltms_g(i), &
               gavgscms_g(i), totcmass_g(i), &
               gleafmas_g(i), bleafmas_g(i), stemmass_g(i), &
               rootmass_g(i), litrmass_g(i), soilcmas_g(i)

!           write to file .CT04D_G
    write(751,8500)iday,iyear, ailcg_g(i),  &
                ailcb_g(i),(rmatctem_g(i,k),k=1,3), &
                veghght_g(i),rootdpth_g(i),roottemp_g(i),slai_g(i)

!           write fire and luc results to file *.CT06D_G

    if (dofire .or. lnduseon) then
    write(781,8800)iday,iyear,  &
          emit_co2_g(i), emit_co_g(i), emit_ch4_g(i), &
          emit_nmhc_g(i), emit_h2_g(i), emit_nox_g(i), &
          emit_n2o_g(i), emit_pm25_g(i), emit_tpm_g(i), &
          emit_tc_g(i), emit_oc_g(i), emit_bc_g(i), &
          burnfrac_g(i)*100., probfire_g(i),lucemcom_g(i), & 
          lucltrin_g(i), lucsocin_g(i), &
          grclarea(i), bterm_g(i), lterm_g(i), mterm_g(i)
    endif
      
!      write CH4 variables to file *.CT08D_G

    if (dowetlands .or. obswetf) then
    write(762,8810)iday,iyear, ch4wet1_g(i),  &
                 ch4wet2_g(i), wetfdyn_g(i),  &
                 ch4dyn1_g(i), ch4dyn2_g(i)
    endif  

    if (compete .or. lnduseon) then 
        sumfare=0.0
        if (mosaic) then
        do m=1,nmos
            sumfare=sumfare+FAREROT(i,m)
        enddo
        write(761,8200)iday,iyear,(FAREROT(i,m)*100.,m=1,nmos),sumfare
        else !composite
        do j=1,icc  !m = 1
            sumfare=sumfare+fcancmxrow(i,1,j)
        enddo
        write(761,8200)iday,iyear,(fcancmxrow(i,1,j)*100.,j=1,icc),(1.0-sumfare)*100.,sumfare
        endif !mosaic/composite
    endif !compete/lnduseon

    end if
    endif !if write daily

851     continue

end subroutine ctem_daily_aw

!==============================================================================================================

subroutine ctem_monthly_aw(nltest,nmtest,iday,FAREROT,iyear,nday)

use ctem_statevars,     only : ctem_tile_mo, vrot, ctem_grd_mo, c_switch
use ctem_params, only : icc,iccp1,nmon,mmday,monthend,monthdays,seed,nmos

implicit none

! arguments
integer, intent(in) :: nltest
integer, intent(in) :: nmtest
integer, intent(in) :: iday
real, intent(in), dimension(:,:) :: FAREROT
integer, intent(in) :: iyear
integer, intent(in) :: nday

! pointers

logical, pointer :: dofire
logical, pointer :: lnduseon
logical, pointer :: compete
logical, pointer :: dowetlands
logical, pointer :: obswetf
logical, pointer :: mosaic
real, pointer, dimension(:,:,:) :: fcancmxrow
real, pointer, dimension(:,:,:) :: laimaxg_mo_m
real, pointer, dimension(:,:,:) :: stemmass_mo_m
real, pointer, dimension(:,:,:) :: rootmass_mo_m
real, pointer, dimension(:,:,:) :: npp_mo_m
real, pointer, dimension(:,:,:) :: gpp_mo_m
real, pointer, dimension(:,:,:) :: vgbiomas_mo_m
real, pointer, dimension(:,:,:) :: autores_mo_m
real, pointer, dimension(:,:,:) :: totcmass_mo_m
real, pointer, dimension(:,:,:) :: litrmass_mo_m
real, pointer, dimension(:,:,:) :: soilcmas_mo_m
real, pointer, dimension(:,:,:) :: nep_mo_m
real, pointer, dimension(:,:,:) :: litres_mo_m
real, pointer, dimension(:,:,:) :: soilcres_mo_m
real, pointer, dimension(:,:,:) :: hetrores_mo_m
real, pointer, dimension(:,:,:) :: nbp_mo_m
real, pointer, dimension(:,:,:) :: emit_co2_mo_m
real, pointer, dimension(:,:,:) :: emit_co_mo_m
real, pointer, dimension(:,:,:) :: emit_ch4_mo_m
real, pointer, dimension(:,:,:) :: emit_nmhc_mo_m
real, pointer, dimension(:,:,:) :: emit_h2_mo_m
real, pointer, dimension(:,:,:) :: emit_nox_mo_m
real, pointer, dimension(:,:,:) :: emit_n2o_mo_m
real, pointer, dimension(:,:,:) :: emit_pm25_mo_m
real, pointer, dimension(:,:,:) :: emit_tpm_mo_m
real, pointer, dimension(:,:,:) :: emit_tc_mo_m
real, pointer, dimension(:,:,:) :: emit_oc_mo_m
real, pointer, dimension(:,:,:) :: emit_bc_mo_m
real, pointer, dimension(:,:,:) :: burnfrac_mo_m
real, pointer, dimension(:,:) :: probfire_mo_m
real, pointer, dimension(:,:) :: bterm_mo_m
real, pointer, dimension(:,:) :: luc_emc_mo_m
real, pointer, dimension(:,:) :: lterm_mo_m
real, pointer, dimension(:,:) :: lucsocin_mo_m
real, pointer, dimension(:,:) :: mterm_mo_m
real, pointer, dimension(:,:) :: lucltrin_mo_m
real, pointer, dimension(:,:) :: ch4wet1_mo_m
real, pointer, dimension(:,:) :: ch4wet2_mo_m
real, pointer, dimension(:,:) :: wetfdyn_mo_m
real, pointer, dimension(:,:) :: ch4dyn1_mo_m
real, pointer, dimension(:,:) :: ch4dyn2_mo_m
logical, pointer, dimension(:,:,:) :: pftexistrow
real, pointer, dimension(:,:,:) :: gppvegrow
real, pointer, dimension(:,:,:) :: nepvegrow
real, pointer, dimension(:,:,:) :: nbpvegrow
real, pointer, dimension(:,:,:) :: nppvegrow
real, pointer, dimension(:,:,:) :: hetroresvegrow
real, pointer, dimension(:,:,:) :: autoresvegrow
real, pointer, dimension(:,:,:) :: litresvegrow
real, pointer, dimension(:,:,:) :: soilcresvegrow
real, pointer, dimension(:,:,:) :: rmlvegaccrow
real, pointer, dimension(:,:,:) :: rmsvegrow
real, pointer, dimension(:,:,:) :: rmrvegrow
real, pointer, dimension(:,:,:) :: rgvegrow
real, pointer, dimension(:,:,:) :: ailcgrow
real, pointer, dimension(:,:,:) :: emit_co2row
real, pointer, dimension(:,:,:) :: emit_corow
real, pointer, dimension(:,:,:) :: emit_ch4row
real, pointer, dimension(:,:,:) :: emit_nmhcrow
real, pointer, dimension(:,:,:) :: emit_h2row
real, pointer, dimension(:,:,:) :: emit_noxrow
real, pointer, dimension(:,:,:) :: emit_n2orow
real, pointer, dimension(:,:,:) :: emit_pm25row
real, pointer, dimension(:,:,:) :: emit_tpmrow
real, pointer, dimension(:,:,:) :: emit_tcrow
real, pointer, dimension(:,:,:) :: emit_ocrow
real, pointer, dimension(:,:,:) :: emit_bcrow
real, pointer, dimension(:,:) :: burnfracrow
real, pointer, dimension(:,:,:) :: burnvegfrow
real, pointer, dimension(:,:) :: probfirerow
real, pointer, dimension(:,:) :: btermrow
real, pointer, dimension(:,:) :: ltermrow
real, pointer, dimension(:,:) :: mtermrow
real, pointer, dimension(:,:) :: lucemcomrow
real, pointer, dimension(:,:) :: lucltrinrow
real, pointer, dimension(:,:) :: lucsocinrow
real, pointer, dimension(:,:) :: ch4wet1row
real, pointer, dimension(:,:) :: ch4wet2row
real, pointer, dimension(:,:) :: wetfdynrow
real, pointer, dimension(:,:) :: ch4dyn1row
real, pointer, dimension(:,:) :: ch4dyn2row
real, pointer, dimension(:,:,:) :: litrmassrow
real, pointer, dimension(:,:,:) :: soilcmasrow
real, pointer, dimension(:,:,:) :: vgbiomas_vegrow
real, pointer, dimension(:,:,:) :: stemmassrow
real, pointer, dimension(:,:,:) :: rootmassrow

real, pointer, dimension(:) :: laimaxg_mo_g
real, pointer, dimension(:) :: stemmass_mo_g
real, pointer, dimension(:) :: rootmass_mo_g
real, pointer, dimension(:) :: litrmass_mo_g
real, pointer, dimension(:) :: soilcmas_mo_g
real, pointer, dimension(:) :: npp_mo_g
real, pointer, dimension(:) :: gpp_mo_g
real, pointer, dimension(:) :: nep_mo_g
real, pointer, dimension(:) :: nbp_mo_g
real, pointer, dimension(:) :: hetrores_mo_g
real, pointer, dimension(:) :: autores_mo_g
real, pointer, dimension(:) :: litres_mo_g
real, pointer, dimension(:) :: soilcres_mo_g
real, pointer, dimension(:) :: vgbiomas_mo_g
real, pointer, dimension(:) :: totcmass_mo_g
real, pointer, dimension(:) :: emit_co2_mo_g
real, pointer, dimension(:) :: emit_co_mo_g
real, pointer, dimension(:) :: emit_ch4_mo_g
real, pointer, dimension(:) :: emit_nmhc_mo_g
real, pointer, dimension(:) :: emit_h2_mo_g
real, pointer, dimension(:) :: emit_nox_mo_g
real, pointer, dimension(:) :: emit_n2o_mo_g
real, pointer, dimension(:) :: emit_pm25_mo_g
real, pointer, dimension(:) :: emit_tpm_mo_g
real, pointer, dimension(:) :: emit_tc_mo_g
real, pointer, dimension(:) :: emit_oc_mo_g
real, pointer, dimension(:) :: emit_bc_mo_g
real, pointer, dimension(:) :: probfire_mo_g
real, pointer, dimension(:) :: luc_emc_mo_g
real, pointer, dimension(:) :: lucltrin_mo_g
real, pointer, dimension(:) :: lucsocin_mo_g
real, pointer, dimension(:) :: burnfrac_mo_g
real, pointer, dimension(:) :: bterm_mo_g
real, pointer, dimension(:) :: lterm_mo_g
real, pointer, dimension(:) :: mterm_mo_g
real, pointer, dimension(:) :: ch4wet1_mo_g
real, pointer, dimension(:) :: ch4wet2_mo_g
real, pointer, dimension(:) :: wetfdyn_mo_g
real, pointer, dimension(:) :: ch4dyn1_mo_g
real, pointer, dimension(:) :: ch4dyn2_mo_g

! local
integer :: i,m,j,nt
real :: barefrac
real :: sumfare
integer :: NDMONTH
integer :: imonth

! point pointers

dofire                => c_switch%dofire
lnduseon              => c_switch%lnduseon
compete               => c_switch%compete
dowetlands            => c_switch%dowetlands
obswetf               => c_switch%obswetf
mosaic                => c_switch%mosaic
pftexistrow           => vrot%pftexist
fcancmxrow            => vrot%fcancmx
laimaxg_mo_m          =>ctem_tile_mo%laimaxg_mo_m
stemmass_mo_m         =>ctem_tile_mo%stemmass_mo_m
rootmass_mo_m         =>ctem_tile_mo%rootmass_mo_m
npp_mo_m              =>ctem_tile_mo%npp_mo_m
gpp_mo_m              =>ctem_tile_mo%gpp_mo_m
vgbiomas_mo_m         =>ctem_tile_mo%vgbiomas_mo_m
autores_mo_m          =>ctem_tile_mo%autores_mo_m
totcmass_mo_m         =>ctem_tile_mo%totcmass_mo_m
litrmass_mo_m         =>ctem_tile_mo%litrmass_mo_m
soilcmas_mo_m         =>ctem_tile_mo%soilcmas_mo_m
nep_mo_m              =>ctem_tile_mo%nep_mo_m
litres_mo_m           =>ctem_tile_mo%litres_mo_m
soilcres_mo_m         =>ctem_tile_mo%soilcres_mo_m
hetrores_mo_m         =>ctem_tile_mo%hetrores_mo_m
nbp_mo_m              =>ctem_tile_mo%nbp_mo_m
emit_co2_mo_m         =>ctem_tile_mo%emit_co2_mo_m
emit_co_mo_m          =>ctem_tile_mo%emit_co_mo_m
emit_ch4_mo_m         =>ctem_tile_mo%emit_ch4_mo_m
emit_nmhc_mo_m        =>ctem_tile_mo%emit_nmhc_mo_m
emit_h2_mo_m          =>ctem_tile_mo%emit_h2_mo_m
emit_nox_mo_m         =>ctem_tile_mo%emit_nox_mo_m
emit_n2o_mo_m         =>ctem_tile_mo%emit_n2o_mo_m
emit_pm25_mo_m        =>ctem_tile_mo%emit_pm25_mo_m
emit_tpm_mo_m         =>ctem_tile_mo%emit_tpm_mo_m
emit_tc_mo_m          =>ctem_tile_mo%emit_tc_mo_m
emit_oc_mo_m          =>ctem_tile_mo%emit_oc_mo_m
emit_bc_mo_m          =>ctem_tile_mo%emit_bc_mo_m
burnfrac_mo_m         =>ctem_tile_mo%burnfrac_mo_m
probfire_mo_m         =>ctem_tile_mo%probfire_mo_m
bterm_mo_m            =>ctem_tile_mo%bterm_mo_m
luc_emc_mo_m          =>ctem_tile_mo%luc_emc_mo_m
lterm_mo_m            =>ctem_tile_mo%lterm_mo_m
lucsocin_mo_m         =>ctem_tile_mo%lucsocin_mo_m
mterm_mo_m            =>ctem_tile_mo%mterm_mo_m
lucltrin_mo_m         =>ctem_tile_mo%lucltrin_mo_m
ch4wet1_mo_m          =>ctem_tile_mo%ch4wet1_mo_m
ch4wet2_mo_m          =>ctem_tile_mo%ch4wet2_mo_m
wetfdyn_mo_m          =>ctem_tile_mo%wetfdyn_mo_m
ch4dyn1_mo_m          =>ctem_tile_mo%ch4dyn1_mo_m
ch4dyn2_mo_m          =>ctem_tile_mo%ch4dyn2_mo_m

gppvegrow         => vrot%gppveg
nepvegrow         => vrot%nepveg
nbpvegrow         => vrot%nbpveg
nppvegrow         => vrot%nppveg
hetroresvegrow    => vrot%hetroresveg
autoresvegrow     => vrot%autoresveg
litresvegrow      => vrot%litresveg
soilcresvegrow    => vrot%soilcresveg
rmlvegaccrow      => vrot%rmlvegacc
rmsvegrow         => vrot%rmsveg
rmrvegrow         => vrot%rmrveg
rgvegrow          => vrot%rgveg
ailcgrow          => vrot%ailcg
emit_co2row       => vrot%emit_co2
emit_corow        => vrot%emit_co
emit_ch4row       => vrot%emit_ch4
emit_nmhcrow      => vrot%emit_nmhc
emit_h2row        => vrot%emit_h2
emit_noxrow       => vrot%emit_nox
emit_n2orow       => vrot%emit_n2o
emit_pm25row      => vrot%emit_pm25
emit_tpmrow       => vrot%emit_tpm
emit_tcrow        => vrot%emit_tc
emit_ocrow        => vrot%emit_oc
emit_bcrow        => vrot%emit_bc
burnfracrow       => vrot%burnfrac
burnvegfrow       => vrot%burnvegf
probfirerow       => vrot%probfire
btermrow          => vrot%bterm
ltermrow          => vrot%lterm
mtermrow          => vrot%mterm
lucemcomrow       => vrot%lucemcom
lucltrinrow       => vrot%lucltrin
lucsocinrow       => vrot%lucsocin
ch4wet1row        => vrot%ch4wet1
ch4wet2row        => vrot%ch4wet2
wetfdynrow        => vrot%wetfdyn
ch4dyn1row        => vrot%ch4dyn1
ch4dyn2row        => vrot%ch4dyn2
litrmassrow       => vrot%litrmass
soilcmasrow       => vrot%soilcmas
vgbiomas_vegrow   => vrot%vgbiomas_veg
stemmassrow       => vrot%stemmass
rootmassrow       => vrot%rootmass

laimaxg_mo_g        =>ctem_grd_mo%laimaxg_mo_g
stemmass_mo_g       =>ctem_grd_mo%stemmass_mo_g
rootmass_mo_g       =>ctem_grd_mo%rootmass_mo_g
litrmass_mo_g       =>ctem_grd_mo%litrmass_mo_g
soilcmas_mo_g       =>ctem_grd_mo%soilcmas_mo_g
npp_mo_g            =>ctem_grd_mo%npp_mo_g
gpp_mo_g            =>ctem_grd_mo%gpp_mo_g
nep_mo_g            =>ctem_grd_mo%nep_mo_g
nbp_mo_g            =>ctem_grd_mo%nbp_mo_g
hetrores_mo_g       =>ctem_grd_mo%hetrores_mo_g
autores_mo_g        =>ctem_grd_mo%autores_mo_g
litres_mo_g         =>ctem_grd_mo%litres_mo_g
soilcres_mo_g       =>ctem_grd_mo%soilcres_mo_g
vgbiomas_mo_g       =>ctem_grd_mo%vgbiomas_mo_g
totcmass_mo_g       =>ctem_grd_mo%totcmass_mo_g
emit_co2_mo_g       =>ctem_grd_mo%emit_co2_mo_g
emit_co_mo_g        =>ctem_grd_mo%emit_co_mo_g
emit_ch4_mo_g       =>ctem_grd_mo%emit_ch4_mo_g
emit_nmhc_mo_g      =>ctem_grd_mo%emit_nmhc_mo_g
emit_h2_mo_g        =>ctem_grd_mo%emit_h2_mo_g
emit_nox_mo_g       =>ctem_grd_mo%emit_nox_mo_g
emit_n2o_mo_g       =>ctem_grd_mo%emit_n2o_mo_g
emit_pm25_mo_g      =>ctem_grd_mo%emit_pm25_mo_g
emit_tpm_mo_g       =>ctem_grd_mo%emit_tpm_mo_g
emit_tc_mo_g        =>ctem_grd_mo%emit_tc_mo_g
emit_oc_mo_g        =>ctem_grd_mo%emit_oc_mo_g
emit_bc_mo_g        =>ctem_grd_mo%emit_bc_mo_g
probfire_mo_g       =>ctem_grd_mo%probfire_mo_g
luc_emc_mo_g        =>ctem_grd_mo%luc_emc_mo_g
lucltrin_mo_g       =>ctem_grd_mo%lucltrin_mo_g
lucsocin_mo_g       =>ctem_grd_mo%lucsocin_mo_g
burnfrac_mo_g       =>ctem_grd_mo%burnfrac_mo_g
bterm_mo_g          =>ctem_grd_mo%bterm_mo_g
lterm_mo_g          =>ctem_grd_mo%lterm_mo_g
mterm_mo_g          =>ctem_grd_mo%mterm_mo_g
ch4wet1_mo_g        =>ctem_grd_mo%ch4wet1_mo_g
ch4wet2_mo_g        =>ctem_grd_mo%ch4wet2_mo_g
wetfdyn_mo_g        =>ctem_grd_mo%wetfdyn_mo_g
ch4dyn1_mo_g        =>ctem_grd_mo%ch4dyn1_mo_g
ch4dyn2_mo_g        =>ctem_grd_mo%ch4dyn2_mo_g


! ------------

!       accumulate monthly outputs
!
        do 862 i=1,nltest

         do 863 m=1,nmtest

          do j=1,icc

           if (ailcgrow(i,m,j) .gt. laimaxg_mo_m(i,m,j)) then
            laimaxg_mo_m(i,m,j)=ailcgrow(i,m,j)
           end if

           npp_mo_m(i,m,j)=npp_mo_m(i,m,j)+nppvegrow(i,m,j)
           gpp_mo_m(i,m,j)=gpp_mo_m(i,m,j)+gppvegrow(i,m,j)
           nep_mo_m(i,m,j)=nep_mo_m(i,m,j)+nepvegrow(i,m,j)
           nbp_mo_m(i,m,j)=nbp_mo_m(i,m,j)+nbpvegrow(i,m,j)
           hetrores_mo_m(i,m,j)=hetrores_mo_m(i,m,j)+hetroresvegrow(i,m,j)
           autores_mo_m(i,m,j) =autores_mo_m(i,m,j)+autoresvegrow(i,m,j)
           litres_mo_m(i,m,j)  =litres_mo_m(i,m,j) +litresvegrow(i,m,j)
           soilcres_mo_m(i,m,j) =soilcres_mo_m(i,m,j) +soilcresvegrow(i,m,j)
           emit_co2_mo_m(i,m,j)=emit_co2_mo_m(i,m,j)+emit_co2row(i,m,j)
           emit_co_mo_m(i,m,j) =emit_co_mo_m(i,m,j)+emit_corow(i,m,j)
           emit_ch4_mo_m(i,m,j) =emit_ch4_mo_m(i,m,j)+emit_ch4row(i,m,j)
           emit_nmhc_mo_m(i,m,j)=emit_nmhc_mo_m(i,m,j)+emit_nmhcrow(i,m,j)
           emit_h2_mo_m(i,m,j) =emit_h2_mo_m(i,m,j)+emit_h2row(i,m,j)
           emit_nox_mo_m(i,m,j) =emit_nox_mo_m(i,m,j)+emit_noxrow(i,m,j)
           emit_n2o_mo_m(i,m,j) =emit_n2o_mo_m(i,m,j)+emit_n2orow(i,m,j)
           emit_pm25_mo_m(i,m,j)=emit_pm25_mo_m(i,m,j)+emit_pm25row(i,m,j)
           emit_tpm_mo_m(i,m,j) =emit_tpm_mo_m(i,m,j)+emit_tpmrow(i,m,j)
           emit_tc_mo_m(i,m,j) =emit_tc_mo_m(i,m,j)+emit_tcrow(i,m,j)
           emit_oc_mo_m(i,m,j) =emit_oc_mo_m(i,m,j)+emit_ocrow(i,m,j)
           emit_bc_mo_m(i,m,j) =emit_bc_mo_m(i,m,j)+emit_bcrow(i,m,j)
           burnfrac_mo_m(i,m,j) =burnfrac_mo_m(i,m,j)+burnvegfrow(i,m,j)

          end do

           nep_mo_m(i,m,iccp1)=nep_mo_m(i,m,iccp1)+nepvegrow(i,m,iccp1)
           nbp_mo_m(i,m,iccp1)=nbp_mo_m(i,m,iccp1)+nbpvegrow(i,m,iccp1)
           hetrores_mo_m(i,m,iccp1)=hetrores_mo_m(i,m,iccp1)+hetroresvegrow(i,m,iccp1)
           litres_mo_m(i,m,iccp1)  =litres_mo_m(i,m,iccp1)+litresvegrow(i,m,iccp1)
           soilcres_mo_m(i,m,iccp1) =soilcres_mo_m(i,m,iccp1) +soilcresvegrow(i,m,iccp1)

           luc_emc_mo_m(i,m) =luc_emc_mo_m(i,m)+lucemcomrow(i,m)
           lucsocin_mo_m(i,m) =lucsocin_mo_m(i,m)+lucsocinrow(i,m)
           lucltrin_mo_m(i,m) =lucltrin_mo_m(i,m)+lucltrinrow(i,m)
                         !CH4 related variables !Rudra
           ch4wet1_mo_m(i,m) = ch4wet1_mo_m(i,m) + CH4WET1ROW(i,m)
           ch4wet2_mo_m(i,m) = ch4wet2_mo_m(i,m) + CH4WET2ROW(i,m)
           wetfdyn_mo_m(i,m) = wetfdyn_mo_m(i,m) + WETFDYNROW(i,m)
           ch4dyn1_mo_m(i,m) = ch4dyn1_mo_m(i,m) + CH4DYN1ROW(i,m)
           ch4dyn2_mo_m(i,m) = ch4dyn2_mo_m(i,m) + CH4DYN2ROW(i,m)

!          Sum the probfire now, later we will make it a per day value.
           probfire_mo_m(i,m) =probfire_mo_m(i,m) + probfirerow(i,m)
           bterm_mo_m(i,m) = bterm_mo_m(i,m) + btermrow(i,m)
           lterm_mo_m(i,m) = lterm_mo_m(i,m) + ltermrow(i,m)
           mterm_mo_m(i,m) = mterm_mo_m(i,m) + mtermrow(i,m)

           do 865 nt=1,nmon

             if(iday.eq.mmday(nt))then

               do j=1,icc
                vgbiomas_mo_m(i,m,j)=0.0
                litrmass_mo_m(i,m,j)=0.0
                soilcmas_mo_m(i,m,j)=0.0
                totcmass_mo_m(i,m,j)=0.0
                stemmass_mo_m(i,m,j)=0.0
                rootmass_mo_m(i,m,j)=0.0
               end do
                litrmass_mo_m(i,m,iccp1)=0.0
                soilcmas_mo_m(i,m,iccp1)=0.0

                do 867 j=1,icc

                  vgbiomas_mo_m(i,m,j)=vgbiomas_vegrow(i,m,j)
                  litrmass_mo_m(i,m,j)=litrmassrow(i,m,j)
                  soilcmas_mo_m(i,m,j)=soilcmasrow(i,m,j)
                  stemmass_mo_m(i,m,j)=stemmassrow(i,m,j)
                  rootmass_mo_m(i,m,j)=rootmassrow(i,m,j)
                  totcmass_mo_m(i,m,j)=vgbiomas_vegrow(i,m,j) + litrmassrow(i,m,j)+soilcmasrow(i,m,j)

867             continue

                ! Do the bare fraction too
                litrmass_mo_m(i,m,iccp1)=litrmassrow(i,m,iccp1)
                soilcmas_mo_m(i,m,iccp1)=soilcmasrow(i,m,iccp1)

                barefrac=1.0

               do j=1,icc
                vgbiomas_mo_g(i)=vgbiomas_mo_g(i)+vgbiomas_mo_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                litrmass_mo_g(i)=litrmass_mo_g(i)+litrmass_mo_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                soilcmas_mo_g(i)=soilcmas_mo_g(i)+soilcmas_mo_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                stemmass_mo_g(i)=stemmass_mo_g(i)+stemmass_mo_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                rootmass_mo_g(i)=rootmass_mo_g(i)+rootmass_mo_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                totcmass_mo_g(i)=totcmass_mo_g(i)+totcmass_mo_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                barefrac=barefrac-FAREROT(i,m)*fcancmxrow(i,m,j)

               end do

!             Also add in the bare fraction contributions.
              litrmass_mo_g(i)=litrmass_mo_g(i)+litrmass_mo_m(i,m,iccp1)*barefrac
              soilcmas_mo_g(i)=soilcmas_mo_g(i)+soilcmas_mo_m(i,m,iccp1)*barefrac

             endif ! mmday (mid-month instantaneous value)

             if(iday.eq.monthend(nt+1))then

               ndmonth=(monthend(nt+1)-monthend(nt))*nday

               barefrac=1.0

               do j=1,icc

                npp_mo_g(i)=npp_mo_g(i)+npp_mo_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                gpp_mo_g(i)=gpp_mo_g(i)+gpp_mo_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                nep_mo_g(i)=nep_mo_g(i)+nep_mo_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                nbp_mo_g(i)=nbp_mo_g(i)+nbp_mo_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                hetrores_mo_g(i)=hetrores_mo_g(i)+hetrores_mo_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                autores_mo_g(i) =autores_mo_g(i) +autores_mo_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                litres_mo_g(i)  =litres_mo_g(i) +litres_mo_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                soilcres_mo_g(i) =soilcres_mo_g(i)+ soilcres_mo_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                laimaxg_mo_g(i)=laimaxg_mo_g(i)+laimaxg_mo_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                emit_co2_mo_g(i)=emit_co2_mo_g(i)+emit_co2_mo_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                emit_co_mo_g(i) =emit_co_mo_g(i)+emit_co_mo_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                emit_ch4_mo_g(i) =emit_ch4_mo_g(i)+emit_ch4_mo_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                emit_nmhc_mo_g(i)=emit_nmhc_mo_g(i)+emit_nmhc_mo_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                emit_h2_mo_g(i) =emit_h2_mo_g(i)+emit_h2_mo_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                emit_nox_mo_g(i) =emit_nox_mo_g(i)+emit_nox_mo_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                emit_n2o_mo_g(i) =emit_n2o_mo_g(i)+emit_n2o_mo_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                emit_pm25_mo_g(i) =emit_pm25_mo_g(i)+emit_pm25_mo_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                emit_tpm_mo_g(i) =emit_tpm_mo_g(i)+emit_tpm_mo_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                emit_tc_mo_g(i) =emit_tc_mo_g(i)+emit_tc_mo_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                emit_oc_mo_g(i) =emit_oc_mo_g(i)+emit_oc_mo_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                emit_bc_mo_g(i) =emit_bc_mo_g(i)+emit_bc_mo_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                burnfrac_mo_g(i)=burnfrac_mo_g(i)+burnfrac_mo_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                barefrac=barefrac-FAREROT(i,m)*fcancmxrow(i,m,j)

               end do !j

              nep_mo_g(i)=nep_mo_g(i)+nep_mo_m(i,m,iccp1)*barefrac
              nbp_mo_g(i)=nbp_mo_g(i)+nbp_mo_m(i,m,iccp1)*barefrac
              hetrores_mo_g(i)=hetrores_mo_g(i)+hetrores_mo_m(i,m,iccp1)*barefrac
              litres_mo_g(i)  =litres_mo_g(i) +litres_mo_m(i,m,iccp1)*barefrac
              soilcres_mo_g(i)=soilcres_mo_g(i)+soilcres_mo_m(i,m,iccp1)*barefrac

               luc_emc_mo_g(i) =luc_emc_mo_g(i)+luc_emc_mo_m(i,m)*FAREROT(i,m)
               lucsocin_mo_g(i) =lucsocin_mo_g(i)+lucsocin_mo_m(i,m)*FAREROT(i,m)
               lucltrin_mo_g(i) =lucltrin_mo_g(i)+lucltrin_mo_m(i,m)*FAREROT(i,m)

               ch4wet1_mo_g(i) = ch4wet1_mo_g(i) +ch4wet1_mo_m(i,m)*FAREROT(i,m)
               ch4wet2_mo_g(i) = ch4wet2_mo_g(i)+ch4wet2_mo_m(i,m)*FAREROT(i,m)

               wetfdyn_mo_m(i,m)=wetfdyn_mo_m(i,m)*(1./real(monthdays(nt)))

               wetfdyn_mo_g(i) = wetfdyn_mo_g(i)+wetfdyn_mo_m(i,m)*FAREROT(i,m)
               ch4dyn1_mo_g(i) = ch4dyn1_mo_g(i)+ch4dyn1_mo_m(i,m)*FAREROT(i,m)
               ch4dyn2_mo_g(i) = ch4dyn2_mo_g(i)+ch4dyn2_mo_m(i,m)*FAREROT(i,m)

!              Make the probability of fire a per day value
               probfire_mo_m(i,m)=probfire_mo_m(i,m)*(1./real(monthdays(nt)))
               probfire_mo_g(i)=probfire_mo_g(i)+probfire_mo_m(i,m)*FAREROT(i,m)
               bterm_mo_m(i,m)=bterm_mo_m(i,m)*(1./real(monthdays(nt)))
               bterm_mo_g(i) =bterm_mo_g(i)+bterm_mo_m(i,m)*FAREROT(i,m)
               lterm_mo_m(i,m)=lterm_mo_m(i,m)*(1./real(monthdays(nt)))
               lterm_mo_g(i) =lterm_mo_g(i)+lterm_mo_m(i,m)*FAREROT(i,m)
               mterm_mo_m(i,m)=mterm_mo_m(i,m)*(1./real(monthdays(nt)))
               mterm_mo_g(i) =mterm_mo_g(i)+mterm_mo_m(i,m)*FAREROT(i,m)


             endif ! monthend (max lai and accumulated npp/gpp/nep over the whole month)
                   ! if(iday.eq.monthend(nt+1))

865        continue ! nmon
863      continue ! m

         do nt=1,nmon
           if(iday.eq.monthend(nt+1))then
             imonth=nt

                barefrac=1.0

!           Write to file .CT01M_M/.CT01M_G
              do m=1,nmtest
               do j=1,icc

                  barefrac=barefrac-fcancmxrow(i,m,j)*FAREROT(i,m)

               if (FAREROT(i,m)*fcancmxrow(i,m,j) .gt. seed) then
                 write(84,8104)imonth,iyear,laimaxg_mo_m(i,m,j),&
                    vgbiomas_mo_m(i,m,j),litrmass_mo_m(i,m,j),&
                    soilcmas_mo_m(i,m,j),npp_mo_m(i,m,j),&
                    gpp_mo_m(i,m,j),nep_mo_m(i,m,j),&
                    nbp_mo_m(i,m,j),hetrores_mo_m(i,m,j),&
                    autores_mo_m(i,m,j),litres_mo_m(i,m,j),&
                    soilcres_mo_m(i,m,j),&
                    ' TILE ',m,' PFT ',j,' FRAC ',FAREROT(i,m)* &
                    fcancmxrow(i,m,j)
               end if
               end do !icc

               if (m .eq. nmtest) then
                if (barefrac .gt. seed) then
                write(84,8104)imonth,iyear,0.0,  &
                    0.0,litrmass_mo_m(i,m,iccp1), &
                    soilcmas_mo_m(i,m,iccp1),0.0, &
                    0.0,0.0, &
                    0.0,hetrores_mo_m(i,m,iccp1), &
                    0.0,litres_mo_m(i,m,iccp1), &
                    soilcres_mo_m(i,m,iccp1), &
                    ' TILE ',m,' PFT ',iccp1,' FRAC ',barefrac
                end if
               end if
              end do !m

              write(84,8104)imonth,iyear,laimaxg_mo_g(i), &
                     vgbiomas_mo_g(i),litrmass_mo_g(i), &
                    soilcmas_mo_g(i),npp_mo_g(i), &
                    gpp_mo_g(i),nep_mo_g(i), &
                    nbp_mo_g(i),hetrores_mo_g(i),autores_mo_g(i), &
                    litres_mo_g(i),soilcres_mo_g(i),' GRDAV'

            if (dofire .or. lnduseon) then

!            write to file .CT06M_M/.CT06M_G

              do m=1,nmtest
               do j=1,icc
                if (FAREROT(i,m)*fcancmxrow(i,m,j) .gt. seed) then
                 write(85,8109)imonth,iyear,emit_co2_mo_m(i,m,j), &
                    emit_co_mo_m(i,m,j),emit_ch4_mo_m(i,m,j), &
                    emit_nmhc_mo_m(i,m,j),emit_h2_mo_m(i,m,j), &
                    emit_nox_mo_m(i,m,j),emit_n2o_mo_m(i,m,j), &
                    emit_pm25_mo_m(i,m,j),emit_tpm_mo_m(i,m,j), &
                    emit_tc_mo_m(i,m,j),emit_oc_mo_m(i,m,j), &
                    emit_bc_mo_m(i,m,j),probfire_mo_m(i,m), &
                    luc_emc_mo_m(i,m),lucltrin_mo_m(i,m), &
                    lucsocin_mo_m(i,m),burnfrac_mo_m(i,m,j)*100., &
                    bterm_mo_m(i,m),lterm_mo_m(i,m),mterm_mo_m(i,m), &
                    ' TILE ',m,' PFT ',j,' FRAC ',FAREROT(i,m)* &
                    fcancmxrow(i,m,j)
                end if
               end do
              end do

             write(85,8109)imonth,iyear,emit_co2_mo_g(i), &
                    emit_co_mo_g(i),emit_ch4_mo_g(i),emit_nmhc_mo_g(i), &
                    emit_h2_mo_g(i),emit_nox_mo_g(i),emit_n2o_mo_g(i), &
                    emit_pm25_mo_g(i),emit_tpm_mo_g(i),emit_tc_mo_g(i), &
                    emit_oc_mo_g(i),emit_bc_mo_g(i), &
                    probfire_mo_g(i),luc_emc_mo_g(i), &
                    lucltrin_mo_g(i),lucsocin_mo_g(i), &
                    burnfrac_mo_g(i)*100.,bterm_mo_g(i),lterm_mo_g(i), &
                    mterm_mo_g(i),' GRDAV '

            endif  !dofire/lnduseon

!           add fraction of each pft and bare \\

            if (compete .or. lnduseon) then
              sumfare=0.0
              if (mosaic) then
               do m=1,nmos
                 sumfare=sumfare+FAREROT(i,m)
               enddo
               write(88,8106)imonth,iyear,(FAREROT(i,m)*100.,m=1,nmos) &
                           ,sumfare,(pftexistrow(i,j,j),j=1,icc)
              else !composite
               m=1
               do j=1,icc
                 sumfare=sumfare+fcancmxrow(i,m,j)
               enddo
               write(88,8106)imonth,iyear,(fcancmxrow(i,m,j)*100., &
                           j=1,icc),(1.0-sumfare)*100.,sumfare, &
                             (pftexistrow(i,m,j),j=1,icc)
              endif !mosaic/composite
            endif !compete/lnduseon

             if (dowetlands .or. obswetf) then
             write(91,8111)imonth,iyear,ch4wet1_mo_g(i), &
                          ch4wet2_mo_g(i),wetfdyn_mo_g(i), &
                          ch4dyn1_mo_g(i),ch4dyn2_mo_g(i)
             endif


!              initialize monthly accumulated arrays
!              for the next round

             do m=1,nmtest

               probfire_mo_m(i,m) =0.0
               luc_emc_mo_m(i,m) =0.0
               lucsocin_mo_m(i,m) =0.0
               lucltrin_mo_m(i,m) =0.0
               bterm_mo_m(i,m) =0.0
               lterm_mo_m(i,m) =0.0
               mterm_mo_m(i,m) =0.0

               ch4wet1_mo_m(i,m)  =0.0
               ch4wet2_mo_m(i,m)  =0.0
               wetfdyn_mo_m(i,m)  =0.0
               ch4dyn1_mo_m(i,m)  =0.0
               ch4dyn2_mo_m(i,m)  =0.0


             do j=1,icc

              laimaxg_mo_m(i,m,j)=0.0
              hetrores_mo_m(i,m,j)=0.0
              autores_mo_m(i,m,j)=0.0
              litres_mo_m(i,m,j)=0.0
              soilcres_mo_m(i,m,j)=0.0

              npp_mo_m(i,m,j)=0.0
              gpp_mo_m(i,m,j)=0.0
              nep_mo_m(i,m,j)=0.0
              nbp_mo_m(i,m,j)=0.0
              emit_co2_mo_m(i,m,j)=0.0
              emit_co_mo_m(i,m,j) =0.0
              emit_ch4_mo_m(i,m,j) =0.0
              emit_nmhc_mo_m(i,m,j) =0.0
              emit_h2_mo_m(i,m,j) =0.0
              emit_nox_mo_m(i,m,j) =0.0
              emit_n2o_mo_m(i,m,j) =0.0
              emit_pm25_mo_m(i,m,j) =0.0
              emit_tpm_mo_m(i,m,j) =0.0
              emit_tc_mo_m(i,m,j) =0.0
              emit_oc_mo_m(i,m,j) =0.0
              emit_bc_mo_m(i,m,j) =0.0
              burnfrac_mo_m(i,m,j) =0.0
             enddo !j

              hetrores_mo_m(i,m,iccp1)=0.0
              litres_mo_m(i,m,iccp1)=0.0
              soilcres_mo_m(i,m,iccp1)=0.0
              nep_mo_m(i,m,iccp1)=0.0
              nbp_mo_m(i,m,iccp1)=0.0


            enddo !m

           endif ! if(iday.eq.monthend(nt+1))
         enddo ! nt=1,nmon

 862     continue ! i

8104  FORMAT(1X,I4,I5,12(F10.3,1X),2(A6,I2),A6,F8.2)
8105  FORMAT(1X,I5,15(F10.3,1X),2(A6,I2),A6,F8.2)
8106  FORMAT(1X,I4,I5,11(F10.5,1X),9L5,2(A6,I2))
8107  FORMAT(1X,I5,11(F10.5,1X),9L5,2(A6,I2))
8108  FORMAT(1X,I5,20(F10.3,1X),2(A6,I2),A6,F8.2)
8109  FORMAT(1X,I4,I5,20(F10.3,1X),2(A6,I2),A6,F8.2)
8111  FORMAT(1X,I4,I5,5(F10.3,1X),2(A6,I2))
8115  FORMAT(1X,I5,5(F10.3,1X),2(A6,I2))


end subroutine ctem_monthly_aw

!==============================================================================================================

subroutine ctem_annual_aw(nltest,nmtest,iday,FAREROT,iyear)

use ctem_statevars,     only : ctem_tile_yr, vrot, ctem_grd_yr, c_switch
use ctem_params, only : icc,iccp1,seed,nmos

implicit none

! arguments
integer, intent(in) :: nltest
integer, intent(in) :: nmtest
integer, intent(in) :: iday
real, intent(in), dimension(:,:) :: FAREROT
integer, intent(in) :: iyear

! pointers

logical, pointer :: dofire
logical, pointer :: lnduseon
logical, pointer :: compete
logical, pointer :: dowetlands
logical, pointer :: obswetf
logical, pointer :: mosaic

real, pointer, dimension(:,:,:) :: laimaxg_yr_m
real, pointer, dimension(:,:,:) :: stemmass_yr_m
real, pointer, dimension(:,:,:) :: rootmass_yr_m
real, pointer, dimension(:,:,:) :: npp_yr_m
real, pointer, dimension(:,:,:) :: gpp_yr_m
real, pointer, dimension(:,:,:) :: vgbiomas_yr_m
real, pointer, dimension(:,:,:) :: autores_yr_m
real, pointer, dimension(:,:,:) :: totcmass_yr_m
real, pointer, dimension(:,:,:) :: litrmass_yr_m
real, pointer, dimension(:,:,:) :: soilcmas_yr_m
real, pointer, dimension(:,:,:) :: nep_yr_m
real, pointer, dimension(:,:,:) :: litres_yr_m
real, pointer, dimension(:,:,:) :: soilcres_yr_m
real, pointer, dimension(:,:,:) :: hetrores_yr_m
real, pointer, dimension(:,:,:) :: nbp_yr_m
real, pointer, dimension(:,:,:) :: emit_co2_yr_m
real, pointer, dimension(:,:,:) :: emit_co_yr_m
real, pointer, dimension(:,:,:) :: emit_ch4_yr_m
real, pointer, dimension(:,:,:) :: emit_nmhc_yr_m
real, pointer, dimension(:,:,:) :: emit_h2_yr_m
real, pointer, dimension(:,:,:) :: emit_nox_yr_m
real, pointer, dimension(:,:,:) :: emit_n2o_yr_m
real, pointer, dimension(:,:,:) :: emit_pm25_yr_m
real, pointer, dimension(:,:,:) :: emit_tpm_yr_m
real, pointer, dimension(:,:,:) :: emit_tc_yr_m
real, pointer, dimension(:,:,:) :: emit_oc_yr_m
real, pointer, dimension(:,:,:) :: emit_bc_yr_m
real, pointer, dimension(:,:,:) :: burnfrac_yr_m
real, pointer, dimension(:,:) :: probfire_yr_m
real, pointer, dimension(:,:) :: bterm_yr_m
real, pointer, dimension(:,:) :: luc_emc_yr_m
real, pointer, dimension(:,:) :: lterm_yr_m
real, pointer, dimension(:,:) :: lucsocin_yr_m
real, pointer, dimension(:,:) :: mterm_yr_m
real, pointer, dimension(:,:) :: lucltrin_yr_m
real, pointer, dimension(:,:) :: ch4wet1_yr_m
real, pointer, dimension(:,:) :: ch4wet2_yr_m
real, pointer, dimension(:,:) :: wetfdyn_yr_m
real, pointer, dimension(:,:) :: ch4dyn1_yr_m
real, pointer, dimension(:,:) :: ch4dyn2_yr_m

logical, pointer, dimension(:,:,:) :: pftexistrow
real, pointer, dimension(:,:,:) :: gppvegrow
real, pointer, dimension(:,:,:) :: nepvegrow
real, pointer, dimension(:,:,:) :: nbpvegrow
real, pointer, dimension(:,:,:) :: nppvegrow
real, pointer, dimension(:,:,:) :: hetroresvegrow
real, pointer, dimension(:,:,:) :: autoresvegrow
real, pointer, dimension(:,:,:) :: litresvegrow
real, pointer, dimension(:,:,:) :: soilcresvegrow
real, pointer, dimension(:,:,:) :: rmlvegaccrow
real, pointer, dimension(:,:,:) :: rmsvegrow
real, pointer, dimension(:,:,:) :: rmrvegrow
real, pointer, dimension(:,:,:) :: rgvegrow
real, pointer, dimension(:,:,:) :: ailcgrow
real, pointer, dimension(:,:,:) :: emit_co2row
real, pointer, dimension(:,:,:) :: emit_corow
real, pointer, dimension(:,:,:) :: emit_ch4row
real, pointer, dimension(:,:,:) :: emit_nmhcrow
real, pointer, dimension(:,:,:) :: emit_h2row
real, pointer, dimension(:,:,:) :: emit_noxrow
real, pointer, dimension(:,:,:) :: emit_n2orow
real, pointer, dimension(:,:,:) :: emit_pm25row
real, pointer, dimension(:,:,:) :: emit_tpmrow
real, pointer, dimension(:,:,:) :: emit_tcrow
real, pointer, dimension(:,:,:) :: emit_ocrow
real, pointer, dimension(:,:,:) :: emit_bcrow
real, pointer, dimension(:,:) :: burnfracrow
real, pointer, dimension(:,:,:) :: burnvegfrow
real, pointer, dimension(:,:) :: probfirerow
real, pointer, dimension(:,:) :: btermrow
real, pointer, dimension(:,:) :: ltermrow
real, pointer, dimension(:,:) :: mtermrow
real, pointer, dimension(:,:) :: lucemcomrow
real, pointer, dimension(:,:) :: lucltrinrow
real, pointer, dimension(:,:) :: lucsocinrow
real, pointer, dimension(:,:) :: ch4wet1row
real, pointer, dimension(:,:) :: ch4wet2row
real, pointer, dimension(:,:) :: wetfdynrow
real, pointer, dimension(:,:) :: ch4dyn1row
real, pointer, dimension(:,:) :: ch4dyn2row
real, pointer, dimension(:,:,:) :: litrmassrow
real, pointer, dimension(:,:,:) :: soilcmasrow
real, pointer, dimension(:,:,:) :: vgbiomas_vegrow
real, pointer, dimension(:,:,:) :: stemmassrow
real, pointer, dimension(:,:,:) :: rootmassrow
real, pointer, dimension(:,:,:) :: fcancmxrow

real, pointer, dimension(:) :: laimaxg_yr_g
real, pointer, dimension(:) :: stemmass_yr_g
real, pointer, dimension(:) :: rootmass_yr_g
real, pointer, dimension(:) :: litrmass_yr_g
real, pointer, dimension(:) :: soilcmas_yr_g
real, pointer, dimension(:) :: npp_yr_g
real, pointer, dimension(:) :: gpp_yr_g
real, pointer, dimension(:) :: nep_yr_g
real, pointer, dimension(:) :: nbp_yr_g
real, pointer, dimension(:) :: hetrores_yr_g
real, pointer, dimension(:) :: autores_yr_g
real, pointer, dimension(:) :: litres_yr_g
real, pointer, dimension(:) :: soilcres_yr_g
real, pointer, dimension(:) :: vgbiomas_yr_g
real, pointer, dimension(:) :: totcmass_yr_g
real, pointer, dimension(:) :: emit_co2_yr_g
real, pointer, dimension(:) :: emit_co_yr_g
real, pointer, dimension(:) :: emit_ch4_yr_g
real, pointer, dimension(:) :: emit_nmhc_yr_g
real, pointer, dimension(:) :: emit_h2_yr_g
real, pointer, dimension(:) :: emit_nox_yr_g
real, pointer, dimension(:) :: emit_n2o_yr_g
real, pointer, dimension(:) :: emit_pm25_yr_g
real, pointer, dimension(:) :: emit_tpm_yr_g
real, pointer, dimension(:) :: emit_tc_yr_g
real, pointer, dimension(:) :: emit_oc_yr_g
real, pointer, dimension(:) :: emit_bc_yr_g
real, pointer, dimension(:) :: probfire_yr_g
real, pointer, dimension(:) :: luc_emc_yr_g
real, pointer, dimension(:) :: lucltrin_yr_g
real, pointer, dimension(:) :: lucsocin_yr_g
real, pointer, dimension(:) :: burnfrac_yr_g
real, pointer, dimension(:) :: bterm_yr_g
real, pointer, dimension(:) :: lterm_yr_g
real, pointer, dimension(:) :: mterm_yr_g
real, pointer, dimension(:) :: ch4wet1_yr_g
real, pointer, dimension(:) :: ch4wet2_yr_g
real, pointer, dimension(:) :: wetfdyn_yr_g
real, pointer, dimension(:) :: ch4dyn1_yr_g
real, pointer, dimension(:) :: ch4dyn2_yr_g

! local
integer :: i,m,j,nt
real :: barefrac
real :: sumfare
integer :: NDMONTH
integer :: IMONTH

! point pointers

dofire                => c_switch%dofire
lnduseon              => c_switch%lnduseon
compete               => c_switch%compete
dowetlands            => c_switch%dowetlands
obswetf               => c_switch%obswetf
mosaic                => c_switch%mosaic

laimaxg_yr_m          =>ctem_tile_yr%laimaxg_yr_m
stemmass_yr_m         =>ctem_tile_yr%stemmass_yr_m
rootmass_yr_m         =>ctem_tile_yr%rootmass_yr_m
npp_yr_m              =>ctem_tile_yr%npp_yr_m
gpp_yr_m              =>ctem_tile_yr%gpp_yr_m
vgbiomas_yr_m         =>ctem_tile_yr%vgbiomas_yr_m
autores_yr_m          =>ctem_tile_yr%autores_yr_m
totcmass_yr_m         =>ctem_tile_yr%totcmass_yr_m
litrmass_yr_m         =>ctem_tile_yr%litrmass_yr_m
soilcmas_yr_m         =>ctem_tile_yr%soilcmas_yr_m
nep_yr_m              =>ctem_tile_yr%nep_yr_m
litres_yr_m           =>ctem_tile_yr%litres_yr_m
soilcres_yr_m         =>ctem_tile_yr%soilcres_yr_m
hetrores_yr_m         =>ctem_tile_yr%hetrores_yr_m
nbp_yr_m              =>ctem_tile_yr%nbp_yr_m
emit_co2_yr_m         =>ctem_tile_yr%emit_co2_yr_m
emit_co_yr_m          =>ctem_tile_yr%emit_co_yr_m
emit_ch4_yr_m         =>ctem_tile_yr%emit_ch4_yr_m
emit_nmhc_yr_m        =>ctem_tile_yr%emit_nmhc_yr_m
emit_h2_yr_m          =>ctem_tile_yr%emit_h2_yr_m
emit_nox_yr_m         =>ctem_tile_yr%emit_nox_yr_m
emit_n2o_yr_m         =>ctem_tile_yr%emit_n2o_yr_m
emit_pm25_yr_m        =>ctem_tile_yr%emit_pm25_yr_m
emit_tpm_yr_m         =>ctem_tile_yr%emit_tpm_yr_m
emit_tc_yr_m          =>ctem_tile_yr%emit_tc_yr_m
emit_oc_yr_m          =>ctem_tile_yr%emit_oc_yr_m
emit_bc_yr_m          =>ctem_tile_yr%emit_bc_yr_m
burnfrac_yr_m         =>ctem_tile_yr%burnfrac_yr_m
probfire_yr_m         =>ctem_tile_yr%probfire_yr_m
bterm_yr_m            =>ctem_tile_yr%bterm_yr_m
luc_emc_yr_m          =>ctem_tile_yr%luc_emc_yr_m
lterm_yr_m            =>ctem_tile_yr%lterm_yr_m
lucsocin_yr_m         =>ctem_tile_yr%lucsocin_yr_m
mterm_yr_m            =>ctem_tile_yr%mterm_yr_m
lucltrin_yr_m         =>ctem_tile_yr%lucltrin_yr_m
ch4wet1_yr_m          =>ctem_tile_yr%ch4wet1_yr_m
ch4wet2_yr_m          =>ctem_tile_yr%ch4wet2_yr_m
wetfdyn_yr_m          =>ctem_tile_yr%wetfdyn_yr_m
ch4dyn1_yr_m          =>ctem_tile_yr%ch4dyn1_yr_m
ch4dyn2_yr_m          =>ctem_tile_yr%ch4dyn2_yr_m

pftexistrow       => vrot%pftexist
gppvegrow         => vrot%gppveg
nepvegrow         => vrot%nepveg
nbpvegrow         => vrot%nbpveg
nppvegrow         => vrot%nppveg
hetroresvegrow    => vrot%hetroresveg
autoresvegrow     => vrot%autoresveg
litresvegrow      => vrot%litresveg
soilcresvegrow    => vrot%soilcresveg
rmlvegaccrow      => vrot%rmlvegacc
rmsvegrow         => vrot%rmsveg
rmrvegrow         => vrot%rmrveg
rgvegrow          => vrot%rgveg
ailcgrow          => vrot%ailcg
emit_co2row       => vrot%emit_co2
emit_corow        => vrot%emit_co
emit_ch4row       => vrot%emit_ch4
emit_nmhcrow      => vrot%emit_nmhc
emit_h2row        => vrot%emit_h2
emit_noxrow       => vrot%emit_nox
emit_n2orow       => vrot%emit_n2o
emit_pm25row      => vrot%emit_pm25
emit_tpmrow       => vrot%emit_tpm
emit_tcrow        => vrot%emit_tc
emit_ocrow        => vrot%emit_oc
emit_bcrow        => vrot%emit_bc
burnfracrow       => vrot%burnfrac
burnvegfrow       => vrot%burnvegf
probfirerow       => vrot%probfire
btermrow          => vrot%bterm
ltermrow          => vrot%lterm
mtermrow          => vrot%mterm
lucemcomrow       => vrot%lucemcom
lucltrinrow       => vrot%lucltrin
lucsocinrow       => vrot%lucsocin
ch4wet1row        => vrot%ch4wet1
ch4wet2row        => vrot%ch4wet2
wetfdynrow        => vrot%wetfdyn
ch4dyn1row        => vrot%ch4dyn1
ch4dyn2row        => vrot%ch4dyn2
litrmassrow       => vrot%litrmass
soilcmasrow       => vrot%soilcmas
vgbiomas_vegrow   => vrot%vgbiomas_veg
stemmassrow       => vrot%stemmass
rootmassrow       => vrot%rootmass
fcancmxrow        => vrot%fcancmx

laimaxg_yr_g          =>ctem_grd_yr%laimaxg_yr_g
stemmass_yr_g         =>ctem_grd_yr%stemmass_yr_g
rootmass_yr_g         =>ctem_grd_yr%rootmass_yr_g
litrmass_yr_g         =>ctem_grd_yr%litrmass_yr_g
soilcmas_yr_g         =>ctem_grd_yr%soilcmas_yr_g
npp_yr_g              =>ctem_grd_yr%npp_yr_g
gpp_yr_g              =>ctem_grd_yr%gpp_yr_g
nep_yr_g              =>ctem_grd_yr%nep_yr_g
nbp_yr_g              =>ctem_grd_yr%nbp_yr_g
hetrores_yr_g         =>ctem_grd_yr%hetrores_yr_g
autores_yr_g          =>ctem_grd_yr%autores_yr_g
litres_yr_g           =>ctem_grd_yr%litres_yr_g
soilcres_yr_g         =>ctem_grd_yr%soilcres_yr_g
vgbiomas_yr_g         =>ctem_grd_yr%vgbiomas_yr_g
totcmass_yr_g         =>ctem_grd_yr%totcmass_yr_g
emit_co2_yr_g         =>ctem_grd_yr%emit_co2_yr_g
emit_co_yr_g          =>ctem_grd_yr%emit_co_yr_g
emit_ch4_yr_g         =>ctem_grd_yr%emit_ch4_yr_g
emit_nmhc_yr_g        =>ctem_grd_yr%emit_nmhc_yr_g
emit_h2_yr_g          =>ctem_grd_yr%emit_h2_yr_g
emit_nox_yr_g         =>ctem_grd_yr%emit_nox_yr_g
emit_n2o_yr_g         =>ctem_grd_yr%emit_n2o_yr_g
emit_pm25_yr_g        =>ctem_grd_yr%emit_pm25_yr_g
emit_tpm_yr_g         =>ctem_grd_yr%emit_tpm_yr_g
emit_tc_yr_g          =>ctem_grd_yr%emit_tc_yr_g
emit_oc_yr_g          =>ctem_grd_yr%emit_oc_yr_g
emit_bc_yr_g          =>ctem_grd_yr%emit_bc_yr_g
probfire_yr_g         =>ctem_grd_yr%probfire_yr_g
luc_emc_yr_g          =>ctem_grd_yr%luc_emc_yr_g
lucltrin_yr_g         =>ctem_grd_yr%lucltrin_yr_g
lucsocin_yr_g         =>ctem_grd_yr%lucsocin_yr_g
burnfrac_yr_g         =>ctem_grd_yr%burnfrac_yr_g
bterm_yr_g            =>ctem_grd_yr%bterm_yr_g
lterm_yr_g            =>ctem_grd_yr%lterm_yr_g
mterm_yr_g            =>ctem_grd_yr%mterm_yr_g
ch4wet1_yr_g          =>ctem_grd_yr%ch4wet1_yr_g
ch4wet2_yr_g          =>ctem_grd_yr%ch4wet2_yr_g
wetfdyn_yr_g          =>ctem_grd_yr%wetfdyn_yr_g
ch4dyn1_yr_g          =>ctem_grd_yr%ch4dyn1_yr_g
ch4dyn2_yr_g          =>ctem_grd_yr%ch4dyn2_yr_g

!------------

!       accumulate yearly outputs

do 882 i=1,nltest
    do 883 m=1,nmtest
        do 884 j=1,icc

            if (ailcgrow(i,m,j).gt.laimaxg_yr_m(i,m,j)) then
            laimaxg_yr_m(i,m,j)=ailcgrow(i,m,j)
            end if

            npp_yr_m(i,m,j)=npp_yr_m(i,m,j)+nppvegrow(i,m,j)
            gpp_yr_m(i,m,j)=gpp_yr_m(i,m,j)+gppvegrow(i,m,j)
            nep_yr_m(i,m,j)=nep_yr_m(i,m,j)+nepvegrow(i,m,j)
            nbp_yr_m(i,m,j)=nbp_yr_m(i,m,j)+nbpvegrow(i,m,j)
            emit_co2_yr_m(i,m,j)=emit_co2_yr_m(i,m,j)+emit_co2row(i,m,j)
            emit_co_yr_m(i,m,j)=emit_co_yr_m(i,m,j)+emit_corow(i,m,j)
            emit_ch4_yr_m(i,m,j)=emit_ch4_yr_m(i,m,j)+emit_ch4row(i,m,j)
            emit_nmhc_yr_m(i,m,j)=emit_nmhc_yr_m(i,m,j)+emit_nmhcrow(i,m,j)
            emit_h2_yr_m(i,m,j)=emit_h2_yr_m(i,m,j)+emit_h2row(i,m,j)
            emit_nox_yr_m(i,m,j)=emit_nox_yr_m(i,m,j)+emit_noxrow(i,m,j)
            emit_n2o_yr_m(i,m,j)=emit_n2o_yr_m(i,m,j)+emit_n2orow(i,m,j)
            emit_pm25_yr_m(i,m,j)=emit_pm25_yr_m(i,m,j)+emit_pm25row(i,m,j)
            emit_tpm_yr_m(i,m,j)=emit_tpm_yr_m(i,m,j)+emit_tpmrow(i,m,j)
            emit_tc_yr_m(i,m,j)=emit_tc_yr_m(i,m,j)+emit_tcrow(i,m,j)
            emit_oc_yr_m(i,m,j)=emit_oc_yr_m(i,m,j)+emit_ocrow(i,m,j)
            emit_bc_yr_m(i,m,j)=emit_bc_yr_m(i,m,j)+emit_bcrow(i,m,j)

            hetrores_yr_m(i,m,j)=hetrores_yr_m(i,m,j)+hetroresvegrow(i,m,j)
            autores_yr_m(i,m,j)=autores_yr_m(i,m,j)+autoresvegrow(i,m,j)
            litres_yr_m(i,m,j)=litres_yr_m(i,m,j)+litresvegrow(i,m,j)
            soilcres_yr_m(i,m,j)=soilcres_yr_m(i,m,j)+soilcresvegrow(i,m,j)
            burnfrac_yr_m(i,m,j)=burnfrac_yr_m(i,m,j)+burnvegfrow(i,m,j)

884         continue

    !   Also do the bare fraction amounts
        hetrores_yr_m(i,m,iccp1)=hetrores_yr_m(i,m,iccp1)+hetroresvegrow(i,m,iccp1)
        litres_yr_m(i,m,iccp1)=litres_yr_m(i,m,iccp1)+litresvegrow(i,m,iccp1)
        soilcres_yr_m(i,m,iccp1)=soilcres_yr_m(i,m,iccp1)+soilcresvegrow(i,m,iccp1)
        nep_yr_m(i,m,iccp1)=nep_yr_m(i,m,iccp1)+nepvegrow(i,m,iccp1)
        nbp_yr_m(i,m,iccp1)=nbp_yr_m(i,m,iccp1)+nbpvegrow(i,m,iccp1)

        probfire_yr_m(i,m)=probfire_yr_m(i,m)+(probfirerow(i,m) * (1./365.))
        bterm_yr_m(i,m)=bterm_yr_m(i,m)+(btermrow(i,m)*(1./365.))
        lterm_yr_m(i,m)=lterm_yr_m(i,m)+(ltermrow(i,m)*(1./365.))
        mterm_yr_m(i,m)=mterm_yr_m(i,m)+(mtermrow(i,m)*(1./365.))
        luc_emc_yr_m(i,m)=luc_emc_yr_m(i,m)+lucemcomrow(i,m)
        lucsocin_yr_m(i,m)=lucsocin_yr_m(i,m)+lucsocinrow(i,m)
        lucltrin_yr_m(i,m)=lucltrin_yr_m(i,m)+lucltrinrow(i,m)

        ch4wet1_yr_m(i,m) = ch4wet1_yr_m(i,m)+ch4wet1row(i,m)
        ch4wet2_yr_m(i,m) = ch4wet2_yr_m(i,m)+ch4wet2row(i,m)
        wetfdyn_yr_m(i,m) = wetfdyn_yr_m(i,m)+(wetfdynrow(i,m)*(1./365.))
        ch4dyn1_yr_m(i,m) = ch4dyn1_yr_m(i,m)+ch4dyn1row(i,m)
        ch4dyn2_yr_m(i,m) = ch4dyn2_yr_m(i,m)+ch4dyn2row(i,m)


        if (iday.eq.365) then

            do 885 j=1,icc
                stemmass_yr_m(i,m,j)=stemmassrow(i,m,j)
                rootmass_yr_m(i,m,j)=rootmassrow(i,m,j)
                litrmass_yr_m(i,m,j)=litrmassrow(i,m,j)
                soilcmas_yr_m(i,m,j)=soilcmasrow(i,m,j)
                vgbiomas_yr_m(i,m,j)=vgbiomas_vegrow(i,m,j)
                totcmass_yr_m(i,m,j)=vgbiomas_yr_m(i,m,j)+litrmass_yr_m(i,m,j)+soilcmas_yr_m(i,m,j)

885           continue

            litrmass_yr_m(i,m,iccp1)=litrmassrow(i,m,iccp1)
            soilcmas_yr_m(i,m,iccp1)=soilcmasrow(i,m,iccp1)

            barefrac=1.0

            do j=1,icc
                laimaxg_yr_g(i)=laimaxg_yr_g(i)+ laimaxg_yr_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                stemmass_yr_g(i)=stemmass_yr_g(i)+stemmass_yr_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                rootmass_yr_g(i)=rootmass_yr_g(i)+rootmass_yr_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                litrmass_yr_g(i)=litrmass_yr_g(i)+litrmass_yr_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                soilcmas_yr_g(i)=soilcmas_yr_g(i)+soilcmas_yr_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                vgbiomas_yr_g(i)=vgbiomas_yr_g(i)+vgbiomas_yr_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                totcmass_yr_g(i)=totcmass_yr_g(i)+totcmass_yr_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                npp_yr_g(i)=npp_yr_g(i)+npp_yr_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                gpp_yr_g(i)=gpp_yr_g(i)+gpp_yr_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                nep_yr_g(i)=nep_yr_g(i)+nep_yr_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                nbp_yr_g(i)=nbp_yr_g(i)+nbp_yr_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                emit_co2_yr_g(i)=emit_co2_yr_g(i)+emit_co2_yr_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                emit_co_yr_g(i)=emit_co_yr_g(i)+emit_co_yr_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                emit_ch4_yr_g(i)=emit_ch4_yr_g(i)+emit_ch4_yr_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                emit_nmhc_yr_g(i)=emit_nmhc_yr_g(i)+emit_nmhc_yr_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                emit_h2_yr_g(i)=emit_h2_yr_g(i)+emit_h2_yr_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                emit_nox_yr_g(i)=emit_nox_yr_g(i)+emit_nox_yr_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                emit_n2o_yr_g(i)=emit_n2o_yr_g(i)+emit_n2o_yr_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                emit_pm25_yr_g(i)=emit_pm25_yr_g(i)+emit_pm25_yr_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                emit_tpm_yr_g(i)=emit_tpm_yr_g(i)+emit_tpm_yr_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                emit_tc_yr_g(i)=emit_tc_yr_g(i)+emit_tc_yr_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                emit_oc_yr_g(i)=emit_oc_yr_g(i)+emit_oc_yr_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                emit_bc_yr_g(i)=emit_bc_yr_g(i)+emit_bc_yr_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)

                hetrores_yr_g(i)=hetrores_yr_g(i)+hetrores_yr_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                autores_yr_g(i) =autores_yr_g(i) +autores_yr_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                litres_yr_g(i)  =litres_yr_g(i)  +litres_yr_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                soilcres_yr_g(i) =soilcres_yr_g(i) +soilcres_yr_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)
                burnfrac_yr_g(i)=burnfrac_yr_g(i)+burnfrac_yr_m(i,m,j)*FAREROT(i,m)*fcancmxrow(i,m,j)

                barefrac=barefrac-FAREROT(i,m)*fcancmxrow(i,m,j)

            end do !j

            litrmass_yr_g(i)=litrmass_yr_g(i)+litrmass_yr_m(i,m,iccp1)*barefrac
            soilcmas_yr_g(i)=soilcmas_yr_g(i)+soilcmas_yr_m(i,m,iccp1)*barefrac
            hetrores_yr_g(i)=hetrores_yr_g(i)+hetrores_yr_m(i,m,iccp1)*barefrac
            litres_yr_g(i)  =litres_yr_g(i)  +litres_yr_m(i,m,iccp1)*barefrac
            soilcres_yr_g(i)=soilcres_yr_g(i)+soilcres_yr_m(i,m,iccp1)*barefrac
            nep_yr_g(i)=nep_yr_g(i)+nep_yr_m(i,m,j)*barefrac
            nbp_yr_g(i)=nbp_yr_g(i)+nbp_yr_m(i,m,j)*barefrac

            probfire_yr_g(i)=probfire_yr_g(i)+probfire_yr_m(i,m)*FAREROT(i,m)
            bterm_yr_g(i)=bterm_yr_g(i)+bterm_yr_m(i,m)*FAREROT(i,m)
            lterm_yr_g(i)=lterm_yr_g(i)+lterm_yr_m(i,m)*FAREROT(i,m)
            mterm_yr_g(i)=mterm_yr_g(i)+mterm_yr_m(i,m)*FAREROT(i,m)
            luc_emc_yr_g(i)=luc_emc_yr_g(i)+luc_emc_yr_m(i,m)*FAREROT(i,m)
            lucsocin_yr_g(i)=lucsocin_yr_g(i)+lucsocin_yr_m(i,m)*FAREROT(i,m)
            lucltrin_yr_g(i)=lucltrin_yr_g(i)+lucltrin_yr_m(i,m)*FAREROT(i,m)

            ch4wet1_yr_g(i) = ch4wet1_yr_g(i)+ch4wet1_yr_m(i,m)*FAREROT(i,m)
            ch4wet2_yr_g(i) = ch4wet2_yr_g(i)+ch4wet2_yr_m(i,m)*FAREROT(i,m)
            wetfdyn_yr_g(i) = wetfdyn_yr_g(i)+wetfdyn_yr_m(i,m)*FAREROT(i,m)
            ch4dyn1_yr_g(i) = ch4dyn1_yr_g(i)+ch4dyn1_yr_m(i,m)*FAREROT(i,m)
            ch4dyn2_yr_g(i) = ch4dyn2_yr_g(i)+ch4dyn2_yr_m(i,m)*FAREROT(i,m)

    endif ! iday 365

883       continue ! m


    if (iday.eq.365) then

        barefrac=1.0

!   Write to file .CT01Y_M/.CT01Y_G

    do m=1,nmtest
    do j=1,icc

        barefrac=barefrac-fcancmxrow(i,m,j)*FAREROT(i,m)

        if (FAREROT(i,m)*fcancmxrow(i,m,j) .gt. seed) then
        write(86,8105)iyear,laimaxg_yr_m(i,m,j), &
            vgbiomas_yr_m(i,m,j),stemmass_yr_m(i,m,j), &
            rootmass_yr_m(i,m,j),litrmass_yr_m(i,m,j), &
            soilcmas_yr_m(i,m,j),totcmass_yr_m(i,m,j), &
            npp_yr_m(i,m,j),gpp_yr_m(i,m,j),nep_yr_m(i,m,j), &
            nbp_yr_m(i,m,j),hetrores_yr_m(i,m,j), &
            autores_yr_m(i,m,j),litres_yr_m(i,m,j), &
            soilcres_yr_m(i,m,j),' TILE ',m,' PFT ',j,' FRAC ' &
            ,FAREROT(i,m)*fcancmxrow(i,m,j)
        end if
     end do !j

!   Now do the bare fraction of the grid cell. Only soil c, hetres
!   and litter are relevant so the rest are set to 0.

    if (m .eq. nmtest) then
        if (barefrac .gt. seed) then
        write(86,8105)iyear,0.,  &
            0.,  &
            0.,0., &
            litrmassrow(i,m,iccp1),soilcmasrow(i,m,iccp1), &
            0.+soilcmasrow(i,m,iccp1)+ &
            litrmassrow(i,m,iccp1),0., &
            0.,0., &
            0.,hetrores_yr_m(i,m,iccp1), &
            0.,litres_yr_m(i,m,iccp1),soilcres_yr_m(i,m,iccp1), &
            ' TILE ',m,' PFT ',iccp1,' FRAC ',barefrac
        end if
    end if

    end do !m

        write(86,8105)iyear,laimaxg_yr_g(i),vgbiomas_yr_g(i), &
            stemmass_yr_g(i),rootmass_yr_g(i),litrmass_yr_g(i), &
            soilcmas_yr_g(i),totcmass_yr_g(i),npp_yr_g(i), &
            gpp_yr_g(i),nep_yr_g(i), &
            nbp_yr_g(i),hetrores_yr_g(i),autores_yr_g(i), &
            litres_yr_g(i),soilcres_yr_g(i),' GRDAV'

    if (dofire .or. lnduseon) then

!   Write to file .CT06Y_M/.CT06Y_G
    do m=1,nmtest
       do j=1,icc
        if (FAREROT(i,m)*fcancmxrow(i,m,j) .gt. seed) then
        write(87,8108)iyear,emit_co2_yr_m(i,m,j), &
            emit_co_yr_m(i,m,j),emit_ch4_yr_m(i,m,j), &
            emit_nmhc_yr_m(i,m,j),emit_h2_yr_m(i,m,j), &
            emit_nox_yr_m(i,m,j),emit_n2o_yr_m(i,m,j), &
            emit_pm25_yr_m(i,m,j),emit_tpm_yr_m(i,m,j), &
            emit_tc_yr_m(i,m,j),emit_oc_yr_m(i,m,j), &
            emit_bc_yr_m(i,m,j),probfire_yr_m(i,m), &
            luc_emc_yr_m(i,m),lucltrin_yr_m(i,m), &
            lucsocin_yr_m(i,m),burnfrac_yr_m(i,m,j)*100., &
            bterm_yr_m(i,m),lterm_yr_m(i,m),mterm_yr_m(i,m), &
            ' TILE ',m,' PFT ',j,' FRAC ' &
            ,FAREROT(i,m)*fcancmxrow(i,m,j)
        end if
       end do
    end do

        write(87,8108)iyear,emit_co2_yr_g(i), &
            emit_co_yr_g(i),emit_ch4_yr_g(i),emit_nmhc_yr_g(i), &
            emit_h2_yr_g(i),emit_nox_yr_g(i),emit_n2o_yr_g(i), &
            emit_pm25_yr_g(i),emit_tpm_yr_g(i),emit_tc_yr_g(i), &
            emit_oc_yr_g(i),emit_bc_yr_g(i),probfire_yr_g(i), &
            luc_emc_yr_g(i),lucltrin_yr_g(i), &
            lucsocin_yr_g(i),burnfrac_yr_g(i)*100.,bterm_yr_g(i), &
            lterm_yr_g(i),mterm_yr_g(i), ' GRDAV'

    endif !dofire,lnduseon

!       write fraction of each pft and bare

        if (compete .or. lnduseon) then

            sumfare=0.0
            if (mosaic) then
                do m=1,nmos
                    sumfare=sumfare+FAREROT(i,m)
                enddo
                write(89,8107)iyear,(FAREROT(i,m)*100.,m=1,nmos), &
                            sumfare,(pftexistrow(i,j,j),j=1,icc)
            else  !composite
                m=1
                do j=1,icc
                    sumfare=sumfare+fcancmxrow(i,m,j)
                enddo

                write(89,8107)iyear,(fcancmxrow(i,m,j)*100., &
                        j=1,icc),(1.0-sumfare)*100.,sumfare, &
                        (pftexistrow(i,m,j),j=1,icc)

            endif
        endif !compete/lnduseon

        if (dowetlands .or. obswetf) then

        write(92,8115)iyear,ch4wet1_yr_g(i), &
                    ch4wet2_yr_g(i),wetfdyn_yr_g(i), &
                    ch4dyn1_yr_g(i),ch4dyn2_yr_g(i)
        endif

!       initialize yearly accumulated arrays
!       for the next round

        do m=1,nmtest
            probfire_yr_m(i,m)=0.0
            luc_emc_yr_m(i,m)=0.0
            lucsocin_yr_m(i,m)=0.0
            lucltrin_yr_m(i,m)=0.0
            bterm_yr_m(i,m)=0.0
            lterm_yr_m(i,m)=0.0
            mterm_yr_m(i,m)=0.0

            ch4wet1_yr_m(i,m)  =0.0
            ch4wet2_yr_m(i,m)  =0.0
            wetfdyn_yr_m(i,m)  =0.0
            ch4dyn1_yr_m(i,m)  =0.0
            ch4dyn2_yr_m(i,m)  =0.0

            do j = 1, icc
                laimaxg_yr_m(i,m,j)=0.0
                npp_yr_m(i,m,j)=0.0
                gpp_yr_m(i,m,j)=0.0
                nep_yr_m(i,m,j)=0.0
                nbp_yr_m(i,m,j)=0.0
                emit_co2_yr_m(i,m,j)=0.0
                emit_co_yr_m(i,m,j)=0.0
                emit_ch4_yr_m(i,m,j)=0.0
                emit_nmhc_yr_m(i,m,j)=0.0
                emit_h2_yr_m(i,m,j)=0.0
                emit_nox_yr_m(i,m,j)=0.0
                emit_n2o_yr_m(i,m,j)=0.0
                emit_pm25_yr_m(i,m,j)=0.0
                emit_tpm_yr_m(i,m,j)=0.0
                emit_tc_yr_m(i,m,j)=0.0
                emit_oc_yr_m(i,m,j)=0.0
                emit_bc_yr_m(i,m,j)=0.0
                hetrores_yr_m(i,m,j)=0.0
                autores_yr_m(i,m,j)=0.0
                litres_yr_m(i,m,j)=0.0
                soilcres_yr_m(i,m,j)=0.0
                burnfrac_yr_m(i,m,j)=0.0
            enddo

            hetrores_yr_m(i,m,iccp1)=0.0
            litres_yr_m(i,m,iccp1)=0.0
            soilcres_yr_m(i,m,iccp1)=0.0
            nep_yr_m(i,m,iccp1)=0.0
            nbp_yr_m(i,m,iccp1)=0.0
        enddo

    endif ! if iday=365

882     continue ! i

8104  FORMAT(1X,I4,I5,12(F10.3,1X),2(A6,I2),A6,F8.2)
8105  FORMAT(1X,I5,15(F10.3,1X),2(A6,I2),A6,F8.2)
8106  FORMAT(1X,I4,I5,11(F10.5,1X),9L5,2(A6,I2))
8107  FORMAT(1X,I5,11(F10.5,1X),9L5,2(A6,I2))
8108  FORMAT(1X,I5,20(F10.3,1X),2(A6,I2),A6,F8.2)
8109  FORMAT(1X,I4,I5,20(F10.3,1X),2(A6,I2),A6,F8.2)
8111  FORMAT(1X,I4,I5,5(F10.3,1X),2(A6,I2))
8115  FORMAT(1X,I5,5(F10.3,1X),2(A6,I2))

end subroutine ctem_annual_aw

!==============================================================================================================

subroutine close_outfiles()
                     
use ctem_statevars, only : c_switch

implicit none

! pointers

logical, pointer :: dofire
logical, pointer :: lnduseon
logical, pointer :: compete
logical, pointer :: dowetlands
logical, pointer :: obswetf
logical, pointer :: mosaic
logical, pointer :: parallelrun

! point pointers

dofire                => c_switch%dofire
lnduseon              => c_switch%lnduseon
compete               => c_switch%compete
dowetlands            => c_switch%dowetlands
obswetf               => c_switch%obswetf
mosaic                => c_switch%mosaic
parallelrun           => c_switch%parallelrun

! -----------------------

      IF (.NOT. PARALLELRUN) THEN
        close(71)
        close(711)
        close(721)
        close(731)
        close(741)
        close(751)

        if (mosaic) then
         close(72)
         close(73)
         close(74)
         close(75)
         close(76)
         if (dofire .or. lnduseon) then
          close(78)
         end if
        end if

        if (compete .or. lnduseon) then
          close(761)
        endif

        if (dowetlands .or. obswetf) then
        close(762)
        endif 

       if (dofire .or. lnduseon) then
        close(781)
       endif
      endif ! if (.not. parallelrun) 

!     CLOSE CLASS OUTPUT FILES      
      CLOSE(81)
      CLOSE(82)
      CLOSE(83)
!     then the CTEM ones
      close(84)
      close(86)
      if (dofire .or. lnduseon) then
       close(85)
       close(87)
      endif
      if (compete .or. lnduseon) then
       close(88)
       close(89)
      endif
 
      if (dowetlands .or. obswetf) then 
       close(91)
       close(92)
      endif 


end subroutine close_outfiles


end module io_driver      