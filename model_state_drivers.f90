module input_dataset_drivers

! This is the central driver to read in, and write out
! all model state variables (replaces INI and CTM files)
! as well as the model inputs such as MET, population density,
! land use change etc.

! J. Melton
! Nov 2016

implicit none

public :: openmet
!public :: readin_met
public :: read_initialstate

contains

!---

subroutine openmet()

!> The openmet subroutine opens the meteorology netcdf file and prepares it to be
!! read in. The subroutine expects the latitude and longitude dimensions to be called
!! 'lat' and 'lon' in the netcdf, respectively. Additionally the names of the meteorology
!! variables must be the same as those defined in io_driver/vname. The dimensions of this meteorology
!! file is used to determine the number of grid cells that we will running for the
!! model simulation (note - it does not find the number of land cells, only the number of
!! cells that fit into the region of interest). It is vital that all model input files are at the same resolution
!! and grid as the meteorology file since we use it to determine the size of the vectors
!! latvect and lonvect.

! J. Melton
! Nov 2016

use netcdf
use netcdf_drivers, only : check_nc
use io_driver, only : vname,niv,metfid,cntx,cnty,srtx,srty,bounds,lonvect,latvect, &
                      yearmetst
use ctem_statevars,     only : c_switch

implicit none

character(180), pointer :: met_file
!character(180) :: climatefile

integer :: i
integer :: xsize
integer :: ysize
integer :: timelen
integer :: dimid
integer :: varid
integer, dimension(1) :: pos
integer, dimension(2) :: xpos,ypos
real, dimension(1,1,1) :: timestart
character(20) :: unitchar

met_file          => c_switch%met_file

!---
!open climate file

call check_nc(nf90_open(trim(met_file),nf90_nowrite,metfid))

!retrieve dimensions

call check_nc(nf90_inq_dimid(metfid,'lon',dimid))
call check_nc(nf90_inquire_dimension(metfid,dimid,len=xsize))
call check_nc(nf90_inq_dimid(metfid,'lat',dimid))
call check_nc(nf90_inquire_dimension(metfid,dimid,len=ysize))
call check_nc(nf90_inq_dimid(metfid,'time',dimid))
call check_nc(nf90_inquire_dimension(metfid,dimid,len=timelen))

!retrieve scale factors and offsets

do i = 1,niv

  call check_nc(nf90_inq_varid(metfid,vname(i),varid))

    ! If the netcdf has scale factors or offsets, you will need to read those in here
    ! and make some vars to store that info.
end do

!calculate the number and indices of the pixels to be calculated

allocate(lonvect(xsize))
allocate(latvect(ysize))

! Figure out the size of the area to be simulated

call check_nc(nf90_inq_varid(metfid,'lon',varid))
call check_nc(nf90_get_var(metfid,varid,lonvect))
call check_nc(nf90_inq_varid(metfid,'lat',varid))
call check_nc(nf90_get_var(metfid,varid,latvect))

call check_nc(nf90_inq_varid(metfid,'lat',varid))
call check_nc(nf90_get_var(metfid,varid,latvect))

call check_nc(nf90_inq_varid(metfid,'time',varid))
call check_nc(nf90_get_var(metfid,varid,timestart,start=(/1,1,1/),count=(/1,1,1/)))
call check_nc(nf90_get_att(metfid,varid,'units',unitchar))

write(*,*)'The metdata in your netcdf file says it starts:',timestart,'for units:',unitchar,&
         &'. The expected units are: day as %Y%m%d.%f.'

yearmetst = int((timestart(1,1,1) - 101.)/10000.) ! This is based on the MET file starting on the first day of the first month at 0:00!

! Based on the bounds, we make vectors of the cells to be run:

pos = minloc(abs(lonvect - bounds(1)))
xpos(1) = pos(1)

pos = minloc(abs(lonvect - bounds(2)))
xpos(2) = pos(1)

pos = minloc(abs(latvect - bounds(3)))
ypos(1) = pos(1)

pos = minloc(abs(latvect - bounds(4)))
ypos(2) = pos(1)

srtx = minval(xpos)
srty = minval(ypos)

if (lonvect(srtx) < bounds(1) .and. bounds(2) /= bounds(1)) srtx = srtx + 1
 cntx = 1 + abs(maxval(xpos) - srtx)

if (latvect(srty) < bounds(3) .and. bounds(4) /= bounds(3)) srty = srty + 1
 cnty = 1 + abs(maxval(ypos) - srty)

end subroutine openmet

! ! ================================================================================

! subroutine readin_met(start_pos,dlat,dlon)
!
! ! Read in one year of the MET file at once.
!
! ! J. Melton
! ! Nov 2016
!
! use typeSizes
! use netcdf
! use netcdf_drivers, only : check_nc
! use io_driver, only : vname,niv,metfid,cntx,cnty,srtx,srty,bounds,lonvect,latvect, &
!                       yearmetst
! use ctem_statevars,     only : c_switch
! use ctem_params, only : nlat
!
! implicit none
!
! ! Arguments
! integer, intent(in) :: start_pos            !< first running timestep of the met dataset to read in.
! real, dimension(:), intent(in) :: dlat   !< Latitude of grid cell [degrees]
! real, dimension(:), intent(in) :: dlon   !< Longitude of grid cell (east of Greenwich) [degrees]
!
! ! Pointers
! integer, pointer :: met_ts_sec
!
! ! Local variables
! integer :: tlen                                 !< total timesteps to readin at once
! integer :: i,y,x,j,mo,px,py
! integer :: clim_mo_read                         !< total number of timesteps to read in
! real, dimension(:,:,:), allocatable :: var_in
! real :: yr_tot
! integer :: varid
!
! ! Parameters
! real, parameter :: secinyear = 31536000.    !< Number of seconds in a non-leap year. Needs to be adapted if leap years considered.
!
! ! Point pointers
! met_ts_sec        => c_switch%met_ts_sec
!
! !---
!
! !>
! !! First we need to figure out where we start to read in from and how many
! !! timesteps of the met netcdf file we need to read
!
! tlen = int(secinyear / real(met_ts_sec))
!
!   write(0,*)'Reading in ',tlen,' timesteps of met'
!
!  clim_mo_read = tlen + 1 - start_pos
!
! allocate(var_in(cntx,cnty,clim_mo_read))
!
!   do i = 1,niv
!
!     call check_nc(nf90_inq_varid(metfid,vname(i),varid))
!     call check_nc(nf90_get_var(metfid,varid,var_in,start=[srtx,srty,start_pos],count=[cntx,cnty,tlen]))  ! FLAG check the y,x order here.
!
!     do y = 1,cnty
!       do x = 1,cntx
!         j = x + 1 * (y-1)
!
! !         if (.not. sv(j)%valid_cell) cycle  ! do not get climate for this gridcell, it is not valid
! !
!          px = srtx + x-1
!          py = srty + y-1
!          dlon = lonvect(px)
!          dlat = latvect(py)
!
!         ! Fill in the 6 hourly met arrays
!
! !         ! Set up yearly climate data
! !         select case (i)
! !         case(1)
! !           ibuf(j)%cld = va(i)%scale_factor * real(var_in(x,y,:)) + va(i)%add_offset
! !
! !                   ! Convert cloud cover (in percent) to a fraction
! !           ibuf(j)%cld = 0.01_dp * ibuf(j)%cld
! !
! !         case(2)
! !           ibuf(j)%dtr = va(i)%scale_factor * real(var_in(x,y,:)) + va(i)%add_offset
! !
! !                   ! Correct the input values for potential inconsistencies
! !           ibuf(j)%dtr = max(ibuf(j)%dtr,0._dp)
! !
! !         case(3)
! !           ibuf(j)%pre = va(i)%scale_factor * real(var_in(x,y,:)) + va(i)%add_offset
! !
! !                    !FLAG this is a temporary fix for problems with dataset pre amounts
! !                    !JM 09.05.2011
! !                    yr_tot = sum(ibuf(j)%pre)
! !                    do mo = 1,12
! !                    if (ibuf(j)%pre(mo) > 0.9 * yr_tot) ibuf(j)%pre(mo) = ibuf(j)%pre(mo) * 0.01
! !                    end do
! !
! !         case(4)
! !           ibuf(j)%tmp = va(i)%scale_factor * real(var_in(x,y,:)) + va(i)%add_offset
! !         case(5)
! !           ibuf(j)%wet = va(i)%scale_factor * real(var_in(x,y,:)) + va(i)%add_offset
! !         end select
!
!       end do
!
!     end do
!   end do
!
! deallocate(var_in)
!
! end subroutine readin_met

!--------------------------------------------------------------------------------------------

subroutine read_initialstate(nmtest,nltest,ignd,onetile_perPFT) ! FLAG these are just placeholders now. Figure out later what I really want.

! J. Melton
! Nov 2016

use netcdf
use netcdf_drivers, only : check_nc
use io_driver, only : cntx,cnty,srtx,srty !vname,niv,metfid,bounds,lonvect,latvect, yearmetst,rsid
use ctem_statevars,     only : c_switch,vrot
use class_statevars,    only : class_rot
! use ctem_params,        only : icc,iccp1,nmos,seed,ignd,ilg,icp1,nlat,ican,abszero

implicit none

! arguments:
integer, intent(out) :: ignd
integer, intent(out) :: nltest
integer, intent(inout) :: nmtest
logical, intent(in) :: onetile_perPFT

! pointers:
real, pointer, dimension(:,:,:) :: FCANROT
real, pointer, dimension(:,:)   :: FAREROT
real, pointer, dimension(:,:,:) :: RSMNROT
real, pointer, dimension(:,:,:) :: QA50ROT
real, pointer, dimension(:,:,:) :: VPDAROT
real, pointer, dimension(:,:,:) :: VPDBROT
real, pointer, dimension(:,:,:) :: PSGAROT
real, pointer, dimension(:,:,:) :: PSGBROT
real, pointer, dimension(:,:)   :: DRNROT
real, pointer, dimension(:,:)   :: SDEPROT
real, pointer, dimension(:,:)   :: XSLPROT
real, pointer, dimension(:,:)   :: GRKFROT
real, pointer, dimension(:,:)   :: WFSFROT
real, pointer, dimension(:,:)   :: WFCIROT
integer, pointer, dimension(:,:)   :: MIDROT
real, pointer, dimension(:,:,:) :: SANDROT
real, pointer, dimension(:,:,:) :: CLAYROT
real, pointer, dimension(:,:,:) :: ORGMROT
real, pointer, dimension(:,:,:) :: TBARROT
real, pointer, dimension(:,:,:) :: THLQROT
real, pointer, dimension(:,:,:) :: THICROT
real, pointer, dimension(:,:)   :: TCANROT
real, pointer, dimension(:,:)   :: TSNOROT
real, pointer, dimension(:,:)   :: TPNDROT
real, pointer, dimension(:,:)   :: ZPNDROT
real, pointer, dimension(:,:)   :: RCANROT
real, pointer, dimension(:,:)   :: SCANROT
real, pointer, dimension(:,:)   :: SNOROT
real, pointer, dimension(:,:)   :: ALBSROT
real, pointer, dimension(:,:)   :: RHOSROT
real, pointer, dimension(:,:)   :: GROROT

logical, pointer :: dofire
logical, pointer :: compete
logical, pointer :: inibioclim
logical, pointer :: dowetlands
logical, pointer :: start_bare
logical, pointer :: lnduseon
logical, pointer :: obswetf
real, pointer, dimension(:,:,:) :: ailcminrow           !
real, pointer, dimension(:,:,:) :: ailcmaxrow           !
real, pointer, dimension(:,:,:) :: dvdfcanrow           !
real, pointer, dimension(:,:,:) :: gleafmasrow          !
real, pointer, dimension(:,:,:) :: bleafmasrow          !
real, pointer, dimension(:,:,:) :: stemmassrow          !
real, pointer, dimension(:,:,:) :: rootmassrow          !
real, pointer, dimension(:,:,:) :: pstemmassrow         !
real, pointer, dimension(:,:,:) :: pgleafmassrow        !
real, pointer, dimension(:,:) :: twarmm            !< temperature of the warmest month (c)
real, pointer, dimension(:,:) :: tcoldm            !< temperature of the coldest month (c)
real, pointer, dimension(:,:) :: gdd5              !< growing degree days above 5 c
real, pointer, dimension(:,:) :: aridity           !< aridity index, ratio of potential evaporation to precipitation
real, pointer, dimension(:,:) :: srplsmon          !< number of months in a year with surplus water i.e.precipitation more than potential evaporation
real, pointer, dimension(:,:) :: defctmon          !< number of months in a year with water deficit i.e.precipitation less than potential evaporation
real, pointer, dimension(:,:) :: anndefct          !< annual water deficit (mm)
real, pointer, dimension(:,:) :: annsrpls          !< annual water surplus (mm)
real, pointer, dimension(:,:) :: annpcp            !< annual precipitation (mm)
real, pointer, dimension(:,:) :: dry_season_length !< length of dry season (months)
real, pointer, dimension(:,:,:) :: litrmassrow
real, pointer, dimension(:,:,:) :: soilcmasrow
real, pointer, dimension(:,:) :: extnprob
real, pointer, dimension(:,:) :: prbfrhuc
real, pointer, dimension(:,:,:) :: mlightng
integer, pointer, dimension(:,:,:) :: lfstatusrow
integer, pointer, dimension(:,:,:) :: pandaysrow
integer, pointer, dimension(:,:) :: stdaln
real, pointer, dimension(:,:,:) :: slopefrac

! local variables

integer :: i,m,j,strlen
integer :: initid
!real, dimension(ilg,2) :: crop_temp_frac
integer, dimension(cnty,cntx) :: maskofregion

! point pointers:
dofire            => c_switch%dofire
compete           => c_switch%compete
inibioclim        => c_switch%inibioclim
dowetlands        => c_switch%dowetlands
start_bare        => c_switch%start_bare
lnduseon          => c_switch%lnduseon
obswetf           => c_switch%obswetf
ailcminrow        => vrot%ailcmin
ailcmaxrow        => vrot%ailcmax
dvdfcanrow        => vrot%dvdfcan
gleafmasrow       => vrot%gleafmas
bleafmasrow       => vrot%bleafmas
stemmassrow       => vrot%stemmass
rootmassrow       => vrot%rootmass
pstemmassrow      => vrot%pstemmass
pgleafmassrow     => vrot%pgleafmass
twarmm            => vrot%twarmm
tcoldm            => vrot%tcoldm
gdd5              => vrot%gdd5
aridity           => vrot%aridity
srplsmon          => vrot%srplsmon
defctmon          => vrot%defctmon
anndefct          => vrot%anndefct
annsrpls          => vrot%annsrpls
annpcp            => vrot%annpcp
dry_season_length => vrot%dry_season_length
litrmassrow       => vrot%litrmass
soilcmasrow       => vrot%soilcmas
extnprob          => vrot%extnprob
prbfrhuc          => vrot%prbfrhuc
mlightng          => vrot%mlightng
slopefrac         => vrot%slopefrac
stdaln            => vrot%stdaln
lfstatusrow       => vrot%lfstatus
pandaysrow        => vrot%pandays

FCANROT           => class_rot%FCANROT
FAREROT           => class_rot%FAREROT
RSMNROT           => class_rot%RSMNROT
QA50ROT           => class_rot%QA50ROT
VPDAROT           => class_rot%VPDAROT
VPDBROT           => class_rot%VPDBROT
PSGAROT           => class_rot%PSGAROT
PSGBROT           => class_rot%PSGBROT
DRNROT           => class_rot%DRNROT
SDEPROT           => class_rot%SDEPROT
XSLPROT           => class_rot%XSLPROT
GRKFROT           => class_rot%GRKFROT
WFSFROT           => class_rot%WFSFROT
WFCIROT           => class_rot%WFCIROT
MIDROT           => class_rot%MIDROT
SANDROT           => class_rot%SANDROT
CLAYROT           => class_rot%CLAYROT
ORGMROT           => class_rot%ORGMROT
TBARROT           => class_rot%TBARROT
THLQROT           => class_rot%THLQROT
THICROT           => class_rot%THICROT
TCANROT           => class_rot%TCANROT
TSNOROT           => class_rot%TSNOROT
TPNDROT           => class_rot%TPNDROT
ZPNDROT           => class_rot%ZPNDROT
RCANROT           => class_rot%RCANROT
SCANROT           => class_rot%SCANROT
SNOROT           => class_rot%SNOROT
ALBSROT           => class_rot%ALBSROT
RHOSROT           => class_rot%RHOSROT
GROROT           => class_rot%GROROT

! ----------------------------

!open initial conditions file

call check_nc(nf90_open(trim(init_file),nf90_nowrite,initid))

!retrieve the number of soil layers

call check_nc(nf90_inq_dimid(initid,'layer',dimid))
call check_nc(nf90_inquire_dimension(initid,dimid,len=ignd))

!retrieve the number of PFTs ! FUTURE?

! CLASS PFTs
!call check_nc(nf90_inq_dimid(initid,'ic',dimid))
!call check_nc(nf90_inquire_dimension(initid,dimid,len=))

! CTEM PFTs
!call check_nc(nf90_inq_dimid(initid,'icc',dimid))
!call check_nc(nf90_inquire_dimension(initid,dimid,len=))


!Figure out how big our nlat vector is by the number of valid cells in the
!region we are running. Use the Mask, it is -1 for land. Also get the nmtest value
!for each of the land cells.

call check_nc(nf90_inq_varid(initid,'Mask',varid))
call check_nc(nf90_get_var(initid,varid,maskofregion,start=[srty,srtx],count=[cnty,cntx]))

call check_nc(nf90_inq_varid(initid,'nmtest',varid))

nltest=0
do x = 1, cntx
    do y = 1, cnty
        if (maskofregion(y,x) == -1) then ! land cell
            nltest = nltest + 1
            call check_nc(nf90_get_var(initid,varid,nmtest,start=[srty,srtx],count=[cnty,cntx]))

        end if
    end do
end do

!  NLTEST,NMTEST
!

! Allocate arrays to proper dimensions

!
! Read in the first values

call check_nc(nf90_inq_varid(initid,'lat',varid))
call check_nc(nf90_get_var(initid,varid,DLATROW,start=[srty,srtx],count=[cnty,cntx]))

call check_nc(nf90_inq_varid(initid,'lon',varid))
call check_nc(nf90_get_var(initid,varid,DLONROW,start=[srty,srtx],count=[cnty,cntx]))

!       JLAT=NINT(DLATROW(1))
!       RADJROW(1)=DLATROW(1)*PI/180.
!       DLONROW(1)=DEGLON
!       Z0ORROW(1)=0.0
!       GGEOROW(1)=0.0
! !     GGEOROW(1)=-0.035


call check_nc(nf90_inq_varid(initid,'ZRFM',varid))
call check_nc(nf90_get_var(initid,varid,ZRFMROW,start=[srty,srtx],count=[cnty,cntx]))

call check_nc(nf90_inq_varid(initid,'ZRFH',varid))
call check_nc(nf90_get_var(initid,varid,ZRFHROW,start=[srty,srtx],count=[cnty,cntx]))

call check_nc(nf90_inq_varid(initid,'ZBLD',varid))
call check_nc(nf90_get_var(initid,varid,ZBLDROW,start=[srty,srtx],count=[cnty,cntx]))

call check_nc(nf90_inq_varid(initid,'ZBLD',varid))
call check_nc(nf90_get_var(initid,varid,ZBLDROW,start=[srty,srtx],count=[cnty,cntx]))

call check_nc(nf90_inq_varid(initid,'GC',varid))
call check_nc(nf90_get_var(initid,varid,GCROW,start=[srty,srtx],count=[cnty,cntx]))

! Now get the per tile variables:

call check_nc(nf90_inq_varid(initid,'FCAN',varid))
call check_nc(nf90_get_var(initid,varid,FCANROT,start=[1,srty,srtx],count=[icp1,cnty,cntx]))

call check_nc(nf90_inq_varid(initid,'LNZ0',varid))
call check_nc(nf90_get_var(initid,varid,LNZ0ROT,start=[1,srty,srtx],count=[icp1,cnty,cntx]))

call check_nc(nf90_inq_varid(initid,'ALVC',varid))
call check_nc(nf90_get_var(initid,varid,ALVCROT,start=[1,srty,srtx],count=[icp1,cnty,cntx]))

call check_nc(nf90_inq_varid(initid,'PAMN',varid))
call check_nc(nf90_get_var(initid,varid,PAMNROT,start=[1,srty,srtx],count=[ic,cnty,cntx]))

call check_nc(nf90_inq_varid(initid,'PAMX',varid))
call check_nc(nf90_get_var(initid,varid,PAMX,start=[1,srty,srtx],count=[ic,cnty,cntx]))

call check_nc(nf90_inq_varid(initid,'CMAS',varid))
call check_nc(nf90_get_var(initid,varid,CMASROT,start=[1,srty,srtx],count=[ic,cnty,cntx]))

call check_nc(nf90_inq_varid(initid,'ALIC',varid))
call check_nc(nf90_get_var(initid,varid,ALICROT,start=[1,srty,srtx],count=[icp1,cnty,cntx]))

call check_nc(nf90_inq_varid(initid,'ROOT',varid))
call check_nc(nf90_get_var(initid,varid,ROOTROT,start=[1,srty,srtx],count=[ic,cnty,cntx]))

call check_nc(nf90_inq_varid(initid,'DRN',varid))
call check_nc(nf90_get_var(initid,varid,DRNROT,start=[srty,srtx],count=[icnty,cntx]))

call check_nc(nf90_inq_varid(initid,'RSMN',varid))
call check_nc(nf90_get_var(initid,varid,RSMNROT,start=[1,srty,srtx],count=[ic,cnty,cntx]))

call check_nc(nf90_inq_varid(initid,'QA50',varid))
call check_nc(nf90_get_var(initid,varid,QA50,start=[1,srty,srtx],count=[ic,cnty,cntx]))

call check_nc(nf90_inq_varid(initid,'VPDA',varid))
call check_nc(nf90_get_var(initid,varid,VPDAROT,start=[1,srty,srtx],count=[ic,cnty,cntx]))

call check_nc(nf90_inq_varid(initid,'VPDB',varid))
call check_nc(nf90_get_var(initid,varid,VPDBROT,start=[1,srty,srtx],count=[ic,cnty,cntx]))

call check_nc(nf90_inq_varid(initid,'PSGA',varid))
call check_nc(nf90_get_var(initid,varid,PSGAROT,start=[1,srty,srtx],count=[ic,cnty,cntx]))

call check_nc(nf90_inq_varid(initid,'PSGB',varid))
call check_nc(nf90_get_var(initid,varid,PSGBROT,start=[1,srty,srtx],count=[ic,cnty,cntx]))

call check_nc(nf90_inq_varid(initid,'SDEP',varid))
call check_nc(nf90_get_var(initid,varid,SDEPROT,start=[srty,srtx],count=[icnty,cntx]))

call check_nc(nf90_inq_varid(initid,'FARE',varid))
call check_nc(nf90_get_var(initid,varid,FAREROT,start=[srty,srtx],count=[icnty,cntx]))

!           ! Error check:
!           if (FAREROT(I,M) .gt. 1.0) then
!            write(*,*)'FAREROT > 1',FAREROT(I,M)
!            call XIT('runclass36ctem', -2)
!           end if

call check_nc(nf90_inq_varid(initid,'XSLP',varid))
call check_nc(nf90_get_var(initid,varid,XSLPROT,start=[srty,srtx],count=[icnty,cntx]))

call check_nc(nf90_inq_varid(initid,'GRKF',varid))
call check_nc(nf90_get_var(initid,varid,GRKFROT,start=[srty,srtx],count=[icnty,cntx]))

call check_nc(nf90_inq_varid(initid,'WFSF',varid))
call check_nc(nf90_get_var(initid,varid,WFSFROT,start=[srty,srtx],count=[icnty,cntx]))

call check_nc(nf90_inq_varid(initid,'WFCI',varid))
call check_nc(nf90_get_var(initid,varid,WFCIROT,start=[srty,srtx],count=[icnty,cntx]))

call check_nc(nf90_inq_varid(initid,'MID',varid))
call check_nc(nf90_get_var(initid,varid,MIDROT,start=[srty,srtx],count=[icnty,cntx]))

call check_nc(nf90_inq_varid(initid,'ZBOT',varid)) !per layer but not per tile (fix?)
call check_nc(nf90_get_var(initid,varid,ZBOT,start=[1,srty,srtx],count=[ignd,icnty,cntx]))

call check_nc(nf90_inq_varid(initid,'DELZ',varid)) !per layer but not per tile (fix?)
call check_nc(nf90_get_var(initid,varid,DELZ,start=[1,srty,srtx],count=[ignd,icnty,cntx]))

call check_nc(nf90_inq_varid(initid,'SAND',varid))
call check_nc(nf90_get_var(initid,varid,SANDROT,start=[1,srty,srtx],count=[ignd,icnty,cntx]))

call check_nc(nf90_inq_varid(initid,'CLAY',varid))
call check_nc(nf90_get_var(initid,varid,CLAYROT,start=[1,srty,srtx],count=[ignd,icnty,cntx]))

call check_nc(nf90_inq_varid(initid,'ORGM',varid))
call check_nc(nf90_get_var(initid,varid,ORGMROT,start=[1,srty,srtx],count=[ignd,icnty,cntx]))

call check_nc(nf90_inq_varid(initid,'TBAR',varid))
call check_nc(nf90_get_var(initid,varid,TBARROT,start=[1,srty,srtx],count=[ignd,icnty,cntx]))

call check_nc(nf90_inq_varid(initid,'THLQ',varid))
call check_nc(nf90_get_var(initid,varid,THLQROT,start=[1,srty,srtx],count=[ignd,icnty,cntx]))

call check_nc(nf90_inq_varid(initid,'THIC',varid))
call check_nc(nf90_get_var(initid,varid,THICROT,start=[1,srty,srtx],count=[ignd,icnty,cntx]))

call check_nc(nf90_inq_varid(initid,'TCAN',varid))
call check_nc(nf90_get_var(initid,varid,TCANROT,start=[srty,srtx],count=[icnty,cntx]))

call check_nc(nf90_inq_varid(initid,'TSNO',varid))
call check_nc(nf90_get_var(initid,varid,TSNOROT,start=[srty,srtx],count=[icnty,cntx]))

call check_nc(nf90_inq_varid(initid,'TPND',varid))
call check_nc(nf90_get_var(initid,varid,TPNDROT,start=[srty,srtx],count=[icnty,cntx]))

call check_nc(nf90_inq_varid(initid,'ZPND',varid))
call check_nc(nf90_get_var(initid,varid,ZPNDROT,start=[srty,srtx],count=[icnty,cntx]))

call check_nc(nf90_inq_varid(initid,'RCAN',varid))
call check_nc(nf90_get_var(initid,varid,RCANROT,start=[srty,srtx],count=[icnty,cntx]))

call check_nc(nf90_inq_varid(initid,'SCAN',varid))
call check_nc(nf90_get_var(initid,varid,SCANROT,start=[srty,srtx],count=[icnty,cntx]))

call check_nc(nf90_inq_varid(initid,'SNO',varid))
call check_nc(nf90_get_var(initid,varid,SNOROT,start=[srty,srtx],count=[icnty,cntx]))

call check_nc(nf90_inq_varid(initid,'ALBS',varid))
call check_nc(nf90_get_var(initid,varid,ALBSROT,start=[srty,srtx],count=[icnty,cntx]))

call check_nc(nf90_inq_varid(initid,'RHOS',varid))
call check_nc(nf90_get_var(initid,varid,RHOSROT,start=[srty,srtx],count=[icnty,cntx]))

call check_nc(nf90_inq_varid(initid,'GRO',varid))
call check_nc(nf90_get_var(initid,varid,GROROT,start=[srty,srtx],count=[icnty,cntx]))

!>The following three variables are needed to run CTEM. 1) min & 2) max leaf area index are needed to break
!>class lai into dcd and evg for trees (for crops and grasses it doesn't matter much). 3) dvdfcanrow is
!>needed to divide needle & broad leaf into dcd and evg, and crops & grasses into c3 and c4 fractions.

call check_nc(nf90_inq_varid(initid,'ailcmin',varid))
call check_nc(nf90_get_var(initid,varid,ailcminrow,start=[1,srty,srtx],count=[icc,icnty,cntx]))

call check_nc(nf90_inq_varid(initid,'ailcmax',varid))
call check_nc(nf90_get_var(initid,varid,ailcmaxrow,start=[1,srty,srtx],count=[icc,icnty,cntx]))

call check_nc(nf90_inq_varid(initid,'dvdfcan',varid))
call check_nc(nf90_get_var(initid,varid,dvdfcanrow,start=[1,srty,srtx],count=[icc,icnty,cntx]))

!>
!>Rest of the initialization variables are needed to run CTEM but if starting from bare ground initialize all
!>live and dead c pools from zero. suitable values of extnprobgrd and prbfrhucgrd would still be required. set
!>stdaln to 1 for operation in non-gcm stand alone mode, in the CTEM initialization file.
!>
call check_nc(nf90_inq_varid(initid,'gleafmas',varid))
call check_nc(nf90_get_var(initid,varid,gleafmasrow,start=[1,srty,srtx],count=[icc,icnty,cntx]))

call check_nc(nf90_inq_varid(initid,'bleafmas',varid))
call check_nc(nf90_get_var(initid,varid,bleafmasrow,start=[1,srty,srtx],count=[icc,icnty,cntx]))

call check_nc(nf90_inq_varid(initid,'stemmass',varid))
call check_nc(nf90_get_var(initid,varid,stemmassrow,start=[1,srty,srtx],count=[icc,icnty,cntx]))

call check_nc(nf90_inq_varid(initid,'rootmass',varid))
call check_nc(nf90_get_var(initid,varid,rootmassrow,start=[1,srty,srtx],count=[icc,icnty,cntx]))

!>
!>If fire and competition are on, save the stemmass and rootmass for use in burntobare subroutine on the first timestep.
if (dofire .and. compete) then
            do j =1,icc
            pstemmassrow(i,m,j)=stemmassrow(i,m,j)
            pgleafmassrow(i,m,j)=rootmassrow(i,m,j)
            end do
end if
!
call check_nc(nf90_inq_varid(initid,'litrmass',varid))
call check_nc(nf90_get_var(initid,varid,litrmassrow,start=[1,srty,srtx],count=[iccp1,icnty,cntx])) !FLAG per layer?

call check_nc(nf90_inq_varid(initid,'soilcmas',varid))
call check_nc(nf90_get_var(initid,varid,soilcmasrow,start=[1,srty,srtx],count=[iccp1,icnty,cntx])) !FLAG per layer?

call check_nc(nf90_inq_varid(initid,'lfstatus',varid))
call check_nc(nf90_get_var(initid,varid,lfstatusrow,start=[1,srty,srtx],count=[icc,icnty,cntx]))

call check_nc(nf90_inq_varid(initid,'pandays',varid))
call check_nc(nf90_get_var(initid,varid,pandaysrow,start=[1,srty,srtx],count=[icc,icnty,cntx]))

!
! NOT per tile
call check_nc(nf90_inq_varid(initid,'mlightng',varid))
call check_nc(nf90_get_var(initid,varid,mlightng,start=[1,srty,srtx],count=[12,icnty,cntx]))

call check_nc(nf90_inq_varid(initid,'extnprob',varid))
call check_nc(nf90_get_var(initid,varid,extnprob,start=[srty,srtx],count=[icnty,cntx]))

call check_nc(nf90_inq_varid(initid,'prbfrhuc',varid))
call check_nc(nf90_get_var(initid,varid,prbfrhuc,start=[srty,srtx],count=[icnty,cntx]))

call check_nc(nf90_inq_varid(initid,'stdaln',varid)) !FLAG - get rid of this? I don't see the use stdaln
call check_nc(nf90_get_var(initid,varid,stdaln,start=[srty,srtx],count=[icnty,cntx]))

if (compete .and. inibioclim) then  !read in the bioclimatic parameters

    call check_nc(nf90_inq_varid(initid,'twarmm',varid))
    call check_nc(nf90_get_var(initid,varid,twarmm,start=[srty,srtx],count=[icnty,cntx]))

    call check_nc(nf90_inq_varid(initid,'tcoldm',varid))
    call check_nc(nf90_get_var(initid,varid,tcoldm,start=[srty,srtx],count=[icnty,cntx]))

    call check_nc(nf90_inq_varid(initid,'gdd5',varid))
    call check_nc(nf90_get_var(initid,varid,gdd5,start=[srty,srtx],count=[icnty,cntx]))

    call check_nc(nf90_inq_varid(initid,'aridity',varid))
    call check_nc(nf90_get_var(initid,varid,aridity,start=[srty,srtx],count=[icnty,cntx]))

    call check_nc(nf90_inq_varid(initid,'srplsmon',varid))
    call check_nc(nf90_get_var(initid,varid,srplsmon,start=[srty,srtx],count=[icnty,cntx]))

    call check_nc(nf90_inq_varid(initid,'defctmon',varid))
    call check_nc(nf90_get_var(initid,varid,defctmon,start=[srty,srtx],count=[icnty,cntx]))

    call check_nc(nf90_inq_varid(initid,'anndefct',varid))
    call check_nc(nf90_get_var(initid,varid,anndefct,start=[srty,srtx],count=[icnty,cntx]))

    call check_nc(nf90_inq_varid(initid,'annsrpls',varid))
    call check_nc(nf90_get_var(initid,varid,annsrpls,start=[srty,srtx],count=[icnty,cntx]))

    call check_nc(nf90_inq_varid(initid,'annpcp',varid))
    call check_nc(nf90_get_var(initid,varid,annpcp,start=[srty,srtx],count=[icnty,cntx]))

    call check_nc(nf90_inq_varid(initid,'dry_season_length',varid))
    call check_nc(nf90_get_var(initid,varid,dry_season_length,start=[srty,srtx],count=[icnty,cntx]))

else if (compete .and. .not. inibioclim) then ! set them to zero
    twarmm(i,1)=0.0
    tcoldm(i,1)=0.0
    gdd5(i,1)=0.0
    aridity(i,1)=0.0
    srplsmon(i,1)=0.0
    defctmon(i,1)=0.0
    anndefct(i,1)=0.0
    annsrpls(i,1)=0.0
    annpcp(i,1)=0.0
    dry_season_length(i,1) = 0.0
endif
!
! !>Take the first tile value now and put it over the other tiles !FLAG
!         if (nmtest > 1) then
!             do m = 2,nmtest
!                 twarmm(i,m)=twarmm(i,1)
!                 tcoldm(i,m)=tcoldm(i,1)
!                 gdd5(i,m)=gdd5(i,1)
!                 aridity(i,m)=aridity(i,1)
!                 srplsmon(i,m)=srplsmon(i,1)
!                 defctmon(i,m)=defctmon(i,1)
!                 anndefct(i,m)=anndefct(i,1)
!                 annsrpls(i,m)=annsrpls(i,1)
!                 annpcp(i,m)=annpcp(i,1)
!                 dry_season_length(i,m) =dry_season_length(i,1)
!                 mlightng(i,m,:) = mlightng(i,1,:)
!                 extnprob(i,m) = extnprob(i,1)
!                 prbfrhuc(i,m) = prbfrhuc(i,1)
!                 stdaln(i,m) = stdaln(i,1)
!             end do
!         end if
if (dowetlands) then !if true then read wetland fractions into the first tile position

    call check_nc(nf90_inq_varid(initid,'slopefrac',varid))
    call check_nc(nf90_get_var(initid,varid,slopefrac,start=[1,srty,srtx],count=[8,icnty,cntx]))

    if (nmtest > 1) then ! if more tiles then just put the first tile across the rest FLAG
    do m = 2,nmtest
        slopefrac(i,m,:) = slopefrac(i,1,:)
    end do
    end if
endif
! 71    continue
!
! close(11)
!
!
! !>Check that a competition or luc run has the correct number of mosaics. if it is not a start_bare run, then nmtest should equal nmos
!       if (onetile_perPFT .and. (compete .or. lnduseon) .and. .not. start_bare) then
!         if (nmtest .ne. nmos) then
!            write(6,*)'compete or luc runs that do not start from bare'
!            write(6,*)'ground need the number of mosaics to equal icc+1'
!            write(6,*)'nmtest = ',nmtest,' nmos = ',nmos
!             call xit('runclass36ctem', -2)
!         endif
!       endif
! !>
! !>if this run uses the competition or lnduseon parameterization and starts from bare ground, set up the model state here. this
! !>overwrites what was read in from the .ini and .ctm files. for composite runs (the composite set up is after this one for mosaics)
!       if ((compete .or. lnduseon) .and. start_bare) then
!
!        if (onetile_perPFT) then
! !>
! !!store the read-in crop fractions as we keep them even when we start bare.
! !!FLAG: this is setup assuming that crops are in mosaics 6 and 7. JM Apr 9 2014.
!          do i=1,nltest
!           crop_temp_frac(i,1)=FAREROT(i,6)
!           crop_temp_frac(i,2)=FAREROT(i,7)
!          end do
!
! !>check the number of mosaics that came from the .ini file
!         if (nmtest .ne. nmos) then
!
! !>we need to transfer some initial parameterization info to all mosaics, so set all values to that of the first mosaic.
!          do i=1,nltest
!           do m=nmtest+1,nmos
!
!            do j=1,ican
!              RSMNROT(i,m,j)=RSMNROT(i,1,j)
!              QA50ROT(i,m,j)=QA50ROT(i,1,j)
!              VPDAROT(i,m,j)=VPDAROT(i,1,j)
!              VPDBROT(i,m,j)=VPDBROT(i,1,j)
!              PSGAROT(i,m,j)=PSGAROT(i,1,j)
!              PSGBROT(i,m,j)=PSGBROT(i,1,j)
!            enddo
!
!            DRNROT(i,m)=DRNROT(i,1)
!            SDEPROT(i,m)=SDEPROT(i,1)
!            FAREROT(i,m)=FAREROT(i,1)
!            XSLPROT(i,m)=XSLPROT(i,1)
!            GRKFROT(i,m)=GRKFROT(i,1)
!            WFSFROT(i,m)=WFSFROT(i,1)
!            WFCIROT(i,m)=WFCIROT(i,1)
!            MIDROT(i,m)=MIDROT(i,1)
!
!            do j=1,3
!             SANDROT(i,m,j)=SANDROT(i,1,j)
!             CLAYROT(i,m,j)=CLAYROT(i,1,j)
!             ORGMROT(i,m,j)=ORGMROT(i,1,j)
!             TBARROT(i,m,j)=TBARROT(i,1,j)
!             THLQROT(i,m,j)=THLQROT(i,1,j)
!             THICROT(i,m,j)=THICROT(i,1,j)
!            enddo
!
!            TCANROT(i,m)=TCANROT(i,1)
!            TSNOROT(i,m)=TSNOROT(i,1)
!            TPNDROT(i,m)=TPNDROT(i,1)
!            ZPNDROT(i,m)=ZPNDROT(i,1)
!            RCANROT(i,m)=RCANROT(i,1)
!            SCANROT(i,m)=SCANROT(i,1)
!            SNOROT(i,m)=SNOROT(i,1)
!            ALBSROT(i,m)=ALBSROT(i,1)
!            RHOSROT(i,m)=RHOSROT(i,1)
!            GROROT(i,m)=GROROT(i,1)
!            do j=1,icc
!              lfstatusrow(i,m,j) = 4
!            enddo !j
!
!           enddo !m
!          enddo !i
!
! !>set the number of mosaics to icc+1
!         nmtest=nmos
!
!         endif  !>if (nmtest .ne. nmos)
!
! !>set the initial conditions for the pfts
! ! (bah, this is such an inelegant way to do this, but oh well...)
!
! !>initalize to zero
!         FCANROT=0.0
!         dvdfcanrow=0.0
!         FAREROT=0.0
!
!         do i=1,nltest
!          do m=1,nmtest
!
! !>set the seed amount for each pft in its mosaic
!           if (compete .or. lnduseon) then
!             if (m .lt. icc+1) then
!              FAREROT(i,m)=seed
!             else
!              FAREROT(i,m)=1.0 - (real(icc) * seed)
!             endif
!           endif
!
!           do j = 1,icc
!             ailcminrow(i,m,j)=0.0
!             ailcmaxrow(i,m,j)=0.0
!             gleafmasrow(i,m,j)=0.0
!             bleafmasrow(i,m,j)=0.0
!             stemmassrow(i,m,j)=0.0
!             rootmassrow(i,m,j)=0.0
!             lfstatusrow(i,m,j)=4
!             pandaysrow(i,m,j)=0
!           enddo
!
!           lfstatusrow(i,m,1)=2
!
!           do j = 1,iccp1
!             litrmassrow(i,m,j)=0.
!             soilcmasrow(i,m,j)=0.
!           enddo
!
! !>initial conditions always required
!           dvdfcanrow(i,m,1)=1.0  !ndl
!           dvdfcanrow(i,m,3)=1.0  !bdl
!           dvdfcanrow(i,m,6)=1.0  !crop
!           dvdfcanrow(i,m,8)=1.0  !grasses
!
! !>then adjusted below for the actual mosaic makeup
!           if (m .le. 2) then                     !ndl
!            FCANROT(i,m,1)=1.0
!            if (m .eq. 2) then
!              dvdfcanrow(i,m,1)=0.0
!              dvdfcanrow(i,m,2)=1.0
!            endif
!           elseif (m .ge. 3 .and. m .le. 5) then  !bdl
!            FCANROT(i,m,2)=1.0
!            if (m .eq. 4) then
!              dvdfcanrow(i,m,3)=0.0
!              dvdfcanrow(i,m,4)=1.0
!            endif
!            if (m .eq. 5) then
!              dvdfcanrow(i,m,3)=0.0
!              dvdfcanrow(i,m,5)=1.0
!            endif
!           elseif (m .eq. 6 .or. m .eq. 7) then  !crop
!            FCANROT(i,m,3)=1.0
!            if (m .eq. 7) then
!              dvdfcanrow(i,m,6)=0.0
!              dvdfcanrow(i,m,7)=1.0
!            endif
!           elseif (m .eq. 8 .or. m .eq. 9) then  !grasses
!            FCANROT(i,m,4)=1.0
!            if (m .eq. 9) then
!              dvdfcanrow(i,m,8)=0.0
!              dvdfcanrow(i,m,9)=1.0
!            endif
!           else                                  !bare/urban?
!            FCANROT(i,m,5)=1.0
!            endif !mosaic adjustments
!          enddo  !m
!         enddo  !i
!
!
!          do i=1,nltest
!           FAREROT(i,6)=crop_temp_frac(i,1)
!           FAREROT(i,7)=crop_temp_frac(i,2)
!          end do
!
!       else if (.not. onetile_perPFT) then
! !>set up for composite runs when start_bare is on and compete or landuseon
!
! !>store the read-in crop fractions as we keep them even when we start bare.
! !!FLAG: this is setup assuming that crops are in pft number 6 and 7.
! !!and the first tile contains the information for the grid cell (assumes we have crops in
! !!every tile too! JM Apr 9 2014.
!          do i=1,nltest
!           crop_temp_frac(i,1)=FCANROT(i,1,3)*dvdfcanrow(i,1,6)
!           crop_temp_frac(i,2)=FCANROT(i,1,3)*dvdfcanrow(i,1,7)
!          end do
! !>
! !>initalize to zero, these will be filled in by the luc or competition subroutines.
!        FCANROT=0.0
!        dvdfcanrow=0.0
!
!        ! Added this as start_bare runs were not properly assigning
!        ! a TCAN on the very first day since the FCANROT was 0. JM Jan 14 2014.
! !        do i=1,nltest
! !         do m = 1,nmtest
! !          do j=1,icp1
! !            if (j .lt. icp1) then
! !             FCANROT(i,m,j)=seed
! !            else
! !             FCANROT(i,m,j)=1.0 - (real(ican) * seed)
! !            endif
! !          end do
! !         end do
! !        end do
!
!        do i=1,nltest
!          do m = 1,nmtest
!
!     !>initial conditions always required
!             dvdfcanrow(i,m,1)=1.0  !ndl
!             dvdfcanrow(i,m,3)=1.0  !bdl
!             dvdfcanrow(i,m,6)=1.0  !crop
!             dvdfcanrow(i,m,8)=1.0  !grasses
!
!             do j = 1,icc
!             ailcminrow(i,m,j)=0.0
!             ailcmaxrow(i,m,j)=0.0
!             gleafmasrow(i,m,j)=0.0
!             bleafmasrow(i,m,j)=0.0
!             stemmassrow(i,m,j)=0.0
!             rootmassrow(i,m,j)=0.0
!             lfstatusrow(i,m,j)=4
!             pandaysrow(i,m,j)=0
!             enddo
!
!             lfstatusrow(i,m,1)=2
!
!             do j = 1,iccp1
!             litrmassrow(i,m,j)=0.0
!             soilcmasrow(i,m,j)=0.0
!             enddo
!          end do ! nmtest
!        enddo !nltest
!
!          do i=1,nltest
!             do m = 1,nmtest
!                 FCANROT(i,m,3) = crop_temp_frac(i,1) + crop_temp_frac(i,2)
!                 if (FCANROT(i,m,3) .gt. abszero) then
!                 dvdfcanrow(i,m,6) = crop_temp_frac(i,1) / FCANROT(i,m,3)
!                 dvdfcanrow(i,m,7) = crop_temp_frac(i,2) / FCANROT(i,m,3)
!                 else
!                 dvdfcanrow(i,m,6) = 1.0
!                 dvdfcanrow(i,m,7) = 0.0
!                 end if
!             end do !nmtest
!          end do !nltest
!
!       end if ! mosaic / composite
!       end if !if (compete/landuseon .and. start_bare)


end subroutine read_initialstate

end module input_dataset_drivers
