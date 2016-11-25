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

!public :: soildriver
!public :: slopedriver
!public :: openGHGfile
!public :: readGHGfile


contains

!---

subroutine openmet(path)

!> The openmet subroutine opens the meteorology netcdf file and prepares it to be
!! read in. The subroutine expects the latitude and longitude dimensions to be called
!! 'lat' and 'lon' in the netcdf, respectively. Additionally the names of the meteorology
!! variables must be the same as those defined in io_driver/vname. The dimensions of this meteorology
!! file is used to determine the number of grid cells that we will running for the
!! model simulation. It is vital that all model input files are at the same resolution
!! and grid as the meteorology file since we use it to determine the size of the vectors
!! latvect and lonvect.

! J. Melton
! Nov 2016

use netcdf
use netcdf_drivers, only : check_nc
use io_driver, only : vname,niv,metfid,cntx,cnty,srtx,srty,bounds,lonvect,latvect

implicit none

character(180), intent(in) :: path
character(180) :: climatefile

integer :: i

integer :: xsize
integer :: ysize
integer :: timelen

integer :: dimid
integer :: varid

integer, dimension(1) :: pos
integer, dimension(2) :: xpos,ypos

!---

 climatefile = trim(path)

!open climate file

call check_nc(nf90_open(climatefile,nf90_nowrite,metfid))

!retrieve dimensions

call check_nc(nf90_inq_dimid(metfid,'lon',dimid))
call check_nc(nf90_inquire_dimension(metfid,dimid,len=xsize))
call check_nc(nf90_inq_dimid(metfid,'lat',dimid))
call check_nc(nf90_inquire_dimension(metfid,dimid,len=ysize))
call check_nc(nf90_inq_dimid(metfid,'time',dimid))
call check_nc(nf90_inquire_dimension(metfid,dimid,len=timelen))

!retrieve scale factors and offsets
write(*,*)xsize,ysize,timelen
do i = 1,niv
    write(*,*)vname(i)
  call check_nc(nf90_inq_varid(metfid,vname(i),varid))

!  status = nf90_get_att(metfid,varid,'scale_factor',va(i)%scale_factor)
!  if (status /= nf90_noerr) call handle_err(status)

!  status = nf90_get_att(metfid,varid,'add_offset',va(i)%add_offset)
!  if (status /= nf90_noerr) call handle_err(status)

end do

!calculate the number and indices of the pixels to be calculated

allocate(lonvect(xsize))
allocate(latvect(ysize))

! Figure out the size of the area to be simulated

call check_nc(nf90_inq_varid(metfid,'lon',varid))
call check_nc(nf90_get_var(metfid,varid,lonvect))
call check_nc(nf90_inq_varid(metfid,'lat',varid))
call check_nc(nf90_get_var(metfid,varid,latvect))

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

! subroutine readin_met(start_pos,tlen)
!
! ! Read in one year of the MET file at once.
!
!
! use typeSizes
! use netcdf
! use netcdf_error
! !use iovariables,         only : ibuf,timebuflen,latvect,lonvect,cntx,cnty,srtx,srty,niv,vname,metfid,va,climatemonths,inputclimlen
! !use statevars,           only : sv
! !use arveparams, only : dp
!
! implicit none
!
! integer, intent(in) :: start_pos                !first running month of the climate dataset to read in.
! integer, intent(in) :: tlen                     !total months to readin
!
! integer :: i,y,x,j,mo,px,py
! integer :: clim_mo_read                         !total number of months to read in
! integer(2), dimension(:,:,:), allocatable :: var_in
! real(dp) :: yr_tot
! integer :: varid
!
! !---
!
!   write(0,*)'Reading in ',tlen,' months of climate'
!
! clim_mo_read = tlen + 1 - start_pos
!
! do i=1,inputclimlen
!  !new climate variables
!   allocate(ibuf(i)%cld(clim_mo_read))
!   allocate(ibuf(i)%dtr(clim_mo_read))
!   allocate(ibuf(i)%pre(clim_mo_read))
!   allocate(ibuf(i)%tmp(clim_mo_read))
!   allocate(ibuf(i)%wet(clim_mo_read))
! end do
!
! allocate(var_in(cntx,cnty,clim_mo_read))
!
!   do i = 1,niv
!
!     status = nf90_inq_varid(metfid,vname(i),varid)
!     if (status /= nf90_noerr) call handle_err(status)
!
!     status = nf90_get_var(metfid,varid,var_in,start=[srtx,srty,start_pos],count=[cntx,cnty,tlen])
!     if (status /= nf90_noerr) call handle_err(status)
!
!     do y = 1,cnty
!       do x = 1,cntx
!         j = x + cntx * (y-1)
!
!         if (.not. sv(j)%valid_cell) cycle  ! do not get climate for this gridcell, it is not valid
!
!         px = srtx + x-1
!         py = srty + y-1
!         sv(j)%lon = lonvect(px)
!         sv(j)%lat = latvect(py)
!
!         ! Set up yearly climate data
!         select case (i)
!         case(1)
!           ibuf(j)%cld = va(i)%scale_factor * real(var_in(x,y,:)) + va(i)%add_offset
!
!                   ! Convert cloud cover (in percent) to a fraction
!           ibuf(j)%cld = 0.01_dp * ibuf(j)%cld
!
!         case(2)
!           ibuf(j)%dtr = va(i)%scale_factor * real(var_in(x,y,:)) + va(i)%add_offset
!
!                   ! Correct the input values for potential inconsistencies
!           ibuf(j)%dtr = max(ibuf(j)%dtr,0._dp)
!
!         case(3)
!           ibuf(j)%pre = va(i)%scale_factor * real(var_in(x,y,:)) + va(i)%add_offset
!
!                    !FLAG this is a temporary fix for problems with dataset pre amounts
!                    !JM 09.05.2011
!                    yr_tot = sum(ibuf(j)%pre)
!                    do mo = 1,12
!                    if (ibuf(j)%pre(mo) > 0.9 * yr_tot) ibuf(j)%pre(mo) = ibuf(j)%pre(mo) * 0.01
!                    end do
!
!         case(4)
!           ibuf(j)%tmp = va(i)%scale_factor * real(var_in(x,y,:)) + va(i)%add_offset
!         case(5)
!           ibuf(j)%wet = va(i)%scale_factor * real(var_in(x,y,:)) + va(i)%add_offset
!         end select
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
! subroutine soildriver(soil_path)
! !new soil driver to take in FAO soil data. Coded Feb 08 JM
!
! use netcdf_error
! use iovariables, only : soil,soilfid,latvect,lonvect,inputsoillen,inputlonlen,inputlatlen,bounds, &
!                         cntx,cnty,srtx,srty
! use netcdf
! use typesizes
! use arveparams, only : dp
!
! implicit none
!
! character(180), intent(in) :: soil_path
! character(100) :: soilfile
!
! integer :: i
! integer :: dimid
! integer :: varid
! integer :: xsize
! integer :: ysize
! integer :: depthint
! integer, dimension(1) :: pos
! integer, dimension(2) :: xpos,ypos
!
! !-------------------------
! !generate file names. Regardless of the path name the input files must always have these names.
! soilfile  = trim(soil_path)
!
! !open soil files
! status = nf90_open(soilfile,nf90_nowrite,soilfid)
! if (status /= nf90_noerr) call handle_err(status)
!
! !find dimensions
! status = nf90_inq_dimid(soilfid,'lon',dimid)
! if (status /= nf90_noerr) call handle_err(status)
!
! status = nf90_inquire_dimension(soilfid,dimid,len=xsize)
! if (status /= nf90_noerr) call handle_err(status)
!
! status = nf90_inq_dimid(soilfid,'lat',dimid)
! if (status /= nf90_noerr) call handle_err(status)
!
! status = nf90_inquire_dimension(soilfid,dimid,len=ysize)
! if (status /= nf90_noerr) call handle_err(status)
!
! status = nf90_inq_dimid(soilfid,'depth',dimid)
! if (status /= nf90_noerr) call handle_err(status)
!
! status = nf90_inquire_dimension(soilfid,dimid,len=depthint)
! if (status /= nf90_noerr) call handle_err(status)
!
! inputsoillen = xsize * ysize
! inputlonlen = xsize
! inputlatlen = ysize
!
! !calculate the number and indices of the pixels to be calculated
! allocate(lonvect(xsize))
! allocate(latvect(ysize))
!
! status = nf90_inq_varid(soilfid,'lon',varid)
! if (status /= nf90_noerr) call handle_err(status)
!
! status = nf90_get_var(soilfid,varid,lonvect)
! if (status /= nf90_noerr) call handle_err(status)
!
! status = nf90_inq_varid(soilfid,'lat',varid)
! if (status /= nf90_noerr) call handle_err(status)
!
! status = nf90_get_var(soilfid,varid,latvect)
! if (status /= nf90_noerr) call handle_err(status)
!
! pos = minloc(abs(lonvect - bounds(1)))
! xpos(1) = pos(1)
!
! pos = minloc(abs(lonvect - bounds(2)))
! xpos(2) = pos(1)
!
! pos = minloc(abs(latvect - bounds(3)))
! ypos(1) = pos(1)
!
! pos = minloc(abs(latvect - bounds(4)))
! ypos(2) = pos(1)
!
! srtx = minval(xpos)
! srty = minval(ypos)
!
! if (lonvect(srtx) < bounds(1) .and. bounds(2) /= bounds(1)) srtx = srtx + 1
! cntx = 1 + abs(maxval(xpos) - srtx)
!
! if (latvect(srty) < bounds(3) .and. bounds(4) /= bounds(3)) srty = srty + 1
! cnty = 1 + abs(maxval(ypos) - srty)
!
! !allocate the vector input buffer for soil data
!
! allocate(soil(inputsoillen))
!
! do i=1,inputsoillen
!   allocate(soil(i)%sdto(depthint))
!   allocate(soil(i)%stpc(depthint))
!   allocate(soil(i)%clpc(depthint))
!   allocate(soil(i)%totc(depthint))
!   allocate(soil(i)%bulk(depthint))
!   allocate(soil(i)%cfrag(depthint))
!   allocate(soil(i)%tawc(depthint))
! end do
!
! !create soil diagnostic file here if needed.
!
! end subroutine soildriver
!
! ! ================================================================================
!
! subroutine openGHGfile(infile)
!
! implicit none
!
! character(80), intent(in) :: infile
! !---------
!
! open(20,file=infile,status='old') ! statue old ensures that it uses an exisiting file
!
! !read header row
! read(20,*) !header
!
! end subroutine openGHGfile
!
! ! ================================================================================
! subroutine readGHGfile
!
! use statevars, only : CO2
!
! implicit none
!
! integer :: year
!
! !------
!
! read(20,*,end=7)year,CO2(1),CO2(2),CO2(3)
!
! return
!
! !this part handles the end-of-file condition
!
! 7 continue
!
! rewind(20)
!
! read(20,*) !header
!
! end subroutine readGHGfile
!
! ! ================================================================================



! !--------------------------------------------

end module input_dataset_drivers
