module netcdf_drivers

! This module creates, writes to, and closes netcdf files of the model outputs

! J. Melton
! Nov 2016

implicit none

!public :: netcdf_create
!public :: netcdf_output
!public :: netcdf_close
public :: check_nc
public :: parsecoords

contains

!--------------------------------------------------------------------------
! subroutine netcdf_create()
!
! use netcdf
! use io_driver, only : srtx,srty,cntx,cnty,endx,endy
! !use netcdf_error
! !use netcdf_drivers, only : handle_err
!
! implicit none
!
! integer :: lon,lat,layer,pft,month,time
! integer :: varid
! integer :: i
!
! character(8)  :: today
! character(10) :: now
!
! real, dimension(2) :: xrange
! real, dimension(2) :: yrange
!
! integer, dimension(1) :: pos
! integer, dimension(2) :: xpos,ypos
!
! integer, allocatable, dimension(:) :: pftnum
!
! !----------
!
! ! Find the range of this particular netcdf file.
!
! ! Set up boundaries
!
! srtx = max(1,srtx)
! srty = max(1,srty)
!
! cntx = min(cntx,inputlonlen)
! cnty = min(cnty,inputlatlen)
!
! endx = srtx + cntx - 1
! endy = srty + cnty - 1
!
! if (lonvect(srtx) < bounds(1) .and. bounds(2) /= bounds(1)) srtx = srtx + 1
! cntx = 1 + abs(maxval(xpos) - srtx)
!
! if (latvect(srty) < bounds(3) .and. bounds(4) /= bounds(3)) srty = srty + 1
! cnty = 1 + abs(maxval(ypos) - srty)
!
! !xrange(1) = minval(lonvect(srtx:endx))
! !xrange(2) = maxval(lonvect(srtx:endx)) + 0.5
!
! !yrange(1) = minval(latvect(srty:endy)) - 0.5
! !yrange(2) = maxval(latvect(srty:endy))
!
! !----1
!
! status = nf90_create(outputfile,cmode=nf90_clobber,ncid=ncid)
! if (status/=nf90_noerr) call handle_err(status)
!
! !write(*,*)'created output file ncid:',ncid
!
! status = nf90_put_att(ncid,nf90_global,'title','CLASS-CTEM netCDF output file')
! if (status/=nf90_noerr) call handle_err(status)
!
! call date_and_time(today,now)
!
! status = nf90_put_att(ncid,nf90_global,'timestamp',today//' '//now(1:4))
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_put_att(ncid,nf90_global,'Conventions','COARDS')
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_put_att(ncid,nf90_global,'node_offset',1)
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_put_att(ncid,nf90_global,'_Format',"netCDF-4")
! if (status/=nf90_noerr) call handle_err(status)
!
! !----2
!
! status = nf90_def_dim(ncid,'lon',cntx,lon)
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_def_var(ncid,'lon',nf90_float,lon,varid)
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_put_att(ncid,varid,'long_name','longitude')
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_put_att(ncid,varid,'units','degrees_east')
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_put_att(ncid,varid,'actual_range',xrange)
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_put_att(ncid,varid,'_Storage',"contiguous")
! if (status/=nf90_noerr) call handle_err(status)
!
!
! !----3
!
! status = nf90_def_dim(ncid,'lat',cnty,lat)
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_def_var(ncid,'lat',nf90_float,lat,varid)
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_put_att(ncid,varid,'long_name','latitude')
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_put_att(ncid,varid,'units','degrees_north')
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_put_att(ncid,varid,'actual_range',yrange)
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_put_att(ncid,varid,'_Storage',"contiguous")
! if (status/=nf90_noerr) call handle_err(status)
!
! !----4
!
! status = nf90_def_dim(ncid,'layer',nl,layer)
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_def_var(ncid,'layer',nf90_short,layer,varid)
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_put_att(ncid,varid,'long_name','soil layer')
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_put_att(ncid,varid,'units','layer')
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_put_att(ncid,varid,'_Storage',"contiguous")
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_put_att(ncid,varid,'_Endianness',"little")
! if (status/=nf90_noerr) call handle_err(status)
!
! !----5
!
! status = nf90_def_dim(ncid,'pft',npft,pft)
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_def_var(ncid,'pft',nf90_short,pft,varid)
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_put_att(ncid,varid,'long_name','Plant Functional Type')
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_put_att(ncid,varid,'units','PFT')
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_put_att(ncid,varid,'_Storage',"contiguous")
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_put_att(ncid,varid,'_Endianness',"little")
! if (status/=nf90_noerr) call handle_err(status)
!
! !----6
!
! status = nf90_def_dim(ncid,'month',12,month)
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_def_var(ncid,'month',nf90_short,month,varid)
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_put_att(ncid,varid,'long_name','month')
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_put_att(ncid,varid,'units','month')
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_put_att(ncid,varid,'_Storage',"contiguous")
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_put_att(ncid,varid,'_Endianness',"little")
! if (status/=nf90_noerr) call handle_err(status)
!
! !----7
!
! status = nf90_def_dim(ncid,'time',nf90_unlimited,time)
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_def_var(ncid,'time',nf90_int,time,varid)
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_put_att(ncid,varid,'long_name','time')
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_put_att(ncid,varid,'units','years since 1900')
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_put_att(ncid,varid,'_Storage',"chunked")
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_put_att(ncid,varid,'_Chunksizes',1)
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_put_att(ncid,varid,'_Endianness',"little")
! if (status/=nf90_noerr) call handle_err(status)
!
! !----------------------------------------------------------------------
! !7
! status = nf90_def_var(ncid,'mw1',nf90_float,[lon,lat,month,time],varid)
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_put_att(ncid,varid,'long_name','mean total soil column moisture (liq+ice) fraction of total porosity')
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_put_att(ncid,varid,'units','fraction')
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_put_att(ncid,varid,'_FillValue',fill_value)
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_put_att(ncid,varid,'missing_value',fill_value)
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_put_att(ncid,varid,'_Storage',"chunked")
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_put_att(ncid,varid,'_Chunksizes',[300,720,12,1])
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_put_att(ncid,varid,'_DeflateLevel',1)
! if (status/=nf90_noerr) call handle_err(status)
!
!
!
! !---
!
! status = nf90_enddef(ncid)
! if (status/=nf90_noerr) call handle_err(status)
!
! !----
!
! !write the dimension variables (except for time)
!
! !lonvect = lonvect + 0.25
! !latvect = latvect - 0.25
!
! allocate(pftnum(npft))
!
! forall (i=1:npft)
!   pftnum(i) = i
! end forall
!
! status = nf90_put_var(ncid,1,lonvect(srtx:endx))
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_put_var(ncid,2,latvect(srty:endy))
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_put_var(ncid,3,[1,2])
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_put_var(ncid,4,pftnum)
! if (status/=nf90_noerr) call handle_err(status)
!
! status = nf90_put_var(ncid,5,[1,2,3,4,5,6,7,8,9,10,11,12])
! if (status/=nf90_noerr) call handle_err(status)
!
! !----------------
!
! deallocate(pftnum)
!
! end subroutine netcdf_create

! !=======================================================================
! subroutine netcdf_output(yr)
!
! use netcdf
! use netcdf_error
!
! use iovariables
! use statevars,                  only : sv,veg
! use arveparams
!
! implicit none
!
! integer :: yr
! integer, dimension(1) :: lyear
! integer :: i,x,y
! real, allocatable, dimension(:,:,:) :: ovar
!
!
! !----
! lyear = yr
!
! !write the time variable
!
! status = nf90_put_var(ncid,6,lyear,start=[lyear],count=[1])
! if (status/=nf90_noerr) call handle_err(status)
!
! !----
! !write other variables
!
! !----
!
! allocate(ovar(cntx,cnty,12))
!
! ovar = fill_value
!
! i = 1
! do y = 1,cnty
!   do x = 1,cntx
!     if (sv(i)%valid_cell) ovar(x,y,:) = ov(i)%mw1
!     i = i + 1
!   end do
! end do
!
! status = nf90_put_var(ncid,7,ovar,start=[1,1,1,lyear],count=[cntx,cnty,12,1])
! if (status/=nf90_noerr) call handle_err(status)
!
! deallocate(ovar)
!
!
! end subroutine netcdf_output
!
!  !=============================================================================
!
! subroutine netcdf_close()
!
! use netcdf
! use netcdf_error
!
! implicit none
!
! status = nf90_close(ncid)
! if (status/=nf90_noerr) call handle_err(status)
!
! end subroutine netcdf_close
!
!=======================================================================

subroutine check_nc(status)

  use netcdf
  use typesizes

  implicit none

  !Internal subroutine - checks error status after each netcdf call,
  !prints out text message each time an error code is returned.

  integer, intent(in) :: status

  if(status /= nf90_noerr) then
    write(0,*)'netCDF error: ',trim(nf90_strerror(status))
    stop
  end if

end subroutine check_nc

!=======================================================================

subroutine parsecoords(coordstring,val)

!subroutine to parse a coordinate string

implicit none

character(45),      intent(in)  :: coordstring
real, dimension(4), intent(out) :: val

character(10), dimension(4) :: cval = '0'

integer :: i
integer :: lasti = 1
integer :: part  = 1

do i=1,len_trim(coordstring)
  if (coordstring(i:i) == '/') then
    cval(part) = coordstring(lasti:i-1)
    lasti=i+1
    part = part + 1
  end if
end do

cval(part) = coordstring(lasti:i-1)

read(cval,*)val

if (part < 4) then
  val(3)=val(2)
  val(4)=val(3)
  val(2)=val(1)
end if

end subroutine parsecoords

end module netcdf_drivers
