program netcdf_create_bf

! This creates a NetCDF files that CTEMCLASS output will be written
! into (using netcdf_writer). The creator_module.f90 contains a module 
! where any parameters must be defined. A simple joboptions file is also required.

! If MOSAIC=.true., you will write both the composite and mosaic variables to the 
! netcdf, in this case MAKEMONTHLY=.false. is likely a good idea to avoid a massive
! file size.

! If the run has many years, also MAKEMONTHLY=.false. could be a good idea

! DOFIRE=.true. turns on the disturbance vars, while COMPETE_LNDUSE=.true. will
! write the competition and luc PFT fractions

! This netcdf generations scheme is designed to be used with master_parallel.sh
! or quick_remake_netcdf.sh

! J. Melton and S. Griffith - May 14, 2012

! Joe Melton - June 19, 2013
!       Massively overhauled to allow easy addition of vars and the mosaic/composite,
!       disturbance and competition/land use vars. Also to make netcdf4 files, but
!       presently cview can not handle them.

use netcdf
use creator_module_bf

implicit none

integer :: totyrs
integer :: yrst
integer :: realyrst
character(120) :: file_to_write
character(120) :: long_path
character(180) :: file_to_write_extended
character(120) :: jobfile
logical :: CTEM
logical :: MAKEMONTHLY
logical :: DOFIRE
logical :: MOSAIC
logical :: COMPETE_LNDUSE
logical :: PARALLELRUN

!name list set up
namelist /joboptions/ &
  PARALLELRUN,        &
  CTEM, 	      &
  MAKEMONTHLY,        &
  DOFIRE,             &
  MOSAIC,             &
  COMPETE_LNDUSE,     &
  totyrs,             &
  yrst,               &
  realyrst,           &
  long_path,          &
  file_to_write                

!----------
! set up the number of longitudes and latitudes,
! as well as the ranges of each.

 call getarg(1,jobfile)

 open(10,file=jobfile,status='old')

 read(10,nml = joboptions)

 close(10)

!fill lonvect and latvect
allocate(lonvect(cntx))
allocate(latvect(cnty))

!pass the lon/lat values to the vector
lonvect=valslons
latvect=valslats

! these are not set accurately!
xrange(1) = minval(lonvect)
xrange(2) = maxval(lonvect)

yrange(1) = minval(latvect)
yrange(2) = maxval(latvect)

!----
! Create the NetCDF file

file_to_write_extended = trim(file_to_write)//'_CLASSCTEM.nc'
call create_netcdf(totyrs,yrst,realyrst,file_to_write_extended,CTEM,MOSAIC,MAKEMONTHLY,DOFIRE,COMPETE_LNDUSE)

deallocate(lonvect)
deallocate(latvect)

end program netcdf_create_bf

!=======================================================================

subroutine create_netcdf(totyrs,yrst,realyrst,file_to_write,CTEM,MOSAIC,MAKEMONTHLY,DOFIRE,COMPETE_LNDUSE)

use creator_module_bf
use netcdf

implicit none

integer :: i
integer,intent(in) :: totyrs
integer,intent(in) :: yrst
integer,intent(in) :: realyrst
character(120),intent(in) :: file_to_write
logical,intent(in) :: CTEM
logical,intent(in) :: MOSAIC
logical,intent(in) :: MAKEMONTHLY
logical,intent(in) :: DOFIRE
logical,intent(in) :: COMPETE_LNDUSE

integer :: grpid_ann_ctem
integer :: grpid_ann_class
integer :: grpid_mon_ctem
integer :: grpid_mon_class
integer :: grpid_mon_dist
integer :: grpid_ann_dist

integer :: grpid_ann_ctem_t
integer :: grpid_mon_ctem_t
integer :: grpid_mon_dist_t
integer :: grpid_ann_dist_t


real, allocatable, dimension(:,:,:) :: threevar
real, allocatable, dimension(:,:,:,:) :: fourvar
real, allocatable, dimension(:,:,:,:,:) :: fivevar
integer, allocatable, dimension(:) :: layersnum
integer, allocatable, dimension(:) :: pftnum
integer, allocatable, dimension(:) :: tilesnum
integer, allocatable, dimension(:) :: timenum

character(8)  :: today
character(10) :: now
character(20) :: yrssince
integer :: z
character(4) :: zchar

if (net4) then
  status = nf90_create(file_to_write,cmode=nf90_netcdf4,ncid=ncid)
  if (status/=nf90_noerr) call handle_err(status)
else
  status = nf90_create(file_to_write,cmode=nf90_64bit_offset,ncid=ncid)
  if (status/=nf90_noerr) call handle_err(status)
end if

status = nf90_put_att(ncid,nf90_global,'title','CLASS-CTEM netCDF output file')
if (status/=nf90_noerr) call handle_err(status)

call date_and_time(today,now)

status = nf90_put_att(ncid,nf90_global,'timestamp',today//' '//now(1:4))
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,nf90_global,'Conventions','COARDS')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,nf90_global,'node_offset',1)
if (status/=nf90_noerr) call handle_err(status)

if (net4) then
  status = nf90_put_att(ncid,nf90_global,'_Format',"netCDF-4")
  if (status/=nf90_noerr) call handle_err(status)
else
  status = nf90_put_att(ncid,nf90_global,'_Format',"netCDF-3")
  if (status/=nf90_noerr) call handle_err(status)
end if

!----1
status = nf90_def_dim(ncid,'lon',cntx,lon)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_def_var(ncid,'lon',nf90_float,lon,varid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','longitude')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units','degrees_east')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'actual_range',xrange)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Storage',"contiguous")
if (status/=nf90_noerr) call handle_err(status)


!----2
status = nf90_def_dim(ncid,'lat',cnty,lat)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_def_var(ncid,'lat',nf90_float,lat,varid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','latitude')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units','degrees_north')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'actual_range',yrange)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Storage',"contiguous")
if (status/=nf90_noerr) call handle_err(status)

!----3
status = nf90_def_dim(ncid,'tile',ntile,tile)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_def_var(ncid,'tile',nf90_short,tile,varid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','mosaic tile')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units','tile number')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Storage',"contiguous")
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Endianness',"little")
if (status/=nf90_noerr) call handle_err(status)

!----4
status = nf90_def_dim(ncid,'pft',ctemnpft,pft)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_def_var(ncid,'pft',nf90_short,pft,varid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','Plant Functional Type')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units','PFT')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Storage',"contiguous")
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Endianness',"little")
if (status/=nf90_noerr) call handle_err(status)

!----5
status = nf90_def_dim(ncid,'layer',nl,layer)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_def_var(ncid,'layer',nf90_short,layer,varid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','soil layer')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units','layer')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Storage',"contiguous")
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Endianness',"little")
if (status/=nf90_noerr) call handle_err(status)

!----6
if (MAKEMONTHLY) then
 status = nf90_def_dim(ncid,'month',12,month)
 if (status/=nf90_noerr) call handle_err(status)

 status = nf90_def_var(ncid,'month',nf90_short,month,varid)
 if (status/=nf90_noerr) call handle_err(status)

 status = nf90_put_att(ncid,varid,'long_name','month')
 if (status/=nf90_noerr) call handle_err(status)

 status = nf90_put_att(ncid,varid,'units','month')
 if (status/=nf90_noerr) call handle_err(status)

 status = nf90_put_att(ncid,varid,'_Storage',"contiguous")
 if (status/=nf90_noerr) call handle_err(status)

 status = nf90_put_att(ncid,varid,'_Endianness',"little")
 if (status/=nf90_noerr) call handle_err(status)
end if

!----7
z=realyrst-1
write (zchar, '(I4)') z
yrssince='years since '//zchar//' AD'

status = nf90_def_dim(ncid,'time',totyrs,time)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_def_var(ncid,'time',nf90_int,time,varid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','time')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units',yrssince)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Storage',"chunked")
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Chunksizes',1)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Endianness',"little")
if (status/=nf90_noerr) call handle_err(status)

!==============Create the groups====================

if (net4) then

  status = nf90_def_grp(ncid,'CLASS-Annual', grpid_ann_class)
  if (status /= nf90_noerr) call handle_err(status)

  if (MAKEMONTHLY) then
   status = nf90_def_grp(ncid,'CLASS-Monthly', grpid_mon_class)
   if (status /= nf90_noerr) call handle_err(status)
  end if

  if (CTEM) then

    status = nf90_def_grp(ncid,'CTEM-Annual GridAvg', grpid_ann_ctem)
    if (status /= nf90_noerr) call handle_err(status)

    if (DOFIRE) then
        status = nf90_def_grp(grpid_ann_ctem,'Annual-Disturbance GridAvg', grpid_ann_dist)
        if (status /= nf90_noerr) call handle_err(status)
    end if

    if (MAKEMONTHLY) then
     status = nf90_def_grp(ncid,'CTEM-Monthly GridAvg', grpid_mon_ctem)
     if (status /= nf90_noerr) call handle_err(status)

     if (DOFIRE) then
        status = nf90_def_grp(grpid_mon_ctem,'Monthly-Disturbance GridAvg', grpid_mon_dist)
        if (status /= nf90_noerr) call handle_err(status)
     end if
    end if !makemonthly

   if (MOSAIC) then

    status = nf90_def_grp(ncid,'CTEM-Annual Tiled', grpid_ann_ctem_t)
    if (status /= nf90_noerr) call handle_err(status)

    if (DOFIRE) then
        status = nf90_def_grp(grpid_ann_ctem_t,'Annual-Disturbance Tiled', grpid_ann_dist_t)
        if (status /= nf90_noerr) call handle_err(status)
    end if

    if (MAKEMONTHLY) then
     status = nf90_def_grp(ncid,'CTEM-Monthly Tiled', grpid_mon_ctem_t)
     if (status /= nf90_noerr) call handle_err(status)

     if (DOFIRE) then
        status = nf90_def_grp(grpid_mon_ctem_t,'Monthly-Disturbance Tiled', grpid_mon_dist_t)
        if (status /= nf90_noerr) call handle_err(status)
     end if
    end if !makemonthly

   end if ! mosaic
  end if  !ctem

else ! netcf3

! No groups allowed in netcdf3 so set all to ncid value
        grpid_ann_ctem=ncid
        grpid_ann_class=ncid
        grpid_mon_ctem=ncid
        grpid_mon_class=ncid
        grpid_mon_dist=ncid
        grpid_ann_dist=ncid
    if (MOSAIC) then
        grpid_ann_ctem_t=ncid
        grpid_mon_ctem_t=ncid
        grpid_mon_dist_t=ncid
        grpid_ann_dist_t=ncid
    end if
end if

status = nf90_enddef(ncid) 
if (status/=nf90_noerr) call handle_err(status)

!========= Fill in the parameters =======================

allocate(tilesnum(ntile))
allocate(pftnum(ctemnpft))
allocate(layersnum(nl))
allocate(timenum(totyrs))

forall (i=1:ntile)
  tilesnum(i) = i
end forall

forall (i=1:ctemnpft)
  pftnum(i) = i
end forall

forall (i=1:nl)
  layersnum(i) = i
end forall

forall (i=yrst:totyrs)
   timenum(i) = i
end forall

status = nf90_put_var(ncid,lon,lonvect)  !longitudes
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_var(ncid,lat,latvect)  !latitudes
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_var(ncid,tile,tilesnum)  !tiles
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_var(ncid,pft,pftnum)    !number of pfts
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_var(ncid,layer,layersnum)    
if (status/=nf90_noerr) call handle_err(status)

if (MAKEMONTHLY) then
status = nf90_put_var(ncid,month,[1,2,3,4,5,6,7,8,9,10,11,12]) !months
if (status/=nf90_noerr) call handle_err(status)
end if

status = nf90_put_var(ncid,time,timenum)    
if (status/=nf90_noerr) call handle_err(status)

deallocate(tilesnum)
deallocate(pftnum)
deallocate(layersnum)
deallocate(timenum)

!================ Now define the variables===========
! CTEM first
if (CTEM) then
 
 !Monthly CTEM MOSAIC==============================
 if (MOSAIC) then 

   allocate(fourvar(cntx,cnty,ntile,totyrs))
   fourvar=fill_value

  if (MAKEMONTHLY) then

   allocate(fivevar(cntx,cnty,ntile,12,totyrs))
   fivevar=fill_value

  do i=lbound(CTEM_M_VAR,1), ubound(CTEM_M_VAR,1)
  
   status = nf90_redef(grpid_mon_ctem_t) 
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_def_var(grpid_mon_ctem_t,trim(CTEM_M_VAR(i)),nf90_float,[lon,lat,tile,month,time],varid)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_att(grpid_mon_ctem_t,varid,'long_name',trim(CTEM_M_NAME(i)))
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_att(grpid_mon_ctem_t,varid,'units',trim(CTEM_M_UNIT(i)))
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_att(grpid_mon_ctem_t,varid,'_FillValue',fill_value)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_att(grpid_mon_ctem_t,varid,'missing_value',fill_value)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_att(grpid_mon_ctem_t,varid,'_Storage',"chunked")
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_att(grpid_mon_ctem_t,varid,'_Chunksizes',[cnty,cntx,ntile,12,1])
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_att(grpid_mon_ctem_t,varid,'_DeflateLevel',1)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_enddef(grpid_mon_ctem_t)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_var(grpid_mon_ctem_t,varid,fivevar,start=[1,1,1,1,1],count=[cntx,cnty,ntile,12,totyrs])
   if (status/=nf90_noerr) call handle_err(status)

  end do

 ! MONTHLY MOSAIC DISTURBANCE VARIABLES

  if (DOFIRE) then

    do i=lbound(CTEM_M_D_VAR,1), ubound(CTEM_M_D_VAR,1)

     status = nf90_redef(grpid_mon_dist_t) 
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_def_var(grpid_mon_dist_t,trim(CTEM_M_D_VAR(i)),nf90_float,[lon,lat,tile,month,time],varid)
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_mon_dist_t,varid,'long_name',trim(CTEM_M_D_NAME(i)))
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_mon_dist_t,varid,'units',trim(CTEM_M_D_UNIT(i)))
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_mon_dist_t,varid,'_FillValue',fill_value)
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_mon_dist_t,varid,'missing_value',fill_value)
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_mon_dist_t,varid,'_Storage',"chunked")
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_mon_dist_t,varid,'_Chunksizes',[cnty,cntx,ntile,12,1])
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_mon_dist_t,varid,'_DeflateLevel',1)
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_enddef(grpid_mon_dist_t)
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_var(grpid_mon_dist_t,varid,fivevar,start=[1,1,1,1,1],count=[cntx,cnty,ntile,12,totyrs])
     if (status/=nf90_noerr) call handle_err(status)

    end do

  end if !dofire
 
  end if !makemonthly
  
  !Annual CTEM MOSAIC====================================

  do i=lbound(CTEM_Y_VAR,1), ubound(CTEM_Y_VAR,1)

   status = nf90_redef(grpid_ann_ctem_t) 
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_def_var(grpid_ann_ctem_t,trim(CTEM_Y_VAR(i)),nf90_float,[lon,lat,tile,time],varid)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_att(grpid_ann_ctem_t,varid,'long_name',trim(CTEM_Y_NAME(i)))
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_att(grpid_ann_ctem_t,varid,'units',trim(CTEM_Y_UNIT(i)))
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_att(grpid_ann_ctem_t,varid,'_FillValue',fill_value)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_att(grpid_ann_ctem_t,varid,'missing_value',fill_value)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_att(grpid_ann_ctem_t,varid,'_Storage',"chunked")
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_att(grpid_ann_ctem_t,varid,'_Chunksizes',[cnty,cntx,ntile,1])
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_att(grpid_ann_ctem_t,varid,'_DeflateLevel',1)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_enddef(grpid_ann_ctem_t)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_var(grpid_ann_ctem_t,varid,fourvar,start=[1,1,1,1],count=[cntx,cnty,ntile,totyrs])
   if (status/=nf90_noerr) call handle_err(status)

  enddo

 ! Annual DISTURBANCE MOSAIC VARIABLES

   if (DOFIRE) then

    do i=lbound(CTEM_Y_D_VAR,1), ubound(CTEM_Y_D_VAR,1)

     status = nf90_redef(grpid_ann_dist_t) 
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_def_var(grpid_ann_dist_t,trim(CTEM_Y_D_VAR(i)),nf90_float,[lon,lat,tile,time],varid)
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_ann_dist_t,varid,'long_name',trim(CTEM_Y_D_NAME(i)))
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_ann_dist_t,varid,'units',trim(CTEM_Y_D_UNIT(i)))
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_ann_dist_t,varid,'_FillValue',fill_value)
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_ann_dist_t,varid,'missing_value',fill_value)
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_ann_dist_t,varid,'_Storage',"chunked")
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_ann_dist_t,varid,'_Chunksizes',[cnty,cntx,ntile,1])
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_ann_dist_t,varid,'_DeflateLevel',1)
     if (status/=nf90_noerr) call handle_err(status)

    status = nf90_enddef(grpid_ann_dist_t)
    if (status/=nf90_noerr) call handle_err(status)
 
    status = nf90_put_var(grpid_ann_dist_t,varid,fourvar,start=[1,1,1,1],count=[cntx,cnty,ntile,totyrs])
    if (status/=nf90_noerr) call handle_err(status)

    end do

   end if !dofire

   deallocate(fourvar)
   if (MAKEMONTHLY) then
    deallocate(fivevar)
   end if

!else if (.NOT. MOSAIC) then
 end if !mosaic
 
   allocate(threevar(cntx,cnty,totyrs))
   threevar=fill_value

 if (MAKEMONTHLY) then

   allocate(fourvar(cntx,cnty,12,totyrs))
   fourvar=fill_value

 !Monthly CTEM COMPOSITE===================================================================================

  do i=lbound(CTEM_M_VAR,1), ubound(CTEM_M_VAR,1)

   status = nf90_redef(grpid_mon_ctem) 
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_def_var(grpid_mon_ctem,trim(CTEM_M_VAR(i)),nf90_float,[lon,lat,month,time],varid)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_att(grpid_mon_ctem,varid,'long_name',trim(CTEM_M_NAME(i)))
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_att(grpid_mon_ctem,varid,'units',trim(CTEM_M_UNIT(i)))
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_att(grpid_mon_ctem,varid,'_FillValue',fill_value)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_att(grpid_mon_ctem,varid,'missing_value',fill_value)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_att(grpid_mon_ctem,varid,'_Storage',"chunked")
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_att(grpid_mon_ctem,varid,'_Chunksizes',[cnty,cntx,12,1])
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_att(grpid_mon_ctem,varid,'_DeflateLevel',1)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_enddef(grpid_mon_ctem)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_var(grpid_mon_ctem,varid,fourvar,start=[1,1,1,1],count=[cntx,cnty,12,totyrs])
   if (status/=nf90_noerr) call handle_err(status)
 
  end do

 ! DISTURBANCE VARIABLES

   if (DOFIRE) then

    do i=lbound(CTEM_M_D_VAR,1), ubound(CTEM_M_D_VAR,1)

     status = nf90_redef(grpid_mon_dist) 
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_def_var(grpid_mon_dist,trim(CTEM_M_D_VAR(i)),nf90_float,[lon,lat,month,time],varid)
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_mon_dist,varid,'long_name',trim(CTEM_M_D_NAME(i)))
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_mon_dist,varid,'units',trim(CTEM_M_D_UNIT(i)))
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_mon_dist,varid,'_FillValue',fill_value)
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_mon_dist,varid,'missing_value',fill_value)
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_mon_dist,varid,'_Storage',"chunked")
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_mon_dist,varid,'_Chunksizes',[cnty,cntx,12,1])
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_mon_dist,varid,'_DeflateLevel',1)
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_enddef(grpid_mon_dist)
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_var(grpid_mon_dist,varid,fourvar,start=[1,1,1,1],count=[cntx,cnty,12,totyrs])
     if (status/=nf90_noerr) call handle_err(status)

    end do
   end if !dofire
  end if ! makemonthly

 !=============Annual CTEM COMPOSITE=================================================

  do i=lbound(CTEM_Y_VAR,1), ubound(CTEM_Y_VAR,1)

   status = nf90_redef(grpid_ann_ctem) 
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_def_var(grpid_ann_ctem,trim(CTEM_Y_VAR(i)),nf90_float,[lon,lat,time],varid)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_att(grpid_ann_ctem,varid,'long_name',trim(CTEM_Y_NAME(i)))
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_att(grpid_ann_ctem,varid,'units',trim(CTEM_Y_UNIT(i)))
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_att(grpid_ann_ctem,varid,'_FillValue',fill_value)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_att(grpid_ann_ctem,varid,'missing_value',fill_value)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_att(grpid_ann_ctem,varid,'_Storage',"chunked")
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_att(grpid_ann_ctem,varid,'_Chunksizes',[cnty,cntx,1])
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_att(grpid_ann_ctem,varid,'_DeflateLevel',1)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_enddef(grpid_ann_ctem)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_var(grpid_ann_ctem,varid,threevar,start=[1,1,1],count=[cntx,cnty,totyrs])
   if (status/=nf90_noerr) call handle_err(status)

  enddo

 ! DISTURBANCE VARIABLES

   if (DOFIRE) then

    do i=lbound(CTEM_Y_D_VAR,1), ubound(CTEM_Y_D_VAR,1)

     status = nf90_redef(grpid_ann_dist) 
     if (status/=nf90_noerr) call handle_err(status)
 
     status = nf90_def_var(grpid_ann_dist,trim(CTEM_Y_D_VAR(i)),nf90_float,[lon,lat,time],varid)
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_ann_dist,varid,'long_name',trim(CTEM_Y_D_NAME(i)))
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_ann_dist,varid,'units',trim(CTEM_Y_D_UNIT(i)))
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_ann_dist,varid,'_FillValue',fill_value)
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_ann_dist,varid,'missing_value',fill_value)
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_ann_dist,varid,'_Storage',"chunked")
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_ann_dist,varid,'_Chunksizes',[cnty,cntx,1])
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_ann_dist,varid,'_DeflateLevel',1)
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_enddef(grpid_ann_dist)
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_var(grpid_ann_dist,varid,threevar,start=[1,1,1],count=[cntx,cnty,totyrs])
     if (status/=nf90_noerr) call handle_err(status)

    end do

   end if !dofire

   deallocate(threevar)
   if (MAKEMONTHLY) then
    deallocate(fourvar)
   end if

!end if !mosaic/composite

if (COMPETE_LNDUSE) then

 if (MAKEMONTHLY) then

   allocate(fourvar(cntx,cnty,12,totyrs))
   fourvar=fill_value

   allocate(fivevar(cntx,cnty,ctemnpft,12,totyrs))
   fivevar=fill_value

!============Monthly Competition/Land Use CTEM====================

! Here the first one is different than the rest

    do i=lbound(CTEM_M_C_VAR,1)+1, ubound(CTEM_M_C_VAR,1)

     status = nf90_redef(grpid_mon_ctem) 
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_def_var(grpid_mon_ctem,trim(CTEM_M_C_VAR(i)),nf90_float,[lon,lat,pft,month,time],varid)
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_mon_ctem,varid,'long_name',trim(CTEM_M_C_NAME(1)))
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_mon_ctem,varid,'units',trim(CTEM_M_C_UNIT(1)))
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_mon_ctem,varid,'_FillValue',fill_value)
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_mon_ctem,varid,'missing_value',fill_value)
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_mon_ctem,varid,'_Storage',"chunked")
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_mon_ctem,varid,'_Chunksizes',[cnty,cntx,ctemnpft,12,1])
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_mon_ctem,varid,'_DeflateLevel',1)
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_enddef(grpid_mon_ctem)
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_var(grpid_mon_ctem,varid,fivevar,start=[1,1,1,1,1],count=[cntx,cnty,ctemnpft,12,totyrs])
     if (status/=nf90_noerr) call handle_err(status)

    end do

     ! Total plant cover

     status = nf90_redef(grpid_mon_ctem) 
     if (status/=nf90_noerr) call handle_err(status)
 
     status = nf90_def_var(grpid_mon_ctem,trim(CTEM_M_C_VAR(1)),nf90_float,[lon,lat,month,time],varid)
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_mon_ctem,varid,'long_name',trim(CTEM_M_C_NAME(2)))
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_mon_ctem,varid,'units',trim(CTEM_M_C_UNIT(2)))
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_mon_ctem,varid,'_FillValue',fill_value)
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_mon_ctem,varid,'missing_value',fill_value)
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_mon_ctem,varid,'_Storage',"chunked")
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_mon_ctem,varid,'_Chunksizes',[cnty,cntx,12,1])
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_mon_ctem,varid,'_DeflateLevel',1)
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_enddef(grpid_mon_ctem)
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_var(grpid_mon_ctem,varid,fourvar,start=[1,1,1,1],count=[cntx,cnty,12,totyrs])
     if (status/=nf90_noerr) call handle_err(status)

     deallocate(fivevar)
     deallocate(fourvar)

   end if ! makemonthly

!============Annual Competition/ Land Use CTEM====================

     allocate(fourvar(cntx,cnty,ctemnpft,totyrs))
     fourvar=fill_value

     allocate(threevar(cntx,cnty,totyrs))
     threevar=fill_value

    do i=lbound(CTEM_Y_C_VAR,1)+1, ubound(CTEM_Y_C_VAR,1)

     status = nf90_redef(grpid_ann_ctem) 
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_def_var(grpid_ann_ctem,trim(CTEM_Y_C_VAR(i)),nf90_float,[lon,lat,pft,time],varid)
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_ann_ctem,varid,'long_name',trim(CTEM_Y_C_NAME(1)))
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_ann_ctem,varid,'units',trim(CTEM_Y_C_UNIT(1)))
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_ann_ctem,varid,'_FillValue',fill_value)
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_ann_ctem,varid,'missing_value',fill_value)
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_ann_ctem,varid,'_Storage',"chunked")
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_ann_ctem,varid,'_Chunksizes',[cnty,cntx,ctemnpft,1])
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_ann_ctem,varid,'_DeflateLevel',1)
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_enddef(grpid_ann_ctem)
     if (status/=nf90_noerr) call handle_err(status)
  
     status = nf90_put_var(grpid_ann_ctem,varid,fourvar,start=[1,1,1,1],count=[cntx,cnty,ctemnpft,totyrs])
     if (status/=nf90_noerr) call handle_err(status)

    end do

     ! Total Plant cover
     status = nf90_redef(grpid_ann_ctem) 
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_def_var(grpid_ann_ctem,trim(CTEM_Y_C_VAR(1)),nf90_float,[lon,lat,time],varid)
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_ann_ctem,varid,'long_name',trim(CTEM_Y_C_NAME(2)))
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_ann_ctem,varid,'units',trim(CTEM_Y_C_UNIT(2)))
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_ann_ctem,varid,'_FillValue',fill_value)
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_ann_ctem,varid,'missing_value',fill_value)
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_ann_ctem,varid,'_Storage',"chunked")
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_ann_ctem,varid,'_Chunksizes',[cnty,cntx,1])
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_att(grpid_ann_ctem,varid,'_DeflateLevel',1)
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_enddef(grpid_ann_ctem)
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_var(grpid_ann_ctem,varid,threevar,start=[1,1,1],count=[cntx,cnty,totyrs])
     if (status/=nf90_noerr) call handle_err(status)

     deallocate(fourvar)
     deallocate(threevar)

end if

end if !CTEMboolean

!============Monthly CLASS MOSAIC====================

! In the future we might want to output the mosaic output of CLASS, right
! now we do not so this can be commented out.

!if (MOSAIC) then

! do i=lbound(CLASS_M_VAR,1), ubound(CLASS_M_VAR,1)
!  status = nf90_def_var(grpid_mon_class,trim(CLASS_M_VAR(i)),nf90_float,[lon,lat,tile,month,time],varid)
!  if (status/=nf90_noerr) call handle_err(status)

!  status = nf90_put_att(grpid_mon_class,varid,'long_name',trim(CLASS_M_NAME(i)))
!  if (status/=nf90_noerr) call handle_err(status)

!  status = nf90_put_att(grpid_mon_class,varid,'units',trim(CLASS_M_UNIT(i)))
!  if (status/=nf90_noerr) call handle_err(status)

!  status = nf90_put_att(grpid_mon_class,varid,'_FillValue',fill_value)
!  if (status/=nf90_noerr) call handle_err(status)

!  status = nf90_put_att(grpid_mon_class,varid,'missing_value',fill_value)
!  if (status/=nf90_noerr) call handle_err(status)

!  status = nf90_put_att(grpid_mon_class,varid,'_Storage',"chunked")
!  if (status/=nf90_noerr) call handle_err(status)

!  status = nf90_put_att(grpid_mon_class,varid,'_Chunksizes',[cnty,cntx,ntile,12,1])
!  if (status/=nf90_noerr) call handle_err(status)

!  status = nf90_put_att(grpid_mon_class,varid,'_DeflateLevel',1)
!  if (status/=nf90_noerr) call handle_err(status)

! enddo

!!============Monthly CLASS COMPOSITE====================

!else if (.NOT. MOSAIC) then

 if (MAKEMONTHLY) then

   allocate(fourvar(cntx,cnty,12,totyrs))
   fourvar=fill_value

   allocate(fivevar(cntx,cnty,nl,12,totyrs))
   fivevar=fill_value

 do i=lbound(CLASS_M_VAR,1), ubound(CLASS_M_VAR,1)

   status = nf90_redef(grpid_mon_class) 
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_def_var(grpid_mon_class,trim(CLASS_M_VAR(i)),nf90_float,[lon,lat,month,time],varid)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_att(grpid_mon_class,varid,'_Chunksizes',[cnty,cntx,12,1])
   if (status/=nf90_noerr) call handle_err(status)

  status = nf90_put_att(grpid_mon_class,varid,'long_name',trim(CLASS_M_NAME(i)))
  if (status/=nf90_noerr) call handle_err(status)

  status = nf90_put_att(grpid_mon_class,varid,'units',trim(CLASS_M_UNIT(i)))
  if (status/=nf90_noerr) call handle_err(status)

  status = nf90_put_att(grpid_mon_class,varid,'_FillValue',fill_value)
  if (status/=nf90_noerr) call handle_err(status)

  status = nf90_put_att(grpid_mon_class,varid,'missing_value',fill_value)
  if (status/=nf90_noerr) call handle_err(status)

  status = nf90_put_att(grpid_mon_class,varid,'_Storage',"chunked")
  if (status/=nf90_noerr) call handle_err(status)

  status = nf90_put_att(grpid_mon_class,varid,'_DeflateLevel',1)
  if (status/=nf90_noerr) call handle_err(status)

   status = nf90_enddef(grpid_mon_class)
   if (status/=nf90_noerr) call handle_err(status)

  status = nf90_put_var(grpid_mon_class,varid,fourvar,start=[1,1,1,1],count=[cntx,cnty,12,totyrs])
  if (status/=nf90_noerr) call handle_err(status)

 end do

  ! Now do the soil layer variables
 do i=lbound(CLASS_M_S_VAR,1), ubound(CLASS_M_S_VAR,1)

  status = nf90_redef(grpid_mon_class) 
  if (status/=nf90_noerr) call handle_err(status)

  status = nf90_def_var(grpid_mon_class,trim(CLASS_M_S_VAR(i)),nf90_float,[lon,lat,layer,month,time],varid)
  if (status/=nf90_noerr) call handle_err(status)

  status = nf90_put_att(grpid_mon_class,varid,'_Chunksizes',[cnty,cntx,nl,12,1])
  if (status/=nf90_noerr) call handle_err(status)

  status = nf90_put_att(grpid_mon_class,varid,'long_name',trim(CLASS_M_S_NAME(i)))
  if (status/=nf90_noerr) call handle_err(status)

  status = nf90_put_att(grpid_mon_class,varid,'units',trim(CLASS_M_S_UNIT(i)))
  if (status/=nf90_noerr) call handle_err(status)

  status = nf90_put_att(grpid_mon_class,varid,'_FillValue',fill_value)
  if (status/=nf90_noerr) call handle_err(status)

  status = nf90_put_att(grpid_mon_class,varid,'missing_value',fill_value)
  if (status/=nf90_noerr) call handle_err(status)

  status = nf90_put_att(grpid_mon_class,varid,'_Storage',"chunked")
  if (status/=nf90_noerr) call handle_err(status)

  status = nf90_put_att(grpid_mon_class,varid,'_DeflateLevel',1)
  if (status/=nf90_noerr) call handle_err(status)

  status = nf90_enddef(grpid_mon_class)
  if (status/=nf90_noerr) call handle_err(status)

  status = nf90_put_var(grpid_mon_class,varid,fivevar,start=[1,1,1,1,1],count=[cntx,cnty,nl,12,totyrs])
  if (status/=nf90_noerr) call handle_err(status)

 end do

  deallocate(fivevar)
  deallocate(fourvar)

  end if !makemonthly

!end if !mosaic/composite

!=============Annual CLASS MOSAIC==================

! In the future we might want to output the mosaic output of CLASS, right
! now we do not so this can be commented out.

!if (MOSAIC) then

! do i=lbound(CLASS_Y_VAR,1), ubound(CLASS_Y_VAR,1)

!  status = nf90_def_var(grpid_ann_class,trim(CLASS_Y_VAR(i)),nf90_float,[lon,lat,tile,time],varid)
!  if (status/=nf90_noerr) call handle_err(status)

!  status = nf90_put_att(grpid_ann_class,varid,'long_name',trim(CLASS_Y_NAME(i)))
!  if (status/=nf90_noerr) call handle_err(status)

!  status = nf90_put_att(grpid_ann_class,varid,'units',trim(CLASS_Y_UNIT(i)))
!  if (status/=nf90_noerr) call handle_err(status)

!  status = nf90_put_att(grpid_ann_class,varid,'_FillValue',fill_value)
!  if (status/=nf90_noerr) call handle_err(status)

!  status = nf90_put_att(grpid_ann_class,varid,'missing_value',fill_value)
!  if (status/=nf90_noerr) call handle_err(status)

!  status = nf90_put_att(grpid_ann_class,varid,'_Storage',"chunked")
!  if (status/=nf90_noerr) call handle_err(status)

!  status = nf90_put_att(grpid_ann_class,varid,'_Chunksizes',[cnty,cntx,ntile,1])
!  if (status/=nf90_noerr) call handle_err(status)

!  status = nf90_put_att(grpid_ann_class,varid,'_DeflateLevel',1)
!  if (status/=nf90_noerr) call handle_err(status)

! enddo

!!=============Annual CLASS COMPOSITE==================

!else if (.NOT. MOSAIC) then

   allocate(threevar(cntx,cnty,totyrs))
   threevar=fill_value

 do i=lbound(CLASS_Y_VAR,1), ubound(CLASS_Y_VAR,1)

  status = nf90_redef(grpid_ann_class) 
  if (status/=nf90_noerr) call handle_err(status)

  status = nf90_def_var(grpid_ann_class,trim(CLASS_Y_VAR(i)),nf90_float,[lon,lat,time],varid)
  if (status/=nf90_noerr) call handle_err(status)

  status = nf90_put_att(grpid_ann_class,varid,'long_name',trim(CLASS_Y_NAME(i)))
  if (status/=nf90_noerr) call handle_err(status)

  status = nf90_put_att(grpid_ann_class,varid,'units',trim(CLASS_Y_UNIT(i)))
  if (status/=nf90_noerr) call handle_err(status)

  status = nf90_put_att(grpid_ann_class,varid,'_FillValue',fill_value)
  if (status/=nf90_noerr) call handle_err(status)

  status = nf90_put_att(grpid_ann_class,varid,'missing_value',fill_value)
  if (status/=nf90_noerr) call handle_err(status)

  status = nf90_put_att(grpid_ann_class,varid,'_Storage',"chunked")
  if (status/=nf90_noerr) call handle_err(status)

  status = nf90_put_att(grpid_ann_class,varid,'_Chunksizes',[cnty,cntx,1])
  if (status/=nf90_noerr) call handle_err(status)

  status = nf90_put_att(grpid_ann_class,varid,'_DeflateLevel',1)
  if (status/=nf90_noerr) call handle_err(status)

  status = nf90_enddef(grpid_ann_class)
  if (status/=nf90_noerr) call handle_err(status)

  status = nf90_put_var(grpid_ann_class,varid,threevar,start=[1,1,1],count=[cntx,cnty,totyrs])
  if (status/=nf90_noerr) call handle_err(status)

 enddo

  deallocate(threevar)

!endif !mosaic/composite

!close the netcdf
status = nf90_close(ncid)
if (status/=nf90_noerr) call handle_err(status)

end subroutine create_netcdf
!------------------------------
subroutine handle_err(status)

  use netcdf

  implicit none

  !Internal subroutine - checks error status after each netcdf call,
  !prints out text message each time an error code is returned.

  integer, intent(in) :: status

  if(status /= nf90_noerr) then
    write(0,*)'netCDF error: ',trim(nf90_strerror(status))
    stop
  end if

end subroutine handle_err

!----------------


