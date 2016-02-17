program netcdf_create_bf

! This creates a NetCDF files that CTEMCLASS output will be written
! into (using netcdf_writer). The creator_module.f90 contains a module 
! where any parameters must be defined. A simple joboptions file is also required.

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

! Joe Melton - July 17, 2013
!       Changed how composite files are handled. Now I output them from CLASS-CTEM
!       with per PFT values as well as grid averaged. This code now can deal with that.
!
! Joe Melton - Feb 26 2014
!       Add in lat and lon bounds to the netcdf files
!
! Joe Melton - Jul 28 2014
!       Split creation into an annual and a monthly file
!
! Joe Melton - Feb 15 2016
!       Adapt so it can deal with tiling (different form than before) and push the
!       month into the years so I can reduce the dimension of the monthly files.
!       Also made it so the number of soil layers comes in from the job options.
!       Added the dothevars subroutine to try and streamline.

!----------------

use netcdf
use creator_module_bf

implicit none

integer :: totyrs
integer :: monyrs
integer :: yrst
integer :: realyrst
integer :: nl ! number of soil layers
integer :: i
integer :: adjustyr
character(120) :: file_to_write
character(120) :: long_path
character(180) :: file_to_write_extended
character(120) :: jobfile
logical :: CTEM
logical :: MAKEMONTHLY
logical :: DOFIRE
logical :: TILED
logical :: DOWETLANDS
logical :: COMPETE_LNDUSE
logical :: PARALLELRUN
logical :: DOPFTS

!name list set up
namelist /joboptions/ &
  PARALLELRUN,        &
  CTEM,               &
  MAKEMONTHLY,        &
  DOFIRE,             &
  TILED,              &
  DOPFTS,             &
  DOWETLANDS,         &
  COMPETE_LNDUSE,     &
  nl,                 &
  totyrs,             &
  monyrs,             &
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

allocate(latboundsvect(cnty+1))
allocate(lonboundsvect(cntx+1))

allocate(latbound(cnty+1,2))
allocate(lonbound(cntx+1,2))

!pass the lon/lat values to the vector
lonvect=valslons
latvect=valslats

!create the bounds vectors
latboundsvect(1)=-90.
do i = 1,cnty-1
  latboundsvect(i+1)=valslats(i)-(valslats(i)-valslats(i+1))*0.5   
end do
latboundsvect(cnty+1)=90.

do i = 2,cnty+1
   latbound(i-1,1)=latboundsvect(i-1)
   latbound(i-1,2)=latboundsvect(i)
end do


lonboundsvect(1)=0.
do i = 1,cntx-1
  lonboundsvect(i+1)=valslons(i)-(valslons(i)-valslons(i+1))*0.5
end do
lonboundsvect(cntx+1)=360.

do i = 2,cntx+1
   lonbound(i-1,1)=lonboundsvect(i-1)
   lonbound(i-1,2)=lonboundsvect(i)
end do

! these are not set accurately!
xrange(1) = minval(lonvect)
xrange(2) = maxval(lonvect)

yrange(1) = minval(latvect)
yrange(2) = maxval(latvect)

!----
! Create the NetCDF file

write(*,*)'Starting annual file'

! First create the annual file
file_to_write_extended = trim(file_to_write)//'_CLASSCTEM_A.nc'
call create_netcdf(totyrs,yrst,realyrst,file_to_write_extended,CTEM,.FALSE.,DOFIRE,DOWETLANDS,COMPETE_LNDUSE,nl,TILED,DOPFTS)

write(*,*)'Done annual file, starting monthly (if needed)'

if (MAKEMONTHLY) then
  ! Then create the monthly file if you are making one:
  file_to_write_extended = trim(file_to_write)//'_CLASSCTEM_M.nc'
  adjustyr=realyrst+(totyrs - monyrs - 1) 
  call create_netcdf(monyrs,yrst,adjustyr,file_to_write_extended,CTEM,.TRUE.,DOFIRE,DOWETLANDS,COMPETE_LNDUSE,nl,TILED,DOPFTS)
end if

deallocate(lonvect)
deallocate(latvect)

deallocate(latboundsvect)
deallocate(lonboundsvect)

deallocate(latbound)
deallocate(lonbound)

end program netcdf_create_bf

!=======================================================================

subroutine create_netcdf(totyrs,yrst,realyrst,file_to_write,CTEM,MONTHFILE,DOFIRE,DOWETLANDS,COMPETE_LNDUSE,nl,TILED,DOPFTS)

use creator_module_bf
use netcdf

implicit none

integer :: i,j
integer,intent(in) :: totyrs
integer,intent(in) :: yrst
integer,intent(in) :: realyrst
character(120),intent(in) :: file_to_write
logical,intent(in) :: CTEM
logical,intent(in) :: MONTHFILE
logical,intent(in) :: DOFIRE
logical,intent(in)  :: DOWETLANDS
logical,intent(in) :: COMPETE_LNDUSE
logical,intent(in) :: TILED
logical,intent(in) :: DOPFTS
integer, intent(in) :: nl

integer :: grpid_ann_ctem
integer :: grpid_ann_class
integer :: grpid_mon_ctem
integer :: grpid_mon_class
integer :: grpid_mon_dist
integer :: grpid_ann_dist
integer :: grpid_ann_wet
integer :: grpid_mon_wet

integer :: grpid_ann_ctem_t
integer :: grpid_mon_ctem_t
integer :: grpid_mon_dist_t
integer :: grpid_ann_dist_t

integer :: tottime

integer, allocatable, dimension(:) :: layersnum
integer, allocatable, dimension(:) :: pftnum
integer, allocatable, dimension(:) :: tilesnum
integer, allocatable, dimension(:) :: timenum

character(8)  :: today
character(10) :: now
character(20) :: yrssince
character(22) :: daysince
integer :: z
character(4) :: zchar
integer :: counter

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

!status = nf90_put_att(ncid,varid,'bounds','lon_bnds')
!if (status/=nf90_noerr) call handle_err(status)

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

!status = nf90_put_att(ncid,varid,'bounds','lat_bnds')
!if (status/=nf90_noerr) call handle_err(status)

!----3

!status = nf90_def_dim(ncid,'bnds',2,bnds)
!if (status/=nf90_noerr) call handle_err(status)

!----3
status = nf90_def_dim(ncid,'tile',ntile,tile)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_def_var(ncid,'tile',nf90_short,tile,varid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','tile')
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

!----7

z=realyrst-1
write (zchar, '(I4)') z
daysince='days since '//zchar//'-01-01'

tottime = totyrs * 12

if (MONTHFILE) then
    status = nf90_def_dim(ncid,'time',tottime,time)
    if (status/=nf90_noerr) call handle_err(status)
else
    status = nf90_def_dim(ncid,'time',totyrs,time)
    if (status/=nf90_noerr) call handle_err(status)
endif

status = nf90_def_var(ncid,'time',nf90_int,time,varid)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'long_name','time')
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units',daysince)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'calendar',"365_day")
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Storage',"chunked")
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Chunksizes',1)
if (status/=nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_Endianness',"little")
if (status/=nf90_noerr) call handle_err(status)

!==============Create the groups====================


if (net4) then

  if (.not. MONTHFILE) then
    status = nf90_def_grp(ncid,'CLASS-Annual', grpid_ann_class)
    if (status /= nf90_noerr) call handle_err(status)
  else 
     status = nf90_def_grp(ncid,'CLASS-Monthly', grpid_mon_class)
     if (status /= nf90_noerr) call handle_err(status)
  end if

  if (CTEM) then

    if (.not. MONTHFILE) then
      status = nf90_def_grp(ncid,'CTEM-Annual GridAvg', grpid_ann_ctem)
      if (status /= nf90_noerr) call handle_err(status)

      if (DOFIRE) then

          status = nf90_def_grp(grpid_ann_ctem,'Annual-Disturbance GridAvg', grpid_ann_dist)
          if (status /= nf90_noerr) call handle_err(status)

      end if

      if (DOWETLANDS) then
          status = nf90_def_grp(grpid_ann_ctem,'Annual-Methane flux GridAvg', grpid_ann_wet)
          if (status /= nf90_noerr) call handle_err(status)
      end if

    else !monthly

      status = nf90_def_grp(ncid,'CTEM-Monthly GridAvg', grpid_mon_ctem)
      if (status /= nf90_noerr) call handle_err(status)

      if (DOFIRE) then
         status = nf90_def_grp(grpid_mon_ctem,'Monthly-Disturbance GridAvg', grpid_mon_dist)
         if (status /= nf90_noerr) call handle_err(status)
      end if

      if (DOWETLANDS) then
         status = nf90_def_grp(grpid_mon_ctem,'Monthly-Methane flux GridAvg', grpid_mon_wet)
         if (status /= nf90_noerr) call handle_err(status)
      end if 

    end if !makemonthly

    if (.not. MONTHFILE) then
        status = nf90_def_grp(ncid,'CTEM-Annual Tiled', grpid_ann_ctem_t)
        if (status /= nf90_noerr) call handle_err(status)

        if (DOFIRE) then
          status = nf90_def_grp(grpid_ann_ctem_t,'Annual-Disturbance Tiled', grpid_ann_dist_t)
          if (status /= nf90_noerr) call handle_err(status)
        end if
  
     else !monthly

        status = nf90_def_grp(ncid,'CTEM-Monthly Tiled', grpid_mon_ctem_t)
        if (status /= nf90_noerr) call handle_err(status)

        if (DOFIRE) then
          status = nf90_def_grp(grpid_mon_ctem_t,'Monthly-Disturbance Tiled', grpid_mon_dist_t)
          if (status /= nf90_noerr) call handle_err(status)
        end if
     end if !makemonthly

  end if  !ctem

else ! netcf3

! No groups allowed in netcdf3 so set all to ncid value
        grpid_ann_ctem=ncid
        grpid_ann_class=ncid
        grpid_mon_ctem=ncid
        grpid_mon_class=ncid
        grpid_mon_dist=ncid
        grpid_ann_dist=ncid
        grpid_ann_wet=ncid   
        grpid_mon_wet=ncid   

        grpid_ann_ctem_t=ncid
        grpid_mon_ctem_t=ncid
        grpid_mon_dist_t=ncid
        grpid_ann_dist_t=ncid

end if

status = nf90_enddef(ncid) 
if (status/=nf90_noerr) call handle_err(status)

!========= Fill in the parameters =======================

allocate(tilesnum(ntile))
allocate(pftnum(ctemnpft))
allocate(layersnum(nl))

if (MONTHFILE) then
    allocate(timenum(tottime))
else
    allocate(timenum(totyrs))
end if

forall (i=1:ntile)
  tilesnum(i) = i
end forall

forall (i=1:ctemnpft)
  pftnum(i) = i
end forall

forall (i=1:nl)
  layersnum(i) = i
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


if (MONTHFILE) then
    counter = 1
    do i=yrst,totyrs
      do j=1,12
        timenum(counter) = i * 365 + monthend(j+1)
        counter = counter + 1
      enddo
    enddo
else
    forall (i=yrst:totyrs)
        timenum(i) = i * 365

    end forall
end if

status = nf90_put_var(ncid,time,timenum)    
if (status/=nf90_noerr) call handle_err(status)

deallocate(tilesnum)
deallocate(pftnum)
deallocate(layersnum)
deallocate(timenum)


! Put in the bounds

!status = nf90_redef(ncid) 
!if (status/=nf90_noerr) call handle_err(status)

!status = nf90_def_var(ncid,'lon_bnds',nf90_float,[lon,bnds],varid)
!if (status/=nf90_noerr) call handle_err(status)

!status = nf90_enddef(ncid) 
!if (status/=nf90_noerr) call handle_err(status)

!status = nf90_put_var(ncid,lon_bnds,lonbound,start=[1,1],count=[cntx+1,2])  !lon bounds
!if (status/=nf90_noerr) call handle_err(status)


!status = nf90_redef(ncid) 
!if (status/=nf90_noerr) call handle_err(status)

!status = nf90_def_var(ncid,'lat_bnds',nf90_float,[lat,bnds],varid)
!if (status/=nf90_noerr) call handle_err(status)

!status = nf90_enddef(ncid) 
!if (status/=nf90_noerr) call handle_err(status)

!status = nf90_put_var(ncid,lat_bnds,latbound,start=[1,1],count=[cnty+1,2])  !lat bounds
!if (status/=nf90_noerr) call handle_err(status)

!================ Now define the variables===========

! CTEM first

 if (CTEM) then

    if (MONTHFILE) then
        ! Annual GRID average:
            call dothevars(3,tottime,grpid_mon_ctem,CTEM_M_VAR_GA,CTEM_M_NAME,CTEM_M_UNIT,numctemvars_m,0,nl)
        ! Fire
            call dothevars(3,tottime,grpid_mon_dist,CTEM_M_D_VAR_GA,CTEM_M_D_NAME,CTEM_M_D_UNIT,nctemdistvars_m,0,nl)
        ! Wetlands
            call dothevars(3,tottime,grpid_mon_wet,CTEM_M_W_VAR,CTEM_M_W_NAME,CTEM_M_W_UNIT,nctemwetvars_m,0,nl)
        if (tiled) then
            ! Annual TILE average:
                call dothevars(4,tottime,grpid_mon_ctem,CTEM_M_VAR_TA,CTEM_M_NAME,CTEM_M_UNIT,numctemvars_m,0,nl)
            ! Fire
                call dothevars(4,tottime,grpid_mon_dist,CTEM_M_D_VAR_TA,CTEM_M_D_NAME,CTEM_M_D_UNIT,nctemdistvars_m,0,nl)
            ! Wetlands
                call dothevars(4,tottime,grpid_mon_wet,CTEM_M_W_T_VAR,CTEM_M_W_NAME,CTEM_M_W_UNIT,nctemwetvars_m,0,nl)
        end if
        if (DOPFTS) then
        ! Annual Per PFT values:
            call dothevars(5,tottime,grpid_mon_ctem,CTEM_M_VAR,CTEM_M_NAME,CTEM_M_UNIT,numctemvars_m,0,nl)
        ! Fire
            call dothevars(5,tottime,grpid_mon_dist,CTEM_M_D_VAR,CTEM_M_D_NAME,CTEM_M_D_UNIT,nctemdistvars_m,0,nl)
        end if
        if (COMPETE_LNDUSE) then
            call dothevars(3,tottime,grpid_mon_ctem,CTEM_M_C_VAR(1),CTEM_M_C_NAME(1),CTEM_M_C_UNIT(1),1,0,nl)
            call dothevars(4,tottime,grpid_mon_ctem,CTEM_M_C_VAR(2:nctemcompvars_m),CTEM_M_C_NAME(2:nctemcompvars_m),CTEM_M_C_UNIT(2:nctemcompvars_m),nctemcompvars_m-1,1,nl)! specialdim = 1 for COMPETE_LNDUSE
        end if
    else ! ANNUAL file
        ! Annual GRID average:
            call dothevars(3,totyrs,grpid_ann_ctem,CTEM_Y_VAR_GA,CTEM_Y_NAME,CTEM_Y_UNIT,numctemvars_a,0,nl)
        ! Fire
            call dothevars(3,totyrs,grpid_ann_dist,CTEM_Y_D_VAR_GA,CTEM_Y_D_NAME,CTEM_Y_D_UNIT,nctemdistvars_a,0,nl)
        ! Wetlands
            call dothevars(3,totyrs,grpid_ann_wet,CTEM_Y_W_VAR,CTEM_Y_W_NAME,CTEM_Y_W_UNIT,nctemwetvars_a,0,nl)
        if (tiled) then
            ! Annual TILE average:
                call dothevars(4,totyrs,grpid_ann_ctem,CTEM_Y_VAR_TA,CTEM_Y_NAME,CTEM_Y_UNIT,numctemvars_a,0,nl)
            ! Fire
                call dothevars(4,totyrs,grpid_ann_dist,CTEM_Y_D_VAR_TA,CTEM_Y_D_NAME,CTEM_Y_D_UNIT,nctemdistvars_a,0,nl)
            ! Wetlands
                call dothevars(4,totyrs,grpid_ann_wet,CTEM_Y_W_T_VAR,CTEM_Y_W_NAME,CTEM_Y_W_UNIT,nctemwetvars_a,0,nl)
        end if
        if (DOPFTS) then
        ! Annual Per PFT values:
            call dothevars(5,totyrs,grpid_ann_ctem,CTEM_Y_VAR,CTEM_Y_NAME,CTEM_Y_UNIT,numctemvars_a,0,nl)
        ! Fire
            call dothevars(5,totyrs,grpid_ann_dist,CTEM_Y_D_VAR,CTEM_Y_D_NAME,CTEM_Y_D_UNIT,nctemdistvars_a,0,nl)
        end if

        if (COMPETE_LNDUSE) then
            call dothevars(3,totyrs,grpid_ann_ctem,CTEM_Y_C_VAR(1),CTEM_Y_C_NAME(1),CTEM_Y_C_UNIT(1),1,0,nl)
            call dothevars(4,totyrs,grpid_ann_ctem,CTEM_Y_C_VAR(2:nctemcompvars_a),CTEM_Y_C_NAME(2:nctemcompvars_a),CTEM_Y_C_UNIT(2:nctemcompvars_a),nctemcompvars_a-1,1,nl)! specialdim = 1 for COMPETE_LNDUSE
        end if

    end  if
 end if

 ! Now do CLASS:
 if (MONTHFILE) then
    call dothevars(3,tottime,grpid_ann_class,CLASS_M_VAR,CLASS_M_NAME,CLASS_M_UNIT,numclasvars_m,0,nl)
    ! Now the soil layer vars
    call dothevars(3,tottime,grpid_ann_class,CLASS_M_S_VAR,CLASS_M_S_NAME,CLASS_M_S_UNIT,nclassoilvars_m,2,nl) ! specialdim = 2 for the CLASS soil layers
 else ! Annual file
    call dothevars(3,totyrs,grpid_mon_class,CLASS_Y_VAR,CLASS_Y_NAME,CLASS_Y_UNIT,numclasvars_a,0,nl)
 end if


!close the netcdf
status = nf90_close(ncid)
if (status/=nf90_noerr) call handle_err(status)

end subroutine create_netcdf

! ==============================================================================================
subroutine dothevars(numdims,tottime,grpid,inarray,namearray,unitarray,sizear,specialdim,nl)

  use netcdf
  use creator_module_bf

  implicit none

  integer, intent(in) :: sizear
  integer, intent(in) :: grpid
  integer, intent(in) :: numdims
  integer, intent(in) :: tottime
  integer, intent(in) :: nl
  real, allocatable, dimension(:,:,:) :: threevar
  real, allocatable, dimension(:,:,:,:) :: fourvar
  real, allocatable, dimension(:,:,:,:,:) :: fivevar
  character(100), intent(in), dimension(sizear) :: inarray
  character(100), intent(in), dimension(sizear) :: namearray
  character(100), intent(in), dimension(sizear) :: unitarray
  integer, intent(in) :: specialdim
  integer :: i

   if (numdims == 3) then
        allocate(threevar(cntx,cnty,tottime))
        threevar=fill_value  !comes in via the module
   else if (numdims == 4) then
        allocate(fourvar(cntx,cnty,ntile,tottime))
        fourvar=fill_value
   else if (numdims == 5) then
        allocate(fivevar(cntx,cnty,ctemnpft,ntile,tottime))
        fivevar=fill_value
   else
        write(*,*)'Incorrect number of dims!'
    end if

  do i=lbound(inarray,1), ubound(inarray,1)

   status = nf90_redef(grpid)
   if (status/=nf90_noerr) call handle_err(status)

   if (numdims == 3) then
        status = nf90_def_var(grpid,trim(inarray(i)),nf90_float,[lon,lat,time],varid)
        if (status/=nf90_noerr) call handle_err(status)
   else if (numdims == 4) then
        if (specialdim == 1) then  !COMPETE_LNDUSE is TRUE
            status = nf90_def_var(grpid,trim(inarray(i)),nf90_float,[lon,lat,pft,time],varid)
            if (status/=nf90_noerr) call handle_err(status)
        else if (specialdim == 2) then ! CLASS SOIL
            status = nf90_def_var(grpid,trim(inarray(i)),nf90_float,[lon,lat,layer,time],varid)
            if (status/=nf90_noerr) call handle_err(status)
        else
            status = nf90_def_var(grpid,trim(inarray(i)),nf90_float,[lon,lat,tile,time],varid)
            if (status/=nf90_noerr) call handle_err(status)
        end if
   else if (numdims == 5) then
        status = nf90_def_var(grpid,trim(inarray(i)),nf90_float,[lon,lat,pft,tile,time],varid)
        if (status/=nf90_noerr) call handle_err(status)
   end if

   status = nf90_enddef(grpid)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_redef(grpid)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_att(grpid,varid,'long_name',trim(namearray(i)))
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_att(grpid,varid,'units',trim(unitarray(i)))
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_att(grpid,varid,'_FillValue',fill_value)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_att(grpid,varid,'missing_value',fill_value)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_att(grpid,varid,'_Storage',"chunked")
   if (status/=nf90_noerr) call handle_err(status)

   if (numdims == 3) then
        status = nf90_put_att(grpid,varid,'_Chunksizes',[cnty,cntx,1])
        if (status/=nf90_noerr) call handle_err(status)
   else if (numdims == 4) then
        if (specialdim == 1) then  !COMPETE_LNDUSE is TRUE
            status = nf90_put_att(grpid,varid,'_Chunksizes',[cnty,cntx,ctemnpft,1])
            if (status/=nf90_noerr) call handle_err(status)
        else if (specialdim == 2) then ! CLASS SOIL
            status = nf90_put_att(grpid,varid,'_Chunksizes',[cnty,cntx,nl,1])
            if (status/=nf90_noerr) call handle_err(status)
        else
            status = nf90_put_att(grpid,varid,'_Chunksizes',[cnty,cntx,ntile,1])
            if (status/=nf90_noerr) call handle_err(status)
        end if
   else if (numdims == 5) then
        status = nf90_put_att(grpid,varid,'_Chunksizes',[cnty,cntx,ctemnpft,ntile,1])
        if (status/=nf90_noerr) call handle_err(status)
   end if

   status = nf90_put_att(grpid,varid,'_DeflateLevel',1)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_enddef(grpid)
   if (status/=nf90_noerr) call handle_err(status)

   if (numdims == 3) then
        status = nf90_put_var(grpid,varid,threevar,start=[1,1,1],count=[cntx,cnty,tottime])
        if (status/=nf90_noerr) call handle_err(status)
   else if (numdims == 4) then
        if (specialdim == 1) then  !COMPETE_LNDUSE is TRUE
            status = nf90_put_var(grpid,varid,fourvar,start=[1,1,1,1],count=[cntx,cnty,ctemnpft,tottime])
            if (status/=nf90_noerr) call handle_err(status)
        else if (specialdim == 2) then ! CLASS SOIL
            status = nf90_put_var(grpid,varid,fourvar,start=[1,1,1,1],count=[cntx,cnty,nl,tottime])
            if (status/=nf90_noerr) call handle_err(status)
        else
            status = nf90_put_var(grpid,varid,fourvar,start=[1,1,1,1],count=[cntx,cnty,ntile,tottime])
            if (status/=nf90_noerr) call handle_err(status)
        end if
   else if (numdims == 5) then
        status = nf90_put_var(grpid,varid,fivevar,start=[1,1,1,1,1],count=[cntx,cnty,ctemnpft,ntile,tottime])
        if (status/=nf90_noerr) call handle_err(status)
    end if

  end do

   if (numdims == 3) then
        deallocate(threevar)
   else if (numdims == 4) then
        deallocate(fourvar)
   else if (numdims == 5) then
        deallocate(fivevar)
    end if

end subroutine dothevars
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


