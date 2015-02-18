program netcdf_writer_fast

!Program to write CLASS and CTEM output variables to NETCDF files

use netcdf
use creator_module_bf

implicit none

! Annual COMPOSITE 
real, allocatable, dimension(:,:) :: ctem_a
real, allocatable, dimension(:,:) :: ctem_d_a
real, allocatable, dimension(:,:) :: ctem_w_a  !Rudra

! Annual MOSAIC
real, allocatable, dimension(:,:,:) :: ctem_a_mos
real, allocatable, dimension(:,:,:) :: ctem_d_a_mos

! Annual Competition/luc
character(1), allocatable, dimension(:,:) :: pftexist
real, allocatable, dimension(:,:) :: pftf
real, allocatable, dimension(:) :: pft_tot

! Monthly COMPOSITE
real, allocatable, dimension(:,:,:) :: ctem_m
real, allocatable, dimension(:,:,:) :: ctem_d_m
real, allocatable, dimension(:,:,:) :: ctem_w_m  !Rudra

! Monthly MOSAIC
real, allocatable, dimension(:,:,:,:) :: ctem_m_mos
real, allocatable, dimension(:,:,:,:) :: ctem_d_m_mos

! Monthly Competition/luc
character(1), allocatable, dimension(:,:,:) :: mpftexist
real, allocatable, dimension(:,:,:) :: mpftf
real, allocatable, dimension(:,:) :: mpft_tot

! Monthly variables for CLASS
real, allocatable, dimension(:,:,:) :: class_m
real, allocatable, dimension(:,:,:,:) :: class_s_m

! Annual variables for CLASS
real, allocatable, dimension(:,:) :: class_a

integer :: y,dummy_year,dummy_month
real :: dummy_var
integer, dimension(1) :: lyear
integer :: i,m,xlon,ylat,yrst,h,l,cell
integer :: totyrs
integer :: monyrs
integer :: yrin
integer :: yr_now
integer :: mo
integer :: tilnum
integer :: dummynum
integer :: var_id,v,p,grpid
character(120) :: ARGBUFF 
character(120) :: file_to_write
character(140) :: file_to_write_extended
character(120) :: folder
character(120) :: long_path
character(120) :: jobfile
character(512) :: command
character(160) :: infile
character(160) :: cellsfile
character(4) :: dummy
character(6) :: lon_in, lat_in
logical :: CTEM
logical :: MAKEMONTHLY
logical :: DOFIRE
logical :: DOWETLANDS   !Rudra
logical :: MOSAIC
logical :: PARALLELRUN
logical :: COMPETE_LNDUSE
integer :: realyrst
integer :: tailyrs
character(1) :: tic
character(len=8) :: fmt ! format descriptor
character(len=10) :: x1
integer :: io_set
logical :: lexist

integer, allocatable, dimension(:) :: tiles_to_write

real, allocatable, dimension(:) :: tmp
real, allocatable, dimension(:) :: tmpd

real, allocatable, dimension(:,:) :: tmpa
real, allocatable, dimension(:,:,:) :: tmpm

namelist /joboptions/ &
  PARALLELRUN,        &
  CTEM, 	      &
  MAKEMONTHLY,        &
  DOFIRE,             &
  DOWETLANDS,         &   
  MOSAIC,             &
  COMPETE_LNDUSE,     &
  totyrs,             &
  monyrs,             &
  yrst,               &
  realyrst ,          &
  long_path ,         &
  file_to_write                
!----

tic=char(39)

allocate(tiles_to_write(ntile))

! Open the grid cells file
 call getarg(1,cellsfile)
 
 call getarg(2,jobfile)

 open(10,file=jobfile,status='old')

 read(10,nml = joboptions)

 close(10)

!********************************************************************
! First do all of the annual writes:

! open the list of coordinates
 open(11,FILE=cellsfile,status='old')

! Open the netcdf file for writing
! Annual:
file_to_write_extended = trim(file_to_write)//'_CLASSCTEM_A.nc'

write(*,*)'annual: writing to',file_to_write_extended

status = nf90_open(file_to_write_extended,nf90_write,ncid)
if (status /= nf90_noerr) call handle_major_err(status)

do cell = 1,num_land_cells

    write(*,*)'annual:writing ',cell,' of ',num_land_cells
    read(11,*)lon_in,lat_in

    ARGBUFF = trim(lon_in)//'_'//trim(lat_in)

    call parsecoords(ARGBUFF,bounds)

    folder=trim(long_path)//'/'//trim(ARGBUFF)//'/'

    ! get the coordinates
    xlon=bounds(1)
    ylat=bounds(3)
    
    tiles_to_write = -1 ! initialize to -1 now.

!==CLASS ANNUAL =====================================================

! OF1Y_G 

if (.NOT. net4) then
   grpid=ncid
end if

inquire(file=trim(folder)//trim(ARGBUFF)//'.OF1Y_G',exist=lexist)
if (lexist) then

    OPEN(83,FILE=trim(folder)//trim(ARGBUFF)//'.OF1Y_G',status='old',form='formatted') ! YEARLY OUTPUT FOR CLASS

    !  first throw out header
       do h = 1,5
    	read(83,*,iostat=io_set)
    	if (io_set .ne. 0) then
    	    write(*,*)'Missing/truncated file',trim(ARGBUFF)//'.OF1Y_G'
    	    close(83)
    	    goto 111
    	end if     
       end do

    !Allocate arrays
    allocate(class_a(numclasvars_a,totyrs))
    
       do y = 1,totyrs
    
         read(83,*,iostat=io_set)dummy_year,class_a(1:numclasvars_a,y)
         if (io_set .ne. 0 .and. y < totyrs) then
    	    write(*,*)'Missing/truncated file',trim(ARGBUFF)//'.OF1Y_G'         
            close(83)
            goto 111
         end if   
       end do

     close(83)


! Find the group id of annual Class
if (net4) then
  status = nf90_inq_ncid(ncid,'CLASS-Annual', grpid)
  if (status /= nf90_noerr) call handle_err(status)
end if

do v = 1,numclasvars_a ! begin vars loop

   status = nf90_inq_varid(grpid,trim(CLASS_Y_VAR(v)), var_id)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_var(grpid,var_id,class_a(v,:),start=[xlon,ylat,yrst],count=[1,1,totyrs])
   if (status/=nf90_noerr) call handle_err(status)

end do ! vars loop

111 continue !error thrown by file

! deallocate arrays
deallocate(class_a)

end if !lexist

!===CTEM COMPETITION ANNUAL FILES======================================================

IF (CTEM) THEN

! Do the COMPETE_LNDUSE ones first since they are independent of mosaic vs. composite.

if (COMPETE_LNDUSE) then

inquire(file=trim(folder)//trim(ARGBUFF)//'.CT07Y_GM',exist=lexist)

if (lexist) then
    OPEN(92,FILE=trim(folder)//trim(ARGBUFF)//'.CT07Y_GM',status='old',form='formatted') ! YEARLY PFR FRAC OUTPUT FOR CTEM

    allocate(pftexist(ctemnpft,totyrs))
    allocate(tmpa(ctemnpft,totyrs))
    allocate(pftf(ctemnpft,totyrs))
    allocate(pft_tot(totyrs))

   do h = 1,6
	read(92,*,iostat=io_set)
	if (io_set .ne. 0) then
      write(*,*)'Missing/truncated file',trim(ARGBUFF)//'.CT07Y_GM'
      close(92)
      goto 112
    end if     
   end do

   do y = 1,totyrs

    read(92,*,iostat=io_set)dummy_year,pftf(1:ctemnpft,y),pft_tot(y),dummy_var,pftexist(1:ctemnpft,y)
    if (io_set .ne. 0 .and. y < totyrs) then
       write(*,*)'Missing/truncated file',trim(ARGBUFF)//'.CT07Y_GM'    
       close(92)
       goto 112
    end if   

   end do

   close(92)

if (net4) then
  status = nf90_inq_ncid(ncid,'CTEM-Annual GridAvg', grpid)
  if (status /= nf90_noerr) call handle_err(status)
end if

! First write the pft fractions and pft exist
do p = 1,ctemnpft ! begin pft loop

   status = nf90_inq_varid(grpid,trim(CTEM_Y_C_VAR(2)), var_id)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_var(grpid,var_id,pftf(p,:),start=[xlon,ylat,p,yrst],count=[1,1,1,totyrs])
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_inq_varid(grpid,trim(CTEM_Y_C_VAR(3)), var_id)
   if (status/=nf90_noerr) call handle_err(status)

     do y = 1,totyrs  
      if (pftexist(p,y) == 'F') then
         tmpa(p,y)=0.0
      else
         tmpa(p,y)=1.0       
      end if
     end do

   status = nf90_put_var(grpid,var_id,tmpa(p,:),start=[xlon,ylat,p,yrst],count=[1,1,1,totyrs])
   if (status/=nf90_noerr) call handle_err(status)

end do ! pfts loop

! Now write the total plant fraction
   status = nf90_inq_varid(grpid,trim(CTEM_Y_C_VAR(1)), var_id)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_var(grpid,var_id,pft_tot(:),start=[xlon,ylat,yrst],count=[1,1,totyrs])
   if (status/=nf90_noerr) call handle_err(status)

112 continue !error thrown by file

deallocate(pftf)
deallocate(pft_tot)
deallocate(pftexist)
deallocate(tmpa)

end if !lexist

end if  !compete/landuse

!====CTEM Grid Average ANNUAL FILES=========================================

! Allocate the arrays in preparation for CT01Y_G, CT06Y_G, CT01Y_GM

  ! Annual CTEM
  if (MOSAIC) then
  infile=trim(folder)//trim(ARGBUFF)//'.CT01Y_M'
  else
  infile=trim(folder)//trim(ARGBUFF)//'.CT01Y_G'
  end if

inquire(file=infile,exist=lexist)
if (lexist) then

  command='sed -n '//tic//'/GRDAV/p'//tic//' '//trim(infile)//' > tmp_a.dat'  !this removes any header!

  call system(command)

  ! Check to make sure the file was created
  !inquire(file='tmp_a.dat',exist=lexist)
  
  !if (lexist) then
  OPEN(75,FILE='tmp_a.dat',status='old',form='formatted') ! YEARLY OUTPUT FOR CTEM
  
  !allocate the size of the arrays for the output data
  allocate(ctem_a(numctemvars_a,totyrs))

  do y = 1,totyrs
     read(75,*,iostat=io_set)dummy_year,ctem_a(1:numctemvars_a,y)
     if (io_set .ne. 0 .and. y < totyrs) then
        write(*,*)'Missing/truncated file',trim(ARGBUFF)//'.CT01Y_G/M'
        close(75)    
        goto 113
     end if   
  end do
  
  close(75)
    
  if (net4) then
    status = nf90_inq_ncid(ncid,'CTEM-Annual GridAvg', grpid)
    if (status /= nf90_noerr) call handle_err(status)
  end if

  do v = 1,numctemvars_a ! begin vars loop

   status = nf90_inq_varid(grpid,trim(CTEM_Y_VAR_GA(v)), var_id)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_var(grpid,var_id,ctem_a(v,:),start=[xlon,ylat,yrst],count=[1,1,totyrs])
   if (status/=nf90_noerr) call handle_err(status)

  end do ! vars loop

113 continue !error thrown by file

  ! deallocate arrays
  deallocate(ctem_a)
  
! endif !lexist (tmp_a file)
end if !lexist

! -----Fire --------------

if (DOFIRE) then
  
    if (MOSAIC) then
     infile=trim(folder)//trim(ARGBUFF)//'.CT06Y_M'
    else
     infile=trim(folder)//trim(ARGBUFF)//'.CT06Y_G'
    end if

    inquire(file=infile,exist=lexist)
    if (lexist) then

     command='sed -n '//tic//'/GRDAV/p'//tic//' '//trim(infile)//' > tmp_a_d.dat'
     call system(command)
     
     OPEN(85,FILE='tmp_a_d.dat',status='old',form='formatted') 

  allocate(ctem_d_a(nctemdistvars_a,totyrs))
  
  do y = 1,totyrs
     read(85,*,iostat=io_set)dummy_year,ctem_d_a(1:nctemdistvars_a,y)
     if (io_set .ne. 0 .and. y < totyrs) then
         write(*,*)'Missing/truncated file',trim(ARGBUFF)//'.CT06Y_G/M'
         close(85)
         goto 114
     end if   
  end do      
  
  close(85)
  
   if (net4) then
    status = nf90_inq_ncid(ncid,'Annual-Disturbance GridAvg', grpid)
    if (status /= nf90_noerr) call handle_err(status)
  end if

  do v = 1,nctemdistvars_a ! begin vars loop

    status = nf90_inq_varid(grpid,trim(CTEM_Y_D_VAR_GA(v)), var_id)
    if (status/=nf90_noerr) call handle_err(status)

    status = nf90_put_var(grpid,var_id,ctem_d_a(v,:),start=[xlon,ylat,yrst],count=[1,1,totyrs])
    if (status/=nf90_noerr) call handle_err(status)

  end do ! vars loop

114 continue !error thrown by file

 deallocate(ctem_d_a)

 end if !lexist
 
endif ! if fire.

! -----Wetlands --------------
if (DOWETLANDS) then 

  inquire(file=trim(folder)//trim(ARGBUFF)//'.CT08Y_G',exist=lexist)
  if (lexist) then

  OPEN(900,FILE=trim(folder)//trim(ARGBUFF)//'.CT08Y_G',status='old',form='formatted') ! YEARLY CTEM wetlands CH4 emission

  allocate(ctem_w_a(nctemwetvars_a,totyrs))         
  !  first throw out header
      do h = 1,6
        read(900,*,iostat=io_set)
        if (io_set .ne. 0) then
          write(*,*)'Missing/truncated file ',trim(ARGBUFF)//'.CT08Y_G/M'
          close(900)
          goto 115
        end if     
      end do
  do y = 1,totyrs   
      read(900,*,iostat=io_set) dummy_year,ctem_w_a(1:nctemwetvars_a,y)
      if (io_set .ne. 0 .and. y < totyrs) then
         write(*,*)'Missing/truncated file ',trim(ARGBUFF)//'.CT08Y_G/M'      
         close(900)
         goto 115
      end if   

  end do    

  close(900)
  
   if (net4) then
    status = nf90_inq_ncid(ncid,'Annual-Methane flux GridAvg', grpid)
    if (status /= nf90_noerr) call handle_err(status)
   end if
   
  do v = 1,nctemwetvars_a ! begin vars loop

    status = nf90_inq_varid(grpid,trim(CTEM_Y_W_VAR(v)), var_id)
    if (status/=nf90_noerr) call handle_err(status)

    status = nf90_put_var(grpid,var_id,ctem_w_a(v,:),start=[xlon,ylat,yrst],count=[1,1,totyrs])
    if (status/=nf90_noerr) call handle_err(status)

  end do ! vars loop

115 continue !error thrown by file

  deallocate(ctem_w_a)
  
  end if !lexist
  
end if  !if wetlands. 

!now per PFT/Tile MODE ====================================================

! Read in from the ascii file

  if (MOSAIC) then
  infile=trim(folder)//trim(ARGBUFF)//'.CT01Y_M'
  else
  infile=trim(folder)//trim(ARGBUFF)//'.CT01Y_G'
  end if

  inquire(file=infile,exist=lexist)
  if (lexist) then

! Allocate the size of the arrays for the output data
  allocate(ctem_a_mos(numctemvars_a,ntile,totyrs))

   ! Make a file of the PFT info, removing the grdav.
   command='sed '//tic//'/GRDAV/d'//tic//' '//trim(infile)//' > tmp_a_p.dat'
   call system(command)
   OPEN(751,FILE='tmp_a_p.dat',status='old',form='formatted') ! YEARLY OUTPUT FOR CTEM

!  first throw out header
   do h = 1,6
        read(751,*,iostat=io_set)
        if (io_set .ne. 0) then
          write(*,*)'Missing/truncated file ',trim(ARGBUFF)//'.CT01Y_G/M tiled'
          close(751)
          goto 116
        end if             
   end do

  allocate(tmp(numctemvars_a))

! We have to keep track of the year that is read in as it is the only way we know that we are done the tiles for a gridcell.
    yrin=realyrst

    do y = 1,totyrs

      yr_now = realyrst + y - 1

      do while (yrin == yr_now)

          read(751,*,iostat=io_set,end=90) yrin,tmp(1:numctemvars_a),dummy,dummynum,dummy,tilnum
          if (io_set .ne. 0 .and. y < totyrs) then
            write(*,*)'Missing/truncated file ',trim(ARGBUFF)//'.CT01Y_G/M tiled'          
            close(751)
            goto 116
          end if   

          if (yrin == realyrst) then
             tiles_to_write = cshift(tiles_to_write,1)
             tiles_to_write(1) = tilnum
          end if     

        if (yrin == yr_now) then
          ! Assign that just read in to its vars 
          do v = 1,numctemvars_a ! begin vars loop
            ctem_a_mos(v,tilnum,y)=tmp(v)
          end do
        else
            backspace(751)
        end if
       end do ! while loop
     end do ! y loop

90     continue
  
! done with ascii file, close it.
  close(751)
  deallocate(tmp)

if (net4) then
  status = nf90_inq_ncid(ncid,'CTEM-Annual Tiled', grpid)
  if (status /= nf90_noerr) call handle_err(status)
end if

 do l=1,ntile
  do v = 1,numctemvars_a ! begin vars loop

   if (tiles_to_write(l) .ne. -1) then
 
     status = nf90_inq_varid(grpid,trim(CTEM_Y_VAR(v)), var_id)
     if (status/=nf90_noerr) call handle_err(status)
 
     status = nf90_put_var(grpid,var_id,ctem_a_mos(v,tiles_to_write(l),:),start=[xlon,ylat,tiles_to_write(l),yrst],count=[1,1,1,totyrs])
     if (status/=nf90_noerr) call handle_err(status)
 
   end if

  end do ! vars loop
 end do !tiles loop

116 continue !error thrown by file

! deallocate arrays
deallocate(ctem_a_mos)

end if !lexist

! --Fire----------------------------------------------
if (DOFIRE) then
 
    if (MOSAIC) then
     infile=trim(folder)//trim(ARGBUFF)//'.CT06Y_M'
    else
     infile=trim(folder)//trim(ARGBUFF)//'.CT06Y_G'
    end if

  inquire(file=infile,exist=lexist)
  if (lexist) then

     allocate(ctem_d_a_mos(nctemdistvars_a,ntile,totyrs))
     allocate(tmpd(nctemdistvars_a))

     command='sed '//tic//'/GRDAV/d'//tic//' '//trim(infile)//' > tmp_a_p_d.dat'
     call system(command)
     OPEN(851,FILE='tmp_a_p_d.dat',status='old',form='formatted')

 !  first throw out header
   do h = 1,6
	read(851,*,iostat=io_set)
    if (io_set .ne. 0) then
      write(*,*)'Missing/truncated file ',trim(ARGBUFF)//'.CT06Y_G/M tiled'
      close(851)
      goto 117
    end if             
   end do

   yrin=realyrst
   do y = 1,totyrs
     yr_now = realyrst + y - 1
     do while (yrin == yr_now)

       read(851,*,iostat=io_set,end=91) yrin,tmpd(1:nctemdistvars_a),dummy,dummynum,dummy,tilnum
       if (io_set .ne. 0 .and. y < totyrs) then
          write(*,*)'Missing/truncated file ',trim(ARGBUFF)//'.CT06Y_G/M tiled'       
          close(851)
          goto 117
       end if   
 
       if (yrin == yr_now) then
         do v = 1,nctemdistvars_a ! begin vars loop
           ctem_d_a_mos(v,tilnum,y)=tmpd(v)
         end do
       else
         backspace(851)  ! go back one line in the file
       end if
     end do ! while loop
   end do ! y loop
91     continue

  close(851)
  deallocate(tmpd)

  if (net4) then
    status = nf90_inq_ncid(ncid,'Annual-Disturbance Tiled', grpid)
    if (status /= nf90_noerr) call handle_err(status)
  end if

  do l=1,ntile
   do v = 1,nctemdistvars_a ! begin vars loop

    if (tiles_to_write(l) .ne. -1) then
      status = nf90_inq_varid(grpid,trim(CTEM_Y_D_VAR(v)), var_id)
      if (status/=nf90_noerr) call handle_err(status)
 
      status = nf90_put_var(grpid,var_id,ctem_d_a_mos(v,tiles_to_write(l),:),start=[xlon,ylat,tiles_to_write(l),yrst],count=[1,1,1,totyrs])
      if (status/=nf90_noerr) call handle_err(status)
    end if

   end do ! vars loop
  end do !tiles loop

117 continue !error thrown by file

  deallocate(ctem_d_a_mos)
  
  end if !lexist
  
endif !if fire.

end if !if CTEM.

end do !loop over the gridcells

!!close the annual netcdf
status = nf90_close(ncid)
if (status/=nf90_noerr) call handle_err(status)

! Close the gridcells file
 close(11)

! *******************************************************************
!=====================CLASS_MONTHLY files============================
! *******************************************************************


if (MAKEMONTHLY) then

! Open the netcdf file for writing

!! Monthly file:
file_to_write_extended = trim(file_to_write)//'_CLASSCTEM_M.nc'
write(*,*)'monthly: writing to',file_to_write_extended
status = nf90_open(file_to_write_extended,nf90_write,ncid_m)
if (status /= nf90_noerr) call handle_major_err(status)

!! open the list of coordinates
 open(12,FILE=cellsfile,status='old')

do cell = 1,num_land_cells

    write(*,*)'monthly: starting',cell,' of ',num_land_cells
    read(12,*)lon_in,lat_in
    ARGBUFF = trim(lon_in)//'_'//trim(lat_in)
    call parsecoords(ARGBUFF,bounds)
    folder=trim(long_path)//'/'//trim(ARGBUFF)//'/'

    ! get the coordinates
    xlon=bounds(1)
    ylat=bounds(3)

    tiles_to_write = -1 ! initialize to -1 now.


if (.NOT. net4) then
   grpid=ncid_m
end if

! Do OF1M_G first
if (net4) then
  status = nf90_inq_ncid(ncid_m,'CLASS-Monthly',grpid) 
  if (status /= nf90_noerr) call handle_err(status)
end if

  inquire(file=trim(folder)//trim(ARGBUFF)//'.OF1M_G',exist=lexist)
  if (lexist) then

! Read in the ascii file
   OPEN(81,FILE=trim(folder)//trim(ARGBUFF)//'.OF1M_G',status='old',form='formatted') ! MONTHLY OUTPUT FOR CLASS

    !Allocate Arrays
    allocate(class_m(numclasvars_m,monyrs,12))

 !  first throw out header (5 lines)
   do h = 1,5
	read(81,*,iostat=io_set)
	if (io_set .ne. 0) then
      write(*,*)'Missing/truncated file ',trim(ARGBUFF)//'OF1M_G'
      close(81)
      goto 118
    end if             
   end do
   
   do y = 1,monyrs
    do m=1,12
     read(81,*,iostat=io_set)dummy_month,dummy_year,class_m(1:numclasvars_m,y,m)
     if (io_set .ne. 0 .and. y < monyrs .and. m < 12) then
         write(*,*)'Missing/truncated file ',trim(ARGBUFF)//'OF1M_G'     
         close(81)
         goto 118
     end if   

    end do
   end do

   close(81)

! Write the inputs to the netcdf
do m=1,12  !Begin Month Loop
  do v = 1,numclasvars_m

   status = nf90_inq_varid(grpid,trim(CLASS_M_VAR(v)), var_id)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_var(grpid,var_id,class_m(v,:,m),start=[xlon,ylat,m,yrst],count=[1,1,1,monyrs])
   if (status/=nf90_noerr) call handle_err(status)

  end do !Var loop
end do   !End month loop

118 continue !error thrown by file

! deallocate arrays
deallocate(class_m)

end if !lexist

!=======================================================================
! Do OF2M_G

  inquire(file=trim(folder)//trim(ARGBUFF)//'.OF2M_G',exist=lexist)
  if (lexist) then

   OPEN(82,FILE=trim(folder)//trim(ARGBUFF)//'.OF2M_G',status='old',form='formatted')

   !Allocate Arrays
   allocate(class_s_m(nclassoilvars_m,nl,monyrs,12))

   !  first throw out header (5 lines)
   do h = 1,5
	read(82,*,iostat=io_set)
	if (io_set .ne. 0) then
      write(*,*)'Missing/truncated file ',trim(ARGBUFF)//'OF2M_G'
      close(82)
      goto 119
    end if             	
   end do
        
   do y = 1,monyrs
    do m=1,12
     read(82,*,iostat=io_set)dummy_month,dummy_year,class_s_m(1,1,y,m),class_s_m(2,1,y,m),class_s_m(3,1,y,m),class_s_m(1,2,y,m),class_s_m(2,2,y,m),class_s_m(3,2,y,m),class_s_m(1,3,y,m),class_s_m(2,3,y,m),class_s_m(3,3,y,m)
     if (io_set .ne. 0 .and. y < monyrs .and. m < 12) then
         write(*,*)'Missing/truncated file ',trim(ARGBUFF)//'OF2M_G'     
         close(82)
         goto 119
     end if   

    end do
   end do

   close(82)

!----
! Write to netcdf file
do m=1,12   !begin month loop
 do l=1,nl   ! begin soil layer loop
  do v = 1,nclassoilvars_m ! begin vars loop

   status = nf90_inq_varid(grpid,trim(CLASS_M_S_VAR(v)), var_id)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_var(grpid,var_id,class_s_m(v,l,:,m),start=[xlon,ylat,l,m,yrst],count=[1,1,1,1,monyrs])
   if (status/=nf90_noerr) call handle_err(status)

  end do ! vars loop
 end do ! soil layer loop
end do  !End Month Loop

119 continue !error thrown by file

! deallocate arrays
deallocate(class_s_m)

end if !lexist

! %%%%%%Now MONTHLY CTEM outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IF (CTEM) THEN

! Start with LNDUSECOMPETE vars
if (COMPETE_LNDUSE) then

  inquire(file=trim(folder)//trim(ARGBUFF)//'.CT07M_GM',exist=lexist)
  if (lexist) then

   OPEN(91,FILE=trim(folder)//trim(ARGBUFF)//'.CT07M_GM',status='old',form='formatted') ! MONTHLY PFT FRAC OUTPUT FOR CTEM

    allocate(mpftexist(ctemnpft,monyrs,12))
    allocate(tmpm(ctemnpft,monyrs,12))
    allocate(mpftf(ctemnpft,monyrs,12))
    allocate(mpft_tot(monyrs,12))

   !  first throw out header (6 lines)
   do h = 1,6
	read(91,*,iostat=io_set)
	if (io_set .ne. 0) then
      write(*,*)'Missing/truncated file ',trim(ARGBUFF)//'.CT07M_GM'
      close(91)
      goto 120
    end if   
   end do
   
   do y = 1,monyrs
    do m=1,12
      read(91,*,iostat=io_set)dummy_month,dummy_year,mpftf(1:ctemnpft,y,m),mpft_tot(y,m),dummy_var,mpftexist(1:ctemnpft,y,m)
      if (io_set .ne. 0 .and. y < monyrs .and. m < 12) then
         write(*,*)'Missing/truncated file ',trim(ARGBUFF)//'.CT07M_GM'
         close(91)
         goto 120
      end if   
    end do
   end do

   close(91)

if (net4) then
  status = nf90_inq_ncid(ncid_m,'CTEM-Monthly GridAvg', grpid)
  if (status /= nf90_noerr) call handle_err(status)
end if

! First write the pft fractions and pftexist
 do m=1,12  !Begin Month Loop
  do p = 1,ctemnpft ! begin pft loop

   status = nf90_inq_varid(grpid,trim(CTEM_M_C_VAR(2)), var_id)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_var(grpid,var_id,mpftf(p,:,m),start=[xlon,ylat,p,m,yrst],count=[1,1,1,1,monyrs])
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_inq_varid(grpid,trim(CTEM_M_C_VAR(3)), var_id)
   if (status/=nf90_noerr) call handle_err(status)

     do y = 1,monyrs 
      if (mpftexist(p,y,m) == 'F') then
         tmpm(p,y,m)=0.0
      else
         tmpm(p,y,m)=1.0       
      end if
     end do
  
   status = nf90_put_var(grpid,var_id,tmpm(p,:,m),start=[xlon,ylat,p,m,yrst],count=[1,1,1,1,monyrs])
   if (status/=nf90_noerr) call handle_err(status)


  end do ! pfts loop
 end do !months

 ! Now write the total plant fraction
 do m=1,12  !Begin Month Loop
   status = nf90_inq_varid(grpid,trim(CTEM_M_C_VAR(1)), var_id)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_var(grpid,var_id,mpft_tot(:,m),start=[xlon,ylat,m,yrst],count=[1,1,1,monyrs])
   if (status/=nf90_noerr) call handle_err(status)
 end do

120 continue !error thrown by file

deallocate(mpftexist)
deallocate(mpftf)
deallocate(tmpm)
deallocate(mpft_tot)

end if !lexist

end if !compete/lnduse

!==============================Start Monthly Grid Average Vals CTEM=============================


if (MOSAIC) then
   infile=trim(folder)//trim(ARGBUFF)//'.CT01M_M'
else
   infile=trim(folder)//trim(ARGBUFF)//'.CT01M_G'
end if

inquire(file=infile,exist=lexist)
if (lexist) then

command='sed -n '//tic//'/GRDAV/p'//tic//' '//trim(infile)//' > tmp_m.dat'
call system(command)

! Open the newly trimmed file
OPEN(74,FILE='tmp_m.dat',status='old',form='formatted') ! MONTHLY OUTPUT FOR CTEM

!Allocate Arrays
allocate(ctem_m(numctemvars_m,monyrs,12))

!---Read in Variables
   do y = 1,monyrs
    do m=1,12
     read(74,*,iostat=io_set)dummy_month,dummy_year,ctem_m(1:numctemvars_m,y,m)
     if (io_set .ne. 0 .and. y < monyrs .and. m < 12) then
         write(*,*)'Missing/truncated file ',trim(ARGBUFF)//'.CT01M_G/M'
         close(74)
         goto 121
     end if   
    end do !m loop
   end do !y loop
 close(74)

! MONTHLY
if (net4) then
  status = nf90_inq_ncid(ncid_m,'CTEM-Monthly GridAvg', grpid)
  if (status /= nf90_noerr) call handle_err(status)
end if

 do v = 1,numctemvars_m ! begin vars loop
  do m=1,12  !Begin Month Loop
   status = nf90_inq_varid(grpid,trim(CTEM_M_VAR_GA(v)), var_id)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_var(grpid,var_id,ctem_m(v,:,m),start=[xlon,ylat,m,yrst],count=[1,1,1,monyrs])
   if (status/=nf90_noerr) call handle_err(status)

  end do !months
 end do ! vars loop

121 continue !error thrown by file

! deallocate arrays

deallocate(ctem_m)

end if

if (DOFIRE) then
  
      if (MOSAIC) then
        infile=trim(folder)//trim(ARGBUFF)//'.CT06M_M'
      else
        infile=trim(folder)//trim(ARGBUFF)//'.CT06M_G'
      end if

inquire(file=infile,exist=lexist)
if (lexist) then

     command='sed -n '//tic//'/GRDAV/p'//tic//' '//trim(infile)//' > tmp_m_d.dat'
     call system(command)
     OPEN(84,FILE='tmp_m_d.dat',status='old',form='formatted') 

  allocate(ctem_d_m(nctemdistvars_m,monyrs,12))

!---Read in Variables
    do y = 1,monyrs
     do m=1,12
      read(84,*,iostat=io_set)dummy_month,dummy_year,ctem_d_m(1:nctemdistvars_m,y,m)
      if (io_set .ne. 0 .and. y < monyrs .and. m < 12) then
         write(*,*)'Missing/truncated file ',trim(ARGBUFF)//'.CT06M_G/M'
         close(84)
         goto 122
      end if   
     end do !m loop
    end do !y loop

    close(84)

 if (net4) then
   status = nf90_inq_ncid(ncid_m,'Monthly-Disturbance GridAvg', grpid)
   if (status /= nf90_noerr) call handle_err(status)
 end if

 do v = 1,nctemdistvars_m ! begin vars loop
  do m=1,12  !Begin Month Loop
  
   status = nf90_inq_varid(grpid,trim(CTEM_M_D_VAR_GA(v)), var_id)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_var(grpid,var_id,ctem_d_m(v,:,m),start=[xlon,ylat,m,yrst],count=[1,1,1,monyrs])
   if (status/=nf90_noerr) call handle_err(status)

  end do !months
 end do ! vars loop

122 continue !error thrown by file

 deallocate(ctem_d_m) 
 
 end if !lexist 
 
end if !if fire.

if (DOWETLANDS) then 
  
   inquire(file=trim(folder)//trim(ARGBUFF)//'.CT08M_G',exist=lexist)
   if (lexist) then

   OPEN(910,FILE=trim(folder)//trim(ARGBUFF)//'.CT08M_G',status='old',form='formatted') ! MONTHLY CTEM wetlands
   
   allocate(ctem_w_m(nctemwetvars_m,monyrs,12))
   
   !  first throw out header (6 lines)
   do h = 1,6
	read(910,*,iostat=io_set)
	if (io_set .ne. 0) then
      write(*,*)'Missing/truncated file ',trim(ARGBUFF)//'.CT08M_G'
      close(910)
      goto 123
    end if   
   end do

!---Read in Variables
   do y = 1,monyrs
    do m=1,12
        read(910,*,iostat=io_set) dummy_month,dummy_year,ctem_w_m(1:nctemwetvars_m,y,m)
        if (io_set .ne. 0 .and. y < monyrs .and. m < 12) then
         write(*,*)'Missing/truncated file ',trim(ARGBUFF)//'.CT08M_G'        
         close(910)
         goto 123
        end if   
    end do !m loop
   end do !y loop
   
   close(910)

 if (net4) then
   status = nf90_inq_ncid(ncid_m,'Monthly-Methane flux GridAvg', grpid)
   if (status /= nf90_noerr) call handle_err(status)
 end if
 
 do v = 1,nctemwetvars_m ! begin vars loop
  do m=1,12  !Begin Month Loop
   status = nf90_inq_varid(grpid,trim(CTEM_M_W_VAR(v)), var_id)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_var(grpid,var_id,ctem_w_m(v,:,m),start=[xlon,ylat,m,yrst],count=[1,1,1,monyrs])
   if (status/=nf90_noerr) call handle_err(status)

  end do !months
 end do ! vars loop

123 continue !error thrown by file

 deallocate(ctem_w_m) 
 
 end if !lexist
   
end if !if wetlands

!===================per PFT/Tile===CTEM Monthly==========================

!---Read in Variables

! Monthly CTEM
if (MOSAIC) then
  infile=trim(folder)//trim(ARGBUFF)//'.CT01M_M'
else
  infile=trim(folder)//trim(ARGBUFF)//'.CT01M_G'
end if

inquire(file=infile,exist=lexist)
if (lexist) then

! Make a file of the PFT info, removing the grdav.
command='sed '//tic//'/GRDAV/d'//tic//' '//trim(infile)//' > tmp_m_p.dat'
call system(command)
  
! Open the newly trimmed file
OPEN(741,FILE='tmp_m_p.dat',status='old',form='formatted') ! MONTHLY OUTPUT FOR CTEM

!Allocate Arrays
allocate(ctem_m_mos(numctemvars_m,ntile,monyrs,12))
allocate(tmp(numctemvars_m))

!  first throw out header (6 lines) as it is there still
do h = 1,6
 read(741,*,iostat=io_set)
	if (io_set .ne. 0) then
      write(*,*)'Missing/truncated file ',trim(ARGBUFF)//'.CT01M_G/M'
      close(741)
      deallocate(tmp)
      goto 124
    end if   
end do

! We have to keep track of the month that is read in as it is the only way we know that we are done the tiles for a gridcell.

   do y = 1,monyrs
    do m=1,12

      mo=m  !initialize the month counter to the month of the outer loop

      do while (mo == m) 

          read(741,*,iostat=io_set,END=11) mo,yrin,tmp(1:numctemvars_m),dummy,dummynum,dummy,tilnum
          if (io_set .ne. 0 .and. y < monyrs .and. m < 12) then
            write(*,*)'Missing/truncated file ',trim(ARGBUFF)//'.CT01M_G/M'          
            close(741)
            deallocate(tmp)
            goto 124
          end if   

          ! Figure out what tiles to write
          if (y == 1 .and. mo == 1) then
             tiles_to_write = cshift(tiles_to_write,1)
             tiles_to_write(1) = tilnum
          end if     

        if (mo == m) then
        ! Assign that just read in to its vars
          do v = 1,numctemvars_m ! begin vars loop
            ctem_m_mos(v,tilnum,y,m)=tmp(v)
          end do
        else
            backspace(741)  ! go back one line in the file
        end if

     end do !mo=m loop
    end do !months
   end do ! years

11 continue

 close(741)
 deallocate(tmp)
 
!----
! MONTHLY
if (net4) then
  status = nf90_inq_ncid(ncid_m,'CTEM-Monthly Tiled', grpid)
  if (status /= nf90_noerr) call handle_err(status)
end if

 do v = 1,numctemvars_m ! begin vars loop
  do m=1,12  !Begin Month Loop
   do l = 1,ntile ! begin tile loop

    if (tiles_to_write(l) .ne. -1) then 
      status = nf90_inq_varid(grpid,trim(CTEM_M_VAR(v)), var_id)
      if (status/=nf90_noerr) call handle_err(status)
 
      status = nf90_put_var(grpid,var_id,ctem_m_mos(v,tiles_to_write(l),:,m),start=[xlon,ylat,tiles_to_write(l),m,yrst],count=[1,1,1,1,monyrs])
      if (status/=nf90_noerr) call handle_err(status)
    end if
    
   end do
  end do
 end do

124 continue !error thrown by file

! deallocate arrays

deallocate(ctem_m_mos)

end if !lexist

! ----Fire ----------------------------------
if (DOFIRE) then
  
      if (MOSAIC) then
        infile=trim(folder)//trim(ARGBUFF)//'.CT06M_M'
      else
        infile=trim(folder)//trim(ARGBUFF)//'.CT06M_G'
      end if

    inquire(file=infile,exist=lexist)
    if (lexist) then

  command='sed '//tic//'/GRDAV/d'//tic//' '//trim(infile)//' > tmp_m_p_d.dat'
  call system(command)
  OPEN(841,FILE='tmp_m_p_d.dat',status='old',form='formatted') 

  allocate(ctem_d_m_mos(nctemdistvars_m,ntile,monyrs,12))
  allocate(tmpd(nctemdistvars_m))

  !  first throw out header (6 lines) as it is there still
  do h = 1,6
   read(841,*,iostat=io_set)
	if (io_set .ne. 0) then
      write(*,*)'Missing/truncated file ',trim(ARGBUFF)//'.CT06M_G/M'
      close(841)
      deallocate(tmpd)
      goto 125
    end if 
  end do

  do y = 1,monyrs
   do m=1,12
      mo=m  !initialize the month counter to the month of the outer loop
      do while (mo == m) 
          read(841,*,iostat=io_set,END=12) mo,yrin,tmpd(1:nctemdistvars_m),dummy,dummynum,dummy,tilnum
          if (io_set .ne. 0 .and. y < monyrs .and. m < 12) then
             write(*,*)'Missing/truncated file ',trim(ARGBUFF)//'.CT06M_G/M'          
             close(841)
             deallocate(tmpd)
             goto 125
          end if   
          if (mo == m) then
            do v = 1,nctemdistvars_m ! begin vars loop
              ctem_d_m_mos(v,tilnum,y,m)=tmpd(v)
            end do
          else
           backspace(841)  ! go back one line in the file
          end if
       end do !mo=m loop
    end do !months
   end do ! years
12 continue
   close(841)
   deallocate(tmpd)

  if (net4) then
    status = nf90_inq_ncid(ncid_m,'Monthly-Disturbance Tiled', grpid)
    if (status /= nf90_noerr) call handle_err(status)
  end if

 do v = 1,nctemdistvars_m ! begin vars loop
  do m=1,12  !Begin Month Loop
   do l = 1,ntile ! begin tile loop

    if (tiles_to_write(l) .ne. -1) then

      status = nf90_inq_varid(grpid,trim(CTEM_M_D_VAR(v)), var_id)
      if (status/=nf90_noerr) call handle_err(status)
 
      status = nf90_put_var(grpid,var_id,ctem_d_m_mos(v,tiles_to_write(l),:,m),start=[xlon,ylat,tiles_to_write(l),m,yrst],count=[1,1,1,1,monyrs])
      if (status/=nf90_noerr) call handle_err(status)
     end if

   end do
  end do
 end do
 
 125 continue !error thrown by file
 
 deallocate(ctem_d_m_mos)
 
 end if ! lexist
 
end if !dofire

end if !if ctem. 

end do !loop over the gridcells

!close the monthly netcdf
status = nf90_close(ncid_m)
if (status/=nf90_noerr) call handle_err(status)

! Close the gridcells file
 close(12)

end if ! makemonthly

deallocate(tiles_to_write)

! remove the tmp files
command='rm tmp*.dat'
call system(command)

end program netcdf_writer_fast

!=====================

subroutine handle_err(status)

  use netcdf

  implicit none

  !Internal subroutine - checks error status after each netcdf call,
  !prints out text message each time an error code is returned.

  integer, intent(in) :: status

  if(status /= nf90_noerr) then
    write(0,*)'netCDF error: ',trim(nf90_strerror(status))
    !stop
  end if

end subroutine handle_err

!====================================

subroutine handle_major_err(status)

  use netcdf

  implicit none

  !Internal subroutine - checks error status after each netcdf call,
  !prints out text message each time an error code is returned.

  integer, intent(in) :: status

  if(status /= nf90_noerr) then
    write(0,*)'netCDF error: ',trim(nf90_strerror(status))
    stop
  end if

end subroutine handle_major_err

!------------------------------
subroutine parsecoords(coordstring,val)

!subroutine to parse a coordinate string

implicit none

character(45),      intent(in)  :: coordstring
integer, dimension(4), intent(out) :: val

character(10), dimension(4) :: cval

integer :: i
integer :: lasti = 1
integer :: part  = 1


cval(:) = '0'
lasti = 1
part = 1

do i=1,len_trim(coordstring)
  if (coordstring(i:i) == '_') then
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


