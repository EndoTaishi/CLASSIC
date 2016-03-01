program netcdf_writer_fast

!Program to write CLASS and CTEM output variables to NETCDF files

use netcdf
use creator_module_bf

implicit none

! Annual COMPOSITE 
real, allocatable, dimension(:,:) :: ctem_a
real, allocatable, dimension(:,:) :: ctem_d_a
real, allocatable, dimension(:,:) :: ctem_w_a

! Annual tile avg
real, allocatable, dimension(:,:,:) :: ctem_a_mos
real, allocatable, dimension(:,:,:) :: ctem_d_a_mos
real, allocatable, dimension(:,:,:) :: ctem_w_a_mos

! Annual pft level
real, allocatable, dimension(:,:,:,:) :: ctem_a_pft
real, allocatable, dimension(:,:,:,:) :: ctem_d_a_pft

! Annual Competition/luc
character(1), allocatable, dimension(:,:) :: pftexist
real, allocatable, dimension(:,:) :: pftf
real, allocatable, dimension(:) :: pft_tot

! Monthly COMPOSITE
real, allocatable, dimension(:,:) :: ctem_m
real, allocatable, dimension(:,:) :: ctem_d_m
real, allocatable, dimension(:,:) :: ctem_w_m

! Monthly tile avg
real, allocatable, dimension(:,:,:) :: ctem_m_mos
real, allocatable, dimension(:,:,:) :: ctem_d_m_mos
real, allocatable, dimension(:,:,:) :: ctem_w_m_mos

! Monthly pft level
real, allocatable, dimension(:,:,:,:) :: ctem_m_pft
real, allocatable, dimension(:,:,:,:) :: ctem_d_m_pft

! Monthly Competition/luc
character(1), allocatable, dimension(:,:) :: mpftexist
real, allocatable, dimension(:,:) :: mpftf
real, allocatable, dimension(:) :: mpft_tot

! Monthly variables for CLASS
real, allocatable, dimension(:,:) :: class_m
real, allocatable, dimension(:,:,:) :: class_s_m

! Annual variables for CLASS
real, allocatable, dimension(:,:) :: class_a

integer :: y,dummy_year,dummy_month
real :: dummy_var
integer, dimension(1) :: lyear
integer :: i,m,xlon,ylat,yrst,h,l,cell
integer :: totyrs
integer :: monyrs
integer :: totmons
integer :: numtiles
integer :: yrin
integer :: yr_now
integer :: mo
integer :: tilnum
integer :: pftnum
integer :: dummynum
integer :: var_id,v,p,grpid,t
character(120) :: ARGBUFF 
character(120) :: file_to_write
character(140) :: file_to_write_extended
character(120) :: folder
character(120) :: long_path
character(120) :: jobfile
character(512) :: command
character(160) :: infile
character(160) :: cellsfile
character(4) :: dummy,dummy2
character(6) :: lon_in, lat_in
logical :: CTEM
logical :: MAKEMONTHLY
logical :: DOFIRE
logical :: DOWETLANDS
logical :: TILED
logical :: DOPFTS
logical :: PARALLELRUN
logical :: COMPETE_LNDUSE
integer :: nl ! number of soil layers
integer :: realyrst
integer :: tailyrs
character(1) :: tic
character(len=8) :: fmt ! format descriptor
character(len=10) :: x1
integer :: io_set
logical :: lexist

logical, allocatable, dimension(:,:) :: pfts_to_write

real, allocatable, dimension(:) :: tmp
real, allocatable, dimension(:) :: tmpd

real, allocatable, dimension(:,:) :: tmpa
real, allocatable, dimension(:,:) :: tmpm

namelist /joboptions/ &
  PARALLELRUN,        &
  CTEM,               &
  MAKEMONTHLY,        &
  DOFIRE,             &
  DOWETLANDS,         &   
  TILED,              &
  DOPFTS,             &
  COMPETE_LNDUSE,     &
  nl,                 &
  totyrs,             &
  monyrs,             &
  yrst,               &
  realyrst ,          &
  long_path ,         &
  file_to_write                
!----

tic=char(39)

allocate(pfts_to_write(ntile,ctemnpft))

! Open the grid cells file
 call getarg(1,cellsfile)
 
 call getarg(2,jobfile)

 open(10,file=jobfile,status='old')

 read(10,nml = joboptions)

 close(10)

!********************************************************************

totmons = monyrs * 12

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


    read(11,*)lon_in,lat_in

    ARGBUFF = trim(lon_in)//'_'//trim(lat_in)

    write(*,'(a9,i5,a5,i5,a8,a10)')'annual: ',cell,' of ',num_land_cells,' cell: ',ARGBUFF

    call parsecoords(ARGBUFF,bounds)

    folder=trim(long_path)//'/'//trim(ARGBUFF)//'/'

    ! get the coordinates
    xlon=bounds(1)
    ylat=bounds(3)

    pfts_to_write = .False. ! initialize to -1 now.
    numtiles = 0

!==CLASS ANNUAL =====================================================

! OF1Y

if (.NOT. net4) then
   grpid=ncid
end if

inquire(file=trim(folder)//trim(ARGBUFF)//'.OF1Y',exist=lexist)
if (lexist) then

    OPEN(83,FILE=trim(folder)//trim(ARGBUFF)//'.OF1Y',status='old',form='formatted') ! YEARLY OUTPUT FOR CLASS

    !  first throw out header
       do h = 1,5
    	read(83,*,iostat=io_set)
    	if (io_set .ne. 0) then
    	    write(*,*)'Missing/truncated file',trim(ARGBUFF)//'.OF1Y'
    	    close(83)
    	    goto 111
    	end if     
       end do

    !Allocate arrays
    allocate(class_a(numclasvars_a,totyrs))
    
       do y = 1,totyrs
    
         read(83,*,iostat=io_set)dummy_year,class_a(1:numclasvars_a,y)
         if (io_set .ne. 0 .and. y < totyrs) then
    	    write(*,*)'Missing/truncated file',trim(ARGBUFF)//'.OF1Y'
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

inquire(file=trim(folder)//trim(ARGBUFF)//'.CT07Y',exist=lexist)

if (lexist) then
    OPEN(92,FILE=trim(folder)//trim(ARGBUFF)//'.CT07Y',status='old',form='formatted') ! YEARLY PFT FRAC OUTPUT FOR CTEM

    allocate(pftexist(ctemnpft,totyrs))
    allocate(tmpa(ctemnpft,totyrs))
    allocate(pftf(ctemnpft,totyrs))
    allocate(pft_tot(totyrs))

   do h = 1,6
	read(92,*,iostat=io_set)
	if (io_set .ne. 0) then
      write(*,*)'Missing/truncated file',trim(ARGBUFF)//'.CT07Y'
      close(92)
      goto 112
    end if     
   end do

   do y = 1,totyrs

    read(92,*,iostat=io_set)dummy_year,pftf(1:ctemnpft,y),pft_tot(y),dummy_var,pftexist(1:ctemnpft,y)
    if (io_set .ne. 0 .and. y < totyrs) then
       write(*,*)'Missing/truncated file',trim(ARGBUFF)//'.CT07Y'
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

  ! Annual CTEM
  infile=trim(folder)//trim(ARGBUFF)//'.CT01Y'

inquire(file=infile,exist=lexist)
if (lexist) then

  command='sed -n '//tic//'/GRDAV/p'//tic//' '//trim(infile)//' > tmp_a.dat'  !this removes any header!

  call system(command)

  OPEN(75,FILE='tmp_a.dat',status='old',form='formatted') ! YEARLY OUTPUT FOR CTEM
  
  !allocate the size of the arrays for the output data
  allocate(ctem_a(numctemvars_a,totyrs))

  do y = 1,totyrs
     read(75,*,iostat=io_set)dummy_year,ctem_a(1:numctemvars_a,y)
     if (io_set .ne. 0 .and. y < totyrs) then
        write(*,*)'Missing/truncated file',trim(ARGBUFF)//'.CT01Y'
        close(75)
        deallocate(ctem_a)
        goto 115
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

  ! deallocate arrays
  deallocate(ctem_a)

!now per Tile MODE  for CTEM Annual ====================================================

  if (TILED) then
! Read in from the ascii file

! Allocate the size of the arrays for the output data
  allocate(ctem_a_mos(numctemvars_a,ntile,totyrs))
  allocate(tmp(numctemvars_a))

   ! Make a file of the TILE AVG info
   command='sed -n '//tic//'/TFRAC/p'//tic//' '//trim(infile)//' > tmp_a_t.dat'  !this removes any header!
   call system(command)

   ! Check if there was only one tile and if so then use the grid avg value and sub in the format of the tile avg so it works well
   command='[ -s tmp_a_t.dat ] || sed -n '//tic//'/GRDAV/p'//tic//' '//trim(infile)//' > tmp_a_t.dat ; sed -i '//tic//'s/GRDAV/TILE  1   OF  1 TFRAC/g'//tic//' tmp_a_t.dat'
   call system(command)

   OPEN(750,FILE='tmp_a_t.dat',status='old',form='formatted') ! YEARLY OUTPUT FOR CTEM

   read(750,*,iostat=io_set,end=9000) yrin,tmp(1:numctemvars_a),dummy,tilnum,dummy,dummynum
   backspace(750)

! We have to keep track of the year that is read in as it is the only way we know that we are done the tiles for a gridcell.
    yrin=realyrst

    do y = 1,totyrs

      yr_now = realyrst + y - 1

      do while (yrin == yr_now)

          read(750,*,iostat=io_set,end=90) yrin,tmp(1:numctemvars_a),dummy,tilnum,dummy,dummynum
          if (io_set .ne. 0 .and. y < totyrs) then
            write(*,*)'Missing/truncated file ',trim(ARGBUFF)//'.CT01Y tiled'
            close(750)
            deallocate(tmp)
            deallocate(ctem_a_mos)
            goto 115
          end if

          if (yrin == realyrst) then
            numtiles=max(tilnum,numtiles)
          end if

        if (yrin == yr_now) then
          ! Assign that just read in to its vars
          do v = 1,numctemvars_a ! begin vars loop
            ctem_a_mos(v,tilnum,y)=tmp(v)
          end do
        else
            backspace(750)
        end if
       end do ! while loop
     end do ! y loop

90     continue

! done with ascii file, close it.
  close(750)
  deallocate(tmp)

if (net4) then
  status = nf90_inq_ncid(ncid,'CTEM-Annual Tiled', grpid)
  if (status /= nf90_noerr) call handle_err(status)
end if

 do l=1,numtiles
  do v = 1,numctemvars_a ! begin vars loop

     status = nf90_inq_varid(grpid,trim(CTEM_Y_VAR_TA(v)), var_id)  !tiled vars
     if (status/=nf90_noerr) call handle_err(status)

     status = nf90_put_var(grpid,var_id,ctem_a_mos(v,l,:),start=[xlon,ylat,l,yrst],count=[1,1,1,totyrs])
     if (status/=nf90_noerr) call handle_err(status)

  end do ! vars loop
 end do !tiles loop

9000 continue ! This is a composite run so just move on.

! deallocate arrays
deallocate(ctem_a_mos)

  end if !TILED

!now per PFT MODE for CTEM Annual ====================================================

  if (DOPFTS) then

  numtiles = 0 ! in case you didn't run with TILED, we can redetermine this here.

! Read in from the ascii file

! Allocate the size of the arrays for the output data
  allocate(ctem_a_pft(numctemvars_a,ctemnpft,ntile,totyrs))

   ! Make a file of the PFT info, removing the grdav and TFRAC and the header (first 6 lines)
   command='sed '//tic//'/GRDAV/d'//tic//' '//trim(infile)//' | sed '//tic//'/TFRAC/d'//tic//' | sed '//tic//'1,6d'//tic//' > tmp_a_p.dat'
   call system(command)
   OPEN(751,FILE='tmp_a_p.dat',status='old',form='formatted') ! YEARLY OUTPUT FOR CTEM

  allocate(tmp(numctemvars_a))

! We have to keep track of the year that is read in as it is the only way we know that we are done the PFTs for a gridcell.
    yrin=realyrst

    do y = 1,totyrs

      yr_now = realyrst + y - 1

      do while (yrin == yr_now)

          read(751,*,iostat=io_set,end=95) yrin,tmp(1:numctemvars_a),dummy,tilnum,dummy2,pftnum

          if (io_set .ne. 0 .and. y < totyrs) then
            write(*,*)'Missing/truncated file ',trim(ARGBUFF)//'.CT01Y pftlevel'
            close(751)
            deallocate(ctem_a_pft)
            deallocate(tmp)
            goto 115
          end if

          if (yrin == realyrst) then

             pfts_to_write(tilnum,pftnum) = .True.
             numtiles=max(tilnum,numtiles)


          end if

        if (yrin == yr_now) then

          ! Assign that just read in to its vars
          do v = 1,numctemvars_a ! begin vars loop
            ctem_a_pft(v,pftnum,tilnum,y)=tmp(v)
          end do
        else
            backspace(751)
        end if
       end do ! while yr loop
     end do ! y loop

95     continue

! done with ascii file, close it.
  close(751)
  deallocate(tmp)

if (net4) then
  status = nf90_inq_ncid(ncid,'CTEM-Annual Tiled', grpid)
  if (status /= nf90_noerr) call handle_err(status)
end if

 do l=1,numtiles
  do p = 1, ctemnpft
   if (pfts_to_write(l,p)) then
   do v = 1,numctemvars_a ! begin vars loop
     if (TILED .and. DOPFTS) then

      status = nf90_inq_varid(grpid,trim(CTEM_Y_VAR(v)), var_id) !pft level vars
      if (status/=nf90_noerr) call handle_err(status)

      status = nf90_put_var(grpid,var_id,ctem_a_pft(v,p,l,:),start=[xlon,ylat,p,l,yrst],count=[1,1,1,1,totyrs])
      if (status/=nf90_noerr) call handle_err(status)

    else if (DOPFTS .and. .not. TILED) then

      status = nf90_inq_varid(grpid,trim(CTEM_Y_VAR(v)), var_id) !pft level vars
      if (status/=nf90_noerr) call handle_err(status)

      status = nf90_put_var(grpid,var_id,ctem_a_pft(v,p,1,:),start=[xlon,ylat,p,yrst],count=[1,1,1,totyrs])
      if (status/=nf90_noerr) call handle_err(status)

    end if
   end do ! vars loop
   end if
  end do ! pft loop
 end do !tiles loop

! deallocate arrays
deallocate(ctem_a_pft)

115 continue !error thrown by file

end if !lexist

end if ! DOPFTS

!=====Fire Annual Grid Average=============================================

if (DOFIRE) then

     infile=trim(folder)//trim(ARGBUFF)//'.CT06Y'

    inquire(file=infile,exist=lexist)
    if (lexist) then

     command='sed -n '//tic//'/GRDAV/p'//tic//' '//trim(infile)//' > tmp_a_d.dat'
     call system(command)
     
     OPEN(85,FILE='tmp_a_d.dat',status='old',form='formatted') 

  allocate(ctem_d_a(nctemdistvars_a,totyrs))
  
  do y = 1,totyrs
     read(85,*,iostat=io_set)dummy_year,ctem_d_a(1:nctemdistvars_a,y)
     if (io_set .ne. 0 .and. y < totyrs) then
         write(*,*)'Missing/truncated file',trim(ARGBUFF)//'.CT06Y'
         close(85)
         deallocate(ctem_d_a)
         goto 215
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

 deallocate(ctem_d_a)

 if (TILED) then

! Per TILE AVG CTEM Annual Fire ==========================================

     allocate(ctem_d_a_mos(nctemdistvars_a,ntile,totyrs))
     allocate(tmpd(nctemdistvars_a))

     ! Grab only the tile avg values
     command='sed -n '//tic//'/TFRAC/p'//tic//' '//trim(infile)//' > tmp_a_d_t.dat'
     call system(command)

     ! Check if there was only one tile and if so then use the grid avg value and sub in the format of the tile avg so it works well
     command='[ -s tmp_a_d_t.dat ] || sed -n '//tic//'/GRDAV/p'//tic//' '//trim(infile)//' > tmp_a_d_t.dat ; sed -i '//tic//'s/GRDAV/TILE  1   OF  1 TFRAC/g'//tic//' tmp_a_d_t.dat'
     call system(command)

     OPEN(851,FILE='tmp_a_d_t.dat',status='old',form='formatted')
     ! test the file
     read(851,*,iostat=io_set,end=9150) yrin,tmpd(1:nctemdistvars_a),dummy,tilnum,dummy,dummynum
     backspace(851)

   yrin=realyrst
   do y = 1,totyrs
     yr_now = realyrst + y - 1
     do while (yrin == yr_now)

       read(851,*,iostat=io_set,end=915) yrin,tmpd(1:nctemdistvars_a),dummy,tilnum,dummy,dummynum

       if (io_set .ne. 0 .and. y < totyrs) then
          write(*,*)'Missing/truncated file ',trim(ARGBUFF)//'.CT06Y tiled',y,totyrs
          close(851)
          deallocate(ctem_d_a_mos)
          deallocate(tmpd)
          goto 215
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
915     continue

  close(851)
  deallocate(tmpd)

  if (net4) then
    status = nf90_inq_ncid(ncid,'Annual-Disturbance Tiled', grpid)
    if (status /= nf90_noerr) call handle_err(status)
  end if

  do l=1,numtiles
   do v = 1,nctemdistvars_a ! begin vars loop

      status = nf90_inq_varid(grpid,trim(CTEM_Y_D_VAR_TA(v)), var_id) !tile avg vars
      if (status/=nf90_noerr) call handle_err(status)

      status = nf90_put_var(grpid,var_id,ctem_d_a_mos(v,l,:),start=[xlon,ylat,l,yrst],count=[1,1,1,totyrs])
      if (status/=nf90_noerr) call handle_err(status)

   end do ! vars loop
  end do !tiles loop

9150 continue ! This was a composite run so just move on.

  deallocate(ctem_d_a_mos)

  end if !TILED

  if (DOPFTS) then

! Per PFT CTEM Annual Fire ==========================================

     allocate(ctem_d_a_pft(nctemdistvars_a,ctemnpft,ntile,totyrs))
     allocate(tmpd(nctemdistvars_a))

     ! Remove the header, grdavg and tfrac values
     command='sed '//tic//'/GRDAV/d'//tic//' '//trim(infile)//' | sed '//tic//'/TFRAC/d'//tic//' | sed '//tic//'1,6d'//tic//' > tmp_a_d_p.dat'
     call system(command)
     OPEN(850,FILE='tmp_a_d_p.dat',status='old',form='formatted')

   yrin=realyrst
   do y = 1,totyrs
     yr_now = realyrst + y - 1
     do while (yrin == yr_now)

       read(850,*,iostat=io_set,end=925) yrin,tmpd(1:nctemdistvars_a),dummy,tilnum,dummy,pftnum
       if (io_set .ne. 0 .and. y < totyrs) then
          write(*,*)'Missing/truncated file ',trim(ARGBUFF)//'.CT06Y pft level'
          close(850)
          deallocate(ctem_d_a_pft)
          deallocate(tmpd)
          goto 215
       end if

       if (yrin == yr_now) then
         do v = 1,nctemdistvars_a ! begin vars loop
           ctem_d_a_pft(v,pftnum,tilnum,y)=tmpd(v)
         end do
       else
         backspace(850)  ! go back one line in the file
       end if
     end do ! while loop
   end do ! y loop
925     continue

  close(850)
  deallocate(tmpd)

  if (net4) then
    status = nf90_inq_ncid(ncid,'Annual-Disturbance PerPFT', grpid)
    if (status /= nf90_noerr) call handle_err(status)
  end if

 do l=1,numtiles
  do p = 1, ctemnpft
    if (pfts_to_write(l,p)) then ! these were determined by the CTO1Y file.
     do v = 1,nctemdistvars_a  ! begin vars loop
      if (TILED .and. DOPFTS) then

        status = nf90_inq_varid(grpid,trim(CTEM_Y_D_VAR(v)), var_id) !pft level vars
        if (status/=nf90_noerr) call handle_err(status)

        status = nf90_put_var(grpid,var_id,ctem_d_a_pft(v,p,l,:),start=[xlon,ylat,p,l,yrst],count=[1,1,1,1,totyrs])
        if (status/=nf90_noerr) call handle_err(status)

      else if (DOPFTS .and. .not. TILED) then

        status = nf90_inq_varid(grpid,trim(CTEM_Y_D_VAR(v)), var_id) !pft level vars
        if (status/=nf90_noerr) call handle_err(status)

        status = nf90_put_var(grpid,var_id,ctem_d_a_pft(v,p,1,:),start=[xlon,ylat,p,yrst],count=[1,1,1,totyrs])
        if (status/=nf90_noerr) call handle_err(status)

      end if

     end do ! vars loop
    end if
  end do ! pft loop
 end do !tiles loop

  deallocate(ctem_d_a_pft)

215 continue !error thrown by file

 end if !lexist

 end if !DOPFTS
 
endif ! if fire.

! =====Annual Wetlands==Grid avg=====================================================================

if (DOWETLANDS) then 

  inquire(file=trim(folder)//trim(ARGBUFF)//'.CT08Y',exist=lexist)
  if (lexist) then

     command='sed -n '//tic//'/GRDAV/p'//tic//' '//trim(infile)//' > tmp_a_w.dat'
     call system(command)

     OPEN(900,FILE='tmp_a_w.dat',status='old',form='formatted') ! YEARLY CTEM wetlands CH4 emission

  allocate(ctem_w_a(nctemwetvars_a,totyrs))

  do y = 1,totyrs   
      read(900,*,iostat=io_set) dummy_year,ctem_w_a(1:nctemwetvars_a,y)
      if (io_set .ne. 0 .and. y < totyrs) then
         write(*,*)'Missing/truncated file ',trim(ARGBUFF)//'.CT08Y'
         close(900)
         deallocate(ctem_w_a)
         goto 315
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

  deallocate(ctem_w_a)

  if (TILED) then
! =====Annual Wetlands==Tile avg=====================================================================

    allocate(ctem_w_a_mos(nctemwetvars_a,ntile,totyrs))
    allocate(tmpd(nctemdistvars_a))

     ! Grab only the tile avg values and remove the header
     command='sed -n '//tic//'/TILE/p'//tic//' '//trim(infile)//' > tmp_a_w_t.dat'
     call system(command)
     ! Check if there was only one tile and if so then use the grid avg value and sub in the format of the tile avg so it works well
     command='[ -s tmp_a_w_t.dat ] || sed -n '//tic//'/GRDAV/p'//tic//' '//trim(infile)//' > tmp_a_w_t.dat ; sed -i '//tic//'s/GRDAV/TILE  1/g'//tic//' tmp_a_w_t.dat'
     call system(command)

     OPEN(901,FILE='tmp_a_w_t.dat',status='old',form='formatted')
     ! test the file
     read(901,*,iostat=io_set,end=10150) yrin,tmpd(1:nctemwetvars_a),dummy,tilnum
     backspace(901)

   yrin=realyrst
   do y = 1,totyrs
     yr_now = realyrst + y - 1
     do while (yrin == yr_now)

       read(901,*,iostat=io_set,end=1015) yrin,tmpd(1:nctemwetvars_a),dummy,tilnum
       if (io_set .ne. 0 .and. y < totyrs) then
          write(*,*)'Missing/truncated file ',trim(ARGBUFF)//'.CT08Y tiled'
          close(901)
          deallocate(ctem_w_a_mos)
          deallocate(tmpd)
          goto 315
       end if

       if (yrin == yr_now) then
         do v = 1,nctemwetvars_a ! begin vars loop
           ctem_w_a_mos(v,tilnum,y)=tmpd(v)
         end do
       else
         backspace(901)  ! go back one line in the file
       end if
     end do ! while loop
   end do ! y loop
1015     continue

  close(901)
  deallocate(tmpd)

  if (net4) then
    status = nf90_inq_ncid(ncid,'Annual-Wetland Tiled', grpid)
    if (status /= nf90_noerr) call handle_err(status)
  end if

  do l=1,numtiles
   do v = 1,nctemwetvars_a ! begin vars loop

      status = nf90_inq_varid(grpid,trim(CTEM_Y_W_T_VAR(v)), var_id) !tile avg vars
      if (status/=nf90_noerr) call handle_err(status)

      status = nf90_put_var(grpid,var_id,ctem_w_a_mos(v,l,:),start=[xlon,ylat,l,yrst],count=[1,1,1,totyrs])
      if (status/=nf90_noerr) call handle_err(status)

   end do ! vars loop
  end do !tiles loop

10150 continue ! composite run so move on.

  deallocate(ctem_w_a_mos)

315 continue !error thrown by file
  
  end if !lexist for wetland file

  end if !TILED
  
end if  !if wetlands. 

end if !if CTEM.

! remove the tmp files
command='rm tmp*.dat'
call system(command)


end do !loop over the gridcells

!close the annual netcdf
status = nf90_close(ncid)
if (status/=nf90_noerr) call handle_err(status)

! Close the gridcells file
 close(11)

! *******************************************************************
!=====================CLASS_MONTHLY files============================
! *******************************************************************

if (MAKEMONTHLY) then

! Open the netcdf file for writing

! Monthly file:
file_to_write_extended = trim(file_to_write)//'_CLASSCTEM_M.nc'

write(*,*)'monthly: writing to',file_to_write_extended

status = nf90_open(file_to_write_extended,nf90_write,ncid_m)
if (status /= nf90_noerr) call handle_major_err(status)

! open the list of coordinates
 open(12,FILE=cellsfile,status='old')

do cell = 1,num_land_cells


    read(12,*)lon_in,lat_in
    ARGBUFF = trim(lon_in)//'_'//trim(lat_in)

    write(*,'(a9,i5,a5,i5,a8,a10)')'monthly:',cell,' of ',num_land_cells,'   cell: ',ARGBUFF

    call parsecoords(ARGBUFF,bounds)
    folder=trim(long_path)//'/'//trim(ARGBUFF)//'/'

    ! get the coordinates
    xlon=bounds(1)
    ylat=bounds(3)

    pfts_to_write = .False.  !initialize to False now.
    numtiles=0


if (.NOT. net4) then
   grpid=ncid_m
end if

! Do OF1M first
if (net4) then
  status = nf90_inq_ncid(ncid_m,'CLASS-Monthly',grpid) 
  if (status /= nf90_noerr) call handle_err(status)
end if

  inquire(file=trim(folder)//trim(ARGBUFF)//'.OF1M',exist=lexist)
  if (lexist) then

! Read in the ascii file
   OPEN(81,FILE=trim(folder)//trim(ARGBUFF)//'.OF1M',status='old',form='formatted') ! MONTHLY OUTPUT FOR CLASS

    !Allocate Arrays
    allocate(class_m(numclasvars_m,totmons))

 !  first throw out header (5 lines)
   do h = 1,5
	read(81,*,iostat=io_set)
	if (io_set .ne. 0) then
      write(*,*)'Missing/truncated file ',trim(ARGBUFF)//'OF1M'
      close(81)
      deallocate(class_m)
      goto 118
    end if             
   end do
   
   do y = 1,totmons
     read(81,*,iostat=io_set)dummy_month,dummy_year,class_m(1:numclasvars_m,y)
     if (io_set .ne. 0 .and. y < totmons) then
         write(*,*)'Missing/truncated file ',trim(ARGBUFF)//'OF1M'
         close(81)
         deallocate(class_m)
         goto 118
     end if   
   end do

   close(81)

! Write the inputs to the netcdf
  do v = 1,numclasvars_m

   status = nf90_inq_varid(grpid,trim(CLASS_M_VAR(v)), var_id)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_var(grpid,var_id,class_m(v,:),start=[xlon,ylat,yrst],count=[1,1,totmons])
   if (status/=nf90_noerr) call handle_err(status)

  end do !Var loop

118 continue !error thrown by file

end if !lexist

!=======================================================================
! Do OF2M_G

  inquire(file=trim(folder)//trim(ARGBUFF)//'.OF2M',exist=lexist)
  if (lexist) then

   OPEN(82,FILE=trim(folder)//trim(ARGBUFF)//'.OF2M',status='old',form='formatted')

   !Allocate Arrays
   allocate(class_s_m(nclassoilvars_m,nl,totmons))

   !  first throw out header (5 lines)
   do h = 1,5
	read(82,*,iostat=io_set)
	if (io_set .ne. 0) then
      write(*,*)'Missing/truncated file ',trim(ARGBUFF)//'OF2M'
      close(82)
      deallocate(class_s_m)
      goto 119
    end if             	
   end do
        
   do y = 1,totmons
     read(82,*,iostat=io_set)dummy_month,dummy_year,class_s_m(1,1,y),class_s_m(2,1,y),class_s_m(3,1,y),class_s_m(1,2,y),class_s_m(2,2,y),class_s_m(3,2,y),class_s_m(1,3,y),class_s_m(2,3,y),class_s_m(3,3,y)
     if (io_set .ne. 0 .and. y < totmons) then
         write(*,*)'Missing/truncated file ',trim(ARGBUFF)//'OF2M'
         close(82)
         deallocate(class_s_m)
         goto 119
     end if   
   end do

   close(82)

!----
! Write to netcdf file
 do l=1,nl   ! begin soil layer loop
  do v = 1,nclassoilvars_m ! begin vars loop

   status = nf90_inq_varid(grpid,trim(CLASS_M_S_VAR(v)), var_id)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_var(grpid,var_id,class_s_m(v,l,:),start=[xlon,ylat,l,yrst],count=[1,1,1,totmons])
   if (status/=nf90_noerr) call handle_err(status)

  end do ! vars loop
 end do ! soil layer loop

! deallocate arrays
deallocate(class_s_m)

119 continue !error thrown by file

end if !lexist

! %%%%%%Now MONTHLY CTEM outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IF (CTEM) THEN

! Start with LNDUSECOMPETE vars
if (COMPETE_LNDUSE) then

  inquire(file=trim(folder)//trim(ARGBUFF)//'.CT07M',exist=lexist)
  if (lexist) then

   OPEN(91,FILE=trim(folder)//trim(ARGBUFF)//'.CT07M',status='old',form='formatted') ! MONTHLY PFT FRAC OUTPUT FOR CTEM

    allocate(mpftexist(ctemnpft,totmons))
    allocate(tmpm(ctemnpft,totmons))
    allocate(mpftf(ctemnpft,totmons))
    allocate(mpft_tot(totmons))

   !  first throw out header (6 lines)
   do h = 1,6
	read(91,*,iostat=io_set)
	if (io_set .ne. 0) then
      write(*,*)'Missing/truncated file ',trim(ARGBUFF)//'.CT07M'
      close(91)
      goto 120
    end if   
   end do
   
   do y = 1,totmons
      read(91,*,iostat=io_set)dummy_month,dummy_year,mpftf(1:ctemnpft,y),mpft_tot(y),dummy_var,mpftexist(1:ctemnpft,y)
      if (io_set .ne. 0 .and. y < totmons) then
         write(*,*)'Missing/truncated file ',trim(ARGBUFF)//'.CT07M'
         close(91)
         goto 120
      end if   
   end do

   close(91)

if (net4) then
  status = nf90_inq_ncid(ncid_m,'CTEM-Monthly GridAvg', grpid)
  if (status /= nf90_noerr) call handle_err(status)
end if

! First write the pft fractions and pftexist
  do p = 1,ctemnpft ! begin pft loop

   status = nf90_inq_varid(grpid,trim(CTEM_M_C_VAR(2)), var_id)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_var(grpid,var_id,mpftf(p,:),start=[xlon,ylat,p,yrst],count=[1,1,1,totmons])
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_inq_varid(grpid,trim(CTEM_M_C_VAR(3)), var_id)
   if (status/=nf90_noerr) call handle_err(status)

     do y = 1,totmons
      if (mpftexist(p,y) == 'F') then
         tmpm(p,y)=0.0
      else
         tmpm(p,y)=1.0
      end if
     end do
  
   status = nf90_put_var(grpid,var_id,tmpm(p,:),start=[xlon,ylat,p,yrst],count=[1,1,1,totmons])
   if (status/=nf90_noerr) call handle_err(status)

  end do ! pfts loop

 ! Now write the total plant fraction
   status = nf90_inq_varid(grpid,trim(CTEM_M_C_VAR(1)), var_id)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_var(grpid,var_id,mpft_tot(:),start=[xlon,ylat,yrst],count=[1,1,totmons])
   if (status/=nf90_noerr) call handle_err(status)

120 continue !error thrown by file

deallocate(mpftexist)
deallocate(mpftf)
deallocate(tmpm)
deallocate(mpft_tot)

end if !lexist

end if !compete/lnduse

!==============================Start Monthly Grid Average Vals CTEM=============================

   infile=trim(folder)//trim(ARGBUFF)//'.CT01M'

inquire(file=infile,exist=lexist)
if (lexist) then

 command='sed -n '//tic//'/GRDAV/p'//tic//' '//trim(infile)//' > tmp_m.dat'
call system(command)

! Open the newly trimmed file
OPEN(74,FILE='tmp_m.dat',status='old',form='formatted') ! MONTHLY OUTPUT FOR CTEM

!Allocate Arrays
allocate(ctem_m(numctemvars_m,totmons))

!---Read in Variables
   do y = 1,totmons
     read(74,*,iostat=io_set)dummy_month,dummy_year,ctem_m(1:numctemvars_m,y)
     if (io_set .ne. 0 .and. y < totmons) then
         write(*,*)'Missing/truncated file ',trim(ARGBUFF)//'.CT01M'
         close(74)
         deallocate(ctem_m)
         goto 1210
     end if   
   end do !y loop
 close(74)

! MONTHLY
if (net4) then
  status = nf90_inq_ncid(ncid_m,'CTEM-Monthly GridAvg', grpid)
  if (status /= nf90_noerr) call handle_err(status)
end if

 do v = 1,numctemvars_m ! begin vars loop

   status = nf90_inq_varid(grpid,trim(CTEM_M_VAR_GA(v)), var_id)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_var(grpid,var_id,ctem_m(v,:),start=[xlon,ylat,yrst],count=[1,1,totmons])
   if (status/=nf90_noerr) call handle_err(status)

 end do ! vars loop

 if (TILED) then
! =================CTEM Monthly Tile Avg Values ===============

! Get only the tile avg values.
 command='sed -n '//tic//'/TFRAC/p'//tic//' '//trim(infile)//' > tmp_m_t.dat'
call system(command)

! Check if there was only one tile and if so then use the grid avg value and sub in the format of the tile avg so it works well
command='[ -s tmp_m_t.dat ] || sed -n '//tic//'/GRDAV/p'//tic//' '//trim(infile)//' > tmp_m_t.dat ; sed -i '//tic//'s/GRDAV/TILE  1   OF  1 TFRAC/g'//tic//' tmp_m_t.dat'
call system(command)


!Allocate Arrays
allocate(ctem_m_mos(numctemvars_m,ntile,totmons))
allocate(tmp(numctemvars_m))

! Open the newly trimmed file
OPEN(741,FILE='tmp_m_t.dat',status='old',form='formatted') ! MONTHLY OUTPUT FOR CTEM
! Test this
read(741,*,iostat=io_set,END=1123) mo,yrin,tmp(1:numctemvars_m),dummy,tilnum,dummy,pftnum
backspace(741)

! We have to keep track of the month that is read in as it is the only way we know that we are done the tiles for a gridcell.

   do y = 1,monyrs
    do m=1,12

      mo=m  !initialize the month counter to the month of the outer loop

      do while (mo == m)

          read(741,*,iostat=io_set,END=11) mo,yrin,tmp(1:numctemvars_m),dummy,tilnum,dummy,pftnum

          if (io_set .ne. 0 .and. y < monyrs .and. m < 12) then
            write(*,*)'Missing/truncated file ',trim(ARGBUFF)//'.CT01M'
            close(741)
            deallocate(tmp)
            deallocate(ctem_m_mos)
            goto 1210
          end if

          ! Figure out what tiles to write
          if (y == 1 .and. mo == 1) then
             numtiles = max(tilnum, numtiles)

          end if

        if (mo == m) then
        ! Assign that just read in to its vars
          t= int((y - 1)* 12 + m)
          do v = 1,numctemvars_m ! begin vars loop
            ctem_m_mos(v,tilnum,t)=tmp(v)
          end do
        else
            backspace(741)  ! go back one line in the file
        end if

     end do !mo=m loop
    end do !months
   end do ! years
11 continue

 close(741)


!----
! MONTHLY
if (net4) then
  status = nf90_inq_ncid(ncid_m,'CTEM-Monthly Tiled', grpid)
  if (status /= nf90_noerr) call handle_err(status)
end if

 do v = 1,numctemvars_m ! begin vars loop
   do l = 1,numtiles ! begin tile loop

      status = nf90_inq_varid(grpid,trim(CTEM_M_VAR_TA(v)), var_id)
      if (status/=nf90_noerr) call handle_err(status)

      status = nf90_put_var(grpid,var_id,ctem_m_mos(v,l,:),start=[xlon,ylat,l,yrst],count=[1,1,1,totmons])
      if (status/=nf90_noerr) call handle_err(status)

   end do
 end do

1123 continue

deallocate(tmp)
deallocate(ctem_m_mos)

end if ! TILED

if (DOPFTS) then

numtiles=0

! ======= CTEM Monthly per PFT value ==================

! Make a file of the PFT info, removing the grdav, tfrac and header
 command='sed '//tic//'/GRDAV/d'//tic//' '//trim(infile)//' | sed '//tic//'/TFRAC/d'//tic//' | sed '//tic//'1,6d'//tic//' > tmp_m_p.dat'
 call system(command)

! Open the newly trimmed file
OPEN(740,FILE='tmp_m_p.dat',status='old',form='formatted') ! MONTHLY OUTPUT FOR CTEM

!Allocate Arrays
allocate(ctem_m_pft(numctemvars_m,ctemnpft,ntile,totmons))
allocate(tmp(numctemvars_m))

! We have to keep track of the month that is read in as it is the only way we know that we are done the tiles for a gridcell.

   do y = 1,monyrs
    do m=1,12

      mo=m  !initialize the month counter to the month of the outer loop

      do while (mo == m)

          read(740,*,iostat=io_set,END=1100) mo,yrin,tmp(1:numctemvars_m),dummy,tilnum,dummy,pftnum
          if (io_set .ne. 0 .and. y < monyrs .and. m < 12) then
            write(*,*)'Missing/truncated file ',trim(ARGBUFF)//'.CT01M pftlevel'
            close(740)
            deallocate(tmp)
            deallocate(ctem_m_pft)
            goto 1210
          end if

          ! Figure out what pfts to write
          if (y == 1 .and. mo == 1) then

             pfts_to_write(tilnum,pftnum) = .True.
             numtiles = max(tilnum, numtiles)

          end if

        if (mo == m) then
        ! Assign that just read in to its vars
          t= int((y - 1)* 12 + m)
          do v = 1,numctemvars_m ! begin vars loop
            ctem_m_pft(v,pftnum,tilnum,t)=tmp(v)
          end do
        else
            backspace(740)  ! go back one line in the file
        end if

     end do !mo=m loop
    end do !months
   end do ! years

1100 continue

 close(740)
 deallocate(tmp)

!----
! MONTHLY
if (net4) then
  status = nf90_inq_ncid(ncid_m,'CTEM-Monthly PFTlevel', grpid)
  if (status /= nf90_noerr) call handle_err(status)
end if

  do p = 1,ctemnpft
   do l = 1,numtiles ! begin tile loop
    if (pfts_to_write(l,p)) then
     do v = 1,numctemvars_m ! begin vars loop

      if (TILED .and. DOPFTS) then

        status = nf90_inq_varid(grpid,trim(CTEM_M_VAR(v)), var_id)
        if (status/=nf90_noerr) call handle_err(status)

        status = nf90_put_var(grpid,var_id,ctem_m_pft(v,p,l,:),start=[xlon,ylat,p,l,yrst],count=[1,1,1,1,totmons])
        if (status/=nf90_noerr) call handle_err(status)

      else if (DOPFTS .and. .not. TILED) then

        status = nf90_inq_varid(grpid,trim(CTEM_M_VAR(v)), var_id)
        if (status/=nf90_noerr) call handle_err(status)

        status = nf90_put_var(grpid,var_id,ctem_m_pft(v,p,1,:),start=[xlon,ylat,p,yrst],count=[1,1,1,totmons])
        if (status/=nf90_noerr) call handle_err(status)

      end if

    end do
   end if
  end do
 end do

deallocate(ctem_m_pft)

1210 continue !error thrown by file

end if  !lexist

end if !DOPFTS

! Fire Monthly Grid Avg ============================================================

if (DOFIRE) then
  
        infile=trim(folder)//trim(ARGBUFF)//'.CT06M'

inquire(file=infile,exist=lexist)
if (lexist) then

     command='sed -n '//tic//'/GRDAV/p'//tic//' '//trim(infile)//' > tmp_m_d.dat'
     call system(command)
     OPEN(84,FILE='tmp_m_d.dat',status='old',form='formatted') 

  allocate(ctem_d_m(nctemdistvars_m,totmons))

!---Read in Variables
    do y = 1,totmons
      read(84,*,iostat=io_set)dummy_month,dummy_year,ctem_d_m(1:nctemdistvars_m,y)
      if (io_set .ne. 0 .and. y < monyrs .and. m < 12) then
         write(*,*)'Missing/truncated file ',trim(ARGBUFF)//'.CT06M'
         close(84)
         deallocate(ctem_d_m)
         goto 1220
      end if   
    end do !y loop

    close(84)

 if (net4) then
   status = nf90_inq_ncid(ncid_m,'Monthly-Disturbance GridAvg', grpid)
   if (status /= nf90_noerr) call handle_err(status)
 end if

 do v = 1,nctemdistvars_m ! begin vars loop
  
   status = nf90_inq_varid(grpid,trim(CTEM_M_D_VAR_GA(v)), var_id)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_var(grpid,var_id,ctem_d_m(v,:),start=[xlon,ylat,yrst],count=[1,1,totmons])
   if (status/=nf90_noerr) call handle_err(status)

 end do ! vars loop

 deallocate(ctem_d_m)

 ! Fire Monthly Tile Avg ============================================================

if (TILED) then

  command='sed -n '//tic//'/TFRAC/p'//tic//' '//trim(infile)//' > tmp_m_t_d.dat'
  call system(command)

  ! Check if there was only one tile and if so then use the grid avg value and sub in the format of the tile avg so it works well
  command='[ -s tmp_m_t_d.dat ] || sed -n '//tic//'/GRDAV/p'//tic//' '//trim(infile)//' > tmp_m_t_d.dat ; sed -i '//tic//'s/GRDAV/TILE  1   OF  1 TFRAC/g'//tic//' tmp_m_t_d.dat'
  call system(command)


  OPEN(841,FILE='tmp_m_t_d.dat',status='old',form='formatted')

  allocate(ctem_d_m_mos(nctemdistvars_m,ntile,totmons))
  allocate(tmpd(nctemdistvars_m))

  do y = 1,monyrs
   do m=1,12
      mo=m  !initialize the month counter to the month of the outer loop
      do while (mo == m)
          read(841,*,iostat=io_set,END=12) mo,yrin,tmpd(1:nctemdistvars_m),dummy,tilnum,dummy,pftnum
          if (io_set .ne. 0 .and. y < monyrs .and. m < 12) then
             write(*,*)'Missing/truncated file ',trim(ARGBUFF)//'.CT06M'
             close(841)
             deallocate(tmpd)
             deallocate(ctem_d_m_mos)
             goto 1220
          end if
          if (mo == m) then
            t= int((y - 1)* 12 + m)
            do v = 1,nctemdistvars_m ! begin vars loop
              ctem_d_m_mos(v,tilnum,t)=tmpd(v)
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
   do l = 1,numtiles ! begin tile loop

      status = nf90_inq_varid(grpid,trim(CTEM_M_D_VAR(v)), var_id)
      if (status/=nf90_noerr) call handle_err(status)

      status = nf90_put_var(grpid,var_id,ctem_d_m_mos(v,l,:),start=[xlon,ylat,l,yrst],count=[1,1,1,totmons])
      if (status/=nf90_noerr) call handle_err(status)


   end do
 end do

 deallocate(ctem_d_m_mos)

 end if !TILED

 if (DOPFTS) then

 ! Fire Monthly Per PFT ============================================================

  command='sed '//tic//'/GRDAV/d'//tic//' '//trim(infile)//' | sed '//tic//'/TFRAC/d'//tic//' | sed '//tic//'1,6d'//tic//' > tmp_m_p_d.dat'
  call system(command)
  OPEN(840,FILE='tmp_m_p_d.dat',status='old',form='formatted')

  allocate(ctem_d_m_pft(nctemdistvars_m,ctemnpft,ntile,totmons))
  allocate(tmpd(nctemdistvars_m))

  do y = 1,monyrs
   do m=1,12
      mo=m  !initialize the month counter to the month of the outer loop
      do while (mo == m)
          read(840,*,iostat=io_set,END=121) mo,yrin,tmpd(1:nctemdistvars_m),dummy,tilnum,dummy,pftnum
          if (io_set .ne. 0 .and. y < monyrs .and. m < 12) then
             write(*,*)'Missing/truncated file ',trim(ARGBUFF)//'.CT06M pftlevel'
             close(840)
             deallocate(tmpd)
             deallocate(ctem_d_m_pft)
             goto 1220
          end if
          if (mo == m) then
            t= int((y - 1)* 12 + m)
            do v = 1,nctemdistvars_m ! begin vars loop
              ctem_d_m_pft(v,pftnum,tilnum,t)=tmpd(v)
            end do
          else
           backspace(840)  ! go back one line in the file
          end if
       end do !mo=m loop
    end do !months
   end do ! years
121 continue
   close(840)
   deallocate(tmpd)

  if (net4) then
    status = nf90_inq_ncid(ncid_m,'Monthly-Disturbance PFTlevel', grpid)
    if (status /= nf90_noerr) call handle_err(status)
  end if

  do p = 1,ctemnpft
   do l = 1,numtiles ! begin tile loop
    if (pfts_to_write(l,p)) then
     do v = 1,nctemdistvars_m ! begin vars loop

      if (TILED .and. DOPFTS) then

        status = nf90_inq_varid(grpid,trim(CTEM_M_D_VAR(v)), var_id)
        if (status/=nf90_noerr) call handle_err(status)

        status = nf90_put_var(grpid,var_id,ctem_d_m_pft(v,p,l,:),start=[xlon,ylat,p,l,yrst],count=[1,1,1,1,totmons])
        if (status/=nf90_noerr) call handle_err(status)

      else if (DOPFTS .and. .not. TILED) then

        status = nf90_inq_varid(grpid,trim(CTEM_M_D_VAR(v)), var_id)
        if (status/=nf90_noerr) call handle_err(status)

        status = nf90_put_var(grpid,var_id,ctem_d_m_pft(v,p,1,:),start=[xlon,ylat,p,yrst],count=[1,1,1,totmons])
        if (status/=nf90_noerr) call handle_err(status)

      end if

    end do
   end if
  end do
 end do

 deallocate(ctem_d_m_pft)

1220 continue !error thrown by file

 end if !lexist 

 end if ! DOPFTS
 
end if !if fire.

! Monthly grid avg wetlands ====================================================================

if (DOWETLANDS) then 
  
    infile=trim(folder)//trim(ARGBUFF)//'.CT08M'
    inquire(file=infile,exist=lexist)

   if (lexist) then

     command='sed -n '//tic//'/GRDAV/p'//tic//' '//trim(infile)//' > tmp_m_w.dat'
     call system(command)
     OPEN(910,FILE='tmp_m_w.dat',status='old',form='formatted') ! MONTHLY CTEM wetlands
   
   allocate(ctem_w_m(nctemwetvars_m,totmons))
   
!---Read in Variables
   do y = 1,totmons
        read(910,*,iostat=io_set) dummy_month,dummy_year,ctem_w_m(1:nctemwetvars_m,y)
        if (io_set .ne. 0 .and. y < totmons) then
         write(*,*)'Missing/truncated file ',trim(ARGBUFF)//'.CT08M'
         close(910)
         deallocate(ctem_w_m)
         goto 123
        end if   
   end do !y loop
   
   close(910)

 if (net4) then
   status = nf90_inq_ncid(ncid_m,'Monthly-Methane flux GridAvg', grpid)
   if (status /= nf90_noerr) call handle_err(status)
 end if
 
 do v = 1,nctemwetvars_m ! begin vars loop
   status = nf90_inq_varid(grpid,trim(CTEM_M_W_VAR(v)), var_id)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_var(grpid,var_id,ctem_w_m(v,:),start=[xlon,ylat,yrst],count=[1,1,totmons])
   if (status/=nf90_noerr) call handle_err(status)

 end do ! vars loop

deallocate(ctem_w_m)

if (TILED) then
!===================per Tile===CTEM Monthly wetlands==========================

 command='sed -n '//tic//'/TILE/p'//tic//' '//trim(infile)//' > tmp_m_t_w.dat'
  call system(command)

  ! Check if there was only one tile and if so then use the grid avg value and sub in the format of the tile avg so it works well
  command='[ -s tmp_m_t_w.dat ] || sed -n '//tic//'/GRDAV/p'//tic//' '//trim(infile)//' > tmp_m_t_w.dat ; sed -i '//tic//'s/GRDAV/TILE  1/g'//tic//' tmp_m_t_w.dat'
  call system(command)

  OPEN(911,FILE='tmp_m_t_w.dat',status='old',form='formatted')

  allocate(ctem_w_m_mos(nctemwetvars_m,ntile,totmons))
  allocate(tmpd(nctemwetvars_m))

  do y = 1,monyrs
   do m=1,12
      mo=m  !initialize the month counter to the month of the outer loop
      do while (mo == m)
          read(911,*,iostat=io_set,END=129) mo,yrin,tmpd(1:nctemwetvars_m),dummy,tilnum,dummy,pftnum
          if (io_set .ne. 0 .and. y < monyrs .and. m < 12) then
             write(*,*)'Missing/truncated file ',trim(ARGBUFF)//'.CT06M'
             close(911)
             deallocate(tmpd)
             deallocate(ctem_w_m_mos)
             goto 123
          end if
          if (mo == m) then
            t= int((y - 1)* 12 + m)
            do v = 1,nctemwetvars_m ! begin vars loop
              ctem_w_m_mos(v,tilnum,t)=tmpd(v)
            end do
          else
           backspace(911)  ! go back one line in the file
          end if
       end do !mo=m loop
    end do !months
   end do ! years
129 continue
   close(911)
   deallocate(tmpd)

  if (net4) then
    status = nf90_inq_ncid(ncid_m,'Monthly-Wetland Tiled', grpid)
    if (status /= nf90_noerr) call handle_err(status)
  end if

 do v = 1,nctemwetvars_m ! begin vars loop
   do l = 1,numtiles ! begin tile loop

      status = nf90_inq_varid(grpid,trim(CTEM_M_W_T_VAR(v)), var_id)
      if (status/=nf90_noerr) call handle_err(status)

      status = nf90_put_var(grpid,var_id,ctem_w_m_mos(v,l,:),start=[xlon,ylat,l,yrst],count=[1,1,1,totmons])
      if (status/=nf90_noerr) call handle_err(status)

   end do
 end do

 deallocate(ctem_w_m_mos)

123 continue !error thrown by monthly wetland file

 end if !lexist

 end if !TILED
   
end if !if wetlands

end if !if ctem. 

end do !loop over the gridcells

!close the monthly netcdf
status = nf90_close(ncid_m)
if (status/=nf90_noerr) call handle_err(status)

! Close the gridcells file
 close(12)

end if ! makemonthly

deallocate(pfts_to_write)

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


