program netcdf_writer_bf

!Program to write CLASS and CTEM output variables to NETCDF files

use netcdf
use creator_module_bf

implicit none

! Annual COMPOSITE 
real, allocatable, dimension(:,:) :: ctem_a
real, allocatable, dimension(:,:) :: ctem_d_a

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
integer :: i,m,xlon,ylat,yrst,h,l
integer :: totyrs
integer :: yrin
integer :: yr_now
integer :: mo
integer :: tilnum
integer :: var_id,v,p,grpid
character(120) :: ARGBUFF 
character(120) :: file_to_write
character(140) :: file_to_write_extended
character(120) :: folder
character(120) :: long_path
character(120) :: jobfile
character(512) :: command
character(160) :: infile
character(4) :: dummy
logical :: CTEM
logical :: MAKEMONTHLY
logical :: DOFIRE
logical :: MOSAIC
logical :: PARALLELRUN
logical :: COMPETE_LNDUSE
integer :: realyrst
character(1) :: tic

real, allocatable, dimension(:) :: tmp
real, allocatable, dimension(:) :: tmpd

real, allocatable, dimension(:,:) :: tmpa
real, allocatable, dimension(:,:,:) :: tmpm

namelist /joboptions/ &
  PARALLELRUN,        &
  CTEM, 	      &
  MAKEMONTHLY,        &
  DOFIRE,             &
  MOSAIC,             &
  COMPETE_LNDUSE,     &
  totyrs,             &
  yrst,               &
  realyrst ,          &
  long_path ,         &
  file_to_write                
!----

!read the file coordinates

 call getarg(1,ARGBUFF)

 call parsecoords(ARGBUFF,bounds)

! get the coordinates
 xlon=bounds(1)
 ylat=bounds(3)

 call getarg(2,jobfile)

 open(10,file=jobfile,status='old')

 read(10,nml = joboptions)

 close(10)

folder=trim(long_path)//'/'//trim(ARGBUFF)//'/'

! open the input files
        OPEN(81,FILE=trim(folder)//trim(ARGBUFF)//'.OF1M_G',status='old',form='formatted') ! MONTHLY OUTPUT FOR CLASS
        OPEN(82,FILE=trim(folder)//trim(ARGBUFF)//'.OF2M_G',status='old',form='formatted')
        OPEN(83,FILE=trim(folder)//trim(ARGBUFF)//'.OF1Y_G',status='old',form='formatted') ! YEARLY OUTPUT FOR CLASS

       if (CTEM) then

        if (DOFIRE) then
         OPEN(86,FILE=trim(folder)//trim(ARGBUFF)//'.CT06M_G',status='old',form='formatted') ! MONTHLY disturbance OUTPUT FOR CTEM
         OPEN(87,FILE=trim(folder)//trim(ARGBUFF)//'.CT06Y_G',status='old',form='formatted') ! YEARLY disturbance OUTPUT FOR CTEM
        end if

       if (MOSAIC) then
        if (DOFIRE) then
         OPEN(861,FILE=trim(folder)//trim(ARGBUFF)//'.CT06M_M',status='old',form='formatted') ! MONTHLY disturbance OUTPUT FOR CTEM
         OPEN(871,FILE=trim(folder)//trim(ARGBUFF)//'.CT06Y_M',status='old',form='formatted') ! YEARLY disturbance OUTPUT FOR CTEM
        end if
       end if
        
        if (COMPETE_LNDUSE) then
         OPEN(91,FILE=trim(folder)//trim(ARGBUFF)//'.CT07M_GM',status='old',form='formatted') ! MONTHLY PFT FRAC OUTPUT FOR CTEM
         OPEN(92,FILE=trim(folder)//trim(ARGBUFF)//'.CT07Y_GM',status='old',form='formatted') ! YEARLY PFR FRAC OUTPUT FOR CTEM
        end if

       END IF


! Prepare the composite/mosaic CTEM files for read-in. These files have a value per PFT for 
! each variable (plus one for the bare fraction) and then one for the grid
! cell average across all fractions (GRDAV). It seems to be easiest for the 
! netcdf creation to create separate files for the grid-averaged values.
! Do this using sed.

  ! Annual CTEM
  if (MOSAIC) then
  infile=trim(folder)//trim(ARGBUFF)//'.CT01Y_M'
  else
  infile=trim(folder)//trim(ARGBUFF)//'.CT01Y_G'
  end if

  tic=char(39)
  command='sed -n '//tic//'/GRDAV/p'//tic//' '//trim(infile)//' > tmp_a.dat'

  call system(command)

   ! Make a file of the PFT info, removing the grdav.

   command='sed '//tic//'/GRDAV/d'//tic//' '//trim(infile)//' > tmp_a_p.dat'

   call system(command)

   ! Monthly CTEM
  if (MOSAIC) then
  infile=trim(folder)//trim(ARGBUFF)//'.CT01M_M'
  else
  infile=trim(folder)//trim(ARGBUFF)//'.CT01M_G'
  end if

  command='sed -n '//tic//'/GRDAV/p'//tic//' '//trim(infile)//' > tmp_m.dat'

  call system(command)

   ! Make a file of the PFT info, removing the grdav.
   command='sed '//tic//'/GRDAV/d'//tic//' '//trim(infile)//' > tmp_m_p.dat'

   call system(command)

  OPEN(74,FILE='tmp_m.dat',status='old',form='formatted') ! MONTHLY OUTPUT FOR CTEM
  OPEN(75,FILE='tmp_a.dat',status='old',form='formatted') ! YEARLY OUTPUT FOR CTEM

  OPEN(741,FILE='tmp_m_p.dat',status='old',form='formatted') ! MONTHLY OUTPUT FOR CTEM
  OPEN(751,FILE='tmp_a_p.dat',status='old',form='formatted') ! YEARLY OUTPUT FOR CTEM
 
! Open the netcdf file for writing

file_to_write_extended = trim(file_to_write)//'_CLASSCTEM.nc'
status = nf90_open(file_to_write_extended,nf90_write,ncid)
if (status /= nf90_noerr) call handle_err(status)

if (.NOT. net4) then
   grpid=ncid
end if

!=====================CLASS_MONTHLY files============================

if (MAKEMONTHLY) then

! Do OF1M_G first
if (net4) then
  status = nf90_inq_ncid(ncid,'CLASS-Monthly',grpid) 
  if (status /= nf90_noerr) call handle_err(status)
end if

!  first throw out header
   do h = 1,5
	read(81,*)
   end do

!Allocate Arrays
allocate(class_m(numclasvars_m,totyrs,12))

! Read in the ascii file
   do y = 1,totyrs
    do m=1,12
     read(81,*)dummy_month,dummy_year,class_m(1:numclasvars_m,y,m)
    end do
   end do

    close(81)

! Write the inputs to the netcdf
do m=1,12  !Begin Month Loop
  do v = 1,numclasvars_m

   status = nf90_inq_varid(grpid,trim(CLASS_M_VAR(v)), var_id)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_var(grpid,var_id,class_m(v,:,m),start=[xlon,ylat,m,yrst],count=[1,1,1,totyrs])
   if (status/=nf90_noerr) call handle_err(status)

  end do !Var loop
end do   !End month loop

! deallocate arrays
deallocate(class_m)

!=======================================================================
! Do OF2M_G
!  first throw out header
   do h = 1,5
	read(82,*)
   end do

        !Allocate Arrays
        allocate(class_s_m(nclassoilvars_m,nl,totyrs,12))

   do y = 1,totyrs
    do m=1,12
     read(82,*)dummy_month,dummy_year,class_s_m(1,1,y,m),class_s_m(2,1,y,m),class_s_m(3,1,y,m),class_s_m(1,2,y,m),class_s_m(2,2,y,m),class_s_m(3,2,y,m),class_s_m(1,3,y,m),class_s_m(2,3,y,m),class_s_m(3,3,y,m)
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

   status = nf90_put_var(grpid,var_id,class_s_m(v,l,:,m),start=[xlon,ylat,l,m,yrst],count=[1,1,1,1,totyrs])
   if (status/=nf90_noerr) call handle_err(status)

  end do ! vars loop
 end do ! soil layer loop
end do  !End Month Loop

! deallocate arrays
deallocate(class_s_m)

end if ! makemonthly

!===========================CLASS ANNUAL =========================================
! Do OF1Y_G 

!  first throw out header
   do h = 1,5
	read(83,*)
   end do

        !Allocate arrays
        allocate(class_a(numclasvars_a,totyrs))

   do y = 1,totyrs
    
     read(83,*)dummy_year,class_a(1:numclasvars_a,y)

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

! deallocate arrays
deallocate(class_a)

!============================CTEM COMPETITION ANNUAL FILES===============================

IF (CTEM) THEN

! Do the COMPETE_LNDUSE ones first since they are independent of mosaic vs. composite.

if (COMPETE_LNDUSE) then

allocate(pftexist(ctemnpft,totyrs))
allocate(tmpa(ctemnpft,totyrs))
allocate(pftf(ctemnpft,totyrs))
allocate(pft_tot(totyrs))

   do h = 1,6
	read(92,*)
   end do

   do y = 1,totyrs

    read(92,*)dummy_year,pftf(1:ctemnpft,y),pft_tot(y),dummy_var,pftexist(1:ctemnpft,y)

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

deallocate(pftf)
deallocate(pft_tot)
deallocate(pftexist)
deallocate(tmpa)

end if

!============================CTEM Grid Average ANNUAL FILES=========================================\

! Allocate the arrays in preparation for CT01Y_G, CT06Y_G, CT01Y_GM

!allocate the size of the arrays for the output data
allocate(ctem_a(numctemvars_a,totyrs))

if (DOFIRE) then
 allocate(ctem_d_a(nctemdistvars_a,totyrs))
end if

! Read in from the ascii file

!  first throw out header
    if (DOFIRE) then 
     do h = 1,6
	read(87,*)
     end do
    end if

  do y = 1,totyrs

     read(75,*)dummy_year,ctem_a(1:numctemvars_a,y)

    if (DOFIRE) then

     read(87,*)dummy_year,ctem_d_a(1:nctemdistvars_a,y) 

    endif

   end do

    close(75)
    if (DOFIRE) then
        close(87)
    endif

!----
! ANNUAL 

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


 if (DOFIRE) then

   if (net4) then
    status = nf90_inq_ncid(ncid,'Annual-Disturbance GridAvg', grpid)
    if (status /= nf90_noerr) call handle_err(status)
  end if

  do v = 1,nctemdistvars_a ! begin vars loop

    status = nf90_inq_varid(grpid,trim(CTEM_Y_D_VAR(v)), var_id)
    if (status/=nf90_noerr) call handle_err(status)

    status = nf90_put_var(grpid,var_id,ctem_d_a(v,:),start=[xlon,ylat,yrst],count=[1,1,totyrs])
    if (status/=nf90_noerr) call handle_err(status)

  end do ! vars loop

 end if !dofire

! deallocate arrays
deallocate(ctem_a)

if (DOFIRE) then
 deallocate(ctem_d_a)
endif


!now per PFT/Tile MODE ====================================================

! Allocate the arrays in preparation for CT01Y_M and CT06Y_M

! Allocate the size of the arrays for the output data
allocate(ctem_a_mos(numctemvars_a,ntile,totyrs))

if (DOFIRE) then
 allocate(ctem_d_a_mos(nctemdistvars_a,ntile,totyrs))
 allocate(tmpd(nctemdistvars_a))
end if

! Read in from the ascii file

!  first throw out header
   do h = 1,6
        read(751,*)
    if (DOFIRE) then 
	read(871,*)
    endif
   end do

  allocate(tmp(numctemvars_a))

! We have to keep track of the year that is read in as it is the only way we know that we are done the tiles for a gridcell.
    yrin=realyrst

    do y = 1,totyrs

      yr_now = realyrst + y - 1

      do while (yrin == yr_now)

          read(751,*,end=90) yrin,tmp(1:numctemvars_a),dummy,tilnum

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
    
      if (DOFIRE) then
         yrin=realyrst
       do y = 1,totyrs
        yr_now = realyrst + y - 1
        do while (yrin == yr_now)

          read(871,*,end=91) yrin,tmpd(1:nctemdistvars_a),dummy,tilnum
 
          if (yrin == yr_now) then
           do v = 1,nctemdistvars_a ! begin vars loop
             ctem_d_a_mos(v,tilnum,y)=tmpd(v)
           end do
          else
           backspace(871)  ! go back one line in the file
          end if
        end do ! while loop
       end do ! y loop
      end if !dofire
91     continue


! done with ascii file, close it.

      close(751)

    if (DOFIRE) then
        close(871)
        deallocate(tmpd)
    endif

  deallocate(tmp)
!----

! ANNUAL
if (net4) then
  status = nf90_inq_ncid(ncid,'CTEM-Annual Tiled', grpid)
  if (status /= nf90_noerr) call handle_err(status)
end if

 do l=1,ntile
  do v = 1,numctemvars_a ! begin vars loop

   status = nf90_inq_varid(grpid,trim(CTEM_Y_VAR(v)), var_id)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_var(grpid,var_id,ctem_a_mos(v,l,:),start=[xlon,ylat,l,yrst],count=[1,1,1,totyrs])
   if (status/=nf90_noerr) call handle_err(status)

  end do ! vars loop
 end do !tiles loop


if (DOFIRE) then

  if (net4) then
    status = nf90_inq_ncid(ncid,'Annual-Disturbance Tiled', grpid)
    if (status /= nf90_noerr) call handle_err(status)
  end if

 do l=1,ntile
  do v = 1,nctemdistvars_a ! begin vars loop

    status = nf90_inq_varid(grpid,trim(CTEM_Y_D_VAR(v)), var_id)
    if (status/=nf90_noerr) call handle_err(status)

    status = nf90_put_var(grpid,var_id,ctem_d_a_mos(v,l,:),start=[xlon,ylat,l,yrst],count=[1,1,1,totyrs])
    if (status/=nf90_noerr) call handle_err(status)

  end do ! vars loop
 end do !tiles loop
end if !do fire

! deallocate arrays
deallocate(ctem_a_mos)

if (DOFIRE) then
  deallocate(ctem_d_a_mos)
endif


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Now MONTHLY CTEM outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (MAKEMONTHLY) then

! Start with LNDUSECOMPETE vars
if (COMPETE_LNDUSE) then

allocate(mpftexist(ctemnpft,totyrs,12))
allocate(tmpm(ctemnpft,totyrs,12))
allocate(mpftf(ctemnpft,totyrs,12))
allocate(mpft_tot(totyrs,12))

!  first throw out header
   do h = 1,6
        read(91,*)
   end do

   do y = 1,totyrs
    do m=1,12
      read(91,*)dummy_month,dummy_year,mpftf(1:ctemnpft,y,m),mpft_tot(y,m),dummy_var,mpftexist(1:ctemnpft,y,m)
    end do
   end do

   close(91)

if (net4) then
  status = nf90_inq_ncid(ncid,'CTEM-Monthly GridAvg', grpid)
  if (status /= nf90_noerr) call handle_err(status)
end if

! First write the pft fractions and pftexist
 do m=1,12  !Begin Month Loop
  do p = 1,ctemnpft ! begin pft loop

   status = nf90_inq_varid(grpid,trim(CTEM_M_C_VAR(2)), var_id)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_var(grpid,var_id,mpftf(p,:,m),start=[xlon,ylat,p,m,yrst],count=[1,1,1,1,totyrs])
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_inq_varid(grpid,trim(CTEM_M_C_VAR(3)), var_id)
   if (status/=nf90_noerr) call handle_err(status)

     do y = 1,totyrs 
      if (mpftexist(p,y,m) == 'F') then
         tmpm(p,y,m)=0.0
      else
         tmpm(p,y,m)=1.0       
      end if
     end do
  
   status = nf90_put_var(grpid,var_id,tmpm(p,:,m),start=[xlon,ylat,p,m,yrst],count=[1,1,1,1,totyrs])
   if (status/=nf90_noerr) call handle_err(status)


  end do ! pfts loop
 end do !months

 ! Now write the total plant fraction
 do m=1,12  !Begin Month Loop
   status = nf90_inq_varid(grpid,trim(CTEM_M_C_VAR(1)), var_id)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_var(grpid,var_id,mpft_tot(:,m),start=[xlon,ylat,m,yrst],count=[1,1,1,totyrs])
   if (status/=nf90_noerr) call handle_err(status)
 end do

deallocate(mpftexist)
deallocate(mpftf)
deallocate(tmpm)
deallocate(mpft_tot)

end if !compete/lnduse

!==============================Start Monthly Grid Average Vals CTEM=============================

!  first throw out header
    if (DOFIRE) then
     do h = 1,6
	read(86,*)
     end do
    end if

!Allocate Arrays
allocate(ctem_m(numctemvars_m,totyrs,12))

if (DOFIRE) then
  allocate(ctem_d_m(nctemdistvars_m,totyrs,12))
end if

!---Read in Variables
   do y = 1,totyrs
    do m=1,12

     read(74,*)dummy_month,dummy_year,ctem_m(1:numctemvars_m,y,m)

     if (DOFIRE) then

      read(86,*)dummy_month,dummy_year,ctem_d_m(1:nctemdistvars_m,y,m)

     end if

    end do !m loop
   end do !y loop

    close(74)
    if (DOFIRE) then
      close(86)
    end if
!----

! MONTHLY
if (net4) then
  status = nf90_inq_ncid(ncid,'CTEM-Monthly GridAvg', grpid)
  if (status /= nf90_noerr) call handle_err(status)
end if

 do v = 1,numctemvars_m ! begin vars loop
  do m=1,12  !Begin Month Loop
   status = nf90_inq_varid(grpid,trim(CTEM_M_VAR_GA(v)), var_id)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_var(grpid,var_id,ctem_m(v,:,m),start=[xlon,ylat,m,yrst],count=[1,1,1,totyrs])
   if (status/=nf90_noerr) call handle_err(status)

  end do !months
 end do ! vars loop

if (DOFIRE) then

 if (net4) then
   status = nf90_inq_ncid(ncid,'Monthly-Disturbance GridAvg', grpid)
   if (status /= nf90_noerr) call handle_err(status)
 end if

 do v = 1,nctemdistvars_m ! begin vars loop
  do m=1,12  !Begin Month Loop
   status = nf90_inq_varid(grpid,trim(CTEM_M_D_VAR(v)), var_id)
   if (status/=nf90_noerr) call handle_err(status)

   status = nf90_put_var(grpid,var_id,ctem_d_m(v,:,m),start=[xlon,ylat,m,yrst],count=[1,1,1,totyrs])
   if (status/=nf90_noerr) call handle_err(status)

  end do !months
 end do ! vars loop

end if !dofire

! deallocate arrays

deallocate(ctem_m)

if (DOFIRE) then
  deallocate(ctem_d_m)
end if

!===================per PFT/Tile===CTEM Monthly==========================

!  first throw out header
   do h = 1,6
        read(741,*)
    if (DOFIRE) then
	read(861,*)
    end if
   end do

!Allocate Arrays
allocate(ctem_m_mos(numctemvars_m,ntile,totyrs,12))
allocate(tmp(numctemvars_m))

if (DOFIRE) then
  allocate(ctem_d_m_mos(nctemdistvars_m,ntile,totyrs,12))
  allocate(tmpd(nctemdistvars_m))
end if

!---Read in Variables

! We have to keep track of the month that is read in as it is the only way we know that we are done the tiles for a gridcell.

   do y = 1,totyrs
    do m=1,12

      mo=m  !initialize the month counter to the month of the outer loop

      do while (mo == m) 
    
          read(741,*) mo,yrin,tmp(1:numctemvars_m),dummy,tilnum
  
        if (mo == m) then
        ! Assign that just read in to its vars
          do v = 1,numctemvars_a ! begin vars loop
            ctem_m_mos(v,tilnum,y,m)=tmp(v)
          end do
        else
            backspace(741)  ! go back one line in the file
        end if

       if (DOFIRE) then

        read(861,*) mo,yrin,tmpd(1:nctemdistvars_m),dummy,tilnum

        if (mo == m) then
          do v = 1,nctemdistvars_a ! begin vars loop
            ctem_d_m_mos(v,tilnum,y,m)=tmpd(v)
          end do
        else
         backspace(861)  ! go back one line in the file
        end if

      end if !dofire
     end do !mo=m loop
    end do !months
   end do ! years

     close(741)

    if (DOFIRE) then
      close(861)
      deallocate(tmpd)
    end if
  deallocate(tmp)
!----

! MONTHLY
if (net4) then
  status = nf90_inq_ncid(ncid,'CTEM-Monthly Tiled', grpid)
  if (status /= nf90_noerr) call handle_err(status)
end if

 do v = 1,numctemvars_a ! begin vars loop
  do m=1,12  !Begin Month Loop
   do l = 1,ntile ! begin tile loop

    status = nf90_inq_varid(grpid,trim(CTEM_M_VAR(v)), var_id)
    if (status/=nf90_noerr) call handle_err(status)

    status = nf90_put_var(grpid,var_id,ctem_m_mos(v,l,:,m),start=[xlon,ylat,l,m,yrst],count=[1,1,1,1,totyrs])
    if (status/=nf90_noerr) call handle_err(status)
   end do
  end do
 end do

if (DOFIRE) then

  if (net4) then
    status = nf90_inq_ncid(ncid,'Monthly-Disturbance Tiled', grpid)
    if (status /= nf90_noerr) call handle_err(status)
end if

 do v = 1,numctemvars_a ! begin vars loop
  do m=1,12  !Begin Month Loop
   do l = 1,ntile ! begin tile loop

    status = nf90_inq_varid(grpid,trim(CTEM_M_D_VAR(v)), var_id)
    if (status/=nf90_noerr) call handle_err(status)

    status = nf90_put_var(grpid,var_id,ctem_d_m_mos(v,l,:,m),start=[xlon,ylat,l,m,yrst],count=[1,1,1,1,totyrs])
    if (status/=nf90_noerr) call handle_err(status)
   end do
  end do
 end do

end if !dofire

! deallocate arrays

deallocate(ctem_m_mos)

if (DOFIRE) then
 deallocate(ctem_d_m_mos)
end if

end if ! makemonthly

end if ! check CTEM

!close the netcdf
status = nf90_close(ncid)
if (status/=nf90_noerr) call handle_err(status)

! remove the tmp files

command='rm tmp_a.dat tmp_m.dat tmp_a_p.dat tmp_m_p.dat'

call system(command)

end program netcdf_writer_bf

!=====================

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

!------------------------------
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


