module met_forcing

implicit none

public :: openmet_nc
public :: readin_met_nc

public :: disagg_6hr_to_hh
private :: linear_dist
private :: rad_diurnal_dist
private :: uniform_dist
private :: precip_dist

contains

subroutine openmet_nc()

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
write(*,*)met_file
call check_nc(nf90_open(trim(met_file),nf90_nowrite,metfid))

!retrieve dimensions

!call check_nc(nf90_inq_dimid(metfid,'lon',dimid))
!call check_nc(nf90_inquire_dimension(metfid,dimid,len=xsize))
!call check_nc(nf90_inq_dimid(metfid,'lat',dimid))
!call check_nc(nf90_inquire_dimension(metfid,dimid,len=ysize))
call check_nc(nf90_inq_dimid(metfid,'time',dimid))
call check_nc(nf90_inquire_dimension(metfid,dimid,len=timelen))

!retrieve scale factors and offsets

do i = 1,niv

  call check_nc(nf90_inq_varid(metfid,vname(i),varid))

    ! If the netcdf has scale factors or offsets, you will need to read those in here
    ! and make some vars to store that info.
end do

!calculate the number and indices of the pixels to be calculated

!allocate(lonvect(xsize))
!allocate(latvect(ysize))

! Figure out the size of the area to be simulated

! call check_nc(nf90_inq_varid(metfid,'lon',varid))
! call check_nc(nf90_get_var(metfid,varid,lonvect))
!
! call check_nc(nf90_inq_varid(metfid,'lat',varid))
! call check_nc(nf90_get_var(metfid,varid,latvect))

call check_nc(nf90_inq_varid(metfid,'time',varid))
call check_nc(nf90_get_var(metfid,varid,timestart,start=(/1,1,1/),ct=(/1,1,1/)))
call check_nc(nf90_get_att(metfid,varid,'units',unitchar))

write(*,*)'The metdata in your netcdf file says it starts:',timestart,'for units:',unitchar,&
         &'. The expected units are: day as %Y%m%d.%f.'

yearmetst = int((timestart(1,1,1) - 101.)/10000.) ! This is based on the MET file starting on the first day of the first month at 0:00!


end subroutine openmet_nc

! ! ================================================================================

! subroutine readin_met_nc(start_pos,dlat,dlon)
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
!     call check_nc(nf90_get_var(metfid,varid,var_in,start=[srtx,srty,start_pos],ct=[cntx,cnty,tlen]))  ! FLAG check the y,x order here.
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
! end subroutine readin_met_nc

!--------------------------------------------------------------------------------------------

subroutine disagg_6hr_to_hh(invar, reallat, reallong, day, fss_hh, fdl_hh, pre_hh,&
                        ta_hh, qa_hh, uv_hh, pres_hh)

    implicit none

    ! Arguments
    real, dimension(:,:), intent(in)  :: invar  ! expected dimension 5,7
    integer, intent(in)            :: day       !< Day of year
    integer, intent(in)            :: lcltmind  !< ?
    real, intent(in)               :: reallat   !< Latitude in radians
    real, intent(in)               :: reallong  !< Longitude in radians

    real, dimension(48), intent(out) :: fss_hh  !< Shortwave
    real, dimension(48), intent(out) :: fdl_hh
    real, dimension(48), intent(out) :: pre_hh
    real, dimension(48), intent(out) :: ta_hh
    real, dimension(48), intent(out) :: qa_hh
    real, dimension(48), intent(out) :: uv_hh
    real, dimension(48), intent(out) :: pres_hh

    ! Local vars
    real, dimension(4) :: SW4

    ! --------

    ! Short Wave Radiation
    ! 21600 is the number of seconds in 6 hours.
    SW4 = invar(:,1)/21600.00

    call rad_diurnal_dist(SW4, reallat, day, reallong, fss_hh, lcltmind)

    ! Long-Wave Radiation

    call uniform_dist(invar(1,2),fdl_hh(1:12))
    call uniform_dist(invar(2,2),fdl_hh(13:24))
    call uniform_dist(invar(3,2),fdl_hh(25:36))
    call uniform_dist(invar(4,2),fdl_hh(37:48))

    ! Precipitation

    call precip_dist(invar(1,3),pre_hh(1:12))
    call precip_dist(invar(2,3),pre_hh(13:24))
    call precip_dist(invar(3,3),pre_hh(25:36))
    call precip_dist(invar(4,3),pre_hh(37:48))

    ! Temperature

    call linear_dist(invar(1,4)-273.16,invar(2,4)-273.16,ta_hh(1:12))  !FLAG check if I need to do the shift to C.
    call linear_dist(invar(2,4)-273.16,invar(3,4)-273.16,ta_hh(13:24))
    call linear_dist(invar(3,4)-273.16,invar(4,4)-273.16,ta_hh(25:36))
    call linear_dist(invar(4,4)-273.16,invar(5,4)-273.16,ta_hh(37:48))

    ! Specific Humidity

    call linear_dist(invar(1,5),invar(2,5),qa_hh(1:12))
    call linear_dist(invar(2,5),invar(3,5),qa_hh(13:24))
    call linear_dist(invar(3,5),invar(4,5),qa_hh(25:36))
    call linear_dist(invar(4,5),invar(5,5),qa_hh(37:48))

    ! Wind Speed

    call linear_dist(invar(1,6),invar(2,6),uv_hh(1:12))
    call linear_dist(invar(2,6),invar(3,6),uv_hh(13:24))
    call linear_dist(invar(3,6),invar(4,6),uv_hh(25:36))
    call linear_dist(invar(4,6),invar(5,6),uv_hh(37:48))

    ! Atmospheric Pressure

    call linear_dist(invar(1,7),invar(2,7),pres_hh(1:12))
    call linear_dist(invar(2,7),invar(3,7),pres_hh(13:24))
    call linear_dist(invar(3,7),invar(4,7),pres_hh(25:36))
    call linear_dist(invar(4,7),invar(5,7),pres_hh(37:48))

! !!!!!!!!!!!!!!!!!!!!!!!! SHIFTING VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!              do 140 count = 1, quantity
!               shift = count + (48.0 - lcltmind - 6.0)
!               if(shift.gt.quantity)shift = shift - quantity
!                 LW_shifted_num(count) = LW_interpolated_num(shift)
!                 PCP_shifted_num(count) = PCP_interpolated_num(shift)
!                 T_shifted_num(count) = T_interpolated_num(shift)
!                 H_shifted_num(count) = H_interpolated_num(shift)
!                 W_shifted_num(count) = W_interpolated_num(shift)
!                 P_shifted_num(count) = P_interpolated_num(shift)
!
! 140          continue

    return

end subroutine disagg_6hr_to_hh

! -------------------------------------------

subroutine linear_dist(val1, val2, outval)

    ! Linearly distribute the quantity over the timesteps

    implicit none

    !Arguments
    real,intent(in) :: val1
    real,intent(in) :: val2
    real,dimension(12), intent(out) :: outval

    ! Local vars
    integer        :: z

    !----

    do z = 1,11
        outval(z) = val1 + ((val2 - val1) * (x(z) / 12))
    end do

    return

end subroutine linear_dist

! -------------------------------------------

subroutine rad_diurnal_dist(values,rad_lat,day,longt,zen_s_fixed,lcltmind)

    implicit none

    ! Arguments
    real, dimension(4), intent(in)  :: values
    real, intent(in)    :: rad_lat
    real, intent(in)    :: longt
    integer, intent(in) :: day
    integer, intent(out) :: lcltmind
    real, allocatable(:), intent(out) :: zen_s_fixed

    ! Local vars
    real :: sum1
    real :: mean_s
    integer :: ct
    integer :: daytim
    integer :: day_length
    integer :: time, K
    integer :: ts_tot   !< Number of timesteps in one day, if model running on 30 min timestep then this is 48.
    real, dimension(1) :: zenith_noon       !< Zenith angle at noon

    ! Allocatable vars
    integer, allocatable(:) :: daytim_indx
    real, allocatable(:) :: times
    real, allocatable(:) :: zenith_cos
    real, allocatable(:) :: zenith_values

    ! Parameters
    integer, parameter :: day_minutes = 1440  !< 24 hours in minutes
    integer, parameter :: time_step = 30  !INPUT !< Minutes in a CLASS timestep
    integer, parameter :: num_of_values = 4
    real, parameter    :: noon(1) = 720     !< Noon in minutes


    ! ---------------

    ! Finds the mean of the input values

    ts_tot = day_minutes / time_step

    allocate(times(ts_tot))
    allocate(zenith_cos(ts_tot))
    allocate(zenith_values(ts_tot))
    allocate(zen_s_fixed(ts_tot))
    allocate(daytim_indx(ts_tot))

    sum1 = sum(values)
    mean_s = sum1/num_of_values

    do t = 0, ts_tot-1
        times(t) = time_step*t
    end do

    !> Get zenith angles
    call get_zenith(times, ts_tot, day, rad_lat, zenith_values)

    !> Do I need this or can I just use the one from above?
    call get_zenith(noon, 1, day, rad_lat, zenith_noon)

    ! Find the day length

    daytim = 0
    do ct = 1, ts_tot
        zenith_cos(ct) = cos(zenith_values(ct))
        if (zenith_cos(ct) .ge. 0.0) then
            daytim = daytim + 1
            daytim_indx(ct) = ct
        else
            daytim_indx(ct)=-1
        endif
    end do
    day_length = daytim * time_step  !int or real?

    ! Find diurnal distribution of sw radiation at specified time step

    call get_zen_s(zenith_values, ts_tot, daytim, daytim_indx, &
                zenith_noon, mean_s, longt,time_step, zen_s_fixed,lcltmind)

    deallocate(times) !Should I deallocate or will they automatically upon return?
    deallocate(zenith_cos)
    deallocate(zenith_values)
    deallocate(zen_s_fixed)
    deallocate(daytim_indx)

    return

end subroutine rad_diurnal_dist

! -------------------------------------------

subroutine get_zenith(times, qty, day, rad_lat,out_zenith)

    use ctem_params, only :: pi

    implicit none

    ! Arguments
    real, dimension(:), intent(in) :: times
    integer, intent(in) :: qty
    integer, intent(in) :: day
    real, intent(in)    :: rad_lat
    real, allocatable, intent(out) :: out_zenith

    ! Local vars
    integer :: ct
    real :: psi
    real :: rad_hour_angle
    real :: deg_hour_angle
    real :: dec

    ! Parameters
    real, dimension(4), parameter :: A = [0.006918, -0.399912, -0.006758, -0.002697]
    real, dimension(4), parameter :: B = [0.0,  0.070257,  0.000907,  0.001480]
    real, dimension(4), parameter :: N = [0.,1.,2.,3.]
    integer, parameter :: minutes_in_hour = 60
    integer, parameter :: days_in_year = 365  !FLAG

    ! ---------

    allocate(out_zenith(qty))

    psi = (2 * pi * (day - 1)) / days_in_year
    dec = 0
    do t = 1,4
        dec = dec + A(t) * cos(N(t) * psi) + (B(t) * sin(N(t) * psi)
    end do

!   Find the hour angle, convert it to radians then find the zenith
!   angle(s).

    do ct = 1, qty
        deg_hour_angle = 15. * (12. - (times(ct) / minutes_in_hour))
        rad_hour_angle = (deg_hour_angle/ 360.) * 2. * pi
        out_zenith(ct) = acos((sin(rad_lat) * sin(dec)) + &
                         (cos(rad_lat) * cos(dec) * cos(rad_hour_angle)))
    end do

    deallocate(out_zenith)

    return

end subroutine get_zenith

! -------------------------------------------

subroutine get_zen_s(zenith_values,qty,raw_size,daytim_indx, &
                    zenith_noon, mean_s, longt, time_step, zen_s,lcltmind)

    use ctem_params, only :: pi

    implicit none

    ! Arguments
    real, dimension(:), intent(in) :: zenith_values
    integer, dimension(:), intent(in) :: daytim_indx
    integer, intent(in) :: qty
    real, intent(in) :: mean_s
    real, intent(in) :: longt
    integer, intent(in) :: time_step
    real, dimension(1), intent(in) :: zenith_noon
    real, dimension(:), allocatable, intent(out) :: zen_s
    integer, intent(out) :: lcltmind

    ! Allocatable vars
    real, dimension(:), allocatable  :: s_raw

    ! Local vars
    integer :: ct,K
    integer :: shift
    real :: mean_s_raw
    real :: sum1
    real :: correction
    integer :: local_time

    ! Parameters
    integer, parameter :: intzero = 0

    ! ----------

    allocate(s_raw(qty))
    allocate(zen_s(qty))

    do ct = 1, qty

        ! Check with daytim_indx to see if it there is daylight at
        ! the time, and if there is, then find the value for s, setting
        ! it to zero otherwise.
        if (daytim_indx(ct) .gt. intzero) then
            s_raw(ct) = ((mean_s * pi) / 2.) * cos(((zenith_values(ct)  &
                        - zenith_noon(1)) / ((pi / 2.) - zenith_noon(1))) * (pi / 2.))
        else
            s_raw(ct) = 0.
        end if
    end do

    sum1 = sum(s_raw)
    mean_s_raw = sum1/qty

    if (mean_s_raw .ne. 0.) then
        correction = mean_s / mean_s_raw
    else
        correction = 0.
    end if

    ! Find what's the local time when GMT equals midnight, assuming that
    ! every 15 degree longitude equals 1 hour.

    local_time = nint(longt*4.0) ! Local time in minutes


    ! Find which time slot index does this time fits in

    lcltmind = 0 ! Local time index
    do K = 0, qty-1
        if(local_time .ge. K*time_step .and. local_time .lt. (K+1)*time_step )then
            lcltmind = K+1
        end if
    end do

    if(lcltmind.eq.0)then
        write(6,*)'LOCAL TIME INDEX ZERO'
!        PAUSE
    end if

    ! Apply the correction and shift so that numbers start at GMT midnight
! FLAG this can't be right below!
    do ct = 1, qty
        shift = ct + lcltmind
        if (shift .gt. qty) shift = shift - qty
        zen_s(ct) = s_raw(ct) * correction
    end do

    deallocate(s_raw)
    deallocate(zen_s)

    return

end subroutine get_zen_s

! -------------------------------------------

subroutine uniform_dist(val1, outvals)

    ! Uniformly distribute the quantity over the timesteps so all get the same value

    implicit none

    ! Arguments
    real,intent(in) :: val1
    real,dimension(12),intent(out) :: outvals

    ! Local vars
    integer :: z

    ! ----

    do z = 1,12
        outvals(z) = val1
    end do

    return

end subroutine uniform_dist

! -------------------------------------------

subroutine precip_dist(totpre,ppt)

    !> Distribute the precipitation over the timesteps randomly but conservatively following
    !! \cite Arora1997-ll.

    implicit none

    ! Arguments
    real,intent(in) :: totpre
    real,dimension(12),intent(out) :: ppt

    ! Local vars
    real :: wet_hh
    real :: temp,sum1
    real,dimension(12) :: random
    integer :: i,k,t
    integer :: wet
    integer,dimension(12) :: sort_ind

    ! ------------

    !P = totpre*6.*60.*60. ! convert mm/s to mm/6hr FLAG check if needed.
    !P = max(P,0.0)

    if (totpre .gt. 0.0) then
        wet_hh = 2.6 * log10(6.93 * totpre)
        wet_hh = max(1.0, min(wet_hh, 12.0))
        wet = nint(wet_hh)
    else
        wet_hh = 0.0
        wet = 0
    end if

    do k = 1,12
        random(k) = rand()
        sort_ind(k) = k
    end do

    do i = 1,12
        do k = i,12
            if(random(i) .lt. random(k)) then
                ! swap the values
                temp = random(i)
                random(i) = random(k)
                random(k) = temp
                t = sort_ind(i)
                sort_ind(i) = sort_ind(k)
                sort_ind(k) = t
            end if
        end do
    end do

    ! reset the values of the random array
    random(:) = 0.0

    sum1 = 0.0
    do k = 1,wet
        random(sort_ind(k)) = rand()
        sum1 = sum1 + random(sort_ind(k))
    end do

    do k = 1,12
        i = mod(k,12)
        if(i.eq.0) i = 12
        if(sum1.gt.0.0) then
            ppt(k) = (random(i)/sum1)*totpre*(1./(0.5*60*60)) !why 0.5? this due to timestep?
        else
            ppt(k)=0.0
        end if
    end do

    return

end subroutine precip_dist

! -------------------------------------------

end module met_forcing