!>Central module for all general utilities
module generalUtils

    implicit none

    public :: abandonCell
    public :: findDaylength
    public :: findCloudiness
    public :: findLeapYears
    public :: findActiveLayerDepth
    public :: parseTimeStamp
    public :: closeEnough

    logical :: run_model           !< Simple logical switch to either keep run going or finish

    contains

    !---------------------------------------------------------------------------------------

    subroutine abandonCell

        use class_statevars,    only : class_rot

        implicit none

        real, pointer, dimension(:) :: DLONROW !<
        real, pointer, dimension(:) :: DLATROW !<

        DLATROW => class_rot%DLATROW
        DLONROW => class_rot%DLONROW

        run_model = .false.
        print*,'died on',DLONROW,DLATROW
        return

    end subroutine abandonCell

    !---------------------------------------------------------------------------------------

    !> Calculate the daylength based on the latitude and day of year
    real function findDaylength(solday,radl)

        ! Joe Melton Dec 18 2015 (taken from phenlogy.f)

        use ctem_params, only : pi

        implicit none

        real, intent(in) :: solday  !day of year
        real, intent(in) :: radl    ! latitude
        real :: theta               ! temp var
        real :: decli               ! temp var
        real :: term                ! temp var

            theta=0.2163108 + 2.0*atan(0.9671396*tan(0.0086*(solday-186.0)))
            decli=asin(0.39795*cos(theta))      !declination !note I see that CLASS does this also but with different formula...
            term=(sin(radl)*sin(decli))  /(cos(radl)*cos(decli))
            term=max(-1.0,min(term,1.0))
            findDaylength=24.0-(24.0/pi)*acos(term)

    end function findDaylength

    !---------------------------------------------------------------------------------------

    subroutine findLeapYears(iyear,leapnow,lastDOY)

        !Check if this year is a leap year

        use ctem_params,        only : monthend, mmday,monthdays

        implicit none

        logical, intent(out) :: leapnow
        integer, intent(in) :: iyear
        integer, intent(inout) :: lastDOY

        if (mod(iyear,4).ne.0) then !it is a common year
            leapnow = .false.
        else if (mod(iyear,100).ne.0) then !it is a leap year
            leapnow = .true.
        else if (mod(iyear,400).ne.0) then !it is a common year
            leapnow = .false.
        else !it is a leap year
            leapnow = .true.
        end if

        if (leapnow) then
            lastDOY = 366
        else
            lastDOY = 365
        end if

        ! We do not check the MET files to make sure the incoming MET is in fact
        ! 366 days if leapnow. You must verify this in your own input files. Later
        ! in the code it will fail and print an error message to screen warning you
        ! that your file is not correct.
        if (leapnow) then ! adjust the calendar and set the error check.
            monthdays = (/ 31,29,31,30,31,30,31,31,30,31,30,31 /)
            monthend = (/ 0,31,60,91,121,152,182,213,244,274,305,335,366 /)
            mmday = (/ 16,46,76,107,137,168,198,229,260,290,321,351 /)

        else
            monthdays = [ 31,28,31,30,31,30,31,31,30,31,30,31 ] !< days in each month
            monthend  = [ 0,31,59,90,120,151,181,212,243,273,304,334,365 ] !< calender day at end of each month
            mmday     = [ 16,46,75,106,136,167,197,228,259,289,320,350 ] !<mid-month day
        end if

    end subroutine findLeapYears

    !---------------------------------------------------------------------------------------
    !> The cosine of the solar zenith angle COSZ is calculated from the day of
    !> the year, the hour, the minute and the latitude using basic radiation geometry,
    !> and (avoiding vanishingly small numbers) is assigned to CSZROW.  The fractional
    !> cloud cover FCLOROW is commonly not available so a rough estimate is
    !> obtained by setting it to 1 when precipitation is occurring, and to the fraction
    !> of incoming diffuse radiation XDIFFUS otherwise (assumed to be 1 when the sun
    !> is at the horizon, and 0.10 when it is at the zenith).

    subroutine findCloudiness(nltest,imin,ihour,iday,lastDOY)

        use ctem_params, only : pi
        use class_statevars, only : class_rot

        implicit none

        integer, intent(in) :: nltest
        integer, intent(in) :: imin
        integer, intent(in) :: ihour
        integer, intent(in) :: iday
        integer, intent(in) :: lastDOY

        integer :: i
        real :: day
        real :: decl    !< Declination
        real :: hour
        real :: cosz    !< Cosine of the zenith angle

        real, pointer, dimension(:) :: RADJROW !<Latitude of grid cell (positive north of equator) [rad]
        real, pointer, dimension(:) :: CSZROW  !<Cosine of solar zenith angle [ ]
        real, pointer, dimension(:) :: PREROW  !< Surface precipitation rate \f$[kg m^{-2} s^{-1} ]\f$
        real, pointer, dimension(:) :: XDIFFUS !< Fraction of diffused radiation
        real, pointer, dimension(:) :: FCLOROW !<Fractional cloud cover [ ]

        RADJROW => class_rot%RADJROW
        CSZROW => class_rot%CSZROW
        PREROW => class_rot%PREROW
        XDIFFUS => class_rot%XDIFFUS
        FCLOROW => class_rot%FCLOROW

        day=real(iday)+(real(ihour)+real(imin)/60.)/24.
        decl=sin(2.*pi*(284.+day)/real(lastDOY))*23.45*pi/180.
        hour=(real(ihour)+real(imin)/60.)*pi/12.-pi

        do i=1,nltest

            cosz=sin(radjrow(i))*sin(decl)+cos(radjrow(i))*cos(decl)*cos(hour)

            cszrow(i)=sign(max(abs(cosz),1.0e-3),cosz)
            if(prerow(i).gt.0.) then
                xdiffus(i)=1.0
            else
                xdiffus(i)=max(0.0,min(1.0-0.9*cosz,1.0))
            endif
            fclorow(i)=xdiffus(i)
        end do

    end subroutine findCloudiness

    !---------------------------------------------------------------------------------------
    !> Parses a time stamp in the expected form "day as %Y%m%d.%f"
    !! Returns an array with 1) year, 2) month, 3) day, 4) fraction of day
    !! 5) day of year
    function parseTimeStamp(timeStamp)

    use ctem_params, only : monthdays

    implicit none

    real, dimension(5) :: parseTimeStamp
    real, intent(in) :: timeStamp
    real :: date, moment
    integer :: intdate, day, month, year,totdays,t

    date = floor(timeStamp) ! remove the part days
    parseTimeStamp(4) = timeStamp - date ! save the part days
    intdate=int(date)
    day = mod(intdate, 100);   intdate = intdate / 100
    month = mod(intdate, 100); intdate = intdate / 100
    year = intdate
    parseTimeStamp(1) = real(year)
    parseTimeStamp(2) = real(month)
    parseTimeStamp(3) = real(day)
    totdays=0
    if (month > 1) then
        do t = 1, month-1
            totdays = totdays + monthdays(t)
        end do
    end if
    parseTimeStamp(5) = real(totdays + day)

    end function parseTimeStamp

    !---------------------------------------------------------------------------------------

    subroutine findActiveLayerDepth


    end subroutine findActiveLayerDepth

    !---------------------------------------------------------------------------------------

    !> As real numbers are not precise, this is a simple way to compare two reals
    logical function closeEnough(num1, num2,error)
        real, intent(in)    :: num1, num2
        real, intent(in)     :: error
        if (abs(num1 - num2) < error) then
            closeEnough = .true.
        else
            closeEnough = .false.
        endif
    end function closeEnough

end module generalUtils
