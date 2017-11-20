!>Central module for all general utilities
module generalUtils

    implicit none

    public :: abandonCell
    public :: findDaylength
    public :: findLeapYears
    public :: parseTimeStamp

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
    
    function parseTimeStamp(timeStamp)

    implicit none

    real, dimension(4) :: parseTimeStamp
    real, intent(in) :: timeStamp
    real :: date, moment
    integer :: intdate, day, month, year

    date = floor(timeStamp) ! remove the part days
    parseTimeStamp(4) = timeStamp - date ! save the part days
    intdate=int(date)
    day = mod(intdate, 100);   intdate = intdate / 100
    month = mod(intdate, 100); intdate = intdate / 100
    year = intdate
    parseTimeStamp(1) = real(year)
    parseTimeStamp(2) = real(month)
    parseTimeStamp(3) = real(day)

    end function parseTimeStamp

end module generalUtils


