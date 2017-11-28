module metDisaggModule

    use model_state_drivers, only : metInputTimeStep,metTime,metFss,metFdl,metPre,metTa,metQa,metUv,metPres
    use generalUtils, only : closeEnough,parseTimeStamp
    use ctem_params, only : pi

    implicit none

    integer, allocatable    :: dayIndices(:)                        ! Used by diurnal distribution of shortwave
    integer                 :: tslongStart, tslongEnd, tslongCount
    integer                 :: firstDayIndex, lastDayIndex

contains

    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    ! Processes one single grid cell
    subroutine disaggGridCell(longitude, latitude,delt,lastDOY,nday) ! longitude, latitude

        real, intent(in) :: longitude, latitude
        real, intent(in) :: delt         !< Simulation physics timestep
        integer, intent(in) :: lastDOY
        integer, intent(in) :: nday

        real, dimension(5) :: timestamp
        integer :: numberPhysInMet  !< Number of physics timesteps that fit into one MET input file timestep
        integer :: numberMetInput   !< Number of MET input file timesteps that fit into one day
        integer :: days,i,vcount

        !> First check that we should be doing the disaggregation. If we are already
        !! at the needed timestep (delt) then we can return the main.
        if (closeEnough(metInputTimeStep,delt)) then
           return
        end if

        !> Determine how many physics timesteps fit into a metInputTimeStep
        numberPhysInMet = int(metInputTimeStep / delt)

        !> Determine how many metInputTimeStep fit into a day
        numberMetInput = int(86400./metInputTimeStep)

        !> Next make sure the MET data goes from day 1 to the last day of the year.
        !! This subroutine is setup assuming full years
        timestamp = parseTimeStamp(metTime(1))
        firstDayIndex = timestamp(5)

        !> Find how many days are in the MET file time array
        days = size(metTime) / numberMetInput

        !> Set the lastDayIndex
        lastDayIndex = firstDayIndex + days - 1

        !> Allocate an array with the day indices for the shortwave diurnal distribution
        allocate(dayIndices(days + 2))
        do i = firstDayIndex - 1, lastDayIndex + 1
            dayIndices(i - firstDayIndex + 2) = i
        enddo

        !> Compute timestep range from first and last day
       tslongStart = (firstDayIndex - 1) * numberMetInput + 1
       tslongEnd = lastDayIndex * numberMetInput

        !> Compute the timestep count (both original timestep and the needed one)
        tslongCount = size(metTime) !tslongEnd - tslongStart + 1     ! original timesteps
        vcount = tslongCount * numberPhysInMet

        !> Add two extra days to the vcount since we need them for interpolation of the
        !! first and last days of our time period
        vcount = vcount + 2 * numberPhysInMet * numberMetInput

        !> Expand original MET data to physics timestep. This increases the array size
        !! and puts the old values within the larger array. It also duplicates the first
        !! and last days for the interpolation.
        call makebig(vcount)

        call stepInterpolation(metFdl)
        call linearInterpolation(metTa)
        call linearInterpolation(metQa)
        call linearInterpolation(metUv)
        call linearInterpolation(metPres)
        call precipDistribution(metPre,delt)
        call diurnalDistribution(metFss, latitude, dayIndices)

!         call timeShift(0.25)                    ! 6h correction for faulty NetCDF input data
!         call timeShift(timeZone(lon(lonIndex)))
!         call refreshTimespan
!         call extractValidTime !>> Just remove the added two days of values.

    contains

        !-----------------------------------------------------------------------------------------------------------------------------------------------------
        !> Expands the meteorological variable arrays so they can accomodate the amount of
        !! timesteps on the physics timestep. Also copies the first day and last day values
        !! into the array to give a startday -1 and endday + 1 values for the interpolations.
        subroutine makebig(vcount)

            integer, intent(in) :: vcount
            real, allocatable, dimension(:) :: tmpFss,tmpFdl,tmpPre,tmpTa,tmpQa,tmpUv,tmpPres
            integer :: i,j,timePlusTwoDays

            !> Transfer the present contents of the met*** arrays to temp ones.
            call move_alloc(metFss,tmpFss)
            call move_alloc(metFdl,tmpFdl)
            call move_alloc(metPre,tmpPre)
            call move_alloc(metTa,tmpTa)
            call move_alloc(metQa,tmpQa)
            call move_alloc(metUv,tmpUv)
            call move_alloc(metPres,tmpPres)

            !> Allocate the met*** arrays in preparation to fill with old contents
            !! the allocated size is now including the met on the physics timestep
            allocate(metFss(vcount),metFdl(vcount),metPre(vcount),metTa(vcount), &
                    metQa(vcount),metUv(vcount),metPres(vcount))

            timePlusTwoDays = size(metTime) + 2 * numberMetInput
            do j = 1, timePlusTwoDays
                if (j < numberMetInput+1) then ! Then just copy the first day into that time
                    i = j
                elseif (j > (timePlusTwoDays - numberMetInput)) then ! Copy the last day into that time
                    i = j-(2* numberMetInput)
                else ! use the value for that time
                    i = j - numberMetInput
                end if
                metFss ((j-1) * numberPhysInMet + 1) = tmpFss(i)
                metFdl ((j-1) * numberPhysInMet + 1) = tmpFdl(i)
                metPre ((j-1) * numberPhysInMet + 1) = tmpPre(i)
                metTa  ((j-1) * numberPhysInMet + 1) = tmpTa(i)
                metQa  ((j-1) * numberPhysInMet + 1) = tmpQa(i)
                metUv  ((j-1) * numberPhysInMet + 1) = tmpUv(i)
                metPres((j-1) * numberPhysInMet + 1) = tmpPres(i)
            enddo

            deallocate(tmpFss,tmpFdl,tmpPre,tmpTa,tmpQa,tmpUv,tmpPres)

        end subroutine makebig

        !-----------------------------------------------------------------------------------------------------------------------------------------------------
        ! Step interpolation
        subroutine stepInterpolation(var)
            real, intent(inout)                         :: var(:)
            integer                                     :: i, j, countr

            vcount = size(var)

            do i = 1, countr / numberPhysInMet
                do j = 2, numberPhysInMet
                    var((i - 1) * numberPhysInMet + j) = var((i - 1) * numberPhysInMet + 1)
                enddo
            enddo

        end subroutine stepInterpolation

        !-----------------------------------------------------------------------------------------------------------------------------------------------------
        ! Linear interpolation
        subroutine linearInterpolation(var)
            real, intent(inout)                         :: var(:)
            integer                                     :: i, j, shortIntervals, start, endpt, countr
            real                                        :: t

            countr = size(var)

            ! TODO: Explain this thoroughly FLAG EDUARD
            shortIntervals = countr / numberPhysInMet
            do i = 1, shortIntervals - 1
                var(i * numberPhysInMet) = var(i * numberPhysInMet + 1)
            enddo
            var(countr) = var(countr - numberPhysInMet)
            ! up to here!!!
            do i = 1, shortIntervals
                do j = 0, numberPhysInMet - 1
                    start = (i - 1) * numberPhysInMet + 1
                    endpt = i * numberPhysInMet
                    t = real(j) / numberPhysInMet
                    var(start + j) = var(start) + (var(endpt) - var(start)) * t
                enddo
            enddo

        end subroutine linearInterpolation

        !-----------------------------------------------------------------------------------------------------------------------------------------------------
        ! Precipitation distribution occurs randomly, but conservatively, over the number of wet timesteps.
        subroutine precipDistribution(var,delt)
            real, intent(inout)                         :: var(:)
            integer, intent(in)                         :: delt
            integer                                     :: i, j, k, T,start, endpt, countr,wetpds,sort_ind(numberPhysInMet)
            real                                        :: random(numberPhysInMet),temp
            real                                        :: totP  ! Total precipitation over the metInputTimeStep timestep

            countr = size(var)

            ! Loop through the metInputTimeStep timesteps (commonly 6 hr)
            do i = 1, countr/numberPhysInMet

                start = (i - 1) * numberPhysInMet + 1
                endpt = i * numberPhysInMet

                ! Convert the incoming precipitation from mm/s to mm.
                totP = var(start) * delt

                !> If precipitation for this metInputTimeStep is greater than 0, use formula
                !! to produce number of wet half hours

                if (totP > 0.) then

                    ! We expect precipitation to be in kg/m2/s (or mm/s) for all of this
                    wetpds = nint( &
                                  max( &
                                      min( &
                                          abs(2.6 * log10(6.93 * totP ) ) &
                                          , real(numberPhysInMet)) &
                                                , 1.0) &
                                      )

                    ! Create an array of random numbers
                    call random_number(random)

                    ! Create an array of integers from 1 to numberPhysInMet
                    do k = 1,numberPhysInMet
                        sort_ind(k) = k
                    end do

                    ! Now use the random numbers to create a list of randomly sorted
                    ! integers that will be used as indices later. Here the integers are
                    ! sorted such that the highest integer ends up in the index of the largest
                    ! random number. The sorting continues until integer 1 is in the index position
                    ! of the smallest random number generated.

                    do k = 1,numberPhysInMet
                        do j = k,numberPhysInMet
                            if (random(k) < random(j)) then
                                temp = random(k)
                                random(k) = random(j)
                                random(j) = temp
                                T = sort_ind(k)
                                sort_ind(k) = sort_ind(j)
                                sort_ind(j) = T
                            end if
                        end do
                    end do

                    ! Set all values in random to 0 in preparation for reassignment
                    random = 0.

                    ! Produces random list of 1s and 0s
                    do j = 1, wetpds
                        call random_number(random(sort_ind(j)))
                    end do
                    random = random / sum(random)

                    ! Disperse precipitiation randomly across the random indice
                    k = 1
                    do j = start,endpt
                        ! Assign this time period its precipitation and convert back to mm/s
                        var(j) = random(k) * totP / delt
                        k = k + 1
                    end do

                else !No precip, move on.
                    wetpds = 0
                    var(start:endpt) = 0.
                end if

            enddo

        end subroutine precipDistribution

        !-----------------------------------------------------------------------------------------------------------------------------------------------------
        ! Diurnal distribution over the entire timespan
        subroutine diurnalDistribution(shortWave, latitude, daysOfYear)
            real, intent(inout)                         :: shortWave(:)
            real, intent(in)                            :: latitude
            integer, intent(in)                         :: daysOfYear(:)
            real, allocatable                           :: swDiurnalDistributed(:)
            real, allocatable                           :: correction
            integer                                     :: i, start, endpt, vcount, shortSteps

            vcount = size(shortWave)

            shortSteps = numberPhysInMet * numberMetInput

            allocate(daylightIndices(shortSteps),daylightIndYear(365,shortSteps),&
                     zenithAngYear(365,shortSteps),zenithAngles(shortSteps))

            ! Assume we are going to run the model for at least one full year so just find the
            ! zenith angles and daylight timesteps for one year and store the values.
            do d = 1,365 !FLAG not set up for leap yet!
                 call zenithDaylight(d, latitude, zenithAngles, daylightIndices)
                 daylightIndYear(d,:) = daylightIndices
                 zenithAngYear(d,:) = zenithAngles
            end do
!FLAG here.
            allocate(swDiurnalDistributed(vcount))
            swDiurnalDistributed = 0
            do i = 1, size(shortWave) / 48      ! Once, every day
                start = (i - 1) * 48 + 1
                endpt = i * 48
                            !swMean = sum(swDayValues(1:vcount)) / 4          ! Finds the mean of this day's values
                call applyDiurnalDistribution(zenithAngles, daylightIndices, zenithNoon, swMean, timestep, swDiurnalDistributed(start:endpt))
                !call diurnalDistributionOneDay(shortWave(start:endpt), daysOfYear(i), latitude, swDiurnalDistributed(start:endpt))
            enddo
            shortWave = swDiurnalDistributed
            deallocate(swDiurnalDistributed, dayIndices)

        end subroutine diurnalDistribution

        !-----------------------------------------------------------------------------------------------------------------------------------------------------
        ! Diurnal distribution over just one day
        subroutine zenithDaylight(day, latitude, zenithAngles, daylightIndices)
            !real,       intent(in)  :: swDayValues(:), latitude
            real,       intent(in)  :: latitude
            integer,    intent(in)  :: day
            real,       intent(inout) :: daylightIndices, zenithAngles
            !real                    :: swDiurnalDistributed(:)
            integer                 :: i, dayLength, daylightCount, vcount
            !real                    :: zenithNoon, swMean, latRad
            real                    :: zenithNoon, latRad, zenithCos
            integer,    parameter   :: timestep = 30!, vcount = 48
            !real                    :: daylightIndices(vcount), zenithCos(vcount), timesteps(vcount), zenithAngles(vcount)
            real, allocatable, dimension(:) :: timesteps

            vcount = numberPhysInMet * numberMetInput
            allocate(timesteps(vcount))

            latRad = latitude * pi / 180.00
!             swMean = sum(swDayValues(1:vcount)) / 4          ! Finds the mean of this day's values

            do i = 1, vcount                                    ! Creates the time steps in 30 minutes steps
                timesteps(i) = (i - 1) * 30
            enddo

            ! Find the zenith angles at every specified time step
            call getZenithAngles(timesteps, day, latRad, zenithAngles)

!             ! Find the zenith angle at noon
!             zenithNoon = zenithAngles(25)

            ! Find day length
            daylightCount = 0
            do i = 1, vcount
                zenithCos(i) = cos(zenithAngles(i))
                if (zenithCos(i) >= 0) then
                    daylightCount = daylightCount + 1
                    daylightIndices(i) = i
                else
                    daylightIndices(i) = -1
                endif
            enddo
            dayLength = daylightCount * timestep !purely diagnostic.

            deallocate(timesteps)

            !call applyDiurnalDistribution(zenithAngles, daylightIndices, zenithNoon, swMean, timestep, swDiurnalDistributed)
        end subroutine zenithDaylight

            !-----------------------------------------------------------------------------------------------------------------------------------------------------
        ! Compute zenith angles
        subroutine getZenithAngles(timesteps, day, latitude, zenithAngles)
            integer,    intent(in)          :: day
            real,       intent(in)          :: timesteps(:), latitude
            real,       intent(out)         :: zenithAngles(:)
            integer                         :: i, n(4), vcount
            real                            :: radHourAngle, degHourAngle, psi, dec, a(4), b(4)
            integer,    parameter           :: minutesInHour = 60, daysInYear = 365
            a = (/0.006918, -0.399912, -0.006758, -0.002697/)
            b = (/0.0, 0.070257, 0.000907, 0.001480/)
            n = (/0, 1, 2, 3/)
            vcount = size(timesteps)
            !!!! Does it matter if day > 365?!?! The original code was day = mod(day, 365) or some such thing
            psi = 2 * pi * (day - 1) / daysInYear
            dec = sum((a * cos(n * psi)) + (b * sin(n * psi)))
            !FIND THE HOUR ANGLE, CONVERT IT TO RADIANS THEN FIND THE ZENITH ANGLE(S).
            do i = 1, vcount
                degHourAngle = 15 * (12 - (timesteps(i) / minutesInHour))
                radHourAngle = (degHourAngle / 360) * 2 * pi
                zenithAngles(i) = acos((sin(latitude) * sin(dec)) + &
                (cos(latitude) * cos(dec) * cos(radHourAngle)))
            enddo
        end subroutine getZenithAngles

        !-----------------------------------------------------------------------------------------------------------------------------------------------------
        ! Apply the diurnal distribution on one day
        subroutine applyDiurnalDistribution(zenithAngles, daylightIndices, zenithNoon, &
        swMean, timestep, swDiurnalDistributed)
            integer,    intent(in)          :: timestep
            real,       intent(in)          :: zenithAngles(:), daylightIndices(:), zenithNoon, swMean
            real,       intent(out)         :: swDiurnalDistributed(:)
            real,       allocatable         :: diurnalDistrib(:)
            integer                         :: i, vcount, daylightIndex
            real                            :: vsum, correction, diurnalMean
            real                    :: zenithNoon

            vcount = size(zenithAngles)
            allocate(diurnalDistrib(vcount))

            ! Find the zenith angle at noon
            zenithNoon = zenithAngles(25)

            vsum = 0.0
            do i = 1, vcount
                daylightIndex = daylightIndices(i)              ! CHECK WITH DAYLIGHT_INDICES TO SEE IF IT THERE IS DAYLIGHT AT THE TIME,
                if (daylightIndex > 0) then                     ! IF THERE IS FIND THE VALUE FOR S
                    diurnalDistrib(i) = swMean * pi / 2 * &
                    cos((zenithAngles(daylightIndex) - zenithNoon) / (pi / 2 - zenithNoon) * pi / 2)
                else                                            ! ELSE SET IT TO 0
                    diurnalDistrib(i) = 0
                endif
                vsum = vsum + diurnalDistrib(i)
            enddo

            diurnalMean = vsum / vcount
            if (diurnalMean /= 0) then
                correction = swMean / diurnalMean
            else
                correction = 0
            endif
            swDiurnalDistributed = diurnalDistrib * correction / 21600 ! Joe needs to explain this one

            deallocate(diurnalDistrib)
        end subroutine applyDiurnalDistribution

    end subroutine disaggGridCell

!
!     !-----------------------------------------------------------------------------------------------------------------------------------------------------
!     ! Perform a literal time shift (i.e. shift the time scale only)
!     subroutine timeShift(timeZoneOffset)
!         real, intent(in)                            :: timeZoneOffset
!         character(200)                              :: message
!
!         write(message, "(a,f10.3)") " Performing time shift with offset: ", timeZoneOffset; call notify(message)
!
!         time = time + timeZoneOffset
!
!         call confirm
!
!         call shiftShortWave(timeZoneOffset)
!
!     contains
!         !-----------------------------------------------------------------------------------------------------------------------------------------------------
!         ! Perform time shift on the shortWave variable
!         subroutine shiftShortWave(timeZoneOffset)
!             real,   intent(in)                          :: timeZoneOffset
!             integer                                     :: timeSlots, i
!             character(200)                              :: message
!
!             write(message, "(a,f10.3)") " Performing shortWave shift with offset: ", timeZoneOffset; call notify(message)
!
!             timeSlots = real(timeZoneOffset * 48)
!             if (timeSlots /= 0) then
!                 if (timeSlots > 0) then
!                     do i = 1, size(shortWave) - timeSlots
!                         shortWave(i) = shortWave(i + timeSlots)
!                     enddo
!                 else
!                     timeSlots = -timeSlots
!                     do i = size(shortWave), timeSlots + 1, -1
!                         shortWave(i) = shortWave(i - timeSlots)
!                     enddo
!                 endif
!             endif
!
!             call confirm
!         end subroutine shiftShortWave
!     end subroutine timeShift
!
!     !-----------------------------------------------------------------------------------------------------------------------------------------------------
!     ! Get time zone
!     real function timeZone(longitude)
!         implicit none
!         real,       intent(in)      :: longitude
!         timeZone = real(nint(longitude * 4) / 30 ) / 48
!         if (longitude >= 180) timeZone = timeZone - 0.99999999
!     end function timeZone
!
!     !-----------------------------------------------------------------------------------------------------------------------------------------------------
!     ! Check the timespan for wonky timesteps
!     subroutine refreshTimespan
!         integer :: i
!         real :: date, moment
!
!         do i = 1, size(time)
!             date = floor(time(i))
!             moment = time(i) - date
!             date = refreshDate(date)
!             time(i) = date + moment
!         enddo
!
!     contains
!         !-----------------------------------------------------------------------------------------------------------------------------------------------------
!         ! Refresh date - check for invalid days e.g. 19011232 or 19010001
!         real function refreshDate(dateIn)
!             real, intent(in)                            :: dateIn
!             integer                                     :: i, year = 0, month = 0, day = 0, date = 0
!
!             date = int(dateIn)
!             day = mod(date, 100);   date = date / 100
!             month = mod(date, 100); date = date / 100
!             year = date
!
!             if (day > daysInMonth(month)) then
!                 month = month + 1
!                 call checkMonth(year, month)
!                 day = 1
!             endif
!             if (day < 1) then
!                 month = month - 1
!                 call checkMonth(year, month)
!                 day = daysInMonth(month)
!             endif
!
!             call checkMonth(year, month)
!
!             date = year * 10000 + month * 100 + day
!             refreshDate = date
!             !print*, date, day, month, year
!         end function refreshDate
!
!         !-----------------------------------------------------------------------------------------------------------------------------------------------------
!         ! Refresh the month
!         subroutine checkMonth(year, month)
!             integer, intent(inout) :: year, month
!             if (month > 12) then
!                 month = 1
!                 year = year + 1
!             endif
!             if (month < 1) then
!                 month = 12
!                 year = year - 1
!             endif
!         end subroutine checkMonth
!     end subroutine refreshTimespan
!

!
!     !-----------------------------------------------------------------------------------------------------------------------------------------------------
!     ! Find timestep by key
!     integer function findTimestep(key)
!         real, intent(in)    :: key
!         integer                         :: i
!         findTimestep = -1
!         do i = 1, size(time)
!             if (closeEnough(key, time(i))) then
!                 findTimestep = i
!             endif
!         end do
!
!
!     end function findTimestep
!
!     !-----------------------------------------------------------------------------------------------------------------------------------------------------
!     ! Deallocates all variables
!     subroutine unloadVariables
!         if (allocated(time)) deallocate(time)
!     end subroutine unloadVariables
!
!     !-----------------------------------------------------------------------------------------------------------------------------------------------------
!     ! Sets the current timespan
!     subroutine setTimespan(startDate, endDate, firstDayIndex, lastDayIndex)
!         integer, optional, intent(in)       :: startDate, endDate, firstDayIndex, lastDayIndex
!
!         timespan%firstDayIndex = firstDayIndex
!         timespan%lastDayIndex = lastDayIndex
!
!         ! afterwards, update startDate and endDate from the NetCDF accordingly
!     end subroutine setTimespan
!
!
!     !-----------------------------------------------------------------------------------------------------------------------------------------------------
!     ! Logs and/or notifies the status of the processes
!     subroutine notify(message, override)
!         character(*), intent(in)        :: message
!         logical, optional, intent(in)   :: override
!
!         if (verbose) then
!             print "(a,a$)", trim(message), "..."
!         else
!             if (present(override) .and. override) print"(a,a$)", trim(message), "..."
!         endif
!     end subroutine notify
!
!     !-----------------------------------------------------------------------------------------------------------------------------------------------------
!     ! Prints [Done] at the end of individual the process
!     subroutine confirm(override)
!         logical, optional, intent(in)   :: override
!
!         if (verbose) then
!             print *, '['//achar(27)//'[32mDone'//achar(27)//'[0m]'
!         else
!             if (present(override) .and. override) print *, '['//achar(27)//'[32mDone'//achar(27)//'[0m]'
!         endif
!     end subroutine confirm

end module metDisaggModule
