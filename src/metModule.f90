module metDisaggModule

    use model_state_drivers, only : metInputTimeStep,metTime
    use generalUtils, only : closeEnough,parseTimeStamp

    implicit none


    integer, allocatable    :: dayIndices(:)                    ! Used by diurnal distribution of shortwave
    integer                 :: tslongStart, tslongEnd, tslongCount    ! 6h timesteps start, end and count
    integer                     :: firstDayIndex, lastDayIndex
!     ! All the following variables are at a grid level

!     type(timespanType)          :: timespan
!     logical                     :: verbose = .true.
!     real, parameter             :: pi = 3.14159265
!     integer, parameter          :: daysInMonth(*) = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
!
!     ! All the following variables are at a grid cell level only
!     integer                 :: lonIndex, latIndex
!     real, allocatable       :: lon(:), lat(:), time(:), timeStart, timeEnd
!     real, allocatable       :: shortWave(:), longWave(:), precipitation(:), temperature(:), humidity(:), wind(:), pressure(:)
!     integer                 :: daysInInputDataset = 0, days, count
!     logical                 :: includesFirstDay, includesLastDay! Are we using the first or last day of the data?
!     real, parameter         :: fillValue = 1.E38
!
!     ! NetCDF variables
!     integer                 :: inputId

contains
    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    ! Processes one single grid cell
    subroutine disaggGridCell(lonIndex, latIndex,delt,lastDOY,nday) ! longitude, latitude

        integer, intent(in) :: lonIndex, latIndex
        real, intent(in) :: delt         !< Simulation physics timestep
        integer, intent(in) :: lastDOY
        integer, intent(in) :: nday

        real, dimension(5) :: timestamp
        integer :: numberPhysInMet  !< Number of physics timesteps that fit into one MET input file timestep
        integer :: numberMetInput   !< Number of MET input file timesteps that fit into one day
        integer :: days,i

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
        if (closeEnough(timestamp(5),1.)) then
            firstDayIndex = 1
        else
            print*,'Problem in metDisaggModule, expecting MET to start on first day of year'
            stop
        end if

        timestamp = parseTimeStamp(metTime(ubound(metTime,1)))
        if (closeEnough(timestamp(5),real(lastDOY))) then
            lastDayIndex = lastDOY
        else
            print*,'Problem in metDisaggModule, last day of MET ',timestamp(5),', is not ',lastDOY
            stop
        end if

        !> Find how many days are in the MET file time array
        days = size(metTime) / numberMetInput
        print*,days
        ! Allocate an array with the day indices for the shortwave diurnal distribution
        allocate(dayIndices(days + 2))
        do i = firstDayIndex - 1, lastDayIndex + 1
            dayIndices(i - firstDayIndex + 2) = i
        enddo

        ! Compute timestep range from first and last day
       tslongStart = (firstDayIndex - 1) * numberMetInput + 1
       tslongEnd = lastDayIndex * numberMetInput

        ! Compute the timestep count (both 6h and hh)
        tslongCount = tslongEnd - tslongStart + 1     ! 6h timesteps
        count = tslongCount * 12                  ! x12 half hour intervals
!
!         ! Load dimensions
!         call loadTime
!
!         call processVariables
!         call timeShift(0.25)                    ! 6h correction for faulty NetCDF input data
!         call timeShift(timeZone(lon(lonIndex)))
!         call refreshTimespan
!         call extractValidTime


    end subroutine disaggGridCell

!     !-----------------------------------------------------------------------------------------------------------------------------------------------------
!     ! Loads all the relevant dimensions and variables for the specified timespan inside a grid cell
!     subroutine loadTimespan
!
!         ! Check that the first and last day are valid
!         call checkTimespan
!
!         ! Sets up the timestep environment
!         call initializeTimesteps
!
!         ! Load dimensions
!         call loadTime
!         call loadDimension('lon')
!         call loadDimension('lat')
!
!         ! Load variables
!         shortWave = loadVariable("sw")
!         longWave = loadVariable("lw")
!         precipitation = loadVariable("pr")
!         humidity = loadVariable("qa")
!         wind = loadVariable("wi")
!         pressure = loadVariable("ap")
!         temperature = loadVariable("ta")
!
!         ! Convert from Kelvin to Celsius
!         temperature = temperature - 273.15
!
!     contains
!         !-----------------------------------------------------------------------------------------------------------------------------------------------------
!         ! Checks to make sure that the desired timespan is valid
! !         subroutine checkTimespan
! !             integer                     :: firstDayIndex, lastDayIndex
! !             character(200)              :: message
! !             integer                     :: inputTimeSteps
! !
! !             ! Load time span boundaries
! !             firstDayIndex = timespan%firstDayIndex; lastDayIndex = timespan%lastDayIndex
! !
! !             write(message, "(a$)") " Verifying the validity of the timespan "
! !             call notify(message)
! !
! !             ! Check that firstDayIndex is valid
! !             if (firstDayIndex < 1) stop("First day must be 1 or greater")
! !
! !             ! Stop if the last day is before the first day
! !             if (lastDayIndex < firstDayIndex) stop("The last day of the interval can't be earlier than the first day")
! !
! !             ! Find out how much data there is
! !             inputTimeSteps = ncGetDimLen(inputId,'time')
! !
! !             ! Find the size of the dataset in days
! !             daysInInputDataset = inputTimeSteps / 4
! !
! !             ! Check that the last day is available
! !             if (lastDayIndex > daysInInputDataset) stop("The last day of the timespan is later than the existing data")
! !
! !             ! Set first day flag
! !             if (firstDayIndex == 1) includesfirstDay = .true.
! !
! !             ! Set last day flag
! !             if (lastDayIndex == daysInInputDataset) includesLastDay = .true.
! !
! !             ! Set the number of days
! !             days = lastDayIndex - firstDayIndex + 1
! !
! !             call confirm
! !         end subroutine checkTimespan
!
!         !-----------------------------------------------------------------------------------------------------------------------------------------------------
!         ! Initializes the timesteps
!         subroutine initializeTimesteps
!             integer                 :: firstDayIndex, lastDayIndex
!             character(200)          :: message
!             integer                 :: i
!
!             ! Load time span boundaries
!             firstDayIndex = timespan%firstDayIndex; lastDayIndex = timespan%lastDayIndex
!
!             ! Allocate an array with the day indices for the shortwave diurnal distribution
!             allocate(dayIndices(days + 2))
!             do i = firstDayIndex - 1, lastDayIndex + 1
!                 dayIndices(i - firstDayIndex + 2) = i
!             enddo
!
!             ! Compute timestep range from first and last day
!             tslongStart = (firstDayIndex - 1) * 4 + 1
!             tslongEnd = lastDayIndex * 4
!
!             ! If there's enough extra data, get the real thing (needed for timezone shift)
!             if (firstDayIndex > 1) tslongStart = tslongStart - 4                    ! grab an extra 4 x 6h intervals from before
!             if (lastDayIndex <= daysInInputDataset - 1) tslongEnd = tslongEnd + 4   ! grab an extra 4 x 6h intervals from after
!
!             ! Compute the timestep count (both 6h and hh)
!             tslongCount = tslongEnd - tslongStart + 1     ! 6h timesteps
!             count = tslongCount * 12                  ! x12 half hour intervals
!
!             write(message, "(a,i7,a,i7,a)") " Computing", count, " 1/2h timesteps from ", tslongCount, " 6h timesteps"
!             call notify(message); call confirm
!         end subroutine initializeTimesteps
!
        !-----------------------------------------------------------------------------------------------------------------------------------------------------
        ! Load the time dimension
!         subroutine loadTime
!             character(200)                              :: message
!             real, allocatable                           :: time6h(:), rawNCData(:)
!             integer                                     :: timeCount, varId, i
!             real                                        :: day, moment
!
!             ! Load existing data
!             timeCount = ncGetDimLen(inputId, 'time')
!             write(message, "(a,i10,a,i8,a)") " Loading", tslongCount, " out of ", timeCount, " time values from NetCDF file..."
!             call notify(message)
!
!             allocate(rawNCData(tslongCount))
!             rawNCData = ncGetDimValues(inputId, 'time',  start=[tslongStart], count=[tslongCount])
!             call confirm
!
!             write(message, "(a)") " Filling in missing days (first/last)"; call notify(message)
!             ! Fill in missing days
!             if (includesFirstDay) then
!                 if (includesLastDay) then                    ! Need to fill in both the first and last days
!                     allocate(time6h(tslongCount + 8))
!                     time6h(1:4) = rawNCData(1:4) - 1
!                     time6h(5:tslongCount + 4) = rawNCData
!                     time6h(tslongCount + 5 : tslongCount + 8) = rawNCData(tslongCount - 3 : tslongCount) + 1
!                     count = count + 48 * 2
!                 else                                        ! Need to fill in just the first day
!                     allocate(time6h(tslongCount + 4))
!                     time6h(5:tslongCount + 4) = rawNCData
!                     time6h(1:4) = rawNCData(1:4) - 1
!                     count = count + 48
!                 endif
!             else
!                 if (includesLastDay) then                    ! Need to fill in just the last day
!                     allocate(time6h(tslongCount + 4))
!                     time6h(1 : tslongCount) = rawNCData      ! Fetch the existing data
!                     time6h(tslongCount + 1 : tslongCount + 4) = rawNCData(tslongCount - 3 : tslongCount) + 1
!                     count = count + 48                      ! Update count
!                 else
!                     time6h = rawNCData
!                 endif
!             endif
!             deallocate(rawNCData)
!             call confirm
!
!             ! Compute range
!             timeStart = time6h(5)
!             timeEnd = time6h(days * 4 + 5)
!
!             ! Expand 6 hourly to half hourly
!             write(message, "(a,i8,a)") " Expanding", count, " time values"; call notify(message)
!             allocate(time(count))
!             time = 0
!             do i = 1, count
!                 day = floor(time6h((i - 1) / 12 + 1))
!                 moment = real(mod(i - 1, 48)) / 48
!                 time(i) = day + moment
!                 !print*, i, day, moment, time(i)
!             enddo
!
!             deallocate(time6h)
!             call confirm
!         end subroutine loadTime
! !
!         !-----------------------------------------------------------------------------------------------------------------------------------------------------
!         ! Load lon/lat dimensions; don't use this function for loading the time dimension
!         subroutine loadDimension(label)
!             character(*), intent(in)    :: label
!             character(200)              :: message
!             integer                     :: count
!             real, allocatable           :: dimension(:)
!
!             count = ncGetDimLen(inputId, label)
!
!             write(message, "(a,i10,a,a,a)") " Loading", count, " ", label, " values from NetCDF file"
!             call notify(message)
!
!             allocate(dimension(count))
!             dimension = ncGetDimValues(inputId, label, count = [count])
!
!             if (label == 'lon') lon = dimension
!             if (label == 'lat') lat = dimension
!
!             call confirm
!         end subroutine loadDimension
!
!         !-----------------------------------------------------------------------------------------------------------------------------------------------------
!         ! Load a variable
!         function loadVariable(label)
!             character(len=2), intent(in)                :: label
!             character(200)                              :: message
!             real, allocatable                           :: loadVariable(:), inputData(:), rawNCData(:) !! inputData is 6h
!             integer                                     :: i, varId
!
!             write(message, "(a,i10,a,a,a)") " Loading", tslongCount, " ", label, " values from NetCDF file"; call notify(message)
!
!             allocate(rawNCData(tslongCount))
!
!             ! Get the existing data
!             rawNCData = ncGet1DVar(inputId, label, start = [lonIndex, latIndex, tslongStart], count = [1, 1, tslongCount])
!
!             ! Fill in missing days
!             if (includesFirstDay) then
!                 if (includesLastDay) then                       ! Need to fill in both the first and last days
!                     allocate(inputData(tslongCount + 8))
!                     inputData(1:4) = rawNCData(1:4) - 1
!                     inputData(5:tslongCount + 4) = rawNCData           ! Fetch the existing data
!                     inputData(tslongCount + 5 : tslongCount + 8) = rawNCData(tslongCount - 3 : tslongCount) + 1
!                 else                                            ! Need to fill in just the first day
!                     allocate(inputData(tslongCount + 4))
!                     inputData(5:tslongCount + 4) = rawNCData           ! Fetch the existing data
!                     inputData(1:4) = rawNCData(1:4) - 1
!                 endif
!             else
!                 if (includesLastDay) then                       ! Need to fill in just the last day
!                     allocate(inputData(tslongCount + 4))
!                     inputData(1 : tslongCount) = rawNCData             ! Fetch the existing data
!                     inputData(tslongCount + 1 : tslongCount + 4) = rawNCData(tslongCount - 3 : tslongCount) + 1
!                 else
!                     inputData = rawNCData
!                 endif
!             endif
!
!             call confirm
!
!             ! Expand 6 hourly to half hourly
!             write(message, "(a,a,a,i10,a)") " Expanding ", label, " to ", count, " values"; call notify(message)
!             allocate(loadVariable(count))
!             loadVariable = 0
!             do i = 1, size(inputData)
!                 loadVariable((i-1) * 12 + 1) = inputData(i)
!             enddo
!             deallocate(inputData, rawNCData)
!
!             call confirm
!         end function loadVariable
!     end subroutine loadTimespan
!
!     !-----------------------------------------------------------------------------------------------------------------------------------------------------
!     ! Performs interpolation and redistribution of data
!     subroutine processVariables
!         call stepInterpolation(longWave)
!         call linearInterpolation(temperature)
!         call linearInterpolation(humidity)
!         call linearInterpolation(wind)
!         call linearInterpolation(pressure)
!         call randomDistribution(precipitation)
!         call diurnalDistribution(shortWave, lat(latIndex), dayIndices)
!
!     contains
!         !-----------------------------------------------------------------------------------------------------------------------------------------------------
!         ! Step interpolation
!         subroutine stepInterpolation(var)
!             real, intent(inout)                         :: var(:)
!             character(200)                              :: message
!             integer                                     :: i, j, current, count
!
!             count = size(var)
!             write(message, "(a,i8,a)") " Performing step interpolation on", count, " values"; call notify(message)
!
!             do i = 1, count/12
!                 do j = 2, 12
!                     var((i - 1) * 12 + j) = var((i - 1) * 12 + 1)
!                 enddo
!             enddo
!
!             call confirm
!         end subroutine stepInterpolation
!
!         !-----------------------------------------------------------------------------------------------------------------------------------------------------
!         ! Linear interpolation
!         subroutine linearInterpolation(var)
!             real, intent(inout)                         :: var(:)
!             character(200)                              :: message
!             integer                                     :: i, j, shortIntervals, start, end, count
!             real                                        :: t
!
!             count = size(var)
!             write(message, "(a,i8,a)") " Performing linear interpolation on", count, " values"; call notify(message)
!
!             ! TODO: Explain this thoroughly
!             shortIntervals = count / 12
!             do i = 1, shortIntervals - 1
!                 var(i * 12) = var(i * 12 + 1)
!             enddo
!             var(count) = var(count - 12)
!             ! up to here!!!
!             do i = 1, shortIntervals
!                 do j = 0, 11
!                     start = (i - 1) * 12 + 1
!                     end = i * 12
!                     t = real(j) / 12
!                     var(start + j) = var(start) + (var(end) - var(start)) * t
!                 enddo
!             enddo
!
!             call confirm
!         end subroutine linearInterpolation
!
!         !-----------------------------------------------------------------------------------------------------------------------------------------------------
!         ! Random distribution
!         subroutine randomDistribution(var)
!             real, intent(inout)                         :: var(:)
!             character(200)                              :: message
!             integer                                     :: i, j, start, end, count
!             real                                        :: t, random(12)
!
!             count = size(var)
!             write(message, "(a,i8,a)") " Performing random interpolation on", count, " values"; call notify(message)
!
!             call random_number(random)
!             random = random / sum(random)
!             do i = 1, count/12
!                 start = (i - 1) * 12 + 1
!                 end = i * 12
!                 var(start:end) = random * var(start)
!             enddo
!
!             call confirm
!         end subroutine randomDistribution
!
!         !-----------------------------------------------------------------------------------------------------------------------------------------------------
!         ! Diurnal distribution over the entire timespan
!         subroutine diurnalDistribution(shortWave, latitude, daysOfYear)
!             real, intent(inout)                         :: shortWave(:), latitude
!             integer, intent(in)                         :: daysOfYear(:)
!             character(200)                              :: message
!             real, allocatable                           :: swDiurnalDistributed(:)
!             integer                                     :: i, start, end, count
!
!             count = size(shortWave)
!             write(message, "(a,i8,a)") " Initializing diurnal distribution on", count, " values"; call notify(message)
!
!             allocate(swDiurnalDistributed(count))
!             swDiurnalDistributed = 0
!             do i = 1, size(shortWave) / 48      ! Once, every day
!                 start = (i - 1) * 48 + 1
!                 end = i * 48
!                 call diurnalDistributionOneDay(shortWave(start:end), daysOfYear(i), latitude, swDiurnalDistributed(start:end))
!             enddo
!             shortWave = swDiurnalDistributed
!             deallocate(swDiurnalDistributed, dayIndices)
!
!             call confirm
!         end subroutine diurnalDistribution
!
!         !-----------------------------------------------------------------------------------------------------------------------------------------------------
!         ! Diurnal distribution over just one day
!         subroutine diurnalDistributionOneDay(swDayValues, day, latitude, swDiurnalDistributed)
!             real,       intent(in)  :: swDayValues(:), latitude
!             integer,    intent(in)  :: day
!             real                    :: swDiurnalDistributed(:)
!             integer                 :: i, time, rawSize, quantity, counter, dayLength, daylightCount
!             real                    :: zenithNoon, swMean, latRad
!             integer,    parameter   :: dayMinutes = 1440, timestep = 30, count = 48
!             real                    :: daylightIndices(count), zenithCos(count), timesteps(count), zenithAngles(count)
!
!             latRad = latitude * pi / 180.00
!             swMean = sum(swDayValues(1:count)) / 4          ! Finds the mean of this day's values
!
!             do i = 1, count                                    ! Creates the time steps in 30 minutes steps
!                 timesteps(i) = (i - 1) * 30
!             enddo
!
!             ! Find the zenith angles at every specified time step
!             call getZenithAngles(timesteps, day, latRad, zenithAngles)
!
!             ! Find the zenith angle at noon
!             zenithNoon = zenithAngles(25)
!
!             ! Find day length
!             daylightCount = 0
!             do i = 1, count
!                 zenithCos(i) = cos(zenithAngles(i))
!                 if (zenithCos(i) >= 0) then
!                     daylightCount = daylightCount + 1
!                     daylightIndices(i) = i
!                 else
!                     daylightIndices(i) = -1
!                 endif
!             enddo
!             dayLength = daylightCount * timestep
!
!             call applyDiurnalDistribution(zenithAngles, daylightIndices, zenithNoon, swMean, timestep, swDiurnalDistributed)
!         end subroutine diurnalDistributionOneDay
!
!             !-----------------------------------------------------------------------------------------------------------------------------------------------------
!         ! Compute zenith angles
!         subroutine getZenithAngles(timesteps, day, latitude, zenithAngles)
!             integer,    intent(in)          :: day
!             real,       intent(in)          :: timesteps(:), latitude
!             real,       intent(out)         :: zenithAngles(:)
!             integer                         :: i, n(4), count
!             real                            :: radHourAngle, degHourAngle, psi, dec, a(4), b(4)
!             integer,    parameter           :: minutesInHour = 60, daysInYear = 365
!             a = (/0.006918, -0.399912, -0.006758, -0.002697/)
!             b = (/0.0, 0.070257, 0.000907, 0.001480/)
!             n = (/0, 1, 2, 3/)
!             count = size(timesteps)
!             !!!! Does it matter if day > 365?!?! The original code was day = mod(day, 365) or some such thing
!             psi = 2 * pi * (day - 1) / daysInYear
!             dec = sum((a * cos(n * psi)) + (b * sin(n * psi)))
!             !FIND THE HOUR ANGLE, CONVERT IT TO RADIANS THEN FIND THE ZENITH ANGLE(S).
!             do i = 1, count
!                 degHourAngle = 15 * (12 - (timesteps(i) / minutesInHour))
!                 radHourAngle = (degHourAngle / 360) * 2 * pi
!                 zenithAngles(i) = acos((sin(latitude) * sin(dec)) + &
!                 (cos(latitude) * cos(dec) * cos(radHourAngle)))
!             enddo
!         end subroutine getZenithAngles
!
!         !-----------------------------------------------------------------------------------------------------------------------------------------------------
!         ! Apply the diurnal distribution on one day
!         subroutine applyDiurnalDistribution(zenithAngles, daylightIndices, zenithNoon, &
!         swMean, timestep, swDiurnalDistributed)
!             integer,    intent(in)          :: timestep
!             real,       intent(in)          :: zenithAngles(:), daylightIndices(:), zenithNoon, swMean
!             real,       intent(out)         :: swDiurnalDistributed(:)
!             real,       allocatable         :: diurnalDistrib(:)
!             integer                         :: i, count, daylightIndex
!             real                            :: sum, correction, diurnalMean
!
!             count = size(zenithAngles)
!             allocate(diurnalDistrib(count))
!
!             sum = 0.0
!             do i = 1, count
!                 daylightIndex = daylightIndices(i)              ! CHECK WITH DAYLIGHT_INDICES TO SEE IF IT THERE IS DAYLIGHT AT THE TIME,
!                 if (daylightIndex > 0) then                     ! IF THERE IS FIND THE VALUE FOR S
!                     diurnalDistrib(i) = swMean * pi / 2 * &
!                     cos((zenithAngles(daylightIndex) - zenithNoon) / (pi / 2 - zenithNoon) * pi / 2)
!                 else                                            ! ELSE SET IT TO 0
!                     diurnalDistrib(i) = 0
!                 endif
!                 sum = sum + diurnalDistrib(i)
!             enddo
!
!             diurnalMean = sum / count
!             if (diurnalMean /= 0) then
!                 correction = swMean / diurnalMean
!             else
!                 correction = 0
!             endif
!             swDiurnalDistributed = diurnalDistrib * correction / 21600 ! Joe needs to explain this one
!
!             deallocate(diurnalDistrib)
!         end subroutine applyDiurnalDistribution
!     end subroutine processVariables
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
!     !-----------------------------------------------------------------------------------------------------------------------------------------------------
!     ! As the initial processing time range is 2 days longer, the desired timespan must be extracted
!     subroutine extractValidTime
!         integer                                     :: startIndex, stopIndex
!         character(200)                              :: message
!
!         write(message, "(a)") " Extracting the desired timespan"; call notify(message)
!
!         ! Find the index of the first day
!         startIndex = findTimestep(timeStart)
!         if (startIndex < 0) stop(" Uh-oh, something sure is wrong with the startIndex...")
!
!         ! Find the index of the last day
!         stopIndex = findTimestep(timeEnd)
!         if (stopIndex < 0) stop(" Uh-oh, something sure is wrong with the stopIndex...")
!
!         ! The last value should be the one before the boundary
!         stopIndex = stopIndex - 1
!
!         ! Extract the time dimension
!         time = time(startIndex : stopIndex)
!
!         ! Extract all the variables
!         longWave = longWave(startIndex : stopIndex)
!         temperature = temperature(startIndex : stopIndex)
!         humidity = humidity(startIndex : stopIndex)
!         wind = wind(startIndex : stopIndex)
!         pressure = pressure(startIndex : stopIndex)
!         precipitation = precipitation(startIndex : stopIndex)
!         shortWave = shortWave(startIndex : stopIndex)
!
!         ! Rectify the new interval count and such
!         count = stopIndex - startIndex + 1
!
!         call confirm
!     end subroutine extractValidTime
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
!     contains
!         !-----------------------------------------------------------------------------------------------------------------------------------------------------
!         ! As real numbers are not precise, this is a way to compare two reals
!         logical function closeEnough(num1, num2)
!             real, intent(in)    :: num1, num2
!             real, parameter     :: error = 0.0001
!             if (abs(num1 - num2) < error) then
!                 closeEnough = .true.
!             else
!                 closeEnough = .false.
!             endif
!         end function closeEnough
!
!     end function findTimestep
!
!     !-----------------------------------------------------------------------------------------------------------------------------------------------------
!     ! Export data to output NetCDF file
!     subroutine exportData
!         character(200)              :: message
!
!         write(message, "(a)") " Starting to export data"; call notify(message)
!
!         !        call ncPut1DVar(fileId, 'sw', shortWave, start = [lonIndex, latIndex, 1], count = [1, 1, size(shortWave)])
!         !        call ncPut1DVar(fileId, 'lw', longWave, start = [lonIndex, latIndex, 1], count = [1, 1, size(longWave)])
!         !        call ncPut1DVar(fileId, 'pr', precipitation, start = [lonIndex, latIndex, 1], count = [1, 1, size(precipitation)])
!         !        call ncPut1DVar(fileId, 'ta', temperature, start = [lonIndex, latIndex, 1], count = [1, 1, size(temperature)])
!         !        call ncPut1DVar(fileId, 'qa', humidity, start = [lonIndex, latIndex, 1], count = [1, 1, size(humidity)])
!         !        call ncPut1DVar(fileId, 'wi', wind, start = [lonIndex, latIndex, 1], count = [1, 1, size(wind)])
!         !        call ncPut1DVar(fileId, 'ap', pressure, start = [lonIndex, latIndex, 1], count = [1, 1, size(pressure)])
!
!         call exportVariable('sw', shortWave)
!         call exportVariable('lw', longWave)
!         call exportVariable('pr', precipitation)
!         call exportVariable('ta', temperature)
!         call exportVariable('qa', humidity)
!         call exportVariable('wi', wind)
!         call exportVariable('ap', pressure)
!
!         call confirm
!     end subroutine exportData
!
!     !-----------------------------------------------------------------------------------------------------------------------------------------------------
!     subroutine exportVariable(label, variable)
!         character(*), intent(in)    :: label
!         real, intent(inout)         :: variable(:)
!         integer                     :: fileId
!         character(200)              :: message, filename, units, gridType, title, name
!         integer                     :: lonDimId, latDimId, timeDimId, varId, err
!
!         ! Make the filename
!         filename = 'metVar_' // label // '.nc'
!
!         ! If the file doesn't already exist, then initialize the file
!         if (.not.fileExists(filename)) then
!             write(message, "(a,a)") ' Creating ', trim(filename); call notify(message)
!
!             ! Create the variable file
!             fileId = ncCreate(filename, NF90_CLOBBER)
!
!             ! Populate the variable with attribute content
!             select case(label)
!             case('sw')
!                 units = 'J/m2'
!                 gridType = 'gaussian'
!                 title = 'Incoming short wave radiation'
!                 name = 'short wave'
!             case('ap')
!                 units = 'Pa'
!                 gridType = 'gaussian'
!                 title = 'Atmospheric pressure'
!                 name = 'athmospheric pressure'
!             case('pr')
!                 units = 'mm/hh'
!                 gridType = 'gaussian'
!                 title = 'Total precipitation'
!                 name = 'precipitation'
!             case('qa')
!                 units = 'g/g'
!                 gridType = 'gaussian'
!                 title = 'Air specific humidity'
!                 name = 'humidity'
!             case('ta')
!                 units = 'K'
!                 gridType = 'gaussian'
!                 title = 'Temperature'
!                 name = 'temperature'
!             case('wi')
!                 units = 'm/s'
!                 gridType = 'gaussian'
!                 title = 'Wind'
!                 name = 'wind'
!             case('lw')
!                 units = ''
!                 gridType = 'gaussian'
!                 title = 'Incoming long wave radiation'
!                 name = 'long wave'
!             case default
!                 stop('Unrecognized label')
!             end select
!
!             ! Define file attributes
!             call ncPutAtt(fileId, nf90_global, 'title', charValue = title)
!             call ncPutAtt(fileId, nf90_global, 'units', charValue = units)
!             call ncPutAtt(fileId, nf90_global, 'grid_type', charValue = gridType)
!             call ncPutAtt(fileId, nf90_global, '_FillValue', realValue = fillValue)
!
!             ! Define the longitude dimension
!             lonDimId = ncDefDim(fileId, 'lon', size(lon))
!             varid = ncDefVar(fileId, 'lon', nf90_double, [lonDimId])
!             call ncPutAtt(fileId, varId, 'standard_name', charValue = 'Longitude')
!             call ncPutAtt(fileId, varId, 'units', charValue = 'degrees_east')
!             call ncPutAtt(fileId, varId, 'axis', charValue = 'X')
!
!             call ncEndDef(fileId)
!
!             call ncPutDimValues(fileId, 'lon', realValues = lon, count = [size(lon)])
!
!             call ncRedef(fileId)
!
!             ! Define the latitude dimension
!             latDimId = ncDefDim(fileId, 'lat', size(lat))
!             varid = ncDefVar(fileId, 'lat', nf90_double, [latDimId])
!             call ncPutAtt(fileId, varId, 'standard_name', charValue = 'Latitude')
!             call ncPutAtt(fileId, varId, 'units', charValue = 'degrees_north')
!             call ncPutAtt(fileId, varId, 'axis', charValue = 'Y')
!
!             call ncEndDef(fileId)
!
!             call ncPutDimValues(fileId, 'lat', lat, count = [size(lat)])
!
!             call ncRedef(fileId)
!
!             ! Define the time dimension
!             timeDimId = ncDefDim(fileId, 'time', size(time))
!             varid = ncDefVar(fileId, 'time', nf90_double, [timeDimId])
!             call ncPutAtt(fileId, varId, 'standard_name', charValue = 'time')
!             call ncPutAtt(fileId, varId, 'units', charValue = 'day as YYYYMMDD.FFFF')
!             call ncPutAtt(fileId, varId, 'calendar', charValue = 'proleptic_gregorian')
!             call ncEndDef(fileId)
!             call ncPutDimValues(fileId, 'time', time, count = [size(time)])
!
!             call ncRedef(fileId)
!
!             ! Define the variable
!             varid = ncDefVar(fileId, label, nf90_double, [lonDimId, latDimId, timeDimId])
!             call ncPutAtt(fileId, varId, '_FillValue', realValue = fillValue)
!             call ncEndDef(fileId)
!
!             call confirm
!         else
!             fileId = ncOpen(filename, nf90_write)
!         endif
!
!         write(message, "(a,a,a)") ' Exporting the ', trim(label), ' data'; call notify(message)
!
!         ! Put in data
!         call ncPut1DVar(fileId, label, variable, start = [lonIndex, latIndex, 1], count = [1, 1, size(variable)])
!
!         call ncClose(fileId)
!
!         call confirm
!     end subroutine exportVariable
!
!     !-----------------------------------------------------------------------------------------------------------------------------------------------------
!     ! Close the NetCDF input file
!     subroutine closeInputFiles
!         character(200)                  :: message
!
!         write(message, "(a)") " Closing the input NetCDF file"; call notify(message)
!
!         call ncClose(inputId)
!
!         call confirm
!     end subroutine closeInputFiles
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
