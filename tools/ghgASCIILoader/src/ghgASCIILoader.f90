program ghgASCIILoader
    use fileIOModule
    implicit none
    real, allocatable   :: moleFrac(:), time(:)
    integer, allocatable    :: hour(:), minute(:), day(:), year(:)
    character(300)          :: inputFile, ghgName
    integer                 :: timesteps, i, currentYear
    integer, dimension(12)  :: daysInMonth = [ 31,28,31,30,31,30,31,31,30,31,30,31 ]
    real, parameter         :: fillValue = 1.E38
    real                    :: lon, lat

    call processArguments
    open(unit = 10, file = inputFile, form = 'formatted', status = 'old', action = 'read')
    call initializeData
    call loadData
    call exportVariable(trim(ghgName), moleFrac)
    close(unit = 10)
contains
    subroutine processArguments
        if (iargc() .ne. 2) then
            stop('Usage is: ghgASCIILoader [ghg name, e.g. CO2] [input file]')
        endif
        call getarg(1, ghgName)
        call getarg(2, inputFile)
    end subroutine processArguments

    subroutine initializeData
        integer                 :: fileSize

        inquire(file = inputFile, size = fileSize)
        timesteps = (fileSize + 1) / 91; ! FLAG check on this. why 91?
        print*,timesteps
        allocate(year(timesteps))
        allocate(moleFrac(timesteps))
        allocate(time(timesteps))
    end subroutine initializeData

    subroutine loadData()
        integer                 :: i
        character(*), parameter :: format = '(I4, 2X, F5.2)'

        do i = 1, timesteps
            read(unit = 10, fmt = format) year(i), moleFrac(i)
        enddo
    end subroutine loadData

    real function buildTimestep(hour, minute, dayIn, year)
        integer, intent(in) :: hour, minute, dayIn, year
        real                :: time, day
        integer             :: i, month

        day = dayIn
        time = (real(hour) * 60 + real(minute)) / 1440  ! 1440 minutes in a day

        month = 0
        do i = 1, 12
            if (day > daysInMonth(i)) then
                day = day - daysInMonth(i)
            else
                month = i
                exit
            endif
        enddo

        buildTimestep = year * 10000 + month * 100 + day + time
    end function buildTimestep

    real function charToReal(input)
        character(len=*), intent(in)    :: input    !< Char input
        read(input,*) charToReal
    end function charToReal

    subroutine exportVariable(label, variable)
        character(*), intent(in)    :: label
        real, intent(inout)         :: variable(:)
        integer                     :: fileId
        character(200)              :: filename, units, gridType, title, name
        integer                     :: varId, timeDimId, lonDimId, latDimId

        ! Make the filename
        filename = 'metVar_' // label // '.nc'

        ! If the file doesn't already exist, then initialize the file
        if (.not.fileExists(filename)) then
            ! Create the variable file
            fileId = ncCreate(filename, NF90_CLOBBER)

            ! Populate the variable with attribute content
            select case(label)
            case('sw')
                units = 'J/m2'
                gridType = 'gaussian'
                title = 'Incoming short wave radiation'
                name = 'short wave'
            case default
                stop('Unrecognized label')
            end select

            ! Define file attributes
            call ncPutAtt(fileId, nf90_global, 'title', charValues = title)
            call ncPutAtt(fileId, nf90_global, 'units', charValues = units)
            call ncPutAtt(fileId, nf90_global, 'grid_type', charValues = gridType)
            call ncPutAtt(fileId, nf90_global, '_FillValue', realValues = fillValue)

            ! Define the time dimension
            timeDimId = ncDefDim(fileId, 'time', size(time))
            varid = ncDefVar(fileId, 'time', nf90_double, [timeDimId])
            call ncPutAtt(fileId, varId, 'standard_name', charValues = 'time')
            call ncPutAtt(fileId, varId, 'units', charValues = 'day as YYYYMMDD.FFFF')
            call ncPutAtt(fileId, varId, 'calendar', charValues = 'proleptic_gregorian')
            call ncEndDef(fileId)
            call ncPutDimValues(fileId, 'time', time, count = [size(time)])

            call ncRedef(fileId)

            ! Define the variable
            varid = ncDefVar(fileId, label, nf90_double, [lonDimId, latDimId, timeDimId])
            call ncPutAtt(fileId, varId, '_FillValue', realValues = fillValue)
            call ncEndDef(fileId)

        else
            fileId = ncOpen(filename, nf90_write)
        endif

        ! Put in data
        call ncPutVar(fileId, label, realValues = variable, start = [1, 1, 1], count = [1, 1, size(variable)])

        call ncClose(fileId)
    end subroutine exportVariable

    logical function fileExists(filename)
        character(*), intent(in) :: filename
        inquire(file = filename, exist = fileExists)
    end function fileExists

end program
