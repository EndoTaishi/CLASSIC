!>Central module for all netcdf output file operations
module outputManager

    use ctem_statevars, only : c_switch

    implicit none

    public :: generateOutputFiles
    private :: generateNetCDFFile
    private :: addVariable
    private :: validGroup
    private :: validTime
    private :: generateFilename
    private :: getDescriptor
    public :: getIdByKey
    public  :: createNetCDF
    private :: determineTime
    private :: identityVector
    public :: writeOutput1D
    public  :: closeNCFiles
    public :: checkForTime


    type simulationDomain
        real, dimension(:), allocatable     :: lonLandCell, latLandCell     !< Long/Lat values of only the land cells in our model domain
        integer, dimension(:), allocatable  :: lonLandIndex, latLandIndex   !< Indexes of only the land cells in our model domain for our resolution
        real, dimension(:), allocatable     :: allLonValues, allLatValues   !< All long/Lat values in our model domain (including ocean/non-land)
        integer, dimension(:), allocatable  :: lonLocalIndex,latLocalIndex  !< The index for only the region that is being simulated
        real, dimension(:), allocatable     :: latUnique,lonUnique          !< The index for only the region that is being simulated with each value only once
        integer                             :: LandCellCount    !> number of land cells that the model will run over
        real, dimension(4) :: domainBounds                      !> Corners of the domain to be simulated (netcdfs)
        integer :: srtx                                         !> starting index for this simulation for longitudes
        integer :: srty                                         !> starting index for this simulation for latitudes
        integer :: cntx                                         !> number of grid cells for this simulation in the longitude direction
        integer :: cnty                                         !> number of grid cells for this simulation in the latitude direction
    end type

    type(simulationDomain) :: myDomain

    !> This data structure is used to set up the output netcdf files.
    type outputDescriptor
        character(80)   :: group                = ''
        character(30)   :: shortName            = ''
        character(30)   :: standardName         = ''
        character(400)  :: longName             = ''        !< Long name of the variable
        character(30)   :: units                = ''        !< Units of the variable
        character(30)   :: timeFreq             = ''        !< Time frequency of variable: half-hourly, daily, monthly, annually
        logical         :: includeBareGround    = .false.   !< If true then expand the PFT list for a final position that is the bare ground.
    end type

    type(outputDescriptor), allocatable     :: outputDescriptors(:)

    !> Contains characteristic information about the output netcdf files and is used in their creation
    type variant
        character(80)   :: nameInCode           = ''
        character(80)   :: timeFrequency        = ''
        character(80)   :: outputForm           = ''
        character(80)   :: shortName            = ''
    end type

    type(variant), allocatable              :: variants(:)

    !> Stores the information needed to retrieve the output netcdf files to write to them
    type netcdfVar
        integer         :: ncid
        character(30)   :: key
        character(350)   :: filename
    end type

    integer, parameter  :: maxncVariableNumber = 300        !< Maximum number of netcdf output files to make (can be adjusted)
    type(netcdfVar)     :: netcdfVars(maxncVariableNumber)

    integer :: variableCount = 0, descriptorCount = 0, variantCount = 0 !< Initializes the counter variables

    integer         :: refyr = 1850                     !< Time reference for netcdf output files
    character(30)   :: timestart = "days since 1850-01-01 00:00" !< Time reference for netcdf output files
    real   :: fill_value = 1.E38             !< Default fill value for missing values in netcdf output files
    real, dimension(:), allocatable :: timeVect  !< Array of the timesteps in days since refyr for this model run and output file

contains

    !---------------------------------------------------------------------------------------
    !>\ingroup output_generateOutputFiles
    !>@{
    !> Runs through all possible variants and calls for the generation of the required output files.
    subroutine generateOutputFiles

        implicit none

        integer             :: i

        do i = 1, variantCount
            call generateNetCDFFile(variants(i)%nameInCode, variants(i)%timeFrequency,&
                variants(i)%outputForm, variants(i)%shortName)
        enddo

        return

    end subroutine generateOutputFiles

    !<@}
    !---------------------------------------------------------------------------------------

    !>\ingroup output_generateNetCDFFile
    !>@{
    !> Generates a new (netcdf) variable
    subroutine generateNetCDFFile(nameInCode, timeFreq, outputForm, descriptorLabel)

        implicit none

        character(*), intent(in)                :: nameInCode, timeFreq, outputForm, descriptorLabel
        type(outputDescriptor)                  :: descriptor
        character(350)                           :: filename = ''
        logical                                 :: isTimeValid, isGroupValid, fileCreatedOk
        integer                                 :: id

        ! Get variable descriptor
        descriptor = getDescriptor(descriptorLabel)

        ! If the requested timeFreq matches the project config, all's good
        isTimeValid = validTime(timeFreq, descriptor)

        ! If the group property of descriptor matches the project config, all's good
        isGroupValid = validGroup(descriptor,outputForm)

        ! If the project config and variable descriptor match the request, process current variable
        if (isTimeValid .and. isGroupValid) then

            ! Generate the filename
            filename = generateFilename(outputForm, descriptor)

            ! Allocate a new variable (ncid, filename, key etc.)
            id = addVariable(nameInCode, filename)

            ! Make the netcdf file for the new variable (mostly definitions)
            call createNetCDF(filename, id, outputForm, descriptor, timeFreq)

            ! Now make sure the file was properly created
            fileCreatedOk = checkFileExists(filename)

            if (.not. fileCreatedOk) then
                print*,'Failed to create',filename
                print*,'Aborting'
                stop
            end if
        endif

    end subroutine generateNetCDFFile

    !<@}
    !---------------------------------------------------------------------------------------

    !>\ingroup output_checkFileExists
    !>@{
    !> Checks if file exists
    logical function checkFileExists(filename)

        character(*), intent(in) :: filename

        inquire(file=filename, exist=checkFileExists)

    end function

    !<@}
    !---------------------------------------------------------------------------------------

    !>\ingroup output_addVariable
    !>@{
    !> Adds the new variable to the list of variables (see the type "netcdfVars")
    integer function addVariable(key, filename)
        use fileIOModule
        implicit none
        character(*), intent(in)    :: key, filename
        integer                     :: ncid

        ncid = ncCreate(fileName, cmode=NF90_CLOBBER)

        variableCount = variableCount + 1
        netcdfVars(variableCount)%ncid = ncid
        netcdfVars(variableCount)%key = key
        netcdfVars(variableCount)%filename = filename
        addVariable = variableCount

    end function addVariable

    !<@}
    !---------------------------------------------------------------------------------------

    !>\ingroup output_validGroup
    !>@{
    !> Determines if the current variable matches the project configuration
    logical function validGroup(descriptor,outputForm)

        implicit none

        character(*), intent(in)              :: outputForm
        type(outputDescriptor), intent(in)    :: descriptor

        ! First weed out the per pft, per tile options if they should not be used.
        if (.not. c_switch%doperpftoutput .and. trim(outputForm) == "pft") then
            validGroup = .false.
            return
        end if

        if (.not. c_switch%dopertileoutput .and. trim(outputForm) == "tile") then
            validGroup = .false.
            return
        end if

        ! Now check over the remaining options
        if (trim(descriptor%group) == "class") then !CLASS outputs always are valid
            validGroup = .true.
        elseif (c_switch%ctem_on .and. trim(descriptor%group) == "ctem") then
            validGroup = .true.
        elseif (c_switch%ctem_on) then ! check the CTEM sub-switches
            if (c_switch%dofire .and. trim(descriptor%group) == "fire") then
                validGroup = .true.
            elseif (c_switch%lnduseon .and. trim(descriptor%group) == "land") then
                validGroup = .true.
            elseif (c_switch%dowetlands .and. trim(descriptor%group) == "methane") then
                validGroup = .true.
            elseif (c_switch%PFTCompetition .and. trim(descriptor%group) == "PFTCompetition") then
                validGroup = .true.
            else
                validGroup = .false.
            end if
        else
            validGroup = .false.
        endif
    end function validGroup

    !<@}
    !---------------------------------------------------------------------------------------

    !>\ingroup output_validTime
    !>@{
    !> Determines wether the current variable matches the project configuration
    logical function validTime(timeFreq, descriptor)

        implicit none

        !type(projectConfiguration), intent(in)  :: config
        type(outputDescriptor), intent(inout) :: descriptor
        character(*), intent(in)                :: timeFreq
        logical                                 :: valid

        valid = .true.
        if (timeFreq == 'annually') then
            descriptor%timeFreq = timeFreq
        elseif (c_switch%domonthoutput .and. timeFreq == 'monthly') then
            descriptor%timeFreq = timeFreq
        elseif (c_switch%dodayoutput .and. timeFreq == 'daily') then
            descriptor%timeFreq = timeFreq
        elseif (c_switch%dohhoutput .and. timeFreq == 'halfhourly') then
            descriptor%timeFreq = timeFreq
        else
            valid = .false.
        endif
        validTime = valid
    end function validTime

    !<@}
    !---------------------------------------------------------------------------------------

    !>\ingroup output_generateFilename
    !>@{
    !> Generates the filename for the current variable
    character(350) function generateFilename(outputForm, descriptor)

        implicit none

        !type(projectConfiguration), intent(in)  :: config
        type(outputDescriptor), intent(in)    :: descriptor
        character(*), intent(in)                :: outputForm
        character(80)                           :: suffix = ''

        select case(trim(outputForm))
            case("pft")
                suffix = '_perpft'
            case ("tile")
                suffix = '_pertil'
            case default
                suffix = ''
        end select
        generateFilename = trim(c_switch%output_directory) // '/' // &
        trim(descriptor%shortName) // '_' // &
        trim(descriptor%timeFreq) // trim(suffix) // '.nc'
    end function generateFilename

    !<@}
    !---------------------------------------------------------------------------------------

    !>\ingroup output_getDescriptor
    !>@{
    !> Retrieve a variable descriptor based on a given key (e.g. shortName)
    type (outputDescriptor) function getDescriptor(key)

        implicit none

        character(len=*), intent(in)       :: key
        integer                            :: i

        do i = 1, descriptorCount
            if (outputDescriptors(i)%shortName == key) then
                getDescriptor = outputDescriptors(i)
                return
            endif
        enddo
        print*, "something went awry with the getDescriptor function"
    end function getDescriptor

    !<@}
    !---------------------------------------------------------------------------------------

    !>\ingroup output_getIdByKey
    !>@{
    !> Find the id of the variable with the following key
    integer function getIdByKey(key)

        implicit none
        character(*), intent(in)   :: key
        integer i
        do i=1, variableCount
            if (netcdfVars(i)%key == key) then
                getIdByKey = i
                return
            endif
        enddo
        getIdByKey = 0
    end function getIdByKey

    !<@}
    !---------------------------------------------------------------------------------------

    !>\ingroup output_createNetCDF
    !>@{
    ! Create the output netcdf files
    subroutine createNetCDF(fileName, id, outputForm, descriptor,timeFreq)

        use fileIOModule
        use ctem_statevars,     only : c_switch
        use ctem_params,        only : ignd,icc,nmos,iccp1

        implicit none

        character(*), intent(in)              :: fileName
        character(*), intent(in)              :: outputForm
        character(*), intent(in)              :: timeFreq
        type(outputDescriptor), intent(in)    :: descriptor
        integer, intent(in)                   :: id

        character(8)  :: today
        character(10) :: now
        integer                     :: ncid, varid, suffix,i,timeLength
        integer                     :: DimId,lonDimId,latDimId,tileDimId,pftDimId,layerDimId,timeDimId
        real, dimension(2)          :: xrange, yrange
        integer, dimension(:), allocatable :: intArray

        character(:), pointer :: Comment   !< Comment about the run that will be written to the output netcdfs
        logical, pointer :: leap           !< set to true if all/some leap years in the .MET file have data for 366 days
                                           !< also accounts for leap years in .MET when cycling over meteorology (metLoop > 1)

        ! Point pointers
        leap => c_switch%leap
        Comment => c_switch%Comment

        ncid=netcdfVars(variableCount)%ncid

        ! ---------
        call ncPutAtt(ncid, nf90_global, 'title',charvalues='CLASSIC output file')
        call date_and_time(today,now)

        call ncPutAtt(ncid,nf90_global,'timestamp',charvalues=today//' '//now(1:4))
        call ncPutAtt(ncid,nf90_global,'Conventions',charvalues='COARDS')
        !call ncPutAtt(ncid,nf90_global,'node_offset',intvalues=1)

        !----1 - Longitude

        lonDimId = ncDefDim(ncid,'lon',myDomain%cntx)
        varid = ncDefVar(ncid,'lon',nf90_double,[lonDimId])
        call ncPutAtt(ncid,varid,'long_name',charvalues='longitude')
        call ncPutAtt(ncid,varid,'units',charvalues='degrees_east')
        !call ncPutAtt(ncid,varid,'actual_range',xrange) #FLAG need to find the xrange from all_lon.
        !call ncPutAtt(ncid,varid,'_Storage',charvalues="contiguous")


        !----2 - Latitude
        latDimId = ncDefDim(ncid,'lat',myDomain%cnty)
        varid = ncDefVar(ncid,'lat',nf90_double,[latDimId])
        call ncPutAtt(ncid,varid,'long_name',charvalues='latitude')
        call ncPutAtt(ncid,varid,'units',charvalues='degrees_north')
        !call ncPutAtt(ncid,varid,'actual_range',yrange) #FLAG need to find the xrange from all_lon.
        !call ncPutAtt(ncid,varid,'_Storage',charvalues="contiguous")

        select case(trim(outputForm))

            case ("tile")       ! Per tile outputs

                tileDimId = ncDefDim(ncid,'tile',nmos)
                varid = ncDefVar(ncid,'tile',nf90_short,[tileDimId])
                call ncPutAtt(ncid,varid,'long_name',charvalues='tile')
                call ncPutAtt(ncid,varid,'units',charvalues='tile number')
                !call ncPutAtt(ncid,varid,'_Storage',charvalues="contiguous")
                !call ncPutAtt(ncid,varid,'_Endianness',charvalues="little")

            case("pft")         ! Per PFT outputs

                if (descriptor%includeBareGround) then
                    pftDimId = ncDefDim(ncid,'pft',iccp1)
                else
                    pftDimId = ncDefDim(ncid,'pft',icc)
                end if

                varid = ncDefVar(ncid,'pft',nf90_short,[pftDimId])
                call ncPutAtt(ncid,varid,'long_name',charvalues='Plant Functional Type')
                call ncPutAtt(ncid,varid,'units',charvalues='PFT')
                !call ncPutAtt(ncid,varid,'_Storage',charvalues="contiguous")
                !call ncPutAtt(ncid,varid,'_Endianness',charvalues="little")

            case ("layer")      ! Per layer outputs

                layerDimId = ncDefDim(ncid,'layer',ignd)
                varid = ncDefVar(ncid,'layer',nf90_short,[layerDimId])
                call ncPutAtt(ncid,varid,'long_name',charvalues='soil layer')
                call ncPutAtt(ncid,varid,'units',charvalues='layer')
                !call ncPutAtt(ncid,varid,'_Storage',charvalues="contiguous")
                !call ncPutAtt(ncid,varid,'_Endianness',charvalues="little")

        end select

        ! Figure out the total run length, make a time vector and add to file.
        call determineTime(timeFreq)

        timeLength = size(timeVect)

        ! Set up the time dimension
        timeDimId = ncDefDim(ncid,'time',timeLength)

        varid = ncDefVar(ncid,'time',nf90_double,[timeDimId])

        call ncPutAtt(ncid,varid,'long_name',charvalues='time')
        call ncPutAtt(ncid,varid,'units',charvalues=trim(timestart))

        if (leap) then
            call ncPutAtt(ncid,varid,'calendar',charvalues="standard")
        else
            call ncPutAtt(ncid,varid,'calendar',charvalues="365_day")
        end if

        !call ncPutAtt(ncid,varid,'_Storage',charvalues="chunked")
        !call ncPutAtt(ncid,varid,'_Chunksizes',intvalues=1)
        !call ncPutAtt(ncid,varid,'_Endianness',charvalues="little")

        call ncEndDef(ncid)

        call ncPutDimValues(ncid, 'time', realValues=timeVect, count=(/timelength/))

        deallocate(timeVect) !needs to be deallocated so the next file can allocate it.

        ! Fill in the dimension variables and define the model output vars

        call ncPutDimValues(ncid, 'lon', realValues=myDomain%lonUnique, count=(/myDomain%cntx/))
        call ncPutDimValues(ncid, 'lat', realValues=myDomain%latUnique, count=(/myDomain%cnty/))
        
        select case(trim(outputForm))
            case ("tile")       ! Per tile outputs
                allocate(intArray(nmos))
                intArray=identityVector(nmos)
                call ncPutDimValues(ncid, 'tile', intValues=intArray, count=(/nmos/))
                call ncReDef(ncid)
                varid = ncDefVar(ncid, trim(descriptor%shortName), nf90_double, [lonDimId,latDimId,tileDimId,timeDimId])

            case("pft")         ! Per PFT outputs

                if (descriptor%includeBareGround) then
                    allocate(intArray(iccp1))
                    intArray=identityVector(iccp1)
                    call ncPutDimValues(ncid, 'pft', intValues=intArray, count=(/iccp1/))
                else
                    allocate(intArray(icc))
                    intArray=identityVector(icc)
                    call ncPutDimValues(ncid, 'pft', intValues=intArray, count=(/icc/))
                end if


                call ncReDef(ncid)

                if (descriptor%includeBareGround) then
                    ! do something for cells that have bare ground
                    suffix = suffix + 1
                else
                    ! do something for cells that don't have bare ground
                    suffix = suffix
                endif
                varid = ncDefVar(ncid, trim(descriptor%shortName), nf90_double, [lonDimId,latDimId,pftDimId,timeDimId])

            case ("layer")      ! Per layer outputs

                allocate(intArray(ignd))
                intArray=identityVector(ignd)
                call ncPutDimValues(ncid, 'layer',intValues=intArray, count=(/ignd/)) 
                call ncReDef(ncid)
                varid = ncDefVar(ncid, trim(descriptor%shortName), nf90_double, [lonDimId,latDimId,layerDimId,timeDimId])

            case default        ! Grid average outputs

                call ncReDef(ncid)
                varid = ncDefVar(ncid, trim(descriptor%shortName), nf90_double, [lonDimId,latDimId,timeDimId])

        end select

        call ncPutAtt(ncid,varid,'long_name',charvalues=descriptor%longName)
        call ncPutAtt(ncid,varid,'units',charvalues=descriptor%units)
        call ncPutAtt(ncid,varid,'_FillValue',realvalues=fill_value)
        !call ncPutAtt(ncid,varid,'missing_value',realvalues=fill_value)
        !call ncPutAtt(ncid,varid,'_Storage',charvalues="chunked")
        !call ncPutAtt(ncid,varid,'_DeflateLevel',intvalues=1)
        call ncPutAtt(ncid,nf90_global,'Comment',c_switch%Comment)

        call ncEndDef(ncid)

    end subroutine createNetCDF

    !<@}
    !---------------------------------------------------------------------------------------
    !>\ingroup output_determineTime
    !>@{
    !> Determine the time vector for this run. This implictly
    !! assumes that leap year meteorological forcing is used for runs with metLoop = 1, otherwise
    !! the timing of the leap years will be off in the output files.
    subroutine determineTime(timeFreq)

        use fileIOModule
        use ctem_statevars,     only : c_switch
        use generalUtils,       only : findLeapYears
        use ctem_params,        only : monthend

        implicit none

        character(*), intent(in)              :: timeFreq
        real, allocatable, dimension(:) :: temptime
        integer :: totsteps, totyrs, i, st, en, j, m, cnt
        integer :: lastDOY, length
        logical :: leapnow

        logical, pointer :: leap           !< set to true if all/some leap years in the .MET file have data for 366 days
                                           !< also accounts for leap years in .MET when cycling over meteorology (metLoop > 1)
        integer, pointer :: metLoop        !< no. of times the meteorological data is to be looped over. this
                                           !< option is useful to equilibrate CTEM's C pools
        integer, pointer :: readMetStartYear !< First year of meteorological forcing to read in from the met file
        integer, pointer :: readMetEndYear   !< Last year of meteorological forcing to read in from the met file
        integer, pointer :: jhhstd  !< day of the year to start writing the half-hourly output
        integer, pointer :: jhhendd !< day of the year to stop writing the half-hourly output
        integer, pointer :: jdstd   !< day of the year to start writing the daily output
        integer, pointer :: jdendd  !< day of the year to stop writing the daily output
        integer, pointer :: jhhsty  !< simulation year (iyear) to start writing the half-hourly output
        integer, pointer :: jhhendy !< simulation year (iyear) to stop writing the half-hourly output
        integer, pointer :: jdsty   !< simulation year (iyear) to start writing the daily output
        integer, pointer :: jdendy  !< simulation year (iyear) to stop writing the daily output
        integer, pointer :: jmosty    !< Year to start writing out the monthly output files. If you want to write monthly outputs right

        leap            => c_switch%leap
        metLoop         => c_switch%metLoop
        readMetStartYear=> c_switch%readMetStartYear
        readMetEndYear  => c_switch%readMetEndYear
        jhhstd          => c_switch%jhhstd
        jhhendd         => c_switch%jhhendd
        jdstd           => c_switch%jdstd
        jdendd          => c_switch%jdendd
        jhhsty          => c_switch%jhhsty
        jhhendy         => c_switch%jhhendy
        jdsty           => c_switch%jdsty
        jdendy          => c_switch%jdendy
        jmosty          => c_switch%jmosty

        lastDOY = 365
        select case(trim(timeFreq))

            case("annually")
                totyrs = (readMetEndYear - readMetStartYear + 1) * metLoop
                totsteps = totyrs
                allocate(timeVect(totsteps))
                do i = 1, totsteps
                    if (leap) call findLeapYears(readMetStartYear + i - 1,leapnow,lastDOY)
                    timeVect(i) = (readMetStartYear + i - 1 - refyr) * lastDOY
                end do
 
            case("monthly")
                ! Monthly may start writing later (after jmosty) so make sure to account for that.
                ! Also if leap years are on, it changes the timestamps
                if (readMetStartYear < jmosty) then
                    totyrs = (readMetEndYear - jmosty + 1) * metLoop
                else
                    totyrs = (readMetEndYear - readMetStartYear + 1) * metLoop
                end if
                totsteps = totyrs*12
                allocate(timeVect(totsteps))
                do i = 1, totyrs
                    ! Find out if this year is a leap year. It adjusts the monthend array.
                    if (leap) call findLeapYears(readMetStartYear + i - 1,leapnow,lastDOY)
                    do m = 1, 12
                        j = ((i - 1) * 12) + m
                        timeVect(j) = (readMetStartYear + i - 1 - refyr) * lastDOY + monthend(m+1) - 1
                    end do
                end do
             
            case("daily")
                ! Daily may start writing later (after jdsty) and end earlier (jdendy) so make sure to account for that.
                ! Also likely doesn't do all days of the year. Lastly if leap years are on, it changes the timestamps
                ! First determine the number of years

                ! Sanity check on jdsty and jdendy
                if (readMetEndYear < jdsty .or. readMetStartYear > jdendy) then
                    print*,'**addTime says: Check your daily output file start and end points, they are outside the range of this run'
                    stop
                end if
                if (readMetStartYear < jdsty) then
                    st = jdsty
                else
                    st = readMetStartYear
                end if
                if (readMetEndYear > jdendy) then
                    en = jdendy
                else
                    en = readMetEndYear
                end if
                totyrs = (en - st + 1) * metLoop
                ! Now determine the total number of timesteps (days) across all years
                totsteps = 0
                allocate(timeVect(0))
                cnt=0
                ! Create the time vector to write to the file
                do i = 1, totyrs
                    if (leap) call findLeapYears(readMetStartYear + i - 1,leapnow,lastDOY)
                    st = max(1, jdstd)
                    en = min(jdendd, lastDOY)
                    totsteps = totsteps + (en - st + 1)
                    allocate(temptime(totsteps))
                    length = size(timeVect)
                    if (i > 1) then
                        temptime(1 : length) = timeVect
                    end if
                    do j = st,en
                        cnt=cnt+1
                        temptime(cnt) = (readMetStartYear + i - 1 - refyr) * lastDOY + j
                end do
                    call move_alloc(temptime,timeVect)
                end do

            case("halfhourly")
                ! Similar to daily in that it may start writing later (after jhhsty) and end earlier (jhhendy) so make sure to account for that.
                ! Also likely doesn't do all days of the year. Lastly if leap years are on, it changes the timestamps

                ! Sanity check on jhhsty and jhhendy
                if (readMetEndYear < jhhsty .or. readMetStartYear > jhhendy) then
                    print*,'**addTime says: Check your half-hourly output file start and end points, they are outside the range of this run'
                    stop
                end if

                if (readMetStartYear < jhhsty) then
                    st = jhhsty
                else
                    st = readMetStartYear
                end if
                if (readMetEndYear > jhhendy) then
                    en = jhhendy
                else
                    en = readMetEndYear
                end if
                totyrs = (en - st + 1) * metLoop

                ! Now determine the total number of timesteps (halfhours) across all years
                totsteps = 0
                allocate(timeVect(0))
                cnt=0
                ! Create the time vector to write to the file
                do i = 1, totyrs
                    if (leap) call findLeapYears(readMetStartYear + i - 1,leapnow,lastDOY)
                    st = max(1, jhhstd)
                    en = min(jhhendd, lastDOY)
                    totsteps = totsteps + (en - st + 1) * 48 ! 48 half hours in a day. !FLAG change this if delt is not half-hour

                    allocate(temptime(totsteps))
                    length = size(timeVect)
                    if (i > 1) then
                        temptime(1 : length) = timeVect
                    end if
                    do j = st,en
                        do m = 1,48
                            cnt=cnt+1
                            temptime(cnt) = (readMetStartYear + i - 1 - refyr) * lastDOY + j + (m - 1) / 48.
                        end do
                    end do
                    call move_alloc(temptime,timeVect)
                end do
            case default
                print*,'addTime says - Unknown timeFreq: ',timeFreq
                stop
        end select

    end subroutine determineTime

    !<@}
    !---------------------------------------------------------------------------------------

    !>\ingroup output_writeOutput1D
    !>@{
    ! Write model outputs to already created netcdf files
    subroutine writeOutput1D(lonLocalIndex,latLocalIndex,key,timeStamp,label,data,specStart)

        use fileIOModule

        implicit none

        integer, intent(in) :: lonLocalIndex,latLocalIndex
        character(*), intent(in) :: key
        integer :: ncid, timeIndex, id, length
        real, dimension(:), intent(in) :: timeStamp
        real, dimension(:), intent(in)  :: data
        character(*), intent(in)        :: label
        integer,  optional :: specStart
        integer :: start
        integer :: posTimeWanted
        real, dimension(:), allocatable :: localData
        real, dimension(1) :: localStamp
        real, allocatable, dimension(:) :: timeWritten
        character(90) :: errmsg

        allocate(localData(size(data)))
        localStamp = timeStamp
        localData = data

        !print*,key,timeStamp,label,lonLocalIndex,latLocalIndex
        id= getIdByKey(key)

        if (id == 0) then
            errmsg = 'writeOutput1D says: Your requested key does not exist (' // trim(key) // ') in netcdfVars.'
            print*,trim(errmsg)
            stop
        end if

        ncid = netcdfVars(id)%ncid

        length = size(data)
        if (present(specStart)) then
            start = specStart
        else
            start = 1
        end if

        ! Check if the time period has already been added to the file
        timeIndex = ncGetDimLen(ncid, "time")

        allocate(timeWritten(timeIndex))
        timeWritten= ncGetDimValues(ncid, "time", count = (/timeIndex/))
        posTimeWanted = checkForTime(timeIndex,timeWritten,localStamp(1))

            if (posTimeWanted == 0) then
            print*,'missing timestep in output file',key,localStamp
                stop
            else
                timeIndex = posTimeWanted
            end if

        if (length > 1) then
           call ncPutVar(ncid, label, localData, start=[lonLocalIndex,latLocalIndex,start,timeIndex], count=[1,1,length,1])
        else
           call ncPutVar(ncid, label, localData, start=[lonLocalIndex,latLocalIndex,timeIndex], count=[1,1,1])
        end if

    end subroutine writeOutput1D

    !<@}
    !---------------------------------------------------------------------------------------

    !>\ingroup output_checkForTime
    !>@{
    ! Find if a time period is already in the timeIndex of the file
    integer function checkForTime(timeIndex,timeWritten,timeStamp)

        implicit none

        real, intent(in)   :: timeStamp
        integer, intent(in) :: timeIndex
        real, dimension(:), intent(in) :: timeWritten
        integer i
        do i=1, timeIndex
            if (timeWritten(i) == timeStamp) then
                checkForTime = i
                return
            endif
        enddo
        checkForTime = 0
    end function checkForTime

    !<@}
    !---------------------------------------------------------------------------------------

    !>\ingroup output_closeNCFiles
    !>@{
    ! Close all output netcdfs or just a select file
    subroutine closeNCFiles(incid)

        use fileIOModule

        implicit none

        integer, optional   :: incid
        integer i

        if (present(incid)) then
            call ncClose(incid)
        else
            do i=1, variableCount
                call ncClose(netcdfVars(i)%ncid)
            enddo
        end if

    end subroutine closeNCFiles
    !<@}
    !---------------------------------------------------------------------------------------

    !>\ingroup output_identityVector
    !>@{
    pure function identityVector(n) result(res)
        integer, allocatable ::res(:)
        integer, intent(in) :: n
        integer             :: i
        allocate(res(n))
        forall (i=1:n)
            res(i) = i
        end forall
    end function identityVector

    !<@}
    !---------------------------------------------------------------------------------------

!>\namespace output

end module outputManager
