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

    integer :: metid                               !> netcdf file id for the meteorology file
    integer :: initid                               !> netcdf file id for the model initialization file
    integer :: rsid                                 !> netcdf file id for the model restart file
    integer :: co2id                                !> netcdf file id for the CO2 input file
    integer :: ch4id                                !> netcdf file id for the CH4 input file
    integer :: popid                                !> netcdf file id for the population density input file
    integer :: lghtid                               !> netcdf file id for the lightning density input file
    integer :: lucid                                !> netcdf file id for the land use change input file

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

    type variant
        character(80)   :: nameInCode           = ''
        character(80)   :: timeFrequency        = ''
        character(80)   :: outputForm           = ''
        character(80)   :: shortName            = ''
    end type

    type(variant), allocatable              :: variants(:)

    type netcdfVar
        integer         :: ncid
        character(30)   :: key
        character(80)   :: filename
    end type

    integer, parameter  :: maxncVariableNumber = 300
    type(netcdfVar)     :: netcdfVars(maxncVariableNumber)

    integer :: variableCount = 0, descriptorCount = 0, variantCount = 0

    integer         :: refyr = 1850                     !< Time reference for netcdf output files
    character(30)   :: timestart = "days since 1850-01-01 00:00" !< Time reference for netcdf output files
    real   :: fill_value = 1.E38             !< Default fill value for missing values in netcdf output files

contains

    subroutine generateOutputFiles

        implicit none

        integer             :: i

        do i = 1, variantCount
            call generateNetCDFFile(variants(i)%nameInCode, variants(i)%timeFrequency,&
                variants(i)%outputForm, variants(i)%shortName)
        enddo

        return

    end subroutine generateOutputFiles

    subroutine generateNetCDFFile(nameInCode, timeFreq, outputForm, descriptorLabel)
    ! Generates a new (netcdf) variable

        implicit none

        character(*), intent(in)                :: nameInCode, timeFreq, outputForm, descriptorLabel
        type(outputDescriptor)                  :: descriptor
        character(80)                           :: filename = ''
        logical                                 :: isTimeValid, isGroupValid, fileCreatedOk
        integer                                 :: id

        ! Get variable descriptor
        descriptor = getDescriptor(descriptorLabel)

        ! If the requested timeFreq matches the project config, all's good
        isTimeValid = validTime(timeFreq, descriptor)

        ! If the group property of descriptor matches the project config, all's good
        isGroupValid = validGroup(descriptor)

        ! If the project config and variable descriptor match the request, process current variable
        if (isTimeValid .and. isGroupValid) then

            ! Generate the filename
            filename = generateFilename(outputForm, descriptor)

            ! Allocate a new variable (ncid, filename, key etc.)
            id = addVariable(nameInCode, filename)

            ! Make the netcdf file for the new variable (mostly definitions)
            call createNetCDF(filename, id, outputForm, descriptor)

            ! Now make sure the file was properly created
            fileCreatedOk = checkFileExists(filename)

            if (.not. fileCreatedOk) then
                print*,'Failed to create',filename
                print*,'Aborting'
                stop
            end if
        endif

    end subroutine generateNetCDFFile

    logical function checkFileExists(filename)

        character(*), intent(in) :: filename

        inquire(file=filename, exist=checkFileExists)

    end function

    ! Adds the new variable to the list of variables (see the type "netcdfVars")
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

    ! Determines if the current variable matches the project configuration
    logical function validGroup(descriptor)

        implicit none

        type(outputDescriptor), intent(in)    :: descriptor

        if (trim(descriptor%group) == "class") then !CLASS outputs always are valid
            validGroup = .true.
        elseif (c_switch%ctem_on .and. trim(descriptor%group) == "ctem") then
            validGroup = .true.
        elseif (c_switch%dofire .and. trim(descriptor%group) == "fire") then
            validGroup = .true.
        elseif (c_switch%lnduseon .and. trim(descriptor%group) == "land") then
            validGroup = .true.
        elseif (c_switch%dowetlands .and. trim(descriptor%group) == "methane") then
            validGroup = .true.
        elseif (c_switch%PFTCompetition .and. trim(descriptor%group) == "PFTCompetition") then
            validGroup = .true.
        else
            validGroup = .false.
        endif
    end function validGroup

    ! Determines wether the current variable matches the project configuration
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

    ! Generates the filename for the current variable
    character(80) function generateFilename(outputForm, descriptor)

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

    ! Retrieve a variable descriptor based on a given key (e.g. shortName)
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

    ! Find the id of the variable with the following key
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

    ! Create the output netcdf files
    subroutine createNetCDF(fileName, id, outputForm, descriptor)

        use fileIOModule
        use ctem_statevars,     only : c_switch
        use ctem_params,        only : ignd,icc,nmos,iccp1

        implicit none

        character(*), intent(in)              :: fileName
        character(*), intent(in)              :: outputForm
        type(outputDescriptor), intent(in)    :: descriptor
        integer, intent(in)                   :: id

        character(8)  :: today
        character(10) :: now
        integer                     :: ncid, varid, suffix,i
        integer                     :: DimId,lonDimId,latDimId,tileDimId,pftDimId,layerDimId,timeDimId
        real, dimension(2)          :: xrange, yrange
        real, dimension(1)          :: dummyArray = [ 0. ]
        integer, dimension(:), allocatable :: intArray

        character(:), pointer :: Comment   !< Comment about the run that will be written to the output netcdfs
        logical, pointer :: leap           !< set to true if all/some leap years in the .MET file have data for 366 days
                                           !< also accounts for leap years in .MET when cycling over meteorology (cyclemet)

        ! Point pointers
        leap => c_switch%leap
        Comment => c_switch%Comment

        ncid=netcdfVars(variableCount)%ncid

        ! ---------
        call ncPutAtt(ncid, nf90_global, 'title',charvalues='CLASSIC output file')
        call date_and_time(today,now)

        call ncPutAtt(ncid,nf90_global,'timestamp',charvalues=today//' '//now(1:4))
        call ncPutAtt(ncid,nf90_global,'Conventions',charvalues='COARDS')
        call ncPutAtt(ncid,nf90_global,'node_offset',intvalues=1)

        !----1 - Longitude

        lonDimId = ncDefDim(ncid,'lon',myDomain%cntx)
        varid = ncDefVar(ncid,'lon',nf90_double,[lonDimId])
        call ncPutAtt(ncid,varid,'long_name',charvalues='longitude')
        call ncPutAtt(ncid,varid,'units',charvalues='degrees_east')
        !call ncPutAtt(ncid,varid,'actual_range',xrange) #FLAG need to find the xrange from all_lon.
        call ncPutAtt(ncid,varid,'_Storage',charvalues="contiguous")


        !----2 - Latitude
        latDimId = ncDefDim(ncid,'lat',myDomain%cnty)
        varid = ncDefVar(ncid,'lat',nf90_double,[latDimId])
        call ncPutAtt(ncid,varid,'long_name',charvalues='latitude')
        call ncPutAtt(ncid,varid,'units',charvalues='degrees_north')
        !call ncPutAtt(ncid,varid,'actual_range',yrange) #FLAG need to find the xrange from all_lon.
        call ncPutAtt(ncid,varid,'_Storage',charvalues="contiguous")

        select case(trim(outputForm))

            case ("tile")       ! Per tile outputs

                tileDimId = ncDefDim(ncid,'tile',nmos)
                varid = ncDefVar(ncid,'tile',nf90_short,[tileDimId])
                call ncPutAtt(ncid,varid,'long_name',charvalues='tile')
                call ncPutAtt(ncid,varid,'units',charvalues='tile number')
                call ncPutAtt(ncid,varid,'_Storage',charvalues="contiguous")
                call ncPutAtt(ncid,varid,'_Endianness',charvalues="little")

            case("pft")         ! Per PFT outputs

                if (descriptor%includeBareGround) then
                    pftDimId = ncDefDim(ncid,'pft',iccp1)
                else
                    pftDimId = ncDefDim(ncid,'pft',icc)
                end if

                varid = ncDefVar(ncid,'pft',nf90_short,[pftDimId])
                call ncPutAtt(ncid,varid,'long_name',charvalues='Plant Functional Type')
                call ncPutAtt(ncid,varid,'units',charvalues='PFT')
                call ncPutAtt(ncid,varid,'_Storage',charvalues="contiguous")
                call ncPutAtt(ncid,varid,'_Endianness',charvalues="little")

            case ("layer")      ! Per layer outputs

                layerDimId = ncDefDim(ncid,'layer',ignd)
                varid = ncDefVar(ncid,'layer',nf90_short,[layerDimId])
                call ncPutAtt(ncid,varid,'long_name',charvalues='soil layer')
                call ncPutAtt(ncid,varid,'units',charvalues='layer')
                call ncPutAtt(ncid,varid,'_Storage',charvalues="contiguous")
                call ncPutAtt(ncid,varid,'_Endianness',charvalues="little")

        end select

        ! Set up the time dimension
        timeDimId = ncDefDim(ncid,'time',nf90_unlimited)
        varid = ncDefVar(ncid,'time',nf90_double,[timeDimId])

        call ncPutAtt(ncid,varid,'long_name',charvalues='time')
        call ncPutAtt(ncid,varid,'units',charvalues=trim(timestart))

        if (leap) then
            call ncPutAtt(ncid,varid,'calendar',charvalues="standard")
        else
            call ncPutAtt(ncid,varid,'calendar',charvalues="365_day")
        end if

        call ncPutAtt(ncid,varid,'_Storage',charvalues="chunked")
        call ncPutAtt(ncid,varid,'_Chunksizes',intvalues=1)
        call ncPutAtt(ncid,varid,'_Endianness',charvalues="little")

        call ncEndDef(ncid)

        ! Fill in the dimension variables and define the model output vars
        call ncPutDimValues(ncid, 'lon', myDomain%lonUnique, count=(/myDomain%cntx/))
        call ncPutDimValues(ncid, 'lat', myDomain%latUnique, count=(/myDomain%cnty/))
        
        select case(trim(outputForm))
            case ("tile")       ! Per tile outputs
                allocate(intArray(nmos))
                intArray=identityVector(nmos)
                call ncPutDimValues(ncid, 'tile', dummyArray,intValues=intArray, count=(/nmos/)) ! pass dummyArray to allow integers
                call ncReDef(ncid)
                varid = ncDefVar(ncid, trim(descriptor%shortName), nf90_double, [lonDimId,latDimId,tileDimId,timeDimId])

            case("pft")         ! Per PFT outputs

                if (descriptor%includeBareGround) then
                    allocate(intArray(iccp1))
                    intArray=identityVector(iccp1)
                    call ncPutDimValues(ncid, 'pft', dummyArray,intValues=intArray, count=(/iccp1/)) ! pass dummyArray to allow integers
                else
                    allocate(intArray(icc))
                    intArray=identityVector(icc)
                    call ncPutDimValues(ncid, 'pft', dummyArray,intValues=intArray, count=(/icc/)) ! pass dummyArray to allow integers
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
                call ncPutDimValues(ncid, 'layer', dummyArray,intValues=intArray, count=(/ignd/)) ! pass dummyArray to allow integers
                call ncReDef(ncid)
                varid = ncDefVar(ncid, trim(descriptor%shortName), nf90_double, [lonDimId,latDimId,layerDimId,timeDimId])

            case default        ! Grid average outputs

                call ncReDef(ncid)
                varid = ncDefVar(ncid, trim(descriptor%shortName), nf90_double, [lonDimId,latDimId,timeDimId])

        end select

        call ncPutAtt(ncid,varid,'long_name',charvalues=descriptor%longName)
        call ncPutAtt(ncid,varid,'units',charvalues=descriptor%units)
        call ncPutAtt(ncid,varid,'_FillValue',realvalues=fill_value)
        call ncPutAtt(ncid,varid,'missing_value',realvalues=fill_value)
        call ncPutAtt(ncid,varid,'_Storage',charvalues="chunked")
        call ncPutAtt(ncid,varid,'_DeflateLevel',intvalues=1)
        call ncPutAtt(ncid,nf90_global,'Comment',c_switch%Comment)
        call ncEndDef(ncid)

    end subroutine createNetCDF

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

        !print*,key,timeStamp,label
        id= getIdByKey(key)

        if (id == 0) then
            errmsg = 'writeOutput1D says: Your requested key does not exist (' // trim(key) // ') in netcdfVars.'
            print*,trim(errmsg)
            stop
        end if

        ncid = netcdfVars(id)%ncid

        ! Check if the time period has already been added to the file
        timeIndex = ncGetDimLen(ncid, "time")

        if (timeIndex == 0) then
            ! This is the first time step so add it.
            timeIndex = timeIndex + 1
            call ncPutDimValues(ncid, "time", localStamp, start=(/timeIndex/), count=(/1/))
        else
            ! This is a subsequent time step so need to check if it has already been added by
            ! another grid cell.
            allocate(timeWritten(timeIndex))
            timeWritten= ncGetDimValues(ncid, "time", count = (/timeIndex/))
            posTimeWanted = checkForTime(timeIndex,timeWritten,localStamp(1))

            if (posTimeWanted == 0) then ! Need to add this time step
                timeIndex = timeIndex + 1
                call ncPutDimValues(ncid, "time", localStamp, start=(/timeIndex/), count=(/1/))
            else
                ! timeStamp already added so just use it.
                timeIndex = posTimeWanted
            end if
        end if

        length = size(data)
        if (present(specStart)) then
            start = specStart
        else
            start = 1
        end if
        if (length > 1) then
            call ncPutVar(ncid, label, localData, start=[lonLocalIndex,latLocalIndex,start,timeIndex], count=[1,1,length,1])
        else
            call ncPutVar(ncid, label, localData, start=[lonLocalIndex,latLocalIndex,timeIndex], count=[1,1,1])
        end if
    end subroutine writeOutput1D

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

    pure function identityVector(n) result(res)
        integer, allocatable ::res(:)
        integer, intent(in) :: n
        integer             :: i
        allocate(res(n))
        forall (i=1:n)
            res(i) = i
        end forall
    end function identityVector


end module outputManager
