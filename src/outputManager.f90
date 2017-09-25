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
    public  :: closeNCFiles


    type simulationDomain
        real, dimension(:), allocatable     :: lonLandCell, latLandCell     ! Long/Lat values of only the land cells in our model domain
        integer, dimension(:), allocatable  :: lonLandIndex, latLandIndex   ! Indexes of only the land cells in our model domain for our resolution
        real, dimension(:), allocatable     :: allLonValues, allLatValues   ! All long/Lat values in our model domain (including ocean/non-land)
        integer                             :: LandCellCount    !> number of land cells that the model will run over
        real, dimension(4) :: domainBounds                      !> Corners of the domain to be simulated (netcdfs)
        integer :: srtx                                         !> starting index for this simulation for longitudes
        integer :: srty                                         !> starting index for this simulation for latitudes
        integer :: cntx                                         !> number of grid cells for this simulation in the longitude direction
        integer :: cnty                                         !> number of grid cells for this simulation in the latitude direction
    end type

    type(simulationDomain) :: myDomain

    integer :: metfid                               !> netcdf file id for the meteorology file
    integer :: initid                               !> netcdf file id for the model initialization file
    integer :: rsid                                 !> netcdf file id for the model restart file

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

    type netcdfVar
        integer         :: ncid
        character(30)   :: key
        character(80)   :: filename
    end type

    integer, parameter  :: maxncVariableNumber = 300
    type(netcdfVar)     :: netcdfVars(maxncVariableNumber)

    integer :: variableCount = 0, descriptorCount = 0

    integer         :: refyr = 1850                     !< !< Time reference for netcdf output files
    character(30)   :: timestart = "days since 1850-01-01 00:00" !< Time reference for netcdf output files
    character(30)   :: fill_value = "1.e38"             !< Default fill value for missing values in netcdf output files

contains

    subroutine generateOutputFiles

        implicit none

                            !nameInCode, timeFreq, outputForm, descriptor
        call generateNetCDFFile("fsstar_yr", "annually", "grid", "rss")
        call generateNetCDFFile("fsstar_mo", "monthly", "grid", "rss")
        call generateNetCDFFile("flstar_yr", "annually", "grid", "rls")
        call generateNetCDFFile("flstar_mo", "monthly", "grid", "rls")
        call generateNetCDFFile("qh_yr", "annually", "grid", "hfss")
        call generateNetCDFFile("qh_mo", "monthly", "grid", "hfss")
        call generateNetCDFFile("qe_yr", "annually", "grid", "hfls")
        call generateNetCDFFile("qe_mo", "monthly", "grid", "hfls")
        call generateNetCDFFile("snoacc_mo", "monthly", "grid", "snw")
        call generateNetCDFFile("wsnoacc_mo", "monthly", "grid", "wsnw")
        call generateNetCDFFile("rofacc_yr", "annually", "grid", "mrro")
        call generateNetCDFFile("rofacc_mo", "monthly", "grid", "mrro")
        call generateNetCDFFile("preacc_yr", "annually", "grid", "pr")
        call generateNetCDFFile("preacc_mo", "monthly", "grid", "pr")
        call generateNetCDFFile("evapacc_yr", "annually", "grid", "evspsbl")
        call generateNetCDFFile("evapacc_mo", "monthly", "grid", "evspsbl")
        call generateNetCDFFile("groundevap", "monthly", "grid", "evspsblsoi")
        call generateNetCDFFile("canopyevap", "monthly", "grid", "evspsblveg")
        call generateNetCDFFile("transpacc_yr", "annually", "grid", "tran")
        call generateNetCDFFile("transpacc_mo", "monthly", "grid", "tran")
        call generateNetCDFFile("taacc_mo", "monthly", "grid", "ts")
        call generateNetCDFFile("altotacc_yr", "annually", "grid", "albs")
        call generateNetCDFFile("altotacc_mo", "monthly", "grid", "albs")
        call generateNetCDFFile("TBARACC_MO", "monthly", "layer", "tsl")
        call generateNetCDFFile("THLQACC_MO", "monthly", "layer", "mrsll")
        call generateNetCDFFile("THICACC_MO", "monthly", "layer", "mrsfl")
        call generateNetCDFFile("ACTLYR_MO", "monthly", "grid", "actlyr")
        call generateNetCDFFile("ACTLYR_MAX_MO", "monthly", "grid", "actlyrmax")
        call generateNetCDFFile("ACTLYR_MIN_MO", "monthly", "grid", "actlyrmin")
        call generateNetCDFFile("FTABLE_MO", "monthly", "grid", "ftable")
        call generateNetCDFFile("FTABLE_MAX_MO", "monthly", "grid", "ftablemax")
        call generateNetCDFFile("FTABLE_MIN_MO", "monthly", "grid", "ftablemin")
        call generateNetCDFFile("laimaxg_yr_g", "annually", "grid", "lai")
        call generateNetCDFFile("laimaxg_mo_g", "monthly", "grid", "lai")
        call generateNetCDFFile("laimaxg_yr", "annually", "pft", "lai")
        call generateNetCDFFile("laimaxg_mo", "monthly", "pft", "lai")
        call generateNetCDFFile("laimaxg_yr_t", "annually", "tile", "lai")
        call generateNetCDFFile("laimaxg_mo_t", "monthly", "tile", "lai")
        call generateNetCDFFile("vgbiomas_yr_g", "annually", "grid", "cVeg")
        call generateNetCDFFile("vgbiomas_mo_g", "monthly", "grid", "cVeg")
        call generateNetCDFFile("vgbiomas_yr", "annually", "pft", "cVeg")
        call generateNetCDFFile("vgbiomas_mo", "monthly", "pft", "cVeg")
        call generateNetCDFFile("vgbiomas_yr_t", "annually", "tile", "cVeg")
        call generateNetCDFFile("vgbiomas_mo_t", "monthly", "tile", "cVeg")
        call generateNetCDFFile("stemmass_yr_g", "annually", "grid", "cStem")
        call generateNetCDFFile("stemmass_yr", "annually", "pft", "cStem")
        call generateNetCDFFile("stemmass_yr_t", "annually", "tile", "cStem")
        call generateNetCDFFile("rootmass_yr_g", "annually", "grid", "cRoot")
         call generateNetCDFFile("rootmass_yr", "annually", "pft", "cRoot")
        call generateNetCDFFile("rootmass_yr_t", "annually", "tile", "cRoot")
        call generateNetCDFFile("litrmass_yr_g", "annually", "grid", "cLitter")
        call generateNetCDFFile("litrmass_mo_t", "monthly", "grid", "cLitter")
        call generateNetCDFFile("litrmass_yr", "annually", "pft", "cLitter")
        call generateNetCDFFile("litrmass_mo", "monthly", "pft", "cLitter")
        call generateNetCDFFile("litrmass_yr_t", "annually", "tile", "cLitter")
        call generateNetCDFFile("litrmass_mo_t", "monthly", "tile", "cLitter")
        call generateNetCDFFile("soilcmas_yr_g", "annually", "grid", "cSoil")
        call generateNetCDFFile("soilcmas_mo_g", "monthly", "grid", "cSoil")
        call generateNetCDFFile("soilcmas_yr", "annually", "pft", "cSoil")
        call generateNetCDFFile("soilcmas_mo", "monthly", "pft", "cSoil")
        call generateNetCDFFile("soilcmas_yr_t", "annually", "tile", "cSoil")
        call generateNetCDFFile("soilcmas_mo_t", "monthly", "tile", "cSoil")
        call generateNetCDFFile("totcmass_yr_g", "annually", "grid", "cLand")
        call generateNetCDFFile("totcmass_yr", "annually", "pft", "cLand")
        call generateNetCDFFile("totcmass_yr_t", "annually", "tile", "cLand")
        call generateNetCDFFile("veghght_yr_g", "annually", "grid", "vegHeight")
        call generateNetCDFFile("veghght_yr", "annually", "pft", "vegHeight")
        call generateNetCDFFile("veghght_yr_t", "annually", "tile", "vegHeight")
        call generateNetCDFFile("npp_yr_g", "annually", "grid", "npp")
        call generateNetCDFFile("npp_mo_g", "monthly", "grid", "npp")
        call generateNetCDFFile("npp_yr", "annually", "pft", "npp")
        call generateNetCDFFile("npp_mo", "monthly", "pft", "npp")
        call generateNetCDFFile("npp_yr_t", "annually", "tile", "npp")
        call generateNetCDFFile("npp_mo_t", "monthly", "tile", "npp")
        call generateNetCDFFile("gpp_yr_g", "annually", "grid", "gpp")
        call generateNetCDFFile("gpp_mo_g", "monthly", "grid", "gpp")
        call generateNetCDFFile("gpp_yr", "annually", "pft", "gpp")
        call generateNetCDFFile("gpp_mo", "monthly", "pft", "gpp")
        call generateNetCDFFile("gpp_yr_t", "annually", "tile", "gpp")
        call generateNetCDFFile("gpp_mo_t", "monthly", "tile", "gpp")
        call generateNetCDFFile("nep_yr_g", "annually", "grid", "nep")
        call generateNetCDFFile("nep_mo_g", "monthly", "grid", "nep")
        call generateNetCDFFile("nep_yr", "annually", "pft", "nep")
        call generateNetCDFFile("nep_mo", "monthly", "pft", "nep")
        call generateNetCDFFile("nep_yr_t", "annually", "tile", "nep")
        call generateNetCDFFile("nep_mo_t", "monthly", "tile", "nep")
        call generateNetCDFFile("nbp_yr_g", "annually", "grid", "nbp")
        call generateNetCDFFile("nbp_mo_g", "monthly", "grid", "nbp")
        call generateNetCDFFile("nbp_yr", "annually", "pft", "nbp")
        call generateNetCDFFile("nbp_mo", "monthly", "pft", "nbp")
        call generateNetCDFFile("nbp_yr_t", "annually", "tile", "nbp")
        call generateNetCDFFile("nbp_mo_t", "monthly", "tile", "nbp")
        call generateNetCDFFile("hetrores_yr_g", "annually", "grid", "rh")
        call generateNetCDFFile("hetrores_mo_g", "monthly", "grid", "rh")
        call generateNetCDFFile("hetrores_yr", "annually", "pft", "rh")
        call generateNetCDFFile("hetrores_mo", "monthly", "pft", "rh")
        call generateNetCDFFile("hetrores_yr_t", "annually", "tile", "rh")
        call generateNetCDFFile("hetrores_mo_t", "monthly", "tile", "rh")
        call generateNetCDFFile("autores_yr_g", "annually", "grid", "ra")
        call generateNetCDFFile("autores_mo_g", "monthly", "grid", "ra")
        call generateNetCDFFile("autores_yr", "annually", "pft", "ra")
        call generateNetCDFFile("autores_mo", "monthly", "pft", "ra")
        call generateNetCDFFile("autores_yr_t", "annually", "tile", "ra")
        call generateNetCDFFile("autores_mo_t", "monthly", "tile", "ra")
        call generateNetCDFFile("litres_yr_g", "annually", "grid", "rhLitter")
        call generateNetCDFFile("litres_mo_g", "monthly", "grid", "rhLitter")
        call generateNetCDFFile("litres_yr", "annually", "pft", "rhLitter")
        call generateNetCDFFile("litres_mo", "monthly", "pft", "rhLitter")
        call generateNetCDFFile("litres_yr_t", "annually", "tile", "rhLitter")
        call generateNetCDFFile("litres_mo_t", "monthly", "tile", "rhLitter")
        call generateNetCDFFile("soilcres_yr_g", "annually", "grid", "rhSoil")
        call generateNetCDFFile("soilcres_mo_g", "monthly", "grid", "rhSoil")
        call generateNetCDFFile("soilcres_yr", "annually", "pft", "rhSoil")
        call generateNetCDFFile("soilcres_mo", "monthly", "pft", "rhSoil")
        call generateNetCDFFile("soilcres_yr_t", "annually", "tile", "rhSoil")
        call generateNetCDFFile("soilcres_mo_t", "monthly", "tile", "rhSoil")
        call generateNetCDFFile("litrfallveg_mo_g", "monthly", "grid", "fVegLitter")
        call generateNetCDFFile("litrfallveg_mo", "monthly", "pft", "fVegLitter")
        call generateNetCDFFile("litrfallveg_mo_t", "monthly", "tile", "fVegLitter")
        call generateNetCDFFile("humiftrsveg_mo_g", "monthly", "grid", "fLitterSoil")
        call generateNetCDFFile("humiftrsveg_mo", "monthly", "pft", "fLitterSoil")
        call generateNetCDFFile("humiftrsveg_mo_t", "monthly", "tile", "fLitterSoil")
        call generateNetCDFFile("emit_co2_yr_g", "annually", "grid", "fFire")
        call generateNetCDFFile("emit_co2_mo_g", "monthly", "grid", "fFire")
        call generateNetCDFFile("emit_co2_yr", "annually", "pft", "fFire")
        call generateNetCDFFile("emit_co2_mo", "monthly", "pft", "fFire")
        call generateNetCDFFile("emit_co2_yr_t", "annually", "tile", "fFire")
        call generateNetCDFFile("emit_co2_mo_t", "monthly", "tile", "fFire")
        call generateNetCDFFile("emit_ch4_yr_g", "annually", "grid", "fFireCH4")
        call generateNetCDFFile("emit_ch4_mo_g", "monthly", "grid", "fFireCH4")
        call generateNetCDFFile("emit_ch4_yr", "annually", "pft", "fFireCH4")
        call generateNetCDFFile("emit_ch4_mo", "monthly", "pft", "fFireCH4")
        call generateNetCDFFile("emit_ch4_yr_t", "annually", "tile", "fFireCH4")
        call generateNetCDFFile("emit_ch4_mo_t", "monthly", "tile", "fFireCH4")
        call generateNetCDFFile("burnfrac_yr_g", "annually", "grid", "burntFractionAll")
        call generateNetCDFFile("burnfrac_mo_g", "monthly", "grid", "burntFractionAll")
        call generateNetCDFFile("burnfrac_yr", "annually", "pft", "burntFractionAll")
        call generateNetCDFFile("burnfrac_mo", "monthly", "pft", "burntFractionAll")
        call generateNetCDFFile("burnfrac_yr_t", "annually", "tile", "burntFractionAll")
        call generateNetCDFFile("burnfrac_mo_t", "monthly", "tile", "burntFractionAll")
        call generateNetCDFFile("luc_emc_yr_g", "annually", "grid", "fDeforestToAtmos")
        call generateNetCDFFile("luc_emc_mo_g", "monthly", "grid", "fDeforestToAtmos")
        call generateNetCDFFile("luc_emc_yr", "annually", "pft", "fDeforestToAtmos")
        call generateNetCDFFile("luc_emc_mo", "monthly", "pft", "fDeforestToAtmos")
        call generateNetCDFFile("luc_emc_yr_t", "annually", "tile", "fDeforestToAtmos")
        call generateNetCDFFile("luc_emc_mo_t", "monthly", "tile", "fDeforestToAtmos")
        call generateNetCDFFile("lucltrin_yr_g", "annually", "grid", "fDeforestToLitter")
        call generateNetCDFFile("lucltrin_mo_g", "monthly", "grid", "fDeforestToLitter")
        call generateNetCDFFile("lucltrin_yr", "annually", "pft", "fDeforestToLitter")
        call generateNetCDFFile("lucltrin_mo", "monthly", "pft", "fDeforestToLitter")
        call generateNetCDFFile("lucltrin_yr_t", "annually", "tile", "fDeforestToLitter")
        call generateNetCDFFile("lucltrin_mo_t", "monthly", "tile", "fDeforestToLitter")
        call generateNetCDFFile("lucsocin_yr_g", "annually", "grid", "fDeforestToSoil")
        call generateNetCDFFile("lucsocin_mo_g", "monthly", "grid", "fDeforestToSoil")
        call generateNetCDFFile("lucsocin_yr", "annually", "pft", "fDeforestToSoil")
        call generateNetCDFFile("lucsocin_mo", "monthly", "pft", "fDeforestToSoil")
        call generateNetCDFFile("lucsocin_yr_t", "annually", "tile", "fDeforestToSoil")
        call generateNetCDFFile("lucsocin_mo_t", "monthly", "tile", "fDeforestToSoil")
        call generateNetCDFFile("fcancmxrow_yr_g", "annually", "grid", "landCoverFrac")
        call generateNetCDFFile("fcancmxrow_mo_g", "monthly", "grid", "landCoverFrac")
        call generateNetCDFFile("pftexistrow_yr_g", "annually", "grid", "landCoverExist")
        call generateNetCDFFile("pftexistrow_mo_g", "monthly", "grid", "landCoverExist")
        call generateNetCDFFile("ch4wet1_yr_g", "annually", "grid", "wetlandCH4spec")
        call generateNetCDFFile("ch4wet1_mo_g", "monthly", "grid", "wetlandCH4spec")
        call generateNetCDFFile("ch4wet1_yr_t", "annually", "tile", "wetlandCH4spec")
        call generateNetCDFFile("ch4wet1_mo_t", "monthly", "tile", "wetlandCH4spec")
        call generateNetCDFFile("ch4dyn1_yr_g", "annually", "grid", "wetlandCH4dyn")
        call generateNetCDFFile("ch4dyn1_mo_g", "monthly", "grid", "wetlandCH4dyn")
        call generateNetCDFFile("ch4dyn1_yr_t", "annually", "tile", "wetlandCH4dyn")
        call generateNetCDFFile("ch4dyn1_mo_t", "monthly", "tile", "wetlandCH4dyn")
        call generateNetCDFFile("ch4soills_yr_g", "annually", "grid", "soilCH4cons")
        call generateNetCDFFile("ch4soills_mo_g", "monthly", "grid", "soilCH4cons")
        call generateNetCDFFile("ch4soills_yr_t", "annually", "tile", "soilCH4cons")
        call generateNetCDFFile("ch4soills_mo_t", "monthly", "tile", "soilCH4cons")
        call generateNetCDFFile("wetfdyn_yr_g", "annually", "grid", "wetlandFrac")
        call generateNetCDFFile("wetfdyn_mo_g", "monthly", "grid", "wetlandFrac")
        call generateNetCDFFile("wetfdyn_yr_t", "annually", "tile", "wetlandFrac")
        call generateNetCDFFile("wetfdyn_mo_t", "monthly", "tile", "wetlandFrac")

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
        elseif (c_switch%compete .and. trim(descriptor%group) == "compete") then
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

    subroutine createNetCDF(fileName,id, outputForm, descriptor)

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
        varid = ncDefVar(ncid,'lon',nf90_float,[lonDimId])
        call ncPutAtt(ncid,varid,'long_name',charvalues='longitude')
        call ncPutAtt(ncid,varid,'units',charvalues='degrees_east')
        !call ncPutAtt(ncid,varid,'actual_range',xrange) #FLAG need to find the xrange from all_lon.
        call ncPutAtt(ncid,varid,'_Storage',charvalues="contiguous")


        !----2 - Latitude
        latDimId = ncDefDim(ncid,'lat',myDomain%cnty)
        varid = ncDefVar(ncid,'lat',nf90_float,[latDimId])
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
        varid = ncDefVar(ncid,'time',nf90_float,[timeDimId])

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
        call ncPutDimValues(ncid, 'lon', myDomain%lonLandCell, count=myDomain%cntx)
        call ncPutDimValues(ncid, 'lat', myDomain%latLandCell, count=myDomain%cnty)
        
        select case(trim(outputForm))
            case ("tile")       ! Per tile outputs
                allocate(intArray(nmos))
                intArray=identityVector(nmos)
                call ncPutDimValues(ncid, 'tile', dummyArray,intdata=intArray, count=nmos) ! pass dummyArray to allow integers
                call ncReDef(ncid)
                varid = ncDefVar(ncid, trim(descriptor%shortName), nf90_float, [lonDimId,latDimId,tileDimId,timeDimId])

            case("pft")         ! Per PFT outputs

                if (descriptor%includeBareGround) then
                    allocate(intArray(iccp1))
                    intArray=identityVector(iccp1)
                    call ncPutDimValues(ncid, 'pft', dummyArray,intdata=intArray, count=iccp1) ! pass dummyArray to allow integers
                else
                    allocate(intArray(icc))
                    intArray=identityVector(icc)
                    call ncPutDimValues(ncid, 'pft', dummyArray,intdata=intArray, count=icc) ! pass dummyArray to allow integers
                end if


                call ncReDef(ncid)

                if (descriptor%includeBareGround) then
                    ! do something for cells that have bare ground
                    suffix = suffix + 1
                else
                    ! do something for cells that don't have bare ground
                    suffix = suffix
                endif
                varid = ncDefVar(ncid, trim(descriptor%shortName), nf90_float, [lonDimId,latDimId,pftDimId,timeDimId])

            case ("layer")      ! Per layer outputs

                allocate(intArray(ignd))
                intArray=identityVector(ignd)
                call ncPutDimValues(ncid, 'layer', dummyArray,intdata=intArray, count=ignd) ! pass dummyArray to allow integers
                call ncReDef(ncid)
                varid = ncDefVar(ncid, trim(descriptor%shortName), nf90_float, [lonDimId,latDimId,layerDimId,timeDimId])

            case default        ! Grid average outputs

                call ncReDef(ncid)
                varid = ncDefVar(ncid, trim(descriptor%shortName), nf90_float, [lonDimId,latDimId,timeDimId])

        end select

        call ncPutAtt(ncid,varid,'long_name',charvalues=descriptor%longName)
        call ncPutAtt(ncid,varid,'units',charvalues=descriptor%units)

        !call ncPutAtt(ncid,varid,'_FillValue',fill_value) !FLAG this is not working at present.

        call ncPutAtt(ncid,varid,'missing_value',charvalues=fill_value)
        call ncPutAtt(ncid,varid,'_Storage',charvalues="chunked")
        call ncPutAtt(ncid,varid,'_DeflateLevel',intvalues=1)
        !call ncPutAtt(ncid,varid,'name_in_code',varinfo%nameincode)
        !call ncPutAtt(ncid,varid,'Comment',varinfo%Comment)
        call ncEndDef(ncid)

    end subroutine createNetCDF

    subroutine writeOutput1D(key,timeStamp,label,data)
        use fileIOModule
        implicit none

        character(*), intent(in) :: key
        integer :: ncid, timeIndex, id, length
        real, intent(in) :: timeStamp
        real, dimension(:), intent(in)  :: data
        character(*), intent(in)        :: label

        !print*,key,timeStamp,label
        id= getIdByKey(key)
        ncid = netcdfVars(id)%ncid
        timeIndex = ncGetDimLen(ncid, "time") + 1
        call ncPutDimValues(ncid, "time", [timeStamp], start=timeIndex, count=1)
        length = size(data)
        if (length > 1) then
            call ncPut1DVar(ncid, label, data, start=[1,1,1,timeIndex], count=[1,1,length,1])
        else
            call ncPut1DVar(ncid, label, data, start=[1,1,timeIndex], count=[1,1,1])
        end if
    end subroutine writeOutput1D

    subroutine closeNCFiles(incid)
        use fileIOModule
        implicit none
        integer, optional   :: incid
        integer i
        ! close all output netcdfs or just a select file
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
