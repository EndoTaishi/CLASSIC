module outputManager

    use io_driver, only : outputDescriptor
    use netcdfFileManager, only :

    implicit none

contains

    subroutine generateOutputVariables

        implicit none

                           ! nameInCode, timeFreq, outputForm, descriptor
        call generateVariable("fsstar_yr", "annually", "grid", "rss")
        call generateVariable("fsstar_mo", "monthly", "grid", "rss")
        call generateVariable("flstar_yr", "annually", "grid", "rls")
        call generateVariable("flstar_mo", "monthly", "grid", "rls")
        call generateVariable("qh_yr", "annually", "grid", "hfss")
        call generateVariable("qh_mo", "monthly", "grid", "hfss")
        call generateVariable("qe_yr", "annually", "grid", "hfls")
        call generateVariable("qe_mo", "monthly", "grid", "hfls")
        call generateVariable("snoacc_mo", "monthly", "grid", "snw")
        call generateVariable("wsnoacc_mo", "monthly", "grid", "wsnw")
        call generateVariable("rofacc_yr", "annually", "grid", "mrro")
        call generateVariable("rofacc_mo", "monthly", "grid", "mrro")
        call generateVariable("preacc_yr", "annually", "grid", "pr")
        call generateVariable("preacc_mo", "monthly", "grid", "pr")
        call generateVariable("evapacc_yr", "annually", "grid", "evspsbl")
        call generateVariable("evapacc_mo", "monthly", "grid", "evspsbl")
        call generateVariable("groundevap", "monthly", "grid", "evspsblsoi")
        call generateVariable("canopyevap", "monthly", "grid", "evspsblveg")
        call generateVariable("transpacc_yr", "annually", "grid", "tran")
        call generateVariable("transpacc_mo", "monthly", "grid", "tran")
        call generateVariable("taacc_mo", "monthly", "grid", "ts")
        call generateVariable("altotacc_yr", "annually", "grid", "albs")
        call generateVariable("altotacc_mo", "monthly", "grid", "albs")
        call generateVariable("TBARACC_MO", "monthly", "layer", "tsl")
        call generateVariable("THLQACC_MO", "monthly", "layer", "mrsll")
        call generateVariable("THICACC_MO", "monthly", "layer", "mrsfl")
        call generateVariable("ACTLYR_MO", "monthly", "grid", "actlyr")
        call generateVariable("ACTLYR_MAX_MO", "monthly", "grid", "actlyrmax")
        call generateVariable("ACTLYR_MIN_MO", "monthly", "grid", "actlyrmin")
        call generateVariable("FTABLE_MO", "monthly", "grid", "ftable")
        call generateVariable("FTABLE_MAX_MO", "monthly", "grid", "ftablemax")
        call generateVariable("FTABLE_MIN_MO", "monthly", "grid", "ftablemin")
        call generateVariable("laimaxg_yr_g", "annually", "grid", "lai")
        call generateVariable("laimaxg_mo_g", "monthly", "grid", "lai")
        call generateVariable("laimaxg_yr", "annually", "pft", "lai")
        call generateVariable("laimaxg_mo", "monthly", "pft", "lai")
        call generateVariable("laimaxg_yr_t", "annually", "tile", "lai")
        call generateVariable("laimaxg_mo_t", "monthly", "tile", "lai")
        call generateVariable("vgbiomas_yr_g", "annually", "grid", "cVeg")
        call generateVariable("vgbiomas_mo_g", "monthly", "grid", "cVeg")
        call generateVariable("vgbiomas_yr", "annually", "pft", "cVeg")
        call generateVariable("vgbiomas_mo", "monthly", "pft", "cVeg")
        call generateVariable("vgbiomas_yr_t", "annually", "tile", "cVeg")
        call generateVariable("vgbiomas_mo_t", "monthly", "tile", "cVeg")
        call generateVariable("stemmass_yr_g", "annually", "grid", "cStem")
        call generateVariable("stemmass_yr", "annually", "pft", "cStem")
        call generateVariable("stemmass_yr_t", "annually", "tile", "cStem")
        call generateVariable("rootmass_yr_g", "annually", "grid", "cRoot")
        call generateVariable("rootmass_yr", "annually", "pft", "cRoot")
        call generateVariable("rootmass_yr_t", "annually", "tile", "cRoot")
        call generateVariable("litrmass_yr_g", "annually", "grid", "cLitter")
        call generateVariable("litrmass_mo_t", "monthly", "grid", "cLitter")
        call generateVariable("litrmass_yr", "annually", "pft", "cLitter")
        call generateVariable("litrmass_mo", "monthly", "pft", "cLitter")
        call generateVariable("litrmass_yr_t", "annually", "tile", "cLitter")
        call generateVariable("litrmass_mo_t", "monthly", "tile", "cLitter")
        call generateVariable("soilcmas_yr_g", "annually", "grid", "cSoil")
        call generateVariable("soilcmas_mo_g", "monthly", "grid", "cSoil")
        call generateVariable("soilcmas_yr", "annually", "pft", "cSoil")
        call generateVariable("soilcmas_mo", "monthly", "pft", "cSoil")
        call generateVariable("soilcmas_yr_t", "annually", "tile", "cSoil")
        call generateVariable("soilcmas_mo_t", "monthly", "tile", "cSoil")
        call generateVariable("totcmass_yr_g", "annually", "grid", "cLand")
        call generateVariable("totcmass_yr", "annually", "pft", "cLand")
        call generateVariable("totcmass_yr_t", "annually", "tile", "cLand")
        call generateVariable("veghght_yr_g", "annually", "grid", "vegHeight")
        call generateVariable("veghght_yr", "annually", "pft", "vegHeight")
        call generateVariable("veghght_yr_t", "annually", "tile", "vegHeight")
        call generateVariable("npp_yr_g", "annually", "grid", "npp")
        call generateVariable("npp_mo_g", "monthly", "grid", "npp")
        call generateVariable("npp_yr", "annually", "pft", "npp")
        call generateVariable("npp_mo", "monthly", "pft", "npp")
        call generateVariable("npp_yr_t", "annually", "tile", "npp")
        call generateVariable("npp_mo_t", "monthly", "tile", "npp")
        call generateVariable("gpp_yr_g", "annually", "grid", "gpp")
        call generateVariable("gpp_mo_g", "monthly", "grid", "gpp")
        call generateVariable("gpp_yr", "annually", "pft", "gpp")
        call generateVariable("gpp_mo", "monthly", "pft", "gpp")
        call generateVariable("gpp_yr_t", "annually", "tile", "gpp")
        call generateVariable("gpp_mo_t", "monthly", "tile", "gpp")
        call generateVariable("nep_yr_g", "annually", "grid", "nep")
        call generateVariable("nep_mo_g", "monthly", "grid", "nep")
        call generateVariable("nep_yr", "annually", "pft", "nep")
        call generateVariable("nep_mo", "monthly", "pft", "nep")
        call generateVariable("nep_yr_t", "annually", "tile", "nep")
        call generateVariable("nep_mo_t", "monthly", "tile", "nep")
        call generateVariable("nbp_yr_g", "annually", "grid", "nbp")
        call generateVariable("nbp_mo_g", "monthly", "grid", "nbp")
        call generateVariable("nbp_yr", "annually", "pft", "nbp")
        call generateVariable("nbp_mo", "monthly", "pft", "nbp")
        call generateVariable("nbp_yr_t", "annually", "tile", "nbp")
        call generateVariable("nbp_mo_t", "monthly", "tile", "nbp")
        call generateVariable("hetrores_yr_g", "annually", "grid", "rh")
        call generateVariable("hetrores_mo_g", "monthly", "grid", "rh")
        call generateVariable("hetrores_yr", "annually", "pft", "rh")
        call generateVariable("hetrores_mo", "monthly", "pft", "rh")
        call generateVariable("hetrores_yr_t", "annually", "tile", "rh")
        call generateVariable("hetrores_mo_t", "monthly", "tile", "rh")
        call generateVariable("autores_yr_g", "annually", "grid", "ra")
        call generateVariable("autores_mo_g", "monthly", "grid", "ra")
        call generateVariable("autores_yr", "annually", "pft", "ra")
        call generateVariable("autores_mo", "monthly", "pft", "ra")
        call generateVariable("autores_yr_t", "annually", "tile", "ra")
        call generateVariable("autores_mo_t", "monthly", "tile", "ra")
        call generateVariable("litres_yr_g", "annually", "grid", "rhLitter")
        call generateVariable("litres_mo_g", "monthly", "grid", "rhLitter")
        call generateVariable("litres_yr", "annually", "pft", "rhLitter")
        call generateVariable("litres_mo", "monthly", "pft", "rhLitter")
        call generateVariable("litres_yr_t", "annually", "tile", "rhLitter")
        call generateVariable("litres_mo_t", "monthly", "tile", "rhLitter")
        call generateVariable("soilcres_yr_g", "annually", "grid", "rhSoil")
        call generateVariable("soilcres_mo_g", "monthly", "grid", "rhSoil")
        call generateVariable("soilcres_yr", "annually", "pft", "rhSoil")
        call generateVariable("soilcres_mo", "monthly", "pft", "rhSoil")
        call generateVariable("soilcres_yr_t", "annually", "tile", "rhSoil")
        call generateVariable("soilcres_mo_t", "monthly", "tile", "rhSoil")
        call generateVariable("litrfallveg_mo_g", "monthly", "grid", "fVegLitter")
        call generateVariable("litrfallveg_mo", "monthly", "pft", "fVegLitter")
        call generateVariable("litrfallveg_mo_t", "monthly", "tile", "fVegLitter")
        call generateVariable("humiftrsveg_mo_g", "monthly", "grid", "fLitterSoil")
        call generateVariable("humiftrsveg_mo", "monthly", "pft", "fLitterSoil")
        call generateVariable("humiftrsveg_mo_t", "monthly", "tile", "fLitterSoil")
        call generateVariable("emit_co2_yr_g", "annually", "grid", "fFire")
        call generateVariable("emit_co2_mo_g", "monthly", "grid", "fFire")
        call generateVariable("emit_co2_yr", "annually", "pft", "fFire")
        call generateVariable("emit_co2_mo", "monthly", "pft", "fFire")
        call generateVariable("emit_co2_yr_t", "annually", "tile", "fFire")
        call generateVariable("emit_co2_mo_t", "monthly", "tile", "fFire")
        call generateVariable("emit_ch4_yr_g", "annually", "grid", "fFireCH4")
        call generateVariable("emit_ch4_mo_g", "monthly", "grid", "fFireCH4")
        call generateVariable("emit_ch4_yr", "annually", "pft", "fFireCH4")
        call generateVariable("emit_ch4_mo", "monthly", "pft", "fFireCH4")
        call generateVariable("emit_ch4_yr_t", "annually", "tile", "fFireCH4")
        call generateVariable("emit_ch4_mo_t", "monthly", "tile", "fFireCH4")
        call generateVariable("burnfrac_yr_g", "annually", "grid", "burntFractionAll")
        call generateVariable("burnfrac_mo_g", "monthly", "grid", "burntFractionAll")
        call generateVariable("burnfrac_yr", "annually", "pft", "burntFractionAll")
        call generateVariable("burnfrac_mo", "monthly", "pft", "burntFractionAll")
        call generateVariable("burnfrac_yr_t", "annually", "tile", "burntFractionAll")
        call generateVariable("burnfrac_mo_t", "monthly", "tile", "burntFractionAll")
        call generateVariable("luc_emc_yr_g", "annually", "grid", "fDeforestToAtmos")
        call generateVariable("luc_emc_mo_g", "monthly", "grid", "fDeforestToAtmos")
        call generateVariable("luc_emc_yr", "annually", "pft", "fDeforestToAtmos")
        call generateVariable("luc_emc_mo", "monthly", "pft", "fDeforestToAtmos")
        call generateVariable("luc_emc_yr_t", "annually", "tile", "fDeforestToAtmos")
        call generateVariable("luc_emc_mo_t", "monthly", "tile", "fDeforestToAtmos")
        call generateVariable("lucltrin_yr_g", "annually", "grid", "fDeforestToLitter")
        call generateVariable("lucltrin_mo_g", "monthly", "grid", "fDeforestToLitter")
        call generateVariable("lucltrin_yr", "annually", "pft", "fDeforestToLitter")
        call generateVariable("lucltrin_mo", "monthly", "pft", "fDeforestToLitter")
        call generateVariable("lucltrin_yr_t", "annually", "tile", "fDeforestToLitter")
        call generateVariable("lucltrin_mo_t", "monthly", "tile", "fDeforestToLitter")
        call generateVariable("lucsocin_yr_g", "annually", "grid", "fDeforestToSoil")
        call generateVariable("lucsocin_mo_g", "monthly", "grid", "fDeforestToSoil")
        call generateVariable("lucsocin_yr", "annually", "pft", "fDeforestToSoil")
        call generateVariable("lucsocin_mo", "monthly", "pft", "fDeforestToSoil")
        call generateVariable("lucsocin_yr_t", "annually", "tile", "fDeforestToSoil")
        call generateVariable("lucsocin_mo_t", "monthly", "tile", "fDeforestToSoil")
        call generateVariable("fcancmxrow_yr_g", "annually", "grid", "landCoverFrac")
        call generateVariable("fcancmxrow_mo_g", "monthly", "grid", "landCoverFrac")
        call generateVariable("pftexistrow_yr_g", "annually", "grid", "landCoverExist")
        call generateVariable("pftexistrow_mo_g", "monthly", "grid", "landCoverExist")
        call generateVariable("ch4wet1_yr_g", "annually", "grid", "wetlandCH4spec")
        call generateVariable("ch4wet1_mo_g", "monthly", "grid", "wetlandCH4spec")
        call generateVariable("ch4wet1_yr_t", "annually", "tile", "wetlandCH4spec")
        call generateVariable("ch4wet1_mo_t", "monthly", "tile", "wetlandCH4spec")
        call generateVariable("ch4dyn1_yr_g", "annually", "grid", "wetlandCH4dyn")
        call generateVariable("ch4dyn1_mo_g", "monthly", "grid", "wetlandCH4dyn")
        call generateVariable("ch4dyn1_yr_t", "annually", "tile", "wetlandCH4dyn")
        call generateVariable("ch4dyn1_mo_t", "monthly", "tile", "wetlandCH4dyn")
        call generateVariable("ch4soills_yr_g", "annually", "grid", "soilCH4cons")
        call generateVariable("ch4soills_mo_g", "monthly", "grid", "soilCH4cons")
        call generateVariable("ch4soills_yr_t", "annually", "tile", "soilCH4cons")
        call generateVariable("ch4soills_mo_t", "monthly", "tile", "soilCH4cons")
        call generateVariable("wetfdyn_yr_g", "annually", "grid", "wetlandFrac")
        call generateVariable("wetfdyn_mo_g", "monthly", "grid", "wetlandFrac")
        call generateVariable("wetfdyn_yr_t", "annually", "tile", "wetlandFrac")
        call generateVariable("wetfdyn_mo_t", "monthly", "tile", "wetlandFrac")

    end subroutine generateOutputVariables

    subroutine generateVariable(nameInCode, timeFreq, outputForm, descriptorLabel)
    ! Generates a new (netcdf) variable

        implicit none

        character(*), intent(in)                :: nameInCode, timeFreq, outputForm, descriptorLabel
        type(outputDescriptor)                :: descriptor
        character(80)                           :: filename = ''
        logical                                 :: isTimeValid, isGroupValid
        integer                                 :: id

        !call printConfig(config)
        !call printDescriptor(descriptor)

        ! Get variable descriptor
        descriptor = getDescriptor(descriptorLabel)
        ! If the requested timeFreq matches the project config, all's good
        isTimeValid = validTime(config, timeFreq, descriptor)
        ! If the group property of descriptor matches the project config, all's good
        isGroupValid = validGroup(config, descriptor)
        ! If the project config and variable descriptor match the request, process current variable
        if (isTimeValid .and. isGroupValid) then
            ! Generate the filename
            filename = generateFilename(outputForm, config, descriptor)
            ! Allocate a new variable (ncid, filename, key etc.)
            id = addVariable(nameInCode, filename)
            ! Make the netcdf file for the new variable (mostly definitions)
            call createNetCDF(id, outputForm, descriptor)
        endif

    end subroutine generateVariable

    ! Adds the new variable to the list of variables (see the type "variable")
    integer function addVariable(key, filename)
        implicit none
        character(*), intent(in)    :: key, filename
        integer                     :: ncid
        ncid = createFile(filename)
        variableCount = variableCount + 1
        variables(variableCount)%ncid = ncid
        variables(variableCount)%key = key
        variables(variableCount)%filename = filename
        addVariable = variableCount
    end function addVariable

    ! Determines if the current variable matches the project configuration
    logical function validGroup(config, descriptor)
        implicit none
        type(projectConfiguration), intent(in)  :: config
        type(outputDescriptor), intent(in)    :: descriptor
        if (config%class .and. trim(descriptor%group) == "class") then
            validGroup = .true.
        elseif (config%ctem .and. trim(descriptor%group) == "ctem") then
            validGroup = .true.
        elseif (config%fire .and. trim(descriptor%group) == "fire") then
            validGroup = .true.
        elseif (config%land .and. trim(descriptor%group) == "land") then
            validGroup = .true.
        elseif (config%methane .and. trim(descriptor%group) == "methane") then
            validGroup = .true.
        elseif (config%compete .and. trim(descriptor%group) == "compete") then
            validGroup = .true.
        else
            validGroup = .false.
        endif
    end function validGroup

    ! Determines wether the current variable matches the project configuration
    logical function validTime(config, timeFreq, descriptor)
        implicit none
        type(projectConfiguration), intent(in)  :: config
        type(outputDescriptor), intent(inout) :: descriptor
        character(*), intent(in)                :: timeFreq
        logical                                 :: valid

        valid = .true.
        if (config%annually .and. timeFreq == 'annually') then
            descriptor%timeFreq = timeFreq
        elseif (config%monthly .and. timeFreq == 'monthly') then
            descriptor%timeFreq = timeFreq
        elseif (config%daily .and. timeFreq == 'daily') then
            descriptor%timeFreq = timeFreq
        elseif (config%halfhourly .and. timeFreq == 'halfhourly') then
            descriptor%timeFreq = timeFreq
        else
            valid = .false.
        endif
        validTime = valid
    end function validTime

    ! Generates the filename for the current variable
    character(80) function generateFilename(outputForm, config, descriptor)
        implicit none
        type(projectConfiguration), intent(in)  :: config
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
        generateFilename = trim(config%path) // '/' // &
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
        character(30), intent(in)   :: key
        integer i
        do i=1, variableCount
            if (variables(i)%key == key) then
                getIdByKey = i
                return
            endif
        enddo
        getIdByKey = 0
    end function getIdByKey

    ! Prints out the project config
    subroutine printConfig(config)
        implicit none
        type(projectConfiguration), intent(in)  :: config
        print*, "Project configuration:"
        print*, "   halfhourly", config%halfhourly
        print*, "   daily", config%daily
        print*, "   monthly", config%monthly
        print*, "   annually", config%annually
    end subroutine printConfig

    ! Prints out a variable descriptor
    subroutine printDescriptor(descriptor)
        implicit none
        type(outputDescriptor), intent(in)    :: descriptor
        print*, "Variable descriptor:"
        print*, "   Group: ", descriptor%group
        print*, "   Short Name: ", descriptor%shortName
        print*, "   Standard Name: ", descriptor%standardName
        print*, "   Long Name: ", trim(descriptor%longName)
        print*, "   Units: ", descriptor%units
        print*, "   Time Frequency: ", descriptor%timeFreq
    end subroutine printDescriptor

    ! Prints out some information about the variables
    subroutine printVariables
        implicit none
        integer         :: i
        do i = 1, variableCount
            !print*, i, getIdByKey(variables(i)%key), variables(i)%ncid, variables(i)%key, trim(variables(i)%filename)
            print*, i, variables(i)%ncid, variables(i)%key, trim(variables(i)%filename)
        enddo
    end subroutine printVariables
end module outputManager
