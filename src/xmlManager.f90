!> The XML Manager Module handles the loading of variable descriptors and variable variants from an XML file.
!! Once the variables are generated, corresponding NetCDF output files are created.
module xmlManager
    use outputManager,  only : outputDescriptor, outputDescriptors, descriptorCount, &
                                variant, variants, variantCount
    use xmlParser,      only : xml_process
    implicit none
    public :: loadoutputDescriptor
    public :: startfunc
    public :: datafunc
    public :: endfunc
    private :: charToLogical
    private :: charToInt

    character(len=80), dimension(2,10)      :: attribs
    character(len=400), dimension(100)      :: data
    logical                                 :: error
    character(len=80)                       :: currentGroup, currentVariableName, variableSetType, variableSetDate, variableSetVersion
    real                                    :: xmlVersion

contains
    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    !> The [change my name] function parses the XML file using the startfunc, datafunc and endfunc for starting tags, data tags and end tags, respectively.
    subroutine loadoutputDescriptor()
        use ctem_statevars, only : c_switch
        character(:), pointer    :: xmlFile
        xmlFile         => c_switch%xmlFile
        call xml_process(xmlFile, attribs, data, startfunc, datafunc, endfunc, 0, error)
    end subroutine loadoutputDescriptor

    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    !> The startfunc function parses opening tags or start tags. It gets triggered for each and every opening tag of the XML document.
    subroutine startfunc(tag, attribs, error)
        character(len=*)                    :: tag, attribs(:,:)
        character(len=40)                   :: attribute
        integer                             :: id
        logical                             :: error
        type(outputDescriptor), allocatable :: tempDescriptors(:)
        type(variant), allocatable          :: tempVariants(:)

        ! Examine the tag
        select case( tag )
            ! If the tag is <variableSet>, allocate the variable descriptors and the variants.
            case( 'variableSet' )
                xmlVersion = 0;
                variableSetType = trim(attribs(2,1))
                variableSetVersion = trim(attribs(2,2));
                if (variableSetVersion == '') then
                    stop("The input XML document doesn't feature the required version field of the <variableSet> node");
                else
                    xmlVersion = charToReal(variableSetVersion);
                    if (xmlVersion < 1.1) stop('Older XML document found, please upgrade to a more recent version');
                endif

                variableSetDate = trim(attribs(2,3))
                allocate(outputDescriptors(0))
                allocate(variants(0))
            ! If the tag is <group>, remember the current group.
            case( 'group' )
                currentGroup = trim(attribs(2,1))
            ! If the tag is <variable>, increment the descriptor count and allocate the temporary descriptors.
            ! Then, add a new descriptor to the array, by copying and extending the array, followed by setting some values.
            case( 'variable' )
                descriptorCount = descriptorCount + 1
                allocate(tempDescriptors(descriptorCount))
                tempDescriptors(1 : descriptorCount - 1) = outputDescriptors(1 : descriptorCount - 1)
                call move_alloc(tempDescriptors, outputDescriptors)
                attribute = attribs(2,1)
                outputDescriptors(descriptorCount)%includeBareGround = charToLogical(attribute)
                outputDescriptors(descriptorCount)%group = currentGroup
            ! If the tag is <variant>, then similar to the above section, extend the array and append a new variant.
            case('variant')
                variantCount = variantCount + 1
                allocate(tempVariants(variantCount))
                tempVariants(1 : variantCount - 1) = variants(1 : variantCount - 1)
                call move_alloc(tempVariants, variants)
                variants(variantCount)%shortName = currentVariableName
        end select
    end subroutine

    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    !> The datafunc function parses the content of a tag.
    subroutine datafunc(tag, data, error)
        character(len=*)               :: tag, data(:)
        character(len=400)             :: info
        integer                        :: i
        logical                        :: error

        info = trim(data(1))
        ! Examine the tag and store the tag content in the descriptors array.
        select case( tag )
            case( 'shortName' )
                outputDescriptors(descriptorCount)%shortName = info
                currentVariableName = info
            case( 'longName' )
                outputDescriptors(descriptorCount)%longName = info
            case( 'standardName' )
                outputDescriptors(descriptorCount)%standardName = info
            case( 'units' )
                !outputDescriptors(descriptorCount)%units = info
                variants(variantCount)%units = info
            case( 'nameInCode' )
                variants(variantCount)%nameInCode = info
            case( 'timeFrequency' )
                variants(variantCount)%timeFrequency = info
            case( 'outputForm' )
                variants(variantCount)%outputForm = info
        end select
    end subroutine

    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    !> The endfunc function gets triggered when the closing XML element is encountered. E.g. if we wanted to deallocate at </variableSet>.
    subroutine endfunc(tag, error)
        character(len=*)               :: tag
        logical                        :: error
        ! Place a select case(tag) in here if you need to do anything on a closing tag.
    end subroutine

    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    !> The charToLogical function returns the logical value of a given char input
    logical function charToLogical(input)
        character(len=*), intent(in)    :: input    !< Char input
        if (input == "true") then
            charToLogical = .true.
        else
            charToLogical = .false.
        endif
    end function charToLogical

    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    !> The charToInt function returns the integer value of a given char input
    integer function charToInt(input)
        character(len=*), intent(in)    :: input    !< Char input
        read(input,*) charToInt
    end function charToInt

    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    !> The charToInt function returns the integer value of a given char input
    real function charToReal(input)
        character(len=*), intent(in)    :: input    !< Char input
        read(input,*) charToReal
    end function charToReal
end module xmlManager
