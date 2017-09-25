module xmlManager

    use outputManager,   only : outputDescriptors, descriptorCount
    use xmlParser,   only : xml_process

    implicit none

    public :: loadoutputDescriptor
    public :: startfunc
    public :: datafunc
    public :: endfunc
    private :: charToLogical
    private :: charToInt

    character(len=80), dimension(2,10)      :: attribs
    character(len=400), dimension(100)      :: data
    character(len=80)                       :: currentGroup
    logical                                 :: error

contains

!> --------------
    subroutine loadoutputDescriptor()

        use ctem_statevars, only : c_switch

        implicit none

        character(:), pointer    :: xmlFile
        xmlFile         => c_switch%xmlFile

        call xml_process(xmlFile, attribs, data, startfunc, datafunc, endfunc, 0, error)

    end subroutine loadoutputDescriptor

!> --------------

    subroutine startfunc(tag, attribs, error)

        implicit none

        character(len=*)    :: tag, attribs(:,:)
        character(len=40)   :: attribute
        integer             :: id
        logical             :: error

        select case( tag )
            case( 'variableSet' )
                attribute = attribs(2,1)
                allocate(outputDescriptors(charToInt(attribute)))
            case( 'group' )
                currentGroup = trim(attribs(2,1))
            case( 'variable' )
                descriptorCount = descriptorCount + 1
                attribute = attribs(2,1)
                outputDescriptors(descriptorCount)%includeBareGround = charToLogical(attribute)
                outputDescriptors(descriptorCount)%group = currentGroup
        end select

    end subroutine

!> --------------

    subroutine datafunc(tag, data, error)

        implicit none

        character(len=*)               :: tag, data(:)
        character(len=400)             :: info
        integer                        :: i
        logical                        :: error

        info = trim(data(1))

        select case( tag )
            case( 'shortName' )
                outputDescriptors(descriptorCount)%shortName = info
            case( 'longName' )
                outputDescriptors(descriptorCount)%longName = info
            case( 'standardName' )
                outputDescriptors(descriptorCount)%standardName = info
            case( 'units' )
                outputDescriptors(descriptorCount)%units = info
        end select

    end subroutine

!> --------------

    subroutine endfunc(tag, error)
        implicit none
        character(len=*)               :: tag
        logical                        :: error
    end subroutine

!> --------------

    logical function charToLogical(input)
        implicit none
        character(len=*), intent(in)    :: input
        if (input == "true") then
            charToLogical = .true.
        else
            charToLogical = .false.
        endif
    end function charToLogical

!> --------------
    integer function charToInt(input)
        implicit none
        character(len=*), intent(in)    :: input
        read(input,*) charToInt
    end function charToInt

end module xmlManager
