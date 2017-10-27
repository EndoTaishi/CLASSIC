!> The Data Transfer Module contains methods for transferring data between functions, by serializing them into 1D real arrays
module dataTransferModule
    implicit none
    type dataPackage
        real, allocatable       :: values(:)        !< the actual data stored
    end type
contains
    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    !> Retrieves a 1D real array out of a data package
    function inflateTo1D(p, dim)

        type(dataPackage), intent(in)           :: p            !< Data package
        integer, optional, intent(in)           :: dim          !< dim is not used, but left optional for compatibility with inflateTo2D and higher dimensions
        real, dimension(:), allocatable         :: inflateTo1D  !< function type
        allocate(inflateTo1D(size(p%values)))
        inflateTo1D = p%values
    end function inflateTo1D

    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    !> Retrieves a 2D real array from a data package
    function inflateTo2D(p, dim)

        type(dataPackage), intent(in)           :: p            !< data package variable
        integer, intent(in)                     :: dim(2)       !< 2 dimensions
        real, dimension(:,:), allocatable       :: inflateTo2D  !< function type declaration
        allocate(inflateTo2D(dim(1), dim(2)))
        inflateTo2D = reshape(p%values, dim)
    end function inflateTo2D

    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    !> Retrieves a 3D real array from a data package
    function inflateTo3D(p, dim)

        type(dataPackage), intent(in)           :: p            !< data package variable
        integer, intent(in)                     :: dim(3)       !< 3 dimensions
        real, dimension(:,:,:), allocatable     :: inflateTo3D  !< function type declaration
        allocate(inflateTo3D(dim(1), dim(2), dim(3)))
        inflateTo3D = reshape(p%values, dim)
    end function inflateTo3D

    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    !> Retrieves a 4D real array from a data package
    function inflateTo4D(p, dim)

        type(dataPackage), intent(in)           :: p            !< data package variable
        integer, intent(in)                     :: dim(4)       !< 4 dimensions
        real, dimension(:,:,:,:), allocatable   :: inflateTo4D  !< function type declaration
        allocate(inflateTo4D(dim(1), dim(2), dim(3), dim(4)))
        inflateTo4D = reshape(p%values, dim)
    end function inflateTo4D

    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    !> Retrieves a data package from a 1D real array
    type(dataPackage) function deflateFrom1D(m)

        real, intent(in)    :: m(:)                             !< input array
        allocate(deflateFrom1D%values(size(m)))
        deflateFrom1D%values = m
    end function deflateFrom1D

    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    !> Retrieves a data package from a 2D real array
    type(dataPackage) function deflateFrom2D(m)

        real, intent(in)    :: m(:,:)                           !< input array
        allocate(deflateFrom2D%values(size(m)))
        deflateFrom2D%values = reshape(m, (/size(m)/))
    end function deflateFrom2D

    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    !> Retrieves a data package from a 3D real array
    type(dataPackage) function deflateFrom3D(m)

        real, intent(in)    :: m(:,:,:)                         !< input array
        allocate(deflateFrom3D%values(size(m)))
        deflateFrom3D%values = reshape(m, (/size(m)/))
    end function deflateFrom3D

    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    !> Retrieves a data package from a 4D real array
    type(dataPackage) function deflateFrom4D(m)

        real, intent(in)    :: m(:,:,:,:)                       !< input array
        allocate(deflateFrom4D%values(size(m)))
        deflateFrom4D%values = reshape(m, (/size(m)/))
    end function deflateFrom4D

    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    !> Retrieves a data package from a 5D real array
    type(dataPackage) function deflateFrom5D(m)

        real, intent(in)    :: m(:,:,:,:,:)                     !< input array
        allocate(deflateFrom5D%values(size(m)))
        deflateFrom5D%values = reshape(m, (/size(m)/))
    end function deflateFrom5D

    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    !> Retrieves a data package variable with all the elements numbered from 1 to count
    type(dataPackage) function generateIdentityPackage(count)

        integer, intent(in) :: count                            !< count of elements to generate
        integer             :: i
        allocate(generateIdentityPackage%values(count))
        do i = 1, count
            generateIdentityPackage%values(i) = i
        enddo
    end function generateIdentityPackage
end module dataTransferModule
