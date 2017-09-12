#define dataPack type(dataPackage)
module dataTransferModule
    implicit none
    type dataPackage
        real, allocatable       :: values(:)
        integer                 :: dimensionCount
        integer, dimension(3)   :: dimension = 0
    end type
contains
    dataPack function subset(inData, from, to)
        implicit none
        dataPack, intent(in)    :: inData
        integer, intent(in)     :: from, to
        dataPack                :: newData
        allocate(newData%values(to - from + 1))
        newData%values = inData%values(from:to)
        subset = newData
    end function subset

    function inflateTo1D(p, dim)
        implicit none
        dataPack, intent(in)                    :: p
        integer, optional, intent(in)           :: dim
        real, dimension(:), allocatable         :: inflateTo1D
        allocate(inflateTo1D(size(p%values)))
        inflateTo1D = p%values
    end function inflateTo1D

    function inflateTo2D(p, dim)
        implicit none
        dataPack, intent(in)                    :: p
        integer, intent(in)                     :: dim(2)
        real, dimension(:,:), allocatable       :: inflateTo2D
        allocate(inflateTo2D(dim(1), dim(2)))
        inflateTo2D = reshape(p%values, dim)
    end function inflateTo2D

    function inflateTo3D(p, dim)
        implicit none
        dataPack, intent(in)                    :: p
        integer, intent(in)                     :: dim(3)
        real, dimension(:,:,:), allocatable     :: inflateTo3D
        allocate(inflateTo3D(dim(1), dim(2), dim(3)))
        inflateTo3D = reshape(p%values, dim)
    end function inflateTo3D

    function inflateTo4D(p, dim)
        implicit none
        dataPack, intent(in)                    :: p
        integer, intent(in)                     :: dim(4)
        real, dimension(:,:,:,:), allocatable   :: inflateTo4D
        allocate(inflateTo4D(dim(1), dim(2), dim(3), dim(4)))
        inflateTo4D = reshape(p%values, dim)
    end function inflateTo4D

    dataPack function deflateFrom1D(m)
        implicit none
        real, intent(in)    :: m(:)
        allocate(deflateFrom1D%values(size(m)))
        deflateFrom1D%values = m
    end function deflateFrom1D

    dataPack function deflateFrom2D(m)
        implicit none
        real, intent(in)    :: m(:,:)
        allocate(deflateFrom2D%values(size(m)))
        deflateFrom2D%values = reshape(m, (/size(m)/))
    end function deflateFrom2D

    dataPack function deflateFrom3D(m)
        implicit none
        real, intent(in)    :: m(:,:,:)
        allocate(deflateFrom3D%values(size(m)))
        deflateFrom3D%values = reshape(m, (/size(m)/))
    end function deflateFrom3D

    dataPack function deflateFrom4D(m)
        implicit none
        real, intent(in)    :: m(:,:,:,:)
        allocate(deflateFrom4D%values(size(m)))
        deflateFrom4D%values = reshape(m, (/size(m)/))
    end function deflateFrom4D

    dataPack function deflateFrom5D(m)
        implicit none
        real, intent(in)    :: m(:,:,:,:,:)
        allocate(deflateFrom5D%values(size(m)))
        deflateFrom5D%values = reshape(m, (/size(m)/))
    end function deflateFrom5D

    dataPack function generateIdentityPackage(count)
        implicit none
        integer, intent(in) :: count
        integer             :: i
        allocate(generateIdentityPackage%values(count))
        do i = 1, count
            generateIdentityPackage%values(i) = i
        enddo
    end function generateIdentityPackage
end module dataTransferModule
