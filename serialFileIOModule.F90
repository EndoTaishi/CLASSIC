#define dataPack type(dataPackage)
module serialFileIOModule
    use netcdf
    use dataTransferModule
    implicit none
contains
    integer function ncOpen(fileName, omode)
        implicit none
        character(*), intent(in)    :: fileName
        integer, intent(in)         :: omode
        call check_nc(nf90_open(trim(fileName), omode, ncOpen))
    end function ncOpen

    integer function ncGetVarId(fileId, label)
        implicit none
        integer, intent(in)         :: fileId
        character(*), intent(in)    :: label
        call check_nc(nf90_inq_varid(fileId, label, ncGetVarId))
    end function ncGetVarId

    integer function ncGetVarDimensions(fileId, varId)
        implicit none
        integer, intent(in)         :: fileId, varId
        call check_nc(nf90_inquire_variable(fileId, varId, ndims = ncGetVarDimensions))
    end function ncGetVarDimensions

    dataPack function ncGetVar(fileId, label, start, count)
        implicit none
        integer, intent(in)                           :: fileId
        character(*), intent(in)                      :: label
        integer, dimension(:), intent(in)             :: start, count
        integer                                       :: varId, ndims, len
        integer, dimension(:), allocatable            :: dimIds(:), dimLens(:)
        real, dimension(:), allocatable               :: temp1D
        real, dimension(:,:), allocatable             :: temp2D
        real, dimension(:,:,:), allocatable           :: temp3D
        real, dimension(:,:,:,:), allocatable         :: temp4D
        real, dimension(:,:,:,:,:), allocatable       :: temp5D

        varId = ncGetVarId(fileId, label)
        ndims = ncGetVarDimensions(fileId, varId)

        select case(ndims)
        case(1)
            allocate(ncGetVar%values(count(1)))
            allocate(temp1D(count(1)))
            call check_nc(nf90_get_var(fileId, varId, temp1D, start = start, count = count))
            ncGetVar = deflateFrom1D(temp1D)
        case(2)
            allocate(ncGetVar%values(count(1) * count(2)))
            allocate(temp2D(count(1), count(2)))
            call check_nc(nf90_get_var(fileId, varId, temp2D, start = start, count = count))
            ncGetVar = deflateFrom2D(temp2D)
        case(3)
            allocate(ncGetVar%values(count(1) * count(2) * count(3)))
            allocate(temp3D(count(1), count(2), count(3)))
            call check_nc(nf90_get_var(fileId, varId, temp3D, start = start, count = count))
            ncGetVar = deflateFrom3D(temp3D)
        case(4)
            allocate(ncGetVar%values(count(1) * count(2) * count(3) * count(4)))
            allocate(temp4D(count(1), count(2), count(3), count(4)))
            call check_nc(nf90_get_var(fileId, varId, temp4D, start = start, count = count))
            ncGetVar = deflateFrom4D(temp4D)
        case(5)
            allocate(ncGetVar%values(count(1) * count(2) * count(3) * count(4) * count(5)))
            allocate(temp5D(count(1), count(2), count(3), count(4), count(5)))
            call check_nc(nf90_get_var(fileId, varId, temp5D, start = start, count = count))
            ncGetVar = deflateFrom5D(temp5D)
        case default
            stop("Only up to 5 dimensions have been implemented!")
        end select
    end function ncGetVar

    subroutine ncPut1DVar(fileId, label, data, start, count)
        implicit none
        integer, intent(in)                                     :: fileId
        character(*), intent(in)                                :: label
        real, dimension(:), intent(in)                          :: data
        integer, dimension(:), intent(in)                       :: start
        integer, dimension(:), optional, intent(in)             :: count
        integer                                                 :: varId
        real, dimension(:,:), allocatable                       :: temp2D
        integer, dimension(2)                                   :: localFormat, localCount
        if (present(count)) then
            localCount = count
        else
            localCount = [1, 1]
        endif
        localFormat = localCount
        varId = ncGetVarId(fileId, label)
        temp2D = reshape(data, localFormat)
        call check_nc(nf90_put_var(fileId, varId, temp2D, start, localCount))
    end subroutine ncPut1DVar

    subroutine ncPut2DVar(fileId, label, data, start, count)
        implicit none
        integer, intent(in)                                     :: fileId
        character(*), intent(in)                                :: label
        real, dimension(:,:), intent(in)                        :: data
        integer, dimension(:), intent(in)                       :: start
        integer, dimension(:), optional, intent(in)             :: count
        integer                                                 :: varId
        real, dimension(:,:,:), allocatable                     :: temp3D
        integer, dimension(3)                                   :: localFormat, localCount
        if (present(count)) then
            localCount = count
        else
            localCount = [1, 1, 1]
        endif
        localFormat = count
        varId = ncGetVarId(fileId, label)
        temp3D = reshape(data, localFormat)
        call check_nc(nf90_put_var(fileId, varId, temp3D, start, localCount))
    end subroutine ncPut2DVar

    subroutine ncPut3DVar(fileId, label, data, start, count)
        implicit none
        integer, intent(in)                                     :: fileId
        character(*), intent(in)                                :: label
        real, dimension(:,:,:), intent(in)                      :: data
        integer, dimension(:), intent(in)                       :: start, count
        integer                                                 :: varId
        real, dimension(:,:,:,:), allocatable                   :: temp4D
        integer, dimension(4)                                   :: format = 0
        format = count
        varId = ncGetVarId(fileId, label)
        temp4D = reshape(data, format)
        call check_nc(nf90_put_var(fileId, varId, temp4D, start, count))
    end subroutine ncPut3DVar

    subroutine ncPut4DVar(fileId, label, data, start, count)
        implicit none
        integer, intent(in)                                     :: fileId
        character(*), intent(in)                                :: label
        real, dimension(:,:,:,:), intent(in)                      :: data
        integer, dimension(:), intent(in)                       :: start, count
        integer                                                 :: varId
        real, dimension(:,:,:,:,:), allocatable                   :: temp5D
        integer, dimension(5)                                   :: format = 0
        format = count
        varId = ncGetVarId(fileId, label)
        temp5D = reshape(data, format)
        call check_nc(nf90_put_var(fileId, varId, temp5D, start, count))
    end subroutine ncPut4DVar

    subroutine ncClose(fileId)
        implicit none
        integer, intent(in)         :: fileId
        call check_nc(nf90_close(fileId))
    end subroutine ncClose

    subroutine check_nc(nc_status)
        use netcdf
        implicit none
        integer, intent(in) :: nc_status
        if(nc_status /= nf90_noerr) then
            write(0,*)'netCDF error: ',trim(nf90_strerror(nc_status))
            stop
        end if
    end subroutine check_nc

    function ncGetDimValues(fileId, label, start, count)
        implicit none
        integer, intent(in)                         :: fileId
        character(*), intent(in)                    :: label
        integer, dimension(:), intent(in)           :: count
        integer, dimension(:), optional, intent(in) :: start
        dataPack                                    :: data
        real, dimension(:), allocatable             :: ncGetDimValues
        integer, dimension(1)                       :: localStart
        if (present(start)) then
            localStart = start
        else
            localStart = [1]
        endif

        data = ncGetVar(fileId, label, localStart, count)
        ncGetDimValues = inflateTo1D(data)
    end function ncGetDimValues

    function ncGet1DVar(fileId, label, start, count, format)
        implicit none
        integer, intent(in)                         :: fileId
        character(*), intent(in)                    :: label
        integer, dimension(:), optional, intent(in) :: start, count, format
        dataPack                                    :: data
        real, dimension(:), allocatable             :: ncGet1DVar
        integer, dimension(2)                       :: localCount
        if (present(count)) then
            localCount = count
        else
            localCount = [1, 1]
        endif
        data = ncGetVar(fileId, label, start, localCount)
        ncGet1DVar = inflateTo1D(data)
    end function ncGet1DVar

    function ncGet2DVar(fileId, label, start, count, format)
        implicit none
        integer, intent(in)                         :: fileId
        character(*), intent(in)                    :: label
        integer, dimension(:), intent(in)           :: start
        integer, dimension(:), optional, intent(in) :: count, format
        dataPack                                    :: data
        real, dimension(:,:), allocatable           :: ncGet2DVar
        integer, dimension(3)                       :: localCount
        integer, dimension(2)                       :: localFormat
        if (present(format)) then
            localFormat = format
        else
            localFormat = [1, 1]
        endif
        if (present(count)) then
            localCount = count
        else
            localCount = [1, 1, 1]
        endif
        data = ncGetVar(fileId, label, start, localCount)
        ncGet2DVar = inflateTo2D(data, localFormat)
    end function ncGet2DVar

    function ncGet3DVar(fileId, label, start, count, format)
        implicit none
        integer, intent(in)                             :: fileId
        character(*), intent(in)                        :: label
        integer, dimension(:), optional, intent(in)     :: start, count, format
        dataPack                                        :: data
        real, dimension(:,:,:), allocatable             :: ncGet3DVar
        integer, dimension(4)                           :: localCount
        integer, dimension(3)                           :: localFormat
        if (present(format)) then
            localFormat = format
        else
            localFormat = [1, 1, 1]
        endif
        if (present(count)) then
            localCount = count
        else
            localCount = [1, 1, 1, 1]
        endif
        data = ncGetVar(fileId, label, start, localCount)
        ncGet3DVar = inflateTo3D(data, localFormat)
    end function ncGet3DVar

    function ncGet4DVar(fileId, label, start, count, format)
        implicit none
        integer, intent(in)                             :: fileId
        character(*), intent(in)                        :: label
        integer, dimension(:), optional, intent(in)     :: start, count, format
        dataPack                                        :: data
        real, dimension(:,:,:,:), allocatable           :: ncGet4DVar
        integer, dimension(4)                           :: localCount, localFormat
        if (present(format)) then
            localFormat = format
        else
            localFormat = [1, 1, 1, 1]
        endif
        if (present(count)) then
            localCount = count
        else
            localCount = [1, 1, 1, 1]
        endif
        data = ncGetVar(fileId, label, start, localCount)
        ncGet4DVar = inflateTo4D(data, localFormat)
    end function ncGet4DVar

    integer function ncGetDimId(fileId, label)
        implicit none
        integer, intent(in)         :: fileId
        character(*), intent(in)    :: label
        call check_nc(nf90_inq_dimid(fileId, label, ncGetDimId))
    end function ncGetDimId

    integer function ncGetDimLen(fileId, label)
        implicit none
        integer, intent(in)             :: fileId
        character(*), intent(in)        :: label
        integer                         :: dimId
        integer                         :: dimLen
        dimId = ncGetDimId(fileId, label)
        call check_nc(nf90_inquire_dimension(fileId, dimId, len=dimLen))
        ncGetDimLen = dimLen
    end function ncGetDimLen

end module serialFileIOModule
