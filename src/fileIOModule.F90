#define dataPack type(dataPackage)
module fileIOModule
#if PARALLEL
    use mpi
    use pnetcdf
#else
    use netcdf
#endif
    use dataTransferModule
    implicit none
contains

! Create new file
    integer function ncCreate(fileName, cmode)
        character(*), intent(in)    :: fileName
        integer, intent(in)         :: cmode
#if PARALLEL
        ! we assume MPI_COMM_WORLD and MPI_INFO_NULL are common
        call check_nc(nf90mpi_create(MPI_COMM_WORLD, trim(fileName), cmode, MPI_INFO_NULL, ncCreate))
#else
        call check_nc(nf90_create(trim(fileName), cmode, ncCreate))
#endif
    end function ncCreate

! Open existing file
    integer function ncOpen(fileName, omode)
        character(*), intent(in)    :: fileName
        integer, intent(in)         :: omode
#if PARALLEL
        ! we assume MPI_COMM_WORLD and MPI_INFO_NULL are common
        call check_nc(nf90mpi_open(MPI_COMM_WORLD, trim(fileName), omode, MPI_INFO_NULL, ncOpen))
#else
        call check_nc(nf90_open(trim(fileName), omode, ncOpen))
#endif
    end function ncOpen

! Get variable ID
    integer function ncGetVarId(fileId, label)
        integer, intent(in)         :: fileId
        character(*), intent(in)    :: label
#if PARALLEL
        call check_nc(nf90mpi_inq_varid(fileId, label, ncGetVarId))
#else
        call check_nc(nf90_inq_varid(fileId, label, ncGetVarId))
#endif
    end function ncGetVarId


! Get variable dimensions
    integer function ncGetVarDimensions(fileId, varId)
        integer, intent(in)         :: fileId, varId
#if PARALLEL
        call check_nc(nf90mpi_inquire_variable(fileId, varId, ndims = ncGetVarDimensions))
#else
        call check_nc(nf90_inquire_variable(fileId, varId, ndims = ncGetVarDimensions))
#endif
    end function ncGetVarDimensions


! Get dimension ID
    integer function ncGetDimId(fileId, label)
        integer, intent(in)         :: fileId
        character(*), intent(in)    :: label
#if PARALLEL
        call check_nc(nf90mpi_inq_dimid(fileId, label, ncGetDimId))
#else
        call check_nc(nf90_inq_dimid(fileId, label, ncGetDimId))
#endif
    end function ncGetDimId


! Get dimension length
    integer function ncGetDimLen(fileId, label)
        integer, intent(in)             :: fileId
        character(*), intent(in)        :: label
        integer                         :: dimId
#if PARALLEL
        integer(kind=MPI_OFFSET_KIND)   :: dimLen
#else
        integer                         :: dimLen
#endif
        dimId = ncGetDimId(fileId, label)
#if PARALLEL
        call check_nc(nf90mpi_inquire_dimension(fileId, dimId, len=dimLen))
#else
        call check_nc(nf90_inquire_dimension(fileId, dimId, len=dimLen))
#endif
        ncGetDimLen = dimLen
    end function ncGetDimLen

! Define dimension
    integer function ncDefDim(fileId, label, length)
        integer, intent(in)                         :: fileId
        integer, intent(in)                         :: length
        character(*), intent(in)                    :: label
#if PARALLEL
        call check_nc(nf90mpi_def_dim(fileId, label, int(length,8), ncDefDim))
#else
        call check_nc(nf90_def_dim(fileId, label, length, ncDefDim))
#endif
    end function ncDefDim

! Define variable
    integer function ncDefVar(fileId, label, type, dimIds)
        integer, intent(in)                         :: fileId, dimIds(:), type
        character(*), intent(in)                    :: label
#if PARALLEL
        call check_nc(nf90mpi_def_var(fileId, label, type, dimIds, ncDefVar))
#else
        call check_nc(nf90_def_var(fileId, label, type, dimIds, ncDefVar))
#endif
    end function ncDefVar

! End definition
    subroutine ncEndDef(fileId)
        integer, intent(in)                 :: fileId
#if PARALLEL
        call check_nc(nf90mpi_enddef(fileId))
#else
        call check_nc(nf90_enddef(fileId))
#endif
    end subroutine ncEndDef

! Redef
    subroutine ncReDef(fileId)
        integer, intent(in)                 :: fileId
#if PARALLEL
        call check_nc(nf90mpi_redef(fileId))
#else
        call check_nc(nf90_redef(fileId))
#endif
    end subroutine ncReDef

! Put attribute
    subroutine ncPutAtt(fileId, varId, label, charvalues, intvalues,realvalues)
        integer, intent(in)     :: fileId, varId
        character(*)            :: label
        character(*), optional  :: charvalues
        integer, optional       :: intvalues
        real, optional          :: realvalues
#if PARALLEL
        if (present(charvalues)) then
            call check_nc(nf90mpi_put_att(fileId, varId, label, charvalues))
        else if (present(intvalues)) then
            call check_nc(nf90mpi_put_att(fileId, varId, label, intvalues))
        else
            call check_nc(nf90mpi_put_att(fileId, varId, label, realvalues))
        end if
#else
        if (present(charvalues)) then
            call check_nc(nf90_put_att(fileId, varId, label, charvalues))
        else if (present(intvalues)) then
            call check_nc(nf90_put_att(fileId, varId, label, intvalues))
        else
            call check_nc(nf90_put_att(fileId, varId, label, realvalues))
        end if
#endif
    end subroutine ncPutAtt

! Close current file
    subroutine ncClose(fileId)
        integer, intent(in)                 :: fileId
#if PARALLEL
        call check_nc(nf90mpi_close(fileId))
#else
        call check_nc(nf90_close(fileId))
#endif
    end subroutine ncClose

! Check for errors in the NetCDF data retrieval process
    subroutine check_nc(nc_status)
        integer, intent(in) :: nc_status
        integer             :: err
        if(nc_status /= nf90_noerr) then
#if PARALLEL
            write(0,*)'netCDF error: ',trim(nf90mpi_strerror(nc_status))
            call MPI_Abort(MPI_COMM_WORLD, -1, err)
#else
            write(0,*)'netCDF error: ',trim(nf90_strerror(nc_status))
            stop
#endif
        end if
    end subroutine check_nc

! Get variable content in the form of a dataPack
#if PARALLEL
    dataPack function ncGetVar(fileId, label, start, count)
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
            call check_nc(nf90mpi_get_var_all(fileId, varId, temp1D, start = int(start,8), count = int(count,8)))
            ncGetVar = deflateFrom1D(temp1D)
        case(2)
            allocate(ncGetVar%values(count(1) * count(2)))
            allocate(temp2D(count(1), count(2)))
            call check_nc(nf90mpi_get_var_all(fileId, varId, temp2D, start = int(start,8), count = int(count,8)))
            ncGetVar = deflateFrom2D(temp2D)
        case(3)
            allocate(ncGetVar%values(count(1) * count(2) * count(3)))
            allocate(temp3D(count(1), count(2), count(3)))
            call check_nc(nf90mpi_get_var_all(fileId, varId, temp3D, start = int(start,8), count = int(count,8)))
            ncGetVar = deflateFrom3D(temp3D)
        case(4)
            allocate(ncGetVar%values(count(1) * count(2) * count(3) * count(4)))
            allocate(temp4D(count(1), count(2), count(3), count(4)))
            call check_nc(nf90mpi_get_var_all(fileId, varId, temp4D, start = int(start,8), count = int(count,8)))
            ncGetVar = deflateFrom4D(temp4D)
        case(5)
            allocate(ncGetVar%values(count(1) * count(2) * count(3) * count(4) * count(5)))
            allocate(temp5D(count(1), count(2), count(3), count(4), count(5)))
            call check_nc(nf90mpi_get_var_all(fileId, varId, temp5D, start = int(start,8), count = int(count,8)))
            ncGetVar = deflateFrom5D(temp5D)
        case default
            stop("Only up to 5 dimensions have been implemented!")
        end select
    end function ncGetVar
#else
    dataPack function ncGetVar(fileId, label, start, count)
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
#endif

! Write the 1D dimension values to the NetCDF file
#if PARALLEL
    subroutine ncPutDimValues(fileId, label, data, intdata, start, count)
        integer, intent(in)                                     :: fileId
        character(*), intent(in)                                :: label
        real, dimension(:), intent(inout)                       :: data
        integer, optional, intent(in)                           :: start
        integer, intent(in)                                     :: count
        integer, dimension(:), optional, intent(inout)          :: intdata
        integer                                                 :: varId,localStart
        if (present(start)) then
            localStart = start
        else
            localStart = 1
        endif
        varId = ncGetVarId(fileId, label)
        if (present(intdata)) then  !if the data is integer
            call check_nc(nf90mpi_put_var_all(fileId, varId, intdata, int((/localStart/),8), int((/count/),8)))
        else
            call check_nc(nf90mpi_put_var_all(fileId, varId, data, int((/localStart/),8), int((/count/),8)))
        end if

    end subroutine ncPutDimValues
#else
    subroutine ncPutDimValues(fileId, label, data, intdata, start, count)
        integer, intent(in)                                     :: fileId
        character(*), intent(in)                                :: label
        real, dimension(:), intent(in)                          :: data
        integer, optional, intent(in)                           :: start
        integer, intent(in)                                     :: count
        integer, dimension(:), optional, intent(in)             :: intdata
        integer                                                 :: varId,localStart
        if (present(start)) then
            localStart = start
        else
            localStart = 1
        endif
        varId = ncGetVarId(fileId, label)
        if (present(intdata)) then  !if the data is integer
            call check_nc(nf90_put_var(fileId, varId, intdata, (/localStart/),(/count/)))
        else
            call check_nc(nf90_put_var(fileId, varId, data, (/localStart/), (/count/)))
        end if

    end subroutine ncPutDimValues
#endif

! Write a local 1D variable (2D in NetCDF file)

    subroutine ncPut1DVar(fileId, label, data, intdata, start, count)
        integer, intent(in)                                     :: fileId
        character(*), intent(in)                                :: label
        real, dimension(:), intent(inout)                       :: data
        integer, dimension(:), intent(in)                       :: start
        integer, dimension(:), optional, intent(in)             :: count
        integer, dimension(:), optional, intent(inout)          :: intdata
        integer                                                 :: varId
        integer, dimension(:), allocatable                      :: localFormat, localCount
        if (present(count)) then
            allocate(localCount(size(count)))
            localCount = count
        else
            allocate(localCount(2))
            localCount = [1, 1]
        endif
        localFormat = localCount
        varId = ncGetVarId(fileId, label)
        if (present(intdata)) then  !if the data is integer
#if PARALLEL
            call check_nc(nf90mpi_put_var_all(fileId, varId, intdata, int(start,8), int(localcount,8)))
#else
            call check_nc(nf90_put_var(fileId, varId, intdata, start, localCount))
#endif
        else
#if PARALLEL
            call check_nc(nf90mpi_put_var_all(fileId, varId, data, int(start,8), int(localcount,8)))
#else
            call check_nc(nf90_put_var(fileId, varId, data, start, localCount))
#endif
        end if

    end subroutine ncPut1DVar


! Write a local 2D variable (3D in NetCDF file)
#if PARALLEL
    subroutine ncPut2DVar(fileId, label, data, start, count)
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
        call check_nc(nf90mpi_put_var_all(fileId, varId, temp3D, int(start,8), int(localCount,8)))
    end subroutine ncPut2DVar
#else
    subroutine ncPut2DVar(fileId, label, data, start, count)
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
#endif

! Write a local 3D variable (4D in NetCDF file)
#if PARALLEL
    subroutine ncPut3DVar(fileId, label, data, start, count)
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
        call check_nc(nf90mpi_put_var_all(fileId, varId, temp4D, int(start,8), int(count,8)))
    end subroutine ncPut3DVar
#else
    subroutine ncPut3DVar(fileId, label, data, start, count)
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
#endif

    ! Get the values stored in a dimension (e.g. get Lon, Lat or Time values)
    function ncGetDimValues(fileId, label, start, count)
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

!     function ncGet1DVar(fileId, label, start, count, format)  !FLAG - Ed looks we have some stuff optional that shouldn't be. Pls look over.
!         integer, intent(in)                         :: fileId
!         character(*), intent(in)                    :: label
!         integer, dimension(:), optional, intent(in) :: start, count, format
!         dataPack                                    :: data
!         real, dimension(:), allocatable             :: ncGet1DVar
!         integer, dimension(2)                       :: localCount
!         if (present(count)) then
!             localCount = count
!         else
!             localCount = [1, 1]
!         endif
!         data = ncGetVar(fileId, label, start, localCount)
!         ncGet1DVar = inflateTo1D(data)
!     end function ncGet1DVar

    function ncGet1DVar(fileId, label, start, count, format)  !FLAG - Ed looks we have some stuff optional that shouldn't be. Pls look over.
        integer, intent(in)                         :: fileId
        character(*), intent(in)                    :: label
        integer, dimension(:), optional, intent(in) :: start, count, format
        dataPack                                    :: data
        real, dimension(:), allocatable             :: ncGet1DVar
        integer, allocatable                        :: localCount(:)
        if (present(count)) then
            allocate(localCount(size(count)))
            localCount = count
        else
            allocate(localCount(2))
            localCount = [1, 1]
        endif
        data = ncGetVar(fileId, label, start, localCount)
        ncGet1DVar = inflateTo1D(data)
    end function ncGet1DVar

!     function ncGet2DVar(fileId, label, start, count, format)
!         integer, intent(in)                         :: fileId
!         character(*), intent(in)                    :: label
!         integer, dimension(:), intent(in)           :: start
!         integer, dimension(:), optional, intent(in) :: count, format
!         dataPack                                    :: data
!         real, dimension(:,:), allocatable           :: ncGet2DVar
!         integer, dimension(3)                       :: localCount
!         integer, dimension(2)                       :: localFormat
!         if (present(format)) then
!             localFormat = format
!         else
!             localFormat = [1, 1]
!         endif
!         if (present(count)) then
!             if (size(count) == 3) then
!                 localCount = count
!                 data = ncGetVar(fileId, label, start, localCount)
!                 ncGet2DVar = inflateTo2D(data, localFormat)
!             else ! this can be called to just get a simple 2D field so do that here.
!                 data = ncGetVar(fileId, label, start, count)
!                 ncGet2DVar = inflateTo2D(data, count)
!             end if
!         else
!             localCount = [1, 1, 1]
!             data = ncGetVar(fileId, label, start, localCount)
!             ncGet2DVar = inflateTo2D(data, localFormat)
!         endif
!     end function ncGet2DVar

    function ncGet2DVar(fileId, label, start, count, format)
        integer, intent(in)                         :: fileId
        character(*), intent(in)                    :: label
        integer, intent(in)                         :: start(:)
        integer, intent(in), optional               :: count(:), format(:)
        real, allocatable                           :: ncGet2DVar(:,:)
        integer, allocatable                        :: localCount(:), localFormat(:)
        dataPack                                    :: data

        if (present(format)) then
            allocate(localFormat(size(format)))
            localFormat = format
        else
            allocate(localFormat(3))
            localFormat = [1, 1, 1]
        endif

        if (present(count)) then
            allocate(localCount(size(count)))
            localCount = count
        else
            allocate(localCount(3))
            localCount = [1, 1, 1]
        endif

        data = ncGetVar(fileId, label, start, localCount)
        ncGet2DVar = inflateTo2D(data, localFormat)

    end function ncGet2DVar

    function ncGet3DVar(fileId, label, start, count, format)
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
        integer, intent(in)                             :: fileId
        character(*), intent(in)                        :: label
        integer, dimension(:), optional, intent(in)     :: start, count, format
        dataPack                                        :: data
        real, dimension(:,:,:), allocatable             :: ncGet4DVar
        integer, dimension(5)                           :: localCount
        integer, dimension(4)                           :: localFormat
        if (present(format)) then
            localFormat = format
        else
            localFormat = [1, 1, 1, 1]
        endif
        if (present(count)) then
            localCount = count
        else
            localCount = [1, 1, 1, 1, 1]
        endif
        data = ncGetVar(fileId, label, start, localCount)
        ncGet4DVar = inflateTo3D(data, localFormat)
    end function ncGet4DVar

end module fileIOModule
