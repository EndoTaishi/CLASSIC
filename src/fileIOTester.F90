program fileIOTester
#if PARALLEL
    use mpi
    use pnetcdf
#else
    use netcdf
#endif
    use fileIOModule
    implicit none
    integer             :: i, fileId, dimId, varId
    real                :: values(10)
    real, allocatable   :: temp(:)
#if PARALLEL
    call MPI_INIT(ierr)
#endif

    ! Create a new file
    fileId = ncCreate('test.nc', NF90_CLOBBER)

    ! Define a new dimension
    dimId = ncDefDim(fileId, 'testDimension', 0)

    ! Define a new variable
    varId = ncDefVar(fileId, 'testVariable', nf90_double, [dimId])

    ! Compute a unit array
    do i = 1, 10
        values(i) = i
    enddo

    ! Close define mode and open data mode
    call ncEndDef(fileId)

    ! Write some data to the file
    call ncPutVar(fileId, 'testVariable', realValues = values, start = [1], count = [10])

    ! Close the file
    call ncClose(fileId)

    ! Open existing file
    fileId = ncOpen('test.nc', NF90_NOWRITE)

    print*, "Use manual allocation on the server!"
    allocate(temp(5))
    temp = ncGet1DVar(fileId, 'testVariable', start = [6], count = [5])
    print*, temp
    deallocate(temp)

    ! Close the file
    call ncClose(fileId)

#if PARALLEL
    call MPI_FINALIZE(ierr)
#endif

end program fileIOTester
