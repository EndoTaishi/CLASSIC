!>\file
!! Principle driver program to run CLASSIC in stand-alone mode using specified boundary
!! conditions and atmospheric forcing. Depending upon the compiler options chosen, this
!! program will either run in MPI mode for use on parallel computing environments or
!! in serial mode for use in running at single sites.
!!
program CLASSIC

    ! Joe Melton and Ed Wisernig @ 2017

    use mpi
    use io_driver,              only : bounds,lonvect,latvect, &
                                       validCount,validLon,validLat
    use model_state_drivers,    only : read_modelsetup
    use netcdf_drivers,         only : create_out_netcdf
    use readjobopts,            only : read_from_job_options
    use main,                   only : main_driver
    use ctem_statevars,         only : alloc_ctem_vars
    use class_statevars,        only : alloc_class_vars
    use ctem_params,            only : nlat
    use supportFunctions,       only : isMainProcess

    implicit none

    double precision                :: time
    integer                         :: ierr, rank, size, i, cell, blocks, remainder

    ! MAIN PROGRAM

    ! Initialize the MPI and PnetCDF session
    call initializeParallelEnvironment

    ! Load the project config file. This parses the command line arguments.
    ! All model switches are read in from a namelist file. This sets up the
    ! run options and points to input files as needed.
    call read_from_job_options()

    ! Load the model setup information based on the metadata in the
    ! initialization netcdf file. The bounds given as an argument to
    ! CLASSIC are used to find the start points (srtx and srty)
    ! in the netcdf file, placing the gridcell on the domain of the
    ! input/output netcdfs. In read_modelsetup we use the netcdf to set
    ! the nmos, ignd,and ilg constants. It also opens the initial conditions
    ! file that is used below in read_initialstate.
    call read_modelsetup()

    ! Execute the following only on the main thread (see supportFunctions.f90)
    if (isMainProcess(rank)) then

        ! Generate the output files based on options in the joboptions file
        ! and the parameters of the initilization netcdf file.
        call create_out_netcdf

    endif

    ! Process the land grid cells, in parallel
    call processLandCells

    ! Close PnetCDF and shut down the MPI session
    call finalizeParallelEnvironment

    ! END MAIN PROGRAM

    !------------------

    contains

    subroutine processLandCells

        ! PROCESS LAND CELLS
        ! This section processes all of the land cells. There are validCount valid(i.e. land) cells, stored in validLon and validLat

        ! Since we know the nlat, nmos, ignd, and ilg we can allocate the CLASS and
        ! CTEM variable structures. This has to be done before call to main_driver.
        call alloc_class_vars()
        call alloc_ctem_vars()

        blocks = validCount / size + 1          ! The number of processing blocks
        remainder = mod(validCount, size)       ! The number of cells for the last block

        do i = 1, blocks - 1                    ! Go through every block except for the last one
            cell = (i - 1) * size + rank + 1
            call main_driver(validLon(cell),validLat(cell))
        enddo

        cell = (blocks - 1) * size + rank + 1   ! In the last block, process only the existing cells (NEEDS BETTER DESCRIPTION)
        if (rank < remainder) call main_driver(validLon(cell),validLat(cell))

    end subroutine processLandCells

    !------------------

    subroutine initializeParallelEnvironment

        ! INITIALIZE MPI and PnetCDf
        ! This section initializes the MPI environment and PnetCDF files

        call MPI_INIT(ierr)
        time = MPI_WTIME()
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
        ! HERE - PnetCDF init code comes here
        
    end subroutine initializeParallelEnvironment

    subroutine finalizeParallelEnvironment

        ! FINALIZE MPI and PnetCDF
        ! This section wraps up the whole MPI and PnetCDF megillah

        ! HERE - PnetCDF close session code comes here
        call MPI_FINALIZE(ierr)

    end subroutine finalizeParallelEnvironment

end program CLASSIC
