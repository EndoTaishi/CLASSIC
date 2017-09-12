!>\file
!! Principle driver program to run CLASSIC in stand-alone mode using specified boundary
!! conditions and atmospheric forcing. Depending upon the compiler options chosen, this
!! program will either run in MPI mode for use on parallel computing environments or
!! in serial mode for use in running at single sites.
!!
program CLASSIC

    ! Joe Melton and Ed Wisernig @ 2017

    use mpi
    use io_driver,              only : bounds,validCount,validLon,validLat, &
                                       validLonIndex,validLatIndex
    use model_state_drivers,    only : read_modelsetup
    !use netcdf_drivers,         only : create_out_netcdf
    use xmlManager,             only : loadoutputDescriptor
    !use outputManager,          only : generateOutputVariables
    use readjobopts,            only : read_from_job_options
    use main,                   only : main_driver
    use ctem_statevars,         only : alloc_ctem_vars
    use class_statevars,        only : alloc_class_vars

    implicit none

    double precision                :: time
    integer                         :: ierr, rank, size, i, cell, blocks, remainder
    logical                         :: ready_to_go = .false.

    ! MAIN PROGRAM

    ! Initialize the MPI and PnetCDF session
    call initializeParallelEnvironment

    ! Load the job options file. This first parses the command line arguments.
    ! Then all model switches are read in from a namelist file. This sets up the
    ! run options and points to input files as needed.
    call read_from_job_options()

    ! Load the run setup information based on the metadata in the
    ! initialization netcdf file. The bounds given as an argument to
    ! CLASSIC are used to find the start points (srtx and srty)
    ! in the netcdf file, placing the gridcell on the domain of the
    ! input/output netcdfs. In read_modelsetup we use the netcdf to set
    ! the nmos, ignd,and ilg constants. It also opens the initial conditions
    ! file that is used below in read_initialstate as well as the restart file
    ! that is written to later.
    call read_modelsetup()

    ! Based on the (FLAG add more here)
    call loadoutputDescriptor()


    ! Execute the following only on the main thread (see supportFunctions.f90)
    if (isMainProcess(rank)) then

        ! Generate the output files based on options in the joboptions file
        ! and the parameters of the initilization netcdf file.
        !call create_out_netcdf
        !call generateOutputVariables

        ready_to_go = .true.

    endif

    ! The other threads will wait until the output files are all finished being
    ! created above.
    call MPI_BCAST(ready_to_go, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

    ! Run model over the land grid cells, in parallel
    call processLandCells

    ! Shut down the MPI session
    call MPI_FINALIZE(ierr)

    ! END MAIN PROGRAM

    !------------------

    contains

    subroutine processLandCells

        ! PROCESS LAND CELLS
        ! This section runs the model over all of the land cells. There are validCount valid(i.e. land) cells, stored in validLon and validLat

        ! Since we know the nlat, nmos, ignd, and ilg we can allocate the CLASS and
        ! CTEM variable structures. This has to be done before call to main_driver.
        call alloc_class_vars()
        call alloc_ctem_vars()

        blocks = validCount / size + 1          ! The number of processing blocks
        remainder = mod(validCount, size)       ! The number of cells for the last block

        do i = 1, blocks - 1                    ! Go through every block except for the last one
            cell = (i - 1) * size + rank + 1
            print*,validLon(cell),validLat(cell),validLonIndex(cell),validLatIndex(cell)
            call main_driver(validLon(cell),validLat(cell),validLonIndex(cell),validLatIndex(cell))
        enddo

        cell = (blocks - 1) * size + rank + 1   ! In the last block, process only the existing cells (NEEDS BETTER DESCRIPTION)
        if (rank < remainder) call main_driver(validLon(cell),validLat(cell),&
                                               validLonIndex(cell),validLatIndex(cell))

    end subroutine processLandCells

    !------------------

    subroutine initializeParallelEnvironment

        ! INITIALIZE MPI
        ! This section initializes the MPI environment

        call MPI_INIT(ierr)
        time = MPI_WTIME()
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
        
    end subroutine initializeParallelEnvironment

    !------------------

    logical function isMainProcess(currentProcess)

        ! IS MAIN PROCESS?
        ! This checks to see if we're on the main process of not

        implicit none

        integer, parameter              :: mainProcess = 0
        integer, intent(in)             :: currentProcess

        if (currentProcess == mainProcess) then
            isMainProcess = .true.
        else
            isMainProcess = .false.
        endif

    end function isMainProcess

end program CLASSIC
