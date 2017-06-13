! Authors: Joe Melton and Ed Wisernig

program CLASSIC
    use mpi
    use io_driver,              only : bounds,lonvect,latvect
    use io_driver,              only : validCount,validLon,validLat
    use model_state_drivers,    only : readModelSetup
    use netcdf_drivers,         only : createOutputNetcdf
    use readjobopts,            only : readFromJobOptions
    use main_driver,            only : CLASSIC_driver
    use ctem_statevars,         only : alloc_ctem_vars
    use class_statevars,        only : alloc_class_vars
    implicit none
    integer, parameter              :: mainProcess = 0
    double precision                :: time
    integer                         :: ierr, rank, size, i, cell, blocks, remainder

                                                ! MAIN PROGRAM
    call initializeParallelEnvironment          ! Initialize the MPI and PnetCDF session
    call readFromJobOptions                     ! Load the project config file
    call readModelSetup                         ! Load the model setup information
    if (isMainProcess()) then                   ! Execute only on the main thread
        call createOutputNetcdf                 ! Generate the output files
    endif
    call processLandCells                       ! Process the land grid cells, in parallel
    call finalizeParallelEnvironment            ! Close PnetCDF and shut down the MPI session
                                                ! END MAIN PROGRAM

contains

    ! PROCESS LAND CELLS
    ! This section processes all of the land cells. There are validCount valid(i.e. land) cells, stored in validLon and validLat
    subroutine processLandCells
        call alloc_class_vars()
        call alloc_ctem_vars()

        blocks = validCount / size + 1          ! The number of processing blocks
        remainder = mod(validCount, size)       ! The number of cells for the last block
        do i = 1, blocks - 1                    ! Go through every block except for the last one
            cell = (i - 1) * size + rank + 1
            call CLASSIC_driver(cell)
        enddo
        cell = (blocks - 1) * size + rank + 1   ! In the last block, process only the existing cells (NEEDS BETTER DESCRIPTION)
        if (rank < remainder) call CLASSIC_driver(cell)
    end subroutine processLandCells

    ! INITIALIZE MPI and PnetCDf
    ! This section initializes the MPI environment and PnetCDF files
    subroutine initializeParallelEnvironment
        call MPI_INIT(ierr)
        time = MPI_WTIME()
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
        ! HERE - PnetCDF init code comes here
    end subroutine initializeParallelEnvironment

    ! FINALIZE MPI and PnetCDF
    ! This section wraps up the whole MPI and PnetCDF megillah
    subroutine finalizeParallelEnvironment
        ! HERE - PnetCDF close session code comes here
        call MPI_FINALIZE(ierr)
    end subroutine finalizeParallelEnvironment

    ! IS MAIN PROCESS
    ! This section checks to see if we're on the main process of not
    logical function isMainProcess()
        implicit none
        if (rank == mainProcess) then
            isMainProcess = .true.
        else
            isMainProcess = .false.
        endif
    end function isMainProcess
end program CLASSIC
