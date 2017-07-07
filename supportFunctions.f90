module supportFunctions

    implicit none

    integer, parameter              :: mainProcess = 0

contains

    logical function isMainProcess(currentProcess)

        ! IS MAIN PROCESS?
        ! This section checks to see if we're on the main process of not

        implicit none

        integer, intent(in)             :: currentProcess
        if (currentProcess == mainProcess) then
            isMainProcess = .true.
        else
            isMainProcess = .false.
        endif

    end function isMainProcess

end module supportFunctions
