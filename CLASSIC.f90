!>\file
!! Principle driver program to run CLASSIC in stand-alone mode using specified boundary
!! conditions and atmospheric forcing.
!!
program CLASSIC

use io_driver,          only : bounds
use input_dataset_drivers, only : read_modelsetup, read_initialstate
!use main_driver, only : CLASSIC_driver

implicit none

! Local variables
real :: longitude, latitude

! ------------

! All model switches are read in from a namelist file
call read_from_job_options()

!> First we set up the run boundaries based on the metadata in the initialization netcdf file.
!! The bounds are used to find the srtx and srty in the netcdf file, placing the gridcell on the
!! domain of the input/output netcdfs. In read_modelsetup we use the netcdf to set the nmos, ignd,
!!and ilg constants. It also opens the initial conditions file that is used in read_initialstate.
call read_modelsetup()
!
!> Set up the longitude and latitude of this gridcell based on the bounds
longitude = bounds(1)
latitude = bounds(3)

!> Next we create all the output files for the model run based on options in the joboptions file
!! and the parameters of the initilization netcdf file.

! Then in parallel we call the main model:
!call CLASSIC_driver()


end program