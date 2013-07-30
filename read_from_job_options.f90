subroutine read_from_job_options(argbuff,mosaic,ctemloop,ctem_on,ncyear,lnduseon,spinfast,cyclemet, &
                  nummetcylyrs,metcylyrst,co2on,setco2conc,popdon,popcycleyr, &
                  parallelrun,dofire,compete,inibioclim,start_bare,rsfile,start_from_rs,idisp,izref,islfd,ipcp,itc,itcg, &
                  itg,iwf,ipai,ihgt,ialc,ials,ialg,jhhstd,jhhendd,jdstd, & 
                  jdendd,jhhsty,jhhendy,jdsty,jdendy)

!#ifdef nagf95
!use f90_unix
!#endif

!           Canadian Terrestrial Ecosystem Model (CTEM) 
!                    Joboptions Read-In Subroutine 
!
!     2   Jul. 2013 - Removed ctem1 and ctem2, replaced with ctem_on
!     J. Melton
!
!     25  Jun. 2013 - Added inibioclim switch for compete runs
!     J. Melton
!
!     17  Oct. 2012 - Added the start_bare switch for compete runs
!     J. Melton
!
!     25  Apr. 2012 - This subroutine takes in model switches from
!     J. Melton       a job file and pushes them to RUNCLASSCTEM
!		      

implicit none

! -------------
! ctem model switches

character(80), intent(out) :: argbuff !prefix of file names

logical, intent(out) :: mosaic   ! true if the run is in mosaic mode, otherwise it
                                 ! is a composite run

integer, intent(out) :: ctemloop ! no. of times the .met file is to be read. this
                    	         ! option is useful to see how ctem's c pools
                    	         ! equilibrate when driven with same climate data
                    	         ! over and over again.

logical, intent(out) :: ctem_on  ! set this to true for using ctem simulated dynamic
 				 ! lai and canopy mass, else class simulated specified
 				 ! lai and canopy mass are used. with this switch on,
 				 ! all ctem subroutines are run.

integer, intent(out) :: ncyear   ! no. of years in the .met file. 

logical, intent(out) :: lnduseon ! set this to 1 if land use change is to be
 				 ! implimented by reading in the fractions of 9 ctem
 				 ! pfts from a file. keep in mind that once on, luc read-in is
                                 ! also influenced by the cyclemet and popcycleyr
                                 ! switches

integer, intent(out) :: spinfast ! set this to a higher number up to 10 to spin up
 				 ! soil carbon pool faster

logical, intent(out) :: cyclemet ! to cycle over only a fixed number of years 
 				 ! (nummetcylyrs) starting at a certain year (metcylyrst)
 				 ! if cyclemet, then put co2on = false and set an appopriate setco2conc, also
 				 ! if popdon is true, it will choose the popn and luc data for year
 				 ! metcylyrst and cycle on that.

integer, intent(out) :: nummetcylyrs ! years of the climate file to spin up on repeatedly
 				 ! ignored if cyclemet is false

integer, intent(out) :: metcylyrst   ! climate year to start the spin up on
 				 ! ignored if cyclemet is false

logical, intent(out) :: co2on    ! use co2 time series, set to false if cyclemet is true

real, intent(out) :: setco2conc  ! set the value of atmospheric co2 if co2on is false.

logical, intent(out) :: popdon   ! if set true use population density data to calculate fire extinguishing 
 				 ! probability and probability of fire due to human causes, 
 				 ! or if false, read directly from .ctm file

integer, intent(out) :: popcycleyr ! popd and luc year to cycle on when cyclemet is true, set to -9999
				 ! to cycle on metcylyrst for both popd and luc. if cyclemet is false
                                 ! this defaults to -9999, which will then cause the model to cycle on
                                 ! whatever is the first year in the popd and luc datasets

logical, intent(out) :: parallelrun ! set this to be true if model is run in parallel mode for 
 				 ! multiple grid cells, output is limited to monthly & yearly 
 				 ! grid-mean only. else the run is in stand alone mode, in which 
 				 ! output includes half-hourly and daily and mosaic-mean as well.

logical, intent(out) :: dofire   ! if true the fire/disturbance subroutine will be used.

logical, intent(out) :: compete  ! set this to true if competition between pfts is
 				 ! to be implimented

logical, intent(out) :: inibioclim  ! set this to true if competition between pfts is
 				    ! to be implimented and you have the mean climate values
                                    ! in the ctm files.

logical, intent(out) :: start_bare !set this to true if competition is true, and if you wish
                                 ! to start from bare ground. if this is set to false, the 
                                 ! ini and ctm file info will be used to set up the run.

logical, intent(out) :: rsfile   ! set this to true if restart files (.ini_rs and .ctm_rs)   
 				 ! are written at the end of each year. these files are  
 				 ! necessary for checking whether the model reaches 
 				 ! equilibrium after running for a certain years. 
 				 ! set this to false if restart files are not needed 
 				 ! (known how many years the model will run)
logical, intent(out) :: start_from_rs ! if true, this option copies the _RS INI and CTM files
                                 ! to be the .INI and .CTM files and then starts the run as per normal.
                                 ! it is handy when spinning up so you don't have to do a complicated copying of the
                                 ! RS files to restart from them.

! -------------
! class model switches

integer, intent(out) :: idisp    ! if idisp=0, vegetation displacement heights are ignored,
				 ! because the atmospheric model considers these to be part
				 ! of the "terrain".
				 ! if idisp=1, vegetation displacement heights are calculated.

integer, intent(out) :: izref    ! if izref=1, the bottom of the atmospheric model is taken
				 ! to lie at the ground surface.
				 ! if izref=2, the bottom of the atmospheric model is taken
				 ! to lie at the local roughness height.

integer, intent(out) :: islfd    ! if islfd=0, drcoef is called for surface stability corrections
				 ! and the original gcm set of screen-level diagnostic calculations 
				 ! is done.
				 ! if islfd=1, drcoef is called for surface stability corrections
				 ! and sldiag is called for screen-level diagnostic calculations. 
				 ! if islfd=2, flxsurfz is called for surface stability corrections
				 ! and diasurf is called for screen-level diagnostic calculations. 

integer, intent(out) :: ipcp     ! if ipcp=1, the rainfall-snowfall cutoff is taken to lie at 0 c.
				 ! if ipcp=2, a linear partitioning of precipitation betweeen 
				 ! rainfall and snowfall is done between 0 c and 2 c.
				 ! if ipcp=3, rainfall and snowfall are partitioned according to
				 ! a polynomial curve between 0 c and 6 c.

integer, intent(out) :: iwf     ! if iwf=0, only overland flow and baseflow are modelled, and
				! the ground surface slope is not modelled.
				! if iwf=n (0<n<4), the watflood calculations of overland flow 
				! and interflow are performed; interflow is drawn from the top 
				! n soil layers.

! itc, itcg and itg are switches to choose the iteration scheme to
! be used in calculating the canopy or ground surface temperature
! respectively.  if the switch is set to 1, a bisection method is
! used; if to 2, the newton-raphson method is used.
INTEGER, INTENT(OUT) :: ITC
INTEGER, INTENT(OUT) :: ITCG
INTEGER, INTENT(OUT) :: ITG

! if ipai, ihgt, ialc, ials and ialg are zero, the values of 
! plant area index, vegetation height, canopy albedo, snow albedo
! and soil albedo respectively calculated by class are used.
! if any of these switches is set to 1, the value of the
! corresponding parameter calculated by class is overridden by
! a user-supplied input value.
!      
INTEGER, INTENT(OUT) :: IPAI
INTEGER, INTENT(OUT) :: IHGT
INTEGER, INTENT(OUT) :: IALC
INTEGER, INTENT(OUT) :: IALS
INTEGER, INTENT(OUT) :: IALG

! -------------
! classctem output switches

! >>>> note: if you wish to use the values in the .ini file, set all to -9999 in the job options file
!            and the .ini file will be used.

integer, intent(out) :: jhhstd    ! day of the year to start writing the half-hourly output
integer, intent(out) :: jhhendd   ! day of the year to stop writing the half-hourly output
integer, intent(out) :: jdstd     ! day of the year to start writing the daily output
integer, intent(out) :: jdendd    ! day of the year to stop writing the daily output
integer, intent(out) :: jhhsty    ! simulation year (iyear) to start writing the half-hourly output
integer, intent(out) :: jhhendy   ! simulation year (iyear) to stop writing the half-hourly output
integer, intent(out) :: jdsty     ! simulation year (iyear) to start writing the daily output
integer, intent(out) :: jdendy    ! simulation year (iyear) to stop writing the daily output

! -------------

namelist /joboptions/ &
  mosaic,             &
  ctemloop,           &
  ctem_on,            &
  ncyear,             &
  lnduseon,           &
  spinfast,           &
  cyclemet,           &
  nummetcylyrs,       &
  metcylyrst,         &
  co2on,              &
  setco2conc,         &
  popdon,             &
  popcycleyr,         &
  parallelrun,        &
  dofire,             &
  compete,            &
  inibioclim,         &
  start_bare,         &
  rsfile,             &
  start_from_rs,      &
  IDISP,              &
  IZREF,              &
  ISLFD,              &
  IPCP,               &
  ITC,                &
  ITCG,               &
  ITG,                &
  IWF,                &
  IPAI,               &
  IHGT,               &
  IALC,               &
  IALS,               &
  IALG,               &
  jhhstd,             &
  jhhendd,            &
  jdstd,              &
  jdendd,             &
  jhhsty,             &
  jhhendy,            &
  jdsty,              &
  jdendy

character(140) :: jobfile
integer :: argcount, iargc 

!-------------------------
!read the joboptions

argcount = iargc()

       if(argcount .ne. 2)then
         write(*,*)'usage is as follows'
         write(*,*)' '
         write(*,*)'CLASS36CTEM joboptions_file site_name'
         write(*,*)' '
         write(*,*)'- joboptions_file - an example can be found '
         write(*,*)'  in the src folder - template_job_options_file.txt.'
         write(*,*)'  descriptions of the various variables '
         write(*,*)'  can be found in read_from_job_options.f90 '
         write(*,*)' '
         write(*,*)'- site_name is the prefix of your input files '
         write(*,*)' '
         stop
      end if


call getarg(1,jobfile)

open(10,file=jobfile,status='old')

read(10,nml = joboptions)

close(10)

call getarg(2,argbuff)


end subroutine read_from_job_options

