!>\file
!> Module for parsing command line arguments to program and reading joboptions files.
module readjobopts

    implicit none

    public :: read_from_job_options
    public :: parsecoords

contains

    ! ---------------------------------------------------

    subroutine read_from_job_options

        !       History:
        !
        !     31 May 2017    - Convert to module, getting ready for new driver.
        !     J. Melton
        !
        !     10 Jan 2017    - igralb no longer supported so removed
        !     J. Melton
        !
        !     9 Nov 2016     - Add the "leap" switch for leap years (.TRUE. if leap years
        !     J.-S. Landry     in the .MET file have data for 366 days, .FALSE. if not)
        !
        !     28  Jul  2016  - Add ability to have changing CO2 but cycling climate
        !     J. Melton        this was for the TRENDY project but generally useful so
        !                      keep in.
        !     3   Feb  2016 - Remove mosaic flag. It is no longer required.
        !     J. Melton

        !     20  Mar. 2015 - Add in new CLASS flags for snow albedos -igralb & isnoalb
        !     J. Melton

        !     4   Sep. 2014 - Add in the transient_run flag.
        !     J. Melton
        !
        !     2   Jul. 2013 - Removed ctem1 and ctem2, replaced with ctem_on
        !     J. Melton
        !
        !     25  Jun. 2013 - Added inibioclim switch for PFTCompetition runs
        !     J. Melton
        !
        !     17  Oct. 2012 - Added the start_bare switch for PFTCompetition runs
        !     J. Melton
        !
        !     25  Apr. 2012 - This subroutine takes in model switches from
        !     J. Melton       a job file and pushes them to RUNCLASSCTEM

        use outputManager, only : myDomain
        use ctem_statevars,     only : c_switch
        use ctem_params, only : icc, ican, l2max, runParamsFile,PFTCompetitionSwitch

        implicit none

        ! -------------

        ! All switches are described in the ConfigurationFiles/template_job_options_file.txt file
        ! and also the user manual.

        ! ctem model switches

        logical, pointer :: transient_run
        integer, pointer :: trans_startyr
        integer, pointer :: metLoop
        logical, pointer :: ctem_on
        integer, pointer :: ncyear
        logical, pointer :: lnduseon
        integer, pointer :: spinfast
        logical, pointer :: cyclemet
        integer, pointer :: nummetcylyrs
        integer, pointer :: metcylyrst
        logical, pointer :: transientCO2
        character(:), pointer :: CO2File
        integer, pointer :: fixedYearCO2
        logical, pointer :: transientCH4
        character(:), pointer :: CH4File
        integer, pointer :: fixedYearCH4
        logical, pointer :: transientPOPD
        integer, pointer :: popcycleyr
        logical, pointer :: dofire
        logical, pointer :: dowetlands
        logical, pointer :: obswetf
        logical, pointer :: PFTCompetition
        logical, pointer :: inibioclim
        logical, pointer :: start_bare
        logical, pointer :: use_netcdf
        character(:), pointer :: met_file
        character(:), pointer :: init_file
        character(:), pointer :: rs_file_to_overwrite
        character(:), pointer :: runparams_file
        character(:), pointer :: output_directory
        character(:), pointer :: xmlFile
        logical, pointer :: leap
        ! -------------
        ! class model switches

        integer, pointer :: idisp
        integer, pointer :: izref
        integer, pointer :: islfd
        integer, pointer :: ipcp
        integer, pointer :: iwf
        INTEGER, pointer :: ITC
        INTEGER, pointer :: ITCG
        INTEGER, pointer :: ITG
        INTEGER, pointer :: IPAI
        INTEGER, pointer :: IHGT
        INTEGER, pointer :: IALC
        INTEGER, pointer :: IALS
        INTEGER, pointer :: IALG
        integer, pointer :: isnoalb

        ! -------------
        ! Output switches

        integer, pointer :: jhhstd  !< day of the year to start writing the half-hourly output
        integer, pointer :: jhhendd !< day of the year to stop writing the half-hourly output
        integer, pointer :: jdstd   !< day of the year to start writing the daily output
        integer, pointer :: jdendd  !< day of the year to stop writing the daily output
        integer, pointer :: jhhsty  !< simulation year (iyear) to start writing the half-hourly output
        integer, pointer :: jhhendy !< simulation year (iyear) to stop writing the half-hourly output
        integer, pointer :: jdsty   !< simulation year (iyear) to start writing the daily output
        integer, pointer :: jdendy  !< simulation year (iyear) to stop writing the daily output
        integer, pointer :: jmosty    !< Year to start writing out the monthly output files. If you want to write monthly outputs right
        logical, pointer :: doperpftoutput    !< Switch for making extra output files that are at the per PFT level
        logical, pointer :: dopertileoutput    !< Switch for making extra output files that are at the per tile level
        logical, pointer :: domonthoutput    !< Switch for making monthly output files (annual are always outputted)
        logical, pointer :: dodayoutput    !< Switch for making daily output files (annual are always outputted)
        logical, pointer :: dohhoutput    !< Switch for making half hourly output files (annual are always outputted)
        character(:), pointer :: Comment   !< Comment about the run that will be written to the output netcdfs

        character(140) :: jobfile
        character(80) :: argbuff
        integer :: argcount, iargc

        ! Order of the namelist and order in the file don't have to match.

        namelist /joboptions/ &
        transient_run, &
        trans_startyr, &
        metLoop, &
        ncyear, &
        cyclemet, &
        nummetcylyrs, &
        metcylyrst, &
        leap, &
        ctem_on, &
        icc, &
        spinfast, &
        lnduseon, &
        transientCO2, &
        CO2File, &
        fixedYearCO2, &
        transientCH4, &
        CH4File, &
        fixedYearCH4, &
        transientPOPD, &
        popcycleyr, &
        PFTCompetition, &
        inibioclim, &
        start_bare, &
        dofire, &
        dowetlands, &
        obswetf, &
        use_netcdf, &
        met_file, &
        init_file, &
        rs_file_to_overwrite, &
        runparams_file, &
        ican, &
        l2max, &
        IDISP, &
        IZREF, &
        ISLFD, &
        IPCP, &
        ITC, &
        ITCG, &
        ITG, &
        IWF, &
        IPAI, &
        IHGT, &
        IALC, &
        IALS, &
        IALG, &
        isnoalb, &
        output_directory, &
        xmlFile, &
        doperpftoutput, &
        dopertileoutput, &
        dohhoutput, &
        JHHSTD, &
        JHHENDD, &
        JHHSTY, &
        JHHENDY, &
        dodayoutput, &
        JDSTD, &
        JDENDD, &
        JDSTY, &
        JDENDY, &
        domonthoutput, &
        JMOSTY, &
        Comment

        ! Point pointers:
        transient_run   => c_switch%transient_run
        trans_startyr   => c_switch%trans_startyr
        metLoop         => c_switch%metLoop
        ctem_on         => c_switch%ctem_on
        ncyear          => c_switch%ncyear
        lnduseon        => c_switch%lnduseon
        spinfast        => c_switch%spinfast
        cyclemet        => c_switch%cyclemet
        nummetcylyrs    => c_switch%nummetcylyrs
        metcylyrst      => c_switch%metcylyrst
        transientCO2    => c_switch%transientCO2
        CO2File         => c_switch%CO2File
        fixedYearCO2    => c_switch%fixedYearCO2
        transientCH4    => c_switch%transientCH4
        CH4File         => c_switch%CH4File
        fixedYearCH4    => c_switch%fixedYearCH4
        transientPOPD   => c_switch%transientPOPD
        popcycleyr      => c_switch%popcycleyr
        dofire          => c_switch%dofire
        dowetlands      => c_switch%dowetlands
        obswetf         => c_switch%obswetf
        PFTCompetition  => c_switch%PFTCompetition
        inibioclim      => c_switch%inibioclim
        start_bare      => c_switch%start_bare
        rs_file_to_overwrite => c_switch%rs_file_to_overwrite
        output_directory => c_switch%output_directory
        xmlFile         => c_switch%xmlFile
        use_netcdf      => c_switch%use_netcdf
        met_file        => c_switch%met_file
        runparams_file  => c_switch%runparams_file
        init_file       => c_switch%init_file
        IDISP           => c_switch%IDISP
        IZREF           => c_switch%IZREF
        ISLFD           => c_switch%ISLFD
        IPCP            => c_switch%IPCP
        ITC             => c_switch%ITC
        ITCG            => c_switch%ITCG
        ITG             => c_switch%ITG
        IWF             => c_switch%IWF
        IPAI            => c_switch%IPAI
        IHGT            => c_switch%IHGT
        IALC            => c_switch%IALC
        IALS            => c_switch%IALS
        IALG            => c_switch%IALG
        isnoalb         => c_switch%isnoalb
        leap            => c_switch%leap
        jhhstd          => c_switch%jhhstd
        jhhendd         => c_switch%jhhendd
        jdstd           => c_switch%jdstd
        jdendd          => c_switch%jdendd
        jhhsty          => c_switch%jhhsty
        jhhendy         => c_switch%jhhendy
        jdsty           => c_switch%jdsty
        jdendy          => c_switch%jdendy
        jmosty          => c_switch%jmosty
        doperpftoutput  => c_switch%doperpftoutput
        dopertileoutput => c_switch%dopertileoutput
        domonthoutput   => c_switch%domonthoutput
        dodayoutput     => c_switch%dodayoutput
        dohhoutput      => c_switch%dohhoutput
        Comment         => c_switch%Comment

        !-------------------------
        !read the joboptions

        argcount = iargc()

        if(argcount .ne. 2)then
            write(*,*)'Usage is as follows'
            write(*,*)' '
            write(*,*)'CLASSIC joboptions_file longitude/latitude'
            write(*,*)' '
            write(*,*)'- joboptions_file - an example can be found '
            write(*,*)'  in the src folder - template_job_options_file.txt.'
            write(*,*)'  descriptions of the various variables '
            write(*,*)'  can be found in read_from_job_options.f90 '
            write(*,*)' '
            write(*,*)'- longitude/latitude '
            write(*,*)'  e.g. 105.23/40.91 '
            write(*,*)' '
            stop
        end if

        call getarg(1,jobfile)

        open(10,file=jobfile,action='read',status='old')

        read(10,nml = joboptions)

        close(10)

        call getarg(2,argbuff)

        if (use_netcdf) then
            call parsecoords(argbuff,myDomain%domainBounds)
        end if

        ! Assign some vars that are passed out
        runParamsFile = runparams_file
        PFTCompetitionSwitch = PFTCompetition

        end subroutine read_from_job_options

! ----------------------------------------------------------------------------------

subroutine parsecoords(coordstring,val)

!subroutine to parse a coordinate string

implicit none

character(45),      intent(in)  :: coordstring
real, dimension(4), intent(out) :: val

character(10), dimension(4) :: cval = '0'

integer :: i
integer :: lasti = 1
integer :: part  = 1

do i=1,len_trim(coordstring)
  if (coordstring(i:i) == '/') then
    cval(part) = coordstring(lasti:i-1)
    lasti=i+1
    part = part + 1
  end if
end do

cval(part) = coordstring(lasti:i-1)

read(cval,*)val

if (part < 4) then
  val(3)=val(2)
  val(4)=val(3)
  val(2)=val(1)
end if

end subroutine parsecoords

!>\defgroup read_from_job_options
!!
!!Joboptions Read-In Subroutine
!!

!>\defgroup parsecoords
!! Parses the bounds info given during the call for the model from the command line.

end module readjobopts
