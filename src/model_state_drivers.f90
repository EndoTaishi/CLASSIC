!> Central driver to read in, and write out all model state variables (replacing INI and CTM files)
!! as well as the model inputs such as MET, population density, land use change, CO2 etc.

module model_state_drivers

    ! J. Melton
    ! Nov 2016

    use fileIOModule

    implicit none

    public  :: read_modelsetup
    public  :: read_initialstate
    public  :: write_restart
    public  :: getInput
    public  :: updateInput
    public  :: getMet
    public  :: updateMet
    public  :: deallocInput
    private :: closestCell

    integer, dimension(:), allocatable :: CO2Time           !< The time (years) from the CO2File
    real, dimension(:), allocatable :: CO2FromFile          !< The array of CO2 values (ppm) from the CO2File
    integer, dimension(:), allocatable :: CH4Time           !< The time (years) from the CH4File
    real, dimension(:), allocatable :: CH4FromFile          !< The array of CH4 values (ppm) from the CH4File
    integer, dimension(:), allocatable :: POPDTime          !< The time (years) from the population density file
    real, dimension(:), allocatable :: POPDFromFile         !< The array of CH4 values (ppm) from the POPDFile
    real, dimension(:), allocatable :: LGHTTime             !< The time from the lightning density file (usually months)
    real, dimension(:), allocatable :: LGHTFromFile         !< The array of lightning density from the LGHTFile
    integer, dimension(:), allocatable :: LUCTime           !< The time from the LUC file
    real, dimension(:,:), allocatable :: LUCFromFile        !< The array of LUC from the LUCFile

    real, dimension(:), allocatable :: metTime              !< The time from the Met file
    real, dimension(:), allocatable :: metFss               !< Incoming shortwave radiation from metFile \f$[W m^{-2} ]\f$
    real, dimension(:), allocatable :: metFdl               !< Incoming longwave radiation from metFile \f$[W m^{-2} ]\f$
    real, dimension(:), allocatable :: metPre               !< Precipitation from metFile \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(:), allocatable :: metTa                !< Air temperature from metFile (Celsius)
    real, dimension(:), allocatable :: metQa                !< Specific humidity from metFile
    real, dimension(:), allocatable :: metUv                !< Wind speed from metFile
    real, dimension(:), allocatable :: metPres              !< Atmospheric pressure from metFile

    integer :: metFssId                             !> netcdf file id for the incoming shortwave radiation meteorology file
    integer :: metFdlId                             !> netcdf file id for the incoming longwave radiation meteorology file
    integer :: metPreId                             !> netcdf file id for the precipitation meteorology file
    integer :: metTaId                              !> netcdf file id for the air temperature meteorology file
    integer :: metQaId                              !> netcdf file id for the specific humidity meteorology file
    integer :: metUvId                              !> netcdf file id for the wind speed meteorology file
    integer :: metPresId                            !> netcdf file id for the atmospheric pressure meteorology file
    integer :: initid                               !> netcdf file id for the model initialization file
    integer :: rsid                                 !> netcdf file id for the model restart file
    integer :: co2id                                !> netcdf file id for the CO2 input file
    integer :: ch4id                                !> netcdf file id for the CH4 input file
    integer :: popid                                !> netcdf file id for the population density input file
    integer :: lghtid                               !> netcdf file id for the lightning density input file
    integer :: lucid                                !> netcdf file id for the land use change input file

    real :: metInputTimeStep                        !> The timestep of the read in meteorology (hours)

contains

    !---

    !>\ingroup model_state_drivers_read_modelsetup
    !!@{
    !> Reads in the model setup from the netcdf initialization file.
    !> The number of latitudes is always 1 offline while the maximum number of
    !> mosaics (nmos), the number of soil layers (ignd), are read from the netcdf.
    !> ilg is then calculated from nlat and nmos.

    subroutine read_modelsetup

        ! J. Melton
        ! Feb 2017

        use ctem_statevars,     only : c_switch
        use ctem_params, only : nmos,nlat,ignd,ilg  ! These are set in this subroutine!
        use outputManager, only : myDomain

        implicit none

        ! pointers:
        character(350), pointer          :: init_file
        character(350), pointer          :: rs_file_to_overwrite
        character(350), pointer          :: CO2File
        character(350), pointer          :: CH4File
        character(350), pointer          :: POPDFile
        character(350), pointer          :: LGHTFile
        character(350), pointer          :: LUCFile
        character(350), pointer          :: metFileFss
        character(350), pointer          :: metFileFdl
        character(350), pointer          :: metFilePre
        character(350), pointer          :: metFileTa
        character(350), pointer          :: metFileQa
        character(350), pointer          :: metFileUv
        character(350), pointer          :: metFilePres
        logical, pointer                 :: ctem_on
        logical, pointer                 :: dofire
        logical, pointer                 :: lnduseon
        integer, pointer                 :: fixedYearLUC

        ! Local vars
        integer, allocatable, dimension(:,:) :: mask
        integer :: i, j
        integer :: totlon,totlat,totsize
        integer, dimension(1) :: pos
        integer, dimension(2) :: xpos,ypos
        integer, dimension(:,:), allocatable :: nmarray
        integer :: lonloc,latloc

        ! point pointers:
        init_file               => c_switch%init_file
        rs_file_to_overwrite    => c_switch%rs_file_to_overwrite
        CO2File                 => c_switch%CO2File
        CH4File                 => c_switch%CH4File
        POPDFile                => c_switch%POPDFile
        LGHTFile                => c_switch%LGHTFile
        LUCFile                 => c_switch%LUCFile
        ctem_on                 => c_switch%ctem_on
        dofire                  => c_switch%dofire
        lnduseon                => c_switch%lnduseon
        fixedYearLUC            => c_switch%fixedYearLUC
        metFileFss              => c_switch%metFileFss
        metFileFdl              => c_switch%metFileFdl
        metFilePre              => c_switch%metFilePre
        metFileTa               => c_switch%metFileTa
        metFileQa               => c_switch%metFileQa
        metFileUv               => c_switch%metFileUv
        metFilePres             => c_switch%metFilePres

        ! ------------

        !> First, open initial conditions file.

        initid = ncOpen(init_file, NF90_NOWRITE)

        !> Next, retrieve dimensions. We assume the file has 'lon' and 'lat' for
        !! names of longitude and latitude.

        totlon = ncGetDimLen(initid,'lon')
        totlat = ncGetDimLen(initid,'lat')

        !>calculate the number and indices of the pixels to be calculated
        allocate(myDomain%allLonValues(totlon), myDomain%allLatValues(totlat))

        myDomain%allLonValues = ncGetDimValues(initid, 'lon', count = (/totlon/))
        myDomain%allLatValues = ncGetDimValues(initid, 'lat', count = (/totlat/))

        !> Try and catch if the user has put in lon values from -180 to 180 or 0 to 360
        !! when the input file expects the opposite.
        if (myDomain%domainBounds(1) < 0. .and. myDomain%allLonValues(1) >= 0.) then
            myDomain%domainBounds(1) = 360. + myDomain%domainBounds(1)
            print *,'Based on init_file, adjusted your domain to',myDomain%domainBounds(1)
        end if
        if (myDomain%domainBounds(2) < 0. .and. myDomain%allLonValues(1) >= 0.) then
            myDomain%domainBounds(2) = 360. + myDomain%domainBounds(2)
            print *,'Based on init_file, adjusted your domain to',myDomain%domainBounds(2)
        end if
        if (myDomain%domainBounds(1) > 180. .and. myDomain%allLonValues(1) < 0.) then
            myDomain%domainBounds(1) = myDomain%domainBounds(1) - 360.
            print *,'Based on init_file, adjusted your domain to',myDomain%domainBounds(1)
        end if
        if (myDomain%domainBounds(2) > 180. .and. myDomain%allLonValues(1) < 0.) then
            myDomain%domainBounds(2) = myDomain%domainBounds(2) - 360.
            print *,'Based on init_file, adjusted your domain to',myDomain%domainBounds(2)
        end if

        !> Check that our domain is within the longitude and latitude limits of
        !! the input files. Otherwise print a warning. Primarily we are trying to
        !! catch instances where the input file runs from 0 to 360 longitude while
        !! the user expects -180 to 180.
        if (myDomain%domainBounds(1) < myDomain%allLonValues(1)) then !W most lon
            print*,'=>Your domain bound ', myDomain%domainBounds(1),' is outside of',&
                ' the limits of the init_file ',myDomain%allLonValues(1)
        else if (myDomain%domainBounds(2) > myDomain%allLonValues(ubound(myDomain%allLonValues,1))) then ! E most lon
            print*,'=>Your domain bound ', myDomain%domainBounds(2),' is outside of',&
                ' the limits of the init_file ',myDomain%allLonValues(ubound(myDomain%allLonValues,1))
        else if (myDomain%domainBounds(3) < myDomain%allLatValues(1)) then !S most lat
            print*,'=>Your domain bound ', myDomain%domainBounds(3),' is outside of',&
                ' the limits of the init_file ',myDomain%allLatValues(1)
        else if (myDomain%domainBounds(4) > myDomain%allLatValues(ubound(myDomain%allLatValues,1))) then !N most lat
            print*,'=>Your domain bound ', myDomain%domainBounds(4),' is outside of',&
                ' the limits of the init_file ',myDomain%allLatValues(ubound(myDomain%allLatValues,1))
        end if

        !> Based on the domainBounds, we make vectors of the cells to be run.
        pos = minloc(abs(myDomain%allLonValues - myDomain%domainBounds(1)))
        xpos(1) = pos(1)

        pos = minloc(abs(myDomain%allLonValues - myDomain%domainBounds(2)))
        xpos(2) = pos(1)

        pos = minloc(abs(myDomain%allLatValues - myDomain%domainBounds(3)))
        ypos(1) = pos(1)

        pos = minloc(abs(myDomain%allLatValues - myDomain%domainBounds(4)))
        ypos(2) = pos(1)

        myDomain%srtx = minval(xpos)
        myDomain%srty = minval(ypos)

        if (myDomain%allLonValues(myDomain%srtx) < myDomain%domainBounds(1) .and.&
            myDomain%domainBounds(2) /= myDomain%domainBounds(1)) myDomain%srtx = myDomain%srtx + 1
        myDomain%cntx = 1 + abs(maxval(xpos) - myDomain%srtx)

        if (myDomain%allLatValues(myDomain%srty) < myDomain%domainBounds(3) .and.&
            myDomain%domainBounds(4) /= myDomain%domainBounds(3)) myDomain%srty = myDomain%srty + 1
        myDomain%cnty = 1 + abs(maxval(ypos) - myDomain%srty)

        !> Save the longitudes and latitudes over the region of interest for making the
        !! output files.
        totsize = myDomain%cntx * myDomain%cnty
        allocate(myDomain%latLandCell(totsize),&
                 myDomain%lonLandCell(totsize),&
                 myDomain%latLandIndex(totsize),&
                 myDomain%lonLandIndex(totsize),&
                 myDomain%latLocalIndex(totsize),&
                 myDomain%lonLocalIndex(totsize),&
                 myDomain%latUnique(myDomain%cnty),&
                 myDomain%lonUnique(myDomain%cntx))

        !> Retrieve the number of soil layers (set ignd!)
        ignd = ncGetDimLen(initid, 'layer')

        !> Grab the model domain. We use GC since it is the land cells we want to run the model over.
        !! the 'Mask' variable is all land (but we don't run over Antarctica).
        allocate(mask(myDomain%cntx, myDomain%cnty))
        mask = ncGet2DVar(initid, 'GC', start = [myDomain%srtx, myDomain%srty],&
                          count = [myDomain%cntx, myDomain%cnty],format = [myDomain%cntx, myDomain%cnty])

        myDomain%LandCellCount = 0
        do i = 1, myDomain%cntx
            do j = 1, myDomain%cnty
                if (mask(i,j) .eq. -1) then
                    !print*, "(", i, ",", j, ") or (", myDomain%allLonValues(i + myDomain%srtx - 1)&
                    !, ",", myDomain%allLatValues(j + myDomain%srty - 1), ") is land"
                    myDomain%LandCellCount = myDomain%LandCellCount + 1
                    myDomain%lonLandCell(myDomain%LandCellCount) = myDomain%allLonValues(i + myDomain%srtx - 1)
                    myDomain%lonLandIndex(myDomain%LandCellCount) = i + myDomain%srtx - 1
                    myDomain%lonLocalIndex(myDomain%LandCellCount) = i
                    myDomain%lonUnique(i) = myDomain%allLonValues(i + myDomain%srtx - 1)
                    myDomain%latLandCell(myDomain%LandCellCount) = myDomain%allLatValues(j + myDomain%srty - 1)
                    myDomain%latLandIndex(myDomain%LandCellCount) = j + myDomain%srty - 1
                    myDomain%latLocalIndex(myDomain%LandCellCount) = j
                    myDomain%latUnique(j) = myDomain%allLatValues(j + myDomain%srty - 1)
                else !keep track of the non-land too for the making of the output files.
                    myDomain%lonUnique(i) = myDomain%allLonValues(i + myDomain%srtx - 1)
                    myDomain%latUnique(j) = myDomain%allLatValues(j + myDomain%srty - 1)
                endif
            enddo
        enddo

        if (myDomain%LandCellCount == 0) then
            print*,'=>Your domain is not land my friend.'
            if (myDomain%domainBounds(1) == myDomain%domainBounds(2)) then !point run
                lonloc = closestCell(initid,'lon',myDomain%domainBounds(1))
                latloc = closestCell(initid,'lat',myDomain%domainBounds(3))
                print*,'Closest grid cell is ',myDomain%allLonValues(lonloc),'/',myDomain%allLatValues(latloc)
                print*,'but that may not be land. Check your input files to be sure'
            end if
        end if
        
        nlat = 1

        !> To determine nmos, we use the largest number in the input file variable nmtest
        !! for the region we are running.
        allocate(nmarray(myDomain%cntx, myDomain%cnty))
        nmarray = ncGet2DVar(initid, 'nmtest', start = [myDomain%srtx, myDomain%srty],&
                             count = [myDomain%cntx, myDomain%cnty],format = [myDomain%cntx, myDomain%cnty])
        nmos= maxval(nmarray)

        !> Determine the size of ilg which is nlat times nmos

        ilg = nlat * nmos

        !> Lastly, open some files so they are ready
        rsid = ncOpen(rs_file_to_overwrite, nf90_write)

        if (ctem_on) then
            co2id = ncOpen(CO2File, nf90_nowrite)
            ch4id = ncOpen(CH4File, nf90_nowrite)
            if (dofire) then
                popid = ncOpen(POPDFile, nf90_nowrite)
                lghtid = ncOpen(LGHTFile, nf90_nowrite)
            end if
            if (lnduseon .or. (fixedYearLUC .ne. -9999)) then
                lucid = ncOpen(LUCFile, nf90_nowrite)
            end if
        end if

        !> Open the meteorological forcing files
        metFssId    = ncOpen(metFileFss, nf90_nowrite)
        metFdlId    = ncOpen(metFileFdl, nf90_nowrite)
        metPreId    = ncOpen(metFilePre, nf90_nowrite)
        metTaId     = ncOpen(metFileTa, nf90_nowrite)
        metQaId     = ncOpen(metFileQa, nf90_nowrite)
        metUvId     = ncOpen(metFileUv, nf90_nowrite)
        metPresId   = ncOpen(metFilePres, nf90_nowrite)

    end subroutine read_modelsetup

    !>@}
    ! ------------------------------------------------------------------------------------

    !>\ingroup model_state_drivers_read_initialstate
    !!@{
    !> Reads in the model initial conditions for both physics and biogeochemistry (if CTEM on)

    subroutine read_initialstate(lonIndex,latIndex)

        ! J. Melton
        ! Nov 2016

        use ctem_statevars,     only : c_switch,vrot,vgat
        use class_statevars,    only : class_rot,class_gat
        use ctem_params,        only : icc,iccp1,nmos,ignd,ilg,icp1,nlat,ican,abszero,pi,crop

        implicit none

        ! arguments
        integer, intent(in) :: lonIndex,latIndex

        ! pointers:
        real, pointer, dimension(:,:,:) :: FCANROT      !<Maximum fractional coverage of modelled
        real, pointer, dimension(:,:)   :: FAREROT
        real, pointer, dimension(:,:,:) :: RSMNROT
        real, pointer, dimension(:,:,:) :: QA50ROT
        real, pointer, dimension(:,:,:) :: VPDAROT
        real, pointer, dimension(:,:,:) :: VPDBROT
        real, pointer, dimension(:,:,:) :: PSGAROT
        real, pointer, dimension(:,:,:) :: PSGBROT
        real, pointer, dimension(:,:,:) :: ALVCROT
        real, pointer, dimension(:,:,:) :: ALICROT
        real, pointer, dimension(:,:,:) :: PAMNROT
        real, pointer, dimension(:,:,:) :: PAMXROT
        real, pointer, dimension(:,:,:) :: LNZ0ROT
        real, pointer, dimension(:,:,:) :: CMASROT
        real, pointer, dimension(:,:,:) :: ROOTROT
        real, pointer, dimension(:,:)   :: DRNROT
        real, pointer, dimension(:,:)   :: SDEPROT
        real, pointer, dimension(:,:)   :: XSLPROT
        real, pointer, dimension(:,:)   :: GRKFROT
        real, pointer, dimension(:,:)   :: WFSFROT
        real, pointer, dimension(:,:)   :: WFCIROT
        integer, pointer, dimension(:,:)   :: MIDROT
        real, pointer, dimension(:,:,:) :: SANDROT
        real, pointer, dimension(:,:,:) :: CLAYROT
        real, pointer, dimension(:,:,:) :: ORGMROT
        real, pointer, dimension(:,:,:) :: TBARROT
        real, pointer, dimension(:,:,:) :: THLQROT
        real, pointer, dimension(:,:,:) :: THICROT
        real, pointer, dimension(:)   :: DELZ
        real, pointer, dimension(:)   :: ZBOT
        real, pointer, dimension(:,:)   :: TCANROT
        real, pointer, dimension(:,:)   :: TSNOROT
        real, pointer, dimension(:,:)   :: TPNDROT
        real, pointer, dimension(:,:)   :: ZPNDROT
        real, pointer, dimension(:,:)   :: RCANROT
        real, pointer, dimension(:,:)   :: SCANROT
        real, pointer, dimension(:,:)   :: SNOROT
        real, pointer, dimension(:,:)   :: ALBSROT
        real, pointer, dimension(:,:)   :: RHOSROT
        real, pointer, dimension(:,:)   :: GROROT
        real, pointer, dimension(:)     :: ZRFHROW !<
        real, pointer, dimension(:)     :: ZRFMROW !<
        real, pointer, dimension(:)     :: DLATROW !<
        real, pointer, dimension(:)     :: DLONROW !<
        real, pointer, dimension(:)     :: GCROW   !<Type identifier for grid cell (1 = sea ice, 0 = ocean, -1 = land)
        real, pointer, dimension(:)     :: ZBLDROW !<
        real, pointer, dimension(:)     :: RADJROW !<Latitude of grid cell (positive north of equator) [rad]
        real, pointer, dimension(:)     :: Z0ORROW !<
        real, pointer, dimension(:)     :: GGEOROW !<Geothermal heat flux at bottom of soil profile \f$[W m^{-2} ]\f$
        real, pointer, dimension(:,:)   :: SOCIROT

        real, pointer, dimension(:,:)   :: TBASROT !<
        logical, pointer :: ctem_on
        logical, pointer :: dofire
        logical, pointer :: PFTCompetition
        logical, pointer :: inibioclim
        logical, pointer :: dowetlands
        logical, pointer :: start_bare
        logical, pointer :: lnduseon
        logical, pointer :: obswetf
        real, pointer, dimension(:,:,:) :: ailcminrow           !
        real, pointer, dimension(:,:,:) :: ailcmaxrow           !
        real, pointer, dimension(:,:,:) :: fcancmxrow           !
        real, pointer, dimension(:,:,:) :: gleafmasrow          !
        real, pointer, dimension(:,:,:) :: bleafmasrow          !
        real, pointer, dimension(:,:,:) :: stemmassrow          !
        real, pointer, dimension(:,:,:) :: rootmassrow          !
        real, pointer, dimension(:,:,:) :: pstemmassrow         !
        real, pointer, dimension(:,:,:) :: pgleafmassrow        !
        real, pointer, dimension(:,:) :: twarmm            !< temperature of the warmest month (c)
        real, pointer, dimension(:,:) :: tcoldm            !< temperature of the coldest month (c)
        real, pointer, dimension(:,:) :: gdd5              !< growing degree days above 5 c
        real, pointer, dimension(:,:) :: aridity           !< aridity index, ratio of potential evaporation to precipitation
        real, pointer, dimension(:,:) :: srplsmon          !< number of months in a year with surplus water i.e.precipitation more than potential evaporation
        real, pointer, dimension(:,:) :: defctmon          !< number of months in a year with water deficit i.e.precipitation less than potential evaporation
        real, pointer, dimension(:,:) :: anndefct          !< annual water deficit (mm)
        real, pointer, dimension(:,:) :: annsrpls          !< annual water surplus (mm)
        real, pointer, dimension(:,:) :: annpcp            !< annual precipitation (mm)
        real, pointer, dimension(:,:) :: dry_season_length !< length of dry season (months)
        real, pointer, dimension(:,:,:) :: litrmassrow
        real, pointer, dimension(:,:,:) :: soilcmasrow
        real, pointer, dimension(:,:) :: extnprob
        real, pointer, dimension(:,:) :: prbfrhuc
        real, pointer, dimension(:,:,:) :: mlightng
        integer, pointer, dimension(:,:,:) :: lfstatusrow
        integer, pointer, dimension(:,:,:) :: pandaysrow
        integer, pointer, dimension(:,:) :: stdaln
        real, pointer, dimension(:,:,:) :: slopefrac
        integer, pointer, dimension(:,:) :: ipeatlandrow   !<Peatland switch: 0 = not a peatland, 1= bog, 2 = fen
        real, pointer, dimension(:,:) :: Cmossmas          !<C in moss biomass, \f$kg C/m^2\f$
        real, pointer, dimension(:,:) :: litrmsmoss        !<moss litter mass, \f$kg C/m^2\f$
        real, pointer, dimension(:,:) :: dmoss             !<depth of living moss (m)
        real, pointer, dimension(:) :: grclarea            !<area of the grid cell, \f$km^2\f$

        ! local variables

        integer :: i,m,j
        real, dimension(ilg,2) :: crop_temp_frac
        real, parameter :: TFREZ = 273.16

        ! point pointers:
        ctem_on           => c_switch%ctem_on
        dofire            => c_switch%dofire
        PFTCompetition    => c_switch%PFTCompetition
        inibioclim        => c_switch%inibioclim
        dowetlands        => c_switch%dowetlands
        start_bare        => c_switch%start_bare
        lnduseon          => c_switch%lnduseon
        obswetf           => c_switch%obswetf
        ailcminrow        => vrot%ailcmin
        ailcmaxrow        => vrot%ailcmax
        fcancmxrow        => vrot%fcancmx
        gleafmasrow       => vrot%gleafmas
        bleafmasrow       => vrot%bleafmas
        stemmassrow       => vrot%stemmass
        rootmassrow       => vrot%rootmass
        pstemmassrow      => vrot%pstemmass
        pgleafmassrow     => vrot%pgleafmass
        twarmm            => vrot%twarmm
        tcoldm            => vrot%tcoldm
        gdd5              => vrot%gdd5
        aridity           => vrot%aridity
        srplsmon          => vrot%srplsmon
        defctmon          => vrot%defctmon
        anndefct          => vrot%anndefct
        annsrpls          => vrot%annsrpls
        annpcp            => vrot%annpcp
        dry_season_length => vrot%dry_season_length
        litrmassrow       => vrot%litrmass
        soilcmasrow       => vrot%soilcmas
        extnprob          => vrot%extnprob
        prbfrhuc          => vrot%prbfrhuc
        mlightng          => vrot%mlightng
        slopefrac         => vrot%slopefrac
        stdaln            => vrot%stdaln
        lfstatusrow       => vrot%lfstatus
        pandaysrow        => vrot%pandays
        ipeatlandrow      => vrot%ipeatland
        Cmossmas          => vrot%Cmossmas
        litrmsmoss        => vrot%litrmsmoss
        dmoss             => vrot%dmoss

        grclarea          => vgat%grclarea

        FCANROT           => class_rot%FCANROT
        FAREROT           => class_rot%FAREROT
        RSMNROT           => class_rot%RSMNROT
        QA50ROT           => class_rot%QA50ROT
        VPDAROT           => class_rot%VPDAROT
        VPDBROT           => class_rot%VPDBROT
        PSGAROT           => class_rot%PSGAROT
        PSGBROT           => class_rot%PSGBROT
        DRNROT            => class_rot%DRNROT
        SDEPROT           => class_rot%SDEPROT
        XSLPROT           => class_rot%XSLPROT
        GRKFROT           => class_rot%GRKFROT
        WFSFROT           => class_rot%WFSFROT
        WFCIROT           => class_rot%WFCIROT
        MIDROT            => class_rot%MIDROT
        DELZ              => class_gat%DELZ
        ZBOT              => class_gat%ZBOT
        SANDROT           => class_rot%SANDROT
        CLAYROT           => class_rot%CLAYROT
        ORGMROT           => class_rot%ORGMROT
        TBARROT           => class_rot%TBARROT
        THLQROT           => class_rot%THLQROT
        THICROT           => class_rot%THICROT
        TCANROT           => class_rot%TCANROT
        TSNOROT           => class_rot%TSNOROT
        TPNDROT           => class_rot%TPNDROT
        ZPNDROT           => class_rot%ZPNDROT
        RCANROT           => class_rot%RCANROT
        SCANROT           => class_rot%SCANROT
        SNOROT            => class_rot%SNOROT
        ALBSROT           => class_rot%ALBSROT
        RHOSROT           => class_rot%RHOSROT
        GROROT            => class_rot%GROROT
        ZRFHROW           => class_rot%ZRFHROW
        ZRFMROW           => class_rot%ZRFMROW
        GCROW             => class_rot%GCROW
        ZBLDROW           => class_rot%ZBLDROW
        ALVCROT           => class_rot%ALVCROT
        ALICROT           => class_rot%ALICROT
        PAMNROT           => class_rot%PAMNROT
        PAMXROT           => class_rot%PAMXROT
        LNZ0ROT           => class_rot%LNZ0ROT
        CMASROT           => class_rot%CMASROT
        ROOTROT           => class_rot%ROOTROT
        DLATROW           => class_rot%DLATROW
        DLONROW           => class_rot%DLONROW
        RADJROW           => class_rot%RADJROW
        Z0ORROW           => class_rot%Z0ORROW
        GGEOROW           => class_rot%GGEOROW
        SOCIROT           => class_rot%SOCIROT
        TBASROT           => class_rot%TBASROT

        ! ----------------------------

        do i = 1, nlat
         RADJROW(i)=DLATROW(i)*PI/180.
         Z0ORROW(i)=0.0
         GGEOROW(i)=0.0
        end do

    !> ZRFMROW and ZRFHROW, the reference heights at which the momentum variables (wind speed) and energy variables
    !> (temperature and specific humidity) are provided.  In a run using atmospheric model forcing data, these heights
    !> would vary by time step, but since this version of the driver is set up to use field data, ZRFMROW and ZRFHROW
    !> refer to the measurement height of these variables, which is fixed.

        ZRFMROW = ncGet1DVar(initid, 'ZRFM', start = [lonIndex, latIndex], count = [1, 1])
        ZRFHROW = ncGet1DVar(initid, 'ZRFH', start = [lonIndex, latIndex], count = [1, 1])

    !> ZBLDROW, the atmospheric blending height.  Technically this variable depends on the length scale of the
    !> patches of roughness elements on the land surface, but this is difficult to ascertain.  Usually it is assigned a value of 50 m.

        ZBLDROW = ncGet1DVar(initid, 'ZBLD', start = [lonIndex, latIndex], count = [1, 1])

    !> GCROW, the GCM surface descriptor variable.  For land surfaces (including inland water) it has a value of -1.

        GCROW = ncGet1DVar(initid, 'GC', start = [lonIndex, latIndex], count = [1, 1])
        DRNROT = ncGet2DVar(initid, 'DRN', start = [lonIndex, latIndex, 1], count = [1, 1, nmos], format = [nlat, nmos])
        SDEPROT = ncGet2DVar(initid, 'SDEP', start = [lonIndex, latIndex, 1], count = [1, 1, nmos], format = [nlat, nmos])
        SOCIROT = ncGet2DVar(initid, 'SOCI', start = [lonIndex, latIndex, 1], count = [1, 1, nmos], format = [nlat, nmos])
        FAREROT = ncGet2DVar(initid, 'FARE', start = [lonIndex, latIndex, 1], count = [1, 1, nmos], format = [nlat, nmos])
        XSLPROT = ncGet2DVar(initid, 'XSLP', start = [lonIndex, latIndex, 1], count = [1, 1, nmos], format = [nlat, nmos])
        GRKFROT = ncGet2DVar(initid, 'GRKF', start = [lonIndex, latIndex, 1], count = [1, 1, nmos], format = [nlat, nmos])
        WFSFROT = ncGet2DVar(initid, 'WFSF', start = [lonIndex, latIndex, 1], count = [1, 1, nmos], format = [nlat, nmos])
        WFCIROT = ncGet2DVar(initid, 'WFCI', start = [lonIndex, latIndex, 1], count = [1, 1, nmos], format = [nlat, nmos])
        TCANROT = ncGet2DVar(initid, 'TCAN', start = [lonIndex, latIndex, 1], count = [1, 1, nmos], format = [nlat, nmos])
        TSNOROT = ncGet2DVar(initid, 'TSNO', start = [lonIndex, latIndex, 1], count = [1, 1, nmos], format = [nlat, nmos])
        TPNDROT = ncGet2DVar(initid, 'TPND', start = [lonIndex, latIndex, 1], count = [1, 1, nmos], format = [nlat, nmos])
        ZPNDROT = ncGet2DVar(initid, 'ZPND', start = [lonIndex, latIndex, 1], count = [1, 1, nmos], format = [nlat, nmos])
        RCANROT = ncGet2DVar(initid, 'RCAN', start = [lonIndex, latIndex, 1], count = [1, 1, nmos], format = [nlat, nmos])
        SCANROT = ncGet2DVar(initid, 'SCAN', start = [lonIndex, latIndex, 1], count = [1, 1, nmos], format = [nlat, nmos])
        SNOROT = ncGet2DVar(initid, 'SNO', start = [lonIndex, latIndex, 1], count = [1, 1, nmos], format = [nlat, nmos])
        ALBSROT = ncGet2DVar(initid, 'ALBS', start = [lonIndex, latIndex, 1], count = [1, 1, nmos], format = [nlat, nmos])
        RHOSROT = ncGet2DVar(initid, 'RHOS', start = [lonIndex, latIndex, 1], count = [1, 1, nmos], format = [nlat, nmos])
        GROROT = ncGet2DVar(initid, 'GRO', start = [lonIndex, latIndex, 1], count = [1, 1, nmos], format = [nlat, nmos])
        MIDROT = ncGet2DVar(initid, 'MID', start = [lonIndex, latIndex, 1], count = [1, 1, nmos], format=[nlat, nmos])
        LNZ0ROT = ncGet3DVar(initid, 'LNZ0', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icp1, nmos], format = [nlat, nmos, icp1])
        ALVCROT = ncGet3DVar(initid, 'ALVC', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icp1, nmos], format = [nlat, nmos, icp1])
        ALICROT = ncGet3DVar(initid, 'ALIC', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icp1, nmos], format = [nlat, nmos, icp1])
        PAMNROT = ncGet3DVar(initid, 'PAMN', start = [lonIndex, latIndex, 1, 1], count = [1, 1, ican, nmos], format = [nlat, nmos, ican])
        PAMXROT = ncGet3DVar(initid, 'PAMX', start = [lonIndex, latIndex, 1, 1], count = [1, 1, ican, nmos], format = [nlat, nmos, ican])
        CMASROT = ncGet3DVar(initid, 'CMAS', start = [lonIndex, latIndex, 1, 1], count = [1, 1, ican, nmos], format = [nlat, nmos, ican])
        ROOTROT = ncGet3DVar(initid, 'ROOT', start = [lonIndex, latIndex, 1, 1], count = [1, 1, ican, nmos], format = [nlat, nmos, ican])
        RSMNROT = ncGet3DVar(initid, 'RSMN', start = [lonIndex, latIndex, 1, 1], count = [1, 1, ican, nmos], format = [nlat, nmos, ican])
        QA50ROT = ncGet3DVar(initid, 'QA50', start = [lonIndex, latIndex, 1, 1], count = [1, 1, ican, nmos], format = [nlat, nmos, ican])
        VPDAROT = ncGet3DVar(initid, 'VPDA', start = [lonIndex, latIndex, 1, 1], count = [1, 1, ican, nmos], format = [nlat, nmos, ican])
        VPDBROT = ncGet3DVar(initid, 'VPDB', start = [lonIndex, latIndex, 1, 1], count = [1, 1, ican, nmos], format = [nlat, nmos, ican])
        PSGAROT = ncGet3DVar(initid, 'PSGA', start = [lonIndex, latIndex, 1, 1], count = [1, 1, ican, nmos], format = [nlat, nmos, ican])
        PSGBROT = ncGet3DVar(initid, 'PSGB', start = [lonIndex, latIndex, 1, 1], count = [1, 1, ican, nmos], format = [nlat, nmos, ican])
        SANDROT = ncGet3DVar(initid, 'SAND', start = [lonIndex, latIndex, 1, 1], count = [1, 1, ignd, nmos], format = [nlat, nmos, ignd])
        CLAYROT = ncGet3DVar(initid, 'CLAY', start = [lonIndex, latIndex, 1, 1], count = [1, 1, ignd, nmos], format = [nlat, nmos, ignd])
        ORGMROT = ncGet3DVar(initid, 'ORGM', start = [lonIndex, latIndex, 1, 1], count = [1, 1, ignd, nmos], format = [nlat, nmos, ignd])
        TBARROT = ncGet3DVar(initid, 'TBAR', start = [lonIndex, latIndex, 1, 1], count = [1, 1, ignd, nmos], format = [nlat, nmos, ignd])
        THLQROT = ncGet3DVar(initid, 'THLQ', start = [lonIndex, latIndex, 1, 1], count = [1, 1, ignd, nmos], format = [nlat, nmos, ignd])
        THICROT = ncGet3DVar(initid, 'THIC', start = [lonIndex, latIndex, 1, 1], count = [1, 1, ignd, nmos], format = [nlat, nmos, ignd])
        ZBOT = reshape(ncGet3DVar(initid, 'ZBOT', start = [lonIndex, latIndex, 1, 1], count = [1, 1, ignd, 1], format = [1, 1, ignd]), [ignd])
        DELZ = reshape(ncGet3DVar(initid, 'DELZ', start = [lonIndex, latIndex, 1, 1], count = [1, 1, ignd, 1], format = [1, 1, ignd]), [ignd])
        ipeatlandrow = ncGet2DVar(initid, 'ipeatland', start = [lonIndex, latIndex, 1], count = [1, 1, nmos], format = [nlat, nmos])

!         SCANROT=0. !FLAG test! This is just until we have a 'normal' initfile... JM Oct30 2017.
!         RCANROT=0.
!         RHOSROT=0.
!         ALBSROT = 0.
!         SNOROT = 0.
!         TBARROT = 258.
!         THICROT = 0.
!         TCANROT = 258.
!
!

        if (.not. ctem_on) then
            FCANROT = ncGet3DVar(initid, 'FCAN', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icp1, nmos], format = [nlat, nmos, icp1])
            ! Error check:
            do i = 1,nlat
                do m = 1,nmos
                    if (FAREROT(i,m) .gt. 1.0) then
                        print *,'FAREROT > 1',FAREROT(I,M)
                        call XIT('read_initialstate', -1)
                    end if
                enddo
            enddo
            !else fcancmx is read in instead and fcanrot is derived later.
        end if

!     Complete some initial set up work:

        DO 100 I=1,nlat
            DO 100 M=1,nmos

                DO J=1,IGND
                    TBARROT(I,M,J)=TBARROT(I,M,J)+TFREZ
                ENDDO
                TSNOROT(I,M)=TSNOROT(I,M)+TFREZ
                TCANROT(I,M)=TCANROT(I,M)+TFREZ

                TPNDROT(I,M)=TPNDROT(I,M)+TFREZ
                TBASROT(I,M)=TBARROT(I,M,IGND)
                !CMAIROT(I,M)=0.
                !WSNOROT(I,M)=0.
                !ZSNLROT(I,M)=0.10
                !TSFSROT(I,M,1)=TFREZ
                !TSFSROT(I,M,2)=TFREZ

                !TSFSROT(I,M,3)=TBARROT(I,M,1)
                !TSFSROT(I,M,4)=TBARROT(I,M,1)
                !TACROT (I,M)=TCANROT(I,M)
                !QACROT (I,M)=0.5E-2

                !DO 75 K=1,6
                !    DO 75 L=1,50
                !        ITCTROT(I,M,K,L)=0
!75              CONTINUE
100     CONTINUE

        ! Check that the THIC and THLQ values are set to zero for soil layers
        ! that are non-permeable (bedrock).
        do i = 1,nlat
            do j = 1,nmos
                do m = 1,ignd-1
                    if (zbot(m) > SDEPROT(i,j) .and. zbot(m+1) > SDEPROT(i,j)) then
                        THLQROT(i,j,m:ignd) = 0.
                        THICROT(i,j,m:ignd) = 0.
                        exit
                    end if
                end do
            end do
        end do

        if (ctem_on) then

            grclarea = ncGet1DVar(initid, 'grclarea', start = [lonIndex, latIndex], count = [1, 1])
            extnprob(:,1) = ncGet1DVar(initid, 'extnprob', start = [lonIndex, latIndex], count = [1, 1])
            prbfrhuc(:,1) = ncGet1DVar(initid, 'prbfrhuc', start = [lonIndex, latIndex], count = [1, 1])

            do i = 1,nmos
                grclarea(i) = grclarea(1)  !grclarea is ilg, but offline nlat is always 1 so ilg = nmos.
                extnprob(:,i) = extnprob(:,1)
                prbfrhuc(:,i) = prbfrhuc(:,1)
            end do

            mlightng(:,1,:) = ncGet2DVar(initid, 'mlightng', start = [lonIndex, latIndex, 1], count = [1, 1, nmos], format = [nlat, nmos])

            do i = 1, nmos
                mlightng(:,i,:) = mlightng(:,1,:)
            end do

            slopefrac = ncGet3DVar(initid, 'slopefrac', start = [lonIndex, latIndex, 1, 1], count = [1, 1, nmos, 8], format = [nlat, nmos, 8])
            Cmossmas = ncGet2DVar(initid, 'Cmossmas', start = [lonIndex, latIndex, 1], count = [1, 1, nmos], format = [nlat, nmos])
            litrmsmoss = ncGet2DVar(initid, 'litrmsmoss', start = [lonIndex, latIndex, 1], count = [1, 1, nmos], format = [nlat, nmos])
            dmoss = ncGet2DVar(initid, 'dmoss', start = [lonIndex, latIndex, 1], count = [1, 1, nmos], format = [nlat, nmos])
            ailcminrow = ncGet3DVar(initid, 'ailcmin', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos], format = [nlat, nmos,icc])
            ailcmaxrow = ncGet3DVar(initid, 'ailcmax', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos], format = [nlat, nmos,icc])
            fcancmxrow = ncGet3DVar(initid, 'fcancmx', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos], format = [nlat, nmos,icc])

            !>Rest of the initialization variables are needed to run CTEM but if starting from bare ground initialize all
            !>live and dead c pools from zero. suitable values of extnprobgrd and prbfrhucgrd would still be required. set
            !>stdaln to 1 for operation in non-gcm stand alone mode, in the CTEM initialization file.

            gleafmasrow = ncGet3DVar(initid, 'gleafmas', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos], format = [nlat, nmos,icc])
            bleafmasrow = ncGet3DVar(initid, 'bleafmas', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos], format = [nlat, nmos,icc])
            stemmassrow = ncGet3DVar(initid, 'stemmass', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos], format = [nlat, nmos,icc])
            rootmassrow = ncGet3DVar(initid, 'rootmass', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos], format = [nlat, nmos,icc])

            !>If fire and competition are on, save the stemmass and rootmass for use in burntobare subroutine on the first timestep.
            if (dofire .and. PFTCompetition) then
                        do j =1,icc
                            pstemmassrow(i,m,j)=stemmassrow(i,m,j)
                            pgleafmassrow(i,m,j)=rootmassrow(i,m,j)
                        end do
            end if

            !litrmassrow = ncGet4DVar(initid, 'litrmass', start = [lonIndex, latIndex, 1, 1, 1], count = [1, 1, iccp1, ignd, nmos], format = [nlat, nmos, iccp1, ignd])
            !soilcmasrow = ncGet4DVar(initid, 'soilcmas', start = [lonIndex, latIndex, 1, 1,1], count = [1, 1, iccp1, ignd,nmos], format = [nlat, nmos,iccp1, ignd])
            litrmassrow = ncGet3DVar(initid, 'litrmass', start = [lonIndex, latIndex, 1, 1, 1], count = [1, 1, iccp1, nmos], format = [nlat, nmos, iccp1])
            soilcmasrow = ncGet3DVar(initid, 'soilcmas', start = [lonIndex, latIndex, 1, 1,1], count = [1, 1, iccp1, nmos], format = [nlat, nmos,iccp1])
            lfstatusrow = ncGet3DVar(initid, 'lfstatus', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos], format = [nlat, nmos,icc])
            pandaysrow = ncGet3DVar(initid, 'pandays', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos], format = [nlat, nmos,icc])

            if (PFTCompetition .and. inibioclim) then  !read in the bioclimatic parameters

                twarmm(:,1) = ncGet1DVar(initid, 'twarmm', start = [lonIndex, latIndex], count = [1, 1])!, format = [nlat])
                tcoldm(:,1) = ncGet1DVar(initid, 'tcoldm', start = [lonIndex, latIndex], count = [1, 1])!, format = [nlat])
                gdd5(:,1) = ncGet1DVar(initid, 'gdd5', start = [lonIndex, latIndex], count = [1, 1])!, format = [nlat])
                aridity(:,1) = ncGet1DVar(initid, 'aridity', start = [lonIndex, latIndex], count = [1, 1])!, format = [nlat])
                srplsmon(:,1) = ncGet1DVar(initid, 'srplsmon', start = [lonIndex, latIndex], count = [1, 1])!, format = [nlat])
                defctmon(:,1) = ncGet1DVar(initid, 'defctmon', start = [lonIndex, latIndex], count = [1, 1])!, format = [nlat])
                anndefct(:,1) = ncGet1DVar(initid, 'anndefct', start = [lonIndex, latIndex], count = [1, 1])!, format = [nlat])
                annsrpls(:,1) = ncGet1DVar(initid, 'annsrpls', start = [lonIndex, latIndex], count = [1, 1])!, format = [nlat])
                annpcp(:,1) = ncGet1DVar(initid, 'annpcp', start = [lonIndex, latIndex], count = [1, 1])!, format = [nlat])
                dry_season_length(:,1) = ncGet1DVar(initid, 'dry_season_length', start = [lonIndex, latIndex], count = [1, 1])!, format = [nlat])

                !>Take the first tile value now and put it over the other tiles
                do m = 1,nmos
                    twarmm(:,m)=twarmm(:,1)
                    tcoldm(:,m)=tcoldm(:,1)
                    gdd5(:,m)=gdd5(:,1)
                    aridity(:,m)=aridity(:,1)
                    srplsmon(:,m)=srplsmon(:,1)
                    defctmon(:,m)=defctmon(:,1)
                    anndefct(:,m)=anndefct(:,1)
                    annsrpls(:,m)=annsrpls(:,1)
                    annpcp(:,m)=annpcp(:,1)
                    dry_season_length(:,m) =dry_season_length(:,1)
                end do

            else if (PFTCompetition .and. .not. inibioclim) then ! set them to zero

                twarmm=0.0
                tcoldm=0.0
                gdd5=0.0
                aridity=0.0
                srplsmon=0.0
                defctmon=0.0
                anndefct=0.0
                annsrpls=0.0
                annpcp=0.0
                dry_season_length = 0.0

            endif

            !>if this run uses the competition and starts from bare ground, set up the model state here. this
            !>overwrites what was read in from the initialization file.
            if (PFTCompetition .and. start_bare) then

                do i=1,nlat
                    do m = 1,nmos

                        do j = 1,icc
                            if (.not. crop(j)) fcancmxrow(i,m,j) = 0.0
                            ailcminrow(i,m,j)=0.0
                            ailcmaxrow(i,m,j)=0.0
                            gleafmasrow(i,m,j)=0.0
                            bleafmasrow(i,m,j)=0.0
                            stemmassrow(i,m,j)=0.0
                            rootmassrow(i,m,j)=0.0
                            lfstatusrow(i,m,j)=4
                            pandaysrow(i,m,j)=0
                        enddo

                        lfstatusrow(i,m,1)=2

                        do j = 1,iccp1
                            litrmassrow(i,m,j)=0.0
                            soilcmasrow(i,m,j)=0.0
                        enddo
                    end do ! nmtest
                enddo !nltest

            end if !if (PFTCompetition .and. start_bare)

        end if !ctem_on

    end subroutine read_initialstate

    !>@}
    ! ------------------------------------------------------------------------------------

    !>\ingroup model_state_drivers_write_restart
    !!@{
    !> Write out the model restart file to netcdf. We only write out the variables that the model
    !! influences. This overwrites a pre-existing netcdf file.

    subroutine write_restart(lonIndex,latIndex)

        ! J. Melton
        ! Jun 2017

        use ctem_statevars,     only : c_switch,vrot
        use class_statevars,    only : class_rot
        use ctem_params,        only : icc,nmos,ignd,icp1,nlat,ican,l2max,modelpft

        implicit none

        ! arguments
        integer, intent(in) :: lonIndex,latIndex

        ! pointers:
        real, pointer, dimension(:,:,:) :: FCANROT
        real, pointer, dimension(:,:)   :: FAREROT
        real, pointer, dimension(:,:,:) :: TBARROT
        real, pointer, dimension(:,:,:) :: THLQROT
        real, pointer, dimension(:,:,:) :: THICROT
        real, pointer, dimension(:,:)   :: TCANROT
        real, pointer, dimension(:,:)   :: TSNOROT
        real, pointer, dimension(:,:)   :: TPNDROT
        real, pointer, dimension(:,:)   :: ZPNDROT
        real, pointer, dimension(:,:)   :: RCANROT
        real, pointer, dimension(:,:)   :: SCANROT
        real, pointer, dimension(:,:)   :: SNOROT
        real, pointer, dimension(:,:)   :: ALBSROT
        real, pointer, dimension(:,:)   :: RHOSROT
        real, pointer, dimension(:,:)   :: GROROT

        logical, pointer :: ctem_on
        logical, pointer :: PFTCompetition
        logical, pointer :: lnduseon
        real, pointer, dimension(:,:,:) :: ailcminrow           !
        real, pointer, dimension(:,:,:) :: ailcmaxrow           !
        real, pointer, dimension(:,:,:) :: fcancmxrow           !
        real, pointer, dimension(:,:,:) :: gleafmasrow          !
        real, pointer, dimension(:,:,:) :: bleafmasrow          !
        real, pointer, dimension(:,:,:) :: stemmassrow          !
        real, pointer, dimension(:,:,:) :: rootmassrow          !
        real, pointer, dimension(:,:) :: twarmm            !< temperature of the warmest month (c)
        real, pointer, dimension(:,:) :: tcoldm            !< temperature of the coldest month (c)
        real, pointer, dimension(:,:) :: gdd5              !< growing degree days above 5 c
        real, pointer, dimension(:,:) :: aridity           !< aridity index, ratio of potential evaporation to precipitation
        real, pointer, dimension(:,:) :: srplsmon          !< number of months in a year with surplus water i.e.precipitation more than potential evaporation
        real, pointer, dimension(:,:) :: defctmon          !< number of months in a year with water deficit i.e.precipitation less than potential evaporation
        real, pointer, dimension(:,:) :: anndefct          !< annual water deficit (mm)
        real, pointer, dimension(:,:) :: annsrpls          !< annual water surplus (mm)
        real, pointer, dimension(:,:) :: annpcp            !< annual precipitation (mm)
        real, pointer, dimension(:,:) :: dry_season_length !< length of dry season (months)
        real, pointer, dimension(:,:,:) :: litrmassrow
        real, pointer, dimension(:,:,:) :: soilcmasrow
        integer, pointer, dimension(:,:,:) :: lfstatusrow
        integer, pointer, dimension(:,:,:) :: pandaysrow
        real, pointer, dimension(:,:) :: Cmossmas          !<C in moss biomass, \f$kg C/m^2\f$
        real, pointer, dimension(:,:) :: litrmsmoss        !<moss litter mass, \f$kg C/m^2\f$
        real, pointer, dimension(:,:) :: dmoss             !<depth of living moss (m)

        ! local variables
        !integer :: i,m,j,k1c,k2c,n
        !integer, dimension(nlat,nmos) :: icountrow
        !real, dimension(icc) :: rnded_pft
        real, parameter :: TFREZ = 273.16

        ! point pointers:
        ctem_on           => c_switch%ctem_on
        PFTCompetition    => c_switch%PFTCompetition
        lnduseon          => c_switch%lnduseon
        ailcminrow        => vrot%ailcmin
        ailcmaxrow        => vrot%ailcmax
        fcancmxrow        => vrot%fcancmx
        gleafmasrow       => vrot%gleafmas
        bleafmasrow       => vrot%bleafmas
        stemmassrow       => vrot%stemmass
        rootmassrow       => vrot%rootmass
        twarmm            => vrot%twarmm
        tcoldm            => vrot%tcoldm
        gdd5              => vrot%gdd5
        aridity           => vrot%aridity
        srplsmon          => vrot%srplsmon
        defctmon          => vrot%defctmon
        anndefct          => vrot%anndefct
        annsrpls          => vrot%annsrpls
        annpcp            => vrot%annpcp
        dry_season_length => vrot%dry_season_length
        litrmassrow       => vrot%litrmass
        soilcmasrow       => vrot%soilcmas
        lfstatusrow       => vrot%lfstatus
        pandaysrow        => vrot%pandays
        Cmossmas          => vrot%Cmossmas
        litrmsmoss        => vrot%litrmsmoss
        dmoss             => vrot%dmoss
        FCANROT           => class_rot%FCANROT
        FAREROT           => class_rot%FAREROT
        TBARROT           => class_rot%TBARROT
        THLQROT           => class_rot%THLQROT
        THICROT           => class_rot%THICROT
        TCANROT           => class_rot%TCANROT
        TSNOROT           => class_rot%TSNOROT
        TPNDROT           => class_rot%TPNDROT
        ZPNDROT           => class_rot%ZPNDROT
        RCANROT           => class_rot%RCANROT
        SCANROT           => class_rot%SCANROT
        SNOROT            => class_rot%SNOROT
        ALBSROT           => class_rot%ALBSROT
        RHOSROT           => class_rot%RHOSROT
        GROROT            => class_rot%GROROT

        call ncPut2DVar(rsid, 'FARE',   realValues = FAREROT, start = [lonIndex, latIndex, 1], count = [1, 1, nmos])
        call ncPut2DVar(rsid, 'TCAN',   realValues = TCANROT-TFREZ, start = [lonIndex, latIndex, 1], count = [1, 1, nmos])
        call ncPut2DVar(rsid, 'TSNO',   realValues = TSNOROT-TFREZ, start = [lonIndex, latIndex, 1], count = [1, 1, nmos])
        call ncPut2DVar(rsid, 'TPND',   realValues = TPNDROT-TFREZ, start = [lonIndex, latIndex, 1], count = [1, 1, nmos])
        call ncPut2DVar(rsid, 'ZPND',   realValues = ZPNDROT, start = [lonIndex, latIndex, 1], count = [1, 1, nmos])
        call ncPut2DVar(rsid, 'RCAN',   realValues = RCANROT, start = [lonIndex, latIndex, 1], count = [1, 1, nmos])
        call ncPut2DVar(rsid, 'SCAN',   realValues = SCANROT, start = [lonIndex, latIndex, 1], count = [1, 1, nmos])
        call ncPut2DVar(rsid, 'SNO',    realValues = SNOROT, start = [lonIndex, latIndex, 1], count = [1, 1, nmos])
        call ncPut2DVar(rsid, 'ALBS',   realValues = ALBSROT, start = [lonIndex, latIndex, 1], count = [1, 1, nmos])
        call ncPut2DVar(rsid, 'RHOS',   realValues = RHOSROT, start = [lonIndex, latIndex, 1], count = [1, 1, nmos])
        call ncPut2DVar(rsid, 'GRO',    realValues = GROROT, start = [lonIndex, latIndex, 1], count = [1, 1, nmos])

        call ncPut3DVar(rsid, 'FCAN',   realValues = FCANROT, start = [lonIndex, latIndex, 1, 1], count = [1, 1, icp1, nmos])
        call ncPut3DVar(rsid, 'THLQ',   realValues = THLQROT, start = [lonIndex, latIndex, 1, 1], count = [1, 1, ignd, nmos])
        call ncPut3DVar(rsid, 'THIC',   realValues = THICROT, start = [lonIndex, latIndex, 1, 1], count = [1, 1, ignd, nmos])
        call ncPut3DVar(rsid, 'TBAR',   realValues = TBARROT-TFREZ, start = [lonIndex, latIndex, 1, 1], count = [1, 1, ignd, nmos])

        if (ctem_on) then

            call ncPut2DVar(rsid, 'Cmossmas',   realValues = Cmossmas, start = [lonIndex, latIndex, 1], count = [1, 1, nmos])
            call ncPut2DVar(rsid, 'litrmsmoss', realValues = litrmsmoss, start = [lonIndex, latIndex, 1], count = [1, 1, nmos])
            call ncPut2DVar(rsid, 'dmoss',      realValues = dmoss, start = [lonIndex, latIndex, 1], count = [1, 1, nmos])

            call ncPut3DVar(rsid, 'ailcmin',    realValues = ailcminrow, start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos])
            call ncPut3DVar(rsid, 'ailcmax',    realValues = ailcmaxrow, start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos])
            call ncPut3DVar(rsid, 'fcancmx',    realValues = fcancmxrow, start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos])
            call ncPut3DVar(rsid, 'gleafmas',   realValues = gleafmasrow, start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos])
            call ncPut3DVar(rsid, 'bleafmas',   realValues = bleafmasrow, start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos])
            call ncPut3DVar(rsid, 'stemmass',   realValues = stemmassrow, start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos])
            call ncPut3DVar(rsid, 'rootmass',   realValues = rootmassrow, start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos])
            call ncPut3DVar(rsid, 'litrmass',   realValues = litrmassrow, start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos])
            call ncPut3DVar(rsid, 'soilcmas',   realValues = soilcmasrow, start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos])
            call ncPut3DVar(rsid, 'lfstatus',   intValues = lfstatusrow, start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos])
            call ncPut3DVar(rsid, 'pandays',    intValues = pandaysrow, start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos])

            if (PFTCompetition) then

                call ncPutVar(rsid, 'twarmm', realValues = reshape(twarmm(1:1,1:1), [1]), start = [lonIndex, latIndex], count = [1, 1])
! JOE
! After checking that writing twarmm works as desired, please uncomment the lines below
!                call ncPutVar(rsid, 'tcoldm', tcoldm, start = [lonIndex, latIndex])
!                call ncPutVar(rsid, 'gdd5', gdd5, start = [lonIndex, latIndex])
!                call ncPutVar(rsid, 'aridity', aridity, start = [lonIndex, latIndex])
!                call ncPutVar(rsid, 'srplsmon', srplsmon, start = [lonIndex, latIndex])
!                call ncPutVar(rsid, 'defctmon', defctmon, start = [lonIndex, latIndex])
!                call ncPutVar(rsid, 'anndefct', anndefct, start = [lonIndex, latIndex])
!                call ncPutVar(rsid, 'annsrpls', annsrpls, start = [lonIndex, latIndex])
!                call ncPutVar(rsid, 'annpcp', annpcp, start = [lonIndex, latIndex])
!                call ncPutVar(rsid, 'dry_season_length', dry_season_length, start = [lonIndex, latIndex])

            end if ! PFTCompetition

        end if !ctem_on

    end subroutine write_restart

    !>@}
    ! ------------------------------------------------------------------------------------

    !>\ingroup model_state_drivers_getInput
    !!@{
    !>  Read in a model input from a netcdf file and store the file's time array
    !! as well as the input values into memory.

    subroutine getInput(inputRequested,longitude,latitude)



        use fileIOModule
        use generalUtils, only : parseTimeStamp
        use ctem_statevars, only : c_switch,vrot
        use ctem_params, only : icc,nmos
        use outputManager, only : checkForTime

        implicit none

        character(*), intent(in) :: inputRequested
        real, intent(in), optional :: longitude
        real, intent(in), optional :: latitude
        integer :: lengthOfFile
        integer :: lonloc,latloc
        integer :: i,arrindex,m,numPFTsinFile
        real, dimension(:), allocatable :: fileTime
        logical, pointer :: transientCO2
        integer, pointer :: fixedYearCO2
        logical, pointer :: transientCH4
        integer, pointer :: fixedYearCH4
        logical, pointer :: transientPOPD
        integer, pointer :: fixedYearPOPD
        logical, pointer :: transientLGHT
        integer, pointer :: fixedYearLGHT
        logical, pointer :: lnduseon
        integer, pointer :: fixedYearLUC

        real, dimension(5) :: dateTime
        real :: startLGHTTime
        real, pointer, dimension(:,:) :: co2concrow
        real, pointer, dimension(:,:) :: ch4concrow
        real, pointer, dimension(:,:) :: popdinrow
        real, pointer, dimension(:,:,:) :: fcancmxrow
        real, pointer, dimension(:,:,:) :: mlightngrow

        transientCO2    => c_switch%transientCO2
        fixedYearCO2    => c_switch%fixedYearCO2
        transientCH4    => c_switch%transientCH4
        fixedYearCH4    => c_switch%fixedYearCH4
        transientPOPD   => c_switch%transientPOPD
        fixedYearPOPD   => c_switch%fixedYearPOPD
        transientLGHT   => c_switch%transientLGHT
        fixedYearLGHT   => c_switch%fixedYearLGHT
        lnduseon        => c_switch%lnduseon
        fixedYearLUC    => c_switch%fixedYearLUC
        co2concrow      => vrot%co2conc
        ch4concrow      => vrot%ch4conc
        popdinrow       => vrot%popdin
        fcancmxrow      => vrot%fcancmx
        mlightngrow       => vrot%mlightng

        select case (trim(inputRequested))

        !> For each of the time varying inputs in this subroutine, we take in the whole dataset
        !! and later determine the year we need (in updateInput). The general approach is that these
        !! files are light enough on memory demands to make this acceptable.

        !! It is important that the files have time as the fastest varying dimension.

        case ('CO2') ! Carbon dioxide concentration

            lengthOfFile = ncGetDimLen(co2id, 'time')
            allocate(fileTime(lengthOfFile))
            allocate(CO2Time(lengthOfFile))

            fileTime = ncGet1DVar(CO2id, 'time', start = [1], count = [lengthOfFile])

            ! Parse these into just years (expected format is "day as %Y%m%d.%f")
            do i = 1, lengthOfFile
                dateTime = parseTimeStamp(fileTime(i))
                CO2Time(i) = dateTime(1) ! Rewrite putting in the year
            end do

            if (transientCO2) then
                ! We read in the whole CO2 times series and store it.
                allocate(CO2FromFile(lengthOfFile))
                CO2FromFile = ncGet1DVar(CO2id, 'mole_fraction_of_carbon_dioxide_in_air', start = [1], count = [lengthOfFile])
            else
                ! Find the requested year in the file.
                arrindex = checkForTime(lengthOfFile,real(CO2Time),real(fixedYearCO2))
                if (arrindex == 0) stop('getInput says: The CO2 file does not contain requested year')

                ! We read in only the suggested year
                i = 1 ! offline nlat is always 1 so just set
                co2concrow(i,:) = ncGet1DVar(CO2id, 'mole_fraction_of_carbon_dioxide_in_air', start = [arrindex], count = [1])
            end if

        case ('CH4') ! Methane concentration

            lengthOfFile = ncGetDimLen(ch4id, 'time')
            allocate(fileTime(lengthOfFile))
            allocate(CH4Time(lengthOfFile))

            fileTime = ncGet1DVar(ch4id, 'time', start = [1], count = [lengthOfFile])

            ! Parse these into just years (expected format is "day as %Y%m%d.%f")
            do i = 1, lengthOfFile
                dateTime = parseTimeStamp(fileTime(i))
                CH4Time(i) = dateTime(1) ! Rewrite putting in the year
            end do

            if (transientCH4) then
                ! We read in the whole CH3 times series and store it.
                allocate(CH4FromFile(lengthOfFile))
                CH4FromFile = ncGet1DVar(ch4id, 'mole_fraction_of_methane_in_air', start = [1], count = [lengthOfFile])
            else
                ! Find the requested year in the file.
                arrindex = checkForTime(lengthOfFile,real(CH4Time),real(fixedYearCH4))
                if (arrindex == 0) stop('getInput says: The CH4 file does not contain requested year')

                ! We read in only the suggested year
                i = 1 ! offline nlat is always 1 so just set
                ch4concrow(i,:) = ncGet1DVar(ch4id, 'mole_fraction_of_methane_in_air', start = [arrindex], count = [1])
            end if

        case ('POPD') ! Population density

            lengthOfFile = ncGetDimLen(popid, 'time')
            allocate(fileTime(lengthOfFile))
            allocate(POPDTime(lengthOfFile))

            fileTime = ncGet1DVar(popid, 'time', start = [1], count = [lengthOfFile])

            ! Parse these into just years (expected format is "day as %Y%m%d.%f")
            do i = 1, lengthOfFile
                dateTime = parseTimeStamp(fileTime(i))
                POPDTime(i) = dateTime(1) ! Rewrite putting in the year
            end do
            lonloc = closestCell(popid,'lon',longitude)
            latloc = closestCell(popid,'lat',latitude)

            if (transientPOPD) then
                ! We read in the whole POPD times series and store it.
                allocate(POPDFromFile(lengthOfFile))
                POPDFromFile = ncGet1DVar(popid, 'popd', start = [1,lonloc,latloc], count = [lengthOfFile,1,1])

            else
                ! Find the requested year in the file.
                arrindex = checkForTime(lengthOfFile,real(POPDTime),real(fixedYearPOPD))
                if (arrindex == 0) stop('getInput says: The POPD file does not contain requested year')

                ! We read in only the suggested year
                i = 1 ! offline nlat is always 1 so just set
                popdinrow(i,:) = ncGet1DVar(popid, 'popd', start = [arrindex,lonloc,latloc], count = [1,1,1])

            end if

         case ('LGHT') ! Lightning strikes

            lengthOfFile = ncGetDimLen(lghtid, 'time')
            allocate(fileTime(lengthOfFile))
            allocate(LGHTTime(lengthOfFile))

            fileTime = ncGet1DVar(lghtid, 'time', start = [1], count = [lengthOfFile])

            ! The lightning file is daily (expected format is "day as %Y%m%d.%f")
            ! We want to retain all except the partial day.
            do i = 1, lengthOfFile
                dateTime = parseTimeStamp(fileTime(i))
                LGHTTime(i) = dateTime(1) * 10000. + dateTime(2) * 100. + dateTime(3)
            end do
            lonloc = closestCell(lghtid,'lon',longitude)
            latloc = closestCell(lghtid,'lat',latitude)

            ! Units expected are "strikes km-2 yr-1"

            if (transientLGHT) then
                ! We read in the whole LGHT times series and store it.
                allocate(LGHTFromFile(lengthOfFile))
                LGHTFromFile = ncGet1DVar(lghtid, 'lght', start = [1,lonloc,latloc], count = [lengthOfFile,1,1])

            else
                ! Find the requested day and year in the file.
                ! Assume we are grabbing from day 1
                startLGHTTime = real(fixedYearLGHT) * 10000. + 1. * 100. + 1.
                arrindex = checkForTime(lengthOfFile,LGHTTime,startLGHTTime)
                if (arrindex == 0) stop('getInput says: The LGHT file does not contain requested year')

                ! We read in only the suggested year of daily inputs
                ! FLAG Not presently set up for leap years!
                allocate(LGHTFromFile(365))
                LGHTFromFile = ncGet1DVar(lghtid, 'lght', start = [arrindex,lonloc,latloc], count = [365,1,1])

            end if

        case ('LUC') ! Land use change

            lengthOfFile = ncGetDimLen(lucid, 'time')
            allocate(fileTime(lengthOfFile))
            allocate(LUCTime(lengthOfFile))

            fileTime = ncGet1DVar(lucid, 'time', start = [1], count = [lengthOfFile])

            ! Parse these into just years (expected format is "day as %Y%m%d.%f")
            do i = 1, lengthOfFile
                dateTime = parseTimeStamp(fileTime(i))
                LUCTime(i) = dateTime(1) ! Rewrite putting in only the year
            end do
            lonloc = closestCell(lucid,'lon',longitude)
            latloc = closestCell(lucid,'lat',latitude)

            ! Ensure the file has the expected number of PFTs
            numPFTsinFile = ncGetDimLen(lucid, 'lev')
            if (numPFTsinFile .ne. icc) stop('getInput says: LUC file does not have expected number of PFTs')

            if (lnduseon) then
                ! We read in the whole LUC times series and store it.
                allocate(LUCFromFile(lengthOfFile,icc))
                LUCFromFile = ncGet2DVar(lucid, 'frac', start = [1,1,lonloc,latloc], count = [lengthOfFile,icc,1,1])
            else
                ! Find the requested year in the file.
                arrindex = checkForTime(lengthOfFile,real(LUCTime),real(fixedYearLUC))
                if (arrindex == 0) stop('getInput says: The LUC file does not contain requested year')

                ! We read in only the suggested year
                i = 1 ! offline nlat is always 1 so just set
                m = 1 ! FLAG this is set up only for 1 tile at PRESENT! JM

                if (nmos .ne. 1) stop('getInput for LUC is not setup for more than one tile at present!')

                fcancmxrow(i,m,:) = ncGet1DVar(lucid, 'frac', start = [arrindex,1,lonloc,latloc], count = [1,icc,1,1])

            end if

        case ('OBSWETF') ! Observed wetland fractions

            stop('Not implemented yet')
            !wetfrac_monrow(i,m,:
            !  wetfrac_presgat(i)=wetfrac_mongat(i,month1)+(real(xday)/30.0)*&
            !            (wetfrac_mongat(i,month2)-wetfrac_mongat(i,month1))

        case default
            stop('Specify an input kind for getInput')

        end select

        deallocate(fileTime)

    end subroutine getInput

    !>@}
    ! ------------------------------------------------------------------------------------

    !>\ingroup model_state_drivers_updateInput
    !!@{
    !> Update the input field variable based on the present model timestep

    subroutine updateInput(inputRequested,iyear,imonth,iday)

        use outputManager, only : checkForTime
        use ctem_statevars, only : vrot,c_switch,vgat
        use ctem_params, only : nmos

        implicit none

        character(*), intent(in) :: inputRequested
        integer, intent(in) :: iyear
        integer, intent(in), optional :: imonth
        integer, intent(in), optional :: iday
        integer :: arrindex,lengthTime,i,m
        real :: LGHTTimeNow
        real, pointer, dimension(:,:) :: co2concrow
        real, pointer, dimension(:,:) :: ch4concrow
        real, pointer, dimension(:,:) :: popdinrow
        real, pointer, dimension(:,:,:) :: nfcancmxrow
        real, pointer, dimension(:) :: lightng       !<total \f$lightning, flashes/(km^2 . year)\f$ it is assumed that cloud
                                                     !<to ground lightning is some fixed fraction of total lightning.

        logical, pointer :: transientLGHT
        integer, pointer :: fixedYearLGHT

        co2concrow      => vrot%co2conc
        ch4concrow      => vrot%ch4conc
        popdinrow       => vrot%popdin
        nfcancmxrow     => vrot%nfcancmx
        transientLGHT   => c_switch%transientLGHT
        fixedYearLGHT   => c_switch%fixedYearLGHT
        lightng         => vgat%lightng

        select case (trim(inputRequested))

        case ('CO2')

            lengthTime = size(CO2Time)

            ! Find the requested year in the file.
            arrindex = checkForTime(lengthTime,real(CO2Time),real(iyear))
            i = 1 ! offline nlat is always 1 so just set
            co2concrow(i,:) = CO2FromFile(arrindex)

        case ('CH4')

            lengthTime = size(CH4Time)

            ! Find the requested year in the file.
            arrindex = checkForTime(lengthTime,real(CH4Time),real(iyear))
            i = 1 ! offline nlat is always 1 so just set
            ch4concrow(i,:) = CH4FromFile(arrindex)

        case ('POPD')

            lengthTime = size(POPDTime)

            ! Find the requested year in the file.
            arrindex = checkForTime(lengthTime,real(POPDTime),real(iyear))
            i = 1 ! offline nlat is always 1 so just set
            popdinrow(i,:) = POPDFromFile(arrindex)

        case ('LUC')

            lengthTime = size(LUCTime)

            ! Find the requested year in the file.
            arrindex = checkForTime(lengthTime,real(LUCTime),real(iyear))
            i = 1 ! offline nlat is always 1 so just set
            m = 1 ! FLAG this is set up only for 1 tile at PRESENT! JM
            if (nmos > 1) stop('updateInput for LUC only set up for 1 tile at present')
            nfcancmxrow(i,m,:) = LUCFromFile(arrindex,:)

       case('LGHT')

            ! This file is daily so we need to find the day we are looking for.
            ! imonth is starting at 0 so add 1 always.

            lengthTime = size(LGHTTime)

            if (transientLGHT) then
                LGHTTimeNow = real(iyear) * 10000. + real(imonth+1) * 100. + real(iday)
            else ! we only draw from the fixedYearLGHT
                LGHTTimeNow = fixedYearLGHT * 10000. + real(imonth+1) * 100. + real(iday)
            end if

            ! Find the requested year in the file.
            arrindex = checkForTime(lengthTime,LGHTTime,LGHTTimeNow)

            lightng(1)= LGHTFromFile(arrindex)

            ! Since lighning is the same for all tiles, and nlat is always 1 offline, then we
            ! can just pass the same values across all ilg.
            do m = 1, size(lightng)
                lightng(m) = lightng(1)
            end do

        case default
            stop('specify an input kind for updateInput')
        end select

    end subroutine updateInput

    !>@}
    ! ------------------------------------------------------------------------------------

    !>\ingroup model_state_drivers_getMet
    !!@{
    !> Read in the meteorological input from a netcdf file
    !! It is **very** important that the files have time as the fastest varying dimension.
    !! There is an orders of magnitude slow-up if the dimensions are out of order!
    subroutine getMet(longitude,latitude,nday,delt)

        use fileIOModule
        use ctem_statevars, only : c_switch
        use generalUtils, only : parseTimeStamp,closeEnough

        implicit none

        real, intent(in) :: longitude       !< Longitude of grid cell of interest
        real, intent(in) :: latitude        !< Latitude of grid cell of interest
        integer, intent(in) :: nday         !< Maximum number of physics timesteps in one day
        real, intent(in) :: delt            !< Physics timestep (s)

        integer, pointer :: readMetStartYear !< First year of meteorological forcing to read in from the met file
        integer, pointer :: readMetEndYear   !< Last year of meteorological forcing to read in from the met file

        real :: moStart,moEnd,domStart,domEnd !< Assumed start and end months and days of month
        real :: timeStart, timeEnd            !< Calculated start and end in the format:%Y%m%d.%f
        integer :: lengthOfFile
        integer :: lonloc,latloc,i
        real, dimension(:), allocatable :: fileTime
        real, dimension(:), allocatable :: tempTime
        integer :: validTimestep
        integer :: firstIndex
        real, dimension(5) :: firstTime,secondTime

        readMetStartYear  => c_switch%readMetStartYear
        readMetEndYear    => c_switch%readMetEndYear

        !! It is very important that the files have time as the fastest varying dimension.
        !! There is a orders of magnitude slow-up if the dimensions are out of order.

        ! Grab the length of time dimension from the SW met file and write it to an array.
        ! NOTE: We assume the user is careful enough to ensure the time array is the same
        ! across all met files!
        lengthOfFile = ncGetDimLen(metFssId, 'time')
        allocate(fileTime(lengthOfFile))
        fileTime = ncGet1DVar(metFssId, 'time', start = [1], count = [lengthOfFile])

        ! Construct the time bounds that we will look for in the file.
        ! We assume that you will start on the first timestep of the day.
        ! Further the default is to start on (or at least look for) Jan 1
        ! of the yrStart year.
        moStart=1.
        domStart=1.
        ! The first time is considered to be the first physics timestep so given a fractional day of 0.
        timeStart = readMetStartYear * 10000. + moStart * 100. + domStart
        moEnd=12.
        domEnd=31.
        ! The last time is considered to be the last physics timestep of the day
        timeEnd =  readMetEndYear * 10000. + moEnd * 100. + domEnd + (real(nday - 1) * delt / 86400.)

        ! Now we read in and append the metTime the timesteps from the time variable of the met file. This
        ! uses the intrinsic move_alloc, but it simply appends to the array.
        allocate(metTime(0))
        validTimestep=0
        firstIndex=999999999 ! set to large value

        do i = 1, lengthOfFile
            if (fileTime(i) >= timeStart .and. fileTime(i) <= timeEnd) then
                validTimestep = validTimestep + 1
                allocate(tempTime(validTimestep))
                tempTime(1:validTimestep - 1) = metTime(1 : validTimestep - 1)
                call move_alloc(tempTime,metTime)
                metTime(validTimestep) = fileTime(i)
                firstIndex = min(firstIndex,i)
            end if
        end do

        ! Check that the first day is Jan 1, otherwise warn the user
        firstTime =  parseTimeStamp(metTime(1))
        if (.not. closeEnough(firstTime(5),1.)) then
            print*,'Warning, your met file does not start on Jan 1.'
        end if

        ! Determine the time step of the met data and
        ! convert from fraction of day to period in seconds
        secondTime =  parseTimeStamp(metTime(2))
        metInputTimeStep = (secondTime(4) - firstTime(4)) * 86400.

        ! Find the closest cell to our lon and lat
        lonloc = closestCell(metFssId,'lon',longitude)
        latloc = closestCell(metFssId,'lat',latitude)

        ! Now read in the whole MET times series and store it for each variable
        allocate(metFss(validTimestep),metFdl(validTimestep),metPre(validTimestep),&
                 metTa(validTimestep),metQa(validTimestep),metUv(validTimestep),metPres(validTimestep))

        ! NOTE: Carefully check that your incoming inputs are in the expected units!

        ! Also take care here. If you use ncdump on a file it will show the opposite order for the
        ! dimensions of a variable than how fortran reads them in. So var(lat,lon,time) is actually
        ! var(time,lon,lat) from the perspective of fortran. Pay careful attention!
        
!         metFss = ncGet1DVar(metFssId, 'sw', start = [firstIndex,lonloc,latloc], count = [validTimestep,1,1])
!         metFdl = ncGet1DVar(metFdlId, 'lw', start = [firstIndex,lonloc,latloc], count = [validTimestep,1,1])
!         metPre = ncGet1DVar(metPreId, 'pr', start = [firstIndex,lonloc,latloc], count = [validTimestep,1,1])
!         metTa = ncGet1DVar(metTaId, 'ta', start = [firstIndex,lonloc,latloc], count = [validTimestep,1,1])
!         metQa = ncGet1DVar(metQaId, 'qa', start = [firstIndex,lonloc,latloc], count = [validTimestep,1,1])
!         metUv = ncGet1DVar(metUvId, 'wi', start = [firstIndex,lonloc,latloc], count = [validTimestep,1,1])
!         metPres = ncGet1DVar(metPresId, 'ap', start = [firstIndex,lonloc,latloc], count = [validTimestep,1,1])

        metFss = ncGet1DVar(metFssId, 'Incoming_Short_Wave_Radiation', start = [firstIndex,lonloc,latloc], count = [validTimestep,1,1])
        metFdl = ncGet1DVar(metFdlId, 'Incoming_Long_Wave_Radiation', start = [firstIndex,lonloc,latloc], count = [validTimestep,1,1])
        metPre = ncGet1DVar(metPreId, 'Total_Precipitation', start = [firstIndex,lonloc,latloc], count = [validTimestep,1,1])
        metTa = ncGet1DVar(metTaId, 'Temperature', start = [firstIndex,lonloc,latloc], count = [validTimestep,1,1])
        metQa = ncGet1DVar(metQaId, 'Air_Specific_Humidity', start = [firstIndex,lonloc,latloc], count = [validTimestep,1,1])
        metUv = ncGet1DVar(metUvId, 'U_wind_component', start = [firstIndex,lonloc,latloc], count = [validTimestep,1,1])
        metPres = ncGet1DVar(metPresId, 'Pression', start = [firstIndex,lonloc,latloc], count = [validTimestep,1,1])

    end subroutine getMet

    !>@}
    ! ------------------------------------------------------------------------------------

    !>\ingroup model_state_drivers_updateMet
    !!@{
    !> This transfers the met data of this time step from the read-in array to the
    !! instantaneous variables. This also sets iyear to the present year of MET being read in.

    subroutine updateMet(metTimeIndex,delt,iyear,iday,ihour,imin,metDone)

        use ctem_params, only : monthend
        use class_statevars, only : class_rot
        use generalUtils, only : parseTimeStamp

        implicit none

        integer, intent(inout) :: metTimeIndex      !< Index to read from met file
        real,    intent(in) :: delt                 !< Physics timestep (s)
        integer, intent(out) :: iyear               !< Present year of simulation
        integer, intent(out) :: iday                !< Present day of simulation
        integer, intent(out) :: ihour               !< Present hour of simulation
        integer, intent(out) :: imin                !< Present minute of simulation
        logical, intent(out) :: metDone             !< Switch signalling end of met data

        real, pointer, dimension(:) :: FDLROW       !< Downwelling longwave sky radiation \f$[W m^{-2} ]\f$
        real, pointer, dimension(:) :: FSSROW       !< Shortwave radiation \f$[W m^{-2} ]\f$
        real, pointer, dimension(:) :: PREROW       !< Surface precipitation rate \f$[kg m^{-2} s^{-1} ]\f$
        real, pointer, dimension(:) :: TAROW        !< Air temperature at reference height [K]
        real, pointer, dimension(:) :: QAROW        !< Specific humidity at reference height \f$[kg kg^{-1}]\f$
        real, pointer, dimension(:) :: UVROW        !< Wind speed at reference height \f$[m s^{-1} ]\f$
        real, pointer, dimension(:) :: PRESROW      !< Surface air pressure \f$[P_a]\f$

        integer :: i,numsteps
        real, dimension(5) :: theTime
        real :: dayfrac, month, dom, hour, minute

        FSSROW => class_rot%FSSROW
        FDLROW => class_rot%FDLROW
        PREROW => class_rot%PREROW
        TAROW => class_rot%TAROW
        QAROW => class_rot%QAROW
        UVROW => class_rot%UVROW
        PRESROW => class_rot%PRESROW

        metDone = .false.

        ! Find the timestep info from the array already read in.
        theTime =  parseTimeStamp(metTime(metTimeIndex))

        iyear = int(theTime(1))
        month = theTime(2)
        dom = theTime(3)
        dayfrac = theTime(4)
        iday = int(theTime(5))

        !> The dayfrac can then be parsed to give the hour and minute.
        numsteps = nint(dayfrac * 24. / (delt / 3600.))
        ihour = floor(real(numsteps) / 2.)
        minute = mod(numsteps,2)
        imin = nint(minute) * (delt / 60.)

        !> The meteorological data is then passed to the instantaneous variables
        !! from the larger variables that store the run's met data read in earlier.
        i = 1 ! always 1 offline
        FSSROW(I)   = metFss(metTimeIndex)
        FDLROW(i)   = metFdl(metTimeIndex)
        PREROW(i)   = metPre(metTimeIndex)
        TAROW(i)    = metTa(metTimeIndex) ! This is converted from the read-in degree C to K in main_driver!
        QAROW(i)    = metQa(metTimeIndex)
        UVROW(i)    = metUv(metTimeIndex)
        PRESROW(i)  = metPres(metTimeIndex)
        
        !> If the end of the timeseries is reached, change the metDone switch to true.
        if (metTimeIndex ==  size(metTime)) metDone = .true.

        return

    end subroutine updateMet

    !>@}
    ! ------------------------------------------------------------------------------------

    !>\ingroup model_state_drivers_closestCell
    !!@{
    !> Finds the closest grid cell in the file

    integer function closestCell(ncid,label,gridPoint)

        use fileIOModule

        implicit none

        integer, intent(in) :: ncid
        character(*), intent(in) :: label
        real, intent(in) :: gridPoint
        integer :: lengthdim
        real, dimension(:), allocatable :: filevals
        integer, dimension(1) :: tempintarr

        lengthdim = ncGetDimLen(ncid, label)
        allocate(filevals(lengthdim))
        filevals = ncGet1DVar(ncid, label, start = [1], count = [lengthdim])
        filevals = filevals - gridPoint
        tempintarr = minloc(abs(filevals))
        closestCell = tempintarr(1)

    end function closestCell
    !>@}
    ! ------------------------------------------------------------------------------------

    !>\ingroup model_state_drivers_deallocInput
    !!@{
    !> Deallocates the input files arrays

    subroutine deallocInput

        implicit none

        if (allocated(CO2Time))       deallocate(CO2Time)
        if (allocated(CO2FromFile))   deallocate(CO2FromFile)
        if (allocated(CH4Time))       deallocate(CH4Time)
        if (allocated(CH4FromFile))   deallocate(CH4FromFile)
        if (allocated(POPDTime))      deallocate(POPDTime)
        if (allocated(POPDFromFile))  deallocate(POPDFromFile)
        if (allocated(LGHTTime))      deallocate(LGHTTime)
        if (allocated(LGHTFromFile))  deallocate(LGHTFromFile)
        if (allocated(LUCTime))       deallocate(LUCTime)
        if (allocated(LUCFromFile))   deallocate(LUCFromFile)
        deallocate(metTime,metFss,metFdl,metPre,metPres,metQa,metTa,metUv)

    end subroutine deallocInput
!!@}
!>\file
!> Central driver to read in, and write out all model state variables (replacing INI and CTM files)
!! as well as the model inputs such as MET, population density, land use change, CO2 etc.

end module model_state_drivers
