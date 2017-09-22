module model_state_drivers

    !> This is the central driver to read in, and write out
    !! all model state variables (replacing INI and CTM files)
    !! as well as the model inputs such as MET, population density,
    !! land use change, CO2 etc.

    ! J. Melton
    ! Nov 2016

    use fileIOModule !, only : ncOpen, ncGetDimLen, ncGet1DVar, ncGet2DVar, ncGetDimValues,&
                     !        ncPut2DVar,ncPut3DVar,ncClose,ncGet3DVar

    implicit none

    public  :: read_modelsetup
    public  :: read_initialstate
    public  :: write_restart

contains

    !---

    subroutine read_modelsetup

        !> This reads in the model setup from the netcdf initialization file.
        !> The number of latitudes is always 1 offline while the maximum number of
        !> mosaics (nmos), the number of soil layers (ignd), are read from the netcdf.
        !> ilg is then calculated from nlat and nmos.

        ! J. Melton
        ! Feb 2017

        use io_driver, only : initid,rsid,myDomain
        use ctem_statevars,     only : c_switch
        use ctem_params, only : nmos,nlat,ignd,ilg  ! These are set in this subroutine!

        implicit none

        ! pointers:
        character(180), pointer         :: init_file
        character(180), pointer         :: rs_file_to_overwrite

        ! Local vars
        integer, allocatable, dimension(:,:) :: mask
        integer :: i, j
        integer :: totlon,totlat,totsize
        integer, dimension(1) :: pos
        integer, dimension(2) :: xpos,ypos
        integer, dimension(:,:), allocatable :: nmarray

        ! point pointers:
        init_file               => c_switch%init_file
        rs_file_to_overwrite    => c_switch%rs_file_to_overwrite

        ! ------------

        !> First, open initial conditions file.

        initid = ncOpen(init_file, NF90_NOWRITE)

        !> Next, retrieve dimensions. We assume the file has 'lon' and 'lat' for
        !! names of longitude and latitude.

        totlon = ncGetDimLen(initid,'lon')
        totlat = ncGetDimLen(initid,'lat')

        !calculate the number and indices of the pixels to be calculated
        !allocate(allLonValues(totlon), allLatValues(totlat))

        myDomain%allLonValues = ncGetDimValues(initid, 'lon', count = (/totlon/))
        myDomain%allLatValues = ncGetDimValues(initid, 'lat', count = (/totlat/))

        ! Based on the domainBounds, we make vectors of the cells to be run:
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
                 myDomain%lonLandIndex(totsize))

        !> Retrieve the number of soil layers (set ignd!)

        ignd = ncGetDimLen(initid, 'layer')

        !> Grab the model domain. We use GC since it is the land cells we want to run the model over.
        !! the 'Mask' variable is all land (we don't run over Antarctica).
        mask = ncGet2DVar(initid, 'GC', start = [myDomain%srtx, myDomain%srty],&
                          count = [myDomain%cntx, myDomain%cnty],format = [myDomain%cntx, myDomain%cnty])

        myDomain%LandCellCount = 0
        do i = 1, myDomain%cntx
            do j = 1, myDomain%cnty
                if (mask(i,j) .eq. -1) then
                    !print*, "(", i, ",", j, ") or (", allLonValues(i + srtx - 1), ",", allLatValues(j + srty - 1), ") is land"
                    myDomain%LandCellCount = myDomain%LandCellCount + 1
                    myDomain%lonLandCell(myDomain%LandCellCount) = myDomain%allLonValues(i + myDomain%srtx - 1)
                    myDomain%lonLandIndex(myDomain%LandCellCount) = i + myDomain%srtx - 1
                    myDomain%latLandCell(myDomain%LandCellCount) = myDomain%allLatValues(j + myDomain%srty - 1)
                    myDomain%latLandIndex(myDomain%LandCellCount) = j + myDomain%srty - 1
                endif
            enddo
        enddo

        nlat = 1

        !> To determine nmos, we use the largest number in the input file variable nmtest
        !! for the region we are running.

        nmarray = ncGet2DVar(initid, 'nmtest', start = [myDomain%srtx, myDomain%srty],&
                             count = [myDomain%cntx, myDomain%cnty],format = [myDomain%cntx, myDomain%cnty])
        nmos= maxval(nmarray)

        !> Determine the size of ilg which is nlat times nmos

        ilg = nlat * nmos

        !> Lastly, open the restart file so it is ready to be written to.
        rsid = ncOpen(rs_file_to_overwrite, nf90_write)

    end subroutine read_modelsetup

    !--------------------------------------------------------------------------------------------

    subroutine read_initialstate(lonIndex,latIndex)

        ! J. Melton
        ! Nov 2016

        use io_driver, only : initid
        use ctem_statevars,     only : c_switch,vrot,vgat
        use class_statevars,    only : class_rot,class_gat
        use ctem_params,        only : icc,iccp1,nmos,ignd,ilg,icp1,nlat,ican,abszero,pi

        implicit none

        ! arguments
        integer, intent(in) :: lonIndex,latIndex

        ! pointers:
        real, pointer, dimension(:,:,:) :: FCANROT
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

        logical, pointer :: ctem_on
        logical, pointer :: dofire
        logical, pointer :: compete
        logical, pointer :: inibioclim
        logical, pointer :: dowetlands
        logical, pointer :: start_bare
        logical, pointer :: lnduseon
        logical, pointer :: obswetf
        real, pointer, dimension(:,:,:) :: ailcminrow           !
        real, pointer, dimension(:,:,:) :: ailcmaxrow           !
        real, pointer, dimension(:,:,:) :: dvdfcanrow           !
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
        integer, pointer, dimension(:,:) :: ipeatlandrow   !<Peatland flag: 0 = not a peatland, 1= bog, 2 = fen
        real, pointer, dimension(:,:) :: Cmossmas          !<C in moss biomass, \f$kg C/m^2\f$
        real, pointer, dimension(:,:) :: litrmsmoss        !<moss litter mass, \f$kg C/m^2\f$
        real, pointer, dimension(:,:) :: dmoss             !<depth of living moss (m)
        real, pointer, dimension(:) :: grclarea            !<area of the grid cell, \f$km^2\f$

        ! local variables

        integer :: i,m,j
        real, dimension(ilg,2) :: crop_temp_frac
        real, allocatable, dimension(:,:) :: temptwod

        ! point pointers:
        ctem_on           => c_switch%ctem_on
        dofire            => c_switch%dofire
        compete           => c_switch%compete
        inibioclim        => c_switch%inibioclim
        dowetlands        => c_switch%dowetlands
        start_bare        => c_switch%start_bare
        lnduseon          => c_switch%lnduseon
        obswetf           => c_switch%obswetf
        ailcminrow        => vrot%ailcmin
        ailcmaxrow        => vrot%ailcmax
        dvdfcanrow        => vrot%dvdfcan
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

        ! ----------------------------

        do i = 1, nlat
         RADJROW(i)=DLATROW(i)*PI/180.
         Z0ORROW(i)=0.0
         GGEOROW(i)=0.0
        end do

        ZRFMROW = ncGet1DVar(initid, 'ZRFM', start = [lonIndex, latIndex], count = [1, 1])
        ZRFHROW = ncGet1DVar(initid, 'ZRFH', start = [lonIndex, latIndex], count = [1, 1])
        ZBLDROW = ncGet1DVar(initid, 'ZBLD', start = [lonIndex, latIndex], count = [1, 1])
        GCROW = ncGet1DVar(initid, 'GC', start = [lonIndex, latIndex], count = [1, 1])
        DRNROT = ncGet2DVar(initid, 'DRN', start = [lonIndex, latIndex, 1], count = [1, 1, nmos], format = [nlat, nmos])
        SDEPROT = ncGet2DVar(initid, 'SDEP', start = [lonIndex, latIndex, 1], count = [1, 1, nmos], format = [nlat, nmos])
        SOCIROT = ncGet2DVar(initid, 'SOCI', start = [lonIndex, latIndex, 1], count = [1, 1, nmos], format = [nlat, nmos])
        FAREROT = ncGet2DVar(initid, 'FARE', start = [lonIndex, latIndex, 1], count = [1, 1, nmos], format = [nlat, nmos])

        ! Error check:
        do i = 1,nlat
            do m = 1,nmos
                if (FAREROT(i,m) .gt. 1.0) then
                    print *,'FAREROT > 1',FAREROT(I,M)
                    call XIT('read_initialstate', -1)
                end if
            enddo
        enddo

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
        FCANROT = ncGet3DVar(initid, 'FCAN', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icp1, nmos], format = [nlat, nmos, icp1])
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

            ipeatlandrow = ncGet2DVar(initid, 'ipeatland', start = [lonIndex, latIndex, 1], count = [1, 1, nmos], format = [nlat, nmos])
            slopefrac = ncGet3DVar(initid, 'slopefrac', start = [lonIndex, latIndex, 1, 1], count = [1, 1, nmos, 8], format = [nlat, nmos, 8])
            Cmossmas = ncGet2DVar(initid, 'Cmossmas', start = [lonIndex, latIndex, 1], count = [1, 1, nmos], format = [nlat, nmos])
            litrmsmoss = ncGet2DVar(initid, 'litrmsmoss', start = [lonIndex, latIndex, 1], count = [1, 1, nmos], format = [nlat, nmos])
            dmoss = ncGet2DVar(initid, 'dmoss', start = [lonIndex, latIndex, 1], count = [1, 1, nmos], format = [nlat, nmos])
            ailcminrow = ncGet3DVar(initid, 'ailcmin', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos], format = [nlat, nmos,icc])
            ailcmaxrow = ncGet3DVar(initid, 'ailcmax', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos], format = [nlat, nmos,icc])
            dvdfcanrow = ncGet3DVar(initid, 'dvdfcan', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos], format = [nlat, nmos,icc])

            !>Rest of the initialization variables are needed to run CTEM but if starting from bare ground initialize all
            !>live and dead c pools from zero. suitable values of extnprobgrd and prbfrhucgrd would still be required. set
            !>stdaln to 1 for operation in non-gcm stand alone mode, in the CTEM initialization file.

            gleafmasrow = ncGet3DVar(initid, 'gleafmas', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos], format = [nlat, nmos,icc])
            bleafmasrow = ncGet3DVar(initid, 'bleafmas', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos], format = [nlat, nmos,icc])
            stemmassrow = ncGet3DVar(initid, 'stemmass', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos], format = [nlat, nmos,icc])
            rootmassrow = ncGet3DVar(initid, 'rootmass', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos], format = [nlat, nmos,icc])

            !>If fire and competition are on, save the stemmass and rootmass for use in burntobare subroutine on the first timestep.
            if (dofire .and. compete) then
                        do j =1,icc
                            pstemmassrow(i,m,j)=stemmassrow(i,m,j)
                            pgleafmassrow(i,m,j)=rootmassrow(i,m,j)
                        end do
            end if

            litrmassrow = ncGet4DVar(initid, 'litrmass', start = [lonIndex, latIndex, 1, 1, 1], count = [1, 1, icc, ignd, nmos], format = [nlat, nmos, icc, ignd])
            soilcmasrow = ncGet4DVar(initid, 'soilcmas', start = [lonIndex, latIndex, 1, 1,1], count = [1, 1, icc, ignd,nmos], format = [nlat, nmos,icc, ignd])
            lfstatusrow = ncGet3DVar(initid, 'lfstatus', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos], format = [nlat, nmos,icc])
            pandaysrow = ncGet3DVar(initid, 'pandays', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos], format = [nlat, nmos,icc])

            if (compete .and. inibioclim) then  !read in the bioclimatic parameters

                twarmm(:,1) = ncGet1DVar(initid, 'twarmm', start = [lonIndex, latIndex], count = [1, 1], format = [nlat])
                tcoldm(:,1) = ncGet1DVar(initid, 'tcoldm', start = [lonIndex, latIndex], count = [1, 1], format = [nlat])
                gdd5(:,1) = ncGet1DVar(initid, 'gdd5', start = [lonIndex, latIndex], count = [1, 1], format = [nlat])
                aridity(:,1) = ncGet1DVar(initid, 'aridity', start = [lonIndex, latIndex], count = [1, 1], format = [nlat])
                srplsmon(:,1) = ncGet1DVar(initid, 'srplsmon', start = [lonIndex, latIndex], count = [1, 1], format = [nlat])
                defctmon(:,1) = ncGet1DVar(initid, 'defctmon', start = [lonIndex, latIndex], count = [1, 1], format = [nlat])
                anndefct(:,1) = ncGet1DVar(initid, 'anndefct', start = [lonIndex, latIndex], count = [1, 1], format = [nlat])
                annsrpls(:,1) = ncGet1DVar(initid, 'annsrpls', start = [lonIndex, latIndex], count = [1, 1], format = [nlat])
                annpcp(:,1) = ncGet1DVar(initid, 'annpcp', start = [lonIndex, latIndex], count = [1, 1], format = [nlat])
                dry_season_length(:,1) = ncGet1DVar(initid, 'dry_season_length', start = [lonIndex, latIndex], count = [1, 1], format = [nlat])

                deallocate(temptwod)

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

            else if (compete .and. .not. inibioclim) then ! set them to zero

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

            !>if this run uses the competition or lnduseon parameterization and starts from bare ground, set up the model state here. this
            !>overwrites what was read in from the initialization file. for composite runs (the composite set up is after this one for mosaics)
            if ((compete .or. lnduseon) .and. start_bare) then

                !>set up for composite runs when start_bare is on and compete or landuseon

                !>store the read-in crop fractions as we keep them even when we start bare.
                !!FLAG: this is setup assuming that crops are in pft number 6 and 7.
                !!and the first tile contains the information for the grid cell (assumes we have crops in
                !!every tile too! JM Apr 9 2014.
                do i=1,nlat
                    crop_temp_frac(i,1)=FCANROT(i,1,3)*dvdfcanrow(i,1,6)
                    crop_temp_frac(i,2)=FCANROT(i,1,3)*dvdfcanrow(i,1,7)
                end do

                !>initalize to zero, these will be filled in by the luc or competition subroutines.
                FCANROT=0.0
                dvdfcanrow=0.0

                ! FLAG- needed anymore?
                ! Added this as start_bare runs were not properly assigning
                ! a TCAN on the very first day since the FCANROT was 0. JM Jan 14 2014.
!                 do i=1,nltest
!                     do m = 1,nmtest
!                         do j=1,icp1
!                         if (j .lt. icp1) then
!                             FCANROT(i,m,j)=seed
!                         else
!                             FCANROT(i,m,j)=1.0 - (real(ican) * seed)
!                         endif
!                         end do
!                     end do
!                 end do

                do i=1,nlat
                    do m = 1,nmos

                        !>initial conditions always required
                        dvdfcanrow(i,m,1)=1.0  !ndl
                        dvdfcanrow(i,m,3)=1.0  !bdl
                        dvdfcanrow(i,m,6)=1.0  !crop
                        dvdfcanrow(i,m,8)=1.0  !grasses

                        do j = 1,icc
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

                do i=1,nlat
                    do m = 1,nmos
                        FCANROT(i,m,3) = crop_temp_frac(i,1) + crop_temp_frac(i,2)
                        if (FCANROT(i,m,3) .gt. abszero) then
                            dvdfcanrow(i,m,6) = crop_temp_frac(i,1) / FCANROT(i,m,3)
                            dvdfcanrow(i,m,7) = crop_temp_frac(i,2) / FCANROT(i,m,3)
                        else
                            dvdfcanrow(i,m,6) = 1.0
                            dvdfcanrow(i,m,7) = 0.0
                        end if
                    end do !nmtest
                end do !nltest

            end if !if (compete/landuseon .and. start_bare)

        end if !ctem_on

    end subroutine read_initialstate

    !---------------------------------------------------------------------------------------------

    subroutine write_restart(lonIndex,latIndex)

        !! Write out the model restart file to netcdf. We only write out the variables that the model
        !! influences

        ! J. Melton
        ! Jun 2017

        !use netcdf
        !use netcdf_drivers, only : check_nc
        !use serialFileIOModule
        use io_driver, only : rsid
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
        logical, pointer :: compete
        logical, pointer :: lnduseon
        real, pointer, dimension(:,:,:) :: ailcminrow           !
        real, pointer, dimension(:,:,:) :: ailcmaxrow           !
        real, pointer, dimension(:,:,:) :: dvdfcanrow           !
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
        integer :: i,m,j,k1c,k2c,n
        integer, dimension(nlat,nmos) :: icountrow
        real, dimension(icc) :: rnded_pft

        ! point pointers:
        ctem_on           => c_switch%ctem_on
        compete           => c_switch%compete
        lnduseon          => c_switch%lnduseon
        ailcminrow        => vrot%ailcmin
        ailcmaxrow        => vrot%ailcmax
        dvdfcanrow        => vrot%dvdfcan
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


        call ncPut2DVar(rsid, 'FARE', FAREROT, start = [lonIndex, latIndex, 1], count = [1, 1, nmos])
        call ncPut3DVar(rsid, 'FCAN', FCANROT, start = [lonIndex, latIndex, 1, 1], count = [1, 1, icp1, nmos])
        call ncPut3DVar(rsid, 'THLQ', THLQROT, start = [lonIndex, latIndex, 1, 1], count = [1, 1, ignd, nmos])
        call ncPut3DVar(rsid, 'THIC', THICROT, start = [lonIndex, latIndex, 1, 1], count = [1, 1, ignd, nmos])
        call ncPut3DVar(rsid, 'TBAR', TBARROT-273.16, start = [lonIndex, latIndex, 1, 1], count = [1, 1, ignd, nmos])
        call ncPut2DVar(rsid, 'TCAN', TCANROT-273.16, start = [lonIndex, latIndex, 1], count = [1, 1, nmos])
        call ncPut2DVar(rsid, 'TSNO', TSNOROT-273.16, start = [lonIndex, latIndex, 1], count = [1, 1, nmos])
        call ncPut2DVar(rsid, 'TPND', TPNDROT-273.16, start = [lonIndex, latIndex, 1], count = [1, 1, nmos])
        call ncPut2DVar(rsid, 'ZPND', ZPNDROT, start = [lonIndex, latIndex, 1], count = [1, 1, nmos])
        call ncPut2DVar(rsid, 'RCAN', RCANROT, start = [lonIndex, latIndex, 1], count = [1, 1, nmos])
        call ncPut2DVar(rsid, 'SCAN', SCANROT, start = [lonIndex, latIndex, 1], count = [1, 1, nmos])
        call ncPut2DVar(rsid, 'SNO', SNOROT, start = [lonIndex, latIndex, 1], count = [1, 1, nmos])
        call ncPut2DVar(rsid, 'ALBS', ALBSROT, start = [lonIndex, latIndex, 1], count = [1, 1, nmos])
        call ncPut2DVar(rsid, 'RHOS', RHOSROT, start = [lonIndex, latIndex, 1], count = [1, 1, nmos])
        call ncPut2DVar(rsid, 'GRO', GROROT, start = [lonIndex, latIndex, 1], count = [1, 1, nmos])

        if (ctem_on) then

            !> if landuseon or competition, then we need to recreate the dvdfcanrow so do so now
            if (lnduseon .or. compete ) then
                icountrow=0
                do j = 1, ican
                    do i = 1,nlat
                        do m = 1,nmos
                            k1c = (j-1)*l2max + 1
                            k2c = k1c + (l2max - 1)
                            do n = k1c, k2c
                                if (modelpft(n) .eq. 1) then
                                    icountrow(i,m) = icountrow(i,m) + 1
                                    if (FCANROT(i,m,j) .gt. 0.) then
                                        dvdfcanrow(i,m,icountrow(i,m)) = fcancmxrow(i,m,icountrow(i,m))/FCANROT(i,m,j)
                                    else
                                        dvdfcanrow(i,m,icountrow(i,m)) = 0.
                                    end if
                                end if !modelpft
                            end do !n
                            !> check to ensure that the dvdfcanrow's add up to 1 across a class-level pft
                            if (dvdfcanrow(i,m,1) .eq. 0. .and. dvdfcanrow(i,m,2) .eq. 0.) then
                                dvdfcanrow(i,m,1)=1.0
                            else if (dvdfcanrow(i,m,3) .eq. 0. .and. dvdfcanrow(i,m,4) .eq. 0. .and. dvdfcanrow(i,m,5) .eq. 0.) then
                                dvdfcanrow(i,m,3)=1.0
                            else if (dvdfcanrow(i,m,6) .eq. 0. .and. dvdfcanrow(i,m,7) .eq. 0.) then
                                dvdfcanrow(i,m,6)=1.0
                            else if (dvdfcanrow(i,m,8) .eq. 0. .and. dvdfcanrow(i,m,9) .eq. 0.) then
                                dvdfcanrow(i,m,8)=1.0
                            end if
                        end do !m
                    enddo !i
                enddo !j

                do i=1,nlat
                    do m=1,nmos
                        do j = 1, icc
                            !>Lastly check if the different pfts accidently add up > 1.0 after rounding to the number of sig figs used in the output
                            !>this rounds to 3 decimal places. if you are found to be over or under, arbitrarily reduce one of the pfts. the amount of
                            !>the change will be inconsequential.
                            rnded_pft(j) =real(int(dvdfcanrow(i,m,j) * 1000.0))/ 1000.0
                            dvdfcanrow(i,m,j) = rnded_pft(j)
                        end do

                        if (dvdfcanrow(i,m,1) + dvdfcanrow(i,m,2) .ne. 1.0) then
                            dvdfcanrow(i,m,1) = 1.0 - rnded_pft(2)
                            dvdfcanrow(i,m,2) = rnded_pft(2)
                        end if
                        if (dvdfcanrow(i,m,3) + dvdfcanrow(i,m,4) +  dvdfcanrow(i,m,5) .ne. 1.0) then
                            dvdfcanrow(i,m,3) = 1.0 - rnded_pft(4) - rnded_pft(5)
                            dvdfcanrow(i,m,4) = rnded_pft(4)
                            dvdfcanrow(i,m,5) = rnded_pft(5)
                        end if
                        if (dvdfcanrow(i,m,6) + dvdfcanrow(i,m,7) .ne. 1.0) then
                            dvdfcanrow(i,m,6) = 1.0 - rnded_pft(7)
                            dvdfcanrow(i,m,7) = rnded_pft(7)
                        end if
                        if (dvdfcanrow(i,m,8) + dvdfcanrow(i,m,9) .ne. 1.0) then
                            dvdfcanrow(i,m,8) = 1.0 - rnded_pft(9)
                            dvdfcanrow(i,m,9) = rnded_pft(9)
                        end if
                    enddo
                enddo
            end if !lnuse/compete

            call ncPut3DVar(rsid, 'ailcmin', ailcminrow, start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos])
            call ncPut3DVar(rsid, 'ailcmax', ailcmaxrow, start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos])
            call ncPut3DVar(rsid, 'dvdfcan', dvdfcanrow, start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos])
            call ncPut3DVar(rsid, 'gleafmas', gleafmasrow, start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos])
            call ncPut3DVar(rsid, 'bleafmas', bleafmasrow, start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos])
            call ncPut3DVar(rsid, 'stemmass', stemmassrow, start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos])
            call ncPut3DVar(rsid, 'rootmass', rootmassrow, start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos])
            call ncPut3DVar(rsid, 'litrmass', litrmassrow, start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos])
            call ncPut3DVar(rsid, 'soilcmas', soilcmasrow, start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos])
            call ncPut3DVar(rsid, 'lfstatus', real(lfstatusrow), start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos])
            call ncPut3DVar(rsid, 'pandays', real(pandaysrow), start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos])
            call ncPut2DVar(rsid, 'Cmossmas', Cmossmas, start = [lonIndex, latIndex, 1], count = [1, 1, nmos])
            call ncPut2DVar(rsid, 'litrmsmoss', litrmsmoss, start = [lonIndex, latIndex, 1], count = [1, 1, nmos])
            call ncPut2DVar(rsid, 'dmoss', dmoss, start = [lonIndex, latIndex, 1], count = [1, 1, nmos])

            if (compete) then

                call ncPut2DVar(rsid, 'twarmm', twarmm, start = [lonIndex, latIndex])
                call ncPut2DVar(rsid, 'tcoldm', tcoldm, start = [lonIndex, latIndex])
                call ncPut2DVar(rsid, 'gdd5', gdd5, start = [lonIndex, latIndex])
                call ncPut2DVar(rsid, 'aridity', aridity, start = [lonIndex, latIndex])
                call ncPut2DVar(rsid, 'srplsmon', srplsmon, start = [lonIndex, latIndex])
                call ncPut2DVar(rsid, 'defctmon', defctmon, start = [lonIndex, latIndex])
                call ncPut2DVar(rsid, 'anndefct', anndefct, start = [lonIndex, latIndex])
                call ncPut2DVar(rsid, 'annsrpls', annsrpls, start = [lonIndex, latIndex])
                call ncPut2DVar(rsid, 'annpcp', annpcp, start = [lonIndex, latIndex])
                call ncPut2DVar(rsid, 'dry_season_length', dry_season_length, start = [lonIndex, latIndex])

            end if ! compete

        end if !ctem_on

    end subroutine write_restart

end module model_state_drivers
