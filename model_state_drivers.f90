module model_state_drivers

    !> This is the central driver to read in, and write out
    !! all model state variables (replacing INI and CTM files)
    !! as well as the model inputs such as MET, population density,
    !! land use change, CO2 etc.

    ! J. Melton
    ! Nov 2016

    implicit none

    private
    public :: read_modelsetup
    public :: read_initialstate

contains

    !---

    subroutine read_modelsetup()

        !> This reads in the model setup from the netcdf initialization file.
        !> The number of latitudes is always 1 offline while the maximum number of
        !> mosaics (nmos), the number of soil layers (ignd), are read from the netcdf.
        !> ilg is then calculated from nlat and nmos.

        ! J. Melton
        ! Feb 2017

        use netcdf
        use netcdf_drivers, only : check_nc
        use io_driver, only : initid,cntx,cnty,srtx,srty,bounds,lonvect,latvect
        use ctem_statevars,     only : c_switch
        use ctem_params, only : nmos,nlat,ignd,ilg  ! These are set in this subroutine!

        implicit none

        ! pointers:
        character(180), pointer         :: init_file

        integer :: dimid
        integer :: varid
        integer :: xsize, ysize
        integer, dimension(1) :: pos
        integer, dimension(2) :: xpos,ypos
        integer, dimension(:,:), allocatable :: nmarray
        real, dimension(:), allocatable :: all_lon,all_lat

        ! point pointers:
        init_file         => c_switch%init_file

        ! ------------

        !> First, open initial conditions file.

        call check_nc(nf90_open(trim(init_file),nf90_nowrite,initid))

        !> Next, retrieve dimensions. We assume the file has 'lon' and 'lat' for
        !! names of longitude and latitude.

        call check_nc(nf90_inq_dimid(initid,'lon',dimid))
        call check_nc(nf90_inquire_dimension(initid,dimid,len=xsize))

        call check_nc(nf90_inq_dimid(initid,'lat',dimid))
        call check_nc(nf90_inquire_dimension(initid,dimid,len=ysize))

        !calculate the number and indices of the pixels to be calculated

        allocate(all_lon(xsize),&
            all_lat(ysize))

        call check_nc(nf90_inq_varid(initid,'lon',varid))
        call check_nc(nf90_get_var(initid,varid,all_lon))

        call check_nc(nf90_inq_varid(initid,'lat',varid))
        call check_nc(nf90_get_var(initid,varid,all_lat))

        ! Based on the bounds, we make vectors of the cells to be run:
        pos = minloc(abs(all_lon - bounds(1)))
        xpos(1) = pos(1)

        pos = minloc(abs(all_lon - bounds(2)))
        xpos(2) = pos(1)

        pos = minloc(abs(all_lat - bounds(3)))
        ypos(1) = pos(1)

        pos = minloc(abs(all_lat - bounds(4)))
        ypos(2) = pos(1)

        srtx = minval(xpos)
        srty = minval(ypos)

        if (all_lon(srtx) < bounds(1) .and. bounds(2) /= bounds(1)) srtx = srtx + 1
        cntx = 1 + abs(maxval(xpos) - srtx)

        if (all_lat(srty) < bounds(3) .and. bounds(4) /= bounds(3)) srty = srty + 1
        cnty = 1 + abs(maxval(ypos) - srty)

        !> The size of nlat should then be cntx x cnty. We later take in GC to determine
        !! which cells are valid land cells and use that in GATPREP to make it so we only
        !! do computations over the valid land cells.

        nlat = cntx * cnty

        !> Save the longitudes and latitudes over the region of interest for making the
        !! output files.
        allocate(lonvect(cntx),&
            latvect(cnty),&
            nmarray(cnty,cntx))

        lonvect = all_lon(srtx:srtx+cntx-1)
        latvect = all_lat(srty:srty+cnty-1)

        !> Retrieve the number of soil layers (set ignd!)

        call check_nc(nf90_inq_dimid(initid,'layer',dimid))
        call check_nc(nf90_inquire_dimension(initid,dimid,len=ignd))

        !> To determine nmos, we use the largest number in the input file variable nmtest
        !! for the region we are running.

        call check_nc(nf90_inq_varid(initid,'nmtest',varid))
        call check_nc(nf90_get_var(initid,varid,nmarray,start=[srtx,srty],count=[cntx,cnty]))

        nmos= maxval(nmarray)

        deallocate(nmarray,&
            all_lat,&
            all_lon)

        !> Lastly we determine the size of ilg which is nlat times nmos

        ilg = nlat * nmos

    end subroutine read_modelsetup


    !--------------------------------------------------------------------------------------------

    subroutine read_initialstate()

        ! J. Melton
        ! Nov 2016

        use netcdf
        use netcdf_drivers, only : check_nc
        use io_driver, only : initid,cntx,cnty,srtx,srty
        use ctem_statevars,     only : c_switch,vrot
        use class_statevars,    only : class_rot
        use ctem_params,        only : icc,iccp1,nmos,seed,ignd,ilg,icp1,nlat,ican,abszero

        implicit none

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
!        real, pointer, dimension(:,:,:) :: DELZ
!        real, pointer, dimension(:,:,:) :: ZBOT
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

        ! local variables

        integer :: i,m,j,x,y,varid,k
        real, dimension(ilg,2) :: crop_temp_frac
        real, allocatable, dimension(:,:) :: temptwod
        real, allocatable, dimension(:,:,:) :: temp3d
        integer, allocatable, dimension(:,:,:) :: temp3di
        real, allocatable, dimension(:,:,:,:) :: temp4d

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
!        DELZ              => class_rot%DELZ
!        ZBOT              => class_rot%ZBOT
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

        ! ----------------------------

        JLAT=NINT(DLATROW(1))
        RADJROW(1)=DLATROW(1)*PI/180.
        DLONROW(1)=DEGLON
        Z0ORROW(1)=0.0
        GGEOROW(1)=0.0
        !     GGEOROW(1)=-0.035  !commented out in original code.

        allocate(temptwod(cntx,cnty))

        call check_nc(nf90_inq_varid(initid,'ZRFM',varid))
        call check_nc(nf90_get_var(initid,varid,temptwod,start=[srtx,srty],count=[cntx,cnty]))
        ZRFMROW = reshape(temptwod,[ nlat ] )

        call check_nc(nf90_inq_varid(initid,'ZRFH',varid))
        call check_nc(nf90_get_var(initid,varid,temptwod,start=[srtx,srty],count=[cntx,cnty]))
        ZRFHROW = reshape(temptwod,[ nlat ] )

        call check_nc(nf90_inq_varid(initid,'ZBLD',varid))
        call check_nc(nf90_get_var(initid,varid,temptwod,start=[srtx,srty],count=[cntx,cnty]))
        ZBLDROW = reshape(temptwod,[ nlat ] )

        call check_nc(nf90_inq_varid(initid,'GC',varid))
        call check_nc(nf90_get_var(initid,varid,temptwod,start=[srtx,srty],count=[cntx,cnty]))
        GCROW = reshape(temptwod,[ nlat ] )

        deallocate(temptwod)

        ! Now get the 3D variables:
        allocate(temp3d(cntx,cnty,nmos))

        call check_nc(nf90_inq_varid(initid,'DRN',varid))
        call check_nc(nf90_get_var(initid,varid,temp3d,start=[srtx,srty,1],count=[cntx,cnty,nmos]))
        DRNROT = reshape(temp3d,[ nlat,nmos ] )

        call check_nc(nf90_inq_varid(initid,'SDEP',varid))
        call check_nc(nf90_get_var(initid,varid,temp3d,start=[srtx,srty,1],count=[cntx,cnty,nmos]))
        SDEPROT = reshape(temp3d,[ nlat,nmos ] )
        !
        !call check_nc(nf90_inq_varid(initid,'FARE',varid))
        !call check_nc(nf90_get_var(initid,varid,temp3d,start=[srtx,srty,1],count=[cntx,cnty,nmos]))
        !FAREROT = reshape(temp3d,[ nlat,nmos ] )
        !
                !  write(*,*)FAREROT
        ! !
        ! !           ! Error check:
        ! !           if (FAREROT(I,M) .gt. 1.0) then
        ! !            write(*,*)'FAREROT > 1',FAREROT(I,M)
        ! !            call XIT('runclass36ctem', -2)
        ! !           end if
        !
        !
        call check_nc(nf90_inq_varid(initid,'XSLP',varid))
        call check_nc(nf90_get_var(initid,varid,temp3d,start=[srtx,srty,1],count=[cntx,cnty,nmos]))
        XSLPROT = reshape(temp3d,[ nlat,nmos ] )

        call check_nc(nf90_inq_varid(initid,'GRKF',varid))
        call check_nc(nf90_get_var(initid,varid,temp3d,start=[srtx,srty,1],count=[cntx,cnty,nmos]))
        GRKFROT = reshape(temp3d,[ nlat,nmos ] )

        call check_nc(nf90_inq_varid(initid,'WFSF',varid))
        call check_nc(nf90_get_var(initid,varid,temp3d,start=[srtx,srty,1],count=[cntx,cnty,nmos]))
        WFSFROT = reshape(temp3d,[ nlat,nmos ] )

        call check_nc(nf90_inq_varid(initid,'WFCI',varid))
        call check_nc(nf90_get_var(initid,varid,temp3d,start=[srtx,srty,1],count=[cntx,cnty,nmos]))
        WFCIROT = reshape(temp3d,[ nlat,nmos ] )

        call check_nc(nf90_inq_varid(initid,'TCAN',varid))
        call check_nc(nf90_get_var(initid,varid,temp3d,start=[srtx,srty,1],count=[cntx,cnty,nmos]))
        TCANROT = reshape(temp3d,[ nlat,nmos ] )

        call check_nc(nf90_inq_varid(initid,'TSNO',varid))
        call check_nc(nf90_get_var(initid,varid,temp3d,start=[srtx,srty,1],count=[cntx,cnty,nmos]))
        TSNOROT = reshape(temp3d,[ nlat,nmos ] )

        call check_nc(nf90_inq_varid(initid,'TPND',varid))
        call check_nc(nf90_get_var(initid,varid,temp3d,start=[srtx,srty,1],count=[cntx,cnty,nmos]))
        TPNDROT = reshape(temp3d,[ nlat,nmos ] )

        call check_nc(nf90_inq_varid(initid,'ZPND',varid))
        call check_nc(nf90_get_var(initid,varid,temp3d,start=[srtx,srty,1],count=[cntx,cnty,nmos]))
        ZPNDROT = reshape(temp3d,[ nlat,nmos ] )

        call check_nc(nf90_inq_varid(initid,'RCAN',varid))
        call check_nc(nf90_get_var(initid,varid,temp3d,start=[srtx,srty,1],count=[cntx,cnty,nmos]))
        RCANROT = reshape(temp3d,[ nlat,nmos ] )

        call check_nc(nf90_inq_varid(initid,'SCAN',varid))
        call check_nc(nf90_get_var(initid,varid,temp3d,start=[srtx,srty,1],count=[cntx,cnty,nmos]))
        SCANROT = reshape(temp3d,[ nlat,nmos ] )

        call check_nc(nf90_inq_varid(initid,'SNO',varid))
        call check_nc(nf90_get_var(initid,varid,temp3d,start=[srtx,srty,1],count=[cntx,cnty,nmos]))
        SNOROT = reshape(temp3d,[ nlat,nmos ] )

        call check_nc(nf90_inq_varid(initid,'ALBS',varid))
        call check_nc(nf90_get_var(initid,varid,temp3d,start=[srtx,srty,1],count=[cntx,cnty,nmos]))
        ALBSROT = reshape(temp3d,[ nlat,nmos ] )

        call check_nc(nf90_inq_varid(initid,'RHOS',varid))
        call check_nc(nf90_get_var(initid,varid,temp3d,start=[srtx,srty,1],count=[cntx,cnty,nmos]))
        RHOSROT = reshape(temp3d,[ nlat,nmos ] )

        call check_nc(nf90_inq_varid(initid,'GRO',varid))
        call check_nc(nf90_get_var(initid,varid,temp3d,start=[srtx,srty,1],count=[cntx,cnty,nmos]))
        GROROT = reshape(temp3d,[ nlat,nmos ] )

        allocate(temp3di(cntx,cnty,nmos)) !this is an integer

        call check_nc(nf90_inq_varid(initid,'MID',varid))
        call check_nc(nf90_get_var(initid,varid,temp3di,start=[srtx,srty,1],count=[cntx,cnty,nmos]))
        MIDROT = reshape(temp3di,[ nlat,nmos ] )

        deallocate(temp3d)
        deallocate(temp3di)

!         allocate(temp3d(cntx,cnty,ignd))
!
!         call check_nc(nf90_inq_varid(initid,'ZBOT',varid)) !per layer but not per tile (fix?)
!         call check_nc(nf90_get_var(initid,varid,temp4d,start=[srtx,srty,1,1],count=[cntx,cnty,ignd]))
!         ZBOT = reshape(temp4d,[nlat,ignd])
!
!         call check_nc(nf90_inq_varid(initid,'DELZ',varid)) !per layer but not per tile (fix?)
!         call check_nc(nf90_get_var(initid,varid,temp4d,start=[srtx,srty,1,1],count=[cntx,cnty,ignd]))
!         DELZ = reshape(temp4d,[nlat,ignd])
!
!         deallocate(temp3d)

        ! Now get the 4D variables:

        allocate(temp4d(cntx,cnty,icp1,nmos))

        call check_nc(nf90_inq_varid(initid,'FCAN',varid))
        call check_nc(nf90_get_var(initid,varid,temp4d,start=[srtx,srty,1,1],count=[cntx,cnty,icp1,nmos]))
        FCANROT = reshape(temp4d,[nlat,nmos,icp1])

        call check_nc(nf90_inq_varid(initid,'LNZ0',varid))
        call check_nc(nf90_get_var(initid,varid,temp4d,start=[srtx,srty,1,1],count=[cntx,cnty,icp1,nmos]))
        LNZ0ROT = reshape(temp4d,[nlat,nmos,icp1])

        call check_nc(nf90_inq_varid(initid,'ALVC',varid))
        call check_nc(nf90_get_var(initid,varid,temp4d,start=[srtx,srty,1,1],count=[cntx,cnty,icp1,nmos]))
        ALVCROT = reshape(temp4d,[nlat,nmos,icp1])

        call check_nc(nf90_inq_varid(initid,'ALIC',varid))
        call check_nc(nf90_get_var(initid,varid,temp4d,start=[srtx,srty,1,1],count=[cntx,cnty,icp1,nmos]))
        ALICROT = reshape(temp4d,[nlat,nmos,icp1])

        deallocate(temp4d)

        allocate(temp4d(cntx,cnty,ican,nmos))

        call check_nc(nf90_inq_varid(initid,'PAMN',varid))
        call check_nc(nf90_get_var(initid,varid,temp4d,start=[srtx,srty,1,1],count=[cntx,cnty,ican,nmos]))
        PAMNROT = reshape(temp4d,[nlat,nmos,ican])

        call check_nc(nf90_inq_varid(initid,'PAMX',varid))
        call check_nc(nf90_get_var(initid,varid,temp4d,start=[srtx,srty,1,1],count=[cntx,cnty,ican,nmos]))
        PAMXROT = reshape(temp4d,[nlat,nmos,ican])

        call check_nc(nf90_inq_varid(initid,'PAMX',varid))
        call check_nc(nf90_get_var(initid,varid,temp4d,start=[srtx,srty,1,1],count=[cntx,cnty,ican,nmos]))
        PAMXROT = reshape(temp4d,[nlat,nmos,ican])

        call check_nc(nf90_inq_varid(initid,'CMAS',varid))
        call check_nc(nf90_get_var(initid,varid,temp4d,start=[srtx,srty,1,1],count=[cntx,cnty,ican,nmos]))
        CMASROT = reshape(temp4d,[nlat,nmos,ican])

        call check_nc(nf90_inq_varid(initid,'ROOT',varid))
        call check_nc(nf90_get_var(initid,varid,temp4d,start=[srtx,srty,1,1],count=[cntx,cnty,ican,nmos]))
        ROOTROT = reshape(temp4d,[nlat,nmos,ican])

        call check_nc(nf90_inq_varid(initid,'RSMN',varid))
        call check_nc(nf90_get_var(initid,varid,temp4d,start=[srtx,srty,1,1],count=[cntx,cnty,ican,nmos]))
        RSMNROT = reshape(temp4d,[nlat,nmos,ican])

        call check_nc(nf90_inq_varid(initid,'QA50',varid))
        call check_nc(nf90_get_var(initid,varid,temp4d,start=[srtx,srty,1,1],count=[cntx,cnty,ican,nmos]))
        QA50ROT = reshape(temp4d,[nlat,nmos,ican])

        call check_nc(nf90_inq_varid(initid,'VPDA',varid))
        call check_nc(nf90_get_var(initid,varid,temp4d,start=[srtx,srty,1,1],count=[cntx,cnty,ican,nmos]))
        VPDAROT = reshape(temp4d,[nlat,nmos,ican])

        call check_nc(nf90_inq_varid(initid,'VPDB',varid))
        call check_nc(nf90_get_var(initid,varid,temp4d,start=[srtx,srty,1,1],count=[cntx,cnty,ican,nmos]))
        VPDBROT = reshape(temp4d,[nlat,nmos,ican])

        call check_nc(nf90_inq_varid(initid,'PSGA',varid))
        call check_nc(nf90_get_var(initid,varid,temp4d,start=[srtx,srty,1,1],count=[cntx,cnty,ican,nmos]))
        PSGAROT = reshape(temp4d,[nlat,nmos,ican])

        call check_nc(nf90_inq_varid(initid,'PSGB',varid))
        call check_nc(nf90_get_var(initid,varid,temp4d,start=[srtx,srty,1,1],count=[cntx,cnty,ican,nmos]))
        PSGBROT = reshape(temp4d,[nlat,nmos,ican])

        deallocate(temp4d)

        allocate(temp4d(cntx,cnty,ignd,nmos))

        call check_nc(nf90_inq_varid(initid,'SAND',varid))
        call check_nc(nf90_get_var(initid,varid,temp4d,start=[srtx,srty,1,1],count=[cntx,cnty,ignd,nmos]))
        SANDROT = reshape(temp4d,[nlat,nmos,ignd])

        call check_nc(nf90_inq_varid(initid,'CLAY',varid))
        call check_nc(nf90_get_var(initid,varid,temp4d,start=[srtx,srty,1,1],count=[cntx,cnty,ignd,nmos]))
        CLAYROT = reshape(temp4d,[nlat,nmos,ignd])

        call check_nc(nf90_inq_varid(initid,'ORGM',varid))
        call check_nc(nf90_get_var(initid,varid,temp4d,start=[srtx,srty,1,1],count=[cntx,cnty,ignd,nmos]))
        ORGMROT = reshape(temp4d,[nlat,nmos,ignd])

        call check_nc(nf90_inq_varid(initid,'TBAR',varid))
        call check_nc(nf90_get_var(initid,varid,temp4d,start=[srtx,srty,1,1],count=[cntx,cnty,ignd,nmos]))
        TBARROT = reshape(temp4d,[nlat,nmos,ignd])

        call check_nc(nf90_inq_varid(initid,'THLQ',varid))
        call check_nc(nf90_get_var(initid,varid,temp4d,start=[srtx,srty,1,1],count=[cntx,cnty,ignd,nmos]))
        THLQROT = reshape(temp4d,[nlat,nmos,ignd])

        call check_nc(nf90_inq_varid(initid,'THIC',varid))
        call check_nc(nf90_get_var(initid,varid,temp4d,start=[srtx,srty,1,1],count=[cntx,cnty,ignd,nmos]))
        THICROT = reshape(temp4d,[nlat,nmos,ignd])

        deallocate(temp4d)

        if (ctem_on) then

            allocate(temptwod(cntx,cnty))

            call check_nc(nf90_inq_varid(initid,'extnprob',varid))
            call check_nc(nf90_get_var(initid,varid,temptwod,start=[srtx,srty],count=[cntx,cnty]))
            extnprob = reshape(temptwod,[ nlat ] )

            call check_nc(nf90_inq_varid(initid,'prbfrhuc',varid))
            call check_nc(nf90_get_var(initid,varid,temptwod,start=[srtx,srty],count=[cntx,cnty]))
            prbfrhuc = reshape(temptwod,[ nlat ] )

            deallocate(temptwod)

            allocate(temp3d(cntx,cnty,12))

            call check_nc(nf90_inq_varid(initid,'mlightng',varid))
            call check_nc(nf90_get_var(initid,varid,temp3d,start=[srtx,srty,1],count=[cntx,cnty,12]))
            mlightng = reshape(temp3d,[ nlat, 12 ] )

            deallocate(temp3d)
            allocate(temp4d(cntx,cnty,icc,nmos))

            !>The following three variables are needed to run CTEM. 1) min & 2) max leaf area index are needed to break
            !>class lai into dcd and evg for trees (for crops and grasses it doesn't matter much). 3) dvdfcanrow is
            !>needed to divide needle & broad leaf into dcd and evg, and crops & grasses into c3 and c4 fractions.

            call check_nc(nf90_inq_varid(initid,'ailcmin',varid))
            call check_nc(nf90_get_var(initid,varid,temp4d,start=[srtx,srty,1,1],count=[cntx,cnty,icc,nmos]))
            ailcminrow = reshape(temp4d,[nlat,nmos,icc])

            call check_nc(nf90_inq_varid(initid,'ailcmax',varid))
            call check_nc(nf90_get_var(initid,varid,temp4d,start=[srtx,srty,1,1],count=[cntx,cnty,icc,nmos]))
            ailcmaxrow = reshape(temp4d,[nlat,nmos,icc])

            call check_nc(nf90_inq_varid(initid,'dvdfcan',varid))
            call check_nc(nf90_get_var(initid,varid,temp4d,start=[srtx,srty,1,1],count=[cntx,cnty,icc,nmos]))
            dvdfcanrow = reshape(temp4d,[nlat,nmos,icc])

            !>Rest of the initialization variables are needed to run CTEM but if starting from bare ground initialize all
            !>live and dead c pools from zero. suitable values of extnprobgrd and prbfrhucgrd would still be required. set
            !>stdaln to 1 for operation in non-gcm stand alone mode, in the CTEM initialization file.

            call check_nc(nf90_inq_varid(initid,'gleafmas',varid))
            call check_nc(nf90_get_var(initid,varid,temp4d,start=[srtx,srty,1,1],count=[cntx,cnty,icc,nmos]))
            gleafmasrow = reshape(temp4d,[nlat,nmos,icc])

            call check_nc(nf90_inq_varid(initid,'bleafmas',varid))
            call check_nc(nf90_get_var(initid,varid,temp4d,start=[srtx,srty,1,1],count=[cntx,cnty,icc,nmos]))
            bleafmasrow = reshape(temp4d,[nlat,nmos,icc])

            call check_nc(nf90_inq_varid(initid,'stemmass',varid))
            call check_nc(nf90_get_var(initid,varid,temp4d,start=[srtx,srty,1,1],count=[cntx,cnty,icc,nmos]))
            stemmassrow = reshape(temp4d,[nlat,nmos,icc])

            call check_nc(nf90_inq_varid(initid,'rootmass',varid))
            call check_nc(nf90_get_var(initid,varid,temp4d,start=[srtx,srty,1,1],count=[cntx,cnty,icc,nmos]))
            rootmassrow = reshape(temp4d,[nlat,nmos,icc])

            !>If fire and competition are on, save the stemmass and rootmass for use in burntobare subroutine on the first timestep.
            if (dofire .and. compete) then
                        do j =1,icc
                            pstemmassrow(i,m,j)=stemmassrow(i,m,j)
                            pgleafmassrow(i,m,j)=rootmassrow(i,m,j)
                        end do
            end if

            call check_nc(nf90_inq_varid(initid,'litrmass',varid))
            call check_nc(nf90_get_var(initid,varid,temp4d,start=[srtx,srty,1,1],count=[cntx,cnty,icc,nmos]))
            litrmassrow = reshape(temp4d,[nlat,nmos,icc])

            call check_nc(nf90_inq_varid(initid,'soilcmas',varid))
            call check_nc(nf90_get_var(initid,varid,temp4d,start=[srtx,srty,1,1],count=[cntx,cnty,icc,nmos]))
            soilcmasrow = reshape(temp4d,[nlat,nmos,icc])

            call check_nc(nf90_inq_varid(initid,'lfstatus',varid))
            call check_nc(nf90_get_var(initid,varid,temp4d,start=[srtx,srty,1,1],count=[cntx,cnty,icc,nmos]))
            lfstatusrow = reshape(temp4d,[nlat,nmos,icc])

            call check_nc(nf90_inq_varid(initid,'pandays',varid))
            call check_nc(nf90_get_var(initid,varid,temp4d,start=[srtx,srty,1,1],count=[cntx,cnty,icc,nmos]))
            pandaysrow = reshape(temp4d,[nlat,nmos,icc])

            !call check_nc(nf90_inq_varid(initid,'stdaln',varid)) !FLAG - get rid of this? I don't see the use of stdaln
            !call check_nc(nf90_get_var(initid,varid,stdaln,start=[srty,srtx],count=[icnty,cntx]))

            deallocate(temp4d)

            if (compete .and. inibioclim) then  !read in the bioclimatic parameters

                allocate(temptwod(cntx,cnty))

                call check_nc(nf90_inq_varid(initid,'twarmm',varid))
                call check_nc(nf90_get_var(initid,varid,temptwod,start=[srtx,srty],count=[cntx,cnty]))
                twarmm = reshape(temptwod,[ nlat ] )

                call check_nc(nf90_inq_varid(initid,'tcoldm',varid))
                call check_nc(nf90_get_var(initid,varid,temptwod,start=[srtx,srty],count=[cntx,cnty]))
                tcoldm = reshape(temptwod,[ nlat ] )

                call check_nc(nf90_inq_varid(initid,'gdd5',varid))
                call check_nc(nf90_get_var(initid,varid,temptwod,start=[srtx,srty],count=[cntx,cnty]))
                gdd5 = reshape(temptwod,[ nlat ] )

                call check_nc(nf90_inq_varid(initid,'aridity',varid))
                call check_nc(nf90_get_var(initid,varid,temptwod,start=[srtx,srty],count=[cntx,cnty]))
                aridity = reshape(temptwod,[ nlat ] )

                call check_nc(nf90_inq_varid(initid,'srplsmon',varid))
                call check_nc(nf90_get_var(initid,varid,temptwod,start=[srtx,srty],count=[cntx,cnty]))
                srplsmon = reshape(temptwod,[ nlat ] )

                call check_nc(nf90_inq_varid(initid,'defctmon',varid))
                call check_nc(nf90_get_var(initid,varid,temptwod,start=[srtx,srty],count=[cntx,cnty]))
                defctmon = reshape(temptwod,[ nlat ] )

                call check_nc(nf90_inq_varid(initid,'anndefct',varid))
                call check_nc(nf90_get_var(initid,varid,temptwod,start=[srtx,srty],count=[cntx,cnty]))
                anndefct = reshape(temptwod,[ nlat ] )

                call check_nc(nf90_inq_varid(initid,'annsrpls',varid))
                call check_nc(nf90_get_var(initid,varid,temptwod,start=[srtx,srty],count=[cntx,cnty]))
                annsrpls = reshape(temptwod,[ nlat ] )

                call check_nc(nf90_inq_varid(initid,'annpcp',varid))
                call check_nc(nf90_get_var(initid,varid,temptwod,start=[srtx,srty],count=[cntx,cnty]))
                annpcp = reshape(temptwod,[ nlat ] )

                call check_nc(nf90_inq_varid(initid,'dry_season_length',varid))
                call check_nc(nf90_get_var(initid,varid,temptwod,start=[srtx,srty],count=[cntx,cnty]))
                dry_season_length = reshape(temptwod,[ nlat ] )

                deallocate(temptwod)

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

            !FLAG why are these vars per tile??!
            ! This below should now not be needed.
            ! !>Take the first tile value now and put it over the other tiles
            !         if (nmtest > 1) then
            !             do m = 2,nmtest
            !                 twarmm(i,m)=twarmm(i,1)
            !                 tcoldm(i,m)=tcoldm(i,1)
            !                 gdd5(i,m)=gdd5(i,1)
            !                 aridity(i,m)=aridity(i,1)
            !                 srplsmon(i,m)=srplsmon(i,1)
            !                 defctmon(i,m)=defctmon(i,1)
            !                 anndefct(i,m)=anndefct(i,1)
            !                 annsrpls(i,m)=annsrpls(i,1)
            !                 annpcp(i,m)=annpcp(i,1)
            !                 dry_season_length(i,m) =dry_season_length(i,1)
            !                 mlightng(i,m,:) = mlightng(i,1,:)
            !                 extnprob(i,m) = extnprob(i,1)
            !                 prbfrhuc(i,m) = prbfrhuc(i,1)
            !                 stdaln(i,m) = stdaln(i,1)
            !             end do
            !         end if

            if (dowetlands) then !if true then read wetland fractions into the first tile position

                allocate(temp3d(cntx,cnty,8))

                call check_nc(nf90_inq_varid(initid,'slopefrac',varid))
                call check_nc(nf90_get_var(initid,varid,temp3d,start=[srtx,srty,1],count=[cntx,cnty,8]))
                slopefrac = reshape(temp3d,[ nlat, 8 ] )

                deallocate(temp3d)

                ! FLAG assumedly this is no longer required.
!                 if (nmtest > 1) then ! if more tiles then just put the first tile across the rest FLAG
!                     do m = 2,nmtest
!                         slopefrac(i,m,:) = slopefrac(i,1,:)
!                     end do
!                 end if
            endif

            !>if this run uses the competition or lnduseon parameterization and starts from bare ground, set up the model state here. this
            !>overwrites what was read in from the initialization file. for composite runs (the composite set up is after this one for mosaics)
            if ((compete .or. lnduseon) .and. start_bare) then

                !>set up for composite runs when start_bare is on and compete or landuseon

                !>store the read-in crop fractions as we keep them even when we start bare.
                !!FLAG: this is setup assuming that crops are in pft number 6 and 7.
                !!and the first tile contains the information for the grid cell (assumes we have crops in
                !!every tile too! JM Apr 9 2014.
                do i=1,nltest
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

                do i=1,nlat  !FLAG change to nlat
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

end module model_state_drivers
