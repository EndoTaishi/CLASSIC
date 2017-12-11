!> Converts old format INI and CTM files into netcdf format initialization files.
!! Also takes in namelist format. Useful for tests at the site-level.

program initFileConverter

    use fileIOModule

    implicit none

    real                    :: lon, lat
    integer                 :: NLTEST, NMTEST
    character(300)          :: INIFile, CTMFile,charLon, charLat
    character(3)            :: fileType
    integer                 :: fileId,slopeDimId,iccp1DimId,icctemDimId,icDimId,monthsDimId,layerDimId, &
                                icp1DimId,lonDimId,latDimId,tileDimId
    character(50)           :: units
    character(350)          :: long_name
    integer, allocatable    :: dimArray(:),start(:),count(:)
    real :: dummy       ! dummyvar

    logical                 :: inclCTEM = .false.   ! True if a .CTM file was provided to the program
    integer                 :: ignd = 3             ! Number of soil layers in INI file.
    integer                 :: ican = 4
    integer                 :: icc = 9
    real, parameter         :: fillValue = -999.
    real                    :: grclarea = 100. !?

    ! All variable descriptions are listed in exportData
    real                    :: ZRFHROW
    real                    :: ZRFMROW
    real                    :: ZBLDROW
    real                    :: GCROW
    real, allocatable, dimension(:,:):: FCANROT
    real, allocatable, dimension(:)  :: FAREROT
    real, allocatable, dimension(:,:):: RSMNROT
    real, allocatable, dimension(:,:):: QA50ROT
    real, allocatable, dimension(:,:):: VPDAROT
    real, allocatable, dimension(:,:):: VPDBROT
    real, allocatable, dimension(:,:):: PSGAROT
    real, allocatable, dimension(:,:):: PSGBROT
    real, allocatable, dimension(:,:):: ALVCROT
    real, allocatable, dimension(:,:):: ALICROT
    real, allocatable, dimension(:,:):: PAMNROT
    real, allocatable, dimension(:,:):: PAMXROT
    real, allocatable, dimension(:,:):: LNZ0ROT
    real, allocatable, dimension(:,:):: CMASROT
    real, allocatable, dimension(:,:):: ROOTROT
    real, allocatable, dimension(:)  :: DRNROT
    real, allocatable, dimension(:)  :: SDEPROT
    real, allocatable, dimension(:)  :: XSLPROT
    real, allocatable, dimension(:)  :: GRKFROT
    real, allocatable, dimension(:)  :: WFSFROT
    real, allocatable, dimension(:)  :: WFCIROT
    integer, allocatable, dimension(:)  :: MIDROT
    real, allocatable, dimension(:,:):: SANDROT
    real, allocatable, dimension(:,:):: CLAYROT
    real, allocatable, dimension(:,:):: ORGMROT
    real, allocatable, dimension(:,:):: TBARROT
    real, allocatable, dimension(:,:):: THLQROT
    real, allocatable, dimension(:,:):: THICROT
    real, allocatable, dimension(:,:)   :: DELZ
    real, allocatable, dimension(:,:)   :: ZBOT
    real, allocatable, dimension(:)  :: TCANROT
    real, allocatable, dimension(:)  :: TSNOROT
    real, allocatable, dimension(:)  :: TPNDROT
    real, allocatable, dimension(:)  :: ZPNDROT
    real, allocatable, dimension(:)  :: RCANROT
    real, allocatable, dimension(:)  :: SCANROT
    real, allocatable, dimension(:)  :: SNOROT
    real, allocatable, dimension(:)  :: ALBSROT
    real, allocatable, dimension(:)  :: RHOSROT
    real, allocatable, dimension(:)  :: GROROT
    integer, allocatable, dimension(:)  :: SOCIROT

    real :: extnprob
    real :: prbfrhuc
    real, allocatable, dimension(:,:):: ailcminrow
    real, allocatable, dimension(:,:):: ailcmaxrow
    real, allocatable, dimension(:,:):: dvdfcanrow
    real, allocatable, dimension(:,:):: fcancmxrow
    real, allocatable, dimension(:,:):: gleafmasrow
    real, allocatable, dimension(:,:):: bleafmasrow
    real, allocatable, dimension(:,:):: stemmassrow
    real, allocatable, dimension(:,:):: rootmassrow
    real, allocatable, dimension(:,:):: pstemmassrow
    real, allocatable, dimension(:,:):: pgleafmassrow
    real, allocatable, dimension(:):: twarmm
    real, allocatable, dimension(:):: tcoldm
    real, allocatable, dimension(:):: gdd5
    real, allocatable, dimension(:):: aridity
    real, allocatable, dimension(:):: srplsmon
    real, allocatable, dimension(:):: defctmon
    real, allocatable, dimension(:):: anndefct
    real, allocatable, dimension(:):: annsrpls
    real, allocatable, dimension(:):: annpcp
    real, allocatable, dimension(:):: dry_season_length
    real, allocatable, dimension(:,:):: litrmassrow
    real, allocatable, dimension(:,:):: soilcmasrow
    real, allocatable, dimension(:):: mlightng
    integer, allocatable, dimension(:,:):: lfstatusrow
    integer, allocatable, dimension(:,:):: pandaysrow
    real, allocatable, dimension(:,:):: slopefrac
    integer, allocatable, dimension(:):: ipeatlandrow   !<Peatland switch: 0 = not a peatland, 1= bog, 2 = fen
    real, allocatable, dimension(:):: Cmossmas          !<C in moss biomass, \f$kg C/m^2\f$
    real, allocatable, dimension(:):: litrmsmoss        !<moss litter mass, \f$kg C/m^2\f$
    real, allocatable, dimension(:):: dmoss             !<depth of living moss (m)
    real, allocatable, dimension(:) :: grclarea         !<area of the grid cell, \f$km^2\f$

    !----------

    ! Parse the arguments and determine if the input file is INI format or namelist
    call processArguments

    ! Read the file headers (INI or nml) so we can allocate arrays
    if (fileType == 'ini') then
        ! Open the INI file, or the namelist file
        open(unit = 10, file = INIFile, form = 'formatted', status = 'old', action = 'read')

        read(10,*) ! Throw out first three lines.
        read(10,*)
        read(10,*)
        READ(10,5020) lat,lon,ZRFMROW,ZRFHROW,ZBLDROW,GCROW,NLTEST,NMTEST
    else ! the file is a namelist
        call readNamelistHeader
    end if

    ! Allocate the arrays
    call setupArrays

    ! Read in the rest of the input from the INI/nml file
    if (fileType == 'ini') then
        call loadINIData
    else
        call readNamelist
    end if

    ! Read in the CTM file
    if (inclCTEM) then
        call loadCTMData
    end if

    ! Make the file
    call makeNetCDF

    ! Write the data out to netcdf
    call exportData

    ! Close the input and output files
    close(unit = 10)
    call ncClose(fileId)
    if (inclCTEM) close(unit = 20)

    5020  FORMAT(5F10.2,F7.1,3I5)

contains

    ! ------------------------------------------------------------------------------------------------------

    subroutine processArguments

        call getarg(1, INIFile)

        if (iargc() .eq. 2) then
            call getarg(2, CTMFile)
            if (index(CTMFile,'CTM') > 0 .or. index(CTMFile,'ctm') > 0) then
                inclCTEM = .true.
            end if
        elseif (iargc() .gt. 2 .or. iargc() == 0) then
            print*
            print*,'Expecting only one or two arguments: First in the INI or namelist file'
            print*,'the second file is the CTM file and only required if you want to run with'
            print*,'CTEM on. The possible suffixes of the input files are:.INI, .ini, .CTM, .ctm, .nml, .txt'
            print*,'Output netcdf is given the same name as the INI/nml file with the suffix .nc'
        endif

        ! Parse the INIfile to help decide if it is a namelist or INI file
        if (index(INIFile,'INI') > 0 .or. index(INIFile,'ini') > 0) then
            fileType = 'ini'
        elseif (index(INIFile,'txt') > 0 .or. index(INIFile,'nml') > 0) then
            fileType = 'nml'
        else
            print*,'***Unknown file given for first file. Expected suffixes: INI,ini,txt,nml'
            call exit
        end if

    end subroutine processArguments

    ! ------------------------------------------------------------------------------------------------------

    subroutine setupArrays

        allocate(FCANROT(NMTEST,ICAN+1))
        allocate(PAMXROT(NMTEST,ICAN))
        allocate(LNZ0ROT(NMTEST,ICAN+1))
        allocate(PAMNROT(NMTEST,ICAN))
        allocate(ALVCROT(NMTEST,ICAN+1))
        allocate(CMASROT(NMTEST,ICAN))
        allocate(ALICROT(NMTEST,ICAN+1))
        allocate(ROOTROT(NMTEST,ICAN))
        allocate(RSMNROT(NMTEST,ICAN))
        allocate(QA50ROT(NMTEST,ICAN))
        allocate(VPDAROT(NMTEST,ICAN))
        allocate(VPDBROT(NMTEST,ICAN))
        allocate(PSGAROT(NMTEST,ICAN))
        allocate(PSGBROT(NMTEST,ICAN))
        allocate(DRNROT(NMTEST))
        allocate(SDEPROT(NMTEST))
        allocate(FAREROT(NMTEST))
        allocate(XSLPROT(NMTEST))
        allocate(GRKFROT(NMTEST))
        allocate(WFSFROT(NMTEST))
        allocate(WFCIROT(NMTEST))
        allocate(MIDROT(NMTEST))
        allocate(SOCIROT(NMTEST))
        allocate(SANDROT(NMTEST,ignd))
        allocate(CLAYROT(NMTEST,ignd))
        allocate(ORGMROT(NMTEST,ignd))
        allocate(TBARROT(NMTEST,ignd))
        allocate(TCANROT(NMTEST))
        allocate(TSNOROT(NMTEST))
        allocate(TPNDROT(NMTEST))
        allocate(THLQROT(NMTEST,ignd))
        allocate(THICROT(NMTEST,ignd))
        allocate(ZPNDROT(NMTEST))
        allocate(RCANROT(NMTEST))
        allocate(SCANROT(NMTEST))
        allocate(SNOROT(NMTEST))
        allocate(ALBSROT(NMTEST))
        allocate(RHOSROT(NMTEST))
        allocate(GROROT(NMTEST))
        allocate(DELZ(nmtest,ignd))
        allocate(ZBOT(nmtest,ignd))

        allocate(ailcminrow(NMTEST,icc))
        allocate(ailcmaxrow(NMTEST,icc))
        allocate(dvdfcanrow(NMTEST,icc))
        allocate(gleafmasrow(NMTEST,icc))
        allocate(bleafmasrow(NMTEST,icc))
        allocate(stemmassrow(NMTEST,icc))
        allocate(rootmassrow(NMTEST,icc))
        allocate(litrmassrow(NMTEST,icc+1))
        allocate(soilcmasrow(NMTEST,icc+1))
        allocate(lfstatusrow(NMTEST,icc))
        allocate(pandaysrow(NMTEST,icc))
        allocate(fcancmxrow(nmtest,icc))
        allocate(mlightng(12))

    end subroutine setupArrays

    ! ------------------------------------------------------------------------------------------------------

    subroutine loadINIData()

        integer                 :: m,j

         DO 50 M=1,NMTEST
          READ(10,5040) (FCANROT(M,J),J=1,ICAN+1),(PAMXROT(M,J),J=1,ICAN)
          READ(10,5040) (LNZ0ROT(M,J),J=1,ICAN+1),(PAMNROT(M,J),J=1,ICAN)
          READ(10,5040) (ALVCROT(M,J),J=1,ICAN+1),(CMASROT(M,J),J=1,ICAN)
          READ(10,5040) (ALICROT(M,J),J=1,ICAN+1),(ROOTROT(M,J),J=1,ICAN)
          READ(10,5030) (RSMNROT(M,J),J=1,ICAN),(QA50ROT(M,J),J=1,ICAN)
          READ(10,5030) (VPDAROT(M,J),J=1,ICAN),(VPDBROT(M,J),J=1,ICAN)
          READ(10,5030) (PSGAROT(M,J),J=1,ICAN),(PSGBROT(M,J),J=1,ICAN)
          READ(10,5040) DRNROT(M),SDEPROT(M),FAREROT(M)
          READ(10,5090) XSLPROT(M),GRKFROT(M),WFSFROT(M),WFCIROT(M),MIDROT(M),SOCIROT(M)
          READ(10,5080) (SANDROT(M,J),J=1,ignd)
          READ(10,5080) (CLAYROT(M,J),J=1,ignd)
          READ(10,5080) (ORGMROT(M,J),J=1,ignd)
          READ(10,5050) (TBARROT(M,J),J=1,ignd),TCANROT(M),TSNOROT(M),TPNDROT(M)
          READ(10,5060) (THLQROT(M,J),J=1,ignd),(THICROT(M,J),J=1,ignd),ZPNDROT(M)
          READ(10,5070) RCANROT(M),SCANROT(M),SNOROT(M),ALBSROT(M),RHOSROT(M),GROROT(M)
50    CONTINUE
     ! In CLASS 3.6.2, we include this soil info in the INI file.
      DO 25 J=1,IGND
          READ(10,*) DELZ(1,J),ZBOT(1,J)
 25   CONTINUE

    ! DELZ and ZBOT get put the same for all tiles.
    if (nmtest > 1) then
        do m = 2,nmtest
            delz(m,:) = delz(1,:)
            zbot(m,:) = zbot(1,:)
        end do
    end if

5040  FORMAT(9F8.3)
5030  FORMAT(4F8.3,8X,4F8.3)
5050  FORMAT(6F10.2)
5060  FORMAT(7F10.3)
5070  FORMAT(2F10.4,F10.2,F10.3,F10.4,F10.3)
5080  FORMAT(3F10.1)
5090  FORMAT(4E8.1,2I8)

    end subroutine loadINIData

    ! ------------------------------------------------------------------------------------------------------

    subroutine loadCTMData()

        integer :: m,j

        !>Read from CTEM initialization file (.CTM)
        open(unit = 11, file = CTMFile, form = 'formatted', status = 'old', action = 'read')

        ! dump header
        read(11,*)
        read(11,*)
        read(11,*)

        do m=1,nmtest
            read(11,*) (ailcminrow(m,j),j=1,icc)
            read(11,*) (ailcmaxrow(m,j),j=1,icc)
            read(11,*) (dvdfcanrow(m,j),j=1,icc)
            read(11,*) (gleafmasrow(m,j),j=1,icc)
            read(11,*) (bleafmasrow(m,j),j=1,icc)
            read(11,*) (stemmassrow(m,j),j=1,icc)
            read(11,*) (rootmassrow(m,j),j=1,icc)
            read(11,*) (litrmassrow(m,j),j=1,icc+1)
            read(11,*) (soilcmasrow(m,j),j=1,icc+1)
            read(11,*) (lfstatusrow(m,j),j=1,icc)
            read(11,*) (pandaysrow(m,j),j=1,icc)
        end do
            read(11,*) (mlightng(j),j=1,6)  !mean monthly lightning frequency
            read(11,*) (mlightng(j),j=7,12) !flashes/km2.year, this is spread over other tiles below
            read(11,*) extnprob
            read(11,*) prbfrhuc
            read(11,*) dummy ! was stdaln but not used so dump.
            !if (compete .and. inibioclim) then  !read in the bioclimatic parameters
            !! read them into the first tile of each grid cell.
            !read(11,*) twarmm(i,1), tcoldm(i,1), gdd5(i,1), aridity(i,1),srplsmon(i,1)
            !read(11,*) defctmon(i,1), anndefct(i,1), annsrpls(i,1), annpcp(i,1), dry_season_length(i,1)


    end subroutine loadCTMData

    ! ------------------------------------------------------------------------------------------------------

    subroutine readNamelistHeader

        namelist /header/ &
        lat,&
        lon,&
        nmtest,&
        ignd,&
        ican,&
        icc

        read(unit=10,nml = header)
        rewind(10)

    end subroutine readNamelistHeader

    ! ------------------------------------------------------------------------------------------------------

    subroutine readNamelist

        namelist /classicvars/ &
            ZRFHROW,&
            ZRFMROW,&
            ZBLDROW,&
            GCROW,&
            FCANROT,&
            FAREROT,&
            RSMNROT,&
            QA50ROT,&
            VPDAROT,&
            VPDBROT,&
            PSGAROT,&
            PSGBROT,&
            ALVCROT,&
            ALICROT,&
            PAMNROT,&
            PAMXROT,&
            LNZ0ROT,&
            CMASROT,&
            ROOTROT,&
            DRNROT,&
            SDEPROT,&
            XSLPROT,&
            GRKFROT,&
            WFSFROT,&
            WFCIROT,&
            MIDROT,&
            SANDROT,&
            CLAYROT,&
            ORGMROT,&
            TBARROT,&
            THLQROT,&
            THICROT,&
            DELZ,&
            ZBOT,&
            TCANROT,&
            TSNOROT,&
            TPNDROT,&
            ZPNDROT,&
            RCANROT,&
            SCANROT,&
            SNOROT,&
            ALBSROT,&
            RHOSROT,&
            GROROT,&
            SOCIROT

        read(unit=10,nml = classicvars)
        !write(*,nml = classicvars)

    end subroutine readNamelist

    ! ------------------------------------------------------------------------------------------------------

    subroutine makeNetCDF

        character(400) :: filename,title
        character(8)   :: date
        integer :: i,varId
        integer :: tile(nmtest)
        integer :: icp1(ican+1)
        integer :: layer(ignd)
        integer :: ic(ican)
        integer :: icctem(icc)
        integer :: iccp1(icc+1)
        integer :: months(12)
        integer :: slope(8)
        integer :: indexend

        tile = (/(i, i=1,nmtest, 1)/)
        icp1 = (/(i, i=1,ican+1, 1)/)
        layer = (/(i, i=1,ignd, 1)/)
        ic = (/(i, i=1,ican, 1)/)
        icctem = (/(i, i=1,icc, 1)/)
        iccp1 = (/(i, i=1,icc+1, 1)/)
        months = (/(i, i=1,12, 1)/)
        slope = (/(i, i=1,8, 1)/)

        ! Filename is going to be the INI file name with .nc as a suffix.
        if (fileType == 'ini') then
            indexend = max(index(INIFile,'INI'),index(INIFile,'ini'))
        else
            indexend = max(index(INIFile,'nml'),index(INIFile,'txt'))
        end if
        filename = INIfile(1:indexend-1)//'nc'

        ! If the file doesn't already exist, then initialize the file
        !if (.not.fileExists(trim(filename))) then
            ! Create the values file
            fileId = ncCreate(filename, NF90_CLOBBER)
        !else
        !    fileId = ncOpen(filename, NF90_WRITE)
        !end if

        ! Add in the metadata for the file

        title = 'CLASSIC initialization file created from: '//INIfile
        if (inclCTEM) then
            title = title//' and '//CTMFile
        end if

        call ncPutAtt(fileId, nf90_global, 'title', charValues = trim(title))
        call date_and_time(date=date)
        call ncPutAtt(fileId, nf90_global, 'creation_date', charValues = date)

        ! Define the longitude dimension
        lonDimId = ncDefDim(fileId, 'lon', 1)
        varid = ncDefVar(fileId, 'lon', nf90_double, [lonDimId])
        call ncPutAtt(fileId, varId, 'standard_name', charValues = 'Longitude')
        call ncPutAtt(fileId, varId, 'units', charValues = 'degrees_east')
        call ncPutAtt(fileId, varId, 'axis', charValues = 'X')
        call ncEndDef(fileId)
        call ncPutDimValues(fileId, 'lon', realValues = [lon], count = [1])

        call ncRedef(fileId)

        ! Define the latitude dimension
        latDimId = ncDefDim(fileId, 'lat', 1)
        varid = ncDefVar(fileId, 'lat', nf90_double, [latDimId])
        call ncPutAtt(fileId, varId, 'standard_name', charValues = 'Latitude')
        call ncPutAtt(fileId, varId, 'units', charValues = 'degrees_north')
        call ncPutAtt(fileId, varId, 'axis', charValues = 'Y')
        call ncEndDef(fileId)
        call ncPutDimValues(fileId, 'lat', [lat], count = [1])

        call ncRedef(fileId)

        ! Define the tile dimension
        tileDimId = ncDefDim(fileId, 'tile', size(tile))
        varid = ncDefVar(fileId, 'tile', nf90_int, [tileDimId])
        call ncPutAtt(fileId, varId, 'standard_name', charValues = 'land surface tile')
        call ncEndDef(fileId)
        call ncPutDimValues(fileId, 'tile', intValues=tile, count = [size(tile)])

        call ncRedef(fileId)

        ! Define the icp1 dimension
        icp1DimId = ncDefDim(fileId, 'icp1', size(icp1))
        varid = ncDefVar(fileId, 'icp1', nf90_int, [icp1DimId])
        call ncPutAtt(fileId, varId, 'standard_name', charValues = 'physics (CLASS) PFTs plus bareground')
        call ncEndDef(fileId)
        call ncPutDimValues(fileId, 'icp1', intValues=icp1, count = [size(icp1)])

        call ncRedef(fileId)

        ! Define the layer dimension
        layerDimId = ncDefDim(fileId, 'layer', size(layer))
        varid = ncDefVar(fileId, 'layer', nf90_int, [layerDimId])
        call ncPutAtt(fileId, varId, 'standard_name', charValues = 'ground layers')
        call ncEndDef(fileId)
        call ncPutDimValues(fileId, 'layer', intValues=layer, count = [size(layer)])

        call ncRedef(fileId)

        ! Define the ic dimension
        icDimId = ncDefDim(fileId, 'ic', size(ic))
        varid = ncDefVar(fileId, 'ic', nf90_int, [icDimId])
        call ncPutAtt(fileId, varId, 'standard_name', charValues = 'physics (CLASS) PFTs')
        call ncEndDef(fileId)
        call ncPutDimValues(fileId, 'ic', intValues=ic, count = [size(ic)])

        call ncRedef(fileId)

        ! Define the icctem dimension
        icctemDimId = ncDefDim(fileId, 'icctem', size(icctem))
        varid = ncDefVar(fileId, 'icctem', nf90_int, [icctemDimId])
        call ncPutAtt(fileId, varId, 'standard_name', charValues = 'biogeochemical (CTEM) PFTs')
        call ncEndDef(fileId)
        call ncPutDimValues(fileId, 'icctem', intValues=icctem, count = [size(icctem)])

        call ncRedef(fileId)

        ! Define the iccp1 dimension
        iccp1DimId = ncDefDim(fileId, 'iccp1', size(iccp1))
        varid = ncDefVar(fileId, 'iccp1', nf90_int, [iccp1DimId])
        call ncPutAtt(fileId, varId, 'standard_name', charValues = 'biogeochemical (CTEM) PFTs plus bareground')
        call ncEndDef(fileId)
        call ncPutDimValues(fileId, 'iccp1', intValues=iccp1, count = [size(iccp1)])

        call ncRedef(fileId)

        ! Define the months dimension
        monthsDimId = ncDefDim(fileId, 'months', size(months))
        varid = ncDefVar(fileId, 'months', nf90_int, [monthsDimId])
        call ncPutAtt(fileId, varId, 'standard_name', charValues = 'months')
        call ncEndDef(fileId)
        call ncPutDimValues(fileId, 'months', intValues=months, count = [size(months)])

        call ncRedef(fileId)

        ! Define the slope dimension
        slopeDimId = ncDefDim(fileId, 'slope', size(slope))
        varid = ncDefVar(fileId, 'slope', nf90_int, [slopeDimId])
        call ncPutAtt(fileId, varId, 'standard_name', charValues = 'wetland slope fractions for 0.025, 0.05, 0.1, 0.15, 0.20, 0.25, 0.3 and 0.35 percent slope threshold')
        call ncEndDef(fileId)
        call ncPutDimValues(fileId, 'slope', intValues=slope, count = [size(slope)])

        ! Add some generic variables that CLASSIC will look for

        ! per grid variables
        allocate(dimArray(2),start(2),count(2))
        dimArray = (/lonDimId,latDimId/)
        start = (/1, 1 /)
        count = (/1, 1 /)
        call exportVariable('MID',units='-',long_name='Mosaic tile type identifier (1 for land surface, 0 for inland lake)',intvalues=MIDROT)
        call exportVariable('GC',units='-',long_name='GCM surface descriptor - land surfaces (inc. inland water) is -1',intvalues=(/-1/))
        call exportVariable('nmtest',units='-',long_name='Number of tiles in each grid cell',intvalues=tile)
        deallocate(dimArray,start,count)

    end subroutine makeNetCDF

    ! ------------------------------------------------------------------------------------------------------

    subroutine exportData

        integer :: m

        ! icp1 variables:
        allocate(dimArray(4),start(4),count(4))
        dimArray = (/lonDimId,latDimId,icp1DimId,tileDimId/)
        start = (/1, 1, 1 ,1/)
        count = (/1, 1, ican+1, 1/)
        call exportVariable('FCAN',units='-',long_name='Annual maximum fractional coverage of modelled area (read in for CLASS only runs)',values2D=FCANROT)
        call exportVariable('LNZ0',units='-',long_name='Natural logarithm of maximum vegetation roughness length',values2D=LNZ0ROT)
        call exportVariable('ALIC',units='-',long_name='Average near-IR albedo of vegetation category when fully-leafed',values2D=ALICROT)
        call exportVariable('ALVC',units='-',long_name='Average visible albedo of vegetation category when fully-leafed',values2D=ALVCROT)


        ! ic variables:
        dimArray = (/lonDimId,latDimId,icDimId,tileDimId/)
        count = (/1, 1, ican, 1/)
        call exportVariable('PAMX',units='m2/m2',long_name='Annual maximum plant area index of vegetation category',values2D=PAMXROT)
        call exportVariable('PAMN',units='m2/m2',long_name='Annual minimum plant area index of vegetation category',values2D=PAMNROT)
        call exportVariable('CMAS',units='$[kg m^{-2} ]$',long_name='Annual maximum canopy mass for vegetation category',values2D=CMASROT)
        call exportVariable('ROOT',units='m',long_name='Annual maximum rooting depth of vegetation category',values2D=ROOTROT)
        call exportVariable('RSMN',units='s/m',long_name='Minimum stomatal resistance of vegetation category',values2D=RSMNROT)
        call exportVariable('QA50',units='W/m2',long_name='Reference value of incoming shortwave radiation (used in stomatal resistance calculation)',values2D=QA50ROT)
        call exportVariable('VPDA',units='-',long_name='Vapour pressure deficit coefficient (used in stomatal resistance calculation)',values2D=VPDAROT)
        call exportVariable('VPDB',units='-',long_name='Vapour pressure deficit coefficient (used in stomatal resistance calculation)',values2D=VPDBROT)
        call exportVariable('PSGA',units='-',long_name='Soil moisture suction coefficient (used in stomatal resistance calculation)',values2D=PSGAROT)
        call exportVariable('PSGB',units='-',long_name='Soil moisture suction coefficient (used in stomatal resistance calculation)',values2D=PSGBROT)

        ! ignd variables:
        dimArray = (/lonDimId,latDimId,layerDimId,tileDimId/)
        count = (/1, 1, ignd, 1/)
        call exportVariable('SAND',units='%',long_name='Percentage sand content',values2D=SANDROT)
        call exportVariable('CLAY',units='%',long_name='Percentage clay content',values2D=CLAYROT)
        call exportVariable('ORGM',units='%',long_name='Percentage organic matter content',values2D=ORGMROT)
        call exportVariable('TBAR',units='C',long_name='Temperature of soil layers',values2D=TBARROT)
        call exportVariable('THIC',units='m3/m3',long_name='Volumetric frozen water content of soil layers',values2D=THICROT)
        call exportVariable('THLQ',units='m3/m3',long_name='Volumetric liquid water content of soil layers',values2D=THLQROT)
        call exportVariable('DELZ',units='m',long_name='Ground layer thickness',values2D=DELZ)
        call exportVariable('ZBOT',units='m',long_name='Depth of bottom of ground layer',values2D=ZBOT)

        deallocate(dimArray,start,count)

        ! nmtest only variables:
        allocate(dimArray(3),start(3),count(3))
        dimArray = (/lonDimId,latDimId,tileDimId/)
        start = (/1, 1, 1 /)
        count = (/1, 1, nmtest/)
        call exportVariable('DRN',units='-',long_name='Soil drainage index',values=DRNROT)
        call exportVariable('FARE',units='fraction',long_name='Tile fractional area of gridcell',values=FAREROT)
        call exportVariable('SDEP',units='m',long_name='Soil permeable depth',values=SDEPROT)
        call exportVariable('XSLP',units='-',long_name='Not in Use: parameters lateral movement of soil water',values=XSLPROT)
        call exportVariable('GRKF',units='-',long_name='Not in Use: parameters lateral movement of soil water',values=GRKFROT)
        call exportVariable('WFCI',units='-',long_name='Not in Use: parameters lateral movement of soil water',values=WFCIROT)
        call exportVariable('WFSF',units='-',long_name='Not in Use: parameters lateral movement of soil water',values=WFSFROT)
        call exportVariable('SOCI',units='index',long_name='Soil colour index',intvalues=SOCIROT)
        call exportVariable('TCAN',units='C',long_name='Vegetation canopy temperature',values=TCANROT)
        call exportVariable('ALBS',units='-',long_name='Soil drainage index',values=ALBSROT)
        call exportVariable('GRO',units='-',long_name='Vegetation growth index',values=GROROT)
        call exportVariable('RCAN',units='-',long_name='Intercepted liquid water stored on canopy',values=RCANROT)
        call exportVariable('RHOS',units='kg/m3',long_name='Density of snow',values=RHOSROT)
        call exportVariable('SCAN',units='kg/m2',long_name='Intercepted frozen water stored on canopy',values=SCANROT)
        call exportVariable('SNO',units='kg/m2',long_name='Mass of snow pack',values=SNOROT)
        call exportVariable('TPND',units='C',long_name='Temperature of ponded water',values=TPNDROT)
        call exportVariable('TSNO',units='C',long_name='Snowpack temperature',values=TSNOROT)
        call exportVariable('ZPND',units='m',long_name='Depth of ponded water on surface',values=ZPNDROT)

        deallocate(dimArray,start,count)

        ! per grid variables
        allocate(dimArray(2),start(2),count(2))
        dimArray = (/lonDimId,latDimId/)
        start = (/1, 1 /)
        count = (/1, 1 /)
        call exportVariable('ZBLD',units='m',long_name='Atmospheric blending height for surface roughness length averaging',values=(/ZBLDROW/))
        call exportVariable('ZRFH',units='m',long_name='Reference height associated with forcing air temperature and humidity',values=(/ZRFHROW/))
        call exportVariable('ZRFM',units='m',long_name='Reference height associated with forcing wind speed',values=(/ZRFMROW/))

        deallocate(dimArray,start,count)

        if (inclCTEM) then

            ! icc variables
            allocate(dimArray(4),start(4),count(4))
            dimArray = (/lonDimId,latDimId,icctemDimId,tileDimId/)
            start = (/1, 1, 1 ,1/)
            count = (/1, 1, icc, 1/)
            call exportVariable('ailcmin',units='m2/m2',long_name='Min. LAI for use with CTEM1 option only. Obsolete',values2D=ailcminrow)
            call exportVariable('ailcmax',units='m2/m2',long_name='Max. LAI for use with CTEM1 option only. Obsolete',values2D=ailcmaxrow)
            call exportVariable('bleafmas',units='kgC/m2',long_name='Brown leaf mass',values2D=bleafmasrow)
            call exportVariable('gleafmas',units='kgC/m2',long_name='Green leaf mass',values2D=gleafmasrow)
            call exportVariable('stemmass',units='kgC/m2',long_name='Stem mass',values2D=stemmassrow)
            call exportVariable('rootmass',units='kgC/m2',long_name='Root mass',values2D=rootmassrow)
            call exportVariable('lfstatus',units='-',long_name='Leaf status, see Phenology',intvalues2D=lfstatusrow)
            call exportVariable('pandays',units='-',long_name='Days with +ve new photosynthesis, see Phenology',intvalues2D=pandaysrow)

            if (fileType == 'ini') then
                if (icc .ne. 9 and ican .ne. 4) print*,'Warning - expected ICC =9 and ICAN = 4'
                do m = 1,nmtest
                    fcancmxrow(m,1) = FCANROT(m,1) * dvdfcanrow(m,1)
                    fcancmxrow(m,2) = FCANROT(m,1) * dvdfcanrow(m,2)
                    fcancmxrow(m,3) = FCANROT(m,2) * dvdfcanrow(m,3)
                    fcancmxrow(m,4) = FCANROT(m,2) * dvdfcanrow(m,4)
                    fcancmxrow(m,5) = FCANROT(m,2) * dvdfcanrow(m,5)
                    fcancmxrow(m,6) = FCANROT(m,3) * dvdfcanrow(m,6)
                    fcancmxrow(m,7) = FCANROT(m,3) * dvdfcanrow(m,7)
                    fcancmxrow(m,8) = FCANROT(m,4) * dvdfcanrow(m,8)
                    fcancmxrow(m,9) = FCANROT(m,4) * dvdfcanrow(m,9)
                end do
            !else - namelist reads in the fcancmx directly without this conversion
            end if
            call exportVariable('fcancmx',units='-',long_name='PFT fractional coverage per grid cell',intvalues2D=fcancmxrow)

            ! iccp1 variables
            dimArray = (/lonDimId,latDimId,iccp1DimId,tileDimId/)
            count = (/1, 1, icc+1, 1/)
            call exportVariable('litrmass',units='kgC/m2',long_name='Litter mass per soil layer',values2D=litrmassrow)
            call exportVariable('soilcmas',units='kgC/m2',long_name='Soil C mass per soil layer',values2D=soilcmasrow)

            deallocate(dimArray,start,count)

            ! per month vars
            allocate(dimArray(3),start(3),count(3))
            dimArray = (/lonDimId,latDimId,monthsDimId/)
            start = (/1, 1, 1 /)
            count = (/1, 1, 12/)
            call exportVariable('mlightng',units='flashes/km2.year',long_name='mean monthly lightning freq. (total flashes)',values=mlightng)

            deallocate(dimArray,start,count)

            ! read(11,*) extnprob
            ! read(11,*) prbfrhuc
    !     float extnprob(lat, lon) ;
    !         extnprob:_FillValue = -999.f ;
    !         extnprob:units = "-" ;
    !         extnprob:long_name = "Fire extinguishing probability (overwritten if POPD true)" ;
    !     float prbfrhuc(lat, lon) ;
    !         prbfrhuc:_FillValue = -999.f ;
    !         prbfrhuc:units = "-" ;
    !         prbfrhuc:long_name = "Probability of fire due to human causes (overwritten if POPD true)" ;

    !     float grclarea(lat, lon) ;
    !         grclarea:_FillValue = -999.f ;
    !         grclarea:units = "km2" ;
    !         grclarea:long_name = "Area of grid cell" ;

    !     float ipeatland(tile, lat, lon) ;
    !         ipeatland:_FillValue = -999.f ;
    !         ipeatland:units = "-" ;
    !         ipeatland:long_name = "Peatland flag: 0 = not a peatland, 1= bog, 2 = fen" ;
    !     float dmoss(tile, lat, lon) ;
    !         dmoss:_FillValue = -999.f ;
    !         dmoss:units = "m" ;
    !         dmoss:long_name = "Depth of living moss" ;
    !     float litrmsmoss(tile, lat, lon) ;
    !         litrmsmoss:_FillValue = -999.f ;
    !         litrmsmoss:units = "kgC/m2" ;
    !         litrmsmoss:long_name = "Moss litter mass" ;
    !     float rice(months, lat, lon) ;
    !         rice:_FillValue = -999.f ;
    !         rice:units = "-" ;
    !         rice:long_name = "Monthly irrigated rice ag. gridcell fraction" ;
    !     float slopefrac(slope, tile, lat, lon) ;
    !         slopefrac:_FillValue = -999.f ;
    !         slopefrac:units = "-" ;
    !         slopefrac:long_name = "Slope-based fraction for dynamic wetlands" ;
    !     float Cmossmas(tile, lat, lon) ;
    !         Cmossmas:_FillValue = -999.f ;
    !         Cmossmas:units = "kgC/m2" ;
    !         Cmossmas:long_name = "C in moss biomass" ;

            end if

    end subroutine exportData

    ! ------------------------------------------------------------------------------------------------------

    subroutine exportVariable(name,units,long_name,values,values2D,intvalues,intvalues2D)
        character(*), intent(in)    :: name
        real, intent(in),optional   :: values(:),values2D(:,:)
        integer, intent(in), optional :: intvalues(:),intvalues2D(:,:)
        character(*), intent(in)    :: units
        character(*), intent(in)    :: long_name
        integer                     :: varId, m

        call ncRedef(fileId)

        varid = ncDefVar(fileId, name, nf90_double, dimArray)
        call ncPutAtt(fileId, varid, 'units', charValues = units)
        call ncPutAtt(fileId, varid, '_FillValue', realValues = fillValue)
        call ncPutAtt(fileId, varId, 'long_name', charValues = long_name)
        call ncEndDef(fileId)

        ! Put in data
        if (present(values)) then
            call ncPutVar(fileId, name, realValues = values, start= start, count = count)
        else if (present(values2D)) then
            do m = 1,nmtest
                call ncPutVar(fileId, name, realValues = values2D(m,:), start= start, count = count)
            end do
        else if (present(intvalues)) then
            call ncPutVar(fileId, name, intValues = intvalues, start= start, count = count)
        else if (present(intvalues2D)) then
            do m = 1,nmtest
                call ncPutVar(fileId, name, intValues = intvalues2D(m,:), start= start, count = count)
            end do
        else
            print*,'Problem in exportVariable'
        end if


    end subroutine exportVariable

    ! ------------------------------------------------------------------------------------------------------

    logical function fileExists(filename)
        character(*), intent(in) :: filename
        inquire(file = filename, exist = fileExists)
    end function fileExists

    ! ------------------------------------------------------------------------------------------------------

    real function charToReal(input)
        character(len=*), intent(in)    :: input    !< Char input
        read(input,*) charToReal
    end function charToReal
    ! ------------------------------------------------------------------------------------------------------

end program initFileConverter
