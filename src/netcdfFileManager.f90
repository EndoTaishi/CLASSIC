module netcdfFileManager
    use types
    use netcdf
    implicit none
    integer, allocatable    :: ncid(:)

contains
    ! Create the definitions for this variable
    subroutine createNetCDF(id, outputForm, descriptor)
        implicit none
        character(*), intent(in)                :: outputForm
        type(variableDescriptor), intent(in)    :: descriptor
        integer, intent(in)                     :: id
        character(8)                            :: today
        character(10)                           :: now
        integer                     :: ncid, varid, lat, lon, cntx = 10, cnty = 20, suffix
        integer                     :: nmos = 3, tile, icc = 10, pft, ignd = 21, layer, time, i
        real, dimension(2)          :: xrange
        real, dimension(2)          :: yrange
        integer, allocatable, dimension(:) :: pftnum
        character(30)               :: timestart = "seconds since 1801-1-1"
        logical                     :: leap = .true.

        ncid = variables(id)%ncid

        call check_nc(nf90_put_att(ncid,nf90_global,'title','CLASSIC output file'))

        call date_and_time(today,now)

        call check_nc(nf90_put_att(ncid,nf90_global,'timestamp',today//' '//now(1:4)))
        call check_nc(nf90_put_att(ncid,nf90_global,'Conventions','COARDS'))
        call check_nc(nf90_put_att(ncid,nf90_global,'node_offset',1))

        ! Longitude
        call check_nc(nf90_def_dim(ncid,'lon',cntx,lon))
        call check_nc(nf90_def_var(ncid,'lon',nf90_float,lon,varid))
        call check_nc(nf90_put_att(ncid,varid,'long_name','longitude'))
        call check_nc(nf90_put_att(ncid,varid,'units','degrees_east'))
        call check_nc(nf90_put_att(ncid,varid,'actual_range',xrange))
        call check_nc(nf90_put_att(ncid,varid,'_Storage',"contiguous"))

        ! Latitude
        call check_nc(nf90_def_dim(ncid,'lat',cnty,lat))
        call check_nc(nf90_def_var(ncid,'lat',nf90_float,lat,varid))
        call check_nc(nf90_put_att(ncid,varid,'long_name','latitude'))
        call check_nc(nf90_put_att(ncid,varid,'units','degrees_north'))
        call check_nc(nf90_put_att(ncid,varid,'actual_range',yrange))
        call check_nc(nf90_put_att(ncid,varid,'_Storage',"contiguous"))

        select case(trim(outputForm))
            case ("tile")       ! Per tile outputs
                call check_nc(nf90_def_dim(ncid,'tile',nmos,tile))
                call check_nc(nf90_def_var(ncid,'tile',nf90_short,tile,varid))
                call check_nc(nf90_put_att(ncid,varid,'long_name','tile'))
                call check_nc(nf90_put_att(ncid,varid,'units','tile number'))
                call check_nc(nf90_put_att(ncid,varid,'_Storage',"contiguous"))
                call check_nc(nf90_put_att(ncid,varid,'_Endianness',"little"))
            case("pft")         ! Per PFT outputs
                call check_nc(nf90_def_dim(ncid,'pft',icc,pft))
                call check_nc(nf90_def_var(ncid,'pft',nf90_short,pft,varid))
                call check_nc(nf90_put_att(ncid,varid,'long_name','Plant Functional Type'))
                call check_nc(nf90_put_att(ncid,varid,'units','PFT'))
                call check_nc(nf90_put_att(ncid,varid,'_Storage',"contiguous"))
                call check_nc(nf90_put_att(ncid,varid,'_Endianness',"little"))
            case ("layer")      ! Per layer outputs
                call check_nc(nf90_def_dim(ncid,'layer',ignd,layer))
                call check_nc(nf90_def_var(ncid,'layer',nf90_short,layer,varid))
                call check_nc(nf90_put_att(ncid,varid,'long_name','soil layer'))
                call check_nc(nf90_put_att(ncid,varid,'units','layer'))
                call check_nc(nf90_put_att(ncid,varid,'_Storage',"contiguous"))
                call check_nc(nf90_put_att(ncid,varid,'_Endianness',"little"))
        end select

        ! Set up the time dimension
        call check_nc(nf90_def_dim(ncid,'time',nf90_unlimited,time))
        call check_nc(nf90_def_var(ncid,'time',nf90_int,time,varid))
        call check_nc(nf90_put_att(ncid,varid,'long_name','time'))
        call check_nc(nf90_put_att(ncid,varid,'units', trim(timestart)))

        if (leap) then
            call check_nc(nf90_put_att(ncid,varid,'calendar',"standard"))
        else
            call check_nc(nf90_put_att(ncid,varid,'calendar',"365_day"))
        end if

        call check_nc(nf90_put_att(ncid,varid,'_Storage',"chunked"))
        call check_nc(nf90_put_att(ncid,varid,'_Chunksizes',1))
        call check_nc(nf90_put_att(ncid,varid,'_Endianness',"little"))

        !call check_nc(nf90_enddef(ncid))
        !call check_nc(nf90_put_var(ncid,lon,lonvect))
        !call check_nc(nf90_put_var(ncid,lat,latvect))
        !call check_nc(nf90_redef(ncid))

        select case(trim(outputForm))
            case ("tile")       ! Per tile outputs
                call check_nc(nf90_enddef(ncid))
                call check_nc(nf90_put_var(ncid,tile,identityVector(nmos)))
                call check_nc(nf90_redef(ncid))
                call check_nc(nf90_def_var(ncid,trim(descriptor%shortName),nf90_float,[lon,lat,tile,time],varid))
                call check_nc(nf90_put_att(ncid,varid,'_Chunksizes',[cnty,cntx,icc,1]))
            case("pft")         ! Per PFT outputs
                call check_nc(nf90_enddef(ncid))
                call check_nc(nf90_put_var(ncid,pft,identityVector(icc)))
                call check_nc(nf90_redef(ncid))
                if (descriptor%includeBareGround) then
                    ! do something for cells that have bare ground
                    suffix = suffix + 1
                else
                    ! do something for cells that don't have bare ground
                    suffix = suffix
                endif
                call check_nc(nf90_def_var(ncid,trim(descriptor%shortName),nf90_float,[lon,lat,pft,time],varid))
            case ("layer")      ! Per layer outputs
                call check_nc(nf90_enddef(ncid))
                call check_nc(nf90_put_var(ncid,layer,identityVector(ignd)))
                call check_nc(nf90_redef(ncid))
                call check_nc(nf90_def_var(ncid,trim(descriptor%shortName),nf90_float,[lon,lat,layer,time],varid))
            case default        ! Grid average outputs
                call check_nc(nf90_def_var(ncid,trim(descriptor%shortName),nf90_float,[lon,lat,time],varid))
        end select

        call check_nc(nf90_put_att(ncid,varid,'long_name',descriptor%longName))
        call check_nc(nf90_put_att(ncid,varid,'units',descriptor%units))

        !call check_nc(nf90_put_att(ncid,varid,'missing_value',fill_vncdum  alue))
        call check_nc(nf90_put_att(ncid,varid,'_Storage',"chunked"))
        call check_nc(nf90_put_att(ncid,varid,'_DeflateLevel',1))
        !call check_nc(nf90_put_att(ncid,varid,'name_in_code', variables(id)%key))
    end subroutine createNetCDF

    ! Create new netcdf file
    integer function createFile(filename)
        implicit none
        character(*), intent(in)    :: filename
        integer                     :: ncid
        !print*, "Creating filename ", filename
        call check_nc(nf90_create(filename,cmode=NF90_CLOBBER,ncid=ncid))
        createFile = ncid
    end function createFile

    ! Close the file associated with id
    subroutine closeFile(id)
        implicit none
        integer, intent(in) :: id
        integer             :: ncid
        ncid = variables(id)%ncid
        call check_nc(nf90_close(ncid))
    end subroutine closeFile

    subroutine check_nc(status)
        use netcdf
        implicit none
        integer, intent(in) :: status
        if(status /= nf90_noerr) then
            write(0,*)'netCDF error: ',trim(nf90_strerror(status))
            stop
        end if
    end subroutine check_nc

    pure function identityVector(n) result(res)
        integer, allocatable ::res(:)
        integer, intent(in) :: n
        integer             :: i
        allocate(res(n))
        forall (i=1:n)
            res(i) = i
        end forall
    end function identityVector
end module netcdfFileManager
