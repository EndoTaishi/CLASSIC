module netcdf_drivers

    ! This module creates, writes to, and closes netcdf files of the model outputs

    ! J. Melton
    ! Nov 2016

    implicit none

    public :: create_out_netcdf
    public :: setup_nc
    !public :: write_nc
    !public :: netcdf_output
    public :: netcdf_close
    public :: check_nc

    contains

    !--------------------------------------------------------------------------

    subroutine create_out_netcdf()

        use netcdf
        use io_driver, only :  winfo
        use ctem_statevars, only : c_switch

        implicit none

        logical, pointer  :: doperpft
        logical, pointer  :: dopertile
        logical, pointer  :: doday
        logical, pointer  :: dohh
        logical, pointer  :: domonth
        logical, pointer  :: ctem_on
        logical, pointer  :: dofire
        logical, pointer  :: lnduseon
        logical, pointer  :: compete
        logical, pointer  :: domethane

        ! Point pointers
        doperpft        =>c_switch%doperpftoutput
        dopertile       =>c_switch%dopertileoutput
        domonth         =>c_switch%domonthoutput
        doday           =>c_switch%dodayoutput
        dohh            =>c_switch%dohhoutput
        ctem_on         =>c_switch%ctem_on
        dofire          =>c_switch%dofire
        lnduseon        =>c_switch%lnduseon
        compete         =>c_switch%compete
        domethane       =>c_switch%dowetlands

        ! ------------

        ! The doday and dohh are not setup yet so check if true and if so, exit
        if (doday .or. dohh) call xit('dodayoutput and dohhoutput not setup yet', -1)

        ! Give the needed info

        ! Physics variables:
        winfo%shortname = 'rss'
        winfo%standard_name = 'Net Shortwave Surface Radiation'
        winfo%long_name = 'Net downward shortwave radiation at the surface'
        winfo%units = 'W/$m^2$'
        winfo%time_freq = 'annually'
        winfo%nameincode = 'fsstar_yr'
        call setup_nc('pergrd')

        if (domonth) then
            winfo%time_freq = 'monthly'
            winfo%nameincode = 'fsstar_mo'
            call setup_nc('pergrd')
        end if

        winfo%shortname = 'rls'
        winfo%standard_name = 'Net longwave surface radiation'
        winfo%long_name = 'Net longwave surface radiation'
        winfo%units = 'W/$m^2$'
        winfo%time_freq = 'annually'
        winfo%nameincode = 'flstar_yr'
        call setup_nc('pergrd')

        if (domonth) then
            winfo%time_freq = 'monthly'
            winfo%nameincode = 'flstar_mo'
                call setup_nc('pergrd')
        end if

        winfo%shortname = 'hfss'
        winfo%standard_name = 'Surface Net Sensible Heat Flux'
        winfo%long_name = 'Surface Net Sensible Heat Flux'
        winfo%units = 'W/$m^2$'
        winfo%time_freq = 'annually'
        winfo%nameincode = 'qh_yr'
        call setup_nc('pergrd')

        if (domonth) then
            winfo%time_freq = 'monthly'
            winfo%nameincode = 'qh_mo'
            call setup_nc('pergrd')
        end if

        winfo%shortname = 'hfls'
        winfo%standard_name = 'Surface Net Latent Heat Flux'
        winfo%long_name = 'Surface Net Latent Heat Flux'
        winfo%units = 'W/$m^2$'
        winfo%time_freq = 'annually'
        winfo%nameincode = 'qe_yr'
        call setup_nc('pergrd')

        if (domonth) then
            winfo%time_freq = 'monthly'
            winfo%nameincode = 'qe_mo'
            call setup_nc('pergrd')
        end if

        if (domonth) then
            winfo%shortname = 'snw'
            winfo%standard_name = 'Surface Snow Amount'
            winfo%time_freq = 'monthly'
            winfo%long_name ='The mass of surface snow on the land portion of the grid cell'//&
                            ' divided by the land area in the grid cell;'// &
                            'reported as missing where the land fraction is 0;'// &
                            'excludes snow on vegetation canopy or on sea ice.'
            winfo%units = 'kg/$m^2$'
            winfo%nameincode = 'snoacc_mo'
            call setup_nc('pergrd')
        end if

        if (domonth) then
            winfo%shortname = 'wsnw'
            winfo%standard_name = 'Liquid water content of snow pack'
            winfo%time_freq = 'monthly'
            winfo%long_name = 'The total mass of liquid water contained interstitially within the'//&
                            ' whole depth of the snow pack of the land portion of a grid cell'// &
                            ' divided by the area of the land portion of the cell.'
            winfo%units = 'kg/$m^2$'
            winfo%nameincode = 'wsnoacc_mo'
            call setup_nc('pergrd')
        end if

        winfo%shortname = 'mrro'
        winfo%standard_name = 'Total Run-off'
        winfo%long_name = 'The total run-off (including drainage through the base of the soil model)'//&
                        ' per unit area leaving the land portion of the grid cell.'
        winfo%time_freq = 'annually'
        winfo%units = 'mm/year'
        winfo%nameincode = 'rofacc_yr'
        call setup_nc('pergrd')

        if (domonth) then
            winfo%time_freq = 'monthly'
            winfo%units = 'mm/month'
            winfo%nameincode = 'rofacc_mo'
            call setup_nc('pergrd')
        end if

        winfo%shortname = 'pr'
        winfo%standard_name = 'Precipitation'
        winfo%long_name = 'Preciptiation, includes both solid and liquid phases'
        winfo%time_freq = 'annually'
        winfo%units = 'mm/year'
        winfo%nameincode = 'preacc_yr'
        call setup_nc('pergrd')

        if (domonth) then
            winfo%time_freq = 'monthly'
            winfo%units = 'mm/month'
            winfo%nameincode = 'preacc_mo'
            call setup_nc('pergrd')
        end if

        winfo%shortname = 'evspsbl'
        winfo%standard_name = 'Evaporation'
        winfo%long_name = 'Evaporation at surface: flux of water into the atmosphere due to'//&
                        ' conversion of both liquid and solid phases to vapor (from underlying surface and vegetation)'
        winfo%time_freq = 'annually'
        winfo%units = 'mm/year'
        winfo%nameincode = 'evapacc_yr'
        call setup_nc('pergrd')

        if (domonth) then
            winfo%time_freq = 'monthly'
            winfo%units = 'mm/month'
            winfo%nameincode = 'evapacc_mo'
            call setup_nc('pergrd')
        end if

        if (domonth) then
            winfo%shortname = 'evspsblsoi'
            winfo%standard_name = 'Water Evaporation from Soil'
            winfo%time_freq = 'monthly'
            winfo%long_name = 'Water evaporation from soil (including sublimation).'
            winfo%units = 'kg/$m^2$'
            winfo%nameincode = 'groundevap'
            call setup_nc('pergrd')
        end if

        if (domonth) then
            winfo%shortname = 'evspsblveg'
            winfo%standard_name = 'Water Evaporation from Canopy'
            winfo%time_freq = 'monthly'
            winfo%long_name = 'The canopy evaporation and sublimation; may include dew formation as a negative flux.'
            winfo%units = 'mm/month'
            winfo%nameincode = 'canopyevap'
            call setup_nc('pergrd')
        end if

        winfo%shortname = 'tran'
        winfo%standard_name = 'Transpiration'
        winfo%long_name = 'Transpiration'
        winfo%time_freq = 'annually'
        winfo%units = 'mm/year'
        winfo%nameincode = 'transpacc_yr'
        call setup_nc('pergrd')

        if (domonth) then
            winfo%time_freq = 'monthly'
            winfo%units = 'mm/month'
            winfo%nameincode = 'transpacc_mo'
            call setup_nc('pergrd')
        end if

        winfo%shortname = 'ts'
        winfo%standard_name = 'Surface Temperature'
        winfo%time_freq = 'monthly'
        winfo%long_name = 'Temperature of the lower boundary of the atmosphere'
        winfo%units = '$\circ$C'
        winfo%nameincode = 'taacc_mo'
        call setup_nc('pergrd')


        winfo%shortname = 'albs'
        winfo%standard_name = 'Surface Albedo'
        winfo%long_name = 'Grid cell average albedo for all wavelengths.'
        winfo%units = 'fraction'
        winfo%time_freq = 'annually'
        winfo%nameincode = 'altotacc_yr'
        call setup_nc('pergrd')

        if (domonth) then
            winfo%time_freq = 'monthly'
            winfo%nameincode = 'altotacc_mo'
            call setup_nc('pergrd')
        end if

        if (domonth) then
            winfo%shortname = 'tsl'
            winfo%standard_name = 'Temperature of Soil'
            winfo%time_freq = 'monthly'
            winfo%long_name = 'Temperature of each soil layer. Reported as missing for grid cells with no land.'
            winfo%units = '$\circ$C'
            winfo%nameincode = 'TBARACC_MO'
            call setup_nc('perlay')  !per soil layer
        end if

        if (domonth) then
            winfo%shortname = 'mrsll'
            winfo%standard_name = 'Liquid water content of soil layer'
            winfo%time_freq = 'monthly'
            winfo%long_name = 'In each soil layer, the mass of water in liquid phase.'//&
                            ' Reported as "missing" for grid cells occupied entirely by "sea"'
            winfo%units = '$m^3$/$m^3$'
            winfo%nameincode = 'THLQACC_MO'
            call setup_nc('perlay')  !per soil layer
        end if

        if (domonth) then
            winfo%shortname = 'mrsfl'
            winfo%standard_name = 'Frozen water content of soil layer'
            winfo%time_freq = 'monthly'
            winfo%long_name = 'In each soil layer, the mass of water in ice phase.'//&
                            ' Reported as "missing" for grid cells occupied entirely by "sea"'
            winfo%units = '$m^3$/$m^3$'
            winfo%nameincode = 'THICACC_MO'
            call setup_nc('perlay')  !per soil layer
        end if

        if (domonth) then
            winfo%shortname = 'actlyr'
            winfo%standard_name = 'Active Layer Depth'
            winfo%time_freq = 'monthly'
            winfo%long_name = 'Mean active layer depth'
            winfo%units = 'm'
            winfo%nameincode = 'ACTLYR_MO'
            call setup_nc('pergrd')
        end if

        if (domonth) then
            winfo%shortname = 'actlyrmax'
            winfo%standard_name = 'Maximum Active Layer Depth'
            winfo%time_freq = 'monthly'
            winfo%long_name = 'Maximum active layer depth'
            winfo%units = 'm'
            winfo%nameincode = 'ACTLYR_MAX_MO'
            call setup_nc('pergrd')
        end if

        if (domonth) then
            winfo%shortname = 'actlyrmin'
            winfo%standard_name = 'Minimum Active Layer Depth'
            winfo%time_freq = 'monthly'
            winfo%long_name = 'Minimum active layer depth'
            winfo%units = 'm'
            winfo%nameincode = 'ACTLYR_MIN_MO'
            call setup_nc('pergrd')
        end if

        if (domonth) then
            winfo%shortname = 'ftable'
            winfo%standard_name = 'Mean depth to frozen water table'
            winfo%time_freq = 'monthly'
            winfo%long_name = 'Mean depth to frozen water table'
            winfo%units = 'm'
            winfo%nameincode = 'FTABLE_MO'
            call setup_nc('pergrd')
        end if

        if (domonth) then
            winfo%shortname = 'ftablemax'
            winfo%standard_name = 'Maximum depth to frozen water table'
            winfo%time_freq = 'monthly'
            winfo%long_name = 'Maximum depth to frozen water table'
            winfo%units = 'm'
            winfo%nameincode = 'FTABLE_MAX_MO'
            call setup_nc('pergrd')
        end if

        if (domonth) then
            winfo%shortname = 'ftablemin'
            winfo%standard_name = 'Minimum depth to frozen water table'
            winfo%time_freq = 'monthly'
            winfo%long_name = 'Minimum depth to frozen water table'
            winfo%units = 'm'
            winfo%nameincode = 'FTABLE_MIN_MO'
            call setup_nc('pergrd')
        end if

        ! Now the biogeochemical variables ('CTEM')

        if (ctem_on) then

            winfo%shortname = 'lai'
            winfo%standard_name = 'Leaf area index'
            winfo%long_name = 'Leaf Area Index'
            winfo%units = '$m^2$/$m^2$'
            winfo%time_freq = 'annually'
            winfo%nameincode = 'laimaxg_yr_g'
            call setup_nc('pergrd')
            if (domonth) then
                winfo%time_freq = 'monthly'
                winfo%nameincode = 'laimaxg_mo_g'
                call setup_nc('pergrd')
            end if
            if (doperpft) then ! per PFT
                winfo%nameincode = 'laimaxg_yr'
                call setup_nc('perpft')
                if (domonth) then
                    winfo%time_freq = 'monthly'
                    winfo%nameincode = 'laimaxg_mo'
                    call setup_nc('perpft')
                end if
            end if
            if (dopertile) then ! per tile
                winfo%nameincode = 'laimaxg_yr_t'
                call setup_nc('pertil')
                if (domonth) then
                    winfo%time_freq = 'monthly'
                    winfo%nameincode = 'laimaxg_mo_t'
                    call setup_nc('pertil')
                end if
            end if

            winfo%shortname = 'cVeg'
            winfo%standard_name = 'Carbon Mass in Vegetation'
            winfo%time_freq = 'annually'
            winfo%long_name = 'Carbon mass per unit area in vegetation'
            winfo%units = 'kg C/$m^2$'
            winfo%nameincode = 'vgbiomas_yr_g'
            call setup_nc('pergrd')
            if (domonth) then
                winfo%nameincode = 'vgbiomas_mo_g'
                winfo%time_freq = 'monthly'
                call setup_nc('pergrd')
            end if
            if (doperpft) then ! per PFT
                winfo%nameincode = 'vgbiomas_yr'
                call setup_nc('perpft')
                if (domonth) then
                    winfo%nameincode = 'vgbiomas_mo'
                    winfo%time_freq = 'monthly'
                    call setup_nc('perpft')
                end if
            end if
            if (dopertile) then ! per tile
                winfo%nameincode = 'vgbiomas_yr_t'
                call setup_nc('pertil')
                if (domonth) then
                    winfo%nameincode = 'vgbiomas_mo_t'
                    winfo%time_freq = 'monthly'
                    call setup_nc('pertil')
                end if
            end if

            winfo%shortname = 'cStem'
            winfo%standard_name = 'Carbon Mass in Stem'
            winfo%time_freq = 'annually'
            winfo%long_name = 'Carbon Mass in Stem including sapwood and hardwood.'
            winfo%units = 'kg C/$m^2$'
            winfo%nameincode = 'stemmass_yr_g'
            call setup_nc('pergrd')
            if (doperpft) then ! per PFT
                winfo%nameincode = 'stemmass_yr'
                call setup_nc('perpft')
            end if
            if (dopertile) then ! per tile
                winfo%nameincode = 'stemmass_yr_t'
                call setup_nc('pertil')
            end if

            winfo%shortname = 'cRoot'
            winfo%standard_name = 'Carbon Mass in Roots'
            winfo%time_freq = 'annually'
            winfo%long_name = 'Carbon mass per unit area in roots, including fine and coarse roots.'
            winfo%units = 'kg C/$m^2$'
            winfo%nameincode = 'rootmass_yr_g'
            call setup_nc('pergrd')
            if (doperpft) then ! per PFT
                winfo%nameincode = 'rootmass_yr'
                call setup_nc('perpft')
            end if
            if (dopertile) then ! per tile
                winfo%nameincode = 'rootmass_yr_t'
                call setup_nc('pertil')
            end if

            winfo%shortname = 'cLitter'
            winfo%standard_name = 'Carbon Mass in Litter Pool'
            winfo%time_freq = 'annually'
            winfo%long_name = 'Carbon Mass in Litter Pool'
            winfo%units = 'kg C/$m^2$'
            winfo%nameincode = 'litrmass_yr_g'
            call setup_nc('pergrd')
            if (domonth) then
                winfo%nameincode = 'litrmass_mo_t'
                winfo%time_freq = 'monthly'
                call setup_nc('pergrd')
            end if
            if (doperpft) then ! per PFT
                winfo%nameincode = 'litrmass_yr'
                call setup_nc('perpft')
                if (domonth) then
                    winfo%nameincode = 'litrmass_mo'
                    winfo%time_freq = 'monthly'
                    call setup_nc('perpft')
                end if
            end if
            if (dopertile) then ! per tile
                winfo%nameincode = 'litrmass_yr_t'
                call setup_nc('pertil')
                if (domonth) then
                    winfo%nameincode = 'litrmass_mo_t'
                    winfo%time_freq = 'monthly'
                    call setup_nc('pertil')
                end if
            end if

            winfo%shortname = 'cSoil'
            winfo%standard_name = 'Carbon Mass in Model Soil Pool'
            winfo%time_freq = 'annually'
            winfo%long_name = 'Carbon mass in the full depth of the soil model.'
            winfo%units = 'kg C/$m^2$'
            winfo%nameincode = 'soilcmas_yr_g'
            call setup_nc('pergrd')
            if (domonth) then
                winfo%nameincode = 'soilcmas_mo_g'
                winfo%time_freq = 'monthly'
                call setup_nc('pergrd')
            end if
            if (doperpft) then ! per PFT
                winfo%nameincode = 'soilcmas_yr'
                call setup_nc('perpft')
                if (domonth) then
                    winfo%nameincode = 'soilcmas_mo'
                    winfo%time_freq = 'monthly'
                    call setup_nc('perpft')
                end if
            end if
            if (dopertile) then ! per tile
                winfo%nameincode = 'soilcmas_yr_t'
                call setup_nc('pertil')
                if (domonth) then
                    winfo%nameincode = 'soilcmas_mo_t'
                    winfo%time_freq = 'monthly'
                    call setup_nc('pertil')
                end if
            end if

            winfo%shortname = 'cLand'
            winfo%standard_name = 'Total Carbon in All Terrestrial Carbon Pools'
            winfo%time_freq = 'annually'
            winfo%long_name = ' Report missing data over ocean grid cells. For fractional'//&
                            ' land report value averaged over the land fraction.'
            winfo%units = 'kg C/$m^2$'
            winfo%nameincode = 'totcmass_yr_g'
            call setup_nc('pergrd')
            if (doperpft) then ! per PFT
                winfo%nameincode = 'totcmass_yr'
                call setup_nc('perpft')
            end if
            if (dopertile) then ! per tile
                winfo%nameincode = 'totcmass_yr_t'
                call setup_nc('pertil')
            end if

            winfo%shortname = 'vegHeight'
            winfo%standard_name = 'Vegetation Height'
            winfo%time_freq = 'annually'
            winfo%long_name = ' Vegetation height over the vegetated fraction of the gridcell.'
            winfo%units = 'm'
            winfo%nameincode = 'veghght_yr_g'
            call setup_nc('pergrd')
            if (doperpft) then ! per PFT
                winfo%nameincode = 'veghght_yr'
                call setup_nc('perpft')
            end if
            if (dopertile) then ! per tile
                winfo%nameincode = 'veghght_yr_t'
                call setup_nc('pertil')
            end if

            winfo%shortname = 'npp'
            winfo%standard_name = 'Net Primary Productivity of Biomass Expressed As Carbon'
            winfo%time_freq = 'annually'
            winfo%long_name = 'Carbon Mass Flux out of Atmosphere due to Net Primary Production on Land'
            winfo%units = 'kg C/$m^2$'
            winfo%nameincode = 'npp_yr_g'
            call setup_nc('pergrd')
            if (domonth) then
                winfo%nameincode = 'npp_mo_g'
                winfo%time_freq = 'monthly'
                call setup_nc('pergrd')
            end if
            if (doperpft) then ! per PFT
                winfo%nameincode = 'npp_yr'
                call setup_nc('perpft')
                if (domonth) then
                    winfo%nameincode = 'npp_mo'
                    winfo%time_freq = 'monthly'
                    call setup_nc('perpft')
                end if
            end if
            if (dopertile) then ! per tile
                winfo%nameincode = 'npp_yr_t'
                call setup_nc('pertil')
                if (domonth) then
                    winfo%nameincode = 'npp_mo_t'
                    winfo%time_freq = 'monthly'
                    call setup_nc('pertil')
                end if
            end if

            winfo%shortname = 'gpp'
            winfo%standard_name = 'Gross Primary Productivity of Biomass Expressed As Carbon'
            winfo%time_freq = 'annually'
            winfo%long_name = 'Carbon Mass Flux out of Atmosphere due to Gross Primary Production on Land'
            winfo%units = 'kg C/$m^2$'
            winfo%nameincode = 'gpp_yr_g'
            call setup_nc('pergrd')
            if (domonth) then
                winfo%nameincode = 'gpp_mo_g'
                winfo%time_freq = 'monthly'
                call setup_nc('pergrd')
            end if
            if (doperpft) then ! per PFT
                winfo%nameincode = 'gpp_yr'
                call setup_nc('perpft')
                if (domonth) then
                    winfo%nameincode = 'gpp_mo'
                    winfo%time_freq = 'monthly'
                    call setup_nc('perpft')
                end if
            end if
            if (dopertile) then ! per tile
                winfo%nameincode = 'gpp_yr_t'
                call setup_nc('pertil')
                if (domonth) then
                    winfo%nameincode = 'gpp_mo_t'
                    winfo%time_freq = 'monthly'
                    call setup_nc('pertil')
                end if
            end if

            winfo%shortname = 'nep'
            winfo%standard_name = 'Net Carbon Mass Flux out of Atmosphere due to Net Ecosystem Productivity on Land'
            winfo%time_freq = 'annually'
            winfo%long_name = 'Natural flux of CO2 (expressed as a mass flux of carbon) from the atmosphere'//&
                            ' to the land calculated as the difference between uptake associated will photosynthesis'//&
                            ' and the release of CO2 from the sum of plant and soil respiration and fire. Positive'//&
                            ' flux is into the land. emissions from natural fires and human ignition fires as '//&
                            'calculated by the fire module of the dynamic vegetation model, but excluding any'//&
                            ' CO2 flux from fire included in fLuc (CO2 Flux to Atmosphere from Land Use Change).'
            winfo%units = 'kg C/$m^2$'
            winfo%nameincode = 'nep_yr_g'
            call setup_nc('pergrd')
            if (domonth) then
                winfo%nameincode = 'nep_mo_g'
                winfo%time_freq = 'monthly'
                call setup_nc('pergrd')
            end if
            if (doperpft) then ! per PFT
                winfo%nameincode = 'nep_yr'
                call setup_nc('perpft')
                if (domonth) then
                    winfo%nameincode = 'nep_mo'
                    winfo%time_freq = 'monthly'
                    call setup_nc('perpft')
                end if
            end if
            if (dopertile) then ! per tile
                winfo%nameincode = 'nep_yr_t'
                call setup_nc('pertil')
                if (domonth) then
                    winfo%nameincode = 'nep_mo_t'
                    winfo%time_freq = 'monthly'
                    call setup_nc('pertil')
                end if
            end if

            winfo%shortname = 'nbp'
            winfo%standard_name = 'Net Carbon Mass Flux out of Atmosphere due to Net Biospheric Productivity on Land'
            winfo%time_freq = 'annually'
            winfo%long_name = 'This is the net mass flux of carbon from atmosphere into land, calculated as photosynthesis'//&
                            ' MINUS the sum of plant and soil respiration, carbon fluxes from fire, harvest, grazing'//&
                            ' and land use change. Positive flux is into the land.'
            winfo%units = 'kg C/$m^2$'
            winfo%nameincode = 'nbp_yr_g'
            call setup_nc('pergrd')
            if (domonth) then
                winfo%nameincode = 'nbp_mo_g'
                winfo%time_freq = 'monthly'
                call setup_nc('pergrd')
            end if
            if (doperpft) then ! per PFT
                winfo%nameincode = 'nbp_yr'
                call setup_nc('perpft')
                if (domonth) then
                    winfo%nameincode = 'nbp_mo'
                    winfo%time_freq = 'monthly'
                    call setup_nc('perpft')
                end if
            end if
            if (dopertile) then ! per tile
                winfo%nameincode = 'nbp_yr_t'
                call setup_nc('pertil')
                if (domonth) then
                    winfo%nameincode = 'nbp_mo_t'
                    winfo%time_freq = 'monthly'
                    call setup_nc('pertil')
                end if
            end if

            winfo%shortname = 'rh'
            winfo%standard_name = 'Carbon Mass Flux into Atmosphere due to Heterotrophic Respiration on Land'
            winfo%time_freq = 'annually'
            winfo%long_name = 'Carbon mass flux per unit area into atmosphere due to heterotrophic'//&
                            ' respiration on land (respiration by consumers)'
            winfo%units = 'kg C/$m^2$'
            winfo%nameincode = 'hetrores_yr_g'
            call setup_nc('pergrd')
            if (domonth) then
                winfo%nameincode = 'hetrores_mo_g'
                winfo%time_freq = 'monthly'
                call setup_nc('pergrd')
            end if
            if (doperpft) then ! per PFT
                winfo%nameincode = 'hetrores_yr'
                call setup_nc('perpft')
                if (domonth) then
                    winfo%nameincode = 'hetrores_mo'
                    winfo%time_freq = 'monthly'
                    call setup_nc('perpft')
                end if
            end if
            if (dopertile) then ! per tile
                winfo%nameincode = 'hetrores_yr_t'
                call setup_nc('pertil')
                if (domonth) then
                    winfo%nameincode = 'hetrores_mo_t'
                    winfo%time_freq = 'monthly'
                    call setup_nc('pertil')
                end if
            end if

            winfo%shortname = 'ra'
            winfo%standard_name = 'Total autotrophic respiration'
            winfo%time_freq = 'annually'
            winfo%long_name = 'Autotrophic Respiration from All Vegetation Pools'
            winfo%units = 'kg C/$m^2$'
            winfo%nameincode = 'autores_yr_g'
            call setup_nc('pergrd')
            if (domonth) then
                winfo%nameincode = 'autores_mo_g'
                winfo%time_freq = 'monthly'
                call setup_nc('pergrd')
            end if
            if (doperpft) then ! per PFT
                winfo%nameincode = 'autores_yr'
                call setup_nc('perpft')
                if (domonth) then
                    winfo%nameincode = 'autores_mo'
                    winfo%time_freq = 'monthly'
                    call setup_nc('perpft')
                end if
            end if
            if (dopertile) then ! per tile
                winfo%nameincode = 'autores_yr_t'
                call setup_nc('pertil')
                if (domonth) then
                    winfo%nameincode = 'autores_mo_t'
                    winfo%time_freq = 'monthly'
                    call setup_nc('pertil')
                end if
            end if

            winfo%shortname = 'rhLitter'
            winfo%standard_name = 'Litter heterotrophic respiration'
            winfo%time_freq = 'annually'
            winfo%long_name = 'Carbon Mass Flux into Atmosphere due to Heterotrophic Respiration from Litter on Land'
            winfo%units = 'kg C/$m^2$'
            winfo%nameincode = 'litres_yr_g'
            call setup_nc('pergrd')
            if (domonth) then
                winfo%nameincode = 'litres_mo_g'
                winfo%time_freq = 'monthly'
                call setup_nc('pergrd')
            end if
            if (doperpft) then ! per PFT
                winfo%nameincode = 'litres_yr'
                call setup_nc('perpft')
                if (domonth) then
                    winfo%nameincode = 'litres_mo'
                    winfo%time_freq = 'monthly'
                    call setup_nc('perpft')
                end if
            end if
            if (dopertile) then ! per tile
                winfo%nameincode = 'litres_yr_t'
                call setup_nc('pertil')
                if (domonth) then
                    winfo%nameincode = 'litres_mo_t'
                    winfo%time_freq = 'monthly'
                    call setup_nc('pertil')
                end if
            end if

            winfo%shortname = 'rhSoil'
            winfo%standard_name = 'Soil heterotrophic respiration'
            winfo%time_freq = 'annually'
            winfo%long_name = 'Carbon Mass Flux into Atmosphere due to Heterotrophic Respiration from Soil on Land'
            winfo%units = 'kg C/$m^2$'
            winfo%nameincode = 'soilcres_yr_g'
            call setup_nc('pergrd')
            if (domonth) then
                winfo%nameincode = 'soilcres_mo_g'
                winfo%time_freq = 'monthly'
                call setup_nc('pergrd')
            end if
            if (doperpft) then ! per PFT
                winfo%nameincode = 'soilcres_yr'
                call setup_nc('perpft')
                if (domonth) then
                    winfo%nameincode = 'soilcres_mo'
                    winfo%time_freq = 'monthly'
                    call setup_nc('perpft')
                end if
            end if
            if (dopertile) then ! per tile
                winfo%nameincode = 'soilcres_yr_t'
                call setup_nc('pertil')
                if (domonth) then
                    winfo%nameincode = 'soilcres_mo_t'
                    winfo%time_freq = 'monthly'
                    call setup_nc('pertil')
                end if
            end if

            if (domonth) then
                winfo%shortname = 'fVegLitter'
                winfo%standard_name = 'Total Carbon Mass Flux from Vegetation to Litter'
                winfo%time_freq = 'monthly'
                winfo%long_name = 'Total Carbon Mass Flux from Vegetation to Litter'
                winfo%units = 'kg C/$m^2$'
                winfo%nameincode = 'litrfallveg_mo_g'
                call setup_nc('pergrd')
                if (doperpft) then ! per PFT
                    winfo%nameincode = 'litrfallveg_mo'
                    call setup_nc('perpft')
                end if
                if (dopertile) then ! per tile
                    winfo%nameincode = 'litrfallveg_mo_t'
                    call setup_nc('pertil')
                end if
            end if

            if (domonth) then
                winfo%shortname = 'fLitterSoil'
                winfo%standard_name = 'Total Carbon Mass Flux from Litter to Soil'
                winfo%time_freq = 'monthly'
                winfo%long_name = 'Total Carbon Mass Flux from Litter to Soil'
                winfo%units = 'kg C/$m^2$'
                winfo%nameincode = 'humiftrsveg_mo_g'
                call setup_nc('pergrd')
                if (doperpft) then ! per PFT
                    winfo%nameincode = 'humiftrsveg_mo'
                    call setup_nc('perpft')
                end if
                if (dopertile) then ! per tile
                    winfo%nameincode = 'humiftrsveg_mo_t'
                    call setup_nc('pertil')
                end if
            end if

            if (dofire) then

                winfo%shortname = 'fFire'
                winfo%standard_name = 'Carbon Mass Flux into Atmosphere due to CO2 Emission from Fire'
                winfo%time_freq = 'annually'
                winfo%long_name = 'CO2 emissions (expressed as a CO2 mass flux per unit area) from'//&
                                'natural fires and human ignition fires as calculated by the fire module'//&
                                ' of the dynamic vegetation model, but excluding any CO2 flux from fire'//&
                                ' included in fLuc (CO2 Flux to Atmosphere from Land Use Change)'
                winfo%units = 'g CO$_2$ $m^{-2}$ year$^{-1}$'
                winfo%nameincode = 'emit_co2_yr_g'
                call setup_nc('pergrd')
                if (domonth) then
                    winfo%nameincode = 'emit_co2_mo_g'
                    winfo%units = 'g CO$_2$ $m^{-2}$ month$^{-1}$'
                    winfo%time_freq = 'monthly'
                    call setup_nc('pergrd')
                end if
                if (doperpft) then ! per PFT
                    winfo%nameincode = 'emit_co2_yr'
                    call setup_nc('perpft')
                    if (domonth) then
                        winfo%nameincode = 'emit_co2_mo'
                        winfo%units = 'g CO$_2$ $m^{-2}$ month$^{-1}$'
                        winfo%time_freq = 'monthly'
                        call setup_nc('perpft')
                    end if
                end if
                if (dopertile) then ! per tile
                    winfo%nameincode = 'emit_co2_yr_t'
                    call setup_nc('pertil')
                    if (domonth) then
                        winfo%nameincode = 'emit_co2_mo_t'
                        winfo%units = 'g CO$_2$ $m^{-2}$ month$^{-1}$'
                        winfo%time_freq = 'monthly'
                        call setup_nc('pertil')
                    end if
                end if

                winfo%shortname = 'fFireCH4'
                winfo%standard_name = 'Carbon Mass Flux into Atmosphere due to CH4 Emission from Fire'
                winfo%time_freq = 'annually'
                winfo%long_name = 'CH4 emissions (expressed as a CH4 mass flux per unit area) from'//&
                                'natural fires and human ignition fires as calculated by the fire module'//&
                                ' of the dynamic vegetation model'
                winfo%units = 'g CH$_4$ $m^{-2}$ year$^{-1}$'
                winfo%nameincode = 'emit_ch4_yr_g'
                call setup_nc('pergrd')
                if (domonth) then
                    winfo%nameincode = 'emit_ch4_mo_g'
                    winfo%units = 'g CH$_4$ $m^{-2}$ month$^{-1}$'
                    winfo%time_freq = 'monthly'
                    call setup_nc('pergrd')
                end if
                if (doperpft) then ! per PFT
                    winfo%nameincode = 'emit_ch4_yr'
                    call setup_nc('perpft')
                    if (domonth) then
                        winfo%nameincode = 'emit_ch4_mo'
                        winfo%units = 'g CH$_4$ $m^{-2}$ month$^{-1}$'
                        winfo%time_freq = 'monthly'
                        call setup_nc('perpft')
                    end if
                end if
                if (dopertile) then ! per tile
                    winfo%nameincode = 'emit_ch4_yr_t'
                    call setup_nc('pertil')
                    if (domonth) then
                        winfo%nameincode = 'emit_ch4_mo_t'
                        winfo%units = 'g CH$_4$ $m^{-2}$ month$^{-1}$'
                        winfo%time_freq = 'monthly'
                        call setup_nc('pertil')
                    end if
                end if

                winfo%shortname = 'burntFractionAll'
                winfo%standard_name = 'Percentage of Entire Grid cell that is Covered by Burnt Vegetation (All Classes)'
                winfo%time_freq = 'annually'
                winfo%long_name = 'Percentage of grid cell burned due to all fires including natural'//&
                                ' and anthropogenic fires'
                winfo%units = 'fraction'
                winfo%nameincode = 'burnfrac_yr_g'
                call setup_nc('pergrd')
                if (domonth) then
                    winfo%nameincode = 'burnfrac_mo_g'
                    winfo%time_freq = 'monthly'
                    call setup_nc('pergrd')
                end if
                if (doperpft) then ! per PFT
                    winfo%nameincode = 'burnfrac_yr'
                    call setup_nc('perpft')
                    if (domonth) then
                        winfo%nameincode = 'burnfrac_mo'
                        winfo%time_freq = 'monthly'
                        call setup_nc('perpft')
                    end if
                end if
                if (dopertile) then ! per tile
                    winfo%nameincode = 'burnfrac_yr_t'
                    call setup_nc('pertil')
                    if (domonth) then
                        winfo%nameincode = 'burnfrac_mo_t'
                        winfo%time_freq = 'monthly'
                        call setup_nc('pertil')
                    end if
                end if

            end if !dofire

            if (lnduseon) then

                winfo%shortname = 'fDeforestToAtmos'
                winfo%standard_name = 'Deforested biomass that goes into atmosphere as a result of anthropogenic land use change'
                winfo%time_freq = 'annually'
                winfo%long_name = 'When land use change results in deforestation of natural vegetation '//&
                                '(trees or grasslands) then natural biomass is removed.'
                winfo%units = 'g C $m^{-2}$ year$^{-1}$'
                winfo%nameincode = 'luc_emc_yr_g'
                call setup_nc('pergrd')
                if (domonth) then
                    winfo%nameincode = 'luc_emc_mo_g'
                    winfo%units = 'g C $m^{-2}$ month$^{-1}$'
                    winfo%time_freq = 'monthly'
                    call setup_nc('pergrd')
                end if
                if (doperpft) then ! per PFT
                    winfo%nameincode = 'luc_emc_yr'
                    call setup_nc('perpft')
                    if (domonth) then
                        winfo%nameincode = 'luc_emc_mo'
                        winfo%units = 'g C $m^{-2}$ month$^{-1}$'
                        winfo%time_freq = 'monthly'
                        call setup_nc('perpft')
                    end if
                end if
                if (dopertile) then ! per tile
                    winfo%nameincode = 'luc_emc_yr_t'
                    call setup_nc('pertil')
                    if (domonth) then
                        winfo%nameincode = 'luc_emc_mo_t'
                        winfo%units = 'g C $m^{-2}$ month$^{-1}$'
                        winfo%time_freq = 'monthly'
                        call setup_nc('pertil')
                    end if
                end if

                winfo%shortname = 'fDeforestToLitter'
                winfo%standard_name = 'Deforested biomass that goes into litter as a result of anthropogenic land use change'
                winfo%time_freq = 'annually'
                winfo%long_name = 'When land use change results in deforestation of natural vegetation '//&
                                '(trees or grasslands) then natural biomass is removed and some transferred to litter.'
                winfo%units = 'g C $m^{-2}$ year$^{-1}$'
                winfo%nameincode = 'lucltrin_yr_g'
                call setup_nc('pergrd')
                if (domonth) then
                    winfo%nameincode = 'lucltrin_mo_g'
                    winfo%units = 'g C $m^{-2}$ month$^{-1}$'
                    winfo%time_freq = 'monthly'
                    call setup_nc('pergrd')
                end if
                if (doperpft) then ! per PFT
                    winfo%nameincode = 'lucltrin_yr'
                    call setup_nc('perpft')
                    if (domonth) then
                        winfo%nameincode = 'lucltrin_mo'
                        winfo%units = 'g C $m^{-2}$ month$^{-1}$'
                        winfo%time_freq = 'monthly'
                        call setup_nc('perpft')
                    end if
                end if
                if (dopertile) then ! per tile
                    winfo%nameincode = 'lucltrin_yr_t'
                    call setup_nc('pertil')
                    if (domonth) then
                        winfo%nameincode = 'lucltrin_mo_t'
                        winfo%units = 'g C $m^{-2}$ month$^{-1}$'
                        winfo%time_freq = 'monthly'
                        call setup_nc('pertil')
                    end if
                end if

                winfo%shortname = 'fDeforestToSoil'
                winfo%standard_name = 'Deforested biomass that goes into soil as a result of anthropogenic land use change'
                winfo%time_freq = 'annually'
                winfo%long_name = 'When land use change results in deforestation of natural vegetation '//&
                                '(trees or grasslands) then natural biomass is removed and some transferred to soil.'
                winfo%units = 'g C $m^{-2}$ year$^{-1}$'
                winfo%nameincode = 'lucsocin_yr_g'
                call setup_nc('pergrd')
                if (domonth) then
                    winfo%nameincode = 'lucsocin_mo_g'
                    winfo%units = 'g C $m^{-2}$ month$^{-1}$'
                    winfo%time_freq = 'monthly'
                    call setup_nc('pergrd')
                end if
                if (doperpft) then ! per PFT
                    winfo%nameincode = 'lucsocin_yr'
                    call setup_nc('perpft')
                    if (domonth) then
                        winfo%nameincode = 'lucsocin_mo'
                        winfo%units = 'g C $m^{-2}$ month$^{-1}$'
                        winfo%time_freq = 'monthly'
                        call setup_nc('perpft')
                    end if
                end if
                if (dopertile) then ! per tile
                    winfo%nameincode = 'lucsocin_yr_t'
                    call setup_nc('pertil')
                    if (domonth) then
                        winfo%nameincode = 'lucsocin_mo_t'
                        winfo%units = 'g C $m^{-2}$ month$^{-1}$'
                        winfo%time_freq = 'monthly'
                        call setup_nc('pertil')
                    end if
                end if

            end if !lnduseon

            if (lnduseon .or. compete) then

                winfo%shortname = 'landCoverFrac'
                winfo%standard_name = 'Percentage of Area by Vegetation/Land Cover Category'
                winfo%time_freq = 'annually'
                winfo%long_name = ' Percentage of grid cell area occupied by different model vegetation/land cover categories.'
                winfo%units = '%'
                winfo%nameincode = 'fcancmxrow'
                call setup_nc('pergrd')
                if (domonth) then
                    winfo%time_freq = 'monthly'
                    call setup_nc('pergrd')
                end if

                winfo%shortname = 'landCoverExist'
                winfo%standard_name = 'Boolean for existance of Vegetation/Land Cover Category based on bioclimatic indices'
                winfo%time_freq = 'annually'
                winfo%long_name = 'Boolean for existance of Vegetation/Land Cover Category based on bioclimatic indices'
                winfo%units = 'T/F'
                winfo%nameincode = 'pftexistrow'
                call setup_nc('pergrd')
                if (domonth) then
                    winfo%time_freq = 'monthly'
                    call setup_nc('pergrd')
                end if

            end if !lnduseon or compete

            if (domethane) then

                winfo%shortname = 'wetlandCH4spec'
                winfo%standard_name = 'Methane emissions from wetlands'
                winfo%time_freq = 'annually'
                winfo%long_name = 'Methane emissions from specified wetland areas based on heterotrophic respiration'
                winfo%units = 'g CH$_4$ $m^{-2}$ year$^{-1}$'
                winfo%nameincode = 'ch4wet1_yr_g'
                call setup_nc('pergrd')
                if (domonth) then
                    winfo%nameincode = 'ch4wet1_mo_g'
                    winfo%units = 'g CH$_4$ $m^{-2}$ month$^{-1}$'
                    winfo%time_freq = 'monthly'
                    call setup_nc('pergrd')
                end if
                if (dopertile) then ! per tile
                    winfo%nameincode = 'ch4wet1_yr_t'
                    call setup_nc('pertil')
                    if (domonth) then
                        winfo%nameincode = 'ch4wet1_mo_t'
                        winfo%units = 'g CH$_4$ $m^{-2}$ month$^{-1}$'
                        winfo%time_freq = 'monthly'
                        call setup_nc('pertil')
                    end if
                end if

                winfo%shortname = 'wetlandCH4dyn'
                winfo%standard_name = 'Methane emissions from wetlands'
                winfo%time_freq = 'annually'
                winfo%long_name = 'Methane emissions from dynamically determined wetland areas based on heterotrophic respiration'
                winfo%units = 'g CH$_4$ $m^{-2}$ year$^{-1}$'
                winfo%nameincode = 'ch4dyn1_yr_g'
                call setup_nc('pergrd')
                if (domonth) then
                    winfo%nameincode = 'ch4dyn1_mo_g'
                    winfo%units = 'g CH$_4$ $m^{-2}$ month$^{-1}$'
                    winfo%time_freq = 'monthly'
                    call setup_nc('pergrd')
                end if
                if (dopertile) then ! per tile
                    winfo%nameincode = 'ch4dyn1_yr_t'
                    call setup_nc('pertil')
                    if (domonth) then
                        winfo%nameincode = 'ch4dyn1_mo_t'
                        winfo%units = 'g CH$_4$ $m^{-2}$ month$^{-1}$'
                        winfo%time_freq = 'monthly'
                        call setup_nc('pertil')
                    end if
                end if

                winfo%shortname = 'soilCH4cons'
                winfo%standard_name = 'Methane uptake (methanotrophy) by oxic soils'
                winfo%time_freq = 'annually'
                winfo%long_name = 'Methane uptake (methanotrophy) by oxic soils'
                winfo%units = 'g CH$_4$ $m^{-2}$ year$^{-1}$'
                winfo%nameincode = 'ch4soills_yr_g'
                call setup_nc('pergrd')
                if (domonth) then
                    winfo%nameincode = 'ch4soills_mo_g'
                    winfo%units = 'g CH$_4$ $m^{-2}$ month$^{-1}$'
                    winfo%time_freq = 'monthly'
                    call setup_nc('pergrd')
                end if
                if (dopertile) then ! per tile
                    winfo%nameincode = 'ch4soills_yr_t'
                    call setup_nc('pertil')
                    if (domonth) then
                        winfo%nameincode = 'ch4soills_mo_t'
                        winfo%units = 'g CH$_4$ $m^{-2}$ month$^{-1}$'
                        winfo%time_freq = 'monthly'
                        call setup_nc('pertil')
                    end if
                end if

                winfo%shortname = 'wetlandFrac'
                winfo%standard_name = 'Fraction of a grid cell covered by wetland'
                winfo%time_freq = 'annually'
                winfo%long_name = 'Fraction of a grid cell covered by wetland. Prognostically determined.'
                winfo%units = 'Fraction'
                winfo%nameincode = 'wetfdyn_yr_g'
                call setup_nc('pergrd')
                if (domonth) then
                    winfo%nameincode = 'wetfdyn_mo_g'
                    winfo%time_freq = 'monthly'
                    call setup_nc('pergrd')
                end if
                if (dopertile) then ! per tile
                    winfo%nameincode = 'wetfdyn_yr_t'
                    call setup_nc('pertil')
                    if (domonth) then
                        winfo%nameincode = 'wetfdyn_mo_t'
                        winfo%time_freq = 'monthly'
                        call setup_nc('pertil')
                    end if
                end if

            end if ! domethane

        end if !ctem_on


    end subroutine create_out_netcdf

    ! !=======================================================================

    subroutine setup_nc(kindout)

        use netcdf
        use io_driver, only :  infopak,winfo,lonvect,latvect,srtx,srty,cntx,cnty,&
                         yrange,xrange, fill_value, timestart
        use ctem_statevars,     only : c_switch
        use ctem_params,        only : ignd,icc,nmos

        implicit none

        character(6), intent(in) :: kindout    !either perpft, pertil, perlay, or pergrd

        type(infopak) :: varinfo

        character(8)  :: today
        character(10) :: now
        character(80) :: filename
        integer :: ncid, lon, lat, varid, tile, pft, layer, time, i

        integer, allocatable, dimension(:) :: layersnum
        integer, allocatable, dimension(:) :: pftnum
        integer, allocatable, dimension(:) :: tilesnum

        character(:), pointer :: Comment   !< Comment about the run that will be written to the output netcdfs
        logical, pointer :: leap           !< set to true if all/some leap years in the .MET file have data for 366 days
                                        !< also accounts for leap years in .MET when cycling over meteorology (cyclemet)
        character(:), pointer :: output_directory !< Directory where the output netcdfs will be placed

        ! Point pointers
        leap => c_switch%leap
        Comment => c_switch%Comment
        output_directory => c_switch%output_directory

        ! ---------

        varinfo = winfo
        varinfo%comment = Comment ! Use the comment from the joboptions file.

        ! Make the output file name
        if (kindout == 'pertil') then
            filename = trim(output_directory)//'/'// trim(varinfo%shortname) // '_' // trim(varinfo%time_freq) // '_pertil.nc'
        else if (kindout == 'perpft') then
            filename = trim(output_directory)//'/'// trim(varinfo%shortname) // '_' // trim(varinfo%time_freq) // '_perpft.nc'
        else
            filename = trim(output_directory)//'/'// trim(varinfo%shortname) // '_' // trim(varinfo%time_freq) // '.nc'
        end if


        call check_nc(nf90_create(filename,cmode=nf90_netcdf4,ncid=ncid))
        call check_nc(nf90_put_att(ncid,nf90_global,'title','CLASSIC output file'))

        call date_and_time(today,now)

        call check_nc(nf90_put_att(ncid,nf90_global,'timestamp',today//' '//now(1:4)))
        call check_nc(nf90_put_att(ncid,nf90_global,'Conventions','COARDS'))
        call check_nc(nf90_put_att(ncid,nf90_global,'node_offset',1))

        !----1 - Longitude
        call check_nc(nf90_def_dim(ncid,'lon',cntx,lon))
        call check_nc(nf90_def_var(ncid,'lon',nf90_float,lon,varid))
        call check_nc(nf90_put_att(ncid,varid,'long_name','longitude'))
        call check_nc(nf90_put_att(ncid,varid,'units','degrees_east'))
        call check_nc(nf90_put_att(ncid,varid,'actual_range',xrange))
        call check_nc(nf90_put_att(ncid,varid,'_Storage',"contiguous"))

        !----2 - Latitude
        call check_nc(nf90_def_dim(ncid,'lat',cnty,lat))
        call check_nc(nf90_def_var(ncid,'lat',nf90_float,lat,varid))
        call check_nc(nf90_put_att(ncid,varid,'long_name','latitude'))
        call check_nc(nf90_put_att(ncid,varid,'units','degrees_north'))
        call check_nc(nf90_put_att(ncid,varid,'actual_range',yrange))
        call check_nc(nf90_put_att(ncid,varid,'_Storage',"contiguous"))

        if (kindout == 'pertil') then   ! Per tile outputs

            call check_nc(nf90_def_dim(ncid,'tile',nmos,tile))
            call check_nc(nf90_def_var(ncid,'tile',nf90_short,tile,varid))
            call check_nc(nf90_put_att(ncid,varid,'long_name','tile'))
            call check_nc(nf90_put_att(ncid,varid,'units','tile number'))
            call check_nc(nf90_put_att(ncid,varid,'_Storage',"contiguous"))
            call check_nc(nf90_put_att(ncid,varid,'_Endianness',"little"))

        else if (kindout == 'perpft') then ! Per PFT outputs

            call check_nc(nf90_def_dim(ncid,'pft',icc,pft))
            call check_nc(nf90_def_var(ncid,'pft',nf90_short,pft,varid))
            call check_nc(nf90_put_att(ncid,varid,'long_name','Plant Functional Type'))
            call check_nc(nf90_put_att(ncid,varid,'units','PFT'))
            call check_nc(nf90_put_att(ncid,varid,'_Storage',"contiguous"))
            call check_nc(nf90_put_att(ncid,varid,'_Endianness',"little"))

        else if (kindout == 'perlay') then ! Per soil layer outputs

            call check_nc(nf90_def_dim(ncid,'layer',ignd,layer))
            call check_nc(nf90_def_var(ncid,'layer',nf90_short,layer,varid))
            call check_nc(nf90_put_att(ncid,varid,'long_name','soil layer'))
            call check_nc(nf90_put_att(ncid,varid,'units','layer'))
            call check_nc(nf90_put_att(ncid,varid,'_Storage',"contiguous"))
            call check_nc(nf90_put_att(ncid,varid,'_Endianness',"little"))

        end if

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

        call check_nc(nf90_enddef(ncid))

        ! Fill in the dimension variables and define the model output vars

        call check_nc(nf90_put_var(ncid,lon,lonvect))
        call check_nc(nf90_put_var(ncid,lat,latvect))

        if (kindout == 'pertil') then   ! Per tile outputs

            allocate(tilesnum(nmos))
            forall (i=1:nmos)
                tilesnum(i) = i
            end forall
            call check_nc(nf90_put_var(ncid,tile,tilesnum))
            call check_nc(nf90_def_var(ncid,trim(varinfo%shortname),nf90_float,[lon,lat,tile,time],varid))
            call check_nc(nf90_put_att(ncid,varid,'_Chunksizes',[cnty,cntx,icc,1]))
            deallocate(tilesnum)

        else if (kindout == 'perpft') then ! Per PFT outputs

            allocate(pftnum(icc))
            forall (i=1:icc)
                pftnum(i) = i
            end forall
            call check_nc(nf90_put_var(ncid,pft,pftnum))
            call check_nc(nf90_def_var(ncid,trim(varinfo%shortname),nf90_float,[lon,lat,pft,time],varid))
            deallocate(pftnum)

        else if (kindout == 'perlay') then ! Per soil layer outputs

            allocate(layersnum(ignd))
            forall (i=1:ignd)
                layersnum(i) = i
            end forall
            call check_nc(nf90_put_var(ncid,layer,layersnum))
            call check_nc(nf90_def_var(ncid,trim(varinfo%shortname),nf90_float,[lon,lat,layer,time],varid))
            deallocate(layersnum)

        else  ! Grid average outputs

            call check_nc(nf90_def_var(ncid,trim(varinfo%shortname),nf90_float,[lon,lat,time],varid))

        end if

        call check_nc(nf90_put_att(ncid,varid,'long_name',varinfo%long_name))
        call check_nc(nf90_put_att(ncid,varid,'units',varinfo%units))

        !call check_nc(nf90_put_att(ncid,varid,'_FillValue',fill_value)) !FLAG this is not working at present.

        call check_nc(nf90_put_att(ncid,varid,'missing_value',fill_value))
        call check_nc(nf90_put_att(ncid,varid,'_Storage',"chunked"))
        call check_nc(nf90_put_att(ncid,varid,'_DeflateLevel',1))
        call check_nc(nf90_put_att(ncid,varid,'name_in_code',varinfo%nameincode))
        call check_nc(nf90_put_att(ncid,varid,'Comment',varinfo%Comment))
        call check_nc(nf90_enddef(ncid))

        !close the netcdf
        call check_nc(nf90_close(ncid))

    end subroutine setup_nc

!  !=============================================================================
!
! subroutine write_nc(towrite,kindout)
!
!     use netcdf
!     use outinfo, only :  infopak,winfo,lonvect,latvect,srtx,srty,cntx,cnty,&
!                          yrange,xrange,nmos,icc, fill_value
!     use ctem_params, only : ignd,icc,
!
!     implicit none
!
!     character(6), intent(in) :: kindout    !either perpft, pertil, perlay, or pergrd
!     real, dimension(:,:,:), intent(in) :: towrite
!
!     type(infopak) :: varinfo
!     integer :: ncid, varid, tile, pft, layer, time, i
!     character(80) :: filename
!     ! ------
!
!     varinfo = winfo
!
!     if (kindout == 'pertil') then
!         filename = trim(varinfo%shortname) // '_' // trim(varinfo%time_freq) // '_pertil.nc'
!     else if (kindout == 'perpft') then
!         filename = trim(varinfo%shortname) // '_' // trim(varinfo%time_freq) // '_perpft.nc'
!     else
!         filename = trim(varinfo%shortname) // '_' // trim(varinfo%time_freq) // '.nc'
!     end if
!
!     call check_nc(nf90_open(filename,nf90_write,ncid))
!     call check_nc(nf90_inq_varid(ncid,trim(varinfo%shortname), varid))
!     call check_nc(nf90_put_var(ncid,varid,towrite))!,start=[1,1,yrst],count=[1,1,totyrs]))
        !call check_nc(nf90_close(ncid))
!
!
!
! end subroutine write_nc

    !=============================================================================

    subroutine netcdf_close(ncid)

        use netcdf
        use netcdf_error

        implicit none

        integer, intent(in) :: ncid

        status = nf90_close(ncid)
        if (status/=nf90_noerr) call handle_err(status)

    end subroutine netcdf_close

    !=======================================================================


    subroutine check_nc(status)

        use netcdf
        use typesizes

        implicit none

        !Internal subroutine - checks error status after each netcdf call,
        !prints out text message each time an error code is returned.

        integer, intent(in) :: status

        if(status /= nf90_noerr) then
            write(0,*)'netCDF error: ',trim(nf90_strerror(status))
            stop
        end if

    end subroutine check_nc

    !=======================================================================

end module netcdf_drivers
