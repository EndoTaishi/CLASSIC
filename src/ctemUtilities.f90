!> Central module for CTEM (biogeochemical)-related utilities
module ctemUtilities

implicit none

public :: dayEndCTEMPreparation
public :: accumulateForCTEM
public :: ctemInit

contains

! --------------------------------------------------------------------------------------------------------------------
!> Prepare the CTEM input (physics) variables at the end of the day.
subroutine dayEndCTEMPreparation(nml,nday)

    use ctem_params, only : icc,ignd
    use ctem_statevars, only : vgat,ctem_tile

    implicit none

    integer, intent(in) :: nml     !<Counter representing number of mosaic tiles on modelled domain that are land
    integer, intent(in) :: nday    !<Number of short (physics) timesteps in one day. e.g., if physics timestep is 15 min this is 48.

    integer, pointer, dimension(:) :: altotcount_ctm !nlat  !< Counter used for calculating total albedo
    real, pointer, dimension(:)    :: fsinacc_gat !(ilg)    !<
    real, pointer, dimension(:)    :: flutacc_gat !(ilg)    !<
    real, pointer, dimension(:)    :: flinacc_gat !(ilg)    !<
    real, pointer, dimension(:)    :: altotacc_gat !(ilg)   !<
    !real, pointer, dimension(:)    :: pregacc_gat !(ilg)    !<
    real, pointer, dimension(:)    :: netrad_gat !(ilg)     !<
!     real, pointer, dimension(:)    :: preacc_gat !(ilg)     !<
    real, pointer, dimension(:) :: anmossac_t
    real, pointer, dimension(:) :: rmlmossac_t
    real, pointer, dimension(:) :: gppmossac_t
    real, pointer, dimension(:) :: fsnowacc_t       !<
    real, pointer, dimension(:) :: tcansacc_t       !<
    real, pointer, dimension(:) :: tcanoaccgat_t    !<
    real, pointer, dimension(:) :: taaccgat_t       !<
    real, pointer, dimension(:) :: uvaccgat_t       !<
    real, pointer, dimension(:) :: vvaccgat_t       !<
    real, pointer, dimension(:,:) :: tbaraccgat_t!<
    real, pointer, dimension(:,:) :: tbarcacc_t  !<
    real, pointer, dimension(:,:) :: tbarcsacc_t !<
    real, pointer, dimension(:,:) :: tbargacc_t  !<
    real, pointer, dimension(:,:) :: tbargsacc_t !<
    real, pointer, dimension(:,:) :: thliqcacc_t !<
    real, pointer, dimension(:,:) :: thliqgacc_t !<
    real, pointer, dimension(:,:) :: thliqacc_t  !<
    real, pointer, dimension(:,:) :: thicecacc_t !<
    real, pointer, dimension(:,:) :: thicegacc_t !<
    real, pointer, dimension(:,:) :: thiceacc_t  !< Added in place of YW's thicaccgat_m. EC Dec 23 2016.
    real, pointer, dimension(:,:) :: ancsvgac_t
    real, pointer, dimension(:,:) :: ancgvgac_t
    real, pointer, dimension(:,:) :: rmlcsvga_t
    real, pointer, dimension(:,:) :: rmlcgvga_t
    integer, pointer, dimension(:) :: ipeatlandgat

    integer :: i,j
    real :: fsstar_gat
    real :: flstar_gat

    anmossac_t        => ctem_tile%anmossac_t
    rmlmossac_t       => ctem_tile%rmlmossac_t
    gppmossac_t       => ctem_tile%gppmossac_t
    altotcount_ctm    => vgat%altotcount_ctm
    fsinacc_gat       => vgat%fsinacc_gat
    flutacc_gat       => vgat%flutacc_gat
    flinacc_gat       => vgat%flinacc_gat
!     pregacc_gat       => vgat%pregacc_gat
    altotacc_gat      => vgat%altotacc_gat
    netrad_gat        => vgat%netrad_gat
!     preacc_gat        => vgat%preacc_gat
    ipeatlandgat      => vgat%ipeatland
    tbaraccgat_t      => ctem_tile%tbaraccgat_t
    tbarcacc_t        => ctem_tile%tbarcacc_t
    tbarcsacc_t       => ctem_tile%tbarcsacc_t
    tbargacc_t        => ctem_tile%tbargacc_t
    tbargsacc_t       => ctem_tile%tbargsacc_t
    thliqcacc_t       => ctem_tile%thliqcacc_t
    thliqgacc_t       => ctem_tile%thliqgacc_t
    thliqacc_t        => ctem_tile%thliqacc_t
    thiceacc_t        => ctem_tile%thiceacc_t  ! Added in place of YW's thicaccgat_m. EC Dec 23 2016.
    thicecacc_t       => ctem_tile%thicecacc_t
    thicegacc_t       => ctem_tile%thicegacc_t
    ancsvgac_t        => ctem_tile%ancsvgac_t
    ancgvgac_t        => ctem_tile%ancgvgac_t
    rmlcsvga_t        => ctem_tile%rmlcsvga_t
    rmlcgvga_t        => ctem_tile%rmlcgvga_t
    fsnowacc_t        => ctem_tile%fsnowacc_t
    tcansacc_t        => ctem_tile%tcansacc_t
    tcanoaccgat_t     => ctem_tile%tcanoaccgat_t
    taaccgat_t        => ctem_tile%taaccgat_t
    uvaccgat_t        => ctem_tile%uvaccgat_t
    vvaccgat_t        => ctem_tile%vvaccgat_t


    do i=1,nml

        !net radiation and precipitation estimates for ctem's bioclim

        uvaccgat_t(i)=uvaccgat_t(i)/real(nday)
        vvaccgat_t(i)=vvaccgat_t(i)/real(nday)
        fsinacc_gat(i)=fsinacc_gat(i)/real(nday)
        flinacc_gat(i)=flinacc_gat(i)/real(nday)
        flutacc_gat(i)=flutacc_gat(i)/real(nday)

        if (altotcount_ctm(i) > 0) then
            altotacc_gat(i)=altotacc_gat(i)/real(altotcount_ctm(i))
        else
            altotacc_gat(i)=0.
        end if

        fsstar_gat=fsinacc_gat(i)*(1.-altotacc_gat(i))
        flstar_gat=flinacc_gat(i)-flutacc_gat(i)
        netrad_gat(i)=fsstar_gat+flstar_gat
        !preacc_gat(i)=pregacc_gat(i)

        fsnowacc_t(i)=fsnowacc_t(i)/real(nday)
        tcanoaccgat_t(i)=tcanoaccgat_t(i)/real(nday)
        tcansacc_t(i)=tcansacc_t(i)/real(nday)
        taaccgat_t(i)=taaccgat_t(i)/real(nday)

        do 831 j=1,ignd
            tbaraccgat_t(i,j)=tbaraccgat_t(i,j)/real(nday)
            tbarcacc_t(i,j) = tbaraccgat_t(i,j)
            tbarcsacc_t(i,j) = tbaraccgat_t(i,j)
            tbargacc_t(i,j) = tbaraccgat_t(i,j)
            tbargsacc_t(i,j) = tbaraccgat_t(i,j)
            !
            thliqcacc_t(i,j)=thliqcacc_t(i,j)/real(nday)
            thliqgacc_t(i,j)=thliqgacc_t(i,j)/real(nday)
            thicecacc_t(i,j)=thicecacc_t(i,j)/real(nday)
            thicegacc_t(i,j)=thicegacc_t(i,j)/real(nday) ! EC Jan 31 2017.
            thliqacc_t(i,j)=thliqacc_t(i,j)/real(nday) ! Assume this replaces YW's thlqaccgat_m.
            thiceacc_t(i,j)=thiceacc_t(i,j)/real(nday) ! Added in place of YW's thicaccgat_m. EC Dec 23 2016.
    831 continue

        do 832 j = 1, icc
            ancsvgac_t(i,j)=ancsvgac_t(i,j)/real(nday)
            ancgvgac_t(i,j)=ancgvgac_t(i,j)/real(nday)
            rmlcsvga_t(i,j)=rmlcsvga_t(i,j)/real(nday)
            rmlcgvga_t(i,j)=rmlcgvga_t(i,j)/real(nday)
    832 continue

        !     -daily average moss C fluxes for ctem.f-------------------\
        !     Capitulum biomass = 0.22 kg/m2 in hummock, 0.1 kg/m2 in lawn
        !     stem biomass = 1.65 kg/m2 in hummock , 0.77 kg/m2 in lawn (Bragazza et al.2004)
        !     the ratio between stem and capitulum = 7.5 and 7.7
        if (ipeatlandgat(i) > 0) then
            anmossac_t(i) = anmossac_t(i)/real(nday)
            rmlmossac_t(i)= rmlmossac_t(i)/real(nday)
            gppmossac_t(i) = gppmossac_t(i)/real(nday)
        endif

    end do !nml loop

end subroutine dayEndCTEMPreparation

! --------------------------------------------------------------------------------------------------------------------
!> Accumulate the CTEM input (physics) variables at the end of each physics timestep
subroutine accumulateForCTEM(nml)

    use ctem_params, only : icc,ignd
    use class_statevars, only : class_gat,class_rot
    use ctem_statevars, only : vgat,ctem_tile

    implicit none

    integer, intent(in) :: nml     !<Counter representing number of mosaic tiles on modelled domain that are land
    !integer, intent(in) :: delt    !<Timestep in seconds

    integer :: i,j

    real :: delt,sbc

    COMMON /CLASS1/ DELT
    COMMON /CLASS2/ SBC

    real, pointer, dimension(:) :: FSIHGAT !<Near-infrared radiation incident on horizontal surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: FSVHGAT !<Visible radiation incident on horizontal surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: ancsmoss
    real, pointer, dimension(:) :: angsmoss
    real, pointer, dimension(:) :: ancmoss
    real, pointer, dimension(:) :: angmoss
    real, pointer, dimension(:) :: rmlcsmoss
    real, pointer, dimension(:) :: rmlgsmoss
    real, pointer, dimension(:) :: rmlcmoss
    real, pointer, dimension(:) :: rmlgmoss
    real, pointer, dimension(:) :: ALIRGAT !<Diagnosed total near-infrared albedo of land surface [ ]
    real, pointer, dimension(:) :: ALVSGAT !<Diagnosed total visible albedo of land surface [ ]
    real, pointer, dimension(:) :: FSSROW  !< Shortwave radiation \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: GTGAT   !<Diagnosed effective surface black-body temperature [K]
    real, pointer, dimension(:) :: FDLGAT  !<Downwelling longwave radiation at bottom of atmosphere (i.e. incident on modelled land surface elements \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: PREGAT  !<Surface precipitation rate \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: FSNOGAT !<Diagnosed fractional snow coverage [ ]
    real, pointer, dimension(:) :: TCANO   !<
    real, pointer, dimension(:) :: TCANS   !<
    real, pointer, dimension(:) :: TAGAT   !<Air temperature at reference height [K]
    real, pointer, dimension(:) :: VLGAT   !<Meridional component of wind velocity \f$[m s^{-1} ]\f$
    real, pointer, dimension(:) :: ULGAT   !<Zonal component of wind velocity \f$[m s^{-1} ]\f$
    real, pointer, dimension(:) :: FSGGGAT !<Diagnosed net shortwave radiation at soil surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:) :: TBARGAT !<Temperature of soil layers [K]
    real, pointer, dimension(:) :: FSGSGAT !<Diagnosed net shortwave radiation at snow surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:) :: TBARC  !<
    real, pointer, dimension(:,:) :: TBARG  !<
    real, pointer, dimension(:,:) :: TBARCS !<
    real, pointer, dimension(:,:) :: TBARGS !<
    real, pointer, dimension(:,:) :: THLIQC !<
    real, pointer, dimension(:,:) :: THLIQG !<
    real, pointer, dimension(:,:) :: THICEC !<
    real, pointer, dimension(:,:) :: THICEG !<
    real, pointer, dimension(:) :: FSGVGAT !<Diagnosed net shortwave radiation on vegetation canopy \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:) :: THICGAT !<Volumetric frozen water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, pointer, dimension(:,:) :: THLQGAT !<Volumetric liquid water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, pointer, dimension(:,:) :: ancsveggat
    real, pointer, dimension(:,:) :: ancgveggat
    real, pointer, dimension(:,:) :: rmlcsveggat
    real, pointer, dimension(:,:) :: rmlcgveggat
    real, pointer, dimension(:) :: gppmossgat
    real, pointer, dimension(:) :: anmossgat
    real, pointer, dimension(:) :: rmlmossgat
    real, pointer, dimension(:) :: FC      !<
    real, pointer, dimension(:) :: FG      !<
    real, pointer, dimension(:) :: FCS     !<
    real, pointer, dimension(:) :: FGS     !<
    integer, pointer, dimension(:) :: altotcount_ctm !nlat  !< Counter used for calculating total albedo
    real, pointer, dimension(:)    :: fsinacc_gat !(ilg)    !<
    real, pointer, dimension(:)    :: flutacc_gat !(ilg)    !<
    real, pointer, dimension(:)    :: flinacc_gat !(ilg)    !<
    real, pointer, dimension(:)    :: altotacc_gat !(ilg)   !<
    !real, pointer, dimension(:)    :: pregacc_gat !(ilg)    !<
    real, pointer, dimension(:)    :: netrad_gat !(ilg)     !<
    real, pointer, dimension(:)    :: preacc_gat !(ilg)     !<
    real, pointer, dimension(:) :: anmossac_t
    real, pointer, dimension(:) :: rmlmossac_t
    real, pointer, dimension(:) :: gppmossac_t
    real, pointer, dimension(:) :: fsnowacc_t       !<
    real, pointer, dimension(:) :: tcansacc_t       !<
    real, pointer, dimension(:) :: tcanoaccgat_t    !<
    real, pointer, dimension(:) :: taaccgat_t       !<
    real, pointer, dimension(:) :: uvaccgat_t       !<
    real, pointer, dimension(:) :: vvaccgat_t       !<
    real, pointer, dimension(:,:) :: tbaraccgat_t!<
    real, pointer, dimension(:,:) :: tbarcacc_t  !<
    real, pointer, dimension(:,:) :: tbarcsacc_t !<
    real, pointer, dimension(:,:) :: tbargacc_t  !<
    real, pointer, dimension(:,:) :: tbargsacc_t !<
    real, pointer, dimension(:,:) :: thliqcacc_t !<
    real, pointer, dimension(:,:) :: thliqgacc_t !<
    real, pointer, dimension(:,:) :: thliqacc_t  !<
    real, pointer, dimension(:,:) :: thicecacc_t !<
    real, pointer, dimension(:,:) :: thicegacc_t !<
    real, pointer, dimension(:,:) :: thiceacc_t  !< Added in place of YW's thicaccgat_m. EC Dec 23 2016.
    real, pointer, dimension(:,:) :: ancsvgac_t
    real, pointer, dimension(:,:) :: ancgvgac_t
    real, pointer, dimension(:,:) :: rmlcsvga_t
    real, pointer, dimension(:,:) :: rmlcgvga_t
    integer, pointer, dimension(:) :: ipeatlandgat


    anmossac_t        => ctem_tile%anmossac_t
    rmlmossac_t       => ctem_tile%rmlmossac_t
    gppmossac_t       => ctem_tile%gppmossac_t
    altotcount_ctm    => vgat%altotcount_ctm
    fsinacc_gat       => vgat%fsinacc_gat
    flutacc_gat       => vgat%flutacc_gat
    flinacc_gat       => vgat%flinacc_gat
    !pregacc_gat       => vgat%pregacc_gat
    altotacc_gat      => vgat%altotacc_gat
    netrad_gat        => vgat%netrad_gat
    preacc_gat        => vgat%preacc_gat
    ipeatlandgat      => vgat%ipeatland
    tbaraccgat_t      => ctem_tile%tbaraccgat_t
    tbarcacc_t        => ctem_tile%tbarcacc_t
    tbarcsacc_t       => ctem_tile%tbarcsacc_t
    tbargacc_t        => ctem_tile%tbargacc_t
    tbargsacc_t       => ctem_tile%tbargsacc_t
    thliqcacc_t       => ctem_tile%thliqcacc_t
    thliqgacc_t       => ctem_tile%thliqgacc_t
    thliqacc_t        => ctem_tile%thliqacc_t
    thiceacc_t        => ctem_tile%thiceacc_t  ! Added in place of YW's thicaccgat_m. EC Dec 23 2016.
    thicecacc_t       => ctem_tile%thicecacc_t
    thicegacc_t       => ctem_tile%thicegacc_t
    ancsvgac_t        => ctem_tile%ancsvgac_t
    ancgvgac_t        => ctem_tile%ancgvgac_t
    rmlcsvga_t        => ctem_tile%rmlcsvga_t
    rmlcgvga_t        => ctem_tile%rmlcgvga_t
    fsnowacc_t        => ctem_tile%fsnowacc_t
    tcansacc_t        => ctem_tile%tcansacc_t
    tcanoaccgat_t     => ctem_tile%tcanoaccgat_t
    taaccgat_t        => ctem_tile%taaccgat_t
    uvaccgat_t        => ctem_tile%uvaccgat_t
    vvaccgat_t        => ctem_tile%vvaccgat_t
    ancsmoss         => vgat%ancsmoss
    angsmoss         => vgat%angsmoss
    ancmoss          => vgat%ancmoss
    angmoss          => vgat%angmoss
    rmlcsmoss        => vgat%rmlcsmoss
    rmlgsmoss        => vgat%rmlgsmoss
    rmlcmoss         => vgat%rmlcmoss
    rmlgmoss         => vgat%rmlgmoss
    FC => class_gat%FC
    FG => class_gat%FG
    FCS => class_gat%FCS
    FGS => class_gat%FGS
    FSIHGAT => class_gat%FSIHGAT
    FSVHGAT => class_gat%FSVHGAT
    ALIRGAT => class_gat%ALIRGAT
    ALVSGAT => class_gat%ALVSGAT
    FSSROW => class_rot%FSSROW
    GTGAT => class_gat%GTGAT
    FDLGAT => class_gat%FDLGAT
    PREGAT => class_gat%PREGAT
    FSNOGAT => class_gat%FSNOGAT
    TCANO => class_gat%TCANO
    TCANS => class_gat%TCANS
    TAGAT => class_gat%TAGAT
    VLGAT => class_gat%VLGAT
    ULGAT => class_gat%ULGAT
    FSGGGAT => class_gat%FSGGGAT
    TBARGAT => class_gat%TBARGAT
    FSGSGAT => class_gat%FSGSGAT
    TBARC => class_gat%TBARC
    TBARG => class_gat%TBARG
    TBARCS => class_gat%TBARCS
    TBARGS => class_gat%TBARGS
    THLIQC => class_gat%THLIQC
    THLIQG => class_gat%THLIQG
    THICEC => class_gat%THICEC
    THICEG => class_gat%THICEG
    FSGVGAT => class_gat%FSGVGAT
    THICGAT => class_gat%THICGAT
    THLQGAT => class_gat%THLQGAT
    ancsveggat        => vgat%ancsveg
    ancgveggat        => vgat%ancgveg
    rmlcsveggat       => vgat%rmlcsveg
    rmlcgveggat       => vgat%rmlcgveg
    anmossgat        => vgat%anmoss
    rmlmossgat       => vgat%rmlmoss
    gppmossgat       => vgat%gppmoss

    do i = 1, nml

        fsinacc_gat(i)=fsinacc_gat(i)+FSSROW(1) ! FLAG! Do this offline only (since all tiles are the same in a gridcell and we run
                                                ! only one gridcell at a time. JM Feb 42016.
        flinacc_gat(i)=flinacc_gat(i)+fdlgat(i)
        flutacc_gat(i)=flutacc_gat(i)+sbc*gtgat(i)**4
        preacc_gat(i)=preacc_gat(i)+pregat(i)*delt
!         pregacc_gat(i)=pregacc_gat(i)+pregat(i)*delt
        fsnowacc_t(i)=fsnowacc_t(i)+fsnogat(i)
        tcanoaccgat_t(i)=tcanoaccgat_t(i)+tcano(i)
        tcansacc_t(i)=tcansacc_t(i)+tcans(i)
        taaccgat_t(i)=taaccgat_t(i)+tagat(i)
        vvaccgat_t(i)=vvaccgat_t(i)+ vlgat(i)
        uvaccgat_t(i)=uvaccgat_t(i)+ulgat(i)
        if (FSSROW(I) .gt. 0.) then
            altotacc_gat(i) = altotacc_gat(i) + (FSSROW(I)-&
                &                (FSGVGAT(I)+FSGSGAT(I)+FSGGGAT(I)))&
                &                /FSSROW(I)
            altotcount_ctm = altotcount_ctm + 1
        end if

        do 710 j=1,ignd
            tbaraccgat_t(i,j)=tbaraccgat_t(i,j)+tbargat(i,j)
            tbarcacc_t(i,j)=tbarcacc_t(i,j)+tbarc(i,j)
            tbarcsacc_t(i,j)=tbarcsacc_t(i,j)+tbarcs(i,j)
            tbargacc_t(i,j)=tbargacc_t(i,j)+tbarg(i,j)
            tbargsacc_t(i,j)=tbargsacc_t(i,j)+tbargs(i,j)
            thliqcacc_t(i,j)=thliqcacc_t(i,j)+thliqc(i,j)
            thliqgacc_t(i,j)=thliqgacc_t(i,j)+thliqg(i,j)
            thicecacc_t(i,j)=thicecacc_t(i,j)+thicec(i,j)
            ! FLAG: Needs to be reviewed. EC Dec 23 2016.
            ! YW's original variables were thlqaccgat_m/thicaccgat_m.
            ! The following 2 variables needs to be passed via ctem to hetres_peat
            ! (otherwise causes a floating exception at line 863 of peatlands_mod).
            thliqacc_t(i,j) = thliqacc_t(i,j) + THLQGAT(i,j) ! Not used elsewhere, so assume replacement for thlqaccgat_m
            thiceacc_t(i,j) = thiceacc_t(i,j) + THICGAT(i,j) ! New.
            thicegacc_t(i,j)=thicegacc_t(i,j)+thiceg(i,j)    ! EC Jan 31 2017.
710                 continue

        do 713 j = 1, icc
            ancsvgac_t(i,j)=ancsvgac_t(i,j)+ancsveggat(i,j)
            ancgvgac_t(i,j)=ancgvgac_t(i,j)+ancgveggat(i,j)
            rmlcsvga_t(i,j)=rmlcsvga_t(i,j)+rmlcsveggat(i,j)
            rmlcgvga_t(i,j)=rmlcgvga_t(i,j)+rmlcgveggat(i,j)
713                 continue

        !    -accumulate moss C fluxes to tile level then daily----
        if (ipeatlandgat(i) > 0) then
            anmossgat(i) = fcs(i)*ancsmoss(i)+fgs(i)*angsmoss(i)+fc(i)*ancmoss(i)+fg(i)*angmoss(i)
            rmlmossgat(i)= fcs(i)*rmlcsmoss(i)+fgs(i)*rmlgsmoss(i)+fc(i)*rmlcmoss(i)+fg(i)*rmlgmoss(i)
            gppmossgat(i) = anmossgat(i) +rmlmossgat(i)

            anmossac_t(i) = anmossac_t(i)   + anmossgat(i)
            rmlmossac_t(i)= rmlmossac_t(i)  + rmlmossgat(i)
            gppmossac_t(i)= gppmossac_t(i)  + gppmossgat(i)

        endif
    end do

end subroutine accumulateForCTEM
!
! --------------------------------------------------------------------------------------------------------------------
!> Find mosaic tile (grid) average vegetation biomass, litter mass, and soil c mass.
!! Also initialize additional variables which are used by CTEM (biogeochemical processes).
subroutine ctemInit(nltest,nmtest)

    use ctem_params, only : icc,ilg
    use ctem_statevars, only : vrot,ctem_tile,vgat
    use class_statevars,only : class_rot

    implicit none

    integer, intent(in) :: nltest
    integer, intent(in) :: nmtest

    integer :: i,m,j

    real, pointer, dimension(:,:,:) :: co2i1cgrow
    real, pointer, dimension(:,:,:) :: co2i1csrow
    real, pointer, dimension(:,:,:) :: co2i2cgrow
    real, pointer, dimension(:,:,:) :: co2i2csrow
    real, pointer, dimension(:,:,:) :: slairow
    real, pointer, dimension(:) :: fsnowacc_t
    real, pointer, dimension(:) :: tcansacc_t
    real, pointer, dimension(:) :: taaccgat_t
    real, pointer, dimension(:,:) :: ancsvgac_t
    real, pointer, dimension(:,:) :: ancgvgac_t
    real, pointer, dimension(:,:) :: rmlcsvga_t
    real, pointer, dimension(:,:) :: rmlcgvga_t
    real, pointer, dimension(:,:)  :: todfrac
    real, pointer, dimension(:,:) :: tbarcacc_t
    real, pointer, dimension(:,:) :: tbarcsacc_t
    real, pointer, dimension(:,:) :: tbargacc_t
    real, pointer, dimension(:,:) :: tbargsacc_t
    real, pointer, dimension(:,:) :: thliqcacc_t
    real, pointer, dimension(:,:) :: thliqgacc_t
    real, pointer, dimension(:,:) :: thliqacc_t
    real, pointer, dimension(:,:) :: thiceacc_t  ! Added in place of YW's thicaccgat_m. EC Dec 23 2016.
    real, pointer, dimension(:,:) :: thicecacc_t
    real, pointer, dimension(:,:) :: thicegacc_t
    real, pointer, dimension(:,:) :: vgbiomasrow
    real, pointer, dimension(:,:) :: gavgltmsrow
    real, pointer, dimension(:,:) :: gavgscmsrow
    real, pointer, dimension(:,:) :: gavglairow
    real, pointer, dimension(:,:) :: lucemcomrow
    real, pointer, dimension(:,:) :: lucltrinrow
    real, pointer, dimension(:,:) :: lucsocinrow
    integer, pointer, dimension(:,:,:) :: colddaysrow
    real, pointer, dimension(:,:,:) :: flhrlossrow
    real, pointer, dimension(:,:,:) :: stmhrlosrow
    real, pointer, dimension(:,:,:) :: rothrlosrow
    real, pointer, dimension(:,:,:) :: tymaxlairow
    real, pointer, dimension(:,:,:) :: grwtheffrow
    real, pointer, dimension(:,:,:) :: lystmmasrow
    real, pointer, dimension(:,:,:) :: lyrotmasrow
    integer, pointer, dimension(:,:) :: ipeatlandrow !This is first set in read_from_ctm.
    real, pointer, dimension(:,:,:) :: fcanrot !<
    real, pointer, dimension(:,:,:) :: gleafmasrow        !
    real, pointer, dimension(:,:,:) :: bleafmasrow        !
    real, pointer, dimension(:,:,:) :: stemmassrow        !
    real, pointer, dimension(:,:,:) :: rootmassrow        !
    real, pointer, dimension(:,:,:) :: litrmassrow
    real, pointer, dimension(:,:,:) :: soilcmasrow
    real, pointer, dimension(:,:,:) :: fcancmxrow
    real, pointer, dimension(:,:) :: peatdeprow
    real, pointer, dimension(:) :: anmossac_t
    real, pointer, dimension(:) :: rmlmossac_t
    real, pointer, dimension(:) :: gppmossac_t
    real, pointer, dimension(:,:) :: litrmsmossrow
    real, pointer, dimension(:,:) :: Cmossmasrow
    real, pointer, dimension(:,:) :: sdeprot !<Depth to bedrock in the soil profile

    ipeatlandrow     => vrot%ipeatland
    co2i1cgrow        => vrot%co2i1cg
    co2i1csrow        => vrot%co2i1cs
    co2i2cgrow        => vrot%co2i2cg
    co2i2csrow        => vrot%co2i2cs
    slairow           => vrot%slai
    fsnowacc_t        => ctem_tile%fsnowacc_t
    tcansacc_t        => ctem_tile%tcansacc_t
    taaccgat_t        => ctem_tile%taaccgat_t
    ancsvgac_t        => ctem_tile%ancsvgac_t
    ancgvgac_t        => ctem_tile%ancgvgac_t
    rmlcsvga_t        => ctem_tile%rmlcsvga_t
    rmlcgvga_t        => ctem_tile%rmlcgvga_t
    todfrac           => vgat%todfrac
    tbarcacc_t        => ctem_tile%tbarcacc_t
    tbarcsacc_t       => ctem_tile%tbarcsacc_t
    tbargacc_t        => ctem_tile%tbargacc_t
    tbargsacc_t       => ctem_tile%tbargsacc_t
    thliqcacc_t       => ctem_tile%thliqcacc_t
    thliqgacc_t       => ctem_tile%thliqgacc_t
    thliqacc_t        => ctem_tile%thliqacc_t
    thiceacc_t        => ctem_tile%thiceacc_t
    thicecacc_t       => ctem_tile%thicecacc_t
    thicegacc_t       => ctem_tile%thicegacc_t
    vgbiomasrow       => vrot%vgbiomas
    gavglairow        => vrot%gavglai
    gavgltmsrow       => vrot%gavgltms
    gavgscmsrow       => vrot%gavgscms
    lucemcomrow       => vrot%lucemcom
    lucltrinrow       => vrot%lucltrin
    lucsocinrow       => vrot%lucsocin
    colddaysrow       => vrot%colddays
    flhrlossrow       => vrot%flhrloss
    stmhrlosrow       => vrot%stmhrlos
    rothrlosrow       => vrot%rothrlos
    tymaxlairow       => vrot%tymaxlai
    grwtheffrow       => vrot%grwtheff
    lystmmasrow       => vrot%lystmmas
    lyrotmasrow       => vrot%lyrotmas
    fcanrot           => class_rot%fcanrot
    gleafmasrow       => vrot%gleafmas
    bleafmasrow       => vrot%bleafmas
    stemmassrow       => vrot%stemmass
    rootmassrow       => vrot%rootmass
    litrmassrow       => vrot%litrmass
    soilcmasrow       => vrot%soilcmas
    fcancmxrow        => vrot%fcancmx
    peatdeprow        => vrot%peatdep
    anmossac_t        => ctem_tile%anmossac_t
    rmlmossac_t       => ctem_tile%rmlmossac_t
    gppmossac_t       => ctem_tile%gppmossac_t
    litrmsmossrow    => vrot%litrmsmoss
    Cmossmasrow      => vrot%Cmossmas
    sdeprot => class_rot%sdeprot
    ! ------

    ! Initialize to zero:
    co2i1csrow(:,:,:)=0.0     !intercellular co2 concentrations
    co2i1cgrow(:,:,:)=0.0
    co2i2csrow(:,:,:)=0.0
    co2i2cgrow(:,:,:)=0.0
    slairow(:,:,:)=0.0        !if bio2str is not called we need to initialize this to zero
    fsnowacc_t(:)=0.0         !daily accu. fraction of snow
    tcansacc_t(:)=0.0         !daily accu. canopy temp. over snow
    taaccgat_t(:)=0.0
    ancsvgac_t(:,:)=0.0    !daily accu. net photosyn. for canopy over snow subarea
    ancgvgac_t(:,:)=0.0    !daily accu. net photosyn. for canopy over ground subarea
    rmlcsvga_t(:,:)=0.0    !daily accu. leaf respiration for canopy over snow subarea
    rmlcgvga_t(:,:)=0.0    !daily accu. leaf respiration for canopy over ground subarea
    todfrac(:,:)=0.0
    tbarcacc_t (:,:)=0.0
    tbarcsacc_t(:,:)=0.0
    tbargacc_t (:,:)=0.0
    tbargsacc_t(:,:)=0.0
    thliqcacc_t(:,:)=0.0
    thliqgacc_t(:,:)=0.0
    thliqacc_t(:,:)=0.0
    thiceacc_t(:,:)=0.0
    thicecacc_t(:,:)=0.0
    thicegacc_t(:,:)=0.0
    vgbiomasrow(:,:)=0.0
    gavglairow(:,:)=0.0
    gavgltmsrow(:,:)=0.0
    gavgscmsrow(:,:)=0.0
    lucemcomrow(:,:)=0.0      !land use change combustion emission losses
    lucltrinrow(:,:)=0.0      !land use change inputs to litter pool
    lucsocinrow(:,:)=0.0      !land use change inputs to soil c pool
    colddaysrow(:,:,:)=0      !cold days counter for ndl dcd and crops
    flhrlossrow(:,:,:)=0.0     !fall/harvest loss
    stmhrlosrow(:,:,:)=0.0     !stem harvest loss for crops
    rothrlosrow(:,:,:)=0.0     !root death for crops
    tymaxlairow(:,:,:)=0.0

    do 115 i = 1,nltest
        do 115 m = 1,nmtest
            do 116 j = 1, icc
                vgbiomasrow(i,m)=vgbiomasrow(i,m)+fcancmxrow(i,m,j)*&
                    &        (gleafmasrow(i,m,j)+stemmassrow(i,m,j)+&
                    &         rootmassrow(i,m,j)+bleafmasrow(i,m,j))
                gavgltmsrow(i,m)=gavgltmsrow(i,m)+fcancmxrow(i,m,j)*&
                    &                       litrmassrow(i,m,j)
                gavgscmsrow(i,m)=gavgscmsrow(i,m)+fcancmxrow(i,m,j)*&
                    &         soilcmasrow(i,m,j)
                grwtheffrow(i,m,j)=100.0   !set growth efficiency to some large number
                                                        !so that no growth related mortality occurs in
                                                        !first year
                lystmmasrow(i,m,j)=stemmassrow(i,m,j)
                lyrotmasrow(i,m,j)=rootmassrow(i,m,j)

116                 continue
115         continue

    do 117 i = 1,nltest
        do 117 m = 1,nmtest
            if (ipeatlandrow(i,m)==0) then ! NON-peatland tile
                gavgltmsrow(i,m)=gavgltmsrow(i,m)+ (1.0-fcanrot(i,m,1)-&
                        fcanrot(i,m,2)-fcanrot(i,m,3)-&
                        fcanrot(i,m,4))*litrmassrow(i,m,icc+1)
                gavgscmsrow(i,m)=gavgscmsrow(i,m)+ (1.0-fcanrot(i,m,1)-&
                        fcanrot(i,m,2)-fcanrot(i,m,3)-&
                        fcanrot(i,m,4))*soilcmasrow(i,m,icc+1)
            else !peatland tile
                gavgltmsrow(i,m)= gavgltmsrow(i,m)+litrmsmossrow(i,m)
                peatdeprow(i,m) = sdeprot(i,m) !the peatdepth is set to the soil depth
                ! The soil carbon on the peatland tiles is assigned based on depth. This
                ! is the same relation as found in decp subroutine.
                gavgscmsrow(i,m) = 0.487*(4056.6*peatdeprow(i,m)**2+&
                            72067.0*peatdeprow(i,m))/1000
                vgbiomasrow(i,m)=vgbiomasrow(i,m)+Cmossmasrow(i,m)
            endif
117         continue

        !    Also initialize the accumulators for moss daily C fluxes.
        !    FLAG perhaps move? JM.
        !
        !    Moved out of 117 loop and ipeatlandrow > 0 block where only 1 to nltest
        !    elements were actually initialized since the loop is really for row
        !    variables. Initialization is now done on all elements, which shouldn't be
        !    a problem. EC - Feb 16, 2016.

        anmossac_t  = 0.0
        rmlmossac_t = 0.0
        gppmossac_t = 0.0

            !write(6,6990) 'peatdeprow=', peatdeprow
            !write(6,6990) 'gavgscms=', gavgscmsrow
            !write(6,6990) 'vgbiomas=', vgbiomasrow
            !6990   format(A15,12f6.2)
        !    ----------------------------YW March 25, 2015 --------------------/
end subroutine ctemInit

end module ctemUtilities
