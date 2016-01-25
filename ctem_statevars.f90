module ctem_statevars

!this module contains the variable type structures:
! 1) c_switch - switches for running CTEM, read from the joboptions file
! 2) vrot - CTEM's 'rot' vars
! 3) vgat - CTEM's 'gat' vars
! 4) class_out - CLASS's monthly outputs
! 5) ctem_grd - CTEM's grid average variables
! 6) ctem_tile - CTEM's variables per tile (mosaic)
! 7) ctem_grd_mo - CTEM's grid average monthly values
! 8) ctem_tile_mo - CTEM's variables per tile (mosaic) monthly values
! 9) ctem_grd_yr - CTEM's grid average annual values
! 10) ctem_tile_yr - CTEM's variables per tile (mosaic) annual values

! J. Melton Apr 2015

use ctem_params,  only : initpftpars, nlat, nmos, ilg, nmon,ican, ignd,icp1, icc, iccp1, &
                    monthend, mmday,modelpft, l2max,deltat, abszero, monthdays,seed, crop, NBS

implicit none

public :: initrowvars
public :: resetclassmon
public :: resetclassyr
public :: resetmidmonth
public :: resetmonthend_g
public :: resetmonthend_m
public :: resetyearend_g
public :: resetyearend_m
public :: resetclassaccum
public :: resetgridavg
public :: resetctem_g
public :: finddaylength

!=================================================================================
type ctem_switches

    logical :: ctem_on
    logical :: parallelrun
    logical :: mosaic
    logical :: cyclemet
    logical :: dofire
    logical :: run_model
    logical :: met_rewound
    logical :: reach_eof
    logical :: compete
    logical :: start_bare
    logical :: rsfile
    logical :: lnduseon
    logical :: co2on
    logical :: ch4on
    logical :: popdon
    logical :: inibioclim
    logical :: start_from_rs
    logical :: dowetlands
    logical :: obswetf
    logical :: transient_run

    character(80) :: titlec1
    character(80) :: titlec2
    character(80) :: titlec3

end type ctem_switches

type (ctem_switches), save, target :: c_switch

!=================================================================================

type veg_rot

    ! This is the basic data structure that contains the state variables
    ! for the Plant Functional Type (PFT). The dimensions are nlat,nmos,{icc,iccp1}

    real, dimension(nlat,nmos,icc) :: ailcmin         !
    real, dimension(nlat,nmos,icc) :: ailcmax         !
    real, dimension(nlat,nmos,icc) :: dvdfcan         !
    real, dimension(nlat,nmos,icc) :: gleafmas        !
    real, dimension(nlat,nmos,icc) :: bleafmas        !
    real, dimension(nlat,nmos,icc) :: stemmass        !
    real, dimension(nlat,nmos,icc) :: rootmass        !
    real, dimension(nlat,nmos,icc) :: pstemmass       !
    real, dimension(nlat,nmos,icc) :: pgleafmass      !
    real, dimension(nlat,nmos,icc) :: fcancmx

    real, dimension(nlat,nmos) :: gavglai

    real, dimension(nlat,nmos,ican) :: zolnc
    real, dimension(nlat,nmos,ican) :: ailc

    real, dimension(nlat,nmos,icc) :: ailcg
    real, dimension(nlat,nmos,icc) :: ailcgs
    real, dimension(nlat,nmos,icc) :: fcancs
    real, dimension(nlat,nmos,icc) :: fcanc

    real, dimension(nlat,nmos) :: co2conc
    real, dimension(nlat,nmos) :: ch4conc

    real, dimension(nlat,nmos,icc) :: co2i1cg
    real, dimension(nlat,nmos,icc) :: co2i1cs
    real, dimension(nlat,nmos,icc) :: co2i2cg
    real, dimension(nlat,nmos,icc) :: co2i2cs
    real, dimension(nlat,nmos,icc) :: ancsveg
    real, dimension(nlat,nmos,icc) :: ancgveg
    real, dimension(nlat,nmos,icc) :: rmlcsveg
    real, dimension(nlat,nmos,icc) :: rmlcgveg
    real, dimension(nlat,nmos,icc) :: slai
    real, dimension(nlat,nmos,icc) :: ailcb
    real, dimension(nlat,nmos) :: canres
    real, dimension(nlat,nmos,icc) :: flhrloss

    real, dimension(nlat,nmos,icc) :: grwtheff
    real, dimension(nlat,nmos,icc) :: lystmmas
    real, dimension(nlat,nmos,icc) :: lyrotmas
    real, dimension(nlat,nmos,icc) :: tymaxlai
    real, dimension(nlat,nmos) :: vgbiomas
    real, dimension(nlat,nmos) :: gavgltms
    real, dimension(nlat,nmos) :: gavgscms
    real, dimension(nlat,nmos,icc) :: stmhrlos
    real, dimension(nlat,nmos,ican,ignd) :: rmatc
    real, dimension(nlat,nmos,icc,ignd) :: rmatctem
    real, dimension(nlat,nmos,iccp1) :: litrmass
    real, dimension(nlat,nmos,iccp1) :: soilcmas
    real, dimension(nlat,nmos,icc) :: vgbiomas_veg

    ! c     Fire-related variables

    real, dimension(nlat,nmos,icc) :: emit_co2
    real, dimension(nlat,nmos,icc) :: emit_co
    real, dimension(nlat,nmos,icc) :: emit_ch4
    real, dimension(nlat,nmos,icc) :: emit_nmhc
    real, dimension(nlat,nmos,icc) :: emit_h2
    real, dimension(nlat,nmos,icc) :: emit_nox
    real, dimension(nlat,nmos,icc) :: emit_n2o
    real, dimension(nlat,nmos,icc) :: emit_pm25
    real, dimension(nlat,nmos,icc) :: emit_tpm
    real, dimension(nlat,nmos,icc) :: emit_tc
    real, dimension(nlat,nmos,icc) :: emit_oc
    real, dimension(nlat,nmos,icc) :: emit_bc
    real, dimension(nlat,nmos) :: burnfrac
    real, dimension(nlat,nmos,icc) :: burnvegf
    real, dimension(nlat,nmos) :: probfire
    real, dimension(nlat,nmos) :: bterm
    real, dimension(nlat,nmos) :: lterm
    real, dimension(nlat,nmos) :: mterm

    real, dimension(nlat) :: extnprob
    real, dimension(nlat) :: prbfrhuc
    real, dimension(nlat,12) :: mlightng

    real, dimension(nlat,nmos,icc) :: bmasveg
    real, dimension(nlat,nmos,ican) :: cmasvegc
    real, dimension(nlat,nmos,icc) :: veghght
    real, dimension(nlat,nmos,icc) :: rootdpth
    real, dimension(nlat,nmos) :: rml
    real, dimension(nlat,nmos) :: rms
    real, dimension(nlat,nmos,icc) :: tltrleaf
    real, dimension(nlat,nmos,icc) :: tltrstem
    real, dimension(nlat,nmos,icc) :: tltrroot
    real, dimension(nlat,nmos,icc) :: leaflitr
    real, dimension(nlat,nmos,icc) :: roottemp
    real, dimension(nlat,nmos,icc) :: afrleaf
    real, dimension(nlat,nmos,icc) :: afrstem
    real, dimension(nlat,nmos,icc) :: afrroot
    real, dimension(nlat,nmos,icc) :: wtstatus
    real, dimension(nlat,nmos,icc) :: ltstatus
    real, dimension(nlat,nmos) :: rmr

    real, dimension(nlat) :: wetfrac
    real, dimension(nlat,nmos) :: ch4wet1
    real, dimension(nlat,nmos) :: ch4wet2
    real, dimension(nlat,nmos) :: wetfdyn
    real, dimension(nlat,nmos) :: ch4dyn1
    real, dimension(nlat,nmos) :: ch4dyn2
    real, dimension(nlat,12) :: wetfrac_mon
    real, dimension(nlat,nmos) :: ch4_soills

    real, dimension(nlat,nmos) :: lucemcom
    real, dimension(nlat,nmos) :: lucltrin
    real, dimension(nlat,nmos) :: lucsocin

    real, dimension(nlat,nmos) :: npp
    real, dimension(nlat,nmos) :: nep
    real, dimension(nlat,nmos) :: nbp
    real, dimension(nlat,nmos) :: gpp
    real, dimension(nlat,nmos) :: hetrores
    real, dimension(nlat,nmos) :: autores
    real, dimension(nlat,nmos) :: soilcresp
    real, dimension(nlat,nmos) :: rm
    real, dimension(nlat,nmos) :: rg
    real, dimension(nlat,nmos) :: litres
    real, dimension(nlat,nmos) :: socres
    real, dimension(nlat,nmos) :: dstcemls
    real, dimension(nlat,nmos) :: litrfall
    real, dimension(nlat,nmos) :: humiftrs

    real, dimension(nlat,nmos,icc) :: gppveg
    real, dimension(nlat,nmos,iccp1) :: nepveg
    real, dimension(nlat,nmos,iccp1) :: nbpveg
    real, dimension(nlat,nmos,icc) :: nppveg
    real, dimension(nlat,nmos,iccp1) :: hetroresveg
    real, dimension(nlat,nmos,icc) :: autoresveg
    real, dimension(nlat,nmos,iccp1) :: litresveg
    real, dimension(nlat,nmos,iccp1) :: soilcresveg
    real, dimension(nlat,nmos,icc) :: rmlvegacc
    real, dimension(nlat,nmos,icc) :: rmsveg
    real, dimension(nlat,nmos,icc) :: rmrveg
    real, dimension(nlat,nmos,icc) :: rgveg

    real, dimension(nlat,nmos,icc) :: rothrlos
    real, dimension(nlat,nmos,icc) :: pfcancmx
    real, dimension(nlat,nmos,icc) :: nfcancmx
    real, dimension(nlat,nmos,ican) :: alvsctm
    real, dimension(nlat,nmos,ican) :: paic
    real, dimension(nlat,nmos,ican) :: slaic
    real, dimension(nlat,nmos,ican) :: alirctm
    real, dimension(nlat,nmos) :: cfluxcg
    real, dimension(nlat,nmos) :: cfluxcs
    real, dimension(nlat,nmos) :: dstcemls3
    real, dimension(nlat,nmos,icc) :: anveg
    real, dimension(nlat,nmos,icc) :: rmlveg

    logical, dimension(nlat,nmos,icc) :: pftexist
    integer, dimension(nlat,nmos,2) :: colddays
    integer, dimension(nlat,nmos) :: icount
    integer, dimension(nlat,nmos,icc) :: lfstatus
    integer, dimension(nlat,nmos,icc) :: pandays
    integer, dimension(nlat) :: stdaln

    real, dimension(nlat,nmos) :: PREACC_M
    real, dimension(nlat,nmos) :: GTACC_M
    real, dimension(nlat,nmos) :: QEVPACC_M
    real, dimension(nlat,nmos) :: HFSACC_M
    real, dimension(nlat,nmos) :: HMFNACC_M
    real, dimension(nlat,nmos) :: ROFACC_M
    real, dimension(nlat,nmos) :: SNOACC_M
    real, dimension(nlat,nmos) :: OVRACC_M
    real, dimension(nlat,nmos) :: WTBLACC_M
    real, dimension(nlat,nmos,ignd) :: TBARACC_M
    real, dimension(nlat,nmos,ignd) :: THLQACC_M
    real, dimension(nlat,nmos,ignd) :: THICACC_M
    real, dimension(nlat,nmos,ignd) :: THALACC_M
    real, dimension(nlat,nmos) :: ALVSACC_M
    real, dimension(nlat,nmos) :: ALIRACC_M
    real, dimension(nlat,nmos) :: RHOSACC_M
    real, dimension(nlat,nmos) :: TSNOACC_M
    real, dimension(nlat,nmos) :: WSNOACC_M
    real, dimension(nlat,nmos) :: SNOARE_M
    real, dimension(nlat,nmos) :: TCANACC_M
    real, dimension(nlat,nmos) :: RCANACC_M
    real, dimension(nlat,nmos) :: SCANACC_M
    real, dimension(nlat,nmos) :: GROACC_M
    real, dimension(nlat,nmos) :: FSINACC_M
    real, dimension(nlat,nmos) :: FLINACC_M
    real, dimension(nlat,nmos) :: TAACC_M
    real, dimension(nlat,nmos) :: UVACC_M
    real, dimension(nlat,nmos) :: PRESACC_M
    real, dimension(nlat,nmos) :: QAACC_M
    real, dimension(nlat,nmos) :: EVAPACC_M
    real, dimension(nlat,nmos) :: FLUTACC_M

    real, dimension(nlat,nmos) :: tcanrs
    real, dimension(nlat,nmos) :: tsnors
    real, dimension(nlat,nmos) :: tpndrs
    real, dimension(nlat,nmos,ican) :: csum
    real, dimension(nlat,nmos,ignd) :: tbaraccrow_m
    real, dimension(nlat,nmos) :: tcanoaccrow_m
    real, dimension(nlat,nmos) :: uvaccrow_m
    real, dimension(nlat,nmos) :: vvaccrow_m

    real, dimension(nlat,nmos) :: tcanoaccrow_out
    real, dimension(nlat,nmos) :: qevpacc_m_save


end type veg_rot

type (veg_rot), save, target :: vrot

!=================================================================================
type veg_gat

    ! This is the basic data structure that contains the state variables
    ! for the Plant Functional Type (PFT). The dimensions are ilg,{icc,iccp1}

    real, dimension(ilg,icc) :: ailcmin           !
    real, dimension(ilg,icc) :: ailcmax           !
    real, dimension(ilg,icc) :: dvdfcan           !
    real, dimension(ilg,icc) :: gleafmas          !
    real, dimension(ilg,icc) :: bleafmas          !
    real, dimension(ilg,icc) :: stemmass          !
    real, dimension(ilg,icc) :: rootmass          !
    real, dimension(ilg,icc) :: pstemmass         !
    real, dimension(ilg,icc) :: pgleafmass        !
    real, dimension(ilg,icc) :: fcancmx

    real, dimension(ilg) :: gavglai

    real, dimension(ilg) :: lightng
    real, dimension(ilg) :: tcanoaccgat_out

    real, dimension(ilg,ican) :: zolnc
    real, dimension(ilg,ican) :: ailc

    real, dimension(ilg,icc) :: ailcg
    real, dimension(ilg,icc) :: ailcgs
    real, dimension(ilg,icc) :: fcancs
    real, dimension(ilg,icc) :: fcanc

    real, dimension(ilg) :: co2conc
    real, dimension(ilg) :: ch4conc

    real, dimension(ilg,icc) :: co2i1cg
    real, dimension(ilg,icc) :: co2i1cs
    real, dimension(ilg,icc) :: co2i2cg
    real, dimension(ilg,icc) :: co2i2cs
    real, dimension(ilg,icc) :: ancsveg
    real, dimension(ilg,icc) :: ancgveg
    real, dimension(ilg,icc) :: rmlcsveg
    real, dimension(ilg,icc) :: rmlcgveg
    real, dimension(ilg,icc) :: slai
    real, dimension(ilg,icc) :: ailcb
    real, dimension(ilg) :: canres
    real, dimension(ilg,icc) :: flhrloss

    real, dimension(ilg,icc) :: grwtheff
    real, dimension(ilg,icc) :: lystmmas
    real, dimension(ilg,icc) :: lyrotmas
    real, dimension(ilg,icc) :: tymaxlai
    real, dimension(ilg) :: vgbiomas
    real, dimension(ilg) :: gavgltms
    real, dimension(ilg) :: gavgscms
    real, dimension(ilg,icc) :: stmhrlos
    real, dimension(ilg,ican,ignd) :: rmatc
    real, dimension(ilg,icc,ignd) :: rmatctem
    real, dimension(ilg,iccp1) :: litrmass
    real, dimension(ilg,iccp1) :: soilcmas
    real, dimension(ilg,icc) :: vgbiomas_veg

    real, dimension(ilg,icc) :: emit_co2
    real, dimension(ilg,icc) :: emit_co
    real, dimension(ilg,icc) :: emit_ch4
    real, dimension(ilg,icc) :: emit_nmhc
    real, dimension(ilg,icc) :: emit_h2
    real, dimension(ilg,icc) :: emit_nox
    real, dimension(ilg,icc) :: emit_n2o
    real, dimension(ilg,icc) :: emit_pm25
    real, dimension(ilg,icc) :: emit_tpm
    real, dimension(ilg,icc) :: emit_tc
    real, dimension(ilg,icc) :: emit_oc
    real, dimension(ilg,icc) :: emit_bc
    real, dimension(ilg) :: burnfrac
    real, dimension(ilg,icc) :: burnvegf
    real, dimension(ilg) :: probfire
    real, dimension(ilg) :: bterm
    real, dimension(ilg) :: lterm
    real, dimension(ilg) :: mterm

    real, dimension(ilg) :: extnprob
    real, dimension(ilg) :: prbfrhuc
    real, dimension(ilg,12) :: mlightng

    real, dimension(ilg,icc) :: bmasveg
    real, dimension(ilg,ican) :: cmasvegc
    real, dimension(ilg,icc) :: veghght
    real, dimension(ilg,icc) :: rootdpth
    real, dimension(ilg) :: rml
    real, dimension(ilg) :: rms
    real, dimension(ilg,icc) :: tltrleaf
    real, dimension(ilg,icc) :: tltrstem
    real, dimension(ilg,icc) :: tltrroot
    real, dimension(ilg,icc) :: leaflitr
    real, dimension(ilg,icc) :: roottemp
    real, dimension(ilg,icc) :: afrleaf
    real, dimension(ilg,icc) :: afrstem
    real, dimension(ilg,icc) :: afrroot
    real, dimension(ilg,icc) :: wtstatus
    real, dimension(ilg,icc) :: ltstatus
    real, dimension(ilg) :: rmr

    real, dimension(ilg,8) :: wetfrac_s
    real, dimension(ilg) :: ch4wet1
    real, dimension(ilg) :: ch4wet2
    real, dimension(ilg) :: wetfdyn
    real, dimension(ilg) :: ch4dyn1
    real, dimension(ilg) :: ch4dyn2
    real, dimension(ilg) :: ch4_soills

    real, dimension(ilg) :: lucemcom
    real, dimension(ilg) :: lucltrin
    real, dimension(ilg) :: lucsocin

    real, dimension(ilg) :: npp
    real, dimension(ilg) :: nep
    real, dimension(ilg) :: nbp
    real, dimension(ilg) :: gpp
    real, dimension(ilg) :: hetrores
    real, dimension(ilg) :: autores
    real, dimension(ilg) :: soilcresp
    real, dimension(ilg) :: rm
    real, dimension(ilg) :: rg
    real, dimension(ilg) :: litres
    real, dimension(ilg) :: socres
    real, dimension(ilg) :: dstcemls
    real, dimension(ilg) :: litrfall
    real, dimension(ilg) :: humiftrs

    real, dimension(ilg,icc) :: gppveg
    real, dimension(ilg,iccp1) :: nepveg
    real, dimension(ilg,iccp1) :: nbpveg
    real, dimension(ilg,icc) :: nppveg
    real, dimension(ilg,iccp1) :: hetroresveg
    real, dimension(ilg,icc) :: autoresveg
    real, dimension(ilg,iccp1) :: litresveg
    real, dimension(ilg,iccp1) :: soilcresveg
    real, dimension(ilg,icc) :: rmlvegacc
    real, dimension(ilg,icc) :: rmsveg
    real, dimension(ilg,icc) :: rmrveg
    real, dimension(ilg,icc) :: rgveg

    real, dimension(ilg,icc) :: rothrlos
    real, dimension(ilg,icc) :: pfcancmx
    real, dimension(ilg,icc) :: nfcancmx
    real, dimension(ilg,ican) :: alvsctm
    real, dimension(ilg,ican) :: paic
    real, dimension(ilg,ican) :: slaic
    real, dimension(ilg,ican) :: alirctm
    real, dimension(ilg) :: cfluxcg
    real, dimension(ilg) :: cfluxcs
    real, dimension(ilg) :: dstcemls3
    real, dimension(ilg,icc) :: anveg
    real, dimension(ilg,icc) :: rmlveg

    real, dimension(ilg) :: twarmm              ! temperature of the warmest month (c)
    real, dimension(ilg) :: tcoldm              ! temperature of the coldest month (c)
    real, dimension(ilg) :: gdd5                ! growing degree days above 5 c
    real, dimension(ilg) :: aridity             ! aridity index, ratio of potential evaporation to precipitation
    real, dimension(ilg) :: srplsmon            ! number of months in a year with surplus water i.e. precipitation more than potential evaporation
    real, dimension(ilg) :: defctmon            ! number of months in a year with water deficit i.e. precipitation less than potential evaporation
    real, dimension(ilg) :: anndefct            ! annual water deficit (mm)
    real, dimension(ilg) :: annsrpls            ! annual water surplus (mm)
    real, dimension(ilg) :: annpcp              ! annual precipitation (mm)
    real, dimension(ilg) :: dry_season_length   ! length of dry season (months)

    logical, dimension(ilg,icc) :: pftexist
    integer, dimension(ilg,2) :: colddays
    integer, dimension(ilg) :: icount
    integer, dimension(ilg,icc) :: lfstatus
    integer, dimension(ilg,icc) :: pandays
    integer, dimension(ilg) :: stdaln

end type veg_gat

type (veg_gat), save, target :: vgat
!=================================================================================

type class_moyr_output

!   MONTHLY OUTPUT FOR CLASS GRID-MEAN

    real, dimension(nlat) :: ALVSACC_MO
    real, dimension(nlat) :: ALIRACC_MO
    real, dimension(nlat) :: FLUTACC_MO
    real, dimension(nlat) :: FSINACC_MO
    real, dimension(nlat) :: FLINACC_MO
    real, dimension(nlat) :: HFSACC_MO
    real, dimension(nlat) :: QEVPACC_MO
    real, dimension(nlat) :: SNOACC_MO
    real, dimension(nlat) :: WSNOACC_MO
    real, dimension(nlat) :: ROFACC_MO
    real, dimension(nlat) :: PREACC_MO
    real, dimension(nlat) :: EVAPACC_MO
    real, dimension(nlat) :: TAACC_MO

    real :: FSSTAR_MO
    real :: FLSTAR_MO
    real :: QH_MO
    real :: QE_MO

    real, dimension(nlat,ignd) :: TBARACC_MO
    real, dimension(nlat,ignd) :: THLQACC_MO
    real, dimension(nlat,ignd) :: THICACC_MO

!   YEARLY OUTPUT FOR CLASS GRID-MEAN

    real, dimension(nlat) :: ALVSACC_YR
    real, dimension(nlat) :: ALIRACC_YR
    real, dimension(nlat) :: FLUTACC_YR
    real, dimension(nlat) :: FSINACC_YR
    real, dimension(nlat) :: FLINACC_YR
    real, dimension(nlat) :: HFSACC_YR
    real, dimension(nlat) :: QEVPACC_YR
    real, dimension(nlat) :: ROFACC_YR
    real, dimension(nlat) :: PREACC_YR
    real, dimension(nlat) :: EVAPACC_YR
    real, dimension(nlat) :: TAACC_YR

    real :: FSSTAR_YR
    real :: FLSTAR_YR
    real :: QH_YR
    real :: QE_YR

end type class_moyr_output

type (class_moyr_output), save, target :: class_out

!=================================================================================

type ctem_gridavg

! Grid-averaged variables (denoted with an ending of "_g")

      real, dimension(nlat) :: WSNOROT_g
      real, dimension(nlat) :: ROFSROT_g
      real, dimension(nlat) :: SNOROT_g
      real, dimension(nlat) :: RHOSROT_g
      real, dimension(nlat) :: ROFROT_g
      real, dimension(nlat) :: ZPNDROT_g
      real, dimension(nlat) :: RCANROT_g
      real, dimension(nlat) :: SCANROT_g
      real, dimension(nlat) :: TROFROT_g
      real, dimension(nlat) :: TROOROT_g
      real, dimension(nlat) :: TROBROT_g
      real, dimension(nlat) :: ROFOROT_g
      real, dimension(nlat) :: ROFBROT_g
      real, dimension(nlat) :: TROSROT_g
      real, dimension(nlat) :: FSGVROT_g
      real, dimension(nlat) :: FSGSROT_g
      real, dimension(nlat) :: FLGVROT_g
      real, dimension(nlat) :: FLGSROT_g
      real, dimension(nlat) :: HFSCROT_g
      real, dimension(nlat) :: HFSSROT_g
      real, dimension(nlat) :: HEVCROT_g
      real, dimension(nlat) :: HEVSROT_g
      real, dimension(nlat) :: HMFCROT_g
      real, dimension(nlat) :: HMFNROT_g
      real, dimension(nlat) :: HTCSROT_g
      real, dimension(nlat) :: HTCCROT_g
      real, dimension(nlat) :: FSGGROT_g
      real, dimension(nlat) :: FLGGROT_g
      real, dimension(nlat) :: HFSGROT_g
      real, dimension(nlat) :: HEVGROT_g
      real, dimension(nlat) :: CDHROT_g
      real, dimension(nlat) :: CDMROT_g
      real, dimension(nlat) :: SFCUROT_g
      real, dimension(nlat) :: SFCVROT_g
      real, dimension(nlat) :: fc_g
      real, dimension(nlat) :: fg_g
      real, dimension(nlat) :: fcs_g
      real, dimension(nlat) :: fgs_g
      real, dimension(nlat) :: PCFCROT_g
      real, dimension(nlat) :: PCLCROT_g
      real, dimension(nlat) :: PCPGROT_g
      real, dimension(nlat) :: QFCFROT_g
      real, dimension(nlat) :: QFGROT_g
      real, dimension(nlat,ignd) :: QFCROT_g
      real, dimension(nlat) :: ROFCROT_g
      real, dimension(nlat) :: ROFNROT_g
      real, dimension(nlat) :: WTRSROT_g
      real, dimension(nlat) :: WTRGROT_g
      real, dimension(nlat) :: PCPNROT_g
      real, dimension(nlat) :: QFCLROT_g
      real, dimension(nlat) :: QFNROT_g
      real, dimension(nlat) :: WTRCROT_g
      real, dimension(nlat) :: gpp_g
      real, dimension(nlat) :: npp_g
      real, dimension(nlat) :: nbp_g
      real, dimension(nlat) :: socres_g
      real, dimension(nlat) :: autores_g
      real, dimension(nlat) :: litres_g
      real, dimension(nlat) :: dstcemls3_g
      real, dimension(nlat) :: litrfall_g
      real, dimension(nlat) :: rml_g
      real, dimension(nlat) :: rms_g
      real, dimension(nlat) :: rg_g
      real, dimension(nlat) :: leaflitr_g
      real, dimension(nlat) :: tltrstem_g
      real, dimension(nlat) :: tltrroot_g
      real, dimension(nlat) :: nep_g
      real, dimension(nlat) :: hetrores_g
      real, dimension(nlat) :: dstcemls_g
      real, dimension(nlat) :: humiftrs_g
      real, dimension(nlat) :: rmr_g
      real, dimension(nlat) :: tltrleaf_g
      real, dimension(nlat) :: gavgltms_g

      real, dimension(nlat) :: vgbiomas_g
      real, dimension(nlat) :: gavglai_g
      real, dimension(nlat) :: gavgscms_g
      real, dimension(nlat) :: gleafmas_g
      real, dimension(nlat) :: bleafmas_g
      real, dimension(nlat) :: stemmass_g
      real, dimension(nlat) :: rootmass_g
      real, dimension(nlat) :: litrmass_g
      real, dimension(nlat) :: soilcmas_g
      real, dimension(nlat) :: slai_g
      real, dimension(nlat) :: ailcg_g
      real, dimension(nlat) :: ailcb_g
      real, dimension(nlat) :: veghght_g
      real, dimension(nlat) :: rootdpth_g
      real, dimension(nlat) :: roottemp_g
      real, dimension(nlat) :: totcmass_g
      real, dimension(nlat) :: tcanoacc_out_g
      real, dimension(nlat) :: burnfrac_g
      real, dimension(nlat) :: probfire_g
      real, dimension(nlat) :: lucemcom_g
      real, dimension(nlat) :: lucltrin_g
      real, dimension(nlat) :: lucsocin_g
      real, dimension(nlat) :: emit_co2_g
      real, dimension(nlat) :: emit_co_g
      real, dimension(nlat) :: emit_ch4_g
      real, dimension(nlat) :: emit_nmhc_g
      real, dimension(nlat) :: emit_h2_g
      real, dimension(nlat) :: emit_nox_g
      real, dimension(nlat) :: emit_n2o_g
      real, dimension(nlat) :: emit_pm25_g
      real, dimension(nlat) :: emit_tpm_g
      real, dimension(nlat) :: emit_tc_g
      real, dimension(nlat) :: emit_oc_g
      real, dimension(nlat) :: emit_bc_g
      real, dimension(nlat) :: bterm_g
      real, dimension(nlat) :: lterm_g
      real, dimension(nlat) :: mterm_g
      real, dimension(nlat) :: ch4wet1_g
      real, dimension(nlat) :: ch4wet2_g
      real, dimension(nlat) :: wetfdyn_g
      real, dimension(nlat) :: ch4dyn1_g
      real, dimension(nlat) :: ch4dyn2_g
      real, dimension(nlat) :: ch4_soills_g
      real, dimension(nlat,icc) :: afrleaf_g
      real, dimension(nlat,icc) :: afrstem_g
      real, dimension(nlat,icc) :: afrroot_g
      real, dimension(nlat,icc) :: lfstatus_g
      real, dimension(nlat,icc) :: rmlvegrow_g
      real, dimension(nlat,icc) :: anvegrow_g
      real, dimension(nlat,ignd) :: rmatctem_g
      real, dimension(nlat,ignd) :: HMFGROT_g
      real, dimension(nlat,ignd) :: HTCROT_g
      real, dimension(nlat,ignd) :: TBARROT_g
      real, dimension(nlat,ignd) :: THLQROT_g
      real, dimension(nlat,ignd) :: THICROT_g
      real, dimension(nlat,ignd) :: GFLXROT_g

      real :: fsstar_g
      real :: flstar_g
      real :: qh_g
      real :: qe_g
      real :: snomlt_g
      real :: beg_g
      real :: gtout_g
      real :: tpn_g
      real :: altot_g
      real :: tcn_g
      real :: tsn_g
      real :: zsn_g

      ! Variables that are the same for the entire gridcell
      real, dimension(nlat) :: dayl_max ! maximum daylength for that location (hours)
      real, dimension(nlat) :: dayl     ! daylength for that location (hours)

end type ctem_gridavg

type (ctem_gridavg), save, target :: ctem_grd

!=================================================================================

type ctem_mosaic_level

!   Mosaic-level variables (denoted by an ending of "_m")

      real, dimension(nlat,nmos) :: leaflitr_m
      real, dimension(nlat,nmos) :: tltrleaf_m
      real, dimension(nlat,nmos) :: tltrstem_m
      real, dimension(nlat,nmos) :: tltrroot_m
      real, dimension(nlat,nmos) :: ailcg_m
      real, dimension(nlat,nmos) :: ailcb_m
      real, dimension(nlat,nmos,ignd) :: rmatctem_m
      real, dimension(nlat,nmos) :: veghght_m
      real, dimension(nlat,nmos) :: rootdpth_m
      real, dimension(nlat,nmos) :: roottemp_m
      real, dimension(nlat,nmos) :: slai_m
      real, dimension(nlat,nmos) :: afrroot_m
      real, dimension(nlat,nmos) :: afrleaf_m
      real, dimension(nlat,nmos) :: afrstem_m
      real, dimension(nlat,nmos) :: laimaxg_m
      real, dimension(nlat,nmos) :: stemmass_m
      real, dimension(nlat,nmos) :: rootmass_m
      real, dimension(nlat,nmos) :: litrmass_m
      real, dimension(nlat,nmos) :: gleafmas_m
      real, dimension(nlat,nmos) :: bleafmas_m
      real, dimension(nlat,nmos) :: soilcmas_m

      real, dimension(ilg) :: fsnowacc_m
      real, dimension(ilg) :: tcansacc_m
      real, dimension(ilg) :: tcanoaccgat_m
      real, dimension(ilg) :: taaccgat_m
      real, dimension(ilg) :: uvaccgat_m
      real, dimension(ilg) :: vvaccgat_m
      real, dimension(ilg,ignd) :: tbaraccgat_m
      real, dimension(ilg,ignd) :: tbarcacc_m
      real, dimension(ilg,ignd) :: tbarcsacc_m
      real, dimension(ilg,ignd) :: tbargacc_m
      real, dimension(ilg,ignd) :: tbargsacc_m
      real, dimension(ilg,ignd) :: thliqcacc_m
      real, dimension(ilg,ignd) :: thliqgacc_m
      real, dimension(ilg,ignd) :: thliqacc_m
      real, dimension(ilg,ignd) :: thicecacc_m
      real, dimension(ilg,icc) :: ancsvgac_m
      real, dimension(ilg,icc) :: ancgvgac_m
      real, dimension(ilg,icc) :: rmlcsvga_m
      real, dimension(ilg,icc) :: rmlcgvga_m
      integer, dimension(nlat,nmos) :: ifcancmx_m

end type ctem_mosaic_level

type (ctem_mosaic_level), save, target :: ctem_tile

!=================================================================================

type ctem_gridavg_monthly

!  Grid averaged monthly variables (denoted by name ending in "_mo_g")

    real, dimension(nlat) :: laimaxg_mo_g
    real, dimension(nlat) :: stemmass_mo_g
    real, dimension(nlat) :: rootmass_mo_g
    real, dimension(nlat) :: litrmass_mo_g
    real, dimension(nlat) :: soilcmas_mo_g
    real, dimension(nlat) :: npp_mo_g
    real, dimension(nlat) :: gpp_mo_g
    real, dimension(nlat) :: nep_mo_g
    real, dimension(nlat) :: nbp_mo_g
    real, dimension(nlat) :: hetrores_mo_g
    real, dimension(nlat) :: autores_mo_g
    real, dimension(nlat) :: litres_mo_g
    real, dimension(nlat) :: soilcres_mo_g
    real, dimension(nlat) :: vgbiomas_mo_g
    real, dimension(nlat) :: totcmass_mo_g
    real, dimension(nlat) :: emit_co2_mo_g
    real, dimension(nlat) :: emit_co_mo_g
    real, dimension(nlat) :: emit_ch4_mo_g
    real, dimension(nlat) :: emit_nmhc_mo_g
    real, dimension(nlat) :: emit_h2_mo_g
    real, dimension(nlat) :: emit_nox_mo_g
    real, dimension(nlat) :: emit_n2o_mo_g
    real, dimension(nlat) :: emit_pm25_mo_g
    real, dimension(nlat) :: emit_tpm_mo_g
    real, dimension(nlat) :: emit_tc_mo_g
    real, dimension(nlat) :: emit_oc_mo_g
    real, dimension(nlat) :: emit_bc_mo_g
    real, dimension(nlat) :: probfire_mo_g
    real, dimension(nlat) :: luc_emc_mo_g
    real, dimension(nlat) :: lucltrin_mo_g
    real, dimension(nlat) :: lucsocin_mo_g
    real, dimension(nlat) :: burnfrac_mo_g
    real, dimension(nlat) :: bterm_mo_g
    real, dimension(nlat) :: lterm_mo_g
    real, dimension(nlat) :: mterm_mo_g
    real, dimension(nlat) :: ch4wet1_mo_g
    real, dimension(nlat) :: ch4wet2_mo_g
    real, dimension(nlat) :: wetfdyn_mo_g
    real, dimension(nlat) :: ch4dyn1_mo_g
    real, dimension(nlat) :: ch4dyn2_mo_g
    real, dimension(nlat) :: ch4soills_mo_g

end type ctem_gridavg_monthly

type (ctem_gridavg_monthly), save, target :: ctem_grd_mo

!=================================================================================

type ctem_mosaic_monthly

!     Mosaic monthly variables (denoted by name ending in "_mo_m")

      real, dimension(nlat,nmos,icc) :: laimaxg_mo_m
      real, dimension(nlat,nmos,icc) :: stemmass_mo_m
      real, dimension(nlat,nmos,icc) :: rootmass_mo_m
      real, dimension(nlat,nmos,icc) :: npp_mo_m
      real, dimension(nlat,nmos,icc) :: gpp_mo_m
      real, dimension(nlat,nmos,icc) :: vgbiomas_mo_m
      real, dimension(nlat,nmos,icc) :: autores_mo_m
      real, dimension(nlat,nmos,icc) :: totcmass_mo_m
      real, dimension(nlat,nmos,iccp1) :: litrmass_mo_m
      real, dimension(nlat,nmos,iccp1) :: soilcmas_mo_m
      real, dimension(nlat,nmos,iccp1) :: nep_mo_m
      real, dimension(nlat,nmos,iccp1) :: litres_mo_m
      real, dimension(nlat,nmos,iccp1) :: soilcres_mo_m
      real, dimension(nlat,nmos,iccp1) :: hetrores_mo_m
      real, dimension(nlat,nmos,iccp1) :: nbp_mo_m
      real, dimension(nlat,nmos,icc) :: emit_co2_mo_m
      real, dimension(nlat,nmos,icc) :: emit_co_mo_m
      real, dimension(nlat,nmos,icc) :: emit_ch4_mo_m
      real, dimension(nlat,nmos,icc) :: emit_nmhc_mo_m
      real, dimension(nlat,nmos,icc) :: emit_h2_mo_m
      real, dimension(nlat,nmos,icc) :: emit_nox_mo_m
      real, dimension(nlat,nmos,icc) :: emit_n2o_mo_m
      real, dimension(nlat,nmos,icc) :: emit_pm25_mo_m
      real, dimension(nlat,nmos,icc) :: emit_tpm_mo_m
      real, dimension(nlat,nmos,icc) :: emit_tc_mo_m
      real, dimension(nlat,nmos,icc) :: emit_oc_mo_m
      real, dimension(nlat,nmos,icc) :: emit_bc_mo_m
      real, dimension(nlat,nmos,icc) :: burnfrac_mo_m
      real, dimension(nlat,nmos) :: probfire_mo_m
      real, dimension(nlat,nmos) :: bterm_mo_m
      real, dimension(nlat,nmos) :: luc_emc_mo_m
      real, dimension(nlat,nmos) :: lterm_mo_m
      real, dimension(nlat,nmos) :: lucsocin_mo_m
      real, dimension(nlat,nmos) :: mterm_mo_m
      real, dimension(nlat,nmos) :: lucltrin_mo_m
      real, dimension(nlat,nmos) :: ch4wet1_mo_m
      real, dimension(nlat,nmos) :: ch4wet2_mo_m
      real, dimension(nlat,nmos) :: wetfdyn_mo_m
      real, dimension(nlat,nmos) :: ch4dyn1_mo_m
      real, dimension(nlat,nmos) :: ch4dyn2_mo_m
      real, dimension(nlat,nmos) :: ch4soills_mo_m

end type ctem_mosaic_monthly

type (ctem_mosaic_monthly), save, target :: ctem_tile_mo


!=================================================================================

type ctem_gridavg_annual

! Annual output for CTEM grid-averaged variables:
! (denoted by name ending in "_yr_g")

    real, dimension(nlat) :: laimaxg_yr_g
    real, dimension(nlat) :: stemmass_yr_g
    real, dimension(nlat) :: rootmass_yr_g
    real, dimension(nlat) :: litrmass_yr_g
    real, dimension(nlat) :: soilcmas_yr_g
    real, dimension(nlat) :: npp_yr_g
    real, dimension(nlat) :: gpp_yr_g
    real, dimension(nlat) :: nep_yr_g
    real, dimension(nlat) :: nbp_yr_g
    real, dimension(nlat) :: hetrores_yr_g
    real, dimension(nlat) :: autores_yr_g
    real, dimension(nlat) :: litres_yr_g
    real, dimension(nlat) :: soilcres_yr_g
    real, dimension(nlat) :: vgbiomas_yr_g
    real, dimension(nlat) :: totcmass_yr_g
    real, dimension(nlat) :: emit_co2_yr_g
    real, dimension(nlat) :: emit_co_yr_g
    real, dimension(nlat) :: emit_ch4_yr_g
    real, dimension(nlat) :: emit_nmhc_yr_g
    real, dimension(nlat) :: emit_h2_yr_g
    real, dimension(nlat) :: emit_nox_yr_g
    real, dimension(nlat) :: emit_n2o_yr_g
    real, dimension(nlat) :: emit_pm25_yr_g
    real, dimension(nlat) :: emit_tpm_yr_g
    real, dimension(nlat) :: emit_tc_yr_g
    real, dimension(nlat) :: emit_oc_yr_g
    real, dimension(nlat) :: emit_bc_yr_g
    real, dimension(nlat) :: probfire_yr_g
    real, dimension(nlat) :: luc_emc_yr_g
    real, dimension(nlat) :: lucltrin_yr_g
    real, dimension(nlat) :: lucsocin_yr_g
    real, dimension(nlat) :: burnfrac_yr_g
    real, dimension(nlat) :: bterm_yr_g
    real, dimension(nlat) :: lterm_yr_g
    real, dimension(nlat) :: mterm_yr_g
    real, dimension(nlat) :: ch4wet1_yr_g
    real, dimension(nlat) :: ch4wet2_yr_g
    real, dimension(nlat) :: wetfdyn_yr_g
    real, dimension(nlat) :: ch4dyn1_yr_g
    real, dimension(nlat) :: ch4dyn2_yr_g
    real, dimension(nlat) :: ch4soills_yr_g

end type ctem_gridavg_annual

type (ctem_gridavg_annual), save, target :: ctem_grd_yr

!=================================================================================

type ctem_mosaic_annual

! c      Annual output for CTEM mosaic variables:
! c      (denoted by name ending in "_yr_m")
!
      real, dimension(nlat,nmos,icc) :: laimaxg_yr_m
      real, dimension(nlat,nmos,icc) :: stemmass_yr_m
      real, dimension(nlat,nmos,icc) :: rootmass_yr_m
      real, dimension(nlat,nmos,icc) :: npp_yr_m
      real, dimension(nlat,nmos,icc) :: gpp_yr_m
      real, dimension(nlat,nmos,icc) :: vgbiomas_yr_m
      real, dimension(nlat,nmos,icc) :: autores_yr_m
      real, dimension(nlat,nmos,icc) :: totcmass_yr_m
      real, dimension(nlat,nmos,iccp1) :: litrmass_yr_m
      real, dimension(nlat,nmos,iccp1) :: soilcmas_yr_m
      real, dimension(nlat,nmos,iccp1) :: nep_yr_m
      real, dimension(nlat,nmos,iccp1) :: litres_yr_m
      real, dimension(nlat,nmos,iccp1) :: soilcres_yr_m
      real, dimension(nlat,nmos,iccp1) :: hetrores_yr_m
      real, dimension(nlat,nmos,iccp1) :: nbp_yr_m
      real, dimension(nlat,nmos,icc) :: emit_co2_yr_m
      real, dimension(nlat,nmos,icc) :: emit_co_yr_m
      real, dimension(nlat,nmos,icc) :: emit_ch4_yr_m
      real, dimension(nlat,nmos,icc) :: emit_nmhc_yr_m
      real, dimension(nlat,nmos,icc) :: emit_h2_yr_m
      real, dimension(nlat,nmos,icc) :: emit_nox_yr_m
      real, dimension(nlat,nmos,icc) :: emit_n2o_yr_m
      real, dimension(nlat,nmos,icc) :: emit_pm25_yr_m
      real, dimension(nlat,nmos,icc) :: emit_tpm_yr_m
      real, dimension(nlat,nmos,icc) :: emit_tc_yr_m
      real, dimension(nlat,nmos,icc) :: emit_oc_yr_m
      real, dimension(nlat,nmos,icc) :: emit_bc_yr_m
      real, dimension(nlat,nmos,icc) :: burnfrac_yr_m
      real, dimension(nlat,nmos) :: probfire_yr_m
      real, dimension(nlat,nmos) :: bterm_yr_m
      real, dimension(nlat,nmos) :: luc_emc_yr_m
      real, dimension(nlat,nmos) :: lterm_yr_m
      real, dimension(nlat,nmos) :: lucsocin_yr_m
      real, dimension(nlat,nmos) :: mterm_yr_m
      real, dimension(nlat,nmos) :: lucltrin_yr_m
      real, dimension(nlat,nmos) :: ch4wet1_yr_m
      real, dimension(nlat,nmos) :: ch4wet2_yr_m
      real, dimension(nlat,nmos) :: wetfdyn_yr_m
      real, dimension(nlat,nmos) :: ch4dyn1_yr_m
      real, dimension(nlat,nmos) :: ch4dyn2_yr_m
      real, dimension(nlat,nmos) :: ch4soills_yr_m

end type ctem_mosaic_annual

type (ctem_mosaic_annual), save, target :: ctem_tile_yr


contains

!=================================================================================

subroutine initrowvars()

use ctem_params, only : nlat, nmos, ican, ignd ,icc, iccp1

implicit none

integer :: j,k,l,m

 do j = 1,nlat

   vrot%prbfrhuc(j)         = 0.0
   vrot%extnprob(j)         = 0.0
   !vrot%barf(j)                = 1.0

   do k =1,12
     vrot%mlightng(j,k) = 0.0
   end do

   do k = 1,nmos

!         vrot%PREACC_M(j,k) = 0.
!         vrot%GTACC_M(j,k) = 0.
!         vrot%QEVPACC_M(j,k) = 0.
!         vrot%HFSACC_M(j,k) = 0.
!         vrot%HMFNACC_M(j,k) = 0.
!         vrot%ROFACC_M(j,k) = 0.
!         vrot%SNOACC_M(j,k) = 0.
!         vrot%OVRACC_M(j,k) = 0.
!         vrot%WTBLACC_M(j,k) = 0.
!
!         vrot%ALVSACC_M(j,k) = 0.
!         vrot%ALIRACC_M(j,k) = 0.
!         vrot%RHOSACC_M(j,k) = 0.
!         vrot%TSNOACC_M(j,k) = 0.
!         vrot%WSNOACC_M(j,k) = 0.
!         vrot%TCANACC_M(j,k) = 0.
!         vrot%RCANACC_M(j,k) = 0.
!         vrot%SCANACC_M(j,k) = 0.
!         vrot%GROACC_M(j,k) = 0.
!         vrot%FSINACC_M(j,k) = 0.
!         vrot%FLINACC_M(j,k) = 0.
!         vrot%TAACC_M(j,k) = 0.
!         vrot%UVACC_M(j,k) = 0.
!         vrot%PRESACC_M(j,k) = 0.
!         vrot%QAACC_M(j,k) = 0.
!         vrot%EVAPACC_M(j,k) = 0.
!         vrot%FLUTACC_M(j,k) = 0.

        vrot%icount(j,k)           = 0
        vrot%co2conc(j,k)          = 0.0
        vrot%npp(j,k)              = 0.0
        vrot%nep(j,k)              = 0.0
        vrot%hetrores(j,k)         = 0.0
        vrot%autores(j,k)          = 0.0
        vrot%soilcresp(j,k)        = 0.0
        vrot%rm(j,k)               = 0.0
        vrot%rg(j,k)               = 0.0
        vrot%nbp(j,k)              = 0.0
        vrot%litres(j,k)           = 0.0
        vrot%socres(j,k)           = 0.0
        vrot%gpp(j,k)              = 0.0
        vrot%dstcemls(j,k)         = 0.0
        vrot%dstcemls3(j,k)        = 0.0
        vrot%litrfall(j,k)         = 0.0
        vrot%humiftrs(j,k)         = 0.0
        vrot%canres(j,k)           = 0.0
        vrot%rml(j,k)              = 0.0
        vrot%rms(j,k)              = 0.0
        vrot%rmr(j,k)              = 0.0
        vrot%lucemcom(j,k)         = 0.0
        vrot%lucltrin(j,k)         = 0.0
        vrot%lucsocin(j,k)         = 0.0
        vrot%burnfrac(j,k)         = 0.0
        vrot%probfire(j,k)         = 0.0
        vrot%bterm(j,k)            = 0.0
        vrot%lterm(j,k)            = 0.0
        vrot%mterm(j,k)            = 0.0
        vrot%cfluxcg(j,k)          = 0.0
        vrot%cfluxcs(j,k)          = 0.0
        !vrot%TCANOACC_M(j,k)       = 0.0
        !vrot%UVACC_M(j,k)          = 0.0
        !vrot%VVACC_M(j,k)          = 0.0
        !vrot%TCANOACC_OUT(j,k)     = 0.0
        vrot%ch4wet1(j,k)          = 0.0
        vrot%ch4wet2(j,k)          = 0.0
        vrot%wetfdyn(j,k)          = 0.0
        vrot%ch4dyn1(j,k)          = 0.0
        vrot%ch4dyn2(j,k)          = 0.0
        vrot%ch4_soills(j,k)       = 0.0


        do l=1,ignd
            vrot%tbaraccrow_m(j,k,l)  = 0.0
!             vrot%TBARACC_M(j,k,l) = 0.
!             vrot%THLQACC_M(j,k,l) = 0.
!             vrot%THICACC_M(j,k,l) = 0.
!             vrot%THALACC_M(j,k,l) = 0.
        end do

        do l=1,ican
            vrot%ZOLNC(j,k,l)        = 0.0
            vrot%AILC(j,k,l)         = 0.0
            vrot%CMASVEGC(j,k,l)     = 0.0
            vrot%ALVSCTM(j,k,l)      = 0.0
            vrot%ALIRCTM(j,k,l)      = 0.0
            vrot%CSUM(j,k,l)            = 0.0
            vrot%PAIC(j,k,l)         = 0.0
            vrot%SLAIC(j,k,l)        = 0.0

            do m = 1, 3
                vrot%RMATC(j,k,l,m)    = 0.0
            end do

        end do

        do l = 1,icc

            vrot%ailcmin(j,k,l) = 0.
            vrot%ailcmax(j,k,l) = 0.
            vrot%dvdfcan(j,k,l) = 0.
            vrot%gleafmas(j,k,l) = 0.
            vrot%bleafmas(j,k,l) = 0.
            vrot%stemmass(j,k,l) = 0.
            vrot%rootmass(j,k,l) = 0.
            vrot%pstemmass(j,k,l) = 0.
            vrot%pgleafmass(j,k,l) = 0.
            vrot%ailcg(j,k,l)        = 0.0
            vrot%ailcgs(j,k,l)       = 0.0
            vrot%fcancs(j,k,l)       = 0.0
            vrot%fcanc(j,k,l)        = 0.0
            vrot%fcancmx(j,k,l)      = 0.0
            vrot%co2i1cg(j,k,l)      = 0.0
            vrot%co2i1cs(j,k,l)      = 0.0
            vrot%co2i2cg(j,k,l)      = 0.0
            vrot%co2i2cs(j,k,l)      = 0.0
            vrot%ancsveg(j,k,l)      = 0.0
            vrot%ancgveg(j,k,l)      = 0.0
            vrot%rmlcsveg(j,k,l)     = 0.0
            vrot%rmlcgveg(j,k,l)     = 0.0
            vrot%stemmass(j,k,l)     = 0.0
            vrot%rootmass(j,k,l)     = 0.0
            vrot%ailcb(j,k,l)        = 0.0
            vrot%grwtheff(j,k,l)     = 0.0
            vrot%dvdfcan(j,k,l)      = 0.0
            vrot%bmasveg(j,k,l)      = 0.0
            vrot%tltrleaf(j,k,l)     = 0.0
            vrot%tltrstem(j,k,l)     = 0.0
            vrot%tltrroot(j,k,l)     = 0.0
            vrot%leaflitr(j,k,l)     = 0.0
            vrot%roottemp(j,k,l)     = 0.0
            vrot%afrleaf(j,k,l)      = 0.0
            vrot%afrstem(j,k,l)      = 0.0
            vrot%afrroot(j,k,l)      = 0.0
            vrot%wtstatus(j,k,l)     = 0.0
            vrot%ltstatus(j,k,l)     = 0.0
            vrot%ailcmin(j,k,l)      = 0.0
            vrot%ailcmax(j,k,l)      = 0.0
            vrot%pfcancmx(j,k,l)     = 0.0
            vrot%nfcancmx(j,k,l)     = 0.0
            vrot%nppveg(j,k,l)       = 0.0
            vrot%veghght(j,k,l)      = 0.0
            vrot%rootdpth(j,k,l)     = 0.0
            vrot%gleafmas(j,k,l)     = 0.0
            vrot%bleafmas(j,k,l)     = 0.0
            vrot%anveg(j,k,l)        = 0.0
            vrot%rmlveg(j,k,l)       = 0.0
            vrot%rmlvegacc(j,k,l)    = 0.0
            vrot%rmsveg(j,k,l)       = 0.0
            vrot%rmrveg(j,k,l)       = 0.0
            vrot%rgveg(j,k,l)        = 0.0
            vrot%vgbiomas_veg(j,k,l) = 0.0
            vrot%gppveg(j,k,l) = 0.0
            vrot%autoresveg(j,k,l) = 0.0
            vrot%emit_co2(j,k,l)         =0.0
            vrot%emit_co(j,k,l)          =0.0
            vrot%emit_ch4(j,k,l)         =0.0
            vrot%emit_nmhc(j,k,l)        =0.0
            vrot%emit_h2(j,k,l)          =0.0
            vrot%emit_nox(j,k,l)         =0.0
            vrot%emit_n2o(j,k,l)         =0.0
            vrot%emit_pm25(j,k,l)        =0.0
            vrot%emit_tpm(j,k,l)         =0.0
            vrot%emit_tc(j,k,l)          =0.0
            vrot%emit_oc(j,k,l)          =0.0
            vrot%emit_bc(j,k,l)          =0.0
            vrot%burnvegf(j,k,l)         =0.0

                do m = 1, ignd
                    vrot%rmatctem(j,k,l,m) = 0.0
                end do

            end do !icc
            !
            do l = 1, iccp1
                vrot%litrmass(j,k,l)    = 0.0
                vrot%soilcmas(j,k,l)    = 0.0
                vrot%hetroresveg(j,k,l) = 0.0
                vrot%litresveg(j,k,l) = 0.0
                vrot%soilcresveg(j,k,l) = 0.0
                vrot%nepveg(j,k,l) = 0.0
                vrot%nbpveg(j,k,l) = 0.0
            end do !iccp1

   end do !nmos
 end do !nlat

end subroutine initrowvars

!==================================================

subroutine resetclassmon(nltest)

use ctem_params, only : ignd

implicit none

integer, intent(in) :: nltest

integer :: i,j

do i=1,nltest
    class_out%ALVSACC_MO(I)=0.
    class_out%ALIRACC_MO(I)=0.
    class_out%FLUTACC_MO(I)=0.
    class_out%FSINACC_MO(I)=0.
    class_out%FLINACC_MO(I)=0.
    class_out%HFSACC_MO(I) =0.
    class_out%QEVPACC_MO(I)=0.
    class_out%SNOACC_MO(I) =0.
    class_out%WSNOACC_MO(I)=0.
    class_out%ROFACC_MO(I) =0.
    class_out%PREACC_MO(I) =0.
    class_out%EVAPACC_MO(I)=0.
    class_out%TAACC_MO(I)=0.

    DO J=1,IGND
        class_out%TBARACC_MO(I,J)=0.
        class_out%THLQACC_MO(I,J)=0.
        class_out%THICACC_MO(I,J)=0.
    end do
end do

end subroutine resetclassmon

!==================================================

subroutine resetclassyr(nltest)

implicit none

integer, intent(in) :: nltest

integer :: i

do i=1,nltest
          class_out%ALVSACC_YR(I)=0.
          class_out%ALIRACC_YR(I)=0.
          class_out%FLUTACC_YR(I)=0.
          class_out%FSINACC_YR(I)=0.
          class_out%FLINACC_YR(I)=0.
          class_out%HFSACC_YR(I) =0.
          class_out%QEVPACC_YR(I)=0.
          class_out%ROFACC_YR(I) =0.
          class_out%PREACC_YR(I) =0.
          class_out%EVAPACC_YR(I)=0.
          class_out%TAACC_YR(I)=0.
end do

end subroutine resetclassyr

!==================================================

subroutine resetmidmonth(nltest)

implicit none

integer, intent(in) :: nltest

integer :: i

do i=1,nltest
    ctem_grd_mo%stemmass_mo_g(i)=0.0
    ctem_grd_mo%rootmass_mo_g(i)=0.0
    ctem_grd_mo%litrmass_mo_g(i)=0.0
    ctem_grd_mo%soilcmas_mo_g(i)=0.0
    ctem_grd_mo%vgbiomas_mo_g(i)=0.0
    ctem_grd_mo%totcmass_mo_g(i)=0.0
end do

end subroutine resetmidmonth

!==================================================

subroutine resetmonthend_g(nltest)

implicit none

integer, intent(in) :: nltest

integer :: i

do i=1,nltest
    ctem_grd_mo%laimaxg_mo_g(i)=0.0
    ctem_grd_mo%npp_mo_g(i)=0.0
    ctem_grd_mo%gpp_mo_g(i)=0.0
    ctem_grd_mo%nep_mo_g(i)=0.0
    ctem_grd_mo%nbp_mo_g(i)=0.0
    ctem_grd_mo%hetrores_mo_g(i)=0.0
    ctem_grd_mo%autores_mo_g(i)=0.0
    ctem_grd_mo%litres_mo_g(i)=0.0
    ctem_grd_mo%soilcres_mo_g(i)=0.0
    ctem_grd_mo%emit_co2_mo_g(i)=0.0
    ctem_grd_mo%emit_co_mo_g(i) =0.0
    ctem_grd_mo%emit_ch4_mo_g(i) =0.0
    ctem_grd_mo%emit_nmhc_mo_g(i) =0.0
    ctem_grd_mo%emit_h2_mo_g(i) =0.0
    ctem_grd_mo%emit_nox_mo_g(i) =0.0
    ctem_grd_mo%emit_n2o_mo_g(i) =0.0
    ctem_grd_mo%emit_pm25_mo_g(i) =0.0
    ctem_grd_mo%emit_tpm_mo_g(i) =0.0
    ctem_grd_mo%emit_tc_mo_g(i) =0.0
    ctem_grd_mo%emit_oc_mo_g(i) =0.0
    ctem_grd_mo%emit_bc_mo_g(i) =0.0
    ctem_grd_mo%probfire_mo_g(i) =0.0
    ctem_grd_mo%luc_emc_mo_g(i) =0.0
    ctem_grd_mo%lucsocin_mo_g(i) =0.0
    ctem_grd_mo%lucltrin_mo_g(i) =0.0
    ctem_grd_mo%burnfrac_mo_g(i) =0.0
    ctem_grd_mo%bterm_mo_g(i)    =0.0
    ctem_grd_mo%lterm_mo_g(i)    =0.0
    ctem_grd_mo%mterm_mo_g(i)    =0.0
    ctem_grd_mo%ch4wet1_mo_g(i)  =0.0
    ctem_grd_mo%ch4wet2_mo_g(i)  =0.0
    ctem_grd_mo%wetfdyn_mo_g(i)  =0.0
    ctem_grd_mo%ch4dyn1_mo_g(i)  =0.0
    ctem_grd_mo%ch4dyn2_mo_g(i)  =0.0
    ctem_grd_mo%ch4soills_mo_g(i)  =0.0
end do

end subroutine resetmonthend_g

!==================================================

subroutine resetmonthend_m(nltest,nmtest)

use ctem_params, only : iccp1,icc

implicit none

integer, intent(in) :: nltest
integer, intent(in) :: nmtest

integer :: i,m,j

do i=1,nltest
 do m = 1,nmtest
   do j=1,icc
    ctem_tile_mo%npp_mo_m(i,m,j)=0.0
    ctem_tile_mo%gpp_mo_m(i,m,j)=0.0
    ctem_tile_mo%nep_mo_m(i,m,j)=0.0
    ctem_tile_mo%nbp_mo_m(i,m,j)=0.0
    ctem_tile_mo%laimaxg_mo_m(i,m,j)=0.0
    ctem_tile_mo%emit_co2_mo_m(i,m,j)=0.0
    ctem_tile_mo%emit_co_mo_m(i,m,j) =0.0
    ctem_tile_mo%emit_ch4_mo_m(i,m,j) =0.0
    ctem_tile_mo%emit_nmhc_mo_m(i,m,j) =0.0
    ctem_tile_mo%emit_h2_mo_m(i,m,j) =0.0
    ctem_tile_mo%emit_nox_mo_m(i,m,j) =0.0
    ctem_tile_mo%emit_n2o_mo_m(i,m,j) =0.0
    ctem_tile_mo%emit_pm25_mo_m(i,m,j) =0.0
    ctem_tile_mo%emit_tpm_mo_m(i,m,j) =0.0
    ctem_tile_mo%emit_tc_mo_m(i,m,j) =0.0
    ctem_tile_mo%emit_oc_mo_m(i,m,j) =0.0
    ctem_tile_mo%emit_bc_mo_m(i,m,j) =0.0
    ctem_tile_mo%burnfrac_mo_m(i,m,j) =0.0

   end do

    ctem_tile_mo%nep_mo_m(i,m,iccp1)=0.0
    ctem_tile_mo%nbp_mo_m(i,m,iccp1)=0.0

    ctem_tile_mo%probfire_mo_m(i,m) =0.0
    ctem_tile_mo%luc_emc_mo_m(i,m) =0.0
    ctem_tile_mo%lucsocin_mo_m(i,m) =0.0
    ctem_tile_mo%lucltrin_mo_m(i,m) =0.0
    ctem_tile_mo%bterm_mo_m(i,m)=0.0
    ctem_tile_mo%lterm_mo_m(i,m)=0.0
    ctem_tile_mo%mterm_mo_m(i,m)=0.0

    ctem_tile_mo%ch4wet1_mo_m(i,m)  =0.0
    ctem_tile_mo%ch4wet2_mo_m(i,m)  =0.0
    ctem_tile_mo%wetfdyn_mo_m(i,m)  =0.0
    ctem_tile_mo%ch4dyn1_mo_m(i,m)  =0.0
    ctem_tile_mo%ch4dyn2_mo_m(i,m)  =0.0
    ctem_tile_mo%ch4soills_mo_m(i,m)  =0.0

  end do !nmtest
end do ! nltest

end subroutine resetmonthend_m

!==================================================

subroutine resetyearend_g(nltest)

implicit none

integer, intent(in) :: nltest

integer :: i

do i=1,nltest

    ctem_grd_yr%laimaxg_yr_g(i)=0.0
    ctem_grd_yr%stemmass_yr_g(i)=0.0
    ctem_grd_yr%rootmass_yr_g(i)=0.0
    ctem_grd_yr%litrmass_yr_g(i)=0.0
    ctem_grd_yr%soilcmas_yr_g(i)=0.0
    ctem_grd_yr%vgbiomas_yr_g(i)=0.0
    ctem_grd_yr%totcmass_yr_g(i)=0.0
    ctem_grd_yr%npp_yr_g(i)=0.0
    ctem_grd_yr%gpp_yr_g(i)=0.0
    ctem_grd_yr%nep_yr_g(i)=0.0
    ctem_grd_yr%nbp_yr_g(i)=0.0
    ctem_grd_yr%hetrores_yr_g(i)=0.0
    ctem_grd_yr%autores_yr_g(i)=0.0
    ctem_grd_yr%litres_yr_g(i)=0.0
    ctem_grd_yr%soilcres_yr_g(i)=0.0
    ctem_grd_yr%emit_co2_yr_g(i)=0.0
    ctem_grd_yr%emit_co_yr_g(i)=0.0
    ctem_grd_yr%emit_ch4_yr_g(i)=0.0
    ctem_grd_yr%emit_nmhc_yr_g(i)=0.0
    ctem_grd_yr%emit_h2_yr_g(i)=0.0
    ctem_grd_yr%emit_nox_yr_g(i)=0.0
    ctem_grd_yr%emit_n2o_yr_g(i)=0.0
    ctem_grd_yr%emit_pm25_yr_g(i)=0.0
    ctem_grd_yr%emit_tpm_yr_g(i)=0.0
    ctem_grd_yr%emit_tc_yr_g(i)=0.0
    ctem_grd_yr%emit_oc_yr_g(i)=0.0
    ctem_grd_yr%emit_bc_yr_g(i)=0.0
    ctem_grd_yr%probfire_yr_g(i)=0.0
    ctem_grd_yr%luc_emc_yr_g(i)=0.0
    ctem_grd_yr%lucsocin_yr_g(i)=0.0
    ctem_grd_yr%lucltrin_yr_g(i)=0.0
    ctem_grd_yr%burnfrac_yr_g(i)=0.0
    ctem_grd_yr%bterm_yr_g(i)=0.0
    ctem_grd_yr%lterm_yr_g(i)=0.0
    ctem_grd_yr%mterm_yr_g(i)=0.0
    ctem_grd_yr%ch4wet1_yr_g(i)  =0.0
    ctem_grd_yr%ch4wet2_yr_g(i)  =0.0
    ctem_grd_yr%wetfdyn_yr_g(i)  =0.0
    ctem_grd_yr%ch4dyn1_yr_g(i)  =0.0
    ctem_grd_yr%ch4dyn2_yr_g(i)  =0.0
    ctem_grd_yr%ch4soills_yr_g(i)  =0.0

end do

end subroutine resetyearend_g

!==================================================

subroutine resetyearend_m(nltest,nmtest)

use ctem_params, only : iccp1,icc

implicit none

integer, intent(in) :: nltest
integer, intent(in) :: nmtest

integer :: i,m,j

do i=1,nltest
 do m = 1,nmtest
   do j=1,icc

    ctem_tile_yr%laimaxg_yr_m(i,m,j)=0.0
    ctem_tile_yr%npp_yr_m(i,m,j)=0.0
    ctem_tile_yr%gpp_yr_m(i,m,j)=0.0
    ctem_tile_yr%nep_yr_m(i,m,j)=0.0
    ctem_tile_yr%nbp_yr_m(i,m,j)=0.0
    ctem_tile_yr%hetrores_yr_m(i,m,j)=0.0
    ctem_tile_yr%autores_yr_m(i,m,j)=0.0
    ctem_tile_yr%litres_yr_m(i,m,j)=0.0
    ctem_tile_yr%soilcres_yr_m(i,m,j)=0.0

    ctem_tile_yr%emit_co2_yr_m(i,m,j)=0.0
    ctem_tile_yr%emit_co_yr_m(i,m,j)=0.0
    ctem_tile_yr%emit_ch4_yr_m(i,m,j)=0.0
    ctem_tile_yr%emit_nmhc_yr_m(i,m,j)=0.0
    ctem_tile_yr%emit_h2_yr_m(i,m,j)=0.0
    ctem_tile_yr%emit_nox_yr_m(i,m,j)=0.0
    ctem_tile_yr%emit_n2o_yr_m(i,m,j)=0.0
    ctem_tile_yr%emit_pm25_yr_m(i,m,j)=0.0
    ctem_tile_yr%emit_tpm_yr_m(i,m,j)=0.0
    ctem_tile_yr%emit_tc_yr_m(i,m,j)=0.0
    ctem_tile_yr%emit_oc_yr_m(i,m,j)=0.0
    ctem_tile_yr%emit_bc_yr_m(i,m,j)=0.0
    ctem_tile_yr%burnfrac_yr_m(i,m,j)=0.0

   end do

    ctem_tile_yr%hetrores_yr_m(i,m,iccp1)=0.0
    ctem_tile_yr%litres_yr_m(i,m,iccp1)=0.0
    ctem_tile_yr%soilcres_yr_m(i,m,iccp1)=0.0
    ctem_tile_yr%nep_yr_m(i,m,iccp1)=0.0
    ctem_tile_yr%nbp_yr_m(i,m,iccp1)=0.0

    ctem_tile_yr%probfire_yr_m(i,m)=0.0
    ctem_tile_yr%luc_emc_yr_m(i,m)=0.0
    ctem_tile_yr%lucsocin_yr_m(i,m)=0.0
    ctem_tile_yr%lucltrin_yr_m(i,m)=0.0
    ctem_tile_yr%bterm_yr_m(i,m)=0.0
    ctem_tile_yr%lterm_yr_m(i,m)=0.0
    ctem_tile_yr%mterm_yr_m(i,m)=0.0

    ctem_tile_yr%ch4wet1_yr_m(i,m)  =0.0
    ctem_tile_yr%ch4wet2_yr_m(i,m)  =0.0
    ctem_tile_yr%wetfdyn_yr_m(i,m)  =0.0
    ctem_tile_yr%ch4dyn1_yr_m(i,m)  =0.0
    ctem_tile_yr%ch4dyn2_yr_m(i,m)  =0.0
    ctem_tile_yr%ch4soills_yr_m(i,m)  =0.0

  end do !nmtest
end do ! nltest

end subroutine resetyearend_m

!==================================================
subroutine resetclassaccum(nltest,nmtest)

use ctem_params, only : ignd

implicit none

integer, intent(in) :: nltest
integer, intent(in) :: nmtest

integer :: i,m,j


DO I=1,NLTEST
  DO M=1,NMTEST

        vrot%PREACC_M(i,m) = 0.
        vrot%GTACC_M(i,m) = 0.
        vrot%QEVPACC_M(i,m) = 0.
        vrot%HFSACC_M(i,m) = 0.
        vrot%HMFNACC_M(i,m) = 0.
        vrot%ROFACC_M(i,m) = 0.
        vrot%SNOACC_M(i,m) = 0.
        vrot%OVRACC_M(i,m) = 0.
        vrot%WTBLACC_M(i,m) = 0.
        vrot%ALVSACC_M(i,m) = 0.
        vrot%ALIRACC_M(i,m) = 0.
        vrot%RHOSACC_M(i,m) = 0.
        vrot%TSNOACC_M(i,m) = 0.
        vrot%WSNOACC_M(i,m) = 0.
        vrot%SNOARE_M(i,m) = 0.
        vrot%TCANACC_M(i,m) = 0.
        vrot%RCANACC_M(i,m) = 0.
        vrot%SCANACC_M(i,m) = 0.
        vrot%GROACC_M(i,m) = 0.
        vrot%FSINACC_M(i,m) = 0.
        vrot%FLINACC_M(i,m) = 0.
        vrot%TAACC_M(i,m) = 0.
        vrot%UVACC_M(i,m) = 0.
        vrot%PRESACC_M(i,m) = 0.
        vrot%QAACC_M(i,m) = 0.
        vrot%EVAPACC_M(i,m) = 0.
        vrot%FLUTACC_M(i,m) = 0.

    DO J=1,IGND
        vrot%TBARACC_M(I,M,J)=0.
        vrot%THLQACC_M(I,M,J)=0.
        vrot%THICACC_M(I,M,J)=0.
        vrot%THALACC_M(I,M,J)=0.
    end do
  end do
end do

end subroutine resetclassaccum


!==================================================

subroutine resetgridavg(nltest)

use ctem_params, only : ignd,icc

implicit none

integer, intent(in) :: nltest

integer :: i,j

        do i = 1, nltest

        ctem_grd%fsstar_g    =0.0 !flag should some of these be per i that are presently not?   JM Jun 2015.
        ctem_grd%flstar_g    =0.0
        ctem_grd%qh_g        =0.0
        ctem_grd%qe_g        =0.0
        ctem_grd%snomlt_g    =0.0
        ctem_grd%beg_g       =0.0
        ctem_grd%gtout_g     =0.0
        ctem_grd%SNOROT_g(i) =0.0
        ctem_grd%RHOSROT_g(i)=0.0
        ctem_grd%WSNOROT_g(i)=0.0
        ctem_grd%altot_g     =0.0
        ctem_grd%ROFROT_g(i) =0.0
        ctem_grd%tpn_g       =0.0
        ctem_grd%ZPNDROT_g(i)=0.0

        do j=1,ignd
         ctem_grd%TBARROT_g(i,j)=0.0
         ctem_grd%THLQROT_g(i,j)=0.0
         ctem_grd%THICROT_g(i,j)=0.0
         ctem_grd%GFLXROT_g(i,j)=0.0
         ctem_grd%HMFGROT_g(i,j)=0.0
         ctem_grd%HTCROT_g(i,j)=0.0
         ctem_grd%QFCROT_g(i,j)=0.0
        end do

        ctem_grd%tcn_g=0.0
        ctem_grd%RCANROT_g(i) =0.0
        ctem_grd%SCANROT_g(i) =0.0
        ctem_grd%tsn_g=0.0
        ctem_grd%zsn_g=0.0
        ctem_grd%TROFROT_g(i)=0.0
        ctem_grd%TROOROT_g(i)=0.0
        ctem_grd%TROBROT_g(i)=0.0
        ctem_grd%TROSROT_g(i)=0.0
        ctem_grd%ROFOROT_g(i)=0.0
        ctem_grd%ROFSROT_g(i)=0.0
        ctem_grd%ROFBROT_g(i)=0.0
        ctem_grd%FSGVROT_g(i)=0.0
        ctem_grd%FSGSROT_g(i)=0.0
        ctem_grd%FSGGROT_g(i)=0.0
        ctem_grd%FLGVROT_g(i)=0.0
        ctem_grd%FLGSROT_g(i)=0.0
        ctem_grd%FLGGROT_g(i)=0.0
        ctem_grd%HFSCROT_g(i)=0.0
        ctem_grd%HFSSROT_g(i)=0.0
        ctem_grd%HFSGROT_g(i)=0.0
        ctem_grd%HEVCROT_g(i)=0.0
        ctem_grd%HEVSROT_g(i)=0.0
        ctem_grd%HEVGROT_g(i)=0.0
        ctem_grd%HMFCROT_g(i)=0.0
        ctem_grd%HMFNROT_g(i)=0.0
        ctem_grd%HTCCROT_g(i)=0.0
        ctem_grd%HTCSROT_g(i)=0.0
        ctem_grd%PCFCROT_g(i)=0.0
        ctem_grd%PCLCROT_g(i)=0.0
        ctem_grd%PCPNROT_g(i)=0.0
        ctem_grd%PCPGROT_g(i)=0.0
        ctem_grd%QFCFROT_g(i)=0.0
        ctem_grd%QFCLROT_g(i)=0.0
        ctem_grd%QFNROT_g(i)=0.0
        ctem_grd%QFGROT_g(i)=0.0
        ctem_grd%ROFCROT_g(i)=0.0
        ctem_grd%ROFNROT_g(i)=0.0
        ctem_grd%WTRCROT_g(i)=0.0
        ctem_grd%WTRSROT_g(i)=0.0
        ctem_grd%WTRGROT_g(i)=0.0
        ctem_grd%CDHROT_g(i)=0.0
        ctem_grd%CDMROT_g(i)=0.0
        ctem_grd%SFCUROT_g(i)=0.0
        ctem_grd%SFCVROT_g(i)=0.0

       if (c_switch%ctem_on) then
          do j=1,icc
            ctem_grd%anvegrow_g(i,j)=0.0
            ctem_grd%rmlvegrow_g(i,j)=0.0
          end do
       end if

       end do

end subroutine resetgridavg

!==================================================
subroutine resetctem_g(nltest)

use ctem_params, only : ignd,icc

implicit none

integer, intent(in) :: nltest

integer :: i,j,k

do i=1,nltest
    ctem_grd%gpp_g(i) =0.0
    ctem_grd%npp_g(i) =0.0
    ctem_grd%nep_g(i) =0.0
    ctem_grd%nbp_g(i) =0.0
    ctem_grd%autores_g(i) =0.0
    ctem_grd%hetrores_g(i)=0.0
    ctem_grd%litres_g(i) =0.0
    ctem_grd%socres_g(i) =0.0
    ctem_grd%dstcemls_g(i)=0.0
    ctem_grd%dstcemls3_g(i)=0.0
    ctem_grd%litrfall_g(i)=0.0
    ctem_grd%humiftrs_g(i)=0.0
    ctem_grd%rml_g(i) =0.0
    ctem_grd%rms_g(i) =0.0
    ctem_grd%rmr_g(i) =0.0
    ctem_grd%rg_g(i) =0.0
    ctem_grd%vgbiomas_g(i) =0.0
    ctem_grd%totcmass_g(i) =0.0
    ctem_grd%gavglai_g(i) =0.0
    ctem_grd%gavgltms_g(i) =0.0
    ctem_grd%gavgscms_g(i) =0.0
    ctem_grd%ailcg_g(i)=0.0
    ctem_grd%ailcb_g(i)=0.0
    ctem_grd%tcanoacc_out_g(i) =0.0
    ctem_grd%burnfrac_g(i) =0.0
    ctem_grd%probfire_g(i) =0.0
    ctem_grd%lucemcom_g(i) =0.0
    ctem_grd%lucltrin_g(i) =0.0
    ctem_grd%lucsocin_g(i) =0.0
    ctem_grd%emit_co2_g(i) =0.0
    ctem_grd%emit_co_g(i)  =0.0
    ctem_grd%emit_ch4_g(i) =0.0
    ctem_grd%emit_nmhc_g(i) =0.0
    ctem_grd%emit_h2_g(i) =0.0
    ctem_grd%emit_nox_g(i) =0.0
    ctem_grd%emit_n2o_g(i) =0.0
    ctem_grd%emit_pm25_g(i) =0.0
    ctem_grd%emit_tpm_g(i) =0.0
    ctem_grd%emit_tc_g(i) =0.0
    ctem_grd%emit_oc_g(i) =0.0
    ctem_grd%emit_bc_g(i) =0.0
    ctem_grd%bterm_g(i)   =0.0
    ctem_grd%lterm_g(i)   =0.0
    ctem_grd%mterm_g(i)   =0.0
    ctem_grd%leaflitr_g(i)=0.0
    ctem_grd%tltrleaf_g(i)=0.0
    ctem_grd%tltrstem_g(i)=0.0
    ctem_grd%tltrroot_g(i)=0.0
    ctem_grd%gleafmas_g(i)=0.0
    ctem_grd%bleafmas_g(i)=0.0
    ctem_grd%stemmass_g(i)=0.0
    ctem_grd%rootmass_g(i)=0.0
    ctem_grd%litrmass_g(i)=0.0
    ctem_grd%soilcmas_g(i)=0.0
    ctem_grd%veghght_g(i)=0.0
    ctem_grd%rootdpth_g(i)=0.0
    ctem_grd%roottemp_g(i)=0.0
    ctem_grd%slai_g(i)=0.0

    ctem_grd%CH4WET1_G(i) = 0.0
    ctem_grd%CH4WET2_G(i) = 0.0
    ctem_grd%WETFDYN_G(i) = 0.0
    ctem_grd%CH4DYN1_G(i) = 0.0
    ctem_grd%CH4DYN2_G(i) = 0.0
    ctem_grd%ch4_soills_g(i) = 0.0

    do k=1,ignd
      ctem_grd%rmatctem_g(i,k)=0.0
    enddo

    do j=1,icc
      ctem_grd%afrleaf_g(i,j)=0.0
      ctem_grd%afrstem_g(i,j)=0.0
      ctem_grd%afrroot_g(i,j)=0.0
    enddo
end do

end subroutine resetctem_g

!==================================================

subroutine finddaylength(solday,radl,daylength)

! Calculate the daylength based on the latitude and day of year

! Joe Melton Dec 18 2015 (taken from phenlogy.f)

use ctem_params, only : pi

implicit none

real, intent(in) :: solday  !day of year
real, intent(in) :: radl    ! latitude
real, intent(out) :: daylength  ! calculated daylength
real :: theta               ! temp var
real :: decli               ! temp var
real :: term                ! temp var

    theta=0.2163108 + 2.0*atan(0.9671396*tan(0.0086*(solday-186.0)))
    decli=asin(0.39795*cos(theta))    !declination !note I see that CLASS does this also but with different formula...
    term=(sin(radl)*sin(decli))/(cos(radl)*cos(decli))
    term=max(-1.0,min(term,1.0))
    daylength=24.0-(24.0/pi)*acos(term)

end subroutine finddaylength

!==================================================


! separate one:
!
! c     reset mosaic accumulator arrays.
! c
!       do 655 i=1,nml
!          uvaccgat_m(i)=0.0
! 655   continue
! c
!       if (ctem_on) then
!         do 705 i = 1, nml
! c
! c         competitition related variables added by y. peng \\
!           fsinacc_gat(i)=0.
!           flinacc_gat(i)=0.
!           flutacc_gat(i)=0.
!           alswacc_gat(i)=0.
!           allwacc_gat(i)=0.
!           pregacc_gat(i)=0.
! c         competitition related variables added by y. peng //
! c
!           fsnowacc_m(i)=0.0
!           tcanoaccgat_out(i)=tcanoaccgat_m(i)
!           tcanoaccgat_m(i)=0.0
! c
!           tcansacc_m(i)=0.0
!           taaccgat_m(i)=0.0
!           vvaccgat_m(i)=0.0
! c
!           do 715 j=1,ignd
!              tbaraccgat_m(i,j)=0.0
!              tbarcacc_m(i,j)=0.0
!              tbarcsacc_m(i,j)=0.0
!              tbargacc_m(i,j)=0.0
!              tbargsacc_m(i,j)=0.0
!              thliqcacc_m(i,j)=0.0
!              thliqgacc_m(i,j)=0.0
!              thicecacc_m(i,j)=0.0
! 715       continue
! c
!           do 716 j = 1, icc
!             ancsvgac_m(i,j)=0.0
!             ancgvgac_m(i,j)=0.0
!             rmlcsvga_m(i,j)=0.0
!             rmlcgvga_m(i,j)=0.0
! 716       continue
! c
! 705     continue
!       endif  ! if(ctem_on)
!       endif  ! if(ncount.eq.nday)
!=================================================================================

end module ctem_statevars
