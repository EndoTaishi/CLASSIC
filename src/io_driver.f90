!>Central module that handles all CTEM preparation and writing of output files
module io_driver

! J. Melton Mar 30 2015

implicit none

! subroutines contained in this module:

!public  :: class_hh_w           ! Accumulates and writes the CLASS half hourly file
!public  :: class_daily_aw     ! Accumulates and writes the CLASS daily file
public  :: class_monthly_aw     ! Accumulates and writes the CLASS monthly file
public  :: class_annual_aw      ! Accumulates and writes the CLASS annual file
public  :: ctem_daily_aw        ! Accumulates and writes the CTEM daily file
public  :: ctem_monthly_aw      ! Accumulates and writes the CTEM monthly file
public  :: ctem_annual_aw       ! Accumulates and writes the CTEM annual file


contains

        !     * CLASS daily and half-hourly output files (monthly and annual are done in io_driver)
!         !
!         if (.not. parallelrun) then ! stand alone mode, includes half-hourly and daily output
!             OPEN(UNIT=61,FILE='test.OF1_G')  ! GRID-LEVEL DAILY OUTPUT FROM CLASS
!             OPEN(UNIT=62,FILE='test.OF2_G')
!             OPEN(UNIT=63,FILE='test.OF3_G')
!
!             OPEN(UNIT=611,FILE='test.OF1_M') ! MOSAIC DAILY OUTPUT FROM CLASS
!             OPEN(UNIT=621,FILE='test.OF2_M')
!             OPEN(UNIT=631,FILE='test.OF3_M')
!
!             OPEN(UNIT=64,FILE='test.OF4_M')  ! MOSAIC HALF-HOURLY OUTPUT FROM CLASS
!             OPEN(UNIT=65,FILE='test.OF5_M')
!             OPEN(UNIT=66,FILE='test.OF6_M')
!             OPEN(UNIT=67,FILE='test.OF7_M')
!             OPEN(UNIT=68,FILE='test.OF8_M')
!             OPEN(UNIT=69,FILE='test.OF9_M')
!
!             OPEN(UNIT=641,FILE='test.OF4_G') ! GRID-LEVEL HALF-HOURLY OUTPUT FROM CLASS
!             OPEN(UNIT=651,FILE='test.OF5_G')
!             OPEN(UNIT=661,FILE='test.OF6_G')
!             OPEN(UNIT=671,FILE='test.OF7_G')
!             OPEN(UNIT=681,FILE='test.OF8_G')
!             OPEN(UNIT=691,FILE='test.OF9_G')
!         end if
!
! !    peatland output-FLAG!!---------------------------------------------\
!       open(unit=93,file='test.CT11D_G') !peatland GPP components
!       open(unit=94,file='test.CT12D_G') !peatland vegetation height & lai
!       open(unit=95,file='test.CT13D_G') !peatland C pools (in ctem.f)
!       open(unit=96,file='test.CT14D_G') !peatland soil resp (in ctem.f)
!       open(unit=97,file='test.CT15D_G') !peatland decompositon components (in decp.f)
!       open(unit=98,file='test.CT16D_G') !peatland moss photosynthesis midday sub-areas(mosspht.f)
!       open(unit=99,file='test.CT17D_G') !peatland water balance
!       open(unit=90,file='test.CT18Y_G') !peatland depth information

!    YW March 25, 2015 ----------------------------------------------/
        !
!         IF (.NOT. PARALLELRUN) THEN ! STAND ALONE MODE, INCLUDES HALF-HOURLY AND DAILY OUTPUT
!             !
!             WRITE(61,6001) 'test'
!             WRITE(61,6002) 'test'
!             WRITE(61,6011)
! 6011            FORMAT(2X,'DAY  YEAR  K*  L*  QH  QE  SM  QG  ',&
!                 &          'TR  SWE  DS  WS  AL  ROF  CUMS')
!
!             WRITE(62,6001) 'test'
!             WRITE(62,6002) 'test'
!
!             IF(IGND.GT.3) THEN
!                 WRITE(62,6012)
! 6012                FORMAT(2X,'DAY  YEAR  TG1  THL1  THI1  TG2  THL2  THI2  ',&
!                     &    'TG3  THL3  THI3  TG4  THL4  THI4  TG5  THL5  ',&
!                     &    'THI5 TG6  THL6  THI6 TG7  THL7  THI7',&
!                     &    'TG8  THL8  THI8 TG9  THL9  THI9 TG10  THL10  THI10',&
!                     &    'TG11  THL11  THI11 TG12  THL12  THI12 TG13  THL13  THI13',&
!                     &    'TG14  THL14  THI14 TG15  THL15  THI15 TG16  THL16  THI16',&
!                     &    'TG17  THL17  THI17 TG18  THL18  THI18 TG19  THL19  THI19',&
!                     &    'TG20  THL20  THI20 ACTLYR FTABLE')
!
!             ELSE
!                 WRITE(62,6212)
! 6212                FORMAT(2X,'DAY  YEAR  TG1  THL1  THI1  TG2  THL2  THI2  ',&
!                     &              'TG3  THL3  THI3  TCN  RCAN  SCAN  TSN  ZSN',&
!                     &              'ACTLYR FTABLE')
!
!             ENDIF
!
!             WRITE(63,6001) 'test'
!             WRITE(63,6002) 'test'
!             WRITE(63,6313)
! 6313            FORMAT(2X,'DAY YEAR KIN LIN TA UV PRES QA PCP EVAP')
!             !
!             WRITE(64,6001) 'test'
!             WRITE(64,6002) 'test'
!             WRITE(64,6014)
! 6014            FORMAT(2X,'HOUR  MIN  DAY  YEAR  K*  L*  QH  QE  SM  QG  ',&
!                 &          'TR  SWE  DS  WS  AL  ROF  TPN  ZPN  CDH  CDM  ',&
!                 &          'SFCU  SFCV  UV')
!
!             WRITE(65,6001) 'test'
!             WRITE(65,6002) 'test'
!
!             WRITE(65,6515)
! 6515            FORMAT(2X,'HOUR  MIN  DAY  YEAR  TG1  THL1  THI1  TG2  ',&
!                 &           'THL2  THI2  TG3  THL3  THI3  TCN  RCAN  SCAN  ',&
!                 &           'TSN  ZSN  TCN-TA  TCANO  TAC  ACTLYR  FTABLE')
!
!             WRITE(66,6001) 'test'
!             WRITE(66,6002) 'test'
!
!             IF(IGND.GT.3) THEN
!                 WRITE(66,6016)
! 6016                FORMAT(2X,'HOUR  MIN  DAY  YEAR  TG6  THL6  THI6  TG7  ',&
!                     &          'THL7  THI7  TG8  THL8  THI8  TG9  THL9  THI9  ',&
!                     &          'TG10  THL10  THI10  TG11  THL11  THI11  TG12  ',&
!                     &          'THL12  THI12  TG13  THL13  THI13  TG14  THL14  ',&
!                     &          'THI14  TG15  THL15  THI15  TG16  THL16  THI16  ',&
!                     &          'TG17  THL17  THI17  TG18  THL18  THI18  TG19  ' ,&
!                     &          'THL19  THI19  TG20  THL20  THI20  TG21  THL21  ',&
!                     &          'THI21  TG22  THL22 THI22  TG23  THL23  THI23  ',&
!                     &          'TG24  THL24  THI24  TG25  THL25  THI25  TG26  ',&
!                     &          'THL26  THI26  G0  G1  G2  G3  G4  G5  G6  ',&
!                     &          'G7  G8  G9')
!
!             ELSE
!                 WRITE(66,6616)
!                 WRITE(66,6615)
! 6616                FORMAT(2X,'HOUR  MIN  DAY  SWIN  LWIN  PCP  TA  VA  PA  QA')
! 6615                FORMAT(2X,'IF IGND <= 3, THIS FILE IS EMPTY')
!             ENDIF
!
!             WRITE(67,6001) 'test'
!             WRITE(67,6002) 'test'
!             WRITE(67,6017)
!             !     6017  FORMAT(2X,'WCAN SCAN CWLCAP CWFCAP FC FG FCS FGS CDH ', !runclass formatted.
!             !     1          'TCANO TCANS ALBS')
! 6017            FORMAT(2X,'HOUR  MIN  DAY  YEAR  ',&
!                 &  'TROF     TROO     TROS     TROB      ROF     ROFO   ',&
!                 &  '  ROFS        ROFB         FCS        FGS        FC       FG')
!
!             WRITE(68,6001) 'test'
!             WRITE(68,6002) 'test'
!             WRITE(68,6018)
! 6018            FORMAT(2X,'HOUR  MIN  DAY  YEAR  ',&
!                 &          'FSGV FSGS FSGG FLGV FLGS FLGG HFSC HFSS HFSG ',&
!                 &          'HEVC HEVS HEVG HMFC HMFS HMFG1 HMFG2 HMFG3 ',&
!                 &          'HTCC HTCS HTC1 HTC2 HTC3')
!
!             WRITE(69,6001) 'test'
!             WRITE(69,6002) 'test'
!             WRITE(69,6019)
! 6019            FORMAT(2X,'HOUR  MIN  DAY  YEAR  ',&
!                 &   'PCFC PCLC PCPN PCPG QFCF QFCL QFN QFG QFC1 ',&
!                 &          'QFC2 QFC3 ROFC ROFN ROFO ROF WTRC WTRS WTRG')
!             !       runclass also has: EVDF ','CTV CTS CT1 CT2 CT3')
!             !
!             WRITE(611,6001) 'test'
!             WRITE(611,6002) 'test'
!             WRITE(611,6011)
!             WRITE(621,6001) 'test'
!             WRITE(621,6002) 'test'
!             !
!             IF(IGND.GT.3) THEN
!                 WRITE(621,6012)
!             ELSE
!                 WRITE(621,6212)
!             ENDIF
!             !
!             WRITE(631,6001) 'test'
!             WRITE(631,6002) 'test'
!             WRITE(631,6313)
!             !
!             WRITE(641,6001) 'test'
!             WRITE(641,6002) 'test'
!             WRITE(641,6008)
!
!             WRITE(651,6001) 'test'
!             WRITE(651,6002) 'test'
!             WRITE(651,6515)
!             !
!             WRITE(661,6001) 'test'
!             WRITE(661,6002) 'test'
!             !
!             IF(IGND.GT.3) THEN
!                 WRITE(661,6016)
!             ELSE
!                 WRITE(661,6616)
!             ENDIF
!             !
!             WRITE(671,6001) 'test'
!             WRITE(671,6002) 'test'
!             WRITE(671,6017)
!             WRITE(681,6001) 'test'
!             WRITE(681,6002) 'test'
!             WRITE(681,6018)
!             WRITE(691,6001) 'test'
!             WRITE(691,6002) 'test'
!             WRITE(691,6019)
!             !
! 6008            FORMAT(2X,'HOUR  MIN  DAY  YEAR  K*  L*  QH  QE  SM  QG  ',&
!                 &          'TR  SWE  DS  WS  AL  ROF  TPN  ZPN  CDH  CDM  ',&
!                 &          'SFCU  SFCV  UV')
!
! !    --------write peatland output-------------------------------------\
!
!           write(93,6903)
!           write(94,6904)
!           write(95,6905)
!           write(96,6906)
!           write(97,6907)
!           write(98,6908)
!           write(99,6909)
!
! 6903  format (2X,'iday iyear nppmoss   armoss   gppmoss   ',&
!          'gppveg1  gppveg2  gppveg3  gppveg4   gppveg5   gppveg6  ',&
!          'gppveg7  gppveg8  gppveg9  gppveg10  gppveg11  gppveg12',&
!          'nppveg1  nppveg2  nppveg3  nppveg4   nppveg5   nppveg6  ',&
!          'nppveg7  nppveg8  nppveg9  nppveg10  nppveg11  nppveg12',&
!          'autoresp1   autoresp2   autoresp3   autoresp4  autoresp5  ',&
!          'autoresp6   autoresp7   autoresp8   autoresp9  autoresp10 ',&
!          'autoresp11  autoresp12  heteresp1   heteresp2  heteresp3  ',&
!          'heteresp4   heteresp5   heteresp6   heteresp7  heteresp8  ',&
!          'heteresp9   heteresp10  heteresp11  heteresp12   ',&
!          'fcancmx1  fcancmx2  fcancmx3  fcancmx4  fcancmx5  fcancmx6 ',&
!          'fcancmx7  fcancmx8  fcancmx9  fcancmx10 fcancmx11 fcancmx12')
! 6904  format (2X,'iday   iyear   veghght1   veghght2   veghght3   ',&
!          'veghght4   veghght5   veghght6   veghght7   veghght8   ',&
!          'veghght9   veghght10  veghght11  veghght12  rootdpt1   ',&
!          'rootdpt2   rootdpt3   rootdpt4   rootdpt5   rootdpt6   ',&
!          'rootdpt7   rootdpt8   rootdpt9   rootdpt10  rootdpt11  ',&
!          'rootdpt12  ailcg1   ailcg2   ailcg3   ailcg4   ailcg5  ',&
!          'ailcg6  ailcg7   ailcg8   ailcg9    ailcg10   ailcg11   ',&
!          'ailcg12     stemmas1    stemmas2   stemmas3   stemmas4   ',&
!          'stemmas5    stemmas6    stemmas7   stemmas8   stemmas9   ',&
!          'stemmas10   stemmas11   stemmas12  rootmas1   rootmas2   ',&
!          'rootmas3   rootmas4   rootmas5   rootmas6    rootmas7    ',&
!          'rootmas8   rootmas9   rootmas10   rootmas11  rootmas12   ',&
!          'litrmas1   litrmas2   litrmas3   litrmas4    litrmas5    ',&
!          'litrmas6   litrmas7   litrmas8   litrmas9    litrmas10   ',&
!          'litrmas11  litrmas12  gleafmas1  gleafmas2   gleafmas3   ',&
!          'gleafmas4  gleafmas5  gleafmas6  gleafmas7   gleafmas8   ',&
!          'gleafmas9  gleafmas10 gleafmas11 gleafmas12  bleafmas1   ',&
!          'bleafmas2  bleafmas3   bleafmas4  bleafmas5  bleafmas6   ',&
!          'bleafmas7  bleafmas8   bleafmas9  bleafmas10  bleafmas11 ',&
!          'bleafmas12')
! 6905  format (2X, 'litrmass6  tlreleaf6  tltrstem6  tltrroot6  ',&
!          'ltresveg6  humtrsvg6  litrmass7  tltrleaf7  tltrstem7  ',&
!          'tltrroot7  ltresveg7  humtrsvg7  plitrmassms  litrmassms  ',&
!          'litrfallms  ltrestepms  humicmstep  nppmosstep  nppmoss  ',&
!          'anmoss  rgmoss  rmlmoss  gppmoss  Cmossmas  pCmossmas ')
! 6906  format (2X, 'hpd  gavgscms  hutrstep_g  socrestep  resoxic  ',&
!         'resanoxic  socresp(umol/m2/s)  resoxic(umol/m2/s)  ',&
!          'resanoxic(umol/m2/s)')
! 6907  format (2X, 'litresms  litpsims  psisat1  ltrmosclms   ',&
!          'litrmassms  tbar1  q10funcms litrtempms  ratescpo  ',&
!          'ratescpa  Cso  Csa  fto  fta  resoxic  resanoxic   ',&
!          'frac  tsoila  toilo  ewtable  lewtable  tbar1  tbar2  tbar3')
! 6908  format(2X, 'iday  tmoss  cevapms  fwmoss  thliq1  dsmoss  ',&
!          'g_moss  wmoss  rmlmoss  mwce  q10rmlmos  wmosmax  wmosmin')
! 6909  format(2X,'WTBLACC ZSN PREACC EVAPACC ROFACC g12acc g23acc')
! !    --------------YW March 30, 2015 ---------------------------------/
!
!         ENDIF !IF NOT PARALLELRUN
!-------------------------------------------------------------------------------------------------------------

!subroutine class_hh_w(nltest,nmtest,iday,FAREROT,iyear,jdstd,jdsty,jdendd,jdendy,grclarea,onetile_perPFT)

!         real, pointer, dimension(:) ::  fsstar_g
!         real, pointer, dimension(:) ::  flstar_g
!         real, pointer, dimension(:) ::  qh_g
!         real, pointer, dimension(:) ::  qe_g
!         real, pointer, dimension(:) ::  snomlt_g
!         real, pointer, dimension(:) ::  beg_g
!         real, pointer, dimension(:) ::  gtout_g
!         real, pointer, dimension(:) ::  tpn_g
!         real, pointer, dimension(:) ::  altot_g
!         real, pointer, dimension(:) ::  tcn_g
!         real, pointer, dimension(:) ::  tsn_g
!         real, pointer, dimension(:) ::  zsn_g
!
!         real, pointer, dimension(:) :: WSNOROT_g
!         real, pointer, dimension(:) :: ROFSROT_g
!         real, pointer, dimension(:) :: SNOROT_g
!         real, pointer, dimension(:) :: RHOSROT_g
!         real, pointer, dimension(:) :: ROFROT_g
!         real, pointer, dimension(:) :: ZPNDROT_g
!         real, pointer, dimension(:) :: RCANROT_g
!         real, pointer, dimension(:) :: SCANROT_g
!         real, pointer, dimension(:) :: TROFROT_g
!         real, pointer, dimension(:) :: TROOROT_g
!         real, pointer, dimension(:) :: TROBROT_g
!         real, pointer, dimension(:) :: ROFOROT_g
!         real, pointer, dimension(:) :: ROFBROT_g
!         real, pointer, dimension(:) :: TROSROT_g
!         real, pointer, dimension(:) :: FSGVROT_g
!         real, pointer, dimension(:) :: FSGSROT_g
!         real, pointer, dimension(:) :: FLGVROT_g
!         real, pointer, dimension(:) :: FLGSROT_g
!         real, pointer, dimension(:) :: HFSCROT_g
!         real, pointer, dimension(:) :: HFSSROT_g
!         real, pointer, dimension(:) :: HEVCROT_g
!         real, pointer, dimension(:) :: HEVSROT_g
!         real, pointer, dimension(:) :: HMFCROT_g
!         real, pointer, dimension(:) :: HMFNROT_g
!         real, pointer, dimension(:) :: HTCSROT_g
!         real, pointer, dimension(:) :: HTCCROT_g
!         real, pointer, dimension(:) :: FSGGROT_g
!         real, pointer, dimension(:) :: FLGGROT_g
!         real, pointer, dimension(:) :: HFSGROT_g
!         real, pointer, dimension(:) :: HEVGROT_g
!         real, pointer, dimension(:) :: CDHROT_g
!         real, pointer, dimension(:) :: CDMROT_g
!         real, pointer, dimension(:) :: SFCUROT_g
!         real, pointer, dimension(:) :: SFCVROT_g
!         real, pointer, dimension(:) :: ACTLYR_g
!         real, pointer, dimension(:) :: FTABLE_g
!         real, pointer, dimension(:) :: fc_g
!         real, pointer, dimension(:) :: fg_g
!         real, pointer, dimension(:) :: fcs_g
!         real, pointer, dimension(:) :: fgs_g
!         real, pointer, dimension(:) :: PCFCROT_g
!         real, pointer, dimension(:) :: PCLCROT_g
!         real, pointer, dimension(:) :: PCPGROT_g
!         real, pointer, dimension(:) :: QFCFROT_g
!         real, pointer, dimension(:) :: QFGROT_g
!         real, pointer, dimension(:,:) :: QFCROT_g
!         real, pointer, dimension(:) :: ROFCROT_g
!         real, pointer, dimension(:) :: ROFNROT_g
!         real, pointer, dimension(:) :: WTRSROT_g
!         real, pointer, dimension(:) :: WTRGROT_g
!         real, pointer, dimension(:) :: PCPNROT_g
!         real, pointer, dimension(:) :: QFCLROT_g
!         real, pointer, dimension(:) :: QFNROT_g
!         real, pointer, dimension(:) :: WTRCROT_g
!         real, pointer, dimension(:,:) :: rmlvegrow_g
!         real, pointer, dimension(:,:) :: anvegrow_g
!         real, pointer, dimension(:,:) :: HMFGROT_g
!         real, pointer, dimension(:,:) :: HTCROT_g
!         real, pointer, dimension(:,:) :: TBARROT_g
!         real, pointer, dimension(:,:) :: THLQROT_g
!         real, pointer, dimension(:,:) :: THICROT_g
!         real, pointer, dimension(:,:) :: GFLXROT_g
!
!         WSNOROT_g         => ctem_grd%WSNOROT_g
!         ROFSROT_g         => ctem_grd%ROFSROT_g
!         SNOROT_g          => ctem_grd%SNOROT_g
!         RHOSROT_g         => ctem_grd%RHOSROT_g
!         ROFROT_g          => ctem_grd%ROFROT_g
!         ZPNDROT_g         => ctem_grd%ZPNDROT_g
!         RCANROT_g         => ctem_grd%RCANROT_g
!         SCANROT_g         => ctem_grd%SCANROT_g
!         TROFROT_g         => ctem_grd%TROFROT_g
!         TROOROT_g         => ctem_grd%TROOROT_g
!         TROBROT_g         => ctem_grd%TROBROT_g
!         ROFOROT_g         => ctem_grd%ROFOROT_g
!         ROFBROT_g         => ctem_grd%ROFBROT_g
!         TROSROT_g         => ctem_grd%TROSROT_g
!         FSGVROT_g         => ctem_grd%FSGVROT_g
!         FSGSROT_g         => ctem_grd%FSGSROT_g
!         FLGVROT_g         => ctem_grd%FLGVROT_g
!         FLGSROT_g         => ctem_grd%FLGSROT_g
!         HFSCROT_g         => ctem_grd%HFSCROT_g
!         HFSSROT_g         => ctem_grd%HFSSROT_g
!         HEVCROT_g         => ctem_grd%HEVCROT_g
!         HEVSROT_g         => ctem_grd%HEVSROT_g
!         HMFCROT_g         => ctem_grd%HMFCROT_g
!         HMFNROT_g         => ctem_grd%HMFNROT_g
!         HTCSROT_g         => ctem_grd%HTCSROT_g
!         HTCCROT_g         => ctem_grd%HTCCROT_g
!         FSGGROT_g         => ctem_grd%FSGGROT_g
!         FLGGROT_g         => ctem_grd%FLGGROT_g
!         HFSGROT_g         => ctem_grd%HFSGROT_g
!         HEVGROT_g         => ctem_grd%HEVGROT_g
!         CDHROT_g          => ctem_grd%CDHROT_g
!         CDMROT_g          => ctem_grd%CDMROT_g
!         SFCUROT_g         => ctem_grd%SFCUROT_g
!         SFCVROT_g         => ctem_grd%SFCVROT_g
!         ACTLYR_g          => ctem_grd%ACTLYR_g
!         FTABLE_g          => ctem_grd%FTABLE_g
!         fc_g              => ctem_grd%fc_g
!         fg_g              => ctem_grd%fg_g
!         fcs_g             => ctem_grd%fcs_g
!         fgs_g             => ctem_grd%fgs_g
!         PCFCROT_g         => ctem_grd%PCFCROT_g
!         PCLCROT_g         => ctem_grd%PCLCROT_g
!         PCPGROT_g         => ctem_grd%PCPGROT_g
!         QFCFROT_g         => ctem_grd%QFCFROT_g
!         QFGROT_g          => ctem_grd%QFGROT_g
!         QFCROT_g          => ctem_grd%QFCROT_g
!         ROFCROT_g         => ctem_grd%ROFCROT_g
!         ROFNROT_g         => ctem_grd%ROFNROT_g
!         WTRSROT_g         => ctem_grd%WTRSROT_g
!         WTRGROT_g         => ctem_grd%WTRGROT_g
!         PCPNROT_g         => ctem_grd%PCPNROT_g
!         QFCLROT_g         => ctem_grd%QFCLROT_g
!         QFNROT_g          => ctem_grd%QFNROT_g
!         WTRCROT_g         => ctem_grd%WTRCROT_g
!         rmlvegrow_g       => ctem_grd%rmlvegrow_g
!         anvegrow_g        => ctem_grd%anvegrow_g
!         HMFGROT_g         => ctem_grd%HMFGROT_g
!         HTCROT_g          => ctem_grd%HTCROT_g
!         TBARROT_g         => ctem_grd%TBARROT_g
!         THLQROT_g         => ctem_grd%THLQROT_g
!         THICROT_g         => ctem_grd%THICROT_g
!         GFLXROT_g         => ctem_grd%GFLXROT_g
!
!         fsstar_g         => ctem_grd%fsstar_g
!         flstar_g         => ctem_grd%flstar_g
!         qh_g             => ctem_grd%qh_g
!         qe_g             => ctem_grd%qe_g
!         snomlt_g         => ctem_grd%snomlt_g
!         beg_g            => ctem_grd%beg_g
!         gtout_g          => ctem_grd%gtout_g
!         tpn_g            => ctem_grd%tpn_g
!         altot_g          => ctem_grd%altot_g
!         tcn_g            => ctem_grd%tcn_g
!         tsn_g            => ctem_grd%tsn_g
!         zsn_g            => ctem_grd%zsn_g
!
!         ! mosaic level variables (CLASS):
!
!         fsnowacc_t        => ctem_tile%fsnowacc_t
!         tcansacc_t        => ctem_tile%tcansacc_t
!         tcanoaccgat_t     => ctem_tile%tcanoaccgat_t
!         taaccgat_t        => ctem_tile%taaccgat_t
!         uvaccgat_t        => ctem_tile%uvaccgat_t
!         vvaccgat_t        => ctem_tile%vvaccgat_t
!         tbaraccgat_t      => ctem_tile%tbaraccgat_t
!         tbarcacc_t        => ctem_tile%tbarcacc_t
!         tbarcsacc_t       => ctem_tile%tbarcsacc_t
!         tbargacc_t        => ctem_tile%tbargacc_t
!         tbargsacc_t       => ctem_tile%tbargsacc_t
!         thliqcacc_t       => ctem_tile%thliqcacc_t
!         thliqgacc_t       => ctem_tile%thliqgacc_t
!         thliqacc_t        => ctem_tile%thliqacc_t
!         thiceacc_t        => ctem_tile%thiceacc_t  ! Added in place of YW's thicaccgat_m. EC Dec 23 2016.
!         thicecacc_t       => ctem_tile%thicecacc_t
!         thicegacc_t       => ctem_tile%thicegacc_t
!         ancsvgac_t        => ctem_tile%ancsvgac_t
!         ancgvgac_t        => ctem_tile%ancgvgac_t
!         rmlcsvga_t        => ctem_tile%rmlcsvga_t
!         rmlcgvga_t        => ctem_tile%rmlcgvga_t
!         anmossac_t        => ctem_tile%anmossac_t
!         rmlmossac_t       => ctem_tile%rmlmossac_t
!         gppmossac_t       => ctem_tile%gppmossac_t

!             if (dohhoutput) then ! stand alone mode, include half-hourly output for CLASS & CTEM
!
!                 DO 450 I=1,NLTEST
!
!                     !       initialization of various grid-averaged variables
!                     call resetgridavg(nltest)
!
!                     DO 425 M=1,NMTEST
!                         IF(FSSROW(I).GT.0.0) THEN
!                             ALTOT=(FSSROW(I)-(FSGVROT(I,M)+FSGSROT(I,M)&
!                                 &              +FSGGROT(I,M)))/FSSROW(I)
!                         ELSE
!                             ALTOT=0.0
!                         ENDIF
!                         FSSTAR=FSSROW(I)*(1.0-ALTOT)
!                         FLSTAR=FDLROW(I)-SBC*GTROT(I,M)**4
!                         QH=HFSROT(I,M)
!                         QE=QEVPROT(I,M)
!                         !          BEG=FSSTAR+FLSTAR-QH-QE !(commented out in runclass.fieldsite)
!                         BEG=GFLXGAT(1,1)  !FLAG!
!                         !          USTARBS=UVROW(1)*SQRT(CDMROT(I,M)) !FLAG (commented out in runclass.fieldsite)
!                         SNOMLT=HMFNROT(I,M)
!                         IF(RHOSROT(I,M).GT.0.0) THEN
!                             ZSN=SNOROT(I,M)/RHOSROT(I,M)
!                         ELSE
!                             ZSN=0.0
!                         ENDIF
!                         IF(TCANROT(I,M).GT.0.01) THEN
!                             TCN=TCANROT(I,M)-TFREZ
!                         ELSE
!                             TCN=0.0
!                         ENDIF
!                         TSURF=FCS(I)*TSFSGAT(I,1)+FGS(I)*TSFSGAT(I,2)+&
!                             &           FC(I)*TSFSGAT(I,3)+FG(I)*TSFSGAT(I,4)
!                         !          IF(FSSROW(I).GT.0.0 .AND. (FCS(I)+FC(I)).GT.0.0) THEN
!                         !          IF(FSSROW(I).GT.0.0) THEN
!                         NFS=NFS+1
!                         ITA=NINT(TAROW(I)-TFREZ)
!                         ITCAN=NINT(TCN)
!                         ITAC=NINT(TACGAT(I)-TFREZ)
!                         ITSCR=NINT(SFCTGAT(I)-TFREZ)
!                         ITS=NINT(TSURF-TFREZ)
!                         !              ITD=ITS-ITA
!                         ITD=ITCAN-ITA
!                         ITD2=ITCAN-ITSCR
!                         ITD3=ITCAN-ITAC
!                         ITD4=ITAC-ITA
!                         !          ENDIF
!                         IF(FC(I).GT.0.1 .AND. RC(I).GT.1.0E5) NDRY=NDRY+1
!                         !           IF((ITCAN-ITA).GE.10) THEN
!                         !               WRITE(6,6070) IHOUR,IMIN,IDAY,IYEAR,FSSTAR,FLSTAR,QH,QE,
!                         !      1                      BEG,TAROW(I)-TFREZ,TCN,TCN-(TAROW(I)-TFREZ),
!                         !      2                      PAICAN(I),FSVF(I),UVROW(I),RC(I)
!                         ! 6070          FORMAT(2X,2I2,I4,I5,9F6.1,F6.3,F6.1,F8.1)
!                         !           ENDIF
!                         !
!                         IF(TSNOROT(I,M).GT.0.01) THEN
!                             TSN=TSNOROT(I,M)-TFREZ
!                         ELSE
!                             TSN=0.0
!                         ENDIF
!                         IF(TPNDROT(I,M).GT.0.01) THEN
!                             TPN=TPNDROT(I,M)-TFREZ
!                         ELSE
!                             TPN=0.0
!                         ENDIF
!                         GTOUT=GTROT(I,M)-TFREZ
!                         EVAPSUM=QFCFROT(I,M)+QFCLROT(I,M)+QFNROT(I,M)+QFGROT(I,M)+&
!                             &                   QFCROT(I,M,1)+QFCROT(I,M,2)+QFCROT(I,M,3)
!
!                         ! start writing output
!
!                         if ((iyear .ge. jhhsty) .and. (iyear .le. jhhendy)) then
!                             if ((iday .ge. jhhstd) .and. (iday .le. jhhendd)) then
!
!                                 WRITE(64,6400) IHOUR,IMIN,IDAY,IYEAR,FSSTAR,FLSTAR,QH,QE,&
!                                     &                   SNOMLT,BEG,GTOUT,SNOROT(I,M),RHOSROT(I,M),&
!                                     &                   WSNOROT(I,M),ALTOT,ROFROT(I,M),&
!                                     &                   TPN,ZPNDROT(I,M),CDHROT(I,M),CDMROT(I,M),&
!                                     &                   SFCUROT(I,M),SFCVROT(I,M),UVROW(I),' TILE ',m
!                                 IF(IGND.GT.3) THEN
!                                     write(66,6601) ihour,imin,iday,iyear,(TBARROT(i,m,j)-&
!                                         &                 tfrez,THLQROT(i,m,j),THICROT(i,m,j),j=1,IGND),&
!                                         &                 (GFLXROT(i,m,j),j=1,IGND),' TILE ',m
!                                 end if
!
!                                 write(65,6500) ihour,imin,iday,iyear,(TBARROT(i,m,j)-&
!                                     &                   tfrez,THLQROT(i,m,j),THICROT(i,m,j),j=1,3),&
!                                     &                  tcn,RCANROT(i,m),SCANROT(i,m),tsn,zsn,&
!                                     &                   TCN-(TAROW(I)-TFREZ),TCANO(I)-TFREZ,&
!                                     &                   TACGAT(I)-TFREZ,' TILE ',m
!                                 !
!                                 WRITE(67,6700) IHOUR,IMIN,IDAY,IYEAR,&
!                                     &                   TROFROT(I,M),TROOROT(I,M),TROSROT(I,M),&
!                                     &                   TROBROT(I,M),ROFROT(I,M),ROFOROT(I,M),&
!                                     &                   ROFSROT(I,M),ROFBROT(I,M),&
!                                     &                   FCS(M),FGS(M),FC(M),FG(M),' TILE ',M
!                                 WRITE(68,6800) IHOUR,IMIN,IDAY,IYEAR,&
!                                     &                   FSGVROT(I,M),FSGSROT(I,M),FSGGROT(I,M),&
!                                     &                   FLGVROT(I,M),FLGSROT(I,M),FLGGROT(I,M),&
!                                     &                   HFSCROT(I,M),HFSSROT(I,M),HFSGROT(I,M),&
!                                     &                   HEVCROT(I,M),HEVSROT(I,M),HEVGROT(I,M),&
!                                     &                   HMFCROT(I,M),HMFNROT(I,M),&
!                                     &                   (HMFGROT(I,M,J),J=1,3),&
!                                     &                   HTCCROT(I,M),HTCSROT(I,M),&
!                                     &                   (HTCROT(I,M,J),J=1,3),' TILE ',M
!                                 WRITE(69,6900) IHOUR,IMIN,IDAY,IYEAR,&
!                                     &                   PCFCROT(I,M),PCLCROT(I,M),PCPNROT(I,M),&
!                                     &                   PCPGROT(I,M),QFCFROT(I,M),QFCLROT(I,M),&
!                                     &                   QFNROT(I,M),QFGROT(I,M),(QFCROT(I,M,J),J=1,3),&
!                                     &                   ROFCROT(I,M),ROFNROT(I,M),ROFOROT(I,M),&
!                                     &                   ROFROT(I,M),WTRCROT(I,M),WTRSROT(I,M),&
!                                     &                   WTRGROT(I,M),' TILE ',M
!                             endif
!                         endif ! half hourly output loop.
!                         !
!                         ! Write half-hourly CTEM results to file *.CT01H
!                         !
!                         ! Net photosynthetic rates and leaf maintenance respiration for
!                         ! each pft. however, if ctem_on then physyn subroutine
!                         ! is using storage lai while actual lai is zero. if actual lai is
!                         ! zero then we make anveg and rmlveg zero as well because these
!                         ! are imaginary just like storage lai. note that anveg and rmlveg
!                         ! are not passed to ctem. rather ancsveg, ancgveg, rmlcsveg, and
!                         ! rmlcgveg are passed.
!                         !
!                         if (ctem_on) then
!
!                             do 760 j = 1,icc
!                                 if(ailcgrow(i,m,j).le.0.0) then
!                                     anvegrow(i,m,j)=0.0
!                                     rmlvegrow(i,m,j)=0.0
!                                 else
!                                     anvegrow(i,m,j)=ancsvegrow(i,m,j)*FSNOROT(i,m) +&
!                                         &                          ancgvegrow(i,m,j)*(1. - FSNOROT(i,m))
!                                     rmlvegrow(i,m,j)=rmlcsvegrow(i,m,j)*FSNOROT(i,m) +&
!                                         &                         rmlcgvegrow(i,m,j)*(1. - FSNOROT(i,m))
!                                 endif
! 760                         continue
!
!                             if ((iyear .ge. jhhsty) .and. (iyear .le. jhhendy)) then
!                                 if ((iday .ge. jhhstd) .and. (iday .le. jhhendd)) then
!
!                                     write(71,7200)ihour,imin,iday,iyear,(anvegrow(i,m,j),&
!                                         &                    j=1,icc),(rmlvegrow(i,m,j),j=1,icc),&
!                                         &                    ' TILE ',m
!                                 endif
!                             end if
!
!                             do j = 1,icc
!                                 anvegrow_g(i,j)=anvegrow_g(i,j)+anvegrow(i,m,j)&
!                                     &                                        *FAREROT(i,m)
!                                 rmlvegrow_g(i,j)=rmlvegrow_g(i,j)+rmlvegrow(i,m,j)&
!                                     &                                         *FAREROT(i,m)
!                             enddo
!                         endif   ! ctem_on
!
! 7200      format(1x,i2,1x,i2,i5,i5,12f11.3,12f11.3,2(a6,i2))
!
!                         fsstar_g(i)    =fsstar_g(i) + fsstar*FAREROT(i,m)
!                         flstar_g(i)    =flstar_g(i) + flstar*FAREROT(i,m)
!                         qh_g(i)        =qh_g(i)     + qh*FAREROT(i,m)
!                         qe_g(i)        =qe_g(i)     + qe*FAREROT(i,m)
!                         snomlt_g(i)    =snomlt_g(i) + snomlt*FAREROT(i,m)
!                         beg_g(i)       =beg_g(i)    + beg*FAREROT(i,m)
!                         gtout_g(i)     =gtout_g(i)  + gtout*FAREROT(i,m)
!                         tcn_g(i)       =tcn_g(i)    + tcn*FAREROT(i,m)
!                         tsn_g(i)       =tsn_g(i)    + tsn*FAREROT(i,m)
!                         zsn_g(i)       =zsn_g(i)    + zsn*FAREROT(i,m)
!                         altot_g(i)     =altot_g(i)  + altot*FAREROT(i,m)
!                         tpn_g(i)       =tpn_g(i)    + tpn*FAREROT(i,m)
!
!                         do j=1,ignd
!                             TBARROT_g(i,j)=TBARROT_g(i,j) + TBARROT(i,m,j)*FAREROT(i,m)
!                             THLQROT_g(i,j)=THLQROT_g(i,j) + THLQROT(i,m,j)*FAREROT(i,m)
!                             THICROT_g(i,j)=THICROT_g(i,j) + THICROT(i,m,j)*FAREROT(i,m)
!                             GFLXROT_g(i,j)=GFLXROT_g(i,j) + GFLXROT(i,m,j)*FAREROT(i,m)
!                             HMFGROT_g(i,j)=HMFGROT_g(i,j) + HMFGROT(i,m,j)*FAREROT(i,m)
!                             HTCROT_g(i,j)=HTCROT_g(i,j) + HTCROT(i,m,j)*FAREROT(i,m)
!                             QFCROT_g(i,j)=QFCROT_g(i,j) + QFCROT(i,m,j)*FAREROT(i,m)
!                         enddo
!
!                         ZPNDROT_g(i)=ZPNDROT_g(i) + ZPNDROT(i,m)*FAREROT(i,m)
!                         RHOSROT_g(i)=RHOSROT_g(i) + RHOSROT(i,m)*FAREROT(i,m)
!                         WSNOROT_g(i)=WSNOROT_g(i) + WSNOROT(i,m)*FAREROT(i,m)
!                         RCANROT_g(i)=RCANROT_g(i) + RCANROT(i,m)*FAREROT(i,m)
!                         SCANROT_g(i)=SCANROT_g(i) + SCANROT(i,m)*FAREROT(i,m)
!                         TROFROT_g(i)=TROFROT_g(i) + TROFROT(i,m)*FAREROT(i,m)
!                         TROOROT_g(i)=TROOROT_g(i) + TROOROT(i,m)*FAREROT(i,m)
!                         TROSROT_g(i)=TROSROT_g(i) + TROSROT(i,m)*FAREROT(i,m)
!                         TROBROT_g(i)=TROBROT_g(i) + TROBROT(i,m)*FAREROT(i,m)
!                         ROFOROT_g(i)=ROFOROT_g(i) + ROFOROT(i,m)*FAREROT(i,m)
!                         ROFSROT_g(i)=ROFSROT_g(i) + ROFSROT(i,m)*FAREROT(i,m)
!                         ROFBROT_g(i)=ROFBROT_g(i) + ROFBROT(i,m)*FAREROT(i,m)
!                         FSGVROT_g(i)=FSGVROT_g(i) + FSGVROT(i,m)*FAREROT(i,m)
!                         FSGSROT_g(i)=FSGSROT_g(i) + FSGSROT(i,m)*FAREROT(i,m)
!                         FSGGROT_g(i)=FSGGROT_g(i) + FSGGROT(i,m)*FAREROT(i,m)
!                         FLGVROT_g(i)=FLGVROT_g(i) + FLGVROT(i,m)*FAREROT(i,m)
!                         FLGSROT_g(i)=FLGSROT_g(i) + FLGSROT(i,m)*FAREROT(i,m)
!                         FLGGROT_g(i)=FLGGROT_g(i) + FLGGROT(i,m)*FAREROT(i,m)
!                         HFSCROT_g(i)=HFSCROT_g(i) + HFSCROT(i,m)*FAREROT(i,m)
!                         HFSSROT_g(i)=HFSSROT_g(i) + HFSSROT(i,m)*FAREROT(i,m)
!                         HFSGROT_g(i)=HFSGROT_g(i) + HFSGROT(i,m)*FAREROT(i,m)
!                         HEVCROT_g(i)=HEVCROT_g(i) + HEVCROT(i,m)*FAREROT(i,m)
!                         HEVSROT_g(i)=HEVSROT_g(i) + HEVSROT(i,m)*FAREROT(i,m)
!                         HEVGROT_g(i)=HEVGROT_g(i) + HEVGROT(i,m)*FAREROT(i,m)
!                         HMFCROT_g(i)=HMFCROT_g(i) + HMFCROT(i,m)*FAREROT(i,m)
!                         HMFNROT_g(i)=HMFNROT_g(i) + HMFNROT(i,m)*FAREROT(i,m)
!                         HTCCROT_g(i)=HTCCROT_g(i) + HTCCROT(i,m)*FAREROT(i,m)
!                         HTCSROT_g(i)=HTCSROT_g(i) + HTCSROT(i,m)*FAREROT(i,m)
!                         PCFCROT_g(i)=PCFCROT_g(i) + PCFCROT(i,m)*FAREROT(i,m)
!                         PCLCROT_g(i)=PCLCROT_g(i) + PCLCROT(i,m)*FAREROT(i,m)
!                         PCPNROT_g(i)=PCPNROT_g(i) + PCPNROT(i,m)*FAREROT(i,m)
!                         PCPGROT_g(i)=PCPGROT_g(i) + PCPGROT(i,m)*FAREROT(i,m)
!                         QFCFROT_g(i)=QFCFROT_g(i) + QFCFROT(i,m)*FAREROT(i,m)
!                         QFCLROT_g(i)=QFCLROT_g(i) + QFCLROT(i,m)*FAREROT(i,m)
!                         ROFCROT_g(i)=ROFCROT_g(i) + ROFCROT(i,m)*FAREROT(i,m)
!                         ROFNROT_g(i)=ROFNROT_g(i) + ROFNROT(i,m)*FAREROT(i,m)
!                         WTRCROT_g(i)=WTRCROT_g(i) + WTRCROT(i,m)*FAREROT(i,m)
!                         WTRSROT_g(i)=WTRSROT_g(i) + WTRSROT(i,m)*FAREROT(i,m)
!                         WTRGROT_g(i)=WTRGROT_g(i) + WTRGROT(i,m)*FAREROT(i,m)
!                         QFNROT_g(i) =QFNROT_g(i) + QFNROT(i,m)*FAREROT(i,m)
!                         QFGROT_g(i) =QFGROT_g(i) + QFGROT(i,m)*FAREROT(i,m)
!                         ROFROT_g(i) =ROFROT_g(i) + ROFROT(i,m)*FAREROT(i,m)
!                         SNOROT_g(i) =SNOROT_g(i) + SNOROT(i,m)*FAREROT(i,m)
!                         CDHROT_g(i) =CDHROT_g(i) + CDHROT(i,m)*FAREROT(i,m)
!                         CDMROT_g(i) =CDMROT_g(i) + CDMROT(i,m)*FAREROT(i,m)
!                         SFCUROT_g(i) =SFCUROT_g(i) + SFCUROT(i,m)*FAREROT(i,m)
!                         SFCVROT_g(i) =SFCVROT_g(i) + SFCVROT(i,m)*FAREROT(i,m)
!                         ACTLYR_g(i) = ACTLYR_g(i) + ACTLYR(i,m) * FAREROT(i,m)
!                         FTABLE_g(i) = FTABLE_g(i) + FTABLE(i,m) * FAREROT(i,m)
!
! 425                 CONTINUE
!
!                     if ((iyear .ge. jhhsty) .and. (iyear .le. jhhendy)) then
!                         if ((iday .ge. jhhstd) .and. (iday .le. jhhendd)) then
!
!                             IF (CTEM_ON) THEN
!                                 WRITE(711,7200)IHOUR,IMIN,IDAY,IYEAR,(ANVEGROW_G(I,J),&
!                                     &                 J=1,ICC),(RMLVEGROW_G(I,J),J=1,ICC)
!                             ENDIF !CTEM_ON
!
!                             WRITE(641,6400) IHOUR,IMIN,IDAY,IYEAR,FSSTAR_G(i),FLSTAR_G(i),&
!                                 &                  QH_G(i),QE_G(i),SNOMLT_G(i),BEG_G(i),GTOUT_G(i),&
!                                 &                  SNOROT_G(I),RHOSROT_G(I),WSNOROT_G(I),&
!                                 &                  ALTOT_G(i),ROFROT_G(I),TPN_G(i),ZPNDROT_G(I),&
!                                 &                  CDHROT_G(I),CDMROT_G(I),SFCUROT_G(I),&
!                                 &                  SFCVROT_G(I),UVROW(I)
!                             WRITE(651,6500) IHOUR,IMIN,IDAY,IYEAR,(TBARROT_G(I,J)-&
!                                 &                   TFREZ,THLQROT_G(I,J),THICROT_G(I,J),J=1,3),&
!                                 &                   TCN_G(i),RCANROT_G(I),SCANROT_G(I),TSN_G(i),&
!                                 &                   ZSN_G(i),TCN_G(i)-(TAROW(I)-TFREZ),&
!                                 &                   TCANO(I)-TFREZ,TACGAT(I)-TFREZ
!
!                             IF(IGND.GT.3) THEN
!                                 WRITE(661,6601) IHOUR,IMIN,IDAY,IYEAR,(TBARROT_G(I,J)-&
!                                     &                   TFREZ,THLQROT_G(I,J),THICROT_G(I,J),J=1,IGND),&
!                                     &                   (GFLXROT_G(I,J),J=1,IGND)
!                             ELSE
!                                 WRITE(661,6600) IHOUR,IMIN,IDAY,FSSROW(I),FDLROW(I),PREROW(I),&
!                                     &                   TAROW(I)-TFREZ,UVROW(I),PRESROW(I),QAROW(I)
!                             ENDIF
!
!                             WRITE(671,6700) IHOUR,IMIN,IDAY,IYEAR,&
!                                 &                   TROFROT_G(I),TROOROT_G(I),TROSROT_G(I),&
!                                 &                   TROBROT_G(I),ROFROT_G(I),ROFOROT_G(I),&
!                                 &                   ROFSROT_G(I),ROFBROT_G(I),&
!                                 &                   FCS_G(I),FGS_G(I),FC_G(I),FG_G(I)
!                             WRITE(681,6800) IHOUR,IMIN,IDAY,IYEAR,&
!                                 &                   FSGVROT_G(I),FSGSROT_G(I),FSGGROT_G(I),&
!                                 &                   FLGVROT_G(I),FLGSROT_G(I),FLGGROT_G(I),&
!                                 &                   HFSCROT_G(I),HFSSROT_G(I),HFSGROT_G(I),&
!                                 &                   HEVCROT_G(I),HEVSROT_G(I),HEVGROT_G(I),&
!                                 &                   HMFCROT_G(I),HMFNROT_G(I),&
!                                 &                   (HMFGROT_G(I,J),J=1,3),&
!                                 &                   HTCCROT_G(I),HTCSROT_G(I),&
!                                 &                   (HTCROT_G(I,J),J=1,3)
!                             WRITE(691,6900) IHOUR,IMIN,IDAY,IYEAR,&
!                                 &                   PCFCROT_G(I),PCLCROT_G(I),PCPNROT_G(I),&
!                                 &                   PCPGROT_G(I),QFCFROT_G(I),QFCLROT_G(I),&
!                                 &                   QFNROT_G(I),QFGROT_G(I),(QFCROT_G(I,J),J=1,3),&
!                                 &                   ROFCROT_G(I),ROFNROT_G(I),ROFOROT_G(I),&
!                                 &                   ROFROT_G(I),WTRCROT_G(I),WTRSROT_G(I),&
!                                 &                   WTRGROT_G(I)
!                         endif
!                     ENDIF ! if write half-hourly
!
! 450             CONTINUE
!
!             endif ! dohhoutput
!subroutine class_hh_w

!==============================================================================================================

    !==============================================================================================================
    !>\ingroup io_driver_class_daily_aw
    !>@{

    !class_daily_aw
!
!         real, pointer, dimension(:,:) :: PREACC_M
!         real, pointer, dimension(:,:) :: GTACC_M
!         real, pointer, dimension(:,:) :: QEVPACC_M
!         real, pointer, dimension(:,:) :: HFSACC_M
!         real, pointer, dimension(:,:) :: HMFNACC_M
!         real, pointer, dimension(:,:) :: ROFACC_M
!         real, pointer, dimension(:,:) :: SNOACC_M
!         real, pointer, dimension(:,:) :: OVRACC_M
!         real, pointer, dimension(:,:) :: WTBLACC_M
!         real, pointer, dimension(:,:,:) :: TBARACC_M
!         real, pointer, dimension(:,:,:) :: THLQACC_M
!         real, pointer, dimension(:,:,:) :: THICACC_M
!         real, pointer, dimension(:,:,:) :: THALACC_M
!         real, pointer, dimension(:,:) :: ALVSACC_M
!         real, pointer, dimension(:,:) :: ALIRACC_M
!         real, pointer, dimension(:,:) :: RHOSACC_M
!         real, pointer, dimension(:,:) :: TSNOACC_M
!         real, pointer, dimension(:,:) :: WSNOACC_M
!         real, pointer, dimension(:,:) :: SNOARE_M
!         real, pointer, dimension(:,:) :: TCANACC_M
!         real, pointer, dimension(:,:) :: RCANACC_M
!         real, pointer, dimension(:,:) :: SCANACC_M
!         real, pointer, dimension(:,:) :: GROACC_M
!         real, pointer, dimension(:,:) :: FSINACC_M
!         real, pointer, dimension(:,:) :: FLINACC_M
!         real, pointer, dimension(:,:) :: TAACC_M
!         real, pointer, dimension(:,:) :: UVACC_M
!         real, pointer, dimension(:,:) :: PRESACC_M
!         real, pointer, dimension(:,:) :: QAACC_M
!         real, pointer, dimension(:,:) :: ALTOTACC_M
!         real, pointer, dimension(:,:) :: EVAPACC_M
!         real, pointer, dimension(:,:) :: FLUTACC_M

!         PREACC_M          => class_rot%PREACC_M
!         GTACC_M           => class_rot%GTACC_M
!         QEVPACC_M         => class_rot%QEVPACC_M
!         HFSACC_M          => class_rot%HFSACC_M
!         HMFNACC_M         => class_rot%HMFNACC_M
!         ROFACC_M          => class_rot%ROFACC_M
!         SNOACC_M          => class_rot%SNOACC_M
!         OVRACC_M          => class_rot%OVRACC_M
!         WTBLACC_M         => class_rot%WTBLACC_M
!         TBARACC_M         => class_rot%TBARACC_M
!         THLQACC_M         => class_rot%THLQACC_M
!         THICACC_M         => class_rot%THICACC_M
!         THALACC_M         => class_rot%THALACC_M
!         ALVSACC_M         => class_rot%ALVSACC_M
!         ALIRACC_M         => class_rot%ALIRACC_M
!         RHOSACC_M         => class_rot%RHOSACC_M
!         TSNOACC_M         => class_rot%TSNOACC_M
!         WSNOACC_M         => class_rot%WSNOACC_M
!         SNOARE_M          => class_rot%SNOARE_M
!         TCANACC_M         => class_rot%TCANACC_M
!         RCANACC_M         => class_rot%RCANACC_M
!         SCANACC_M         => class_rot%SCANACC_M
!         GROACC_M          => class_rot%GROACC_M
!         FSINACC_M         => class_rot%FSINACC_M
!         FLINACC_M         => class_rot%FLINACC_M
!         TAACC_M           => class_rot%TAACC_M
!         UVACC_M           => class_rot%UVACC_M
!         PRESACC_M         => class_rot%PRESACC_M
!         QAACC_M           => class_rot%QAACC_M
!         ALTOTACC_M        => class_rot%ALTOTACC_M
!         EVAPACC_M         => class_rot%EVAPACC_M
!         FLUTACC_M         => class_rot%FLUTACC_M

            !     * CALCULATE GRID CELL AVERAGE DIAGNOSTIC FIELDS.

!             if(dohhoutput) then !or dodayoutput???          ! stand alone mode, includes diagnostic fields
!
!                 DO 525 I=1,NLTEST
!                     CDHROW(I)=0.
!                     CDMROW(I)=0.
!                     HFSROW(I)=0.
!                     TFXROW(I)=0.
!                     QEVPROW(I)=0.
!                     QFSROW(I)=0.
!                     QFXROW(I)=0.
!                     PETROW(I)=0.
!                     GAROW(I)=0.
!                     EFROW(I)=0.
!                     GTROW(I)=0.
!                     QGROW(I)=0.
!                     ALVSROW(I)=0.
!                     ALIRROW(I)=0.
!                     SFCTROW(I)=0.
!                     SFCUROW(I)=0.
!                     SFCVROW(I)=0.
!                     SFCQROW(I)=0.
!                     SFRHROW(I)=0.
!                     FSNOROW(I)=0.
!                     FSGVROW(I)=0.
!                     FSGSROW(I)=0.
!                     FSGGROW(I)=0.
!                     FLGVROW(I)=0.
!                     FLGSROW(I)=0.
!                     FLGGROW(I)=0.
!                     HFSCROW(I)=0.
!                     HFSSROW(I)=0.
!                     HFSGROW(I)=0.
!                     HEVCROW(I)=0.
!                     HEVSROW(I)=0.
!                     HEVGROW(I)=0.
!                     HMFCROW(I)=0.
!                     HMFNROW(I)=0.
!                     HTCCROW(I)=0.
!                     HTCSROW(I)=0.
!                     PCFCROW(I)=0.
!                     PCLCROW(I)=0.
!                     PCPNROW(I)=0.
!                     PCPGROW(I)=0.
!                     QFGROW(I)=0.
!                     QFNROW(I)=0.
!                     QFCLROW(I)=0.
!                     QFCFROW(I)=0.
!                     ROFROW(I)=0.
!                     ROFOROW(I)=0.
!                     ROFSROW(I)=0.
!                     ROFBROW(I)=0.
!                     ROFCROW(I)=0.
!                     ROFNROW(I)=0.
!                     ROVGROW(I)=0.
!                     WTRCROW(I)=0.
!                     WTRSROW(I)=0.
!                     WTRGROW(I)=0.
!                     DRROW(I)=0.
!                     wtableROW(I)=0.
!                     ILMOROW(I)=0.
!                     UEROW(I)=0.
!                     HBLROW(I)=0.
!                     G12GRD(I)= 0.       !YW March 27, 2015
!                     G23GRD(I)= 0.       !YW March 27, 2015
!                     DO 500 J=1,IGND
!                         HMFGROW(I,J)=0.
!                         HTCROW(I,J)=0.
!                         QFCROW(I,J)=0.
!                         GFLXROW(I,J)=0.
! 500                 CONTINUE
! 525             CONTINUE
!
!                 DO 600 I=1,NLTEST
!                     DO 575 M=1,NMTEST
!                         CDHROW(I)=CDHROW(I)+CDHROT(I,M)*FAREROT(I,M)
!                         CDMROW(I)=CDMROW(I)+CDMROT(I,M)*FAREROT(I,M)
!                         HFSROW(I)=HFSROW(I)+HFSROT(I,M)*FAREROT(I,M)
!                         TFXROW(I)=TFXROW(I)+TFXROT(I,M)*FAREROT(I,M)
!                         QEVPROW(I)=QEVPROW(I)+QEVPROT(I,M)*FAREROT(I,M)
!                         QFSROW(I)=QFSROW(I)+QFSROT(I,M)*FAREROT(I,M)
!                         QFXROW(I)=QFXROW(I)+QFXROT(I,M)*FAREROT(I,M)
!                         PETROW(I)=PETROW(I)+PETROT(I,M)*FAREROT(I,M)
!                         GAROW(I)=GAROW(I)+GAROT(I,M)*FAREROT(I,M)
!                         EFROW(I)=EFROW(I)+EFROT(I,M)*FAREROT(I,M)
!                         GTROW(I)=GTROW(I)+GTROT(I,M)*FAREROT(I,M)
!                         QGROW(I)=QGROW(I)+QGROT(I,M)*FAREROT(I,M)
!                         ALVSROW(I)=ALVSROW(I)+ALVSROT(I,M)*FAREROT(I,M)
!                         ALIRROW(I)=ALIRROW(I)+ALIRROT(I,M)*FAREROT(I,M)
!                         SFCTROW(I)=SFCTROW(I)+SFCTROT(I,M)*FAREROT(I,M)
!                         SFCUROW(I)=SFCUROW(I)+SFCUROT(I,M)*FAREROT(I,M)
!                         SFCVROW(I)=SFCVROW(I)+SFCVROT(I,M)*FAREROT(I,M)
!                         SFCQROW(I)=SFCQROW(I)+SFCQROT(I,M)*FAREROT(I,M)
!                         SFRHROW(I)=SFRHROW(I)+SFRHROT(I,M)*FAREROT(I,M)
!                         FSNOROW(I)=FSNOROW(I)+FSNOROT(I,M)*FAREROT(I,M)
!                         FSGVROW(I)=FSGVROW(I)+FSGVROT(I,M)*FAREROT(I,M)
!                         FSGSROW(I)=FSGSROW(I)+FSGSROT(I,M)*FAREROT(I,M)
!                         FSGGROW(I)=FSGGROW(I)+FSGGROT(I,M)*FAREROT(I,M)
!                         FLGVROW(I)=FLGVROW(I)+FLGVROT(I,M)*FAREROT(I,M)
!                         FLGSROW(I)=FLGSROW(I)+FLGSROT(I,M)*FAREROT(I,M)
!                         FLGGROW(I)=FLGGROW(I)+FLGGROT(I,M)*FAREROT(I,M)
!                         HFSCROW(I)=HFSCROW(I)+HFSCROT(I,M)*FAREROT(I,M)
!                         HFSSROW(I)=HFSSROW(I)+HFSSROT(I,M)*FAREROT(I,M)
!                         HFSGROW(I)=HFSGROW(I)+HFSGROT(I,M)*FAREROT(I,M)
!                         HEVCROW(I)=HEVCROW(I)+HEVCROT(I,M)*FAREROT(I,M)
!                         HEVSROW(I)=HEVSROW(I)+HEVSROT(I,M)*FAREROT(I,M)
!                         HEVGROW(I)=HEVGROW(I)+HEVGROT(I,M)*FAREROT(I,M)
!                         HMFCROW(I)=HMFCROW(I)+HMFCROT(I,M)*FAREROT(I,M)
!                         HMFNROW(I)=HMFNROW(I)+HMFNROT(I,M)*FAREROT(I,M)
!                         HTCCROW(I)=HTCCROW(I)+HTCCROT(I,M)*FAREROT(I,M)
!                         HTCSROW(I)=HTCSROW(I)+HTCSROT(I,M)*FAREROT(I,M)
!                         PCFCROW(I)=PCFCROW(I)+PCFCROT(I,M)*FAREROT(I,M)
!                         PCLCROW(I)=PCLCROW(I)+PCLCROT(I,M)*FAREROT(I,M)
!                         PCPNROW(I)=PCPNROW(I)+PCPNROT(I,M)*FAREROT(I,M)
!                         PCPGROW(I)=PCPGROW(I)+PCPGROT(I,M)*FAREROT(I,M)
!                         QFGROW(I)=QFGROW(I)+QFGROT(I,M)*FAREROT(I,M)
!                         QFNROW(I)=QFNROW(I)+QFNROT(I,M)*FAREROT(I,M)
!                         QFCLROW(I)=QFCLROW(I)+QFCLROT(I,M)*FAREROT(I,M)
!                         QFCFROW(I)=QFCFROW(I)+QFCFROT(I,M)*FAREROT(I,M)
!                         ROFROW(I)=ROFROW(I)+ROFROT(I,M)*FAREROT(I,M)
!                         ROFOROW(I)=ROFOROW(I)+ROFOROT(I,M)*FAREROT(I,M)
!                         ROFSROW(I)=ROFSROW(I)+ROFSROT(I,M)*FAREROT(I,M)
!                         ROFBROW(I)=ROFBROW(I)+ROFBROT(I,M)*FAREROT(I,M)
!                         ROFCROW(I)=ROFCROW(I)+ROFCROT(I,M)*FAREROT(I,M)
!                         ROFNROW(I)=ROFNROW(I)+ROFNROT(I,M)*FAREROT(I,M)
!                         ROVGROW(I)=ROVGROW(I)+ROVGROT(I,M)*FAREROT(I,M)
!                         WTRCROW(I)=WTRCROW(I)+WTRCROT(I,M)*FAREROT(I,M)
!                         WTRSROW(I)=WTRSROW(I)+WTRSROT(I,M)*FAREROT(I,M)
!                         WTRGROW(I)=WTRGROW(I)+WTRGROT(I,M)*FAREROT(I,M)
!                         DRROW(I)=DRROW(I)+DRROT(I,M)*FAREROT(I,M)
!                         wtableROW(I)=wtableROW(I)+wtableROT(I,M)*FAREROT(I,M)
!                         ILMOROW(I)=ILMOROW(I)+ILMOROT(I,M)*FAREROT(I,M)
!                         UEROW(I)=UEROW(I)+UEROT(I,M)*FAREROT(I,M)
!                         HBLROW(I)=HBLROW(I)+HBLROT(I,M)*FAREROT(I,M)
!                         G12GRD(I)=G12C(I)*FC(I)+G12G(I)*FG(I)+G12CS(I)*FCS(I)+&
!                         G12GS(I)*FGS(I)     !YW March 27, 2015
!                         G23GRD(I)=G23C(I)*FC(I)+G23G(I)*FG(I)+G23CS(I)*FCS(I)+&
!                         G23GS(I)*FGS(I)     !YW March 27, 2015
!                         DO 550 J=1,IGND
!                             HMFGROW(I,J)=HMFGROW(I,J)+HMFGROT(I,M,J)*FAREROT(I,M)
!                             HTCROW(I,J)=HTCROW(I,J)+HTCROT(I,M,J)*FAREROT(I,M)
!                             QFCROW(I,J)=QFCROW(I,J)+QFCROT(I,M,J)*FAREROT(I,M)
!                             GFLXROW(I,J)=GFLXROW(I,J)+GFLXROT(I,M,J)*FAREROT(I,M)
! 550                     CONTINUE
! 575                 CONTINUE
! 600             CONTINUE
!
!             endif ! dodayoutput, for diagnostic fields
!
!             if(dodayoutput) then ! stand alone mode, includes daily output for class
!
!                 !     * ACCUMULATE OUTPUT DATA FOR DIURNALLY AVERAGED FIELDS. BOTH GRID
!                 !       MEAN AND MOSAIC MEAN
!                 !
!                 DO 675 I=1,NLTEST
!
!                     IF (FSSROW(I) .gt. 0.) then
!                         ALTOTACC(I)=ALTOTACC(I) + (FSSROW(I)-(FSGVROW(I)&
!                             &                   +FSGSROW(I)+FSGGROW(I)))/FSSROW(I)
!                         altotcntr_d(i)=altotcntr_d(i) + 1
!                     END IF
!
!                     DO 650 M=1,NMTEST
!                         PREACC(I)=PREACC(I)+PREROW(I)*FAREROT(I,M)*DELT
!                         GTACC(I)=GTACC(I)+GTROT(I,M)*FAREROT(I,M)
!                         QEVPACC(I)=QEVPACC(I)+QEVPROT(I,M)*FAREROT(I,M)
!                         EVAPACC(I)=EVAPACC(I)+QFSROT(I,M)*FAREROT(I,M)*DELT
!                         HFSACC(I)=HFSACC(I)+HFSROT(I,M)*FAREROT(I,M)
!                         HMFNACC(I)=HMFNACC(I)+HMFNROT(I,M)*FAREROT(I,M)
!                         ROFACC(I)=ROFACC(I)+ROFROT(I,M)*FAREROT(I,M)*DELT
!                         OVRACC(I)=OVRACC(I)+ROFOROT(I,M)*FAREROT(I,M)*DELT
!                         WTBLACC(I)=WTBLACC(I)+wtableROT(I,M)*FAREROT(I,M)
!                         IF (FSSROW(I) .gt. 0.) then
!                             ALTOTACC(I)=ALTOTACC(I) + (FSSROW(I)-(FSGVROW(I)&
!                                         +FSGSROW(I)+FSGGROW(I)))/FSSROW(I)
!                         END IF
!                         DO 625 J=1,IGND
!                             TBARACC(I,J)=TBARACC(I,J)+TBARROT(I,M,J)*FAREROT(I,M)
!                             THLQACC(I,J)=THLQACC(I,J)+THLQROT(I,M,J)*FAREROT(I,M)
!                             THICACC(I,J)=THICACC(I,J)+THICROT(I,M,J)*FAREROT(I,M)
!                             THALACC(I,J)=THALACC(I,J)+(THLQROT(I,M,J)+THICROT(I,M,J))&
!                                 &                    *FAREROT(I,M)
! 625                     CONTINUE
!                         ALVSACC(I)=ALVSACC(I)+ALVSROT(I,M)*FAREROT(I,M)*FSVHROW(I)
!                         ALIRACC(I)=ALIRACC(I)+ALIRROT(I,M)*FAREROT(I,M)*FSIHROW(I)
!                         IF(SNOROT(I,M).GT.0.0) THEN
!                             RHOSACC(I)=RHOSACC(I)+RHOSROT(I,M)*FAREROT(I,M)
!                             TSNOACC(I)=TSNOACC(I)+TSNOROT(I,M)*FAREROT(I,M)
!                             WSNOACC(I)=WSNOACC(I)+WSNOROT(I,M)*FAREROT(I,M)
!                             SNOARE(I)=SNOARE(I)+FAREROT(I,M)
!                         ENDIF
!                         IF(TCANROT(I,M).GT.0.5) THEN
!                             TCANACC(I)=TCANACC(I)+TCANROT(I,M)*FAREROT(I,M)
!                             CANARE(I)=CANARE(I)+FAREROT(I,M)
!                         ENDIF
!                         SNOACC(I)=SNOACC(I)+SNOROT(I,M)*FAREROT(I,M)
!                         RCANACC(I)=RCANACC(I)+RCANROT(I,M)*FAREROT(I,M)
!                         SCANACC(I)=SCANACC(I)+SCANROT(I,M)*FAREROT(I,M)
!                         GROACC(I)=GROACC(I)+GROROT(I,M)*FAREROT(I,M)
!                         FSINACC(I)=FSINACC(I)+FSSROW(I)*FAREROT(I,M)
!                         FLINACC(I)=FLINACC(I)+FDLROW(I)*FAREROT(I,M)
!                         FLUTACC(I)=FLUTACC(I)+SBC*GTROT(I,M)**4*FAREROT(I,M)
!                         TAACC(I)=TAACC(I)+TAROW(I)*FAREROT(I,M)
!                         UVACC(I)=UVACC(I)+UVROW(I)*FAREROT(I,M)
!                         PRESACC(I)=PRESACC(I)+PRESROW(I)*FAREROT(I,M)
!                         QAACC(I)=QAACC(I)+QAROW(I)*FAREROT(I,M)
!                         G12ACC(I)=G12ACC(I)+G12GRD(I)*FAREROT(I,M)  !YW March 23, 2015
!                         G23ACC(I)=G23ACC(I)+G23GRD(I)*FAREROT(I,M)  !YW March 23, 2015
! 650                 CONTINUE
! 675             CONTINUE
!
!                 ! * CALCULATE AND PRINT DAILY AVERAGES.
!
!                 IF(NCOUNT.EQ.NDAY) THEN
!
!                     DO 800 I=1,NLTEST
!                         PREACC(I)=PREACC(I)
!                         GTACC(I)=GTACC(I)/REAL(NDAY)
!                         QEVPACC(I)=QEVPACC(I)/REAL(NDAY)
!                         EVAPACC(I)=EVAPACC(I)
!                         HFSACC(I)=HFSACC(I)/REAL(NDAY)
!                         HMFNACC(I)=HMFNACC(I)/REAL(NDAY)
!                         ROFACC(I)=ROFACC(I)
!                         OVRACC(I)=OVRACC(I)
!                         WTBLACC(I)=WTBLACC(I)/REAL(NDAY)
!                         DO 725 J=1,IGND
!                             TBARACC(I,J)=TBARACC(I,J)/REAL(NDAY)
!                             THLQACC(I,J)=THLQACC(I,J)/REAL(NDAY)
!                             THICACC(I,J)=THICACC(I,J)/REAL(NDAY)
!                             THALACC(I,J)=THALACC(I,J)/REAL(NDAY)
! 725                     CONTINUE
!                         IF(FSINACC(I).GT.0.0) THEN
!                             ALVSACC(I)=ALVSACC(I)/(FSINACC(I)*0.5)
!                             ALIRACC(I)=ALIRACC(I)/(FSINACC(I)*0.5)
!                         ELSE
!                             ALVSACC(I)=0.0
!                             ALIRACC(I)=0.0
!                         ENDIF
!                         IF(SNOARE(I).GT.0.0) THEN
!                             RHOSACC(I)=RHOSACC(I)/SNOARE(I)
!                             TSNOACC(I)=TSNOACC(I)/SNOARE(I)
!                             WSNOACC(I)=WSNOACC(I)/SNOARE(I)
!                         ENDIF
!                         IF(CANARE(I).GT.0.0) THEN
!                             TCANACC(I)=TCANACC(I)/CANARE(I)
!                         ENDIF
!                         SNOACC(I)=SNOACC(I)/REAL(NDAY)
!                         RCANACC(I)=RCANACC(I)/REAL(NDAY)
!                         SCANACC(I)=SCANACC(I)/REAL(NDAY)
!                         GROACC(I)=GROACC(I)/REAL(NDAY)
!                         FSINACC(I)=FSINACC(I)/REAL(NDAY)
!                         FLINACC(I)=FLINACC(I)/REAL(NDAY)
!                         FLUTACC(I)=FLUTACC(I)/REAL(NDAY)
!                         TAACC(I)=TAACC(I)/REAL(NDAY)
!                         UVACC(I)=UVACC(I)/REAL(NDAY)
!                         PRESACC(I)=PRESACC(I)/REAL(NDAY)
!                         QAACC(I)=QAACC(I)/REAL(NDAY)
!                         if (altotcntr_d(i) > 0) then
!                             ALTOTACC(I)=ALTOTACC(I)/REAL(altotcntr_d(i))
!                         else
!                             ALTOTACC(I)=0.
!                         end if
!                         FSSTAR=FSINACC(I)*(1.-ALTOTACC(I))
!                         FLSTAR=FLINACC(I)-FLUTACC(I)
!                         QH=HFSACC(I)
!                         QE=QEVPACC(I)
!                         BEG=FSSTAR+FLSTAR-QH-QE
!                         SNOMLT=HMFNACC(I)
!                         IF(RHOSACC(I).GT.0.0) THEN
!                             ZSN=SNOACC(I)/RHOSACC(I)
!                         ELSE
!                             ZSN=0.0
!                         ENDIF
!                         IF(TCANACC(I).GT.0.01) THEN
!                             TCN=TCANACC(I)-TFREZ
!                         ELSE
!                             TCN=0.0
!                         ENDIF
!                         IF(TSNOACC(I).GT.0.01) THEN
!                             TSN=TSNOACC(I)-TFREZ
!                         ELSE
!                             TSN=0.0
!                         ENDIF
!                         GTOUT=GTACC(I)-TFREZ
!
!                         if ((iyear .ge. jdsty) .and. (iyear .le. jdendy)) then
!                             if ((iday .ge. jdstd) .and. (iday .le. jdendd)) then
!
!                                 WRITE(61,6100) IDAY,IYEAR,FSSTAR,FLSTAR,QH,QE,SNOMLT,&
!                                     &                       BEG,GTOUT,SNOACC(I),RHOSACC(I),&
!                                     &                       WSNOACC(I),ALTOTACC(I),ROFACC(I),CUMSNO
!                                 IF(IGND.GT.3) THEN
!                                     WRITE(62,6201) IDAY,IYEAR,(TBARACC(I,J)-TFREZ,&
!                                         &                       THLQACC(I,J),THICACC(I,J),J=1,IGND),&
!                                         &                       ACTLYR_G(I),FTABLE_g(I)
!                                 ELSE
!                                     WRITE(62,6200) IDAY,IYEAR,(TBARACC(I,J)-TFREZ,&
!                                         &                       THLQACC(I,J),THICACC(I,J),J=1,3),&
!                                         &                       TCN,RCANACC(I),SCANACC(I),TSN,ZSN,&
!                                         &                       ACTLYR_G(I),FTABLE_g(I)
!                                 ENDIF
!                                 WRITE(63,6300) IDAY,IYEAR,FSINACC(I),FLINACC(I),&
!                                     &                       TAACC(I)-TFREZ,UVACC(I),PRESACC(I),&
!                                     &                       QAACC(I),PREACC(I),EVAPACC(I)
!
!                             endif
!                         ENDIF
!
!                         !    ----peatland output-----------------------------------------------\
!
!                         write(99,6999)  IDAY,IYEAR,WTBLACC(i), ZSN,PREACC(i),EVAPACC(i),ROFACC(i),g12acc(i),g23acc(i)
! 6999    format(1X,I4,I5,10f12.3)
!                         !    ----YW March 23, 2015 --------------------------------------------/
!
!                         !* RESET ACCUMULATOR ARRAYS.

! FLAG! Here call : resetAccVars from class_statevars. It covers down to the 800 continue line.
!
!                         PREACC(I)=0.
!                         GTACC(I)=0.
!                         QEVPACC(I)=0.
!                         HFSACC(I)=0.
!                         HMFNACC(I)=0.
!                         ROFACC(I)=0.
!                         SNOACC(I)=0.
!                         CANARE(I)=0.
!                         SNOARE(I)=0.
!                         OVRACC(I)=0.
!                         WTBLACC(I)=0.
!                         DO 750 J=1,IGND
!                             TBARACC(I,J)=0.
!                             THLQACC(I,J)=0.
!                             THICACC(I,J)=0.
!                             THALACC(I,J)=0.
! 750                     CONTINUE
!                         ALVSACC(I)=0.
!                         ALIRACC(I)=0.
!                         RHOSACC(I)=0.
!                         TSNOACC(I)=0.
!                         WSNOACC(I)=0.
!                         TCANACC(I)=0.
!                         RCANACC(I)=0.
!                         SCANACC(I)=0.
!                         GROACC(I)=0.
!                         FSINACC(I)=0.
!                         FLINACC(I)=0.
!                         TAACC(I)=0.
!                         UVACC(I)=0.
!                         PRESACC(I)=0.
!                         QAACC(I)=0.
!                         ALTOTACC(I) = 0.
!                         EVAPACC(I)=0.
!                         FLUTACC(I)=0.
! 800                 CONTINUE
!
!                 ENDIF ! IF(NCOUNT.EQ.NDAY)
!
!                 !     CALCULATE AND PRINT MOSAIC DAILY AVERAGES.
!
!                 !       start -> FLAG JM
!                 DO 676 I=1,NLTEST
!                     DO 658 M=1,NMTEST
!                         PREACC_M(I,M)=PREACC_M(I,M)+PREROW(I)*DELT
!                         GTACC_M(I,M)=GTACC_M(I,M)+GTROT(I,M)
!                         QEVPACC_M(I,M)=QEVPACC_M(I,M)+QEVPROT(I,M)
!                         EVAPACC_M(I,M)=EVAPACC_M(I,M)+QFSROT(I,M)*DELT
!                         HFSACC_M(I,M)=HFSACC_M(I,M)+HFSROT(I,M)
!                         HMFNACC_M(I,M)=HMFNACC_M(I,M)+HMFNROT(I,M)
!                         ROFACC_M(I,M)=ROFACC_M(I,M)+ROFROT(I,M)*DELT
!                         OVRACC_M(I,M)=OVRACC_M(I,M)+ROFOROT(I,M)*DELT
!                         WTBLACC_M(I,M)=WTBLACC_M(I,M)+wtableROT(I,M)
!                         DO 626 J=1,IGND
!                             TBARACC_M(I,M,J)=TBARACC_M(I,M,J)+TBARROT(I,M,J)
!                             THLQACC_M(I,M,J)=THLQACC_M(I,M,J)+THLQROT(I,M,J)
!                             THICACC_M(I,M,J)=THICACC_M(I,M,J)+THICROT(I,M,J)
!                             THALACC_M(I,M,J)=THALACC_M(I,M,J)+(THLQROT(I,M,J)+THICROT(I,M,J))
! 626                     CONTINUE
!                         ALVSACC_M(I,M)=ALVSACC_M(I,M)+ALVSROT(I,M)*FSVHROW(I)
!                         ALIRACC_M(I,M)=ALIRACC_M(I,M)+ALIRROT(I,M)*FSIHROW(I)
!                         IF(SNOROT(I,M).GT.0.0) THEN
!                             RHOSACC_M(I,M)=RHOSACC_M(I,M)+RHOSROT(I,M)
!                             TSNOACC_M(I,M)=TSNOACC_M(I,M)+TSNOROT(I,M)
!                             WSNOACC_M(I,M)=WSNOACC_M(I,M)+WSNOROT(I,M)
!                             SNOARE_M(I,M) = SNOARE_M(I,M) + 1.0 !FLAG test.
!                         ENDIF
!                         IF(TCANROT(I,M).GT.0.5) THEN
!                             TCANACC_M(I,M)=TCANACC_M(I,M)+TCANROT(I,M)
!                         !              CANARE(I)=CANARE(I)+FAREROT(I,M)
!                         ENDIF
!                         SNOACC_M(I,M)=SNOACC_M(I,M)+SNOROT(I,M)
!                         RCANACC_M(I,M)=RCANACC_M(I,M)+RCANROT(I,M)
!                         SCANACC_M(I,M)=SCANACC_M(I,M)+SCANROT(I,M)
!                         GROACC_M(I,M)=GROACC_M(I,M)+GROROT(I,M)
!                         IF (FSSROW(I) .gt. 0.) THEN ! we will reuse the altotcntr_d counter values so don't need to do again.
!                             ALTOTACC_M(I,M)=ALTOTACC_M(I,M) + (FSSROW(I)-&
!                                 &                    (FSGVROT(I,M)+FSGSROT(I,M)+&
!                                 &                     FSGGROT(I,M)))/FSSROW(I)
!                         END IF
!                         FSINACC_M(I,M)=FSINACC_M(I,M)+FSSROW(I)
!                         FLINACC_M(I,M)=FLINACC_M(I,M)+FDLROW(I)
!                         FLUTACC_M(I,M)=FLUTACC_M(I,M)+SBC*GTROT(I,M)**4
!                         TAACC_M(I,M)=TAACC_M(I,M)+TAROW(I)
!                         UVACC_M(I,M)=UVACC_M(I,M)+UVROW(I)
!                         PRESACC_M(I,M)=PRESACC_M(I,M)+PRESROW(I)
!                         QAACC_M(I,M)=QAACC_M(I,M)+QAROW(I)
! 658                 CONTINUE
! 676             CONTINUE
!
!                 !CALCULATE AND PRINT DAILY AVERAGES.
!
!                 IF(NCOUNT.EQ.NDAY) THEN
!
!                     DO 808 I=1,NLTEST
!                         DO 809 M=1,NMTEST
!                             PREACC_M(I,M)=PREACC_M(I,M)     !became [kg m-2 day-1] instead of [kg m-2 s-1]
!                             GTACC_M(I,M)=GTACC_M(I,M)/REAL(NDAY)
!                             QEVPACC_M(I,M)=QEVPACC_M(I,M)/REAL(NDAY)
!                             EVAPACC_M(I,M)=EVAPACC_M(I,M)   !became [kg m-2 day-1] instead of [kg m-2 s-1]
!                             HFSACC_M(I,M)=HFSACC_M(I,M)/REAL(NDAY)
!                             HMFNACC_M(I,M)=HMFNACC_M(I,M)/REAL(NDAY)
!                             ROFACC_M(I,M)=ROFACC_M(I,M)   !became [kg m-2 day-1] instead of [kg m-2 s-1
!                             OVRACC_M(I,M)=OVRACC_M(I,M)   !became [kg m-2 day-1] instead of [kg m-2 s-1]
!                             WTBLACC_M(I,M)=WTBLACC_M(I,M)/REAL(NDAY)
!                             DO 726 J=1,IGND
!                                 TBARACC_M(I,M,J)=TBARACC_M(I,M,J)/REAL(NDAY)
!                                 THLQACC_M(I,M,J)=THLQACC_M(I,M,J)/REAL(NDAY)
!                                 THICACC_M(I,M,J)=THICACC_M(I,M,J)/REAL(NDAY)
!                                 THALACC_M(I,M,J)=THALACC_M(I,M,J)/REAL(NDAY)
! 726                         CONTINUE
!
!                             IF(FSINACC_M(I,M).GT.0.0) THEN
!                                 ALVSACC_M(I,M)=ALVSACC_M(I,M)/(FSINACC_M(I,M)*0.5)
!                                 ALIRACC_M(I,M)=ALIRACC_M(I,M)/(FSINACC_M(I,M)*0.5)
!                             ELSE
!                                 ALVSACC_M(I,M)=0.0
!                                 ALIRACC_M(I,M)=0.0
!                             ENDIF
!
!                             SNOACC_M(I,M)=SNOACC_M(I,M)/REAL(NDAY)
!                             if (SNOARE_M(I,M) .GT. 0.) THEN
!                                 RHOSACC_M(I,M)=RHOSACC_M(I,M)/SNOARE_M(I,M)
!                                 TSNOACC_M(I,M)=TSNOACC_M(I,M)/SNOARE_M(I,M)
!                                 WSNOACC_M(I,M)=WSNOACC_M(I,M)/SNOARE_M(I,M)
!                             END IF
!                             TCANACC_M(I,M)=TCANACC_M(I,M)/REAL(NDAY)
!                             RCANACC_M(I,M)=RCANACC_M(I,M)/REAL(NDAY)
!                             SCANACC_M(I,M)=SCANACC_M(I,M)/REAL(NDAY)
!                             GROACC_M(I,M)=GROACC_M(I,M)/REAL(NDAY)
!                             FSINACC_M(I,M)=FSINACC_M(I,M)/REAL(NDAY)
!                             FLINACC_M(I,M)=FLINACC_M(I,M)/REAL(NDAY)
!                             FLUTACC_M(I,M)=FLUTACC_M(I,M)/REAL(NDAY)
!                             TAACC_M(I,M)=TAACC_M(I,M)/REAL(NDAY)
!                             UVACC_M(I,M)=UVACC_M(I,M)/REAL(NDAY)
!                             PRESACC_M(I,M)=PRESACC_M(I,M)/REAL(NDAY)
!                             QAACC_M(I,M)=QAACC_M(I,M)/REAL(NDAY)
!                             if (altotcntr_d(i) > 0) then ! altotcntr_d(i) could be 0
!                                 ALTOTACC_M(I,M)=ALTOTACC_M(I,M)/REAL(altotcntr_d(i))
!                             else
!                                 ALTOTACC_M(I,M)=0.
!                             endif
!                             FSSTAR=FSINACC_M(I,M)*(1.-ALTOTACC_M(I,M))
!                             FLSTAR=FLINACC_M(I,M)-FLUTACC_M(I,M)
!                             QH=HFSACC_M(I,M)
!                             QE=QEVPACC_M(I,M)
!                             QEVPACC_M_SAVE(I,M)=QEVPACC_M(I,M)   !FLAG! What is the point of this? JM Apr 12015
!                             BEG=FSSTAR+FLSTAR-QH-QE
!                             SNOMLT=HMFNACC_M(I,M)
!
!                             IF(RHOSACC_M(I,M).GT.0.0) THEN
!                                 ZSN=SNOACC_M(I,M)/RHOSACC_M(I,M)
!                             ELSE
!                                 ZSN=0.0
!                             ENDIF
!
!                             IF(TCANACC_M(I,M).GT.0.01) THEN
!                                 TCN=TCANACC_M(I,M)-TFREZ
!                             ELSE
!                                 TCN=0.0
!                             ENDIF
!
!                             IF(TSNOACC_M(I,M).GT.0.01) THEN
!                                 TSN=TSNOACC_M(I,M)-TFREZ
!                             ELSE
!                                 TSN=0.0
!                             ENDIF
!
!                             GTOUT=GTACC_M(I,M)-TFREZ
!
!                             if ((iyear .ge. jdsty) .and. (iyear .le. jdendy)) then
!                                 if ((iday .ge. jdstd) .and. (iday .le. jdendd)) then
!                                     !
!                                     !         WRITE TO OUTPUT FILES
!                                     !
!                                     WRITE(611,6100) IDAY,IYEAR,FSSTAR,FLSTAR,QH,QE,SNOMLT,&
!                                         &                    BEG,GTOUT,SNOACC_M(I,M),RHOSACC_M(I,M),&
!                                         &                    WSNOACC_M(I,M),ALTOTACC_M(I,M),ROFACC_M(I,M),&
!                                         &                    CUMSNO,' TILE ',M
!                                     IF(IGND.GT.3) THEN
!                                         WRITE(621,6201) IDAY,IYEAR,(TBARACC_M(I,M,J)-TFREZ,&
!                                             &                  THLQACC_M(I,M,J),THICACC_M(I,M,J),J=1,IGND)&
!                                             &                  ,ACTLYR(I,M), FTABLE(I,M),' TILE ',M
!                                     ELSE
!                                         WRITE(621,6200) IDAY,IYEAR,(TBARACC_M(I,M,J)-TFREZ,&
!                                             &                  THLQACC_M(I,M,J),THICACC_M(I,M,J),J=1,3),&
!                                             &                  TCN,RCANACC_M(I,M),SCANACC_M(I,M),TSN,ZSN,&
!                                             &                  ' TILE ',M
!                                     ENDIF
!                                     WRITE(631,6300) IDAY,IYEAR,FSINACC_M(I,M),FLINACC_M(I,M),&
!                                         &                  TAACC_M(I,M)-TFREZ,UVACC_M(I,M),PRESACC_M(I,M),&
!                                         &                  QAACC_M(I,M),PREACC_M(I,M),EVAPACC_M(I,M),&
!                                         &                  ' TILE ',M
!                                 endif
!                             ENDIF ! IF write daily
!
!                             ! INITIALIZATION FOR MOSAIC TILE AND GRID VARIABLES
!
!                             call resetclassaccum(nltest,nmtest)
!
! 809                     CONTINUE
! 808                 CONTINUE
!                 ENDIF ! IF(NCOUNT.EQ.NDAY)
!             ENDIF !  IFdodayoutput
!


    !==============================================================================================================

    !>\ingroup io_driver_class_monthly_aw
    !>@{
    !>Accumulate and write out the monthly physics outputs

    subroutine class_monthly_aw(lonLocalIndex,latLocalIndex,IDAY,IYEAR,NCOUNT,NDAY,SBC,DELT,nltest,nmtest,TFREZ,&
                                ACTLYR,FTABLE,lastDOY)
                           
        use class_statevars, only : class_out,resetclassmon,class_rot
        use ctem_params, only : nmon, monthend, nlat, nmos, ignd
        use outputManager, only : writeOutput1D,refyr

        implicit none

        ! arguments
        integer, intent(in) :: lonLocalIndex,latLocalIndex
        integer, intent(in) :: IDAY
        integer, intent(in) :: IYEAR
        integer, intent(in) :: NCOUNT
        integer, intent(in) :: NDAY
        integer, intent(in) :: lastDOY
        integer, intent(in) :: nltest
        integer, intent(in) :: nmtest
        real, intent(in) :: SBC  !CLASS common block items,
        real, intent(in) :: DELT !CLASS common block items,
        real, intent(in) :: TFREZ !CLASS common block items,
        real, dimension(nlat,nmos), intent(in) :: ACTLYR          ! Active layer depth (m)
        real, dimension(nlat,nmos), intent(in) :: FTABLE          ! Depth to frozen water table (m)

        ! pointers
        real, dimension(:,:,:), pointer :: TBARROT
        real, dimension(:,:,:), pointer :: THLQROT
        real, dimension(:,:,:), pointer :: THICROT
        real, dimension(:,:,:), pointer :: QFCROT
        real, dimension(:), pointer :: FSSROW
        real, dimension(:), pointer :: FDLROW
        real, dimension(:), pointer :: FSVHROW
        real, dimension(:), pointer :: FSIHROW
        real, dimension(:), pointer :: TAROW
        real, dimension(:), pointer :: PREROW
        real, dimension(:,:), pointer :: ALVSROT
        real, dimension(:,:), pointer :: FAREROT
        real, dimension(:,:), pointer :: ALIRROT
        real, dimension(:,:), pointer :: GTROT
        real, dimension(:,:), pointer :: HFSROT
        real, dimension(:,:), pointer :: QEVPROT
        real, dimension(:,:), pointer :: SNOROT
        real, dimension(:,:), pointer :: WSNOROT
        real, dimension(:,:), pointer :: ROFROT
        real, dimension(:,:), pointer :: QFSROT
        real, dimension(:,:), pointer :: QFGROT
        real, dimension(:,:), pointer :: QFNROT
        real, dimension(:,:), pointer :: QFCLROT
        real, dimension(:,:), pointer :: QFCFROT
        real, dimension(:,:), pointer :: FSGVROT           !< Diagnosed net shortwave radiation on vegetation canopy
        real, dimension(:,:), pointer :: FSGSROT           !< Diagnosed net shortwave radiation on ground snow surface
        real, dimension(:,:), pointer :: FSGGROT           !< Diagnosed net shortwave radiation on ground surface
        real, pointer, dimension(:) :: ALVSACC_MO
        real, pointer, dimension(:) :: ALIRACC_MO
        real, pointer, dimension(:) :: FLUTACC_MO
        real, pointer, dimension(:) :: FSINACC_MO
        real, pointer, dimension(:) :: FLINACC_MO
        real, pointer, dimension(:) :: HFSACC_MO
        real, pointer, dimension(:) :: QEVPACC_MO
        real, pointer, dimension(:) :: SNOACC_MO
        real, pointer, dimension(:) :: WSNOACC_MO
        real, pointer, dimension(:) :: ROFACC_MO
        real, pointer, dimension(:) :: PREACC_MO
        real, pointer, dimension(:) :: EVAPACC_MO
        real, pointer, dimension(:) :: TRANSPACC_MO
        real, pointer, dimension(:) :: TAACC_MO
        real, pointer, dimension(:) :: ACTLYR_MO
        real, pointer, dimension(:) :: FTABLE_MO
        real, pointer, dimension(:) :: ACTLYR_MIN_MO
        real, pointer, dimension(:) :: FTABLE_MIN_MO
        real, pointer, dimension(:) :: ACTLYR_MAX_MO
        real, pointer, dimension(:) :: FTABLE_MAX_MO
        real, pointer, dimension(:) :: ALTOTACC_MO
        real, pointer, dimension(:) :: GROUNDEVAP
        real, pointer, dimension(:) :: CANOPYEVAP
        real, pointer, dimension(:,:) :: TBARACC_MO
        real, pointer, dimension(:,:) :: THLQACC_MO
        real, pointer, dimension(:,:) :: THICACC_MO
        integer, pointer, dimension(:) :: altotcntr_m
   
        ! local

        real :: ALTOT_MO
        integer :: NT
        integer :: NDMONTH
        integer :: i,m,j
        integer :: IMONTH
        real :: tovere
        real, dimension(nlat) :: ACTLYR_tmp
        real, dimension(nlat) :: FTABLE_tmp
        real :: FSSTAR_MO
        real :: FLSTAR_MO
        real :: QH_MO
        real :: QE_MO
        real, dimension(1) :: timeStamp

        ! point pointers
        TBARROT         => class_rot%TBARROT
        THLQROT         => class_rot%THLQROT
        THICROT         => class_rot%THICROT
        QFCROT          => class_rot%QFCROT
        ALVSROT         => class_rot%ALVSROT
        FAREROT         => class_rot%FAREROT
        ALIRROT         => class_rot%ALIRROT
        GTROT           => class_rot%GTROT
        HFSROT          => class_rot%HFSROT
        QEVPROT         => class_rot%QEVPROT
        SNOROT          => class_rot%SNOROT
        WSNOROT         => class_rot%WSNOROT
        ROFROT          => class_rot%ROFROT
        QFSROT          => class_rot%QFSROT
        QFGROT          => class_rot%QFGROT
        QFNROT          => class_rot%QFNROT
        QFCLROT         => class_rot%QFCLROT
        QFCFROT         => class_rot%QFCFROT
        FSGVROT         => class_rot%FSGVROT
        FSGSROT         => class_rot%FSGSROT
        FSGGROT         => class_rot%FSGGROT
        FSSROW          => class_rot%FSSROW
        FDLROW          => class_rot%FDLROW
        FSVHROW         => class_rot%FSVHROW
        FSIHROW         => class_rot%FSIHROW
        TAROW           => class_rot%TAROW
        PREROW          => class_rot%PREROW
        ALVSACC_MO        => class_out%ALVSACC_MO
        ALIRACC_MO        => class_out%ALIRACC_MO
        FLUTACC_MO        => class_out%FLUTACC_MO
        FSINACC_MO        => class_out%FSINACC_MO
        FLINACC_MO        => class_out%FLINACC_MO
        HFSACC_MO         => class_out%HFSACC_MO
        QEVPACC_MO        => class_out%QEVPACC_MO
        SNOACC_MO         => class_out%SNOACC_MO
        WSNOACC_MO        => class_out%WSNOACC_MO
        ROFACC_MO         => class_out%ROFACC_MO
        PREACC_MO         => class_out%PREACC_MO
        EVAPACC_MO        => class_out%EVAPACC_MO
        TRANSPACC_MO      => class_out%TRANSPACC_MO
        TAACC_MO          => class_out%TAACC_MO
        TBARACC_MO        => class_out%TBARACC_MO
        THLQACC_MO        => class_out%THLQACC_MO
        THICACC_MO        => class_out%THICACC_MO
        ACTLYR_MO         => class_out%ACTLYR_MO
        FTABLE_MO         => class_out%FTABLE_MO
        ACTLYR_MIN_MO     => class_out%ACTLYR_MIN_MO
        FTABLE_MIN_MO     => class_out%FTABLE_MIN_MO
        ACTLYR_MAX_MO     => class_out%ACTLYR_MAX_MO
        FTABLE_MAX_MO     => class_out%FTABLE_MAX_MO
        GROUNDEVAP        => class_out%GROUNDEVAP
        CANOPYEVAP        => class_out%CANOPYEVAP
        ALTOTACC_MO       => class_out%ALTOTACC_MO
        altotcntr_m       => class_out%altotcntr_m

        ! ------------

        !> Accumulate output data for monthly averaged fields for class grid-mean.
        !> for both parallel mode and stand alone mode

        FSSTAR_MO   =0.0
        FLSTAR_MO   =0.0
        QH_MO       =0.0
        QE_MO       =0.0
        ACTLYR_tmp  =0.0
        FTABLE_tmp  =0.0

        i = 1 ! offline nlat is always 1 so this array position is always 1.
        DO 821 M=1,NMTEST

            ! These are presently not being outputted but the code is kept in place if the need arises.
            !     ALVSACC_MO(I)=ALVSACC_MO(I)+ALVSROT(I,M)*FAREROT(I,M)*FSVHROW(I)
            !     ALIRACC_MO(I)=ALIRACC_MO(I)+ALIRROT(I,M)*FAREROT(I,M)*FSIHROW(I)
            FLUTACC_MO(I)=FLUTACC_MO(I)+SBC*GTROT(I,M)**4*FAREROT(I,M)
            FSINACC_MO(I)=FSINACC_MO(I)+FSSROW(I)*FAREROT(I,M)
            FLINACC_MO(I)=FLINACC_MO(I)+FDLROW(I)*FAREROT(I,M)
            HFSACC_MO(I) =HFSACC_MO(I)+HFSROT(I,M)*FAREROT(I,M)
            QEVPACC_MO(I)=QEVPACC_MO(I)+QEVPROT(I,M)*FAREROT(I,M)
            SNOACC_MO(I) =SNOACC_MO(I)+SNOROT(I,M)*FAREROT(I,M)
            TAACC_MO(I)=TAACC_MO(I)+TAROW(I)*FAREROT(I,M)
            ACTLYR_MO(I) = ACTLYR_MO(I) + ACTLYR(I,M) * FAREROT(I,M)
            FTABLE_MO(I) = FTABLE_MO(I) + FTABLE(I,M) * FAREROT(I,M)
            ACTLYR_tmp(I) = ACTLYR_tmp(I) + ACTLYR(I,M) * FAREROT(I,M)
            FTABLE_tmp(I) = FTABLE_tmp(I) + FTABLE(I,M) * FAREROT(I,M)
            GROUNDEVAP(I)=GROUNDEVAP(I)+(QFGROT(I,M)+QFNROT(I,M))*FAREROT(I,M)*DELT !ground evap includes both evap and sublimation from snow
            CANOPYEVAP(I)=CANOPYEVAP(I)+(QFCLROT(I,M)+QFCFROT(I,M))*FAREROT(I,M)*DELT !canopy evap includes both evap and sublimation

            IF(SNOROT(I,M).GT.0.0) THEN
                WSNOACC_MO(I)=WSNOACC_MO(I)+WSNOROT(I,M)*FAREROT(I,M)
            ENDIF

            ROFACC_MO(I) =ROFACC_MO(I)+ROFROT(I,M)*FAREROT(I,M)*DELT
            PREACC_MO(I) =PREACC_MO(I)+PREROW(I)*FAREROT(I,M)*DELT
            EVAPACC_MO(I)=EVAPACC_MO(I)+QFSROT(I,M)*FAREROT(I,M)*DELT

            IF(FSSROW(I).GT.0.0) THEN
                ALTOTACC_MO(I)=ALTOTACC_MO(I) + ( (FSSROW(I)-(FSGVROT(I,M)+FSGSROT(I,M)+FSGGROT(I,M))) &
                /FSSROW(I) )*FAREROT(I,M)
                altotcntr_m(i) = altotcntr_m(i) + 1
            ENDIF

            DO 823 J=1,IGND
                TBARACC_MO(I,J)=TBARACC_MO(I,J)+TBARROT(I,M,J)*FAREROT(I,M)
                THLQACC_MO(I,J)=THLQACC_MO(I,J)+THLQROT(I,M,J)*FAREROT(I,M)
                THICACC_MO(I,J)=THICACC_MO(I,J)+THICROT(I,M,J)*FAREROT(I,M)
                TRANSPACC_MO(I)=TRANSPACC_MO(I)+QFCROT(I,M,J)*FAREROT(I,M)*DELT
823          CONTINUE

821     CONTINUE

        ! Check if the active layer has become more shallow or deepended.
        ACTLYR_MAX_MO(I) = max(ACTLYR_MAX_MO(I), ACTLYR_tmp(I))
        ACTLYR_MIN_MO(I) = min(ACTLYR_MIN_MO(I), ACTLYR_tmp(I))
        FTABLE_MAX_MO(I) = max(FTABLE_MAX_MO(I), FTABLE_tmp(I))
        FTABLE_MIN_MO(I) = min(FTABLE_MIN_MO(I), FTABLE_tmp(I))

        DO NT=1,NMON
            IF(IDAY.EQ.monthend(NT+1).AND.NCOUNT.EQ.NDAY)THEN
    
                IMONTH=NT
                NDMONTH=(monthend(NT+1)-monthend(NT))*NDAY

                    ! These are presently not being outputted but the code is kept in place if the need arises.
                !             IF(FSINACC_MO(I).GT.0.0) THEN
                !                 ALVSACC_MO(I)=ALVSACC_MO(I)/(FSINACC_MO(I)*0.5)
                !                 ALIRACC_MO(I)=ALIRACC_MO(I)/(FSINACC_MO(I)*0.5)
                !             ELSE
                !                 ALVSACC_MO(I)=0.0
                !                 ALIRACC_MO(I)=0.0
                !             ENDIF

                ! Albedo is only counted when sun is above horizon so it uses its own counter.\

                if (altotcntr_m(i) > 0) then
                    ALTOTACC_MO(I) = ALTOTACC_MO(I)/REAL(altotcntr_m(i))
                else
                    ALTOTACC_MO(I) = 0.
                end if

                FLUTACC_MO(I)=FLUTACC_MO(I)/REAL(NDMONTH)
                FSINACC_MO(I)=FSINACC_MO(I)/REAL(NDMONTH)
                FLINACC_MO(I)=FLINACC_MO(I)/REAL(NDMONTH)
                HFSACC_MO(I) =HFSACC_MO(I)/REAL(NDMONTH)
                QEVPACC_MO(I)=QEVPACC_MO(I)/REAL(NDMONTH)
                SNOACC_MO(I) =SNOACC_MO(I)/REAL(NDMONTH)
                WSNOACC_MO(I)=WSNOACC_MO(I)/REAL(NDMONTH)
                TAACC_MO(I)=TAACC_MO(I)/REAL(NDMONTH)

                ACTLYR_MO(I) = ACTLYR_MO(I)/REAL(NDMONTH)
                FTABLE_MO(I) = FTABLE_MO(I)/REAL(NDMONTH)
                ! The accumulated quantities don't change.
                !ROFACC_MO,PREACC_MO(I),EVAPACC_MO(I),TRANSPACC_MO(I),GROUNDEVAP,CANOPYEVAP
                DO J=1,IGND
                    TBARACC_MO(I,J)=TBARACC_MO(I,J)/REAL(NDMONTH)
                    THLQACC_MO(I,J)=THLQACC_MO(I,J)/REAL(NDMONTH)
                    THICACC_MO(I,J)=THICACC_MO(I,J)/REAL(NDMONTH)
                ENDDO

                FSSTAR_MO=FSINACC_MO(I)*(1.-ALTOTACC_MO(I))
                FLSTAR_MO=FLINACC_MO(I)-FLUTACC_MO(I)
                QH_MO=HFSACC_MO(I)
                QE_MO=QEVPACC_MO(I)

                if (EVAPACC_MO(I) > 0.) then
                    tovere = TRANSPACC_MO(I)/EVAPACC_MO(I)
                else
                    tovere = 0.
                end if


                ! Prepare the timestamp for this month. Take one day off so it is the last day of the month
                ! rather than the first day of the next month.
                timeStamp = (iyear - 1 - refyr) * lastDOY + monthend(imonth+1) - 1

                call writeOutput1D(lonLocalIndex,latLocalIndex,'fsstar_mo' ,timeStamp,'rss', [FSSTAR_MO])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'flstar_mo' ,timeStamp,'rls', [FLSTAR_MO])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'qh_mo'     ,timeStamp,'hfss', [QH_MO])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'qe_mo'     ,timeStamp,'hfls', [QE_MO])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'snoacc_mo' ,timeStamp,'snw', [SNOACC_MO(I)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'wsnoacc_mo',timeStamp,'wsnw', [WSNOACC_MO(I)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'taacc_mo'  ,timeStamp,'ts', [TAACC_MO(I)-TFREZ])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'groundevap',timeStamp,'evspsblsoi', [GROUNDEVAP(I)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'canopyevap',timeStamp,'evspsblveg', [CANOPYEVAP(I)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'rofacc_mo' ,timeStamp,'mrro', [ROFACC_MO(I)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'preacc_mo' ,timeStamp,'pr', [PREACC_MO(I)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'evapacc_mo',timeStamp,'evspsbl', [EVAPACC_MO(I)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'transpacc_mo',timeStamp,'tran', [TRANSPACC_MO(I)])

                call writeOutput1D(lonLocalIndex,latLocalIndex,'altotacc_mo',timeStamp,'albs', [ALTOTACC_MO(I)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'tbaracc_mo',timeStamp,'tsl', [TBARACC_MO(I,:)-TFREZ])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'thlqacc_mo',timeStamp,'mrsll', [THLQACC_MO(I,:)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'thicacc_mo',timeStamp,'mrsfl', [THICACC_MO(I,:)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'actlyr_mo',timeStamp,'actlyr', [ACTLYR_MO(I)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'actlyr_max_mo',timeStamp,'actlyrmax', [ACTLYR_MAX_MO(I)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'actlyr_min_mo',timeStamp,'actlyrmin', [ACTLYR_MIN_MO(I)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'ftable_mo',timeStamp,'ftable', [FTABLE_MO(I)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'ftable_max_mo',timeStamp,'ftablemax', [FTABLE_MAX_MO(I)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'ftable_min_mo',timeStamp,'ftablemin', [FTABLE_MIN_MO(I)])
                !tovere

                call resetclassmon(nltest)
          
            END IF ! IF(IDAY.EQ.monthend(NT+1).AND.NCOUNT.EQ.NDAY)
        END DO ! NMON

end subroutine class_monthly_aw
!>@}
!==============================================================================================================

!>\ingroup io_driver_class_annual_aw
!>@{
!>Accumulate and write out the annual physics outputs

subroutine class_annual_aw(lonLocalIndex,latLocalIndex,IDAY,IYEAR,NCOUNT,NDAY,SBC,DELT, &
                            nltest,nmtest,ACTLYR,FTABLE,lastDOY)

    use class_statevars,     only : class_out,resetclassyr,class_rot
    use ctem_params, only : nmon, monthend, nlat, nmos, ignd
    use outputManager, only : writeOutput1D,refyr

    implicit none

    ! arguments
    integer, intent(in) :: lonLocalIndex,latLocalIndex
    integer, intent(in) :: IDAY
    integer, intent(in) :: IYEAR
    integer, intent(in) :: NCOUNT
    integer, intent(in) :: NDAY
    integer, intent(in) :: nltest
    integer, intent(in) :: nmtest
    integer, intent(in) :: lastDOY
    real, intent(in) :: SBC
    real, intent(in) :: DELT
    real, dimension(nlat,nmos), intent(in) :: ACTLYR          ! Active layer depth (m)
    real, dimension(nlat,nmos), intent(in) :: FTABLE          ! Depth to frozen water table (m)

    ! pointers
    real, dimension(:), pointer :: FSSROW
    real, dimension(:), pointer :: FDLROW
    real, dimension(:), pointer :: FSVHROW
    real, dimension(:), pointer :: FSIHROW
    real, dimension(:), pointer :: TAROW
    real, dimension(:), pointer :: PREROW
    real, dimension(:,:), pointer :: ALVSROT
    real, dimension(:,:), pointer :: FAREROT
    real, dimension(:,:), pointer :: ALIRROT
    real, dimension(:,:), pointer :: GTROT
    real, dimension(:,:), pointer :: HFSROT
    real, dimension(:,:), pointer :: QEVPROT
    real, dimension(:,:), pointer :: ROFROT
    real, dimension(:,:), pointer :: QFSROT
    real, dimension(:,:,:), pointer :: QFCROT
    real, dimension(:,:), pointer :: FSGVROT           !< Diagnosed net shortwave radiation on vegetation canopy
    real, dimension(:,:), pointer :: FSGSROT           !< Diagnosed net shortwave radiation on ground snow surface
    real, dimension(:,:), pointer :: FSGGROT           !< Diagnosed net shortwave radiation on ground surface
    integer, pointer, dimension(:) :: altotcntr_yr
    real, pointer, dimension(:) :: ALVSACC_YR
    real, pointer, dimension(:) :: ALIRACC_YR
    real, pointer, dimension(:) :: FLUTACC_YR
    real, pointer, dimension(:) :: FSINACC_YR
    real, pointer, dimension(:) :: FLINACC_YR
    real, pointer, dimension(:) :: HFSACC_YR
    real, pointer, dimension(:) :: QEVPACC_YR
    real, pointer, dimension(:) :: ROFACC_YR
    real, pointer, dimension(:) :: PREACC_YR
    real, pointer, dimension(:) :: EVAPACC_YR
    real, pointer, dimension(:) :: TRANSPACC_YR
    real, pointer, dimension(:) :: TAACC_YR
    real, pointer, dimension(:) :: ACTLYR_YR
    real, pointer, dimension(:) :: FTABLE_YR
    real, pointer, dimension(:) :: ALTOTACC_YR

    !local
    integer :: i,m,j
    real :: tovere
    real :: FSSTAR_YR
    real :: FLSTAR_YR
    real :: QH_YR
    real :: QE_YR
    real, dimension(1) :: timeStamp

    !point pointers
    ALVSROT         => class_rot%ALVSROT
    FAREROT         => class_rot%FAREROT
    ALIRROT         => class_rot%ALIRROT
    GTROT           => class_rot%GTROT
    HFSROT          => class_rot%HFSROT
    QEVPROT         => class_rot%QEVPROT
    ROFROT          => class_rot%ROFROT
    QFSROT          => class_rot%QFSROT
    QFCROT          => class_rot%QFCROT
    FSGVROT         => class_rot%FSGVROT
    FSGSROT         => class_rot%FSGSROT
    FSGGROT         => class_rot%FSGGROT
    FSSROW          => class_rot%FSSROW
    FDLROW          => class_rot%FDLROW
    FSVHROW         => class_rot%FSVHROW
    FSIHROW         => class_rot%FSIHROW
    TAROW           => class_rot%TAROW
    PREROW          => class_rot%PREROW
    ALVSACC_YR        => class_out%ALVSACC_YR
    ALIRACC_YR        => class_out%ALIRACC_YR
    FLUTACC_YR        => class_out%FLUTACC_YR
    FSINACC_YR        => class_out%FSINACC_YR
    FLINACC_YR        => class_out%FLINACC_YR
    HFSACC_YR         => class_out%HFSACC_YR
    QEVPACC_YR        => class_out%QEVPACC_YR
    ROFACC_YR         => class_out%ROFACC_YR
    PREACC_YR         => class_out%PREACC_YR
    EVAPACC_YR        => class_out%EVAPACC_YR
    TRANSPACC_YR      => class_out%TRANSPACC_YR
    TAACC_YR          => class_out%TAACC_YR
    ACTLYR_YR         => class_out%ACTLYR_YR
    FTABLE_YR         => class_out%FTABLE_YR
    ALTOTACC_YR       => class_out%ALTOTACC_YR
    altotcntr_yr      => class_out%altotcntr_yr

    !> Accumulate output data for yearly averaged fields for class grid-mean.
    !> for both parallel mode and stand alone mode
    FSSTAR_YR   =0.0
    FLSTAR_YR   =0.0
    QH_YR       =0.0
    QE_YR       =0.0

    i = 1 ! offline nlat is always 1 so this array position is always 1.
    DO 828 M=1,NMTEST

            ! These are presently not being outputted but the code is kept in place if the need arises.
        !         ALVSACC_YR(I)=ALVSACC_YR(I)+ALVSROT(I,M)*FAREROT(I,M)*FSVHROW(I)
        !         ALIRACC_YR(I)=ALIRACC_YR(I)+ALIRROT(I,M)*FAREROT(I,M)*FSIHROW(I)

        FLUTACC_YR(I)=FLUTACC_YR(I)+SBC*GTROT(I,M)**4*FAREROT(I,M)
        FSINACC_YR(I)=FSINACC_YR(I)+FSSROW(I)*FAREROT(I,M)
        FLINACC_YR(I)=FLINACC_YR(I)+FDLROW(I)*FAREROT(I,M)
        HFSACC_YR(I) =HFSACC_YR(I)+HFSROT(I,M)*FAREROT(I,M)
        QEVPACC_YR(I)=QEVPACC_YR(I)+QEVPROT(I,M)*FAREROT(I,M)
        TAACC_YR(I)=TAACC_YR(I)+TAROW(I)*FAREROT(I,M)
        ROFACC_YR(I) =ROFACC_YR(I)+ROFROT(I,M)*FAREROT(I,M)*DELT
        PREACC_YR(I) =PREACC_YR(I)+PREROW(I)*FAREROT(I,M)*DELT
        EVAPACC_YR(I)=EVAPACC_YR(I)+QFSROT(I,M)*FAREROT(I,M)*DELT
        DO J = 1,IGND
            TRANSPACC_YR(I)=TRANSPACC_YR(I)+QFCROT(I,M,J)*FAREROT(I,M)*DELT
        END DO

        !ACTLYR_MO(I) = ACTLYR_MO(I) + ACTLYR(I,M) * FAREROT(I,M)
        !FTABLE_MO(I) = FTABLE_MO(I) + FTABLE(I,M) * FAREROT(I,M)

        IF(FSSROW(I).GT.0.0) THEN
            ALTOTACC_YR(I)=ALTOTACC_YR(I) + ((FSSROW(I)-(FSGVROT(I,M)+FSGSROT(I,M)+FSGGROT(I,M))) &
            /FSSROW(I) )*FAREROT(I,M)
            altotcntr_yr(i) = altotcntr_yr(i) + 1
        ENDIF

828     CONTINUE

    IF (IDAY .EQ. lastDOY .AND. NCOUNT .EQ. NDAY)  THEN

        ! These are presently not being outputted but the code is kept in place if the need arises.
        !             IF(FSINACC_YR(I).GT.0.0) THEN
        !                 ALVSACC_YR(I)=ALVSACC_YR(I)/(FSINACC_YR(I)*0.5)
        !                 ALIRACC_YR(I)=ALIRACC_YR(I)/(FSINACC_YR(I)*0.5)
        !             ELSE
        !                 ALVSACC_YR(I)=0.0
        !                 ALIRACC_YR(I)=0.0
        !             ENDIF

        FLUTACC_YR(I)=FLUTACC_YR(I)/(REAL(NDAY)*real(lastDOY))
        FSINACC_YR(I)=FSINACC_YR(I)/(REAL(NDAY)*real(lastDOY))
        FLINACC_YR(I)=FLINACC_YR(I)/(REAL(NDAY)*real(lastDOY))
        HFSACC_YR(I) =HFSACC_YR(I)/(REAL(NDAY)*real(lastDOY))
        QEVPACC_YR(I)=QEVPACC_YR(I)/(REAL(NDAY)*real(lastDOY))
        ROFACC_YR(I) =ROFACC_YR(I)
        PREACC_YR(I) =PREACC_YR(I)
        EVAPACC_YR(I)=EVAPACC_YR(I)
        TRANSPACC_YR(I)=TRANSPACC_YR(I)
        TAACC_YR(I)=TAACC_YR(I)/(REAL(NDAY)*real(lastDOY))

        ! Albedo is only counted when sun is above horizon so it uses its own counter.
        if (altotcntr_yr(i) > 0) then
            ALTOTACC_YR(I)=ALTOTACC_YR(I)/(REAL(altotcntr_yr(i)))
        else
            ALTOTACC_YR(I)= 0.
        end if

        FSSTAR_YR=FSINACC_YR(I)*(1.-ALTOTACC_YR(I))
        FLSTAR_YR=FLINACC_YR(I)-FLUTACC_YR(I)
        QH_YR=HFSACC_YR(I)
        QE_YR=QEVPACC_YR(I)

        if (EVAPACC_YR(I) > 0.) then
            tovere = TRANSPACC_YR(I)/EVAPACC_YR(I)
        else
            tovere = 0.
        end if

        ! Prepare the timestamp for this year
        timeStamp = (iyear - refyr) * lastDOY

        call writeOutput1D(lonLocalIndex,latLocalIndex,'fsstar_yr' ,timeStamp,'rss', [FSSTAR_YR])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'flstar_yr' ,timeStamp,'rls', [FLSTAR_YR])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'qh_yr'     ,timeStamp,'hfss', [QH_YR])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'qe_yr'     ,timeStamp,'hfls', [QE_YR])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'rofacc_yr' ,timeStamp,'mrro', [ROFACC_YR(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'preacc_yr' ,timeStamp,'pr', [PREACC_YR(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'evapacc_yr' ,timeStamp,'evspsbl', [EVAPACC_YR(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'transpacc_yr' ,timeStamp,'tran', [TRANSPACC_YR(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'altotacc_yr' ,timeStamp,'albs', [ALTOTACC_YR(i)])
        ! tovere

        !> ADD INITIALIZTION FOR YEARLY ACCUMULATED ARRAYS
        call resetclassyr(nltest)

    ENDIF !> IDAY.EQ.365/366 .AND. NDAY

end subroutine class_annual_aw
!>@}
!==============================================================================================================

!>\ingroup io_driver_ctem_daily_aw
!>@{
!> Accumulate and write the daily biogeochemical outputs

subroutine ctem_daily_aw(nltest,nmtest,iday,FAREROT,iyear,jdstd,jdsty,jdendd,jdendy,grclarea,onetile_perPFT,ipeatlandrow)

    ! Accumulate and write out the daily CTEM outputs

    ! J. Melton Feb 2016.

    use ctem_statevars,     only : ctem_tile, vrot, c_switch, &
    resetdaily, ctem_grd
    use ctem_params, only : icc,ignd,nmos,iccp1,wtCH4,seed

    implicit none

    ! arguments
    integer, intent(in) :: nltest
    integer, intent(in) :: nmtest
    integer, intent(in) :: iday
    real, intent(in), dimension(:,:) :: FAREROT
    integer, intent(in) :: iyear
    integer, intent(in) :: jdstd
    integer, intent(in) :: jdsty
    integer, intent(in) :: jdendd
    integer, intent(in) :: jdendy
    real, intent(in), dimension(:) :: grclarea
    logical, intent(in) :: onetile_perPFT

    integer, intent(in), dimension(:,:) :: ipeatlandrow

    ! pointers

    logical, pointer :: dofire
    logical, pointer :: lnduseon
    logical, pointer :: PFTCompetition
    logical, pointer :: dowetlands
    logical, pointer :: obswetf
    real, pointer, dimension(:,:,:) :: fcancmxrow
    real, pointer, dimension(:,:,:) :: gppvegrow
    real, pointer, dimension(:,:,:) :: nepvegrow
    real, pointer, dimension(:,:,:) :: nbpvegrow
    real, pointer, dimension(:,:,:) :: nppvegrow
    real, pointer, dimension(:,:,:) :: hetroresvegrow
    real, pointer, dimension(:,:,:) :: autoresvegrow
    real, pointer, dimension(:,:,:) :: litresvegrow
    real, pointer, dimension(:,:,:) :: soilcresvegrow
    real, pointer, dimension(:,:,:) :: rmlvegaccrow
    real, pointer, dimension(:,:,:) :: rmsvegrow
    real, pointer, dimension(:,:,:) :: rmrvegrow
    real, pointer, dimension(:,:,:) :: rgvegrow
    real, pointer, dimension(:,:,:) :: ailcgrow
    real, pointer, dimension(:,:,:) :: emit_co2row
    real, pointer, dimension(:,:,:) :: emit_corow
    real, pointer, dimension(:,:,:) :: emit_ch4row
    real, pointer, dimension(:,:,:) :: emit_nmhcrow
    real, pointer, dimension(:,:,:) :: emit_h2row
    real, pointer, dimension(:,:,:) :: emit_noxrow
    real, pointer, dimension(:,:,:) :: emit_n2orow
    real, pointer, dimension(:,:,:) :: emit_pm25row
    real, pointer, dimension(:,:,:) :: emit_tpmrow
    real, pointer, dimension(:,:,:) :: emit_tcrow
    real, pointer, dimension(:,:,:) :: emit_ocrow
    real, pointer, dimension(:,:,:) :: emit_bcrow
    real, pointer, dimension(:,:) :: burnfracrow
    real, pointer, dimension(:,:,:) :: burnvegfrow
    real, pointer, dimension(:,:,:) :: smfuncvegrow
    real, pointer, dimension(:,:,:) :: btermrow
    real, pointer, dimension(:,:) :: ltermrow
    real, pointer, dimension(:,:,:) :: mtermrow
    real, pointer, dimension(:,:) :: lucemcomrow
    real, pointer, dimension(:,:) :: lucltrinrow
    real, pointer, dimension(:,:) :: lucsocinrow
    real, pointer, dimension(:,:) :: ch4wet1row
    real, pointer, dimension(:,:) :: ch4wet2row
    real, pointer, dimension(:,:) :: wetfdynrow
    real, pointer, dimension(:,:) :: ch4dyn1row
    real, pointer, dimension(:,:) :: ch4dyn2row
    real, pointer, dimension(:,:) :: ch4soillsrow
    real, pointer, dimension(:,:,:) :: litrmassrow
    real, pointer, dimension(:,:,:) :: soilcmasrow
    real, pointer, dimension(:,:,:) :: vgbiomas_vegrow
    real, pointer, dimension(:,:,:) :: stemmassrow
    real, pointer, dimension(:,:,:) :: rootmassrow
    real, pointer, dimension(:,:,:) :: gleafmasrow        !
    real, pointer, dimension(:,:,:) :: bleafmasrow        !
    real, pointer, dimension(:,:) :: gavglairow
    real, pointer, dimension(:,:,:) :: slairow
    real, pointer, dimension(:,:,:) :: ailcbrow
    real, pointer, dimension(:,:,:) :: flhrlossrow
    real, pointer, dimension(:,:) :: dstcemls3row
    integer, pointer, dimension(:,:,:) :: lfstatusrow
    real, pointer, dimension(:,:) :: vgbiomasrow
    real, pointer, dimension(:,:) :: gavgltmsrow
    real, pointer, dimension(:,:) :: gavgscmsrow
      
    real, pointer, dimension(:,:) :: leaflitr_t
    real, pointer, dimension(:,:) :: tltrleaf_t
    real, pointer, dimension(:,:) :: tltrstem_t
    real, pointer, dimension(:,:) :: tltrroot_t
    real, pointer, dimension(:,:) :: ailcg_t
    real, pointer, dimension(:,:) :: ailcb_t
    real, pointer, dimension(:,:,:) :: rmatctem_t
    real, pointer, dimension(:,:) :: veghght_t
    real, pointer, dimension(:,:) :: rootdpth_t
    real, pointer, dimension(:,:) :: roottemp_t
    real, pointer, dimension(:,:) :: slai_t
    real, pointer, dimension(:,:) :: afrroot_t
    real, pointer, dimension(:,:) :: afrleaf_t
    real, pointer, dimension(:,:) :: afrstem_t
    real, pointer, dimension(:,:) :: laimaxg_t
    real, pointer, dimension(:,:) :: stemmass_t
    real, pointer, dimension(:,:) :: rootmass_t
    real, pointer, dimension(:,:) :: litrmass_t
    real, pointer, dimension(:,:) :: gleafmas_t
    real, pointer, dimension(:,:) :: bleafmas_t
    real, pointer, dimension(:,:) :: soilcmas_t
    real, pointer, dimension(:,:) :: emit_co2_t
    real, pointer, dimension(:,:) :: emit_co_t
    real, pointer, dimension(:,:) :: emit_ch4_t
    real, pointer, dimension(:,:) :: emit_nmhc_t
    real, pointer, dimension(:,:) :: emit_h2_t
    real, pointer, dimension(:,:) :: emit_nox_t
    real, pointer, dimension(:,:) :: emit_n2o_t
    real, pointer, dimension(:,:) :: emit_pm25_t
    real, pointer, dimension(:,:) :: emit_tpm_t
    real, pointer, dimension(:,:) :: emit_tc_t
    real, pointer, dimension(:,:) :: emit_oc_t
    real, pointer, dimension(:,:) :: emit_bc_t
    real, pointer, dimension(:,:) :: bterm_t
    real, pointer, dimension(:,:) :: mterm_t
    real, pointer, dimension(:,:) :: smfuncveg_t
    real, pointer, dimension(:,:) :: tcanoaccrow_out

    real, pointer, dimension(:,:) :: npprow
    real, pointer, dimension(:,:) :: neprow
    real, pointer, dimension(:,:) :: nbprow
    real, pointer, dimension(:,:) :: gpprow
    real, pointer, dimension(:,:) :: hetroresrow
    real, pointer, dimension(:,:) :: autoresrow
    real, pointer, dimension(:,:) :: soilcresprow
    real, pointer, dimension(:,:) :: rgrow
    real, pointer, dimension(:,:) :: litresrow
    real, pointer, dimension(:,:) :: socresrow

    real, pointer, dimension(:,:) :: nppmossrow
    real, pointer, dimension(:,:) :: armossrow

    real, pointer, dimension(:) :: gpp_g
    real, pointer, dimension(:) :: npp_g
    real, pointer, dimension(:) :: nbp_g
    real, pointer, dimension(:) :: autores_g
    real, pointer, dimension(:) :: socres_g
    real, pointer, dimension(:) :: litres_g
    real, pointer, dimension(:) :: dstcemls3_g
    real, pointer, dimension(:) :: litrfall_g
    real, pointer, dimension(:) :: rml_g
    real, pointer, dimension(:) :: rms_g
    real, pointer, dimension(:) :: rg_g
    real, pointer, dimension(:) :: leaflitr_g
    real, pointer, dimension(:) :: tltrstem_g
    real, pointer, dimension(:) :: tltrroot_g
    real, pointer, dimension(:) :: nep_g
    real, pointer, dimension(:) :: hetrores_g
    real, pointer, dimension(:) :: dstcemls_g
    real, pointer, dimension(:) :: humiftrs_g
    real, pointer, dimension(:) :: rmr_g
    real, pointer, dimension(:) :: tltrleaf_g
    real, pointer, dimension(:) :: gavgltms_g
    real, pointer, dimension(:) :: vgbiomas_g
    real, pointer, dimension(:) :: gavglai_g
    real, pointer, dimension(:) :: gavgscms_g
    real, pointer, dimension(:) :: gleafmas_g
    real, pointer, dimension(:) :: bleafmas_g
    real, pointer, dimension(:) :: stemmass_g
    real, pointer, dimension(:) :: rootmass_g
    real, pointer, dimension(:) :: litrmass_g
    real, pointer, dimension(:) :: soilcmas_g
    real, pointer, dimension(:) :: slai_g
    real, pointer, dimension(:) :: ailcg_g
    real, pointer, dimension(:) :: ailcb_g
    real, pointer, dimension(:) :: veghght_g
    real, pointer, dimension(:) :: rootdpth_g
    real, pointer, dimension(:) :: roottemp_g
    real, pointer, dimension(:) :: totcmass_g
    real, pointer, dimension(:) :: tcanoacc_out_g
    real, pointer, dimension(:) :: burnfrac_g
    real, pointer, dimension(:) :: smfuncveg_g
    real, pointer, dimension(:) :: lucemcom_g
    real, pointer, dimension(:) :: lucltrin_g
    real, pointer, dimension(:) :: lucsocin_g
    real, pointer, dimension(:) :: emit_co2_g
    real, pointer, dimension(:) :: emit_co_g
    real, pointer, dimension(:) :: emit_ch4_g
    real, pointer, dimension(:) :: emit_nmhc_g
    real, pointer, dimension(:) :: emit_h2_g
    real, pointer, dimension(:) :: emit_nox_g
    real, pointer, dimension(:) :: emit_n2o_g
    real, pointer, dimension(:) :: emit_pm25_g
    real, pointer, dimension(:) :: emit_tpm_g
    real, pointer, dimension(:) :: emit_tc_g
    real, pointer, dimension(:) :: emit_oc_g
    real, pointer, dimension(:) :: emit_bc_g
    real, pointer, dimension(:) :: bterm_g
    real, pointer, dimension(:) :: lterm_g
    real, pointer, dimension(:) :: mterm_g
    real, pointer, dimension(:) :: ch4wet1_g
    real, pointer, dimension(:) :: ch4wet2_g
    real, pointer, dimension(:) :: wetfdyn_g
    real, pointer, dimension(:) :: ch4dyn1_g
    real, pointer, dimension(:) :: ch4dyn2_g
    real, pointer, dimension(:) :: ch4soills_g
    real, pointer, dimension(:,:) :: afrleaf_g
    real, pointer, dimension(:,:) :: afrstem_g
    real, pointer, dimension(:,:) :: afrroot_g
    real, pointer, dimension(:,:) :: lfstatus_g
    real, pointer, dimension(:,:) :: rmlvegrow_g
    real, pointer, dimension(:,:) :: anvegrow_g
    real, pointer, dimension(:,:) :: rmatctem_g

    !real, pointer, dimension(:) :: gppmosac_g

    real, pointer, dimension(:,:,:) :: bmasvegrow
    real, pointer, dimension(:,:,:) :: cmasvegcrow
    real, pointer, dimension(:,:,:) :: veghghtrow
    real, pointer, dimension(:,:,:) :: rootdpthrow
    real, pointer, dimension(:,:) :: rmlrow
    real, pointer, dimension(:,:) :: rmsrow
    real, pointer, dimension(:,:,:) :: tltrleafrow
    real, pointer, dimension(:,:,:) :: tltrstemrow
    real, pointer, dimension(:,:,:) :: tltrrootrow
    real, pointer, dimension(:,:,:) :: leaflitrrow
    real, pointer, dimension(:,:,:) :: roottemprow
    real, pointer, dimension(:,:,:) :: afrleafrow
    real, pointer, dimension(:,:,:) :: afrstemrow
    real, pointer, dimension(:,:,:) :: afrrootrow
    real, pointer, dimension(:,:,:) :: wtstatusrow
    real, pointer, dimension(:,:,:) :: ltstatusrow
    real, pointer, dimension(:,:) :: rmrrow

    real, pointer, dimension(:,:,:,:) :: rmatctemrow
    real, pointer, dimension(:,:) :: dstcemlsrow
    real, pointer, dimension(:,:) :: litrfallrow
    real, pointer, dimension(:,:) :: humiftrsrow
   
    ! local
    integer :: i,m,j,nt,k
    real :: barefrac
    real :: sumfare

    ! point pointers

    dofire                => c_switch%dofire
    lnduseon              => c_switch%lnduseon
    PFTCompetition        => c_switch%PFTCompetition
    dowetlands            => c_switch%dowetlands
    obswetf               => c_switch%obswetf

    fcancmxrow        => vrot%fcancmx
    gppvegrow         => vrot%gppveg
    nepvegrow         => vrot%nepveg
    nbpvegrow         => vrot%nbpveg
    nppvegrow         => vrot%nppveg
    hetroresvegrow    => vrot%hetroresveg
    autoresvegrow     => vrot%autoresveg
    litresvegrow      => vrot%litresveg
    soilcresvegrow    => vrot%soilcresveg
    rmlvegaccrow      => vrot%rmlvegacc
    rmsvegrow         => vrot%rmsveg
    rmrvegrow         => vrot%rmrveg
    rgvegrow          => vrot%rgveg
    ailcgrow          => vrot%ailcg
    emit_co2row       => vrot%emit_co2
    emit_corow        => vrot%emit_co
    emit_ch4row       => vrot%emit_ch4
    emit_nmhcrow      => vrot%emit_nmhc
    emit_h2row        => vrot%emit_h2
    emit_noxrow       => vrot%emit_nox
    emit_n2orow       => vrot%emit_n2o
    emit_pm25row      => vrot%emit_pm25
    emit_tpmrow       => vrot%emit_tpm
    emit_tcrow        => vrot%emit_tc
    emit_ocrow        => vrot%emit_oc
    emit_bcrow        => vrot%emit_bc
    burnfracrow       => vrot%burnfrac
    burnvegfrow       => vrot%burnvegf
    smfuncvegrow      => vrot%smfuncveg
    btermrow          => vrot%bterm
    ltermrow          => vrot%lterm
    mtermrow          => vrot%mterm
    lucemcomrow       => vrot%lucemcom
    lucltrinrow       => vrot%lucltrin
    lucsocinrow       => vrot%lucsocin
    ch4wet1row        => vrot%ch4wet1
    ch4wet2row        => vrot%ch4wet2
    wetfdynrow        => vrot%wetfdyn
    ch4dyn1row        => vrot%ch4dyn1
    ch4dyn2row        => vrot%ch4dyn2
    ch4soillsrow      => vrot%ch4_soills
    litrmassrow       => vrot%litrmass
    soilcmasrow       => vrot%soilcmas
    vgbiomas_vegrow   => vrot%vgbiomas_veg
    stemmassrow       => vrot%stemmass
    rootmassrow       => vrot%rootmass
    flhrlossrow       => vrot%flhrloss
    dstcemls3row      => vrot%dstcemls3
    lfstatusrow       => vrot%lfstatus

    !gppmosac_g        => ctem_tile%gppmosac_g

    tcanoaccrow_out   => vrot%tcanoaccrow_out
    npprow            => vrot%npp
    neprow            => vrot%nep
    nbprow            => vrot%nbp
    gpprow            => vrot%gpp
    hetroresrow       => vrot%hetrores
    autoresrow        => vrot%autores
    soilcresprow      => vrot%soilcresp
    rgrow             => vrot%rg
    litresrow         => vrot%litres
    socresrow         => vrot%socres
    vgbiomasrow       => vrot%vgbiomas
    gavgltmsrow       => vrot%gavgltms
    gavgscmsrow       => vrot%gavgscms

    nppmossrow         => vrot%nppmoss
    armossrow          => vrot%armoss
      
    bmasvegrow        => vrot%bmasveg
    cmasvegcrow       => vrot%cmasvegc
    veghghtrow        => vrot%veghght
    rootdpthrow       => vrot%rootdpth
    rmlrow            => vrot%rml
    rmsrow            => vrot%rms
    tltrleafrow       => vrot%tltrleaf
    tltrstemrow       => vrot%tltrstem
    tltrrootrow       => vrot%tltrroot
    leaflitrrow       => vrot%leaflitr
    roottemprow       => vrot%roottemp
    afrleafrow        => vrot%afrleaf
    afrstemrow        => vrot%afrstem
    afrrootrow        => vrot%afrroot
    wtstatusrow       => vrot%wtstatus
    ltstatusrow       => vrot%ltstatus
    rmrrow            => vrot%rmr
    gleafmasrow       => vrot%gleafmas
    bleafmasrow       => vrot%bleafmas
    gavglairow        => vrot%gavglai
    slairow           => vrot%slai
    ailcbrow          => vrot%ailcb
    flhrlossrow       => vrot%flhrloss
    rmatctemrow       => vrot%rmatctem
    dstcemlsrow       => vrot%dstcemls
    litrfallrow       => vrot%litrfall
    humiftrsrow       => vrot%humiftrs

    leaflitr_t        => ctem_tile%leaflitr_t
    tltrleaf_t        => ctem_tile%tltrleaf_t
    tltrstem_t        => ctem_tile%tltrstem_t
    tltrroot_t        => ctem_tile%tltrroot_t
    ailcg_t           => ctem_tile%ailcg_t
    ailcb_t           => ctem_tile%ailcb_t
    rmatctem_t        => ctem_tile%rmatctem_t
    veghght_t         => ctem_tile%veghght_t
    rootdpth_t        => ctem_tile%rootdpth_t
    roottemp_t        => ctem_tile%roottemp_t
    slai_t            => ctem_tile%slai_t
    afrroot_t         => ctem_tile%afrroot_t
    afrleaf_t         => ctem_tile%afrleaf_t
    afrstem_t         => ctem_tile%afrstem_t
    laimaxg_t         => ctem_tile%laimaxg_t
    stemmass_t        => ctem_tile%stemmass_t
    rootmass_t        => ctem_tile%rootmass_t
    litrmass_t        => ctem_tile%litrmass_t
    gleafmas_t        => ctem_tile%gleafmas_t
    bleafmas_t        => ctem_tile%bleafmas_t
    soilcmas_t        => ctem_tile%soilcmas_t
    emit_co2_t        => ctem_tile%emit_co2_t
    emit_co_t         => ctem_tile%emit_co_t
    emit_ch4_t        => ctem_tile%emit_ch4_t
    emit_nmhc_t       => ctem_tile%emit_nmhc_t
    emit_h2_t         => ctem_tile%emit_h2_t
    emit_nox_t        => ctem_tile%emit_nox_t
    emit_n2o_t        => ctem_tile%emit_n2o_t
    emit_pm25_t       => ctem_tile%emit_pm25_t
    emit_tpm_t        => ctem_tile%emit_tpm_t
    emit_tc_t         => ctem_tile%emit_tc_t
    emit_oc_t         => ctem_tile%emit_oc_t
    emit_bc_t         => ctem_tile%emit_bc_t
    bterm_t           => ctem_tile%bterm_t
    mterm_t           => ctem_tile%mterm_t
    smfuncveg_t       => ctem_tile%smfuncveg_t

    gpp_g             => ctem_grd%gpp_g
    npp_g             => ctem_grd%npp_g
    nbp_g             => ctem_grd%nbp_g
    autores_g         => ctem_grd%autores_g
    socres_g          => ctem_grd%socres_g
    litres_g          => ctem_grd%litres_g
    dstcemls3_g       => ctem_grd%dstcemls3_g
    litrfall_g        => ctem_grd%litrfall_g
    rml_g             => ctem_grd%rml_g
    rms_g             => ctem_grd%rms_g
    rg_g              => ctem_grd%rg_g
    leaflitr_g        => ctem_grd%leaflitr_g
    tltrstem_g        => ctem_grd%tltrstem_g
    tltrroot_g        => ctem_grd%tltrroot_g
    nep_g             => ctem_grd%nep_g
    hetrores_g        => ctem_grd%hetrores_g
    dstcemls_g        => ctem_grd%dstcemls_g
    humiftrs_g        => ctem_grd%humiftrs_g
    rmr_g             => ctem_grd%rmr_g
    tltrleaf_g        => ctem_grd%tltrleaf_g
    gavgltms_g        => ctem_grd%gavgltms_g
    vgbiomas_g        => ctem_grd%vgbiomas_g
    gavglai_g         => ctem_grd%gavglai_g
    gavgscms_g        => ctem_grd%gavgscms_g
    gleafmas_g        => ctem_grd%gleafmas_g
    bleafmas_g        => ctem_grd%bleafmas_g
    stemmass_g        => ctem_grd%stemmass_g
    rootmass_g        => ctem_grd%rootmass_g
    litrmass_g        => ctem_grd%litrmass_g
    soilcmas_g        => ctem_grd%soilcmas_g
    slai_g            => ctem_grd%slai_g
    ailcg_g           => ctem_grd%ailcg_g
    ailcb_g           => ctem_grd%ailcb_g
    veghght_g         => ctem_grd%veghght_g
    rootdpth_g        => ctem_grd%rootdpth_g
    roottemp_g        => ctem_grd%roottemp_g
    totcmass_g        => ctem_grd%totcmass_g
    tcanoacc_out_g    => ctem_grd%tcanoacc_out_g
    burnfrac_g        => ctem_grd%burnfrac_g
    smfuncveg_g       => ctem_grd%smfuncveg_g
    lucemcom_g        => ctem_grd%lucemcom_g
    lucltrin_g        => ctem_grd%lucltrin_g
    lucsocin_g        => ctem_grd%lucsocin_g
    emit_co2_g        => ctem_grd%emit_co2_g
    emit_co_g         => ctem_grd%emit_co_g
    emit_ch4_g        => ctem_grd%emit_ch4_g
    emit_nmhc_g       => ctem_grd%emit_nmhc_g
    emit_h2_g         => ctem_grd%emit_h2_g
    emit_nox_g        => ctem_grd%emit_nox_g
    emit_n2o_g        => ctem_grd%emit_n2o_g
    emit_pm25_g       => ctem_grd%emit_pm25_g
    emit_tpm_g        => ctem_grd%emit_tpm_g
    emit_tc_g         => ctem_grd%emit_tc_g
    emit_oc_g         => ctem_grd%emit_oc_g
    emit_bc_g         => ctem_grd%emit_bc_g
    bterm_g           => ctem_grd%bterm_g
    lterm_g           => ctem_grd%lterm_g
    mterm_g           => ctem_grd%mterm_g
    ch4wet1_g         => ctem_grd%ch4wet1_g
    ch4wet2_g         => ctem_grd%ch4wet2_g
    wetfdyn_g         => ctem_grd%wetfdyn_g
    ch4dyn1_g         => ctem_grd%ch4dyn1_g
    ch4dyn2_g         => ctem_grd%ch4dyn2_g
    ch4soills_g       => ctem_grd%ch4_soills_g
    afrleaf_g         => ctem_grd%afrleaf_g
    afrstem_g         => ctem_grd%afrstem_g
    afrroot_g         => ctem_grd%afrroot_g
    lfstatus_g        => ctem_grd%lfstatus_g
    rmlvegrow_g       => ctem_grd%rmlvegrow_g
    anvegrow_g        => ctem_grd%anvegrow_g
    rmatctem_g        => ctem_grd%rmatctem_g


    !       ---------------------------------------------------------

    !>write daily ctem results
    if ((iyear .ge. jdsty).and.(iyear.le.jdendy))then
        if ((iday .ge. jdstd).and.(iday .le.jdendd))then

            !>Reset the grid and tile average variables.
            call resetdaily(nltest,nmtest)

            !>First some unit conversions:

            do 10 i = 1,nltest
                do 20 m = 1 , nmtest
                    !   ------convert peatland C fluxes to gC/m2/day for output-----------\
                    nppmossrow(i,m)=nppmossrow(i,m)*1.0377 ! convert to gc/m2.day
                    armossrow(i,m)=armossrow(i,m)*1.0377 ! convert to gc/m2.day

                    do 30 j=1,icc
                        if (fcancmxrow(i,m,j) .gt.0.0) then

                            gppvegrow(i,m,j)=gppvegrow(i,m,j)*1.0377 ! convert to gc/m2.day
                            nppvegrow(i,m,j)=nppvegrow(i,m,j)*1.0377 ! convert to gc/m2.day
                            nepvegrow(i,m,j)=nepvegrow(i,m,j)*1.0377 ! convert to gc/m2.day
                            nbpvegrow(i,m,j)=nbpvegrow(i,m,j)*1.0377 ! convert to gc/m2.day
                            hetroresvegrow(i,m,j)=hetroresvegrow(i,m,j)*1.0377 ! convert to gc/m2.day
                            autoresvegrow(i,m,j)=autoresvegrow(i,m,j)*1.0377 ! convert to gc/m2.day
                            litresvegrow(i,m,j)=litresvegrow(i,m,j)*1.0377 ! convert to gc/m2.day
                            soilcresvegrow(i,m,j)=soilcresvegrow(i,m,j)*1.0377 ! convert to gc/m2.day

                        end if

30                  continue ! icc

                    !>Now for the bare fraction of the grid cell.
                    hetroresvegrow(i,m,iccp1)=hetroresvegrow(i,m,iccp1)*1.0377 ! convert to gc/m2.day
                    litresvegrow(i,m,iccp1)=litresvegrow(i,m,iccp1)*1.0377 ! convert to gc/m2.day
                    soilcresvegrow(i,m,iccp1)=soilcresvegrow(i,m,iccp1)*1.0377 ! convert to gc/m2.day
                    nepvegrow(i,m,iccp1)=nepvegrow(i,m,iccp1)*1.0377 ! convert to gc/m2.day
                    nbpvegrow(i,m,iccp1)=nbpvegrow(i,m,iccp1)*1.0377 ! convert to gc/m2.day

                    npprow(i,m)     =npprow(i,m)*1.0377 ! convert to gc/m2.day
                    gpprow(i,m)     =gpprow(i,m)*1.0377 ! convert to gc/m2.day
                    neprow(i,m)     =neprow(i,m)*1.0377 ! convert to gc/m2.day
                    nbprow(i,m)     =nbprow(i,m)*1.0377 ! convert to gc/m2.day
                    lucemcomrow(i,m)=lucemcomrow(i,m)*1.0377 ! convert to gc/m2.day
                    lucltrinrow(i,m)=lucltrinrow(i,m)*1.0377 ! convert to gc/m2.day
                    lucsocinrow(i,m)=lucsocinrow(i,m)*1.0377 ! convert to gc/m2.day
                    hetroresrow(i,m)=hetroresrow(i,m)*1.0377 ! convert to gc/m2.day
                    autoresrow(i,m) =autoresrow(i,m)*1.0377  ! convert to gc/m2.day
                    litresrow(i,m)  =litresrow(i,m)*1.0377   ! convert to gc/m2.day
                    socresrow(i,m)  =socresrow(i,m)*1.0377   ! convert to gc/m2.day
                    ch4wet1row(i,m) = ch4wet1row(i,m)*1.0377 * wtCH4 / 12. ! convert from umolch4/m2/s to gch4/m2.day
                    ch4wet2row(i,m) = ch4wet2row(i,m)*1.0377 * wtCH4 / 12. ! convert from umolch4/m2/s to gch4/m2.day
                    ch4dyn1row(i,m) = ch4dyn1row(i,m)*1.0377 * wtCH4 / 12. ! convert from umolch4/m2/s to gch4/m2.day
                    ch4dyn2row(i,m) = ch4dyn2row(i,m)*1.0377 * wtCH4 / 12. ! convert from umolch4/m2/s to gch4/m2.day
                    ch4soillsrow(i,m) = ch4soillsrow(i,m)*1.0377 * wtCH4 / 12. ! convert from umolch4/m2/s to gch4/m2.day

20              continue
10          continue


            !>Aggregate to the tile avg vars:
            do 60 i=1,nltest
                do 70 m=1,nmtest
                    barefrac = 1.0
                    do j=1,icc
                        barefrac = barefrac - fcancmxrow(i,m,j)
                        leaflitr_t(i,m)=leaflitr_t(i,m)+leaflitrrow(i,m,j)*fcancmxrow(i,m,j)
                        tltrleaf_t(i,m)=tltrleaf_t(i,m)+tltrleafrow(i,m,j)*fcancmxrow(i,m,j)
                        tltrstem_t(i,m)=tltrstem_t(i,m)+tltrstemrow(i,m,j)*fcancmxrow(i,m,j)
                        tltrroot_t(i,m)=tltrroot_t(i,m)+tltrrootrow(i,m,j)*fcancmxrow(i,m,j)
                        veghght_t(i,m)=veghght_t(i,m)+veghghtrow(i,m,j)*fcancmxrow(i,m,j)
                        rootdpth_t(i,m)=rootdpth_t(i,m)+rootdpthrow(i,m,j)*fcancmxrow(i,m,j)
                        roottemp_t(i,m)=roottemp_t(i,m)+roottemprow(i,m,j)*fcancmxrow(i,m,j)
                        slai_t(i,m)=slai_t(i,m)+slairow(i,m,j)*fcancmxrow(i,m,j)
                        afrleaf_t(i,m)=afrleaf_t(i,m)+afrleafrow(i,m,j)*fcancmxrow(i,m,j)
                        afrstem_t(i,m)=afrstem_t(i,m)+afrstemrow(i,m,j)*fcancmxrow(i,m,j)
                        afrroot_t(i,m)=afrroot_t(i,m)+afrrootrow(i,m,j)*fcancmxrow(i,m,j)
                        ailcg_t(i,m)=ailcg_t(i,m)+ailcgrow(i,m,j)*fcancmxrow(i,m,j)
                        ailcb_t(i,m)=ailcb_t(i,m)+ailcbrow(i,m,j)*fcancmxrow(i,m,j)
                        gleafmas_t(i,m) = gleafmas_t(i,m) + gleafmasrow(i,m,j)*fcancmxrow(i,m,j)
                        bleafmas_t(i,m) = bleafmas_t(i,m) + bleafmasrow(i,m,j)*fcancmxrow(i,m,j)
                        stemmass_t(i,m) = stemmass_t(i,m) + stemmassrow(i,m,j)*fcancmxrow(i,m,j)
                        rootmass_t(i,m) = rootmass_t(i,m) + rootmassrow(i,m,j)*fcancmxrow(i,m,j)
                        litrmass_t(i,m) = litrmass_t(i,m) + litrmassrow(i,m,j)*fcancmxrow(i,m,j)
                        soilcmas_t(i,m) = soilcmas_t(i,m) + soilcmasrow(i,m,j)*fcancmxrow(i,m,j)
                        emit_co2_t(i,m) =emit_co2_t(i,m)+ emit_co2row(i,m,j)*fcancmxrow(i,m,j)
                        emit_co_t(i,m)  =emit_co_t(i,m) + emit_corow(i,m,j)*fcancmxrow(i,m,j)
                        emit_ch4_t(i,m) =emit_ch4_t(i,m)+ emit_ch4row(i,m,j)*fcancmxrow(i,m,j)
                        emit_nmhc_t(i,m)=emit_nmhc_t(i,m)+emit_nmhcrow(i,m,j)*fcancmxrow(i,m,j)
                        emit_h2_t(i,m)  =emit_h2_t(i,m) + emit_h2row(i,m,j)*fcancmxrow(i,m,j)
                        emit_nox_t(i,m) =emit_nox_t(i,m)+ emit_noxrow(i,m,j)*fcancmxrow(i,m,j)
                        emit_n2o_t(i,m) =emit_n2o_t(i,m)+ emit_n2orow(i,m,j)*fcancmxrow(i,m,j)
                        emit_pm25_t(i,m)=emit_pm25_t(i,m)+emit_pm25row(i,m,j)*fcancmxrow(i,m,j)
                        emit_tpm_t(i,m) =emit_tpm_t(i,m)+ emit_tpmrow(i,m,j)*fcancmxrow(i,m,j)
                        emit_tc_t(i,m)  =emit_tc_t(i,m) + emit_tcrow(i,m,j)*fcancmxrow(i,m,j)
                        emit_oc_t(i,m)  =emit_oc_t(i,m) + emit_ocrow(i,m,j)*fcancmxrow(i,m,j)
                        emit_bc_t(i,m)  =emit_bc_t(i,m) + emit_bcrow(i,m,j)*fcancmxrow(i,m,j)
                        bterm_t(i,m)  =bterm_t(i,m) + btermrow(i,m,j)*fcancmxrow(i,m,j)
                        mterm_t(i,m)  =mterm_t(i,m) + mtermrow(i,m,j)*fcancmxrow(i,m,j)
                        smfuncveg_t(i,m) = smfuncveg_t(i,m) + smfuncvegrow(i,m,j)*fcancmxrow(i,m,j)

                        do k=1,ignd
                            rmatctem_t(i,m,k)=rmatctem_t(i,m,k)+rmatctemrow(i,m,j,k)*fcancmxrow(i,m,j)
                        enddo
                    enddo !icc

                    !>Do the bare ground also:
                    litrmass_t(i,m) = litrmass_t(i,m) + litrmassrow(i,m,iccp1)*barefrac
                    soilcmas_t(i,m) = soilcmas_t(i,m) + soilcmasrow(i,m,iccp1)*barefrac

                    !>Calculation of grid averaged variables

                    gpp_g(i) =gpp_g(i) + gpprow(i,m)*FAREROT(i,m)
                    npp_g(i) =npp_g(i) + npprow(i,m)*FAREROT(i,m)
                    nep_g(i) =nep_g(i) + neprow(i,m)*FAREROT(i,m)
                    nbp_g(i) =nbp_g(i) + nbprow(i,m)*FAREROT(i,m)
                    autores_g(i) =autores_g(i) +autoresrow(i,m)*FAREROT(i,m)
                    hetrores_g(i)=hetrores_g(i)+hetroresrow(i,m)*FAREROT(i,m)
                    litres_g(i) =litres_g(i) + litresrow(i,m)*FAREROT(i,m)
                    socres_g(i) =socres_g(i) + socresrow(i,m)*FAREROT(i,m)
                    dstcemls_g(i)=dstcemls_g(i)+dstcemlsrow(i,m)*FAREROT(i,m)
                    dstcemls3_g(i)=dstcemls3_g(i)+dstcemls3row(i,m)*FAREROT(i,m)
                    litrfall_g(i)=litrfall_g(i)+litrfallrow(i,m)*FAREROT(i,m)
                    humiftrs_g(i)=humiftrs_g(i)+humiftrsrow(i,m)*FAREROT(i,m)
                    rml_g(i) =rml_g(i) + rmlrow(i,m)*FAREROT(i,m)
                    rms_g(i) =rms_g(i) + rmsrow(i,m)*FAREROT(i,m)
                    rmr_g(i) =rmr_g(i) + rmrrow(i,m)*FAREROT(i,m)
                    rg_g(i) =rg_g(i) + rgrow(i,m)*FAREROT(i,m)
                    leaflitr_g(i) = leaflitr_g(i) + leaflitr_t(i,m)*FAREROT(i,m)
                    tltrleaf_g(i) = tltrleaf_g(i) + tltrleaf_t(i,m)*FAREROT(i,m)
                    tltrstem_g(i) = tltrstem_g(i) + tltrstem_t(i,m)*FAREROT(i,m)
                    tltrroot_g(i) = tltrroot_g(i) + tltrroot_t(i,m)*FAREROT(i,m)
                    slai_g(i) = slai_g(i) + slai_t(i,m)*FAREROT(i,m)
                    ailcg_g(i)=ailcg_g(i)+ailcg_t(i,m)*FAREROT(i,m)
                    ailcb_g(i)=ailcb_g(i)+ailcb_t(i,m)*FAREROT(i,m)
                    vgbiomas_g(i) =vgbiomas_g(i) + vgbiomasrow(i,m)*FAREROT(i,m)
                    veghght_g(i) = veghght_g(i) + veghght_t(i,m)*FAREROT(i,m)
                    gavglai_g(i) =gavglai_g(i) + gavglairow(i,m)*FAREROT(i,m)
                    gavgltms_g(i) =gavgltms_g(i) + gavgltmsrow(i,m)*FAREROT(i,m)
                    gavgscms_g(i) =gavgscms_g(i) + gavgscmsrow(i,m)*FAREROT(i,m)
                    tcanoacc_out_g(i) =tcanoacc_out_g(i)+tcanoaccrow_out(i,m)*FAREROT(i,m)
                    totcmass_g(i) =vgbiomas_g(i) + gavgltms_g(i) + gavgscms_g(i)
                    gleafmas_g(i) = gleafmas_g(i) + gleafmas_t(i,m)*FAREROT(i,m)
                    bleafmas_g(i) = bleafmas_g(i) + bleafmas_t(i,m)*FAREROT(i,m)
                    stemmass_g(i) = stemmass_g(i) + stemmass_t(i,m)*FAREROT(i,m)
                    rootmass_g(i) = rootmass_g(i) + rootmass_t(i,m)*FAREROT(i,m)
                    rootdpth_g(i) = rootdpth_g(i) + rootdpth_t(i,m)*FAREROT(i,m)
                    roottemp_g(i) = roottemp_g(i) + roottemp_t(i,m)*FAREROT(i,m)
                    litrmass_g(i) = litrmass_g(i) + litrmass_t(i,m)*FAREROT(i,m)
                    soilcmas_g(i) = soilcmas_g(i) + soilcmas_t(i,m)*FAREROT(i,m)
                    burnfrac_g(i) =burnfrac_g(i)+ burnfracrow(i,m)*FAREROT(i,m)
                    smfuncveg_g(i) =smfuncveg_g(i)+smfuncveg_t(i,m)*FAREROT(i,m)
                    lucemcom_g(i) =lucemcom_g(i)+lucemcomrow(i,m)*FAREROT(i,m)
                    lucltrin_g(i) =lucltrin_g(i)+lucltrinrow(i,m)*FAREROT(i,m)
                    lucsocin_g(i) =lucsocin_g(i)+lucsocinrow(i,m)*FAREROT(i,m)
                    bterm_g(i)    =bterm_g(i)   +bterm_t(i,m)*FAREROT(i,m)
                    lterm_g(i)    =lterm_g(i)   +ltermrow(i,m)*FAREROT(i,m)
                    mterm_g(i)    =mterm_g(i)   +mterm_t(i,m)*FAREROT(i,m)
                    ch4wet1_g(i) = ch4wet1_g(i) + ch4wet1row(i,m)*farerot(i,m)
                    ch4wet2_g(i) = ch4wet2_g(i) + ch4wet2row(i,m)*farerot(i,m)
                    wetfdyn_g(i) = wetfdyn_g(i) + wetfdynrow(i,m)*farerot(i,m)
                    ch4dyn1_g(i) = ch4dyn1_g(i) + ch4dyn1row(i,m)*farerot(i,m)
                    ch4dyn2_g(i) = ch4dyn2_g(i) + ch4dyn2row(i,m)*farerot(i,m)
                    ch4soills_g(i) = ch4soills_g(i) + ch4soillsrow(i,m)*farerot(i,m)
                    emit_co2_g(i) =emit_co2_g(i)+ emit_co2_t(i,m)*FAREROT(i,m)
                    emit_co_g(i)  =emit_co_g(i) + emit_co_t(i,m)*FAREROT(i,m)
                    emit_ch4_g(i) =emit_ch4_g(i)+ emit_ch4_t(i,m)*FAREROT(i,m)
                    emit_nmhc_g(i)=emit_nmhc_g(i)+emit_nmhc_t(i,m)*FAREROT(i,m)
                    emit_h2_g(i)  =emit_h2_g(i) + emit_h2_t(i,m)*FAREROT(i,m)
                    emit_nox_g(i) =emit_nox_g(i)+ emit_nox_t(i,m)*FAREROT(i,m)
                    emit_n2o_g(i) =emit_n2o_g(i)+ emit_n2o_t(i,m)*FAREROT(i,m)
                    emit_pm25_g(i)=emit_pm25_g(i)+emit_pm25_t(i,m)*FAREROT(i,m)
                    emit_tpm_g(i) =emit_tpm_g(i)+ emit_tpm_t(i,m)*FAREROT(i,m)
                    emit_tc_g(i)  =emit_tc_g(i) + emit_tc_t(i,m)*FAREROT(i,m)
                    emit_oc_g(i)  =emit_oc_g(i) + emit_oc_t(i,m)*FAREROT(i,m)
                    emit_bc_g(i)  =emit_bc_g(i) + emit_bc_t(i,m)*FAREROT(i,m)
                    ! nppmoss_g(i)  = nppmoss_g(i) +nppmossrow(i,m)*FAREROT(i,m)
                    ! armoss_g(i)   = armoss_g(i) + armossrow(i,m)*FAREROT(i,m)

                    do k=1,ignd
                        rmatctem_g(i,k)=rmatctem_g(i,k)+rmatctem_t(i,m,k)*FAREROT(i,m)
                    end do


70              continue !nmtest
60          continue !nltest

            !>Write daily ctem results

            do 80 i=1,nltest
                do 90 m=1,nmtest

                    barefrac = 1.0

                    !>First the per PFT values to file .CT01D
                    do j=1,icc

                        if (fcancmxrow(i,m,j) .gt. seed) then

                            barefrac = barefrac - fcancmxrow(i,m,j)

                            !>File: .CT01D
                            write(72,8200)iday,iyear,gppvegrow(i,m,j),nppvegrow(i,m,j), &
                            nepvegrow(i,m,j),nbpvegrow(i,m,j),autoresvegrow(i,m,j), &
                            hetroresvegrow(i,m,j),litresvegrow(i,m,j),soilcresvegrow(i,m,j), &
                            (dstcemlsrow(i,m)+dstcemls3row(i,m)), &   ! FLAG at present dstcemls are only per tile values
                            litrfallrow(i,m),humiftrsrow(i,m), & ! same with litrfall and humiftrs.
                            ' TILE ',m,' PFT ',j,' FRAC ',fcancmxrow(i,m,j)

                            !>File .CT02D
                            write(73,8300)iday,iyear,rmlvegaccrow(i,m,j), &
                            rmsvegrow(i,m,j),rmrvegrow(i,m,j),rgvegrow(i,m,j), &
                            leaflitrrow(i,m,j),tltrleafrow(i,m,j), &
                            tltrstemrow(i,m,j),tltrrootrow(i,m,j), &
                            ' TILE ',m,' PFT ',j,' FRAC ',fcancmxrow(i,m,j)

                            !>File *.CT03D
                            write(74,8401)iday,iyear,vgbiomas_vegrow(i,m,j), &
                            ailcgrow(i,m,j),gleafmasrow(i,m,j), &
                            bleafmasrow(i,m,j), stemmassrow(i,m,j), &
                            rootmassrow(i,m,j), litrmassrow(i,m,j),  &
                            soilcmasrow(i,m,j), &
                            ' TILE ',m,' PFT ',j,' FRAC ',fcancmxrow(i,m,j)

                            !>File .CT04D
                            write(75,8500)iday,iyear, ailcgrow(i,m,j),  &
                            ailcbrow(i,m,j),(rmatctemrow(i,m,j,k),k=1,3), &
                            veghghtrow(i,m,j),rootdpthrow(i,m,j), &
                            roottemprow(i,m,j),slairow(i,m,j), &
                            ' TILE ',m,' PFT ',j,' FRAC ',fcancmxrow(i,m,j)

                            ! File .CT05D
                            !write(76,8600)iday,iyear, afrleafrow(i,m,j),  &
                            !afrstemrow(i,m,j),afrrootrow(i,m,j),  &
                            !tcanoaccrow_out(i,m), lfstatusrow(i,m,j), &
                            !' TILE ',m,' PFT ',j,' FRAC ',fcancmxrow(i,m,j)

                            !>File *.CT06D
                            if (dofire .or. lnduseon) then
                                write(77,8800)iday,iyear, &
                                emit_co2row(i,m,j),emit_corow(i,m,j),emit_ch4row(i,m,j), &
                                emit_nmhcrow(i,m,j),emit_h2row(i,m,j),emit_noxrow(i,m,j), &
                                emit_n2orow(i,m,j),emit_pm25row(i,m,j), &
                                emit_tpmrow(i,m,j),emit_tcrow(i,m,j),emit_ocrow(i,m,j), &
                                emit_bcrow(i,m,j),burnvegfrow(i,m,j)*100., &
                                smfuncvegrow(i,m,j),lucemcom_g(i), &  !FLAG only per grid values for these last ones.
                                lucltrin_g(i), lucsocin_g(i), &
                                grclarea(i), btermrow(i,m,j), lterm_g(i), mtermrow(i,m,j), &
                                ' TILE ',m,' PFT ',j,' FRAC ',fcancmxrow(i,m,j)
                            endif

                        end if !fcancmx

                    end do !icc

                    !>Now write out the bare fraction values ( only needed if you have vars that are affected by barefrac values)
                    if (barefrac .gt. seed) then

                        !>File: .CT01D
                        write(72,8200)iday,iyear,0.0,0.0, &
                        nepvegrow(i,m,iccp1),nbpvegrow(i,m,iccp1),0.0, &
                        hetroresvegrow(i,m,iccp1),litresvegrow(i,m,iccp1),soilcresvegrow(i,m,iccp1), &
                        (dstcemlsrow(i,m)+dstcemls3row(i,m)), &   ! FLAG at present dstcemls are only per tile values
                        litrfallrow(i,m),humiftrsrow(i,m), & ! same with litrfall and humiftrs.
                        ' TILE ',m,' PFT ',iccp1,' FRAC ',barefrac

                        !>File *.CT03D
                        write(74,8401)iday,iyear,0.0, &
                        0.0,0.0, &
                        0.0, 0.0, &
                        0.0, litrmassrow(i,m,iccp1),  &
                        soilcmasrow(i,m,iccp1), &
                        ' TILE ',m,' PFT ',iccp1,' FRAC ',barefrac

                    end if

                    !>Now write out the tile average values for each tile if the tile number
                    !>is greater than 1 (nmtest > 1).
                    if (nmtest > 1) then

                        !>File: .CT01D
                        write(72,8200)iday,iyear,gpprow(i,m),npprow(i,m), &
                        neprow(i,m),nbprow(i,m),autoresrow(i,m), &
                        hetroresrow(i,m),litresrow(i,m),socresrow(i,m), &
                        (dstcemlsrow(i,m)+dstcemls3row(i,m)), &
                        litrfallrow(i,m),humiftrsrow(i,m), &
                        ' TILE ',m,' OF ',nmtest,' TFRAC ',FAREROT(i,m)

                        !>File .CT02D
                        write(73,8300)iday,iyear,rmlrow(i,m),rmsrow(i,m), &
                        rmrrow(i,m),rgrow(i,m),leaflitr_t(i,m),tltrleaf_t(i,m), &
                        tltrstem_t(i,m),tltrroot_t(i,m), &
                        ' TILE ',m,' OF ',nmtest,' TFRAC ',FAREROT(i,m)

                        !>File .CT03D
                        write(74,8401)iday,iyear,vgbiomasrow(i,m), &
                        ailcg_t(i,m), gleafmas_t(i,m), &
                        bleafmas_t(i,m), stemmass_t(i,m), &
                        rootmass_t(i,m), litrmass_t(i,m), &
                        soilcmas_t(i,m),&
                        ' TILE ',m,' OF ',nmtest,' TFRAC ',FAREROT(i,m)

                        !>File .CT04D
                        write(75,8500)iday,iyear,ailcg_t(i,m), &
                        ailcb_t(i,m),(rmatctem_t(i,m,k),k=1,3), &
                        veghght_t(i,m),rootdpth_t(i,m), &
                        roottemp_t(i,m),slai_t(i,m), &
                        ' TILE ',m,' OF ',nmtest,' TFRAC ',FAREROT(i,m)

                        ! File .CT05D
                        !write(76,8601)iday,iyear, afrleaf_t(i,m), &
                        !    afrstem_t(i,m),afrroot_t(i,m),  &
                        !    tcanoaccrow_out(i,m), -999,   & ! lfstatus is kinda meaningless grid avg so set to -999
                        !    ' TILE ',m,' OF ',nmtest,' TFRAC ',FAREROT(i,m)

                        if (dofire .or. lnduseon) then
                            write(77,8800)iday,iyear,  &
                            emit_co2_t(i,m), emit_co_t(i,m), emit_ch4_t(i,m), &
                            emit_nmhc_t(i,m), emit_h2_t(i,m), emit_nox_t(i,m), &
                            emit_n2o_t(i,m), emit_pm25_t(i,m), emit_tpm_t(i,m), &
                            emit_tc_t(i,m), emit_oc_t(i,m), emit_bc_t(i,m), &
                            burnfrac_g(i)*100., smfuncveg_t(i,m),lucemcom_g(i), & !FLAG only per grid values for these last ones.
                            lucltrin_g(i), lucsocin_g(i), &
                            grclarea(i), bterm_t(i,m), lterm_g(i), mterm_t(i,m), &
                            ' TILE ',m,' OF ',nmtest,' TFRAC ',FAREROT(i,m)
                        endif

                    end if !nmtest>1

90              continue !nmtest

                !>Finally do the grid avg values:

                !>File: .CT01D
                write(72,8200)iday,iyear,gpp_g(i),npp_g(i), &
                nep_g(i),nbp_g(i),autores_g(i), &
                hetrores_g(i),litres_g(i),socres_g(i), &
                (dstcemls_g(i)+dstcemls3_g(i)), &
                litrfall_g(i),humiftrs_g(i),' GRDAV'

                !>File .CT02D
                write(73,8300)iday,iyear,rml_g(i),rms_g(i), &
                rmr_g(i),rg_g(i),leaflitr_g(i),tltrleaf_g(i), &
                tltrstem_g(i),tltrroot_g(i),' GRDAV'

                !>File .CT03D
                write(74,8401)iday,iyear,vgbiomas_g(i), &
                gavglai_g(i), &
                gleafmas_g(i), bleafmas_g(i), stemmass_g(i), &
                rootmass_g(i), litrmass_g(i), soilcmas_g(i),' GRDAV'

                !>File .CT04D
                write(75,8500)iday,iyear, ailcg_g(i),  &
                ailcb_g(i),(rmatctem_g(i,k),k=1,3), &
                veghght_g(i),rootdpth_g(i),roottemp_g(i),&
                slai_g(i),' GRDAV'

                ! File .CT05D
                !write(76,8601)iday,iyear, afrleaf_t(i,m), &
                !    afrstem_t(i,m),afrroot_t(i,m),  &
                !    tcanoaccrow_out(i,m), -999,   & ! lfstatus is kinda meaningless grid avg so set to -999
                !    ' GRDAV'

                !>File *.CT06D
                if (dofire .or. lnduseon) then
                    write(77,8800)iday,iyear,  &
                    emit_co2_g(i), emit_co_g(i), emit_ch4_g(i), &
                    emit_nmhc_g(i), emit_h2_g(i), emit_nox_g(i), &
                    emit_n2o_g(i), emit_pm25_g(i), emit_tpm_g(i), &
                    emit_tc_g(i), emit_oc_g(i), emit_bc_g(i), &
                    burnfrac_g(i)*100., smfuncveg_g(i),lucemcom_g(i), &
                    lucltrin_g(i), lucsocin_g(i), &
                    grclarea(i), bterm_g(i), lterm_g(i), mterm_g(i),' GRDAV'
                endif


                !>File .CT08D
                if (dowetlands .or. obswetf) then
                    write(79,8810)iday,iyear, ch4wet1_g(i),  &
                    ch4wet2_g(i), wetfdyn_g(i),  &
                    ch4dyn1_g(i), ch4dyn2_g(i),  &
                    ch4soills_g(i),' GRDAV'
                endif

                ! FLAG FLAG FLAG
                !   ---------------peatland outputs-----------------------------------\
                !   - Note that YW's original code used a mixture of gat/row variables.
                !     The gat variables are replaced with row versions except for
                !     gppmosac_g which has no row version equivalent. Needs to be scattered out?
                !   - Also note that in YW's original code the arrays were hard-coded
                !     for just the 1st tile of the 1st grid point.
                !   - Leave gppmosac_g hard-coded to 1 for now, but should be changed!
                !   EC - Feb 2106.
                ! I think I will just remove these unless they are needed. JM Nov 2016.
                do m=1,nmtest

                    !   CT11D_G   convert moss gpp from umol/m2/s to g/m2/day
                    write (93,6993) iday,iyear, &
                    nppmossrow(i,m),armossrow(i,m), &!,gppmosac_g(1)*1.0377, &
                    (fcancmxrow(i,m,j)*gppvegrow(i,m,j),j=1,icc),      &
                    (fcancmxrow(i,m,j)*nppvegrow(i,m,j),j=1,icc),      &
                    (fcancmxrow(i,m,j)*autoresvegrow(i,m,j),j=1,icc),  &
                    (fcancmxrow(i,m,j)*hetroresvegrow(i,m,j),j=1,icc), &
                    (fcancmxrow(i,m,j),j=1,icc)
                    !   CT12D_G
                    write (94,6993) iday,iyear,(veghghtrow(i,m,j),j=1,icc), &
                    (rootdpthrow(i,m,j),j=1,icc),(ailcgrow(i,m,j),j=1,icc), &
                    (fcancmxrow(i,m,j)*stemmassrow(i,m,j),j=1,icc), &
                    (fcancmxrow(i,m,j)*rootmassrow(i,m,j),j=1,icc), &
                    (fcancmxrow(i,m,j)*litrmassrow(i,m,j),j=1,icc), &
                    (fcancmxrow(i,m,j)*gleafmasrow(i,m,j),j=1,icc), &
                    (fcancmxrow(i,m,j)*bleafmasrow(i,m,j),j=1,icc)

                enddo

6993            format(2i5,100f12.6)

                !   ----------------YW March 27, 2015 -------------------------------/

                if (PFTCompetition .or. lnduseon) then

                    sumfare=0.0
                    if (onetile_perPFT) then
                        do m=1,nmos
                            sumfare=sumfare+FAREROT(i,m)
                        enddo
                        write(78,8200)iday,iyear,(FAREROT(i,m)*100.,m=1,nmos),sumfare
                    else !composite
                        do m=1,nmos
                            sumfare=0.0
                            do j=1,icc  !m = 1
                                sumfare=sumfare+fcancmxrow(i,m,j)
                            enddo !j
                            write(78,8200)iday,iyear,(fcancmxrow(i,m,j)*100.,j=1,icc),(1.0-sumfare)*100.,sumfare,' TILE ',m
                        end do !m
                    endif !mosaic/composite
                endif !PFTCompetition/lnduseon

80          continue ! nltest

        end if !if write daily
    endif !if write daily

8200 format(1x,i4,i5,11f10.5,2(a6,i2),a8,f5.3)
8201 format(1x,i4,i5,3f10.5,80x,2(a6,i2),a8,f5.3)
8300 format(1x,i4,i5,8f10.5,2(a6,i2),a8,f5.3)
8301 format(1x,i4,i5,4f10.5,40x,2(a6,i2),a8,f5.3)
    !8400       format(1x,i4,i5,11f10.5,2(a6,i2),a6,f5.3)
8401 format(1x,i4,i5,2f10.5,6f10.5,2(a6,i2),a8,f5.3)
    !8402       format(1x,i4,i5,10f10.5,2(a6,i2),a6,f5.3)
    ! 8402       format(1x,i4,i5,2f10.5,40x,2f10.5,2(a6,i2))
8500 format(1x,i4,i5,9f10.5,2(a6,i2),a8,f5.3)
8600 format(1x,i4,i5,4f10.5,i8,2(a6,i2),a8,f5.3)
8601 format(1x,i4,i5,4f10.5,8x,2(a6,i2),a8,f5.3)
8800 format(1x,i4,i5,20f11.4,2x,f9.2,2(a6,i2),a8,f5.3)
8810 format(1x,i4,i5,6f11.4,2(a6,i2),a8,f5.3)

end subroutine ctem_daily_aw
!>@}
!==============================================================================================================

!>\ingroup io_driver_ctem_monthly_aw
!>@{
!> Accumulate and write out the monthly CTEM outputs

subroutine ctem_monthly_aw(lonLocalIndex,latLocalIndex,nltest,nmtest,iday,FAREROT,iyear,nday,lastDOY)

    ! J. Melton Feb 2016.

    use ctem_statevars,     only : ctem_tile_mo, vrot, ctem_grd_mo, c_switch, &
                                   resetmonthend,ctem_mo
    use ctem_params, only : icc,iccp1,nmon,mmday,monthend,monthdays,seed,nmos,nlat
    use outputManager, only : writeOutput1D,refyr

    implicit none

    ! arguments
    integer, intent(in) :: lonLocalIndex,latLocalIndex
    integer, intent(in) :: nltest
    integer, intent(in) :: nmtest
    integer, intent(in) :: iday
    real, intent(in), dimension(:,:) :: FAREROT
    integer, intent(in) :: iyear
    integer, intent(in) :: nday
    integer, intent(in) :: lastDOY

    ! pointers

    logical, pointer :: dofire
    logical, pointer :: lnduseon
    logical, pointer :: PFTCompetition
    logical, pointer :: dowetlands
    logical, pointer :: obswetf
    logical, pointer :: doperpftoutput
    logical, pointer :: dopertileoutput

    real, pointer, dimension(:,:,:) :: fcancmxrow
    real, pointer, dimension(:,:,:) :: laimaxg_mo
    real, pointer, dimension(:,:,:) :: stemmass_mo
    real, pointer, dimension(:,:,:) :: rootmass_mo
    real, pointer, dimension(:,:,:) :: litrfallveg_mo
    real, pointer, dimension(:,:,:) :: humiftrsveg_mo
    real, pointer, dimension(:,:,:) :: npp_mo
    real, pointer, dimension(:,:,:) :: gpp_mo
    real, pointer, dimension(:,:,:) :: vgbiomas_mo
    real, pointer, dimension(:,:,:) :: autores_mo
    real, pointer, dimension(:,:,:) :: totcmass_mo
    real, pointer, dimension(:,:,:) :: litrmass_mo
    real, pointer, dimension(:,:,:) :: soilcmas_mo
    real, pointer, dimension(:,:,:) :: nep_mo
    real, pointer, dimension(:,:,:) :: litres_mo
    real, pointer, dimension(:,:,:) :: soilcres_mo
    real, pointer, dimension(:,:,:) :: hetrores_mo
    real, pointer, dimension(:,:,:) :: nbp_mo
    real, pointer, dimension(:,:,:) :: emit_co2_mo
    real, pointer, dimension(:,:,:) :: emit_co_mo
    real, pointer, dimension(:,:,:) :: emit_ch4_mo
    real, pointer, dimension(:,:,:) :: emit_nmhc_mo
    real, pointer, dimension(:,:,:) :: emit_h2_mo
    real, pointer, dimension(:,:,:) :: emit_nox_mo
    real, pointer, dimension(:,:,:) :: emit_n2o_mo
    real, pointer, dimension(:,:,:) :: emit_pm25_mo
    real, pointer, dimension(:,:,:) :: emit_tpm_mo
    real, pointer, dimension(:,:,:) :: emit_tc_mo
    real, pointer, dimension(:,:,:) :: emit_oc_mo
    real, pointer, dimension(:,:,:) :: emit_bc_mo
    real, pointer, dimension(:,:,:) :: bterm_mo
    real, pointer, dimension(:,:,:) :: mterm_mo
    real, pointer, dimension(:,:,:) :: burnfrac_mo
    real, pointer, dimension(:,:,:) :: smfuncveg_mo

    real, pointer, dimension(:,:) :: laimaxg_mo_t
    real, pointer, dimension(:,:) :: stemmass_mo_t
    real, pointer, dimension(:,:) :: rootmass_mo_t
    real, pointer, dimension(:,:) :: litrfall_mo_t
    real, pointer, dimension(:,:) :: humiftrs_mo_t
    real, pointer, dimension(:,:) :: npp_mo_t
    real, pointer, dimension(:,:) :: gpp_mo_t
    real, pointer, dimension(:,:) :: vgbiomas_mo_t
    real, pointer, dimension(:,:) :: autores_mo_t
    real, pointer, dimension(:,:) :: totcmass_mo_t
    real, pointer, dimension(:,:) :: litrmass_mo_t
    real, pointer, dimension(:,:) :: soilcmas_mo_t
    real, pointer, dimension(:,:) :: nep_mo_t
    real, pointer, dimension(:,:) :: litres_mo_t
    real, pointer, dimension(:,:) :: soilcres_mo_t
    real, pointer, dimension(:,:) :: hetrores_mo_t
    real, pointer, dimension(:,:) :: nbp_mo_t
    real, pointer, dimension(:,:) :: emit_co2_mo_t
    real, pointer, dimension(:,:) :: emit_co_mo_t
    real, pointer, dimension(:,:) :: emit_ch4_mo_t
    real, pointer, dimension(:,:) :: emit_nmhc_mo_t
    real, pointer, dimension(:,:) :: emit_h2_mo_t
    real, pointer, dimension(:,:) :: emit_nox_mo_t
    real, pointer, dimension(:,:) :: emit_n2o_mo_t
    real, pointer, dimension(:,:) :: emit_pm25_mo_t
    real, pointer, dimension(:,:) :: emit_tpm_mo_t
    real, pointer, dimension(:,:) :: emit_tc_mo_t
    real, pointer, dimension(:,:) :: emit_oc_mo_t
    real, pointer, dimension(:,:) :: emit_bc_mo_t
    real, pointer, dimension(:,:) :: burnfrac_mo_t
    real, pointer, dimension(:,:) :: smfuncveg_mo_t
    real, pointer, dimension(:,:) :: bterm_mo_t
    real, pointer, dimension(:,:) :: luc_emc_mo_t
    real, pointer, dimension(:,:) :: lterm_mo_t
    real, pointer, dimension(:,:) :: lucsocin_mo_t
    real, pointer, dimension(:,:) :: mterm_mo_t
    real, pointer, dimension(:,:) :: lucltrin_mo_t
    real, pointer, dimension(:,:) :: ch4wet1_mo_t
    real, pointer, dimension(:,:) :: ch4wet2_mo_t
    real, pointer, dimension(:,:) :: wetfdyn_mo_t
    real, pointer, dimension(:,:) :: ch4dyn1_mo_t
    real, pointer, dimension(:,:) :: ch4dyn2_mo_t
    real, pointer, dimension(:,:) :: ch4soills_mo_t
    real, pointer, dimension(:,:) :: wind_mo_t

    logical, pointer, dimension(:,:,:) :: pftexistrow
    real, pointer, dimension(:,:,:) :: gppvegrow
    real, pointer, dimension(:,:,:) :: nepvegrow
    real, pointer, dimension(:,:,:) :: nbpvegrow
    real, pointer, dimension(:,:,:) :: nppvegrow
    real, pointer, dimension(:,:,:) :: hetroresvegrow
    real, pointer, dimension(:,:,:) :: autoresvegrow
    real, pointer, dimension(:,:,:) :: litresvegrow
    real, pointer, dimension(:,:,:) :: soilcresvegrow
    real, pointer, dimension(:,:,:) :: rmlvegaccrow
    real, pointer, dimension(:,:,:) :: rmsvegrow
    real, pointer, dimension(:,:,:) :: rmrvegrow
    real, pointer, dimension(:,:,:) :: rgvegrow
    real, pointer, dimension(:,:,:) :: ailcgrow
    real, pointer, dimension(:,:,:) :: emit_co2row
    real, pointer, dimension(:,:,:) :: emit_corow
    real, pointer, dimension(:,:,:) :: emit_ch4row
    real, pointer, dimension(:,:,:) :: emit_nmhcrow
    real, pointer, dimension(:,:,:) :: emit_h2row
    real, pointer, dimension(:,:,:) :: emit_noxrow
    real, pointer, dimension(:,:,:) :: emit_n2orow
    real, pointer, dimension(:,:,:) :: emit_pm25row
    real, pointer, dimension(:,:,:) :: emit_tpmrow
    real, pointer, dimension(:,:,:) :: emit_tcrow
    real, pointer, dimension(:,:,:) :: emit_ocrow
    real, pointer, dimension(:,:,:) :: emit_bcrow
    real, pointer, dimension(:,:) :: burnfracrow
    real, pointer, dimension(:,:,:) :: burnvegfrow
    real, pointer, dimension(:,:,:) :: smfuncvegrow
    real, pointer, dimension(:,:,:) :: btermrow
    real, pointer, dimension(:,:) :: ltermrow
    real, pointer, dimension(:,:,:) :: mtermrow
    real, pointer, dimension(:,:) :: lucemcomrow
    real, pointer, dimension(:,:) :: lucltrinrow
    real, pointer, dimension(:,:) :: lucsocinrow
    real, pointer, dimension(:,:) :: ch4wet1row
    real, pointer, dimension(:,:) :: ch4wet2row
    real, pointer, dimension(:,:) :: wetfdynrow
    real, pointer, dimension(:,:) :: ch4dyn1row
    real, pointer, dimension(:,:) :: ch4dyn2row
    real, pointer, dimension(:,:) :: ch4soillsrow
    real, pointer, dimension(:,:,:) :: litrmassrow
    real, pointer, dimension(:,:,:) :: soilcmasrow
    real, pointer, dimension(:,:,:) :: vgbiomas_vegrow
    real, pointer, dimension(:,:,:) :: stemmassrow
    real, pointer, dimension(:,:,:) :: rootmassrow
    real, pointer, dimension(:,:,:) :: litrfallvegrow
    real, pointer, dimension(:,:,:) :: humiftrsvegrow
    real, pointer, dimension(:,:) ::uvaccrow_m
    real, pointer, dimension(:,:) ::vvaccrow_m

    real, pointer, dimension(:) :: laimaxg_mo_g
    real, pointer, dimension(:) :: stemmass_mo_g
    real, pointer, dimension(:) :: rootmass_mo_g
    real, pointer, dimension(:) :: litrmass_mo_g
    real, pointer, dimension(:) :: soilcmas_mo_g
    real, pointer, dimension(:) :: litrfall_mo_g
    real, pointer, dimension(:) :: humiftrs_mo_g
    real, pointer, dimension(:) :: npp_mo_g
    real, pointer, dimension(:) :: gpp_mo_g
    real, pointer, dimension(:) :: nep_mo_g
    real, pointer, dimension(:) :: nbp_mo_g
    real, pointer, dimension(:) :: hetrores_mo_g
    real, pointer, dimension(:) :: autores_mo_g
    real, pointer, dimension(:) :: litres_mo_g
    real, pointer, dimension(:) :: soilcres_mo_g
    real, pointer, dimension(:) :: vgbiomas_mo_g
    real, pointer, dimension(:) :: totcmass_mo_g
    real, pointer, dimension(:) :: emit_co2_mo_g
    real, pointer, dimension(:) :: emit_co_mo_g
    real, pointer, dimension(:) :: emit_ch4_mo_g
    real, pointer, dimension(:) :: emit_nmhc_mo_g
    real, pointer, dimension(:) :: emit_h2_mo_g
    real, pointer, dimension(:) :: emit_nox_mo_g
    real, pointer, dimension(:) :: emit_n2o_mo_g
    real, pointer, dimension(:) :: emit_pm25_mo_g
    real, pointer, dimension(:) :: emit_tpm_mo_g
    real, pointer, dimension(:) :: emit_tc_mo_g
    real, pointer, dimension(:) :: emit_oc_mo_g
    real, pointer, dimension(:) :: emit_bc_mo_g
    real, pointer, dimension(:) :: smfuncveg_mo_g
    real, pointer, dimension(:) :: luc_emc_mo_g
    real, pointer, dimension(:) :: lucltrin_mo_g
    real, pointer, dimension(:) :: lucsocin_mo_g
    real, pointer, dimension(:) :: burnfrac_mo_g
    real, pointer, dimension(:) :: bterm_mo_g
    real, pointer, dimension(:) :: lterm_mo_g
    real, pointer, dimension(:) :: mterm_mo_g
    real, pointer, dimension(:) :: ch4wet1_mo_g
    real, pointer, dimension(:) :: ch4wet2_mo_g
    real, pointer, dimension(:) :: wetfdyn_mo_g
    real, pointer, dimension(:) :: ch4dyn1_mo_g
    real, pointer, dimension(:) :: ch4dyn2_mo_g
    real, pointer, dimension(:) :: ch4soills_mo_g

    ! local
    integer :: i,m,j,nt
    real :: barefrac
    real :: sumfare
    integer :: NDMONTH
    integer :: imonth
    real, dimension(1) :: timeStamp
    real, dimension(icc) :: pftExist

    ! point pointers

    dofire                => c_switch%dofire
    lnduseon              => c_switch%lnduseon
    PFTCompetition        => c_switch%PFTCompetition
    dowetlands            => c_switch%dowetlands
    obswetf               => c_switch%obswetf
    doperpftoutput        => c_switch%doperpftoutput
    dopertileoutput       => c_switch%dopertileoutput
    pftexistrow           => vrot%pftexist
    fcancmxrow            => vrot%fcancmx
    laimaxg_mo            =>ctem_mo%laimaxg_mo
    stemmass_mo           =>ctem_mo%stemmass_mo
    rootmass_mo           =>ctem_mo%rootmass_mo
    litrfallveg_mo           =>ctem_mo%litrfallveg_mo
    humiftrsveg_mo           =>ctem_mo%humiftrsveg_mo
    npp_mo                =>ctem_mo%npp_mo
    gpp_mo                =>ctem_mo%gpp_mo
    vgbiomas_mo           =>ctem_mo%vgbiomas_mo
    autores_mo            =>ctem_mo%autores_mo
    totcmass_mo           =>ctem_mo%totcmass_mo
    litrmass_mo           =>ctem_mo%litrmass_mo
    soilcmas_mo           =>ctem_mo%soilcmas_mo
    nep_mo                =>ctem_mo%nep_mo
    litres_mo             =>ctem_mo%litres_mo
    soilcres_mo           =>ctem_mo%soilcres_mo
    hetrores_mo           =>ctem_mo%hetrores_mo
    nbp_mo                =>ctem_mo%nbp_mo
    emit_co2_mo           =>ctem_mo%emit_co2_mo
    emit_co_mo            =>ctem_mo%emit_co_mo
    emit_ch4_mo           =>ctem_mo%emit_ch4_mo
    emit_nmhc_mo          =>ctem_mo%emit_nmhc_mo
    emit_h2_mo            =>ctem_mo%emit_h2_mo
    emit_nox_mo           =>ctem_mo%emit_nox_mo
    emit_n2o_mo           =>ctem_mo%emit_n2o_mo
    emit_pm25_mo          =>ctem_mo%emit_pm25_mo
    emit_tpm_mo           =>ctem_mo%emit_tpm_mo
    emit_tc_mo            =>ctem_mo%emit_tc_mo
    emit_oc_mo            =>ctem_mo%emit_oc_mo
    emit_bc_mo            =>ctem_mo%emit_bc_mo
    bterm_mo              =>ctem_mo%bterm_mo
    mterm_mo              =>ctem_mo%mterm_mo
    burnfrac_mo           =>ctem_mo%burnfrac_mo
    smfuncveg_mo          =>ctem_mo%smfuncveg_mo

    laimaxg_mo_t          =>ctem_tile_mo%laimaxg_mo_t
    stemmass_mo_t         =>ctem_tile_mo%stemmass_mo_t
    rootmass_mo_t         =>ctem_tile_mo%rootmass_mo_t
    litrfall_mo_t         =>ctem_tile_mo%litrfall_mo_t
    humiftrs_mo_t         =>ctem_tile_mo%humiftrs_mo_t
    npp_mo_t              =>ctem_tile_mo%npp_mo_t
    gpp_mo_t              =>ctem_tile_mo%gpp_mo_t
    vgbiomas_mo_t         =>ctem_tile_mo%vgbiomas_mo_t
    autores_mo_t          =>ctem_tile_mo%autores_mo_t
    totcmass_mo_t         =>ctem_tile_mo%totcmass_mo_t
    litrmass_mo_t         =>ctem_tile_mo%litrmass_mo_t
    soilcmas_mo_t         =>ctem_tile_mo%soilcmas_mo_t
    nep_mo_t              =>ctem_tile_mo%nep_mo_t
    litres_mo_t           =>ctem_tile_mo%litres_mo_t
    soilcres_mo_t         =>ctem_tile_mo%soilcres_mo_t
    hetrores_mo_t         =>ctem_tile_mo%hetrores_mo_t
    nbp_mo_t              =>ctem_tile_mo%nbp_mo_t
    emit_co2_mo_t         =>ctem_tile_mo%emit_co2_mo_t
    emit_co_mo_t          =>ctem_tile_mo%emit_co_mo_t
    emit_ch4_mo_t         =>ctem_tile_mo%emit_ch4_mo_t
    emit_nmhc_mo_t        =>ctem_tile_mo%emit_nmhc_mo_t
    emit_h2_mo_t          =>ctem_tile_mo%emit_h2_mo_t
    emit_nox_mo_t         =>ctem_tile_mo%emit_nox_mo_t
    emit_n2o_mo_t         =>ctem_tile_mo%emit_n2o_mo_t
    emit_pm25_mo_t        =>ctem_tile_mo%emit_pm25_mo_t
    emit_tpm_mo_t         =>ctem_tile_mo%emit_tpm_mo_t
    emit_tc_mo_t          =>ctem_tile_mo%emit_tc_mo_t
    emit_oc_mo_t          =>ctem_tile_mo%emit_oc_mo_t
    emit_bc_mo_t          =>ctem_tile_mo%emit_bc_mo_t
    burnfrac_mo_t         =>ctem_tile_mo%burnfrac_mo_t
    smfuncveg_mo_t        =>ctem_tile_mo%smfuncveg_mo_t
    bterm_mo_t            =>ctem_tile_mo%bterm_mo_t
    luc_emc_mo_t          =>ctem_tile_mo%luc_emc_mo_t
    lterm_mo_t            =>ctem_tile_mo%lterm_mo_t
    lucsocin_mo_t         =>ctem_tile_mo%lucsocin_mo_t
    mterm_mo_t            =>ctem_tile_mo%mterm_mo_t
    lucltrin_mo_t         =>ctem_tile_mo%lucltrin_mo_t
    ch4wet1_mo_t          =>ctem_tile_mo%ch4wet1_mo_t
    ch4wet2_mo_t          =>ctem_tile_mo%ch4wet2_mo_t
    wetfdyn_mo_t          =>ctem_tile_mo%wetfdyn_mo_t
    ch4dyn1_mo_t          =>ctem_tile_mo%ch4dyn1_mo_t
    ch4dyn2_mo_t          =>ctem_tile_mo%ch4dyn2_mo_t
    ch4soills_mo_t        =>ctem_tile_mo%ch4soills_mo_t
    wind_mo_t             =>ctem_tile_mo%wind_mo_t

    gppvegrow         => vrot%gppveg
    nepvegrow         => vrot%nepveg
    nbpvegrow         => vrot%nbpveg
    nppvegrow         => vrot%nppveg
    hetroresvegrow    => vrot%hetroresveg
    autoresvegrow     => vrot%autoresveg
    litresvegrow      => vrot%litresveg
    soilcresvegrow    => vrot%soilcresveg
    rmlvegaccrow      => vrot%rmlvegacc
    rmsvegrow         => vrot%rmsveg
    rmrvegrow         => vrot%rmrveg
    rgvegrow          => vrot%rgveg
    ailcgrow          => vrot%ailcg
    emit_co2row       => vrot%emit_co2
    emit_corow        => vrot%emit_co
    emit_ch4row       => vrot%emit_ch4
    emit_nmhcrow      => vrot%emit_nmhc
    emit_h2row        => vrot%emit_h2
    emit_noxrow       => vrot%emit_nox
    emit_n2orow       => vrot%emit_n2o
    emit_pm25row      => vrot%emit_pm25
    emit_tpmrow       => vrot%emit_tpm
    emit_tcrow        => vrot%emit_tc
    emit_ocrow        => vrot%emit_oc
    emit_bcrow        => vrot%emit_bc
    burnfracrow       => vrot%burnfrac
    burnvegfrow       => vrot%burnvegf
    smfuncvegrow      => vrot%smfuncveg
    btermrow          => vrot%bterm
    ltermrow          => vrot%lterm
    mtermrow          => vrot%mterm
    lucemcomrow       => vrot%lucemcom
    lucltrinrow       => vrot%lucltrin
    lucsocinrow       => vrot%lucsocin
    ch4wet1row        => vrot%ch4wet1
    ch4wet2row        => vrot%ch4wet2
    wetfdynrow        => vrot%wetfdyn
    ch4dyn1row        => vrot%ch4dyn1
    ch4dyn2row        => vrot%ch4dyn2
    ch4soillsrow      => vrot%ch4_soills
    litrmassrow       => vrot%litrmass
    soilcmasrow       => vrot%soilcmas
    vgbiomas_vegrow   => vrot%vgbiomas_veg
    stemmassrow       => vrot%stemmass
    rootmassrow       => vrot%rootmass
    uvaccrow_m        => vrot%uvaccrow_m
    vvaccrow_m        => vrot%vvaccrow_m
    litrfallvegrow    => vrot%litrfallveg
    humiftrsvegrow    => vrot%humiftrsveg

    laimaxg_mo_g        =>ctem_grd_mo%laimaxg_mo_g
    stemmass_mo_g       =>ctem_grd_mo%stemmass_mo_g
    rootmass_mo_g       =>ctem_grd_mo%rootmass_mo_g
    litrmass_mo_g       =>ctem_grd_mo%litrmass_mo_g
    soilcmas_mo_g       =>ctem_grd_mo%soilcmas_mo_g
    litrfall_mo_g       =>ctem_grd_mo%litrfall_mo_g
    humiftrs_mo_g       =>ctem_grd_mo%humiftrs_mo_g
    npp_mo_g            =>ctem_grd_mo%npp_mo_g
    gpp_mo_g            =>ctem_grd_mo%gpp_mo_g
    nep_mo_g            =>ctem_grd_mo%nep_mo_g
    nbp_mo_g            =>ctem_grd_mo%nbp_mo_g
    hetrores_mo_g       =>ctem_grd_mo%hetrores_mo_g
    autores_mo_g        =>ctem_grd_mo%autores_mo_g
    litres_mo_g         =>ctem_grd_mo%litres_mo_g
    soilcres_mo_g       =>ctem_grd_mo%soilcres_mo_g
    vgbiomas_mo_g       =>ctem_grd_mo%vgbiomas_mo_g
    totcmass_mo_g       =>ctem_grd_mo%totcmass_mo_g
    emit_co2_mo_g       =>ctem_grd_mo%emit_co2_mo_g
    emit_co_mo_g        =>ctem_grd_mo%emit_co_mo_g
    emit_ch4_mo_g       =>ctem_grd_mo%emit_ch4_mo_g
    emit_nmhc_mo_g      =>ctem_grd_mo%emit_nmhc_mo_g
    emit_h2_mo_g        =>ctem_grd_mo%emit_h2_mo_g
    emit_nox_mo_g       =>ctem_grd_mo%emit_nox_mo_g
    emit_n2o_mo_g       =>ctem_grd_mo%emit_n2o_mo_g
    emit_pm25_mo_g      =>ctem_grd_mo%emit_pm25_mo_g
    emit_tpm_mo_g       =>ctem_grd_mo%emit_tpm_mo_g
    emit_tc_mo_g        =>ctem_grd_mo%emit_tc_mo_g
    emit_oc_mo_g        =>ctem_grd_mo%emit_oc_mo_g
    emit_bc_mo_g        =>ctem_grd_mo%emit_bc_mo_g
    smfuncveg_mo_g      =>ctem_grd_mo%smfuncveg_mo_g
    luc_emc_mo_g        =>ctem_grd_mo%luc_emc_mo_g
    lucltrin_mo_g       =>ctem_grd_mo%lucltrin_mo_g
    lucsocin_mo_g       =>ctem_grd_mo%lucsocin_mo_g
    burnfrac_mo_g       =>ctem_grd_mo%burnfrac_mo_g
    bterm_mo_g          =>ctem_grd_mo%bterm_mo_g
    lterm_mo_g          =>ctem_grd_mo%lterm_mo_g
    mterm_mo_g          =>ctem_grd_mo%mterm_mo_g
    ch4wet1_mo_g        =>ctem_grd_mo%ch4wet1_mo_g
    ch4wet2_mo_g        =>ctem_grd_mo%ch4wet2_mo_g
    wetfdyn_mo_g        =>ctem_grd_mo%wetfdyn_mo_g
    ch4dyn1_mo_g        =>ctem_grd_mo%ch4dyn1_mo_g
    ch4dyn2_mo_g        =>ctem_grd_mo%ch4dyn2_mo_g
    ch4soills_mo_g      =>ctem_grd_mo%ch4soills_mo_g

    !> ------------

    !> Accumulate monthly outputs

    i = 1 ! offline nlat is always 1 so this array position is always 1.
    do 863 m = 1,nmtest
        do j=1,icc

            !> Accumulate monthly outputs at the per PFT level.
            if (ailcgrow(i,m,j) .gt. laimaxg_mo(i,m,j)) then
                laimaxg_mo(i,m,j)=ailcgrow(i,m,j)
            end if

            npp_mo(i,m,j)=npp_mo(i,m,j)+nppvegrow(i,m,j)
            gpp_mo(i,m,j)=gpp_mo(i,m,j)+gppvegrow(i,m,j)
            nep_mo(i,m,j)=nep_mo(i,m,j)+nepvegrow(i,m,j)
            nbp_mo(i,m,j)=nbp_mo(i,m,j)+nbpvegrow(i,m,j)
            hetrores_mo(i,m,j)=hetrores_mo(i,m,j)+hetroresvegrow(i,m,j)
            autores_mo(i,m,j) =autores_mo(i,m,j)+autoresvegrow(i,m,j)
            litres_mo(i,m,j)  =litres_mo(i,m,j) +litresvegrow(i,m,j)
            soilcres_mo(i,m,j) =soilcres_mo(i,m,j) +soilcresvegrow(i,m,j)
            emit_co2_mo(i,m,j)=emit_co2_mo(i,m,j)+emit_co2row(i,m,j)
            emit_co_mo(i,m,j) =emit_co_mo(i,m,j)+emit_corow(i,m,j)
            emit_ch4_mo(i,m,j) =emit_ch4_mo(i,m,j)+emit_ch4row(i,m,j)
            emit_nmhc_mo(i,m,j)=emit_nmhc_mo(i,m,j)+emit_nmhcrow(i,m,j)
            emit_h2_mo(i,m,j) =emit_h2_mo(i,m,j)+emit_h2row(i,m,j)
            emit_nox_mo(i,m,j) =emit_nox_mo(i,m,j)+emit_noxrow(i,m,j)
            emit_n2o_mo(i,m,j) =emit_n2o_mo(i,m,j)+emit_n2orow(i,m,j)
            emit_pm25_mo(i,m,j)=emit_pm25_mo(i,m,j)+emit_pm25row(i,m,j)
            emit_tpm_mo(i,m,j) =emit_tpm_mo(i,m,j)+emit_tpmrow(i,m,j)
            emit_tc_mo(i,m,j) =emit_tc_mo(i,m,j)+emit_tcrow(i,m,j)
            emit_oc_mo(i,m,j) =emit_oc_mo(i,m,j)+emit_ocrow(i,m,j)
            emit_bc_mo(i,m,j) =emit_bc_mo(i,m,j)+emit_bcrow(i,m,j)
            bterm_mo(i,m,j) = bterm_mo(i,m,j) + btermrow(i,m,j)
            mterm_mo(i,m,j) = mterm_mo(i,m,j) + mtermrow(i,m,j)
            burnfrac_mo(i,m,j) =burnfrac_mo(i,m,j)+burnvegfrow(i,m,j)
            smfuncveg_mo(i,m,j) =smfuncveg_mo(i,m,j) + smfuncvegrow(i,m,j)
            litrfallveg_mo(i,m,j) = litrfallveg_mo(i,m,j) + litrfallvegrow(i,m,j)
            humiftrsveg_mo(i,m,j) = humiftrsveg_mo(i,m,j) + humiftrsvegrow(i,m,j)

        end do !j

        !> Also do the bare ground
        nep_mo(i,m,iccp1)=nep_mo(i,m,iccp1)+nepvegrow(i,m,iccp1)
        nbp_mo(i,m,iccp1)=nbp_mo(i,m,iccp1)+nbpvegrow(i,m,iccp1)
        hetrores_mo(i,m,iccp1)=hetrores_mo(i,m,iccp1)+hetroresvegrow(i,m,iccp1)
        litres_mo(i,m,iccp1)  =litres_mo(i,m,iccp1)+litresvegrow(i,m,iccp1)
        soilcres_mo(i,m,iccp1) =soilcres_mo(i,m,iccp1) +soilcresvegrow(i,m,iccp1)

        !> Accumulate monthly outputs at the per tile level.
        luc_emc_mo_t(i,m) =luc_emc_mo_t(i,m)+lucemcomrow(i,m)
        lucsocin_mo_t(i,m) =lucsocin_mo_t(i,m)+lucsocinrow(i,m)
        lucltrin_mo_t(i,m) =lucltrin_mo_t(i,m)+lucltrinrow(i,m)
        ch4wet1_mo_t(i,m) = ch4wet1_mo_t(i,m) + ch4wet1row(i,m)
        ch4wet2_mo_t(i,m) = ch4wet2_mo_t(i,m) + ch4wet2row(i,m)
        wetfdyn_mo_t(i,m) = wetfdyn_mo_t(i,m) + wetfdynrow(i,m)
        ch4dyn1_mo_t(i,m) = ch4dyn1_mo_t(i,m) + ch4dyn1row(i,m)
        ch4dyn2_mo_t(i,m) = ch4dyn2_mo_t(i,m) + ch4dyn2row(i,m)
        ch4soills_mo_t(i,m) = ch4soills_mo_t(i,m) + ch4soillsrow(i,m)
        lterm_mo_t(i,m) = lterm_mo_t(i,m) + ltermrow(i,m)
        wind_mo_t(i,m) = wind_mo_t(i,m) + (sqrt(uvaccrow_m(i,m)**2.0 + vvaccrow_m(i,m)**2.0))*3.6 !>take mean wind speed and convert to km/h

863     continue ! m

    do 865 nt=1,nmon

        if(iday.eq.mmday(nt))then

            !> Do the mid-month variables (these are not accumulated, we just keep the mid month value for printing in the monthly file)

            do 866 m=1,nmtest

                do 867 j=1,icc

                    vgbiomas_mo(i,m,j)=vgbiomas_vegrow(i,m,j)
                    litrmass_mo(i,m,j)=litrmassrow(i,m,j)
                    soilcmas_mo(i,m,j)=soilcmasrow(i,m,j)
                    stemmass_mo(i,m,j)=stemmassrow(i,m,j)
                    rootmass_mo(i,m,j)=rootmassrow(i,m,j)
                    totcmass_mo(i,m,j)=vgbiomas_vegrow(i,m,j) + litrmassrow(i,m,j) + soilcmasrow(i,m,j)

867                 continue

                !> Do the bare fraction too
                litrmass_mo(i,m,iccp1)=litrmassrow(i,m,iccp1)
                soilcmas_mo(i,m,iccp1)=soilcmasrow(i,m,iccp1)
                totcmass_mo(i,m,iccp1)=soilcmasrow(i,m,iccp1) + litrmassrow(i,m,iccp1)

                barefrac=1.0

                !> Now find the per tile values:
                do j=1,icc
                    vgbiomas_mo_t(i,m)=vgbiomas_mo_t(i,m)+vgbiomas_mo(i,m,j)*fcancmxrow(i,m,j)
                    litrmass_mo_t(i,m)=litrmass_mo_t(i,m)+litrmass_mo(i,m,j)*fcancmxrow(i,m,j)
                    soilcmas_mo_t(i,m)=soilcmas_mo_t(i,m)+soilcmas_mo(i,m,j)*fcancmxrow(i,m,j)
                    stemmass_mo_t(i,m)=stemmass_mo_t(i,m)+stemmass_mo(i,m,j)*fcancmxrow(i,m,j)
                    rootmass_mo_t(i,m)=rootmass_mo_t(i,m)+rootmass_mo(i,m,j)*fcancmxrow(i,m,j)
                    totcmass_mo_t(i,m)=totcmass_mo_t(i,m)+totcmass_mo(i,m,j)*fcancmxrow(i,m,j)
                    barefrac=barefrac-fcancmxrow(i,m,j)
                end do

                !>Also add in the bare fraction contributions.
                litrmass_mo_t(i,m)=litrmass_mo_t(i,m)+litrmass_mo(i,m,iccp1)*barefrac
                soilcmas_mo_t(i,m)=soilcmas_mo_t(i,m)+soilcmas_mo(i,m,iccp1)*barefrac
                totcmass_mo_t(i,m)=totcmass_mo_t(i,m)+(litrmass_mo(i,m,iccp1)+soilcmas_mo(i,m,iccp1))*barefrac

                !> Now find the gridcell level values:
                vgbiomas_mo_g(i)=vgbiomas_mo_g(i)+vgbiomas_mo_t(i,m)*FAREROT(i,m)
                litrmass_mo_g(i)=litrmass_mo_g(i)+litrmass_mo_t(i,m)*FAREROT(i,m)
                soilcmas_mo_g(i)=soilcmas_mo_g(i)+soilcmas_mo_t(i,m)*FAREROT(i,m)
                stemmass_mo_g(i)=stemmass_mo_g(i)+stemmass_mo_t(i,m)*FAREROT(i,m)
                rootmass_mo_g(i)=rootmass_mo_g(i)+rootmass_mo_t(i,m)*FAREROT(i,m)
                totcmass_mo_g(i)=totcmass_mo_g(i)+totcmass_mo_t(i,m)*FAREROT(i,m)

866             continue  !nmtest loop.

        endif ! mmday (mid-month instantaneous value)

        if(iday.eq.monthend(nt+1))then

            !> Do the end of month variables
            ndmonth=(monthend(nt+1)-monthend(nt))*nday

            do 900 m = 1,nmtest


                !> Convert some quantities into per day values
                wetfdyn_mo_t(i,m)=wetfdyn_mo_t(i,m)*(1./real(monthdays(nt)))
                lterm_mo_t(i,m)=lterm_mo_t(i,m)*(1./real(monthdays(nt)))
                wind_mo_t(i,m) = wind_mo_t(i,m)*(1./real(monthdays(nt)))
                do j = 1, icc
                    bterm_mo(i,m,j)=bterm_mo(i,m,j)*(1./real(monthdays(nt)))
                    mterm_mo(i,m,j)=mterm_mo(i,m,j)*(1./real(monthdays(nt)))
                    smfuncveg_mo(i,m,j) =smfuncveg_mo(i,m,j) *(1./real(monthdays(nt)))
                end do

                barefrac=1.0

                do j=1,icc

                    !> Find the monthly outputs at the per tile level from the outputs at the per PFT level
                    npp_mo_t(i,m)=npp_mo_t(i,m)+npp_mo(i,m,j)*fcancmxrow(i,m,j)
                    gpp_mo_t(i,m)=gpp_mo_t(i,m)+gpp_mo(i,m,j)*fcancmxrow(i,m,j)
                    nep_mo_t(i,m)=nep_mo_t(i,m)+nep_mo(i,m,j)*fcancmxrow(i,m,j)
                    nbp_mo_t(i,m)=nbp_mo_t(i,m)+nbp_mo(i,m,j)*fcancmxrow(i,m,j)
                    hetrores_mo_t(i,m)=hetrores_mo_t(i,m)+hetrores_mo(i,m,j)*fcancmxrow(i,m,j)
                    autores_mo_t(i,m) =autores_mo_t(i,m)+autores_mo(i,m,j)*fcancmxrow(i,m,j)
                    litres_mo_t(i,m)  =litres_mo_t(i,m) +litres_mo(i,m,j)*fcancmxrow(i,m,j)
                    soilcres_mo_t(i,m) =soilcres_mo_t(i,m) +soilcres_mo(i,m,j)*fcancmxrow(i,m,j)
                    emit_co2_mo_t(i,m)=emit_co2_mo_t(i,m)+emit_co2_mo(i,m,j)*fcancmxrow(i,m,j)
                    emit_co_mo_t(i,m) =emit_co_mo_t(i,m)+emit_co_mo(i,m,j)*fcancmxrow(i,m,j)
                    emit_ch4_mo_t(i,m) =emit_ch4_mo_t(i,m)+emit_ch4_mo(i,m,j)*fcancmxrow(i,m,j)
                    emit_nmhc_mo_t(i,m)=emit_nmhc_mo_t(i,m)+emit_nmhc_mo(i,m,j)*fcancmxrow(i,m,j)
                    emit_h2_mo_t(i,m) =emit_h2_mo_t(i,m)+emit_h2_mo(i,m,j)*fcancmxrow(i,m,j)
                    emit_nox_mo_t(i,m) =emit_nox_mo_t(i,m)+emit_nox_mo(i,m,j)*fcancmxrow(i,m,j)
                    emit_n2o_mo_t(i,m) =emit_n2o_mo_t(i,m)+emit_n2o_mo(i,m,j)*fcancmxrow(i,m,j)
                    emit_pm25_mo_t(i,m)=emit_pm25_mo_t(i,m)+emit_pm25_mo(i,m,j)*fcancmxrow(i,m,j)
                    emit_tpm_mo_t(i,m) =emit_tpm_mo_t(i,m)+emit_tpm_mo(i,m,j)*fcancmxrow(i,m,j)
                    emit_tc_mo_t(i,m) =emit_tc_mo_t(i,m)+emit_tc_mo(i,m,j)*fcancmxrow(i,m,j)
                    emit_oc_mo_t(i,m) =emit_oc_mo_t(i,m)+emit_oc_mo(i,m,j)*fcancmxrow(i,m,j)
                    emit_bc_mo_t(i,m) =emit_bc_mo_t(i,m)+emit_bc_mo(i,m,j)*fcancmxrow(i,m,j)
                    bterm_mo_t(i,m) =bterm_mo_t(i,m)+bterm_mo(i,m,j)*fcancmxrow(i,m,j)
                    mterm_mo_t(i,m) =mterm_mo_t(i,m)+mterm_mo(i,m,j)*fcancmxrow(i,m,j)
                    smfuncveg_mo_t(i,m) =smfuncveg_mo_t(i,m)+smfuncveg_mo(i,m,j)*fcancmxrow(i,m,j)
                    burnfrac_mo_t(i,m) =burnfrac_mo_t(i,m)+burnfrac_mo(i,m,j)*fcancmxrow(i,m,j)
                    laimaxg_mo_t(i,m)=laimaxg_mo_t(i,m)+laimaxg_mo(i,m,j)*fcancmxrow(i,m,j)
                    litrfall_mo_t(i,m)=litrfall_mo_t(i,m)+litrfallveg_mo(i,m,j)*fcancmxrow(i,m,j)
                    humiftrs_mo_t(i,m)=humiftrs_mo_t(i,m)+humiftrsveg_mo(i,m,j)*fcancmxrow(i,m,j)
                    barefrac=barefrac-fcancmxrow(i,m,j)

                end do !j

                nep_mo_t(i,m)=nep_mo_t(i,m)+nep_mo(i,m,iccp1)*barefrac
                nbp_mo_t(i,m)=nbp_mo_t(i,m)+nbp_mo(i,m,iccp1)*barefrac
                hetrores_mo_t(i,m)=hetrores_mo_t(i,m)+hetrores_mo(i,m,iccp1)*barefrac
                litres_mo_t(i,m)  =litres_mo_t(i,m) +litres_mo(i,m,iccp1)*barefrac
                soilcres_mo_t(i,m)=soilcres_mo_t(i,m)+soilcres_mo(i,m,iccp1)*barefrac
                humiftrs_mo_t(i,m)=humiftrs_mo_t(i,m)+humiftrsveg_mo(i,m,iccp1)*barefrac

                !> Find the monthly outputs at the per grid cell level from the outputs at the per tile level
                npp_mo_g(i)=npp_mo_g(i)+npp_mo_t(i,m)*FAREROT(i,m)
                gpp_mo_g(i)=gpp_mo_g(i)+gpp_mo_t(i,m)*FAREROT(i,m)
                nep_mo_g(i)=nep_mo_g(i)+nep_mo_t(i,m)*FAREROT(i,m)
                nbp_mo_g(i)=nbp_mo_g(i)+nbp_mo_t(i,m)*FAREROT(i,m)
                hetrores_mo_g(i)=hetrores_mo_g(i)+hetrores_mo_t(i,m)*FAREROT(i,m)
                autores_mo_g(i) =autores_mo_g(i) +autores_mo_t(i,m)*FAREROT(i,m)
                litres_mo_g(i)  =litres_mo_g(i) +litres_mo_t(i,m)*FAREROT(i,m)
                soilcres_mo_g(i) =soilcres_mo_g(i)+ soilcres_mo_t(i,m)*FAREROT(i,m)
                laimaxg_mo_g(i)=laimaxg_mo_g(i)+laimaxg_mo_t(i,m)*FAREROT(i,m)
                emit_co2_mo_g(i)=emit_co2_mo_g(i)+emit_co2_mo_t(i,m)*FAREROT(i,m)
                emit_co_mo_g(i) =emit_co_mo_g(i)+emit_co_mo_t(i,m)*FAREROT(i,m)
                emit_ch4_mo_g(i) =emit_ch4_mo_g(i)+emit_ch4_mo_t(i,m)*FAREROT(i,m)
                emit_nmhc_mo_g(i)=emit_nmhc_mo_g(i)+emit_nmhc_mo_t(i,m)*FAREROT(i,m)
                emit_h2_mo_g(i) =emit_h2_mo_g(i)+emit_h2_mo_t(i,m)*FAREROT(i,m)
                emit_nox_mo_g(i) =emit_nox_mo_g(i)+emit_nox_mo_t(i,m)*FAREROT(i,m)
                emit_n2o_mo_g(i) =emit_n2o_mo_g(i)+emit_n2o_mo_t(i,m)*FAREROT(i,m)
                emit_pm25_mo_g(i) =emit_pm25_mo_g(i)+emit_pm25_mo_t(i,m)*FAREROT(i,m)
                emit_tpm_mo_g(i) =emit_tpm_mo_g(i)+emit_tpm_mo_t(i,m)*FAREROT(i,m)
                emit_tc_mo_g(i) =emit_tc_mo_g(i)+emit_tc_mo_t(i,m)*FAREROT(i,m)
                emit_oc_mo_g(i) =emit_oc_mo_g(i)+emit_oc_mo_t(i,m)*FAREROT(i,m)
                emit_bc_mo_g(i) =emit_bc_mo_g(i)+emit_bc_mo_t(i,m)*FAREROT(i,m)
                burnfrac_mo_g(i)=burnfrac_mo_g(i)+burnfrac_mo_t(i,m)*FAREROT(i,m)
                luc_emc_mo_g(i) =luc_emc_mo_g(i)+luc_emc_mo_t(i,m)*FAREROT(i,m)
                lucsocin_mo_g(i) =lucsocin_mo_g(i)+lucsocin_mo_t(i,m)*FAREROT(i,m)
                lucltrin_mo_g(i) =lucltrin_mo_g(i)+lucltrin_mo_t(i,m)*FAREROT(i,m)
                ch4wet1_mo_g(i) = ch4wet1_mo_g(i) +ch4wet1_mo_t(i,m)*FAREROT(i,m)
                ch4wet2_mo_g(i) = ch4wet2_mo_g(i)+ch4wet2_mo_t(i,m)*FAREROT(i,m)
                wetfdyn_mo_g(i) = wetfdyn_mo_g(i)+wetfdyn_mo_t(i,m)*FAREROT(i,m)
                ch4dyn1_mo_g(i) = ch4dyn1_mo_g(i)+ch4dyn1_mo_t(i,m)*FAREROT(i,m)
                ch4dyn2_mo_g(i) = ch4dyn2_mo_g(i)+ch4dyn2_mo_t(i,m)*FAREROT(i,m)
                ch4soills_mo_g(i) = ch4soills_mo_g(i)+ch4soills_mo_t(i,m)*FAREROT(i,m)
                smfuncveg_mo_g(i)=smfuncveg_mo_g(i)+smfuncveg_mo_t(i,m)*FAREROT(i,m)
                bterm_mo_g(i) =bterm_mo_g(i)+bterm_mo_t(i,m)*FAREROT(i,m)
                lterm_mo_g(i) =lterm_mo_g(i)+lterm_mo_t(i,m)*FAREROT(i,m)
                mterm_mo_g(i) =mterm_mo_g(i)+mterm_mo_t(i,m)*FAREROT(i,m)
                litrfall_mo_g(i)=litrfall_mo_g(i)+litrfall_mo_t(i,m)*FAREROT(i,m)
                humiftrs_mo_g(i)=humiftrs_mo_g(i)+humiftrs_mo_t(i,m)*FAREROT(i,m)

900             continue

            imonth=nt

            ! Prepare the timestamp for this month !FLAG this isn't correct for leap yet. Need to look at yrs before iyear too!
            !Take one day off so it is the last day of the month rather than the first day of the next month.
            timeStamp(1) = (iyear - 1 - refyr) * lastDOY + monthend(imonth+1) - 1

            call writeOutput1D(lonLocalIndex,latLocalIndex,'laimaxg_mo_g' ,timeStamp,'lai', [laimaxg_mo_g(i)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'vgbiomas_mo_g',timeStamp,'cVeg',[vgbiomas_mo_g(i)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'litrmass_mo_g',timeStamp,'cLitter',[litrmass_mo_g(i)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcmas_mo_g',timeStamp,'cSoil',[soilcmas_mo_g(i)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'npp_mo_g'     ,timeStamp,'npp',[npp_mo_g(i)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'gpp_mo_g'     ,timeStamp,'gpp',[gpp_mo_g(i)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nep_mo_g'     ,timeStamp,'nep',[nep_mo_g(i)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nbp_mo_g'     ,timeStamp,'nbp',[nbp_mo_g(i)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'hetrores_mo_g',timeStamp,'rh',[hetrores_mo_g(i)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'autores_mo_g' ,timeStamp,'ra',[autores_mo_g(i)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'litres_mo_g'  ,timeStamp,'rhLitter',[litres_mo_g(i)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcres_mo_g',timeStamp,'rhSoil',[soilcres_mo_g(i)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'litrfall_mo_g' ,timeStamp,'fVegLitter',[litrfall_mo_g(i)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'humiftrs_mo_g' ,timeStamp,'fLitterSoil',[humiftrs_mo_g(i)])
            do m=1,nmtest
                sumfare = 0.0
                do j=1,icc
                    sumfare=sumfare+fcancmxrow(i,m,j)
                end do !j
            end do !m
            call writeOutput1D(lonLocalIndex,latLocalIndex,'fcancmxrow_mo_g' ,timeStamp,'landCoverFrac',[fcancmxrow(i,1,1:icc),1.-sumfare]) !flag only set up for one tile!

            if (dofire) then
                call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_co2_mo_g' ,timeStamp,'fFire',[emit_co2_mo_g(i)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'burnfrac_mo_g' ,timeStamp,'burntFractionAll',[burnfrac_mo_g(i)*100.])
            end if
            if (lnduseon) then
                call writeOutput1D(lonLocalIndex,latLocalIndex,'luc_emc_mo_g' ,timeStamp,'fDeforestToAtmos',[luc_emc_mo_g(i)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'lucltrin_mo_g' ,timeStamp,'fDeforestToLitter',[lucltrin_mo_g(i)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'lucsocin_mo_g' ,timeStamp,'fDeforestToSoil',[lucsocin_mo_g(i)])
            end if

            if (PFTCompetition) then
                pftExist = 0.0
                do j=1,icc
                    if (pftexistrow(i,1,j)) then
                        pftExist(j) = 1.0
                    end if
                end do
                call writeOutput1D(lonLocalIndex,latLocalIndex,'pftexistrow_yr_g' ,timeStamp,'landCoverExist',[pftExist]) !flag only set up for one tile!
            end if
            if (dowetlands .or. obswetf) then
                call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4wet1_mo_g' ,timeStamp,'wetlandCH4spec',[ch4wet1_mo_g(i)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4dyn1_mo_g' ,timeStamp,'wetlandCH4dyn',[ch4dyn1_mo_g(i)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'wetfdyn_mo_g' ,timeStamp,'wetlandFrac',[wetfdyn_mo_g(i)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4soills_mo_g' ,timeStamp,'soilCH4cons',[ch4soills_mo_g(i)])
            end if

            if (doperpftoutput) then
                if (nmtest > 1) then
                    print*,'Per PFT and per tile outputs together not implemented yet'
                else
                    m = 1
                    call writeOutput1D(lonLocalIndex,latLocalIndex,'laimaxg_mo' ,timeStamp,'lai', [laimaxg_mo(i,m,:)])
                    call writeOutput1D(lonLocalIndex,latLocalIndex,'vgbiomas_mo',timeStamp,'cVeg',[vgbiomas_mo(i,m,:)])
                    call writeOutput1D(lonLocalIndex,latLocalIndex,'litrmass_mo',timeStamp,'cLitter',[litrmass_mo(i,m,:)])
                    call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcmas_mo',timeStamp,'cSoil',[soilcmas_mo(i,m,:)])
                    call writeOutput1D(lonLocalIndex,latLocalIndex,'npp_mo'     ,timeStamp,'npp',[npp_mo(i,m,:)])
                    call writeOutput1D(lonLocalIndex,latLocalIndex,'gpp_mo'     ,timeStamp,'gpp',[gpp_mo(i,m,:)])
                    call writeOutput1D(lonLocalIndex,latLocalIndex,'nep_mo'     ,timeStamp,'nep',[nep_mo(i,m,:)])
                    call writeOutput1D(lonLocalIndex,latLocalIndex,'nbp_mo'     ,timeStamp,'nbp',[nbp_mo(i,m,:)])
                    call writeOutput1D(lonLocalIndex,latLocalIndex,'hetrores_mo',timeStamp,'rh',[hetrores_mo(i,m,:)])
                    call writeOutput1D(lonLocalIndex,latLocalIndex,'autores_mo' ,timeStamp,'ra',[autores_mo(i,m,:)])
                    call writeOutput1D(lonLocalIndex,latLocalIndex,'litres_mo'  ,timeStamp,'rhLitter',[litres_mo(i,m,:)])
                    call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcres_mo',timeStamp,'rhSoil',[soilcres_mo(i,m,:)])
                    call writeOutput1D(lonLocalIndex,latLocalIndex,'litrfallveg_mo' ,timeStamp,'fVegLitter',[litrfallveg_mo(i,m,:)])
                    call writeOutput1D(lonLocalIndex,latLocalIndex,'humiftrsveg_mo' ,timeStamp,'fLitterSoil',[humiftrsveg_mo(i,m,:)])
                    if (dofire .or. lnduseon) then
                        call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_co2_mo' ,timeStamp,'fFire',[emit_co2_mo(i,m,:)])
                        call writeOutput1D(lonLocalIndex,latLocalIndex,'burnfrac_mo' ,timeStamp,'burntFractionAll',[burnfrac_mo(i,m,:)*100.])
!                            smfuncveg_mo(i,m,j), &
!                            bterm_mo(i,m,j),lterm_mo_t(i,m),mterm_mo(i,m,j),wind_mo_t &
                    end if
                end if
            end if

            if (dopertileoutput) then

                if (nmtest > 1) then
                    call writeOutput1D(lonLocalIndex,latLocalIndex,'laimaxg_mo_t' ,timeStamp,'lai', [laimaxg_mo_t(i,:)])
                    call writeOutput1D(lonLocalIndex,latLocalIndex,'vgbiomas_mo_t',timeStamp,'cVeg',[vgbiomas_mo_t(i,:)])
                    call writeOutput1D(lonLocalIndex,latLocalIndex,'litrmass_mo_t',timeStamp,'cLitter',[litrmass_mo_t(i,:)])
                    call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcmas_mo_t',timeStamp,'cSoil',[soilcmas_mo_t(i,:)])
                    call writeOutput1D(lonLocalIndex,latLocalIndex,'npp_mo_t'     ,timeStamp,'npp',[npp_mo_t(i,:)])
                    call writeOutput1D(lonLocalIndex,latLocalIndex,'gpp_mo_t'     ,timeStamp,'gpp',[gpp_mo_t(i,:)])
                    call writeOutput1D(lonLocalIndex,latLocalIndex,'nep_mo_t'     ,timeStamp,'nep',[nep_mo_t(i,:)])
                    call writeOutput1D(lonLocalIndex,latLocalIndex,'nbp_mo_t'     ,timeStamp,'nbp',[nbp_mo_t(i,:)])
                    call writeOutput1D(lonLocalIndex,latLocalIndex,'hetrores_mo_t',timeStamp,'rh',[hetrores_mo_t(i,:)])
                    call writeOutput1D(lonLocalIndex,latLocalIndex,'autores_mo_t' ,timeStamp,'ra',[autores_mo_t(i,:)])
                    call writeOutput1D(lonLocalIndex,latLocalIndex,'litres_mo_t'  ,timeStamp,'rhLitter',[litres_mo_t(i,:)])
                    call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcres_mo_t',timeStamp,'rhSoil',[soilcres_mo_t(i,:)])
                    call writeOutput1D(lonLocalIndex,latLocalIndex,'litrfall_mo_t' ,timeStamp,'fVegLitter',[litrfall_mo_t(i,:)])
                    call writeOutput1D(lonLocalIndex,latLocalIndex,'humiftrs_mo_t' ,timeStamp,'fLitterSoil',[humiftrs_mo_t(i,:)])
                    if (dofire .or. lnduseon) then
                        call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_co2_mo_t' ,timeStamp,'fFire',[emit_co2_mo_t(i,:)])
                        call writeOutput1D(lonLocalIndex,latLocalIndex,'burnfrac_mo_t' ,timeStamp,'burntFractionAll',[burnfrac_mo_t(i,:)*100.])
                        call writeOutput1D(lonLocalIndex,latLocalIndex,'luc_emc_mo_t' ,timeStamp,'fDeforestToAtmos',[luc_emc_mo_t(i,:)])
                        call writeOutput1D(lonLocalIndex,latLocalIndex,'lucltrin_mo_t' ,timeStamp,'fDeforestToLitter',[lucltrin_mo_t(i,:)])
                        call writeOutput1D(lonLocalIndex,latLocalIndex,'lucsocin_mo_t' ,timeStamp,'fDeforestToSoil',[lucsocin_mo_t(i,:)])
                    end if
                    if (dowetlands .or. obswetf) then
                        call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4wet1_mo_t' ,timeStamp,'wetlandCH4spec',[ch4wet1_mo_t(i,:)])
                        call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4dyn1_mo_t' ,timeStamp,'wetlandCH4dyn',[ch4dyn1_mo_t(i,:)])
                        call writeOutput1D(lonLocalIndex,latLocalIndex,'wetfdyn_mo_t' ,timeStamp,'wetlandFrac',[wetfdyn_mo_t(i,:)])
                        call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4soills_mo_t' ,timeStamp,'soilCH4cons',[ch4soills_mo_t(i,:)])
                    end if
                end if
            end if

            !> Reset all end of month accumulated arrays
            call resetmonthend(nltest,nmtest)

        end if ! end of month

865     continue ! nmon

end subroutine ctem_monthly_aw
!>@}
!==============================================================================================================

!>\ingroup io_driver_ctem_annual_aw
!>@{
!> Accumulate and write out the annual biogeochemical (CTEM) outputs

subroutine ctem_annual_aw(lonLocalIndex,latLocalIndex,iday,imonth,iyear,nltest,nmtest,FAREROT,lastDOY)

    use ctem_statevars,     only : ctem_tile_yr, vrot, ctem_grd_yr, c_switch, ctem_yr, &
    resetyearend
    use ctem_params, only : icc,iccp1,seed,nmos
    !use fileIOModule
    use outputManager, only : writeOutput1D,refyr

    implicit none

    ! arguments
    integer, intent(in) :: lonLocalIndex,latLocalIndex
    integer, intent(in) :: nltest
    integer, intent(in) :: nmtest
    integer, intent(in) :: iday
    integer, intent(in) :: imonth
    real, intent(in), dimension(:,:) :: FAREROT
    integer, intent(in) :: iyear
    integer, intent(in) :: lastDOY

    ! pointers

    logical, pointer :: dofire
    logical, pointer :: lnduseon
    logical, pointer :: PFTCompetition
    logical, pointer :: dowetlands
    logical, pointer :: obswetf
    logical, pointer :: doperpftoutput
    logical, pointer :: dopertileoutput

    real, pointer, dimension(:,:,:) :: laimaxg_yr
    real, pointer, dimension(:,:,:) :: stemmass_yr
    real, pointer, dimension(:,:,:) :: rootmass_yr
    real, pointer, dimension(:,:,:) :: npp_yr
    real, pointer, dimension(:,:,:) :: gpp_yr
    real, pointer, dimension(:,:,:) :: vgbiomas_yr
    real, pointer, dimension(:,:,:) :: autores_yr
    real, pointer, dimension(:,:,:) :: totcmass_yr
    real, pointer, dimension(:,:,:) :: litrmass_yr
    real, pointer, dimension(:,:,:) :: soilcmas_yr
    real, pointer, dimension(:,:,:) :: nep_yr
    real, pointer, dimension(:,:,:) :: litres_yr
    real, pointer, dimension(:,:,:) :: soilcres_yr
    real, pointer, dimension(:,:,:) :: hetrores_yr
    real, pointer, dimension(:,:,:) :: nbp_yr
    real, pointer, dimension(:,:,:) :: emit_co2_yr
    real, pointer, dimension(:,:,:) :: emit_co_yr
    real, pointer, dimension(:,:,:) :: emit_ch4_yr
    real, pointer, dimension(:,:,:) :: emit_nmhc_yr
    real, pointer, dimension(:,:,:) :: emit_h2_yr
    real, pointer, dimension(:,:,:) :: emit_nox_yr
    real, pointer, dimension(:,:,:) :: emit_n2o_yr
    real, pointer, dimension(:,:,:) :: emit_pm25_yr
    real, pointer, dimension(:,:,:) :: emit_tpm_yr
    real, pointer, dimension(:,:,:) :: emit_tc_yr
    real, pointer, dimension(:,:,:) :: emit_oc_yr
    real, pointer, dimension(:,:,:) :: emit_bc_yr
    real, pointer, dimension(:,:,:) :: bterm_yr
    real, pointer, dimension(:,:,:) :: mterm_yr
    real, pointer, dimension(:,:,:) :: smfuncveg_yr
    real, pointer, dimension(:,:,:) :: burnfrac_yr
    real, pointer, dimension(:,:,:) :: veghght_yr

    real, pointer, dimension(:,:) :: laimaxg_yr_t
    real, pointer, dimension(:,:) :: stemmass_yr_t
    real, pointer, dimension(:,:) :: rootmass_yr_t
    real, pointer, dimension(:,:) :: npp_yr_t
    real, pointer, dimension(:,:) :: gpp_yr_t
    real, pointer, dimension(:,:) :: vgbiomas_yr_t
    real, pointer, dimension(:,:) :: autores_yr_t
    real, pointer, dimension(:,:) :: totcmass_yr_t
    real, pointer, dimension(:,:) :: litrmass_yr_t
    real, pointer, dimension(:,:) :: soilcmas_yr_t
    real, pointer, dimension(:,:) :: nep_yr_t
    real, pointer, dimension(:,:) :: litres_yr_t
    real, pointer, dimension(:,:) :: soilcres_yr_t
    real, pointer, dimension(:,:) :: hetrores_yr_t
    real, pointer, dimension(:,:) :: nbp_yr_t
    real, pointer, dimension(:,:) :: emit_co2_yr_t
    real, pointer, dimension(:,:) :: emit_co_yr_t
    real, pointer, dimension(:,:) :: emit_ch4_yr_t
    real, pointer, dimension(:,:) :: emit_nmhc_yr_t
    real, pointer, dimension(:,:) :: emit_h2_yr_t
    real, pointer, dimension(:,:) :: emit_nox_yr_t
    real, pointer, dimension(:,:) :: emit_n2o_yr_t
    real, pointer, dimension(:,:) :: emit_pm25_yr_t
    real, pointer, dimension(:,:) :: emit_tpm_yr_t
    real, pointer, dimension(:,:) :: emit_tc_yr_t
    real, pointer, dimension(:,:) :: emit_oc_yr_t
    real, pointer, dimension(:,:) :: emit_bc_yr_t
    real, pointer, dimension(:,:) :: burnfrac_yr_t
    real, pointer, dimension(:,:) :: smfuncveg_yr_t
    real, pointer, dimension(:,:) :: bterm_yr_t
    real, pointer, dimension(:,:) :: luc_emc_yr_t
    real, pointer, dimension(:,:) :: lterm_yr_t
    real, pointer, dimension(:,:) :: lucsocin_yr_t
    real, pointer, dimension(:,:) :: mterm_yr_t
    real, pointer, dimension(:,:) :: lucltrin_yr_t
    real, pointer, dimension(:,:) :: ch4wet1_yr_t
    real, pointer, dimension(:,:) :: ch4wet2_yr_t
    real, pointer, dimension(:,:) :: wetfdyn_yr_t
    real, pointer, dimension(:,:) :: ch4dyn1_yr_t
    real, pointer, dimension(:,:) :: ch4dyn2_yr_t
    real, pointer, dimension(:,:) :: ch4soills_yr_t
    real, pointer, dimension(:,:) :: veghght_yr_t
    real, pointer, dimension(:,:) :: peatdep_yr_t


    logical, pointer, dimension(:,:,:) :: pftexistrow
    real, pointer, dimension(:,:,:) :: gppvegrow
    real, pointer, dimension(:,:,:) :: nepvegrow
    real, pointer, dimension(:,:,:) :: nbpvegrow
    real, pointer, dimension(:,:,:) :: nppvegrow
    real, pointer, dimension(:,:,:) :: hetroresvegrow
    real, pointer, dimension(:,:,:) :: autoresvegrow
    real, pointer, dimension(:,:,:) :: litresvegrow
    real, pointer, dimension(:,:,:) :: soilcresvegrow
    real, pointer, dimension(:,:,:) :: rmlvegaccrow
    real, pointer, dimension(:,:,:) :: rmsvegrow
    real, pointer, dimension(:,:,:) :: rmrvegrow
    real, pointer, dimension(:,:,:) :: rgvegrow
    real, pointer, dimension(:,:,:) :: ailcgrow
    real, pointer, dimension(:,:,:) :: emit_co2row
    real, pointer, dimension(:,:,:) :: emit_corow
    real, pointer, dimension(:,:,:) :: emit_ch4row
    real, pointer, dimension(:,:,:) :: emit_nmhcrow
    real, pointer, dimension(:,:,:) :: emit_h2row
    real, pointer, dimension(:,:,:) :: emit_noxrow
    real, pointer, dimension(:,:,:) :: emit_n2orow
    real, pointer, dimension(:,:,:) :: emit_pm25row
    real, pointer, dimension(:,:,:) :: emit_tpmrow
    real, pointer, dimension(:,:,:) :: emit_tcrow
    real, pointer, dimension(:,:,:) :: emit_ocrow
    real, pointer, dimension(:,:,:) :: emit_bcrow
    real, pointer, dimension(:,:) :: burnfracrow
    real, pointer, dimension(:,:,:) :: burnvegfrow
    real, pointer, dimension(:,:,:) :: smfuncvegrow
    real, pointer, dimension(:,:,:) :: btermrow
    real, pointer, dimension(:,:) :: ltermrow
    real, pointer, dimension(:,:,:) :: mtermrow
    real, pointer, dimension(:,:) :: lucemcomrow
    real, pointer, dimension(:,:) :: lucltrinrow
    real, pointer, dimension(:,:) :: lucsocinrow
    real, pointer, dimension(:,:) :: ch4wet1row
    real, pointer, dimension(:,:) :: ch4wet2row
    real, pointer, dimension(:,:) :: wetfdynrow
    real, pointer, dimension(:,:) :: ch4dyn1row
    real, pointer, dimension(:,:) :: ch4dyn2row
    real, pointer, dimension(:,:) :: ch4soillsrow
    real, pointer, dimension(:,:,:) :: litrmassrow
    real, pointer, dimension(:,:,:) :: soilcmasrow
    real, pointer, dimension(:,:,:) :: vgbiomas_vegrow
    real, pointer, dimension(:,:,:) :: stemmassrow
    real, pointer, dimension(:,:,:) :: rootmassrow
    real, pointer, dimension(:,:,:) :: fcancmxrow
    real, pointer, dimension(:,:,:) :: veghghtrow

    real, pointer, dimension(:,:) :: peatdeprow

    real, pointer, dimension(:) :: laimaxg_yr_g
    real, pointer, dimension(:) :: stemmass_yr_g
    real, pointer, dimension(:) :: rootmass_yr_g
    real, pointer, dimension(:) :: litrmass_yr_g
    real, pointer, dimension(:) :: soilcmas_yr_g
    real, pointer, dimension(:) :: npp_yr_g
    real, pointer, dimension(:) :: gpp_yr_g
    real, pointer, dimension(:) :: nep_yr_g
    real, pointer, dimension(:) :: nbp_yr_g
    real, pointer, dimension(:) :: hetrores_yr_g
    real, pointer, dimension(:) :: autores_yr_g
    real, pointer, dimension(:) :: litres_yr_g
    real, pointer, dimension(:) :: soilcres_yr_g
    real, pointer, dimension(:) :: vgbiomas_yr_g
    real, pointer, dimension(:) :: totcmass_yr_g
    real, pointer, dimension(:) :: emit_co2_yr_g
    real, pointer, dimension(:) :: emit_co_yr_g
    real, pointer, dimension(:) :: emit_ch4_yr_g
    real, pointer, dimension(:) :: emit_nmhc_yr_g
    real, pointer, dimension(:) :: emit_h2_yr_g
    real, pointer, dimension(:) :: emit_nox_yr_g
    real, pointer, dimension(:) :: emit_n2o_yr_g
    real, pointer, dimension(:) :: emit_pm25_yr_g
    real, pointer, dimension(:) :: emit_tpm_yr_g
    real, pointer, dimension(:) :: emit_tc_yr_g
    real, pointer, dimension(:) :: emit_oc_yr_g
    real, pointer, dimension(:) :: emit_bc_yr_g
    real, pointer, dimension(:) :: smfuncveg_yr_g
    real, pointer, dimension(:) :: luc_emc_yr_g
    real, pointer, dimension(:) :: lucltrin_yr_g
    real, pointer, dimension(:) :: lucsocin_yr_g
    real, pointer, dimension(:) :: burnfrac_yr_g
    real, pointer, dimension(:) :: bterm_yr_g
    real, pointer, dimension(:) :: lterm_yr_g
    real, pointer, dimension(:) :: mterm_yr_g
    real, pointer, dimension(:) :: ch4wet1_yr_g
    real, pointer, dimension(:) :: ch4wet2_yr_g
    real, pointer, dimension(:) :: wetfdyn_yr_g
    real, pointer, dimension(:) :: ch4dyn1_yr_g
    real, pointer, dimension(:) :: ch4dyn2_yr_g
    real, pointer, dimension(:) :: ch4soills_yr_g
    real, pointer, dimension(:) :: veghght_yr_g
    real, pointer, dimension(:) :: peatdep_yr_g

    ! local
    integer :: i,m,j,nt
    real :: barefrac
    real :: sumfare
    real, dimension(1) :: timeStamp
    real, dimension(icc) :: pftExist

    ! point pointers

    dofire                => c_switch%dofire
    lnduseon              => c_switch%lnduseon
    PFTCompetition        => c_switch%PFTCompetition
    dowetlands            => c_switch%dowetlands
    obswetf               => c_switch%obswetf
    doperpftoutput        => c_switch%doperpftoutput
    dopertileoutput       => c_switch%dopertileoutput

    laimaxg_yr          =>ctem_yr%laimaxg_yr
    stemmass_yr         =>ctem_yr%stemmass_yr
    rootmass_yr         =>ctem_yr%rootmass_yr
    npp_yr              =>ctem_yr%npp_yr
    gpp_yr              =>ctem_yr%gpp_yr
    vgbiomas_yr         =>ctem_yr%vgbiomas_yr
    autores_yr          =>ctem_yr%autores_yr
    totcmass_yr         =>ctem_yr%totcmass_yr
    litrmass_yr         =>ctem_yr%litrmass_yr
    soilcmas_yr         =>ctem_yr%soilcmas_yr
    nep_yr              =>ctem_yr%nep_yr
    litres_yr           =>ctem_yr%litres_yr
    soilcres_yr         =>ctem_yr%soilcres_yr
    hetrores_yr         =>ctem_yr%hetrores_yr
    nbp_yr              =>ctem_yr%nbp_yr
    emit_co2_yr         =>ctem_yr%emit_co2_yr
    emit_co_yr          =>ctem_yr%emit_co_yr
    emit_ch4_yr         =>ctem_yr%emit_ch4_yr
    emit_nmhc_yr        =>ctem_yr%emit_nmhc_yr
    emit_h2_yr          =>ctem_yr%emit_h2_yr
    emit_nox_yr         =>ctem_yr%emit_nox_yr
    emit_n2o_yr         =>ctem_yr%emit_n2o_yr
    emit_pm25_yr        =>ctem_yr%emit_pm25_yr
    emit_tpm_yr         =>ctem_yr%emit_tpm_yr
    emit_tc_yr          =>ctem_yr%emit_tc_yr
    emit_oc_yr          =>ctem_yr%emit_oc_yr
    emit_bc_yr          =>ctem_yr%emit_bc_yr
    bterm_yr            =>ctem_yr%bterm_yr
    mterm_yr            =>ctem_yr%mterm_yr
    burnfrac_yr         =>ctem_yr%burnfrac_yr
    smfuncveg_yr        =>ctem_yr%smfuncveg_yr
    veghght_yr          =>ctem_yr%veghght_yr

    laimaxg_yr_t          =>ctem_tile_yr%laimaxg_yr_t
    stemmass_yr_t         =>ctem_tile_yr%stemmass_yr_t
    rootmass_yr_t         =>ctem_tile_yr%rootmass_yr_t
    npp_yr_t              =>ctem_tile_yr%npp_yr_t
    gpp_yr_t              =>ctem_tile_yr%gpp_yr_t
    vgbiomas_yr_t         =>ctem_tile_yr%vgbiomas_yr_t
    autores_yr_t          =>ctem_tile_yr%autores_yr_t
    totcmass_yr_t         =>ctem_tile_yr%totcmass_yr_t
    litrmass_yr_t         =>ctem_tile_yr%litrmass_yr_t
    soilcmas_yr_t         =>ctem_tile_yr%soilcmas_yr_t
    nep_yr_t              =>ctem_tile_yr%nep_yr_t
    litres_yr_t           =>ctem_tile_yr%litres_yr_t
    soilcres_yr_t         =>ctem_tile_yr%soilcres_yr_t
    hetrores_yr_t         =>ctem_tile_yr%hetrores_yr_t
    nbp_yr_t              =>ctem_tile_yr%nbp_yr_t
    emit_co2_yr_t         =>ctem_tile_yr%emit_co2_yr_t
    emit_co_yr_t          =>ctem_tile_yr%emit_co_yr_t
    emit_ch4_yr_t         =>ctem_tile_yr%emit_ch4_yr_t
    emit_nmhc_yr_t        =>ctem_tile_yr%emit_nmhc_yr_t
    emit_h2_yr_t          =>ctem_tile_yr%emit_h2_yr_t
    emit_nox_yr_t         =>ctem_tile_yr%emit_nox_yr_t
    emit_n2o_yr_t         =>ctem_tile_yr%emit_n2o_yr_t
    emit_pm25_yr_t        =>ctem_tile_yr%emit_pm25_yr_t
    emit_tpm_yr_t         =>ctem_tile_yr%emit_tpm_yr_t
    emit_tc_yr_t          =>ctem_tile_yr%emit_tc_yr_t
    emit_oc_yr_t          =>ctem_tile_yr%emit_oc_yr_t
    emit_bc_yr_t          =>ctem_tile_yr%emit_bc_yr_t
    burnfrac_yr_t         =>ctem_tile_yr%burnfrac_yr_t
    smfuncveg_yr_t        =>ctem_tile_yr%smfuncveg_yr_t
    bterm_yr_t            =>ctem_tile_yr%bterm_yr_t
    luc_emc_yr_t          =>ctem_tile_yr%luc_emc_yr_t
    lterm_yr_t            =>ctem_tile_yr%lterm_yr_t
    lucsocin_yr_t         =>ctem_tile_yr%lucsocin_yr_t
    mterm_yr_t            =>ctem_tile_yr%mterm_yr_t
    lucltrin_yr_t         =>ctem_tile_yr%lucltrin_yr_t
    ch4wet1_yr_t          =>ctem_tile_yr%ch4wet1_yr_t
    ch4wet2_yr_t          =>ctem_tile_yr%ch4wet2_yr_t
    wetfdyn_yr_t          =>ctem_tile_yr%wetfdyn_yr_t
    ch4dyn1_yr_t          =>ctem_tile_yr%ch4dyn1_yr_t
    ch4dyn2_yr_t          =>ctem_tile_yr%ch4dyn2_yr_t
    ch4soills_yr_t        =>ctem_tile_yr%ch4soills_yr_t
    veghght_yr_t          =>ctem_tile_yr%veghght_yr_t
    peatdep_yr_t          =>ctem_tile_yr%peatdep_yr_t


    pftexistrow       => vrot%pftexist
    gppvegrow         => vrot%gppveg
    nepvegrow         => vrot%nepveg
    nbpvegrow         => vrot%nbpveg
    nppvegrow         => vrot%nppveg
    hetroresvegrow    => vrot%hetroresveg
    autoresvegrow     => vrot%autoresveg
    litresvegrow      => vrot%litresveg
    soilcresvegrow    => vrot%soilcresveg
    rmlvegaccrow      => vrot%rmlvegacc
    rmsvegrow         => vrot%rmsveg
    rmrvegrow         => vrot%rmrveg
    rgvegrow          => vrot%rgveg
    ailcgrow          => vrot%ailcg
    emit_co2row       => vrot%emit_co2
    emit_corow        => vrot%emit_co
    emit_ch4row       => vrot%emit_ch4
    emit_nmhcrow      => vrot%emit_nmhc
    emit_h2row        => vrot%emit_h2
    emit_noxrow       => vrot%emit_nox
    emit_n2orow       => vrot%emit_n2o
    emit_pm25row      => vrot%emit_pm25
    emit_tpmrow       => vrot%emit_tpm
    emit_tcrow        => vrot%emit_tc
    emit_ocrow        => vrot%emit_oc
    emit_bcrow        => vrot%emit_bc
    burnfracrow       => vrot%burnfrac
    burnvegfrow       => vrot%burnvegf
    smfuncvegrow      => vrot%smfuncveg
    btermrow          => vrot%bterm
    ltermrow          => vrot%lterm
    mtermrow          => vrot%mterm
    lucemcomrow       => vrot%lucemcom
    lucltrinrow       => vrot%lucltrin
    lucsocinrow       => vrot%lucsocin
    ch4wet1row        => vrot%ch4wet1
    ch4wet2row        => vrot%ch4wet2
    wetfdynrow        => vrot%wetfdyn
    ch4dyn1row        => vrot%ch4dyn1
    ch4dyn2row        => vrot%ch4dyn2
    ch4soillsrow      => vrot%ch4_soills
    litrmassrow       => vrot%litrmass
    soilcmasrow       => vrot%soilcmas
    vgbiomas_vegrow   => vrot%vgbiomas_veg
    stemmassrow       => vrot%stemmass
    rootmassrow       => vrot%rootmass
    fcancmxrow        => vrot%fcancmx
    veghghtrow        => vrot%veghght

    peatdeprow            => vrot%peatdep

    laimaxg_yr_g          =>ctem_grd_yr%laimaxg_yr_g
    stemmass_yr_g         =>ctem_grd_yr%stemmass_yr_g
    rootmass_yr_g         =>ctem_grd_yr%rootmass_yr_g
    litrmass_yr_g         =>ctem_grd_yr%litrmass_yr_g
    soilcmas_yr_g         =>ctem_grd_yr%soilcmas_yr_g
    npp_yr_g              =>ctem_grd_yr%npp_yr_g
    gpp_yr_g              =>ctem_grd_yr%gpp_yr_g
    nep_yr_g              =>ctem_grd_yr%nep_yr_g
    nbp_yr_g              =>ctem_grd_yr%nbp_yr_g
    hetrores_yr_g         =>ctem_grd_yr%hetrores_yr_g
    autores_yr_g          =>ctem_grd_yr%autores_yr_g
    litres_yr_g           =>ctem_grd_yr%litres_yr_g
    soilcres_yr_g         =>ctem_grd_yr%soilcres_yr_g
    vgbiomas_yr_g         =>ctem_grd_yr%vgbiomas_yr_g
    totcmass_yr_g         =>ctem_grd_yr%totcmass_yr_g
    emit_co2_yr_g         =>ctem_grd_yr%emit_co2_yr_g
    emit_co_yr_g          =>ctem_grd_yr%emit_co_yr_g
    emit_ch4_yr_g         =>ctem_grd_yr%emit_ch4_yr_g
    emit_nmhc_yr_g        =>ctem_grd_yr%emit_nmhc_yr_g
    emit_h2_yr_g          =>ctem_grd_yr%emit_h2_yr_g
    emit_nox_yr_g         =>ctem_grd_yr%emit_nox_yr_g
    emit_n2o_yr_g         =>ctem_grd_yr%emit_n2o_yr_g
    emit_pm25_yr_g        =>ctem_grd_yr%emit_pm25_yr_g
    emit_tpm_yr_g         =>ctem_grd_yr%emit_tpm_yr_g
    emit_tc_yr_g          =>ctem_grd_yr%emit_tc_yr_g
    emit_oc_yr_g          =>ctem_grd_yr%emit_oc_yr_g
    emit_bc_yr_g          =>ctem_grd_yr%emit_bc_yr_g
    smfuncveg_yr_g        =>ctem_grd_yr%smfuncveg_yr_g
    luc_emc_yr_g          =>ctem_grd_yr%luc_emc_yr_g
    lucltrin_yr_g         =>ctem_grd_yr%lucltrin_yr_g
    lucsocin_yr_g         =>ctem_grd_yr%lucsocin_yr_g
    burnfrac_yr_g         =>ctem_grd_yr%burnfrac_yr_g
    bterm_yr_g            =>ctem_grd_yr%bterm_yr_g
    lterm_yr_g            =>ctem_grd_yr%lterm_yr_g
    mterm_yr_g            =>ctem_grd_yr%mterm_yr_g
    ch4wet1_yr_g          =>ctem_grd_yr%ch4wet1_yr_g
    ch4wet2_yr_g          =>ctem_grd_yr%ch4wet2_yr_g
    wetfdyn_yr_g          =>ctem_grd_yr%wetfdyn_yr_g
    ch4dyn1_yr_g          =>ctem_grd_yr%ch4dyn1_yr_g
    ch4dyn2_yr_g          =>ctem_grd_yr%ch4dyn2_yr_g
    ch4soills_yr_g        =>ctem_grd_yr%ch4soills_yr_g
    veghght_yr_g          =>ctem_grd_yr%veghght_yr_g
    peatdep_yr_g          =>ctem_grd_yr%peatdep_yr_g
    !------------

    !> Accumulate yearly outputs

    i = 1 ! offline nlat is always 1 so this array position is always 1.
    do 883 m=1,nmtest
        do 884 j=1,icc

            !> Accumulate the variables at the per PFT level

            if (ailcgrow(i,m,j).gt.laimaxg_yr(i,m,j)) then
                laimaxg_yr(i,m,j)=ailcgrow(i,m,j)
            end if

            npp_yr(i,m,j)=npp_yr(i,m,j)+nppvegrow(i,m,j)
            gpp_yr(i,m,j)=gpp_yr(i,m,j)+gppvegrow(i,m,j)
            nep_yr(i,m,j)=nep_yr(i,m,j)+nepvegrow(i,m,j)
            nbp_yr(i,m,j)=nbp_yr(i,m,j)+nbpvegrow(i,m,j)
            emit_co2_yr(i,m,j)=emit_co2_yr(i,m,j)+emit_co2row(i,m,j)
            emit_co_yr(i,m,j)=emit_co_yr(i,m,j)+emit_corow(i,m,j)
            emit_ch4_yr(i,m,j)=emit_ch4_yr(i,m,j)+emit_ch4row(i,m,j)
            emit_nmhc_yr(i,m,j)=emit_nmhc_yr(i,m,j)+emit_nmhcrow(i,m,j)
            emit_h2_yr(i,m,j)=emit_h2_yr(i,m,j)+emit_h2row(i,m,j)
            emit_nox_yr(i,m,j)=emit_nox_yr(i,m,j)+emit_noxrow(i,m,j)
            emit_n2o_yr(i,m,j)=emit_n2o_yr(i,m,j)+emit_n2orow(i,m,j)
            emit_pm25_yr(i,m,j)=emit_pm25_yr(i,m,j)+emit_pm25row(i,m,j)
            emit_tpm_yr(i,m,j)=emit_tpm_yr(i,m,j)+emit_tpmrow(i,m,j)
            emit_tc_yr(i,m,j)=emit_tc_yr(i,m,j)+emit_tcrow(i,m,j)
            emit_oc_yr(i,m,j)=emit_oc_yr(i,m,j)+emit_ocrow(i,m,j)
            emit_bc_yr(i,m,j)=emit_bc_yr(i,m,j)+emit_bcrow(i,m,j)

            bterm_yr(i,m,j)=bterm_yr(i,m,j)+(btermrow(i,m,j)*(1./real(lastDOY)))
            mterm_yr(i,m,j)=mterm_yr(i,m,j)+(mtermrow(i,m,j)*(1./real(lastDOY)))
            smfuncveg_yr(i,m,j)=smfuncveg_yr(i,m,j)+(smfuncvegrow(i,m,j) * (1./real(lastDOY)))
            hetrores_yr(i,m,j)=hetrores_yr(i,m,j)+hetroresvegrow(i,m,j)
            autores_yr(i,m,j)=autores_yr(i,m,j)+autoresvegrow(i,m,j)
            litres_yr(i,m,j)=litres_yr(i,m,j)+litresvegrow(i,m,j)
            soilcres_yr(i,m,j)=soilcres_yr(i,m,j)+soilcresvegrow(i,m,j)
            burnfrac_yr(i,m,j)=burnfrac_yr(i,m,j)+burnvegfrow(i,m,j)

884         continue

    !>   Also do the bare fraction amounts
    hetrores_yr(i,m,iccp1)=hetrores_yr(i,m,iccp1)+hetroresvegrow(i,m,iccp1)
    litres_yr(i,m,iccp1)=litres_yr(i,m,iccp1)+litresvegrow(i,m,iccp1)
    soilcres_yr(i,m,iccp1)=soilcres_yr(i,m,iccp1)+soilcresvegrow(i,m,iccp1)
    nep_yr(i,m,iccp1)=nep_yr(i,m,iccp1)+nepvegrow(i,m,iccp1)
    nbp_yr(i,m,iccp1)=nbp_yr(i,m,iccp1)+nbpvegrow(i,m,iccp1)

    peatdep_yr_t(i,m)=peatdeprow(i,m)      !YW September 04, 2015


    !> Accumulate the variables at the per tile level
    lterm_yr_t(i,m)=lterm_yr_t(i,m)+(ltermrow(i,m)*(1./real(lastDOY)))
    wetfdyn_yr_t(i,m) = wetfdyn_yr_t(i,m)+(wetfdynrow(i,m)*(1./real(lastDOY)))
    luc_emc_yr_t(i,m)=luc_emc_yr_t(i,m)+lucemcomrow(i,m)
    lucsocin_yr_t(i,m)=lucsocin_yr_t(i,m)+lucsocinrow(i,m)
    lucltrin_yr_t(i,m)=lucltrin_yr_t(i,m)+lucltrinrow(i,m)
    ch4wet1_yr_t(i,m) = ch4wet1_yr_t(i,m)+ch4wet1row(i,m)
    ch4wet2_yr_t(i,m) = ch4wet2_yr_t(i,m)+ch4wet2row(i,m)
    ch4dyn1_yr_t(i,m) = ch4dyn1_yr_t(i,m)+ch4dyn1row(i,m)
    ch4dyn2_yr_t(i,m) = ch4dyn2_yr_t(i,m)+ch4dyn2row(i,m)
    ch4soills_yr_t(i,m) = ch4soills_yr_t(i,m)+ch4soillsrow(i,m)

883     continue ! m

    if (iday.eq.lastDOY) then

        do 900 m = 1, nmtest

            do 925 j=1,icc

                !> The pools are looked at just at the end of the year.
                stemmass_yr(i,m,j)=stemmassrow(i,m,j)
                rootmass_yr(i,m,j)=rootmassrow(i,m,j)
                litrmass_yr(i,m,j)=litrmassrow(i,m,j)
                soilcmas_yr(i,m,j)=soilcmasrow(i,m,j)
                vgbiomas_yr(i,m,j)=vgbiomas_vegrow(i,m,j)
                totcmass_yr(i,m,j)=vgbiomas_yr(i,m,j)+litrmass_yr(i,m,j)+soilcmas_yr(i,m,j)
                veghght_yr(i,m,j)=veghghtrow(i,m,j)

925             continue


            peatdep_yr_g(i)=peatdep_yr_g(i)+peatdep_yr_t(i,m)*farerot(i,m)    !YW September 04, 2015
            litrmass_yr(i,m,iccp1)=litrmassrow(i,m,iccp1)
            soilcmas_yr(i,m,iccp1)=soilcmasrow(i,m,iccp1)
            totcmass_yr(i,m,iccp1)=litrmassrow(i,m,iccp1) + soilcmasrow(i,m,iccp1)


            barefrac=1.0

            !> Add values to the per tile vars
            do j=1,icc

                laimaxg_yr_t(i,m)=laimaxg_yr_t(i,m)+ laimaxg_yr(i,m,j)*fcancmxrow(i,m,j)
                stemmass_yr_t(i,m)=stemmass_yr_t(i,m)+stemmass_yr(i,m,j)*fcancmxrow(i,m,j)
                rootmass_yr_t(i,m)=rootmass_yr_t(i,m)+rootmass_yr(i,m,j)*fcancmxrow(i,m,j)
                litrmass_yr_t(i,m)=litrmass_yr_t(i,m)+litrmass_yr(i,m,j)*fcancmxrow(i,m,j)
                soilcmas_yr_t(i,m)=soilcmas_yr_t(i,m)+soilcmas_yr(i,m,j)*fcancmxrow(i,m,j)
                vgbiomas_yr_t(i,m)=vgbiomas_yr_t(i,m)+vgbiomas_yr(i,m,j)*fcancmxrow(i,m,j)
                totcmass_yr_t(i,m)=totcmass_yr_t(i,m)+totcmass_yr(i,m,j)*fcancmxrow(i,m,j)
                npp_yr_t(i,m)=npp_yr_t(i,m)+npp_yr(i,m,j)*fcancmxrow(i,m,j)
                gpp_yr_t(i,m)=gpp_yr_t(i,m)+gpp_yr(i,m,j)*fcancmxrow(i,m,j)
                nep_yr_t(i,m)=nep_yr_t(i,m)+nep_yr(i,m,j)*fcancmxrow(i,m,j)
                nbp_yr_t(i,m)=nbp_yr_t(i,m)+nbp_yr(i,m,j)*fcancmxrow(i,m,j)
                emit_co2_yr_t(i,m)=emit_co2_yr_t(i,m)+emit_co2_yr(i,m,j)*fcancmxrow(i,m,j)
                emit_co_yr_t(i,m)=emit_co_yr_t(i,m)+emit_co_yr(i,m,j)*fcancmxrow(i,m,j)
                emit_ch4_yr_t(i,m)=emit_ch4_yr_t(i,m)+emit_ch4_yr(i,m,j)*fcancmxrow(i,m,j)
                emit_nmhc_yr_t(i,m)=emit_nmhc_yr_t(i,m)+emit_nmhc_yr(i,m,j)*fcancmxrow(i,m,j)
                emit_h2_yr_t(i,m)=emit_h2_yr_t(i,m)+emit_h2_yr(i,m,j)*fcancmxrow(i,m,j)
                emit_nox_yr_t(i,m)=emit_nox_yr_t(i,m)+emit_nox_yr(i,m,j)*fcancmxrow(i,m,j)
                emit_n2o_yr_t(i,m)=emit_n2o_yr_t(i,m)+emit_n2o_yr(i,m,j)*fcancmxrow(i,m,j)
                emit_pm25_yr_t(i,m)=emit_pm25_yr_t(i,m)+emit_pm25_yr(i,m,j)*fcancmxrow(i,m,j)
                emit_tpm_yr_t(i,m)=emit_tpm_yr_t(i,m)+emit_tpm_yr(i,m,j)*fcancmxrow(i,m,j)
                emit_tc_yr_t(i,m)=emit_tc_yr_t(i,m)+emit_tc_yr(i,m,j)*fcancmxrow(i,m,j)
                emit_oc_yr_t(i,m)=emit_oc_yr_t(i,m)+emit_oc_yr(i,m,j)*fcancmxrow(i,m,j)
                emit_bc_yr_t(i,m)=emit_bc_yr_t(i,m)+emit_bc_yr(i,m,j)*fcancmxrow(i,m,j)
                bterm_yr_t(i,m)=bterm_yr_t(i,m)+bterm_yr(i,m,j)*fcancmxrow(i,m,j)
                mterm_yr_t(i,m)=mterm_yr_t(i,m)+mterm_yr(i,m,j)*fcancmxrow(i,m,j)
                smfuncveg_yr_t(i,m)=smfuncveg_yr_t(i,m)+smfuncveg_yr(i,m,j)*fcancmxrow(i,m,j)
                hetrores_yr_t(i,m)=hetrores_yr_t(i,m)+hetrores_yr(i,m,j)*fcancmxrow(i,m,j)
                autores_yr_t(i,m) =autores_yr_t(i,m) +autores_yr(i,m,j)*fcancmxrow(i,m,j)
                litres_yr_t(i,m)  =litres_yr_t(i,m)  +litres_yr(i,m,j)*fcancmxrow(i,m,j)
                soilcres_yr_t(i,m) =soilcres_yr_t(i,m) +soilcres_yr(i,m,j)*fcancmxrow(i,m,j)
                burnfrac_yr_t(i,m)=burnfrac_yr_t(i,m)+burnfrac_yr(i,m,j)*fcancmxrow(i,m,j)
                veghght_yr_t(i,m) = veghght_yr_t(i,m)+veghght_yr(i,m,j)*fcancmxrow(i,m,j)

                barefrac=barefrac-fcancmxrow(i,m,j)

            end do !j

            litrmass_yr_t(i,m)=litrmass_yr_t(i,m)+litrmass_yr(i,m,iccp1)*barefrac
            soilcmas_yr_t(i,m)=soilcmas_yr_t(i,m)+soilcmas_yr(i,m,iccp1)*barefrac
            hetrores_yr_t(i,m)=hetrores_yr_t(i,m)+hetrores_yr(i,m,iccp1)*barefrac
            litres_yr_t(i,m)  =litres_yr_t(i,m)  +litres_yr(i,m,iccp1)*barefrac
            soilcres_yr_t(i,m)=soilcres_yr_t(i,m)+soilcres_yr(i,m,iccp1)*barefrac
            nep_yr_t(i,m)=nep_yr_t(i,m)+nep_yr(i,m,iccp1)*barefrac
            nbp_yr_t(i,m)=nbp_yr_t(i,m)+nbp_yr(i,m,iccp1)*barefrac
            totcmass_yr_t(i,m) = totcmass_yr_t(i,m)+(litrmass_yr(i,m,iccp1) + soilcmas_yr(i,m,iccp1))*barefrac

            !> Add values to the per gridcell vars
            laimaxg_yr_g(i)=laimaxg_yr_g(i)+ laimaxg_yr_t(i,m)*FAREROT(i,m)
            stemmass_yr_g(i)=stemmass_yr_g(i)+stemmass_yr_t(i,m)*FAREROT(i,m)
            rootmass_yr_g(i)=rootmass_yr_g(i)+rootmass_yr_t(i,m)*FAREROT(i,m)
            litrmass_yr_g(i)=litrmass_yr_g(i)+litrmass_yr_t(i,m)*FAREROT(i,m)
            soilcmas_yr_g(i)=soilcmas_yr_g(i)+soilcmas_yr_t(i,m)*FAREROT(i,m)
            vgbiomas_yr_g(i)=vgbiomas_yr_g(i)+vgbiomas_yr_t(i,m)*FAREROT(i,m)
            totcmass_yr_g(i)=totcmass_yr_g(i)+totcmass_yr_t(i,m)*FAREROT(i,m)
            npp_yr_g(i)=npp_yr_g(i)+npp_yr_t(i,m)*FAREROT(i,m)
            gpp_yr_g(i)=gpp_yr_g(i)+gpp_yr_t(i,m)*FAREROT(i,m)
            nep_yr_g(i)=nep_yr_g(i)+nep_yr_t(i,m)*FAREROT(i,m)
            nbp_yr_g(i)=nbp_yr_g(i)+nbp_yr_t(i,m)*FAREROT(i,m)
            emit_co2_yr_g(i)=emit_co2_yr_g(i)+emit_co2_yr_t(i,m)*FAREROT(i,m)
            emit_co_yr_g(i)=emit_co_yr_g(i)+emit_co_yr_t(i,m)*FAREROT(i,m)
            emit_ch4_yr_g(i)=emit_ch4_yr_g(i)+emit_ch4_yr_t(i,m)*FAREROT(i,m)
            emit_nmhc_yr_g(i)=emit_nmhc_yr_g(i)+emit_nmhc_yr_t(i,m)*FAREROT(i,m)
            emit_h2_yr_g(i)=emit_h2_yr_g(i)+emit_h2_yr_t(i,m)*FAREROT(i,m)
            emit_nox_yr_g(i)=emit_nox_yr_g(i)+emit_nox_yr_t(i,m)*FAREROT(i,m)
            emit_n2o_yr_g(i)=emit_n2o_yr_g(i)+emit_n2o_yr_t(i,m)*FAREROT(i,m)
            emit_pm25_yr_g(i)=emit_pm25_yr_g(i)+emit_pm25_yr_t(i,m)*FAREROT(i,m)
            emit_tpm_yr_g(i)=emit_tpm_yr_g(i)+emit_tpm_yr_t(i,m)*FAREROT(i,m)
            emit_tc_yr_g(i)=emit_tc_yr_g(i)+emit_tc_yr_t(i,m)*FAREROT(i,m)
            emit_oc_yr_g(i)=emit_oc_yr_g(i)+emit_oc_yr_t(i,m)*FAREROT(i,m)
            emit_bc_yr_g(i)=emit_bc_yr_g(i)+emit_bc_yr_t(i,m)*FAREROT(i,m)
            hetrores_yr_g(i)=hetrores_yr_g(i)+hetrores_yr_t(i,m)*FAREROT(i,m)
            autores_yr_g(i) =autores_yr_g(i) +autores_yr_t(i,m)*FAREROT(i,m)
            litres_yr_g(i)  =litres_yr_g(i)  +litres_yr_t(i,m)*FAREROT(i,m)
            soilcres_yr_g(i) =soilcres_yr_g(i) +soilcres_yr_t(i,m)*FAREROT(i,m)
            burnfrac_yr_g(i)=burnfrac_yr_g(i)+burnfrac_yr_t(i,m)*FAREROT(i,m)
            smfuncveg_yr_g(i)=smfuncveg_yr_g(i)+smfuncveg_yr_t(i,m)*FAREROT(i,m)
            bterm_yr_g(i)=bterm_yr_g(i)+bterm_yr_t(i,m)*FAREROT(i,m)
            lterm_yr_g(i)=lterm_yr_g(i)+lterm_yr_t(i,m)*FAREROT(i,m)
            mterm_yr_g(i)=mterm_yr_g(i)+mterm_yr_t(i,m)*FAREROT(i,m)
            luc_emc_yr_g(i)=luc_emc_yr_g(i)+luc_emc_yr_t(i,m)*FAREROT(i,m)
            lucsocin_yr_g(i)=lucsocin_yr_g(i)+lucsocin_yr_t(i,m)*FAREROT(i,m)
            lucltrin_yr_g(i)=lucltrin_yr_g(i)+lucltrin_yr_t(i,m)*FAREROT(i,m)
            ch4wet1_yr_g(i) = ch4wet1_yr_g(i)+ch4wet1_yr_t(i,m)*FAREROT(i,m)
            ch4wet2_yr_g(i) = ch4wet2_yr_g(i)+ch4wet2_yr_t(i,m)*FAREROT(i,m)
            wetfdyn_yr_g(i) = wetfdyn_yr_g(i)+wetfdyn_yr_t(i,m)*FAREROT(i,m)
            ch4dyn1_yr_g(i) = ch4dyn1_yr_g(i)+ch4dyn1_yr_t(i,m)*FAREROT(i,m)
            ch4dyn2_yr_g(i) = ch4dyn2_yr_g(i)+ch4dyn2_yr_t(i,m)*FAREROT(i,m)
            ch4soills_yr_g(i) = ch4soills_yr_g(i)+ch4soills_yr_t(i,m)*FAREROT(i,m)
            veghght_yr_g(i) = veghght_yr_g(i) + veghght_yr_t(i,m)*FAREROT(i,m)

900         continue !m

        !>Write to annual output files:

        ! Prepare the timestamp for this year  !FLAG this isn't correct for leap yet. Need to look at yrs before iyear too!
        timeStamp = (iyear - refyr) * lastDOY

        !> First write out the per gridcell values
        call writeOutput1D(lonLocalIndex,latLocalIndex,'laimaxg_yr_g' ,timeStamp,'lai', [laimaxg_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'vgbiomas_yr_g',timeStamp,'cVeg',[vgbiomas_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'stemmass_yr_g',timeStamp,'cStem',[stemmass_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'rootmass_yr_g',timeStamp,'cRoot',[rootmass_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'litrmass_yr_g',timeStamp,'cLitter',[litrmass_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcmas_yr_g',timeStamp,'cSoil',[soilcmas_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'totcmass_yr_g',timeStamp,'cLand',[totcmass_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'npp_yr_g'     ,timeStamp,'npp',[npp_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'gpp_yr_g'     ,timeStamp,'gpp',[gpp_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nep_yr_g'     ,timeStamp,'nep',[nep_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nbp_yr_g'     ,timeStamp,'nbp',[nbp_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'hetrores_yr_g',timeStamp,'rh',[hetrores_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'autores_yr_g' ,timeStamp,'ra',[autores_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'litres_yr_g'  ,timeStamp,'rhLitter',[litres_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcres_yr_g',timeStamp,'rhSoil',[soilcres_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'veghght_yr_g' ,timeStamp,'vegHeight',[veghght_yr_g(i)])
        do m=1,nmtest
            sumfare = 0.0
            do j=1,icc
                sumfare=sumfare+fcancmxrow(i,m,j)
            end do !j
        end do !m
        call writeOutput1D(lonLocalIndex,latLocalIndex,'fcancmxrow_yr_g' ,timeStamp,'landCoverFrac',[fcancmxrow(i,1,1:icc), 1- sumfare]) !flag only set up for one tile!

        if (dofire .or. lnduseon) then
            call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_co2_yr_g' ,timeStamp,'fFire',[emit_co2_yr_g(i)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'burnfrac_yr_g' ,timeStamp,'burntFractionAll',[burnfrac_yr_g(i)*100.])
        end if

        if (PFTCompetition) then
            pftExist = 0.0
            do j=1,icc
                if (pftexistrow(i,1,j)) then
                    pftExist(j) = 1.0
                end if
            end do
            call writeOutput1D(lonLocalIndex,latLocalIndex,'pftexistrow_yr_g' ,timeStamp,'landCoverExist',[pftExist]) !flag only set up for one tile!
        end if
        if (dowetlands .or. obswetf) then
            call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4wet1_yr_g' ,timeStamp,'wetlandCH4spec',[ch4wet1_yr_g(i)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4dyn1_yr_g' ,timeStamp,'wetlandCH4dyn',[ch4dyn1_yr_g(i)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'wetfdyn_yr_g' ,timeStamp,'wetlandFrac',[wetfdyn_yr_g(i)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4soills_yr_g' ,timeStamp,'soilCH4cons',[ch4soills_yr_g(i)])
        end if

        if (doperpftoutput) then
            if (nmtest > 1) then
                print*,'Per PFT and per tile outputs not implemented yet'
            else
                m = 1 !FLAG
                call writeOutput1D(lonLocalIndex,latLocalIndex,'laimaxg_yr' ,timeStamp,'lai', [laimaxg_yr(i,m,:)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'vgbiomas_yr',timeStamp,'cVeg',[vgbiomas_yr(i,m,:)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'stemmass_yr',timeStamp,'cStem',[stemmass_yr(i,m,:)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'rootmass_yr',timeStamp,'cRoot',[rootmass_yr(i,m,:)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'litrmass_yr',timeStamp,'cLitter',[litrmass_yr(i,m,:)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcmas_yr',timeStamp,'cSoil',[soilcmas_yr(i,m,:)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'totcmass_yr',timeStamp,'cLand',[totcmass_yr(i,m,:)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'npp_yr'     ,timeStamp,'npp',[npp_yr(i,m,:)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'gpp_yr'     ,timeStamp,'gpp',[gpp_yr(i,m,:)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'nep_yr'     ,timeStamp,'nep',[nep_yr(i,m,:)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'nbp_yr'     ,timeStamp,'nbp',[nbp_yr(i,m,:)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'hetrores_yr',timeStamp,'rh',[hetrores_yr(i,m,:)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'autores_yr' ,timeStamp,'ra',[autores_yr(i,m,:)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'litres_yr'  ,timeStamp,'rhLitter',[litres_yr(i,m,:)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcres_yr',timeStamp,'rhSoil',[soilcres_yr(i,m,:)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'veghght_yr' ,timeStamp,'vegHeight',[veghght_yr(i,m,:)])

                if (dofire .or. lnduseon) then

                    call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_co2_yr' ,timeStamp,'fFire',[emit_co2_yr(i,m,:)])
                    call writeOutput1D(lonLocalIndex,latLocalIndex,'burnfrac_yr' ,timeStamp,'burntFractionAll',[burnfrac_yr(i,m,:)*100.])
!                            smfuncveg_yr(i,m,j), &
!                            bterm_yr(i,m,j),lterm_yr_t(i,m),mterm_yr(i,m,j), &
                end if

            end if
        end if

        if (dopertileoutput) then

            if (nmtest == 1) then
                print*,'Switch selected for per tile output but number of tiles is only one.'
            else
                !> Write out the per tile values
                call writeOutput1D(lonLocalIndex,latLocalIndex,'laimaxg_yr_t' ,timeStamp,'lai', [laimaxg_yr_t(i,:)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'vgbiomas_yr_t',timeStamp,'cVeg',[vgbiomas_yr_t(i,:)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'stemmass_yr_t',timeStamp,'cStem',[stemmass_yr_t(i,:)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'rootmass_yr_t',timeStamp,'cRoot',[rootmass_yr_t(i,:)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'litrmass_yr_t',timeStamp,'cLitter',[litrmass_yr_t(i,:)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcmas_yr_t',timeStamp,'cSoil',[soilcmas_yr_t(i,:)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'totcmass_yr_t',timeStamp,'cLand',[totcmass_yr_t(i,:)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'npp_yr_t'     ,timeStamp,'npp',[npp_yr_t(i,:)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'gpp_yr_t'     ,timeStamp,'gpp',[gpp_yr_t(i,:)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'nep_yr_t'     ,timeStamp,'nep',[nep_yr_t(i,:)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'nbp_yr_t'     ,timeStamp,'nbp',[nbp_yr_t(i,:)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'hetrores_yr_t',timeStamp,'rh',[hetrores_yr_t(i,:)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'autores_yr_t' ,timeStamp,'ra',[autores_yr_t(i,:)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'litres_yr_t'  ,timeStamp,'rhLitter',[litres_yr_t(i,:)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcres_yr_t',timeStamp,'rhSoil',[soilcres_yr_t(i,:)])
                call writeOutput1D(lonLocalIndex,latLocalIndex,'veghght_yr_t' ,timeStamp,'vegHeight',[veghght_yr_t(i,:)])

                if (dofire .or. lnduseon) then
                    call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_co2_yr_t' ,timeStamp,'fFire',[emit_co2_yr_t(i,:)])
                    call writeOutput1D(lonLocalIndex,latLocalIndex,'burnfrac_yr_t' ,timeStamp,'burntFractionAll',[burnfrac_yr_t(i,:)*100.])
                end if
                if (dowetlands .or. obswetf) then
                    call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4wet1_yr_t' ,timeStamp,'wetlandCH4spec',[ch4wet1_yr_t(i,:)])
                    call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4dyn1_yr_t' ,timeStamp,'wetlandCH4dyn',[ch4dyn1_yr_t(i,:)])
                    call writeOutput1D(lonLocalIndex,latLocalIndex,'wetfdyn_yr_t' ,timeStamp,'wetlandFrac',[wetfdyn_yr_t(i,:)])
                    call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4soills_yr_t' ,timeStamp,'soilCH4cons',[ch4soills_yr_t(i,:)])
                end if
            end if
        end if

        !> Reset all annual vars in preparation for next year
        call resetyearend(nltest,nmtest)

    endif ! if iday=365/366

end subroutine ctem_annual_aw
!>@}
!>\file
!>Central module that handles all CTEM preparation and writing of output files

end module io_driver      
