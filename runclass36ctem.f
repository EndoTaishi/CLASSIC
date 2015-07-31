      PROGRAM RUNCLASS36CTEM
C
C     * SHELL PROGRAM TO RUN "CLASS" ("CANADIAN LAND SURFACE SCHEME")
C     * VERSION 3.6 IN STAND-ALONE MODE USING SPECIFIED BOUNDARY
C     * CONDITIONS AND ATMOSPHERIC FORCING, COUPLED TO CTEM (CANADIAN TERRESTRIAL
C     * ECOSYSTEM MODEL).
C
C     REVISION HISTORY:
C
C     * JUL 2 2015
C     * JOE MELTON : Took many calculations out of this driver and into subroutines. Introduced
C                    modular structure and made CTEM vars into pointers. Harmonized CLASS v. 3.6.1
C                    code with CTEM code and into this driver.
C
C     * JAN 14 2014
C     * JOE MELTON : Harmonized the field capacity and wilting point calculations between CLASS and CTEM.
C                    took the code out of runclassctem and it is now fully done within CLASSB. Harmonized names too.
C
C     * JUN 2014
C     * RUDRA SHRESTHA : ADD IN WETLAND CODE
C
C     * JUL 2013
C     * JOE MELTON : REMOVED CTEM1 AND CTEM2 OPTIONS, REPLACED WITH CTEM_ON. INTRODUCE
C                                  MODULES. RESTRUCTURE OUTPUTS AND CTEM VARIABLE DECLARATIONS
C                                  OUTPUTS ARE NOW THE SAME FOR BOTH MOSAIC AND COMPOSITE MODES
C
C     * DEC 2012
C     * JOE MELTON : REMOVED GOTO STATEMENTS, CLEANED UP AND FIXED INCONSISTENCIES
C                                  IN HOW INPUT DATA READ IN. ALSO MADE LUC WORK FOR
C                                  BOTH COMPOSITE AND MOSAIC APPROACHES.
C
C     * OCT 2012
C     * YIRAN PENG AND JOE MELTON: BRING IN COMPETITION TO 3.6 AND MAKE IT
C                                  SO THE MODEL CAN START FROM BARE GROUND
C                                  OR FROM THE INI AND CTM FILES INPUTS
C
C     * SEP 2012
C     * JOE MELTON: COUPLED CLASS3.6 AND CTEM
C
C     * NOV 2011
C     * YIRAN PENG AND VIVEK ARORA: COUPLED CLASS3.5 AND CTEM
C
C     * SEPT 8, 2009
C     * RONG LI AND VIVEK ARORA: COUPLED CLASS3.4 AND CTEM
C
C=======================================================================

C     * DIMENSION STATEMENTS.

C     * FIRST SET OF DEFINITIONS:
C     * BACKGROUND VARIABLES, AND PROGNOSTIC AND DIAGNOSTIC
C     * VARIABLES NORMALLY PROVIDED BY AND/OR USED BY THE GCM.
C     * THE SUFFIX "ROT" REFERS TO VARIABLES EXISTING ON THE
C     * MOSAIC GRID ON THE CURRENT LATITUDE CIRCLE.  THE SUFFIX
C     * "GAT" REFERS TO THE SAME VARIABLES AFTER THEY HAVE UNDERGONE
C     * A "GATHER" OPERATION IN WHICH THE TWO MOSAIC DIMENSIONS
C     * ARE COLLAPSED INTO ONE.  THE SUFFIX "ROW" REFERS BOTH TO
C     * GRID-CONSTANT INPUT VARIABLES. AND TO GRID-AVERAGED
C     * DIAGNOSTIC VARIABLES.
C
C     * THE FIRST DIMENSION ELEMENT OF THE "ROT" VARIABLES
C     * REFERS TO THE NUMBER OF GRID CELLS ON THE CURRENT
C     * LATITUDE CIRCLE.  IN THIS STAND-ALONE VERSION, THIS
C     * NUMBER IS ARBITRARILY SET TO THREE, TO ALLOW UP TO THREE
C     * SIMULTANEOUS TESTS TO BE RUN.  THE SECOND DIMENSION
C     * ELEMENT OF THE "ROT" VARIABLES REFERS TO THE MAXIMUM
C     * NUMBER OF TILES IN THE MOSAIC.  IN THIS STAND-ALONE
C     * VERSION, THIS NUMBER IS SET TO EIGHT.  THE FIRST
C     * DIMENSION ELEMENT IN THE "GAT" VARIABLES IS GIVEN BY
C     * THE PRODUCT OF THE FIRST TWO DIMENSION ELEMENTS IN THE
C     * "ROT" VARIABLES.

C     The majority of CTEM parameters are stored in ctem_params.f90. We access them
c     through use statements for modules:
      use ctem_params,        only : initpftpars,nlat,nmos,ilg,nmon,
     1                               ican, ignd,icp1, icc, iccp1,
     2                               monthend, mmday,modelpft, l2max,
     3                                deltat, abszero, monthdays,seed,
     4                                crop, NBS

      use landuse_change,     only : initialize_luc, readin_luc

      use ctem_statevars,     only : vrot,vgat,c_switch,initrowvars,
     1                               class_out,resetclassmon,
     2                               resetclassyr,resetmidmonth,
     3                               resetmonthend_g,ctem_grd_mo,
     4                               resetyearend_g,ctem_grd_yr,
     5                               resetclassaccum,ctem_grd,
     6                               ctem_tile, ctem_tile_mo,
     7                               ctem_tile_yr,resetgridavg,
     8                               resetmonthend_m,resetyearend_g,
     9                               resetyearend_m

      use io_driver,          only : read_from_ctm, create_outfiles,
     1                               write_ctm_rs, class_monthly_aw,
     2                               ctem_annual_aw,ctem_monthly_aw,
     3                               close_outfiles,ctem_daily_aw

c
      implicit none
C
C     * INTEGER CONSTANTS.
C
      INTEGER IDISP,IZREF,ISLFD,IPCP,IWF,IPAI,IHGT,IALC,
     1        IALS,IALG,N,ITG,ITC,ITCG,isnoalb,igralb

      INTEGER NLTEST,NMTEST,NCOUNT,NDAY,
     1        IMONTH,NDMONTH,NT,
     2        IHOUR,IMIN,IDAY,IYEAR,NML,NMW,NWAT,NICE,JLAT,
     3        NLANDCS,NLANDGS,NLANDC,NLANDG,NLANDI,I,J,K,L,M,
     4        NTLD
C
      INTEGER K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,ITA,ITCAN,ITD,
     1        ITAC,ITS,ITSCR,ITD2,ITD3,ITD4,ISTEPS,NFS,NDRY,NAL,NFT
      REAL TAHIST(200), TCHIST(200), TACHIST(200), TDHIST(200),
     1     TSHIST(200),TSCRHIST(200),TD2HIST(200),TD3HIST(200),
     2     TD4HIST(200),PAICAN(ILG)

      INTEGER*4 TODAY(3), NOW(3)
C
C     * LAND SURFACE PROGNOSTIC VARIABLES.
C
      REAL,DIMENSION(NLAT,NMOS,IGND) ::
     1        TBARROT,   THLQROT,   THICROT
C
      REAL,DIMENSION(NLAT,NMOS) ::
     1        TPNDROT,   ZPNDROT,   TBASROT,
     2        ALBSROT,   TSNOROT,   RHOSROT,
     3        SNOROT ,   TCANROT,   RCANROT,
     4        SCANROT,   GROROT ,   CMAIROT,
     5        TACROT ,   QACROT ,   WSNOROT,
     6        REFROT,    BCSNROT
C
      REAL    TSFSROT(NLAT,NMOS,4)
C
      REAL,DIMENSION(ILG,IGND) ::
     1        TBARGAT, THLQGAT, THICGAT
C
      REAL,DIMENSION(ILG) ::
     1        TPNDGAT,   ZPNDGAT,   TBASGAT,
     2        ALBSGAT,   TSNOGAT,   RHOSGAT,
     3        SNOGAT ,   TCANGAT,   RCANGAT,
     4        SCANGAT,   GROGAT ,   CMAIGAT,
     5        TACGAT ,   QACGAT ,   WSNOGAT,
     6        REFGAT,    BCSNGAT

C
      REAL    TSFSGAT(ILG,4)
C
C     * GATHER-SCATTER INDEX ARRAYS.
C
      INTEGER  ILMOS (ILG),  JLMOS  (ILG),  IWMOS  (ILG),  JWMOS (ILG)
C
C     * CANOPY AND SOIL INFORMATION ARRAYS.
C     * (THE LAST DIMENSION OF MOST OF THESE ARRAYS IS GIVEN BY
C     * THE NUMBER OF SOIL LAYERS (IGND), THE NUMBER OF BROAD
C     * VEGETATION CATEGORIES (ICAN), OR ICAN+1.
C
      REAL,DIMENSION(NLAT,NMOS,ICP1) ::
     1              FCANROT,  LNZ0ROT,
     2              ALVCROT,  ALICROT
C
      REAL,DIMENSION(NLAT,NMOS,ICAN) ::
     1              PAMXROT,  PAMNROT,
     2              CMASROT,  ROOTROT,
     3              RSMNROT,  QA50ROT,
     4              VPDAROT,  VPDBROT,
     5              PSGAROT,  PSGBROT,
     6              PAIDROT,  HGTDROT,
     7              ACVDROT,  ACIDROT
C
      REAL,DIMENSION(ILG,ICP1) ::
     1              FCANGAT,  LNZ0GAT,
     2              ALVCGAT,  ALICGAT
C
      REAL,DIMENSION(ILG,ICAN) ::
     1              PAMXGAT,  PAMNGAT,
     2              CMASGAT,  ROOTGAT,
     3              RSMNGAT,  QA50GAT,
     4              VPDAGAT,  VPDBGAT,
     5              PSGAGAT,  PSGBGAT,
     6              PAIDGAT,  HGTDGAT,
     7              ACVDGAT,  ACIDGAT
C
      REAL,DIMENSION(NLAT,NMOS,IGND) ::
     1        THPROT ,  THRROT ,  THMROT ,
     2        BIROT  ,  PSISROT,  GRKSROT,
     3        THRAROT,  HCPSROT,
     4        TCSROT ,  THFCROT,  PSIWROT,
     5        THLWROT,  DLZWROT,  ZBTWROT
C
      REAL,DIMENSION(NLAT,NMOS) ::
     1        DRNROT ,   XSLPROT,   GRKFROT,
     2        WFSFROT,   WFCIROT,   ALGWROT,
     3        ALGDROT,   ASVDROT,   ASIDROT,
     4        AGVDROT,   AGIDROT,   ZSNLROT,
     5        ZPLGROT,   ZPLSROT,   ZSNOROT,
     6        ALGWVROT,  ALGWNROT,  ALGDVROT,
     7        ALGDNROT,  EMISROT

C
      REAL,DIMENSION(NLAT,NMOS,NBS) ::
     1        SALBROT,   CSALROT
C
      REAL,DIMENSION(NLAT,NBS) ::
     1        FSDBROL,   FSFBROL,   FSSBROL

C
      REAL,DIMENSION(ILG,IGND) ::
     1        THPGAT ,  THRGAT ,  THMGAT ,
     2        BIGAT  ,  PSISGAT,  GRKSGAT,
     3        THRAGAT,  HCPSGAT,
     4        TCSGAT ,  THFCGAT,  PSIWGAT,
     5        THLWGAT,  DLZWGAT,  ZBTWGAT
C
      REAL,DIMENSION(ILG) ::
     1        DRNGAT ,   XSLPGAT,   GRKFGAT,
     2        WFSFGAT,   WFCIGAT,   ALGWGAT,
     3        ALGDGAT,   ASVDGAT,   ASIDGAT,
     4        AGVDGAT,   AGIDGAT,   ZSNLGAT,
     5        ZPLGGAT,   ZPLSGAT,   ALGWVGAT,
     6        ALGWNGAT,  ALGDVGAT,  ALGDNGAT,
     7        EMISGAT
C
      REAL    SANDROT(NLAT,NMOS,IGND), CLAYROT(NLAT,NMOS,IGND),
     1        ORGMROT(NLAT,NMOS,IGND), SOCIROT(NLAT,NMOS),
     2        SDEPROT(NLAT,NMOS),      FAREROT(NLAT,NMOS)
C
      INTEGER MIDROT (NLAT,NMOS),     ISNDROT(NLAT,NMOS,IGND),
     1        ISNDGAT( ILG,IGND),     IGDRROT(NLAT,NMOS),
     2        IGDRGAT( ILG)
C
      REAL,DIMENSION(ILG,NBS) ::
     1        FSDBGAT,   FSFBGAT,   FSSBGAT,
     2        SALBGAT,   CSALGAT
C
C     * ARRAYS ASSOCIATED WITH COMMON BLOCKS.
C
      REAL  THPORG (  3), THRORG (  3), THMORG (  3), BORG   (  3),
     1      PSISORG(  3), GRKSORG(  3)
C
      REAL  CANEXT(ICAN), XLEAF (ICAN), ZORAT (ICAN),
     1      DELZ  (IGND), ZBOT  (IGND),
     2      GROWYR (  18,4,2)
C
C     * ATMOSPHERIC AND GRID-CONSTANT INPUT VARIABLES.
C
      REAL,DIMENSION(NLAT) ::
     1      ZRFMROW,   ZRFHROW,   ZDMROW ,   ZDHROW ,
     2      ZBLDROW,   FSVHROW,   FSIHROW,   RADJROW,
     3      CSZROW ,   FDLROW ,   ULROW  ,   VLROW  ,
     4      TAROW  ,   QAROW  ,   PRESROW,   PREROW ,
     5      PADRROW,   VPDROW ,   TADPROW,   RHOAROW,
     6      RPCPROW,   TRPCROW,   SPCPROW,   TSPCROW,
     7      RHSIROW,   FCLOROW,   DLONROW,   UVROW  ,
     8      XDIFFUS,   GCROW  ,   Z0ORROW,   GGEOROW,
     9      RPREROW,   SPREROW,   VMODROW,   DLATROW,
     A      FSSROW,    PRENROW,   CLDTROW,   FSGROL,
     B      FLGROL,    GUSTROL,   DEPBROW
C
      REAL,DIMENSION(ILG) ::
     1      ZRFMGAT,   ZRFHGAT,   ZDMGAT ,   ZDHGAT ,
     2      ZBLDGAT,   FSVHGAT,   FSIHGAT,   RADJGAT,
     3      CSZGAT ,   FDLGAT ,   ULGAT  ,   VLGAT  ,
     4      TAGAT  ,   QAGAT  ,   PRESGAT,   PREGAT ,
     5      PADRGAT,   VPDGAT ,   TADPGAT,   RHOAGAT,
     6      RPCPGAT,   TRPCGAT,   SPCPGAT,   TSPCGAT,
     7      RHSIGAT,   FCLOGAT,   DLONGAT,   Z0ORGAT,
     8      GGEOGAT,   VMODGAT,   FSGGAT,    FLGGAT,
     9      GUSTGAT,   DEPBGAT,   GTBS,      SFCUBS,
     +      SFCVBS,    USTARBS,   TCSNOW,    GSNOW

C
C     * LAND SURFACE DIAGNOSTIC VARIABLES.
C
      REAL,DIMENSION(NLAT,NMOS) ::
     1      CDHROT ,   CDMROT ,   HFSROT ,   TFXROT ,
     2      QEVPROT,   QFSROT ,   QFXROT ,   PETROT ,
     3      GAROT  ,   EFROT  ,   GTROT  ,   QGROT  ,
     4      ALVSROT,   ALIRROT,   FSNOROT,   SFRHROT,
     5      SFCTROT,   SFCUROT,   SFCVROT,   SFCQROT,
     6      FSGVROT,   FSGSROT,   FSGGROT,   FLGVROT,
     7      FLGSROT,   FLGGROT,   HFSCROT,   HFSSROT,
     8      HFSGROT,   HEVCROT,   HEVSROT,   HEVGROT,
     9      HMFCROT,   HMFNROT,   HTCCROT,   HTCSROT,
     A      PCFCROT,   PCLCROT,   PCPNROT,   PCPGROT,
     B      QFGROT ,   QFNROT ,   QFCLROT,   QFCFROT,
     C      ROFROT ,   ROFOROT,   ROFSROT,   ROFBROT,
     D      TROFROT,   TROOROT,   TROSROT,   TROBROT,
     E      ROFCROT,   ROFNROT,   ROVGROT,   WTRCROT,
     F      WTRSROT,   WTRGROT,   DRROT  ,   WTABROT,
     G      ILMOROT,   UEROT  ,   HBLROT
C
      REAL,DIMENSION(ILG) ::
     1      CDHGAT ,   CDMGAT ,   HFSGAT ,   TFXGAT ,
     2      QEVPGAT,   QFSGAT ,   QFXGAT ,   PETGAT ,
     3      GAGAT  ,   EFGAT  ,   GTGAT  ,   QGGAT  ,
     4      ALVSGAT,   ALIRGAT,   FSNOGAT,   SFRHGAT,
     5      SFCTGAT,   SFCUGAT,   SFCVGAT,   SFCQGAT,
     6      FSGVGAT,   FSGSGAT,   FSGGGAT,   FLGVGAT,
     7      FLGSGAT,   FLGGGAT,   HFSCGAT,   HFSSGAT,
     8      HFSGGAT,   HEVCGAT,   HEVSGAT,   HEVGGAT,
     9      HMFCGAT,   HMFNGAT,   HTCCGAT,   HTCSGAT,
     A      PCFCGAT,   PCLCGAT,   PCPNGAT,   PCPGGAT,
     B      QFGGAT ,   QFNGAT ,   QFCLGAT,   QFCFGAT,
     C      ROFGAT ,   ROFOGAT,   ROFSGAT,   ROFBGAT,
     D      TROFGAT,   TROOGAT,   TROSGAT,   TROBGAT,
     E      ROFCGAT,   ROFNGAT,   ROVGGAT,   WTRCGAT,
     F      WTRSGAT,   WTRGGAT,   DRGAT  ,   WTABGAT,
     G      ILMOGAT,   UEGAT  ,   HBLGAT ,   QLWOGAT,
     H      FTEMP  ,   FVAP   ,   RIB
C
      REAL,DIMENSION(NLAT) ::
     1      CDHROW ,   CDMROW ,   HFSROW ,   TFXROW ,
     2      QEVPROW,   QFSROW ,   QFXROW ,   PETROW ,
     3      GAROW  ,   EFROW  ,   GTROW  ,   QGROW  ,
     4      ALVSROW,   ALIRROW,   FSNOROW,   SFRHROW,
     5      SFCTROW,   SFCUROW,   SFCVROW,   SFCQROW,
     6      FSGVROW,   FSGSROW,   FSGGROW,   FLGVROW,
     7      FLGSROW,   FLGGROW,   HFSCROW,   HFSSROW,
     8      HFSGROW,   HEVCROW,   HEVSROW,   HEVGROW,
     9      HMFCROW,   HMFNROW,   HTCCROW,   HTCSROW,
     A      PCFCROW,   PCLCROW,   PCPNROW,   PCPGROW,
     B      QFGROW ,   QFNROW ,   QFCLROW,   QFCFROW,
     C      ROFROW ,   ROFOROW,   ROFSROW,   ROFBROW,
     D      ROFCROW,   ROFNROW,   ROVGROW,   WTRCROW,
     E      WTRSROW,   WTRGROW,   DRROW  ,   WTABROW,
     F      ILMOROW,   UEROW  ,   HBLROW
C
      REAL    HMFGROT(NLAT,NMOS,IGND),   HTCROT (NLAT,NMOS,IGND),
     1        QFCROT (NLAT,NMOS,IGND),   GFLXROT(NLAT,NMOS,IGND),
     2        HMFGGAT(ILG,IGND),         HTCGAT (ILG,IGND),
     3        QFCGAT (ILG,IGND),         GFLXGAT(ILG,IGND),
     4        HMFGROW(NLAT,IGND),        HTCROW (NLAT,IGND),
     5        QFCROW (NLAT,IGND),        GFLXROW(NLAT,IGND)
C
      INTEGER     ITCTROT(NLAT,NMOS,6,50),  ITCTGAT(ILG,6,50)
      INTEGER     ISUM(6)

C     * ARRAYS USED FOR OUTPUT AND DISPLAY PURPOSES.
C     * (THE SUFFIX "ACC" REFERS TO ACCUMULATOR ARRAYS USED IN
C     * CALCULATING TIME AVERAGES.)

      CHARACTER     TITLE1*4,     TITLE2*4,     TITLE3*4,
     1              TITLE4*4,     TITLE5*4,     TITLE6*4
      CHARACTER     NAME1*4,      NAME2*4,      NAME3*4,
     1              NAME4*4,      NAME5*4,      NAME6*4
      CHARACTER     PLACE1*4,     PLACE2*4,     PLACE3*4,
     1              PLACE4*4,     PLACE5*4,     PLACE6*4

      REAL,DIMENSION(NLAT) ::
     1              PREACC ,   GTACC  ,   QEVPACC,
     2              HFSACC ,   ROFACC ,   SNOACC ,
     3              ALVSACC,   ALIRACC,   FSINACC,
     4              FLINACC,   TAACC  ,   UVACC  ,
     5              PRESACC,   QAACC  ,
     6              EVAPACC,   FLUTACC,   OVRACC ,
     7              HMFNACC,   WTBLACC,   WSNOACC,
     8              RHOSACC,   TSNOACC,   TCANACC,
     9              RCANACC,   SCANACC,   GROACC ,
     A              CANARE ,   SNOARE

      REAL          TBARACC(NLAT,IGND), THLQACC(NLAT,IGND),
     1              THICACC(NLAT,IGND), THALACC(NLAT,IGND)
C
! C     * MONTHLY OUTPUT FOR CLASS GRID-MEAN
! C
!       REAL,DIMENSION(NLAT) ::
!      1              ALVSACC_MO,ALIRACC_MO,FLUTACC_MO,
!      2              FSINACC_MO,FLINACC_MO,HFSACC_MO,
!      3              QEVPACC_MO,SNOACC_MO, WSNOACC_MO,
!      4              ROFACC_MO, PREACC_MO, EVAPACC_MO,
!      5              TAACC_MO
!       REAL ::       FSSTAR_MO,FLSTAR_MO,QH_MO,QE_MO
!
!       REAL TBARACC_MO(NLAT,IGND), THLQACC_MO(NLAT,IGND),
!      1     THICACC_MO(NLAT,IGND)
! C
! C     * YEARLY OUTPUT FOR CLASS GRID-MEAN
! C
!       REAL,DIMENSION(NLAT) ::
!      1              ALVSACC_YR,ALIRACC_YR,FLUTACC_YR,
!      2              FSINACC_YR,FLINACC_YR,HFSACC_YR,
!      3              QEVPACC_YR,ROFACC_YR, PREACC_YR,
!      4              EVAPACC_YR,TAACC_YR
!       REAL ::       FSSTAR_YR,FLSTAR_YR,QH_YR,QE_YR
C
C     * ARRAYS DEFINED TO PASS INFORMATION BETWEEN THE THREE MAJOR
C     * SUBSECTIONS OF CLASS ("CLASSA", "CLASST" AND "CLASSW").

      REAL,DIMENSION(ILG,IGND) ::
     1        TBARC  ,     TBARG  ,     TBARCS ,
     2        TBARGS ,     THLIQC ,     THLIQG ,
     3        THICEC ,     THICEG ,     FROOT  ,
     4        HCPC   ,     HCPG   ,
     5        TCTOPC ,     TCBOTC ,
     6        TCTOPG ,     TCBOTG
C
      REAL  FC     (ILG), FG     (ILG), FCS    (ILG), FGS    (ILG),
     1      RBCOEF (ILG), ZSNOW  (ILG),
     2      FSVF   (ILG), FSVFS  (ILG),
     3      ALVSCN (ILG), ALIRCN (ILG), ALVSG  (ILG), ALIRG  (ILG),
     4      ALVSCS (ILG), ALIRCS (ILG), ALVSSN (ILG), ALIRSN (ILG),
     5      ALVSGC (ILG), ALIRGC (ILG), ALVSSC (ILG), ALIRSC (ILG),
     6      TRVSCN (ILG), TRIRCN (ILG), TRVSCS (ILG), TRIRCS (ILG),
     7      RC     (ILG), RCS    (ILG), FRAINC (ILG), FSNOWC (ILG),
     8      FRAICS (ILG), FSNOCS (ILG),
     9      CMASSC (ILG), CMASCS (ILG), DISP   (ILG), DISPS  (ILG),
     A      ZOMLNC (ILG), ZOELNC (ILG), ZOMLNG (ILG), ZOELNG (ILG),
     B      ZOMLCS (ILG), ZOELCS (ILG), ZOMLNS (ILG), ZOELNS (ILG),
     C      TRSNOWC (ILG), CHCAP  (ILG), CHCAPS (ILG),
     D      GZEROC (ILG), GZEROG (ILG), GZROCS (ILG), GZROGS (ILG),
     E      G12C   (ILG), G12G   (ILG), G12CS  (ILG), G12GS  (ILG),
     F      G23C   (ILG), G23G   (ILG), G23CS  (ILG), G23GS  (ILG),
     G      QFREZC (ILG), QFREZG (ILG), QMELTC (ILG), QMELTG (ILG),
     I      EVAPC  (ILG), EVAPCG (ILG), EVAPG  (ILG), EVAPCS (ILG),
     J      EVPCSG (ILG), EVAPGS (ILG), TCANO  (ILG), TCANS  (ILG),
     K      RAICAN (ILG), SNOCAN (ILG), RAICNS (ILG), SNOCNS (ILG),
     L      CWLCAP (ILG), CWFCAP (ILG), CWLCPS (ILG), CWFCPS (ILG),
     M      TSNOCS (ILG), TSNOGS (ILG), RHOSCS (ILG), RHOSGS (ILG),
     N      WSNOCS (ILG), WSNOGS (ILG),
     O      TPONDC (ILG), TPONDG (ILG), TPNDCS (ILG), TPNDGS (ILG),
     P      ZPLMCS (ILG), ZPLMGS (ILG), ZPLIMC (ILG), ZPLIMG (ILG)
C
      REAL  ALTG(ILG,NBS),ALSNO(ILG,NBS),TRSNOWG(ILG,NBS)

C
C     * DIAGNOSTIC ARRAYS USED FOR CHECKING ENERGY AND WATER
C     * BALANCES.
C
      REAL CTVSTP(ILG),   CTSSTP(ILG),   CT1STP(ILG),   CT2STP(ILG),
     1     CT3STP(ILG),   WTVSTP(ILG),   WTSSTP(ILG),   WTGSTP(ILG)
C
C     * CONSTANTS AND TEMPORARY VARIABLES.
C
      REAL DEGLON,DAY,DECL,HOUR,COSZ,CUMSNO,EVAPSUM,
     1     QSUMV,QSUMS,QSUM1,QSUM2,QSUM3,WSUMV,WSUMS,WSUMG,ALTOT,
     2     FSSTAR,FLSTAR,QH,QE,BEG,SNOMLT,ZSN,TCN,TSN,TPN,GTOUT,TAC,
     3     ALTOT_YR,TSURF,ALAVG,ALMAX,ACTLYR,FTAVG,FTMAX,FTABLE
C
C     * COMMON BLOCK PARAMETERS.
C
      REAL X1,X2,X3,X4,G,GAS,X5,X6,CPRES,GASV,X7,CPI,X8,CELZRO,X9,
     1     X10,X11,X12,X13,X14,X15,SIGMA,X16,DELTIM,DELT,TFREZ,
     2     RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN,TCW,TCICE,TCSAND,TCCLAY,
     3     TCOM,TCDRYS,RHOSOL,RHOOM,HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,
     4     HCPCLY,SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,TCGLAC,CLHMLT,
     5     CLHVAP,PI,ZOLNG,ZOLNS,ZOLNI,ZORATG,ALVSI,ALIRI,ALVSO,ALIRO,
     6     ALBRCK,DELTA,CGRAV,CKARM,CPD,AS,ASX,CI,BS,BETA,FACTN,HMIN,
     7     ANGMAX,A,B
C
c================= CTEM array declaration ===============================\
c
c     Local variables for coupling CLASS and CTEM
c
      integer ictemmod

      integer strlen !strlen can go in arg reader subroutine.
      character*80   titlec1 !, titlec2, titlec3
      character*80   argbuff
      character*160  command
c
!       logical ctem_on,    parallelrun,    mosaic,  !NOW POINTERS (BELOW)
!      1        cyclemet,   dofire,         run_model,
!      2     met_rewound,   reach_eof,      compete,
!      3      start_bare,   rsfile,         lnduseon,
!      4           co2on,   popdon,         inibioclim,
!      5    start_from_rs, dowetlands,      obswetf,
!      6    transient_run
c
       integer   lopcount,  isumc,
     1           k1c,       k2c,     jhhstd,
     2           jhhendd,   jdstd,   jdendd,      jhhsty,
     3           jhhendy,   jdsty,   jdendy,
     4           month1,
     5           month2,      xday,  ctemloop,nummetcylyrs,
     6           ncyear,  co2yr,
     5           spinfast,   nol2pfts(4),
     6           popyr,
     7           metcylyrst, metcycendyr, climiyear, popcycleyr,
     8           cypopyr, lucyr, cylucyr, endyr,bigpftc(2),
     9           obswetyr, cywetldyr, trans_startyr, jmosty
c
!        real      fsstar_g,
!      1           flstar_g,  qh_g,    qe_g,        snomlt_g,
!      2           beg_g,     gtout_g, tpn_g,       altot_g,
!      3           tcn_g,     tsn_g,   zsn_g

      real, pointer ::  fsstar_g
      real, pointer ::  flstar_g
      real, pointer ::  qh_g
      real, pointer ::  qe_g
      real, pointer ::  snomlt_g
      real, pointer ::  beg_g
      real, pointer ::  gtout_g
      real, pointer ::  tpn_g
      real, pointer ::  altot_g
      real, pointer ::  tcn_g
      real, pointer ::  tsn_g
      real, pointer ::  zsn_g


       real      co2concin,  popdin,    setco2conc, sumfare,
     1           temp_var, barefrac,  todfrac(ilg,icc), barf(nlat)

      real grclarea(ilg), crop_temp_frac(ilg,2)
c
!       real tcanrs(nlat,nmos), tsnors(nlat,nmos), tpndrs(nlat,nmos),
!      1     csum(nlat,nmos,ican),       tbaraccrow_m(nlat,nmos,ignd),
!      2     tcanoaccrow_m(nlat,nmos),   lightng(ilg),
!      3     uvaccrow_m(nlat,nmos),      vvaccrow_m(nlat,nmos)

c     Competition related variables

       real fsinacc_gat(ilg), flutacc_gat(ilg), flinacc_gat(ilg),
     1      alswacc_gat(ilg), allwacc_gat(ilg), pregacc_gat(ilg),
     2      altot_gat,        fsstar_gat,       flstar_gat,
     3      netrad_gat(ilg),  preacc_gat(ilg)
c
!       These go into CTEM but that is about it...
       real tcurm(ilg),       srpcuryr   (ilg), dftcuryr(ilg),
     1      tmonth(12,ilg),      anpcpcur(ilg),  anpecur(ilg),
     2      gdd5cur(ilg),        surmncur(ilg), defmncur(ilg),
     3      srplscur(ilg),       defctcur(ilg)

!        real, dimension(ilg) :: twarmm    ! temperature of the warmest month (c)
!        real, dimension(ilg) :: tcoldm    ! temperature of the coldest month (c)
!        real, dimension(ilg) :: gdd5      ! growing degree days above 5 c
!        real, dimension(ilg) :: aridity   ! aridity index, ratio of potential evaporation to precipitation
!        real, dimension(ilg) :: srplsmon  ! number of months in a year with surplus water i.e.
!                                                   !  precipitation more than potential evaporation
!        real, dimension(ilg) :: defctmon  ! number of months in a year with water deficit i.e.
!                                                   ! precipitation less than potential evaporation
!        real, dimension(ilg) :: anndefct  ! annual water deficit (mm)
!        real, dimension(ilg) :: annsrpls  ! annual water surplus (mm)
!        real, dimension(ilg) :: annpcp    ! annual precipitation (mm)
!        real, dimension(ilg) :: dry_season_length  ! length of dry season (months)
c
!       These go into CTEM but that is about it...
       real lyglfmasgat(ilg,icc),   geremortgat(ilg,icc),
     1      intrmortgat(ilg,icc),     lambdagat(ilg,icc),
         !2            rnded_pft(icc),
     3            ccgat(ilg,icc),         mmgat(ilg,icc) !,
!     4            temparray(icc),                  temp

      real  xdiffusgat(ilg) ! the corresponding ROW is CLASS's XDIFFUS

!     For these below, the corresponding ROWs are defined by CLASS

      real  sdepgat(ilg),       orgmgat(ilg,ignd),
     1      sandgat(ilg,ignd),  claygat(ilg,ignd)

!     Set up the variables that have both row and gat such that
!     they are beside each other in the declaration. This makes it
!     more apparent if the naming is inconsistent or if there is a problem
!     with how the variable is defined. (JM Jul 2013)

!     ||      This column is ROW       ||  This column is GAT    ||

!       logical pftexistrow(nlat,nmos,icc), pftexistgat(ilg,icc)
!
!       integer colddaysrow(nlat,nmos,2), colddaysgat(ilg,2),
!      1     icountrow(nlat,nmos),        icount(ilg),
!      2     lfstatusrow(nlat,nmos,icc),  lfstatusgat(ilg,icc),
!      3     pandaysrow(nlat,nmos,icc),   pandaysgat(ilg,icc),
!      4     stdalngrd(nlat),             stdalngat(ilg)
c
      !FLAG testing area!
!       real ailcminrow(nlat,nmos,icc),   ailcmingat(ilg,icc),
!      1     ailcmaxrow(nlat,nmos,icc),   ailcmaxgat(ilg,icc),
!      2     dvdfcanrow(nlat,nmos,icc),   dvdfcan(ilg,icc),
!      3     gleafmasrow(nlat,nmos,icc),  gleafmasgat(ilg,icc),
!      4     bleafmasrow(nlat,nmos,icc),  bleafmasgat(ilg,icc),
!      5     stemmassrow(nlat,nmos,icc),  stemmassgat(ilg,icc),
!      6     rootmassrow(nlat,nmos,icc),  rootmassgat(ilg,icc),
!      7     pstemmassrow(nlat,nmos,icc), pstemmassgat(ilg,icc),
!      8     pgleafmassrow(nlat,nmos,icc), pgleafmassgat(ilg,icc)
!
!       real fcancmxrow(nlat,nmos,icc),  fcancmxgat(ilg,icc),
!      1     gavglairow(nlat,nmos),      gavglaigat(ilg),
!      2     zolncrow(nlat,nmos,ican),   zolncgat(ilg,ican),
!      3     ailcrow(nlat,nmos,ican),    ailcgat(ilg,ican),
!      4     ailcgrow(nlat,nmos,icc),    ailcggat(ilg,icc),
!      5     ailcgsrow(nlat,nmos,icc),   ailcgsgat(ilg,icc),
!      6     fcancsrow(nlat,nmos,icc),   fcancsgat(ilg,icc),
!      7     fcancrow(nlat,nmos,icc),    fcancgat(ilg,icc),
!      8     co2concrow(nlat,nmos),      co2concgat(ilg),
!      9     co2i1cgrow(nlat,nmos,icc),  co2i1cggat(ilg,icc),
!      a     co2i1csrow(nlat,nmos,icc),  co2i1csgat(ilg,icc),
!      b     co2i2cgrow(nlat,nmos,icc),  co2i2cggat(ilg,icc),
!      c     co2i2csrow(nlat,nmos,icc),  co2i2csgat(ilg,icc),
!      d     ancsvegrow(nlat,nmos,icc),  ancsveggat(ilg,icc),
!      e     ancgvegrow(nlat,nmos,icc),  ancgveggat(ilg,icc),
!      f     rmlcsvegrow(nlat,nmos,icc), rmlcsveggat(ilg,icc),
!      g     rmlcgvegrow(nlat,nmos,icc), rmlcgveggat(ilg,icc),
!      h     slairow(nlat,nmos,icc),     slaigat(ilg,icc),
!      i     ailcbrow(nlat,nmos,icc),    ailcbgat(ilg,icc),
!      j     canresrow(nlat,nmos),       canresgat(ilg),
!      a     flhrlossrow(nlat,nmos,icc), flhrlossgat(ilg,icc)
! c
!       real grwtheffrow(nlat,nmos,icc), grwtheffgat(ilg,icc),
!      1     lystmmasrow(nlat,nmos,icc), lystmmasgat(ilg,icc),
!      2     lyrotmasrow(nlat,nmos,icc), lyrotmasgat(ilg,icc),
!      3     tymaxlairow(nlat,nmos,icc), tymaxlaigat(ilg,icc),
!      4     vgbiomasrow(nlat,nmos),     vgbiomasgat(ilg),
!      5     gavgltmsrow(nlat,nmos),     gavgltmsgat(ilg),
!      6     gavgscmsrow(nlat,nmos),     gavgscmsgat(ilg),
!      7     stmhrlosrow(nlat,nmos,icc), stmhrlosgat(ilg,icc)

!       real rmatcrow(nlat,nmos,ican,ignd),  rmatcgat(ilg,ican,ignd),
!      1     rmatctemrow(nlat,nmos,icc,ignd),rmatctemgat(ilg,icc,ignd),
!      2     litrmassrow(nlat,nmos,iccp1),   litrmassgat(ilg,iccp1),
!      3     soilcmasrow(nlat,nmos,iccp1),   soilcmasgat(ilg,iccp1),
!      4     vgbiomas_vegrow(nlat,nmos,icc), vgbiomas_veggat(ilg,icc)

! c     Fire-related variables
!
!       real emit_co2row(nlat,nmos,icc), emit_co2gat(ilg,icc),
!      1     emit_corow(nlat,nmos,icc),  emit_cogat(ilg,icc),
!      2     emit_ch4row(nlat,nmos,icc), emit_ch4gat(ilg,icc),
!      3     emit_nmhcrow(nlat,nmos,icc),emit_nmhcgat(ilg,icc),
!      4     emit_h2row(nlat,nmos,icc),  emit_h2gat(ilg,icc),
!      5     emit_noxrow(nlat,nmos,icc), emit_noxgat(ilg,icc),
!      6     emit_n2orow(nlat,nmos,icc), emit_n2ogat(ilg,icc),
!      7     emit_pm25row(nlat,nmos,icc),emit_pm25gat(ilg,icc),
!      8     emit_tpmrow(nlat,nmos,icc), emit_tpmgat(ilg,icc),
!      9     emit_tcrow(nlat,nmos,icc),  emit_tcgat(ilg,icc),
!      a     emit_ocrow(nlat,nmos,icc),  emit_ocgat(ilg,icc),
!      b     emit_bcrow(nlat,nmos,icc),  emit_bcgat(ilg,icc),
!      c     burnfracrow(nlat,nmos),     burnfracgat(ilg),
!      d     burnvegfrow(nlat,nmos,icc), burnvegfgat(ilg,icc),
!      e     probfirerow(nlat,nmos),     probfiregat(ilg),
!      f     btermrow(nlat,nmos),        btermgat(ilg),
!      g     ltermrow(nlat,nmos),        ltermgat(ilg),
!      h     mtermrow(nlat,nmos),        mtermgat(ilg)
!
!        real extnprobgrd(nlat),          extnprobgat(ilg),
!      1      prbfrhucgrd(nlat),          prbfrhucgat(ilg),
!      1      mlightnggrd(nlat,12),       mlightnggat(ilg,12)

! c      Methane(wetland) related variables    !Rudra added on 03/12/2013
!
!        real  WETFRACGRD(nlat),              wetfrac_sgrd(ilg,8),
!      1       CH4WET1ROW(nlat,nmos),         CH4WET1GAT(ILG),
!      2       CH4WET2ROW(nlat,nmos),         CH4WET2GAT(ILG),
!      3       WETFDYNROW(nlat,nmos),         WETFDYNGAT(ILG),
!      4       CH4DYN1ROW(nlat,nmos),         CH4DYN1GAT(ILG),
!      5       CH4DYN2ROW(nlat,nmos),         CH4DYN2GAT(ILG),
!      6       wetfrac_mon(nlat,12)
!
! !      Land-use related variables
!
!       real lucemcomrow(nlat,nmos),     lucemcomgat(ilg),
!      1     lucltrinrow(nlat,nmos),     lucltringat(ilg),
!      2     lucsocinrow(nlat,nmos),     lucsocingat(ilg)
! c
!       real bmasvegrow(nlat,nmos,icc),  bmasveggat(ilg,icc),
!      1     cmasvegcrow(nlat,nmos,ican),cmasvegcgat(ilg,ican),
!      2     veghghtrow(nlat,nmos,icc),  veghghtgat(ilg,icc),
!      3     rootdpthrow(nlat,nmos,icc), rootdpthgat(ilg,icc),
!      4     rmlrow(nlat,nmos),          rmlgat(ilg),
!      5     rmsrow(nlat,nmos),          rmsgat(ilg),
!      6     tltrleafrow(nlat,nmos,icc), tltrleafgat(ilg,icc),
!      7     tltrstemrow(nlat,nmos,icc), tltrstemgat(ilg,icc),
!      8     tltrrootrow(nlat,nmos,icc), tltrrootgat(ilg,icc),
!      9     leaflitrrow(nlat,nmos,icc), leaflitrgat(ilg,icc),
!      a     roottemprow(nlat,nmos,icc), roottempgat(ilg,icc),
!      b     afrleafrow(nlat,nmos,icc),  afrleafgat(ilg,icc),
!      c     afrstemrow(nlat,nmos,icc),  afrstemgat(ilg,icc),
!      d     afrrootrow(nlat,nmos,icc),  afrrootgat(ilg,icc),
!      e     wtstatusrow(nlat,nmos,icc), wtstatusgat(ilg,icc),
!      f     ltstatusrow(nlat,nmos,icc), ltstatusgat(ilg,icc),
!      g     rmrrow(nlat,nmos),          rmrgat(ilg)

!       real npprow(nlat,nmos),          nppgat(ilg),
!      1     neprow(nlat,nmos),          nepgat(ilg),
!      2     nbprow(nlat,nmos),          nbpgat(ilg),
!      3     gpprow(nlat,nmos),          gppgat(ilg),
!      4     hetroresrow(nlat,nmos),     hetroresgat(ilg),
!      5     autoresrow(nlat,nmos),      autoresgat(ilg),
!      6     soilcresprow(nlat,nmos),    soilcrespgat(ilg),
!      7     rmrow(nlat,nmos),           rmgat(ilg),
!      8     rgrow(nlat,nmos),           rggat(ilg),
!      9     litresrow(nlat,nmos),       litresgat(ilg),
!      a     socresrow(nlat,nmos),       socresgat(ilg),
!      b     dstcemlsrow(nlat,nmos),     dstcemlsgat(ilg),
!      c     litrfallrow(nlat,nmos),     litrfallgat(ilg),
!      d     humiftrsrow(nlat,nmos),     humiftrsgat(ilg)
!
!       real  gppvegrow(nlat,nmos,icc),   gppveggat(ilg,icc),
!      1      nepvegrow(nlat,nmos,iccp1),   nepveggat(ilg,iccp1),
!      2      nbpvegrow(nlat,nmos,iccp1),   nbpveggat(ilg,iccp1),
!      3      nppvegrow(nlat,nmos,icc),   nppveggat(ilg,icc),
!      4      hetroresvegrow(nlat,nmos,iccp1),hetroresveggat(ilg,iccp1),
!      5      autoresvegrow(nlat,nmos,icc),autoresveggat(ilg,icc),
!      6      litresvegrow(nlat,nmos,iccp1),litresveggat(ilg,iccp1),
!      7      soilcresvegrow(nlat,nmos,iccp1),soilcresveggat(ilg,iccp1),
!      8      rmlvegaccrow(nlat,nmos,icc),rmlvegaccgat(ilg,icc),
!      9      rmsvegrow(nlat,nmos,icc),   rmsveggat(ilg,icc),
!      a      rmrvegrow(nlat,nmos,icc),   rmrveggat(ilg,icc),
!      b      rgvegrow(nlat,nmos,icc),    rgveggat(ilg,icc)
! c
!       real
!      1     rothrlosrow(nlat,nmos,icc), rothrlosgat(ilg,icc),
!      2     pfcancmxrow(nlat,nmos,icc), pfcancmxgat(ilg,icc),
!      3     nfcancmxrow(nlat,nmos,icc), nfcancmxgat(ilg,icc),
!      4     alvsctmrow(nlat,nmos,ican), alvsctmgat(ilg,ican),
!      5     paicrow(nlat,nmos,ican),    paicgat(ilg,ican),
!      6     slaicrow(nlat,nmos,ican),   slaicgat(ilg,ican),
!      7     alirctmrow(nlat,nmos,ican), alirctmgat(ilg,ican),
!      8     cfluxcgrow(nlat,nmos),      cfluxcggat(ilg),
!      9     cfluxcsrow(nlat,nmos),      cfluxcsgat(ilg),
!      a     dstcemls3row(nlat,nmos),    dstcemls3gat(ilg),
!      b     anvegrow(nlat,nmos,icc),    anveggat(ilg,icc),
!      c     rmlvegrow(nlat,nmos,icc),   rmlveggat(ilg,icc)

      ! CLASS Monthly Outputs:

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
      real, pointer, dimension(:) :: TAACC_MO

      real, pointer :: FSSTAR_MO
      real, pointer :: FLSTAR_MO
      real, pointer :: QH_MO
      real, pointer :: QE_MO

      real, pointer, dimension(:,:) :: TBARACC_MO
      real, pointer, dimension(:,:) :: THLQACC_MO
      real, pointer, dimension(:,:) :: THICACC_MO

    ! CLASS yearly output for class grid-mean

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
      real, pointer, dimension(:) :: TAACC_YR

      real, pointer :: FSSTAR_YR
      real, pointer :: FLSTAR_YR
      real, pointer :: QH_YR
      real, pointer :: QE_YR

      logical, pointer :: ctem_on
      logical, pointer :: parallelrun
      logical, pointer :: mosaic
      logical, pointer :: cyclemet
      logical, pointer :: dofire
      logical, pointer :: run_model
      logical, pointer :: met_rewound
      logical, pointer :: reach_eof
      logical, pointer :: compete
      logical, pointer :: start_bare
      logical, pointer :: rsfile
      logical, pointer :: lnduseon
      logical, pointer :: co2on
      logical, pointer :: popdon
      logical, pointer :: inibioclim
      logical, pointer :: start_from_rs
      logical, pointer :: dowetlands
      logical, pointer :: obswetf
      logical, pointer :: transient_run

      ! ROW vars:
      logical, pointer, dimension(:,:,:) :: pftexistrow
      integer, pointer, dimension(:,:,:) :: colddaysrow
      integer, pointer, dimension(:,:) :: icountrow
      integer, pointer, dimension(:,:,:) :: lfstatusrow
      integer, pointer, dimension(:,:,:) :: pandaysrow

      integer, pointer, dimension(:) :: stdalngrd

      real, pointer, dimension(:,:) :: tcanrs
      real, pointer, dimension(:,:) :: tsnors
      real, pointer, dimension(:,:) :: tpndrs
      real, pointer, dimension(:,:,:) :: csum
      real, pointer, dimension(:,:,:) :: tbaraccrow_m
      real, pointer, dimension(:,:) :: tcanoaccrow_m
      real, pointer, dimension(:,:) :: uvaccrow_m
      real, pointer, dimension(:,:) :: vvaccrow_m

      real, pointer, dimension(:,:,:) :: ailcminrow         !
      real, pointer, dimension(:,:,:) :: ailcmaxrow         !
      real, pointer, dimension(:,:,:) :: dvdfcanrow         !
      real, pointer, dimension(:,:,:) :: gleafmasrow        !
      real, pointer, dimension(:,:,:) :: bleafmasrow        !
      real, pointer, dimension(:,:,:) :: stemmassrow        !
      real, pointer, dimension(:,:,:) :: rootmassrow        !
      real, pointer, dimension(:,:,:) :: pstemmassrow       !
      real, pointer, dimension(:,:,:) :: pgleafmassrow      !
      real, pointer, dimension(:,:,:) :: fcancmxrow
      real, pointer, dimension(:,:) :: gavglairow
      real, pointer, dimension(:,:,:) :: zolncrow
      real, pointer, dimension(:,:,:) :: ailcrow
      real, pointer, dimension(:,:,:) :: ailcgrow
      real, pointer, dimension(:,:,:) :: ailcgsrow
      real, pointer, dimension(:,:,:) :: fcancsrow
      real, pointer, dimension(:,:,:) :: fcancrow
      real, pointer, dimension(:,:) :: co2concrow
      real, pointer, dimension(:,:,:) :: co2i1cgrow
      real, pointer, dimension(:,:,:) :: co2i1csrow
      real, pointer, dimension(:,:,:) :: co2i2cgrow
      real, pointer, dimension(:,:,:) :: co2i2csrow
      real, pointer, dimension(:,:,:) :: ancsvegrow
      real, pointer, dimension(:,:,:) :: ancgvegrow
      real, pointer, dimension(:,:,:) :: rmlcsvegrow
      real, pointer, dimension(:,:,:) :: rmlcgvegrow
      real, pointer, dimension(:,:,:) :: slairow
      real, pointer, dimension(:,:,:) :: ailcbrow
      real, pointer, dimension(:,:) :: canresrow
      real, pointer, dimension(:,:,:) :: flhrlossrow

      real, pointer, dimension(:,:,:) :: grwtheffrow
      real, pointer, dimension(:,:,:) :: lystmmasrow
      real, pointer, dimension(:,:,:) :: lyrotmasrow
      real, pointer, dimension(:,:,:) :: tymaxlairow
      real, pointer, dimension(:,:) :: vgbiomasrow
      real, pointer, dimension(:,:) :: gavgltmsrow
      real, pointer, dimension(:,:) :: gavgscmsrow
      real, pointer, dimension(:,:,:) :: stmhrlosrow
      real, pointer, dimension(:,:,:,:) :: rmatcrow
      real, pointer, dimension(:,:,:,:) :: rmatctemrow
      real, pointer, dimension(:,:,:) :: litrmassrow
      real, pointer, dimension(:,:,:) :: soilcmasrow
      real, pointer, dimension(:,:,:) :: vgbiomas_vegrow

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
      real, pointer, dimension(:,:) :: probfirerow
      real, pointer, dimension(:,:) :: btermrow
      real, pointer, dimension(:,:) :: ltermrow
      real, pointer, dimension(:,:) :: mtermrow

      real, pointer, dimension(:) :: extnprobgrd
      real, pointer, dimension(:) :: prbfrhucgrd
      real, pointer, dimension(:,:) :: mlightnggrd

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

      real, pointer, dimension(:) :: wetfracgrd
      real, pointer, dimension(:,:) :: ch4wet1row
      real, pointer, dimension(:,:) :: ch4wet2row
      real, pointer, dimension(:,:) :: wetfdynrow
      real, pointer, dimension(:,:) :: ch4dyn1row
      real, pointer, dimension(:,:) :: ch4dyn2row
      real, pointer, dimension(:,:) :: wetfrac_mon

      real, pointer, dimension(:,:) :: lucemcomrow
      real, pointer, dimension(:,:) :: lucltrinrow
      real, pointer, dimension(:,:) :: lucsocinrow

      real, pointer, dimension(:,:) :: npprow
      real, pointer, dimension(:,:) :: neprow
      real, pointer, dimension(:,:) :: nbprow
      real, pointer, dimension(:,:) :: gpprow
      real, pointer, dimension(:,:) :: hetroresrow
      real, pointer, dimension(:,:) :: autoresrow
      real, pointer, dimension(:,:) :: soilcresprow
      real, pointer, dimension(:,:) :: rmrow
      real, pointer, dimension(:,:) :: rgrow
      real, pointer, dimension(:,:) :: litresrow
      real, pointer, dimension(:,:) :: socresrow
      real, pointer, dimension(:,:) :: dstcemlsrow
      real, pointer, dimension(:,:) :: litrfallrow
      real, pointer, dimension(:,:) :: humiftrsrow

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

      real, pointer, dimension(:,:,:) :: rothrlosrow
      real, pointer, dimension(:,:,:) :: pfcancmxrow
      real, pointer, dimension(:,:,:) :: nfcancmxrow
      real, pointer, dimension(:,:,:) :: alvsctmrow
      real, pointer, dimension(:,:,:) :: paicrow
      real, pointer, dimension(:,:,:) :: slaicrow
      real, pointer, dimension(:,:,:) :: alirctmrow
      real, pointer, dimension(:,:) :: cfluxcgrow
      real, pointer, dimension(:,:) :: cfluxcsrow
      real, pointer, dimension(:,:) :: dstcemls3row
      real, pointer, dimension(:,:,:) :: anvegrow
      real, pointer, dimension(:,:,:) :: rmlvegrow

      ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.
      ! GAT version:

      logical, pointer, dimension(:,:) :: pftexistgat
      integer, pointer, dimension(:,:) :: colddaysgat
      integer, pointer, dimension(:) :: icountgat
      integer, pointer, dimension(:,:) :: lfstatusgat
      integer, pointer, dimension(:,:) :: pandaysgat

      integer, pointer, dimension(:) :: stdalngat

      real, pointer, dimension(:) :: lightng

      real, pointer, dimension(:,:) :: ailcmingat         !
      real, pointer, dimension(:,:) :: ailcmaxgat         !
      real, pointer, dimension(:,:) :: dvdfcangat         !
      real, pointer, dimension(:,:) :: gleafmasgat        !
      real, pointer, dimension(:,:) :: bleafmasgat        !
      real, pointer, dimension(:,:) :: stemmassgat        !
      real, pointer, dimension(:,:) :: rootmassgat        !
      real, pointer, dimension(:,:) :: pstemmassgat       !
      real, pointer, dimension(:,:) :: pgleafmassgat      !
      real, pointer, dimension(:,:) :: fcancmxgat
      real, pointer, dimension(:) :: gavglaigat
      real, pointer, dimension(:,:) :: zolncgat
      real, pointer, dimension(:,:) :: ailcgat
      real, pointer, dimension(:,:) :: ailcggat
      real, pointer, dimension(:,:) :: ailcgsgat
      real, pointer, dimension(:,:) :: fcancsgat
      real, pointer, dimension(:,:) :: fcancgat
      real, pointer, dimension(:) :: co2concgat
      real, pointer, dimension(:,:) :: co2i1cggat
      real, pointer, dimension(:,:) :: co2i1csgat
      real, pointer, dimension(:,:) :: co2i2cggat
      real, pointer, dimension(:,:) :: co2i2csgat
      real, pointer, dimension(:,:) :: ancsveggat
      real, pointer, dimension(:,:) :: ancgveggat
      real, pointer, dimension(:,:) :: rmlcsveggat
      real, pointer, dimension(:,:) :: rmlcgveggat
      real, pointer, dimension(:,:) :: slaigat
      real, pointer, dimension(:,:) :: ailcbgat
      real, pointer, dimension(:) :: canresgat
      real, pointer, dimension(:,:) :: flhrlossgat

      real, pointer, dimension(:,:) :: grwtheffgat
      real, pointer, dimension(:,:) :: lystmmasgat
      real, pointer, dimension(:,:) :: lyrotmasgat
      real, pointer, dimension(:,:) :: tymaxlaigat
      real, pointer, dimension(:) :: vgbiomasgat
      real, pointer, dimension(:) :: gavgltmsgat
      real, pointer, dimension(:) :: gavgscmsgat
      real, pointer, dimension(:,:) :: stmhrlosgat
      real, pointer, dimension(:,:,:) :: rmatcgat
      real, pointer, dimension(:,:,:) :: rmatctemgat
      real, pointer, dimension(:,:) :: litrmassgat
      real, pointer, dimension(:,:) :: soilcmasgat
      real, pointer, dimension(:,:) :: vgbiomas_veggat

      real, pointer, dimension(:,:) :: emit_co2gat
      real, pointer, dimension(:,:) :: emit_cogat
      real, pointer, dimension(:,:) :: emit_ch4gat
      real, pointer, dimension(:,:) :: emit_nmhcgat
      real, pointer, dimension(:,:) :: emit_h2gat
      real, pointer, dimension(:,:) :: emit_noxgat
      real, pointer, dimension(:,:) :: emit_n2ogat
      real, pointer, dimension(:,:) :: emit_pm25gat
      real, pointer, dimension(:,:) :: emit_tpmgat
      real, pointer, dimension(:,:) :: emit_tcgat
      real, pointer, dimension(:,:) :: emit_ocgat
      real, pointer, dimension(:,:) :: emit_bcgat
      real, pointer, dimension(:) :: burnfracgat
      real, pointer, dimension(:,:) :: burnvegfgat
      real, pointer, dimension(:) :: probfiregat
      real, pointer, dimension(:) :: btermgat
      real, pointer, dimension(:) :: ltermgat
      real, pointer, dimension(:) :: mtermgat

      real, pointer, dimension(:) :: extnprobgat
      real, pointer, dimension(:) :: prbfrhucgat
      real, pointer, dimension(:,:) :: mlightnggat

      real, pointer, dimension(:,:) :: bmasveggat
      real, pointer, dimension(:,:) :: cmasvegcgat
      real, pointer, dimension(:,:) :: veghghtgat
      real, pointer, dimension(:,:) :: rootdpthgat
      real, pointer, dimension(:) :: rmlgat
      real, pointer, dimension(:) :: rmsgat
      real, pointer, dimension(:,:) :: tltrleafgat
      real, pointer, dimension(:,:) :: tltrstemgat
      real, pointer, dimension(:,:) :: tltrrootgat
      real, pointer, dimension(:,:) :: leaflitrgat
      real, pointer, dimension(:,:) :: roottempgat
      real, pointer, dimension(:,:) :: afrleafgat
      real, pointer, dimension(:,:) :: afrstemgat
      real, pointer, dimension(:,:) :: afrrootgat
      real, pointer, dimension(:,:) :: wtstatusgat
      real, pointer, dimension(:,:) :: ltstatusgat
      real, pointer, dimension(:) :: rmrgat

      real, pointer, dimension(:,:) :: wetfrac_sgrd
      real, pointer, dimension(:) :: ch4wet1gat
      real, pointer, dimension(:) :: ch4wet2gat
      real, pointer, dimension(:) :: wetfdyngat
      real, pointer, dimension(:) :: ch4dyn1gat
      real, pointer, dimension(:) :: ch4dyn2gat

      real, pointer, dimension(:) :: lucemcomgat
      real, pointer, dimension(:) :: lucltringat
      real, pointer, dimension(:) :: lucsocingat

      real, pointer, dimension(:) :: nppgat
      real, pointer, dimension(:) :: nepgat
      real, pointer, dimension(:) :: nbpgat
      real, pointer, dimension(:) :: gppgat
      real, pointer, dimension(:) :: hetroresgat
      real, pointer, dimension(:) :: autoresgat
      real, pointer, dimension(:) :: soilcrespgat
      real, pointer, dimension(:) :: rmgat
      real, pointer, dimension(:) :: rggat
      real, pointer, dimension(:) :: litresgat
      real, pointer, dimension(:) :: socresgat
      real, pointer, dimension(:) :: dstcemlsgat
      real, pointer, dimension(:) :: litrfallgat
      real, pointer, dimension(:) :: humiftrsgat

      real, pointer, dimension(:,:) :: gppveggat
      real, pointer, dimension(:,:) :: nepveggat
      real, pointer, dimension(:,:) :: nbpveggat
      real, pointer, dimension(:,:) :: nppveggat
      real, pointer, dimension(:,:) :: hetroresveggat
      real, pointer, dimension(:,:) :: autoresveggat
      real, pointer, dimension(:,:) :: litresveggat
      real, pointer, dimension(:,:) :: soilcresveggat
      real, pointer, dimension(:,:) :: rmlvegaccgat
      real, pointer, dimension(:,:) :: rmsveggat
      real, pointer, dimension(:,:) :: rmrveggat
      real, pointer, dimension(:,:) :: rgveggat

      real, pointer, dimension(:,:) :: rothrlosgat
      real, pointer, dimension(:,:) :: pfcancmxgat
      real, pointer, dimension(:,:) :: nfcancmxgat
      real, pointer, dimension(:,:) :: alvsctmgat
      real, pointer, dimension(:,:) :: paicgat
      real, pointer, dimension(:,:) :: slaicgat
      real, pointer, dimension(:,:) :: alirctmgat
      real, pointer, dimension(:) :: cfluxcggat
      real, pointer, dimension(:) :: cfluxcsgat
      real, pointer, dimension(:) :: dstcemls3gat
      real, pointer, dimension(:,:) :: anveggat
      real, pointer, dimension(:,:) :: rmlveggat

      real, pointer, dimension(:) :: twarmm    ! temperature of the warmest month (c)
      real, pointer, dimension(:) :: tcoldm    ! temperature of the coldest month (c)
      real, pointer, dimension(:) :: gdd5      ! growing degree days above 5 c
      real, pointer, dimension(:) :: aridity   ! aridity index, ratio of potential evaporation to precipitation
      real, pointer, dimension(:) :: srplsmon  ! number of months in a year with surplus water i.e.
                                                  !  precipitation more than potential evaporation
      real, pointer, dimension(:) :: defctmon  ! number of months in a year with water deficit i.e.
                                                  ! precipitation less than potential evaporation
      real, pointer, dimension(:) :: anndefct  ! annual water deficit (mm)
      real, pointer, dimension(:) :: annsrpls  ! annual water surplus (mm)
      real, pointer, dimension(:) :: annpcp    ! annual precipitation (mm)
      real, pointer, dimension(:) :: dry_season_length  ! length of dry season (months)

      ! Mosaic level:

      real, pointer, dimension(:,:) :: PREACC_M
      real, pointer, dimension(:,:) :: GTACC_M
      real, pointer, dimension(:,:) :: QEVPACC_M
      real, pointer, dimension(:,:) :: HFSACC_M
      real, pointer, dimension(:,:) :: HMFNACC_M
      real, pointer, dimension(:,:) :: ROFACC_M
      real, pointer, dimension(:,:) :: SNOACC_M
      real, pointer, dimension(:,:) :: OVRACC_M
      real, pointer, dimension(:,:) :: WTBLACC_M
      real, pointer, dimension(:,:,:) :: TBARACC_M
      real, pointer, dimension(:,:,:) :: THLQACC_M
      real, pointer, dimension(:,:,:) :: THICACC_M
      real, pointer, dimension(:,:,:) :: THALACC_M
      real, pointer, dimension(:,:) :: ALVSACC_M
      real, pointer, dimension(:,:) :: ALIRACC_M
      real, pointer, dimension(:,:) :: RHOSACC_M
      real, pointer, dimension(:,:) :: TSNOACC_M
      real, pointer, dimension(:,:) :: WSNOACC_M
      real, pointer, dimension(:,:) :: SNOARE_M
      real, pointer, dimension(:,:) :: TCANACC_M
      real, pointer, dimension(:,:) :: RCANACC_M
      real, pointer, dimension(:,:) :: SCANACC_M
      real, pointer, dimension(:,:) :: GROACC_M
      real, pointer, dimension(:,:) :: FSINACC_M
      real, pointer, dimension(:,:) :: FLINACC_M
      real, pointer, dimension(:,:) :: TAACC_M
      real, pointer, dimension(:,:) :: UVACC_M
      real, pointer, dimension(:,:) :: PRESACC_M
      real, pointer, dimension(:,:) :: QAACC_M
      real, pointer, dimension(:,:) :: EVAPACC_M
      real, pointer, dimension(:,:) :: FLUTACC_M

!      Outputs

!        real tcanoaccrow_out(nlat,nmos),tcanoaccgat_out(ilg),
!      1     qevpacc_m_save(nlat,nmos)
       real, pointer, dimension(:,:) :: tcanoaccrow_out
       real, pointer, dimension(:) :: tcanoaccgat_out
       real, pointer, dimension(:,:) :: qevpacc_m_save
c
!       integer ifcancmx_g(nlat,icc)  !FLAG!! What is the point of this one? JM Apr 1 2015.

!     -----------------------
!      Mosaic-level variables (denoted by an ending of "_m")
c
!        real   PREACC_M(NLAT,NMOS),     GTACC_M(NLAT,NMOS),
!      1     QEVPACC_M(NLAT,NMOS),       HFSACC_M(NLAT,NMOS),
!      2     HMFNACC_M(NLAT,NMOS),       ROFACC_M(NLAT,NMOS),
!      3     SNOACC_M(NLAT,NMOS),        OVRACC_M(NLAT,NMOS),
!      3     WTBLACC_M(NLAT,NMOS),       TBARACC_M(NLAT,NMOS,IGND),
!      4     THLQACC_M(NLAT,NMOS,IGND),  THICACC_M(NLAT,NMOS,IGND),
!      5     THALACC_M(NLAT,NMOS,IGND),  ALVSACC_M(NLAT,NMOS),
!      6     ALIRACC_M(NLAT,NMOS),       RHOSACC_M(NLAT,NMOS),
!      7     TSNOACC_M(NLAT,NMOS),       WSNOACC_M(NLAT,NMOS),
!      8     TCANACC_M(NLAT,NMOS),       RCANACC_M(NLAT,NMOS),
!      9     SCANACC_M(NLAT,NMOS),       GROACC_M(NLAT,NMOS),
!      A     FSINACC_M(NLAT,NMOS),       FLINACC_M(NLAT,NMOS),
!      B     TAACC_M(NLAT,NMOS),         UVACC_M(NLAT,NMOS),
!      C     PRESACC_M(NLAT,NMOS),       QAACC_M(NLAT,NMOS),
!      D     EVAPACC_M(NLAT,NMOS),       FLUTACC_M(NLAT,NMOS)

      real faregat(ilg)


!       real  fsnowacc_m(ilg),           tcansacc_m(ilg),
!      1      tbarcacc_m(ilg,ignd),      tbarcsacc_m(ilg,ignd),
!      2      tbargacc_m(ilg,ignd),      tbargsacc_m(ilg,ignd),
!      3      thliqcacc_m(ilg,ignd),     thliqgacc_m(ilg,ignd),
!      4      thicecacc_m(ilg,ignd),
!      5      tcanoaccgat_m(ilg),        taaccgat_m(ilg),
!      6      uvaccgat_m(ilg),           vvaccgat_m(ilg),
!      7      tbaraccgat_m(ilg,ignd)
!
!       real ancsvgac_m(ilg,icc),        ancgvgac_m(ilg,icc),
!      b     rmlcsvga_m(ilg,icc),        rmlcgvga_m(ilg,icc)
!
!        integer ifcancmx_m(nlat,nmos)
!
!        real leaflitr_m(nlat,nmos),     tltrleaf_m(nlat,nmos),
!      1     tltrstem_m(nlat,nmos),      tltrroot_m(nlat,nmos),
!      2     ailcg_m(nlat,nmos),         ailcb_m(nlat,nmos),
!      3     rmatctem_m(nlat,nmos,ignd), veghght_m(nlat,nmos),
!      4     rootdpth_m(nlat,nmos),      roottemp_m(nlat,nmos),
!      5     slai_m(nlat,nmos),          afrroot_m(nlat,nmos),
!      6     afrleaf_m(nlat,nmos),       afrstem_m(nlat,nmos),
!      7     laimaxg_m(nlat,nmos),       stemmass_m(nlat,nmos),
!      8     rootmass_m(nlat,nmos),      litrmass_m(nlat,nmos),
!      9     gleafmas_m(nlat,nmos),      bleafmas_m(nlat,nmos),
!      a     soilcmas_m(nlat,nmos)
c
      real, pointer, dimension(:,:) :: leaflitr_m
      real, pointer, dimension(:,:) :: tltrleaf_m
      real, pointer, dimension(:,:) :: tltrstem_m
      real, pointer, dimension(:,:) :: tltrroot_m
      real, pointer, dimension(:,:) :: ailcg_m
      real, pointer, dimension(:,:) :: ailcb_m
      real, pointer, dimension(:,:,:) :: rmatctem_m
      real, pointer, dimension(:,:) :: veghght_m
      real, pointer, dimension(:,:) :: rootdpth_m
      real, pointer, dimension(:,:) :: roottemp_m
      real, pointer, dimension(:,:) :: slai_m
      real, pointer, dimension(:,:) :: afrroot_m
      real, pointer, dimension(:,:) :: afrleaf_m
      real, pointer, dimension(:,:) :: afrstem_m
      real, pointer, dimension(:,:) :: laimaxg_m
      real, pointer, dimension(:,:) :: stemmass_m
      real, pointer, dimension(:,:) :: rootmass_m
      real, pointer, dimension(:,:) :: litrmass_m
      real, pointer, dimension(:,:) :: gleafmas_m
      real, pointer, dimension(:,:) :: bleafmas_m
      real, pointer, dimension(:,:) :: soilcmas_m

      real, pointer, dimension(:) :: fsnowacc_m
      real, pointer, dimension(:) :: tcansacc_m
      real, pointer, dimension(:) :: tcanoaccgat_m
      real, pointer, dimension(:) :: taaccgat_m
      real, pointer, dimension(:) :: uvaccgat_m
      real, pointer, dimension(:) :: vvaccgat_m
      real, pointer, dimension(:,:) :: tbaraccgat_m
      real, pointer, dimension(:,:) :: tbarcacc_m
      real, pointer, dimension(:,:) :: tbarcsacc_m
      real, pointer, dimension(:,:) :: tbargacc_m
      real, pointer, dimension(:,:) :: tbargsacc_m
      real, pointer, dimension(:,:) :: thliqcacc_m
      real, pointer, dimension(:,:) :: thliqgacc_m
      real, pointer, dimension(:,:) :: thicecacc_m
      real, pointer, dimension(:,:) :: ancsvgac_m
      real, pointer, dimension(:,:) :: ancgvgac_m
      real, pointer, dimension(:,:) :: rmlcsvga_m
      real, pointer, dimension(:,:) :: rmlcgvga_m
      integer, pointer, dimension(:,:) :: ifcancmx_m

!     -----------------------
!     Grid-averaged variables (denoted with an ending of "_g")

!       real  WSNOROT_g(nlat),           ROFSROT_g(nlat),
!      1     SNOROT_g(nlat),             RHOSROT_g(nlat),
!      2     ROFROT_g(nlat),             ZPNDROT_g(nlat),
!      3     RCANROT_g(nlat),            SCANROT_g(nlat),
!      4     TROFROT_g(nlat),            TROOROT_g(nlat),
!      5     TROBROT_g(nlat),            ROFOROT_g(nlat),
!      6     ROFBROT_g(nlat),            TROSROT_g(nlat),
!      7     FSGVROT_g(nlat),            FSGSROT_g(nlat),
!      8     FLGVROT_g(nlat),            FLGSROT_g(nlat),
!      9     HFSCROT_g(nlat),            HFSSROT_g(nlat),
!      a     HEVCROT_g(nlat),            HEVSROT_g(nlat),
!      b     HMFCROT_g(nlat),            HMFNROT_g(nlat),
!      c     HTCSROT_g(nlat),            HTCCROT_g(nlat),
!      d     FSGGROT_g(nlat),            FLGGROT_g(nlat),
!      e     HFSGROT_g(nlat),            HEVGROT_g(nlat),
!      f     fc_g(nlat),                 fg_g(nlat),
!      g     fcs_g(nlat),                fgs_g(nlat),
!      h     PCFCROT_g(nlat),            PCLCROT_g(nlat),
!      1     PCPGROT_g(nlat),            QFCFROT_g(nlat),
!      2     QFGROT_g(nlat),             QFCROT_g(nlat,ignd),
!      3     ROFCROT_g(nlat),            ROFNROT_g(nlat),
!      4     WTRSROT_g(nlat),            WTRGROT_g(nlat),
!      5     PCPNROT_g(nlat),            QFCLROT_g(nlat),
!      6     QFNROT_g(nlat),             WTRCROT_g(nlat),
!      8     gpp_g(nlat),                npp_g(nlat),
!      9     nbp_g(nlat),                autores_g(nlat),
!      a     litres_g(nlat),             socres_g(nlat),
!      b     dstcemls3_g(nlat),          litrfall_g(nlat),
!      c     rml_g(nlat),                rms_g(nlat),
!      d     rg_g(nlat),                 leaflitr_g(nlat),
!      e     tltrstem_g(nlat),           tltrroot_g(nlat),
!      f     nep_g(nlat),                hetrores_g(nlat),
!      g     dstcemls_g(nlat),           humiftrs_g(nlat),
!      h     rmr_g(nlat),                tltrleaf_g(nlat),
!      i     gavgltms_g(nlat)

      real, pointer, dimension(:) :: WSNOROT_g
      real, pointer, dimension(:) :: ROFSROT_g
      real, pointer, dimension(:) :: SNOROT_g
      real, pointer, dimension(:) :: RHOSROT_g
      real, pointer, dimension(:) :: ROFROT_g
      real, pointer, dimension(:) :: ZPNDROT_g
      real, pointer, dimension(:) :: RCANROT_g
      real, pointer, dimension(:) :: SCANROT_g
      real, pointer, dimension(:) :: TROFROT_g
      real, pointer, dimension(:) :: TROOROT_g
      real, pointer, dimension(:) :: TROBROT_g
      real, pointer, dimension(:) :: ROFOROT_g
      real, pointer, dimension(:) :: ROFBROT_g
      real, pointer, dimension(:) :: TROSROT_g
      real, pointer, dimension(:) :: FSGVROT_g
      real, pointer, dimension(:) :: FSGSROT_g
      real, pointer, dimension(:) :: FLGVROT_g
      real, pointer, dimension(:) :: FLGSROT_g
      real, pointer, dimension(:) :: HFSCROT_g
      real, pointer, dimension(:) :: HFSSROT_g
      real, pointer, dimension(:) :: HEVCROT_g
      real, pointer, dimension(:) :: HEVSROT_g
      real, pointer, dimension(:) :: HMFCROT_g
      real, pointer, dimension(:) :: HMFNROT_g
      real, pointer, dimension(:) :: HTCSROT_g
      real, pointer, dimension(:) :: HTCCROT_g
      real, pointer, dimension(:) :: FSGGROT_g
      real, pointer, dimension(:) :: FLGGROT_g
      real, pointer, dimension(:) :: HFSGROT_g
      real, pointer, dimension(:) :: HEVGROT_g
      real, pointer, dimension(:) :: CDHROT_g
      real, pointer, dimension(:) :: CDMROT_g
      real, pointer, dimension(:) :: SFCUROT_g
      real, pointer, dimension(:) :: SFCVROT_g
      real, pointer, dimension(:) :: fc_g
      real, pointer, dimension(:) :: fg_g
      real, pointer, dimension(:) :: fcs_g
      real, pointer, dimension(:) :: fgs_g
      real, pointer, dimension(:) :: PCFCROT_g
      real, pointer, dimension(:) :: PCLCROT_g
      real, pointer, dimension(:) :: PCPGROT_g
      real, pointer, dimension(:) :: QFCFROT_g
      real, pointer, dimension(:) :: QFGROT_g
      real, pointer, dimension(:,:) :: QFCROT_g
      real, pointer, dimension(:) :: ROFCROT_g
      real, pointer, dimension(:) :: ROFNROT_g
      real, pointer, dimension(:) :: WTRSROT_g
      real, pointer, dimension(:) :: WTRGROT_g
      real, pointer, dimension(:) :: PCPNROT_g
      real, pointer, dimension(:) :: QFCLROT_g
      real, pointer, dimension(:) :: QFNROT_g
      real, pointer, dimension(:) :: WTRCROT_g
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

! c
!        real  HMFGROT_g(nlat,ignd),     HTCROT_g(nlat,ignd),
!      1     TBARROT_g(nlat,ignd),       THLQROT_g(nlat,ignd),
!      2     THICROT_g(nlat,ignd),       GFLXROT_g(nlat,ignd),
!      3     anvegrow_g(nlat,icc),       rmlvegrow_g(nlat,icc)
! c
!        real vgbiomas_g(nlat),          gavglai_g(nlat),
!      1     gavgscms_g(nlat),           gleafmas_g(nlat),
!      2     bleafmas_g(nlat),           stemmass_g(nlat),
!      3     rootmass_g(nlat),           litrmass_g(nlat),
!      4     soilcmas_g(nlat),           slai_g(nlat),
!      5     ailcg_g(nlat),              ailcb_g(nlat),
!      6     rmatctem_g(nlat,ignd),      veghght_g(nlat),
!      7     rootdpth_g(nlat),           roottemp_g(nlat),
!      a     totcmass_g(nlat)
! c
!        real afrleaf_g(nlat,icc),       afrstem_g(nlat,icc),
!      1     afrroot_g(nlat,icc),        lfstatus_g(nlat,icc),
!      2     tcanoacc_out_g(nlat),
!      3     burnfrac_g(nlat),           probfire_g(nlat),
!      4     lucemcom_g(nlat),           lucltrin_g(nlat),
!      5     lucsocin_g(nlat),
!      6     emit_co2_g(nlat),           emit_co_g(nlat),
!      7     emit_ch4_g(nlat),           emit_nmhc_g(nlat),
!      8     emit_h2_g(nlat),            emit_nox_g(nlat),
!      9     emit_n2o_g(nlat),           emit_pm25_g(nlat),
!      a     emit_tpm_g(nlat),           emit_tc_g(nlat),
!      b     emit_oc_g(nlat),            emit_bc_g(nlat),
!      c     bterm_g(nlat),              lterm_g(nlat),
!      d     mterm_g(nlat)
!
!        real    CH4WET1_G(nlat),            CH4WET2_G(nlat),        !Rudra addes on 03/12/2013
!      &     WETFDYN_G(nlat),            CH4DYN1_G(nlat),
!      &     CH4DYN2_G(nlat)

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
      real, pointer, dimension(:) :: probfire_g
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
      real, pointer, dimension(:,:) :: afrleaf_g
      real, pointer, dimension(:,:) :: afrstem_g
      real, pointer, dimension(:,:) :: afrroot_g
      real, pointer, dimension(:,:) :: lfstatus_g
      real, pointer, dimension(:,:) :: rmlvegrow_g
      real, pointer, dimension(:,:) :: anvegrow_g
      real, pointer, dimension(:,:) :: rmatctem_g
      real, pointer, dimension(:,:) :: HMFGROT_g
      real, pointer, dimension(:,:) :: HTCROT_g
      real, pointer, dimension(:,:) :: TBARROT_g
      real, pointer, dimension(:,:) :: THLQROT_g
      real, pointer, dimension(:,:) :: THICROT_g
      real, pointer, dimension(:,:) :: GFLXROT_g

!     -----------------------
!      Grid averaged monthly variables (denoted by name ending in "_mo_g")

!        real laimaxg_mo_g(nlat),          stemmass_mo_g(nlat),
!      1      rootmass_mo_g(nlat),         litrmass_mo_g(nlat),
!      2      soilcmas_mo_g(nlat),         npp_mo_g(nlat),
!      3      gpp_mo_g(nlat),              nep_mo_g(nlat),
!      4      nbp_mo_g(nlat),              hetrores_mo_g(nlat),
!      5      autores_mo_g(nlat),          litres_mo_g(nlat),
!      6      soilcres_mo_g(nlat),
!      7      vgbiomas_mo_g(nlat),         totcmass_mo_g(nlat),
! !
!      8      emit_co2_mo_g(nlat),       emit_co_mo_g(nlat),
!      9      emit_ch4_mo_g(nlat),       emit_nmhc_mo_g(nlat),
!      a      emit_h2_mo_g(nlat),        emit_nox_mo_g(nlat),
!      b      emit_n2o_mo_g(nlat),       emit_pm25_mo_g(nlat),
!      c      emit_tpm_mo_g(nlat),       emit_tc_mo_g(nlat),
!      d      emit_oc_mo_g(nlat),        emit_bc_mo_g(nlat),
!      e      probfire_mo_g(nlat),       luc_emc_mo_g(nlat),
!      f      lucltrin_mo_g(nlat),       lucsocin_mo_g(nlat),
!      g      burnfrac_mo_g(nlat),       bterm_mo_g(nlat),
!      h      lterm_mo_g(nlat),          mterm_mo_g(nlat),
! !
!      &      ch4wet1_mo_g(nlat),        ch4wet2_mo_g(nlat),     !Rudra added on 03/12/2013
!      &      wetfdyn_mo_g(nlat),        ch4dyn1_mo_g(nlat),
!      &      ch4dyn2_mo_g(nlat)

        real, pointer, dimension(:) :: laimaxg_mo_g
        real, pointer, dimension(:) :: stemmass_mo_g
        real, pointer, dimension(:) :: rootmass_mo_g
        real, pointer, dimension(:) :: litrmass_mo_g
        real, pointer, dimension(:) :: soilcmas_mo_g
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
        real, pointer, dimension(:) :: probfire_mo_g
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

!      Mosaic monthly variables (denoted by name ending in "_mo_m")

!        real laimaxg_mo_m(nlat,nmos,icc),   stemmass_mo_m(nlat,nmos,icc),
!      1      rootmass_mo_m(nlat,nmos,icc),litrmass_mo_m(nlat,nmos,iccp1),
!      2      soilcmas_mo_m(nlat,nmos,iccp1),  npp_mo_m(nlat,nmos,icc),
!      3      gpp_mo_m(nlat,nmos,icc),       nep_mo_m(nlat,nmos,iccp1),
!      4      nbp_mo_m(nlat,nmos,iccp1),     vgbiomas_mo_m(nlat,nmos,icc),
!      a      hetrores_mo_m(nlat,nmos,iccp1), autores_mo_m(nlat,nmos,icc),
!      b      litres_mo_m(nlat,nmos,iccp1),soilcres_mo_m(nlat,nmos,iccp1),
!      c      totcmass_mo_m(nlat,nmos,icc)
! c
!        real emit_co2_mo_m(nlat,nmos,icc),emit_co_mo_m(nlat,nmos,icc),
!      1      emit_ch4_mo_m(nlat,nmos,icc),emit_nmhc_mo_m(nlat,nmos,icc),
!      2      emit_h2_mo_m(nlat,nmos,icc), emit_nox_mo_m(nlat,nmos,icc),
!      3      emit_n2o_mo_m(nlat,nmos,icc),emit_pm25_mo_m(nlat,nmos,icc),
!      4      emit_tpm_mo_m(nlat,nmos,icc),emit_tc_mo_m(nlat,nmos,icc),
!      5      emit_oc_mo_m(nlat,nmos,icc), emit_bc_mo_m(nlat,nmos,icc),
!      6      probfire_mo_m(nlat,nmos),    bterm_mo_m(nlat,nmos),
!      7      luc_emc_mo_m(nlat,nmos),     lterm_mo_m(nlat,nmos),
!      8      lucsocin_mo_m(nlat,nmos),    mterm_mo_m(nlat,nmos),
!      9      lucltrin_mo_m(nlat,nmos),
!      a      burnfrac_mo_m(nlat,nmos,icc),
!      &      ch4wet1_mo_m(nlat,nmos), ch4wet2_mo_m(nlat,nmos),
!      &      wetfdyn_mo_m(nlat,nmos),
!      &      ch4dyn1_mo_m(nlat,nmos), ch4dyn2_mo_m(nlat,nmos)
c
      real, pointer, dimension(:,:,:) :: laimaxg_mo_m
      real, pointer, dimension(:,:,:) :: stemmass_mo_m
      real, pointer, dimension(:,:,:) :: rootmass_mo_m
      real, pointer, dimension(:,:,:) :: npp_mo_m
      real, pointer, dimension(:,:,:) :: gpp_mo_m
      real, pointer, dimension(:,:,:) :: vgbiomas_mo_m
      real, pointer, dimension(:,:,:) :: autores_mo_m
      real, pointer, dimension(:,:,:) :: totcmass_mo_m
      real, pointer, dimension(:,:,:) :: litrmass_mo_m
      real, pointer, dimension(:,:,:) :: soilcmas_mo_m
      real, pointer, dimension(:,:,:) :: nep_mo_m
      real, pointer, dimension(:,:,:) :: litres_mo_m
      real, pointer, dimension(:,:,:) :: soilcres_mo_m
      real, pointer, dimension(:,:,:) :: hetrores_mo_m
      real, pointer, dimension(:,:,:) :: nbp_mo_m
      real, pointer, dimension(:,:,:) :: emit_co2_mo_m
      real, pointer, dimension(:,:,:) :: emit_co_mo_m
      real, pointer, dimension(:,:,:) :: emit_ch4_mo_m
      real, pointer, dimension(:,:,:) :: emit_nmhc_mo_m
      real, pointer, dimension(:,:,:) :: emit_h2_mo_m
      real, pointer, dimension(:,:,:) :: emit_nox_mo_m
      real, pointer, dimension(:,:,:) :: emit_n2o_mo_m
      real, pointer, dimension(:,:,:) :: emit_pm25_mo_m
      real, pointer, dimension(:,:,:) :: emit_tpm_mo_m
      real, pointer, dimension(:,:,:) :: emit_tc_mo_m
      real, pointer, dimension(:,:,:) :: emit_oc_mo_m
      real, pointer, dimension(:,:,:) :: emit_bc_mo_m
      real, pointer, dimension(:,:,:) :: burnfrac_mo_m
      real, pointer, dimension(:,:) :: probfire_mo_m
      real, pointer, dimension(:,:) :: bterm_mo_m
      real, pointer, dimension(:,:) :: luc_emc_mo_m
      real, pointer, dimension(:,:) :: lterm_mo_m
      real, pointer, dimension(:,:) :: lucsocin_mo_m
      real, pointer, dimension(:,:) :: mterm_mo_m
      real, pointer, dimension(:,:) :: lucltrin_mo_m
      real, pointer, dimension(:,:) :: ch4wet1_mo_m
      real, pointer, dimension(:,:) :: ch4wet2_mo_m
      real, pointer, dimension(:,:) :: wetfdyn_mo_m
      real, pointer, dimension(:,:) :: ch4dyn1_mo_m
      real, pointer, dimension(:,:) :: ch4dyn2_mo_m

!     -----------------------
c      Annual output for CTEM grid-averaged variables:
c      (denoted by name ending in "_yr_g")

!        real laimaxg_yr_g(nlat),        stemmass_yr_g(nlat),
!      1      rootmass_yr_g(nlat),       litrmass_yr_g(nlat),
!      2      soilcmas_yr_g(nlat),       npp_yr_g(nlat),
!      3      gpp_yr_g(nlat),            nep_yr_g(nlat),
!      a      nbp_yr_g(nlat),            hetrores_yr_g(nlat),
!      b      autores_yr_g(nlat),        litres_yr_g(nlat),
!      c      soilcres_yr_g(nlat),
!      4      vgbiomas_yr_g(nlat),       totcmass_yr_g(nlat),
!      5      emit_co2_yr_g(nlat),       emit_co_yr_g(nlat),
!      6      emit_ch4_yr_g(nlat),       emit_nmhc_yr_g(nlat),
!      7      emit_h2_yr_g(nlat),        emit_nox_yr_g(nlat),
!      8      emit_n2o_yr_g(nlat),       emit_pm25_yr_g(nlat),
!      9      emit_tpm_yr_g(nlat),       emit_tc_yr_g(nlat),
!      a      emit_oc_yr_g(nlat),        emit_bc_yr_g(nlat),
!      b      probfire_yr_g(nlat),       luc_emc_yr_g(nlat),
!      c      lucltrin_yr_g(nlat),       lucsocin_yr_g(nlat),
!      d      burnfrac_yr_g(nlat),       bterm_yr_g(nlat),
!      e      lterm_yr_g(nlat),          mterm_yr_g(nlat),
!      &      ch4wet1_yr_g(nlat),        ch4wet2_yr_g(nlat),
!      &      wetfdyn_yr_g(nlat),        ch4dyn1_yr_g(nlat),
!      &      ch4dyn2_yr_g(nlat)

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
      real, pointer, dimension(:) :: probfire_yr_g
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


! c      Annual output for CTEM mosaic variables:
! c      (denoted by name ending in "_yr_m")
!
!        real laimaxg_yr_m(nlat,nmos,icc), stemmass_yr_m(nlat,nmos,icc),
!      1      rootmass_yr_m(nlat,nmos,icc),litrmass_yr_m(nlat,nmos,iccp1),
!      2      soilcmas_yr_m(nlat,nmos,iccp1),npp_yr_m(nlat,nmos,icc),
!      3      gpp_yr_m(nlat,nmos,icc),       nep_yr_m(nlat,nmos,iccp1),
!      4      nbp_yr_m(nlat,nmos,iccp1),     vgbiomas_yr_m(nlat,nmos,icc),
!      a      hetrores_yr_m(nlat,nmos,iccp1),autores_yr_m(nlat,nmos,icc),
!      b      litres_yr_m(nlat,nmos,iccp1),soilcres_yr_m(nlat,nmos,iccp1),
!      5      emit_co2_yr_m(nlat,nmos,icc), emit_co_yr_m(nlat,nmos,icc),
!      6      emit_ch4_yr_m(nlat,nmos,icc), emit_nmhc_yr_m(nlat,nmos,icc),
!      7      emit_h2_yr_m(nlat,nmos,icc), emit_nox_yr_m(nlat,nmos,icc),
!      8      emit_n2o_yr_m(nlat,nmos,icc),emit_pm25_yr_m(nlat,nmos,icc),
!      9      emit_tpm_yr_m(nlat,nmos,icc),emit_tc_yr_m(nlat,nmos,icc),
!      a      emit_oc_yr_m(nlat,nmos,icc), emit_bc_yr_m(nlat,nmos,icc),
!      b      probfire_yr_m(nlat,nmos),    bterm_yr_m(nlat,nmos),
!      c      luc_emc_yr_m(nlat,nmos),     totcmass_yr_m(nlat,nmos,icc),
!      d      lucsocin_yr_m(nlat,nmos),    lterm_yr_m(nlat,nmos),
!      e      lucltrin_yr_m(nlat,nmos),    mterm_yr_m(nlat,nmos),
!      f      burnfrac_yr_m(nlat,nmos,icc),
!      &      ch4wet1_yr_m(nlat,nmos),   ch4wet2_yr_m(nlat,nmos),
!      &      wetfdyn_yr_m(nlat,nmos),   ch4dyn1_yr_m(nlat,nmos),
!      &      ch4dyn2_yr_m(nlat,nmos)

      real, pointer, dimension(:,:,:) :: laimaxg_yr_m
      real, pointer, dimension(:,:,:) :: stemmass_yr_m
      real, pointer, dimension(:,:,:) :: rootmass_yr_m
      real, pointer, dimension(:,:,:) :: npp_yr_m
      real, pointer, dimension(:,:,:) :: gpp_yr_m
      real, pointer, dimension(:,:,:) :: vgbiomas_yr_m
      real, pointer, dimension(:,:,:) :: autores_yr_m
      real, pointer, dimension(:,:,:) :: totcmass_yr_m
      real, pointer, dimension(:,:,:) :: litrmass_yr_m
      real, pointer, dimension(:,:,:) :: soilcmas_yr_m
      real, pointer, dimension(:,:,:) :: nep_yr_m
      real, pointer, dimension(:,:,:) :: litres_yr_m
      real, pointer, dimension(:,:,:) :: soilcres_yr_m
      real, pointer, dimension(:,:,:) :: hetrores_yr_m
      real, pointer, dimension(:,:,:) :: nbp_yr_m
      real, pointer, dimension(:,:,:) :: emit_co2_yr_m
      real, pointer, dimension(:,:,:) :: emit_co_yr_m
      real, pointer, dimension(:,:,:) :: emit_ch4_yr_m
      real, pointer, dimension(:,:,:) :: emit_nmhc_yr_m
      real, pointer, dimension(:,:,:) :: emit_h2_yr_m
      real, pointer, dimension(:,:,:) :: emit_nox_yr_m
      real, pointer, dimension(:,:,:) :: emit_n2o_yr_m
      real, pointer, dimension(:,:,:) :: emit_pm25_yr_m
      real, pointer, dimension(:,:,:) :: emit_tpm_yr_m
      real, pointer, dimension(:,:,:) :: emit_tc_yr_m
      real, pointer, dimension(:,:,:) :: emit_oc_yr_m
      real, pointer, dimension(:,:,:) :: emit_bc_yr_m
      real, pointer, dimension(:,:,:) :: burnfrac_yr_m
      real, pointer, dimension(:,:) :: probfire_yr_m
      real, pointer, dimension(:,:) :: bterm_yr_m
      real, pointer, dimension(:,:) :: luc_emc_yr_m
      real, pointer, dimension(:,:) :: lterm_yr_m
      real, pointer, dimension(:,:) :: lucsocin_yr_m
      real, pointer, dimension(:,:) :: mterm_yr_m
      real, pointer, dimension(:,:) :: lucltrin_yr_m
      real, pointer, dimension(:,:) :: ch4wet1_yr_m
      real, pointer, dimension(:,:) :: ch4wet2_yr_m
      real, pointer, dimension(:,:) :: wetfdyn_yr_m
      real, pointer, dimension(:,:) :: ch4dyn1_yr_m
      real, pointer, dimension(:,:) :: ch4dyn2_yr_m

c
c============= CTEM array declaration done =============================/
C
C=======================================================================
C     * PHYSICAL CONSTANTS.
C     * PARAMETERS IN THE FOLLOWING COMMON BLOCKS ARE NORMALLY DEFINED
C     * WITHIN THE GCM.

      COMMON /PARAMS/ X1,    X2,    X3,    X4,   G,GAS,   X5,
     1                X6,    CPRES, GASV,  X7
      COMMON /PARAM1/ CPI,   X8,    CELZRO,X9,    X10,    X11
      COMMON /PARAM3/ X12,   X13,   X14,   X15,   SIGMA,  X16
      COMMON  /TIMES/ DELTIM,K1,    K2,    K3,    K4,     K5,
     1                K6,    K7,    K8,    K9,    K10,    K11
!      COMMON /EPS/    A,B,EPS1,EPS2   !FLAG not in new CLASS
!      COMMON /HTCP/   T1S,T2S,AI,BI,AW,BW,SLP !FLAG not in new CLASS
C
C     * THE FOLLOWING COMMON BLOCKS ARE DEFINED SPECIFICALLY FOR USE
C     * IN CLASS, VIA BLOCK DATA AND THE SUBROUTINE "CLASSD".
C
      COMMON /CLASS1/ DELT,TFREZ
      COMMON /CLASS2/ RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN
      COMMON /CLASS3/ TCW,TCICE,TCSAND,TCCLAY,TCOM,TCDRYS,
     1                RHOSOL,RHOOM
      COMMON /CLASS4/ HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,
     1                SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,
     2                TCGLAC,CLHMLT,CLHVAP
      COMMON /CLASS5/ THPORG,THRORG,THMORG,BORG,PSISORG,GRKSORG
      COMMON /CLASS6/ PI,GROWYR,ZOLNG,ZOLNS,ZOLNI,ZORAT,ZORATG
      COMMON /CLASS7/ CANEXT,XLEAF
      COMMON /CLASS8/ ALVSI,ALIRI,ALVSO,ALIRO,ALBRCK
      COMMON /PHYCON/ DELTA,CGRAV,CKARM,CPD
      COMMON /CLASSD2/ AS,ASX,CI,BS,BETA,FACTN,HMIN,ANGMAX
!      COMMON /ESTWI/  RW1,RW2,RW3,RI1,RI2,RI3 !FLAG not in new CLASS
C
      CALL CLASSD
C
      ZDMROW(1)=10.0
      ZDHROW(1)=2.0
      NTLD=NMOS
      CUMSNO = 0.0
C
C===================== CTEM ==============================================\

!       FLAG testing!:

    ! Point pointers

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
      TAACC_MO          => class_out%TAACC_MO
      FSSTAR_MO         => class_out%FSSTAR_MO
      FLSTAR_MO         => class_out%FLSTAR_MO
      QH_MO             => class_out%QH_MO
      QE_MO             => class_out%QE_MO
      TBARACC_MO        => class_out%TBARACC_MO
      THLQACC_MO        => class_out%THLQACC_MO
      THICACC_MO        => class_out%THICACC_MO
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
      TAACC_YR          => class_out%TAACC_YR
      FSSTAR_YR         => class_out%FSSTAR_YR
      FLSTAR_YR         => class_out%FLSTAR_YR
      QH_YR             => class_out%QH_YR
      QE_YR             => class_out%QE_YR

      ctem_on           => c_switch%ctem_on
      parallelrun       => c_switch%parallelrun
      mosaic            => c_switch%mosaic
      cyclemet          => c_switch%cyclemet
      dofire            => c_switch%dofire
      run_model         => c_switch%run_model
      met_rewound       => c_switch%met_rewound
      reach_eof         => c_switch%reach_eof
      compete           => c_switch%compete
      start_bare        => c_switch%start_bare
      rsfile            => c_switch%rsfile
      lnduseon          => c_switch%lnduseon
      co2on             => c_switch%co2on
      popdon            => c_switch%popdon
      inibioclim        => c_switch%inibioclim
      start_from_rs     => c_switch%start_from_rs
      dowetlands        => c_switch%dowetlands
      obswetf           => c_switch%obswetf
      transient_run     => c_switch%transient_run

      tcanrs            => vrot%tcanrs
      tsnors            => vrot%tsnors
      tpndrs            => vrot%tpndrs
      csum              => vrot%csum
      tbaraccrow_m      => vrot%tbaraccrow_m
      tcanoaccrow_m     => vrot%tcanoaccrow_m
      uvaccrow_m        => vrot%uvaccrow_m
      vvaccrow_m        => vrot%vvaccrow_m

      ! ROW:
      ailcminrow        => vrot%ailcmin
      ailcmaxrow        => vrot%ailcmax
      dvdfcanrow        => vrot%dvdfcan
      gleafmasrow       => vrot%gleafmas
      bleafmasrow       => vrot%bleafmas
      stemmassrow       => vrot%stemmass
      rootmassrow       => vrot%rootmass
      pstemmassrow      => vrot%pstemmass
      pgleafmassrow     => vrot%pgleafmass
      fcancmxrow        => vrot%fcancmx
      gavglairow        => vrot%gavglai
      zolncrow          => vrot%zolnc
      ailcrow           => vrot%ailc
      ailcgrow          => vrot%ailcg
      ailcgsrow         => vrot%ailcgs
      fcancsrow         => vrot%fcancs
      fcancrow          => vrot%fcanc
      co2concrow        => vrot%co2conc
      co2i1cgrow        => vrot%co2i1cg
      co2i1csrow        => vrot%co2i1cs
      co2i2cgrow        => vrot%co2i2cg
      co2i2csrow        => vrot%co2i2cs
      ancsvegrow        => vrot%ancsveg
      ancgvegrow        => vrot%ancgveg
      rmlcsvegrow       => vrot%rmlcsveg
      rmlcgvegrow       => vrot%rmlcgveg
      slairow           => vrot%slai
      ailcbrow          => vrot%ailcb
      canresrow         => vrot%canres
      flhrlossrow       => vrot%flhrloss

      tcanoaccrow_out   => vrot%tcanoaccrow_out
      qevpacc_m_save    => vrot%qevpacc_m_save

      grwtheffrow       => vrot%grwtheff
      lystmmasrow       => vrot%lystmmas
      lyrotmasrow       => vrot%lyrotmas
      tymaxlairow       => vrot%tymaxlai
      vgbiomasrow       => vrot%vgbiomas
      gavgltmsrow       => vrot%gavgltms
      gavgscmsrow       => vrot%gavgscms
      stmhrlosrow       => vrot%stmhrlos
      rmatcrow          => vrot%rmatc
      rmatctemrow       => vrot%rmatctem
      litrmassrow       => vrot%litrmass
      soilcmasrow       => vrot%soilcmas
      vgbiomas_vegrow   => vrot%vgbiomas_veg

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
      probfirerow       => vrot%probfire
      btermrow          => vrot%bterm
      ltermrow          => vrot%lterm
      mtermrow          => vrot%mterm

      extnprobgrd       => vrot%extnprob
      prbfrhucgrd       => vrot%prbfrhuc
      mlightnggrd       => vrot%mlightng

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

      wetfracgrd        => vrot%wetfrac
      ch4wet1row        => vrot%ch4wet1
      ch4wet2row        => vrot%ch4wet2
      wetfdynrow        => vrot%wetfdyn
      ch4dyn1row        => vrot%ch4dyn1
      ch4dyn2row        => vrot%ch4dyn2
      wetfrac_mon       => vrot%wetfrac_mon

      lucemcomrow       => vrot%lucemcom
      lucltrinrow       => vrot%lucltrin
      lucsocinrow       => vrot%lucsocin

      npprow            => vrot%npp
      neprow            => vrot%nep
      nbprow            => vrot%nbp
      gpprow            => vrot%gpp
      hetroresrow       => vrot%hetrores
      autoresrow        => vrot%autores
      soilcresprow      => vrot%soilcresp
      rmrow             => vrot%rm
      rgrow             => vrot%rg
      litresrow         => vrot%litres
      socresrow         => vrot%socres
      dstcemlsrow       => vrot%dstcemls
      litrfallrow       => vrot%litrfall
      humiftrsrow       => vrot%humiftrs

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

      rothrlosrow       => vrot%rothrlos
      pfcancmxrow       => vrot%pfcancmx
      nfcancmxrow       => vrot%nfcancmx
      alvsctmrow        => vrot%alvsctm
      paicrow           => vrot%paic
      slaicrow          => vrot%slaic
      alirctmrow        => vrot%alirctm
      cfluxcgrow        => vrot%cfluxcg
      cfluxcsrow        => vrot%cfluxcs
      dstcemls3row      => vrot%dstcemls3
      anvegrow          => vrot%anveg
      rmlvegrow         => vrot%rmlveg

      pftexistrow       => vrot%pftexist
      colddaysrow       => vrot%colddays
      icountrow         => vrot%icount
      lfstatusrow       => vrot%lfstatus
      pandaysrow        => vrot%pandays
      stdalngrd         => vrot%stdaln


      ! >>>>>>>>>>>>>>>>>>>>>>>>>>
      ! GAT:

      lightng           => vgat%lightng
      tcanoaccgat_out   => vgat%tcanoaccgat_out

      ailcmingat        => vgat%ailcmin
      ailcmaxgat        => vgat%ailcmax
      dvdfcangat        => vgat%dvdfcan
      gleafmasgat       => vgat%gleafmas
      bleafmasgat       => vgat%bleafmas
      stemmassgat       => vgat%stemmass
      rootmassgat       => vgat%rootmass
      pstemmassgat      => vgat%pstemmass
      pgleafmassgat     => vgat%pgleafmass
      fcancmxgat        => vgat%fcancmx
      gavglaigat        => vgat%gavglai
      zolncgat          => vgat%zolnc
      ailcgat           => vgat%ailc
      ailcggat          => vgat%ailcg
      ailcgsgat         => vgat%ailcgs
      fcancsgat         => vgat%fcancs
      fcancgat          => vgat%fcanc
      co2concgat        => vgat%co2conc
      co2i1cggat        => vgat%co2i1cg
      co2i1csgat        => vgat%co2i1cs
      co2i2cggat        => vgat%co2i2cg
      co2i2csgat        => vgat%co2i2cs
      ancsveggat        => vgat%ancsveg
      ancgveggat        => vgat%ancgveg
      rmlcsveggat       => vgat%rmlcsveg
      rmlcgveggat       => vgat%rmlcgveg
      slaigat           => vgat%slai
      ailcbgat          => vgat%ailcb
      canresgat         => vgat%canres
      flhrlossgat       => vgat%flhrloss

      grwtheffgat       => vgat%grwtheff
      lystmmasgat       => vgat%lystmmas
      lyrotmasgat       => vgat%lyrotmas
      tymaxlaigat       => vgat%tymaxlai
      vgbiomasgat       => vgat%vgbiomas
      gavgltmsgat       => vgat%gavgltms
      gavgscmsgat       => vgat%gavgscms
      stmhrlosgat       => vgat%stmhrlos
      rmatcgat          => vgat%rmatc
      rmatctemgat       => vgat%rmatctem
      litrmassgat       => vgat%litrmass
      soilcmasgat       => vgat%soilcmas
      vgbiomas_veggat   => vgat%vgbiomas_veg

      emit_co2gat       => vgat%emit_co2
      emit_cogat        => vgat%emit_co
      emit_ch4gat       => vgat%emit_ch4
      emit_nmhcgat      => vgat%emit_nmhc
      emit_h2gat        => vgat%emit_h2
      emit_noxgat       => vgat%emit_nox
      emit_n2ogat       => vgat%emit_n2o
      emit_pm25gat      => vgat%emit_pm25
      emit_tpmgat       => vgat%emit_tpm
      emit_tcgat        => vgat%emit_tc
      emit_ocgat        => vgat%emit_oc
      emit_bcgat        => vgat%emit_bc
      burnfracgat       => vgat%burnfrac
      burnvegfgat       => vgat%burnvegf
      probfiregat       => vgat%probfire
      btermgat          => vgat%bterm
      ltermgat          => vgat%lterm
      mtermgat          => vgat%mterm

      extnprobgat       => vgat%extnprob
      prbfrhucgat       => vgat%prbfrhuc
      mlightnggat       => vgat%mlightng

      bmasveggat        => vgat%bmasveg
      cmasvegcgat       => vgat%cmasvegc
      veghghtgat        => vgat%veghght
      rootdpthgat       => vgat%rootdpth
      rmlgat            => vgat%rml
      rmsgat            => vgat%rms
      tltrleafgat       => vgat%tltrleaf
      tltrstemgat       => vgat%tltrstem
      tltrrootgat       => vgat%tltrroot
      leaflitrgat       => vgat%leaflitr
      roottempgat       => vgat%roottemp
      afrleafgat        => vgat%afrleaf
      afrstemgat        => vgat%afrstem
      afrrootgat        => vgat%afrroot
      wtstatusgat       => vgat%wtstatus
      ltstatusgat       => vgat%ltstatus
      rmrgat            => vgat%rmr

      wetfrac_sgrd      => vgat%wetfrac_s
      ch4wet1gat        => vgat%ch4wet1
      ch4wet2gat        => vgat%ch4wet2
      wetfdyngat        => vgat%wetfdyn
      ch4dyn1gat        => vgat%ch4dyn1
      ch4dyn2gat        => vgat%ch4dyn2

      lucemcomgat       => vgat%lucemcom
      lucltringat       => vgat%lucltrin
      lucsocingat       => vgat%lucsocin

      nppgat            => vgat%npp
      nepgat            => vgat%nep
      nbpgat            => vgat%nbp
      gppgat            => vgat%gpp
      hetroresgat       => vgat%hetrores
      autoresgat        => vgat%autores
      soilcrespgat      => vgat%soilcresp
      rmgat             => vgat%rm
      rggat             => vgat%rg
      litresgat         => vgat%litres
      socresgat         => vgat%socres
      dstcemlsgat       => vgat%dstcemls
      litrfallgat       => vgat%litrfall
      humiftrsgat       => vgat%humiftrs

      gppveggat         => vgat%gppveg
      nepveggat         => vgat%nepveg
      nbpveggat         => vgat%nbpveg
      nppveggat         => vgat%nppveg
      hetroresveggat    => vgat%hetroresveg
      autoresveggat     => vgat%autoresveg
      litresveggat      => vgat%litresveg
      soilcresveggat    => vgat%soilcresveg
      rmlvegaccgat      => vgat%rmlvegacc
      rmsveggat         => vgat%rmsveg
      rmrveggat         => vgat%rmrveg
      rgveggat          => vgat%rgveg

      rothrlosgat       => vgat%rothrlos
      pfcancmxgat       => vgat%pfcancmx
      nfcancmxgat       => vgat%nfcancmx
      alvsctmgat        => vgat%alvsctm
      paicgat           => vgat%paic
      slaicgat          => vgat%slaic
      alirctmgat        => vgat%alirctm
      cfluxcggat        => vgat%cfluxcg
      cfluxcsgat        => vgat%cfluxcs
      dstcemls3gat      => vgat%dstcemls3
      anveggat          => vgat%anveg
      rmlveggat         => vgat%rmlveg

      twarmm            => vgat%twarmm
      tcoldm            => vgat%tcoldm
      gdd5              => vgat%gdd5
      aridity           => vgat%aridity
      srplsmon          => vgat%srplsmon
      defctmon          => vgat%defctmon
      anndefct          => vgat%anndefct
      annsrpls          => vgat%annsrpls
      annpcp            => vgat%annpcp
      dry_season_length => vgat%dry_season_length

      pftexistgat       => vgat%pftexist
      colddaysgat       => vgat%colddays
      icountgat         => vgat%icount
      lfstatusgat       => vgat%lfstatus
      pandaysgat        => vgat%pandays
      stdalngat         => vgat%stdaln

      ! Mosaic-level:

      PREACC_M          => vrot%PREACC_M
      GTACC_M           => vrot%GTACC_M
      QEVPACC_M         => vrot%QEVPACC_M
      HFSACC_M          => vrot%HFSACC_M
      HMFNACC_M         => vrot%HMFNACC_M
      ROFACC_M          => vrot%ROFACC_M
      SNOACC_M          => vrot%SNOACC_M
      OVRACC_M          => vrot%OVRACC_M
      WTBLACC_M         => vrot%WTBLACC_M
      TBARACC_M         => vrot%TBARACC_M
      THLQACC_M         => vrot%THLQACC_M
      THICACC_M         => vrot%THICACC_M
      THALACC_M         => vrot%THALACC_M
      ALVSACC_M         => vrot%ALVSACC_M
      ALIRACC_M         => vrot%ALIRACC_M
      RHOSACC_M         => vrot%RHOSACC_M
      TSNOACC_M         => vrot%TSNOACC_M
      WSNOACC_M         => vrot%WSNOACC_M
      SNOARE_M          => vrot%SNOARE_M
      TCANACC_M         => vrot%TCANACC_M
      RCANACC_M         => vrot%RCANACC_M
      SCANACC_M         => vrot%SCANACC_M
      GROACC_M          => vrot%GROACC_M
      FSINACC_M         => vrot%FSINACC_M
      FLINACC_M         => vrot%FLINACC_M
      TAACC_M           => vrot%TAACC_M
      UVACC_M           => vrot%UVACC_M
      PRESACC_M         => vrot%PRESACC_M
      QAACC_M           => vrot%QAACC_M
      EVAPACC_M         => vrot%EVAPACC_M
      FLUTACC_M         => vrot%FLUTACC_M

      ! grid-averaged

      WSNOROT_g         => ctem_grd%WSNOROT_g
      ROFSROT_g         => ctem_grd%ROFSROT_g
      SNOROT_g          => ctem_grd%SNOROT_g
      RHOSROT_g         => ctem_grd%RHOSROT_g
      ROFROT_g          => ctem_grd%ROFROT_g
      ZPNDROT_g         => ctem_grd%ZPNDROT_g
      RCANROT_g         => ctem_grd%RCANROT_g
      SCANROT_g         => ctem_grd%SCANROT_g
      TROFROT_g         => ctem_grd%TROFROT_g
      TROOROT_g         => ctem_grd%TROOROT_g
      TROBROT_g         => ctem_grd%TROBROT_g
      ROFOROT_g         => ctem_grd%ROFOROT_g
      ROFBROT_g         => ctem_grd%ROFBROT_g
      TROSROT_g         => ctem_grd%TROSROT_g
      FSGVROT_g         => ctem_grd%FSGVROT_g
      FSGSROT_g         => ctem_grd%FSGSROT_g
      FLGVROT_g         => ctem_grd%FLGVROT_g
      FLGSROT_g         => ctem_grd%FLGSROT_g
      HFSCROT_g         => ctem_grd%HFSCROT_g
      HFSSROT_g         => ctem_grd%HFSSROT_g
      HEVCROT_g         => ctem_grd%HEVCROT_g
      HEVSROT_g         => ctem_grd%HEVSROT_g
      HMFCROT_g         => ctem_grd%HMFCROT_g
      HMFNROT_g         => ctem_grd%HMFNROT_g
      HTCSROT_g         => ctem_grd%HTCSROT_g
      HTCCROT_g         => ctem_grd%HTCCROT_g
      FSGGROT_g         => ctem_grd%FSGGROT_g
      FLGGROT_g         => ctem_grd%FLGGROT_g
      HFSGROT_g         => ctem_grd%HFSGROT_g
      HEVGROT_g         => ctem_grd%HEVGROT_g
      CDHROT_g          => ctem_grd%CDHROT_g
      CDMROT_g          => ctem_grd%CDMROT_g
      SFCUROT_g         => ctem_grd%SFCUROT_g
      SFCVROT_g         => ctem_grd%SFCVROT_g
      fc_g              => ctem_grd%fc_g
      fg_g              => ctem_grd%fg_g
      fcs_g             => ctem_grd%fcs_g
      fgs_g             => ctem_grd%fgs_g
      PCFCROT_g         => ctem_grd%PCFCROT_g
      PCLCROT_g         => ctem_grd%PCLCROT_g
      PCPGROT_g         => ctem_grd%PCPGROT_g
      QFCFROT_g         => ctem_grd%QFCFROT_g
      QFGROT_g          => ctem_grd%QFGROT_g
      QFCROT_g          => ctem_grd%QFCROT_g
      ROFCROT_g         => ctem_grd%ROFCROT_g
      ROFNROT_g         => ctem_grd%ROFNROT_g
      WTRSROT_g         => ctem_grd%WTRSROT_g
      WTRGROT_g         => ctem_grd%WTRGROT_g
      PCPNROT_g         => ctem_grd%PCPNROT_g
      QFCLROT_g         => ctem_grd%QFCLROT_g
      QFNROT_g          => ctem_grd%QFNROT_g
      WTRCROT_g         => ctem_grd%WTRCROT_g
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
      probfire_g        => ctem_grd%probfire_g
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
      afrleaf_g         => ctem_grd%afrleaf_g
      afrstem_g         => ctem_grd%afrstem_g
      afrroot_g         => ctem_grd%afrroot_g
      lfstatus_g        => ctem_grd%lfstatus_g
      rmlvegrow_g       => ctem_grd%rmlvegrow_g
      anvegrow_g        => ctem_grd%anvegrow_g
      rmatctem_g        => ctem_grd%rmatctem_g
      HMFGROT_g         => ctem_grd%HMFGROT_g
      HTCROT_g          => ctem_grd%HTCROT_g
      TBARROT_g         => ctem_grd%TBARROT_g
      THLQROT_g         => ctem_grd%THLQROT_g
      THICROT_g         => ctem_grd%THICROT_g
      GFLXROT_g         => ctem_grd%GFLXROT_g

       fsstar_g         => ctem_grd%fsstar_g
       flstar_g         => ctem_grd%flstar_g
       qh_g             => ctem_grd%qh_g
       qe_g             => ctem_grd%qe_g
       snomlt_g         => ctem_grd%snomlt_g
       beg_g            => ctem_grd%beg_g
       gtout_g          => ctem_grd%gtout_g
       tpn_g            => ctem_grd%tpn_g
       altot_g          => ctem_grd%altot_g
       tcn_g            => ctem_grd%tcn_g
       tsn_g            => ctem_grd%tsn_g
       zsn_g            => ctem_grd%zsn_g

      ! mosaic level variables:

      leaflitr_m        => ctem_tile%leaflitr_m
      tltrleaf_m        => ctem_tile%tltrleaf_m
      tltrstem_m        => ctem_tile%tltrstem_m
      tltrroot_m        => ctem_tile%tltrroot_m
      ailcg_m           => ctem_tile%ailcg_m
      ailcb_m           => ctem_tile%ailcb_m
      rmatctem_m        => ctem_tile%rmatctem_m
      veghght_m         => ctem_tile%veghght_m
      rootdpth_m        => ctem_tile%rootdpth_m
      roottemp_m        => ctem_tile%roottemp_m
      slai_m            => ctem_tile%slai_m
      afrroot_m         => ctem_tile%afrroot_m
      afrleaf_m         => ctem_tile%afrleaf_m
      afrstem_m         => ctem_tile%afrstem_m
      laimaxg_m         => ctem_tile%laimaxg_m
      stemmass_m        => ctem_tile%stemmass_m
      rootmass_m        => ctem_tile%rootmass_m
      litrmass_m        => ctem_tile%litrmass_m
      gleafmas_m        => ctem_tile%gleafmas_m
      bleafmas_m        => ctem_tile%bleafmas_m
      soilcmas_m        => ctem_tile%soilcmas_m
      fsnowacc_m        => ctem_tile%fsnowacc_m
      tcansacc_m        => ctem_tile%tcansacc_m
      tcanoaccgat_m     => ctem_tile%tcanoaccgat_m
      taaccgat_m        => ctem_tile%taaccgat_m
      uvaccgat_m        => ctem_tile%uvaccgat_m
      vvaccgat_m        => ctem_tile%vvaccgat_m
      tbaraccgat_m      => ctem_tile%tbaraccgat_m
      tbarcacc_m        => ctem_tile%tbarcacc_m
      tbarcsacc_m       => ctem_tile%tbarcsacc_m
      tbargacc_m        => ctem_tile%tbargacc_m
      tbargsacc_m       => ctem_tile%tbargsacc_m
      thliqcacc_m       => ctem_tile%thliqcacc_m
      thliqgacc_m       => ctem_tile%thliqgacc_m
      thicecacc_m       => ctem_tile%thicecacc_m
      ancsvgac_m        => ctem_tile%ancsvgac_m
      ancgvgac_m        => ctem_tile%ancgvgac_m
      rmlcsvga_m        => ctem_tile%rmlcsvga_m
      rmlcgvga_m        => ctem_tile%rmlcgvga_m
      ifcancmx_m        => ctem_tile%ifcancmx_m


      ! grid level monthly outputs

        laimaxg_mo_g        =>ctem_grd_mo%laimaxg_mo_g
        stemmass_mo_g       =>ctem_grd_mo%stemmass_mo_g
        rootmass_mo_g       =>ctem_grd_mo%rootmass_mo_g
        litrmass_mo_g       =>ctem_grd_mo%litrmass_mo_g
        soilcmas_mo_g       =>ctem_grd_mo%soilcmas_mo_g
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
        probfire_mo_g       =>ctem_grd_mo%probfire_mo_g
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

      ! mosaic monthly outputs

      laimaxg_mo_m          =>ctem_tile_mo%laimaxg_mo_m
      stemmass_mo_m         =>ctem_tile_mo%stemmass_mo_m
      rootmass_mo_m         =>ctem_tile_mo%rootmass_mo_m
      npp_mo_m              =>ctem_tile_mo%npp_mo_m
      gpp_mo_m              =>ctem_tile_mo%gpp_mo_m
      vgbiomas_mo_m         =>ctem_tile_mo%vgbiomas_mo_m
      autores_mo_m          =>ctem_tile_mo%autores_mo_m
      totcmass_mo_m         =>ctem_tile_mo%totcmass_mo_m
      litrmass_mo_m         =>ctem_tile_mo%litrmass_mo_m
      soilcmas_mo_m         =>ctem_tile_mo%soilcmas_mo_m
      nep_mo_m              =>ctem_tile_mo%nep_mo_m
      litres_mo_m           =>ctem_tile_mo%litres_mo_m
      soilcres_mo_m         =>ctem_tile_mo%soilcres_mo_m
      hetrores_mo_m         =>ctem_tile_mo%hetrores_mo_m
      nbp_mo_m              =>ctem_tile_mo%nbp_mo_m
      emit_co2_mo_m         =>ctem_tile_mo%emit_co2_mo_m
      emit_co_mo_m          =>ctem_tile_mo%emit_co_mo_m
      emit_ch4_mo_m         =>ctem_tile_mo%emit_ch4_mo_m
      emit_nmhc_mo_m        =>ctem_tile_mo%emit_nmhc_mo_m
      emit_h2_mo_m          =>ctem_tile_mo%emit_h2_mo_m
      emit_nox_mo_m         =>ctem_tile_mo%emit_nox_mo_m
      emit_n2o_mo_m         =>ctem_tile_mo%emit_n2o_mo_m
      emit_pm25_mo_m        =>ctem_tile_mo%emit_pm25_mo_m
      emit_tpm_mo_m         =>ctem_tile_mo%emit_tpm_mo_m
      emit_tc_mo_m          =>ctem_tile_mo%emit_tc_mo_m
      emit_oc_mo_m          =>ctem_tile_mo%emit_oc_mo_m
      emit_bc_mo_m          =>ctem_tile_mo%emit_bc_mo_m
      burnfrac_mo_m         =>ctem_tile_mo%burnfrac_mo_m
      probfire_mo_m         =>ctem_tile_mo%probfire_mo_m
      bterm_mo_m            =>ctem_tile_mo%bterm_mo_m
      luc_emc_mo_m          =>ctem_tile_mo%luc_emc_mo_m
      lterm_mo_m            =>ctem_tile_mo%lterm_mo_m
      lucsocin_mo_m         =>ctem_tile_mo%lucsocin_mo_m
      mterm_mo_m            =>ctem_tile_mo%mterm_mo_m
      lucltrin_mo_m         =>ctem_tile_mo%lucltrin_mo_m
      ch4wet1_mo_m          =>ctem_tile_mo%ch4wet1_mo_m
      ch4wet2_mo_m          =>ctem_tile_mo%ch4wet2_mo_m
      wetfdyn_mo_m          =>ctem_tile_mo%wetfdyn_mo_m
      ch4dyn1_mo_m          =>ctem_tile_mo%ch4dyn1_mo_m
      ch4dyn2_mo_m          =>ctem_tile_mo%ch4dyn2_mo_m

      ! grid level annual outputs
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
      probfire_yr_g         =>ctem_grd_yr%probfire_yr_g
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

      ! mosaic annual outputs

      laimaxg_yr_m          =>ctem_tile_yr%laimaxg_yr_m
      stemmass_yr_m         =>ctem_tile_yr%stemmass_yr_m
      rootmass_yr_m         =>ctem_tile_yr%rootmass_yr_m
      npp_yr_m              =>ctem_tile_yr%npp_yr_m
      gpp_yr_m              =>ctem_tile_yr%gpp_yr_m
      vgbiomas_yr_m         =>ctem_tile_yr%vgbiomas_yr_m
      autores_yr_m          =>ctem_tile_yr%autores_yr_m
      totcmass_yr_m         =>ctem_tile_yr%totcmass_yr_m
      litrmass_yr_m         =>ctem_tile_yr%litrmass_yr_m
      soilcmas_yr_m         =>ctem_tile_yr%soilcmas_yr_m
      nep_yr_m              =>ctem_tile_yr%nep_yr_m
      litres_yr_m           =>ctem_tile_yr%litres_yr_m
      soilcres_yr_m         =>ctem_tile_yr%soilcres_yr_m
      hetrores_yr_m         =>ctem_tile_yr%hetrores_yr_m
      nbp_yr_m              =>ctem_tile_yr%nbp_yr_m
      emit_co2_yr_m         =>ctem_tile_yr%emit_co2_yr_m
      emit_co_yr_m          =>ctem_tile_yr%emit_co_yr_m
      emit_ch4_yr_m         =>ctem_tile_yr%emit_ch4_yr_m
      emit_nmhc_yr_m        =>ctem_tile_yr%emit_nmhc_yr_m
      emit_h2_yr_m          =>ctem_tile_yr%emit_h2_yr_m
      emit_nox_yr_m         =>ctem_tile_yr%emit_nox_yr_m
      emit_n2o_yr_m         =>ctem_tile_yr%emit_n2o_yr_m
      emit_pm25_yr_m        =>ctem_tile_yr%emit_pm25_yr_m
      emit_tpm_yr_m         =>ctem_tile_yr%emit_tpm_yr_m
      emit_tc_yr_m          =>ctem_tile_yr%emit_tc_yr_m
      emit_oc_yr_m          =>ctem_tile_yr%emit_oc_yr_m
      emit_bc_yr_m          =>ctem_tile_yr%emit_bc_yr_m
      burnfrac_yr_m         =>ctem_tile_yr%burnfrac_yr_m
      probfire_yr_m         =>ctem_tile_yr%probfire_yr_m
      bterm_yr_m            =>ctem_tile_yr%bterm_yr_m
      luc_emc_yr_m          =>ctem_tile_yr%luc_emc_yr_m
      lterm_yr_m            =>ctem_tile_yr%lterm_yr_m
      lucsocin_yr_m         =>ctem_tile_yr%lucsocin_yr_m
      mterm_yr_m            =>ctem_tile_yr%mterm_yr_m
      lucltrin_yr_m         =>ctem_tile_yr%lucltrin_yr_m
      ch4wet1_yr_m          =>ctem_tile_yr%ch4wet1_yr_m
      ch4wet2_yr_m          =>ctem_tile_yr%ch4wet2_yr_m
      wetfdyn_yr_m          =>ctem_tile_yr%wetfdyn_yr_m
      ch4dyn1_yr_m          =>ctem_tile_yr%ch4dyn1_yr_m
      ch4dyn2_yr_m          =>ctem_tile_yr%ch4dyn2_yr_m


c     all model switches are read in from a namelist file
      call read_from_job_options(argbuff,mosaic,transient_run,
     1             trans_startyr,ctemloop,ctem_on,ncyear,lnduseon,
     2             spinfast,cyclemet,nummetcylyrs,metcylyrst,co2on,
     3             setco2conc,popdon,popcycleyr,parallelrun,dofire,
     4             dowetlands,obswetf,compete,inibioclim,start_bare,
     5             rsfile,start_from_rs,jmosty,idisp,izref,islfd,ipcp,
     6             itc,itcg,itg,iwf,ipai,ihgt,ialc,ials,ialg,isnoalb,
     7             igralb,jhhstd,jhhendd,jdstd,jdendd,jhhsty,jhhendy,
     8             jdsty,jdendy)

c     Initialize the CTEM parameters
      call initpftpars(compete)
c
c     set ictemmod, which is the class switch for coupling to ctem
c     either to 1 (ctem is coupled to class) or 0 (class runs alone)
c     this switch is set based on ctem_on that was set by read_from_job_options
c
      if (ctem_on) then
        ictemmod = 1
      else  !ctem_on is false
        ictemmod = 0
      end if
c
      lopcount = 1   ! initialize loop count to 1.
c
c     checking the time spent for running model
c
c      call idate(today)
c      call itime(now)
c      write(*,1000)   today(2), today(1), 2000+today(3), now
c 1000 format( 'start date: ', i2.2, '/', i2.2, '/', i4.4,
c     &      '; start time: ', i2.2, ':', i2.2, ':', i2.2 )
c
C     INITIALIZATION FOR COUPLING CLASS AND CTEM
C
       call initrowvars()
       call resetclassaccum(nltest,nmtest)

       IMONTH = 0
!       DO I=1,NLAT
!         DO M=1,NMOS
!           PREACC_M(I,M)=0.
!           GTACC_M(I,M)=0.
!           QEVPACC_M(I,M)=0.
!           HFSACC_M(I,M)=0.
!           HMFNACC_M(I,M)=0.
!           ROFACC_M(I,M)=0.
!           SNOACC_M(I,M)=0.
!           OVRACC_M(I,M)=0.
!           WTBLACC_M(I,M)=0.
!               DO J=1,IGND
!                 TBARACC_M(I,M,J)=0.
!                 THLQACC_M(I,M,J)=0.
!                 THICACC_M(I,M,J)=0.
!                 THALACC_M(I,M,J)=0.
!               ENDDO
!           ALVSACC_M(I,M)=0.
!           ALIRACC_M(I,M)=0.
!           RHOSACC_M(I,M)=0.
!           TSNOACC_M(I,M)=0.
!           WSNOACC_M(I,M)=0.
!           TCANACC_M(I,M)=0.
!           RCANACC_M(I,M)=0.
!           SCANACC_M(I,M)=0.
!           GROACC_M(I,M)=0.
!           FSINACC_M(I,M)=0.
!           FLINACC_M(I,M)=0.
!           TAACC_M(I,M)=0.
!           UVACC_M(I,M)=0.
!           PRESACC_M(I,M)=0.
!           QAACC_M(I,M)=0.
!           EVAPACC_M(I,M)=0.
!           FLUTACC_M(I,M)=0.
!         ENDDO
!       ENDDO
! C
       do 11 i=1,nlat
!        prbfrhucgrd(i)         = 0.0
!        extnprobgrd(i)         = 0.0
        barf(i)                = 1.0
! c
!        do j =1,12
!          mlightnggrd(i,j)=0.0
!        enddo
! c
        do 11 m=1,nmos
!         icountrow(i,m)           = 0
!         co2concrow(i,m)          = 0.0
!         npprow(i,m)              = 0.0
!         neprow(i,m)              = 0.0
!         hetroresrow(i,m)         = 0.0
!         autoresrow(i,m)          = 0.0
!         soilcresprow(i,m)         = 0.0
!         rmrow(i,m)               = 0.0
!         rgrow(i,m)               = 0.0
!         nbprow(i,m)              = 0.0
!         litresrow(i,m)           = 0.0
!         socresrow(i,m)           = 0.0
!         gpprow(i,m)              = 0.0
!         dstcemlsrow(i,m)         = 0.0
!         dstcemls3row(i,m)        = 0.0
!         litrfallrow(i,m)         = 0.0
!         humiftrsrow(i,m)         = 0.0
!         canresrow(i,m)           = 0.0
!         rmlrow(i,m)              = 0.0
!         rmsrow(i,m)              = 0.0
!         rmrrow(i,m)              = 0.0
!         lucemcomrow(i,m)         = 0.0
!         lucltrinrow(i,m)         = 0.0
!         lucsocinrow(i,m)         = 0.0
!         burnfracrow(i,m)         = 0.0
!         probfirerow(i,m)         = 0.0
!         btermrow(i,m)            = 0.0
!         ltermrow(i,m)            = 0.0
!         mtermrow(i,m)            = 0.0
!         cfluxcgrow(i,m)          = 0.0
!         cfluxcsrow(i,m)          = 0.0
! c
         TCANOACCROW_M(I,M)       = 0.0
         UVACCROW_M(I,M)          = 0.0
         VVACCROW_M(I,M)          = 0.0
         TCANOACCROW_OUT(I,M)     = 0.0
! c                                         !Rudra addes CH4 related variables on 03/12/2013
!         ch4wet1row(i,m)          = 0.0
!         ch4wet2row(i,m)          = 0.0
!         wetfdynrow(i,m)          = 0.0
!         ch4dyn1row(i,m)          = 0.0
!         ch4dyn2row(i,m)          = 0.0
!
! c
!         do j = 1, ignd
!            tbaraccrow_m(i,m,j)  = 0.0
!         enddo
! C
!         DO J = 1, ICAN
!           ZOLNCROW(I,M,J)        = 0.0
!           AILCROW(I,M,J)         = 0.0
!           CMASVEGCROW(I,M,J)     = 0.0
!           ALVSCTMROW(I,M,J)      = 0.0
!           ALIRCTMROW(I,M,J)      = 0.0
!           CSUM(I,M,J)            = 0.0
!           PAICROW(I,M,J)         = 0.0
!           SLAICROW(I,M,J)        = 0.0
!           DO K = 1, 3
!             RMATCROW(I,M,J,K)    = 0.0
!           ENDDO
!         ENDDO
!
!         do j = 1, icc
!           ailcgrow(i,m,j)        = 0.0
!           ailcgsrow(i,m,j)       = 0.0
!           fcancsrow(i,m,j)       = 0.0
!           fcancrow(i,m,j)        = 0.0
!           fcancmxrow(i,m,j)      = 0.0
!           co2i1cgrow(i,m,j)      = 0.0
!           co2i1csrow(i,m,j)      = 0.0
!           co2i2cgrow(i,m,j)      = 0.0
!           co2i2csrow(i,m,j)      = 0.0
!           ancsvegrow(i,m,j)      = 0.0
!           ancgvegrow(i,m,j)      = 0.0
!           rmlcsvegrow(i,m,j)     = 0.0
!           rmlcgvegrow(i,m,j)     = 0.0
!           stemmassrow(i,m,j)     = 0.0
!           rootmassrow(i,m,j)     = 0.0
!           ailcbrow(i,m,j)        = 0.0
!           grwtheffrow(i,m,j)     = 0.0
!           dvdfcanrow(i,m,j)      = 0.0
!           bmasvegrow(i,m,j)      = 0.0
!           tltrleafrow(i,m,j)     = 0.0
!           tltrstemrow(i,m,j)     = 0.0
!           tltrrootrow(i,m,j)     = 0.0
!           leaflitrrow(i,m,j)     = 0.0
!           roottemprow(i,m,j)     = 0.0
!           afrleafrow(i,m,j)      = 0.0
!           afrstemrow(i,m,j)      = 0.0
!           afrrootrow(i,m,j)      = 0.0
!           wtstatusrow(i,m,j)     = 0.0
!           ltstatusrow(i,m,j)     = 0.0
!           ailcminrow(i,m,j)      = 0.0
!           ailcmaxrow(i,m,j)      = 0.0
!           pfcancmxrow(i,m,j)     = 0.0
!           nfcancmxrow(i,m,j)     = 0.0
!           nppvegrow(i,m,j)       = 0.0
!           veghghtrow(i,m,j)      = 0.0
!           rootdpthrow(i,m,j)     = 0.0
!           gleafmasrow(i,m,j)     = 0.0
!           bleafmasrow(i,m,j)     = 0.0
!           anvegrow(i,m,j)        = 0.0
!           rmlvegrow(i,m,j)       = 0.0
! c
!           rmlvegaccrow(i,m,j)    = 0.0
!           rmsvegrow(i,m,j)       = 0.0
!           rmrvegrow(i,m,j)       = 0.0
!           rgvegrow(i,m,j)        = 0.0
! c
!           vgbiomas_vegrow(i,m,j) = 0.0
! c
!           gppvegrow(i,m,j) = 0.0
!           autoresvegrow(i,m,j) = 0.0
!
!           emit_co2row(i,m,j)         =0.0
!           emit_corow(i,m,j)          =0.0
!           emit_ch4row(i,m,j)         =0.0
!           emit_nmhcrow(i,m,j)        =0.0
!           emit_h2row(i,m,j)          =0.0
!           emit_noxrow(i,m,j)         =0.0
!           emit_n2orow(i,m,j)         =0.0
!           emit_pm25row(i,m,j)        =0.0
!           emit_tpmrow(i,m,j)         =0.0
!           emit_tcrow(i,m,j)          =0.0
!           emit_ocrow(i,m,j)          =0.0
!           emit_bcrow(i,m,j)          =0.0
!           burnvegfrow(i,m,j)         =0.0
!
! c
!           do k = 1, ignd
!             rmatctemrow(i,m,j,k) = 0.0
!           enddo
!         enddo
!
!         do j = 1, iccp1
!           litrmassrow(i,m,j)    = 0.0
!           soilcmasrow(i,m,j)    = 0.0
!           hetroresvegrow(i,m,j) = 0.0
!           litresvegrow(i,m,j) = 0.0
!           soilcresvegrow(i,m,j) = 0.0
!           nepvegrow(i,m,j) = 0.0
!           nbpvegrow(i,m,j) = 0.0
!
!         enddo
!
 11     continue
c
!     ==================================================================================================

c     do some initializations for the reading in of data from files. these
c     initializations primarily affect how the model does a spinup or transient
c     simulation and which years of the input data are being read.

      if (.not. cyclemet .and. transient_run) then !transient simulation, set to dummy values
        metcylyrst=-9999
        metcycendyr=9999
      else
c       find the final year of the cycling met
c       metcylyrst is defined in the joboptions file
        metcycendyr = metcylyrst + nummetcylyrs - 1
      endif

c     if cycling met (and not doing a transient run), find the popd and luc year to cycle with.
c     it is assumed that you always want to cycle the popd and luc
c     on the same year to be consistent. so if you are cycling the
c     met data, you can set a popd year (first case), or if cycling
c     the met data you can let the popcycleyr default to the met cycling
c     year by setting popcycleyr to -9999 (second case). if not cycling
c     the met data or you are doing a transient run that cycles the MET
c     at the start, cypopyr and cylucyr will default to a dummy value
c     (last case). (See example at bottom of read_from_job_options.f90
c     if confused)
c
      if (cyclemet .and. popcycleyr .ne. -9999 .and.
     &                                .not. transient_run) then
        cypopyr = popcycleyr
        cylucyr = popcycleyr
!        cywetldyr = popcycleyr
      else if (cyclemet .and. .not. transient_run) then
        cypopyr = metcylyrst
        cylucyr = metcylyrst
!        cywetldyr = metcylyrst
      else  ! give dummy value
        cypopyr = -9999
        cylucyr = -9999
!        cywetldyr = -9999
      end if

c     ctem initialization done
c
c     open files for reading and writing.
c     these are for coupled model (class_ctem)
c     we added both grid and mosaic output files
c
c     * input files

c         If we wish to restart from the .CTM_RS and .INI_RS files, then
c         we move the original RS files into place and start from them.
          if (start_from_rs) then
             command='mv '//argbuff(1:strlen(argbuff))//'.INI_RS '
     &                    //argbuff(1:strlen(argbuff))//'.INI'
             call system(command)
             command='mv '//argbuff(1:strlen(argbuff))//'.CTM_RS '
     &                    //argbuff(1:strlen(argbuff))//'.CTM'
             call system(command)
          end if


        open(unit=10,file=argbuff(1:strlen(argbuff))//'.INI',
     &       status='old')

!       if (ctem_on) then
!       open(unit=11,file=argbuff(1:strlen(argbuff))//'.CTM',
!     &       status='old')
!       endif
        open(unit=12,file=argbuff(1:strlen(argbuff))//'.MET',
     &      status='old')

c     luc file is opened in initialize_luc subroutine

      if (popdon) then
        open(unit=13,file=argbuff(1:strlen(argbuff))//'.POPD',
     &       status='old')
        read(13,*)  !Skip 3 lines of header
        read(13,*)
        read(13,*)
      endif
      if (co2on) then
        open(unit=14,file=argbuff(1:strlen(argbuff))//'.CO2',
     &         status='old')
      endif

c
      if (obswetf) then
        open(unit=16,file=argbuff(1:strlen(argbuff))//'.WET',
     &         status='old')
      endif
c
c     * CLASS output files
c
      if (.not. parallelrun) then ! stand alone mode, includes half-hourly and daily output
       OPEN(UNIT=61,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF1_G')  ! GRID-LEVEL DAILY OUTPUT FROM CLASS
       OPEN(UNIT=62,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF2_G')
       OPEN(UNIT=63,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF3_G')

       OPEN(UNIT=611,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF1_M') ! MOSAIC DAILY OUTPUT FROM CLASS
       OPEN(UNIT=621,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF2_M')
       OPEN(UNIT=631,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF3_M')

       OPEN(UNIT=64,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF4_M')  ! MOSAIC HALF-HOURLY OUTPUT FROM CLASS
       OPEN(UNIT=65,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF5_M')
       OPEN(UNIT=66,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF6_M')
       OPEN(UNIT=67,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF7_M')
       OPEN(UNIT=68,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF8_M')
       OPEN(UNIT=69,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF9_M')

       OPEN(UNIT=641,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF4_G') ! GRID-LEVEL HALF-HOURLY OUTPUT FROM CLASS
       OPEN(UNIT=651,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF5_G')
       OPEN(UNIT=661,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF6_G')
       OPEN(UNIT=671,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF7_G')
       OPEN(UNIT=681,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF8_G')
       OPEN(UNIT=691,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF9_G')
       end if
C
c      the ctem output file suffix naming convention is as follows:
c                       ".CT##{time}_{mosaic/grid}"
c      where the ## is a numerical identifier, {time} is any of H, D, M,
c      or Y for half hourly, daily, monthly, or yearly, respectively.
c      after the underscore M or G is used to denote mosaic or grid
c      -averaged values, respectively. also possible is GM for competition
c      outputs since they are the same format in either composite or
c      mosaic modes.
c
       ! Set up the CTEM half-hourly, daily, monthly and yearly files (if all needed), also
       ! setup the CLASS monthly and annual output files:
! !       JRM: comment below:
!        if (ctem_on) then
! c       ctem half hourly output files
!         open(unit=71, file=argbuff(1:strlen(argbuff))//'.CT01H_M')
!         open(unit=711,file=argbuff(1:strlen(argbuff))//'.CT01H_G')
!
!         if (mosaic) then
! c        ctem daily output files (mosaic)
!          open(unit=72,file=argbuff(1:strlen(argbuff))//'.CT01D_M')
!          open(unit=73,file=argbuff(1:strlen(argbuff))//'.CT02D_M')
!          open(unit=74,file=argbuff(1:strlen(argbuff))//'.CT03D_M')
!          open(unit=75,file=argbuff(1:strlen(argbuff))//'.CT04D_M')
!          open(unit=76,file=argbuff(1:strlen(argbuff))//'.CT05D_M')
!
!          if (dofire .or. lnduseon) then
!           open(unit=78,file=argbuff(1:strlen(argbuff))//'.CT06D_M') ! disturbance vars
!          endif
!
!         end if ! mosaic
!
! c        ctem daily output files (grid-average)
!          open(unit=721,file=argbuff(1:strlen(argbuff))//'.CT01D_G')
!          open(unit=731,file=argbuff(1:strlen(argbuff))//'.CT02D_G')
!          open(unit=741,file=argbuff(1:strlen(argbuff))//'.CT03D_G')
!          open(unit=751,file=argbuff(1:strlen(argbuff))//'.CT04D_G')
! c
!          if (dofire .or. lnduseon) then
!           open(unit=781,file=argbuff(1:strlen(argbuff))//'.CT06D_G') ! disturbance vars
!          endif
!
!          if (compete .or. lnduseon) then
!           open(unit=761,file=argbuff(1:strlen(argbuff))//'.CT07D_G') ! competition
!          end if
! c
!           if (dowetlands .or. obswetf) then
!           open(unit=762,file=argbuff(1:strlen(argbuff))//'.CT08D_G') ! Methane(Wetland)
!           endif
! c
!        endif ! ctem_on
!        endif ! parallelrun JRM: Keep this line uncommented!
!
! c     monthly & yearly output for both parallel mode and stand alone mode
!
! C     CLASS MONTHLY OUTPUT FILES
!       OPEN(UNIT=81,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF1M_G')
!       OPEN(UNIT=82,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF2M_G')
!
! C     CLASS YEARLY OUTPUT FILES
!       OPEN(UNIT=83,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF1Y_G')
!
!       if (ctem_on) then
!
!        if (mosaic) then
!
! c       Mosaic file names:
!
! c       CTEM monthly output files
!         open(unit=84,file=argbuff(1:strlen(argbuff))//'.CT01M_M')
!
! c       CTEM yearly output files
!         open(unit=86,file=argbuff(1:strlen(argbuff))//'.CT01Y_M')
!
!         if (dofire .or. lnduseon) then
!          open(unit=85,file=argbuff(1:strlen(argbuff))//'.CT06M_M') ! Monthly disturbance
!          open(unit=87,file=argbuff(1:strlen(argbuff))//'.CT06Y_M') ! Annual disturbance
!         endif
! c
!        else
!
! c       Composite file names:
!
! c       CTEM monthly output files
!         open(unit=84,file=argbuff(1:strlen(argbuff))//'.CT01M_G')
!
! c       CTEM yearly output files
!         open(unit=86,file=argbuff(1:strlen(argbuff))//'.CT01Y_G')
!
!         if (dofire .or. lnduseon) then
!          open(unit=85,file=argbuff(1:strlen(argbuff))//'.CT06M_G') ! Monthly disturbance
!          open(unit=87,file=argbuff(1:strlen(argbuff))//'.CT06Y_G') ! Annual disturbance
!         endif
!
!       end if !mosiac/composite
!
!         if (compete .or. lnduseon) then
!          open(unit=88,file=argbuff(1:strlen(argbuff))//'.CT07M_GM')! ctem pft fractions MONTHLY
!
!          open(unit=89,file=argbuff(1:strlen(argbuff))//'.CT07Y_GM')! ctem pft fractions YEARLY
!         endif
! c
!        if (dowetlands .or. obswetf) then
!         open(unit=91,file=argbuff(1:strlen(argbuff))//'.CT08M_G')  !Methane(wetland) MONTHLY
! c
!         open(unit=92,file=argbuff(1:strlen(argbuff))//'.CT08Y_G')  !Methane(wetland) YEARLY
!        endif !dowetlands
!       end if !ctem_on
! c
! !        JRM to here!
! C=======================================================================
! C
C     * READ AND PROCESS INITIALIZATION AND BACKGROUND INFORMATION.
C     * FIRST, MODEL RUN SPECIFICATIONS.

      READ (10,5010) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
      READ (10,5010) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
      READ (10,5010) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
C

       ! Set up the CTEM half-hourly, daily, monthly and yearly files (if all needed), also
       ! setup the CLASS monthly and annual output files:

       call create_outfiles(argbuff,title1, title2, title3, title4,
     1                     title5,title6,name1, name2, name3, name4,
     2                     name5, name6, place1,place2, place3,
     3                     place4, place5, place6)

      IF(CTEM_ON) THEN

!        READ (11,7010) TITLEC1
!        READ (11,7010) TITLEC2
!        READ (11,7010) TITLEC3
!
        if(obswetf) then
         read(16,*) TITLEC1
        end if
       ENDIF
C
      IF (.NOT. PARALLELRUN) THEN ! STAND ALONE MODE, INCLUDES HALF-HOURLY AND DAILY OUTPUT
C
       WRITE(61,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(61,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
       WRITE(61,6011)
6011  FORMAT(2X,'DAY  YEAR  K*  L*  QH  QE  SM  QG  ',
     1          'TR  SWE  DS  WS  AL  ROF  CUMS')

       WRITE(62,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(62,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6

       IF(IGND.GT.3) THEN
          WRITE(62,6012)
6012      FORMAT(2X,'DAY  YEAR  TG1  THL1  THI1  TG2  THL2  THI2  ',
     1              'TG3  THL3  THI3  TG4  THL4  THI4  TG5  THL5  ',
     2              'THI5')

       ELSE
          WRITE(62,6212)
6212      FORMAT(2X,'DAY  YEAR  TG1  THL1  THI1  TG2  THL2  THI2  ',
     1              'TG3  THL3  THI3  TCN  RCAN  SCAN  TSN  ZSN')

       ENDIF

       WRITE(63,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(63,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6

       IF(IGND.GT.3) THEN
          WRITE(63,6013)
6013      FORMAT(2X,'DAY  YEAR  TG6  THL6  THI6  TG7  THL7  THI7  ',
     1              'TG8  THL8  THI8  TG9  THL9  THI9  TG10'  ,
     2              'THL10  THI10')

       ELSE
          WRITE(63,6313)
6313      FORMAT(2X,'DAY YEAR KIN LIN TA UV PRES QA PCP EVAP')

       ENDIF
C
       WRITE(64,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(64,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
       WRITE(64,6014)
6014  FORMAT(2X,'HOUR  MIN  DAY  YEAR  K*  L*  QH  QE  SM  QG  ',
     1          'TR  SWE  DS  WS  AL  ROF  TPN  ZPN  CDH  CDM  ',
     2          'SFCU  SFCV  UV')

       WRITE(65,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(65,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6

       IF(IGND.GT.3) THEN
          WRITE(65,6015)
6015      FORMAT(2X,'HOUR  MIN  DAY  YEAR  TG1  THL1  THI1  TG2  ',
     1          'THL2  THI2  TG3  THL3  THI3  TG4  THL4  THI4  ',
     2          'TG5  THL5  THI5')

       ELSE
          WRITE(65,6515)
6515      FORMAT(2X,'HOUR  MIN  DAY  YEAR  TG1  THL1  THI1  TG2  ',
     1           'THL2  THI2  TG3  THL3  THI3  TCN  RCAN  SCAN  ',
     2           'TSN  ZSN  TCN-TA  TCANO  TAC  ACTLYR  FTABLE')

       ENDIF

       WRITE(66,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(66,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6

       IF(IGND.GT.3) THEN
          WRITE(66,6016)
6016      FORMAT(2X,'HOUR  MIN  DAY  YEAR  TG6  THL6  THI6  TG7  ',
     1          'THL7  THI7  TG8  THL8  THI8  TG9  THL9  THI9  ',
     2          'TG10  THL10  THI10  G0  G1  G2  G3  G4  G5  G6  ',
     3          'G7  G8  G9')

       ELSE
          WRITE(66,6616)
          WRITE(66,6615)
6616  FORMAT(2X,'HOUR  MIN  DAY  SWIN  LWIN  PCP  TA  VA  PA  QA')
6615  FORMAT(2X,'IF IGND <= 3, THIS FILE IS EMPTY')
       ENDIF

       WRITE(67,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(67,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
       WRITE(67,6017)
!     6017  FORMAT(2X,'WCAN SCAN CWLCAP CWFCAP FC FG FCS FGS CDH ', !runclass formatted.
!     1          'TCANO TCANS ALBS')
6017  FORMAT(2X,'HOUR  MIN  DAY  YEAR  ',
     1  'TROF     TROO     TROS     TROB      ROF     ROFO   ',
     2  '  ROFS        ROFB         FCS        FGS        FC       FG')

       WRITE(68,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(68,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
       WRITE(68,6018)
6018  FORMAT(2X,'HOUR  MIN  DAY  YEAR  ',
     1          'FSGV FSGS FSGG FLGV FLGS FLGG HFSC HFSS HFSG ',
     2          'HEVC HEVS HEVG HMFC HMFS HMFG1 HMFG2 HMFG3 ',
     3          'HTCC HTCS HTC1 HTC2 HTC3')

       WRITE(69,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(69,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
       WRITE(69,6019)
6019  FORMAT(2X,'HOUR  MIN  DAY  YEAR  ',
     1   'PCFC PCLC PCPN PCPG QFCF QFCL QFN QFG QFC1 ',
     2          'QFC2 QFC3 ROFC ROFN ROFO ROF WTRC WTRS WTRG')
!       runclass also has: EVDF ','CTV CTS CT1 CT2 CT3')
C
       WRITE(611,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(611,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
       WRITE(611,6011)
       WRITE(621,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(621,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
C
       IF(IGND.GT.3) THEN
           WRITE(621,6012)
       ELSE
           WRITE(621,6212)
       ENDIF
C
       WRITE(631,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(631,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
C
       IF(IGND.GT.3) THEN
          WRITE(631,6013)
       ELSE
          WRITE(631,6313)
       ENDIF
C
       WRITE(641,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(641,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
       WRITE(641,6008)
       WRITE(651,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(651,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
C
       IF(IGND.GT.3) THEN
          WRITE(651,6015)
       ELSE
          WRITE(651,6515)
       ENDIF
C
       WRITE(661,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(661,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
C
       IF(IGND.GT.3) THEN
          WRITE(661,6016)
       ELSE
          WRITE(661,6616)
       ENDIF
C
       WRITE(671,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(671,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
       WRITE(671,6017)
       WRITE(681,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(681,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
       WRITE(681,6018)
       WRITE(691,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(691,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
       WRITE(691,6019)
C
6008  FORMAT(2X,'HOUR  MIN  DAY  YEAR  K*  L*  QH  QE  SM  QG  ',
     1          'TR  SWE  DS  WS  AL  ROF  TPN  ZPN  CDH  CDM  ',
     2          'SFCU  SFCV  UV')

C
! !          JRM comment out below:
! C     CTEM FILE TITLES
! C
!       IF (CTEM_ON) THEN
!         WRITE(71,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
!         WRITE(71,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
!         WRITE(71,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
!         WRITE(71,7020)
!         WRITE(71,7030)
! C
!        IF (MOSAIC) THEN
!         WRITE(72,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
!         WRITE(72,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
!         WRITE(72,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
!         WRITE(72,7020)
!         WRITE(72,7040)
! C
!         WRITE(73,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
!         WRITE(73,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
!         WRITE(73,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
!         WRITE(73,7020)
!         WRITE(73,7050)
! C
!         WRITE(74,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
!         WRITE(74,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
!         WRITE(74,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
!         WRITE(74,7020)
!         WRITE(74,7061)
! C
!         WRITE(75,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
!         WRITE(75,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
!         WRITE(75,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
!         WRITE(75,7020)
!         WRITE(75,7070)
! C
!         WRITE(76,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
!         WRITE(76,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
!         WRITE(76,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
!         WRITE(76,7020)
!         WRITE(76,7080)
! C
!        IF (DOFIRE .OR. LNDUSEON) THEN
!         WRITE(78,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
!         WRITE(78,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
!         WRITE(78,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
!         WRITE(78,7021)
!         WRITE(78,7110)
!         WRITE(78,7111)
!        ENDIF
!       END IF !mosaic
! C
! 7010  FORMAT(A80)
! 7020  FORMAT('CANADIAN TERRESTRIAL ECOSYSTEM MODEL (CTEM) DAILY RESULTS'
!      &)
! 7021  FORMAT('CANADIAN TERRESTRIAL ECOSYSTEM MODEL (CTEM) DAILY ',
!      &' DISTURBANCE RESULTS')
! 7030  FORMAT('HOUR MIN  DAY, An FOR 9 PFTs, RmL FOR 9 PFTs')
! 7040  FORMAT('  DAY YEAR       GPP       NPP       NEP       NBP',
!      &'   AUTORES  HETRORES    LITRES    SOCRES  DSTCEMLS  LITRFALL',
!      &'  HUMIFTRS')
! 7050  FORMAT('  DAY YEAR       RML       RMS       RMR        RG',
!      &'  LEAFLITR  TLTRLEAF  TLTRSTEM  TLTRROOT ')
! 7060  FORMAT('  DAY YEAR  VGBIOMAS   GAVGLAI  GAVGLTMS  GAVGSCMS  ',
!      &'TOTCMASS  GLEAFMAS   BLEAFMAS STEMMASS   ROOTMASS  LITRMASS ',
!      &' SOILCMAS')
! 7061  FORMAT('  DAY YEAR  VGBIOMAS   GAVGLAI  GLEAFMAS   BLEAFMAS ',
!      & 'STEMMASS   ROOTMASS  LITRMASS SOILCMAS')
! 7070  FORMAT('  DAY YEAR     AILCG     AILCB    RMATCTEM ',
!      &'LAYER 1,2, & 3     VEGHGHT  ROOTDPTH  ROOTTEMP      SLAI')
! 7075  FORMAT('  DAY YEAR   FRAC #1   FRAC #2   FRAC #3   FRAC #4   ',
!      &'FRAC #5   FRAC #6   FRAC #7   FRAC #8   FRAC #9  ',
!      &'FRAC #10[%] SUMCHECK')
! 7080  FORMAT('  DAY YEAR   AFRLEAF   AFRSTEM   AFRROOT  TCANOACC',
!      &'  LFSTATUS')
! 7110  FORMAT('  DAY YEAR   EMIT_CO2',
!      &'    EMIT_CO   EMIT_CH4  EMIT_NMHC    EMIT_H2   EMIT_NOX',
!      &'   EMIT_N2O  EMIT_PM25   EMIT_TPM    EMIT_TC    EMIT_OC',
!      &'    EMIT_BC   BURNFRAC   PROBFIRE   LUCEMCOM   LUCLTRIN',
!      &'   LUCSOCIN   GRCLAREA   BTERM   LTERM   MTERM')
! 7111  FORMAT('               g/m2.D     g/m2.d',
!      &'     g/m2.d     g/m2.d     g/m2.d     g/m2.d     g/m2.d',
!      &'     g/m2.d     g/m2.d     g/m2.d     g/m2.d     g/m2.d   ',
!      &'       %  avgprob/d uMOL-CO2/M2.S KgC/M2.D',
!      &'   KgC/M2.D      KM^2    prob/d       prob/d       prob/d')
! 7112  FORMAT(' DAY  YEAR   CH4WET1    CH4WET2    WETFDYN   CH4DYN1
!      & CH4DYN2 ')
! 7113  FORMAT('          umolCH4/M2.S    umolCH4/M2.S          umolCH4/M2.S
!      & umolCH4/M2.S')
!
! C
!         WRITE(711,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
!         WRITE(711,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
!         WRITE(711,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
!         WRITE(711,7020)
!         WRITE(711,7030)
! C
!         WRITE(721,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
!         WRITE(721,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
!         WRITE(721,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
!         WRITE(721,7020)
!         WRITE(721,7040)
! C
!         WRITE(731,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
!         WRITE(731,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
!         WRITE(731,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
!         WRITE(731,7020)
!         WRITE(731,7050)
! C
!         WRITE(741,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
!         WRITE(741,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
!         WRITE(741,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
!         WRITE(741,7020)
!         WRITE(741,7060)
! C
!         WRITE(751,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
!         WRITE(751,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
!         WRITE(751,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
!         WRITE(751,7020)
!         WRITE(751,7070)
! C
!         IF (COMPETE .OR. LNDUSEON) THEN
!          WRITE(761,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
!          WRITE(761,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
!          WRITE(761,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
!          WRITE(761,7020)
!          WRITE(761,7075)
!         ENDIF
! C
!        IF (DOFIRE .OR. LNDUSEON) THEN
!         WRITE(781,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
!         WRITE(781,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
!         WRITE(781,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
!         WRITE(781,7020)
!         WRITE(781,7110)
!         WRITE(781,7111)
!        ENDIF
! c             Methane(Wetland) variables !Rudra
!        IF (DOWETLANDS .OR. OBSWETF) THEN
!         WRITE(762,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
!         WRITE(762,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
!         WRITE(762,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
!         WRITE(762,7020)
!         WRITE(762,7112)
!         WRITE(762,7113)
!        ENDIF
!
!       ENDIF !CTEM_ON
C
       ENDIF !IF NOT PARALLELRUN  ! JRM keep this line uncommented

! C     MONTHLY & YEARLY OUTPUT FOR BOTH PARALLEL MODE AND STAND ALONE MODE
!       WRITE(81,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
!       WRITE(81,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
!       WRITE(81,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
!       WRITE(81,6021)
!       WRITE(81,6121)
! C
!       WRITE(82,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
!       WRITE(82,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
!       WRITE(82,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
!       WRITE(82,6022)
!       WRITE(82,6122)
! C
!       WRITE(83,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
!       WRITE(83,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
!       WRITE(83,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
!       WRITE(83,6023)
!       WRITE(83,6123)
! C
!       IF (CTEM_ON) THEN
!          WRITE(84,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
!          WRITE(84,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
!          WRITE(84,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
!          WRITE(84,6024)
!          WRITE(84,6124)
!          WRITE(84,6224)
! C
!        IF (DOFIRE .OR. LNDUSEON) THEN
!         WRITE(85,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
!         WRITE(85,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
!         WRITE(85,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
!         WRITE(85,6025)
!         WRITE(85,6125)
!         WRITE(85,6225)
!        ENDIF
! C
!          WRITE(86,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
!          WRITE(86,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
!          WRITE(86,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
!          WRITE(86,6026)
!          WRITE(86,6126)
!          WRITE(86,6226)
! C
!         IF (DOFIRE .OR. LNDUSEON) THEN
!          WRITE(87,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
!          WRITE(87,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
!          WRITE(87,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
!          WRITE(87,6027)
!          WRITE(87,6127)
!          WRITE(87,6227)
!         ENDIF
! C
!         IF (COMPETE .OR. LNDUSEON) THEN
!           WRITE(88,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
!           WRITE(88,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
!           WRITE(88,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
!           WRITE(88,6028)
!           WRITE(88,6128)
!           WRITE(88,6228)
! C
!           WRITE(89,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
!           WRITE(89,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
!           WRITE(89,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
!           WRITE(89,6029)
!           WRITE(89,6129)
!           WRITE(89,6229)
!         ENDIF !COMPETE
!
!          IF (DOWETLANDS .OR. OBSWETF) THEN
!           WRITE(91,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
!           WRITE(91,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
!           WRITE(91,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
!           WRITE(91,6028)
!           WRITE(91,6230)
!           WRITE(91,6231)
!
!           WRITE(92,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
!           WRITE(92,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
!           WRITE(92,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
!           WRITE(92,6028)
!           WRITE(92,6232)
!           WRITE(92,6233)
!          ENDIF
!
!
! C
!       ENDIF !CTEM_ON
! C
! 6021  FORMAT(2X,'MONTH YEAR  SW     LW      QH      QE    SNOACC    ',
!      &        'WSNOACC    ROFACC      PCP      EVAP       TAIR')
! 6121  FORMAT(2X,'           W/m2    W/m2    W/m2    W/m2    kg/m2   ',
!      &        'kg/m2      mm.mon    mm.mon    mm.mon      degC')
! 6022  FORMAT(2X,'MONTH  YEAR  TG1  THL1  THI1     TG2  THL2  THI2',
!      &      '     TG3  THL3  THI3')
! 6122  FORMAT(2X,'             deg  m3/m3  m3/m3   deg  m3/m3  ',
!      &      'm3/m3   deg  m3/m3  m3/m3')
! 6023  FORMAT(2X,'YEAR   SW     LW      QH      QE     ROFACC   ',
!      &     ' PCP     EVAP  ')
! 6123  FORMAT(2X,'      W/m2   W/m2    W/m2    W/m2    mm.yr    ',
!      &     'mm.yr    mm.yr')
! 6024  FORMAT('CANADIAN TERRESTRIAL ECOSYSTEM MODEL (CTEM) MONTHLY ',
!      &'RESULTS')
! 6124  FORMAT('  MONTH  YEAR  LAIMAXG  VGBIOMAS  LITTER    SOIL_C  ',
!      &'  NPP       GPP        NEP       NBP    HETRES',
!      &'   AUTORES    LITRES   SOILCRES')
! 6224  FORMAT('                 m2/m2  Kg C/m2  Kg C/m2   Kg C/m2  ',
!      &       'gC/m2.mon  gC/m2.mon  gC/m2.mon  g/m2.mon   g/m2.mon ',
!      &       'gC/m2.mon  gC/m2.mon  gC/m2.mon')
! 6025  FORMAT('CANADIAN TERRESTRIAL ECOSYSTEM MODEL (CTEM) MONTHLY ',
!      &'RESULTS FOR DISTURBANCES')
! 6125  FORMAT('  MONTH  YEAR  CO2',
!      &'        CO        CH4      NMHC       H2       NOX       N2O',
!      &'       PM25       TPM        TC        OC        BC  ',
!      &' PROBFIRE  LUC_CO2_E  LUC_LTRIN  LUC_SOCIN   BURNFRAC    BTERM',
!      &' LTERM   MTERM')
! 6225  FORMAT('            g/m2.mon  g/m2.mon',
!      &'  g/m2.mon  g/m2.mon  g/m2.mon  g/m2.mon  g/m2.mon',
!      &'  g/m2.mon  g/m2.mon  g/m2.mon  g/m2.mon  g/m2.mon',
!      &'  prob/mon    g C/m2    g C/m2    g C/m2         %  prob/mon',
!      &'  prob/mon  prob/mon')
! 6026  FORMAT('CANADIAN TERRESTRIAL ECOSYSTEM MODEL (CTEM) YEARLY ',
!      &'RESULTS')
! 6126  FORMAT('  YEAR   LAIMAXG  VGBIOMAS  STEMMASS  ROOTMASS  LITRMASS',
!      &'  SOILCMAS  TOTCMASS  ANNUALNPP ANNUALGPP ANNUALNEP ANNUALNBP',
!      &' ANNHETRSP ANAUTORSP ANNLITRES ANSOILCRES')
! 6226  FORMAT('          m2/m2   Kg C/m2   Kg C/m2   Kg C/m2    Kg C/m2',
!      &'  Kg C/m2   Kg C/m2   gC/m2.yr  gC/m2.yr  gC/m2.yr  gC/m2.yr',
!      &'  gC/m2.yr  gC/m2.yr  gC/m2.yr  gC/m2.yr')
! 6027  FORMAT('CANADIAN TERRESTRIAL ECOSYSTEM MODEL (CTEM) YEARLY ',
!      &'RESULTS FOR DISTURBANCES')
!
! 6127  FORMAT('  YEAR   ANNUALCO2',
!      &'  ANNUALCO  ANNUALCH4  ANN_NMHC ANNUAL_H2 ANNUALNOX ANNUALN2O',
!      &'  ANN_PM25  ANNUALTPM ANNUAL_TC ANNUAL_OC ANNUAL_BC APROBFIRE',
!      &' ANNLUCCO2  ANNLUCLTR ANNLUCSOC ABURNFRAC ANNBTERM ANNLTERM',
!      &' ANNMTERM')
! 6227  FORMAT('         g/m2.yr',
!      &'  g/m2.yr  g/m2.yr  g/m2.yr  g/m2.yr  g/m2.yr  g/m2.yr',
!      &'  g/m2.yr  g/m2.yr  g/m2.yr  g/m2.yr  g/m2.yr  prob/yr ',
!      &'  g/m2.yr  g/m2.yr  g/m2.yr    %     prob/yr  prob/yr',
!      &'  prob/yr')
! 6028  FORMAT('CANADIAN TERRESTRIAL ECOSYSTEM MODEL (CTEM) MONTHLY ',
!      &'RESULTS')
! 6128  FORMAT(' MONTH YEAR  FRAC #1   FRAC #2   FRAC #3   FRAC #4   ',
!      &'FRAC #5   FRAC #6   FRAC #7   FRAC #8   FRAC #9   FRAC #10   ',
!      &'SUMCHECK, PFT existence for each of the 9 pfts')
! 6228  FORMAT('             %         %         %         %         ',
!      &'%         %         %         %         %         %          ',
!      &'     ')
! 6029  FORMAT('CANADIAN TERRESTRIAL ECOSYSTEM MODEL (CTEM) YEARLY ',
!      &'RESULTS')
! 6129  FORMAT('  YEAR   FRAC #1   FRAC #2   FRAC #3   FRAC #4   ',
!      &'FRAC #5   FRAC #6   FRAC #7   FRAC #8   FRAC #9   FRAC #10   ',
!      &'SUMCHECK, PFT existence for each of the 9 pfts')
! 6229  FORMAT('         %         %         %         %         ',
!      &'%         %         %         %         %         %          ',
!      &'     ')
!
! 6230  FORMAT('MONTH  YEAR   CH4WET1    CH4WET2    WETFDYN   CH4DYN1
!      & CH4DYN2 ')
! 6231  FORMAT('       gCH4/M2.MON     gCH4/M2.MON        gCH4/M2.MON
!      & gCH4/M2.MON')
! 6232  FORMAT('  YEAR   CH4WET1    CH4WET2    WETFDYN   CH4DYN1
!      & CH4DYN2 ')
! 6233  FORMAT('      gCH4/M2.YR      gCH4/M2.YR         gCH4/M2.YR
!      & gCH4/M2.YR ')
!
! C           JRM: To here!

C     CTEM FILE TITLES DONE
C======================= CTEM ========================================== /
C
C=======================================================================

C     BEGIN READ IN OF THE .INI FILE

      READ(10,5020) DLATROW(1),DEGLON,ZRFMROW(1),ZRFHROW(1),ZBLDROW(1),
     1              GCROW(1),NLTEST,NMTEST

      JLAT=NINT(DLATROW(1))
      RADJROW(1)=DLATROW(1)*PI/180.
      DLONROW(1)=DEGLON
      Z0ORROW(1)=0.0
      GGEOROW(1)=0.0
C     GGEOROW(1)=-0.035

      DO 50 I=1,NLTEST
      DO 50 M=1,NMTEST
          READ(10,5040) (FCANROT(I,M,J),J=1,ICAN+1),(PAMXROT(I,M,J),
     1                  J=1,ICAN)
          READ(10,5040) (LNZ0ROT(I,M,J),J=1,ICAN+1),(PAMNROT(I,M,J),
     1                  J=1,ICAN)
          READ(10,5040) (ALVCROT(I,M,J),J=1,ICAN+1),(CMASROT(I,M,J),
     1                  J=1,ICAN)
          READ(10,5040) (ALICROT(I,M,J),J=1,ICAN+1),(ROOTROT(I,M,J),
     1                  J=1,ICAN)
          READ(10,5030) (RSMNROT(I,M,J),J=1,ICAN),
     1                  (QA50ROT(I,M,J),J=1,ICAN)
          READ(10,5030) (VPDAROT(I,M,J),J=1,ICAN),
     1                  (VPDBROT(I,M,J),J=1,ICAN)
          READ(10,5030) (PSGAROT(I,M,J),J=1,ICAN),
     1                  (PSGBROT(I,M,J),J=1,ICAN)
          READ(10,5040) DRNROT(I,M),SDEPROT(I,M),FAREROT(I,M)
          READ(10,5090) XSLPROT(I,M),GRKFROT(I,M),WFSFROT(I,M),
     1                  WFCIROT(I,M),MIDROT(I,M)
          READ(10,5080) (SANDROT(I,M,J),J=1,3)
          READ(10,5080) (CLAYROT(I,M,J),J=1,3)
          READ(10,5080) (ORGMROT(I,M,J),J=1,3)
          READ(10,5050) (TBARROT(I,M,J),J=1,3),TCANROT(I,M),
     1                  TSNOROT(I,M),TPNDROT(I,M)
          READ(10,5060) (THLQROT(I,M,J),J=1,3),(THICROT(I,M,J),
     1                  J=1,3),ZPNDROT(I,M)
          READ(10,5070) RCANROT(I,M),SCANROT(I,M),SNOROT(I,M),
     1                  ALBSROT(I,M),RHOSROT(I,M),GROROT(I,M)

50    CONTINUE
C
      DO 25 J=1,IGND
          READ(10,*) DELZ(J),ZBOT(J) !was format 5002 but was causing problems so let be free-format. JM Jun 2015
 25   CONTINUE
 5002 FORMAT(2X,2F8.2)
C
C======================= CTEM ========================================== \

c     the output year ranges can be read in from the job options file or not.
c     if the values should be read in from the .ini file, and not
c     from the job options file, the job options file values are set to
c     -9999 thus triggering the read in of the .ini file values below
      if (jhhstd .eq. -9999) then
        read(10,5200) jhhstd,jhhendd,jdstd,jdendd
       read(10,5200) jhhsty,jhhendy,jdsty,jdendy
      end if
C======================= CTEM ========================================== /

      !JHHST=JHHSTY*1000+JHHSTD
      !JHHEND=JHHENDY*1000+JHHENDD
      !JDST=JDSTY*1000+JDSTD
      !JDEND=JDENDY*1000+JDENDD

      CLOSE(10)
C
C====================== CTEM =========================================== \
C
c     read from ctem initialization file (.CTM)

      if (ctem_on) then
      call read_from_ctm(nltest,nmtest,FCANROT,FAREROT,
     1                   RSMNROT,QA50ROT,VPDAROT,VPDBROT,PSGAROT,
     2                   PSGBROT,DRNROT,SDEPROT, XSLPROT,GRKFROT,
     3                   WFSFROT,WFCIROT,MIDROT,SANDROT, CLAYROT,
     4                   ORGMROT,TBARROT,THLQROT,THICROT,TCANROT,
     5                   TSNOROT,TPNDROT,ZPNDROT,RCANROT,SCANROT,
     6                   SNOROT, ALBSROT,RHOSROT,GROROT,argbuff)
      end if
c
!       if (ctem_on) then
!       do 71 i=1,nltest
!       do 72 m=1,nmtest
! c
! c         following three variables are needed to run ctem.
! c         min & max leaf area index are needed to break
! c         class lai into dcd and evg for trees (for crops and grasses it
! c         doesn't matter much).
! c
! c         dvdfcanrow is needed to divide needle & broad leaf into dcd and evg,
! c         and crops & grasses into c3 and c4 fractions.
! c
!           read(11,*) (ailcminrow(i,m,j),j=1,icc)
!           read(11,*) (ailcmaxrow(i,m,j),j=1,icc)
!           read(11,*) (dvdfcanrow(i,m,j),j=1,icc)
! c
! c         rest of the initialization variables are needed to run ctem.
! c         if starting from bare ground initialize all live and dead c pools from zero. suitable values
! c         of extnprobgrd and prbfrhucgrd would still be required. set stdalngrd to
! c         1 for operation in non-gcm stand alone mode, in the ctem
! c         initialization file.
! c
!           read(11,*) (gleafmasrow(i,m,j),j=1,icc)
!           read(11,*) (bleafmasrow(i,m,j),j=1,icc)
!           read(11,*) (stemmassrow(i,m,j),j=1,icc)
!           read(11,*) (rootmassrow(i,m,j),j=1,icc)
!
! c           If fire and competition are on, save the stemmass and rootmass for
! c           use in burntobare subroutine on the first timestep.
!             if (dofire .and. compete) then
!              do j =1,icc
!               pstemmassrow(i,m,j)=stemmassrow(i,m,j)
!               pgleafmassrow(i,m,j)=rootmassrow(i,m,j)
!              end do
!             end if
!
!           read(11,*) (litrmassrow(i,m,j),j=1,iccp1)
!           read(11,*) (soilcmasrow(i,m,j),j=1,iccp1)
!           read(11,*) (lfstatusrow(i,m,j),j=1,icc)
!           read(11,*) (pandaysrow(i,m,j),j=1,icc)
!
! 72      continue
!
!          read(11,*) (mlightnggrd(i,j),j=1,6)  !mean monthly lightning frequency
!          read(11,*) (mlightnggrd(i,j),j=7,12) !flashes/km2.year
!          read(11,*) extnprobgrd(i)
!          read(11,*) prbfrhucgrd(i)
!          read(11,*) stdalngrd(i)
!
!          if (compete .and. inibioclim) then  !read in the bioclimatic parameters
!           read(11,*) twarmm(i), tcoldm(i), gdd5(i), aridity(i),
!      1              srplsmon(i)
!           read(11,*) defctmon(i), anndefct(i), annsrpls(i),
!      1              annpcp(i), dry_season_length(i)
!          else if (compete .and. .not. inibioclim) then ! set them to zero
!            twarmm(i)=0.0
!            tcoldm(i)=0.0
!            gdd5(i)=0.0
!            aridity(i)=0.0
!            srplsmon(i)=0.0
!            defctmon(i)=0.0
!            anndefct(i)=0.0
!            annsrpls(i)=0.0
!            annpcp(i)=0.0
!            dry_season_length(i) = 0.0
!          endif
!
!          if (dowetlands) then      ! Rudra !if true then read wetland fractions
!              read(11,*) (wetfrac_sgrd(i,j),j=1,8)
!          endif
! 71    continue
!       close(11)
!       endif
!
! c     check that a competition or luc run has the correct number of mosaics
! c     if it is not a start_bare run, then nmtest should equal nmos
!       if (mosaic .and. (compete .or. lnduseon)
!      &                 .and. .not. start_bare) then
!         if (nmtest .ne. nmos) then
!            write(6,*)'compete or luc runs that do not start from bare'
!            write(6,*)'ground need the number of mosaics to equal icc+1'
!            write(6,*)'nmtest = ',nmtest,' nmos = ',nmos
!             call xit('runclass36ctem', -2)
!         endif
!       endif
! c
! c     if this run uses the competition or lnduseon parameterization and
! c     starts from bare ground, set up the model state here. this
! c     overwrites what was read in from the .ini and .ctm files.
! c     for composite runs (the composite set up is after this one for mosaics)
!       if ((compete .or. lnduseon) .and. start_bare) then
!
!        if (mosaic) then
!
! c       store the read-in crop fractions as we keep them even when we start bare.
! c       FLAG: this is setup assuming that crops are in mosaics 6 and 7. JM Apr 9 2014.
!          do i=1,nltest
!           crop_temp_frac(i,1)=FAREROT(i,6)
!           crop_temp_frac(i,2)=FAREROT(i,7)
!          end do
!
! c       check the number of mosaics that came from the .ini file
!         if (nmtest .ne. nmos) then
!
! c        we need to transfer some initial parameterization info to all
! c        mosaics, so set all values to that of the first mosaic.
!          do i=1,nltest
!           do m=nmtest+1,nmos
! c
!            do j=1,ican
!              RSMNROT(i,m,j)=RSMNROT(i,1,j)
!              QA50ROT(i,m,j)=QA50ROT(i,1,j)
!              VPDAROT(i,m,j)=VPDAROT(i,1,j)
!              VPDBROT(i,m,j)=VPDBROT(i,1,j)
!              PSGAROT(i,m,j)=PSGAROT(i,1,j)
!              PSGBROT(i,m,j)=PSGBROT(i,1,j)
!            enddo
! c
!            DRNROT(i,m)=DRNROT(i,1)
!            SDEPROT(i,m)=SDEPROT(i,1)
!            FAREROT(i,m)=FAREROT(i,1)
!            XSLPROT(i,m)=XSLPROT(i,1)
!            GRKFROT(i,m)=GRKFROT(i,1)
!            WFSFROT(i,m)=WFSFROT(i,1)
!            WFCIROT(i,m)=WFCIROT(i,1)
!            MIDROT(i,m)=MIDROT(i,1)
! c
!            do j=1,3
!             SANDROT(i,m,j)=SANDROT(i,1,j)
!             CLAYROT(i,m,j)=CLAYROT(i,1,j)
!             ORGMROT(i,m,j)=ORGMROT(i,1,j)
!             TBARROT(i,m,j)=TBARROT(i,1,j)
!             THLQROT(i,m,j)=THLQROT(i,1,j)
!             THICROT(i,m,j)=THICROT(i,1,j)
!            enddo
! c
!            TCANROT(i,m)=TCANROT(i,1)
!            TSNOROT(i,m)=TSNOROT(i,1)
!            TPNDROT(i,m)=TPNDROT(i,1)
!            ZPNDROT(i,m)=ZPNDROT(i,1)
!            RCANROT(i,m)=RCANROT(i,1)
!            SCANROT(i,m)=SCANROT(i,1)
!            SNOROT(i,m)=SNOROT(i,1)
!            ALBSROT(i,m)=ALBSROT(i,1)
!            RHOSROT(i,m)=RHOSROT(i,1)
!            GROROT(i,m)=GROROT(i,1)
!            do j=1,icc
!              lfstatusrow(i,m,j) = 4
!            enddo !j
! c
!           enddo !m
!          enddo !i
!
! c       set the number of mosaics to icc+1
!         nmtest=nmos
!
!         endif  !if (nmtest .ne. nmos)
! c
! c       set the initial conditions for the pfts
! c       (bah, this is such an inelegant way to do this, but oh well...)
! c
! c       initalize to zero
!         FCANROT=0.0
!         dvdfcanrow=0.0
!         FAREROT=0.0
!
!         do i=1,nltest
!          do m=1,nmtest
! c
! c         set the seed amount for each pft in its mosaic
!           if (compete .or. lnduseon) then
!             if (m .lt. icc+1) then
!              FAREROT(i,m)=seed
!             else
!              FAREROT(i,m)=1.0 - (real(icc) * seed)
!             endif
!           endif
!
!           do j = 1,icc
!             ailcminrow(i,m,j)=0.0
!             ailcmaxrow(i,m,j)=0.0
!             gleafmasrow(i,m,j)=0.0
!             bleafmasrow(i,m,j)=0.0
!             stemmassrow(i,m,j)=0.0
!             rootmassrow(i,m,j)=0.0
!             lfstatusrow(i,m,j)=4
!             pandaysrow(i,m,j)=0
!           enddo
!
!           lfstatusrow(i,m,1)=2
!
!           do j = 1,iccp1
!             litrmassrow(i,m,j)=0.
!             soilcmasrow(i,m,j)=0.
!           enddo
!
! c         initial conditions always required
!           dvdfcanrow(i,m,1)=1.0  !ndl
!           dvdfcanrow(i,m,3)=1.0  !bdl
!           dvdfcanrow(i,m,6)=1.0  !crop
!           dvdfcanrow(i,m,8)=1.0  !grasses
!
! c         then adjusted below for the actual mosaic makeup
!           if (m .le. 2) then                     !ndl
!            FCANROT(i,m,1)=1.0
!            if (m .eq. 2) then
!              dvdfcanrow(i,m,1)=0.0
!              dvdfcanrow(i,m,2)=1.0
!            endif
!           elseif (m .ge. 3 .and. m .le. 5) then  !bdl
!            FCANROT(i,m,2)=1.0
!            if (m .eq. 4) then
!              dvdfcanrow(i,m,3)=0.0
!              dvdfcanrow(i,m,4)=1.0
!            endif
!            if (m .eq. 5) then
!              dvdfcanrow(i,m,3)=0.0
!              dvdfcanrow(i,m,5)=1.0
!            endif
!           elseif (m .eq. 6 .or. m .eq. 7) then  !crop
!            FCANROT(i,m,3)=1.0
!            if (m .eq. 7) then
!              dvdfcanrow(i,m,6)=0.0
!              dvdfcanrow(i,m,7)=1.0
!            endif
!           elseif (m .eq. 8 .or. m .eq. 9) then  !grasses
!            FCANROT(i,m,4)=1.0
!            if (m .eq. 9) then
!              dvdfcanrow(i,m,8)=0.0
!              dvdfcanrow(i,m,9)=1.0
!            endif
!           else                                  !bare/urban?
!            FCANROT(i,m,5)=1.0
!            endif !mosaic adjustments
!          enddo  !m
!         enddo  !i
!
!
!          do i=1,nltest
!           FAREROT(i,6)=crop_temp_frac(i,1)
!           FAREROT(i,7)=crop_temp_frac(i,2)
!          end do
!
!       else if (.not. mosaic) then  !set up for composite runs when start_bare is on and compete or landuseon
!
! c       store the read-in crop fractions as we keep them even when we start bare.
! c       FLAG: this is setup assuming that crops are in pft number 6 and 7. JM Apr 9 2014.
!          do i=1,nltest
!           crop_temp_frac(i,1)=FCANROT(i,1,3)*dvdfcanrow(i,1,6)
!           crop_temp_frac(i,2)=FCANROT(i,1,3)*dvdfcanrow(i,1,7)
!          end do
!
! c      initalize to zero, these will be filled in by the luc or
! c      competition subroutines.
!        FCANROT=0.0
!        dvdfcanrow=0.0
!
!        ! Added this as start_bare runs were not properly assigning
!        ! a TCAN on the very first day since the FCANROT was 0. JM Jan 14 2014.
!        do i=1,nltest
!         do j=1,iccp1
!            if (j .lt. icc+1) then
!             FCANROT(i,1,j)=seed
!            else
!             FCANROT(i,1,j)=1.0 - (real(icc) * seed)
!            endif
!         end do
!        end do
!
!        do i=1,nltest
!
! c      initial conditions always required
!          dvdfcanrow(i,1,1)=1.0  !ndl
!          dvdfcanrow(i,1,3)=1.0  !bdl
!          dvdfcanrow(i,1,6)=1.0  !crop
!          dvdfcanrow(i,1,8)=1.0  !grasses
!
!          do j = 1,icc
!            ailcminrow(i,1,j)=0.0
!            ailcmaxrow(i,1,j)=0.0
!            gleafmasrow(i,1,j)=0.0
!            bleafmasrow(i,1,j)=0.0
!            stemmassrow(i,1,j)=0.0
!            rootmassrow(i,1,j)=0.0
!            lfstatusrow(i,1,j)=4
!            pandaysrow(i,1,j)=0
!          enddo
!
!          lfstatusrow(i,1,1)=2
!
!          do j = 1,iccp1
!            litrmassrow(i,1,j)=0.0
!            soilcmasrow(i,1,j)=0.0
!          enddo
!        enddo !nltest
!
!          do i=1,nltest
!           FCANROT(i,1,3) = crop_temp_frac(i,1) + crop_temp_frac(i,2)
!           if (FCANROT(i,1,3) .gt. abszero) then
!            dvdfcanrow(i,1,6) = crop_temp_frac(i,1) / FCANROT(i,1,3)
!            dvdfcanrow(i,1,7) = crop_temp_frac(i,2) / FCANROT(i,1,3)
!           else
!            dvdfcanrow(i,1,6) = 1.0
!            dvdfcanrow(i,1,7) = 0.0
!           end if
!          end do
!
!       end if ! mosaic / composite
!       end if !if (compete/landuseon .and. start_bare)

C===================== CTEM =============================================== /
C
      DO 100 I=1,NLTEST
      DO 100 M=1,NMTEST

          TBARROT(I,M,1)=TBARROT(I,M,1)+TFREZ
          TBARROT(I,M,2)=TBARROT(I,M,2)+TFREZ
          TBARROT(I,M,3)=TBARROT(I,M,3)+TFREZ
          TSNOROT(I,M)=TSNOROT(I,M)+TFREZ
          TCANROT(I,M)=TCANROT(I,M)+TFREZ

C         FLAG! TEMP FIX. ON RESTART A SMALL NUMBER OF CELLS (<10) WILL HAVE
C         A RELATIVELY LARGE DIFFERENCE BETWEEN THE CANOPY TEMP AND SNOW TEMP
C         THIS WILL CAUSE THE MODEL TO FAIL RIGHT AWAY. TO PREVENT THIS CHECK
C         IF THE CANOPY TEMP IS WITHIN 5 DEGREES OF THE SNOW TEMP (IF THERE IS
C         SNOW), AND IF SO THEN OVERWRITE THE TCAN WITH 1 DEGREE COLDER THAN THE SNOW TEMP.
C         JM FEB 5 2013
!          FLAG commented out for 3.6.2 JM Jun 2015.
!          IF (SNOROT(I,M) .GT. 0.0) then
!           IF ( ABS(TCANROT(I,M)-TSNOROT(I,M)) .GT. 5. ) THEN
!              TCANROT(I,M)=TSNOROT(I,M) - 0.5
!           ENDIF
!          ENDIF

          TPNDROT(I,M)=TPNDROT(I,M)+TFREZ
          TBASROT(I,M)=TBARROT(I,M,3)
          CMAIROT(I,M)=0.
          WSNOROT(I,M)=0.
          ZSNLROT(I,M)=0.10

C         THIS FIX BELOW IS TO CORRECT A BUG THAT CAUSES A CRASH DUE
C         TO UNREASONABLE CANOPY TEMPERATURES IN THE FIRST YEAR OF A RESTART
C         WITH SNOW ON THE GROUND. NOTE: RUNCLASS.f HAS THIS SAME PROBLEM. JM JAN 2013
!          FLAG commented out for 3.6.2 JM Jun 2015.
!          IF (SNOROT(I,M) .GT. 0.) THEN !THERE IS SNOW ON THE GROUND
!           TSFSROT(I,M,1)=TBARROT(I,M,1)
!           TSFSROT(I,M,2)=TBARROT(I,M,1)
!          ELSE ! NO SNOW SO JUST SET THESE TO FREEZING POINT
           TSFSROT(I,M,1)=TFREZ
           TSFSROT(I,M,2)=TFREZ
!          ENDIF

          TSFSROT(I,M,3)=TBARROT(I,M,1)
          TSFSROT(I,M,4)=TBARROT(I,M,1)
          TACROT (I,M)=TCANROT(I,M)
          QACROT (I,M)=0.5E-2

          IF(IGND.GT.3)                                 THEN
              DO 65 J=4,IGND,-1 ! FLAG: Yuanqiao said this caused problems. Check. JM Mar 31 2015.
                  TBARROT(I,M,J)=TBARROT(I,M,3)
                  IF(SDEPROT(I,M).LT.(ZBOT(J-1)+0.001) .AND.
     1                  SANDROT(I,M,3).GT.-2.5)     THEN
                      SANDROT(I,M,J)=-3.0
                      CLAYROT(I,M,J)=-3.0
                      ORGMROT(I,M,J)=-3.0
                      THLQROT(I,M,J)=0.0
                      THICROT(I,M,J)=0.0
                  ELSE
                      SANDROT(I,M,J)=SANDROT(I,M,3)
                      CLAYROT(I,M,J)=CLAYROT(I,M,3)
                      ORGMROT(I,M,J)=ORGMROT(I,M,3)
                      THLQROT(I,M,J)=THLQROT(I,M,3)
                      THICROT(I,M,J)=THICROT(I,M,3)
                  ENDIF
65            CONTINUE
          ENDIF

          DO 75 K=1,6
          DO 75 L=1,50
              ITCTROT(I,M,K,L)=0
75        CONTINUE
100   CONTINUE

      DO 150 I=1,NLTEST
          PREACC(I)=0.
          GTACC(I)=0.
          QEVPACC(I)=0.
          EVAPACC(I)=0.
          HFSACC(I)=0.
          HMFNACC(I)=0.
          ROFACC(I)=0.
          OVRACC(I)=0.
          WTBLACC(I)=0.
          ALVSACC(I)=0.
          ALIRACC(I)=0.
          RHOSACC(I)=0.
          SNOACC(I)=0.
          WSNOACC(I)=0.
          CANARE(I)=0.
          SNOARE(I)=0.
          TSNOACC(I)=0.
          TCANACC(I)=0.
          RCANACC(I)=0.
          SCANACC(I)=0.
          GROACC(I)=0.
          FSINACC(I)=0.
          FLINACC(I)=0.
          FLUTACC(I)=0.
          TAACC(I)=0.
          UVACC(I)=0.
          PRESACC(I)=0.
          QAACC(I)=0.
          DO 125 J=1,IGND
              TBARACC(I,J)=0.
              THLQACC(I,J)=0.
              THICACC(I,J)=0.
              THALACC(I,J)=0.
125       CONTINUE
150   CONTINUE
C
C===================== CTEM =============================================== \
c
c     initialize accumulated array for monthly & yearly output for class
c
      call resetclassmon(nltest)
      call resetclassyr(nltest)

!       do 151 i=1,nltest
!           ALVSACC_MO(I)=0.
!           ALIRACC_MO(I)=0.
!           FLUTACC_MO(I)=0.
!           FSINACC_MO(I)=0.
!           FLINACC_MO(I)=0.
!           HFSACC_MO(I) =0.
!           QEVPACC_MO(I)=0.
!           SNOACC_MO(I) =0.
!           WSNOACC_MO(I)=0.
!           ROFACC_MO(I) =0.
!           PREACC_MO(I) =0.
!           EVAPACC_MO(I)=0.
!           TAACC_MO(I)=0.
!
!           DO 152 J=1,IGND
!               TBARACC_MO(I,J)=0.
!               THLQACC_MO(I,J)=0.
!               THICACC_MO(I,J)=0.
! 152       CONTINUE
C

!           ALVSACC_YR(I)=0.
!           ALIRACC_YR(I)=0.
!           FLUTACC_YR(I)=0.
!           FSINACC_YR(I)=0.
!           FLINACC_YR(I)=0.
!           HFSACC_YR(I) =0.
!           QEVPACC_YR(I)=0.
!           ROFACC_YR(I) =0.
!           PREACC_YR(I) =0.
!           EVAPACC_YR(I)=0.
!           TAACC_YR(I)=0.
! C
! 151   CONTINUE
C===================== CTEM =============================================== /

      DO 175 I=1,200
          TAHIST(I)=0.0
          TCHIST(I)=0.0
          TACHIST(I)=0.0
          TDHIST(I)=0.0
          TD2HIST(I)=0.0
          TD3HIST(I)=0.0
          TD4HIST(I)=0.0
          TSHIST(I)=0.0
          TSCRHIST(I)=0.0
175   CONTINUE
      ALAVG=0.0
      ALMAX=0.0
      ACTLYR=0.0
      FTAVG=0.0
      FTMAX=0.0
      FTABLE=0.0

      CALL CLASSB(THPROT,THRROT,THMROT,BIROT,PSISROT,GRKSROT,
     1            THRAROT,HCPSROT,TCSROT,THFCROT,THLWROT,PSIWROT,
     2            DLZWROT,ZBTWROT,ALGWROT,ALGDROT,
     +            ALGWVROT,ALGWNROT,ALGDVROT,ALGDNROT,
     3            SANDROT,CLAYROT,ORGMROT,SOCIROT,DELZ,ZBOT,
     4            SDEPROT,ISNDROT,IGDRROT,
     5            NLAT,NMOS,1,NLTEST,NMTEST,IGND,IGRALB)

5010  FORMAT(2X,6A4)
5020  FORMAT(5F10.2,F7.1,3I5)
5030  FORMAT(4F8.3,8X,4F8.3)
5040  FORMAT(9F8.3)
5050  FORMAT(6F10.2)
5060  FORMAT(7F10.3)
5070  FORMAT(2F10.4,F10.2,F10.3,F10.4,F10.3)
5080  FORMAT(3F10.1)
5090  FORMAT(4E8.1,I8)
5200  FORMAT(4I10)
5300  FORMAT(1X,I2,I3,I5,I6,2F9.2,E14.4,F9.2,E12.3,F8.2,F12.2,3F9.2,
     1       F9.4)
5301  FORMAT(I5,F10.4)
6001  FORMAT('CLASS TEST RUN:     ',6A4)
6002  FORMAT('RESEARCHER:         ',6A4)
6003  FORMAT('INSTITUTION:        ',6A4)
C
C===================== CTEM =============================================== \
C
c     ctem initializations.
c
      if (ctem_on) then
c
c     calculate number of level 2 pfts using modelpft
c
      do 101 j = 1, ican
        isumc = 0
        k1c = (j-1)*l2max + 1
        k2c = k1c + (l2max - 1)
        do n = k1c, k2c
          if(modelpft(n).eq.1) isumc = isumc + 1
        enddo
        nol2pfts(j)=isumc  ! number of level 2 pfts
101   continue
c
      do 110 i=1,nltest
       do 110 m=1,nmtest
        do 111 j = 1, icc
          co2i1csrow(i,m,j)=0.0     !intercellular co2 concentrations
          co2i1cgrow(i,m,j)=0.0
          co2i2csrow(i,m,j)=0.0
          co2i2cgrow(i,m,j)=0.0
          slairow(i,m,j)=0.0        !if bio2str is not called we need to initialize this to zero
111     continue
110   continue
c
      do 123 i =1, ilg
         fsnowacc_m(i)=0.0         !daily accu. fraction of snow
         tcansacc_m(i)=0.0         !daily accu. canopy temp. over snow
         taaccgat_m(i)=0.0
c
         do 128 j = 1, icc
           ancsvgac_m(i,j)=0.0    !daily accu. net photosyn. for canopy over snow subarea
           ancgvgac_m(i,j)=0.0    !daily accu. net photosyn. for canopy over ground subarea
           rmlcsvga_m(i,j)=0.0    !daily accu. leaf respiration for canopy over snow subarea
           rmlcgvga_m(i,j)=0.0    !daily accu. leaf respiration for canopy over ground subarea
           todfrac(i,j)=0.0
128      continue
c
         do 112 j = 1,ignd       !soil temperature and moisture over different subareas
            tbarcacc_m (i,j)=0.0
            tbarcsacc_m(i,j)=0.0
            tbargacc_m (i,j)=0.0
            tbargsacc_m(i,j)=0.0
            thliqcacc_m(i,j)=0.0
            thliqgacc_m(i,j)=0.0
            thicecacc_m(i,j)=0.0
112      continue
123    continue
c
c     find fcancmx with class' fcanmxs and dvdfcans read from ctem's
c     initialization file. this is to divide needle leaf and broad leaf
c     into dcd and evg, and crops and grasses into c3 and c4.
c
      do 113 j = 1, ican
        do 114 i=1,nltest
        do 114 m=1,nmtest
c
          k1c = (j-1)*l2max + 1
          k2c = k1c + (l2max - 1)
c
          do n = k1c, k2c
            if(modelpft(n).eq.1)then
              icountrow(i,m) = icountrow(i,m) + 1
              csum(i,m,j) = csum(i,m,j) +
     &         dvdfcanrow(i,m,icountrow(i,m))

!              Added in seed here to prevent competition from getting
!              pfts with no seed fraction.  JM Feb 20 2014.
              if (compete .and. .not. mosaic) then
               fcancmxrow(i,m,icountrow(i,m))=max(seed,FCANROT(i,m,j)*
     &         dvdfcanrow(i,m,icountrow(i,m)))
               barf(i) = barf(i) - fcancmxrow(i,m,icountrow(i,m))
              else
               fcancmxrow(i,m,icountrow(i,m))=FCANROT(i,m,j)*
     &         dvdfcanrow(i,m,icountrow(i,m))
              end if
            endif
          enddo
c
          if( abs(csum(i,m,j)-1.0).gt.abszero ) then
           write(6,1130)i,m,j
1130       format('dvdfcans for (',i1,',',i1,',',i1,') must add to 1.0')
            call xit('runclass36ctem', -3)
          endif
c
114     continue
113   continue

!     Now make sure that you aren´t over 1.0 for a grid cell (i.e. with a negative
!     bare ground fraction due to the seed fractions being added in.) JM Mar 27 2014
      do i=1,nltest
       if (barf(i) .lt. 0.) then
        bigpftc=maxloc(fcancmxrow(i,:,:))
        ! reduce the most predominant PFT by barf and 1.0e-5,
        ! which ensures that our barefraction is non-zero to avoid
        ! problems later.
        fcancmxrow(i,bigpftc(1),bigpftc(2))=fcancmxrow
     &                (i,bigpftc(1),bigpftc(2))+barf(i) - 1.0e-5
       end if
      end do
c
c     ----------

c     preparation with the input datasets prior to launching run:

      iyear=-99999  ! initialization, forces entry to loop below
      obswetyr=-99999

c     find the first year of met data

       do while (iyear .lt. metcylyrst)
c
        do i=1,nltest  ! formatting was 5300
          read(12,*) ihour,imin,iday,iyear,FSSROW(I),FDLROW(i),
     1         PREROW(i),TAROW(i),QAROW(i),UVROW(i),PRESROW(i)
        enddo
       enddo

c      back up one space in the met file so it is ready for the next readin
       backspace(12)

c  /--------------Rudra-------------/

       if(obswetf) then
         do while (obswetyr .lt. metcylyrst)
            do i=1,nltest
              read(16,*) obswetyr,(wetfrac_mon(i,j),j=1,12)
            end do
         end do
         backspace(16)
       else
           do i=1,nltest
             do j = 1,12
               wetfrac_mon(i,j) = 0.0
             enddo
           enddo

       end if


c    \---------------Rudra----------\

c      If you are not cycling over the MET, you can still specify to end on a
c      year that is shorter than the total climate file length.
       if (.not. cyclemet) endyr = iyear + ncyear

c      find the popd data to cycle over, popd is only cycled over when the met is cycled.
       popyr=-99999  ! initialization, forces entry to loop below

       if (cyclemet .and. popdon) then
        do while (popyr .lt. cypopyr)
         do i = 1, nltest
          read(13,5301) popyr,popdin
         enddo
        enddo
       endif
c
c     if land use change switch is on then read the fractional coverages
c     of ctem's 9 pfts for the first year.
c
      if (lnduseon .and. transient_run) then

         reach_eof=.false.  !flag for when read to end of luc input file

         call initialize_luc(iyear,argbuff,nmtest,nltest,
     1                     mosaic,nol2pfts,cyclemet,
     2                     cylucyr,lucyr,FCANROT,FAREROT,nfcancmxrow,
     3                     pfcancmxrow,fcancmxrow,reach_eof,start_bare,
     4                     compete)

         if (reach_eof) goto 999

      endif ! if (lnduseon)
c
c     with fcancmx calculated above and initialized values of all ctem pools,
c     find mosaic tile (grid) average vegetation biomass, litter mass, and soil c mass.
c     also initialize additional variables which are used by ctem.
c
      do 115 i = 1,nltest
        do 115 m = 1,nmtest
          vgbiomasrow(i,m)=0.0
          gavglairow(i,m)=0.0
          gavgltmsrow(i,m)=0.0
          gavgscmsrow(i,m)=0.0
          lucemcomrow(i,m)=0.0      !land use change combustion emission losses
          lucltrinrow(i,m)=0.0      !land use change inputs to litter pool
          lucsocinrow(i,m)=0.0      !land use change inputs to soil c pool
          colddaysrow(i,m,1)=0      !cold days counter for ndl dcd
          colddaysrow(i,m,2)=0      !cold days counter for crops

          do 116 j = 1, icc
            vgbiomasrow(i,m)=vgbiomasrow(i,m)+fcancmxrow(i,m,j)*
     &        (gleafmasrow(i,m,j)+stemmassrow(i,m,j)+
     &         rootmassrow(i,m,j)+bleafmasrow(i,m,j))
            gavgltmsrow(i,m)=gavgltmsrow(i,m)+fcancmxrow(i,m,j)*
     &                       litrmassrow(i,m,j)
            gavgscmsrow(i,m)=gavgscmsrow(i,m)+fcancmxrow(i,m,j)*
     &         soilcmasrow(i,m,j)
            grwtheffrow(i,m,j)=100.0   !set growth efficiency to some large number
c                                      !so that no growth related mortality occurs in
c                                      !first year
            flhrlossrow(i,m,j)=0.0     !fall/harvest loss
            stmhrlosrow(i,m,j)=0.0     !stem harvest loss for crops
            rothrlosrow(i,m,j)=0.0     !root death for crops
            lystmmasrow(i,m,j)=stemmassrow(i,m,j)
            lyrotmasrow(i,m,j)=rootmassrow(i,m,j)
            tymaxlairow(i,m,j)=0.0

116      continue

c
c *     initialize accumulated array for monthly and yearly output for ctem
c
        if(ctem_on) then

         call resetmonthend_m(nltest,nmtest)
         call resetyearend_m(nltest,nmtest)

!          do j=1,icc
!           npp_mo_m(i,m,j)=0.0
!           gpp_mo_m(i,m,j)=0.0
!           nep_mo_m(i,m,j)=0.0
!           nbp_mo_m(i,m,j)=0.0
!           laimaxg_mo_m(i,m,j)=0.0
!           emit_co2_mo_m(i,m,j)=0.0
!           emit_co_mo_m(i,m,j) =0.0
!           emit_ch4_mo_m(i,m,j) =0.0
!           emit_nmhc_mo_m(i,m,j) =0.0
!           emit_h2_mo_m(i,m,j) =0.0
!           emit_nox_mo_m(i,m,j) =0.0
!           emit_n2o_mo_m(i,m,j) =0.0
!           emit_pm25_mo_m(i,m,j) =0.0
!           emit_tpm_mo_m(i,m,j) =0.0
!           emit_tc_mo_m(i,m,j) =0.0
!           emit_oc_mo_m(i,m,j) =0.0
!           emit_bc_mo_m(i,m,j) =0.0
!           burnfrac_mo_m(i,m,j) =0.0
!
!           laimaxg_yr_m(i,m,j)=0.0
!           npp_yr_m(i,m,j)=0.0
!           gpp_yr_m(i,m,j)=0.0
!           nep_yr_m(i,m,j)=0.0
!           nbp_yr_m(i,m,j)=0.0
!           hetrores_yr_m(i,m,j)=0.0
!           autores_yr_m(i,m,j)=0.0
!           litres_yr_m(i,m,j)=0.0
!           soilcres_yr_m(i,m,j)=0.0
!
!           emit_co2_yr_m(i,m,j)=0.0
!           emit_co_yr_m(i,m,j)=0.0
!           emit_ch4_yr_m(i,m,j)=0.0
!           emit_nmhc_yr_m(i,m,j)=0.0
!           emit_h2_yr_m(i,m,j)=0.0
!           emit_nox_yr_m(i,m,j)=0.0
!           emit_n2o_yr_m(i,m,j)=0.0
!           emit_pm25_yr_m(i,m,j)=0.0
!           emit_tpm_yr_m(i,m,j)=0.0
!           emit_tc_yr_m(i,m,j)=0.0
!           emit_oc_yr_m(i,m,j)=0.0
!           emit_bc_yr_m(i,m,j)=0.0
!           burnfrac_yr_m(i,m,j)=0.0
!
!          end do
!
!           nep_mo_m(i,m,iccp1)=0.0
!           nbp_mo_m(i,m,iccp1)=0.0
!           hetrores_yr_m(i,m,iccp1)=0.0
!           litres_yr_m(i,m,iccp1)=0.0
!           soilcres_yr_m(i,m,iccp1)=0.0
!           nep_yr_m(i,m,iccp1)=0.0
!           nbp_yr_m(i,m,iccp1)=0.0
!
!           probfire_mo_m(i,m) =0.0
!           luc_emc_mo_m(i,m) =0.0
!           lucsocin_mo_m(i,m) =0.0
!           lucltrin_mo_m(i,m) =0.0
!           bterm_mo_m(i,m)=0.0
!           lterm_mo_m(i,m)=0.0
!           mterm_mo_m(i,m)=0.0
!
!           probfire_yr_m(i,m)=0.0
!           luc_emc_yr_m(i,m)=0.0
!           lucsocin_yr_m(i,m)=0.0
!           lucltrin_yr_m(i,m)=0.0
!           bterm_yr_m(i,m)=0.0
!           lterm_yr_m(i,m)=0.0
!           mterm_yr_m(i,m)=0.0
! c                                      !CH4(wetland) related variables !Rudra 04/12/2013
!           ch4wet1_mo_m(i,m)  =0.0
!           ch4wet2_mo_m(i,m)  =0.0
!           wetfdyn_mo_m(i,m)  =0.0
!           ch4dyn1_mo_m(i,m)  =0.0
!           ch4dyn2_mo_m(i,m)  =0.0
!
!           ch4wet1_yr_m(i,m)  =0.0
!           ch4wet2_yr_m(i,m)  =0.0
!           wetfdyn_yr_m(i,m)  =0.0
!           ch4dyn1_yr_m(i,m)  =0.0
!           ch4dyn2_yr_m(i,m)  =0.0
! c
        endif ! ctem_on
c
115   continue
c
      do 117 i = 1,nltest
        do 117 m = 1,nmtest
         gavgltmsrow(i,m)=gavgltmsrow(i,m)+ (1.0-FCANROT(i,m,1)-
     &       FCANROT(i,m,2)-
     &    FCANROT(i,m,3)-FCANROT(i,m,4))*litrmassrow(i,m,icc+1)
         gavgscmsrow(i,m)=gavgscmsrow(i,m)+ (1.0-FCANROT(i,m,1)-
     &   FCANROT(i,m,2)-
     &    FCANROT(i,m,3)-FCANROT(i,m,4))*soilcmasrow(i,m,icc+1)
c
117   continue
c
C===================== CTEM =============================================== /

      CALL GATPREP(ILMOS,JLMOS,IWMOS,JWMOS,
     1             NML,NMW,GCROW,FAREROT,MIDROT,
     2             NLAT,NMOS,ILG,1,NLTEST,NMTEST)

C===================== CTEM =============================================== \

c
      call ctemg1(gleafmasgat,bleafmasgat,stemmassgat,rootmassgat,
     1      fcancmxgat,zbtwgat,dlzwgat,sdepgat,ailcggat,ailcbgat,
     2      ailcgat,zolncgat,rmatcgat,rmatctemgat,slaigat,
     3      bmasveggat,cmasvegcgat,veghghtgat,
     4      rootdpthgat,alvsctmgat,alirctmgat,
     5      paicgat,    slaicgat,
     6      ilmos,jlmos,iwmos,jwmos,
     7      nml,
     8      gleafmasrow,bleafmasrow,stemmassrow,rootmassrow,
     9      fcancmxrow,ZBTWROT,DLZWROT,SDEPROT,ailcgrow,ailcbrow,
     a      ailcrow,zolncrow,rmatcrow,rmatctemrow,slairow,
     b      bmasvegrow,cmasvegcrow,veghghtrow,
     c      rootdpthrow,alvsctmrow,alirctmrow,
     d      paicrow,    slaicrow)
c
c
      call bio2str( gleafmasgat,bleafmasgat,stemmassgat,rootmassgat,
     1                           1,      nml,    fcancmxgat, zbtwgat,
     2                        dlzwgat, nol2pfts,   sdepgat,
     4                       ailcggat, ailcbgat,  ailcgat, zolncgat,
     5                       rmatcgat, rmatctemgat,slaigat,bmasveggat,
     6                 cmasvegcgat,veghghtgat, rootdpthgat,alvsctmgat,
     7                     alirctmgat, paicgat,  slaicgat )
c
      call ctems1(gleafmasrow,bleafmasrow,stemmassrow,rootmassrow,
     1      fcancmxrow,ZBTWROT,DLZWROT,SDEPROT,ailcgrow,ailcbrow,
     2      ailcrow,zolncrow,rmatcrow,rmatctemrow,slairow,
     3      bmasvegrow,cmasvegcrow,veghghtrow,
     4      rootdpthrow,alvsctmrow,alirctmrow,
     5      paicrow,    slaicrow,
     6      ilmos,jlmos,iwmos,jwmos,
     7      nml,
     8      gleafmasgat,bleafmasgat,stemmassgat,rootmassgat,
     9      fcancmxgat,zbtwgat,dlzwgat,sdepgat,ailcggat,ailcbgat,
     a      ailcgat,zolncgat,rmatcgat,rmatctemgat,slaigat,
     b      bmasveggat,cmasvegcgat,veghghtgat,
     c      rootdpthgat,alvsctmgat,alirctmgat,
     d      paicgat,    slaicgat)
c
      endif   ! if (ctem_on)
c
c     ctem initial preparation done

C===================== CTEM ============================================ /
C
C     **** LAUNCH RUN. ****

      N=0
      NCOUNT=1
      NDAY=86400/NINT(DELT)

C===================== CTEM ============================================ \

      run_model=.true.
      met_rewound=.false.

200   continue

c     start up the main model loop

      do while (run_model)


c     if the met file has been rewound (due to cycling over the met data)
c     then we need to find the proper year in the file before we continue
c     on with the run
      if (met_rewound) then
        do while (iyear .lt. metcylyrst)
         do i=1,nltest
c         this reads in one 30 min slice of met data, when it reaches
c         the end of file it will go to label 999.  !formatting was 5300
          read(12,*,end=999) ihour,imin,iday,iyear,FSSROW(I),FDLROW(i),
     1         PREROW(i),TAROW(i),QAROW(i),UVROW(i),PRESROW(i)

         enddo
        enddo

c       back up one space in the met file so it is ready for the next readin
c       but only if it was read in during the loop above.
        if (metcylyrst .ne. -9999) backspace(12)

      if (ctem_on) then
        if (obswetf) then
          do while (obswetyr .lt. metcylyrst)
              do i = 1,nltest
                read(16,*) obswetyr,(wetfrac_mon(i,j),j=1,12)
              enddo
          enddo
         if (metcylyrst .ne. -9999) backspace(16)
        else
           do i=1,nltest
             do j = 1,12
               wetfrac_mon(i,j) = 0.0
             enddo
           enddo
        endif !obswetf
       endif ! ctem_on

      met_rewound = .false.

      endif

C===================== CTEM ============================================ /
C
C========================================================================
C     * READ IN METEOROLOGICAL FORCING DATA FOR CURRENT TIME STEP;
C     * CALCULATE SOLAR ZENITH ANGLE AND COMPONENTS OF INCOMING SHORT-
C     * WAVE RADIATION FLUX; ESTIMATE FLUX PARTITIONS IF NECESSARY.
C
      N=N+1

      DO 250 I=1,NLTEST
C         THIS READS IN ONE 30 MIN SLICE OF MET DATA, WHEN IT REACHES
C         THE END OF FILE IT WILL GO TO 999. !formatting was 5300
          READ(12,*,END=999) IHOUR,IMIN,IDAY,IYEAR,FSSROW(I),FDLROW(I),
     1         PREROW(I),TAROW(I),QAROW(I),UVROW(I),PRESROW(I)

          if (iday.eq.1.and.ihour.eq.0.and.imin.eq.0) then
            if (ctem_on) then
              if (obswetf) then
                  read(16,*,end=1001) obswetyr,(wetfrac_mon(i,j),j=1,12)
              else
                   do j = 1,12
                     wetfrac_mon(i,j) = 0.0
                   enddo
              endif !obswetf
            endif ! ctem_on

          endif


C===================== CTEM ============================================ \

c         assign the met climate year to climiyear
          climiyear = iyear

!         If in a transient_run that has to cycle over MET then change
!         the iyear here:
          if (transient_run .and. cyclemet) then
            iyear = iyear - (metcylyrst - trans_startyr)
          end if
c
          if(lopcount .gt. 1) then
            if (cyclemet) then
              iyear=iyear + nummetcylyrs*(lopcount-1)
            else
              iyear=iyear + ncyear*(lopcount-1)
            end if
          endif   ! lopcount .gt. 1

c
!         write(*,*)'year=',iyear,'day=',iday,' hour=',ihour,' min=',imin
c
C===================== CTEM ============================================ /
          FSVHROW(I)=0.5*FSSROW(I)
          FSIHROW(I)=0.5*FSSROW(I)
          TAROW(I)=TAROW(I)+TFREZ
          ULROW(I)=UVROW(I)
          VLROW(I)=0.0
          VMODROW(I)=UVROW(I)
250   CONTINUE

C
      DAY=REAL(IDAY)+(REAL(IHOUR)+REAL(IMIN)/60.)/24.
      DECL=SIN(2.*PI*(284.+DAY)/365.)*23.45*PI/180.
      HOUR=(REAL(IHOUR)+REAL(IMIN)/60.)*PI/12.-PI
      COSZ=SIN(RADJROW(1))*SIN(DECL)+COS(RADJROW(1))*COS(DECL)*COS(HOUR)

      DO 300 I=1,NLTEST
          CSZROW(I)=SIGN(MAX(ABS(COSZ),1.0E-3),COSZ)
          IF(PREROW(I).GT.0.) THEN
              XDIFFUS(I)=1.0
          ELSE
              XDIFFUS(I)=MAX(0.0,MIN(1.0-0.9*COSZ,1.0))
          ENDIF
          FCLOROW(I)=XDIFFUS(I)
300   CONTINUE
C
C===================== CTEM ============================================ \
C
      if (iday.eq.1.and.ihour.eq.0.and.imin.eq.0) then
c
c      if popdon=true
c      calculate fire extinguishing probability and
c      probability of fire due to human causes
c      from population density input data. In disturb.f90 this will
c      overwrite extnprobgrd(i) and prbfrhucgrd(i) that are
c      read in from the .ctm file. Set
c      cypopyr = -9999 when we don't want to cycle over the popd data
c      so this allows us to grab a new value each year.

       if(popdon .and. cypopyr .eq. -9999) then
         do while (popyr .lt. iyear)
          do i=1,nltest
           read(13,5301,end=999) popyr,popdin
          enddo
         enddo
       endif
c
c      if co2on is true
c      read co2concin from input datafile and
c      overwrite co2concrow, otherwise set to constant value
c
       if(co2on) then

        do while (co2yr .lt. iyear)
          do i=1,nltest
           read(14,*,end=999) co2yr,co2concin
           do m=1,nmtest
            co2concrow(i,m)=co2concin
           enddo !nmtest
          enddo !nltest
        enddo !co2yr < iyear

       else !constant co2

         do i=1,nltest
          do m=1,nmtest
           co2concrow(i,m)=setco2conc
          enddo
         enddo

       endif !co2on

c      if lnduseon is true, read in the luc data now

       if (ctem_on .and. lnduseon .and. transient_run) then

         call readin_luc(iyear,nmtest,nltest,mosaic,lucyr,
     &                   nfcancmxrow,pfcancmxrow,reach_eof,compete)
         if (reach_eof) goto 999

       else ! lnduseon = false or met is cycling in a spin up run

c          land use is not on or the met data is being cycled, so the
c          pfcancmx value is also the nfcancmx value.
c
           nfcancmxrow=pfcancmxrow

       endif ! lnduseon/cyclemet
c
      endif   ! at the first day of each year i.e.
c             ! if (iday.eq.1.and.ihour.eq.0.and.imin.eq.0)
c

C===================== CTEM ============================================ /
C
      CALL CLASSI(VPDROW,TADPROW,PADRROW,RHOAROW,RHSIROW,
     1            RPCPROW,TRPCROW,SPCPROW,TSPCROW,TAROW,QAROW,
     2            PREROW,RPREROW,SPREROW,PRESROW,
     3            IPCP,NLAT,1,NLTEST)

C
      CUMSNO=CUMSNO+SPCPROW(1)*RHSIROW(1)*DELT
C
      CALL GATPREP(ILMOS,JLMOS,IWMOS,JWMOS,
     1             NML,NMW,GCROW,FAREROT,MIDROT,
     2             NLAT,NMOS,ILG,1,NLTEST,NMTEST)
C
      CALL CLASSG (TBARGAT,THLQGAT,THICGAT,TPNDGAT,ZPNDGAT,
     1             TBASGAT,ALBSGAT,TSNOGAT,RHOSGAT,SNOGAT,
     2             TCANGAT,RCANGAT,SCANGAT,GROGAT, CMAIGAT,
     3             FCANGAT,LNZ0GAT,ALVCGAT,ALICGAT,PAMXGAT,
     4             PAMNGAT,CMASGAT,ROOTGAT,RSMNGAT,QA50GAT,
     5             VPDAGAT,VPDBGAT,PSGAGAT,PSGBGAT,PAIDGAT,
     6             HGTDGAT,ACVDGAT,ACIDGAT,TSFSGAT,WSNOGAT,
     7             THPGAT, THRGAT, THMGAT, BIGAT,  PSISGAT,
     8             GRKSGAT,THRAGAT,HCPSGAT,TCSGAT, IGDRGAT,
     9             THFCGAT,THLWGAT,PSIWGAT,DLZWGAT,ZBTWGAT,
     A             VMODGAT,ZSNLGAT,ZPLGGAT,ZPLSGAT,TACGAT,
     B             QACGAT,DRNGAT, XSLPGAT,GRKFGAT,WFSFGAT,
     C             WFCIGAT,ALGWVGAT,ALGWNGAT,ALGDVGAT,ALGDNGAT,
     +             ALGWGAT,ALGDGAT,ASVDGAT,ASIDGAT,AGVDGAT,
     D             AGIDGAT,ISNDGAT,RADJGAT,ZBLDGAT,Z0ORGAT,
     E             ZRFMGAT,ZRFHGAT,ZDMGAT, ZDHGAT, FSVHGAT,
     F             FSIHGAT,FSDBGAT,FSFBGAT,FSSBGAT,CSZGAT,
     +             FSGGAT, FLGGAT, FDLGAT, ULGAT,  VLGAT,
     G             TAGAT,  QAGAT,  PRESGAT,PREGAT, PADRGAT,
     H             VPDGAT, TADPGAT,RHOAGAT,RPCPGAT,TRPCGAT,
     I             SPCPGAT,TSPCGAT,RHSIGAT,FCLOGAT,DLONGAT,
     J             GGEOGAT,GUSTGAT,REFGAT, BCSNGAT,DEPBGAT,
     K             ILMOS,JLMOS,
     L             NML,NLAT,NTLD,NMOS,ILG,IGND,ICAN,ICAN+1,NBS,
     M             TBARROT,THLQROT,THICROT,TPNDROT,ZPNDROT,
     N             TBASROT,ALBSROT,TSNOROT,RHOSROT,SNOROT,
     O             TCANROT,RCANROT,SCANROT,GROROT, CMAIROT,
     P             FCANROT,LNZ0ROT,ALVCROT,ALICROT,PAMXROT,
     Q             PAMNROT,CMASROT,ROOTROT,RSMNROT,QA50ROT,
     R             VPDAROT,VPDBROT,PSGAROT,PSGBROT,PAIDROT,
     S             HGTDROT,ACVDROT,ACIDROT,TSFSROT,WSNOROT,
     T             THPROT, THRROT, THMROT, BIROT,  PSISROT,
     U             GRKSROT,THRAROT,HCPSROT,TCSROT, IGDRROT,
     V             THFCROT,THLWROT,PSIWROT,DLZWROT,ZBTWROT,
     W             VMODROW,ZSNLROT,ZPLGROT,ZPLSROT,TACROT,
     X             QACROT,DRNROT, XSLPROT,GRKFROT,WFSFROT,
     Y             WFCIROT,ALGWVROT,ALGWNROT,ALGDVROT,ALGDNROT,
     +             ALGWROT,ALGDROT,ASVDROT,ASIDROT,AGVDROT,
     Z             AGIDROT,ISNDROT,RADJROW,ZBLDROW,Z0ORROW,
     +             ZRFMROW,ZRFHROW,ZDMROW, ZDHROW, FSVHROW,
     +             FSIHROW,FSDBROL,FSFBROL,FSSBROL,CSZROW,
     +             FSGROL, FLGROL, FDLROW, ULROW,  VLROW,
     +             TAROW,  QAROW,  PRESROW,PREROW, PADRROW,
     +             VPDROW, TADPROW,RHOAROW,RPCPROW,TRPCROW,
     +             SPCPROW,TSPCROW,RHSIROW,FCLOROW,DLONROW,
     +             GGEOROW,GUSTROL,REFROT, BCSNROT,DEPBROW )
C
C    * INITIALIZATION OF DIAGNOSTIC VARIABLES SPLIT OUT OF CLASSG
C    * FOR CONSISTENCY WITH GCM APPLICATIONS.
C
      DO 330 K=1,ILG
          CDHGAT (K)=0.0
          CDMGAT (K)=0.0
          HFSGAT (K)=0.0
          TFXGAT (K)=0.0
          QEVPGAT(K)=0.0
          QFSGAT (K)=0.0
          QFXGAT (K)=0.0
          PETGAT (K)=0.0
          GAGAT  (K)=0.0
          EFGAT  (K)=0.0
          GTGAT  (K)=0.0
          QGGAT  (K)=0.0
          ALVSGAT(K)=0.0
          ALIRGAT(K)=0.0
          SFCTGAT(K)=0.0
          SFCUGAT(K)=0.0
          SFCVGAT(K)=0.0
          SFCQGAT(K)=0.0
          FSNOGAT(K)=0.0
          FSGVGAT(K)=0.0
          FSGSGAT(K)=0.0
          FSGGGAT(K)=0.0
          FLGVGAT(K)=0.0
          FLGSGAT(K)=0.0
          FLGGGAT(K)=0.0
          HFSCGAT(K)=0.0
          HFSSGAT(K)=0.0
          HFSGGAT(K)=0.0
          HEVCGAT(K)=0.0
          HEVSGAT(K)=0.0
          HEVGGAT(K)=0.0
          HMFCGAT(K)=0.0
          HMFNGAT(K)=0.0
          HTCCGAT(K)=0.0
          HTCSGAT(K)=0.0
          PCFCGAT(K)=0.0
          PCLCGAT(K)=0.0
          PCPNGAT(K)=0.0
          PCPGGAT(K)=0.0
          QFGGAT (K)=0.0
          QFNGAT (K)=0.0
          QFCFGAT(K)=0.0
          QFCLGAT(K)=0.0
          ROFGAT (K)=0.0
          ROFOGAT(K)=0.0
          ROFSGAT(K)=0.0
          ROFBGAT(K)=0.0
          TROFGAT(K)=0.0
          TROOGAT(K)=0.0
          TROSGAT(K)=0.0
          TROBGAT(K)=0.0
          ROFCGAT(K)=0.0
          ROFNGAT(K)=0.0
          ROVGGAT(K)=0.0
          WTRCGAT(K)=0.0
          WTRSGAT(K)=0.0
          WTRGGAT(K)=0.0
          DRGAT  (K)=0.0
330   CONTINUE
C
      DO 334 L=1,IGND
      DO 332 K=1,ILG
          HMFGGAT(K,L)=0.0
          HTCGAT (K,L)=0.0
          QFCGAT (K,L)=0.0
          GFLXGAT(K,L)=0.0
332   CONTINUE
334   CONTINUE
C
      DO 340 M=1,50
          DO 338 L=1,6
              DO 336 K=1,NML
                  ITCTGAT(K,L,M)=0
336           CONTINUE
338       CONTINUE
340   CONTINUE
C
C========================================================================
C
      CALL CLASSZ (0,      CTVSTP, CTSSTP, CT1STP, CT2STP, CT3STP,
     1             WTVSTP, WTSSTP, WTGSTP,
     2             FSGVGAT,FLGVGAT,HFSCGAT,HEVCGAT,HMFCGAT,HTCCGAT,
     3             FSGSGAT,FLGSGAT,HFSSGAT,HEVSGAT,HMFNGAT,HTCSGAT,
     4             FSGGGAT,FLGGGAT,HFSGGAT,HEVGGAT,HMFGGAT,HTCGAT,
     5             PCFCGAT,PCLCGAT,QFCFGAT,QFCLGAT,ROFCGAT,WTRCGAT,
     6             PCPNGAT,QFNGAT, ROFNGAT,WTRSGAT,PCPGGAT,QFGGAT,
     7             QFCGAT, ROFGAT, WTRGGAT,CMAIGAT,RCANGAT,SCANGAT,
     8             TCANGAT,SNOGAT, WSNOGAT,TSNOGAT,THLQGAT,THICGAT,
     9             HCPSGAT,THPGAT, DLZWGAT,TBARGAT,ZPNDGAT,TPNDGAT,
     A             DELZ,   FCS,    FGS,    FC,     FG,
     B             1,      NML,    ILG,    IGND,   N    )
C
C========================================================================
C
C===================== CTEM ============================================ \
C
      call ctemg2(fcancmxgat,rmatcgat,zolncgat,paicgat,
     1      ailcgat,     ailcggat,    cmasvegcgat,  slaicgat,
     2      ailcgsgat,   fcancsgat,   fcancgat,     rmatctemgat,
     3      co2concgat,  co2i1cggat,  co2i1csgat,   co2i2cggat,
     4      co2i2csgat,  xdiffusgat,  slaigat,      cfluxcggat,
     5      cfluxcsgat,  ancsveggat,  ancgveggat,   rmlcsveggat,
     6      rmlcgveggat, canresgat,   sdepgat,
     7      sandgat,     claygat,     orgmgat,
     8      anveggat,    rmlveggat,   tcanoaccgat_m,tbaraccgat_m,
     9      uvaccgat_m,  vvaccgat_m,  mlightnggat,  prbfrhucgat,
     a      extnprobgat, stdalngat,   pfcancmxgat,  nfcancmxgat,
     b      stemmassgat, rootmassgat, litrmassgat,  gleafmasgat,
     c      bleafmasgat, soilcmasgat, ailcbgat,     flhrlossgat,
     d      pandaysgat,  lfstatusgat, grwtheffgat,  lystmmasgat,
     e      lyrotmasgat, tymaxlaigat, vgbiomasgat,  gavgltmsgat,
     f      stmhrlosgat, bmasveggat,  colddaysgat,  rothrlosgat,
     g      alvsctmgat,  alirctmgat,  gavglaigat,   nppgat,
     h      nepgat,      hetroresgat, autoresgat,   soilcrespgat,
     i      rmgat,       rggat,       nbpgat,       litresgat,
     j      socresgat,   gppgat,      dstcemlsgat,  litrfallgat,
     k      humiftrsgat, veghghtgat,  rootdpthgat,  rmlgat,
     l      rmsgat,      rmrgat,      tltrleafgat,  tltrstemgat,
     m      tltrrootgat, leaflitrgat, roottempgat,  afrleafgat,
     n      afrstemgat,  afrrootgat,  wtstatusgat,  ltstatusgat,
     o      burnfracgat, probfiregat, lucemcomgat,  lucltringat,
     p      lucsocingat, nppveggat,   dstcemls3gat,
     q      faregat,     gavgscmsgat, rmlvegaccgat, pftexistgat,
     &      rmsveggat,   rmrveggat,   rgveggat,    vgbiomas_veggat,
     &      gppveggat,   nepveggat,   ailcmingat,   ailcmaxgat,
     &      emit_co2gat,  emit_cogat, emit_ch4gat,  emit_nmhcgat,
     &      emit_h2gat,   emit_noxgat,emit_n2ogat,  emit_pm25gat,
     &      emit_tpmgat,  emit_tcgat, emit_ocgat,   emit_bcgat,
     &      btermgat,     ltermgat,   mtermgat,
     &      nbpveggat,    hetroresveggat, autoresveggat,litresveggat,
     &      soilcresveggat, burnvegfgat, pstemmassgat, pgleafmassgat,
     &      CH4WET1GAT, CH4WET2GAT,
     &      WETFDYNGAT, CH4DYN1GAT,  CH4DYN2GAT,
c
     r      ilmos,       jlmos,       iwmos,        jwmos,
     s      nml,      fcancmxrow,  rmatcrow,    zolncrow,  paicrow,
     v      ailcrow,     ailcgrow,    cmasvegcrow,  slaicrow,
     w      ailcgsrow,   fcancsrow,   fcancrow,     rmatctemrow,
     x      co2concrow,  co2i1cgrow,  co2i1csrow,   co2i2cgrow,
     y      co2i2csrow,  xdiffus,     slairow,      cfluxcgrow,
     z      cfluxcsrow,  ancsvegrow,  ancgvegrow,   rmlcsvegrow,
     1      rmlcgvegrow, canresrow,   SDEPROT,
     2      SANDROT,     CLAYROT,     ORGMROT,
     3      anvegrow,    rmlvegrow,   tcanoaccrow_m,tbaraccrow_m,
     4      uvaccrow_m,  vvaccrow_m,  mlightnggrd,  prbfrhucgrd,
     5      extnprobgrd, stdalngrd,   pfcancmxrow,  nfcancmxrow,
     6      stemmassrow, rootmassrow, litrmassrow,  gleafmasrow,
     7      bleafmasrow, soilcmasrow, ailcbrow,     flhrlossrow,
     8      pandaysrow,  lfstatusrow, grwtheffrow,  lystmmasrow,
     9      lyrotmasrow, tymaxlairow, vgbiomasrow,  gavgltmsrow,
     a      stmhrlosrow, bmasvegrow,  colddaysrow,  rothrlosrow,
     b      alvsctmrow,  alirctmrow,  gavglairow,   npprow,
     c      neprow,      hetroresrow, autoresrow,   soilcresprow,
     d      rmrow,       rgrow,       nbprow,       litresrow,
     e      socresrow,   gpprow,      dstcemlsrow,  litrfallrow,
     f      humiftrsrow, veghghtrow,  rootdpthrow,  rmlrow,
     g      rmsrow,      rmrrow,      tltrleafrow,  tltrstemrow,
     h      tltrrootrow, leaflitrrow, roottemprow,  afrleafrow,
     i      afrstemrow,  afrrootrow,  wtstatusrow,  ltstatusrow,
     j      burnfracrow, probfirerow, lucemcomrow,  lucltrinrow,
     k      lucsocinrow, nppvegrow,   dstcemls3row,
     l      FAREROT,     gavgscmsrow, rmlvegaccrow, pftexistrow,
     &      rmsvegrow,   rmrvegrow,   rgvegrow,    vgbiomas_vegrow,
     &      gppvegrow,   nepvegrow,   ailcminrow,   ailcmaxrow,
     &      emit_co2row,  emit_corow, emit_ch4row,  emit_nmhcrow,
     &      emit_h2row,   emit_noxrow,emit_n2orow,  emit_pm25row,
     &      emit_tpmrow,  emit_tcrow, emit_ocrow,   emit_bcrow,
     &      btermrow,     ltermrow,   mtermrow,
     &      nbpvegrow,    hetroresvegrow, autoresvegrow,litresvegrow,
     &      soilcresvegrow, burnvegfrow, pstemmassrow, pgleafmassrow,
     &      CH4WET1ROW, CH4WET2ROW,
     &      WETFDYNROW, CH4DYN1ROW, CH4DYN2ROW)
c
C===================== CTEM ============================================ /
C
C-----------------------------------------------------------------------
C     * ALBEDO AND TRANSMISSIVITY CALCULATIONS; GENERAL VEGETATION
C     * CHARACTERISTICS.
C     * ADAPTED TO COUPLING OF CLASS3.6 AND CTEM
C
      CALL CLASSA    (FC,     FG,     FCS,    FGS,    ALVSCN, ALIRCN,
     1                ALVSG,  ALIRG,  ALVSCS, ALIRCS, ALVSSN, ALIRSN,
     2                ALVSGC, ALIRGC, ALVSSC, ALIRSC, TRVSCN, TRIRCN,
     3                TRVSCS, TRIRCS, FSVF,   FSVFS,
     4                RAICAN, RAICNS, SNOCAN, SNOCNS, FRAINC, FSNOWC,
     5                FRAICS, FSNOCS, DISP,   DISPS,  ZOMLNC, ZOMLCS,
     6                ZOELNC, ZOELCS, ZOMLNG, ZOMLNS, ZOELNG, ZOELNS,
     7                CHCAP,  CHCAPS, CMASSC, CMASCS, CWLCAP, CWFCAP,
     8                CWLCPS, CWFCPS, RC,     RCS,    RBCOEF, FROOT,
     9                ZPLIMC, ZPLIMG, ZPLMCS, ZPLMGS, ZSNOW,
     A                WSNOGAT,ALVSGAT,ALIRGAT,HTCCGAT,HTCSGAT,HTCGAT,
     +                ALTG,   ALSNO,  TRSNOWC,TRSNOWG,
     B                WTRCGAT,WTRSGAT,WTRGGAT,CMAIGAT,FSNOGAT,
     C                FCANGAT,LNZ0GAT,ALVCGAT,ALICGAT,PAMXGAT,PAMNGAT,
     D                CMASGAT,ROOTGAT,RSMNGAT,QA50GAT,VPDAGAT,VPDBGAT,
     E                PSGAGAT,PSGBGAT,PAIDGAT,HGTDGAT,ACVDGAT,ACIDGAT,
     F                ASVDGAT,ASIDGAT,AGVDGAT,AGIDGAT,ALGWGAT,ALGDGAT,
     +                ALGWVGAT,ALGWNGAT,ALGDVGAT,ALGDNGAT,
     G                THLQGAT,THICGAT,TBARGAT,RCANGAT,SCANGAT,TCANGAT,
     H                GROGAT, SNOGAT, TSNOGAT,RHOSGAT,ALBSGAT,ZBLDGAT,
     I                Z0ORGAT,ZSNLGAT,ZPLGGAT,ZPLSGAT,
     J                FCLOGAT,TAGAT,  VPDGAT, RHOAGAT,CSZGAT,
     +                FSDBGAT,FSFBGAT,REFGAT, BCSNGAT,
     K                FSVHGAT,RADJGAT,DLONGAT,RHSIGAT,DELZ,   DLZWGAT,
     L                ZBTWGAT,THPGAT, THMGAT, PSISGAT,BIGAT,  PSIWGAT,
     M                HCPSGAT,ISNDGAT,
     P                FCANCMXGAT,ICC,ICTEMMOD,RMATCGAT,ZOLNCGAT,
     Q                CMASVEGCGAT,AILCGAT,PAICGAT,L2MAX, NOL2PFTS,
     R                SLAICGAT,AILCGGAT,AILCGSGAT,FCANCGAT,FCANCSGAT,
     R                IDAY,   ILG,    1,      NML,  NBS,
     N                JLAT,N, ICAN,   ICAN+1, IGND,   IDISP,  IZREF,
     O                IWF,    IPAI,   IHGT,   IALC,   IALS,   IALG,
     P                ISNOALB,IGRALB, alvsctmgat,alirctmgat )
C
C-----------------------------------------------------------------------
C          * SURFACE TEMPERATURE AND FLUX CALCULATIONS.
C          * ADAPTED TO COUPLING OF CLASS3.6 AND CTEM
C
      CALL CLASST     (TBARC,  TBARG,  TBARCS, TBARGS, THLIQC, THLIQG,
     1  THICEC, THICEG, HCPC,   HCPG,   TCTOPC, TCBOTC, TCTOPG, TCBOTG,
     2  GZEROC, GZEROG, GZROCS, GZROGS, G12C,   G12G,   G12CS,  G12GS,
     3  G23C,   G23G,   G23CS,  G23GS,  QFREZC, QFREZG, QMELTC, QMELTG,
     4  EVAPC,  EVAPCG, EVAPG,  EVAPCS, EVPCSG, EVAPGS, TCANO,  TCANS,
     5  RAICAN, SNOCAN, RAICNS, SNOCNS, CHCAP,  CHCAPS, TPONDC, TPONDG,
     6  TPNDCS, TPNDGS, TSNOCS, TSNOGS, WSNOCS, WSNOGS, RHOSCS, RHOSGS,
     7  ITCTGAT,CDHGAT, CDMGAT, HFSGAT, TFXGAT, QEVPGAT,QFSGAT,
     8  PETGAT, GAGAT,  EFGAT,  GTGAT,  QGGAT,
     +  SFCTGAT,SFCUGAT,SFCVGAT,SFCQGAT,SFRHGAT,
     +  GTBS,   SFCUBS, SFCVBS, USTARBS,
     9  FSGVGAT,FSGSGAT,FSGGGAT,FLGVGAT,FLGSGAT,FLGGGAT,
     A  HFSCGAT,HFSSGAT,HFSGGAT,HEVCGAT,HEVSGAT,HEVGGAT,HMFCGAT,HMFNGAT,
     B  HTCCGAT,HTCSGAT,HTCGAT, QFCFGAT,QFCLGAT,DRGAT,  WTABGAT,ILMOGAT,
     C  UEGAT,  HBLGAT, TACGAT, QACGAT, ZRFMGAT,ZRFHGAT,ZDMGAT, ZDHGAT,
     D  VPDGAT, TADPGAT,RHOAGAT,FSVHGAT,FSIHGAT,FDLGAT, ULGAT,  VLGAT,
     E  TAGAT,  QAGAT,  PADRGAT,FC,     FG,     FCS,    FGS,    RBCOEF,
     F  FSVF,   FSVFS,  PRESGAT,VMODGAT,ALVSCN, ALIRCN, ALVSG,  ALIRG,
     G  ALVSCS, ALIRCS, ALVSSN, ALIRSN, ALVSGC, ALIRGC, ALVSSC, ALIRSC,
     H  TRVSCN, TRIRCN, TRVSCS, TRIRCS, RC,     RCS,    WTRGGAT,QLWOGAT,
     I  FRAINC, FSNOWC, FRAICS, FSNOCS, CMASSC, CMASCS, DISP,   DISPS,
     J  ZOMLNC, ZOELNC, ZOMLNG, ZOELNG, ZOMLCS, ZOELCS, ZOMLNS, ZOELNS,
     K  TBARGAT,THLQGAT,THICGAT,TPNDGAT,ZPNDGAT,TBASGAT,TCANGAT,TSNOGAT,
     L  ZSNOW,  RHOSGAT,WSNOGAT,THPGAT, THRGAT, THMGAT, THFCGAT,THLWGAT,
     +  TRSNOWC,TRSNOWG,ALSNO,  FSSBGAT,
     M  RADJGAT,PREGAT, HCPSGAT,TCSGAT, TSFSGAT,DELZ,   DLZWGAT,ZBTWGAT,
     N  FTEMP,  FVAP,   RIB,    ISNDGAT,
     O  AILCGGAT,  AILCGSGAT, FCANCGAT,FCANCSGAT,CO2CONCGAT,CO2I1CGGAT,
     P  CO2I1CSGAT,CO2I2CGGAT,CO2I2CSGAT,CSZGAT,XDIFFUSGAT,SLAIGAT,ICC,
     Q  ICTEMMOD,RMATCTEMGAT,FCANCMXGAT,L2MAX,  NOL2PFTS,CFLUXCGGAT,
     R  CFLUXCSGAT,ANCSVEGGAT,ANCGVEGGAT,RMLCSVEGGAT,RMLCGVEGGAT,
     S  TCSNOW,GSNOW,ITC,ITCG,ITG,    ILG,    1,NML,  JLAT,N, ICAN,
     T  IGND,   IZREF,  ISLFD,  NLANDCS,NLANDGS,NLANDC, NLANDG, NLANDI,
     U  NBS,    ISNOALB,lfstatusgat)
C
C-----------------------------------------------------------------------
C          * WATER BUDGET CALCULATIONS.
C

          CALL CLASSW  (THLQGAT,THICGAT,TBARGAT,TCANGAT,RCANGAT,SCANGAT,
     1                  ROFGAT, TROFGAT,SNOGAT, TSNOGAT,RHOSGAT,ALBSGAT,
     2                  WSNOGAT,ZPNDGAT,TPNDGAT,GROGAT, TBASGAT,GFLXGAT,
     3                  PCFCGAT,PCLCGAT,PCPNGAT,PCPGGAT,QFCFGAT,QFCLGAT,
     4                  QFNGAT, QFGGAT, QFCGAT, HMFCGAT,HMFGGAT,HMFNGAT,
     5                  HTCCGAT,HTCSGAT,HTCGAT, ROFCGAT,ROFNGAT,ROVGGAT,
     6                  WTRSGAT,WTRGGAT,ROFOGAT,ROFSGAT,ROFBGAT,
     7                  TROOGAT,TROSGAT,TROBGAT,QFSGAT, QFXGAT, RHOAGAT,
     8                  TBARC,  TBARG,  TBARCS, TBARGS, THLIQC, THLIQG,
     9                  THICEC, THICEG, HCPC,   HCPG,   RPCPGAT,TRPCGAT,
     A                  SPCPGAT,TSPCGAT,PREGAT, TAGAT,  RHSIGAT,GGEOGAT,
     B                  FC,     FG,     FCS,    FGS,    TPONDC, TPONDG,
     C                  TPNDCS, TPNDGS, EVAPC,  EVAPCG, EVAPG,  EVAPCS,
     D                  EVPCSG, EVAPGS, QFREZC, QFREZG, QMELTC, QMELTG,
     E                  RAICAN, SNOCAN, RAICNS, SNOCNS, FROOT,  FSVF,
     F                  FSVFS,  CWLCAP, CWFCAP, CWLCPS, CWFCPS, TCANO,
     G                  TCANS,  CHCAP,  CHCAPS, CMASSC, CMASCS, ZSNOW,
     H                  GZEROC, GZEROG, GZROCS, GZROGS, G12C,   G12G,
     I                  G12CS,  G12GS,  G23C,   G23G,   G23CS,  G23GS,
     J                  TSNOCS, TSNOGS, WSNOCS, WSNOGS, RHOSCS, RHOSGS,
     K                  ZPLIMC, ZPLIMG, ZPLMCS, ZPLMGS, TSFSGAT,
     L                  TCTOPC, TCBOTC, TCTOPG, TCBOTG,
     M                  THPGAT, THRGAT, THMGAT, BIGAT,  PSISGAT,GRKSGAT,
     N                  THRAGAT,THFCGAT,DRNGAT, HCPSGAT,DELZ,
     O                  DLZWGAT,ZBTWGAT,XSLPGAT,GRKFGAT,WFSFGAT,WFCIGAT,
     P                  ISNDGAT,IGDRGAT,
     Q                  IWF,    ILG,    1,      NML,    N,
     R                  JLAT,   ICAN,   IGND,   IGND+1, IGND+2,
     S                  NLANDCS,NLANDGS,NLANDC, NLANDG, NLANDI )
C

C========================================================================
C
      CALL CLASSZ (1,      CTVSTP, CTSSTP, CT1STP, CT2STP, CT3STP,
     1             WTVSTP, WTSSTP, WTGSTP,
     2             FSGVGAT,FLGVGAT,HFSCGAT,HEVCGAT,HMFCGAT,HTCCGAT,
     3             FSGSGAT,FLGSGAT,HFSSGAT,HEVSGAT,HMFNGAT,HTCSGAT,
     4             FSGGGAT,FLGGGAT,HFSGGAT,HEVGGAT,HMFGGAT,HTCGAT,
     5             PCFCGAT,PCLCGAT,QFCFGAT,QFCLGAT,ROFCGAT,WTRCGAT,
     6             PCPNGAT,QFNGAT, ROFNGAT,WTRSGAT,PCPGGAT,QFGGAT,
     7             QFCGAT, ROFGAT, WTRGGAT,CMAIGAT,RCANGAT,SCANGAT,
     8             TCANGAT,SNOGAT, WSNOGAT,TSNOGAT,THLQGAT,THICGAT,
     9             HCPSGAT,THPGAT, DLZWGAT,TBARGAT,ZPNDGAT,TPNDGAT,
     A             DELZ,   FCS,    FGS,    FC,     FG,
     B             1,      NML,    ILG,    IGND,   N    )
C
C=======================================================================

c     * accumulate output data for running ctem.
c
      do 660 i=1,nml
         uvaccgat_m(i)=uvaccgat_m(i)+ulgat(i)
660   continue
c
c     accumulate variables not already accumulated but which are required by
c     ctem.
c
      if (ctem_on) then
        do 700 i = 1, nml
c
          alswacc_gat(i)=alswacc_gat(i)+alvsgat(i)*fsvhgat(i)
          allwacc_gat(i)=allwacc_gat(i)+alirgat(i)*fsihgat(i)
          fsinacc_gat(i)=fsinacc_gat(i)+FSSROW(I)
          flinacc_gat(i)=flinacc_gat(i)+fdlgat(i)
          flutacc_gat(i)=flutacc_gat(i)+sbc*gtgat(i)**4
          pregacc_gat(i)=pregacc_gat(i)+pregat(i)*delt
c
          fsnowacc_m(i)=fsnowacc_m(i)+fsnogat(i)
          tcanoaccgat_m(i)=tcanoaccgat_m(i)+tcano(i)
          tcansacc_m(i)=tcansacc_m(i)+tcans(i)
          taaccgat_m(i)=taaccgat_m(i)+tagat(i)
          vvaccgat_m(i)=vvaccgat_m(i)+ vlgat(i)
c
          do 710 j=1,ignd
             tbaraccgat_m(i,j)=tbaraccgat_m(i,j)+tbargat(i,j)
             tbarcacc_m(i,j)=tbarcacc_m(i,j)+tbarc(i,j)
             tbarcsacc_m(i,j)=tbarcsacc_m(i,j)+tbarcs(i,j)
             tbargacc_m(i,j)=tbargacc_m(i,j)+tbarg(i,j)
             tbargsacc_m(i,j)=tbargsacc_m(i,j)+tbargs(i,j)
             thliqcacc_m(i,j)=thliqcacc_m(i,j)+thliqc(i,j)
             thliqgacc_m(i,j)=thliqgacc_m(i,j)+thliqg(i,j)
             thicecacc_m(i,j)=thicecacc_m(i,j)+thicec(i,j)
710       continue
c
          do 713 j = 1, icc
            ancsvgac_m(i,j)=ancsvgac_m(i,j)+ancsveggat(i,j)
            ancgvgac_m(i,j)=ancgvgac_m(i,j)+ancgveggat(i,j)
            rmlcsvga_m(i,j)=rmlcsvga_m(i,j)+rmlcsveggat(i,j)
            rmlcgvga_m(i,j)=rmlcgvga_m(i,j)+rmlcgveggat(i,j)
713       continue
c
700     continue
      endif !if (ctem_on)
c
      if(ncount.eq.nday) then
c
      do 855 i=1,nml
          uvaccgat_m(i)=uvaccgat_m(i)/real(nday)
          vvaccgat_m(i)=vvaccgat_m(i)/real(nday)
c
c         daily averages of accumulated variables for ctem
c
          if (ctem_on) then
c
c           net radiation and precipitation estimates for ctem's bioclim
c
            if(fsinacc_gat(i).gt.0.0) then
              alswacc_gat(i)=alswacc_gat(i)/(fsinacc_gat(i)*0.5)
              allwacc_gat(i)=allwacc_gat(i)/(fsinacc_gat(i)*0.5)
            else
              alswacc_gat(i)=0.0
              allwacc_gat(i)=0.0
            endif
c
            fsinacc_gat(i)=fsinacc_gat(i)/real(nday)
            flinacc_gat(i)=flinacc_gat(i)/real(nday)
            flutacc_gat(i)=flutacc_gat(i)/real(nday)
c
            altot_gat=(alswacc_gat(i)+allwacc_gat(i))/2.0
            fsstar_gat=fsinacc_gat(i)*(1.-altot_gat)
            flstar_gat=flinacc_gat(i)-flutacc_gat(i)
            netrad_gat(i)=fsstar_gat+flstar_gat
            preacc_gat(i)=pregacc_gat(i)
c
            fsnowacc_m(i)=fsnowacc_m(i)/real(nday)
            tcanoaccgat_m(i)=tcanoaccgat_m(i)/real(nday)
            tcansacc_m(i)=tcansacc_m(i)/real(nday)
            taaccgat_m(i)=taaccgat_m(i)/real(nday)
c
            do 831 j=1,ignd
              tbaraccgat_m(i,j)=tbaraccgat_m(i,j)/real(nday)
              tbarcacc_m(i,j) = tbaraccgat_m(i,j)
              tbarcsacc_m(i,j) = tbaraccgat_m(i,j)
              tbargacc_m(i,j) = tbaraccgat_m(i,j)
              tbargsacc_m(i,j) = tbaraccgat_m(i,j)
c
              thliqcacc_m(i,j)=thliqcacc_m(i,j)/real(nday)
              thliqgacc_m(i,j)=thliqgacc_m(i,j)/real(nday)
              thicecacc_m(i,j)=thicecacc_m(i,j)/real(nday)
831         continue
c
            do 832 j = 1, icc
              ancsvgac_m(i,j)=ancsvgac_m(i,j)/real(nday)
              ancgvgac_m(i,j)=ancgvgac_m(i,j)/real(nday)
              rmlcsvga_m(i,j)=rmlcsvga_m(i,j)/real(nday)
              rmlcgvga_m(i,j)=rmlcgvga_m(i,j)/real(nday)
832         continue
c
c           pass on mean monthly lightning for the current month to ctem
c           lightng(i)=mlightng(i,month)
c
c           in a very simple way try to interpolate monthly lightning to
c           daily lightning
c
            if(iday.ge.15.and.iday.le.45)then ! mid jan - mid feb
              month1=1
              month2=2
              xday=iday-15
            else if(iday.ge.46.and.iday.le.74)then ! mid feb - mid mar
              month1=2
              month2=3
              xday=iday-46
            else if(iday.ge.75.and.iday.le.105)then ! mid mar - mid apr
              month1=3
              month2=4
              xday=iday-75
            else if(iday.ge.106.and.iday.le.135)then ! mid apr - mid may
              month1=4
              month2=5
              xday=iday-106
            else if(iday.ge.136.and.iday.le.165)then ! mid may - mid june
              month1=5
              month2=6
              xday=iday-136
            else if(iday.ge.166.and.iday.le.196)then ! mid june - mid july
              month1=6
              month2=7
              xday=iday-166
            else if(iday.ge.197.and.iday.le.227)then ! mid july - mid aug
              month1=7
              month2=8
              xday=iday-197
            else if(iday.ge.228.and.iday.le.258)then ! mid aug - mid sep
              month1=8
              month2=9
              xday=iday-228
            else if(iday.ge.259.and.iday.le.288)then ! mid sep - mid oct
              month1=9
              month2=10
              xday=iday-259
            else if(iday.ge.289.and.iday.le.319)then ! mid oct - mid nov
              month1=10
              month2=11
              xday=iday-289
            else if(iday.ge.320.and.iday.le.349)then ! mid nov - mid dec
              month1=11
              month2=12
              xday=iday-320
            else if(iday.ge.350.or.iday.lt.14)then ! mid dec - mid jan
              month1=12
              month2=1
              xday=iday-350
              if(xday.lt.0)xday=iday
            endif
c
            lightng(i)=mlightnggat(i,month1)+(real(xday)/30.0)*
     &                 (mlightnggat(i,month2)-mlightnggat(i,month1))
c
            if (obswetf) then
              wetfracgrd(i)=wetfrac_mon(i,month1)+(real(xday)/30.0)*
     &                 (wetfrac_mon(i,month2)-wetfrac_mon(i,month1))
            endif !obswetf

          endif ! if(ctem_on)
c
855   continue
c
c     call canadian terrestrial ecosystem model which operates at a
c     daily time step, and uses daily accumulated values of variables
c     simulated by class.
c
      if (ctem_on) then
c
        call ctem ( fcancmxgat, fsnowacc_m,    sandgat,    claygat,
     2                      1,        nml,        iday,    radjgat,
     4          tcanoaccgat_m,  tcansacc_m, tbarcacc_m,tbarcsacc_m,
     5             tbargacc_m, tbargsacc_m, taaccgat_m,    dlzwgat,
     6             ancsvgac_m,  ancgvgac_m, rmlcsvga_m, rmlcgvga_m,
     7                zbtwgat, thliqcacc_m,thliqgacc_m,     deltat,
     8             uvaccgat_m,  vvaccgat_m,    lightng,prbfrhucgat,
     9            extnprobgat,   stdalngat,tbaraccgat_m,  popdon,
     a               nol2pfts, pfcancmxgat, nfcancmxgat,  lnduseon,
     b            thicecacc_m,     sdepgat,    spinfast,   todfrac,
     &                compete,  netrad_gat,  preacc_gat,
     &                 popdin,  dofire, dowetlands,obswetf, isndgat,
     &                faregat,      mosaic, WETFRACGRD, wetfrac_sgrd,
c    -------------- inputs used by ctem are above this line ---------
     c            stemmassgat, rootmassgat, litrmassgat, gleafmasgat,
     d            bleafmasgat, soilcmasgat,    ailcggat,    ailcgat,
     e               zolncgat,  rmatctemgat,   rmatcgat,  ailcbgat,
     f            flhrlossgat,  pandaysgat, lfstatusgat, grwtheffgat,
     g            lystmmasgat, lyrotmasgat, tymaxlaigat, vgbiomasgat,
     h            gavgltmsgat, gavgscmsgat, stmhrlosgat,     slaigat,
     i             bmasveggat, cmasvegcgat,  colddaysgat, rothrlosgat,
     j                fcangat,  alvsctmgat,   alirctmgat,  gavglaigat,
     &                  tcurm,    srpcuryr,     dftcuryr,  inibioclim,
     &                 tmonth,    anpcpcur,      anpecur,     gdd5cur,
     &               surmncur,    defmncur,     srplscur,    defctcur,
     &            geremortgat, intrmortgat,    lambdagat, lyglfmasgat,
     &            pftexistgat,      twarmm,       tcoldm,        gdd5,
     1                aridity,    srplsmon,     defctmon,    anndefct,
     2               annsrpls,      annpcp,  dry_season_length,
     &              burnvegfgat, pstemmassgat, pgleafmassgat,
c    -------------- inputs updated by ctem are above this line ------
     k                 nppgat,      nepgat, hetroresgat, autoresgat,
     l            soilcrespgat,       rmgat,       rggat,      nbpgat,
     m              litresgat,    socresgat,     gppgat, dstcemlsgat,
     n            litrfallgat,  humiftrsgat, veghghtgat, rootdpthgat,
     o                 rmlgat,      rmsgat,     rmrgat,  tltrleafgat,
     p            tltrstemgat, tltrrootgat, leaflitrgat, roottempgat,
     q             afrleafgat,  afrstemgat,  afrrootgat, wtstatusgat,
     r            ltstatusgat, burnfracgat, probfiregat, lucemcomgat,
     s            lucltringat, lucsocingat,   nppveggat, grclarea,
     t            dstcemls3gat,    paicgat,    slaicgat,
     u            emit_co2gat,  emit_cogat,  emit_ch4gat, emit_nmhcgat,
     v             emit_h2gat, emit_noxgat,  emit_n2ogat, emit_pm25gat,
     w            emit_tpmgat,  emit_tcgat,   emit_ocgat,   emit_bcgat,
     &               btermgat,    ltermgat,     mtermgat,
     &            ccgat,             mmgat,
     &          rmlvegaccgat,    rmsveggat,  rmrveggat,  rgveggat,
     &       vgbiomas_veggat, gppveggat,  nepveggat, nbpveggat,
     &        hetroresveggat, autoresveggat, litresveggat,
     &           soilcresveggat, nml, ilmos, jlmos, CH4WET1GAT,
     &          CH4WET2GAT, WETFDYNGAT, CH4DYN1GAT, CH4DYN2GAT)
c    ---------------- outputs are listed above this line ------------
c
      endif  !if(ctem_on)
c
c     reset mosaic accumulator arrays.
c
      do 655 i=1,nml
         uvaccgat_m(i)=0.0
655   continue
c
      if (ctem_on) then
        do 705 i = 1, nml
c
c         competitition related variables added by y. peng \\
          fsinacc_gat(i)=0.
          flinacc_gat(i)=0.
          flutacc_gat(i)=0.
          alswacc_gat(i)=0.
          allwacc_gat(i)=0.
          pregacc_gat(i)=0.
c         competitition related variables added by y. peng //
c
          fsnowacc_m(i)=0.0
          tcanoaccgat_out(i)=tcanoaccgat_m(i)
          tcanoaccgat_m(i)=0.0
c
          tcansacc_m(i)=0.0
          taaccgat_m(i)=0.0
          vvaccgat_m(i)=0.0
c
          do 715 j=1,ignd
             tbaraccgat_m(i,j)=0.0
             tbarcacc_m(i,j)=0.0
             tbarcsacc_m(i,j)=0.0
             tbargacc_m(i,j)=0.0
             tbargsacc_m(i,j)=0.0
             thliqcacc_m(i,j)=0.0
             thliqgacc_m(i,j)=0.0
             thicecacc_m(i,j)=0.0
715       continue
c
          do 716 j = 1, icc
            ancsvgac_m(i,j)=0.0
            ancgvgac_m(i,j)=0.0
            rmlcsvga_m(i,j)=0.0
            rmlcgvga_m(i,j)=0.0
716       continue
c
705     continue
      endif  ! if(ctem_on)
      endif  ! if(ncount.eq.nday)
C===================== CTEM ============================================ /
C
      CALL CLASSS (TBARROT,THLQROT,THICROT,TSFSROT,TPNDROT,
     1             ZPNDROT,TBASROT,ALBSROT,TSNOROT,RHOSROT,
     2             SNOROT, GTROT, TCANROT,RCANROT,SCANROT,
     3             GROROT, CMAIROT,TACROT, QACROT, WSNOROT,
     +             REFROT, BCSNROT,EMISROT,SALBROT,CSALROT,
     4             ILMOS,JLMOS,NML,NLAT,NTLD,NMOS,
     5             ILG,IGND,ICAN,ICAN+1,NBS,
     6             TBARGAT,THLQGAT,THICGAT,TSFSGAT,TPNDGAT,
     7             ZPNDGAT,TBASGAT,ALBSGAT,TSNOGAT,RHOSGAT,
     8             SNOGAT, GTGAT, TCANGAT,RCANGAT,SCANGAT,
     9             GROGAT, CMAIGAT,TACGAT, QACGAT, WSNOGAT,
     +             REFGAT, BCSNGAT,EMISGAT,SALBGAT,CSALGAT)

C
C    * SCATTER OPERATION ON DIAGNOSTIC VARIABLES SPLIT OUT OF
C    * CLASSS FOR CONSISTENCY WITH GCM APPLICATIONS.
C
      DO 380 K=1,NML
          CDHROT (ILMOS(K),JLMOS(K))=CDHGAT (K)
          CDMROT (ILMOS(K),JLMOS(K))=CDMGAT (K)
          HFSROT (ILMOS(K),JLMOS(K))=HFSGAT (K)
          TFXROT (ILMOS(K),JLMOS(K))=TFXGAT (K)
          QEVPROT(ILMOS(K),JLMOS(K))=QEVPGAT(K)
          QFSROT (ILMOS(K),JLMOS(K))=QFSGAT (K)
          QFXROT (ILMOS(K),JLMOS(K))=QFXGAT (K)
          PETROT (ILMOS(K),JLMOS(K))=PETGAT (K)
          GAROT  (ILMOS(K),JLMOS(K))=GAGAT  (K)
          EFROT  (ILMOS(K),JLMOS(K))=EFGAT  (K)
          QGROT  (ILMOS(K),JLMOS(K))=QGGAT  (K)
          ALVSROT(ILMOS(K),JLMOS(K))=ALVSGAT(K)
          ALIRROT(ILMOS(K),JLMOS(K))=ALIRGAT(K)
          SFCTROT(ILMOS(K),JLMOS(K))=SFCTGAT(K)
          SFCUROT(ILMOS(K),JLMOS(K))=SFCUGAT(K)
          SFCVROT(ILMOS(K),JLMOS(K))=SFCVGAT(K)
          SFCQROT(ILMOS(K),JLMOS(K))=SFCQGAT(K)
          FSNOROT(ILMOS(K),JLMOS(K))=FSNOGAT(K)
          FSGVROT(ILMOS(K),JLMOS(K))=FSGVGAT(K)
          FSGSROT(ILMOS(K),JLMOS(K))=FSGSGAT(K)
          FSGGROT(ILMOS(K),JLMOS(K))=FSGGGAT(K)
          FLGVROT(ILMOS(K),JLMOS(K))=FLGVGAT(K)
          FLGSROT(ILMOS(K),JLMOS(K))=FLGSGAT(K)
          FLGGROT(ILMOS(K),JLMOS(K))=FLGGGAT(K)
          HFSCROT(ILMOS(K),JLMOS(K))=HFSCGAT(K)
          HFSSROT(ILMOS(K),JLMOS(K))=HFSSGAT(K)
          HFSGROT(ILMOS(K),JLMOS(K))=HFSGGAT(K)
          HEVCROT(ILMOS(K),JLMOS(K))=HEVCGAT(K)
          HEVSROT(ILMOS(K),JLMOS(K))=HEVSGAT(K)
          HEVGROT(ILMOS(K),JLMOS(K))=HEVGGAT(K)
          HMFCROT(ILMOS(K),JLMOS(K))=HMFCGAT(K)
          HMFNROT(ILMOS(K),JLMOS(K))=HMFNGAT(K)
          HTCCROT(ILMOS(K),JLMOS(K))=HTCCGAT(K)
          HTCSROT(ILMOS(K),JLMOS(K))=HTCSGAT(K)
          PCFCROT(ILMOS(K),JLMOS(K))=PCFCGAT(K)
          PCLCROT(ILMOS(K),JLMOS(K))=PCLCGAT(K)
          PCPNROT(ILMOS(K),JLMOS(K))=PCPNGAT(K)
          PCPGROT(ILMOS(K),JLMOS(K))=PCPGGAT(K)
          QFGROT (ILMOS(K),JLMOS(K))=QFGGAT (K)
          QFNROT (ILMOS(K),JLMOS(K))=QFNGAT (K)
          QFCLROT(ILMOS(K),JLMOS(K))=QFCLGAT(K)
          QFCFROT(ILMOS(K),JLMOS(K))=QFCFGAT(K)
          ROFROT (ILMOS(K),JLMOS(K))=ROFGAT (K)
          ROFOROT(ILMOS(K),JLMOS(K))=ROFOGAT(K)
          ROFSROT(ILMOS(K),JLMOS(K))=ROFSGAT(K)
          ROFBROT(ILMOS(K),JLMOS(K))=ROFBGAT(K)
          TROFROT(ILMOS(K),JLMOS(K))=TROFGAT(K)
          TROOROT(ILMOS(K),JLMOS(K))=TROOGAT(K)
          TROSROT(ILMOS(K),JLMOS(K))=TROSGAT(K)
          TROBROT(ILMOS(K),JLMOS(K))=TROBGAT(K)
          ROFCROT(ILMOS(K),JLMOS(K))=ROFCGAT(K)
          ROFNROT(ILMOS(K),JLMOS(K))=ROFNGAT(K)
          ROVGROT(ILMOS(K),JLMOS(K))=ROVGGAT(K)
          WTRCROT(ILMOS(K),JLMOS(K))=WTRCGAT(K)
          WTRSROT(ILMOS(K),JLMOS(K))=WTRSGAT(K)
          WTRGROT(ILMOS(K),JLMOS(K))=WTRGGAT(K)
          DRROT  (ILMOS(K),JLMOS(K))=DRGAT  (K)
          WTABROT(ILMOS(K),JLMOS(K))=WTABGAT(K)
          ILMOROT(ILMOS(K),JLMOS(K))=ILMOGAT(K)
          UEROT  (ILMOS(K),JLMOS(K))=UEGAT(K)
          HBLROT (ILMOS(K),JLMOS(K))=HBLGAT(K)
380   CONTINUE
C
      DO 390 L=1,IGND
      DO 390 K=1,NML
          HMFGROT(ILMOS(K),JLMOS(K),L)=HMFGGAT(K,L)
          HTCROT (ILMOS(K),JLMOS(K),L)=HTCGAT (K,L)
          QFCROT (ILMOS(K),JLMOS(K),L)=QFCGAT (K,L)
          GFLXROT(ILMOS(K),JLMOS(K),L)=GFLXGAT(K,L)
390   CONTINUE
C
      DO 430 M=1,50
          DO 420 L=1,6
              DO 410 K=1,NML
                  ITCTROT(ILMOS(K),JLMOS(K),L,M)=ITCTGAT(K,L,M)
410           CONTINUE
420       CONTINUE
430   CONTINUE
C

C
C===================== CTEM ============================================ \
C
      call ctems2(fcancmxrow,rmatcrow,zolncrow,paicrow,
     1      ailcrow,     ailcgrow,    cmasvegcrow,  slaicrow,
     2      ailcgsrow,   fcancsrow,   fcancrow,     rmatctemrow,
     3      co2concrow,  co2i1cgrow,  co2i1csrow,   co2i2cgrow,
     4      co2i2csrow,  xdiffus,     slairow,      cfluxcgrow,
     5      cfluxcsrow,  ancsvegrow,  ancgvegrow,   rmlcsvegrow,
     6      rmlcgvegrow, canresrow,   SDEPROT,
     7      SANDROT,     CLAYROT,     ORGMROT,
     8      anvegrow,    rmlvegrow,   tcanoaccrow_m,tbaraccrow_m,
     9      uvaccrow_m,  vvaccrow_m,  mlightnggrd,  prbfrhucgrd,
     a      extnprobgrd, stdalngrd,   pfcancmxrow,  nfcancmxrow,
     b      stemmassrow, rootmassrow, litrmassrow,  gleafmasrow,
     c      bleafmasrow, soilcmasrow, ailcbrow,     flhrlossrow,
     d      pandaysrow,  lfstatusrow, grwtheffrow,  lystmmasrow,
     e      lyrotmasrow, tymaxlairow, vgbiomasrow,  gavgltmsrow,
     f      stmhrlosrow, bmasvegrow,  colddaysrow,  rothrlosrow,
     g      alvsctmrow,  alirctmrow,  gavglairow,   npprow,
     h      neprow,      hetroresrow, autoresrow,   soilcresprow,
     i      rmrow,       rgrow,       nbprow,       litresrow,
     j      socresrow,   gpprow,      dstcemlsrow,  litrfallrow,
     k      humiftrsrow, veghghtrow,  rootdpthrow,  rmlrow,
     l      rmsrow,      rmrrow,      tltrleafrow,  tltrstemrow,
     m      tltrrootrow, leaflitrrow, roottemprow,  afrleafrow,
     n      afrstemrow,  afrrootrow,  wtstatusrow,  ltstatusrow,
     o      burnfracrow, probfirerow, lucemcomrow,  lucltrinrow,
     p      lucsocinrow, nppvegrow,   dstcemls3row,
     q      FAREROT,     gavgscmsrow, tcanoaccrow_out,
     &      rmlvegaccrow, rmsvegrow,  rmrvegrow,    rgvegrow,
     &      vgbiomas_vegrow,gppvegrow,nepvegrow,ailcminrow,ailcmaxrow,
     &      FCANROT,      pftexistrow,
     &      emit_co2row,  emit_corow, emit_ch4row,  emit_nmhcrow,
     &      emit_h2row,   emit_noxrow,emit_n2orow,  emit_pm25row,
     &      emit_tpmrow,  emit_tcrow, emit_ocrow,   emit_bcrow,
     &      btermrow,     ltermrow,   mtermrow,
     &      nbpvegrow,   hetroresvegrow, autoresvegrow,litresvegrow,
     &      soilcresvegrow, burnvegfrow, pstemmassrow, pgleafmassrow,
!     &      WETFRACROW, WETFRAC_SROW,
     &      CH4WET1ROW, CH4WET2ROW,
     &      WETFDYNROW, CH4DYN1ROW, CH4DYN2ROW,
c    ----
     r      ilmos,       jlmos,       iwmos,        jwmos,
     s      nml,     fcancmxgat,  rmatcgat,    zolncgat,     paicgat,
     v      ailcgat,     ailcggat,    cmasvegcgat,  slaicgat,
     w      ailcgsgat,   fcancsgat,   fcancgat,     rmatctemgat,
     x      co2concgat,  co2i1cggat,  co2i1csgat,   co2i2cggat,
     y      co2i2csgat,  xdiffusgat,  slaigat,      cfluxcggat,
     z      cfluxcsgat,  ancsveggat,  ancgveggat,   rmlcsveggat,
     1      rmlcgveggat, canresgat,   sdepgat,
     2      sandgat,     claygat,     orgmgat,
     3      anveggat,    rmlveggat,   tcanoaccgat_m,tbaraccgat_m,
     4      uvaccgat_m,  vvaccgat_m,  mlightnggat,  prbfrhucgat,
     5      extnprobgat, stdalngat,   pfcancmxgat,  nfcancmxgat,
     6      stemmassgat, rootmassgat, litrmassgat,  gleafmasgat,
     7      bleafmasgat, soilcmasgat, ailcbgat,     flhrlossgat,
     8      pandaysgat,  lfstatusgat, grwtheffgat,  lystmmasgat,
     9      lyrotmasgat, tymaxlaigat, vgbiomasgat,  gavgltmsgat,
     a      stmhrlosgat, bmasveggat,  colddaysgat,  rothrlosgat,
     b      alvsctmgat,  alirctmgat,  gavglaigat,   nppgat,
     c      nepgat,      hetroresgat, autoresgat,   soilcrespgat,
     d      rmgat,       rggat,       nbpgat,       litresgat,
     e      socresgat,   gppgat,      dstcemlsgat,  litrfallgat,
     f      humiftrsgat, veghghtgat,  rootdpthgat,  rmlgat,
     g      rmsgat,      rmrgat,      tltrleafgat,  tltrstemgat,
     h      tltrrootgat, leaflitrgat, roottempgat,  afrleafgat,
     i      afrstemgat,  afrrootgat,  wtstatusgat,  ltstatusgat,
     j      burnfracgat, probfiregat, lucemcomgat,  lucltringat,
     k      lucsocingat, nppveggat,   dstcemls3gat,
     l      faregat,     gavgscmsgat, tcanoaccgat_out,
     &      rmlvegaccgat, rmsveggat,  rmrveggat,    rgveggat,
     &      vgbiomas_veggat,gppveggat,nepveggat,ailcmingat,ailcmaxgat,
     &      fcangat,      pftexistgat,
     &      emit_co2gat,  emit_cogat, emit_ch4gat,  emit_nmhcgat,
     &      emit_h2gat,   emit_noxgat,emit_n2ogat,  emit_pm25gat,
     &      emit_tpmgat,  emit_tcgat, emit_ocgat,   emit_bcgat,
     &      btermgat,     ltermgat,   mtermgat,
     &      nbpveggat, hetroresveggat, autoresveggat,litresveggat,
     &      soilcresveggat, burnvegfgat, pstemmassgat, pgleafmassgat,
!     &      WETFRACGAT, WETFRAC_SGAT,
     &      CH4WET1GAT, CH4WET2GAT,
     &      WETFDYNGAT, CH4DYN1GAT, CH4DYN2GAT)
c
C===================== CTEM ============================================ /
C
C=======================================================================
C     * WRITE FIELDS FROM CURRENT TIME STEP TO OUTPUT FILES.

6100  FORMAT(1X,I4,I5,9F8.2,2F8.3,F12.4,F8.2,2(A6,I2))
6200  FORMAT(1X,I4,I5,3(F8.2,2F6.3),F8.2,2F8.4,F8.2,F8.3,2(A6,I2))
6201  FORMAT(1X,I4,I5,5(F7.2,2F6.3),2(A6,I2))
6300  FORMAT(1X,I4,I5,3F9.2,F8.2,F10.2,E12.3,2F12.3,A6,I2)
6400  FORMAT(1X,I2,I3,I5,I6,9F8.2,2F7.3,E11.3,F8.2,F12.4,5F9.5,2(A6,I2))
6500  FORMAT(1X,I2,I3,I5,I6,3(F7.2,2F6.3),F8.2,2F8.4,F8.2,4F8.3,
     &       2F7.3,2(A6,I2))
6600  FORMAT(1X,I2,I3,I5,2F10.2,E12.3,F10.2,F8.2,F10.2,E12.3,2(A6,I2))
6501  FORMAT(1X,I2,I3,I5,I6,5(F7.2,2F6.3),2(A6,I2))
6601  FORMAT(1X,I2,I3,I5,I6,7(F7.2,2F6.3),10F9.4,2(A6,I2))
6700  FORMAT(1X,I2,I3,I5,I6,2X,12E11.4,2(A6,I2))
6800  FORMAT(1X,I2,I3,I5,I6,2X,22(F10.4,2X),2(A6,I2))
6900  FORMAT(1X,I2,I3,I5,I6,2X,18(E12.4,2X),2(A6,I2))
C
C===================== CTEM ============================================ \
c
c  fc,fg,fcs and fgs are one_dimensional in class subroutines
c  the transformations here to grid_cell mean fc_g,fg_g,fcs_g and fgs_g
c  are only applicable when nltest=1 (e.g., one grid cell)
c
      do i=1,nltest
        fc_g(i)=0.0
        fg_g(i)=0.0
        fcs_g(i)=0.0
        fgs_g(i)=0.0
        do m=1,nmtest
          fc_g(i)=fc_g(i)+fc(m)
          fg_g(i)=fg_g(i)+fg(m)
          fcs_g(i)=fcs_g(i)+fcs(m)
          fgs_g(i)=fgs_g(i)+fgs(m)
        enddo
      enddo
c
      if (.not. parallelrun) then ! stand alone mode, include half-hourly
c                                 ! output for class & ctem
C
C===================== CTEM =====================================/
      DO 450 I=1,NLTEST
C===================== CTEM =====================================\
c       initialization of various grid-averaged variables
        call resetgridavg(nltest)
!         fsstar_g    =0.0
!         flstar_g    =0.0
!         qh_g        =0.0
!         qe_g        =0.0
!         snomlt_g    =0.0
!         beg_g       =0.0
!         gtout_g     =0.0
! !        SNOROT_g(i) =0.0
! !        RHOSROT_g(i)=0.0
! !        WSNOROT_g(i)=0.0
!         altot_g     =0.0
! !        ROFROT_g(i) =0.0
!         tpn_g       =0.0
! !        ZPNDROT_g(i)=0.0
! c
! !         do j=1,ignd
! !          TBARROT_g(i,j)=0.0
! !          THLQROT_g(i,j)=0.0
! !          THICROT_g(i,j)=0.0
! !          GFLXROT_g(i,j)=0.0
! !         enddo
! c
!         tcn_g=0.0
! !        RCANROT_g(i) =0.0
! !        SCANROT_g(i) =0.0
!         tsn_g=0.0
!         zsn_g=0.0
! !         TROFROT_g(i)=0.0
!         TROOROT_g(i)=0.0
!         TROSROT_g(i)=0.0
!         TROBROT_g(i)=0.0
!         ROFOROT_g(i)=0.0
!         ROFSROT_g(i)=0.0
!         ROFBROT_g(i)=0.0
!         FSGVROT_g(i)=0.0
!         FSGSROT_g(i)=0.0
!         FSGGROT_g(i)=0.0
!         FLGVROT_g(i)=0.0
!         FLGSROT_g(i)=0.0
!         FLGGROT_g(i)=0.0
!         HFSCROT_g(i)=0.0
!         HFSSROT_g(i)=0.0
!         HFSGROT_g(i)=0.0
!         HEVCROT_g(i)=0.0
!         HEVSROT_g(i)=0.0
!         HEVGROT_g(i)=0.0
!         HMFCROT_g(i)=0.0
!         HMFNROT_g(i)=0.0
! c
!         do j=1,ignd
!          HMFGROT_g(i,j)=0.0
!          HTCROT_g(i,j)=0.0
!         enddo
! c
!         HTCCROT_g(i)=0.0
!         HTCSROT_g(i)=0.0
!         PCFCROT_g(i)=0.0
!         PCLCROT_g(i)=0.0
!         PCPNROT_g(i)=0.0
!         PCPGROT_g(i)=0.0
!         QFCFROT_g(i)=0.0
!         QFCLROT_g(i)=0.0
!         QFNROT_g(i)=0.0
!         QFGROT_g(i)=0.0
! c
!         do j=1,ignd
!          QFCROT_g(i,j)=0.0
!         enddo
! c
!         ROFCROT_g(i)=0.0
!         ROFNROT_g(i)=0.0
!         WTRCROT_g(i)=0.0
!         WTRSROT_g(i)=0.0
!         WTRGROT_g(i)=0.0
! c
!        if (ctem_on) then
!           do j=1,icc
!             anvegrow_g(i,j)=0.0
!             rmlvegrow_g(i,j)=0.0
!           enddo
!        endif   !if (ctem_on)
c
C===================== CTEM =====================================/
C

      ACTLYR=0.0
      FTABLE=0.0
      DO 440 J=1,IGND
          IF(ABS(TBARGAT(1,J)-TFREZ).LT.0.0001) THEN
              IF(ISNDGAT(1,J).GT.-3) THEN
                  ACTLYR=ACTLYR+(THLQGAT(1,J)/(THLQGAT(1,J)+
     1                THICGAT(1,J)))*DLZWGAT(1,J)
              ELSEIF(ISNDGAT(1,J).EQ.-3) THEN
                  ACTLYR=ACTLYR+DELZ(J)
              ENDIF
          ELSEIF(TBARGAT(1,J).GT.TFREZ) THEN
              ACTLYR=ACTLYR+DELZ(J)
          ENDIF
          IF(ABS(TBARGAT(1,J)-TFREZ).LT.0.0001) THEN
              IF(ISNDGAT(1,J).GT.-3) THEN
                  FTABLE=FTABLE+(THICGAT(1,J)/(THLQGAT(1,J)+
     1                THICGAT(1,J)-THMGAT(1,J)))*DLZWGAT(1,J)
              ELSE
                  FTABLE=FTABLE+DELZ(J)
              ENDIF
          ELSEIF(TBARGAT(1,J).LT.TFREZ) THEN
              FTABLE=FTABLE+DELZ(J)
          ENDIF
440   CONTINUE
C
      IF(IDAY.GE.182 .AND. IDAY.LE.243)  THEN
          ALAVG=ALAVG+ACTLYR
          NAL=NAL+1
          IF(ACTLYR.GT.ALMAX) ALMAX=ACTLYR
      ENDIF
C
      IF(IDAY.GE.1 .AND. IDAY.LE.59)   THEN
          FTAVG=FTAVG+FTABLE
          NFT=NFT+1
          IF(FTABLE.GT.FTMAX) FTMAX=FTABLE
      ENDIF

       DO 425 M=1,NMTEST
          IF(FSSROW(I).GT.0.0) THEN
C              ALTOT=(ALVSROT(I,M)+ALIRROT(I,M))/2.0
              ALTOT=(FSSROW(I)-FSGGGAT(1))/FSSROW(I)  !FLAG I adopt the runclass approach of using 1 for index here. JM Jul 2015.
          ELSE
              ALTOT=0.0
          ENDIF
          FSSTAR=FSSROW(I)*(1.0-ALTOT)
          FLSTAR=FDLROW(I)-SBC*GTROT(I,M)**4
          QH=HFSROT(I,M)
          QE=QEVPROT(I,M)
C          BEG=FSSTAR+FLSTAR-QH-QE !(commented out in runclass.fieldsite)
          BEG=GFLXGAT(1,1)  !FLAG!
C          USTARBS=UVROW(1)*SQRT(CDMROT(I,M)) !FLAG (commented out in runclass.fieldsite)
          SNOMLT=HMFNROT(I,M)
          IF(RHOSROT(I,M).GT.0.0) THEN
              ZSN=SNOROT(I,M)/RHOSROT(I,M)
          ELSE
              ZSN=0.0
          ENDIF
          IF(TCANROT(I,M).GT.0.01) THEN
              TCN=TCANROT(I,M)-TFREZ
          ELSE
              TCN=0.0
          ENDIF
          TSURF=FCS(I)*TSFSGAT(I,1)+FGS(I)*TSFSGAT(I,2)+
     1           FC(I)*TSFSGAT(I,3)+FG(I)*TSFSGAT(I,4)
C          IF(FSSROW(I).GT.0.0 .AND. (FCS(I)+FC(I)).GT.0.0) THEN
C          IF(FSSROW(I).GT.0.0) THEN
              NFS=NFS+1
              ITA=NINT(TAROW(I)-TFREZ)
              ITCAN=NINT(TCN)
              ITAC=NINT(TACGAT(I)-TFREZ)
              ITSCR=NINT(SFCTGAT(I)-TFREZ)
              ITS=NINT(TSURF-TFREZ)
C              ITD=ITS-ITA
              ITD=ITCAN-ITA
              ITD2=ITCAN-ITSCR
              ITD3=ITCAN-ITAC
              ITD4=ITAC-ITA
C              IF(ITA.GT.0.0) THEN
                  TAHIST(ITA+100)=TAHIST(ITA+100)+1.0
                  TCHIST(ITCAN+100)=TCHIST(ITCAN+100)+1.0
                  TSHIST(ITS+100)=TSHIST(ITS+100)+1.0
                  TACHIST(ITAC+100)=TACHIST(ITAC+100)+1.0
                  TDHIST(ITD+100)=TDHIST(ITD+100)+1.0
                  TD2HIST(ITD2+100)=TD2HIST(ITD2+100)+1.0
                  TD3HIST(ITD3+100)=TD3HIST(ITD3+100)+1.0
                  TD4HIST(ITD4+100)=TD4HIST(ITD4+100)+1.0
                  TSCRHIST(ITSCR+100)=TSCRHIST(ITSCR+100)+1.0
C              ENDIF
C          ENDIF     
          IF(FC(I).GT.0.1 .AND. RC(I).GT.1.0E5) NDRY=NDRY+1
!           IF((ITCAN-ITA).GE.10) THEN
!               WRITE(6,6070) IHOUR,IMIN,IDAY,IYEAR,FSSTAR,FLSTAR,QH,QE,
!      1                      BEG,TAROW(I)-TFREZ,TCN,TCN-(TAROW(I)-TFREZ),
!      2                      PAICAN(I),FSVF(I),UVROW(I),RC(I)
! 6070          FORMAT(2X,2I2,I4,I5,9F6.1,F6.3,F6.1,F8.1)
!           ENDIF
C
          IF(TSNOROT(I,M).GT.0.01) THEN
              TSN=TSNOROT(I,M)-TFREZ
          ELSE
              TSN=0.0
          ENDIF
          IF(TPNDROT(I,M).GT.0.01) THEN
              TPN=TPNDROT(I,M)-TFREZ
          ELSE
              TPN=0.0
          ENDIF
          GTOUT=GTROT(I,M)-TFREZ
          EVAPSUM=QFCFROT(I,M)+QFCLROT(I,M)+QFNROT(I,M)+QFGROT(I,M)+
     1                   QFCROT(I,M,1)+QFCROT(I,M,2)+QFCROT(I,M,3)
C
C===================== CTEM =====================================\
c         start writing output
c
          if ((iyear .ge. jhhsty) .and. (iyear .le. jhhendy)) then
           if ((iday .ge. jhhstd) .and. (iday .le. jhhendd)) then
C===================== CTEM =====================================/
          WRITE(64,6400) IHOUR,IMIN,IDAY,IYEAR,FSSTAR,FLSTAR,QH,QE,
     1                   SNOMLT,BEG,GTOUT,SNOROT(I,M),RHOSROT(I,M),
     2                   WSNOROT(I,M),ALTOT,ROFROT(I,M),
     3                   TPN,ZPNDROT(I,M),CDHROT(I,M),CDMROT(I,M),
     4                   SFCUROT(I,M),SFCVROT(I,M),UVROW(I),' TILE ',m
          IF(IGND.GT.3) THEN
C===================== CTEM =====================================\

              write(65,6500) ihour,imin,iday,iyear,(TBARROT(i,m,j)-
     1                   tfrez,THLQROT(i,m,j),THICROT(i,m,j),j=1,3),
     2                  tcn,RCANROT(i,m),SCANROT(i,m),tsn,zsn,
     3                   TCN-(TAROW(I)-TFREZ),TCANO(I)-TFREZ,
     4                   TACGAT(I)-TFREZ,ACTLYR,FTABLE,' TILE ',m
              write(66,6601) ihour,imin,iday,iyear,(TBARROT(i,m,j)-
     1                   tfrez,THLQROT(i,m,j),THICROT(i,m,j),j=4,10),
     2                   (GFLXROT(i,m,j),j=1,10),
     3                   ' TILE ',m
          else
              write(65,6500) ihour,imin,iday,iyear,(TBARROT(i,m,j)-
     1                   tfrez,THLQROT(i,m,j),THICROT(i,m,j),j=1,3),
     2                  tcn,RCANROT(i,m),SCANROT(i,m),tsn,zsn,
     3                   TCN-(TAROW(I)-TFREZ),TCANO(I)-TFREZ,
     4                   TACGAT(I)-TFREZ,ACTLYR,FTABLE,' TILE ',m

C===================== CTEM =====================================/
          ENDIF
C
          WRITE(67,6700) IHOUR,IMIN,IDAY,IYEAR,
     1                   TROFROT(I,M),TROOROT(I,M),TROSROT(I,M),
     2                   TROBROT(I,M),ROFROT(I,M),ROFOROT(I,M),
     3                   ROFSROT(I,M),ROFBROT(I,M),
     4                   FCS(M),FGS(M),FC(M),FG(M),' TILE ',M
          WRITE(68,6800) IHOUR,IMIN,IDAY,IYEAR,
     1                   FSGVROT(I,M),FSGSROT(I,M),FSGGROT(I,M),
     2                   FLGVROT(I,M),FLGSROT(I,M),FLGGROT(I,M),
     3                   HFSCROT(I,M),HFSSROT(I,M),HFSGROT(I,M),
     4                   HEVCROT(I,M),HEVSROT(I,M),HEVGROT(I,M),
     5                   HMFCROT(I,M),HMFNROT(I,M),
     6                   (HMFGROT(I,M,J),J=1,3),
     7                   HTCCROT(I,M),HTCSROT(I,M),
     8                   (HTCROT(I,M,J),J=1,3),' TILE ',M
          WRITE(69,6900) IHOUR,IMIN,IDAY,IYEAR,
     1                   PCFCROT(I,M),PCLCROT(I,M),PCPNROT(I,M),
     2                   PCPGROT(I,M),QFCFROT(I,M),QFCLROT(I,M),
     3                   QFNROT(I,M),QFGROT(I,M),(QFCROT(I,M,J),J=1,3),
     4                   ROFCROT(I,M),ROFNROT(I,M),ROFOROT(I,M),
     5                   ROFROT(I,M),WTRCROT(I,M),WTRSROT(I,M),
     6                   WTRGROT(I,M),' TILE ',M
C===================== CTEM =====================================\
C
         endif
        endif ! half hourly output loop.
c
c         Write half-hourly CTEM results to file *.CT01H
c
c         Net photosynthetic rates and leaf maintenance respiration for
c         each pft. however, if ctem_on then physyn subroutine
c         is using storage lai while actual lai is zero. if actual lai is
c         zero then we make anveg and rmlveg zero as well because these
c         are imaginary just like storage lai. note that anveg and rmlveg
c         are not passed to ctem. rather ancsveg, ancgveg, rmlcsveg, and
c         rmlcgveg are passed.
c
          if (ctem_on) then

            do 760 j = 1,icc
             if(ailcgrow(i,m,j).le.0.0) then
                anvegrow(i,m,j)=0.0
                rmlvegrow(i,m,j)=0.0
              else
                anvegrow(i,m,j)=anvegrow(i,m,j)
                rmlvegrow(i,m,j)=rmlvegrow(i,m,j)
              endif
760         continue
c
          if ((iyear .ge. jhhsty) .and. (iyear .le. jhhendy)) then
           if ((iday .ge. jhhstd) .and. (iday .le. jhhendd)) then

              write(71,7200)ihour,imin,iday,iyear,(anvegrow(i,m,j),
     1                    j=1,icc),(rmlvegrow(i,m,j),j=1,icc),
     2                    ' TILE ',m
            endif
           end if

           do j = 1,icc
              anvegrow_g(i,j)=anvegrow_g(i,j)+anvegrow(i,m,j)
     1                                        *FAREROT(i,m)
              rmlvegrow_g(i,j)=rmlvegrow_g(i,j)+rmlvegrow(i,m,j)
     1                                         *FAREROT(i,m)
            enddo

          endif   ! ctem_on

7200      format(1x,i2,1x,i2,i5,i5,9f11.3,9f11.3,2(a6,i2))
c
          fsstar_g    =fsstar_g + fsstar*FAREROT(i,m)
          flstar_g    =flstar_g + flstar*FAREROT(i,m)
          qh_g        =qh_g     + qh*FAREROT(i,m)
          qe_g        =qe_g     + qe*FAREROT(i,m)
          snomlt_g    =snomlt_g + snomlt*FAREROT(i,m)
          beg_g       =beg_g    + beg*FAREROT(i,m)
          gtout_g     =gtout_g  + gtout*FAREROT(i,m)
          tcn_g=tcn_g + tcn*FAREROT(i,m)
          tsn_g=tsn_g + tsn*FAREROT(i,m)
          zsn_g=zsn_g + zsn*FAREROT(i,m)
          altot_g     =altot_g   + altot*FAREROT(i,m)
          tpn_g       =tpn_g       + tpn*FAREROT(i,m)
c
          do j=1,ignd
            TBARROT_g(i,j)=TBARROT_g(i,j) + TBARROT(i,m,j)*FAREROT(i,m)
            THLQROT_g(i,j)=THLQROT_g(i,j) + THLQROT(i,m,j)*FAREROT(i,m)
            THICROT_g(i,j)=THICROT_g(i,j) + THICROT(i,m,j)*FAREROT(i,m)
            GFLXROT_g(i,j)=GFLXROT_g(i,j) + GFLXROT(i,m,j)*FAREROT(i,m)
            HMFGROT_g(i,j)=HMFGROT_g(i,j) + HMFGROT(i,m,j)*FAREROT(i,m)
            HTCROT_g(i,j)=HTCROT_g(i,j) + HTCROT(i,m,j)*FAREROT(i,m)
            QFCROT_g(i,j)=QFCROT_g(i,j) + QFCROT(i,m,j)*FAREROT(i,m)
          enddo
c
          ZPNDROT_g(i)=ZPNDROT_g(i) + ZPNDROT(i,m)*FAREROT(i,m)
          RHOSROT_g(i)=RHOSROT_g(i) + RHOSROT(i,m)*FAREROT(i,m)
          WSNOROT_g(i)=WSNOROT_g(i) + WSNOROT(i,m)*FAREROT(i,m)
          RCANROT_g(i)=RCANROT_g(i) + RCANROT(i,m)*FAREROT(i,m)
          SCANROT_g(i)=SCANROT_g(i) + SCANROT(i,m)*FAREROT(i,m)
          TROFROT_g(i)=TROFROT_g(i) + TROFROT(i,m)*FAREROT(i,m)
          TROOROT_g(i)=TROOROT_g(i) + TROOROT(i,m)*FAREROT(i,m)
          TROSROT_g(i)=TROSROT_g(i) + TROSROT(i,m)*FAREROT(i,m)
          TROBROT_g(i)=TROBROT_g(i) + TROBROT(i,m)*FAREROT(i,m)
          ROFOROT_g(i)=ROFOROT_g(i) + ROFOROT(i,m)*FAREROT(i,m)
          ROFSROT_g(i)=ROFSROT_g(i) + ROFSROT(i,m)*FAREROT(i,m)
          ROFBROT_g(i)=ROFBROT_g(i) + ROFBROT(i,m)*FAREROT(i,m)
          FSGVROT_g(i)=FSGVROT_g(i) + FSGVROT(i,m)*FAREROT(i,m)
          FSGSROT_g(i)=FSGSROT_g(i) + FSGSROT(i,m)*FAREROT(i,m)
          FSGGROT_g(i)=FSGGROT_g(i) + FSGGROT(i,m)*FAREROT(i,m)
          FLGVROT_g(i)=FLGVROT_g(i) + FLGVROT(i,m)*FAREROT(i,m)
          FLGSROT_g(i)=FLGSROT_g(i) + FLGSROT(i,m)*FAREROT(i,m)
          FLGGROT_g(i)=FLGGROT_g(i) + FLGGROT(i,m)*FAREROT(i,m)
          HFSCROT_g(i)=HFSCROT_g(i) + HFSCROT(i,m)*FAREROT(i,m)
          HFSSROT_g(i)=HFSSROT_g(i) + HFSSROT(i,m)*FAREROT(i,m)
          HFSGROT_g(i)=HFSGROT_g(i) + HFSGROT(i,m)*FAREROT(i,m)
          HEVCROT_g(i)=HEVCROT_g(i) + HEVCROT(i,m)*FAREROT(i,m)
          HEVSROT_g(i)=HEVSROT_g(i) + HEVSROT(i,m)*FAREROT(i,m)
          HEVGROT_g(i)=HEVGROT_g(i) + HEVGROT(i,m)*FAREROT(i,m)
          HMFCROT_g(i)=HMFCROT_g(i) + HMFCROT(i,m)*FAREROT(i,m)
          HMFNROT_g(i)=HMFNROT_g(i) + HMFNROT(i,m)*FAREROT(i,m)
          HTCCROT_g(i)=HTCCROT_g(i) + HTCCROT(i,m)*FAREROT(i,m)
          HTCSROT_g(i)=HTCSROT_g(i) + HTCSROT(i,m)*FAREROT(i,m)
          PCFCROT_g(i)=PCFCROT_g(i) + PCFCROT(i,m)*FAREROT(i,m)
          PCLCROT_g(i)=PCLCROT_g(i) + PCLCROT(i,m)*FAREROT(i,m)
          PCPNROT_g(i)=PCPNROT_g(i) + PCPNROT(i,m)*FAREROT(i,m)
          PCPGROT_g(i)=PCPGROT_g(i) + PCPGROT(i,m)*FAREROT(i,m)
          QFCFROT_g(i)=QFCFROT_g(i) + QFCFROT(i,m)*FAREROT(i,m)
          QFCLROT_g(i)=QFCLROT_g(i) + QFCLROT(i,m)*FAREROT(i,m)
          ROFCROT_g(i)=ROFCROT_g(i) + ROFCROT(i,m)*FAREROT(i,m)
          ROFNROT_g(i)=ROFNROT_g(i) + ROFNROT(i,m)*FAREROT(i,m)
          WTRCROT_g(i)=WTRCROT_g(i) + WTRCROT(i,m)*FAREROT(i,m)
          WTRSROT_g(i)=WTRSROT_g(i) + WTRSROT(i,m)*FAREROT(i,m)
          WTRGROT_g(i)=WTRGROT_g(i) + WTRGROT(i,m)*FAREROT(i,m)
          QFNROT_g(i) =QFNROT_g(i) + QFNROT(i,m)*FAREROT(i,m)
          QFGROT_g(i) =QFGROT_g(i) + QFGROT(i,m)*FAREROT(i,m)
          ROFROT_g(i) =ROFROT_g(i) + ROFROT(i,m)*FAREROT(i,m)
          SNOROT_g(i) =SNOROT_g(i) + SNOROT(i,m)*FAREROT(i,m)
          CDHROT_g(i) =CDHROT_g(i) + CDHROT(i,m)*FAREROT(i,m)
          CDMROT_g(i) =CDMROT_g(i) + CDMROT(i,m)*FAREROT(i,m)
          SFCUROT_g(i) =SFCUROT_g(i) + SFCUROT(i,m)*FAREROT(i,m)
          SFCVROT_g(i) =SFCVROT_g(i) + SFCVROT(i,m)*FAREROT(i,m)
C
C======================== CTEM =====================================/
425    CONTINUE
C===================== CTEM =====================================\
C      WRITE CTEM OUTPUT FILES
C
      if ((iyear .ge. jhhsty) .and. (iyear .le. jhhendy)) then
       if ((iday .ge. jhhstd) .and. (iday .le. jhhendd)) then

       IF (CTEM_ON) THEN
           WRITE(711,7200)IHOUR,IMIN,IDAY,IYEAR,(ANVEGROW_G(I,J),
     1                 J=1,ICC),(RMLVEGROW_G(I,J),J=1,ICC)
       ENDIF !CTEM_ON

       WRITE(641,6400) IHOUR,IMIN,IDAY,IYEAR,FSSTAR_G,FLSTAR_G,QH_G,
     1      QE_G,SNOMLT_G,BEG_G,GTOUT_G,SNOROT_G(I),RHOSROT_G(I),
     2                   WSNOROT_G(I),ALTOT_G,ROFROT_G(I),
     3                   TPN_G,ZPNDROT_G(I),CDHROT_G(I),CDMROT_G(I),
     4                   SFCUROT_G(I),SFCVROT_G(I),UVROW(I)
         WRITE(651,6500) IHOUR,IMIN,IDAY,IYEAR,(TBARROT_G(I,J)-
     1                   TFREZ,THLQROT_G(I,J),THICROT_G(I,J),J=1,3),
     2                   TCN_G,RCANROT_G(I),SCANROT_G(I),TSN_G,ZSN_G,
     3                   TCN_G-(TAROW(I)-TFREZ),TCANO(I)-TFREZ,
     4                   TACGAT(I)-TFREZ,ACTLYR,FTABLE
C
         IF(IGND.GT.3) THEN
          WRITE(661,6601) IHOUR,IMIN,IDAY,IYEAR,(TBARROT_G(I,J)-
     1                   TFREZ,THLQROT_G(I,J),THICROT_G(I,J),J=4,10),
     2                   (GFLXROT_G(I,J),J=1,10)
         ELSE
          WRITE(661,6600) IHOUR,IMIN,IDAY,FSSROW(I),FDLROW(I),PREROW(I),
     1                   TAROW(I)-TFREZ,UVROW(I),PRESROW(I),QAROW(I)
         ENDIF
C
         WRITE(671,6700) IHOUR,IMIN,IDAY,IYEAR,
     &                   TROFROT_G(I),TROOROT_G(I),TROSROT_G(I),
     1                   TROBROT_G(I),ROFROT_G(I),ROFOROT_G(I),
     2                   ROFSROT_G(I),ROFBROT_G(I),
     3                   FCS_G(I),FGS_G(I),FC_G(I),FG_G(I)
         WRITE(681,6800) IHOUR,IMIN,IDAY,IYEAR,
     &                   FSGVROT_G(I),FSGSROT_G(I),FSGGROT_G(I),
     1                   FLGVROT_G(I),FLGSROT_G(I),FLGGROT_G(I),
     2                   HFSCROT_G(I),HFSSROT_G(I),HFSGROT_G(I),
     3                   HEVCROT_G(I),HEVSROT_G(I),HEVGROT_G(I),
     4                   HMFCROT_G(I),HMFNROT_G(I),
     5                   (HMFGROT_G(I,J),J=1,3),
     6                   HTCCROT_G(I),HTCSROT_G(I),
     7                   (HTCROT_G(I,J),J=1,3)
         WRITE(691,6900) IHOUR,IMIN,IDAY,IYEAR,
     &                   PCFCROT_G(I),PCLCROT_G(I),PCPNROT_G(I),
     1                   PCPGROT_G(I),QFCFROT_G(I),QFCLROT_G(I),
     2                   QFNROT_G(I),QFGROT_G(I),(QFCROT_G(I,J),J=1,3),
     3                   ROFCROT_G(I),ROFNROT_G(I),ROFOROT_G(I),
     4                   ROFROT_G(I),WTRCROT_G(I),WTRSROT_G(I),
     5                   WTRGROT_G(I)
C
        endif
       ENDIF ! if write half-hourly
C===================== CTEM =====================================/
450   CONTINUE
C
C===================== CTEM =====================================\

      endif ! not parallelrun
C===================== CTEM =====================================/
C
C=======================================================================
C     * CALCULATE GRID CELL AVERAGE DIAGNOSTIC FIELDS.
C
C===================== CTEM =====================================\

      if(.not.parallelrun) then ! stand alone mode, includes
c                               ! diagnostic fields
C===================== CTEM =====================================/
C
      DO 525 I=1,NLTEST
          CDHROW(I)=0.
          CDMROW(I)=0.
          HFSROW(I)=0.
          TFXROW(I)=0.
          QEVPROW(I)=0.
          QFSROW(I)=0.
          QFXROW(I)=0.
          PETROW(I)=0.
          GAROW(I)=0.
          EFROW(I)=0.
          GTROW(I)=0.
          QGROW(I)=0.
          ALVSROW(I)=0.
          ALIRROW(I)=0.
          SFCTROW(I)=0.
          SFCUROW(I)=0.
          SFCVROW(I)=0.
          SFCQROW(I)=0.
          SFRHROW(I)=0.
          FSNOROW(I)=0.
          FSGVROW(I)=0.
          FSGSROW(I)=0.
          FSGGROW(I)=0.
          FLGVROW(I)=0.
          FLGSROW(I)=0.
          FLGGROW(I)=0.
          HFSCROW(I)=0.
          HFSSROW(I)=0.
          HFSGROW(I)=0.
          HEVCROW(I)=0.
          HEVSROW(I)=0.
          HEVGROW(I)=0.
          HMFCROW(I)=0.
          HMFNROW(I)=0.
          HTCCROW(I)=0.
          HTCSROW(I)=0.
          PCFCROW(I)=0.
          PCLCROW(I)=0.
          PCPNROW(I)=0.
          PCPGROW(I)=0.
          QFGROW(I)=0.
          QFNROW(I)=0.
          QFCLROW(I)=0.
          QFCFROW(I)=0.
          ROFROW(I)=0.
          ROFOROW(I)=0.
          ROFSROW(I)=0.
          ROFBROW(I)=0.
          ROFCROW(I)=0.
          ROFNROW(I)=0.
          ROVGROW(I)=0.
          WTRCROW(I)=0.
          WTRSROW(I)=0.
          WTRGROW(I)=0.
          DRROW(I)=0.
          WTABROW(I)=0.
          ILMOROW(I)=0.
          UEROW(I)=0.
          HBLROW(I)=0.
          DO 500 J=1,IGND
              HMFGROW(I,J)=0.
              HTCROW(I,J)=0.
              QFCROW(I,J)=0.
              GFLXROW(I,J)=0.
500       CONTINUE
525   CONTINUE
C
      DO 600 I=1,NLTEST
      DO 575 M=1,NMTEST
          CDHROW(I)=CDHROW(I)+CDHROT(I,M)*FAREROT(I,M)
          CDMROW(I)=CDMROW(I)+CDMROT(I,M)*FAREROT(I,M)
          HFSROW(I)=HFSROW(I)+HFSROT(I,M)*FAREROT(I,M)
          TFXROW(I)=TFXROW(I)+TFXROT(I,M)*FAREROT(I,M)
          QEVPROW(I)=QEVPROW(I)+QEVPROT(I,M)*FAREROT(I,M)
          QFSROW(I)=QFSROW(I)+QFSROT(I,M)*FAREROT(I,M)
          QFXROW(I)=QFXROW(I)+QFXROT(I,M)*FAREROT(I,M)
          PETROW(I)=PETROW(I)+PETROT(I,M)*FAREROT(I,M)
          GAROW(I)=GAROW(I)+GAROT(I,M)*FAREROT(I,M)
          EFROW(I)=EFROW(I)+EFROT(I,M)*FAREROT(I,M)
          GTROW(I)=GTROW(I)+GTROT(I,M)*FAREROT(I,M)
          QGROW(I)=QGROW(I)+QGROT(I,M)*FAREROT(I,M)
          ALVSROW(I)=ALVSROW(I)+ALVSROT(I,M)*FAREROT(I,M)
          ALIRROW(I)=ALIRROW(I)+ALIRROT(I,M)*FAREROT(I,M)
          SFCTROW(I)=SFCTROW(I)+SFCTROT(I,M)*FAREROT(I,M)
          SFCUROW(I)=SFCUROW(I)+SFCUROT(I,M)*FAREROT(I,M)
          SFCVROW(I)=SFCVROW(I)+SFCVROT(I,M)*FAREROT(I,M)
          SFCQROW(I)=SFCQROW(I)+SFCQROT(I,M)*FAREROT(I,M)
          SFRHROW(I)=SFRHROW(I)+SFRHROT(I,M)*FAREROT(I,M)
          FSNOROW(I)=FSNOROW(I)+FSNOROT(I,M)*FAREROT(I,M)
          FSGVROW(I)=FSGVROW(I)+FSGVROT(I,M)*FAREROT(I,M)
          FSGSROW(I)=FSGSROW(I)+FSGSROT(I,M)*FAREROT(I,M)
          FSGGROW(I)=FSGGROW(I)+FSGGROT(I,M)*FAREROT(I,M)
          FLGVROW(I)=FLGVROW(I)+FLGVROT(I,M)*FAREROT(I,M)
          FLGSROW(I)=FLGSROW(I)+FLGSROT(I,M)*FAREROT(I,M)
          FLGGROW(I)=FLGGROW(I)+FLGGROT(I,M)*FAREROT(I,M)
          HFSCROW(I)=HFSCROW(I)+HFSCROT(I,M)*FAREROT(I,M)
          HFSSROW(I)=HFSSROW(I)+HFSSROT(I,M)*FAREROT(I,M)
          HFSGROW(I)=HFSGROW(I)+HFSGROT(I,M)*FAREROT(I,M)
          HEVCROW(I)=HEVCROW(I)+HEVCROT(I,M)*FAREROT(I,M)
          HEVSROW(I)=HEVSROW(I)+HEVSROT(I,M)*FAREROT(I,M)
          HEVGROW(I)=HEVGROW(I)+HEVGROT(I,M)*FAREROT(I,M)
          HMFCROW(I)=HMFCROW(I)+HMFCROT(I,M)*FAREROT(I,M)
          HMFNROW(I)=HMFNROW(I)+HMFNROT(I,M)*FAREROT(I,M)
          HTCCROW(I)=HTCCROW(I)+HTCCROT(I,M)*FAREROT(I,M)
          HTCSROW(I)=HTCSROW(I)+HTCSROT(I,M)*FAREROT(I,M)
          PCFCROW(I)=PCFCROW(I)+PCFCROT(I,M)*FAREROT(I,M)
          PCLCROW(I)=PCLCROW(I)+PCLCROT(I,M)*FAREROT(I,M)
          PCPNROW(I)=PCPNROW(I)+PCPNROT(I,M)*FAREROT(I,M)
          PCPGROW(I)=PCPGROW(I)+PCPGROT(I,M)*FAREROT(I,M)
          QFGROW(I)=QFGROW(I)+QFGROT(I,M)*FAREROT(I,M)
          QFNROW(I)=QFNROW(I)+QFNROT(I,M)*FAREROT(I,M)
          QFCLROW(I)=QFCLROW(I)+QFCLROT(I,M)*FAREROT(I,M)
          QFCFROW(I)=QFCFROW(I)+QFCFROT(I,M)*FAREROT(I,M)
          ROFROW(I)=ROFROW(I)+ROFROT(I,M)*FAREROT(I,M)
          ROFOROW(I)=ROFOROW(I)+ROFOROT(I,M)*FAREROT(I,M)
          ROFSROW(I)=ROFSROW(I)+ROFSROT(I,M)*FAREROT(I,M)
          ROFBROW(I)=ROFBROW(I)+ROFBROT(I,M)*FAREROT(I,M)
          ROFCROW(I)=ROFCROW(I)+ROFCROT(I,M)*FAREROT(I,M)
          ROFNROW(I)=ROFNROW(I)+ROFNROT(I,M)*FAREROT(I,M)
          ROVGROW(I)=ROVGROW(I)+ROVGROT(I,M)*FAREROT(I,M)
          WTRCROW(I)=WTRCROW(I)+WTRCROT(I,M)*FAREROT(I,M)
          WTRSROW(I)=WTRSROW(I)+WTRSROT(I,M)*FAREROT(I,M)
          WTRGROW(I)=WTRGROW(I)+WTRGROT(I,M)*FAREROT(I,M)
          DRROW(I)=DRROW(I)+DRROT(I,M)*FAREROT(I,M)
          WTABROW(I)=WTABROW(I)+WTABROT(I,M)*FAREROT(I,M)
          ILMOROW(I)=ILMOROW(I)+ILMOROT(I,M)*FAREROT(I,M)
          UEROW(I)=UEROW(I)+UEROT(I,M)*FAREROT(I,M)
          HBLROW(I)=HBLROW(I)+HBLROT(I,M)*FAREROT(I,M)
          DO 550 J=1,IGND
              HMFGROW(I,J)=HMFGROW(I,J)+HMFGROT(I,M,J)*FAREROT(I,M)
              HTCROW(I,J)=HTCROW(I,J)+HTCROT(I,M,J)*FAREROT(I,M)
              QFCROW(I,J)=QFCROW(I,J)+QFCROT(I,M,J)*FAREROT(I,M)
              GFLXROW(I,J)=GFLXROW(I,J)+GFLXROT(I,M,J)*FAREROT(I,M)
550       CONTINUE
575   CONTINUE
600   CONTINUE
C
C===================== CTEM =====================================\

      endif ! not parallelrun, for diagnostic fields
c
      if(.not.parallelrun) then ! stand alone mode, includes daily output for class
C===================== CTEM =====================================/
C
C     * ACCUMULATE OUTPUT DATA FOR DIURNALLY AVERAGED FIELDS. BOTH GRID
C       MEAN AND MOSAIC MEAN
C
      DO 675 I=1,NLTEST
      DO 650 M=1,NMTEST
          PREACC(I)=PREACC(I)+PREROW(I)*FAREROT(I,M)*DELT
          GTACC(I)=GTACC(I)+GTROT(I,M)*FAREROT(I,M)
          QEVPACC(I)=QEVPACC(I)+QEVPROT(I,M)*FAREROT(I,M)
          EVAPACC(I)=EVAPACC(I)+QFSROT(I,M)*FAREROT(I,M)*DELT
          HFSACC(I)=HFSACC(I)+HFSROT(I,M)*FAREROT(I,M)
          HMFNACC(I)=HMFNACC(I)+HMFNROT(I,M)*FAREROT(I,M)
          ROFACC(I)=ROFACC(I)+ROFROT(I,M)*FAREROT(I,M)*DELT
          OVRACC(I)=OVRACC(I)+ROFOROT(I,M)*FAREROT(I,M)*DELT
          WTBLACC(I)=WTBLACC(I)+WTABROT(I,M)*FAREROT(I,M)
          DO 625 J=1,IGND
              TBARACC(I,J)=TBARACC(I,J)+TBARROT(I,M,J)*FAREROT(I,M)
              THLQACC(I,J)=THLQACC(I,J)+THLQROT(I,M,J)*FAREROT(I,M)
              THICACC(I,J)=THICACC(I,J)+THICROT(I,M,J)*FAREROT(I,M)
              THALACC(I,J)=THALACC(I,J)+(THLQROT(I,M,J)+THICROT(I,M,J))
     1                    *FAREROT(I,M)
625       CONTINUE
          ALVSACC(I)=ALVSACC(I)+ALVSROT(I,M)*FAREROT(I,M)*FSVHROW(I)
          ALIRACC(I)=ALIRACC(I)+ALIRROT(I,M)*FAREROT(I,M)*FSIHROW(I)
          IF(SNOROT(I,M).GT.0.0) THEN
              RHOSACC(I)=RHOSACC(I)+RHOSROT(I,M)*FAREROT(I,M)
              TSNOACC(I)=TSNOACC(I)+TSNOROT(I,M)*FAREROT(I,M)
              WSNOACC(I)=WSNOACC(I)+WSNOROT(I,M)*FAREROT(I,M)
              SNOARE(I)=SNOARE(I)+FAREROT(I,M)
          ENDIF
          IF(TCANROT(I,M).GT.0.5) THEN
              TCANACC(I)=TCANACC(I)+TCANROT(I,M)*FAREROT(I,M)
              CANARE(I)=CANARE(I)+FAREROT(I,M)
          ENDIF
          SNOACC(I)=SNOACC(I)+SNOROT(I,M)*FAREROT(I,M)
          RCANACC(I)=RCANACC(I)+RCANROT(I,M)*FAREROT(I,M)
          SCANACC(I)=SCANACC(I)+SCANROT(I,M)*FAREROT(I,M)
          GROACC(I)=GROACC(I)+GROROT(I,M)*FAREROT(I,M)
          FSINACC(I)=FSINACC(I)+FSSROW(I)*FAREROT(I,M)
          FLINACC(I)=FLINACC(I)+FDLROW(I)*FAREROT(I,M)
          FLUTACC(I)=FLUTACC(I)+SBC*GTROT(I,M)**4*FAREROT(I,M)
          TAACC(I)=TAACC(I)+TAROW(I)*FAREROT(I,M)
          UVACC(I)=UVACC(I)+UVROW(I)*FAREROT(I,M)
          PRESACC(I)=PRESACC(I)+PRESROW(I)*FAREROT(I,M)
          QAACC(I)=QAACC(I)+QAROW(I)*FAREROT(I,M)
650   CONTINUE
675   CONTINUE
C
C     * CALCULATE AND PRINT DAILY AVERAGES.
C
      IF(NCOUNT.EQ.NDAY) THEN

      DO 800 I=1,NLTEST
          PREACC(I)=PREACC(I)
          GTACC(I)=GTACC(I)/REAL(NDAY)
          QEVPACC(I)=QEVPACC(I)/REAL(NDAY)
          EVAPACC(I)=EVAPACC(I)
          HFSACC(I)=HFSACC(I)/REAL(NDAY)
          HMFNACC(I)=HMFNACC(I)/REAL(NDAY)
          ROFACC(I)=ROFACC(I)
          OVRACC(I)=OVRACC(I)
          WTBLACC(I)=WTBLACC(I)/REAL(NDAY)
          DO 725 J=1,IGND
              TBARACC(I,J)=TBARACC(I,J)/REAL(NDAY)
              THLQACC(I,J)=THLQACC(I,J)/REAL(NDAY)
              THICACC(I,J)=THICACC(I,J)/REAL(NDAY)
              THALACC(I,J)=THALACC(I,J)/REAL(NDAY)
725       CONTINUE
          IF(FSINACC(I).GT.0.0) THEN
              ALVSACC(I)=ALVSACC(I)/(FSINACC(I)*0.5)
              ALIRACC(I)=ALIRACC(I)/(FSINACC(I)*0.5)
          ELSE
              ALVSACC(I)=0.0
              ALIRACC(I)=0.0
          ENDIF
          IF(SNOARE(I).GT.0.0) THEN
              RHOSACC(I)=RHOSACC(I)/SNOARE(I)
              TSNOACC(I)=TSNOACC(I)/SNOARE(I)
              WSNOACC(I)=WSNOACC(I)/SNOARE(I)
          ENDIF
          IF(CANARE(I).GT.0.0) THEN
              TCANACC(I)=TCANACC(I)/CANARE(I)
          ENDIF
          SNOACC(I)=SNOACC(I)/REAL(NDAY)
          RCANACC(I)=RCANACC(I)/REAL(NDAY)
          SCANACC(I)=SCANACC(I)/REAL(NDAY)
          GROACC(I)=GROACC(I)/REAL(NDAY)
          FSINACC(I)=FSINACC(I)/REAL(NDAY)
          FLINACC(I)=FLINACC(I)/REAL(NDAY)
          FLUTACC(I)=FLUTACC(I)/REAL(NDAY)
          TAACC(I)=TAACC(I)/REAL(NDAY)
          UVACC(I)=UVACC(I)/REAL(NDAY)
          PRESACC(I)=PRESACC(I)/REAL(NDAY)
          QAACC(I)=QAACC(I)/REAL(NDAY)

              ALTOT=(ALVSACC(I)+ALIRACC(I))/2.0
              FSSTAR=FSINACC(I)*(1.-ALTOT)
              FLSTAR=FLINACC(I)-FLUTACC(I)
              QH=HFSACC(I)
              QE=QEVPACC(I)
              BEG=FSSTAR+FLSTAR-QH-QE
              SNOMLT=HMFNACC(I)
              IF(RHOSACC(I).GT.0.0) THEN
                  ZSN=SNOACC(I)/RHOSACC(I)
              ELSE
                  ZSN=0.0
              ENDIF
              IF(TCANACC(I).GT.0.01) THEN
                  TCN=TCANACC(I)-TFREZ
              ELSE
                  TCN=0.0
              ENDIF
              IF(TSNOACC(I).GT.0.01) THEN
                  TSN=TSNOACC(I)-TFREZ
              ELSE
                  TSN=0.0
              ENDIF
              GTOUT=GTACC(I)-TFREZ
C
             if ((iyear .ge. jdsty) .and. (iyear .le. jdendy)) then
              if ((iday .ge. jdstd) .and. (iday .le. jdendd)) then

              WRITE(61,6100) IDAY,IYEAR,FSSTAR,FLSTAR,QH,QE,SNOMLT,
     1                       BEG,GTOUT,SNOACC(I),RHOSACC(I),
     2                       WSNOACC(I),ALTOT,ROFACC(I),CUMSNO
              IF(IGND.GT.3) THEN
                  WRITE(62,6201) IDAY,IYEAR,(TBARACC(I,J)-TFREZ,
     1                       THLQACC(I,J),THICACC(I,J),J=1,5)
                  WRITE(63,6201) IDAY,IYEAR,(TBARACC(I,J)-TFREZ,
     1                       THLQACC(I,J),THICACC(I,J),J=6,10)
              ELSE
                  WRITE(62,6200) IDAY,IYEAR,(TBARACC(I,J)-TFREZ,
     1                       THLQACC(I,J),THICACC(I,J),J=1,3),
     2                       TCN,RCANACC(I),SCANACC(I),TSN,ZSN
                  WRITE(63,6300) IDAY,IYEAR,FSINACC(I),FLINACC(I),
     1                       TAACC(I)-TFREZ,UVACC(I),PRESACC(I),
     2                       QAACC(I),PREACC(I),EVAPACC(I)
              ENDIF
             endif
            ENDIF
C
C     * RESET ACCUMULATOR ARRAYS.
C
          PREACC(I)=0.
          GTACC(I)=0.
          QEVPACC(I)=0.
          HFSACC(I)=0.
          HMFNACC(I)=0.
          ROFACC(I)=0.
          SNOACC(I)=0.
          CANARE(I)=0.
          SNOARE(I)=0.
          OVRACC(I)=0.
          WTBLACC(I)=0.
          DO 750 J=1,IGND
              TBARACC(I,J)=0.
              THLQACC(I,J)=0.
              THICACC(I,J)=0.
              THALACC(I,J)=0.
750       CONTINUE
          ALVSACC(I)=0.
          ALIRACC(I)=0.
          RHOSACC(I)=0.
          TSNOACC(I)=0.
          WSNOACC(I)=0.
          TCANACC(I)=0.
          RCANACC(I)=0.
          SCANACC(I)=0.
          GROACC(I)=0.
          FSINACC(I)=0.
          FLINACC(I)=0.
          TAACC(I)=0.
          UVACC(I)=0.
          PRESACC(I)=0.
          QAACC(I)=0.
          EVAPACC(I)=0.
          FLUTACC(I)=0.
800   CONTINUE

      ENDIF ! IF(NCOUNT.EQ.NDAY)

C===================== CTEM =====================================\
C
C     CALCULATE AND PRINT MOSAIC DAILY AVERAGES.
C
!       start -> FLAG JM
      DO 676 I=1,NLTEST
      DO 658 M=1,NMTEST
          PREACC_M(I,M)=PREACC_M(I,M)+PREROW(I)*DELT
          GTACC_M(I,M)=GTACC_M(I,M)+GTROT(I,M)
          QEVPACC_M(I,M)=QEVPACC_M(I,M)+QEVPROT(I,M)
          EVAPACC_M(I,M)=EVAPACC_M(I,M)+QFSROT(I,M)*DELT
          HFSACC_M(I,M)=HFSACC_M(I,M)+HFSROT(I,M)
          HMFNACC_M(I,M)=HMFNACC_M(I,M)+HMFNROT(I,M)
          ROFACC_M(I,M)=ROFACC_M(I,M)+ROFROT(I,M)*DELT
          OVRACC_M(I,M)=OVRACC_M(I,M)+ROFOROT(I,M)*DELT
          WTBLACC_M(I,M)=WTBLACC_M(I,M)+WTABROT(I,M)
          DO 626 J=1,IGND
              TBARACC_M(I,M,J)=TBARACC_M(I,M,J)+TBARROT(I,M,J)
              THLQACC_M(I,M,J)=THLQACC_M(I,M,J)+THLQROT(I,M,J)
              THICACC_M(I,M,J)=THICACC_M(I,M,J)+THICROT(I,M,J)
              THALACC_M(I,M,J)=THALACC_M(I,M,J)+(THLQROT(I,M,J)+
     1           THICROT(I,M,J))
626       CONTINUE
          ALVSACC_M(I,M)=ALVSACC_M(I,M)+ALVSROT(I,M)*FSVHROW(I)
          ALIRACC_M(I,M)=ALIRACC_M(I,M)+ALIRROT(I,M)*FSIHROW(I)
          IF(SNOROT(I,M).GT.0.0) THEN
              RHOSACC_M(I,M)=RHOSACC_M(I,M)+RHOSROT(I,M)
              TSNOACC_M(I,M)=TSNOACC_M(I,M)+TSNOROT(I,M)
              WSNOACC_M(I,M)=WSNOACC_M(I,M)+WSNOROT(I,M)
              SNOARE_M(I,M) = SNOARE_M(I,M) + 1.0 !FLAG test.
          ENDIF
          IF(TCANROT(I,M).GT.0.5) THEN
              TCANACC_M(I,M)=TCANACC_M(I,M)+TCANROT(I,M)
C              CANARE(I)=CANARE(I)+FAREROT(I,M)
          ENDIF
          SNOACC_M(I,M)=SNOACC_M(I,M)+SNOROT(I,M)
          RCANACC_M(I,M)=RCANACC_M(I,M)+RCANROT(I,M)
          SCANACC_M(I,M)=SCANACC_M(I,M)+SCANROT(I,M)
          GROACC_M(I,M)=GROACC_M(I,M)+GROROT(I,M)
          FSINACC_M(I,M)=FSINACC_M(I,M)+FSSROW(I)
          FLINACC_M(I,M)=FLINACC_M(I,M)+FDLROW(I)
          FLUTACC_M(I,M)=FLUTACC_M(I,M)+SBC*GTROT(I,M)**4
          TAACC_M(I,M)=TAACC_M(I,M)+TAROW(I)
          UVACC_M(I,M)=UVACC_M(I,M)+UVROW(I)
          PRESACC_M(I,M)=PRESACC_M(I,M)+PRESROW(I)
          QAACC_M(I,M)=QAACC_M(I,M)+QAROW(I)
658   CONTINUE
676   CONTINUE
C
C     CALCULATE AND PRINT DAILY AVERAGES.
C
      IF(NCOUNT.EQ.NDAY) THEN

      DO 808 I=1,NLTEST
        DO 809 M=1,NMTEST
          PREACC_M(I,M)=PREACC_M(I,M)     !became [kg m-2 day-1] instead of [kg m-2 s-1]
          GTACC_M(I,M)=GTACC_M(I,M)/REAL(NDAY)
          QEVPACC_M(I,M)=QEVPACC_M(I,M)/REAL(NDAY)
          EVAPACC_M(I,M)=EVAPACC_M(I,M)   !became [kg m-2 day-1] instead of [kg m-2 s-1]
          HFSACC_M(I,M)=HFSACC_M(I,M)/REAL(NDAY)
          HMFNACC_M(I,M)=HMFNACC_M(I,M)/REAL(NDAY)
          ROFACC_M(I,M)=ROFACC_M(I,M)   !became [kg m-2 day-1] instead of [kg m-2 s-1
          OVRACC_M(I,M)=OVRACC_M(I,M)   !became [kg m-2 day-1] instead of [kg m-2 s-1]
          WTBLACC_M(I,M)=WTBLACC_M(I,M)/REAL(NDAY)
          DO 726 J=1,IGND
            TBARACC_M(I,M,J)=TBARACC_M(I,M,J)/REAL(NDAY)
            THLQACC_M(I,M,J)=THLQACC_M(I,M,J)/REAL(NDAY)
            THICACC_M(I,M,J)=THICACC_M(I,M,J)/REAL(NDAY)
            THALACC_M(I,M,J)=THALACC_M(I,M,J)/REAL(NDAY)
726       CONTINUE
C
          IF(FSINACC_M(I,M).GT.0.0) THEN
            ALVSACC_M(I,M)=ALVSACC_M(I,M)/(FSINACC_M(I,M)*0.5)
            ALIRACC_M(I,M)=ALIRACC_M(I,M)/(FSINACC_M(I,M)*0.5)
          ELSE
            ALVSACC_M(I,M)=0.0
            ALIRACC_M(I,M)=0.0
          ENDIF
C
          SNOACC_M(I,M)=SNOACC_M(I,M)/REAL(NDAY)
          if (SNOARE_M(I,M) .GT. 0.) THEN
             RHOSACC_M(I,M)=RHOSACC_M(I,M)/SNOARE_M(I,M)
             TSNOACC_M(I,M)=TSNOACC_M(I,M)/SNOARE_M(I,M)
             WSNOACC_M(I,M)=WSNOACC_M(I,M)/SNOARE_M(I,M)
          END IF
          TCANACC_M(I,M)=TCANACC_M(I,M)/REAL(NDAY)
          RCANACC_M(I,M)=RCANACC_M(I,M)/REAL(NDAY)
          SCANACC_M(I,M)=SCANACC_M(I,M)/REAL(NDAY)
          GROACC_M(I,M)=GROACC_M(I,M)/REAL(NDAY)
          FSINACC_M(I,M)=FSINACC_M(I,M)/REAL(NDAY)
          FLINACC_M(I,M)=FLINACC_M(I,M)/REAL(NDAY)
          FLUTACC_M(I,M)=FLUTACC_M(I,M)/REAL(NDAY)
          TAACC_M(I,M)=TAACC_M(I,M)/REAL(NDAY)
          UVACC_M(I,M)=UVACC_M(I,M)/REAL(NDAY)
          PRESACC_M(I,M)=PRESACC_M(I,M)/REAL(NDAY)
          QAACC_M(I,M)=QAACC_M(I,M)/REAL(NDAY)
          ALTOT=(ALVSACC_M(I,M)+ALIRACC_M(I,M))/2.0
          FSSTAR=FSINACC_M(I,M)*(1.-ALTOT)
          FLSTAR=FLINACC_M(I,M)-FLUTACC_M(I,M)
          QH=HFSACC_M(I,M)
          QE=QEVPACC_M(I,M)
          QEVPACC_M_SAVE(I,M)=QEVPACC_M(I,M)   !FLAG!! What is the point of this? JM Apr 1 2015
          BEG=FSSTAR+FLSTAR-QH-QE
          SNOMLT=HMFNACC_M(I,M)
C
          IF(RHOSACC_M(I,M).GT.0.0) THEN
              ZSN=SNOACC_M(I,M)/RHOSACC_M(I,M)
          ELSE
              ZSN=0.0
          ENDIF
C
          IF(TCANACC_M(I,M).GT.0.01) THEN
              TCN=TCANACC_M(I,M)-TFREZ
          ELSE
              TCN=0.0
          ENDIF
C
          IF(TSNOACC_M(I,M).GT.0.01) THEN
              TSN=TSNOACC_M(I,M)-TFREZ
          ELSE
              TSN=0.0
          ENDIF
C
          GTOUT=GTACC_M(I,M)-TFREZ
C 
          if ((iyear .ge. jdsty) .and. (iyear .le. jdendy)) then
           if ((iday .ge. jdstd) .and. (iday .le. jdendd)) then
C
C         WRITE TO OUTPUT FILES
C
          WRITE(611,6100) IDAY,IYEAR,FSSTAR,FLSTAR,QH,QE,SNOMLT,
     1                    BEG,GTOUT,SNOACC_M(I,M),RHOSACC_M(I,M),
     2                    WSNOACC_M(I,M),ALTOT,ROFACC_M(I,M),
     3                    CUMSNO,' TILE ',M
            IF(IGND.GT.3) THEN
               WRITE(621,6201) IDAY,IYEAR,(TBARACC_M(I,M,J)-TFREZ,
     1                  THLQACC_M(I,M,J),THICACC_M(I,M,J),J=1,5),
     2                  ' TILE ',M
               WRITE(631,6201) IDAY,IYEAR,(TBARACC_M(I,M,J)-TFREZ,
     1                  THLQACC_M(I,M,J),THICACC_M(I,M,J),J=6,10),
     2                  ' TILE ',M
            ELSE
               WRITE(621,6200) IDAY,IYEAR,(TBARACC_M(I,M,J)-TFREZ,
     1                  THLQACC_M(I,M,J),THICACC_M(I,M,J),J=1,3),
     2                  TCN,RCANACC_M(I,M),SCANACC_M(I,M),TSN,ZSN,
     3                  ' TILE ',M
               WRITE(631,6300) IDAY,IYEAR,FSINACC_M(I,M),FLINACC_M(I,M),
     1                  TAACC_M(I,M)-TFREZ,UVACC_M(I,M),PRESACC_M(I,M),
     2                  QAACC_M(I,M),PREACC_M(I,M),EVAPACC_M(I,M),
     3                  ' TILE ',M
            ENDIF
C
           endif
          ENDIF ! IF write daily
C
C          INITIALIZTION FOR MOSAIC TILE AND GRID VARIABLES
C
            call resetclassaccum(nltest,nmtest)
!            PREACC_M(I,M)=0.
!            GTACC_M(I,M)=0.
!            QEVPACC_M(I,M)=0.
!            HFSACC_M(I,M)=0.
!            HMFNACC_M(I,M)=0.
!            ROFACC_M(I,M)=0.
!            SNOACC_M(I,M)=0.
!            OVRACC_M(I,M)=0.
!            WTBLACC_M(I,M)=0.
!            ALVSACC_M(I,M)=0.
!            ALIRACC_M(I,M)=0.
!            RHOSACC_M(I,M)=0.
!            TSNOACC_M(I,M)=0.
!            WSNOACC_M(I,M)=0.
!            TCANACC_M(I,M)=0.
!            RCANACC_M(I,M)=0.
!            SCANACC_M(I,M)=0.
!            GROACC_M(I,M)=0.
!            FSINACC_M(I,M)=0.
!            FLINACC_M(I,M)=0.
!            TAACC_M(I,M)=0.
!            UVACC_M(I,M)=0.
!            PRESACC_M(I,M)=0.
!            QAACC_M(I,M)=0.
!            EVAPACC_M(I,M)=0.
!            FLUTACC_M(I,M)=0.
! C
!            DO 759 J=1,IGND
!              TBARACC_M(I,M,J)=0.
!              THLQACC_M(I,M,J)=0.
!              THICACC_M(I,M,J)=0.
!              THALACC_M(I,M,J)=0.
! 759        CONTINUE
C
809   CONTINUE
808   CONTINUE

      ENDIF ! IF(NCOUNT.EQ.NDAY)
C
      ENDIF !  IF(.NOT.PARALLELRUN)
C
C=======================================================================
C
!       CTEM--------------\
!     Only bother with monthly calculations if we desire those outputs to be written out.
      if (iyear .ge. jmosty) then

        call class_monthly_aw(IDAY,IYEAR,NCOUNT,NDAY,SBC,DELT,
     1                       nltest,nmtest,ALVSROT,FAREROT,FSVHROW,
     2                       ALIRROT,FSIHROW,GTROT,FSSROW,FDLROW,
     3                       HFSROT,ROFROT,PREROW,QFSROT,QEVPROT,
     4                       SNOROT,TAROW,WSNOROT,TBARROT,THLQROT,
     5                       THICROT,TFREZ)
!       CTEM--------------/

! C     ACCUMULATE OUTPUT DATA FOR MONTHLY AVERAGED FIELDS FOR CLASS GRID-MEAN.
! C     FOR BOTH PARALLEL MODE AND STAND ALONE MODE
! C
!       FSSTAR_MO   =0.0
!       FLSTAR_MO   =0.0
!       QH_MO       =0.0
!       QE_MO       =0.0
!       ALTOT_MO    =0.0
! C
!       DO 820 I=1,NLTEST
!        DO 821 M=1,NMTEST
!           ALVSACC_MO(I)=ALVSACC_MO(I)+ALVSROT(I,M)*FAREROT(I,M)
!      1                  *FSVHROW(I)
!           ALIRACC_MO(I)=ALIRACC_MO(I)+ALIRROT(I,M)*FAREROT(I,M)
!      1                  *FSIHROW(I)
!           FLUTACC_MO(I)=FLUTACC_MO(I)+SBC*GTROT(I,M)**4*FAREROT(I,M)
!           FSINACC_MO(I)=FSINACC_MO(I)+FSSROW(I)*FAREROT(I,M)
!           FLINACC_MO(I)=FLINACC_MO(I)+FDLROW(I)*FAREROT(I,M)
!           HFSACC_MO(I) =HFSACC_MO(I)+HFSROT(I,M)*FAREROT(I,M)
!           QEVPACC_MO(I)=QEVPACC_MO(I)+QEVPROT(I,M)*FAREROT(I,M)
!           SNOACC_MO(I) =SNOACC_MO(I)+SNOROT(I,M)*FAREROT(I,M)
!           TAACC_MO(I)=TAACC_MO(I)+TAROW(I)*FAREROT(I,M)
! C
!           IF(SNOROT(I,M).GT.0.0) THEN
!            WSNOACC_MO(I)=WSNOACC_MO(I)+WSNOROT(I,M)*FAREROT(I,M)
!           ENDIF
! C
!           ROFACC_MO(I) =ROFACC_MO(I)+ROFROT(I,M)*FAREROT(I,M)*DELT
!           PREACC_MO(I) =PREACC_MO(I)+PREROW(I)*FAREROT(I,M)*DELT
!           EVAPACC_MO(I)=EVAPACC_MO(I)+QFSROT(I,M)*FAREROT(I,M)*DELT
! C
!           DO 823 J=1,IGND
!            TBARACC_MO(I,J)=TBARACC_MO(I,J)+TBARROT(I,M,J)*FAREROT(I,M)
!            THLQACC_MO(I,J)=THLQACC_MO(I,J)+THLQROT(I,M,J)*FAREROT(I,M)
!            THICACC_MO(I,J)=THICACC_MO(I,J)+THICROT(I,M,J)*FAREROT(I,M)
! 823       CONTINUE
! C
! 821    CONTINUE
! 820   CONTINUE
! C
       DO NT=1,NMON
        IF(IDAY.EQ.monthend(NT+1).AND.NCOUNT.EQ.NDAY)THEN
         IMONTH=NT
!         NDMONTH=(monthend(NT+1)-monthend(NT))*NDAY
! C
!         DO 824 I=1,NLTEST
!          IF(FSINACC_MO(I).GT.0.0) THEN
!           ALVSACC_MO(I)=ALVSACC_MO(I)/(FSINACC_MO(I)*0.5)
!           ALIRACC_MO(I)=ALIRACC_MO(I)/(FSINACC_MO(I)*0.5)
!          ELSE
!           ALVSACC_MO(I)=0.0
!           ALIRACC_MO(I)=0.0
!          ENDIF
!          FLUTACC_MO(I)=FLUTACC_MO(I)/REAL(NDMONTH)
!          FSINACC_MO(I)=FSINACC_MO(I)/REAL(NDMONTH)
!          FLINACC_MO(I)=FLINACC_MO(I)/REAL(NDMONTH)
!          HFSACC_MO(I) =HFSACC_MO(I)/REAL(NDMONTH)
!          QEVPACC_MO(I)=QEVPACC_MO(I)/REAL(NDMONTH)
!          SNOACC_MO(I) =SNOACC_MO(I)/REAL(NDMONTH)
!          WSNOACC_MO(I)=WSNOACC_MO(I)/REAL(NDMONTH)
!          ROFACC_MO(I) =ROFACC_MO(I)
!          PREACC_MO(I) =PREACC_MO(I)
!          EVAPACC_MO(I)=EVAPACC_MO(I)
!          TAACC_MO(I)=TAACC_MO(I)/REAL(NDMONTH)
!          DO J=1,IGND
!           TBARACC_MO(I,J)=TBARACC_MO(I,J)/REAL(NDMONTH)
!           THLQACC_MO(I,J)=THLQACC_MO(I,J)/REAL(NDMONTH)
!           THICACC_MO(I,J)=THICACC_MO(I,J)/REAL(NDMONTH)
!          ENDDO
! C
!          ALTOT_MO=(ALVSACC_MO(I)+ALIRACC_MO(I))/2.0
!          FSSTAR_MO=FSINACC_MO(I)*(1.-ALTOT_MO)
!          FLSTAR_MO=FLINACC_MO(I)-FLUTACC_MO(I)
!          QH_MO=HFSACC_MO(I)
!          QE_MO=QEVPACC_MO(I)
! C
!           WRITE(81,8100)IMONTH,IYEAR,FSSTAR_MO,FLSTAR_MO,QH_MO,
!      1                 QE_MO,SNOACC_MO(I),WSNOACC_MO(I),
!      2                 ROFACC_MO(I),PREACC_MO(I),EVAPACC_MO(I),
!      3                 TAACC_MO(I)-TFREZ
!           IF (IGND.GT.3) THEN
!            WRITE(82,8101)IMONTH,IYEAR,(TBARACC_MO(I,J)-TFREZ,
!      1                  THLQACC_MO(I,J),THICACC_MO(I,J),J=1,5)
!            WRITE(82,8101)IMONTH,IYEAR,(TBARACC_MO(I,J)-TFREZ,
!      1                  THLQACC_MO(I,J),THICACC_MO(I,J),J=6,10)
!           ELSE
!            WRITE(82,8102)IMONTH,IYEAR,(TBARACC_MO(I,J)-TFREZ,
!      1                  THLQACC_MO(I,J),THICACC_MO(I,J),J=1,3)
!           ENDIF
C
C ADD INITIALIZTION FOR MONTHLY ACCUMULATED ARRAYS
C
!      call resetclassmon(nltest)
!          ALVSACC_MO(I)=0.0
!          ALIRACC_MO(I)=0.0
!          FLUTACC_MO(I)=0.0
!          FSINACC_MO(I)=0.0
!          FLINACC_MO(I)=0.0
!          HFSACC_MO(I)=0.0
!          QEVPACC_MO(I)=0.0
!          SNOACC_MO(I)=0.0
!          WSNOACC_MO(I)=0.0
!          ROFACC_MO(I)=0.0
!          PREACC_MO(I)=0.0
!          EVAPACC_MO(I)=0.0
!          TAACC_MO(I)=0.0
!          DO 826 J=1,IGND
!            TBARACC_MO(I,J)=0.
!            THLQACC_MO(I,J)=0.
!            THICACC_MO(I,J)=0.
! 826      CONTINUE
C
! 824     CONTINUE ! I
! C
        ENDIF ! IF(IDAY.EQ.monthend(NT+1).AND.NCOUNT.EQ.NDAY)
       ENDDO ! NMON

!       CTEM--------------\

      end if !skip the monthly calculations/writing unless iyear>=jmosty
!       CTEM--------------/

C
! 8100  FORMAT(1X,I4,I5,5(F8.2,1X),F8.3,F12.4,3(E12.3,1X),2(A6,I2))
! 8101  FORMAT(1X,I4,I5,5(F7.2,1X,2F6.3,1X),2(A6,I2))
! 8102  FORMAT(1X,I4,I5,3(F8.2,1X,2F6.3,1X),2(A6,I2))


C
C     ACCUMULATE OUTPUT DATA FOR YEARLY AVERAGED FIELDS FOR CLASS GRID-MEAN.
C     FOR BOTH PARALLEL MODE AND STAND ALONE MODE
C
      FSSTAR_YR   =0.0
      FLSTAR_YR   =0.0
      QH_YR       =0.0
      QE_YR       =0.0
      ALTOT_YR    =0.0
C
      DO 827 I=1,NLTEST
       DO 828 M=1,NMTEST
          ALVSACC_YR(I)=ALVSACC_YR(I)+ALVSROT(I,M)*FAREROT(I,M)
     1                  *FSVHROW(I)
          ALIRACC_YR(I)=ALIRACC_YR(I)+ALIRROT(I,M)*FAREROT(I,M)
     1                  *FSIHROW(I)
          FLUTACC_YR(I)=FLUTACC_YR(I)+SBC*GTROT(I,M)**4*FAREROT(I,M)
          FSINACC_YR(I)=FSINACC_YR(I)+FSSROW(I)*FAREROT(I,M)
          FLINACC_YR(I)=FLINACC_YR(I)+FDLROW(I)*FAREROT(I,M)
          HFSACC_YR(I) =HFSACC_YR(I)+HFSROT(I,M)*FAREROT(I,M)
          QEVPACC_YR(I)=QEVPACC_YR(I)+QEVPROT(I,M)*FAREROT(I,M)
          TAACC_YR(I)=TAACC_YR(I)+TAROW(I)*FAREROT(I,M)
          ROFACC_YR(I) =ROFACC_YR(I)+ROFROT(I,M)*FAREROT(I,M)*DELT
          PREACC_YR(I) =PREACC_YR(I)+PREROW(I)*FAREROT(I,M)*DELT
          EVAPACC_YR(I)=EVAPACC_YR(I)+QFSROT(I,M)*FAREROT(I,M)*DELT
828    CONTINUE
827   CONTINUE
C
      IF (IDAY.EQ.365.AND.NCOUNT.EQ.NDAY) THEN
C
       DO 829 I=1,NLTEST
         IF(FSINACC_YR(I).GT.0.0) THEN
          ALVSACC_YR(I)=ALVSACC_YR(I)/(FSINACC_YR(I)*0.5)
          ALIRACC_YR(I)=ALIRACC_YR(I)/(FSINACC_YR(I)*0.5)
         ELSE
          ALVSACC_YR(I)=0.0
          ALIRACC_YR(I)=0.0
         ENDIF
         FLUTACC_YR(I)=FLUTACC_YR(I)/(REAL(NDAY)*365.)
         FSINACC_YR(I)=FSINACC_YR(I)/(REAL(NDAY)*365.)
         FLINACC_YR(I)=FLINACC_YR(I)/(REAL(NDAY)*365.)
         HFSACC_YR(I) =HFSACC_YR(I)/(REAL(NDAY)*365.)
         QEVPACC_YR(I)=QEVPACC_YR(I)/(REAL(NDAY)*365.)
         ROFACC_YR(I) =ROFACC_YR(I)
         PREACC_YR(I) =PREACC_YR(I)
         EVAPACC_YR(I)=EVAPACC_YR(I)
         TAACC_YR(I)=TAACC_YR(I)/(REAL(NDAY)*365.)
C
         ALTOT_YR=(ALVSACC_YR(I)+ALIRACC_YR(I))/2.0
         FSSTAR_YR=FSINACC_YR(I)*(1.-ALTOT_YR)
         FLSTAR_YR=FLINACC_YR(I)-FLUTACC_YR(I)
         QH_YR=HFSACC_YR(I)
         QE_YR=QEVPACC_YR(I)
C
         WRITE(*,*) 'IYEAR=',IYEAR,' CLIMATE YEAR=',CLIMIYEAR

         WRITE(83,8103)IYEAR,FSSTAR_YR,FLSTAR_YR,QH_YR,
     1                  QE_YR,ROFACC_YR(I),PREACC_YR(I),
     2                  EVAPACC_YR(I)
C
C ADD INITIALIZTION FOR YEARLY ACCUMULATED ARRAYS
C
      call resetclassyr(nltest)
!          ALVSACC_YR(I)=0.0
!          ALIRACC_YR(I)=0.0
!          FLUTACC_YR(I)=0.0
!          FSINACC_YR(I)=0.0
!          FLINACC_YR(I)=0.0
!          HFSACC_YR(I)=0.0
!          QEVPACC_YR(I)=0.0
!          ROFACC_YR(I)=0.0
!          PREACC_YR(I)=0.0
!          EVAPACC_YR(I)=0.0
!          TAACC_YR(I)=0.0
C
829    CONTINUE ! I
C
      ENDIF ! IDAY.EQ.365 .AND. NDAY
C
8103  FORMAT(1X,I5,4(F8.2,1X),F12.4,1X,2(F12.3,1X),2(A5,I1))
C
c     CTEM output and write out
c
      if(.not.parallelrun) then ! stand alone mode, includes daily and yearly mosaic-mean output for ctem
c
c     calculate daily outputs from ctem
c
       if (ctem_on) then
        if(ncount.eq.nday) then
         call ctem_daily_aw(nltest,nmtest,iday,FAREROT,
     1                      iyear,jdstd,jdsty,jdendd,jdendy,grclarea)
! c
!         do i=1,nltest
!            do j=1,icc
!              ifcancmx_g(i,j)=0
!            enddo
!         enddo
! c
!         do i=1,nltest
!            do m=1,nmtest
!               ifcancmx_m(i,m)=0   !0=bare soil tile; 1=tile with vegetation
!               leaflitr_m(i,m)=0.0
!               tltrleaf_m(i,m)=0.0
!               tltrstem_m(i,m)=0.0
!               tltrroot_m(i,m)=0.0
!               ailcg_m(i,m)=0.0
!               ailcb_m(i,m)=0.0
!               afrleaf_m(i,m)=0.0
!               afrstem_m(i,m)=0.0
!               afrroot_m(i,m)=0.0
!               veghght_m(i,m)=0.0
!               rootdpth_m(i,m)=0.0
!               roottemp_m(i,m)=0.0
!               slai_m(i,m)=0.0
!               gleafmas_m(i,m) = 0.0
!               bleafmas_m(i,m) = 0.0
!               stemmass_m(i,m) = 0.0
!               rootmass_m(i,m) = 0.0
!               litrmass_m(i,m) = 0.0
!               soilcmas_m(i,m) = 0.0
!
! c
!               do j=1,icc
!                 if (fcancmxrow(i,m,j) .gt.0.0) then
!                 ifcancmx_g(i,j)=1
!                 ifcancmx_m(i,m)=1
!                 endif
!               enddo
! c
!               do k=1,ignd
!                 rmatctem_m(i,m,k)=0.0
!               enddo
! c
!            enddo ! m
!         enddo  ! i
! c
! c       ---------------------------------------------------------
! c
!         call resetctem_g(nltest)
!
!          do 851 i=1,nltest
! !           gpp_g(i) =0.0
! !           npp_g(i) =0.0
! !           nep_g(i) =0.0
! !           nbp_g(i) =0.0
! !           autores_g(i) =0.0
! !           hetrores_g(i)=0.0
! !           litres_g(i) =0.0
! !           socres_g(i) =0.0
! !           dstcemls_g(i)=0.0
! !           dstcemls3_g(i)=0.0
! !           litrfall_g(i)=0.0
! !           humiftrs_g(i)=0.0
! !           rml_g(i) =0.0
! !           rms_g(i) =0.0
! !           rmr_g(i) =0.0
! !           rg_g(i) =0.0
! !           vgbiomas_g(i) =0.0
! !           totcmass_g(i) =0.0
! !           gavglai_g(i) =0.0
! !           gavgltms_g(i) =0.0
! !           gavgscms_g(i) =0.0
! !           ailcg_g(i)=0.0
! !           ailcb_g(i)=0.0
! !           tcanoacc_out_g(i) =0.0
! !           burnfrac_g(i) =0.0
! !           probfire_g(i) =0.0
! !           lucemcom_g(i) =0.0
! !           lucltrin_g(i) =0.0
! !           lucsocin_g(i) =0.0
! !           emit_co2_g(i) =0.0
! !           emit_co_g(i)  =0.0
! !           emit_ch4_g(i) =0.0
! !           emit_nmhc_g(i) =0.0
! !           emit_h2_g(i) =0.0
! !           emit_nox_g(i) =0.0
! !           emit_n2o_g(i) =0.0
! !           emit_pm25_g(i) =0.0
! !           emit_tpm_g(i) =0.0
! !           emit_tc_g(i) =0.0
! !           emit_oc_g(i) =0.0
! !           emit_bc_g(i) =0.0
! !           bterm_g(i)   =0.0
! !           lterm_g(i)   =0.0
! !           mterm_g(i)   =0.0
! !           leaflitr_g(i)=0.0
! !           tltrleaf_g(i)=0.0
! !           tltrstem_g(i)=0.0
! !           tltrroot_g(i)=0.0
! !           gleafmas_g(i)=0.0
! !           bleafmas_g(i)=0.0
! !           stemmass_g(i)=0.0
! !           rootmass_g(i)=0.0
! !           litrmass_g(i)=0.0
! !           soilcmas_g(i)=0.0
! !           veghght_g(i)=0.0
! !           rootdpth_g(i)=0.0
! !           roottemp_g(i)=0.0
! !           slai_g(i)=0.0
! ! c                         !Rudra added CH4 realted variables on 03/12/2013
! !           CH4WET1_G(i) = 0.0
! !           CH4WET2_G(i) = 0.0
! !           WETFDYN_G(i) = 0.0
! !           CH4DYN1_G(i) = 0.0
! !           CH4DYN2_G(i) = 0.0
! !
! !           do k=1,ignd
! !            rmatctem_g(i,k)=0.0
! !           enddo
! ! c
! !           do j=1,icc
! !             afrleaf_g(i,j)=0.0
! !             afrstem_g(i,j)=0.0
! !             afrroot_g(i,j)=0.0
! !           enddo
! c
!           do 852 m=1,nmtest
! c
!            do j=1,icc
!              leaflitr_m(i,m)=leaflitr_m(i,m)+
!      &                       leaflitrrow(i,m,j)*fcancmxrow(i,m,j)
!              tltrleaf_m(i,m)=tltrleaf_m(i,m)+
!      &                       tltrleafrow(i,m,j)*fcancmxrow(i,m,j)
!              tltrstem_m(i,m)=tltrstem_m(i,m)+
!      &                       tltrstemrow(i,m,j)*fcancmxrow(i,m,j)
!              tltrroot_m(i,m)=tltrroot_m(i,m)+
!      &                       tltrrootrow(i,m,j)*fcancmxrow(i,m,j)
!              veghght_m(i,m)=veghght_m(i,m)+
!      &                            veghghtrow(i,m,j)*fcancmxrow(i,m,j)
!              rootdpth_m(i,m)=rootdpth_m(i,m)+
!      &                            rootdpthrow(i,m,j)*fcancmxrow(i,m,j)
!              roottemp_m(i,m)=roottemp_m(i,m)+
!      &                            roottemprow(i,m,j)*fcancmxrow(i,m,j)
!              slai_m(i,m)=slai_m(i,m)+slairow(i,m,j)*fcancmxrow(i,m,j)
! c
!              afrleaf_m(i,m)=afrleaf_m(i,m)+
!      &                              afrleafrow(i,m,j)*fcancmxrow(i,m,j)
!              afrstem_m(i,m)=afrstem_m(i,m)+
!      &                              afrstemrow(i,m,j)*fcancmxrow(i,m,j)
!              afrroot_m(i,m)=afrroot_m(i,m)+
!      &                              afrrootrow(i,m,j)*fcancmxrow(i,m,j)
! c
!              ailcg_m(i,m)=ailcg_m(i,m)+ailcgrow(i,m,j)*fcancmxrow(i,m,j)
!              ailcb_m(i,m)=ailcb_m(i,m)+ailcbrow(i,m,j)*fcancmxrow(i,m,j)
!
!              gleafmas_m(i,m) = gleafmas_m(i,m) + gleafmasrow(i,m,j)
!      &                                          *fcancmxrow(i,m,j)
!              bleafmas_m(i,m) = bleafmas_m(i,m) + bleafmasrow(i,m,j)
!      &                                          *fcancmxrow(i,m,j)
!              stemmass_m(i,m) = stemmass_m(i,m) + stemmassrow(i,m,j)
!      &                                          *fcancmxrow(i,m,j)
!              rootmass_m(i,m) = rootmass_m(i,m) + rootmassrow(i,m,j)
!      &                                          *fcancmxrow(i,m,j)
!              litrmass_m(i,m) = litrmass_m(i,m) + litrmassrow(i,m,j)
!      &                                          *fcancmxrow(i,m,j)
!              soilcmas_m(i,m) = soilcmas_m(i,m) + soilcmasrow(i,m,j)
!      &                                          *fcancmxrow(i,m,j)
!
! c
!              do k=1,ignd
!                 rmatctem_m(i,m,k)=rmatctem_m(i,m,k)+
!      &                            rmatctemrow(i,m,j,k)*fcancmxrow(i,m,j)
!              enddo
!            enddo
! c
!            npprow(i,m)     =npprow(i,m)*1.0377 ! convert to gc/m2.day
!            gpprow(i,m)     =gpprow(i,m)*1.0377 ! convert to gc/m2.day
!            neprow(i,m)     =neprow(i,m)*1.0377 ! convert to gc/m2.day
!            nbprow(i,m)     =nbprow(i,m)*1.0377 ! convert to gc/m2.day
!            lucemcomrow(i,m)=lucemcomrow(i,m)*1.0377 ! convert to gc/m2.day
!            lucltrinrow(i,m)=lucltrinrow(i,m)*1.0377 ! convert to gc/m2.day
!            lucsocinrow(i,m)=lucsocinrow(i,m)*1.0377 ! convert to gc/m2.day
! c
!            hetroresrow(i,m)=hetroresrow(i,m)*1.0377 ! convert to gc/m2.day
!            autoresrow(i,m) =autoresrow(i,m)*1.0377  ! convert to gc/m2.day
!            litresrow(i,m)  =litresrow(i,m)*1.0377   ! convert to gc/m2.day
!            socresrow(i,m)  =socresrow(i,m)*1.0377   ! convert to gc/m2.day
! c
!            CH4WET1ROW(i,m) = CH4WET1ROW(i,m)*1.0377 * 16.044 / 12. ! convert from umolCH4/m2/s to gCH4/m2.day
!            CH4WET2ROW(i,m) = CH4WET2ROW(i,m)*1.0377 * 16.044 / 12. ! convert from umolCH4/m2/s to gCH4/m2.day
!            CH4DYN1ROW(i,m) = CH4DYN1ROW(i,m)*1.0377 * 16.044 / 12. ! convert from umolCH4/m2/s to gCH4/m2.day
!            CH4DYN2ROW(i,m) = CH4DYN2ROW(i,m)*1.0377 * 16.044 / 12. ! convert from umolCH4/m2/s to gCH4/m2.day
! c
! c          write daily ctem results
! c
!            if ((iyd.ge.jdst).and.(iyd.le.jdend)) then
! c
! c             write grid-averaged fluxes of basic quantities to
! c             file *.CT01D_M
! c
!              if (mosaic) then
!               write(72,8200)iday,iyear,gpprow(i,m),npprow(i,m),
!      1                neprow(i,m),nbprow(i,m),autoresrow(i,m),
!      2                hetroresrow(i,m),litresrow(i,m),socresrow(i,m),
!      3                (dstcemlsrow(i,m)+dstcemls3row(i,m)),
!      4               litrfallrow(i,m),humiftrsrow(i,m),' TILE ',m,'AVGE'
!              end if
!
! c             write breakdown of some of basic fluxes to file *.CT3
! c             and selected litter fluxes for selected pft
!
! !              First for the bare fraction of the grid cell.
!              hetroresvegrow(i,m,iccp1)=hetroresvegrow(i,m,iccp1)*1.0377 ! convert to gc/m2.day
!              litresvegrow(i,m,iccp1)=litresvegrow(i,m,iccp1)*1.0377 ! convert to gc/m2.day
!              soilcresvegrow(i,m,iccp1)=soilcresvegrow(i,m,iccp1)*1.0377 ! convert to gc/m2.day
!
! c
!               do 853 j=1,icc
! c
!                 if (fcancmxrow(i,m,j) .gt.0.0) then
! c
!                  gppvegrow(i,m,j)=gppvegrow(i,m,j)*1.0377 ! convert to gc/m2.day
!                  nppvegrow(i,m,j)=nppvegrow(i,m,j)*1.0377 ! convert to gc/m2.day
!                  nepvegrow(i,m,j)=nepvegrow(i,m,j)*1.0377 ! convert to gc/m2.day
!                  nbpvegrow(i,m,j)=nbpvegrow(i,m,j)*1.0377 ! convert to gc/m2.day
!                  hetroresvegrow(i,m,j)=hetroresvegrow(i,m,j)*1.0377 ! convert to gc/m2.day
!                  autoresvegrow(i,m,j)=autoresvegrow(i,m,j)*1.0377 ! convert to gc/m2.day
!                  litresvegrow(i,m,j)=litresvegrow(i,m,j)*1.0377 ! convert to gc/m2.day
!                  soilcresvegrow(i,m,j)=soilcresvegrow(i,m,j)*1.0377 ! convert to gc/m2.day
!
!
! c                write to file .CT01D_M
!                  if (mosaic) then
!                   write(72,8201)iday,iyear,gppvegrow(i,m,j),
!      1            nppvegrow(i,m,j),nepvegrow(i,m,j),
!      2            ' TILE ',m,'PFT',j
!
!
! c                write to file .CT02D_M
!                  write(73,8300)iday,iyear,rmlvegaccrow(i,m,j),
!      1           rmsvegrow(i,m,j),rmrvegrow(i,m,j),rgvegrow(i,m,j),
!      2           leaflitrrow(i,m,j),tltrleafrow(i,m,j),
!      3          tltrstemrow(i,m,j),tltrrootrow(i,m,j),' TILE ',m,'PFT',j
! c
! c
! c                write grid-averaged pool sizes and component sizes for
! c                seleced ctem pft to file *.CT03D_M
! c
!                  write(74,8401)iday,iyear,vgbiomas_vegrow(i,m,j),
!      1               ailcgrow(i,m,j),gleafmasrow(i,m,j),
!      3               bleafmasrow(i,m,j), stemmassrow(i,m,j),
!      4               rootmassrow(i,m,j), litrmassrow(i,m,j),
!      5               soilcmasrow(i,m,j),' TILE ',m,'PFT',j
! c
! c                write lai, rmatctem, & structural attributes for selected
! c                pft to file *.CT04D_M
! c
!                  write(75,8500)iday,iyear, ailcgrow(i,m,j),
!      1                ailcbrow(i,m,j),(rmatctemrow(i,m,j,k),k=1,3),
!      2                veghghtrow(i,m,j),rootdpthrow(i,m,j),
!      3              roottemprow(i,m,j),slairow(i,m,j),' TILE ',m,'PFT',j
! c
! c                write allocation fractions for selected pft to
! c                file *.CT05D_M
! c
!                  write(76,8600)iday,iyear, afrleafrow(i,m,j),
!      1                afrstemrow(i,m,j),afrrootrow(i,m,j),
!      2                tcanoaccrow_out(i,m), lfstatusrow(i,m,j),
!      3                ' TILE ',m,'PFT',j
!
! c             write fire and luc results to file *.CT06D_M
! c
!               if (dofire .or. lnduseon) then   !FLAG FIX THIS
!                write(78,8800)iday,iyear,
!      1         emit_co2row(i,m,j),emit_corow(i,m,j),emit_ch4row(i,m,j),
!      2         emit_nmhcrow(i,m,j),emit_h2row(i,m,j),emit_noxrow(i,m,j),
!      3         emit_n2orow(i,m,j),emit_pm25row(i,m,j),
!      4         emit_tpmrow(i,m,j),emit_tcrow(i,m,j),emit_ocrow(i,m,j),
!      5         emit_bcrow(i,m,j),burnvegfrow(i,m,j),probfirerow(i,m),
!      6         lucemcomrow(i,m),lucltrinrow(i,m), lucsocinrow(i,m),
!      7         grclarea(i),btermrow(i,m),ltermrow(i,m),mtermrow(i,m),
!      8         ' TILE ',m,'PFT',j
!                endif
!
!               end if !mosaic
!
!               endif  !if (fcancmxrow(i,m,j) .gt.0.0) then
! c
! 853           continue
! c
!               if (mosaic) then
!                if (ifcancmx_m(i,m) .gt. 0) then
!
!
! c               write to file .CT02D_M
!                 write(73,8300)iday,iyear,rmlrow(i,m),rmsrow(i,m),
!      1          rmrrow(i,m),rgrow(i,m),leaflitr_m(i,m),tltrleaf_m(i,m),
!      2          tltrstem_m(i,m),tltrroot_m(i,m),' TILE ',m,'AVGE'
! c
! c               write to file .CT03D_M
!                 write(74,8402)iday,iyear,vgbiomasrow(i,m),
!      1               gavglairow(i,m),gavgltmsrow(i,m),
!      2               gavgscmsrow(i,m),gleafmasrow(i,m,j),
!      3               bleafmasrow(i,m,j), stemmassrow(i,m,j),
!      4               rootmassrow(i,m,j), litrmassrow(i,m,j),
!      5               soilcmasrow(i,m,j),' TILE ',m, 'AVGE'
! c
! c               write to file .CT04D_M
!                 write(75,8500)iday,iyear,ailcg_m(i,m),
!      1                ailcb_m(i,m),(rmatctem_m(i,m,k),k=1,3),
!      2                veghght_m(i,m),rootdpth_m(i,m),
!      3                roottemp_m(i,m),slai_m(i,m),' TILE ',m, 'AVGE'
! c
! c               write to file .CT05D_M
!                 write(76,8601)iday,iyear, afrleaf_m(i,m),
!      1                afrstem_m(i,m),afrroot_m(i,m),
!      2                tcanoaccrow_out(i,m),
!      3                ' TILE ',m,'AVGE'
!
!                end if !if (ifcancmx_m(i,m) .gt.0.0) then
!               endif !mosaic
! c
!            endif ! if ((iyd.ge.jdst).and.(iyd.le.jdend))
! c
! 8200       format(1x,i4,i5,11f10.5,2(a6,i2))
! 8201       format(1x,i4,i5,3f10.5,80x,2(a6,i2))
! 8300       format(1x,i4,i5,8f10.5,2(a6,i2))
! 8301       format(1x,i4,i5,4f10.5,40x,2(a6,i2))
! 8400       format(1x,i4,i5,11f10.5,2(a6,i2))
! 8401       format(1x,i4,i5,2f10.5,6f10.5,2(a6,i2))
! 8402       format(1x,i4,i5,10f10.5,2(a6,i2))
! !                   8402       format(1x,i4,i5,2f10.5,40x,2f10.5,2(a6,i2))
! 8500       format(1x,i4,i5,9f10.5,2(a6,i2))
! 8600       format(1x,i4,i5,4f10.5,i8,2(a6,i2))
! 8601       format(1x,i4,i5,4f10.5,8x,2(a6,i2))
! 8800       format(1x,i4,i5,20f11.4,2x,f9.2,2(a6,i2))
! 8810       format(1x,i4,i5,5f11.4,2(a6,i2))
! c
! c          Calculation of grid averaged variables
! c
!            gpp_g(i) =gpp_g(i) + gpprow(i,m)*FAREROT(i,m)
!            npp_g(i) =npp_g(i) + npprow(i,m)*FAREROT(i,m)
!            nep_g(i) =nep_g(i) + neprow(i,m)*FAREROT(i,m)
!            nbp_g(i) =nbp_g(i) + nbprow(i,m)*FAREROT(i,m)
!            autores_g(i) =autores_g(i) +autoresrow(i,m)*FAREROT(i,m)
!            hetrores_g(i)=hetrores_g(i)+hetroresrow(i,m)*FAREROT(i,m)
!            litres_g(i) =litres_g(i) + litresrow(i,m)*FAREROT(i,m)
!            socres_g(i) =socres_g(i) + socresrow(i,m)*FAREROT(i,m)
!            dstcemls_g(i)=dstcemls_g(i)+dstcemlsrow(i,m)*FAREROT(i,m)
!            dstcemls3_g(i)=dstcemls3_g(i)
!      &                      +dstcemls3row(i,m)*FAREROT(i,m)
!
!            litrfall_g(i)=litrfall_g(i)+litrfallrow(i,m)*FAREROT(i,m)
!            humiftrs_g(i)=humiftrs_g(i)+humiftrsrow(i,m)*FAREROT(i,m)
!            rml_g(i) =rml_g(i) + rmlrow(i,m)*FAREROT(i,m)
!            rms_g(i) =rms_g(i) + rmsrow(i,m)*FAREROT(i,m)
!            rmr_g(i) =rmr_g(i) + rmrrow(i,m)*FAREROT(i,m)
!            rg_g(i) =rg_g(i) + rgrow(i,m)*FAREROT(i,m)
!            leaflitr_g(i) = leaflitr_g(i) + leaflitr_m(i,m)*FAREROT(i,m)
!            tltrleaf_g(i) = tltrleaf_g(i) + tltrleaf_m(i,m)*FAREROT(i,m)
!            tltrstem_g(i) = tltrstem_g(i) + tltrstem_m(i,m)*FAREROT(i,m)
!            tltrroot_g(i) = tltrroot_g(i) + tltrroot_m(i,m)*FAREROT(i,m)
!            vgbiomas_g(i) =vgbiomas_g(i) + vgbiomasrow(i,m)*FAREROT(i,m)
!            gavglai_g(i) =gavglai_g(i) + gavglairow(i,m)*FAREROT(i,m)
!            gavgltms_g(i) =gavgltms_g(i) + gavgltmsrow(i,m)*FAREROT(i,m)
!            gavgscms_g(i) =gavgscms_g(i) + gavgscmsrow(i,m)*FAREROT(i,m)
!            tcanoacc_out_g(i) =tcanoacc_out_g(i)+
!      1                        tcanoaccrow_out(i,m)*FAREROT(i,m)
!            totcmass_g(i) =vgbiomas_g(i) + gavgltms_g(i) + gavgscms_g(i)
!            gleafmas_g(i) = gleafmas_g(i) + gleafmas_m(i,m)*FAREROT(i,m)
!            bleafmas_g(i) = bleafmas_g(i) + bleafmas_m(i,m)*FAREROT(i,m)
!            stemmass_g(i) = stemmass_g(i) + stemmass_m(i,m)*FAREROT(i,m)
!            rootmass_g(i) = rootmass_g(i) + rootmass_m(i,m)*FAREROT(i,m)
!            litrmass_g(i) = litrmass_g(i) + litrmass_m(i,m)*FAREROT(i,m)
!            soilcmas_g(i) = soilcmas_g(i) + soilcmas_m(i,m)*FAREROT(i,m)
! c
!            burnfrac_g(i) =burnfrac_g(i)+ burnfracrow(i,m)*FAREROT(i,m)
!
!            probfire_g(i) =probfire_g(i)+probfirerow(i,m)*FAREROT(i,m)
!            lucemcom_g(i) =lucemcom_g(i)+lucemcomrow(i,m)*FAREROT(i,m)
!            lucltrin_g(i) =lucltrin_g(i)+lucltrinrow(i,m)*FAREROT(i,m)
!            lucsocin_g(i) =lucsocin_g(i)+lucsocinrow(i,m)*FAREROT(i,m)
!            bterm_g(i)    =bterm_g(i)   +btermrow(i,m)*FAREROT(i,m)
!            lterm_g(i)    =lterm_g(i)   +ltermrow(i,m)*FAREROT(i,m)
!            mterm_g(i)    =mterm_g(i)   +mtermrow(i,m)*FAREROT(i,m)
! c                                                   !Rudra added CH4 related variables on 03/12/2013
!            CH4WET1_G(i) = CH4WET1_G(i) + CH4WET1ROW(i,m)*FAREROT(i,m)
!            CH4WET2_G(i) = CH4WET2_G(i) + CH4WET2ROW(i,m)*FAREROT(i,m)
!            WETFDYN_G(i) = WETFDYN_G(i) + WETFDYNROW(i,m)*FAREROT(i,m)
!            CH4DYN1_G(i) = CH4DYN1_G(i) + CH4DYN1ROW(i,m)*FAREROT(i,m)
!            CH4DYN2_G(i) = CH4DYN2_G(i) + CH4DYN2ROW(i,m)*FAREROT(i,m)
!
!            do j=1,icc
!
!             do k=1,ignd
!              rmatctem_g(i,k)=rmatctem_g(i,k)+rmatctemrow(i,m,j,k)
!      1                            *FAREROT(i,m)*fcancmxrow(i,m,j)
!             end do
!
!            veghght_g(i) = veghght_g(i) + veghghtrow(i,m,j)
!      1                            *FAREROT(i,m)*fcancmxrow(i,m,j)
!            rootdpth_g(i) = rootdpth_g(i) + rootdpthrow(i,m,j)
!      1                            *FAREROT(i,m)*fcancmxrow(i,m,j)
!            roottemp_g(i) = roottemp_g(i) + roottemprow(i,m,j)
!      1                            *FAREROT(i,m)*fcancmxrow(i,m,j)
!            slai_g(i) = slai_g(i) + slairow(i,m,j)
!      1                            *FAREROT(i,m)*fcancmxrow(i,m,j)
!
!            emit_co2_g(i) =emit_co2_g(i)+ emit_co2row(i,m,j)
!      1                            *FAREROT(i,m)*fcancmxrow(i,m,j)
!            emit_co_g(i)  =emit_co_g(i) + emit_corow(i,m,j)
!      1                            *FAREROT(i,m)*fcancmxrow(i,m,j)
!            emit_ch4_g(i) =emit_ch4_g(i)+ emit_ch4row(i,m,j)
!      1                            *FAREROT(i,m)*fcancmxrow(i,m,j)
!            emit_nmhc_g(i)=emit_nmhc_g(i)+emit_nmhcrow(i,m,j)
!      1                            *FAREROT(i,m)*fcancmxrow(i,m,j)
!            emit_h2_g(i)  =emit_h2_g(i) + emit_h2row(i,m,j)
!      1                            *FAREROT(i,m)*fcancmxrow(i,m,j)
!            emit_nox_g(i) =emit_nox_g(i)+ emit_noxrow(i,m,j)
!      1                            *FAREROT(i,m)*fcancmxrow(i,m,j)
!            emit_n2o_g(i) =emit_n2o_g(i)+ emit_n2orow(i,m,j)
!      1                            *FAREROT(i,m)*fcancmxrow(i,m,j)
!            emit_pm25_g(i)=emit_pm25_g(i)+emit_pm25row(i,m,j)
!      1                            *FAREROT(i,m)*fcancmxrow(i,m,j)
!            emit_tpm_g(i) =emit_tpm_g(i)+ emit_tpmrow(i,m,j)
!      1                            *FAREROT(i,m)*fcancmxrow(i,m,j)
!            emit_tc_g(i)  =emit_tc_g(i) + emit_tcrow(i,m,j)
!      1                            *FAREROT(i,m)*fcancmxrow(i,m,j)
!            emit_oc_g(i)  =emit_oc_g(i) + emit_ocrow(i,m,j)
!      1                            *FAREROT(i,m)*fcancmxrow(i,m,j)
!            emit_bc_g(i)  =emit_bc_g(i) + emit_bcrow(i,m,j)
!      1                            *FAREROT(i,m)*fcancmxrow(i,m,j)
!
!              ailcg_g(i)=ailcg_g(i)+ailcgrow(i,m,j)*fcancmxrow(i,m,j)
!      &                   *FAREROT(i,m)
!              ailcb_g(i)=ailcb_g(i)+ailcbrow(i,m,j)*fcancmxrow(i,m,j)
!      &                   *FAREROT(i,m)
!
!            enddo
! c
! 852       continue
! c
!           if ((iyd.ge.jdst).and.(iyd.le.jdend)) then
! c           write to file .CT01D_G
!             write(721,8200)iday,iyear,gpp_g(i),npp_g(i),
!      1                nep_g(i),nbp_g(i),autores_g(i),
!      2                hetrores_g(i),litres_g(i),socres_g(i),
!      3                (dstcemls_g(i)+dstcemls3_g(i)),
!      4                litrfall_g(i),humiftrs_g(i)
!
! c           write breakdown of some of basic fluxes to file
! c           *.CT02D_G and selected litter fluxes for selected pft
!             write(731,8300)iday,iyear,rml_g(i),rms_g(i),
!      1          rmr_g(i),rg_g(i),leaflitr_g(i),tltrleaf_g(i),
!      2          tltrstem_g(i),tltrroot_g(i)
! c
! c           write to file .CT03D_G
!             write(741,8400)iday,iyear,vgbiomas_g(i),
!      1               gavglai_g(i),gavgltms_g(i),
!      2               gavgscms_g(i), totcmass_g(i),
!      3               gleafmas_g(i), bleafmas_g(i), stemmass_g(i),
!      4               rootmass_g(i), litrmass_g(i), soilcmas_g(i)
! c
! c           write to file .CT04D_G
!             write(751,8500)iday,iyear, ailcg_g(i),
!      1                ailcb_g(i),(rmatctem_g(i,k),k=1,3),
!      2                veghght_g(i),rootdpth_g(i),roottemp_g(i),slai_g(i)
! c
! c           write fire and luc results to file *.CT06D_G
! c
!            if (dofire .or. lnduseon) then
!             write(781,8800)iday,iyear,
!      1          emit_co2_g(i), emit_co_g(i), emit_ch4_g(i),
!      2          emit_nmhc_g(i), emit_h2_g(i), emit_nox_g(i),
!      3          emit_n2o_g(i), emit_pm25_g(i), emit_tpm_g(i),
!      4          emit_tc_g(i), emit_oc_g(i), emit_bc_g(i),
!      5          burnfrac_g(i)*100., probfire_g(i),lucemcom_g(i),
!      6          lucltrin_g(i), lucsocin_g(i),
!      7          grclarea(i), bterm_g(i), lterm_g(i), mterm_g(i)
!            endif
! c
! c      write CH4 variables to file *.CT08D_G
! c
!            if (dowetlands .or. obswetf) then
!             write(762,8810)iday,iyear, ch4wet1_g(i),
!      1                 ch4wet2_g(i), wetfdyn_g(i),
!      2                 ch4dyn1_g(i), ch4dyn2_g(i)
!            endif
!
! c
!
!             if (compete .or. lnduseon) then
!               sumfare=0.0
!               if (mosaic) then
!                do m=1,nmos
!                  sumfare=sumfare+FAREROT(i,m)
!                enddo
!                write(761,8200)iday,iyear,(FAREROT(i,m)*100.,m=1,nmos),
!      1                      sumfare
!               else !composite
!                do j=1,icc  !m = 1
!                  sumfare=sumfare+fcancmxrow(i,1,j)
!                enddo
!                write(761,8200)iday,iyear,(fcancmxrow(i,1,j)*100.,
!      1                      j=1,icc),(1.0-sumfare)*100.,sumfare
!               endif !mosaic/composite
!             endif !compete/lnduseon
! c
!           endif !if ((iyd.ge.jdst).and.(iyd.le.jdend)) then
! c
! 851     continue
! c
       endif ! if(ncount.eq.nday)
       endif ! if(ctem_on)
! c
       endif ! if(not.parallelrun)
c
c=======================================================================
c     Calculate monthly & yearly output for ctem
c

c     First initialize some output variables
c     initialization is done just before use.

      if (ctem_on) then
       if(ncount.eq.nday) then
c
        !do 861 i=1,nltest

c
          do nt=1,nmon
           if (iday.eq.mmday(nt)) then

            call resetmidmonth(nltest)

!             stemmass_mo_g(i)=0.0
!             rootmass_mo_g(i)=0.0
!             litrmass_mo_g(i)=0.0
!             soilcmas_mo_g(i)=0.0
!             vgbiomas_mo_g(i)=0.0
!             totcmass_mo_g(i)=0.0
           endif
          enddo
c

          if(iday.eq.monthend(imonth+1))then

            call resetmonthend_g(nltest)

!            laimaxg_mo_g(i)=0.0
!            npp_mo_g(i)=0.0
!            gpp_mo_g(i)=0.0
!            nep_mo_g(i)=0.0
!            nbp_mo_g(i)=0.0
!            hetrores_mo_g(i)=0.0
!            autores_mo_g(i)=0.0
!            litres_mo_g(i)=0.0
!            soilcres_mo_g(i)=0.0
!
!            emit_co2_mo_g(i)=0.0
!            emit_co_mo_g(i) =0.0
!            emit_ch4_mo_g(i) =0.0
!            emit_nmhc_mo_g(i) =0.0
!            emit_h2_mo_g(i) =0.0
!            emit_nox_mo_g(i) =0.0
!            emit_n2o_mo_g(i) =0.0
!            emit_pm25_mo_g(i) =0.0
!            emit_tpm_mo_g(i) =0.0
!            emit_tc_mo_g(i) =0.0
!            emit_oc_mo_g(i) =0.0
!            emit_bc_mo_g(i) =0.0
!            probfire_mo_g(i) =0.0
!            luc_emc_mo_g(i) =0.0
!            lucsocin_mo_g(i) =0.0
!            lucltrin_mo_g(i) =0.0
!            burnfrac_mo_g(i) =0.0
!            bterm_mo_g(i)    =0.0
!            lterm_mo_g(i)    =0.0
!            mterm_mo_g(i)    =0.0
! c          CH4(wetland) related variables !Rudra 04/12/2013
!            ch4wet1_mo_g(i)  =0.0
!            ch4wet2_mo_g(i)  =0.0
!            wetfdyn_mo_g(i)  =0.0
!            ch4dyn1_mo_g(i)  =0.0
!            ch4dyn2_mo_g(i)  =0.0

          endif !mid-month
c
          if (iday .eq. 365) then

            call resetyearend_g(nltest)

!            laimaxg_yr_g(i)=0.0
!            stemmass_yr_g(i)=0.0
!            rootmass_yr_g(i)=0.0
!            litrmass_yr_g(i)=0.0
!            soilcmas_yr_g(i)=0.0
!            vgbiomas_yr_g(i)=0.0
!            totcmass_yr_g(i)=0.0
!            npp_yr_g(i)=0.0
!            gpp_yr_g(i)=0.0
!            nep_yr_g(i)=0.0
!            nbp_yr_g(i)=0.0
!            hetrores_yr_g(i)=0.0
!            autores_yr_g(i)=0.0
!            litres_yr_g(i)=0.0
!            soilcres_yr_g(i)=0.0
!            emit_co2_yr_g(i)=0.0
!            emit_co_yr_g(i)=0.0
!            emit_ch4_yr_g(i)=0.0
!            emit_nmhc_yr_g(i)=0.0
!            emit_h2_yr_g(i)=0.0
!            emit_nox_yr_g(i)=0.0
!            emit_n2o_yr_g(i)=0.0
!            emit_pm25_yr_g(i)=0.0
!            emit_tpm_yr_g(i)=0.0
!            emit_tc_yr_g(i)=0.0
!            emit_oc_yr_g(i)=0.0
!            emit_bc_yr_g(i)=0.0
!            probfire_yr_g(i)=0.0
!            luc_emc_yr_g(i)=0.0
!            lucsocin_yr_g(i)=0.0
!            lucltrin_yr_g(i)=0.0
!            burnfrac_yr_g(i)=0.0
!            bterm_yr_g(i)=0.0
!            lterm_yr_g(i)=0.0
!            mterm_yr_g(i)=0.0
! c          CH4(wetland) related variables !Rudra 04/12/2013
!            ch4wet1_yr_g(i)  =0.0
!            ch4wet2_yr_g(i)  =0.0
!            wetfdyn_yr_g(i)  =0.0
!            ch4dyn1_yr_g(i)  =0.0
!            ch4dyn2_yr_g(i)  =0.0

          endif

861     continue
c

!       CTEM--------------\
!     Only bother with monthly calculations if we desire those outputs to be written out.
      if (iyear .ge. jmosty) then
!       CTEM--------------/

        call ctem_monthly_aw(nltest,nmtest,iday,FAREROT,iyear,nday)

! c       accumulate monthly outputs
! c
!         do 862 i=1,nltest
!
!          do 863 m=1,nmtest
!
!           do j=1,icc
!
!            if (ailcgrow(i,m,j) .gt. laimaxg_mo_m(i,m,j)) then
!             laimaxg_mo_m(i,m,j)=ailcgrow(i,m,j)
!            end if
!
!            npp_mo_m(i,m,j)=npp_mo_m(i,m,j)+nppvegrow(i,m,j)
!            gpp_mo_m(i,m,j)=gpp_mo_m(i,m,j)+gppvegrow(i,m,j)
!            nep_mo_m(i,m,j)=nep_mo_m(i,m,j)+nepvegrow(i,m,j)
!            nbp_mo_m(i,m,j)=nbp_mo_m(i,m,j)+nbpvegrow(i,m,j)
!            hetrores_mo_m(i,m,j)=hetrores_mo_m(i,m,j)
!      1                               +hetroresvegrow(i,m,j)
!            autores_mo_m(i,m,j) =autores_mo_m(i,m,j)+autoresvegrow(i,m,j)
!            litres_mo_m(i,m,j)  =litres_mo_m(i,m,j) +litresvegrow(i,m,j)
!            soilcres_mo_m(i,m,j) =soilcres_mo_m(i,m,j)
!      1                               +soilcresvegrow(i,m,j)
!            emit_co2_mo_m(i,m,j)=emit_co2_mo_m(i,m,j)+emit_co2row(i,m,j)
!            emit_co_mo_m(i,m,j) =emit_co_mo_m(i,m,j)+emit_corow(i,m,j)
!            emit_ch4_mo_m(i,m,j) =emit_ch4_mo_m(i,m,j)+emit_ch4row(i,m,j)
!            emit_nmhc_mo_m(i,m,j)=emit_nmhc_mo_m(i,m,j)
!      1                                     +emit_nmhcrow(i,m,j)
!            emit_h2_mo_m(i,m,j) =emit_h2_mo_m(i,m,j)+emit_h2row(i,m,j)
!            emit_nox_mo_m(i,m,j) =emit_nox_mo_m(i,m,j)+emit_noxrow(i,m,j)
!            emit_n2o_mo_m(i,m,j) =emit_n2o_mo_m(i,m,j)+emit_n2orow(i,m,j)
!            emit_pm25_mo_m(i,m,j)=emit_pm25_mo_m(i,m,j)
!      1                                     +emit_pm25row(i,m,j)
!            emit_tpm_mo_m(i,m,j) =emit_tpm_mo_m(i,m,j)+emit_tpmrow(i,m,j)
!            emit_tc_mo_m(i,m,j) =emit_tc_mo_m(i,m,j)+emit_tcrow(i,m,j)
!            emit_oc_mo_m(i,m,j) =emit_oc_mo_m(i,m,j)+emit_ocrow(i,m,j)
!            emit_bc_mo_m(i,m,j) =emit_bc_mo_m(i,m,j)+emit_bcrow(i,m,j)
!            burnfrac_mo_m(i,m,j) =burnfrac_mo_m(i,m,j)+burnvegfrow(i,m,j)
!
!           end do
!
!            nep_mo_m(i,m,iccp1)=nep_mo_m(i,m,iccp1)+nepvegrow(i,m,iccp1)
!            nbp_mo_m(i,m,iccp1)=nbp_mo_m(i,m,iccp1)+nbpvegrow(i,m,iccp1)
!            hetrores_mo_m(i,m,iccp1)=hetrores_mo_m(i,m,iccp1)
!      1                               +hetroresvegrow(i,m,iccp1)
!            litres_mo_m(i,m,iccp1)  =litres_mo_m(i,m,iccp1)
!      1                               +litresvegrow(i,m,iccp1)
!            soilcres_mo_m(i,m,iccp1) =soilcres_mo_m(i,m,iccp1)
!      1                              +soilcresvegrow(i,m,iccp1)
!
!            luc_emc_mo_m(i,m) =luc_emc_mo_m(i,m)
!      &                             +lucemcomrow(i,m)
!            lucsocin_mo_m(i,m) =lucsocin_mo_m(i,m)
!      &                             +lucsocinrow(i,m)
!            lucltrin_mo_m(i,m) =lucltrin_mo_m(i,m)
!      &                             +lucltrinrow(i,m)
! C                         !CH4 related variables !Rudra
!            ch4wet1_mo_m(i,m) = ch4wet1_mo_m(i,m) + CH4WET1ROW(i,m)
!            ch4wet2_mo_m(i,m) = ch4wet2_mo_m(i,m) + CH4WET2ROW(i,m)
!            wetfdyn_mo_m(i,m) = wetfdyn_mo_m(i,m) + WETFDYNROW(i,m)
!            ch4dyn1_mo_m(i,m) = ch4dyn1_mo_m(i,m) + CH4DYN1ROW(i,m)
!            ch4dyn2_mo_m(i,m) = ch4dyn2_mo_m(i,m) + CH4DYN2ROW(i,m)
!
! !          Sum the probfire now, later we will make it a per day value.
!            probfire_mo_m(i,m) =probfire_mo_m(i,m) + probfirerow(i,m)
!            bterm_mo_m(i,m) = bterm_mo_m(i,m) + btermrow(i,m)
!            lterm_mo_m(i,m) = lterm_mo_m(i,m) + ltermrow(i,m)
!            mterm_mo_m(i,m) = mterm_mo_m(i,m) + mtermrow(i,m)
! c
!            do 865 nt=1,nmon
! c
!              if(iday.eq.mmday(nt))then
!
!                do j=1,icc
!                 vgbiomas_mo_m(i,m,j)=0.0
!                 litrmass_mo_m(i,m,j)=0.0
!                 soilcmas_mo_m(i,m,j)=0.0
!                 totcmass_mo_m(i,m,j)=0.0
!                 stemmass_mo_m(i,m,j)=0.0
!                 rootmass_mo_m(i,m,j)=0.0
!                end do
!                 litrmass_mo_m(i,m,iccp1)=0.0
!                 soilcmas_mo_m(i,m,iccp1)=0.0
!
!                 do 867 j=1,icc
!
!                   vgbiomas_mo_m(i,m,j)=vgbiomas_vegrow(i,m,j)
!                   litrmass_mo_m(i,m,j)=litrmassrow(i,m,j)
!                   soilcmas_mo_m(i,m,j)=soilcmasrow(i,m,j)
!                   stemmass_mo_m(i,m,j)=stemmassrow(i,m,j)
!                   rootmass_mo_m(i,m,j)=rootmassrow(i,m,j)
!                   totcmass_mo_m(i,m,j)=vgbiomas_vegrow(i,m,j) +
!      1                        litrmassrow(i,m,j)+soilcmasrow(i,m,j)
!
! 867             continue
!
!                 ! Do the bare fraction too
!                 litrmass_mo_m(i,m,iccp1)=litrmassrow(i,m,iccp1)
!                 soilcmas_mo_m(i,m,iccp1)=soilcmasrow(i,m,iccp1)
!
!                 barefrac=1.0
!
!                do j=1,icc
!                 vgbiomas_mo_g(i)=vgbiomas_mo_g(i)+vgbiomas_mo_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 litrmass_mo_g(i)=litrmass_mo_g(i)+litrmass_mo_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 soilcmas_mo_g(i)=soilcmas_mo_g(i)+soilcmas_mo_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 stemmass_mo_g(i)=stemmass_mo_g(i)+stemmass_mo_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 rootmass_mo_g(i)=rootmass_mo_g(i)+rootmass_mo_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 totcmass_mo_g(i)=totcmass_mo_g(i)+totcmass_mo_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 barefrac=barefrac-FAREROT(i,m)*fcancmxrow(i,m,j)
!
!                end do
!
! !             Also add in the bare fraction contributions.
!               litrmass_mo_g(i)=litrmass_mo_g(i)+litrmass_mo_m(i,m,iccp1)
!      &                          *barefrac
!               soilcmas_mo_g(i)=soilcmas_mo_g(i)+soilcmas_mo_m(i,m,iccp1)
!      &                          *barefrac
!
!              endif ! mmday (mid-month instantaneous value)
! c
!              if(iday.eq.monthend(nt+1))then
!
!                ndmonth=(monthend(nt+1)-monthend(nt))*nday
!
!                barefrac=1.0
! c
!                do j=1,icc
!
!                 npp_mo_g(i)=npp_mo_g(i)+npp_mo_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 gpp_mo_g(i)=gpp_mo_g(i)+gpp_mo_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 nep_mo_g(i)=nep_mo_g(i)+nep_mo_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 nbp_mo_g(i)=nbp_mo_g(i)+nbp_mo_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 hetrores_mo_g(i)=hetrores_mo_g(i)+hetrores_mo_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 autores_mo_g(i) =autores_mo_g(i) +autores_mo_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 litres_mo_g(i)  =litres_mo_g(i) +litres_mo_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 soilcres_mo_g(i) =soilcres_mo_g(i)+ soilcres_mo_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 laimaxg_mo_g(i)=laimaxg_mo_g(i)+laimaxg_mo_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 emit_co2_mo_g(i)=emit_co2_mo_g(i)+emit_co2_mo_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 emit_co_mo_g(i) =emit_co_mo_g(i)+emit_co_mo_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 emit_ch4_mo_g(i) =emit_ch4_mo_g(i)+emit_ch4_mo_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 emit_nmhc_mo_g(i)=emit_nmhc_mo_g(i)+
!      &                            emit_nmhc_mo_m(i,m,j)*FAREROT(i,m)
!      &                                            *fcancmxrow(i,m,j)
!                 emit_h2_mo_g(i) =emit_h2_mo_g(i)+emit_h2_mo_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 emit_nox_mo_g(i) =emit_nox_mo_g(i)+emit_nox_mo_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 emit_n2o_mo_g(i) =emit_n2o_mo_g(i)+emit_n2o_mo_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 emit_pm25_mo_g(i) =emit_pm25_mo_g(i)+
!      &                            emit_pm25_mo_m(i,m,j)*FAREROT(i,m)
!      &                                           *fcancmxrow(i,m,j)
!                 emit_tpm_mo_g(i) =emit_tpm_mo_g(i)+emit_tpm_mo_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 emit_tc_mo_g(i) =emit_tc_mo_g(i)+emit_tc_mo_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 emit_oc_mo_g(i) =emit_oc_mo_g(i)+emit_oc_mo_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 emit_bc_mo_g(i) =emit_bc_mo_g(i)+emit_bc_mo_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                burnfrac_mo_g(i)=burnfrac_mo_g(i)+burnfrac_mo_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 barefrac=barefrac-FAREROT(i,m)*fcancmxrow(i,m,j)
!
!                end do !j
!
!                 nep_mo_g(i)=nep_mo_g(i)+nep_mo_m(i,m,iccp1)
!      &                          *barefrac
!                 nbp_mo_g(i)=nbp_mo_g(i)+nbp_mo_m(i,m,iccp1)
!      &                          *barefrac
!               hetrores_mo_g(i)=hetrores_mo_g(i)+hetrores_mo_m(i,m,iccp1)
!      &                          *barefrac
!               litres_mo_g(i)  =litres_mo_g(i) +litres_mo_m(i,m,iccp1)
!      &                          *barefrac
!               soilcres_mo_g(i)=soilcres_mo_g(i)+soilcres_mo_m(i,m,iccp1)
!      &                          *barefrac
!
!                luc_emc_mo_g(i) =luc_emc_mo_g(i)
!      &                          +luc_emc_mo_m(i,m)*FAREROT(i,m)
!                lucsocin_mo_g(i) =lucsocin_mo_g(i)
!      &                          +lucsocin_mo_m(i,m)*FAREROT(i,m)
!                lucltrin_mo_g(i) =lucltrin_mo_g(i)
!      &                          +lucltrin_mo_m(i,m)*FAREROT(i,m)
! c    CH4(wetland) variables !Rudra
!
!                ch4wet1_mo_g(i) = ch4wet1_mo_g(i)
!      &                           +ch4wet1_mo_m(i,m)*FAREROT(i,m)
!                ch4wet2_mo_g(i) = ch4wet2_mo_g(i)
!      &                           +ch4wet2_mo_m(i,m)*FAREROT(i,m)
!
!                wetfdyn_mo_m(i,m)=wetfdyn_mo_m(i,m)*(1./
!      &                                      real(monthdays(nt)))
!
!                wetfdyn_mo_g(i) = wetfdyn_mo_g(i)
!      &                           +wetfdyn_mo_m(i,m)*FAREROT(i,m)
!                ch4dyn1_mo_g(i) = ch4dyn1_mo_g(i)
!      &                           +ch4dyn1_mo_m(i,m)*FAREROT(i,m)
!                ch4dyn2_mo_g(i) = ch4dyn2_mo_g(i)
!      &                           +ch4dyn2_mo_m(i,m)*FAREROT(i,m)
!
! !              Make the probability of fire a per day value
!                probfire_mo_m(i,m)=probfire_mo_m(i,m)*
!      &                                (1./real(monthdays(nt)))
!                probfire_mo_g(i)=probfire_mo_g(i)
!      &                          +probfire_mo_m(i,m)*FAREROT(i,m)
!                bterm_mo_m(i,m)=bterm_mo_m(i,m)*(1./real(monthdays(nt)))
!                bterm_mo_g(i) =bterm_mo_g(i)+bterm_mo_m(i,m)*FAREROT(i,m)
!                lterm_mo_m(i,m)=lterm_mo_m(i,m)*(1./real(monthdays(nt)))
!                lterm_mo_g(i) =lterm_mo_g(i)+lterm_mo_m(i,m)*FAREROT(i,m)
!                mterm_mo_m(i,m)=mterm_mo_m(i,m)*(1./real(monthdays(nt)))
!                mterm_mo_g(i) =mterm_mo_g(i)+mterm_mo_m(i,m)*FAREROT(i,m)
!
!
!              endif ! monthend (max lai and accumulated npp/gpp/nep over the whole month)
! c                  ! if(iday.eq.monthend(nt+1))
!
! 865        continue ! nmon
! 863      continue ! m
! c
!          do nt=1,nmon
!            if(iday.eq.monthend(nt+1))then
!              imonth=nt
!
!                 barefrac=1.0
!
! c            Write to file .CT01M_M/.CT01M_G
!               do m=1,nmtest
!                do j=1,icc
!
!                   barefrac=barefrac-fcancmxrow(i,m,j)*FAREROT(i,m)
!
!                 if (FAREROT(i,m)*fcancmxrow(i,m,j) .gt. seed) then
!                  write(84,8104)imonth,iyear,laimaxg_mo_m(i,m,j),
!      1               vgbiomas_mo_m(i,m,j),litrmass_mo_m(i,m,j),
!      2               soilcmas_mo_m(i,m,j),npp_mo_m(i,m,j),
!      3               gpp_mo_m(i,m,j),nep_mo_m(i,m,j),
!      4               nbp_mo_m(i,m,j),hetrores_mo_m(i,m,j),
!      5               autores_mo_m(i,m,j),litres_mo_m(i,m,j),
!      6               soilcres_mo_m(i,m,j),
!      9               ' TILE ',m,' PFT ',j,' FRAC ',FAREROT(i,m)*
!      a               fcancmxrow(i,m,j)
!                 end if
!                end do !icc
!
!                if (m .eq. nmtest) then
!                 if (barefrac .gt. seed) then
!                 write(84,8104)imonth,iyear,0.0,
!      1               0.0,litrmass_mo_m(i,m,iccp1),
!      2               soilcmas_mo_m(i,m,iccp1),0.0,
!      3               0.0,0.0,
!      4               0.0,hetrores_mo_m(i,m,iccp1),
!      5               0.0,litres_mo_m(i,m,iccp1),
!      6               soilcres_mo_m(i,m,iccp1),
!      7               ' TILE ',m,' PFT ',iccp1,' FRAC ',barefrac
!                 end if
!                end if
!               end do !m
!
!               write(84,8104)imonth,iyear,laimaxg_mo_g(i),
!      1                vgbiomas_mo_g(i),litrmass_mo_g(i),
!      2               soilcmas_mo_g(i),npp_mo_g(i),
!      3               gpp_mo_g(i),nep_mo_g(i),
!      4               nbp_mo_g(i),hetrores_mo_g(i),autores_mo_g(i),
!      5               litres_mo_g(i),soilcres_mo_g(i),' GRDAV'
!
!             if (dofire .or. lnduseon) then
!
! c            write to file .CT06M_M/.CT06M_G
!
!               do m=1,nmtest
!                do j=1,icc
!                 if (FAREROT(i,m)*fcancmxrow(i,m,j) .gt. seed) then
!                  write(85,8109)imonth,iyear,emit_co2_mo_m(i,m,j),
!      1               emit_co_mo_m(i,m,j),emit_ch4_mo_m(i,m,j),
!      2               emit_nmhc_mo_m(i,m,j),emit_h2_mo_m(i,m,j),
!      3               emit_nox_mo_m(i,m,j),emit_n2o_mo_m(i,m,j),
!      4               emit_pm25_mo_m(i,m,j),emit_tpm_mo_m(i,m,j),
!      5               emit_tc_mo_m(i,m,j),emit_oc_mo_m(i,m,j),
!      6               emit_bc_mo_m(i,m,j),probfire_mo_m(i,m),
!      7               luc_emc_mo_m(i,m),lucltrin_mo_m(i,m),
!      8               lucsocin_mo_m(i,m),burnfrac_mo_m(i,m,j)*100.,
!      9               bterm_mo_m(i,m),lterm_mo_m(i,m),mterm_mo_m(i,m),
!      &               ' TILE ',m,' PFT ',j,' FRAC ',FAREROT(i,m)*
!      &               fcancmxrow(i,m,j)
!                 end if
!                end do
!               end do
!
!              write(85,8109)imonth,iyear,emit_co2_mo_g(i),
!      3               emit_co_mo_g(i),emit_ch4_mo_g(i),emit_nmhc_mo_g(i),
!      4               emit_h2_mo_g(i),emit_nox_mo_g(i),emit_n2o_mo_g(i),
!      5               emit_pm25_mo_g(i),emit_tpm_mo_g(i),emit_tc_mo_g(i),
!      6               emit_oc_mo_g(i),emit_bc_mo_g(i),
!      7               probfire_mo_g(i),luc_emc_mo_g(i),
!      8               lucltrin_mo_g(i),lucsocin_mo_g(i),
!      8               burnfrac_mo_g(i)*100.,bterm_mo_g(i),lterm_mo_g(i),
!      9               mterm_mo_g(i),' GRDAV '
!
!             endif  !dofire/lnduseon
!
! c           add fraction of each pft and bare \\
! c
!             if (compete .or. lnduseon) then
!               sumfare=0.0
!               if (mosaic) then
!                do m=1,nmos
!                  sumfare=sumfare+FAREROT(i,m)
!                enddo
!                write(88,8106)imonth,iyear,(FAREROT(i,m)*100.,m=1,nmos)
!      1                      ,sumfare,(pftexistrow(i,j,j),j=1,icc)
!               else !composite
!                m=1
!                do j=1,icc
!                  sumfare=sumfare+fcancmxrow(i,m,j)
!                enddo
!                write(88,8106)imonth,iyear,(fcancmxrow(i,m,j)*100.,
!      1                      j=1,icc),(1.0-sumfare)*100.,sumfare,
!      2                        (pftexistrow(i,m,j),j=1,icc)
!               endif !mosaic/composite
!             endif !compete/lnduseon
!
!              if (dowetlands .or. obswetf) then
!              write(91,8111)imonth,iyear,ch4wet1_mo_g(i),
!      1                     ch4wet2_mo_g(i),wetfdyn_mo_g(i),
!      2                     ch4dyn1_mo_g(i),ch4dyn2_mo_g(i)
!              endif
!
! c
! c              initialize monthly accumulated arrays
! c              for the next round
!
!              do m=1,nmtest
!
!                probfire_mo_m(i,m) =0.0
!                luc_emc_mo_m(i,m) =0.0
!                lucsocin_mo_m(i,m) =0.0
!                lucltrin_mo_m(i,m) =0.0
!                bterm_mo_m(i,m) =0.0
!                lterm_mo_m(i,m) =0.0
!                mterm_mo_m(i,m) =0.0
! C       !Rudra
!                ch4wet1_mo_m(i,m)  =0.0
!                ch4wet2_mo_m(i,m)  =0.0
!                wetfdyn_mo_m(i,m)  =0.0
!                ch4dyn1_mo_m(i,m)  =0.0
!                ch4dyn2_mo_m(i,m)  =0.0
!
!
!              do j=1,icc
!
!               laimaxg_mo_m(i,m,j)=0.0
!               hetrores_mo_m(i,m,j)=0.0
!               autores_mo_m(i,m,j)=0.0
!               litres_mo_m(i,m,j)=0.0
!               soilcres_mo_m(i,m,j)=0.0
!
!               npp_mo_m(i,m,j)=0.0
!               gpp_mo_m(i,m,j)=0.0
!               nep_mo_m(i,m,j)=0.0
!               nbp_mo_m(i,m,j)=0.0
!               emit_co2_mo_m(i,m,j)=0.0
!               emit_co_mo_m(i,m,j) =0.0
!               emit_ch4_mo_m(i,m,j) =0.0
!               emit_nmhc_mo_m(i,m,j) =0.0
!               emit_h2_mo_m(i,m,j) =0.0
!               emit_nox_mo_m(i,m,j) =0.0
!               emit_n2o_mo_m(i,m,j) =0.0
!               emit_pm25_mo_m(i,m,j) =0.0
!               emit_tpm_mo_m(i,m,j) =0.0
!               emit_tc_mo_m(i,m,j) =0.0
!               emit_oc_mo_m(i,m,j) =0.0
!               emit_bc_mo_m(i,m,j) =0.0
!               burnfrac_mo_m(i,m,j) =0.0
!              enddo !j
!
!               hetrores_mo_m(i,m,iccp1)=0.0
!               litres_mo_m(i,m,iccp1)=0.0
!               soilcres_mo_m(i,m,iccp1)=0.0
!               nep_mo_m(i,m,iccp1)=0.0
!               nbp_mo_m(i,m,iccp1)=0.0
!
!
!             enddo !m
!
!            endif ! if(iday.eq.monthend(nt+1))
!          enddo ! nt=1,nmon
! c
! 862     continue ! i

        end if !to write out the monthly outputs or not
c
c       accumulate yearly outputs
            call ctem_annual_aw(nltest,nmtest,iday,FAREROT,iyear)
c
!         do 882 i=1,nltest
!           do 883 m=1,nmtest
!             do 884 j=1,icc
!
!              if (ailcgrow(i,m,j).gt.laimaxg_yr_m(i,m,j)) then
!                laimaxg_yr_m(i,m,j)=ailcgrow(i,m,j)
!              end if
!
!             npp_yr_m(i,m,j)=npp_yr_m(i,m,j)+nppvegrow(i,m,j)
!             gpp_yr_m(i,m,j)=gpp_yr_m(i,m,j)+gppvegrow(i,m,j)
!             nep_yr_m(i,m,j)=nep_yr_m(i,m,j)+nepvegrow(i,m,j)
!             nbp_yr_m(i,m,j)=nbp_yr_m(i,m,j)+nbpvegrow(i,m,j)
!             emit_co2_yr_m(i,m,j)=emit_co2_yr_m(i,m,j)+emit_co2row(i,m,j)
!             emit_co_yr_m(i,m,j)=emit_co_yr_m(i,m,j)+emit_corow(i,m,j)
!             emit_ch4_yr_m(i,m,j)=emit_ch4_yr_m(i,m,j)+emit_ch4row(i,m,j)
!             emit_nmhc_yr_m(i,m,j)=emit_nmhc_yr_m(i,m,j)+
!      1                            emit_nmhcrow(i,m,j)
!             emit_h2_yr_m(i,m,j)=emit_h2_yr_m(i,m,j)+emit_h2row(i,m,j)
!             emit_nox_yr_m(i,m,j)=emit_nox_yr_m(i,m,j)+emit_noxrow(i,m,j)
!             emit_n2o_yr_m(i,m,j)=emit_n2o_yr_m(i,m,j)+emit_n2orow(i,m,j)
!             emit_pm25_yr_m(i,m,j)=emit_pm25_yr_m(i,m,j)+
!      1                            emit_pm25row(i,m,j)
!             emit_tpm_yr_m(i,m,j)=emit_tpm_yr_m(i,m,j)+emit_tpmrow(i,m,j)
!             emit_tc_yr_m(i,m,j)=emit_tc_yr_m(i,m,j)+emit_tcrow(i,m,j)
!             emit_oc_yr_m(i,m,j)=emit_oc_yr_m(i,m,j)+emit_ocrow(i,m,j)
!             emit_bc_yr_m(i,m,j)=emit_bc_yr_m(i,m,j)+emit_bcrow(i,m,j)
!
!             hetrores_yr_m(i,m,j)=hetrores_yr_m(i,m,j)
!      &                                +hetroresvegrow(i,m,j)
!             autores_yr_m(i,m,j)=autores_yr_m(i,m,j)
!      &                                +autoresvegrow(i,m,j)
!             litres_yr_m(i,m,j)=litres_yr_m(i,m,j)+litresvegrow(i,m,j)
!             soilcres_yr_m(i,m,j)=soilcres_yr_m(i,m,j)
!      &                                +soilcresvegrow(i,m,j)
!             burnfrac_yr_m(i,m,j)=burnfrac_yr_m(i,m,j)+burnvegfrow(i,m,j)
!
! 884         continue
!
! !           Also do the bare fraction amounts
!             hetrores_yr_m(i,m,iccp1)=hetrores_yr_m(i,m,iccp1)+
!      &                                  hetroresvegrow(i,m,iccp1)
!             litres_yr_m(i,m,iccp1)=litres_yr_m(i,m,iccp1)+
!      &                                  litresvegrow(i,m,iccp1)
!             soilcres_yr_m(i,m,iccp1)=soilcres_yr_m(i,m,iccp1)+
!      &                                  soilcresvegrow(i,m,iccp1)
!             nep_yr_m(i,m,iccp1)=nep_yr_m(i,m,iccp1)+nepvegrow(i,m,iccp1)
!             nbp_yr_m(i,m,iccp1)=nbp_yr_m(i,m,iccp1)+nbpvegrow(i,m,iccp1)
!
!             probfire_yr_m(i,m)=probfire_yr_m(i,m)
!      &                         +(probfirerow(i,m) * (1./365.))
!             bterm_yr_m(i,m)=bterm_yr_m(i,m)+(btermrow(i,m)*(1./365.))
!             lterm_yr_m(i,m)=lterm_yr_m(i,m)+(ltermrow(i,m)*(1./365.))
!             mterm_yr_m(i,m)=mterm_yr_m(i,m)+(mtermrow(i,m)*(1./365.))
!             luc_emc_yr_m(i,m)=luc_emc_yr_m(i,m)+lucemcomrow(i,m)
!             lucsocin_yr_m(i,m)=lucsocin_yr_m(i,m)+lucsocinrow(i,m)
!             lucltrin_yr_m(i,m)=lucltrin_yr_m(i,m)+lucltrinrow(i,m)
! c             CH4(wetland) variables !Rudra
!
!                ch4wet1_yr_m(i,m) = ch4wet1_yr_m(i,m)
!      &                           +ch4wet1row(i,m)
!                ch4wet2_yr_m(i,m) = ch4wet2_yr_m(i,m)
!      &                           +ch4wet2row(i,m)
!                wetfdyn_yr_m(i,m) = wetfdyn_yr_m(i,m)
!      &                           +(wetfdynrow(i,m)*(1./365.))
!                ch4dyn1_yr_m(i,m) = ch4dyn1_yr_m(i,m)
!      &                           +ch4dyn1row(i,m)
!                ch4dyn2_yr_m(i,m) = ch4dyn2_yr_m(i,m)
!      &                           +ch4dyn2row(i,m)
!
!
!             if (iday.eq.365) then
!
!               do 885 j=1,icc
!                 stemmass_yr_m(i,m,j)=stemmassrow(i,m,j)
!                 rootmass_yr_m(i,m,j)=rootmassrow(i,m,j)
!                 litrmass_yr_m(i,m,j)=litrmassrow(i,m,j)
!                 soilcmas_yr_m(i,m,j)=soilcmasrow(i,m,j)
!                 vgbiomas_yr_m(i,m,j)=vgbiomas_vegrow(i,m,j)
!                 totcmass_yr_m(i,m,j)=vgbiomas_yr_m(i,m,j)+
!      &                               litrmass_yr_m(i,m,j)+
!      &                               soilcmas_yr_m(i,m,j)
!
! 885           continue
!
!                 litrmass_yr_m(i,m,iccp1)=litrmassrow(i,m,iccp1)
!                 soilcmas_yr_m(i,m,iccp1)=soilcmasrow(i,m,iccp1)
!
!                 barefrac=1.0
!
!               do j=1,icc
!                 laimaxg_yr_g(i)=laimaxg_yr_g(i)+ laimaxg_yr_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 stemmass_yr_g(i)=stemmass_yr_g(i)+stemmass_yr_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 rootmass_yr_g(i)=rootmass_yr_g(i)+rootmass_yr_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 litrmass_yr_g(i)=litrmass_yr_g(i)+litrmass_yr_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 soilcmas_yr_g(i)=soilcmas_yr_g(i)+soilcmas_yr_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 vgbiomas_yr_g(i)=vgbiomas_yr_g(i)+vgbiomas_yr_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 totcmass_yr_g(i)=totcmass_yr_g(i)+totcmass_yr_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 npp_yr_g(i)=npp_yr_g(i)+npp_yr_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 gpp_yr_g(i)=gpp_yr_g(i)+gpp_yr_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 nep_yr_g(i)=nep_yr_g(i)+nep_yr_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 nbp_yr_g(i)=nbp_yr_g(i)+nbp_yr_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 emit_co2_yr_g(i)=emit_co2_yr_g(i)+emit_co2_yr_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 emit_co_yr_g(i)=emit_co_yr_g(i)+emit_co_yr_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 emit_ch4_yr_g(i)=emit_ch4_yr_g(i)+emit_ch4_yr_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                emit_nmhc_yr_g(i)=emit_nmhc_yr_g(i)+emit_nmhc_yr_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 emit_h2_yr_g(i)=emit_h2_yr_g(i)+emit_h2_yr_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 emit_nox_yr_g(i)=emit_nox_yr_g(i)+emit_nox_yr_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 emit_n2o_yr_g(i)=emit_n2o_yr_g(i)+emit_n2o_yr_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                emit_pm25_yr_g(i)=emit_pm25_yr_g(i)+emit_pm25_yr_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 emit_tpm_yr_g(i)=emit_tpm_yr_g(i)+emit_tpm_yr_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 emit_tc_yr_g(i)=emit_tc_yr_g(i)+emit_tc_yr_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 emit_oc_yr_g(i)=emit_oc_yr_g(i)+emit_oc_yr_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 emit_bc_yr_g(i)=emit_bc_yr_g(i)+emit_bc_yr_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!
!                 hetrores_yr_g(i)=hetrores_yr_g(i)+hetrores_yr_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 autores_yr_g(i) =autores_yr_g(i) +autores_yr_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 litres_yr_g(i)  =litres_yr_g(i)  +litres_yr_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 soilcres_yr_g(i) =soilcres_yr_g(i) +soilcres_yr_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!                 burnfrac_yr_g(i)=burnfrac_yr_g(i)+burnfrac_yr_m(i,m,j)
!      &                          *FAREROT(i,m)*fcancmxrow(i,m,j)
!
!                 barefrac=barefrac-FAREROT(i,m)*fcancmxrow(i,m,j)
!
!               end do !j
!
!               litrmass_yr_g(i)=litrmass_yr_g(i)+litrmass_yr_m(i,m,iccp1)
!      &                          *barefrac
!               soilcmas_yr_g(i)=soilcmas_yr_g(i)+soilcmas_yr_m(i,m,iccp1)
!      &                          *barefrac
!               hetrores_yr_g(i)=hetrores_yr_g(i)+hetrores_yr_m(i,m,iccp1)
!      &                          *barefrac
!               litres_yr_g(i)  =litres_yr_g(i)  +litres_yr_m(i,m,iccp1)
!      &                          *barefrac
!               soilcres_yr_g(i)=soilcres_yr_g(i)+soilcres_yr_m(i,m,iccp1)
!      &                          *barefrac
!                 nep_yr_g(i)=nep_yr_g(i)+nep_yr_m(i,m,j)
!      &                          *barefrac
!                 nbp_yr_g(i)=nbp_yr_g(i)+nbp_yr_m(i,m,j)
!      &                          *barefrac
!
!
!               probfire_yr_g(i)=probfire_yr_g(i)
!      &                         +probfire_yr_m(i,m)*FAREROT(i,m)
!               bterm_yr_g(i)=bterm_yr_g(i)+bterm_yr_m(i,m)*FAREROT(i,m)
!               lterm_yr_g(i)=lterm_yr_g(i)+lterm_yr_m(i,m)*FAREROT(i,m)
!               mterm_yr_g(i)=mterm_yr_g(i)+mterm_yr_m(i,m)*FAREROT(i,m)
!               luc_emc_yr_g(i)=luc_emc_yr_g(i)
!      &                         +luc_emc_yr_m(i,m)*FAREROT(i,m)
!               lucsocin_yr_g(i)=lucsocin_yr_g(i)
!      &                         +lucsocin_yr_m(i,m)*FAREROT(i,m)
!               lucltrin_yr_g(i)=lucltrin_yr_g(i)
!      &                         +lucltrin_yr_m(i,m)*FAREROT(i,m)
! c    CH4(wetland) variables !Rudra
!
!                ch4wet1_yr_g(i) = ch4wet1_yr_g(i)
!      &                           +ch4wet1_yr_m(i,m)*FAREROT(i,m)
!                ch4wet2_yr_g(i) = ch4wet2_yr_g(i)
!      &                           +ch4wet2_yr_m(i,m)*FAREROT(i,m)
!                wetfdyn_yr_g(i) = wetfdyn_yr_g(i)
!      &                           +wetfdyn_yr_m(i,m)*FAREROT(i,m)
!                ch4dyn1_yr_g(i) = ch4dyn1_yr_g(i)
!      &                           +ch4dyn1_yr_m(i,m)*FAREROT(i,m)
!                ch4dyn2_yr_g(i) = ch4dyn2_yr_g(i)
!      &                           +ch4dyn2_yr_m(i,m)*FAREROT(i,m)
!
! c
!             endif ! iday 365
! c
! 883       continue ! m
! c
!           if (iday.eq.365) then
!
!               barefrac=1.0
!
! c            Write to file .CT01Y_M/.CT01Y_G
!
!            do m=1,nmtest
!             do j=1,icc
!
!                 barefrac=barefrac-fcancmxrow(i,m,j)*FAREROT(i,m)
!
!              if (FAREROT(i,m)*fcancmxrow(i,m,j) .gt. seed) then
!               write(86,8105)iyear,laimaxg_yr_m(i,m,j),
!      1            vgbiomas_yr_m(i,m,j),stemmass_yr_m(i,m,j),
!      2            rootmass_yr_m(i,m,j),litrmass_yr_m(i,m,j),
!      3            soilcmas_yr_m(i,m,j),totcmass_yr_m(i,m,j),
!      4            npp_yr_m(i,m,j),gpp_yr_m(i,m,j),nep_yr_m(i,m,j),
!      5            nbp_yr_m(i,m,j),hetrores_yr_m(i,m,j),
!      6            autores_yr_m(i,m,j),litres_yr_m(i,m,j),
!      9            soilcres_yr_m(i,m,j),' TILE ',m,' PFT ',j,' FRAC '
!      a            ,FAREROT(i,m)*fcancmxrow(i,m,j)
!              end if
!             end do !j
!
! !           Now do the bare fraction of the grid cell. Only soil c, hetres
! !           and litter are relevant so the rest are set to 0.
!             if (m .eq. nmtest) then
!              if (barefrac .gt. seed) then
!               write(86,8105)iyear,0.,
!      1            0.,
!      1            0.,0.,
!      2            litrmassrow(i,m,iccp1),soilcmasrow(i,m,iccp1),
!      3            0.+soilcmasrow(i,m,iccp1)+
!      2            litrmassrow(i,m,iccp1),0.,
!      4            0.,0.,
!      5            0.,hetrores_yr_m(i,m,iccp1),
!      6            0.,litres_yr_m(i,m,iccp1),soilcres_yr_m(i,m,iccp1),
!      7            ' TILE ',m,' PFT ',iccp1,' FRAC ',barefrac
!              end if
!             end if
!
!            end do !m
!
!              write(86,8105)iyear,laimaxg_yr_g(i),vgbiomas_yr_g(i),
!      1            stemmass_yr_g(i),rootmass_yr_g(i),litrmass_yr_g(i),
!      2            soilcmas_yr_g(i),totcmass_yr_g(i),npp_yr_g(i),
!      3            gpp_yr_g(i),nep_yr_g(i),
!      4            nbp_yr_g(i),hetrores_yr_g(i),autores_yr_g(i),
!      5            litres_yr_g(i),soilcres_yr_g(i),' GRDAV'
!
!            if (dofire .or. lnduseon) then
! c            Write to file .CT06Y_M/.CT06Y_G
!             do m=1,nmtest
!              do j=1,icc
!               if (FAREROT(i,m)*fcancmxrow(i,m,j) .gt. seed) then
!                write(87,8108)iyear,emit_co2_yr_m(i,m,j),
!      1            emit_co_yr_m(i,m,j),emit_ch4_yr_m(i,m,j),
!      2            emit_nmhc_yr_m(i,m,j),emit_h2_yr_m(i,m,j),
!      3            emit_nox_yr_m(i,m,j),emit_n2o_yr_m(i,m,j),
!      4            emit_pm25_yr_m(i,m,j),emit_tpm_yr_m(i,m,j),
!      5            emit_tc_yr_m(i,m,j),emit_oc_yr_m(i,m,j),
!      6            emit_bc_yr_m(i,m,j),probfire_yr_m(i,m),
!      7            luc_emc_yr_m(i,m),lucltrin_yr_m(i,m),
!      8            lucsocin_yr_m(i,m),burnfrac_yr_m(i,m,j)*100.,
!      9            bterm_yr_m(i,m),lterm_yr_m(i,m),mterm_yr_m(i,m),
!      9            ' TILE ',m,' PFT ',j,' FRAC '
!      a            ,FAREROT(i,m)*fcancmxrow(i,m,j)
!              end if
!             end do
!            end do
!
!              write(87,8108)iyear,emit_co2_yr_g(i),
!      4            emit_co_yr_g(i),emit_ch4_yr_g(i),emit_nmhc_yr_g(i),
!      5            emit_h2_yr_g(i),emit_nox_yr_g(i),emit_n2o_yr_g(i),
!      6            emit_pm25_yr_g(i),emit_tpm_yr_g(i),emit_tc_yr_g(i),
!      7            emit_oc_yr_g(i),emit_bc_yr_g(i),probfire_yr_g(i),
!      8            luc_emc_yr_g(i),lucltrin_yr_g(i),
!      9            lucsocin_yr_g(i),burnfrac_yr_g(i)*100.,bterm_yr_g(i),
!      a            lterm_yr_g(i),mterm_yr_g(i), ' GRDAV'
!
!            endif !dofire,lnduseon
!
! c           write fraction of each pft and bare
! c
!              if (compete .or. lnduseon) then
!                  sumfare=0.0
!                if (mosaic) then
!                  do m=1,nmos
!                     sumfare=sumfare+FAREROT(i,m)
!                  enddo
!                  write(89,8107)iyear,(FAREROT(i,m)*100.,m=1,nmos),
!      &                         sumfare,(pftexistrow(i,j,j),j=1,icc)
!                else  !composite
!                  m=1
!                  do j=1,icc
!                     sumfare=sumfare+fcancmxrow(i,m,j)
!                  enddo
!                 write(89,8107)iyear,(fcancmxrow(i,m,j)*100.,
!      1                       j=1,icc),(1.0-sumfare)*100.,sumfare,
!      2                      (pftexistrow(i,m,j),j=1,icc)
!
!                endif
!              endif !compete/lnduseon
! C
!               if (dowetlands .or. obswetf) then
!                 write(92,8115)iyear,ch4wet1_yr_g(i),
!      1                     ch4wet2_yr_g(i),wetfdyn_yr_g(i),
!      2                     ch4dyn1_yr_g(i),ch4dyn2_yr_g(i)
!               endif
!
!
!
! c             initialize yearly accumulated arrays
! c             for the next round
!
!              do m=1,nmtest
!               probfire_yr_m(i,m)=0.0
!               luc_emc_yr_m(i,m)=0.0
!               lucsocin_yr_m(i,m)=0.0
!               lucltrin_yr_m(i,m)=0.0
!               bterm_yr_m(i,m)=0.0
!               lterm_yr_m(i,m)=0.0
!               mterm_yr_m(i,m)=0.0
! C       !Rudra
!                ch4wet1_yr_m(i,m)  =0.0
!                ch4wet2_yr_m(i,m)  =0.0
!                wetfdyn_yr_m(i,m)  =0.0
!                ch4dyn1_yr_m(i,m)  =0.0
!                ch4dyn2_yr_m(i,m)  =0.0
!
!
!                do j = 1, icc
!                 laimaxg_yr_m(i,m,j)=0.0
!                 npp_yr_m(i,m,j)=0.0
!                 gpp_yr_m(i,m,j)=0.0
!                 nep_yr_m(i,m,j)=0.0
!                 nbp_yr_m(i,m,j)=0.0
!                 emit_co2_yr_m(i,m,j)=0.0
!                 emit_co_yr_m(i,m,j)=0.0
!                 emit_ch4_yr_m(i,m,j)=0.0
!                 emit_nmhc_yr_m(i,m,j)=0.0
!                 emit_h2_yr_m(i,m,j)=0.0
!                 emit_nox_yr_m(i,m,j)=0.0
!                 emit_n2o_yr_m(i,m,j)=0.0
!                 emit_pm25_yr_m(i,m,j)=0.0
!                 emit_tpm_yr_m(i,m,j)=0.0
!                 emit_tc_yr_m(i,m,j)=0.0
!                 emit_oc_yr_m(i,m,j)=0.0
!                 emit_bc_yr_m(i,m,j)=0.0
!                 hetrores_yr_m(i,m,j)=0.0
!                 autores_yr_m(i,m,j)=0.0
!                 litres_yr_m(i,m,j)=0.0
!                 soilcres_yr_m(i,m,j)=0.0
!                 burnfrac_yr_m(i,m,j)=0.0
!
!                enddo
!                 hetrores_yr_m(i,m,iccp1)=0.0
!                 litres_yr_m(i,m,iccp1)=0.0
!                 soilcres_yr_m(i,m,iccp1)=0.0
!                 nep_yr_m(i,m,iccp1)=0.0
!                 nbp_yr_m(i,m,iccp1)=0.0
!              enddo

!            endif ! if iday=365
!c
!882     continue ! i
c
      endif ! if(ncount.eq.nday)
      endif ! if(ctem_on)
C
! 8104  FORMAT(1X,I4,I5,12(F10.3,1X),2(A6,I2),A6,F8.2)
! 8105  FORMAT(1X,I5,15(F10.3,1X),2(A6,I2),A6,F8.2)
! 8106  FORMAT(1X,I4,I5,11(F10.5,1X),9L5,2(A6,I2))
! 8107  FORMAT(1X,I5,11(F10.5,1X),9L5,2(A6,I2))
! 8108  FORMAT(1X,I5,20(F10.3,1X),2(A6,I2),A6,F8.2)
! 8109  FORMAT(1X,I4,I5,20(F10.3,1X),2(A6,I2),A6,F8.2)
! 8111  FORMAT(1X,I4,I5,5(F10.3,1X),2(A6,I2))
! 8115  FORMAT(1X,I5,5(F10.3,1X),2(A6,I2))
C
C     OPEN AND WRITE TO THE RESTART FILES
C
      IF (RSFILE) THEN
       IF (IDAY.EQ.365.AND.NCOUNT.EQ.NDAY) THEN
C
C       WRITE .INI_RS FOR CLASS RESTART DATA
C
        OPEN(UNIT=100,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.INI_RS')
C
        WRITE(100,5010) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
        WRITE(100,5010) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
        WRITE(100,5010) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
C
        WRITE(100,5020)DLATROW(1),DEGLON,ZRFMROW(1),ZRFHROW(1),
     1                 ZBLDROW(1),GCROW(1),NLTEST,NMTEST
        DO I=1,NLTEST
          DO M=1,NMTEST

C         IF START_BARE (SO EITHER COMPETE OR LNDUSEON), THEN WE NEED TO CREATE
C         THE FCANROT FOR THE RS FILE.
          IF (START_BARE .AND. MOSAIC) THEN
           IF (M .LE. 2) THEN                     !NDL
            FCANROT(I,M,1)=1.0
           ELSEIF (M .GE. 3 .AND. M .LE. 5) THEN  !BDL
            FCANROT(I,M,2)=1.0
           ELSEIF (M .EQ. 6 .OR. M .EQ. 7) THEN  !CROP
            FCANROT(I,M,3)=1.0
           ELSEIF (M .EQ. 8 .OR. M .EQ. 9) THEN  !GRASSES
            FCANROT(I,M,4)=1.0
           ELSE                                  !BARE
            FCANROT(I,M,5)=1.0
           ENDIF
          ENDIF !START_BARE/MOSAIC

            WRITE(100,5040) (FCANROT(I,M,J),J=1,ICAN+1),(PAMXROT(I,M,J),
     1                      J=1,ICAN)
            WRITE(100,5040) (LNZ0ROT(I,M,J),J=1,ICAN+1),(PAMNROT(I,M,J),
     1                      J=1,ICAN)
            WRITE(100,5040) (ALVCROT(I,M,J),J=1,ICAN+1),(CMASROT(I,M,J),
     1                      J=1,ICAN)
            WRITE(100,5040) (ALICROT(I,M,J),J=1,ICAN+1),(ROOTROT(I,M,J),
     1                      J=1,ICAN)
            WRITE(100,5030) (RSMNROT(I,M,J),J=1,ICAN),
     1                      (QA50ROT(I,M,J),J=1,ICAN)
            WRITE(100,5030) (VPDAROT(I,M,J),J=1,ICAN),
     1                      (VPDBROT(I,M,J),J=1,ICAN)
            WRITE(100,5030) (PSGAROT(I,M,J),J=1,ICAN),
     1                      (PSGBROT(I,M,J),J=1,ICAN)
            WRITE(100,5040) DRNROT(I,M),SDEPROT(I,M),FAREROT(I,M)
            WRITE(100,5090) XSLPROT(I,M),GRKFROT(I,M),WFSFROT(I,M),
     1                      WFCIROT(I,M),MIDROT(I,M)
            WRITE(100,5080) (SANDROT(I,M,J),J=1,3)
            WRITE(100,5080) (CLAYROT(I,M,J),J=1,3)
            WRITE(100,5080) (ORGMROT(I,M,J),J=1,3)
C           Temperatures are in degree C
            IF (TCANROT(I,M).NE.0.0) TCANRS(I,M)=TCANROT(I,M)-273.16
            IF (TSNOROT(I,M).NE.0.0) TSNORS(I,M)=TSNOROT(I,M)-273.16
            IF (TPNDROT(I,M).NE.0.0) TPNDRS(I,M)=TPNDROT(I,M)-273.16
            WRITE(100,5050) (TBARROT(I,M,J)-273.16,J=1,3),TCANRS(I,M),
     2                      TSNORS(I,M),TPNDRS(I,M)
            WRITE(100,5060) (THLQROT(I,M,J),J=1,3),(THICROT(I,M,J),
     1                      J=1,3),ZPNDROT(I,M)
C
            WRITE(100,5070) RCANROT(I,M),SCANROT(I,M),SNOROT(I,M),
     1                      ALBSROT(I,M),RHOSROT(I,M),GROROT(I,M)
C           WRITE(100,5070) 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
          ENDDO
        ENDDO
C
        DO J=1,IGND
          WRITE(100,5002) DELZ(J),ZBOT(J)
        ENDDO
C
        WRITE(100,5200) JHHSTD,JHHENDD,JDSTD,JDENDD
        WRITE(100,5200) JHHSTY,JHHENDY,JDSTY,JDENDY
        CLOSE(100)
C
c       write .CTM_RS for ctem restart data
c
         if (ctem_on) then
            call write_ctm_rs(nltest,nmtest,FCANROT,argbuff)
!           open(unit=101,file=argbuff(1:strlen(argbuff))//'.CTM_RS')
! c
!           write(101,7010) titlec1
!           write(101,7010) titlec2
!           write(101,7010) titlec3
!
! c         if landuseon or competition, then we need to recreate the
! c         dvdfcanrow so do so now
!           if (lnduseon .or. compete ) then
!            icountrow=0
!            do j = 1, ican
!             do i=1,nltest
!              do m=1,nmtest
!                k1c = (j-1)*l2max + 1
!                k2c = k1c + (l2max - 1)
! c
!                 do n = k1c, k2c
!                  if(modelpft(n).eq.1)then
!                   icountrow(i,m) = icountrow(i,m) + 1
!                   if (FCANROT(i,m,j) .gt. 0.) then
!                    dvdfcanrow(i,m,icountrow(i,m))=
!      1                 fcancmxrow(i,m,icountrow(i,m))/FCANROT(i,m,j)
!                   else
!                    dvdfcanrow(i,m,icountrow(i,m))=0.
!                   end if
!                  endif !modelpft
!                 enddo !n
!
! c            check to ensure that the dvdfcanrow's add up to 1 across a
! c            class-level pft
!              if(dvdfcanrow(i,m,1) .eq. 0. .and. dvdfcanrow(i,m,2)
!      1                                              .eq. 0.) then
!                   dvdfcanrow(i,m,1)=1.0
!              elseif(dvdfcanrow(i,m,3) .eq. 0. .and. dvdfcanrow(i,m,4)
!      1                     .eq. 0. .and. dvdfcanrow(i,m,5) .eq. 0.) then
!                   dvdfcanrow(i,m,3)=1.0
!              elseif(dvdfcanrow(i,m,6) .eq. 0. .and. dvdfcanrow(i,m,7)
!      1                                              .eq. 0.) then
!                   dvdfcanrow(i,m,6)=1.0
!              elseif(dvdfcanrow(i,m,8) .eq. 0. .and. dvdfcanrow(i,m,9)
!      1                                              .eq. 0.) then
!                   dvdfcanrow(i,m,8)=1.0
!              endif
!
!              enddo !m
!             enddo !i
!            enddo !j
!
!            do i=1,nltest
!             do m=1,nmtest
!              do j = 1, icc
! c            lastly check if the different pfts accidently add up > 1.0
! c            after rounding to the number of sig figs used in the output
! c            this rounds to 3 decimal places. if you are found to be over
! c            or under, arbitrarily reduce one of the pfts. the amount of
! c            the change will be inconsequential.
!               rnded_pft(j) =real(int(dvdfcanrow(i,m,j) * 1000.0))
!      1                                                         / 1000.0
!
!               dvdfcanrow(i,m,j) = rnded_pft(j)
!              enddo
!
!              if (dvdfcanrow(i,m,1) + dvdfcanrow(i,m,2) .ne. 1.0) then
!               dvdfcanrow(i,m,1) = 1.0 - rnded_pft(2)
!               dvdfcanrow(i,m,2) = rnded_pft(2)
!              end if
!              if (dvdfcanrow(i,m,3) + dvdfcanrow(i,m,4) +
!      1                                 dvdfcanrow(i,m,5) .ne. 1.0) then
!               dvdfcanrow(i,m,3) = 1.0 - rnded_pft(4) - rnded_pft(5)
!               dvdfcanrow(i,m,4) = rnded_pft(4)
!               dvdfcanrow(i,m,5) = rnded_pft(5)
!              end if
!              if (dvdfcanrow(i,m,6) + dvdfcanrow(i,m,7) .ne. 1.0) then
!               dvdfcanrow(i,m,6) = 1.0 - rnded_pft(7)
!               dvdfcanrow(i,m,7) = rnded_pft(7)
!              end if
!              if (dvdfcanrow(i,m,8) + dvdfcanrow(i,m,9) .ne. 1.0) then
!               dvdfcanrow(i,m,8) = 1.0 - rnded_pft(9)
!               dvdfcanrow(i,m,9) = rnded_pft(9)
!              endif
!             enddo
!            enddo
!
!           endif !lnuse/compete
! c
!           do i=1,nltest
!             do m=1,nmtest
!               write(101,7011) (ailcminrow(i,m,j),j=1,icc)
!               write(101,7011) (ailcmaxrow(i,m,j),j=1,icc)
!               write(101,'(9f8.3)') (dvdfcanrow(i,m,j),j=1,icc)
! c
!               write(101,7011) (gleafmasrow(i,m,j),j=1,icc)
!               write(101,7011) (bleafmasrow(i,m,j),j=1,icc)
!               write(101,7011) (stemmassrow(i,m,j),j=1,icc)
!               write(101,7011) (rootmassrow(i,m,j),j=1,icc)
!               write(101,7013) (litrmassrow(i,m,j),j=1,iccp1)
!               write(101,7013) (soilcmasrow(i,m,j),j=1,iccp1)
!               write(101,7012) (lfstatusrow(i,m,j),j=1,icc)
!               write(101,7012) (pandaysrow(i,m,j),j=1,icc)
!             enddo
! c
!             write(101,"(6f8.3)") (mlightnggrd(i,j),j=1,6)  !mean monthly lightning frequency
!             write(101,"(6f8.3)") (mlightnggrd(i,j),j=7,12) !flashes/km2.year
!             write(101,"(f8.2)") extnprobgrd(i)
!             write(101,"(f8.2)") prbfrhucgrd(i)
!             write(101,"(i4)") stdalngrd(i)
!
!             if (compete) then
!              write(101,"(5f8.2)")twarmm(i),tcoldm(i),gdd5(i),
!      1                            aridity(i),srplsmon(i)
!              write(101,"(5f8.2)")defctmon(i),anndefct(i),annsrpls(i),
!      1                        annpcp(i),dry_season_length(i)
!             endif
!
!             if (dowetlands) then
!               write(101,"(8f9.5)")(wetfrac_sgrd(i,j),j=1,8)
!             endif
!
!           enddo
! c
!
!           close(101)
        endif ! ctem_on
c
       endif ! if iday=365
      endif ! if generate restart files
c
! 7011  format(9f8.2)
! 7012  format(9i8)
! 7013  format(10f8.2)
c
c      check if the model is done running.
       if (iday.eq.365.and.ncount.eq.nday) then

          if (cyclemet .and. climiyear .ge. metcycendyr) then

            lopcount = lopcount+1

             if(lopcount.le.ctemloop .and. .not. transient_run)then

              rewind(12)   ! rewind met file
c /---------------------Rudra----------------/
               if(obswetf) then
                rewind(16) !rewind obswetf file
                read(16,*) ! read in the header
               endif
c\----------------------Rudra---------------\
              met_rewound = .true.
              iyear=-9999
              obswetyr=-9999     !Rudra

               if(popdon) then
                 rewind(13) !rewind popd file
                 read(13,*) ! skip header (3 lines)
                 read(13,*) ! skip header (3 lines)
                 read(13,*) ! skip header (3 lines)
               endif
               if(co2on) then
                 rewind(14) !rewind co2 file
               endif

             else if (lopcount.le.ctemloop .and. transient_run)then
             ! rewind only the MET file (since we are looping over the MET  while
             ! the other inputs continue on.
               rewind(12)   ! rewind met file

             else

              if (transient_run .and. cyclemet) then
              ! Now switch from cycling over the MET to running through the file
              rewind(12)   ! rewind met file
              cyclemet = .false.
              lopcount = 1
              endyr = iyear + ncyear  !set the new end year

              else
               run_model = .false.
              endif

             endif

          else if (iyear .eq. endyr .and. .not. cyclemet) then

             run_model = .false.

          endif !if cyclemet and iyear > metcycendyr
       endif !last day of year check

C===================== CTEM =====================================/
C
        NCOUNT=NCOUNT+1
        IF(NCOUNT.GT.NDAY) THEN
            NCOUNT=1
        ENDIF

      ENDDO !MAIN MODEL LOOP

C     MODEL RUN HAS COMPLETED SO NOW CLOSE OUTPUT FILES AND EXIT
C==================================================================
c
c      checking the time spent for running model
c
c      call idate(today)
c      call itime(now)
c      write(*,1001) today(2), today(1), 2000+today(3), now
c 1001 format( 'end date: ', i2.2, '/', i2.2, '/', i4.4,
c     &      '; end time: ', i2.2, ':', i2.2, ':', i2.2 )
c
      IF (.NOT. PARALLELRUN) THEN
C       FIRST ANY CLASS OUTPUT FILES
        CLOSE(61)
        CLOSE(62)
        CLOSE(63)
        CLOSE(64)
        CLOSE(65)
        CLOSE(66)
        CLOSE(67)
        CLOSE(68)
        CLOSE(69)
        CLOSE(611)
        CLOSE(621)
        CLOSE(631)
        CLOSE(641)
        CLOSE(651)
        CLOSE(661)
        CLOSE(671)
        CLOSE(681)
        CLOSE(691)
        end if ! moved this up from below so it calls the close subroutine. JRM.

c       then ctem ones

        call close_outfiles()

!         close(71)
!         close(711)
!         close(721)
!         close(731)
!         close(741)
!         close(751)
!
!         if (mosaic) then
!          close(72)
!          close(73)
!          close(74)
!          close(75)
!          close(76)
!          if (dofire .or. lnduseon) then
!           close(78)
!          end if
!         end if
! c
!         if (compete .or. lnduseon) then
!           close(761)
!         endif
! c
!         if (dowetlands .or. obswetf) then
!         close(762)
!         endif
! c
!        if (dofire .or. lnduseon) then
!         close(781)
!        endif
!       endif ! if (.not. parallelrun)
! C
! C     CLOSE CLASS OUTPUT FILES
!       CLOSE(81)
!       CLOSE(82)
!       CLOSE(83)
! c     then the CTEM ones
!       close(84)
!       close(86)
!       if (dofire .or. lnduseon) then
!        close(85)
!        close(87)
!       endif
!       if (compete .or. lnduseon) then
!        close(88)
!        close(89)
!       endif
!
!       if (dowetlands .or. obswetf) then
!        close(91)
!        close(92)
!       endif
c
c     close the input files too
      close(12)
      close(13)
      close(14)
      if (obswetf) then
        close(16)  !*.WET
      end if
      call exit
C
c         the 999 label below is hit when an input file reaches its end.
999       continue

            lopcount = lopcount+1

             if(lopcount.le.ctemloop)then

              rewind(12)   ! rewind met file
c /-----------Rudra-----------------/
                if(obswetf) then
                  rewind(16) !rewind obswetf file
                  read(16,*) ! read in the header
                endif
c \------------Rudra---------------\
              met_rewound = .true.
              iyear=-9999
              obswetyr=-9999   !Rudra

               if(popdon) then
                 rewind(13) !rewind popd file
                 read(13,*) ! skip header (3 lines)
                 read(13,*) ! skip header
                 read(13,*) ! skip header
               endif
               if(co2on) then
                 rewind(14) !rewind co2 file
               endif

             else

              run_model = .false.

             endif

c     return to the time stepping loop
      if (run_model) then
         goto 200
      else

c     close the output files
C
      IF (.NOT. PARALLELRUN) THEN
C       FIRST ANY CLASS OUTPUT FILES
        CLOSE(61)
        CLOSE(62)
        CLOSE(63)
        CLOSE(64)
        CLOSE(65)
        CLOSE(66)
        CLOSE(67)
        CLOSE(68)
        CLOSE(69)
        CLOSE(611)
        CLOSE(621)
        CLOSE(631)
        CLOSE(641)
        CLOSE(651)
        CLOSE(661)
        CLOSE(671)
        CLOSE(681)
        CLOSE(691)
        end if ! moved this up from below so it calls the close subroutine. JRM.
c       then ctem ones
        call close_outfiles()
!         close(71)
!         close(711)
!         close(721)
!         close(731)
!         close(741)
!         close(751)
!
!         if (mosaic) then
!          close(72)
!          close(73)
!          close(74)
!          close(75)
!          close(76)
!          if (dofire .or. lnduseon) then
!           close(78)
!          end if
!         end if
!
! c
!         if (compete .or. lnduseon) then
!           close(761)
!         endif
! c
!        if (dofire .or. lnduseon) then
!         close(781)
!        endif
!       endif ! if (.not. parallelrun)
! C
! C     CLOSE CLASS OUTPUT FILES
!       CLOSE(81)
!       CLOSE(82)
!       CLOSE(83)
! c     then the CTEM ones
!       close(84)
!       close(86)
!       if (dofire .or. lnduseon) then
!        close(85)
!        close(87)
!       endif
!       if (compete .or. lnduseon) then
!        close(88)
!        close(89)
!       endif
!       if (dowetlands .or. obswetf) then
!        close(91)
!        close(92)
!       endif

C
C     CLOSE THE INPUT FILES TOO
      CLOSE(12)
      CLOSE(13)
      CLOSE(14)

         CALL EXIT
      END IF

1001  continue

       write(*,*)'Error while reading WETF file'
       run_model=.false.



C ============================= CTEM =========================/

      END

      INTEGER FUNCTION STRLEN(ST)
      INTEGER       I
      CHARACTER     ST*(*)
      I = LEN(ST)
      DO WHILE (ST(I:I) .EQ. ' ')
        I = I - 1
      ENDDO
      STRLEN = I
      RETURN
      END
