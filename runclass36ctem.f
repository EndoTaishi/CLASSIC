      PROGRAM RUNCLASS36CTEM
C
C     * SHELL PROGRAM TO RUN "CLASS" ("CANADIAN LAND SURFACE SCHEME")
C     * VERSION 3.6 IN STAND-ALONE MODE USING SPECIFIED BOUNDARY 
C     * CONDITIONS AND ATMOSPHERIC FORCING, COUPLED TO CTEM (CANADIAN TERRESTRIAL 
C     * ECOSYSTEM MODEL.
C
C     REVISION HISTORY:
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
C     * THE SUFFIX "ROW" REFERS TO VARIABLES EXISTING ON THE 
C     * MOSAIC GRID ON THE CURRENT LATITUDE CIRCLE.  THE SUFFIX 
C     * "GAT" REFERS TO THE SAME VARIABLES AFTER THEY HAVE UNDERGONE
C     * A "GATHER" OPERATION IN WHICH THE TWO MOSAIC DIMENSIONS
C     * ARE COLLAPSED INTO ONE.  THE SUFFIX "GRD" REFERS BOTH TO 
C     * GRID-CONSTANT INPUT VARIABLES. AND TO GRID-AVERAGED
C     * DIAGNOSTIC VARIABLES.
C     
C     * THE FIRST DIMENSION ELEMENT OF THE "ROW" VARIABLES 
C     * REFERS TO THE NUMBER OF GRID CELLS ON THE CURRENT 
C     * LATITUDE CIRCLE.  IN THIS STAND-ALONE VERSION, THIS 
C     * NUMBER IS ARBITRARILY SET TO THREE, TO ALLOW UP TO THREE
C     * SIMULTANEOUS TESTS TO BE RUN.  THE SECOND DIMENSION 
C     * ELEMENT OF THE "ROW" VARIABLES REFERS TO THE MAXIMUM
C     * NUMBER OF TILES IN THE MOSAIC.  IN THIS STAND-ALONE
C     * VERSION, THIS NUMBER IS SET TO EIGHT.  THE FIRST 
C     * DIMENSION ELEMENT IN THE "GAT" VARIABLES IS GIVEN BY
C     * THE PRODUCT OF THE FIRST TWO DIMENSION ELEMENTS IN THE
C     * "ROW" VARIABLES.

c     use statements for modules:
      use ctem_params,        only : initpftpars, nlat, nmos, ilg, nmon, 
     1                               ican, ignd,icp1, icc, iccp1, 
     2                               monthend, mmday,modelpft, l2max,
     3                                deltat, abszero, monthdays,seed,
     4                                crop
     
      use landuse_change,     only : initialize_luc, readin_luc
      
c
      implicit none
C
C     * INTEGER CONSTANTS.
C
      INTEGER IDISP,IZREF,ISLFD,IPCP,IWF,IPAI,IHGT,IALC,
     1        IALS,IALG,N,ITG,ITC,ITCG

      INTEGER NLTEST,NMTEST,NCOUNT,NDAY,
     1        IMONTH,NDMONTH,NT,
     2        IHOUR,IMIN,IDAY,IYEAR,NML,NMW,NWAT,NICE,JLAT,
     3        NLANDCS,NLANDGS,NLANDC,NLANDG,NLANDI,I,J,K,L,M
C
      INTEGER K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11
C
      INTEGER*4 TODAY(3), NOW(3)
C
C     * LAND SURFACE PROGNOSTIC VARIABLES.
C
      REAL,DIMENSION(NLAT,NMOS,IGND) ::
     1        TBARROW,   THLQROW,   THICROW
C
      REAL,DIMENSION(NLAT,NMOS) ::
     1        TPNDROW,   ZPNDROW,   TBASROW,   
     2        ALBSROW,   TSNOROW,   RHOSROW,   
     3        SNOROW ,   TCANROW,   RCANROW,   
     4        SCANROW,   GROROW ,   CMAIROW,
     5        TACROW ,   QACROW ,   WSNOROW
C
      REAL    TSFSROW(NLAT,NMOS,4)
C
      REAL,DIMENSION(ILG,IGND) ::
     1        TBARGAT, THLQGAT, THICGAT 
C
      REAL,DIMENSION(ILG) ::
     1        TPNDGAT,   ZPNDGAT,   TBASGAT,   
     2        ALBSGAT,   TSNOGAT,   RHOSGAT,   
     3        SNOGAT ,   TCANGAT,   RCANGAT,   
     4        SCANGAT,   GROGAT ,   CMAIGAT,
     5        TACGAT ,   QACGAT ,   WSNOGAT
C
      REAL    TSFSGAT(ILG,4)
C
C     * GATHER-SCATTER INDEX ARRAYS.
C
      INTEGER  ILMOS (ILG),  JLMOS  (ILG),  IWMOS  (ILG),  JWMOS (ILG),
     1         IWAT  (NLAT), IICE   (NLAT)
C
C     * CANOPY AND SOIL INFORMATION ARRAYS.
C     * (THE LAST DIMENSION OF MOST OF THESE ARRAYS IS GIVEN BY
C     * THE NUMBER OF SOIL LAYERS (IGND), THE NUMBER OF BROAD 
C     * VEGETATION CATEGORIES (ICAN), OR ICAN+1.
C
      REAL,DIMENSION(NLAT,NMOS,ICP1) ::
     1              FCANROW,  LNZ0ROW,
     2              ALVCROW,  ALICROW
C
      REAL,DIMENSION(NLAT,NMOS,ICAN) ::
     1              PAMXROW,  PAMNROW,
     2              CMASROW,  ROOTROW,
     3              RSMNROW,  QA50ROW,
     4              VPDAROW,  VPDBROW,
     5              PSGAROW,  PSGBROW,
     6              PAIDROW,  HGTDROW,
     7              ACVDROW,  ACIDROW
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
     1        THPROW ,  THRROW ,  THMROW , 
     2        BIROW  ,  PSISROW,  GRKSROW,    
     3        THRAROW,  HCPSROW,  
     4        TCSROW ,  THFCROW,  PSIWROW,   
     5        DLZWROW,  ZBTWROW,  THLWROW 
C
      REAL,DIMENSION(NLAT,NMOS) ::
     1        DRNROW ,   XSLPROW,   GRKFROW,
     2        WFSFROW,   WFCIROW,   ALGWROW,   
     3        ALGDROW,   ASVDROW,   ASIDROW,   
     4        AGVDROW,   AGIDROW,   ZSNLROW,
     5        ZPLGROW,   ZPLSROW
C
      REAL,DIMENSION(ILG,IGND) ::
     1        THPGAT ,  THRGAT ,  THMGAT , 
     2        BIGAT  ,  PSISGAT,  GRKSGAT,    
     3        THRAGAT,  HCPSGAT,  
     4        TCSGAT ,  THFCGAT,  PSIWGAT,   
     5        DLZWGAT,  ZBTWGAT,  THLWGAT
C
      REAL,DIMENSION(ILG) ::
     1        DRNGAT ,   XSLPGAT,   GRKFGAT,
     2        WFSFGAT,   WFCIGAT,   ALGWGAT,     
     3        ALGDGAT,   ASVDGAT,   ASIDGAT,     
     4        AGVDGAT,   AGIDGAT,   ZSNLGAT,
     5        ZPLGGAT,   ZPLSGAT
C
      REAL    SANDROW(NLAT,NMOS,IGND), CLAYROW(NLAT,NMOS,IGND), 
     1        ORGMROW(NLAT,NMOS,IGND),
     2        SDEPROW(NLAT,NMOS),      FAREROW(NLAT,NMOS)
C
      INTEGER MIDROW (NLAT,NMOS),     ISNDROW(NLAT,NMOS,IGND), 
     1        ISNDGAT( ILG,IGND),     IORG   (NLAT,NMOS,IGND),
     2        IGDRROW(NLAT,NMOS),     IGDRGAT( ILG)       
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
     1      ZRFMGRD,   ZRFHGRD,   ZDMGRD ,   ZDHGRD ,  
     2      ZBLDGRD,   FSVHGRD,   FSIHGRD,   RADJGRD,
     3      CSZGRD ,   FDLGRD ,   ULGRD  ,   VLGRD  ,   
     4      TAGRD  ,   QAGRD  ,   PRESGRD,   PREGRD ,  
     5      PADRGRD,   VPDGRD ,   TADPGRD,   RHOAGRD,  
     6      RPCPGRD,   TRPCGRD,   SPCPGRD,   TSPCGRD,  
     7      RHSIGRD,   FCLOGRD,   DLONGRD,   UVGRD  ,   
     8      XDIFFUS,   GCGRD  ,   Z0ORGRD,   GGEOGRD,
     9      RPREGRD,   SPREGRD,   VMODGRD
C
      REAL,DIMENSION(ILG) ::
     1      ZRFMGAT,   ZRFHGAT,   ZDMGAT ,   ZDHGAT ,  
     2      ZBLDGAT,   FSVHGAT,   FSIHGAT,   RADJGAT,
     3      CSZGAT ,   FDLGAT ,   ULGAT  ,   VLGAT  ,   
     4      TAGAT  ,   QAGAT  ,   PRESGAT,   PREGAT ,  
     5      PADRGAT,   VPDGAT ,   TADPGAT,   RHOAGAT,  
     6      RPCPGAT,   TRPCGAT,   SPCPGAT,   TSPCGAT,  
     7      RHSIGAT,   FCLOGAT,   DLONGAT,   Z0ORGAT,
     8      GGEOGAT,   VMODGAT
C
C     * LAND SURFACE DIAGNOSTIC VARIABLES.
C
      REAL,DIMENSION(NLAT,NMOS) ::
     1      CDHROW ,   CDMROW ,   HFSROW ,   TFXROW ,  
     2      QEVPROW,   QFSROW ,   QFXROW ,   PETROW ,  
     3      GAROW  ,   EFROW  ,   GTROW  ,   QGROW  ,   
     4      TSFROW ,   ALVSROW,   ALIRROW,   FSNOROW,  
     5      SFCTROW,   SFCUROW,   SFCVROW,   SFCQROW,   
     6      FSGVROW,   FSGSROW,   FSGGROW,   FLGVROW,   
     7      FLGSROW,   FLGGROW,   HFSCROW,   HFSSROW,  
     8      HFSGROW,   HEVCROW,   HEVSROW,   HEVGROW,   
     9      HMFCROW,   HMFNROW,   HTCCROW,   HTCSROW,   
     A      PCFCROW,   PCLCROW,   PCPNROW,   PCPGROW,   
     B      QFGROW ,   QFNROW ,   QFCLROW,   QFCFROW,   
     C      ROFROW ,   ROFOROW,   ROFSROW,   ROFBROW,  
     D      TROFROW,   TROOROW,   TROSROW,   TROBROW,  
     E      ROFCROW,   ROFNROW,   ROVGROW,   WTRCROW,   
     F      WTRSROW,   WTRGROW,   DRROW  ,   WTABROW,  
     G      ILMOROW,   UEROW  ,   HBLROW 
C
      REAL,DIMENSION(ILG) ::
     1      CDHGAT ,   CDMGAT ,   HFSGAT ,   TFXGAT ,  
     2      QEVPGAT,   QFSGAT ,   QFXGAT ,   PETGAT ,  
     3      GAGAT  ,   EFGAT  ,   GTGAT  ,   QGGAT  ,   
     4      TSFGAT ,   ALVSGAT,   ALIRGAT,   FSNOGAT,  
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
     G      ILMOGAT,   UEGAT  ,   HBLGAT ,   SFRHGAT,
     I      FTEMP,     FVAP,      RIB,       QLWOGAT
C
      REAL,DIMENSION(NLAT) ::
     1      CDHGRD ,   CDMGRD ,   HFSGRD ,   TFXGRD ,  
     2      QEVPGRD,   QFSGRD ,   QFXGRD ,   PETGRD ,  
     3      GAGRD  ,   EFGRD  ,   GTGRD  ,   QGGRD  ,   
     4      TSFGRD ,   ALVSGRD,   ALIRGRD,   FSNOGRD,  
     5      SFCTGRD,   SFCUGRD,   SFCVGRD,   SFCQGRD,   
     6      FSGVGRD,   FSGSGRD,   FSGGGRD,   FLGVGRD,   
     7      FLGSGRD,   FLGGGRD,   HFSCGRD,   HFSSGRD,  
     8      HFSGGRD,   HEVCGRD,   HEVSGRD,   HEVGGRD,   
     9      HMFCGRD,   HMFNGRD,   HTCCGRD,   HTCSGRD,   
     A      PCFCGRD,   PCLCGRD,   PCPNGRD,   PCPGGRD,   
     B      QFGGRD ,   QFNGRD ,   QFCLGRD,   QFCFGRD,   
     C      ROFGRD ,   ROFOGRD,   ROFSGRD,   ROFBGRD,  
     D      ROFCGRD,   ROFNGRD,   ROVGGRD,   WTRCGRD,   
     E      WTRSGRD,   WTRGGRD,   DRGRD  ,   WTABGRD,  
     F      ILMOGRD,   UEGRD  ,   HBLGRD
C
      REAL    HMFGROW(NLAT,NMOS,IGND),   HTCROW (NLAT,NMOS,IGND),
     1        QFCROW (NLAT,NMOS,IGND),   GFLXROW(NLAT,NMOS,IGND),
     2        HMFGGAT(ILG,IGND),         HTCGAT (ILG,IGND), 
     3        QFCGAT (ILG,IGND),         GFLXGAT(ILG,IGND),
     4        HMFGGRD(NLAT,IGND),        HTCGRD (NLAT,IGND),
     5        QFCGRD (NLAT,IGND),        GFLXGRD(NLAT,IGND)
C
      INTEGER     ITCTROW(NLAT,NMOS,6,50),  ITCTGAT(ILG,6,50)
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
C     * MONTHLY OUTPUT FOR CLASS GRID-MEAN
C 
      REAL,DIMENSION(NLAT) ::
     1              ALVSACC_MO,ALIRACC_MO,FLUTACC_MO, 
     2              FSINACC_MO,FLINACC_MO,HFSACC_MO,
     3              QEVPACC_MO,SNOACC_MO, WSNOACC_MO,
     4              ROFACC_MO, PREACC_MO, EVAPACC_MO,
     5              TAACC_MO
      REAL ::       FSSTAR_MO,FLSTAR_MO,QH_MO,QE_MO

      REAL TBARACC_MO(NLAT,IGND), THLQACC_MO(NLAT,IGND),
     1     THICACC_MO(NLAT,IGND) 
C
C     * YEARLY OUTPUT FOR CLASS GRID-MEAN
C
      REAL,DIMENSION(NLAT) ::
     1              ALVSACC_YR,ALIRACC_YR,FLUTACC_YR, 
     2              FSINACC_YR,FLINACC_YR,HFSACC_YR,
     3              QEVPACC_YR,ROFACC_YR, PREACC_YR, 
     4              EVAPACC_YR,TAACC_YR
      REAL ::       FSSTAR_YR,FLSTAR_YR,QH_YR,QE_YR
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
     C      TRSNOW (ILG), CHCAP  (ILG), CHCAPS (ILG),
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
C     * DIAGNOSTIC ARRAYS USED FOR CHECKING ENERGY AND WATER 
C     * BALANCES.
C
      REAL CTVSTP(ILG),   CTSSTP(ILG),   CT1STP(ILG),   CT2STP(ILG),
     1     CT3STP(ILG),   WTVSTP(ILG),   WTSSTP(ILG),   WTGSTP(ILG)
C
C     * CONSTANTS AND TEMPORARY VARIABLES.
C
      REAL DEGLAT,DEGLON,FSDOWN,DAY,DECL,HOUR,COSZ,CUMSNO,
     1     QSUMV,QSUMS,QSUM1,QSUM2,QSUM3,WSUMV,WSUMS,WSUMG,ALTOT,
     2     FSSTAR,FLSTAR,QH,QE,BEG,SNOMLT,ZSN,TCN,TSN,TPN,GTOUT,
     3     ALTOT_MO,ALTOT_YR
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
     7     ANGMAX
C
c================= CTEM array declaration ===============================\
c
c     variables for coupling CLASS and CTEM
c
      integer argcount,iargc,strlen,ictemmod
c
      character*80   titlec1, titlec2, titlec3
      character*80   argbuff
      character*160  command
c
      logical ctem_on,    parallelrun,    mosaic,
     1        cyclemet,   dofire,         run_model,
     2     met_rewound,   reach_eof,      compete, 
     3      start_bare,   rsfile,         lnduseon, 
     4           co2on,   popdon,         inibioclim,
     5    start_from_rs, dowetlands,      obswetf,
     6    transient_run
c
       integer   lopcount,  isumc,   nol2pfts(4),  
     1           k1c,       k2c,     iyd,         jhhstd,
     2           jhhendd,   jdstd,   jdendd,      jhhsty,
     3           jhhendy,   jdsty,   jdendy,      jhhst, 
     4           jhhend,    jdst,    jdend,       ctemloop,  
     5           spinfast,  month1,  month2,      xday, 
     6           ncyear, co2yr, popyr, nummetcylyrs,
     7           metcylyrst, metcycendyr, climiyear, popcycleyr,
     8           cypopyr, lucyr, cylucyr, endyr,bigpftc(2),
     9           obswetyr, cywetldyr, trans_startyr  
c
       real      rlim,        fsstar_g,
     1           flstar_g,  qh_g,    qe_g,        snomlt_g,
     2           beg_g,     gtout_g, tpn_g,       altot_g,
     3           tcn_g,     tsn_g,   zsn_g       

       real      co2concin,  popdin,    setco2conc, sumfare,
     1           temp_var, barefrac,  todfrac(ilg,icc), barf(nlat)

      real grclarea(ilg), crop_temp_frac(ilg,2)
c
      real tcanrs(nlat,nmos), tsnors(nlat,nmos), tpndrs(nlat,nmos),
     1     csum(nlat,nmos,ican),       tbaraccrow_m(nlat,nmos,ignd),
     2     tcanoaccrow_m(nlat,nmos),   lightng(ilg),
     3     uvaccrow_m(nlat,nmos),      vvaccrow_m(nlat,nmos)
        
c     wilting and field capacities vars
       ! FLAG these can be removed soon - as long as the PSIWLT limit gets decided between CLASS and CTEM - JM Jan 14 2015
c      real     fieldsm(ilg,ignd),     wiltsm(ilg,ignd)
c      real     psisat(ilg,ignd),      grksat(ilg,ignd)
c      real     thpor(ilg,ignd),       bterm(ilg,ignd)

c     Competition related variables

       real fsinacc_gat(ilg), flutacc_gat(ilg), flinacc_gat(ilg),
     1      alswacc_gat(ilg), allwacc_gat(ilg), pregacc_gat(ilg),
     2      altot_gat,        fsstar_gat,       flstar_gat,
     3      netrad_gat(ilg),  preacc_gat(ilg)
c
       real tcurm(ilg),       srpcuryr   (ilg), dftcuryr(ilg),
     1      tmonth(12,ilg),      anpcpcur(ilg),  anpecur(ilg),
     2      gdd5cur(ilg),        surmncur(ilg), defmncur(ilg),
     3      srplscur(ilg),       defctcur(ilg)
  
       real, dimension(ilg) :: twarmm    ! temperature of the warmest month (c)
       real, dimension(ilg) :: tcoldm    ! temperature of the coldest month (c)
       real, dimension(ilg) :: gdd5      ! growing degree days above 5 c
       real, dimension(ilg) :: aridity   ! aridity index, ratio of potential evaporation to precipitation
       real, dimension(ilg) :: srplsmon  ! number of months in a year with surplus water i.e.
                                                  !  precipitation more than potential evaporation
       real, dimension(ilg) :: defctmon  ! number of months in a year with water deficit i.e.
                                                  ! precipitation less than potential evaporation
       real, dimension(ilg) :: anndefct  ! annual water deficit (mm) 
       real, dimension(ilg) :: annsrpls  ! annual water surplus (mm)
       real, dimension(ilg) :: annpcp    ! annual precipitation (mm)
       real, dimension(ilg) :: dry_season_length  ! length of dry season (months)
c
       real lyglfmasgat(ilg,icc),   geremortgat(ilg,icc),
     1      intrmortgat(ilg,icc),     lambdagat(ilg,icc),
     2            rnded_pft(icc),
     3            ccgat(ilg,icc),         mmgat(ilg,icc),
     4            temparray(icc),                   temp

      real  xdiffusgat(ilg) ! the corresponding ROW is CLASS's XDIFFUS    
     
!     For these below, the corresponding ROWs are defined by CLASS  
 
      real  sdepgat(ilg),       orgmgat(ilg,ignd), 
     1      sandgat(ilg,ignd),  claygat(ilg,ignd)

!     Set up the variables that have both row and gat such that
!     they are beside each other in the declaration. This makes it 
!     more apparent if the naming is inconsistent or if there is a problem
!     with how the variable is defined. (JM Jul 2013)

!     ||      This column is ROW       ||  This column is GAT    ||

      logical pftexistrow(nlat,nmos,icc), pftexistgat(ilg,icc)

      integer colddaysrow(nlat,nmos,2), colddaysgat(ilg,2),
     1     icountrow(nlat,nmos),        icount(ilg),
     2     lfstatusrow(nlat,nmos,icc),  lfstatusgat(ilg,icc),
     3     pandaysrow(nlat,nmos,icc),   pandaysgat(ilg,icc),
     4     stdalngrd(nlat),             stdalngat(ilg)
c
      real ailcminrow(nlat,nmos,icc),   ailcmingat(ilg,icc),
     1     ailcmaxrow(nlat,nmos,icc),   ailcmaxgat(ilg,icc),
     2     dvdfcanrow(nlat,nmos,icc),   dvdfcan(ilg,icc),
     3     gleafmasrow(nlat,nmos,icc),  gleafmasgat(ilg,icc),
     4     bleafmasrow(nlat,nmos,icc),  bleafmasgat(ilg,icc),
     5     stemmassrow(nlat,nmos,icc),  stemmassgat(ilg,icc),
     6     rootmassrow(nlat,nmos,icc),  rootmassgat(ilg,icc),
     7     pstemmassrow(nlat,nmos,icc), pstemmassgat(ilg,icc),
     8     pgleafmassrow(nlat,nmos,icc), pgleafmassgat(ilg,icc)
c
      real fcancmxrow(nlat,nmos,icc),  fcancmxgat(ilg,icc),
     1     gavglairow(nlat,nmos),      gavglaigat(ilg),
     2     zolncrow(nlat,nmos,ican),   zolncgat(ilg,ican),
     3     ailcrow(nlat,nmos,ican),    ailcgat(ilg,ican),
     4     ailcgrow(nlat,nmos,icc),    ailcggat(ilg,icc),
     5     ailcgsrow(nlat,nmos,icc),   ailcgsgat(ilg,icc),
     6     fcancsrow(nlat,nmos,icc),   fcancsgat(ilg,icc),
     7     fcancrow(nlat,nmos,icc),    fcancgat(ilg,icc),
     8     co2concrow(nlat,nmos),      co2concgat(ilg),
     9     co2i1cgrow(nlat,nmos,icc),  co2i1cggat(ilg,icc),
     a     co2i1csrow(nlat,nmos,icc),  co2i1csgat(ilg,icc), 
     b     co2i2cgrow(nlat,nmos,icc),  co2i2cggat(ilg,icc),
     c     co2i2csrow(nlat,nmos,icc),  co2i2csgat(ilg,icc),
     d     ancsvegrow(nlat,nmos,icc),  ancsveggat(ilg,icc),
     e     ancgvegrow(nlat,nmos,icc),  ancgveggat(ilg,icc),
     f     rmlcsvegrow(nlat,nmos,icc), rmlcsveggat(ilg,icc),
     g     rmlcgvegrow(nlat,nmos,icc), rmlcgveggat(ilg,icc),
     h     slairow(nlat,nmos,icc),     slaigat(ilg,icc),
     i     ailcbrow(nlat,nmos,icc),    ailcbgat(ilg,icc),
     j     canresrow(nlat,nmos),       canresgat(ilg),
     a     flhrlossrow(nlat,nmos,icc), flhrlossgat(ilg,icc)
c
      real grwtheffrow(nlat,nmos,icc), grwtheffgat(ilg,icc),
     1     lystmmasrow(nlat,nmos,icc), lystmmasgat(ilg,icc),
     2     lyrotmasrow(nlat,nmos,icc), lyrotmasgat(ilg,icc),
     3     tymaxlairow(nlat,nmos,icc), tymaxlaigat(ilg,icc),
     4     vgbiomasrow(nlat,nmos),     vgbiomasgat(ilg),
     5     gavgltmsrow(nlat,nmos),     gavgltmsgat(ilg), 
     6     gavgscmsrow(nlat,nmos),     gavgscmsgat(ilg),
     7     stmhrlosrow(nlat,nmos,icc), stmhrlosgat(ilg,icc)

      real rmatcrow(nlat,nmos,ican,ignd),  rmatcgat(ilg,ican,ignd),
     1     rmatctemrow(nlat,nmos,icc,ignd),rmatctemgat(ilg,icc,ignd),
     2     litrmassrow(nlat,nmos,iccp1),   litrmassgat(ilg,iccp1),     
     3     soilcmasrow(nlat,nmos,iccp1),   soilcmasgat(ilg,iccp1),
     4     vgbiomas_vegrow(nlat,nmos,icc), vgbiomas_veggat(ilg,icc)

c     Fire-related variables

      real emit_co2row(nlat,nmos,icc), emit_co2gat(ilg,icc),
     1     emit_corow(nlat,nmos,icc),  emit_cogat(ilg,icc),
     2     emit_ch4row(nlat,nmos,icc), emit_ch4gat(ilg,icc),
     3     emit_nmhcrow(nlat,nmos,icc),emit_nmhcgat(ilg,icc),
     4     emit_h2row(nlat,nmos,icc),  emit_h2gat(ilg,icc), 
     5     emit_noxrow(nlat,nmos,icc), emit_noxgat(ilg,icc),
     6     emit_n2orow(nlat,nmos,icc), emit_n2ogat(ilg,icc),
     7     emit_pm25row(nlat,nmos,icc),emit_pm25gat(ilg,icc),
     8     emit_tpmrow(nlat,nmos,icc), emit_tpmgat(ilg,icc),
     9     emit_tcrow(nlat,nmos,icc),  emit_tcgat(ilg,icc),
     a     emit_ocrow(nlat,nmos,icc),  emit_ocgat(ilg,icc),
     b     emit_bcrow(nlat,nmos,icc),  emit_bcgat(ilg,icc),
     c     burnfracrow(nlat,nmos),     burnfracgat(ilg),
     d     burnvegfrow(nlat,nmos,icc), burnvegfgat(ilg,icc),
     e     probfirerow(nlat,nmos),     probfiregat(ilg),
     f     btermrow(nlat,nmos),        btermgat(ilg),
     g     ltermrow(nlat,nmos),        ltermgat(ilg),
     h     mtermrow(nlat,nmos),        mtermgat(ilg)

       real extnprobgrd(nlat),          extnprobgat(ilg),
     1      prbfrhucgrd(nlat),          prbfrhucgat(ilg), 
     1      mlightnggrd(nlat,12),       mlightnggat(ilg,12)

c      Methane(wetland) related variables    !Rudra added on 03/12/2013

       real  WETFRACGRD(nlat),              wetfrac_sgrd(ilg,8), 
!     1       WETFRAC_SROW(nlat),            WETFRAC_SGAT(ILG),
     1       CH4WET1ROW(nlat,nmos),         CH4WET1GAT(ILG),
     2       CH4WET2ROW(nlat,nmos),         CH4WET2GAT(ILG),
     3       WETFDYNROW(nlat,nmos),         WETFDYNGAT(ILG),
     4       CH4DYN1ROW(nlat,nmos),         CH4DYN1GAT(ILG),
     5       CH4DYN2ROW(nlat,nmos),         CH4DYN2GAT(ILG),
     6       wetfrac_mon(nlat,12)

!      Land-use related variables

      real lucemcomrow(nlat,nmos),     lucemcomgat(ilg), 
     1     lucltrinrow(nlat,nmos),     lucltringat(ilg),
     2     lucsocinrow(nlat,nmos),     lucsocingat(ilg)
c
      real bmasvegrow(nlat,nmos,icc),  bmasveggat(ilg,icc),
     1     cmasvegcrow(nlat,nmos,ican),cmasvegcgat(ilg,ican),
     2     veghghtrow(nlat,nmos,icc),  veghghtgat(ilg,icc),
     3     rootdpthrow(nlat,nmos,icc), rootdpthgat(ilg,icc),
     4     rmlrow(nlat,nmos),          rmlgat(ilg),
     5     rmsrow(nlat,nmos),          rmsgat(ilg),
     6     tltrleafrow(nlat,nmos,icc), tltrleafgat(ilg,icc),
     7     tltrstemrow(nlat,nmos,icc), tltrstemgat(ilg,icc),
     8     tltrrootrow(nlat,nmos,icc), tltrrootgat(ilg,icc), 
     9     leaflitrrow(nlat,nmos,icc), leaflitrgat(ilg,icc),
     a     roottemprow(nlat,nmos,icc), roottempgat(ilg,icc),
     b     afrleafrow(nlat,nmos,icc),  afrleafgat(ilg,icc),
     c     afrstemrow(nlat,nmos,icc),  afrstemgat(ilg,icc),
     d     afrrootrow(nlat,nmos,icc),  afrrootgat(ilg,icc),
     e     wtstatusrow(nlat,nmos,icc), wtstatusgat(ilg,icc),
     f     ltstatusrow(nlat,nmos,icc), ltstatusgat(ilg,icc),
     g     rmrrow(nlat,nmos),          rmrgat(ilg)

      real npprow(nlat,nmos),          nppgat(ilg),
     1     neprow(nlat,nmos),          nepgat(ilg),
     2     nbprow(nlat,nmos),          nbpgat(ilg),
     3     gpprow(nlat,nmos),          gppgat(ilg),
     4     hetroresrow(nlat,nmos),     hetroresgat(ilg),
     5     autoresrow(nlat,nmos),      autoresgat(ilg),
     6     soilcresprow(nlat,nmos),    soilcrespgat(ilg),
     7     rmrow(nlat,nmos),           rmgat(ilg),
     8     rgrow(nlat,nmos),           rggat(ilg),
     9     litresrow(nlat,nmos),       litresgat(ilg),
     a     socresrow(nlat,nmos),       socresgat(ilg),
     b     dstcemlsrow(nlat,nmos),     dstcemlsgat(ilg),
     c     litrfallrow(nlat,nmos),     litrfallgat(ilg),
     d     humiftrsrow(nlat,nmos),     humiftrsgat(ilg)

      real  gppvegrow(nlat,nmos,icc),   gppveggat(ilg,icc),
     1      nepvegrow(nlat,nmos,iccp1),   nepveggat(ilg,iccp1),
     2      nbpvegrow(nlat,nmos,iccp1),   nbpveggat(ilg,iccp1),
     3      nppvegrow(nlat,nmos,icc),   nppveggat(ilg,icc), 
     4      hetroresvegrow(nlat,nmos,iccp1),hetroresveggat(ilg,iccp1),
     5      autoresvegrow(nlat,nmos,icc),autoresveggat(ilg,icc),
     6      litresvegrow(nlat,nmos,iccp1),litresveggat(ilg,iccp1),
     7      soilcresvegrow(nlat,nmos,iccp1),soilcresveggat(ilg,iccp1),
     8      rmlvegaccrow(nlat,nmos,icc),rmlvegaccgat(ilg,icc),     
     9      rmsvegrow(nlat,nmos,icc),   rmsveggat(ilg,icc),
     a      rmrvegrow(nlat,nmos,icc),   rmrveggat(ilg,icc),        
     b      rgvegrow(nlat,nmos,icc),    rgveggat(ilg,icc) 
c
      real 
     1     rothrlosrow(nlat,nmos,icc), rothrlosgat(ilg,icc),
     2     pfcancmxrow(nlat,nmos,icc), pfcancmxgat(ilg,icc),
     3     nfcancmxrow(nlat,nmos,icc), nfcancmxgat(ilg,icc),
     4     alvsctmrow(nlat,nmos,ican), alvsctmgat(ilg,ican),
     5     paicrow(nlat,nmos,ican),    paicgat(ilg,ican),
     6     slaicrow(nlat,nmos,ican),   slaicgat(ilg,ican),
     7     alirctmrow(nlat,nmos,ican), alirctmgat(ilg,ican),
     8     cfluxcgrow(nlat,nmos),      cfluxcggat(ilg),
     9     cfluxcsrow(nlat,nmos),      cfluxcsgat(ilg),
     a     dstcemls3row(nlat,nmos),    dstcemls3gat(ilg),
     b     anvegrow(nlat,nmos,icc),    anveggat(ilg,icc),                           
     c     rmlvegrow(nlat,nmos,icc),   rmlveggat(ilg,icc)
 
!      Outputs

       real tcanoaccrow_out(nlat,nmos),tcanoaccgat_out(ilg),
     1     qevpacc_m_save(nlat,nmos)
c
       integer ifcancmx_g(nlat,icc)

!     -----------------------
!      Mosaic-level variables (denoted by an ending of "_m")  
c
       real   PREACC_M(NLAT,NMOS),     GTACC_M(NLAT,NMOS),
     1     QEVPACC_M(NLAT,NMOS),       HFSACC_M(NLAT,NMOS),
     2     HMFNACC_M(NLAT,NMOS),       ROFACC_M(NLAT,NMOS),
     3     SNOACC_M(NLAT,NMOS),        OVRACC_M(NLAT,NMOS),
     3     WTBLACC_M(NLAT,NMOS),       TBARACC_M(NLAT,NMOS,IGND),
     4     THLQACC_M(NLAT,NMOS,IGND),  THICACC_M(NLAT,NMOS,IGND),
     5     THALACC_M(NLAT,NMOS,IGND),  ALVSACC_M(NLAT,NMOS),
     6     ALIRACC_M(NLAT,NMOS),       RHOSACC_M(NLAT,NMOS),
     7     TSNOACC_M(NLAT,NMOS),       WSNOACC_M(NLAT,NMOS),
     8     TCANACC_M(NLAT,NMOS),       RCANACC_M(NLAT,NMOS),
     9     SCANACC_M(NLAT,NMOS),       GROACC_M(NLAT,NMOS),
     A     FSINACC_M(NLAT,NMOS),       FLINACC_M(NLAT,NMOS),
     B     TAACC_M(NLAT,NMOS),         UVACC_M(NLAT,NMOS),
     C     PRESACC_M(NLAT,NMOS),       QAACC_M(NLAT,NMOS),
     D     EVAPACC_M(NLAT,NMOS),       FLUTACC_M(NLAT,NMOS)

      real  fsnowacc_m(ilg),           tcansacc_m(ilg),
     1      tbarcacc_m(ilg,ignd),      tbarcsacc_m(ilg,ignd),
     2      tbargacc_m(ilg,ignd),      tbargsacc_m(ilg,ignd),
     3      thliqcacc_m(ilg,ignd),     thliqgacc_m(ilg,ignd),
     4      thicecacc_m(ilg,ignd),     faregat(ilg), 
     5      tcanoaccgat_m(ilg),        taaccgat_m(ilg),
     6      uvaccgat_m(ilg),           vvaccgat_m(ilg),
     7      tbaraccgat_m(ilg,ignd)            

      real ancsvgac_m(ilg,icc),        ancgvgac_m(ilg,icc), 
     b     rmlcsvga_m(ilg,icc),        rmlcgvga_m(ilg,icc)

       integer ifcancmx_m(nlat,nmos)

       real leaflitr_m(nlat,nmos),     tltrleaf_m(nlat,nmos),
     1     tltrstem_m(nlat,nmos),      tltrroot_m(nlat,nmos),
     2     ailcg_m(nlat,nmos),         ailcb_m(nlat,nmos),
     3     rmatctem_m(nlat,nmos,ignd), veghght_m(nlat,nmos),
     4     rootdpth_m(nlat,nmos),      roottemp_m(nlat,nmos),
     5     slai_m(nlat,nmos),          afrroot_m(nlat,nmos), 
     6     afrleaf_m(nlat,nmos),       afrstem_m(nlat,nmos),
     7     laimaxg_m(nlat,nmos),       stemmass_m(nlat,nmos),
     8     rootmass_m(nlat,nmos),      litrmass_m(nlat,nmos),
     9     gleafmas_m(nlat,nmos),      bleafmas_m(nlat,nmos),
     a     soilcmas_m(nlat,nmos)
c
!     -----------------------
!     Grid-averaged variables (denoted with an ending of "_g")

      real  wsnorow_g(nlat),           rofsrow_g(nlat),  
     1     snorow_g(nlat),             rhosrow_g(nlat),
     2     rofrow_g(nlat),             zpndrow_g(nlat),
     3     rcanrow_g(nlat),            scanrow_g(nlat), 
     4     trofrow_g(nlat),            troorow_g(nlat), 
     5     trobrow_g(nlat),            roforow_g(nlat),
     6     rofbrow_g(nlat),            trosrow_g(nlat),
     7     fsgvrow_g(nlat),            fsgsrow_g(nlat),
     8     flgvrow_g(nlat),            flgsrow_g(nlat),
     9     hfscrow_g(nlat),            hfssrow_g(nlat),
     a     hevcrow_g(nlat),            hevsrow_g(nlat),
     b     hmfcrow_g(nlat),            hmfnrow_g(nlat), 
     c     htcsrow_g(nlat),            htccrow_g(nlat),
     d     fsggrow_g(nlat),            flggrow_g(nlat),
     e     hfsgrow_g(nlat),            hevgrow_g(nlat),
     f     fc_g(nlat),                 fg_g(nlat),
     g     fcs_g(nlat),                fgs_g(nlat),
     h     pcfcrow_g(nlat),            pclcrow_g(nlat),
     1     pcpgrow_g(nlat),            qfcfrow_g(nlat),
     2     qfgrow_g(nlat),             qfcrow_g(nlat,ignd),
     3     rofcrow_g(nlat),            rofnrow_g(nlat),
     4     wtrsrow_g(nlat),            wtrgrow_g(nlat),
     5     pcpnrow_g(nlat),            qfclrow_g(nlat),
     6     qfnrow_g(nlat),             wtrcrow_g(nlat),
     8     gpp_g(nlat),                npp_g(nlat),     
     9     nbp_g(nlat),                autores_g(nlat), 
     a     litres_g(nlat),             socres_g(nlat),  
     b     dstcemls3_g(nlat),          litrfall_g(nlat), 
     c     rml_g(nlat),                rms_g(nlat),         
     d     rg_g(nlat),                 leaflitr_g(nlat), 
     e     tltrstem_g(nlat),           tltrroot_g(nlat),
     f     nep_g(nlat),                hetrores_g(nlat),
     g     dstcemls_g(nlat),           humiftrs_g(nlat),
     h     rmr_g(nlat),                tltrleaf_g(nlat),
     i     gavgltms_g(nlat)
c
       real  hmfgrow_g(nlat,ignd),     htcrow_g(nlat,ignd), 
     1     tbarrow_g(nlat,ignd),       thlqrow_g(nlat,ignd),
     2     thicrow_g(nlat,ignd),       gflxrow_g(nlat,ignd),
     3     anvegrow_g(nlat,icc),       rmlvegrow_g(nlat,icc)
c
       real vgbiomas_g(nlat),          gavglai_g(nlat),
     1     gavgscms_g(nlat),           gleafmas_g(nlat),
     2     bleafmas_g(nlat),           stemmass_g(nlat),
     3     rootmass_g(nlat),           litrmass_g(nlat),
     4     soilcmas_g(nlat),           slai_g(nlat),
     5     ailcg_g(nlat),              ailcb_g(nlat),
     6     rmatctem_g(nlat,ignd),      veghght_g(nlat),
     7     rootdpth_g(nlat),           roottemp_g(nlat),
     a     totcmass_g(nlat)
c
       real afrleaf_g(nlat,icc),       afrstem_g(nlat,icc),
     1     afrroot_g(nlat,icc),        lfstatus_g(nlat,icc),
     2     tcanoacc_out_g(nlat),      
     3     burnfrac_g(nlat),           probfire_g(nlat),
     4     lucemcom_g(nlat),           lucltrin_g(nlat),
     5     lucsocin_g(nlat),
     6     emit_co2_g(nlat),           emit_co_g(nlat),
     7     emit_ch4_g(nlat),           emit_nmhc_g(nlat),
     8     emit_h2_g(nlat),            emit_nox_g(nlat),
     9     emit_n2o_g(nlat),           emit_pm25_g(nlat),
     a     emit_tpm_g(nlat),           emit_tc_g(nlat),
     b     emit_oc_g(nlat),            emit_bc_g(nlat),
     c     bterm_g(nlat),              lterm_g(nlat),
     d     mterm_g(nlat)

       real    CH4WET1_G(nlat),            CH4WET2_G(nlat),        !Rudra addes on 03/12/2013
     &     WETFDYN_G(nlat),            CH4DYN1_G(nlat),
     &     CH4DYN2_G(nlat)  

!     -----------------------
!      Grid averaged monthly variables (denoted by name ending in "_mo_g")

       real laimaxg_mo_g(nlat),          stemmass_mo_g(nlat),
     1      rootmass_mo_g(nlat),         litrmass_mo_g(nlat),
     2      soilcmas_mo_g(nlat),         npp_mo_g(nlat),
     3      gpp_mo_g(nlat),              nep_mo_g(nlat),
     4      nbp_mo_g(nlat),              hetrores_mo_g(nlat),
     5      autores_mo_g(nlat),          litres_mo_g(nlat),
     6      soilcres_mo_g(nlat),
     7      vgbiomas_mo_g(nlat),         totcmass_mo_g(nlat),
!
     8      emit_co2_mo_g(nlat),       emit_co_mo_g(nlat),
     9      emit_ch4_mo_g(nlat),       emit_nmhc_mo_g(nlat),
     a      emit_h2_mo_g(nlat),        emit_nox_mo_g(nlat),
     b      emit_n2o_mo_g(nlat),       emit_pm25_mo_g(nlat),
     c      emit_tpm_mo_g(nlat),       emit_tc_mo_g(nlat), 
     d      emit_oc_mo_g(nlat),        emit_bc_mo_g(nlat),
     e      probfire_mo_g(nlat),       luc_emc_mo_g(nlat),
     f      lucltrin_mo_g(nlat),       lucsocin_mo_g(nlat),
     g      burnfrac_mo_g(nlat),       bterm_mo_g(nlat),
     h      lterm_mo_g(nlat),          mterm_mo_g(nlat),
!    
     &      ch4wet1_mo_g(nlat),        ch4wet2_mo_g(nlat),     !Rudra added on 03/12/2013
     &      wetfdyn_mo_g(nlat),        ch4dyn1_mo_g(nlat),
     &      ch4dyn2_mo_g(nlat)  

!      Mosaic monthly variables (denoted by name ending in "_mo_m")

       real laimaxg_mo_m(nlat,nmos,icc),   stemmass_mo_m(nlat,nmos,icc),
     1      rootmass_mo_m(nlat,nmos,icc),litrmass_mo_m(nlat,nmos,iccp1),
     2      soilcmas_mo_m(nlat,nmos,iccp1),  npp_mo_m(nlat,nmos,icc),
     3      gpp_mo_m(nlat,nmos,icc),       nep_mo_m(nlat,nmos,iccp1), 
     4      nbp_mo_m(nlat,nmos,iccp1),     vgbiomas_mo_m(nlat,nmos,icc),  
     a      hetrores_mo_m(nlat,nmos,iccp1), autores_mo_m(nlat,nmos,icc), 
     b      litres_mo_m(nlat,nmos,iccp1),soilcres_mo_m(nlat,nmos,iccp1),
     c      totcmass_mo_m(nlat,nmos,icc)
c
       real emit_co2_mo_m(nlat,nmos,icc),emit_co_mo_m(nlat,nmos,icc),
     1      emit_ch4_mo_m(nlat,nmos,icc),emit_nmhc_mo_m(nlat,nmos,icc),
     2      emit_h2_mo_m(nlat,nmos,icc), emit_nox_mo_m(nlat,nmos,icc),
     3      emit_n2o_mo_m(nlat,nmos,icc),emit_pm25_mo_m(nlat,nmos,icc),
     4      emit_tpm_mo_m(nlat,nmos,icc),emit_tc_mo_m(nlat,nmos,icc),
     5      emit_oc_mo_m(nlat,nmos,icc), emit_bc_mo_m(nlat,nmos,icc),
     6      probfire_mo_m(nlat,nmos),    bterm_mo_m(nlat,nmos),
     7      luc_emc_mo_m(nlat,nmos),     lterm_mo_m(nlat,nmos),
     8      lucsocin_mo_m(nlat,nmos),    mterm_mo_m(nlat,nmos),
     9      lucltrin_mo_m(nlat,nmos),
     a      burnfrac_mo_m(nlat,nmos,icc),
c                                                            !Rudra CH4 variable
     &      ch4wet1_mo_m(nlat,nmos), ch4wet2_mo_m(nlat,nmos),
     &      wetfdyn_mo_m(nlat,nmos), 
     &      ch4dyn1_mo_m(nlat,nmos), ch4dyn2_mo_m(nlat,nmos) 
c
!     -----------------------
c      Annual output for CTEM grid-averaged variables:
c      (denoted by name ending in "_yr_g")

       real laimaxg_yr_g(nlat),        stemmass_yr_g(nlat),            
     1      rootmass_yr_g(nlat),       litrmass_yr_g(nlat),
     2      soilcmas_yr_g(nlat),       npp_yr_g(nlat),   
     3      gpp_yr_g(nlat),            nep_yr_g(nlat), 
     a      nbp_yr_g(nlat),            hetrores_yr_g(nlat),
     b      autores_yr_g(nlat),        litres_yr_g(nlat),
     c      soilcres_yr_g(nlat),
     4      vgbiomas_yr_g(nlat),       totcmass_yr_g(nlat),
     5      emit_co2_yr_g(nlat),       emit_co_yr_g(nlat),
     6      emit_ch4_yr_g(nlat),       emit_nmhc_yr_g(nlat),
     7      emit_h2_yr_g(nlat),        emit_nox_yr_g(nlat),
     8      emit_n2o_yr_g(nlat),       emit_pm25_yr_g(nlat),
     9      emit_tpm_yr_g(nlat),       emit_tc_yr_g(nlat),
     a      emit_oc_yr_g(nlat),        emit_bc_yr_g(nlat),
     b      probfire_yr_g(nlat),       luc_emc_yr_g(nlat),
     c      lucltrin_yr_g(nlat),       lucsocin_yr_g(nlat),
     d      burnfrac_yr_g(nlat),       bterm_yr_g(nlat),
     e      lterm_yr_g(nlat),          mterm_yr_g(nlat),
     &      ch4wet1_yr_g(nlat),        ch4wet2_yr_g(nlat),
     &      wetfdyn_yr_g(nlat),        ch4dyn1_yr_g(nlat),
     &      ch4dyn2_yr_g(nlat)

c      Annual output for CTEM mosaic variables:
c      (denoted by name ending in "_yr_m")

       real laimaxg_yr_m(nlat,nmos,icc), stemmass_yr_m(nlat,nmos,icc),
     1      rootmass_yr_m(nlat,nmos,icc),litrmass_yr_m(nlat,nmos,iccp1),
     2      soilcmas_yr_m(nlat,nmos,iccp1),npp_yr_m(nlat,nmos,icc),
     3      gpp_yr_m(nlat,nmos,icc),       nep_yr_m(nlat,nmos,iccp1),
     4      nbp_yr_m(nlat,nmos,iccp1),     vgbiomas_yr_m(nlat,nmos,icc),
     a      hetrores_yr_m(nlat,nmos,iccp1),autores_yr_m(nlat,nmos,icc),
     b      litres_yr_m(nlat,nmos,iccp1),soilcres_yr_m(nlat,nmos,iccp1),  
     5      emit_co2_yr_m(nlat,nmos,icc), emit_co_yr_m(nlat,nmos,icc),
     6      emit_ch4_yr_m(nlat,nmos,icc), emit_nmhc_yr_m(nlat,nmos,icc),
     7      emit_h2_yr_m(nlat,nmos,icc), emit_nox_yr_m(nlat,nmos,icc),
     8      emit_n2o_yr_m(nlat,nmos,icc),emit_pm25_yr_m(nlat,nmos,icc),
     9      emit_tpm_yr_m(nlat,nmos,icc),emit_tc_yr_m(nlat,nmos,icc),
     a      emit_oc_yr_m(nlat,nmos,icc), emit_bc_yr_m(nlat,nmos,icc),
     b      probfire_yr_m(nlat,nmos),    bterm_yr_m(nlat,nmos),
     c      luc_emc_yr_m(nlat,nmos),     totcmass_yr_m(nlat,nmos,icc),
     d      lucsocin_yr_m(nlat,nmos),    lterm_yr_m(nlat,nmos),   
     e      lucltrin_yr_m(nlat,nmos),    mterm_yr_m(nlat,nmos),
     f      burnfrac_yr_m(nlat,nmos,icc),
     &      ch4wet1_yr_m(nlat,nmos),   ch4wet2_yr_m(nlat,nmos),
     &      wetfdyn_yr_m(nlat,nmos),   ch4dyn1_yr_m(nlat,nmos),
     &      ch4dyn2_yr_m(nlat,nmos)
  
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
C
      CALL CLASSD
C
      ZDMGRD(1)=10.0
      ZDHGRD(1)=2.0
      CUMSNO=0.0
C
C===================== CTEM ==============================================\

c     all model switches are read in from a namelist file
      call read_from_job_options(argbuff,mosaic,transient_run,
     1             trans_startyr,ctemloop,ctem_on,ncyear,lnduseon,
     2             spinfast,cyclemet,nummetcylyrs,metcylyrst,co2on,
     3             setco2conc,popdon,popcycleyr,parallelrun,dofire,
     4             dowetlands,obswetf,compete,inibioclim,start_bare,
     5             rsfile,start_from_rs,idisp,izref,islfd,ipcp,itc,
     6             itcg,itg,iwf,ipai,ihgt,ialc,ials,ialg,jhhstd,
     7             jhhendd,jdstd,jdendd,jhhsty,jhhendy,jdsty,jdendy)

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
      IMONTH = 0   
      DO I=1,NLAT
        DO M=1,NMOS
          PREACC_M(I,M)=0.
          GTACC_M(I,M)=0.
          QEVPACC_M(I,M)=0.
          HFSACC_M(I,M)=0.
          HMFNACC_M(I,M)=0.
          ROFACC_M(I,M)=0.
          SNOACC_M(I,M)=0.
C         CANARE(I,M)=0.
C         SNOARE(I,M)=0.
          OVRACC_M(I,M)=0.
          WTBLACC_M(I,M)=0.
              DO J=1,IGND
                TBARACC_M(I,M,J)=0.
                THLQACC_M(I,M,J)=0.
                THICACC_M(I,M,J)=0.
                THALACC_M(I,M,J)=0.
              ENDDO
          ALVSACC_M(I,M)=0.
          ALIRACC_M(I,M)=0.
          RHOSACC_M(I,M)=0.
          TSNOACC_M(I,M)=0.
          WSNOACC_M(I,M)=0.
          TCANACC_M(I,M)=0.
          RCANACC_M(I,M)=0.
          SCANACC_M(I,M)=0.
          GROACC_M(I,M)=0.
          FSINACC_M(I,M)=0.
          FLINACC_M(I,M)=0.
          TAACC_M(I,M)=0.
          UVACC_M(I,M)=0.
          PRESACC_M(I,M)=0.
          QAACC_M(I,M)=0.
          EVAPACC_M(I,M)=0.
          FLUTACC_M(I,M)=0.
        ENDDO
      ENDDO
C
      do 11 i=1,nlat  
       prbfrhucgrd(i)         = 0.0
       extnprobgrd(i)         = 0.0
       barf(i)                = 1.0
c
       do j =1,12  
         mlightnggrd(i,j)=0.0
       enddo
c
       do 11 m=1,nmos    
        icountrow(i,m)           = 0
        co2concrow(i,m)          = 0.0 
        npprow(i,m)              = 0.0
        neprow(i,m)              = 0.0
        hetroresrow(i,m)         = 0.0
        autoresrow(i,m)          = 0.0
        soilcresprow(i,m)         = 0.0
        rmrow(i,m)               = 0.0
        rgrow(i,m)               = 0.0
        nbprow(i,m)              = 0.0
        litresrow(i,m)           = 0.0
        socresrow(i,m)           = 0.0
        gpprow(i,m)              = 0.0
        dstcemlsrow(i,m)         = 0.0
        dstcemls3row(i,m)        = 0.0
        litrfallrow(i,m)         = 0.0
        humiftrsrow(i,m)         = 0.0
        canresrow(i,m)           = 0.0
        rmlrow(i,m)              = 0.0
        rmsrow(i,m)              = 0.0
        rmrrow(i,m)              = 0.0  
        lucemcomrow(i,m)         = 0.0
        lucltrinrow(i,m)         = 0.0
        lucsocinrow(i,m)         = 0.0
        burnfracrow(i,m)         = 0.0
        probfirerow(i,m)         = 0.0
        btermrow(i,m)            = 0.0
        ltermrow(i,m)            = 0.0
        mtermrow(i,m)            = 0.0
        cfluxcgrow(i,m)          = 0.0
        cfluxcsrow(i,m)          = 0.0 
c 
        TCANOACCROW_M(I,M)       = 0.0
        UVACCROW_M(I,M)          = 0.0
        VVACCROW_M(I,M)          = 0.0
        TCANOACCROW_OUT(I,M)     = 0.0
c                                         !Rudra addes CH4 related variables on 03/12/2013
        CH4WET1ROW(i,m)          = 0.0
        CH4WET2ROW(i,m)          = 0.0
        WETFDYNROW(i,m)          = 0.0
        CH4DYN1ROW(i,m)          = 0.0
        CH4DYN2ROW(i,m)          = 0.0
        
c
        do j = 1, ignd
           tbaraccrow_m(i,m,j)  = 0.0
        enddo
C
        DO J = 1, ICAN
          ZOLNCROW(I,M,J)        = 0.0
          AILCROW(I,M,J)         = 0.0
          CMASVEGCROW(I,M,J)     = 0.0
          ALVSCTMROW(I,M,J)      = 0.0
          ALIRCTMROW(I,M,J)      = 0.0
          CSUM(I,M,J)            = 0.0
          PAICROW(I,M,J)         = 0.0
          SLAICROW(I,M,J)        = 0.0
          DO K = 1, 3
            RMATCROW(I,M,J,K)    = 0.0
          ENDDO
        ENDDO

        do j = 1, icc
          ailcgrow(i,m,j)        = 0.0
          ailcgsrow(i,m,j)       = 0.0   
          fcancsrow(i,m,j)       = 0.0 
          fcancrow(i,m,j)        = 0.0
          fcancmxrow(i,m,j)      = 0.0
          co2i1cgrow(i,m,j)      = 0.0  
          co2i1csrow(i,m,j)      = 0.0
          co2i2cgrow(i,m,j)      = 0.0
          co2i2csrow(i,m,j)      = 0.0
          ancsvegrow(i,m,j)      = 0.0
          ancgvegrow(i,m,j)      = 0.0 
          rmlcsvegrow(i,m,j)     = 0.0
          rmlcgvegrow(i,m,j)     = 0.0
          stemmassrow(i,m,j)     = 0.0
          rootmassrow(i,m,j)     = 0.0  
          ailcbrow(i,m,j)        = 0.0
          grwtheffrow(i,m,j)     = 0.0
          dvdfcanrow(i,m,j)      = 0.0
          bmasvegrow(i,m,j)      = 0.0
          tltrleafrow(i,m,j)     = 0.0
          tltrstemrow(i,m,j)     = 0.0
          tltrrootrow(i,m,j)     = 0.0
          leaflitrrow(i,m,j)     = 0.0
          roottemprow(i,m,j)     = 0.0
          afrleafrow(i,m,j)      = 0.0
          afrstemrow(i,m,j)      = 0.0
          afrrootrow(i,m,j)      = 0.0
          wtstatusrow(i,m,j)     = 0.0
          ltstatusrow(i,m,j)     = 0.0
          ailcminrow(i,m,j)      = 0.0
          ailcmaxrow(i,m,j)      = 0.0
          pfcancmxrow(i,m,j)     = 0.0
          nfcancmxrow(i,m,j)     = 0.0
          nppvegrow(i,m,j)       = 0.0
          veghghtrow(i,m,j)      = 0.0
          rootdpthrow(i,m,j)     = 0.0
          gleafmasrow(i,m,j)     = 0.0
          bleafmasrow(i,m,j)     = 0.0
          anvegrow(i,m,j)        = 0.0
          rmlvegrow(i,m,j)       = 0.0
c
          rmlvegaccrow(i,m,j)    = 0.0
          rmsvegrow(i,m,j)       = 0.0
          rmrvegrow(i,m,j)       = 0.0
          rgvegrow(i,m,j)        = 0.0
c
          vgbiomas_vegrow(i,m,j) = 0.0
c
          gppvegrow(i,m,j) = 0.0 
          autoresvegrow(i,m,j) = 0.0

          emit_co2row(i,m,j)         =0.0
          emit_corow(i,m,j)          =0.0
          emit_ch4row(i,m,j)         =0.0
          emit_nmhcrow(i,m,j)        =0.0
          emit_h2row(i,m,j)          =0.0
          emit_noxrow(i,m,j)         =0.0
          emit_n2orow(i,m,j)         =0.0
          emit_pm25row(i,m,j)        =0.0
          emit_tpmrow(i,m,j)         =0.0
          emit_tcrow(i,m,j)          =0.0
          emit_ocrow(i,m,j)          =0.0
          emit_bcrow(i,m,j)          =0.0
          burnvegfrow(i,m,j)         =0.0

c
          do k = 1, ignd
            rmatctemrow(i,m,j,k) = 0.0     
          enddo
        enddo
     
        do j = 1, iccp1 
          litrmassrow(i,m,j)    = 0.0
          soilcmasrow(i,m,j)    = 0.0
          hetroresvegrow(i,m,j) = 0.0
          litresvegrow(i,m,j) = 0.0
          soilcresvegrow(i,m,j) = 0.0
          nepvegrow(i,m,j) = 0.0 
          nbpvegrow(i,m,j) = 0.0

        enddo

11     continue
c
      rlim                = abszero
c
c     do some initializations for the reading in of data from files. these
c     initializations primarily affect how the model does a spinup or transient
c     simulation and which years of the input data are being read.

      if (.not. cyclemet) then !transient simulation, set to dummy values
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
        if (ctem_on) then
        open(unit=11,file=argbuff(1:strlen(argbuff))//'.CTM',
     &       status='old')
        endif
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
c     * output files
c
      if (.not. parallelrun) then ! stand alone mode, includes half-hourly and daily output
       OPEN(UNIT=61,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF1_G')  ! DAILY OUTPUT FROM CLASS
       OPEN(UNIT=62,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF2_G')  
       OPEN(UNIT=63,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF3_G')
       OPEN(UNIT=611,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF1_M') ! DAILY OUTPUT FROM CLASS
       OPEN(UNIT=621,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF2_M')
       OPEN(UNIT=631,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF3_M')   
           
       OPEN(UNIT=64,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF4_M')  ! HALF-HOURLY OUTPUT FROM CLASS  
       OPEN(UNIT=65,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF5_M')
       OPEN(UNIT=66,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF6_M')
       OPEN(UNIT=67,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF7_M')
       OPEN(UNIT=68,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF8_M')   
       OPEN(UNIT=69,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF9_M')
       OPEN(UNIT=641,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF4_G') ! HALF-HOURLY OUTPUT FROM CLASS
       OPEN(UNIT=651,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF5_G')
       OPEN(UNIT=661,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF6_G')
       OPEN(UNIT=671,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF7_G')
       OPEN(UNIT=681,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF8_G')   
       OPEN(UNIT=691,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF9_G')
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
       if (ctem_on) then
c       ctem half hourly output files
        open(unit=71, file=argbuff(1:strlen(argbuff))//'.CT01H_M')  
        open(unit=711,file=argbuff(1:strlen(argbuff))//'.CT01H_G')

        if (mosaic) then
c        ctem daily output files (mosaic)
         open(unit=72,file=argbuff(1:strlen(argbuff))//'.CT01D_M')  
         open(unit=73,file=argbuff(1:strlen(argbuff))//'.CT02D_M')
         open(unit=74,file=argbuff(1:strlen(argbuff))//'.CT03D_M')
         open(unit=75,file=argbuff(1:strlen(argbuff))//'.CT04D_M')
         open(unit=76,file=argbuff(1:strlen(argbuff))//'.CT05D_M')
         
         if (dofire .or. lnduseon) then
          open(unit=78,file=argbuff(1:strlen(argbuff))//'.CT06D_M') ! disturbance vars
         endif

        end if ! mosaic

c        ctem daily output files (grid-average)
         open(unit=721,file=argbuff(1:strlen(argbuff))//'.CT01D_G') 
         open(unit=731,file=argbuff(1:strlen(argbuff))//'.CT02D_G')
         open(unit=741,file=argbuff(1:strlen(argbuff))//'.CT03D_G')
         open(unit=751,file=argbuff(1:strlen(argbuff))//'.CT04D_G')
c
         if (dofire .or. lnduseon) then
          open(unit=781,file=argbuff(1:strlen(argbuff))//'.CT06D_G') ! disturbance vars
         endif

         if (compete .or. lnduseon) then
          open(unit=761,file=argbuff(1:strlen(argbuff))//'.CT07D_G') ! competition
         end if
c
          if (dowetlands .or. obswetf) then
          open(unit=762,file=argbuff(1:strlen(argbuff))//'.CT08D_G') ! Methane(Wetland)
          endif 
c
       endif ! ctem_on
      endif ! parallelrun

c     monthly & yearly output for both parallel mode and stand alone mode

C     CLASS MONTHLY OUTPUT FILES 
      OPEN(UNIT=81,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF1M_G') 
      OPEN(UNIT=82,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF2M_G')

C     CLASS YEARLY OUTPUT FILES
      OPEN(UNIT=83,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF1Y_G')

      if (ctem_on) then

       if (mosaic) then

c       Mosaic file names:

c       CTEM monthly output files
        open(unit=84,file=argbuff(1:strlen(argbuff))//'.CT01M_M')

c       CTEM yearly output files
        open(unit=86,file=argbuff(1:strlen(argbuff))//'.CT01Y_M')

        if (dofire .or. lnduseon) then
         open(unit=85,file=argbuff(1:strlen(argbuff))//'.CT06M_M') ! Monthly disturbance
         open(unit=87,file=argbuff(1:strlen(argbuff))//'.CT06Y_M') ! Annual disturbance
        endif 
c
       else

c       Composite file names:

c       CTEM monthly output files
        open(unit=84,file=argbuff(1:strlen(argbuff))//'.CT01M_G')

c       CTEM yearly output files
        open(unit=86,file=argbuff(1:strlen(argbuff))//'.CT01Y_G')

        if (dofire .or. lnduseon) then
         open(unit=85,file=argbuff(1:strlen(argbuff))//'.CT06M_G') ! Monthly disturbance
         open(unit=87,file=argbuff(1:strlen(argbuff))//'.CT06Y_G') ! Annual disturbance
        endif 

      end if !mosiac/composite 
         
        if (compete .or. lnduseon) then
         open(unit=88,file=argbuff(1:strlen(argbuff))//'.CT07M_GM')! ctem pft fractions MONTHLY

         open(unit=89,file=argbuff(1:strlen(argbuff))//'.CT07Y_GM')! ctem pft fractions YEARLY
        endif
c
       if (dowetlands .or. obswetf) then
        open(unit=91,file=argbuff(1:strlen(argbuff))//'.CT08M_G')  !Methane(wetland) MONTHLY
c
        open(unit=92,file=argbuff(1:strlen(argbuff))//'.CT08Y_G')  !Methane(wetland) YEARLY
       endif !dowetlands
      end if !ctem_on
c
C=======================================================================
C
C     * READ AND PROCESS INITIALIZATION AND BACKGROUND INFORMATION.
C     * FIRST, MODEL RUN SPECIFICATIONS.

      READ (10,5010) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
      READ (10,5010) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
      READ (10,5010) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
C
      IF(CTEM_ON) THEN
        READ (11,7010) TITLEC1
        READ (11,7010) TITLEC2
        READ (11,7010) TITLEC3

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
       WRITE(62,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(62,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
C
       IF(IGND.GT.3) THEN
          WRITE(62,6012)
       ELSE
          WRITE(62,6212)
       ENDIF
C
       WRITE(63,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(63,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
C
       IF(IGND.GT.3) THEN
          WRITE(63,6013)
       ELSE
          WRITE(63,6313)
       ENDIF
C
       WRITE(64,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(64,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
       WRITE(64,6014) 
       WRITE(65,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(65,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
C
       IF(IGND.GT.3) THEN
          WRITE(65,6015) 
       ELSE
          WRITE(65,6515)
       ENDIF
C
       WRITE(66,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(66,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
C
       IF(IGND.GT.3) THEN
          WRITE(66,6016)
       ELSE
          WRITE(66,6616)
          WRITE(66,6615)
       ENDIF
C
       WRITE(67,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(67,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
       WRITE(67,6017)
       WRITE(68,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(68,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
       WRITE(68,6018) 
       WRITE(69,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(69,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
       WRITE(69,6019) 
C
6011  FORMAT(2X,'DAY  YEAR  K*  L*  QH  QE  SM  QG  ',
     1          'TR  SWE  DS  WS  AL  ROF')
6012  FORMAT(2X,'DAY  YEAR  TG1  THL1  THI1  TG2  THL2  THI2  ',
     1              'TG3  THL3  THI3  TG4  THL4  THI4  TG5  THL5  ',
     2              'THI5')
6212  FORMAT(2X,'DAY  YEAR  TG1  THL1  THI1  TG2  THL2  THI2  ',
     1              'TG3  THL3  THI3  TCN  RCAN  SCAN  TSN  ZSN')
6013  FORMAT(2X,'DAY  YEAR  TG6  THL6  THI6  TG7  THL7  THI7  ',
     1              'TG8  THL8  THI8  TG9  THL9  THI9  TG10'  ,
     2              'THL10  THI10')
6313  FORMAT(2X,'DAY YEAR KIN LIN TA UV PRES QA PCP EVAP')
6014  FORMAT(2X,'HOUR  MIN  DAY  YEAR  K*  L*  QH  QE  SM  QG  ',
     1          'TR  SWE  DS  WS  AL  ROF  TPN  ZPN CANRES')
6015  FORMAT(2X,'HOUR  MIN  DAY  YEAR  TG1  THL1  THI1  TG2  ',
     1          'THL2  THI2  TG3  THL3  THI3')
6515  FORMAT(2X,'HOUR  MIN  DAY  YEAR  TG1  THL1  THI1  TG2  ',
     1           'THL2  THI2  TG3  THL3  THI3  TCN  RCAN  SCAN  ',
     2           'TSN  ZSN')
6016  FORMAT(2X,'HOUR  MIN  DAY  YEAR  TG4  THL4  THI4  TG5  ',
     1          'THL5  THI5  TG6  THL6  THI6  TG7  ',
     2          'THL7  THI7  TG8  THL8  THI8  TG9  THL9  THI9  ',
     3          'TG10  THL10  THI10  G0  G1  G2  G3  G4  G5  G6  ',
     4          'G7  G8  G9')
6616  FORMAT(2X,'HOUR  MIN  DAY  SWIN  LWIN  PCP  TA  VA  PA  QA')
6615  FORMAT(2X,'IF IGND <= 3, THIS FILE IS EMPTY')
6017  FORMAT(2X,'HOUR  MIN  DAY  YEAR  ',
     1  'TROF     TROO     TROS     TROB      ROF     ROFO   ',
     2  '  ROFS        ROFB         FCS        FGS        FC       FG')
6018  FORMAT(2X,'HOUR  MIN  DAY  YEAR  ',
     1   'FSGV FSGS FSGG FLGV FLGS FLGG HFSC HFSS HFSG ',
     2          'HEVC HEVS HEVG HMFC HMFS HMFG1 HMFG2 HMFG3 ',
     3          'HTCC HTCS HTC1 HTC2 HTC3')
6019  FORMAT(2X,'HOUR  MIN  DAY  YEAR  ',
     1   'PCFC PCLC PCPN PCPG QFCF QFCL QFN QFG QFC1 ',
     2          'QFC2 QFC3 ROFC ROFN ROFO ROF WTRC WTRS WTRG')
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
     1          'TR  SWE  DS  WS  AL  ROF  TPN  ZPN ')
C
C     CTEM FILE TITLES
C
      IF (CTEM_ON) THEN
        WRITE(71,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
        WRITE(71,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
        WRITE(71,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
        WRITE(71,7020)
        WRITE(71,7030)
C
       IF (MOSAIC) THEN
        WRITE(72,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
        WRITE(72,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
        WRITE(72,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
        WRITE(72,7020)
        WRITE(72,7040)
C
        WRITE(73,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
        WRITE(73,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
        WRITE(73,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
        WRITE(73,7020)
        WRITE(73,7050)
C
        WRITE(74,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
        WRITE(74,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
        WRITE(74,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
        WRITE(74,7020)
        WRITE(74,7061)
C
        WRITE(75,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
        WRITE(75,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
        WRITE(75,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
        WRITE(75,7020)
        WRITE(75,7070)
C
        WRITE(76,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
        WRITE(76,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
        WRITE(76,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
        WRITE(76,7020)
        WRITE(76,7080)
C
       IF (DOFIRE .OR. LNDUSEON) THEN
        WRITE(78,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
        WRITE(78,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
        WRITE(78,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
        WRITE(78,7021)
        WRITE(78,7110)
        WRITE(78,7111)
       ENDIF
      END IF !mosaic
C
7010  FORMAT(A80)
7020  FORMAT('CANADIAN TERRESTRIAL ECOSYSTEM MODEL (CTEM) DAILY RESULTS'
     &)
7021  FORMAT('CANADIAN TERRESTRIAL ECOSYSTEM MODEL (CTEM) DAILY ',
     &' DISTURBANCE RESULTS')
7030  FORMAT('HOUR MIN  DAY, An FOR 9 PFTs, RmL FOR 9 PFTs')
7040  FORMAT('  DAY YEAR       GPP       NPP       NEP       NBP',
     &'   AUTORES  HETRORES    LITRES    SOCRES  DSTCEMLS  LITRFALL',
     &'  HUMIFTRS')
7050  FORMAT('  DAY YEAR       RML       RMS       RMR        RG',
     &'  LEAFLITR  TLTRLEAF  TLTRSTEM  TLTRROOT ')
7060  FORMAT('  DAY YEAR  VGBIOMAS   GAVGLAI  GAVGLTMS  GAVGSCMS  ',
     &'TOTCMASS  GLEAFMAS   BLEAFMAS STEMMASS   ROOTMASS  LITRMASS ',
     &' SOILCMAS')
7061  FORMAT('  DAY YEAR  VGBIOMAS   GAVGLAI  GLEAFMAS   BLEAFMAS ',
     & 'STEMMASS   ROOTMASS  LITRMASS SOILCMAS')
7070  FORMAT('  DAY YEAR     AILCG     AILCB    RMATCTEM ',
     &'LAYER 1,2, & 3     VEGHGHT  ROOTDPTH  ROOTTEMP      SLAI')
7075  FORMAT('  DAY YEAR   FRAC #1   FRAC #2   FRAC #3   FRAC #4   ',
     &'FRAC #5   FRAC #6   FRAC #7   FRAC #8   FRAC #9  ',
     &'FRAC #10[%] SUMCHECK')
7080  FORMAT('  DAY YEAR   AFRLEAF   AFRSTEM   AFRROOT  TCANOACC',
     &'  LFSTATUS')
7110  FORMAT('  DAY YEAR   EMIT_CO2',
     &'    EMIT_CO   EMIT_CH4  EMIT_NMHC    EMIT_H2   EMIT_NOX',
     &'   EMIT_N2O  EMIT_PM25   EMIT_TPM    EMIT_TC    EMIT_OC',
     &'    EMIT_BC   BURNFRAC   PROBFIRE   LUCEMCOM   LUCLTRIN',
     &'   LUCSOCIN   GRCLAREA   BTERM   LTERM   MTERM')
7111  FORMAT('               g/m2.D     g/m2.d',
     &'     g/m2.d     g/m2.d     g/m2.d     g/m2.d     g/m2.d',
     &'     g/m2.d     g/m2.d     g/m2.d     g/m2.d     g/m2.d   ',
     &'       %  avgprob/d uMOL-CO2/M2.S KgC/M2.D',
     &'   KgC/M2.D      KM^2    prob/d       prob/d       prob/d')
7112  FORMAT(' DAY  YEAR   CH4WET1    CH4WET2    WETFDYN   CH4DYN1 
     & CH4DYN2 ')
7113  FORMAT('          umolCH4/M2.S    umolCH4/M2.S          umolCH4/M2.S 
     & umolCH4/M2.S')

C
        WRITE(711,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
        WRITE(711,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
        WRITE(711,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
        WRITE(711,7020)
        WRITE(711,7030)
C
        WRITE(721,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
        WRITE(721,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
        WRITE(721,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
        WRITE(721,7020)
        WRITE(721,7040)
C
        WRITE(731,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
        WRITE(731,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
        WRITE(731,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
        WRITE(731,7020)
        WRITE(731,7050)
C
        WRITE(741,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
        WRITE(741,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
        WRITE(741,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
        WRITE(741,7020)
        WRITE(741,7060)
C
        WRITE(751,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
        WRITE(751,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
        WRITE(751,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
        WRITE(751,7020)
        WRITE(751,7070)
C
        IF (COMPETE .OR. LNDUSEON) THEN
         WRITE(761,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
         WRITE(761,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
         WRITE(761,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
         WRITE(761,7020)
         WRITE(761,7075)
        ENDIF
C
       IF (DOFIRE .OR. LNDUSEON) THEN
        WRITE(781,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
        WRITE(781,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
        WRITE(781,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
        WRITE(781,7020)
        WRITE(781,7110)
        WRITE(781,7111)
       ENDIF
c             Methane(Wetland) variables !Rudra
       IF (DOWETLANDS .OR. OBSWETF) THEN     
        WRITE(762,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
        WRITE(762,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
        WRITE(762,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
        WRITE(762,7020)
        WRITE(762,7112)
        WRITE(762,7113)
       ENDIF 

      ENDIF !CTEM_ON
C
      ENDIF !IF NOT PARALLELRUN
 
C     MONTHLY & YEARLY OUTPUT FOR BOTH PARALLEL MODE AND STAND ALONE MODE
      WRITE(81,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
      WRITE(81,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
      WRITE(81,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
      WRITE(81,6021) 
      WRITE(81,6121)
C
      WRITE(82,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
      WRITE(82,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
      WRITE(82,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
      WRITE(82,6022)
      WRITE(82,6122)
C
      WRITE(83,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
      WRITE(83,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
      WRITE(83,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
      WRITE(83,6023)
      WRITE(83,6123)
C
      IF (CTEM_ON) THEN 
         WRITE(84,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
         WRITE(84,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
         WRITE(84,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
         WRITE(84,6024)
         WRITE(84,6124)
         WRITE(84,6224)
C
       IF (DOFIRE .OR. LNDUSEON) THEN
        WRITE(85,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
        WRITE(85,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
        WRITE(85,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
        WRITE(85,6025)
        WRITE(85,6125)
        WRITE(85,6225)
       ENDIF
C
         WRITE(86,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
         WRITE(86,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
         WRITE(86,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
         WRITE(86,6026)
         WRITE(86,6126)
         WRITE(86,6226)
C
        IF (DOFIRE .OR. LNDUSEON) THEN        
         WRITE(87,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
         WRITE(87,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
         WRITE(87,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
         WRITE(87,6027)
         WRITE(87,6127)
         WRITE(87,6227)
        ENDIF
C
        IF (COMPETE .OR. LNDUSEON) THEN
          WRITE(88,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
          WRITE(88,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
          WRITE(88,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
          WRITE(88,6028)
          WRITE(88,6128)
          WRITE(88,6228)
C
          WRITE(89,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
          WRITE(89,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
          WRITE(89,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
          WRITE(89,6029)
          WRITE(89,6129)
          WRITE(89,6229)
        ENDIF !COMPETE
   
         IF (DOWETLANDS .OR. OBSWETF) THEN
          WRITE(91,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
          WRITE(91,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
          WRITE(91,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
          WRITE(91,6028)
          WRITE(91,6230)
          WRITE(91,6231)

          WRITE(92,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
          WRITE(92,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
          WRITE(92,6003) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
          WRITE(92,6028)
          WRITE(92,6232)
          WRITE(92,6233)
         ENDIF 


C
      ENDIF !CTEM_ON
C
6021  FORMAT(2X,'MONTH YEAR  SW     LW      QH      QE    SNOACC    ',  
     &        'WSNOACC    ROFACC      PCP      EVAP       TAIR')
6121  FORMAT(2X,'           W/m2    W/m2    W/m2    W/m2    kg/m2   ',
     &        'kg/m2      mm.mon    mm.mon    mm.mon      degC') 
6022  FORMAT(2X,'MONTH  YEAR  TG1  THL1  THI1     TG2  THL2  THI2',
     &      '     TG3  THL3  THI3')
6122  FORMAT(2X,'             deg  m3/m3  m3/m3   deg  m3/m3  ',
     &      'm3/m3   deg  m3/m3  m3/m3')
6023  FORMAT(2X,'YEAR   SW     LW      QH      QE     ROFACC   ',
     &     ' PCP     EVAP  ')
6123  FORMAT(2X,'      W/m2   W/m2    W/m2    W/m2    mm.yr    ',
     &     'mm.yr    mm.yr')
6024  FORMAT('CANADIAN TERRESTRIAL ECOSYSTEM MODEL (CTEM) MONTHLY ',
     &'RESULTS')
6124  FORMAT('  MONTH  YEAR  LAIMAXG  VGBIOMAS  LITTER    SOIL_C  ', 
     &'  NPP       GPP        NEP       NBP    HETRES',
     &'   AUTORES    LITRES   SOILCRES')
6224  FORMAT('                 m2/m2  Kg C/m2  Kg C/m2   Kg C/m2  ',
     &       'gC/m2.mon  gC/m2.mon  gC/m2.mon  g/m2.mon   g/m2.mon ',
     &       'gC/m2.mon  gC/m2.mon  gC/m2.mon')   
6025  FORMAT('CANADIAN TERRESTRIAL ECOSYSTEM MODEL (CTEM) MONTHLY ',
     &'RESULTS FOR DISTURBANCES')
6125  FORMAT('  MONTH  YEAR  CO2',
     &'        CO        CH4      NMHC       H2       NOX       N2O',
     &'       PM25       TPM        TC        OC        BC  ',
     &' PROBFIRE  LUC_CO2_E  LUC_LTRIN  LUC_SOCIN   BURNFRAC    BTERM',
     &' LTERM   MTERM')
6225  FORMAT('            g/m2.mon  g/m2.mon',
     &'  g/m2.mon  g/m2.mon  g/m2.mon  g/m2.mon  g/m2.mon',
     &'  g/m2.mon  g/m2.mon  g/m2.mon  g/m2.mon  g/m2.mon',
     &'  prob/mon    g C/m2    g C/m2    g C/m2         %  prob/mon',
     &'  prob/mon  prob/mon')  
6026  FORMAT('CANADIAN TERRESTRIAL ECOSYSTEM MODEL (CTEM) YEARLY ',
     &'RESULTS')
6126  FORMAT('  YEAR   LAIMAXG  VGBIOMAS  STEMMASS  ROOTMASS  LITRMASS', 
     &'  SOILCMAS  TOTCMASS  ANNUALNPP ANNUALGPP ANNUALNEP ANNUALNBP',
     &' ANNHETRSP ANAUTORSP ANNLITRES ANSOILCRES')
6226  FORMAT('          m2/m2   Kg C/m2   Kg C/m2   Kg C/m2    Kg C/m2',
     &'  Kg C/m2   Kg C/m2   gC/m2.yr  gC/m2.yr  gC/m2.yr  gC/m2.yr',
     &'  gC/m2.yr  gC/m2.yr  gC/m2.yr  gC/m2.yr')
6027  FORMAT('CANADIAN TERRESTRIAL ECOSYSTEM MODEL (CTEM) YEARLY ',
     &'RESULTS FOR DISTURBANCES')

6127  FORMAT('  YEAR   ANNUALCO2',
     &'  ANNUALCO  ANNUALCH4  ANN_NMHC ANNUAL_H2 ANNUALNOX ANNUALN2O',
     &'  ANN_PM25  ANNUALTPM ANNUAL_TC ANNUAL_OC ANNUAL_BC APROBFIRE',
     &' ANNLUCCO2  ANNLUCLTR ANNLUCSOC ABURNFRAC ANNBTERM ANNLTERM',
     &' ANNMTERM')
6227  FORMAT('         g/m2.yr',
     &'  g/m2.yr  g/m2.yr  g/m2.yr  g/m2.yr  g/m2.yr  g/m2.yr',
     &'  g/m2.yr  g/m2.yr  g/m2.yr  g/m2.yr  g/m2.yr  prob/yr ',
     &'  g/m2.yr  g/m2.yr  g/m2.yr    %     prob/yr  prob/yr',
     &'  prob/yr')
6028  FORMAT('CANADIAN TERRESTRIAL ECOSYSTEM MODEL (CTEM) MONTHLY ',
     &'RESULTS')
6128  FORMAT(' MONTH YEAR  FRAC #1   FRAC #2   FRAC #3   FRAC #4   ',
     &'FRAC #5   FRAC #6   FRAC #7   FRAC #8   FRAC #9   FRAC #10   ',
     &'SUMCHECK, PFT existence for each of the 9 pfts') 
6228  FORMAT('             %         %         %         %         ',
     &'%         %         %         %         %         %          ',
     &'     ')   
6029  FORMAT('CANADIAN TERRESTRIAL ECOSYSTEM MODEL (CTEM) YEARLY ',
     &'RESULTS')
6129  FORMAT('  YEAR   FRAC #1   FRAC #2   FRAC #3   FRAC #4   ',
     &'FRAC #5   FRAC #6   FRAC #7   FRAC #8   FRAC #9   FRAC #10   ',
     &'SUMCHECK, PFT existence for each of the 9 pfts')
6229  FORMAT('         %         %         %         %         ',
     &'%         %         %         %         %         %          ',
     &'     ')
     
6230  FORMAT('MONTH  YEAR   CH4WET1    CH4WET2    WETFDYN   CH4DYN1 
     & CH4DYN2 ')
6231  FORMAT('       gCH4/M2.MON     gCH4/M2.MON        gCH4/M2.MON 
     & gCH4/M2.MON')
6232  FORMAT('  YEAR   CH4WET1    CH4WET2    WETFDYN   CH4DYN1 
     & CH4DYN2 ')
6233  FORMAT('      gCH4/M2.YR      gCH4/M2.YR         gCH4/M2.YR  
     & gCH4/M2.YR ')
 
C
C     CTEM FILE TITLES DONE
C======================= CTEM ========================================== /
C
C=======================================================================

C     BEGIN READ IN OF THE .INI FILE
      
      READ(10,5020) DEGLAT,DEGLON,ZRFMGRD(1),ZRFHGRD(1),ZBLDGRD(1),
     1              GCGRD(1),NLTEST,NMTEST

      JLAT=NINT(DEGLAT)
      RADJGRD(1)=DEGLAT*PI/180.
      DLONGRD(1)=DEGLON
      Z0ORGRD(1)=0.0
      GGEOGRD(1)=0.0
C     GGEOGRD(1)=-0.035

      DO 50 I=1,NLTEST
      DO 50 M=1,NMTEST
          READ(10,5040) (FCANROW(I,M,J),J=1,ICAN+1),(PAMXROW(I,M,J),
     1                  J=1,ICAN)
          READ(10,5040) (LNZ0ROW(I,M,J),J=1,ICAN+1),(PAMNROW(I,M,J),
     1                  J=1,ICAN)
          READ(10,5040) (ALVCROW(I,M,J),J=1,ICAN+1),(CMASROW(I,M,J),
     1                  J=1,ICAN)
          READ(10,5040) (ALICROW(I,M,J),J=1,ICAN+1),(ROOTROW(I,M,J),
     1                  J=1,ICAN)
          READ(10,5030) (RSMNROW(I,M,J),J=1,ICAN),
     1                  (QA50ROW(I,M,J),J=1,ICAN)
          READ(10,5030) (VPDAROW(I,M,J),J=1,ICAN),
     1                  (VPDBROW(I,M,J),J=1,ICAN)
          READ(10,5030) (PSGAROW(I,M,J),J=1,ICAN),
     1                  (PSGBROW(I,M,J),J=1,ICAN)
          READ(10,5040) DRNROW(I,M),SDEPROW(I,M),FAREROW(I,M)
          READ(10,5090) XSLPROW(I,M),GRKFROW(I,M),WFSFROW(I,M),
     1                  WFCIROW(I,M),MIDROW(I,M)
          READ(10,5080) (SANDROW(I,M,J),J=1,3)
          READ(10,5080) (CLAYROW(I,M,J),J=1,3)
          READ(10,5080) (ORGMROW(I,M,J),J=1,3)
          READ(10,5050) (TBARROW(I,M,J),J=1,3),TCANROW(I,M),
     1                  TSNOROW(I,M),TPNDROW(I,M)
          READ(10,5060) (THLQROW(I,M,J),J=1,3),(THICROW(I,M,J),
     1                  J=1,3),ZPNDROW(I,M)
          READ(10,5070) RCANROW(I,M),SCANROW(I,M),SNOROW(I,M),
     1                  ALBSROW(I,M),RHOSROW(I,M),GROROW(I,M)

50    CONTINUE
C
      DO 25 J=1,IGND                     
          READ(10,5002) DELZ(J),ZBOT(J)  
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

      JHHST=JHHSTY*1000+JHHSTD
      JHHEND=JHHENDY*1000+JHHENDD
      JDST=JDSTY*1000+JDSTD
      JDEND=JDENDY*1000+JDENDD

      CLOSE(10)
C
C====================== CTEM =========================================== \
C
c     read from ctem initialization file (.CTM)
c
      if (ctem_on) then
      do 71 i=1,nltest
      do 72 m=1,nmtest
c
c         following three variables are needed to run ctem. 
c         min & max leaf area index are needed to break
c         class lai into dcd and evg for trees (for crops and grasses it
c         doesn't matter much).
c
c         dvdfcanrow is needed to divide needle & broad leaf into dcd and evg,
c         and crops & grasses into c3 and c4 fractions.
c
          read(11,*) (ailcminrow(i,m,j),j=1,icc)
          read(11,*) (ailcmaxrow(i,m,j),j=1,icc)
          read(11,*) (dvdfcanrow(i,m,j),j=1,icc)
c
c         rest of the initialization variables are needed to run ctem.
c         if starting from bare ground initialize all live and dead c pools from zero. suitable values
c         of extnprobgrd and prbfrhucgrd would still be required. set stdalngrd to
c         1 for operation in non-gcm stand alone mode, in the ctem
c         initialization file.
c
          read(11,*) (gleafmasrow(i,m,j),j=1,icc)
          read(11,*) (bleafmasrow(i,m,j),j=1,icc)
          read(11,*) (stemmassrow(i,m,j),j=1,icc)
          read(11,*) (rootmassrow(i,m,j),j=1,icc)

c           If fire and competition are on, save the stemmass and rootmass for
c           use in burntobare subroutine on the first timestep.
            if (dofire .and. compete) then
             do j =1,icc
              pstemmassrow(i,m,j)=stemmassrow(i,m,j)
              pgleafmassrow(i,m,j)=rootmassrow(i,m,j)    
             end do           
            end if

          read(11,*) (litrmassrow(i,m,j),j=1,iccp1)
          read(11,*) (soilcmasrow(i,m,j),j=1,iccp1)
          read(11,*) (lfstatusrow(i,m,j),j=1,icc)
          read(11,*) (pandaysrow(i,m,j),j=1,icc)

72      continue

         read(11,*) (mlightnggrd(i,j),j=1,6)  !mean monthly lightning frequency
         read(11,*) (mlightnggrd(i,j),j=7,12) !flashes/km2.year
         read(11,*) extnprobgrd(i)
         read(11,*) prbfrhucgrd(i)
         read(11,*) stdalngrd(i)

         if (compete .and. inibioclim) then  !read in the bioclimatic parameters
          read(11,*) twarmm(i), tcoldm(i), gdd5(i), aridity(i),
     1              srplsmon(i)
          read(11,*) defctmon(i), anndefct(i), annsrpls(i), 
     1              annpcp(i), dry_season_length(i)
         else if (compete .and. .not. inibioclim) then ! set them to zero
           twarmm(i)=0.0
           tcoldm(i)=0.0
           gdd5(i)=0.0
           aridity(i)=0.0
           srplsmon(i)=0.0
           defctmon(i)=0.0
           anndefct(i)=0.0
           annsrpls(i)=0.0
           annpcp(i)=0.0
           dry_season_length(i) = 0.0
         endif

         if (dowetlands) then      ! Rudra !if true then read wetland fractions
             read(11,*) (wetfrac_sgrd(i,j),j=1,8)
         endif   
71    continue
      close(11)
      endif

c     check that a competition or luc run has the correct number of mosaics
c     if it is not a start_bare run, then nmtest should equal nmos
      if (mosaic .and. (compete .or. lnduseon) 
     &                 .and. .not. start_bare) then
        if (nmtest .ne. nmos) then
           write(6,*)'compete or luc runs that do not start from bare'
           write(6,*)'ground need the number of mosaics to equal icc+1'
           write(6,*)'nmtest = ',nmtest,' nmos = ',nmos
            call xit('runclass36ctem', -2)
        endif
      endif
c
c     if this run uses the competition or lnduseon parameterization and
c     starts from bare ground, set up the model state here. this 
c     overwrites what was read in from the .ini and .ctm files. 
c     for composite runs (the composite set up is after this one for mosaics)
      if ((compete .or. lnduseon) .and. start_bare) then

       if (mosaic) then 

c       store the read-in crop fractions as we keep them even when we start bare. 
c       FLAG: this is setup assuming that crops are in mosaics 6 and 7. JM Apr 9 2014.
         do i=1,nltest
          crop_temp_frac(i,1)=farerow(i,6)
          crop_temp_frac(i,2)=farerow(i,7)
         end do

c       check the number of mosaics that came from the .ini file
        if (nmtest .ne. nmos) then

c        we need to transfer some initial parameterization info to all
c        mosaics, so set all values to that of the first mosaic.
         do i=1,nltest
          do m=nmtest+1,nmos
c
           do j=1,ican
             rsmnrow(i,m,j)=rsmnrow(i,1,j)
             qa50row(i,m,j)=qa50row(i,1,j)
             vpdarow(i,m,j)=vpdarow(i,1,j)
             vpdbrow(i,m,j)=vpdbrow(i,1,j)
             psgarow(i,m,j)=psgarow(i,1,j)
             psgbrow(i,m,j)=psgbrow(i,1,j)
           enddo
c
           drnrow(i,m)=drnrow(i,1)
           sdeprow(i,m)=sdeprow(i,1)
           farerow(i,m)=farerow(i,1)
           xslprow(i,m)=xslprow(i,1)
           grkfrow(i,m)=grkfrow(i,1)
           wfsfrow(i,m)=wfsfrow(i,1)
           wfcirow(i,m)=wfcirow(i,1)
           midrow(i,m)=midrow(i,1)
c
           do j=1,3
            sandrow(i,m,j)=sandrow(i,1,j)
            clayrow(i,m,j)=clayrow(i,1,j)
            orgmrow(i,m,j)=orgmrow(i,1,j)
            tbarrow(i,m,j)=tbarrow(i,1,j)
            thlqrow(i,m,j)=thlqrow(i,1,j)
            thicrow(i,m,j)=thicrow(i,1,j)
           enddo
c
           tcanrow(i,m)=tcanrow(i,1)
           tsnorow(i,m)=tsnorow(i,1)
           tpndrow(i,m)=tpndrow(i,1)
           zpndrow(i,m)=zpndrow(i,1)
           rcanrow(i,m)=rcanrow(i,1)
           scanrow(i,m)=scanrow(i,1)
           snorow(i,m)=snorow(i,1)
           albsrow(i,m)=albsrow(i,1)
           rhosrow(i,m)=rhosrow(i,1)
           grorow(i,m)=grorow(i,1)
           do j=1,icc
             lfstatusrow(i,m,j) = 4
           enddo !j
c
          enddo !m
         enddo !i

c       set the number of mosaics to icc+1        
        nmtest=nmos

        endif  !if (nmtest .ne. nmos)
c
c       set the initial conditions for the pfts
c       (bah, this is such an inelegant way to do this, but oh well...)
c
c       initalize to zero
        fcanrow=0.0
        dvdfcanrow=0.0
        farerow=0.0

        do i=1,nltest
         do m=1,nmtest
c
c         set the seed amount for each pft in its mosaic
          if (compete .or. lnduseon) then
            if (m .lt. icc+1) then
             farerow(i,m)=seed
            else
             farerow(i,m)=1.0 - (real(icc) * seed)
            endif
          endif

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
            litrmassrow(i,m,j)=0. 
            soilcmasrow(i,m,j)=0. 
          enddo

c         initial conditions always required
          dvdfcanrow(i,m,1)=1.0  !ndl
          dvdfcanrow(i,m,3)=1.0  !bdl
          dvdfcanrow(i,m,6)=1.0  !crop
          dvdfcanrow(i,m,8)=1.0  !grasses

c         then adjusted below for the actual mosaic makeup
          if (m .le. 2) then                     !ndl
           fcanrow(i,m,1)=1.0
           if (m .eq. 2) then
             dvdfcanrow(i,m,1)=0.0
             dvdfcanrow(i,m,2)=1.0        
           endif
          elseif (m .ge. 3 .and. m .le. 5) then  !bdl
           fcanrow(i,m,2)=1.0
           if (m .eq. 4) then
             dvdfcanrow(i,m,3)=0.0
             dvdfcanrow(i,m,4)=1.0        
           endif
           if (m .eq. 5) then
             dvdfcanrow(i,m,3)=0.0
             dvdfcanrow(i,m,5)=1.0        
           endif
          elseif (m .eq. 6 .or. m .eq. 7) then  !crop
           fcanrow(i,m,3)=1.0
           if (m .eq. 7) then
             dvdfcanrow(i,m,6)=0.0
             dvdfcanrow(i,m,7)=1.0        
           endif
          elseif (m .eq. 8 .or. m .eq. 9) then  !grasses
           fcanrow(i,m,4)=1.0
           if (m .eq. 9) then
             dvdfcanrow(i,m,8)=0.0
             dvdfcanrow(i,m,9)=1.0        
           endif
          else                                  !bare/urban? 
           fcanrow(i,m,5)=1.0
           endif !mosaic adjustments
         enddo  !m
        enddo  !i


         do i=1,nltest
          farerow(i,6)=crop_temp_frac(i,1)
          farerow(i,7)=crop_temp_frac(i,2)
         end do

      else if (.not. mosaic) then  !set up for composite runs when start_bare is on and compete or landuseon

c       store the read-in crop fractions as we keep them even when we start bare. 
c       FLAG: this is setup assuming that crops are in pft number 6 and 7. JM Apr 9 2014.
         do i=1,nltest
          crop_temp_frac(i,1)=fcanrow(i,1,3)*dvdfcanrow(i,1,6)
          crop_temp_frac(i,2)=fcanrow(i,1,3)*dvdfcanrow(i,1,7)
         end do

c      initalize to zero, these will be filled in by the luc or 
c      competition subroutines.
       fcanrow=0.0
       dvdfcanrow=0.0

       ! Added this as start_bare runs were not properly assigning 
       ! a TCAN on the very first day since the fcanrow was 0. JM Jan 14 2014. 
       do i=1,nltest
        do j=1,iccp1       
           if (j .lt. icc+1) then
            fcanrow(i,1,j)=seed
           else
            fcanrow(i,1,j)=1.0 - (real(icc) * seed)
           endif
        end do
       end do

       do i=1,nltest

c      initial conditions always required
         dvdfcanrow(i,1,1)=1.0  !ndl
         dvdfcanrow(i,1,3)=1.0  !bdl
         dvdfcanrow(i,1,6)=1.0  !crop
         dvdfcanrow(i,1,8)=1.0  !grasses

         do j = 1,icc
           ailcminrow(i,1,j)=0.0
           ailcmaxrow(i,1,j)=0.0
           gleafmasrow(i,1,j)=0.0
           bleafmasrow(i,1,j)=0.0
           stemmassrow(i,1,j)=0.0
           rootmassrow(i,1,j)=0.0
           lfstatusrow(i,1,j)=4
           pandaysrow(i,1,j)=0
         enddo

         lfstatusrow(i,1,1)=2

         do j = 1,iccp1
           litrmassrow(i,1,j)=0.0 
           soilcmasrow(i,1,j)=0.0 
         enddo
       enddo !nltest

         do i=1,nltest
          fcanrow(i,1,3) = crop_temp_frac(i,1) + crop_temp_frac(i,2)
          if (fcanrow(i,1,3) .gt. abszero) then
           dvdfcanrow(i,1,6) = crop_temp_frac(i,1) / fcanrow(i,1,3)
           dvdfcanrow(i,1,7) = crop_temp_frac(i,2) / fcanrow(i,1,3)
          else
           dvdfcanrow(i,1,6) = 1.0
           dvdfcanrow(i,1,7) = 0.0
          end if
         end do

      end if ! mosaic / composite
      end if !if (compete/landuseon .and. start_bare) 

C===================== CTEM =============================================== /
C
      DO 100 I=1,NLTEST
      DO 100 M=1,NMTEST

          TBARROW(I,M,1)=TBARROW(I,M,1)+TFREZ
          TBARROW(I,M,2)=TBARROW(I,M,2)+TFREZ
          TBARROW(I,M,3)=TBARROW(I,M,3)+TFREZ
          TSNOROW(I,M)=TSNOROW(I,M)+TFREZ
          TCANROW(I,M)=TCANROW(I,M)+TFREZ

C         FLAG! TEMP FIX. ON RESTART A SMALL NUMBER OF CELLS (<10) WILL HAVE
C         A RELATIVELY LARGE DIFFERENCE BETWEEN THE CANOPY TEMP AND SNOW TEMP
C         THIS WILL CAUSE THE MODEL TO FAIL RIGHT AWAY. TO PREVENT THIS CHECK
C         IF THE CANOPY TEMP IS WITHIN 5 DEGREES OF THE SNOW TEMP (IF THERE IS
C         SNOW), AND IF SO THEN OVERWRITE THE TCAN WITH 1 DEGREE COLDER THAN THE SNOW TEMP.
C         JM FEB 5 2013
          IF (SNOROW(I,M) .GT. 0.0) then
           IF ( ABS(TCANROW(I,M)-TSNOROW(I,M)) .GT. 5. ) THEN
              TCANROW(I,M)=TSNOROW(I,M) - 0.5
           ENDIF
          ENDIF
 
          TPNDROW(I,M)=TPNDROW(I,M)+TFREZ
          TBASROW(I,M)=TBARROW(I,M,3)
          CMAIROW(I,M)=0.
          WSNOROW(I,M)=0.
          ZSNLROW(I,M)=0.10

C         THIS FIX BELOW IS TO CORRECT A BUG THAT CAUSES A CRASH DUE
C         TO UNREASONABLE CANOPY TEMPERATURES IN THE FIRST YEAR OF A RESTART
C         WITH SNOW ON THE GROUND. NOTE: RUNCLASS.f HAS THIS SAME PROBLEM. JM JAN 2013
          IF (SNOROW(I,M) .GT. 0.) THEN !THERE IS SNOW ON THE GROUND
           TSFSROW(I,M,1)=TBARROW(I,M,1)
           TSFSROW(I,M,2)=TBARROW(I,M,1)
          ELSE ! NO SNOW SO JUST SET THESE TO FREEZING POINT
           TSFSROW(I,M,1)=TFREZ  
           TSFSROW(I,M,2)=TFREZ
          ENDIF

          TSFSROW(I,M,3)=TBARROW(I,M,1)
          TSFSROW(I,M,4)=TBARROW(I,M,1)
          TACROW (I,M)=TCANROW(I,M)
          QACROW (I,M)=0.5E-2

          IF(IGND.GT.3)                                 THEN
              DO 65 J=4,IGND,-1
                  TBARROW(I,M,J)=TBARROW(I,M,3)
                  IF(SDEPROW(I,M).LT.(ZBOT(J-1)+0.001) .AND.
     1                  SANDROW(I,M,3).GT.-2.5)     THEN
                      SANDROW(I,M,J)=-3.0
                      CLAYROW(I,M,J)=-3.0
                      ORGMROW(I,M,J)=-3.0
                      THLQROW(I,M,J)=0.0
                      THICROW(I,M,J)=0.0
                  ELSE
                      SANDROW(I,M,J)=SANDROW(I,M,3)
                      CLAYROW(I,M,J)=CLAYROW(I,M,3)
                      ORGMROW(I,M,J)=ORGMROW(I,M,3)
                      THLQROW(I,M,J)=THLQROW(I,M,3)
                      THICROW(I,M,J)=THICROW(I,M,3)
                  ENDIF
65            CONTINUE
          ENDIF

          DO 75 K=1,6
          DO 75 L=1,50
              ITCTROW(I,M,K,L)=0
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
      do 151 i=1,nltest
          ALVSACC_MO(I)=0.
          ALIRACC_MO(I)=0. 
          FLUTACC_MO(I)=0.
          FSINACC_MO(I)=0.
          FLINACC_MO(I)=0.
          HFSACC_MO(I) =0.
          QEVPACC_MO(I)=0.
          SNOACC_MO(I) =0.
          WSNOACC_MO(I)=0.
          ROFACC_MO(I) =0.
          PREACC_MO(I) =0.
          EVAPACC_MO(I)=0.
          TAACC_MO(I)=0.

          DO 152 J=1,IGND
              TBARACC_MO(I,J)=0.
              THLQACC_MO(I,J)=0.
              THICACC_MO(I,J)=0.
152       CONTINUE  
C
          ALVSACC_YR(I)=0.
          ALIRACC_YR(I)=0. 
          FLUTACC_YR(I)=0.
          FSINACC_YR(I)=0.
          FLINACC_YR(I)=0.
          HFSACC_YR(I) =0.
          QEVPACC_YR(I)=0.
          ROFACC_YR(I) =0.
          PREACC_YR(I) =0.
          EVAPACC_YR(I)=0.
          TAACC_YR(I)=0.
C       
151   CONTINUE
C===================== CTEM =============================================== /

      CALL CLASSB(THPROW,THRROW,THMROW,BIROW,PSISROW,GRKSROW,
     1            THRAROW,HCPSROW,TCSROW,THFCROW,PSIWROW,THLWROW,
     2            DLZWROW,ZBTWROW,ALGWROW,ALGDROW,
     3            SANDROW,CLAYROW,ORGMROW,DELZ,ZBOT,
     4            SDEPROW,ISNDROW,IGDRROW,
     5            NLAT,NMOS,1,NLTEST,NMTEST,IGND,ICTEMMOD)

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
               fcancmxrow(i,m,icountrow(i,m))=max(seed,fcanrow(i,m,j)*
     &         dvdfcanrow(i,m,icountrow(i,m)))
               barf(i) = barf(i) - fcancmxrow(i,m,icountrow(i,m))
              else
               fcancmxrow(i,m,icountrow(i,m))=fcanrow(i,m,j)*
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
          read(12,*) ihour,imin,iday,iyear,fsdown,fdlgrd(i),
     1         pregrd(i),tagrd(i),qagrd(i),uvgrd(i),presgrd(i)
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
     2                     cylucyr,lucyr,fcanrow,farerow,nfcancmxrow,     
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

         do j=1,icc 
          npp_mo_m(i,m,j)=0.0
          gpp_mo_m(i,m,j)=0.0
          nep_mo_m(i,m,j)=0.0
          nbp_mo_m(i,m,j)=0.0
          laimaxg_mo_m(i,m,j)=0.0 
          emit_co2_mo_m(i,m,j)=0.0
          emit_co_mo_m(i,m,j) =0.0
          emit_ch4_mo_m(i,m,j) =0.0
          emit_nmhc_mo_m(i,m,j) =0.0
          emit_h2_mo_m(i,m,j) =0.0
          emit_nox_mo_m(i,m,j) =0.0
          emit_n2o_mo_m(i,m,j) =0.0
          emit_pm25_mo_m(i,m,j) =0.0
          emit_tpm_mo_m(i,m,j) =0.0
          emit_tc_mo_m(i,m,j) =0.0
          emit_oc_mo_m(i,m,j) =0.0
          emit_bc_mo_m(i,m,j) =0.0
          burnfrac_mo_m(i,m,j) =0.0

          laimaxg_yr_m(i,m,j)=0.0    
          npp_yr_m(i,m,j)=0.0
          gpp_yr_m(i,m,j)=0.0
          nep_yr_m(i,m,j)=0.0
          nbp_yr_m(i,m,j)=0.0
          hetrores_yr_m(i,m,j)=0.0
          autores_yr_m(i,m,j)=0.0
          litres_yr_m(i,m,j)=0.0
          soilcres_yr_m(i,m,j)=0.0

          emit_co2_yr_m(i,m,j)=0.0
          emit_co_yr_m(i,m,j)=0.0
          emit_ch4_yr_m(i,m,j)=0.0
          emit_nmhc_yr_m(i,m,j)=0.0
          emit_h2_yr_m(i,m,j)=0.0
          emit_nox_yr_m(i,m,j)=0.0
          emit_n2o_yr_m(i,m,j)=0.0
          emit_pm25_yr_m(i,m,j)=0.0
          emit_tpm_yr_m(i,m,j)=0.0
          emit_tc_yr_m(i,m,j)=0.0
          emit_oc_yr_m(i,m,j)=0.0
          emit_bc_yr_m(i,m,j)=0.0
          burnfrac_yr_m(i,m,j)=0.0

         end do

          nep_mo_m(i,m,iccp1)=0.0
          nbp_mo_m(i,m,iccp1)=0.0
          hetrores_yr_m(i,m,iccp1)=0.0
          litres_yr_m(i,m,iccp1)=0.0
          soilcres_yr_m(i,m,iccp1)=0.0
          nep_yr_m(i,m,iccp1)=0.0
          nbp_yr_m(i,m,iccp1)=0.0

          probfire_mo_m(i,m) =0.0
          luc_emc_mo_m(i,m) =0.0
          lucsocin_mo_m(i,m) =0.0
          lucltrin_mo_m(i,m) =0.0
          bterm_mo_m(i,m)=0.0
          lterm_mo_m(i,m)=0.0
          mterm_mo_m(i,m)=0.0

          probfire_yr_m(i,m)=0.0
          luc_emc_yr_m(i,m)=0.0
          lucsocin_yr_m(i,m)=0.0
          lucltrin_yr_m(i,m)=0.0
          bterm_yr_m(i,m)=0.0 
          lterm_yr_m(i,m)=0.0
          mterm_yr_m(i,m)=0.0
c                                      !CH4(wetland) related variables !Rudra 04/12/2013
          ch4wet1_mo_m(i,m)  =0.0
          ch4wet2_mo_m(i,m)  =0.0
          wetfdyn_mo_m(i,m)  =0.0
          ch4dyn1_mo_m(i,m)  =0.0
          ch4dyn2_mo_m(i,m)  =0.0

          ch4wet1_yr_m(i,m)  =0.0
          ch4wet2_yr_m(i,m)  =0.0
          wetfdyn_yr_m(i,m)  =0.0
          ch4dyn1_yr_m(i,m)  =0.0 
          ch4dyn2_yr_m(i,m)  =0.0     
c
        endif ! ctem_on
c
115   continue
c
      do 117 i = 1,nltest
        do 117 m = 1,nmtest
         gavgltmsrow(i,m)=gavgltmsrow(i,m)+ (1.0-fcanrow(i,m,1)-
     &       fcanrow(i,m,2)-
     &    fcanrow(i,m,3)-fcanrow(i,m,4))*litrmassrow(i,m,icc+1)
         gavgscmsrow(i,m)=gavgscmsrow(i,m)+ (1.0-fcanrow(i,m,1)-
     &   fcanrow(i,m,2)-
     &    fcanrow(i,m,3)-fcanrow(i,m,4))*soilcmasrow(i,m,icc+1)
c
117   continue
c
C===================== CTEM =============================================== /

      CALL GATPREP(ILMOS,JLMOS,IWMOS,JWMOS,
     1             NML,NMW,GCGRD,FAREROW,MIDROW,
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
     9      fcancmxrow,zbtwrow,dlzwrow,sdeprow,ailcgrow,ailcbrow,
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
c    find the wilting point and field capacity for classt
c    (it would be preferable to have this in a subroutine 
c    rather than here. jm sep 06/12)
c
!       FLAG this can be removed once the wilting point matric pot limit is decided. JM Jan 14 2015.
!        do 119 i = 1,ilg
!         do 119 j = 1,ignd

!           psisat(i,j)= (10.0**(-0.0131*sandgat(i,j)+1.88))/100.0
!           grksat(i,j)= (10.0**(0.0153*sandgat(i,j)-0.884))*7.0556e-6
!           thpor(i,j) = (-0.126*sandgat(i,j)+48.9)/100.0
!           bterm(i,j)     = 0.159*claygat(i,j)+2.91

!           wiltsm(i,j) = (150./psisat(i,j))**(-1.0/bterm(i,j))
!           wiltsm(i,j) = thpor(i,j) * wiltsm(i,j)

!           fieldsm(i,j) = (1.157e-09/grksat(i,j))**
!     &      (1.0/(2.0*bterm(i,j)+3.0))
!           fieldsm(i,j) = thpor(i,j) *  fieldsm(i,j)

!119    continue
c
      call ctems1(gleafmasrow,bleafmasrow,stemmassrow,rootmassrow,
     1      fcancmxrow,zbtwrow,dlzwrow,sdeprow,ailcgrow,ailcbrow,
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
          read(12,*,end=999) ihour,imin,iday,iyear,fsdown,fdlgrd(i),
     1         pregrd(i),tagrd(i),qagrd(i),uvgrd(i),presgrd(i)

         enddo
        enddo

c       back up one space in the met file so it is ready for the next readin
c       but only if it was read in during the loop above.    
        if (metcylyrst .ne. -9999) backspace(12)

c  /------------------Rudra----------------/

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


c  \------------------Rudra---------------\     


      met_rewound = .false.

      endif

C===================== CTEM ============================================ /
C
C     * READ IN METEOROLOGICAL FORCING DATA FOR CURRENT TIME STEP;
C     * CALCULATE SOLAR ZENITH ANGLE AND COMPONENTS OF INCOMING SHORT-
C     * WAVE RADIATION FLUX; ESTIMATE FLUX PARTITIONS IF NECESSARY.
C
      N=N+1
C
      DO 250 I=1,NLTEST
C         THIS READS IN ONE 30 MIN SLICE OF MET DATA, WHEN IT REACHES 
C         THE END OF FILE IT WILL GO TO 999. !formatting was 5300
          READ(12,*,END=999) IHOUR,IMIN,IDAY,IYEAR,FSDOWN,FDLGRD(I),
     1         PREGRD(I),TAGRD(I),QAGRD(I),UVGRD(I),PRESGRD(I)

c         /---------------Rudra-----------------/

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
         
c         \----------------Rudra---------------\ 

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
          FSVHGRD(I)=0.5*FSDOWN
          FSIHGRD(I)=0.5*FSDOWN
          TAGRD(I)=TAGRD(I)+TFREZ
          ULGRD(I)=UVGRD(I)
          VLGRD(I)=0.0
          VMODGRD(I)=UVGRD(I) 
250   CONTINUE

C
      DAY=REAL(IDAY)+(REAL(IHOUR)+REAL(IMIN)/60.)/24.
      DECL=SIN(2.*PI*(284.+DAY)/365.)*23.45*PI/180.
      HOUR=(REAL(IHOUR)+REAL(IMIN)/60.)*PI/12.-PI
      COSZ=SIN(RADJGRD(1))*SIN(DECL)+COS(RADJGRD(1))*COS(DECL)*COS(HOUR)

      DO 300 I=1,NLTEST
          CSZGRD(I)=SIGN(MAX(ABS(COSZ),1.0E-3),COSZ)
          IF(PREGRD(I).GT.0.) THEN
              XDIFFUS(I)=1.0
          ELSE
              XDIFFUS(I)=MAX(0.0,MIN(1.0-0.9*COSZ,1.0))
          ENDIF
          FCLOGRD(I)=XDIFFUS(I)
300   CONTINUE
C
C===================== CTEM ============================================ \
C
      if (iday.eq.1.and.ihour.eq.0.and.imin.eq.0) then
c
c      if popdon=true
c      calculate fire extinguishing probability and 
c      probability of fire due to human causes
c      from population density input data
c      overwrite extnprobgrd(i) and prbfrhucgrd(i) 
c      read from .ctm file
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
      CALL CLASSI(VPDGRD,TADPGRD,PADRGRD,RHOAGRD,RHSIGRD,
     1            RPCPGRD,TRPCGRD,SPCPGRD,TSPCGRD,TAGRD,QAGRD,
     2            PREGRD,RPREGRD,SPREGRD,PRESGRD,
     3            IPCP,NLAT,1,NLTEST)
C
      CUMSNO=CUMSNO+SPCPGRD(1)*RHSIGRD(1)*DELT
C
      CALL GATPREP(ILMOS,JLMOS,IWMOS,JWMOS,
     1             NML,NMW,GCGRD,FAREROW,MIDROW,
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
     8             GRKSGAT,THRAGAT,HCPSGAT,TCSGAT,IGDRGAT,
     9             THFCGAT,PSIWGAT,DLZWGAT,ZBTWGAT,VMODGAT,
     A             ZSNLGAT,ZPLGGAT,ZPLSGAT,TACGAT, QACGAT,
     B             DRNGAT, XSLPGAT,GRKFGAT,WFSFGAT,WFCIGAT,
     C             ALGWGAT,ALGDGAT,ASVDGAT,ASIDGAT,AGVDGAT,
     D             AGIDGAT,ISNDGAT,RADJGAT,ZBLDGAT,Z0ORGAT,
     E             ZRFMGAT,ZRFHGAT,ZDMGAT, ZDHGAT, FSVHGAT,
     F             FSIHGAT,CSZGAT, FDLGAT, ULGAT,  VLGAT,  
     G             TAGAT,  QAGAT,  PRESGAT,PREGAT, PADRGAT,
     H             VPDGAT, TADPGAT,RHOAGAT,RPCPGAT,TRPCGAT,
     I             SPCPGAT,TSPCGAT,RHSIGAT,FCLOGAT,DLONGAT,
     J             GGEOGAT,THLWGAT,
     K             ILMOS,JLMOS,IWMOS,JWMOS,
     L             NML,NLAT,NMOS,ILG,IGND,ICAN,ICAN+1,
     M             TBARROW,THLQROW,THICROW,TPNDROW,ZPNDROW,
     N             TBASROW,ALBSROW,TSNOROW,RHOSROW,SNOROW, 
     O             TCANROW,RCANROW,SCANROW,GROROW, CMAIROW,
     P             FCANROW,LNZ0ROW,ALVCROW,ALICROW,PAMXROW,
     Q             PAMNROW,CMASROW,ROOTROW,RSMNROW,QA50ROW,
     R             VPDAROW,VPDBROW,PSGAROW,PSGBROW,PAIDROW,
     S             HGTDROW,ACVDROW,ACIDROW,TSFSROW,WSNOROW,
     T             THPROW, THRROW, THMROW, BIROW,  PSISROW,
     U             GRKSROW,THRAROW,HCPSROW,TCSROW, IGDRROW,
     V             THFCROW,PSIWROW,DLZWROW,ZBTWROW,VMODGRD,
     W             ZSNLROW,ZPLGROW,ZPLSROW,TACROW, QACROW,
     X             DRNROW, XSLPROW,GRKFROW,WFSFROW,WFCIROW,
     Y             ALGWROW,ALGDROW,ASVDROW,ASIDROW,AGVDROW,
     Z             AGIDROW,ISNDROW,RADJGRD,ZBLDGRD,Z0ORGRD,
     +             ZRFMGRD,ZRFHGRD,ZDMGRD, ZDHGRD, FSVHGRD,
     +             FSIHGRD,CSZGRD, FDLGRD, ULGRD,  VLGRD,  
     +             TAGRD,  QAGRD,  PRESGRD,PREGRD, PADRGRD,
     +             VPDGRD, TADPGRD,RHOAGRD,RPCPGRD,TRPCGRD,
     +             SPCPGRD,TSPCGRD,RHSIGRD,FCLOGRD,DLONGRD,
     +             GGEOGRD, THLWROW  )
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
!     &      WETFRACGAT, WETFRAC_SGAT,
     &      CH4WET1GAT, CH4WET2GAT, 
     &      WETFDYNGAT, CH4DYN1GAT,  CH4DYN2GAT,
c
     r      ilmos,       jlmos,       iwmos,        jwmos,
     s      nml,      fcancmxrow,  rmatcrow,    zolncrow,     paicrow,
     v      ailcrow,     ailcgrow,    cmasvegcrow,  slaicrow,
     w      ailcgsrow,   fcancsrow,   fcancrow,     rmatctemrow,
     x      co2concrow,  co2i1cgrow,  co2i1csrow,   co2i2cgrow,
     y      co2i2csrow,  xdiffus,     slairow,      cfluxcgrow,
     z      cfluxcsrow,  ancsvegrow,  ancgvegrow,   rmlcsvegrow,
     1      rmlcgvegrow, canresrow,   sdeprow,
     2      sandrow,     clayrow,     orgmrow,
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
     l      farerow,     gavgscmsrow, rmlvegaccrow, pftexistrow,
     &      rmsvegrow,   rmrvegrow,   rgvegrow,    vgbiomas_vegrow,
     &      gppvegrow,   nepvegrow,   ailcminrow,   ailcmaxrow,
     &      emit_co2row,  emit_corow, emit_ch4row,  emit_nmhcrow,
     &      emit_h2row,   emit_noxrow,emit_n2orow,  emit_pm25row,
     &      emit_tpmrow,  emit_tcrow, emit_ocrow,   emit_bcrow,
     &      btermrow,     ltermrow,   mtermrow,
     &      nbpvegrow,    hetroresvegrow, autoresvegrow,litresvegrow,
     &      soilcresvegrow, burnvegfrow, pstemmassrow, pgleafmassrow,
!     &      WETFRACROW, WETFRAC_SROW,
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
     9                ZPLIMC, ZPLIMG, ZPLMCS, ZPLMGS, TRSNOW, ZSNOW,  
     A                WSNOGAT,ALVSGAT,ALIRGAT,HTCCGAT,HTCSGAT,HTCGAT, 
     B                WTRCGAT,WTRSGAT,WTRGGAT,CMAIGAT,FSNOGAT,
     C                FCANGAT,LNZ0GAT,ALVCGAT,ALICGAT,PAMXGAT,PAMNGAT,
     D                CMASGAT,ROOTGAT,RSMNGAT,QA50GAT,VPDAGAT,VPDBGAT,
     E                PSGAGAT,PSGBGAT,PAIDGAT,HGTDGAT,ACVDGAT,ACIDGAT, 
     F                ASVDGAT,ASIDGAT,AGVDGAT,AGIDGAT,ALGWGAT,ALGDGAT, 
     G                THLQGAT,THICGAT,TBARGAT,RCANGAT,SCANGAT,TCANGAT,   
     H                GROGAT, SNOGAT, TSNOGAT,RHOSGAT,ALBSGAT,ZBLDGAT,
     I                Z0ORGAT,ZSNLGAT,ZPLGGAT,ZPLSGAT,
     J                FCLOGAT,TAGAT,  VPDGAT, RHOAGAT,CSZGAT, 
     K                FSVHGAT,RADJGAT,DLONGAT,RHSIGAT,DELZ,   DLZWGAT,
     L                ZBTWGAT,THPGAT, THMGAT, PSISGAT,BIGAT,  PSIWGAT,
     M                HCPSGAT,ISNDGAT,
     P                FCANCMXGAT,ICC,ICTEMMOD,RMATCGAT,ZOLNCGAT, 
     Q                CMASVEGCGAT,AILCGAT,PAICGAT,L2MAX, NOL2PFTS,
     R                SLAICGAT,AILCGGAT,AILCGSGAT,FCANCGAT,FCANCSGAT,
     R                IDAY,   ILG,    1,      NML,      
     N                JLAT,N, ICAN,   ICAN+1, IGND,   IDISP,  IZREF,
     O                IWF,    IPAI,   IHGT,   IALC,   IALS,   IALG,
     P                ALVSCTMGAT, ALIRCTMGAT )
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
     7  ITCTGAT,CDHGAT, CDMGAT, HFSGAT, TFXGAT, QEVPGAT,QFSGAT, QFXGAT, 
     8  PETGAT, GAGAT,  EFGAT,  GTGAT,  QGGAT,  SFCTGAT,SFCUGAT,SFCVGAT,
     9  SFCQGAT,SFRHGAT,FSGVGAT,FSGSGAT,FSGGGAT,FLGVGAT,FLGSGAT,FLGGGAT,
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
     L  ZSNOW,  TRSNOW, RHOSGAT,WSNOGAT,THPGAT, THRGAT, THMGAT, THFCGAT,
     M  RADJGAT,PREGAT, HCPSGAT,TCSGAT, TSFSGAT,DELZ,   DLZWGAT,ZBTWGAT,
     N  FTEMP,  FVAP,   RIB,    ISNDGAT,
     O  AILCGGAT,  AILCGSGAT, FCANCGAT,FCANCSGAT,CO2CONCGAT,CO2I1CGGAT,
     P  CO2I1CSGAT,CO2I2CGGAT,CO2I2CSGAT,CSZGAT,XDIFFUSGAT,SLAIGAT,ICC,
     Q  ICTEMMOD,RMATCTEMGAT,FCANCMXGAT,L2MAX,  NOL2PFTS,CFLUXCGGAT,
     R  CFLUXCSGAT,ANCSVEGGAT,ANCGVEGGAT,RMLCSVEGGAT,RMLCGVEGGAT,
     S  THLWGAT,ITC,ITCG,ITG,    ILG,    1,NML,  JLAT,N, ICAN,   
     T  IGND,   IZREF,  ISLFD,  NLANDCS,NLANDGS,NLANDC, NLANDG, NLANDI) 
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
     7                  TROOGAT,TROSGAT,TROBGAT,QFSGAT, 
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

C-----------------------------------------------------------------------
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
C-----------------------------------------------------------------------
C
C===================== CTEM ============================================ \
C
c     use net photosynthetic rates from canopy over snow and canopy over 
c     ground sub-areas to find average net photosynthetic rate for each
c     pft. and similarly for leaf respiration.
c
      if (ctem_on) then
        do 605 j = 1, icc
          do 610 i = 1, nml
            if ( (fcancgat(i,j)+fcancsgat(i,j)).gt.abszero) then
              anveggat(i,j)=(fcancgat(i,j)*ancgveggat(i,j) + 
     &                   fcancsgat(i,j)*ancsveggat(i,j)) / 
     &                   (fcancgat(i,j)+fcancsgat(i,j))   
              rmlveggat(i,j)=(fcancgat(i,j)*rmlcgveggat(i,j) + 
     &                    fcancsgat(i,j)*rmlcsveggat(i,j)) / 
     &                    (fcancgat(i,j)+fcancsgat(i,j))   
            else
              anveggat(i,j)=0.0
              rmlveggat(i,j)=0.0
            endif
610       continue
605     continue
      endif 

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
          fsinacc_gat(i)=fsinacc_gat(i)+fsdown
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
            else if(iday.ge.350.or.iday.lt.14)then ! mid nov - mid dec
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
     9            extnprobgat,   stdalngat,tbaraccgat_m,   
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
      CALL CLASSS (TBARROW,THLQROW,THICROW,GFLXROW,TSFSROW,
     1             TPNDROW,ZPNDROW,TBASROW,ALBSROW,TSNOROW,
     2             RHOSROW,SNOROW, TCANROW,RCANROW,SCANROW,
     3             GROROW, CMAIROW,TACROW, QACROW, WSNOROW,
     4             ILMOS,JLMOS,IWMOS,JWMOS,
     5             NML,NLAT,NMOS,ILG,IGND,ICAN,ICAN+1,
     6             TBARGAT,THLQGAT,THICGAT,GFLXGAT,TSFSGAT,
     7             TPNDGAT,ZPNDGAT,TBASGAT,ALBSGAT,TSNOGAT,
     8             RHOSGAT,SNOGAT, TCANGAT,RCANGAT,SCANGAT,
     9             GROGAT, CMAIGAT,TACGAT, QACGAT, WSNOGAT)

C    * SCATTER OPERATION ON DIAGNOSTIC VARIABLES SPLIT OUT OF 
C    * CLASSS FOR CONSISTENCY WITH GCM APPLICATIONS.
C
      DO 380 K=1,NML
          CDHROW (ILMOS(K),JLMOS(K))=CDHGAT (K)  
          CDMROW (ILMOS(K),JLMOS(K))=CDMGAT (K)  
          HFSROW (ILMOS(K),JLMOS(K))=HFSGAT (K)  
          TFXROW (ILMOS(K),JLMOS(K))=TFXGAT (K)  
          QEVPROW(ILMOS(K),JLMOS(K))=QEVPGAT(K)  
          QFSROW (ILMOS(K),JLMOS(K))=QFSGAT (K)  
          QFXROW (ILMOS(K),JLMOS(K))=QFXGAT (K)  
          PETROW (ILMOS(K),JLMOS(K))=PETGAT (K)  
          GAROW  (ILMOS(K),JLMOS(K))=GAGAT  (K)  
          EFROW  (ILMOS(K),JLMOS(K))=EFGAT  (K)  
          GTROW  (ILMOS(K),JLMOS(K))=GTGAT  (K)  
          QGROW  (ILMOS(K),JLMOS(K))=QGGAT  (K)  
          ALVSROW(ILMOS(K),JLMOS(K))=ALVSGAT(K)  
          ALIRROW(ILMOS(K),JLMOS(K))=ALIRGAT(K)  
          SFCTROW(ILMOS(K),JLMOS(K))=SFCTGAT(K)  
          SFCUROW(ILMOS(K),JLMOS(K))=SFCUGAT(K)  
          SFCVROW(ILMOS(K),JLMOS(K))=SFCVGAT(K)  
          SFCQROW(ILMOS(K),JLMOS(K))=SFCQGAT(K)  
          FSNOROW(ILMOS(K),JLMOS(K))=FSNOGAT(K)  
          FSGVROW(ILMOS(K),JLMOS(K))=FSGVGAT(K)  
          FSGSROW(ILMOS(K),JLMOS(K))=FSGSGAT(K)  
          FSGGROW(ILMOS(K),JLMOS(K))=FSGGGAT(K)  
          FLGVROW(ILMOS(K),JLMOS(K))=FLGVGAT(K)  
          FLGSROW(ILMOS(K),JLMOS(K))=FLGSGAT(K)  
          FLGGROW(ILMOS(K),JLMOS(K))=FLGGGAT(K)  
          HFSCROW(ILMOS(K),JLMOS(K))=HFSCGAT(K)  
          HFSSROW(ILMOS(K),JLMOS(K))=HFSSGAT(K)  
          HFSGROW(ILMOS(K),JLMOS(K))=HFSGGAT(K)  
          HEVCROW(ILMOS(K),JLMOS(K))=HEVCGAT(K)  
          HEVSROW(ILMOS(K),JLMOS(K))=HEVSGAT(K)  
          HEVGROW(ILMOS(K),JLMOS(K))=HEVGGAT(K)  
          HMFCROW(ILMOS(K),JLMOS(K))=HMFCGAT(K)  
          HMFNROW(ILMOS(K),JLMOS(K))=HMFNGAT(K)  
          HTCCROW(ILMOS(K),JLMOS(K))=HTCCGAT(K)  
          HTCSROW(ILMOS(K),JLMOS(K))=HTCSGAT(K)  
          PCFCROW(ILMOS(K),JLMOS(K))=PCFCGAT(K)  
          PCLCROW(ILMOS(K),JLMOS(K))=PCLCGAT(K)  
          PCPNROW(ILMOS(K),JLMOS(K))=PCPNGAT(K)  
          PCPGROW(ILMOS(K),JLMOS(K))=PCPGGAT(K)  
          QFGROW (ILMOS(K),JLMOS(K))=QFGGAT (K)  
          QFNROW (ILMOS(K),JLMOS(K))=QFNGAT (K)  
          QFCLROW(ILMOS(K),JLMOS(K))=QFCLGAT(K)  
          QFCFROW(ILMOS(K),JLMOS(K))=QFCFGAT(K)  
          ROFROW (ILMOS(K),JLMOS(K))=ROFGAT (K)  
          ROFOROW(ILMOS(K),JLMOS(K))=ROFOGAT(K)  
          ROFSROW(ILMOS(K),JLMOS(K))=ROFSGAT(K)  
          ROFBROW(ILMOS(K),JLMOS(K))=ROFBGAT(K)  
          TROFROW(ILMOS(K),JLMOS(K))=TROFGAT(K)  
          TROOROW(ILMOS(K),JLMOS(K))=TROOGAT(K)  
          TROSROW(ILMOS(K),JLMOS(K))=TROSGAT(K)  
          TROBROW(ILMOS(K),JLMOS(K))=TROBGAT(K)  
          ROFCROW(ILMOS(K),JLMOS(K))=ROFCGAT(K)  
          ROFNROW(ILMOS(K),JLMOS(K))=ROFNGAT(K)  
          ROVGROW(ILMOS(K),JLMOS(K))=ROVGGAT(K)  
          WTRCROW(ILMOS(K),JLMOS(K))=WTRCGAT(K)  
          WTRSROW(ILMOS(K),JLMOS(K))=WTRSGAT(K)  
          WTRGROW(ILMOS(K),JLMOS(K))=WTRGGAT(K)  
          DRROW  (ILMOS(K),JLMOS(K))=DRGAT  (K)  
          WTABROW(ILMOS(K),JLMOS(K))=WTABGAT(K)  
          ILMOROW(ILMOS(K),JLMOS(K))=ILMOGAT(K)  
          UEROW  (ILMOS(K),JLMOS(K))=UEGAT(K)  
          HBLROW (ILMOS(K),JLMOS(K))=HBLGAT(K)  
380   CONTINUE
C
      DO 390 L=1,IGND
      DO 390 K=1,NML
          HMFGROW(ILMOS(K),JLMOS(K),L)=HMFGGAT(K,L)
          HTCROW (ILMOS(K),JLMOS(K),L)=HTCGAT (K,L)
          QFCROW (ILMOS(K),JLMOS(K),L)=QFCGAT (K,L)
390   CONTINUE
C
      DO 430 M=1,50
          DO 420 L=1,6
              DO 410 K=1,NML
                  ITCTROW(ILMOS(K),JLMOS(K),L,M)=ITCTGAT(K,L,M)
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
     6      rmlcgvegrow, canresrow,   sdeprow,
     7      sandrow,     clayrow,     orgmrow,
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
     q      farerow,     gavgscmsrow, tcanoaccrow_out,
     &      rmlvegaccrow, rmsvegrow,  rmrvegrow,    rgvegrow,
     &      vgbiomas_vegrow,gppvegrow,nepvegrow,ailcminrow,ailcmaxrow,
     &      fcanrow,      pftexistrow,
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

6100  FORMAT(1X,I4,I5,9F8.2,2F8.3,F12.4,2(A6,I2))
6200  FORMAT(1X,I4,I5,3(F8.2,2F6.3),F8.2,2F8.4,F8.2,F8.3,2(A6,I2))
6201  FORMAT(1X,I4,I5,5(F7.2,2F6.3),2(A6,I2))
6300  FORMAT(1X,I4,I5,3F9.2,F8.2,F10.2,E12.3,2F12.3,A6,I2)   
6400  FORMAT(1X,I2,I3,I5,I6,9F8.2,2F7.3,E11.3,F8.2,F12.4,F9.2,2(A6,I2))
6500  FORMAT(1X,I2,I3,I5,I6,3(F7.2,2F6.3),F8.2,2F8.4,F8.2,F8.3,
     &    2(A6,I2))
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
        fsstar_g    =0.0
        flstar_g    =0.0
        qh_g        =0.0
        qe_g        =0.0
        snomlt_g    =0.0
        beg_g       =0.0
        gtout_g     =0.0
        snorow_g(i) =0.0
        rhosrow_g(i)=0.0
        wsnorow_g(i)=0.0
        altot_g     =0.0
        rofrow_g(i) =0.0
        tpn_g       =0.0
        zpndrow_g(i)=0.0
c
        do j=1,ignd
         tbarrow_g(i,j)=0.0
         thlqrow_g(i,j)=0.0
         thicrow_g(i,j)=0.0
         gflxrow_g(i,j)=0.0
        enddo
c
        tcn_g=0.0
        rcanrow_g(i) =0.0
        scanrow_g(i) =0.0
        tsn_g=0.0
        zsn_g=0.0
        trofrow_g(i)=0.0
        troorow_g(i)=0.0
        trosrow_g(i)=0.0
        trobrow_g(i)=0.0
        roforow_g(i)=0.0
        rofsrow_g(i)=0.0
        rofbrow_g(i)=0.0
        fsgvrow_g(i)=0.0
        fsgsrow_g(i)=0.0
        fsggrow_g(i)=0.0
        flgvrow_g(i)=0.0
        flgsrow_g(i)=0.0
        flggrow_g(i)=0.0
        hfscrow_g(i)=0.0
        hfssrow_g(i)=0.0
        hfsgrow_g(i)=0.0
        hevcrow_g(i)=0.0
        hevsrow_g(i)=0.0
        hevgrow_g(i)=0.0
        hmfcrow_g(i)=0.0
        hmfnrow_g(i)=0.0
c
        do j=1,ignd    
         hmfgrow_g(i,j)=0.0
         htcrow_g(i,j)=0.0
        enddo
c               
        htccrow_g(i)=0.0
        htcsrow_g(i)=0.0
        pcfcrow_g(i)=0.0
        pclcrow_g(i)=0.0
        pcpnrow_g(i)=0.0
        pcpgrow_g(i)=0.0
        qfcfrow_g(i)=0.0
        qfclrow_g(i)=0.0
        qfnrow_g(i)=0.0
        qfgrow_g(i)=0.0
c       
        do j=1,ignd    
         qfcrow_g(i,j)=0.0
        enddo        
c       
        rofcrow_g(i)=0.0
        rofnrow_g(i)=0.0
        wtrcrow_g(i)=0.0
        wtrsrow_g(i)=0.0
        wtrgrow_g(i)=0.0
c
       if (ctem_on) then
          do j=1,icc
            anvegrow_g(i,j)=0.0
            rmlvegrow_g(i,j)=0.0
          enddo
       endif   !if (ctem_on) 
c
C===================== CTEM =====================================/
C
       DO 425 M=1,NMTEST
          IF(FSDOWN.GT.0.0) THEN
              ALTOT=(ALVSROW(I,M)+ALIRROW(I,M))/2.0
          ELSE
              ALTOT=0.0
          ENDIF
          FSSTAR=FSDOWN*(1.0-ALTOT)
          FLSTAR=FDLGRD(I)-SBC*GTROW(I,M)**4
          QH=HFSROW(I,M)
          QE=QEVPROW(I,M)
          BEG=FSSTAR+FLSTAR-QH-QE
C         BEG=GFLXGAT(1,1)
          SNOMLT=HMFNROW(I,M)
          IF(RHOSROW(I,M).GT.0.0) THEN
              ZSN=SNOROW(I,M)/RHOSROW(I,M)
          ELSE
              ZSN=0.0
          ENDIF
          IF(TCANROW(I,M).GT.0.01) THEN
              TCN=TCANROW(I,M)-TFREZ
          ELSE
              TCN=0.0
          ENDIF
          IF(TSNOROW(I,M).GT.0.01) THEN
              TSN=TSNOROW(I,M)-TFREZ
          ELSE
              TSN=0.0
          ENDIF
          IF(TPNDROW(I,M).GT.0.01) THEN
              TPN=TPNDROW(I,M)-TFREZ
          ELSE
              TPN=0.0
          ENDIF
          GTOUT=GTROW(I,M)-TFREZ
C 
C===================== CTEM =====================================\
c         start writing output
c
          iyd=iyear*1000+iday                         
          if ((iyd.ge.jhhst).and.(iyd.le.jhhend)) then 
C===================== CTEM =====================================/
          WRITE(64,6400) IHOUR,IMIN,IDAY,IYEAR,FSSTAR,FLSTAR,QH,QE,
     1                   SNOMLT,BEG,GTOUT,SNOROW(I,M),RHOSROW(I,M),
     2                   WSNOROW(I,M),ALTOT,ROFROW(I,M),
     3                   TPN,ZPNDROW(I,M),CANRESROW(I,M),' TILE ',M
          IF(IGND.GT.3) THEN
C===================== CTEM =====================================\

              write(65,6500) ihour,imin,iday,iyear,(tbarrow(i,m,j)-
     1                   tfrez,thlqrow(i,m,j),thicrow(i,m,j),j=1,3),
     2                  tcn,rcanrow(i,m),scanrow(i,m),tsn,zsn,' TILE ',m
              write(66,6601) ihour,imin,iday,iyear,(tbarrow(i,m,j)-
     1                   tfrez,thlqrow(i,m,j),thicrow(i,m,j),j=4,10),
     2                   (gflxrow(i,m,j),j=1,10),
     3                   ' TILE ',m
          else
              write(65,6500) ihour,imin,iday,iyear,(tbarrow(i,m,j)-
     1                   tfrez,thlqrow(i,m,j),thicrow(i,m,j),j=1,3),
     2                  tcn,rcanrow(i,m),scanrow(i,m),tsn,zsn,' TILE ',m
C===================== CTEM =====================================/
          ENDIF
C
          WRITE(67,6700) IHOUR,IMIN,IDAY,IYEAR,                    
     1                   TROFROW(I,M),TROOROW(I,M),TROSROW(I,M),
     2                   TROBROW(I,M),ROFROW(I,M),ROFOROW(I,M),
     3                   ROFSROW(I,M),ROFBROW(I,M),
     4                   FCS(M),FGS(M),FC(M),FG(M),' TILE ',M
          WRITE(68,6800) IHOUR,IMIN,IDAY,IYEAR,                    
     1                   FSGVROW(I,M),FSGSROW(I,M),FSGGROW(I,M),
     2                   FLGVROW(I,M),FLGSROW(I,M),FLGGROW(I,M),
     3                   HFSCROW(I,M),HFSSROW(I,M),HFSGROW(I,M),
     4                   HEVCROW(I,M),HEVSROW(I,M),HEVGROW(I,M),
     5                   HMFCROW(I,M),HMFNROW(I,M),
     6                   (HMFGROW(I,M,J),J=1,3),
     7                   HTCCROW(I,M),HTCSROW(I,M),
     8                   (HTCROW(I,M,J),J=1,3),' TILE ',M
          WRITE(69,6900) IHOUR,IMIN,IDAY,IYEAR,                   
     1                   PCFCROW(I,M),PCLCROW(I,M),PCPNROW(I,M),
     2                   PCPGROW(I,M),QFCFROW(I,M),QFCLROW(I,M),
     3                   QFNROW(I,M),QFGROW(I,M),(QFCROW(I,M,J),J=1,3),
     4                   ROFCROW(I,M),ROFNROW(I,M),ROFOROW(I,M),
     5                   ROFROW(I,M),WTRCROW(I,M),WTRSROW(I,M),
     6                   WTRGROW(I,M),' TILE ',M
C===================== CTEM =====================================\
C
          endif ! if ((iyd.ge.jhhst).and.(iyd.le.jhhend))
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
              iyd=iyear*1000+iday                  
              if ((iyd.ge.jhhst).and.(iyd.le.jhhend)) then 
              write(71,7200)ihour,imin,iday,(anvegrow(i,m,j),j=1,icc),
     1                    (rmlvegrow(i,m,j),j=1,icc),' TILE ',m
              endif
          endif  ! if(ctem_on) 
c
7200      format(1x,i2,1x,i2,i5,9f11.3,9f11.3,2(a6,i2))
c
          if (ctem_on) then
            do j = 1,icc
              anvegrow_g(i,j)=anvegrow_g(i,j)+anvegrow(i,m,j)
     1                                        *farerow(i,m)
              rmlvegrow_g(i,j)=rmlvegrow_g(i,j)+rmlvegrow(i,m,j)
     1                                         *farerow(i,m)
            enddo
          endif   ! ctem_on
c
          fsstar_g    =fsstar_g + fsstar*farerow(i,m)
          flstar_g    =flstar_g + flstar*farerow(i,m)
          qh_g        =qh_g     + qh*farerow(i,m)
          qe_g        =qe_g     + qe*farerow(i,m)
          snomlt_g    =snomlt_g + snomlt*farerow(i,m)
          beg_g       =beg_g    + beg*farerow(i,m)
          gtout_g     =gtout_g  + gtout*farerow(i,m)
          tcn_g=tcn_g + tcn*farerow(i,m)
          tsn_g=tsn_g + tsn*farerow(i,m)
          zsn_g=zsn_g + zsn*farerow(i,m)
          altot_g     =altot_g   + altot*farerow(i,m)
          tpn_g       =tpn_g       + tpn*farerow(i,m)
c
          do j=1,ignd
            tbarrow_g(i,j)=tbarrow_g(i,j) + tbarrow(i,m,j)*farerow(i,m)
            thlqrow_g(i,j)=thlqrow_g(i,j) + thlqrow(i,m,j)*farerow(i,m)
            thicrow_g(i,j)=thicrow_g(i,j) + thicrow(i,m,j)*farerow(i,m)
            gflxrow_g(i,j)=gflxrow_g(i,j) + gflxrow(i,m,j)*farerow(i,m)
            hmfgrow_g(i,j)=hmfgrow_g(i,j) + hmfgrow(i,m,j)*farerow(i,m)
            htcrow_g(i,j)=htcrow_g(i,j) + htcrow(i,m,j)*farerow(i,m)
            qfcrow_g(i,j)=qfcrow_g(i,j) + qfcrow(i,m,j)*farerow(i,m)
          enddo
c
          zpndrow_g(i)=zpndrow_g(i) + zpndrow(i,m)*farerow(i,m) 
          rhosrow_g(i)=rhosrow_g(i) + rhosrow(i,m)*farerow(i,m)
          wsnorow_g(i)=wsnorow_g(i) + wsnorow(i,m)*farerow(i,m)
          rcanrow_g(i)=rcanrow_g(i) + rcanrow(i,m)*farerow(i,m)
          scanrow_g(i)=scanrow_g(i) + scanrow(i,m)*farerow(i,m)
          trofrow_g(i)=trofrow_g(i) + trofrow(i,m)*farerow(i,m)
          troorow_g(i)=troorow_g(i) + troorow(i,m)*farerow(i,m)
          trosrow_g(i)=trosrow_g(i) + trosrow(i,m)*farerow(i,m)
          trobrow_g(i)=trobrow_g(i) + trobrow(i,m)*farerow(i,m)
          roforow_g(i)=roforow_g(i) + roforow(i,m)*farerow(i,m)
          rofsrow_g(i)=rofsrow_g(i) + rofsrow(i,m)*farerow(i,m)
          rofbrow_g(i)=rofbrow_g(i) + rofbrow(i,m)*farerow(i,m)
          fsgvrow_g(i)=fsgvrow_g(i) + fsgvrow(i,m)*farerow(i,m)
          fsgsrow_g(i)=fsgsrow_g(i) + fsgsrow(i,m)*farerow(i,m)
          fsggrow_g(i)=fsggrow_g(i) + fsggrow(i,m)*farerow(i,m)
          flgvrow_g(i)=flgvrow_g(i) + flgvrow(i,m)*farerow(i,m)
          flgsrow_g(i)=flgsrow_g(i) + flgsrow(i,m)*farerow(i,m)
          flggrow_g(i)=flggrow_g(i) + flggrow(i,m)*farerow(i,m)
          hfscrow_g(i)=hfscrow_g(i) + hfscrow(i,m)*farerow(i,m)
          hfssrow_g(i)=hfssrow_g(i) + hfssrow(i,m)*farerow(i,m)
          hfsgrow_g(i)=hfsgrow_g(i) + hfsgrow(i,m)*farerow(i,m)
          hevcrow_g(i)=hevcrow_g(i) + hevcrow(i,m)*farerow(i,m)
          hevsrow_g(i)=hevsrow_g(i) + hevsrow(i,m)*farerow(i,m)
          hevgrow_g(i)=hevgrow_g(i) + hevgrow(i,m)*farerow(i,m)
          hmfcrow_g(i)=hmfcrow_g(i) + hmfcrow(i,m)*farerow(i,m)
          hmfnrow_g(i)=hmfnrow_g(i) + hmfnrow(i,m)*farerow(i,m)               
          htccrow_g(i)=htccrow_g(i) + htccrow(i,m)*farerow(i,m)
          htcsrow_g(i)=htcsrow_g(i) + htcsrow(i,m)*farerow(i,m)
          pcfcrow_g(i)=pcfcrow_g(i) + pcfcrow(i,m)*farerow(i,m)
          pclcrow_g(i)=pclcrow_g(i) + pclcrow(i,m)*farerow(i,m)
          pcpnrow_g(i)=pcpnrow_g(i) + pcpnrow(i,m)*farerow(i,m)
          pcpgrow_g(i)=pcpgrow_g(i) + pcpgrow(i,m)*farerow(i,m)
          qfcfrow_g(i)=qfcfrow_g(i) + qfcfrow(i,m)*farerow(i,m)
          qfclrow_g(i)=qfclrow_g(i) + qfclrow(i,m)*farerow(i,m)
          rofcrow_g(i)=rofcrow_g(i) + rofcrow(i,m)*farerow(i,m)
          rofnrow_g(i)=rofnrow_g(i) + rofnrow(i,m)*farerow(i,m)
          wtrcrow_g(i)=wtrcrow_g(i) + wtrcrow(i,m)*farerow(i,m)
          wtrsrow_g(i)=wtrsrow_g(i) + wtrsrow(i,m)*farerow(i,m)
          wtrgrow_g(i)=wtrgrow_g(i) + wtrgrow(i,m)*farerow(i,m)
          qfnrow_g(i) =qfnrow_g(i) + qfnrow(i,m)*farerow(i,m)
          qfgrow_g(i) =qfgrow_g(i) + qfgrow(i,m)*farerow(i,m)
          rofrow_g(i) =rofrow_g(i) + rofrow(i,m)*farerow(i,m)
          snorow_g(i) =snorow_g(i) + snorow(i,m)*farerow(i,m)
C
C======================== CTEM =====================================/
425    CONTINUE
C===================== CTEM =====================================\
C      WRITE CTEM OUTPUT FILES
C
       IF (CTEM_ON) THEN
           IF ((IYD.GE.JHHST).AND.(IYD.LE.JHHEND)) THEN  
           WRITE(711,7200)IHOUR,IMIN,IDAY,(ANVEGROW_G(I,J),J=1,ICC),
     1                 (RMLVEGROW_G(I,J),J=1,ICC)
           ENDIF  
       ENDIF !CTEM_ON

       IF ((IYD.GE.JHHST).AND.(IYD.LE.JHHEND)) THEN 
         WRITE(641,6400) IHOUR,IMIN,IDAY,IYEAR,FSSTAR_G,FLSTAR_G,QH_G,
     1      QE_G,SNOMLT_G,BEG_G,GTOUT_G,SNOROW_G(I),RHOSROW_G(I),
     2                   WSNOROW_G(I),ALTOT_G,ROFROW_G(I),
     3                   TPN_G,ZPNDROW_G(I)
         WRITE(651,6500) IHOUR,IMIN,IDAY,IYEAR,(TBARROW_G(I,J)-
     1                   TFREZ,THLQROW_G(I,J),THICROW_G(I,J),J=1,3),
     2                   TCN_G,RCANROW_G(I),SCANROW_G(I),TSN_G,ZSN_G
C
         IF(IGND.GT.3) THEN
          WRITE(661,6601) IHOUR,IMIN,IDAY,IYEAR,(TBARROW_G(I,J)-
     1                   TFREZ,THLQROW_G(I,J),THICROW_G(I,J),J=4,10),
     2                   (GFLXROW_G(I,J),J=1,10)
         ELSE
          WRITE(661,6600) IHOUR,IMIN,IDAY,FSDOWN,FDLGRD(I),PREGRD(I),
     1                   TAGRD(I)-TFREZ,UVGRD(I),PRESGRD(I),QAGRD(I)
         ENDIF
C
         WRITE(671,6700) IHOUR,IMIN,IDAY,IYEAR,                    
     &                   TROFROW_G(I),TROOROW_G(I),TROSROW_G(I),
     1                   TROBROW_G(I),ROFROW_G(I),ROFOROW_G(I),
     2                   ROFSROW_G(I),ROFBROW_G(I),
     3                   FCS_G(I),FGS_G(I),FC_G(I),FG_G(I)
         WRITE(681,6800) IHOUR,IMIN,IDAY,IYEAR,                   
     &                   FSGVROW_G(I),FSGSROW_G(I),FSGGROW_G(I),
     1                   FLGVROW_G(I),FLGSROW_G(I),FLGGROW_G(I),
     2                   HFSCROW_G(I),HFSSROW_G(I),HFSGROW_G(I),
     3                   HEVCROW_G(I),HEVSROW_G(I),HEVGROW_G(I),
     4                   HMFCROW_G(I),HMFNROW_G(I),
     5                   (HMFGROW_G(I,J),J=1,3),
     6                   HTCCROW_G(I),HTCSROW_G(I),
     7                   (HTCROW_G(I,J),J=1,3)
         WRITE(691,6900) IHOUR,IMIN,IDAY,IYEAR,                   
     &                   PCFCROW_G(I),PCLCROW_G(I),PCPNROW_G(I),
     1                   PCPGROW_G(I),QFCFROW_G(I),QFCLROW_G(I),
     2                   QFNROW_G(I),QFGROW_G(I),(QFCROW_G(I,J),J=1,3),
     3                   ROFCROW_G(I),ROFNROW_G(I),ROFOROW_G(I),
     4                   ROFROW_G(I),WTRCROW_G(I),WTRSROW_G(I),
     5                   WTRGROW_G(I)
C
       ENDIF ! IF ((IYD.GE.JHHST).AND.(IYD.LE.JHHEND))
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
          CDHGRD(I)=0.
          CDMGRD(I)=0.
          HFSGRD(I)=0.
          TFXGRD(I)=0.
          QEVPGRD(I)=0.
          QFSGRD(I)=0.
          QFXGRD(I)=0.
          PETGRD(I)=0.
          GAGRD(I)=0.
          EFGRD(I)=0.
          GTGRD(I)=0.
          QGGRD(I)=0.
          TSFGRD(I)=0.
          ALVSGRD(I)=0.
          ALIRGRD(I)=0.
          SFCTGRD(I)=0.
          SFCUGRD(I)=0.
          SFCVGRD(I)=0.
          SFCQGRD(I)=0.
          FSNOGRD(I)=0.
          FSGVGRD(I)=0.
          FSGSGRD(I)=0.
          FSGGGRD(I)=0.
          FLGVGRD(I)=0.
          FLGSGRD(I)=0.
          FLGGGRD(I)=0.
          HFSCGRD(I)=0.
          HFSSGRD(I)=0.
          HFSGGRD(I)=0.
          HEVCGRD(I)=0.
          HEVSGRD(I)=0.
          HEVGGRD(I)=0.
          HMFCGRD(I)=0.
          HMFNGRD(I)=0.
          HTCCGRD(I)=0.
          HTCSGRD(I)=0.
          PCFCGRD(I)=0.
          PCLCGRD(I)=0.
          PCPNGRD(I)=0.
          PCPGGRD(I)=0.
          QFGGRD(I)=0.
          QFNGRD(I)=0.
          QFCLGRD(I)=0.
          QFCFGRD(I)=0.
          ROFGRD(I)=0.
          ROFOGRD(I)=0.
          ROFSGRD(I)=0.
          ROFBGRD(I)=0.
          ROFCGRD(I)=0.
          ROFNGRD(I)=0.
          ROVGGRD(I)=0.
          WTRCGRD(I)=0.
          WTRSGRD(I)=0.
          WTRGGRD(I)=0.
          DRGRD(I)=0.
          WTABGRD(I)=0.
          ILMOGRD(I)=0.
          UEGRD(I)=0.
          HBLGRD(I)=0.
          DO 500 J=1,IGND
              HMFGGRD(I,J)=0.
              HTCGRD(I,J)=0.
              QFCGRD(I,J)=0.
              GFLXGRD(I,J)=0.
500       CONTINUE
525   CONTINUE
C
      DO 600 I=1,NLTEST
      DO 575 M=1,NMTEST
          CDHGRD(I)=CDHGRD(I)+CDHROW(I,M)*FAREROW(I,M)
          CDMGRD(I)=CDMGRD(I)+CDMROW(I,M)*FAREROW(I,M)
          HFSGRD(I)=HFSGRD(I)+HFSROW(I,M)*FAREROW(I,M)
          TFXGRD(I)=TFXGRD(I)+TFXROW(I,M)*FAREROW(I,M)
          QEVPGRD(I)=QEVPGRD(I)+QEVPROW(I,M)*FAREROW(I,M)
          QFSGRD(I)=QFSGRD(I)+QFSROW(I,M)*FAREROW(I,M)
          QFXGRD(I)=QFXGRD(I)+QFXROW(I,M)*FAREROW(I,M)
          PETGRD(I)=PETGRD(I)+PETROW(I,M)*FAREROW(I,M)
          GAGRD(I)=GAGRD(I)+GAROW(I,M)*FAREROW(I,M)
          EFGRD(I)=EFGRD(I)+EFROW(I,M)*FAREROW(I,M)
          GTGRD(I)=GTGRD(I)+GTROW(I,M)*FAREROW(I,M)
          QGGRD(I)=QGGRD(I)+QGROW(I,M)*FAREROW(I,M)
          TSFGRD(I)=TSFGRD(I)+TSFROW(I,M)*FAREROW(I,M)
          ALVSGRD(I)=ALVSGRD(I)+ALVSROW(I,M)*FAREROW(I,M)
          ALIRGRD(I)=ALIRGRD(I)+ALIRROW(I,M)*FAREROW(I,M)
          SFCTGRD(I)=SFCTGRD(I)+SFCTROW(I,M)*FAREROW(I,M)
          SFCUGRD(I)=SFCUGRD(I)+SFCUROW(I,M)*FAREROW(I,M)
          SFCVGRD(I)=SFCVGRD(I)+SFCVROW(I,M)*FAREROW(I,M)
          SFCQGRD(I)=SFCQGRD(I)+SFCQROW(I,M)*FAREROW(I,M)
          FSNOGRD(I)=FSNOGRD(I)+FSNOROW(I,M)*FAREROW(I,M)
          FSGVGRD(I)=FSGVGRD(I)+FSGVROW(I,M)*FAREROW(I,M)
          FSGSGRD(I)=FSGSGRD(I)+FSGSROW(I,M)*FAREROW(I,M)
          FSGGGRD(I)=FSGGGRD(I)+FSGGROW(I,M)*FAREROW(I,M)
          FLGVGRD(I)=FLGVGRD(I)+FLGVROW(I,M)*FAREROW(I,M)
          FLGSGRD(I)=FLGSGRD(I)+FLGSROW(I,M)*FAREROW(I,M)
          FLGGGRD(I)=FLGGGRD(I)+FLGGROW(I,M)*FAREROW(I,M)
          HFSCGRD(I)=HFSCGRD(I)+HFSCROW(I,M)*FAREROW(I,M)
          HFSSGRD(I)=HFSSGRD(I)+HFSSROW(I,M)*FAREROW(I,M)
          HFSGGRD(I)=HFSGGRD(I)+HFSGROW(I,M)*FAREROW(I,M)
          HEVCGRD(I)=HEVCGRD(I)+HEVCROW(I,M)*FAREROW(I,M)
          HEVSGRD(I)=HEVSGRD(I)+HEVSROW(I,M)*FAREROW(I,M)
          HEVGGRD(I)=HEVGGRD(I)+HEVGROW(I,M)*FAREROW(I,M)
          HMFCGRD(I)=HMFCGRD(I)+HMFCROW(I,M)*FAREROW(I,M)
          HMFNGRD(I)=HMFNGRD(I)+HMFNROW(I,M)*FAREROW(I,M)
          HTCCGRD(I)=HTCCGRD(I)+HTCCROW(I,M)*FAREROW(I,M)
          HTCSGRD(I)=HTCSGRD(I)+HTCSROW(I,M)*FAREROW(I,M)
          PCFCGRD(I)=PCFCGRD(I)+PCFCROW(I,M)*FAREROW(I,M)
          PCLCGRD(I)=PCLCGRD(I)+PCLCROW(I,M)*FAREROW(I,M)
          PCPNGRD(I)=PCPNGRD(I)+PCPNROW(I,M)*FAREROW(I,M)
          PCPGGRD(I)=PCPGGRD(I)+PCPGROW(I,M)*FAREROW(I,M)
          QFGGRD(I)=QFGGRD(I)+QFGROW(I,M)*FAREROW(I,M)
          QFNGRD(I)=QFNGRD(I)+QFNROW(I,M)*FAREROW(I,M)
          QFCLGRD(I)=QFCLGRD(I)+QFCLROW(I,M)*FAREROW(I,M)
          QFCFGRD(I)=QFCFGRD(I)+QFCFROW(I,M)*FAREROW(I,M)
          ROFGRD(I)=ROFGRD(I)+ROFROW(I,M)*FAREROW(I,M)
          ROFOGRD(I)=ROFOGRD(I)+ROFOROW(I,M)*FAREROW(I,M)
          ROFSGRD(I)=ROFSGRD(I)+ROFSROW(I,M)*FAREROW(I,M)
          ROFBGRD(I)=ROFBGRD(I)+ROFBROW(I,M)*FAREROW(I,M)
          ROFCGRD(I)=ROFCGRD(I)+ROFCROW(I,M)*FAREROW(I,M)
          ROFNGRD(I)=ROFNGRD(I)+ROFNROW(I,M)*FAREROW(I,M)
          ROVGGRD(I)=ROVGGRD(I)+ROVGROW(I,M)*FAREROW(I,M)
          WTRCGRD(I)=WTRCGRD(I)+WTRCROW(I,M)*FAREROW(I,M)
          WTRSGRD(I)=WTRSGRD(I)+WTRSROW(I,M)*FAREROW(I,M)
          WTRGGRD(I)=WTRGGRD(I)+WTRGROW(I,M)*FAREROW(I,M)
          DRGRD(I)=DRGRD(I)+DRROW(I,M)*FAREROW(I,M)
          WTABGRD(I)=WTABGRD(I)+WTABROW(I,M)*FAREROW(I,M)
          ILMOGRD(I)=ILMOGRD(I)+ILMOROW(I,M)*FAREROW(I,M)
          UEGRD(I)=UEGRD(I)+UEROW(I,M)*FAREROW(I,M)
          HBLGRD(I)=HBLGRD(I)+HBLROW(I,M)*FAREROW(I,M)
          DO 550 J=1,IGND
              HMFGGRD(I,J)=HMFGGRD(I,J)+HMFGROW(I,M,J)*FAREROW(I,M)        
              HTCGRD(I,J)=HTCGRD(I,J)+HTCROW(I,M,J)*FAREROW(I,M)        
              QFCGRD(I,J)=QFCGRD(I,J)+QFCROW(I,M,J)*FAREROW(I,M)        
              GFLXGRD(I,J)=GFLXGRD(I,J)+GFLXROW(I,M,J)*FAREROW(I,M)        
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
          PREACC(I)=PREACC(I)+PREGRD(I)*FAREROW(I,M)*DELT
          GTACC(I)=GTACC(I)+GTROW(I,M)*FAREROW(I,M)
          QEVPACC(I)=QEVPACC(I)+QEVPROW(I,M)*FAREROW(I,M)
          EVAPACC(I)=EVAPACC(I)+QFSROW(I,M)*FAREROW(I,M)*DELT
          HFSACC(I)=HFSACC(I)+HFSROW(I,M)*FAREROW(I,M)
          HMFNACC(I)=HMFNACC(I)+HMFNROW(I,M)*FAREROW(I,M)
          ROFACC(I)=ROFACC(I)+ROFROW(I,M)*FAREROW(I,M)*DELT
          OVRACC(I)=OVRACC(I)+ROFOROW(I,M)*FAREROW(I,M)*DELT
          WTBLACC(I)=WTBLACC(I)+WTABROW(I,M)*FAREROW(I,M)
          DO 625 J=1,IGND
              TBARACC(I,J)=TBARACC(I,J)+TBARROW(I,M,J)*FAREROW(I,M)
              THLQACC(I,J)=THLQACC(I,J)+THLQROW(I,M,J)*FAREROW(I,M)
              THICACC(I,J)=THICACC(I,J)+THICROW(I,M,J)*FAREROW(I,M)
              THALACC(I,J)=THALACC(I,J)+(THLQROW(I,M,J)+THICROW(I,M,J))
     1                    *FAREROW(I,M)
625       CONTINUE
          ALVSACC(I)=ALVSACC(I)+ALVSROW(I,M)*FAREROW(I,M)*FSVHGRD(I)
          ALIRACC(I)=ALIRACC(I)+ALIRROW(I,M)*FAREROW(I,M)*FSIHGRD(I)
          IF(SNOROW(I,M).GT.0.0) THEN
              RHOSACC(I)=RHOSACC(I)+RHOSROW(I,M)*FAREROW(I,M)
              TSNOACC(I)=TSNOACC(I)+TSNOROW(I,M)*FAREROW(I,M)
              WSNOACC(I)=WSNOACC(I)+WSNOROW(I,M)*FAREROW(I,M)
              SNOARE(I)=SNOARE(I)+FAREROW(I,M)
          ENDIF
          IF(TCANROW(I,M).GT.0.5) THEN
              TCANACC(I)=TCANACC(I)+TCANROW(I,M)*FAREROW(I,M)
              CANARE(I)=CANARE(I)+FAREROW(I,M)
          ENDIF
          SNOACC(I)=SNOACC(I)+SNOROW(I,M)*FAREROW(I,M)
          RCANACC(I)=RCANACC(I)+RCANROW(I,M)*FAREROW(I,M)
          SCANACC(I)=SCANACC(I)+SCANROW(I,M)*FAREROW(I,M)
          GROACC(I)=GROACC(I)+GROROW(I,M)*FAREROW(I,M)
          FSINACC(I)=FSINACC(I)+FSDOWN*FAREROW(I,M)
          FLINACC(I)=FLINACC(I)+FDLGRD(I)*FAREROW(I,M)
          FLUTACC(I)=FLUTACC(I)+SBC*GTROW(I,M)**4*FAREROW(I,M)
          TAACC(I)=TAACC(I)+TAGRD(I)*FAREROW(I,M)
          UVACC(I)=UVACC(I)+UVGRD(I)*FAREROW(I,M)
          PRESACC(I)=PRESACC(I)+PRESGRD(I)*FAREROW(I,M)
          QAACC(I)=QAACC(I)+QAGRD(I)*FAREROW(I,M)
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
             IYD=IYEAR*1000+IDAY                         
             IF ((IYD.GE.JDST).AND.(IYD.LE.JDEND)) THEN  
              WRITE(61,6100) IDAY,IYEAR,FSSTAR,FLSTAR,QH,QE,SNOMLT,
     1                       BEG,GTOUT,SNOACC(I),RHOSACC(I),
     2                       WSNOACC(I),ALTOT,ROFACC(I)
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
      DO 676 I=1,NLTEST
      DO 658 M=1,NMTEST
          PREACC_M(I,M)=PREACC_M(I,M)+PREGRD(I)*DELT
          GTACC_M(I,M)=GTACC_M(I,M)+GTROW(I,M)
          QEVPACC_M(I,M)=QEVPACC_M(I,M)+QEVPROW(I,M)
          EVAPACC_M(I,M)=EVAPACC_M(I,M)+QFSROW(I,M)*DELT
          HFSACC_M(I,M)=HFSACC_M(I,M)+HFSROW(I,M)
          HMFNACC_M(I,M)=HMFNACC_M(I,M)+HMFNROW(I,M)
          ROFACC_M(I,M)=ROFACC_M(I,M)+ROFROW(I,M)*DELT
          OVRACC_M(I,M)=OVRACC_M(I,M)+ROFOROW(I,M)*DELT
          WTBLACC_M(I,M)=WTBLACC_M(I,M)+WTABROW(I,M)
          DO 626 J=1,IGND
              TBARACC_M(I,M,J)=TBARACC_M(I,M,J)+TBARROW(I,M,J)
              THLQACC_M(I,M,J)=THLQACC_M(I,M,J)+THLQROW(I,M,J)
              THICACC_M(I,M,J)=THICACC_M(I,M,J)+THICROW(I,M,J)
              THALACC_M(I,M,J)=THALACC_M(I,M,J)+(THLQROW(I,M,J)+
     1           THICROW(I,M,J))
626       CONTINUE
          ALVSACC_M(I,M)=ALVSACC_M(I,M)+ALVSROW(I,M)*FSVHGRD(I)
          ALIRACC_M(I,M)=ALIRACC_M(I,M)+ALIRROW(I,M)*FSIHGRD(I)
          IF(SNOROW(I,M).GT.0.0) THEN
              RHOSACC_M(I,M)=RHOSACC_M(I,M)+RHOSROW(I,M)
              TSNOACC_M(I,M)=TSNOACC_M(I,M)+TSNOROW(I,M)
              WSNOACC_M(I,M)=WSNOACC_M(I,M)+WSNOROW(I,M)
          ENDIF
          IF(TCANROW(I,M).GT.0.5) THEN
              TCANACC_M(I,M)=TCANACC_M(I,M)+TCANROW(I,M)
C              CANARE(I)=CANARE(I)+FAREROW(I,M)
          ENDIF
          SNOACC_M(I,M)=SNOACC_M(I,M)+SNOROW(I,M)
          RCANACC_M(I,M)=RCANACC_M(I,M)+RCANROW(I,M)
          SCANACC_M(I,M)=SCANACC_M(I,M)+SCANROW(I,M)
          GROACC_M(I,M)=GROACC_M(I,M)+GROROW(I,M)
          FSINACC_M(I,M)=FSINACC_M(I,M)+FSDOWN
          FLINACC_M(I,M)=FLINACC_M(I,M)+FDLGRD(I)
          FLUTACC_M(I,M)=FLUTACC_M(I,M)+SBC*GTROW(I,M)**4
          TAACC_M(I,M)=TAACC_M(I,M)+TAGRD(I)
          UVACC_M(I,M)=UVACC_M(I,M)+UVGRD(I)
          PRESACC_M(I,M)=PRESACC_M(I,M)+PRESGRD(I)
          QAACC_M(I,M)=QAACC_M(I,M)+QAGRD(I)
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
          RHOSACC_M(I,M)=RHOSACC_M(I,M)/REAL(NDAY)  
          TSNOACC_M(I,M)=TSNOACC_M(I,M)/REAL(NDAY)  
          WSNOACC_M(I,M)=WSNOACC_M(I,M)/REAL(NDAY) 
          TCANACC_M(I,M)=TCANACC_M(I,M)/REAL(NDAY)  
          SNOACC_M(I,M)=SNOACC_M(I,M)/REAL(NDAY)
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
          QEVPACC_M_SAVE(I,M)=QEVPACC_M(I,M)  
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
          IYD=IYEAR*1000+IDAY                        
          IF ((IYD.GE.JDST).AND.(IYD.LE.JDEND)) THEN
C
C         WRITE TO OUTPUT FILES
C
          WRITE(611,6100) IDAY,IYEAR,FSSTAR,FLSTAR,QH,QE,SNOMLT,
     1                    BEG,GTOUT,SNOACC_M(I,M),RHOSACC_M(I,M),
     2                    WSNOACC_M(I,M),ALTOT,ROFACC_M(I,M),' TILE ',M
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
           ENDIF ! IF ((IYD.GE.JDST).AND.(IYD.LE.JDEND))
C
C          INITIALIZTION FOR MOSAIC TILE AND GRID VARIABLES
C
           PREACC_M(I,M)=0.
           GTACC_M(I,M)=0.
           QEVPACC_M(I,M)=0.
           HFSACC_M(I,M)=0.
           HMFNACC_M(I,M)=0.
           ROFACC_M(I,M)=0.
           SNOACC_M(I,M)=0.
           OVRACC_M(I,M)=0.
           WTBLACC_M(I,M)=0.
           ALVSACC_M(I,M)=0.
           ALIRACC_M(I,M)=0.
           RHOSACC_M(I,M)=0.
           TSNOACC_M(I,M)=0.
           WSNOACC_M(I,M)=0.
           TCANACC_M(I,M)=0.
           RCANACC_M(I,M)=0.
           SCANACC_M(I,M)=0.
           GROACC_M(I,M)=0.
           FSINACC_M(I,M)=0.
           FLINACC_M(I,M)=0.
           TAACC_M(I,M)=0.
           UVACC_M(I,M)=0.
           PRESACC_M(I,M)=0.
           QAACC_M(I,M)=0.
           EVAPACC_M(I,M)=0.
           FLUTACC_M(I,M)=0.
C
           DO 759 J=1,IGND
             TBARACC_M(I,M,J)=0.
             THLQACC_M(I,M,J)=0.
             THICACC_M(I,M,J)=0.
             THALACC_M(I,M,J)=0.
759        CONTINUE
C
809   CONTINUE
808   CONTINUE

      ENDIF ! IF(NCOUNT.EQ.NDAY)
C
      ENDIF !  IF(.NOT.PARALLELRUN)
C
C=======================================================================
C
C     ACCUMULATE OUTPUT DATA FOR MONTHLY AVERAGED FIELDS FOR CLASS GRID-MEAN.
C     FOR BOTH PARALLEL MODE AND STAND ALONE MODE
C
      FSSTAR_MO   =0.0
      FLSTAR_MO   =0.0
      QH_MO       =0.0
      QE_MO       =0.0
      ALTOT_MO    =0.0
C
      DO 820 I=1,NLTEST
       DO 821 M=1,NMTEST
          ALVSACC_MO(I)=ALVSACC_MO(I)+ALVSROW(I,M)*FAREROW(I,M)
     1                  *FSVHGRD(I)
          ALIRACC_MO(I)=ALIRACC_MO(I)+ALIRROW(I,M)*FAREROW(I,M)
     1                  *FSIHGRD(I) 
          FLUTACC_MO(I)=FLUTACC_MO(I)+SBC*GTROW(I,M)**4*FAREROW(I,M)
          FSINACC_MO(I)=FSINACC_MO(I)+FSDOWN*FAREROW(I,M)
          FLINACC_MO(I)=FLINACC_MO(I)+FDLGRD(I)*FAREROW(I,M)
          HFSACC_MO(I) =HFSACC_MO(I)+HFSROW(I,M)*FAREROW(I,M)
          QEVPACC_MO(I)=QEVPACC_MO(I)+QEVPROW(I,M)*FAREROW(I,M)
          SNOACC_MO(I) =SNOACC_MO(I)+SNOROW(I,M)*FAREROW(I,M)
          TAACC_MO(I)=TAACC_MO(I)+TAGRD(I)*FAREROW(I,M)
C
          IF(SNOROW(I,M).GT.0.0) THEN
           WSNOACC_MO(I)=WSNOACC_MO(I)+WSNOROW(I,M)*FAREROW(I,M)
          ENDIF
C
          ROFACC_MO(I) =ROFACC_MO(I)+ROFROW(I,M)*FAREROW(I,M)*DELT
          PREACC_MO(I) =PREACC_MO(I)+PREGRD(I)*FAREROW(I,M)*DELT
          EVAPACC_MO(I)=EVAPACC_MO(I)+QFSROW(I,M)*FAREROW(I,M)*DELT
C
          DO 823 J=1,IGND
           TBARACC_MO(I,J)=TBARACC_MO(I,J)+TBARROW(I,M,J)*FAREROW(I,M)
           THLQACC_MO(I,J)=THLQACC_MO(I,J)+THLQROW(I,M,J)*FAREROW(I,M)
           THICACC_MO(I,J)=THICACC_MO(I,J)+THICROW(I,M,J)*FAREROW(I,M)
823       CONTINUE
C
821    CONTINUE
820   CONTINUE 
C
      DO NT=1,NMON
       IF(IDAY.EQ.monthend(NT+1).AND.NCOUNT.EQ.NDAY)THEN
        IMONTH=NT
        NDMONTH=(monthend(NT+1)-monthend(NT))*NDAY
C
        DO 824 I=1,NLTEST
         IF(FSINACC_MO(I).GT.0.0) THEN
          ALVSACC_MO(I)=ALVSACC_MO(I)/(FSINACC_MO(I)*0.5)
          ALIRACC_MO(I)=ALIRACC_MO(I)/(FSINACC_MO(I)*0.5)
         ELSE
          ALVSACC_MO(I)=0.0
          ALIRACC_MO(I)=0.0
         ENDIF
         FLUTACC_MO(I)=FLUTACC_MO(I)/REAL(NDMONTH)
         FSINACC_MO(I)=FSINACC_MO(I)/REAL(NDMONTH)
         FLINACC_MO(I)=FLINACC_MO(I)/REAL(NDMONTH)
         HFSACC_MO(I) =HFSACC_MO(I)/REAL(NDMONTH)
         QEVPACC_MO(I)=QEVPACC_MO(I)/REAL(NDMONTH)
         SNOACC_MO(I) =SNOACC_MO(I)/REAL(NDMONTH)
         WSNOACC_MO(I)=WSNOACC_MO(I)/REAL(NDMONTH)
         ROFACC_MO(I) =ROFACC_MO(I)
         PREACC_MO(I) =PREACC_MO(I)
         EVAPACC_MO(I)=EVAPACC_MO(I)
         TAACC_MO(I)=TAACC_MO(I)/REAL(NDMONTH)
         DO J=1,IGND
          TBARACC_MO(I,J)=TBARACC_MO(I,J)/REAL(NDMONTH)
          THLQACC_MO(I,J)=THLQACC_MO(I,J)/REAL(NDMONTH)
          THICACC_MO(I,J)=THICACC_MO(I,J)/REAL(NDMONTH)
         ENDDO
C
         ALTOT_MO=(ALVSACC_MO(I)+ALIRACC_MO(I))/2.0
         FSSTAR_MO=FSINACC_MO(I)*(1.-ALTOT_MO)
         FLSTAR_MO=FLINACC_MO(I)-FLUTACC_MO(I)
         QH_MO=HFSACC_MO(I)
         QE_MO=QEVPACC_MO(I)
C
         WRITE(81,8100)IMONTH,IYEAR,FSSTAR_MO,FLSTAR_MO,QH_MO,
     1                 QE_MO,SNOACC_MO(I),WSNOACC_MO(I),
     2                 ROFACC_MO(I),PREACC_MO(I),EVAPACC_MO(I),
     3                 TAACC_MO(I)-TFREZ
         IF (IGND.GT.3) THEN
          WRITE(82,8101)IMONTH,IYEAR,(TBARACC_MO(I,J)-TFREZ,
     1                  THLQACC_MO(I,J),THICACC_MO(I,J),J=1,5)
          WRITE(82,8101)IMONTH,IYEAR,(TBARACC_MO(I,J)-TFREZ,
     1                  THLQACC_MO(I,J),THICACC_MO(I,J),J=6,10)
         ELSE
          WRITE(82,8102)IMONTH,IYEAR,(TBARACC_MO(I,J)-TFREZ,
     1                  THLQACC_MO(I,J),THICACC_MO(I,J),J=1,3)
         ENDIF     
C
C ADD INITIALIZTION FOR MONTHLY ACCUMULATED ARRAYS
C
         ALVSACC_MO(I)=0.0
         ALIRACC_MO(I)=0.0  
         FLUTACC_MO(I)=0.0
         FSINACC_MO(I)=0.0
         FLINACC_MO(I)=0.0
         HFSACC_MO(I)=0.0
         QEVPACC_MO(I)=0.0
         SNOACC_MO(I)=0.0
         WSNOACC_MO(I)=0.0
         ROFACC_MO(I)=0.0
         PREACC_MO(I)=0.0
         EVAPACC_MO(I)=0.0
         TAACC_MO(I)=0.0
         DO 826 J=1,IGND
           TBARACC_MO(I,J)=0.
           THLQACC_MO(I,J)=0.
           THICACC_MO(I,J)=0.
826      CONTINUE   
C
824     CONTINUE ! I
C               
       ENDIF ! IF(IDAY.EQ.monthend(NT+1).AND.NCOUNT.EQ.NDAY)
      ENDDO ! NMON
C
8100  FORMAT(1X,I4,I5,5(F8.2,1X),F8.3,F12.4,3(E12.3,1X),2(A6,I2))
8101  FORMAT(1X,I4,I5,5(F7.2,1X,2F6.3,1X),2(A6,I2))
8102  FORMAT(1X,I4,I5,3(F8.2,1X,2F6.3,1X),2(A6,I2))
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
          ALVSACC_YR(I)=ALVSACC_YR(I)+ALVSROW(I,M)*FAREROW(I,M)
     1                  *FSVHGRD(I)
          ALIRACC_YR(I)=ALIRACC_YR(I)+ALIRROW(I,M)*FAREROW(I,M)
     1                  *FSIHGRD(I) 
          FLUTACC_YR(I)=FLUTACC_YR(I)+SBC*GTROW(I,M)**4*FAREROW(I,M)
          FSINACC_YR(I)=FSINACC_YR(I)+FSDOWN*FAREROW(I,M)
          FLINACC_YR(I)=FLINACC_YR(I)+FDLGRD(I)*FAREROW(I,M)
          HFSACC_YR(I) =HFSACC_YR(I)+HFSROW(I,M)*FAREROW(I,M)
          QEVPACC_YR(I)=QEVPACC_YR(I)+QEVPROW(I,M)*FAREROW(I,M)
          TAACC_YR(I)=TAACC_YR(I)+TAGRD(I)*FAREROW(I,M)
          ROFACC_YR(I) =ROFACC_YR(I)+ROFROW(I,M)*FAREROW(I,M)*DELT
          PREACC_YR(I) =PREACC_YR(I)+PREGRD(I)*FAREROW(I,M)*DELT
          EVAPACC_YR(I)=EVAPACC_YR(I)+QFSROW(I,M)*FAREROW(I,M)*DELT
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
         ALVSACC_YR(I)=0.0
         ALIRACC_YR(I)=0.0  
         FLUTACC_YR(I)=0.0
         FSINACC_YR(I)=0.0
         FLINACC_YR(I)=0.0
         HFSACC_YR(I)=0.0
         QEVPACC_YR(I)=0.0
         ROFACC_YR(I)=0.0
         PREACC_YR(I)=0.0
         EVAPACC_YR(I)=0.0
         TAACC_YR(I)=0.0
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
c
        do i=1,nltest
           do j=1,icc
             ifcancmx_g(i,j)=0   
           enddo
        enddo
c
        do i=1,nltest
           do m=1,nmtest
              ifcancmx_m(i,m)=0   !0=bare soil tile; 1=tile with vegetation
              leaflitr_m(i,m)=0.0
              tltrleaf_m(i,m)=0.0
              tltrstem_m(i,m)=0.0
              tltrroot_m(i,m)=0.0
              ailcg_m(i,m)=0.0
              ailcb_m(i,m)=0.0
              afrleaf_m(i,m)=0.0
              afrstem_m(i,m)=0.0
              afrroot_m(i,m)=0.0
              veghght_m(i,m)=0.0
              rootdpth_m(i,m)=0.0
              roottemp_m(i,m)=0.0
              slai_m(i,m)=0.0
              gleafmas_m(i,m) = 0.0
              bleafmas_m(i,m) = 0.0
              stemmass_m(i,m) = 0.0
              rootmass_m(i,m) = 0.0
              litrmass_m(i,m) = 0.0
              soilcmas_m(i,m) = 0.0

c
              do j=1,icc
                if (fcancmxrow(i,m,j) .gt.0.0) then
                ifcancmx_g(i,j)=1
                ifcancmx_m(i,m)=1
                endif 
              enddo
c
              do k=1,ignd
                rmatctem_m(i,m,k)=0.0
              enddo
c
           enddo ! m
        enddo  ! i
c
c       ---------------------------------------------------------
c
        do 851 i=1,nltest
          gpp_g(i) =0.0
          npp_g(i) =0.0
          nep_g(i) =0.0
          nbp_g(i) =0.0
          autores_g(i) =0.0
          hetrores_g(i)=0.0
          litres_g(i) =0.0
          socres_g(i) =0.0
          dstcemls_g(i)=0.0
          dstcemls3_g(i)=0.0
          litrfall_g(i)=0.0
          humiftrs_g(i)=0.0
          rml_g(i) =0.0
          rms_g(i) =0.0
          rmr_g(i) =0.0
          rg_g(i) =0.0
          vgbiomas_g(i) =0.0
          totcmass_g(i) =0.0
          gavglai_g(i) =0.0
          gavgltms_g(i) =0.0
          gavgscms_g(i) =0.0
          ailcg_g(i)=0.0
          ailcb_g(i)=0.0
          tcanoacc_out_g(i) =0.0
          burnfrac_g(i) =0.0
          probfire_g(i) =0.0
          lucemcom_g(i) =0.0
          lucltrin_g(i) =0.0
          lucsocin_g(i) =0.0
          emit_co2_g(i) =0.0
          emit_co_g(i)  =0.0
          emit_ch4_g(i) =0.0  
          emit_nmhc_g(i) =0.0
          emit_h2_g(i) =0.0
          emit_nox_g(i) =0.0
          emit_n2o_g(i) =0.0 
          emit_pm25_g(i) =0.0
          emit_tpm_g(i) =0.0 
          emit_tc_g(i) =0.0
          emit_oc_g(i) =0.0  
          emit_bc_g(i) =0.0
          bterm_g(i)   =0.0
          lterm_g(i)   =0.0
          mterm_g(i)   =0.0
          leaflitr_g(i)=0.0  
          tltrleaf_g(i)=0.0
          tltrstem_g(i)=0.0
          tltrroot_g(i)=0.0
          gleafmas_g(i)=0.0
          bleafmas_g(i)=0.0
          stemmass_g(i)=0.0
          rootmass_g(i)=0.0
          litrmass_g(i)=0.0
          soilcmas_g(i)=0.0
          veghght_g(i)=0.0
          rootdpth_g(i)=0.0
          roottemp_g(i)=0.0
          slai_g(i)=0.0
c                         !Rudra added CH4 realted variables on 03/12/2013
          CH4WET1_G(i) = 0.0
          CH4WET2_G(i) = 0.0
          WETFDYN_G(i) = 0.0
          CH4DYN1_G(i) = 0.0
          CH4DYN2_G(i) = 0.0

          do k=1,ignd
           rmatctem_g(i,k)=0.0
          enddo
c
          do j=1,icc        
            afrleaf_g(i,j)=0.0
            afrstem_g(i,j)=0.0
            afrroot_g(i,j)=0.0
          enddo
c
          do 852 m=1,nmtest
c
           do j=1,icc
             leaflitr_m(i,m)=leaflitr_m(i,m)+
     &                       leaflitrrow(i,m,j)*fcancmxrow(i,m,j)
             tltrleaf_m(i,m)=tltrleaf_m(i,m)+
     &                       tltrleafrow(i,m,j)*fcancmxrow(i,m,j)
             tltrstem_m(i,m)=tltrstem_m(i,m)+
     &                       tltrstemrow(i,m,j)*fcancmxrow(i,m,j)
             tltrroot_m(i,m)=tltrroot_m(i,m)+
     &                       tltrrootrow(i,m,j)*fcancmxrow(i,m,j)
             veghght_m(i,m)=veghght_m(i,m)+
     &                            veghghtrow(i,m,j)*fcancmxrow(i,m,j)
             rootdpth_m(i,m)=rootdpth_m(i,m)+
     &                            rootdpthrow(i,m,j)*fcancmxrow(i,m,j)
             roottemp_m(i,m)=roottemp_m(i,m)+
     &                            roottemprow(i,m,j)*fcancmxrow(i,m,j)
             slai_m(i,m)=slai_m(i,m)+slairow(i,m,j)*fcancmxrow(i,m,j)
c
             afrleaf_m(i,m)=afrleaf_m(i,m)+
     &                              afrleafrow(i,m,j)*fcancmxrow(i,m,j)
             afrstem_m(i,m)=afrstem_m(i,m)+
     &                              afrstemrow(i,m,j)*fcancmxrow(i,m,j)
             afrroot_m(i,m)=afrroot_m(i,m)+
     &                              afrrootrow(i,m,j)*fcancmxrow(i,m,j)
c
             ailcg_m(i,m)=ailcg_m(i,m)+ailcgrow(i,m,j)*fcancmxrow(i,m,j)
             ailcb_m(i,m)=ailcb_m(i,m)+ailcbrow(i,m,j)*fcancmxrow(i,m,j)

             gleafmas_m(i,m) = gleafmas_m(i,m) + gleafmasrow(i,m,j)
     &                                          *fcancmxrow(i,m,j)
             bleafmas_m(i,m) = bleafmas_m(i,m) + bleafmasrow(i,m,j)
     &                                          *fcancmxrow(i,m,j)
             stemmass_m(i,m) = stemmass_m(i,m) + stemmassrow(i,m,j)
     &                                          *fcancmxrow(i,m,j)
             rootmass_m(i,m) = rootmass_m(i,m) + rootmassrow(i,m,j)
     &                                          *fcancmxrow(i,m,j)
             litrmass_m(i,m) = litrmass_m(i,m) + litrmassrow(i,m,j)
     &                                          *fcancmxrow(i,m,j)
             soilcmas_m(i,m) = soilcmas_m(i,m) + soilcmasrow(i,m,j)
     &                                          *fcancmxrow(i,m,j)

c
             do k=1,ignd
                rmatctem_m(i,m,k)=rmatctem_m(i,m,k)+
     &                            rmatctemrow(i,m,j,k)*fcancmxrow(i,m,j)
             enddo
           enddo
c
           npprow(i,m)     =npprow(i,m)*1.0377 ! convert to gc/m2.day
           gpprow(i,m)     =gpprow(i,m)*1.0377 ! convert to gc/m2.day
           neprow(i,m)     =neprow(i,m)*1.0377 ! convert to gc/m2.day
           nbprow(i,m)     =nbprow(i,m)*1.0377 ! convert to gc/m2.day
           lucemcomrow(i,m)=lucemcomrow(i,m)*1.0377 ! convert to gc/m2.day
           lucltrinrow(i,m)=lucltrinrow(i,m)*1.0377 ! convert to gc/m2.day
           lucsocinrow(i,m)=lucsocinrow(i,m)*1.0377 ! convert to gc/m2.day
c
           hetroresrow(i,m)=hetroresrow(i,m)*1.0377 ! convert to gc/m2.day
           autoresrow(i,m) =autoresrow(i,m)*1.0377  ! convert to gc/m2.day
           litresrow(i,m)  =litresrow(i,m)*1.0377   ! convert to gc/m2.day
           socresrow(i,m)  =socresrow(i,m)*1.0377   ! convert to gc/m2.day
c
           CH4WET1ROW(i,m) = CH4WET1ROW(i,m)*1.0377 * 16.044 / 12. ! convert from umolCH4/m2/s to gCH4/m2.day 
           CH4WET2ROW(i,m) = CH4WET2ROW(i,m)*1.0377 * 16.044 / 12. ! convert from umolCH4/m2/s to gCH4/m2.day
           CH4DYN1ROW(i,m) = CH4DYN1ROW(i,m)*1.0377 * 16.044 / 12. ! convert from umolCH4/m2/s to gCH4/m2.day
           CH4DYN2ROW(i,m) = CH4DYN2ROW(i,m)*1.0377 * 16.044 / 12. ! convert from umolCH4/m2/s to gCH4/m2.day 
c
c          write daily ctem results
c
           if ((iyd.ge.jdst).and.(iyd.le.jdend)) then   
c
c             write grid-averaged fluxes of basic quantities to 
c             file *.CT01D_M
c
             if (mosaic) then
              write(72,8200)iday,iyear,gpprow(i,m),npprow(i,m),
     1                neprow(i,m),nbprow(i,m),autoresrow(i,m),
     2                hetroresrow(i,m),litresrow(i,m),socresrow(i,m),
     3                (dstcemlsrow(i,m)+dstcemls3row(i,m)),
     4               litrfallrow(i,m),humiftrsrow(i,m),' TILE ',m,'AVGE'
             end if

c             write breakdown of some of basic fluxes to file *.CT3 
c             and selected litter fluxes for selected pft

!              First for the bare fraction of the grid cell.
             hetroresvegrow(i,m,iccp1)=hetroresvegrow(i,m,iccp1)*1.0377 ! convert to gc/m2.day
             litresvegrow(i,m,iccp1)=litresvegrow(i,m,iccp1)*1.0377 ! convert to gc/m2.day
             soilcresvegrow(i,m,iccp1)=soilcresvegrow(i,m,iccp1)*1.0377 ! convert to gc/m2.day

c
              do 853 j=1,icc
c
                if (fcancmxrow(i,m,j) .gt.0.0) then
c
                 gppvegrow(i,m,j)=gppvegrow(i,m,j)*1.0377 ! convert to gc/m2.day
                 nppvegrow(i,m,j)=nppvegrow(i,m,j)*1.0377 ! convert to gc/m2.day
                 nepvegrow(i,m,j)=nepvegrow(i,m,j)*1.0377 ! convert to gc/m2.day
                 nbpvegrow(i,m,j)=nbpvegrow(i,m,j)*1.0377 ! convert to gc/m2.day
                 hetroresvegrow(i,m,j)=hetroresvegrow(i,m,j)*1.0377 ! convert to gc/m2.day
                 autoresvegrow(i,m,j)=autoresvegrow(i,m,j)*1.0377 ! convert to gc/m2.day
                 litresvegrow(i,m,j)=litresvegrow(i,m,j)*1.0377 ! convert to gc/m2.day
                 soilcresvegrow(i,m,j)=soilcresvegrow(i,m,j)*1.0377 ! convert to gc/m2.day

 
c                write to file .CT01D_M 
                 if (mosaic) then
                  write(72,8201)iday,iyear,gppvegrow(i,m,j),
     1            nppvegrow(i,m,j),nepvegrow(i,m,j),
     2            ' TILE ',m,'PFT',j
                 

c                write to file .CT02D_M 
                 write(73,8300)iday,iyear,rmlvegaccrow(i,m,j), 
     1           rmsvegrow(i,m,j),rmrvegrow(i,m,j),rgvegrow(i,m,j),
     2           leaflitrrow(i,m,j),tltrleafrow(i,m,j),
     3          tltrstemrow(i,m,j),tltrrootrow(i,m,j),' TILE ',m,'PFT',j
c
c
c                write grid-averaged pool sizes and component sizes for
c                seleced ctem pft to file *.CT03D_M
c
                 write(74,8401)iday,iyear,vgbiomas_vegrow(i,m,j),
     1               ailcgrow(i,m,j),gleafmasrow(i,m,j),
     3               bleafmasrow(i,m,j), stemmassrow(i,m,j),
     4               rootmassrow(i,m,j), litrmassrow(i,m,j), 
     5               soilcmasrow(i,m,j),' TILE ',m,'PFT',j
c
c                write lai, rmatctem, & structural attributes for selected 
c                pft to file *.CT04D_M
c
                 write(75,8500)iday,iyear, ailcgrow(i,m,j), 
     1                ailcbrow(i,m,j),(rmatctemrow(i,m,j,k),k=1,3),
     2                veghghtrow(i,m,j),rootdpthrow(i,m,j),
     3              roottemprow(i,m,j),slairow(i,m,j),' TILE ',m,'PFT',j
c
c                write allocation fractions for selected pft to 
c                file *.CT05D_M
c
                 write(76,8600)iday,iyear, afrleafrow(i,m,j), 
     1                afrstemrow(i,m,j),afrrootrow(i,m,j), 
     2                tcanoaccrow_out(i,m), lfstatusrow(i,m,j),
     3                ' TILE ',m,'PFT',j

c             write fire and luc results to file *.CT06D_M
c
              if (dofire .or. lnduseon) then   !FLAG FIX THIS
               write(78,8800)iday,iyear,
     1         emit_co2row(i,m,j),emit_corow(i,m,j),emit_ch4row(i,m,j),
     2         emit_nmhcrow(i,m,j),emit_h2row(i,m,j),emit_noxrow(i,m,j),
     3         emit_n2orow(i,m,j),emit_pm25row(i,m,j),
     4         emit_tpmrow(i,m,j),emit_tcrow(i,m,j),emit_ocrow(i,m,j),
     5         emit_bcrow(i,m,j),burnvegfrow(i,m,j),probfirerow(i,m), 
     6         lucemcomrow(i,m),lucltrinrow(i,m), lucsocinrow(i,m),
     7         grclarea(i),btermrow(i,m),ltermrow(i,m),mtermrow(i,m),
     8         ' TILE ',m,'PFT',j
               endif

              end if !mosaic

              endif  !if (fcancmxrow(i,m,j) .gt.0.0) then
c
853           continue
c
              if (mosaic) then
               if (ifcancmx_m(i,m) .gt. 0) then


c               write to file .CT02D_M 
                write(73,8300)iday,iyear,rmlrow(i,m),rmsrow(i,m),
     1          rmrrow(i,m),rgrow(i,m),leaflitr_m(i,m),tltrleaf_m(i,m),
     2          tltrstem_m(i,m),tltrroot_m(i,m),' TILE ',m,'AVGE'
c
c               write to file .CT03D_M 
                write(74,8402)iday,iyear,vgbiomasrow(i,m),
     1               gavglairow(i,m),gavgltmsrow(i,m),
     2               gavgscmsrow(i,m),gleafmasrow(i,m,j),
     3               bleafmasrow(i,m,j), stemmassrow(i,m,j),
     4               rootmassrow(i,m,j), litrmassrow(i,m,j), 
     5               soilcmasrow(i,m,j),' TILE ',m, 'AVGE'
c
c               write to file .CT04D_M
                write(75,8500)iday,iyear,ailcg_m(i,m),
     1                ailcb_m(i,m),(rmatctem_m(i,m,k),k=1,3),
     2                veghght_m(i,m),rootdpth_m(i,m),
     3                roottemp_m(i,m),slai_m(i,m),' TILE ',m, 'AVGE'
c
c               write to file .CT05D_M
                write(76,8601)iday,iyear, afrleaf_m(i,m), 
     1                afrstem_m(i,m),afrroot_m(i,m), 
     2                tcanoaccrow_out(i,m), 
     3                ' TILE ',m,'AVGE'

               end if !if (ifcancmx_m(i,m) .gt.0.0) then
              endif !mosaic
c
           endif ! if ((iyd.ge.jdst).and.(iyd.le.jdend))
c
8200       format(1x,i4,i5,11f10.5,2(a6,i2))
8201       format(1x,i4,i5,3f10.5,80x,2(a6,i2))
8300       format(1x,i4,i5,8f10.5,2(a6,i2))
8301       format(1x,i4,i5,4f10.5,40x,2(a6,i2)) 
8400       format(1x,i4,i5,11f10.5,2(a6,i2))
8401       format(1x,i4,i5,2f10.5,6f10.5,2(a6,i2))
8402       format(1x,i4,i5,10f10.5,2(a6,i2))
!                   8402       format(1x,i4,i5,2f10.5,40x,2f10.5,2(a6,i2))
8500       format(1x,i4,i5,9f10.5,2(a6,i2))
8600       format(1x,i4,i5,4f10.5,i8,2(a6,i2))
8601       format(1x,i4,i5,4f10.5,8x,2(a6,i2))   
8800       format(1x,i4,i5,20f11.4,2x,f9.2,2(a6,i2))
8810       format(1x,i4,i5,5f11.4,2(a6,i2))
c
c          Calculation of grid averaged variables
c
           gpp_g(i) =gpp_g(i) + gpprow(i,m)*farerow(i,m)
           npp_g(i) =npp_g(i) + npprow(i,m)*farerow(i,m)
           nep_g(i) =nep_g(i) + neprow(i,m)*farerow(i,m)
           nbp_g(i) =nbp_g(i) + nbprow(i,m)*farerow(i,m)
           autores_g(i) =autores_g(i) +autoresrow(i,m)*farerow(i,m)
           hetrores_g(i)=hetrores_g(i)+hetroresrow(i,m)*farerow(i,m)
           litres_g(i) =litres_g(i) + litresrow(i,m)*farerow(i,m)
           socres_g(i) =socres_g(i) + socresrow(i,m)*farerow(i,m)
           dstcemls_g(i)=dstcemls_g(i)+dstcemlsrow(i,m)*farerow(i,m)
           dstcemls3_g(i)=dstcemls3_g(i)
     &                      +dstcemls3row(i,m)*farerow(i,m)

           litrfall_g(i)=litrfall_g(i)+litrfallrow(i,m)*farerow(i,m)
           humiftrs_g(i)=humiftrs_g(i)+humiftrsrow(i,m)*farerow(i,m)
           rml_g(i) =rml_g(i) + rmlrow(i,m)*farerow(i,m)
           rms_g(i) =rms_g(i) + rmsrow(i,m)*farerow(i,m)
           rmr_g(i) =rmr_g(i) + rmrrow(i,m)*farerow(i,m)
           rg_g(i) =rg_g(i) + rgrow(i,m)*farerow(i,m)
           leaflitr_g(i) = leaflitr_g(i) + leaflitr_m(i,m)*farerow(i,m)
           tltrleaf_g(i) = tltrleaf_g(i) + tltrleaf_m(i,m)*farerow(i,m)
           tltrstem_g(i) = tltrstem_g(i) + tltrstem_m(i,m)*farerow(i,m)
           tltrroot_g(i) = tltrroot_g(i) + tltrroot_m(i,m)*farerow(i,m)
           vgbiomas_g(i) =vgbiomas_g(i) + vgbiomasrow(i,m)*farerow(i,m)
           gavglai_g(i) =gavglai_g(i) + gavglairow(i,m)*farerow(i,m)
           gavgltms_g(i) =gavgltms_g(i) + gavgltmsrow(i,m)*farerow(i,m)
           gavgscms_g(i) =gavgscms_g(i) + gavgscmsrow(i,m)*farerow(i,m)
           tcanoacc_out_g(i) =tcanoacc_out_g(i)+
     1                        tcanoaccrow_out(i,m)*farerow(i,m)
           totcmass_g(i) =vgbiomas_g(i) + gavgltms_g(i) + gavgscms_g(i)
           gleafmas_g(i) = gleafmas_g(i) + gleafmas_m(i,m)*farerow(i,m)
           bleafmas_g(i) = bleafmas_g(i) + bleafmas_m(i,m)*farerow(i,m)
           stemmass_g(i) = stemmass_g(i) + stemmass_m(i,m)*farerow(i,m)
           rootmass_g(i) = rootmass_g(i) + rootmass_m(i,m)*farerow(i,m)
           litrmass_g(i) = litrmass_g(i) + litrmass_m(i,m)*farerow(i,m)
           soilcmas_g(i) = soilcmas_g(i) + soilcmas_m(i,m)*farerow(i,m)
c
           burnfrac_g(i) =burnfrac_g(i)+ burnfracrow(i,m)*farerow(i,m) 

           probfire_g(i) =probfire_g(i)+probfirerow(i,m)*farerow(i,m)
           lucemcom_g(i) =lucemcom_g(i)+lucemcomrow(i,m)*farerow(i,m)
           lucltrin_g(i) =lucltrin_g(i)+lucltrinrow(i,m)*farerow(i,m)
           lucsocin_g(i) =lucsocin_g(i)+lucsocinrow(i,m)*farerow(i,m) 
           bterm_g(i)    =bterm_g(i)   +btermrow(i,m)*farerow(i,m) 
           lterm_g(i)    =lterm_g(i)   +ltermrow(i,m)*farerow(i,m) 
           mterm_g(i)    =mterm_g(i)   +mtermrow(i,m)*farerow(i,m)
c                                                   !Rudra added CH4 related variables on 03/12/2013
           CH4WET1_G(i) = CH4WET1_G(i) + CH4WET1ROW(i,m)*farerow(i,m)
           CH4WET2_G(i) = CH4WET2_G(i) + CH4WET2ROW(i,m)*farerow(i,m)
           WETFDYN_G(i) = WETFDYN_G(i) + WETFDYNROW(i,m)*farerow(i,m)
           CH4DYN1_G(i) = CH4DYN1_G(i) + CH4DYN1ROW(i,m)*farerow(i,m)
           CH4DYN2_G(i) = CH4DYN2_G(i) + CH4DYN2ROW(i,m)*farerow(i,m)

           do j=1,icc  

            do k=1,ignd
             rmatctem_g(i,k)=rmatctem_g(i,k)+rmatctemrow(i,m,j,k)
     1                            *farerow(i,m)*fcancmxrow(i,m,j)
            end do

           veghght_g(i) = veghght_g(i) + veghghtrow(i,m,j)
     1                            *farerow(i,m)*fcancmxrow(i,m,j)
           rootdpth_g(i) = rootdpth_g(i) + rootdpthrow(i,m,j)
     1                            *farerow(i,m)*fcancmxrow(i,m,j)
           roottemp_g(i) = roottemp_g(i) + roottemprow(i,m,j)
     1                            *farerow(i,m)*fcancmxrow(i,m,j)
           slai_g(i) = slai_g(i) + slairow(i,m,j)
     1                            *farerow(i,m)*fcancmxrow(i,m,j)

           emit_co2_g(i) =emit_co2_g(i)+ emit_co2row(i,m,j)
     1                            *farerow(i,m)*fcancmxrow(i,m,j)
           emit_co_g(i)  =emit_co_g(i) + emit_corow(i,m,j)
     1                            *farerow(i,m)*fcancmxrow(i,m,j)
           emit_ch4_g(i) =emit_ch4_g(i)+ emit_ch4row(i,m,j)
     1                            *farerow(i,m)*fcancmxrow(i,m,j)
           emit_nmhc_g(i)=emit_nmhc_g(i)+emit_nmhcrow(i,m,j)
     1                            *farerow(i,m)*fcancmxrow(i,m,j)
           emit_h2_g(i)  =emit_h2_g(i) + emit_h2row(i,m,j)
     1                            *farerow(i,m)*fcancmxrow(i,m,j)
           emit_nox_g(i) =emit_nox_g(i)+ emit_noxrow(i,m,j)
     1                            *farerow(i,m)*fcancmxrow(i,m,j)
           emit_n2o_g(i) =emit_n2o_g(i)+ emit_n2orow(i,m,j)
     1                            *farerow(i,m)*fcancmxrow(i,m,j)
           emit_pm25_g(i)=emit_pm25_g(i)+emit_pm25row(i,m,j)
     1                            *farerow(i,m)*fcancmxrow(i,m,j)
           emit_tpm_g(i) =emit_tpm_g(i)+ emit_tpmrow(i,m,j)
     1                            *farerow(i,m)*fcancmxrow(i,m,j)
           emit_tc_g(i)  =emit_tc_g(i) + emit_tcrow(i,m,j)
     1                            *farerow(i,m)*fcancmxrow(i,m,j)
           emit_oc_g(i)  =emit_oc_g(i) + emit_ocrow(i,m,j)
     1                            *farerow(i,m)*fcancmxrow(i,m,j)
           emit_bc_g(i)  =emit_bc_g(i) + emit_bcrow(i,m,j)
     1                            *farerow(i,m)*fcancmxrow(i,m,j)

             ailcg_g(i)=ailcg_g(i)+ailcgrow(i,m,j)*fcancmxrow(i,m,j)
     &                   *farerow(i,m)
             ailcb_g(i)=ailcb_g(i)+ailcbrow(i,m,j)*fcancmxrow(i,m,j)
     &                   *farerow(i,m)

           enddo 
c
852       continue
c
          if ((iyd.ge.jdst).and.(iyd.le.jdend)) then   
c           write to file .CT01D_G               
            write(721,8200)iday,iyear,gpp_g(i),npp_g(i),
     1                nep_g(i),nbp_g(i),autores_g(i),
     2                hetrores_g(i),litres_g(i),socres_g(i),
     3                (dstcemls_g(i)+dstcemls3_g(i)),
     4                litrfall_g(i),humiftrs_g(i)

c           write breakdown of some of basic fluxes to file  
c           *.CT02D_G and selected litter fluxes for selected pft
            write(731,8300)iday,iyear,rml_g(i),rms_g(i),
     1          rmr_g(i),rg_g(i),leaflitr_g(i),tltrleaf_g(i),
     2          tltrstem_g(i),tltrroot_g(i)
c
c           write to file .CT03D_G
            write(741,8400)iday,iyear,vgbiomas_g(i),
     1               gavglai_g(i),gavgltms_g(i),
     2               gavgscms_g(i), totcmass_g(i),
     3               gleafmas_g(i), bleafmas_g(i), stemmass_g(i),
     4               rootmass_g(i), litrmass_g(i), soilcmas_g(i)
c
c           write to file .CT04D_G
            write(751,8500)iday,iyear, ailcg_g(i), 
     1                ailcb_g(i),(rmatctem_g(i,k),k=1,3),
     2                veghght_g(i),rootdpth_g(i),roottemp_g(i),slai_g(i)
c
c           write fire and luc results to file *.CT06D_G
c
           if (dofire .or. lnduseon) then
            write(781,8800)iday,iyear, 
     1          emit_co2_g(i), emit_co_g(i), emit_ch4_g(i),
     2          emit_nmhc_g(i), emit_h2_g(i), emit_nox_g(i),
     3          emit_n2o_g(i), emit_pm25_g(i), emit_tpm_g(i),
     4          emit_tc_g(i), emit_oc_g(i), emit_bc_g(i),
     5          burnfrac_g(i)*100., probfire_g(i),lucemcom_g(i), 
     6          lucltrin_g(i), lucsocin_g(i),
     7          grclarea(i), bterm_g(i), lterm_g(i), mterm_g(i)
           endif
c      
c      write CH4 variables to file *.CT08D_G
c
           if (dowetlands .or. obswetf) then
            write(762,8810)iday,iyear, ch4wet1_g(i), 
     1                 ch4wet2_g(i), wetfdyn_g(i), 
     2                 ch4dyn1_g(i), ch4dyn2_g(i)
           endif  
    
c

            if (compete .or. lnduseon) then 
              sumfare=0.0
              if (mosaic) then
               do m=1,nmos
                 sumfare=sumfare+farerow(i,m)
               enddo
               write(761,8200)iday,iyear,(farerow(i,m)*100.,m=1,nmos),
     1                      sumfare
              else !composite
               do j=1,icc  !m = 1
                 sumfare=sumfare+fcancmxrow(i,1,j)
               enddo
               write(761,8200)iday,iyear,(fcancmxrow(i,1,j)*100.,
     1                      j=1,icc),(1.0-sumfare)*100.,sumfare
              endif !mosaic/composite
            endif !compete/lnduseon
c
          endif !if ((iyd.ge.jdst).and.(iyd.le.jdend)) then  
c
851     continue
c
      endif ! if(ncount.eq.nday) 
      endif ! if(ctem_on)
c
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
        do 861 i=1,nltest

c
          do nt=1,nmon
           if (iday.eq.mmday(nt)) then

            stemmass_mo_g(i)=0.0
            rootmass_mo_g(i)=0.0
            litrmass_mo_g(i)=0.0
            soilcmas_mo_g(i)=0.0
            vgbiomas_mo_g(i)=0.0
            totcmass_mo_g(i)=0.0
           endif
          enddo
c
          if(iday.eq.monthend(imonth+1))then

           laimaxg_mo_g(i)=0.0
           npp_mo_g(i)=0.0
           gpp_mo_g(i)=0.0
           nep_mo_g(i)=0.0
           nbp_mo_g(i)=0.0
           hetrores_mo_g(i)=0.0
           autores_mo_g(i)=0.0
           litres_mo_g(i)=0.0
           soilcres_mo_g(i)=0.0

           emit_co2_mo_g(i)=0.0
           emit_co_mo_g(i) =0.0
           emit_ch4_mo_g(i) =0.0
           emit_nmhc_mo_g(i) =0.0
           emit_h2_mo_g(i) =0.0
           emit_nox_mo_g(i) =0.0
           emit_n2o_mo_g(i) =0.0
           emit_pm25_mo_g(i) =0.0
           emit_tpm_mo_g(i) =0.0
           emit_tc_mo_g(i) =0.0
           emit_oc_mo_g(i) =0.0
           emit_bc_mo_g(i) =0.0
           probfire_mo_g(i) =0.0
           luc_emc_mo_g(i) =0.0
           lucsocin_mo_g(i) =0.0
           lucltrin_mo_g(i) =0.0
           burnfrac_mo_g(i) =0.0
           bterm_mo_g(i)    =0.0
           lterm_mo_g(i)    =0.0
           mterm_mo_g(i)    =0.0
c          CH4(wetland) related variables !Rudra 04/12/2013
           ch4wet1_mo_g(i)  =0.0
           ch4wet2_mo_g(i)  =0.0
           wetfdyn_mo_g(i)  =0.0
           ch4dyn1_mo_g(i)  =0.0
           ch4dyn2_mo_g(i)  =0.0

          endif !mid-month
c
          if (iday .eq. 365) then
           laimaxg_yr_g(i)=0.0
           stemmass_yr_g(i)=0.0
           rootmass_yr_g(i)=0.0
           litrmass_yr_g(i)=0.0
           soilcmas_yr_g(i)=0.0 
           vgbiomas_yr_g(i)=0.0 
           totcmass_yr_g(i)=0.0 
           npp_yr_g(i)=0.0
           gpp_yr_g(i)=0.0
           nep_yr_g(i)=0.0 
           nbp_yr_g(i)=0.0
           hetrores_yr_g(i)=0.0
           autores_yr_g(i)=0.0
           litres_yr_g(i)=0.0
           soilcres_yr_g(i)=0.0
           emit_co2_yr_g(i)=0.0
           emit_co_yr_g(i)=0.0
           emit_ch4_yr_g(i)=0.0
           emit_nmhc_yr_g(i)=0.0
           emit_h2_yr_g(i)=0.0
           emit_nox_yr_g(i)=0.0
           emit_n2o_yr_g(i)=0.0
           emit_pm25_yr_g(i)=0.0
           emit_tpm_yr_g(i)=0.0
           emit_tc_yr_g(i)=0.0
           emit_oc_yr_g(i)=0.0
           emit_bc_yr_g(i)=0.0
           probfire_yr_g(i)=0.0
           luc_emc_yr_g(i)=0.0
           lucsocin_yr_g(i)=0.0
           lucltrin_yr_g(i)=0.0
           burnfrac_yr_g(i)=0.0
           bterm_yr_g(i)=0.0 
           lterm_yr_g(i)=0.0
           mterm_yr_g(i)=0.0
c          CH4(wetland) related variables !Rudra 04/12/2013
           ch4wet1_yr_g(i)  =0.0
           ch4wet2_yr_g(i)  =0.0
           wetfdyn_yr_g(i)  =0.0
           ch4dyn1_yr_g(i)  =0.0
           ch4dyn2_yr_g(i)  =0.0


          endif

861     continue
c
c       accumulate monthly outputs
c
        do 862 i=1,nltest

         do 863 m=1,nmtest

          do j=1,icc

           if (ailcgrow(i,m,j) .gt. laimaxg_mo_m(i,m,j)) then
            laimaxg_mo_m(i,m,j)=ailcgrow(i,m,j)
           end if

           npp_mo_m(i,m,j)=npp_mo_m(i,m,j)+nppvegrow(i,m,j)
           gpp_mo_m(i,m,j)=gpp_mo_m(i,m,j)+gppvegrow(i,m,j) 
           nep_mo_m(i,m,j)=nep_mo_m(i,m,j)+nepvegrow(i,m,j) 
           nbp_mo_m(i,m,j)=nbp_mo_m(i,m,j)+nbpvegrow(i,m,j) 
           hetrores_mo_m(i,m,j)=hetrores_mo_m(i,m,j)
     1                               +hetroresvegrow(i,m,j)
           autores_mo_m(i,m,j) =autores_mo_m(i,m,j)+autoresvegrow(i,m,j)
           litres_mo_m(i,m,j)  =litres_mo_m(i,m,j) +litresvegrow(i,m,j)
           soilcres_mo_m(i,m,j) =soilcres_mo_m(i,m,j) 
     1                               +soilcresvegrow(i,m,j)
           emit_co2_mo_m(i,m,j)=emit_co2_mo_m(i,m,j)+emit_co2row(i,m,j)
           emit_co_mo_m(i,m,j) =emit_co_mo_m(i,m,j)+emit_corow(i,m,j)
           emit_ch4_mo_m(i,m,j) =emit_ch4_mo_m(i,m,j)+emit_ch4row(i,m,j)
           emit_nmhc_mo_m(i,m,j)=emit_nmhc_mo_m(i,m,j)
     1                                     +emit_nmhcrow(i,m,j)
           emit_h2_mo_m(i,m,j) =emit_h2_mo_m(i,m,j)+emit_h2row(i,m,j)
           emit_nox_mo_m(i,m,j) =emit_nox_mo_m(i,m,j)+emit_noxrow(i,m,j)
           emit_n2o_mo_m(i,m,j) =emit_n2o_mo_m(i,m,j)+emit_n2orow(i,m,j)
           emit_pm25_mo_m(i,m,j)=emit_pm25_mo_m(i,m,j)
     1                                     +emit_pm25row(i,m,j)
           emit_tpm_mo_m(i,m,j) =emit_tpm_mo_m(i,m,j)+emit_tpmrow(i,m,j)
           emit_tc_mo_m(i,m,j) =emit_tc_mo_m(i,m,j)+emit_tcrow(i,m,j)
           emit_oc_mo_m(i,m,j) =emit_oc_mo_m(i,m,j)+emit_ocrow(i,m,j)
           emit_bc_mo_m(i,m,j) =emit_bc_mo_m(i,m,j)+emit_bcrow(i,m,j)
           burnfrac_mo_m(i,m,j) =burnfrac_mo_m(i,m,j)+burnvegfrow(i,m,j)

          end do

           nep_mo_m(i,m,iccp1)=nep_mo_m(i,m,iccp1)+nepvegrow(i,m,iccp1) 
           nbp_mo_m(i,m,iccp1)=nbp_mo_m(i,m,iccp1)+nbpvegrow(i,m,iccp1) 
           hetrores_mo_m(i,m,iccp1)=hetrores_mo_m(i,m,iccp1)
     1                               +hetroresvegrow(i,m,iccp1)
           litres_mo_m(i,m,iccp1)  =litres_mo_m(i,m,iccp1) 
     1                               +litresvegrow(i,m,iccp1)
           soilcres_mo_m(i,m,iccp1) =soilcres_mo_m(i,m,iccp1) 
     1                              +soilcresvegrow(i,m,iccp1)

           luc_emc_mo_m(i,m) =luc_emc_mo_m(i,m)
     &                             +lucemcomrow(i,m)
           lucsocin_mo_m(i,m) =lucsocin_mo_m(i,m)
     &                             +lucsocinrow(i,m)
           lucltrin_mo_m(i,m) =lucltrin_mo_m(i,m)
     &                             +lucltrinrow(i,m)
C                         !CH4 related variables !Rudra
           ch4wet1_mo_m(i,m) = ch4wet1_mo_m(i,m) + CH4WET1ROW(i,m)
           ch4wet2_mo_m(i,m) = ch4wet2_mo_m(i,m) + CH4WET2ROW(i,m)
           wetfdyn_mo_m(i,m) = wetfdyn_mo_m(i,m) + WETFDYNROW(i,m)
           ch4dyn1_mo_m(i,m) = ch4dyn1_mo_m(i,m) + CH4DYN1ROW(i,m) 
           ch4dyn2_mo_m(i,m) = ch4dyn2_mo_m(i,m) + CH4DYN2ROW(i,m) 

!          Sum the probfire now, later we will make it a per day value. 
           probfire_mo_m(i,m) =probfire_mo_m(i,m) + probfirerow(i,m) 
           bterm_mo_m(i,m) = bterm_mo_m(i,m) + btermrow(i,m)
           lterm_mo_m(i,m) = lterm_mo_m(i,m) + ltermrow(i,m)
           mterm_mo_m(i,m) = mterm_mo_m(i,m) + mtermrow(i,m)
c
           do 865 nt=1,nmon
c
             if(iday.eq.mmday(nt))then

               do j=1,icc 
                vgbiomas_mo_m(i,m,j)=0.0
                litrmass_mo_m(i,m,j)=0.0
                soilcmas_mo_m(i,m,j)=0.0
                totcmass_mo_m(i,m,j)=0.0
                stemmass_mo_m(i,m,j)=0.0
                rootmass_mo_m(i,m,j)=0.0
               end do
                litrmass_mo_m(i,m,iccp1)=0.0
                soilcmas_mo_m(i,m,iccp1)=0.0

                do 867 j=1,icc

                  vgbiomas_mo_m(i,m,j)=vgbiomas_vegrow(i,m,j)
                  litrmass_mo_m(i,m,j)=litrmassrow(i,m,j)
                  soilcmas_mo_m(i,m,j)=soilcmasrow(i,m,j)
                  stemmass_mo_m(i,m,j)=stemmassrow(i,m,j)
                  rootmass_mo_m(i,m,j)=rootmassrow(i,m,j)
                  totcmass_mo_m(i,m,j)=vgbiomas_vegrow(i,m,j) + 
     1                        litrmassrow(i,m,j)+soilcmasrow(i,m,j)
  
867             continue

                ! Do the bare fraction too
                litrmass_mo_m(i,m,iccp1)=litrmassrow(i,m,iccp1)
                soilcmas_mo_m(i,m,iccp1)=soilcmasrow(i,m,iccp1)

                barefrac=1.0

               do j=1,icc
                vgbiomas_mo_g(i)=vgbiomas_mo_g(i)+vgbiomas_mo_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j) 
                litrmass_mo_g(i)=litrmass_mo_g(i)+litrmass_mo_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j) 
                soilcmas_mo_g(i)=soilcmas_mo_g(i)+soilcmas_mo_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j) 
                stemmass_mo_g(i)=stemmass_mo_g(i)+stemmass_mo_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j) 
                rootmass_mo_g(i)=rootmass_mo_g(i)+rootmass_mo_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j) 
                totcmass_mo_g(i)=totcmass_mo_g(i)+totcmass_mo_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j) 
                barefrac=barefrac-farerow(i,m)*fcancmxrow(i,m,j) 

               end do

!             Also add in the bare fraction contributions.        
              litrmass_mo_g(i)=litrmass_mo_g(i)+litrmass_mo_m(i,m,iccp1)
     &                          *barefrac
              soilcmas_mo_g(i)=soilcmas_mo_g(i)+soilcmas_mo_m(i,m,iccp1)
     &                          *barefrac

             endif ! mmday (mid-month instantaneous value)
c
             if(iday.eq.monthend(nt+1))then

               ndmonth=(monthend(nt+1)-monthend(nt))*nday

               barefrac=1.0
c
               do j=1,icc

                npp_mo_g(i)=npp_mo_g(i)+npp_mo_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j)
                gpp_mo_g(i)=gpp_mo_g(i)+gpp_mo_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j)
                nep_mo_g(i)=nep_mo_g(i)+nep_mo_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j)
                nbp_mo_g(i)=nbp_mo_g(i)+nbp_mo_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j)
                hetrores_mo_g(i)=hetrores_mo_g(i)+hetrores_mo_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j)
                autores_mo_g(i) =autores_mo_g(i) +autores_mo_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j)
                litres_mo_g(i)  =litres_mo_g(i) +litres_mo_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j)
                soilcres_mo_g(i) =soilcres_mo_g(i)+ soilcres_mo_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j)
                laimaxg_mo_g(i)=laimaxg_mo_g(i)+laimaxg_mo_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j)  
                emit_co2_mo_g(i)=emit_co2_mo_g(i)+emit_co2_mo_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j)   
                emit_co_mo_g(i) =emit_co_mo_g(i)+emit_co_mo_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j)   
                emit_ch4_mo_g(i) =emit_ch4_mo_g(i)+emit_ch4_mo_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j)   
                emit_nmhc_mo_g(i)=emit_nmhc_mo_g(i)+
     &                            emit_nmhc_mo_m(i,m,j)*farerow(i,m)
     &                                            *fcancmxrow(i,m,j)   
                emit_h2_mo_g(i) =emit_h2_mo_g(i)+emit_h2_mo_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j)   
                emit_nox_mo_g(i) =emit_nox_mo_g(i)+emit_nox_mo_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j)   
                emit_n2o_mo_g(i) =emit_n2o_mo_g(i)+emit_n2o_mo_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j)   
                emit_pm25_mo_g(i) =emit_pm25_mo_g(i)+
     &                            emit_pm25_mo_m(i,m,j)*farerow(i,m)
     &                                           *fcancmxrow(i,m,j)   
                emit_tpm_mo_g(i) =emit_tpm_mo_g(i)+emit_tpm_mo_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j)   
                emit_tc_mo_g(i) =emit_tc_mo_g(i)+emit_tc_mo_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j)   
                emit_oc_mo_g(i) =emit_oc_mo_g(i)+emit_oc_mo_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j)   
                emit_bc_mo_g(i) =emit_bc_mo_g(i)+emit_bc_mo_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j)  
               burnfrac_mo_g(i)=burnfrac_mo_g(i)+burnfrac_mo_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j)
                barefrac=barefrac-farerow(i,m)*fcancmxrow(i,m,j) 
 
               end do !j

                nep_mo_g(i)=nep_mo_g(i)+nep_mo_m(i,m,iccp1)
     &                          *barefrac
                nbp_mo_g(i)=nbp_mo_g(i)+nbp_mo_m(i,m,iccp1)
     &                          *barefrac
              hetrores_mo_g(i)=hetrores_mo_g(i)+hetrores_mo_m(i,m,iccp1)
     &                          *barefrac
              litres_mo_g(i)  =litres_mo_g(i) +litres_mo_m(i,m,iccp1)
     &                          *barefrac
              soilcres_mo_g(i)=soilcres_mo_g(i)+soilcres_mo_m(i,m,iccp1)
     &                          *barefrac

               luc_emc_mo_g(i) =luc_emc_mo_g(i)
     &                          +luc_emc_mo_m(i,m)*farerow(i,m)   
               lucsocin_mo_g(i) =lucsocin_mo_g(i)
     &                          +lucsocin_mo_m(i,m)*farerow(i,m)   
               lucltrin_mo_g(i) =lucltrin_mo_g(i)
     &                          +lucltrin_mo_m(i,m)*farerow(i,m)
c    CH4(wetland) variables !Rudra 

               ch4wet1_mo_g(i) = ch4wet1_mo_g(i) 
     &                           +ch4wet1_mo_m(i,m)*farerow(i,m)
               ch4wet2_mo_g(i) = ch4wet2_mo_g(i)
     &                           +ch4wet2_mo_m(i,m)*farerow(i,m)

               wetfdyn_mo_m(i,m)=wetfdyn_mo_m(i,m)*(1./
     &                                      real(monthdays(nt))) 

               wetfdyn_mo_g(i) = wetfdyn_mo_g(i)
     &                           +wetfdyn_mo_m(i,m)*farerow(i,m)
               ch4dyn1_mo_g(i) = ch4dyn1_mo_g(i)
     &                           +ch4dyn1_mo_m(i,m)*farerow(i,m)
               ch4dyn2_mo_g(i) = ch4dyn2_mo_g(i)
     &                           +ch4dyn2_mo_m(i,m)*farerow(i,m)

!              Make the probability of fire a per day value
               probfire_mo_m(i,m)=probfire_mo_m(i,m)*
     &                                (1./real(monthdays(nt)))
               probfire_mo_g(i)=probfire_mo_g(i)
     &                          +probfire_mo_m(i,m)*farerow(i,m)   
               bterm_mo_m(i,m)=bterm_mo_m(i,m)*(1./real(monthdays(nt)))
               bterm_mo_g(i) =bterm_mo_g(i)+bterm_mo_m(i,m)*farerow(i,m)  
               lterm_mo_m(i,m)=lterm_mo_m(i,m)*(1./real(monthdays(nt)))
               lterm_mo_g(i) =lterm_mo_g(i)+lterm_mo_m(i,m)*farerow(i,m)  
               mterm_mo_m(i,m)=mterm_mo_m(i,m)*(1./real(monthdays(nt)))
               mterm_mo_g(i) =mterm_mo_g(i)+mterm_mo_m(i,m)*farerow(i,m)  


             endif ! monthend (max lai and accumulated npp/gpp/nep over the whole month)
c                  ! if(iday.eq.monthend(nt+1))

865        continue ! nmon
863      continue ! m
c
         do nt=1,nmon
           if(iday.eq.monthend(nt+1))then
             imonth=nt

                barefrac=1.0  

c            Write to file .CT01M_M/.CT01M_G
              do m=1,nmtest
               do j=1,icc

                  barefrac=barefrac-fcancmxrow(i,m,j)*farerow(i,m)

                if (farerow(i,m)*fcancmxrow(i,m,j) .gt. seed) then
                 write(84,8104)imonth,iyear,laimaxg_mo_m(i,m,j),
     1               vgbiomas_mo_m(i,m,j),litrmass_mo_m(i,m,j),
     2               soilcmas_mo_m(i,m,j),npp_mo_m(i,m,j),
     3               gpp_mo_m(i,m,j),nep_mo_m(i,m,j),
     4               nbp_mo_m(i,m,j),hetrores_mo_m(i,m,j),
     5               autores_mo_m(i,m,j),litres_mo_m(i,m,j),
     6               soilcres_mo_m(i,m,j),
     9               ' TILE ',m,' PFT ',j,' FRAC ',farerow(i,m)*
     a               fcancmxrow(i,m,j)
                end if
               end do !icc

               if (m .eq. nmtest) then
                if (barefrac .gt. seed) then
                write(84,8104)imonth,iyear,0.0,  
     1               0.0,litrmass_mo_m(i,m,iccp1),
     2               soilcmas_mo_m(i,m,iccp1),0.0,
     3               0.0,0.0,
     4               0.0,hetrores_mo_m(i,m,iccp1),
     5               0.0,litres_mo_m(i,m,iccp1),
     6               soilcres_mo_m(i,m,iccp1),
     7               ' TILE ',m,' PFT ',iccp1,' FRAC ',barefrac
                end if
               end if
              end do !m

              write(84,8104)imonth,iyear,laimaxg_mo_g(i),
     1                vgbiomas_mo_g(i),litrmass_mo_g(i),
     2               soilcmas_mo_g(i),npp_mo_g(i),
     3               gpp_mo_g(i),nep_mo_g(i),
     4               nbp_mo_g(i),hetrores_mo_g(i),autores_mo_g(i),
     5               litres_mo_g(i),soilcres_mo_g(i),' GRDAV'

            if (dofire .or. lnduseon) then

c            write to file .CT06M_M/.CT06M_G

              do m=1,nmtest
               do j=1,icc
                if (farerow(i,m)*fcancmxrow(i,m,j) .gt. seed) then
                 write(85,8109)imonth,iyear,emit_co2_mo_m(i,m,j),
     1               emit_co_mo_m(i,m,j),emit_ch4_mo_m(i,m,j),
     2               emit_nmhc_mo_m(i,m,j),emit_h2_mo_m(i,m,j),
     3               emit_nox_mo_m(i,m,j),emit_n2o_mo_m(i,m,j),
     4               emit_pm25_mo_m(i,m,j),emit_tpm_mo_m(i,m,j),
     5               emit_tc_mo_m(i,m,j),emit_oc_mo_m(i,m,j),
     6               emit_bc_mo_m(i,m,j),probfire_mo_m(i,m),
     7               luc_emc_mo_m(i,m),lucltrin_mo_m(i,m),
     8               lucsocin_mo_m(i,m),burnfrac_mo_m(i,m,j)*100.,
     9               bterm_mo_m(i,m),lterm_mo_m(i,m),mterm_mo_m(i,m),
     &               ' TILE ',m,' PFT ',j,' FRAC ',farerow(i,m)*
     &               fcancmxrow(i,m,j)
                end if
               end do
              end do

             write(85,8109)imonth,iyear,emit_co2_mo_g(i),
     3               emit_co_mo_g(i),emit_ch4_mo_g(i),emit_nmhc_mo_g(i),
     4               emit_h2_mo_g(i),emit_nox_mo_g(i),emit_n2o_mo_g(i),
     5               emit_pm25_mo_g(i),emit_tpm_mo_g(i),emit_tc_mo_g(i),
     6               emit_oc_mo_g(i),emit_bc_mo_g(i),
     7               probfire_mo_g(i),luc_emc_mo_g(i),
     8               lucltrin_mo_g(i),lucsocin_mo_g(i),
     8               burnfrac_mo_g(i)*100.,bterm_mo_g(i),lterm_mo_g(i),
     9               mterm_mo_g(i),' GRDAV '

            endif  !dofire/lnduseon

c           add fraction of each pft and bare \\
c
            if (compete .or. lnduseon) then
              sumfare=0.0
              if (mosaic) then
               do m=1,nmos
                 sumfare=sumfare+farerow(i,m)
               enddo
               write(88,8106)imonth,iyear,(farerow(i,m)*100.,m=1,nmos)
     1                      ,sumfare,(pftexistrow(i,j,j),j=1,icc)
              else !composite
               m=1 
               do j=1,icc  
                 sumfare=sumfare+fcancmxrow(i,m,j)
               enddo
               write(88,8106)imonth,iyear,(fcancmxrow(i,m,j)*100.,
     1                      j=1,icc),(1.0-sumfare)*100.,sumfare,
     2                        (pftexistrow(i,m,j),j=1,icc)
              endif !mosaic/composite
            endif !compete/lnduseon
             
             if (dowetlands .or. obswetf) then
             write(91,8111)imonth,iyear,ch4wet1_mo_g(i),
     1                     ch4wet2_mo_g(i),wetfdyn_mo_g(i),
     2                     ch4dyn1_mo_g(i),ch4dyn2_mo_g(i)
             endif 

c
c              initialize monthly accumulated arrays
c              for the next round
 
             do m=1,nmtest

               probfire_mo_m(i,m) =0.0
               luc_emc_mo_m(i,m) =0.0
               lucsocin_mo_m(i,m) =0.0
               lucltrin_mo_m(i,m) =0.0
               bterm_mo_m(i,m) =0.0
               lterm_mo_m(i,m) =0.0
               mterm_mo_m(i,m) =0.0
C       !Rudra
               ch4wet1_mo_m(i,m)  =0.0
               ch4wet2_mo_m(i,m)  =0.0
               wetfdyn_mo_m(i,m)  =0.0
               ch4dyn1_mo_m(i,m)  =0.0
               ch4dyn2_mo_m(i,m)  =0.0
 

             do j=1,icc

              laimaxg_mo_m(i,m,j)=0.0
              hetrores_mo_m(i,m,j)=0.0
              autores_mo_m(i,m,j)=0.0
              litres_mo_m(i,m,j)=0.0
              soilcres_mo_m(i,m,j)=0.0

              npp_mo_m(i,m,j)=0.0
              gpp_mo_m(i,m,j)=0.0
              nep_mo_m(i,m,j)=0.0
              nbp_mo_m(i,m,j)=0.0
              emit_co2_mo_m(i,m,j)=0.0
              emit_co_mo_m(i,m,j) =0.0
              emit_ch4_mo_m(i,m,j) =0.0
              emit_nmhc_mo_m(i,m,j) =0.0
              emit_h2_mo_m(i,m,j) =0.0
              emit_nox_mo_m(i,m,j) =0.0
              emit_n2o_mo_m(i,m,j) =0.0
              emit_pm25_mo_m(i,m,j) =0.0
              emit_tpm_mo_m(i,m,j) =0.0
              emit_tc_mo_m(i,m,j) =0.0
              emit_oc_mo_m(i,m,j) =0.0
              emit_bc_mo_m(i,m,j) =0.0
              burnfrac_mo_m(i,m,j) =0.0
             enddo !j

              hetrores_mo_m(i,m,iccp1)=0.0
              litres_mo_m(i,m,iccp1)=0.0
              soilcres_mo_m(i,m,iccp1)=0.0
              nep_mo_m(i,m,iccp1)=0.0
              nbp_mo_m(i,m,iccp1)=0.0


            enddo !m

           endif ! if(iday.eq.monthend(nt+1))
         enddo ! nt=1,nmon
c
862     continue ! i
c
c       accumulate yearly outputs
c
        do 882 i=1,nltest
          do 883 m=1,nmtest
            do 884 j=1,icc          

             if (ailcgrow(i,m,j).gt.laimaxg_yr_m(i,m,j)) then
               laimaxg_yr_m(i,m,j)=ailcgrow(i,m,j)
             end if

            npp_yr_m(i,m,j)=npp_yr_m(i,m,j)+nppvegrow(i,m,j)
            gpp_yr_m(i,m,j)=gpp_yr_m(i,m,j)+gppvegrow(i,m,j) 
            nep_yr_m(i,m,j)=nep_yr_m(i,m,j)+nepvegrow(i,m,j) 
            nbp_yr_m(i,m,j)=nbp_yr_m(i,m,j)+nbpvegrow(i,m,j) 
            emit_co2_yr_m(i,m,j)=emit_co2_yr_m(i,m,j)+emit_co2row(i,m,j)
            emit_co_yr_m(i,m,j)=emit_co_yr_m(i,m,j)+emit_corow(i,m,j)
            emit_ch4_yr_m(i,m,j)=emit_ch4_yr_m(i,m,j)+emit_ch4row(i,m,j)
            emit_nmhc_yr_m(i,m,j)=emit_nmhc_yr_m(i,m,j)+
     1                            emit_nmhcrow(i,m,j)
            emit_h2_yr_m(i,m,j)=emit_h2_yr_m(i,m,j)+emit_h2row(i,m,j)
            emit_nox_yr_m(i,m,j)=emit_nox_yr_m(i,m,j)+emit_noxrow(i,m,j)
            emit_n2o_yr_m(i,m,j)=emit_n2o_yr_m(i,m,j)+emit_n2orow(i,m,j)
            emit_pm25_yr_m(i,m,j)=emit_pm25_yr_m(i,m,j)+
     1                            emit_pm25row(i,m,j)
            emit_tpm_yr_m(i,m,j)=emit_tpm_yr_m(i,m,j)+emit_tpmrow(i,m,j)
            emit_tc_yr_m(i,m,j)=emit_tc_yr_m(i,m,j)+emit_tcrow(i,m,j)
            emit_oc_yr_m(i,m,j)=emit_oc_yr_m(i,m,j)+emit_ocrow(i,m,j)
            emit_bc_yr_m(i,m,j)=emit_bc_yr_m(i,m,j)+emit_bcrow(i,m,j)

            hetrores_yr_m(i,m,j)=hetrores_yr_m(i,m,j)
     &                                +hetroresvegrow(i,m,j) 
            autores_yr_m(i,m,j)=autores_yr_m(i,m,j)
     &                                +autoresvegrow(i,m,j) 
            litres_yr_m(i,m,j)=litres_yr_m(i,m,j)+litresvegrow(i,m,j) 
            soilcres_yr_m(i,m,j)=soilcres_yr_m(i,m,j)
     &                                +soilcresvegrow(i,m,j) 
            burnfrac_yr_m(i,m,j)=burnfrac_yr_m(i,m,j)+burnvegfrow(i,m,j)

884         continue

!           Also do the bare fraction amounts
            hetrores_yr_m(i,m,iccp1)=hetrores_yr_m(i,m,iccp1)+
     &                                  hetroresvegrow(i,m,iccp1) 
            litres_yr_m(i,m,iccp1)=litres_yr_m(i,m,iccp1)+
     &                                  litresvegrow(i,m,iccp1) 
            soilcres_yr_m(i,m,iccp1)=soilcres_yr_m(i,m,iccp1)+
     &                                  soilcresvegrow(i,m,iccp1) 
            nep_yr_m(i,m,iccp1)=nep_yr_m(i,m,iccp1)+nepvegrow(i,m,iccp1) 
            nbp_yr_m(i,m,iccp1)=nbp_yr_m(i,m,iccp1)+nbpvegrow(i,m,iccp1) 

            probfire_yr_m(i,m)=probfire_yr_m(i,m)
     &                         +(probfirerow(i,m) * (1./365.))
            bterm_yr_m(i,m)=bterm_yr_m(i,m)+(btermrow(i,m)*(1./365.))  
            lterm_yr_m(i,m)=lterm_yr_m(i,m)+(ltermrow(i,m)*(1./365.))
            mterm_yr_m(i,m)=mterm_yr_m(i,m)+(mtermrow(i,m)*(1./365.))
            luc_emc_yr_m(i,m)=luc_emc_yr_m(i,m)+lucemcomrow(i,m)
            lucsocin_yr_m(i,m)=lucsocin_yr_m(i,m)+lucsocinrow(i,m)
            lucltrin_yr_m(i,m)=lucltrin_yr_m(i,m)+lucltrinrow(i,m)
c             CH4(wetland) variables !Rudra 

               ch4wet1_yr_m(i,m) = ch4wet1_yr_m(i,m)
     &                           +ch4wet1row(i,m)
               ch4wet2_yr_m(i,m) = ch4wet2_yr_m(i,m)
     &                           +ch4wet2row(i,m)
               wetfdyn_yr_m(i,m) = wetfdyn_yr_m(i,m)
     &                           +(wetfdynrow(i,m)*(1./365.))
               ch4dyn1_yr_m(i,m) = ch4dyn1_yr_m(i,m)
     &                           +ch4dyn1row(i,m)
               ch4dyn2_yr_m(i,m) = ch4dyn2_yr_m(i,m)
     &                           +ch4dyn2row(i,m)


            if (iday.eq.365) then

              do 885 j=1,icc
                stemmass_yr_m(i,m,j)=stemmassrow(i,m,j)
                rootmass_yr_m(i,m,j)=rootmassrow(i,m,j)
                litrmass_yr_m(i,m,j)=litrmassrow(i,m,j)
                soilcmas_yr_m(i,m,j)=soilcmasrow(i,m,j)
                vgbiomas_yr_m(i,m,j)=vgbiomas_vegrow(i,m,j)
                totcmass_yr_m(i,m,j)=vgbiomas_yr_m(i,m,j)+
     &                               litrmass_yr_m(i,m,j)+
     &                               soilcmas_yr_m(i,m,j)

885           continue

                litrmass_yr_m(i,m,iccp1)=litrmassrow(i,m,iccp1)
                soilcmas_yr_m(i,m,iccp1)=soilcmasrow(i,m,iccp1)

                barefrac=1.0

              do j=1,icc
                laimaxg_yr_g(i)=laimaxg_yr_g(i)+ laimaxg_yr_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j) 
                stemmass_yr_g(i)=stemmass_yr_g(i)+stemmass_yr_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j) 
                rootmass_yr_g(i)=rootmass_yr_g(i)+rootmass_yr_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j) 
                litrmass_yr_g(i)=litrmass_yr_g(i)+litrmass_yr_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j) 
                soilcmas_yr_g(i)=soilcmas_yr_g(i)+soilcmas_yr_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j)   
                vgbiomas_yr_g(i)=vgbiomas_yr_g(i)+vgbiomas_yr_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j) 
                totcmass_yr_g(i)=totcmass_yr_g(i)+totcmass_yr_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j) 
                npp_yr_g(i)=npp_yr_g(i)+npp_yr_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j)
                gpp_yr_g(i)=gpp_yr_g(i)+gpp_yr_m(i,m,j)    
     &                          *farerow(i,m)*fcancmxrow(i,m,j)
                nep_yr_g(i)=nep_yr_g(i)+nep_yr_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j)
                nbp_yr_g(i)=nbp_yr_g(i)+nbp_yr_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j)
                emit_co2_yr_g(i)=emit_co2_yr_g(i)+emit_co2_yr_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j)
                emit_co_yr_g(i)=emit_co_yr_g(i)+emit_co_yr_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j)
                emit_ch4_yr_g(i)=emit_ch4_yr_g(i)+emit_ch4_yr_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j)
               emit_nmhc_yr_g(i)=emit_nmhc_yr_g(i)+emit_nmhc_yr_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j)
                emit_h2_yr_g(i)=emit_h2_yr_g(i)+emit_h2_yr_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j)
                emit_nox_yr_g(i)=emit_nox_yr_g(i)+emit_nox_yr_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j)
                emit_n2o_yr_g(i)=emit_n2o_yr_g(i)+emit_n2o_yr_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j)
               emit_pm25_yr_g(i)=emit_pm25_yr_g(i)+emit_pm25_yr_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j)
                emit_tpm_yr_g(i)=emit_tpm_yr_g(i)+emit_tpm_yr_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j)
                emit_tc_yr_g(i)=emit_tc_yr_g(i)+emit_tc_yr_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j)
                emit_oc_yr_g(i)=emit_oc_yr_g(i)+emit_oc_yr_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j)
                emit_bc_yr_g(i)=emit_bc_yr_g(i)+emit_bc_yr_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j)

                hetrores_yr_g(i)=hetrores_yr_g(i)+hetrores_yr_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j)
                autores_yr_g(i) =autores_yr_g(i) +autores_yr_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j)
                litres_yr_g(i)  =litres_yr_g(i)  +litres_yr_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j)
                soilcres_yr_g(i) =soilcres_yr_g(i) +soilcres_yr_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j)
                burnfrac_yr_g(i)=burnfrac_yr_g(i)+burnfrac_yr_m(i,m,j)
     &                          *farerow(i,m)*fcancmxrow(i,m,j)

                barefrac=barefrac-farerow(i,m)*fcancmxrow(i,m,j) 

              end do !j

              litrmass_yr_g(i)=litrmass_yr_g(i)+litrmass_yr_m(i,m,iccp1)
     &                          *barefrac
              soilcmas_yr_g(i)=soilcmas_yr_g(i)+soilcmas_yr_m(i,m,iccp1)
     &                          *barefrac
              hetrores_yr_g(i)=hetrores_yr_g(i)+hetrores_yr_m(i,m,iccp1)
     &                          *barefrac
              litres_yr_g(i)  =litres_yr_g(i)  +litres_yr_m(i,m,iccp1)
     &                          *barefrac
              soilcres_yr_g(i)=soilcres_yr_g(i)+soilcres_yr_m(i,m,iccp1)
     &                          *barefrac
                nep_yr_g(i)=nep_yr_g(i)+nep_yr_m(i,m,j)
     &                          *barefrac
                nbp_yr_g(i)=nbp_yr_g(i)+nbp_yr_m(i,m,j)
     &                          *barefrac
   

              probfire_yr_g(i)=probfire_yr_g(i)
     &                         +probfire_yr_m(i,m)*farerow(i,m)   
              bterm_yr_g(i)=bterm_yr_g(i)+bterm_yr_m(i,m)*farerow(i,m)   
              lterm_yr_g(i)=lterm_yr_g(i)+lterm_yr_m(i,m)*farerow(i,m)   
              mterm_yr_g(i)=mterm_yr_g(i)+mterm_yr_m(i,m)*farerow(i,m)   
              luc_emc_yr_g(i)=luc_emc_yr_g(i)
     &                         +luc_emc_yr_m(i,m)*farerow(i,m)   
              lucsocin_yr_g(i)=lucsocin_yr_g(i)
     &                         +lucsocin_yr_m(i,m)*farerow(i,m)   
              lucltrin_yr_g(i)=lucltrin_yr_g(i)
     &                         +lucltrin_yr_m(i,m)*farerow(i,m) 
c    CH4(wetland) variables !Rudra 

               ch4wet1_yr_g(i) = ch4wet1_yr_g(i)
     &                           +ch4wet1_yr_m(i,m)*farerow(i,m)
               ch4wet2_yr_g(i) = ch4wet2_yr_g(i)
     &                           +ch4wet2_yr_m(i,m)*farerow(i,m)
               wetfdyn_yr_g(i) = wetfdyn_yr_g(i)
     &                           +wetfdyn_yr_m(i,m)*farerow(i,m)
               ch4dyn1_yr_g(i) = ch4dyn1_yr_g(i)
     &                           +ch4dyn1_yr_m(i,m)*farerow(i,m)
               ch4dyn2_yr_g(i) = ch4dyn2_yr_g(i)
     &                           +ch4dyn2_yr_m(i,m)*farerow(i,m)
  
c
            endif ! iday 365
c
883       continue ! m
c
          if (iday.eq.365) then

              barefrac=1.0

c            Write to file .CT01Y_M/.CT01Y_G

           do m=1,nmtest
            do j=1,icc

                barefrac=barefrac-fcancmxrow(i,m,j)*farerow(i,m)

             if (farerow(i,m)*fcancmxrow(i,m,j) .gt. seed) then
              write(86,8105)iyear,laimaxg_yr_m(i,m,j),
     1            vgbiomas_yr_m(i,m,j),stemmass_yr_m(i,m,j),
     2            rootmass_yr_m(i,m,j),litrmass_yr_m(i,m,j),
     3            soilcmas_yr_m(i,m,j),totcmass_yr_m(i,m,j),
     4            npp_yr_m(i,m,j),gpp_yr_m(i,m,j),nep_yr_m(i,m,j),
     5            nbp_yr_m(i,m,j),hetrores_yr_m(i,m,j),
     6            autores_yr_m(i,m,j),litres_yr_m(i,m,j),
     9            soilcres_yr_m(i,m,j),' TILE ',m,' PFT ',j,' FRAC '
     a            ,farerow(i,m)*fcancmxrow(i,m,j)
             end if
            end do !j

!           Now do the bare fraction of the grid cell. Only soil c, hetres
!           and litter are relevant so the rest are set to 0. 
            if (m .eq. nmtest) then
             if (barefrac .gt. seed) then
              write(86,8105)iyear,0.,
     1            0., 
     1            0.,0.,
     2            litrmassrow(i,m,iccp1),soilcmasrow(i,m,iccp1),
     3            0.+soilcmasrow(i,m,iccp1)+
     2            litrmassrow(i,m,iccp1),0.,
     4            0.,0.,
     5            0.,hetrores_yr_m(i,m,iccp1),
     6            0.,litres_yr_m(i,m,iccp1),soilcres_yr_m(i,m,iccp1),
     7            ' TILE ',m,' PFT ',iccp1,' FRAC ',barefrac
             end if
            end if

           end do !m

             write(86,8105)iyear,laimaxg_yr_g(i),vgbiomas_yr_g(i),
     1            stemmass_yr_g(i),rootmass_yr_g(i),litrmass_yr_g(i),
     2            soilcmas_yr_g(i),totcmass_yr_g(i),npp_yr_g(i),
     3            gpp_yr_g(i),nep_yr_g(i),
     4            nbp_yr_g(i),hetrores_yr_g(i),autores_yr_g(i),
     5            litres_yr_g(i),soilcres_yr_g(i),' GRDAV'

           if (dofire .or. lnduseon) then
c            Write to file .CT06Y_M/.CT06Y_G
            do m=1,nmtest
             do j=1,icc
              if (farerow(i,m)*fcancmxrow(i,m,j) .gt. seed) then
               write(87,8108)iyear,emit_co2_yr_m(i,m,j),
     1            emit_co_yr_m(i,m,j),emit_ch4_yr_m(i,m,j),
     2            emit_nmhc_yr_m(i,m,j),emit_h2_yr_m(i,m,j),
     3            emit_nox_yr_m(i,m,j),emit_n2o_yr_m(i,m,j),
     4            emit_pm25_yr_m(i,m,j),emit_tpm_yr_m(i,m,j),
     5            emit_tc_yr_m(i,m,j),emit_oc_yr_m(i,m,j),
     6            emit_bc_yr_m(i,m,j),probfire_yr_m(i,m),
     7            luc_emc_yr_m(i,m),lucltrin_yr_m(i,m),
     8            lucsocin_yr_m(i,m),burnfrac_yr_m(i,m,j)*100.,
     9            bterm_yr_m(i,m),lterm_yr_m(i,m),mterm_yr_m(i,m),
     9            ' TILE ',m,' PFT ',j,' FRAC '
     a            ,farerow(i,m)*fcancmxrow(i,m,j)
             end if
            end do
           end do

             write(87,8108)iyear,emit_co2_yr_g(i),
     4            emit_co_yr_g(i),emit_ch4_yr_g(i),emit_nmhc_yr_g(i),
     5            emit_h2_yr_g(i),emit_nox_yr_g(i),emit_n2o_yr_g(i),
     6            emit_pm25_yr_g(i),emit_tpm_yr_g(i),emit_tc_yr_g(i),
     7            emit_oc_yr_g(i),emit_bc_yr_g(i),probfire_yr_g(i),
     8            luc_emc_yr_g(i),lucltrin_yr_g(i),
     9            lucsocin_yr_g(i),burnfrac_yr_g(i)*100.,bterm_yr_g(i),
     a            lterm_yr_g(i),mterm_yr_g(i), ' GRDAV'

           endif !dofire,lnduseon

c           write fraction of each pft and bare 
c
             if (compete .or. lnduseon) then
                 sumfare=0.0
               if (mosaic) then
                 do m=1,nmos
                    sumfare=sumfare+farerow(i,m)
                 enddo
                 write(89,8107)iyear,(farerow(i,m)*100.,m=1,nmos),
     &                         sumfare,(pftexistrow(i,j,j),j=1,icc)
               else  !composite
                 m=1 
                 do j=1,icc  
                    sumfare=sumfare+fcancmxrow(i,m,j)
                 enddo
                write(89,8107)iyear,(fcancmxrow(i,m,j)*100.,
     1                       j=1,icc),(1.0-sumfare)*100.,sumfare,
     2                      (pftexistrow(i,m,j),j=1,icc)

               endif
             endif !compete/lnduseon
C            
              if (dowetlands .or. obswetf) then 
                write(92,8115)iyear,ch4wet1_yr_g(i),
     1                     ch4wet2_yr_g(i),wetfdyn_yr_g(i),
     2                     ch4dyn1_yr_g(i),ch4dyn2_yr_g(i)
              endif 



c             initialize yearly accumulated arrays
c             for the next round

             do m=1,nmtest
              probfire_yr_m(i,m)=0.0
              luc_emc_yr_m(i,m)=0.0
              lucsocin_yr_m(i,m)=0.0
              lucltrin_yr_m(i,m)=0.0
              bterm_yr_m(i,m)=0.0
              lterm_yr_m(i,m)=0.0
              mterm_yr_m(i,m)=0.0
C       !Rudra
               ch4wet1_yr_m(i,m)  =0.0
               ch4wet2_yr_m(i,m)  =0.0
               wetfdyn_yr_m(i,m)  =0.0
               ch4dyn1_yr_m(i,m)  =0.0
               ch4dyn2_yr_m(i,m)  =0.0


               do j = 1, icc 
                laimaxg_yr_m(i,m,j)=0.0
                npp_yr_m(i,m,j)=0.0
                gpp_yr_m(i,m,j)=0.0
                nep_yr_m(i,m,j)=0.0 
                nbp_yr_m(i,m,j)=0.0 
                emit_co2_yr_m(i,m,j)=0.0
                emit_co_yr_m(i,m,j)=0.0
                emit_ch4_yr_m(i,m,j)=0.0
                emit_nmhc_yr_m(i,m,j)=0.0
                emit_h2_yr_m(i,m,j)=0.0
                emit_nox_yr_m(i,m,j)=0.0
                emit_n2o_yr_m(i,m,j)=0.0
                emit_pm25_yr_m(i,m,j)=0.0
                emit_tpm_yr_m(i,m,j)=0.0
                emit_tc_yr_m(i,m,j)=0.0
                emit_oc_yr_m(i,m,j)=0.0
                emit_bc_yr_m(i,m,j)=0.0
                hetrores_yr_m(i,m,j)=0.0 
                autores_yr_m(i,m,j)=0.0 
                litres_yr_m(i,m,j)=0.0 
                soilcres_yr_m(i,m,j)=0.0 
                burnfrac_yr_m(i,m,j)=0.0
              
               enddo
                hetrores_yr_m(i,m,iccp1)=0.0  
                litres_yr_m(i,m,iccp1)=0.0 
                soilcres_yr_m(i,m,iccp1)=0.0 
                nep_yr_m(i,m,iccp1)=0.0 
                nbp_yr_m(i,m,iccp1)=0.0 
             enddo

            endif ! if iday=365
c
882     continue ! i
c
      endif ! if(ncount.eq.nday)
      endif ! if(ctem_on) 
C
8104  FORMAT(1X,I4,I5,12(F10.3,1X),2(A6,I2),A6,F8.2)
8105  FORMAT(1X,I5,15(F10.3,1X),2(A6,I2),A6,F8.2)
8106  FORMAT(1X,I4,I5,11(F10.5,1X),9L5,2(A6,I2))
8107  FORMAT(1X,I5,11(F10.5,1X),9L5,2(A6,I2))
8108  FORMAT(1X,I5,20(F10.3,1X),2(A6,I2),A6,F8.2)
8109  FORMAT(1X,I4,I5,20(F10.3,1X),2(A6,I2),A6,F8.2)
8111  FORMAT(1X,I4,I5,5(F10.3,1X),2(A6,I2))
8115  FORMAT(1X,I5,5(F10.3,1X),2(A6,I2))
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
        WRITE(100,5020)DEGLAT,DEGLON,ZRFMGRD(1),ZRFHGRD(1),ZBLDGRD(1),
     1                 GCGRD(1),NLTEST,NMTEST
        DO I=1,NLTEST
          DO M=1,NMTEST

C         IF START_BARE (SO EITHER COMPETE OR LNDUSEON), THEN WE NEED TO CREATE
C         THE FCANROW FOR THE RS FILE.
          IF (START_BARE .AND. MOSAIC) THEN
           IF (M .LE. 2) THEN                     !NDL
            FCANROW(I,M,1)=1.0
           ELSEIF (M .GE. 3 .AND. M .LE. 5) THEN  !BDL
            FCANROW(I,M,2)=1.0
           ELSEIF (M .EQ. 6 .OR. M .EQ. 7) THEN  !CROP
            FCANROW(I,M,3)=1.0
           ELSEIF (M .EQ. 8 .OR. M .EQ. 9) THEN  !GRASSES
            FCANROW(I,M,4)=1.0
           ELSE                                  !BARE 
            FCANROW(I,M,5)=1.0
           ENDIF
          ENDIF !START_BARE/MOSAIC
             
            WRITE(100,5040) (FCANROW(I,M,J),J=1,ICAN+1),(PAMXROW(I,M,J),
     1                      J=1,ICAN)
            WRITE(100,5040) (LNZ0ROW(I,M,J),J=1,ICAN+1),(PAMNROW(I,M,J),
     1                      J=1,ICAN)
            WRITE(100,5040) (ALVCROW(I,M,J),J=1,ICAN+1),(CMASROW(I,M,J),
     1                      J=1,ICAN)
            WRITE(100,5040) (ALICROW(I,M,J),J=1,ICAN+1),(ROOTROW(I,M,J),
     1                      J=1,ICAN)
            WRITE(100,5030) (RSMNROW(I,M,J),J=1,ICAN),
     1                      (QA50ROW(I,M,J),J=1,ICAN)
            WRITE(100,5030) (VPDAROW(I,M,J),J=1,ICAN),
     1                      (VPDBROW(I,M,J),J=1,ICAN)
            WRITE(100,5030) (PSGAROW(I,M,J),J=1,ICAN),
     1                      (PSGBROW(I,M,J),J=1,ICAN)
            WRITE(100,5040) DRNROW(I,M),SDEPROW(I,M),FAREROW(I,M)
            WRITE(100,5090) XSLPROW(I,M),GRKFROW(I,M),WFSFROW(I,M),
     1                      WFCIROW(I,M),MIDROW(I,M)
            WRITE(100,5080) (SANDROW(I,M,J),J=1,3)
            WRITE(100,5080) (CLAYROW(I,M,J),J=1,3)
            WRITE(100,5080) (ORGMROW(I,M,J),J=1,3)
C           Temperatures are in degree C
            IF (TCANROW(I,M).NE.0.0) TCANRS(I,M)=TCANROW(I,M)-273.16
            IF (TSNOROW(I,M).NE.0.0) TSNORS(I,M)=TSNOROW(I,M)-273.16
            IF (TPNDROW(I,M).NE.0.0) TPNDRS(I,M)=TPNDROW(I,M)-273.16
            WRITE(100,5050) (TBARROW(I,M,J)-273.16,J=1,3),TCANRS(I,M),
     2                      TSNORS(I,M),TPNDRS(I,M)
            WRITE(100,5060) (THLQROW(I,M,J),J=1,3),(THICROW(I,M,J),
     1                      J=1,3),ZPNDROW(I,M)
C
            WRITE(100,5070) RCANROW(I,M),SCANROW(I,M),SNOROW(I,M),
     1                      ALBSROW(I,M),RHOSROW(I,M),GROROW(I,M)
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
          open(unit=101,file=argbuff(1:strlen(argbuff))//'.CTM_RS')
c
          write(101,7010) titlec1
          write(101,7010) titlec2
          write(101,7010) titlec3

c         if landuseon or competition, then we need to recreate the 
c         dvdfcanrow so do so now
          if (lnduseon .or. compete ) then
           icountrow=0
           do j = 1, ican
            do i=1,nltest
             do m=1,nmtest 
               k1c = (j-1)*l2max + 1
               k2c = k1c + (l2max - 1)
c
                do n = k1c, k2c
                 if(modelpft(n).eq.1)then
                  icountrow(i,m) = icountrow(i,m) + 1
                  if (fcanrow(i,m,j) .gt. 0.) then
                   dvdfcanrow(i,m,icountrow(i,m))=
     1                 fcancmxrow(i,m,icountrow(i,m))/fcanrow(i,m,j)
                  else
                   dvdfcanrow(i,m,icountrow(i,m))=0.
                  end if
                 endif !modelpft
                enddo !n
              
c            check to ensure that the dvdfcanrow's add up to 1 across a 
c            class-level pft
             if(dvdfcanrow(i,m,1) .eq. 0. .and. dvdfcanrow(i,m,2) 
     1                                              .eq. 0.) then
                  dvdfcanrow(i,m,1)=1.0
             elseif(dvdfcanrow(i,m,3) .eq. 0. .and. dvdfcanrow(i,m,4) 
     1                     .eq. 0. .and. dvdfcanrow(i,m,5) .eq. 0.) then
                  dvdfcanrow(i,m,3)=1.0
             elseif(dvdfcanrow(i,m,6) .eq. 0. .and. dvdfcanrow(i,m,7) 
     1                                              .eq. 0.) then
                  dvdfcanrow(i,m,6)=1.0
             elseif(dvdfcanrow(i,m,8) .eq. 0. .and. dvdfcanrow(i,m,9) 
     1                                              .eq. 0.) then
                  dvdfcanrow(i,m,8)=1.0
             endif
  
             enddo !m
            enddo !i 
           enddo !j

           do i=1,nltest
            do m=1,nmtest 
             do j = 1, icc
c            lastly check if the different pfts accidently add up > 1.0
c            after rounding to the number of sig figs used in the output
c            this rounds to 2 decimal places. if you are found to be over
c            or under, arbitrarily reduce one of the pfts. the amount of
c            the change will be inconsequential. 
              rnded_pft(j) =real(int(dvdfcanrow(i,m,j) * 1000.0 + 5.0))
     1                                                         / 1000.0
             enddo

             if (rnded_pft(1) + rnded_pft(2) .ne. 1.0) then
              dvdfcanrow(i,m,1) = 1.0 - rnded_pft(2)
              dvdfcanrow(i,m,2) = rnded_pft(2)
             else if (rnded_pft(3) + rnded_pft(4) + rnded_pft(5) 
     1                                                 .ne. 1.0) then
              dvdfcanrow(i,m,3) = 1.0 - rnded_pft(4) - rnded_pft(5)
              dvdfcanrow(i,m,4) = rnded_pft(4)
              dvdfcanrow(i,m,5) = rnded_pft(5)
             else if (rnded_pft(6) + rnded_pft(7) .ne. 1.0) then
              dvdfcanrow(i,m,6) = 1.0 - rnded_pft(7)
              dvdfcanrow(i,m,7) = rnded_pft(7)
             else if (rnded_pft(8) + rnded_pft(9) .ne. 1.0) then
              dvdfcanrow(i,m,8) = 1.0 - rnded_pft(9)
              dvdfcanrow(i,m,9) = rnded_pft(9)
             endif
            enddo
           enddo


          endif !lnuse/compete
c
          do i=1,nltest
            do m=1,nmtest
              write(101,7011) (ailcminrow(i,m,j),j=1,icc)
              write(101,7011) (ailcmaxrow(i,m,j),j=1,icc)
              write(101,'(9f8.3)') (dvdfcanrow(i,m,j),j=1,icc)
c
              write(101,7011) (gleafmasrow(i,m,j),j=1,icc)
              write(101,7011) (bleafmasrow(i,m,j),j=1,icc)
              write(101,7011) (stemmassrow(i,m,j),j=1,icc)
              write(101,7011) (rootmassrow(i,m,j),j=1,icc)
              write(101,7013) (litrmassrow(i,m,j),j=1,iccp1)
              write(101,7013) (soilcmasrow(i,m,j),j=1,iccp1)
              write(101,7012) (lfstatusrow(i,m,j),j=1,icc)
              write(101,7012) (pandaysrow(i,m,j),j=1,icc)
            enddo
c
            write(101,"(6f8.3)") (mlightnggrd(i,j),j=1,6)  !mean monthly lightning frequency
            write(101,"(6f8.3)") (mlightnggrd(i,j),j=7,12) !flashes/km2.year
            write(101,"(f8.2)") extnprobgrd(i)
            write(101,"(f8.2)") prbfrhucgrd(i)
            write(101,"(i4)") stdalngrd(i)

            if (compete) then
             write(101,"(5f8.2)")twarmm(i),tcoldm(i),gdd5(i),
     1                            aridity(i),srplsmon(i)
             write(101,"(5f8.2)")defctmon(i),anndefct(i),annsrpls(i),
     1                        annpcp(i),dry_season_length(i)
            endif

            if (dowetlands) then     
              write(101,"(8f9.5)")(wetfrac_sgrd(i,j),j=1,8)
            endif   

          enddo
c

          close(101)
        endif ! ctem_on
c
       endif ! if iday=365
      endif ! if generate restart files
c
7011  format(9f8.2)
7012  format(9i8)
7013  format(10f8.2)
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
C
C     FLAG! I CAN FIND NO USE OF THIS ISUM VAR. JM DEC 11 2012

C      DO 825 J=1,6
C          ISUM(J)=0
C825   CONTINUE
C
C      DO 900 I=1,50
C        DO 810 J=1,6
C          ISUM(J)=ISUM(J)+ITCTGAT(1,J,I)
C810     CONTINUE
C900   CONTINUE
C
C===================== CTEM =====================================\
c
c      checking the time spent for running model
c
c      call idate(today) 
c      call itime(now)
c      write(*,1001) today(2), today(1), 2000+today(3), now
c 1001 format( 'end date: ', i2.2, '/', i2.2, '/', i4.4, 
c     &      '; end time: ', i2.2, ':', i2.2, ':', i2.2 )
c
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
c       then ctem ones
        close(71)
        close(711)
        close(721)
        close(731)
        close(741)
        close(751)

        if (mosaic) then
         close(72)
         close(73)
         close(74)
         close(75)
         close(76)
         if (dofire .or. lnduseon) then
          close(78)
         end if
        end if
c
        if (compete .or. lnduseon) then
          close(761)
        endif
c
        if (dowetlands .or. obswetf) then
        close(762)
        endif 
c
       if (dofire .or. lnduseon) then
        close(781)
       endif
      endif ! if (.not. parallelrun) 
C
C     CLOSE CLASS OUTPUT FILES      
      CLOSE(81)
      CLOSE(82)
      CLOSE(83)
c     then the CTEM ones
      close(84)
      close(86)
      if (dofire .or. lnduseon) then
       close(85)
       close(87)
      endif
      if (compete .or. lnduseon) then
       close(88)
       close(89)
      endif
 
      if (dowetlands .or. obswetf) then 
       close(91)
       close(92)
      endif 
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
c       then ctem ones
        close(71)
        close(711)
        close(721)
        close(731)
        close(741)
        close(751)

        if (mosaic) then
         close(72)
         close(73)
         close(74)
         close(75)
         close(76)
         if (dofire .or. lnduseon) then
          close(78)
         end if
        end if

c
        if (compete .or. lnduseon) then
          close(761)
        endif
c
       if (dofire .or. lnduseon) then
        close(781)
       endif
      endif ! if (.not. parallelrun) 
C
C     CLOSE CLASS OUTPUT FILES      
      CLOSE(81)
      CLOSE(82)
      CLOSE(83)
c     then the CTEM ones
      close(84)
      close(86)
      if (dofire .or. lnduseon) then
       close(85)
       close(87)
      endif
      if (compete .or. lnduseon) then
       close(88)
       close(89)
      endif
      if (dowetlands .or. obswetf) then 
       close(91)
       close(92)
      endif 

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
